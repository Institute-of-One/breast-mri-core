"""Reference-ROI-conditioned FTVser reconstruction and frozen-model evaluation."""
import argparse
import hashlib
import json
import numpy as np
import pydicom
from theme3_start import ROOT
from theme3_lock_development import OUT


def load_case(pid,visit):
    extraction=json.loads((OUT/f'{pid}_extraction.json').read_text())
    v=extraction['visits'][visit]
    pair=v['pair_file']
    if extraction['manufacturer']=='Philips Medical Systems':
        correction=json.loads((OUT/f'{pid}_{visit}_native_precision.json').read_text())
        assert correction['status']=='pass';pair=correction['pair_file']
    a=np.load(ROOT/pair)
    masksource=next(r for r in json.loads((OUT/'ftv_mask_sources.json').read_text()) if r['patient_id']==pid and r['visit']==visit)
    mask=pydicom.dcmread(ROOT/masksource['path'])
    old=json.loads((ROOT/f'data/derived/training_cohort/patients/{pid}.json').read_text())
    ref=pydicom.dcmread(ROOT/old['visits'][visit]['inputs']['ser'],stop_before_pixels=True)
    assert mask.PatientID==ref.PatientID and mask.StudyInstanceUID==ref.StudyInstanceUID
    assert mask.FrameOfReferenceUID==ref.FrameOfReferenceUID
    assert str(mask.ClinicalTrialTimePointID)==visit and str(mask.PatientID)==pid
    assert len(mask.SegmentSequence)==1 and mask.SegmentSequence[0].SegmentLabel=='VOLSER Analysis Mask'
    shared=mask.SharedFunctionalGroupsSequence[0]
    assert np.allclose(shared.PlaneOrientationSequence[0].ImageOrientationPatient,a['orientation'])
    assert np.allclose(shared.PixelMeasuresSequence[0].PixelSpacing,a['pixel_spacing'])
    data=mask.pixel_array
    assert data.shape[1:]==a['pre'].shape[1:]
    grid=np.empty(a['pre'].shape,dtype=np.uint8);seen=set()
    assert len(mask.PerFrameFunctionalGroupsSequence)==len(data)==len(grid)
    for frame,plane in zip(mask.PerFrameFunctionalGroupsSequence,data):
        assert int(frame.SegmentIdentificationSequence[0].ReferencedSegmentNumber)==1
        pos=frame.PlanePositionSequence[0].ImagePositionPatient
        match=np.flatnonzero(np.all(np.isclose(a['positions'],pos,atol=.001,rtol=0),axis=1))
        assert len(match)==1 and int(match[0]) not in seen
        i=int(match[0]);seen.add(i);grid[i]=plane
    assert len(seen)==len(grid)
    normal=np.cross(a['orientation'][:3],a['orientation'][3:])
    gaps=np.diff(a['positions']@normal)
    assert np.all(gaps>0) and np.allclose(gaps,gaps[0],atol=.001,rtol=0)
    voxel_cc=float(np.prod(a['pixel_spacing'])*gaps[0]/1000)
    values={}
    for item in ref[(0x0117,0x10b0)].value:
        label=str(item[(0x0117,0x10b5)].value)
        kind='ser' if label.startswith('FTVser(') else 'pe' if label.startswith('FTVpe(') else None
        if kind:
            assert kind not in values
            values[kind]={'label':label,'lower':float(item[(0x0117,0x10b1)].value),
                'upper':float(item[(0x0117,0x10b2)].value),'count':int(item[(0x0117,0x10b3)].value),
                'volume_cc':float(item[(0x0117,0x10b4)].value)}
    assert set(values)=={'pe','ser'}
    pre,early,late=[a[k].astype(float) for k in ['pre','early','late']]
    return pre,early,late,grid==0,voxel_cc,values,v,sorted(int(x) for x in np.unique(grid))


def volume(pre,early,late,eligible,voxel_cc,lower,upper):
    denominator=late-pre
    ratio=np.divide(early-pre,denominator,out=np.full_like(pre,np.nan),where=denominator!=0)
    capped=np.clip(ratio,0,upper)
    selected=eligible & np.isfinite(capped) & (capped>lower)
    finite_values=capped[eligible & np.isfinite(capped)]
    return {'count':int(selected.sum()),'volume_cc':float(selected.sum()*voxel_cc),
            'eligible_voxels':int(eligible.sum()),
            'fraction_above_threshold':float(selected.sum()/eligible.sum()) if eligible.any() else None,
            'capped_SER_quantiles_0_25_50_75_100':np.quantile(finite_values,[0,.25,.5,.75,1]).tolist() if finite_values.size else None,
            'zero_denominator_in_mask':int(np.count_nonzero(eligible & (denominator==0))),
            'negative_denominator_in_mask':int(np.count_nonzero(eligible & (denominator<0)))}


def reproduce():
    reports=[]
    patients=json.loads((OUT/'extraction_summary.json').read_text())['patients']
    for r in patients:
        for visit in ['T0','T1']:
            pre,early,late,mask,voxel_cc,stored,v,maskvalues=load_case(r['patient_id'],visit)
            calculated={kind:volume(pre,early,late,mask,voxel_cc,x['lower'],x['upper']) for kind,x in stored.items()}
            passed=all(calculated[k]['count']==stored[k]['count'] and
                       abs(calculated[k]['volume_cc']-stored[k]['volume_cc'])<1e-5 for k in stored)
            report={'patient_id':r['patient_id'],'visit':visit,'split':r['split'],'status':'pass' if passed else 'failed',
                'stored':stored,'recalculated':calculated,'eligible_mask_volume_cc':float(mask.sum()*voxel_cc),
                'eligible_voxels':int(mask.sum()),'mask_values':maskvalues,'voxel_cc':voxel_cc}
            reports.append(report)
            (OUT/'ftv_reference_reproduction.json').write_text(json.dumps(reports,indent=2),encoding='utf-8')
            print(r['patient_id'],visit,report['status'],'counts',
                  {k:(calculated[k]['count'],stored[k]['count']) for k in stored},flush=True)
    if any(r['status']!='pass' for r in reports):raise SystemExit(1)


def evaluate(diagnostic=False):
    reference=json.loads((OUT/'ftv_reference_reproduction.json').read_text())
    assert len(reference)==18
    if not diagnostic:
        assert all(r['status']=='pass' for r in reference), 'Stored FTV reproduction failed; primary evaluation remains blocked'
    train=[r for r in reference if r['split']=='development']
    validation=[r for r in reference if r['split']=='internal_validation']
    assert len(train)==12 and len(validation)==6
    target_key='recalculated' if diagnostic else 'stored'
    x=np.array([r['eligible_mask_volume_cc'] for r in train]);y=np.array([r[target_key]['ser']['volume_cc'] for r in train])
    coefficient=float(np.clip(x@y/(x@x),0,1))
    direct={'coefficient':coefficient,'training_patients':sorted({r['patient_id'] for r in train}),
        'target':'recalculated reference-ROI SER-volume' if diagnostic else 'stored FTVser',
        'input':'volume of published inverse mask == 0; reference-ROI-conditioned',
        'fitting':'least squares through origin, coefficient restricted to [0,1]',
        'plan_sha256':hashlib.sha256((ROOT/'docs/theme3_ftv_evaluation_plan_v1.md').read_bytes()).hexdigest()}
    directpath=OUT/('direct_roi_ser_diagnostic_frozen.json' if diagnostic else 'direct_ftv_model_frozen_before_validation.json')
    if directpath.exists():assert json.loads(directpath.read_text())==direct
    else:directpath.write_text(json.dumps(direct,indent=2),encoding='utf-8')
    modelpath=OUT/'ridge_frozen_before_validation.json'
    model=json.loads(modelpath.read_text());coeff=np.array(model['coef'])
    assert set(direct['training_patients']).isdisjoint(r['patient_id'] for r in validation)
    assert set(model['training_patients'])==set(direct['training_patients'])
    results=[]
    for r in validation:
        pre,early,late,mask,voxel_cc,stored,v,maskvalues=load_case(r['patient_id'],r['visit'])
        scale=float(np.percentile(pre[pre>0],95))
        prediction=early+scale*(coeff[0]*pre/scale+coeff[1]*(early-pre)/scale+
            coeff[2]*v['delays_seconds'][1]/600+coeff[3]*v['delays_seconds'][2]/600+model['intercept'])
        methods={'copy_early':volume(pre,early,early,mask,voxel_cc,stored['ser']['lower'],stored['ser']['upper']),
                 'ridge_image':volume(pre,early,prediction,mask,voxel_cc,stored['ser']['lower'],stored['ser']['upper']),
                 'direct_volume':{'volume_cc':coefficient*r['eligible_mask_volume_cc']}}
        truth=r[target_key]['ser']['volume_cc']
        for value in methods.values():
            value['error_cc']=value['volume_cc']-truth
            value['absolute_error_cc']=abs(value['error_cc'])
        results.append({'patient_id':r['patient_id'],'visit':r['visit'],'real_reference_volume_cc':truth,'methods':methods,
            'real_reference_diagnostics':volume(pre,early,late,mask,voxel_cc,stored['ser']['lower'],stored['ser']['upper']),
            'stored_FTVser_cc':stored['ser']['volume_cc']})
    patients=[]
    for pid in sorted({r['patient_id'] for r in results}):
        visits={r['visit']:r for r in results if r['patient_id']==pid}
        real0,real1=[visits[v]['real_reference_volume_cc'] for v in ['T0','T1']]
        metrics={}
        for name in ['copy_early','ridge_image','direct_volume']:
            p0,p1=[visits[v]['methods'][name]['volume_cc'] for v in ['T0','T1']]
            metrics[name]={'mean_absolute_volume_error_cc':(abs(p0-real0)+abs(p1-real1))/2,
                'predicted_change_cc':p1-p0,'change_error_cc':(p1-p0)-(real1-real0),
                'real_percent_change':100*(real1-real0)/real0 if real0>0 else None,
                'predicted_percent_change':100*(p1-p0)/p0 if p0>0 else None,
                'false_disappearance':bool(real1>0 and p1==0)}
        patients.append({'patient_id':pid,'real_change_cc':real1-real0,'metrics':metrics})
        print(pid,json.dumps(metrics),flush=True)
    report={'scope':'Exploratory extension on 3 consumed validation patients; fixed reference ROI; not autonomous clinical validation',
        'endpoint':'recalculated reference-ROI SER-volume' if diagnostic else 'stored FTVser',
        'stored_FTV_reproduction_passed':all(r['status']=='pass' for r in reference),
        'model_sha256':hashlib.sha256(modelpath.read_bytes()).hexdigest(),
        'direct_model_sha256':hashlib.sha256(directpath.read_bytes()).hexdigest(),
        'visits':results,'patients':patients}
    name='roi_ser_diagnostic.json' if diagnostic else 'ftv_predictor_evaluation.json'
    (OUT/name).write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--evaluate',action='store_true');p.add_argument('--diagnostic-recomputed',action='store_true');args=p.parse_args()
    if args.diagnostic_recomputed and not args.evaluate:p.error('--diagnostic-recomputed requires --evaluate')
    evaluate(args.diagnostic_recomputed) if args.evaluate else reproduce()
