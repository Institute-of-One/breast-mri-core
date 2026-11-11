"""Resumable 234-exam public-data acquisition and fixed-mask S2 ablation."""
import argparse
import hashlib
import json
import shutil
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
import numpy as np
import pydicom
import theme3_extract_development as extraction
import theme3_restore_philips as native
import theme3_ftv as ftv
from theme3_start import ROOT, BASE, get

OLD=ROOT/'data/derived/theme3_development_v1'
OUT=ROOT/'data/derived/theme4_full_v2'
RAW=ROOT/'data/raw/theme4_full_v2'
SEED='theme4-q3-20260917-v2'


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024**2),b''):h.update(chunk)
    return h.hexdigest()


def lock():
    OUT.mkdir(parents=True,exist_ok=True);RAW.mkdir(parents=True,exist_ok=True)
    records=json.loads((OLD/'locked_plan.json').read_text())['records']
    pilots={r['patient_id'] for r in records if r['pilot']}
    groups=defaultdict(list)
    for r in records:groups[r['manufacturer']].append(r['patient_id'])
    split={pid:'development' for pid in pilots}
    for group in groups.values():
        candidates=sorted(set(group)-pilots,key=lambda pid:hashlib.sha256((SEED+pid).encode()).hexdigest())
        n=max(1,round(len(group)*.25))
        for i,pid in enumerate(candidates):split[pid]='test' if i<n else 'development'
    plan={'seed':SEED,'patients':split,'exposed_pilots':sorted(pilots),
          'gates_sha256':digest(ROOT/'docs/theme4_gates_v2.md'),
          'script_sha256':digest(__import__('pathlib').Path(__file__)),
          'Q3_status':'not_fitted','target':'T1 real/sub retention; T0 secondary'}
    path=OUT/'locked_plan.json'
    if path.exists():assert json.loads(path.read_text())==plan
    else:path.write_text(json.dumps(plan,indent=2),encoding='utf-8')
    sources=json.loads((OLD/'ftv_mask_sources.json').read_text())
    previous={(r['patient_id'],r['visit']):r for r in sources}
    masks=[]
    for r in records:
        rows=json.loads((ROOT/f"data/derived/training_cohort/series/{r['patient_id']}.json").read_text())
        for visit in ['T0','T1']:
            candidates=[x for x in rows if x['StudyDesc'].endswith('_'+visit) and x.get('SeriesDescription','').endswith('Analysis Mask')]
            assert len(candidates)==1
            old=previous.get((r['patient_id'],visit))
            masks.append({'patient_id':r['patient_id'],'visit':visit,'series_uid':candidates[0]['SeriesInstanceUID'],
                          'path':old['path'] if old else str((RAW/f"{r['patient_id']}_{visit}_mask.dcm").relative_to(ROOT))})
    (OUT/'ftv_mask_sources.json').write_text(json.dumps(masks,indent=2),encoding='utf-8')
    return records,masks


def download(path,endpoint,params,iszip=False):
    if path.exists():return path
    for attempt in range(3):
        partial=path.with_suffix('.part')
        try:
            assert shutil.disk_usage(ROOT).free>10*1024**3, 'Less than 10 GiB free'
            with get(BASE+endpoint,params=params,stream=True) as r,partial.open('wb') as f:
                size=0
                for chunk in r.iter_content(1024**2):
                    size+=len(chunk)
                    if size>2*1024**3:raise RuntimeError('Single-series safety bound 2 GiB')
                    f.write(chunk)
            if iszip:
                with zipfile.ZipFile(partial) as z:assert z.testzip() is None
            else:pydicom.dcmread(partial,stop_before_pixels=True)
            partial.replace(path);return path
        except Exception:
            if partial.exists():partial.unlink()
            if attempt==2:raise
            time.sleep(2)


def tolerance(calc,ref):
    return calc==0 if ref==0 else abs(calc-ref)<=.01 and abs(calc-ref)/abs(ref)<=.001


def process(record,masks):
    pid=record['patient_id'];path=OUT/f'{pid}_extraction.json'
    oldpath=OLD/f'{pid}_extraction.json'
    result=json.loads(path.read_text()) if path.exists() else (json.loads(oldpath.read_text()) if oldpath.exists() else
           {'patient_id':pid,'manufacturer':record['manufacturer'],'visits':{}})
    rows=json.loads((ROOT/f'data/derived/training_cohort/series/{pid}.json').read_text())
    for visit in ['T0','T1']:
        output=OUT/f'{pid}_{visit}_q2.json'
        if output.exists() and json.loads(output.read_text()).get('status')=='complete':continue
        try:
            if result['visits'].get(visit,{}).get('status')!='pass':
                row=record['visits'][visit]['series']
                source=download(RAW/f'{pid}_{visit}.zip','getImage',{'SeriesInstanceUID':row['SeriesInstanceUID']},True)
                result['visits'][visit]=extraction.prepare(source,record,visit)
            path.write_text(json.dumps(result,indent=2),encoding='utf-8')
            precision='published_crop'
            if record['manufacturer']=='Philips Medical Systems':
                correction=OUT/f'{pid}_{visit}_native_precision.json'
                old=OLD/correction.name
                if old.exists() and not correction.exists():shutil.copyfile(old,correction)
                if not correction.exists():
                    candidates=[r for r in rows if r['StudyDesc'].endswith('_'+visit) and
                                int(r['ImageCount'])==int(record['visits'][visit]['series']['ImageCount']) and
                                'VOLSER' not in r.get('SeriesDescription','') and r.get('Modality')=='MR']
                    assert len(candidates)==1, 'Native series selection not unique'
                    row=candidates[0]
                    source=download(RAW/f'{pid}_{visit}_native.zip','getImage',{'SeriesInstanceUID':row['SeriesInstanceUID']},True)
                    corrected=native.prepare(pid,visit,row,source)
                    correction.write_text(json.dumps(corrected,indent=2),encoding='utf-8')
                precision='native_philips_exact_grid'
            source=next(r for r in masks if r['patient_id']==pid and r['visit']==visit)
            maskpath=ROOT/source['path']
            if not maskpath.exists():
                sops=get(BASE+'getSOPInstanceUIDs',params={'SeriesInstanceUID':source['series_uid']}).json()
                assert len(sops)==1
                download(maskpath,'getSingleImage',{'SeriesInstanceUID':source['series_uid'],'SOPInstanceUID':sops[0]['SOPInstanceUID']})
            pre,early,late,mask,cc,stored,_,bits=ftv.load_case(pid,visit)
            calculated={k:ftv.volume(pre,early,late,mask,cc,x['lower'],x['upper']) for k,x in stored.items()}
            sub=ftv.volume(pre,early,early,mask,cc,stored['ser']['lower'],stored['ser']['upper'])
            real=calculated['ser'];vr=real['volume_cc'];vs=sub['volume_cc']
            assert 0<stored['ser']['lower']<1<stored['ser']['upper']
            # Under this fixed support S1 substitution must include every nonzero early increment.
            assert sub['count']==int(np.count_nonzero(mask & (early!=pre)))
            report={'patient_id':pid,'visit':visit,'status':'complete','precision':precision,
                    'endpoint':'reference-mask-conditioned SER volume','mask_sha256':digest(maskpath),
                    'stored':stored,'calculated':calculated,'substitution':sub,
                    'exact_match':all(calculated[k]['count']==stored[k]['count'] and abs(calculated[k]['volume_cc']-stored[k]['volume_cc'])<1e-5 for k in stored),
                    'within_tolerance':all(tolerance(calculated[k]['volume_cc'],stored[k]['volume_cc']) for k in stored),
                    'real_cc':vr,'sub_cc':vs,'difference_cc':vs-vr,
                    'retention_real_over_sub':vr/vs if vs else None,
                    'relative_difference':(vs-vr)/vr if vr else None,
                    'abs_ge10':abs(vs-vr)/vr>=.1 if vr else vs>0,
                    'over_ge10':(vs-vr)/vr>=.1 if vr else vs>0,
                    'under_ge10':(vs-vr)/vr<=-.1 if vr else False,
                    'unknown_bit16':any(b&16 for b in bits),
                    'zero_reference':vr==0}
        except Exception as exc:
            report={'patient_id':pid,'visit':visit,'status':'unresolved','error':repr(exc)}
        output.write_text(json.dumps(report,indent=2),encoding='utf-8')
        print(pid,visit,report['status'],report.get('error',''),flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run',action='store_true');args=parser.parse_args()
    records,masks=lock()
    print('Locked patients',len(records),'exams',len(masks),flush=True)
    if not args.run:return
    extraction.OUT=OUT;native.OUT=OUT;ftv.OUT=OUT
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(process,r,masks) for r in records]
        for f in as_completed(futures):f.result()
    reports=[json.loads(p.read_text()) for p in sorted(OUT.glob('*_q2.json'))]
    assert len(reports)==234
    summary={'planned':234,'complete':sum(r['status']=='complete' for r in reports),
             'unresolved':[r for r in reports if r['status']!='complete'],
             'scope':'all planned exams attempted; conditional-mask endpoint; Q3 not fitted'}
    (OUT/'flow.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
