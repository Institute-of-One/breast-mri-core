"""Frozen-direction paired AUC bootstrap for an explicitly post-hoc endpoint pivot."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score
from theme4_full_cohort import ROOT,OUT as SOURCE

OUT=ROOT/'data/derived/theme4_pcr_v1'


def auc_boot(score,pos,neg):
    values=np.concatenate([score[pos],score[neg]],axis=1)
    ranks=rankdata(values,axis=1,method='average')
    npos=pos.shape[1];nneg=neg.shape[1]
    return (ranks[:,:npos].sum(axis=1)-npos*(npos+1)/2)/(npos*nneg)


def compare(y,real,sub):
    positives=np.flatnonzero(y==1);negatives=np.flatnonzero(y==0)
    base={'n':len(y),'pcr':len(positives),'non_pcr':len(negatives),
          'small_class_warning':min(len(positives),len(negatives))<10}
    if not len(positives) or not len(negatives):return dict(base,status='not_estimable')
    ar=roc_auc_score(y,real);ass=roc_auc_score(y,sub)
    # Independent pairwise tie-credit check of the observed AUCs.
    for s,a in [(real,ar),(sub,ass)]:
        d=s[positives,None]-s[negatives]
        assert np.isclose(np.mean((d>0)+.5*(d==0)),a)
    rng=np.random.default_rng(20260918)
    pos=rng.choice(positives,(10000,len(positives)),replace=True)
    neg=rng.choice(negatives,(10000,len(negatives)),replace=True)
    br,bs=auc_boot(real,pos,neg),auc_boot(sub,pos,neg)
    # Check the vectorized bootstrap against sklearn on the first replicates.
    for i in range(3):
        ix=np.r_[pos[i],neg[i]]
        assert np.isclose(br[i],roc_auc_score(y[ix],real[ix]))
        assert np.isclose(bs[i],roc_auc_score(y[ix],sub[ix]))
    return dict(base,status='estimated',real_auc=float(ar),sub_auc=float(ass),
                real_ci=np.quantile(br,[.025,.975]).tolist(),sub_ci=np.quantile(bs,[.025,.975]).tolist(),
                difference_sub_minus_real=float(ass-ar),difference_ci=np.quantile(bs-br,[.025,.975]).tolist())


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    assert not (OUT/'results.json').exists(),'Keep results immutable; use a new version'
    label=next((ROOT/'data/raw/acrin6698').glob('BMMR2-Training*.csv'))
    sources=[Path(__file__),ROOT/'docs/theme4_pcr_plan_v1.md',label,*sorted(SOURCE.glob('*_q2.json'))]
    lock={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    (OUT/'before_auc.json').write_text(json.dumps(lock,indent=2),encoding='utf-8')
    with label.open(encoding='utf-8-sig') as f:clinical=list(csv.DictReader(f))
    assert len(clinical)==len({r['Patient ID DICOM'] for r in clinical})==117
    assert {r['Split'] for r in clinical}=={'Train'}
    assert {r['pcr'] for r in clinical}=={'pCR','Non-pCR'}
    rows=[]
    for c in sorted(clinical,key=lambda r:r['Patient ID DICOM']):
        pid=c['Patient ID DICOM'];v=[json.loads((SOURCE/f'{pid}_{t}_q2.json').read_text()) for t in ['T0','T1']]
        assert all(r['patient_id']==pid and r['status']=='complete' for r in v)
        r0,r1=[r['real_cc'] for r in v];s0,s1=[r['sub_cc'] for r in v]
        rows.append({'patient_id':pid,'pcr':int(c['pcr']=='pCR'),'subtype':c['hrher4g'],
                     'real_T0':-r0,'sub_T0':-s0,'real_T1':-r1,'sub_T1':-s1,
                     'real_change':r0-r1,'sub_change':s0-s1,
                     'real_relative':(r0-r1)/r0 if r0>0 else None,
                     'sub_relative':(s0-s1)/s0 if s0>0 else None})
    with (OUT/'paired_scores.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    results=[]
    for group in ['all']+sorted({r['subtype'] for r in rows}):
        cohort=[r for r in rows if group=='all' or r['subtype']==group]
        for metric in ['T0','T1','change','relative']:
            rs=[r for r in cohort if r['real_'+metric] is not None and r['sub_'+metric] is not None]
            res=compare(np.array([r['pcr'] for r in rs]),np.array([r['real_'+metric] for r in rs]),np.array([r['sub_'+metric] for r in rs]))
            res.update(group=group,metric=metric,excluded=[r['patient_id'] for r in cohort if r not in rs])
            results.append(res);print(json.dumps(res),flush=True)
    report={'scope':'Post-hoc single-biomarker discrimination, fixed public mask, previously used training117 cohort',
            'difference_direction':'substitution minus measured late','primary':'all/change','bootstrap':10000,'results':results}
    (OUT/'results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')


if __name__=='__main__':main()
