"""Reproduce the published-score AUC calculations; no image download or fitting."""
import csv
import json
from pathlib import Path
import numpy as np
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score
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
    base=Path(__file__).resolve().parent
    rows=list(csv.DictReader((base/'paired_scores.csv').open(encoding='utf-8-sig')))
    expected=json.loads((base/'auc_results.json').read_text())['results']
    assert len(rows)==117 and len(expected)==20
    verified=[]
    for e in expected:
        group=e['group'];metric=e['metric']
        selected=[r for r in rows if group=='all' or r['subtype']==group]
        result=compare(np.array([int(r['pcr']) for r in selected]),
                       np.array([float(r['real_'+metric]) for r in selected]),
                       np.array([float(r['sub_'+metric]) for r in selected]))
        for key in ['real_auc','sub_auc','difference_sub_minus_real','real_ci','sub_ci','difference_ci']:
            assert np.allclose(result[key],e[key],rtol=0,atol=1e-12),(group,metric,key,result[key],e[key])
        verified.append({'group':group,'metric':metric,'matches_frozen_results':True})
    report={'comparisons':verified,'all_match':True,'absolute_tolerance':1e-12}
    (base/'verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('All 20 AUC comparisons and confidence intervals reproduced within 1e-12.')

if __name__=='__main__':main()
