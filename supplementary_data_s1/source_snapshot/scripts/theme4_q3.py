"""Locked crop-conditioned exploratory retention regression; no image synthesis."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error,mean_squared_error,r2_score
import joblib
from theme4_full_cohort import ROOT,OUT as SOURCE

OUT=ROOT/'data/derived/theme4_q3_v1'


def dump(name,x):
    (OUT/name).write_text(json.dumps(x,indent=2,allow_nan=False),encoding='utf-8')


def features(pid,visit):
    meta=json.loads((SOURCE/f'{pid}_extraction.json').read_text())
    v=meta['visits'][visit];path=ROOT/v['pair_file']
    if meta['manufacturer']=='Philips Medical Systems':
        path=ROOT/json.loads((SOURCE/f'{pid}_{visit}_native_precision.json').read_text())['pair_file']
    with np.load(path) as data:
        pre=data['pre'].astype(float);early=data['early'].astype(float)
    scale=np.percentile(pre[pre>0],95)
    support=pre>.1*scale
    pe=np.divide(100*(early-pre),pre,out=np.full_like(pre,np.nan),where=pre>0)
    selected=support & np.isfinite(pe) & (pe>70)
    delay=(v['delays_seconds'][1]-v['delays_seconds'][0])/60
    assert delay>0
    quant=np.quantile(pe[selected],[.25,.5,.75,.9]) if selected.any() else [np.nan]*4
    return [quant[1],quant[2]-quant[0],quant[3],selected.sum()/support.sum(),quant[1]/delay,delay]


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    assert not (OUT/'results.json').exists(),'Retain completed evaluation; new analysis needs a new version'
    paths=[ROOT/'docs/theme4_q3_plan_v1.md',Path(__file__),SOURCE/'locked_plan.json']
    freeze={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    lock=OUT/'before_fit_undefined_target_handling.json'
    if lock.exists():assert json.loads(lock.read_text())==freeze
    else:dump(lock.name,freeze)
    split=json.loads((SOURCE/'locked_plan.json').read_text())['patients']
    with next((ROOT/'data/raw/acrin6698').glob('*Training*.csv')).open(encoding='utf-8-sig') as f:
        subtype={r['Patient ID DICOM']:r['hrher4g'] for r in csv.DictReader(f)}
    ids=sorted(split);train=np.array([split[p]=='development' for p in ids]);test=~train
    assert train.sum()==88 and test.sum()==29
    output={}
    for visit in ['T1','T0']:
        featurefile=OUT/f'{visit}_features.json'
        if featurefile.exists():rows=json.loads(featurefile.read_text())
        else:
            rows=[]
            for i,pid in enumerate(ids):
                values=features(pid,visit)
                rows.append({'patient_id':pid,'values':[float(x) if np.isfinite(x) else None for x in values],'subtype':subtype[pid]})
                if (i+1)%20==0:print(visit,'features',i+1,flush=True)
            dump(featurefile.name,rows)
        assert [r['patient_id'] for r in rows]==ids
        eligible=np.array([json.loads((SOURCE/f'{p}_{visit}_q2.json').read_text())['sub_cc']>0 for p in ids])
        train=np.array([split[p]=='development' for p in ids]) & eligible
        test=np.array([split[p]=='test' for p in ids]) & eligible
        dump(f'{visit}_eligibility.json',{'excluded_undefined_denominator':np.array(ids)[~eligible].tolist(),
             'development':int(train.sum()),'test':int(test.sum()),'split_changed':False})
        X=np.array([[np.nan if x is None else x for x in r['values']]+[r['subtype']] for r in rows],dtype=object)
        ytrain=np.array([json.loads((SOURCE/f'{p}_{visit}_q2.json').read_text())['retention_real_over_sub'] for p in np.array(ids)[train]],float)
        assert np.isfinite(ytrain).all()
        preds={'constant':np.repeat(ytrain.mean(),test.sum())};fitted={}
        for name,nums,cat in [('subtype',[],True),('early',list(range(6)),False),('combined',list(range(6)),True)]:
            transforms=[]
            if nums:transforms.append(('numeric',make_pipeline(SimpleImputer(strategy='median'),StandardScaler()),nums))
            if cat:transforms.append(('subtype',OneHotEncoder(handle_unknown='ignore',sparse_output=False),[6]))
            model=make_pipeline(ColumnTransformer(transforms),Ridge(alpha=10.))
            model.fit(X[train],ytrain);preds[name]=model.predict(X[test])
            modelpath=OUT/f'{visit}_{name}.joblib';joblib.dump(model,modelpath)
            fitted[name]=hashlib.sha256(modelpath.read_bytes()).hexdigest()
        dump(f'{visit}_predictions_before_test_target_read.json',{'patients':np.array(ids)[test].tolist(),
             'predictions':{k:v.tolist() for k,v in preds.items()},'model_hashes':fitted})
        y=np.array([json.loads((SOURCE/f'{p}_{visit}_q2.json').read_text())['retention_real_over_sub'] for p in np.array(ids)[test]],float)
        assert np.isfinite(y).all()
        rng=np.random.default_rng(20260918);indices=rng.integers(0,len(y),(10000,len(y)))
        results={}
        for name,p in preds.items():
            loss=abs(p-y)
            results[name]={'MAE':mean_absolute_error(y,p),'RMSE':float(np.sqrt(mean_squared_error(y,p))),
                'R2':r2_score(y,p),'out_of_range':int(((p<0)|(p>1)).sum()),'comparisons':{}}
            for baseline in ['constant','subtype']:
                d=loss-abs(preds[baseline]-y)
                results[name]['comparisons'][baseline]={'MAE_difference':float(d.mean()),
                    'bootstrap95':np.quantile(d[indices].mean(axis=1),[.025,.975]).tolist()}
        output[visit]=results;print(visit,json.dumps(results),flush=True)
    dump('results.json',{'scope':'exploratory, prelocalized public crop, exposed Q2 cohort; no independent clinical validation','visits':output})


if __name__=='__main__':main()
