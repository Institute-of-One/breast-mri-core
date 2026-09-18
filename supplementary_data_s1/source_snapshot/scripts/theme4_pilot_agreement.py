"""Descriptive agreement on exposed pilot patients; no population claim."""
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from theme4_audit import OUT


def wilson(k, n):
    z=1.959963984540054
    p=k/n; den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [float(center-half), float(center+half)]


def agreement(real, sub):
    d=sub-real
    def estimates(a):
        mean=a.mean(axis=-1); sd=a.std(axis=-1,ddof=1)
        return np.stack([mean,mean-1.96*sd,mean+1.96*sd],axis=-1)
    rng=np.random.default_rng(20260917)
    boot=estimates(d[rng.integers(0,len(d),size=(10000,len(d)))])
    point=estimates(d)
    return {'n_patients':len(d),'bias_and_limits_cc':point.tolist(),
            'bootstrap_95_intervals_for_bias_lower_limit_upper_limit':np.quantile(boot,[.025,.975],axis=0).T.tolist()}


def main():
    with (OUT/'pilot_exam_audit.csv').open(encoding='utf-8-sig') as f: rows=list(csv.DictReader(f))
    with (OUT/'pilot_longitudinal_audit.csv').open(encoding='utf-8-sig') as f: changes=list(csv.DictReader(f))
    report={'scope':'9 selected exposed pilot patients; fixed reference masks; not full FTV or independent validation',
            'proportions':{},'agreement':{}}
    fig, axes=plt.subplots(1,3,figsize=(13,4))
    for ax,name in zip(axes,['T0','T1','Change']):
        selected=changes if name=='Change' else [r for r in rows if r['visit']==name]
        keys=('real_change_cc','sub_change_cc') if name=='Change' else ('Vreal_mask_conditioned_cc','Vsub_mask_conditioned_cc')
        real,sub=[np.array([float(r[k]) for r in selected]) for k in keys]
        stats=agreement(real,sub); report['agreement'][name]=stats
        ax.scatter((real+sub)/2, sub-real)
        for val in stats['bias_and_limits_cc']:ax.axhline(val,color='grey',linestyle='--')
        ax.set(title=name+' (n=9)',xlabel='Mean of paired volumes / changes (cc)',ylabel='S1 substitution minus real (cc)')
        if name!='Change':
            k=sum(r['at_least_10_percent']=='True' for r in selected)
            report['proportions'][name]={'count':k,'n':len(selected),'wilson_95':wilson(k,len(selected))}
    ids={r['patient_id'] for r in rows}
    k=sum(any(r['at_least_10_percent']=='True' for r in rows if r['patient_id']==pid) for pid in ids)
    report['proportions']['either_visit']={'count':k,'n':len(ids),'wilson_95':wilson(k,len(ids))}
    fig.suptitle('Exploratory pilot: reference-mask-conditioned SER volume')
    fig.tight_layout();fig.savefig(OUT/'pilot_bland_altman.png',dpi=160);plt.close(fig)
    (OUT/'pilot_agreement.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
