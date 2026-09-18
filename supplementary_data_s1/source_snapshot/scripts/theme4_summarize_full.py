"""Aggregate only after every planned examination has a terminal record."""
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from theme4_full_cohort import OUT
from theme4_pilot_agreement import wilson, agreement


def csvout(name, rows):
    if not rows:return
    with (OUT/name).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def proportion(k,n,planned):
    return {'k':int(k),'evaluable':n,'planned':planned,'fraction':k/n if n else None,
            'wilson95':wilson(k,n) if n else None,
            'missing_worst_case_bounds':[k/planned,(k+planned-n)/planned]}


def main():
    plan=json.loads((OUT/'locked_plan.json').read_text())
    reports=[json.loads(p.read_text()) for p in sorted(OUT.glob('*_q2.json'))]
    expected={(pid,v) for pid in plan['patients'] for v in ['T0','T1']}
    assert len(reports)==234 and {(r['patient_id'],r['visit']) for r in reports}==expected,'Wait for all 234 terminal records'
    complete=[r for r in reports if r['status']=='complete']
    visits=[];audit=[]
    for r in reports:
        row={'patient_id':r['patient_id'],'visit':r['visit'],'status':r['status'],'error':r.get('error','')}
        for key in ['precision','real_cc','sub_cc','difference_cc','relative_difference','retention_real_over_sub','abs_ge10','over_ge10','under_ge10','zero_reference']:
            row[key]=r.get(key)
        row['excluded_fraction']=1-r['retention_real_over_sub'] if r.get('retention_real_over_sub') is not None else None
        row['changed_fraction_eligible']=None
        if r['status']=='complete':
            # Fixed support, substitution includes all nonzero early increments;
            # real support is a subset when eligible PE is positive. Check rather
            # than silently assuming inclusion for atypical source data.
            row['eligible_voxels']=r['substitution']['eligible_voxels']
            row['removed_count_real_vs_sub']=r['substitution']['count']-r['calculated']['ser']['count']
            row['changed_fraction_eligible']=row['removed_count_real_vs_sub']/row['eligible_voxels'] if row['eligible_voxels'] else None
        else:row['eligible_voxels']=None;row['removed_count_real_vs_sub']=None
        visits.append(row)
        for kind in ['pe','ser']:
            a={'patient_id':r['patient_id'],'visit':r['visit'],'endpoint':kind,'status':r['status'],
               'stored_cc':None,'recomputed_cc':None,'count_difference':None,'difference_cc':None,
               'relative_difference':None,'exam_exact_match':r.get('exact_match'),
               'exam_within_tolerance':r.get('within_tolerance'),'cause':'unresolved_not_attributed'}
            if r['status']=='complete':
                s,c=r['stored'][kind],r['calculated'][kind]
                a.update(stored_cc=s['volume_cc'],recomputed_cc=c['volume_cc'],
                         count_difference=c['count']-s['count'],difference_cc=c['volume_cc']-s['volume_cc'],
                         relative_difference=(c['volume_cc']-s['volume_cc'])/s['volume_cc'] if s['volume_cc'] else None)
            audit.append(a)
    csvout('all_examinations.csv',visits);csvout('Q1_all_endpoints.csv',audit)
    paired=[]
    for pid in sorted(plan['patients']):
        rs={r['visit']:r for r in complete if r['patient_id']==pid}
        if set(rs)!={'T0','T1'}:continue
        r0,r1=[rs[v]['real_cc'] for v in ['T0','T1']];s0,s1=[rs[v]['sub_cc'] for v in ['T0','T1']]
        paired.append({'patient_id':pid,'q3_split':plan['patients'][pid],
                       'real_change_cc':r1-r0,'sub_change_cc':s1-s0,'difference_change_cc':(s1-s0)-(r1-r0),
                       'percentage_change_difference_points':100*((s1-s0)/s0-(r1-r0)/r0) if s0>0 and r0>0 else None,
                       'either_abs_ge10':rs['T0']['abs_ge10'] or rs['T1']['abs_ge10'],
                       'either_over_ge10':rs['T0']['over_ge10'] or rs['T1']['over_ge10']})
    csvout('all_patient_changes.csv',paired)
    summary={'planned_exams':234,'evaluable_exams':len(complete),'paired_patients':len(paired),
             'exact_match_exams':sum(r['exact_match'] for r in complete),
             'tolerance_match_exams':sum(r['within_tolerance'] for r in complete),
             'unresolved':[r for r in reports if r['status']!='complete'],
             'endpoint':'reference-mask-conditioned SER volume; not autonomous early-only FTV',
             'Q3_fitted':False,'by_visit':{},'agreement':{}}
    fig,axes=plt.subplots(1,3,figsize=(15,4))
    for ax,name in zip(axes,['T0','T1','Change']):
        selected=paired if name=='Change' else [r for r in complete if r['visit']==name]
        if name!='Change':
            summary['by_visit'][name]={key:proportion(sum(r[key] for r in selected),len(selected),117) for key in ['abs_ge10','over_ge10','under_ge10']}
        if len(selected)<2:continue
        keys=('real_change_cc','sub_change_cc') if name=='Change' else ('real_cc','sub_cc')
        real,sub=[np.array([r[k] for r in selected]) for k in keys]
        stats=agreement(real,sub);summary['agreement'][name]=stats
        ax.scatter((real+sub)/2,sub-real,s=12,alpha=.65)
        for value in stats['bias_and_limits_cc']:ax.axhline(value,color='grey',ls='--')
        ax.set(title=f'{name}: n={len(selected)} patients',xlabel='Mean (cc)',ylabel='S1 substitution minus measured S2 (cc)')
    summary['either_visit_abs_ge10']=proportion(sum(r['either_abs_ge10'] for r in paired),len(paired),117)
    fig.suptitle('Full cohort: reference-mask-conditioned SER-volume ablation')
    fig.tight_layout();fig.savefig(OUT/'bland_altman.png',dpi=160);plt.close(fig)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
