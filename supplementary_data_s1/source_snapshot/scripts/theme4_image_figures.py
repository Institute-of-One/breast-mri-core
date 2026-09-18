"""Real-image manuscript figures; no fitting or revision of frozen estimates."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

import theme3_ftv as ftv
from theme3_start import ROOT

SOURCE = ROOT / 'data/derived/theme4_full_v2'
PCR = ROOT / 'data/derived/theme4_pcr_v1'
OUT = ROOT / 'data/derived/theme4_figures_v1'
CYAN = '#30c7d5'
ORANGE = '#f1a340'
BLUE = '#276b96'
PURPLE = '#8b5b9a'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(fig, name):
    for ext in ['png', 'svg']:
        fig.savefig(OUT / f'{name}.{ext}', dpi=250, facecolor='white')
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.fonttype': 'none'})
    rows = [json.loads(p.read_text()) for p in sorted(SOURCE.glob('*_q2.json'))]
    assert len(rows) == 234 and all(r['status'] == 'complete' for r in rows)
    by = {(r['patient_id'], r['visit']): r for r in rows}
    ids = sorted({r['patient_id'] for r in rows})
    assert len(ids) == 117
    medians = {v: float(np.median([r['retention_real_over_sub'] for r in rows
                    if r['visit'] == v and r['retention_real_over_sub'] is not None]))
               for v in ['T0', 'T1']}
    candidates = [p for p in ids if all(by[p, v]['retention_real_over_sub'] is not None
                                       for v in ['T0', 'T1'])]
    distance = lambda p: sum((by[p, v]['retention_real_over_sub']-medians[v])**2
                             for v in ['T0', 'T1'])
    pid = min(candidates, key=lambda p: (distance(p), p))
    ftv.OUT = SOURCE
    cases = []
    checks = []
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in sorted(SOURCE.glob('*_q2.json'))}
    source_list = json.loads((SOURCE/'ftv_mask_sources.json').read_text())
    extraction = json.loads((SOURCE/f'{pid}_extraction.json').read_text())
    for visit in ['T0', 'T1']:
        pre, early, late, eligible, cc, stored, meta, _ = ftv.load_case(pid, visit)
        ratio = np.divide(early-pre, late-pre, out=np.full_like(pre, np.nan), where=late!=pre)
        real = eligible & np.isfinite(ratio) & (np.clip(ratio, 0, stored['ser']['upper']) > stored['ser']['lower'])
        sub = eligible & (early != pre)
        r = by[pid, visit]
        assert int(real.sum()) == r['calculated']['ser']['count']
        assert int(sub.sum()) == r['substitution']['count']
        assert np.isclose(real.sum()*cc, r['real_cc'])
        assert np.isclose(sub.sum()*cc, r['sub_cc'])
        z = int(np.argmax(eligible.sum(axis=(1, 2))))
        pair = ROOT / meta['pair_file']
        if extraction['manufacturer'] == 'Philips Medical Systems':
            correction = json.loads((SOURCE/f'{pid}_{visit}_native_precision.json').read_text())
            pair = ROOT/correction['pair_file']
        with np.load(pair) as a:
            spacing = a['pixel_spacing'].copy()
        maskpath = ROOT/next(x['path'] for x in source_list if x['patient_id']==pid and x['visit']==visit)
        for p in [pair, maskpath, SOURCE/f'{pid}_extraction.json']:
            hashes[str(p.relative_to(ROOT))] = sha(p)
        cases.append((pre, early, late, real, sub, z, spacing, meta, r))
        checks.append({'visit': visit, 'slice_index_zero_based': z,
                       'real_voxels': int(real.sum()), 'sub_voxels': int(sub.sum()),
                       'real_only_voxels': int((real & ~sub).sum()),
                       'real_cc': r['real_cc'], 'sub_cc': r['sub_cc'],
                       'pixel_spacing_mm': spacing.tolist(),
                       'early_delay_seconds': meta['delays_seconds'][1],
                       'late_delay_seconds': meta['delays_seconds'][2]})
        assert not np.any(real & ~sub), 'Unexpected real-only voxels: revise legend'

    # Common physical field width, centered on the analysis support in each visit.
    field_mm = max(max(c[0].shape[1]*c[6][0], c[0].shape[2]*c[6][1]) for c in cases)
    fig, axes = plt.subplots(2, 4, figsize=(14, 8.2))
    fig.subplots_adjust(left=.035, right=.995, top=.86, bottom=.16, hspace=.26, wspace=.035)
    fig.suptitle('Late-phase removal changes the measured tumor volume', x=.04, ha='left', y=.98,
                 fontsize=20, fontweight='bold')
    fig.text(.04, .93, 'Actual MRI and computed masks | representative case selected without pCR labels',
             fontsize=12, color='#52616b')
    titles = ['Precontrast  S0', 'Early phase  S1', 'Measured late phase  S2', 'Volume membership on S1']
    for row, c in enumerate(cases):
        pre, early, late, real, sub, z, spacing, meta, r = c
        planes = [pre[z], early[z], late[z], early[z]]
        lo = 0
        hi = np.percentile(np.concatenate([a.ravel() for a in planes[:3]]), 99.5)
        assert hi > lo
        h, w = early[z].shape
        ext = (0, w*spacing[1], h*spacing[0], 0)
        yy, xx = np.where(sub[z])
        cx, cy = (xx.mean()+.5)*spacing[1], (yy.mean()+.5)*spacing[0]
        for col, ax in enumerate(axes[row]):
            ax.set_facecolor('black')
            ax.imshow(planes[col], cmap='gray', vmin=lo, vmax=hi, extent=ext, interpolation='nearest')
            ax.set_xlim(cx-field_mm/2, cx+field_mm/2)
            ax.set_ylim(cy+field_mm/2, cy-field_mm/2)
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values(): spine.set_visible(False)
            label = chr(65+row*4+col)
            ax.set_title(f'{label}  {titles[col]}', fontsize=11, loc='left', pad=8)
            if col in [1, 2]:
                ax.text(.97,.04,f"{meta['delays_seconds'][col]:.0f} s", transform=ax.transAxes,
                        color='white', ha='right', fontsize=10)
            if col == 0:
                ax.text(.04,.94,['T0 | baseline', 'T1 | early treatment'][row],
                        transform=ax.transAxes, color='white', va='top', fontsize=11, weight='bold')
                x0, y0 = cx-field_mm*.42, cy+field_mm*.39
                ax.plot([x0, x0+20],[y0,y0],color='white',lw=3)
                ax.text(x0, y0-field_mm*.035, '20 mm', color='white', fontsize=9)
            if col == 3:
                rgba = np.zeros((*real[z].shape,4))
                rgba[real[z]] = matplotlib.colors.to_rgba(CYAN, .72)
                rgba[sub[z] & ~real[z]] = matplotlib.colors.to_rgba(ORANGE, .72)
                ax.imshow(rgba, extent=ext, interpolation='nearest')
                ax.text(.04,.04,f"3D volume: {r['real_cc']:.2f} -> {r['sub_cc']:.2f} mL",
                        transform=ax.transAxes, color='white', fontsize=10,
                        bbox={'facecolor':'black','alpha':.72,'edgecolor':'none','pad':3})
    fig.legend(handles=[Patch(color=CYAN,label='Included with measured S2 and S1 substitution'),
                        Patch(color=ORANGE,label='Additional volume with S2 := S1')],
               loc='lower left',bbox_to_anchor=(.035,.073),ncol=2,frameon=False,fontsize=11)
    fig.text(.04,.045,'Same published analysis support in both calculations. One source slice per visit; volumes use all slices.',fontsize=10)
    fig.text(.04,.018,'Shared grayscale window within each visit; no cross-visit registration. Colors show calculation membership, not pathology.',fontsize=10,color='#52616b')
    save(fig,'figure1_actual_mri')

    fig, axs = plt.subplots(1,3,figsize=(16,5.6),gridspec_kw={'width_ratios':[1,1,1.4]})
    fig.subplots_adjust(left=.06,right=.975,top=.77,bottom=.21,wspace=.39)
    fig.suptitle('Volume disagreement and pCR discrimination are different endpoints',x=.06,ha='left',y=.97,
                 fontsize=19,fontweight='bold')
    fig.text(.06,.89,'117 patients | 234 examinations | fixed published masks | paired estimates',fontsize=12,color='#52616b')
    ax=axs[0]
    for v,color in [('T0',BLUE),('T1',ORANGE)]:
        rr=[by[p,v] for p in ids]
        ax.scatter([r['real_cc'] for r in rr],[r['sub_cc'] for r in rr],s=22,color=color,alpha=.7,label=v,edgecolor='white',linewidth=.25)
    maximum=max(r['sub_cc'] for r in rows)*1.06
    ax.plot([0,maximum],[0,maximum],ls='--',color='#7b8790',lw=1)
    ax.set(xlim=(0,maximum),ylim=(0,maximum),xlabel='Measured-S2 volume (mL)',ylabel='S1-substitution volume (mL)')
    ax.set_title('A  Volume comparison',loc='left',fontweight='bold',pad=14)
    ax.legend(frameon=False);ax.grid(alpha=.15)
    ax=axs[1]
    real_relative=np.array([(by[p,'T0']['real_cc']-by[p,'T1']['real_cc'])/by[p,'T0']['real_cc'] for p in ids])
    sub_relative=np.array([(by[p,'T0']['sub_cc']-by[p,'T1']['sub_cc'])/by[p,'T0']['sub_cc'] for p in ids])
    scores=list(csv.DictReader((PCR/'paired_scores.csv').open(encoding='utf-8-sig')))
    scores={r['patient_id']:r for r in scores}
    for i,p in enumerate(ids):
        assert np.isclose(real_relative[i],float(scores[p]['real_relative']))
        assert np.isclose(sub_relative[i],float(scores[p]['sub_relative']))
    low=min(real_relative.min(),sub_relative.min())*100
    low=np.floor(low/20)*20-5; high=105
    ax.scatter(real_relative*100,sub_relative*100,s=25,c=PURPLE,alpha=.75,edgecolor='white',linewidth=.3)
    ax.plot([low,high],[low,high],ls='--',color='#7b8790',lw=1)
    ax.set(xlim=(low,high),ylim=(low,high),xlabel='Measured-S2 reduction (%)',ylabel='S1-substitution reduction (%)')
    ax.set_title('B  Relative reduction: T0 to T1',loc='left',fontweight='bold',pad=14)
    ax.grid(alpha=.15)
    ax=axs[2]
    auc=json.loads((PCR/'results.json').read_text())['results']
    selected=[next(r for r in auc if r['group']=='all' and r['metric']==m) for m in ['T0','relative','T1']]
    labels=['T0 volume','Relative reduction','T1 volume']
    for y,label,r in zip([2,1,0],labels,selected):
        d=r['difference_sub_minus_real'];l,h=r['difference_ci']
        ax.plot([l,h],[y,y],color=BLUE,lw=2.4)
        ax.scatter(d,y,s=45,c=BLUE,zorder=3)
        ax.text(-.078,y+.37,label,fontsize=11,fontweight='bold')
        ax.text(.083,y,f'{d:+.3f}\n[{l:+.3f}, {h:+.3f}]',va='center',fontsize=10)
    ax.axvline(0,ls='--',lw=1,color='#7b8790')
    ax.set(xlim=(-.08,.17),ylim=(-.55,2.7),yticks=[],xlabel='AUC difference: substitution minus measured')
    ax.set_xticks([-.05,0,.05]);ax.grid(axis='x',alpha=.15)
    ax.set_title('C  Paired AUC differences',loc='left',fontweight='bold',pad=14)
    ax.spines['left'].set_visible(False)
    fig.text(.06,.08,'Dashed lines: equality (A, B) or zero AUC difference (C). Negative reduction denotes volume increase.',fontsize=10)
    fig.text(.06,.035,'Intervals: pointwise 95% paired patient-bootstrap intervals (10,000 resamples). No prespecified noninferiority margin.',fontsize=10,color='#52616b')
    save(fig,'figure2_cohort_results')
    hashes[str((PCR/'results.json').relative_to(ROOT))]=sha(PCR/'results.json')
    hashes[str((PCR/'paired_scores.csv').relative_to(ROOT))]=sha(PCR/'paired_scores.csv')
    provenance={'selected_patient':pid,'selection':'Minimum squared distance to T0 and T1 cohort median real/sub retention; positive substitution volumes at both visits; ties by ID; no pCR labels used',
        'selection_candidates':len(candidates),'cohort_retention_medians':medians,
        'selected_distance_squared':distance(pid),'slice_selection':'Largest eligible-mask cross-sectional area at each visit independently',
        'display':'Source-plane orientation; no anatomical side labels, resampling, smoothing, or registration. One common physical field width; shared S0/S1/S2 window within each visit (0 to pooled 99.5th percentile).',
        'checks':checks,'all_234_complete':True,'all_117_relative_scores_match_frozen_inputs':True,
        'input_sha256':hashes,'script_sha256':sha(Path(__file__))}
    (OUT/'provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    caption=f'''# Theme 4: actual-image and cohort figures

## Figure 1. Effect of late-phase substitution on SER-filtered tumor-volume membership

Precontrast (S0), early postcontrast (S1), and measured late postcontrast (S2) images are shown at baseline (T0; A-D) and early treatment (T1; E-H). The overlay on S1 distinguishes voxels retained with both measured S2 and S2:=S1 (cyan) from additional voxels included after substitution (orange). The public analysis support is held fixed. Displayed volumes are three-dimensional recalculated SER-filtered volumes in mL; the overlays depict only the displayed slice. They are not histologic tumor boundaries. Grayscale windowing is shared across phases within each visit. A common physical field width is used, with 20-mm scale bars. No cross-visit registration was performed; the displayed plane is selected separately at each visit by maximum eligible-mask area. The collection supplies tumor-side cropped images; this is not a validation of autonomous early-only localization.

The illustrative patient ({pid}) was selected without pCR labels by minimum squared distance to the cohort median real/substitution volume ratio at T0 and T1, among {len(candidates)} patients with defined ratios at both visits. Median ratios: T0={medians['T0']:.4f}, T1={medians['T1']:.4f}. This is representative by this specified volume-ratio criterion, not by every clinical characteristic. Full data provenance, phase delays, slice indices, voxel counts, and hashes are in provenance.json.

## Figure 2. Cohort-level volume disagreement and paired pCR discrimination

(A) Measured-late-phase and S1-substitution volumes at T0 and T1 (234 examinations in 117 patients). (B) Relative volume reduction, 100*(T0-T1)/T0, calculated separately for each method (117 patients); negative values denote an increase. Every observation is retained, including zero T1 volumes. Dashed diagonals indicate equality. (C) Differences in pCR discrimination AUC (substitution minus measured late), with pointwise 95% confidence intervals from the frozen 10,000 paired outcome-stratified patient-bootstrap resamples. Negative values favor measured S2. The plotted intervals are not noninferiority tests or guarantees of future loss. No models, thresholds, or scores were retuned for these figures. Main-endpoint prioritization was exploratory and followed inspection of the AUC results, as documented in the manuscript.

## 日本語での読み方

Figure 1は実MRIと計算した領域の重ね合わせ。水色は両方法で含まれる領域、橙色はS1置換で追加される領域です。病理学的な正解領域を示すものではありません。T0/T1の各画像に記載した体積は全スライスから計算した値です。

Figure 2は、体積・相対縮小率の一致とpCR識別AUCが別の評価軸であることを示します。AUC差の区間が0を含むことだけで同等性・非劣性を宣言しません。

## Files

- figure1_actual_mri.png / .svg
- figure2_cohort_results.png / .svg
- provenance.json
- Rebuild: .venv\\Scripts\\python.exe -X utf8 scripts/theme4_image_figures.py

Source images: TCIA ACRIN-6698; retain the dataset acknowledgement and citation required by the source collection in any submitted manuscript. This figure generation does not alter source images or frozen analyses.
'''
    (OUT/'captions.md').write_text(caption,encoding='utf-8')
    print(json.dumps({'selected_patient':pid,'checks':checks,'output':str(OUT)},indent=2))


if __name__ == '__main__':
    main()
