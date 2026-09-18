"""Theme 4 stage 1: immutable source audit, not source-pipeline reconstruction."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pydicom

from theme3_ftv import load_case, volume
from theme3_start import ROOT
from theme3_lock_development import OUT as SOURCE

OUT = ROOT / 'data/derived/theme4_v1'


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    protocol = ROOT / 'docs/theme4_protocol_v1.md'
    lock = {'protocol_sha256': hashlib.sha256(protocol.read_bytes()).hexdigest(),
            'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'source_report_sha256': hashlib.sha256((SOURCE / 'ftv_reference_reproduction.json').read_bytes()).hexdigest(),
            'scope': '18 exposed pilot exams; inventory of 234 cached headers; no Q3'}
    lockpath = OUT / 'stage1_lock.json'
    if lockpath.exists():
        assert json.loads(lockpath.read_text()) == lock, 'Changed analysis requires a new version'
    else:
        lockpath.write_text(json.dumps(lock, indent=2), encoding='utf-8')
    inventory = []
    for path in sorted((ROOT / 'data/derived/training_cohort/patients').glob('*.json')):
        patient = json.loads(path.read_text())
        for visit in ['T0', 'T1']:
            row = {'patient_id': path.stem, 'visit': visit}
            try:
                d = pydicom.dcmread(ROOT / patient['visits'][visit]['inputs']['ser'], stop_before_pixels=True)
                params = {str(x[(0x0117,0x1014)].value): x for x in d[(0x0117,0x1010)].value}
                for key in ['PE_threshold', 'pre_contrast_threshold', 'minimum_neighbor_count']:
                    row[key] = float(params[key][(0x0117,0x1019)].value)
                row['status'] = 'header_available'
                row['error'] = ''
            except Exception as exc:
                for key in ['PE_threshold', 'pre_contrast_threshold', 'minimum_neighbor_count']:
                    row.setdefault(key, None)
                row['status'], row['error'] = 'unresolved', str(exc)
            inventory.append(row)
    write_csv(OUT / 'reference_header_inventory.csv', inventory)
    references = json.loads((SOURCE / 'ftv_reference_reproduction.json').read_text())
    rows = []
    for ref in references:
        pid, visit = ref['patient_id'], ref['visit']
        pre, early, late, mask, cc, stored, _, bits = load_case(pid, visit)
        threshold, upper = stored['ser']['lower'], stored['ser']['upper']
        real = volume(pre, early, late, mask, cc, threshold, upper)
        sub = volume(pre, early, early, mask, cc, threshold, upper)
        ratio = np.divide(early-pre, late-pre, out=np.full_like(pre, np.nan), where=late!=pre)
        clipped = np.clip(ratio, 0, upper)
        real_support = mask & np.isfinite(clipped) & (clipped > threshold)
        sub_support = mask & (early != pre)  # SER=1, with validated threshold below.
        assert 0 < threshold < 1 < upper
        assert int(sub_support.sum()) == sub['count']
        pe_error = ref['recalculated']['pe']['count'] - stored['pe']['count']
        ser_error = real['count'] - stored['ser']['count']
        inclusive = int(np.count_nonzero(mask & np.isfinite(clipped) & (clipped >= threshold)))
        rounded = int(np.count_nonzero(mask & np.isfinite(clipped) & (clipped > .9)))
        vreal, vsub = real['volume_cc'], sub['volume_cc']
        rows.append({'patient_id': pid, 'visit': visit, 'original_exact_match': ref['status'],
                     'pe_count_difference': pe_error, 'ser_count_difference': ser_error,
                     'ser_cc_difference': ser_error*cc,
                     'inclusive_threshold_count_change': inclusive-real['count'],
                     'rounded_0_9_count_change': rounded-real['count'],
                     'voxels_within_0_001_of_threshold': int(np.count_nonzero(mask & (abs(clipped-threshold)<=.001))),
                     'undocumented_bit16_present': any(b & 16 for b in bits),
                     'geometry': 'physical_grid_matched_without_resampling',
                     'VOI_connectivity_intensity_cause': 'not_identified' if pe_error else 'not_established_by_count_agreement',
                     'threshold_cause': 'candidate_only' if ser_error else 'no_count_discrepancy',
                     'Vreal_mask_conditioned_cc': vreal, 'Vsub_mask_conditioned_cc': vsub,
                     'sub_minus_real_cc': vsub-vreal,
                     'relative_difference_real_denominator': (vsub-vreal)/vreal if vreal else None,
                     'at_least_10_percent': bool(abs(vsub-vreal)/vreal>=.1) if vreal else bool(vsub>0),
                     'discordant_zero_reference': bool(vreal==0 and vsub>0),
                     'changed_voxel_fraction_eligible': float(np.count_nonzero(real_support ^ sub_support)/mask.sum()) if mask.any() else None})
        print(pid, visit, pe_error, ser_error, flush=True)
    write_csv(OUT / 'pilot_exam_audit.csv', rows)
    changes=[]
    for pid in sorted({r['patient_id'] for r in rows}):
        visits = {r['visit']: r for r in rows if r['patient_id']==pid}
        r0,r1 = [visits[v]['Vreal_mask_conditioned_cc'] for v in ['T0','T1']]
        s0,s1 = [visits[v]['Vsub_mask_conditioned_cc'] for v in ['T0','T1']]
        changes.append({'patient_id':pid,'real_change_cc':r1-r0,'sub_change_cc':s1-s0,
                        'change_difference_cc':(s1-s0)-(r1-r0),
                        'percentage_change_difference_points':100*((s1-s0)/s0-(r1-r0)/r0) if r0>0 and s0>0 else None})
    write_csv(OUT / 'pilot_longitudinal_audit.csv', changes)
    summary = {'inventory_exams':len(inventory), 'header_available':sum(r['status']=='header_available' for r in inventory),
               'image_mask_exams':len(rows), 'exact_matches':sum(r['original_exact_match']=='pass' for r in rows),
               'nonzero_PE_count_discrepancies':sum(r['pe_count_difference']!=0 for r in rows),
               'nonzero_SER_count_discrepancies':sum(r['ser_count_difference']!=0 for r in rows),
               'causal_decomposition_complete':False, 'full_cohort_Q2_complete':False,
               'warning':'Pilot only. Fixed source masks. Full FTV and late-independent support not established.'}
    (OUT / 'stage1_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
