"""Freeze an outcome-free, manufacturer-stratified development split and pilot."""
import hashlib
import json
from collections import Counter, defaultdict
import pydicom
from theme3_start import ROOT, PID

OUT = ROOT / 'data/derived/theme3_development_v1'
SEED = 'theme3-development-20260917-v1'


def rank(pid):
    return hashlib.sha256((SEED + ':' + pid).encode()).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / 'locked_plan.json'
    if path.exists():
        print('Existing plan retained:', path)
        return
    records = []
    for p in sorted((ROOT / 'data/derived/training_cohort/series').glob('*.json')):
        rows = json.loads(p.read_text())
        extraction = json.loads((ROOT / 'data/derived/training_cohort/patients' / p.name).read_text())
        visits = {}
        for visit in ['T0', 'T1']:
            candidates = [r for r in rows if r['StudyDesc'].endswith('_' + visit)
                          and 'original DCE' in r.get('SeriesDescription', '')]
            assert len(candidates) == 1
            ref = extraction['visits'][visit]['inputs']['ser']
            d = pydicom.dcmread(ROOT / ref, stop_before_pixels=True)
            assert str(d.PatientID) == p.stem and str(d.ClinicalTrialTimePointID) == visit
            visits[visit] = {'series': candidates[0], 'reference': ref,
                'ser_indices': [int(x) for x in d[(0x0117, 0x1035)].value],
                'site': str(getattr(d, 'ClinicalTrialSiteID', 'unknown'))}
        records.append({'patient_id': p.stem, 'manufacturer': visits['T1']['series']['Manufacturer'],
                        'visits': visits})
    groups = defaultdict(list)
    for r in records:
        groups[r['manufacturer']].append(r)
    for vendor, group in groups.items():
        ordered = sorted(group, key=lambda r: (r['patient_id'] != PID, rank(r['patient_id'])))
        nvalidation = max(1, round(len(group) * .25))
        for i, r in enumerate(ordered):
            r['split'] = 'development' if i < len(group) - nvalidation else 'internal_validation'
            r['pilot'] = i < 2 or i == len(group) - nvalidation
    assert next(r for r in records if r['patient_id'] == PID)['split'] == 'development'
    assert len(records) == 117 and len({r['patient_id'] for r in records}) == 117
    plan = {'version': SEED, 'records': records,
        'counts': dict(Counter(r['split'] for r in records)),
        'pilot_counts': dict(Counter(r['split'] for r in records if r['pilot'])),
        'scope': 'internal exploratory development; original 117 previously used for pCR, not a new confirmatory cohort',
        'unit': 'patient; all visits and slices together',
        'selection': 'SHA256 within manufacturer; engineering patient forced first in development; no outcome reads',
        'preprocessing': 'published crop engineering pilot only; pre-positive p95 scaling; no target normalization or registration',
        'native_input_gate': 'native full field or crop determined from allowed inputs required before main synthesis claims',
        'model': {'name': 'ridge voxelwise residual', 'alpha': 1.0, 'sample_voxels_per_visit': 12000,
                  'features': ['pre/scale', '(early-pre)/scale', 'early_seconds/600', 'late_seconds/600'],
                  'target': '(late-early)/scale', 'seed': 20260917},
        'evaluation': 'same input-derived masks as pilot_v1; patient means over visits; compare copy early; no tuning on validation',
        'stop_rules': ['Do not claim paper Go from a nine-patient pilot.',
            'Stop main-model claims if native input policy cannot be implemented without hidden target information.',
            'Do not expand a failed simple baseline as evidence of success; review patient-level harm and technical failures.',
            'Stop simple synthesis-paper framing if prior art already answers the proposed question.'],
        'download_cap_bytes': 1024**3}
    text = json.dumps(plan, indent=2)
    path.write_text(text, encoding='utf-8')
    (OUT / 'locked_plan.sha256').write_text(hashlib.sha256(text.encode()).hexdigest(), encoding='ascii')
    print(json.dumps({k:v for k,v in plan.items() if k != 'records'}, indent=2))
    print('Pilot:', [(r['patient_id'], r['manufacturer'], r['split']) for r in records if r['pilot']])


if __name__ == '__main__':
    main()
