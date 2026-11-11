"""Outcome-free DCE inventory and bounded, training-only temporal parsing pilot.

Published tumor-side crops are engineering inputs, not deployment-ready inputs.
No clinical labels, original model configuration, or held-out data are read.
"""
import argparse
import hashlib
import io
import json
import zipfile
from collections import Counter
from pathlib import Path

import numpy as np
import pydicom
from probe_sources import ROOT, BASE, get

OUT = ROOT / 'data/derived/theme3_pilot_v1'
RAW = ROOT / 'data/raw/theme3_pilot_v1'
PID = 'ACRIN-6698-107700'
CAP = 256 * 1024 * 1024


def inventory():
    records = []
    for path in sorted((ROOT / 'data/derived/training_cohort/series').glob('*.json')):
        rows = json.loads(path.read_text())
        assert rows and all(r['PatientID'] == path.stem for r in rows)
        for visit in ['T0', 'T1']:
            matching = [r for r in rows if r.get('StudyDesc', '').endswith('_' + visit)
                        and 'original DCE' in r.get('SeriesDescription', '')]
            records.append({'patient_id': path.stem, 'visit': visit,
                            'status': 'unique_crop_candidate' if len(matching) == 1 else 'review_required',
                            'candidates': matching})
    result = {'scope': 'cached original training metadata; no outcome reads',
              'patients': len({r['patient_id'] for r in records}),
              'status_counts': dict(Counter(r['status'] for r in records)),
              'both_visits_unique': sum(all(any(r['patient_id'] == pid and r['visit'] == v
                    and r['status'] == 'unique_crop_candidate' for r in records) for v in ['T0', 'T1'])
                    for pid in {r['patient_id'] for r in records}),
              'warning': 'Series presence does not establish usable temporal pairs. Published crops may use full DCE knowledge.',
              'records': records}
    (OUT / 'inventory.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    return result


def fetch(row, visit):
    path = RAW / f'{PID}_{visit}_published_crop.zip'
    if not path.exists():
        used = sum(p.stat().st_size for p in RAW.glob('*.zip'))
        partial = path.with_suffix('.part')
        try:
            with get(BASE + 'getImage', params={'SeriesInstanceUID': row['SeriesInstanceUID']}, stream=True) as response:
                size = 0
                with partial.open('wb') as f:
                    for chunk in response.iter_content(1024 * 1024):
                        size += len(chunk)
                        if used + size > CAP:
                            raise RuntimeError('256 MiB compressed pilot cap exceeded')
                        f.write(chunk)
            with zipfile.ZipFile(partial) as z:
                assert z.testzip() is None
            partial.replace(path)
        finally:
            if partial.exists():
                partial.unlink()
    return path


def inspect(path, row, visit):
    with zipfile.ZipFile(path) as z:
        datasets = [pydicom.dcmread(io.BytesIO(z.read(n))) for n in z.namelist() if n.lower().endswith('.dcm')]
    assert len(datasets) == int(row['ImageCount'])
    assert all(str(d.PatientID) == PID and str(d.ClinicalTrialTimePointID) == visit
               and str(d.SeriesInstanceUID) == row['SeriesInstanceUID']
               and str(d.StudyInstanceUID) == row['StudyInstanceUID'] for d in datasets)
    assert len({str(d.SOPInstanceUID) for d in datasets}) == len(datasets)
    keys = ['TemporalPositionIdentifier', 'NumberOfTemporalPositions', 'AcquisitionNumber',
            'AcquisitionTime', 'ContentTime', 'TriggerTime', 'ContrastBolusStartTime',
            'InstanceNumber', 'ClinicalTrialSiteID']
    first = datasets[0]
    report = {'patient_id': PID, 'visit': visit, 'instances': len(datasets),
              'file': str(path.relative_to(ROOT)), 'bytes': path.stat().st_size,
              'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
              'series_uid': row['SeriesInstanceUID'], 'license': row.get('LicenseURI'),
              'shape': [int(first.Rows), int(first.Columns)],
              'field_unique_values': {k: sorted({str(getattr(d, k, '')) for d in datasets}) for k in keys},
              'positions': len({tuple(d.ImagePositionPatient) for d in datasets}),
              'private_tags_first_instance': [str(e) for e in first if e.tag.group == 0x0117],
              'status': 'headers_inspected_temporal_assignment_pending',
              'scope': 'published crop engineering audit only; no synthesis accuracy or clinical result'}
    (OUT / f'{visit}_headers.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    inv = inventory()
    print(json.dumps({k: v for k, v in inv.items() if k != 'records'}, indent=2), flush=True)
    if args.download:
        for visit in ['T0', 'T1']:
            record = next(r for r in inv['records'] if r['patient_id'] == PID and r['visit'] == visit)
            assert record['status'] == 'unique_crop_candidate'
            row = record['candidates'][0]
            report = inspect(fetch(row, visit), row, visit)
            print(visit, report['instances'], 'instances;', report['positions'], 'positions', flush=True)


if __name__ == '__main__':
    main()
