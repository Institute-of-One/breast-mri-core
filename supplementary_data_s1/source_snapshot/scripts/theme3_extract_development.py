"""Bounded nine-patient engineering extraction; published crop limitations apply."""
import argparse
import hashlib
import io
import json
import threading
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pydicom
from theme3_start import ROOT, BASE, get
from theme3_lock_development import OUT

RAW = ROOT / 'data/raw/theme3_development_v1'
LOCK = threading.Lock()
used = 0


def fetch(row, pid, visit, online, cap):
    global used
    old = ROOT / f'data/raw/theme3_pilot_v1/{pid}_{visit}_published_crop.zip'
    path = old if old.exists() else RAW / f'{pid}_{visit}.zip'
    if path.exists():
        return path
    if not online:
        raise RuntimeError('Missing source; use --download')
    partial = path.with_suffix('.part')
    reserved = 0
    success = False
    try:
        with get(BASE + 'getImage', params={'SeriesInstanceUID': row['SeriesInstanceUID']}, stream=True) as r:
            with partial.open('wb') as f:
                for chunk in r.iter_content(1024 * 1024):
                    with LOCK:
                        if used + len(chunk) > cap:
                            raise RuntimeError('1 GiB compressed new-download cap exceeded')
                        used += len(chunk)
                        reserved += len(chunk)
                    f.write(chunk)
        with zipfile.ZipFile(partial) as z:
            assert z.testzip() is None
        partial.replace(path)
        success = True
    finally:
        if partial.exists():
            partial.unlink()
        if not success:
            with LOCK:
                used -= reserved
    return path


def prepare(path, record, visit):
    pid = record['patient_id']
    row = record['visits'][visit]['series']
    refpath = ROOT / record['visits'][visit]['reference']
    ref = pydicom.dcmread(refpath, stop_before_pixels=True)
    indices = record['visits'][visit]['ser_indices']
    with zipfile.ZipFile(path) as z:
        datasets = [pydicom.dcmread(io.BytesIO(z.read(n))) for n in z.namelist() if n.lower().endswith('.dcm')]
    assert len(datasets) == int(row['ImageCount'])
    assert len({str(d.SOPInstanceUID) for d in datasets}) == len(datasets)
    first = datasets[0]
    assert 'UCSF BIRP' in str(first[(0x0117, 0x0010)].value)
    nphase = int(first[(0x0117, 0x1030)].value)
    delays = np.array(first[(0x0117, 0x1034)].value, float)
    assert len(delays) == nphase and np.all(np.isfinite(delays))
    assert indices[0] == 0 and 0 < indices[1] < indices[2] < nphase
    assert 0 < delays[indices[1]] < delays[indices[2]]
    assert np.allclose(delays, ref[(0x0117, 0x1034)].value)
    assert str(ref.StudyInstanceUID) == row['StudyInstanceUID']
    assert str(ref.PatientID) == pid
    groups = {}
    for d in datasets:
        assert str(d.PatientID) == pid and str(d.ClinicalTrialTimePointID) == visit
        assert str(d.SeriesInstanceUID) == row['SeriesInstanceUID'] and str(d.StudyInstanceUID) == row['StudyInstanceUID']
        assert d.FrameOfReferenceUID == first.FrameOfReferenceUID
        assert np.allclose(d.ImageOrientationPatient, first.ImageOrientationPatient)
        assert np.allclose(d.PixelSpacing, first.PixelSpacing)
        assert (d.Rows, d.Columns) == (first.Rows, first.Columns)
        groups.setdefault(int(d.TemporalPositionIdentifier), []).append(d)
    assert sorted(groups) == list(range(nphase)), 'Unexpected phase identifiers'
    normal = np.cross(np.array(first.ImageOrientationPatient[:3], float), np.array(first.ImageOrientationPatient[3:], float))
    for group in groups.values():
        group.sort(key=lambda d: float(np.dot(normal, np.array(d.ImagePositionPatient, float))))
    positions = np.array([d.ImagePositionPatient for d in groups[0]], float)
    assert len({tuple(p) for p in positions}) == len(positions)
    assert all(len(g) == len(positions) and np.allclose(
        [d.ImagePositionPatient for d in g], positions, atol=.001, rtol=0) for g in groups.values())
    arrays = {role: np.stack([d.pixel_array.astype(np.float32) * float(getattr(d, 'RescaleSlope', 1))
                             + float(getattr(d, 'RescaleIntercept', 0)) for d in groups[index]])
              for role, index in zip(['pre', 'early', 'late'], indices)}
    assert all(np.isfinite(a).all() for a in arrays.values())
    target = OUT / f'{pid}_{visit}.npz'
    np.savez_compressed(target, **arrays, positions=positions,
                        orientation=np.array(first.ImageOrientationPatient, float),
                        pixel_spacing=np.array(first.PixelSpacing, float))
    return {'status': 'pass', 'pair_file': str(target.relative_to(ROOT)),
            'source': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'phases': nphase, 'shape': list(arrays['pre'].shape), 'indices': indices,
            'delays_seconds': delays[indices].tolist(), 'site': str(getattr(first, 'ClinicalTrialSiteID', 'unknown')),
            'native_input_policy_passed': False, 'motion_registered': False}


def main():
    global used
    parser = argparse.ArgumentParser()
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args()
    # Plan hash uses canonical LF text, independent of Windows newline encoding.
    payload = (OUT / 'locked_plan.json').read_text(encoding='utf-8').encode('utf-8')
    assert hashlib.sha256(payload).hexdigest() == (OUT / 'locked_plan.sha256').read_text()
    plan = json.loads(payload)
    RAW.mkdir(parents=True, exist_ok=True)
    used = sum(p.stat().st_size for p in RAW.glob('*.zip'))
    def one(record):
        pid = record['patient_id']
        output = OUT / f'{pid}_extraction.json'
        if output.exists():
            previous = json.loads(output.read_text())
            if all(x['status'] == 'pass' for x in previous['visits'].values()):
                return previous
        result = {'patient_id': pid, 'split': record['split'], 'manufacturer': record['manufacturer'], 'visits': {}}
        for visit in ['T0', 'T1']:
            try:
                source = fetch(record['visits'][visit]['series'], pid, visit, args.download, plan['download_cap_bytes'])
                result['visits'][visit] = prepare(source, record, visit)
            except Exception as e:
                result['visits'][visit] = {'status': 'failed', 'error': repr(e)}
            print(pid, visit, result['visits'][visit]['status'], flush=True)
        output.write_text(json.dumps(result, indent=2), encoding='utf-8')
        return result
    results = []
    with ThreadPoolExecutor(max_workers=2) as pool:
        tasks = [pool.submit(one, r) for r in plan['records'] if r['pilot']]
        for future in as_completed(tasks):
            results.append(future.result())
    summary = {'patients': sorted(results, key=lambda r: r['patient_id']),
               'complete_patients': sum(all(v['status'] == 'pass' for v in r['visits'].values()) for r in results),
               'new_folder_compressed_bytes': used}
    (OUT / 'extraction_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print('Complete patients:', summary['complete_patients'], '/', len(results), flush=True)


if __name__ == '__main__':
    main()
