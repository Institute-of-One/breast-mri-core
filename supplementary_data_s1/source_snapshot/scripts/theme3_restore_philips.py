"""Preserve native Philips rescaled intensities on the existing engineering grid.

No inferred quantization correction: sample actual native pixels at matching
physical coordinates. Original crop pairs are retained unchanged.
"""
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

RAW = ROOT / 'data/raw/theme3_native_probe'
SERIES = {'ACRIN-6698-342959': [701,1001], 'ACRIN-6698-618860': [1301,801],
          'ACRIN-6698-861745': [1101,801]}
LOCK = threading.Lock()
used = 0


def download(row,pid,visit):
    global used
    path=RAW/f'{pid}_{visit}_native.zip'
    if path.exists():return path
    partial=path.with_suffix('.part'); reserved=0; success=False
    try:
        with get(BASE+'getImage',params={'SeriesInstanceUID':row['SeriesInstanceUID']},stream=True) as response:
            with partial.open('wb') as f:
                for chunk in response.iter_content(1024**2):
                    with LOCK:
                        if used+len(chunk)>2*1024**3:raise RuntimeError('Native folder 2 GiB cap')
                        used+=len(chunk);reserved+=len(chunk)
                    f.write(chunk)
        with zipfile.ZipFile(partial) as z:assert z.testzip() is None
        partial.replace(path);success=True
    finally:
        if partial.exists():partial.unlink()
        if not success:
            with LOCK:used-=reserved
    return path


def prepare(pid,visit,row,path):
    record=json.loads((OUT/f'{pid}_extraction.json').read_text())['visits'][visit]
    crop=np.load(ROOT/record['pair_file'])
    result={};native_counts={}
    with zipfile.ZipFile(path) as z:
        groups={}
        for name in z.namelist():
            if not name.lower().endswith('.dcm'):continue
            d=pydicom.dcmread(io.BytesIO(z.read(name)),stop_before_pixels=True)
            assert str(d.PatientID)==pid and str(d.SeriesInstanceUID)==row['SeriesInstanceUID']
            assert str(d.StudyInstanceUID)==row['StudyInstanceUID'] and str(d.ClinicalTrialTimePointID)==visit
            groups.setdefault(int(d.TemporalPositionIdentifier)-1,[]).append((name,d))
        assert sum(map(len,groups.values()))==int(row['ImageCount'])
        assert sorted(groups)==list(range(record['phases']))
        first=groups[0][0][1]
        orientation=np.array(first.ImageOrientationPatient,float)
        normal=np.cross(orientation[:3],orientation[3:])
        for g in groups.values():g.sort(key=lambda pair:float(np.dot(normal,np.array(pair[1].ImagePositionPatient,float))))
        pos=np.array([d.ImagePositionPatient for _,d in groups[0]],float)
        spacing=float(np.diff(pos@normal)[0])
        assert spacing>0 and np.allclose(np.diff(pos@normal),spacing,atol=.001,rtol=0)
        assert all(np.allclose([d.ImagePositionPatient for _,d in g],pos,atol=.001,rtol=0) for g in groups.values())
        assert len({str(d.SOPInstanceUID) for g in groups.values() for _,d in g})==int(row['ImageCount'])
        for g in groups.values():
            for _,d in g:
                assert d.FrameOfReferenceUID==first.FrameOfReferenceUID
                assert np.allclose(d.ImageOrientationPatient,orientation)
                assert np.allclose(d.PixelSpacing,first.PixelSpacing)
        height,width=crop['pre'].shape[1:]
        yy,xx=np.mgrid[:height,:width]
        plane_offset=(yy[...,None]*crop['orientation'][3:]*crop['pixel_spacing'][0]
                      +xx[...,None]*crop['orientation'][:3]*crop['pixel_spacing'][1])
        indices=[]
        for origin in crop['positions']:
            world=origin+plane_offset-pos[0]
            coordinates=np.stack([world@normal/spacing,world@orientation[3:]/float(first.PixelSpacing[0]),
                                  world@orientation[:3]/float(first.PixelSpacing[1])])
            assert np.max(np.abs(coordinates-np.rint(coordinates)))<.002, 'Native/crop interpolation required: do not guess'
            indices.append(np.rint(coordinates).astype(int))
        for role,phase in zip(['pre','early','late'],record['indices']):
            volume=np.stack([pydicom.dcmread(io.BytesIO(z.read(name))).pixel_array.astype(float)*float(getattr(d,'RescaleSlope',1))
                             +float(getattr(d,'RescaleIntercept',0)) for name,d in groups[phase]])
            planes=[]
            for coordinate in indices:
                assert all(coordinate[ax].min()>=0 and coordinate[ax].max()<volume.shape[ax] for ax in range(3))
                planes.append(volume[tuple(coordinate)])
            result[role]=np.stack(planes)  # Float64 intentionally retains rescale precision.
            native_counts[role]={'native_shape':list(volume.shape),
                'rescale_slopes':sorted({float(getattr(d,'RescaleSlope',1)) for _,d in groups[phase]})}
    target=OUT/f'{pid}_{visit}_native_precision.npz'
    np.savez_compressed(target,**result,positions=crop['positions'],orientation=crop['orientation'],pixel_spacing=crop['pixel_spacing'])
    return {'patient_id':pid,'visit':visit,'status':'pass','pair_file':str(target.relative_to(ROOT)),
        'source':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'native_counts':native_counts,'scope':'Native intensities on reference crop for engineering only; not input-only crop'}


def main():
    global used
    used=sum(p.stat().st_size for p in RAW.glob('*.zip'))
    def one(pid):
        rows=json.loads((ROOT/f'data/derived/training_cohort/series/{pid}.json').read_text())
        output=[]
        for visit,number in zip(['T0','T1'],SERIES[pid]):
            report_path=OUT/f'{pid}_{visit}_native_precision.json'
            if report_path.exists() and json.loads(report_path.read_text()).get('status')=='pass':
                output.append(json.loads(report_path.read_text()));continue
            try:
                candidates=[r for r in rows if r['StudyDesc'].endswith('_'+visit) and int(r['SeriesNumber'])==number]
                assert len(candidates)==1
                path=download(candidates[0],pid,visit)
                report=prepare(pid,visit,candidates[0],path)
            except Exception as e:
                report={'patient_id':pid,'visit':visit,'status':'failed','error':repr(e)}
            report_path.write_text(json.dumps(report,indent=2),encoding='utf-8')
            output.append(report)
            print(pid,visit,report['status'],report.get('error',''),flush=True)
        return output
    reports=[]
    with ThreadPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(one,pid) for pid in SERIES]):reports.extend(future.result())
    (OUT/'native_precision_summary.json').write_text(json.dumps(reports,indent=2),encoding='utf-8')


if __name__=='__main__':main()
