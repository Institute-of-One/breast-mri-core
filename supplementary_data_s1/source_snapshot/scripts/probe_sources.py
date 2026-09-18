"""Download public documentation and metadata only for two predefined training cases."""
from pathlib import Path
import json
import requests
import csv
import hashlib
from datetime import datetime, timezone
from pypdf import PdfReader
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'data/derived/probe'
BASE = 'https://services.cancerimagingarchive.net/nbia-api/services/v1/'

def get(url, **kw):
    r=requests.get(url,timeout=60,**kw); r.raise_for_status(); return r

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    errors=[]
    sources = {
        'data_descriptions.pdf':'https://www.cancerimagingarchive.net/wp-content/uploads/ACRIN-6698-ISPY2-DWI-and-DCE-MRI-Data-Descriptions_20210520.pdf',
        'bmmr2.html':'https://pmc.ncbi.nlm.nih.gov/articles/PMC10825718/',
        'acrin_primary.html':'https://pmc.ncbi.nlm.nih.gov/articles/PMC6283325/',
        'collection.html':'https://www.cancerimagingarchive.net/collection/acrin-6698/'
    }
    for name,url in sources.items():
        try:
            path=OUT/name
            if not path.exists(): path.write_bytes(get(url).content)
            if name.endswith('.pdf'):
                (OUT/(name+'.txt')).write_text('\n'.join(f'PAGE {i+1}\n'+(p.extract_text() or '') for i,p in enumerate(PdfReader(path).pages)),encoding='utf-8')
            print(name,path.stat().st_size,flush=True)
        except Exception as e:
            errors.append({'resource':name,'error':str(e)})
            print(name,type(e).__name__,str(e)[:180],flush=True)
    with next((ROOT/'data/raw/acrin6698').glob('*Training*.csv')).open() as f: train={r['Patient ID DICOM'] for r in csv.DictReader(f)}
    for pid in ['ACRIN-6698-107700','ACRIN-6698-373346']:
        assert pid in train
        for endpoint in ['getPatientStudy','getSeries']:
            try:
                path=OUT/(pid+'_'+endpoint+'.json')
                if not path.exists(): path.write_text(json.dumps(get(BASE+endpoint,params={'Collection':'ACRIN-6698','PatientID':pid}).json(),indent=2),encoding='utf-8')
                data=json.loads(path.read_text()); print(pid,endpoint,len(data),data[:1],flush=True)
            except Exception as e:
                errors.append({'resource':pid+'_'+endpoint,'error':str(e)})
                print(pid,endpoint,type(e).__name__,str(e)[:180],flush=True)
    manifest={'checked_utc':datetime.now(timezone.utc).isoformat(),'sources':sources,'api_base':BASE,'errors':errors,'file_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file() and (p.name in sources or p.name.endswith(('_getPatientStudy.json','_getSeries.json')))}}
    (OUT/'retrieval_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    if errors: raise SystemExit(1)

if __name__=='__main__': main()
