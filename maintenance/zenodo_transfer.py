"""Transfer verified public release assets into the existing Zenodo draft.

One-off archival maintenance only; no analysis, metadata edits or publication.
ZENODO_UPLOAD_TOKEN is supplied by the hosting platform's encrypted secret store.
Only the authorized draft is touched. Existing matching files are retained.
"""
import hashlib
import http.client
import json
import os
from pathlib import Path
import tempfile
import time
import urllib.request
import urllib.error
from concurrent.futures import ThreadPoolExecutor,as_completed
from urllib.parse import urlparse

BASE='https://zenodo.org/api/records/24073205/draft/files'
TOKEN=os.environ.pop('ZENODO_UPLOAD_TOKEN')
EXPECTED=json.loads(Path(__file__).with_name('zenodo_expected_files.json').read_text())
assert len(EXPECTED)==29

def api(method,url,payload=None):
    assert url.startswith(BASE)
    data=None if payload is None else json.dumps(payload).encode()
    request=urllib.request.Request(url,data=data,method=method,headers={
        'Authorization':'Bearer '+TOKEN,
        'User-Agent':'ChinaEHI-ArchivalTransfer/1.0',
        'Content-Type':'application/json'})
    with urllib.request.urlopen(request,timeout=180) as response:
        body=response.read()
        return json.loads(body) if body else None

def files():return {f['key']:f for f in api('GET',BASE)['entries']}
def match(f,e):
    return (f['status']=='completed' and f['size']==e['bytes']
            and f.get('checksum','').removeprefix('md5:').zfill(32)==e['md5'].zfill(32))

remote=files()
assert not set(remote)-set(EXPECTED)
buckets={f['bucket_id'] for f in remote.values()}
assert len(buckets)==1
BUCKET='https://zenodo.org/api/files/'+buckets.pop()
for name,f in remote.items():
    if f['status']=='completed':
        assert match(f,EXPECTED[name]),name
    else:
        # Stop the local uploader before dispatch. Only empty partial slots
        # created by that uploader are eligible for removal and reconstruction.
        assert f['status']=='pending' and f.get('size') in (None,0),name
        api('DELETE',f['links']['self'])
        print('Removed incomplete slot:',name,flush=True)

def transfer_one(item):
    name,e=item
    remote=files()
    if name in remote:
        assert match(remote[name],e),name
        print('Verified existing:',name,flush=True)
        return
    assert e['source_url'].startswith('https://github.com/lizexuan581-stack/china-ehi-continuous-domain-2000-2025/releases/download/v1.0.0/')
    begin=time.monotonic()
    with tempfile.TemporaryDirectory() as tmp:
        path=Path(tmp)/name
        for attempt in range(5):
            md5=hashlib.md5();sha256=hashlib.sha256();size=0
            try:
                request=urllib.request.Request(e['source_url'],headers={'User-Agent':'ChinaEHI-ArchivalTransfer/1.0'})
                with urllib.request.urlopen(request,timeout=180) as response,path.open('wb') as out:
                    while chunk:=response.read(1024*1024):
                        out.write(chunk);md5.update(chunk);sha256.update(chunk);size+=len(chunk)
                break
            except (urllib.error.URLError,TimeoutError) as exc:
                if isinstance(exc,urllib.error.HTTPError) and exc.code not in [429,500,502,503,504]:raise
                if attempt==4:raise
                print('Retrying public download:',name,'attempt:',attempt+2,flush=True)
                time.sleep(10*(attempt+1))
        assert size==e['bytes'] and md5.hexdigest()==e['md5'] and sha256.hexdigest()==e['sha256'],name
        print('Downloaded and verified:',name,flush=True)
        # Official Zenodo large-file bucket API initializes and commits the
        # draft entry itself. A small checksum-file probe verified this route.
        url=BUCKET+'/'+name
        parsed=urlparse(url)
        connection=http.client.HTTPSConnection(parsed.hostname,timeout=600)
        try:
            connection.putrequest('PUT',parsed.path)
            connection.putheader('Authorization','Bearer '+TOKEN)
            connection.putheader('Content-Type','application/octet-stream')
            connection.putheader('Content-Length',str(size))
            connection.putheader('User-Agent','ChinaEHI-ArchivalTransfer/1.0')
            connection.endheaders()
            transferred=0;reported=0
            with path.open('rb') as source:
                while chunk:=source.read(1024*1024):
                    connection.send(chunk);transferred+=len(chunk)
                    if transferred-reported>=64*1024*1024:
                        reported=transferred
                        print('Upload progress:',name,transferred,'of',size,flush=True)
            response=connection.getresponse();body=response.read()
            assert 200<=response.status<300,(name,response.status)
            uploaded=json.loads(body)
            assert uploaded['size']==e['bytes'] and uploaded['checksum'].removeprefix('md5:').zfill(32)==e['md5'].zfill(32),name
        finally:connection.close()
        assert match(files()[name],e),name
        print('Uploaded and verified:',name,'seconds:',round(time.monotonic()-begin,1),flush=True)

with ThreadPoolExecutor(max_workers=3) as pool:
    futures=[pool.submit(transfer_one,item) for item in EXPECTED.items()]
    for future in as_completed(futures):future.result()
remote=files()
assert set(remote)==set(EXPECTED)
assert all(match(remote[name],e) for name,e in EXPECTED.items())
print('ALL_29_FILES_VERIFIED. Publication remains a separate browser action.',flush=True)
