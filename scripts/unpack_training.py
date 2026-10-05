#!/usr/bin/env python3
"""Validate compressed and extracted training inputs using the committed manifest."""
import hashlib, json, pathlib, zipfile
root=pathlib.Path(__file__).resolve().parents[1]
records=json.loads((root/'data/manifest.json').read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()
for r in records:
    archive=root/r['archive']
    if sha(archive)!=r['archive_sha256']: raise SystemExit(f'Archive checksum mismatch: {archive.name}')
    member=pathlib.PurePosixPath(r['member'])
    if member.is_absolute() or '..' in member.parts: raise SystemExit('Unsafe archive member')
    output=root/'data/extracted'/member
    if not output.exists() or sha(output)!=r['sha256']:
        with zipfile.ZipFile(archive) as z:
            if z.namelist()!=[r['member']]: raise SystemExit('Unexpected archive contents')
            output.parent.mkdir(parents=True,exist_ok=True)
            with z.open(r['member']) as incoming,output.open('wb') as out:
                import shutil
                shutil.copyfileobj(incoming,out)
    if sha(output)!=r['sha256']: raise SystemExit(f'Extracted checksum mismatch: {member}')
    print(f'OK {member}')
