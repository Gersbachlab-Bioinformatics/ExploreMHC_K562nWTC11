#!/usr/bin/env python3
"""Fetch only WTC11 lowQ BED members via verified-TLS curl byte-range requests."""
import io
import json
import hashlib
from pathlib import Path
import subprocess
import tarfile

URL = 'https://s3.kopah.uw.edu/sjn/web/public/project/g.crawford/asm.tar'
ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'data/wtc11_lowq'
TARGETS = {f'asm/assembly.asm.bp.hap{h}.p_ctg.lowQ.bed' for h in [1, 2]}

class Remote(io.RawIOBase):
    def __init__(self):
        self.pos = 0
        self.transferred = 0
    def tell(self): return self.pos
    def seek(self, offset, whence=0):
        if whence == 0: self.pos = offset
        elif whence == 1: self.pos += offset
        else: raise ValueError('SEEK_END not supported')
        return self.pos
    def readable(self): return True
    def seekable(self): return True
    def read(self, n=-1):
        if n < 0 or n > 16*1024*1024: raise ValueError(f'Unexpected large read: {n}')
        if n == 0: return b''
        data = subprocess.check_output(['curl','--fail','--silent','--show-error','--max-time','40',
            '--max-filesize',str(n),'--header','Accept-Encoding: identity',
            '--range',f'{self.pos}-{self.pos+n-1}',URL])
        if len(data) != n: raise ValueError(f'Short range response: {len(data)} vs {n}')
        self.pos += n
        self.transferred += n
        return data

if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    remote = Remote()
    found, listing = [], []
    with tarfile.open(fileobj=remote, mode='r:') as tar:
        for member in tar:
            listing.append(dict(name=member.name,size=member.size,offset=member.offset_data))
            print(member.name,member.size,flush=True)
            if member.name in TARGETS:
                if not member.isfile(): raise ValueError('Expected a regular BED file')
                data=tar.extractfile(member).read()
                dest=OUT / Path(member.name).name
                dest.write_bytes(data)
                found.append(dict(member=member.name,offset=member.offset_data,size=member.size,
                                  local=str(dest),sha256=hashlib.sha256(data).hexdigest()))
    (OUT/'archive_members.json').write_text(json.dumps(listing,indent=2)+'\n')
    if {r['member'] for r in found} != TARGETS: raise ValueError('Missing expected BED members')
    (OUT/'download_manifest.json').write_text(json.dumps(dict(url=URL,members=found,bytes_transferred=remote.transferred),indent=2)+'\n')
    print(f'Fetched {len(found)} BED files; total bytes transferred: {remote.transferred:,}')
