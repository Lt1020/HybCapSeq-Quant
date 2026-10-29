#!/usr/bin/env python3
"""Re-download only the frozen diagnostic windows from NCBI and compare bases."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import time
import urllib.request

def sequence(data):
    lines=data.decode('ascii').splitlines()
    assert lines and lines[0].startswith('>'),'NCBI did not return FASTA'
    return lines[0][1:].split()[0],''.join(s.strip().upper() for s in lines[1:])

def refresh(inputs,output):
    if output.exists():raise ValueError('Choose a new output directory')
    output.mkdir(parents=True);(output/'source_windows').mkdir()
    specs=json.loads((inputs/'source_catalog.json').read_text())
    for s in specs:
        source=inputs/'source_windows'/s['source_file']
        assert hashlib.sha256(source.read_bytes()).hexdigest()==s['source_fasta_sha256']
        expected_header,expected=sequence(source.read_bytes())
        with urllib.request.urlopen(s['source_fetch_url'],timeout=45) as response:data=response.read()
        header,bases=sequence(data)
        assert header.startswith(s['source_accession_version']),'Unexpected accession version'
        assert bases==expected,'NCBI reference bases differ from the saved source; review before updating'
        # Keep original frozen file bytes and provenance after verifying the bases.
        shutil.copy2(source,output/'source_windows'/s['source_file'])
        print(s['target'],s['assay'],s['source_accession_version'],'MATCH',flush=True)
        time.sleep(0.4)
    for p in inputs.iterdir():
        if p.is_file():shutil.copy2(p,output/p.name)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();refresh(args.inputs.resolve(),args.output.resolve())
