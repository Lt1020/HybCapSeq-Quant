#!/usr/bin/env python3
"""Independently check the emitted diagnostic-window files (standard library)."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import urllib.parse
import xml.etree.ElementTree as ET

COLUMNS='Accession description Name Segment kingdom phylum tclass order family genus species subspecies Length TaxID Assembly Asm_length'.split()
BASES={'A':'A','C':'C','G':'G','T':'T','R':'AG','Y':'CT','S':'CG','W':'AT','K':'GT','M':'AC','B':'CGT','D':'AGT','H':'ACT','V':'ACG','N':'ACGT','I':'ACGT'}
COMPLEMENT=dict(zip('ACGTRYSWKMBDHVNI','TGCAYRSWMKVHDBNI'))

def read_tsv(path):
    with path.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def fasta(path):
    rows={};key=None
    with path.open(encoding='ascii') as f:
        for line in f:
            if line.startswith('>'):
                key=line[1:].strip().split()[0]
                assert key not in rows,(path,'duplicate ID',key)
                rows[key]=''
            elif line.strip():
                assert key is not None
                rows[key]+=line.strip().upper()
    assert all(rows.values())
    for seq in rows.values():assert set(seq)<=set(BASES)-{'I'}
    return rows

def reverse_complement(seq):return ''.join(COMPLEMENT[b] for b in reversed(seq))

def validate_package(root):
    catalog=json.loads((root/'sources/source_catalog.json').read_text())
    specs={(s['target'],s['assay']):s for s in catalog}
    assert len(specs)==len(catalog)
    allrefs=fasta(root/'all_candidates/virus_pathogen_database.fna')
    taxonomy={}
    for node in ET.parse(root/'sources/taxonomy_snapshot.xml').getroot().findall('Taxon'):
        lineage={x.findtext('Rank'):x.findtext('ScientificName') for x in node.findall('./LineageEx/Taxon')}
        lineage[node.findtext('Rank')]=node.findtext('ScientificName')
        taxonomy[int(node.findtext('TaxId'))]=lineage
    all_metadata={m['Accession']:m for m in read_tsv(root/'all_candidates/virus_pathogen_database.all_metadata.tsv')}
    groupids={}
    for sub in ['', 'all_candidates', 'exceptions/boundary_limited', 'exceptions/oligo_mismatch']:
        base=root/sub;seqs=fasta(base/'virus_pathogen_database.fna')
        with (base/'virus_pathogen_database.all_metadata.tsv').open(encoding='utf-8',newline='') as f:
            reader=csv.DictReader(f,delimiter='\t');assert reader.fieldnames==COLUMNS;meta=list(reader)
        assert len(meta)==len(seqs)
        assert {m['Accession'] for m in meta}==set(seqs)
        for m in meta:
            ident=m['Accession'];n=len(seqs[ident])
            assert re.fullmatch(r'[A-Za-z0-9_.-]+',ident)
            assert int(m['Length'])==n==int(m['Asm_length'])
            assert m['Name']==m['Assembly']
            assert allrefs[ident]==seqs[ident]
        groupids[sub]=set(seqs)
    primary=groupids[''];boundary=groupids['exceptions/boundary_limited'];mismatch=groupids['exceptions/oligo_mismatch']
    assert not (primary&boundary or primary&mismatch or boundary&mismatch)
    assert primary|boundary|mismatch==set(allrefs)
    regions=read_tsv(root/'reference_regions.tsv')
    assert len(regions)==len(catalog)
    byassay={}
    expectedseq={}
    for r in regions:
        s=specs[(r['target'],r['assay'])];rid=r['reference_id']
        metadata=all_metadata[rid]
        assert int(metadata['TaxID'])==s['taxid']
        for rank in ['kingdom','phylum','order','family','genus','species','subspecies']:
            assert metadata[rank]==taxonomy[s['taxid']].get(rank,'NA')
        assert metadata['tclass']==taxonomy[s['taxid']].get('class','NA')
        src=root/'sources/source_windows'/s['source_file']
        assert hashlib.sha256(src.read_bytes()).hexdigest()==s['source_fasta_sha256']
        parsed=fasta(src);assert len(parsed)==1
        header,seq=next(iter(parsed.items()));assert header.startswith(s['source_accession_version'])
        assert len(seq)==s['source_window_end']-s['source_window_start']+1
        start,end=int(r['reference_start_1based']),int(r['reference_end_1based'])
        left,right=int(r['amplicon_start_1based']),int(r['amplicon_end_1based'])
        assert start==max(1,left-100) and end==min(s['source_length_nt'],right+100)
        assert int(r['reference_length_nt'])==end-start+1
        assert int(r['amplicon_length_nt'])==right-left+1
        assert int(r['left_flank_nt'])==left-start
        assert int(r['right_flank_nt'])==end-right
        assert int(r['amplicon_start_in_output_1based'])==left-start+1
        assert int(r['amplicon_end_in_output_1based'])==right-start+1
        fragment=seq[start-s['source_window_start']:end-s['source_window_start']+1]
        assert fragment==allrefs[rid]
        assert hashlib.sha256(fragment.encode('ascii')).hexdigest()==r['output_sequence_sha256']
        q=urllib.parse.parse_qs(urllib.parse.urlparse(s['source_fetch_url']).query)
        assert q['id']==[s['source_accession_version']]
        assert int(q['seq_start'][0])==s['source_window_start'] and int(q['seq_stop'][0])==s['source_window_end']
        if rid in primary:
            assert int(r['left_flank_nt'])==int(r['right_flank_nt'])==100
            assert r['selected_sites_compatible']=='true'
            assert r['selected_sites_unambiguous_ACGT']=='true'
        byassay[(r['target'],r['assay'])]=(s,r,seq)
        expectedseq[rid]=fragment
    oligos=read_tsv(root/'oligo_alignment_audit.tsv')
    assert len(oligos)==sum(len(s['oligos']) for s in catalog)
    selected={}
    for a in oligos:
        s,r,seq=byassay[(a['target'],a['assay'])]
        assert a['reference_id']==r['reference_id']
        assert any(o['name']==a['oligo_name'] and o['sequence_5to3']==a['sequence_5to3'] and o['direction']==a['direction'] for o in s['oligos'])
        if a['source_start_1based']=='NA':
            assert a['selected_for_representative_reference']=='false';continue
        lo,hi=int(a['source_start_1based']),int(a['source_end_1based'])
        bases=seq[lo-s['source_window_start']:hi-s['source_window_start']+1]
        assert bases==a['reference_bases_plus_strand']
        oriented=bases if a['strand']=='+' else reverse_complement(bases)
        assert oriented==a['reference_bases_in_oligo_direction']
        assert len(oriented)==len(a['sequence_5to3'])
        mm=[i+1 for i,(x,y) in enumerate(zip(a['sequence_5to3'],oriented)) if not set(BASES[x]).intersection(BASES[y])]
        assert len(mm)==int(a['mismatches'])
        assert (','.join(map(str,mm)) or 'NA')==a['mismatch_positions_1based_in_oligo']
        assert int(a['output_start_1based'])==lo-int(r['reference_start_1based'])+1
        assert int(a['output_end_1based'])==hi-int(r['reference_start_1based'])+1
        if a['selected_for_representative_reference']=='true':
            key=(a['target'],a['assay']);selected.setdefault(key,{})[a['direction']]=a
    for key,(s,r,seq) in byassay.items():
        chosen=selected[key]
        assert {'F','R'}<=set(chosen)
        if any(o['direction']=='P' for o in s['oligos']):assert 'P' in chosen
        assert int(chosen['F']['source_start_1based'])==int(r['amplicon_start_1based'])
        assert int(chosen['R']['source_end_1based'])==int(r['amplicon_end_1based'])
        assert int(chosen['F']['source_end_1based'])<int(chosen['R']['source_start_1based'])
        if 'P' in chosen:
            p=chosen['P'];assert int(r['amplicon_start_1based'])<=int(p['source_start_1based'])<=int(p['source_end_1based'])<=int(r['amplicon_end_1based'])
        if r['reference_id'] in primary:assert all(a['mismatches']=='0' for a in chosen.values())
    mapping=read_tsv(root/'target_record_mapping.tsv')
    assert len(mapping)==len(catalog) and {m['reference_id'] for m in mapping}==set(allrefs)
    lit=read_tsv(root/'literature_assays.tsv');assert len(lit)==len(catalog)
    for r in lit:assert r['reference_id'] in allrefs and r['doi'].startswith('10.')
    report=json.loads((root/'validation_summary.json').read_text())
    assert report['primary_record_count']==len(primary)
    assert report['candidate_record_count']==len(allrefs)
    assert report['assay_record_count']==len(catalog)
    assert report['boundary_limited_count']==len(boundary) and report['oligo_mismatch_count']==len(mismatch)
    return dict(result='PASS',target_count=len({s['target'] for s in catalog}),assay_count=len(catalog),
                unique_reference_count=len(allrefs),primary_reference_count=len(primary),boundary_reference_count=len(boundary),
                mismatch_reference_count=len(mismatch),oligo_count=len(oligos),metadata_column_count=len(COLUMNS),
                checks=['fixed_accession_versions','source_checksums','unchanged_source_substrings','one_based_inclusive_coordinates',
                        'forward_reverse_orientation','oligo_sequence_and_mismatch_counts','primary_full100_flanks',
                        'exact_sequence_aliases','FASTA_metadata_lengths_and_identifiers','NCBI_taxonomy_fields',
                        'disjoint_output_groups','literature_assay_mapping'])

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--package',type=Path,required=True)
    args=parser.parse_args();print(json.dumps(validate_package(args.package.resolve()),indent=2))
