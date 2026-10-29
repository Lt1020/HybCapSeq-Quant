#!/usr/bin/env python3
"""Build an offline, provenance-linked EsViritu diagnostic-reference increment.

Only unchanged subsequences of the supplied public reference windows are emitted.
The output separates records that meet the stated sequence checks from exceptions.
"""
import argparse
import collections
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
import sys
import textwrap
import urllib.parse
import xml.etree.ElementTree as ET

COLS = ['Accession','description','Name','Segment','kingdom','phylum','tclass',
        'order','family','genus','species','subspecies','Length','TaxID','Assembly','Asm_length']
IUPAC = {k:set(v) for k,v in zip('ACGTRYSWKMBDHVN',[
    'A','C','G','T','AG','CT','GC','AT','GT','AC','CGT','AGT','ACT','ACG','ACGT'])}
IUPAC['I']=set('ACGT')  # Published inosine, for oligo placement only.
COMP = str.maketrans('ACGTRYSWKMBDHVNI','TGCAYRSWMKVHDBNI')

def require(condition, message):
    if not condition: raise ValueError(message)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def rc(sequence):
    return sequence.translate(COMP)[::-1]

def read_fasta(path):
    result={}; key=None
    for line in path.read_text(encoding='ascii').splitlines():
        if line.startswith('>'):
            key=line[1:].split()[0]
            require(key not in result,'Duplicate FASTA ID: '+key)
            result[key]=''
        elif line.strip():
            require(key is not None,'Sequence before FASTA header')
            result[key]+=line.strip().upper()
    require(result and all(result.values()),'Empty FASTA')
    require(all(set(s)<=set(IUPAC)-{'I'} for s in result.values()),'Invalid reference nucleotide alphabet')
    return result

def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)

def mismatch_count(sequence, pattern):
    require(len(sequence)==len(pattern),'Oligo/reference length mismatch')
    return sum(not bool(IUPAC[a]&IUPAC[b]) for a,b in zip(sequence,pattern))

def locate(window, oligo, direction):
    patterns=[('+',oligo)] if direction=='F' else [('-',rc(oligo))]
    if direction=='P':patterns=[('+',oligo),('-',rc(oligo))]
    hits=[]
    for strand,pattern in patterns:
        for i in range(len(window)-len(pattern)+1):
            ref=window[i:i+len(pattern)]
            n=mismatch_count(ref,pattern)
            if n<=3:
                hits.append(dict(local_start=i+1,local_end=i+len(pattern),strand=strand,
                                 mismatches=n,reference_bases=ref,pattern=pattern))
    return sorted(hits,key=lambda x:(x['mismatches'],x['local_start'],x['strand']))

def taxonomy_map(path):
    out={}
    for node in ET.parse(path).getroot().findall('Taxon'):
        fields={x.findtext('Rank'):x.findtext('ScientificName') for x in node.findall('./LineageEx/Taxon')}
        fields[node.findtext('Rank')]=node.findtext('ScientificName')
        out[int(node.findtext('TaxId'))]=fields
    return out

def write_dataset(directory, records):
    directory.mkdir(parents=True,exist_ok=True)
    fasta=directory/'virus_pathogen_database.fna'
    with fasta.open('w',encoding='ascii',newline='\n') as stream:
        for r in records:
            stream.write('>'+r['reference_id']+'\n')
            for block in textwrap.wrap(r['sequence'],80):stream.write(block+'\n')
    write_tsv(directory/'virus_pathogen_database.all_metadata.tsv',[r['metadata'] for r in records],COLS)
    if records:
        check=read_fasta(fasta)
        require(set(check)=={r['reference_id'] for r in records},'FASTA/metadata IDs differ')
        require(len(check)==len(records),'Duplicate emitted identifier')
        with (directory/'virus_pathogen_database.all_metadata.tsv').open(encoding='utf-8',newline='') as stream:
            reader=csv.DictReader(stream,delimiter='\t')
            require(reader.fieldnames==COLS,'Metadata schema changed')
            for row in reader:
                length=len(check[row['Accession']])
                require(int(row['Length'])==length==int(row['Asm_length']),'Metadata length mismatch')
                require(row['Assembly']==row['Name'],'Independent grouping mismatch')

def build(inputs, output):
    require(not output.exists(),'Output directory already exists')
    output.mkdir(parents=True)
    sources=output/'sources';sources.mkdir()
    catalog=json.loads((inputs/'source_catalog.json').read_text())
    tax=taxonomy_map(inputs/'taxonomy_snapshot.xml')
    records=[]; alignments=[];issues=[]
    for spec in catalog:
        src=inputs/'source_windows'/spec['source_file']
        require(sha(src)==spec['source_fasta_sha256'],'Source checksum mismatch')
        parsed=read_fasta(src)
        require(len(parsed)==1,'One source window per record is required')
        source_header,window=next(iter(parsed.items()))
        accession=spec['source_accession_version']
        require(source_header.startswith(accession),'Source accession/version mismatch')
        offset=spec['source_window_start']-1
        require(len(window)==spec['source_window_end']-offset,'Source window coordinate mismatch')
        selected={};oligo_rows=[]
        for o in spec['oligos']:
            hits=locate(window,o['sequence_5to3'],o['direction'])
            if hits:
                h=hits[0]
                a=dict(target=spec['target'],assay=spec['assay'],oligo_name=o['name'],
                       direction=o['direction'],sequence_5to3=o['sequence_5to3'],
                       source_accession_version=accession,source_start_1based=h['local_start']+offset,
                       source_end_1based=h['local_end']+offset,strand=h['strand'],
                       mismatches=h['mismatches'],reference_bases_plus_strand=h['reference_bases'],
                       match_status='compatible' if h['mismatches']==0 else 'mismatch',
                       uniqueness_within_window=sum(v['mismatches']==h['mismatches'] for v in hits))
                oriented=h['reference_bases'] if h['strand']=='+' else rc(h['reference_bases'])
                a['reference_bases_in_oligo_direction']=oriented
                a['mismatch_positions_1based_in_oligo']=','.join(str(i+1) for i,(x,y) in enumerate(zip(o['sequence_5to3'],oriented)) if not IUPAC[x]&IUPAC[y]) or 'NA'
                key=(h['mismatches'],h['local_start'],o['name'])
                if o['direction'] not in selected or key<selected[o['direction']][0]:
                    selected[o['direction']]=(key,a,h)
            else:
                a=dict(target=spec['target'],assay=spec['assay'],oligo_name=o['name'],
                       direction=o['direction'],sequence_5to3=o['sequence_5to3'],
                       source_accession_version=accession,source_start_1based='NA',source_end_1based='NA',
                       strand='NA',mismatches='NA',reference_bases_plus_strand='NA',
                       match_status='no_site_within_3_mismatches_in_bounded_window',uniqueness_within_window=0)
                a['reference_bases_in_oligo_direction']='NA'
                a['mismatch_positions_1based_in_oligo']='NA'
            oligo_rows.append(a)
        classes=['F','R']+(['P'] if any(o['direction']=='P' for o in spec['oligos']) else [])
        require(set(selected)==set(classes),'No bounded candidate site for an oligo class')
        chosen=[selected[k][1] for k in classes]
        f,r=chosen[:2]
        p=selected['P'][1] if 'P' in selected else None
        left,right=f['source_start_1based'],r['source_end_1based']
        require(f['source_end_1based']<r['source_start_1based'],'Primers do not form a convergent pair')
        if p:require(left<=p['source_start_1based']<=p['source_end_1based']<=right,'Probe outside amplicon')
        start=max(1,left-100);end=min(spec['source_length_nt'],right+100)
        require(spec['source_window_start']<=start<=end<=spec['source_window_end'],'Source window lacks required available flanks')
        sequence=window[start-1-offset:end-offset]
        require(len(sequence)==end-start+1,'Extracted length mismatch')
        lf,rf=left-start,end-right
        full=lf==100 and rf==100
        compatible=all(x['mismatches']==0 for x in chosen)
        concrete=all(set(x['reference_bases_plus_strand'])<=set('ACGT') for x in chosen)
        if not full:status='boundary_limited'
        elif not compatible or not concrete:status='oligo_mismatch_review_required'
        else:status='full100_compatible_sites'
        require(all(x['uniqueness_within_window']==1 for x in chosen),'Ambiguous selected oligo site')
        rid='PCRROI_'+spec['target']+'_'+spec['assay']+'_'+accession+'_'+str(start)+'_'+str(end)
        require(re.fullmatch(r'[A-Za-z0-9_.-]+',rid) is not None,'Invalid reference ID')
        name='PCRROI_'+spec['target']+'_'+spec['assay']
        lineage=tax[spec['taxid']]
        metadata={c:'NA' for c in COLS}
        for rank in ['kingdom','phylum','order','family','genus','species','subspecies']:
            metadata[rank]=lineage.get(rank,'NA')
        metadata.update(Accession=rid,Name=name,Assembly=name,Length=len(sequence),Asm_length=len(sequence),
                        TaxID=spec['taxid'],Segment=spec['segment'],tclass=lineage.get('class','NA'),
                        description=f"Literature diagnostic reference window; target={spec['target']}; source={accession}:{start}-{end}; assay={spec['assay']}; status={status}")
        warnings=[]
        if not full:warnings.append(f'available_flanks_left={lf}_right={rf}; no_padding')
        if not compatible:warnings.append('selected_site_mismatches='+','.join(x['oligo_name']+':'+str(x['mismatches']) for x in chosen if x['mismatches']))
        if spec['target']=='Influenza_B' and spec['assay']=='Chen2011':warnings.append('paper_labels_target_NP_but_accession_is_segment_8_NS1_NEP')
        if spec['target']=='SARS-CoV-2' and spec['assay'].startswith('CDC2020'):warnings.append('actual_sequence_match_coordinates_are_16_nt_before_paper_table_coordinates')
        if spec['target']=='Influenza_A' and spec['assay']=='Chen2011':warnings.append('forward_candidate_at_146_167_differs_from_paper_156_177_and_has_1_mismatch')
        if spec['assay']=='Johnson2005':warnings.append('actual_sequence_match_coordinates_are_17_nt_before_paper_table_coordinates')
        if spec['assay']=='Ward2004':warnings.append('NCBI_accession_strain_name_differs_from_original_article_label')
        if spec['assay']=='Li2010_G2R_G':warnings.append('inverted_terminal_repeat_target_one_physical_locus_retained_not_copy_number_corrected')
        if any('I' in o['sequence_5to3'] for o in spec['oligos']):warnings.append('inosine_is_a_location_wildcard_not_a_reference_DNA_base_or_binding_efficiency_claim')
        if spec['target'] in ['Human_astrovirus','Sapovirus','Norovirus_GI','Rotavirus_A']:
            warnings.append('representative_reference_does_not_require_every_mixture_oligo_to_match')
        region=dict(reference_id=rid,target=spec['target'],assay=spec['assay'],source_accession_version=accession,
                    source_title=spec['source_title'],source_length_nt=spec['source_length_nt'],source_taxid=spec['taxid'],
                    reference_start_1based=start,reference_end_1based=end,reference_length_nt=len(sequence),
                    amplicon_start_1based=left,amplicon_end_1based=right,amplicon_length_nt=right-left+1,
                    amplicon_start_in_output_1based=lf+1,amplicon_end_in_output_1based=len(sequence)-rf,
                    left_flank_nt=lf,right_flank_nt=rf,full_100bp_flanks=str(full).lower(),
                    forward_oligo=f['oligo_name'],reverse_oligo=r['oligo_name'],probe_oligo=p['oligo_name'] if p else 'NA_SYBR',
                    forward_mismatches=f['mismatches'],reverse_mismatches=r['mismatches'],probe_mismatches=p['mismatches'] if p else 'NA_SYBR',
                    selected_sites_compatible=str(compatible).lower(),selected_sites_unambiguous_ACGT=str(concrete).lower(),
                    sequence_check_status=status,actual_commercial_assay_confirmed='false',
                    paper_amplicon_start_1based=spec.get('paper_amplicon_start') or 'NA',
                    paper_amplicon_end_1based=spec.get('paper_amplicon_end') or 'NA',doi=spec['doi'],
                    oligo_sequence_source_url=spec['oligo_sequence_source_url'],
                    source_record_url='https://www.ncbi.nlm.nih.gov/nuccore/'+accession,
                    source_fetch_url=spec['source_fetch_url'],retrieved_utc=spec['retrieved_utc'],
                    source_window_start_1based=spec['source_window_start'],source_window_end_1based=spec['source_window_end'],
                    source_file='sources/source_windows/'+spec['source_file'],source_fasta_sha256=spec['source_fasta_sha256'],
                    output_sequence_sha256=hashlib.sha256(sequence.encode('ascii')).hexdigest(),notes='; '.join(warnings),
                    assay_type=spec['assay_type'],gene=spec.get('gene','NA'),
                    design_reference_accessions=spec['design_reference_accessions'],published_reference_accessions=spec['published_reference_accessions'],
                    reference_role=spec['reference_role'],reference_evidence_url=spec['reference_evidence_url'],literature_notes=spec['literature_notes'])
        for a in oligo_rows:
            a['reference_id']=rid
            a['selected_for_representative_reference']=str(a['oligo_name'] in {x['oligo_name'] for x in chosen}).lower()
            if a['source_start_1based']!='NA':
                a['output_start_1based']=a['source_start_1based']-start+1
                a['output_end_1based']=a['source_end_1based']-start+1
            else:a['output_start_1based']='NA';a['output_end_1based']='NA'
            a['doi']=spec['doi'];a['source_url']=spec['oligo_sequence_source_url']
        records.append(dict(reference_id=rid,sequence=sequence,metadata=metadata,region=region,status=status))
        alignments.extend(oligo_rows)
        for warning in warnings:issues.append(dict(reference_id=rid,target=spec['target'],status=status,issue=warning))

    require(len({r['reference_id'] for r in records})==len(records),'Duplicate reference ID')
    # Exact sequence duplicates for the same target share one FASTA record.
    # Keep every assay and its provenance in the mapping and audit tables.
    canonical={};aliases={};unique=[]
    for rec in records:
        key=(rec['region']['target'],rec['sequence'])
        if key in canonical:
            first=canonical[key]
            require(first['status']==rec['status'],'Duplicate sequence has differing assay QC status; resolve explicitly')
            aliases[rec['reference_id']]=first['reference_id']
            rec['reference_id']=first['reference_id']
            rec['region']['reference_id']=first['reference_id']
            rec['metadata']=first['metadata']
        else:
            canonical[key]=rec;unique.append(rec)
    for a in alignments:a['reference_id']=aliases.get(a['reference_id'],a['reference_id'])
    for a in issues:a['reference_id']=aliases.get(a['reference_id'],a['reference_id'])
    passed=[r for r in unique if r['status']=='full100_compatible_sites']
    boundary=[r for r in unique if r['status']=='boundary_limited']
    mismatch=[r for r in unique if r['status']=='oligo_mismatch_review_required']
    write_dataset(output,passed)
    write_dataset(output/'all_candidates',unique)
    write_dataset(output/'exceptions/boundary_limited',boundary)
    write_dataset(output/'exceptions/oligo_mismatch',mismatch)
    write_tsv(output/'reference_regions.tsv',[r['region'] for r in records],list(records[0]['region']))
    write_tsv(output/'oligo_alignment_audit.tsv',alignments,list(alignments[0]))
    write_tsv(output/'record_issues.tsv',issues,['reference_id','target','status','issue'])
    write_tsv(output/'target_record_mapping.tsv',[
        dict(target=r['region']['target'],assay=r['region']['assay'],reference_id=r['reference_id'],
             Name=r['metadata']['Name'],Assembly=r['metadata']['Assembly'],sequence_check_status=r['status'],
             in_primary_files=str(any(x['reference_id']==r['reference_id'] for x in passed)).lower(),species=r['metadata']['species'],TaxID=r['metadata']['TaxID'])
        for r in records],['target','assay','reference_id','Name','Assembly','sequence_check_status','in_primary_files','species','TaxID'])
    write_tsv(output/'sequence_aliases.tsv',[dict(alias_reference_id=k,canonical_reference_id=v,reason='same_target_exact_sequence_duplicate') for k,v in aliases.items()],['alias_reference_id','canonical_reference_id','reason'])
    bibliography=json.loads((inputs/'bibliography.json').read_text())
    literature=[]
    for spec,rec in zip(catalog,records):
        b=bibliography.get(spec['doi'].lower(),{})
        require(b.get('doi','').lower()==spec['doi'].lower(),'Missing exact-DOI bibliography entry')
        z=rec['region']
        line=dict(target=spec['target'],assay=spec['assay'],assay_type=spec['assay_type'],gene=spec.get('gene','NA'),
                  year=b.get('pubYear','NA'),authors=b.get('authorString','NA'),title=b.get('title','NA'),journal=b.get('journalTitle','NA'),
                  doi=spec['doi'],doi_url='https://doi.org/'+spec['doi'],PMID=b.get('id','NA') if b.get('source')=='MED' else 'NA',
                  public_oligo_source_url=spec['oligo_sequence_source_url'],reference_evidence_url=spec['reference_evidence_url'],
                  design_reference_accessions=spec['design_reference_accessions'],published_reference_accessions=spec['published_reference_accessions'],
                  extraction_accession_version=spec['source_accession_version'],extraction_reference_role=spec['reference_role'],
                  ncbi_record_title=spec['source_title'],reference_id=rec['reference_id'],
                  amplicon_start_1based=z['amplicon_start_1based'],amplicon_end_1based=z['amplicon_end_1based'],
                  amplicon_length_nt=z['amplicon_length_nt'],output_start_1based=z['reference_start_1based'],output_end_1based=z['reference_end_1based'],
                  left_flank_nt=z['left_flank_nt'],right_flank_nt=z['right_flank_nt'],selected_sites_compatible=z['selected_sites_compatible'],
                  sequence_check_status=rec['status'],in_primary_files=str(rec['status']=='full100_compatible_sites').lower(),
                  selection_basis=spec['selection_basis'],notes=spec['literature_notes'])
        for d,label in [('F','forward_oligos_5to3'),('R','reverse_oligos_5to3'),('P','probe_oligos_5to3')]:
            line[label]='; '.join(o['name']+':'+o['sequence_5to3'] for o in spec['oligos'] if o['direction']==d) or 'NA_SYBR'
        literature.append(line)
    write_tsv(output/'literature_assays.tsv',literature,list(literature[0]))
    # Exact records that overlap in one reference accession are retained as
    # alternative diagnostic windows, never asserted to be independent counts.
    overlaps=[]
    for i,a in enumerate(unique):
        for b in unique[i+1:]:
            x,y=a['region'],b['region']
            if x['source_accession_version']!=y['source_accession_version']:continue
            amount=min(x['reference_end_1based'],y['reference_end_1based'])-max(x['reference_start_1based'],y['reference_start_1based'])+1
            if amount>0:overlaps.append(dict(reference_id_a=a['reference_id'],reference_id_b=b['reference_id'],
                accession_version=x['source_accession_version'],overlap_nt=amount,reason='overlapping_reference_windows_do_not_sum_read_counts'))
    write_tsv(output/'overlapping_reference_windows.tsv',overlaps,['reference_id_a','reference_id_b','accession_version','overlap_nt','reason'])
    shutil.copytree(inputs/'source_windows',sources/'source_windows')
    for name in ['source_catalog.json','taxonomy_snapshot.xml','fluB_feature_audit.json','bibliography.json']:
        shutil.copy2(inputs/name,sources/name)
    code=output/'scripts';code.mkdir()
    shutil.copy2(Path(__file__),code/'build_reference_package.py')
    report=dict(package_type='literature_diagnostic_ROI_reference_increment',created_utc=max(s['retrieved_utc'] for s in catalog),
                source_target_count=len({r['region']['target'] for r in records}),assay_record_count=len(records),candidate_record_count=len(unique),exact_sequence_alias_count=len(aliases),
                primary_record_count=len(passed),boundary_limited_count=len(boundary),oligo_mismatch_count=len(mismatch),
                metadata_columns=COLS,coordinate_system='one-based inclusive on accession.version; plus strand',
                exact_match_definition='At least one IUPAC-compatible forward and reverse member and internal probe member (probe omitted for SYBR); selected reference sites unambiguous ACGT. Published inosine is a placement wildcard only.',
                validation=dict(fasta_metadata_ids_and_lengths_match=True,unique_ids_within_package=True,
                                emitted_sequences_are_unchanged_source_substrings=True,no_padding_or_primer_replacement=True,
                                primary_records_have_exactly_100_bases_each_flank=True,source_checksums_verified=True,
                                collision_check_against_full_base_database_performed=False,minimap2_index_built=False,
                                experimental_validation_performed=False,actual_commercial_assay_confirmed=False),
                release_ready=False,notes='Sequence matching is a computational check, not validation of assay performance or quantitative accuracy.')
    (output/'validation_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    summary=[]
    labels={'full100_compatible_sites':'完整100 bp且所选引物/探针兼容',
            'boundary_limited':'参考边界不足',
            'oligo_mismatch_review_required':'所选参考存在引物/探针失配'}
    for target in dict.fromkeys(x['target'] for x in catalog):
        specs=[x for x in catalog if x['target']==target]
        refs=[x for x in unique if x['region']['target']==target]
        summary.append(dict(target=target,assay_count=len(specs),unique_reference_count=len(refs),
            primary_reference_count=sum(x['status']=='full100_compatible_sites' for x in refs),
            boundary_limited_count=sum(x['status']=='boundary_limited' for x in refs),
            mismatch_count=sum(x['status']=='oligo_mismatch_review_required' for x in refs)))
    write_tsv(output/'target_summary.tsv',summary,list(summary[0]))
    readme=f'''# 文献诊断靶区参考片段合集（2026-10-29扩展版）

本包覆盖{report['source_target_count']}个靶标，整理{len(records)}套公开检测体系。本次按用户要求，新冠保留N1、N2、E和ORF1ab共4条；其他靶标收录数及原有序列保持上一版不变。去除{len(aliases)}组完全相同的序列后，共{len(unique)}条真实NCBI参考片段。其中{len(passed)}条满足双侧各100 bp及所选引物/探针兼容条件，放在根目录主FASTA中。{len(boundary)}条边界不足、{len(mismatch)}条存在所选参考上的引物/探针失配，分别列入异常目录。`all_candidates/`包含全部{len(unique)}条片段。Sapovirus Oka2006与Oka2019在AY237420.2上提取得到完全相同的序列，FASTA仅保存一份，各套体系的文献和比对记录全部保留。

## 文件

| 文件/目录 | 内容 |
|---|---|
| `virus_pathogen_database.fna` | {len(passed)}条完整100 bp且所选位点兼容的参考片段 |
| `virus_pathogen_database.all_metadata.tsv` | 与主FASTA逐条对应的16列EsViritu格式metadata |
| `all_candidates/` | 全部{len(unique)}条去重候选及对应metadata，包含下列异常 |
| `exceptions/boundary_limited/` | {len(boundary)}条边界不足的真实片段；其中部分同时存在失配 |
| `exceptions/oligo_mismatch/` | {len(mismatch)}条双侧100 bp齐全但存在所选位点失配的候选 |
| `literature_assays.tsv` | {len(records)}套体系的完整引物/探针、文献、DOI、参考accession及其角色 |
| `reference_regions.tsv` | 各套体系的参考坐标、局部坐标、两端实际延伸长度和QC状态 |
| `oligo_alignment_audit.tsv` | 每条已发表引物/探针在所选参考中的定位、方向、失配数及代表组分 |
| `target_record_mapping.tsv` | 检测体系与FASTA记录的对应关系 |
| `sequence_aliases.tsv` | 完全相同片段的合并记录 |
| `overlapping_reference_windows.tsv` | 同一accession上重叠的不同参考片段 |
| `literature_review_zh.md` | 按靶标组织的检索说明和选取结果 |
| `sources/` | 原始NCBI短片段、固定版本accession、下载URL、分类快照和来源目录 |
| `scripts/` | Python标准库构建、重新下载及独立核验代码 |
| `SHA256SUMS` | 文件校验和 |

## 提取规则

所有坐标均为固定`accession.version`上的**正链、1起始、两端包含**坐标。扩增子外边界定义为正向引物结合位点的最左端a和反向引物结合位点的最右端b；反向引物按反向互补序列定位。提取`[a-100,b+100]`。例如a=203、b=508时，提取103–608，共506 nt；引物已包含在a–b内，不额外叠加一次引物长度。

源记录边界不足时，仅保存`[max(1,a-100),min(source_length,b+100)]`，并从主FASTA中分出。**有accession并不保证该记录覆盖所需的全部延伸范围。** 本包不补N、不拼接其他毒株、不将参考碱基替换为引物序列。全长片段仍逐碱基保持NCBI参考记录原样。

主FASTA的位点兼容条件：同一参考中至少一个已发表正向引物、一个已发表反向引物，以及至少一个内部探针与参考IUPAC兼容且对应参考位点为ACGT。SYBR体系没有探针要求。混合引物中的其他组分可不匹配该代表毒株，逐条结果均在审计表中。已发表脱氧肌苷I仅作为定位通配符，不能据此推断结合效率。失配候选按既有文献靶区附近的最佳比对定位，不标为精确引物匹配记录。

## accession的含义

`design_reference_accessions`只记录能在原始论文中核实的设计模板；未能核实单一设计accession时明确为NA。`published_reference_accessions`记录文献明确给出的坐标、比对或标准物质参考；`reference_role`说明其角色。`source_accession_version`/`extraction_accession_version`是本次实际下载和提取的固定版本。三者不能相互替代。后续应用论文采用的标准物质参考不被写成原始设计模板。

公开文献体系未必等于本次实际数字PCR试剂盒的检测体系；本包没有对商业试剂的引物身份作推断。若用于比较相同物理靶区，应先核对实际数字PCR试剂/引物信息。

## EsViritu字段

metadata严格保持以下16列，使用制表符分隔：

`{' '.join(COLS)}`

FASTA第一段标识与`Accession`完全一致。`Length=Asm_length=实际片段长度`，不是原基因组长度。每个唯一参考片段使用独立的`Name`与`Assembly`，二者相同。`Name`是PCRROI记录标签；实际病毒/分类信息见`target_record_mapping.tsv`及分类列。分类使用NCBI Taxonomy快照，未知字段为NA；`tclass`对应NCBI的class。不同文献体系可能取自同一病毒同一基因区域，不能按记录数统计病毒种数。

本包是短靶区**增量参考文件**，未与v3.3全库合并，也未生成MMI索引。与原有全基因组或其他重叠ROI同时比对会出现多重映射；必须定义去重/归属规则，不能直接把各ROI的read_count或浓度相加。RPKMF若针对ROI计算，应使用ROI有效长度，并与原全基因组RPKMF区分。计算定位核验不能代替数据库特异性、定量性能或实验验证。

## 离线重建

在本包解压目录运行：

```bash
python3 scripts/build_reference_package.py --inputs sources --output rebuilt_reference
python3 scripts/validate_reference_package.py --package rebuilt_reference
```

构建程序只读取本包保存的短参考片段和来源文件，不联网、不设计或修改引物。重新向NCBI核对固定版本源片段，可运行`python3 scripts/refresh_source_windows.py --inputs sources --output refreshed_sources`，随后将`--inputs`改为`refreshed_sources`。下载程序检查原版本和已保存的碱基序列，变化时停止，不自动替换为最新版本。

## 各靶标收录情况

| 靶标 | 文献体系数 | 去重片段数 | 主FASTA | 边界不足 | 失配候选 |
|---|---:|---:|---:|---:|---:|
'''
    for s in summary:readme+='| '+s['target']+' | '+' | '.join(str(s[k]) for k in ['assay_count','unique_reference_count','primary_reference_count','boundary_limited_count','mismatch_count'])+' |\n'
    readme+='''
## 关键核查结果

- SARS-CoV-2 CDC N1/N2：本次在MN908947.3上的实际匹配坐标较原文表1提前16 nt。提取采用真实序列匹配坐标；原文坐标另存，N1实际扩增子长度为72 nt。
- SARS-CoV-2 China CDC ORF1ab：新增公开ORF1ab体系，保留原有N1、N2和E。MN908947.3上的扩增子为13342–13460（119 nt），双侧各延伸100 bp后为13242–13560（319 nt），三个引物/探针位点均精确匹配。MN908947.3作为本次实际定位参考，未写成原作者唯一设计模板。
- DENV-1 Johnson2005：M87512.1上实际区间8956–9067，较原表坐标提前17 nt。
- Influenza B Chen2011：原文将目标标作NP，但CY040733.1是segment 8 NS1/NEP。本包按实际NCBI记录标注。
- Ward2004：原文毒株表述与所引accession的现有记录名称不同；U52940.1实际为A/Texas/36/1991的340 nt部分M1记录，AF100376.1实际为B/Victoria/2/87。实际记录名称、失配及边界不足均已保留。
- Le Cann2004：L06802.1在2484结束，正好到达反向引物外边界，下游没有100 bp。公开复用表的AV2与探针位置与其长度顺序不符，本次按实际匹配位置纠正。
- Rotavirus A：本次保留的3套参考均有一侧不足100 bp，故均位于边界异常目录，主FASTA中未放入RVA记录。Zeng2008与Jothikumar2009虽均为NSP3体系，但位于不同区域。
- MPXV G2R位于倒置末端重复区，本包只保存一个实际位点；基因拷贝数与基因组拷贝数不可直接混同。

“常用”的筛选依据是公开诊断用途、独立原始研究中的复用及可核查程度，未进行全球引用量排名。PMMoV本次保留1套已核验常用定量体系；这不表示相关文献只有1篇。Zhang2006/Haramoto2013同一靶区的历史沿用未重复计数。详细检索取舍见文献说明文件。
'''
    (output/'README_zh.md').write_text(readme,encoding='utf-8')
    review='''# 公开定量检测文献与参考片段核查

检索和NCBI核对日期：2026-10-29。范围为前文12个数字PCR靶标；本次纳入的公开体系为qPCR/RT-qPCR，按用户要求可替代dPCR文献来源。不同体系不是本次商业数字PCR试剂盒身份的证明。序列依据优先采用原始论文表格；原文不能直接取得完整表格时，使用公开的独立原始应用研究复用表，并分别给出原始方法DOI和引物公开来源URL。

本次新冠保留N1、N2、E和新增ORF1ab共4套体系，其他靶标保持上一版收录数量。原版依据公共卫生/污水诊断使用、独立研究复用和序列/参考可追溯性选取，未实施引用量排名，不能称为全球“最常用前三名”。混合引物组分和同一体系的重复引用不单独扩充体系数。

## 收录结果

| 靶标 | 体系 | 实际提取accession.version | 参考区间 | 左/右延伸nt | 主文件/异常 |
|---|---|---|---|---|---|
'''
    for z in literature:
        review+='| '+z['target']+' | '+z['assay']+' | '+z['extraction_accession_version']+' | '+str(z['output_start_1based'])+'–'+str(z['output_end_1based'])+' | '+str(z['left_flank_nt'])+'/'+str(z['right_flank_nt'])+' | '+labels[z['sequence_check_status']]+' |\n'
    review+='\n## 文献及accession角色\n\n'
    for z in literature:
        review+='### '+z['target']+' / '+z['assay']+'\n\n'
        review+=z['authors']+' ('+str(z['year'])+'). '+z['title']+' '+z['journal']+'. DOI: ['+z['doi']+']('+z['doi_url']+').\n\n'
        review+='引物公开来源：'+z['public_oligo_source_url']+'；参考序列证据：'+z['reference_evidence_url']+'。\n\n'
        review+='设计参考：`'+z['design_reference_accessions']+'`；文献所列参考：`'+z['published_reference_accessions']+'`；实际提取：`'+z['extraction_accession_version']+'`；角色：`'+z['extraction_reference_role']+'`。\n\n'
        review+='正向引物（5′→3′）：`'+z['forward_oligos_5to3']+'`。\n\n反向引物（5′→3′）：`'+z['reverse_oligos_5to3']+'`。\n\n探针（5′→3′）：`'+z['probe_oligos_5to3']+'`。\n\n'
        review+=z['notes']+'\n\n'
    review+='''## 检索取舍与适用范围

检索使用靶标名与qPCR、digital PCR、real-time RT-PCR、primer、probe、GenBank/accession等组合，并追溯公开表格的原始方法文献。辅助复用证据包括Mans2014（10.1111/tmi.12251）、CHIKV双体系评估2019（10.3390/v11080755）、Sapovirus原始应用研究2024（10.1128/jvi.00639-24）及RVA的公开原始研究。综述只用来定位原始文献，未将其中未能核实的序列直接作为主库依据。

PMMoV：Haramoto2013对Zhang2006同一位点的正向引物作了修正，重复引用未形成独立靶区。另检到2022年的CAPS/END研究（10.1371/journal.pone.0278784），但其原文说明主要评价存在/不存在，未作为本次常用定量体系补足3条；该文公开表中两个探针相同，本次未进一步将其作为已验证定量参考收录。另一些检索命中为常规终点PCR或特定植物材料中的相对表达研究，本次未据此填满3条。

CHIKV：保留Pastorino2005、Panning2008一般体系和Lanciotti2007的874体系；没有把同一原文中的适配变体或第二套体系再计成新的文献记录。

Sapovirus：保留Oka2006、Chan2006和Oka2019。Oka2006的GI/GII/GV代表株不再作为三篇独立方法计数；本次用GII代表序列。Oka2006与Oka2019输出窗口完全相同，去重后共享一个FASTA标识。

Rotavirus A：3套体系均公开可核查，但所选参考的靶区靠近端点或记录本身只覆盖部分基因。边界不足并不说明方法无效；本包只据此判断其不满足双侧各100 bp的文件要求。其余失配同样只描述文献引物与本次所选固定版本参考的关系，不据此否定原论文的实验结果。

公开文献可能采用多序列比对设计，且并非总是报告单一设计accession。未披露或未核实的原始设计编号保持NA，不能以实际提取accession替代。需要相同实际数字PCR靶区、完整不同基因型覆盖或数据库全局特异性时，应在本合集基础上另行核对，而不能由一个代表片段推出全部覆盖结论。
'''
    (output/'literature_review_zh.md').write_text(review,encoding='utf-8')
    for helper in ['validate_reference_package.py','refresh_source_windows.py']:
        src=Path(__file__).with_name(helper)
        require(src.is_file(),'Missing companion script: '+helper)
        shutil.copy2(src,code/helper)
    validator_spec=importlib.util.spec_from_file_location('independent_roi_validator',code/'validate_reference_package.py')
    validator=importlib.util.module_from_spec(validator_spec)
    previous_bytecode_setting=sys.dont_write_bytecode
    sys.dont_write_bytecode=True
    try:validator_spec.loader.exec_module(validator)
    finally:sys.dont_write_bytecode=previous_bytecode_setting
    independent=validator.validate_package(output)
    (output/'independent_validation.json').write_text(json.dumps(independent,indent=2)+'\n')
    if (inputs/'ncbi_online_validation.json').exists():
        shutil.copy2(inputs/'ncbi_online_validation.json',sources/'ncbi_online_validation.json')
    sums=[sha(p)+'  '+str(p.relative_to(output)) for p in sorted(output.rglob('*')) if p.is_file() and p.name!='SHA256SUMS']
    (output/'SHA256SUMS').write_text('\n'.join(sums)+'\n')
    return report

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True,help='Directory containing source_catalog.json, source_windows/, taxonomy_snapshot.xml, fluB_feature_audit.json')
    parser.add_argument('--output',type=Path,required=True,help='New output directory')
    args=parser.parse_args()
    print(json.dumps(build(args.inputs.resolve(),args.output.resolve()),indent=2))

if __name__=='__main__':main()
