#!/usr/bin/env python3
"""RefSeq first-exon transcript-start proximity for flagged DHSs."""
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data/refseq'
OUT=ROOT/'results/wtc11/refseq_proximity'
RELEASE='GCF_000001405.40-RS_2025_08'
ALLOWED_BIOTYPES={'protein_coding', 'transcribed_pseudogene'}
URL=f'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/{RELEASE}/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz'

def tss(exons,strand):
    """First exonic base in transcript orientation; 0-based coordinate."""
    e=sorted(set(exons))
    return e[0][0] if strand=='+' else e[-1][1]-1

def distance(a,b,pos):
    return max(a-pos,pos-(b-1),0)

def table(d):
    return '\n'.join(['| '+' | '.join(d.columns)+' |','| '+' | '.join(['---']*len(d.columns))+' |']+
                     ['| '+' | '.join(str(v) for v in row)+' |' for row in d.itertuples(index=False,name=None)])

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    # Strand, BED boundaries and single-exon regression checks.
    assert tss([(100,150),(200,250)],'+')==100
    assert tss([(100,150),(200,250)],'-')==249
    assert tss([(100,150)],'+')==100 and tss([(100,150)],'-')==149
    assert distance(100,200,199)==0 and distance(100,200,200)==1
    assert distance(100,200,1100)==901 and distance(100,200,1199)==1000
    path=DATA/f'{RELEASE}_genomic.gtf.gz'
    digest=hashlib.md5(path.read_bytes()).hexdigest()
    lines=(DATA/'md5checksums.txt').read_text().splitlines()
    expected=[x.split()[0] for x in lines if x.split()[-1]=='./GCF_000001405.40_GRCh38.p14_genomic.gtf.gz']
    assert expected==[digest],f'MD5 mismatch: {expected} vs {digest}'
    exons=defaultdict(list); meta={}; headers=[]; gene_biotypes={}
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith('#'):
                headers.append(line.strip());continue
            v=line.rstrip().split('\t')
            if v[0]!='NC_000006.12':continue
            attrs=dict(re.findall(r'(\w+) "([^"]*)"',v[8]))
            if v[2]=='gene':
                gene_biotypes[attrs['gene_id']]=attrs.get('gene_biotype','unknown')
                continue
            if v[2]!='exon':continue
            tx=attrs.get('transcript_id','')
            if not tx:continue
            key=(attrs['gene_id'],tx,v[6])
            exons[key].append((int(v[3])-1,int(v[4])))
            geneid=re.search(r'GeneID:(\d+)',v[8])
            meta[key]=dict(gene=attrs.get('gene',attrs['gene_id']),gene_id=geneid.group(1) if geneid else '',
                           transcript=tx,strand=v[6],source=v[1])
    rows=[]
    for key,e in exons.items():
        assert key[0] in gene_biotypes, f'Missing gene biotype: {key[0]}'
        if gene_biotypes[key[0]] not in ALLOWED_BIOTYPES:continue
        e=sorted(set(e)); m=meta[key]
        rows.append(dict(**m,chrom='chr6',gene_biotype=gene_biotypes[key[0]],exon_count=len(e),annotated_tss0=tss(e,m['strand'])))
    tx=pd.DataFrame(rows)
    tx.to_csv(OUT/'refseq_chr6_transcript_sites.tsv.gz',sep='\t',index=False)
    dhs=pd.read_csv(ROOT/'results/wtc11/flagged_hit_DHSs.tsv',sep='\t')
    assert len(dhs)==17 and dhs.dhs.is_unique
    matches=[]; summaries=[]
    for r in dhs.itertuples():
        chrom,span=r.dhs.split(':'); a,b=map(int,span.split('-'))
        item=dict(dhs=r.dhs,library_genomic_element=r.genomic_element,screen_response_genes=r.genes)
        for kind,col in [('annotated_tss','annotated_tss0')]:
            valid=tx[tx[col].notna()].copy()
            valid['site0']=valid[col].astype(int)
            valid['distance_bp']=[distance(a,b,pos) for pos in valid.site0]
            near=valid[valid.distance_bp.le(1000)]
            nearest=valid.sort_values(['distance_bp','gene','transcript']).iloc[0]
            item[f'{kind}_within_1kb']=len(near)>0
            item[f'{kind}_genes']=';'.join(sorted(near.gene.unique()))
            item[f'{kind}_n_transcripts']=len(near)
            item[f'{kind}_nearest_gene']=nearest.gene
            item[f'{kind}_nearest_distance_bp']=int(nearest.distance_bp)
            for n in near.itertuples():
                matches.append(dict(dhs=r.dhs,definition=kind,gene=n.gene,gene_id=n.gene_id,transcript=n.transcript,
                    source=n.source,gene_biotype=n.gene_biotype,strand=n.strand,exon_count=n.exon_count,site0=n.site0,site1=n.site0+1,
                    distance_bp=n.distance_bp,window_start0=max(0,n.site0-1000),window_end0=n.site0+1001))
        summaries.append(item)
    s=pd.DataFrame(summaries); m=pd.DataFrame(matches)
    s.to_csv(OUT/'flagged_DHS_refseq_proximity.tsv',sep='\t',index=False)
    m.to_csv(OUT/'flagged_DHS_refseq_transcript_matches.tsv',sep='\t',index=False)
    info=dict(release=RELEASE,url=URL,md5=digest,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
              gene_biotype_filter=sorted(ALLOWED_BIOTYPES),gtf_headers=headers,chr6_transcripts=len(tx),single_exon_transcripts=int(tx.exon_count.eq(1).sum()),
              flagged_DHSs=len(s),annotated_tss_near_DHSs=int(s.annotated_tss_within_1kb.sum()),distance_rule='any DHS base within inclusive +/-1000 bp of site')
    (OUT/'provenance.json').write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps(info,indent=2))
    print(s.to_string(index=False))
    report=f"""# RefSeq TSS proximity of flagged WTC11 DHSs

**{info['annotated_tss_near_DHSs']} of 17 flagged DHSs ({100*info['annotated_tss_near_DHSs']/len(s):.1f}%) overlap a ±1-kb window
around at least one protein-coding gene or transcribed pseudogene's RefSeq TSS.** The remaining
{len(s)-info['annotated_tss_near_DHSs']} do not meet this positional
criterion. All 17 retain their original library `enhancer` labels in a separate
column; those labels are not a reliable substitute for this TSS-distance analysis.

## Definition and annotation

TSS means the **first base of the first exon in transcript orientation**:
lowest exonic coordinate on the positive strand, highest exonic coordinate on
the negative strand. Single-exon transcripts are included. Any base of the DHS
within an inclusive ±1,000 bp of any TSS qualifies; the DHS midpoint is not used.
A distance of zero means the TSS lies inside the DHS. Each DHS counts once,
regardless of how many transcripts match.

The annotation is **{RELEASE} on GRCh38.p14**, the latest GRCh38 release listed
in NCBI's live directory when checked September 10, 2026. All curated and
predicted RefSeq transcript isoforms with exon records belonging to genes with
RefSeq `gene_biotype=protein_coding` or `transcribed_pseudogene` are included on primary chr6 (NC_000006.12).
The gene-level filter applies to both nearest-TSS searches and ±1-kb calls.
All isoforms of those genes are retained, including noncoding transcript variants;
Other noncoding gene classes and other pseudogene classes are excluded.
“Transcribed pseudogene” follows the explicit RefSeq gene biotype, rather than
inferring transcription from any pseudogene with an exon record. Alternative contigs are excluded.
Genes without an annotated transcript/exon structure cannot contribute a TSS.

The nearby RefSeq gene is not necessarily the screen response gene. This is
a positional promoter-proximity label, not proof of a promoter's function or
of the DHS's regulatory target. In particular, a DHS more than 1 kb from every
included TSS is not thereby proven to be a functional distal enhancer.

## Per-DHS results

"""
    view=s[['dhs','screen_response_genes','annotated_tss_genes','annotated_tss_nearest_gene','annotated_tss_nearest_distance_bp']].copy()
    view=view.rename(columns={'dhs':'DHS','screen_response_genes':'Screen response genes',
        'annotated_tss_genes':'RefSeq genes within 1 kb','annotated_tss_nearest_gene':'Nearest TSS gene',
        'annotated_tss_nearest_distance_bp':'Nearest TSS distance (bp)'})
    report+=table(view.replace('', 'None'))
    report+=f"""

## Provenance, validation and reproduction

[NCBI release directory](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/{RELEASE}/)
provides the [GTF]({URL}). Its header states annotation date 08/01/2025 and
GRCh38.p14. Its MD5 matches NCBI's checksum file. Saved directory listings
record the available release versions at retrieval.

Retained {len(tx):,} chr6 transcripts of protein-coding genes and transcribed pseudogenes with exons, including
{info['single_exon_transcripts']:,} single-exon transcripts. Synthetic checks
cover both strands, single-exon transcripts, DHS interval boundaries and
inclusive 1-kb distance. Complete transcript-level matches preserve RefSeq
accessions, gene IDs, strands and exact site coordinates for inspection.

Reproduce with `mamba activate mambaforge` followed by
`python annotate_flagged_dhs_refseq.py`. Inputs are in `data/refseq/`.
Outputs under `results/wtc11/refseq_proximity/`:

- `flagged_DHS_refseq_proximity.tsv`: all 17 DHSs, proximity calls, matching
  genes/transcript counts and nearest-TSS distances.
- `flagged_DHS_refseq_transcript_matches.tsv`: all matching transcript accessions,
  strand, 0-based and 1-based TSS coordinates, and distances.
- `refseq_chr6_transcript_sites.tsv.gz`: derived chr6 transcript starts for the two included gene biotypes.
- `provenance.json`: release, source URL, hashes, GTF header and counts.
"""
    (ROOT/'WTC11_REFSEQ_PROXIMITY.md').write_text(report)

if __name__=='__main__':main()
