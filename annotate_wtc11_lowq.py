#!/usr/bin/env python3
"""Annotate exact and reference-projected WTC11 guide sites with hifiasm lowQ BEDs.

Run after fetch_wtc11_lowq.py and minimap2 -x asm20 -c --eqx per haplotype,
using the local mhc_ref.fa as query and the WTC11 haplotype FASTA as target.
Unmapped, deleted, ambiguous or low-MAPQ projections remain unassessed.
"""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd
from run_wtc11_analysis import read_hits, digest

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/wtc11/lowq'
DATA=ROOT/'data/wtc11_lowq'

def fasta(path):
    with open(path) as f:
        name=f.readline().strip()[1:]
        seq=''.join(line.strip() for line in f).upper()
    return name,seq

def overlap(a,b,intervals):
    return any(a < end and b > start for start,end in intervals)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    name,ref=fasta(ROOT/'data/mhc_ref.fa')
    nominal=int(name.split(':')[1].split('-')[0])
    summary=pd.read_csv(ROOT/'results/wtc11/guide_exact_match_summary_trim1bp.tsv.gz',sep='\t')
    gr=read_hits(ROOT/'data/grch38_exact_alignments_trim1bp.tsv.gz')
    gr=gr[gr.contig.eq('chr6') & gr.start0.ge(nominal) & gr.end0.le(nominal+len(ref))].copy()
    seqs=summary.set_index('spacer_id').trimmed_spacer.to_dict()
    rc=str.maketrans('ACGT','TGCA')
    # Resolve the FASTA header's coordinate convention against known exact hits.
    checks={}
    for offset in [nominal,nominal-1]:
        matches=0
        for r in gr.itertuples():
            seq=seqs[r.spacer_id]
            if r.strand=='-': seq=seq.translate(rc)[::-1]
            matches += ref[r.start0-offset:r.end0-offset]==seq
        checks[offset]=matches
    offset=max(checks,key=checks.get)
    assert checks[offset]==len(gr),f'Reference sequence does not reproduce cached GRCh38 hits: {checks}'
    print('Reference offset:',offset,'verified records:',len(gr),flush=True)
    coords=summary[['guide_id','spacer_id','guide_start','intended_target_name','WTC11_MHC_exact_match_count']].merge(gr[['spacer_id','start0','end0']],on='spacer_id')
    coords['distance']=(coords.start0-pd.to_numeric(coords.guide_start,errors='coerce')).abs()
    coords=coords.sort_values(['guide_id','distance','start0']).drop_duplicates('guide_id')
    records=[]; provenance=[]; bed_stats=[]
    members=json.loads((DATA/'archive_members.json').read_text())
    for hap in ['hap1','hap2']:
        ctg,seq=fasta(ROOT/f'data/wtc11.{hap}.mhc_contigs.fa')
        original=ctg.removeprefix(hap+'_')
        gfa=ROOT/f'{hap}.gfa'
        archived=next(m for m in members if m['name']==f'asm/assembly.asm.bp.{hap}.p_ctg.gfa')
        assert gfa.stat().st_size==archived['size']
        same=False
        with open(gfa) as f:
            for line in f:
                if line.startswith('S\t'+original+'\t'):
                    same=line.split('\t',3)[2].strip().upper()==seq
                    break
        assert same,'Selected FASTA differs from local GFA'
        bed=pd.read_csv(DATA/f'assembly.asm.bp.{hap}.p_ctg.lowQ.bed',sep='\t',header=None,usecols=[0,1,2],names=['contig','start0','end0'])
        bed=bed[bed.contig.eq(original)].copy()
        assert (bed.start0.ge(0)&bed.end0.le(len(seq))&bed.end0.gt(bed.start0)).all()
        bed.to_csv(OUT/f'{hap}_MHC_contig_lowQ.bed',sep='\t',header=False,index=False)
        intervals=list(bed[['start0','end0']].itertuples(index=False,name=None))
        merged=[]
        for a,b in sorted(intervals):
            if merged and a<=merged[-1][1]: merged[-1][1]=max(merged[-1][1],b)
            else: merged.append([a,b])
        bed_stats.append(dict(haplotype=hap,contig=ctg,contig_bp=len(seq),lowq_intervals=len(bed),lowq_union_bp=sum(b-a for a,b in merged)))
        mapping=np.full(len(ref),-1,dtype=np.int64)
        ambiguous=np.zeros(len(ref),dtype=bool)
        with open(OUT/f'{hap}_reference_to_hap.paf') as f:
            for line in f:
                v=line.rstrip().split('\t'); tags={t.split(':',2)[0]:t.split(':',2)[2] for t in v[12:]}
                if int(v[11])<20 or tags.get('tp')!='P': continue
                assert v[5]==ctg
                strand=v[4]; q=int(v[2]) if strand=='+' else int(v[3]); t=int(v[7])
                for length,op in re.findall(r'(\d+)([MID=X])',tags['cg']):
                    n=int(length)
                    if op in 'M=X':
                        a,b=(q,q+n) if strand=='+' else (q-n,q)
                        values=np.arange(t,t+n) if strand=='+' else np.arange(t+n-1,t-1,-1)
                        if op=='=':
                            expected=ref[a:b] if strand=='+' else ref[a:b].translate(rc)[::-1]
                            assert expected==seq[t:t+n], 'PAF exact block disagrees with FASTAs'
                        old=mapping[a:b]
                        ambiguous[a:b] |= (old>=0)&(old!=values)
                        mapping[a:b]=values
                    if op in 'MI=X': q+=n if strand=='+' else -n
                    if op in 'MD=X': t+=n
                assert q==(int(v[3]) if strand=='+' else int(v[2])) and t==int(v[8])
        hits=read_hits(ROOT/f'results/wtc11/wtc11_{hap}_exact_alignments_trim1bp.tsv.gz')
        grouped={sid:g for sid,g in hits.groupby('spacer_id')}
        consistent=checked=0
        for r in coords.itertuples():
            a,b=r.start0-offset,r.end0-offset
            vals=mapping[a:b]
            status='projected'
            if ambiguous[a:b].any(): status='ambiguous_projection'
            elif len(vals)!=19 or (vals<0).any(): status='unmapped_or_deleted_bases'
            elif len(vals)>1 and not (np.diff(vals)>0).all() and not (np.diff(vals)<0).all(): status='noncolinear_projection'
            start=end=lowq=distance=None
            if status=='projected':
                start,end=int(vals.min()),int(vals.max())+1
                if end-start!=19: status='projected_spanning_indel'
                lowq=overlap(start,end,merged)
                distance=min((max(lo-end,start-hi,0) for lo,hi in merged),default=None)
            h=grouped.get(r.spacer_id)
            count=0 if h is None else len(h)
            exact_lowq=None if h is None else any(overlap(int(z.start0),int(z.end0),merged) for z in h.itertuples())
            if count==1 and status=='projected':
                checked+=1
                z=h.iloc[0]
                consistent+=int(start==z.start0 and end==z.end0)
            records.append(dict(guide_id=r.guide_id,spacer_id=r.spacer_id,dhs=r.intended_target_name,
                grch38_start0=r.start0,grch38_end0=r.end0,haplotype=hap,contig=ctg,
                exact_match_count=count,missing_both=r.WTC11_MHC_exact_match_count==0,
                projection_status=status,projected_start0=start,projected_end0=end,
                hifiasm_lowQ_overlap=lowq,distance_to_lowQ_bp=distance,any_exact_hit_lowQ_overlap=exact_lowq))
        provenance.append(dict(haplotype=hap,selected_sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
            local_gfa_sequence_identical=True,local_gfa_size_matches_archive=True,
            projection_checked_unique_exact_guides=checked,projection_agrees_with_exact_site=consistent,
            bed_sha256=digest(DATA/f'assembly.asm.bp.{hap}.p_ctg.lowQ.bed'),paf_sha256=digest(OUT/f'{hap}_reference_to_hap.paf')))
        print(hap,bed_stats[-1], 'unique exact projection agreement',consistent,'/',checked,flush=True)
    d=pd.DataFrame(records)
    d.to_csv(OUT/'guide_lowQ_annotations.tsv.gz',sep='\t',index=False)
    pd.DataFrame(bed_stats).to_csv(OUT/'lowQ_contig_summary.tsv',sep='\t',index=False)
    assessed=d.hifiasm_lowQ_overlap.notna()
    stats=[]
    for (hap,missing),g in d.groupby(['haplotype','missing_both']):
        a=g.hifiasm_lowQ_overlap.notna()
        stats.append(dict(haplotype=hap,missing_both=missing,guides=len(g),assessed=int(a.sum()),
                          lowQ_overlap=int(g.hifiasm_lowQ_overlap.eq(True).sum()),unassessed=int((~a).sum())))
    pd.DataFrame(stats).to_csv(OUT/'lowQ_missing_vs_retained.tsv',sep='\t',index=False)
    flagged=pd.read_csv(ROOT/'results/wtc11/flagged_hit_unique_guides.tsv',sep='\t')
    flagged.merge(d,on=['guide_id','dhs'],validate='one_to_many').to_csv(OUT/'flagged_hit_guides_lowQ.tsv',sep='\t',index=False)
    d[d.dhs.eq('chr6:29934071-29934293')].to_csv(OUT/'chr6_29934071_29934293_lowQ.tsv',sep='\t',index=False)
    (OUT/'annotation_provenance.json').write_text(json.dumps(dict(reference_offset=offset,reference_sha256=digest(ROOT/'data/mhc_ref.fa'),
        reference_verified_exact_records=len(gr),minimum_paf_mapq=20,primary_paf_only=True,haplotypes=provenance),indent=2)+'\n')
    print(pd.DataFrame(stats).to_string(index=False))
    selected=d[d.guide_id.isin(flagged.loc[flagged.dhs.eq('chr6:29934071-29934293'),'guide_id'])]
    print(selected[['guide_id','haplotype','projection_status','projected_start0','projected_end0','hifiasm_lowQ_overlap','distance_to_lowQ_bp']].to_string(index=False))

if __name__=='__main__': main()
