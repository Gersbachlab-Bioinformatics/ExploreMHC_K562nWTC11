#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
src,map_path,grch_path,k562_path,results_dir,plots_dir=sys.argv[1:7]
results_dir,plots_dir=Path(results_dir),Path(plots_dir)
results_dir.mkdir(parents=True,exist_ok=True); plots_dir.mkdir(parents=True,exist_ok=True)
df=pd.read_csv(src,sep='\t',compression='infer',dtype=str).reset_index().rename(columns={'index':'row_index'})
df['spacer']=df['spacer'].fillna('').str.upper()
mp=pd.read_csv(map_path,sep='\t',compression='infer',dtype=str); mp['row_index']=mp['row_index'].astype(int)
df=df.merge(mp[['row_index','spacer_id','original_spacer','trimmed_spacer']],on='row_index',how='left')
valid_mask=df['spacer_id'].notna()
def read_hits(path):
    x=pd.read_csv(path,sep='\t',compression='infer')
    if x.empty:return x
    return x.loc[(x.NM==0)&(~x.secondary.astype(bool))&(~x.supplementary.astype(bool))].copy()
gr,k5=read_hits(grch_path),read_hits(k562_path)
def agg(x,prefix):
    if x.empty:return pd.DataFrame(columns=['spacer_id',f'{prefix}_exact_match_count',f'{prefix}_alignments'])
    x=x.sort_values(['spacer_id','contig','start0','strand']).copy()
    x['formatted']=x.contig.astype(str)+':'+x.start0.astype(str)+'-'+x.end0.astype(str)+':'+x.strand.astype(str)+':'+x.CIGAR.astype(str)+':MAPQ'+x.MAPQ.astype(str)+':NM'+x.NM.astype(str)
    return x.groupby('spacer_id',as_index=False).agg(**{f'{prefix}_exact_match_count':('formatted','size'),f'{prefix}_alignments':('formatted',lambda s:';'.join(s))})
summary=df.merge(agg(gr,'GRCh38'),on='spacer_id',how='left').merge(agg(k5,'K562_hap2_MHC'),on='spacer_id',how='left')
for c in ['GRCh38_exact_match_count','K562_hap2_MHC_exact_match_count']: summary[c]=summary[c].fillna(0).astype(int)
for c in ['GRCh38_alignments','K562_hap2_MHC_alignments']: summary[c]=summary[c].fillna('')
g=summary.GRCh38_exact_match_count>0; k=summary.K562_hap2_MHC_exact_match_count>0
summary['exact_match_category']=np.select([g&k,g&~k,~g&k,~g&~k],['Both','GRCh38 only','K562 MHC only','Neither'],default='Invalid spacer')
summary.to_csv(results_dir/'guide_exact_match_summary_trim1bp.tsv.gz',sep='\t',index=False,compression='gzip')
pd.concat([gr,k5],ignore_index=True).to_csv(results_dir/'all_exact_alignments_trim1bp.tsv.gz',sep='\t',index=False,compression='gzip')
counts=summary.loc[valid_mask].groupby('exact_match_category').size().rename('guide_row_count').reset_index(); counts.to_csv(results_dir/'exact_match_category_counts_trim1bp.tsv',sep='\t',index=False)
membership=mp[['spacer_id','trimmed_spacer']].drop_duplicates().sort_values('spacer_id').reset_index(drop=True)
gset=set(gr.spacer_id) if not gr.empty else set(); kset=set(k5.spacer_id) if not k5.empty else set()
membership['GRCh38_exact']=membership.spacer_id.isin(gset).astype(int); membership['K562_hap2_MHC_exact']=membership.spacer_id.isin(kset).astype(int)
membership['category']=np.select([(membership.GRCh38_exact==1)&(membership.K562_hap2_MHC_exact==1),(membership.GRCh38_exact==1)&(membership.K562_hap2_MHC_exact==0),(membership.GRCh38_exact==0)&(membership.K562_hap2_MHC_exact==1)],['Both','GRCh38 only','K562 MHC only'],default='Neither')
membership.to_csv(results_dir/'unique_trimmed_spacer_membership.tsv.gz',sep='\t',index=False,compression='gzip')
membership[['GRCh38_exact','K562_hap2_MHC_exact']].corr().to_csv(results_dir/'reference_membership_correlation_trim1bp.tsv',sep='\t')
km=membership.loc[membership.K562_hap2_MHC_exact==1].copy()
with open(results_dir/'k562_matching_unique_trimmed_spacers.fa','w') as h:
    for r in km.itertuples(index=False): h.write(f'>{r.spacer_id}\n{r.trimmed_spacer}\n')
if len(km)>=2:
    b={'A':0,'C':1,'G':2,'T':3}; mat=np.array([[b[x] for x in s] for s in km.trimmed_spacer]); d=pdist(mat,metric='hamming')*19.0; Z=linkage(d,method='average')
    plt.figure(figsize=(12,max(6,min(40,len(km)*0.12)))); show=len(km)<=250
    dendrogram(Z,labels=km.spacer_id.tolist() if show else None,orientation='right',leaf_font_size=6,no_labels=not show)
    plt.xlabel('Hamming distance (number of differing bases; 19-nt spacers)'); plt.title(f'K562 MHC exact-matching trimmed spacers (n={len(km):,})'); plt.tight_layout(); plt.savefig(plots_dir/'k562_trim1bp_spacer_hamming_dendrogram.png',dpi=200); plt.close()
if not k5.empty:
    pos=k5.sort_values(['contig','start0']).copy(); contigs=list(pos.contig.drop_duplicates()); offsets={}; running=0
    for c in contigs: offsets[c]=running; running += int(pos.loc[pos.contig==c,'end0'].max())+1000
    pos['plot_position']=pos.apply(lambda r: offsets[r.contig]+r.start0,axis=1)
    plt.figure(figsize=(13,4)); plt.scatter(pos.plot_position,np.zeros(len(pos)),s=8,alpha=.55)
    for c in contigs: plt.axvline(offsets[c],linewidth=.5)
    plt.yticks([]); plt.xlabel('Position along concatenated K562 hap2 MHC contig(s)'); plt.title(f'Exact 19-nt trimmed spacer placements in K562 hap2 MHC ({len(pos):,} alignment records)'); plt.tight_layout(); plt.savefig(plots_dir/'k562_mhc_trim1bp_spacer_positions.png',dpi=200); plt.close()
print(counts.to_string(index=False)); print(f'Unique trimmed 19-nt spacers matching K562 MHC exactly: {len(km):,}')
