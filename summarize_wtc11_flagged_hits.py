#!/usr/bin/env python3
"""Resolve flagged WTC11 hit slots to unique guides, DHSs, genes and library annotations."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results/wtc11'
meta = pd.read_csv(ROOT / 'data/IGVFFI2404DYFG.tsv.gz', sep='\t')
assert meta.guide_id.is_unique
slots = pd.read_csv(OUT / 'screen_top_guides_mapping.tsv.gz', sep='\t', low_memory=False)
hits = slots[slots.fractel_hit].copy()
hits = hits.merge(meta[['guide_id', 'intended_target_name', 'genomic_element', 'type', 'description']],
                  on='guide_id', how='left', validate='many_to_one', indicator=True)
assert hits._merge.eq('both').all()
assert hits.dhs.eq(hits.intended_target_name).all()
# Require a consistent annotation across ALL library guides in each hit DHS.
annotations = meta[meta.intended_target_name.isin(hits.dhs)].groupby('intended_target_name').genomic_element.agg(lambda x: sorted(set(x.fillna('unknown'))))
assert annotations.map(len).eq(1).all(), 'Mixed DHS annotations require review'
hits['genomic_element'] = hits.genomic_element.fillna('unknown')
f = hits[hits.relevant_problem].copy()
old = pd.read_csv(OUT / 'fractel_hits_with_problem_guides.tsv', sep='\t')
keys=['screen','dhs','ensembl_id','guide_id']
assert set(map(tuple,f[keys].values)) == set(map(tuple,old[keys].values))

def joined(s):
    return ';'.join(sorted(set(s.dropna().astype(str))))

def metrics(d):
    return dict(slots=len(d), unique_guides=d.guide_id.nunique(), unique_DHSs=d.dhs.nunique(),
                unique_genes=d.ensembl_id.nunique(), unique_DHS_gene_pairs=len(d[['dhs','ensembl_id']].drop_duplicates()),
                screen_specific_pairs=len(d[['screen','dhs','ensembl_id']].drop_duplicates()))

summary=pd.DataFrame([dict(scope='all_hit_slots',**metrics(hits)),dict(scope='flagged_hit_slots',**metrics(f))])
classes=pd.DataFrame([dict(mapping_class=c,**metrics(d)) for c,d in f.groupby('mapping_class')])
elements=[]
for c,d in hits.groupby('genomic_element'):
    flagged=d[d.relevant_problem]
    elements.append(dict(genomic_element=c, total_hit_DHSs=d.dhs.nunique(), **metrics(flagged)))
elements=pd.DataFrame(elements)
pairs=hits.groupby(['screen','dhs','ensembl_id','gene_symbol','genomic_element'],as_index=False).agg(
    top_slots=('guide_id','size'), flagged_slots=('relevant_problem','sum'))
pairs=pairs[pairs.flagged_slots.gt(0)].copy()
pairs['all_top_guides_flagged']=pairs.top_slots.eq(pairs.flagged_slots)
dhs=f.groupby(['dhs','genomic_element'],as_index=False).agg(
    flagged_slots=('guide_id','size'), unique_flagged_guides=('guide_id','nunique'),
    genes=('gene_symbol',joined), gene_ids=('ensembl_id',joined), screens=('screen',joined),
    mapping_classes=('mapping_class',joined), low_coverage_flagged_slots=('absence_in_assembly_gap','sum'))
for flag, mask in [('missing',f.mapping_class.eq('missing_in_WTC11')),('multi',f.mapping_class.isin(['multi_hg38','multi_WTC11']))]:
    dhs[f'unique_{flag}_guides']=dhs.dhs.map(f[mask].groupby('dhs').guide_id.nunique()).fillna(0).astype(int)
genes=f.groupby(['ensembl_id','gene_symbol'],as_index=False).agg(unique_DHSs=('dhs','nunique'),
    unique_flagged_guides=('guide_id','nunique'), flagged_slots=('guide_id','size'),
    DHSs=('dhs',joined), screens=('screen',joined), mapping_classes=('mapping_class',joined))
guides=f.groupby(['guide_id','dhs','genomic_element','mapping_class'],as_index=False).agg(
    flagged_slots=('screen','size'), genes=('gene_symbol',joined),screens=('screen',joined),
    best_rank=('guide_rank','min'), absence_in_assembly_gap=('absence_in_assembly_gap','first'))
for name,d in [('flagged_hit_overview',summary),('flagged_hit_mapping_classes',classes),
               ('flagged_hit_element_classes',elements),('flagged_hit_DHSs',dhs),
               ('flagged_hit_genes',genes),('flagged_hit_unique_guides',guides),
               ('flagged_hit_pairs',pairs),('flagged_hit_slots_with_metadata',f.drop(columns='_merge'))]:
    d.to_csv(OUT/f'{name}.tsv',sep='\t',index=False)
print(summary.to_string(index=False)); print(classes.to_string(index=False)); print(elements.to_string(index=False))
print('\nDHSs\n',dhs.to_string(index=False))
print('\nAll top slots flagged\n',pairs[pairs.all_top_guides_flagged].to_string(index=False))
print('\nGenes\n',genes[['gene_symbol','unique_DHSs','unique_flagged_guides','flagged_slots']].to_string(index=False))
print('\nMissing without low coverage:',metrics(f[f.mapping_class.eq('missing_in_WTC11') & ~f.absence_in_assembly_gap]))
print('Missing with low coverage:',metrics(f[f.mapping_class.eq('missing_in_WTC11') & f.absence_in_assembly_gap]))
