#!/usr/bin/env python3
"""Audit all library guides in DHSs containing a missing, low-coverage WTC11 guide."""
from pathlib import Path
import numpy as np
import pandas as pd
from run_wtc11_analysis import read_hits, digest
import json

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results/wtc11'
manifest = json.loads((OUT / 'alignment_manifest.json').read_text())
for item in manifest['inputs'].values():
    assert digest(item['path']) == item['sha256'], f"Input changed: {item['path']}"
lookup = pd.read_csv(OUT / 'guide_mapping_lookup.tsv.gz', sep='\t')
assert lookup.guide_id.is_unique
# Independently recalculate the inherited +/-25 kb guide-match fractions.
gr = read_hits(manifest['inputs']['grch38']['path'])
gr = gr[gr.contig.eq('chr6') & gr.start0.ge(28_000_000) & gr.start0.lt(34_000_000)]
anchors = gr[~gr.spacer_id.duplicated(keep=False)].sort_values('start0')
positions = gr.sort_values(['spacer_id', 'start0']).drop_duplicates('spacer_id').set_index('spacer_id').start0
for hap in ['hap1', 'hap2']:
    hits = read_hits(OUT / f'wtc11_{hap}_exact_alignments_trim1bp.tsv.gz')
    present = set(hits.spacer_id)
    fractions = []
    for sid in lookup.spacer_id:
        pos = positions.get(sid, np.nan)
        local = anchors.loc[anchors.start0.between(pos-25_000, pos+25_000)]
        fractions.append(local.spacer_id.isin(present).mean())
    np.testing.assert_allclose(fractions, lookup[f'{hap}_local_coverage'], equal_nan=True)
seed = lookup.mapping_class.eq('missing_in_WTC11') & lookup.absence_in_assembly_gap
assert (lookup.loc[seed, ['hap1_local_coverage', 'hap2_local_coverage']] < .5).all().all()
dhs = set(lookup.loc[seed, 'intended_target_name'])
d = lookup.loc[lookup.intended_target_name.isin(dhs)].copy()
d['seed_missing_and_gap_flag'] = seed.loc[d.index]
d['other_guide'] = ~d.seed_missing_and_gap_flag
# Unique mapping here is a sequence-mapping criterion, not functional validation.
d['unique_in_both_no_recorded_offtarget'] = d.mapping_class.eq('unique_both') & d.n_hg38_offtarget_loci.eq(0)
top = pd.read_csv(OUT / 'screen_top_guides_mapping.tsv.gz', sep='\t', low_memory=False)
seed_ids = set(d.loc[d.seed_missing_and_gap_flag, 'guide_id'])
hit_seed = top[top.fractel_hit & top.guide_id.isin(seed_ids)]
rows = []
for name, group in d.groupby('intended_target_name'):
    others = group[group.other_guide]
    rows.append(dict(dhs=name, total_guides=len(group), missing_and_gap_flag=int(group.seed_missing_and_gap_flag.sum()),
                     other_guides=len(others), other_unique_both=int(others.mapping_class.eq('unique_both').sum()),
                     other_missing=int(others.mapping_class.eq('missing_in_WTC11').sum()),
                     other_multi=int(others.mapping_class.isin(['multi_hg38','multi_WTC11']).sum()),
                     other_single_hap=int(others.mapping_class.eq('hemizygous_WTC11').sum()),
                     other_single_hap_low_coverage=int(others.mapping_class.eq('hemizygous_assembly_gap').sum()),
                     other_is_problematic=int(others.is_problematic.sum()),
                     all_other_unique_both=bool(others.unique_in_both_no_recorded_offtarget.all()) if len(others) else False,
                     seed_top_slots_at_hits=int(hit_seed.dhs.eq(name).sum())))
summary = pd.DataFrame(rows)
assert (summary.other_guides == summary[['other_unique_both','other_missing','other_multi','other_single_hap','other_single_hap_low_coverage']].sum(axis=1)).all()
summary.to_csv(OUT / 'gap_flagged_DHS_summary.tsv', sep='\t', index=False)
d.sort_values(['intended_target_name','guide_start']).to_csv(OUT / 'gap_flagged_DHS_all_guides.tsv', sep='\t', index=False)
print(summary.to_string(index=False))
print('\nAll other guides:', d.loc[d.other_guide, 'mapping_class'].value_counts().to_dict())
print('Unique_both guides with recorded hg38 off-target loci:', int((d.other_guide & d.mapping_class.eq('unique_both') & d.n_hg38_offtarget_loci.gt(0)).sum()))
print('Validated input hashes and independently reproduced both WTC11 local-coverage columns.')
