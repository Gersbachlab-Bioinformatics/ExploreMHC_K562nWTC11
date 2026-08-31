#!/usr/bin/env python3
import re, sys
import pandas as pd
src, fasta_out, map_out = sys.argv[1:4]
df = pd.read_csv(src, sep='\t', compression='infer', dtype=str)
if 'spacer' not in df.columns:
    raise SystemExit("ERROR: input table has no 'spacer' column")
df['spacer'] = df['spacer'].fillna('').str.upper()
valid = df['spacer'].str.fullmatch(r'[ACGT]{20}', na=False)
dv = df.loc[valid].copy()
dv['trimmed_spacer'] = dv['spacer'].str[1:]
uniq = sorted(dv['trimmed_spacer'].unique())
seq_to_id = {s:f'trimmed_spacer_{i:05d}' for i,s in enumerate(uniq,1)}
with open(fasta_out,'w') as h:
    for s in uniq:
        h.write(f'>{seq_to_id[s]}\n{s}\n')
mapdf = pd.DataFrame({
    'row_index': dv.index,
    'guide_id': dv.get('guide_id', pd.Series(index=dv.index, dtype=str)),
    'spacer_id': dv['trimmed_spacer'].map(seq_to_id),
    'original_spacer': dv['spacer'],
    'trimmed_spacer': dv['trimmed_spacer']
})
mapdf.to_csv(map_out, sep='\t', index=False, compression='gzip')
print(f'Rows in source table: {len(df):,}')
print(f'Valid original 20-nt spacer rows: {len(dv):,}')
print(f'Unique original 20-nt spacers: {dv["spacer"].nunique():,}')
print(f'Unique trimmed 19-nt spacers: {len(uniq):,}')
