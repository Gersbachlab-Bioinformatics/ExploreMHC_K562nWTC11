#!/usr/bin/env python3
"""Validate every reported WTC11 placement and downstream accounting."""
from pathlib import Path
import argparse
import json
import tempfile
import unittest

import pandas as pd
from run_wtc11_analysis import read_hits, digest

REPO = Path(__file__).resolve().parent


def validate(out):
    manifest = json.loads((out / 'alignment_manifest.json').read_text())
    for entry in manifest['inputs'].values():
        assert digest(entry['path']) == entry['sha256'], entry['path']
    guides = pd.read_csv(out / 'guide_exact_match_summary_trim1bp.tsv.gz', sep='\t')
    assert guides.row_index.is_unique
    assert len(guides) == len(pd.read_csv(manifest['inputs']['guides']['path'], sep='\t'))
    seqs = guides.set_index('spacer_id').trimmed_spacer.to_dict()
    rc = str.maketrans('ACGT', 'TGCA')
    for hap in ['hap1', 'hap2']:
        contigs, name, parts = {}, None, []
        with open(manifest['inputs'][hap]['path']) as f:
            for line in f:
                if line.startswith('>'):
                    if name is not None:
                        contigs[name] = ''.join(parts).upper()
                    name, parts = line[1:].split()[0], []
                else:
                    parts.append(line.strip())
        if name is not None:
            contigs[name] = ''.join(parts).upper()
        hits = read_hits(out / f'wtc11_{hap}_exact_alignments_trim1bp.tsv.gz')
        for r in hits.itertuples():
            expected = seqs[r.spacer_id]
            if r.strand == '-':
                expected = expected.translate(rc)[::-1]
            assert r.end0 - r.start0 == 19 and r.CIGAR == '19M', r
            assert contigs[r.contig][r.start0:r.end0] == expected, r
        count = guides.spacer_id.map(hits.groupby('spacer_id').size()).fillna(0)
        assert count.eq(guides[f'WTC11_{hap}_MHC_exact_match_count']).all()
        print(f'{hap}: verified all {len(hits):,} exact placements directly against FASTA')
    assert guides.WTC11_MHC_exact_match_count.eq(guides.WTC11_hap1_MHC_exact_match_count + guides.WTC11_hap2_MHC_exact_match_count).all()
    assert guides.exact_match_category.eq('GRCh38 only').equals(guides.GRCh38_exact_match_count.gt(0) & guides.WTC11_MHC_exact_match_count.eq(0))
    lookup = pd.read_csv(out / 'guide_mapping_lookup.tsv.gz', sep='\t')
    assert len(lookup) == len(guides) and lookup.guide_id.is_unique
    assert not any('k562' in c.lower() for c in lookup.columns)
    paired = lookup.merge(guides, on='guide_id', suffixes=('_lookup', '_summary'), validate='one_to_one')
    assert paired.n_wtc11_loci.eq(paired.WTC11_MHC_exact_match_count).all()
    for cohort in ['realigned', 'intended']:
        for label in ['WTC11_hap1_MHC', 'WTC11_hap2_MHC', 'WTC11_MHC']:
            root = out / cohort / label
            frame = pd.read_csv(root / 'guide_summary.tsv.gz', sep='\t')
            for width in [10, 100]:
                bins = pd.read_csv(root / f'loss_by_{width}kb_window.tsv', sep='\t')
                assert bins.total_guides.sum() == len(frame)
                assert bins.missing_guides.sum() == frame.lost.sum()
                tested = bins.dropna(subset=['binom_p']).sort_values('binom_p')
                assert (tested.binom_q.diff().dropna() >= -1e-12).all()
            clusters = pd.read_csv(root / 'missing_guide_clusters.tsv', sep='\t')
            assert clusters.n_missing_guides.sum() == frame.lost.sum()
    screens = pd.read_csv(out / 'screen_top_guides_mapping.tsv.gz', sep='\t', low_memory=False)
    assert screens.screen.str.startswith(('ipsc', 'npc')).all()
    assert screens.screen_is_wtc11.all()
    assert screens.guide_rank.between(1, 3).all()
    assert screens.mapping_class.notna().all()
    assert not screens.mapping_class.str.contains('K562').any()
    print('Passed: input hashes, guide IDs/counts, union logic, lookup joins, loss windows/clusters, BH monotonicity, WTC11 screen scope')


class HitFilteringTest(unittest.TestCase):
    def test_secondary_is_retained_supplementary_and_mismatches_are_removed(self):
        rows = [dict(spacer_id='s1', contig='c', start0=i, end0=i+19,
                     strand='+', NM=nm, secondary=secondary, supplementary=supp)
                for i, nm, secondary, supp in [(0, 0, False, False), (50, 0, True, False),
                                                (100, 0, False, True), (150, 1, False, False)]]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'hits.tsv'
            pd.DataFrame(rows + [rows[1]]).to_csv(path, sep='\t', index=False)
            self.assertEqual(read_hits(path).start0.tolist(), [0, 50])


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outdir', type=Path, default=REPO / 'results/wtc11')
    args = ap.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(HitFilteringTest)
    if not unittest.TextTestRunner().run(suite).wasSuccessful():
        raise SystemExit(1)
    validate(args.outdir)
