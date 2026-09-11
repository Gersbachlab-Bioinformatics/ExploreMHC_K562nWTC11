#!/usr/bin/env python3
"""Replicate the exact trim1bp MHC analysis for both WTC11 haplotypes.

Run in the mambaforge environment. Reuses the existing GRCh38 alignment table;
all WTC11 outputs stay under results/wtc11. No downloads are performed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from scipy.stats import binomtest, false_discovery_control

REPO = Path(__file__).resolve().parent
HAPS = ('hap1', 'hap2')


def run(*command, log=None):
    print(' '.join(map(str, command)), flush=True)
    if log:
        with open(log, 'w') as handle:
            subprocess.run(list(map(str, command)), check=True, stdout=handle, stderr=subprocess.STDOUT)
    else:
        subprocess.run(list(map(str, command)), check=True)


def read_hits(path):
    d = pd.read_csv(path, sep='\t')
    d = d.loc[d.NM.eq(0) & ~d.supplementary.astype(bool)].copy()
    return d.drop_duplicates(['spacer_id', 'contig', 'start0', 'end0', 'strand'])


def save(d, path):
    d.to_csv(path, sep='\t', index=False)


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def fasta_stats(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                rows.append(dict(contig=line[1:].split()[0], length=0, n_bases=0))
            else:
                seq = line.strip().upper()
                if seq:
                    if not rows:
                        raise ValueError(f'FASTA sequence before header: {path}')
                    rows[-1]['length'] += len(seq)
                    rows[-1]['n_bases'] += seq.count('N')
    if not rows or any(r['length'] == 0 for r in rows):
        raise ValueError(f'Empty FASTA/contig: {path}')
    if len({r['contig'] for r in rows}) != len(rows):
        raise ValueError(f'Duplicate FASTA contig IDs: {path}')
    return rows


def category(g, h):
    return np.select([g & h, g & ~h, ~g & h],
                     ['Both', 'GRCh38 only', 'WTC11 MHC only'], default='Neither')


def loss_analysis(frame, pos_col, end_col, out, label):
    """Same 10 kb windows, <=10 kb clusters and 100 kb density baseline."""
    out.mkdir(parents=True, exist_ok=True)
    x = frame.dropna(subset=[pos_col, end_col]).copy()
    x['position'] = x[pos_col].astype(int)
    x['position_end'] = x[end_col].astype(int)
    # Intended-coordinate cohort includes guides without a GRCh38 exact hit;
    # preserve the original definition: lost requires a GRCh38 exact match.
    x['lost'] = x.GRCh38_exact_match_count.gt(0) & x[f'{label}_exact_match_count'].eq(0)
    x = x.sort_values(['position', 'guide_id'])
    save(x, out / 'guide_summary.tsv.gz')
    lost = x.loc[x.lost].copy()
    lost['cluster_id'] = (lost.position.diff().isna() | lost.position.diff().gt(10_000)).cumsum()
    save(lost, out / 'missing_guides.tsv.gz')
    clusters = lost.groupby('cluster_id', as_index=False).agg(
        start0=('position', 'min'), end0=('position_end', 'max'), n_missing_guides=('guide_id', 'size'))
    clusters['span_bp'] = clusters.end0 - clusters.start0
    save(clusters, out / 'missing_guide_clusters.tsv')
    if x.empty:
        return dict(cohort=out.parent.name, reference=label, total=0, missing=0, loss_percent=None)
    for width in [10_000, 100_000]:
        x['bin_start'] = x.position // width * width
        s = x.groupby('bin_start').agg(total_guides=('lost', 'size'), missing_guides=('lost', 'sum'))
        s = s.reindex(np.arange(s.index.min(), s.index.max() + width, width), fill_value=0).reset_index()
        s['bin_end'] = s.bin_start + width
        s['found_guides'] = s.total_guides - s.missing_guides
        s['loss_fraction'] = s.missing_guides / s.total_guides.replace(0, np.nan)
        s['expected_missing'] = s.total_guides * x.lost.mean()
        s['enrichment'] = s.missing_guides / s.expected_missing.replace(0, np.nan)
        s['binom_p'] = [binomtest(int(m), int(n), x.lost.mean(), alternative='greater').pvalue
                        if n else np.nan for m, n in zip(s.missing_guides, s.total_guides)]
        tested = s.binom_p.notna()
        # Full BH adjustment includes the reverse cumulative minimum.
        s.loc[tested, 'binom_q'] = false_discovery_control(s.loc[tested, 'binom_p'].values)
        save(s, out / f'loss_by_{width // 1000}kb_window.tsv')
        centers = (s.bin_start + width / 2) / 1e6
        if width == 10_000:
            for col, ylabel in [('missing_guides', 'Missing guide rows per 10 kb'), ('loss_fraction', 'Fraction of guide rows lost')]:
                fig, ax = plt.subplots(figsize=(14, 4))
                ax.bar(centers, s[col], width=width / 1e6)
                ax.set(xlabel='GRCh38 chr6 position (Mb)', ylabel=ylabel, title=f'{label}: {out.parent.name} coordinates')
                if col == 'loss_fraction':
                    ax.set_ylim(0, 1)
                fig.tight_layout(); fig.savefig(out / f'{col}.png', dpi=160); plt.close(fig)
        else:
            fig, ax = plt.subplots(figsize=(14, 4))
            ax.bar(centers - .021, s.missing_guides, width=.042, label='Observed missing')
            ax.bar(centers + .021, s.expected_missing, width=.042, label='Expected from guide density')
            sig = s.binom_q.lt(.05)
            ax.plot(centers[sig], s.loc[sig, ['missing_guides', 'expected_missing']].max(axis=1) + 1.5,
                    'k*', label='BH q < 0.05')
            ax.set(xlabel='GRCh38 chr6 position (Mb)', ylabel='Missing guide rows per 100 kb', title=label)
            ax.legend(); fig.tight_layout(); fig.savefig(out / 'missing_observed_vs_expected.png', dpi=160); plt.close(fig)
            fig, ax = plt.subplots(figsize=(14, 4))
            right = ax.twinx()
            b1 = ax.bar(centers - .021, s.found_guides, width=.042, color='#2a78d6', label='Found (left)')
            b2 = right.bar(centers + .021, s.missing_guides, width=.042, color='#eb6834', label='Missing (right)')
            ax.set(xlabel='GRCh38 chr6 position (Mb)', ylabel='Found guide rows', title=f'{label}: independent y-scales')
            right.set_ylabel('Missing guide rows'); ax.legend(handles=[b1, b2])
            fig.tight_layout(); fig.savefig(out / 'found_vs_missing_dual_axis.png', dpi=160); plt.close(fig)
    return dict(cohort=out.parent.name, reference=label, total=len(x), missing=int(x.lost.sum()),
                loss_percent=100 * x.lost.mean())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--guides', type=Path, default=REPO / 'data/IGVFFI2404DYFG.tsv.gz')
    ap.add_argument('--grch38', type=Path, default=REPO / 'data/grch38_exact_alignments_trim1bp.tsv.gz')
    for hap in HAPS:
        ap.add_argument(f'--{hap}', type=Path, default=REPO / f'data/wtc11.{hap}.mhc_contigs.fa')
    ap.add_argument('--outdir', type=Path, default=REPO / 'results/wtc11')
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--reuse-alignments', action='store_true', help='reuse existing parsed WTC11 alignments (inputs must match manifest)')
    ap.add_argument('--skip-dendrogram', action='store_true', help='omit quadratic-memory sequence clustering')
    args = ap.parse_args()
    out = args.outdir.resolve()
    for folder in ['', 'tmp', 'index', 'sam', 'logs', 'plots']:
        (out / folder).mkdir(parents=True, exist_ok=True)
    inputs = {key: dict(path=str(getattr(args, key).resolve()), sha256=digest(getattr(args, key)))
              for key in ['guides', 'grch38', *HAPS]}
    manifest = out / 'alignment_manifest.json'
    if args.reuse_alignments:
        if not manifest.exists() or json.loads(manifest.read_text())['inputs'] != inputs:
            raise SystemExit('Cannot reuse alignments: missing manifest or changed inputs. Run without --reuse-alignments.')
    run(sys.executable, REPO / 'extract_spacers_trim1bp.py', args.guides,
        out / 'tmp/spacers_19nt.fa', out / 'tmp/spacer_map.tsv.gz')
    mp = pd.read_csv(out / 'tmp/spacer_map.tsv.gz', sep='\t')
    guides = pd.read_csv(args.guides, sep='\t').reset_index(names='row_index')
    guides = guides.merge(mp[['row_index', 'spacer_id', 'trimmed_spacer']], on='row_index', how='left', validate='one_to_one')
    # The cached reference table uses sorted spacer IDs; verify against the
    # original per-guide summary rather than trusting matching ID ranges alone.
    old_path = REPO / 'data/guide_exact_match_summary_trim1bp.tsv.gz'
    if not old_path.exists():
        raise SystemExit(f'Missing ID provenance table: {old_path}')
    old = pd.read_csv(old_path, sep='\t')[['spacer_id', 'spacer']].drop_duplicates()
    old['trimmed_spacer'] = old.spacer.str.upper().str[1:]
    check = mp.merge(old[['spacer_id', 'trimmed_spacer']], on='spacer_id', how='left', suffixes=('', '_cached'), validate='many_to_one')
    if not check.trimmed_spacer.eq(check.trimmed_spacer_cached).all():
        raise SystemExit('Spacer sequences/IDs disagree with the cached GRCh38 provenance table.')
    hits = {'GRCh38': read_hits(args.grch38)}
    stats = []
    for hap in HAPS:
        stats.extend(dict(haplotype=hap, **r) for r in fasta_stats(getattr(args, hap)))
        dest = out / f'wtc11_{hap}_exact_alignments_trim1bp.tsv.gz'
        if not args.reuse_alignments:
            run('bowtie2-build', '--threads', args.threads, getattr(args, hap), out / f'index/{hap}', log=out / f'logs/{hap}_build.log')
            run('bowtie2', '--end-to-end', '-a', '--no-unal', '-f', '-p', args.threads,
                '-x', out / f'index/{hap}', '-U', out / 'tmp/spacers_19nt.fa', '-S', out / f'sam/{hap}.sam', log=out / f'logs/{hap}_align.log')
            run(sys.executable, REPO / 'parse_exact_alignments_trim1bp.py', out / f'sam/{hap}.sam', f'WTC11_{hap}_MHC', dest)
        hits[f'WTC11_{hap}_MHC'] = read_hits(dest)
    if not args.reuse_alignments:
        manifest.write_text(json.dumps(dict(inputs=inputs, command='bowtie2 --end-to-end -a --no-unal -f',
                                             bowtie2_version=subprocess.check_output(['bowtie2', '--version'], text=True).splitlines()[0]), indent=2) + '\n')
    save(pd.DataFrame(stats), out / 'assembly_contig_summary.tsv')
    for label, h in hits.items():
        if not set(h.spacer_id) <= set(mp.spacer_id):
            raise ValueError(f'{label}: unknown spacer IDs')
        guides[f'{label}_exact_match_count'] = guides.spacer_id.map(h.groupby('spacer_id').size()).fillna(0).astype(int)
        formatted = h.contig + ':' + h.start0.astype(str) + '-' + h.end0.astype(str) + ':' + h.strand
        guides[f'{label}_alignments'] = guides.spacer_id.map(formatted.groupby(h.spacer_id).agg(';'.join)).fillna('')
    guides['WTC11_MHC_exact_match_count'] = guides.WTC11_hap1_MHC_exact_match_count + guides.WTC11_hap2_MHC_exact_match_count
    g = guides.GRCh38_exact_match_count.gt(0)
    labels = ['WTC11_hap1_MHC', 'WTC11_hap2_MHC', 'WTC11_MHC']
    for label in labels:
        guides[f'{label}_category'] = category(g, guides[f'{label}_exact_match_count'].gt(0))
        guides.loc[guides.spacer_id.isna(), f'{label}_category'] = 'Invalid spacer'
    guides['exact_match_category'] = guides.WTC11_MHC_category
    guides['wtc11_haplotypes'] = np.select([
        guides.WTC11_hap1_MHC_exact_match_count.gt(0) & guides.WTC11_hap2_MHC_exact_match_count.gt(0),
        guides.WTC11_hap1_MHC_exact_match_count.gt(0), guides.WTC11_hap2_MHC_exact_match_count.gt(0)],
        ['both', 'hap1_only', 'hap2_only'], default='neither')
    summary = out / 'guide_exact_match_summary_trim1bp.tsv.gz'
    save(guides, summary)
    membership = guides[['spacer_id', 'trimmed_spacer'] + [f'{l}_exact_match_count' for l in ['GRCh38', *labels]]].drop_duplicates()
    save(membership, out / 'unique_trimmed_spacer_membership.tsv.gz')
    membership.filter(like='_exact_match_count').gt(0).corr().to_csv(out / 'reference_membership_correlation.tsv', sep='\t')
    counts = pd.concat([guides[f'{l}_category'].value_counts().rename_axis('category').reset_index(name='guide_rows').assign(reference=l) for l in labels])
    save(counts, out / 'exact_match_category_counts.tsv')
    gr = hits['GRCh38']
    mhc = gr.loc[gr.contig.eq('chr6') & gr.start0.ge(28_000_000) & gr.start0.lt(34_000_000)]
    # Choose coordinates per guide row, not per spacer, so duplicate spacer
    # sequences with different intended loci cannot relocate each other.
    coords = guides[['row_index', 'spacer_id', 'guide_start']].merge(mhc[['spacer_id', 'start0', 'end0']], on='spacer_id')
    coords['offset'] = (coords.start0 - pd.to_numeric(coords.guide_start, errors='coerce')).abs()
    coords = coords.sort_values(['row_index', 'offset', 'start0']).drop_duplicates('row_index')
    actual = guides.merge(coords[['row_index', 'start0', 'end0']], on='row_index', validate='one_to_one')
    intended = guides.loc[guides.intended_target_chr.eq('chr6') &
                         pd.to_numeric(guides.intended_target_start, errors='coerce').ge(28_000_000) &
                         pd.to_numeric(guides.intended_target_end, errors='coerce').le(34_000_000)].copy()
    for c in ['guide_start', 'guide_end']:
        intended[c] = pd.to_numeric(intended[c], errors='coerce')
    losses = []
    for label in labels:
        losses.append(loss_analysis(actual, 'start0', 'end0', out / 'realigned' / label, label))
        losses.append(loss_analysis(intended, 'guide_start', 'guide_end', out / 'intended' / label, label))
    save(pd.DataFrame(losses), out / 'loss_summary.tsv')
    for label in labels[:2]:
        h = hits[label]
        fig, axes = plt.subplots(max(1, h.contig.nunique()), 1, figsize=(14, max(3, h.contig.nunique() * 2)), squeeze=False)
        for ax, (contig, group) in zip(axes.flat, h.groupby('contig')):
            ax.scatter(group.start0 / 1e6, np.zeros(len(group)), s=5, alpha=.5)
            ax.set(xlabel=f'{contig} position (Mb)', yticks=[])
        fig.suptitle(f'{label}: all exact spacer placements'); fig.tight_layout()
        fig.savefig(out / f'plots/{label}_positions.png', dpi=160); plt.close(fig)
    matching = membership.loc[membership.WTC11_MHC_exact_match_count.gt(0)]
    with open(out / 'wtc11_matching_unique_trimmed_spacers.fa', 'w') as f:
        for r in matching.itertuples():
            f.write(f'>{r.spacer_id}\n{r.trimmed_spacer}\n')
    if not args.skip_dendrogram and len(matching) >= 2:
        print(f'Clustering {len(matching):,} matching 19-mers (quadratic memory).', flush=True)
        encoded = np.array([list(s.encode()) for s in matching.trimmed_spacer], dtype=np.uint8)
        tree = linkage(pdist(encoded, metric='hamming') * 19, method='average')
        np.save(out / 'wtc11_spacer_linkage.npy', tree)
        fig, ax = plt.subplots(figsize=(12, 10))
        dendrogram(tree, no_labels=True, orientation='right', ax=ax)
        ax.set(xlabel='Hamming distance (differing bases in 19 nt)', title=f'WTC11 either-haplotype matching spacers (n={len(matching):,})')
        fig.tight_layout(); fig.savefig(out / 'plots/wtc11_spacer_hamming_dendrogram.png', dpi=160); plt.close(fig)
    run(sys.executable, REPO / 'build_guide_mapping_lookup.py', '--cell-line', 'WTC11', '--summary', summary,
        '--grch38', args.grch38, '--hap1', out / 'wtc11_hap1_exact_alignments_trim1bp.tsv.gz',
        '--hap2', out / 'wtc11_hap2_exact_alignments_trim1bp.tsv.gz', '--outdir', out, log=out / 'logs/lookup.log')
    run(sys.executable, REPO / 'check_problem_guides_in_screens.py', '--cell-line', 'WTC11',
        '--lookup', out / 'guide_mapping_lookup.tsv.gz', '--screens', 'ipsc', 'npc', '--outdir', out, log=out / 'logs/screens.log')
    # Compare to both cached K562 haplotypes on the identical GRCh38 MHC cohort.
    k1 = read_hits(REPO / 'data/k562_hap1_exact_alignments_trim1bp.tsv.gz')
    k2 = read_hits(REPO / 'data/k562_mhc_exact_alignments_trim1bp.tsv.gz')
    comparison = actual[['guide_id', 'spacer_id', 'start0', 'wtc11_haplotypes']].copy()
    comparison['present_WTC11'] = actual.WTC11_MHC_exact_match_count.gt(0).values
    comparison['present_K562'] = comparison.spacer_id.isin(set(k1.spacer_id) | set(k2.spacer_id))
    save(comparison, out / 'k562_wtc11_guide_comparison.tsv.gz')
    save(comparison.groupby(['present_K562', 'present_WTC11']).size().reset_index(name='guide_rows'), out / 'k562_wtc11_presence_counts.tsv')
    print(pd.DataFrame(losses).to_string(index=False))
    print(f'Completed: {out}')


if __name__ == '__main__':
    main()
