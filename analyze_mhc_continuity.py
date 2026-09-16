"""Per-haplotype MHC continuity check: pick the candidate contig (largest
merged reference span at MAPQ >= --min-mapq) for hap1 and hap2, report flank
proximity and reference-span coverage against the extracted GRCh38 interval,
and write a gap table (both reference-relative and chromosome coordinates).

Reads data/{hap}.mhc.paf (from run_mhc_alignment.sh / the HPC transfer).
"""

import argparse
from pathlib import Path

import pandas as pd

from mhc_paf_utils import (
    candidate_contig, gap_table, load_paf, parse_reference_interval,
)

REPO = Path(__file__).resolve().parent
DATA = REPO / "data"


def data_file(name):
    """Prefer data/<name>; fall back to the repo root for the older layout."""
    candidate = DATA / name
    return candidate if candidate.exists() else REPO / name


ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--hap1-paf", type=Path, default=data_file("hap1.mhc.paf"))
ap.add_argument("--hap2-paf", type=Path, default=data_file("hap2.mhc.paf"))
ap.add_argument("--min-mapq", type=int, default=30,
                 help="minimum MAPQ to treat an alignment as confident, "
                      "excluding multi-mapped paralogous noise (default: %(default)s)")
ap.add_argument("--outdir", type=Path, default=REPO / "results" / "wtc11")
args = ap.parse_args()

args.outdir.mkdir(parents=True, exist_ok=True)

HAP_PAFS = {"hap1": args.hap1_paf, "hap2": args.hap2_paf}

summary_rows = []

for hap, paf_path in HAP_PAFS.items():
    print(f"\n=== {hap} ===")
    all_records, confident, excluded = load_paf(paf_path, min_mapq=args.min_mapq)

    chrom, interval_start, interval_end = parse_reference_interval(
        all_records["reference"].iloc[0]
    )
    # samtools faidx region syntax is 1-based inclusive, so the extracted
    # sequence is (interval_end - interval_start + 1) bp; trust the PAF's own
    # ref_len (column 7) as authoritative rather than re-deriving it, since
    # it's what alignment coordinates are actually measured against.
    ref_length = int(all_records["ref_len"].iloc[0])
    expected_length = interval_end - interval_start + 1
    if ref_length != expected_length:
        print(f"  WARNING: PAF ref_len {ref_length} != interval length "
              f"{expected_length} derived from '{chrom}:{interval_start}-{interval_end}'")

    contig, span, merged = candidate_contig(confident)
    span_pct = 100 * span / ref_length

    expected_file = data_file(f"{hap}.mhc_contigs.txt")
    expected = expected_file.read_text().split() if expected_file.exists() else []
    if expected and contig not in expected:
        print(f"  WARNING: candidate contig {contig} not in {expected_file.name} "
              f"({expected})")

    flank_start_dist = merged[0]["start"] - 0
    flank_end_dist = ref_length - merged[-1]["end"]

    gaps = gap_table(merged, interval_start)
    gaps.insert(0, "haplotype", hap)
    gaps.insert(1, "contig", contig)

    gap_path = args.outdir / f"mhc_continuity_{hap}.tsv"
    gaps.to_csv(gap_path, sep="\t", index=False)

    excluded_summary = (
        excluded.groupby("contig")
        .agg(n_records=("mapq", "size"), total_block_len=("block_len", "sum"))
        .sort_values("total_block_len", ascending=False)
        .reset_index()
    )
    excluded_path = args.outdir / f"mhc_lowmapq_excluded_summary_{hap}.tsv"
    excluded_summary.to_csv(excluded_path, sep="\t", index=False)

    print(f"  candidate contig: {contig}")
    print(f"  reference span coverage: {span} / {ref_length} bp ({span_pct:.1f}%)")
    print(f"  flank distances: start={flank_start_dist} bp, end={flank_end_dist} bp")
    print(f"  gaps: {len(gaps)}, total gap bp: {gaps['ref_gap_bp'].sum() if len(gaps) else 0}")
    print(f"  -> {gap_path}")
    print(f"  -> {excluded_path}")

    summary_rows.append({
        "haplotype": hap,
        "candidate_contig": contig,
        "reference_length_bp": ref_length,
        "merged_span_bp": span,
        "span_coverage_pct": span_pct,
        "flank_start_dist_bp": flank_start_dist,
        "flank_end_dist_bp": flank_end_dist,
        "n_gaps": len(gaps),
        "total_gap_bp": int(gaps["ref_gap_bp"].sum()) if len(gaps) else 0,
    })

summary = pd.DataFrame(summary_rows)
summary_path = args.outdir / "mhc_continuity_summary.tsv"
summary.to_csv(summary_path, sep="\t", index=False)
print(f"\n=== summary -> {summary_path} ===")
print(summary.to_string(index=False))
