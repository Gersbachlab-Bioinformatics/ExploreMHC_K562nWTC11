"""HLA/MHC gene-content footprint check: for every gene annotated in the
extracted GRCh38 interval, report what fraction of its genomic span (not
just exons -- this is a presence/continuity check, not allele typing) is
covered by each haplotype's candidate contig, the mean alignment identity
across that span, and whether a continuity gap (from analyze_mhc_continuity.py)
falls inside it.

This is genomic-footprint/identity evidence of preservation, not
confirmation of an intact coding sequence or a specific HLA allele call --
that distinction, and comparison against IPD-IMGT/HLA, is a later,
out-of-scope step.
"""

import argparse
from pathlib import Path

import pandas as pd

from mhc_paf_utils import (
    CLASSIC_HLA_GENES, candidate_contig, load_genes, load_paf,
    parse_reference_interval,
)

REPO = Path(__file__).resolve().parent
DATA = REPO / "data"


def data_file(name):
    candidate = DATA / name
    return candidate if candidate.exists() else REPO / name


ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--hap1-paf", type=Path, default=data_file("hap1.mhc.paf"))
ap.add_argument("--hap2-paf", type=Path, default=data_file("hap2.mhc.paf"))
ap.add_argument("--gtf", type=Path,
                 default=data_file("refseq/GCF_000001405.40-RS_2025_08_genomic.gtf.gz"))
ap.add_argument("--min-mapq", type=int, default=30)
ap.add_argument("--min-coverage-flag", type=float, default=95.0,
                 help="genes below this footprint coverage %% (either haplotype) "
                      "are highlighted in the filtered output (default: %(default)s)")
ap.add_argument("--outdir", type=Path, default=REPO / "results" / "wtc11")
args = ap.parse_args()

args.outdir.mkdir(parents=True, exist_ok=True)

HAP_PAFS = {"hap1": args.hap1_paf, "hap2": args.hap2_paf}


def interval_overlap(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def footprint(gene_ref_start, gene_ref_end, merged, on_contig):
    """gene_ref_start/end: gene span in ref-relative (0-based) coordinates.
    merged: candidate contig's merged alignment intervals.
    on_contig: candidate contig's confident alignment records (for identity).
    Returns (coverage_bp, coverage_pct, mean_identity_pct)."""
    gene_len = gene_ref_end - gene_ref_start
    covered = sum(
        interval_overlap(gene_ref_start, gene_ref_end, m["start"], m["end"])
        for m in merged
    )
    overlapping = on_contig[
        (on_contig["ref_end"] > gene_ref_start) & (on_contig["ref_start"] < gene_ref_end)
    ]
    mean_identity = overlapping["identity_pct"].mean() if len(overlapping) else None
    return covered, 100 * covered / gene_len if gene_len else 0.0, mean_identity


per_hap = {}
gap_tables = {}

for hap, paf_path in HAP_PAFS.items():
    print(f"\n=== {hap} ===")
    all_records, confident, _ = load_paf(paf_path, min_mapq=args.min_mapq)
    chrom, interval_start, interval_end = parse_reference_interval(
        all_records["reference"].iloc[0]
    )
    contig, span, merged = candidate_contig(confident)
    on_contig = confident[confident["contig"] == contig]
    per_hap[hap] = {
        "interval_start": interval_start,
        "interval_end": interval_end,
        "merged": merged,
        "on_contig": on_contig,
    }

    gap_path = args.outdir / f"mhc_continuity_{hap}.tsv"
    gap_tables[hap] = pd.read_csv(gap_path, sep="\t") if gap_path.exists() else pd.DataFrame(
        columns=["chrom_start_1based", "chrom_end_1based"]
    )

interval_start_1based = per_hap["hap1"]["interval_start"]
interval_end_1based = per_hap["hap1"]["interval_end"]
ref_length = interval_end_1based - interval_start_1based + 1
genes = load_genes(args.gtf, interval_start_1based, interval_end_1based)
print(f"\n{len(genes)} genes annotated in the interval")

rows = []
for _, g in genes.iterrows():
    gene_ref_start = g["start"] - interval_start_1based  # 1-based chrom -> 0-based ref-relative
    gene_ref_end = g["end"] - interval_start_1based + 1
    # Genes straddling the extraction boundary (e.g. GPX6, which starts
    # ~6.7 kb before this interval) can never reach 100% footprint coverage
    # regardless of assembly quality -- the missing part isn't in
    # mhc_ref.fa at all. Clip to the interval and score coverage against
    # the visible portion, flagging when that clip happened.
    truncated = gene_ref_start < 0 or gene_ref_end > ref_length
    gene_ref_start_clipped = max(gene_ref_start, 0)
    gene_ref_end_clipped = min(gene_ref_end, ref_length)
    row = {
        "gene": g["gene"],
        "gene_biotype": g["gene_biotype"],
        "chrom_start_1based": g["start"],
        "chrom_end_1based": g["end"],
        "strand": g["strand"],
        "is_classic_hla": g["gene"] in CLASSIC_HLA_GENES,
        "truncated_by_interval": truncated,
    }
    for hap in HAP_PAFS:
        covered, pct, identity = footprint(
            gene_ref_start_clipped, gene_ref_end_clipped,
            per_hap[hap]["merged"], per_hap[hap]["on_contig"],
        )
        gaps = gap_tables[hap]
        gap_overlap = bool(len(gaps)) and (
            (gaps["chrom_start_1based"] < g["end"]) & (gaps["chrom_end_1based"] > g["start"])
        ).any()
        row[f"{hap}_footprint_coverage_pct"] = pct
        row[f"{hap}_mean_identity_pct"] = identity
        row[f"{hap}_gap_overlap"] = gap_overlap
    rows.append(row)

table = pd.DataFrame(rows).sort_values("chrom_start_1based")
out_path = args.outdir / "mhc_gene_content_footprint.tsv"
table.to_csv(out_path, sep="\t", index=False)
print(f"-> {out_path}")

flagged = table[
    table["is_classic_hla"]
    | (table["hap1_footprint_coverage_pct"] < args.min_coverage_flag)
    | (table["hap2_footprint_coverage_pct"] < args.min_coverage_flag)
    | table["hap1_gap_overlap"]
    | table["hap2_gap_overlap"]
]
flagged_path = args.outdir / "mhc_gene_content_footprint_flagged.tsv"
flagged.to_csv(flagged_path, sep="\t", index=False)
print(f"-> {flagged_path} ({len(flagged)} of {len(table)} genes: classic HLA, "
      f"<{args.min_coverage_flag}% coverage, or gap overlap)")

print("\nClassic HLA genes:")
print(
    table[table["is_classic_hla"]]
    [["gene", "chrom_start_1based", "chrom_end_1based",
      "hap1_footprint_coverage_pct", "hap1_mean_identity_pct", "hap1_gap_overlap",
      "hap2_footprint_coverage_pct", "hap2_mean_identity_pct", "hap2_gap_overlap"]]
    .to_string(index=False, float_format=lambda v: f"{v:.1f}")
)
