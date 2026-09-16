"""Dotplot comparing the extracted GRCh38 MHC interval (x-axis, chromosome
coordinates) against the candidate contig from each WTC11 haplotype
(y-axis, position along that contig) -- two side-by-side panels, each with
strand-colored alignment segments, gap shading from
analyze_mhc_continuity.py's output, and a classic-HLA-gene tick track.

Run analyze_mhc_continuity.py first (writes the gap tables this reads).
"""

import argparse
import collections
from pathlib import Path

import matplotlib.pyplot as plt

from mhc_paf_utils import (
    CLASSIC_HLA_GENES, candidate_contig, load_genes, load_paf,
    parse_reference_interval, to_chrom_pos,
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
ap.add_argument("--gapdir", type=Path, default=REPO / "results" / "wtc11",
                 help="directory containing mhc_continuity_{hap}.tsv (default: %(default)s)")
ap.add_argument("--outdir", type=Path, default=REPO / "results" / "wtc11")
args = ap.parse_args()

args.outdir.mkdir(parents=True, exist_ok=True)

HAP_PAFS = {"hap1": args.hap1_paf, "hap2": args.hap2_paf}
STRAND_COLOR = {"+": "#2a78d6", "-": "#eb6834"}
INK = "#0b0b0b"
INK_MUTED = "#52514e"
GAP_COLOR = "#e6e5e1"

fig, axes = plt.subplots(1, len(HAP_PAFS), figsize=(20, 11))

interval_start_1based = None
interval_end_1based = None

for ax, (hap, paf_path) in zip(axes, HAP_PAFS.items()):
    all_records, confident, _ = load_paf(paf_path, min_mapq=args.min_mapq)
    chrom, interval_start_1based, interval_end_1based = parse_reference_interval(
        all_records["reference"].iloc[0]
    )

    contig, span, merged = candidate_contig(confident)

    # Plot one segment per merged reference interval (i.e. the same
    # backbone used for the gap table), not one per raw alignment record.
    # The MHC region contains dispersed repeats that independently align
    # back to the same narrow reference window at MAPQ>=30 from dozens of
    # unrelated query positions across the contig (confirmed by inspecting
    # data/hap1.mhc.paf around chr6:30.08-30.09 Mb) -- real signal, but it
    # doesn't represent the assembly's colinear path through the interval,
    # and plotting every record drowns that path in scatter.
    for m in merged:
        open_r, close_r = m["open_rec"], m["close_rec"]
        x0 = to_chrom_pos(m["start"], interval_start_1based) / 1e6
        x1 = to_chrom_pos(m["end"], interval_start_1based) / 1e6
        y0 = (open_r["query_start"] if open_r["strand"] == "+" else open_r["query_end"]) / 1e6
        y1 = (close_r["query_end"] if close_r["strand"] == "+" else close_r["query_start"]) / 1e6
        ax.plot([x0, x1], [y0, y1], color=STRAND_COLOR[close_r["strand"]],
                 linewidth=1.5, alpha=0.85, marker="o", markersize=2)

    gap_path = args.gapdir / f"mhc_continuity_{hap}.tsv"
    if gap_path.exists():
        import pandas as pd
        gaps = pd.read_csv(gap_path, sep="\t")
        for _, g in gaps.iterrows():
            ax.axvspan(g["chrom_start_1based"] / 1e6, g["chrom_end_1based"] / 1e6,
                       color=GAP_COLOR, zorder=0)

    # hifiasm doesn't orient contigs to match the reference strand -- a
    # contig assembled on the reference's minus strand plots as a mirrored
    # (decreasing) diagonal, which isn't a biological inversion, just an
    # arbitrary assembly-output orientation. Flip the axis so both panels
    # read left-to-right the same way, and label it so the flip is explicit
    # rather than silently changing what "up" means between panels.
    backbone_strands = collections.Counter(m["close_rec"]["strand"] for m in merged)
    flipped = backbone_strands["-"] > backbone_strands["+"]
    if flipped:
        ax.invert_yaxis()

    ax.set_ylabel(f"{hap}\n{contig}\nposition (Mb)", color=INK)
    ax.set_title(
        f"{hap}: {contig}  |  reference span coverage "
        f"{100 * span / int(all_records['ref_len'].iloc[0]):.1f}%"
        + ("  |  axis flipped: contig is minus-strand vs. reference" if flipped else ""),
        loc="left", color=INK_MUTED, fontsize=10,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="both", color="#f0efec", linewidth=0.6)
    ax.set_axisbelow(True)
    # Equal aspect: 1 Mb on the reference axis == 1 Mb on the contig axis,
    # so a true 1:1 colinear alignment reads as a 45-degree diagonal.
    ax.set_aspect("equal", adjustable="box")

    genes = load_genes(args.gtf, interval_start_1based, interval_end_1based)
    hla = genes[genes["gene"].isin(CLASSIC_HLA_GENES)].sort_values("start")

    for _, g in hla.iterrows():
        mid_mb = (g["start"] + g["end"]) / 2 / 1e6
        ax.axvline(mid_mb, color="#9a9890", linewidth=0.5, linestyle=":", zorder=0)

    for i, (_, g) in enumerate(hla.iterrows()):
        mid_mb = (g["start"] + g["end"]) / 2 / 1e6
        stagger = 1 + 0.10 * (i % 3)
        ax.annotate(
            g["gene"], xy=(mid_mb, stagger), xycoords=("data", "axes fraction"),
            rotation=90, fontsize=6, ha="center", va="bottom", color=INK_MUTED,
        )

    ax.set_xlabel(f"{chrom} position (Mb)", color=INK)
fig.suptitle("WTC11 MHC assembly vs. GRCh38: reference vs. hap1 / hap2 candidate contigs",
             color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.90))

out_path = args.outdir / "mhc_ref_vs_hap1_hap2_dotplot.png"
fig.savefig(out_path, dpi=300)
plt.close(fig)
print(f"-> {out_path}")
print(f"HLA genes plotted: {len(hla)} of {len(CLASSIC_HLA_GENES)} classic HLA genes found in interval")
