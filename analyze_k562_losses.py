
import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Files
#
# Each input defaults to the copy checked into the repo, so the script runs
# as-is from anywhere; every path can be overridden independently.

REPO = Path(__file__).resolve().parent
DATA = REPO / "data"


def data_file(name):
    """Prefer data/<name>; fall back to the repo root for the older layout."""
    candidate = DATA / name
    return candidate if candidate.exists() else REPO / name

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--summary", type=Path,
                default=data_file("guide_exact_match_summary_trim1bp.tsv.gz"),
                help="per-guide exact-match summary (default: %(default)s)")
ap.add_argument("--grch38", type=Path,
                default=data_file("grch38_exact_alignments_trim1bp.tsv.gz"),
                help="exact GRCh38 alignment records (default: %(default)s)")
ap.add_argument("--outdir", type=Path, default=REPO / "results",
                help="where to write tables and plots (default: %(default)s)")
ap.add_argument("--plot-bin-size", type=int, default=100_000,
                help="bin size for the found-vs-missing comparison figures; the "
                     "10-kb analysis bins are too dense to draw as grouped bars "
                     "(default: %(default)s)")
args = ap.parse_args()

SUMMARY, GRCH38, OUTDIR = args.summary, args.grch38, args.outdir
PLOT_BIN_SIZE = args.plot_bin_size
OUTDIR.mkdir(parents=True, exist_ok=True)

# Broad GRCh38 MHC interval
MHC_CHR = "chr6"
MHC_START = 28_000_000
MHC_END = 34_000_000

BIN_SIZE = 10_000       # 10 kb
CLUSTER_GAP = 10_000    # lost guides <=10 kb apart = same cluster

# Read files

summary = pd.read_csv(
    SUMMARY,
    sep="\t",
    compression="infer"
)

gr = pd.read_csv(
    GRCH38,
    sep="\t",
    compression="infer"
)

# Keep exact GRCh38 alignments.
#
# Do NOT filter on the `secondary` flag here. bowtie2 was run with -a, and an
# MHC guide typically has ~7 equally-good exact hits: chr6 plus the seven alt
# MHC haplotype contigs (chr6_GL000250v2_alt ... chr6_GL000256v2_alt). Exactly
# one of those is flagged primary, chosen arbitrarily among equal-scoring hits,
# and it lands on chr6 for only ~18% of guides. Dropping secondary alignments
# therefore discards ~82% of MHC-targeting guides at random rather than
# selecting the best ones. Restricting to contig == chr6 below is what picks
# the primary-assembly locus.

gr = gr[gr["NM"] == 0].copy()
if "supplementary" in gr.columns:
    gr = gr[~gr["supplementary"].astype(bool)]

# Restrict GRCh38 alignments to MHC

gr_mhc = gr[
    (gr["contig"] == MHC_CHR) &
    (gr["start0"] >= MHC_START) &
    (gr["start0"] < MHC_END)
].copy()

print("Exact GRCh38 MHC alignment records:", len(gr_mhc))
print("Unique spacers in GRCh38 MHC:",
      gr_mhc["spacer_id"].nunique())

# Get one MHC coordinate per spacer.
#
# 253 guides have more than one exact 19-nt locus in the MHC on the chr6
# primary assembly -- overwhelmingly the ~32.7 kb RCCX/C4 segmental duplication
# near 32.0 Mb, plus smaller paralogous pairs. Simply taking the leftmost locus
# misplaces 33 of them by a median 119 kb (24 land in the wrong 100-kb bin).
# The library's own guide_start pins the intended copy to within 1 bp while the
# runner-up sits ~32.7 kb away, so use it to pick the right locus and fall back
# to the leftmost only when guide_start is unavailable.

cand = gr_mhc.merge(
    summary[["spacer_id", "guide_start"]],
    on="spacer_id",
    how="left"
)

cand["locus_offset"] = (
    cand["start0"] -
    pd.to_numeric(cand["guide_start"], errors="coerce")
).abs()

coords = (
    cand
    .sort_values(["spacer_id", "locus_offset", "start0"], na_position="last")
    .drop_duplicates("spacer_id")
    [["spacer_id", "contig", "start0", "end0", "strand"]]
)

# Add coordinates to guide-level summary
x = summary.merge(
    coords,
    on="spacer_id",
    how="left"
)

# Keep guides with a GRCh38 MHC target

x = x[x["start0"].notna()].copy()

x["start0"] = x["start0"].astype(int)
x["end0"] = x["end0"].astype(int)

# Define lost versus retained

x["lost_in_K562"] = (
    x["exact_match_category"] == "GRCh38 only"
)

lost = x[x["lost_in_K562"]].copy()

print()
print("GRCh38 MHC-targeting guide rows:", len(x))
print("Missing from K562:", len(lost))

# Save actual MHC losses

lost.to_csv(
    OUTDIR / "K562_missing_MHC_guides.tsv",
    sep="\t",
    index=False
)


# ANALYSIS 1: 10-kb distribution

x["bin_start"] = (
    x["start0"] // BIN_SIZE
) * BIN_SIZE

stats = (
    x.groupby("bin_start")
     .agg(
         total_GRCh38_guides=("lost_in_K562", "size"),
         missing_guides=("lost_in_K562", "sum")
     )
     .reset_index()
)

stats["found_guides"] = (
    stats["total_GRCh38_guides"] -
    stats["missing_guides"]
)

stats["loss_fraction"] = (
    stats["missing_guides"] /
    stats["total_GRCh38_guides"]
)

stats["loss_percent"] = (
    stats["loss_fraction"] * 100
)

stats.to_csv(
    OUTDIR / "K562_loss_by_10kb_window.tsv",
    sep="\t",
    index=False
)


# Figure 1: number of missing guides

plt.figure(figsize=(14, 4))

plt.bar(
    stats["bin_start"] / 1_000_000,
    stats["missing_guides"],
    width=BIN_SIZE / 1_000_000
)

plt.xlabel("GRCh38 chr6 position (Mb)")
plt.ylabel("Missing guides per 10 kb")
plt.title(
    "Distribution of GRCh38 MHC targets absent from K562 hap2 MHC"
)

plt.tight_layout()
plt.savefig(
    OUTDIR / "K562_missing_guides_distribution.png",
    dpi=300
)
plt.close()

# Figure 2: fraction lost

plt.figure(figsize=(14, 4))

plt.bar(
    stats["bin_start"] / 1_000_000,
    stats["loss_percent"],
    width=BIN_SIZE / 1_000_000
)

plt.xlabel("GRCh38 chr6 position (Mb)")
plt.ylabel("Targets absent from K562 (%)")
plt.title(
    "Fraction of GRCh38 MHC targets absent from K562"
)

plt.ylim(0, 100)

plt.tight_layout()
plt.savefig(
    OUTDIR / "K562_missing_guides_fraction.png",
    dpi=300
)
plt.close()

# ANALYSIS 1b: found versus missing, side by side
#
# The question these two figures answer: are the apparent hotspots of missing
# guides real jackpots, or do they just track where the library is dense?
#
# At 10 kb there are ~270 populated bins, which is far too many to draw as
# grouped bars, so these use PLOT_BIN_SIZE (100 kb by default). Empty bins
# inside the range are kept so the x-axis stays proportional to the genome.

SERIES_FOUND = "#2a78d6"    # categorical slot 1
SERIES_MISSING = "#eb6834"  # categorical slot 2
INK = "#0b0b0b"
INK_MUTED = "#52514e"

x["plot_bin"] = (
    x["start0"] // PLOT_BIN_SIZE
) * PLOT_BIN_SIZE

pstats = (
    x.groupby("plot_bin")
     .agg(
         total_guides=("lost_in_K562", "size"),
         missing_guides=("lost_in_K562", "sum")
     )
)

# Keep empty bins so gaps in the library read as gaps, not as compression
full_index = np.arange(
    pstats.index.min(),
    pstats.index.max() + PLOT_BIN_SIZE,
    PLOT_BIN_SIZE
)
pstats = pstats.reindex(full_index, fill_value=0)
pstats.index.name = "bin_start"
pstats = pstats.reset_index()

pstats["found_guides"] = (
    pstats["total_guides"] -
    pstats["missing_guides"]
)

# Expected losses if the loss rate were uniform across the MHC. Departures
# from this are what "jackpot" means; agreement means the missing guides are
# just tracking library density.
global_loss_rate = (
    x["lost_in_K562"].sum() /
    len(x)
)

pstats["expected_missing"] = (
    pstats["total_guides"] *
    global_loss_rate
)

pstats["enrichment"] = np.where(
    pstats["expected_missing"] > 0,
    pstats["missing_guides"] / pstats["expected_missing"],
    np.nan
)

# One-sided binomial test per bin: more losses than the global rate predicts?
from scipy.stats import binomtest

pstats["binom_p"] = [
    binomtest(int(m), int(n), global_loss_rate, alternative="greater").pvalue
    if n > 0 else np.nan
    for m, n in zip(pstats["missing_guides"], pstats["total_guides"])
]

# Benjamini-Hochberg across the tested bins
tested = pstats["binom_p"].notna()
ranked = pstats.loc[tested, "binom_p"].rank(method="first")
pstats.loc[tested, "binom_q"] = (
    pstats.loc[tested, "binom_p"] *
    tested.sum() /
    ranked
).clip(upper=1.0)

pstats["bin_mb"] = pstats["bin_start"] / 1e6
pstats.to_csv(
    OUTDIR / "K562_found_vs_missing_by_plot_window.tsv",
    sep="\t",
    index=False
)

bin_mb = PLOT_BIN_SIZE / 1e6
bar_w = bin_mb * 0.42          # leaves a visible gap between the paired bars
centers = pstats["bin_mb"] + bin_mb / 2

# Figure 3: grouped bars, found on the left axis and missing on the right.
#
# NOTE: two independent y-scales. This is the layout requested, but the two
# heights are NOT directly comparable -- the axis ranges alone decide how tall
# each series looks. Read Figure 4 for the actual enrichment question.

fig, ax_found = plt.subplots(figsize=(15, 5))
ax_missing = ax_found.twinx()

b1 = ax_found.bar(
    centers - bar_w / 2,
    pstats["found_guides"],
    width=bar_w,
    color=SERIES_FOUND,
    label="Found in K562 (left axis)"
)

b2 = ax_missing.bar(
    centers + bar_w / 2,
    pstats["missing_guides"],
    width=bar_w,
    color=SERIES_MISSING,
    label="Missing from K562 (right axis)"
)

ax_found.set_xlabel("GRCh38 chr6 position (Mb)", color=INK)
ax_found.set_ylabel(
    f"Guides found in K562 per {PLOT_BIN_SIZE // 1000} kb",
    color=SERIES_FOUND
)
ax_missing.set_ylabel(
    f"Guides missing from K562 per {PLOT_BIN_SIZE // 1000} kb",
    color=SERIES_MISSING
)

ax_found.tick_params(axis="y", colors=SERIES_FOUND)
ax_missing.tick_params(axis="y", colors=SERIES_MISSING)

ax_found.set_axisbelow(True)
ax_found.grid(axis="y", color="#e6e5e1", linewidth=0.8)
for spine in ("top",):
    ax_found.spines[spine].set_visible(False)
    ax_missing.spines[spine].set_visible(False)

ax_found.legend(
    handles=[b1, b2],
    loc="upper left",
    frameon=False,
    labelcolor=INK
)

ax_found.set_title(
    "Guides found vs. missing in K562 hap2 MHC "
    f"(independent y-scales; {PLOT_BIN_SIZE // 1000}-kb bins)",
    color=INK
)

fig.tight_layout()
fig.savefig(
    OUTDIR / "K562_found_vs_missing_grouped_dual_axis.png",
    dpi=300
)
plt.close(fig)

# Figure 4: the same bins on ONE scale -- observed losses against the number
# expected from library density alone. This is the figure that answers whether
# a peak is a jackpot.

fig, ax = plt.subplots(figsize=(15, 5))

ax.bar(
    centers - bar_w / 2,
    pstats["missing_guides"],
    width=bar_w,
    color=SERIES_MISSING,
    label="Missing (observed)"
)

ax.bar(
    centers + bar_w / 2,
    pstats["expected_missing"],
    width=bar_w,
    color=SERIES_FOUND,
    label=f"Missing expected from library density "
          f"(uniform {100 * global_loss_rate:.1f}% loss)"
)

# Mark bins with significantly more losses than density predicts
sig = pstats[pstats["binom_q"] < 0.05]
if len(sig):
    ax.plot(
        sig["bin_mb"] + bin_mb / 2,
        sig[["missing_guides", "expected_missing"]].max(axis=1) + 1.5,
        marker="*",
        linestyle="none",
        color=INK,
        markersize=9,
        label="Excess loss (BH q < 0.05)"
    )

ax.set_xlabel("GRCh38 chr6 position (Mb)", color=INK)
ax.set_ylabel(
    f"Guides missing per {PLOT_BIN_SIZE // 1000} kb",
    color=INK
)
ax.set_axisbelow(True)
ax.grid(axis="y", color="#e6e5e1", linewidth=0.8)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(loc="upper left", frameon=False, labelcolor=INK)
ax.set_title(
    "Are the missing-guide peaks jackpots, or just library density?",
    color=INK
)

fig.tight_layout()
fig.savefig(
    OUTDIR / "K562_missing_observed_vs_expected.png",
    dpi=300
)
plt.close(fig)

print(f"\nGlobal loss rate: {100 * global_loss_rate:.2f}%")
print(f"Bins with excess loss (BH q < 0.05): {len(sig)} of {int(tested.sum())}")
if len(sig):
    print(
        sig[["bin_mb", "total_guides", "missing_guides",
             "expected_missing", "enrichment", "binom_q"]]
        .sort_values("enrichment", ascending=False)
        .to_string(index=False, float_format=lambda v: f"{v:.3g}")
    )

# ANALYSIS 2: clusters of missing guides

lost = lost.sort_values("start0").reset_index(drop=True)

lost["gap_from_previous"] = (
    lost["start0"].diff()
)

lost["new_cluster"] = (
    lost["gap_from_previous"].isna() |
    (lost["gap_from_previous"] > CLUSTER_GAP)
)

lost["cluster_id"] = (
    lost["new_cluster"].cumsum()
)

clusters = (
    lost.groupby("cluster_id")
        .agg(
            start=("start0", "min"),
            end=("end0", "max"),
            n_missing_guides=("spacer_id", "size")
        )
        .reset_index()
)

clusters["span_bp"] = (
    clusters["end"] -
    clusters["start"]
)

clusters["start_mb"] = clusters["start"] / 1e6
clusters["end_mb"] = clusters["end"] / 1e6

clusters.to_csv(
    OUTDIR / "K562_missing_guide_clusters.tsv",
    sep="\t",
    index=False
)


# Print most interesting regions

print("\nHighest-loss 10-kb windows:")

print(
    stats[
        stats["total_GRCh38_guides"] >= 5
    ]
    .sort_values(
        "loss_percent",
        ascending=False
    )
    .head(20)
    .to_string(index=False)
)


print("\nLargest clusters of missing guides:")

print(
    clusters
    .sort_values(
        "n_missing_guides",
        ascending=False
    )
    .head(20)
    .to_string(index=False)
)
