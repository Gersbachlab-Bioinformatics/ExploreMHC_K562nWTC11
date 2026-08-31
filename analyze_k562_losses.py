
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Files

SUMMARY = "/hpc/group/gersbachlab/rr151/Lexi_MHC/k562_mhc_spacer_analysis_19bp_try2/results/guide_exact_match_summary_trim1bp.tsv.gz"
GRCH38 = "/hpc/group/gersbachlab/rr151/Lexi_MHC/k562_mhc_spacer_analysis_19bp_try2/results/grch38_exact_alignments_trim1bp.tsv.gz"

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

# Keep exact primary GRCh38 alignments

gr = gr[gr["NM"] == 0].copy()
if "secondary" in gr.columns:
    gr = gr[~gr["secondary"].astype(bool)]
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

# Get one MHC coordinate per spacer

coords = (
    gr_mhc
    .sort_values(["spacer_id", "start0"])
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
    "K562_missing_MHC_guides.tsv",
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

stats["loss_fraction"] = (
    stats["missing_guides"] /
    stats["total_GRCh38_guides"]
)

stats["loss_percent"] = (
    stats["loss_fraction"] * 100
)

stats.to_csv(
    "K562_loss_by_10kb_window.tsv",
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
    "K562_missing_guides_distribution.png",
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
    "K562_missing_guides_fraction.png",
    dpi=300
)
plt.close()

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
    "K562_missing_guide_clusters.tsv",
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
