#!/usr/bin/env python3
"""Ask whether badly-mapping guides are driving the manuscript's screen hits.

For each screen (cell line x effector) the guides targeting a DHS are ranked by
their gRNA-level p-value and the top 3 are taken. Those top-3 slots are then
cross-referenced against the guide mapping lookup, with particular attention to
the DHS-gene pairs FRACTEL calls significant at FDR-corrected p == 0.

Ranking happens per (dhs, gene) pair, because that is the unit both tables
share: grna_glm_nb tests each guide against each gene, and FRACTEL reports one
p-value per DHS-gene pair.
"""

import argparse
import re
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent
DATA = REPO / "data"

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--cell-line", choices=["K562", "WTC11"], default="K562",
                help="host assembly represented by --lookup; WTC11 applies to ipsc/npc")
ap.add_argument("--manuscript", type=Path, default=DATA / "manuscript",
                help="directory of manuscript result tables (default: %(default)s)")
ap.add_argument("--lookup", type=Path, default=REPO / "results" / "guide_mapping_lookup.tsv.gz",
                help="guide mapping lookup (default: %(default)s)")
ap.add_argument("--outdir", type=Path, default=REPO / "results",
                help="where to write results (default: %(default)s)")
ap.add_argument("--top", type=int, default=3,
                help="how many top guides per DHS-gene pair (default: %(default)s)")
ap.add_argument("--screens", nargs="*", default=None, metavar="PREFIX",
                help="only analyse screens whose name starts with one of these "
                     "(e.g. --screens k562). Default: all screens found.")
args = ap.parse_args()
args.outdir.mkdir(parents=True, exist_ok=True)

lookup = pd.read_csv(args.lookup, sep="\t")
if args.cell_line == "WTC11":
    lookup = lookup.rename(columns=lambda c: c.replace("wtc11", "k562").replace("WTC11", "K562"))
    for c in ["mapping_class", "missing_reason"]:
        lookup[c] = lookup[c].str.replace("WTC11", "K562", regex=False)
keep = ["guide_id", "mapping_class", "is_problematic", "missing_reason",
        "missing_reason_ambiguous", "n_hg38_mhc_loci",
        "n_k562_hap1_loci", "n_k562_hap2_loci", "k562_haplotypes",
        "is_hemizygous_K562", "absence_in_assembly_gap"]
lookup = lookup[keep]

screens = sorted(args.manuscript.glob("*.grna_glm_nb_results.txt.gz"))
if not screens:
    raise SystemExit(f"no grna_glm_nb tables under {args.manuscript}")

if args.screens:
    def screen_name(path):
        return re.sub(r"^mhc\.|\.grna_glm_nb_results\.txt\.gz$", "", path.name)
    wanted = tuple(args.screens)
    available = [screen_name(p) for p in screens]
    screens = [p for p in screens if screen_name(p).startswith(wanted)]
    if not screens:
        raise SystemExit(f"no screens match {wanted}; available: {available}")
    missing = [w for w in wanted
               if not any(n.startswith(w) for n in available)]
    if missing:
        print(f"NOTE: no screen matches {missing}; available: {available}\n")

all_top, summary_rows = [], []

for gpath in screens:
    screen = re.sub(r"^mhc\.|\.grna_glm_nb_results\.txt\.gz$", "", gpath.name)
    fpath = gpath.with_name(gpath.name.replace("grna_glm_nb_results",
                                               "dhs_fractel_results"))
    g = pd.read_csv(gpath, sep="\t")
    f = pd.read_csv(fpath, sep="\t")

    g = g.merge(lookup, on="guide_id", how="left")
    g["is_problematic"] = g.is_problematic.fillna(False).astype(bool)

    # Scope each label to the screen it can actually speak to. A guide that is
    # multi-mapping in GRCh38 cuts several sites in any cell line, so that flag
    # always applies. "Missing from K562" only means the guide has no target in
    # a K562 screen -- in the iPSC and NPC screens the host genome is a
    # different one (WTC11), which this repo has not yet compared, so the flag
    # says nothing about whether the guide had a target there.
    is_k562 = screen.startswith("k562") if args.cell_line == "K562" else screen.startswith(("ipsc", "npc"))
    g["screen_is_k562"] = is_k562
    g["relevant_problem"] = g.mapping_class.eq("multi_hg38")
    if is_k562:
        # Absent from both K562 haplotypes = no target at all; multiple sites
        # on one haplotype = more than one cut on that chromosome copy.
        g["relevant_problem"] |= g.mapping_class.isin(["missing_in_K562",
                                                       "multi_K562"])
    # Hemizygous guides do have a target, just on one allele only, so they are
    # reported separately rather than folded into `relevant_problem`.
    # Only count hemizygosity where the other haplotype's assembly actually
    # covers the region -- otherwise the guide is merely in an assembly gap.
    g["hemizygous_here"] = (is_k562
                            & g.is_hemizygous_K562.fillna(False).astype(bool)
                            & ~g.absence_in_assembly_gap.fillna(False).astype(bool))

    # Rank guides within each DHS-gene pair; break p-value ties on effect size.
    g["abs_lfc"] = g.log2fc.abs()
    g = g.sort_values(["dhs", "ensembl_id", "p_value", "abs_lfc"],
                      ascending=[True, True, True, False])
    grp = g.groupby(["dhs", "ensembl_id"], sort=False)
    g["guide_rank"] = grp.cumcount() + 1
    g["n_guides_in_group"] = grp.guide_id.transform("size")

    top = g[g.guide_rank <= args.top].copy()

    # Which DHS-gene pairs does FRACTEL call significant?
    hits = f.loc[f.FRACTEL_pval_fdr_corr == 0, ["dhs", "ensembl_id"]].drop_duplicates()
    hits["fractel_hit"] = True
    top = top.merge(hits, on=["dhs", "ensembl_id"], how="left")
    top["fractel_hit"] = top.fractel_hit.eq(True)
    top["screen"] = screen
    all_top.append(top)

    hit_top = top[top.fractel_hit]
    summary_rows.append(dict(
        screen=screen,
        dhs_gene_pairs=g[["dhs", "ensembl_id"]].drop_duplicates().shape[0],
        fractel_hits=len(hits),
        top_slots=len(top),
        bg_pct_relevant=100 * top.loc[~top.fractel_hit, "relevant_problem"].mean(),
        hit_top_slots=len(hit_top),
        hit_relevant_problem=int(hit_top.relevant_problem.sum()),
        hit_pct_relevant=(100 * hit_top.relevant_problem.mean()
                          if len(hit_top) else float("nan")),
        hit_pairs_with_problem=int(hit_top.groupby(["dhs", "ensembl_id"])
                                   .relevant_problem.any().sum()),
        hit_hemizygous=int(hit_top.hemizygous_here.sum()),
        hit_pairs_all_hemizygous=int(hit_top.groupby(["dhs", "ensembl_id"])
                                     .hemizygous_here.all().sum()),
    ))

top = pd.concat(all_top, ignore_index=True)
summary = pd.DataFrame(summary_rows)

cols = ["screen", "screen_is_k562", "dhs", "gene_symbol", "ensembl_id", "guide_id", "guide_rank",
        "n_guides_in_group", "p_value", "log2fc", "fractel_hit",
        "mapping_class", "is_problematic", "relevant_problem", "missing_reason",
        "missing_reason_ambiguous", "n_hg38_mhc_loci",
        "n_k562_hap1_loci", "n_k562_hap2_loci", "k562_haplotypes",
        "absence_in_assembly_gap", "hemizygous_here"]
def export_labels(frame):
    frame = frame.rename(columns=lambda c: c.replace("k562", args.cell_line.lower()).replace("K562", args.cell_line)).copy()
    for c in ["mapping_class", "missing_reason"]:
        frame[c] = frame[c].str.replace("K562", args.cell_line, regex=False)
    return frame

export_labels(top[cols]).to_csv(args.outdir / "screen_top_guides_mapping.tsv.gz",
                 sep="\t", index=False, compression="gzip")

flagged = top[top.fractel_hit & top.relevant_problem]
export_labels(flagged[cols]).sort_values(["screen", "dhs", "p_value"]).to_csv(
    args.outdir / "fractel_hits_with_problem_guides.tsv", sep="\t", index=False)

summary.to_csv(args.outdir / "screen_problem_guide_summary.tsv",
               sep="\t", index=False)

pd.set_option("display.width", 200)
print("=== per screen ===")
print(summary.to_string(index=False, float_format=lambda v: f"{v:.1f}"))

print(f"\n=== top-{args.top} slots at FRACTEL-significant (p_fdr == 0) pairs ===")
hit_top = top[top.fractel_hit]
print(f"total slots: {len(hit_top):,}   screen-relevant problems: "
      f"{hit_top.relevant_problem.sum()} "
      f"({100 * hit_top.relevant_problem.mean():.1f}%)")
print("\nby mapping_class:")
print(export_labels(hit_top).mapping_class.value_counts().to_string())

if len(flagged):
    print(f"\n=== the {len(flagged)} flagged guide/hit combinations ===")
    print(export_labels(flagged).groupby(["screen", "mapping_class"]).size()
          .rename("n").reset_index().to_string(index=False))
    print("\nDHS-gene pairs affected, by rank of the problematic guide:")
    print(flagged.guide_rank.value_counts().sort_index().rename("n_guides").to_string())
    print(f"\naffected FRACTEL hits: "
          f"{flagged.groupby(['screen','dhs','ensembl_id']).ngroups} of "
          f"{hit_top.groupby(['screen','dhs','ensembl_id']).ngroups}")

print(f"\nwrote {args.outdir}/screen_top_guides_mapping.tsv.gz")
print(f"wrote {args.outdir}/fractel_hits_with_problem_guides.tsv")
print(f"wrote {args.outdir}/screen_problem_guide_summary.tsv")
