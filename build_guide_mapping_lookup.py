#!/usr/bin/env python3
"""Build a per-guide lookup of how cleanly each guide maps in GRCh38 and K562.

K562 is represented by BOTH haplotype assemblies. A guide counts as present in
K562 if it matches either haplotype, so "missing" means absent from both; a
guide found in only one haplotype is flagged hemizygous, since it can only cut
one of the two alleles.

Every guide gets a `mapping_class`; guides absent from both haplotypes also get
a `missing_reason` inferred from the local colinearity of GRCh38 and each
haplotype. See the README for the method and its caveats.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
DATA = REPO / "data"


def data_file(name):
    """Prefer data/<name>; fall back to the repo root for the older layout."""
    candidate = DATA / name
    return candidate if candidate.exists() else REPO / name


ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--cell-line", choices=["K562", "WTC11"], default="K562",
                help="label used in output columns and mapping classes")
ap.add_argument("--summary", type=Path,
                default=data_file("guide_exact_match_summary_trim1bp.tsv.gz"),
                help="per-guide exact-match summary (default: %(default)s)")
ap.add_argument("--grch38", type=Path,
                default=data_file("grch38_exact_alignments_trim1bp.tsv.gz"),
                help="GRCh38 exact alignments, all hits (default: %(default)s)")
ap.add_argument("--k562-hap2", "--hap2", type=Path,
                default=None,
                help="K562 hap2 exact alignments (default: %(default)s)")
ap.add_argument("--k562-hap1", "--hap1", type=Path,
                default=None,
                help="K562 hap1 exact alignments (default: %(default)s)")
ap.add_argument("--outdir", type=Path, default=REPO / "results",
                help="where to write the lookup (default: %(default)s)")
args = ap.parse_args()
if args.cell_line != "K562" and (args.k562_hap1 is None or args.k562_hap2 is None):
    ap.error("--cell-line WTC11 requires explicit --hap1 and --hap2 alignment tables")
args.k562_hap1 = args.k562_hap1 or data_file("k562_hap1_exact_alignments_trim1bp.tsv.gz")
args.k562_hap2 = args.k562_hap2 or data_file("k562_mhc_exact_alignments_trim1bp.tsv.gz")
args.outdir.mkdir(parents=True, exist_ok=True)

MHC_CHR, MHC_START, MHC_END = "chr6", 28_000_000, 34_000_000
INDEL_TOL = 500          # bp; smaller length differences count as "small"


def read_hits(path):
    """Exact, non-supplementary alignments.

    The `secondary` flag is deliberately ignored: with bowtie2 -a a guide has
    one primary hit chosen arbitrarily among equally-good ones, so filtering on
    it discards real loci at random.
    """
    d = pd.read_csv(path, sep="\t")
    return d[(d.NM == 0) & (~d.supplementary.astype(bool))]


summary = pd.read_csv(args.summary, sep="\t")
gr = read_hits(args.grch38)
haps = {"hap1": read_hits(args.k562_hap1), "hap2": read_hits(args.k562_hap2)}

# GRCh38 loci on the primary assembly inside the MHC. Alt contigs
# (chr6_GL0002*_alt) are excluded: they are alternate representations of the
# same locus, not additional target sites.
gr_mhc = gr[(gr.contig == MHC_CHR) &
            (gr.start0 >= MHC_START) & (gr.start0 < MHC_END)]

d = summary.copy()
counts = {
    "n_hg38_mhc_loci": gr_mhc.groupby("spacer_id").size(),
    "n_hg38_offtarget_loci": (gr[(gr.contig != MHC_CHR) &
                                 (~gr.contig.str.contains("_alt"))]
                              .groupby("spacer_id").size()),
    "n_k562_hap1_loci": haps["hap1"].groupby("spacer_id").size(),
    "n_k562_hap2_loci": haps["hap2"].groupby("spacer_id").size(),
}
for name, col in counts.items():
    d = d.merge(col.rename(name), left_on="spacer_id", right_index=True, how="left")
    d[name] = d[name].fillna(0).astype(int)

d["n_k562_loci"] = d.n_k562_hap1_loci + d.n_k562_hap2_loci
d["n_k562_haplotypes_hit"] = ((d.n_k562_hap1_loci > 0).astype(int) +
                              (d.n_k562_hap2_loci > 0).astype(int))
# Multiplicity that matters for cutting is per-haplotype: two sites on one
# haplotype means two cuts on that chromosome copy, whereas one site on each
# haplotype is simply the normal diploid state.
d["max_loci_per_haplotype"] = d[["n_k562_hap1_loci", "n_k562_hap2_loci"]].max(axis=1)

has_hg38 = d.GRCh38_exact_match_count > 0

d["mapping_class"] = np.select(
    [
        ~d.targeting.astype(bool),
        ~has_hg38,
        (d.n_k562_haplotypes_hit == 0) & (d.n_hg38_mhc_loci == 0),
        d.n_k562_haplotypes_hit == 0,
        d.n_hg38_mhc_loci > 1,
        d.max_loci_per_haplotype > 1,
        d.n_k562_haplotypes_hit == 1,
        d.n_hg38_mhc_loci == 1,
    ],
    [
        "non_targeting",
        "no_hg38_match",
        "outside_K562_assembly",
        "missing_in_K562",
        "multi_hg38",
        "multi_K562",
        "hemizygous_K562",
        "unique_both",
    ],
    default="other",
)

d["is_multi_hg38"] = d.n_hg38_mhc_loci > 1
d["is_multi_K562"] = d.max_loci_per_haplotype > 1
d["is_hemizygous_K562"] = (d.n_k562_haplotypes_hit == 1) & (d.n_hg38_mhc_loci >= 1)
d["k562_haplotypes"] = np.select(
    [d.n_k562_haplotypes_hit == 2,
     (d.n_k562_hap1_loci > 0) & (d.n_k562_hap2_loci == 0),
     (d.n_k562_hap2_loci > 0) & (d.n_k562_hap1_loci == 0)],
    ["both", "hap1_only", "hap2_only"], default="neither")

# ---------------------------------------------------------------- why missing?
#
# Anchors map GRCh38 onto one haplotype: guides with exactly one locus in each.
# For a guide absent from that haplotype, compare the GRCh38 distance between
# its flanking anchors to the corresponding haplotype distance. Equal distances
# mean the sequence is present and colinear, so the guide was lost to a
# base-level change; a much shorter haplotype distance means a deletion.

g1 = gr_mhc.groupby("spacer_id").filter(lambda x: len(x) == 1)
gr_first = gr_mhc.sort_values(["spacer_id", "start0"]).drop_duplicates("spacer_id")
miss_pos = dict(zip(gr_first.spacer_id, gr_first.start0))

# ------------------------------------------------- assembly gaps vs. real loss
#
# Measure local exact-guide sequence presence, not assembly/base coverage.
# A low fraction can reflect divergence, deletion, or incomplete assembly;
# it does not establish a physical gap. This is recomputed from the supplied
# haplotype tables (including WTC11); no K562 coverage values are reused.
# The historical output name `absence_in_assembly_gap` is retained for compatibility.

COV_WINDOW = 25_000      # bp either side
COV_MIN = 0.5            # threshold for low local exact-guide match fraction

anchored = g1[["spacer_id", "start0"]].sort_values("start0")
anchor_pos = anchored.start0.values


def local_coverage(hits):
    """Fraction of nearby GRCh38-anchored guides that this haplotype matches."""
    present = anchored.spacer_id.isin(set(hits.spacer_id)).to_numpy().astype(float)
    csum = np.concatenate([[0.0], np.cumsum(present)])

    def frac(pos):
        lo = np.searchsorted(anchor_pos, pos - COV_WINDOW, "left")
        hi = np.searchsorted(anchor_pos, pos + COV_WINDOW, "right")
        n = hi - lo
        return (csum[hi] - csum[lo]) / n if n else np.nan

    return frac


cov_fn = {h: local_coverage(hits) for h, hits in haps.items()}


def classify(hits, spacer_ids, coverage_fn):
    """Explain, per guide, why it is absent from this haplotype."""
    k1 = hits.groupby("spacer_id").filter(lambda x: len(x) == 1)
    anchors = (g1[["spacer_id", "start0", "strand"]]
               .rename(columns={"start0": "gpos", "strand": "gstr"})
               .merge(k1[["spacer_id", "contig", "start0", "strand"]]
                      .rename(columns={"contig": "kctg", "start0": "kpos",
                                       "strand": "kstr"}), on="spacer_id")
               .sort_values("gpos").reset_index(drop=True))

    # GRCh38 stretches covered by more than one contig of this haplotype.
    # Deletion calls inside them are unreliable, because the flanking anchors
    # can come from two contigs that both represent the region. This must use
    # ALL hits: a guide in a redundant region matches two contigs and is
    # therefore absent from the anchor set that would reveal the overlap.
    coverage = (g1[["spacer_id", "start0"]]
                .merge(hits[["spacer_id", "contig"]], on="spacer_id")
                .rename(columns={"start0": "gpos", "contig": "kctg"}))
    spans = coverage.groupby("kctg").gpos.agg(["min", "max"])
    overlaps = [(max(ra["min"], rb["min"]), min(ra["max"], rb["max"]))
                for a, ra in spans.iterrows() for b, rb in spans.iterrows()
                if a < b and max(ra["min"], rb["min"]) < min(ra["max"], rb["max"])]

    gpos_arr = anchors.gpos.values
    out = {}
    for sid in spacer_ids:
        pos = miss_pos.get(sid)
        if pos is None:
            out[sid] = ("no_hg38_mhc_locus", "", False, False)
            continue
        i = np.searchsorted(gpos_arr, pos)
        if i == 0 or i >= len(gpos_arr):
            # No anchor on one side -- the guide sits at the very edge of the
            # library's span -- so the two-sided span comparison cannot run.
            # Fall back to local coverage: if this haplotype assembles the
            # surrounding region, the sequence is there and the guide was lost
            # to a base-level change. Flagged `unbracketed` because it rests on
            # coverage rather than on a measured span.
            cov = coverage_fn(pos)
            if cov is not None and not np.isnan(cov) and cov >= COV_MIN:
                out[sid] = ("snv_or_small_indel",
                            f"no flanking anchor; local coverage {cov:.2f}",
                            False, True)
            else:
                out[sid] = ("no_flanking_anchor",
                            f"local coverage {cov:.2f}" if cov == cov else "",
                            False, True)
            continue
        L, R = anchors.iloc[i - 1], anchors.iloc[i]
        gspan = int(R.gpos - L.gpos)
        amb = any(lo <= pos <= hi for lo, hi in overlaps)
        if L.kctg != R.kctg:
            out[sid] = ("contig_break", f"flanked by {L.kctg} and {R.kctg}",
                        True, False)
            continue
        kspan = int(abs(R.kpos - L.kpos))
        delta = gspan - kspan
        if (L.gstr == L.kstr) != (R.gstr == R.kstr):
            reason = "orientation_change"
        elif delta > INDEL_TOL:
            reason = "deletion_in_K562"
        elif delta < -INDEL_TOL:
            reason = "insertion_in_K562"
        else:
            reason = "snv_or_small_indel"
        out[sid] = (reason,
                    f"GRCh38 span {gspan} bp, hap span {kspan} bp, delta {delta} bp",
                    amb, False)
    return out


missing_ids = list(d.loc[d.mapping_class == "missing_in_K562", "spacer_id"])
per_hap = {h: classify(hits, missing_ids, cov_fn[h])
           for h, hits in haps.items()}

for h in haps:
    d[f"missing_reason_{h}"] = d.spacer_id.map(
        {k: v[0] for k, v in per_hap[h].items()}).fillna("")

# Combined call across the two haplotypes. If either still carries the sequence
# colinearly, the guide was lost to base-level divergence rather than structure;
# a deletion only counts as structural when both haplotypes agree.
r1, r2 = d["missing_reason_hap1"], d["missing_reason_hap2"]
d["missing_reason"] = np.select(
    [(r1 == "") & (r2 == ""),
     (r1 == "snv_or_small_indel") | (r2 == "snv_or_small_indel"),
     (r1 == "deletion_in_K562") & (r2 == "deletion_in_K562"),
     (r1 == "deletion_in_K562") | (r2 == "deletion_in_K562"),
     (r1 == "insertion_in_K562") | (r2 == "insertion_in_K562"),
     (r1 == "orientation_change") | (r2 == "orientation_change"),
     (r1 == "contig_break") | (r2 == "contig_break")],
    ["", "snv_or_small_indel", "deletion_in_both_haplotypes",
     "deletion_in_one_haplotype", "insertion_in_K562", "orientation_change",
     "contig_break"],
    default=r2.where(r2 != "", r1))

d["missing_reason_unbracketed"] = d.spacer_id.map(
    {k: (per_hap["hap1"][k][3] and per_hap["hap2"][k][3])
     for k in missing_ids}).eq(True)
d["missing_reason_ambiguous"] = d.spacer_id.map(
    {k: (per_hap["hap1"][k][2] or per_hap["hap2"][k][2])
     for k in missing_ids}).eq(True)
d["missing_reason_detail"] = d.spacer_id.map(
    {k: f"hap1: {per_hap['hap1'][k][1]} | hap2: {per_hap['hap2'][k][1]}"
     for k in missing_ids}).fillna("")

guide_pos = d.spacer_id.map(miss_pos)
for h in haps:
    d[f"{h}_local_coverage"] = [
        cov_fn[h](p) if pd.notna(p) else np.nan for p in guide_pos]

# Coverage of the haplotype(s) the guide is absent from. When it is absent from
# both, the question is whether the BETTER-covered haplotype covers the region:
# use the better fraction; low fractions in both trigger the heuristic flag.
h1c, h2c = d.hap1_local_coverage, d.hap2_local_coverage
absent_cov = np.select(
    [(d.n_k562_hap1_loci == 0) & (d.n_k562_hap2_loci == 0),
     d.n_k562_hap1_loci == 0,
     d.n_k562_hap2_loci == 0],
    [np.fmax(h1c, h2c), h1c, h2c],
    default=np.nan)
d["absent_haplotype_coverage"] = absent_cov

# Flag low local exact-guide presence; this does not distinguish assembly gaps
# from allelic divergence or biological deletions.
d["absence_in_assembly_gap"] = (
    (d.n_k562_haplotypes_hit < 2) &
    (d.n_hg38_mhc_loci >= 1) &
    (pd.Series(absent_cov, index=d.index) < COV_MIN)
).fillna(False)

d.loc[d.mapping_class.eq("hemizygous_K562") & d.absence_in_assembly_gap,
      "mapping_class"] = "hemizygous_assembly_gap"

d["is_problematic"] = d.mapping_class.isin(
    ["missing_in_K562", "multi_hg38", "multi_K562", "hemizygous_K562"])

cols = ["guide_id", "spacer_id", "spacer", "targeting",
        "guide_chr", "guide_start", "guide_end", "strand",
        "intended_target_name", "exact_match_category",
        "n_hg38_mhc_loci", "n_hg38_offtarget_loci",
        "n_k562_hap1_loci", "n_k562_hap2_loci", "n_k562_loci",
        "n_k562_haplotypes_hit", "max_loci_per_haplotype", "k562_haplotypes",
        "mapping_class", "is_multi_hg38", "is_multi_K562",
        "is_hemizygous_K562", "is_problematic",
        "hap1_local_coverage", "hap2_local_coverage",
        "absent_haplotype_coverage", "absence_in_assembly_gap",
        "missing_reason", "missing_reason_hap1", "missing_reason_hap2",
        "missing_reason_ambiguous", "missing_reason_unbracketed",
        "missing_reason_detail"]
out = d[cols].sort_values(["guide_chr", "guide_start"])
dest = args.outdir / "guide_mapping_lookup.tsv.gz"
# Keep the historical K562 schema by default; label other assemblies at export.
export = out.rename(columns=lambda c: c.replace("k562", args.cell_line.lower())
                    .replace("K562", args.cell_line)).copy()
for c in ["mapping_class", "missing_reason", "missing_reason_hap1", "missing_reason_hap2"]:
    export[c] = export[c].str.replace("K562", args.cell_line, regex=False)
export.to_csv(dest, sep="\t", index=False, compression="gzip")

print(f"wrote {dest}  ({len(out):,} guides)\n")
print("mapping_class:")
print(export.mapping_class.value_counts().to_string())
print(f"\n{args.cell_line} haplotype presence (guides with a GRCh38 MHC locus):")
print(export.loc[export.n_hg38_mhc_loci >= 1, f"{args.cell_line.lower()}_haplotypes"].value_counts().to_string())
print(f"\nproblematic guides: {out.is_problematic.sum():,}")
print("\nmissing_reason (absent from BOTH haplotypes):")
mr = out[out.mapping_class == "missing_in_K562"]
print(mr.missing_reason.str.replace("K562", args.cell_line, regex=False).value_counts().to_string())
print(f"\n  flagged ambiguous (inside a contig overlap): "
      f"{mr.missing_reason_ambiguous.sum()}")
print(f"  absent where the assembly does not cover the region: "
      f"{mr.absence_in_assembly_gap.sum()}")
print("\nsingle-haplotype guide matches: local sequence-coverage heuristic")
hz = out[out.is_hemizygous_K562]
print(f"  total single-haplotype: {len(hz):,}")
print(f"    low local coverage : {hz.absence_in_assembly_gap.sum():,}")
print(f"    higher coverage    : {(~hz.absence_in_assembly_gap).sum():,}")
