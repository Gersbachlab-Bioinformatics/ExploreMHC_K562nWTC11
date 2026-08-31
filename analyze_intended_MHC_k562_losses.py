#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

src = Path(sys.argv[1])
out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(".")
out.mkdir(parents=True, exist_ok=True)

MHC_CHR="chr6"; MHC_START=28_000_000; MHC_END=34_000_000
BIN_SIZE=10_000; CLUSTER_GAP=10_000

df=pd.read_csv(src, sep="\t", compression="infer")
mhc=((df.intended_target_chr==MHC_CHR) &
     (pd.to_numeric(df.intended_target_start,errors="coerce")>=MHC_START) &
     (pd.to_numeric(df.intended_target_end,errors="coerce")<=MHC_END))
x=df.loc[mhc].copy()
x["guide_start"]=pd.to_numeric(x.guide_start,errors="coerce")
x["guide_end"]=pd.to_numeric(x.guide_end,errors="coerce")
x["lost_in_K562"]=x.exact_match_category.eq("GRCh38 only")
x.to_csv(out/"intended_MHC_guide_summary_trim1bp.tsv.gz",sep="\t",index=False,compression="gzip")

lost=x.loc[x.lost_in_K562].copy().sort_values(["guide_chr","guide_start","guide_end"])
lost.to_csv(out/"K562_missing_intended_MHC_guides_trim1bp.tsv",sep="\t",index=False)

v=x.dropna(subset=["guide_start"]).copy()
v["bin_start"]=(v.guide_start.astype(int)//BIN_SIZE)*BIN_SIZE
s=v.groupby("bin_start",as_index=False).agg(
    total_intended_MHC_guides=("lost_in_K562","size"),
    missing_in_K562=("lost_in_K562","sum"))
s["retained_in_K562"]=s.total_intended_MHC_guides-s.missing_in_K562
s["loss_fraction"]=s.missing_in_K562/s.total_intended_MHC_guides
s["loss_percent"]=100*s.loss_fraction
s["bin_end"]=s.bin_start+BIN_SIZE
s=s[["bin_start","bin_end","total_intended_MHC_guides","retained_in_K562",
     "missing_in_K562","loss_fraction","loss_percent"]]
s.to_csv(out/"K562_intended_MHC_loss_by_10kb_window_trim1bp.tsv",sep="\t",index=False)

c=lost.dropna(subset=["guide_start","guide_end"]).copy()
c["guide_start"]=c.guide_start.astype(int); c["guide_end"]=c.guide_end.astype(int)
c=c.sort_values("guide_start").reset_index(drop=True)
c["gap_from_previous_bp"]=c.guide_start.diff()
c["new_cluster"]=c.gap_from_previous_bp.isna()|(c.gap_from_previous_bp>CLUSTER_GAP)
c["cluster_id"]=c.new_cluster.cumsum().astype(int)
z=c.groupby("cluster_id",as_index=False).agg(
    start0=("guide_start","min"),end0=("guide_end","max"),
    n_missing_guides=("spacer_id","size"),
    first_guide_id=("guide_id","first"),last_guide_id=("guide_id","last"))
z["span_bp"]=z.end0-z.start0
z["start_Mb"]=z.start0/1e6; z["end_Mb"]=z.end0/1e6
z["classification"]=np.where(z.n_missing_guides==1,"isolated","cluster")
z.to_csv(out/"K562_missing_intended_MHC_guide_clusters_trim1bp.tsv",sep="\t",index=False)
c.to_csv(out/"K562_missing_intended_MHC_guides_with_clusters_trim1bp.tsv",sep="\t",index=False)

plt.figure(figsize=(14,4.5))
plt.bar(s.bin_start/1e6,s.missing_in_K562,width=BIN_SIZE/1e6)
plt.xlabel("GRCh38 chr6 position (Mb)"); plt.ylabel("Missing intended MHC guides per 10 kb")
plt.title("Distribution of intended GRCh38 MHC targets absent from K562 hap2 MHC")
plt.tight_layout(); plt.savefig(out/"K562_intended_MHC_missing_guides_distribution_trim1bp.png",dpi=300); plt.close()

plt.figure(figsize=(14,4.5))
plt.bar(s.bin_start/1e6,s.loss_percent,width=BIN_SIZE/1e6)
plt.xlabel("GRCh38 chr6 position (Mb)"); plt.ylabel("Intended MHC targets absent from K562 (%)")
plt.ylim(0,100); plt.title("Fraction of intended GRCh38 MHC targets absent from K562 hap2 MHC")
plt.tight_layout(); plt.savefig(out/"K562_intended_MHC_loss_fraction_trim1bp.png",dpi=300); plt.close()

print(x.exact_match_category.value_counts().to_string())
print(f"\nIntended MHC guides: {len(x):,}")
print(f"Missing in K562: {len(lost):,} ({100*len(lost)/len(x):.2f}%)")
