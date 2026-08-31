import pandas as pd

df = pd.read_csv("/hpc/group/gersbachlab/rr151/Lexi_MHC/IGVFFI2404DYFG.tsv.gz", sep="\t")

with open("/hpc/group/gersbachlab/rr151/Lexi_MHC/spacers.fa","w") as out:
    for i,s in enumerate(df["spacer"]):
        out.write(f">guide_{i}\n{s}\n")
