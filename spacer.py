import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent
DATA = REPO / "data"


def data_file(name):
    """Prefer data/<name>; fall back to the repo root for the older layout."""
    candidate = DATA / name
    return candidate if candidate.exists() else REPO / name

src = Path(sys.argv[1]) if len(sys.argv) > 1 else data_file("IGVFFI2404DYFG.tsv.gz")
out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else (DATA if DATA.is_dir() else REPO) / "spacers.fa"

df = pd.read_csv(src, sep="\t")

with open(out_path, "w") as out:
    for i, s in enumerate(df["spacer"]):
        out.write(f">guide_{i}\n{s}\n")
