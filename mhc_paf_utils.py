"""Shared helpers for reading data/{hap}.mhc.paf and computing MHC-interval
continuity metrics. Used by analyze_mhc_continuity.py, plot_mhc_dotplot.py,
and analyze_mhc_gene_content.py so the candidate-contig / gap logic stays in
one place.
"""

import gzip
import hashlib
import re

import pandas as pd

# RefSeq accession for chr6 in the GTF (data/refseq/GCF_000001405.40-RS_2025_08_genomic.gtf.gz),
# GRCh38.p14. The extracted mhc_ref.fa interval is on this chromosome.
REFSEQ_CHR6_ACCESSION = "NC_000006.12"

CLASSIC_HLA_GENES = [
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
    "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
    "HLA-DQA1", "HLA-DQB1", "HLA-DQA2", "HLA-DQB2",
    "HLA-DPA1", "HLA-DPB1",
    "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
    "TAP1", "TAP2", "C4A", "C4B", "MICA", "MICB",
]

PAF_COLUMNS = [
    "contig", "contig_len", "query_start", "query_end", "strand",
    "reference", "ref_len", "ref_start", "ref_end",
    "n_match", "block_len", "mapq",
]


def file_stats(path):
    """Return (n_lines, md5) for a file, read once. Printed by load_paf so a
    partial/mid-transfer read is visible instead of silently producing a
    truncated result (bit us once while building this: a PAF read mid-copy
    looked like a real "assembly missing" finding until re-checked)."""
    md5 = hashlib.md5()
    n_lines = 0
    with open(path, "rb") as fh:
        for line in fh:
            md5.update(line)
            n_lines += 1
    return n_lines, md5.hexdigest()


def load_paf(path, min_mapq=30):
    """Load a minimap2 PAF (first 12 columns), split into confident
    (mapq >= min_mapq) and excluded subsets. Returns (all, confident, excluded).
    """
    n_lines, md5 = file_stats(path)
    print(f"{path}: {n_lines} lines, md5={md5}")

    df = pd.read_csv(path, sep="\t", header=None, usecols=range(12), names=PAF_COLUMNS)
    df["identity_pct"] = 100 * df["n_match"] / df["block_len"]

    confident = df[df["mapq"] >= min_mapq].copy()
    excluded = df[df["mapq"] < min_mapq].copy()
    print(f"  {len(df)} records total: {len(confident)} at MAPQ >= {min_mapq}, "
          f"{len(excluded)} excluded")
    return df, confident, excluded


def parse_reference_interval(reference_name):
    """'chr6:28510000-33480000' -> ('chr6', 28510000, 33480000)."""
    chrom, span = reference_name.split(":")
    start, end = span.split("-")
    return chrom, int(start), int(end)


def to_chrom_pos(ref_pos_0based, interval_start_1based):
    """0-based position within the extracted reference -> 1-based chromosome
    position. interval_start_1based is the extraction interval's start (e.g.
    28510000 for chr6:28510000-33480000)."""
    return ref_pos_0based + interval_start_1based


def merge_with_records(records):
    """records: DataFrame with ref_start/ref_end (and other columns).
    Returns a list of merged-interval dicts, each carrying the record that
    opened it and the record that most recently extended it -- needed later
    to report query coordinates/strand/identity at each gap boundary, which a
    plain (start, end) merge would lose.
    """
    recs = records.sort_values("ref_start").to_dict("records")
    merged = []
    for r in recs:
        s, e = r["ref_start"], r["ref_end"]
        if merged and s <= merged[-1]["end"]:
            if e > merged[-1]["end"]:
                merged[-1]["end"] = e
                merged[-1]["close_rec"] = r
        else:
            merged.append({"start": s, "end": e, "open_rec": r, "close_rec": r})
    return merged


def merged_span(records):
    """Total merged (deduplicated) reference span covered by records."""
    merged = merge_with_records(records)
    return sum(m["end"] - m["start"] for m in merged)


def candidate_contig(confident):
    """Pick the contig with the largest merged reference span among
    confident (already MAPQ-filtered) records. Returns
    (contig_name, span_bp, merged_intervals)."""
    best_contig, best_span, best_merged = None, -1, []
    for contig, sub in confident.groupby("contig"):
        merged = merge_with_records(sub)
        span = sum(m["end"] - m["start"] for m in merged)
        if span > best_span:
            best_contig, best_span, best_merged = contig, span, merged
    return best_contig, best_span, best_merged


def gap_table(merged, interval_start_1based, strand_aware=True):
    """Given merged intervals (from merge_with_records, sorted by start),
    return one row per gap between consecutive merged intervals, with both
    reference-relative and chromosome coordinates."""
    rows = []
    for i, (prev, curr) in enumerate(zip(merged, merged[1:]), start=1):
        prev_r, curr_r = prev["close_rec"], curr["open_rec"]
        ref_gap = curr["start"] - prev["end"]

        query_gap = None
        if strand_aware and prev_r["strand"] == curr_r["strand"]:
            if prev_r["strand"] == "+":
                query_gap = curr_r["query_start"] - prev_r["query_end"]
            else:
                query_gap = prev_r["query_start"] - curr_r["query_end"]

        rows.append({
            "gap_index": i,
            "ref_start_0based": prev["end"],
            "ref_end_0based": curr["start"],
            "chrom_start_1based": to_chrom_pos(prev["end"], interval_start_1based),
            "chrom_end_1based": to_chrom_pos(curr["start"], interval_start_1based),
            "ref_gap_bp": ref_gap,
            "query_gap_bp": query_gap,
            "prev_strand": prev_r["strand"],
            "curr_strand": curr_r["strand"],
            "prev_identity_pct": prev_r["identity_pct"],
            "curr_identity_pct": curr_r["identity_pct"],
            "prev_mapq": prev_r["mapq"],
            "curr_mapq": curr_r["mapq"],
        })
    return pd.DataFrame(rows)


def load_genes(gtf_path, start_1based, end_1based, chrom_accession=REFSEQ_CHR6_ACCESSION):
    """Parse RefSeq GTF 'gene' records overlapping [start_1based, end_1based]
    on chrom_accession. Returns a DataFrame with gene, gene_biotype, start,
    end, strand (1-based, GTF/chromosome coordinates)."""
    rows = []
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != chrom_accession or f[2] != "gene":
                continue
            gstart, gend = int(f[3]), int(f[4])
            if gend < start_1based or gstart > end_1based:
                continue
            attrs = f[8]
            m_gene = re.search(r'gene "([^"]*)"', attrs)
            m_biotype = re.search(r'gene_biotype "([^"]*)"', attrs)
            rows.append({
                "gene": m_gene.group(1) if m_gene else None,
                "gene_biotype": m_biotype.group(1) if m_biotype else None,
                "start": gstart,
                "end": gend,
                "strand": f[6],
            })
    return pd.DataFrame(rows)
