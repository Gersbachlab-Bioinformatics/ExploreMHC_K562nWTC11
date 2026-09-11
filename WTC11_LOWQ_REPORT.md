# WTC11 hifiasm lowQ interval audit

## Result

**No assessed guide position overlaps a hifiasm lowQ interval**, including the
expected positions of all three top iPSC–KRAB guides for
`chr6:29934071-29934293`, in both WTC11 haplotypes.

The BEDs therefore provide **no positive evidence** that hifiasm-flagged assembly
inconsistency explains these missing exact matches. Absence of a lowQ flag does
not establish an error-free assembly, and unassessed projections remain unresolved.

## Files and assembly correspondence

Retrieved from the user-supplied archive:
`https://s3.kopah.uw.edu/sjn/web/public/project/g.crawford/asm.tar`.

The archive is 254,164,664,320 bytes. HTTP byte-range extraction transferred
7,431,628 bytes in total, rather than downloading the full archive.

Files saved in `data/wtc11_lowq/`:

- `assembly.asm.bp.hap1.p_ctg.lowQ.bed` (3,657,104 bytes)
- `assembly.asm.bp.hap2.p_ctg.lowQ.bed` (3,764,266 bytes)
- `download_manifest.json`: member offsets, sizes and SHA-256 hashes
- `archive_members.json`: archive listing and offsets

The BED contig names `h1tg000020l` and `h2tg000025l` correspond to FASTA names
`hap1_h1tg000020l` and `hap2_h2tg000025l`. The selected FASTA sequences are
identical to their S records in the existing local haplotype GFAs. The local GFA
file sizes also match their respective archive members. The remote multi-GB GFA
members were not downloaded again for a full checksum comparison.

| Selected contig | Length | lowQ intervals | Union of lowQ bases |
| --- | ---: | ---: | ---: |
| WTC11 hap1 | 61,671,305 | 18 | 23,400 |
| WTC11 hap2 | 32,322,129 | 10 | 15,761 |

Hifiasm documents lowQ BEDs as contig intervals flagged for inconsistency;
`--lowQ` has a documented default threshold of 70%. These are not raw-read
Phred scores or a direct measure of depth. The original run's full version and
command were not established by this download.
[Hifiasm parameter reference](https://hifiasm.readthedocs.io/en/latest/parameter-reference.html)

## How missing guides were located

A missing exact match has no Bowtie2 coordinate to intersect with the BED.
The local GRCh38 MHC sequence was therefore aligned separately to each full
WTC11 contig with minimap2 2.31-r1302:

```bash
minimap2 -x asm20 -c --eqx -t 8 data/wtc11.hap1.mhc_contigs.fa mhc_ref.fa
minimap2 -x asm20 -c --eqx -t 8 data/wtc11.hap2.mhc_contigs.fa mhc_ref.fa
```

The FASTA header `chr6:28510000-33480000` uses a 1-based start: the correct
0-based offset is **28,509,999**, independently confirmed by reproducing all
11,520 cached GRCh38 exact alignment records in this interval from its sequence.
The analysis has 11,253 distinct MHC guide rows, as before.

For each guide, the GRCh38 exact locus nearest its intended guide coordinate
was selected. Primary minimap2 PAF alignments with MAPQ >=20 were used to
project its 19 reference bases through the CIGAR into haplotype coordinates.
Conflicting projections, missing/deleted bases, and noncolinear projections
remain unassessed. If all bases project monotonically across an insertion in
the haplotype, the full projected span is assessed and labeled
`projected_spanning_indel`. BED overlap uses 0-based half-open coordinates.

All `=` CIGAR blocks were checked against the two FASTA sequences, including
reverse-strand handling. As an independent positional check, projections agree
with the unique exact spacer placement for 10,093/10,106 guides in hap1 and
10,273/10,285 in hap2 among guides eligible for that comparison. The remaining
13/12 discordances show that projection is not infallible, particularly when
reference copy choice and the surviving matching copy differ. A projected site
is an inferred corresponding locus, not proof of the correct allele or paralog.

## Missing versus retained guide positions

“Missing” here means absent from **both** WTC11 haplotypes in the existing
exact-match analysis. “Retained” means found in at least one; each haplotype's
projected reference site is assessed independently.

| Cohort | Haplotype | Guide rows | Assessed | lowQ overlaps | Unassessed |
| --- | --- | ---: | ---: | ---: | ---: |
| Missing from both | hap1 | 475 | 424 | 0 | 51 |
| Missing from both | hap2 | 475 | 451 | 0 | 24 |
| Retained in either | hap1 | 10,778 | 10,740 | 0 | 38 |
| Retained in either | hap2 | 10,778 | 10,660 | 0 | 118 |

Of the 475 missing guides, 413 are assessed in both haplotypes, 49 in just one,
and 13 in neither. No exact spacer placements for this MHC cohort overlap lowQ
intervals either. There is consequently no observed lowQ-overlap enrichment
among missing guides; this flag does not discriminate missing from retained
sites in the assessed cohort.

Among the 17 unique missing guides contributing to flagged screen hits, 16 are
assessed in hap1 and all 17 in hap2; none overlaps lowQ.

## The HLA-A DHS chr6:29934071-29934293

These are the three top iPSC–KRAB guides, with inferred haplotype target spans.
Coordinates in the two haplotype columns refer to their respective contigs,
not GRCh38.

| GRCh38 guide coordinates | Strand | Hap1 projected span | Hap2 projected span | lowQ overlap |
| --- | --- | --- | --- | --- |
| chr6:29934202-29934222 | + | 29943917–29943936 | 29196181–29196200 | Neither |
| chr6:29934191-29934211 | − | 29943905–29943924 | 29196193–29196212 | Neither |
| chr6:29934206-29934226 | + | 29943921–29943940 | 29196177–29196196 | Neither |

The nearest lowQ interval is approximately 7.07 Mb away in hap1 and 1.88 Mb
away in hap2. Their absence as exact 19-mers is therefore not supported by a
nearby lowQ annotation. This analysis does not establish the causal sequence
changes or test the guides' activity.

## Outputs and reproduction

The original match counts, `missing_in_WTC11`, and
`absence_in_assembly_gap` annotations were preserved. New annotations are
stored separately under `results/wtc11/lowq/`:

- `guide_lowQ_annotations.tsv.gz`: one row per guide per haplotype; projection
  status, expected position, overlap, distance to lowQ, and overlap at any actual
  exact-match site. Blank overlap means unassessed, not False.
- `flagged_hit_guides_lowQ.tsv`: unique flagged hit guides joined to annotations.
- `chr6_29934071_29934293_lowQ.tsv`: all 20 guides for the HLA-A DHS, both haplotypes.
- `lowQ_missing_vs_retained.tsv`, `lowQ_contig_summary.tsv`.
- `hap1_MHC_contig_lowQ.bed`, `hap2_MHC_contig_lowQ.bed`: selected-contig intervals.
- `hap1_reference_to_hap.paf`, `hap2_reference_to_hap.paf` and minimap2 logs.
- `annotation_provenance.json`: reference offset, hashes and validation counts.

```bash
mamba activate mambaforge
python fetch_wtc11_lowq.py
mkdir -p results/wtc11/lowq
minimap2 -x asm20 -c --eqx -t 8 data/wtc11.hap1.mhc_contigs.fa mhc_ref.fa > results/wtc11/lowq/hap1_reference_to_hap.paf
minimap2 -x asm20 -c --eqx -t 8 data/wtc11.hap2.mhc_contigs.fa mhc_ref.fa > results/wtc11/lowq/hap2_reference_to_hap.paf
python annotate_wtc11_lowq.py
```

Use existing `mambaforge` dependencies; no packages were installed.
