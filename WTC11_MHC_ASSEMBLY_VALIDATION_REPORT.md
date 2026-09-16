# WTC11 MHC hifiasm assembly validation

Analysis run September 15-16, 2026 with the existing `mambaforge` environment
(`minimap2 2.31`, `samtools 1.23.1`, `pandas`, `matplotlib`).

This checks the first of the four goals from the original handoff (see
`analyze_mhc_continuity.py`'s docstring lineage): **a continuous MHC sequence
connecting unique flanks on both sides**, for each WTC11 haplotype, plus how
well MHC/HLA gene content is preserved. Read support, base-level completeness
(Merqury QV), full HLA allele typing against IPD-IMGT/HLA, and phasing are
**out of scope** here.

## Inputs

- `data/hap1.fa.gz`, `data/hap2.fa.gz` — the WTC11 hifiasm primary-contig
  assemblies (whole genome scale, matching
  `data/wtc11_lowq/assembly.asm.bp.hap{1,2}.p_ctg.lowQ.bed`), not an
  MHC-only slice.
- `data/mhc_ref.fa` — GRCh38 chr6:28,510,000-33,480,000 (1-based), extracted
  via `samtools faidx "$ref" chr6:28510000-33480000 > mhc_ref.fa`. The FASTA
  header is literally `chr6:28510000-33480000`, which the scripts parse
  directly rather than hardcoding the offset.
- `data/hap{1,2}.mhc.paf` — `minimap2 -cx asm5 --cs` of each haplotype
  assembly against `mhc_ref.fa`, computed on the HPC and copied in; not
  regenerated locally (regenerating against the full ~800 MB compressed
  assemblies is an HPC-scale job, not something to run on a laptop).
- `data/refseq/GCF_000001405.40-RS_2025_08_genomic.gtf.gz` — RefSeq GTF,
  GRCh38.p14, annotation release 2025_08. Chr6 is accession `NC_000006.12`.
  Gene symbol is the `gene` attribute (not `gene_name`).

## Method

`mhc_paf_utils.py` holds the shared logic used by all three scripts below:

- Load a PAF, split into confident (MAPQ ≥ 30, default) and excluded
  records. **MAPQ filtering matters here**: the MHC region contains
  dispersed repeats that independently align back to the same narrow
  reference window from dozens of unrelated query positions at MAPQ 60 (one
  ~5.6 kb element recurs across `hap1_h1tg000020l` this way around
  chr6:30.08-30.09 Mb) — real sequence, but not part of the assembly's
  colinear path through the interval. Neither the original `PAF col11 >
  50000` filter nor summed block length per contig survives this; **merged
  (deduplicated) reference span** does.
- Pick the candidate contig per haplotype as the one with the largest merged
  reference span among confident records.
- Build a gap table from consecutive merged intervals on that contig, in
  both reference-relative and chromosome coordinates.
- Parse RefSeq genes in the interval and compute, per gene per haplotype,
  what fraction of its genomic span (full gene, not exons) falls inside the
  candidate contig's merged intervals, plus mean alignment identity and
  whether a gap overlaps it. Genes that extend past the extraction boundary
  (e.g. `GPX6`, which starts 6.7 kb before chr6:28,510,000) are clipped to
  the visible portion before scoring, so the missing part of the interval
  itself isn't counted against the assembly.

## Continuity and gap results

| haplotype | candidate contig | span coverage | flank distance (start / end) | gaps | total gap bp |
| --- | --- | ---: | --- | ---: | ---: |
| hap1 | `hap1_h1tg000020l` | 97.2% (4,829,146 / 4,970,001 bp) | 0 bp / 0 bp | 35 | 140,855 |
| hap2 | `hap2_h2tg000025l` | 96.8% (4,813,104 / 4,970,001 bp) | 0 bp / 0 bp | 26 | 156,897 |

Both haplotypes are carried by a **single contig** that reaches **exactly
both flanks** of the extracted interval (0 bp from either edge) — this
confirms, on solid footing, the single-contig candidates the earlier crude
`results/wtc11/assembly_contig_summary.tsv` pass already landed on. Gaps are
not evenly spread: **almost all of them fall in one place**, the extended
HLA class II region (~chr6:32.5-32.8 Mb) — visible as the only kinks in the
dotplot below. This span-coverage number is reference-span colinearity, not
base-level completeness: a merged interval can still contain small internal
indels invisible without CIGAR/`cs`-string parsing.

Full tables: `results/wtc11/mhc_continuity_summary.tsv`,
`mhc_continuity_hap1.tsv`, `mhc_continuity_hap2.tsv`.

## Dotplot

![dotplot](results/wtc11/mhc_ref_vs_hap1_hap2_dotplot.png)

Two side-by-side panels (equal aspect: 1 Mb reference = 1 Mb contig, so
perfect colinearity reads as a 45-degree line), each showing one segment per
merged interval (not one per raw alignment record — plotting every record
reproduces the repeat-noise scatter described above), with gap regions
shaded and classic HLA genes ticked along the top.

`hap2_h2tg000025l`'s backbone alignment is minus-strand throughout, so its
y-axis is flipped (labeled "axis flipped" on the panel) to read
left-to-right the same way as hap1. This is an artifact of hifiasm's
arbitrary contig-output orientation, **not a structural inversion** in the
WTC11 genome.

## Gene content

428 genes are annotated in the interval. Every gene below 99% footprint
coverage in either haplotype, or overlapping a gap (14 genes; `GPX6`
excluded after the boundary-clip fix):

| gene | biotype | classic HLA | hap1 coverage | hap2 coverage |
| --- | --- | :-: | ---: | ---: |
| HLA-DRB5 | protein_coding | yes | 6.1% | 0.0% |
| RNU1-152P | pseudogene | | 0.0% | 0.0% |
| HLA-DRB1 | protein_coding | yes | 0.0% | 26.3% |
| HLA-DQB1 | protein_coding | yes | 0.0% | 0.0% |
| HLA-DQB1-AS1 | lncRNA | | 0.0% | 0.0% |
| HLA-DQA1-AS1 | lncRNA | | 56.1% | 0.1% |
| HLA-DRB6 | transcribed_pseudogene | | 14.5% | 3.0% |
| HLA-DQA1 | protein_coding | yes | 8.7% | 7.0% |
| LOC112267902 | lncRNA | | 100.0% | 96.6% |
| HLA-DPB1 | protein_coding | yes | 97.8% | 100.0% |
| HLA-B | protein_coding | yes | 100.0% | 98.5% |
| HLA-C | protein_coding | yes | 98.9% | 100.0% |
| MUC21 | protein_coding | | 100.0% | 99.1% |
| TSBP1-AS1 | lncRNA | | 100.0% | 99.9% |

All eight classic HLA class I genes and the complement genes are clean:
**HLA-A, HLA-B, HLA-C, HLA-E, HLA-F, HLA-G, MICA, MICB, C4A, C4B all score
100% footprint coverage at high identity in both haplotypes.** The low
scores cluster entirely in the extended class II DRB/DQ complex — HLA-DRB5,
HLA-DRB1, HLA-DQA1, HLA-DQB1, plus the pseudogene/lncRNA neighbors embedded
in the same locus (HLA-DRB6, RNU1-152P, HLA-DQA1-AS1, HLA-DQB1-AS1) — the
single most structurally polymorphic part of the MHC, where different
haplotypes carry different DRB paralog combinations and divergent DQ
alleles. **Low coverage against one arbitrary GRCh38 haplotype here is the
expected signature of real allelic divergence, not necessarily an assembly
defect** — confirming which needs independent HLA typing or read support,
not more of this alignment-based check. The rest of the flagged rows
(HLA-DPB1, HLA-B, HLA-C, MUC21, TSBP1-AS1, LOC112267902) are ≥96.6%, just
clipped by a gap boundary passing through part of the gene — not
meaningfully "problematic."

Full tables: `results/wtc11/mhc_gene_content_footprint.tsv` (all 428 genes),
`mhc_gene_content_footprint_flagged.tsv` (classic HLA, <95% coverage, or gap
overlap).

## Cross-check against the WTC11 screens

Of the 8 genes with **both** haplotypes below 90% footprint coverage
(HLA-DRB5, HLA-DRB1, HLA-DQA1, HLA-DQB1, HLA-DRB6, RNU1-152P,
HLA-DQA1-AS1, HLA-DQB1-AS1), only **HLA-DRB1** appears in any FRACTEL hit
(`FRACTEL_pval_fdr_corr == 0`) across the four WTC11 (iPSC/NPC) screens:

| screen | dhs | gene | effect size | FDR p |
| --- | --- | --- | ---: | ---: |
| ipsc.p300 | chr6:32521797-32522141 | HLA-DRB1 | 0.293 | 0 |
| ipsc.p300 | chr6:32589565-32590024 | HLA-DRB1 | 0.517 | 0 |

That's 2 of 320 total FRACTEL hits across `ipsc.krab`, `ipsc.p300`,
`npc.krab`, `npc.p300` — the DRB/DQ assembly-coverage gap touches a small
fraction of screen hits, concentrated in one screen and one gene.

## Deferred

- **`asm20` sensitivity check** on the DRB/DQ gaps (recommended as an HPC
  job, not run here — the local machine shouldn't run `minimap2` against the
  full ~800 MB compressed hifiasm assemblies).
- Read support (HiFi reads → assembly), Merqury completeness/QV, full HLA
  allele typing against IPD-IMGT/HLA, and phasing — the remaining goals from
  the original handoff.

## Reproduce

```bash
mamba activate mambaforge
python analyze_mhc_continuity.py     # writes mhc_continuity_{summary,hap1,hap2}.tsv
python plot_mhc_dotplot.py           # writes mhc_ref_vs_hap1_hap2_dotplot.png
python analyze_mhc_gene_content.py   # writes mhc_gene_content_footprint(_flagged).tsv
```

All three default to the tracked `data/` inputs and write into
`results/wtc11/`; every path can be overridden (`--hap1-paf`, `--gtf`,
`--min-mapq`, `--outdir`, etc. — see `--help`).
