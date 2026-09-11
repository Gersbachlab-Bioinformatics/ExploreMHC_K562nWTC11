# RefSeq TSS proximity of flagged WTC11 DHSs

**8 of 17 flagged DHSs (47.1%) overlap a ±1-kb window
around at least one protein-coding gene or transcribed pseudogene's RefSeq TSS.** The remaining
9 do not meet this positional
criterion. All 17 retain their original library `enhancer` labels in a separate
column; those labels are not a reliable substitute for this TSS-distance analysis.

## Definition and annotation

TSS means the **first base of the first exon in transcript orientation**:
lowest exonic coordinate on the positive strand, highest exonic coordinate on
the negative strand. Single-exon transcripts are included. Any base of the DHS
within an inclusive ±1,000 bp of any TSS qualifies; the DHS midpoint is not used.
A distance of zero means the TSS lies inside the DHS. Each DHS counts once,
regardless of how many transcripts match.

The annotation is **GCF_000001405.40-RS_2025_08 on GRCh38.p14**, the latest GRCh38 release listed
in NCBI's live directory when checked September 10, 2026. All curated and
predicted RefSeq transcript isoforms with exon records belonging to genes with
RefSeq `gene_biotype=protein_coding` or `transcribed_pseudogene` are included on primary chr6 (NC_000006.12).
The gene-level filter applies to both nearest-TSS searches and ±1-kb calls.
All isoforms of those genes are retained, including noncoding transcript variants;
Other noncoding gene classes and other pseudogene classes are excluded.
“Transcribed pseudogene” follows the explicit RefSeq gene biotype, rather than
inferring transcription from any pseudogene with an exon record. Alternative contigs are excluded.
Genes without an annotated transcript/exon structure cannot contribute a TSS.

The nearby RefSeq gene is not necessarily the screen response gene. This is
a positional promoter-proximity label, not proof of a promoter's function or
of the DHS's regulatory target. In particular, a DHS more than 1 kb from every
included TSS is not thereby proven to be a functional distal enhancer.

## Per-DHS results

| DHS | Screen response genes | RefSeq genes within 1 kb | Nearest TSS gene | Nearest TSS distance (bp) |
| --- | --- | --- | --- | --- |
| chr6:29887280-29888879 | HLA-A | HLA-H | HLA-H | 0 |
| chr6:29917733-29918203 | POU5F1 | None | HLA-A | 24329 |
| chr6:29924848-29925041 | HLA-A;HLA-B | None | HLA-A | 17491 |
| chr6:29934071-29934293 | HLA-A | None | HLA-A | 8239 |
| chr6:29941886-29943609 | HLA-A | HLA-A | HLA-A | 0 |
| chr6:30006203-30006475 | HLA-A | HLA-J | HLA-J | 233 |
| chr6:31270985-31272710 | HLA-C | HLA-C | HLA-C | 0 |
| chr6:31395558-31396207 | HSPA1B | None | MICA | 4504 |
| chr6:31396975-31397224 | MICA | None | MICA | 3487 |
| chr6:31497604-31498618 | MICB | MICB | MICB | 0 |
| chr6:31613382-31615931 | AIF1 | AIF1 | AIF1 | 0 |
| chr6:31817239-31817580 | HSPA1A;HSPA1B | None | HSPA1A | 1697 |
| chr6:32521797-32522141 | HLA-DRB1 | None | HLA-DRB5 | 8146 |
| chr6:32589565-32590024 | HLA-DRB1 | HLA-DRB1 | HLA-DRB1 | 0 |
| chr6:32664796-32665308 | SOX2 | None | HLA-DQB1 | 1349 |
| chr6:32763294-32763473 | HSPA1B | HLA-DQB2 | HLA-DQB2 | 59 |
| chr6:33051226-33051462 | HSPA1A;HSPA1B | None | HLA-DPA1 | 22187 |

## Provenance, validation and reproduction

[NCBI release directory](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/GCF_000001405.40-RS_2025_08/)
provides the [GTF](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/GCF_000001405.40-RS_2025_08/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz). Its header states annotation date 08/01/2025 and
GRCh38.p14. Its MD5 matches NCBI's checksum file. Saved directory listings
record the available release versions at retrieval.

Retained 6,988 chr6 transcripts of protein-coding genes and transcribed pseudogenes with exons, including
131 single-exon transcripts. Synthetic checks
cover both strands, single-exon transcripts, DHS interval boundaries and
inclusive 1-kb distance. Complete transcript-level matches preserve RefSeq
accessions, gene IDs, strands and exact site coordinates for inspection.

Reproduce with `mamba activate mambaforge` followed by
`python annotate_flagged_dhs_refseq.py`. Inputs are in `data/refseq/`.
Outputs under `results/wtc11/refseq_proximity/`:

- `flagged_DHS_refseq_proximity.tsv`: all 17 DHSs, proximity calls, matching
  genes/transcript counts and nearest-TSS distances.
- `flagged_DHS_refseq_transcript_matches.tsv`: all matching transcript accessions,
  strand, 0-based and 1-based TSS coordinates, and distances.
- `refseq_chr6_transcript_sites.tsv.gz`: derived chr6 transcript starts for the two included gene biotypes.
- `provenance.json`: release, source URL, hashes, GTF header and counts.
