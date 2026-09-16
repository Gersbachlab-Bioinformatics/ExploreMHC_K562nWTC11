# Guide mapping in the K562 MHC, and what it means for the screen hits

How cleanly each CRISPR guide in the IGVF MHC library maps to the K562 genome
versus GRCh38, and whether badly-mapping guides are driving the manuscript's
DHS–gene hits.

**Short answer: they are not.** Of 418 DHS–gene pairs called significant by
FRACTEL across the five screens, 13 have a badly-mapping guide among their top
three, and in the K562 screen specifically it is 3 of 98. The guides that cannot
work in K562 at all are *depleted* among top hits, not enriched — which is what
should happen if the screen is measuring real biology.

Everything below is reproducible from the scripts in this repo; see
[Reproducing](#reproducing).

---

## Method

Guides are compared as their 3′ 19-nt portion (the first base is trimmed), and
only **exact** matches count. Alignment is `bowtie2 --end-to-end -a`, so every
equally-good locus is reported.

Presence/absence is one question; *why* a guide is absent is another. To answer
the second, guides that map to exactly one place in each genome are used as
**anchors**. Walking the anchors in GRCh38 order builds a coordinate map, and
for a missing guide the GRCh38 distance between its flanking anchors is compared
to the corresponding distance in the K562 haplotype:

- **equal distances** → the sequence is present and colinear, so the guide was
  lost to a base change inside the 19-mer;
- **much shorter in K562** → the segment is deleted.

K562 is diploid, so **both haplotype assemblies are used**. A guide is missing
only if it is absent from *both*; one found in a single haplotype is hemizygous
and can cut only one allele.

---

## Two upstream bugs, found and fixed

Both were in inherited code and both silently distorted the results.

**1. Dropping secondary alignments discarded 82% of MHC guides.** With
`bowtie2 -a`, an MHC guide has ~7 equally-good exact hits — chr6 plus the seven
GRCh38 alt MHC haplotype contigs. bowtie2 flags exactly one as primary, chosen
*arbitrarily* among equal scorers, and it lands on chr6 for only 2,008 of 11,368
guides. Filtering on `~secondary` therefore kept a random 18% subsample.

| | Before | After |
|---|---|---|
| MHC-targeting guides | 2,002 | **11,253** |
| Missing from K562 (hap2) | 347 | **658** |
| Loss rate | 17.33% | **5.85%** |

Restricting to `contig == "chr6"` is what selects the primary-assembly locus;
the flag is not a quality signal. Fixed in `analyze_k562_losses.py`. The same
filter remains in `summarize_and_cluster_exact_trim1bp.py`, which is why that
script's `*_exact_match_count` columns are really 0/1 flags — presence/absence
is unaffected, so downstream category calls were always correct.

**2. Multi-locus guides were placed at the wrong copy.** Resolving them with
`drop_duplicates` on a position sort always takes the leftmost copy, which put
33 of 253 guides at the wrong locus (median error 119 kb; 24 landed in the wrong
100-kb bin). The library's own `guide_start` identifies the intended copy to
within 0–1 bp while the runner-up sits ~32.7 kb away, so the choice is
unambiguous. Fixed.

---

## The guide lookup

`results/guide_mapping_lookup.tsv.gz`, one row per guide (12,625):

| `mapping_class` | Guides | Meaning |
|---|---|---|
| `unique_both` | 9,205 | one locus in GRCh38 and in each K562 haplotype |
| `non_targeting` | 1,262 | library controls |
| `hemizygous_assembly_gap` | 703 | single-haplotype, but the other haplotype has no assembly there |
| `hemizygous_K562` | 701 | genuinely present on one haplotype only |
| `missing_in_K562` | 330 | absent from **both** haplotypes |
| `multi_hg38` | 252 | more than one exact locus in the GRCh38 MHC |
| `outside_K562_assembly` | 110 | positive controls on other chromosomes |
| `multi_K562` | 62 | more than one locus within a single haplotype |

`outside_K562_assembly` matters: those 110 guides sit on chr1/2/3/4/11/12 and
are "absent" from K562 only because the assemblies cover the MHC alone. They are
not defective and are excluded from the problem counts.

---

## Why guides are missing: divergence, not deletion

Of the 330 guides absent from both haplotypes:

| `missing_reason` | Guides | Share |
|---|---|---|
| `snv_or_small_indel` | 252 | **76%** |
| `deletion_in_one_haplotype` | 55 | 17% |
| `contig_break` | 13 | 4% |
| `deletion_in_both_haplotypes` | 10 | 3% |

The divergence call is well supported: across the 467 per-haplotype colinear
calls the flanking anchors sit a median 147 bp apart, the median
GRCh38-versus-haplotype length difference is 0 bp, and 83% agree within ±10 bp.
The DNA is present and colinear; the guide simply no longer matches it.

This distinction is actionable. A guide lost to divergence still has a target in
K562 and could be redesigned against the K562 haplotype. Nothing rescues a guide
whose target is deleted.

---

## Multi-mapping guides are real paralogy

253 guides have more than one exact locus in the GRCh38 MHC (243 with two, 7
with three, 3 with five or more). These are genuine segmental duplications, not
reference artifacts: the closest pair is 12.2 kb apart, median separation 32.7
kb, and alt contigs were excluded before counting.

One duplication dominates. About 134 guides form pairs separated by almost
exactly 32,738 bp around 31.98–32.04 Mb — the spacing and position of the RCCX
module (the C4A/C4B tandem repeat). Others include a 12.2 kb pair at
31.81↔31.82 Mb and a 62 kb pair at 32.52↔32.58 Mb.

**They are conserved in K562**: 246 of 253 (97%) still multi-map on at least
one haplotype, so the duplications are ancestral rather than reference-specific.
Only 6 collapse to a single site per haplotype and 1 is absent from both.

They are also the *least* likely guides to be lost — 1.98% versus 5.94% for
single-locus guides — which makes sense, since only one of two copies has to
survive.

---

## The two haplotypes disagree, and it matters

Using hap2 alone (as the previous analysis did) overstates loss by roughly a
factor of two.

| K562 haplotype presence | Guides |
|---|---|
| both | 9,491 |
| hap2 only | 1,104 |
| **neither (truly missing)** | **330** |
| hap1 only | 328 |

**328 guides match hap1 but not hap2** and were previously misclassified as
missing from K562. Counting both haplotypes moves the missing set from 658 to
330.

**But hap1 is materially less complete** — 11 contigs versus 4, 9,819 spacers
matched versus 10,595, and a **~157 kb assembly gap at chr6:31.751–31.907 Mb**
where hap1 coverage falls to 0.00–0.08 while hap2 holds 0.96–0.99. So a
single-haplotype call can mean a real allelic difference *or* missing assembly.
The lookup measures each haplotype's local coverage and sets
`absence_in_assembly_gap`, which splits the 1,432 single-haplotype guides into
728 assembly gaps and 704 plausibly real, and reclassifies 101 of the
"missing from both" the same way.

**Always filter on `absence_in_assembly_gap` before reading a hemizygous or
missing call as biology.**

---

## Do bad guides drive the screen hits?

Guides targeting each DHS are ranked by their gRNA-level p-value within each
(DHS, gene) pair — the unit both manuscript tables share — and the top three are
cross-referenced against the lookup.

Labels are scoped to the screen they can speak to. `multi_hg38` applies
everywhere, since a multi-mapping guide cuts several sites in any cell line.
`missing_in_K562` and `multi_K562` apply **only to the K562 screen**: the iPSC
and NPC screens run on a different genome (WTC11), which this repo has not yet
compared, so those flags say nothing there.

| Screen | FRACTEL hits | Top-3 slots at hits | Flagged | Affected pairs |
|---|---|---|---|---|
| ipsc.krab | 153 | 459 | 9 | 5 |
| ipsc.p300 | 26 | 78 | 2 | 2 |
| **k562.krab** | **98** | **294** | **9** | **3** |
| npc.krab | 80 | 240 | 4 | 2 |
| npc.p300 | 61 | 183 | 1 | 1 |
| **total** | **418** | **1,254** | **25 (2.0%)** | **13** |

2.0% at hits versus a 2.1–2.2% background: no enrichment. Badly-mapping guides
are not driving the hits.

### The K562 screen validates itself

Among the 294 top-3 slots at K562 FRACTEL hits, guides absent from both K562
haplotypes appear **zero times**, against 8.5 expected from their 2.90%
background rate (binomial p = 1.8e-04). A guide with no target in the host
genome never rises to the top at a real hit — exactly as it should not. This
independently corroborates both the mapping calls and the screen.

### The 13 pairs that are affected

Six have **all three** top guides multi-mapping, so the DHS→gene assignment
cannot distinguish the paralogs:

| Screen | DHS | Gene | Flagged |
|---|---|---|---|
| ipsc.krab | chr6:31817239-31817580 | HSPA1A | 3/3 |
| ipsc.krab | chr6:31817239-31817580 | HSPA1B | 3/3 |
| k562.krab | chr6:31817239-31817580 | HSPA1A | 3/3 |
| npc.krab | chr6:31817239-31817580 | HSPA1B | 3/3 |
| k562.krab | chr6:32011727-32012212 | TNXB | 3/3 |
| k562.krab | chr6:32012313-32013264 | TNXB | 3/3 |

The HSPA1A/HSPA1B DHS is the 12.2 kb duplication above, and it appears assigned
to *both* paralogs in different screens — the signature of guides that cut both.
The TNXB DHSs sit in the RCCX module. The remaining seven pairs (HLA-A ×3,
HLA-B, HLA-DRB1, MICA, MICB) have one flagged guide out of three, which is much
less concerning.

Two further K562 hits — EGFL8 and PPT2, both at chr6:32152952-32155097 — have
all three top guides hemizygous after gap filtering, but hap1 coverage there is
0.90, so treat them as uncertain rather than established.

---

## Caveats

- **Deletion sizes are not trustworthy.** Every deletion interval falls where
  two contigs of the same haplotype both represent the region (Seq3/Seq4 overlap
  32.476–32.665 Mb in hap2). Flanking anchors can then be drawn from different
  contigs and manufacture a spurious length difference. These carry
  `missing_reason_ambiguous`; 270 of 330 missing guides are so flagged. Settling
  them needs the assemblies' own contig/gap records.
- **Δ ≈ 0 shows no *net* indel** between flanking anchors, not that the guide's
  exact 19 bp is unchanged; a deletion plus compensating insertion would cancel.
  Point divergence is the parsimonious reading, not a proof.
- **Gene identities are inferred, not derived.** `putative_target_genes` is empty
  throughout — the library targets coordinate-named CREs — so RCCX/C4, HSPA1A/B
  and class II assignments come from MHC architecture and spacing, not from
  these tables.
- **One missing guide is unbracketed.** The leftmost MHC guide in the library has
  no anchor on its left, so its reason is inferred from local coverage (0.81 /
  0.83, sequence present) rather than a measured span. Flagged with
  `missing_reason_unbracketed`.
- **`analyze_k562_losses.py` is still hap2-only.** Its loss-distribution figures
  and the 5.85% rate predate the hap1 alignment; the lookup's 330 supersedes its
  658. The two are not contradictory — they answer different questions — but do
  not quote them side by side.
- **WTC11 is analysed separately.** See [WTC11_REPORT.md](WTC11_REPORT.md) for
  the same guide-mapping question against both WTC11 haplotypes (iPSC/NPC
  screens), and
  [WTC11_MHC_ASSEMBLY_VALIDATION_REPORT.md](WTC11_MHC_ASSEMBLY_VALIDATION_REPORT.md)
  for whether the WTC11 hifiasm assembly itself reconstructs the MHC well.
  Neither result should be assumed to transfer to K562 or vice versa.

---

## Reproducing

```bash
python3 build_guide_mapping_lookup.py          # -> results/guide_mapping_lookup.tsv.gz
python3 check_problem_guides_in_screens.py     # all five screens
python3 check_problem_guides_in_screens.py --screens k562
```

See [README.md](README.md) for the data layout and for adding a further
haplotype assembly.

## Outputs

| File | Contents |
|---|---|
| `results/guide_mapping_lookup.tsv.gz` | per-guide mapping class, haplotype counts, missing reason, gap flags |
| `results/screen_top_guides_mapping.tsv.gz` | every top-3 guide in every screen with its class and FRACTEL status |
| `results/fractel_hits_with_problem_guides.tsv` | just the affected FRACTEL hits |
| `results/screen_problem_guide_summary.tsv` | one row per screen |
