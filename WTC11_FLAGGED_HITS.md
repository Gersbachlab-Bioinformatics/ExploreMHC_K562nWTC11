# WTC11 flagged hit guides: DHSs, genes and target annotations

This audit expands the 42 flagged top-guide slots in the four WTC11 screens.
It uses the existing definition: FRACTEL FDR-corrected p = 0, top three guides
ranked per screen/DHS/Ensembl-gene pair, and `relevant_problem=True`.
A flag identifies a mapping concern for review, not a demonstrated false hit.

**Terminology:** “Slots” means **screen-specific guide–gene entries**. Each slot
is one guide in the top-three list for a particular screen × DHS × gene
combination. The same guide can occupy multiple slots across genes or screens;
slot counts therefore differ from unique-guide counts.

## How many distinct targets are affected?

- **31 unique guide IDs** among 742 unique guides occupying top slots at hits.
- **17 unique DHSs** among 169 hit DHSs (10.1%).
- **11 unique Ensembl genes** among 64 hit genes (17.2%).
- **20 unique DHS–gene pairs** among 243, after collapsing screens (8.2%).
- **27 screen-specific DHS–gene pairs** among 320 (8.4%).
- **42 flagged guide slots** among 960 (4.38%).

An affected gene has at least one flagged guide in a top-three set for a
significant DHS–gene pair. This does not invalidate its other regulatory links.
Guide and gene counts across categories or DHS rows can overlap.

## Promoters versus enhancer-labeled elements

The source is `data/IGVFFI2404DYFG.tsv.gz`. The intended-target classification
is stored in **`genomic_element`**; `intended_target_name` contains the target
identifier/coordinates, and `type` distinguishes targeting/control guides.

**All 17 flagged DHSs are labeled `enhancer`; zero are labeled `promoter`.**
All 31 unique flagged guides carry that enhancer label. There are no missing or
mixed classifications for these DHSs. The join was checked by `guide_id`, and
`intended_target_name` agrees with the screen's `dhs` for every hit slot.
Classification consistency was also checked across all library guides assigned
to each hit DHS, not only the flagged guides.

The metadata does not provide a `distal` category or distance to a transcript's
TSS. Thus these are **17 enhancer-labeled targets**, not 17 independently
confirmed distal elements. A subsequent RefSeq analysis using the first base of
the first exon as the TSS, restricted to RefSeq protein-coding genes and explicitly annotated transcribed
pseudogenes, found **8/17 DHSs within ±1 kb of at least one TSS** and
**9/17 outside that distance**. See [WTC11_REFSEQ_PROXIMITY.md](WTC11_REFSEQ_PROXIMITY.md)
for the GRCh38 RS_2025_08 annotation, per-DHS gene matches and distances.
This positional annotation is separate from the original library labels.

All 169 DHSs in the current significant-hit/top-guide cohort are enhancer-labeled,
so this cohort does not support comparing promoter versus enhancer risk rates.
The full guide library does include promoter-labeled guides (149), but those
are not represented among these hit slots.

## Mapping concerns

| Class | Slots | Unique guides | DHSs | Genes | DHS–gene pairs | Screen-specific pairs |
| --- | --- | --- | --- | --- | --- | --- |
| missing_in_WTC11 | 24 | 17 | 9 | 8 | 10 | 15 |
| multi_WTC11 | 2 | 2 | 2 | 2 | 2 | 2 |
| multi_hg38 | 16 | 12 | 6 | 7 | 8 | 10 |

Of the 17 unique exact-missing guides, four have the low-coverage flag. These
four account for seven slots across three DHSs, three genes and six
screen-specific pairs. The remaining 13 unique missing guides account for
17 slots across six DHSs, five genes and nine screen-specific pairs.

`multi_hg38` means multiple reference MHC loci, not necessarily multiple
functional WTC11 targets. Mapping classes use the existing precedence rules;
these are exclusive class counts, not counts of every independent flag.
The low-coverage flag is based on nearby exact guide matches and does not
establish an assembly gap. One multi_hg38 slot also has that flag; it is not
included in the seven low-coverage *missing* slots.

## Every affected DHS

The gene column names screen response genes, not an independently annotated
promoter owner. Every row below has `genomic_element=enhancer`.

| DHS | Unique flagged guides | Slots | Response genes | Screens | Class |
| --- | --- | --- | --- | --- | --- |
| chr6:29887280-29888879 | 1 | 1 | HLA-A | ipsc.krab | multi_WTC11 |
| chr6:29917733-29918203 | 2 | 2 | POU5F1 | ipsc.krab | missing_in_WTC11 |
| chr6:29924848-29925041 | 2 | 2 | HLA-A;HLA-B | ipsc.krab;ipsc.p300 | multi_hg38 |
| chr6:29934071-29934293 | 3 | 3 | HLA-A | ipsc.krab | missing_in_WTC11 |
| chr6:29941886-29943609 | 4 | 6 | HLA-A | ipsc.krab;ipsc.p300;npc.p300 | missing_in_WTC11 |
| chr6:30006203-30006475 | 1 | 2 | HLA-A | ipsc.krab;npc.krab | multi_hg38 |
| chr6:31270985-31272710 | 2 | 5 | HLA-C | ipsc.krab;ipsc.p300;npc.krab;npc.p300 | missing_in_WTC11 |
| chr6:31395558-31396207 | 1 | 1 | HSPA1B | npc.p300 | missing_in_WTC11 |
| chr6:31396975-31397224 | 1 | 1 | MICA | npc.p300 | multi_hg38 |
| chr6:31497604-31498618 | 1 | 1 | MICB | ipsc.krab | multi_hg38 |
| chr6:31613382-31615931 | 1 | 1 | AIF1 | ipsc.krab | missing_in_WTC11 |
| chr6:31817239-31817580 | 6 | 9 | HSPA1A;HSPA1B | ipsc.krab;npc.krab | multi_hg38 |
| chr6:32521797-32522141 | 1 | 1 | HLA-DRB1 | ipsc.p300 | multi_hg38 |
| chr6:32589565-32590024 | 1 | 1 | HLA-DRB1 | ipsc.p300 | missing_in_WTC11 |
| chr6:32664796-32665308 | 1 | 1 | SOX2 | npc.p300 | missing_in_WTC11 |
| chr6:32763294-32763473 | 1 | 1 | HSPA1B | ipsc.krab | multi_WTC11 |
| chr6:33051226-33051462 | 2 | 4 | HSPA1A;HSPA1B | ipsc.krab | missing_in_WTC11 |

## Every affected gene

| Gene | Affected DHSs | Unique flagged guides | Slots |
| --- | --- | --- | --- |
| SOX2 | 1 | 1 | 1 |
| HLA-DRB1 | 2 | 2 | 2 |
| HSPA1B | 4 | 10 | 10 |
| HSPA1A | 2 | 5 | 5 |
| AIF1 | 1 | 1 | 1 |
| MICB | 1 | 1 | 1 |
| MICA | 1 | 1 | 1 |
| HLA-C | 1 | 2 | 5 |
| POU5F1 | 1 | 2 | 2 |
| HLA-A | 5 | 10 | 13 |
| HLA-B | 1 | 1 | 1 |

## Pairs with all three top guides flagged

These five screen-specific pairs are particularly useful starting points for
review because none of their top three guides escapes the mapping flags.
They represent three DHSs and three response genes (four distinct DHS–gene
pairs when screens are collapsed).

| Screen | DHS | Gene | Flagged | Top slots |
| --- | --- | --- | --- | --- |
| ipsc.krab | chr6:29934071-29934293 | HLA-A | 3 | 3 |
| ipsc.krab | chr6:31817239-31817580 | HSPA1B | 3 | 3 |
| ipsc.krab | chr6:31817239-31817580 | HSPA1A | 3 | 3 |
| ipsc.p300 | chr6:29941886-29943609 | HLA-A | 3 | 3 |
| npc.krab | chr6:31817239-31817580 | HSPA1B | 3 | 3 |

The two HLA-A rows involve exact-missing guides. The three HSPA1A/HSPA1B rows
involve reference multi-mapping. Other guides outside the top-three set may
still map cleanly; this is not an all-guides-per-DHS assessment.

## Reproduce and inspect individual guides

```bash
mamba activate mambaforge
python summarize_wtc11_flagged_hits.py
```

Saved tables under `results/wtc11/`:

- `flagged_hit_unique_guides.tsv`: all 31 guide IDs, target annotation, mapping
  class, low-coverage flag, affected genes/screens, slot count and best rank.
- `flagged_hit_slots_with_metadata.tsv`: all 42 slots, with original guide
  metadata added, including `intended_target_name`, `genomic_element`, `type`
  and `description`.
- `flagged_hit_DHSs.tsv` and `flagged_hit_genes.tsv`: distinct-target summaries.
- `flagged_hit_pairs.tsv`: all 27 affected screen-specific pairs, with flagged
  count out of the actual top-guide count.
- `flagged_hit_overview.tsv`, `flagged_hit_mapping_classes.tsv` and
  `flagged_hit_element_classes.tsv`: aggregate counts and denominators.

Validation checked unique metadata IDs, complete joins, matching DHS identifiers,
consistent element annotations, and exact agreement with the previously saved
42 flagged slots. The existing alignment and screen results were not changed.
