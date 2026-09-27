# False-Positive Spliced Alignments Under Extreme gDNA Contamination

**Status:** investigation complete, fixes not yet implemented
**Date:** 2026-08-03
**Index case:** `mctp_LBX0077_SI_44263_HKWNYDRX7` (human, GRCh38 / gencode v46)
**Companion documents:**
[GDNA_SPLICE_ARTIFACTS_ALIGNABLE.md](GDNA_SPLICE_ARTIFACTS_ALIGNABLE.md) ·
[GDNA_SPLICE_ARTIFACTS_STAR.md](GDNA_SPLICE_ARTIFACTS_STAR.md) ·
[GDNA_SPLICE_ARTIFACTS_RIGEL.md](GDNA_SPLICE_ARTIFACTS_RIGEL.md)

---

## 1. Summary

In libraries with extreme genomic-DNA contamination, unspliced gDNA fragments are
being given spliced alignments by STAR and are surviving rigel's splice-junction
blacklist. In the index case (96.4% gDNA) **roughly 43% of the spliced fragments
rigel treats as RNA carry an alignment-level artifact signature.**

Two independent estimators agree:

| estimator | method | result |
|---|---|---|
| strand mixture | observed strand specificity 0.7856 vs 0.991–0.999 in all 20 same-flowcell siblings; two-component mixture with `s_real = 0.995`, `s_artifact = 0.5` | **42.3%** |
| read-level forensics | fraction of mRNA-called fragments whose best surviving junction has a short anchor, no annotated junction, or a mismatch within 5 bp | **43.4%** |

A third, non-circular estimator (novel/annotated junction mass ratio, tightly
conserved at 0.0143 ± 0.001 across clean controls, 0.308 here — 21.6×) puts 95.4%
of *novel*-junction read mass in this library in the artifact category.

**The problem is real, it is quantified, and its causes are mechanical rather than
statistical.** No single component is at fault; the failure is distributed across
`alignable` (blacklist construction), STAR (permissive stitching), and `rigel`
(blacklist matching semantics).

## 2. The consequential reframe

Escaping spliced artifacts are **not** the dominant cause of mRNA over-estimation
in this library. Splicing evidence supports only ~0.19% RNA content, while rigel
reports `mrna_fraction + nrna_fraction = 3.61%` — a 19× gap. Most of the
over-estimate comes from *unspliced* fragments the EM assigns to exons.

But the escaping spliced artifacts corrupt the two priors that the EM uses to make
those unspliced assignments:

- `fragment_length.rna.n_observations = 10,643` — exactly `spliced_annotated`
- `strand_model.n_training_fragments = 11,135` — trained only on exonic spliced
  fragments (`strand_model.py:391-473`)

Both are estimated from a population that is ~43% gDNA. rigel's own diagnostics
already flag the contradiction: `exonic_all_specificity = 0.5555` (near-random)
with `contamination_gap = 0.2301`.

**So the fix matters for a second-order reason: it de-contaminates the priors, not
primarily the spliced counts.** Any evaluation that measures only spliced-fragment
counts will understate the benefit.

## 3. Evidence base

### 3.1 Reproduction

The original run deleted every BAM (all declared `temp()`). The library was
re-run in an isolated tree with references symlinked, retaining all intermediates:

```
/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_gdna_debug/
  refs -> /scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs   (symlink)
  sample_sheets/gdna_debug_LBX0077.xlsx
  logs/rerun_LBX0077.log
  runs/human/mctp_LBX0077_SI_44263_HKWNYDRX7/
    fastq/k2_prefilter_{1,2}.dedup.fastq.gz    exact STAR input
    bam/star.srt.rmdup.collate.bam             direct rigel input
    bam/star.srt.rmdup.bam(.bai)               coord-sorted, indexed
    bam/spliced_only.bam                       11 MB, all spliced records
    rigel/annotated.bam                        rigel per-record ZI/ZJ/ZF/ZS/ZB tags
```

Command used (12 jobs; no reference rebuild, no Globus upload, no
salmon/kallisto/arriba/ciri3/bcftools/stringtie/zna):

```bash
export PATH=/home/mkiyer/sw/miniforge3/envs/hulkrna/bin:$PATH
export HULKRNA_PREFLIGHT_OK=1
BASE=/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_gdna_debug
snakemake --snakefile workflow/Snakefile --directory "$BASE" \
  --profile profiles/arc_armis2_slurm \
  --conda-prefix /scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/.snakemake/conda \
  --config sample_sheet="$BASE/sample_sheets/gdna_debug_LBX0077.xlsx" \
  --notemp \
  runs/human/mctp_LBX0077_SI_44263_HKWNYDRX7/rigel/annotated.bam
```

**Reproducibility is exact where it matters.** Every spliced counter is identical
to the original run: `spliced_annotated` 10,643, `spliced_unannotated` 2,126,
`spliced_implicit` 707, `with_annotated_sj` 12,492, `with_unannotated_sj` 5,715,
`n_training_fragments` 11,135. Drift is +213 alignment records / +2 fragments
(`fastq-dupaway --compare-seq loose` is memory-bounded and order-dependent);
`sj_blacklisted` 50,907 vs 50,912; `strand_specificity` 0.7859 vs 0.785631.

### 3.2 Evidence bundle

`/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_gdna_debug/evidence/`

| file | contents |
|---|---|
| `spliced_only.bam` | 126,813 spliced alignment records / 52,545 fragments, with rigel tags |
| `spliced_junctions.parquet` | 135,186 junction rows: rigel-semantics anchors, intron length, NM/AS/NH/MD, mismatch distance to junction, all Z tags |
| `spliced_annotated_joined.parquet` | above, joined to `splice_blacklist.feather` + `sj.feather`, with the reject rule re-applied |
| `bug1_halfblind_escapes.tsv` | 8,617 rows / 2,028 junctions / 4,317 fragments |
| `bug2_bothsites_known_but_pair_missing.tsv` | 21,617 rows / 10,479 junctions / 2,851 fragments |
| `bug3_anchor_ceiling_62.tsv` | 538 rows / 307 junctions / 249 fragments |
| `bug4_sjdb_stitch_no_entry.tsv` | 7,706 rows / 1,398 junctions / 5,163 fragments |
| `validation_targets.json` | pre-registered expected counts per bug class |
| `extract_spliced.py`, `analyze.py`, `discriminate.py`, `sitesets.py`, `make_targets.py` | reproducible analysis scripts |
| `rigel_summary_{original,rerun}.json` | reproducibility comparison |

> **Note:** this lives on `/scratch`, which is subject to purge policy. Copy the
> bundle somewhere durable before relying on it long-term.

### 3.3 Validation of the analysis itself

The analysis re-implements rigel's reject rule from source and **reproduces
`sj_blacklisted = 50,907` exactly**. This confirms two semantics that are easy to
get wrong and that any future work must respect:

1. **Anchors are cumulative, not per-block.** From `bam_scanner.cpp:430-470`:
   `anchor_left` = *all* ref-advancing bases (M/D/=/X, excluding N) before the
   junction; `anchor_right = total - anchor_left`. For multi-junction reads the
   anchors span earlier introns, making the reject test systematically **more
   lenient** on multi-junction reads.
2. **Coordinates are 0-based half-open** `[start, end)` in both
   `splice_blacklist.feather` and `sj.feather`. STAR `SJ.out.tab` (1-based
   inclusive) maps as `bl_start = SJ_start - 1`, `bl_end = SJ_end`. Verified two
   ways: offset scan (6,080 hits at this offset, 0–2 at all others) and 100.0000%
   agreement between `sj.feather` membership and STAR's `annotated_flag` over
   5,741 junctions.

## 4. Where the alignments actually go

Of 135,186 junction rows in the index case:

```
rejected by blacklist                        50,907  (37.7%)
survived, in blacklist but escaped           14,554  (10.8%)
survived, not in blacklist at all            69,725  (51.6%)
                                            -------
surviving junction rows                      84,279
```

Escape routes among the 14,554 that *were* in the blacklist:

| route | rows | share |
|---|---|---|
| half-blind entry (one stored anchor `= 0`) | 8,617 | **59.2%** |
| min-side margin ≤ 5 | 2,340 | 16.1% |
| min-side margin ≤ 10 | 4,543 | 31.2% |
| stored anchor saturated at the 62 ceiling | 538 | 3.7% |

**The dominant failure is coverage** (51.6% of surviving rows have no blacklist
entry at all), followed by the half-blind-entry bookkeeping flaw. The 62-bp anchor
ceiling — which looked like the primary cause from aggregate statistics — accounts
for only 538 rows at the read level.

## 5. Two distinct artifact populations

They have **orthogonal** signatures and require **different** filters. Any fix
that targets only one will leave the other untouched.

| class | n rows | anchor ≤10 | mismatch ≤5bp | NM≥3 | multimapping |
|---|---|---|---|---|---|
| `not_in_bl_annot` | 21,608 | **35.7%** | 1.6% | 0.9% | 13.0% |
| `not_in_bl_novel` | 48,117 | 5.8% | **20.1%** | **26.1%** | **82.8%** |
| `bl_escape` | 14,554 | 4.7% | 8.7% | 16.2% | 52.5% |

**Class A — sjdb-stitch.** A gDNA fragment matches ≥3 bases by chance past an
annotated junction and STAR stitches it (`alignSJDBoverhangMin = 3`). The
sequence is *pristine*: the `anchor_min ≤ 3` bin has mean NM **0.14** and only
**0.24%** carry a mismatch within 5 bp of the junction. **Mismatch-based
detection cannot work on this class** — three bases match perfectly far too
often. Anchor length is the only discriminator.

**Class B — novel multimapper.** Long anchors, high mismatch load, 20% with a
mismatch within 5 bp, and **82.8% multimapping**. Anchor filters do not touch it;
mismatch and multimap filters do.

Anchor-length distribution by class confirms the split:

```
class              n        5%    25%   50%   75%   95%
rejected        50,907      3      3     3    12    47
not_in_bl_annot 21,608      3      3    23    47    69     <- bimodal
not_in_bl_novel 48,117     10     24    40    57    70
bl_escape       14,554     11     35    50    64    72
```

A revealing detail: in the `anchor_min` 4–5 bin, **14.0%** carry a mismatch
within 5 bp, versus 0.24% in the 1–3 bin. This is `alignSJstitchMismatchNmax
5 -1 5 5` permitting mismatched stitches once the anchor is long enough to
require them.

## 6. Cohort context

The index case is an outlier even among high-gDNA libraries. Across 344 libraries
with both `rigel/summary.json` and `star/SJ.out.tab.gz`:

- novel-junction blacklist capture correlates with `gdna_fraction` at Pearson
  **r = +0.857** / Spearman **+0.890** — the most gDNA-responsive junction-level
  statistic available
- `LBX0077` ranks **1 of 344** (worst) on strand specificity, at 0.7856
- the best on-disk predictor of leakage among the 19 libraries with
  `gdna_fraction ≥ 0.90` is `calibration.gdna_density_global`
  (r = +0.758), **not** `gdna_fraction` (r = +0.347)
- two same-flowcell siblings at 74% and 82% gDNA show clean strand specificity
  (0.9911, 0.9931). They have **104 bp** mates; `LBX0077` has **138 bp**.
  `104/2 = 52 < 62`, `138/2 = 69 > 62` — consistent with the anchor ceiling being
  a contributing mechanism, and implying **the problem worsens as read lengths
  grow**.

Second-worst candidate for follow-up: `mctp_MI_1104_SI_43016_HNMVYDMXY`
(98.6% gDNA, strand specificity 0.9088, 19.1M fragments).

## 7. Design principle

> **The system must be robust to aligner settings.**

STAR parameter tightening is a mitigation, not a fix: it is fragile (any future
parameter change, STAR version bump, or new library chemistry re-opens the hole),
it is shared with arriba fusion calling, and it trades true positives directly.
The durable defenses live in `alignable` (what artifacts are catalogued) and
`rigel` (how the catalogue is matched). STAR changes should be limited to
choices that are defensible on their own merits.

A corollary that the evidence supports: **the blacklist should be built and
matched so that it does not depend on having simulated the exact junction.** The
splice-site decomposition in
[GDNA_SPLICE_ARTIFACTS_RIGEL.md](GDNA_SPLICE_ARTIFACTS_RIGEL.md) §3 is the
cheapest realisation of that principle.

## 8. Recommended priority order

Detail and specifications are in the companion documents. Ranked by
(benefit × confidence) / cost:

| # | change | owner | benefit | cost | needs zarr regen? |
|---|---|---|---|---|---|
| 1 | splice-**site** blacklist: reject when left site AND right site both known | rigel | novel capture 36.4% → **64.6%**; cost 0.19–0.59% of genuine junctions | low | no |
| 2 | fix half-blind aggregation (one anchor stored as 0) | alignable | 8,617 rows / 4,317 fragments become rejectable | low | rebuild from existing zarr |
| 3 | expose and lower `--splice-blacklist-min-count` | rigel/hulkrna | annotated capture 53.6% → 70.5%, novel 36.4% → 49.2% | trivial | no |
| 4 | raise `--read-lengths` past production read length | alignable | removes the 62-bp ceiling; future-proofs against longer reads | high (re-simulate) | yes |
| 5 | raise `--coverage` above 5 | alignable | more junction *pairings*; addresses the 51.6% no-entry population | high (re-simulate) | yes |
| 6 | `alignSJDBoverhangMin` 3 → 8–10 | hulkrna/STAR | suppresses class A at the source | low | ideally yes, for consistency |
| 7 | multimap-aware spliced-evidence weighting | rigel | targets class B (82.8% multimapping) | medium | no |

Items 1–3 are cheap, do not require re-simulating the 26 GB zarr, and are
independently testable against the pre-registered targets in
`validation_targets.json`.
