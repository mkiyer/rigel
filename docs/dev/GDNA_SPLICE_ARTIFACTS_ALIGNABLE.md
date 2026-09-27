# `alignable`: Defects and Specifications for Splice-Artifact Cataloguing

**Audience:** a dedicated session working in `/home/mkiyer/proj/alignable`
**Status:** specification; no changes made to `alignable` by this investigation
**Parent:** [GDNA_SPLICE_ARTIFACTS_OVERVIEW.md](GDNA_SPLICE_ARTIFACTS_OVERVIEW.md)
**Version investigated:** `alignable` 0.1.0, zarr built 2026-04-18

---

## 0. What this document is

`alignable` produces the splice-junction artifact catalogue that `rigel` uses to
reject false-positive spliced alignments of genomic DNA. This document specifies
**three defects and two coverage limitations** in that catalogue, each with the
exact code location, the mechanism, the measured impact on a real library, and a
pre-registered validation target.

**Two hypotheses were tested and REJECTED.** Do not spend time on them:

1. ~~`alignable` aligns with different STAR parameters than production~~ —
   **false.** The zarr's `aligner_config` attribute is byte-identical to
   `hulkrna/resources/star_arriba_parameters.txt`; `diff` returns identical,
   including `alignSJstitchMismatchNmax 5 -1 5 5`.
2. ~~`alignable`'s STAR index lacked the gencode sjdb, so annotated-junction
   stitching was never simulated~~ — **false.** All 56,166 blacklist entries with
   `max_anchor` of 3–4 bp are **100.000%** annotated junctions (69.75×
   enrichment over the 1.434% background). A 3-bp anchor is reachable *only* via
   `alignSJDBoverhangMin`, so the annotation was loaded and the stitch mechanism
   *was* simulated.

## 1. The production artefact under discussion

```
/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/alignable.zarr.zip   26 GB
```

Recorded root attributes (ground truth for what the current blacklist covers):

```json
{
  "aligner": "star",
  "aligner_version": "2.7.11b",
  "aligner_config": "<byte-identical to hulkrna star_arriba_parameters.txt>",
  "read_length_bins": [50, 75, 100, 125],
  "frag_len_mean": 250, "frag_len_sd": 50,
  "frag_len_min": 50,   "frag_len_max": 1000,
  "error_rate": 0.01,
  "tolerance": 3,
  "coverage": 5,
  "genome_fasta": "/scratch/.../shared_data/alignable/refs/genome_controls.fasta.bgz"
}
```

The generating config survives at
`/scratch/mkiyer_root/mkiyer0/shared_data/alignable/config.yaml`. Note the STAR
index it referenced
(`/scratch/.../shared_data/alignable/refs/star_index`) **has been deleted**;
regeneration will need a fresh index built with `--sjdbGTFfile` and
`--sjdbOverhang` matching production (currently 149).

Raw splice table inside the zip: **34,894,972** `(chrom, intron, read_length)`
rows spanning **27,018,296** unique junctions. After `rigel index` applies
`--splice-blacklist-min-count 2`, the shipped blacklist is **5,833,092** rows.

## 2. DEFECT 1 — half-blind aggregation stores an unusable anchor of `0`

**Severity: high. Cheapest high-value fix. Does not require re-simulation.**

### Location

`src/alignable/cli.py`, `_aggregate_splice_artifacts` (defined at line 920):

```python
960:  left_mask  = pc.less_equal(table["anchor_left"],  table["anchor_right"])
961:  right_mask = pc.less_equal(table["anchor_right"], table["anchor_left"])
964:      table.filter(left_mask)      # -> max_anchor_left  computed ONLY here
970:      table.filter(right_mask)     # -> max_anchor_right computed ONLY here
999:      pc.fill_null(out["max_anchor_left"],  0),
1004:     pc.fill_null(out["max_anchor_right"], 0),
```

Note `left_mask` and `right_mask` both use `<=`, so ties (`anchor_left ==
anchor_right`) fall into both — the sides are not partitioned, they overlap.
That is not the defect, but it is worth understanding before rewriting.

### Mechanism

`max_anchor_left` is updated **only** from records where
`anchor_left <= anchor_right`, and symmetrically for the right side. The
aggregation then left-outer-joins and fills missing values with `0`.

A junction never observed with its *left* side limiting therefore ships with
`max_anchor_left = 0`.

`rigel` rejects a junction when
(`bam_scanner.cpp:493-497`):

```cpp
return sj.anchor_left  <= hit->max_anchor_left
    || sj.anchor_right <= hit->max_anchor_right;
```

A real CIGAR junction always has both anchors ≥ 1, so a stored `0` **can never
satisfy the comparison**. The entry is inert on that side.

### Measured impact

- **3,530,708 of 5,833,092** blacklist rows (**60.53%**) have exactly one anchor
  equal to `0`. (Zero rows have both at `0`.)
- Of the 6,080 blacklist entries actually hit by the index-case library's
  junctions, **3,631 (59.7%)** are half-blind.
- **This is the single largest escape route** among junctions that *are* in the
  blacklist: 8,617 of 14,554 escaping junction rows (**59.2%**).

### Specification

Store the true per-side maximum, computed independently for each side over **all**
records touching that junction — not only the records where that side happened to
be the shorter one.

```
max_anchor_left[j]  = max over all observed artifact records r at junction j of r.anchor_left
max_anchor_right[j] = max over all observed artifact records r at junction j of r.anchor_right
```

Then a stored value is `0` only if the junction was never observed at all, in
which case no row should be emitted.

**Semantic note for the fix session:** think carefully about which envelope you
intend. Recording the per-side max over all records makes each side's threshold
*larger*, so the OR-based reject rule becomes **more aggressive**. That is
desirable for sensitivity but changes the specificity profile, so it must be
evaluated jointly on a clean library (see §7). An alternative, more conservative
reading is to record `max(min(anchor_left, anchor_right))` as a single scalar
"limiting anchor" per junction and have `rigel` compare against
`min(observed_left, observed_right)`. **Coordinate this choice with the `rigel`
side** — see [GDNA_SPLICE_ARTIFACTS_RIGEL.md](GDNA_SPLICE_ARTIFACTS_RIGEL.md) §2.

### Validation target

`evidence/bug1_halfblind_escapes.tsv` — **8,617 junction rows / 2,028 distinct
junctions / 4,317 fragments** in `spliced_only.bam` that currently escape solely
because one stored anchor is `0`. After the fix, re-running `rigel quant` on the
retained BAM must reject a large majority of these. Each row carries the observed
`anchor_left`/`anchor_right` and the stored `max_anchor_left`/`max_anchor_right`,
so the expected outcome is checkable without re-running anything.

## 3. DEFECT 2 — the anchor ceiling is `floor(max_read_length / 2)`

**Severity: medium at current read lengths, growing. Requires re-simulation.**

### Mechanism

Because aggregation records the **limiting** (shorter) anchor, the largest anchor
observable in read-length bin `L` is exactly `floor(L / 2)`. With
`read_length_bins = [50, 75, 100, 125]` the per-bin ceilings are **25, 37, 50,
62**, and the observed maximum across all 5,833,092 rows is exactly **62**.

The distribution is saturated against that ceiling:

```
max(max_anchor_left, max_anchor_right):  min=3  p50=45  p90=60  p99=62  max=62
```

Any junction whose crossing reads have **both** anchors > 62 is unrejectable
regardless of blacklist membership.

### Measured impact

- The index-case library has **138 bp** mates (151 bp raw; STAR
  `Average input read length` 278 after PE merge) and STAR reports a maximum
  spliced overhang of **75**.
- 2,398 of 13,988 junctions (17.1%), carrying **36.5%** of unique crossing-read
  mass, have `max_spliced_overhang > 62`.
- At the read level, however, only **538 junction rows** escape *specifically*
  at the 62 ceiling. The ceiling is a real but **secondary** route — the
  aggregate statistic overstates it because most >62 junctions have no blacklist
  entry at all (defect 3 / limitation).
- **Trend evidence:** two same-flowcell siblings at 74% and 82% gDNA have 104 bp
  mates (`104/2 = 52 < 62`) and show clean strand specificity (0.9911, 0.9931).
  The index case has 138 bp mates (`138/2 = 69 > 62`) and leaks. **The ceiling
  becomes the binding constraint as read length grows.**

### Specification

- `--read-lengths` must cover the production read length. Current production:
  `fastp.max_len = 150`, mates ~138 bp post-trim, merged reads longer still.
  Simulate at least up to **150**, and prefer a bin set that brackets future
  chemistry (e.g. `50,75,100,125,150,175,200`).
- Record in the zarr attributes the **maximum representable anchor**, so `rigel`
  can detect at match time that an observed anchor exceeds the catalogue's
  dynamic range and decline to claim the junction is clean. This is the
  robustness property that matters: **the catalogue should be able to say
  "unknown" rather than silently implying "not an artifact."**

### Validation target

`evidence/bug3_anchor_ceiling_62.tsv` — 538 rows / 307 junctions / 249 fragments.
A necessary (not sufficient) check: after regeneration at
`--read-lengths` ≥ 150, `max(max_anchor_*)` in the emitted blacklist must exceed
62, and these junctions must become rejectable.

## 4. DEFECT 3 — `coverage: 5` under-samples junction *pairings*

**Severity: high — this is the dominant failure. Requires re-simulation.**

### Mechanism

At `coverage: 5`, each genomic position contributes 5 fragments per read-length
bin. For any given spurious junction, a **singleton** observation is the modal
outcome: of 34,894,972 raw rows, **26,821,287 (76.9%)** have `count == 1`.
`rigel index` then applies `--splice-blacklist-min-count 2` (default), discarding
**21.2 M unique junctions — 78.4% of the catalogue**.

The result is that most artifact junctions were either never sampled or sampled
once and dropped.

### Measured impact

- **51.6%** of surviving junction rows in the index case (69,725 of 135,186) have
  **no blacklist entry at all**. This is larger than every other escape route
  combined.
- Realised blacklist sensitivity on this library's novel junctions is **36.4%**
  (3,000 of 8,247). At 96.4% gDNA essentially all novel junctions are DNA-derived,
  so this is a direct sensitivity measurement.
- Only **20.69%** of gencode v46's 404,168 annotated junctions have any blacklist
  entry — so **79.31%** are unprotected against a 3-bp `alignSJDBoverhangMin`
  stitch.
- Relaxing to `min_count = 1` raises novel capture 36.4% → 49.2% and annotated
  capture 53.6% → 70.5%. **Still leaves 50.8% of novel junctions uncovered** — so
  the threshold is not the whole story; sampling depth is.

### The near-miss result — why more depth pays off disproportionately

Of the index case's 5,247 non-blacklisted novel junctions:

| | index case | clean controls | 74–82% gDNA |
|---|---|---|---|
| shares donor **or** acceptor with a blacklist entry | **71.4%** | 20.7–24.3% | 54.7–60.3% |
| shares **both** sites (with different partners) | **44.4%** | 2.8–4.9% | 25.3–27.0% |

**The simulation already knows both splice sites — it never observed them
paired.** This is the signature of insufficient sampling of *combinations*, not
of ignorance about which loci are dangerous.

It also means there is a much cheaper route to most of the benefit: match on
sites rather than exact pairs. See
[GDNA_SPLICE_ARTIFACTS_RIGEL.md](GDNA_SPLICE_ARTIFACTS_RIGEL.md) §3 — measured at
**36.4% → 64.6%** capture with a 0.19–0.59% specificity cost and **no
re-simulation**. Regenerating at higher coverage should be evaluated *after* that
change, since it may become redundant.

### Specification

- Raise `--coverage` (5 → 20+ recommended for evaluation) and re-measure. Cost is
  approximately linear in coverage × number of read-length bins; the current zarr
  took >24 h on 36 cores, so a 4× coverage and 7 bins (vs 4) is roughly a 7×
  increase. Budget accordingly and consider whether a splice-only product is
  worth building (see §6).
- Record the achieved per-junction observation count distribution in the zarr
  attributes so downstream consumers can reason about statistical support instead
  of guessing.
- Emit `count` per junction in the shipped blacklist (it exists in the raw table)
  so `rigel` can threshold at match time rather than at index time. This turns an
  irreversible index-time decision into a tunable runtime parameter.

### Validation target

`evidence/bug2_bothsites_known_but_pair_missing.tsv` — **21,617 junction rows /
10,479 distinct junctions / 2,851 fragments** where both sites are already in the
catalogue but the pair is absent. 21,549 of the 21,617 rows are novel junctions.
Higher coverage should convert a substantial share of these into real entries.

## 5. LIMITATION — fragment-length distribution

**Severity: low-to-medium. Frequently misdiagnosed — read this before acting.**

Fragment lengths are **already varied**: `frag_len_mean = 250`, `frag_len_sd =
50`, clamped to `[50, 1000]`. This is *not* a single fixed length.

However it does not match production. The index-case library:

```
global   fragment length: mean 165.76  sd 61.14  median 155  mode 161
gdna     fragment length: mean 189.45  sd 77.95  median 165  mode 162
```

The simulation is centred **85 bp longer** than the real gDNA population
(250 vs 165). Consequences:

- Blacklist intron lengths skew long by construction (median **194,310**;
  1st percentile 87; max 4,819,393), because 250 bp fragments spanning the genome
  produce long-range stitches.
- The index case's escaping novel junctions have **longer** introns (median
  5,433) than its blacklisted ones (median 2,037) — i.e. long-range spurious
  stitches are *under*-sampled, consistent with the fragment model being centred
  too long **and** too narrow to populate the tails densely.
- By contrast, in `LBX0054` (74% gDNA, 104 bp reads) the artifact classes are much
  shorter (novel median 603, blacklisted 234), so the short-intron signature is
  **read-length dependent** and is not a general discriminator.

### Specification

Match the empirical distribution rather than guessing. `rigel` already writes
`fragment_lengths.feather` per library and
`summary.json → fragment_length.gdna` (`mean`, `std`, `median`, `mode`). Either:

- parameterise from a representative high-gDNA library
  (`frag_len_mean ≈ 165`, `frag_len_sd ≈ 78`), or better
- accept an **empirical fragment-length histogram** file so the simulation can
  reproduce the real distribution including its tails, and record the source in
  the zarr attributes.

Widening `frag_len_sd` matters as much as re-centring `frag_len_mean`: the
escaping junctions live in the tails.

## 6. Cost control — is a splice-only product feasible?

The full zarr is a **mappability** product; the splice blacklist is a by-product.
Regeneration at higher coverage and more read lengths is expensive (>24 h × ~7).

Worth investigating in the fix session:

- `scripts/reaggregate_splice.py`, `scripts/migrate_splice_to_feather.py`,
  `scripts/reaggregate_splice.py` — whether re-aggregation can be done from
  retained intermediates without re-alignment.
- Whether a CLI mode can emit **only** the splice table (skipping per-base
  mappability arrays), which would cut both runtime and the 26 GB footprint
  dramatically and make iteration on coverage/read-length practical.
- **DEFECT 1 does not need re-simulation at all** — it is a re-aggregation of the
  existing raw table inside the zip. Do that first and ship it.

## 7. Joint sensitivity/specificity evaluation — required for every change

Every change must be evaluated on **both** an extreme-gDNA library and a clean
library. Loosening the catalogue has a real specificity cost:

Moving to `min_count = 1` would newly blacklist **3,315 annotated junctions**
carrying **185,542 unique crossing reads** in clean control `LBX0629`
(annotated capture 25.1% → 43.6%). Under the current anchor rule most still
escape via long anchors, so immediate damage is limited — but raising coverage
raises `max_anchor` values, which starts rejecting genuine RNA.

### Recommended evaluation panel

| role | library | gDNA | note |
|---|---|---|---|
| positive (artifact-rich) | `mctp_LBX0077_SI_44263_HKWNYDRX7` | 0.964 | index case, BAM retained |
| positive (second) | `mctp_MI_1104_SI_43016_HNMVYDMXY` | 0.986 | strand specificity 0.9088 |
| negative (clean) | `mctp_MI_1089_SI_40085_HKF7JDRX5` | 0.001 | strand 0.9982 |
| negative (clean) | `mctp_MI_1083_SI_40080_HKF7JDRX5` | 0.003 | strand 0.9975 |
| intermediate | `mctp_LBX0054_SI_43886_HKWNYDRX7` | 0.742 | 104 bp mates, does **not** leak |

### Primary metrics

1. **Sensitivity** — fraction of the index case's novel junctions captured
   (currently 36.4%).
2. **Specificity cost** — count of newly-rejected junctions in the clean controls
   that are annotated with `max_spliced_overhang > 20` and `n_uniquely_mapped ≥ 5`
   (presumed genuine). Baselines: 1,441 and 1,087 such junctions respectively.
3. **Strand specificity of the index case** — `summary.json →
   strand_model.strand_specificity`. Currently **0.7856**; clean siblings are
   0.991–0.999. **This is the best single end-to-end readout**, because it is
   computed from the surviving spliced population and is independent of the
   blacklist itself. Moving it toward 0.99 is the goal.
4. **`splice_artifact / all_spliced` ratio** — currently 0.491 in the index case
   versus 0.010–0.046 in every other library.

### Reproduction commands

Re-running `rigel quant` against a new index using the retained BAM takes ~25 min
and requires no re-alignment:

```bash
BAM=/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_gdna_debug/runs/human/\
mctp_LBX0077_SI_44263_HKWNYDRX7/bam/star.srt.rmdup.collate.bam
rigel quant --index <NEW_INDEX> --bam "$BAM" --threads 8 -o <OUTDIR> \
            --annotated-bam <OUTDIR>/annotated.bam
```

Then re-run `evidence/extract_spliced.py` and `evidence/analyze.py` against the
new `annotated.bam` and compare to `evidence/validation_targets.json`.
