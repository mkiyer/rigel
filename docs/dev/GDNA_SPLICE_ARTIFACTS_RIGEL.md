# `rigel`: Improving Splice-Artifact Detection

**Audience:** a session working in `/home/mkiyer/proj/rigel`
**Status:** specification; no changes made to `rigel` by this investigation
**Parent:** [GDNA_SPLICE_ARTIFACTS_OVERVIEW.md](GDNA_SPLICE_ARTIFACTS_OVERVIEW.md)
**Version investigated:** rigel 0.7.1, index format_version 7

---

## 0. Why the highest-leverage fixes live here

`rigel` receives the BAM and decides which spliced alignments are real. It is the
right place for the durable defence, because:

- it can be changed and re-evaluated in **~25 minutes** against a retained BAM,
  with no re-alignment and no 26 GB re-simulation;
- it sees per-fragment evidence (anchors, mismatches, multimap structure, strand,
  the gDNA density model) that neither STAR nor `alignable` can combine;
- it is the component that makes the system **robust to aligner settings**, which
  is the stated design goal.

The single best-measured intervention in this whole investigation is a `rigel`
matching change requiring no new data (§3).

## 1. Current mechanism (as implemented in 0.7.1)

### 1.1 Anchor computation — `bam_scanner.cpp:430-470`

```cpp
if (op == CIG_REF_SKIP) {
    sj.start = pos; sj.end = pos + len;
    sj.anchor_left = ref_advanced;      // cumulative, NOT the adjacent block
    ...
}
// after the CIGAR loop:
sjs[k].anchor_right = ref_advanced - sjs[k].anchor_left;
```

`anchor_left` is the **cumulative** count of ref-advancing bases (M/D/=/X,
excluding N) from the alignment start to the junction; `anchor_right` is the
remainder. I/S/H/P do not contribute.

**Implication:** for multi-junction reads the anchors span earlier introns, so
both anchors are larger than the adjacent exon blocks and the reject test is
**systematically more lenient on multi-junction reads.** This is a latent
soundness issue independent of everything else — see §6.

### 1.2 Reject rule — `bam_scanner.cpp:493-497`

```cpp
const auto* hit = resolver.sj_blacklist_lookup(ref_id, sj.start, sj.end);
if (!hit) return false;
return sj.anchor_left  <= hit->max_anchor_left
    || sj.anchor_right <= hit->max_anchor_right;
```

Two properties worth stating explicitly:

1. **Exact-triple match required.** `(ref, start, end)` must be present. A
   junction one base away, or a novel pairing of two known artifact sites, is
   invisible.
2. **A stored anchor of `0` is inert.** A real CIGAR junction always has both
   anchors ≥ 1, so `anchor <= 0` never fires. 60.53% of blacklist rows have
   exactly one anchor at `0`
   (see [GDNA_SPLICE_ARTIFACTS_ALIGNABLE.md](GDNA_SPLICE_ARTIFACTS_ALIGNABLE.md) §2).

Coordinates are **0-based half-open** `[start, end)`, matching
`bam_scanner.cpp:449-450` (`sj.start = pos` with `pos = b->core.pos`).

### 1.3 Classification precedence — `resolve_context.h:1311-1320, 1379-1384`

A fragment is relabelled `SPLICE_ARTIFACT` only when `n_sj_blacklisted > 0`
**AND** its classification is not already `SPLICE_SPLICED_ANNOT` /
`SPLICE_SPLICED_UNANNOT`. `SPLICED_IMPLICIT` detection additionally requires
`splice_type == SPLICE_UNSPLICED`.

**Consequence:** a fragment with two junctions — one blacklisted, one surviving —
is counted `spliced_annotated`, enters the strand model, and enters the RNA
fragment-length model. Measured: **3,423 junction rows** were individually
rejected while their record still reads `spliced_annot`/`spliced_unannot`.

### 1.4 Counter semantics — `bam_scanner.cpp:1305-1312`

`stats_.n_sj_blacklisted` is incremented **per alignment record** inside
`parse_bam_record`, including secondary and supplementary alignments. The index
case's BAM has 39,189,259 records for 17,821,747 read-name groups (2,530,267
secondary + 1,033,751 supplementary).

**Do not compare `sj_blacklisted` (50,907, per-record) against
`spliced_annotated` (10,643, per-fragment).** They are different units. The
per-fragment picture: 52,545 spliced fragments, of which 135,186 junction rows →
50,907 rejected, 84,279 surviving.

## 2. FIX 1 — make the anchor semantics coherent and defensible

The current `max_anchor_left` / `max_anchor_right` pair combined with an `OR`
reject rule is hard to reason about, and interacts badly with `alignable`'s
half-blind aggregation.

### Options

**(a) Single limiting-anchor scalar (recommended).** Have `alignable` store one
value per junction — the maximum observed **limiting** anchor,
`max over records of min(anchor_left, anchor_right)` — and have `rigel` reject
when `min(observed_left, observed_right) <= max_limiting_anchor`. This is a single
well-defined quantity, is exactly what the current aggregation *tried* to compute,
and eliminates the `0`-inert failure mode by construction.

**(b) True per-side maxima.** Keep two values but compute each over all records
(the `alignable` §2 fix). More aggressive than (a) because each side's threshold
rises independently. Requires joint specificity evaluation.

**Either way, `rigel` must treat a stored `0` as "no information" rather than as a
threshold**, and must not silently accept a junction because one side was never
observed.

### Dynamic-range guard

`alignable` should record the maximum representable anchor (currently **62** =
`floor(125/2)`). When an observed anchor exceeds it, `rigel` cannot conclude the
junction is clean — the catalogue simply has no reach there. Expose this as a
distinct outcome (e.g. `sj_blacklist_out_of_range`) rather than folding it into
"passed", and surface the count in `summary.json`. **This is the single most
important robustness change**: it converts a silent false negative into a
reportable one.

## 3. FIX 2 — splice-**site** blacklist (highest-value change, no new data)

This is the proposal to build left-side and right-side site sets and reject a
junction when **both** of its sites are known artifact sites. It was measured and
it is close to free.

### Rationale

Of the index case's 5,247 non-blacklisted novel junctions, **71.4%** share a donor
or acceptor with a blacklist entry and **44.4%** share **both** (with different
partners); clean controls sit at 20.7–24.3% and 2.8–4.9%. The simulation already
knows which loci are dangerous — it never observed the specific pairing.

### Measured performance

Site sets derived from the existing 5,833,092-row blacklist: **1,172,288** distinct
left sites, **1,150,976** distinct right sites.

Junction capture, exact-triple rule vs both-sites rule:

| library | gDNA | exact | both-sites | newly-caught *genuine* junctions¹ |
|---|---|---|---|---|
| `mctp_LBX0077_SI_44263_HKWNYDRX7` | 0.964 | 43.5% | **60.2%** | 1 of 523 (0.19%) |
| `mctp_MI_1104_SI_43016_HNMVYDMXY` | 0.986 | 37.6% | **64.1%** | 6 of 1,024 (0.59%) |
| `mctp_MI_1089_SI_40085_HKF7JDRX5` | 0.001 | 20.6% | 21.1% | 4 of 1,441 (0.28%) |
| `mctp_MI_1083_SI_40080_HKF7JDRX5` | 0.003 | 19.6% | 19.9% | 2 of 1,087 (0.18%) |

¹ newly-rejected junctions that are annotated with `max_spliced_overhang > 20`
and `n_uniquely_mapped_crossing ≥ 5` — the presumed-genuine population.

**Novel-junction capture in the index case: 36.38% → 64.62%**, adding 2,329
junctions carrying 478 unique reads.

**The rule is self-targeting.** It adds +16.7 and +26.5 percentage points of
capture in the two extreme-gDNA libraries and only +0.3 to +0.5 points in clean
libraries — because clean libraries' real junctions rarely have *both* sites in
the artifact catalogue. Ceiling check: across all of gencode v46, the rule moves
annotation exposure from 20.69% to only **21.50%** (+0.81 points), so it cannot
become catastrophic even in the worst case.

### Specification

At index time, additionally emit from the alignable splice table:

```
sj_blacklist_left_sites  : (ref, start) -> max_anchor_left_observed
sj_blacklist_right_sites : (ref, end)   -> max_anchor_right_observed
```

At match time, a junction is a blacklist candidate if **either**:

- the exact triple `(ref, start, end)` is present (current behaviour), **or**
- `(ref, start)` is in the left-site set **AND** `(ref, end)` is in the
  right-site set.

Then apply the anchor test as usual, using the per-site anchor maxima for the
site-derived path.

### Cautions

- **Anchor test must still apply.** Do not reject on site membership alone; the
  anchor test is what preserves specificity. The measured 0.19–0.59% cost above
  is *with* the anchor test in force.
- **Report the two paths separately** in `summary.json` (e.g.
  `sj_blacklisted_exact` vs `sj_blacklisted_site_pair`) so the site rule's
  contribution and cost are auditable per library.
- Consider gating the site rule on estimated gDNA content. It is nearly free in
  clean libraries, so gating is probably unnecessary — but making it a documented
  toggle allows a clean A/B.
- **Do not** relax to "either site is in the blacklist." At 71.4% one-site
  sharing in the index case but 20.7–24.3% in clean controls, the specificity cost
  would be an order of magnitude larger and was not measured here.

### Also worth testing: should the blacklist include ANY junction?

The related question — drop `--splice-blacklist-min-count` to 1 — was measured
separately: annotated capture 53.6% → 70.5%, novel 36.4% → 49.2%. Cost: 3,315
newly-blacklisted annotated junctions carrying 185,542 unique crossing reads in
clean control `LBX0629`.

**The site-pair rule dominates it** (+28 points of novel capture vs +13) at a much
lower specificity cost, because requiring *both* sites is a stronger condition
than the min-count relaxation. Evaluate the site rule first; `min_count` may become
unnecessary. If both are adopted, evaluate them jointly rather than additively.

**Recommendation:** expose `--splice-blacklist-min-count` through
`hulkrna` `config.yaml → rigel.index_options` (currently only
`--nrna-tolerance 20` is passed) so it is tunable without editing rigel, and
better, **ship `count` per junction in the blacklist** so the threshold becomes a
`rigel quant` runtime parameter rather than an irreversible index-time decision.

## 4. FIX 3 — `--splicing-anchor-tolerance` is not the knob you want

Currently `3` (`summary.json → scan.splicing_anchor_tolerance`). This is a
**separate** mechanism from the blacklist test and does **not** relax or tighten
it. Raising it will not address the artifacts described here.

It is worth documenting the distinction in the rigel manual, because the name
invites exactly this confusion.

If a *global* anchor floor is wanted — reject any spliced fragment whose
`min(anchor_left, anchor_right)` is below a threshold, regardless of blacklist
membership — that should be a **new, explicitly-named parameter**, e.g.
`--min-splice-anchor`. Measured basis for such a floor:

| threshold | index-case surviving annotated rows removed | signature |
|---|---|---|
| `< 4` | 5,472 | mean NM 0.14 — pristine chance matches |
| `< 6` | 6,319 | + the 14% near-junction-mismatch bin |
| `< 11` | 7,706 | 35.7% of the class |

This is the `rigel`-side equivalent of STAR's `alignSJDBoverhangMin` and is
**strictly preferable to changing STAR**, because it is aligner-independent,
tunable per run, and does not perturb arriba. It directly serves the robustness
goal. **The true-positive cost has not been measured** and must be, on the clean
controls, before adopting a value.

## 5. FIX 4 — use the evidence rigel already has but ignores

Two artifact populations with orthogonal signatures (parent document §5). The
class-B population is currently unaddressed:

| class | n rows | anchor ≤10 | mismatch ≤5bp | NM≥3 | multimapping |
|---|---|---|---|---|---|
| `not_in_bl_annot` (class A) | 21,608 | 35.7% | 1.6% | 0.9% | 13.0% |
| `not_in_bl_novel` (class B) | 48,117 | 5.8% | **20.1%** | **26.1%** | **82.8%** |

### 5a. Multimap-aware spliced evidence

**82.8%** of surviving novel-junction rows are multimapping, and 78.7% of
surviving novel junctions have **zero** uniquely-mapped crossing reads (clean
controls: 16.3–19.2%). `rigel` runs with `include_multimap = true` and uses these
placements.

A spliced fragment whose *only* support is multimapping placements is weak
evidence for a novel junction. Options, in increasing sophistication:

1. exclude multimapper-only novel junctions from the **strand model** and the
   **RNA fragment-length model** training sets (cheapest; directly de-contaminates
   the two corrupted priors);
2. down-weight their contribution by posterior `ZW`;
3. require ≥1 unique placement before a novel junction can promote a fragment to
   `spliced_unannotated`.

Option 1 alone should move `strand_specificity` materially and is low-risk,
because it changes only what the model *trains on*, not what gets counted.

### 5b. Mismatch proximity to the junction

`MD` is present in the BAM (`outSAMattributes ... NM MD MC`) and is currently
unused for splice adjudication. Measured discrimination on surviving rows:

```
class              mean NM   mismatch within 5bp   within 10bp
not_in_bl_annot       0.21          1.62%             2.20%
not_in_bl_novel       2.01         20.08%            29.02%
bl_escape             1.35          8.75%            14.12%
rejected              0.37          2.32%             3.65%
```

A mismatch within 5 bp of a junction is a **12× enrichment** in class B over
class A. Note this is useless for class A (1.62%, and 0.24% in the
`anchor_min ≤ 3` bin) — mismatch and anchor signals are complementary, not
redundant, and a combined score should use both.

### 5c. Strand agreement per surviving spliced fragment

The library-level strand mixture already implies ~42.3% artifacts. `rigel` can
compute strand agreement **per fragment** against the (contamination-corrected)
protocol strand and use disagreement as artifact evidence. This is a within-`rigel`
signal requiring no new inputs, and it is the same statistic that produced the
independent 42.3% estimate — so its discriminative power is already demonstrated
at the aggregate level.

**Circularity warning:** the strand model is currently *trained* on exonic spliced
fragments (`strand_model.py:391-473`), i.e. on the contaminated population. Using
strand to filter the same population needs either an iterative scheme or a
training set restricted to high-confidence junctions (long anchor, annotated,
unique placement). Fix 5a helps here directly.

## 6. FIX 5 — multi-junction leniency (soundness)

Because anchors are cumulative (§1.1), a fragment with junctions J1 and J2 has,
for J2, an `anchor_left` that includes all bases aligned before J1. The blacklist
test for J2 is therefore evaluated against a larger anchor than the local exon
block, making rejection less likely.

Combined with the precedence rule (§1.3), a single surviving junction shields
every blacklisted junction in the same fragment from the `SPLICE_ARTIFACT` label.
**3,423 junction rows** are affected in the index case.

### Specification

- Evaluate the blacklist test per junction using **locally-bounded** anchors —
  the ref-advancing bases between this junction and the adjacent junction (or
  alignment end) — since that is the quantity `alignable` actually measured from
  single-junction simulated reads. Keep the cumulative value if it is deliberate,
  but then `alignable` must record cumulative anchors too. **Currently the two
  sides disagree about what "anchor" means**, which is the underlying defect.
- Reconsider the precedence rule: a fragment with any blacklisted junction is
  suspect even if another junction survives. At minimum, report these fragments
  in a separate counter so their contribution to the strand and fragment-length
  models is visible.

## 7. Verification harness

The analysis in this investigation re-implements the reject rule in Python and
**reproduces `sj_blacklisted = 50,907` exactly**. Reuse it as a regression check.

```bash
EV=/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_gdna_debug/evidence
BAM=/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_gdna_debug/runs/human/\
mctp_LBX0077_SI_44263_HKWNYDRX7/bam/star.srt.rmdup.collate.bam
PY=/home/mkiyer/sw/miniforge3/envs/rigel/bin/python

# 1. re-quant against a modified rigel / new index  (~25 min, no re-alignment)
rigel quant --index <INDEX> --bam "$BAM" --threads 8 -o <OUT> \
            --annotated-bam <OUT>/annotated.bam

# 2. rebuild the per-junction table and re-measure
$PY $EV/extract_spliced.py <OUT>/annotated.bam <OUT>/spliced_only.bam <OUT>/sj.parquet
$PY $EV/analyze.py          # crosstabs, escape routes, anchor and mismatch profiles
$PY $EV/discriminate.py     # artifact-signature fraction among mRNA-called fragments
```

### Acceptance criteria

| metric | current | target |
|---|---|---|
| `strand_model.strand_specificity` (index case) | **0.7856** | → 0.95+ (clean siblings 0.991–0.999) |
| artifact-signature fraction among mRNA-called spliced fragments | **43.4%** | → < 15% |
| `splice_artifact / all_spliced` (index case) | **0.491** | (rises as detection improves — not a failure) |
| clean-control `mrna_fraction` | 0.995–0.997 | unchanged within 1% |
| clean-control annotated junctions with `maxoh>20 & uniq≥5` | 1,441 / 1,087 | ≥ 99% retained |
| `bug1_halfblind_escapes.tsv` rows now rejected | 0 of 8,617 | majority |
| `bug2_bothsites…tsv` rows now rejected | 0 of 21,617 | majority (via §3) |

`evidence/validation_targets.json` holds these counts machine-readably.

**Primary metric is `strand_specificity`.** It is computed from the surviving
spliced population, is independent of the blacklist, and has a tight
cross-cohort reference range (0.991–0.999 across 20 same-flowcell siblings,
including two at 74% and 82% gDNA). It is the closest thing to ground truth
available without orthogonal experimental data.

## 8. Priority within `rigel`

1. **§3 site-pair blacklist** — measured, +28 points sensitivity, 0.19–0.59% cost,
   needs no new data. Do this first.
2. **§2 dynamic-range guard + treat stored `0` as unknown** — converts silent
   false negatives into reported ones; small, high-value, enables the `alignable`
   fix to land safely.
3. **§5a multimap-aware model training** — directly de-contaminates the strand and
   fragment-length priors, which is where the real downstream damage is.
4. **§4 `--min-splice-anchor`** as the aligner-independent alternative to changing
   STAR's `alignSJDBoverhangMin`. Measure the true-positive cost first.
5. **§6 anchor-semantics alignment with `alignable`** — soundness; coordinate the
   definition across both repositories before either changes.
6. **§5b/5c mismatch and strand scoring** — largest modelling effort, most
   sophistication, addresses class B.
