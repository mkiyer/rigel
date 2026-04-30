# SRD v2 — Handoff for Phases 2–5 (with refined implicit-splice approach)

> **Audience.** A fresh Claude Code instance running on the HPC node, taking
> over after Phases 0 and 1 have shipped. This document is self-contained —
> read it end to end before touching code. It supersedes the Phase 2+ scope
> in `docs/calibration/srd_v2_final_plan.md` in one place: the implicit-splice
> approach (Phase 2c below) is **no longer deferred to v2.1** because the
> elegant detection method makes it almost free.
>
> **Source of truth for what's already done.**
> - `docs/calibration/srd_v2_audit.md` — Phase 0 audit; eight items verified.
> - `docs/calibration/srd_v2_final_plan.md` — the consolidated plan that
>   Phase 1 was implemented against. Pillars 1, 2, 3, 5, 7 still apply
>   verbatim. Pillar 6 (implicit splicing) is **upgraded** by this doc.
> - Latest commit on `main` (`updates to bam scanner and calibration`)
>   ships Phase 1.

---

## 1. Status

| Phase | Status | Reference |
|---|---|---|
| 0 — Audit | ✅ Done | `srd_v2_audit.md` |
| 1 — Calibration-only safety patch | ✅ Done & tested | latest commit on `main` |
| 2 — Scanner overhaul (consolidated) | ⬜ Not started | this doc, §3 |
| 3 — Categorization rewrite (Python) | ⬜ Not started | this doc, §4 |
| 4 — Validation | ⬜ Not started | this doc, §5 |
| 5 — Cleanup | ⬜ Not started | this doc, §6 |

Phase 1 alone eliminates the bin-0 artifact (~33% of EXON_INCOMPATIBLE pool
mass in `mctp_vcap_rna20m_dna00m`). It does **not** address the structural
issue that `INTERGENIC` and `INTRONIC` Python categories are essentially
always 0 because the resolver drops zero-candidate fragments before the
buffer. Phase 2 fixes that.

## 2. Big-picture architecture (recap)

The plan in two sentences:

> Augment the BAM scanner's existing per-fragment cgranges hit-walk to emit
> **strand-aware collapsed overlap counts** and **let zero-candidate fragments
> through the buffer**, so calibration becomes a pure column-arithmetic
> consumer of `FragmentBuffer`. No region partition, no separate region
> tables, no annotation queries downstream of the scanner.

The four new buffer fields per fragment (Pillar 2a from `srd_v2_final_plan.md`,
unchanged):

```
exon_bp_pos    int32   bp of fragment overlapping ANY (+)-strand transcript's exon
exon_bp_neg    int32   bp of fragment overlapping ANY (−)-strand transcript's exon
tx_bp_pos      int32   bp of fragment overlapping ANY (+)-strand transcript's span
tx_bp_neg      int32   bp of fragment overlapping ANY (−)-strand transcript's span
```

Note: we store `tx_bp_*` (transcript-span overlap, including exonic positions),
**not** `intron_bp_*`. The cgranges index has only `ITYPE_EXON` and
`ITYPE_TRANSCRIPT` interval types — no explicit `ITYPE_INTRON` — so
`intron_bp_strand` is **derived in Python** as `max(tx_bp_strand − exon_bp_strand, 0)`.

Categories computed in Python from these four fields:

```
fp                         = genomic_footprint
exon_bp_either_lower_bound = max(exon_bp_pos, exon_bp_neg)
tol                        = exon_fit_tolerance_bp  (default 3)

INTERGENIC       : tx_bp_pos == 0  AND  tx_bp_neg == 0
INTRONIC         : exon_bp_pos == 0  AND  exon_bp_neg == 0  AND  not INTERGENIC
EXON_CONTAINED   : (fp - exon_bp_either_lower_bound) <= tol
EXON_INCOMPATIBLE: otherwise
```

Strand sub-label per category (none / pos / neg / ambig).

## 3. Phase 2 — Scanner overhaul (one consolidated C++ build)

This phase touches `src/rigel/native/resolve_context.h` and
`src/rigel/native/bam_scanner.cpp`. **One recompile**, three sub-changes
folded together:

### 3a. Strand-aware collapsed overlap counts

`src/rigel/native/resolve_context.h` — augment `RawResolveResult` (struct
near line 60ish) with the four new int32 fields, and the inner loop at
**lines ~880–915** with the strand-routing logic:

```cpp
// Existing per-cgranges-hit code already computes:
//   bp = max(0, min(bend, h_end) - max(bstart, h_start))
//   t_set = t_set_data_[off : off + cnt]
// Cgranges intervals are non-overlapping by construction (collapsed
// exon/transcript tiling), so simple addition is correct.
bool any_pos = false, any_neg = false;
for (int32_t k = 0; k < cnt; k++) {
    int8_t s = t_strand_arr_[t_set_data_[off + k]];
    any_pos |= (s == STRAND_POS);
    any_neg |= (s == STRAND_NEG);
}
if (itype == ITYPE_EXON) {
    if (any_pos) cr.exon_bp_pos += bp;
    if (any_neg) cr.exon_bp_neg += bp;
}
if (itype == ITYPE_TRANSCRIPT) {  // NB: TRANSCRIPT intervals already
    if (any_pos) cr.tx_bp_pos += bp;   //     subsume the exon footprint
    if (any_neg) cr.tx_bp_neg += bp;   //     of the same strand.
}
```

Strand-overlap zone semantics: when a single cgranges interval has both
+ and − transcripts in its `t_set`, `bp` is added to both pos and neg.
Therefore `exon_bp_pos + exon_bp_neg ≥ exon_bp_either` always, and Phase 3
treats positions covered by both strands as `strand_ambig`.

### 3b. Stop dropping zero-candidate fragments

Currently `resolve_context.h:1014–1016` early-exits when `cr.t_inds.empty()`,
preventing INTERGENIC fragments from ever reaching the buffer. Remove that
early exit. A zero-candidate fragment now flows through the buffer with:

- empty `t_inds` (already handled downstream by every consumer that checks
  `n_cands > 0`)
- `genomic_footprint > 0`
- `exon_bp_pos = exon_bp_neg = tx_bp_pos = tx_bp_neg = 0`
- categorized as INTERGENIC in Phase 3

The existing `intergenic_obs / intergenic_truth / intergenic_lengths`
accumulators in `bam_scanner.cpp` (~lines 227–238 and 1320–1326) become
**redundant** but leave them in place for Phase 2 — Phase 5 cleanup
removes them. Cross-check: `n_intergenic_unspliced` from the accumulator
should match the new INTERGENIC category count (modulo the multimap filter).

**Verify before merging:** `locus.py` and the per-locus EM never iterate a
fragment with empty `t_inds` (they only project to transcripts via
`t_indices`). `scan.py FragmentRouter` may need a one-line skip-on-zero
guard; check it. The audit (`srd_v2_audit.md` §5) already swept these
consumers.

### 3c. SPLICE_ARTIFACT classification (consolidated here)

Add `SpliceType.SPLICE_ARTIFACT = 4` in `src/rigel/types.py`. In the
resolver where splice_type is currently set on the blacklist-rejection path
(`resolve_context.h` ~lines 1006–1008, the `has_unannotated_sj ? ... : SPLICE_UNSPLICED`
ternary), intercept the case where `n_sj_blacklisted > 0` was set by the
upstream `filter_blacklisted_sjs` (per-read counter at
`bam_scanner.cpp:178–181`) and assign `cr.splice_type = SPLICE_ARTIFACT`
instead of downgrading to `UNSPLICED`.

The `n_sj_blacklisted` counter is already threaded per-read; just propagate
it into the resolver context. **No buffer schema change needed** — the
audit (`srd_v2_audit.md` §6) confirmed there is no separate "block merging"
step to remove. The merged span in `genomic_footprint` is the natural
side effect of the existing footprint computation; for SPLICE_ARTIFACT
fragments it's inflated, but since they're held out of calibration this
is harmless.

### 3d. **NEW: SPLICED_IMPLICIT — the elegant approach**

> **This supersedes the "deferred to v2.1" disposition in `srd_v2_final_plan.md` §Pillar 6.**

The original final-plan deferral assumed implicit-splice detection required
a new annotated-introns range index (because `sj_map_` is keyed by exact
intron coordinates and doesn't support range queries). That assumption was
wrong: the existing `compute_frag_lengths` machinery (resolve_context.h:725–757)
**already implicitly handles intron-in-gap** through transcript-coordinate
projection.

**The mechanism.** For a multi-block fragment (paired-end with gap), `compute_frag_lengths`
takes the outer extents `gstart`/`gend` and projects them onto each candidate
transcript's spliced coordinate space via `genomic_to_tx_pos`
(resolve_context.h:678–703). When the unsequenced gap fully contains an
annotated intron `I` of transcript `T`, the projection naturally skips that
intron and `tx_distance(T) ≈ genomic_span − len(I)`. The shortening **is**
the implicit-splicing signal.

**Worked example.** Transcript T with exons (500,1000), (5000,5100), (7000,7400).
Multi-block fragment with R1=(700,900) and R2=(5020,5100):
- `gstart=700`, `gend=5100`, `genomic_span = 4400`
- `tx_s = genomic_to_tx_pos(700)  = 200`  (200bp into exon 1)
- `tx_e = genomic_to_tx_pos(5100) = 600`  (cumsum[1] + 100 = end of exon 2)
- `tx_distance = 400`
- `genomic_span − tx_distance = 4000` ≡ length of intron 1 → SPLICED_IMPLICIT.

**Implementation (~5 lines, no new data structures):**

In the resolver, after `compute_frag_lengths` populates `frag_length_map`
for a multi-block fragment, before the final splice_type assignment:

```cpp
// SRD v2: implicit-splice detection via existing transcript projection.
// For multi-block (paired-end gap) fragments, if the gap fully spans
// an annotated intron of any candidate transcript, the projected
// tx_distance for that transcript will be < genomic_span by the
// intron's length. We use this as a "free" implicit-splice detector.
if (cr.splice_type == SPLICE_UNSPLICED && exons.size() >= 2) {
    int32_t genomic_span = gend - gstart;
    constexpr int32_t MIN_INTRON_THRESHOLD = 30;  // smaller than any real intron
    for (const auto& [t, fl] : cr.frag_length_map) {
        if (fl > 0 && (genomic_span - fl) >= MIN_INTRON_THRESHOLD) {
            cr.splice_type = SPLICE_IMPLICIT;
            break;
        }
    }
}
```

Add `SpliceType.SPLICED_IMPLICIT = 3` in `src/rigel/types.py`.

**Conservative "ANY candidate" rule.** If any candidate transcript shows
shortening, classify as SPLICED_IMPLICIT — this means a fragment consistent
with mature mRNA on T1 stays out of the gDNA pool even if T2 has different
geometry. This is the right policy: fragments that *could* be RNA must not
be treated as gDNA.

**Single-block fragments** never have implicit splicing — there's no gap
to contain an intron. The single-block path (resolve_context.h:735–739)
sets `fl = block_length` for all candidates; no projection runs. The
`exons.size() >= 2` guard above skips them cheaply.

**Optional but clean:** for SPLICED_IMPLICIT fragments, recompute
`genomic_footprint` to be `sum(block_length)` rather than `max_end − min_start`,
so downstream consumers see a physically meaningful span. Not strictly
required for calibration since SPLICED_IMPLICIT is held out of the pool,
but the diagnostic counters become more interpretable.

### 3e. Buffer schema additions

`src/rigel/buffer.py` (`_FinalizedChunk` and `from_raw`): add the four new
columns. Update `memory_bytes`. Cost: 16 bytes/fragment.

`src/rigel/types.py`: add `SpliceType.SPLICED_IMPLICIT = 3` and
`SpliceType.SPLICE_ARTIFACT = 4`.

### 3f. Phase 2 tests

- `tests/test_resolve_strand_overlap.py` — synthetic mini-annotation with
  fragments crossing exon/intron/intergenic boundaries on +/− strands;
  assert exact bp counts; include the strand-overlap zone case (a position
  covered by both + and − transcripts); verify the invariants:
  - `exon_bp_pos ≤ tx_bp_pos`
  - `exon_bp_neg ≤ tx_bp_neg`
  - `max(tx_bp_pos, tx_bp_neg) ≤ genomic_footprint`

- `tests/test_zero_candidate_fragments.py` — synthetic intergenic fragments
  now reach the buffer with empty `t_inds` and zeroed strand counts;
  downstream consumers (FragmentRouter, locus builder) tolerate them.

- `tests/test_splice_artifact.py` — blacklisted SJ → `SPLICE_ARTIFACT`;
  not downgraded to `UNSPLICED`.

- `tests/test_implicit_splicing.py` — synthetic paired-end fragments:
  - Gap containing an intron of length 4000 → SPLICED_IMPLICIT (positive case).
  - Gap containing an intron of length 20 (below MIN_INTRON_THRESHOLD) →
    stays UNSPLICED (negative case, threshold edge).
  - Gap not containing any intron of any candidate (e.g., gap is fully
    intergenic) → stays UNSPLICED.
  - Gap containing intron of T1 but not T2 (T2 is exonic in the gap region)
    → SPLICED_IMPLICIT (conservative "ANY candidate" rule).
  - Single-block fragment overhanging an exon-intron boundary → stays
    UNSPLICED (no gap to contain an intron; takes the single-block path).

### 3g. Phase 2 exit criteria

- `pytest tests/ -v` green.
- Smoke run on `mctp_vcap_rna20m_dna00m`: scanner runtime delta < 5%; new
  buffer columns populated; INTERGENIC and INTRONIC category counts now
  non-zero (cross-check with `n_intergenic_unspliced` from the legacy
  accumulator — should match modulo multimap filter).

## 4. Phase 3 — Categorization rewrite (Python only)

No recompile. Single PR. Files:

### 4a. `src/rigel/calibration/_categorize.py` — full rewrite

Stop reading `cand_frag_lengths`, `t_offsets`, `t_indices`, `intron_bp`,
`exon_bp` (the per-candidate arrays). Read only:
- `splice_type`
- `num_hits`
- `genomic_footprint`
- `exon_bp_pos`, `exon_bp_neg`, `tx_bp_pos`, `tx_bp_neg` (new)

Categories per Pillar 3 of `srd_v2_final_plan.md`:

```python
# Filter: unique-mappers AND truly-unspliced
keep = (chunk.num_hits == 1) & (chunk.splice_type == SpliceType.UNSPLICED)

fp = chunk.genomic_footprint
exon_bp_pos = chunk.exon_bp_pos
exon_bp_neg = chunk.exon_bp_neg
tx_bp_pos   = chunk.tx_bp_pos
tx_bp_neg   = chunk.tx_bp_neg

intron_bp_pos = np.maximum(tx_bp_pos - exon_bp_pos, 0)
intron_bp_neg = np.maximum(tx_bp_neg - exon_bp_neg, 0)

exon_either_lower = np.maximum(exon_bp_pos, exon_bp_neg)
tol = exon_fit_tolerance_bp  # default 3

intergenic = (tx_bp_pos == 0) & (tx_bp_neg == 0)
intronic   = (exon_bp_pos == 0) & (exon_bp_neg == 0) & ~intergenic
exon_contained = (fp - exon_either_lower) <= tol
exon_incompat  = ~(intergenic | intronic | exon_contained)
```

Strand sub-label (per category — for `INTRONIC` use `intron_bp_*`; for the
exon classes use `exon_bp_*`):

```python
def strand_label(bp_pos, bp_neg):
    pos = (bp_pos > 0) & (bp_neg == 0)
    neg = (bp_neg > 0) & (bp_pos == 0)
    ambig = (bp_pos > 0) & (bp_neg > 0)
    return where(pos, POS, where(neg, NEG, where(ambig, AMBIG, NONE)))
```

Output `category × strand` 2-D counts.

### 4b. `src/rigel/calibration/_result.py`

Widen `category_counts` to shape `(N_CATEGORIES, 4)` with last axis
`(none, pos, neg, ambig)`. JSON-serialize as a flat list with documented
stride.

Add diagnostic fields:
- `n_pool_intronic_strand_pos`, `n_pool_intronic_strand_neg` (for nRNA
  pollution surfacing — see §7 below).

### 4c. `src/rigel/calibration/_simple.py`

Pool: `Pool = INTERGENIC ∪ INTRONIC ∪ EXON_INCOMPATIBLE` (unchanged from
v1 in spirit; the categories are now genuinely populated).

**Remove** the `frag_length_models.intergenic` merge (~lines 138–157 in
the current `_simple.py`). Phase 2b (zero-candidate-drop removal) means
those fragments now flow through the buffer directly and are categorized
as INTERGENIC. Pulling them in twice would double-count.

Length: `chunk.genomic_footprint` (already wired by Phase 1; nothing to
change there).

### 4d. Phase 3 tests

- `tests/test_categorize_v2.py` — synthetic per-category and per-strand
  cases (4 categories × 4 strand sub-labels = 16 base cases) plus
  tolerance edges (overhang of exactly `tol` and `tol + 1`) plus the
  strand-overlap zone case.
- Update `tests/test_calibration_simple.py` for the new schema.
- Update goldens via `pytest --update-golden` after eyeballing diffs.

### 4e. Phase 3 exit criteria

`pytest tests/ -v` green; smoke run on `mctp_vcap_rna20m_dna00m` shows:
- INTERGENIC and INTRONIC counts now substantial (10× to 100× larger than
  before, since they were near-zero pre-fix).
- Pool size grows correspondingly.
- `pi_pool` may shift meaningfully — record before/after.

## 5. Phase 4 — Validation

Rerun the spike series end-to-end:

```bash
bash /scratch/.../rigel/srd_v1_vcap_mixture/run_all.sh
```

Output to a sibling `srd_v2_vcap_mixture/default/` directory and produce
`docs/calibration/srd_v2_results.md` with side-by-side tables vs the
SRD v1 baseline (`srd_v1_results.md`).

### 5a. Quantitative acceptance

| Metric | Library | Threshold |
|---|---|---|
| `gDNA_FL` bin-0 mass | all | < 1% |
| `gDNA_FL` bin-max mass | all | < 1% |
| `gDNA_FL` mode | all | in `[100, 400]` bp |
| `gDNA_FL` mode consistency | spiked libs | within ±50 bp across `dna01m..dna80m` |
| Headline `gdna_fraction` monotonicity | dna00m → dna80m | strictly increasing |
| Headline `gdna_fraction` at dna80m (~80% nominal) | spiked | ≥ 0.75 (improvement vs SRD v1's ~0.69) |
| Headline `gdna_fraction` at dna00m | "pure" RNA | report only; no pass/fail. 0.5–1% is acceptable. |
| Pool size growth from Pillar 2b | dna00m | report (expect 10×–100× with INTERGENIC + INTRONIC populated) |
| `n_pool_dropped_out_of_range` | all | < 0.1% of pool |
| Scanner runtime delta vs SRD v1 | full sample | < 5% |
| `n_pool_implicit_splice_caught` (new diagnostic) | all | report — expect a few % of multi-block UNSPLICED |

### 5b. Framing reminders

- `mctp_vcap_rna20m_dna00m` is the real-world floor; no ground-truth 0%
  exists. 0.5–1% headline `gdna_fraction` is acceptable. We are testing
  *shape* and *trend*, not zero.
- If dna80m headline doesn't improve, the structural nRNA limitation
  (§7) is the suspect. Drill in with `INTRONIC × strand` diagnostic.
- `gdna_fl_model_quality` may legitimately be `"weak"` or `"fallback"` on
  the unspiked library — that's correct behavior, not a bug.

## 6. Phase 5 — Cleanup

Pure deletion. Sequence after Phase 4 acceptance.

- Remove the C++ intergenic accumulator path in `bam_scanner.cpp`
  (`intergenic_obs / intergenic_truth / intergenic_lengths`,
  ~lines 227–238 and 1320–1326). Pillar 2b makes it redundant.
- Remove the `frag_length_models.intergenic` plumbing in Python
  (`pipeline._replay_fraglen_observations` and its caller).
- Remove the `n_intergenic_unspliced` cross-check counter once parity is
  confirmed (or keep as a permanent telemetry counter — your call).
- Update `docs/MANUAL.md`, `docs/METHODS.md`, `CHANGELOG.md` to describe
  SRD v2 as shipped. Cite this handoff doc and `srd_v2_results.md`.
- Bump version (minor; no breaking API change for end users; new
  diagnostic fields in `CalibrationResult`).

## 7. Accepted v2 limitations (document, don't try to fix)

### 7a. nRNA pool pollution

> **UPDATE 2026-04-28 (SRD v3 Phase 1 landed):** the diagnostic
> referenced below was renamed to transcript frame —
> `n_pool_intronic_strand_sense / _antisense / _ambig`. The old
> `n_pool_intronic_strand_pos / _neg` reflected which genomic strand
> the local gene happened to live on (a coordinate artifact); the new
> labels compare the read's own genomic strand against the
> overlapping transcript's strand and are therefore directly
> interpretable as a nRNA-vs-gDNA strand-asymmetry signal. The full
> resolution of this limitation (a joint FL × strand mixture) is
> scheduled for v3 Phase 2; see `srd_v3_early_plan.md` and
> `srd_v3_phase1_strand_relabel.md`.

The 1-D mixture `pool ∝ π·gDNA + (1−π)·RNA_FL` cannot distinguish gDNA
from nascent RNA when their fragment-length distributions overlap, which
they typically do (both span 100–600 bp from sonication). nRNA fragments
enter the pool (almost entirely INTRONIC; partly EXON_INCOMPATIBLE),
inflating fitted `π`.

This is a **structural identifiability limit**, not a bug — fixing it
requires a strand-aware mixture or a third component, both out of scope
for v2. Surface it via the new `n_pool_intronic_strand_sense` /
`n_pool_intronic_strand_antisense` diagnostic so the user can spot
nRNA-heavy libraries (extreme strand asymmetry inside the INTRONIC
bucket given high SS).

### 7b. SPLICE_ARTIFACT block-disjunction not preserved in buffer

Per audit §6, there is no separate block-merging step to remove. For v2
the merged `genomic_footprint` for SPLICE_ARTIFACT fragments is acceptable
because they are held out of calibration. v2.1 may add per-block
preservation if downstream needs it; for now, no plumbing change.

## 8. Risk register

| Risk | Mitigation |
|---|---|
| Phase 2 strand-bp accounting wrong at strand-overlap zones | `test_resolve_strand_overlap.py` invariant tests; the `any_pos`/`any_neg` flags are conservative |
| Removing zero-candidate drop breaks downstream | Audit §5 enumerated consumers; one-line guard at most in `scan.py FragmentRouter` |
| Implicit-splice detector false-positives at small introns | `MIN_INTRON_THRESHOLD = 30` is below any real intron; the conservative "ANY candidate" rule is intentional |
| Implicit-splice detector misses cases where no candidate has the relevant intron | Acceptable — those fragments still take the per-locus EM path, where they may classify correctly via FL likelihood |
| Scanner overhead from 4 new fields | Profile in Phase 4; expect <5% (work is constant-factor on existing hit-walk) |
| nRNA contamination still pollutes gDNA_FL | Accepted limitation §7a; new diagnostic surfaces it |
| `n_pool_dropped_out_of_range` materially > 0.1% | Phase 1 already validates this is small for VCaP libraries; investigate in Phase 4 if it spikes |
| Pure-RNA library has no real gDNA signal | Expected; `gdna_fl_quality = "weak"` or `"fallback"` correctly fires |
| SPLICED_IMPLICIT detection misses gap-hidden introns where compute_frag_lengths returned 0 (sentinel skip) | Worth checking: if a candidate's projection collapses to 0, the implicit-splice check skips it. Conservative — those fragments go to per-locus EM. May warrant follow-up if Phase 4 shows residual gDNA pool contamination from this path. |

## 9. Out of scope for v2 (deliberate)

- SS-gated antisense-exonic pool inclusion (largest open lever for
  pure-library sensitivity per `srd_v2_final_plan.md` §2 Pillar 4 note).
  v2.1 follow-up.
- Strand-aware 2D mixture for nRNA / gDNA separation.
- Cross-region density signals; mappability-corrected effective length;
  region partition resurrection.
- Strand-model / RNA-FL training changes.
- Per-block preservation for SPLICE_ARTIFACT.

## 10. References

- `docs/calibration/srd_v2_audit.md` — Phase 0 audit results.
- `docs/calibration/srd_v2_final_plan.md` — Pillars 1–7 (Pillar 6 superseded by §3d above).
- `docs/calibration/srd_v2_bugs_and_new_plan.md` — original bug analysis.
- `docs/calibration/srd_v1_results.md` — SRD v1 baseline for Phase 4 comparison.
- Latest commit on `main` (`updates to bam scanner and calibration`) — Phase 1 implementation.
