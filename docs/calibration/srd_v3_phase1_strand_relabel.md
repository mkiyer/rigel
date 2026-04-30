# SRD v3 Phase 1 — Strand Relabeling: Implementation Plan

**Date:** 2026-04-28
**Scope:** Replace genomic-frame POS/NEG strand labels in calibration
with transcript-frame SENSE/ANTISENSE/AMBIG labels. Bookkeeping +
diagnostic only — no model change.
**Roadmap context:** `srd_v3_early_plan.md` §4.1 Phase 1.

---

## 1. Problem (one paragraph)

`_categorize._strand_label()` reports which **transcript strand** has
overlap (`exon_bp_pos > 0` etc.) without ever consulting the read's
own strand (`chunk.exon_strand`). The resulting `StrandLabel.POS / NEG`
on, e.g., the INTRONIC bucket therefore reflects which genomic strand
the local gene happens to live on — a coordinate-system artifact, not
a biological signal. To diagnose nRNA pollution of the INTRONIC pool
(and to enable the v3 Phase 2 joint mixture later) we need labels in
the **transcript frame**: sense ↔ read agrees with overlapping
transcript strand, antisense ↔ read disagrees.

## 2. New enum

`src/rigel/calibration/_categorize.py`:

```python
class FragmentStrand(IntEnum):
    """SRD v3 transcript-frame strand label.

    Defined relative to the strand of the transcript(s) overlapping the
    fragment. For INTERGENIC fragments (no transcript overlap), the
    enum is reused as a pure naming convention: SENSE = read on (+)
    genomic strand, ANTISENSE = read on (−) genomic strand. The
    INTERGENIC labeling carries no transcript-frame meaning and downstream
    consumers should treat the INTERGENIC row of category_counts as
    aggregate-only.
    """
    SENSE     = 0   # read strand == overlapping-transcript strand
                    # (or "+" by convention for INTERGENIC)
    ANTISENSE = 1   # read strand != overlapping-transcript strand
                    # (or "−" by convention for INTERGENIC)
    AMBIG     = 2   # transcripts on BOTH strands overlap (sense to one,
                    # antisense to another); cannot classify
```

`StrandLabel` (the old enum, 4 values incl. `NONE`) is **removed**.
The `NONE` slot disappears: every unique-mapper UNSPLICED fragment with
a known read strand receives one of the 3 new values. There is no
`UNKNOWN` value because read strand is always observed in the BAM —
what is probabilistic is the inference from genomic strand to RNA
template strand, and that is a downstream (per-source likelihood)
operation, not a labeling decision.

`N_STRAND_LABELS` is renamed `N_FRAGMENT_STRANDS` and becomes 3.

## 3. Classification rule (vectorized)

For each unique-mapper UNSPLICED fragment:

```text
read_strand  = chunk.exon_strand        # POS=1 / NEG=2 / UNKNOWN=0 / AMBIG=3
ebp_pos, ebp_neg                        # transcript-strand exon bp counts
ibp_pos, ibp_neg = max(tbp_* − ebp_*, 0)  # derived per-strand intron bp

For exon classes (EXON_CONTAINED, EXON_INCOMPATIBLE):
    bp_pos, bp_neg = ebp_pos, ebp_neg
For INTRONIC:
    bp_pos, bp_neg = ibp_pos, ibp_neg
For INTERGENIC:
    bp_pos = bp_neg = 0 (by definition)

label = AMBIG     if (bp_pos > 0 AND bp_neg > 0)                # tx on both strands
        SENSE     if (read_strand == POS AND bp_pos > 0 AND bp_neg == 0)
                  or (read_strand == NEG AND bp_neg > 0 AND bp_pos == 0)
        ANTISENSE if (read_strand == POS AND bp_neg > 0 AND bp_pos == 0)
                  or (read_strand == NEG AND bp_pos > 0 AND bp_neg == 0)
        SENSE     if (INTERGENIC AND read_strand == POS)        # convention
        ANTISENSE if (INTERGENIC AND read_strand == NEG)        # convention
```

### 3.1 Edge cases

**Read strand AMBIGUOUS (R1/R2 disagree).** `chunk.exon_strand` carries
`STRAND_AMBIGUOUS = 3` in this case (already excluded from
`StrandModel` training). For calibration, these are dropped from the
per-strand label assignment entirely — they fail the `keep` mask,
landing in the same exclusion bin as multi-mappers and non-UNSPLICED
fragments. Population is small in real data (≤ 1% of unique-mapper
unspliced fragments per spot-check); no information loss.

> Implementation: extend the `keep` predicate in `categorize_chunk`
> from `(num_hits == 1) & (splice_type == UNSPLICED)` to also require
> `(exon_strand == STRAND_POS) | (exon_strand == STRAND_NEG)`.

**Read strand UNKNOWN (`exon_strand == 0`).** Possible if the C++
scanner couldn't infer strand at all (rare; mostly a programming
sentinel). Treated identically to AMBIGUOUS: dropped from `keep`.

**INTERGENIC fragment with AMBIG read strand.** Filtered out by `keep`
above. This means INTERGENIC fragments that survive labeling all carry
SENSE or ANTISENSE; the row therefore acts as a 50/50 sanity check on
genomic strand assignment for gDNA fragments (which is, in practice, a
useful diagnostic — a gross deviation from 50/50 would indicate an
alignment-strand bias).

**Strand-overlap zone (transcripts on both genomic strands at the same
position).** Both `bp_pos > 0` and `bp_neg > 0`. Fragment is, by
definition, sense to one transcript and antisense to another →
AMBIG. Already handled by the rule above.

**INTRONIC with derived `ibp_pos == ibp_neg == 0`.** Cannot happen by
construction: INTRONIC requires `exon_bp_*` both zero AND `tx_bp_*` not
both zero, so at least one of `ibp_*` is positive. (Sanity-asserted in
test.)

## 4. File-by-file edits

| File | Change |
|---|---|
| `src/rigel/calibration/_categorize.py` | Replace `StrandLabel` → `FragmentStrand` (3 values). Rename `N_STRAND_LABELS` → `N_FRAGMENT_STRANDS`. Rewrite `_strand_label()` → `_compute_fragment_strand()` taking `(bp_pos, bp_neg, read_strand, intergenic_mask)`. Tighten `keep` mask to require known read strand. Update module docstring + `CategorizedChunk.strand` docstring. |
| `src/rigel/calibration/_simple.py` | Update import. Replace `n_pool_intronic_strand_pos / _neg` accumulators with `_sense / _antisense`. Update logger format string. Pass renamed values to `CalibrationResult`. |
| `src/rigel/calibration/_result.py` | Rename `n_pool_intronic_strand_pos / _neg` → `_sense / _antisense`. Update `category_counts` shape comment to `(N_CATEGORIES, 3)`. Update docstring to spell out the transcript-frame semantics and the INTERGENIC convention reuse. Update `as_dict()` JSON keys. |
| `tests/test_categorize.py` | Replace `StrandLabel` import. Rewrite the parameterised `test_*_strand_labels` cases to drive `exon_strand` explicitly and assert SENSE/ANTISENSE/AMBIG outputs. Add new test cases for the read-strand-aware logic (sense vs antisense flip, AMBIG read strand drops fragment). Drop `test_intergenic_zero_strand_bp` strand-NONE assertion; replace with the SENSE/ANTISENSE convention test. |
| `tests/test_calibration_simple.py` | Update `StrandLabel` import. The `n_pool` denominator tests already subtract `n_pool_dropped_out_of_range`; if the new `keep` predicate further reduces the pool, these tests need to track it through `category_counts.sum()` (which already excludes filtered-out fragments). Verify and adjust. |
| `docs/calibration/srd_v2_phase2plus_handoff.md` §7a | Add an inline note that v3 Phase 1 makes the INTRONIC strand diagnostic interpretable (transcript frame); v3 Phase 2 will retire the limitation entirely. |

**No C++ changes.** All required inputs (`exon_strand`, `exon_bp_pos/neg`,
`tx_bp_pos/neg`) are already in the buffer.

## 5. Schema changes (visible to users)

| Surface | Old | New |
|---|---|---|
| `CalibrationResult.category_counts` | shape `(4, 4)` | shape `(4, 3)` |
| `CalibrationResult.n_pool_intronic_strand_pos` | int | **removed** |
| `CalibrationResult.n_pool_intronic_strand_neg` | int | **removed** |
| `CalibrationResult.n_pool_intronic_strand_sense` | — | int (new) |
| `CalibrationResult.n_pool_intronic_strand_antisense` | — | int (new) |
| `CalibrationResult.n_pool_intronic_strand_ambig` | — | int (new) |
| `summary.json` keys | as above | as above |

These are deliberate breaks; downstream tools should either be updated
or pin to the previous version. Acceptable while v3 is in flight.

## 6. Validation

1. **Unit tests** for the new classification rule — small parametrised
   table per (`category`, `read_strand`, `bp_pos`, `bp_neg`) → expected
   `FragmentStrand`.
2. **Existing tests** must pass: golden output, pipeline smoke,
   calibration regression. The categorical aggregate counts may shift
   (some POS/NEG fragments redistribute into SENSE/ANTISENSE) but the
   total `n_pool` and `pi_pool` should be **identical** because no FL
   data has changed and the mixture EM uses neither the strand label
   nor the `keep` predicate beyond category counts. Any drift is
   attributable to the slightly tighter `keep` mask (drops AMBIG-read
   fragments); expected to be sub-percent.
3. **Field-level golden updates** if needed: regenerate via
   `pytest tests/ --update-golden` and inspect diffs.
4. **Smoke test on a real BAM** (one of the VCaP libraries on
   armis2) — verify that on a stranded library, the INTRONIC
   sense-fraction is dramatically asymmetric in transcript frame
   (expected) vs. ~50/50 in the old genomic frame.

## 7. Acceptance criteria

- All `pytest tests/` pass (with goldens regenerated where the strand
  re-labeling shifts category-count distributions; mixture-EM-derived
  fields like `pi_pool` should be unchanged within numerical noise).
- `n_pool` is unchanged ± the small AMBIG-read drop (typically < 1%).
- `pi_pool`, `gdna_fl_*`, `mixture_iterations` **unchanged within
  ε = 1e-6** on synthetic test scenarios (these never produce AMBIG
  read strands).
- A new diagnostic test confirms that on a stranded synthetic
  scenario, INTRONIC SENSE > INTRONIC ANTISENSE.

## 8. Implementation order

1. Edit `_categorize.py` — new enum, classification rule, tightened
   `keep`. Run `pytest tests/test_categorize.py` (will fail) and update
   tests to match new semantics.
2. Edit `_simple.py` — accumulators, logger, result construction.
3. Edit `_result.py` — schema rename + docstrings.
4. Run full `pytest tests/`. Regenerate goldens if needed.
5. Update `srd_v2_phase2plus_handoff.md` §7a footnote.
6. Hand-spot-check on one real BAM (optional but recommended).

## 9. Out of scope for Phase 1

- The actual joint (FL × strand) mixture (Phase 2).
- Any change to `StrandModel`, `FragmentLengthModels`, or
  `_fl_mixture.py`.
- Any C++ change.
- Schema additions for the future SS-weighted mixture.

## 10. References

- `srd_v3_early_plan.md` — roadmap, §3 theoretical foundation.
- `src/rigel/calibration/_categorize.py` — current code under
  modification.
- Conversation log capturing the design decision to use Python-side
  frame transformation (vs. pushing into C++) on the basis of code
  simplicity and future evolutionary flexibility.
