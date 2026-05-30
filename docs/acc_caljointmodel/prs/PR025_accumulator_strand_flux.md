# PR 2.5 — Accumulator strand & spliced-flux correction

**Parent plan:** [`../00_implementation_plan.md`](../00_implementation_plan.md) §4 D2/D7, §7.
**Type:** **C++ (accumulator + scanner) + Python (payload, substrate).**
**Build required:** **yes** (`pip install --no-build-isolation -e .`).
**Status:** Design settled (PR03 §0 review). Magic-number budget: **0 new
constants** (channel encoding only).

This PR corrects the accumulator's **strand** and **boundary-flux** semantics
so the RNA strand model (PR 3) and the AMBIG-region handling (D7) rest on
correct sufficient statistics. It precedes PR 3. No calibration inference here.

## Why (from the PR 3 review)

1. **Spliced strand is motif-defined and region-independent.** A spliced
   fragment's RNA strand is fixed by the splice motif (GT/AG → XS/ts tag), so it
   is valid even inside strand-ambiguous regions. Orienting at deposit
   (`sense = align_strand == sj_strand`) lets annotated spliced fragments train
   the RNA strand model anywhere (D2/D7).
2. **Boundary flux must be per-side for spliced.** A spliced fragment skipping
   an intron lands on two *different* boundaries (exon→intron, intron→exon); a
   single shared flux falsely credits both sides of each. Per-side flux
   (deposit on the side carrying the fragment's block) fixes it — and, because
   it just follows the slice geometry, an unspliced contiguous crossing
   naturally credits **both** sides of its one boundary (the old shared
   behaviour) while a spliced jump credits **one** side of each.
3. **The accumulator emits genome strand, not transcript-relative.** Unspliced
   orientation (→ sense) belongs downstream (region strand + library mode);
   AMBIG unspliced stays unrecoverable by strand (→ density/sweep, D7).
4. **`exon_strand` is a poor name** for the read's genome alignment strand.

## Changes

### T1 — Rename `exon_strand` → `align_strand` (mechanical, no behaviour change)

Across C++ (`constants.h`, `resolve_context.h`, `resolve.cpp`,
`bam_scanner.cpp`) and Python (`buffer.py`, `strand_model.py`, the resolve
result key, the buffer Arrow column) + tests. Pure refactor; its own commit.
`align_strand` = the read's genome alignment strand (distinct from `sj_strand`
= splice-motif strand). Buffers are transient, so the column rename has no
on-disk-compat cost.

### T2 — Channel encoding: spliced = motif-relative sense/antisense

Channel layout (`kNChannels = 4`, `ch = (spliced?2:0) + (primary?0:1)`):

```
ch0 = unspliced & align_strand +     (genome strand)
ch1 = unspliced & align_strand −
ch2 = spliced   & SENSE              (align_strand == sj_strand)
ch3 = spliced   & ANTISENSE          (align_strand != sj_strand)
```

The accumulator becomes channel-agnostic: `deposit(..., int channel)` (drops
the `spliced`/`strand_pos` bools; `channel_idx(spliced, primary)` stays as a
helper). **The scanner computes the channel**:

- unspliced: `primary = (align_strand == STRAND_POS)`  → genome strand.
- spliced:   `primary = (align_strand == sj_strand)`   → motif-relative sense.

**Deposit gate (documented):** unspliced deposits iff `align_strand ∈ {POS,NEG}`
(unchanged); spliced deposits iff `align_strand ∈ {POS,NEG}` **and**
`sj_strand ∈ {POS,NEG}` (motif present) — orienting by the motif, valid in
AMBIG regions. Spliced fragments without a motif strand are not deposited (they
cannot inform the RNA strand model; a negligible fraction). The
*annotated-only* restriction for the strand **fit** is enforced in PR 3 by
which boundary observations it pools (the partition's exon/intron boundaries
are annotated junctions by construction).

### T3 — Per-side boundary flux

`Boundary.flux[4]` → `flux_left[4]` + `flux_right[4]` (struct 48 B → 64 B). In
`deposit`, the crossing path increments `flux_left[ch]` on `b_out` (the
left-region's right boundary, matching its `mass_left`) and `flux_right[ch]` on
`b_in` (the right-region's left boundary, matching its `mass_right`) — replacing
the old dedup-and-shared-`flux` step. `merge_from` sums both. No `spliced` flag
needed; the slice geometry distinguishes contiguous (both sides) from
intron-skip (one side each).

### T4 — Payload schema (`scan_payload.py`)

`boundary_flux` → `boundary_flux_left` + `boundary_flux_right` (both
`uint32[B_obj, 4]`). Update the scanner payload build (cal dict) and
`AccumulatorPayload`. Document the new channel semantics in the schema header.
(No index/on-disk change — this is the scan payload only.)

### T5 — Substrate rework (`substrate.py`, raw genome strand; orient downstream)

`SubstrateView` exposes **raw** counts, no pre-oriented `k_plus`:

```
n_unspliced_pos, n_unspliced_neg   # int64[R] — genome strand (orient downstream)
n_spliced_sense, n_spliced_antisense  # int64[R] — motif-relative (ready to use)
mass_unspliced, mass_spliced       # float64[R] — strand-agnostic magnitude
```

Boundary views use per-side flux: left view ← `flux_right[left_boundary(r)]`,
right view ← `flux_left[right_boundary(r)]` (the side facing region `r`), with
mass from `mass_right`/`mass_left` as before. The E-step (PR 4) orients
unspliced via region strand + library mode; AMBIG unspliced has no valid
orientation (D7). The spliced sense/antisense are already transcript-relative.

### T6 — Tests + rebuild

- `tests/native/test_accumulator_spec.py` — per-side flux (contiguous ⇒ both
  sides; intron-skip ⇒ one side each, no false intron flux); channel-by-index
  deposit. Re-baseline the reference if needed.
- `test_accumulator_payload.py` — `boundary_flux_left/right` shapes + invariants.
- `tests/calibration/test_substrate*.py` — raw pos/neg + sense/antisense
  reductions; boundary projection over per-side flux.
- `test_scanner_accumulator_integration.py` — re-baseline payload keys.
- Rename touch-ups in `test_buffer.py`, `test_pipeline_routing.py`,
  `scoring_helpers.py`.
- Rebuild the extension; full suite failure mode unchanged (post-calibrate
  `NotImplementedError`).

## Acceptance gate

- Builds clean; `ruff` clean; `pytest --collect-only` clean.
- Accumulator/substrate/integration tests green; no new full-suite failure mode.
- A spliced fragment crossing an annotated junction deposits sense/antisense
  correctly and one-sided flux (no false intron flux); unspliced unchanged.
- Magic-number budget: 0 new constants.

## Rollback

Revert the PR; rebuild. Channel semantics + flux revert to shared genome-strand.
