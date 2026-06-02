# PR 10 — Spliced junction reads are double-counted in the ρ_r_bb pool

> **SUBSUMED BY PR 9 (no separate PR needed).** PR 9 replaced the method-of-moments
> ρ_r_bb fit (which consumed the double-counted substrate spliced pool) with the
> posterior-predictive Beta-Binomial drawn from the `StrandModel`'s **2×2 fragment
> counts** — which are per-fragment from the scanner and never double-counted. With
> the substrate pool no longer used for the strand fit, the double-count described
> below cannot occur. This doc is retained as the rationale for *why* the MoM pool
> was unsafe. No code action remains.

Third of the 3-PR robustness split (PR8 AMBIG → PR9 strand knife-edge → **PR10
spliced double-count**). Smallest of the three; a correctness fix to the
overdispersion fit.

## Summary

`fit_strand_balance` pools its `(k_sense, n)` spliced observations from all three
substrate views (contained + left + right). A spliced junction read crosses a
boundary, so it is deposited on **both** sides of that boundary and therefore
appears in **both** flanking regions' boundary views — pooled twice (more for
multi-junction reads). On single_exon this pooled 124 reads from only **39**
spliced fragments. The inflated, duplicated observation set biases the
method-of-moments `ρ_r_bb` fit.

## Root cause (confirmed)

The accumulator deposits a boundary-crossing fragment's flux on both sides of the
boundary (`boundary_flux_left[B]` and `boundary_flux_right[B]`). The substrate maps
those to the flanking regions
([substrate.py:113](../../../src/rigel/calibration/substrate.py#L113)):

- region `r`'s **right** view ← `boundary_flux_left[B]` (r is B's left region)
- region `r+1`'s **left** view ← `boundary_flux_right[B]` (r+1 is B's right region)

`fit_strand_balance`
([strand_balance.py:92](../../../src/rigel/calibration/strand_balance.py#L92))
then loops `for view in (contained, left, right)` and concatenates all non-empty
`(n_spliced_sense, n_spliced)` rows. So the single fragment crossing `B` is counted
once via region `r`'s right view and again via region `r+1`'s left view.

This is *only* a problem for the strand pool: the per-side split is correct and
necessary for **mass** attribution (D1), where each side legitimately carries its
share. But the `ρ_r_bb` MoM wants each junction *observation* once.

**Evidence:** single_exon `fit_strand_balance → n_obs=4, n_reads=124` from 39
spliced fragments (the strand-model 2×2 counts the same 39 once; the substrate pool
double-counts).

## The fix

Pool each spliced junction fragment **once** for the `ρ_r_bb` fit. Since contained
spliced flux is ~0 by construction (a spliced fragment spans a boundary, so it is a
boundary observation, not a contained one), the strand pool is effectively the
boundary flux — and each boundary's flux is the *same* fragments on both sides.

Preferred: pool **each boundary once** rather than both flanking region views — e.g.
read the boundary spliced flux directly (one entry per internal boundary) instead of
concatenating `substrate.left` and `substrate.right`. This keeps integer
`(k_sense, n)` observations (required by the BB MoM) and counts each crossing
fragment exactly once.

Constant-free (Q6-clean) — it changes *which* rows are pooled, adds nothing.

## Alternatives considered

- **Weight each side by 1/2.** Restores the count but breaks the integer
  `(k_sense, n)` the Beta-Binomial MoM consumes. Rejected.
- **Use only the contained view.** Contained spliced ≈ 0 → almost no observations →
  `ρ_r_bb` falls back. Rejected (throws away the real signal).
- **De-dup by fragment id.** The substrate is aggregate flux, not per-fragment, so
  there is no id to de-dup on. The boundary-level fix is the aggregate-correct
  equivalent.

## Test plan

- Unit: `fit_strand_balance` on a known junction layout reports `n_total_reads`
  equal to the true spliced fragment count (no 2×/3× inflation), and `ρ_r_bb`
  matches the hand-computed MoM on the de-duplicated observations.
- Confirm the corrected `ρ_r_bb` does not regress the calibration on the scenario
  suite (it should be neutral-to-better). Regenerate goldens if `ρ_r_bb` shifts.

## Interaction with PR8 / PR9

Smallest and most localized. Orthogonal to PR8. Mild interaction with PR9: both
touch `ρ_r_bb`, but PR9's `ε_CI` floor is computed from the strand-model 2×2 counts
(already de-duplicated), so PR10 mainly cleans the MoM *point estimate* while PR9
sets its *floor*. Recommend landing PR10 last (after PR8, PR9) as a clean-up, or
independently — no ordering constraint.

## Rollback

Single revert point: restore the three-view concatenation in `fit_strand_balance`.

## Open questions (for critique)

1. Which boundary side to read when pooling once — left, right, or the average?
   (They carry the same fragments; the motif-relative `n_spliced_sense` should be
   identical on both sides — confirm in implementation, then pick one canonically.)
2. Do multi-junction reads (crossing >2 boundaries) need special handling, or does
   "one observation per boundary" already count them correctly as one observation
   per junction? (My read: per-junction is the right granularity for a per-junction
   strand model; confirm against the accumulator's deposit semantics.)
