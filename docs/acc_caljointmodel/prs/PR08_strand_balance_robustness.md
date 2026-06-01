# PR 8 — Strand-balance robustness: orientation, perfect-stranding degeneracy, boundary double-count

A separate calibration bug fix surfaced by the zero-gDNA deep dive
([calibration_zero_gdna_diagnosis.md](../calibration_zero_gdna_diagnosis.md)).
Independent of PR 7 (eps_s removal); both improve robustness before validation.

`κ_rna` (RNA sense mean) comes from `StrandModel.p_r1_sense`; `ρ_r_bb` (RNA strand
overdispersion) is fit by `fit_strand_balance` from the substrate's spliced
channels. On the single_exon baseline (ss=1.0, 39 spliced reads) the strand model
is degenerate in three distinct ways.

## Finding 1 — `p_r1_sense` is almost certainly **inverted** (the clear bug)

Observed: `p_r1_sense = 0` exactly (all 39 spliced reads measured antisense),
`strand_specificity = 1.0`.

Expected: the simulator (`sim/reads.py`) uses an **R1-antisense** library prep,
**but documents that the BAM scanner flips R2's strand so both mates report the
correct transcript strand** ("R1-antisense convention … the BAM scanner flips
R2's strand → both report + → correct"). After that flip, spliced reads should
align *sense* to their transcript ⇒ `p_r1_sense ≈ 1`, not 0.

So `p_r1_sense = 0` is strong evidence that the **spliced-read strand is measured
inverted** (sense/antisense swapped) somewhere between the scanner's SJ-strand
assignment and `StrandModel`'s `pos_pos/neg_neg` accumulation. It is internally
self-consistent enough that single_exon's *contained* reads still resolve to RNA
(they're "antisense" in the same flipped frame), which is why it didn't surface
until now — but an inverted `κ_rna` will mis-call on mixed-strand / gDNA-bearing
scenarios.

**Why it matters:** `κ_rna` feeds the E-step strand log-Bayes-factor
(`RNA BB(κ_rna, ρ_r_bb)` vs `gDNA BB(0.5, ρ_d_bb)`). If it's the wrong sense mean,
every strand call is reasoning from a backwards RNA model.

## Finding 2 — perfect stranding collapses the strand channel to a knife-edge

`κ_rna = p_r1_sense` is clamped to `[_BB_FLOOR, 1 − _BB_FLOOR] = [1e-6, 1−1e-6]`.
A perfectly stranded library (`p_r1_sense ∈ {0,1}`) therefore yields a
Beta-Binomial with shape `(a≈1, b≈1e6)` (or its mirror) — razor-sharp: **any**
read even slightly off the expected strand scores as gDNA. And `ρ_r_bb` *also*
floors at `1e-6` here (MoM returns ≈0 when `κ` is at a bound), removing the
overdispersion that is supposed to soften the channel. A perfectly stranded
library is the *easy* case; it should not make the strand discriminator
pathologically brittle.

## Finding 3 — boundary spliced reads are double-counted

`fit_strand_balance` pooled `n_reads = 124` from only **39** spliced fragments. A
junction read crosses two boundaries, so it appears in *both* flanking regions'
boundary views (exon1's right view *and* exon2's left view) and is pooled twice
(more for multi-junction reads). This inflates the observation set feeding the
`ρ_r_bb` method-of-moments fit and biases it.

## Investigation plan (resolve before finalizing fixes)

1. **Confirm the inversion.** Trace a handful of spliced fragments end-to-end:
   transcript strand, the scanner's SJ-motif strand + R2 flip, the resolved
   align strand, and which `StrandModel` counter (`pos_pos/neg_neg` vs
   `pos_neg/neg_pos`) each increments. Determine where sense/antisense flips.
   Cross-check a `−`-strand gene and a real (non-oracle) BAM if available.
2. **Quantify the double-count** and confirm each junction fragment is pooled
   once per boundary side it touches.

## Candidate fixes (refine after investigation)

- **Orientation (Finding 1):** correct the spliced-read sense/antisense
  accounting at its source (scanner SJ-strand vs `StrandModel` convention) so
  `p_r1_sense` reflects the true sense fraction. (Likely a one-line convention
  fix once located — but it must be traced, not guessed.)
- **Perfect-stranding robustness (Finding 2):** keep the strand channel
  well-conditioned at the bounds. Options to weigh (some introduce a constant →
  Q6 discussion): a less-extreme `κ_rna` floor; a small **minimum effective
  `ρ_r_bb`** so the BB is never a knife-edge; or a cap on the per-region strand
  LLR. Prefer the overdispersion floor (it's the model's existing robustness
  knob) if a principled value exists.
- **Double-count (Finding 3):** de-duplicate spliced observations in the pool —
  count each junction fragment once (e.g., attribute each junction to a single
  boundary side, or weight by 1/sides), consistent with how the accumulator
  emits per-side flux.

## Open questions (for critique)

- Is Finding 1 a true inversion or a self-consistent convention we can leave?
  (Investigation decides; my read is it's an inversion given the sim's documented
  scanner-flip intent.)
- Finding 2: which robustness mechanism, and what value — is there a principled
  `ρ_r_bb` floor, or should `κ_rna` simply not be pushed to the `_BB_FLOOR`
  extreme?
- Does the double-count change any conclusion once `κ_rna` is corrected (here
  `ρ_r_bb` is degenerate anyway because `κ_rna` is at a bound)?

## Relationship to PR 7

Independent. PR 7 removes the spliced→gDNA `eps_s` phantom; PR 8 fixes how the
strand channel reasons about (genuine) spliced reads and behaves at perfect
stranding. Either can land first; both are needed for robust low-/mixed-gDNA
behaviour ahead of larger validation.

## Rollback

Each fix is independently revertable (orientation convention, `ρ_r_bb` floor,
pooling de-dup).
