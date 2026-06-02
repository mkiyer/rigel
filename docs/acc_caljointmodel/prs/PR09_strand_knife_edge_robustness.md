# PR 9 — Strand knife-edge: keep the RNA Beta-Binomial well-conditioned at (near-)perfect stranding

Second of the 3-PR robustness split (PR8 AMBIG → **PR9 strand knife-edge** → PR10
spliced double-count). Independent of PR8; either can land first.

## Summary

At perfect or near-perfect strand specificity the RNA strand Beta-Binomial
collapses to a razor edge: the sense mean `κ_rna` is pushed to its clamp bound and
the overdispersion `ρ_r_bb` floors to `1e-6`, giving `BB(α≈1, β≈1e6)`. Then **any**
read even slightly off the expected strand scores as gDNA. A perfectly stranded
library is the *easy* case; it should not make the discriminator pathologically
brittle. This drives RNA→gDNA leakage on the imperfect-`ss` and mixed scenarios.

## Root cause (confirmed)

`fit_strand_balance` ([strand_balance.py](../../../src/rigel/calibration/strand_balance.py)):
- `κ_rna = clamp(p_r1_sense, _BB_FLOOR, 1−_BB_FLOOR)` with `_BB_FLOOR = 1e-6`. A
  perfectly stranded library has `p_r1_sense ∈ {0,1}` ⇒ `κ_rna` sits **at the bound**.
- `ρ_r_bb` is method-of-moments; with no strand spread (perfect `ss`) MoM returns
  ≈0 and is floored to `_BB_FLOOR = 1e-6` — removing the very overdispersion that is
  supposed to soften the channel.

Result on single_exon (ss=1.0): `κ_rna=1e-6, ρ_r_bb=1e-6` →
`a_r = κ(1−ρ)/ρ ≈ 1`, `b_r = (1−κ)(1−ρ)/ρ ≈ 1e6` → the RNA BB puts essentially all
mass on a single sense count. In the E-step strand log-Bayes-factor
([estep.py:89](../../../src/rigel/calibration/estep.py#L89)) a region whose sense
fraction is off by even one read gets `ll_d → −∞` relative to the symmetric gDNA
`BB(0.5, ρ_d)` ⇒ scored gDNA.

(Note: the `κ_rna≈0`-vs-`≈1` orientation is the **self-consistent convention**
established in the PR8 diagnosis — both `κ_rna` and the unspliced `k_sense` live in
the same R1-orientation frame, so the channel still separates RNA from gDNA. The
bug here is the *sharpness*, not the orientation. The only orientation action is to
fix the misleading "both report the correct transcript strand" docstring in
[strand_model.py](../../../src/rigel/strand_model.py).)

## The fix (principled — reuse the existing measurement-uncertainty floor)

`StrandModel` already exposes exactly the right quantity but it is **not wired into
the calibration strand channel** (it is only referenced for a pipeline warning):

> `strand_specificity_ci_epsilon(confidence)` →
> `ε_CI = 1 − UCL(strand_specificity)`, the one-sided upper credible limit on the
> minor-orientation rate `1−ss` under a `Beta(k_minor+1, k_major+1)` posterior
> ([strand_model.py:224](../../../src/rigel/strand_model.py#L224)). Its own
> docstring says it exists "to avoid runaway strand LLRs when the training set is
> small or when ss_est saturates to 1.0."

`ε_CI` is the **statistically justified minimum off-strand rate** consistent with
the finite spliced sample — not a magic number. Use it to keep the RNA BB off the
knife edge. Two coupled levers (decide which, or both, in critique):

- **Mean floor:** clamp `κ_rna` into `[ε_CI, 1−ε_CI]` instead of
  `[_BB_FLOOR, 1−_BB_FLOOR]`. At perfect stranding with finite `n`, `ε_CI > 0`
  pulls `κ_rna` off the hard bound so the BB expects an `ε_CI` minority off-strand.
- **Dispersion floor:** floor `ρ_r_bb` at the overdispersion implied by `ε_CI`
  (the minor-rate posterior's spread) rather than `_BB_FLOOR`. This restores the
  softening the channel is designed to have.

Both derive entirely from the spliced sample size — **no new constant** (Q6-clean):
they replace a fixed `1e-6` floor with the data's own credible interval. As `n → ∞`
with a truly perfect library, `ε_CI → 0` and we recover the sharp channel; with the
small spliced sets these toy scenarios have, `ε_CI` is appreciable and the channel
stays tolerant.

## Alternatives considered

- **A fixed `ρ_r_bb` floor (e.g. 0.01/0.05).** Simple, but a magic number with no
  principled value (Q6 violation). The `ε_CI` route is strictly better — same
  effect, derived from data.
- **Cap the per-region strand LLR at ±C.** Treats the symptom; needs a constant.
- **Leave `κ_rna` clamp, only fix `ρ_r_bb`.** May be insufficient — a razor BB
  centred at the bound still rejects the `ε_CI` minority. Prefer the coupled fix.

## Test plan

- The imperfect-`ss` scenarios (`overlapping_antisense` strand sweeps, the `ss<1`
  arms of the sweeps): RNA no longer leaks to gDNA at `ss∈{0.9,0.65}`.
- Perfect-`ss` regression: single_exon / non_overlapping still resolve cleanly
  (the channel must not become *too* soft and admit gDNA as RNA).
- Unit tests: `fit_strand_balance` returns `κ_rna`, `ρ_r_bb` off the hard bound when
  `n` is small/perfectly stranded; recovers the sharp values as `n → ∞`.
- A gDNA-bearing stranded scenario: confirm gDNA is still rejected from RNA (the
  soften must not blunt real discrimination). Regenerate goldens.

## Interaction with PR8

Complementary. PR8 stops the `ρ₀` runaway (the dominant failure); PR9 hardens the
strand channel so RNA at imperfect `ss` is not lost to gDNA. Some PR8 scenarios run
at ss=1.0 and are fixed by PR8 alone; the `nrna_dc g20_n70_s65`-type low-`ss`
leakage needs PR9. Land order is free.

## Rollback

Single revert point: restore the `_BB_FLOOR` clamps in `fit_strand_balance` and
unwire `ε_CI`.

## Open questions (for critique)

1. Mean floor, dispersion floor, or both? (My read: both — they address the two
   independent degeneracies; either alone leaves a residual edge.)
2. `ε_CI` confidence level — `strand_specificity_ci_epsilon` defaults to 0.99.
   Is 0.99 right here, or should the calibration pick its own? (This *is* a
   parameter, but a meaningful statistical one — worth an explicit decision.)
3. Does the double-count (PR10) materially bias `ρ_r_bb` such that PR9 should land
   after PR10? (My read: PR10 changes the `ρ_r_bb` MoM input but PR9's `ε_CI` floor
   is computed from the strand-model counts, not the substrate pool, so they're
   largely orthogonal.)
