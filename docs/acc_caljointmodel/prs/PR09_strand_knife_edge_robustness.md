# PR 9 — Strand Beta-Binomial: posterior-predictive overdispersion (subsumes PR 10)

Second of the post-PR07 robustness work. Replaces the RNA strand BB's
arbitrary-floor overdispersion with a **posterior-predictive** fit driven by
spliced-read statistical power — removing three magic numbers and folding in PR 10.

## Investigation first: the knife-edge is NOT what fails the nrna_dc cases

We suspected the strand "knife-edge" (`ρ_r_bb` floored to `1e-6` ≡ a Binomial at
`κ_rna`) caused the 3 remaining `nrna_dc g20` failures. **Disproven empirically:**
softening the strand BB (dispersion floor, mean floor, both — across the full range
of candidate values) changes those scenarios' metrics by ≤4 fragments. The reason:
`π_g = expit(logit(prior) + llr_count + llr_strand)` — the **count channel + prior
already saturate π_g**, so the strand LLR is irrelevant there. The 3 `nrna_dc`
failures are count/density/EM issues (a separate investigation), not strand.

So PR 9 is **not** about fixing those scenarios. It is a *principled-robustness +
magic-number removal* step (your framing): real datasets with **sparse or zero
spliced RNA** exist and must degrade gracefully, and the `1e-6` floor is an
arbitrary constant we can replace with statistical power.

## What ρ_r_bb is, and what was wrong

`κ_rna` (RNA sense mean) and `ρ_r_bb` (RNA strand BB overdispersion) are a **fixed
one-time fit** (not EM-trained) — `fit_strand_balance` computes them once and the
M-step never touches them. Previously:
- `κ_rna` = `StrandModel.p_r1_sense` (clean per-fragment 2×2), clamped to `[1e-6, 1−1e-6]`;
- `ρ_r_bb` = method-of-moments over the substrate's per-region spliced pool, floored
  at `_BB_FLOOR=1e-6`, with a `_RHO_R_BB_FALLBACK=0.01` for `<_MIN_OBS_FOR_OVERDISPERSION=2`.

Three problems: (a) the `1e-6` floor is arbitrary; (b) the MoM pool **double-counts
boundary spliced reads** (a junction read lands in both flanking views — this was
PR 10); (c) with 1–2 spliced reads a point-estimate κ feeds a razor-sharp Binomial.

## The fix: posterior-predictive Beta-Binomial

The RNA strand BB's Beta prior **is** the κ_rna posterior `Beta(n_same+1, n_opp+1)`
from the StrandModel's 2×2 fragment counts (Laplace prior `Beta(1,1)`). So:

- `κ_rna = (n_same+1)/(n_obs+2)` — posterior mean (never at 0/1 → no clamp needed);
- `ρ_r_bb = 1/(n_obs+3)` — the posterior-predictive overdispersion, a pure function
  of spliced-read count (statistical power).

The Beta-Binomial collapses to the **Binomial** as `ρ_r_bb → 0` (confirmed:
`max|BB − Binomial| = 2.5e-6` at `1e-6`), so:

```
n_obs = 200,000 → ρ ≈ 0       → Binomial (correct sharp limit, real data)
n_obs = 200     → ρ ≈ 0.005   → ~Binomial
n_obs = 5       → ρ ≈ 0.125   → wide (robust)
n_obs = 0       → Beta(1,1)   → κ=0.5 symmetric → strand channel neutral
```

This **deletes** `_mom_overdispersion`, `_RHO_R_BB_FALLBACK`,
`_MIN_OBS_FOR_OVERDISPERSION`, and the `_BB_FLOOR` clamp on the RNA params
(`_BB_FLOOR` stays only for the gDNA `ρ_d_bb` bound, moved to mstep). And because the
fit now draws from the StrandModel's per-fragment 2×2 — not the substrate pool — the
boundary double-count disappears: **PR 10 is subsumed.** `fit_strand_balance` no
longer takes the substrate.

## Validation

- `nrna_dc g20` cases **unchanged** (dense spliced → ρ≈0.005 → ~Binomial → as before;
  count/EM still dominates them — PR 9 was never going to fix them).
- single_exon false gDNA **0.69 → 8e-5** (the softer sparse-data BB stops mis-flagging
  the clean RNA) — a real improvement.
- 115 calibration units + 21 golden tests pass; ruff clean. Goldens regenerated
  (global strand-param change, like PR 7): the κ_rna posterior mean shifts per-region
  splits across ~all scenarios (1–5%); sanity-checked, no garbage.
- Edge-case behaviour proven by unit tests (sparse → wide, dense → Binomial, n_obs=0
  → symmetric); a **scenario** pool for zero-RNA / sparse-spliced / zero-gDNA is the
  next (separate) PR.

## Rollback

Single revert point: restore the MoM `fit_strand_balance` + the three constants.
