# ARCHIVED — the joint count×strand deconvolution

**Status:** archived 2026-06-10. Retired in favour of the decoupled architecture
(`../decoupled_calibration_design.md`). This document is the **resurrection guide**: enough detail to
rebuild the joint approach if the future re-coupling step (decoupled doc §8) proves worthwhile. The
code is not kept in the live tree — recover it from git.

## Last intact commit

- **`f5da375`** — the joint deconvolution + the count overdispersion model are both present and wired.
  `git show f5da375:src/rigel/calibration/joint_deconv.py`,
  `git show f5da375:src/rigel/calibration/count_dispersion.py`,
  `git show f5da375:src/rigel/calibration/density_model.py`,
  `git show f5da375:src/rigel/calibration/calibrate.py`.
- (The intermediate "Phase 0/1" working state — overdispersion torn down, own-mean floor added — was
  never committed; it is superseded and not worth preserving.)

## What it was

A per-node Bayesian posterior over the gDNA fraction `g`, combining two conditionally-independent
clues as a **product**:

```
posterior(g) ∝ Beta(g ; a_c, b_c)  ×  BetaBinomial(sense | n, p(g), overdispersions)
p(g) = g·½ + (1−g)·κ                              # node sense rate
```

- **Count prior** `Beta(a_c, b_c)`: mean `count_gdna_frac = clip(density·eff_len / mass, 0, 1)` where
  `density` was the per-region gDNA density from `density_model` — **strand-cleaned** by the linear
  unmix `ĝ = (s−κ)/(½−κ)` (`density_model.strand_clean_gdna_frac`). Concentration = the
  **overdispersion-limited effective count** `N_eff = N/(1+α·N)` from `count_dispersion.py` (NB MoM
  per count-type, contained vs crossing), Jeffreys-floored `Beta(½,½)` at zero evidence.
- **Strand likelihood**: `strand_likelihood.strand_loglik` — a two-component (gDNA at ½, RNA at κ)
  Beta-Binomial with fitted gDNA and RNA strand overdispersions (`gdna_strand.py`).
- Reported `g` = posterior **median** on a grid, optionally log-odds-shifted by `gdna_strand_llr_bias`.
- gDNA mass = `g·M`, RNA = `(1−g)·M` + spliced. Mass conserved per node. Lived in
  `joint_deconv._joint_per_node`, applied via `deconv_regions` (contained) and `deconv_sides`
  (boundary sides, `_compute_side`).

Strand was thus used **twice**: to strand-clean the count-prior mean *and* in the strand likelihood
(the "double-use").

## Why it was retired

The 20-scenario benchmark showed the product mis-weights two **unequal** estimators:

- The strand estimator is *unbiased* (noisy at low N/SS); the count estimator is *biased* under
  hybrid capture (exon density imputed from depleted off-target neighbours, ~2× low).
- The count overdispersion `α` exploded under capture (between-region mean heterogeneity booked as
  dispersion: α_crossing 0.0005→86), crushing `N_eff`→~0. This **accidentally** helped — it nullified
  the biased count so the excellent strand could win.
- Removing the overdispersion (to "honour the count") let the biased count prior fight the strand and
  **tripled** the flagship leak (gdna1000 cap_on ss0.99: 3.74% → 10.95%). Titrating the count's
  precision to defer to strand at high SS proved intractable.

Conclusion: don't mix a biased estimator into an unbiased one. Decouple (use strand where it has
signal, count only where strand is absent). See `../decoupled_calibration_design.md` §1.

## When to resurrect / how to do it right

The joint form is *correct in principle* — it is a precision weighting — and would beat the decoupled
handoff in the genuine synergy regimes (low-count nodes, weak-SS libraries). Resurrect **only** once
the count module is made *honest*: an accurate mean (the point-5 unspliced-fraction imputation,
`../count_channel_capture_design.md`) and a **calibrated precision** that reflects imputation
uncertainty rather than the raw count. With an honest count precision, the product weights correctly
on its own (no overdispersion crutch). Re-introduce as the final, optional step of the decoupled plan
(§8), gated to fire only at low strand information, and re-validate on the full benchmark + the nrna
nascent scenarios.
