# Calibration gDNA prior — PRODUCTION reference (M2 background + KDE, summed)

**Status:** authoritative reference for what SHIPS on `main` (M2, commit `9582693a` and later), 2026-07-06.
Use this as the baseline when considering how to overcome the mixture problems. What we ship is **not** the
clean two-model mixture — it is the pre-existing layered architecture with the background variance made
deterministic (M2). The unified mixture was explored and refuted; see §7 + the pointers.

## 1. Shape in one line
A **2-pass** per-node solve. An **M2 background** Gaussian prior on `log f_g` runs in **BOTH passes**; a
population **KDE** is trained on pass-1 densities and **ADDED (summed) in pass 2**. Background and KDE are
combined by **addition in log-space = a product of densities (product-of-experts), NOT a mixture (logsumexp)**.

## 2. The two passes (`calibration/calibrate.py`)
```
belief = _sweep(None)                     # PASS 1: strand + M2 background only
train_sub = build_training_substrate(...) # all non-AMBIG region nodes (intergenic/intron/exon), UNIT weight,
gdna_prior = GdnaDensityPrior.fit(...)    #   AMBIG excluded; a single KDE over EVERYTHING (no bg/enriched split)
belief = _sweep(gdna_prior)               # PASS 2: strand + M2 background + KDE
```
Pass 2 runs only if there are ≥`_MIN_KDE_TEACHERS` (10) teachers; else pass-1 belief stands.

## 3. What a node "sees" — the per-node objective ψ (`simplex_logodds`)
The latent is the log-odds `λ = logit(f_g)`; the solver grids λ on `[−L, L]` and reports the posterior MEDIAN
`f_g = σ(λ)` (interpolated at CDF=0.5) + posterior-mean strand fractions. Over the grid,
```
ψ(λ) =  strand Beta-Binomial mixture log-lik            (the ONLY intrinsic per-node signal; count enters
                                                          only as overdispersed Fisher info — count-zero-info)
      + sided spliced-floor lower bound                 (one-sided mature RNA at boundaries)
      + Jeffreys Beta(½,½) strand reference             (strand-observable nodes)
      + global_lp                                       (§4/§5: the M2 background, +KDE in pass 2)
      + gDNA & per-strand RNA imputation messages       (neighbour density messages, per-edge precision)
      + log σ'(λ) = log(f_g(1−f_g))                     (Jacobian: uniform-λ Riemann sum → ∫·df_g)
```
Posterior ∝ exp(ψ); `f_g` = its median. G1 sinks / empty nodes keep their signature-binary init.

## 4. The M2 background (`bp_solver._global_logprior` + `_LogLinearVarMean`) — BOTH passes
A Gaussian on `log f_g`, per node:
```
target = log( ρ·E/M ) clipped to (0,1] ,   ρ = ρ_global (exon/boundary) | ρ_floor (intergenic/intron override)
N      = min( 1 / (var_mean + σ²_g(ρ_global)), 1 )        # precision, capped at ONE pseudo-observation
glob   = −½ · N · (log f_g − target)²
```
- **`ρ_global`** = exposure-pooled gDNA rate over gDNA-clean seeds; **`ρ_floor`** = pooled rate over
  intergenic+intron REGIONS (the depleted-floor override for those nodes); both fit once, pre-solve.
- **`σ²_g(μ) = max(a + b·log μ, 0)`** — the deterministic closed-form log-linear law (`_LogLinearVarMean`,
  WLS on seed-edge Poisson-corrected excess vs edge-midpoint log-density). This REPLACED the bistable
  MonotoneVarMean P-spline (the cross-process nondeterminism root cause). Evaluated at `ρ_global`.
- **The load-bearing property:** σ²_g is *large under capture* (→ N≈0.03, weak → spares enriched exons) and
  *≈0 at no-gDNA* (→ N=1, strong → suppresses false-positive gDNA). This ONE global density-dependence does
  both jobs. It is GLOBAL (evaluated at ρ_global), not per-node.
- Capped at one pseudo-observation (`_GLOBAL_STAB_PREC=1`) so it can never overrule a node's own strand.

## 5. The KDE (`bp_solver._kde_logprior` + `gdna_density_prior`) — PASS 2 only, ADDED to §4
```
log_rho = log f_g + log(M/E)                             # per node, per grid point
kde     = GdnaDensityPrior.logpdf_kernel(log_rho)        # weighted-Gaussian KDE, REAL quadratic tails
jeff    = −log(1 − f_g)                                  # RNA Jeffreys p(ρ_r)∝1/ρ_r pushed into f_g
_kde_logprior = kde + jeff
```
- Trained on the pass-1 solved densities of all non-AMBIG region nodes, **unit weight**, AMBIG excluded,
  Silverman bandwidth floored at the per-node Poisson noise (continuous `_weighted_median`, deterministic).
- **Real Gaussian tails** (`logpdf_kernel`, NOT the clamped `logpdf`) — a genuine `−½((x−xᵢ)/h)²` tail
  penalizes implausibly-high ρ_g (removes the FP the clamped tail caused: ~585k → ~2k).
- The **Jeffreys** term is load-bearing: without it the cheap explanation of any flat-strand node is "dump mass
  into free RNA, park gDNA at the tall depleted KDE mode" (the cliff + FP). With it, lowering f_g raises the
  penalized ρ_r, so gDNA is the residual after a typical RNA.

## 6. How M2 + KDE combine — SUM, not mixture
`global_lp = _global_logprior + _kde_logprior` (pass 2). This is a **product of densities**: both must be
satisfied. M2 (weak, capped) anchors + suppresses no-gDNA FP + spares enriched via σ²_g; the KDE (uncapped,
real tails) adds the empirical population structure (the enriched mode under capture) + extra FP suppression.
Because M2 is capped weak, the KDE dominates where it is informative while M2 keeps the solve anchored/finite.
This is the "floor + KDE" layering — NOT the clean-slate "gravity to the nearest of a background mode and an
enriched mode" (that would be a logsumexp mixture).

## 7. Are we unified? No — and why the mixture didn't win (for §-considering-fixes)
| | unified two-model (aspiration) | PRODUCTION (M2 + KDE) |
|---|---|---|
| background | its own model, one mode | M2 Gaussian (ρ_global / ρ_floor) |
| enriched | its own KDE (exon+boundary) | KDE over ALL nodes |
| combine | MIXTURE (logsumexp, nearest-mode gravity) | SUM (product of densities) |
| σ²_g | eliminated | KEPT (M2's deterministic line) |

Refuted mixture attempts (do not repeat; details in the linked docs + memory):
- **Seeded mixture / drop the genome baseline** (`unified_two_model_prior_WIP.md`): the background seed fixes the
  no-gDNA FP (validated) but σ²_g's *consistent unstranded-exon anchoring* (ranks) is irreducibly load-bearing;
  net −0.049 spearman.
- **Distance-redescending background** (`unified_prior_redescending_design.md` §9): breaks no-gDNA FP — for an
  UNSTRANDED node the prior alone sets f_g and a redescending prior is weak far out, so no-gDNA exons drift up.
- **The real wall = IDENTIFIABILITY on unstranded capture:** an enriched exon's high M/E is gDNA-vs-RNA
  ambiguous with no strand. No prior shape invents the missing signal. The win must come from a **non-prior
  signal** — structural boundary/spliced RNA residual, or mappability/GC covariates — NOT a cleverer prior.

**Why M2 beats every mixture we tried:** its GLOBAL density-dependent precision is *strong-everywhere at
no-gDNA* (suppress FP) and *weak-everywhere under capture* (spare enriched). A per-node distance/uncertainty
scale cannot be "strong-everywhere at no-gDNA" — that is the property the mixture keeps losing.

## 8. Determinism + tests
Every prior input is a continuous data function (no argmin/GCV/IRLS; interpolated median). Cross-process
2.6% → ~0.05% (the M2 change; residual is the C++ scanner's parallel FP — a separate follow-up). Full suite
green; goldens regenerated on the M2 ship.
