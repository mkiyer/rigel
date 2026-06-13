# Count local-dispersion estimator (PARKED idea)

**Status:** parked design note. 2026-06-10. Not in the v1 plan — recorded so the idea isn't lost.
Phase 4-var (`count_posterior_design.md`) ships a Poisson counting baseline + the variance~mean fit
from paired boundary disagreements; this note keeps the *local count-dispersion* estimator for if/when
the counting-noise term needs to be more than Poisson.

## The idea

RNA-seq counts are NB-overdispersed (`Var = μ + α·μ²`). A **global** `α` is exactly what exploded under
hybrid capture — it booked on/off-target *mean* heterogeneity as dispersion (the count-overdispersion
teardown that triggered the decoupling). The fix is to estimate `α` **locally in count space** so
on-target (high count) and off-target (low count) seeds get their own dispersion:

- **kNN:** for a query count `x`, take the `k` nearest observable-seed counts; `α̂(x)` from their
  excess variance. Knob: `k`.
- **Gaussian kernel over sorted counts (preferred):** convolve a Gaussian (bandwidth `h` in **log-count**
  space) across the sorted seed counts; `α̂(x)` is the proximity-weighted local excess variance —
  smooth, no hard bins, near neighbours weighted more. Knob: bandwidth `h`. A 1-D convolution producing
  a function `α(count)`.

Estimator per kernel window around expected count `μ`:
```
α̂(μ) = [ Σᵢ wᵢ ((Nᵢ − μ)² − Nᵢ) ] / [ Σᵢ wᵢ μ² ],   wᵢ = Gaussian(log Nᵢ − log μ; h)
```
over observable gDNA-pure seed nodes.

## Why it's parked (the caveat)

The kernel pools **unpaired** counts from *different* seed nodes that do **not** share a true rate
(capture gives each its own), so `(Nᵢ−μ)²` still conflates within-window rate-variation with sampling
dispersion — the same conflation as the global `α`, only localized and reduced. It is an *upper bound*
on true dispersion, and the bandwidth `h` trades conflation against noise.

**Phase 4-var's variance~mean fit avoids this** by pooling **paired** same-node anchor disagreements
(`d_L` vs `d_R` of one node — one true rate, so the disagreement is pure variance). That's the
better-founded path, and it already gave a clean `var ∝ mean²` on the benchmark. So this local-count
estimator is only worth revisiting if the *counting* (not imputation) term needs to exceed Poisson —
e.g. for observable nodes that route to count, or to sharpen the propagated anchor count-noise term.

## If revived

- Use the Gaussian kernel in log-count space (bandwidth a data-driven fraction of the log-count IQR).
- Validate against the variance~mean fit: the two should agree where both apply; if they diverge,
  trust the paired (variance~mean) estimate.
- Keep it opt-in behind a config flag; default Poisson.
