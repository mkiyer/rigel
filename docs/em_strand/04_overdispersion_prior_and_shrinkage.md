# Design — replace the overdispersion-fit thresholds with principled shrinkage

**Status: IMPLEMENTED (2026-06-07).** The two Phase-2 gate constants (`_MIN_SEED_NODES = 30`,
`_SIGNIFICANCE_Z = 2`) are **removed**, replaced by precision-weighted (seed-node-count)
shrinkage toward a prior overdispersion. Final design (per the user): the prior is parametrized
as a **symmetric `Beta(a, a)`** with `od₀ = 1/(2a + 1)` — config `gdna_strand_prior_alpha_beta`
(default **3** ⇒ `od₀ ≈ 0.143`, a conservative overdispersed "floor" that works well for sparse
data) — with strength `gdna_strand_prior_weight` (effective seed-node units, default 30). The
fit is clamped to the **`Beta(2, 2)` ceiling `od ≤ 0.2`** ("the most overdispersion we allow").
Sparse/no data ⇒ `od₀`; abundant data ⇒ the fitted MoM. (Implementation note: shrinkage is
weighted by *seed-node count* rather than the raw MoM SE of §3A — the SE underestimates
uncertainty at low node counts, the very regime the prior must dominate; node-count weighting is
the robust choice that retired the hard min-node gate.) The §2 reasoning (default overdispersed,
not Binomial) and the conversion are unchanged.

## 1. The problem with the current gate

Phase 2 trusts a positive gDNA strand overdispersion only above **two arbitrary procedural
cutoffs**: ≥ 30 seed nodes *and* a 2σ excess-variance significance test; below either, od snaps
to 0 (Binomial). These are magic numbers in the worst sense — no physical meaning, a
discontinuous cliff (29 nodes → Binomial, 30 → full MLE), and they silently encode a *Binomial
default* on thin data, which is precisely the assumption this whole effort exists to remove.

**Principle:** estimation should degrade **continuously** with information — the data's signal
should determine how far we move from a sensible default toward the library-specific estimate.
That is regularization / Bayesian shrinkage, and it needs **no thresholds**.

## 2. The two questions, answered

### 2.1 Default: Binomial or overdispersed? → **Overdispersed.**

The default (what we believe *before* this library's seeds inform us) should be a **mild
positive overdispersion `od_0 > 0`**, not Binomial:

1. **We know gDNA strand is overdispersed.** Real RNA-seq gDNA shows substantial
   region-to-region strand skew (the motivating observation). Defaulting to Binomial on weak
   data **reinstates the exact bug** Phase 1–2 fixed.
2. **It is the FP-safe direction.** Overdispersion widens the gDNA strand distribution → a
   skewed node reads as gDNA rather than RNA → fewer RNA false positives. Erring toward `od_0`
   is conservative in the direction we care about.
3. **`od = 0` is a degenerate boundary.** The Beta concentration `a = ½(1−od)/od → ∞` there; a
   default sitting on the boundary is awkward and asymmetric (the data can only push *up*).

So the default is *overdispersed*, and **Binomial becomes a data-driven outcome** (a library
whose seeds clearly show no overdispersion shrinks `od` down toward 0) rather than the
thin-data fallback.

### 2.2 Shrink to what? → **A biological prior mean `od_0`, precision-weighted by the data.**

Shrink the library MoM estimate toward `od_0` — the *typical gDNA strand intra-class
correlation*, a single physical quantity we can elicit from real-data studies — with the weight
set by the data's own information content. Abundant / strong-signal seeds → the library MLE;
thin / weak seeds → `od_0`. Continuous, no cliff.

## 3. Options

**A. Precision-weighted shrinkage (Normal–Normal conjugate) — recommended.**
The MoM already yields a point estimate `od_mom = num/den` *and* a sampling SE
`se = √(Σ 2·(N·μ(1−μ))²)/den` (the same quantity the current gate thresholds on — repurpose it
as a **weight**, not a cutoff). With a prior `od ~ Normal(od_0, σ0²)`:

```
od_post = (od_mom/se² + od_0/σ0²) / (1/se² + 1/σ0²),  clamped to [0, 1)
```

Thin/weak data ⇒ `se` large ⇒ `od_post → od_0`. Abundant data ⇒ `se` small ⇒ `od_post → od_mom`.
Closed-form, O(n_seed_nodes), **removes both magic numbers**. Introduces one biological prior
`od_0` (+ a width `σ0`, see §5). No threshold anywhere.

**B. Beta-Binomial MAP with a prior on od.** Replace the MoM with the exact Beta-Binomial
log-likelihood over seeds + a prior on `od`; report the posterior mode/mean (1-D optimize).
More faithful at low counts than the MoM's Gaussian SE; same shrinkage philosophy; modestly
more compute. A good upgrade if §A's normal approximation proves too coarse on real data.

**C. Empirical Bayes — estimate `od_0` from data.** With a *corpus* of libraries (or grouped
nodes within one), estimate the prior mean/width from the spread of per-group estimates, then
shrink each toward it (James–Stein). Fully data-driven — *eliminates the fixed `od_0`* when
there is enough data; degenerates to needing a prior in the thin limit (unavoidable). Best
long-term once we maintain a library corpus.

**D. Marginalize od.** Carry the full `od` posterior into the decode (integrate the strand
likelihood over `od`) instead of a point estimate. Most principled, most invasive — likely
overkill given the effect size; noted for completeness.

## 4. Is `od_0` just another magic number? — No.

Categorically different from `30` and `2σ`. `od_0` is a **physical quantity** — the typical
gDNA strand intra-class correlation — with units, meaning, and an empirical referent we have
*measured* and can re-measure. `30` and `2σ` are arbitrary procedural cutoffs with no physical
content. A prior grounded in measurement is the opposite of a magic number; it is a calibration
constant. **Action item:** set `od_0` from the real-data gDNA strand-balance studies (please
supply the measured value); until then we use a documented placeholder and flag it.

## 5. Precise "default" semantics + the prior width

- No seed data → `od_post = od_0`.
- Weak/noisy seeds → shrunk toward `od_0`.
- Abundant, clearly-overdispersed seeds → the library-specific MLE.
- Abundant, clearly-Binomial seeds → `od_post → 0` (the data overrides the prior *downward* too
  — Binomial is reachable, just not the default).

The width `σ0` controls how quickly data overrides the prior. To avoid a *second* free number,
tie it to `od_0`: **`σ0 = od_0`** (a weakly-informative CV = 100% prior — the prior says "about
`od_0`, give or take itself," so even modest data dominates). This keeps a single biological
input `od_0`; `σ0` is derived, not chosen. (Or express the prior strength in pseudo-node units
`k0`: the prior is worth `k0` seed nodes of evidence, `k0` ~ a few. Equivalent; pick whichever
reads cleaner in code.)

## 6. Recommendation

Adopt **Option A** now: precision-weighted shrinkage toward a biological `od_0` (default
overdispersed), `σ0 = od_0`. This deletes `_MIN_SEED_NODES` and `_SIGNIFICANCE_Z`, replaces them
with one measured biological prior, and degrades gracefully with no cliffs. Keep Option B in
reserve if the normal approximation is too coarse on real low-count seeds; pursue Option C
(empirical Bayes) once a library corpus exists, to retire even the fixed `od_0`.

## 7. Implementation sketch (when approved)

- In `fit_gdna_strand_overdispersion`: compute `od_mom = num/den` and `se` (already computed),
  then `od_post` per §3A; clamp `[0, ceil]`. Delete the `_MIN_SEED_NODES`/`_SIGNIFICANCE_Z`
  gate and `_OVERDISPERSION_FLOOR`-as-fallback; `fallback_used` becomes "no seeds at all"
  (`den ≤ 0` → `od_0`).
- `CalibrationConfig`: add `gdna_strand_overdispersion_prior` (= `od_0`, default = the measured
  value; documented as a biological prior, not a tunable threshold).
- Tests: thin data → `≈ od_0` (not 0); abundant overdispersed → recovered; abundant Binomial →
  `→ 0`. The existing recovery tests pass unchanged (strong data dominates the prior).

## 8. Open question for the user

**What is `od_0`?** Please supply the typical gDNA strand overdispersion from the real-data
studies (a single number in `[0, 1)`, e.g. an intra-class correlation). That value *is* the
default, and it is the only input the shrinkage needs.
