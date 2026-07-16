# NPMLE gDNA prior — current struggles (for review & to guide the search)

**Status:** prototype findings on branch `calib-ambig-init-wip`. This documents why the count-space
Poisson-rate **NPMLE** (plan of record in `npmle_roadmap.md`, motivated in `gdna_prior_zero_handling.md`)
is producing **spiky / jagged / brittle** fits rather than the smooth, robust prior we want — and lays out
the tradeoffs so we can choose a better estimator. It also precisely describes the prototype implementation
so it can be code-reviewed.

> **RESOLUTION (§9, added after review):** the spikiness is fixed by **not running the NPMLE optimizer at
> all** — replace the EM-to-convergence with the **population average of the per-node likelihood curves**
> (`method="avg"`; a deterministic kernel sum, one EM step from uniform, NO splines/λ/argmin). Every fit is
> now smooth on every substrate. See §9 and `figures/npmle_*_avg.png` vs `figures/npmle_*_npmle.png`.

All evidence is from two real hybrid-capture cfRNA caches, scan-cached once
(`scripts/debug/calib_cache.py`) and fit on the compact per-pass beliefs (`scripts/debug/npmle_fusion.py`
+ `--belief-cache`): **LBX0190** (near-pristine, ~94 % of the confident set is count-0) and **vcap**
(~20–25 % gDNA in-silico mix). Figures are in `docs/calibration/figures/`.

---

## 1. What we are estimating (recap)

Each chain node has a latent gDNA rate `ρ` and observes a gDNA count `g ~ Poisson(ρ·E)` over a gDNA
effective length `E`. We want the **population rate distribution** `P(ρ)` — the prior that, via the
projection (`bp_solver._kde_logprior`), tells the per-node solver "at observed density `d = M/E`, how much
of it is plausibly gDNA." Zero must be representable (a pristine node has `ρ=0`), which the density-space
KDE cannot do (it floors at `1/E`); that is why we moved to a **count/rate-space** estimator.

The **NPMLE** is the nonparametric maximum-likelihood estimator of `P(ρ)`: choose weights `{w_j}` on a fixed
log-ρ grid `{ρ_j}` to maximise `Σ_i log Σ_j w_j · L_i(ρ_j)`, where `L_i(ρ_j)` is node `i`'s likelihood at
rate `ρ_j`.

---

## 2. The prototype implementation, precisely (for scrutiny)

All in `scripts/debug/npmle_fusion.py`. Symbols: `k` = total unspliced count (`u_pos+u_neg`), `f_g` = the
current-belief gDNA fraction, `ĝ = f_g·k` = believed gDNA count, `E = eff_global` (region: contained gDNA
eff-length; boundary: `½(E_left+E_right)`), `τ² = Var(log f_g) = belief.var_gdna`.

**Fit set (`node_obs`).** Every node with **finite `Var(log f_g)`** (solved or structurally locked). Nodes
with `Var=∞` (AMBIG at init, empty nodes) are excluded. A `substrate` switch restricts this to `region`
(drop boundaries) or `structural` (intergenic+intron regions only — belief-independent).

**Per-node likelihood — the Poisson-lognormal (`_loglik_pln`).** The count is `g ~ Poisson(ρ·E)`; the belief
places `log g ~ N(log ĝ, τ²)`. Marginalise `g` by `n_gh=7`-point Gauss-Hermite quadrature `x_q, w_q`:

```
g_q      = ĝ · exp(√2 · τ · x_q)                          # quadrature counts
logL(ρ_j)= logΣ_q [ log(w_q/√π) + g_q·(log ρ_j + log E) − ρ_j·E − lgamma(g_q+1) ]
```

- `ĝ=0` (every count-0 or believed-pure-RNA node) ⇒ `g_q=0` ⇒ `logL = −ρ_j·E` (the exact zero anchor).
- high `ĝ` ⇒ the Poisson sharpens, the lognormal width dominates ⇒ Gaussian-in-log-rate (seam is automatic).
- This is the "honest width" idea: `τ` (belief uncertainty) *broadens* each node's likelihood.

**Cell collapse (`_collapse`).** `L_i` depends only on `(ĝ, E, τ)`, so nodes are binned into cells keyed on
`(log ĝ | a shared zero-bin, log E, τ)` (`dlog=0.1`, `dt=0.25`); the representative is the in-cell mean, the
weight is the node count. ~10⁶ nodes → ~10³–10⁴ cells → the EM is sub-minute. Binning error is a perf knob.

**EM (`fusion_fit`).** Standard mixture EM over cells with cell-count weights: `r_cj ∝ w_j·L_c(ρ_j)`;
`w_j ∝ Σ_c count_c·r_cj`; ~200 iters, tol 1e-6, deterministic (no RNG). Grid: `log ρ ∈ [log(0.05/E_max),
log(3·max ĝ/E)]`, 300 points.

**Projection read (`smooth_logprior` + `projection_curve`).** For the plots only: convolve the discrete `w`
with a Gaussian of width `h_dex=0.15` decades (real tails), add a uniform `eps=0.02` floor **bounded to the
fitted support** `[q0.5%, q99.5%]` (fills interior valleys, leaves tail decay intact). The prior over `f_g`
at density `d` is `logP(f_g·d) + Jeffreys(−log(1−f_g))`; we plot the posterior-mean `f_g`.

**How the beliefs are obtained (`_run_sweep`).** We reproduce production `calibrate` exactly to get the
belief at pass 0 (init), pass 1 (sweep, no KDE), pass 2 (sweep + production KDE), then fit the NPMLE on each.

---

## 3. The core struggle: the NPMLE is ATOMIC (spiky) *by construction*

**This is the answer to "ALL NODES should give the smoothest fit."** It does not — and cannot — because of
a fundamental property of the estimator, not the data.

By the **Kiefer–Wolfowitz / Lindsay** theorem, the NPMLE of a mixing distribution is a **discrete measure
with at most as many atoms as there are distinct likelihood vectors**, and in practice it concentrates on a
*handful* of spikes. Maximum likelihood, with no smoothness constraint, is *maximised* by a spiky discrete
`P(ρ)`. **More data does not smooth it** — it can only add atoms up to the number of distinct `(ĝ, E, τ)`
cells, and the optimiser still piles mass onto a few. Smoothness in a KDE comes entirely from the imposed
bandwidth; the NPMLE imposes nothing, so it is jagged.

The figures confirm it directly: **vcap `all`** (750 k nodes — the most data we have) is still three sharp
peaks (`figures/npmle_vcap_all.png`), and **LBX0190 `all`** is a sharp comb (`figures/npmle_LBX0190_all.png`).
Volume did not help.

⇒ If we want smooth + robust (the goal), the raw NPMLE is the wrong estimator. We need a **regularised /
smooth deconvolution** (see §7). This is the central finding.

---

## 4. The ALL-NODES approach, in detail — why it is spiky *and* biased high

Two distinct problems compound (they are separable — smoothness vs bias):

**(a) Spikiness** — §3, the NPMLE atomicity. Independent of substrate.

**(b) Upward bias / brittleness of the fitted rate**, from three mechanisms:

1. **The discreteness comb.** 28 % of LBX0190's fit set is `ĝ≈1` over near-fixed short boundary
   `E≈115 bp` (96 % boundaries). `Poisson(1 | ρE)` peaks at `ρ=1/E≈10⁻²`, and thousands of these nodes share
   the same `E` ⇒ the NPMLE plants a sharp atom at exactly `10⁻²`. This is quantization masquerading as a
   real gDNA mode (the same comb that defeated the density KDE — `gdna_prior_zero_handling.md` §4).
2. **GIGO.** The fit uses the *current belief* `ĝ = f_g·k`. On near-pristine LBX0190 the belief already
   reports **`Σĝ/Σk = 0.57`** (57 % of reads called gDNA — a gross over-call, inherited from the current
   floor+KDE). Fitting `P(ρ)` on this belief *reproduces the over-call*.
3. **Population concentration defeats the honest width.** We hoped `τ` (belief uncertainty) would down-weight
   the unreliable comb nodes (they have `τ≈1.2`, large). It does not: **width broadens each node's
   likelihood but does not move its centre, and unit-mass observations are never down-weighted.** 60 k broad
   likelihoods all *centred* at `10⁻²` still sum to a mode at `10⁻²`. Width-as-spread ≠ down-weighting. This
   is the empirical refutation of roadmap decision **D2** ("honest width makes fit-all safe").

So `all` is *both* spiky (a) and biased high (b). You said you would accept an overestimate — but note the
bias here is not gentle (a pristine sample called majority-gDNA), and the spikes make it *brittle*, not just
biased: the projection latches onto individual atoms (see §5).

**Why `structural` is clean.** Intergenic+intron regions are gDNA-clean **by structure** (their mass is
gDNA regardless of any belief), so fitting them is **belief-independent** — no GIGO, no comb (they are
long-E regions, not the short-E boundary comb). Result: a clean single mode at `10⁻⁶` (LBX0190) / `10⁻⁴`
(vcap), stable across passes (`figures/npmle_*_structural.png`). This is the count-zero-information
principle reasserting itself: the *structural* nodes are the honest gDNA observations.

**The catch.** `structural` has no exons, so it never sees the **capture-enriched on-target** mode. Its
projection then crushes high-density on-target nodes to `f_g→0` — i.e. calls enriched gDNA "RNA" (the
capture leak). Visible in `figures/npmle_vcap_structural.png` (right panel → 0 at high density). So:

| substrate | pristine | capture | spiky? |
|---|---|---|---|
| `all` | biased high (57 % gDNA), comb at 10⁻² | broad, has an enriched mode | **yes** |
| `region` | clean at pass 0/1; small over-call bump by pass 2 | intermediate | some |
| `structural` | **clean** (mode 10⁻⁶) | **misses enrichment** → leak | least |

---

## 5. The projection is brittle at the ends of the spectrum

The prior only helps a node if the projection is well-behaved. Two failures:

- **Tiny atoms → false positives.** False-positive suppression needs `P(ρ) ≈ 0` *above* the background mode
  (a high-density node must find no plausible high gDNA rate → `f_g→0`). But the NPMLE sprinkles tiny atoms
  (weight ~5e-4) out to `10⁻²`; each becomes a small "acceptable rate," so the projection *rises* again at
  high density (`figures/npmle_LBX0190_structural.png`, right panel oscillates up to 0.85). A single spurious
  atom re-injects gDNA. This is the brittleness.
- **Flat-tail inversion (fixed, but instructive).** An earlier bug applied the uniform floor across the whole
  grid, flattening the high-ρ tail; with no penalty there, the RNA-parsimony Jeffreys drove `f_g→1`
  everywhere — a fully inverted projection. Bounding the floor to the fitted support restored real decay.
  The lesson (same as production's `logpdf_kernel` vs the clamped `logpdf`): the tail penalty is load-bearing
  and *fragile*.
- **Low-density end.** A node whose *total* density equals the background rate is correctly called gDNA
  (`f_g→1`); that is right. But it means the low end is entirely at the mercy of where the fitted mode sits —
  another reason a jagged fit is dangerous.

---

## 6. The underlying tension (independent of estimator)

Even with a perfect estimator, two goals pull apart:

- **Non-circularity** wants a belief-independent substrate (structural) → clean, but no enriched mode.
- **Capture-enrichment coverage** wants exons/boundaries → has the enriched mode, but is belief-dependent
  (GIGO) or quantization-prone (comb).

The enriched mode must come from a source that is **both on-target and belief-independent**. Our earlier
validated idea (memory `ambig_g1emit_stratified_kde`, 53 % recovery): the **RNA-free single-strand
splice-junction-adjacent exons** — on-target by structure, `f_g≈1` by the spliced=0 gate (oracle 0.997).
That is a structural enriched sample, not "all exons via belief."

---

## 7. Tradeoffs & candidate directions (the search)

| estimator | smooth? | zero-native? | discreteness-robust? | circular? | notes |
|---|---|---|---|---|---|
| density KDE (current prod) | ✅ | ❌ (floors 1/E) | ❌ (comb→fake mode) | fit on solved nodes | what we're replacing |
| raw NPMLE (prototype) | ❌ **spiky** | ✅ | ✅ | depends on substrate | this doc |
| **smoothed/penalised deconvolution** | ✅ | ✅ | ✅ | depends on substrate | the likely answer |

**Directions to evaluate:**

1. **Smooth the deconvolution (address §3 directly).** Keep the count/Poisson likelihood (zero-native) but
   estimate a *smooth* `P(ρ)` instead of the atomic MLE:
   - **g-modeling** (Efron 2016): `log P(ρ) = spline(log ρ)`, fit by penalised Poisson likelihood. Smooth by
     construction; a handful of parameters; still a genuine deconvolution. **My lead candidate.**
   - **Penalised NPMLE / P-spline**: NPMLE `+ λ·roughness(log w)`. One knob `λ`; `λ→0` = NPMLE, `λ↑` = smooth.
   - **Fixed-kernel mixture (a "Poisson-space KDE")**: each atom is a kernel of fixed width, not a δ — smooth,
     and honours your intuition that all-nodes should give a stable, robust (if biased) curve.
   - (Cheap fallback: a coarser grid ⇒ fewer possible atoms ⇒ less spiky, but crude.)
2. **Choose the substrate deliberately (address §4b bias / §6 tension).** Likely **structural depleted +
   structural enriched (RNA-free SJ-adjacent exons)** → a non-circular bimodal `P(ρ)`. Fitting *all* nodes is
   robust for smoothness (most data) but reintroduces GIGO — so if we fit all, we must accept (and bound) an
   upward bias, and rely on the strand/messages/EM downstream to correct it.
3. **Explicit zero handling.** Even structural count-0 nodes give a *flat* `e^{−ρE}` (not a spike at 0); an
   explicit zero-atom (`P = π·δ₀ + (1−π)·smooth`) may be needed so the "gDNA absent" mass concentrates rather
   than spreading — the zero-inflation of `gdna_prior_zero_handling.md` §3.2, now motivated by data.

The reconciliation of your intuition: **use as much data as possible (robustness) with a *smooth* estimator
(not the atomic NPMLE), and get the enriched mode from a *structural* (non-circular) source.** Smoothness and
non-circularity are separate levers; the raw NPMLE gets neither for free.

---

## 8. Resolution of the spikiness: average the kernels, don't optimize (no splines)

The reviewer's diagnosis is right — we need smoothness — but their fix (g-modeling / penalized-spline) is
exactly the spline machinery that caused our past `GCV-λ` argmin nondeterminism
(`calibrate_cross_process_nondeterminism.md`). There is a smoother-*and*-more-robust option that needs **no
optimization at all**.

**First principle.** The KDE was smooth not because it had "more data" but because it is a **deterministic
sum of fixed kernels**, not an optimizer. The NPMLE is a maximum-likelihood optimizer, and ML *wants* atoms
(§3). So: stop optimizing. A node's `(count, length)` is not a point — it is a smooth likelihood curve over
the rate ρ (a count-1 node over E≈115 is `Gamma(≈2,115)`, ~a decade wide). **Let each node contribute its
own likelihood curve and average them:**

```
P̂(ρ) = (1/N) Σ_i L_i(ρ)          # L_i = the per-node Poisson-lognormal curve, unit-mass normalized
```

This is a KDE **in rate space**, where the kernel is the node's own likelihood and its **bandwidth is set by
the physics**: wide for low-count/short-E (so the boundary comb *smooths itself away*), narrow for
high-count. No bandwidth knob, no λ, no argmin, no spline — a deterministic weighted sum ⇒ maximally
reproducible. Zero-native (`k=0 → e^{−ρE}`). It is exactly **one EM step from uniform weights**; the NPMLE is
that EM run to convergence. The average is the smooth extreme; convergence is the spiky extreme. (It is
*oversmoothed*, like any KDE — the accepted, safe behaviour.) This realizes the "3-D KDE over (count,
length)" instinct correctly: collapse the (count, length) observation onto the one physical axis (the rate),
with the Poisson supplying the bandwidth.

**Result (`method="avg"` in `npmle_fusion.py`):** smooth on **every** substrate, deterministic, and ~10–50×
faster than the EM (a single sum). `figures/npmle_LBX0190_all_avg.png` — the comb at 10⁻² is gone, dissolved
into a smooth curve; `figures/npmle_*_structural_avg.png` — clean smooth background, stable across passes;
`figures/npmle_vcap_all_avg.png` — smooth, and even *sharpens* pass 0→2 (σ 2.02→1.61: the weak→sharp
behaviour). Contrast the spiky `figures/npmle_*_npmle.png`.

**ADOPTED refinement — the Fixed-Kernel Poisson Mixture (`method="kernel"`, reviewer's synthesis).** `avg`
is *oversmoothed* — it never deconvolves, so it cannot separate the depleted from the enriched mode (which
capture needs). The better estimator parameterises `P(log ρ) = Σ_j w_j·N(log ρ; log ρ_j, h²)` — a mixture of
**fixed-width** Gaussian kernels (`h = kernel_bw` decades) — and fits the weights `w_j` by **standard EM**
(monotone, arithmetic-only, deterministic — no spline, no λ, no matrix inversion; the `GCV-λ` fragility that
caused past nondeterminism is absent). The EM **deconvolves** (recovers the true modes) but the fixed kernel
forbids any peak sharper than `h` ⇒ no bed-of-nails. `h` is precisely the KDE bandwidth knob. It subsumes the
spectrum: `h=0` = the atomic NPMLE, `h→∞`/no-iterate = `avg`. `figures/npmle_*_kernel.png`: sharp *and*
smooth *and* stable (LBX0190/structural: σ=0.62 at 10⁻⁶). **This is the fitter we adopt.**

**Projection ↔ Jeffreys balance — diagnosed (still to calibrate in the solver).** The prior-only projection
(`prior × RNA-Jeffreys`, no strand/messages) pins f_g high because the RNA-parsimony **Jeffreys `−log(1−f_g)`
is UNBOUNDED** (+12 at the log-odds grid top) while the prototype **floors the prior tail** (`max(·,_EPS)`)
— so near f_g=1 the unbounded Jeffreys beats the capped penalty. This is the documented "clamped-tail →
false-positive" failure (production's `logpdf_kernel` uses REAL quadratic tails for exactly this reason). It
is a *consumption* issue, not a fit issue, and the prior-only projection is a worst case: in the real solve
the prior is ADDED to the strand likelihood + neighbour messages, which counteract the Jeffreys. The right
place to calibrate the Jeffreys/tail balance is **in the solver, against oracle** — not this viz.

**What `avg` does NOT fix (orthogonal, still open):**
- **Substrate bias.** `all` is smooth but still biased high (GIGO — the 57 % belief over-call); `structural`
  is smooth *and* clean (near-zero for pristine). Bias is a substrate choice, not an estimator choice — the
  reviewer's **dual-structural** (depleted intergenic+intron **+** structural-enriched RNA-free SJ-adjacent
  exons) is the belief-independent substrate that is smooth, clean, *and* carries the enriched mode. Combine:
  **`avg` estimator (no splines) on the dual-structural substrate.**
- **Projection ↔ Jeffreys balance.** The pure-prior projection (prior × RNA-Jeffreys, no strand/messages) is
  delicate: a soft high-ρ prior tail lets the Jeffreys (which pushes f_g→1) dominate → f_g→1 everywhere; the
  earlier flat-tail-floor bug let it invert the other way. This is a *projection/consumption* question (the
  Jeffreys strength, the tail sharpness, and — crucially — that in the real solve the prior is ADDED to the
  strand likelihood + messages, which counteract the Jeffreys), NOT a property of the fit. The **fit** (left
  panels) is the validated deliverable; the projection needs its own calibration pass.

---

## 9. Files to send for code review

**Core (the NPMLE math + projection — review these first):**
- `scripts/debug/npmle_fusion.py` — `node_obs`, `_collapse`, `_loglik_pln`, `fusion_fit` (both `avg` and
  `npmle`), `smooth_logprior`, `projection_curve`, and `_run_sweep` (how beliefs are produced).
- `scripts/debug/npmle_diag.py` — the fit-set dissection (count-0 share, the comb, `Σĝ/Σk`, `var_g`).

**Plan & motivation:**
- `docs/calibration/npmle_roadmap.md` — the plan of record (§4 records decisions D1/D2/D3, which this
  refutes).
- `docs/calibration/gdna_prior_zero_handling.md` — why count/rate space; the zero & discreteness problem.
- `docs/calibration/npmle_struggles.md` — this document.

**Production seam it will replace (context for the reviewer):**
- `src/rigel/calibration/gdna_density_prior.py` — `GdnaDensityPrior` (the KDE) + `build_training_substrate`
  (the current substrate). The NPMLE replaces both.
- `src/rigel/calibration/bp_solver.py` — `_kde_logprior` (the projection, ~line 400), `_global_logprior` +
  `_floor_estimate` (the current floor), and `node_sweep`'s prior-build block (~lines 664–724).
- `src/rigel/calibration/simplex_logodds.py` — `_solve_nodes_logodds_all`: the `(n_nodes, K)` additive
  `global_logprior` contract the projection must produce (the clean seam; unchanged by any of this).

**Harness (reproduce the figures):**
- `scripts/debug/calib_cache.py` — scan-once/cache-everything.
- caches at `~/Downloads/rigel_runs/cfrna/_calib_cache/{LBX0190,vcap}.pkl`; index at
  `~/Downloads/rigel_runs/rigel_index`.
- reproduce: `python scripts/debug/npmle_fusion.py --cache …/LBX0190.pkl,…/vcap.pkl --out
  docs/calibration/figures --belief-cache <dir> --substrate all,region,structural --method avg,npmle`
