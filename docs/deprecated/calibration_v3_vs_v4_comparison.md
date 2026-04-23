# Calibration v3 vs v4 — A structural comparison

**Purpose.** Compare the current `rigel.calibration` (v4, mixture EM)
against the earlier "excellent calibration results" implementation at
commit `3bb0afd` (v3, *Aggregate-First Empirical EM*).  The goal is to
see what survived, what was re-architected, and which v3 ideas are worth
re-adopting.

**Code references.**
- v3: `git show 3bb0afd:src/rigel/calibration.py` (single file, 1210 lines,
  `calibrate_gdna` at the bottom).  Design doc: `docs/rigel_revision_plan/
  calibration_v3_aggregate_em.md` at that commit.
- v4: [src/rigel/calibration/_em.py](../../src/rigel/calibration/_em.py),
  [_calibrate.py](../../src/rigel/calibration/_calibrate.py),
  [_stats.py](../../src/rigel/calibration/_stats.py).  Design doc:
  [calibration_v4_em_mixture_plan.md](calibration_v4_em_mixture_plan.md).

---

## 1. Model at a glance

Both versions are two-component mixtures over genomic regions (G = gDNA-only,
R = expressed) fused as additive LLRs on logit γ, with a hard constraint
`k^s_i > 0 → γ_i = 0`.  Everything else differs.

| Aspect                                   | v3 (`3bb0afd`)                                                    | v4 (current)                                                                 |
|------------------------------------------|-------------------------------------------------------------------|------------------------------------------------------------------------------|
| **Primary density signal**               | Gaussian on `log(n/L + ε)`                                        | Poisson `k^u ∼ Poi(λ·E)` vs Poisson-LogNormal marginal                       |
| **Rate denominator**                     | Raw region length `L`                                             | `E_i = mappable_bp_i` (index-supplied)                                       |
| **RNA rate distribution**                | Per-class Gaussian `(μ_R, σ_R²)` on log-density                   | `log μ_i ∼ N(μ_R, σ_R²)` (LogNormal prior), 21-node Gauss-HermiteE quadrature |
| **Strand model**                         | Beta-Binomial with **shared κ** (Binomial fallback)               | Plain Binomial                                                               |
| **Strand channel gating**                | Always on; `SS=0.5` makes LLR vanish naturally                    | **Data-driven aggregate z-test** (z ≥ 3 on pooled sense fraction)            |
| **Strand LLR floor**                     | Symmetric `_EPS=1e-12` (same latent bug)                          | Asymmetric `max(ε_bio, ε_CI)` — see companion doc                            |
| **Fragment-length LLR**                  | Shape-normalised histogram, α = 1/M Dirichlet                     | **Identical** (shape-normalised histogram, α = 1/M Dirichlet)                |
| **Seed — RNA**                           | `k^s > 0 → γ=0`                                                   | Same                                                                         |
| **Seed — gDNA**                          | Density-eCDF bottom 10% + strand-symmetric filter                 | **None** (anchor-only; everything else starts γ=0.5)                         |
| **Magic numbers in init**                | `GDNA_INIT_DENSITY_PERCENTILE=0.10`, `abs(sense_frac-0.5)<0.1`    | **None** (pooled-MLE `λ_G₀`, no thresholds)                                  |
| **σ_R M-step update**                    | Over all eligible regions, (1-γ)-weighted                         | **Anchor-only** (hard-spliced regions only) — breaks circularity             |
| **λ_G M-step update**                    | γ-weighted Poisson rate over eligible                             | Same                                                                         |
| **κ estimation**                         | Golden-section maximum of **mixture** marginal log-lik, every iter | N/A                                                                          |
| **Zero-count / ineligible regions**      | Excluded from M-step; γ defaults to π                             | Ineligible if `E_i < mean_frag_len`; γ defaults to π_soft                    |
| **E-step prior**                         | π_soft (conditioned on unspliced)                                 | **Same** (π_soft)                                                            |
| **Per-chromosome gDNA rate**             | `gdna_density_per_ref` in output                                  | Not computed                                                                 |
| **Pristine-sample flag**                 | `init_diag["pristine_sample"]`                                    | Not computed                                                                 |
| **FL model construction**                | `build_gdna_fl_model(γ)` and `build_gdna_fl_model(1-γ)` every iter | Same pattern (plain histogram helper `_build_fl_histogram`)                  |
| **Output FL model smoothing/prior**      | Plain posterior histogram                                         | Blended with intergenic FL via `fl_prior_ess` pseudocount (ESS = 1000)       |

---

## 2. Where v4 is theoretically superior

### 2.1 The density signal is actually a count

v3 computed `log_d = log(n_total/L + ε)` and fit two Gaussians on that
scalar.  That collapses the real observation (a Poisson count `k^u` over
an exposure `E`) into a summary statistic and throws away the sample
size.  A region with `k^u = 1, E = 1 kb` and a region with `k^u = 100,
E = 100 kb` have the same log-density but enormously different
likelihoods for any given λ.  v3 treated them identically.

**v4's fix.** `P(k|z=G) = Poisson(k | λ_G · E)` and `P(k|z=R) =
∫ Poisson(k | (λ_G + μ) E) · LogNormal(μ ; μ_R, σ_R²) dμ` — the full
count model.  Small-`E` regions naturally carry less weight because
their Poisson log-likelihoods are flatter over λ.  The Gauss-HermiteE
quadrature on log μ makes the LogNormal marginal tractable in closed
form (21 nodes is overkill but the cost is trivial).

Practical consequence: low-coverage regions in v3 could be
mis-classified by tiny fluctuations in observed density; v4 correctly
integrates over their uncertainty.

### 2.2 Mappability as the Poisson exposure

v3's λ_G denominator was `L` (raw length), so any genic region
containing repetitive sequence or unmappable hg38 windows systematically
**under-reported** its Poisson exposure, inflating the apparent
per-region density and pulling those regions into class R.  v4 uses
`mappable_bp_i = mappable_effective_length` from the index, which is
the correct Poisson offset.

### 2.3 The anchor-only σ_R update breaks a circular loop

v3's M-step updates μ_R, σ_R² using `(1-γ)`-weighted moments over every
eligible region.  This is the textbook Gaussian-mixture update, but it
has a nasty feedback mode specific to this problem: empty or low-count
regions have log-density near the floor `log(ε)`, which is far below
the gene-expression mean.  If a few of those regions get γ < 1 at any
iteration (for any reason — prior, FL noise, strand LLR blip), their
near-floor log-density enters the R-class variance and **broadens** σ_R.
A broader σ_R makes the R-class PDF wider and flatter, which assigns
more posterior mass to R on low-count regions, which then get γ < 1,
which further broadens σ_R.  The result is a sigma-inflation runaway
that collapses λ_G.

v4 cuts that loop by updating μ_R, σ_R² from **only the hard-spliced
anchor regions**.  Those regions are definitionally RNA; γ never moves
away from 0; their log-rate moments reflect the actual expressed
distribution.  The count LLR then uses that anchor-derived σ_R as a
prior on *all* regions.  This is the structural improvement that made
v4 robust at real-genome scale.

### 2.4 Anchor-only init, no magic numbers

v3 seeded γ = 1 for the "bottom 10% by log-density eCDF" plus all
strand-symmetric regions (|sense_frac - 0.5| < 0.1).  Both thresholds
are magic and both fail in corner cases:
- Pristine samples have no bottom-10% gDNA, so the eCDF seed
  contaminates the gDNA class with genuinely expressed tail regions.
- Stranded libraries with `SS = 0.5` (i.e., unstranded) trigger the
  strand-symmetric clause on every region.

v4 uses one initialisation: `k^s > 0 → γ = 0`, everything else
γ = 0.5.  The first E-step's count LLR discovers the gDNA class from
its Poisson signature directly.  This removes two hyperparameters and
one failure mode.

### 2.5 Data-driven strand channel gating

v3 always computed the strand LLR; v4 pools
`(Σ k_sense, Σ n_unspliced)` across eligible regions and runs a
one-sided Binomial z-test against H0: `p_sense = 0.5`.  Only if
`z ≥ 3` is the strand LLR enabled.  This is a genuinely nice touch —
the question "is this library stranded?" should be answered by the
data, not by a config flag or an `SS` field the user might lie about.
v3's Beta-Binomial naturally suppressed the channel at `SS = 0.5` but
still paid a log-likelihood-evaluation cost and exposed itself to
noise floor bugs.

---

## 3. Where v3 was arguably better — candidates for re-adoption

### 3.1 Beta-Binomial strand model with shared κ (**strong candidate**)

v3's key insight on the strand channel was this: the per-region sense
fraction is **overdispersed** relative to a pure Binomial.  Sources:
gene-level strand-specific biases (some transcripts are longer, some
have secondary structure, some have cleavage bias), per-region depth
differences, PCR duplicates, R1/R2 asymmetry.  A plain Binomial treats
every sense/antisense read as an iid trial, which **inflates the LLR
per region at high counts** and makes the strand channel over-confident
precisely where it has the most leverage.

v3 solved this with a Beta-Binomial sharing a single dispersion
parameter κ:
- `z = G`: `BetaBin(k | n, κ/2, κ/2)` — symmetric around 0.5.
- `z = R`: `BetaBin(k | n, κ·SS, κ·(1−SS))` — centered at SS.
- `κ` estimated each iteration by **golden-section maximum** of the
  exact mixture marginal log-likelihood.

Properties:
- At `κ → ∞` degenerates to v4's Binomial.
- At finite κ automatically tempers overdispersion (the per-region LLR
  caps out at a data-dependent value determined by the mixture fit).
- At `SS = 0.5` LLR vanishes identically (G and R coincide).

**Why this matters to the bug we just fixed.**  Our ε-floor fix caps
per-antisense LLR at `log(0.5/ε)`, which is a *global* symmetric cap.
A Beta-Binomial with a learned κ gives you that cap *per region*,
scaling with the actual overdispersion in that region's counts — which
is statistically the right answer.  At vcap scale, the learned κ would
probably give similar protection to ε = 0.001 without the hand-picked
floor.

**Recommendation.** Seriously evaluate re-adopting the Beta-Binomial
strand model with shared κ.  The machinery already exists in git at
`3bb0afd:src/rigel/calibration.py`:
- `_compute_strand_llr_betabinom` (lines ~333–400)
- `estimate_kappa_marginal` (lines ~415–490, ~50 lines with golden-section)
- Two helper functions for Beta function evaluation.

Cost: ~150 lines added to `_em.py`, one extra golden-section call per
iteration (cheap), one extra output field (κ).  Benefit: principled
dispersion handling and an orthogonal backstop against the LLR blow-up
that the ε-floor addresses mechanically.

### 3.2 Per-chromosome λ_G reporting (**easy win**)

v3 returned `gdna_density_per_ref: dict[str, float]` — the γ-weighted
gDNA density per reference sequence.  Useful for diagnostics:
chromosome-specific mappability problems, mitochondrial/rRNA
contamination, control sequences, etc.  v4 dropped this; 15 lines of
code would bring it back.  Low-priority but cheap.

### 3.3 `pristine_sample` flag (**minor**)

v3 computed `pristine_sample = (pi_init < 0.02 or n_gdna_seed < 50)` as
a diagnostic.  v4 exposes `pi`, `pi_soft`, and `λ_G`, so it's
derivable, but a single boolean summary is nice for log-scanning.

### 3.4 Iteration-level κ re-estimation as a template for other
self-consistent updates (**structural**)

v3 re-estimates κ **inside** the EM loop using the current γ as
responsibilities.  This is the self-consistent pattern that the design
doc in [strand_llr_noise_floor_design.md](strand_llr_noise_floor_design.md)
§6 described as "layer (2)" for coupling the per-region `f_{G,i}` with
λ_G.  v3 already demonstrated the pattern works (and converges) for κ.
If we ever do implement layer (2), the template is `estimate_kappa_
marginal` at `3bb0afd:src/rigel/calibration.py:415`.

---

## 4. Where v3 was worse and we should not go back

### 4.1 Density-percentile gDNA seed

Hard-coded 10th-percentile seed fails on pristine samples and is a
magic number.  Don't reintroduce.

### 4.2 Strand-symmetric seed filter `|sense_frac - 0.5| < 0.1`

Triggers on everything when the library is unstranded.  Don't
reintroduce.

### 4.3 Gaussian density M-step over all eligible regions

The circular broadening loop described in §2.3.  Don't reintroduce.

### 4.4 Raw `L` as the Poisson rate denominator

v3 predates the mappability work; v4's `mappable_bp` is strictly better
and there is no reason to revert.

### 4.5 `_EPS = 1e-12` strand LLR clip

**Both versions had this bug** — we only hit it recently because vcap
runs at `ss_est ≈ 1.0` with modern strand models.  The ε-floor fix
from this week is orthogonal to v4 vs v3 and should be kept regardless
of which strand model (Binomial vs Beta-Binomial) we use going forward.

---

## 5. What carried over unchanged

1. **The hard-splice anchor** (`k^s > 0 → γ = 0`).  Definitionally correct.
2. **The FL channel** (shape-normalised histogram with α = 1/M
   Dirichlet, per-fragment Σ log-ratio scatter-added to regions).  The
   v4 helper `_fl_llr` is almost character-for-character identical to
   v3's `_compute_fl_llr`.
3. **π_soft as the E-step prior** (conditioned on soft/unspliced
   regions only).  Both versions discovered this independently; v3
   pushed it first.
4. **FL-model double-update each iteration** (γ-weighted gDNA FL,
   (1-γ)-weighted RNA FL, both rebuilt each iter).
5. **Logit-sum fusion of independent LLR channels** (count/density,
   strand, FL).
6. **Convergence criterion** `|Δπ_soft| < tol = 1e-4`.
7. **Golden-section search** as the univariate maximiser of choice
   (v4 still imports it implicitly via the architecture, even though it
   no longer uses κ).  If we re-adopt κ, v3's `_golden_section_max`
   is importable as-is.

---

## 6. Net assessment

**v4 is a theoretical improvement over v3 on the count/density channel,
the exposure denominator, the σ_R update, and the initialisation.**
These are the dominant failure modes we saw in practice and v4 closes
all four.

**v3's strand model (Beta-Binomial with shared κ) is theoretically
stronger than v4's plain Binomial** and is the one piece of v3
architecture that deserves a serious look for re-adoption.  The
overdispersion handling it provides is orthogonal to the ε-floor fix
we just shipped and would protect against the same class of
LLR-blow-up failures from a different direction.

**The ε-floor fix is independent of all of the above.**  Both v3 and v4
would have hit the same `_EPS = 1e-12` pathology under modern strand
trainers — we just happened to find it in v4.

---

## 7. Concrete recommendations

1. **Ship nothing new right now.**  v4 + the ε-floor fix works; vcap
   is healed.
2. **Open an issue / design doc** titled *"Re-adopt Beta-Binomial
   strand model with shared κ"* and cross-reference v3's
   `_compute_strand_llr_betabinom` and `estimate_kappa_marginal`.
   Pilot it on the mini-stress grid first (where the strand LLR
   stress-test amplifier is strongest) to see whether it collapses
   the residual 1.45–1.86× error at ss = 1.0 / gdna_frac ≥ 1.0.
3. **Restore `gdna_density_per_ref`** as a low-priority diagnostic
   (15 lines).  Useful for debugging per-chromosome contamination.
4. **Do not** restore the v3 density-Gaussian M-step, the
   density-percentile gDNA seed, or the strand-symmetric seed filter.
5. **Cross-link this document** from `METHODS.md` when the time comes
   to write an external-facing methods section, since it is the
   clearest statement of why the count/density channel was
   re-architected.

---

## 8. Appendix — file-level diff at a glance

| File                                    | v3 state                            | v4 state                                          |
|----------------------------------------|-------------------------------------|---------------------------------------------------|
| `src/rigel/calibration.py`             | 1210 lines, monolithic              | Removed — split into package                      |
| `src/rigel/calibration/_stats.py`      | —                                   | 141 lines, shared stats (unchanged semantics)     |
| `src/rigel/calibration/_em.py`         | —                                   | 765 lines, Poisson-LogNormal mixture EM           |
| `src/rigel/calibration/_calibrate.py`  | —                                   | 126 lines, orchestrator → CalibrationResult       |
| `src/rigel/calibration/_fl_model.py`   | —                                   | 102 lines, gDNA FL model with prior blending      |
| `src/rigel/calibration/_result.py`     | —                                   | 66 lines, frozen result dataclass                 |
| `src/rigel/calibration/__init__.py`    | —                                   | 18 lines, public re-exports                       |

v4 is roughly the same total LoC (~1218) as v3 (~1210) — not a
simplification in size, but a re-architecture in substance.
