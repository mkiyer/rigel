# Calibration V4 — Mappability-Aware EM Mixture Deconvolution

## 1. Synopsis of the Problem

Rigel must estimate a background genomic-DNA contamination rate `λ_G`
(fragments per base) from exonic regions alone — hybrid-capture panels
will not have intergenic data. Chr22 STAR+oracle sweeps show the
current density / strand / pool calibration stack overshoots `λ_G`:

- **Additive leak floor ≈ 4 × 10⁻⁵** at `gdna_fraction = 0`
  (should be 0).
- **Multiplicative residual 2–6 %** across `gdna_fraction ∈ [0.05, 2.0]`.
- Diagnostic decomposition attributes 100 % of the leak to exonic
  regions; the Tukey-biweight IRLS in `_density.py` fails to trim
  moderately-expressed genes and their RNA mass is absorbed into
  `λ_G`.

The structural problem is that **density regression is the wrong
frame**. The data is a **mixture of two populations** — silent exonic
regions (gDNA only) and expressed exonic regions (gDNA + RNA) — and
the estimator should deconvolve the mixture, not fit a robust line
through it.

We already built a mixture-EM for this (commit `3565dfe`, "new
calibration procedure is working reasonably well"). Two structural
defects in that version caused it to underperform and eventually be
replaced by the current regression stack:

1. **No mappability correction.** The exposure was raw region length
   `L_i`. Short or low-mappability exons biased the density estimate
   because `k_i / L_i` is not an unbiased estimator of the local rate
   when fragments cannot be uniquely assigned.
2. **Zero-count regions excluded.** `eligible = (n_total > 0)` dropped
   the very population that carries the strongest evidence for `λ_G`:
   silent genes whose exons legitimately drew zero fragments. This
   made the "silent seed" tiny and biased it upward toward regions
   that happened to catch a stray read.

Both defects are now fixable. BUG #2 (closed) supplies 1/NH-weighted
multimapper credits for `E_i` (`mappable_effective_length`), and the
same 1/NH weighting applies to `k_i`, making `k_i ~ Poisson(λ · E_i)`
an unbiased model even with multimappers and even for zero-count
regions.

## 2. Proposed Solution

A **mappability-aware, splice-anchored, three-signal EM mixture
model** over exonic regions. This is a direct descendant of the
`3565dfe` estimator, with the two structural defects fixed and the
count likelihood reformulated so spliced vs unspliced reads carry
distinct weights.

### 2.1 Generative model (per region *i*)

Let `E_i` be the mappable effective length (1/NH-weighted bases) and
`L_i` the raw span. Each region has a latent class
`z_i ∈ {G, R}` (gDNA-only vs expressed). Observed counts are split:
`k_i^u` = unspliced fragments (1/NH-weighted), `k_i^s` = spliced
fragments.

- **gDNA-only class (`z_i = G`).** Fragments arise from uniform
  Poisson coverage of mappable bases:
  $$k_i^u \mid z_i = G \;\sim\; \mathrm{Poisson}(\lambda_G E_i), \qquad k_i^s \mid z_i = G \;=\; 0 \text{ (hard)}.$$
  Strand split is binomial with success 1/2:
  $$k_{i,+}^u \mid k_i^u, z_i = G \;\sim\; \mathrm{Binom}(k_i^u,\, 0.5).$$

- **Expressed class (`z_i = R`).** Region has a latent RNA rate
  `μ_i` (fragments per base). Total counts combine gDNA background
  and RNA signal:
  $$k_i^u \mid z_i = R \;\sim\; \mathrm{Poisson}\big((\lambda_G + (1-\rho_i)\mu_i) E_i\big),$$
  $$k_i^s \mid z_i = R \;\sim\; \mathrm{Poisson}\big(\rho_i \mu_i E_i\big),$$
  where `ρ_i ∈ [0, 1]` is the region's splicing rate. We do **not**
  fit `μ_i` per region; instead we model the expressed population's
  **per-region rate** as heavy-tailed:
  $$\log\mu_i \mid z_i = R \;\sim\; \mathcal{N}(\mu_R,\, \sigma_R^2).$$
  This LogNormal prior is closed-form compatible with the Poisson
  likelihood under a plug-in MAP for `μ_i` at the E-step. It covers
  the extreme dynamic range of gene expression (orders of magnitude)
  without underflow, and it matches the observed log-density
  distribution of the spliced-evidence subset.

  Strand split has directional bias `p_i = SS` if the region is on a
  `+`-only gene, `1−SS` if `−`-only, `0.5` if ambiguous (two-gene
  overlap).

- **Mixing proportion.** `π = P(z_i = G)`. Estimated from data; we
  do not fix it, but we re-derive it against the soft (no-splice)
  subset (`π_soft`) when computing the E-step prior, exactly as in
  the old implementation.

### 2.2 Three signal channels (LLR decomposition)

The posterior for class `G` given region data factors into independent
log-likelihood ratios:

$$\mathrm{logit}\,P(z_i = G \mid \text{data}) = \log\frac{\pi_{\text{soft}}}{1 - \pi_{\text{soft}}} + \mathrm{LLR}^{(\text{count})}_i + \mathrm{LLR}^{(\text{strand})}_i + \mathrm{LLR}^{(\text{FL})}_i,$$

with the hard constraint `P(z_i = G | k_i^s > 0) = 0` applied first.

1. **Count LLR** uses the full Poisson comparison, not a Gaussian on
   log-density:
   $$\mathrm{LLR}^{(\text{count})}_i = \log \mathrm{Poisson}\big(k_i^u \mid \lambda_G E_i\big) - \log \int \mathrm{Poisson}\big(k_i^u \mid (\lambda_G + \mu) E_i\big)\, \mathrm{LN}(\mu \mid \mu_R, \sigma_R^2)\, d\mu.$$
   We evaluate the R-class marginal once per E-step via a small
   Gauss–Hermite quadrature (21 nodes). `k_i^u = 0` is now a valid
   observation with a meaningful count LLR, which is precisely the
   fix we need.
2. **Strand LLR** — unchanged from the old module (Binomial: G ∼
   Binom(n, 0.5); R ∼ Binom(n, p_i)). Disabled when SS ≤ 0.55 or
   `gene_strand == 0`. When unavailable it contributes exactly zero.
3. **Fragment-length LLR** — unchanged form (shape-normalised
   histograms with shared Dirichlet smoothing), but the FL models
   are trained only after iteration 0 so that the initial E-step
   runs on count + strand only.

### 2.3 M-step

- **`λ_G`** is a pooled Poisson MLE over the γ-weighted unspliced
  counts and exposures:
  $$\hat{\lambda}_G = \frac{\sum_i \gamma_i\, k_i^u}{\sum_i \gamma_i\, E_i}.$$
  Zero-count regions contribute `γ_i · 0` to the numerator and
  `γ_i · E_i` to the denominator — exactly the evidence that was
  missing in the old version.
- **`(μ_R, σ_R^2)`** are weighted MLEs on the log-rate, computed
  from the spliced-evidence subset at initialisation, then
  refined each M-step using the (1-γ)-weighted plug-in MAP rate
  `\hat{\mu}_i = k_i^u / E_i` (clipped at a small ε to avoid
  `log 0` on near-silent expressed regions).
- **`π`** and **`π_soft`** computed as before; `π_soft` conditions
  on `k_i^s = 0`.
- **gDNA FL model** rebuilt each iteration from γ-weighted
  fragments; RNA FL model frozen from the spliced-evidence
  fragments (gold standard — does not need γ-weighting).
- **Per-ref `λ_G^c`** is the γ-weighted Poisson MLE restricted to
  each reference.

### 2.4 Initialisation — splice-anchored two-phase seeding

Unchanged from `3565dfe` in shape, tightened in detail:

1. **Expressed seed:** `γ_i = 0` for every region with `k_i^s > 0`.
2. **gDNA seed:** among the `k_i^s = 0` set, set `γ_i = 1` for the
   bottom-decile of rate `k_i^u / E_i` (lowest-density silent
   candidates). This is a weak, deliberately generous seed — the EM
   will correct it.
3. For all other soft regions, `γ_i = π_init` (mixing proportion
   derived from the seed counts).
4. The initialiser never uses Tukey / IRLS / threshold tuning.

### 2.5 Convergence

Iterate up to 50 rounds; stop when `|Δπ_soft| < 10⁻⁴`. In practice
the old implementation converged in 8–15 iterations; we expect
similar with the new count channel.

## 3. Theoretical Justification

### 3.1 Why a mixture EM beats the regression stack

The regression stack (`_density.py`) treats every exonic region as a
noisy observation of the same underlying rate `λ_G`, and uses Tukey
biweight to "robustly" exclude outliers. This is statistically
misspecified: the outliers are not noise, they are a structurally
different population (expressed genes) drawn from a heavy-tailed
distribution. Tukey downweights proportional to residual size, but a
moderately expressed gene (2–10× `λ_G`) is not heavily downweighted,
and there are thousands of such genes per library. Their aggregate
pull on `λ_G` is systematic and in the direction we observe.

A properly specified two-component model with an explicit expressed
distribution handles this correctly: even the mildly expressed
genes have higher posterior mass on the expressed class than on
gDNA, because the expressed class has a fatter tail above `λ_G`
while the gDNA class has (essentially) a point mass at `λ_G`. The
classification is unambiguous at the 2–10× level.

### 3.2 Why LogNormal for the expressed rate

Real RNA expression spans 5–6 decades. A Normal-on-log-rate model:

- keeps the right tail heavy enough that extreme expressers do not
  drag `λ_G` up;
- keeps the left tail heavy enough that low expressers (ρ_i small,
  nascent-dominant, single-exon silent-ish) are not falsely
  classified as gDNA;
- is closed-form under Gauss–Hermite integration;
- is exactly what the old module used implicitly via the Gaussian
  on `log(n/L + ε)` density.

Any fatter-tailed alternative (Gamma-mixture, NegBinom with learnt
dispersion) would add parameters without data to inform them.

### 3.3 Why including zero-count regions matters

Silent exons with `k_i^u = 0` provide the strongest evidence that
`λ_G` is low. Under `Poisson(λ_G E_i)`, observing `k_i^u = 0`
contributes `-λ_G E_i` to the log-likelihood — a direct pull
**downward** on the `λ_G` estimate. The old module threw this
evidence away. The regression stack also throws it away (it fits in
density space, which is undefined at `k_i = 0 / E_i`). The new
model uses it natively.

### 3.4 Why splice-anchoring is safe

A spliced read is a physical indicator: gDNA cannot produce reads
spanning annotated junctions. Every other signal in this design is
probabilistic; the splice signal is deterministic and non-tunable.
It provides the expressed-class seed that bootstraps the whole EM.
Single-exon genes (no splice evidence) are NOT assumed silent —
they enter the soft classifier on equal terms.

### 3.5 Identifiability

With only one count channel we would be weakly identified because
`λ_G + μ_i` and `λ_G` trade off. Identifiability is rescued by:

- the large silent sub-population (`γ_i = 1`) providing a direct
  estimate of `λ_G`;
- the spliced-read hard carveout pinning the expressed class;
- the FL channel (once trained) discriminating long-ish RNA
  inserts from short-ish gDNA ones;
- the strand channel (when available) resolving the remaining
  ambiguous regions.

## 4. Required Setup

### 4.1 Data dependencies (already in place)

- `region_counts` — per-region unspliced/spliced, +/− counts with
  1/NH weighting (BUG #2).
- `region_df` — includes `length`, `tx_pos/neg`, `exon_pos/neg`,
  `ref`, and `mappable_effective_length`.
- `fl_table` — per-fragment (region_id, frag_len).
- Library `strand_specificity` estimated upstream from unique
  mappers.

### 4.2 Code dependencies

- `FragmentLengthModel` (unchanged).
- No new C++ extensions; the EM is pure NumPy/Python as in
  `3565dfe`.

### 4.3 Things we can reuse from `3565dfe`

| Piece | Reuse as-is? | Notes |
|:-|:-:|:-|
| `compute_region_stats` | adapt | add `exposure` field (already in `_stats.py`) |
| `compute_sense_fraction` | ✓ | unchanged |
| `_compute_strand_llr_binomial` | ✓ | unchanged |
| `estimate_kappa_sym` | ✓ | post-hoc diagnostic only |
| `build_gdna_fl_model` | ✓ | unchanged |
| `_compute_fl_llr` | ✓ | unchanged |
| `_seed_initial_partition` | adapt | drop the non-splice fallback; use rate `k^u / E`, not density `k / L` |
| `_compute_density_llr_gaussian` | **replace** | becomes `_compute_count_llr_poisson_ln` (Gauss–Hermite) |
| `_m_step` (λ_G branch) | adapt | Poisson MLE over `E_i`, not `L_i` |
| `_e_step` | adapt | use new count LLR, keep strand / FL LLRs |
| driver `calibrate_gdna` | adapt | new eligibility rule: any region with `E_i > 0`, zero counts included |

### 4.4 Tunables (all with defensible defaults)

| parameter | default | sensitivity | justification |
|:-|:-:|:-|:-|
| `max_iterations` | 50 | none | generous cap |
| `convergence_tol` | 10⁻⁴ | low | on `Δπ_soft` |
| `min_exposure` | 1 bp | negligible | eligibility floor |
| `init_gdna_percentile` | 0.10 | low | seed size, not output |
| `gh_nodes` | 21 | low | quadrature precision |
| strand-signal threshold | SS ≥ 0.55 | negligible | matches current |

No Tukey constant, no Anscombe, no IRLS, no dispersion prior, no
MAD-based φ init, no Fisher-information pool weights.

## 5. Implementation Details

### 5.1 New module layout

```
src/rigel/calibration/
    __init__.py
    _stats.py          # kept (already has `exposure`)
    _fl_model.py       # kept
    _result.py         # kept (GDNACalibration dataclass)
    _em.py             # NEW — core EM driver (replaces _density.py,
                       #       _strand.py, _pool.py)
    _calibrate.py      # thin public entry, delegates to _em.py
    _density.py        # DELETE
    _strand.py         # DELETE
    _pool.py           # DELETE
```

The existing `GDNACalibration` dataclass already contains all output
fields we need (`region_posteriors`, `gdna_density_global`,
`gdna_density_per_ref`, `kappa_sym`, `gdna_fl_model`,
`mixing_proportion`, `expressed_density`, diagnostics).

### 5.2 `_em.py` structure (~350 lines)

```python
# --- small utilities -------------------------------------------------
_GH21_NODES, _GH21_WEIGHTS = np.polynomial.hermite_e.hermegauss(21)

def _poisson_logpmf(k, mean):
    # stable: k*log(mean) - mean - lgamma(k+1)

def _count_llr_poisson_ln(k_u, E, lam_G, mu_R, sigma_R):
    """Return per-region log[P(k_u | λ_G E) / ∫ P(k_u | (λ_G+μ)E) LN(μ) dμ].
    Gauss-Hermite at nodes x_j with weights w_j on
    log μ ≈ μ_R + √2 σ_R x_j.
    """

# --- phases ----------------------------------------------------------
def _seed(stats, exposure): ...        # splice-anchored two-phase
def _e_step(stats, exposure, params, fl_models): ...
def _m_step(stats, exposure, gamma): ...

# --- driver ----------------------------------------------------------
def calibrate_gdna(region_counts, fl_table, region_df,
                   strand_specificity, *,
                   max_iterations=50, convergence_tol=1e-4,
                   diagnostics=False) -> GDNACalibration: ...
```

### 5.3 Eligibility and exposure

- **Old rule:** `eligible = (n_total > 0) & (L > 0)` — wrong.
- **New rule:** `eligible = (E_i ≥ 1)` — accept zero-count regions.

All M-step aggregations use `E_i`, not `L_i`. The `log(n/L + ε)`
density channel is removed entirely; replaced by the Poisson count
LLR that consumes `E_i` natively.

### 5.4 Per-region rate aliases for the expressed distribution

The expressed M-step needs a per-region rate observation. We use the
plug-in:

```
mu_hat_i = max(k_u_i / E_i - lam_G, eps_mu)
```

with `eps_mu = 1 / median(E_i)`. Then `(μ_R, σ_R²)` are the
(1-γ)-weighted moments of `log μ_hat_i` restricted to eligible
regions. This is the same log-density M-step as before, just with
the gDNA contribution subtracted and the mappable exposure in the
denominator.

### 5.5 Initialisation changes from `3565dfe`

- **Exposure-based rate** everywhere the old code used
  density-from-length.
- **Zero-count regions** are part of the gDNA seed when their
  percentile rank is low; previously they were dropped.
- No change to the spliced-region hard seed.

### 5.6 Test plan

Port the three most-informative calibration tests from the current
suite onto the new module, delete the rest:

- `test_calibration.py::TestStrandLLR::test_biased_toward_ss_favors_rna`
  (known pre-existing failure — re-evaluate against the new model;
  expect to pass).
- `test_calibration_stress.py::*` — keep the stress grids but
  recompute golden lambdas from the new model.
- `test_calibrated_density.py` — replace with
  `test_calibrated_count_llr.py` exercising the Poisson-LN channel.

Add unit tests:

- Zero-count region admissibility and its contribution to `λ_G`.
- Gauss–Hermite integration tail stability for large `σ_R`.
- Unstranded (`SS = 0.5`) hard-off of the strand channel
  (`LLR_strand ≡ 0`).
- Mappability reduction (synthetic `E_i = 0.3 L_i`) — `λ_G` estimate
  scales appropriately.

### 5.7 Benchmark plan

1. Rerun `scripts/benchmark/chr22_calibration_sweep.py` (oracle).
2. Rerun `scripts/benchmark/chr22_star_calibration_sweep.py` (STAR).
3. Compare against golden results in
   `scripts/benchmark/results/calibration_{oracle,star}_chr22/`.

### 5.8 Success criteria

- **Leak floor** at `gdna_fraction = 0`: `λ_G ≤ 5 × 10⁻⁶` on both
  oracle and STAR (currently 2.6 × 10⁻⁵ / 4.1 × 10⁻⁵ — ≥ 5× better).
- **Multiplicative residual** at `gdna_fraction ∈ [0.05, 2.0]`:
  `|λ̂/λ_true − 1| < 2 %` after leak subtraction (currently 2–6 %).
- **Unstranded (SS = 0.5)**: strand channel contributes exactly
  zero, verifiable by unit test.
- **Synthetic benchmarks**: all existing `scripts/benchmark/configs/*`
  sweeps at least as accurate as current.
- **Full pytest**: all previously-passing tests continue to pass;
  `test_biased_toward_ss_favors_rna` expected to flip to pass.

### 5.9 Migration order

1. Port `_em.py` driver + tests; keep `_density.py`, `_strand.py`,
   `_pool.py` in place temporarily.
2. Switch `_calibrate.py` to delegate to `_em.py` behind a feature
   flag (`RIGEL_CALIBRATION_V4=1`).
3. Run full benchmark + regression suites under both paths.
4. Once V4 meets success criteria, remove the flag, delete old
   modules, update callers.
5. Clean up memory / doc references to the old density+strand+pool
   stack.

## 6. What This Plan Is NOT

- Not a rewrite of fragment scoring or EM quantification.
- Not a change to the index or the BAM scanner.
- Not a fallback onto intergenic regions.
- Not a Tukey / Anscombe / IRLS tuning exercise.
- Not a new C++ extension — stays in Python / NumPy.
- Not an attempt to re-derive every statistic; this is a restoration
  of the `3565dfe` mixture EM with the two known defects repaired
  and with proper count-likelihood handling of the `k_i = 0` case.
