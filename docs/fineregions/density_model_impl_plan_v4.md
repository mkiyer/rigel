# gDNA Density Model Implementation Plan v4

Date: 2026-05-24
Status: implementation ready
Companion: `calibration_roadmap_v3.md`
Supersedes: `density_model_impl_plan_v1.md`, `density_model_impl_plan_v2.md`,
`density_model_impl_plan_v3.md`

## 0. What Changed Relative to v3

v3 was architecturally sound but had several open seams flagged in review.
v4 closes them and is the only plan to implement against. Concrete deltas:

1. **Pearson-scored trim (§7.2).** The robust MoM trims by standardized
   squared residual `(C - mu)^2 / max(mu, eps)`, not raw `(C - mu)^2`. Raw
   squared residuals correlate with `mu`, so v3 preferentially trimmed
   long/high-opportunity rows even under the Gamma–Poisson null and biased
   `phi` downward. Pearson scoring removes that bias at the same cost.
2. **`rho_ref` and the depth-2 fallback are specified (§7.2.4, §7.3).**
   `rho_ref = ALL.mean_density` when the all-anchor fit succeeds, otherwise a
   precision-weighted mean of family priors, with an explicit
   `rho_ref <= 0` path that emits a deterministic zero-density evidence
   surface. The depth-2 fallback prior is `Gamma(alpha_fallback,
   alpha_fallback / rho_ref)` with `alpha_fallback = 1` — a genuinely weak
   prior — not `Gamma(rho_ref * beta_cap, beta_cap)`.
3. **Closed-form `expected_tail_count` (§7.4).** The NB stop-loss identity
   `E[(X - c)+] = mu * Sf_{r+1}(floor(c) - 1) - c * Sf_r(c)` replaces the
   branchy truncated-sum / normal-approximation logic. One vectorised SciPy
   call per side.
4. **Fractional-count predictive convention (§7.4).** All discrete CDF/SF
   evaluations consume `floor(C_r)` and `floor(B_tot)`. The continuous
   posterior identities consume the unrounded values. Locked and tested.
5. **`DensityObservation` is transient (§8.1).** It does not live on
   `CalibrationResult`. `DensityEvidence` is shrunk to the arrays downstream
   code actually consumes; full-resolution diagnostics are summarised, not
   stored per-region.
6. **`RegionExposure.from_density` lands inside Phase 4 (§8).** The
   ordering-dependency between v3's Phase 4 and Phase 5 is removed. v4 has
   one fewer phase and a single PR boundary for "density evidence becomes the
   exposure source." The leftover Phase 5 in v4 is exposure cleanup and
   `A_r > 1` audit only.
7. **`RegionExposure` contract update is explicit (§9).** The `(0, 1]` clip
   in `_exposure.bp_weighted_mean_exposure_over_blocks`, the `A_r` docstring
   in `exposure.py`, and the `mode` Literal are all named as code edits with
   acceptance rules for ineligible regions.
8. **`RegionCountLedger` absorbs `PayloadArrays` totals (§5).** There is no
   parallel "typed view" class. `PayloadArrays` keeps its sorted matrix and
   per-channel projections; `RegionCountLedger` is a thin object that owns
   *only* the derived totals (unspliced/spliced/contained/boundary) and the
   `(POS, NEG)` channel views. The duplicated unspliced-total fields on
   `PayloadArrays` are removed in Phase 1, not deferred.

All other locked decisions from v3 §1 stand unchanged.

## 1. Decisions Locked From v3 (Recap, Unchanged)

- Background prior is fit only on contained anchors (INTERGENIC, INTRON-only).
- Boundary counts are local Bayesian evidence, never prior training data.
- Only unspliced channels enter the density likelihood.
- Estimator is robust, exposure-adjusted method of moments with an upper trim.
- Prior precision is bounded by a coefficient-of-variation floor
  (`beta_cap = 1 / (rho * CV_min^2) = 400 / rho`), a depth-invariant scale.
- Density evidence is independent of strand deconvolution and EM.
- `RegionExposure.A_r` from density is a relative scalar that may exceed 1.
- v4 does not assemble locus priors or unblock `quant_from_buffer`.

Notation (unchanged from v3):

```text
C_r   = contained_unspliced_pos + contained_unspliced_neg
B_l   = boundary_left_unspliced_pos + boundary_left_unspliced_neg
B_r   = boundary_right_unspliced_pos + boundary_right_unspliced_neg
B_tot = B_l + B_r
L_c   = contained effective length under gDNA FL
L_b   = L_l + L_r (per-side boundary opportunities)
```

## 2. Phase Map

```text
Phase 0 — Cleanup of stale density and strand scaffolding (v3 §4, unchanged)
Phase 1 — RegionCountLedger as the single totals provider (§5)
Phase 2 — DensityObservation transient builder (§6)
Phase 3 — DensityModel: robust Pearson-trim MoM, NB predictive, closed-form
          stop-loss (§7)
Phase 4 — Result wiring, orchestrator wiring, and RegionExposure.from_density
          all in one PR (§8)
Phase 5 — RegionExposure cleanup: remove A_r upper clip, update docstrings
          and Literal, rewrite test_exposure (§9)
Phase 6 — Config, CLI, docs, dead-script removal (§10)
```

`quant_from_buffer` stays blocked at its current `NotImplementedError`
throughout. EM behavior does not change in any phase.

## 3. Inventory

Identical to v3 §2 with two clarifications:

- §2.2 (move `estimate_strand_balance` to `strand_balance.py`,
  delete `EXON-INTRON`/`exon_contained`): unchanged.
- §2.3 (move `l_eff_contained` into `_exposure.py`): unchanged.
- §2.4: `CalibrationResult.global_densities` and `density_global.py` are both
  deleted at the end of Phase 4 (not retained as a summary stub).
- §2.5 (dead `scripts/debug/*`): unchanged.
- §2.6 (test inventory): unchanged.

## 4. Phase 0 — Cleanup

Identical to v3 §4. Reproduced here only for completeness:

1. Add `calibration/strand_balance.py`, move `StrandBalanceEstimate` and
   `estimate_strand_balance` from `density_global.py`, re-export from
   `calibration/__init__.py`.
2. Move `l_eff_contained` to `_exposure.py`. One-cycle re-export shim in
   `density_global.py`.
3. Trim `GlobalGdnaDensity` to `type`, `rho`, `n_fragments`, `eff_length_bp`,
   `n_regions_used`. Update `to_summary_dict`.
4. Remove `exon_intron` and `exon_contained` slots from `GlobalDensityTable`.
5. Audit imports.
6. Confirm `scripts/debug/*` deletion list with the owner.
7. `pytest tests/ -v` clean except for the pre-existing
   `test_calibration.py::TestStrandLLR::test_biased_toward_ss_favors_rna`.

Exit criteria: `density_global.py` contains only the trimmed
`GlobalGdnaDensity`/`GlobalDensityTable` and `compute_global_densities`.

## 5. Phase 1 — RegionCountLedger

Goal: one owner of derived totals. `RegionCountLedger` does **not** replicate
the per-channel arrays already exposed by `PayloadArrays`; it is a thin
companion that exposes:

```python
@dataclass(frozen=True, slots=True)
class RegionCountLedger:
    # Views into payload_arrays.region_counts_sorted (no copies).
    contained_unspliced_pos: np.ndarray
    contained_unspliced_neg: np.ndarray
    boundary_left_unspliced_pos: np.ndarray
    boundary_left_unspliced_neg: np.ndarray
    boundary_right_unspliced_pos: np.ndarray
    boundary_right_unspliced_neg: np.ndarray
    contained_spliced_pos: np.ndarray
    contained_spliced_neg: np.ndarray
    boundary_left_spliced_pos: np.ndarray
    boundary_left_spliced_neg: np.ndarray
    boundary_right_spliced_pos: np.ndarray
    boundary_right_spliced_neg: np.ndarray

    def contained_unspliced_total(self) -> np.ndarray   # float32[R]
    def boundary_left_unspliced_total(self) -> np.ndarray
    def boundary_right_unspliced_total(self) -> np.ndarray
    def boundary_unspliced_total(self) -> np.ndarray
    def contained_spliced_total(self) -> np.ndarray
    def boundary_spliced_total(self) -> np.ndarray
    def spliced_total(self) -> np.ndarray
```

Totals are returned as freshly-summed `float32[R]` arrays (the consumer is
free to cache). No `out=` parameter — premature optimisation.

Builder:

```python
def build_region_count_ledger(payload_arrays: PayloadArrays) -> RegionCountLedger
```

**PayloadArrays cleanup in this PR (not deferred):**

- Delete `contained_unspliced_total`, `boundary_left_unspliced_total`,
  `boundary_right_unspliced_total` fields from `PayloadArrays`.
- Update every consumer to use `ledger.contained_unspliced_total()` etc.
  Audit shows two consumers: `density_global._contained_density_for_mask`
  (about to be deleted in Phase 4) and `strand_deconv.estimate_kappa_d`
  (uses `payload_arrays.contained_unspliced_pos/neg`, which stay on
  `PayloadArrays` because they are direct channel slices, not totals).
- `PayloadArrays` keeps `region_counts_sorted` and the six per-channel
  unspliced contained/boundary POS/NEG views. Those are direct slices of
  `region_counts_sorted` and used by the strand path; they are not totals
  and do not duplicate ledger state.

Tests: `tests/test_region_count_ledger.py`.

- Channels indexed correctly for every signature.
- Totals equal direct slice sums bitwise.
- The ledger holds views, not copies (`np.shares_memory`).
- TS_NONE / TS_AMBIG rows retain observed counts.
- After removal, `PayloadArrays` exposes no total fields and no strand-path
  test regresses.

Exit criteria: orchestrator and density path both consume `RegionCountLedger`
totals; `PayloadArrays` no longer carries totals.

## 6. Phase 2 — DensityObservation

Transient, sorted, per-region table of gDNA-relevant counts and FL-aware
opportunities. **Not stored on `CalibrationResult`.**

```python
@dataclass(frozen=True, slots=True)
class DensityObservation:
    # counts (float32[R])
    contained_count: np.ndarray
    boundary_left_count: np.ndarray
    boundary_right_count: np.ndarray
    boundary_count: np.ndarray
    observed_compatible_count: np.ndarray  # contained + boundary

    # opportunities (float64[R])
    contained_leff: np.ndarray
    boundary_left_leff: np.ndarray
    boundary_right_leff: np.ndarray
    boundary_leff: np.ndarray

    # classification (bool[R])
    anchor_intergenic: np.ndarray
    anchor_intron: np.ndarray
    is_anchor: np.ndarray

    # diagnostics
    spliced_count: np.ndarray
    region_length: np.ndarray
```

Builder:

```python
def build_density_observation(
    region_arrays: RegionArrays,
    ledger: RegionCountLedger,
    gdna_fl: FragmentLengthModel,
) -> DensityObservation
```

Implementation is unchanged from v3 §6.

Tests: `tests/test_density_observation.py` (per v3 §6, unchanged).

Exit criteria: `DensityObservation` builds on every existing fixture and on
calibration scenario BAMs. No reference to it appears in
`CalibrationResult`.

## 7. Phase 3 — DensityModel

Goal: turn `DensityObservation` into `DensityEvidence` via robust MoM Gamma
prior, local update, and Negative Binomial predictive with exact stop-loss.

New file: `calibration/density_model.py`.

### 7.1 Dataclasses

`GammaRatePrior` is unchanged from v3 §7.1 except for one new field:

```python
    pearson_chi2_trimmed: float          # sum of (C - mu)^2 / mu over keep set
```

`DensityEvidence` is shrunk to the arrays that downstream code consumes:

```python
@dataclass(frozen=True, slots=True)
class DensityEvidence:
    # predictive (kept: read by exposure / locus prior consumers)
    rho_post: np.ndarray                  # float64[R]
    relative_exposure: np.ndarray         # float64[R] = rho_post / rho_ref
    mean_unbounded: np.ndarray            # float64[R]
    upper_unbounded: np.ndarray           # float64[R]

    # provenance and quality
    prior_family: np.ndarray              # uint8[R]: 0=INTERGENIC, 1=INTRON,
                                          # 2=ALL, 3=FALLBACK_BROAD,
                                          # 4=DETERMINISTIC_ZERO
    fallback_depth: np.ndarray            # uint8[R]
    flags: np.ndarray                     # uint8[R]
    confidence: float

    # priors and reference
    priors: dict[str, GammaRatePrior]
    rho_ref: float
    rho_ref_source: str                   # "ALL" | "WEIGHTED_FAMILIES" |
                                          # "ZERO"
```

Removed from the v3 evidence dataclass (computed on demand for summary, not
stored per-region):

- `variance_unbounded`, `precision_unbounded`
- `alpha_post`, `beta_post`
- `w_boundary_opportunity`, `w_boundary_count`
- `tail_probability`, `expected_tail_count`,
  `density_over_observed_ratio`

Summaries that need these (the §8.4 schema below) compute them once during
`to_summary_dict` from `(observation, priors, rho_post, rho_ref)`. This drops
the per-region float64 footprint from ~14 arrays to 4 + 3 uint8 metadata
arrays. At 10^7 regions this is ~360 MB instead of ~1.2 GB.

### 7.2 Background prior fit

Estimator signature is unchanged from v3 §7.2. Internals differ on the trim
selector and on the `rho_ref` / fallback definitions below.

#### 7.2.1 Eligibility

```text
eligible = (opportunities >= min_eff_length) & isfinite(counts)
if eligible.sum() < min_regions:
    raise InsufficientAnchors
```

#### 7.2.2 Pearson-trim MoM (v4 change)

```text
C    = counts[eligible].astype(float64)
L    = opportunities[eligible].astype(float64)
rho  = C.sum() / L.sum()
mu   = rho * L

# Standardised score: Pearson squared residual.
# Under the Gamma-Poisson null this is approximately scale-free in mu,
# so trimming by `score` does not preferentially drop high-mu rows.
score = (C - mu) ** 2 / np.maximum(mu, _SCORE_FLOOR)

k_keep = max(min_regions, int(np.ceil((1.0 - trim_upper) * score.size)))
keep   = np.argpartition(score, k_keep - 1)[:k_keep]

# Closed-form (S - B) / A identity computed on the same keep set so the
# variance estimator stays exposure-adjusted.
S_trim = np.sum((C[keep] - mu[keep]) ** 2)
B_trim = np.sum(mu[keep])
A_trim = np.sum(mu[keep] ** 2)

phi_hat   = max(0.0, (S_trim - B_trim) / max(A_trim, _PHI_EPS))
alpha_raw = 1.0 / max(phi_hat, phi_floor)
beta_raw  = alpha_raw / max(rho, _RHO_EPS)

beta  = min(beta_raw, beta_cap)
alpha = rho * beta
```

Constants:

```python
_SCORE_FLOOR = 1.0e-6   # protects (C - mu)^2 / mu when mu -> 0
_PHI_EPS     = 1.0e-12
_RHO_EPS     = 1.0e-12
```

Why Pearson scoring: `raw_squared_residual = (C - mu)^2` has expectation
`Var(C) = mu * (1 + mu * phi)` under Gamma–Poisson, which scales with `mu`.
Trimming the top quantile of raw residuals therefore selects high-`mu` rows
preferentially, removing variance the estimator was supposed to capture and
biasing `phi_hat` toward zero. Pearson scaling `(C - mu)^2 / mu` has
expectation `1 + mu * phi`, which only weakly increases in `mu`, so the trim
becomes near-uniform under the null and still removes nRNA-bearing intronic
outliers (whose `C/mu` ratio is large by construction).

Diagnostics recorded in `GammaRatePrior` (unchanged + one addition):
`n_trimmed`, `trim_upper`, `trimmed_mu_fraction`, `residual_sum`,
`poisson_variance_sum`, `extra_variance_basis_sum`, and the new
`pearson_chi2_trimmed = np.sum(score[keep])`.

#### 7.2.3 Beta cap (unchanged)

```python
def compute_beta_cap(rho: float) -> float:
    if rho <= 0.0 or not np.isfinite(rho):
        return float("inf")
    return float(DENSITY_PRIOR_MAX_PRECISION / rho)  # 400 / rho
```

#### 7.2.4 Reference density `rho_ref` (v4 change, locked)

`rho_ref` is computed once per `DensityEvidence` from the three family priors:

```python
def select_rho_ref(
    priors: dict[str, GammaRatePrior],
) -> tuple[float, str]:
    """Return (rho_ref, source)."""
    all_prior = priors.get("ALL")
    if all_prior is not None and all_prior.fit_status == "ok" and all_prior.mean_density > 0.0:
        return float(all_prior.mean_density), "ALL"

    weighted_num = 0.0
    weighted_den = 0.0
    for name in ("INTERGENIC", "INTRON"):
        p = priors.get(name)
        if p is None or p.fit_status != "ok" or p.mean_density <= 0.0:
            continue
        # Precision-weight: alpha is the Gamma shape, monotone in confidence.
        weighted_num += p.alpha * p.mean_density
        weighted_den += p.alpha
    if weighted_den > 0.0:
        return float(weighted_num / weighted_den), "WEIGHTED_FAMILIES"

    return 0.0, "ZERO"
```

#### 7.2.5 Deterministic zero-density path

When `rho_ref == 0.0` (no anchor family produced a valid fit) the entire
`DensityEvidence` is constructed deterministically:

```text
rho_post           = zeros(R)
relative_exposure  = ones(R)              # neutral exposure surface
mean_unbounded     = B_tot + C_r          # no prior shrinkage available
upper_unbounded    = mean_unbounded       # cannot quantify upper bound
prior_family       = DETERMINISTIC_ZERO   # uint8 = 4
fallback_depth     = 3
flags              = FLAG_FALLBACK_USED | FLAG_PRIOR_DOMINATED
rho_ref_source     = "ZERO"
priors             = {}                   # empty; summary writes null
```

This branch is fully covered by tests and is the only safe behaviour when no
contained gDNA evidence exists in the library.

### 7.3 Fallback selection (v4 change)

Per-region prior selection:

```text
depth 0: family prior for the anchor class (INTERGENIC or INTRON)
depth 1: all-anchor prior fit across INTERGENIC ∪ INTRON
depth 2: broad weak fallback = Gamma(alpha_fallback,
                                     alpha_fallback / rho_ref)
         with alpha_fallback = DENSITY_FALLBACK_PRIOR_ALPHA = 1.0
depth 3: deterministic-zero (only when rho_ref == 0; see §7.2.5)
```

The depth-2 prior is **weak by construction** — Gamma(1, 1 / rho_ref) has
CV = 1.0 and mean `rho_ref`. It is not capped against `beta_cap`; `beta_cap`
is only a cap on *fitted* priors (depths 0 and 1). Non-anchor regions (exon,
mixed) start at depth 1 because they have no native family.

### 7.4 Local update and predictive (v4 change)

Continuous posterior identities use unrounded `(B_tot, L_b, L_c)`:

```text
alpha_post = alpha_prior + B_tot
beta_post  = beta_prior  + L_b
rho_post   = alpha_post / beta_post

mean_c     = rho_post * L_c
p_nb       = beta_post / (beta_post + L_c)

mean_unbounded  = B_tot + mean_c
upper_unbounded = B_tot + scipy.stats.nbinom.ppf(confidence, alpha_post, p_nb)
```

Discrete CDF/SF evaluations use floored counts. This is the locked
fractional-count convention:

```text
c_int   = np.floor(C_r).astype(np.int64)
b_int   = np.floor(B_tot).astype(np.int64)

# nbinom in SciPy parameterises X = # failures before r successes,
# matching the Gamma-Poisson predictive with r = alpha_post,
# p = p_nb.
sf_c      = scipy.stats.nbinom.sf(c_int, alpha_post, p_nb)  # = P(X > c_int)
sf_cminus = scipy.stats.nbinom.sf(c_int - 1, alpha_post, p_nb)
tail_probability = sf_c
```

Closed-form stop-loss (v4 change). Replace v3 §7.4's truncated-sum +
normal-approximation with the NB stop-loss identity:

```text
# E[X * 1{X > c}] = mu * P(Y >= floor(c)) where Y ~ NB(r + 1, p).
# Equivalently P(Y >= k) = SF_{r+1}(k - 1).
sf_y = scipy.stats.nbinom.sf(c_int - 1, alpha_post + 1.0, p_nb)
expected_x_above_c = mean_c * sf_y
expected_tail_count = expected_x_above_c - c_int.astype(np.float64) * sf_c
```

No branching on mean magnitude. One additional SciPy call per region. The
identity is exact for any `(r, p, c)` and matches the discrete distribution
SciPy uses for `nbinom.cdf`, so the per-region tail summaries are internally
consistent.

For the summary-only diagnostics:

```text
density_over_observed_ratio = mean_unbounded / max(B_tot + C_r, eps)
var_c                = mean_c * (1.0 + L_c / beta_post)
cv                   = sqrt(var_c) / max(mean_c, eps)
precision_unbounded  = 1.0 / (1.0 + cv)
w_boundary_opportunity = L_b / (beta_prior + L_b)
w_boundary_count       = where(B_tot > 0, B_tot / (alpha_prior + B_tot), 0.0)
```

These are computed inside `DensityEvidence.to_summary_dict()` and reduced to
quantile statistics; they are not stored per region.

### 7.5 Flags

`FLAG_NON_ANCHOR`, `FLAG_FALLBACK_USED`, `FLAG_PRIOR_DOMINATED`,
`FLAG_LOW_CONTAINED_OPPORTUNITY`, `FLAG_LOW_BOUNDARY_OPPORTUNITY`,
`FLAG_HIGH_TAIL_TENSION` (rules unchanged from v3 §7.5).

`FLAG_HIGH_TAIL_TENSION` is evaluated against `tail_probability` computed
above and recorded as a bit in `flags`; the underlying float is not stored.

### 7.6 Public entry point

```python
def fit_density_evidence(
    observation: DensityObservation,
    *,
    confidence: float = 0.95,
    min_eff_length: float = DENSITY_MIN_EFF_LENGTH,
) -> DensityEvidence
```

Internal constants (additions to v3 §7.6):

```python
DENSITY_FALLBACK_PRIOR_ALPHA = 1.0   # Gamma shape for depth-2 fallback
```

All other constants are unchanged from v3 §7.6.

### 7.7 Tests

`tests/test_density_model.py`. Required cases (v4-specific deltas marked):

- `fit_gamma_prior` recovers mean density on Poisson anchors.
- `fit_gamma_prior` reports finite `phi` on overdispersed simulated anchors.
- `phi_hat == 0` falls back to `phi_floor` instead of `inf` alpha.
- **(v4) Pearson trim is unbiased under heteroscedastic `mu`:** simulate
  anchors with `L_i` log-uniform over three orders of magnitude and `C_i ~
  Poisson(rho * L_i)`. The Pearson-trimmed `phi_hat` is within tolerance of
  the untrimmed estimator. A raw-residual trim on the same data biases
  `phi_hat` upward by removing low-`mu` rows or downward by removing
  high-`mu` rows, depending on `trim_upper`. Assert the v4 estimator does
  not exhibit either bias.
- nRNA contamination (1–5% anchors injected with `C_i = 50 * mu_i`): trimmed
  `phi_hat` and `beta_raw` remain close to the Poisson fit; an untrimmed
  fitter collapses `beta_raw` toward zero.
- `compute_beta_cap(rho)` returns `400 / rho` for `rho > 0` and `inf` for
  `rho == 0`; depth invariant under common rescaling of `C` and `L`.
- `beta_raw > beta_cap` preserves `alpha / beta == rho` post-cap and yields
  prior CV `>= DENSITY_PRIOR_MIN_CV`.
- Sparse anchors trigger fallback depth 1 then depth 2 deterministically.
- **(v4) `select_rho_ref` priority order:** ALL ok → "ALL"; ALL fails, both
  families ok → precision-weighted; only one family ok → that family;
  nothing ok → `(0.0, "ZERO")`.
- **(v4) Depth-2 fallback prior has the documented shape:**
  `alpha == DENSITY_FALLBACK_PRIOR_ALPHA`, `beta == alpha / rho_ref`, mean
  `== rho_ref`, CV `== 1.0`.
- **(v4) Deterministic zero-density branch:** when no family fits,
  `relative_exposure` is all-ones, `rho_post` is all-zeros, `prior_family`
  is `DETERMINISTIC_ZERO`, every region carries `FLAG_FALLBACK_USED |
  FLAG_PRIOR_DOMINATED`, and `to_summary_dict()['rho_ref_source'] == "ZERO"`.
- Local update is correct: zero boundary leaves posterior equal to prior.
- High boundary count produces `rho_post > rho_ref`.
- NB predictive mean and ppf match `scipy.stats.nbinom` for integer cases.
- **(v4) Fractional-count convention:** non-integer `C_r` and `B_tot`
  produce the same `tail_probability` and `expected_tail_count` as their
  floored integer values, because the continuous posterior identities are
  used for `rho_post`/`mean_c` and the SF evaluations are explicitly floored.
- **(v4) Closed-form stop-loss matches truncated sum:** for small
  `(alpha_post, mean_c)` where direct enumeration is feasible, the
  closed-form `expected_tail_count` matches `sum_{k=c+1..K} (k - c) *
  P(X = k)` to within `1e-9 * max(1, mean_c)`.
- Tension flag fires only when `mean_unbounded` substantially exceeds
  observed.
- Flag bits set correctly on hand-crafted regions.

Exit criteria: `fit_density_evidence` finite on calibration scenarios.

## 8. Phase 4 — Result, Orchestrator, and Exposure Wiring

Goal: one PR that (a) makes `density_evidence` a first-class field on
`CalibrationResult`, (b) wires the orchestrator to build it, (c) constructs
`RegionExposure` from it. v3's Phase 4/5 split is collapsed because Phase 4
already calls `from_density`.

### 8.1 `_result.py`

```python
@dataclass(frozen=True, slots=True)
class CalibrationResult:
    fl_models: FLModels
    diagnostics: Diagnostics
    region_gdna: "RegionGdnaEstimate"
    region_exposure: "RegionExposure"
    density_evidence: "DensityEvidence"
    n_multi_loci: int = 0
    rna_lower_confidence: float = 0.95
```

Removed: `global_densities`, `density_observation`. `DensityObservation` is
transient; it is built inside the orchestrator, consumed by
`fit_density_evidence`, and dropped at function exit. Nothing downstream of
`CalibrationResult` needs the raw observation table; all relevant counts and
opportunities are summarised inside `DensityEvidence` or recomputable from
the original `payload_arrays` + `ledger`.

`to_summary_dict()` drops the `global_densities` block; it emits the
`density_evidence` block from §8.4.

`build_calibration_result` requires `density_evidence` and `region_exposure`.

### 8.2 `_orchestrator.py`

```text
1. Build RegionArrays and PayloadArrays.
2. Build RegionCountLedger.
3. Build FL models.
4. Run strand path (independent).
5. observation = build_density_observation(region_arrays, ledger,
                                           fl_models.gdna)
6. density_evidence = fit_density_evidence(
       observation,
       confidence=config.gdna_density_confidence,
   )
7. region_exposure = RegionExposure.from_density(
       density_evidence,
       max_exposure=config.density_max_exposure,
   )
8. Assemble CalibrationResult (without observation).
```

`compute_global_densities` is no longer called. `density_global.py` is
deleted in this PR.

### 8.3 `RegionExposure.from_density` (lives in `exposure.py`, lands here)

```python
@classmethod
def from_density(
    cls,
    density_evidence: "DensityEvidence",
    *,
    max_exposure: float | None = None,
) -> "RegionExposure":
    A_r       = density_evidence.relative_exposure.astype(np.float32, copy=True)
    rho_r     = density_evidence.rho_post.astype(np.float32, copy=True)
    eligible  = (density_evidence.flags & FLAG_LOW_BOUNDARY_OPPORTUNITY) == 0
    if max_exposure is not None:
        np.minimum(A_r, np.float32(max_exposure), out=A_r)
    return cls(
        mode="density",
        A_r=A_r,
        rho_r=rho_r,
        rho_ref=float(density_evidence.rho_ref),
        reference_quantile=0.0,
        eligible=eligible,
        flags=np.asarray(density_evidence.flags, dtype=np.uint8),
    )
```

**Ineligible-region rule (locked):** ineligible regions keep their posterior
`A_r` and `rho_r` values; `eligible[i] = False` is the consumer's signal not
to trust them. We do not overwrite ineligible `A_r` to 1.0 — doing so would
silently bias downstream consumers that ignore the mask.

The `mode` Literal becomes `Literal["uniform", "density"]`. `"unsupervised"`
is removed (it was never produced).

### 8.4 Summary schema

```text
density_evidence:
  confidence: float
  n_regions: int
  n_anchor_intergenic: int
  n_anchor_intron: int
  rho_ref: float
  rho_ref_source: "ALL" | "WEIGHTED_FAMILIES" | "ZERO"
  priors:
    INTERGENIC: { alpha, beta, mean_density, phi, beta_raw, beta_cap,
                  cap_applied, n_regions, fit_status,
                  pearson_chi2_trimmed, trimmed_mu_fraction }
    INTRON:     { ... }
    ALL:        { ... } | null
  fallback_depth_histogram: { 0: int, 1: int, 2: int, 3: int }
  prior_family_histogram:   { INTERGENIC, INTRON, ALL, FALLBACK_BROAD,
                              DETERMINISTIC_ZERO }
  flags:
    n_low_contained_opportunity: int
    n_low_boundary_opportunity: int
    n_prior_dominated: int
    n_high_tail_tension: int
    n_fallback_used: int
  relative_exposure: { min, p50, p95, max, mean }
  rho_post:           { p50, p95, max }
  mean_unbounded:     { p50, p95, max }
  tail_probability:   { p50, p95, max }   # computed on-demand
```

No per-region arrays serialised.

### 8.5 Downstream consumers (locked actions)

- `src/rigel/sim/locus_sweep.py` — replace
  `cal.global_densities.intergenic.lambda_gdna` with
  `cal.density_evidence.priors["INTERGENIC"].mean_density`. Same column,
  renamed.
- `src/rigel/sim/analysis.py` — the `EXON-INTRON` rho column maps to
  `cal.density_evidence.priors["INTRON"].mean_density` if INTRON is the only
  available family; drop the metric otherwise. Document in the diff.
- `tests/golden/*` regenerate with `pytest --update-golden`.

### 8.6 Tests

- `tests/test_calibrate.py` — assert `density_evidence` and
  `region_exposure` are populated; strand-only changes do not move density
  evidence; ineligible regions retain `A_r` from posterior.
- `tests/test_region_exposure_from_density.py` — alignment with
  `density_evidence.relative_exposure`, `max_exposure` clipping behaviour,
  eligibility wiring (no overwrite), `mode == "density"`.

Exit criteria:

- `density_global.py` deleted.
- `compute_global_densities` gone.
- `CalibrationResult.global_densities` and `.density_observation` do not
  exist.
- All non-debug tests pass.

## 9. Phase 5 — RegionExposure Cleanup and A_r Audit

Goal: remove the lingering "A_r in (0, 1]" assumptions now that density-derived
exposure can exceed 1.

Concrete code edits:

1. `src/rigel/calibration/_exposure.py` — in
   `bp_weighted_mean_exposure_over_blocks`, replace
   `np.clip(weighted_bp / raw_bp, min_weight, 1.0)` with
   `max(weighted_bp / raw_bp, min_weight)`. Update the docstring to state
   that returned weights are non-negative bp-weighted means and may exceed 1
   when the exposure surface is density-derived. The `min_weight` floor
   remains.
2. `src/rigel/calibration/exposure.py` — update `RegionExposure` docstring
   and `A_r` field comment to read `float32[R], non-negative; may exceed 1
   when constructed from density evidence`. Update the `mode` Literal to
   `Literal["uniform", "density"]`.
3. `RegionExposure.to_summary_dict` — emit `A_p99` in addition to the
   existing min/mean/max so the summary surfaces the upper tail without a
   full histogram.
4. Audit consumers of `RegionExposure.A_r` (grep `src/`). Today the only
   non-test reader is `bp_weighted_mean_exposure_over_blocks` itself. Add a
   `# A_r may exceed 1` comment at every read site discovered by the audit
   to make the contract change unmissable.
5. Remove `tests/test_exposure.py` from `tests/conftest.py::collect_ignore`
   and rewrite it: cover the uniform constructor, the density constructor,
   ineligible-row behaviour, and `max_exposure` clipping.

Exit criteria:

- No code path clips `A_r` to 1.
- `RegionExposure.from_density` is the orchestrator's only exposure source
  (the `uniform` constructor remains as a fallback for tests).
- `test_exposure.py` runs in CI.

## 10. Phase 6 — Config, CLI, Docs

### 10.1 Config

`src/rigel/config.py` `CalibrationConfig` adds:

```python
gdna_density_confidence: float = 0.95
density_min_eff_length: float = 1.0
density_max_exposure: float | None = None
```

Validation: `0.5 <= gdna_density_confidence < 1.0`; `density_min_eff_length
>= 0.0`; `density_max_exposure is None or density_max_exposure > 0.0`.

### 10.2 CLI

```text
--gdna-density-confidence FLOAT (default 0.95)
--density-max-exposure   FLOAT (default unset; no clip)
```

Wire via `_ParamSpec`. `density_min_eff_length` stays internal for now.

### 10.3 Docs

- `docs/MANUAL.md` — short subsection on the new CLI flags and the meaning
  of `A_r > 1`.
- `docs/fineregions/calibration_roadmap_v3.md` — mark Phase 1–4 (this plan)
  as complete on landing.
- Move `density_model_impl_plan_v1.md`, `_v2.md`, `_v3.md` into
  `docs/fineregions/archive/`.

## 11. Removal Checklist (after Phase 6)

- `calibration/density_global.py`
- `GlobalGdnaDensity`, `GlobalDensityTable`, `DensityType`,
  `_contained_density_for_mask`, `_empty_density`
- `CalibrationResult.global_densities`
- `CalibrationResult.density_observation` (never lands)
- `DENSITY_PRIOR_EQUIV_BOUNDARIES` (never landed)
- `to_summary_dict()` keys `"global_densities"`, `"EXON-INTRON"`,
  `"exon_contained"`
- Legacy fields `lambda_gdna`, `strand_active`, `rho_uncorrected`,
  `strand_correction_fragments`, `n_strand_informative_regions`,
  `strand_informative_exposure_fraction`, `n_uninf_observed`,
  `n_fragments_estimated`, `n_rows_eligible`
- `PayloadArrays.contained_unspliced_total`,
  `PayloadArrays.boundary_left_unspliced_total`,
  `PayloadArrays.boundary_right_unspliced_total`
- `mode == "unsupervised"` on `RegionExposure`
- Upper clip in `bp_weighted_mean_exposure_over_blocks`
- The `scripts/debug/*` files listed in v3 §2.5

Verification:

```bash
grep -rEn "global_densities|EXON-INTRON|exon_contained|lambda_gdna|\
contained_unspliced_total" src/ tests/
```

Zero non-archive hits.

## 12. Risks and Unexpected Aspects

1. Golden outputs change. Mitigation: explicit `--update-golden` in the
   Phase 4 PR with reviewer-visible diff.
2. **Memory (v4 update).** The shrunk `DensityEvidence` carries 4 float64
   per-region arrays plus 3 uint8 metadata arrays. At 10^7 regions that is
   ~360 MB. `DensityObservation` (transient, larger) is freed before
   `CalibrationResult` is returned. If 360 MB is still too large, switch
   `mean_unbounded` and `upper_unbounded` to float32 after measurement.
3. NB ppf cost for very large `R`. Vectorize `scipy.stats.nbinom.ppf`; add a
   saddlepoint or normal approximation path keyed off `(alpha_post, L_c /
   beta_post)` if it becomes a hotspot. The closed-form stop-loss removes
   one branchy fallback path that v3 had.
4. **`RegionExposure.A_r > 1` (v4 mitigation).** Phase 5 explicitly removes
   the upper clip and audits consumers; the contract change is no longer
   silent.
5. Beta cap behaviour varies between datasets. Cap is a CV floor (depth
   invariant); `cap_applied` exposed in prior diagnostics.
6. Intronic anchors contaminated by nascent RNA. **(v4 strengthening)** the
   Pearson-scored trim is robust to heteroscedastic `mu` so the
   contamination test on a long, high-`mu` intron is no longer trimmed
   preferentially. Per-family diagnostics still record trimmed mass.
7. Capture-mode boundary leakage into INTRON-only anchors. Same as v3:
   `beta_cap` keeps the prior movable; per-family `n_fragments` and
   `eff_length` are recorded.
8. **(v4) Fractional counts.** Continuous posterior identities consume raw
   floats. Discrete CDF/SF evaluations always floor; this is tested. Bias
   per region is bounded by `P(X = floor(C_r))`, negligible for any region
   with non-trivial counts.

## 13. Things v4 Intentionally Does Not Do

Unchanged from v3 §13.

## 14. Done Definition

v4 is done when:

- All six phases are merged.
- `CalibrationResult.density_evidence` populated for every calibration run.
- `RegionExposure.from_density` is the orchestrator's exposure source; the
  upper-clip is gone.
- `density_global.py` and listed legacy fields/scripts no longer exist.
- `PayloadArrays` no longer carries totals.
- `pytest tests/ -v` passes except the pre-existing
  `tests/test_calibration.py::TestStrandLLR::test_biased_toward_ss_favors_rna`.
