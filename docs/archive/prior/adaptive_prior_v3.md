# Adaptive Prior v3 — Quantile/ESS Grouped Prior (Implementation Spec)

Date: 2026-05-26
Status: implementation-ready
Supersedes: `adaptive_prior_v2.md`, `adaptive_prior_v1.md`, `one_knob_parameterization.md`

## 1. Scope

This document specifies the grouped RNA/gDNA EM prior policy. It is a
complete, self-contained build spec: data structures, math, file layout,
phase order, and acceptance tests. No prior reading is required to
implement it.

It does **not** change the EM solver, the calibration pipeline structure,
the strand model, or the unspliced-only contract from prior redesign v3.
It changes only how grouped `(alpha_gdna_add, alpha_rna_add)` are derived
from calibration outputs and how that derivation is exposed to users.

## 2. The contract in one paragraph

For every locus `l`, calibration produces a posterior over the gDNA share
of unspliced mass `phi_l = D_l / U_l`. The user supplies one scalar
`rna_confidence = q ∈ (0, 1)` expressing risk preference. v3 returns a
grouped Beta-equivalent prior whose **share** is `Q_q(phi_l)` (the q-th
posterior quantile) and whose **effective sample size** is the posterior
ESS implied by the posterior's mean and variance, capped by a single
safety constant. No other user knobs control prior strength, edge mass,
or logit bias.

## 3. User surface

### 3.1 CLI

Add one option:

```text
--rna-confidence FLOAT     Posterior quantile of the gDNA share used as the
                           grouped EM prior share. Range (0, 1). Default 0.5.
                           Higher = more conservative RNA reporting.
```

### 3.2 Deprecation

For exactly one release, the old flags are accepted as **identity
aliases** with a `DeprecationWarning`:

| Old flag                      | New flag             | Mapping       |
| ----------------------------- | -------------------- | ------------- |
| `--rna-lower-confidence Q`    | `--rna-confidence Q` | identity      |
| `--gdna-density-confidence Q` | `--rna-confidence Q` | identity      |

Identity is correct: `rna_lower_confidence = 0.95` previously meant
"use the 95% RNA lower credible bound", i.e. the 95% gDNA upper credible
bound, i.e. `Q_{0.95}(phi)`. Both old defaults were `0.95`. The
unified default is changed to `0.5` (Bayes-optimal under symmetric loss);
old configs that set explicit `0.95` keep getting `q = 0.95` with a
warning. Conflicting flag values (differing by more than `1e-12`) fail
fast.

### 3.3 Removed from the user surface

Any value set in YAML or CLI for the following is **rejected** at config
load with an actionable error pointing at `--rna-confidence`:

```text
aggregate_prior_strength
aggregate_prior_edge_count
aggregate_prior_max_count
gdna_prior_logit_bias
```

These names survive only as private safety constants in §6.7.

## 4. Inputs from calibration (per locus `l`)

| Symbol                       | Source                                 | Meaning                                                  |
| ---------------------------- | -------------------------------------- | -------------------------------------------------------- |
| `U_l`                        | `PriorMassDeconvolution.unspliced_total` | total prior-owned unspliced mass                       |
| `mu_D_l`                     | `PriorMassDeconvolution.gdna_unspliced_mean` | posterior mean of gDNA-owned unspliced count       |
| `var_D_l`                    | `PriorMassDeconvolution.gdna_unspliced_var` (NEW) | posterior variance of the same count          |
| `has_gdna_candidate_l`       | native structural flag                 | gates whether the grouped prior may be nonzero           |
| `prior_allocated_fraction_l` | existing projection diagnostic         | share of `U_l` actually allocated to a locus             |

The RNA-side moments are determined by mass conservation
(`mu_R = U - mu_D`, `var_R = var_D` for fixed observed `U`); they are not
stored separately.

All other v3 fields (`mean`, `precision`, `method`, `flags`) remain. No
fields are removed.

## 5. The algorithm

### 5.1 Step 1 — local Beta posterior over `phi_l`

Use moment matching on the bounded ratio `phi = D / U`:

```text
p_l   = clamp(mu_D_l / U_l, 0, 1)                       if U_l > eps else 0
v_max = p_l * (1 - p_l)
v_l   = clamp(var_D_l / U_l^2, 0, v_max - EPS_BETA)
tau_l = p_l * (1 - p_l) / v_l - 1                       if v_l > 0 else 0
tau_l = clamp(tau_l, 0, MAX_ESS)
a_local_l = tau_l * p_l
b_local_l = tau_l * (1 - p_l)
```

Boundary handling:

- `U_l <= eps` → `a_local = b_local = 0`, locus is purely shrinkage-driven.
- `var_D_l` missing → set `tau_l = 0`, mark `FLAG_PRIOR_NO_LOCAL_VARIANCE`,
  rely entirely on global shrinkage (§5.2).
- `p_l ∈ {0, 1}` with `tau_l > 0` is allowed and means "calibration is
  near-certain"; the quantile in §5.3 will collapse to that boundary.

`EPS_BETA = 1e-9`. This is a numerical guard for `betaincinv`, not a
modeling parameter.

### 5.2 Step 2 — empirical-Bayes global Beta

Across all loci with `U_l > 0` and a finite `var_D_l`:

```text
p_i = mu_D_i / U_i
v_i = var_D_i / U_i^2
w_i = 1 / max(v_i, EPS_W)                       (inverse-variance weight)
w_i = min(w_i, U_i + 1)                         (cap by data content)

psi   = sum(w_i * p_i) / sum(w_i)
v_obs = sum(w_i * (p_i - psi)^2) / sum(w_i)
v_mea = sum(w_i * v_i)         / sum(w_i)
v_bet = max(v_obs - v_mea, EPS_BETA)

kappa = clamp( psi * (1 - psi) / v_bet - 1, 0, MAX_ESS )
```

If `sum(w_i) == 0` (no usable loci): set `psi = 0`, `kappa = 0`, record
`FLAG_GLOBAL_PROFILE_EMPTY`. No magic `+0.5` Laplace term: an empty
profile must not silently invent a global gDNA rate.

The global Beta is `Beta(kappa * psi, kappa * (1 - psi))`.

### 5.3 Step 3 — combined posterior and decision summary

```text
a_post_l = a_local_l + kappa * psi
b_post_l = b_local_l + kappa * (1 - psi)

if a_post_l + b_post_l > 0 and has_gdna_candidate_l:
    share_l = Q_q( Beta(a_post_l, b_post_l) )
    ess_l   = min(a_post_l + b_post_l, MAX_ESS)
else:
    share_l = 0
    ess_l   = 0

alpha_gdna_add_l = ess_l * share_l
alpha_rna_add_l  = ess_l * (1 - share_l)
```

`Q_q` is `scipy.special.betaincinv(a, b, q)`. Degenerate `Beta(0, b)` or
`Beta(a, 0)` short-circuit to `share = 0` or `share = 1` without calling
the special function.

### 5.4 What this design deliberately omits

- **No reliability multipliers.** The v2 chain
  `structural * prior_mass * allocation * convergence * variance_completeness`
  is replaced by: (a) `has_gdna_candidate_l` as a hard gate, (b) variance
  itself shrinks `tau_l` automatically when projection or convergence is
  uncertain. Producing low ESS is the *job* of variance; we do not
  multiply variance-derived ESS by hand-picked constants.
- **No variance inflation formulas.** Projection ambiguity is handled in
  §6.2 by the variance projection rule directly. Independence-violating
  cases are reflected as larger `var_D_l`, which is the honest channel.
- **No conflict gating.** Prior-vs-EM conflict is recorded as a
  diagnostic (`prior_conflict_score`, §7) but never used to override the
  prior — that would re-introduce hidden tuning.
- **No alpha floors.** Any numerical floor used to evaluate `betaincinv`
  must not leak into the returned `(alpha_gdna_add, alpha_rna_add)`.

## 6. Required changes to existing data flow

### 6.1 Calibration must emit posterior variance

#### `RegionGdnaEstimate` (`src/rigel/calibration/_result.py`)

Add one field:

```python
var_count: np.ndarray   # Var(D_r | strand evidence); shape == mean_count.shape
```

Population (both strand paths already compute the same SD):

```text
var_count = sd_r ** 2
```

Degenerate cases: `var_count = 0` when `n_total = 0`; `var_count = NaN`
when the strand path is structurally uninformative (consumers treat NaN
as "use density fallback").

#### `RegionGdnaChannelEstimate`

Add three fields:

```python
contained_var:      np.ndarray
boundary_left_var:  np.ndarray
boundary_right_var: np.ndarray
```

The aggregate `Var(D_r)` is the **sum** of the three, treating the
compartments as approximately independent. This approximation is recorded
in `summary.json` under `prior_diagnostics.compartment_independence_assumed`.

#### `PriorMassDeconvolution`

Add one field:

```python
gdna_unspliced_var: np.ndarray
```

Invariants:

```text
gdna_unspliced_mean + rna_unspliced_mean == unspliced_total
gdna_unspliced_var  >= 0  and  finite (or marked via FLAG_PRIOR_NO_LOCAL_VARIANCE)
gdna_unspliced_var  <= gdna_unspliced_mean * (unspliced_total - gdna_unspliced_mean)
```

The upper bound is enforced by clamping at construction (Bernoulli-like
ceiling for a bounded count); a clamp event sets `FLAG_PRIOR_VARIANCE_CLAMPED`.

### 6.2 Variance projection through `_allocate_counts_by_geometry`

The existing allocator computes overlap shares `s_lr` from region `r` to
locus `l` and sums:

```text
U_l    = sum_r s_lr * U_r
mu_D_l = sum_r s_lr * mu_D_r
```

Add the variance projection alongside, **using squared shares**:

```text
var_D_l = sum_r s_lr^2 * var_D_r
```

This is the exact rule for independent regions. When a single region is
split across multiple loci, the per-locus contributions are bounded
above by the region's own variance (since `sum s^2 <= sum s = 1`); when
multiple disjoint regions contribute, variance adds. Both are correct
limits for the modeled assumption.

No extra inflation term. If overlap-induced dependence is a future
concern, model it explicitly in calibration; do not patch it here.

### 6.3 Density-fallback variance

When the density path is used instead of strand deconvolution, the
posterior predictive count is Gamma–Poisson:

```text
rho ~ Gamma(alpha, beta);   D | rho ~ Poisson(rho * L)
E[D]   = alpha * L / beta
Var[D] = E[D] + E[D]^2 / alpha
```

For the mixture used in the current density fallback
(`p_captured * captured + (1 - p_captured) * off_target`):

```text
mu_D    = p * mu_cap + (1 - p) * mu_off
var_D   = p * var_cap + (1 - p) * var_off
        + p * (1 - p) * (mu_cap - mu_off)^2
```

with `alpha, beta` for each component taken from
`background.rho_off_(alpha|beta)` and the sweep excess.

If neither variance source is available for a region (legacy artefact),
the region is excluded from `PriorMassDeconvolution.gdna_unspliced_var`
by setting that entry to `NaN` and flagging
`FLAG_PRIOR_VARIANCE_FROM_INTERVAL`. The EB step in §5.2 already excludes
non-finite variances. No `z_interval` heuristic is introduced.

### 6.4 Config (`src/rigel/config.py`)

Add to `CalibrationConfig`:

```python
rna_confidence: float = 0.5
```

Validation: `0.0 < rna_confidence < 1.0`.

Keep `rna_lower_confidence` and `gdna_density_confidence` as `Optional[float] = None`
**only** to allow the resolver in §6.5 to detect explicit user input.
Their docstrings mark them deprecated. Setting them simultaneously with
different non-None values raises `ValueError`.

Remove from `CalibrationConfig`:

```text
aggregate_prior_strength
aggregate_prior_edge_count    (if present)
aggregate_prior_max_count
gdna_prior_logit_bias
```

YAML loaders that encounter these keys raise:

```text
Configuration option '{name}' was removed in adaptive prior v3.
Use --rna-confidence (or calibration.rna_confidence) instead.
See docs/prior/adaptive_prior_v3.md.
```

### 6.5 CLI (`src/rigel/cli.py`)

Add:

```python
parser.add_argument(
    "--rna-confidence", dest="rna_confidence", type=float, default=None,
    help="Posterior quantile of gDNA share used as the grouped EM prior "
         "share. Range (0, 1). Higher = more conservative RNA reporting. "
         "Default 0.5.",
)
```

Resolver (run after parsing, before constructing `CalibrationConfig`):

```python
def resolve_rna_confidence(args) -> float:
    candidates = {
        name: getattr(args, name)
        for name in ("rna_confidence", "rna_lower_confidence",
                     "gdna_density_confidence")
        if getattr(args, name, None) is not None
    }
    if not candidates:
        return 0.5
    values = list(candidates.values())
    if max(values) - min(values) > 1e-12:
        raise SystemExit(
            f"Conflicting confidence flags: {candidates}. "
            "Use a single --rna-confidence."
        )
    for name in ("rna_lower_confidence", "gdna_density_confidence"):
        if name in candidates:
            warnings.warn(
                f"--{name.replace('_','-')} is deprecated; "
                "use --rna-confidence.",
                DeprecationWarning,
                stacklevel=2,
            )
    return next(iter(candidates.values()))
```

The existing `_ParamSpec` table at `cli.py:592-593` is updated to drop
the two legacy entries and add `rna_confidence`.

### 6.6 New helper module

```text
src/rigel/calibration/adaptive_prior.py
```

Public functions (all pure, no I/O, NumPy-vectorized):

```python
def local_beta_from_moments(
    mu_D: np.ndarray, var_D: np.ndarray, U: np.ndarray,
    *, max_ess: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (a_local, b_local). Handles U==0, missing variance, clamping."""

def global_beta_from_profile(
    mu_D: np.ndarray, var_D: np.ndarray, U: np.ndarray,
    *, max_ess: float,
) -> tuple[float, float, dict]:
    """Return (psi, kappa, diagnostics)."""

def compute_grouped_prior_counts(
    decon: PriorMassDeconvolution,
    has_gdna_candidate: np.ndarray,
    *, rna_confidence: float, max_ess: float = MAX_ESS,
) -> AdaptivePriorResult:
    """Return arrays alpha_gdna_add, alpha_rna_add, ess, share, plus
    diagnostic columns listed in §7."""
```

`AdaptivePriorResult` is a `@dataclass(frozen=True, slots=True)` carrying
exactly the columns listed in §7; consumers select what they need.

### 6.7 Internal safety constants

In `adaptive_prior.py`:

```python
MAX_ESS    = 3000.0      # cap on tau_l and kappa
EPS_BETA   = 1.0e-9      # numerical guard for betaincinv
EPS_W      = 1.0e-12     # guard for 1/var weighting
```

These are not exposed via CLI or YAML. Changes require a code review and
a regression run of the adaptive-prior sentinel suite (§9). They are
documented inline as safety rails, not model knobs.

## 7. Outputs

### 7.1 `loci.feather` — add columns

```text
prior_gdna_share_mean       # p_local_l
prior_gdna_share_var        # v_l
prior_gdna_share_quantile   # share_l (= Q_q of posterior)
prior_ess_local             # tau_l
prior_ess_global            # kappa contribution = kappa
prior_ess_final             # ess_l
alpha_gdna_add              # final prior alpha (gDNA group)
alpha_rna_add               # final prior alpha (RNA group)
prior_flags                 # bitfield (see §7.3)
prior_conflict_score        # filled post-EM (NaN until EM completes)
```

Legacy columns `prior_budget` and `prior_gdna_share_biased` are removed.
This is a one-shot golden regeneration, not an alias migration.

### 7.2 `summary.json` — add block

```json
"prior_policy": {
  "name": "adaptive_quantile_v3",
  "rna_confidence": 0.5,
  "global_gdna_share": 0.04,
  "global_gdna_ess":   17.3,
  "n_loci_total":      12345,
  "n_loci_with_local_variance": 11892,
  "n_loci_structural_gated":     410,
  "max_ess_cap":       3000.0,
  "flag_histogram":    {"PRIOR_NO_LOCAL_VARIANCE": 453, ...}
}
```

### 7.3 `prior_flags` bitfield

```text
0x1  PRIOR_NO_LOCAL_VARIANCE        var_D_l missing; relied on global only
0x2  PRIOR_VARIANCE_CLAMPED         var_D_l clamped to bounded ceiling
0x4  PRIOR_VARIANCE_FROM_INTERVAL   density region lacked variance (excluded)
0x8  PRIOR_STRUCTURAL_GATED         has_gdna_candidate_l == False
0x10 GLOBAL_PROFILE_EMPTY           kappa = 0 (no usable loci)
```

### 7.4 Post-EM conflict diagnostic (filled by EM driver)

```text
em_share_l            = gdna_em_count_l / max(total_em_count_l, eps)
prior_conflict_score  = | logit_clip(em_share_l, 1e-6) -
                          logit_clip(prior_gdna_share_quantile_l, 1e-6) |
```

Diagnostic only. Not used to alter priors.

## 8. Phase plan

Each phase is a self-contained PR with its own green test set. No phase
depends on tests defined in a later phase.

### Phase A — helper module + unit tests (no integration)

Add `src/rigel/calibration/adaptive_prior.py` and
`tests/test_adaptive_prior_helpers.py`. Tests:

1. `test_local_beta_monotone_in_variance` — larger `var_D` → smaller `tau`.
2. `test_local_beta_zero_total` — `U=0` ⇒ `a=b=0`.
3. `test_local_beta_missing_variance` — `NaN var_D` ⇒ `a=b=0`, no exception.
4. `test_local_beta_ess_cap` — `tau <= MAX_ESS`.
5. `test_global_beta_homogeneous` — all `p_i = c` ⇒ `psi = c`, large `kappa`.
6. `test_global_beta_heterogeneous` — disagreement ⇒ small `kappa`.
7. `test_global_beta_empty` — no usable loci ⇒ `psi=0, kappa=0, flag set`.
8. `test_quantile_monotone_in_q` — fixed posterior, `Q_{0.1} <= Q_{0.5} <= Q_{0.9}`.
9. `test_quantile_one_sided_high_confidence` — `p=0.01`, `tau=200`, `q=0.9` ⇒ `share` small.
10. `test_quantile_one_sided_low_confidence` — `p=0.01`, `tau=2`,  `q=0.9` ⇒ `share` larger than (9).
11. `test_structural_gate_zeroes_alphas` — `has_gdna_candidate=False` ⇒ both alphas exactly 0.

### Phase B — propagate variance through calibration

- Add `var_count` to `RegionGdnaEstimate`; populate from both strand
  paths and density paths.
- Add `(contained|boundary_left|boundary_right)_var` to
  `RegionGdnaChannelEstimate`.
- Add `gdna_unspliced_var` to `PriorMassDeconvolution`; populate via
  squared-share projection (§6.2) and density-fallback formula (§6.3).
- Tests: `test_calibration_variance_propagation.py` covers exact strand
  variance vs discrete posterior, normal-path SD round-trip, Gamma–Poisson
  mixture variance, and squared-share projection (disjoint and split
  regions).

### Phase C — wire helper into `assemble_priors`

- Replace v3 grouped-prior construction with
  `compute_grouped_prior_counts`. Keep the call site signature unchanged
  for downstream code; only the body changes.
- Add the new `loci.feather` columns and `summary.json` block (§7).
- Delete the old grouped-prior code path and constants. No `legacy_v3`
  toggle is added; the only way to run the previous policy is to check
  out the previous tag.

### Phase D — CLI + config migration

- Add `--rna-confidence` and resolver (§6.5).
- Delete `aggregate_prior_*` and `gdna_prior_logit_bias` from
  `CalibrationConfig` and from the YAML loader's known-key set; raise on
  unknown keys with the migration error message (§6.4).
- Update `_ParamSpec` in `cli.py`.

### Phase E — golden regeneration + sentinel sweep

```bash
conda activate rigel && pytest tests/ -v
conda activate rigel && pytest tests/ --update-golden
```

Required monotonicity sweep (new test `test_adaptive_prior_q_sweep.py`):
run the synthetic 24-condition harness at `q ∈ {0.1, 0.5, 0.9}`. For
every locus with `ess_final > 0`:

```text
share(q=0.1) <= share(q=0.5) <= share(q=0.9)
alpha_gdna_add monotone non-decreasing in q
alpha_rna_add  monotone non-increasing in q
```

The cap and EB constants are not touched between sweeps; failure means
the sweep, not the constants, is wrong.

## 9. Acceptance criteria

The implementation is accepted when **all** of the following hold:

1. `pytest tests/ -v` is green (modulo the pre-existing known failure
   `test_calibration.py::TestStrandLLR::test_biased_toward_ss_favors_rna`).
2. `--rna-confidence` is the only prior-related option in
   `rigel quant --help`.
3. Setting any removed key in YAML produces the migration error from §6.4.
4. The q-sweep monotonicity test from Phase E passes without touching
   `MAX_ESS`, `EPS_BETA`, or `EPS_W`.
5. Zero-gDNA sentinels (`gdna_none_ss_*`) recover RNA without any logit
   bias being set anywhere in code or config.
6. High-gDNA sentinels (`gdna_high_ss_*`) do not see `alpha_gdna_add`
   collapse to zero in loci with a structural gDNA candidate.
7. `loci.feather` carries the columns in §7.1; `summary.json` carries
   the block in §7.2.
8. No file outside `src/rigel/calibration/adaptive_prior.py` and the
   updated config/CLI/output code references the removed v3 constants.

## 10. Non-goals

- Replacing the EM solver's prior consumption interface.
- Changing the unspliced-only contract.
- Auto-tuning `rna_confidence` from the data. The whole point of the
  knob is that it is the one quantity the data cannot supply.
- Converting calibration self-training damping to ELBO line-search.
  That is a separate workstream (`docs/newcal/`); the present spec does
  not need it and does not block on it.

## 11. Summary

```text
calibration  ->  (mu_D_l, var_D_l, U_l) per locus
             ->  local Beta via moment matching
             ->  empirical-Bayes global Beta from inverse-variance pooling
             ->  combined Beta posterior over phi_l
             ->  share_l = Q_q(posterior)
             ->  ess_l   = min(a + b, MAX_ESS) * 1{has_gdna_candidate}
             ->  (alpha_gdna_add, alpha_rna_add) = ess_l * (share_l, 1 - share_l)
```

One user knob. One safety cap. One algorithm. No reliability multipliers,
no inflation terms, no logit bias, no edge counts, no strength
constants. Variance does the work that magic constants used to do.
