# Critique & Implementation Spec for `new_cal_plan_v5.md`

Date: 2026-05-26
Audience: implementer of the v5 calibration redesign.
Companion to: [new_cal_plan_v5.md](new_cal_plan_v5.md)

This document is in two parts:

1. **Critique** — concrete issues that block implementation as written.
2. **Implementation Spec** — phase-by-phase, module-by-module description
   that can be translated directly into code on top of the current
   `src/rigel/calibration/` tree.

---

## Part 1 — Critique

### Summary verdict

**Not implementation-ready as-is.** Phase I is close. Phase II is
under-specified in the most consequential place (the sweep recurrence
and the "transfer weight"). Phase III mixes two ideas (an evidence
tensor and a small EM) without giving either an identifying anchor or a
convergence story. Phase IV depends on II and III to land first.

Below: issues ordered by how much they will hurt during implementation.

### C1. The boundary sweep recurrence (§4.3) is numerically unsafe

The plan proposes

```
f[i] = local[i] + w_boundary[i] * f[i-1]
cp[i] = prod_{k=1..i} w_boundary[k]
f[i]  = cp[i] * cumsum(local[i] / cp[i])
```

`cp[i]` is monotonically non-increasing in `[0, 1]` and reaches zero
rapidly when even one `w_boundary[k]` is small. Dividing by `cp` is
catastrophic; the proposed "reset where `w_boundary` is effectively
zero" is unspecified.

The fix is to drop the cumulative-product reformulation entirely and
keep a **plain sequential scan in log space**, accumulating
`alpha_excess` and `beta_excess` directly:

```
alpha_fwd[i] = alpha_local[i] + w_boundary[i] * alpha_fwd[i-1]
beta_fwd[i]  = beta_local[i]  + w_boundary[i] * beta_fwd[i-1]
```

This is an O(R) Python loop per reference; with 5–50k reference slices
in a typical genome each ~10–10000 regions long, it remains cheap (a
few hundred ms total). If the loop later becomes a hot spot, vectorize
per-slice with a Numba kernel or a `scipy.signal.lfilter` IIR call — not
by reintroducing the cumulative-product hack.

Also: the recurrence currently attenuates `f[i-1]` by `w_boundary[i]`,
which is the boundary *between* `i-1` and `i`. This is correct, but the
plan never states how `w_boundary` interacts with the prior on
`alpha_excess`. If `w_boundary = 1`, do we double-count `local[i-1]`'s
prior? Define the recurrence on **excess** quantities only (the
contribution above the prior) so the prior remains added exactly once
when computing the final posterior.

### C2. `gamma_global` is not identifiable from data alone (§3.3)

> `gamma_global = captured density / rho_off`

This requires separating "gDNA on-target" from "gDNA off-target". In
stranded non-capture data, "captured" regions don't exist and
`gamma_global = 1`. In stranded capture data, the only true label for
on-target gDNA is the bait/probe set, which the calibrator does not
receive. In unstranded capture data both ratios are confounded with
expressed regions.

The plan needs one of:

- a partial-supervision anchor (bait BED, or annotation-derived
  expected-on-target mask) — preferred, gated by capture mode;
- a hierarchical model where `gamma_r` is shrunk toward 1 with a
  data-driven precision rather than toward `gamma_global`.

As written, Phase II will produce a `gamma_global` estimate that is
implicitly defined by whichever regions happen to score highly in the
initial Phase I exclusion set. That couples the four-state E-step to
the (deliberately conservative) Phase I anchor exclusion and is fragile.

**Fix:** rename the field to `gamma_shrink_target` and document that in
non-capture mode it is constructed as `1.0` with high precision,
collapsing the captured/offtarget distinction.

### C3. Identifiability of the four-state model is not anchored (§5–§6)

With four states and the proposed evidence terms, the joint posterior
has a global symmetry: a region with low contained density, low
boundary mass, and uninformative strand cannot distinguish
`background` from `expressed_offtarget` without a prior. The plan
mentions `log_prior_states` but never specifies how it is computed.

Required additions:

- `log_prior_states[r, :]` must be derived from `region_arrays.signature`
  (annotation): intergenic regions get a strong prior on
  `background ∪ gdna_only_capture`; exonic regions get a flat prior
  across states unless strand RNA evidence is present; intronic regions
  bias toward `background` with intermediate weight.
- The anchor exclusion set used in Phase I is *not* a one-time training
  set. The same regions should carry a `is_background_seed` flag into
  Phase III, and the E-step should clamp those rows to high
  `p_background` mass for the first 1–2 passes.

Without these anchors, the M-step refits described in §6.2 will drift.

### C4. M-step couplings have no monotonicity guarantee (§6.2)

`rho_off`, `kappa_d`, the boundary mixture, and `gamma_global` are
refit simultaneously from posterior weights. Each refit changes the
likelihood for the others; convergence in five passes is asserted, not
proved. In practice:

- assume monotone ELBO is not available;
- enforce damping: `theta_new = (1-eta) * theta_old + eta * theta_hat`
  with `eta = 0.5` initially, growing to 1.0 only after the state
  probabilities are stable;
- report ELBO-like progress (sum of `logsumexp` over rows) and stop
  when its absolute change is below tolerance, *and* `max |dp_states|`
  is below tolerance.

This makes the iteration debuggable without claiming variational EM.

### C5. Boundary "compatibility weight" mixes two distinct things (§4.2)

`w_boundary` is described both as a transfer weight (how much
information should flow across a boundary) and as a gDNA-fraction
estimate at the boundary. These are not the same.

- A boundary with high gDNA fraction is *informative* about gDNA in
  the adjacent regions; transfer weight should be high.
- A boundary with low gDNA fraction (because RNA crossings dominate)
  carries little gDNA evidence; transfer weight should be **low**.

The unstranded formula

```
w_boundary = unspliced_boundary / (unspliced_boundary + spliced_boundary + prior_mass)
```

confounds expression with gDNA. In unstranded mode, the safer default
is `w_boundary = sigmoid(low_information_score)` driven by total
boundary opportunity `L_side` and `B_boundary`, *not* by spliced ratio.

Define this explicitly in code:

```python
def transfer_weight(B_side, L_side, stranded, kappa_d=None, ...):
    """Return per-boundary scalar in [w_floor, 1.0]."""
```

with a hard floor `w_floor = 0.05` to prevent total information loss
and a soft cap `1 - w_floor`.

### C6. "Likelihoods" of mixture posteriors are category errors (§5.1)

The plan composes `log_likelihoods` by adding `logL_capture`, where
`logL_capture[:, captured]` is described as a "boundary mixture
posterior". Posteriors are not likelihoods. Adding them to a tensor
and then renormalizing double-counts the data through the boundary
prior.

The clean abstraction is **log Bayes factors**:

```
logBF_capture[:, captured]      = log p(B | captured)   - log p(B | not captured)
logBF_gdna_density[:, captured] = log p(D | captured)   - log p(D | not captured)
logBF_expression[:, expressed]  = log p(spliced, strand | expressed) - ...
```

Each row is anchored to its `not_*` baseline; the tensor sum becomes a
log-odds tensor over states, and the `log_prior_states[R, 4]` term is
the *only* prior contribution. This eliminates the boundary-prior
double counting and gives a sane semantic for the E-step.

### C7. `BoundaryTable` derivation needs an exactness contract (§2)

The Python `BoundaryTable` must be derivable from the existing
`RegionCountLedger` arrays without recomputing from scan-level data.
Spell out the bijection:

```
local_i  = region_id - ref_region_offsets[ref_id]
left_b   = ref_boundary_offsets[ref_id] + local_i        # boundary i
right_b  = ref_boundary_offsets[ref_id] + local_i + 1    # boundary i+1
```

and the channel mapping:

```
internal boundary b == boundary between local regions i-1 and i:
  left_region_*_pos[b]  = ledger.boundary_right_*_pos[region_id=i-1]
  left_region_*_neg[b]  = ledger.boundary_right_*_neg[region_id=i-1]
  right_region_*_pos[b] = ledger.boundary_left_*_pos[region_id=i]
  right_region_*_neg[b] = ledger.boundary_left_*_neg[region_id=i]
```

Terminal boundaries: zero, `is_terminal=True`. Consumers must check
`is_terminal` before reading. Add a `validate_boundary_table()` helper
that, in `O(R)`, recomputes a few rows and asserts equality with the
ledger.

### C8. `A_r` units and semantics (§1, §3.4, §7.2)

`A_r` is referenced as both a regional *density ratio* (`A_r =
E[lambda_gdna]/E[rho_off]`) and as a *denominator multiplier*
(`gdna_eff_len = unweighted * bp_weighted_mean(A_r)`).

Pin this in code:

```python
# A_r[r] is unitless; A_r == 1 means "this region's gDNA density
# matches the reference off-target density rho_off".
# Locus-level denominator uses bp-weighted mean over the locus's
# region span:
#     gdna_eff_len_corrected = sum_r (L_r * A_r) / sum_r L_r * gdna_eff_len_raw
```

Add a `dataclass` invariant: `A_r.dtype == float32`, `A_r >= 0`, finite.

### C9. Hard-zero vs. smooth modes are inconsistent (§3.2, §7.1)

§3.2 says "no exact zero mode". §7.1 says `enable_gdna` is computed
from evidence and can still be a boolean. These are in tension: if the
solver returns `upper_gdna > 0` everywhere (because the prior is
smooth), `enable_gdna` is always true and the boolean is vestigial.

**Recommendation:** delete `enable_gdna` from the `PriorTable` and
always run the locus EM with a gDNA component whose prior is set by
`mu_gdna` / `upper_gdna`. If a regime requires shutting off gDNA, it
should be a *configuration* switch (`config.disable_gdna_component`)
not a per-locus data-driven boolean. This removes a class of bugs
where a near-zero but positive `upper_gdna` produces a brittle on/off
flip.

### C10. Strand `kappa_d` cannot be refit per-compartment (§3.1)

`RegionGdnaChannelEstimate` returns contained / boundary-left /
boundary-right posteriors but does not say whether they share
`kappa_d`. With three independent fits and limited per-compartment
seed regions, per-compartment `kappa_d` will be noisy.

Decision the plan must make explicit: **`kappa_d` is one scalar per
library**, fit jointly across compartments. The compartment-aware
posteriors share it. If a future need arises to vary it, add
`kappa_d_contained`, `kappa_d_boundary` as separate fields.

### C11. Phase IV migration risk (§8 Phase IV)

`pipeline.py` already calls `fuse_density_and_strand` and reads
`fused_region_gdna`. Removing that conduit before Phase V benchmarks
have run risks regressing the existing golden tests
(`tests/test_golden_output.py`). Mitigation:

- keep `fused_region_gdna` as a back-compat field on
  `CalibrationResult` that is *derived from* the new
  `RegionCalibration` in Phase IV;
- update goldens only after Phase V demonstrates improvement on the
  benchmark fixtures.

### C12. Minor items

- §5.2 strand contrast: existing `strand_log_likelihood_d_grid` is
  exact for `n_total ≤ 200`. The tensor path must call into it without
  resurrecting the old grid double work — vectorize across regions
  with a precomputed grid.
- §9.4 renames `gdna_density_confidence → gdna_confidence`. Keep both
  names with a deprecation alias for one release.
- §10 region-debug output should be opt-in (`--debug-regions`) — it
  scales as `O(R × float64)` and a human genome has `R ≈ 10⁶`.
- `R6` (mappability) — exposing `mappability_correction(L_eff)` as a
  no-op hook now is cheap and de-risks the future addition.

---

## Part 2 — Implementation Spec

The spec below assumes the critique items above are accepted. It is
written so a reviewer can map each subsection 1:1 to a PR.

### Common ground rules

- All new modules live under `src/rigel/calibration/`.
- All public dataclasses are `@dataclass(frozen=True, slots=True)`.
- All hot arrays use `np.float64` for accumulation and `np.float32`
  for storage on `CalibrationResult` fields.
- All flags are `uint16` bitfields; constants live next to the module
  that produces them.
- No new C++ code in Phases I–IV. Phase V may propose a C++
  accumulator if profiling demands.

### Module map (new + edited)

```
src/rigel/calibration/
├── boundaries.py            # NEW  Phase I  BoundaryTable
├── background_model.py      # NEW  Phase I  rho_off + seed selection
├── boundary_model.py        # NEW  Phase II Gamma-Poisson boundary->contained
├── boundary_sweep.py        # NEW  Phase II Forward/backward scan
├── transfer_weight.py       # NEW  Phase II w_boundary computation
├── latent_states.py         # NEW  Phase III [R,4] log-BF tensor + EM
├── _result.py               # EDIT Phase III add RegionCalibration field
├── strand_deconv.py         # EDIT Phase I  add compartment-aware API
├── _orchestrator.py         # EDIT Phase IV rewire pipeline
├── exposure.py              # EDIT Phase IV consume RegionCalibration
├── prior.py                 # EDIT Phase IV consume RegionCalibration
└── integration.py           # KEEP Phase IV deprecate, derive from new
```

### Phase I — Compartment strand + background + BoundaryTable

#### Phase I.1 `boundaries.py`

```python
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ._arrays import RegionArrays
from .region_count_ledger import RegionCountLedger


@dataclass(frozen=True, slots=True)
class BoundaryTable:
    boundary_pos: np.ndarray            # int64[B]   genomic coord of boundary
    ref_id: np.ndarray                  # int32[B]
    ref_region_offsets: np.ndarray      # int32[n_refs+1]
    ref_boundary_offsets: np.ndarray    # int32[n_refs+1]
    is_terminal: np.ndarray             # bool[B]

    # 8 mass channels (4 left-side, 4 right-side):
    left_region_unspliced_pos: np.ndarray   # float32[B]
    left_region_unspliced_neg: np.ndarray
    left_region_spliced_pos:   np.ndarray
    left_region_spliced_neg:   np.ndarray
    right_region_unspliced_pos: np.ndarray
    right_region_unspliced_neg: np.ndarray
    right_region_spliced_pos:   np.ndarray
    right_region_spliced_neg:   np.ndarray

    left_region_boundary_leff:  np.ndarray  # float64[B]
    right_region_boundary_leff: np.ndarray  # float64[B]

    @property
    def n_boundaries(self) -> int: ...
    def left_region_index(self) -> np.ndarray: ...   # int64[B], -1 at left terminal
    def right_region_index(self) -> np.ndarray: ...  # int64[B], -1 at right terminal


def build_boundary_table(
    region_arrays: RegionArrays,
    ledger: RegionCountLedger,
    side_leff: np.ndarray,             # float64[R] (= DensityObservation.boundary_left_leff)
) -> BoundaryTable:
    """Derive B = R + n_refs boundary slots from sorted regions."""


def validate_boundary_table(
    bt: BoundaryTable,
    region_arrays: RegionArrays,
    ledger: RegionCountLedger,
) -> None:
    """O(R) recompute-and-assert helper used only in tests."""
```

Concrete steps inside `build_boundary_table`:

1. `n_refs = region_arrays.n_refs`; per-ref boundary count =
   `(ref_region_offsets[r+1] - ref_region_offsets[r]) + 1`;
   `ref_boundary_offsets` is `cumsum`.
2. Allocate `B = R + n_refs` arrays, all zero-initialised.
3. For each `r` in `0..R-1`, with `i_local = r - ref_region_offsets[ref]`,
   compute `left_b = ref_boundary_offsets[ref] + i_local` and
   `right_b = left_b + 1`. Vectorise across all `r`.
4. Scatter ledger right-side channels into `left_region_*` slots; scatter
   ledger left-side channels into `right_region_*` slots.
5. Scatter `side_leff[r]` into both `left_region_boundary_leff[right_b]`
   and `right_region_boundary_leff[left_b]`.
6. Mark `is_terminal[ref_boundary_offsets[r]]` and
   `is_terminal[ref_boundary_offsets[r+1]-1]` True for every ref.
7. `boundary_pos`: for internal slot `b` between regions `i-1, i`,
   set to `region_arrays.end[i-1]` (equivalently `region_arrays.start[i]`,
   which equals it by adjacency invariant). For terminals use the ref's
   first `start` (left) and last `end` (right) — see assertion below.

Invariants (must be `assert`-checked):

- `ref_boundary_offsets[-1] == B`
- For every internal `b`, `boundary_pos[b] == region_arrays.start[i] == region_arrays.end[i-1]`.

Tests (`tests/test_boundaries.py`):

- empty index → `B == n_refs`
- 3-region single-ref fixture → exact reconstruction of left/right slots
- mass conservation: sum over `left_region_unspliced_pos` equals
  `ledger.boundary_right_unspliced_pos.sum()`.

#### Phase I.2 `strand_deconv.py` — compartment-aware extension

Add a new function (do not delete the existing aggregate path; it
remains used by `kappa_d` estimation):

```python
def build_compartment_strand_counts(
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    *, p_r1_sense: float,
) -> "CompartmentStrandCounts":
    ...


@dataclass(frozen=True, slots=True)
class CompartmentStrandCounts:
    # contained:
    k_sense_c: np.ndarray
    k_anti_c:  np.ndarray
    n_total_c: np.ndarray
    # boundary left:
    k_sense_bl: np.ndarray
    k_anti_bl:  np.ndarray
    n_total_bl: np.ndarray
    # boundary right:
    k_sense_br: np.ndarray
    k_anti_br:  np.ndarray
    n_total_br: np.ndarray
    eligible:  np.ndarray
    p_r1_sense: float
```

Then add:

```python
def deconvolve_compartments_by_strand(
    counts: CompartmentStrandCounts,
    *, kappa_d: float,
    rna_lower_confidence: float,
) -> "RegionGdnaChannelEstimate":
    """Three-compartment posterior reusing strand_log_likelihood_d_grid
    once per compartment. Shares kappa_d (see critique C10)."""
```

Where `RegionGdnaChannelEstimate` has the nine `*_mean / *_upper /
*_rna_lower` arrays plus shared `kappa_d`, `p_r1_sense`, `confidence`,
and `flags[R]` (uint16).

Tests (`tests/test_compartment_strand_deconv.py`):

- aggregate sum of compartment `k_sense + k_anti` equals the v4
  aggregate `n_total`.
- A region with only contained mass → boundary posteriors equal the
  prior; contained posterior matches the v4 aggregate path bit-for-bit.

#### Phase I.3 `background_model.py`

```python
@dataclass(frozen=True, slots=True)
class BackgroundModel:
    rho_off_alpha: float
    rho_off_beta: float
    rho_off_mean: float
    n_seed_regions: int
    seed_mask: np.ndarray         # bool[R]   -- carried forward!
    flags: np.ndarray             # uint16[R]
    fit_status: str               # "ok" | "fallback_floor" | "fallback_uniform"


def fit_background_model(
    region_arrays: RegionArrays,
    observation: DensityObservation,
    strand_channels: RegionGdnaChannelEstimate | None,   # None = unstranded
    *,
    top_t_fraction: float = 0.02,
    min_seed_regions: int = 200,
    alpha_floor: float = 1.0,
    beta_floor: float = 1.0,
) -> BackgroundModel:
    ...
```

Algorithm:

1. Build initial candidate mask: intergenic OR intron-only (use
   `fractional_evidence.is_intergenic`, `is_intron_only`), with
   `observation.contained_leff > min_eff_length` and
   `observation.spliced_count == 0`.
2. Compute initial naive density `d_r = contained / max(contained_leff, eps)`.
3. Trim top `top_t_fraction` of candidates by `d_r` (union exclusion).
4. If `strand_channels is not None`: also exclude any region with
   `strand_channels.contained_rna_lower > 0`.
5. Fit Gamma(α, β) via MoM on the surviving set; floor `α ≥ alpha_floor`,
   `β ≥ beta_floor`.
6. If `n_eligible < min_seed_regions`, mark `fit_status =
   "fallback_floor"` and use prior-only.

The `seed_mask` is the **post-trim, post-RNA-exclusion** mask; this
mask is what Phase III consumes to anchor `p_background` in early
passes (critique C3).

Tests (`tests/test_background_model.py`):

- Synthetic uniform-density fixture: `rho_off_mean ≈ true λ ±10%`.
- Strand-RNA-positive intron is excluded from seed when stranded.
- Empty candidate set → `fit_status == "fallback_floor"`.

#### Phase I integration

`_orchestrator.calibrate` continues to return a `CalibrationResult`,
but now additionally:

- builds `BoundaryTable`;
- builds `CompartmentStrandCounts` and
  `deconvolve_compartments_by_strand` (in stranded mode);
- builds `BackgroundModel`;
- stores these as new optional fields on `CalibrationResult` so Phase
  II/III can read them. `fused_region_gdna` / `region_exposure` paths
  remain.

Gate: existing `pytest tests/` passes unchanged (Phase I is additive).

### Phase II — Boundary imputation + sweep

#### Phase II.1 `transfer_weight.py`

```python
W_FLOOR: float = 0.05
W_CAP:   float = 1.0 - W_FLOOR

def compute_transfer_weights(
    bt: BoundaryTable,
    strand_channels: RegionGdnaChannelEstimate | None,
    background: BackgroundModel,
) -> np.ndarray:   # float64[B]   in [W_FLOOR, W_CAP]
    """Per-boundary gDNA-information transfer weight.

    Stranded:
        gdna_at_boundary = strand_channels.boundary_left_mean[right_region]
                          + strand_channels.boundary_right_mean[left_region]
        opportunity      = bt.left_region_boundary_leff[b]
                          + bt.right_region_boundary_leff[b]
        info             = (gdna_at_boundary + alpha_floor) / (opportunity + beta_floor)
        w[b]             = sigmoid( log(info / rho_off_mean) * gain )
    Unstranded:
        total_b = sum of all 8 channels at b
        L_b     = left+right leff
        w[b]    = 1.0 - exp(-total_b / max(L_b * rho_off_mean, eps))
    Terminal boundaries:
        w[b]    = 0.0  (no transfer)
    """
```

Critical: gain is a fixed internal constant (`gain = 1.0`) — do not
expose as a knob until Phase V benchmarks ask for it.

#### Phase II.2 `boundary_model.py`

```python
@dataclass(frozen=True, slots=True)
class BoundaryLocalPosterior:
    alpha_excess: np.ndarray       # float64[R]   from left+right boundaries
    beta_excess:  np.ndarray       # float64[R]
    mu_local:     np.ndarray       # float32[R]   predictive mean of contained gDNA
    upper_local:  np.ndarray       # float32[R]   predictive upper at confidence c
    flags:        np.ndarray       # uint16[R]


def fit_boundary_local(
    region_arrays: RegionArrays,
    observation: DensityObservation,
    bt: BoundaryTable,
    strand_channels: RegionGdnaChannelEstimate | None,
    background: BackgroundModel,
    *, confidence: float,
) -> BoundaryLocalPosterior:
    """Predict contained gDNA from this region's own two boundaries.

    For each region r with boundaries (left_b, right_b):
        B_side, L_side   = stranded-gDNA counts (or unspliced if unstranded)
        alpha_post = alpha_floor + B_left + B_right
        beta_post  = beta_floor  + L_left + L_right
        Predictive D | lambda ~ Poisson(lambda * contained_leff_r)
        Marginal D ~ NegBin(r = alpha_post,
                            p = beta_post / (beta_post + contained_leff_r))
        mu_local[r], upper_local[r] = NB mean / quantile
        alpha_excess[r] = alpha_post - alpha_floor
        beta_excess[r]  = beta_post  - beta_floor
    """
```

Use `scipy.stats.nbinom`. Vectorize across all R.

#### Phase II.3 `boundary_sweep.py`

```python
@dataclass(frozen=True, slots=True)
class BoundarySweepPosterior:
    alpha_excess_fwd: np.ndarray   # float64[R]
    beta_excess_fwd:  np.ndarray
    alpha_excess_rev: np.ndarray
    beta_excess_rev:  np.ndarray
    alpha_excess_total: np.ndarray
    beta_excess_total:  np.ndarray
    mu_sweep: np.ndarray           # float32[R]   NB predictive mean using combined posterior
    upper_sweep: np.ndarray        # float32[R]
    flags: np.ndarray              # uint16[R]


def run_boundary_sweep(
    region_arrays: RegionArrays,
    bt: BoundaryTable,
    local: BoundaryLocalPosterior,
    w: np.ndarray,                  # transfer weights, float64[B]
    background: BackgroundModel,
    observation: DensityObservation,
    *, confidence: float,
) -> BoundarySweepPosterior:
    ...
```

Implementation (critique C1 — replaces the cumulative-product scheme):

```python
def _scan_ref_slice(alpha_local, beta_local, w_internal):
    """Forward IIR: x[i] = local[i] + w[i] * x[i-1].

    Args:
        alpha_local: float64[N]   per-region local excess in this ref slice.
        beta_local:  float64[N]
        w_internal:  float64[N]   w_internal[0] == 0 (no boundary before first region).
    """
    n = alpha_local.size
    out_a = np.empty(n, dtype=np.float64)
    out_b = np.empty(n, dtype=np.float64)
    a = 0.0
    b = 0.0
    for i in range(n):
        a = alpha_local[i] + w_internal[i] * a
        b = beta_local[i]  + w_internal[i] * b
        out_a[i] = a
        out_b[i] = b
    return out_a, out_b
```

This per-slice scan can be Numba-jitted later. Reverse pass is the
same on reversed slices. Combine:

```
alpha_total = alpha_local + alpha_fwd_prev + alpha_rev_next
beta_total  = beta_local  + beta_fwd_prev  + beta_rev_next
```

where `*_prev/next` exclude the local contribution to avoid double
counting. Implement by passing `alpha_local_prev_only[i] = w[i] *
fwd[i-1]` and similar for reverse.

`mu_sweep`, `upper_sweep` are NB predictive moments from
`alpha = alpha_floor + alpha_total`, `beta = beta_floor + beta_total`,
`leff = observation.contained_leff`.

Tests (`tests/test_boundary_sweep.py`):

- Single isolated region with `w == 0` everywhere → sweep equals
  local.
- Three contiguous regions, middle has zero boundaries, outer two
  have strong boundary mass → middle's `mu_sweep > mu_local`.
- T1/T2 giant-exon fixture: outer-boundary gDNA evidence propagates
  inward; central region's `upper_sweep` is strictly greater than
  `upper_local`.
- Stress: 10⁵ regions in one ref runs in <0.5 s without Numba.

### Phase III — Four-state log-BF tensor + calibration EM

#### Phase III.1 `latent_states.py` — constants and shapes

```python
STATE_BACKGROUND          = 0
STATE_GDNA_ONLY_CAPTURE   = 1
STATE_EXPRESSED_CAPTURE   = 2
STATE_EXPRESSED_OFFTARGET = 3
N_STATES = 4

STATE_IS_EXPRESSED = np.array([0, 0, 1, 1], dtype=np.int8)   # column index
STATE_IS_CAPTURED  = np.array([0, 1, 1, 0], dtype=np.int8)


@dataclass(frozen=True, slots=True)
class RegionCalibration:
    p_states:   np.ndarray   # float32[R, 4]
    mu_gdna:    np.ndarray   # float32[R]
    upper_gdna: np.ndarray   # float32[R]
    rna_lower:  np.ndarray   # float32[R]
    A_r:        np.ndarray   # float32[R]
    gamma_r:    np.ndarray   # float32[R]
    rho_off:    float
    kappa_d:    float | None
    n_passes:   int
    converged:  bool
    flags:      np.ndarray   # uint16[R]

    # convenience properties: p_background, p_captured, p_expressed, ...
```

#### Phase III.2 Likelihood terms (vectorised, log-BF anchored)

```python
def log_bf_expression(
    observation: DensityObservation,
    strand_channels: RegionGdnaChannelEstimate | None,
    background: BackgroundModel,
    sweep: BoundarySweepPosterior,
) -> np.ndarray:   # float64[R, 2]   columns: not_expressed, expressed
    """Spliced mass + RNA lower bound + excess contained vs gDNA upper."""


def log_bf_capture(
    boundary_local: BoundaryLocalPosterior,
    boundary_sweep: BoundarySweepPosterior,
    background: BackgroundModel,
    gamma_shrink_target: float,
) -> np.ndarray:   # float64[R, 2]   columns: not_captured, captured


def log_bf_gdna_density(
    observation: DensityObservation,
    background: BackgroundModel,
    strand_channels: RegionGdnaChannelEstimate | None,
    captured_density: float,    # gamma_shrink_target * rho_off
) -> np.ndarray:   # float64[R, 2]


def log_bf_strand(
    strand_channels: RegionGdnaChannelEstimate | None,
    kappa_d: float | None,
    n_total: np.ndarray,
) -> np.ndarray:   # float64[R, 4]   state-specific
```

Each function is fully vectorized; **no Python loops over regions**.
Sentinels for unstranded (return zeros).

#### Phase III.3 Tensor assembly

```python
def assemble_state_posterior(
    log_prior: np.ndarray,                # float64[R, 4]
    logbf_expression: np.ndarray,         # float64[R, 2]
    logbf_capture: np.ndarray,            # float64[R, 2]
    logbf_gdna_density: np.ndarray,       # float64[R, 2]
    logbf_strand: np.ndarray,             # float64[R, 4]
    seed_mask: np.ndarray,                # bool[R] — clamps to background
    pass_index: int,
) -> tuple[np.ndarray, np.ndarray]:       # (p_states[R,4], log_evidence[R])
    expr_col = STATE_IS_EXPRESSED          # broadcasts via fancy index
    cap_col  = STATE_IS_CAPTURED
    L = log_prior.copy()
    L += logbf_expression[:, expr_col]
    L += logbf_capture[:, cap_col]
    L += logbf_gdna_density[:, cap_col]
    L += logbf_strand
    # Anchor seed regions for the first 2 passes:
    if pass_index < 2:
        L[seed_mask, STATE_BACKGROUND] += 5.0   # ~e^5 odds-on
    log_norm = logsumexp(L, axis=1)
    p = np.exp(L - log_norm[:, None]).astype(np.float32)
    return p, log_norm
```

Note `log_prior_states[r, :]` comes from `region_arrays.signature`:

```python
def initial_state_log_prior(region_arrays, has_capture: bool) -> np.ndarray:
    """Annotation-derived weak prior. Intergenic → background heavy;
    exonic → flat; intronic → background-heavy when no_capture."""
```

#### Phase III.4 The calibration EM

```python
def run_region_calibration(
    region_arrays: RegionArrays,
    observation: DensityObservation,
    bt: BoundaryTable,
    background: BackgroundModel,
    strand_channels: RegionGdnaChannelEstimate | None,
    *,
    max_passes: int = 5,
    p_tol: float = 1e-3,
    log_evidence_tol: float = 1e-4,
    damping: float = 0.5,
) -> RegionCalibration:
    """Drives Phase II/III E/M loop.

    M-step refits (with damping per critique C4):
        rho_off       <- weighted by p_background
        kappa_d       <- weighted by (1 - p_expressed) on stranded data
        gamma_target  <- weighted by p_captured / p_offtarget ratio
        boundary mix  <- weighted by p_captured contributions
    """
```

Loop structure (deliberately explicit so the implementer matches
critique C4):

```
for k in range(max_passes):
    boundary_local = fit_boundary_local(..., background=background_k, ...)
    w              = compute_transfer_weights(..., background=background_k, ...)
    sweep          = run_boundary_sweep(..., w=w, ...)
    L_expr  = log_bf_expression(...)
    L_cap   = log_bf_capture(...)
    L_dens  = log_bf_gdna_density(...)
    L_str   = log_bf_strand(...)
    p_new, log_ev_new = assemble_state_posterior(...)
    delta_p  = max(|p_new - p_old|)
    delta_le = sum(log_ev_new - log_ev_old)
    if delta_p < p_tol and |delta_le| < log_evidence_tol:
        break
    background_k_plus_1 = damped(refit_background(p_new), background_k, damping)
    kappa_d_k_plus_1    = damped(refit_kappa(p_new), kappa_d_k, damping)
    gamma_target_k_plus_1 = damped(refit_gamma(p_new), gamma_target_k, damping)
```

After convergence, derive:

```
mu_gdna    = p_captured * mu_sweep_captured + p_offtarget * mu_sweep_offtarget
upper_gdna = NB upper from combined alpha/beta posterior weighted by states
rna_lower  = strand_channels.contained_rna_lower (when stranded) else 0
A_r        = mu_gdna / max(rho_off * contained_leff, eps)
gamma_r    = mu_gdna_captured / max(rho_off * contained_leff, eps)
```

`flags` collects: ambiguous (max p < 0.5), strand uninformative,
boundary sparse, sweep damped, exact vs approx strand path, seed anchor.

Tests (`tests/test_latent_states.py`,
`tests/test_calibration_iteration.py`):

- Four-mode synthetic mixture (background-only, captured-only,
  expressed-only, mixed): each archetype's posterior `p_*` is the
  largest column.
- Damped M-step: passing pathological pre-fit `rho_off = 100x` true
  converges within 5 passes to the right value.
- Unstranded fixture: `logbf_strand == 0`, but four states still
  identifiable for the corner cases (background vs expressed_offtarget
  remains ambiguous as designed; flag is set).

### Phase IV — Downstream wiring

#### Phase IV.1 `_orchestrator.py`

Drop-in replacement at the bottom of `calibrate(...)`:

```python
boundary_table = build_boundary_table(region_arrays, ledger,
                                      observation.boundary_left_leff)
strand_compartments = (
    build_compartment_strand_counts(region_arrays, payload_arrays,
                                    p_r1_sense=strand_summary.p_r1_sense)
    if strand_summary.is_informative else None
)
region_channels = (
    deconvolve_compartments_by_strand(strand_compartments,
                                      kappa_d=kappa_d.kappa,
                                      rna_lower_confidence=rna_lower_confidence)
    if strand_compartments is not None else None
)
background = fit_background_model(region_arrays, observation,
                                  region_channels, ...)
region_cal = run_region_calibration(region_arrays, observation,
                                    boundary_table, background,
                                    region_channels, ...)
# Back-compat: derive the old types from RegionCalibration:
region_exposure = RegionExposure.from_region_calibration(region_cal)
fused_region_gdna = FusedRegionGdnaEvidence.from_region_calibration(region_cal)
```

`build_calibration_result(...)` gains a new `region_calibration` field
on `CalibrationResult`.

#### Phase IV.2 `exposure.py`

Add a classmethod that consumes `RegionCalibration` and writes the
existing `RegionExposure` schema unchanged:

```python
@classmethod
def from_region_calibration(cls, rc: RegionCalibration) -> "RegionExposure":
    return cls(
        mode="density",
        A_r=np.asarray(rc.A_r, dtype=np.float32),
        rho_r=(rc.A_r * rc.rho_off).astype(np.float32),
        rho_ref=float(rc.rho_off),
        reference_quantile=0.0,
        eligible=(rc.flags & FLAG_INELIGIBLE == 0),
        flags=rc.flags.astype(np.uint8),
    )
```

#### Phase IV.3 `prior.py`

`assemble_priors` already consumes per-region `mu / upper` via
`FusedRegionGdnaEvidence`. Add a parallel constructor:

```python
def assemble_priors_from_region_calibration(
    multi_loci: list[MultiLocus],
    rc: RegionCalibration,
    region_arrays: RegionArrays,
    ...,
) -> PriorTable:
    """Same allocation logic, reading rc.mu_gdna / rc.upper_gdna /
    rc.A_r instead of fused_region_gdna."""
```

Per critique C9, drop `enable_gdna` driving by `upper_count > eps`;
instead set `enable_gdna[:] = not config.disable_gdna_component`.

#### Phase IV.4 `pipeline.py`

Switch the call from `assemble_priors(...)` to
`assemble_priors_from_region_calibration(...)` behind a feature flag
`PipelineConfig.use_v5_calibration` (default `False` until Phase V
benchmarks pass). When flag is `True`, the old `fused_region_gdna`
path is unused but kept compiled.

Tests:

- All `tests/test_golden_output.py` cases pass with the flag `False`
  (unchanged behavior).
- A new `tests/test_golden_v5_output.py` records the v5-flag-on
  outputs and asserts deterministic reproduction.
- `pytest -k "calibration"` continues to pass.

### Phase V — Benchmark + refit

Out of scope for code generation here. Spec:

1. Extend the simulator to emit post-capture truth (existing
   `scripts/sim/` helpers already track per-fragment origin).
2. Run the four-regime grid (`scripts/benchmark/configs/`) with both
   `use_v5_calibration=False` and `=True`.
3. Tabulate mRNA / nRNA / gDNA relative error; require no regression
   vs v4 in any regime and ≥10% improvement in stranded-capture
   high-gDNA fixtures.
4. Only when the above passes, flip `use_v5_calibration` default to
   `True` and update goldens.

---

## Acceptance checklist (per phase)

- [ ] **Phase I** — `BoundaryTable`, `BackgroundModel`,
  `CompartmentStrandCounts` land additively; existing tests pass.
- [ ] **Phase II** — sweep produces sensible values on T1/T2 fixture;
  10⁵-region per-ref scan runs in <0.5 s.
- [ ] **Phase III** — synthetic four-mode fixture identifies the
  correct dominant state per archetype; EM converges in ≤ 5 passes.
- [ ] **Phase IV** — `use_v5_calibration=False` keeps every existing
  golden test byte-for-byte; `use_v5_calibration=True` produces a new
  golden corpus that round-trips deterministically.
- [ ] **Phase V** — benchmarks show no regression and target
  improvement on the agreed fixtures.

