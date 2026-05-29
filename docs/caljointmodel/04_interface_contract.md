# 04 — Interface Contract

This document defines the public API of the new calibration module
and — critically — the **exact algebraic translation** from
calibration output to the per-transcript Dirichlet pseudocounts that
the locus EM consumes. The pseudocount formulas appear here, before
any Phase 4 code is written.

## 1. Module surface

```
src/rigel/calibration/calibrate.py    # the new module (single file, ~300 LoC)
src/rigel/calibration/__init__.py     # re-exports only the symbols below
```

Public surface — three names:

```python
from rigel.calibration import (
    CalibrationConfig,
    CalibrationResult,
    calibrate,
)
```

## 2. Substrate (native)

> **Status (2026-05-29):** the calibrator's substrate is the
> `AccumulatorPayload` produced by accumulator v2 — see
> [`../accumulator/00_design.md`](../accumulator/00_design.md).
> Accumulator v2 lands **before** calibration v6 Phase 3; until then,
> the field names below describe the *target* schema. The
> "extension" framing in the original draft is obsolete: regions and
> boundaries are first-class native objects, not a side-table on top
> of `RegionCountLedger`.

The native BAM scanner emits per-region and per-boundary counters
directly. The calibrator consumes them as the substrate sufficient
statistics defined in [`01_generative_model.md`](01_generative_model.md) §1.

### 2.1 New native struct

```cpp
// src/rigel/native/accumulator.h  (additions)
struct CalibrationAggregates {
  // Per region (parallel to RegionCountLedger row index)
  std::vector<uint32_t> region_n_unspliced;            // (R,)
  std::vector<uint32_t> region_n_spliced;              // (R,)
  std::vector<uint32_t> region_k_plus_unspliced;       // (R,)

  // Per boundary (left/right of each region)
  std::vector<uint32_t> left_n_unspliced;              // (R,)
  std::vector<uint32_t> left_n_spliced;                // (R,)
  std::vector<uint32_t> left_k_plus_unspliced;         // (R,)

  std::vector<uint32_t> right_n_unspliced;             // (R,)
  std::vector<uint32_t> right_n_spliced;               // (R,)
  std::vector<uint32_t> right_k_plus_unspliced;        // (R,)
};
```

No FL histograms — see
[`00_overview.md`](00_overview.md) §1 and the accumulator v2 design
([`../accumulator/audit_phase1.md`](../accumulator/audit_phase1.md)
decision #6).

### 2.2 Memory budget

At $R = 200{,}000$:
- 9 scalar counters × 4 bytes × $R$ = 7.2 MB

Total: **~7 MB** added to the native substrate. An order of magnitude
smaller than previous (FL-bearing) designs.

### 2.3 Native API addition

```cpp
// nanobind binding in _bam_impl
nb::class_<CalibrationAggregates>(m, "CalibrationAggregates")
  .def_ro("region_n_unspliced",        &CalibrationAggregates::region_n_unspliced)
  .def_ro("region_n_spliced",          &CalibrationAggregates::region_n_spliced)
  .def_ro("region_k_plus_unspliced",   &CalibrationAggregates::region_k_plus_unspliced)
  // ... boundaries analogously ...
  ;
```

Python adapter in `scan_payload.py` (KEEP):

```python
@dataclass(frozen=True, slots=True)
class CalibrationSubstrate:
    # Region arrays  (each shape (R,))
    region_n_unspliced:        np.ndarray  # int64
    region_n_spliced:          np.ndarray
    region_k_plus_unspliced:   np.ndarray

    # Left boundary arrays
    left_n_unspliced:          np.ndarray
    left_n_spliced:            np.ndarray
    left_k_plus_unspliced:     np.ndarray

    # Right boundary arrays
    right_n_unspliced:         np.ndarray
    right_n_spliced:           np.ndarray
    right_k_plus_unspliced:    np.ndarray

    # Per-region scalar metadata (re-exposed from RegionArrays for convenience)
    L_eff:                     np.ndarray  # float64
    strand:                    np.ndarray  # int8
    signature:                 np.ndarray  # int8
```

The existing `RegionCountLedger` and 12-channel `region_counts` matrix
**continue to exist**; they are still used elsewhere in the pipeline.
The new substrate is purely additive.

## 3. Configuration

```python
@dataclass(frozen=True, slots=True)
class CalibrationConfig:
    """Five fields. None are decision thresholds."""

    max_outer_iterations: int = 25
    """Hard cap on outer EM iterations. Typical convergence in 3-10."""

    mass_rel_tol: float = 1e-4
    """Relative max-region mass-change termination criterion."""

    phi_floor: float = 1e-6
    """Floor on phi to prevent 1/phi blowup in exposure posterior."""

    boundary_split_factor: float = 0.5
    """Fraction of each boundary's gDNA mass attributed to the inside
    region. Default 0.5 (symmetric split). Only knob that materially
    shifts the exposure estimate."""
```

## 4. Inputs to `calibrate()`

```python
def calibrate(
    substrate: CalibrationSubstrate,
    strand_model: StrandModel,
    config: CalibrationConfig = CalibrationConfig(),
) -> CalibrationResult: ...
```

### 4.1 Substrate invariants enforced at entry

1. All region/boundary arrays have consistent length $R$.
2. `strand_model.p_r1_sense.shape == (R,)`.

Violations raise `CalibrationSubstrateError` immediately.

## 5. Output

```python
@dataclass(frozen=True, slots=True)
class CalibrationResult:
    """Per-region deconvoluted mass + exposure posteriors + library
    hyperparameters. Each output traces to one of G2–G4."""

    # --- G2: contained-fragment mass (length R) ---
    mass_g_contained: np.ndarray   # float64[R]
    mass_d_contained: np.ndarray   # float64[R]

    # --- G3: boundary-flux mass (length R) ---
    mass_g_left:  np.ndarray       # float64[R]
    mass_d_left:  np.ndarray       # float64[R]
    mass_g_right: np.ndarray       # float64[R]
    mass_d_right: np.ndarray       # float64[R]

    # --- G4: per-region exposure posterior (length R) ---
    omega:         np.ndarray      # float64[R]  = E[omega | data]   (Gamma posterior mean)
    log_omega_var: np.ndarray      # float64[R]  = Var(log omega | data)  (delta-method)

    # --- library hyperparameters ---
    rho_0:     float               # > 0
    phi:       float               # > 0  (NB dispersion = Gamma exposure prior shape)
    rho_d_bb:  float               # in (0, 1), gDNA strand BB dispersion (kappa_d = 0.5 fixed, not in schema)
    rho_r_bb:  float               # in (0, 1), RNA strand BB dispersion (kappa_rna pre-computed by StrandModel, not in schema)
    eps_s:     float               # in (0, 1), gDNA splice-artifact rate

    # --- convergence diagnostics ---
    n_iterations:        int
    converged:           bool
    mass_change_history: np.ndarray  # float64[n_iterations]

    # --- provenance ---
    n_regions:   int
    config:      CalibrationConfig
```

### 5.1 Invariants (checked in `__post_init__`)

- All `mass_*` arrays are length $R$ and non-negative.
- For each region: `mass_g_contained[r] + mass_d_contained[r]` equals the contained fragment count (unspliced + spliced) to within $10^{-9}$. Same for each boundary.
- `omega[r] > 0`, `log_omega_var[r] > 0`.
- `0 < rho_d_bb, rho_r_bb, eps_s < 1`; `rho_0 > 0`; `phi > 0`.
- `mass_change_history` is non-increasing (within $10^{-9}$ scale tolerance).

## 6. Downstream consumer contract: locus prior pseudocounts

There is one downstream consumer: `compute_locus_priors_from_partitions` in `src/rigel/locus.py`.

### 6.1 Signature

```python
def compute_locus_priors_from_partitions(
    locus_partitions: Sequence[LocusPartition],
    region_arrays: RegionArrays,
    calibration: CalibrationResult,
    *,
    kappa: float,                  # locus-level prior concentration
) -> LocusPriors: ...
```

### 6.2 Per-transcript pseudocount formulas

For each transcript $t$ in a locus, the EM has two state classes per
transcript ("true RNA" and "gDNA contamination"). The Dirichlet prior
receives two pseudocounts per transcript:

Define the **region-of-transcript fractional allocator** $\phi_{t,r}$
(the existing `LocusPartition` field) and the **inverse-variance
weight** $w_r = 1 / (1 + \hat{\sigma}_{\log\omega_r}^2)$ — the latter
down-weights high-uncertainty regions naturally.

**RNA pseudocount:**

$$
\boxed{\;\alpha_t^{(d)} = \kappa \cdot \sum_{r \in r(t)} \phi_{t,r} \cdot w_r \cdot \bigl[M_r^{(d, \text{cont})} + \tfrac{1}{2}(M_r^{(d,L)} + M_r^{(d,R)})\bigr]\;}
$$

**gDNA pseudocount:**

$$
\boxed{\;\alpha_t^{(g)} = \kappa \cdot \sum_{r \in r(t)} \phi_{t,r} \cdot w_r \cdot \hat{\omega}_r \cdot \rho_0 \cdot L_{t,r}^{\text{eff}}\;}
$$

where:
- $r(t)$ = regions overlapping transcript $t$,
- $\phi_{t,r}$ = fraction of region $r$ inside transcript $t$'s exonic footprint (from `LocusPartition`),
- $w_r = 1 / (1 + \hat{\sigma}_{\log\omega_r}^2)$ = inverse-variance weight,
- $L_{t,r}^{\text{eff}}$ = effective length contribution; v1: $\phi_{t,r} \cdot L_r^{\text{eff}}$,
- $\kappa$ = user knob, inherited from `EMConfig`.

### 6.3 Boundary mass conservation

The half-split in $\alpha_t^{(d)}$ matches the half-split in
$M_r^{(g, \text{tot})}$ used for exposure inference (see
[`03_inference.md`](03_inference.md) §4), preserving total
fragment-mass conservation. For a locus where transcripts fully tile
the regions:

$$
\sum_t \alpha_t^{(d)} + \alpha_t^{(g)} \;=\; \kappa \cdot \sum_r w_r \cdot \bigl[\hat{\omega}_r \rho_0 L_r^{\text{eff}} + n_r^{\text{u}} + n_r^{\text{s}} - M_r^{(g, \text{cont})}\bigr]
$$

This invariant is the basis of the locus-prior consumer unit test.

### 6.4 Degenerate-case behavior

| Case | Behavior |
|---|---|
| Region with 0 fragments | $M = 0$, $\hat{\omega}_r = 1$ (prior mean), $\hat{\sigma}_{\log\omega_r}^2 = \phi$, $w_r = 1/(1+\phi)$. Contributes only $\kappa w_r \rho_0 L_{t,r}^{\text{eff}}$ to $\alpha_t^{(g)}$. |
| Paralog multimap-only region | $M_r^{(d, \text{cont})} \approx 0$ via three-leg rescue (count NB + BB strand + zero spliced; see [`02_failure_audit.md`](02_failure_audit.md) §2.5). Symmetric pseudocounts across paralog transcripts; pathology cannot occur once benchmarks verify the rescue. |
| Transcript not overlapping any region | $\alpha_t^{(d)} = \alpha_t^{(g)} = 0$; handled by EM's existing Dirichlet floor. |
| Region with extreme $\hat{\omega}_r \gg 1$ | gDNA pseudocount scales linearly; intended behavior. |

## 7. Pipeline wiring (post-rewrite sketch)

```python
# src/rigel/pipeline.py::quant_from_buffer (post-rewrite, ~12 lines)

def quant_from_buffer(buffer, index, *, config):
    # ... existing BAM scan ...

    substrate = build_calibration_substrate(buffer, index)            # new helper
    region_arrays = RegionArrays.from_index(index)

    # G2 + G3 + G4: the new calibrator
    calibration = calibrate(
        substrate=substrate,
        strand_model=strand_model,
        config=config.calibration,
    )

    locus_priors = compute_locus_priors_from_partitions(
        locus_partitions=partitions,
        region_arrays=region_arrays,
        calibration=calibration,
        kappa=config.em.locus_prior_kappa,
    )

    # ... EM ...
```

## 8. `PipelineConfig` change

```python
@dataclass(frozen=True, slots=True)
class PipelineConfig:
    em: EMConfig = field(default_factory=EMConfig)
    scan: BamScanConfig = field(default_factory=BamScanConfig)
    calibration: CalibrationConfig = field(default_factory=CalibrationConfig)
    # All legacy calibration sub-configs deleted.
```

## 9. Errors

```python
class CalibrationSubstrateError(ValueError):
    """Substrate inputs violate documented invariants. Indicates a bug
    upstream of the calibrator."""

class CalibrationConvergenceError(RuntimeError):
    """The mass-change diagnostic increased between iterations. Under
    EM theory this can only happen due to an implementation bug."""
```

The calibrator does not raise on data-shape pathologies (too few
regions, unstranded library, etc.) — the Bayesian priors handle those
naturally.

## 10. Logging

Single `logger.info` line at completion:

```
calibration: R=200k iters=6  rho_0=0.0012 phi=0.23 rho_d_BB=0.018 rho_r_BB=0.005 eps_s=8e-4
```

Per-iteration `logger.debug` with mass-change diagnostic.

No flag histograms. No fit_status enums. Rich diagnostics live in
`scripts/debug/dump_calibration_state.py`.

## 11. Backward compatibility

None. Legacy `CalibrationResult` had 28 fields, most diagnostic. New
schema has goal-derived outputs plus library hyperparameters plus
convergence diagnostics. No shim, no adapter. Tests reading legacy
fields are deleted.
