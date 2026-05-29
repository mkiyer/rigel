# Adaptive Prior v5 + v6 — Phased Implementation Plan

Date: 2026-05-26
Sources: `docs/prior/adaptive_prior_v5.md`, `docs/prior/adaptive_prior_v6.md`
Status: execution-ready

This document is the operational roadmap for landing v5 (the parameter-free
entropy-Dirichlet prior) with v6 (the single `rna_call_bias` dial) as a clean
break. It is organized into 8 phases. Each phase is independently runnable but
they must be applied in order: phases 1–2 add the new code, phase 3 swaps it
in, phases 4–6 delete the v3/v4 surface, phase 7 fixes tests, phase 8
regenerates goldens.

## Guiding Principles

- **No backwards compatibility.** No deprecation aliases. No legacy code paths.
- **Default `rna_call_bias = 0.5` ⇒ bit-identical to v5.** All v5 tests pass
  unchanged at the default.
- **Native EM ABI unchanged.** The prior is still
  `(alpha_gdna_add, alpha_rna_add, gdna_eff_len, enable_gdna)` per locus.
- **Calibration tensor is consumed as a contract.** No changes to
  `RegionCalibration`, `PriorMassDeconvolution`, or `latent_states.py`.
- **`_INTERNAL_RNA_LOWER_CI = 0.95`** replaces the user-facing
  `rna_lower_confidence` knob inside `strand_deconv.py`; the strand
  deconvolution algorithm itself is unchanged.
- **`_INTERNAL_GDNA_DENSITY_CI = 0.95`** replaces the user-facing
  `gdna_density_confidence` knob inside the calibration orchestrator; the
  calibration-iteration math is unchanged.
- **No stale user-facing prior/config names.** Any remaining internal use of
  removed confidence levels must be renamed to `internal_*` or hidden from
  `summary.json` so removed knobs do not reappear as apparent configuration.
- **Track deviations** in `/memories/session/adaptive_prior_v6_implementation.md`.

---

## Phase 1 — New helper module + unit tests

**Files added:**
- `src/rigel/calibration/adaptive_prior.py` (new)
- `tests/test_adaptive_prior.py` (new)

**Files touched:** none

**Pseudocode for `adaptive_prior.py`:**

```python
# Module constants (private)
MAX_ESS  = 3000.0
EPS_MASS = 1.0e-12

# Flag bits
PRIOR_NO_UNSPLICED_MASS = 0x1
PRIOR_STRUCTURAL_GATED  = 0x2
PRIOR_ESS_CAPPED        = 0x4
PRIOR_BIAS_APPLIED      = 0x8


@dataclass(frozen=True, slots=True)
class AdaptivePriorResult:
    alpha_gdna_add:    np.ndarray   # float64[L]  (final, after cap/gate/bias)
    alpha_rna_add:     np.ndarray   # float64[L]
    n_local:           np.ndarray   # float64[L, 2]
    n_other:           np.ndarray   # float64[L, 2]
    locus_weight:      np.ndarray   # float64[L]  = W_l
    shrink_weight:     np.ndarray   # float64[L]  = 1 - W_l
    region_weight:     np.ndarray   # float64[R]  = w_r
    global_counts:     np.ndarray   # float64[2]
    locus_unspliced:   np.ndarray   # float64[L]  = U_l
    n_regions_touched: np.ndarray   # int32[L]
    multi_locus_region_mass: np.ndarray  # float64[L], geometry diagnostic on U_r scale
    partial_coverage_region_mass: np.ndarray  # float64[L], geometry diagnostic on U_r scale
    unallocated_unspliced: float
    unallocated_weighted_unspliced: float
    rna_share_v5:      np.ndarray   # float64[L]  = α_r / (α_g + α_r) before bias
    rna_share_final:   np.ndarray   # float64[L]  = α_r' / (α_g' + α_r') after bias
    ess_final:         np.ndarray   # float64[L]  = α_g' + α_r'
    flags:             np.ndarray   # uint16[L]


def _validate_locus_ids(multi_loci: list[MultiLocus], n_loci: int) -> None:
    """Require contiguous `multi_locus_id` values because outputs are id-indexed."""
    ids = sorted(int(locus.multi_locus_id) for locus in multi_loci)
    if ids != list(range(n_loci)):
        raise ValueError("multi_loci must have contiguous ids 0..L-1")


def _validate_inputs(...):
    """Validate helper contracts before arithmetic."""
    # p_states: shape (R, N_STATES), finite, nonnegative within tiny tolerance,
    # row sums positive. Normalize rows defensively after clipping tiny roundoff.
    # unspliced_total: shape (R,), finite; negative values are not accepted.
    # has_gdna_candidate: shape (L,).
    # rna_call_bias: finite and 0.0 < value < 1.0.
    # max_ess: finite and > 0.


def _entropy_weight(p_states: np.ndarray) -> np.ndarray:
    """w_r = 1 - H(P_r)/log(N_STATES), with 0 log 0 := 0. Returns float64[R]."""
    from .latent_states import N_STATES

    p = np.clip(p_states.astype(np.float64, copy=False), 0.0, 1.0)
    row_sum = p.sum(axis=1)
    p = np.divide(p, row_sum[:, None], out=np.zeros_like(p), where=row_sum[:, None] > 0.0)
    # safe x log x: x*log(x) for x>0 else 0
    xlogx = np.zeros_like(p)
    positive = p > 0.0
    xlogx[positive] = p[positive] * np.log(p[positive])
    H = -xlogx.sum(axis=1)
    H_max = math.log(N_STATES)
    w = 1.0 - H / H_max
    return np.clip(w, 0.0, 1.0)


def _project_to_loci(
    *,
    region_arrays: RegionArrays,
    multi_loci: list[MultiLocus],
    region_values: dict[str, np.ndarray],   # name -> float64[R]
    n_loci: int,
) -> tuple[dict[str, np.ndarray], dict[str, float], np.ndarray]:
    """Generic s_lr projector. One geometric pass, multiple value arrays.
    Returns (per_locus, unallocated, n_regions_touched, multi_locus_region_mass,
    partial_coverage_region_mass).
    The s_lr fractions match _allocate_counts_by_geometry exactly:
        share_lr = (overlap_len / r_len) / sum_l (overlap_len / r_len)

    The two retained geometry mass diagnostics are measured on the
    `unspliced_total` scale (`U_r`), not the legacy `gdna_expected` scale.
    This is a deliberate v5 schema change: `U_r` is the mass the adaptive
    prior owns, while `gdna_expected` no longer exists as a prior input.
    """
    out = {name: np.zeros(n_loci, dtype=np.float64) for name in region_values}
    unalloc = {name: 0.0 for name in region_values}
    n_touched = np.zeros(n_loci, dtype=np.int32)
    multi_mass = np.zeros(n_loci, dtype=np.float64)
    partial_mass = np.zeros(n_loci, dtype=np.float64)
    # ... same geometry walk as existing _allocate_counts_by_geometry
    # but allocating every array in region_values in lockstep.
    # When a region touches multiple loci, add share * U_r to multi_mass[lid].
    # When a region is only partially covered by all loci, add share * U_r to
    # partial_mass[lid].


def compute_adaptive_prior(
    *,
    region_arrays: RegionArrays,
    multi_loci: list[MultiLocus],
    p_states: np.ndarray,            # float[R, 4]  (cast to float64 inside)
    unspliced_total: np.ndarray,     # float[R]
    has_gdna_candidate: np.ndarray,  # bool[L]
    rna_call_bias: float = 0.5,
    max_ess: float = MAX_ESS,
) -> AdaptivePriorResult:
    # 0. Validate and normalize contracts
    _validate_locus_ids(multi_loci, len(multi_loci))
    P, U, has_gdna_candidate, rna_call_bias, max_ess = _validate_inputs(...)

    # 1. Entropy weights
    w_r = _entropy_weight(P)                              # [R]

    # 2. q_r via state masks (collapse to 2-group)
    from .latent_states import STATE_IS_EXPRESSED
    rna_mask = np.asarray(STATE_IS_EXPRESSED, dtype=bool) # [4]
    gdna_mask = ~rna_mask
    q_gdna = P[:, gdna_mask].sum(axis=1)
    q_rna  = P[:, rna_mask].sum(axis=1)
    # normalize defensively
    q_sum = q_gdna + q_rna
    q_gdna = np.where(q_sum > 0, q_gdna / q_sum, 0.0)
    q_rna  = np.where(q_sum > 0, q_rna  / q_sum, 0.0)

    # 3. Regional pseudo-counts n_r = U_r * w_r * q_r
    n_gdna_r = U * w_r * q_gdna
    n_rna_r  = U * w_r * q_rna
    weighted_U_r = U * w_r

    # 4. Project to loci in ONE geometric pass
    L = len(multi_loci)
    projected, unalloc, n_touched, multi_mass, partial_mass = _project_to_loci(
        region_arrays=region_arrays,
        multi_loci=multi_loci,
        region_values={
            "U_l":         U,
            "N_l_gdna":    n_gdna_r,
            "N_l_rna":     n_rna_r,
            "U_l_weighted": weighted_U_r,
        },
        n_loci=L,
    )
    U_l           = projected["U_l"]
    N_l_gdna      = projected["N_l_gdna"]
    N_l_rna       = projected["N_l_rna"]
    U_l_weighted  = projected["U_l_weighted"]

    # 5. Locus information weight W_l
    W_l = np.where(U_l > EPS_MASS, U_l_weighted / np.maximum(U_l, EPS_MASS), 0.0)
    W_l = np.clip(W_l, 0.0, 1.0)
    shrink = 1.0 - W_l

    # 6. Global aggregate (one pass)
    N_global_gdna = float(N_l_gdna.sum())
    N_global_rna  = float(N_l_rna.sum())

    # 7. N_other = N_global - N_l (componentwise nonnegative by construction)
    N_other_gdna = np.maximum(N_global_gdna - N_l_gdna, 0.0)
    N_other_rna  = np.maximum(N_global_rna  - N_l_rna,  0.0)

    # 8. Combined alphas (pre-cap, pre-gate, pre-bias)
    alpha_g = N_l_gdna + shrink * N_other_gdna
    alpha_r = N_l_rna  + shrink * N_other_rna

    # 9. Cap by min(U_l, max_ess), preserve direction
    total = alpha_g + alpha_r
    cap = np.minimum(U_l, max_ess)
    capped_mask = total > cap
    scale = np.where(capped_mask & (total > 0), cap / np.maximum(total, EPS_MASS), 1.0)
    alpha_g = alpha_g * scale
    alpha_r = alpha_r * scale

    # 10. Structural gate + no-mass gate
    has_mass = U_l > EPS_MASS
    enable   = np.asarray(has_gdna_candidate, dtype=bool) & has_mass
    alpha_g  = np.where(enable, alpha_g, 0.0)
    alpha_r  = np.where(enable, alpha_r, 0.0)

    # 11. Record v5 share BEFORE the v6 dial
    T_v5 = alpha_g + alpha_r
    rna_share_v5 = np.divide(alpha_r, T_v5, out=np.zeros_like(T_v5), where=T_v5 > 0)

    # 12. v6 mass-conserving logit shift (only if non-default)
    bias_applied = np.zeros(L, dtype=bool)
    if rna_call_bias != 0.5:
        two_w   = 2.0 * rna_call_bias
        two_omw = 2.0 * (1.0 - rna_call_bias)
        s_v = rna_share_v5
        num = two_w * s_v
        den = num + two_omw * (1.0 - s_v)
        s_w = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
        # only apply where T_v5 > 0 (i.e., locus actually has prior mass)
        mass = T_v5 > 0
        alpha_g = np.where(mass, T_v5 * (1.0 - s_w), alpha_g)
        alpha_r = np.where(mass, T_v5 * s_w,         alpha_r)
        bias_applied = mass

    # 13. Final share
    T_final = alpha_g + alpha_r
    rna_share_final = np.divide(alpha_r, T_final, out=np.zeros_like(T_final), where=T_final > 0)

    # 14. Flags
    flags = np.zeros(L, dtype=np.uint16)
    flags |= np.where(~has_mass, PRIOR_NO_UNSPLICED_MASS, 0).astype(np.uint16)
    flags |= np.where(~np.asarray(has_gdna_candidate, dtype=bool),
                      PRIOR_STRUCTURAL_GATED, 0).astype(np.uint16)
    flags |= np.where(capped_mask, PRIOR_ESS_CAPPED, 0).astype(np.uint16)
    flags |= np.where(bias_applied, PRIOR_BIAS_APPLIED, 0).astype(np.uint16)

    return AdaptivePriorResult(
        alpha_gdna_add=alpha_g, alpha_rna_add=alpha_r,
        n_local=np.column_stack([N_l_gdna, N_l_rna]),
        n_other=np.column_stack([N_other_gdna, N_other_rna]),
        locus_weight=W_l, shrink_weight=shrink,
        region_weight=w_r,
        global_counts=np.array([N_global_gdna, N_global_rna]),
        locus_unspliced=U_l,
        n_regions_touched=n_touched,
        multi_locus_region_mass=multi_mass,
        partial_coverage_region_mass=partial_mass,
        unallocated_unspliced=float(unalloc["U_l"]),
        unallocated_weighted_unspliced=float(unalloc["U_l_weighted"]),
        rna_share_v5=rna_share_v5, rna_share_final=rna_share_final,
        ess_final=T_final, flags=flags,
    )
```

**Implementation notes fixed during plan review:**
- The helper output arrays are indexed by `multi_locus_id`, not input-list
  position. Validate ids are contiguous (`0..L-1`) so shuffled-list tests are
  meaningful and production assumptions fail loudly if violated.
- `p_states` is validated against `N_STATES == 4` and row-normalized after
  clipping tiny floating-point roundoff. Bad shape, nonfinite values, or
  materially negative probabilities raise `ValueError`.
- `unspliced_total` is not silently clipped. Negative or nonfinite values
  indicate an upstream conservation bug and raise `ValueError`.
- `PRIOR_ESS_CAPPED` follows the v5 definition exactly: it is set when
  `sum(alpha_l)` before scaling exceeds `min(U_l, MAX_ESS)`, even if the same
  locus is later zeroed by the no-mass or structural gate.

**Unit tests for `tests/test_adaptive_prior.py` (17 cases):**

v5 core (10 from v5 §9.1):
1. `H(P_r)`: one-hot ⇒ `w_r=1`; uniform ⇒ `w_r=0`; monotone in between.
2. `n_r`: `U_r=0` or `w_r=0` ⇒ zero; certain-RNA ⇒ `(0, U_r)`; certain-gDNA ⇒ `(U_r, 0)`.
3. Geometry: disjoint regions sum exactly per locus; split region splits proportionally.
4. `N_other ≥ 0` componentwise.
5. `W_l=1` + nonzero `N_l` ignores `N_other`.
6. `W_l=0` + nonzero `N_other` adopts full `N_other` direction.
7. Cap: rescales to exactly `min(U_l, MAX_ESS)`, preserves direction.
8. Structural gate ⇒ both alphas exactly 0.
9. Single-locus sample: `N_other=0` ⇒ `alpha = N_l`.
10. Locus order invariance (shuffle test).

v6 dial (7 from v6 §6.1):
11. `ω=0.5` ⇒ bit-identical to v5.
12. `ω→0` (e.g. 1e-4) ⇒ `α_r→0`, `α_g→T_v5`.
13. `ω→1` (e.g. 1-1e-4) ⇒ `α_g→0`, `α_r→T_v5`.
14. Mass conservation: `α_g'+α_r' == T_v5` to float tol.
15. Structural gate still wins under any `ω`.
16. Cap still wins: `α_g'+α_r' ≤ min(U_l, MAX_ESS)`.
17. Monotonicity: for fixed locus with `T>0`, `α_r'(ω)` strictly increasing on `(0,1)`.

Phase exit criterion: `pytest tests/test_adaptive_prior.py -v` all green.

---

## Phase 2 — `EMConfig.rna_call_bias` field

**Files touched:**
- `src/rigel/config.py`

**Changes:**
- Add `rna_call_bias: float = 0.5` field to `EMConfig`.
- Extend `__post_init__` to validate `0.0 < rna_call_bias < 1.0`.

(We do this in phase 2 so phase 3 wiring can read it. The deletions of the old
EMConfig fields happen in phase 5, after wiring is in place.)

Phase exit: `pytest tests/test_adaptive_prior.py` still passes; import of
`EMConfig()` works.

---

## Phase 3 — Wire `assemble_priors()` to the new helper

**Files touched:**
- `src/rigel/calibration/prior.py`

**Important review fix:** phase 3 must not leave `n_regions_touched`,
`multi_locus_region_mass`, or `partial_coverage_region_mass` as local
variables in `prior.py`; those diagnostics now come from
`AdaptivePriorResult`. This avoids running two geometry passes and avoids the
undefined-variable bug in the first draft of this plan.

**Pseudocode:**

```python
# Trim PriorTable to the v5+v6 schema
@dataclass(frozen=True, slots=True)
class PriorTable:
    # Native EM inputs (unchanged ABI)
    alpha_gdna_add:           np.ndarray  # float64[L]
    alpha_rna_add:            np.ndarray  # float64[L]
    gdna_eff_len:             np.ndarray  # float64[L]
    gdna_eff_len_unweighted:  np.ndarray  # float64[L]
    gdna_em_exposure_weight:  np.ndarray  # float64[L]
    enable_gdna:              np.ndarray  # uint8[L]
    # v5/v6 diagnostics
    prior_locus_weight:       np.ndarray  # float64[L]  W_l
    prior_shrink_weight:      np.ndarray  # float64[L]  1-W_l
    prior_n_local_gdna:       np.ndarray  # float64[L]
    prior_n_local_rna:        np.ndarray  # float64[L]
    prior_n_other_gdna:       np.ndarray  # float64[L]
    prior_n_other_rna:        np.ndarray  # float64[L]
    prior_ess_final:          np.ndarray  # float64[L]
    prior_rna_share_v5:       np.ndarray  # float64[L]
    prior_rna_share_final:    np.ndarray  # float64[L]
    prior_unspliced_total:    np.ndarray  # float64[L]  U_l (kept; cheap and used)
    prior_flags:              np.ndarray  # uint16[L]
    # Geometry/coverage diagnostics (kept; cheap and informative)
    n_regions_touched:        np.ndarray  # int32[L]
    multi_locus_region_mass:  np.ndarray  # float64[L]
    partial_coverage_region_mass: np.ndarray  # float64[L]
    # Sample-level diagnostics (now from helper's globals)
    global_n_gdna:            float = 0.0
    global_n_rna:             float = 0.0
    prior_policy_name:        str = "entropy_dirichlet_v5_v6"
    prior_max_ess:            float = 3000.0
    rna_call_bias:            float = 0.5

    def to_summary_dict(self) -> dict:
      # Return the docs/prior v5+v6 `prior_policy` payload, including:
      # name, max_ess, rna_call_bias, global_counts, region/locus weight
      # summaries, rna_share_shift_p50/p90, n_loci_total,
      # n_loci_with_prior_mass, n_loci_structural_gated, and flag_histogram.


def assemble_priors(*, multi_loci, em_data, index, calibration, em_config=None) -> PriorTable:
    if em_config is None:
        em_config = EMConfig()
    n_loci = len(multi_loci)
    if n_loci == 0:
        return _empty_prior_table()
    # ... existing validation (region_df present, region_calibration present)
    region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    rc = calibration.region_calibration
    p_states = np.asarray(rc.p_states, dtype=np.float64)               # [R, 4]
    unspliced_total = np.asarray(rc.prior_mass.unspliced_total, dtype=np.float64)  # [R]
    has_gdna_candidate = np.array(
        [enable_gdna_for_multilocus(l, em_data) for l in multi_loci], dtype=bool,
    )

    result = compute_adaptive_prior(
        region_arrays=region_arrays,
        multi_loci=multi_loci,
        p_states=p_states,
        unspliced_total=unspliced_total,
        has_gdna_candidate=has_gdna_candidate,
        rna_call_bias=float(em_config.rna_call_bias),
    )

    # gdna_eff_len + exposure (kept from current implementation)
    gdna_eff_len, gdna_eff_len_unweighted, gdna_em_exposure_weight = _compute_gdna_eff_len(
        multi_loci=multi_loci, index=index, calibration=calibration,
        region_arrays=region_arrays,
    )

    return PriorTable(
        alpha_gdna_add=result.alpha_gdna_add,
        alpha_rna_add=result.alpha_rna_add,
        gdna_eff_len=gdna_eff_len,
        gdna_eff_len_unweighted=gdna_eff_len_unweighted,
        gdna_em_exposure_weight=gdna_em_exposure_weight,
        enable_gdna=has_gdna_candidate.astype(np.uint8),
        prior_locus_weight=result.locus_weight,
        prior_shrink_weight=result.shrink_weight,
        prior_n_local_gdna=result.n_local[:, 0],
        prior_n_local_rna=result.n_local[:, 1],
        prior_n_other_gdna=result.n_other[:, 0],
        prior_n_other_rna=result.n_other[:, 1],
        prior_ess_final=result.ess_final,
        prior_rna_share_v5=result.rna_share_v5,
        prior_rna_share_final=result.rna_share_final,
        prior_unspliced_total=result.locus_unspliced,
        prior_flags=result.flags,
        prior_policy_name="entropy_dirichlet_v5_v6",
        prior_max_ess=MAX_ESS,
        rna_call_bias=float(em_config.rna_call_bias),
        n_regions_touched=result.n_regions_touched,
        multi_locus_region_mass=result.multi_locus_region_mass,
        partial_coverage_region_mass=result.partial_coverage_region_mass,
        global_n_gdna=float(result.global_counts[0]),
        global_n_rna=float(result.global_counts[1]),
    )
```

**Symbols deleted from `prior.py`:**
- `GroupedPriorCounts` dataclass
- `_compute_grouped_prior_counts()`
- All `prior_budget*`, `prior_gdna_share_*`, `gdna_prior_count_em`,
  `gdna_prior_density`, `rna_expected_count`, `gdna_expected_count`,
  `gdna_upper_count`, `prior_mass_conservation_error`,
  `prior_allocated_fraction`, `count_allocation_mode`,
  `count_allocation_fallback`, `n_units_used_for_diagnostics`,
  `unallocated_*` fields from `PriorTable`.
- `COUNT_ALLOC_*` constants if no longer used after pipeline strip.

**Summary wiring:**
- `PriorTable.to_summary_dict()` is the canonical `prior_policy` payload.
- `cli._write_results()` should expose it at top-level
  `summary["prior_policy"]` and may also leave it under
  `summary["calibration"]["prior_table"]` if that is the least invasive path.
  The top-level key is required by the v5/v6 docs.

**Compatibility deferral:** `pipeline.py` and `estimator.py` still reference
the old fields. They WILL break here — we intentionally let them break and
fix them in phase 6. (Phase 3 must be followed quickly by phase 6 to restore
test runnability. The split exists for review clarity.)

Phase exit: `pytest tests/test_adaptive_prior.py` still green; `from
rigel.calibration.prior import PriorTable, assemble_priors` imports cleanly.
The rest of the test suite WILL fail import — that is expected.

---

## Phase 4 — Strand-deconv & calibration knob removal

**Files touched:**
- `src/rigel/calibration/strand_deconv.py`
- `src/rigel/calibration/_orchestrator.py`
- `src/rigel/calibration/_result.py`
- `src/rigel/config.py`

**Strand deconv:**
1. Add module constant `_INTERNAL_RNA_LOWER_CI = 0.95` at top of
   `strand_deconv.py`.
2. Remove `rna_lower_confidence: float` from every function signature inside
   that module (`deconvolve_regions_by_strand` and helpers). All call sites
   inside the module now reference `_INTERNAL_RNA_LOWER_CI` directly.
3. Search for any other module that reads this constant (likely
   `_orchestrator.py`); make their internal references use the same constant.
4. Same treatment for `gdna_density_confidence` if a density module uses it
  internally. In the current code this is consumed by `_orchestrator.py`, so
  define `_INTERNAL_GDNA_DENSITY_CI = 0.95` there and pass it to
  `run_calibration_iteration(confidence=...)`.

**Orchestrator:**
- Remove `rna_lower_confidence` and `gdna_density_confidence` kwargs from the
  orchestrator entry point and from any helper it forwards them to.
- Continue passing explicit internal constants to `estimate_kappa_d()`,
  `deconvolve_compartments_by_strand()`, and `run_calibration_iteration()` so
  defaults do not drift silently.

**Calibration result schema:**
- Remove `CalibrationResult.rna_lower_confidence` and the
  `calibration_config.rna_lower_confidence` summary payload.
- Remove `build_calibration_result(..., rna_lower_confidence=...)`.
- If strand-deconvolution result objects need to retain the internal CI for
  debugging, rename the field and summary key to `internal_rna_lower_ci`.
  Otherwise remove it entirely. Do not leave `rna_lower_confidence` in
  `summary.json`.

**Config:**
- Delete `CalibrationConfig.rna_lower_confidence` and
  `CalibrationConfig.gdna_density_confidence` fields and their
  `__post_init__` validation blocks.

Phase exit: `python -c "from rigel.config import CalibrationConfig;
print(CalibrationConfig())"` works; `from rigel.calibration.strand_deconv
import deconvolve_regions_by_strand` works. Strand deconv unit tests still
pass mathematically (the constant equals the previous default 0.95).

---

## Phase 5 — `EMConfig` deletions

**Files touched:**
- `src/rigel/config.py`

**Changes:**
- Delete `EMConfig.aggregate_prior_strength`,
  `EMConfig.aggregate_prior_edge_count`,
  `EMConfig.aggregate_prior_max_count`,
  `EMConfig.gdna_prior_logit_bias`.
- Delete the corresponding `__post_init__` validation loop.
- Keep `rna_call_bias` (already added in phase 2).

Phase exit: `EMConfig()` instantiates with only the new field set; no
references to removed fields remain in `src/rigel/`.

---

## Phase 6 — Pipeline / estimator / CLI cleanup

**Files touched:**
- `src/rigel/pipeline.py`
- `src/rigel/estimator.py`
- `src/rigel/cli.py`

**`pipeline.py`:**
- `_run_em_for_locus` signature: drop deprecated kwargs
  (`prior_budget_raw`, `prior_budget`, `prior_gdna_share_raw`,
  `prior_gdna_share_biased`, `gdna_prior_count_em`); thread through
  `prior_locus_weight`, `prior_shrink_weight`, `prior_n_local_*`,
  `prior_n_other_*`, `prior_ess_final`, `prior_rna_share_v5`,
  `prior_rna_share_final`, `prior_flags`.
- Locus-results dict: replace deprecated keys with new ones.
- L1029 area: pass new `PriorTable` fields through.
- L1136 area: remove `rna_lower_confidence` / `gdna_density_confidence`
  kwargs to orchestrator call.

**`estimator.py`:**
- `loci_df` builder (L815–928): replace old prior columns with the v5+v6
  schema. New columns:
  ```
  alpha_gdna_add, alpha_rna_add,
  prior_locus_weight, prior_shrink_weight,
  prior_n_local_gdna, prior_n_local_rna,
  prior_n_other_gdna, prior_n_other_rna,
  prior_ess_final, prior_rna_share_v5, prior_rna_share_final,
  prior_unspliced_total, prior_flags,
  gdna_eff_len, gdna_eff_len_unweighted, gdna_em_exposure_weight,
  enable_gdna, n_regions_touched,
  multi_locus_region_mass, partial_coverage_region_mass
  ```
- Remove the derived legacy `gdna_prior` rate column along with
  `gdna_prior_count` / `gdna_prior_count_em`. It is a legacy alias for gDNA
  pseudocount strength and conflicts with the no-alias cutover.

**`cli.py`:**
- Delete `_PARAM_SPECS` entries for `rna_lower_confidence`,
  `gdna_density_confidence`.
- Delete corresponding `--rna-lower-confidence`,
  `--gdna-density-confidence` argparse entries. The aggregate-prior and
  logit-bias fields are currently config-only, not CLI flags, so there is
  nothing to delete from argparse for them unless new references appear.
- Add `--rna-call-bias FLOAT` flag (range `(0, 1)`, default `0.5`) and a
  matching `_ParamSpec` mapped to `EMConfig.rna_call_bias`.
- YAML loader: do **not** overload `_REMOVED_QUANT_CONFIG_KEYS`, which is a
  scan-key replacement map with a different error message. Add a separate
  `_REMOVED_PRIOR_CONFIG_KEYS` set for the seven removed prior/config keys
  (`rna_confidence`, `rna_lower_confidence`, `gdna_density_confidence`, the
  four removed `EMConfig` aggregate/logit fields) and raise the v5 migration
  error:
  `Configuration option '{name}' was removed in adaptive prior v5. The grouped
  EM prior is now parameter-free. See docs/prior/adaptive_prior_v5.md.`
- YAML shape: v6 docs show `quant: {rna_call_bias: 0.5}` while the current CLI
  resolver accepts flat quant keys. Update the resolver to accept both by
  flattening a top-level `quant` mapping before validation, and update
  `_write_config_yaml()` to emit the documented `quant:` form. This is a YAML
  shape cleanup, not a prior-backcompat mode.
- `cli._write_results()`: add top-level `summary["prior_policy"]` from
  `result.calibration.prior_table.to_summary_dict()` when present.

Phase exit: `rigel quant --help` works; CLI shows `--rna-call-bias`; no
deprecated flags surface; `python -c "import rigel"` clean; `pytest
tests/test_cli.py -v` is the gate.

---

## Phase 7 — Test fixes

**Files touched / created / deleted:**

| File | Action |
|------|--------|
| `tests/test_rna_lower_confidence.py` | delete |
| `tests/test_calibration_prior.py` | rewrite for v5+v6 schema |
| `tests/test_bayesian_prior_acceptance.py` | rewrite or delete; replace with adaptive_prior acceptance subset |
| `tests/test_pipeline_wiring.py` | update PriorTable construction + EM kwargs |
| `tests/test_strand_deconv.py` | drop `rna_lower_confidence=` kwarg from all calls |
| `tests/test_calibrate.py` | drop kwargs as above |
| `tests/test_calibration_iteration.py` | drop kwargs as above |
| `tests/test_calibration_result.py` | drop kwargs / removed fields |
| `tests/test_latent_states.py` | check for stale references |
| `tests/test_estimator.py`, `tests/test_pipeline_smoke.py`, etc. | update column expectations |

**Process per test file:**
1. Run the test file alone.
2. Triage failures: import error → rename/remove; assertion on removed
   column → swap to new column or delete assertion.
3. Verify no test masks a removed knob by hardcoding `rna_lower_confidence`.

**New tests added:**
- CLI: `--rna-call-bias 0.5` accepted; `0.0` / `1.0` rejected; YAML
  round-trip preserves the field.
- CLI: removed-keys YAML produces the v5 migration error.

Phase exit: `pytest tests/ -v --ignore=tests/test_golden_output.py` is green.

---

## Phase 8 — Golden regeneration

**Files touched:**
- `tests/golden/**` (TSV / feather / JSON)
- `docs/newcal/new_cal_implementation_log.md` (append v5+v6 entry)

**Process:**
1. Inspect golden TSV headers to enumerate affected files:
   `grep -l 'gdna_prior_count_em\|prior_budget\|prior_gdna_share_biased' tests/golden/`.
2. Run `pytest tests/test_golden_output.py -v --update-golden`.
3. Diff the regenerated files. For each, verify:
   - `alpha_gdna_add`, `alpha_rna_add` are nonzero only where `enable_gdna=1`.
   - `prior_ess_final ≤ min(prior_unspliced_total, 3000)`.
   - `prior_rna_share_v5 == prior_rna_share_final` (default `ω=0.5`).
   - `prior_flags` bit 0x8 unset (default `ω=0.5`).
4. Spot-check a single high-prior locus by hand against the formula:
   `α = N_l + (1-W_l)(N_global - N_l)`, capped.
5. Commit goldens; append a dated entry to
   `docs/newcal/new_cal_implementation_log.md`.

Phase exit: `pytest tests/ -v` is fully green (modulo the pre-existing
`test_calibration.py::TestStrandLLR::test_biased_toward_ss_favors_rna`
failure noted in copilot-instructions).

---

## Execution Order Summary

| Phase | Adds | Removes | Test gate |
|-------|------|---------|-----------|
| 1 | adaptive_prior.py + unit tests | — | `test_adaptive_prior.py` |
| 2 | `EMConfig.rna_call_bias` | — | import smoke |
| 3 | New `PriorTable` schema | v3 grouped prior code | `test_adaptive_prior.py` |
| 4 | `_INTERNAL_*_CI` constants | strand/density confidence knobs | `test_strand_deconv.py` |
| 5 | — | 4 `EMConfig` knobs | import smoke |
| 6 | `--rna-call-bias` | deprecated CLI flags, removed-key rejection | `test_cli.py` |
| 7 | new CLI/YAML tests | obsolete tests | full suite minus goldens |
| 8 | — | — | full suite |

## Rollback Strategy

Phases 1–2 are purely additive and safe to leave in if rollback is needed.
Phases 3–6 are tightly coupled: if any one is reverted, all must be reverted
together (they collectively constitute the v3→v5 swap). Phases 7–8 are
test/golden hygiene and can be reverted independently.

## Deviation Log

Track deviations in
`/memories/session/adaptive_prior_v6_implementation.md` under the
"Deviations from plan" section. Each deviation must record:
- which phase
- which file
- what was supposed to happen
- what actually happened, and why
