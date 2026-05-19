# v4.2 Implementation Plan

**Status**: implementation plan, accepted.
**Date**: 2026-05-19
**Supersedes**: `gdna_exposure_plan_v4.1.md`, `gdna_exposure_plan_v4.1_implementation.md`.

v4.2 is v4.1 with three corrections from the dry-run:

1. The per-class signal attenuator is replaced by **per-region
   empirical-Bayes shrinkage** of `ρ̂_r` toward a global mean
   (Gamma-Poisson conjugate). Class membership no longer enters the
   weight calculation at all.
2. Containment iterates exons in **genomic order regardless of
   strand**. The midpoint-mass integral over valid fragment
   placements is strand-invariant.
3. `ScoredFragments.genomic_midpoint` is **deleted** with the
   numerator path. It has no other consumer.

Everything else from v4.1 stands. Column names follow the existing
`PriorTable` convention (P2 accepted).

---

## 1. Concept summary

Two rules, one engine, one global reference, per-region shrinkage.

- **Engine**: `integrate_exposure_over_midpoints(windows, exposure)` —
  a single piecewise-linear integrator of `A(x) = ρ(x) / ρ_ref` over
  per-ref half-open intervals. Used in every length calculation.
- **Overlap rule** (`weighted_eff_len_overlap`): for gDNA loci. For
  each `ℓ`, shift the merged genomic footprint by `+ℓ/2` and
  integrate A. Sum over `ℓ` with the FL pmf.
- **Containment rule** (`weighted_eff_len_contained`): for
  transcripts and nRNA spans. For each `ℓ`, restrict to valid
  fragment starts in transcript coordinates, project the **midpoint**
  to genome, integrate A over the resulting per-block genomic windows.
  Strand-independent (see §3).
- **Global reference** `ρ_ref`: weighted-Q95 of `ρ̂_r` over **all**
  regions in `RegionArrays`, weighted by `E_r`. One number, period.
- **Shrinkage**: each `ρ̂_r` is replaced by its Gamma-Poisson
  posterior mean before dividing by `ρ_ref`. Regions with small
  `E_r` are pulled toward the global mean. No class structure
  enters this step.
- **Denominator-only EM**: numerator (per-unit `log A`) deleted. EM
  integrates A only through `L_g^A` (gDNA), `L_t^A` (mRNA), and
  `L_s^A` (nRNA span) — all from the same engine.
- **Prior rescaling**: `η_g_em = η_g · L_g^A / L_g` so the calibrated
  Poisson density `η_g / L_g` is preserved.

---

## 2. Per-region shrinkage (replaces v4.1 per-class signal)

### 2.1 Why per-region

Uncertainty in `ρ̂_r = c_r / E_r` is intrinsic to the region: tiny
`E_r` regions have wide credible intervals on `ρ̂_r` no matter what
"class" they happen to belong to. Lumping tens of thousands of
heterogeneous regions under one class-wide attenuator throws away
the natural per-region information about which `ρ̂_r` values to
trust.

### 2.2 Model

Treat region counts as Gamma-Poisson:

$$
\rho_r \sim \mathrm{Gamma}(\alpha, \beta), \qquad
c_r \mid \rho_r \sim \mathrm{Poisson}(\rho_r \cdot E_r).
$$

The prior is shared across all regions in the partition; the
posterior is conjugate:

$$
\rho_r \mid c_r, E_r \sim \mathrm{Gamma}(\alpha + c_r, \beta + E_r),
$$

with posterior mean

$$
\tilde{\rho}_r \;=\; \frac{\alpha + c_r}{\beta + E_r}
\;=\; w_r \cdot \hat{\rho}_r \;+\; (1 - w_r) \cdot \bar{\rho},
\quad w_r = \frac{E_r}{\beta + E_r}, \quad \bar{\rho} = \frac{\alpha}{\beta}.
$$

Interpretation:

- `E_r` large ⇒ `w_r → 1` ⇒ trust `ρ̂_r`.
- `E_r` small ⇒ `w_r → 0` ⇒ shrink toward global mean `ρ̄`.
- `β` is the prior precision in units of exposure; it controls how
  aggressively small-`E_r` regions are pulled in.

### 2.3 Hyperparameter estimation (empirical Bayes)

Estimate `(α, β)` from the marginal `c_r` distribution by
method-of-moments:

- Marginal mean: `E[c_r] = α/β · E_r`.
- Marginal variance: `Var(c_r) = α/β · E_r · (1 + E_r/β)`.

Using `E_r`-weighted aggregates over all regions:

1. `ρ̄ = Σ_r c_r / Σ_r E_r` (global mean rate). Sets `α/β = ρ̄`.
2. Excess Poisson variance:
   `s² = Σ_r E_r · (ρ̂_r − ρ̄)² / Σ_r E_r − ρ̄ / Σ̄_E`,
   where `Σ̄_E = Σ_r E_r / R` is the mean exposure.
3. `β = ρ̄ / max(s², ε)`; `α = ρ̄ · β`.

If `s² ≤ 0` (data are pure Poisson around `ρ̄`, no real exposure
variation): set `β = +∞` ⇒ every `w_r = 0` ⇒ `Ã_r ≡ 1` everywhere.
This is the principled auto-uniform.

### 2.4 Final exposure field

```
Ã_r = ρ̃_r / ρ_ref    (clipped log to ±LOG_RHO_CLIP_NATS)
A_r = exp(clipped log Ã_r)
```

`ρ_ref` is the global weighted-Q95 of the **unshrunk** `ρ̂_r`.
Quantile uses unshrunk values so the reference does not depend on
the prior; this keeps `Ã_r` stable as `β` is re-estimated.

### 2.5 What goes away

- `signal_per_class` dict.
- `rho_ref_per_class` dict and the `max()` collapse.
- The `max(signal · raw, LOG_A_FLOOR)` clamp (replaced by hard
  `±LOG_RHO_CLIP_NATS` clipping after shrinkage).

### 2.6 What gets exposed in `summary.json`

- `rho_ref_global` (the single Q95).
- `eb_prior` dict: `{alpha, beta, mean_rho, n_regions, mean_exposure}`.
- `rho_ref_class` (per-class Q95 of unshrunk `ρ̂_r`) as a backward-
  compat diagnostic only; not used to compute `A_r`.
- Per-region histogram of `w_r` (e.g. 10 bins) for QC.

---

## 3. Containment is strand-independent

Claim: for any transcript with exon blocks `B_1, …, B_N` in genomic
order, exonic length `L_t = Σ |B_k|`, fragment length `ℓ ≤ L_t`,

$$
L_t^A(\ell) \;=\; \int_0^{L_t - \ell} A\!\left(\mathrm{proj}_t(t_m)\right) \, dt_m
$$

is invariant under reversing the projection orientation (i.e. under
strand flip). Proof: the projection is a measure-preserving bijection
between `[0, L_t)` and the disjoint union of `[B_k]`. Reversing it
permutes the same set of genomic positions; the integral of `A` over
that set, weighted by uniform `dt_m`, is identical. Only the
**mapping from `t_m` to genomic position** changes; the **integral**
does not.

Practical consequence: the containment routine iterates exons in
**genomic order** unconditionally. No strand field is needed. The
function signature drops the `strand` parameter that was proposed in
v4.1.

Implementation sketch (Python, vectorized, per `ℓ`):

```
# blocks: int64[N, 2] genomic [start, end), genomic-ordered.
# Cumulative exonic positions C[k] = Σ_{j<k} |B_j|.
# For each ℓ, valid t_m ∈ [0, L_t - ℓ).
# Per block B_k = [g_lo, g_hi), the t_m range that maps into B_k is
#   t_m ∈ [C[k], C[k+1]) ∩ [0, L_t - ℓ).
# But we want the midpoint of the *fragment*, not the start.
# Fragment start s ∈ [0, L_t - ℓ), midpoint t_m = s + ℓ/2.
# Midpoint mass over fragments whose midpoint lies in block B_k:
#   t_m ∈ [C[k], C[k+1)) ∩ [ℓ/2, L_t - ℓ/2),
# which projects to genomic window
#   [g_lo + (max(C[k], ℓ/2) - C[k]),
#    g_lo + (min(C[k+1], L_t - ℓ/2) - C[k])).
# Pass these windows to integrate_exposure_over_midpoints.
```

Integer ℓ/2 with half-open clipping; no floating arithmetic in the
window bounds.

---

## 4. Per-ref prefix-sum exposure (performance)

`RegionalGdnaExposure` lazily builds, per ref:

- `_starts[ref]: int64[R+1]` — region starts plus sentinel at end.
- `_csum[ref]:   float64[R+1]` — `_csum[ref][i] = Σ_{j<i} A_j · (end_j − start_j)`.

`weighted_length_on_ref(ref, lo, hi)` becomes:

1. `i = searchsorted(_starts[ref], lo, side='right') − 1`
2. `j = searchsorted(_starts[ref], hi, side='right') − 1`
3. Interior mass = `_csum[ref][j] − _csum[ref][i+1]`
4. Add partial contributions of region `i` and region `j`.
5. Add identity tail for any `[lo,hi)` portion beyond the last region.

`weighted_length_on_ref_batch(ref, lo_arr, hi_arr)` does the two
searchsorted calls on full arrays and computes interior masses with
vector arithmetic. This is what containment uses.

Total cost for 200k transcripts × ~40 truncated `ℓ` × ~10 blocks:
two vectorized searchsorted per `(ℓ, transcript)`, dominated by
NumPy overhead at ~8M outer calls. With per-transcript batching
(one searchsorted call across all `(ℓ, block)` pairs for that
transcript), total Python-level calls drop to ~200k. Comfortably
inside the G9 +10% budget.

Native port: deferred. Revisit only if profiling shows containment
remains the hot path after the prefix-sum + batching change.

---

## 5. Engine and rule signatures

```python
# src/rigel/calibration/_exposure.py

def integrate_exposure_over_midpoints(
    windows: list[tuple[int, int, int]],   # (ref, g_lo, g_hi)
    exposure: "RegionalGdnaExposure",
) -> float: ...

def weighted_eff_len_overlap(
    footprint: list[tuple[int, int, int]],   # (ref, g_lo, g_hi), no merge required
    ref_lengths: np.ndarray,
    fl: FragmentLengthModel,
    exposure: "RegionalGdnaExposure",
    *,
    min_value: float = 1.0,
) -> float: ...

def weighted_eff_len_contained(
    footprint: list[tuple[int, int, int]],   # genomic-ordered exon blocks
    fl: FragmentLengthModel,
    exposure: "RegionalGdnaExposure",
    *,
    min_value: float = 1.0,
) -> float: ...
```

No `strand` argument (per §3). Existing
`weighted_gdna_eff_len_for_loci` becomes a thin adapter from
`MultiLocus.loci` to `weighted_eff_len_overlap`.

Per-transcript driver (`_exposure_lengths.py` or appended to
`_exposure.py`):

```python
def weighted_transcript_eff_lens(
    index: TranscriptIndex,
    rna_fl: FragmentLengthModel,
    exposure: "RegionalGdnaExposure",
) -> np.ndarray: ...
```

Loops over `index._t_exon_intervals` (which already covers synthetic
nRNA single-block entries), calls `weighted_eff_len_contained` per
transcript. Returns `effective_lengths_em` shape `(n_t,)`.

---

## 6. Estimator EM/output split

(`P1` in dry-run, unchanged.)

- `TranscriptGeometry.effective_lengths_em: np.ndarray | None = None`.
- `AbundanceEstimator` stores `_t_eff_len_em` and `_t_eff_len_output`.
- Public `effective_lengths` property → `_t_eff_len_output`
  (no breaking change for external readers).
- New `effective_lengths_em` property for EM internals.
- Uniform mode: both arrays are the same object.

---

## 7. `PriorTable` and `loci.feather` schema

Accepted naming from dry-run P2:

| Field | Meaning |
|---|---|
| `gdna_eff_len` | `L_g^A` — EM-consumed weighted length. |
| `gdna_eff_len_unweighted` | `L_g` — raw overlap effective length. |
| `gdna_prior_count` | `η_g` — calibration's Poisson prior count. |
| `gdna_prior_count_em` | `η_g · L_g^A / L_g` — what EM actually consumes. |
| `gdna_em_exposure_weight` | `L_g^A / L_g` (diagnostic). |

`gdna_prior_count_regional` and the `_expected_gdna_count_regional_diag`
helper are deleted. They are confusing and not used by anything that
v4.2 ships.

For `quant.feather` and `nrna_quant.feather`:

| New column | Meaning |
|---|---|
| `em_effective_length` | `L_t^A` (`L_s^A` for nRNA). |
| `em_exposure_weight` | `L_t^A / L_t`. |

Existing `effective_length` and `TPM` columns are unchanged
(they continue to use the **unweighted** length).

---

## 8. Numerator deletion

- Delete `pipeline._apply_unit_gdna_weights` and its call site.
- Delete `RegionalWeightApplicationStats` dataclass and its
  `summary.json["regional_exposure"]["application"]` block.
- Delete `RegionalGdnaExposure.log_weights_for_positions` (no
  remaining caller).
- Delete `ScoredFragments.genomic_midpoint` (no remaining
  consumer). Update `bam_scanner.cpp` to stop emitting it.
- Update `test_bam_tag_parsing.py` and friends that read
  `genomic_midpoint` (if any).

This is the right time to drop these — keeping them as stubs
buys nothing and clutters the schema.

---

## 9. Tests

Files to add:

- `tests/test_regional_exposure_global_ref.py` — global Q95 of
  unshrunk `ρ̂_r`; cross-class equivalence; uniform-override sentinel.
- `tests/test_regional_exposure_shrinkage.py` — Gamma-Poisson MoM
  recovery on simulated data; `w_r` monotone in `E_r`;
  `s² ≤ 0` → `β = ∞` ⇒ `A_r ≡ 1`.
- `tests/test_exposure_overlap_engine.py` — parity with
  `gdna_eff_len_for_loci` under uniform `A`; merge stress.
- `tests/test_exposure_contained_engine.py` — parity with
  `compute_all_transcript_eff_lens` under uniform `A`;
  strand-flip invariance (sentinel for §3); constant-`A_0` law
  with full-coverage fixture.
- `tests/test_estimator_em_output_split.py` — public
  `effective_length` column stays unweighted; EM denominator
  differs only when `regional_exposure.mode == "regional"`.
- `tests/test_regional_exposure_locus.py` — synthetic two-region
  locus end-to-end.

Existing pipeline / golden tests must keep passing in uniform mode.

---

## 10. Implementation order

1. **Per-region shrinkage** — `_regional_exposure.py::build`:
   global `ρ_ref` + MoM `(α, β)` + posterior `ρ̃_r`. Drop per-class
   signal and per-class `ρ_ref`. New summary fields. Tests under
   `test_regional_exposure_global_ref.py` and
   `test_regional_exposure_shrinkage.py`.
2. **Prefix-sum exposure** — `_regional_exposure.py`: lazy
   `_starts[ref]`, `_csum[ref]`; new
   `weighted_length_on_ref_batch`. Existing
   `weighted_length_on_ref` switches to the fast path. Internal
   parity test against the old implementation.
3. **Engine and overlap rule** — `_exposure.py`:
   `integrate_exposure_over_midpoints`, `weighted_eff_len_overlap`.
   Refactor `weighted_gdna_eff_len_for_loci` to delegate. Parity
   test.
4. **Containment rule** — `_exposure.py`:
   `weighted_eff_len_contained` (strand-free). Tests for parity,
   strand-flip invariance, constant-A_0.
5. **Per-transcript lengths** — `_exposure_lengths.py`:
   `weighted_transcript_eff_lens`. Includes synthetic nRNA spans
   via existing `_t_exon_intervals`. Test parity.
6. **Estimator EM/output split** — `config.py`, `estimator.py`,
   `pipeline._setup_geometry_and_estimator`. Test split.
7. **Prior rescaling** — `locus_prior.assemble_priors`: compute
   `gdna_prior_count_em = η_g · L_g^A / L_g`. Add to
   `PriorTable`. Wire through `_run_locus_em_partitioned` to
   native EM (no ABI change — same parameter, different value).
   Drop `gdna_prior_count_regional` and
   `_expected_gdna_count_regional_diag`.
8. **Numerator and dead-code deletion** — pipeline
   `_apply_unit_gdna_weights`, `RegionalWeightApplicationStats`,
   `log_weights_for_positions`, `genomic_midpoint` (Python +
   `bam_scanner.cpp`). Update affected tests.
9. **Diagnostic columns** — `quant.feather`, `nrna_quant.feather`,
   `loci.feather` per §7.
10. **Synthetic two-region locus** — `test_regional_exposure_locus.py`.
11. **VCaP gates** and 24-condition regression. Record numbers.

Steps 1–2 unblock everything else and are mutually independent of
the EM/output split. Steps 3–5 depend on 2. Step 6 is independent
and can land in parallel. Steps 7–9 require 1–6. Step 10–11 are
validation.

---

## 11. Surface area

- New files: 4 (`_exposure_lengths.py` + 3 new test modules; the
  rest of new tests go in existing test files).
- Files touched: ~7 (`_regional_exposure.py`, `_exposure.py`,
  `locus_prior.py`, `pipeline.py`, `estimator.py`, `config.py`,
  `bam_scanner.cpp` for `genomic_midpoint` removal).
- Lines deleted: ~250 (per-class signal, numerator path,
  `RegionalWeightApplicationStats`, `genomic_midpoint`
  plumbing, regional diag emitter).
- Lines added: ~450 (shrinkage, prefix sums, two engine
  routines, per-transcript driver, estimator split, tests).
- Native code touched: 1 small change in `bam_scanner.cpp`
  (drop `genomic_midpoint` emission). No ABI change to EM.
- New CLI flags: 0.

---

## 12. Gates

- All existing pytest + golden tests pass in uniform mode (bit
  exact).
- New tests in §9 all pass.
- Synthetic 24-condition harness: no regression in mRNA/nRNA/gDNA
  relative error.
- VCaP rna20m_gdna20m: gDNA contamination estimate within prior
  baseline ± agreed tolerance; intronic spread reduced vs v3
  numerator-era baseline.
- G9 +10% wallclock budget on the existing benchmark sample.

If any gate misses, stop and root-cause before continuing.
