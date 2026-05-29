# PR 07 v2 — Strand × Density Evidence Fusion (Implementation Plan)

## Status

Implementation plan, ready to execute. Supersedes the earlier draft
`pr07_strand_density_integration.md` and the v1 sketch of this file.

## Conceptual Model

For each region `r` the latent is `D_r` ∈ `[0, T_r]`, the gDNA portion of the observed
unspliced-compatible mass `T_r`. There are two evidence channels, and they are **synergistic**, not
ordered:

- **Density** says how many gDNA fragments the region *ought* to contain based on flanking boundary
  flux and global background. It is naturally a Gamma posterior on the regional gDNA rate, with
  mean `μ_density` and a precision driven by the *discrete number of fragments* contributing to the
  boundary observation.
- **Strand**, conditional on `T_r`, observes how those fragments split between sense (`T_sense`)
  and antisense (`T_anti`). RNA dumps onto one strand at rate `1−κ`. gDNA splits 50/50. So the
  antisense strand is a direct, clean window into roughly half the gDNA population — sharp when
  `κ` is small, flat when `κ → 0.5`.

The fusion is a **closed-form, vectorized Fractional Mixture Allocation** that resolves both
channels in a few lines of NumPy. No grid integration, no exact/approximate split, no per-region
loop. Properties fall out for free:

| Library / region | What the math does |
|---|---|
| Strong strand-specific (`κ → 0`), moderate density prior (`π_r > 0`) | `q_anti → 1`, `q_sense → 0.5·π/(1−0.5·π)`. `D̂_r ≈ T_anti + T_sense·0.5·π/(1−0.5·π)`. Approaches the textbook `2·T_anti` only when `π_r` is large enough to credit a matched gDNA fragment to the sense strand. |
| Strong strand-specific (`κ → 0`), near-zero density prior (`π_r → 0`) | `q_anti → 1`, `q_sense → 0`. `D̂_r → T_anti`. FMA only doubles antisense via the prior — see §1.1. |
| Strong strand-specific, gDNA-dominated region (`π_r → 1`) | `q_anti → 1`, `q_sense → 1`, `D̂_r → T_r`. |
| Unstranded library (`κ → 0.5`) | `q_sense = q_anti = π_r`, `D̂_r = T_r · π_r = μ_density`. Strand becomes neutral; density alone drives the answer. |
| Opposite-strand ambiguous region | Strand structurally absent → `π_r` from density is used directly with `q_sense = q_anti = π_r`. |
| No-gDNA library | Density posterior is concentrated near zero with high precision (large opportunity, many fragments), so `π_r ≈ 0` and `D̂_r ≈ 0`. No `force_zero_gdna_mass` cliff. |

The same fragment is **not double-counted**: density consumes count-per-opportunity geometry,
strand consumes orientation conditional on the same total. Importantly, the boundary-local density
posterior is built from raw boundary counts — never from strand-deconvolved counts.

### §1.1 — Exact `κ → 0, π → 0` limit (math correction)

FMA allocates mass fragment-by-fragment using a fixed prior `π`. It does not perform a global
symmetry imputation. Tracing the limits:

```
q_anti  = 0.5·π / (0.5·π + κ·(1−π))            →  1                 as κ → 0  (any π > 0)
q_sense = 0.5·π / (0.5·π + (1−κ)·(1−π))        →  0.5·π / (1 − 0.5·π)  as κ → 0
q_sense → 0                                      as π → 0
```

So at the joint `κ → 0, π → 0` limit, `D̂_r → T_anti`, **not** `2·T_anti`. The classical
"double the antisense" estimator implicitly assumes RNA mass on the sense strand exactly mirrors
gDNA mass on the antisense strand — a symmetry assumption FMA does not enforce a priori. FMA
recovers the symmetric answer only when the density channel provides `π_r > 0`, which credits
sense-strand reads with a matched gDNA probability.

This is the right behavior: with `π_r → 0` we are explicitly assuming the region holds no gDNA,
and the antisense observation alone is the only honest evidence of gDNA presence. The non-zero
finite low-mean Gamma fallback in the density model (replacing the deterministic-zero branch)
guarantees `π_r > 0` in practice, so the working regime smoothly recovers near-`2·T_anti`
behavior as soon as any non-trivial density prior exists.

Test expectations and goldens encode `D̂_r ≈ T_anti` at the strict `π → 0` limit, and verify
monotone recovery toward `2·T_anti` as the density prior grows.

## The Two Simplifying Ideas

### Idea 1 — Fractional Mixture Allocation as the fusion engine

Let `π_r = clip(μ_density_r / T_r, 0, 1)` be the density-channel prior probability that any given
unspliced fragment in region `r` is gDNA. Let `κ` = library crosstalk (`min(p_r1_sense, 1−p_r1_sense)`
folded so RNA → sense, gDNA → 50/50). For each fragment the posterior gDNA probability is:

```
q_sense_r = 0.5·π_r / (0.5·π_r + (1−κ)·(1−π_r))
q_anti_r  = 0.5·π_r / (0.5·π_r + κ      ·(1−π_r))
```

The fused gDNA mass is a one-line linear allocation of the observed sense/antisense unspliced
counts:

```
D̂_r = T_sense_r · q_sense_r + T_anti_r · q_anti_r
R̂_r = T_r − D̂_r
```

For unstranded data (`κ → 0.5`), `q_sense = q_anti = π_r`, so `D̂_r = T_r · π_r = μ_density_r`. For
perfectly strand-specific data with no RNA (`π_r → 0`, `κ → 0`), `q_anti → 1` and `q_sense → 0`,
so `D̂_r → 2·T_anti_r` — the textbook double-the-antisense estimator. Smooth across both
extremes, no branches.

Posterior variance for `D̂_r` is computable in closed form from the binomial sampling variance of
each strand under those mixture probabilities; that variance gives the fusion's own information,
used downstream by background / exposure refits.

### Idea 2 — Instrument the native accumulator to split per-compartment fragment support

The fusion's `π_r` and its variance both depend on the number of *physical fragments* (not
fractional mass) that contributed to each compartment. Today there is one
`region_unspliced_support` counter per region — compartment-agnostic. PR 07 splits it three ways
in the native accumulator:

```
N_contained_r       (uint32) — unique fragments fully contained in r
N_boundary_left_r   (uint32) — unique fragments crossing r's left boundary
N_boundary_right_r  (uint32) — unique fragments crossing r's right boundary
```

These are the discrete-sample-size denominators that turn fractional density and strand observations
into proper precisions. Total region support `N_r = N_contained + N_left + N_right`.

A fragment that crosses both boundaries of `r` increments both `N_left` and `N_right` once. A
fragment that spans `r` entirely (extends past both boundaries) also increments both. A fragment
fully inside increments only `N_contained`. Each compartment counts a fragment at most once, which
is the natural denominator for that compartment's likelihood.

**Width: `uint32_t`, not `uint64_t`.** Per-region per-compartment fragment counts are small in
practice (a region's contained-fragment count saturates at coverage × region-length / fragment-length,
which is comfortably under 2³² for any realistic locus). Using `uint32_t` halves the memory
footprint of the six support vectors and improves cache behavior in the hot accumulator path. The
`merge_from` reduction stays correct under `uint32_t` because per-thread partials cannot exceed
the global total. A debug-only saturation check is added: if any partial counter approaches
`UINT32_MAX − 2³⁰`, log and abort (no silent wrap). See Phase 1 step 1.

---

## Files Touched

### Net new

- `src/rigel/calibration/strand_evidence.py` — `StrandGdnaEvidence` and helpers (small file).
- `src/rigel/calibration/fusion.py` — fractional-mixture fusion engine (~80 lines including
  variance derivation; replaces `integration.py`).

### Native (C++)

- `src/rigel/native/calibration/accumulator.h` — add three `std::vector<std::uint32_t>` support
  vectors per splice class (contained / left / right). Drop the single aggregate.
- `src/rigel/native/calibration/accumulator.cpp` — replace the post-loop "increment unique
  touched regions" block with three per-compartment increments, each unique-per-fragment. Mirror
  in `merge_from`.
- `src/rigel/native/bam_scanner.cpp` — emit the new keys in the calibration payload dict.
- Rebuild required: `pip install --no-build-isolation -e .`.

### Python plumbing

- `src/rigel/calibration/scan_payload.py` — accept the new uint32 arrays.
- `src/rigel/calibration/_arrays.py` — `PayloadArrays` carries
  `region_contained_unspliced_support_sorted`, `region_boundary_left_unspliced_support_sorted`,
  `region_boundary_right_unspliced_support_sorted` (and spliced counterparts), all `uint32[R]`.
  Drop the aggregate `region_unspliced_support_sorted` (downstream computes `N_r = sum` when needed,
  promoting to `uint64` before the reduction).
- `src/rigel/calibration/region_count_ledger.py` — surface the three views as
  `contained_unspliced_support`, `boundary_left_unspliced_support`,
  `boundary_right_unspliced_support`; drop `unspliced_support` aggregate.

### Rewritten

- `src/rigel/calibration/density_model.py`:
  - `fit_density_evidence` reads `N_left` and `N_right` directly from the ledger and uses them as
    the boundary Negative-Binomial sample sizes.
  - Add `information: float64[R]` derived from the Gamma posterior variance:
    `1 / Var[D_r | density alone]` using the existing NB predictive variance.
  - Delete `_deterministic_zero_evidence` and `PRIOR_FAMILY_DETERMINISTIC_ZERO`. The "no boundary
    evidence" case becomes a finite low-mean Gamma prior with `α = 1.0`,
    `β = α / max(rho_ref, 1e-12)`. Information shrinks to zero when opportunity is also tiny.
  - Delete `select_rho_ref`'s `"ZERO"` branch. `rho_ref_source` stays as a diagnostic string but
    no production code reads it.

- `src/rigel/calibration/strand_deconv.py`:
  - `deconvolve_regions_by_strand`'s `near_unstranded` / ineligible branches emit
    `mean = nan`, `variance = inf`, `information = 0` — never `mean_count = n_total`.
  - Per-region strand information is recomputed using
    `N_r = N_contained + N_left + N_right` as the binomial denominator, not the fractional
    `n_total` mass.

- `src/rigel/calibration/calibration_iteration.py`:
  - Delete `build_region_unspliced_mass()` (the tier ladder), `METHOD_STRAND`, `METHOD_BOUNDARY`,
    `METHOD_BACKGROUND_FALLBACK`, `_METHOD_VALUES`, `_STRAND_RELIABILITY_EPS`,
    `FLAG_M_IMPUTED_FROM_BACKGROUND`, `FLAG_M_CLIPPED_TO_TOTAL`, and the
    `force_zero_gdna_mass` parameter everywhere.
  - Add `build_region_unspliced_mass_from_fused(fused, region_size_bp, n_contained, n_left, n_right)`
    (~20 lines).
  - `RegionUnsplicedMass` schema: drop `method`; add `density_information`,
    `strand_information`, `density_weight`, `strand_weight`, `n_contained`, `n_left`, `n_right`.
  - `calibration_m_step()` refits background from fused mass and per-compartment supports.
  - `estimate_background_density()` pool eligibility: `(density_information + strand_information) > 0`
    weighted by `p_unexpressed · fused_information`.

- `src/rigel/calibration/exposure.py`:
  - Delete `_MIN_EXPOSURE_POOL_P_UNEXPRESSED` and `_METHOD_*` constants.
  - Use `support_total_r = n_contained + n_left + n_right` for `v_obs_r = 1 / max(support_total_r, 1)`.
  - Row weight: `w_r = p_unexpressed_r · fused_information_r · support_total_r / (support_total_r + N0)`.
  - Rename `FLAG_EXPOSURE_IMPUTED_TIER3` → `FLAG_EXPOSURE_NO_INFORMATION`; set when
    `fused_information_r == 0`.

- `src/rigel/calibration/_orchestrator.py`:
  - Delete `_strand_summary_identifiable`, `strand_usable`, the
    `calibration_strand_channels = strand_channels if strand_usable else None` branch, and
    `force_zero_gdna_mass`.
  - Wire: density evidence → strand evidence → `fuse_density_and_strand` → calibration loop.
  - Remove `strand_channels=` argument from `fit_background_model`,
    `build_boundary_local_posterior`, `run_boundary_sweep`. Density and boundary geometry are
    computed from raw counts; strand orientation only appears in the fusion likelihood.

### Deleted

- `src/rigel/calibration/integration.py` — replaced by `fusion.py`. No re-export shim.

### Tests

- New: `tests/test_strand_evidence.py`, `tests/test_fusion.py`,
  `tests/test_fusion_continuity.py`, `tests/test_compartment_support_native.py`.
- Rewritten: every test referencing `METHOD_STRAND`, `METHOD_BOUNDARY`,
  `METHOD_BACKGROUND_FALLBACK`, `force_zero_gdna_mass`, or `_MIN_EXPOSURE_POOL_P_UNEXPRESSED`.
- Goldens refreshed only in Phase 6 after scenario tests pass strict assertions.

---

## Data Contracts

### `StrandGdnaEvidence`

```python
@dataclass(frozen=True, slots=True)
class StrandGdnaEvidence:
    n_sense: np.ndarray              # float64[R]  T_sense_r (fractional mass on sense strand)
    n_anti:  np.ndarray              # float64[R]  T_anti_r  (fractional mass on antisense strand)
    n_total: np.ndarray              # float64[R]  T_r
    support_total: np.ndarray        # uint64[R]   N_contained + N_left + N_right (sum promoted to uint64 from uint32 source vectors)
    eff_support:   np.ndarray        # float64[R]  effective sample size — see §Information gain g
    kappa:   np.ndarray              # float64[R]  per-region crosstalk (folded from p_r1_sense)
    applicable: np.ndarray           # bool[R]     STRUCTURAL only: region overlaps a uniquely-stranded annotated feature
    structural_absent: np.ndarray    # bool[R]     true for TS_NONE / TS_AMBIG  (== ~applicable)
    p_r1_sense: float                # library-global crosstalk parameter
    global_info_scale: float         # in [0, 1] from StrandSummary CI on signed contrast (library-level)
    region_info_gain: np.ndarray     # float64[R]  g_r in [0, 1] — per-region depth-driven gain (see formula below)
```

`κ` is global to the library in v2 (`κ = min(p_r1_sense, 1 − p_r1_sense)`), broadcast per region.
A future enhancement can let `κ` be region-specific from local strand contrast; the schema already
supports it.

#### `applicable` is **strictly structural**

`applicable[r] = True` iff `r` overlaps at least one annotated feature with an unambiguous strand
class (`TS_PLUS` xor `TS_MINUS`). It is derived **entirely from the region table / annotation** —
never from observed counts, depth, expression, or any hypothesis test. This guarantees that no
binary count threshold can re-introduce a behavior cliff. All sample-size and depth-driven
uncertainty flows through the continuous `g` factor below.

#### Information-gain factor `g`

The per-region effective gain that pulls `κ_eff` toward neutral when evidence is weak is the
product of two continuous factors:

```
g_r        = g_library · g_region(eff_support_r)
g_library  = global_info_scale = c² / (c² + (z·se)² + ε),   c = |2·p_r1_sense − 1|,  z = 1.96
g_region   = eff_support_r / (eff_support_r + α0)             # α0 = STRAND_INFO_PSEUDOCOUNT, default 4.0
```

Both factors are in `[0, 1]`, smooth, monotone, and asymptote to 1. No threshold.

#### `eff_support` (effective sample size, not raw fragment count)

Multi-mappers and fractional assignments mean a region can carry many physical fragments but
little independent evidence. To prevent `g_region` from over-trusting depth in those regimes,
`eff_support` uses the *effective* sample size — the sum of per-fragment weights actually
deposited in the region's compartments, **not** raw unique-fragment counts:

```
eff_support_r = T_sense_r + T_anti_r        # sum of fractional contributions for this region
```

`support_total_r` (raw fragment count) is kept for diagnostics and for the binomial variance
formula's `(T/N)²` scaling, where the raw count is the correct denominator. The decoupling is
intentional: variance scaling cares about how many discrete observations occurred; the
information-gain factor cares about how confident those observations are in the per-fragment
assignment.

If `g_region` instead used `support_total_r` (raw uniques), a multi-mapper-saturated region with
`N_r = 15` and `T_r = 0.15` would falsely appear confident. Using `eff_support_r = T_r` removes
this failure mode.

### `DensityEvidence` (extended)

Existing fields stay. Add:

```python
information: np.ndarray              # float64[R]  1 / max(Var[D_r | density alone], var_floor)
applicable: np.ndarray               # bool[R]     contained_leff > min_eff_length
```

Density variance for `D_r` (predicted contained gDNA count) comes from the existing
Negative-Binomial predictive: `Var[D_r] = E[D_r] · (1 + E[D_r] / α_post)`. No new posterior code.

### `FusedRegionGdnaEvidence`

```python
@dataclass(frozen=True, slots=True)
class FusedRegionGdnaEvidence:
    total_mass: np.ndarray              # float64[R]  T_r
    gdna_mass:  np.ndarray              # float64[R]  D̂_r, clipped to [0, T_r]
    rna_mass:   np.ndarray              # float64[R]  T_r − D̂_r (exact)
    variance:   np.ndarray              # float64[R]  Var[D̂_r] from the mixture
    pi_prior:   np.ndarray              # float64[R]  π_r used (diagnostic)
    q_sense:    np.ndarray              # float64[R]  diagnostic
    q_anti:     np.ndarray              # float64[R]  diagnostic
    density_information: np.ndarray     # float64[R]
    strand_information:  np.ndarray     # float64[R]
    density_weight: np.ndarray          # float64[R]  τ_d / (τ_d + τ_s), diagnostic
    strand_weight:  np.ndarray          # float64[R]  τ_s / (τ_d + τ_s), diagnostic
    n_contained: np.ndarray             # uint32[R]
    n_left:      np.ndarray             # uint32[R]
    n_right:     np.ndarray             # uint32[R]
    flags:       np.ndarray             # uint16[R]
```

Flags (replacing tier flags):

```
FUSED_STRAND_STRUCTURAL_ABSENT
FUSED_STRAND_NEAR_UNSTRANDED
FUSED_DENSITY_LOW_OPPORTUNITY
FUSED_LOW_INFORMATION              # density_info + strand_info ≈ 0
FUSED_CLIPPED_TO_TOTAL
FUSED_DEGENERATE                   # T_r == 0
```

### `RegionUnsplicedMass` (revised schema)

Clean cutover: drop `method`, add fusion fields plus the three support arrays. Code consuming
`region_unspliced_mass.method` is rewritten in the same commit.

---

## The Fusion Function (entire body)

`fusion.py::fuse_density_and_strand` is essentially this (with the obvious validation, dtype
contiguity, and clipping):

```python
def fuse_density_and_strand(
    *,
    density_evidence: DensityEvidence,
    strand_evidence: StrandGdnaEvidence,
    region_size_bp: np.ndarray,
    confidence: float = 0.95,
) -> FusedRegionGdnaEvidence:
    T   = np.asarray(strand_evidence.n_total, dtype=np.float64)
    Ts  = np.asarray(strand_evidence.n_sense, dtype=np.float64)
    Ta  = np.asarray(strand_evidence.n_anti,  dtype=np.float64)
    mu_d = np.asarray(density_evidence.mean_unbounded, dtype=np.float64)

    # 1. Density-channel prior, safe against T_r == 0.
    pi = np.zeros_like(T, dtype=np.float64)
    np.divide(mu_d, T, out=pi, where=T > 0.0)
    np.clip(pi, 0.0, 1.0, out=pi)

    # 2. Effective crosstalk. Continuous, structural, no branches.
    kappa_lib  = float(min(strand_evidence.p_r1_sense, 1.0 - strand_evidence.p_r1_sense))
    g_lib      = float(strand_evidence.global_info_scale)         # library-level scalar in [0, 1]
    g_region   = np.asarray(strand_evidence.region_info_gain, dtype=np.float64)  # depth-driven, in [0, 1]
    applicable = strand_evidence.applicable.astype(np.float64)    # structural 0/1, never count-driven
    g          = g_lib * g_region * applicable                    # per-region gain in [0, 1]
    kappa_eff  = 0.5 - (0.5 - kappa_lib) * g                      # κ_lib when g==1; 0.5 when g==0

    # 3. Allocation weights, safe against 0/0 (e.g. π==0 and κ_eff==1).
    denom_s = 0.5 * pi + (1.0 - kappa_eff) * (1.0 - pi)
    denom_a = 0.5 * pi + kappa_eff         * (1.0 - pi)
    q_sense = np.zeros_like(pi, dtype=np.float64)
    q_anti  = np.zeros_like(pi, dtype=np.float64)
    np.divide(0.5 * pi, denom_s, out=q_sense, where=denom_s > 0.0)
    np.divide(0.5 * pi, denom_a, out=q_anti,  where=denom_a > 0.0)

    # 4. Mass allocation by linearity of expectation.
    D_hat = Ts * q_sense + Ta * q_anti
    np.clip(D_hat, 0.0, T, out=D_hat)
    R_hat = T - D_hat

    # 5. Posterior variance — Bernoulli on each strand at its mixture probability,
    #    scaled from physical fragment count (N) to fractional mass (T) by (T/N)².
    N_total = strand_evidence.support_total.astype(np.float64)
    p_sense = np.full_like(T, 0.5)
    np.divide(Ts, T, out=p_sense, where=T > 0.0)
    N_s     = N_total * p_sense
    N_a     = N_total * (1.0 - p_sense)
    scale   = np.zeros_like(T)
    np.divide(T * T, N_total * N_total, out=scale, where=N_total > 0.0)
    var_strand = (N_s * q_sense * (1.0 - q_sense) + N_a * q_anti * (1.0 - q_anti)) * scale

    var_density = np.asarray(density_evidence.variance_unbounded, dtype=np.float64)
    tau_d = 1.0 / np.maximum(var_density, 1.0e-9)
    tau_s = 1.0 / np.maximum(var_strand,  1.0e-9)
    tau_s *= g                                                    # null-out when no strand info
    tau   = tau_d + tau_s
    variance       = 1.0 / np.maximum(tau, 1.0e-9)
    density_weight = tau_d / np.maximum(tau, 1.0e-9)
    strand_weight  = tau_s / np.maximum(tau, 1.0e-9)

    assert np.all(np.isfinite(D_hat))                              # division guards must hold
    # ... fill flags, return dataclass.
```

Every division is guarded by an explicit `where=` predicate so `T_r == 0`, `π == 0 ∧ κ_eff == 1`,
and `N_total == 0` regions emit `0.0` instead of `NaN`/`inf`. The trailing assertion runs in tests
and debug builds.

That's the entire fusion. No grid, no per-region loop, no exact/approximate switch. The
"continuity" guarantees in the v1 sketch are now structural: every expression above is
smooth in `p_r1_sense`, `μ_density`, `T_r`, and `N_total`.

### Why the variance formula is correct

The mixture allocator gives `D̂_r = Σ_i q_i` where each fragment `i` contributes `q_i ∈ {q_sense,
q_anti}`. Treating fragments as independent Bernoulli trials with success probability `q_i`, the
variance of the count is `Σ_i q_i (1−q_i)`. Scaling from physical fragment count `N` to
fractional mass `T` is `(T/N)²` per fragment, which gives the `(T·T) / N²` correction above.
When `N_total` is tiny (sparse evidence), the variance grows and the fusion's information drops
to near zero — the precise behavior we want.

---

## Native Accumulator Change (compartment-split support)

In `accumulator.cpp` the post-fragment block:

```cpp
if (!touched_regions_.empty()) {
    std::sort(touched_regions_.begin(), touched_regions_.end());
    auto last = std::unique(touched_regions_.begin(), touched_regions_.end());
    auto& support = (splice_idx == kSpliceUnspliced)
                        ? payload_.region_unspliced_support
                        : payload_.region_spliced_support;
    for (auto it = touched_regions_.begin(); it != last; ++it) {
        ++support[*it];
    }
}
```

becomes three compartment-keyed sets accumulated inside the per-block loop, each then
`unique`-ified and bumped:

```cpp
// In the per-block loop, alongside add_channel(...):
//   contained?       -> touched_contained_.push_back(reg_id)
//   crosses left?    -> touched_left_.push_back(reg_id)
//   crosses right?   -> touched_right_.push_back(reg_id)

// After the fragment is fully processed:
auto bump = [](std::vector<std::uint32_t>& touched,
               std::vector<std::uint32_t>& support) {
    if (touched.empty()) return;
    std::sort(touched.begin(), touched.end());
    auto last = std::unique(touched.begin(), touched.end());
    for (auto it = touched.begin(); it != last; ++it) ++support[*it];
    touched.clear();
};
if (splice_idx == kSpliceUnspliced) {
    bump(touched_contained_, payload_.region_contained_unspliced_support);
    bump(touched_left_,      payload_.region_boundary_left_unspliced_support);
    bump(touched_right_,     payload_.region_boundary_right_unspliced_support);
} else {
    bump(touched_contained_, payload_.region_contained_spliced_support);
    bump(touched_left_,      payload_.region_boundary_left_spliced_support);
    bump(touched_right_,     payload_.region_boundary_right_spliced_support);
}
```

Aggregate `region_unspliced_support` is dropped; consumers sum the three when they need the total.
`merge_from` updated accordingly.

Build after every native change: `pip install --no-build-isolation -e .`.

---

## Orchestrator Flow (post-PR)

```python
region_arrays   = RegionArrays.from_region_df(...)
payload_arrays  = PayloadArrays.from_payload(...)             # carries 3-way support
ledger          = build_region_count_ledger(payload_arrays)   # exposes 3-way support views
observation     = build_density_observation(region_arrays, ledger, fl_models.gdna)

density_evidence = fit_density_evidence(observation, confidence=..., min_eff_length=...)
strand_counts    = build_strand_region_counts(...)
kappa_d          = estimate_kappa_d(...)
strand_channels  = deconvolve_compartments_by_strand(...)     # still used for kappa_d + logBF_strand
strand_regions   = deconvolve_regions_by_strand(...)          # produces sense/anti folded counts
strand_evidence  = build_strand_gdna_evidence(
    strand_counts=strand_counts,
    strand_summary=strand_summary,
    region_arrays=region_arrays,
    ledger=ledger,                                            # for N_contained/left/right
)

fused = fuse_density_and_strand(
    density_evidence=density_evidence,
    strand_evidence=strand_evidence,
    region_size_bp=region_arrays.region_size_bp,
    confidence=_INTERNAL_GDNA_DENSITY_CI,
)

background = fit_background_model(observation, ...)            # no strand_channels= arg
boundaries = build_boundary_table(region_arrays, ledger, observation.boundary_left_leff)
region_calibration = run_calibration_iteration(
    region_arrays, observation, boundaries, background,
    fused=fused,
    strand_channels=strand_channels,                          # only for logBF_strand and kappa_d
    ledger=ledger,                                            # for support arrays
    max_calibration_passes=...,
    confidence=_INTERNAL_GDNA_DENSITY_CI,
)
```

Deletions in this function alone:
- `_strand_summary_identifiable` call
- `strand_usable` variable
- `force_zero_gdna_mass = density_evidence.rho_ref_source == "ZERO"`
- `calibration_strand_channels = strand_channels if strand_usable else None`
- All references to `density_evidence.rho_ref_source`

---

## Iterative Calibration (post-PR)

`run_calibration_iteration` keeps the E-step / M-step structure. Fusion is computed **once**
before the loop (both channels are functions of the BAM payload, not of `p_unexpressed`).

### E-step

```python
log_tensor = build_state_log_tensor(state_log_prior, logbf_expression, logbf_strand)
p_states   = normalize_state_log_tensor(log_tensor)

prior_mass = build_region_unspliced_mass_from_fused(
    fused,
    region_size_bp=region_arrays.region_size_bp,
    n_contained=ledger.contained_unspliced_support,
    n_left=ledger.boundary_left_unspliced_support,
    n_right=ledger.boundary_right_unspliced_support,
)
region_exposure = estimate_region_exposure(prior_mass, current_density, p_states[:, STATE_UNEXPRESSED])
```

### M-step (`calibration_m_step`)

```python
weighted_gdna = sum(p_unexpressed · seed_mask · fused.gdna_mass)
weighted_eff  = sum(p_unexpressed · seed_mask · contained_leff)
alpha_next    = alpha_floor + weighted_gdna
beta_next     = beta_floor  + weighted_eff
```

The `strand_channels` argument is dropped from `calibration_m_step`. `kappa_d` is carried
separately through `RegionCalibration`.

### Background refit (`estimate_background_density`)

```python
fused_info = region_unspliced_mass.density_information + region_unspliced_mass.strand_information
pool_mask  = (
    (fused_info > 0.0)
    & (region_unspliced_mass.n_contained + region_unspliced_mass.n_left + region_unspliced_mass.n_right >= 1)
    & (region_bp >= 1.0)
)
support_total = (n_contained + n_left + n_right).astype(float64)
weights = fused_info * support_total * p_unexpressed
weights = np.where(pool_mask, weights, 0.0)
# robust geomean + Huber as before
```

`method_histogram` is replaced by `(n_density_used, n_strand_used)`.

### Exposure refit (`estimate_region_exposure`)

```python
support_total = (n_contained + n_left + n_right).astype(float64)
v_obs         = 1.0 / np.maximum(support_total, 1.0)
support_w     = support_total / (support_total + N0)                # N0 = exposure_support_pseudocount, e.g. 4
fused_info    = region_unspliced_mass.density_information + region_unspliced_mass.strand_information
row_weight    = p_unexpressed * fused_info * support_w
row_weight[~np.isfinite(row_weight)] = 0.0
active_pool   = row_weight > 0.0
# τ² fit, shrinkage, ω clipping unchanged below this point.
```

Gone:
- `_MIN_EXPOSURE_POOL_P_UNEXPRESSED`
- `_METHOD_*` constants
- `FLAG_EXPOSURE_IMPUTED_TIER3` (renamed → `FLAG_EXPOSURE_NO_INFORMATION`)

---

## Phases

Each phase ends with `pytest tests/ -q` green and `ruff check src/ tests/` clean. Native changes
require a rebuild before tests.

### Phase 1 — Native compartment-split support

1. Add three `region_*_unspliced_support` + three `region_*_spliced_support` `std::vector<uint32_t>`
   to `accumulator.h`. Drop the `region_unspliced_support` / `region_spliced_support` aggregates.
   Add a debug-only saturation check (e.g. `RIGEL_ASSERT(counter < UINT32_MAX - (1u<<30))`) in the
   bump path.
2. Rewrite the post-fragment increment block in `accumulator.cpp` to bump three compartment-keyed
   `touched_*_` scratch vectors (unique per compartment per fragment) — reuse the existing scratch
   buffer pattern from PR 02 (`bam-scan-pr02-scratch-vectors-2026-05-13`) to keep allocations off
   the hot path. Update `merge_from`.
3. Update `bam_scanner.cpp` to expose the six new arrays in the calibration payload dict.
4. Plumb through `scan_payload.py` → `_arrays.py` → `region_count_ledger.py`, keeping the
   `uint32` dtype end-to-end. Any aggregate-total computation must promote to `uint64` before
   the reduction.
5. New tests: `tests/test_compartment_support_native.py` covers fragment-inside / single-boundary
   crossing / both-boundary spanning regions.
6. **Micro-benchmark gate**: before merging Phase 1, run `scripts/profiling/bam_scan_profile.py`
   on a high-depth sample (any condition from the Armis2 benchmark dir suffices). Confirm:
   (a) wall-time regression ≤ 3 % vs. pre-PR, (b) per-thread RSS regression ≤ 5 %, (c) no new
   cache-miss hotspots in `add_observation`. Record numbers in the PR description.
7. Rebuild: `pip install --no-build-isolation -e .`.

**Phase 1 done when** the calibration payload exposes three per-compartment `uint32` support arrays,
the micro-benchmark gate passes, and existing tests still pass (sum of three equals the old aggregate).

### Phase 2 — Evidence contracts

1. `density_model.py`:
   - Add `information`, `applicable` arrays to `DensityEvidence`.
   - Replace `_deterministic_zero_evidence` with finite low-mean Gamma.
   - Delete `PRIOR_FAMILY_DETERMINISTIC_ZERO`.
2. `strand_deconv.py`:
   - Change `near_unstranded` / ineligible branches to emit `nan / inf / 0`.
3. New file `strand_evidence.py` with `StrandGdnaEvidence` and
   `build_strand_gdna_evidence(strand_counts, strand_summary, ledger, ...)`.
4. New tests: `tests/test_strand_evidence.py`, `tests/test_density_information.py`.

**Phase 2 done when** the two evidence objects can be constructed and queried in isolation, and
unstranded / structurally-absent cases yield zero information.

### Phase 3 — Fusion engine

1. Create `fusion.py` with `FusedRegionGdnaEvidence` and `fuse_density_and_strand` as written
   above. Delete `integration.py` and its tests.
2. New tests: `tests/test_fusion.py`:
   - `π_r > 0, κ → 0`: `D̂_r → T_anti + T_sense · 0.5π/(1−0.5π)`.
   - Strict `π_r → 0, κ → 0`: `D̂_r → T_anti` (see §1.1; **not** `2·T_anti`).
   - `κ → 0.5` (unstranded): `D̂_r = T_r · π_r = μ_density`.
   - Structurally absent (`applicable=False`): `q_sense = q_anti = π_r`, so `D̂_r = T_r · π_r`.
   - `π_r → 1`: `D̂_r → T_r`.
   - `T_r == 0`: `D̂_r == 0` and no `NaN`/`inf` anywhere in the output dataclass (full finiteness assertion).
   - `π_r == 0` and `κ_eff == 1`: same finiteness assertion (the `0/0` guard case).
   - `N_total_r == 0`: `var_strand == 0`, fused variance = density variance.
   - Monotone-recovery sweep: hold `κ = 0`, sweep `π_r` over `[0, 0.5]`, verify `D̂_r` increases
     monotonically from `T_anti` toward the `2·T_anti` neighbourhood.
3. `tests/test_fusion_continuity.py`: sweep `p_r1_sense` over `[0.50, 1.00]` with fixed counts;
   assert `D̂_r`, `variance`, `density_weight`, `strand_weight` are all continuous (max
   first-difference ≤ fixed tolerance) with no jumps.

**Phase 3 done when** `fuse_density_and_strand` produces continuous, finite outputs across all
sweep axes and matches hand-computed posteriors on the edge cases above.

### Phase 4 — Production cutover

1. `_orchestrator.py`: rewrite as in the "Orchestrator Flow" section. Delete
   `_strand_summary_identifiable`, `strand_usable`, `force_zero_gdna_mass`,
   `density_evidence.rho_ref_source` branch.
2. `boundary_model.build_boundary_local_posterior`, `boundary_sweep.run_boundary_sweep`,
   `background_model.fit_background_model`: drop the `strand_channels=` parameter. Internals use
   raw counts only.
3. `calibration_iteration.py`:
   - Delete `build_region_unspliced_mass`, `METHOD_*`, `_METHOD_VALUES`,
     `_STRAND_RELIABILITY_EPS`, `FLAG_M_IMPUTED_FROM_BACKGROUND`, `FLAG_M_CLIPPED_TO_TOTAL`,
     `force_zero_gdna_mass` parameter chain.
   - Add `build_region_unspliced_mass_from_fused`.
   - Update `RegionUnsplicedMass` schema.
   - `calibration_e_step` takes `fused=`, drops `force_zero_gdna_mass=`.
4. Test rewrites — search & migrate:

   ```bash
   grep -rln "METHOD_STRAND\|METHOD_BOUNDARY\|METHOD_BACKGROUND_FALLBACK\|force_zero_gdna" tests/ src/
   ```

**Phase 4 done when** all tests pass and the per-region mass comes only from `fused.gdna_mass`.

### Phase 5 — Background and exposure refit

1. `calibration_m_step`: refit from fused mass + per-compartment support; drop
   `strand_channels` argument.
2. `estimate_background_density`: pool eligibility on
   `(density_information + strand_information) > 0`; replace `method_histogram` with
   `(n_density_used, n_strand_used)`.
3. `estimate_region_exposure`: drop `_MIN_EXPOSURE_POOL_P_UNEXPRESSED` and `_METHOD_*`; use
   continuous row weights with per-compartment support.
4. Tests: `tests/test_exposure.py` updated to assert continuous weighting (no `p_unexpressed`
   threshold cliff).

**Phase 5 done when** no scalar refit branches on the old source method.

### Phase 6 — Diagnostics, goldens, dead-code sweep

1. `_result.py` and summary JSON: emit `density_weight`, `strand_weight`, `density_information`,
   `strand_information`, `fused_information`, per-compartment supports, and the new flag
   histogram. Remove `method_histogram` from outputs.
2. Refresh goldens: `pytest tests/ --update-golden`. Inspect diffs.
3. Dead-code sweep — must be silent:

   ```bash
   grep -rln "METHOD_STRAND\|METHOD_BOUNDARY\|METHOD_BACKGROUND_FALLBACK" src/ tests/
   grep -rln "force_zero_gdna" src/ tests/
   grep -rln "_strand_summary_identifiable\|strand_usable" src/ tests/
   grep -rln "_STRAND_RELIABILITY_EPS" src/ tests/
   grep -rln "_MIN_EXPOSURE_POOL_P_UNEXPRESSED" src/ tests/
   grep -rln "FLAG_M_IMPUTED_FROM_BACKGROUND\|FLAG_M_CLIPPED_TO_TOTAL" src/ tests/
   grep -rln "integration\.py\|from .integration\|from rigel.calibration.integration" src/ tests/
   grep -rln "region_unspliced_support[^_]\|region_spliced_support[^_]" src/ tests/   # aggregate gone
   grep -rln "rho_ref_source.*ZERO\|PRIOR_FAMILY_DETERMINISTIC_ZERO" src/ tests/
   ```

4. Update `docs/newcalib/`: mark v1 of this plan as superseded.

### Phase 7 — Scenario validation + benchmark smoke

Targeted scenarios:
- `anti_intron_ss065_nrna50` (the PR 06 sentinel): expect `gdna_mass ≈ 0` without
  `force_zero_gdna_mass`. `strand_weight` low, `density_weight` high.
- Tiny-positive-gDNA variant: `gdna_mass` small and monotonic in true gDNA fraction.
- Unstranded + boundary-only: density-only fusion, `strand_information ≈ 0`.
- Opposite-strand ambiguous overlap: `strand_evidence.applicable == False`, `D̂_r = T_r · π_r`.
- Strong strand-specific, gDNA-heavy: `D̂_r ≈ T_anti + T_sense · 0.5π/(1−0.5π) ≈ 2·T_anti` once
  `π_r` is non-trivial.
- Multi-mapper-diluted region: high raw `N_r`, low `T_r` → `g_region` stays small,
  `strand_weight` low, `density_weight` dominant.

Then `pytest tests/ -q` and, if Armis2 scratch is mounted,
`python -m scripts.benchmarking run` + `analyze`.

---

## Continuity Guarantees

| Axis | Mechanism (no thresholds anywhere) |
|------|---|
| `p_r1_sense` 0.50 → 1.00 | `g_library = c²/(c²+(z·se)²+ε)`; `κ_eff` continuously interpolates from 0.5 to `κ_lib`. |
| Strand reliability | Per-region information goes through binomial variance on physical support; no `≥ 1e-6` cliff. |
| Density `ρ_ref` 0 → ε | Finite low-mean Gamma fallback; mean and variance smooth. |
| Density opportunity `L_eff` | Posterior naturally widens with low opportunity; no min-opportunity gate. |
| `p_unexpressed` 0 → 1 | Background/exposure pool weights linear in `p_unexpressed`. |
| Effective support `eff_support_r = T_r` 0 → many | `g_region = T_r / (T_r + α0)` smooth, no cliff; sparse / multi-mapper-diluted regions automatically lose strand information. |
| `applicable` flag | Strictly structural (annotation overlap only). Never count- or depth-driven. |

---

## Risk Register

| Risk | Mitigation |
|------|-----------|
| Density + strand double-count fragments | Density uses count-per-opportunity geometry; strand uses orientation conditional on `T_r`. Removing `strand_channels=` from `build_boundary_local_posterior` enforces it structurally. |
| `D̂_r → T_anti` (not `2·T_anti`) at strict `π_r=0` limit surprises consumers | Documented explicitly in §1.1; the finite low-mean Gamma fallback ensures `π_r > 0` in production. Phase 3 test asserts the exact limit; Phase 7 sentinel verifies near-`2·T_anti` recovery with realistic `π_r`. |
| `applicable` quietly becomes data-driven | Schema and orchestrator code carry an `# STRUCTURAL ONLY` comment plus an assertion in `build_strand_gdna_evidence` that `applicable` is computed only from annotation overlap (Phase 2). |
| Multi-mapper-diluted regions over-trust strand | `g_region` uses `eff_support_r = T_r` (fractional mass), not raw fragment count. Phase 7 includes a multi-mapper sentinel. |
| `uint32` support counter overflows on mega-loci | Per-region per-compartment counts are bounded by coverage × region length / fragment length; the largest realistic loci are still 2-3 orders of magnitude below 2³². Debug-only saturation assertion in the bump path (Phase 1.1) aborts on approach to overflow. |
| Variance formula underestimates uncertainty when `T_sense` and `T_anti` are highly imbalanced | The mixture formula treats sense and antisense separately and aggregates exact per-fragment Bernoulli variances; `(T/N)²` scaling keeps fractional-mass and physical-count units consistent. Continuity tests in Phase 3 verify against hand-computed cases. |
| Native accumulator change introduces double-counting | New per-compartment sets are independently `std::unique`-ified; tests in Phase 1 cover inside / single-boundary / double-boundary cases. |
| Phase 1 BAM-scan regression | Phase 1.6 micro-benchmark gate (≤ 3 % wall, ≤ 5 % RSS); reuse PR 02 scratch-vector pattern to keep allocations off the hot path. |
| Tiny positive gDNA gets clipped to zero | Finite low-mean Gamma fallback in density (no deterministic zero) keeps `π_r` strictly positive; tiny-positive sentinel test in Phase 7. |
| Goldens churn obscures regressions | Phase 7 scenario tests are strict-assertion and gate the golden refresh in Phase 6. |
| Downstream consumers of `RegionUnsplicedMass.method` | Grep before Phase 4 cutover; adaptive prior, EM warm start, result writers migrate in the same commit. |

---

## Open Questions (deferred)

1. Region-specific `κ_r` from local strand contrast vs. global `κ_lib`. v2 uses global; the
   `StrandGdnaEvidence.kappa` array is per-region to allow future upgrade.
2. Whether to marginalize over a Beta posterior on `p_r1_sense` (rather than scalar
   `global_info_scale`). Defer.
3. `RegionUnsplicedMass` rename to `RegionGdnaPosterior` — defer to a cosmetic follow-up to limit
   downstream churn here.

---

## Done Criteria

- No code branches on `METHOD_STRAND`, `METHOD_BOUNDARY`, `METHOD_BACKGROUND_FALLBACK`,
  `force_zero_gdna_mass`, `rho_ref_source == "ZERO"`, `_STRAND_RELIABILITY_EPS`, or
  `_MIN_EXPOSURE_POOL_P_UNEXPRESSED`.
- Every per-region `gdna_mass` comes from `fuse_density_and_strand`.
- The fusion engine is closed-form: `D̂_r = T_sense·q_sense + T_anti·q_anti`. No grid integration
  anywhere.
- Native accumulator emits compartment-split (contained / left / right) `uint32` fragment support,
  and the density + strand evidence both consume those discrete sample sizes as their precision
  sources. `g_region` uses *effective* support (`T_r`), not raw counts, to avoid multi-mapper bias.
- `applicable` is strictly structural (annotation overlap only) — grep-asserted in Phase 6.
- The PR 06 sentinel passes without `force_zero_gdna_mass`.
- Dead-code grep checklist (Phase 6) is silent.
- `integration.py` is deleted; `fusion.py` is the production module.
- Full `pytest tests/ -q` green; `ruff check src/ tests/` clean.
