# Fixes: Prior–Likelihood Imbalance under High gDNA, Strand-Specific Inputs (v2)

Status: ready-to-implement
Scope: capture-on and capture-off, strand-specific (SS ≥ 0.9), gDNA-high regime
Replaces: `fixes_prior_noise_v1.md`

## 1. Problem

Two synthetic regression sites motivate this work; both have well-instrumented
truth and well-instrumented intermediate outputs. They expose the same
underlying disease in two different forms.

### 1.1 Capture-ON, SS=0.99, gDNA-high — bleed despite a correct prior

Scenario: `gdna_high_ss_0.99_nrna_none_capture_on`, dominant locus is
`Locus 1` (GENE0002 family), 100,992 EM fragments, true mRNA ≈ 72,963 and
true gDNA ≈ 28,029.

The all-region prior-mass channel assembles 94,869 "local gDNA counts" at
Locus 1, which closely matches the 92–100k regional truth. Yet rigel returns
gDNA = 22,529 (an 8,400-fragment shortfall) and mRNA = 81,350 (an
8,400-fragment overcount). Five thousand or more gDNA fragments leak into
RNA at this single locus.

Diagnostic signature:

| field | value |
|---|---|
| `prior_n_local_gdna` (all-region) | 94,869 |
| `prior_n_local_rna` (all-region) | 79,248 |
| `prior_ess_final` (after cap) | **3,000** |
| `n_em_fragments` | 100,992 |
| `gdna_eff_len` | 33,261 bp |
| `gdna_em_exposure_weight` | **1.104** |
| `capture_enrichment_target` (library-wide) | **1.267** |

The prior is approximately correct. The likelihood is not, and the prior is
too weak to overpower it.

### 1.2 Capture-OFF, SS=0.99, gDNA-high — a single locus collapses

Scenario: `gdna_high_ss_0.99_nrna_none_capture_off`, Locus 7 (GENE0008).
Switching from state-label-derived prior mass to `prior_mass` (all-region)
moves the locus prior `n_local_gdna` from 245 → 2,396 (10×) while truth at
this locus is essentially zero gDNA. The EM responds:

| transcript | truth | baseline | all-region |
|---|---|---|---|
| GENE0008.2 | 340 | 824 | **1,972** |
| GENE0008.3 | 4,119 | 1,047 | **0** |
| GENE0008.4 | 21,854 | 10,534 | 10,411 |

GENE0008.3 collapses to zero. Existing analysis attributed this to
"prior strength interacting with isoform allocation". That description is
incomplete: the prior mass itself is wrong at this locus, and the EM
allocator's failure mode is independent and downstream.

### 1.3 The unifying diagnosis

Both failures are about the **fidelity of the channels that feed the EM**,
not about which mass-semantics philosophy (state-label vs. prior-mass) is
"right". Specifically:

1. The `gdna_em_exposure_weight` channel (per-locus `A_r`) is structurally
   under-estimated in captured regions because its inputs are state-labelled
   `(p_captured, captured_mu)`, and state labels assign captured expressed
   mass to RNA components by construction.
2. The `prior_mass` channel inherits a per-region BetaBinomial deconvolution
   that has near-zero shrinkage in high-`n` regions and no precision floor,
   so a few percent of antisense noise in a long expressed intron is
   amplified into thousands of phantom gDNA counts.
3. The locus-level prior is then ESS-capped to a tiny constant (3,000) so a
   correct prior has no power, and the EM solver's RNA-pool re-allocation
   has no per-isoform regularizer so it can drive a true isoform to zero.

Each pathology is local to one module. None of them require redesigning the
latent-state model or the prior-mass handoff itself.

---

## 2. What today's code does (and where it goes wrong)

### 2.1 The `A_r` exposure channel
File: `src/rigel/calibration/calibration_iteration.py` (around L460–L520).

```python
mu_gdna = p_captured * captured_mu + (1 - p_captured) * off_target_mu
denominator = rho_off_mean * contained_leff
A_r = mu_gdna / denominator                # per-region exposure ratio
```

`A_r` is then averaged bp-weighted into a single `gdna_em_exposure_weight`
per locus (see `_exposure.py::bp_weighted_mean_exposure_over_blocks`) and
multiplied into `gdna_eff_len` before EM.

Failure mode: `captured_mu` comes from `boundary_sweep.mu_sweep`, which is
itself state-conditioned via `sweep`. In a captured expressed region the
state E-step puts almost all mass on `expressed_capture` (no gDNA in that
state) so `captured_mu` is small. `p_captured` is high but `captured_mu` is
tiny, so `mu_gdna ≈ off_target_mu` and `A_r ≈ 1`. The reported
`capture_enrichment_target = 1.267` and per-locus `A_r ≈ 1.10` are direct
consequences.

The downstream effect inside the EM is large. For a captured fragment, the
per-component log-likelihood difference is approximately

$$
\log L_{\text{RNA}} - \log L_{\text{gDNA}}
  \approx \log \tilde L_{\text{gDNA}}
        - \log \tilde L_{\text{RNA}}
  \approx \log \frac{A_r \cdot L_{\text{gDNA,raw}}}{L_{\text{RNA}}}.
$$

With `A_r = 1.1`, `L_gdna_raw ≈ 33 kbp`, and `L_rna ≈ 7 kbp` (Locus 1), the
gap is ~+1.55 nats per fragment in favor of RNA. Times 100k fragments this
swamps a 3,000-count prior.

### 2.2 The `prior_mass` deconvolution
File: `src/rigel/calibration/calibration_iteration.py::build_prior_mass_deconvolution` (L317–L379).

```python
gdna = contained.mean_count + boundary_left.mean_count + boundary_right.mean_count
gdna = clip(gdna, 0, unspliced_total)
rna  = unspliced_total - gdna
```

Each `*_mean_count` is `n_total − R̂` from a BetaBinomial posterior with
`kappa_d ≈ 1e6` (very sharp). For SS=0.99 and an expressed region with
`n_total = n` and `k_antisense = k` ≪ n, the gDNA point estimate behaves
like `Ĝ ≈ 2k`. There is no upper sanity bound, no precision floor, and no
weighting by `p_expressed`.

Failure mode: in a long expressed intron, `k` is small as a fraction but
large in absolute count (FL noise, mis-orientation, residual unmodelled
antisense pre-mRNA). The deconvolution returns thousands of "gDNA" counts
that are actually RNA. This is what drives Locus 7 from 245 → 2,396.

### 2.3 The locus prior ESS cap
File: `src/rigel/calibration/adaptive_prior.py` (L30).

```python
MAX_ESS: float = 3000.0
```

`compute_adaptive_prior` constructs `alpha_gdna_add + alpha_rna_add`, then
caps the sum to `min(locus_unspliced, MAX_ESS)`. The result is that no
matter how accurate `prior_mass` is at a locus and no matter how large
`n_em` is, the prior gets at most 3,000 effective counts. Against 100k
fragments this is a 3% effect — not enough to win against the per-fragment
likelihood imbalance described in §2.1.

### 2.4 The EM RNA re-allocator
File: `src/rigel/native/em_solver.cpp::apply_grouped_prior_update`
(around L709–L757).

```cpp
if (n_rna > EM_LOG_EPSILON) {
    double inv = 1.0 / n_rna;
    for (int i = 0; i < n_components; ++i) {
        if (i == gdna_idx) continue;
        out_counts[i] = R * raw_counts[i] * inv;   // R = n_rna + a_r
    }
}
```

The implementation is isoform-neutral in share-space: the RNA pool gets
`R = n_rna + a_r` mass distributed in proportion to the *raw EM*'s current
RNA-component shares. This is correct as written — it does not actively
push toward dominant isoforms.

Failure mode: the implementation is correct, but there is no per-isoform
floor. A small share-pool isoform whose raw EM share underflows to zero
(through equivalence-class redistribution after gDNA absorbs a few thousand
fragments) stays at zero. GENE0008.3 has truth 4,119 fragments, every one
of which is also compatible with GENE0008.2; once GENE0008.2 grows enough
to "win" all the ambiguous reads, GENE0008.3's raw share is zero and the
prior re-allocation amplifies zero into zero.

---

## 3. Theoretical basis for the fix

### 3.1 Principle: every channel must have its own correctness criterion

The pipeline currently treats `prior_mass`, `A_r`, the ESS cap, and the EM
re-allocator as independent knobs. They are not — they multiply. A correct
prior mass at ESS = 3,000 cannot fix an `A_r` that's 50× too small, and a
correct `A_r` cannot fix an over-confident prior in an expressed intron.

The robust answer is not a new global hyperparameter. It is to give each
channel a self-consistency criterion derived from its own evidence:

- `A_r` must reduce to "raw observed mass relative to off-target
  expectation", with no dependence on a latent state that asserts source
  identity.
- `prior_mass` must clip itself to the maximum mass that the observed
  channel could plausibly contribute under its own noise model.
- The locus prior must scale its ESS to a fraction of the locus's data, so
  it is always meaningful but never dominant.
- The EM RNA re-allocator must include a per-isoform smoothing floor so
  isoforms with shared equivalence classes cannot be driven to identifiable
  zero by a small change in the gDNA prior.

### 3.2 The four fixes are each one-line statements in math

1. **`A_r` is a density ratio of observed mass to off-target expectation.**

   $$
   A_r = \frac{\text{observed compatible unspliced mass}(r)}
                {\rho_{\text{off}} \cdot L_{\text{eff}}(r)},
   \quad \text{clipped to } [1, A_{\max}].
   $$

   This is what `A_r` *means*: how many off-target-equivalents of unspliced
   mass do we actually see here? Replace the state-conditioned mixture
   estimator with this single ratio. It is identically correct in
   off-target regions (where it returns ≈ 1) and recovers the true capture
   enrichment in probe-overlap regions (where unspliced mass is 50–1000×
   the off-target rate).

2. **`prior_mass.gdna_unspliced_mean` is clipped by its own noise model.**

   Under the null "all mass is RNA", the expected antisense count is
   `n_total × (1 − p_r1_sense)`, and the gDNA point estimate is
   `Ĝ ≈ 2 k_antisense`. The plausible-cap is therefore

   $$
   G_{\max}(r) = c \cdot 2 \cdot n_{\text{total}}(r) \cdot
                  (1 - p_{r1\_\text{sense}}) \cdot \sqrt{1 + 1/\kappa_d}
   $$

   with c = 1 (or, conservatively, c = 3 to allow 3σ over-shoot).
   Equivalently: any `Ĝ(r)` that exceeds what the strand channel can
   produce *under the null* must be statistical noise. Shrink it.

3. **Locus prior ESS scales with the locus.**

   $$
   \text{ESS}(\ell) = \min\!\Big(
        \pi \cdot n_{\text{em}}(\ell),\;\;
        \text{precision}_\ell \cdot n_{\text{local}}(\ell)
   \Big)
   $$

   with `π = 0.25` (the prior may contribute at most 25% of the locus's
   effective evidence). This is the standard empirical-Bayes shrinkage
   target. It has no free constants beyond a sensible upper-bound fraction.

4. **EM RNA pool has a minimum per-isoform pseudocount.**

   Replace the dynamic isoform-neutral re-allocator with a one-step Dirichlet
   smoothing: every RNA component receives `α_rna_add / n_rna_components` as
   a floor in addition to share-proportional allocation. This is the
   standard Dirichlet symmetric prior; it is automatic, isoform-neutral by
   construction, and prevents collapse.

   $$
   \alpha_i = R \cdot \frac{\max(\text{raw}_i, \epsilon)}{\sum_j \max(\text{raw}_j, \epsilon)}
   $$

   with `ε = α_rna_add / (n_rna · K_rna)` for K_rna RNA components.

None of these requires a tunable knob in production. Each has a single
well-defined statistical interpretation. Each can be unit-tested against a
synthetic input where the right answer is known by construction.

---

## 4. Implementation plan

The four fixes are independent. They can ship in separate PRs and be
benchmark-tested in isolation. Recommended order: §4.1 → §4.3 → §4.2 → §4.4
(largest expected effect first, lowest risk of new failure modes first).

### 4.1 Fix A — Label-independent `A_r` (highest priority)

#### Files

- `src/rigel/calibration/calibration_iteration.py`
- `src/rigel/calibration/_exposure.py` (no change to API; consumer-side
  only validates the new `A_r` shape and floor)
- `tests/test_calibration_iteration.py`
- `tests/test_compartment_strand_deconv.py` (sanity check, not affected)

#### Change

In `calibration_e_step`, around the existing `A_r = exposure` assignment
(L513), replace the `mu_gdna / (rho_off × L_eff)` computation with a
label-independent density ratio.

```python
# Label-independent A_r: observed unspliced mass per region, normalized by
# the off-target expectation. This is what A_r physically represents and
# does not collapse when the four-state E-step assigns expressed-capture
# mass to RNA components.
observed_unspliced = np.asarray(
    observation.observed_compatible_count, dtype=np.float64
)
off_target_expected = np.maximum(
    float(background.rho_off_mean), 0.0
) * np.maximum(contained_leff, 0.0)

A_r_raw = np.divide(
    observed_unspliced,
    np.maximum(off_target_expected, _EPS),
    out=np.zeros_like(observed_unspliced),
    where=off_target_expected > 0.0,
)
A_r = np.clip(A_r_raw, 1.0, A_R_MAX).astype(np.float32)
```

Add a single module-level constant:

```python
# Maximum per-region exposure ratio. A_r > 1000 is implausible for any
# realistic hyb-capture protocol and is almost certainly an L_eff = 0
# artifact. The clip is a safety rail, not a tuning knob.
A_R_MAX: float = 1000.0
```

Keep `gamma_r` (captured exposure) on the old formula for now — it is used
only for diagnostics and the calibration M-step's enrichment target, both
of which we adjust separately.

#### Why this is robust

The numerator and denominator both come from the same fragment scan; there
is no latent variable in the loop. In off-target regions `observed ≈
rho_off × L_eff` so `A_r → 1`. In probe-overlap captured regions
`observed ≫ rho_off × L_eff` so `A_r` recovers the true enrichment whether
the region is expressed or not.

#### Validation

1. Unit test: synthetic region with observed = 100 × off-target expectation
   must produce `A_r ≈ 100` regardless of `p_states`.
2. Golden test refresh: `tests/golden/` — expect `capture_enrichment_target`
   and per-locus `gdna_em_exposure_weight` to increase substantially in the
   capture-on scenarios. Document the new golden values in the
   refresh PR.
3. Benchmark gate: rerun
   `gdna_high_ss_0.99_nrna_none_capture_on` and assert
   `gdna_em_exposure_weight` p50 > 10 (was 1.05).

#### Rollback

The new constant `A_R_MAX = 1000` can be set to `1.0` to recover today's
behavior in a single edit. Keep the old `_safe_exposure` helper in tree
during the transition for easy comparison.

---

### 4.2 Fix B — Self-consistent `prior_mass` clipping

#### Files

- `src/rigel/calibration/calibration_iteration.py::build_prior_mass_deconvolution`
- `src/rigel/calibration/strand_deconv.py` (no API change; we read
  `kappa_d` and `p_r1_sense` already on the `RegionGdnaChannelEstimate`)
- `tests/test_calibration_prior.py` (new test for the cap)

#### Change

Inside `build_prior_mass_deconvolution`, after the `gdna = contained +
left + right` line and before the `clip(gdna, 0, unspliced_total)` line,
add a null-hypothesis cap:

```python
# Null-hypothesis cap: under H0 "all mass is RNA", the expected antisense
# count is n_total * (1 - p_r1_sense), so the strand deconvolution can
# attribute at most 2 * E[k_anti] worth of gDNA before crossing into
# noise-amplification. Allow a 3-sigma over-shoot to preserve sensitivity
# in truly gDNA-rich regions.
if strand_channels is not None:
    p_sense = float(strand_channels.p_r1_sense)
    kappa = float(strand_channels.kappa_d)
    # NULL_SIGMA = 3.0 — three-sigma noise band, statistical constant.
    sigma_inflation = math.sqrt(1.0 + 1.0 / max(kappa, 1.0))
    gdna_null_cap = (
        2.0
        * unspliced_total
        * max(min(1.0 - p_sense, p_sense), _EPS)
        * (1.0 + NULL_SIGMA * sigma_inflation)
    )
    gdna = np.minimum(gdna, gdna_null_cap)

gdna = np.clip(np.nan_to_num(gdna, nan=0.0, posinf=0.0, neginf=0.0),
               0.0, unspliced_total)
```

with `NULL_SIGMA: float = 3.0` added at module top. The cap is a
statistical inevitability, not a tunable hyperparameter: any `Ĝ` larger
than this is more easily explained by sampling noise on RNA antisense than
by real gDNA.

#### Why this is robust

In capture-OFF SS=0.99 expressed introns (the Locus 7 case),
`p_sense ≈ 0.99` so `min(1-p, p) = 0.01`. For `n_total = 12,000` the cap
is `2 × 12,000 × 0.01 × (1 + 3) = 960`. The current code returns ~2,400.
The cap forces the deconvolution to either explain the antisense excess
with another channel (density) or accept the null. In genuinely gDNA-rich
regions (SS=0.5 capture-off off-target), `min(1-p, p) = 0.5` and the cap
is `2 × n_total × 0.5 × 4 = 4 n_total ≫ unspliced_total`, so the cap is
inactive and the density-based estimate is preserved.

#### Validation

1. Unit test: synthetic region with `p_sense = 0.99`, `n_total = 10,000`,
   `k_antisense = 100` (pure noise) — assert
   `gdna_unspliced_mean ≤ 1,000` (not 2 × k = 200, but the noise cap).
2. Unit test: synthetic region with `p_sense = 0.5`, `n_total = 10,000`,
   `k_antisense = 5,000` (genuine gDNA) — assert cap is inactive.
3. Benchmark gate: rerun `gdna_high_ss_0.99_nrna_none_capture_off` with
   all-region `prior_mass` and assert Locus 7 `prior_n_local_gdna < 500`
   (was 2,396).

#### Rollback

Set `NULL_SIGMA = 1e6` to disable.

---

### 4.3 Fix C — ESS that scales with the locus

#### Files

- `src/rigel/calibration/adaptive_prior.py` (`MAX_ESS`, `compute_adaptive_prior`)
- `src/rigel/calibration/prior.py` (`assemble_priors` passes `n_em_per_locus`
  through)
- `src/rigel/pipeline.py` (collect `n_em_per_locus` once after partition;
  it's already available as `len(partition_tuples[li][...])`)
- `tests/test_adaptive_prior.py`

#### Change

Replace the constant `MAX_ESS = 3000.0` cap with a locus-scaled cap.

```python
# Maximum fraction of a locus's data evidence the prior may contribute.
# Standard empirical-Bayes shrinkage: prior gets at most 25% weight even
# at infinite n. This is the only free parameter in the policy.
PRIOR_DATA_FRACTION: float = 0.25
```

Modify `compute_adaptive_prior` to accept `n_em_per_locus: np.ndarray`
(shape `(n_loci,)`, dtype int64) and to compute, for each locus:

```python
ess_data_cap = PRIOR_DATA_FRACTION * np.maximum(
    n_em_per_locus.astype(np.float64), 0.0
)
# locus_unspliced is the mass we have evidence for; precision (mean of
# prior_mass.precision over the locus's regions, bp-weighted) is the
# evidence quality. The final cap is the more conservative of the two.
ess_evidence_cap = locus_unspliced * locus_precision_bp_weighted
cap_max_loc = np.maximum(np.minimum(ess_data_cap, ess_evidence_cap), 1.0)
```

Then use `cap_max_loc` per-locus where `cap_max` (scalar) was used. The
existing `total_before_cap > cap` path becomes per-locus.

Add a new locus diagnostic in `PriorTable`:

```python
prior_ess_cap_data: np.ndarray      # PRIOR_DATA_FRACTION * n_em
prior_ess_cap_evidence: np.ndarray  # locus_unspliced * precision
prior_ess_final: np.ndarray         # the smaller of the two, after capping
```

#### Why this is robust

The locus prior is now self-scaling. At a small locus with 100 EM
fragments the prior can contribute at most 25 effective counts. At
GENE0002's Locus 1 with 100k EM fragments the prior can contribute up to
25,000 effective counts — 8× today's cap and enough to actually move the
posterior. The 25% fraction is the standard "data dominates but prior is
heard" target; it has no scenario-specific tuning.

#### Validation

1. Unit test: at `n_em = 100`, prior ESS ≤ 25; at `n_em = 100,000`,
   prior ESS ≤ 25,000.
2. Benchmark gate: capture-on SS=0.99 Locus 1 `prior_ess_final` should rise
   from 3,000 to ~10,000–25,000 and gDNA share should rise toward 28%.
3. Regression gate: capture-off SS=0.99 Locus 7 — combined with Fix B, the
   final gDNA share at this locus should be ≤ baseline (state-label
   protected) value.

#### Rollback

Set `PRIOR_DATA_FRACTION = 0.0` and add back the constant `MAX_ESS = 3000`
floor. Single-commit revert.

---

### 4.4 Fix D — Per-isoform RNA pool floor in the EM solver

#### Files

- `src/rigel/native/em_solver.cpp::apply_grouped_prior_update` (L709–L757)
- `src/rigel/native/em_solver.cpp` test bindings, if exposed
- `tests/test_em_impl.py` (new test for isoform-floor)

#### Change

Today's re-allocator is:

```cpp
out_counts[i] = R * raw_counts[i] / n_rna;     // n_rna = sum of RNA raw_counts
```

Replace with a Dirichlet-smoothed re-allocator. Let `K_rna` be the number
of RNA components (n_components − 1 when gDNA is enabled, else
n_components). Define:

```cpp
const int    K_rna   = has_gdna ? (n_components - 1) : n_components;
const double floor_i = (K_rna > 0) ? (a_r / static_cast<double>(K_rna)) : 0.0;
double smoothed_total = 0.0;
for (int i = 0; i < n_components; ++i) {
    if (i == gdna_idx) continue;
    smoothed_total += nonnegative_finite(raw_counts[i]) + floor_i;
}
if (smoothed_total > 0.0) {
    const double inv = R / smoothed_total;
    for (int i = 0; i < n_components; ++i) {
        if (i == gdna_idx) continue;
        out_counts[i] = inv * (nonnegative_finite(raw_counts[i]) + floor_i);
    }
}
```

This is the standard symmetric Dirichlet posterior under a per-isoform
pseudocount of `a_r / K_rna` (the additive RNA prior spread uniformly
across RNA components). It is exactly isoform-neutral, never produces
identifiable zeros for ambiguous-equivalence-class isoforms, and reduces
to today's behavior in the limit `a_r → 0`.

Remove the `carried_state` fallback path; it is no longer needed because
the floor ensures every RNA component gets at least `R / K_rna · (floor_i /
(smoothed_total))` mass.

#### Why this is robust

A symmetric Dirichlet prior is the standard regularizer for categorical
allocation under shared evidence. It has the right limits (zero pseudocount
→ no smoothing; large pseudocount → uniform allocation). It cannot push
toward any specific isoform, by symmetry. It removes the GENE0008.3 → 0
collapse by construction: every isoform receives a baseline share
proportional to `a_r / K_rna`.

#### Validation

1. Unit test: 3-isoform locus, 2 isoforms each carrying 100 raw counts and
   1 isoform with 0 raw counts, `a_r = 30` — assert the zero-isoform
   receives ≥ `30 / 3 / (200 + 30) × R` mass (≈ 4.3% of total).
2. Existing EM tests: must continue to pass with `a_r = 0` (no smoothing
   active).
3. Benchmark gate: capture-off SS=0.99 Locus 7, with Fix B applied, must
   produce GENE0008.3 count > 0.

#### Build step

Requires `pip install --no-build-isolation -e .` after the C++ change.

#### Rollback

Wrap in `if (RIGEL_USE_DIRICHLET_FLOOR) {...}` with a Python-level kill
switch. Keep the old branch under an `#else` for one release.

---

## 5. Cross-cutting diagnostics

These are low-cost additions that make this class of failure visible in
`summary.json` going forward. They are not required for correctness but
are required to prevent regressions from going undetected.

### 5.1 Per-locus likelihood balance diagnostic

In `_batch_locus_em_partitioned` (or a thin wrapper post-pass), accumulate
per-locus `mean(log L_RNA - log L_gDNA)` over all EM fragments that have
both components alive, weighted by posterior responsibility. Emit as
`locus_logL_ratio_rna_minus_gdna` in `locus_stats.feather`. A value > +1
nat means the per-fragment likelihood is systematically biased toward RNA
before any prior — diagnostic of an `A_r` problem.

### 5.2 Per-locus prior-vs-posterior gap

In `assemble_priors`, log:

```python
prior_gdna_share_in  = alpha_gdna_add / (alpha_gdna_add + alpha_rna_add)
em_gdna_share_out    = (from locus_results)
prior_posterior_drift = em_gdna_share_out - prior_gdna_share_in
```

into the loci output. Persistently large negative drift at high-quality
loci is the signature of §1.1.

### 5.3 Per-region prior-mass cap activation

Add a `PRIOR_MASS_NULL_CAPPED` flag (uint16 bit) to the per-region
`prior_mass.flags`, set when Fix B's cap was active. Aggregate the count
into `summary.json::region_calibration.prior_mass.flag_histogram`.

---

## 6. Acceptance criteria for the v2 patch (all four fixes applied)

Measured on
`/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb/`.

| scenario | metric | current | v2 target |
|---|---|---|---|
| `gdna_high_ss_0.99_capture_on` | est gDNA fraction | 0.390 | ≥ 0.46 |
| `gdna_high_ss_0.99_capture_on` | mRNA MARD | 6.75 | ≤ 1.5 |
| `gdna_high_ss_0.99_capture_on` | `gdna_em_exposure_weight` p50 | 1.05 | ≥ 10 |
| `gdna_high_ss_0.99_capture_off` | mRNA MARD | 0.47 | ≤ 0.50 (no regression) |
| `gdna_high_ss_0.99_capture_off` | GENE0008.3 count | 1,047 | ≥ 2,000 |
| `gdna_none_ss_*` (all four) | est gDNA fraction | 0.00 | ≤ 0.001 |
| `gdna_high_ss_0.50_*` | mRNA MARD | unchanged ± 5% | unchanged ± 5% |

The fourth row (capture-off non-regression) is the hardest gate: it
ensures Fixes B + C + D collectively recover the protection that today's
state-label-derived prior accidentally provides at Locus 7, without
re-introducing label-based mass semantics.

---

## 7. Out-of-scope (deferred)

These were considered and intentionally deferred:

- Unstranded capture-on (`ss_0.50_capture_on`) deconvolution. This case
  has no strand signal and the per-region prior mass is wrong upstream of
  any fix in this document. It is the subject of Phase D of the
  gDNA-everywhere design doc and needs a fragment-length-informed
  source-split channel that does not yet exist.
- Renaming the four latent states. Cosmetic; tracked in the
  gDNA-everywhere design doc, Phase A.
- Removal of `capture_gated_prior_mass` snapshot logic from the codebase.
  Safe to delete once §4 lands and benchmarks confirm the all-region path
  is uniformly ≥ the gated path.

---

## 8. PR sequencing

1. **PR-1 (Fix A — `A_r`)**: largest expected delta on the capture-on
   target, smallest blast radius (one function, one constant). Refresh
   goldens.
2. **PR-2 (Fix C — locus-scaled ESS)**: depends on PR-1 in spirit
   (otherwise the larger ESS would amplify the bad `A_r`). Refresh
   goldens.
3. **PR-3 (Fix B — prior_mass null cap)**: removes the Q2 regression
   blocker.
4. **PR-4 (Fix D — Dirichlet floor in EM)**: requires C++ rebuild. Last
   so that EM behavior changes only after the inputs are clean.
5. **PR-5 (Diagnostics §5)**: ship together with PR-4 or as a follow-up.

Each PR is independently revertible. Each is gated by a synthetic-suite
benchmark run.