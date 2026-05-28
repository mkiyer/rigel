# PR 2 (v2): Source-Reliable Capture Exposure — Coherent EM Contract

> **Superseded by** `pr02_capture_exposure_v3.md`. Keep this file as the
> detailed critique trail. Do not implement v2 directly.

> **Supersedes** `pr02_empirical_bayes_capture_exposure.md`. v1 had the right
> ambition (kill the latent-state derivation of `A_r`, pool at locus scale,
> emit a continuous unit excess) but contained two foundational errors and
> several statistical artifacts that would cause it to fail on the very
> regression it was designed to fix. v2 keeps v1's organizational stance
> and replaces its mathematical core.

---

## 0. Why v1 Could Not Ship

Two foundational issues, in priority order:

1. **`A_r > 1` decreases the gDNA EM weight, not increases it.**  v1 inherits
   the convention `gdna_eff_len := unweighted * A_r` from
   [`src/rigel/calibration/prior.py`](../../src/rigel/calibration/prior.py).
   In the EM, [`src/rigel/native/em_solver.cpp`](../../src/rigel/native/em_solver.cpp)
   sets `log_eff_len[gdna] = log(gdna_eff_len)` and per-fragment
   `log_weight = log(theta) - log_eff_len`. The scorer
   ([`src/rigel/native/scoring.cpp`](../../src/rigel/native/scoring.cpp))
   emits a position-independent `gdna_log_lik = gdna_fl + gdna_log_sp +
   LOG_HALF + log_nm`. Net: raising `A_r` subtracts `log A_r` from every
   gDNA fragment's score in the locus. A "correctly inferred" `λ_u = 100`
   would make the Q1 bleed *worse*. v1 never wrote down the EM contract,
   so this sign error is invisible at the spec level.

2. **Per-region reliability gating then aggregation collapses on Q1.**  v1
   feeds PR01's per-region `w_strand_r * Ĝ_r` into `G_u`. In an
   expressed-dominated probe-overlap exon, the all-RNA model fits each
   individual region well, so `w_strand_r ≈ 0`, so `G_u ≈ 0`, so
   `λ_u ≈ 1`, so `A_r = 1`. The Q1 failure is exactly the case where
   *pooled* antisense across a locus carries strong gDNA evidence even
   though no single exon does. v1's aggregation order makes the fix
   self-defeating on its target case.

Plus four statistical artifacts (asymmetric clip bias, panel-aware
`ρ_off`, double-counting prior_mass and exposure, intra-unit
heterogeneity), and three engineering gaps (orchestration sequencing,
debug surface, acceptance gates). All addressed below.

---

## 1. Goal

Replace latent-state-derived `A_r` with a coherent, source-reliable,
locus-pooled exposure estimate, with an **explicit EM contract**, a
**shrinkage posterior** on excess exposure, **panel-aware background**,
and **clear ownership of gDNA evidence** between prior mass and exposure.

Works for plain RNA-seq, stranded hybrid capture, and unstranded hybrid
capture without any external panel file, probe BED, or threshold knob.

---

## 2. EM Contract for `A_r` (resolve first)

This is the §1 of the plan. **No code in this PR may be written until this
section is committed to.**

### 2.1 What the EM consumes

The EM per-component log-weight at fragment `f`, component `g` is:

```
log_weight_{f,g} = log p_score(f | g) + log θ_g - log L_eff,g
```

For the gDNA component, `log p_score(f | gDNA)` is position-independent
([`scoring.cpp`](../../src/rigel/native/scoring.cpp) §`score_mm_alignment`),
so all locus geometry enters through `L_eff,gDNA`. To make probe-overlap
fragments more gDNA-favored, the locus must *decrease* `L_eff,gDNA` by
the local enrichment factor.

### 2.2 Contract

Define `λ_u ≥ 1` as the **mean local sampling-opportunity multiplier**
for unit `u` relative to its off-target background. Then:

```
L_eff,gDNA(locus L) := L_eff,gDNA^baseline(L) / λ̄_L
```

where `λ̄_L` is the bp-weighted mean of `λ_{u(r)}` over the regions of
the locus. Equivalently, an `A_r` field on regions:

```
A_r := 1 / (1 + T_r * θ_{u(r)})           with θ_u = λ_u - 1 ≥ 0
```

so that the locus aggregator
[`bp_weighted_mean_exposure_over_blocks`](../../src/rigel/calibration/_exposure.py)
produces `A_L = 1 / λ̄_L`, and
[`prior.py`](../../src/rigel/calibration/prior.py) multiplies:

```
gdna_eff_len(L) := unweighted(L) * A_L = unweighted(L) / λ̄_L
```

Per-fragment `log p(f | gDNA)` then *gains* `log λ̄_L` in the EM
log-weight, which is the desired direction.

### 2.3 Invariants

| Condition | Required `λ_u` | Required `A_r` |
|---|---|---|
| no panel / no source evidence | `1.0` | `1.0` |
| capture-off, ordinary RNA-seq | `1.0` ± shrinkage noise | `≥ 1.0` |
| capture-on, off-target unit | `1.0` | `1.0` |
| capture-on, probe-overlap unit | `≫ 1.0` | `≪ 1.0` for probe regions |

Note the asymmetry vs v1: under v2, `A_r ∈ (0, 1]`. v1's `A_r ≥ 1` is
abandoned. This is a **breaking change** to the `A_r` field semantics
and must be reflected in every downstream consumer and golden file.

### 2.4 Alternative considered, rejected

Adding a per-fragment `+log A(p)` term in `scoring.cpp` is more
expressive (allows true position-dependent enrichment) but requires
per-fragment annotation lookups in the hot scoring loop and breaks the
scoring/EM separation. Deferred to a future PR. v2 stays in EM
denominator-land.

---

## 3. Statistical Model

### 3.1 Unit-level counts pooling (replaces v1 per-region reliability)

For each unit `u` (MultiLocus → gene → none), pool **raw counts** across
all regions/compartments in the unit:

```
n_sense_u      = Σ_{r ∈ u} T_{r,c} * n_sense_{r,c}
n_anti_u       = Σ_{r ∈ u} T_{r,c} * n_anti_{r,c}
n_total_u      = n_sense_u + n_anti_u
contained_leff_u = Σ_{r ∈ u} T_{r,c} * leff_{r,c}
```

Then run the *same* BetaBinomial source deconvolution from
[`strand_deconv.py`](../../src/rigel/calibration/strand_deconv.py) on
the **unit-level counts**, producing:

```
Ĝ_u          = posterior gDNA mass at the unit
w_strand_u   = BetaBinomial log-Bayes-factor for "has gDNA" at the unit
```

This is the central correction over v1. The strand deconvolver already
handles the precision/noise floor; we use its unit-level posterior, not
a weighted sum of region-level posteriors.

### 3.2 Background expectation (panel-aware)

The off-target gDNA rate `ρ_off` is calibrated library-wide today, but in
a capture-on library it is *depleted* (probes pulled reads toward
targets), which inflates every `λ_u = G_u / B_u`. v2 estimates `ρ_off`
from **clearly-off-panel evidence only**:

```
ρ_off_panel-aware =
    (Σ over intronic+intergenic compartments with low source-reliability mass) Ĝ_{r,c}
  / (Σ over the same set) leff_{r,c}
```

When the eligible set has insufficient mass (`Σ leff < MIN_OFFTARGET_LEFF`),
fall back to the existing library-wide `ρ_off` but raise the shrinkage
prior width (§3.4) by 2× to flag the loss of identifiability.

Then:

```
B_u = ρ_off_panel-aware * contained_leff_u
```

### 3.3 Raw excess

```
λ_u_raw = Ĝ_u / B_u
ε_u     = λ_u_raw - 1            # signed excess; can be negative
```

`ε_u` is the data summary; it is *not* the final exposure.

### 3.4 Shrinkage posterior on `θ_u` (replaces `max(λ−1, 0)`)

v1's asymmetric clip is biased upward under noise; on a no-capture sample
every unit drifts to `A_r > 1` from sampling variance alone. v2 places a
half-Normal-on-`θ_u` prior with empirical-Bayes scale:

```
θ_u | data  ~  half-Normal(0, σ_θ)  *  N(ε_u, se_u^2)        truncated to [0, ∞)

σ_θ = MAD-based across-unit scale of {ε_u} after centering at 0
se_u^2 = Var(λ_u_raw | counts)  from delta-method on the BetaBinomial posterior
```

Posterior mean `θ̂_u` shrinks toward `0` when `|ε_u| ≲ se_u`, and
recovers `ε_u` when `ε_u ≫ se_u`. **This replaces v1's clip-at-zero
entirely.** Closed-form posterior is a truncated Normal; no MCMC needed.

Reliability weighting is applied on top:

```
θ̂_u_final = sigmoid(w_strand_u - W_TAU) * θ̂_u
```

with `W_TAU` a fixed-by-derivation BF threshold (log-BF = 0 → 0.5 weight).
This keeps low-reliability units near neutral without a hard cutoff.

### 3.5 Region exposure (with within-unit weighting)

A simple uniform application across the unit understates the bleed at
probe positions. v2 applies `θ̂_u_final` per-compartment with a
data-driven local weight:

```
A_r = 1 / (1 + T_{r,c} * θ̂_u_final * w_local_{r,c})

w_local_{r,c} = clamp( (Ĝ_{r,c} / leff_{r,c}) / (Ĝ_u / Σ leff_u),  [0, w_max] )
```

`w_local_{r,c}` is the region's local gDNA density divided by the unit
mean; it concentrates the exposure boost on regions where source-reliable
gDNA actually accumulates. `w_max = 10` caps a single noisy compartment.
When `Ĝ_{r,c}` is missing or unreliable, fall back to `w_local = 1`
(unit-mean behavior, matching v1's intent).

Boundary compartments get their own `T_{r,boundary}` and so contribute
to *both* fitting and region exposure — the v1 asymmetry (boundaries
fit but don't lift) is removed.

---

## 4. Evidence Ownership: prior_mass vs A_r

Source-reliable `Ĝ` currently flows into both
`PriorMassDeconvolution.gdna_unspliced_mean` (locus prior `α_gdna_add`)
and would now also flow into `λ_u → A_r`. That is the **same evidence
entering the posterior twice**.

**v2 decision: `A_r` owns the source-reliable evidence; `prior_mass`
demotes to a soft seed.**

Concretely, in [`calibration_iteration.py::build_prior_mass_deconvolution`](../../src/rigel/calibration/calibration_iteration.py):

```python
# v1: gdna_unspliced_mean = contained.mean_count + left + right   # full strength
# v2: gdna_unspliced_mean = PRIOR_DEMOTION * (contained + left + right)
#                            where PRIOR_DEMOTION ∈ (0, 1], default 0.25
```

The justification is information-theoretic: per-fragment likelihood (via
`A_r`) is strictly stronger evidence than a Dirichlet pseudo-count on
share, so `A_r` should dominate and the prior is reduced to a
regularizer that keeps the EM from collapsing in zero-coverage loci.

Acceptance test: with `PRIOR_DEMOTION = 1.0` (full double-count) vs
`0.25` (demoted), the latter must match or beat the former on all
benchmark gates in §10.

---

## 5. Configuration

```python
# config.py
capture_exposure_mode: Literal["auto", "off", "legacy"] = "auto"
capture_exposure_prior_demotion: float = 0.25      # §4
capture_exposure_w_max: float = 10.0               # §3.5
capture_exposure_min_offtarget_leff: float = 1e6   # §3.2 fallback threshold
```

```text
--capture-exposure-mode {auto,off,legacy}
--capture-exposure-debug                  # emit per-unit TSV without changing EM
```

`legacy` retains the v1-style `μ_gdna / (ρ_off × leff)` estimator for
one release for A/B benchmarking. Removed in PR03.

No seed quantiles, enrichment thresholds, or panel files. No `ann`
naming. v1's targetability primitives `T_{r,c}` are retained:

```
T_{r,contained} = 1   on targetable exonic unit regions
T_{r,boundary}  ∈ [0,1] when boundary mass attaches to the unit
T_{r,c}         = 0   for unrelated intronic/intergenic compartments
```

---

## 6. Orchestration

v1 underspecified when `fit_capture_exposure()` is called relative to
MultiLocus construction. Concrete sequencing for v2:

1. `_orchestrator.py` builds `region_arrays`, loci, boundaries → produces
   `CaptureUnitMap` keyed by **gene** (always available from index) and
   optionally **MultiLocus** if the partition is already built.
2. `build_prior_mass_deconvolution()` runs (unchanged structure).
3. `fit_capture_exposure()` runs against unit-level pooled counts (§3.1),
   consuming `prior_mass` only for diagnostics, not as evidence.
4. `prior_mass` is then *demoted* per §4 before being handed to
   `prior.py`.
5. `prior.py` builds `gdna_eff_len` from `A_r` per the §2.2 contract.

If MultiLocus is not yet available at step 1, gene is the fallback unit;
this is documented in `CaptureUnitMap.unit_kind`. Calibration must not
depend on MultiLocus being constructed first.

---

## 7. Data Structures

```python
CAPTURE_EXPOSURE_DISABLED              = uint16(0x0001)
CAPTURE_EXPOSURE_NO_UNIT               = uint16(0x0002)
CAPTURE_EXPOSURE_NO_TARGET_WEIGHT      = uint16(0x0004)
CAPTURE_EXPOSURE_NO_SOURCE_RELIABILITY = uint16(0x0008)
CAPTURE_EXPOSURE_NO_BACKGROUND         = uint16(0x0010)
CAPTURE_EXPOSURE_OFFTARGET_FALLBACK    = uint16(0x0020)   # §3.2
CAPTURE_EXPOSURE_SHRUNK_TO_ZERO        = uint16(0x0040)   # §3.4
CAPTURE_EXPOSURE_LOCAL_WEIGHT_CLIPPED  = uint16(0x0080)   # §3.5

@dataclass(frozen=True, slots=True)
class CaptureUnitMap:
    unit_id: np.ndarray                       # int32[R], -1 when no unit
    unit_names: tuple[str, ...]
    unit_kind: np.ndarray                     # uint8[U]: 0=gene, 1=multilocus
    contained_target_weight: np.ndarray       # float32[R]
    boundary_left_target_weight: np.ndarray
    boundary_right_target_weight: np.ndarray

@dataclass(frozen=True, slots=True)
class CaptureExposureFit:
    # Region-level (consumed by prior.py via A_r contract §2.2)
    A_r: np.ndarray                           # float32[R], in (0, 1]
    region_theta_local: np.ndarray            # float32[R], T * θ̂_u * w_local
    region_local_weight: np.ndarray           # float32[R], w_local_{r,c}
    flags: np.ndarray                         # uint16[R]

    # Unit-level (consumed by diagnostics + panel BED)
    unit_eps_raw: np.ndarray                  # float64[U], ε_u (signed)
    unit_eps_se: np.ndarray                   # float64[U], se_u
    unit_theta_posterior: np.ndarray          # float64[U], θ̂_u
    unit_theta_final: np.ndarray              # float64[U], after w_strand gate
    unit_lambda: np.ndarray                   # float64[U], 1 + θ̂_u_final
    unit_background: np.ndarray               # float64[U], B_u
    unit_gdna_source: np.ndarray              # float64[U], Ĝ_u
    unit_strand_logbf: np.ndarray             # float64[U], w_strand_u
    unit_flags: np.ndarray                    # uint16[U]

    # Hyperparameters used (audit)
    rho_off_used: float
    rho_off_panel_aware: bool
    sigma_theta: float                        # §3.4 empirical-Bayes scale
    prior_demotion: float                     # §4
```

---

## 8. Diagnostics

In addition to the v1 panel BED/TSV (kept as a *derived* report, not a
calibration input):

```
inferred_capture_panel.tsv columns:
  unit_id, unit_name, unit_kind, chrom, start, end, strand,
  B_u, G_u, w_strand_u, eps_u, eps_se_u, theta_posterior_u,
  theta_final_u, lambda_u, n_regions, locus_id
```

**New mandatory per-locus diagnostic**: emit the per-locus mean
`log p(f|RNA) − log p(f|gDNA)` over the locus's EM fragments, before
and after `A_r` is applied. This is the only diagnostic that catches a
sign-flip regression of §2 directly. Add to `RegionCalibration`
summary JSON as `loci_em_lr_pre_Ar` and `loci_em_lr_post_Ar`.

---

## 9. Files Touched

| File | Change |
|---|---|
| `src/rigel/config.py` | Add 4 fields from §5. |
| `src/rigel/cli.py` | Add `--capture-exposure-mode {auto,off,legacy}`, `--capture-exposure-debug`. |
| `src/rigel/calibration/capture_exposure.py` | **New.** Unit pooling, BetaBinomial-at-unit, panel-aware ρ_off, shrinkage posterior, region exposure with `w_local`. |
| `src/rigel/calibration/calibration_iteration.py` | Replace latent-state ratio; demote `prior_mass` per §4. Keep `legacy` branch one release. |
| `src/rigel/calibration/_orchestrator.py` | Build `CaptureUnitMap` from gene/MultiLocus; pass to `fit_capture_exposure()`. |
| `src/rigel/calibration/strand_deconv.py` | Expose a `deconvolve(counts)` callable that accepts pre-pooled unit counts (§3.1). |
| `src/rigel/calibration/_result.py` | Add `CaptureExposureFit`, `loci_em_lr_pre_Ar`, `loci_em_lr_post_Ar`. |
| `src/rigel/calibration/prior.py` | Update `assemble_priors` doc + invariant assertion for `A_r ∈ (0, 1]`. Math is unchanged (multiply). |
| `src/rigel/calibration/_exposure.py` | No math change. Verify `bp_weighted_mean_exposure_over_blocks` is correct under `A_r ≤ 1`. |
| `src/rigel/native/em_solver.cpp` | No code change. Add comment at L1847 documenting the §2.2 contract. |
| `src/rigel/pipeline.py` | Emit panel BED/TSV + per-locus likelihood-ratio JSON. |
| `scripts/benchmarking/runner.py` | No probe file path. |
| `tests/test_capture_exposure.py` | **New** (see §11). |
| `tests/test_per_locus_gdna_mass.py` | Update for `A_r ∈ (0, 1]`. |
| `tests/golden/` | Regenerate after benchmark gate (§10) passes. |

---

## 10. Benchmark Gates (quantitative, gate-the-merge)

| condition | metric | gate |
|---|---|---|
| no-gDNA, capture-off | `p95(A_r)` over all regions | `≥ 0.95` (i.e., A_r ≈ 1) |
| no-gDNA, capture-on | `p95(A_r)` over off-target regions | `≥ 0.95` |
| no-gDNA, capture-on | `p5(A_r)` over probe-overlap regions | `≤ 0.30` (boost present even without gDNA truth — calibration knows panel from strand-pool noise floor only) |
| high-gDNA, capture-off, SS=0.99 | `p95(A_r)` | `≥ 0.90` |
| high-gDNA, capture-on, SS=0.99 | median(`λ_u`) over targeted units | `≥ 10` |
| high-gDNA, capture-on, SS=0.99 | gDNA bleed into RNA at locus 1 | `≤ 500 fragments` (was ~5000) |
| capture-off regression | GENE0008.3 quant | within 5% of v1-baseline |
| all conditions | per-locus `loci_em_lr_post_Ar < loci_em_lr_pre_Ar` for elevated-`A_r` loci | sign check: §2 contract holds |

The last row is the **non-negotiable** sign check that detects §0 issue
1 regressions automatically.

---

## 11. Tests

```text
tests/test_capture_exposure.py

  Contract:
    test_A_r_in_unit_interval                 # A_r ∈ (0, 1] always
    test_neutral_when_disabled                # mode='off' ⇒ A_r ≡ 1
    test_neutral_when_no_unit                 # CaptureUnitMap.unit_id = -1
    test_neutral_when_no_source_reliability   # w_strand_u below floor
    test_legacy_mode_matches_v1               # mode='legacy' bit-equal to old

  Unit pooling (§3.1):
    test_pool_counts_then_deconvolve_vs_per_region_collapse
        # the Q1 regression: per-exon ≈ all-RNA, unit ≈ has-gDNA
    test_fine_region_split_invariance         # splitting a region
                                              # doesn't change unit Ĝ_u

  Panel-aware background (§3.2):
    test_rho_off_uses_offpanel_only
    test_rho_off_fallback_raises_sigma        # flag set, σ_θ doubled

  Shrinkage posterior (§3.4):
    test_no_capture_sample_shrinks_to_zero    # no upward bias
    test_strong_signal_recovers_eps           # ε_u ≫ se_u ⇒ θ̂_u ≈ ε_u
    test_symmetric_noise_no_bias              # E[θ̂_u | true=0] ≈ 0

  Region exposure (§3.5):
    test_w_local_concentrates_on_probe_regions
    test_w_local_clipped_at_w_max
    test_boundaries_get_lifted                # not just contained

  Evidence ownership (§4):
    test_prior_mass_is_demoted
    test_demotion_does_not_regress_no_capture # capture-off invariant

  EM consequence:
    test_A_r_below_1_increases_per_fragment_gdna_score   # direct §2.2 check
    test_per_locus_lr_diagnostic_emitted

tests/test_per_locus_gdna_mass.py
    update existing assertions for A_r ∈ (0, 1]
```

---

## 12. Migration / Rollout

1. Land v2 behind `capture_exposure_mode='legacy'` default for one release;
   add the §10 gates to CI as warnings only.
2. Run hyb-capture 500kb benchmark suite with `auto` and verify all §10
   gates pass.
3. Flip default to `auto`; downgrade `legacy` to opt-in.
4. PR03: remove `legacy`, remove demoted `prior_mass` codepath if
   benchmarks confirm A_r alone is sufficient.

---

## 13. Review Checklist (additions over v1)

- [ ] §2 EM contract written down before any code is changed.
- [ ] `A_r ∈ (0, 1]` invariant asserted in `prior.py` and tested.
- [ ] BetaBinomial deconvolution called on **pooled unit counts**, not
      per-region results.
- [ ] `ρ_off` computed from off-panel evidence with documented fallback.
- [ ] `θ_u` produced by truncated-Normal shrinkage posterior, not by
      `max(·, 0)`.
- [ ] Per-region `w_local` applied so probe-overlap regions get the boost.
- [ ] Boundary compartments contribute symmetrically to fit and exposure.
- [ ] `prior_mass` gDNA contribution demoted; double-counting eliminated.
- [ ] `loci_em_lr_pre_Ar` / `loci_em_lr_post_Ar` emitted and gated in §10.
- [ ] `legacy` mode retained for one-release A/B.
- [ ] No external panel file, no seed quantile, no `ann` naming.
