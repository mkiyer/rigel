# New Calibration Plan v2

Date: 2026-05-25
Status: implementation-ready design
Supersedes: `new_cal_plan_v1.md`
Inputs synthesized: `docs/TODO.md` (current selection), `hybrid_capture_failure_research_2026-05-25.md`, `hybrid_capture_exposure_model_2026-05-18.md`, `density_model_impl_plan_v4.md`, `strand_model_impl_plan_v5.md`
Live code audited: [src/rigel/calibration/_orchestrator.py](src/rigel/calibration/_orchestrator.py), [density_model.py](src/rigel/calibration/density_model.py), [density_observation.py](src/rigel/calibration/density_observation.py), [integration.py](src/rigel/calibration/integration.py), [strand_deconv.py](src/rigel/calibration/strand_deconv.py), [region_count_ledger.py](src/rigel/calibration/region_count_ledger.py), [exposure.py](src/rigel/calibration/exposure.py), [prior.py](src/rigel/calibration/prior.py)

---

## 0. North Star

Calibration must answer, **per fine region**, two latent questions:

```text
e_r = P(region r contains any RNA fragments)            "is_expressed"
c_r = P(region r is enriched by the capture assay)      "is_captured"
```

These are continuous because expression is a spectrum and probe efficiency is a spectrum. From them we derive the only two things downstream EM actually needs:

```text
mu_gdna(r)   = expected gDNA fragments in r            (numerator surface)
A_r          = relative gDNA opportunity at r          (denominator surface)
```

Hybrid-capture failure is not a model-shape failure. The 12-channel ledger, FL-aware effective lengths, conjugate Gamma/NB density, and strand deconvolution are correct. The failure is that today's pipeline conflates **off-target background rate** with **local exposure**, then locks the second to the first with an arbitrary precision cap of `1/CV_min^2 = 1/0.05^2 = 400` (see [density_model.py](src/rigel/calibration/density_model.py) line 57). Under capture the off-target rate collapses, the cap stays at 400, and boundary evidence cannot move anything.

v2 fixes this by promoting `e_r`, `c_r`, an off-target rate `rho_off`, and a capture enrichment `gamma_r` to first-class latent variables, then *deriving* `mu_gdna` and `A_r` from them. Strand evidence, when available, accelerates and disambiguates the same fields without changing the contract.

---

## 1. What Has To Go

These are the only things in the current calibration tree that the new design has to remove or repurpose. Everything else stays.

| Live artifact | Status under v2 | Reason |
|---|---|---|
| `DENSITY_PRIOR_MIN_CV = 0.05` and the derived `DENSITY_PRIOR_MAX_PRECISION = 400` | **delete** | Magic precision cap on a rate that is no longer "the local rate". Local exposure has its own concentration in v2. |
| `compute_beta_cap` and `beta_cap`/`cap_applied` on `GammaRatePrior` | **delete** | Same reason. The off-target prior fits as-is from anchor data. |
| `DensityEvidence.relative_exposure` (and the entire "rho_post / rho_ref" pathway as the exposure source) | **delete** | This is the conflation. `A_r` no longer comes from `rho_post / rho_ref`. |
| `RegionExposure.from_density` | **rewrite** | Constructed from the new exposure model, not from density per-region posteriors. |
| `_fit_anchor_priors` family `"ALL"` | **delete** | Background must be fit on *background* regions only, not on the union with exons. Capture data is exactly where "ALL" is wrong. |
| `fuse_density_and_strand` | **rewrite as integration of e_r/c_r** | Today it fuses two posteriors over the *same* observed count `n_total`. In v2 the gDNA posterior is conditioned on `e_r`, `c_r`, `rho_off`, `gamma_r`. |
| `FLAG_PRIOR_DOMINATED`, `n_prior_dominated` | **demote to diagnostic only** | In v2 it should be ~0 for normal use; if it spikes, it indicates the latent model failed to identify on-target regions. |
| `density_max_exposure` CLI/config knob | **delete** | The cap was a workaround for the rho_post/rho_ref pathway exploding. v2 does not need it. |
| `select_rho_ref` `"WEIGHTED_FAMILIES"` path | **delete** | `rho_ref` becomes `rho_off`; we never average exonic densities into background. |
| `density_observation.boundary_left_leff == boundary_right_leff` (both equal `side_leff`) | **keep but rename** | The fact that both sides share the same FL-aware opportunity is correct; the duplicate field is just confusing. |
| `bp_weighted_mean_exposure_over_blocks` in [prior.py](src/rigel/calibration/prior.py) | **keep as fallback only** | A bp-weighted mean of `A_r` is a crude denominator. The proper denominator is the FL-convolved exposure integral. v2 ships the bp mean first and replaces it in a follow-up. See Open Question Q1. |

Nothing in the C++ scanner, signatures, payload, ledger, FL models, strand summary, or strand kappa estimation needs to go.

---

## 2. The Latent Model

### 2.1 Variables

For each fine region `r`:

```text
e_r        in [0, 1]      expression probability
c_r        in [0, 1]      capture probability
gamma_r    >= 1           capture enrichment when captured (rate multiplier)
rho_off                   library off-target gDNA rate (frag / bp)              -- one scalar
kappa_d                   gDNA strand-balance concentration (when stranded)     -- one scalar (already exists)
```

`gamma_r` is per-region because probe efficiency varies. We do **not** need a global on/off ratio as a primary parameter; we report it as a diagnostic (the median or mean of `gamma_r` over `c_r > 0.5` regions).

### 2.2 Per-channel generative model

Let `L_c(r)` = contained gDNA effective length and `L_b(r)` = boundary gDNA opportunity (already produced by `DensityObservation`). Define a per-region local gDNA rate:

```text
lambda_gdna(r) = rho_off * (1 + c_r * (gamma_r - 1))
```

Then the count-generating model the calibration layer fits is:

```text
contained_unspliced_total(r) ~ Poisson( lambda_gdna(r) * L_c(r) + lambda_rna_unspliced(r) * L_c(r) )
boundary_unspliced_total(r)  ~ Poisson( lambda_gdna(r) * L_b(r) + lambda_rna_unspliced(r) * L_b(r) )
contained_spliced_total(r)   ~ Poisson( lambda_rna_spliced(r) * L_c(r) )
strand split of unspliced    ~ BetaBinomial under kappa_d for gDNA, Binomial under SS for RNA   (only when stranded)
```

`lambda_rna_unspliced` and `lambda_rna_spliced` are nuisance rates; we do not need their values, only the indicator `e_r > 0`. They give us two things, though:

1. spliced counts force `e_r -> 1`;
2. strand-deconvolved RNA-lower force `e_r > 0`.

This is the only generative content. Everything below is inference machinery on top of it.

### 2.3 Initialization (data-driven seeds, no magic constants)

```text
SEED 0  --  off-target seed set S_off:
            { r : is_anchor(r)  AND  no spliced mass  AND  (no strand RNA evidence | unstranded) }
            optionally trim top T% by contained density to drop unannotated tx and nRNA hotspots.

SEED 1  --  rho_off = sum(contained[S_off]) / sum(L_c[S_off])

SEED 2  --  e_r seeds:
            spliced_count(r) > 0                                   -> e_r := 1
            strand RNA-lower (R_lower_count(r)) > 0                -> e_r := 1   (only when stranded)
            r in S_off                                              -> e_r := 0
            else                                                    -> e_r := unknown (0.5)

SEED 3  --  c_r seeds:
            r in S_off                                              -> c_r := 0
            r exonic AND no spliced AND no strand RNA              -> c_r candidate (most useful for capture detection)
            else                                                    -> c_r := unknown (0.5)

SEED 4  --  Boundary-flux enrichment test on c_r candidates:
            B_r ~ Poisson(rho_off * L_b(r))    null hypothesis
            p_r = P(X >= B_r | null)
            if p_r < alpha_capture -> seed gamma_r := B_r / max(rho_off * L_b(r), eps)
            else                   -> seed c_r := 0
```

There is **one** tunable threshold here: `alpha_capture` (default 0.01, exposed in `CalibrationConfig`). It controls the false-positive rate for declaring a region as captured. It is not a precision; it is a hypothesis-test alpha and is interpretable.

### 2.4 Iteration (exactly two passes)

Pass 1: refit `rho_off` from `S_off` only, refit `gamma_r` from boundary counts on regions with `c_r > 0.5` and `e_r < 0.5`. No EM yet.

Pass 2: for each region, compute the posterior `p(e_r, c_r | counts, rho_off, gamma_r, kappa_d, SS)` and the bounded posterior over gDNA count `D_r`. The structure is the same exact/approx split that today's `fuse_density_and_strand` already implements, but with the *correct* expected rate `lambda_gdna(r)`.

Two passes are enough because the only feedback loop is `S_off` membership, and that set is dominated by intergenic/intronic regions that the first pass already labels correctly with very high confidence. We deliberately do not iterate further; see Open Question Q3.

### 2.5 Derived outputs for downstream

```text
mu_gdna(r)    = E[ D_r | counts, posteriors ]                      (numerator)
upper_gdna(r) = c-quantile of D_r posterior                        (numerator upper)
A_r           = (1 + c_r * (gamma_r - 1))                          (denominator scalar, dimensionless)
rna_lower(r)  = max(observed_compatible(r) - upper_gdna(r), 0)     (RNA protective lower bound)
```

`A_r` is now an honest per-region scalar describing *how much more gDNA opportunity this region has than the off-target background*, derived from latent variables rather than from a ratio of posterior means.

---

## 3. Module Layout After v2

```
src/rigel/calibration/
  _arrays.py                     # unchanged
  _exposure.py                   # unchanged (FL-aware geometry)
  _fl_sources.py                 # unchanged
  _orchestrator.py               # rewritten body (Phase 4)
  _result.py                     # extended (Phase 4)
  fl.py                          # unchanged
  fractional_evidence.py         # unchanged
  region_count_ledger.py         # unchanged
  regions.py                     # unchanged
  scan_payload.py                # unchanged
  signature.py                   # unchanged
  strand_balance.py              # unchanged
  strand_deconv.py               # unchanged
  strand_summary.py              # unchanged

  density_observation.py         # mild cleanup (Phase 0)
  density_model.py               # gutted: keep rho_off fitting only (Phase 1)
  exposure.py                    # rewritten: A_r from latent c_r, gamma_r (Phase 3)
  integration.py                 # rewritten: posterior over (e_r, c_r, D_r) (Phase 3)
  prior.py                       # mostly unchanged; switch to new A_r (Phase 4)

  # NEW:
  background.py                  # rho_off fit + S_off selection (Phase 1)
  capture_model.py               # boundary enrichment test, gamma_r fit (Phase 2)
  latent_states.py               # e_r, c_r joint posterior + 2-pass refit (Phase 3)
```

No new C++. No changes to the scanner.

---

## 4. PRs (Cleanup → Models → Wiring)

Each PR is independently testable and reviewable. Estimated diff sizes are deliberately small; do not bundle.

### PR 0 — Cleanup of dead state (no behavior change)

1. Delete `DENSITY_PRIOR_MIN_CV`, `DENSITY_PRIOR_MAX_PRECISION`, `compute_beta_cap`. Remove the `beta_cap` / `cap_applied` fields from `GammaRatePrior` and from `to_summary_dict`. Update tests that assert these fields.
2. Delete the `"ALL"` anchor family in `_fit_anchor_priors` and `select_rho_ref`'s `"WEIGHTED_FAMILIES"` path. `rho_ref` is now the `INTRON` prior mean when available, else `INTERGENIC`, else 0 (deterministic-zero path stays).
3. Delete `density_max_exposure` from `calibrate()`, from `CalibrationConfig`, and from the CLI. Delete its only consumer (`RegionExposure.from_density(... max_exposure=)`).
4. Delete `DensityEvidence.relative_exposure`. Anything that reads it (today only `RegionExposure.from_density`) gets the `A_r` field from PR 3.
5. Rename `DensityObservation.boundary_left_leff` and `boundary_right_leff` to a single `boundary_side_leff` field; keep `boundary_leff = 2 * boundary_side_leff` for clarity. The two existing arrays are bit-identical anyway.
6. Move the rho_ref selector logic into a tiny `background.py` placeholder that re-exports it. No new behavior.

**Acceptance**: `pytest tests/ -v` clean except the documented pre-existing strand LLR test. `grep -rE "MIN_CV|MAX_PRECISION|beta_cap|relative_exposure|density_max_exposure"` returns nothing in `src/`.

### PR 1 — `background.py`: off-target rate `rho_off` and seed set `S_off`

```python
@dataclass(frozen=True, slots=True)
class BackgroundModel:
    rho_off: float                       # frag / bp
    n_seed_regions: int                  # |S_off|
    n_fragments: float                   # sum contained on S_off
    eff_length: float                    # sum L_c on S_off
    trim_top_fraction: float             # config knob, default 0.01
    phi: float                           # NB overdispersion of seed-set densities (diagnostic)
    seed_mask: np.ndarray                # bool[R], True where S_off
    fit_status: Literal["ok", "sparse", "fallback_broad", "zero"]
```

Builder: `fit_background_model(observation, strand_counts=None, region_gdna=None, *, trim_top_fraction)`. Inputs:
- contained counts and `L_c` from `DensityObservation`;
- `is_anchor` mask (intergenic + intron-only);
- spliced mass from the ledger (must be zero to remain in `S_off`);
- *optional* strand-deconvolution `rna_lower_count`: if available, drop any anchor with `rna_lower_count > 0`.

Top-K trim removes the top `trim_top_fraction * n_anchors` regions by contained density to defend against unannotated transcripts and nRNA hotspots (matches the user's stated requirement).

`phi` is fit by the existing Pearson-trimmed MoM machinery on `S_off` only. This is purely diagnostic in v2 because we no longer use the prior precision to gate boundary updates.

**Acceptance**:
- Synthetic capture-on/high test: `rho_off` matches the simulator's intergenic density within 10%.
- Non-capture/no-gDNA: `fit_status == "zero"` and downstream gracefully degrades.
- Stranded data: regions with strand-deconvolved `rna_lower_count > 0` are excluded from `S_off`.

### PR 2 — `capture_model.py`: per-region capture enrichment `gamma_r`

```python
@dataclass(frozen=True, slots=True)
class CaptureModel:
    gamma_r: np.ndarray                  # float32[R], >= 1
    c_seed: np.ndarray                   # float32[R], in [0, 1]; capture probability seed
    p_value: np.ndarray                  # float32[R], boundary-Poisson tail p-value
    n_candidates: int
    n_captured_seeded: int
    alpha_capture: float                 # config knob, default 0.01
    fit_status: Literal["ok", "no_signal"]
```

For each region `r` we evaluate the Poisson tail of `B_r` under `rho_off * L_b(r)`:
- `p_value(r) = 1 - PoissonCDF(B_r - 1; rho_off * L_b(r))`
- `c_seed(r) = sigmoid_logp(-log p_value, alpha_capture)` (sigmoid keyed on the chosen alpha; default uses an indicator if you prefer, see Open Question Q4)
- `gamma_r = max(1, B_r / max(rho_off * L_b(r), eps))` for `c_seed > 0.5`, else `1`.

This is the *only* place the boundary enrichment hypothesis test lives. It does not touch contained counts.

When `rho_off == 0` (no-gDNA libraries), `CaptureModel` returns `fit_status == "no_signal"`, `gamma_r == 1`, `c_seed == 0`. This is correct: with no background we cannot detect capture, and we also do not need to.

**Acceptance**:
- Synthetic capture-on/high: `n_captured_seeded` is concentrated on simulator on-target regions; median `gamma_r` of captured regions matches the simulator's enrichment within 30%.
- Capture-off/high: `n_captured_seeded` is at the alpha-rate floor (≤ `alpha_capture * R` false positives).
- Capture-on/no-gDNA: `fit_status == "no_signal"`, all `gamma_r == 1`.

### PR 3 — `latent_states.py`: joint posterior over `(e_r, c_r, D_r)`

This is the only PR that does real Bayesian fusion. Inputs:
- `DensityObservation` (counts and opportunities);
- `BackgroundModel` (rho_off, S_off);
- `CaptureModel` (gamma_r seed, c_seed);
- optional `StrandRegionCounts` + `kappa_d` + `p_r1_sense`;
- optional `RegionGdnaEstimate.rna_lower_count` (strand-derived RNA-lower).

Algorithm (two passes, both vectorized except the per-region exact-posterior inner loop):

```text
Pass 1:
  initialize e_r and c_r from seeds (PR 1 + PR 2 + spliced + strand)
  for each region:
    if observed_compatible_count == 0:
      D_r posterior is delta(0); skip
    if n <= MAX_EXACT_POSTERIOR_N (reuse existing constant = 200):
      enumerate D in [0, n], compute:
        log p(D | e_r, c_r, rho_off, gamma_r) using Poisson(lambda_gdna * L_observed)
        log L_strand(D | k_sense, n, kappa_d, p_r1_sense)         if stranded
        log L_rna(n - D | e_r)                                    weak prior
      normalize, store mean/upper/var
    else:
      bounded Laplace approximation (existing approx path)
  recompute e_r as posterior of (n - D_r > eligibility threshold), capped at 1 when spliced/strand-RNA evidence
  recompute c_r as posterior of (gamma_r > 1) from boundary tail + contained tail

Pass 2:
  refit rho_off using S_off, where S_off membership is re-screened by the updated e_r and c_r posteriors
  refit gamma_r on c_r > 0.5 regions
  repeat the per-region posterior computation once
```

There is no third pass. If two passes are not enough on a fixture, the design has failed and we revisit (see Open Question Q3).

Output:

```python
@dataclass(frozen=True, slots=True)
class LatentStates:
    e_r: np.ndarray                      # float32[R]
    c_r: np.ndarray                      # float32[R]
    gamma_r: np.ndarray                  # float32[R]
    mu_gdna: np.ndarray                  # float32[R]   E[D_r]
    upper_gdna: np.ndarray               # float32[R]
    variance_gdna: np.ndarray            # float32[R]
    rna_lower: np.ndarray                # float32[R]
    A_r: np.ndarray                      # float32[R] = 1 + c_r * (gamma_r - 1)
    flags: np.ndarray                    # uint8[R]
```

Flag bits should reuse what `FusedRegionGdnaEvidence` already has where they still apply. Add `FLAG_EXPRESSED_SEED`, `FLAG_CAPTURED_SEED`, `FLAG_EXPRESSED_POSTERIOR`, `FLAG_CAPTURED_POSTERIOR`.

**Acceptance**:
- Capture-off ablations across (gdna {none, high}) x (ss {0.50, 0.99}): `A_r` ≈ 1 for all regions; no regression on existing tests.
- Capture-on/high/ss=0.99: median `A_r` on simulator on-target exons is in the simulator's enrichment range; gDNA→RNA leak drops at least 5x relative to today.
- Capture-on/no-gDNA/ss=0.99: `mu_gdna ≈ 0` everywhere; we do not invent gDNA from RNA-heavy exons.
- `n_prior_dominated` (now purely diagnostic) is < 5% of regions in the capture suite.

### PR 4 — Rewire `_orchestrator.py` and `exposure.py`

1. `calibrate()` now returns a `CalibrationResult` that carries `BackgroundModel`, `CaptureModel`, `LatentStates`, and the existing `RegionGdnaEstimate` (strand-only, kept for diagnostics).
2. `RegionExposure` becomes a thin reshape of `LatentStates.A_r` plus a flag mask. Its `from_density` constructor is removed; only `RegionExposure.from_latent_states(latent)` remains, plus `uniform()`.
3. `fused_region_gdna` is replaced by `LatentStates` in `CalibrationResult`. Backwards-compatible accessor for one minor release so locus tests do not all churn at once.
4. `_orchestrator.calibrate()` enforces the order: `fl_models → strand_summary → strand_counts → kappa_d → background → capture → latent_states → region_exposure`.
5. The orchestrator must not silently accept `index.region_df is None`; that branch already raises, keep it.

**Acceptance**:
- `summary.json` gains a `latent_states` block and a `capture_model` block.
- All existing tests still pass after the deprecation shim is in place.

### PR 5 — `prior.py`: use new `LatentStates` for numerator and denominator

1. Replace `fused.mean_count` / `fused.upper_count` reads with `latent.mu_gdna` / `latent.upper_gdna`.
2. Replace `bp_weighted_mean_exposure_over_blocks(... region_exposure)` with the same call on the new `region_exposure` (same shape; backed by `A_r` from latent states). Keep the bp-weighted mean as the v2 denominator. The FL-convolved exposure integral is deferred (Open Question Q1).
3. Add `enable_gdna` gating: a locus may only enable gDNA in EM when at least one of its regions has `c_r * gamma_r > 1 + epsilon` *or* the locus contains anchor regions with `mu_gdna > 0`. This drops the current count-based-only enable.

**Acceptance**:
- Synthetic capture suite: on-target locus `gdna_prior_count` recovers (no longer collapses below capture-off levels at high gDNA).
- Capture-off synthetic suite: byte-identical numerical results within tolerance to today.

### PR 6 — Diagnostics, CLI, docs

1. `summary.json` block `capture_detection` with: `mode in {uniform, suspected_capture, capture}`, median `gamma_r`, `n_captured_seeded`, p99 `A_r`.
2. `loci.tsv` and `loci.feather` columns: `n_regions_captured`, `mean_A_r`, `n_regions_expressed_post`.
3. CLI flags: `--alpha-capture`, `--background-trim-top-fraction`. Remove `--gdna-density-confidence` (now part of strand confidence) only if review agrees; otherwise keep as a deprecated alias.
4. Update `docs/calibration/parameters.md` and the manual.

**Acceptance**: the docs and tests reference no removed flags; deprecated aliases warn once.

### PR 7 — Tests and synthetic regression gates

1. Unit tests for `BackgroundModel`, `CaptureModel`, `LatentStates`, each independent of the orchestrator.
2. Add the 8-condition synthetic suite as a `pytest --slow` regression with explicit acceptance thresholds for `mu_gdna` calibration on simulator on-target/off-target masks.
3. Add an oracle-BED smoke test (run only when `tests/scenarios_aligned/capture_oracle.bed` is present) that compares the learned `c_r` against the oracle and reports AUC/precision.

---

## 5. Cross-cutting Implementation Notes

1. **One alpha, one trim fraction.** `alpha_capture` and `background_trim_top_fraction` are the only new tunables. They are statistically interpretable. Do not add precision floors or CV floors anywhere.
2. **No FL re-estimation under capture in v2.** The gDNA FL model is the observed (capture-conditioned) one. Effective lengths use the same observed FL, so numerator and denominator stay consistent. See Open Question Q5.
3. **Determinism.** All RNG (resampling, jitter) must thread the existing `seed` through `calibrate()`. No `np.random.default_rng()` without an explicit seed.
4. **Vectorize first, then approximate.** The exact per-region posterior loop is already O(n^2) capped at `n = 200` and acceptable. Do not pre-optimize.
5. **Flag bits.** Re-use the existing `FusedRegionGdnaEvidence` flag bits where semantics overlap; add new bits only for `e_r` and `c_r` posteriors.
6. **Backwards-compatible `CalibrationResult`.** Keep `fused_region_gdna` as a deprecated property that materializes from `LatentStates` for one minor release.

---

## 6. Open Questions And Unresolved Issues

These are intentionally not deferred to "later". Each one is a real fork in the road that should be answered before, or during, the PR it pertains to.

### Q1. Bp-weighted vs FL-convolved exposure denominator (PR 5)

The current `gdna_eff_len(locus) = unweighted * bp_weighted_mean(A_r)`. A large locus with a small captured island will get a diluted average. The proper denominator is the FL-convolved integral

```text
L_gdna_weighted(locus) = sum_ell h_G(ell) * integral_{valid start s} W(s, ell) ds
```

where `W(s, ell) = mean_{x in [s, s+ell)} A(x)`. Open question: does v2 need this in the first release? The synthetic suite says the bp mean is good enough when capture islands are large (500 kb panel, 1 kb islands), but real exon-scale panels may need the FL integral. **Recommendation**: ship bp mean in PR 5, add FL integral as PR 8 if benchmarks demand it.

### Q2. Local `log W(f)` in fragment scoring

The current pipeline does not add `log W(f)` to per-fragment gDNA log-likelihoods; exposure only enters as a locus denominator. The original hybrid-capture exposure proposal (`hybrid_capture_exposure_model_2026-05-18.md`) argues both must be applied or the EM is inconsistent. Open question: is the locus-level mean enough for synthetic cases? **Recommendation**: instrument first, decide later. Add a `--debug-log-frag-exposure` flag to compute the per-fragment correction without applying it and report the average correction. If the average is < 0.2 nats, skip per-fragment scoring; else add it in PR 9.

### Q3. Convergence of the two-pass scheme

Two passes is a design choice, not a proof. Open question: is there a pathological mixture of expressed-and-captured exons where pass 1 mislabels `S_off` and pass 2 cannot recover? **Recommendation**: write a deliberate adversarial fixture (one locus where every anchor is actually expressed nRNA) and measure. If it breaks, allow up to three passes gated by a convergence test on `|rho_off_new - rho_off_old|`.

### Q4. Hard threshold vs continuous `c_seed`

PR 2 uses a sigmoid on the boundary-Poisson p-value to keep `c_r` continuous. The simpler alternative is a hard indicator `c_seed = 1[p_value < alpha_capture]` and let the joint posterior in PR 3 soften it. Open question: which is more robust under sparse boundary counts? **Recommendation**: hard indicator first; soften only if boundary sparsity causes posterior instability.

### Q5. Capture-conditioned gDNA FL model

The observed gDNA FL mean shifts from ~150 to ~168 bp under capture (length-biased probe overlap). Today's `fl_models.gdna` is trained on observed gDNA fragments and therefore *is* the capture-conditioned distribution. Open question: do we need a "latent" FL for denominator computations? **Recommendation**: no. Keep numerator and denominator both on the same observed FL. The simulator's "pre-capture" mean is not directly observable from a real BAM and is not needed for calibration.

### Q6. Strand evidence and `S_off` selection

PR 1 says: "stranded data drops anchors with `rna_lower_count > 0` from `S_off`". But strand deconvolution runs *after* background fitting in today's orchestrator (look at `_orchestrator.calibrate`). Open question: do we re-order so strand runs first, or do we accept a third "fit background again" sub-pass inside `latent_states.py`? **Recommendation**: re-order. `strand_summary → strand_counts → kappa_d → background → capture → latent_states`. This requires changing exactly two function-call positions in `_orchestrator.py` and one test for ordering.

### Q7. RNA exposure / transcript-level `A_t`

v2 explicitly does **not** add an RNA exposure model. Under capture, uncaptured transcripts will continue to report low TPM. The TODO confirms this is the intended behavior ("Hybrid capture panels are designed intentionally to enrich certain transcripts and deplete others"). Open question: should `quant.feather` carry a `capture_state` column (captured/uncaptured/unknown) so the user can interpret zero-TPM rows? **Recommendation**: yes, add a transcript-level `mean_A_r` and `n_regions_captured` column in PR 6. Do not modify any RNA effective length.

### Q8. Unstranded capture

Without strand evidence, the latent model loses its strongest tool for separating `e_r` from `c_r`. In `latent_states`, we still have boundary flux and spliced counts. Open question: how degraded is the result? **Recommendation**: add unstranded-capture cells to the synthetic suite (gdna_high_ss_0.50_capture_on already exists; add gdna_none_ss_0.50_capture_on). Measure. If gDNA-leak is unacceptable, document that unstranded capture is an unsupported regime in v2 and gate `mode = "capture"` to stranded libraries.

### Q9. Multi-locus regions and the geometry allocator

`prior._allocate_counts_by_geometry` distributes regional counts to overlapping multi-loci by bp. With per-region `mu_gdna` and `A_r` this is still fine for the numerator, but the denominator (`gdna_eff_len`) is computed independently per locus, so a region shared by two multi-loci contributes its `A_r` twice in two denominators. Open question: is double-counting in denominators acceptable? **Recommendation**: yes for v2; the EM does not see the same fragment twice. Document and revisit in PR 8.

### Q10. The "negative control probes deliberately placed in intergenic/intronic space" case

The TODO explicitly calls these out as a confounder for `S_off`: panels deliberately place baseline probes in intergenic regions. The top-K trim handles them statistically but can cause a real off-target probe to be excluded. Open question: do we expose a `--background-exclude-bed` flag? **Recommendation**: not in v2. The top-K trim is enough for the synthetic suite and we cannot test the BED path without ground-truth panels.

---

## 7. What v2 Does *Not* Solve

In keeping with the TODO's directive to stop measuring capture against non-capture:

- v2 does not attempt to recover all-transcript TPM under capture.
- v2 does not attempt to detect capture targets from RNA evidence alone (the `e_r` field can confound `c_r` if `gamma_r ≈ 1`).
- v2 does not implement per-fragment `log W(f)` (deferred to PR 9 contingent on Q2).
- v2 does not add an RNA-side capture exposure.

These are deliberate scope choices. Each has a follow-up PR ticket or an open question above.

---

## 8. Acceptance Summary

The plan is successful when **all** of the following hold on the 8-condition synthetic suite:

| Condition | Acceptance |
|---|---|
| capture-off, gdna=none | byte-identical numerics to today (within float tolerance) |
| capture-off, gdna=high | byte-identical numerics to today |
| capture-on, gdna=none | `mu_gdna ≈ 0` everywhere; no invented gDNA |
| capture-on, gdna=high, ss=0.99 | gDNA→RNA leak drops by ≥5x; on-target loci have `gdna_prior_count` ≥ capture-off levels |
| capture-on, gdna=high, ss=0.50 | gDNA→RNA leak drops by ≥2x (weaker; see Q8) |

And these structural properties:

- `n_prior_dominated` is purely diagnostic and small (< 5% of regions in capture mode).
- No magic constants on probability or precision (no `400`, no `0.05` CV floor).
- All new tunables are statistically interpretable and documented.
- Adding a real probe BED is a strict superset: PR 6's CLI can accept `--capture-bed` to override the seed `c_r` without changing the rest of the pipeline. (We do not ship the flag in v2, but the seam exists.)
