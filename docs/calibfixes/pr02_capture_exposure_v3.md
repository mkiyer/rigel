# PR 2 (v3): Minimal Source-Reliable Capture Exposure

> **Supersedes** `pr02_capture_exposure_v2.md` and `pr02_empirical_bayes_capture_exposure.md`.
> v3 keeps the two useful ideas from v2 and removes the machinery that made the plan hard to ship.

## 0. Verdict On v2

v2 correctly identified two real issues:

1. **The EM sign contract matters.** `prior.py` multiplies `gdna_eff_len` by the regional exposure weight, and the native EM subtracts `log_eff_len`. Therefore capture enrichment must make the regional multiplier smaller, not larger.
2. **Source evidence must be pooled at the unit level.** Per-region strand reliability can miss the Q1 failure because each exon can look RNA-compatible in isolation while the pooled locus has strong gDNA evidence.

Everything else in v2 is too much for PR02:

- no truncated-Normal shrinkage;
- no panel-aware `rho_off` refit;
- no region-level `w_local` jitter weights;
- no prior-mass demotion knob;
- no legacy mode;
- no external panel files;
- no seed thresholds.

The clean PR02 is: fix the sign, pool strand counts by unit, and only lower `A_r` when pooled source evidence clearly supports excess gDNA over background.

## 1. Goal

Replace latent-state-derived `A_r` with a neutral-by-default, unit-pooled, source-reliable exposure multiplier.

PR02 should answer one question:

```text
Does this locus/gene have source-reliable gDNA above off-target expectation?
```

If yes, decrease the gDNA effective length for that unit. If no, leave `A_r = 1`.

## 2. EM Contract

Current EM weight is:

```text
log_weight_{f,g} = log p_score(f | g) + log theta_g - log L_eff,g
```

For gDNA, the scorer emits a position-independent gDNA log-likelihood and the EM applies the per-locus gDNA effective length. Therefore, local capture enrichment must **decrease** the gDNA effective length.

Define `lambda_u >= 1` as local gDNA sampling enrichment for unit `u`. Region exposure is:

```text
A_r = 1 / lambda_{unit(r)}
```

and `prior.py` keeps its existing multiplication:

```text
gdna_eff_len(L) = baseline_gdna_eff_len(L) * mean_A_L
```

So enriched units have `A_r < 1`, giving gDNA a `+log(lambda)` advantage in EM.

Required invariant:

```text
0 < A_r <= 1
A_r = 1 when exposure is not source-identifiable
```

## 3. What Is Broken Today

The current pipeline derives exposure from latent capture states:

```text
mu_gdna = p_captured * captured_mu + (1 - p_captured) * off_target_mu
A_r     = mu_gdna / (rho_off * contained_leff)
```

This is wrong because `p_captured` is a capture/expression stratum probability, not a gDNA source probability. Expressed captured RNA can increase `p_captured`, which then inflates `mu_gdna`, which then changes gDNA EM effective length. That is the feedback loop PR02 must break.

The minimal fixes are:

1. `A_r` must not depend on `p_captured`, `p_expressed_capture`, or `sweep.mu_sweep`.
2. `A_r` must be estimated from source evidence, not latent state evidence.
3. Source evidence must be pooled before strand deconvolution.
4. Non-identifiable units must stay neutral: `A_r = 1`.

## 4. Unit Definition

Use one exposure unit `u` per region:

1. MultiLocus when available;
2. otherwise gene;
3. otherwise no unit.

Fine regions are observations. They are not exposure units.

For PR02, apply a unit-level exposure uniformly to targetable contained regions in the unit. Do not add within-unit local weights yet.

## 5. Target Weight

Use `target_weight`, not `ann`:

```text
T_{r,contained} = 1 for targetable exonic region in a unit
T_{r,boundary}  = 1 when boundary evidence attaches to the same unit
T_{r,c}         = 0 otherwise
```

`target_weight` is not a panel prediction. It is only a deterministic mask saying where capture exposure is allowed to affect gDNA effective length.

## 6. Pooled Source Evidence

Pool raw strand counts over the unit before source deconvolution:

```text
sense_u = sum_{r,c in u} T_{r,c} * sense_{r,c}
anti_u  = sum_{r,c in u} T_{r,c} * anti_{r,c}
n_u     = sense_u + anti_u
L_u     = sum_{r,c in u} T_{r,c} * leff_{r,c}
```

Run the PR01 beta-binomial strand deconvolution once on the pooled unit counts. It should return:

```text
G_u_mean       = posterior mean gDNA count
G_u_lower      = lower credible bound for gDNA count
log_bf_u       = log P(data | mixed source) - log P(data | all RNA)
```

If the current deconvolver does not expose `G_u_lower`, add it. Do not approximate it from per-region outputs.

Background expectation uses the existing background model for PR02:

```text
B_u = rho_off_mean * L_u
```

Do not refit panel-aware `rho_off` in PR02. That can be a follow-up only if benchmark diagnostics show the denominator is the limiting problem.

## 7. Neutral-By-Default Exposure Rule

Use a conservative MAP-style rule, not a truncated-Normal posterior mean:

```text
if mode == "off":
    lambda_u = 1
elif no unit, no target weight, no strand source channel, or L_u <= 0:
    lambda_u = 1
elif log_bf_u <= 0:
    lambda_u = 1
elif G_u_lower <= B_u:
    lambda_u = 1
else:
    lambda_u = max(G_u_mean / B_u, 1)
```

Then:

```text
A_r = 1 / lambda_u          for targetable contained regions in unit u
A_r = 1                     otherwise
```

This rule has the behavior we want:

- ordinary RNA-seq stays neutral;
- RNA-rich capture regions without source evidence stay neutral;
- Q1-like loci can activate because evidence is pooled before deconvolution;
- negative or ambiguous excess cannot lower `A_r`;
- there is no upward truncated-normal bias and no seed threshold.

The only threshold is `log_bf_u > 0`, which is the Bayes-factor MAP boundary where mixed source beats all-RNA. The credible-bound check reuses the existing calibration confidence level; it is not a new capture tuning knob.

## 8. Prior Mass Ownership

Do **not** demote `PriorMassDeconvolution` in PR02.

PR02 changes exposure only. It should not simultaneously retune source-prior strength. The source prior and exposure are different consumers:

- `prior_mass` controls aggregate RNA/gDNA prior mass;
- `A_r` controls gDNA effective length in EM.

If benchmarks later show double counting, make that a separate PR with an isolated switch and benchmark evidence.

## 9. Configuration

Keep config minimal:

```python
capture_exposure_mode: Literal["auto", "off"] = "auto"
```

CLI:

```text
--capture-exposure-mode {auto,off}
```

Do not add:

```text
legacy mode
external panel paths
seed quantiles
prior demotion
w_local / w_max
panel-aware rho_off knobs
```

## 10. Data Structures

```python
CAPTURE_EXPOSURE_DISABLED = np.uint16(0x0001)
CAPTURE_EXPOSURE_NO_UNIT = np.uint16(0x0002)
CAPTURE_EXPOSURE_NO_TARGET_WEIGHT = np.uint16(0x0004)
CAPTURE_EXPOSURE_NO_SOURCE_CHANNEL = np.uint16(0x0008)
CAPTURE_EXPOSURE_NO_BACKGROUND = np.uint16(0x0010)
CAPTURE_EXPOSURE_NOT_SOURCE_SUPPORTED = np.uint16(0x0020)

@dataclass(frozen=True, slots=True)
class CaptureUnitMap:
    unit_id: np.ndarray
    unit_names: tuple[str, ...]
    unit_kind: np.ndarray
    contained_target_weight: np.ndarray
    boundary_left_target_weight: np.ndarray
    boundary_right_target_weight: np.ndarray

@dataclass(frozen=True, slots=True)
class CaptureExposureFit:
    A_r: np.ndarray                    # float32[R], 0 < A_r <= 1
    unit_lambda: np.ndarray            # float64[U]
    unit_background: np.ndarray        # B_u
    unit_gdna_mean: np.ndarray         # G_u_mean
    unit_gdna_lower: np.ndarray        # G_u_lower
    unit_log_bf: np.ndarray            # log_bf_u
    flags: np.ndarray                  # uint16[R]
    unit_flags: np.ndarray             # uint16[U]
```

## 11. Integration Recipe

In `calibration_iteration.py`, stop doing this:

```python
exposure = _safe_exposure(mu_gdna, rho_off * contained_leff)
```

Instead:

1. Build `CaptureUnitMap` after region arrays and gene/locus structure are available.
2. Pool unit strand counts using `target_weight`.
3. Run beta-binomial strand deconvolution on pooled unit counts.
4. Build `CaptureExposureFit` using the neutral-by-default rule.
5. Set `A_r = capture_fit.A_r` and `gamma_r = capture_fit.A_r`.
6. Keep `prior_mass` unchanged.
7. Assert `0 < A_r <= 1` before `prior.py` consumes it.

`prior.py` keeps multiplying baseline gDNA effective length by the bp-weighted mean `A_r`.

## 12. Diagnostics

Emit a TSV derived from fitted units:

```text
inferred_capture_panel.tsv
```

Columns:

```text
unit_id, unit_name, unit_kind, chrom, start, end,
B_u, G_u_mean, G_u_lower, log_bf_u, lambda_u, A_unit, n_regions
```

Optional BED can be derived from the TSV for units with `lambda_u > 1`. It is diagnostic only and never feeds calibration.

Add one sign-contract diagnostic to summaries:

```text
gdna_eff_len_multiplier_summary   # bp-weighted A_r by locus, must be <= 1
```

Do not add per-fragment likelihood-ratio diagnostics in PR02. A focused unit test can prove the EM sign contract more cheaply.

## 13. Tests

Add tests for:

1. `A_r` is always in `(0, 1]`.
2. `mode="off"` gives all-one exposure.
3. no unit gives all-one exposure.
4. high latent capture probability without pooled source support gives all-one exposure.
5. per-region weak evidence but pooled unit evidence activates exposure.
6. `G_u_lower <= B_u` stays neutral even if `G_u_mean > B_u`.
7. strong pooled excess gives `A_r = B_u / G_u_mean` on targetable regions.
8. fine-region splitting does not change unit exposure.
9. `prior.py` multiplying by `A_r < 1` decreases gDNA effective length.
10. native EM sign unit test: smaller gDNA effective length increases gDNA posterior weight for a fixed fragment/locus.
11. `prior_mass` is unchanged by capture exposure fitting.
12. no external panel path/config is accepted by PR02.

## 14. Benchmark Gate

Run the hybrid-capture suite with no probe file passed to `rigel quant`.

Required first-pass gates:

- capture-off: `A_r` is effectively `1` for nearly all regions;
- no-gDNA capture-on: RNA depth alone does not lower `A_r`;
- high-gDNA stranded capture-on: pooled target units get `A_r < 1`;
- Q1/locus bleed: gDNA leakage into RNA is materially reduced;
- unstranded capture-on: no material regression, but full source split remains deferred.

If these fail, diagnose the observed failure before adding machinery. Candidate follow-ups are, in order: background denominator refinement, shrinkage, then local within-unit weighting.

## 15. Review Checklist

- [ ] EM sign contract is documented and tested.
- [ ] `A_r` means effective-length multiplier and satisfies `0 < A_r <= 1`.
- [ ] No latent capture state contributes to `A_r`.
- [ ] Raw counts are pooled before beta-binomial deconvolution.
- [ ] Exposure is neutral unless pooled source evidence beats all-RNA and exceeds background by its lower bound.
- [ ] No truncated-Normal mean, no `w_local`, no panel-aware rho refit, no prior demotion.
- [ ] No external panel file, no seed quantile, no `ann` naming.
