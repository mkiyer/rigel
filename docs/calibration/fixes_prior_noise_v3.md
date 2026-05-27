# Fixes for Prior Noise and gDNA Handoff, v3

Date: 2026-05-27

This document reviews and supersedes `fixes_prior_noise_v2.md`. The main
conclusion is that v2 correctly identified the failure surfaces, but two of the
proposed fixes need to be made more conservative before implementation:

- `A_r` must not be computed from raw observed unspliced mass. Raw observed mass
  is dominated by RNA in expressed captured regions.
- `prior_mass` must not receive a hard null cap. A hard cap removes real gDNA in
  high-gDNA strand-specific scenarios.

The implementation-ready version below preserves the clean semantic model:
latent states stratify expression and capture exposure; they are not source
labels. All four latent states can contain gDNA. Source mass handed to EM comes
only from `PriorMassDeconvolution`.

## 1. What was verified

I verified the v2 hypotheses against the current code paths and the existing
hybrid-capture suite at:

```text
/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb
```

I added and ran:

```text
scripts/debug/review_prior_noise_v2.py
```

The script writes:

```text
/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb/diagnostics/prior_noise_v2_review/
```

Key outputs:

- `prior_noise_v2_review.md`
- `exposure_counterfactuals.tsv`
- `prior_null_cap_counterfactual.tsv`
- `ess_and_floor_counterfactual.tsv`
- `gene0008_quant_focus.tsv`

Important diagnostic readouts:

| check | result | decision |
|---|---:|---|
| Raw observed exposure in no-gDNA SS=0.99 capture-on probe exons | observed/off-target ratio p50 = 161k, p95 = 4.36M | Reject raw observed `A_r` |
| Prior-gDNA exposure in high-gDNA SS=0.99 capture-on probe exons | prior-gDNA/off-target ratio p50 = 972, p95 = 11.7k; true ratio p50 = 977, p95 = 13.3k | Use prior-gDNA source signal |
| Prior-gDNA exposure in no-gDNA SS=0.99 capture-on probe exons | prior-gDNA/off-target ratio p50 = 1.74k, but only 159 total false gDNA counts | Use Gamma/count shrinkage, not raw ratios |
| v2 hard null cap in high-gDNA SS=0.99 capture-on | would reduce prior gDNA from 98.7k to 14.7k | Reject hard cap |
| v2 hard null cap in high-gDNA SS=0.99 capture-off | binds regions containing 97.8k true gDNA | Reject hard cap |
| Locus 1, high-gDNA SS=0.99 capture-on | `n_em_fragments` = 100,992; current ESS cap = 3,000; data-fraction cap at 25% = 25,248 | Accept locus-scaled ESS |
| Locus 7, high-gDNA SS=0.99 capture-off | all-region prior makes GENE0008.3 exactly zero | Add EM RNA floor, but treat it as a stability guard |

A useful exposure shrinkage sanity check was also run outside the persisted TSVs.
For SS=0.99 capture-on probe exons, a per-region Gamma-Poisson exposure estimate
with prior count `tau = 100` gives:

| condition | source prior-gDNA sum | off-target expected sum | shrunken `A_r` p50 | shrunken `A_r` p95 |
|---|---:|---:|---:|---:|
| no gDNA | 158.9 | 0.0419 | 1.017 | 1.227 |
| high gDNA | 92,315 | 101.3 | 26.9 | 69.3 |

That is the shape we want: near 1 in no-gDNA capture-on data, but strongly
elevated when the source split itself carries high absolute gDNA mass.

## 2. v2 verdict

| v2 item | verdict | v3 replacement |
|---|---|---|
| Fix A: label-independent `A_r` from observed unspliced mass | Reject as written | Source-specific, Gamma-shrunk `A_r` from `prior_mass.gdna_unspliced_mean` |
| Fix B: null-hypothesis hard cap on `prior_mass` | Reject as written | Strand-specific excess-over-RNA-noise shrinkage, applied before prior assembly |
| Fix C: ESS scales with locus data | Accept with edits | Data-fraction cap first; precision cap only when the precision channel is meaningful |
| Fix D: per-isoform RNA floor | Accept as stability guard | Symmetric Dirichlet floor in grouped-prior redistribution, not a full isoform-resolution fix |

Recommended PR order changes from v2:

1. Source-specific Gamma-shrunk `A_r`.
2. RNA-strand-noise shrinkage for `prior_mass`.
3. Locus-scaled ESS.
4. Symmetric RNA floor in native EM.
5. Cross-cutting diagnostics.

This order avoids amplifying false strand deconvolution mass before the false-mass
shrinkage is in place.

## 3. Fix A v3: source-specific Gamma-shrunk `A_r`

### Problem

Today `A_r` is derived from the latent-state mixture:

```python
mu_gdna = p_captured * captured_mu + (1.0 - p_captured) * off_target_mu
A_r = mu_gdna / (rho_off * contained_leff)
```

This collapses in expressed capture regions because the latent state can assign
the region to `expressed_capture` while source mass still contains gDNA. The
high-gDNA SS=0.99 capture-on run shows exactly that: `prior_mass` contains about
98.7k gDNA, but current probe-exon `A_r` has p50 only 1.51 and p95 2.10.

v2 proposed replacing this with raw observed unspliced mass. That is unsafe:
RNA-only probe exons in no-gDNA SS=0.99 capture-on have raw observed/off-target
ratios in the hundreds of thousands to millions.

### Change

Compute `A_r` from the explicit source split, not from latent states and not from
raw observed mass:

```python
A_R_MAX: float = 1000.0
A_R_SOURCE_PRIOR_COUNT: float = 100.0

def _source_specific_gdna_exposure(
    gdna_source_mass: np.ndarray,
    off_target_expected: np.ndarray,
) -> np.ndarray:
    numerator = np.maximum(gdna_source_mass, 0.0) + A_R_SOURCE_PRIOR_COUNT
    denominator = np.maximum(off_target_expected, 0.0) + A_R_SOURCE_PRIOR_COUNT
    exposure = numerator / np.maximum(denominator, _EPS)
    return np.clip(exposure, 1.0, A_R_MAX).astype(np.float32)
```

In `calibration_e_step`, build `prior_mass` before `A_r` and then use:

```python
off_target_expected = (
    max(float(background.rho_off_mean), 0.0)
    * np.maximum(contained_leff, 0.0)
)
A_r = _source_specific_gdna_exposure(
    prior_mass.gdna_unspliced_mean,
    off_target_expected,
)
```

`gamma_r` can remain on its current captured-mixture formula for the first PR. It
is an enrichment diagnostic and M-step signal, not the EM gDNA effective-length
weight. If later diagnostics show `capture_enrichment_target` is still misleading,
revise it in the diagnostics PR rather than coupling it to this behavioral fix.

### Why this is theoretically cleaner

`A_r` is a source-specific exposure multiplier for gDNA effective length. The
only compatible numerator is source-specific gDNA mass. The Gamma-Poisson prior
count is a statistical shrinkage prior centered at `A_r = 1`; it prevents tiny
off-target expectations from converting tens or hundreds of false source counts
into enormous exposure multipliers.

With `tau = 100`, the no-gDNA SS=0.99 capture-on probe-exon false signal remains
near 1, while the high-gDNA SS=0.99 capture-on probe-exon exposure rises above
the target p50 > 10 gate.

### Files

- `src/rigel/calibration/calibration_iteration.py`
- `tests/test_calibration_iteration.py`
- `tests/test_per_locus_gdna_mass.py` or a new focused exposure test module

### Tests

Add focused tests for `_source_specific_gdna_exposure`:

1. `gdna_source_mass = 0`, `off_target_expected > 0` gives `A_r = 1`.
2. `gdna_source_mass = off_target_expected` gives `A_r = 1` after clipping.
3. `gdna_source_mass = 1000`, `off_target_expected = 10`, `tau = 100` gives
   `(1100 / 110) = 10`.
4. `gdna_source_mass = 5`, `off_target_expected = 1e-3`, `tau = 100` gives about
   1.05, not thousands.

Benchmark gates after PR-1:

- `gdna_high_ss_0.99_nrna_none_capture_on`: probe-exon or locus
  `gdna_em_exposure_weight` p50 > 10.
- Both `gdna_none_*_capture_on` scenarios: estimated gDNA fraction remains <=
  0.001 or within current no-gDNA tolerance.
- `gdna_high_ss_0.50_nrna_none_capture_on`: no claimed fix yet; no material
  regression beyond noise.

## 4. Fix B v3: strand excess-over-RNA-noise shrinkage

### Problem

The all-region prior-mass path correctly stops treating latent states as source
labels, but it exposes a second issue: strand deconvolution can attribute RNA
strand leakage to gDNA in expressed strand-specific regions. The capture-off
GENE0008/Locus 7 regression is the visible failure.

v2 proposed a hard null cap:

```python
gdna = min(gdna, 2 * n_total * q_error * noise_factor)
```

The diagnostics show this is not safe. In high-gDNA SS=0.99 capture-on, it would
reduce total prior gDNA from 98.7k to 14.7k, even though the true gDNA is 100k.

### Change

Shrink only the amount of gDNA that cannot be explained as excess over the RNA
strand-error null. This belongs in `strand_deconv.py`, before
`PriorMassDeconvolution` consumes the per-region estimates.

For a region with `n` compatible observations, `k_sense` sense-strand
observations, and strand model `p = p_r1_sense`, define the RNA error channel:

```python
if p >= 0.5:
    k_error = n - k_sense
    q_error = 1.0 - p
else:
    k_error = k_sense
    q_error = p
```

Under an all-RNA null:

```python
E_error = n * q_error
sd_error = sqrt(n * q_error * (1.0 - q_error))
error_upper = E_error + RNA_NOISE_Z * sd_error
```

Only counts above `error_upper` support gDNA. Since gDNA contributes 0.5 to the
error channel and RNA contributes `q_error`, the maximum gDNA mass supported by
the excess is:

```python
denom = max(0.5 - q_error, _EPS)
gdna_excess_supported = max(0.0, (k_error - error_upper) / denom)
gdna_mean = min(gdna_mean_from_existing_deconv, gdna_excess_supported)
```

Use:

```python
RNA_NOISE_Z: float = 3.0
STRAND_FLAG_RNA_NOISE_SHRUNK: int = 1 << 5
```

Apply the shrinkage only when the region is eligible for an RNA-strand null:

- transcript-strand information is available;
- `0 <= q_error < 0.5 - eps`;
- `n > 0`.

Do not apply it to ineligible/intergenic regions, where there is no RNA null to
test against. Do not apply it in near-unstranded mode where `q_error` is close to
0.5; the formula becomes uninformative and the existing density/fallback channel
must carry the case.

### Why this preserves true gDNA

The rule is not a global cap. It asks whether the observed minor-strand channel
exceeds what RNA strand leakage can plausibly explain. A truly mixed RNA+gDNA
region has a much larger minor-strand count and therefore survives the shrinkage.
A pure-RNA strand-specific region with ordinary antisense leakage is shrunk.

This directly targets the Locus 7 regression without deleting real high-gDNA
capture-on mass.

### Files

- `src/rigel/calibration/strand_deconv.py`
- `src/rigel/calibration/calibration_iteration.py` only if flags need wiring
- `tests/test_strand_deconv.py` or a new focused test module
- `tests/test_calibration_iteration.py`

### Tests

Add deterministic tests around the helper:

1. Pure RNA, strand-specific: `p = 0.99`, `n = 10000`, `k_error` at the expected
   null plus 1 sigma. Assert supported gDNA is 0 or near 0.
2. True mixed gDNA: `p = 0.99`, `n = 10000`, `k_error = 5000`. Assert supported
   gDNA is near `n`, minus the 3-sigma allowance.
3. Near unstranded: `p = 0.50`. Assert the shrinkage is inactive.
4. Ineligible region: assert existing density/fallback behavior is unchanged.

Benchmark gates after PR-2:

- `gdna_high_ss_0.99_nrna_none_capture_off`: Locus 7 all-region
  `prior_n_local_gdna` should drop substantially from 2,396 without reducing
  true high-gDNA off-target regions globally.
- `gdna_high_ss_0.99_nrna_none_capture_on`: total prior gDNA must remain close
  to 100k, not collapse toward the v2 hard-cap 14.7k counterfactual.

## 5. Fix C v3: locus-scaled adaptive-prior ESS

### Problem

The current global `MAX_ESS = 3000` makes the prior nearly irrelevant in large
loci. In high-gDNA SS=0.99 capture-on Locus 1, there are 100,992 EM fragments;
the prior can contribute only 3,000 effective counts. That is too small to move a
strong but biased likelihood surface.

### Change

Use a data-fraction cap:

```python
PRIOR_DATA_FRACTION: float = 0.25
```

In `assemble_priors`, derive the per-locus data count before partitioning:

```python
n_em_per_locus = np.asarray(
    [len(locus.unit_indices) for locus in multi_loci],
    dtype=np.float64,
)
```

Pass this into `compute_adaptive_prior`. Per locus:

```python
ess_cap_data = PRIOR_DATA_FRACTION * np.maximum(n_em_per_locus, 0.0)
ess_cap_unspliced = np.maximum(locus_unspliced, 0.0)
ess_cap_final = np.minimum(ess_cap_data, ess_cap_unspliced)
```

Use `ess_cap_final` where the scalar `MAX_ESS` cap is currently applied.

Do not add the v2 precision cap in this PR. Current `PriorMassDeconvolution`
precision is meaningful for strand-informative regions but is zero for important
density/fallback paths. A naive precision cap would accidentally disable priors
in the unstranded cases. Add precision scaling only after the precision contract
is made consistent across source-split methods.

### Diagnostics to add

Extend `PriorTable` and locus output with:

```python
prior_ess_cap_data
prior_ess_cap_unspliced
prior_ess_final
```

Keep `prior_n_local_gdna`, `prior_n_local_rna`, `alpha_gdna_add`, and
`alpha_rna_add` as separate diagnostics. They answer different questions.

### Files

- `src/rigel/calibration/adaptive_prior.py`
- `src/rigel/calibration/prior.py`
- `src/rigel/pipeline.py` only if new locus diagnostics need output wiring
- `tests/test_adaptive_prior.py`
- `tests/test_pipeline_wiring.py` if output schema changes

### Tests

1. `n_em = 100`, `locus_unspliced` large: final ESS <= 25.
2. `n_em = 100000`, `locus_unspliced` large: final ESS <= 25000.
3. `locus_unspliced = 500`, `n_em = 100000`: final ESS <= 500.
4. Existing small-locus tests remain numerically stable or are updated with
   explicit expected caps.

Benchmark gates after PR-3:

- High-gDNA SS=0.99 capture-on Locus 1 `prior_ess_final` rises from 3,000 toward
  the 25k data cap, unless limited by available unspliced mass.
- High-gDNA SS=0.99 capture-off does not regress after Fix B is applied.
- No-gDNA scenarios retain near-zero estimated gDNA.

## 6. Fix D v3: symmetric RNA floor in native EM

### Problem

`apply_grouped_prior_update` redistributes aggregate RNA prior according to raw
RNA counts. If a transcript has zero raw count in the update, it can receive zero
after the grouped-prior redistribution even when it is still a plausible isoform.
The GENE0008.3 collapse in the all-region snapshot is the concrete symptom.

### Change

Use a symmetric Dirichlet floor when redistributing the aggregate RNA prior.

Let:

- `R` be the target aggregate RNA mass after the grouped update.
- `a_r` be the additive RNA prior mass.
- `K_rna` be the number of RNA components excluding gDNA.
- `floor_i = a_r / K_rna`.

Then:

```cpp
double smoothed_total = 0.0;
for (int i = 0; i < n_components; ++i) {
    if (i == gdna_idx) continue;
    smoothed_total += nonnegative_finite(raw_counts[i]) + floor_i;
}

if (smoothed_total > 0.0) {
    const double scale = R / smoothed_total;
    for (int i = 0; i < n_components; ++i) {
        if (i == gdna_idx) continue;
        out_counts[i] = scale * (nonnegative_finite(raw_counts[i]) + floor_i);
    }
}
```

Keep the existing carried-state fallback only for the pathological case where
`a_r <= 0` and all raw RNA counts are zero. In the normal adaptive-prior path,
the floor should make the fallback unnecessary.

### Interpretation

This does not solve isoform identifiability. It prevents exact zeros introduced
by the grouped-prior redistribution. The guaranteed zero-isoform share is only:

```text
(a_r / K_rna) / (sum_raw_rna + a_r) * R
```

The diagnostics showed Locus 7 all-region `alpha_rna_add / K_rna` is about 348,
so this is a stability guard, not a reason to expect GENE0008.3 to recover to its
truth by itself.

### Files

- `src/rigel/native/em_solver.cpp`
- `tests/test_batch_em_impl.py`

### Tests

1. Add or replace the current aggregate-prior test with a three-RNA-component
   case where one component has zero raw RNA count and `a_r > 0`. Assert it
   receives positive mass matching the Dirichlet formula.
2. Assert `a_r = 0` preserves current behavior.
3. Assert the grouped RNA and gDNA totals still sum to the intended aggregate
   totals after smoothing.

### Build

After the C++ edit:

```bash
conda activate rigel && pip install --no-build-isolation -e .
```

Then run the focused native tests before broader suite tests.

## 7. Cross-cutting diagnostics

These should land either with the relevant PRs or immediately after Fix D.

### 7.1 Exposure diagnostics

Persist per-locus summaries of the new source-specific `A_r`:

- `gdna_em_exposure_weight`
- `gdna_source_exposure_p50`
- `gdna_source_exposure_p95`
- `gdna_source_exposure_max`

The capture-on high-gDNA target should visibly separate from no-gDNA capture-on.

### 7.2 Prior source-quality diagnostics

Add per-region and per-locus flag histograms for:

- strand deconvolution ineligible;
- RNA-noise shrinkage active;
- Gamma-shrunk exposure clipped at `A_R_MAX`;
- prior ESS data cap active;
- prior ESS unspliced cap active.

### 7.3 Prior-posterior gap

For each locus, emit:

```text
prior_gdna_share_in = alpha_gdna_add / (alpha_gdna_add + alpha_rna_add)
posterior_gdna_share_out = gdna_count / (gdna_count + mrna_count + nrna_count)
prior_posterior_gdna_drift = posterior_gdna_share_out - prior_gdna_share_in
```

Large negative drift after Fix A+C means likelihood/exposure still overwhelms a
high-quality source prior. Large positive drift in no-gDNA scenarios means the
prior or likelihood is overcalling gDNA.

## 8. Acceptance criteria

All criteria should be evaluated on the eight-condition hybrid-capture suite.
The goal is not only to improve the target scenario; it is to preserve all four
data modes: strand-specific yes/no, capture yes/no.

### Required non-regression gates

| scenario group | gate |
|---|---|
| All `gdna_none_*` scenarios | estimated gDNA fraction <= 0.001 or no worse than current tolerance |
| High-gDNA SS=0.99 capture-off | mRNA MARD no worse than baseline/current by more than 5% |
| High-gDNA SS=0.50 capture-off | mRNA MARD no worse than baseline/current by more than 5% |
| High-gDNA SS=0.50 capture-on | no material regression; do not expect a fix from these PRs |

### Required improvement gates

| scenario | gate |
|---|---|
| High-gDNA SS=0.99 capture-on | `gdna_em_exposure_weight` p50 > 10 |
| High-gDNA SS=0.99 capture-on | estimated gDNA fraction moves substantially toward 0.5 |
| High-gDNA SS=0.99 capture-on | mRNA MARD improves relative to current capture-gated/all-region snapshots |
| High-gDNA SS=0.99 capture-off Locus 7 | all-region source prior no longer drives GENE0008.3 to an exact zero |

The old v2 target `GENE0008.3 count >= 2000` should not be used as a hard gate.
The EM floor cannot guarantee that much recovery. The correct gate is no exact
zero plus no scenario-level mRNA MARD regression; if the isoform remains badly
misallocated, investigate isoform-specific likelihood/equivalence-class structure
as a separate issue.

## 9. Validation commands

For Python-only PRs:

```bash
conda activate rigel && ruff check src/rigel/calibration tests scripts/debug/review_prior_noise_v2.py
conda activate rigel && pytest tests/test_calibration_iteration.py tests/test_strand_deconv.py tests/test_adaptive_prior.py -v
```

For the native EM PR:

```bash
conda activate rigel && pip install --no-build-isolation -e .
conda activate rigel && pytest tests/test_batch_em_impl.py -v
```

For full proof after all four fixes:

```bash
conda activate rigel && python scripts/debug/review_prior_noise_v2.py
conda activate rigel && pytest tests/ -v
```

Then rerun the eight-condition hybrid-capture suite and regenerate the condition
report. Golden outputs should be refreshed only after the eight-condition gates
are satisfied and the behavioral changes are accepted.

## 10. Deferred work

Unstranded high-gDNA capture-on remains an upstream source-deconvolution failure.
The diagnostics show `prior_mass` has only about 2.8k gDNA when truth is 100k.
No EM prior handoff, ESS, or grouped-prior smoothing can recover source mass that
never reaches `PriorMassDeconvolution`. That scenario needs a separate
fragment-length or capture-density source channel for near-unstranded data.

Precision-aware ESS scaling is also deferred. It is desirable, but only after
`prior_mass.precision` has a consistent meaning across strand, density, and
fallback methods.