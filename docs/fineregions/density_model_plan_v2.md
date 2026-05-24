# Density Model Plan v2

Date: 2026-05-24
Status: implementation-ready
Supersedes: `density_model_plan_v1.md`
Companion: `calibration_cleanup_plan.md`

## 0. Scope and First Principles

This plan describes the smallest density model that:

1. Produces a per-region gDNA prior count `mu_r` and a conservative upper
   bound `up_r` independently of strand.
2. Handles both **no hybrid capture** and **hybrid capture** modes via a
   single switch chosen automatically from the calibration scan.
3. Plugs into the cleanup-era `pipeline.py` through a single new function
   `build_prior_table(...)` whose return shape is the four arrays EM
   already consumes.
4. Has no per-region flag taxonomy. Mode is global. Class is intrinsic.
   Uncertainty is one scalar.

No NB-GLM. No quantile regression. No exposure surface. No L-dependent
dispersion. No covariates. Just pooled Poisson rates and FL-aware
opportunity, with a capture-mode toggle that changes which anchors and
which opportunities are used.

If the simple version is insufficient on real data, density v3 can layer
NB or covariates on top without changing the consumer interface.

## 1. Inputs

After the cleanup plan lands, the density model consumes:

```text
RegionArrays           # _arrays.py, sorted region geometry + signature + ts_class
PayloadArrays          # _arrays.py, contained and boundary unspliced POS/NEG totals
FLModels.gdna          # fl.py, gDNA fragment-length model
GlobalDensityTable     # density_global.py, pooled INTERGENIC and INTRON-ONLY rho
fl helpers             # _exposure.py: l_eff_contained,
                       #              fractional_boundary_side_exposure,
                       #              gdna_eff_len_for_loci
```

No strand state. No exposure. No CLI knob is required at v2; one optional
override is exposed for testing only.

## 2. Two Modes, Detected Globally

The library is either "depleted background" (typical for hybrid capture)
or "uniform background" (typical for non-capture RNA-seq). The choice is
made once per run from the calibration scan.

### 2.1 Mode definition

```text
mode in {"uniform", "capture"}

uniform mode:
  contained intergenic and intron anchors carry the gDNA signal.
  per-region (N, L_eff) uses contained mass and contained opportunity.

capture mode:
  contained intergenic anchors are depleted by capture.
  per-region (N, L_eff) uses contained + boundary mass and contained +
  boundary opportunity. Class anchors come from regions adjacent to
  captured exons.
```

In both modes, the predicted Poisson rate is class-pooled:

```text
mu_r = rho_class[class(r)] * L_eff_r
```

The capture mode does not learn an "exposure surface". It just changes the
anchor set and includes boundary opportunity. This is the smallest change
that captures the depleted-background phenomenon, and it is reversible.

### 2.2 Detector

```text
rho_intergenic_contained  = pooled rho over INTERGENIC, contained only
rho_intron_contained      = pooled rho over INTRON-ONLY, contained only
rho_intron_boundary       = pooled rho over INTRON-ONLY, boundary only

If rho_intergenic_contained <= 0.25 * rho_intron_contained
   and rho_intron_boundary  >= 2.0 * rho_intron_contained:
    mode = "capture"
Else:
    mode = "uniform"
```

Both thresholds are constants in the new module. They can be promoted to
config later if we have a real reason. The default detection is logged in
`summary.json` along with the four pooled rhos so a human can sanity-check.

Override:

```text
density_mode_override in {None, "uniform", "capture"}
```

Set via `CalibrationConfig.density_mode = None` by default. CLI flag
`--density-mode {auto,uniform,capture}` is acceptable but optional. Auto
is the default.

## 3. Anchors and Observations

The model fits **two** Poisson rates: one for INTERGENIC, one for
INTRON-ONLY. These are the only classes used for fitting. Other classes
(EXON-CONTAINED, EXON-INTRON) only receive predictions.

### 3.1 Per-region observation

For every region `r` we compute one `(N_r, L_eff_r)` pair. The choice of
contained-only vs contained+boundary depends on mode.

Notation:

```text
n_contained_r     = contained_unspliced_total[r]
n_boundary_l_r    = boundary_left_unspliced_total[r]
n_boundary_r_r    = boundary_right_unspliced_total[r]

l_contained_r     = l_eff_contained(end[r] - start[r], gdna_fl)
l_boundary_l_r    = fractional_boundary_side_exposure(
                       left_neighbor_length[r], gdna_fl)
l_boundary_r_r    = fractional_boundary_side_exposure(
                       right_neighbor_length[r], gdna_fl)
```

`left_neighbor_length[r]` is the length of the adjacent region in
`RegionArrays`, capped at the gDNA FL support. It is computed once per
ref-CSR slice using `ref_offsets`. Implementation guidance:

```text
left_neighbor_length[r] = start[r] - end[r-1]    if r is not the first
                                                  region on its ref,
                       else 0
right_neighbor_length[r] = start[r+1] - end[r]   if r is not the last
                                                  region on its ref,
                       else 0
```

These bp values are then passed through
`fractional_boundary_side_exposure` to convert into FL-weighted boundary
opportunity. Side lengths exceeding the FL support are naturally clipped
by the helper.

### 3.2 Mode-specific aggregation

```text
uniform mode:
  N_r     = n_contained_r
  L_eff_r = l_contained_r

capture mode:
  N_r     = n_contained_r + n_boundary_l_r + n_boundary_r_r
  L_eff_r = l_contained_r + l_boundary_l_r + l_boundary_r_r
```

The boundary contribution is always included in **capture** mode, even
for intergenic regions. We do not gate by class: the same `(N, L_eff)`
rule is applied for fitting and for prediction. This keeps the math one
formula.

### 3.3 Anchor masks

```text
anchor_intergenic = is_intergenic(signature)
anchor_intron     = is_intron_only(signature)

both masks require:
  L_eff_r > L_EFF_MIN              (default 1.0)
  no NaNs in N_r or L_eff_r
```

`L_EFF_MIN` is a constant in the new module. Regions failing the mask are
not used for fitting; they still receive predictions.

## 4. The Model

Two Poisson rates, period.

```text
rho_intergenic = sum N_r over anchor_intergenic
               / sum L_eff_r over anchor_intergenic

rho_intron     = sum N_r over anchor_intron
               / sum L_eff_r over anchor_intron
```

Per-region prediction:

```text
class(r) selection:
  if is_intergenic(r):              class = "intergenic"
  elif is_intron_only(r):           class = "intron"
  elif is_exon_any(r):              class = "intron"   # see 4.1
  else:                             class = fallback class

mu_r  = rho_class[class(r)] * L_eff_r
var_r = mu_r                                  # Poisson
up_r  = poisson_upper_quantile(mu_r, alpha)   # see 4.2
```

`alpha` is fixed at 0.95 in v2. Not a CLI knob.

### 4.1 Why exon classes use the intron rate

Exonic regions are not in the anchor set. We have to predict their gDNA
density from something. The intron-only rate is the safest near-exon
proxy because:

- intron-only regions are gDNA-dominated;
- they are physically adjacent to exons;
- in capture mode they share the same boundary-flux opportunity profile.

Using intergenic for exon prediction in capture mode would systematically
underpredict gDNA in the captured part of the genome. Using intron-only
matches the local genomic neighborhood.

In **uniform** mode the two rates should be close; the choice rarely
matters. In **capture** mode the choice matters and intron-only is the
correct one.

### 4.2 Upper quantile

We do not need a real Negative Binomial. For each `mu_r` compute:

```text
up_r = poisson_upper_quantile(mu_r, alpha=0.95)
```

implemented as:

```text
if mu_r < 50:
  up_r = scipy.stats.poisson.ppf(alpha, mu=mu_r)
else:
  # normal approximation
  up_r = mu_r + 1.6449 * sqrt(mu_r)
```

The `1.6449` is the standard normal 95th-percentile. The threshold 50 is
where Poisson quantiles are well-approximated by the normal. Inside the
v2 module both branches return a float; rounding to integer counts is the
caller's choice and not required for EM.

`up_r` is reported as a diagnostic per region. It is **not** used to bound
or replace `mu_r` in the prior. The prior is `mu_r`. The upper bound is a
QA surface. This is intentional: the EM kernel already accepts a single
prior count per locus.

### 4.3 Fallback class

If a region matches none of `is_intergenic`, `is_intron_only`,
`is_exon_any`, fall back to the larger of `rho_intergenic` and
`rho_intron`. Log a one-time count of fallback regions. This is a small
defensive path; in practice nearly all regions match one class.

### 4.4 Failure modes

If either `rho_intergenic` or `rho_intron` cannot be estimated (zero
anchors, zero `L_eff`):

```text
rho_class = max(0.0, the other rho_class)
flag DensityResult.fit_quality = "weak"
```

If both fail:

```text
rho_intergenic = rho_intron = 0.0
flag DensityResult.fit_quality = "no_anchors"
mu_r = 0 for all regions
```

EM still runs; the prior is effectively uniform-zero. This must not raise.
Real-data triage uses `fit_quality` from `summary.json`.

## 5. Module Layout

New module: `src/rigel/calibration/density_model.py`. Single file. Target
under 350 lines including docstring and helpers.

### 5.1 Public surface

```python
@dataclass(frozen=True, slots=True)
class DensityObservation:
    n: np.ndarray              # float64[R]
    l_eff: np.ndarray          # float64[R]
    left_neighbor_len: np.ndarray   # int64[R], bp pre-clip
    right_neighbor_len: np.ndarray  # int64[R], bp pre-clip
    mode: Literal["uniform", "capture"]


@dataclass(frozen=True, slots=True)
class DensityModelFit:
    mode: Literal["uniform", "capture"]
    rho_intergenic: float
    rho_intron: float
    n_anchor_intergenic: int
    n_anchor_intron: int
    fit_quality: Literal["good", "weak", "no_anchors"]
    rho_intergenic_contained: float    # always populated for diagnostics
    rho_intron_contained: float
    rho_intergenic_boundary: float
    rho_intron_boundary: float
    mode_auto_detected: bool


@dataclass(frozen=True, slots=True)
class DensityPrediction:
    mu: np.ndarray             # float64[R]
    up95: np.ndarray           # float64[R]
    class_used: np.ndarray     # uint8[R], 0=intergenic, 1=intron,
                               #          2=fallback


def build_density_observation(
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    *,
    gdna_fl: FragmentLengthModel,
    mode: Literal["uniform", "capture"],
) -> DensityObservation: ...


def detect_density_mode(
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    *,
    gdna_fl: FragmentLengthModel,
    override: Literal["uniform", "capture"] | None = None,
) -> tuple[Literal["uniform", "capture"], DensityModelFit]: ...


def fit_density(
    obs: DensityObservation,
    region_arrays: RegionArrays,
    *,
    mode_info: DensityModelFit,
) -> DensityModelFit: ...


def predict_density(
    fit: DensityModelFit,
    obs: DensityObservation,
    region_arrays: RegionArrays,
) -> DensityPrediction: ...


def build_prior_table(
    loci,                       # sequence of MultiLocus
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    *,
    gdna_fl: FragmentLengthModel,
    fit: DensityModelFit,
    prediction: DensityPrediction,
    ref_lengths,
) -> PriorTable: ...
```

### 5.2 PriorTable

```python
@dataclass(frozen=True, slots=True)
class PriorTable:
    gdna_prior_count_em: np.ndarray   # float64[L]
    gdna_eff_len: np.ndarray          # float64[L]
    enable_gdna: np.ndarray           # uint8[L]
    gdna_prior_count: np.ndarray      # float64[L], same as _em in v2
```

This is the four-array contract `pipeline._run_locus_em_partitioned`
already accepts. No exposure-weighted vs raw split. No precision channel.

### 5.3 Allocation rule

For a `MultiLocus` containing genomic blocks `B` over reference ids:

```text
candidate regions = regions whose start/end intersects any block in B
gdna_prior_count_em[L] = sum_r mu_r * overlap_bp(r, B) / region_length(r)
gdna_eff_len[L]        = gdna_eff_len_for_loci(loci, ref_lengths, gdna_fl)
enable_gdna[L]         = 1 if gdna_prior_count_em[L] > 0 else 0
gdna_prior_count[L]    = gdna_prior_count_em[L]
```

This is base-pair prorating against the model **prediction** `mu_r`, not
against observed counts. Smearing the smooth `mu_r` is safe because
`mu_r` was already constructed to be locally uniform-per-bp under the
class rate. We are not reallocating spiky observed mass; we are
reallocating a per-bp predicted rate.

If a region overlaps multiple MultiLoci, its overlap shares sum to at most
`region_length(r)`; the formula conserves total predicted mass within
rounding.

### 5.4 Orchestrator wiring

`_orchestrator.calibrate(...)` adds the density model as a single block:

```python
mode, mode_info = detect_density_mode(region_arrays, payload_arrays,
                                      gdna_fl=fl_models.gdna,
                                      override=config.density_mode)
obs = build_density_observation(region_arrays, payload_arrays,
                                gdna_fl=fl_models.gdna, mode=mode)
fit = fit_density(obs, region_arrays, mode_info=mode_info)
pred = predict_density(fit, obs, region_arrays)
```

`CalibrationResult` gains two optional fields:

```python
density_fit: DensityModelFit | None
density_prediction: DensityPrediction | None
```

Summary JSON adds a `"density"` block:

```json
"density": {
  "mode": "capture",
  "mode_auto_detected": true,
  "fit_quality": "good",
  "rho_intergenic": 1.23e-05,
  "rho_intron":     1.15e-05,
  "n_anchor_intergenic": 42810,
  "n_anchor_intron":     31987,
  "rho_intergenic_contained": 3.0e-06,
  "rho_intron_contained":     1.0e-05,
  "rho_intergenic_boundary":  2.5e-05,
  "rho_intron_boundary":      2.7e-05
}
```

## 6. Pipeline Wiring

`pipeline.quant_from_buffer` replaces `build_uniform_prior_table` (the
cleanup-era placeholder) with:

```python
prior_table = build_prior_table(
    loci=loci,
    region_arrays=region_arrays,
    payload_arrays=payload_arrays,
    gdna_fl=fl_models.gdna,
    fit=cal.density_fit,
    prediction=cal.density_prediction,
    ref_lengths=ref_lengths,
)

batch_gdna_prior_count_em = prior_table.gdna_prior_count_em
batch_gdna_eff_len        = prior_table.gdna_eff_len
batch_enable_gdna         = prior_table.enable_gdna
batch_gdna_prior_count    = prior_table.gdna_prior_count
```

`_run_locus_em_partitioned` is not touched.

When `cal.density_fit is None` (density failed to initialize for any
reason), fall back to `build_uniform_prior_table`. EM still runs.

## 7. Tests

New test file `tests/test_density_model.py`. Target ~250 lines.

### 7.1 Unit tests

- `test_detect_mode_uniform_synthetic`:
  synthetic payload with equal contained intergenic and intron mass and
  small boundary mass; expect `mode == "uniform"`.
- `test_detect_mode_capture_synthetic`:
  synthetic payload with depleted contained intergenic mass and high
  boundary mass; expect `mode == "capture"`.
- `test_detect_mode_override`:
  forcing `override="uniform"` or `"capture"` bypasses detection and
  records `mode_auto_detected = False`.
- `test_build_observation_uniform`:
  in uniform mode, `obs.n` equals contained totals and
  `obs.l_eff` equals contained `l_eff` for every region.
- `test_build_observation_capture`:
  in capture mode, `obs.n` and `obs.l_eff` include boundary
  contributions.
- `test_left_right_neighbor_length`:
  first and last regions on a ref have zero neighbor length on the
  outside; interior regions get the correct gap-based bp value.
- `test_fit_density_recovers_constant_rho`:
  draw counts from `Poisson(rho * L_eff)` with known `rho`; fit recovers
  `rho` within 1 percent for 1000+ anchors.
- `test_fit_density_zero_anchors`:
  no INTRON-ONLY regions; `fit_quality == "weak"`, `rho_intron == 0`,
  no exception.
- `test_predict_uses_intron_rate_for_exons`:
  exon-class region predictions equal `rho_intron * L_eff`.
- `test_upper_quantile_branches`:
  for small `mu_r` the prediction uses `poisson.ppf`; for `mu_r >= 50`
  the prediction uses the normal approximation; both are within one
  count at the threshold.
- `test_build_prior_table_conservation`:
  a single `MultiLocus` that covers exactly the footprint of N regions
  has `gdna_prior_count_em == sum mu_r` over those regions within
  rounding.
- `test_build_prior_table_partial_overlap_prorates`:
  a `MultiLocus` covering half of one region inherits half of that
  region's `mu_r`.

### 7.2 Integration test

- `test_density_model_end_to_end_uniform`:
  tiny synthetic non-capture BAM, run `rigel quant`, assert the density
  block in `summary.json` reports `mode == "uniform"` and EM converges.
- `test_density_model_end_to_end_capture`:
  tiny synthetic hybrid-capture BAM, run `rigel quant`, assert the
  density block reports `mode == "capture"` and the gDNA prior is
  nonzero on capture-adjacent loci.

The benchmarking harness (`scripts/benchmarking/`) is the next acceptance
layer. Synthetic v2 should not regress gDNA assignment error on
non-capture cases and should improve it on capture cases.

## 8. Acceptance Criteria

```text
- src/rigel/calibration/density_model.py exists, under 350 lines
- DensityModelFit and DensityPrediction are populated on CalibrationResult
- summary.json contains the "density" block as specified
- detect_density_mode picks the right mode on the two synthetic fixtures
- rigel quant runs end-to-end in both modes
- on the lab's non-capture synthetic benchmark: gDNA mean relative error
  is no worse than the cleanup-era uniform-zero baseline
- on the lab's capture synthetic benchmark: gDNA mean relative error
  on captured loci improves measurably vs uniform-zero baseline
- pytest tests/ green
- ruff check src/ tests/ clean
```

## 9. What Density v2 Deliberately Does Not Do

- It does not model overdispersion. Poisson only.
- It does not learn an exposure surface.
- It does not fuse with strand evidence. Strand is diagnostic-only.
- It does not categorize regions beyond {intergenic, intron, exon-like,
  fallback}.
- It does not consume mappability or GC.
- It does not use exact discrete posteriors. Poisson and the normal tail
  are sufficient for the prior count and the QA upper bound.
- It does not introduce any new fragment-length or geometry helpers
  beyond what already exists in `_exposure.py`.

If real data demands more, density v3 layers on top:

```text
v3: NB with class-pooled dispersion          (replaces 4.2)
v4: separate boundary rho from contained rho  (modifies 2.2)
v5: NB-GLM with log L slope                   (replaces 4)
v6: covariates (GC, mappability)              (extends 4)
v7: unsupervised exposure surface             (replaces 2.2)
```

Each step is an isolated change to one section of this module. The
consumer interface (`PriorTable`) does not change.

## 10. Risks

- **Mode misclassification on borderline libraries.** The two-threshold
  rule is conservative. Mitigation: a CLI override and explicit
  `mode_auto_detected` in summary JSON. If misclassification becomes
  frequent on real data, replace the rule with a learned classifier,
  not new flags.
- **Intron-only rate biased by nRNA contamination.** Pre-existing
  problem. The cleanup explicitly leaves intron-only anchors in the
  pool; v2 inherits this. Mitigation: report
  `rho_intergenic_contained` and `rho_intron_contained` together so we
  can spot inflation. Density v3+ can subtract.
- **Boundary opportunity overcounted when adjacent regions are tiny.**
  `fractional_boundary_side_exposure` already clips by FL support, so
  this is bounded. Mitigation: spot-check the integration test fixtures
  and add an assertion that `L_eff_r >= l_contained_r` for every region
  in capture mode.

## 11. First Implementation Sequence

After the cleanup plan is fully merged:

```text
PR 1: scaffold density_model.py with DensityObservation,
      build_density_observation, detect_density_mode, unit tests for
      observation + mode detection.

PR 2: fit_density, predict_density, unit tests for fitting and
      prediction. CalibrationResult fields and summary JSON block.

PR 3: build_prior_table, pipeline wiring, fallback to uniform when
      density is unavailable, integration tests.

PR 4: synthetic benchmarks vs uniform-zero baseline. Document numbers in
      docs/benchmarks/.
```

Each PR is independently reviewable, independently revertable, and
exercises the same four-array `PriorTable` contract.
