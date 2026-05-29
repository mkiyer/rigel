# Strand-Based gDNA / RNA Deconvolution — Implementation Plan v1

Date: 2026-05-23
Status: implementation plan (active)
Scope: strand-specific RNA-seq libraries (capture or no-capture)
Design doc: `rnaseq_mode_aware_gdna_density_plan.md` (§4)
Related: `fine_region_phase4_calibration.md`

## 0. Naming Conventions (Read First)

All strand vocabulary in this plan and in the implementation is
**transcript-relative**. No "major/minor" terms. No "RNA-major" /
"RNA-minor". No "primary/secondary".

| Term | Meaning |
|---|---|
| `sense` | fragment aligned in the same orientation as the region's annotated transcript strand |
| `antisense` | fragment aligned opposite to the region's annotated transcript strand |
| `k_sense` | observed (fractional) count on the sense channel |
| `k_antisense` | observed (fractional) count on the antisense channel |
| `n` | `k_sense + k_antisense` |
| `p` | `p_r1_sense` from `StrandModel`: P(R1 aligns transcript-sense) |
| `kappa_d` | beta-binomial concentration for gDNA strand split |
| `R`, `D` | unknown RNA-fragment count and gDNA-fragment count (`D = n − R`) |

Protocol direction is encoded by `p` alone. We do **not** invent a rotated
frame; the model's RNA-sense emission probability is `p` for R1-sense preps
and small (close to 0) for R1-antisense preps, which is exactly what `p`
already is.  This keeps every variable name pointing at something a
biologist or alignment engineer recognises immediately.

Regions with `TS_AMBIG` (both transcript strands annotated) and `TS_NONE`
(intergenic, no transcript strand) are out of scope for the per-region
strand model. Intergenic regions still contribute to `kappa_d` training using
raw genomic strand counts (see §M2).

## 1. Model Recap (Single Source Of Truth)

For a strand-informative region (`TS_POS` or `TS_NEG`) with unspliced
fractional counts `k_sense`, `k_antisense`, and `n = k_sense + k_antisense`:

```text
R ~ prior on RNA fragment count, default Uniform{0, …, n}
D = n − R
DNA_sense | D ~ BetaBinomial(D, mean=0.5, concentration=kappa_d)
RNA_sense | R ~ Binomial(R, p)
k_sense   = DNA_sense + RNA_sense
```

Closed-form moment estimator (clipped to `[0, n]`):

```text
R_hat = (k_sense − 0.5·n) / (p − 0.5)
```

The denominator `(p − 0.5)` is signed, so this single expression handles
both R1-sense (`p ≈ 1`) and R1-antisense (`p ≈ 0`) libraries without any
rotation.

Conditional variance given `R`:

```text
Var(DNA_sense | D) = 0.25 · D · (D + kappa_d) / (1 + kappa_d)
Var(RNA_sense | R) = R · p · (1 − p)
Var(k_sense   | R) = Var(DNA_sense | D) + Var(RNA_sense | R)
```

Posterior used for confidence intervals:

```text
P(R = r | k_sense, n) ∝ P(k_sense | R = r, n, p, kappa_d) · P(R = r)
```

with the likelihood evaluated as the convolution of the gDNA beta-binomial
(over `D = n − r`) and the RNA binomial (over `r`).

Confidence-level outputs:

```text
rna_lower_α  = posterior_quantile(R; 1 − α)
gdna_upper_α = n − rna_lower_α
rna_mean     = E[R | k_sense, n]
gdna_mean    = n − rna_mean
```

These are the only numerical outputs the rest of the pipeline consumes.

## 2. Target Code Layout

New code lands under `src/rigel/calibration/`, alongside the existing
payload/region modules.

```
src/rigel/calibration/
├── strand_deconv/
│   ├── __init__.py              # public API re-exports only
│   ├── observations.py          # M1: protocol-aware observation builder
│   ├── kappa_d.py               # M2: symmetric beta-binomial MoM
│   ├── seed_training.py         # M2: high-purity training set assembly
│   ├── exon_screen.py           # M3: posterior screen for no-RNA exons
│   ├── posterior.py             # M4: per-region posterior + quantiles
│   ├── tail_approx.py           # M4: vectorized normal-tail fallback
│   └── result.py                # dataclasses: StrandObservations,
│                                #   StrandKappaFit, StrandDeconvResult
└── region_gdna.py               # M5: combined RegionGdnaEstimate builder
                                 #     (consumes StrandDeconvResult)
```

Why a sub-package: the strand deconvolution path has enough internal state
(observations → kappa → screen → posterior) that bundling it keeps
`calibration/` flat and lets us version the module as a unit. Public surface
is re-exported from `strand_deconv/__init__.py`; nothing else outside this
sub-package imports from its private modules.

Existing modules touched (additive only; no breaking changes in M1–M4):

| File | Change |
|---|---|
| `calibration/_arrays.py` | extend `PayloadArrays` with spliced sense/antisense projections (§M1) |
| `calibration/strand_summary.py` | add a `signed_q_sense` helper that just returns `p_r1_sense` |
| `calibration/_orchestrator.py` | wire M1–M5 between FL training and prior assembly |
| `calibration/_result.py` | add `StrandDeconvResult` field, JSON summary |

## 3. Public Dataclasses

All dataclasses are `@dataclass(frozen=True, slots=True)` and numpy-backed.
Sorted-position arrays follow `RegionArrays.order`.

### 3.1 `StrandObservations` (M1)

```python
@dataclass(frozen=True, slots=True)
class StrandObservations:
    """Per-region (sense, antisense) unspliced fractional counts.

    Built once per payload. Sorted in RegionArrays order.
    Ambiguous and no-strand regions are zero-filled and masked.
    """
    k_sense: np.ndarray          # float32[R] — unspliced sense, all compartments summed
    k_antisense: np.ndarray      # float32[R] — unspliced antisense, all compartments summed
    n_total: np.ndarray          # float32[R] — k_sense + k_antisense
    k_spliced_sense: np.ndarray  # float32[R] — diagnostic only
    k_spliced_antisense: np.ndarray  # float32[R] — diagnostic only
    strand_informative: np.ndarray   # bool[R] — TS_POS or TS_NEG, and n_total > 0
    p_r1_sense: float            # signed protocol direction from StrandModel
```

Compartment policy: `k_sense` aggregates `contained + boundary_left +
boundary_right`, all unspliced.  Boundary mass is fractional so this sum is
per-region-anchored and does not double-count fragments. Diagnostic
per-compartment views are produced by helper functions; they are not stored
on `StrandObservations` to keep the hot object small.

### 3.2 `StrandKappaFit` (M2)

```python
@dataclass(frozen=True, slots=True)
class StrandKappaSourceStat:
    source: Literal[
        "intergenic_unstranded",
        "intron_strand_specific",
        "boundary_intron_intergenic",
        "no_rna_exon",
    ]
    n_regions: int
    n_fragments: float
    kappa_d: float
    excess_variance_positive: bool

@dataclass(frozen=True, slots=True)
class StrandKappaFit:
    kappa_d: float               # always finite, clipped to [KAPPA_D_MIN, KAPPA_D_MAX]
    n_training_regions: int
    n_training_fragments: float
    source_breakdown: dict[str, StrandKappaSourceStat]
    fallback_used: bool
    fallback_reason: str
```

### 3.3 `StrandDeconvResult` (M4)

```python
@dataclass(frozen=True, slots=True)
class StrandDeconvResult:
    # sorted in RegionArrays.order
    n_total: np.ndarray          # float32[R]
    rna_mean: np.ndarray         # float32[R]
    gdna_mean: np.ndarray        # float32[R]
    rna_lower_95: np.ndarray     # float32[R]
    rna_lower_99: np.ndarray     # float32[R]
    gdna_upper_95: np.ndarray    # float32[R]
    gdna_upper_99: np.ndarray    # float32[R]
    precision: np.ndarray        # float32[R] — see §M4.4
    flags: np.ndarray            # uint8[R]   — see §3.4
    kappa_d: float
    p_r1_sense: float
    eligible_mask: np.ndarray    # bool[R]
```

### 3.4 Flag bits (`flags`)

```python
FLAG_INELIGIBLE           = 1 << 0  # not strand-informative or no counts
FLAG_LOW_COUNT            = 1 << 1  # n_total below the posterior-grid threshold
FLAG_HIGH_COUNT_TAIL      = 1 << 2  # fell back to normal-tail approximation
FLAG_NEAR_UNSTRANDED      = 1 << 3  # |p − 0.5| below identifiability floor
FLAG_POSTERIOR_DEGENERATE = 1 << 4  # numeric underflow during grid eval
FLAG_TRAINING_SEED        = 1 << 5  # region was used to train kappa_d
FLAG_NO_RNA_EXON          = 1 << 6  # region passed the M3 no-RNA screen
```

## 4. Milestone Details

### M1 — `observations.py`: protocol-aware observation builder

**Goal.** Produce `StrandObservations` from `PayloadArrays + RegionArrays +
StrandSummary` in one vectorized pass.

**Steps.**

1. Extend `PayloadArrays` (additive, no schema break): add
   `contained_spliced_pos/neg`, `boundary_left_spliced_pos/neg`,
   `boundary_right_spliced_pos/neg`. These come straight from the 12-channel
   matrix via `channel_index(..., SPLICE_SPLICED, ...)`.
2. For each compartment, call `sense_antisense_split(...)` once (unspliced)
   and once (spliced).  Sum the three `SenseSplit.sense` arrays into
   `k_sense`; same for antisense and spliced fields.
3. Compute `strand_informative = ((ts_class == TS_POS) | (ts_class == TS_NEG))
   & (n_total > 0)`.
4. Capture `p_r1_sense` from `StrandSummary`.

**Why not store per-compartment views.**  Downstream consumers (M3, M5)
re-use the same `(k_sense, k_antisense)` aggregation. Per-compartment access
during training/diagnostics is fine to recompute on demand because it's
cheap.

**Tests** (`tests/calibration/test_strand_observations.py`):

- single `TS_POS` region with all mass on POS channel → `k_sense == mass`,
  `k_antisense == 0`.
- single `TS_NEG` region with all mass on POS channel → `k_sense == 0`,
  `k_antisense == mass`.
- `TS_AMBIG` region → `strand_informative == False`, sense/antisense both 0.
- conservation: `sum(k_sense + k_antisense) ==
  sum(unspliced mass in strand-informative regions)`.
- R1-antisense library (`p_r1_sense = 0.05`): the *same* counts produce the
  same `StrandObservations`. Protocol direction must enter the model in M4,
  not here.

### M2 — `kappa_d.py` + `seed_training.py`: kappa_d training

**Goal.** Estimate `kappa_d` from high-purity gDNA-compatible regions using
the symmetric beta-binomial MoM already in
`density_global.estimate_strand_balance`. We are not re-deriving the
statistics; we are choosing a defensible training set and emitting a
per-source diagnostic.

**Training sources, in this order:**

1. **`intergenic_unstranded`** — `signature == 0x0`. Use raw genomic POS/NEG
   counts (the existing `estimate_strand_balance` consumes
   `contained_unspliced_pos/neg`). Highest purity.
2. **`intron_strand_specific`** — `signature ∈ {0x4, 0x8}` (`TS_POS` or
   `TS_NEG`, no exon bits). Use `(k_sense, k_antisense)`. Exclude `0xC`
   (`TS_AMBIG` intron). nRNA contamination is possible but small in capture.
3. **`boundary_intron_intergenic`** — strand-informative intronic regions
   neighbouring intergenic (derived from `left_signature`/`right_signature`).
   Boundary-side mass only.
4. **`no_rna_exon`** — exon regions passing the M3 screen. Folded in by M3
   re-running `fit_kappa_d` on the expanded set exactly once.

**Fit function.**

```python
def fit_kappa_d(
    x: np.ndarray,   # sense-side counts per training region
    n: np.ndarray,   # totals per training region
    *, source: str,
) -> StrandKappaSourceStat: ...
```

Uses the same `S, B, A` MoM as `estimate_strand_balance`. Clip with
`KAPPA_D_MIN = 1.0`, `KAPPA_D_MAX = 1.0e6`. Fall back to `KAPPA_D_MAX` and
set `excess_variance_positive=False` when `S ≤ B`.

**Pooling.** Pool sources by summing `S`, `B`, `A` rather than averaging
per-source `kappa_d`. This is the only statistically correct pool for MoM.
The pooled value is `StrandKappaFit.kappa_d`. Per-source values stay in
`source_breakdown` for diagnostics only.

**Conservative override.** After pooling, take
`kappa_d = min(pooled, min_source_kappa_d_with_n_regions_≥ K_MIN)`. This
guards against one very clean source (e.g. tiny intergenic set) pulling
`kappa_d` artificially high.  `K_MIN = 25` for the first cut; revisit after
real-data inspection.

**Tests** (`tests/calibration/test_strand_kappa_d.py`):

- pure binomial simulated data → `kappa_d` saturates to `KAPPA_D_MAX`.
- beta-binomial simulated data with known `kappa_true ∈ {10, 100, 1000}` →
  recovered `kappa_d` within ±20% across 1000 regions.
- empty training set → `fallback_used`, `kappa_d == KAPPA_D_DEFAULT`.
- two sources with very different `kappa_d` → pooled value is between them
  and the conservative override picks the smaller when both have ≥ `K_MIN`.

### M3 — `exon_screen.py`: no-RNA exon selection

**Goal.** Find strand-informative exon regions whose observed
`(k_sense, k_antisense)` is plausibly explained by gDNA alone, and fold them
into `kappa_d` training.

**Eligibility filter.** A candidate region must satisfy:

- `ts_class ∈ {TS_POS, TS_NEG}` (no `TS_AMBIG`);
- exon bit set in `signature`;
- `n_total ≥ N_MIN_EXON` (initial: `N_MIN_EXON = 10`);
- `k_spliced_sense + k_spliced_antisense ≤ S_MAX_SPLICED` (initial:
  `S_MAX_SPLICED = 1`); reject anything with clear RNA evidence;
- no opposite-strand transcript annotation overlap (`TS_AMBIG` neighbours
  via `left_signature`/`right_signature` are OK as neighbours, but the
  region itself must be strand-pure).

**Screen.** Use the seed `kappa_d` (M2 pooled, before exon folding) and the
M4 fast tail path to compute `P(R > R_THRESH | k_sense, n)` with
`R_THRESH = max(1, ceil(0.05 · n))`. Accept the region as `no_rna_exon`
training only if this tail probability is below `P_REJECT = 0.01`.

**Refit.** Call `fit_kappa_d(...)` once more with the no-RNA exon set
appended to the seed set (sufficient statistics simply sum). Do not iterate
further; this avoids the self-training feedback loop the design doc called
out.

**Implementation note.** The screen reuses M4's tail approximation (not the
grid). The grid is reserved for per-region output where higher accuracy
matters.

**Tests** (`tests/calibration/test_strand_exon_screen.py`):

- exon with 50/50 split and `n=100` → accepted as no-RNA.
- exon with 90/10 sense/antisense in an R1-sense `p=0.99` library →
  rejected (clear RNA signal).
- exon with 1 spliced fragment → rejected by spliced filter regardless of
  unspliced split.
- final `kappa_d` after refit is ≥ seed `kappa_d` when added exons are
  clean (more data → tighter estimate).

### M4 — `posterior.py` + `tail_approx.py`: per-region deconvolution

**Goal.** Compute `rna_mean`, `rna_lower_{95,99}`, `gdna_upper_{95,99}`,
`precision`, and `flags` for every strand-informative region.

**Two evaluation paths**, switched per-region:

#### M4.1 Posterior grid (the principled path)

Use when `n_total ≤ N_GRID_MAX` (initial: `N_GRID_MAX = 2000`).

1. Round fractional counts: `k_int = round(k_sense)`,
   `n_int = round(n_total)`. Fractional counts arrive nearly integer for
   unspliced contained mass; boundary mass produces small residues which
   round safely. Set `FLAG_LOW_COUNT` if `n_int < 2`.
2. Enumerate `r ∈ {0, 1, …, n_int}`. For each `r`:
   - DNA component: `BetaBinomial(D = n_int − r, mean=0.5, kappa=kappa_d)`
     PMF over `d_sense ∈ {0, …, D}`.
   - RNA component: `Binomial(r, p)` PMF over `r_sense ∈ {0, …, r}`.
   - Likelihood `P(k_int | r)` = convolution at index `k_int`.
3. Normalize over `r` against `Uniform{0, …, n_int}` prior.
4. Compute `rna_mean = Σ r · P(r | k_int)` and the requested lower quantiles
   from the cumulative posterior.

**Vectorization.** Per-`r` convolution is `O(n)`; full sweep is `O(n²)`.
For `n ≤ 2000` and modern numpy this is ~4M ops per region — acceptable —
but batch regions of similar `n` to share the BetaBinomial PMF table and
use `scipy.signal.fftconvolve` when `n > 256`. Implementation note: do not
prematurely write a native kernel; profile first.

#### M4.2 Normal/saddlepoint tail (the fast path)

Use when `n_total > N_GRID_MAX` or when the M3 screen calls into M4.

1. For each candidate `r`, compute mean and variance from §1.
2. Evaluate the upper-tail probability of observing `k_sense ≥ k_int` under
   the normal approximation `N(μ(r), σ²(r))`.
3. Bisect over `r ∈ [0, n_int]` to find the largest `r` with
   `P(k_sense ≥ k_int | r) ≥ α` — that is `rna_lower_α`.
4. Saddlepoint refinement is optional in v1; the normal approximation is
   conservative (slightly heavier-tailed in absolute terms) when `kappa_d`
   is small, which is the right direction for an RNA lower bound.
5. `rna_mean` from the closed-form moment estimator, `gdna_*` by
   complement.

Set `FLAG_HIGH_COUNT_TAIL` for regions evaluated this way.

#### M4.3 Identifiability guard

If `|p − 0.5| < STRAND_CONTRAST_NUMERICAL_FLOOR` (already defined in
`strand_summary.py` as `1e-3`), the model is non-identifiable. Set
`FLAG_NEAR_UNSTRANDED`, return `rna_mean = 0`, `rna_lower_* = 0`,
`gdna_mean = n_total`, `gdna_upper_* = n_total`, `precision = 0`. This is
the correct conservative output: with no strand contrast we cannot assert
RNA.

#### M4.4 Precision

```text
precision = 1 / (1 + Var(R | k_sense, n_total))
```

computed from the posterior in M4.1 or from the bracket width in M4.2.
Stored on `StrandDeconvResult.precision` and forwarded to EM prior
assembly.

**Tests** (`tests/calibration/test_strand_posterior.py`):

- closed-form check on a tiny grid (`n=4`, `kappa_d=∞`) against analytic
  binomial convolution.
- monotonicity: increasing `k_sense` at fixed `n` produces non-decreasing
  `rna_mean` and `rna_lower_95` (for `p > 0.5`).
- protocol symmetry: `p` and `1 − p` give the same `rna_lower_95` when the
  inputs are reflected (`k_sense ↔ k_antisense`).
- low-`n` calibration: simulate 10000 regions from the model with known
  `R_true`; `rna_lower_95` should cover `R_true` with frequency ≥ 0.95.
- high-`n` fast path agrees with grid path within 1% on regions sized at
  `N_GRID_MAX − 1` vs `N_GRID_MAX + 1`.
- `|p − 0.5| < 1e-3` → all RNA outputs zero, `FLAG_NEAR_UNSTRANDED` set.

### M5 — `region_gdna.py`: integration into `RegionGdnaEstimate`

**Goal.** Wire `StrandDeconvResult` into a `RegionGdnaEstimate` table that
downstream code (per-locus prior assembly) can consume without knowing
about the strand model.

**Mapping.**

| `RegionGdnaEstimate` field | Source from `StrandDeconvResult` |
|---|---|
| `mean_count` | `gdna_mean` |
| `upper_count` | `gdna_upper_95` (default; 99 stored for diagnostics) |
| `rna_lower_count` | `rna_lower_95` |
| `exposure` | placeholder (filled by M6 capture exposure; for now `l_eff_contained(span)` so the field is present and finite) |
| `precision` | `StrandDeconvResult.precision` |
| `flags` | merged with region-class flags |

For regions where `strand_informative == False` (intergenic, ambiguous,
strand-unidentifiable), `RegionGdnaEstimate` falls back to the existing
pooled estimator from `density_global.compute_global_densities`. This keeps
EM happy while the density-based fallback (design doc §M7) is built.

**Orchestrator wiring** (`calibration/_orchestrator.py`):

```text
payload = ...                                 # existing
strand_summary = StrandSummary.from_model(scan_trained.strand_model)
obs            = build_strand_observations(region_arrays, payload_arrays, strand_summary)
kappa_seed     = fit_kappa_d_seed(obs, region_arrays, payload_arrays)
no_rna_exons   = screen_no_rna_exons(obs, kappa_seed, strand_summary)
kappa_fit      = refit_kappa_d(kappa_seed, no_rna_exons)
deconv         = deconvolve_regions(obs, kappa_fit, strand_summary)
region_gdna    = build_region_gdna_estimate(
    region_arrays, payload_arrays, deconv, fallback_densities=global_table
)
result = CalibrationResult(..., strand_deconv=deconv, region_gdna=region_gdna)
```

The existing `FractionalCutoverPending` fail-fast in `quant_from_buffer` is
removed only after M5 lands and `region_gdna` is consumed by the prior
table.

**Tests** (`tests/calibration/test_region_gdna_strand.py`):

- strand-informative regions: `mean_count + rna_mean == n_total` within
  float tolerance.
- ambiguous regions: `mean_count` comes from the pooled fallback, not
  zeros.
- summary JSON contains `kappa_d`, `p_r1_sense`, per-source kappa
  diagnostics, and per-flag region counts.

## 5. Edge Cases & Numerical Safeguards

- **Fractional counts.** Round-to-integer is acceptable because the
  accumulator's fractional residues are small relative to total mass.
  Track `frac_residue = n_total − n_int` and skip the posterior path
  (use `FLAG_LOW_COUNT`) when `n_int < 2`.
- **Zero-count regions.** `strand_informative == False`; carry through as
  fallback. No NaNs allowed in `RegionGdnaEstimate`.
- **Extreme `p`.** `p ∈ {0, 1}` is numerically fine: RNA mass goes
  entirely to one channel, `Var(RNA_sense | R) = 0`, and all uncertainty
  comes from the DNA beta-binomial. The closed-form
  `R_hat = (k_sense − 0.5n) / (p − 0.5)` is well-defined for `p ∈ {0, 1}`.
- **`kappa_d` saturation.** When all training data is consistent with
  binomial (no overdispersion), `kappa_d = KAPPA_D_MAX`. The model behaves
  like the binomial-DNA case; that is correct, not a bug.
- **Determinism.** All numpy operations are deterministic in sorted region
  order. No random sampling in the posterior path.

## 6. Performance Targets

On a 30M-fragment human run with the current fine-region partition
(~1.5M regions, ~600k strand-informative):

- M1 observation build: ≤ 1 s (one pass over `region_counts`).
- M2 seed kappa: ≤ 100 ms (pooled MoM).
- M3 exon screen: ≤ 2 s (normal-tail bisection, vectorized).
- M4 deconvolution: ≤ 30 s total budget, split as
   - grid path: ≤ 25 s for the ≤ 5% of regions with `n_int > 256`
     (FFT conv);
   - tail path: ≤ 5 s for the rest.
- M5 assembly: ≤ 1 s.

If profiling exceeds the M4 budget, hoist the per-`r` PMF tables into
shared arrays before considering a native kernel. We are not writing C++
here unless the Python+numpy path measurably blocks `rigel quant`
end-to-end.

## 7. Implementation Order & Acceptance Gates

Each milestone is a self-contained PR. Don't fold M3+M4 into one PR; the
exon screen depends on M4's tail approximation but should ship behind a
feature flag if M4 is delayed.

| Step | Gate to merge |
|---|---|
| M1 | unit tests above pass; `StrandObservations` round-trips through pickle |
| M2 | simulated kappa recovery within ±20% at `kappa ∈ {10, 100, 1000}`; pooled MoM matches per-source MoM on single-source inputs |
| M3 | on the unstranded-no-capture synthetic benchmark, ≤ 5% of true RNA exons get added to no-RNA training; ≥ 80% of true gDNA-only exons are added |
| M4 | coverage tests above; protocol symmetry; grid/tail agreement |
| M5 | full `rigel quant` succeeds on a strand-specific capture fixture; benchmark deltas reported in PR |

## 8. Out Of Scope For v1

- Capture exposure model (`E_capture_R`, `eta_capture`) — design doc §6.
- Ambiguous-region density fallback — design doc §7 / §M7.
- Beta-binomial RNA component — only if spliced diagnostics show
  overdispersion; current assumption is Binomial RNA.
- Spliced-RNA expression prior on `P(R)` — uniform prior is the v1
  default.
- Mappability-corrected exposure — tracked separately.

## 9. Open Questions Tracked During Implementation

1. `N_GRID_MAX` ceiling: 2000 is a first guess. Profile and revisit after
   M4.
2. `K_MIN = 25` regions per source for the conservative override — is this
   adequate? Inspect real-data distributions during M2 PR review.
3. Should `no_rna_exon` regions be flagged in `RegionGdnaEstimate.flags`?
   Default in plan: yes (`FLAG_NO_RNA_EXON`), because downstream EM may
   want to weight these stronger.
4. Boundary-side training (`boundary_intron_intergenic`) needs
   side-specific `(sense, antisense)` projection — confirm
   `sense_antisense_split` can be called per-compartment without producing
   a redundant `SenseSplit` for the other side. Audit during M2; if it's
   a copy concern, add a `compartments=` parameter.

## 10. Done When

- `rigel quant` runs end-to-end on a strand-specific hybrid capture BAM
  without `FractionalCutoverPending`.
- `summary.json` contains a `strand_deconv` block with `kappa_d`,
  `p_r1_sense`, per-source diagnostics, region-count breakdown by flag.
- Golden fixtures regenerated; `tests/calibration/` strand-deconv suite
  green; existing tests still pass.
- A benchmark PR comparing the unmoved baseline (current pooled estimator
  fallback path) to the strand-deconvolved path on the lab's
  strand-specific capture truth data, with absolute and relative error
  reported per region class.
