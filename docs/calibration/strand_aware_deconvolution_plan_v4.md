# Strand-Aware SRD Calibration - v4 Plan

Status: design, 2026-05-12. Supersedes
[strand_aware_deconvolution_plan_v3.md](strand_aware_deconvolution_plan_v3.md).

v4 keeps the best parts of v3: aggregate raw moments, signed `p_r1_sense`,
canonical orientation routing, and no per-region clipping. It deliberately removes the
high-risk EXON BLUE estimator from the production path and narrows the first
implementation to the two channels that are actually contaminated by nRNA today:
`INTRON` contained counts and `EXON-INTRON` boundary flux.

## Part A - Review Findings From v3

### A.1 Keep the canonical orientation helper

The single most important implementation idea in v3 is still correct. There should be
one C++ helper and one Python helper for classifying a region/fragment strand pair:

```cpp
// src/rigel/native/calibration/orient.h
inline uint8_t classify_orient(uint8_t region_strand,
                               int8_t fragment_strand) noexcept;
```

```python
# src/rigel/calibration/_orient.py
def classify_orient(region_strand: int, fragment_strand: int) -> int: ...
```

Every call site must use these helpers. The value passed from the scanner must be
`result.exon_strand`, the same post-R2-flip convention used to train `StrandModel`.
The accumulator must not infer fragment strand from `obs_exons[0].strand`.

### A.2 Do not put exon-contained strand contrast in the production estimator yet

v3's main statistical overreach is the unified EXON BLUE estimator. It treats the
legacy boundary-flux estimator and a new exon-contained strand-contrast estimator as two
unbiased, independent observations of one exon gDNA density.

That is not a safe production assumption yet:

- Exon-contained counts are usually dominated by mature mRNA. The algebra is unbiased
  only if `p_r1_sense` is exact and globally transferable into every exon context.
  Small systematic errors in `p_r1_sense` are amplified by the very large exonic RNA
  count.
- Region-specific antisense transcription, overlapping annotations, and residual
  alignment-strand quirks become first-order errors in exon-only regions.
- The existing `EXON-INTRON` density is a boundary-flux density that is later projected
  onto exon body exposure. A direct exon-contained estimator is a different evidence
  source with different failure modes.
- v3's inverse-variance BLUE weights use observed-count variance. Zero-count channels
  can get pathological variance unless the implementation adds additional modeling
  guards.

Recommendation: v4 fixes nRNA contamination in the existing boundary channel by
strand-deconvolving `u_left` and `u_right`. Exon-contained strand contrast is kept as a
future diagnostic/experimental path, not part of the first production estimator.

### A.3 Do not use observed-count variances as channel-selection weights

The v2 review correctly removed per-region WLS. v3 reintroduces a related problem at
the channel level by using inverse-variance BLUE with observed-count plug-in variance.
For Poisson ratios, a zero observed count does not mean zero uncertainty; using it as a
weight can make zero-count channels dominate.

In v4, variances are diagnostics. They do not select between boundary and strand
estimators in the production path. The production point estimate is a direct ratio on
the relevant physical channel.

### A.4 Orientation-binned FL histograms are not free

`fl_hist` is one count per observed fragment mask. Per-region orientation arrays are
one count per qualified region hit. A fragment that overlaps multiple regions can have
multiple region orientations. Therefore a global `fl_hist_by_orient` cannot be filled
by simply classifying each region hit without double-counting the fragment.

Recommendation: do not add `fl_hist_by_orient` in Phase 1. If strand-decontaminated FL
is later needed, implement a separate fragment-level orientation histogram with an
agreement rule: a fragment contributes to SAME or OPP only when all informative
qualified regions contributing to that mask agree; otherwise it contributes to UNINF.

### A.5 Store only the orientation arrays needed by v4

Dense `(R, 8, 3)` orientation counts are larger than needed and increase native/Python
ABI surface. v4 needs only:

- intron contained counts by orientation for the `INTRON` density;
- exon boundary flux by orientation for the `EXON-INTRON` density.

This is smaller, clearer, and easier to validate.

### A.6 Configuration belongs under calibration

If a future opt-in strand-aware gDNA FL path is added, its flag belongs in
`CalibrationConfig`, not `EMConfig`. It changes calibration FL construction, not the EM
solver.

## Part B - Statistical Core

For a strand-unambiguous region or boundary side, define:

```text
Y_same = count where fragment_strand == region_strand
Y_opp  = count where fragment_strand != region_strand
T      = Y_same + Y_opp
D      = Y_same - Y_opp
p      = StrandModel.p_r1_sense
c      = 2p - 1
```

The RNA and gDNA observation model is:

```text
E[Y_same] = p       * N_r + 0.5 * N_g
E[Y_opp]  = (1 - p) * N_r + 0.5 * N_g
E[D]      = c * N_r
E[T]      = N_r + N_g
```

For `c != 0`, the raw gDNA moment for a count unit is:

```text
G_i_raw = T_i - D_i / c
```

These are raw signed moments, not local physical assignments. They may be negative or
larger than `T_i` in small regions. That is required for unbiased aggregation under
pure gDNA strand noise.

The density estimator is:

```text
G_raw  = sum_i G_i_raw
E_info = sum_i E_i
rho    = max(0, G_raw / E_info)
```

The only non-negativity clip is the final density clip. There is no per-region clip and
no per-region assignment of fragments to RNA or gDNA.

### Fallback Rules

Fallback is per density channel, not per run:

- If `abs(c) < c_min`, strand information is numerically unidentifiable. Use the legacy
  geometry estimator for that density.
- If strand-informative exposure is zero, use the legacy geometry estimator.
- If strand-informative exposure is positive and the observed count is zero, do not
  fallback. This is a valid zero-count rate observation: numerator zero, denominator
  positive.
- Rows with `RegionStrand.NONE`, `RegionStrand.AMBIG`, or fragment strand `NONE` route
  to `UNINF` and do not enter the strand numerator or denominator.

`c_min = 1e-3` is only a floating-point guard for the `1/c` terms. It is not a
statistical threshold.

### Small Counts and Overdispersion

The estimator should not decide region-by-region whether `3/0` is RNA or gDNA. It
aggregates signed moments across the channel and reports aggregate orientation evidence.

For diagnostics, compute the RNA-major contrast:

```text
delta   = abs(c)
Y_major = Y_same if p >= 0.5 else Y_opp
Y_minor = Y_opp  if p >= 0.5 else Y_same
D_major = sum_i(Y_major_i - Y_minor_i)
T_total = sum_i(Y_major_i + Y_minor_i)
```

Under a pure-gDNA orientation model:

```text
Y_major_total | T_total ~ Binomial(T_total, 0.5)
Var(D_major) ~= phi_orient * T_total
```

`phi_orient` starts at `1.0`. If real data show overdispersed orientation balance, it
can become a beta-binomial or empirical diagnostic later. The diagnostic answers the
small-count question: `3/0` has weak evidence; `100/0` has overwhelming evidence.

Do not use this diagnostic to clip each region. If a conservative production mode is
ever needed, it should be aggregate-only and explicitly named, for example:

```text
R_excess_lower = max(0, (D_major - z_q * sqrt(phi_orient * T_total)) / delta)
G_conservative = T_total - R_excess_lower
```

That is a future bias-variance tradeoff, not the v4 default.

## Part C - Density Channels

### C.1 INTERGENIC: unchanged legacy estimator

Intergenic regions have `RegionStrand.NONE`; there is no transcript strand against
which to classify RNA orientation.

Implementation remains the current contained-fragment ratio:

```text
rho_intergenic = sum(per_region_counts[:, MASK_INTERGENIC] on INTERGENIC rows)
                 / sum(L_eff on INTERGENIC rows)
```

### C.2 INTRON: strand-deconvolved contained-fragment ratio

Use only INTRON regions with `RegionStrand.POS` or `RegionStrand.NEG` and the contained
intron mask `MASK_INTRON`.

For each informative intron row:

```text
same_i = intron_counts_by_orient[i, ORIENT_SAME]
opp_i  = intron_counts_by_orient[i, ORIENT_OPP]
T_i    = same_i + opp_i
D_i    = same_i - opp_i
E_i    = L_eff_i
```

Then:

```text
rho_intron = max(0, sum_i(T_i - D_i / c) / sum_i(E_i))
```

If `abs(c) < c_min` or informative intron exposure is zero, use the current legacy
`_density_simple(..., MASK_INTRON)` result.

The UNINF counts are excluded from the strand estimator. Report both the informative
exposure fraction and the legacy collapsed density so users can see how much data the
strand channel used.

### C.3 EXON-INTRON: strand-deconvolved boundary flux

The existing `EXON-INTRON` density is a boundary-flux estimator on eligible exon sides.
v4 keeps that physical channel and decontaminates its numerator by orientation.

For each eligible exon side on an EXON row with unambiguous `RegionStrand`:

```text
same_side = u_side_by_orient[i, ORIENT_SAME]
opp_side  = u_side_by_orient[i, ORIENT_OPP]
T_side    = same_side + opp_side
D_side    = same_side - opp_side
E_side    = B_cross(splicing_anchor_tolerance)
```

Summing over eligible left and right sides:

```text
rho_exon_intron = max(0, sum_side(T_side - D_side / c) / sum_side(E_side))
```

If `abs(c) < c_min` or informative boundary exposure is zero, use the current legacy
`_density_exon_intron(...)` result.

This directly addresses the nRNA contamination that made boundary flux look like gDNA.
It also preserves the existing projection logic: `rho_exon_intron` is still the density
used to impute exon-contained gDNA mass downstream.

### C.4 Exon-contained strand contrast is deferred

Exon-only orientation counts can be useful diagnostics, but they should not feed the
production global gDNA density in v4. A future experimental estimator must include:

- sensitivity to `p_r1_sense` error, because exon RNA counts are large;
- robust variance that does not give zero-count channels infinite weight;
- covariance with boundary flux on short exons;
- validation on real high-expression genes and antisense-overlap loci.

## Part D - Implementation Phases

### Phase 1 - Native orientation routing with minimal payload

Behavior change: none. The scanner emits new orientation-resolved arrays; Python still
uses legacy collapsed arrays.

Files:

- `src/rigel/native/calibration/orient.h` (new)
- `src/rigel/native/calibration/region_index.h`
- `src/rigel/native/calibration/accumulator.h`
- `src/rigel/native/calibration/accumulator.cpp`
- `src/rigel/native/bam_scanner.cpp`
- `src/rigel/calibration/_orient.py` (new)
- `src/rigel/calibration/scan_payload.py`
- `src/rigel/pipeline.py`
- `tests/test_orient_routing.py` (new)
- `tests/test_calibration_accumulator.py`

Changes:

1. Add `ORIENT_SAME = 0`, `ORIENT_OPP = 1`, `ORIENT_UNINF = 2`, `ORIENT_N = 3` in C++
   and Python.
2. Add `classify_orient(region_strand, fragment_strand)` in C++ and Python with one
   routing matrix test covering all region/fragment pairs.
3. Extend `RegionIndex::set(...)` to accept `strands`; store `strands_`; expose
   `strand(rid)`.
4. Update `BamScanner.set_regions(...)` and `_wire_calibration_regions` to pass
   `region_df["strand"]` through the same sort permutation as starts/ends/type masks.
5. Extend `CalibrationAccumulator::observe(...)` with `int8_t fragment_strand`. The
   caller passes `result.exon_strand` at the current calibration observation site.
6. Add minimal orientation arrays:

```cpp
std::vector<int64_t> intron_counts_by_orient;  // (R, ORIENT_N), only MASK_INTRON
std::vector<int64_t> u_left_by_orient;         // (R, ORIENT_N)
std::vector<int64_t> u_right_by_orient;        // (R, ORIENT_N)
```

7. Increment legacy arrays exactly as today. Additionally:
   - when `obs_mask == MASK_INTRON`, route each intron per-region count by
     `classify_orient(regions.strand(rid), fragment_strand)`;
   - when incrementing `u_left` or `u_right`, route that boundary-side event by the
     exon region's strand and the fragment strand.
8. Validate in `CalibrationScanPayload.from_scan_dict(...)`:
   - `intron_counts_by_orient.sum(axis=1) == per_region_counts[:, MASK_INTRON]`;
   - `u_left_by_orient.sum(axis=1) == u_left`;
   - `u_right_by_orient.sum(axis=1) == u_right`.

Memory cost is about 72 bytes per region for the new dense arrays, compared with about
192 bytes per region for v3's `(R, 8, 3)` design. If this still matters, these arrays can
be compressed later because only INTRON and EXON rows consume them.

Exit criteria:

- legacy arrays are byte-identical to pre-Phase-1 output;
- orientation sums equal legacy arrays;
- one-thread and many-thread scans match on every new array;
- no Python density behavior changes yet.

### Phase 2 - Python strand-aware INTRON and EXON-INTRON densities

Behavior change: calibration densities use strand-aware estimators when available.

Files:

- `src/rigel/calibration/density_global.py`
- `src/rigel/calibration/_orchestrator.py`
- `src/rigel/calibration/_result.py`
- `src/rigel/pipeline.py`
- `tests/test_density_global.py`
- `tests/test_calibrate_orchestrator.py`

Implementation notes:

1. Add a small `StrandSummary` near `density_global.py` or in a new private module:

```python
@dataclass(frozen=True, slots=True)
class StrandSummary:
    p_r1_sense: float
    n_observations: int

    @property
    def c(self) -> float:
        return 2.0 * self.p_r1_sense - 1.0

    @classmethod
    def from_model(cls, model: StrandModel) -> "StrandSummary": ...

    @classmethod
    def uninformative(cls) -> "StrandSummary":
        return cls(p_r1_sense=0.5, n_observations=0)
```

2. Thread `StrandSummary.from_model(strand_models.exonic_spliced)` into
   `calibrate(...)` after strand models are finalized.
3. Keep legacy `_density_simple` and `_density_exon_intron` as fallback helpers. Do not
   require a large helper collapse in the same PR.
4. Add `_density_intron_strand_aware(...)` and
   `_density_exon_intron_strand_aware(...)`. Each returns a `GlobalGdnaDensity` plus
   additive diagnostics.
5. Preserve `GlobalGdnaDensity.n_fragments` as the legacy integer numerator for
   backwards-compatible summaries. Add fields or diagnostics for:

```text
rho_legacy
rho_strand
n_fragments_estimated
strand_source: legacy | strand
strand_informative_regions
strand_informative_sides
strand_informative_eff_length_bp
strand_informative_exposure_fraction
strand_orientation_z
strand_orientation_pvalue
```

6. Kappa: keep the existing kappa estimate based on collapsed legacy per-region counts
   in Phase 2 and label it as such in diagnostics if exposed. A robust kappa estimator
   for signed/fractional strand moments is separate work. The canonical global-only
   prior uses `rho`; kappa is currently locoregional diagnostic/shrinkage support.

Exit criteria:

- unstranded `p_r1_sense = 0.5` is bit-equivalent to legacy densities;
- R1-antisense mirror tests pass by carrying signed `c`;
- pure-gDNA strand simulations are unbiased after aggregate raw moments;
- pure-RNA intron/boundary simulations drive the corresponding gDNA density toward zero
  in stranded libraries;
- no per-region clipping is present.

### Phase 3 - Summary diagnostics and locoregional pass-through

Behavior change: none to priors or EM. Surface the information needed to audit the new
calibration.

Files:

- `src/rigel/calibration/_result.py`
- `src/rigel/calibration/locus_prior.py` (diagnostic fields only, if needed)
- `tests/test_calibration_result.py`
- `tests/test_per_locus_gdna_mass.py`

Additive summary fields:

```text
strand_aware: bool
p_r1_sense: float
strand_specificity: float
strand_n_observations: int

INTRON.rho_legacy
INTRON.rho_strand
INTRON.strand_source
INTRON.strand_informative_eff_length_bp
INTRON.strand_orientation_pvalue

EXON-INTRON.rho_legacy
EXON-INTRON.rho_strand
EXON-INTRON.strand_source
EXON-INTRON.strand_informative_sides
EXON-INTRON.strand_orientation_pvalue
```

`expected_gdna_count_global(...)` remains observation-free. Do not do locus-prior surgery
in this plan.

### Phase 4 - Optional strand-decontaminated gDNA FL, only if needed

This phase is out of the production critical path. It should not block Phase 2.

If implemented, put the flag in `CalibrationConfig`:

```python
strand_aware_gdna_fl: bool = False
```

Do not build this from per-region orientation arrays. Add a fragment-level FL
orientation histogram with exactly-one-bucket-per-fragment semantics:

- SAME/OPP only when all informative qualified regions for the fragment's gDNA-relevant
  mask agree on orientation;
- UNINF otherwise;
- cross-orient sum equals the legacy `fl_hist` for each supported mask.

Only then can `extract_gdna_counts(...)` optionally use raw aggregate FL moments. The
existing EB smoothing in `build_fl_models` can remain the finalizer.

## Part E - Validation Strategy

Unit tests:

- C++ and Python routing matrix for all `{NONE, POS, NEG, AMBIG} x {NONE, POS, NEG}`
  pairs.
- `result.exon_strand` is passed into `observe(...)`; no first-exon fallback remains.
- Cross-orient sums equal legacy arrays for intron counts and boundary flux.
- Worker merge byte equality includes every orientation array.
- `p_r1_sense = 0.5` falls back to legacy densities.
- R1-antisense mirror: `(p, same, opp)` and `(1-p, opp, same)` produce identical
  densities.
- Pure-gDNA many-small-region simulation is unbiased with aggregate raw moments.
- Per-region clipping regression test proves negative local moments survive into the
  aggregate sum.
- Zero-count informative exposure produces a valid zero numerator and positive
  denominator, not an artificial high-precision weight.
- Boundary channel deconvolution uses `u_left_by_orient` and `u_right_by_orient`, not
  exon-only counts.

Scenario tests:

- `tests/scenarios/test_nrna_double_counting.py`: remove xfails for `ss = 0.9` and
  `ss = 1.0` only after Phase 2 passes.
- Keep `ss = 0.65` as a stress case; record outcome before tightening expectations.
- Implicit-splice and clean `nrna_none` gates must not regress.

Benchmark checks:

- Stranded, gDNA-poor, nRNA-rich: `INTRON.rho` and `EXON-INTRON.rho` move toward truth.
- Stranded, gDNA-rich, nRNA-none: no regression versus current golden outputs.
- Unstranded: bit-equivalent to current outputs.
- Armis2 primary target: `gdna_*_ss_1.00_nrna_rand` reduces global gDNA over-allocation
  and improves mRNA recovery downstream.

## Part F - Risks and Mitigations

| Risk | v4 mitigation |
|------|---------------|
| Two orientation conventions diverge | One C++ helper, one Python helper, one routing matrix. |
| Fragment strand inferred incorrectly | `observe(...)` takes `result.exon_strand`; accumulator never derives it. |
| Per-region clipping bias | Raw moments aggregate first; final `max(0, rho)` only. |
| Observed-variance zero-count bias | No observed-count WLS or BLUE in production path. |
| EXON contained counts overreact to `p_r1_sense` error | Exon-only strand contrast deferred to a future experimental phase. |
| Boundary nRNA contamination remains | Boundary flux itself is orientation-deconvolved via `u_left/right_by_orient`. |
| Strand-informative subset is small or nonrepresentative | Report informative exposure fraction and legacy-vs-strand densities. |
| Kappa for signed moments is undefined | Keep legacy kappa in Phase 2; defer strand-aware kappa. |
| FL orientation histogram double-counts fragments | Do not add it in Phase 1; future FL histogram must be fragment-level. |
| Payload memory grows too much | Store channel-specific arrays `(R,3)`, not `(R,8,3)`. |

## Part G - Out of Scope

- Production EXON-contained strand-contrast estimator or EXON BLUE.
- Profiled-Poisson joint estimator with per-region RNA nuisance.
- Beta-binomial orientation overdispersion as a production correction.
- Strand-aware kappa estimation for signed/fractional moments.
- Locus-prior surgery using local strand evidence.
- EM-stage strand-model changes.
