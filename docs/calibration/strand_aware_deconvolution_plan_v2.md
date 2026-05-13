# Strand-Aware SRD Calibration — v2 Plan

Status: design, 2026-05-11. Supersedes
[`strand_aware_srd_plan.md`](strand_aware_srd_plan.md).

This document is a **review + revised implementation plan**. The v1 design is
statistically sound in its core derivation but has several gaps and
opportunities to simplify the existing calibration code. v2 keeps v1's
estimator math and tightens the parts that are easy to get wrong in code.

---

## Part A — Review of the v1 plan

### A.1 What v1 gets right

- Routing through `p_r1_sense` rather than `strand_specificity` (preserves
  R1-antisense protocols).
- Aggregate raw-moment ratio `ρ̂ = Σᵢ(Tᵢ − Dᵢ/c) / Σᵢ Eᵢ`, with the only
  non-negativity clip applied once at the end. The argument against
  per-region clipping is correct: clipping creates a one-sided bias under
  pure-gDNA noise.
- Identifying that the variance collapse `Var(N̂g,ᵢ) ∝ ρEᵢ(1 + 1/c²)`
  algebraically reduces a WLS estimator to the OLS ratio, **so no ad-hoc
  variance floor or `Wss` ramp is needed**.
- Routing `RegionStrand.AMBIG` and `RegionStrand.NONE` to a `UNINF` bucket
  rather than fabricating an answer.
- Treating strand evidence as a *refinement of the same observations* rather
  than an independent likelihood (no double counting).

### A.2 Statistical issues that need explicit treatment

**Issue 1 — variance collapse holds only under the gDNA-only null.**
The derivation in v1 §4 uses `µsame = µopp = 0.5 ρ E`, which is the gDNA-only
mean. Under the full mixture model

$$
Y_{\text{same}} \sim \text{Poisson}(\tfrac12 \rho E + p R), \qquad
Y_{\text{opp}}  \sim \text{Poisson}(\tfrac12 \rho E + q R)
$$

the per-region variance becomes

$$
\mathrm{Var}(\hat N_g) =
  (1 - 1/c)^2 (\tfrac12 \rho E + p R) +
  (1 + 1/c)^2 (\tfrac12 \rho E + q R),
$$

which depends on the per-region nuisance `R`. Genic introns where nRNA
dominates therefore have **higher variance than the WLS-collapse argument
suggests**. The aggregate ratio `ρ̂` remains *consistent and unbiased*
(linearity of expectation), but its standard error and any reported
confidence interval must use this RNA-aware variance, not the gDNA-only
one. v1's stated CI intervals will be too tight in nRNA-heavy runs.

Action in v2: state the consistency result clearly, but compute reported
uncertainty from the **observed mixture variance** (use observed counts as
plug-in for the Poisson means). The point estimator is unchanged.

**Issue 2 — strand-model uncertainty term scaling.**
v1 carries `(2D/c²)² Var(p)`. With production-scale strand training
(`n_obs ≳ 10⁶`, `p ≈ 0.95`), `Var(p) ≈ 5e-8`, while `(2D/c²)²` is bounded
above by `4 (Σ Tᵢ)² / c⁴`. For very stranded data the term is negligible;
for `c → 0` the per-region count term diverges first. v2 keeps the
strand-model variance term in the public API but documents that the
**observation-noise term dominates in all realistic regimes**, so v1
implementations may safely set the strand-model term to zero behind a flag
without affecting decisions.

**Issue 3 — gDNA fragment-length model is trained from contaminated channels.**
The current scanner builds the gDNA FL distribution from
`INTERGENIC + INTRON + EXON-INTRON` histograms. In nRNA-heavy stranded
runs, the intron and boundary histograms contain RNA. That biases `gdna_fl`,
which in turn biases `L_eff_contained` and `b_cross` — the **denominators**
of every density. Strand-aware **numerators** without strand-aware FL
training only half-fixes the problem. v1 defers this to "Phase 4-ish";
v2 promotes orientation-binned FL histograms into Phase 2 (the same
native PR that adds orientation counts) and fits a deconvolved gDNA FL
in Phase 3 alongside the deconvolved densities. Without this, the v1 plan
will likely show *partial* gains on the
`tests/scenarios/test_nrna_double_counting.py` xfailed cases.

**Issue 4 — exon-only strand evidence is excluded but is the strongest signal.**
v1 §5 declines to use exon-region strand counts to estimate gDNA density,
because high-expression exons amplify small `p` errors. That is correct as
stated, but it leaves on the table the most informative channel for **gDNA
detection** (not density) — pure-mRNA exons have ratio `Ymajor/Yminor ≈ p/q`,
and any deviation from that ratio is exon-localized gDNA. v2 keeps the
density estimator off exons (consistent with v1) but exposes
`exon_strand_minor_excess` as a QC diagnostic. Promoting it to a density
channel is a separate, validated phase.

**Issue 5 — gDNA orientation is assumed 50/50, but never measured.**
v1 §10 calls this out as "monitor explicitly". Make it a first-class
diagnostic in v2 Phase 3, computed from intergenic genomic strand balance.
If a real run shows >55/45 imbalance, the equations need `α ≠ 0.5` and we
should warn rather than silently bias.

**Issue 6 — fragment strand provenance.**
v1 correctly says use `result.exon_strand` (the convention already used by
strand-model training) rather than `obs_exons[0].strand`. We must extend
the C++ `CalibrationAccumulator::observe(...)` signature to accept this
strand, **not** infer it inside the accumulator. Doing the inference in two
places is the most likely future bug.

### A.3 Algorithmic / design issues

- v1 splits "carry strand into RegionIndex" (Phase 1) and "add orientation
  buckets" (Phase 2). Phase 1 produces an unused field. v2 merges them into
  one native PR — strand only enters the codebase once it is consumed.
- v1 does not state how the `(R, 8, 3)` orientation array is consumed for
  fragment masks that are not relevant to gDNA channels (e.g. EXON-only
  fragments). v2 either restricts orientation counting to the masks that
  matter (`MASK_INTRON`, `MASK_INTRON|EXON`) or keeps the dense layout but
  documents that other slices are diagnostic only.
- v1's `GlobalGdnaDensity.n_fragments: int` does not accommodate fractional
  deconvolved counts. v2 introduces a parallel `n_fragments_estimated:
  float` and keeps `n_fragments` as the raw integer numerator from the
  legacy collapsed array (so existing summaries are not rewritten).
- The current `density_global._density_simple` /
  `_density_exon_intron` are two near-duplicate functions parameterized
  by exposure and a numerator vector. v2 takes the strand work as the
  forcing function to **collapse them into one helper**, with the
  strand-deconvolution as a pluggable preprocessor over `(numerator,
  exposure)`.
- v1 doesn't specify how the no-strand-info fallback is selected per
  density type. A run with strong strand info on introns but tiny boundary
  counts should use strand-aware introns and legacy boundaries — i.e. the
  decision is per-channel, not per-run. v2 makes this explicit.

### A.4 Items removed in v2

- Variance floor on per-region weights (already removed in v1 — kept
  removed and explicitly regression-tested).
- `(SS - 0.55)/0.35` ramp.
- Per-region clipping of `N̂g,ᵢ` to `[0, T]`.
- Any consumption of strand evidence in the per-locus prior layer in
  Phase 1 of this plan. Local prior changes are deferred to v2 Phase 4
  diagnostic-only and a future plan for prior surgery.

---

## Part B — Revised implementation plan

The v2 plan is **three phases of code change plus one diagnostic-only
phase**. Phase 1 is a single native+payload PR; Phase 2 is a single Python
estimator PR; Phase 3 is a contamination-aware FL PR; Phase 4 wires
diagnostics into `summary.json` and the locoregional prior.

### Phase 1 — Native: strand-carrying regions + orientation-binned counts

Single PR. After this lands, the scanner can produce orientation-resolved
counts but Python still uses the legacy collapsed arrays.

**Files**
- `src/rigel/native/calibration/region_index.h`
- `src/rigel/native/calibration/accumulator.h`
- `src/rigel/native/calibration/accumulator.cpp`
- `src/rigel/native/bam_scanner.cpp` (binding + observation site)
- `src/rigel/calibration/scan_payload.py`
- `src/rigel/pipeline.py` (`_wire_calibration_regions`)
- `tests/test_calibration_accumulator.py`

**Changes**

1. Extend `RegionIndex::set(...)` to take `const uint8_t* strands` lock-step
   with the existing arrays. Store `std::vector<uint8_t> strands_` and add
   `inline uint8_t strand(int32_t rid) const`.
2. In the nanobind binding, accept a `uint8` strand array on
   `BamScanner.set_regions(...)`. Validate length equality.
3. In `_wire_calibration_regions`, pass `region_df["strand"].to_numpy(
   np.uint8)` through the same `lexsort` permutation as the other arrays.
4. Add an orientation enum to `accumulator.h`:

   ```cpp
   namespace orient {
   constexpr uint8_t SAME   = 0;  // fragment strand == region tx strand
   constexpr uint8_t OPP    = 1;
   constexpr uint8_t UNINF  = 2;  // region NONE/AMBIG, fragment NONE, c==0 not handled here
   constexpr int     N      = 3;
   }
   ```

5. Extend `CalibrationPayload`:

   ```cpp
   // Same/Opp/Uninf orientation refinement of the legacy arrays.
   // Sums across the orient axis equal the corresponding legacy array
   // bit-exactly. Worker merge is element-wise additive.
   std::vector<int64_t> per_region_counts_by_orient;     // (R, mask::N_STATES, orient::N)
   std::vector<int64_t> u_left_by_orient;                // (R, orient::N)
   std::vector<int64_t> u_right_by_orient;               // (R, orient::N)
   std::vector<int64_t> fl_hist_by_orient;               // (mask::N_STATES, orient::N, kFlBins)
   std::vector<int64_t> intergenic_genomic_strand_count; // (2,) POS/NEG QC for §A.2-Issue 5
   ```

   Memory: `R × 8 × 3 × 8B ≈ 192 B/region`. For 10⁶ regions that is
   ~190 MiB peak per worker. Acceptable; if it becomes an issue, restrict
   to `MASK_INTRON` and `MASK_INTRON|MASK_EXON` only (the only fragment
   masks that contribute to gDNA channels) — that drops the factor 8 to 2.

6. Extend `CalibrationAccumulator::observe(...)` signature with **one new
   argument**: `int8_t fragment_strand` (POS=1, NEG=2, NONE=0). The caller
   passes `result.exon_strand`. The accumulator never derives this from
   `obs_exons[0].strand`.

7. In `observe(...)`:
   - Compute the legacy mask exactly as today; increment legacy arrays
     bit-exactly (this is the byte-equivalence regression guard).
   - For each per-region increment, look up
     `bucket = orient_of(region_strand, fragment_strand)`:
     - `region_strand ∈ {NONE, AMBIG}` or `fragment_strand == NONE`
       → `UNINF`.
     - else `SAME` if equal, `OPP` otherwise.
   - Add 1 into the orientation-resolved arrays.
   - For boundary `u_left`/`u_right` increments, route the same way using
     the EXON region's strand.
   - For FL histogram, use the fragment-level mask × orientation cell.
   - For each accepted observed fragment, also increment
     `intergenic_genomic_strand_count[fragment_strand - 1]` if the
     fragment's mask is `MASK_INTERGENIC` (gives a POS/NEG balance check
     for the gDNA orientation assumption).

8. Worker merge adds every new array element-wise.

9. `CalibrationScanPayload.from_scan_dict` validates the new arrays and
   asserts `sum(by_orient, axis=-1) == legacy` for every cell. This is
   the central correctness invariant; if a worker merge or routing bug
   creeps in, the scan fails fast.

**Tests**

- Strand round-trip: `region_df["strand"]` survives the lexsort and arrives
  in `RegionIndex::strand(rid)` in the right order.
- Routing matrix: deterministic synthetic fragments cover every
  `(region_strand, fragment_strand) ∈ {NONE, POS, NEG, AMBIG} × {NONE,
  POS, NEG}` pair.
- R1-antisense protocol: `p_r1_sense < 0.5` is not consumed yet but the
  routing is symmetric — verify the orientation arrays for an R1-antisense
  trace are the mirror of an R1-sense trace.
- Worker merge: 1-thread and 8-thread scans of the same BAM produce
  byte-identical orientation arrays.
- Sum-equals-legacy invariant on a real scan.

**Exit:** scan emits orientation-resolved arrays; legacy outputs are
byte-identical; `compute_global_densities` is unchanged.

---

### Phase 2 — Python: strand-aware density estimator + helper unification

Single Python PR. After this lands, calibration uses strand-aware
densities when strand information is available, and falls back to the
existing estimator per-channel when it is not.

**Files**
- `src/rigel/calibration/density_global.py` (rewrite, see below)
- `src/rigel/calibration/_orchestrator.py`
- `src/rigel/calibration/_result.py` (add fields, do not remove)
- `src/rigel/pipeline.py` (thread `StrandModel` summary into `calibrate`)
- `tests/test_density_global.py`
- `tests/test_calibrate_orchestrator.py`

**Refactor**

Collapse `_density_simple` and `_density_exon_intron` into a single
helper that operates on *(numerator, exposure)* arrays and is composed
with an optional strand-deconvolution preprocessor:

```python
@dataclass(frozen=True, slots=True)
class _Channel:
    label: DensityType
    numerator: np.ndarray          # (R,) int64, legacy counts
    exposure:  np.ndarray          # (R,) float64
    eligible:  np.ndarray          # (R,) bool
    # Orientation-resolved counts aligned to numerator. Each is (R,) int64.
    n_same:    np.ndarray | None   # None ⇒ no strand info available
    n_opp:     np.ndarray | None
    region_strand: np.ndarray | None  # (R,) uint8

def _density(channel: _Channel,
             *,
             p_sense: float,
             p_sense_var: float,
             ) -> GlobalGdnaDensity:
    ...
```

`_density` runs the strand-deconvolved estimator iff `channel.n_same is
not None`, the strand model is informative (`|c| ≥ c_min`, default
`c_min = 0.05`), and there is at least one strand-informative eligible
row. Otherwise it returns the legacy ratio. The decision is **per
channel**, not per run.

Estimator math (matches v1 §3):

```text
informative_i = eligible_i AND region_strand_i IN {POS, NEG}
                          AND (n_same_i + n_opp_i) is meaningful
T_i = n_same_i + n_opp_i
D_i = (n_same_i - n_opp_i) if region_strand_i == POS else (n_opp_i - n_same_i)
                          # convention: D_i has sign of "RNA-major minus RNA-minor"
                          # relative to p_sense's "same" definition

# Strand-uninformative rows still contribute their legacy numerator to the
# fallback total, NOT to the deconvolved sum. We choose ONE of two paths
# per channel based on the share of strand-informative exposure.
G_raw  = (T_i - D_i / c).sum() over informative rows
E_raw  = exposure_i.sum()      over informative rows
rho_strand = max(0.0, G_raw / E_raw) if E_raw > 0 else 0.0

# Variance, plug-in mixture form (Issue 1 fix)
mu_same_i = max(0.5 * rho_strand * exposure_i, eps) + ...    # see below
mu_opp_i  = symmetric
Var_per_i = (1 - 1/c)**2 * mu_same_i + (1 + 1/c)**2 * mu_opp_i
Var(G_raw) = Var_per_i.sum()
Var(rho_strand) = Var(G_raw) / E_raw**2
```

For the variance plug-in we use **observed counts** as the unbiased
estimate of `µ`: `mu_same_i ≈ n_same_i`, `mu_opp_i ≈ n_opp_i`. This is
the textbook plug-in for a Poisson rate and avoids re-introducing a
data-dependent weight inside the *point estimator*.

**Channel selection rule (per channel):**

```python
if n_same is None:                         # no orientation arrays
    return legacy
if abs(2*p_sense - 1) < c_min:             # essentially unstranded
    return legacy
informative_exposure = exposure[informative].sum()
if informative_exposure < 0.5 * exposure[eligible].sum():
    return legacy                          # too sparse to trust strand path
return strand_aware
```

The threshold `0.5` is conservative; tune from validation. The selection
is recorded in `GlobalGdnaDensity.method ∈ {"legacy", "strand_aware"}`.

**`GlobalGdnaDensity` schema additions (additive, no breakage):**

```python
@dataclass(frozen=True, slots=True)
class GlobalGdnaDensity:
    type: DensityType
    rho: float
    n_fragments: int                     # legacy integer numerator
    eff_length_bp: float
    n_regions_used: int
    kappa: KappaEstimate
    # New v2 fields
    method: Literal["legacy", "strand_aware"] = "legacy"
    n_fragments_estimated: float = 0.0   # = rho * eff_length_bp
    rho_var: float = 0.0
    n_strand_informative_regions: int = 0
    strand_orientation_pvalue: float | None = None  # H0: gDNA-only orientation
```

**Channel construction:**

| Channel       | numerator                                    | exposure          | eligible              | strand source                                               |
|---------------|----------------------------------------------|-------------------|-----------------------|-------------------------------------------------------------|
| INTERGENIC    | `prc[:, MASK_INTERGENIC]`                    | `L_eff_contained` | `type == INTERGENIC`  | strand always `NONE` ⇒ legacy only (and intergenic POS/NEG QC from genomic_strand_count) |
| INTRON        | `prc[:, MASK_INTRON]`                        | `L_eff_contained` | `type == INTRON`      | `prc_orient[:, MASK_INTRON, SAME/OPP]`                      |
| EXON-INTRON   | `u_left * 1_L + u_right * 1_R`               | `B_cross * sides` | `type == EXON & sides>0` | `u_{left,right}_by_orient`, summed per side                |

For EXON-INTRON, build per-side strand counts. A region with `1_L = 1_R
= True` contributes both sides to the numerator and both side-resolved
strand counts; the convention for "same vs opp" is relative to the EXON
region's transcript strand for both sides (the boundary direction does
not flip strand semantics).

**Fallback:** if every channel falls back to legacy, the function output
is bit-equivalent to today.

**Strand model plumbing:** add a tiny frozen summary struct so
`density_global` does not depend on the full `StrandModel`:

```python
@dataclass(frozen=True, slots=True)
class StrandSummary:
    p_r1_sense: float
    p_r1_sense_variance: float
    n_observations: int
    @classmethod
    def from_model(cls, m: StrandModel) -> "StrandSummary": ...
    @classmethod
    def uninformative(cls) -> "StrandSummary":
        return cls(0.5, 0.25, 0)
```

`scan_and_buffer` finalizes its `StrandModels` already; pipeline passes
`StrandSummary.from_model(strand_models.exonic_spliced)` (the
authoritative model) into `calibrate(...)`. If no model exists, pass
`uninformative()` and every channel falls back to legacy.

**Tests**

- `p_r1_sense=0.5` → every channel chooses legacy; outputs match the
  pre-PR golden bit-for-bit.
- `p_r1_sense=1.0`, pure RNA in introns → `rho_intron ≈ 0`, with
  `rho_intron_var` reported.
- `p_r1_sense=1.0`, pure gDNA in introns → `rho_intron` matches legacy
  within Poisson noise.
- `p_r1_sense=0.05` (R1-antisense) gives the same `rho` as the
  symmetric `p_r1_sense=0.95` case on the mirrored fragment set.
- Sparse strand-informative coverage triggers per-channel fallback.
- AMBIG and NONE region rows never enter the strand sum.
- `n_same.sum() + n_opp.sum() + n_uninf.sum()` equals the legacy
  `n_fragments` for every channel — wired as an assertion in `_density`.
- Per-region clipping regression: synthetic many-small-region pure-gDNA
  trace must give `rhô` within 2σ of truth; v1's per-region clip path
  fails this.

**Scenario tests**

- Target the xfailed cases in `tests/scenarios/test_nrna_double_counting.py`
  at `ss=1.0` and `ss=0.9`. v2 should remove the xfail.
- `ss=0.65` is a stretch — record the result but do not gate the PR on
  it; it is the "intermediate-information" regime where the per-channel
  fallback may engage.

**Exit:** strand-aware densities are produced when strand info exists;
legacy outputs are preserved otherwise; the scenario xfails are
addressed.

---

### Phase 3 — Strand-aware gDNA fragment-length model

Single PR, immediately after Phase 2. Necessary because Phase 2 fixes
numerators while denominators (`L_eff_contained`, `B_cross`) still depend
on a **contaminated** gDNA FL. Without this, the
`ss=1.0, nRNA-heavy` scenario shows a residual bias.

**Files**
- `src/rigel/calibration/_fl_sources.py`
- `src/rigel/calibration/fl.py`
- `src/rigel/calibration/_orchestrator.py`
- `tests/test_fl_models.py` (new orientation-deconvolution case)

**Changes**

1. Add `extract_gdna_counts_strand_aware(payload, *, p_r1_sense, c_min)`
   that, per FL bin, computes
   ```
   G_bin = sum over informative mask cells of: T_bin - D_bin / c
   ```
   using `payload.fl_hist_by_orient`. Sum is over `MASK_INTRON` and
   `MASK_INTRON|MASK_EXON` cells where the routing was strand-informative.
   Add `MASK_INTERGENIC` cells with no orientation deconvolution.
   Apply non-negativity clip per bin **only after** summing across mask
   cells (same rule as densities).
2. If `|c| < c_min` or strand-informative bins account for `< 50%` of
   total gDNA evidence, fall back to the legacy
   `extract_gdna_counts(payload)`.
3. Smoothing: the resulting per-bin counts are sparser than the legacy
   counts. The existing `build_fl_models` empirical-Bayes Dirichlet
   shrinkage handles this (it is what the EB step exists for); no new
   smoothing primitive needed. Verify `prior_ess` is still appropriate
   in tests.
4. Pipe `StrandSummary` through `calibrate(...)` to `build_fl_models`.

**Tests**
- Orientation-deconvolution at `p=1.0` reduces to `2 * minor_orientation_hist`.
- Synthetic trace with intron nRNA shows the deconvolved gDNA FL mean
  matches the true gDNA FL mean within tolerance, while the legacy FL
  mean is biased toward the (longer or shorter) RNA FL.
- All existing FL tests continue to pass when `p=0.5`.

**Exit:** the gDNA FL model used by `compute_global_densities` is
strand-decontaminated; the full nRNA-heavy scenario suite passes.

---

### Phase 4 — Diagnostics and locoregional pass-through (no behavior change)

Smallest PR. Surfaces what we now know.

**Files**
- `src/rigel/calibration/_diagnostics.py`
- `src/rigel/calibration/locus_prior.py` (read-only diagnostics)
- `tests/test_per_locus_gdna_mass.py`

**Changes**

1. Add `summary.json` keys (extending, not replacing, existing keys):
   ```
   strand_aware:                   bool
   p_r1_sense:                     float
   strand_specificity:             float
   strand_n_observations:          int
   intergenic_genomic_strand:      {pos: int, neg: int, p_value: float}
   INTRON.method:                  "legacy" | "strand_aware"
   INTRON.n_strand_informative_regions: int
   INTRON.rho_var:                 float
   INTRON.strand_orientation_pvalue: float | null
   EXON-INTRON.{same fields...}
   ```
2. Locoregional `n_gdna_intron` and `n_gdna_boundary_observed` get
   `_strand_aware` parallel fields when oriented data is available.
   These are **diagnostic only** in this PR — `expected_gdna_count_global`
   continues to drive priors. Surgery on the prior estimator is out of
   scope for this plan.

**Exit:** users can see how much information the strand channel
contributed, and the team has the data to decide whether to wire it
into priors next.

---

## Part C — Validation strategy

### Unit tests (must all pass)

- Routing-matrix tests in Phase 1.
- Single- vs multi-thread byte equality, including new arrays.
- `sum(orient axis) == legacy` invariant baked into payload validation.
- Phase 2 fallback selection at `c → 0`, `n_obs → 0`, sparse coverage.
- Phase 2 R1-antisense / R1-sense symmetry.
- Phase 2 per-region clip regression (forbids reintroduction of v1's
  earlier mistake).
- Phase 3 deconvolution-at-`p=1` closed-form check.
- Phase 4 schema additions are additive (no existing key removed).

### Scenario tests

- `tests/scenarios/test_nrna_double_counting.py`:
  remove xfail for `ss=0.9` and `ss=1.0` after Phase 3; record the
  `ss=0.65` outcome.
- Implicit-splice and clean `nrna_none` golden gates must not regress.

### Synthetic / benchmark

Run the standard sweep
(`scripts/sim/locus_sweep.py` configs under `scripts/benchmark/configs/`):

- Stranded, gDNA-poor, nRNA-rich: must drop INTRON/EXON-INTRON `rho`
  toward truth.
- Stranded, gDNA-rich, nRNA-none: must not regress vs golden.
- Unstranded: must match golden bit-for-bit (per-channel fallback).

Armis2 full-scale benchmarks
(`scripts/benchmarking/configs/default.yaml`):

- `gdna_*_ss_0.50_nrna_*` conditions must be unchanged.
- `gdna_*_ss_1.00_nrna_rand` conditions are the primary target — the
  expected outcome is reduced gDNA over-allocation in the global density
  block and improved mRNA recovery in the EM that consumes it.

### Acceptance gate

Phase 1 PR may merge as soon as the byte-equality and routing tests pass
(it does not change calibration behavior).

Phase 2 + Phase 3 PRs must merge together for the scenario xfails to
clear; merging Phase 2 alone is allowed but should keep the scenario
xfails in place pending Phase 3.

Phase 4 may merge any time after Phase 2.

---

## Part D — Risks carried from v1, with v2 mitigations

| Risk | v2 mitigation |
|------|---------------|
| Using `strand_specificity` instead of `p_r1_sense` | Estimator API takes `p_r1_sense` only; `strand_specificity` is a presentation field. Symmetric R1-antisense test in Phase 2. |
| Inferring fragment strand from `obs_exons[0]` | `observe()` signature requires fragment strand parameter; accumulator never derives it. |
| Per-region clipping bias | Single end-of-pipeline `max(0, ratio)`; explicit regression test. |
| WLS variance under nRNA-rich conditions too tight | Variance computed from observed-count plug-in mixture, not gDNA-only collapse. Documented as Issue 1. |
| Variance floor / data-dependent weights | Not implemented; estimator is the OLS ratio of raw moments. |
| Contaminated gDNA FL bias | Phase 3 deconvolves the FL too; no longer optional. |
| Fractional deconvolved counts | New `n_fragments_estimated: float` field; legacy `n_fragments: int` retained for back-compat. |
| Mega-payload memory for large genomes | Default dense `(R, 8, 3)`; documented escape hatch to restrict to relevant mask bits. |
| Per-run vs per-channel fallback | Decision is per channel; recorded in `GlobalGdnaDensity.method`. |
| Locus prior changes leaking into this work | Phase 4 is diagnostic only; prior surgery deferred. |
| gDNA orientation assumed 50/50 | Intergenic genomic strand QC counter shipped in Phase 1; surfaced in summary in Phase 4. |

---

## Part E — Out of scope

- Profiled-Poisson exact estimator with per-region RNA nuisance — the
  aggregate moment ratio is a consistent estimator and v2 takes it as
  the production path. Promoting to a profiled likelihood is a follow-up
  if validation shows residual bias.
- Strand-aware exonic gDNA channel — diagnostic only in v2.
- Beta-binomial orientation overdispersion — `phi_orient = 1` until QC
  data shows it is needed.
- Locus-prior surgery using strand-aware evidence — separate plan.
- EM-stage strand model changes — this plan is about calibration only.
