# Strand-Aware SRD Calibration — v3 Plan

Status: superseded by
[`strand_aware_deconvolution_plan_v5.md`](strand_aware_deconvolution_plan_v5.md),
2026-05-12. Historical design, 2026-05-11. Supersedes
[`strand_aware_deconvolution_plan_v2.md`](strand_aware_deconvolution_plan_v2.md).

This revision incorporates feedback on v2:

1. Strand-model variance term is dropped from the production estimator
   (Issue 2 confirmed negligible at production training sizes).
2. nRNA-contaminated gDNA FL is real but rare in production data; the
   strand-decontaminated FL build is **moved to an opt-in diagnostic**,
   not a blocker.
3. EXON regions get a **single unified gDNA estimator** that combines
   the existing boundary-flux extrapolation with the new strand-contrast
   evidence via inverse-variance weighting. No magic-number cutoff.
4. The intergenic-orientation QC counter is **removed**. Genomic DNA is
   physically double-stranded; orientation is 50/50 by biology, not by
   measurement.
5. A **single canonical site** for "fragment strand vs region transcript
   strand" comparison is defined and reused everywhere, in C++ and
   Python.

---

## Part A — Updated review summary

### What carries over from v2

- Aggregate raw-moment ratio
  $\hat\rho = \max\!\bigl(0,\;\sum_i (T_i - D_i/c) \,/\, \sum_i E_i\bigr)$
  with the only non-negativity clip applied once at the end.
- Routing through `p_r1_sense` (preserves R1-antisense protocols).
- `RegionStrand.AMBIG` and `RegionStrand.NONE` route to a `UNINF` bucket;
  `UNINF` rows do not enter the strand estimator.
- Per-channel (not per-run) fallback to the legacy estimator when
  $|c| < c_\text{min}$ or strand-informative exposure is too sparse.
- Variance reported is the observed-count plug-in mixture variance,
  not the gDNA-only collapse (so reported CIs stay honest under
  nRNA-rich conditions).
- One-time refactor: collapse `_density_simple` and `_density_exon_intron`
  into a single `_density(channel, …)` helper.

### What v3 changes

- **Issue 4 — unified EXON estimator (NEW).** EXON regions now have
  *two* unbiased estimators of the same gDNA density `rho_g_exon`:
  the existing boundary-flux extrapolation and the new strand-contrast
  estimator. v3 combines them by **inverse-variance weighting** of the
  two `rho` estimates. This is the textbook BLUE for two unbiased
  observations of one parameter; it requires no thresholds and
  degrades gracefully when one channel is uninformative (its variance
  blows up and its weight goes to zero).
- **Issue 5 — drop intergenic POS/NEG QC counter.** Removed from the
  payload and from `summary.json`. gDNA is 50/50 by molecular biology.
- **Issue 3 — strand-aware gDNA FL is opt-in.** The orientation-binned
  FL histogram is still produced in Phase 1 (it is essentially free
  in the same loop). The deconvolved FL **builder** is implemented as
  an opt-in code path behind `EMConfig.strand_aware_gdna_fl: bool =
  False`, with a one-paragraph caveat in the docstring. Production
  default behavior is unchanged.
- **Issue 6 — single canonical strand-comparison site.** Both Python
  and C++ get a one-line helper:
  ```cpp
  // src/rigel/native/calibration/orient.h  (NEW)
  inline uint8_t classify_orient(uint8_t region_strand,
                                 int8_t  fragment_strand) noexcept;
  ```
  ```python
  # src/rigel/calibration/_orient.py  (NEW)
  def classify_orient(region_strand: int, fragment_strand: int) -> int: ...
  ```
  Every call site uses these. The convention is documented at exactly
  one location and tested with a routing matrix.

### What v3 removes from v2

- v2 Phase 3 ("strand-aware gDNA FL") as a *required* PR. v3 ships
  the histogram in Phase 1 and the *consumer* (deconvolved FL build)
  as an opt-in diagnostic in Phase 4.
- `intergenic_genomic_strand_count` from the payload.
- Any per-channel "0.5 of eligible exposure" magic threshold for the
  legacy fallback. v3 falls back **only** when `|c| < c_min` or when
  the channel literally has no strand-informative observations
  (`Σ T_i = 0`). Inverse-variance weighting handles the
  intermediate-sparsity case naturally — the strand-channel weight
  shrinks toward zero as counts shrink.
- The `0.5 * exposure[eligible].sum()` selection rule from v2 Phase 2.

---

## Part B — The unified EXON estimator

This is the substantive new statistics in v3. Read this section before
the implementation phases.

### B.1 Two unbiased estimators of `rho_g_exon`

For an EXON region partition with boundary-flux eligibility, define on
each eligible row $i$ the **exon-area gDNA density**
$\rho_g$ (fragments / bp).

**Channel 1 — boundary-flux extrapolation (legacy).** The current
`_density_exon_intron` estimator measures gDNA flux at exon-intron
edges and extrapolates the rate per bp into the exon body via the
FL-PMF crossing exposure $B_\text{cross}$:

$$
\hat\rho_{g,B} = \frac{\sum_i u_i \cdot 1^\text{eligible}_i}
                      {\sum_i s_i \cdot B_\text{cross}}, \qquad
\mathrm{Var}(\hat\rho_{g,B}) \approx
   \frac{\sum_i u_i \cdot 1^\text{eligible}_i}
        {(\sum_i s_i \cdot B_\text{cross})^2}
$$

where $u_i$ is the per-side flux count and $s_i \in \{1, 2\}$ is the
number of eligible boundary sides on row $i$. This is unchanged from
the legacy code; only the variance is now reported.

**Channel 2 — strand contrast (new).** On EXON regions whose
transcript strand is unambiguous, the existing strand identity holds:

$$
\hat\rho_{g,S} = \frac{\sum_i (T_i - D_i / c) \cdot \mathbf{1}^\text{exon-info}_i}
                      {\sum_i L^\text{eff}_i \cdot \mathbf{1}^\text{exon-info}_i}
$$

where the numerator is the strand-deconvolved exon gDNA mass and the
denominator is the FL-PMF containment effective length on that exon
region. Variance uses the observed-count plug-in mixture form (see
v2 §A.2 / Issue 1):

$$
\mathrm{Var}(\hat\rho_{g,S}) =
   \frac{\sum_i \bigl[(1 - 1/c)^2 \, n_{\text{same},i}
                    + (1 + 1/c)^2 \, n_{\text{opp},i}\bigr]}
        {\bigl(\sum_i L^\text{eff}_i\bigr)^2}.
$$

Both estimators target the same quantity: the gDNA density across the
exon body. They use **disjoint physical observations** (boundary
fragments vs contained-fragment strand counts), so their estimation
errors are approximately independent.

### B.2 Combining the two channels — BLUE

The Best Linear Unbiased Estimator of `rho_g_exon` from two unbiased,
approximately independent estimators is inverse-variance weighting:

$$
w_B = 1/\mathrm{Var}(\hat\rho_{g,B}), \qquad
w_S = 1/\mathrm{Var}(\hat\rho_{g,S}),
$$

$$
\hat\rho_g = \frac{w_B \,\hat\rho_{g,B} + w_S \,\hat\rho_{g,S}}{w_B + w_S}, \qquad
\mathrm{Var}(\hat\rho_g) = \frac{1}{w_B + w_S}.
$$

Properties that make this the right choice:

- **Unbiased**: a convex combination of unbiased estimators is unbiased.
- **No thresholds**: the only "switch" is `w = 0` when a channel has
  zero counts or `c = 0`, which falls out of the formulas.
- **Graceful degradation in both directions**: an unstranded library
  has $c \to 0$, $\mathrm{Var}(\hat\rho_{g,S}) \to \infty$, $w_S \to 0$,
  and the combined estimator collapses to the legacy boundary-flux
  estimator bit-equivalently. Conversely, an exon region without
  eligible boundary sides has $w_B = 0$ and the combined estimator is
  the strand-only estimator.
- **Optimal**: minimum variance among linear unbiased combinations.
- **Honest variance**: the reported $\mathrm{Var}(\hat\rho_g)$ shrinks
  exactly as expected when both channels carry weight.

The independence assumption is satisfied because the two channels
count different fragment populations: boundary flux uses fragments
whose alignment crosses an exon edge, strand contrast uses fragments
contained in the exon body. Edge cases (a fragment that spans an
entire short exon and crosses both edges) contribute to both channels;
in production these are a tiny fraction of fragments and the
correlation they induce is bounded by their mass share — small enough
to ignore in v3 and to revisit only if benchmarks show systematic
under-coverage of $\hat\rho_g$ confidence intervals.

### B.3 INTRON channel — same recipe, one channel only

INTRON has only the strand-contrast channel (there is no boundary
extrapolation for intron). The unified function therefore handles
INTRON as a degenerate case: $w_B = 0$, $w_S$ as above. INTERGENIC
similarly has only the legacy contained-fragment channel: $w_S = 0$
because intergenic regions have `RegionStrand.NONE`. The same
combinator covers all three densities.

### B.4 Falling back when strand information is unavailable

Definition: `c_min` is the minimum value of $|c| = |2p - 1|$ at which
the strand-channel variance is finite enough to contribute. The clean
mathematical condition is simply

```
w_S > 0  ⟺  c != 0 AND Σ (n_same + n_opp) > 0 over informative rows
```

i.e. the channel only weights itself when it has both strand
distinguishability and observations. We clamp `|c| ≥ c_min = 1e-3`
purely to avoid division blow-up in the variance term; this is a
numerical guard, not a statistical threshold (any value far below
production strand specificity values works).

---

## Part C — Convention: fragment strand vs region transcript strand

This is **the** authoritative paragraph. Every implementation site
must defer to it.

- `region_strand ∈ {NONE = 0, POS = 1, NEG = 2, AMBIG = 3}` matches
  `rigel.calibration.regions.RegionStrand`.
- `fragment_strand ∈ {NONE = 0, POS = 1, NEG = 2}` matches
  `rigel.types.Strand`. It is the **post-R2-flip exon-block strand**
  (the convention already used to train `StrandModel`). Specifically,
  the value passed to `CalibrationAccumulator::observe(...)` MUST be
  exactly `result.exon_strand` from the resolver, not derived from
  `obs_exons[0].strand`, not from the BAM flag directly, and not
  recomputed in C++.
- `p_r1_sense = P(fragment_strand == region_strand | RNA, region
  unambiguous)`. This is `StrandModel.p_r1_sense` after finalize.
  An R1-sense library has `p_r1_sense > 0.5`; an R1-antisense
  library has `p_r1_sense < 0.5`. The **same** `p_r1_sense` is used
  by every consumer; we do not invert it anywhere.
- The orientation classification is a single function:

  ```python
  # rigel/calibration/_orient.py
  ORIENT_SAME, ORIENT_OPP, ORIENT_UNINF = 0, 1, 2

  def classify_orient(region_strand: int, fragment_strand: int) -> int:
      """Classify a (region, fragment) strand pair for SRD deconvolution.

      ORIENT_SAME  : region_strand == fragment_strand and both ∈ {POS, NEG}
      ORIENT_OPP   : region_strand != fragment_strand and both ∈ {POS, NEG}
      ORIENT_UNINF : region NONE/AMBIG, or fragment NONE
      """
      if (region_strand == 1 or region_strand == 2) \
              and (fragment_strand == 1 or fragment_strand == 2):
          return ORIENT_SAME if region_strand == fragment_strand else ORIENT_OPP
      return ORIENT_UNINF
  ```

  The C++ analogue lives in `src/rigel/native/calibration/orient.h`
  and uses identical numeric values. The values are referenced
  symbolically (no bare `0/1/2` literals at call sites).
- The strand-deconvolution math always treats `n_same` as
  `ORIENT_SAME` and `n_opp` as `ORIENT_OPP`. The sign convention
  $D = n_\text{same} - n_\text{opp}$ matches `p_r1_sense` directly:

  $$
  \mathbb{E}[D] = (2p_\text{r1\_sense} - 1) \cdot R \;=\; c \cdot R
  $$

  with no per-region sign flip on `region_strand`. Both R1-sense
  (`p > 0.5`, `c > 0`) and R1-antisense (`p < 0.5`, `c < 0`) are
  handled by carrying signed `c` into `D / c`.

A CI test asserts the routing matrix end-to-end (Python helper, C++
helper, full scan) for all
$\{NONE, POS, NEG, AMBIG\} \times \{NONE, POS, NEG\}$ pairs, including
R1-antisense traces.

---

## Part D — Implementation phases

### Phase 1 — Native: strand-carrying regions, orientation arrays, and one helper

Single PR. After this lands, the scanner emits orientation-resolved
arrays; Python still uses the legacy collapsed arrays.

**Files**
- `src/rigel/native/calibration/orient.h` *(new)*
- `src/rigel/native/calibration/region_index.h`
- `src/rigel/native/calibration/accumulator.h`
- `src/rigel/native/calibration/accumulator.cpp`
- `src/rigel/native/bam_scanner.cpp`
- `src/rigel/calibration/_orient.py` *(new)*
- `src/rigel/calibration/scan_payload.py`
- `src/rigel/pipeline.py` (`_wire_calibration_regions`)
- `tests/test_calibration_accumulator.py`
- `tests/test_orient_routing.py` *(new)*

**Changes**

1. New `orient.h` exposing `classify_orient(region_strand,
   fragment_strand)` and the `ORIENT_*` constants. Used at every
   observation site.
2. `RegionIndex::set(...)` gains `const uint8_t* strands`; stores
   `std::vector<uint8_t> strands_`; exposes `strand(rid)`.
3. `BamScanner.set_regions(...)` Python binding accepts a `uint8`
   strand array. `_wire_calibration_regions` passes
   `region_df["strand"].to_numpy(np.uint8)` through the same lexsort
   permutation.
4. `CalibrationAccumulator::observe(...)` gains an `int8_t
   fragment_strand` parameter. Caller passes `result.exon_strand`.
   The accumulator never derives the fragment strand.
5. `CalibrationPayload` gains:
   ```cpp
   // Strand-orientation refinement of the legacy arrays.
   // Sums across the orient axis equal the corresponding legacy array
   // bit-exactly. Worker merge is element-wise additive.
   std::vector<int64_t> per_region_counts_by_orient;  // (R, mask::N_STATES, ORIENT_N)
   std::vector<int64_t> u_left_by_orient;             // (R, ORIENT_N)
   std::vector<int64_t> u_right_by_orient;            // (R, ORIENT_N)
   std::vector<int64_t> fl_hist_by_orient;            // (mask::N_STATES, ORIENT_N, kFlBins)
   ```
   No intergenic POS/NEG counter (Issue 5).
6. `observe(...)` increments legacy arrays exactly as today, then
   calls `classify_orient(...)` once and increments the new arrays.
7. Worker merge adds new arrays element-wise.
8. `CalibrationScanPayload.from_scan_dict(...)` validates the new
   arrays and asserts the cross-orient sum equals each legacy array
   element-wise. **This is the key correctness invariant.**

**Memory note.** Default dense layout costs ~192 B/region. If a
deployment has >5M regions and that becomes a problem, the dense
layout is restricted to mask cells that contribute to gDNA channels
(`MASK_INTERGENIC`, `MASK_INTRON`, `MASK_INTRON|MASK_EXON`), dropping
the factor by ~3×. Not done by default.

**Tests**
- Routing matrix using `classify_orient` (Python and C++).
- Strand round-trip into `RegionIndex`.
- R1-antisense / R1-sense symmetry on a mirrored synthetic trace.
- Single-thread vs multi-thread byte equality on every new array.
- Cross-orient sum equals legacy per-cell.

**Exit:** scanner emits orientation-resolved arrays; legacy outputs
unchanged; `compute_global_densities` unchanged.

---

### Phase 2 — Python: unified BLUE estimator and helper collapse

Single Python PR. After this lands, calibration uses BLUE-combined
densities when strand info exists, and degrades to legacy when it
does not.

**Files**
- `src/rigel/calibration/density_global.py` (rewritten around BLUE)
- `src/rigel/calibration/_orchestrator.py`
- `src/rigel/calibration/_result.py`
- `src/rigel/pipeline.py` (thread `StrandSummary` into `calibrate`)
- `tests/test_density_global.py`
- `tests/test_density_global_blue.py` *(new)*
- `tests/test_calibrate_orchestrator.py`

**Refactor — single density helper**

```python
@dataclass(frozen=True, slots=True)
class _ChannelObs:
    """Per-region observations for ONE density channel.

    A channel is an unbiased estimator of `rho_g` (fragments / bp) on
    a defined subset of regions, computed from a defined subset of
    fragments. Two channels in the same density type combine via
    inverse-variance weighting (BLUE).
    """
    label:          str
    rho_hat:        float          # point estimate
    rho_var:        float          # variance of point estimate (inf if no info)
    n_fragments:    float          # numerator (may be fractional for strand)
    eff_length_bp:  float          # denominator
    n_regions_used: int

@dataclass(frozen=True, slots=True)
class _BlueResult:
    rho_hat:   float
    rho_var:   float
    weights:   dict[str, float]    # channel label -> normalized weight
    channels:  tuple[_ChannelObs, ...]

def blue_combine(channels: Sequence[_ChannelObs]) -> _BlueResult: ...
```

`blue_combine` ignores channels with `w = 0` (i.e. infinite or
non-finite variance, or zero exposure). If every channel has zero
weight, returns `rho = 0, rho_var = inf` and a single
`weights = {label: 0.0}` entry — the caller can detect this and emit
a zero density.

**Channels per density type:**

| Density       | Channel B (boundary)                             | Channel S (strand)                                                   |
|---------------|--------------------------------------------------|-----------------------------------------------------------------------|
| INTERGENIC    | contained-fragment ratio (legacy)                | none (regions are `NONE`)                                             |
| INTRON        | none                                             | strand-deconvolved contained-fragment ratio on `POS`/`NEG` regions    |
| EXON-INTRON   | boundary-flux extrapolation (legacy)             | strand-deconvolved contained-fragment ratio on `POS`/`NEG` exon regions |

EXON-INTRON is the case where BLUE genuinely combines two channels.

**`_density(region_df, payload, type, …, strand: StrandSummary)`**

A single function builds the appropriate channels for the requested
type and calls `blue_combine`. Returns the new `GlobalGdnaDensity`:

```python
@dataclass(frozen=True, slots=True)
class GlobalGdnaDensity:
    type: DensityType
    rho: float
    rho_var: float                       # NEW
    n_fragments: int                     # legacy integer numerator (unchanged)
    n_fragments_estimated: float         # NEW: rho * eff_length_bp
    eff_length_bp: float
    n_regions_used: int
    kappa: KappaEstimate
    channels: tuple[_ChannelObs, ...]    # NEW: full diagnostic record
    weights: dict[str, float]            # NEW
```

`n_fragments` continues to be the **legacy integer numerator** so
existing summary consumers keep working. The new
`n_fragments_estimated` reflects the BLUE-combined point estimate.

**`StrandSummary` plumbed through `calibrate(...)`:**

```python
@dataclass(frozen=True, slots=True)
class StrandSummary:
    p_r1_sense: float          # finalized
    n_observations: int
    @property
    def c(self) -> float:      # 2p - 1, signed
        return 2.0 * self.p_r1_sense - 1.0
    @classmethod
    def from_model(cls, m: StrandModel) -> "StrandSummary": ...
    @classmethod
    def uninformative(cls) -> "StrandSummary":
        return cls(p_r1_sense=0.5, n_observations=0)
```

No `posterior_variance` field on `StrandSummary` — Issue 2 confirms it
is negligible. If we ever need it back, add it without changing the
estimator.

**Numerical guard.** Inside the strand-channel constructor, if
`abs(c) < 1e-3` the channel returns `_ChannelObs(rho_var = inf)`,
which BLUE will weight to zero. This is the only "magic number" in
the file; it is a finite-precision guard (well below the value of
`c` for any real stranded library; far above floating-point epsilon
for the `1/c` and `1/c²` terms).

**Tests**

- `p_r1_sense = 0.5` ⇒ every strand channel has `w = 0` ⇒
  every density falls back bit-equivalent to the pre-PR golden.
- `p_r1_sense = 1.0`, pure RNA in introns ⇒ INTRON `rho ≈ 0` with
  finite `rho_var`.
- `p_r1_sense = 1.0`, pure gDNA in introns ⇒ INTRON `rho` matches
  legacy within Poisson noise.
- EXON-INTRON BLUE combination on synthetic data with known truth:
  - boundary-only: `w_S = 0`, output equals legacy boundary estimator.
  - strand-only (eligibility flag stripped from rows): `w_B = 0`,
    output equals strand-only estimator.
  - both: variance shrinks vs either alone; mean within tolerance.
- R1-antisense mirror test: replace `(p, n_same, n_opp)` with
  `(1-p, n_opp, n_same)` and assert identical `rho_hat`, `rho_var`.
- Per-region clipping forbidden: assert internal raw-moment sum
  preserves negative per-region contributions.
- Cross-orient sum equality (already in Phase 1) re-checked
  end-to-end.

**Scenario tests**

- Remove xfail on `tests/scenarios/test_nrna_double_counting.py` for
  `ss = 1.0` and `ss = 0.9`. Record the `ss = 0.65` outcome.
- Implicit-splice and `nrna_none` golden gates must not regress.

**Exit:** strand-aware densities live; legacy preserved when no
strand info; xfailed scenarios cleared.

---

### Phase 3 — Diagnostics in `summary.json` + locoregional pass-through

Smallest PR. Surfaces information; does not change behavior.

**Files**
- `src/rigel/calibration/_diagnostics.py`
- `src/rigel/calibration/locus_prior.py` (read-only diagnostic fields)
- `tests/test_per_locus_gdna_mass.py`

**`summary.json` schema additions (additive; no removed keys):**

```
strand_aware:           bool
p_r1_sense:             float
strand_specificity:     float
strand_n_observations:  int

INTRON.rho_var:                  float
INTRON.weights:                  {strand: 1.0}    # one channel for INTRON
INTRON.n_strand_informative_regions: int

EXON-INTRON.rho_var:             float
EXON-INTRON.weights:             {boundary: w_B, strand: w_S}
EXON-INTRON.n_strand_informative_regions: int
```

Locoregional `n_gdna_intron` and `n_gdna_boundary_observed` get
parallel `_strand_aware` diagnostic fields when oriented data is
available. `expected_gdna_count_global` is unchanged.

**Exit:** users can see how much information the strand channel
contributed; team has the data to decide on prior surgery as a
follow-up.

---

### Phase 4 — Opt-in strand-deconvolved gDNA FL (diagnostic)

Optional small PR. Production default unchanged.

**Files**
- `src/rigel/calibration/_fl_sources.py`
- `src/rigel/calibration/fl.py`
- `src/rigel/config.py` (new flag)
- `tests/test_fl_models.py`

**Changes**
- New flag `EMConfig.strand_aware_gdna_fl: bool = False`.
- When true, `extract_gdna_counts(payload)` returns the
  per-FL-bin strand-deconvolved gDNA count
  $G_{\text{bin}} = \sum_\text{informative cells} (T_\text{bin} -
  D_\text{bin}/c)$ from `payload.fl_hist_by_orient`, with the same
  end-of-pipeline `max(0, …)` clip applied per bin after summing
  across mask cells. Strand-uninformative mask cells (e.g.
  `MASK_INTERGENIC`) are added in raw.
- When false, behavior is exactly today's.
- The existing EB Dirichlet smoothing in `build_fl_models` handles
  the sparser per-bin counts; no new smoothing primitive.

**Tests**
- Flag off ⇒ all existing FL tests pass unchanged.
- Flag on, `p = 1.0` ⇒ deconvolution reduces per-bin to
  `2 * minor_orientation`.
- Synthetic intron-nRNA trace: deconvolved gDNA FL mean closer to
  truth than legacy.

**Exit:** Issue 3 has a tested, documented opt-in path. If a real
data benchmark later shows the contamination bias matters in
practice, the default flips; until then, production stays on the
simpler path.

---

## Part E — Validation strategy (unchanged from v2 except as noted)

### Unit tests

- Routing-matrix tests (Phase 1) using `classify_orient`.
- Cross-orient-sum equals legacy invariant baked into payload validation.
- Worker-merge byte equality including new arrays.
- Phase 2 BLUE behavior at extremes (`w_B = 0`, `w_S = 0`, both
  positive) and on the R1-antisense mirror.
- Per-region clipping forbidden.
- Phase 4 FL deconvolution-at-`p=1` closed form.

### Scenario tests

- `tests/scenarios/test_nrna_double_counting.py`: remove xfail at
  `ss = 0.9` and `ss = 1.0` after Phase 2.
- Implicit-splice and clean `nrna_none` golden gates do not regress.

### Synthetic / benchmark

`scripts/synthetic_sim_sweep.py` configs:

- Stranded, gDNA-poor, nRNA-rich: INTRON `rho` and EXON-INTRON `rho`
  drop toward truth; EXON-INTRON `rho_var` drops vs Phase-1-only
  baseline (BLUE giving credit to the strand channel).
- Stranded, gDNA-rich, nRNA-none: no regression vs golden.
- Unstranded: bit-equivalent to golden.

Armis2 full-scale (`scripts/benchmarking/configs/default.yaml`):

- `gdna_*_ss_0.50_nrna_*` unchanged.
- `gdna_*_ss_1.00_nrna_rand` is the primary target — expected
  reduction in global gDNA over-allocation and improved mRNA
  recovery in the consuming EM.

### Acceptance gates

- Phase 1 may merge as soon as routing tests and byte-equality tests
  pass (no behavior change).
- Phase 2 must merge with the scenario xfail change.
- Phase 3 may merge any time after Phase 2.
- Phase 4 is opt-in; merges when Phase 4 tests are green.

---

## Part F — Risks and mitigations

| Risk | v3 mitigation |
|------|---------------|
| Two strand-comparison sites silently diverge | One canonical helper (`classify_orient`) in C++ and Python; routing matrix CI covers all pairs. |
| Inferring fragment strand from `obs_exons[0]` | `observe(...)` requires `fragment_strand` parameter; never derived in C++. |
| Per-region clip bias | One end-of-pipeline `max(0, ratio)`; explicit regression test forbids per-region clipping. |
| WLS variance under nRNA-rich conditions too tight | Variance computed from observed-count plug-in mixture, not gDNA-only collapse. |
| Magic-number cutoffs for channel selection | Replaced by inverse-variance weighting; only remaining number is `c_min = 1e-3` numerical guard. |
| Strand-model variance term mis-scaled | Dropped from estimator (Issue 2). API has no slot for it. |
| Contaminated gDNA FL bias | Opt-in Phase 4; orientation histogram already produced in Phase 1 so no second native PR is needed if the flag flips. |
| Locus prior surgery leaking into this work | Phase 3 is diagnostic only; prior surgery is a separate plan. |
| Edge correlation between boundary and strand channels on tiny exons | Acknowledged; small mass share in production; revisit only if benchmark CIs are systematically optimistic. |
| Mega-payload memory for huge genomes | Default dense `(R, 8, 3)` ~192 B/region; documented escape hatch to restrict to gDNA-relevant mask cells. |
| Removing intergenic POS/NEG counter could hide a future bug | Not a real risk: gDNA orientation is 50/50 by molecular biology. The data cannot violate it without an upstream bug that would surface in many other places first. |

---

## Part G — Out of scope

- Profiled-Poisson exact estimator with per-region RNA nuisance.
  BLUE on the moment-ratio is the v3 production path; promoting to a
  joint Poisson MLE is a follow-up if benchmarks demand it.
- Beta-binomial orientation overdispersion. `phi_orient = 1` until QC
  shows otherwise.
- Locus-prior surgery using strand-aware evidence — separate plan.
- EM-stage strand model changes — calibration-only plan.
- Strand-aware exon-only "gDNA detection" classifier (distinct from
  density estimation) — useful future work, not v3.
