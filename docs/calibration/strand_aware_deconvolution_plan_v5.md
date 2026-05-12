# Strand-Aware SRD Calibration — v5 Plan

Status: design, 2026-05-12. Supersedes
[`strand_aware_deconvolution_plan_v4.md`](strand_aware_deconvolution_plan_v4.md).

Implementation note: the final strand-mode decision rule is published in
[`strand_mode_decision_2026-05-12.md`](strand_mode_decision_2026-05-12.md).
That note supersedes the earlier public `min_strand_contrast` framing. The CLI
parameter was removed; `1e-3` remains only as an internal numerical floor.

v5 is a small refinement of v4. It accepts the v4 critique of v3 in full and
fixes two things v4 left underdone:

1. The framing "unstranded vs strand-corrected" must not imply that unstranded
  data are deprecated. They are not. On contaminated stranded channels, the
  uncorrected count/exposure estimator can be biased by RNA; the
  strand-deconvolved estimator is the unbiased version of the same calculation
  on the same physical data. v5 makes this explicit.
2. v4 declines to consolidate `_density_simple` and `_density_exon_intron`
   into a single helper. v5 does the consolidation, because the strand
   work is the natural forcing function and doing the same change twice is
   a code-health regression.

Everything else in v4 is preserved unchanged: minimal `(R, 3)` orientation
arrays, no production BLUE, no observed-count weights, no `fl_hist_by_orient`
in Phase 1, exon-contained strand contrast deferred, kappa deferred, opt-in
strand-aware FL flag (if ever added) lives in `CalibrationConfig`.

Review note: the v5 framing is implementation-ready only after one important
clarification. Strand-uninformative rows (`NONE`, `AMBIG`, or unknown fragment
strand) must not be added as raw RNA+gDNA counts into a strand-corrected
production density when identifiable exposure exists. Doing so would reintroduce
RNA contamination bias through the ambiguous subset. v5 therefore uses
one estimator over the **identifiable exposure subset** when strand information is
available, and uses the unstranded all-eligible ratio only when a channel is
globally unidentifiable or the strand summary is not statistically identifiable.

Finalization note: eight post-review polish items are folded in below. They do
not block implementation, but they lock the audit semantics before coding:
`n_fragments_estimated` projects onto all eligible exposure, each channel reports
its own `strand_active` flag, `n_rows_eligible` is preserved beside
`n_regions_used`, kappa remains a collapsed-count diagnostic, and a 99%
strand-model uncertainty margin prevents unstranded libraries from activating
the signed correction on sampling noise.

---

## Part A — Position relative to v3 and v4

| Question | v3 | v4 | v5 |
|----------|----|----|----|
| Deconvolve INTRON contained counts? | Yes | Yes | Yes |
| Deconvolve EXON-INTRON boundary flux? | Yes | Yes | Yes |
| Use exon-contained counts for global gDNA density? | Yes (BLUE) | No | No |
| Combine boundary + strand BLUE on EXON-INTRON? | Yes | No | No |
| Production observed-count variance weights? | Yes | No | No |
| Payload shape | `(R, 8, 3)` | `(R, 3)` per channel | `(R, 3)` per channel |
| `fl_hist_by_orient` in Phase 1? | Yes | No | No |
| Strand-aware FL flag location | `EMConfig` | `CalibrationConfig` | `CalibrationConfig` |
| Collapse `_density_simple` + `_density_exon_intron`? | Yes | No | **Yes** |
| Frame unstranded vs strand-corrected as two estimators? | n/a | Yes (fallback) | **No (one estimator, two regimes)** |

The two v5-only changes are the last two rows.

---

## Part B — Why "unstranded vs strand-corrected" is the right framing

The existing `INTRON` estimator is

$$\hat\rho_\text{intron} \;=\; \frac{\sum_i (n_{\text{same},i} + n_{\text{opp},i})}{\sum_i L^\text{eff}_i} \;=\; \frac{\sum_i T_i}{\sum_i E_i}.$$

Under contamination, $\mathbb E[T_i] = N_{r,i} + N_{g,i}$, so this estimator
returns "RNA + gDNA per bp", not "gDNA per bp". It is the right answer only
when $N_r = 0$ on intron rows — which is exactly the assumption that breaks
with stranded nRNA.

The strand-deconvolved estimator on the *same data* is

$$\hat\rho_\text{intron} \;=\; \frac{\sum_i (T_i - D_i / \texttt{signed\_strand\_contrast})}{\sum_i E_i}.$$

Under the same observation model,
$\mathbb E[T_i - D_i/\texttt{signed\_strand\_contrast}] = N_{g,i}$, so this
estimator returns "gDNA per bp" — the quantity we actually want — without any
assumption on $N_r$.

The two formulas differ by exactly
$D_i/\texttt{signed\_strand\_contrast}$. As `signed_strand_contrast` approaches
zero (unstranded), the correction blows up because the data carries no
information about the split, and the only sensible thing is to set the
correction to zero, which gives the unstranded count/exposure formula. As the absolute
contrast approaches 1 (perfectly stranded), the correction becomes exact.

So there is one estimator with a strand-correction term that the data turns
on or off. There are not two competing estimators. v5 codes it that way:

```python
def _gdna_count_moment(
    n_same: np.ndarray,        # informative rows only
    n_opp:  np.ndarray,
    *,
    signed_strand_contrast: float,
) -> np.ndarray:
    """Return per-row unbiased gDNA count moments.

    The caller uses this only after the absolute signed strand contrast clears
    both the internal numerical floor and the strand-model uncertainty margin.
    Aggregation and the final non-negativity clip happen in the caller.
    """
    T = n_same + n_opp
    D = (n_same.astype(np.float64) - n_opp.astype(np.float64))
    return T.astype(np.float64) - D / signed_strand_contrast
```

The unstranded path collapses to "set the signed correction term to zero when a
whole channel is statistically or numerically unidentifiable". It does
**not** mean that `UNINF` rows should be added as raw counts to a corrected
stranded estimate. Those raw counts are still RNA+gDNA and can be biased in
exactly the same way as the uncorrected count/exposure estimator.

The only remaining decision is which exposure is identifiable. A row or boundary
side is identifiable iff `region_strand ∈ {POS, NEG}`. The observed fragments on
that identifiable exposure enter `n_same`, `n_opp`, or `n_uninf`; only the
`SAME/OPP` observations carry strand contrast. `UNINF` observations are reported
as diagnostics and included in the unstranded estimator, but are not added to the
corrected production numerator.

When a channel has identifiable exposure and
`abs(signed_strand_contrast) >= max(STRAND_CONTRAST_NUMERICAL_FLOOR, uncertainty_margin_99)`,
estimate density on that identifiable subset:

```python
G_raw = _gdna_count_moment(n_same, n_opp, signed_strand_contrast).sum()
E_raw = E_identifiable.sum()
rho   = max(0.0, G_raw / E_raw)
```

When a channel has no identifiable exposure, or the absolute signed contrast
does not clear that effective guard, the correction is unavailable and the
estimator is the unstranded all-eligible ratio. This is still one estimator with two
identifiability regimes, not two competing estimators.

Small residual bias: on identifiable rows, observations with unknown fragment
strand (`n_uninf`) are excluded from the corrected numerator while the density is
estimated over the full identifiable exposure. If unknown-strand status is
independent of RNA-vs-gDNA state, this slightly under-counts the gDNA share of
`n_uninf`. In practice these counts should be rare; v5 reports them rather than
adding a noisy small-count correction.

---

## Part C — Statistical core (unchanged from v4)

For an informative row $i$:

```text
Y_same_i, Y_opp_i  = orientation-classified counts
T_i = Y_same_i + Y_opp_i
D_i = Y_same_i - Y_opp_i
p   = StrandModel.p_r1_sense  (signed; R1-sense > 0.5, R1-antisense < 0.5)
signed_strand_contrast = 2p - 1
```

Observation model:

```text
E[Y_same_i] = p       * N_r,i + 0.5 * N_g,i
E[Y_opp_i]  = (1 - p) * N_r,i + 0.5 * N_g,i
E[D_i]      = signed_strand_contrast * N_r,i
E[T_i]      = N_r,i + N_g,i
=>   E[T_i - D_i / signed_strand_contrast] = N_g,i
```

Estimator (per density channel):

```text
effective_min = max(STRAND_CONTRAST_NUMERICAL_FLOOR, uncertainty_margin_99)
if identifiable exposure exists and abs(signed_strand_contrast) >= effective_min:
  G_raw = sum_i (T_i - D_i / signed_strand_contrast) over identifiable rows
  E_raw = sum_i E_i            over identifiable rows
else:
  G_raw = sum_i T_i            over all eligible rows
  E_raw = sum_i E_i            over all eligible rows
rho = max(0, G_raw / E_raw)
```

Properties carried from v4 unchanged:

- Per-row signed-strand moments may be negative or exceed $T_i$ under noise. The only
  non-negativity clip is at the final `max(0, ratio)`.
- `STRAND_CONTRAST_NUMERICAL_FLOOR = 1e-3` is an internal numerical floor; the
  effective production guard is the maximum of that floor and a 99% uncertainty
  margin for the signed strand contrast. There is no public CLI/YAML parameter
  for this decision.
- `RegionStrand.AMBIG`, `RegionStrand.NONE`, and unknown fragment strand all
  route to `UNINF`. In strand-active mode, `UNINF` counts are diagnostics;
  they are not used as raw gDNA counts in the corrected numerator.
- Uncertainty is reported through aggregate diagnostics such as the informative
  exposure fraction and UNINF count fraction, and it also gates whether the
  signed correction is identifiable. It does not weight rows in the production
  point estimate.
- Small nonzero absolute contrast values make the signed correction unbiased in
  expectation but very noisy in finite samples; if the 99% interval contains
  zero signed contrast, the unstranded all-eligible ratio is the identifiable regime.

---

## Part D — Density channels (unchanged from v4 in scope)

| Channel        | Rows                                       | Numerator (v5 formula)                                    | Exposure                            |
|----------------|--------------------------------------------|-----------------------------------------------------------|-------------------------------------|
| INTERGENIC     | `type == INTERGENIC`                       | `n_uninf` only (regions are `NONE`, no `n_same/n_opp`)    | `L_eff_contained`                   |
| INTRON         | `type == INTRON`, strand-identifiable subset when active | identifiable: `T - D / signed_strand_contrast`; unstranded: `T` over all eligible rows | `L_eff_contained`                   |
| EXON-INTRON    | eligible EXON boundary sides, strand-identifiable subset when active | identifiable side: `T - D / signed_strand_contrast`; unstranded: raw boundary flux | `B_cross * sides_eligible`         |

INTERGENIC is the degenerate case: there are no `POS/NEG` intergenic regions,
so the channel has no identifiable exposure and uses the unstranded
count/exposure ratio.

INTRON and EXON-INTRON are where the deconvolution actually does work. Both
use the *same* numerator helper. EXON-INTRON applies it once per eligible
boundary side rather than once per region, which is exactly how the collapsed
boundary flux is summed today.

When strand-active, `rho` is estimated on the identifiable POS/NEG subset and
then projected as the channel-wide density for all eligible rows of that type.
This is the same uniform-density assumption the unstranded estimator already makes,
but now estimated from the subset that can separate RNA from gDNA.

Exon-contained strand contrast is **not** in v5 production, for the reasons
v4 lays out (sensitivity to small `p_r1_sense` errors at high mRNA counts,
antisense annotation overlap, etc.). It can be a future diagnostic.

---

## Part E — One density helper (the v5 code-health change)

`density_global.py` currently has two near-duplicate functions:

```python
def _density_simple(*, type_label, region_mask, counts_col, leff): ...
def _density_exon_intron(*, region_mask, bf_left, bf_right,
                         u_left, u_right, b_cross): ...
```

Both build a numerator and a per-row exposure, mask to eligible rows, sum,
divide, and call `estimate_kappa`. The only differences are:

- INTERGENIC/INTRON have one numerator column; EXON-INTRON has two
  (left and right side numerators) and a per-row eligibility mask from the
  boundary flags.
- INTERGENIC/INTRON exposure is `L_eff`; EXON-INTRON exposure is
  `B_cross * sides`.

v5 unifies them into one helper that operates on a small struct of per-row
arrays:

```python
@dataclass(frozen=True, slots=True)
class _ChannelRows:
    """Per-row inputs for a single density channel.

    Each array has length R (the full region count). Rows that are not
    eligible carry zeros and do not affect the ratio.
    """
    label:    DensityType
    eligible: np.ndarray   # (R,) bool: row contributes to unstranded denominator
    n_same:   np.ndarray   # (R,) int64: strand-SAME counts
    n_opp:    np.ndarray   # (R,) int64: strand-OPP counts
    n_uninf:  np.ndarray   # (R,) int64: observed counts without usable orientation
    exposure: np.ndarray   # (R,) float64: per-row effective length / boundary exposure
    identifiable: np.ndarray  # (R,) bool: row/side has POS or NEG region strand


def _compute_density(
    rows: _ChannelRows,
    *,
  strand_summary: StrandSummary,
) -> GlobalGdnaDensity:
    elig = rows.eligible
    ident = elig & rows.identifiable & (rows.exposure > 0)
  uncorrected_counts_all = rows.n_same + rows.n_opp + rows.n_uninf
  total_exposure = float(rows.exposure[elig].sum())
  n_observed = float(uncorrected_counts_all[elig].sum())
  rho_uncorrected = n_observed / total_exposure if total_exposure > 0.0 else 0.0
  signed_strand_contrast = strand_summary.signed_strand_contrast
  effective_min = max(
    STRAND_CONTRAST_NUMERICAL_FLOOR,
    strand_summary.signed_strand_contrast_margin(confidence=0.99),
  )
    use_strand = (
    abs(signed_strand_contrast) >= effective_min
        and bool(np.any(ident))
    )
  # Point estimate: corrected on identifiable exposure only; otherwise unstranded.
    if use_strand:
        D = rows.n_same.astype(np.float64) - rows.n_opp.astype(np.float64)
        T = (rows.n_same + rows.n_opp).astype(np.float64)
        corrected = np.where(ident, T - D / signed_strand_contrast, 0.0)
        G_raw = float(corrected.sum())
        E_raw = float(rows.exposure[ident].sum())
        n_used = int(np.count_nonzero(ident))
    else:
        G_raw = n_observed
        E_raw = total_exposure
        n_used = int(np.count_nonzero(elig & (rows.exposure > 0)))

    rho = max(0.0, G_raw / E_raw) if E_raw > 0 else 0.0
      n_fragments_estimated = float(rho * total_exposure)
      strand_correction_fragments = n_fragments_estimated - n_observed
      # Kappa: keep using collapsed counts as in v4; strand-aware
    # kappa for signed moments is deferred. estimate_kappa is called against
      # the internally consistent count/rho_uncorrected pair, not the corrected rho.
      uncorrected_counts = uncorrected_counts_all[elig].astype(np.int64)
      kappa = estimate_kappa(uncorrected_counts, rows.exposure[elig], rho_uncorrected)
    return GlobalGdnaDensity(
        type=rows.label, rho=rho,
        n_fragments=int(uncorrected_counts.sum()),  # raw integer numerator
        n_fragments_estimated=n_fragments_estimated,  # projected channel total
        eff_length_bp=E_raw, n_regions_used=n_used, kappa=kappa,
        n_rows_eligible=int(np.count_nonzero(elig & (rows.exposure > 0))),
        strand_active=use_strand,
        rho_uncorrected=rho_uncorrected,
        strand_correction_fragments=float(strand_correction_fragments),
        n_strand_informative_regions=int(np.count_nonzero(ident)),
        strand_informative_exposure_fraction=(
          float(rows.exposure[ident].sum() / total_exposure)
          if total_exposure > 0.0 else 0.0
        ),
        n_uninf_observed=int(rows.n_uninf[elig].sum()),
    )
```

Three thin builders construct `_ChannelRows` for each density type. They
hold no statistics — just data shaping:

```python
def _channel_intergenic(region_df, payload) -> _ChannelRows: ...
def _channel_intron    (region_df, payload) -> _ChannelRows: ...
def _channel_exon_intron(region_df, payload, *, b_cross) -> _ChannelRows: ...
```

`compute_global_densities(...)` becomes a 10-line orchestrator that builds
the three channels and calls `_compute_density` three times.

This is a net reduction in `density_global.py`'s size and complexity, in
addition to enabling the strand correction. It is the right time to do this
refactor because the alternative — leaving two near-identical functions and
adding strand handling to each — duplicates the same code change.

---

## Part F — Convention: fragment strand vs region transcript strand

Verbatim from v4 (and v3); this paragraph is the authoritative source.

- `region_strand ∈ {NONE = 0, POS = 1, NEG = 2, AMBIG = 3}` matches
  `rigel.calibration.regions.RegionStrand`.
- `fragment_strand ∈ {NONE = 0, POS = 1, NEG = 2}` matches `rigel.types.Strand`
  and is the post-R2-flip exon-block strand. It is exactly `result.exon_strand`
  from the resolver. The accumulator never derives it from
  `obs_exons[0].strand` or from BAM flags.
- `p_r1_sense = P(fragment_strand == region_strand | RNA, region unambiguous)`,
  taken from `StrandModel.p_r1_sense` after finalize. R1-sense > 0.5;
  R1-antisense < 0.5. `signed_strand_contrast = 2p - 1` is signed; we never
  invert `p` anywhere.
- A single helper `classify_orient(region_strand, fragment_strand)` exists in
  C++ (`src/rigel/native/calibration/orient.h`) and Python
  (`src/rigel/calibration/_orient.py`) with identical numeric values
  `ORIENT_SAME = 0`, `ORIENT_OPP = 1`, `ORIENT_UNINF = 2`. Every routing site
  calls one of these helpers; bare `0/1/2` literals are forbidden at call
  sites.

---

## Part G — Implementation phases

Three phases. Phase 1 lands native plumbing; Phase 2 lands the Python change
and the helper consolidation together; Phase 3 surfaces diagnostics.

### Phase 1 — Native: orientation routing, minimal arrays, single helper

Behavior change: none.

Files:
- `src/rigel/native/calibration/orient.h` *(new)*
- `src/rigel/native/calibration/region_index.h`
- `src/rigel/native/calibration/accumulator.h`
- `src/rigel/native/calibration/accumulator.cpp`
- `src/rigel/native/bam_scanner.cpp`
- `src/rigel/calibration/_orient.py` *(new)*
- `src/rigel/calibration/scan_payload.py`
- `src/rigel/pipeline.py` (`_wire_calibration_regions`)
- `tests/test_orient_routing.py` *(new)*
- `tests/test_calibration_accumulator.py`

Changes:

1. `orient.h` defines `ORIENT_SAME = 0`, `ORIENT_OPP = 1`, `ORIENT_UNINF = 2`,
   `ORIENT_N = 3`, and a single inline `classify_orient(region_strand,
   fragment_strand)`.
2. `RegionIndex::set(...)` takes `const uint8_t* strands`; stores
   `strands_`; exposes `strand(rid)`.
3. `BamScanner.set_regions(...)` Python binding accepts a `uint8` strand
   array. `_wire_calibration_regions` passes
   `region_df["strand"].to_numpy(np.uint8)` through the lexsort permutation.
4. `CalibrationAccumulator::observe(...)` gains `int8_t fragment_strand`.
   The caller passes `result.exon_strand` and only that. The accumulator
   never derives the fragment strand.
5. Add three minimal arrays to `CalibrationPayload`:
   ```cpp
   std::vector<int64_t> intron_counts_by_orient;  // (R, ORIENT_N), populated only on MASK_INTRON observations
   std::vector<int64_t> u_left_by_orient;         // (R, ORIENT_N)
   std::vector<int64_t> u_right_by_orient;        // (R, ORIENT_N)
   ```
   Memory: ~72 B/region for the three arrays combined (vs ~192 B/region in
   v3's dense layout).
6. Increment collapsed arrays exactly as today, then route the matching new
   array slot via `classify_orient`. `u_left/u_right` boundary increments
   route by the EXON region's strand and the fragment strand.
7. Worker merge adds new arrays element-wise.
8. `CalibrationScanPayload.from_scan_dict(...)` validates and asserts the
   correctness invariants element-wise:
   ```text
   intron_counts_by_orient[r, :].sum() == per_region_counts[r, MASK_INTRON]
   u_left_by_orient[r, :].sum()        == u_left[r]
   u_right_by_orient[r, :].sum()       == u_right[r]
   ```
   These three asserts are the regression guard for routing bugs.

Tests:
- C++/Python routing matrix for all 4×3 strand pairs, including R1-antisense
  traces (verified by mirroring `(p, n_same, n_opp) ↔ (1-p, n_opp, n_same)`).
- Strand round-trip into `RegionIndex` after lexsort.
- Cross-orient sums equal collapsed arrays element-wise.
- Single-thread vs many-thread scan byte equality on every new array.

Exit: scan emits orientation arrays; collapsed outputs unchanged; Python
behavior unchanged.

### Phase 2 — Python: helper consolidation + strand-deconvolved densities

Behavior change: INTRON and EXON-INTRON densities use the unbiased
strand-corrected numerator on informative rows when the absolute signed strand
contrast clears both the internal numerical floor and the 99% strand-model
uncertainty margin.

Files:
- `src/rigel/calibration/density_global.py` (rewritten around `_compute_density`)
- `src/rigel/calibration/_orchestrator.py`
- `src/rigel/calibration/_result.py` (additive fields)
- `src/rigel/pipeline.py` (thread `StrandSummary` into `calibrate`)
- `tests/test_density_global.py`
- `tests/test_calibrate_orchestrator.py`

Changes:

1. Add `StrandSummary` to `src/rigel/calibration/_orient.py` beside the Python
  `classify_orient` helper. It carries the MLE plus enough observations to
  compute a normal-approximation standard error for signed strand contrast.
   ```python
   @dataclass(frozen=True, slots=True)
   class StrandSummary:
       p_r1_sense: float
       n_observations: int

       @property
       def signed_strand_contrast(self) -> float:
           return 2.0 * self.p_r1_sense - 1.0

         @property
         def signed_strand_contrast_se(self) -> float: ...

         def signed_strand_contrast_margin(self, *, confidence: float = 0.99) -> float: ...

       @classmethod
       def from_model(cls, m: StrandModel) -> "StrandSummary": ...

       @classmethod
       def uninformative(cls) -> "StrandSummary":
           return cls(p_r1_sense=0.5, n_observations=0)
   ```
2. Replace `_density_simple` and `_density_exon_intron` with `_ChannelRows`,
  three thin channel builders, and one `_compute_density(...)` that takes a
  `StrandSummary`. The pre-existing `compute_global_densities(...)` signature
  gains `strand_summary: StrandSummary | None`; it defaults to
  `StrandSummary.uninformative()`.
3. `GlobalGdnaDensity` schema additions (additive, no removed fields):
  ```python
  n_fragments_estimated: float = 0.0          # = rho * all eligible exposure
  n_rows_eligible: int = 0                    # all eligible rows before correction
  strand_active: bool = False                 # this channel used correction
  rho_uncorrected: float = 0.0                # unstranded all-eligible ratio
  strand_correction_fragments: float = 0.0    # estimated total minus raw observed
  n_strand_informative_regions: int = 0
  strand_informative_exposure_fraction: float = 0.0
  n_uninf_observed: int = 0
  ```
  `n_fragments` continues to be the raw integer numerator across all eligible
  rows so existing summaries keep working. `n_regions_used` continues
  to mean rows that actually entered the point estimate; `n_rows_eligible`
  preserves cross-mode comparability.
4. Pipeline threads a calibration-specific `StrandSummary` into `calibrate(...)`.
  The source is always `exonic_spliced`. If that RNA-pure model is not
  identifiable, calibration does not substitute the exonic diagnostic model;
  channels run unstranded and the run is warned for inspection.

Tests:
- `p_r1_sense = 0.5`, or near-0.5 estimates whose signed contrast is within the
  99% uncertainty margin, ⇒ correction is zero ⇒ output matches the unstranded
  all-eligible ratio on all three densities.
- Sparse-spliced toy libraries must not be rescued with exonic diagnostic
  strand signal. Tests should either provide enough spliced observations for the
  intended stranded assertion or expect unstranded calibration with a warning.
- `p_r1_sense = 1.0`, pure RNA in introns ⇒ `rho_intron ≈ 0`.
- `p_r1_sense = 1.0`, pure gDNA in introns ⇒ `rho_intron` matches the unstranded
  count/exposure density
  within Poisson noise.
- R1-antisense mirror: `(p, n_same, n_opp) → (1-p, n_opp, n_same)` produces
  identical `rho`.
- Many-small-region pure-gDNA simulation: aggregate raw moments give
  unbiased `rho_intron`. Per-region clipping regression test asserts
  negative per-row moments survive into the aggregate sum.
- Mixed informative/AMBIG intron fixture: in strand-active mode, AMBIG/UNINF
  raw counts do not enter the corrected production numerator; they are visible
  only in `n_uninf_observed`, `rho_uncorrected`, and unstranded-mode tests.
- Kappa consistency test: corrected `rho` and collapsed-count `kappa` are not mixed in
  the method-of-moments call; `estimate_kappa` sees raw counts with
  `rho_uncorrected` until signed-moment kappa is designed.
- Cross-orient sum equals collapsed-count invariant exercised end-to-end.
- `StrandSummary.uninformative()` path is bit-equivalent to passing no
  strand summary.

Scenario tests:
- Remove xfail on `tests/scenarios/test_nrna_double_counting.py` for
  `ss = 0.9` and `ss = 1.0`. Record `ss = 0.65` outcome but do not gate
  the PR on it.
- Implicit-splice and `nrna_none` golden gates must not regress.

Exit: INTRON and EXON-INTRON densities are unbiased on contaminated stranded
data; unstranded data is bit-equivalent; `density_global.py` is shorter and
has one density helper instead of two.

### Phase 3 — Diagnostics

Behavior change: none.

`summary.json` additive keys:
```text
strand_aware:               bool             # = correction active after identifiability gate
p_r1_sense:                 float
strand_specificity:         float
strand_n_observations:      int
INTRON.n_strand_informative_regions:        int
INTRON.n_rows_eligible:                     int
INTRON.strand_active:                       bool
INTRON.rho_uncorrected:                     float
INTRON.strand_correction_fragments:         float
INTRON.strand_informative_exposure_fraction: float
INTRON.n_uninf_observed:                    int
INTRON.n_fragments_estimated:               float
EXON-INTRON.n_strand_informative_regions:   int
EXON-INTRON.n_rows_eligible:                int
EXON-INTRON.strand_active:                  bool
EXON-INTRON.rho_uncorrected:                float
EXON-INTRON.strand_correction_fragments:    float
EXON-INTRON.strand_informative_exposure_fraction: float
EXON-INTRON.n_uninf_observed:               int
EXON-INTRON.n_fragments_estimated:          float
```

No locus-prior changes. `expected_gdna_count_global` continues to use the
global `rho` values; those are now unbiased on contaminated channels, which
is the whole point.

Known limitation: until signed-moment kappa is designed, `kappa` is estimated
from collapsed counts with `rho_uncorrected`, while downstream shrinkage
uses the corrected `rho`. This is internally consistent for the kappa diagnostic
itself but can slightly mis-calibrate shrinkage scale in channels where
deconvolution strongly lowers `rho`.

---

## Part H — Out of scope (consolidated)

- Exon-contained strand contrast as a production estimator (v3's BLUE
  channel). Reserved as a future experimental diagnostic; not required
  to fix the contamination problem this plan targets.
- Inverse-variance BLUE combiners with observed-count plug-in variances.
- Profiled-Poisson joint MLE with per-region RNA nuisance.
- Strand-aware kappa for signed/fractional moments.
- Strand-aware gDNA fragment-length model. Out of v5; if added later, the
  flag belongs in `CalibrationConfig`. The orientation FL histogram itself
  also stays out of Phase 1, because a fragment-level FL orientation
  histogram requires an explicit per-fragment agreement rule that v5 does
  not specify.
- Beta-binomial orientation overdispersion as a production correction.
- Locus-prior surgery using local strand evidence.
- EM-stage strand-model changes.

---

## Part I — Risks and mitigations

| Risk | v5 mitigation |
|------|---------------|
| Two strand-comparison sites diverge | One C++ helper, one Python helper, routing-matrix CI. |
| Fragment strand inferred incorrectly | `observe(...)` takes `result.exon_strand` only; never derived in C++. |
| Per-region clipping bias | Single end-of-pipeline `max(0, ratio)`; explicit regression test forbids per-region clipping. |
| Observed-variance zero-count weighting | No production variance weighting at all; `_compute_density` is a direct ratio. |
| Exon-contained mRNA poisons gDNA estimate | Exon-contained strand contrast is not in production. |
| Boundary nRNA contamination remains | Boundary flux numerator is orientation-deconvolved via `u_*_by_orient`. |
| Ambiguous/unknown-strand rows reintroduce RNA contamination bias | Strand-active densities use identifiable exposure only; UNINF counts are diagnostics/fallback data. |
| Unknown-strand counts on identifiable rows cause tiny under-count | `n_uninf_observed` is reported; no noisy multiplicative correction in v5. |
| Identifiable subset may not represent all eligible rows | Documented uniform-density extrapolation; report informative exposure fraction. |
| Mode-flag combinatorics ("unstranded vs strand-corrected") | Eliminated. The estimator has one branch (`abs(signed_strand_contrast) >= max(STRAND_CONTRAST_NUMERICAL_FLOOR, uncertainty_margin_99)`) that toggles whether the correction term is computed; everything else is uniform. |
| Code-health regression from leaving two density helpers | Phase 2 collapses them into `_compute_density` while making the strand change. |
| Mega-payload memory | `(R, 3) × 3` arrays only (~72 B/region). |
| Strand-aware kappa undefined for signed moments | Phase 2 keeps kappa on collapsed counts; deferred. |
| Weakly stranded data amplifies orientation noise | Require the signed contrast to clear both the internal numerical floor and the 99% uncertainty margin; synthetic ss=0.50 benchmarks pin this fallback. |
| FL histogram double-counting | `fl_hist_by_orient` not introduced; future work needs an explicit fragment-level rule. |
| Strand-model variance term mis-scaled | Used only as an identifiability margin, not as a row weight or correction term. |

---

## Part J — Acceptance gates

- Phase 1 may merge once routing-matrix tests, byte-equality tests, and the
  three cross-orient-sum-equals-collapsed-count invariants pass. No behavior change.
- Phase 2 must merge with the scenario xfail removals at `ss = 0.9` and
  `ss = 1.0`, and with the helper consolidation in the same PR.
- Phase 3 may merge any time after Phase 2.

Benchmarks to run after Phase 2:

- `scripts/synthetic_sim_sweep.py` configs:
  - stranded gDNA-poor nRNA-rich: INTRON/EXON-INTRON `rho` move toward truth;
  - stranded gDNA-rich nRNA-none: no regression vs golden;
  - unstranded: bit-equivalent.
- Armis2 `scripts/benchmarking/configs/default.yaml`:
  - `gdna_*_ss_0.50_nrna_*` unchanged;
  - `gdna_*_ss_1.00_nrna_rand` is the primary target — expected to reduce
    global gDNA over-allocation and improve mRNA recovery downstream.
