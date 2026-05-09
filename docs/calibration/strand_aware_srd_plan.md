# Strand-Aware SRD Calibration Plan

Status: design revision, 2026-05-09.

This document replaces the rough `Wss`-weighted proposal. The central change is that
strand information should enter calibration through the observation likelihood and its
uncertainty, not through hard strand-specificity thresholds.

## 1. Goal

Rigel's current SRD calibration estimates gDNA density from geometric categories:

- `INTERGENIC`: contained fragments in intergenic regions.
- `INTRON`: contained fragments in intronic regions.
- `EXON-INTRON`: unspliced fragments crossing eligible exon/intron boundaries.

That works when the intron and boundary channels are mostly gDNA. It breaks when
nascent RNA is abundant: stranded nRNA produces exactly the same intron and boundary
geometry as gDNA, so the geometry-only numerator overestimates gDNA and inflates the
global gDNA prior.

The enhancement is to keep the current geometry foundation, but split genic calibration
observations by orientation whenever the annotation region has a unique transcript
strand. In a stranded library, RNA is enriched in one orientation and gDNA remains
approximately 50/50. That lets SRD separate gDNA from nRNA in the exact channels that
are currently contaminated.

Non-goals for the first implementation:

- Do not require strand-specific data. Unstranded libraries must preserve the current
  calibration behavior.
- Do not add heuristic `SS` cutoffs or a hand-tuned trust ramp.
- Do not change the EM component model in this phase; this only improves the gDNA
  calibration priors feeding EM.

## 2. Current Code Facts

The existing code is already close to the right data model, but the scanner drops the
strand dimension before calibration reaches Python.

- `src/rigel/calibration/regions.py` already emits a `strand` column with
  `RegionStrand.NONE`, `POS`, `NEG`, and `AMBIG`.
- `src/rigel/pipeline.py::_wire_calibration_regions` passes region ids, starts, ends,
  and type masks into native code, but not `region_df["strand"]`.
- `src/rigel/native/calibration/region_index.h` stores only starts, ends, and type
  masks. It needs a `strands_` vector and `strand(rid)` accessor.
- `src/rigel/native/calibration/accumulator.cpp` emits collapsed `global_counts`,
  `per_region_counts`, `fl_hist`, `u_left`, and `u_right` arrays. There is no
  orientation split.
- `src/rigel/calibration/scan_payload.py` validates only the collapsed payload schema.
- `src/rigel/calibration/density_global.py` computes the three geometry-only densities.
- `src/rigel/strand_model.py` already exposes the needed library protocol quantities:
   `p_r1_sense`, `strand_specificity`, `read1_sense`, `n_observations`, and
   `posterior_variance()`.

Important correction: the strand-aware estimator must use `p_r1_sense`, not only
`strand_specificity = max(p_r1_sense, 1 - p_r1_sense)`. Some protocols are
R1-antisense, so the RNA-enriched bucket can be opposite the transcript strand.

## 3. Strand Model

For a region with a unique transcript strand, define two observed orientation counts:

- `Y_same`: fragment strand equals the region transcript strand.
- `Y_opp`: fragment strand is opposite the region transcript strand.

The fragment strand should use the resolver/scanner convention already used for strand
training: R1 on its reference strand, R2 flipped into R1-equivalent orientation. In
practice the calibration observation site should pass `result.exon_strand` or an
explicit fragment-strand value derived from the same logic into the accumulator. It
should not infer the fragment strand from `obs_exons[0].strand`; multi-block fragments
can make that shortcut ambiguous or wrong.

Let:

- `p = p_r1_sense = P(fragment strand == transcript strand | RNA)`.
- `q = 1 - p`.
- `c = p - q = 2p - 1`.
- `T = Y_same + Y_opp`.
- `D = Y_same - Y_opp`.
- `N_r` be the unknown RNA count in the calibration region or boundary side.
- `N_g` be the unknown gDNA count in the same unit.

The expected counts are:

```text
E[Y_same] = p * N_r + 0.5 * N_g
E[Y_opp]  = q * N_r + 0.5 * N_g
```

Therefore:

```text
E[D] = c * N_r
E[T] = N_r + N_g
```

When `c != 0`, the unconstrained moment estimator is:

```text
N_r_hat = D / c
N_g_hat = T - D / c
```

The implementation should use the constrained non-negative form:

```text
N_r_hat = clip(D / c, 0, T)
N_g_hat = T - N_r_hat
```

Equivalent major/minor form, useful for diagnostics:

```text
delta = abs(2p - 1)
Y_major = Y_same if p >= 0.5 else Y_opp
Y_minor = Y_opp  if p >= 0.5 else Y_same

N_r_hat = clip((Y_major - Y_minor) / delta, 0, T)
N_g_hat = T - N_r_hat
```

Boundary checks:

- If `p = 1.0`, then `N_g_hat = 2 * Y_opp`.
- If `p = 0.0`, then `N_g_hat = 2 * Y_same`.
- If `p = 0.5`, then `c = 0`; strand evidence is unidentifiable and contributes no
  precision.

This is the key reason no `Wss` ramp is needed. The math itself becomes uninformative
as the library approaches unstranded.

## 4. Uncertainty and Information

The strand estimator should report a count estimate and a variance. Downstream density
code can then use inverse-variance weighting or a profiled likelihood, both of which
naturally give zero weight to uninformative strand data.

For a single count unit, use a delta-method approximation around the observed or fitted
counts. With signed `c = 2p - 1`:

```text
N_g_hat = T - D / c
        = (1 - 1/c) * Y_same + (1 + 1/c) * Y_opp
```

Assuming independent Poisson orientation counts, the count-noise variance is:

```text
Var_count(N_g_hat) ~= (1 - 1/c)^2 * mu_same
                    + (1 + 1/c)^2 * mu_opp
```

where `mu_same` and `mu_opp` can be plug-in fitted values after constraining
`N_r_hat` and `N_g_hat`.

Strand-model uncertainty should be included when the model has few training reads. If
`Var(p)` comes from `StrandModel.posterior_variance()`, then:

```text
dN_g_hat/dp = 2D / c^2
Var_p(N_g_hat) ~= (2D / c^2)^2 * Var(p)
```

Total variance:

```text
Var(N_g_hat) = Var_count(N_g_hat) + Var_p(N_g_hat) + variance_floor
precision    = 1 / Var(N_g_hat)
```

The floor is only a numerical guard against divide-by-zero in degenerate tiny counts;
it should be fixed in units of fragments, not tuned to performance outcomes.

Consequences:

- As `p -> 0.5`, `c -> 0`, the variance diverges and precision goes to zero.
- Sparse regions get low precision even if `SS` is high.
- Highly stranded, high-count regions get high precision.
- A poorly trained strand model contributes less because `Var(p)` is larger.

This replaces the previous proposal's magic
`clamp((SS - 0.55) / 0.35, 0, 1)` weight.

## 5. Density Estimation Strategy

Do not treat the geometry-only total and the strand split as two independent pieces of
evidence for the same fragments. The strand split refines the same observations.

Recommended first implementation:

1. Keep the existing `INTERGENIC` density unchanged. Intergenic rows have no
   transcript strand and are already the cleanest geometry-only gDNA control.
2. For `INTRON`, compute strand-deconvolved gDNA counts only for intron regions whose
   `RegionStrand` is `POS` or `NEG` and whose fragment orientation is known.
3. For `EXON-INTRON`, compute the same deconvolution independently for each eligible
   boundary side on unambiguous exon regions.
4. Rows or sides with `RegionStrand.NONE`, `AMBIG`, missing fragment strand, or
   `c == 0` are strand-uninformative. They should not enter the strand-corrected
   estimate. If the entire run lacks strand-informative data, fall back to the current
   geometry-only estimator for backwards compatibility.
5. Estimate density from deconvolved rows with inverse-variance weighted least squares:

```text
rho_hat = sum_i(E_i * N_g_hat_i / Var_i) / sum_i(E_i^2 / Var_i)
Var(rho_hat) = 1 / sum_i(E_i^2 / Var_i)
```

where `E_i` is the effective length for that row:

- contained `L_eff` for intron rows;
- `B_cross(splicing_anchor_tolerance)` for each eligible boundary side.

This is the continuous, information-based fusion point. The estimator uses strand
evidence when it is precise, fades it out when it is noisy, and falls back structurally
when it is unavailable.

Longer-term exact version: replace the moment/WLS estimator with the profiled Poisson
likelihood

```text
Y_same_i ~ Poisson(0.5 * rho * E_i + p * R_i)
Y_opp_i  ~ Poisson(0.5 * rho * E_i + q * R_i)
R_i >= 0
rho >= 0
```

where `R_i` is a per-region RNA nuisance parameter. The profiled likelihood has the
same important limit: genic rows contribute no information about `rho` when `p = 0.5`.
The WLS moment estimator is simpler and should be enough for the first pass.

## 6. Fragment-Length Model Impact

The current gDNA FL source path sums histograms from `INTERGENIC`, `INTRON`, and
`EXON-INTRON` masks. In nRNA-heavy stranded runs, the intron and boundary histograms can
also be RNA-contaminated, which biases the gDNA FL model before density estimation.

Implementation should therefore add orientation-aware FL diagnostics from the start,
even if full deconvolved FL fitting is a second PR.

Preferred staged approach:

1. Preserve legacy `fl_hist` exactly.
2. Add `fl_hist_by_orientation` for mask categories where a single informative region
   orientation can be assigned without ambiguity.
3. Initially keep the current FL model path unless tests show FL contamination remains
   a dominant error after density deconvolution.
4. If needed, build deconvolved gDNA FL counts per bin using the same estimator:

```text
G_hat_bin = T_bin - (Y_same_bin - Y_opp_bin) / c
```

with stronger smoothing because per-bin counts are sparse. At perfect stranding this
reduces to `2 * minor_orientation_hist`.

This keeps the first implementation focused while leaving the data needed to debug and
fix FL contamination without another native payload change.

One assumption should be monitored explicitly: gDNA is treated as 50/50 across fragment
orientation. Intergenic fragments do not have a transcript strand, but their genomic
POS/NEG balance can still be reported as a QC check. If real data shows strong gDNA
orientation bias, the equations can be generalized by replacing `0.5` with an estimated
background orientation probability.

## 7. Implementation Plan

### Phase 1: Native Region Strand Plumbing

Files:

- `src/rigel/native/calibration/region_index.h`
- `src/rigel/native/bam_scanner.cpp`
- `src/rigel/pipeline.py`
- `tests/test_calibration_accumulator.py`

Changes:

1. Extend `RegionIndex::set(...)` to accept `const uint8_t* strands` in lockstep with
   `ref_ids`, `starts`, `ends`, and `type_masks`.
2. Store `std::vector<uint8_t> strands_` and expose `uint8_t strand(int32_t rid) const`.
3. Extend the nanobind `BamScanner.set_regions(...)` binding to accept a `uint8`
   strand array.
4. Update `_wire_calibration_regions` to pass `region_df["strand"]` after applying the
   same `(ref_id, start)` sort permutation as the other arrays.
5. Update set-region contract tests for length mismatch, dtype mismatch, and one
   end-to-end scan.

### Phase 2: Orientation-Split Payload

Files:

- `src/rigel/native/calibration/accumulator.h`
- `src/rigel/native/calibration/accumulator.cpp`
- `src/rigel/native/bam_scanner.cpp`
- `src/rigel/calibration/scan_payload.py`
- `src/rigel/calibration/_arrays.py`
- `tests/test_calibration_accumulator.py`

Add orientation buckets while keeping all legacy collapsed arrays:

```text
ORIENT_RNA_SAME = 0
ORIENT_RNA_OPP  = 1
ORIENT_UNINF    = 2
N_ORIENT        = 3
```

Suggested payload additions:

```text
per_region_counts_by_orient: int64[R, 8, 3]
u_left_by_orient:            int64[R, 3]
u_right_by_orient:           int64[R, 3]
fl_hist_by_orient:           int64[8, 3, FL_HIST_N_BINS]   # diagnostic/staged
```

Routing rules:

1. Collapsed arrays are still incremented exactly as today.
2. For each per-region or boundary increment, read `region_index.strand(rid)`.
3. If the region strand is `POS` or `NEG` and the fragment strand is `POS` or `NEG`,
   route to `RNA_SAME` when they match and `RNA_OPP` otherwise.
4. Otherwise route to `UNINF`.
5. The observation site should pass a resolver-derived fragment strand into
   `CalibrationAccumulator::observe(...)`. Prefer the existing `result.exon_strand`
   convention over deriving from the first exon block.
6. Worker merge must add all new arrays and remain byte-identical between one and many
   scan threads.

### Phase 3: Python Strand Estimator

Files:

- `src/rigel/calibration/density_global.py`
- `src/rigel/calibration/_orchestrator.py`
- `src/rigel/calibration/_result.py`
- `src/rigel/pipeline.py`
- `tests/test_density_global.py`
- `tests/test_calibrate_orchestrator.py`

Changes:

1. Thread the finalized `StrandModel` or a compact frozen summary
   (`p_r1_sense`, `p_r1_sense_variance`, `n_observations`) from `scan_and_buffer` into
   `calibrate(...)` and `compute_global_densities(...)`.
2. Implement a helper that takes same/opposite counts plus `p` and returns
   `(n_gdna_hat, variance, n_rna_hat, diagnostics)`.
3. Add `_density_intron_strand_aware(...)` and
   `_density_exon_intron_strand_aware(...)` using inverse-variance WLS.
4. Preserve exact legacy results when orientation arrays are missing, no strand model is
   available, or all genic rows are uninformative.
5. Extend summary diagnostics without breaking existing keys. Suggested additions:

```text
strand_aware: true/false
p_r1_sense
strand_specificity
n_strand_observations
INTRON.n_strand_informative_regions
INTRON.n_strand_uninformative_regions
INTRON.strand_precision
EXON-INTRON.n_strand_informative_sides
EXON-INTRON.strand_precision
```

`GlobalGdnaDensity.n_fragments` is currently an `int`. Deconvolved strand counts can be
fractional. Either broaden that field to `float` with a compatibility note, or add
`n_gdna_estimated` while leaving the legacy integer numerator intact for summaries.

### Phase 4: Local Prior Diagnostics

Files:

- `src/rigel/calibration/locus_prior.py`
- `src/rigel/calibration/_arrays.py`
- `tests/test_per_locus_gdna_mass.py`
- `tests/test_locus_partition.py`

The canonical Phase 2 prior currently uses global-only expected gDNA counts, while the
locoregional path remains diagnostic. After global strand-aware densities are stable:

1. Teach `PayloadArrays` to expose oriented intron and boundary counts.
2. Update diagnostic local `n_gdna_intron` and `n_gdna_boundary_observed` to use
   strand-deconvolved evidence when the local rows are informative.
3. Keep `expected_gdna_count_global(...)` global-only and observation-free unless the
   prior design is explicitly revised.

## 8. Validation Plan

Unit tests:

- Region strand is carried from `region_df` into native `RegionIndex` in sorted order.
- Orientation routing covers `POS/POS`, `POS/NEG`, `NEG/POS`, `NEG/NEG`, `NONE`, and
  `AMBIG` regions.
- R1-sense and R1-antisense protocols both route the RNA-major bucket correctly through
  the Python estimator.
- `p=1.0` gives `N_g=2*minor`.
- `p=0.5` returns zero strand precision and triggers the structural fallback.
- Sparse zero-count rows do not produce negative counts, NaNs, or infinite densities.
- Worker merge equality includes every new orientation array.

Scenario tests:

- Target `tests/scenarios/test_nrna_double_counting.py` first. The xfailed
  `nrna >= 30 and gdna <= 20` cases at `ss=0.9` and `ss=1.0` are the primary success
  criteria.
- Keep `ss=0.65` as an intermediate-information stress case. It should improve when
  counts are sufficient, but the test should not demand perfect separation.
- Preserve the clean `nrna_none` synthetic acceptance gates from the implicit-splice
  validation report.
- Run `tests/test_calibration_accumulator.py`, `tests/test_density_global.py`,
  `tests/test_calibrate_orchestrator.py`, and scenario tests before broad goldens.

Benchmark checks:

- gDNA-poor, nRNA-present stranded conditions should no longer inflate global intron or
  boundary gDNA density.
- Perfectly stranded, zero-gDNA nRNA simulations should estimate near-zero gDNA from
  genic strand-corrected channels.
- Unstranded simulations should remain close to the current baseline.
- Existing high-gDNA `nrna_none` conditions should not regress.

## 9. Key Risks and Mitigations

- Protocol direction bug: using `SS` alone loses whether RNA is same or opposite
  strand. Use `p_r1_sense` for math and `strand_specificity` only for summaries.
- First-exon strand shortcut: deriving fragment orientation from `obs_exons[0]` can
  misroute complex fragments. Pass `result.exon_strand` or an explicitly validated
  fragment strand.
- Double-counting evidence: do not combine raw total counts and orientation-split
  counts as independent observations for the same region.
- Ambiguous antisense loci: `RegionStrand.AMBIG` is not deconvolvable with a single
  RNA nuisance parameter. Route to `UNINF` and keep fallback behavior.
- Fractional deconvolved numerators: schema and summaries currently assume integer
  numerators. Update the schema deliberately instead of rounding away information.
- FL contamination: add orientation histograms early so the gDNA FL model can be
  decontaminated if density-only correction is insufficient.

## 10. Exit Criteria

The first implementation can be considered complete when:

1. Legacy collapsed calibration payloads remain byte-identical when summed across the
   new orientation buckets.
2. No-strand or unstranded runs preserve current density outputs within numerical
   tolerance.
3. Highly stranded nRNA-heavy toy scenarios no longer over-allocate gDNA in the global
   intron and boundary densities.
4. The summary JSON reports how much strand information was available and how much
   precision it contributed.
5. All new native changes are rebuilt with `pip install --no-build-isolation -e .` and
   the focused calibration/scenario tests pass.
