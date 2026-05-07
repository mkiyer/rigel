# Calibration Fix Plan — Post-Sweep Validation

## Summary

Five targeted fixes arising from the gDNA calibration sweep (2026-05-05).
Ordered by dependency (bottom-up): fixes are independent except where noted.

---

## Fix A: Test script uses wrong `sj_strand_tag` (TRIVIAL)

**Root cause**: The sweep script passes `sj_strand_tag="ts"` but the oracle BAM
only sets the `XS` tag. The scanner looks for `ts`, finds nothing, and the strand
model gets zero spliced observations.

**Fix**: Use `sj_strand_tag="auto"` in the test script (and in general for oracle
BAMs). The auto-detection correctly reads spliced records, finds `XS`, and returns
`"XS"`. This is already the `PipelineConfig` default.

**Files changed**: `scripts/debug/test_gdna_calibration_sweep.py`

---

## Fix B: `loci_df["mrna"]` excludes unambiguous counts (REPORTING)

**Root cause**: `get_loci_df()` in `estimator.py` computes `mrna = rna_total - nrna`
where `rna_total` is the EM posterior sum (from `r["rna_total"]`). Unambiguous
(splice-confirmed) fragments bypass the EM and are never reflected in locus-level
mrna/total/gdna_rate columns.

**Design**: Add `unambig_mrna` counts per locus from `self.unambig_counts`:
- During `run_batch_locus_em_partitioned`, accumulate per-locus unambig counts
- In `get_loci_df()`, add unambig_mrna to the mrna column and total column
- Recompute gdna_rate = gdna / (mrna + nrna + gdna) with corrected mrna

This changes the semantics of loci_df to reflect *all* fragments, not just EM fragments.
The `n_em_fragments` column remains unchanged (still EM-only for diagnostic purposes).

**Files changed**: `src/rigel/estimator.py`

**Schema change**: Add column `count_unambig` to loci_df (per-locus unambiguous count).
Existing `mrna` column now includes both EM + unambig. `total` = mrna + nrna + gdna.

---

## Fix C: Global contamination rate includes intergenic gDNA (REPORTING)

**Root cause**: `gdna_contamination_rate` in `estimator.py` uses only `_gdna_em_total`
(EM-assigned gDNA within loci). Intergenic gDNA (which never enters loci) is ignored,
making the reported rate far lower than truth.

**Design**: After calibration is complete (post-`calibrate()`, post-EM), the pipeline
should compute a global gDNA estimate from calibration densities projected over the
entire genome's region partition. This requires:

1. **New function** `compute_global_gdna_estimate()` in `calibration/locus_prior.py`:
   - Takes `global_densities`, `region_df`, `gdna_fl_mean`
   - For each region type (INTERGENIC, INTRON, EXON):
     - Computes `leff_total = sum(region_span + fl_mean - 1)`  
     - INTERGENIC: `n_gdna_ig = rho_intergenic × leff_ig_total`
     - INTRON: `n_gdna_in = rho_intron × leff_in_total`
     - EXON: `n_gdna_ex = rho_exon_intron × leff_ex_total` (projected from boundary flux)
   - Returns `n_gdna_global` (predicted total gDNA fragments across entire genome)

2. **Compute library-wide gDNA rate** = `n_gdna_global / n_total_observed`
   where `n_total_observed = calibration.diagnostics.total()`

3. **Store on CalibrationResult** as `global_gdna_rate` property

4. **Report in pipeline**: Replace `estimator.gdna_contamination_rate` with
   calibration-derived global rate in the log line. Keep the EM-only rate
   available as a secondary diagnostic.

**Files changed**:
- `src/rigel/calibration/locus_prior.py` (new function)
- `src/rigel/calibration/_result.py` (new property)
- `src/rigel/pipeline.py` (use global rate in log)
- `src/rigel/estimator.py` (rename property for clarity)

---

## Fix D: Per-locus π_gdna includes exonic gDNA projection (PRIOR ACCURACY)

**Root cause**: `estimate_locus_gdna()` computes `n_gdna = n_ig + n_intron + n_exon_intron`
but the exon-intron branch uses `rho_loco_ei × leff_ei` where `leff_ei` is only the
**exonic** effective length. The problem is that the *density* (`rho_loco_ei`) comes
from boundary-flux observations (counts from edges), but the *projection* already
applies it to the full exonic L_eff. So the exonic gDNA IS already being predicted!

**Wait** — re-reading the code carefully:
```python
n_gdna_ei = rho_loco_ei * leff_ei  # leff_ei = sum of clipped exon effective lengths
```
This IS projecting boundary flux onto the full exon area. The underestimation comes
from the denominator: `n_obs` includes both gDNA AND mRNA fragments in the locus.
The `n_gdna_total / n_obs` calculation produces a fraction that includes mRNA in the
denominator, which dilutes the prior.

**Actual root cause**: At `gdna_frac=1.0`:
- `n_gdna_total = 1235` (density-predicted from non-exonic + exon boundary flux)
- `n_obs = 9210` (all fragments entering locus EM = 8074 mRNA + ~1136 gDNA)
- `π_gdna = 1235 / 9210 = 0.134`

But the TRUE rate within this locus = 1136 / 9210 = 0.123. So **the estimate is
actually CORRECT** (slightly overestimates, which is fine). The "underestimate"
appearance in the earlier analysis compared π_gdna against the *library-wide* rate
(0.50), which is wrong — π_gdna is a *per-locus* rate that should be compared to
the within-locus gDNA rate.

**Validation**: From the per_locus_gdna.tsv at gdna_1.00:
- Locus 0: π_gdna=0.134, EM result gdna=1068/9210=0.116 — **prior is reasonable**
- Locus 1 (tx_single, no mRNA): π_gdna=1.0, EM result gdna=778/778=1.0 — **correct**

**Conclusion**: The per-locus π_gdna prior is working correctly. No code change needed.
The confusion arose from comparing a per-locus prior to a library-wide rate.

However, we CAN improve by including **neighboring intergenic regions** (your response #1).
This is a separate enhancement that improves robustness for loci with no intronic regions.

### Enhancement: Extend locus region overlap to include neighboring intergenic regions

**Design**:
- In `estimate_locus_gdna()`, after overlapping regions within `[locus.start, locus.end]`,
  optionally extend the overlap query to include the adjacent intergenic region(s)
  on each side (one region left, one region right).
- `RegionIndexPy.overlap()` already returns sorted region IDs.
- Add a parameter `include_neighbors: bool = True` to `estimate_locus_gdna()`.
- The intergenic contribution from neighbors adds a more robust local density estimate.
- EB shrinkage already handles the case where local data is sparse by shrinking toward global.

**Deferred** — this is an improvement for robustness but NOT a correctness fix.
The current system already works well via EB shrinkage to global densities.

---

## Fix E: (No change needed) — RNA FL quality "weak" is expected

The RNA FL pool uses only spliced-annotated fragments (1,926 in this simulation).
The EB-shrinkage system correctly handles this by shrinking the weak RNA FL model
toward the global distribution. No code change needed.

---

## Implementation Order

1. **Fix A** (trivial, test script only) — Use `sj_strand_tag="auto"` in sweep script
2. **Fix B** (estimator.py) — Add unambig counts to loci_df
3. **Fix C** (calibration + pipeline) — Global gDNA rate from calibration densities
4. Re-run sweep to validate all fixes

Fixes B and C are independent and can be developed in parallel.

---

## Validation Plan

After implementing fixes A–C, re-run the gDNA calibration sweep:
- Fix A: Strand model should now train correctly (1,926 spliced-annotated observations)
- Fix B: `loci_df["mrna"]` should show ≈10,000 (not ≈8,074)
- Fix C: Reported contamination rate should match truth (±3% relative error)
- Overall mRNA accuracy should remain <2% error (currently already excellent)
