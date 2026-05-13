# Density Effective-Length Consistency Fix (R1)

**Date:** 2026-05-06
**Scope:** Calibration density estimation (INTERGENIC, INTRON) and
locoregional shrinkage. Replaces the FL-mean overlap effective length
with the FL-PMF-weighted containment effective length, the MLE for the
underlying Poisson process.

---

## 1. The bug, in one sentence

The numerator of $\hat\rho$ counts fragments **contained** in regions
(C++ accumulator: every aligned block must fall in mask-bit
positions), while the denominator uses
$L_{\text{eff,overlap}} = \text{span} + \bar\ell - 1$, the count of
genomic positions where a fragment of *mean* length **overlaps** the
region by at least one base. Numerator and denominator measure
different geometric events; the result is biased low whenever
$\bar\ell$ is non-negligible relative to region span (introns,
exon-intron boundary projection).

## 2. Verified bias on synthetic mini-genome

| Region type | True ρ (simulated) | Estimated ρ | Ratio | Reason |
|-------------|-------------------:|------------:|------:|--------|
| INTERGENIC  | 0.10               | 0.0996      | 0.996 | spans ≫ FL → overlap ≈ containment |
| INTRON      | 0.10               | 0.0787      | 0.787 | spans ~3 kbp, FL ~350 → 21 % low |

Predicted ratio under the bug: $(\overline{\text{span}}-\bar\ell+1) / (\overline{\text{span}}+\bar\ell-1) = 2550/3250 = 0.785$ — matches observation.

## 3. The fix

Replace the per-region scalar denominator with the FL-PMF-weighted
containment expectation:

$$
L_{\text{eff,contained}}(r) \;=\; \sum_{\ell=1}^{\text{span}_r} h_{\text{gdna}}(\ell)\,(\text{span}_r - \ell + 1)
$$

This is the **MLE** denominator for a Poisson process with rate ρ
fragments/bp and FL distribution $h_{\text{gdna}}$:
$E[N_r] = \rho \cdot L_{\text{eff,contained}}(r)$.

### Why Option A (PMF-weighted), not Option B (per-fragment)

Option B (sum of $L_{\text{eff}}(r,\ell_f)$ over observed fragments)
is biased: its expectation reduces to a function of region geometry
that does **not** depend on ρ. See §3 of the previous design dialog
for the algebra.

## 4. Implementation: leverage existing infrastructure

`FragmentLengthModel.compute_all_transcript_eff_lens(lengths)`
**already** computes exactly this quantity for transcripts:

```python
eff = (L + 1) * (CDF[L] - P(0)) - CMOM[L]      # = sum_{l=1..L} P(l)*(L-l+1)
```

Vectorized, eCDF-cached, battle-tested in EM scoring. We reuse it
verbatim — no new arithmetic helper, no new C++ pass, no payload
changes.

The cost: one method call replaces the `l_eff_overlap` line in two
sites (one global, one per-locus). The `gdna_fl_mean: float` parameter
threading through `compute_global_densities`, `assemble_priors`,
`estimate_locus_gdna`, and `_shrink_one_type` is replaced with a
single `gdna_fl: FragmentLengthModel` parameter — narrower to write,
broader in information content.

## 5. Files touched and exact changes

### 5a. `src/rigel/calibration/density_global.py`

- **Delete** `l_eff_overlap` (no longer the canonical formula).
- **Change** `compute_global_densities()` signature:
  `gdna_fl_mean: float` → `gdna_fl: FragmentLengthModel`.
- **Replace** the `leff = l_eff_overlap(spans, gdna_fl_mean)` line with
  `leff = gdna_fl.compute_all_transcript_eff_lens(spans)`.
- **Keep** `_density_exon_intron` unchanged — its denominator
  (`sides * gdna_fl_mean`) reflects a different geometry (a capture
  window of width $\bar\ell$ crossing a boundary), not a containment
  computation. Pass `gdna_fl.mean` through.
- **Update** `estimate_global_gdna_fragments()` similarly: replace
  `leff = l_eff_overlap(...)` with the FL-PMF call.
- **Update** the `GlobalDensityTable` schema: `gdna_fl_mean` field
  remains (used downstream for diagnostics and EXON-INTRON geometry);
  no new fields needed.

### 5b. `src/rigel/calibration/density_loco.py`

- **Delete** the `l_eff_overlap` re-export. No remaining callers.

### 5c. `src/rigel/calibration/locus_prior.py`

- **Change** `estimate_locus_gdna()` signature:
  `gdna_fl_mean: float` → `gdna_fl: FragmentLengthModel`. (The
  function then computes `gdna_fl.mean` locally for the EXON-INTRON
  branch.)
- **Replace** the `leff = cl_len + (gdna_fl_mean - 1.0)` line with
  `leff = gdna_fl.compute_all_transcript_eff_lens(cl_len.astype(np.int64))`.
- **Change** `_shrink_exon_intron()` to take the scalar mean (it's
  already only computing capture-window arithmetic).
- **Change** `assemble_priors()` to take `gdna_fl: FragmentLengthModel | None`.
  The `None` path falls back to `cal.fl_models.gdna`.

### 5d. `src/rigel/calibration/__init__.py`

- **Drop** `l_eff_overlap` from re-exports.

### 5e. `src/rigel/calibration/_orchestrator.py`

- **Change** the `compute_global_densities(...)` call: pass
  `gdna_fl=fl_models.gdna` instead of `gdna_fl_mean=fl_models.gdna.mean`.

### 5f. `src/rigel/pipeline.py`

- The `assemble_priors(...)` call: no signature change needed (it
  reads `gdna_fl` from the `CalibrationResult`).
- The `estimate_global_gdna_fragments(...)` call: no signature change
  (the function reads `gdna_fl` from `global_densities` indirectly via
  the calibration result we already pass).

## 6. Tests touched

| Test file | Action |
|-----------|--------|
| `tests/test_density_global.py`   | Pass a stub `FragmentLengthModel` instead of `gdna_fl_mean=...`. Update expected ρ values for cases where region spans are small. |
| `tests/test_per_locus_gdna_mass.py` | Same: stub model instead of mean. |
| `tests/test_assemble_priors.py`     | Same. |
| `tests/test_density_loco.py`        | Drop any `l_eff_overlap` reference. |
| `tests/golden/`                     | Regenerate via `pytest --update-golden` after the code change (loci_df `gdna_prior` column will shift slightly). |

A small helper in `tests/conftest.py` (or a per-test fixture) builds a
finalized `FragmentLengthModel` from a (mean, std) pair — used by all
calibration tests that previously passed `gdna_fl_mean=...`.

## 7. Validation

- **Unit:** new `tests/test_density_global.py::test_intron_density_unbiased`
  builds 100 identical 1 kbp intron regions, samples 1000 fragments
  from a known $h(\ell)$ via direct simulation, asserts $|\hat\rho/\rho-1| < 0.05$.
- **Integration:** re-run `scripts/sim/evaluate_suite.py` on the
  10-condition synthetic sweep. Expected:
  - $\rho_{\text{intron}} / \rho_{\text{intergenic}}$ moves from 0.79 → 1.00 ± 0.02.
  - Global gDNA rate at high contamination stays accurate (intergenic
    estimate already matches truth).
  - Per-locus $\alpha_\text{gdna}$ priors increase modestly for
    intron-rich loci; mRNA quantification correlation should be
    unchanged or slightly improved.
- **Regression:** golden file diffs reviewed for sane direction (gdna
  estimates rise slightly in intron-rich loci, unchanged in
  intergenic-dominated loci).

## 8. Cleanups bundled in this PR

1. **Drop the `l_eff_overlap` re-export** from `density_loco.py` (was
   only there to spare callers an import — no callers remain after R1).
2. **Promote the FL-PMF eff-length to a single named call site.**
   Previously `l_eff_overlap` was three lines each in three modules;
   now it is one method call to a finalized FL model.
3. **Remove the misleading "forbidden alternative" comment** in
   `l_eff_overlap`'s docstring (the function itself is being deleted,
   resolving the confusion at the source).
4. **Narrow the scalar-mean threading.** `gdna_fl_mean: float`
   appeared in 4 function signatures purely as a denominator
   ingredient. Replace with the full FL model where containment is
   needed; keep the scalar only at the EXON-INTRON site that needs it
   geometrically. Net: −2 parameters from public API.

## 9. Out of scope (deferred)

- Reworking the EXON-INTRON capture-window geometry to use the FL
  PMF (the `sides * gdna_fl_mean` denominator). The geometry there is
  fundamentally different — a probe window crossing a boundary — and
  scalar-mean is a defensible approximation. Revisit if the synthetic
  sweep shows residual bias in the EXON-INTRON ratio after R1.
- Per-region $\kappa$ estimation already uses `sub_leff`; under the
  new denominator $\kappa$ values shift but the estimator's
  invariance is preserved (rate × leff still equals expected count).

## 10. Acceptance criteria

- [ ] All existing tests pass (with golden regen for output schema).
- [ ] New `test_intron_density_unbiased` passes.
- [ ] Synthetic sweep: $\rho_{\text{intron}}/\rho_{\text{intergenic}} \in [0.95, 1.05]$.
- [ ] No new public API surface; net −1 helper function (`l_eff_overlap`
      removed), net −2 parameters in calibration entry points.
