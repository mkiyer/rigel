# gDNA Regional Exposure Plan v4.3

**Status**: implementation plan, not implemented.  
**Date**: 2026-05-19  
**Supersedes**: `gdna_exposure_plan_v4.2.md`, `gdna_exposure_plan_v4.1.md`,
`gdna_exposure_plan_v4.1_implementation.md`,
`gdna_exposure_revised_plan_v3a.md`.

v4.3 is the clean denominator-only implementation plan. It removes the old
class-normalized exposure idea entirely. There is one molecular template, one
capture field, one global normalization factor, and one footprint-to-exposure
rule used for gDNA, nascent RNA, and mature mRNA.

---

## 0. Decisions

1. **One capture field.** Hybrid capture enriches genomic regions without
   knowing whether the molecule came from chromosomal DNA, nascent RNA, or
   mature RNA. Exposure is therefore a single genomic function `A(x)`.
2. **One normalization factor.** Region classes can still help compute the
   correct evidence numerator and opportunity denominator, but they must not
   create separate reference scales. `rho_ref` is global across all regions.
3. **Denominator only.** Production EM does not apply per-fragment
   `log A(x_u)` numerator weights. Exposure enters through EM effective
   lengths.
4. **One footprint rule.** Every component is represented by genomic blocks:
   gDNA locus blocks, nascent RNA spans, or mature mRNA exon blocks. Project
   those blocks onto `A(x)` and compute one component exposure weight.
5. **Raw public output.** Public `effective_length`, TPM, gene TPM, and nRNA
   TPM remain unweighted for cross-sample comparability. Weighted lengths are
   EM-only diagnostics.

The previous class-specific reference quantiles were a relic of the broken
gDNA-only numerator implementation. They tried to protect one component from a
bad asymmetric likelihood. Once the model is denominator-only and every
component receives the same kind of exposure correction, separate class
normalizers are conceptually wrong.

---

## 1. What Region Classes Still Do

Region classes are not part of the exposure scale. They are bookkeeping for how
we measure conservative gDNA evidence:

- INTERGENIC and INTRON regions use contained-fragment evidence.
- EXON-INTRON regions use boundary-crossing evidence.
- Each class therefore has its own way to compute `Y_r` and `E_r`.

After `Y_r` and `E_r` are computed, class membership ends. Every region now has
the same physical quantity:

```text
rho_r = Y_r / E_r
```

interpreted as an estimate of:

```text
rho_gdna_bulk * A_r
```

The plan assumes those evidence channels are comparable after their proper
opportunity denominators are applied. If validation later proves they are not,
the fix is to improve `Y_r`/`E_r` calibration, not to hide the mismatch behind
within-class normalization.

---

## 2. Phase 0: Baselines And Gates

Before implementation, record:

- current `--regional-exposure off` test/golden behavior;
- synthetic 24-condition pool errors;
- VCaP same-version uniform confusion matrix;
- current numerator-era regional confusion matrix as historical context only.

Acceptance gates:

| Gate | Requirement |
|---|---|
| Uniform mode | `--regional-exposure off` unchanged. |
| Auto-uniform mode | numerically identical to `off`. |
| Synthetic 24 | no mRNA/nRNA/gDNA pool regression beyond agreed tolerance. |
| VCaP gDNA -> RNA | <= uniform baseline + 1 percentage point. |
| VCaP gDNA -> nRNA | <= uniform baseline + 0.5 percentage point. |
| VCaP gDNA -> mRNA | <= uniform baseline + 0.5 percentage point. |
| VCaP RNA -> gDNA | within +/- 0.5 percentage point of uniform baseline. |

Do not tune tolerances to pass. If gates fail, inspect the exposure field and
component exposure weights first.

---

## 3. Phase 1: Build One Global Exposure Field

Modify `src/rigel/calibration/_regional_exposure.py`.

### 3.1 Evidence and opportunity

Reuse the existing code that computes per-region `Y` and `E`:

- intergenic `Y` from `payload_arrays.intergenic_per_region`;
- intron `Y` from strand-corrected intronic orientation counts;
- exon-boundary `Y` from strand-corrected boundary evidence;
- intergenic/intron `E` from contained FL-aware opportunity;
- exon-boundary `E` from boundary-crossing opportunity.

This is the only place where region-type-specific evidence logic remains.

### 3.2 Single global rate and shrinkage

For all regions with `E_r > 0`:

```text
rho_global = sum_r Y_r / sum_r E_r
```

Use one global shrinkage estimate, not one per class. The simplest
implementation reuses the existing method-of-moments helper on all eligible
regions:

```python
kappa_global = estimate_kappa(Y[valid], E[valid], rho_global).value
rho_hat_r = (Y_r + kappa_global * rho_global) / (E_r + kappa_global)
```

For `E_r == 0`, set `rho_hat_r = rho_global` and `A_r = 1` unless a later
diagnostic proves this is unsafe.

### 3.3 One global reference

Compute one reference over all eligible regions:

```text
rho_ref = weighted_quantile(rho_hat_r, weights=E_r, q=0.95)
```

Then:

```text
log_A_r = clip(log(rho_hat_r) - log(rho_ref), LOG_A_FLOOR, 0)
A_r = exp(log_A_r)
```

There is no `rho_ref_per_class`. There is no `signal_per_class`. There is no
`max(rho_ref_per_class.values())`.

### 3.4 Auto-uniform

Auto-uniform is global. Compute one spread diagnostic over all eligible
regions. If the observed global spread is no larger than the null/noise spread,
return uniform exposure (`A_r = 1`). Do not partially attenuate by class.

### 3.5 Summary diagnostics

Emit:

- `rho_global`;
- `rho_ref`;
- `kappa_global` and fallback reason;
- global q05/q50/q95 of `rho_hat`;
- per-class q05/q50/q95 diagnostics only, explicitly labeled not used for
  normalization;
- number of regions at floor.

Tests:

- Regions with identical `rho_hat` get identical `A` regardless of class.
- A region at global `rho_ref` gets `A = 1` regardless of class.
- A region 100x below global `rho_ref` gets `A ~= 0.01` regardless of class.
- No code path computes or uses class-specific reference quantiles for weights.

---

## 4. Phase 2: One Footprint Exposure Rule

Add a small footprint helper in `src/rigel/calibration/_exposure.py`:

```python
def footprint_exposure_weight(
    blocks: list[tuple[int, int, int]],  # (ref_id, start, end)
    exposure: "RegionalGdnaExposure",
    *,
    min_weight: float = 1.0e-4,
) -> float:
    """Return bp-weighted mean exposure over a component footprint."""
```

Algorithm:

```text
weighted_bp = sum_blocks exposure.weighted_length_on_ref(ref_id, start, end)
raw_bp      = sum_blocks (end - start)
weight      = weighted_bp / raw_bp
```

Rules:

- merge overlapping blocks per reference before summing, so gDNA loci do not
  double-count overlapping locus intervals;
- if `raw_bp <= 0`, return `1`;
- clamp to `[min_weight, 1]`;
- uniform exposure returns `1`.

This is intentionally a footprint-average correction. It is not exact
fragment-midpoint integration. That exact refinement is deferred until the
simple consistent model has been tested.

Tests:

- single block entirely in `A=0.5` region gives `0.5`;
- two exons of equal length in `A=1` and `A=0.5` regions give `0.75`;
- overlapping gDNA blocks are merged before averaging;
- uniform exposure gives `1` for every footprint.

---

## 5. Phase 3: Apply The Same Rule To Every Component

### 5.1 gDNA locus

For each `MultiLocus`, build footprint blocks from its constituent `Locus`
intervals:

```text
gdna_weight_M = footprint_exposure_weight(multilocus_blocks, exposure)
L_g_em        = L_g_raw * gdna_weight_M
```

`L_g_raw` remains the current FL-aware overlap effective length from
`gdna_eff_len_for_loci(...)`. v4.3 does not refactor that geometry.

### 5.2 Nascent RNA

Nascent RNA is unspliced by definition, so its footprint is the contiguous
genomic span of the nRNA entity. This applies whether the current index stores
the entity as a synthetic nRNA row, an annotated nRNA row, or a shared nRNA
parent row:

```text
nrna_weight_t = footprint_exposure_weight([(ref, start, end)], exposure)
L_n_em        = L_n_raw * nrna_weight_t
```

The raw length `L_n_raw` is the existing nRNA effective length for that entity
row. Do not use mature exon blocks for rows marked `is_nrna`; even if the row
has exon annotations, the nRNA component footprint is its full genomic span.

### 5.3 Mature mRNA

Mature mRNA uses exon blocks and skips introns:

```text
mrna_weight_t = footprint_exposure_weight(exon_blocks_for_t, exposure)
L_m_em        = L_m_raw * mrna_weight_t
```

Single-exon mRNA is the same algorithm with one block.

### 5.4 Implementation helper

Add:

```python
def transcript_exposure_weights(index, exposure) -> np.ndarray:
    """Return per-transcript EM exposure weights for all transcript rows."""
```

Use `index.t_df["is_nrna"]`, transcript reference ids, transcript starts/ends,
and `index.build_exon_csr()`:

- rows marked `is_nrna` use one contiguous `[start, end)` span;
- all other rows use exon blocks from the exon CSR;
- future non-row nRNA span tables should call the same footprint helper on
   their contiguous spans.

Tests:

- all transcript classes use the same helper;
- synthetic and annotated nRNA rows use the contiguous span, not exon blocks;
- multi-exon mRNA skips introns;
- uniform exposure gives all ones.

---

## 6. Phase 4: EM/Output Length Split

Modify `TranscriptGeometry` in `src/rigel/config.py`:

```python
effective_lengths: np.ndarray          # raw public lengths
effective_lengths_em: np.ndarray | None = None
```

Modify `AbundanceEstimator` in `src/rigel/estimator.py`:

- `_t_eff_len_output` stores raw public lengths;
- `_t_eff_len_em` stores EM lengths or aliases output lengths;
- native EM uses `_t_eff_len_em`;
- public `effective_lengths`, `quant.feather`, `nrna_quant.feather`, and gene
  outputs use `_t_eff_len_output`.

Pipeline setup:

```text
t_eff_len_raw = current_effective_lengths(...)
t_weight      = transcript_exposure_weights(index, exposure)
t_eff_len_em  = max(t_eff_len_raw * t_weight, 1.0)
```

Tests:

- with no EM lengths, output and EM behavior are unchanged;
- with synthetic EM lengths, native EM uses them while public output remains
  raw.

---

## 7. Phase 5: gDNA Prior Rescaling

Modify `PriorTable` in `src/rigel/calibration/locus_prior.py`:

```python
gdna_prior_count: np.ndarray       # eta_raw diagnostic
gdna_prior_count_em: np.ndarray    # eta passed to native EM
gdna_eff_len: np.ndarray           # L_g_em, EM denominator
gdna_eff_len_unweighted: np.ndarray
gdna_em_exposure_weight: np.ndarray
```

In `assemble_priors(...)`:

```text
L_g_raw = gdna_eff_len_for_loci(...)
W_g     = footprint_exposure_weight(multilocus_blocks, exposure)
L_g_em  = max(L_g_raw * W_g, 1.0)
eta_raw = expected_gdna_count_global(...)
eta_em  = eta_raw * L_g_em / L_g_raw
```

Pass `gdna_prior_count_em` to native EM. Keep `gdna_prior_count` as the raw
diagnostic. Stop emitting `gdna_prior_count_regional`; it is not part of this
model.

Tests:

- `gdna_prior_count_em / gdna_eff_len == gdna_prior_count / gdna_eff_len_unweighted`;
- uniform exposure gives identical raw and EM values;
- native EM receives `gdna_prior_count_em`.

---

## 8. Phase 6: Remove Numerator Weighting From Production

In `src/rigel/pipeline.py`:

- remove the production call to `_apply_unit_gdna_weights(...)`;
- do not mutate `em_data.gdna_log_liks`;
- do not add nRNA or mRNA numerator mutation;
- keep `ScoredFragments.genomic_midpoint`, `log_weights_for_positions`, and
  `RegionalWeightApplicationStats` as dead scaffolding for now if deleting
  them would cause schema/native churn.

After VCaP gates pass, delete the dead scaffolding in a cleanup PR.

Tests:

- sentinel/monkeypatch proves `_apply_unit_gdna_weights` is not called;
- regional mode changes EM denominators, not fragment log-likelihood arrays.

---

## 9. Phase 7: Diagnostics And Output

Public outputs remain raw:

- `quant.feather effective_length`, `tpm`, `tpm_total_rna`;
- `nrna_quant.feather effective_length`, `tpm`;
- `gene_quant.feather effective_length`, `tpm`.

Add diagnostics:

| Column | File | Meaning |
|---|---|---|
| `em_effective_length` | `quant.feather` | transcript EM denominator |
| `em_exposure_weight` | `quant.feather` | `em_effective_length / effective_length` |
| `em_effective_length` | `nrna_quant.feather` | nRNA EM denominator |
| `em_exposure_weight` | `nrna_quant.feather` | nRNA EM/raw ratio |
| `gdna_eff_len_unweighted` | `loci.feather` | raw gDNA opportunity |
| `gdna_eff_len` | `loci.feather` | EM gDNA opportunity |
| `gdna_em_exposure_weight` | `loci.feather` | gDNA footprint exposure weight |
| `gdna_prior_count` | `loci.feather` | raw prior diagnostic |
| `gdna_prior_count_em` | `loci.feather` | prior passed to EM |

Tests:

- public TPM/effective-length columns are unchanged in uniform mode;
- diagnostic ratios are 1 in uniform mode;
- regional mode can change diagnostics without changing public length columns.

---

## 10. Phase 8: Validation

Run in order:

1. Focused unit tests from Phases 1-7.
2. Full test suite.
3. Synthetic 24-condition benchmark.
4. VCaP same-version uniform vs v4.3 auto confusion matrix.

If VCaP fails:

1. Inspect global exposure field distribution and class diagnostic summaries.
2. Inspect component exposure weights for gDNA, nRNA, and mRNA in failing loci.
3. Check whether the field evidence channels (`Y_r`, `E_r`) are actually on a
   common scale.
4. Decide whether the next experiment is better field calibration, exact
   fragment-midpoint integration, or performance optimization.

Do not reintroduce numerator weighting. Do not return to class-specific
normalization unless the biological model changes; it should be considered a
diagnostic failure mode, not a solution.

---

## 11. Deferred Work

- Gamma-Poisson or other improved field estimators.
- Prefix-sum/batched integration for speed.
- Exact fragment-midpoint effective-length integration.
- Native cleanup of `genomic_midpoint`.
- Deletion of dead numerator scaffolding.

These are follow-ups after the simple globally normalized denominator-only
model is tested.

---

## 12. Success Criterion

v4.3 succeeds if a single globally normalized exposure field plus one
footprint exposure rule fixes the VCaP numerator-era regression without
breaking uniform behavior or synthetic pool accuracy.

The plan is intentionally simple: one field, one reference, one footprint
projection, denominator-only EM.