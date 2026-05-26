I read the v2 plan carefully. It's much improved over v1 — the additive-pseudocount parameterization is the right call, the carried-state edge case is honest, and removing gene-level grouping is correct. **Conditional certification: implement Phase 0–2 as written, but Phase 3 has 3 specification bugs and Phase 1 has 1 formula issue that should be fixed before code lands.** Findings below, in priority order.

---

## Blocking issues (fix before implementation)

### B1. Default `aggregate_prior_strength = 1.0` is dangerously hot

`budget = strength * (m_g + m_r)`. For an expressed locus with `m_r ≈ 1000`, `m_g ≈ 5`, this injects 1005 pseudocounts. With `N_rna ≈ 1000` from the E-step, the multiplier `(1 + a_r / N_rna) ≈ 2.0` — calibration is doubling every RNA component's posterior mass each M-step. For low-expressed loci (`N_rna = 2`, `a_r = 50`), the prior is 25× the data. **This is not a smooth knob; it is a hammer.**

The current `gdna_prior_count_em` ships with calibration prior on the same scale as `m_g` (typically 0.1–10 fragments per locus). The plan grants RNA the *symmetric* counter-prior, but the symmetry is at the scale of raw expected counts, which are radically asymmetric between gDNA and RNA in any expressed locus.

**Fixes (pick one, prefer the first):**

- **Decouple budget magnitude from raw count magnitude.** Use the *smaller* side as the budget: `budget_raw = 2 * min(m_g, m_r)`. This caps prior pull at 2× the weaker side of the calibration evidence, which is the side that actually constrains the share. The Beta MAP semantics are preserved (same `w_raw`) but the budget cannot dominate either side's data.
- Alternative: `budget_raw = sqrt(m_g * m_r)` (geometric mean). Same property, slightly less conservative.
- Default `aggregate_prior_strength = 0.1` and document that `1.0` is "trust calibration completely."

If the first fix is adopted, the default `strength = 1.0` becomes safe and the knob semantics ("one calibrated observation per locus on the constrained side") match user intuition.

### B2. `mu_rna` formula in Phase 1 risks double-counting

The formula:

```
rna_residual = max(contained_total - mu_gdna, 0)
mu_rna       = spliced_total + max(rna_lower, p_expr * rna_residual)
mu_rna       = min(mu_rna, spliced_total + contained_total)
```

`rna_lower` is the strand-derived lower bound on RNA in the region. Per the strand-deconvolution module (`strand_deconv.py::deconvolve_regions_by_strand`), it is computed from **all unspliced sense-vs-antisense balance**, and in some regions it can already include spliced contributions depending on how the accumulator emits counts. Adding `spliced_total` again creates a measurable double-count in expressed regions with strong strand signal — exactly the regions where `mu_rna` matters most.

**Fix:** make the decomposition explicit and non-overlapping. Two clean options:

```
# Option A: spliced is mature, residual is nascent/contained unspliced
mu_rna_spliced   = spliced_total
mu_rna_unspliced = max(rna_lower_unspliced, p_expr * rna_residual)
mu_rna           = mu_rna_spliced + mu_rna_unspliced
mu_rna           = min(mu_rna, spliced_total + contained_total)
```

This requires `rna_lower` to be restricted to the **unspliced** strand contribution; if it is not, factor it that way in `_categorize.py` first. Alternatively:

```
# Option B: ignore strand component entirely in mu_rna; rely on residual + spliced
mu_rna = spliced_total + p_expr * rna_residual
mu_rna = min(mu_rna, spliced_total + contained_total)
```

Option B is simpler and avoids the audit; only adopt Option A if benchmarks show strand signal must contribute. Either way, add a regression test on a synthetic region with known spliced + strand counts to lock the semantics.

### B3. Disabled-gDNA case still applies the RNA counter-pull

The plan zeroes `a_g_effective` when `!enable_gdna` but **leaves `alpha_rna_add` at its full magnitude**. The intent of `a_r` is to balance gDNA's pull; with no gDNA competitor, the residual RNA pull becomes a uniform inflation of RNA components.

In a single-MAP-step it has no effect (since RNA is the only group and the M-step normalization absorbs the scaling). But under SQUAREM with multiple components and length-asymmetric `log_eff_len`, the multiplicative `(1 + a_r/N_rna)` factor *does* shift the EM trajectory slightly because it amplifies before normalization.

**Fix:** when `!enable_gdna`, zero **both** sides:

```cpp
if (!enable_gdna) {
    alpha_gdna_effective = 0.0;
    alpha_rna_effective  = 0.0;
}
```

Document this clearly in Phase 3: the aggregate prior only acts when there is an actual aggregate split to constrain.

### B4. SQUAREM stabilization for grouped VBEM is under-specified

Plan says "Stabilization step calls `grouped_vbem_step()`." Good. But it doesn't specify what happens when SQUAREM's extrapolated `alpha` vector has the wrong gDNA/RNA *ratio* (which it generically will, because extrapolation is linear in component space, not in the Beta-constrained subspace).

The grouped step rebuilds the ratio correctly from raw_counts (which come from a fresh E-step using the extrapolated `alpha`), so this is mathematically fine — but the SQUAREM step-length search compares `‖F(x_p) - F(x_q)‖` to decide step length. If `F` is the grouped step, then off-manifold extrapolations get pulled hard back onto the Beta-MAP manifold, which can make SQUAREM perceive instability and shrink the step. Worth instrumenting.

**Fix:** add to Phase 5 validation an explicit metric: SQUAREM step-length acceptance rate on aggregate-prior loci vs. unpriored loci. If acceptance collapses, fall back to plain EM for those loci (the fallback path is already in the plan; just make sure the trigger is sensitive enough).

---

## Non-blocking issues (fix during implementation)

### N1. `logit_bias` semantics are nonlinear and surprising at extremes

`w = sigmoid(logit(clamp(w_raw, eps, 1-eps)) + bias)`. At `w_raw = 0` this becomes `sigmoid(logit(eps) + bias)`. For `eps = 1e-6` and `bias = +2`, this is ≈ `sigmoid(-13.8 + 2) ≈ 7e-6`. So strong positive bias barely moves loci where calibration says "no gDNA." That is *correct* behavior, but the document advertises bias as "operating point" without explaining that the bias scale is **odds**, not share. A user setting `bias = +1` will see vastly different effects in loci with different `w_raw`.

**Fix:** document that bias is on the log-odds scale and that its practical range is ±2. Optionally, gate the bias-induced share `w` with a calibration-confidence weight (e.g., `confidence = (m_g + m_r) / (m_g + m_r + κ)`) so loci with thin calibration evidence resist bias amplification.

### N2. Conflating mRNA and nRNA under one "RNA" aggregate is risky given known siphon

Per the project's known benchmark findings, the **nRNA siphon** is the dominant remaining failure mode in loci with high NTA1/TA1 ratios. The Phase 1 design pools mRNA and nRNA into a single aggregate and only **scales** their joint mass; it does not affect their *relative* partition. So in principle siphon is unchanged by the prior. Good.

But there is a second-order effect: if calibration's `mu_rna` is dominated by unspliced contained mass (i.e., looks like nRNA), the aggregate pull inflates the joint pool, which through the E-step's effective-length kernel preferentially routes mass to nRNA components (longer eff len). **This may worsen siphon, not improve it.**

**Fix:** Phase 5 must include an explicit siphon benchmark (NTA1/TA1 ratio sweep) before and after the redesign. If siphon worsens by even 5%, accelerate Phase 6 to split into two groups (annotated-RNA vs synthetic-nRNA) before merging Phase 3.

### N3. Cold-start (no carried state) edge case is unsafe

The plan's edge case ordering is:

```
q_t = previous RNA proportions if sum(previous RNA mass) > eps
q_t = warm-start RNA proportions if available
q_t = effective-length-weighted fallback only as a last resort
```

The "last resort" path creates uniform-ish RNA prior across all transcripts — exactly the false-floor pathology the redesign exists to prevent. With `a_r = 50` and `N_rna = 0` at iteration 0, this floors every isoform.

**Fix:** add an additional guard: **when `N_rna == 0` AND the carried state has no RNA mass, set `a_r_effective = 0` for that iteration.** Calibration cannot push expression onto isoforms it has no evidence for; the prior should respect that. The next iteration, once some RNA mass exists, the prior activates.

### N4. `assemble_priors()` signature change is not in the file list

Phase 2 says "Thread `em_config` into `assemble_priors()`". The current signature is `assemble_priors(multi_loci, em_data, index, calibration)`. Add `em_config: EMConfig` (or a smaller `PriorConfig` extracted from `EMConfig`). Update the file list to include pipeline.py where `assemble_priors` is called, and update the call site in `quant_from_buffer`.

### N5. PriorTable schema migration needs explicit feather/golden plan

Phase 4 says "Update goldens intentionally" but doesn't enumerate which golden files. At minimum:
- `tests/golden/loci.feather` — locus-level prior columns change.
- Any golden referencing `gdna_prior_count`.

Add an explicit list to Phase 4 and a one-checkpoint alias retention (`gdna_prior_count_em` aliased to `alpha_gdna_add` for one release).

### N6. The `min(mu_rna, spliced_total + contained_total)` upper bound is per-region but projection is summed

`mu_rna` is bounded per region by `spliced + contained`, but after geometry allocation to a `MultiLocus`, the sum across regions of `mu_rna` is upper-bounded only loosely. Not a bug, just worth a diagnostic: report `sum(mu_rna_alloc) / sum(observed_spliced + observed_contained)` per locus; if > 1 anywhere, the projection has a leak.

### N7. Native warm-start application: explicit math for VBEM state

Phase 3 says "Use that as `theta_init` for MAP or `alpha_init` seed for VBEM." For VBEM, `alpha_init` is unnormalized count mass. The warm start is built in count space (unambig totals + ambiguous shares), so it is already in the right space. Just clarify in the doc that VBEM does **not** re-normalize the warm start before SQUAREM.

### N8. Capture exposure double-counting risk

`gdna_eff_len` already absorbs `A_r` (capture exposure). `mu_gdna` is also calibration's posterior count *given* capture-aware exposure. So `alpha_gdna_add / gdna_eff_len` (the effective density the gDNA component sees in the E-step) divides A_r out of the prior count by an A_r-aware denominator — **fine, no double-count**. But add a sanity assertion: in an off-target locus, `alpha_gdna_add / gdna_eff_len` should match `ρ_off` to within the prior strength scale. Useful diagnostic, not a bug.

---

## Things the plan got right (worth preserving against revision pressure)

- **Additive pseudocount semantics** instead of standard Dirichlet `α - 1`. Correctly aligns with calibration's count-scale outputs.
- **`enable_gdna` strictly technical.** This is the right invariant; do not let benchmarks tempt you into folding calibration confidence back into it.
- **Separating `aggregate_prior_strength` from `gdna_prior_logit_bias`.** This is the v1 plan's biggest weakness fixed. Keep them separate.
- **Carried state for zero-RNA edge case.** Better than length-weighted uniform.
- **Phase 6 ban on gene-level priors.** Correct enforcement of the transcript-first invariant.

---

## Certification

**Certify Phase 0, 1, 2, 4, 5 for implementation as written, after applying B2 (mu_rna formula).**

**Do not certify Phase 3 as written.** Required changes before native ABI lands:
- B1: budget rescaling to `min(m_g, m_r)` or default `strength = 0.1`.
- B3: zero both `a_g` and `a_r` when `!enable_gdna`.
- B4: explicit SQUAREM step-length acceptance instrumentation.
- N3: zero `a_r` when both posterior and carried RNA mass are zero.

These are all small surface-area changes — none require re-architecting the design.

**Suggested order:** patch the plan with the four B-issues, then implement Phase 0 (tests) which will lock the corrected math, then continue through Phase 1–5 as planned.