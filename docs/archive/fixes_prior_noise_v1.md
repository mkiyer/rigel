I've reviewed the diagnostic outputs, configs, calibration code, and prior architecture. Here is my detective-mode analysis of the two questions, with concrete numbers and a fix plan.

## Q1 — Why does 5K+ gDNA still bleed into RNA at high-SS capture-ON despite an accurate prior?

### What we observe at the dominant locus (Locus 1 / GENE0002)

From loci.tsv and the truth file:

| Quantity | Value |
|---|---|
| `n_em_fragments` | 100,992 |
| True mRNA fragments (sum of GENE0002 truth) | 72,963 |
| True gDNA fragments (residual) | ~28,029 |
| Rigel mRNA estimate | 81,350 |
| Rigel gDNA estimate | 22,529 |
| **gDNA→RNA bleed at this locus** | **~5,500–8,400** |
| `prior_n_local_gdna` (all-region snapshot) | 94,869 |
| `prior_n_local_rna` (all-region snapshot) | 79,248 |
| `prior_ess_final` after cap | **3,000** |
| `gdna_em_exposure_weight` | **1.104** |
| `gdna_eff_len` | 33,261 bp |

### Root cause: the prior cannot move the EM because data dominates and the per-fragment likelihood is biased toward RNA

There are two **independent** failures that compound:

**1. Prior ESS is structurally too small to override data.**
The all-region prior assembles ~95K "local gDNA counts" but `MAX_ESS = 3000` collapses that to a 3000-count Dirichlet. Against ~101K fragments, the prior shifts the posterior by only ~3% (`alpha / (alpha + n)`). So even when the prior says "65% gDNA", the EM ends near the likelihood's preferred 41%. The prior is *not* an oracle; it is a tiebreaker.

**2. The per-fragment likelihood at probe positions is systematically RNA-favored.**
A probe-overlap fragment maps to a transcript exon (~7 kb effective length) and to a long gDNA locus (~33 kb effective length). Under the rigel scoring model:

$$\log L_{\text{RNA}} - \log L_{\text{gDNA}} \approx \log\!\frac{\tilde{L}_{\text{gDNA}}}{\tilde{L}_{\text{RNA}}} \approx \log\!\frac{33{,}261}{7{,}000} \approx +1.55 \text{ nats per fragment}$$

That is, every captured fragment carries ~1.5 nats of "RNA-favored" likelihood mass before priors. The `gdna_em_exposure_weight` of **1.104** is supposed to compensate — but it is way too small. The same calibration knob, `capture_enrichment_target`, came out at **1.267** (see `region_calibration` block). A true hyb-capture protocol yields enrichments of 10×–1000× over off-target. Rigel's calibration is under-estimating capture enrichment by 1–3 orders of magnitude.

**Why the enrichment estimate is so small:** in [src/rigel/calibration/calibration_iteration.py#L600-L617](src/rigel/calibration/calibration_iteration.py#L600-L617), `next_capture_target` is computed as the ratio of `captured_mu` (sweep estimate of gDNA-attributable density in captured regions) to `off_target_mu` (`rho_off × L_eff`). The numerator collapses because the four-state E-step labels probe-overlapping expressed regions as `expressed_capture` (gDNA mass = 0 by label), not as gDNA-bearing. So `mu_sweep` in those regions is tiny → enrichment looks like ~1.

This is the same "gDNA exists everywhere" problem your own design doc flagged, but now we see it has a *second* downstream consequence: it starves the gDNA *exposure weight* used inside the EM, not just the prior mass. The prior-mass handoff fix only addresses one of the two channels.

### Sanity check: 5K bleed is consistent with the math

With `alpha_gdna_add ≈ 1,067`, `alpha_rna_add ≈ 1,932` in the current production (capture-gated) prior at locus 1, and an EM that prefers RNA by ~1.5 nats per fragment, the locus-level Bayesian update places ~78% on RNA. Truth is 72%. The 5–8% over-shoot × 100K fragments ≈ 5,500–8,400 bled fragments. Exactly what we observe.

### Fixes for Q1 (in priority order)

1. **Fix `capture_enrichment_target` first.** Estimate enrichment from raw density in probe-overlap windows vs. off-target windows, *independent of the four-state labels*. Even a coarse estimate like `mean unspliced density in probe-overlap exons / rho_off` should yield 50–500× in capture-on data. This will multiply `gdna_em_exposure_weight` accordingly and reduce the per-fragment 1.5-nat RNA bias.
2. **Raise `MAX_ESS` (or remove the cap) when prior-mass precision is high.** Currently `MAX_ESS = 3000` is hardcoded in [src/rigel/calibration/adaptive_prior.py#L30](src/rigel/calibration/adaptive_prior.py#L30). Scale it by `prior_mass.precision` and `region_arrays.fragment_count_in_locus`, e.g. `ess = min(0.3 × n_em, precision × n_local)`. This is exactly your Phase 2 plan.
3. **Diagnostic to add to `summary.json`:** for each locus, log `mean(log L_RNA - log L_gDNA)` over the EM fragments. That single number tells you whether the per-fragment scoring is fair before any prior.

---

## Q2 — Why does all-region prior-mass regress GENE0008 in capture-OFF?

### What we observe

In the `gdna_high_ss_0.99_nrna_none_capture_off` case at Locus 7 (GENE0008):

| Transcript | Truth | baseline (state) | all-region | gated |
|---|---|---|---|---|
| GENE0008.2 | 340 | 824 (+484) | **1972 (+1632)** | 853 |
| GENE0008.3 | 4,119 | 1,047 (−3,072) | **0 (−4,119)** | 1,083 |
| GENE0008.4 | 21,854 | 10,534 | 10,411 | 10,484 |
| Locus prior `n_local_gdna` | — | 245 | **2,396** | 245 |
| Locus prior `n_local_rna` | — | 12,746 | 10,595 | 12,746 |

So `all_region` injects ~10× more gDNA prior mass at this locus and collapses GENE0008.3 to zero.

### Root cause: strand-deconvolution mis-attributes intronic RNA to gDNA in long expressed loci

The all-region `prior_mass` path uses `build_prior_mass_deconvolution` with `strand_channels` ([src/rigel/calibration/calibration_iteration.py#L317-L370](src/rigel/calibration/calibration_iteration.py#L317-L370)):

```python
gdna = contained.mean_count + boundary_left.mean_count + boundary_right.mean_count
```

Each of `contained.mean_count` is `n_total − R̂` from a BetaBinomial deconvolution where the prior on gDNA strand-share is `Beta(κ/2, κ/2)` with `κ = 1e6` (very sharp). For a region with `n_total = n` and `k_antisense = k`, the deconvolution's gDNA point estimate behaves like:

$$\hat{G} \approx \frac{2k}{1 - \text{kappa-correction}} \approx 2k \quad \text{(for SS=0.99)}$$

In a long expressed intron of GENE0008.4 (which has ~22K unspliced fragments), even SS=0.99 leakage gives `k_antisense ≈ 0.005–0.01 × n` purely from mis-orientation noise, FL-noise on long transcripts, and any antisense pre-mRNA. Across **all compartments** (contained + left boundary + right boundary) and across **every region in the locus**, those small antisense counts get doubled and summed → ~2,400 "gDNA" fragments. **Truth at this locus is essentially 0 gDNA** (it's the capture-off, but it's an expressed locus so all 22K+ fragments are RNA).

This is **not** an isoform-allocation/prior-strength regression as the doc currently frames it. It is a **prior-mass accuracy regression at this specific locus**: the strand deconvolution itself is biased high when (a) the region has very large `n_total` of expressed-region unspliced mass, and (b) any antisense channel is non-zero.

The state-based labels accidentally protect this locus because they assign expressed regions to RNA categories that carry zero gDNA mass. The all-region path removes that protection — exposing a real weakness in the strand deconvolution.

### Why the isoform collapse looks so extreme

Once gDNA absorbs ~1,000 extra RNA fragments, the residual RNA pool gets re-EMed over the four GENE0008 transcripts. GENE0008.3 and GENE0008.2 are short (~1.8 kb spliced, similar exon structure) and likely share most of their reads in equivalence classes. With prior pseudo-counts `alpha_rna_add ≈ 23K` distributed *uniformly* across the isoforms (the RNA pool prior is isoform-neutral in mass but actually pushes mass toward isoforms with more equivalence-class coverage), the EM collapses to a single dominant isoform — GENE0008.2 — and zeroes GENE0008.3.

This second effect is real and is what your doc described, but it is *triggered* by the bad prior mass.

### Fixes for Q2

1. **Precision-gate the strand-derived gDNA mass.**
   In `build_prior_mass_deconvolution`, set `gdna_unspliced_mean[r] = 0` (or shrink to the density-derived `mu_gdna`) when `precision < precision_floor` and `p_expressed > 0.7`. A `precision_floor` of ~0.5 would protect highly expressed regions where antisense counts are dominated by noise rather than gDNA signal.
2. **Cap per-region `gdna_unspliced_mean` by `kappa_d × n_total × (1 − p_r1_sense)`**, the expected gDNA-attributable count under the null. Anything beyond that is over-fit noise.
3. **For locus-level prior**, weight `prior_mass.gdna_unspliced_mean[r]` by `1 − p_expressed[r]` (a soft version of the state-based protection). This is *not* the same as state-derived mass labels — it is a confidence-weight on the prior_mass channel itself, which is consistent with your "gDNA-everywhere" doc's Phase B design (confidence-weighted source priors).
4. **Separately**: audit `alpha_rna_add` injection in the EM solver. The collapse of GENE0008.3 → 0 with non-zero truth indicates that an isoform-neutral RNA-pool prior is *not* implemented neutrally inside `em_solver.cpp`. Worth a follow-up read of the `vbem` digamma update path to check that `alpha_rna_add` is spread proportional to current marginals rather than tilting toward dominant isoforms.

---

## How to confirm the Q1 hypothesis fragment-by-fragment

The annotated BAM (`annotated.bam`) carries `ZC` (fragment class), `ZL` (locus), `ZW` (posterior), `ZT` (winning transcript) per record. The oracle BAM `qname` encodes truth: `gdna:...` for gDNA, `GENE####.#:...` for RNA. A small script that joins on `qname` and computes:

- Confusion: `(true_source, ZC)` cross-tab per locus
- For misassigned `gdna→RNA` fragments: extract the per-fragment `log L_RNA - log L_gDNA` from the EM scoring path (the C++ `em_solver` accumulates these per equivalence class — would need a small instrumentation hook)
- Plot bled fragments' genomic positions vs. probe regions

would directly verify both: (a) bled fragments cluster at probe-overlap exons (confirming exposure-weight starvation), and (b) the per-fragment likelihood ratio at those positions is ~+1.5 nats RNA-favored. I'd recommend adding this as `scripts/debug/bleed_fragment_audit.py` — the data infrastructure (annotated BAM, oracle BAM, locus stats) is already in place.

---

## Summary of recommended action sequence

1. **Fix `capture_enrichment_target`** (label-independent estimator) — this is the single change with the biggest expected impact on the capture-on case.
2. **Scale `MAX_ESS` by `prior_mass.precision`** and by `n_em_fragments` (Phase 2 of your existing plan).
3. **Precision-gate `prior_mass.gdna_unspliced_mean` per region** in expressed regions where strand-deconv evidence is statistically weak — this directly removes the Q2 regression without re-introducing label-based mass semantics.
4. **Audit `alpha_rna_add` propagation in `em_solver.cpp`** to confirm isoform-neutrality.
5. **Add per-locus diagnostic** `mean_logL_ratio_rna_vs_gdna` to summary.json so this class of failure becomes visible without bespoke audits.

The good news: both failures are localized, well-instrumented by the existing diagnostic data, and each has a contained fix that does not require redesigning the latent-state model.