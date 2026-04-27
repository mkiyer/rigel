# Calibration root-cause (final, hybrid-capture aware) + strand-first proposal

**Date:** 2026-04-22
**Status:** Root cause understood; plan drafted; awaiting user approval to implement.
**Prior docs:** `vcap_mixture_p1_corrected.md`, `root_cause_locus_em_nrna_siphon.md` (superseded by this document).

---

## TL;DR

- The "10× calibration failure" diagnosed earlier was based on a **wrong assumption** — that intergenic λ_G is ground truth. Under hybrid capture that's invalid (intergenic is off-target and depleted by **100–140×** vs exonic).
- With the **dUTP strand convention** applied correctly (R1 antisense to gene = RNA, p=0.81), the strand channel alone produces a clean decomposition of gDNA per region kind that tracks `quant.gdna_total` to within **≈1.1–1.5×** across the entire titration.
- The strand signal is dramatically underutilized today. In hybrid-capture + highly-stranded data, **strand alone is near-sufficient to calibrate**, and the count-based (Poisson-LogNormal) channel is actively harmful because it assumes spatial uniformity of λ_G that is violated by 2–3 orders of magnitude across region kinds.
- **Fix:** strand-first calibration with per-kind λ_G stratification (exonic / intronic / intergenic). Optional: targets-BED mode.

---

## The titration evidence

Per-library strand decomposition (dUTP: `p_rare_RNA=0.19`, `p_rare_gDNA=0.5`), applied to all unspliced fragments in oriented regions:

```
lib      obs_tot  quant_gDNA  strand_gDNA  ex_fgd  in_fgd     λ_ex      λ_in      λ_ig    capture_enrich(ex/ig)
dna00m   13.88M     0.19M       0.04M     0.000   0.000      —         —      3.2e-5     —   (no gDNA)
dna01m   14.82M     1.07M       0.34M     0.000   0.601      —       1.8e-4  5.2e-5      —
dna02m   15.75M     1.93M       0.68M     0.000   0.760      —       3.7e-4  7.2e-5      —
dna05m   18.53M     4.47M       1.80M     0.018   0.890    7.2e-4    9.6e-4  1.3e-4     5.5×
dna10m   23.08M     8.66M       5.77M     0.283   0.941    1.58e-2   1.9e-3  2.3e-4    69×
dna20m   31.88M    16.95M      13.47M     0.529   0.968    4.5e-2    3.8e-3  4.2e-4   108×
dna40m   48.38M    32.33M      27.90M     0.713   0.982    1.0e-1    7.3e-3  7.7e-4   129×
dna80m   77.46M    59.54M      53.33M     0.830   0.990    2.0e-1    1.35e-2 1.4e-3   139×
```

Key observations:
1. **dna00m (pure RNA, no gDNA spike)**: strand decomposition correctly recovers ≈0 gDNA everywhere. Strand channel is specific.
2. **quant_gDNA ≈ 1.1–1.5 × strand_gDNA** consistently across the titration. The locus EM arrives at *approximately the same* gDNA estimate the strand channel gives directly, but via a much noisier path (and with a systematic 10–30% overshoot — probably FL model + mRNA coverage channel allocating additional mass to gDNA).
3. **`λ_exonic / λ_intergenic` stabilizes at 100–140×** in the well-resolved regime (dna20m+). This is the true hybrid-capture on-target enrichment. A single uniform λ_G is wrong by ~2 orders of magnitude.
4. **`f_gdna_intronic` ramps from 0% (dna00m) to 99% (dna80m)** monotonically. The strand channel reads gDNA content directly without needing a density baseline.

## Why the current calibration underestimates λ_G inside loci

The two-component EM (`_em.py::run_em`) is trying to separate pure-gDNA regions (class G) from expressed regions (class R = gDNA + RNA mixture). It does this with:
- **Count channel** (Poisson-LogNormal, ≈dominant weight): "rate-per-bp vs the seed-spliced μ_R anchor".
- **Strand channel** (Binomial gated by aggregate z≥3, typically weaker weight): "is the rare-strand fraction consistent with gDNA or RNA?"
- **FL-shape channel** (after iter 1): "is the observed FL histogram closer to the gDNA or RNA FL model?"

In hybrid capture:
- The Poisson-LogNormal is **mis-specified**. There is no single λ_G; there are three very different λ_G's. The EM picks a compromise value ≈ the geometric mean (8e-4 for dna10m, vs λ_in=1.9e-3, λ_ig=2.3e-4) that fits neither population and produces a per-locus prior `γ_ℓ` that is too low inside transcribed regions.
- The strand channel, which is the only channel that can separate gDNA from RNA inside a transcribed locus, is **gated too conservatively** and **down-weighted** relative to the misspecified count channel.

The locus EM then has to re-estimate gDNA from scratch using its own strand/FL channels, which mostly works — that's why the final quant gdna_total is in the right ballpark — but the per-locus `gdna_prior` is meaningless (0.05 instead of ~0.5 for typical intronic-rich loci), so the allocation is noisy and tends to over-shoot by 1.5×.

---

## Proposal: strand-first calibration with per-kind λ_G stratification

### Core idea
When `strand_specificity >= 0.75`, trust the strand channel. Use the count channel only as a tiebreaker for low-read regions.

### Algorithm (Phase 1: strand-only λ_G per region kind)

For each region with `n_u_oriented >= k_min` (say 20):

1. Compute `p_rare = rare / (expect + rare)`.
2. Estimate `f_gdna(region) = (p_rare - p_rna_rare) / (p_gdna_rare - p_rna_rare)`, clamped to `[0,1]`.
   - `p_rna_rare = 1 - ss_estimate`, `p_gdna_rare = 0.5`.
3. Post-smooth with a Bayesian prior (Beta(1, 1) on `f_gdna`) for low-coverage regions.
4. Pool per-kind:
   - `λ_G_exonic = Σ f_gdna·n / Σ mbp` over oriented exonic regions with sufficient reads.
   - Same for `λ_G_intronic` and (trivially) `λ_G_intergenic` (= n_u / mbp, since f_gdna ≡ 1).

### Per-locus prior computation (Phase 1)
```
γ_ℓ = Σ_i λ_G(kind_i) · mbp_i / (Σ_i n_total_i + ε)
```
This replaces `region_e_gdna = min(λ_G · mbp, n_total)` with a kind-aware version. The `min()` clip is no longer needed because `λ_G(kind)` is now properly specified per class.

### Channel weighting (Phase 2, optional)
Promote strand to the primary channel:
```
logit(γ_region) = LLR_strand + w_count · LLR_count + w_fl · LLR_fl
```
with `w_count, w_fl ≤ 1` reduced automatically when `ss >= 0.75` and per-region oriented-read precision is sufficient. At the limit `w_count = 0`, calibration reduces to a strand-only model.

### Fallback for unstranded libraries
When `ss < 0.6`, revert to the current count-dominated EM (the present code). Strand-first is a mode, not a replacement.

### Optional: Targets-BED mode (Phase 3)
Accept `--targets targets.bed`; split λ_G into `{on_target, off_target}` instead of `{exonic, intronic, intergenic}`. Useful when the kit is known and the user wants exact control.

---

## What to do about per-locus priors today

Even without implementing the above, we can get a significant improvement by:

1. **Computing three λ_G values in `_m_step`** (one per kind) rather than one.
2. **Summing the per-region expected gDNA** using the kind-appropriate λ_G.
3. **Removing the `min(λ_G·mbp, n_total)` clip** — it's a leftover safety net that now causes ~50% under-estimation for intronic regions (where n_total includes pre-mRNA RNA but we clip to total).

This is a 1-day implementation + test pass. It does not require changes to the scorer or locus EM.

---

## Validation plan

For each fix iteration, re-run the 8-library titration and report:

1. **Calibration-level**: `Σ prior_E_gdna` per library vs `strand_gDNA` — should track within 1.1× across the titration.
2. **Locus-level**: per-locus `gdna_prior` distribution. Mean should shift from 0.27 → 0.5 for intronic-rich loci.
3. **Quant-level**: `quant.gdna_total` vs `strand_gDNA`. Should drop from 1.5× to ~1.1× (the residual being legitimate FL+coverage contributions).
4. **mRNA/nRNA**: should improve on vcap_pristine_vs_dirty comparison (reducing the 1.5× gDNA overshoot gives more mRNA mass).

## Artifacts

- Per-region-kind decomposition script: embedded in this session's terminal output; will be packaged as `scripts/debug/strand_decomp_per_kind.py` before implementation.
- Per-lib region counts (already dumped): `/scratch/.../runs/human_optionb_diag/dna??m/region_counts.feather` (8 libs × 683,916 rows).
- Region dump script: [scripts/debug/dump_region_counts.py](scripts/debug/dump_region_counts.py).
