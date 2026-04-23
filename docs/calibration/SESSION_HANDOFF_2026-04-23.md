# Session handoff: calibration channel-fusion & NB density model

**Date:** 2026-04-23
**Status:** Active design. No code changes yet. Ready to implement pending user sign-off on three questions at the end of this doc.

This file hands off the state of the calibration redesign work so a fresh session can pick up without reloading context.

---

## 1. Where we are

### 1.1 The problem (confirmed)

Rigel's calibration EM ([src/rigel/calibration/_em.py](../src/rigel/calibration/_em.py)) fuses three channels — **density, strand, gDNA-FL** — by **unit-weight LLR summation** in `_e_step`. On hybrid-capture RNA-seq (VCaP MCTP titration), this fails: the density channel is **confidently wrong** because it assumes a single global per-bp gDNA rate `λ_G`, but hybrid-capture data has 100–140× spatial variation between on-target and off-target regions.

- The strand channel works correctly (SS ≈ 0.81, clean gDNA-vs-RNA separation via 0.5 vs 0.19 rare-fraction means).
- The density channel produces per-region LLRs of ±20+ nats (mis-specified), while strand LLRs are bounded by `n · KL(0.5||ss) ≈ 0.23·n` nats. Unit-weight summation effectively down-weights the correct channel.

**Result:** `quant.gdna_total / strand_inferred_gdna ≈ 1.1–1.5×` across VCaP titration; expressed loci mis-classified as class R when strand clearly identifies them as G.

### 1.2 What was ruled out

- **nRNA siphon as root cause** — wrong, based on assumption that intergenic = ground truth. Hybrid capture violates that.
- **Per-kind λ_G stratification** — rejected: "hybrid capture panels are arbitrary and variable".
- **Arbitrary SS ≥ 0.75 threshold** — rejected: no magic numbers.
- **The historical `W = (2·SS−1)²` weighting in [docs/METHODS.md](METHODS.md) §5.1** — verified NOT in current code (removed during v3→v4 rewrite). Retains no authority.

### 1.3 Current design documents

- **[docs/calibration/channel_fusion_hybrid_capture.md](calibration/channel_fusion_hybrid_capture.md)** — 7-section research doc enumerating problem, current implementation, failure analysis, and 4 solution classes (A/B/C/D).
- **[docs/calibration/density_nb_model_plan.md](calibration/density_nb_model_plan.md)** — First-draft implementation plan for Negative Binomial density channel. **Needs revision per §3 below.**
- Superseded: `calibration/strand_first_hybrid_capture_plan.md`, `calibration/root_cause_locus_em_nrna_siphon.md`.

---

## 2. Decisions locked in

From user direction on the 4 open questions in `channel_fusion_hybrid_capture.md` §6:

| # | Decision |
|---|---|
| 1 | **NB density channel — YES**, primary direction. Borrow theory from DESeq2. |
| 2 | **Targets-BED mode — YES** but **deferred to Phase 5**. Required for unstranded hybrid capture. Separate design document to be written later. |
| 3 | **Strand-anchored λ_G re-initialization — NO.** Not pursuing. |
| 4 | **Retire `(2·SS−1)²` from METHODS.md — YES, but only AFTER NB density is validated.** Gate at Phase 4 of rollout. |

---

## 3. Open revision to `density_nb_model_plan.md` (pending user sign-off)

The initial plan put NB on **both** G and R classes with a shared `φ_G`. User flagged feeling it was "somewhat overengineered" and wanted open-ended thinking on R.

I (agent) conceded the simplification in the last exchange. **Proposed revised plan:**

### 3.1 Proposed simplifications

1. **NB on G class only.** R class stays as Poisson-LogNormal (unchanged). Rationale:
   - The G class only had Poisson for flexibility — nothing else to absorb cross-region rate variance. Adding `φ_G` is the minimum fix.
   - The R class already has unlimited overdispersion via σ_R² in the LogNormal prior. Poisson-LogNormal is a compound distribution; `σ_R² → ∞` can match any NB variance. "R is not broken because R cannot be broken."
   - The symmetry argument (both G and R as Gamma-distributed processes) is aesthetically nice but doesn't change behavior — LogNormal puts negligible mass at μ→0 so the μ=0 consistency boundary is never evaluated.
   - **Practical win:** the 21-node Gauss-HermiteE quadrature stays Poisson. No NB PMF in the inner loop. Same per-iter cost as today.

2. **Fit `φ_G` on the hard-G subset** (γ_i > 0.95 over eligible), not γ-weighted over all regions. Rationale:
   - Parallels existing design: μ_R, σ_R are already fit only on hard-spliced R-anchors (`k^s > 0 & eligible`). Making G fit symmetric in spirit.
   - Insulates `φ_G` from R-class heterogeneity leakage. If φ_G were fit γ-weighted over all regions, residual σ_R² tail variance could inflate φ_G, over-shrinking density LLRs even when density should be informative.
   - Fit φ_G represents something cleanly interpretable: "across regions we believe are gDNA-only, what is the cross-region CV² of the per-bp rate?"

3. **Phase 1 regression test becomes trivial.** With `nb_density=False`, full fallback to current code. With `nb_density=True` on synthetic data with <100 hard-G regions, clamp `φ_G=0` → identical to Poisson. No reliance on φ_G "happening to converge near zero".

### 3.2 Numerical cautions (locked in from user review)

- **Poisson branch threshold `_EPS = 1e-4`** (not machine ε). Avoids catastrophic cancellation in `r·log(r/(r+μ))` when `r ≫ μ`. At `φ=1e-4`, NB is biologically indistinguishable from Poisson.
- **`r·log1p(-μ/(r+μ))`** form preferred for numerical safety on the boundary.
- **Warm-started golden section**: bracket `[max(0, φ_prev/2), min(50, 2·φ_prev)]` after iter 0, with fallback to full range if optimum hits endpoint. Drops ~60 → ~20 evaluations per iter.

### 3.3 Additional operational polish

- **Skip φ_G update in iter 0.** No classification → no hard-G subset → no safe fit. Set `phi_G=0` for iter 0, start fitting at iter 1+. Dovetails with FL channel's existing "disabled in iter 0" pattern.
- **Hard-G subset size guard.** If `|{γ_i > 0.95}| < 100`, clamp φ_G=0 and log warning. Protects against weird inputs; doesn't affect typical libraries (VCaP has 100k+ candidate regions).
- **Surface `phi_G` AND empirical CV² of `k_i/E_i` over hard-G subset** in `summary.json`. Disagreement between fitted and empirical values signals different mis-specification (e.g., mappability bias).

### 3.4 Three questions currently open for user

1. **Accept the simplification: NB on G only, R unchanged?** (agent leans yes)
2. **Hard-G subset (γ_i > 0.95) for φ_G fitting?** vs γ-weighted over all eligible? (agent leans hard-G)
3. **Any open-ended R-class concerns?** User flagged wanting to think about this. Agent's current view: R is fine, σ_R² is the load-bearing flexibility there, leave it alone. User may have specific thoughts on LogNormal-vs-alternative as the expression prior.

---

## 4. Where to resume

**Immediate next step in the new session:**

1. User responds to the three questions in §3.4.
2. Agent revises [docs/calibration/density_nb_model_plan.md](calibration/density_nb_model_plan.md) to reflect the simplified design (NB on G only, hard-G subset fit, numerical cautions, iter-0 guard).
3. User approves revised plan.
4. **Phase 1 implementation begins**: `_nb_logpmf`, `_count_llr_nb_g_only` (new name reflecting scope), `_estimate_phi_G_cr` with hard-G subset, add `nb_density` flag to `_config.py`, plumb `phi_G` through `_types.py` + `_summary.py`. Default `nb_density=False` so all existing tests pass unchanged.

**Then phases 2–5 per the plan:**
- Phase 2: VCaP titration validation (8 libs at `/scratch/.../runs/human_optionb_diag/dna*m/`). Primary metrics: φ_G vs gDNA spike monotone; per-region density LLR shrinkage; γ vs strand-only correlation; quant.gdna_total vs strand-inferred ratio.
- Phase 3: Flip `nb_density=True` as default. Regenerate goldens once.
- Phase 4: Update [METHODS.md](METHODS.md), remove `(2·SS−1)²` weighting section.
- Phase 5: Separate work — `--targets BED` mode for unstranded hybrid-capture support.

---

## 5. Key code locations (for fast re-orientation)

- Fusion: [src/rigel/calibration/_em.py](../src/rigel/calibration/_em.py) `_e_step` (line ~530), currently `log_odds = log_prior_odds + llr_count + llr_strand + llr_fl_soft` (unit weights).
- M-step: same file, `_m_step` (line ~598). λ_G weighted Poisson MLE; μ_R, σ_R from hard-spliced R anchors.
- Density kernel to replace: `_count_llr_poisson_ln` (line ~134).
- Strand kernels (unchanged): `_aggregate_strand_z` (~178), `_strand_llr` (~210), `_strand_llr_betabinom` (~295), `_estimate_kappa_marginal` (~400).
- FL kernel (unchanged): `_fl_llr` (~421).
- Orchestration: [src/rigel/calibration/_calibrate.py](../src/rigel/calibration/_calibrate.py); line ~91 has `region_e_gdna = minimum(fit.lam_G * mappable_bp, n_total)` clip.
- Types to extend: [src/rigel/calibration/_types.py](../src/rigel/calibration/_types.py).
- Config to extend: [src/rigel/calibration/_config.py](../src/rigel/calibration/_config.py).
- Summary emitter: [src/rigel/calibration/_summary.py](../src/rigel/calibration/_summary.py).

---

## 6. Data ready for validation

- **VCaP titration BAMs & quant outputs:** `/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human/mctp_vcap_rna20m_dna{00,01,02,05,10,20,40,80}m/`.
- **Per-region diagnostic dumps (8 libs):** `/scratch/.../runs/human_optionb_diag/dna{00,01,02,05,10,20,40,80}m/{region_counts,fl_table}.feather + scan_meta.json`.
- **Dump script:** [scripts/debug/dump_region_counts.py](../scripts/debug/dump_region_counts.py) — works, documented.
- **Rigel index used:** `/scratch/.../refs/human/rigel_index`, `regions.feather` has 683,916 rows, 2.86 Gb mappable.
- **Env:** `source /home/mkiyer/sw/miniforge3/etc/profile.d/conda.sh && conda activate rigel`. Build: `pip install --no-build-isolation -e .`. Tests currently: 1065/1065 pass.

---

## 7. One-paragraph recap for the fresh session

Rigel's three-channel calibration EM (density / strand / gDNA-FL) fuses per-region LLRs with unit weights, which fails on hybrid-capture RNA-seq because the density channel assumes a single global `λ_G` that is violated by ~100× capture-enrichment targeting, producing confident-wrong LLRs that overwhelm the correct strand channel. We've decided to fix it by replacing the G-class Poisson with a Negative Binomial (Gamma-Poisson compound) parameterized DESeq2-style with a shared dispersion `φ_G`, estimated via Cox-Reid-corrected MLE. Latest refinement: apply NB to G only (keep R as Poisson-LogNormal — σ_R² already absorbs any R overdispersion), fit φ_G on the hard-G subset (γ_i > 0.95), with numerical guards (Poisson branch at φ ≤ 1e-4, warm-started golden section, iter-0 skip, small-subset clamp). Targets-BED mode deferred to Phase 5 for unstranded hybrid capture. Next action in the new session: user answers three questions in §3.4, agent revises `density_nb_model_plan.md`, Phase 1 implementation begins.
