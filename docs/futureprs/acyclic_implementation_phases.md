# Acyclic Calibration — Implementation Phases

The calibration is being rebuilt as a **single-pass, acyclic deconvolution**
(`acyclic_deconvolution_design.md`). Each node (region or boundary) is split into RNA and
gDNA from **local** evidence; the global gDNA density `ρ₀` and the per-region exposure are
**derived afterward**. This removes the `ρ₀ → decode → ρ₀` feedback loop that caused the
sparse-data collapse, and dissolves the exposure-dispersion prior.

This document is the phase map. Each phase is built and validated on its own before the
next is layered on.

---

## Phase 1 — Density model (the "count clue")  ← *current*

**Goal.** A per-node gDNA estimate from the **observed counts** in gDNA-representative
nodes, propagated to the rest — never from the global `ρ₀·ω·L`.
**Produces.** Per-node gDNA mass/density from: contained fragments in seed regions
(intergenic/intronic) and crossing fragments at seed boundaries (exon-intron,
exon-intergenic), read through the gDNA FL model; imputed into exons.
**Validates.** On uniform scenarios: `ρ₀` is **bounded** (no 5e-5 collapse) and recovers
in the right range. (Count-only is upper-biased where introns are transcribed — the strand
clue corrects this in Phase 3; that is expected, not a bug.)
**Tears down.** Nothing yet — additive, standalone module + validation script.

## Phase 2 — Strand model (the "BB clue")

**Goal.** A per-node gDNA-vs-RNA split from the sense/antisense imbalance (Beta-Binomial,
regularized so thin counts stay diffuse — the review's Risk A).
**Produces.** Per-node gDNA fraction + its uncertainty, for non-AMBIG nodes.
**Validates.** Correctly calls expressed exons as RNA; diffuse (not over-committed) at low
counts.
**Tears down.** Nothing yet — refactors the existing `_llr_strand` / `fit_strand_balance`
into a standalone decoder.

## Phase 3 — Joint decode (+ the confidence knob)

**Goal.** Integrate the two clues per node, using whichever are valid: count-only for
AMBIG, strand-only for exon-exon boundaries, both elsewhere.
**Produces.** Per-node (RNA mass, gDNA mass) + the single continuous confidence/quantile
knob (specificity↔sensitivity, default = unbiased).
**Validates.** Sparse toy, dense-uniform, dense-capture (incl. silent-but-captured exons).
**Tears down.** The old combined E-step (`estep.estep_view`, `_llr_count` over `μ_g`).

## Phase 4 — Impute + derive `ρ₀` / exposure

**Goal.** Guarantee every node has an answer; then compute the library summaries.
**Produces.** Imputation of non-decodable nodes (neighbors → baseline fallback);
`ρ₀ = ΣĜ/ΣL`; exposure = (region density ÷ global). The unit-information shrinkage of the
*derived* exposure toward 1 (the review's Risk C — protects downstream `gdna_eff_len`).
**Validates.** Every node decoded; exposures finite and sane.
**Tears down.** `update_rho_0` (→ derived aggregate); `exposure.py`'s Gamma posterior;
`update_exposure_dispersion`; `exposure_prior_pseudocount`.

## Phase 5 — Measure stability, regularize only if needed

**Goal.** *Empirically* assess exposure stability across the three regimes. **No
regularization is designed up front.** Only if/where the derived exposures wobble do we add
a simple post-hoc shrinkage — decided by evidence, not in advance.
**Produces.** A stability report; a shrinkage only if warranted.
**Validates.** Graceful degradation sparse→dense; capture enrichment preserved.

## Phase 6 — Integrate + clean up

**Goal.** Wire the new calibrator into `calibrate.py` / the pipeline; remove the dead
feedback machinery; regenerate goldens.
**Produces.** A feed-forward calibration: clues → per-node decode → derive → prior.
**Validates.** Full scenario + golden suites; the 3 nrna_dc failures resolved.
**Tears down.** The outer EM loop in `calibrate.py`; any remaining `μ_g`-feedback paths.

---

**Invariant across all phases:** each node is decoded from its own local data; `ρ₀` and
exposure are *outputs*; nothing global flows back into a node's decode. That one property
is what makes the whole thing stable.
