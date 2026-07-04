# Diagnosis: captured-exon gDNA peak smoothing (Initiative #1)

_2026-07-04. Branch `calib-gdna-accuracy`. The highest-leverage item in
`effective_length_state_and_roadmap.md` §6: the deconvolution's captured-exon gDNA accuracy, whose ~2×
under-read blocks the principled enriched-mode (`kmeans_midpoint`) effective-length reference._

## The instrument

The sim oracle BAM encodes each fragment's true origin + span in its read name. **gDNA read-names are
genomic** (they match the alignment), so per-region **true** gDNA contained mass is exact. (RNA read-names
are **transcript-relative** — do NOT map them as genomic; that was a false start.) Comparing observed
(calibrated) vs true per-region gDNA gives a direct per-node accuracy readout. Tooling:
`scripts/debug/calib_gdna_smoothing.py`, `oracle_gdna_density.py`, `oracle_calibration.py` (+ the
`RIGEL_ORACLE_CALIB` hook in `calibrate.py`).

## The finding: it's peak SMOOTHING, not mass loss

Target (gdna300, ss0.99, capture-on):

- **Total gDNA mass is accurate** (obs 0.998× true); per-node density correlation 0.928.
- **But the density PEAKS are smoothed**: peak obs/true = 0.53×. The top-5 hottest true nodes lose half
  their mass (0.50×); top-20 → 0.65×; top-50 → 0.85×. Mass flows from the sharp peaks down to the
  moderately-enriched shoulder — the enriched half's total stays 1.00× (pure smoothing).
- **Intergenic gDNA is recovered perfectly** (peak 1.00×, corr 1.000) — no RNA there, no ambiguity.
  **Exon peaks are worst** (0.65×). So it is a gDNA-vs-RNA **split** error on high-density exon nodes.
- **Strand-dependent**: exon peak 0.65× stranded vs **0.23× unstranded** (ss0.50); total 0.998× vs 0.903×.
  Unstranded, gDNA (50/50) and unstranded RNA (50/50) are indistinguishable by strand, so the split
  collapses. Capture-off "peaks" are just Poisson noise on a uniform field (correctly not tracked).

## The trace: two stacked causes

Calibration is a **two-pass** solve (`calibrate.py`): PASS-1 `_sweep(None)` → train a nonparametric KDE
density-mixture prior (`GdnaDensityPrior`) on the pass-1 nodes → PASS-2 `_sweep(gdna_prior)`. Ablating
Phase-2 (force the training substrate empty):

| | exon-peak obs/true | total gDNA |
|---|---|---|
| PASS-1 only (no Phase-2 prior) | 0.80× | 0.918× |
| FINAL (2-pass, + Phase-2 KDE prior) | 0.65× | 0.998× |

1. **PASS-1 (strand-solve) under-calls the enriched-exon gDNA fraction** — peaks → 0.80×, total → 0.918×;
   severe unstranded (no strand signal to separate gDNA from RNA).
2. **The Phase-2 KDE density prior *restores* the total (0.918→0.998×) but *smooths* the peaks (0.80→0.65×)**
   — it pulls high-density enriched exon nodes toward the (lower) population density modes. It is the
   **primary peak-smoother**: its regularization, meant for under-determined thin/unstranded nodes, also
   drags down confident high-density exon peaks.

## Fix direction (next cycle)

1. **Peak-preserving density prior** — down-weight the KDE prior where the local evidence is strong /
   high-density, so it regularizes the depleted/thin tail but leaves confident enriched peaks intact.
2. **Strand-independent enriched-gDNA signal in pass-1** — e.g. the component FL model, or anchoring on the
   exon↔intron boundary-crossing gDNA, so the unstranded case does not collapse.

Restoring the exon peaks makes the enriched-mode (`kmeans_midpoint`) reference read correctly, unlocking the
co-design: `kmeans` as the production reference + a milder, more accurate effective-length contraction.
