# VCaP RNA/DNA Mixture — gDNA Underestimation Analysis

**Date:** 2026-04-22
**Scope:** Analysis of 8 in-silico VCaP RNA + exome-DNA mixtures
(`/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human/mctp_vcap_rna20m_dna*m/`)
and root-cause hypotheses for Rigel's gDNA underestimation at high
contamination.

## 1. The data

- 20 M starting RNA fragments (post-pipeline-filter) from VCaP cell line
  (real RNA-seq, assumed near-pure).
- Pure exome DNA mixed in at 1, 2, 5, 10, 20, 40, 80 M fragments.
- Mixture ratios span 0.10× → 4.0× RNA-equivalent (expected gDNA fraction
  0% → 80%).
- Ground-truth *transcripts* are unknown (real data), but ground-truth
  *mixture proportions* are known.

## 2. Observed calibration + quantification (Rigel v0.4.0, VBEM mode)

| Run | expected gDNA | cal_gdna_frac | λ_G (per bp) | quant_gdna_frac | quant_nrna_frac |
|---|---:|---:|---:|---:|---:|
| rna20+dna00 | 0.000 | 0.0005 | 1.0e-5 | 0.007 | 0.025 |
| rna20+dna01 | 0.048 | 0.0069 | 9.4e-5 | 0.056 | 0.030 |
| rna20+dna02 | 0.091 | 0.0125 | 1.8e-4 | 0.102 | 0.033 |
| rna20+dna05 | 0.200 | 0.0250 | 4.2e-4 | 0.209 | 0.041 |
| rna20+dna10 | 0.333 | 0.0384 | 8.0e-4 | 0.330 | 0.050 |
| rna20+dna20 | 0.500 | 0.0535 | 1.5e-3 | 0.463 | 0.063 |
| rna20+dna40 | 0.667 | 0.0671 | 2.9e-3 | 0.582 | 0.076 |
| rna20+dna80 | 0.800 | 0.0774 | 5.4e-3 | 0.666 | 0.088 |

## 3. Corrections to the initial read of these numbers

- **`cal_gdna_fraction` is not comparable to the mixture ratio.** It is
  the expected gDNA fraction *within exonic calibration regions*, where
  the per-bp RNA rate is 10–1000× the gDNA rate. Comparing 0.077 to
  0.80 is not the right test.
- The relevant global quantities are **λ_G (fragments per mappable exonic
  bp)** and the downstream **per-locus γ_ℓ prior**.
- `λ_G` scales near-linearly with input gDNA up through dna20; doubling
  input DNA (20→40→80 M) yields λ_G ratios of 1.91× and 1.85×. The
  global rate estimator is approximately unbiased.

## 4. Where gDNA is actually being lost

### 4.1 Reporting bug (P0) — now fixed

The `gdna_prior` column of `loci.feather` was stuck at 0.0 for every
locus in every run. Root cause:

- `_build_locus_meta` in `pipeline.py` writes keys `alpha_gdna` and
  `alpha_rna` per locus.
- `AbundanceEstimator.get_locus_df()` in `estimator.py` was reading a
  stale key `gdna_init` that nobody sets.

**Fix** computes `gdna_prior = α_gDNA / (α_gDNA + α_RNA)` and also
emits the raw `alpha_gdna` / `alpha_rna` columns. No change to any
solver math — only to what the feather file reports. This unblocks
every downstream per-locus analysis of the calibration prior.

### 4.2 The mega-locus converges to γ = 0.5 regardless of mixture

Size-stratified per-locus `gdna_rate` (fraction of locus fragments
assigned to gDNA) across mixtures:

| Run | <100 kb | 100 kb – 1 Mb | 1–10 Mb | ≥10 Mb (mega) |
|---|---:|---:|---:|---:|
| dna00m | 0.009 | 0.004 | 0.008 | — |
| dna05m | 0.201 | 0.233 | 0.240 | 0.319 |
| dna10m | 0.319 | 0.350 | 0.403 | 0.274 |
| dna20m | 0.455 | 0.461 | 0.432 | **0.410** |
| dna40m | 0.579 | 0.557 | 0.500 | **0.467** |
| dna80m | 0.671 | 0.626 | 0.532 | **0.500** |

At dna80m the single mega-locus (569 Mb span, 70 k transcripts) carries
10.5 M fragments (14% of the total) but classifies only 50% as gDNA.
This single locus accounts for ≈1.9 M of the ~11 M fragments missing
from the overall gDNA count at dna80m. 0.500 is the degenerate
symmetric equilibrium reached when the per-locus Dirichlet prior has
α_gdna ≈ α_rna and the unspliced-fragment likelihoods are
near-degenerate across mRNA/nRNA/gDNA.

### 4.3 nRNA siphon (secondary)

nRNA fraction rises monotonically from 2.5% (pure RNA) → 8.8% (dna80m)
and is 13–17% inside mega-loci. This is the NTA1/TA1 identifiability
issue already catalogued in the repo. It is **not** the dominant
failure mode at high gDNA — it accounts for perhaps 3–5 pp of the 14 pp
underestimate at dna80m, vs the ~9 pp attributable to the mega-locus
prior issue.

### 4.4 Why the exonic-only calibration cannot see genomic gDNA load

The v4 calibration EM fits a single global λ_G from ~400 k exonic
reference regions with mappable bp. At dna80m, fit parameters are
μ_R = −4.13, σ_R = 1.86, giving a LogNormal on μ whose 25th percentile
≈ 4.6e-3 — **numerically identical to the true λ_G at dna80m**. The
two mixture components overlap; the count-LLR channel is near-zero;
γ_ℓ is driven almost entirely by the strand-LLR channel and by π_soft.

For mega-loci whose footprint is ≳99% intronic, `compute_locus_priors`
still computes γ_ℓ from **exonic** regions overlapping the locus —
regions that are RNA-dominated. The result: α_gDNA ≪ α_RNA, and the
per-locus EM has no informative prior to pull it off the symmetric
0.5 equilibrium. The denominators are mismatched: we estimate "gDNA
per exonic bp" and apply it to a predominantly-intronic footprint.

## 5. Proposed directions (not yet committed — awaiting P0 re-run)

P0 → **Fix reporting bug** (done this session).

P1 → **Separate genomic/intronic calibration channel for λ_G.** The
exonic v4 EM stays as-is for gDNA-vs-RNA disambiguation in exonic bp.
For the *mega-locus prior* (and potentially for the global λ_G that
feeds it), estimate a pure "genomic λ_G" from intergenic fragments
(already counted by the scanner: `fragment_stats.intergenic`) and
deep-intronic windows of single-isoform genes. Compare ratio to
exonic-mixture λ_G as a bias diagnostic.

P2 → **Span-weighted per-locus γ_ℓ.** Replace the current
count-weighted formula
`γ_ℓ = Σ E[gDNA]_i / Σ n_total_i` (over exonic regions overlapping the
locus) with a span-aware version
`γ_ℓ = λ_G_gen · gdna_span / (λ_G_gen · gdna_span + E[RNA in locus])`.
This respects that gDNA generates fragments proportional to the full
genomic footprint while RNA generates only from exonic bp.

P3 → **Quantitative nRNA siphon measurement.** Using the indexed
annotated CRAMs, cross-check the dna00m vs dna80m nRNA-assigned
fragments at the same genomic positions. Grows with mixture →
siphoned gDNA; stays flat → true nRNA.

P4 → **Per-region calibration diagnostics TSV** (one row per
calibration region: n_u, n_s, E, γ_final, LLR_count, LLR_strand,
LLR_fl). Trivially cheap to emit; future debugging is then a pandas
query instead of an instrumented rebuild.

## 6. Next step

Rerun all 8 mixtures with the P0 fix in place, then inspect:
- Distribution of `gdna_prior` across loci and its correlation with
  observed `gdna_rate`.
- Specifically: is `gdna_prior` in the mega-locus truly near 0 (as
  predicted by section 4.4), or already near the expected mixture
  value?

The answer to that single question determines whether P1+P2 are
actually needed or whether the EM itself is the culprit.
