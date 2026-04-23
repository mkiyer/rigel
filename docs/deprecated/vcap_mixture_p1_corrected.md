# VCaP P0 + P1 — Corrected Analysis

**Date:** 2026-04-22 (evening revision)
**Status:** Supersedes `vcap_mixture_p0_findings.md` for anything
related to per-locus gdna_rate. See §0 "Retraction" below.

## 0. Retraction

The earlier P0 analysis reported the mega-locus `gdna_rate = 0.500`
and took that as evidence of a "symmetric OVR-EM identifiability
trap". **That was wrong.** It was a reporting artifact of a
**double-count bug in `AbundanceEstimator.get_locus_df`** (now fixed):
the per-locus `total` column computed `r["mrna"] + nrna + r["gdna"]`,
but `r["mrna"]` from the C++ solver *already includes* synthetic nRNA
mass (see `em_solver.cpp::assign_posteriors`: "Transcript (mRNA or
synthetic nRNA — both handled the same way) ... `mrna_total += p`").

Corrected mega-locus numbers at dna80m:

| quantity | value | note |
|---|---:|---|
| `n_em_fragments` | 10,538,621 | |
| `r["mrna"]` (includes nRNA) | 4,201,472 | |
| `nrna` (synthetic only) | 2,134,678 | |
| `mrna_pure = r["mrna"] − nrna` | 2,066,794 | |
| `gdna` | 6,337,149 | |
| true_total = mrna + gdna | 10,538,621 | ✓ matches n_em_fragments |
| **true_rate = gdna / true_total** | **0.6013** | |
| reported rate (double-count) | 0.5000 | was the bug |

Headline summary numbers (`quantification.gdna_fraction`) were NOT
affected — those use `quant_df["count"].sum()` which excludes synthetic
transcripts, so they were computed correctly. The bug was purely in
`loci.feather`.

## 1. Docstring fixes

Two misleading docstrings corrected:
- `src/rigel/calibration/_em.py` — removed "over exonic regions";
  calibration EM operates on the full genome partition (intergenic,
  intronic, exonic). Only the hard-splice RNA anchor can fire on
  exonic regions.
- `src/rigel/locus.py::compute_locus_priors` — now explicitly
  documents what γ is (expected gDNA counts / observed total counts
  over regions overlapping the locus), and how it's used downstream
  (both as Dirichlet pseudocount and as EM warm-start ratio for θ).

## 2. Accounting bug fix

`AbundanceEstimator.get_locus_df` now:
- reports `mrna` column as **annotated-only** (`r["mrna"] − nrna`)
- reports `total = mrna + nrna + gdna = n_em_fragments` (no overlap)
- reports `gdna_rate = gdna / total` (physically meaningful)

Golden tests regenerated (21 scenarios); all 98 tests pass.

## 3. Corrected deficit attribution at dna80m

Expected gDNA fraction in total VCaP mixture: 0.800.
Summary-level reported: `quantification.gdna_fraction = 0.666`.
Locus-rollup (denominator excludes intergenic): `0.685`.

Per-size-bin with true accounting:

| bin | n_loci | frags (M) | true_rate | nrna_frac | share of deficit |
|---|---:|---:|---:|---:|---:|
| <100 kb | 22,030 | 39.6 | 0.708 | 0.055 | 42% |
| 100 k–1 M | 3,334 | 24.5 | 0.687 | 0.097 | 32% |
| 1–10 M | 31 | 0.7 | 0.608 | 0.144 | 2% |
| ≥10 M (mega) | 1 | 10.5 | 0.601 | 0.203 | 24% |

**Key pattern:** gDNA under-estimation is *consistent* across all
locus sizes (all bins at 0.60–0.71 vs expected ~0.78–0.82), not
concentrated at the mega-locus. nRNA fraction grows monotonically
with locus size (5% → 20%), which correlates with the siphon
hypothesis but does not single it out as the sole driver at small
loci.

## 4. Revised top-down story

The deficit is **spread across loci of all sizes** with a magnitude
of ~10–20 percentage points below expected. The mega-locus is not
qualitatively different from the 100 kb–1 M bin — it just carries
more fragments and therefore more absolute deficit.

Candidate mechanisms, not mutually exclusive:

1. **nRNA absorbing unstranded unspliced fragments** (the "siphon").
   Correlates with nrna_frac growing by size bin. Confirmable by
   oracle-BAM fragment tracking.
2. **`region_e_gdna` clipping bias.** `_calibrate.py:91`:
   `region_e_gdna = min(λ_G × mappable_bp, n_total)`. For pure-gDNA
   regions, Poisson noise makes `n_total < λ_G × mappable_bp` roughly
   half the time; the clip systematically pulls down per-region
   expected gDNA, which flows to lower per-locus α_gDNA and thus a
   lower warm-start θ_gDNA. Cheap to test: drop the clip, rerun,
   compare.
3. **Global λ_G under-estimate.** The calibration EM classifies
   regions into expressed / not-expressed; if the expressed class
   absorbs some "low-μ" regions that are actually pure gDNA, the
   estimator of λ_G from the RNA side is diluted. Harder to test
   without ground truth per-region classifications.

## 5. Next steps I recommend

In descending ROI:

**A. Test hypothesis #2 first — remove the `min(λ_G·mbp, n_total)`
clip.** One-line change in `_calibrate.py`. No expected production
regression risk since `region_e_gdna` only feeds the per-locus prior
(warm-start), not any EM posterior directly. Rerun 8 mixtures and
see what the mega-locus and per-bin true_rate do.

**B. Investigate hypothesis #1 (nRNA siphon) via oracle BAM.** For
dna80m, pull fragments where the oracle says "gDNA" and check what
rigel assigns. Quantifies the siphon directly.

**C. Calibration robustness check.** Per your request: pull
high-γ regions from the v4 EM (confident not-expressed), recompute
λ_G from just those with a simple pooled MLE, and compare to the
EM-estimated λ_G. Also compute per-region observed/expected ratio
to measure Poisson overdispersion in the "pure gDNA" class. This
is a standalone script, no production code change.

**Not recommended now:** mega-locus EM instrumentation. The
mega-locus is not anomalous; fixing it won't move the needle
disproportionately.

## 6. Files changed

- `src/rigel/calibration/_em.py` — docstring only
- `src/rigel/locus.py` — docstring only
- `src/rigel/estimator.py::get_locus_df` — accounting fix
- `tests/golden/*_loci_df.tsv` — 21 files regenerated
