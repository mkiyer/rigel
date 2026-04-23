# VCaP RNA/DNA Mixture — P0 Analysis Findings

**Date:** 2026-04-22
**Scope:** Results after the P0 reporting fix: `get_locus_df` now emits
the actual per-locus Dirichlet prior `gdna_prior = α_gDNA / (α_gDNA +
α_RNA)` plus `alpha_gdna` / `alpha_rna`.
**Data:** `/scratch/mkiyer_root/mkiyer0/shared_data/rigel/runs/vcap_mixture_p0fix/`
(8 mixtures rerun with `rigel quant --em-mode vbem --threads 10` on
`annotated.bam` from the original runs).

## 1. P0 is a pure reporting fix — quant unchanged

Summary-level numbers from rerun match the original hulkrna runs bit-for-bit
(modulo RNG noise in sample-mode EM). No calibration or quantification
math was touched. See [docs/calibration/vcap_mixture_gdna_underestimation.md](vcap_mixture_gdna_underestimation.md)
for the original hypothesis before I could see the per-locus priors.

## 2. Priors look sensible in aggregate

Per-locus prior distribution scales monotonically with the true mixture:

| mixture | expected gDNA | prior_mean | prior_med | prior_p90 | prior_max | #priors>0 |
|---|---:|---:|---:|---:|---:|---:|
| dna00m | 0.000 | 0.014 | 0.001 | 0.020 | 1.00 | 15864 / 15864 |
| dna01m | 0.048 | 0.103 | 0.022 | 0.257 | 1.00 | 21777 / 21777 |
| dna02m | 0.091 | 0.144 | 0.034 | 0.495 | 1.00 | 22997 / 22997 |
| dna05m | 0.200 | 0.212 | 0.053 | 0.920 | 1.00 | 24798 / 24798 |
| dna10m | 0.333 | 0.266 | 0.069 | 1.000 | 1.00 | 26091 / 26091 |
| dna20m | 0.500 | 0.291 | 0.077 | 1.000 | 1.00 | 25059 / 25059 |
| dna40m | 0.667 | 0.332 | 0.090 | 1.000 | 1.00 | 25404 / 25404 |
| dna80m | 0.800 | 0.362 | 0.102 | 1.000 | 1.00 | 25396 / 25396 |

So the prior "reporting bug" was **purely cosmetic** — every locus
always had a non-zero prior. That was my first hypothesis and it was
wrong. The prior distribution is roughly lognormal-looking: a long
heavy body at low values with a spike at 1.0 from single-region loci
overlapping pure-gDNA regions. Median is much lower than mean.

## 3. The EM is stronger than the prior

Priors are globally ~2× smaller than the truth (mean 0.36 at dna80m vs
expected 0.80), but the likelihood moves the posterior well above the
prior in most bins. At dna80m, fragment-weighted diagnostics by prior
decile (loci with ≥10 EM fragments):

| prior_bin     | n_loci | frags | prior_med | rate_med | **rate_agg** |
|---|---:|---:|---:|---:|---:|
| 0 – 0.01 | 1321 | 7.6 M | 0.007 | 0.845 | **0.586** |
| 0.01 – 0.05 | 6983 | 32.1 M | 0.027 | 0.775 | **0.641** |
| 0.05 – 0.10 | 4249 | 16.5 M | 0.070 | 0.761 | **0.669** |
| 0.10 – 0.25 | 3385 | 18.3 M | 0.145 | 0.814 | **0.582** |
| 0.25 – 0.50 | 1026 | 0.7 M | 0.330 | 0.859 | **0.769** |
| 0.50 – 0.75 | 401 | 0.05 M | 0.615 | 0.829 | **0.764** |
| 0.75 – 1.00 | 2527 | 0.12 M | 1.000 | 0.896 | **0.865** |

Loci with prior ≤ 0.01 still average 59% gDNA; loci with prior ≥ 0.75
hit 87%. The likelihood is clearly dominant — raising the prior will
probably not help proportionally.

## 4. Where the 14 pp undercount at dna80m actually comes from

Deficit attribution by locus-size bin (dna80m, total deficit = 0.134):

| bin | n_loci | frags | obs_rate | true-rate<br>estimate | **share of deficit** |
|---|---:|---:|---:|---:|---:|
| <100 kb | 22030 | 39.6 M | 0.671 | ~0.69¹ | ~0–38% |
| 100 k – 1 M | 3334 | 24.5 M | 0.626 | ~0.75 | **33%** |
| 1 – 10 M | 31 | 0.7 M | 0.531 | ~0.75 | 2% |
| ≥10 M (mega) | 1 | 10.5 M | **0.500** | ~0.78 | **27%** |

¹ Under a uniform-per-bp gDNA model, the <100 kb bin should have
a *lower* expected gDNA rate than the overall 0.80 because RNA
concentrates there. The observed 0.671 is therefore probably close
to truth.

The **mega-locus is the biggest single source** of the deficit, but it
explains only ~27% of the gap. The remaining ~70% comes from the
**100 kb – 1 M bin** and from nRNA siphon.

## 5. The mega-locus is a likelihood-degeneracy failure, not a prior problem

The 569 Mb mega-locus at dna80m:

| metric | value |
|---|---:|
| locus_span | 569,789,373 bp |
| n_transcripts | 70,539 |
| n_genes | 7,674 |
| n_nrna_entities | 36,695 |
| n_em_fragments | 10,538,621 |
| **gdna_prior** | **0.1145** |
| α_gDNA | 0.572 |
| α_RNA | 4.428 |
| **observed gdna_rate** | **0.500** |

Key observation: the EM **moved from prior 0.11 to posterior 0.50** —
by 39 percentage points — and stopped there. With 10.5 M fragments
against 5.0 total Dirichlet pseudocounts, the prior is numerically
negligible; the 0.500 equilibrium is a property of the likelihood, not
the prior. **Raising `gdna_prior_c_base` from 5 to any reasonable value
will not fix this.**

0.500 is the characteristic equilibrium of a 2-way competition (gDNA
vs "collapsed RNA class"). With 70 k transcripts + 36 k nRNA entities
and VBEM sparsification, the RNA side collapses into a few survivors
sharing ≈50% of the mass; gDNA — the only OVR-class component with
non-shared prior mass — absorbs the other 50%. This is an
**OVR-EM identifiability failure**, not a calibration failure.

## 6. nRNA siphon is real, grows with mixture, and accounts for ~30–40% of the deficit

nRNA fraction per mixture per size-bin:

| mixture | overall | <100 kb | 100 k – 1 M | 1 – 10 M | ≥10 M |
|---|---:|---:|---:|---:|---:|
| dna00m (pure RNA) | 0.027 | 0.023 | 0.034 | 0.036 | — |
| dna05m | 0.042 | 0.033 | 0.055 | 0.085 | 0.137 |
| dna10m | 0.050 | 0.038 | 0.065 | 0.119 | 0.099 |
| dna20m | 0.062 | 0.043 | 0.075 | 0.125 | 0.140 |
| dna40m | 0.073 | 0.048 | 0.085 | 0.131 | 0.155 |
| dna80m | 0.083 | 0.052 | 0.089 | 0.126 | 0.168 |

The dna00m row is the "baseline" nRNA rate (2–4%) under pure RNA.
Everything above that is siphoned gDNA. Excess nRNA at dna80m per bin:

| bin | excess nRNA frac | × bin frags | siphoned (M) |
|---|---:|---:|---:|
| <100 kb | 0.029 | 39.6 M | 1.15 M |
| 100 k – 1 M | 0.055 | 24.5 M | 1.35 M |
| 1 – 10 M | 0.090 | 0.7 M | 0.06 M |
| ≥10 M | 0.130 | 10.5 M | 1.37 M |
| **total** |  |  | **~3.9 M** |

The total missing gDNA at dna80m is ≈ 0.134 × 77 M ≈ 10.3 M fragments.
The nRNA siphon accounts for ~40% of that (≈3.9 M). The remainder is
split between mega-locus EM degeneracy (~2.5 M) and modest
under-estimation in 100 kb – 1 M loci (~3–4 M).

## 7. Revised diagnosis

The initial hypothesis (prior zeroing + mega-locus prior collapse) was
wrong. The actual failure modes, in order of contribution:

1. **nRNA siphon (~40% of deficit).** Structural: unspliced
   unstranded fragments that should go to gDNA are absorbed by
   nRNA entities (which share genomic footprint with their parent
   mRNA and are unspliced by construction). Grows monotonically
   with gDNA load.

2. **Mega-locus EM degeneracy (~25% of deficit).** 70 k
   transcripts, 36 k nRNAs, 10.5 M fragments. Likelihood dominates,
   OVR-EM converges to symmetric ≈0.500 regardless of prior.

3. **Moderate-locus (100 kb – 1 M) underestimation (~30% of
   deficit).** Smaller per-locus effect but spread across thousands
   of loci. Likely a mix of (1) and (2) at smaller scale.

None of these is "the global λ_G is wrong." The calibration module's
global rate is fine. What's wrong is the per-locus EM's ability to
*use* that rate when:
- an alternative explanation exists (nRNA)
- the likelihood is flat over many candidates (mega-loci)

## 8. Next-step priorities (revised)

**Prior-focused ideas (P1, P2 from previous doc) are downgraded** —
they'd help only the ~25% mega-locus share, and only if we can first
break the symmetric-0.5 equilibrium.

**Promoted priorities:**

A. **Investigate why unspliced unstranded fragments siphon into nRNA.**
   Direct inquiry: pull random dna80m fragments tagged as nRNA from
   the annotated CRAM, check how many coincide with dna00m nRNA
   fragments at the same position. The fraction that does *not* is
   the siphon.

B. **Break the OVR-EM symmetric-0.5 trap in mega-loci.** Options:
   - Pre-assign unspliced-fragment evidence to gDNA more aggressively
     when the locus likelihood is flat (identifiability-aware prior
     boost).
   - Split mega-loci by sub-chromosome-arm connectivity and run each
     independently.
   - Cap VBEM sparsification intensity per locus.

C. **Add per-region calibration diagnostics** (still a good idea; now
   tagged P4-equivalent) — but this is diagnostic, not a fix.

The P0 fix itself is landed and required no solver changes. The next
substantive work should focus on (A) and (B).
