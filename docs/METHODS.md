# Rigel: Bayesian RNA-seq Transcript Quantification with Joint mRNA, Nascent RNA, and Genomic DNA Deconvolution

**Matthew K. Iyer**

---

The scientific method behind Rigel's joint mRNA / nascent-RNA / gDNA quantification. This document
follows the code; for the code map see `ARCHITECTURE.md`, and for the calibration internals see
`calibration/CALIBRATION_ARCHITECTURE.md`.

## 1. Model

A sequenced library is a mixture of three molecular origins:

- **mRNA** — mature spliced transcripts; the abundances we want.
- **nascent RNA (nRNA)** — pre-mRNA: unspliced, intron-retaining, strand-matched to its gene.
- **gDNA** — genomic DNA contamination: unspliced, **unstranded** (50/50 sense/antisense),
  distributed genome-wide (and, under hybrid capture, enriched on probed exons).

Rigel estimates per-transcript abundances (mRNA and nRNA spans alike) and the per-locus gDNA
contamination jointly, rather than discarding multi-mappers or ignoring contamination.

## 2. Reference index

`rigel index` builds a region partition and interval trees from a GTF + FASTA. The genome is cut
into **regions** — maximal intervals of constant transcript-feature *signature* (a 4-bit code:
exon±, intron±). Regions are the unit of calibration; transcripts (including unique nRNA spans,
materialized as synthetic transcript rows) are the unit of quantification.

## 3. Single-pass scan and fragment resolution

One C++ htslib pass over the name-sorted BAM: group reads into fragments, resolve each fragment
against the interval tree (which transcripts/regions it is compatible with), train the strand and
fragment-length (FL) models from unique mappers, buffer the resolved fragments, and deposit
fractional fragment mass into the per-region/per-boundary **accumulator** (the calibration
sufficient statistics).

## 4. Fragment likelihood scoring

Each fragment is scored against every compatible component (transcripts + a gDNA component) by a
likelihood combining:
- **fragment length** — RNA vs gDNA FL pmfs (empirical-Bayes smoothed), queried at the fragment's
  genomic span;
- **strand** — sense/antisense consistency with the transcript orientation vs gDNA's 50/50;
- **coverage / position** — within the transcript's effective length;
- **splice** — junction-spanning fragments are compatible only with the spliced transcript (gDNA,
  being genomic, cannot span a junction).

## 5. Calibration: gDNA-vs-RNA deconvolution

Before quantification, the calibrator deconvolves each genomic node's **unspliced** fragment mass
into the 2-simplex **(f_rna₊, f_rna₋, f_g)** — sense-RNA / antisense-RNA / gDNA. Calibration models
only the **RNA-vs-gDNA** split; the nascent-vs-mature separation is left to the per-locus EM
downstream.

**The count-zero-information principle.** A fragment count carries *no intrinsic information* about a
node's gDNA/RNA composition: 50 sense / 50 antisense unspliced fragments are equally consistent with
pure gDNA (genomic strand is symmetric) and pure RNA in an unstranded library. The count enters only
as **statistical power** — the precision of an estimate — never as a vote on composition. A node's
composition is therefore set by exactly three sources:

- the **strand likelihood** — the Beta-Binomial tilt of the per-strand counts (RNA tilts the sense
  rate toward the library sense fraction κ; gDNA stays at ½). This is the only *intrinsic* gDNA/RNA
  signal, and the count enters it only as overdispersed Fisher information (statistical power),
  vanishing to zero at κ=½ (unstranded) or low count;
- **cross-node imputation** — density messages from neighbouring nodes at the belief-free Poisson
  disagreement variance σ²_imp (fit once). gDNA flows genomically; per-strand RNA flows only across
  edges where that strand is transcriptionally continuous (the transcript-structure gate);
- the **global gDNA prior** — the population baseline gDNA density (plus a trained Phase-2 gDNA-density
  KDE) at a robust MAD-spread precision, which sways nodes the strand and imputation leave undetermined.

The genome is modelled as a **bipartite region↔boundary node chain** — regions and boundaries are both
first-class nodes owning fragment mass, an effective length, and a composition pie; boundaries
additionally own the one-sided, motif-stranded **spliced** RNA (a fixed mature-RNA floor). Inference is
a **single forward-backward belief-propagation sweep** over the chain (L→R then R→L), which is exact
because the chain is a forest of linear paths — there is no iterative fixed-point loop. σ²_imp and every
global-prior input are fit once *before* the sweep.

The result is the library hyperparameters (`gdna_density_global`, `rna_sense_frac` = κ, the gDNA and
RNA strand Beta-Binomial overdispersions) plus the per-region / per-boundary-side deconvolved gDNA/RNA
mass, assembled downstream into a **per-locus gDNA-vs-RNA prior**. Full theory: the
count-zero-information architecture in `calibration/CALIBRATION_ARCHITECTURE.md`.

## 6. Per-locus EM

Loci are connected components of overlapping transcripts/regions. Each locus is an independent
subproblem with **`n_t + 1` components** — one per transcript row (annotated mRNA + synthetic nRNA
spans) plus one gDNA component. The calibration prior enters as **two per-locus Dirichlet
pseudocounts** (`gdna_prior_count`, `rna_prior_count`) plus the gDNA component's effective length;
the EM distributes the RNA mass among the compatible transcripts by expectation-maximization
(VBEM-capable; OpenMP-parallel across loci; Kahan-compensated summation).

## 7. Outputs

Per-transcript abundance (counts, TPM) for mRNA and nRNA spans, per-locus gDNA contamination, and
library QC scalars (gDNA fraction, strand specificity, FL means). An optional annotated BAM tags
each fragment with its assignment for downstream auditing.

## 8. Validation

A synthetic simulator (`rigel sim`) generates ground-truth mixtures across gDNA level × strand
specificity × hybrid-capture, with a net-fragment-flow evaluation that separates *systematic*
deconvolution bias from *unrecoverable* sequence-identical ambiguity. See `BENCHMARKING.md` and the
`calibration-benchmark` skill.
