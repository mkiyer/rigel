# Calibration — theory and inference

The authoritative description of Rigel's calibration stage: the generative model, the acyclic
single-pass inference, and the interface to the per-locus EM. (Consolidates the former
`caljointmodel/` theory docs.) For the code map see `docs/ARCHITECTURE.md`; for the scientific
overview see `docs/METHODS.md`; for open work see `CALIBRATION_TODO.md`.

## 1. The problem

An RNA-seq library is a mixture of three molecular origins:
- **mRNA** — mature, spliced transcripts (the quantity of interest).
- **nRNA** — nascent (pre-mRNA), unspliced, intron-retaining, strand-matched to the gene.
- **gDNA** — genomic DNA contamination: unspliced, **unstranded** (double-stranded ⇒ 50/50
  sense/antisense), present genome-wide.

Calibration's job is to estimate the **library hyperparameters** and a **per-locus gDNA-vs-RNA
prior** that the locus EM consumes. It does *not* assign individual fragments — that is the EM.

## 2. Two conditionally-independent clues

Each genomic **region** (a maximal interval of constant transcript-feature signature) and each
**boundary side** between regions is a *node*. A node's unspliced fragment mass is deconvolved into
gDNA vs RNA using two clues that are conditionally independent given the library scalars
`(gdna_density_global, rna_sense_frac)`:

- **Count clue** — *how much* mass is here vs expected gDNA density. gDNA is roughly uniform across
  the genome at density `ρ`; a node whose unspliced mass greatly exceeds `ρ·eff_len` carries RNA.
- **Strand clue** — *which strand*. gDNA is 50/50 (sense rate ½); RNA is stranded (sense rate
  `κ = rna_sense_frac`). A node's observed sense fraction places it between the two.

They are combined per node as a Bayesian posterior over the gDNA fraction `g`:

```
posterior(g) ∝ Beta(g ; count prior)  ·  BetaBinomial(sense | g, κ, overdispersions)
```

When one clue is uninformative the posterior gracefully falls back to the other: unstranded data
(κ→½) ⇒ flat strand likelihood ⇒ count-only; a node with no count evidence ⇒ Jeffreys count prior ⇒
strand-only.

## 3. Count-observability (signature-based, no circularity)

The clues are gated by the **region signature** (4 bits: exon±, intron±). A node is **count-
observable** — i.e. its unspliced mass is gDNA by construction — when no mature RNA can be present:

- **Region** observable ⇔ **no exon bit** (intergenic or intron-only): no mature transcript covers
  it, so contained unspliced mass is gDNA (+ nascent, which the strand clue removes).
- **Boundary** observable ⇔ **no shared exon bit** across its two sides: no single exon continues
  across it, so no unspliced *mature* RNA crosses; the crossing mass is gDNA(+nascent).

Observable nodes anchor the gDNA density **without needing the strand clue or the deconvolution** —
this is what makes the calibration *acyclic* (no density→deconv→density loop).

**Intergenic is a strand wildcard.** A region with no transcript (`TS_NONE`) carries no defined
sense and is compatible with *either* strand. So a gene-edge boundary (a stranded exon abutting
intergenic) is strand-oriented by the gene side — `{POS,NONE}→POS`, `{NEG,NONE}→NEG`. Only an
opposite-strand pairing `{POS,NEG}` or an ambiguous (`TS_AMBIG`) flank leaves the sense undefined.
(Treating gene-edge boundaries as strand-blind was a bug that leaked ~half their pure-gDNA crossing
mass; see CALIBRATION_TODO Issue #3.)

## 4. The acyclic single-pass inference

`calibrate()` runs one feed-forward pass — no EM loop, no convergence test:

```
substrate (per-region 3-view sufficient statistics from the accumulator payload)
  1. RNA strand balance        → κ = rna_sense_frac           (strand_balance.py)
  2. node gDNA density          → per-node density + count-prior mean, strand-cleaned by κ
                                                                (density_model.py)
  3. gDNA & RNA strand overdispersion (Beta-Binomial, pooled MoM + prior shrinkage)
                                                                (gdna_strand.py)
  4. count overdispersion       → per-type NB α (contained / crossing)  (count_dispersion.py)
  5. joint count×strand deconv  → per-node gDNA / RNA mass      (joint_deconv.py)
  6. derive                     → gdna_density_global, geometric gDNA length  (derive.py)
```

### 4.1 Strand balance (κ)
`rna_sense_frac` is the posterior-mean sense fraction of the **spliced** unique-mapper channel
(spliced ⇒ pure RNA). Zero spliced reads ⇒ the orientation is unidentifiable ⇒ `CalibrationStrandError`
(a real RNA-seq library always has spliced reads).

### 4.2 Count clue: per-node gDNA density (`density_model.py`)
The gDNA density is read directly from count-observable nodes and **imputed locally** for the rest
(an exon region anchors from its observable boundary sides; run interiors carry the nearest anchor;
no global sweep). It is **strand-cleaned** by κ first, so the density is *clean gDNA*, not
gDNA+nascent. The strand clean is the linear unmix `ĝ = (s−κ)/(½−κ)`, made robust by precision-
weighted shrinkage (Issue #1). The node returns: the **mean** `count_gdna_frac` (the prior mean)
and the **honest support** `count_support` (the gDNA count behind it).

### 4.3 Strand overdispersion (`gdna_strand.py`)
gDNA is unstranded (mean ½) but real libraries are **overdispersed** about ½; RNA likewise about κ.
Both are fit as Beta-Binomial intra-class correlations by a shared pooled method-of-moments core
(gDNA from count-observable seed regions + boundary sides; RNA from boundary-side spliced counts),
shrunk toward the *same* default prior so that under sparse data both collapse to one distribution
and unstranded data is uninformative.

### 4.4 Count overdispersion (`count_dispersion.py`)
RNA-seq/gDNA counts are NB-overdispersed, so the count-prior **concentration** is not the raw count
(Poisson) but the **overdispersion-limited effective count** `N_eff = N/(1+α·N)` (→ 1/α for large
N). `α` is fit per count-type (contained regions vs crossing boundary sides) by NB MoM against the
global-density expectation, shrunk toward the pooled-seed trend (geometric `α₀=1` fallback). This
replaced an earlier categorical "defer to strand" zeroing with one principled precision.

### 4.5 Joint deconvolution (`joint_deconv.py`)
Per node, the count prior `Beta(N_eff·gf, N_eff·(1−gf))` (Jeffreys-floored) is multiplied by the
strand Beta-Binomial likelihood; the reported gDNA fraction is the posterior median, optionally
shifted in log-odds by the FP-aversion knob `gdna_strand_llr_bias` (default 0). gDNA mass = `g·M`,
RNA mass = `(1−g)·M` + the deterministic spliced mass. Mass is conserved per node.

### 4.6 Derive (`derive.py`)
`gdna_density_global` (a QC scalar) and the geometric gDNA length are derived from the aggregate
deconvolved masses — not looped — which is what dissolved the old EM-loop's sparse-data collapse.

## 5. Interface to the EM: the per-locus prior (`priors.py`)

`assemble_priors` turns the per-region `CalibrationResult` into **two per-locus Dirichlet scalars**
plus the gDNA component's effective length:

- `gdna_prior_count` — Σ over the locus of the deconvolved gDNA mass (boundary-crossing gDNA first
  re-attributed to its origin region by length-bias-free flux transport).
- `rna_prior_count` — the deconvolved RNA mass.
- `gdna_eff_len` — the **inverse participation ratio** of the per-region gDNA mass (Laplace-
  smoothed), so under capture the gDNA support contracts to the exons and the gDNA component
  competes at its true local density.

The prior's only job is to set each locus's gDNA-vs-RNA split; the EM distributes the RNA mass among
the compatible transcripts. The locus EM has **`n_t + 1` components** — one per transcript row
(annotated mRNA + synthetic nRNA spans) plus one gDNA component.

## 6. Design invariants

- **Acyclic:** observable nodes anchor the density strand-free; the global density is *derived* from
  the aggregate, never fed back. No density→deconv→density loop.
- **Graceful degeneracy:** zero gDNA, zero spliced (errors loudly), unstranded (count-only), sparse
  seeds (priors govern) are all valid, bounded outputs — not failures.
- **Conditional independence:** count and strand are combined as independent clues given the library
  scalars; a node's ~1/N self-contribution to the global density is negligible (empirical Bayes).
- **Mass conservation:** every node's gDNA+RNA mass equals its input mass.

## 7. Open work

Tracked in `CALIBRATION_TODO.md`. Active design docs: `strand_cleaning_robustness_design.md` +
`strand_clean_global_target_design.md` (Issue #1), `count_overdispersion_integration_plan.md`
(Issue #2), `density_phase2_dna_fraction_design.md` (the deferred DNA-fraction lever).
