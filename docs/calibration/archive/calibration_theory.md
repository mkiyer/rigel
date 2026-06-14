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

The two clues are **decoupled** and combined by a **precision-weighted gradient**, not a product —
because the strand estimator is *unbiased* (gDNA symmetric at ½, RNA at κ) while the count estimator
is *biased* under hybrid capture, so mixing them by a product re-introduces bias:

```
g = w·g_strand + (1−w)·g_count,   w = I/(I+I₀),  I = N·(2κ−1)²,  I₀ = 10
   g_strand : Beta-Binomial posterior over g, weak Beta(½,½) prior — strand DIRECTION, capture-invariant
   g_count  : clip(ρ·eff_len / mass) — count MAGNITUDE, on raw contained + strand-cleaned boundaries
```

The weight `w = I/(I+I₀)` is a **per-node strand-trust gradient**: `I = N·(2κ−1)²` is the carried strand
**information** — the per-fragment discriminability `(2κ−1)²` (the strand Fisher information / squared
standardized separation) times the node's unspliced count `N`, i.e. the total information the node's own
reads carry about `g`. `I₀` is the half-trust scale. `w → 1` when the strand is genuinely informative
(high κ **and** enough depth — trust the node's own, capture-invariant deconvolution); `w → 0` at `κ=½`,
thin `N`, or `AMBIG`/no-shared-strand (defer to the count imputation). A gradient, **no gate**.

**Why `I/(I+I₀)` and not the κ-only `(2κ−1)²`** — the strand estimate has `SE(g) = 1/√I`, so it is only
meaningful at `I ≫ 1`. The κ-only effect size ignores `N`: it both *under*-trusts a confident strand
(0.99 → w=0.96, dinging an excellent model) and, fatally, *over*-trusts a near-½ strand at high depth
(under capture, `N` huge ⇒ even `κ=0.501` makes `I` non-trivial ⇒ a garbage `g_strand`, `SE≈1`, gets
weight → manufactured gDNA). `I/(I+I₀)` with **`I₀=10`** sets half-trust at `SE≈0.3` (`I=11`): the strand
is trusted only once its own information is real. `g_strand` (direction) and `g_count` (magnitude) use
orthogonal aspects of the data and `g_count`'s *contained* counts are raw (the strand enters once, via
`g_strand`; only the **boundary** crossings are strand-cleaned, for the exon/AMBIG imputation), so the
combine does **not** double-count. The retired joint product is archived in
`archive/joint_deconvolution.md`; the strand-first redesign that arrived at this per-node weight is in
`redesign_phase3_plan.md`.

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
  1. RNA strand balance         → κ = rna_sense_frac                          (strand_balance.py)
  2. node gDNA density (RAW)     → seed for the strand overdispersion fit       (density_model.py)
  3. gDNA & RNA strand overdispersion (Beta-Binomial, pooled MoM + prior)      (gdna_strand.py)
  4. strand_deconvolve → cleaned_gdna_count: clean the BOUNDARY crossings      (strand_deconv.py)
  5. node gDNA density (RAW contained + CLEANED boundaries) → g_count          (density_model.py)
     + splice-junction 3-term gDNA-FRACTION upgrade                            (splice_junction.py)
  6. per-node gradient combine g=w·g_strand+(1−w)·g_count, w=I/(I+I₀)          (strand_deconv.py)
  7. derive                     → gdna_density_global, geometric gDNA length    (derive.py)
```

### 4.1 Strand balance (κ)
`rna_sense_frac` (κ) is the posterior-mean sense fraction of the **spliced** unique-mapper channel
(spliced ⇒ pure RNA). It enters the per-node combine weight `w = I/(I+I₀)` via the per-fragment
discriminability `(2κ−1)²` in `I = N·(2κ−1)²` (§2): there is **no** hard identifiability gate — an
unstranded library has κ≈½ ⇒ I≈0 ⇒ w≈0 ⇒ count governs (and, unlike the κ-only effect size, depth `N`
cannot rescue a near-½ strand — see §2 on the knife-edge). Zero spliced reads ⇒ κ undefined ⇒ not an
RNA-seq library ⇒ `CalibrationStrandError`.
(`strand_summary.strand_contrast_identifiable` survives only as a pipeline QC *warning*, not a gate.)

### 4.2 Count module: per-node gDNA density (`density_model.py`)
The gDNA density is read directly from count-observable nodes and **imputed locally** for the rest (an
exon region anchors from its observable boundary sides; run interiors carry the nearest anchor; a
no-anchor region takes the global count-weighted-mean observable density). It operates on whatever
per-view counts it is handed (`gdna_counts`): the live combine passes **RAW contained + STRAND-CLEANED
boundaries** (step 5) — so an exon/AMBIG region imputed from its crossings inherits a nascent-free
density, while the contained own-density of count-observable regions stays raw (the strand enters the
*region* only via `g_strand`). The node returns the count module's gDNA fraction
`count_gdna_frac = clip(density·eff_len / mass)`, used as `g_count` in the combine and as the gDNA
strand-fit seed weight. Improving this estimate under capture (the point-5 unspliced-fraction
projection) is tracked in `count_channel_capture_design.md`.

### 4.3 Strand overdispersion (`gdna_strand.py`)
gDNA is unstranded (mean ½) but real libraries are **overdispersed** about ½; RNA likewise about κ.
Both are fit as Beta-Binomial intra-class correlations by a shared pooled method-of-moments core
(gDNA from count-observable seed regions + boundary sides; RNA from boundary-side spliced counts),
shrunk toward the *same* default prior so that under sparse data both collapse to one distribution
and unstranded data is uninformative. These parameterise the strand module.

### 4.4 Per-node strand/count combine (`strand_deconv.py`)
Each node (§2): `g = w·g_strand + (1−w)·g_count`, `w = I/(I+I₀)`, `I = N·(2κ−1)²`, `I₀=10` — the per-node
strand-trust gradient (`deconv_regions`/`deconv_sides` → `_deconv_per_node`). `g_strand` is the **strand
module** — the Beta-Binomial posterior over `g` (weak Beta(½,½) prior × the strand likelihood) median;
the exact, clip-free robust strand deconvolution (its MLE is the linear unmix `(s−κ)/(½−κ)` but bounded
and overdispersion-widened — no clip bias). `g_count` is the **count module**'s `count_gdna_frac`,
computed on **raw contained + strand-cleaned boundary** counts (`cleaned_gdna_count` removes the
strand-identified RNA from the crossings before the exon/AMBIG imputation — the nascent the count clue
can't otherwise see; the contained stays raw so the strand enters once, via `g_strand`). A
strand-unobservable node (κ=½ / no defined sense / `AMBIG`) is count-only (`w=0`). The blended point
estimate `center` is then
read at the **FP-rate quantile** `g(q) = clip(center + Φ⁻¹(q)·σ)` (the FP-aversion knob
`gdna_deconv_quantile`, default ½ ⇒ no shift), where `σ` is the combined per-node posterior std
(`√(w²·σ²_strand + (1−w)²·σ²_count)`, or `σ_count` when count-routed). The shift is uncertainty-aware
(wider posterior ⇒ larger shift) and widening-only — it never sharpens the bias-robust blend (§4.5).
gDNA mass = `g·M`, RNA mass = `(1−g)·M` + the deterministic spliced mass. Mass is conserved per node.

### 4.5 Count posterior variance (Phase 1) — the quantile width
The count module also returns a per-node variance `σ²_count = μ²·v_rel` capped at the Bernoulli maximum
`μ(1−μ)`. Every count is **Poisson** at baseline: an observable region/boundary side uses its own count
floor `v_rel = 1/N`; an imputed region (exon/AMBIG) uses a non-parametric variance~mean LOESS over its
two anchors, floored by the anchors' Poisson noise. It feeds **only** the quantile width above — never
the combine weight `w = I/(I+I₀)` (which derives from the *strand* information `I`, not the count σ):
under hybrid capture the count σ is *anti-calibrated* (confident = confidently-biased), so it is
trustworthy for *widening* but not for *sharpening* (`docs/calibration/phase2_design.md`).

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

Tracked in `CALIBRATION_TODO.md`. Active design docs: `count_channel_capture_design.md` (the
count-channel direction under capture) + `phase0_phase1_implementation_plan.md` (the shipped
overdispersion teardown + own-mean floor); `strand_clean_robust_deferred.md` (deferred strand-clean
concept); `density_phase2_dna_fraction_design.md` (the deferred DNA-fraction lever).
