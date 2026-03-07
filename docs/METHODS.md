# Rigel: Bayesian RNA-seq Transcript Quantification with Joint mRNA, Nascent RNA, and Genomic DNA Deconvolution

**Matthew K. Iyer**

---

## Abstract

We present Rigel, a Bayesian transcript quantification method for RNA-seq
data that jointly models three nucleic acid species present in a single
library: mature messenger RNA (mRNA), nascent pre-mRNA (nRNA), and
genomic DNA contamination (gDNA). Rigel introduces a linked kinetic model
that couples mRNA and nRNA abundances through a shared per-transcript
parameter and a nascent fraction grounded in the steady-state kinetics of
transcription, splicing, and degradation. A hierarchical empirical Bayes
framework initializes per-transcript nascent fractions and per-locus gDNA
rates from strand specificity and genomic density signals. The method
operates through a two-stage single-pass pipeline: a C++ BAM scanner that
resolves fragments against a reference index and trains strand and
fragment-length models, followed by locus-level Expectation-Maximization
(EM) with SQUAREM acceleration. We describe the complete mathematical
framework, the fragment likelihood model, the hierarchical prior system,
convergence properties, and the output specification.

---

## 1. Introduction

### 1.1 Motivation

RNA-seq measures transcript abundances by sequencing fragmented cDNA
libraries. The standard quantification problem—estimating how many
fragments originate from each transcript—has been addressed by
methods including Cufflinks (Trapnell et al., 2010), RSEM (Li and Dewey,
2011), Salmon (Patro et al., 2017), and kallisto (Bray et al., 2016).
These methods model fragments as originating from mature (spliced)
transcripts and assign ambiguous fragments via expectation-maximization
or equivalent optimization.

However, total RNA-seq libraries and many stranded library preparations
capture additional nucleic acid species beyond spliced mRNA:

1. **Nascent RNA (nRNA / pre-mRNA):** Unspliced precursor transcripts
   captured mid-transcription. These contain intronic sequence and
   represent the active transcriptional output of the cell. Nascent RNA
   is abundant in total RNA-seq, nuclear RNA-seq, and chromatin-associated
   RNA-seq protocols.

2. **Genomic DNA (gDNA):** Residual genomic DNA contamination that
   survives DNase treatment during library preparation. gDNA fragments
   are unstranded and uniformly distributed across the genome.

Existing quantifiers do not model nRNA or gDNA, leading to systematic
biases: intronic fragments from pre-mRNA inflate apparent expression of
long transcripts, while gDNA contamination creates false-positive
expression signals particularly for long genes with extensive intronic
sequence.

### 1.2 Approach

Rigel addresses these limitations through three innovations:

1. **Three-species generative model.** Each fragment in the library
   originates from one of three species—mRNA, nRNA, or gDNA—with
   species-specific strand, splice, and length distributions.

2. **Linked kinetic model.** Rather than treating mRNA and nRNA as
   independent parameters, Rigel couples them through a shared
   per-transcript total abundance $\theta_t$ and a per-transcript
   nascent fraction $\beta_t$, motivated by the steady-state kinetics
   of the RNA lifecycle.

3. **Hierarchical empirical Bayes initialization.** Per-transcript
   nascent fractions and per-locus gDNA rates are initialized through
   multi-level shrinkage hierarchies that borrow strength across
   transcripts, loci, chromosomes, and the genome.

### 1.3 Notation

Throughout this document we use the following notation:

| Symbol | Description |
|--------|-------------|
| $T$ | Number of transcripts in a locus |
| $G$ | Number of genes |
| $L$ | Number of loci (connected components of overlapping transcripts) |
| $F$ | Total number of fragments |
| $\theta_t$ | Total abundance (mRNA + nRNA) for transcript $t$ on the simplex |
| $\beta_t$ | Nascent RNA fraction for transcript $t$, $\beta_t \in [\varepsilon, 1-\varepsilon]$ |
| $\gamma_\ell$ | gDNA fraction for locus $\ell$ |
| $\text{SS}$ | Strand specificity, $\text{SS} \in [0.5, 1.0]$ |
| $p_{\text{r1s}}$ | Probability that read 1 aligns sense to the gene |
| $\kappa$ | Shrinkage pseudo-count (effective sample size) |
| $\alpha$ | Dirichlet prior pseudo-count per component |
| $\gamma_{\text{ovr}}$ | One Virtual Read (OVR) prior scale factor |

---

## 2. Generative Model

### 2.1 Three-Species Mixture

We model each fragment $f$ in the library as drawn from a mixture of
three species at each genomic locus $\ell$ containing $T$ transcripts:

$$
p(f \mid \Theta) = \sum_{t=1}^{T} \left[
  \theta_t (1 - \beta_t) \cdot \frac{\ell_{\text{mRNA}}(f, t)}{K_{t,\text{mRNA}}}
  + \theta_t \beta_t \cdot \frac{\ell_{\text{nRNA}}(f, t)}{K_{t,\text{nRNA}}}
\right]
+ \gamma_\ell \cdot \frac{\ell_{\text{gDNA}}(f)}{K_{\text{gDNA}}}
$$

where:
- $\theta_t$ is the total abundance of transcript $t$ on the simplex
  (with $\sum_t \theta_t + \gamma_\ell = 1$)
- $\beta_t$ is the nascent fraction of transcript $t$
- $\ell_{\text{mRNA}}(f, t)$, $\ell_{\text{nRNA}}(f, t)$,
  $\ell_{\text{gDNA}}(f)$ are per-component likelihoods
- $K_{t,c}$ are normalization constants (effective lengths)

This defines a **2T + 1 component** mixture model per locus:
$T$ mRNA components, $T$ nRNA components, and 1 gDNA component.

### 2.2 The Linked Kinetic Model

The central modeling innovation is the coupling of mRNA and nRNA through
the RNA lifecycle. At steady state, the kinetic model of transcription,
splicing, and degradation gives:

$$
\text{DNA} \xrightarrow{r_{\text{txn}}} \text{nRNA} \xrightarrow{r_{\text{splice}}}
\text{mRNA} \xrightarrow{r_{\text{deg}}} \varnothing
$$

The steady-state nascent fraction for transcript $t$ is:

$$
\beta_t = \frac{1}{1 + \rho_t}, \quad
\rho_t = \frac{r_{\text{splice},t}}{r_{\text{deg},t}}
$$

where $\rho_t$ is the ratio of splicing rate to degradation rate.

**Key properties:**
- When splicing is fast relative to degradation ($\rho_t \gg 1$):
  $\beta_t \to 0$, and most RNA is mature mRNA.
- When splicing is slow ($\rho_t \ll 1$): $\beta_t \to 1$, and most RNA
  is nascent pre-mRNA.
- $\beta_t$ is per-transcript, allowing different isoforms of the same
  gene to have different splicing kinetics.

The derived abundances are:

$$
\text{mRNA}_t = \theta_t (1 - \beta_t), \quad
\text{nRNA}_t = \theta_t \beta_t
$$

This parameterization halves the effective degrees of freedom compared to
treating mRNA and nRNA as fully independent, encodes the biological
constraint that observing mRNA from a transcript increases the probability
of observing nRNA from the same transcript, and prevents the EM from
freely trading mass between mRNA and nRNA of unrelated transcripts.

---

## 3. Pipeline Architecture

Rigel operates through a two-stage single-pass pipeline.

### 3.1 Stage 1: BAM Scan and Model Training

A high-performance C++ BAM scanner reads the name-sorted input BAM once
using htslib, performing the following operations per fragment:

**Fragment resolution.** Each fragment (mate pair or single read) is
resolved against a cgranges-based interval index to determine:
- Exonic overlap with each candidate transcript (in base pairs)
- Intronic overlap with each candidate transcript
- Intergenic extent (bases outside any annotated gene)
- Splice junction classification: annotated (matching known intron
  boundaries), unannotated (novel boundaries), or unspliced

**Read 2 strand flip.** For paired-end reads, the BAM scanner
normalizes strand representation by flipping R2's alignment strand,
ensuring all exon blocks carry R1's alignment strand regardless of
which mate contributed the alignment. This foundational normalization
enables protocol-agnostic downstream analysis.

**Fragment classification.** Each fragment is classified into one of
the following categories:
- `FRAG_UNIQUE`: Maps to a single transcript unambiguously
- `FRAG_AMBIG_SAME_STRAND`: Maps to multiple transcripts of the same
  gene — resolved by the EM
- `FRAG_MULTIMAPPER`: Multiple globally optimal alignments (NH > 1)
- `FRAG_CHIMERIC`: Inter-chromosome, inter-locus, or mixed-strand mappings

**Strand model training.** The strand model is trained from spliced
fragments with known splice-junction strand (the XS or ts BAM tag). A
Beta posterior on the R1-sense probability is computed:

$$
p_{\text{r1s}} = \frac{n_{\text{same}} + \kappa/2}{n_{\text{same}} + n_{\text{opp}} + \kappa}
$$

where $n_{\text{same}}$ counts fragments where the exon alignment strand
matches the splice-junction strand, $n_{\text{opp}}$ counts mismatches,
and $\kappa$ is a prior pseudo-count (default 2.0, giving a uniform
Beta(1,1) prior). The strand specificity is:

$$
\text{SS} = \max(p_{\text{r1s}}, 1 - p_{\text{r1s}}) \in [0.5, 1.0]
$$

The library protocol is determined as R1-sense ($p_{\text{r1s}} > 0.5$)
or R1-antisense ($p_{\text{r1s}} < 0.5$). Both protocols are handled
identically through the same code path because all downstream
computations use gene-relative sense/antisense classifications rather
than read-relative orientations.

**Fragment length model training.** Separate fragment length histograms
are trained for RNA and gDNA from uniquely mapped fragments:
- RNA model: trained from spliced fragments (which are definitively RNA)
- gDNA model: trained from intergenic fragments

Each model stores per-length bin probabilities with Laplace smoothing for
lengths $[0, L_{\max})$ (default $L_{\max} = 1000$) and an exponential
tail decay beyond $L_{\max}$:

$$
\log p(k) = \log p(L_{\max}) + (k - L_{\max}) \cdot \log(0.99), \quad k > L_{\max}
$$

**Buffering.** Fragments are stored in a memory-efficient columnar format
with configurable memory limits (default 2 GB). When the buffer fills,
chunks are spilled to disk in a temporary directory and recombined during
the quantification stage.

### 3.2 Stage 2: Locus-Level EM Quantification

**Pre-EM accumulation.** Before entering the EM, four accumulator arrays
are computed from single-gene fragments:
- Per-gene sense and antisense counts (for gDNA initialization)
- Per-transcript intronic sense and antisense counts (for nRNA
  initialization)

**Fragment routing.** Fragments are routed into compressed sparse row
(CSR) format. For each fragment $f$, the CSR stores:
- Candidate transcript indices
- Per-candidate log-likelihoods incorporating strand, insert length,
  alignment penalties, and effective length
- Coverage weights
- Splice type and strand classification

**Locus partitioning.** Transcripts are partitioned into independent loci
via union-find with path compression. An edge exists between two
transcripts if any fragment maps to both. Each connected component defines
a locus that is solved independently.

**EM optimization.** For each locus, the EM solves the 2T + 1 component
mixture (Section 4). Post-EM pruning (Section 4.6) removes low-evidence
components and re-optimizes.

**Output assembly.** Per-transcript mRNA and nRNA abundances are computed
from the converged $\theta$ and $\beta$ parameters. TPM values are
computed using mRNA effective lengths. Results are assembled into
per-transcript, per-gene, and per-locus tables.

---

## 4. EM Algorithm

### 4.1 Component Layout

For a locus with $T$ transcripts, the EM operates on a vector of
$2T + 1$ components:

| Index range | Component | Description |
|-------------|-----------|-------------|
| $[0, T)$ | mRNA$_t$ | Mature RNA from transcript $t$ |
| $[T, 2T)$ | nRNA$_t$ | Nascent RNA from transcript $t$ |
| $2T$ | gDNA | Genomic DNA background for the locus |

The abundance vector $\boldsymbol{\theta}$ lies on the simplex:
$\sum_{k=0}^{2T} \theta_k = 1$.

### 4.2 Fragment Likelihood

For each fragment $f$ and candidate component $c$, the log-likelihood
combines four terms:

$$
\log \ell(f, c) = \log p_{\text{strand}}(f, c)
  + \log p_{\text{insert}}(l_f, c)
  + \log p_{\text{align}}(f, c)
  - \log K_c(f)
$$

where:

**Strand likelihood.** The strand probability depends on whether the
fragment aligns sense or antisense to the gene:

$$
p_{\text{strand}}(f, c) = \begin{cases}
p_{\text{r1s}} & \text{if exon strand = gene strand (sense)} \\
1 - p_{\text{r1s}} & \text{if exon strand} \neq \text{gene strand (antisense)}
\end{cases}
$$

For gDNA, which is unstranded, the strand likelihood is 0.5 regardless
of alignment orientation.

**Insert size likelihood.** The probability of observing insert length
$l_f$ under the RNA or gDNA fragment length model:

$$
p_{\text{insert}}(l_f, c) = \begin{cases}
p_{\text{RNA}}(l_f) & \text{if } c \in \{\text{mRNA}, \text{nRNA}\} \\
p_{\text{gDNA}}(l_f) & \text{if } c = \text{gDNA}
\end{cases}
$$

**Alignment penalty.** A multiplicative log-penalty for alignment
imperfections:

$$
\log p_{\text{align}}(f, c) = n_{\text{overhang}} \cdot \log \alpha_o
  + n_{\text{mismatch}} \cdot \log \alpha_m
  + \mathbb{1}[\text{unannot. splice}] \cdot \log \alpha_s
$$

where $\alpha_o = 0.01$ (per-base overhang penalty), $\alpha_m = 0.1$
(per-mismatch penalty from the NM tag), and $\alpha_s = 0.01$ (gDNA
penalty for unannotated splice junctions). A fragment with 5 bp overhang
and 2 mismatches receives a penalty of $5 \times \log(0.01) + 2 \times
\log(0.1) \approx -27.6$ in log space.

**Effective length normalization.** The normalization constant $K_c(f)$
accounts for the probability of observing a fragment of length $l_f$ at
a random position within the transcript:

$$
K_c = \text{eff\_len}(t, c) = (L_t + 1) \cdot (\text{CDF}[k] - P(0))
  - \text{CMOM}[k]
$$

where $k = \min(L_t, L_{\max})$, CDF is the cumulative distribution
function of the fragment length model, and CMOM is the cumulative first
moment. This is the Salmon-style effective length correction (Patro et
al., 2017), computed efficiently via pre-tabulated cumulative arrays.

For mRNA components, $L_t$ is the sum of exon lengths (spliced
transcript length). For nRNA components, $L_t$ is the genomic span of
the transcript (including introns). For the gDNA component, $L_t$ is
the genomic span of the locus.

### 4.3 E-Step

The E-step computes posterior responsibilities for each fragment $f$
over its candidate components:

$$
r(f, c) = \frac{\theta_c \cdot \ell(f, c) / K_c}
  {\sum_{c'} \theta_{c'} \cdot \ell(f, c') / K_{c'}}
$$

Fragments that map unambiguously to a single component (unique mappers
to single-transcript loci) bypass the E-step entirely and contribute
directly to the M-step accumulators.

### 4.4 M-Step

The M-step updates the abundance vector using Maximum A Posteriori (MAP)
estimation with a Dirichlet prior:

$$
\theta_c^{(\text{new})} = \frac{1}{Z} \left(
  U_c + \sum_f r(f, c) + \alpha_c
\right)
$$

where $U_c$ is the unambiguous fragment count for component $c$, the
sum is over ambiguous fragments, and $\alpha_c$ is the Dirichlet prior
for component $c$. The normalization constant $Z$ ensures the simplex
constraint.

**Prior structure.** The prior $\alpha_c$ has two additive components:

$$
\alpha_c = \alpha_{\text{base}} + \gamma_{\text{ovr}} \cdot \text{OVR}_c
$$

where $\alpha_{\text{base}} = 0.01$ is a flat pseudo-count providing
numerical stability, and $\text{OVR}_c$ is the coverage-weighted One
Virtual Read share (Section 5.3). The OVR distributes a single virtual
read's worth of prior mass across components proportionally to their
coverage evidence:

$$
\text{OVR}_c = \frac{\sum_f w(f, c) \cdot r(f, c)}
  {\sum_{c'} \sum_f w(f, c') \cdot r(f, c')}
$$

where $w(f, c)$ is the positional coverage weight of fragment $f$ under
component $c$.

**Nascent fraction update.** The per-transcript nascent fraction $\beta_t$
is updated via MAP estimation with a Beta prior:

$$
\beta_t^{(\text{new})} = \frac{N_{t,\text{nRNA}} + \alpha_t^{(\beta)} - 1}
  {N_{t,\text{total}} + \alpha_t^{(\beta)} + \beta_t^{(\beta)} - 2}
$$

where $N_{t,\text{nRNA}} = \sum_f r(f, t, \text{nRNA})$,
$N_{t,\text{total}} = N_{t,\text{mRNA}} + N_{t,\text{nRNA}}$, and
$\text{Beta}(\alpha_t^{(\beta)}, \beta_t^{(\beta)})$ is the hierarchical
prior (Section 5.1). The result is clamped to $[\varepsilon, 1-\varepsilon]$
with $\varepsilon = 10^{-8}$.

### 4.5 SQUAREM Acceleration

The EM is accelerated using SQUAREM (Squared Iterative Methods; Varadhan
and Roland, 2008), which operates on the abundance vector $\boldsymbol{\theta}$
(the simplex of total component abundances including gDNA):

1. Compute one EM step: $\boldsymbol{\theta}^{(k)} \to \boldsymbol{\theta}^{(k+1)}$
2. Compute extrapolation direction:
   $\Delta = \boldsymbol{\theta}^{(k+1)} - \boldsymbol{\theta}^{(k)}$
3. Evaluate extrapolated point:
   $\boldsymbol{\theta}^{(k+2)} = \boldsymbol{\theta}^{(k+1)} + \lambda \Delta$
   with adaptive step size $\lambda$
4. Project back to the simplex if necessary
5. Revert to the standard EM step if the extrapolation increases the
   objective

The nascent fraction $\beta_t$ is **not** included in the SQUAREM state
vector. It is derived deterministically from the M-step ratio after each
EM iteration and clamped to $[\varepsilon, 1-\varepsilon]$. This avoids
boundary violations from extrapolation and requires no structural changes
to the SQUAREM infrastructure.

In practice, SQUAREM achieves convergence in approximately 10–50× fewer
iterations than standard EM on real RNA-seq data.

### 4.6 Post-EM Pruning

After initial convergence, Rigel identifies components that are likely
false positives—components sustained by prior mass rather than data
evidence. For each component $c$:

1. Compute the evidence ratio: $E_c = D_c / \alpha_c$, where $D_c$ is
   the data-contributed count and $\alpha_c$ is the prior.
2. If $U_c = 0$ (no unambiguous fragments) and $E_c < \tau_{\text{prune}}$
   (default $\tau_{\text{prune}} = 0.1$), zero the component.
3. Re-run the EM to redistribute the pruned mass.

This two-pass strategy allows the initial EM to explore the full
component space, then eliminates ghost components that cannot sustain
themselves without prior support.

### 4.7 Confidence-Based Assignment

For each ambiguous fragment, the RNA-normalized posterior is computed:

$$
p_{\text{RNA}}(f, t) = \frac{r(f, t, \text{mRNA}) + r(f, t, \text{nRNA})}
  {\sum_{t'} [r(f, t', \text{mRNA}) + r(f, t', \text{nRNA})]}
$$

If $\max_t p_{\text{RNA}}(f, t) \geq \tau_{\text{conf}}$ (default 0.95),
the fragment is classified as a high-confidence assignment.

---

## 5. Hierarchical Prior System

### 5.1 Nascent Fraction Prior

The nascent fraction $\beta_t$ for each transcript $t$ is given a
$\text{Beta}(\alpha_t, \beta_t)$ prior computed through a three-level
empirical Bayes hierarchy:

$$
\text{transcript} \leftarrow \text{TSS group} \leftarrow
\text{locus-strand} \leftarrow \text{global}
$$

At each level, the prior is parameterized by a mean estimate $\hat{\beta}$
and a concentration $\kappa = \alpha + \beta$ (effective sample size),
giving $\alpha = \hat{\beta} \kappa$ and $\beta = (1 - \hat{\beta}) \kappa$.

#### Level 0: Global prior

The global nascent fraction is a weakly informative prior of
$\hat{\beta}_{\text{global}} = 0.5$ with $\kappa = 2.0$, corresponding
to a uniform $\text{Beta}(1, 1)$ prior.

#### Level 1: Locus-strand

For each (locus, strand) combination, the nascent fraction is estimated
from the ratio of unspliced to total exonic fragments across all
transcripts on that strand:

$$
\hat{\beta}_{\ell,s} = \frac{\sum_{t \in (\ell,s)} U_t^{\text{unspliced}}}
  {\sum_{t \in (\ell,s)} (U_t^{\text{spliced}} + U_t^{\text{unspliced}})}
$$

If the total evidence exceeds a minimum threshold (default 10 fragments),
this estimate is used with $\kappa$ estimated via Method of Moments
(Section 5.4). Otherwise, the global prior ($\hat{\beta} = 0.5$,
$\kappa = 2$) is used.

The $\kappa$ pulling locus-strand estimates toward the global mean is
auto-estimated or manually set via `--nrna-frac-kappa-global`.

#### Level 2: TSS group

Transcripts sharing a transcription start site (TSS) are expected to have
similar nascent fractions because they share the same promoter and
5' kinetics. Rigel groups transcripts by fuzzy TSS matching using
**single-linkage clustering**: transcripts on the same chromosome and
strand whose 5' positions lie within a configurable window (default
200 bp) are merged into a single TSS group.

The 5' position is defined as the genomic start coordinate for
positive-strand transcripts and the genomic end coordinate for
negative-strand transcripts.

Within each TSS group, the nascent fraction is estimated from the
pooled spliced/unspliced ratio, with $\kappa$ shrinking toward the
locus-strand estimate.

#### Level 3: Transcript

If the transcript has sufficient unique evidence (≥ 10 fragments),
its own spliced/unspliced ratio provides the finest-grained estimate,
with $\kappa$ shrinking toward the TSS group.

#### Cascade logic

The cascade proceeds from the transcript level upward. At each level,
if the evidence exceeds the minimum threshold, the estimate at that
level is used; otherwise, the algorithm falls back to the next coarser
level:

$$
(\hat{\beta}_t, \kappa_t) = \begin{cases}
(\hat{\beta}_t^{\text{obs}}, n_t) & \text{if } n_t \geq n_{\min} \\
(\hat{\beta}_{\text{tss}}, \kappa_{\text{tss}}) & \text{if } n_{\text{tss}} \geq n_{\min} \\
(\hat{\beta}_{\ell,s}, \kappa_{\ell,s}) & \text{if } n_{\ell,s} \geq n_{\min} \\
(0.5, 2.0) & \text{otherwise}
\end{cases}
$$

### 5.2 gDNA Prior

The gDNA rate prior uses a two-level empirical Bayes hierarchy:

$$
\text{locus} \leftarrow \text{chromosome} \leftarrow \text{global}
$$

#### The antisense principle

At perfect strand specificity ($\text{SS} = 1.0$), every fragment aligning
antisense to a gene must originate from gDNA, because both mRNA and nRNA
are stranded (sense-only) while gDNA is unstranded:

$$
\text{gDNA}_{\text{total},g} = 2 \times \text{antisense}_g
$$

At imperfect strand specificity, some sense RNA fragments are flipped to
the antisense strand by library preparation noise. The corrected
antisense count accounts for this:

$$
\text{corrected}_g = \max\!\left(0, \;
  \text{antisense}_g - \text{sense}_g \times (1 - \text{SS})\right)
$$

$$
\text{gDNA}_{\text{init},g} = 2 \times \text{corrected}_g
$$

#### Density estimation

A complementary gDNA signal comes from intergenic fragments—reads mapping
outside all annotated gene bodies, which are unambiguously gDNA:

$$
\text{gDNA}_{\text{density}} = \frac{N_{\text{intergenic}}}
  {L_{\text{intergenic}}}
$$

where $N_{\text{intergenic}}$ is the count of intergenic fragments and
$L_{\text{intergenic}}$ is the total intergenic sequence length. The
expected gDNA count at a locus is then:

$$
\text{gDNA}_{\text{density},\ell} = \text{gDNA}_{\text{density}} \times
  L_{\text{exonic},\ell}
$$

#### Hybrid estimation

The strand-based and density-based estimates are combined via
inverse-variance weighting:

$$
W = (2 \cdot \text{SS} - 1)^2 \in [0, 1]
$$

$$
\hat{\gamma}_\ell = W \cdot \hat{\gamma}_{\text{strand},\ell}
  + (1 - W) \cdot \hat{\gamma}_{\text{density},\ell}
$$

The weight $W$ is the squared reliability of the strand signal:
- At $\text{SS} = 1.0$: $W = 1$, pure strand estimation
- At $\text{SS} = 0.5$: $W = 0$, pure density estimation
- Smooth, continuous transition with no hard cutoffs

#### Hierarchical shrinkage

The per-locus gDNA rate is shrunk toward chromosome-level and
global estimates through two levels of Beta-distribution shrinkage
with auto-estimated $\kappa$ parameters via Method of Moments.

### 5.3 One Virtual Read (OVR) Prior

The OVR prior addresses a fundamental asymmetry in the Dirichlet prior:
a flat pseudo-count $\alpha$ assigns equal prior mass to all components
regardless of whether any fragment evidence supports them. For loci with
many transcripts (e.g., alternative splicing hubs with 40+ isoforms), the
accumulated flat prior mass ($83 \times 0.5 = 41.5$ virtual fragments
for a 41-transcript locus) can overwhelm the actual data and create
self-sustaining false positives.

The OVR replaces the dominant prior source with coverage-weighted
geometric evidence:

$$
\alpha_c = \alpha_{\text{base}} + \gamma_{\text{ovr}} \cdot \text{OVR}_c
$$

where $\alpha_{\text{base}} = 0.01$ provides a numerical stability floor
and $\text{OVR}_c$ distributes a single virtual read across components
proportionally to their accumulated coverage weights. At $\alpha_{\text{base}}
= 0.01$, the total flat prior mass for an 83-component locus is only
$0.83$, while the OVR contributes $1.0$—making coverage evidence the
primary source of prior information.

Components with no geometric evidence receive $\alpha_c \approx 0.01$,
which is insufficient to sustain a false positive through the EM. Only
components with genuine fragment overlap receive meaningful prior support.

### 5.4 Method of Moments for κ Estimation

At each level of the hierarchy, the shrinkage concentration $\kappa$ is
estimated via Method of Moments (MoM) from the observed distribution of
nascent fractions (or gDNA rates) across features at that level.

Given observed proportions $\{p_1, \ldots, p_n\}$ from features with
sufficient evidence:

$$
\bar{p} = \frac{1}{n} \sum_i p_i, \quad
s^2 = \frac{1}{n-1} \sum_i (p_i - \bar{p})^2
$$

$$
\hat{\kappa} = \frac{\bar{p}(1 - \bar{p})}{s^2} - 1
$$

The estimated $\hat{\kappa}$ is clamped to $[\kappa_{\min}, \kappa_{\max}]$
(default $[2.0, 200.0]$). If fewer than $n_{\min}$ features pass the
evidence filter (default 20), a fallback $\kappa$ is used (default 5.0).

**Intuition:** High variance in observed proportions ($s^2$ large) yields
low $\kappa$ (weak shrinkage, allowing heterogeneity). Low variance yields
high $\kappa$ (strong shrinkage toward the group mean). This data-adaptive
behavior avoids the need for manual tuning in most scenarios.

---

## 6. Strand Model

### 6.1 Training

The strand model is trained from spliced fragments whose splice-junction
strand is known (via the XS or ts BAM tag). After R2 strand normalization
(Section 3.1), the BAM scanner accumulates a 2×2 contingency table:

| | SJ strand = exon strand | SJ strand ≠ exon strand |
|---|---|---|
| **Count** | $n_{\text{same}}$ | $n_{\text{opp}}$ |

The R1-sense probability is estimated as:

$$
p_{\text{r1s}} = \frac{n_{\text{same}} + \kappa/2}
  {n_{\text{same}} + n_{\text{opp}} + \kappa}
$$

with $\kappa = 2.0$ (uniform $\text{Beta}(1, 1)$ prior).

### 6.2 Protocol Detection

The library protocol is classified as:
- **R1-sense** (e.g., KAPA Stranded): $p_{\text{r1s}} > 0.5$
- **R1-antisense** (e.g., Illumina dUTP/TruSeq): $p_{\text{r1s}} < 0.5$

### 6.3 Protocol-Agnostic Design

A critical design property is that all downstream computations use
**gene-relative** sense/antisense classifications rather than
read-relative orientations. The classification function:

$$
\text{is\_antisense}(f) =
\begin{cases}
\text{true} & \text{if } p_{\text{strand}}(\text{exon}, \text{gene}) < 0.5 \\
\text{false} & \text{otherwise}
\end{cases}
$$

produces correct gene-relative labels for both R1-sense and R1-antisense
protocols through the same code path:

| Protocol | Exon = gene | $p$ returned | Antisense? |
|----------|------------|--------------|------------|
| R1-sense | Yes | 0.95 | No (correct) |
| R1-sense | No | 0.05 | Yes (correct) |
| R1-antisense | Yes | 0.05 | Yes (correct) |
| R1-antisense | No | 0.95 | No (correct) |

This eliminates the need for protocol-specific code branches and ensures
that all strand-dependent computations—gDNA estimation, nRNA estimation,
fragment scoring—operate identically regardless of the library preparation
protocol.

---

## 7. Fragment Length Model

### 7.1 Histogram Representation

Rigel maintains separate fragment length distributions for RNA and gDNA,
each represented as a smoothed histogram:

- **Bins:** Individual probability bins for lengths $[0, L_{\max})$
  (default $L_{\max} = 1000$ bp)
- **Overflow:** A single bin for lengths $\geq L_{\max}$
- **Smoothing:** Laplace smoothing ensures non-zero probability for all
  lengths
- **Tail decay:** For lengths beyond $L_{\max}$, an exponential tail
  provides graceful probability decay:

$$
\log p(k) = \log p(L_{\max}) + (k - L_{\max}) \cdot \log(0.99),
\quad k > L_{\max}
$$

### 7.2 Effective Length

The effective length of a transcript accounts for the fragment length
distribution using the Salmon-style eCDF correction:

$$
K_t = (L_t + 1) \cdot (\text{CDF}[k] - P(0)) - \text{CMOM}[k]
$$

where $k = \min(L_t, L_{\max})$, CDF is the cumulative distribution
function, and CMOM is the cumulative first moment of the fragment length
distribution. This correction is computed via pre-tabulated cumulative
arrays and vectorized across all transcripts.

The effective length differs by component type:
- **mRNA:** Spliced transcript length (sum of exon lengths)
- **nRNA:** Genomic span of the transcript (exons + introns)
- **gDNA:** Genomic span of the locus

---

## 8. Coverage Weight Model

### 8.1 Positional Coverage Trapezoid

Under uniform coverage, the expected fragment density varies across the
transcript due to boundary effects: fragments cannot start before the
5' end or extend past the 3' end. Rigel models this as a trapezoid:

- **Interior plateau:** Fragments in the middle of the transcript have
  uniform coverage (weight = 1.0)
- **5' ramp:** Linearly increasing coverage from transcript start over
  a distance of $\bar{l}/2$ (half the mean fragment length)
- **3' ramp:** Linearly decreasing coverage toward the transcript end

For each fragment $f$ mapping to transcript $t$, the coverage weight is
computed by analytically integrating the trapezoid function over the
fragment's genomic span:

$$
w(f, t) = \frac{\bar{l}/2}{\text{MeanCapacity}(f, t)}
$$

where MeanCapacity is the average of the trapezoid function across the
fragment interval. Interior fragments receive weight 1.0; edge fragments
receive weight > 1.0, reflecting the lower probability of observing a
fragment at transcript boundaries.

### 8.2 Short Transcript Handling

For transcripts shorter than twice the mean fragment length
($L_t < \bar{l}$), the two ramps meet in the middle and the trapezoid
never reaches its maximum height. The maximum capacity is:

$$
w_{\max} = \min\!\left(\frac{\bar{l}}{2}, \frac{L_t}{2}\right)
$$

This naturally handles the increased positional uncertainty for short
transcripts without special-case logic.

---

## 9. nRNA Initialization

### 9.1 Per-Transcript Intronic Evidence

Nascent RNA is the only stranded source that produces fragments
overlapping introns. The per-transcript nRNA initialization uses the
strand-corrected intronic sense excess:

$$
\text{nrna\_init}_t = \max\!\left(0, \;
  \frac{\text{intronic\_sense}_t - \text{intronic\_antisense}_t}
  {2 \cdot \text{SS} - 1}\right)
$$

This formula is valid when $2 \cdot \text{SS} - 1 > 0.2$ (a numerical
stability threshold). The subtraction removes the gDNA contribution
(symmetric across strands), isolating the stranded nRNA signal.

### 9.2 Single-Exon Transcripts

Transcripts with a single exon have no annotated introns. Since nRNA
from such transcripts is structurally identical to mRNA (no introns to
distinguish), the nRNA initialization is set to zero: $\text{nrna\_init}_t = 0$.
The EM may still assign some mass to the nRNA component if the likelihood
warrants it.

### 9.3 Ambiguous Fragment Handling

For fragments mapping ambiguously to multiple transcripts of the same
gene (`FRAG_AMBIG_SAME_STRAND`), intronic evidence is distributed
fractionally. If a fragment has $n$ candidate transcripts and overlaps
an intron of candidate $k$ (as determined by $\text{intron\_bp}[k] > 0$),
then candidate $k$ receives $1/n$ units of intronic evidence. Candidates
with no intronic overlap receive no count.

---

## 10. gDNA Initialization

### 10.1 The Antisense Signal

Antisense fragments provide the primary gDNA estimator. In a stranded
library, mRNA and nRNA produce exclusively sense fragments relative to
the gene strand. Any antisense fragment must therefore originate from
either gDNA (which is unstranded) or library preparation noise (strand
flip).

The strand-corrected gDNA initialization per gene is:

$$
\text{gDNA\_init}_g = 2 \times \max\!\left(0, \;
  N_g^{\text{anti}} - N_g^{\text{sense}} \times (1 - \text{SS})\right)
$$

The factor of 2 accounts for the sense projection: if $n$ fragments are
observed on the antisense strand from gDNA, an equal number are expected
on the sense strand.

### 10.2 Restriction to Unspliced Fragments

Only unspliced fragments contribute to the gene-level gDNA accumulators.
gDNA cannot produce splice junctions (neither annotated nor unannotated),
so spliced antisense fragments must originate from RNA of overlapping
genes on the opposite strand. Excluding them makes the gDNA
initialization more specific.

### 10.3 Behavior at Limiting Strand Specificity

| SS | Flip rate | Strand reliability | Behavior |
|----|-----------|-------------------|----------|
| 1.0 | 0.0 | Perfect | $\text{corrected} = \text{antisense}$ (exact) |
| 0.95 | 0.05 | High | Small correction; strand estimate dominates |
| 0.7 | 0.3 | Moderate | Strand/density hybrid with $W = 0.16$ |
| 0.5 | 0.5 | None | Strand estimate→0; pure density fallback |

At $\text{SS} = 0.5$ (unstranded libraries), strand-based gDNA estimation
is fundamentally impossible—this is a physics limitation. Rigel falls
back entirely to the density estimator through the smooth weighting
function $W = (2 \cdot \text{SS} - 1)^2$.

---

## 11. Locus Construction

### 11.1 Connected Components

Transcripts are partitioned into loci using union-find with path
compression. Two transcripts are connected if any fragment in the library
maps to both (sharing at least one fragment). This connectivity-based
definition captures biologically meaningful groupings:

- Overlapping transcripts of the same gene
- Convergent or divergent gene pairs sharing fragments in overlapping
  UTRs
- Distant transcripts linked by multimapping reads

Each connected component defines an independent sub-problem for the EM,
allowing the global quantification problem to be decomposed into many
small, tractable optimizations.

### 11.2 CSR Representation

Fragment-to-component mappings are stored in compressed sparse row (CSR)
format for memory efficiency. For each fragment $f$, the CSR stores:

| Field | Description |
|-------|-------------|
| `t_indices` | Candidate transcript indices |
| `log_liks` | Per-candidate log-likelihoods |
| `coverage_weights` | Positional coverage weights |
| `tx_starts`, `tx_ends` | Fragment coordinates on transcript |
| `count_cols` | Splice type × strand column encoding |

mRNA and nRNA candidates are stored during fragment routing. gDNA
candidates are appended per-locus during EM data construction: for each
unspliced fragment in a locus, an edge to the gDNA component ($2T$) is
added with the pre-computed gDNA likelihood.

---

## 12. Output Specification

### 12.1 Transcript-Level Output

The primary output (`quant.feather` / `quant.tsv`) contains one row per
annotated transcript with the following fields:

| Field | Description |
|-------|-------------|
| `transcript_id` | Transcript identifier from the GTF |
| `gene_id` | Parent gene identifier |
| `gene_name` | Gene symbol |
| `mrna` | Estimated mRNA fragment count |
| `nrna` | Estimated nRNA fragment count |
| `length` | Transcript length (bp) |
| `effective_length` | Fragment-length-corrected effective length |
| `tpm` | Transcripts Per Million (mRNA-based) |

TPM is computed from mRNA counts and effective lengths:

$$
\text{TPM}_t = \frac{\text{mRNA}_t / K_t}{\sum_{t'} \text{mRNA}_{t'} / K_{t'}} \times 10^6
$$

### 12.2 Gene-Level Output

Gene-level abundances are computed by summing transcript-level estimates
within each gene.

### 12.3 Locus-Level Output

Per-locus summaries include total fragment counts, mRNA/nRNA/gDNA
totals, number of transcripts and genes, and genomic coordinates.

### 12.4 Summary Statistics

The `summary.json` file reports:
- Library protocol and strand specificity
- Fragment length distribution statistics
- Alignment statistics (total, mapped, unique, multimapping, duplicate reads)
- Fragment classification counts (genic, intergenic, chimeric)
- Global quantification totals (mRNA/nRNA/gDNA fractions)

### 12.5 Annotated BAM

When requested, an annotated BAM is produced with per-fragment tags
encoding the EM assignment:

| Tag | Type | Description |
|-----|------|-------------|
| ZT | int | Assigned transcript index |
| ZG | int | Assigned gene index |
| ZP | float | Posterior probability |
| ZW | float | Coverage weight |
| ZC | int | Count column (splice type × strand) |
| ZH | string | Candidate transcript hash |
| ZN | int | Number of candidates |

---

## 13. Implementation

### 13.1 Software Architecture

Rigel is implemented as a Python package with a C++ native extension:

- **Python layer:** Pipeline orchestration, configuration management, CLI,
  I/O formatting, hierarchical prior computation, locus construction
- **C++ extension:** BAM scanning (htslib), fragment resolution, interval
  indexing (cgranges), EM kernel with SQUAREM acceleration, scoring

The C++ extension is compiled to a stable ABI wheel via scikit-build-core
and nanobind, supporting Python 3.12+.

### 13.2 Parallelism

The BAM scan stage uses multi-threaded htslib decompression. The EM stage
parallelizes across independent loci. Both stages are controlled by the
`--threads` parameter.

### 13.3 Memory Management

Fragment data is buffered in memory-efficient columnar format. When the
buffer exceeds the configured memory limit (default 2 GB), chunks are
spilled to a temporary directory and recombined during quantification.

### 13.4 Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| Python | ≥ 3.12 | Runtime |
| numpy | ≥ 1.26 | Numerical arrays |
| pandas | ≥ 2.1 | DataFrames |
| pyarrow | ≥ 14.0 | Feather I/O |
| pysam | ≥ 0.22 | htslib Python bindings |
| pyyaml | ≥ 6.0 | Configuration parsing |
| nanobind | ≥ 2.0 | C++/Python bindings (build) |
| scikit-build-core | ≥ 0.4.3 | Build system (build) |

---

## 14. Discussion

### 14.1 Relationship to Existing Methods

Rigel builds on the EM-based quantification framework pioneered by RSEM
(Li and Dewey, 2011) and refined by Salmon (Patro et al., 2017).
Key differences include:

1. **Three-species model.** Salmon and kallisto model only mRNA. Rigel
   adds nRNA and gDNA as first-class components, enabling libraries
   with significant nascent or genomic contamination to be analyzed
   without pre-filtering.

2. **Linked kinetic model.** The coupling of mRNA and nRNA through a
   shared abundance eliminates the doubled parameter space and
   incorporates biological knowledge about RNA processing.

3. **Alignment-based.** Unlike Salmon's quasi-mapping and kallisto's
   pseudoalignment, Rigel operates on full BAM alignments, providing
   access to alignment quality information (edit distance, overhang,
   splice junction annotation) that informs the likelihood model.

4. **Hierarchical Bayesian priors.** The multi-level shrinkage
   framework replaces the fixed Dirichlet prior with data-adaptive
   priors that borrow strength across transcripts, loci, and the genome.

### 14.2 Limitations and Future Directions

**Equivalence class collapse.** Transcripts with identical fragment
compatibility patterns (equivalence classes) are indistinguishable to
the EM, which splits counts evenly among them. Positional bias modeling
could partially resolve such degeneracies.

**Overlapping antisense genes.** When genes overlap on opposite strands,
RNA from one gene may appear as antisense evidence for the other,
inflating gDNA estimates. A per-locus gDNA architecture (rather than
per-gene) would address this.

**Mappability.** Rigel does not currently account for regional variation
in read mappability, which can bias both gDNA and transcript abundance
estimates. Integration of pre-computed mappability tracks is a planned
enhancement.

**Positional bias.** While the coverage weight model corrects for
boundary effects, more sophisticated positional bias models (5'/3' bias,
GC content) could improve accuracy, particularly for short transcripts
and competing isoforms.

---

## References

Bray NL, Pimentel H, Melsted P, Pachter L. Near-optimal probabilistic
RNA-seq quantification. *Nat Biotechnol.* 2016;34(5):525-527.

Li B, Dewey CN. RSEM: accurate transcript quantification from RNA-Seq
data with or without a reference genome. *BMC Bioinformatics.*
2011;12:323.

Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides
fast and bias-aware quantification of transcript expression. *Nat
Methods.* 2017;14(4):417-419.

Trapnell C, Williams BA, Pertea G, et al. Transcript assembly and
quantification by RNA-Seq reveals unannotated transcripts and isoform
switching during cell differentiation. *Nat Biotechnol.*
2010;28(5):511-515.

Varadhan R, Roland C. Simple and globally convergent methods for
accelerating the convergence of any EM algorithm. *Scand J Stat.*
2008;35(2):335-353.
