# Rigel: Bayesian RNA-seq Transcript Quantification with Joint mRNA, Nascent RNA, and Genomic DNA Deconvolution

**Matthew K. Iyer**

---

K_c = \text{eff\_len}(t, c) = (L_t + 1) \cdot (\text{CDF}[k] - P(0))

# Rigel Methods

This document describes the model that is currently implemented in the Rigel
codebase. It is intentionally implementation-facing: when the code and earlier
design notes disagree, this document follows the code.

For a theory-first description of the underlying inference problem, separate
from the current implementation, see `docs/THEORETICAL_MODEL.md`.

---

## 1. Scope and motivation

Rigel quantifies RNA-seq alignments while explicitly modeling three signal
sources in the same library:

- mature spliced RNA (`mRNA`)
- nascent RNA (`nRNA`)
- genomic DNA contamination (`gDNA`)

The main motivation is to avoid forcing intronic or contamination-derived signal
into transcript abundance estimates. In total RNA or imperfectly DNase-treated
libraries, that distinction matters.

---

## 2. Reference model

Rigel builds an index from a FASTA and GTF and materializes five core tables:

- reference lengths
- transcripts
- splice junctions
- genomic intervals
- nRNA spans

The nRNA table is central to the current implementation. An nRNA entry is not a
copy of a transcript. It is a unique genomic span:

$$
n = (\mathrm{ref}, \mathrm{strand}, \mathrm{start}, \mathrm{end})
$$

Each transcript $t$ maps to exactly one nRNA span through a transcript-to-nRNA
map $a(t)$. If multiple transcripts share the same genomic span, they share the
same nRNA entry.

For a locus with $T$ transcripts there are therefore two counts that matter:

- $T$: transcript count
- $N$: number of unique nRNA spans among those transcripts

Usually $N \leq T$, often strictly smaller in isoform-rich loci.

---

## 3. Pipeline architecture

Rigel runs in two logical stages.

### 3.1 Native BAM scan

A native htslib-backed scanner reads the BAM once and performs:

- fragment grouping by read name
- annotation overlap resolution against the interval index
- splice classification as `UNSPLICED`, `SPLICED_UNANNOT`, or `SPLICED_ANNOT`
- strand observation collection for the RNA strand model
- fragment-length observation collection for RNA and intergenic gDNA models
- buffering of resolved fragments into a columnar in-memory structure with optional spill-to-disk

The scanner normalizes paired-end strand representation by flipping read 2 so
that downstream logic consistently works in a read-1-relative orientation.

### 3.2 Fragment routing and loci

Buffered fragments are rescored and routed into CSR arrays containing:

- candidate component indices
- per-candidate log-likelihoods
- coverage weights
- transcript-space start and end coordinates
- per-fragment splice metadata

Transcripts are then partitioned into loci using connected components: two
transcripts belong to the same locus if at least one fragment can support both.

### 3.3 Batch locus EM

All loci are solved through a batched native EM call. Before the EM, Rigel
computes:

- per-transcript deterministic counts
- per-nRNA nRNA initialization from intronic strand-aware evidence
- per-locus empirical Bayes gDNA initializations
- per-nRNA Beta prior parameters for nascent fraction updates

---

## 4. Component layout

For a locus with $T$ transcripts and $N$ unique nRNA spans, Rigel solves a
mixture with:

$$
T + N + 1
$$

components:

| Index range | Component |
|-------------|-----------|
| $[0, T)$ | one mRNA component per transcript |
| $[T, T+N)$ | one nRNA component per unique nRNA span |
| $T+N$ | one merged gDNA component for the whole locus |

This is the main architectural difference from the older `2T + 1` design.
nRNA is no longer represented as one shadow per transcript.

Only unspliced ambiguous units receive a gDNA candidate. Spliced units compete
among RNA candidates only.

---

## 5. Fragment likelihood

For a fragment $f$ and candidate component $c$, Rigel combines four terms in
log space:

$$
\log \ell(f, c) =
\log p_{\mathrm{strand}}(f, c) +
\log p_{\mathrm{insert}}(f, c) +
\log p_{\mathrm{align}}(f, c) -
\log K(f, c)
$$

where:

- $p_{\mathrm{strand}}$ uses the trained RNA strand model for RNA candidates and
  a fixed value of $0.5$ for gDNA
- $p_{\mathrm{insert}}$ uses the RNA fragment-length model for RNA candidates and
  the intergenic model for gDNA
- $p_{\mathrm{align}}$ applies mismatch and overhang penalties, plus a gDNA
  penalty for unannotated splice junctions
- $K(f,c)$ is the effective-length term used for length normalization

The user-facing penalty parameters are:

- `overhang_alpha = 0.01`
- `mismatch_alpha = 0.1`
- `gdna_splice_penalty_unannot = 0.01`

Rigel stores these in config space as log penalties and applies them additively
in log space.

---

## 6. Strand model

Rigel retains three strand models, but only one is used for scoring.

### 6.1 Primary RNA strand model

The primary model is trained from annotated spliced fragments with unique gene
assignment. Those fragments are the cleanest RNA-only evidence in the pipeline.

Let $p_{\mathrm{r1s}}$ be the posterior mean probability that read 1 is sense to
the gene. With a Beta prior parameterized by `strand_prior_kappa`, the model
derives:

$$
\mathrm{SS} = \max(p_{\mathrm{r1s}}, 1 - p_{\mathrm{r1s}})
$$

and reports protocol as:

- `R1-sense` if $p_{\mathrm{r1s}} \ge 0.5$
- `R1-antisense` otherwise

### 6.2 Diagnostic strand models

Rigel also tracks:

- an exonic diagnostic model trained on all exonic signal
- an intergenic diagnostic model trained on intergenic fragments

These are reported in `summary.json` but are not used to score gDNA. gDNA is
always treated as unstranded for scoring.

---

## 7. Fragment-length model and effective length

Rigel trains separate fragment-length models for:

- RNA-like fragments
- intergenic gDNA-like fragments

Transcript effective lengths are computed from the RNA model. When no fragment-
length observations are available, Rigel falls back to:

$$
\max(L - \bar{f} + 1, 1)
$$

using a default mean fragment length.

At EM time, the native solver also applies per-fragment effective-length logic
through component-specific bias profiles:

- mRNA uses exonic transcript length
- nRNA uses genomic span length of the shared nRNA entity
- gDNA uses locus span

---

## 8. EM objective

Rigel performs MAP-EM on the component simplex with effective-length
normalization and a coverage-weighted One Virtual Read prior.

### 8.1 E-step

For each ambiguous unit, posterior mass is distributed over candidate
components in proportion to current weight times candidate likelihood.

### 8.2 M-step

The component simplex is updated from:

- deterministic counts
- expected ambiguous counts
- flat Dirichlet prior mass
- OVR prior mass

When a locus contains shared nRNA spans, the native M-step switches into a
hierarchical nRNA-group mode and enforces one nascent-fraction update per nRNA
span rather than one per transcript.

### 8.3 SQUAREM acceleration

The native solver accelerates the simplex updates with SQUAREM. This is the
dominant compute kernel in realistic runs and is where most runtime is spent.

### 8.4 Pruning

After convergence, components with zero deterministic support and insufficient
data-to-prior evidence may be pruned according to `prune_threshold`, then the
EM is rerun once to redistribute mass.

### 8.5 Confidence assignments

Rigel separately tracks high-confidence RNA assignments using the
RNA-normalized posterior and `confidence_threshold`.

---

## 9. nRNA prior system

The current nRNA prior is defined on shared nRNA spans, not transcripts.

### 9.1 Evidence model

Rigel accumulates strand-aware exon and intron evidence, then computes a hybrid
nRNA fraction estimate from two views of the data:

- density subtraction against the global gDNA background
- strand subtraction using the library strand specificity

The strand contribution is weighted by:

$$
W = (2 \cdot \mathrm{SS} - 1)^2
$$

so that density-based estimates dominate when the library is weakly stranded.

### 9.2 Hierarchy

The implemented hierarchy is:

$$
\mathrm{global} \rightarrow \mathrm{locus\text{-}strand} \rightarrow \mathrm{nRNA}
$$

The user-visible parameters are:

- `nrna_frac_kappa_global`
- `nrna_frac_kappa_locus`
- `nrna_frac_kappa_nrna`

The first two can be auto-estimated by Method of Moments subject to minimum-
evidence thresholds and lower/upper clamps. The last is a fixed pseudo-count for
the final Beta prior passed into the EM solver.

### 9.3 Initialization and gating

Rigel computes a separate `nrna_init` signal from intronic sense-excess counts.
nRNA components are disabled when they fail biological or structural checks,
including cases where the shared nRNA span is supported only by single-exon
transcripts.

---

## 10. gDNA prior system

gDNA is modeled as one merged locus-level component.

### 10.1 Background density

Rigel estimates a global intergenic density from total intergenic fragments and
the inferred intergenic territory of the reference.

### 10.2 Hierarchy

The gDNA prior hierarchy is:

$$
\mathrm{global} \rightarrow \mathrm{reference} \rightarrow \mathrm{locus}
$$

with optional auto-estimation of:

- `gdna_kappa_ref`
- `gdna_kappa_locus`

using Method of Moments and evidence thresholds.

### 10.3 gDNA strand penalty

gDNA contamination is strand-symmetric, so a single gDNA component is used per
locus. A constant log(0.5) term is added to the per-fragment gDNA likelihood,
reflecting the a priori expectation that gDNA fragments are equally likely to
originate from either strand. This effectively halves the gDNA likelihood
relative to strand-specific components, preventing gDNA from absorbing
sense-strand fragments that belong to nRNA.

---

## 11. Reporting semantics

Rigel reports several views of the fitted model.

### 11.1 Transcript table

`quant.feather` and `quant.tsv` contain transcript-level mRNA and nRNA values.
The `nrna` column is not a direct per-transcript latent state from the EM. It
is produced by equal-share fan-out of shared nRNA-span counts back to the
transcripts mapping to that span.

### 11.2 Gene table

`gene_quant.feather` and `gene_quant.tsv` sum transcript-level quantities by
gene after transcript fan-out.

### 11.3 Locus table

`loci.feather` and `loci.tsv` report the direct locus-level EM totals:

- `mrna`
- `nrna`
- `gdna`
- `gdna_init`
- `gdna_rate`

### 11.4 Detail table

`quant_detail.feather` and `quant_detail.tsv` are long-format QC tables over:

- transcript
- gene
- splice category
- source (`unambig` or `em`)

### 11.5 Annotated BAM

When enabled, a second BAM pass stamps each record with tags describing the
winning assignment, its pool, posterior, fragment class, number of candidates,
and splice type.

---

## 12. Implementation notes

Rigel is split between Python orchestration and native kernels.

### Python responsibilities

- CLI parsing and YAML merge logic
- index build and load
- pipeline orchestration
- prior estimation
- output writing

### Native responsibilities

- BAM scanning through htslib
- high-volume fragment resolution and routing
- batched locus EM
- SQUAREM acceleration
- annotated BAM writing

The project is packaged with scikit-build-core, nanobind, and CMake, and builds
stable-ABI wheels for Python 3.12+.
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
**single-linkage clustering**: transcripts on the same reference and
strand whose 5' positions lie within a configurable window (default
200 bp) are merged into a single TSS group.

The 5' position is defined as the genomic start coordinate for
pos-strand transcripts and the genomic end coordinate for
neg-strand transcripts.

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
	ext{locus} \leftarrow \text{reference} \leftarrow \text{global}
$$

#### The antisense principle

At perfect strand specificity ($\text{SS} = 1.0$), every fragment aligning
antisense to annotated transcripts must originate from gDNA, because both
mRNA and nRNA are stranded (sense-only) while gDNA is unstranded:

$$
\text{gDNA}_{\text{total},\ell} = 2 \times \text{antisense}_\ell
$$

At imperfect strand specificity, some sense RNA fragments are flipped to
the antisense strand by library preparation noise. The corrected
antisense count accounts for this:

$$
\text{corrected}_\ell = \max\!\left(0, \;
  \text{antisense}_\ell - \text{sense}_\ell \times (1 - \text{SS})\right)
$$

$$
\text{gDNA}_{\text{init},\ell} = 2 \times \text{corrected}_\ell
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

The per-locus gDNA rate is shrunk toward reference-level and
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
the transcript strand. Any antisense fragment must therefore originate
from either gDNA (which is unstranded) or library preparation noise
(strand flip).

The strand-corrected gDNA initialization per locus is:

$$
\text{gDNA\_init}_\ell = 2 \times \max\!\left(0, \;
  N_\ell^{\text{anti}} - N_\ell^{\text{sense}} \times (1 - \text{SS})\right)
$$

where $N_\ell^{\text{anti}}$ and $N_\ell^{\text{sense}}$ are the
unspliced antisense and sense counts summed over all transcripts in
locus $\ell$.

The factor of 2 accounts for the sense projection: if $n$ fragments are
observed on the antisense strand from gDNA, an equal number are expected
on the sense strand.

### 10.2 Restriction to Unspliced Fragments

Only unspliced fragments contribute to the locus-level gDNA accumulators.
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
inflating gDNA estimates. Rigel mitigates this by modeling gDNA at the
locus level (a single shadow per locus) and restricting gDNA accumulators
to unspliced fragments, but residual cross-talk remains possible at
highly overlapping antisense loci.

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
