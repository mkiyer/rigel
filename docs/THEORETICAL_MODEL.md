# Rigel Theoretical Model

This document is intentionally theory-first. It is not a description of the
current code. Its purpose is to define the mathematical problem Rigel is trying
to solve, starting from biology and probability rather than from existing
implementation choices.

The main use of this note is to provide a clean target for redesigning the EM.

For a compact locus-level MAP-EM formulation using the preferred `T + N + 2`
gDNA parameterization, see `docs/THEORETICAL_EM_OBJECTIVE.md`.

## 1. The Problem

We observe aligned sequencing fragments from a library that may contain a mixture
of three biological source classes:

- mature RNA
- nascent RNA
- genomic DNA contamination

The inferential task is to explain each observed fragment as arising from one of
those sources, then estimate:

- transcript-level mature RNA abundance
- transcript-group or span-level nascent RNA abundance
- locus-level genomic DNA abundance

The hard part is that these sources overlap in genomic coordinates and can be
partially confounded by alignment ambiguity, intronic signal, fragment length,
and imperfect library strandedness.

## 2. Observed Data and Latent Variables

For each fragment $f$, we observe a bundle of alignment-derived features:

$$
x_f = (s_f, j_f, a_f, l_f)
$$

where:

- $s_f$ is the observed strand-relative orientation
- $j_f$ is splice class
- $a_f$ is alignment geometry and compatibility with candidate components
- $l_f$ is fragment length or insert information

At a locus $\ell$, the latent source assignment for fragment $f$ is a discrete
variable:

$$
z_f \in \mathcal{C}_\ell
$$

where $\mathcal{C}_\ell$ is the set of components active in that locus.

## 3. Minimal Component Set

The biologically faithful component set for a locus is:

- one mature RNA component per transcript $t$
- one nascent RNA component per unique nascent span $n$
- genomic DNA components for the two DNA strands

So the cleanest first-principles component set is:

$$
\mathcal{C}_\ell = \{m_t\}_{t=1}^{T_\ell}
\cup \{n_k\}_{k=1}^{N_\ell}
\cup \{d_{\ell,+}, d_{\ell,-}\}
$$

This gives $T_\ell + N_\ell + 2$ components.

This is now the preferred redesign target. Rigel should keep one collapsed
public gDNA abundance output per locus, but the internal EM parameterization
should use two strand-specific gDNA components and sum them at reporting time.

## 4. Likelihood Factorization

The natural decomposition is:

$$
p(x_f \mid z_f = c)
= p(s_f \mid c) \; p(j_f \mid c) \; p(l_f \mid c) \; p(a_f \mid c)
$$

The details can vary, but this factorization separates conceptually distinct
signals:

- strand behavior
- splice behavior
- fragment-length behavior
- positional and alignment compatibility

This separation matters because the source of confusion in the current system is
not only parameter count. It is also that several conceptually different ideas
have been mixed together: likelihood design, initialization, priors, and pruning.

## 5. RNA Strand Theory

For mature RNA and nascent RNA, strand is informative because the underlying RNA
molecule is single-stranded and tied to the annotated transcription direction.

Under a stranded RNA-seq protocol, opposite-strand observations arise because of
library-preparation failure or protocol noise, not because the RNA exists on
both strands.

So for RNA components, the strand model should encode:

- high probability for the expected strand
- low but non-zero probability for the opposite strand when strandedness is
  imperfect

This is the right place for a strand-specificity parameter.

## 6. gDNA Strand Theory

For genomic DNA, the biology is different.

- DNA exists on both strands.
- A single contaminating DNA fragment comes from one strand or the other.
- Across many fragments, absent bias, the two strands should be represented
  symmetrically.

This means there are two different questions:

1. What is the per-fragment strand likelihood under gDNA?
2. What is the prior belief about the aggregate strand balance of gDNA at a
   locus?

These should not be conflated.

## 7. Two-Strand gDNA Model

The clean explicit gDNA model introduces two latent components per locus:

$$
d_{\ell,+}, \; d_{\ell,-}
$$

with total gDNA mass:

$$
G_\ell = \theta_{\ell,+}^{(d)} + \theta_{\ell,-}^{(d)}
$$

and strand fraction:

$$
\phi_\ell = \frac{\theta_{\ell,+}^{(d)}}{G_\ell}
$$

Then a strand-observed fragment has deterministic gDNA strand compatibility:

$$
p(s_f = + \mid d_{\ell,+}) = 1, \qquad
p(s_f = - \mid d_{\ell,+}) = 0
$$

$$
p(s_f = - \mid d_{\ell,-}) = 1, \qquad
p(s_f = + \mid d_{\ell,-}) = 0
$$

or, if numerical robustness is needed, those zeroes can be replaced by a tiny
error probability.

The important point is that under this model, a given DNA fragment competes only
with the matching DNA strand, not with both DNA strands equally.

## 8. Collapsed One-gDNA Model

If we collapse the two DNA strands into a single gDNA state $d_\ell$ and do not
retain strand-specific gDNA abundances, then the per-fragment gDNA strand
likelihood becomes the marginal over the hidden DNA strand:

$$
p(s_f \mid d_\ell)
= p(s_f \mid d_{\ell,+}) p(d_{\ell,+} \mid d_\ell)
+ p(s_f \mid d_{\ell,-}) p(d_{\ell,-} \mid d_\ell)
$$

If the two DNA strands are a priori equally likely, then:

$$
p(d_{\ell,+} \mid d_\ell) = p(d_{\ell,-} \mid d_\ell) = \frac{1}{2}
$$

and therefore:

$$
p(s_f \mid d_\ell) = \frac{1}{2}
$$

for either observed strand.

This is the key theoretical conclusion:

### A one-gDNA-state model with strand likelihood $0.5$ is the correct collapsed form of a symmetric two-strand gDNA model.

So changing the current one-gDNA strand term from $\log(0.5)$ to $\log(1.0)$
would not be a harmless reinterpretation. It would define a different and more
gDNA-favorable model.

## 9. Where Symmetry Belongs

The desire for a "gravitational pull" toward equal gDNA signal on the two
strands is mathematically sensible, but that pull belongs in a prior on the
strand split, not in the per-fragment likelihood.

The right abstraction is:

- total gDNA mass at the locus
- strand fraction of that gDNA mass

That is:

$$
G_\ell \ge 0, \qquad \phi_\ell \in [0,1]
$$

with:

$$
\theta_{\ell,+}^{(d)} = G_\ell \phi_\ell
$$

$$
\theta_{\ell,-}^{(d)} = G_\ell (1 - \phi_\ell)
$$

Then symmetry is encoded with a prior centered at $0.5$:

$$
\phi_\ell \sim \operatorname{Beta}(a_{\mathrm{sym}}, a_{\mathrm{sym}})
$$

or equivalently:

$$
\phi_\ell \sim \operatorname{Beta}\left(\frac{\kappa_{\mathrm{sym}}}{2},
\frac{\kappa_{\mathrm{sym}}}{2}\right)
$$

This prior does exactly what the intuition wants:

- mean at $0.5$
- weak spring when $\kappa_{\mathrm{sym}}$ is small
- strong spring when $\kappa_{\mathrm{sym}}$ is large

Crucially, it acts on the locus-level aggregate strand balance, not on each
fragment independently.

It is important to distinguish expected symmetry from realized counts.

If we fix

$$
\phi_\ell = \frac{1}{2}
$$

that does not force a locus with 10 gDNA fragments to appear as exactly 5 plus
and 5 minus. It only says that, conditional on the total number of gDNA
fragments $N_\ell$, the observed positive-strand count is binomial:

$$
K_\ell^+ \mid N_\ell, \phi_\ell = \frac{1}{2}
\sim \operatorname{Binomial}\left(N_\ell, \frac{1}{2}\right)
$$

So imbalances such as $1/9$, $2/8$, or $8/2$ are already allowed by a strict
50/50 expected symmetry model through ordinary sampling variation.

If real loci vary more than the binomial model predicts, then the right
extension is not to abandon symmetry but to let $\phi_\ell$ vary across loci
under the symmetric Beta prior above. The resulting marginal count model is then
Beta-Binomial rather than Binomial.

## 10. Theoretical Recommendation on gDNA

From first principles, the cleanest model is:

1. two explicit gDNA strand components per locus
2. deterministic or near-deterministic strand likelihood for each fragment
3. a symmetric prior on the gDNA strand fraction

This is the most biologically faithful and most interpretable formulation.

It is also now the preferred implementation target. The main reason is
architectural: the `T + N + 2` representation keeps the EM closer to an
ordinary mixture over simplex weights, avoids a separate collapsed-gDNA
nuisance-parameter update, and makes the strand-symmetry prior a coupling prior
on an ordinary pair of components.

However, there is an important caveat:

### The benefit of two gDNA strands is not per-fragment competition. The benefit is a cleaner parameterization of locus-level strand symmetry.

Per fragment, once the observed strand is known, the incompatible gDNA strand is
effectively absent. The competition is not between $d_{\ell,+}$ and $d_{\ell,-}$
for the same fragment. The competition is between RNA and the compatible DNA
strand, while the two DNA strands are coupled only through the prior on
$\phi_\ell$.

That means the case for two gDNA components depends on whether we want an
explicit strand-balance latent variable. For the redesign, the answer is yes:
the preferred internal representation is two gDNA components, while the
preferred public representation remains collapsed total gDNA abundance.

A collapsed single-component gDNA model with a $0.5$ strand term remains
mathematically coherent and remains the correct marginal story. It is simply no
longer the preferred implementation target.

## 11. Why a Single gDNA State with 1.0 Strand Likelihood Is Not Right

It is tempting to say:

- DNA exists on both strands
- therefore any observed strand is fully compatible with DNA
- therefore the strand likelihood should be $1.0$

That reasoning skips an important distinction between compatibility and
probability.

Yes, either observed strand is compatible with DNA. But if the model does not
explicitly represent which DNA strand generated the fragment, then the observed
strand remains a latent mixture outcome. Under symmetric DNA, each observed
strand has probability $0.5$, not $1.0$.

So:

- `compatible with either strand` does not imply
- `probability 1 under a collapsed non-strand-specific DNA state`

This is exactly the same reason that marginalizing over two equally likely
hidden states gives a factor of $1/2$.

## 12. A Minimal First-Principles Mixture Model

If we redesign from scratch, the model should ideally separate four ideas.

### 12.1 Structural likelihood

How well does a fragment fit a component in terms of:

- splice structure
- genomic or transcript compatibility
- fragment length
- strand behavior

This should be a likelihood statement, not a prior.

### 12.2 Source abundances

At each locus, there is a simplex or nonnegative mass vector over active
components.

This is the primary object the EM is estimating.

### 12.3 Biological coupling priors

These are principled priors on interpretable latent quantities, for example:

- nascent fraction within a shared nRNA group
- total gDNA contamination mass
- gDNA strand symmetry fraction

These are the correct place for biological inductive bias.

### 12.4 Initialization

Initialization should be treated as numerical help only. It should not carry
biological meaning that is not also represented in the actual probabilistic
model.

This distinction is important because one of the current system's weaknesses is
that some quantities are halfway between initialization, gating, and prior.

## 13. Suggested Theoretical Parameterization

The following parameterization is conceptually clean.

### Mature RNA

For each transcript $t$:

$$
M_t \ge 0
$$

### Nascent RNA

For each unique nRNA span $n$:

$$
N_n \ge 0
$$

or, if one wants a more structured model, use a total RNA mass plus a nascent
fraction inside each shared group.

### gDNA

For each locus $\ell$:

$$
G_\ell \ge 0, \qquad \phi_\ell \in [0,1]
$$

with:

$$
D_{\ell,+} = G_\ell \phi_\ell, \qquad
D_{\ell,-} = G_\ell (1 - \phi_\ell)
$$

This separates:

- how much gDNA exists at the locus
- how that gDNA is split across strands

That is cleaner than letting one monolithic gDNA state absorb both jobs.

## 14. Splicing and gDNA

Theoretical biology suggests:

- mature RNA can generate annotated spliced fragments and many unspliced ones,
  depending on the alignment and transcript context
- nascent RNA can generate unspliced fragments and some exon-overlapping signal
- gDNA should overwhelmingly generate unspliced genomic fragments

So the clean theoretical stance is that spliced gDNA likelihood should be zero
or extremely small. That is a structural statement, not a prior statement.

## 15. What Should Be Kept Minimal

If the redesign is guided by first principles, the model should aim to keep only
three kinds of prior information:

1. a prior on nascent fraction
2. a prior on total gDNA abundance
3. a prior on gDNA strand symmetry

Everything else should ideally live in the likelihood or in a purely numerical
initializer.

## 16. Recommended Working Conclusion

The best current theoretical position is:

1. The biology supports modeling gDNA as two strand-specific latent components
   per locus.
2. A single gDNA component with strand likelihood $0.5$ is the collapsed form of
   that two-strand model under a symmetric strand prior.
3. A single gDNA component with strand likelihood $1.0$ is not the same model
   and would artificially favor gDNA.
4. The "spring toward symmetry" should be represented as a prior on the gDNA
   strand fraction, not as a change from $0.5$ to $1.0$ in the collapsed
   likelihood.
5. If interpretability, implementation clarity, and solver simplicity are the
   priority, the two-strand gDNA model is the preferred redesign target.
6. Rigel should report total gDNA abundance by collapsing the inferred
   strand-specific gDNA pair at output time.

## 17. Immediate Design Implication

The core theoretical decision is now resolved:

### gDNA should be represented internally as an explicit strand-resolved latent variable pair, and reported externally as collapsed total gDNA abundance.

The clean internal representation is therefore:

$$
(G_\ell, \phi_\ell)
$$

or equivalently:

$$
(D_{\ell,+}, D_{\ell,-})
$$

with the public gDNA output defined by:

$$
G_\ell = D_{\ell,+} + D_{\ell,-}
$$

The collapsed one-component formulation remains useful as a marginal reference
model and as a conceptual cross-check, but it is no longer the preferred design
target for implementation.

## 18. How to Treat $\kappa_{\mathrm{sym}}$

The most principled next step is to treat $\kappa_{\mathrm{sym}}$ as a global
hyperparameter of the gDNA symmetry process and estimate it empirically from the
data.

The intended hierarchy is:

$$
\phi_\ell \mid \kappa_{\mathrm{sym}}
\sim \operatorname{Beta}\left(\frac{\kappa_{\mathrm{sym}}}{2},
\frac{\kappa_{\mathrm{sym}}}{2}\right)
$$

for loci or regions whose observed fragments are believed to be predominantly
gDNA.

This gives a clear interpretation:

- each locus has its own latent strand fraction $\phi_\ell$
- all loci borrow strength through one shared concentration parameter
   $\kappa_{\mathrm{sym}}$
- the global model enforces symmetry in expectation while letting individual
   loci vary substantially

In this view, $\kappa_{\mathrm{sym}}$ is not a hand-tuned penalty coefficient.
It is a data-derived description of how variable gDNA strand balance really is
across loci in that library.

## 19. What Data Should Inform $\kappa_{\mathrm{sym}}$

In theory, $\kappa_{\mathrm{sym}}$ should be estimated from regions where the
observed counts can be attributed to gDNA with high confidence.

The key point is that these regions do not need to be intergenic. They only need
to be RNA-poor.

That means the clean theoretical calibration set is:

- regions with negligible expected mature RNA abundance
- negligible expected nascent RNA abundance
- enough total counts to inform strand imbalance

Such regions may be:

- intergenic intervals
- intronic intervals of unexpressed genes
- exonic intervals of transcripts or genes with essentially zero RNA evidence

This matters because in capture-based or targeted libraries, intergenic space may
be too sparse and too unrepresentative to support stable estimation of gDNA
properties.

So the correct theory is not “estimate from intergenic regions”. The correct
theory is “estimate from regions that are effectively pure gDNA in that sample”.

## 20. A Theory-Level Calibration Model

Suppose we assemble a calibration set of regions

$$
\mathcal{R} = \{r_1, \dots, r_m\}
$$

that are believed to be gDNA-dominant.

For each such region $r$, let:

- $N_r$ be the total number of fragments attributed to gDNA in that region
- $K_r^+$ be the positive-strand count among them

Then the theory-level calibration model is:

$$
\phi_r \mid \kappa_{\mathrm{sym}}
\sim \operatorname{Beta}\left(\frac{\kappa_{\mathrm{sym}}}{2},
\frac{\kappa_{\mathrm{sym}}}{2}\right)
$$

$$
K_r^+ \mid N_r, \phi_r
\sim \operatorname{Binomial}(N_r, \phi_r)
$$

and therefore marginally:

$$
K_r^+ \mid N_r, \kappa_{\mathrm{sym}}
\sim \operatorname{BetaBinomial}\left(
N_r, \frac{\kappa_{\mathrm{sym}}}{2}, \frac{\kappa_{\mathrm{sym}}}{2}
\right)
$$

The role of the calibration step is then simply to estimate the shared
concentration parameter that best explains the observed cross-region dispersion.

## 21. Relationship to Initialization

This is best thought of as an initialization-level goal, not because it is
heuristic, but because it provides a library-level hyperparameter needed before
the main locus EM is fully specified.

The same theoretical issue appears for gDNA fragment-length estimation:

- RNA-rich regions do not give clean gDNA observations
- gDNA-dominant regions do

So the general initialization problem is:

1. identify regions with minimal RNA contamination
2. estimate gDNA-specific nuisance parameters from those regions
3. feed those estimated hyperparameters into the main abundance model

Under that view, $\kappa_{\mathrm{sym}}$ and the gDNA fragment-length
distribution are parallel statistical objects: both should be learned from
gDNA-dominant data, not from arbitrary genomic background.

## 22. What "RNA-Poor" Means

The correct unit of observation is the fragment, not the raw BAM record.

After alignment parsing and read-pair assembly, each sequenced molecule is
represented by one fragment object. A genomic region should therefore be judged
RNA-rich or RNA-poor by the fragment evidence it contains.

At the theory level, a region is RNA-poor if the posterior evidence for RNA is
small under the three strongest RNA-defining signals:

1. splicing evidence
2. exon-versus-intron-or-intergenic density contrast
3. strand asymmetry

These signals have different reliability and should be understood separately.

### 22.1 Splicing evidence

Splicing is the strongest and cleanest indicator of RNA.

A fragment that traverses an intron gap cannot arise from genomic DNA. Therefore,
spliced fragments are effectively pure RNA evidence.

In theory, for a region $r$, define:

$$
S_r = \text{number of spliced fragments overlapping } r
$$

Then large $S_r$ is direct evidence that the region is not gDNA-pure.

This is the most trusted RNA signal.

### 22.2 Density contrast evidence

RNA tends to create non-uniform coverage across annotation classes:

- mature RNA enriches exonic signal relative to intronic and intergenic signal
- nascent RNA increases intronic signal within transcribed spans
- gDNA contributes more uniformly with respect to exon versus intron identity,
  up to capture or assay-specific enrichment effects

So a region is less plausibly gDNA-pure when its observed fragment density is
better explained by an RNA-driven contrast than by background DNA.

At theory level, this can be represented through summary quantities such as:

$$
D_r^{\mathrm{exon}}, \qquad D_r^{\mathrm{intron}}, \qquad D_r^{\mathrm{intergenic}}
$$

or any equivalent density statistics derived from fragment counts per effective
base.

This signal is weaker than splicing because density differences can be affected
by enrichment, baiting, mapping, and annotation geometry.

### 22.3 Strand-asymmetry evidence

In a stranded library:

- RNA produces strand-asymmetric signal
- gDNA is symmetric in expectation

So a region with strong strand asymmetry is less plausibly gDNA-pure.

At theory level, for region $r$, let:

$$
N_r^+ = \text{positive-strand fragment count},
\qquad
N_r^- = \text{negative-strand fragment count}
$$

and total count:

$$
N_r = N_r^+ + N_r^-
$$

Then deviation of $N_r^+ / N_r$ from $1/2$ is evidence for RNA contamination,
with the strength of that evidence depending on:

- total count $N_r$
- library strandedness
- the global gDNA strand-variability model

This signal is especially useful when splicing is sparse but strandedness is
strong.

## 23. A Latent Purity View of Calibration Regions

It is useful to define a latent region-specific RNA contamination fraction:

$$
\omega_r \in [0,1]
$$

where:

- $\omega_r = 0$ means region $r$ is effectively pure gDNA
- $\omega_r > 0$ means region $r$ contains some RNA contribution

Then “RNA-poor” means:

$$
\omega_r \approx 0
$$

not necessarily that the region has no annotation or no expression label.

This is the clean way to state the calibration problem:

- we seek regions whose fragments are overwhelmingly explained by gDNA
- those regions may still overlap annotated exons or introns
- what matters is source purity, not annotation category

In practice, the three evidence channels above are exactly the signals that
should inform whether $\omega_r$ is close to zero.

The preferred theoretical formulation is soft rather than hard. That is, we do
not declare each region to be either usable or unusable for gDNA calibration.
Instead, each region contributes according to how close it appears to the
gDNA-dominant end of the spectrum.

Define a region-specific gDNA purity weight:

$$
w_r \in [0,1]
$$

with the interpretation:

- $w_r \approx 1$: region is close to gDNA-only
- $w_r \approx 0$: region contains substantial RNA contamination

Conceptually, $w_r$ should be a decreasing function of the latent RNA fraction
$\omega_r$:

$$
w_r = h(\omega_r)
$$

for some monotone decreasing map $h$, with $h(0) = 1$.

The statistical purpose of $w_r$ is not to create a new biological parameter.
It is to let calibration borrow information from many mostly-gDNA regions while
downweighting regions with stronger RNA evidence.

This is preferable to a hard-threshold rule because:

1. it respects that RNA-poor versus RNA-contaminated is a continuum
2. it avoids introducing arbitrary cutoffs
3. it lets the many low-expression regions in real data contribute partial
   information rather than being discarded

## 24. Soft Evidence for gDNA Purity

At theory level, the weight $w_r$ should be informed by the same three fragment
evidence channels already identified:

1. splicing evidence
2. density-contrast evidence
3. strand-asymmetry evidence

So one can write abstractly:

$$
w_r = W\left(E_r^{\mathrm{splice}},
E_r^{\mathrm{density}},
E_r^{\mathrm{strand}}\right)
$$

where:

- larger splice evidence should decrease $w_r$
- stronger RNA-like density contrast should decrease $w_r$
- stronger RNA-like strand asymmetry should decrease $w_r$

No specific functional form is needed at the theoretical stage. What matters is
that the calibration system is driven by fragment-level RNA evidence and that it
produces a continuous notion of gDNA purity.

## 25. Calibration Targets Share the Same Purity Problem

The same notion of RNA-poor regions should underlie all gDNA-specific nuisance
parameter estimation, including:

- gDNA strand-symmetry hyperparameter $\kappa_{\mathrm{sym}}$
- gDNA fragment-length distribution
- any future gDNA-specific positional or enrichment corrections

This unifies what might otherwise look like separate heuristics.

The theory is:

1. score regions for gDNA purity using fragment-level RNA evidence
2. let each region contribute to calibration in proportion to that purity
3. estimate gDNA-specific nuisance parameters from the weighted calibration set

So “RNA-poor region selection” is not a side issue. It is the central upstream
statistical problem that makes empirical Bayes calibration of gDNA possible.

## 26. Region Partition and Annotation Context

The calibration problem operates on genomic regions, each defined by:

$$
r = (\mathrm{ref}, \mathrm{start}, \mathrm{end})
$$

These regions are not arbitrary windows. They should arise from a partition of
the genome induced by annotation structure.

At theory level, one reference sequence is partitioned into disjoint intervals
whose annotation context is constant within each interval.

The relevant annotation axes are:

1. genic versus intergenic
2. exonic versus intronic membership
3. positive-strand versus negative-strand versus strand-ambiguous annotation

Because exons and introns may overlap across transcripts, the resulting context
classes are combinatorial. A convenient conceptual partition includes classes
such as:

- exon-only `+`
- exon-only `-`
- exon-only ambiguous
- intron-only `+`
- intron-only `-`
- intron-only ambiguous
- exon-intron `+`
- exon-intron `-`
- exon-intron ambiguous
- intergenic

The exact label set is less important than the principle: every genomic base in
the partition should have a well-defined annotation context that can be used to
interpret fragment density and strand asymmetry.

## 27. Admissible Fragments for Calibration

Calibration should be performed on fragments, not reads, and ideally on the
subset of fragments whose overlap with the annotation partition is itself
unambiguous.

For a fragment $f$, let $\mathcal{O}(f)$ denote the set of partition intervals it
overlaps. In the simplest theory version, retain only fragments for which the
context assignment is effectively unique.

This means holding out fragments whose annotation context is ambiguous because
they:

- overlap opposite-strand annotations simultaneously
- span multiple region classes in a way that prevents a clear regional tally
- lie in regions whose annotation context is itself strand-ambiguous

This is not a statement that ambiguous-context fragments are useless for the
main quantification model. It is only a statement that the calibration problem
benefits from a cleaner subset.

The clean theoretical calibration sample therefore consists of:

- unspliced fragments
- with unambiguous annotation context
- assigned to partition intervals with interpretable strand and density meaning

## 28. Region-Level Evidence Channels

For each calibration region $r$, define three theory-level RNA evidence scores.

### 28.1 Splice evidence score

Let:

$$
S_r = \text{number of spliced fragments associated with region } r
$$

Because splicing is effectively pure RNA evidence, $S_r$ should strongly reduce
the gDNA purity weight.

In a refined model, spliced evidence can be converted into an expected amount of
unspliced RNA in the same transcriptional neighborhood using transcript
geometry, exon structure, and fragment-length considerations. Conceptually,
spliced fragments are not only direct RNA evidence; they are anchors for the
amount of unspliced RNA one should expect nearby.

### 28.2 Density evidence score

For each region or region class, define fragment density as:

$$
D_r = \frac{F_r}{L_r^{\mathrm{eff}}}
$$

where:

- $F_r$ is the count of admissible unspliced fragments overlapping the region
- $L_r^{\mathrm{eff}}$ is the effective genomic length of the interval or region
   class

The density signal is then derived from contrasts between annotation contexts.

At the theory level, a particularly natural contrast is:

$$
E_r^{\mathrm{density}} = \Psi\!
\left(
\frac{D_r^{\mathrm{exonic}}}{D_r^{\mathrm{nonexonic}}}
\right)
$$

where nonexonic means intronic plus intergenic or, in practice, whichever
nonexonic contexts are sufficiently represented and reliable.

The rationale is exactly as you described: intergenic space may be sparse in
modern annotations, whereas intronic space often remains abundant and RNA-poor.
So exonic versus nonexonic density contrast is the more robust theoretical RNA
signal than exonic versus intergenic alone.

### 28.3 Strand evidence score

For a strand-resolved region $r$, let:

$$
N_r^+ = \text{positive-strand fragment count},
\qquad
N_r^- = \text{negative-strand fragment count},
\qquad
N_r = N_r^+ + N_r^-
$$

The strand evidence score should measure deviation from the symmetric gDNA
reference distribution.

At theory level, this should not be a raw imbalance statistic alone. It should
be a divergence from the symmetric Beta-Binomial family induced by the global
hyperparameter $\kappa_{\mathrm{sym}}$.

So one may write abstractly:

$$
E_r^{\mathrm{strand}} = \Delta\!
\left(N_r^+, N_r^- ; \kappa_{\mathrm{sym}}\right)
$$

where larger divergence means less gDNA-like and therefore lower purity.

## 29. Soft Purity Weight Revisited

With these ingredients, the gDNA-dominance weight becomes:

$$
w_r = W\!
\left(
E_r^{\mathrm{splice}},
E_r^{\mathrm{density}},
E_r^{\mathrm{strand}}
\right)
$$

with the qualitative behavior:

- more splice evidence lowers $w_r$
- stronger exonic-versus-nonexonic contrast lowers $w_r$
- stronger strand asymmetry relative to the gDNA symmetry model lowers $w_r$

The exact functional form of $W$ is an implementation question. The theoretical
requirement is only that it be monotone in the appropriate directions and yield
a continuous spectrum of gDNA purity.

## 30. The Circularity in Strand Calibration

There is a genuine apparent circularity:

- to compute the strand evidence score we would like to know
   $\kappa_{\mathrm{sym}}$
- but $\kappa_{\mathrm{sym}}$ itself is estimated from weighted calibration
   regions

This is not a flaw. It simply means the calibration problem is itself a small
hierarchical estimation problem.

The correct theoretical view is self-consistent empirical Bayes.

One clean formulation is alternating optimization:

1. start from a provisional symmetry model
2. compute provisional region weights $w_r$
3. estimate $\kappa_{\mathrm{sym}}$ from the weighted Beta-Binomial calibration
    objective
4. recompute strand evidence and weights under the updated
    $\kappa_{\mathrm{sym}}$
5. iterate to self-consistency

This is the same general logic as many empirical-Bayes procedures: the weights
and the hyperparameter are learned jointly by alternation, not by a one-shot
closed form.

Importantly, the circularity is mild because the splice and density channels do
not depend on $\kappa_{\mathrm{sym}}$. They anchor the purity model even before
the strand component is fully calibrated.

## 31. Probabilistic gDNA Purity Model

The cleanest theoretical formulation is to treat regional gDNA purity itself as
a latent probabilistic quantity rather than as a deterministic score.

For each calibration region $r$, introduce a latent gDNA purity variable:

$$
\pi_r \in [0,1]
$$

interpreted as the fraction of fragments in that region arising from gDNA.

Equivalently, one can introduce a binary latent region state

$$
u_r \in \{\mathrm{gDNA\text{-}dominant},\; \mathrm{RNA\text{-}contaminated}\}
$$

with uncertainty about that state. The continuous version is usually cleaner,
because it matches the biological reality that purity varies along a spectrum.

Under this view, the regional soft weight is not an ad hoc score. It is a
posterior quantity such as:

$$
w_r = \mathbb{E}[\pi_r \mid \text{regional evidence}]
$$

or, in a two-state simplification,

$$
w_r = \Pr(u_r = \mathrm{gDNA\text{-}dominant} \mid \text{regional evidence})
$$

This is the most coherent interpretation of a soft calibration weight.

## 32. Regional Evidence Model

Let the observed evidence for region $r$ be:

$$
y_r = \left(
S_r,
E_r^{\mathrm{density}},
N_r^+,
N_r^-
\right)
$$

where:

- $S_r$ is splice evidence
- $E_r^{\mathrm{density}}$ summarizes RNA-like density contrast
- $N_r^+, N_r^-$ are strand-resolved unspliced fragment counts

The theoretically optimal calibration model factorizes as:

$$
p(y_r \mid \pi_r, \kappa_{\mathrm{sym}}, \Theta_{\mathrm{cal}})
= p(S_r \mid \pi_r, \Theta_{\mathrm{cal}})
\; p(E_r^{\mathrm{density}} \mid \pi_r, \Theta_{\mathrm{cal}})
\; p(N_r^+, N_r^- \mid \pi_r, \kappa_{\mathrm{sym}}, \Theta_{\mathrm{cal}})
$$

where $\Theta_{\mathrm{cal}}$ collects any additional nuisance parameters for
the calibration model.

The exact parametric form of these conditional distributions is not required at
the theory stage. The essential qualitative structure is:

- larger $\pi_r$ should imply less splice evidence
- larger $\pi_r$ should imply weaker RNA-like density contrast
- larger $\pi_r$ should imply strand counts closer to the symmetric gDNA model

## 33. Posterior-Weighted Empirical Bayes

In the fully probabilistic view, the calibration problem is hierarchical:

1. each region has latent purity $\pi_r$
2. the strand behavior of the gDNA part is governed by the global hyperparameter
   $\kappa_{\mathrm{sym}}$
3. the evidence channels jointly update posterior beliefs about each region's
   purity

The ideal calibration weight is then a posterior expectation rather than a
deterministic hand-built rule.

That is,

$$
w_r = \mathbb{E}[\pi_r \mid y_r, \kappa_{\mathrm{sym}}, \Theta_{\mathrm{cal}}]
$$

and the global symmetry hyperparameter is estimated from a weighted or joint
likelihood over all regions.

This is theoretically optimal because it lets the model acknowledge uncertainty
about whether a region is truly gDNA-dominant.

## 34. Practical Relaxation

A full probabilistic purity model may be too expensive or too elaborate for the
first implementation.

So the right theoretical posture is:

1. the ideal object is posterior regional gDNA purity
2. in practice, one may approximate that posterior with a heuristic or
   semiparametric score built from the three evidence channels
3. any such heuristic should be judged by how well it approximates the posterior
   ordering of regions by gDNA purity

This gives a principled role for heuristics: they are approximations to a latent
probabilistic purity model, not replacements for the underlying statistical
concept.

## 35. Comparison of Calibration-Region Models

There are three natural modeling levels for regional gDNA purity.

### 35.1 Continuous latent purity model

Each region has:

$$
\pi_r \in [0,1]
$$

and the calibration weight is a posterior expectation:

$$
w_r = \mathbb{E}[\pi_r \mid y_r]
$$

Advantages:

- matches the biological continuum from pure gDNA to mixed RNA plus gDNA
- uses all regions without hard cutoffs
- gives a natural interpretation to soft weights
- integrates cleanly with empirical-Bayes hyperparameter learning

Disadvantages:

- requires the richest calibration submodel
- likely needs approximation in practice

### 35.2 Two-state mixture model

Each region has a latent discrete state such as:

$$
u_r \in \{\mathrm{gDNA\text{-}dominant}, \mathrm{RNA\text{-}contaminated}\}
$$

and the soft weight is:

$$
w_r = \Pr(u_r = \mathrm{gDNA\text{-}dominant} \mid y_r)
$$

Advantages:

- probabilistic and interpretable
- simpler than a continuous latent purity model
- can often capture most of the practical benefit

Disadvantages:

- imposes an artificial dichotomy on a genuinely continuous phenomenon
- less flexible when many regions are only mildly contaminated

### 35.3 Heuristic scoring model

Define a deterministic score from splice, density, and strand evidence and map
it to a weight.

Advantages:

- simplest to implement initially
- easiest to debug

Disadvantages:

- lacks a proper generative interpretation
- can drift into arbitrary tuning unless anchored to a probabilistic target

## 36. Final Theoretical Recommendation

The preferred target model is the continuous latent purity model.

That is:

1. each calibration region has a latent gDNA purity fraction $\pi_r$
2. the soft calibration weight is a posterior expectation of that purity
3. $\kappa_{\mathrm{sym}}$ and other gDNA nuisance parameters are learned by
   empirical Bayes from those posterior-weighted regions

For first implementation work, approximate inference or semiparametric
shortcuts are acceptable, but they should be viewed as approximations to this
continuous latent-purity target rather than as the model itself.
