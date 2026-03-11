# Rigel Theoretical Model

This document is intentionally theory-first. It is not a description of the
current code. Its purpose is to define the mathematical problem Rigel is trying
to solve, starting from biology and probability rather than from existing
implementation choices.

The main use of this note is to provide a clean target for redesigning the EM.

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

That does not automatically mean the implementation must literally use two gDNA
components. It means the theoretical model should first acknowledge that a DNA
fragment is generated from one strand or the other, never both at once.

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

## 10. Theoretical Recommendation on gDNA

From first principles, the cleanest model is:

1. two explicit gDNA strand components per locus
2. deterministic or near-deterministic strand likelihood for each fragment
3. a symmetric prior on the gDNA strand fraction

This is the most biologically faithful and most interpretable formulation.

However, there is an important caveat:

### The benefit of two gDNA strands is not per-fragment competition. The benefit is a cleaner parameterization of locus-level strand symmetry.

Per fragment, once the observed strand is known, the incompatible gDNA strand is
effectively absent. The competition is not between $d_{\ell,+}$ and $d_{\ell,-}$
for the same fragment. The competition is between RNA and the compatible DNA
strand, while the two DNA strands are coupled only through the prior on
$\phi_\ell$.

That means the case for two gDNA components depends on whether we want an
explicit strand-balance latent variable. If yes, two components are the natural
representation. If no, a collapsed single-component gDNA model with a $0.5$
strand term remains mathematically coherent.

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
5. If interpretability and clean probabilistic structure are the priority, the
   two-strand gDNA model is the better starting point for redesign.

## 17. Immediate Design Implication

Before discussing implementation, the core theoretical decision is this:

### Do we want gDNA to be represented as an explicit strand-resolved latent variable, or only as a collapsed total nuisance component?

If we want explicit strand symmetry as part of the model, then the clean answer
is to represent:

$$
(G_\ell, \phi_\ell)
$$

or equivalently:

$$
(D_{\ell,+}, D_{\ell,-})
$$

If we only care about total gDNA mass, then one collapsed component is still
theoretically sound, but its strand term should remain $0.5$.
