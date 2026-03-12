# Rigel Compact EM Model

This note gives a compact mathematical model for Rigel while staying out of
implementation details.

It adopts the preferred `T + N + 2` modeling choice discussed in
`docs/rigel_revision_plan/THEORETICAL_MODEL.md`:

- two strand-specific gDNA components per locus
- one collapsed public gDNA abundance obtained by summing those two components
- transcript-level mature RNA components
- nRNA-span-level nascent RNA components

The purpose is to define the latent variables, likelihood, and EM objective as
cleanly as possible.

## 1. Locus-Level Latent Variables

Consider one locus $\ell$.

Let:

- $\mathcal{T}_\ell$ be the set of mature RNA transcript components
- $\mathcal{N}_\ell$ be the set of nRNA span components
- $g_{pos}$ and $g_{neg}$ be the two strand-specific gDNA components

The component set is:

$$
\mathcal{C}_\ell = \mathcal{T}_\ell \cup \mathcal{N}_\ell \cup \{g_{pos}, g_{neg}\}
$$

For each fragment $f$ assigned to this locus, define a latent source label:

$$
z_f \in \mathcal{C}_\ell
$$

Let the mixture weights be:

$$
   heta_\ell = \{\theta_c\}_{c \in \mathcal{C}_\ell},
\qquad
   heta_c \ge 0,
\qquad
\sum_{c \in \mathcal{C}_\ell} \theta_c = 1
$$

Define total gDNA mass and the implied strand fraction by:

$$
G_\ell = \theta_{g_{pos}} + \theta_{g_{neg}}
$$

and, when $G_\ell > 0$,

$$
\phi_\ell = \frac{\theta_{g_{pos}}}{G_\ell}
$$

The public gDNA abundance reported for the locus is $G_\ell$, while the EM
internally carries the pair $(\theta_{g_{pos}}, \theta_{g_{neg}})$.

## 2. Observed Fragment Features

For each fragment $f$, let the observed feature bundle be:

$$
x_f = (s_f, j_f, a_f, l_f)
$$

where:

- $s_f \in \{\mathrm{pos}, \mathrm{neg}\}$ is observed genomic strand orientation
- $j_f$ is splice class
- $a_f$ is alignment or compatibility information
- $l_f$ is fragment-length information

## 3. Likelihood Factorization

For any component $c$, define:

$$
p(x_f \mid z_f = c)
= p(s_f \mid c) \, p(j_f \mid c) \, p(a_f \mid c) \, p(l_f \mid c)
$$

It is useful to separate strand from the rest:

$$
\lambda_{fc}
= p(j_f \mid c) \, p(a_f \mid c) \, p(l_f \mid c)
$$

so that:

$$
p(x_f \mid z_f = c) = p(s_f \mid c) \, \lambda_{fc}
$$

### 3.1 RNA strand likelihood

For mature RNA transcript $t$ and nRNA span $n$, the strand term is determined by
the RNA strand model:

$$
p(s_f \mid t) = \rho_t(s_f),
\qquad
p(s_f \mid n) = \rho_n(s_f)
$$

where the expected strand gets high probability and the opposite strand gets the
protocol-error probability.

### 3.2 gDNA strand likelihood

For the two gDNA components:

$$
p(s_f \mid g_{pos}) = \rho_{g_{pos}}(s_f),
\qquad
p(s_f \mid g_{neg}) = \rho_{g_{neg}}(s_f)
$$

where in the idealized model:

$$
\rho_{g_{pos}}(\mathrm{pos}) = 1,
\quad
\rho_{g_{pos}}(\mathrm{neg}) = 0,
\quad
\rho_{g_{neg}}(\mathrm{pos}) = 0,
\quad
\rho_{g_{neg}}(\mathrm{neg}) = 1
$$

Exact zeroes may be replaced by tiny numerical error probabilities if the
implementation requires it.

## 4. Mixture Likelihood

Given locus mixture parameters $\theta_\ell$, the observed-data likelihood for
one fragment is:

$$
p(x_f \mid \theta_\ell)
= \sum_{c \in \mathcal{C}_\ell} \theta_c \, p(x_f \mid z_f = c)
$$

Equivalently:

$$
p(x_f \mid \theta_\ell)
= \sum_{t \in \mathcal{T}_\ell} \theta_t \, \rho_t(s_f) \, \lambda_{ft}
+ \sum_{n \in \mathcal{N}_\ell} \theta_n \, \rho_n(s_f) \, \lambda_{fn}
+ \theta_{g_{pos}} \, \rho_{g_{pos}}(s_f) \, \lambda_{f g_{pos}}
+ \theta_{g_{neg}} \, \rho_{g_{neg}}(s_f) \, \lambda_{f g_{neg}}
$$

## 5. Priors

The compact model keeps only three kinds of regularization.

### 5.1 Weak simplex prior on component weights

Use a Dirichlet prior on the locus mixture:

$$
   heta_\ell \sim \operatorname{Dirichlet}(\alpha_\ell)
$$

For the redesign target, optional components should use a sparse MAP regime
with:

$$
0 < \alpha_c < 1
$$

so that unsupported components can be driven to the simplex boundary rather
than merely retaining diffuse fractional mass.

This is the generic occupancy prior.

### 5.2 nRNA fraction prior

For each nRNA span $n$, let $\mathcal{T}(n)$ be the transcripts sharing that
nRNA span.

Define mature mass in that group:

$$
M_n(\theta) = \sum_{t \in \mathcal{T}(n)} \theta_t
$$

and group total RNA mass:

$$
S_n(\theta) = \theta_n + M_n(\theta)
$$

Define the nascent fraction:

$$
\eta_n(\theta) = \frac{\theta_n}{S_n(\theta)}
$$

Place a Beta prior on this interpretable quantity:

$$
\eta_n \sim \operatorname{Beta}(a_n, b_n)
$$

This is the principled way to couple mature and nascent RNA inside a shared
transcriptional group.

### 5.3 gDNA strand-symmetry prior

Put a symmetric Beta prior on the implied locus gDNA strand fraction:

$$
\phi_\ell \sim \operatorname{Beta}\left(\frac{\kappa_{\mathrm{sym}}}{2},
\frac{\kappa_{\mathrm{sym}}}{2}\right)
$$

with:

$$
\phi_\ell = \frac{\theta_{g_{pos}}}{\theta_{g_{pos}} + \theta_{g_{neg}}}
$$

This gives:

- mean $1/2$
- a weak spring toward symmetry when $\kappa_{\mathrm{sym}}$ is small
- a strong spring when $\kappa_{\mathrm{sym}}$ is large

Integrating over the locus-specific strand split gives the correct overdispersed
count model for gDNA strand balance. If the total gDNA count at locus $\ell$ is
$N_\ell^{(g)}$, then the pos-strand count satisfies:

$$
K_\ell^{pos} \mid N_\ell^{(g)}
\sim \operatorname{BetaBinomial}\left(
N_\ell^{(g)}, \frac{\kappa_{\mathrm{sym}}}{2}, \frac{\kappa_{\mathrm{sym}}}{2}
\right)
$$

with mean:

$$
\mathbb{E}[K_\ell^{pos} \mid N_\ell^{(g)}] = \frac{N_\ell^{(g)}}{2}
$$

and variance:

$$
\operatorname{Var}(K_\ell^{pos} \mid N_\ell^{(g)})
= \frac{N_\ell^{(g)}\left(N_\ell^{(g)} + \kappa_{\mathrm{sym}}\right)}
{4(\kappa_{\mathrm{sym}} + 1)}
$$

As $\kappa_{\mathrm{sym}} \to \infty$, this reduces to the Binomial variance
$N_\ell^{(g)}/4$. Smaller $\kappa_{\mathrm{sym}}$ preserves symmetry in
expectation while allowing more locus-to-locus strand imbalance.

## 6. MAP Objective

For one locus, the MAP objective is:

$$
\mathcal{L}(\theta)
= \sum_f \log \left( \sum_{c \in \mathcal{C}_\ell}
+ \theta_c \, p(x_f \mid z_f = c) \right)
+ \log p(\theta)
+ \sum_{n \in \mathcal{N}_\ell} \log p(\eta_n(\theta))
+ \log p(\phi(\theta))
$$

Expanding the prior terms and dropping constants gives:

$$
\mathcal{L}(\theta)
= \sum_f \log \left( \sum_{c} \theta_c \, p(x_f \mid c) \right)
+ \sum_c (\alpha_c - 1) \log \theta_c
+ \sum_n \left[(a_n - 1) \log \eta_n(\theta) + (b_n - 1) \log (1 - \eta_n(\theta))\right]
+ \left(\frac{\kappa_{\mathrm{sym}}}{2} - 1\right)
\left[\log \phi(\theta) + \log (1 - \phi(\theta))\right]
$$

subject to:

$$
   heta_c \ge 0,
\qquad
\sum_c \theta_c = 1
$$

## 7. Equivalent Form of the nRNA Prior Term

Because:

$$
\eta_n(\theta) = \frac{\theta_n}{S_n(\theta)},
\qquad
1 - \eta_n(\theta) = \frac{M_n(\theta)}{S_n(\theta)}
$$

the nRNA Beta prior contributes:

$$
\log p(\eta_n(\theta))
= (a_n - 1) \log \theta_n
+ (b_n - 1) \log M_n(\theta)
- (a_n + b_n - 2) \log S_n(\theta)
+ \mathrm{const}
$$

This makes clear that the prior acts on the mature-versus-nascent partition of
the shared group, not on the nRNA component in isolation.

## 8. EM Formulation

Introduce the posterior responsibility for fragment $f$ and component $c$:

$$
r_{fc} = p(z_f = c \mid x_f, \theta)
$$

### 8.1 E-step

Given current parameters $\theta^{old}$:

$$
r_{fc}
= \frac{\theta_c^{old} \, p(x_f \mid c)}
{\sum_{c'} \theta_{c'}^{old} \, p(x_f \mid c')}
$$

### 8.2 Expected complete-data objective

Ignoring constants that do not depend on $\theta$, the EM auxiliary function is:

$$
Q(\theta)
= \sum_f \sum_c r_{fc} \log \theta_c
+ \sum_c (\alpha_c - 1) \log \theta_c
+ \sum_n \log p(\eta_n(\theta))
+ \log p(\phi(\theta))
$$

The structural terms $\lambda_{fc}$ do not appear in the M-step because they are
fixed with respect to $\theta$.

Combining terms gives:

$$
Q(\theta)
= \sum_c (R_c + \alpha_c - 1) \log \theta_c
+ \sum_n \log p(\eta_n(\theta))
+ \log p(\phi(\theta))
$$

where:

$$
R_c = \sum_f r_{fc}
$$

## 9. M-step Structure for the gDNA Pair

Define gDNA responsibility totals:

$$
R_{g_{pos}} = \sum_f r_{f g_{pos}},
\qquad
R_{g_{neg}} = \sum_f r_{f g_{neg}}
$$

Then the gDNA-specific behavior is represented inside the same simplex as the
rest of the mixture. There is no separate nuisance-parameter emission update for
collapsed gDNA. After the weight update, the implied strand fraction is:

$$
\phi_\ell^{new}
= \frac{\theta_{g_{pos}}^{new}}{\theta_{g_{pos}}^{new} + \theta_{g_{neg}}^{new}}
$$

If we temporarily ignore the nRNA coupling terms and focus only on the gDNA
pair given total gDNA mass, the Beta prior implies the familiar conditional MAP
update:

$$
\phi_\ell^{new}
= \frac{R_{g_{pos}} + \frac{\kappa_{\mathrm{sym}}}{2} - 1}
{R_{g_{pos}} + R_{g_{neg}} + \kappa_{\mathrm{sym}} - 2}
$$

provided the Beta prior parameters are greater than $1$. The key point is that
this condition is now a property of the gDNA pair inside the ordinary mixture,
not a separate collapsed-gDNA emission update.

If we want the strictly symmetric model with no strand adaptation, we can
constrain:

$$
   heta_{g_{pos}} = \theta_{g_{neg}}
$$

which corresponds to the fixed-symmetry submodel.

### 9.1 Relationship to the current targeted-excess penalty

The current native solver applies a targeted-excess penalty rather than a fully
explicit gDNA pair update. That penalty decomposes gDNA mass into:

- a symmetric protected baseline
- an asymmetric excess subject to discount

The discount factor has the form:

$$
w_{\mathrm{sym}} \propto \left[4 \hat{p}(1 - \hat{p})\right]^{\kappa/2 - 1}
$$

which is the kernel of a symmetric Beta density. So the current implementation
is already an approximation to the same strand-symmetry prior family proposed
here.

What changes in the redesign is not the biological regularizer. What changes is
the parameterization and the explicitness with which that regularizer is
represented.

### 9.2 Identifiability caveat and asymmetric constraint

There is a real identifiability risk at loci with both gDNA and nascent RNA.

- nRNA contributes sense-dominant unspliced signal
- a flexible gDNA pair can absorb part of that same signal if unconstrained

So a purely symmetric Beta prior is not necessarily sufficient as the full
practical story. A principled first implementation may need an asymmetric
constraint or evidence gate that treats sense-heavy gDNA as more suspect than
antisense-heavy gDNA, because sense-heavy excess is more confounded with nRNA.

## 10. M-step for the Mixture Weights

The M-step for $\theta$ is the constrained optimization problem:

$$
\max_{\theta \in \Delta}
\left[
\sum_c (R_c + \alpha_c - 1) \log \theta_c
+ \sum_n \log p(\eta_n(\theta))
+ \log p(\phi(\theta))
\right]
$$

where $\Delta$ is the simplex:

$$
\Delta = \left\{\theta : \theta_c \ge 0, \sum_c \theta_c = 1\right\}
$$

Without the nRNA coupling prior and the gDNA symmetry prior, this reduces to the
usual normalized-count Dirichlet update:

$$
   heta_c^{new}
\propto R_c + \alpha_c - 1
$$

With the nRNA Beta priors included, the M-step is no longer fully separable,
because each nRNA component is coupled to the mature transcript masses in its
shared group. The gDNA symmetry prior similarly couples the pair
$(\theta_{g_{pos}}, \theta_{g_{neg}})$ through the induced strand fraction
$\phi(\theta)$.

That coupling is intentional: it is exactly the biological prior that the model
is trying to encode.

## 11. The Compact Model in Words

The redesign target can be summarized as follows.

For each locus:

1. infer a simplex over mature RNA transcripts, nRNA spans, and two
   strand-specific gDNA components
2. score each fragment by a product of strand, splice, alignment, and
   fragment-length likelihood terms
3. regularize the mature-versus-nascent split through a Beta prior on each nRNA
   group's nascent fraction
4. regularize gDNA strand symmetry through a Beta prior on the implied strand
   fraction of the gDNA pair

This keeps the public output simple:

- transcript-level mature RNA
- nRNA-span-level or transcript-fanned-out nascent RNA
- locus-level total gDNA

while using a cleaner internal parameterization for the EM.

## 12. Recommended Simplest Version

Among all theory-consistent variants, the simplest version is:

1. two internal gDNA components per locus, one per strand
2. report their sum as the public locus-level gDNA abundance
3. optionally fix the gDNA pair to exact symmetry at first
4. if the data show extra locus-level strand variability, free the implied
   strand fraction under the symmetric Beta prior
5. keep the nRNA prior only as a prior on group nascent fractions, not as an
   additional heuristic initializer

That gives a clean progression:

- simplest exact symmetric gDNA-pair model
- then one extra degree of freedom for strand asymmetry if the data justify it

## 13. Hyperparameter Strategy for $\kappa_{\mathrm{sym}}$

The preferred theoretical treatment is empirical Bayes with hierarchical sharing
across loci.

That is, within the main locus model the implied strand fraction obeys:

$$
\phi_\ell \mid \kappa_{\mathrm{sym}}
\sim \operatorname{Beta}\left(\frac{\kappa_{\mathrm{sym}}}{2},
\frac{\kappa_{\mathrm{sym}}}{2}\right)
$$

with:

$$
\phi_\ell = \frac{\theta_{g_{pos}}}{\theta_{g_{pos}} + \theta_{g_{neg}}}
$$

but we do not regard $\kappa_{\mathrm{sym}}$ as fixed a priori. Instead it is a
global library-level hyperparameter estimated from a calibration collection of
gDNA-dominant regions.

The calibration likelihood is:

$$
\prod_{r \in \mathcal{R}}
\Pr\left(K_r^+ \mid N_r, \kappa_{\mathrm{sym}}\right)
$$

where each factor is Beta-Binomial:

$$
K_r^+ \mid N_r, \kappa_{\mathrm{sym}}
\sim \operatorname{BetaBinomial}\left(
N_r, \frac{\kappa_{\mathrm{sym}}}{2}, \frac{\kappa_{\mathrm{sym}}}{2}
\right)
$$

Conceptually, this means:

- each calibration region has its own latent strand fraction
- all such regions share one global dispersion parameter
- the estimated $\kappa_{\mathrm{sym}}$ is then recycled as the prior strength
  on the gDNA pair in the locus EM

The same calibration framework should also estimate the gDNA fragment-length
distribution. These are joint nuisance targets learned from the same
purity-weighted regional evidence, not two unrelated upstream procedures.

## 14. What Counts as a Calibration Region

In theory, a calibration region is any genomic region for which RNA abundance is
negligible and residual fragments can therefore be treated as gDNA-dominant.

This may include:

- intergenic regions
- intronic regions of unexpressed genes
- exonic regions of genes or transcripts with effectively zero expression

The important theoretical point is that purity matters more than annotation
class. In targeted or capture-enriched libraries, intergenic regions may be too
sparse to support stable estimation, so low-RNA exonic or intronic regions are a
better conceptual target.

More precisely, calibration regions should be defined in terms of fragment-level
RNA evidence. The three strongest signals are:

1. spliced fragments, which are effectively pure RNA evidence
2. density contrast between exon, intron, and intergenic space
3. strand asymmetry in stranded libraries

So the calibration problem is not “find intergenic regions”, but rather “find
regions whose fragment populations are overwhelmingly gDNA-like under those
three evidence channels”.

This same calibration set should in principle support both:

- estimation of $\kappa_{\mathrm{sym}}$
- estimation of the gDNA fragment-length model

## 15. Soft Calibration Weights

The preferred theoretical formulation is soft weighting rather than hard region
selection.

For each candidate calibration region $r$, define a gDNA-dominance weight:

$$
w_r \in [0,1]
$$

interpreted as the degree to which the region's fragment population is believed
to be explained by gDNA rather than RNA.

The weight should be informed by the same fragment-level RNA evidence channels:

- spliced-fragment evidence
- exon versus intron or intergenic density contrast
- strand-asymmetry evidence

So in abstract form:

$$
w_r = W\left(E_r^{\mathrm{splice}},
E_r^{\mathrm{density}},
E_r^{\mathrm{strand}}\right)
$$

where stronger RNA evidence drives the weight downward.

The theoretical role of these weights is straightforward:

- regions with $w_r \approx 1$ contribute strongly to gDNA calibration
- regions with intermediate $w_r$ contribute weakly but are not discarded
- regions with $w_r \approx 0$ contribute negligibly

This avoids introducing hard-threshold calibration rules while still exploiting
the fact that many real genomic regions are nearly unexpressed.

## 16. Weighted Empirical-Bayes Calibration for $\kappa_{\mathrm{sym}}$

With soft weights, the calibration objective for the global symmetry
hyperparameter becomes a weighted Beta-Binomial likelihood:

$$
\mathcal{L}_{\mathrm{cal}}(\kappa_{\mathrm{sym}})
= \sum_{r \in \mathcal{R}} w_r \,
\log \Pr\left(K_r^+ \mid N_r, \kappa_{\mathrm{sym}}\right)
$$

with:

$$
K_r^+ \mid N_r, \kappa_{\mathrm{sym}} 
\sim \operatorname{BetaBinomial}\left(
N_r,
\frac{\kappa_{\mathrm{sym}}}{2},
\frac{\kappa_{\mathrm{sym}}}{2}
\right)
$$

This is the cleanest theory-level empirical-Bayes formulation:

1. all regions may contribute
2. mostly gDNA regions contribute most strongly
3. the global hyperparameter is estimated from a continuous purity spectrum

The same weighted calibration principle should also apply to learning the gDNA
fragment-length distribution.

## 17. Region-Based Soft Calibration Model

Let the genome be partitioned into annotation-context intervals:

$$
r = (\mathrm{ref}, \mathrm{start}, \mathrm{end}, \mathrm{context})
$$

where context records exon, intron, intergenic, and strand-resolved annotation
status.

For calibration, focus on admissible unspliced fragments whose overlap with this
partition is itself unambiguous. For each region $r$, define:

- $S_r$: spliced-fragment evidence in the region or neighborhood
- $D_r^{\mathrm{exonic}}$: exonic-context density
- $D_r^{\mathrm{nonexonic}}$: nonexonic density
- $N_r^+, N_r^-$: strand-resolved unspliced fragment counts

Then define a latent RNA contamination fraction $\omega_r$ and a derived
gDNA-dominance weight $w_r = h(\omega_r)$ with $w_r \in [0,1]$.

At theory level, the weight is determined by the evidence channels through a
monotone map:

$$
w_r = W\!
\left(
S_r,
\Psi\!\left(D_r^{\mathrm{exonic}} / D_r^{\mathrm{nonexonic}}\right),
\Delta\!\left(N_r^+, N_r^- ; \kappa_{\mathrm{sym}}\right)
\right)
$$

where:

- larger $S_r$ lowers $w_r$
- larger exonic enrichment lowers $w_r$
- greater divergence from symmetric strand behavior lowers $w_r$

## 18. Weighted Calibration Objective with Self-Consistency

Using these soft weights, the empirical-Bayes calibration objective for the
global symmetry hyperparameter is:

$$
\mathcal{L}_{\mathrm{cal}}(\kappa_{\mathrm{sym}})
= \sum_{r \in \mathcal{R}} w_r(\kappa_{\mathrm{sym}}) \,
\log \Pr\left(K_r^+ \mid N_r, \kappa_{\mathrm{sym}}\right)
$$

where:

$$
K_r^+ = N_r^+,
\qquad
N_r = N_r^+ + N_r^-
$$

and the likelihood factor is Beta-Binomial.

Because the weights themselves may depend on
$\kappa_{\mathrm{sym}}$ through the strand-divergence score, the clean theory is
an alternating empirical-Bayes procedure:

1. initialize a provisional symmetry model
2. compute soft purity weights over regions
3. update $\kappa_{\mathrm{sym}}$ from the weighted Beta-Binomial objective
4. recompute the strand evidence and weights
5. iterate to self-consistency

This avoids hard thresholds while resolving the apparent circularity in a
statistically coherent way.

## 19. Probabilistic Interpretation of the Soft Weights

The soft weights are best understood as posterior gDNA-purity quantities, not as
arbitrary region scores.

For each region $r$, introduce a latent purity variable $\pi_r \in [0,1]$ and
regional evidence:

$$
y_r = (S_r, E_r^{\mathrm{density}}, N_r^+, N_r^-)
$$

The theoretically optimal weight is then:

$$
w_r = \mathbb{E}[\pi_r \mid y_r, \kappa_{\mathrm{sym}}, \Theta_{\mathrm{cal}}]
$$

or, under a two-state approximation,

$$
w_r = \Pr(u_r = \mathrm{gDNA\text{-}dominant} \mid y_r)
$$

This interpretation makes the weighted empirical-Bayes calibration objective a
posterior-weighted likelihood rather than a purely handcrafted scoring rule.

## 20. Theory Versus Approximation

The full joint calibration model would infer simultaneously:

- regional purity variables
- the global symmetry hyperparameter $\kappa_{\mathrm{sym}}$
- any additional nuisance parameters governing splice or density evidence

In practice, that may be approximated by a simpler weight function built from
the three evidence channels.

The correct theoretical standard remains the same:

- the heuristic weight should approximate posterior gDNA purity
- its purpose is to rank and weight regions by likely gDNA dominance
- it should be evaluated against the latent probabilistic target, not viewed as
   an unrelated ad hoc rule

## 21. Recommended Calibration Submodel

Among the candidate calibration formulations, the preferred first
implementation is:

1. a two-state mixture with latent region state
   $u_r \in \{\mathrm{gDNA\text{-}dominant},\mathrm{RNA\text{-}contaminated}\}$
2. posterior calibration weights
   $w_r = \Pr(u_r = \mathrm{gDNA\text{-}dominant} \mid y_r)$
3. joint empirical-Bayes estimation of global gDNA nuisance parameters,
   especially $\kappa_{\mathrm{sym}}$ and the gDNA fragment-length
   distribution, from those posterior-weighted regions

The continuous latent-purity model remains the cleaner long-run theoretical
target. Pure heuristic scoring is acceptable only as an approximation to the
posterior purity ordering implied by the probabilistic model.
