# 01 — Generative Model

The model is stated at the level of **per-region (and per-boundary)
sufficient statistics**. There is no per-fragment latent variable in
the calibrator; that level of resolution belongs downstream, at the
locus EM, where it is already implemented.

The goal of this document is to derive every per-region likelihood
factor from a generative story, then read off the sufficient
statistics, then state the inference target. The result is a model
with two per-region latents ($\pi_r^{(g)}, \omega_r$) and five
library-wide hyperparameters ($\rho_0, \phi, \rho_d^{\text{BB}}, \rho_r^{\text{BB}}, \epsilon_s$).

Note on what is *not* a free hyperparameter: the gDNA sense-strand
mean is **biologically and molecularly fixed at $\kappa_d = 0.5$** —
genomic DNA is unstranded by construction (double-stranded molecule
sheared and converted to dsDNA library), so the sense/antisense
balance is exactly 0.5 in expectation. What varies is the
*overdispersion* around 0.5, captured by $\rho_d^{\text{BB}}$. The
per-region RNA sense-strand mean $\kappa_r^{\text{rna}}$ is also not
fitted in calibration — it is **pre-computed from the initial BAM
scan** by `StrandModel` using spliced unique mappers, where the
gDNA/RNA ambiguity does not exist. Only the RNA strand
overdispersion $\rho_r^{\text{BB}}$ is fitted.

## 1. Notation

| Symbol | Meaning |
|---|---|
| $r$ | Region index, $r \in \{1, \ldots, R\}$. |
| $b$ | Boundary index. Two per region: $b \in \{L_r, R_r\}$. |
| $L_r^{\text{eff}}$ | Effective length of region $r$ for gDNA accessibility. v1: physical length. |
| $n_r^{\text{u}}$ | Unspliced fragment count in region $r$ (or boundary $b$). |
| $n_r^{\text{s}}$ | Spliced fragment count. |
| $k_r^{(+)}$ | Sense-strand count among unspliced fragments. |
| $\kappa_d$ | gDNA sense-strand mean. **Fixed at 0.5** (biology). Not a free parameter. |
| $\rho_d^{\text{BB}}$ | gDNA strand Beta-Binomial dispersion, $\rho_d^{\text{BB}} \in (0, 1)$. Fitted. |
| $\kappa_r^{\text{rna}}$ | RNA sense-strand mean for region $r$. **Pre-computed by `StrandModel`** from spliced uniques during the BAM scan. Not fitted in calibration. |
| $\rho_r^{\text{BB}}$ | RNA strand Beta-Binomial dispersion, $\rho_r^{\text{BB}} \in (0, 1)$. Fitted. |
| $\epsilon_s$ | gDNA splice-artifact rate. Very small. |
| $\rho_0$ | Library-wide gDNA fragment density (fragments per bp). |
| $\phi$ | Library-wide overdispersion. Unifies NB count overdispersion and Gamma exposure spread (§4.2). |
| $\omega_r$ | Per-region gDNA exposure factor, multiplicative on $\rho_0$. |
| $\pi_r^{(g)}$ | Per-region gDNA mixing proportion: $\Pr(\text{fragment is gDNA} \mid r)$. |
| $\mu_r^{(g)}$ | $\omega_r \rho_0 L_r^{\text{eff}}$. |
| $\mu_r^{(d)}$ | Per-region RNA mean count, inferred. |

## 2. Why model per region, not per fragment

The calibrator's job is to inform the downstream locus prior with
deconvoluted mass and exposure. The locus prior consumes per-region
scalars, not per-fragment posteriors. Operating per-fragment at the
calibration stage would (a) duplicate work the downstream fragment
scorer already performs more precisely, (b) require an additional
streaming pass over the entire fragment buffer, and (c) prevent the
natural Beta-Binomial / Negative-Binomial treatment of channel
overdispersion that only makes sense at the region scale.

The sufficient statistics in §1 are exactly what an NB-count +
BB-strand + Poisson-spliced model needs. Everything else is wasted
resolution.

> **No FL channel.** Earlier drafts modelled a per-region unspliced
> FL histogram $h_r[\ell]$ as a two-component mixture multinomial.
> That term has been removed: spliced fragments are FL-ambiguous
> (compatible with many isoforms), and the calibrator's count + strand
> + spliced channels carry enough signal without paying the substrate
> cost of per-region/per-boundary FL histograms. The full-resolution
> FL signal is still used per-fragment by the downstream EM scorer
> ([`scoring.py`](../../src/rigel/scoring.py),
> [`frag_length_model.py`](../../src/rigel/frag_length_model.py)).

## 3. Spliced fragments are deterministic RNA evidence

A spliced fragment is generated only by an RNA molecule traversing a
splice junction — gDNA cannot produce one except via alignment
artifact at rate $\epsilon_s \approx 10^{-4}$ – $10^{-3}$. The
*rate* at which a given RNA region produces spliced fragments depends
on transcript geometry (number of introns, exon sizes, read length) —
this is a **transcript-specific geometric quantity, not a library
biochemical hyperparameter**. There is no scalar "RNA splice rate"
that the calibrator can or should estimate.

Treatment in the calibrator: the spliced count $n_r^{\text{s}}$ is
routed directly into the RNA mass:

$$
M_r^{(d, \text{spliced})} = (1 - \epsilon_s) \cdot n_r^{\text{s}}
\quad\approx\quad n_r^{\text{s}}.
$$

The $\epsilon_s$ correction is a **failsafe**: the upstream BAM
scanner already runs splice-artifact detection and removal on
fragments (`bam_scanner.cpp` discards alignments that pattern-match
known artifact signatures before they reach the substrate aggregates).
By the time spliced counts enter calibration, the residual artifact
rate should already be tiny; modelling $\epsilon_s$ in the calibrator
guards against any residual mis-classifications and lets the
estimator quantify them rather than silently absorb them. In practice
$\epsilon_s$ contributes essentially nothing. The spliced count does
not enter the per-region mixing-proportion update; it is added to
$M_r^{(d)}$ after the unspliced-channel deconvolution. (See §6 for
how this composes.)

## 4. The unspliced-channel generative model

The unspliced channel carries all the ambiguity that calibration
exists to resolve. It has two observed sufficient statistics
($n_r^{\text{u}}$, $k_r^{(+)}$) per region, modelled as follows.

### 4.1 The mixture latent

A fragment in region $r$'s unspliced channel is gDNA with probability
$\pi_r^{(g)}$ and RNA with probability $1 - \pi_r^{(g)}$. The mixing
proportion is itself derived (rather than free) from the per-region
means:

$$
\pi_r^{(g)} = \frac{\mu_r^{(g)}}{\mu_r^{(g)} + \mu_r^{(d)}}.
$$

This identification ties the mixing proportion to the exposure-driven
count means rather than treating it as a free Beta-distributed latent.

### 4.2 Count channel: Negative-Binomial via Gamma-Poisson identity

$$
n_r^{\text{u}} \sim \mathrm{NB}\!\bigl(\mu_r^{(g)} + \mu_r^{(d)},\, \phi\bigr).
$$

Parametrization: mean $\mu$, dispersion $\phi$ such that $\mathrm{Var}(n) = \mu + \phi \mu^2$.

The NB has the equivalent representation as a Gamma-mixed Poisson:

$$
\omega_r \sim \mathrm{Gamma}(1/\phi,\, 1/\phi) \qquad \text{(mean 1, variance } \phi\text{)},
$$

$$
n_r^{(g)} \mid \omega_r \sim \mathrm{Poisson}(\omega_r \cdot \rho_0 \cdot L_r^{\text{eff}}).
$$

So **$\phi$ governs both the NB count overdispersion and the prior
spread of the per-region exposure factor $\omega_r$ around its
library mean of 1**. This is the elegant unification: rather than
estimating an exposure-variance hyperparameter $\tau^2$ separately
from a count-dispersion hyperparameter, both are a single scalar
$\phi$ with one global M-step.

**Consequence (closed-form exposure posterior).** Given the
soft-allocated per-region gDNA mass $M_r^{(g)}$ (computed in the
E-step; §6), the conjugate Gamma-Poisson posterior is

$$
\omega_r \mid M_r^{(g)},\, \rho_0,\, \phi \;\sim\; \mathrm{Gamma}\!\Bigl(\frac{1}{\phi} + M_r^{(g)},\; \frac{1}{\phi} + \rho_0 L_r^{\text{eff}}\Bigr),
$$

with point estimate (posterior mean)

$$
\hat{\omega}_r = \frac{1/\phi + M_r^{(g)}}{1/\phi + \rho_0 L_r^{\text{eff}}},
$$

and delta-method log-variance

$$
\hat{\sigma}_{\log \omega_r}^2 = \frac{1}{1/\phi + M_r^{(g)}}.
$$

Closed form. No Newton, no Laplace, no quadrature.

### 4.3 Strand channel: Beta-Binomial for both gDNA and RNA

**gDNA side.** The gDNA sense-strand mean is biologically fixed at
$\kappa_d = 0.5$ (genomic DNA is double-stranded by construction).
The strand count is overdispersed Beta-Binomial:

$$
K_r^{(g,+)} \mid N_r^{(g)} \sim \mathrm{BetaBinom}(N_r^{(g)},\, \kappa_d = 0.5,\, \rho_d^{\text{BB}})
$$

with intra-class correlation $\rho_d^{\text{BB}} \in (0, 1)$ the only
free gDNA strand parameter. Variance is

$$
\mathrm{Var}\bigl(K_r^{(g,+)}\bigr) = N_r^{(g)} \cdot 0.25 \cdot \bigl[1 + (N_r^{(g)} - 1)\rho_d^{\text{BB}}\bigr].
$$

Equivalent Beta-mixed Binomial representation: $\kappa_r^{(g)} \sim \mathrm{Beta}(\alpha_d, \alpha_d)$ with $\alpha_d = 0.5 \cdot (1 - \rho_d^{\text{BB}}) / \rho_d^{\text{BB}}$ (symmetric about 0.5).

**RNA side.** RNA libraries are also overdispersed around their
strand-specificity mean — perfectly stranded library targets give
mean 1.0 (or 0.0), unstranded targets give 0.5, mixed protocols give
something in between, but in every case there is biological and
technical variability that Binomial under-models. RNA strand is
also Beta-Binomial:

$$
K_r^{(d,+)} \mid N_r^{(d)} \sim \mathrm{BetaBinom}(N_r^{(d)},\, \kappa_r^{\text{rna}},\, \rho_r^{\text{BB}})
$$

where $\kappa_r^{\text{rna}}$ is **pre-computed** by `StrandModel`
during the BAM scan from spliced unique mappers (which are
unambiguous RNA), and $\rho_r^{\text{BB}}$ is the only free RNA
strand parameter (single global scalar, fitted in the calibrator's
M-step). For unstranded libraries $\kappa_r^{\text{rna}} = 0.5$, the
RNA and gDNA strand likelihoods become similar (both have mean 0.5)
and the strand channel contributes near-zero log-evidence-ratio —
handled automatically by the model.

**Observed sense count** decomposes as

$$
k_r^{(+)} = K_r^{(g, +)} + K_r^{(d, +)}
$$

with the per-region soft-allocation E-step computing the expected split (§6).

> **Strand-ambiguous regions (added 2026-05-29; see
> [`../acc_caljointmodel/00_implementation_plan.md`](../acc_caljointmodel/00_implementation_plan.md)
> §4 D7).** This strand factor assumes the region has a single transcript
> strand, so that "sense" is defined. A region with overlapping transcripts on
> **both** strands (signature with both `+` and `−` bits) has **no valid sense
> split** — every read is sense for one transcript and antisense for the other
> — and the strand factor is **omitted** for it. Such regions are deconvolved
> by the count channel + a boundary-density sweep + a global gDNA-density
> fallback (D7). Intergenic (no transcript) regions are different: gDNA is
> unstranded, so an arbitrary sense is harmless and the strand factor stays
> (contributing ~zero log-evidence). This paragraph is the interim contract;
> §4.3/§7 get a full rewrite in Phase 8.

### 4.4 Splice channel (already covered in §3)

The spliced count contributes mass $(1 - \epsilon_s) n_r^{\text{s}}$ directly to $M_r^{(d)}$. There is no Bayes step; the spliced count is a sufficient statistic that is added at the end of the unspliced-channel deconvolution.

## 5. The per-region likelihood (composition)

Combining §4.2–4.4, the per-region log-likelihood factors as

$$
\log p\bigl(n_r^{\text{u}}, k_r^{(+)}, n_r^{\text{s}} \mid \pi_r^{(g)}, \omega_r, \theta\bigr) \;=\; \mathcal{L}_r^{\text{count}} + \mathcal{L}_r^{\text{strand}} + \mathcal{L}_r^{\text{spliced}}
$$

with

$$
\mathcal{L}_r^{\text{count}} = \log \mathrm{NB}\bigl(n_r^{\text{u}} \mid \mu_r^{(g)} + \mu_r^{(d)}, \phi\bigr)
$$

$$
\mathcal{L}_r^{\text{strand}} = \log\bigl[\text{BetaBinom mixture marginal at } k_r^{(+)}\bigr]
$$

$$
\mathcal{L}_r^{\text{spliced}} = \log\bigl[\text{Poisson}(n_r^{\text{s}} \mid (1-\epsilon_s) \mu_r^{(d)} + \epsilon_s \mu_r^{(g)} \cdot \text{(geometric term)})\bigr]
$$

The spliced-channel likelihood term in v1 is **not** used to update
the latents (since the per-transcript splice geometry is unknown to
the calibrator); spliced fragments enter only through the determined
$M_r^{(d, \text{spliced})}$ pathway. The strand channel mixture marginal
is approximated in inference by per-channel BetaBinom evaluations
under the soft allocation (see [`03_inference.md`](03_inference.md)
§3); this is exact in the limit of well-resolved channels and
empirically very tight even at small counts.

## 6. Per-region E-step: soft allocation summary

The per-region E-step computes $(M_r^{(g)}, M_r^{(d)})$ from the
sufficient statistics. The count and strand channels combine
through log Bayes-factor addition into a single per-region scalar
$\pi_r^{(g)}$, then the unspliced count is split:

$$
M_r^{(g, \text{unspl})} = n_r^{\text{u}} \cdot \pi_r^{(g)}
$$

where $\pi_r^{(g)}$ is updated from the count and strand log-evidence
ratios. Details in [`03_inference.md`](03_inference.md) §3.

Total per-region masses (unspliced + spliced):

$$
\boxed{\;M_r^{(d, \text{cont})} = (n_r^{\text{u}} - M_r^{(g, \text{unspl})}) + (1 - \epsilon_s)\, n_r^{\text{s}}, \qquad M_r^{(g, \text{cont})} = M_r^{(g, \text{unspl})} + \epsilon_s \, n_r^{\text{s}}\;}
$$

Boundary equivalents follow the same formulas with boundary sufficient statistics.

## 7. Per-region exposure posterior (G4)

After the soft allocation, define the total per-region gDNA mass used
for exposure inference:

$$
M_r^{(g, \text{tot})} = M_r^{(g, \text{cont})} + \tfrac{1}{2}\bigl(M_r^{(g, L)} + M_r^{(g, R)}\bigr).
$$

The half-split prevents double-counting boundary mass across adjacent
regions. The Gamma posterior on $\omega_r$ from §4.2 is then:

$$
\boxed{\;\omega_r \mid M_r^{(g, \text{tot})} \;\sim\; \mathrm{Gamma}\!\Bigl(\frac{1}{\phi} + M_r^{(g, \text{tot})},\; \frac{1}{\phi} + \rho_0 L_r^{\text{eff}}\Bigr).\;}
$$

Posterior mean and log-variance as in §4.2.

## 8. Library hyperparameters

Five global scalars:

| Param | Role | M-step approach |
|---|---|---|
| $\rho_0$ | Library mean gDNA density | Closed form: $\sum_r M_r^{(g, \text{tot})} / \sum_r \hat{\omega}_r L_r^{\text{eff}}$ |
| $\phi$ | NB count dispersion = Gamma exposure prior shape | 1-D Newton on NB log-likelihood gradient (single global scalar; moment estimator warm start) |
| $\rho_d^{\text{BB}}$ | gDNA strand BB dispersion (mean fixed at 0.5) | 1-D Newton on BB log-likelihood gradient (single global scalar; moment estimator warm start) |
| $\rho_r^{\text{BB}}$ | RNA strand BB dispersion (mean $\kappa_r^{\text{rna}}$ pre-computed by `StrandModel`) | 1-D Newton on BB log-likelihood gradient (single global scalar; moment estimator warm start) |
| $\epsilon_s$ | gDNA splice-artifact rate | **Deferred (PR 7)** — not modelled; spliced fragments are RNA, artifacts filtered upstream by the `alignable` splice blacklist (see 03_inference §5.5) |

All five M-steps in [`03_inference.md`](03_inference.md) §4.

## 9. What the model does NOT include (v1)

| Excluded | Reason |
|---|---|
| Per-fragment latent variables in the calibrator | Calibration operates on sufficient statistics. Fragment-level evidence is downstream EM's job. |
| RNA splice-rate hyperparameter $q$ or $\bar{q}$ | Transcript-geometric, not library-biochemical. Spliced fragments are deterministic RNA evidence. |
| Two-state expression indicator $z_r$ | The mass is the output; the state was a means in earlier drafts. |
| Per-region RNA exposure factor | The RNA mass is the output; per-transcript RNA factors are downstream EM's job. |
| Per-region gDNA FL pmf | FL is not a calibration channel; see \u00a72 note. |
| Per-region strand BB dispersion | Single library-wide $\rho_d^{\text{BB}}$ and $\rho_r^{\text{BB}}$ are sufficient. |
| Free gDNA strand mean $\kappa_d$ | Biologically fixed at 0.5; estimating it would invite identifiability problems with regional RNA contamination. |
| Free per-region RNA strand mean | Pre-computed by `StrandModel` from spliced uniques during the scan; the calibrator must not re-fit it. |
| Mappability-aware $L_r^{\text{eff}}$ | v1 uses physical length. v2 may wire in mappability tracks. |

## 10. Why this kills the paralog phantom-RNA pathway

Recap from [`02_failure_audit.md`](02_failure_audit.md) §2: legacy
calibration manufactured ~26 units of RNA mass on a paralog region
from a noise-driven 126/26 strand asymmetry of pair-anchored gDNA
straddlers, cascading into a `(4, 638)` EM split.

In the new model (no FL channel), the rescue argument rests on three
legs that together — not individually — force the gDNA explanation:

1. **Count channel (NB).** The paralog region's unspliced count
   ($n_r^{\text{u}} \approx 152$) is *consistent with the gDNA-only
   mean* $\mu_r^{(g)} = \hat{\omega}_r \rho_0 L_r^{\text{eff}} \approx 150$
   under the Gamma exposure posterior pooled across the surrounding
   intergenic regions. The NB log-likelihood at $\mu = 150, \phi$
   evaluated at $n = 152$ is essentially at the mode; admitting an
   extra $\mu_r^{(d)}$ would push the mean off-mode for negligible
   gain. Score: gDNA-favored.
2. **Strand channel (BB).** A 126/26 split under
   $\mathrm{BetaBinom}(N=152, \kappa_d = 0.5, \rho_d^{\text{BB}})$
   with the library-fit $\rho_d^{\text{BB}}$ has a finite tail at
   this asymmetry: the BB log-likelihood under fixed $\kappa_d = 0.5$
   has enough overdispersion mass at moderate asymmetry that 126/26
   is not strong evidence for RNA. Score: near-zero log-evidence-ratio
   (the legacy Binomial mis-read this as decisive RNA signal).
3. **Spliced channel.** Zero spliced fragments on the region.
   $n_r^{\text{s}} = 0$ contributes no RNA mass via the deterministic
   pathway. Score: gDNA-favored.

The combined per-region mass allocation is $M_r^{(g, \text{cont})} \approx 145$,
$M_r^{(d, \text{cont})} \approx 5$ — weaker than the prior (count + FL +
strand) rescue but architecturally cleaner. The exposure update is
the conjugate Gamma posterior on $\omega_r$; the locus prior gets
near-symmetric RNA pseudocounts across paralog transcripts. The EM
splits multimapper mass near-symmetrically.

> **Risk flag.** Without the FL channel, paralog rescue depends on
> all three of {count NB matching gDNA mean, BB strand absorbing the
> asymmetry, zero spliced support} firing together. If any one of
> these breaks (e.g. a paralog region happens to have $n_r^{\text{u}}$
> well above the gDNA-only mean, or a single misclassified spliced
> fragment, or unusually low $\rho_d^{\text{BB}}$ from a poorly-mixed
> library) the rescue could fail in ways the FL-bearing model would
> have caught. Synthetic + scenario benchmarks must compare this
> model's paralog behavior against the FL-bearing baseline before
> production cutover. See [`02_failure_audit.md`](02_failure_audit.md)
> §2.5 and [`06_validation_plan.md`](06_validation_plan.md).
