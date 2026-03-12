# Rigel EM Redesign: Design Review

This is a design review of the four interconnected documents proposing a
first-principles redesign of Rigel's EM initialization and prior system:

- `THEORETICAL_MODEL.md` — mathematical foundations
- `THEORETICAL_EM_OBJECTIVE.md` — compact EM formulation
- `EM_PARAMETER_AUDIT.md` — current system audit
- `FIRST_PRINCIPLES_IMPLEMENTATION_PLAN.md` — implementation roadmap

The critique is grounded in the actual codebase, including the C++ EM solver,
Python prior/locus pipeline, fragment-length models, scoring infrastructure,
and annotation index structures.

---

## 1. What the Documents Get Right

### 1.1 The gDNA strand model derivation is airtight

The theoretical derivation of the two-strand gDNA model is the strongest
part of the work. The key conclusions are:

- Biology dictates two gDNA strand sources per locus ($g_{pos}$, $g_{neg}$)
- A fragment from a given DNA strand has deterministic strand likelihood
  ($p = 1$ for the matching strand, $p = 0$ for the opposite)
- Symmetry is expressed as a prior on the strand balance, not in the
  per-fragment likelihood

The collapsed one-component model with strand likelihood 0.5 is correctly
identified as the marginal of the symmetric two-strand model. However, since
we are adopting T+N+2 with two explicit gDNA strand components (see Section
2.1), the collapsed model is no longer the implementation target — the
explicit two-strand formulation with deterministic strand likelihoods is.

### 1.2 Separating likelihood from prior from initialization

The documents correctly identify that the current system conflates three
conceptually distinct roles: likelihood terms, regularization priors, and
numerical initialization. The EM_PARAMETER_AUDIT is particularly valuable —
it maps the actual code to these categories honestly.

### 1.3 Soft purity weighting over hard thresholds

The argument for continuous gDNA-purity weights over binary region selection
is well-motivated. The current system uses intergenic density as a single
gDNA anchor, which is fragile in capture-enriched or annotation-dense
genomes. Soft weighting from multiple evidence channels is a genuine
improvement.

### 1.4 The three-channel evidence framework

Identifying splice evidence, density contrast, and strand asymmetry as the
three independent RNA-contamination signals is correct and complete. These
are the right axes.

### 1.5 The phased build order is pragmatic

Starting with region partition, then evidence extraction, then purity model,
then EB calibration, then EM integration is the right dependency order.

---

## 2. Significant Concerns

### 2.1 T+N+1 vs T+N+2: the architecture favors two explicit gDNA strands

The documents commit to T+N+1 (one collapsed gDNA component with free
$\phi_\ell$). After examining the actual solver architecture, **T+N+2 (two
explicit gDNA strand components) is the better choice for practical reasons
the theoretical docs don't address.**

**The problem with T+N+1 and free $\phi_\ell$:**

The current scoring pre-computes gDNA log-likelihoods as:

    gdna_ll = LOG_HALF + gdna_fl + gdna_log_sp

(see `scoring.cpp:1250`). This `LOG_HALF` is the collapsed symmetric strand
term, baked into the stored log-lik. If $\phi_\ell$ is free, it changes each
EM iteration, so the gDNA log-lik would need to change too. But the
equivalence class architecture (`em_solver.cpp:133-277`) pre-computes and
caches log-liks. To make $\phi_\ell$ free in T+N+1, you'd need to either:

- Re-score gDNA every iteration (expensive, breaks pre-computation)
- Store gDNA log-liks without strand term and apply $\log(\phi)$ or
  $\log(1-\phi)$ dynamically in the E-step kernel (complicates the tight
  inner loop)

**T+N+2 avoids this entirely:**

With two explicit gDNA components $g_{pos}$ and $g_{neg}$:

- $g_{pos}$ log-lik for sense fragments = `gdna_fl + gdna_log_sp` (since
  $\log(1) = 0$)
- $g_{pos}$ log-lik for antisense fragments = $-\infty$
- $g_{neg}$ log-lik: vice versa
- These are ALL pre-computable and iteration-independent
- The E-step kernel stays completely standard
- The $-\infty$ entries get zeroed by the existing early-skip logic (NEON
  path)

**M-step coupling is straightforward:**

After standard responsibility accumulation:

$$
R_{g_{pos}} = \sum_f r_{f,g_{pos}}, \qquad R_{g_{neg}} = \sum_f r_{f,g_{neg}}
$$

$$
R_g = R_{g_{pos}} + R_{g_{neg}}
$$

$$
\phi = \frac{R_{g_{pos}} + \kappa_{\mathrm{sym}}/2 - 1}
{R_g + \kappa_{\mathrm{sym}} - 2}
$$

$$
  heta_{g_{pos}} = \phi \cdot \theta_g^{\mathrm{total}}, \qquad
  heta_{g_{neg}} = (1 - \phi) \cdot \theta_g^{\mathrm{total}}
$$

This is the same closed-form update as the $\phi_\ell$ M-step in the
theoretical docs.

**SQUAREM benefits:**

In T+N+2, both $g_{pos}$ and $g_{neg}$ weights are in the state vector that SQUAREM
accelerates. In T+N+1, $\phi_\ell$ is a separate parameter outside SQUAREM's
view. T+N+2 should converge faster.

**Equivalence class cost:**

T+N+2 splits ECs by gDNA strand (sense units get $g_{pos}$ candidate, antisense
get $g_{neg}$). This creates more ECs but each is smaller. Net computational cost
is similar since total work is $O(NK)$ regardless of grouping.

**The asymmetric $\phi$ constraint maps cleanly to T+N+2:**

Instead of constraining $\phi_\ell$, constrain $g_{pos}$ relative to $g_{neg}$. For
example, cap $\theta_{g_{pos}}$ at some multiple of $\theta_{g_{neg}}$ based on
$\kappa_{\mathrm{sym}}$. This is the same math, just expressed on mixture
weights rather than on a separate parameter.

**Recommendation:** Adopt T+N+2 as the target architecture. It is more
natural for the existing solver, avoids the $\phi_\ell$-in-likelihood
complication, plays better with SQUAREM, and makes the asymmetric constraint
a simple weight inequality. The collapsed one-gDNA output is trivially
recovered as $\theta_g = \theta_{g_{pos}} + \theta_{g_{neg}}$.

### 2.1b The current penalty's relationship to the proposed model

The current targeted-excess penalty (`em_solver.cpp:592-644`) approximates a
MAP update for $\phi_\ell$ under the proposed Beta prior:

- It decomposes gDNA EM counts into symmetric baseline ($2 \min$) +
  asymmetric excess
- The discount $w_{\mathrm{sym}} = \mathrm{balance}^{\kappa/2 - 1}$ with
  $\mathrm{balance} = 4 \hat{p}(1 - \hat{p})$ is the kernel of a
  $\operatorname{Beta}(\kappa/2, \kappa/2)$ density

This correspondence is never acknowledged in the documents. Whether T+N+1
or T+N+2 is adopted, the bridge between the current penalty and the proposed
Beta prior should be documented.

### 2.2 Identifiability risk with free $\phi_\ell$

If $\phi_\ell$ is free to move away from 0.5, there is a serious
identifiability issue at loci with both gDNA and nascent RNA. Both gDNA and
nRNA produce unspliced fragments. nRNA is strand-asymmetric (dominated by
the sense strand). If $\phi_\ell$ can shift toward the sense strand, the EM
can absorb some nRNA signal into gDNA by adjusting $\phi_\ell$, especially
when the nRNA Beta prior is weak.

**This is a confirmed real problem, not just a theoretical concern.**

The documents do not discuss this interaction at all. The current targeted
penalty avoids this by construction: it only discounts excess, never inflates
the gDNA sense fraction. A free $\phi_\ell$ does not have this guarantee.

**Recommendation:** Analyze the $\phi_\ell$/nRNA interaction formally.
Consider whether the $\operatorname{Beta}(\kappa_{\mathrm{sym}}/2,
\kappa_{\mathrm{sym}}/2)$ prior is strong enough to prevent $\phi_\ell$ from
drifting toward the dominant RNA strand at high-expression loci. If not, the
"one-sided" nature of the current penalty may be a feature, not a bug. See
Section 5 for detailed analysis.

### 2.3 Fragment-length estimation shares the same calibration problem

The documents treat FL estimation and purity estimation as separate concerns.
But the current FL model (`frag_length_model.py:470-583`) has a known
fragility: it relies on intergenic fragments for gDNA FL training, and modern
annotations leave extremely sparse intergenic space. The FL model needs
"RNA-poor" regions — the same thing the purity model identifies.

This creates a unification opportunity: the calibration framework should
estimate BOTH $\kappa_{\mathrm{sym}}$ AND the gDNA FL distribution from the
same purity-weighted regions. The bootstrap is:

1. Three-channel purity weights (splice, density, strand) — no FL needed
2. Estimate gDNA FL from purity-weighted regions (both strands)
3. Optionally: FL becomes a fourth purity signal, iterate to convergence

This directly addresses the known FL fragility while also making FL
available as a calibration signal in later iterations.

**Recommendation:** Explicitly include gDNA FL estimation as a co-target
of the calibration framework (alongside $\kappa_{\mathrm{sym}}$), not as a
separate upstream step. The documents' Workstream E mentions FL but doesn't
emphasize that FL estimation and purity estimation are coupled.

### 2.4 The calibration region partition is underspecified

Workstream A proposes a region table with ~12 annotation fields. But the
documents don't address the combinatorial explosion from overlapping genes on
opposite strands. In a typical mammalian genome, antisense transcription
creates many regions where both strands have exonic annotation. The proposed
`strand_context` and `is_context_ambiguous` fields hint at this but don't
resolve it.

More importantly, the existing `TranscriptIndex._gen_genomic_intervals()`
already produces a genome-wide tiling into EXON, TRANSCRIPT, INTERGENIC, and
UNAMBIG_INTRON intervals with per-transcript associations. This is a strong
starting point, but the documents don't reference it or discuss how to extend
it.

**Recommendation:** Ground Workstream A in the existing interval
infrastructure. Define the calibration partition as a refinement of the
existing tiling, not as a new parallel system. Use the four-flag annotation
scheme (see Section 4.2) to handle overlapping antisense annotations cleanly.

**Caveat:** Refactoring the index to a four-flag system would cascade through
the codebase (cgranges queries, fragment resolver, buffer structures, locus
construction). This is a major refactoring effort that should be scoped
carefully. It may be better to build the calibration partition as a separate
post-index-build step initially, using the four-flag scheme for calibration
regions only, and defer the full index refactoring to a later phase.

### 2.5 The continuous purity model may not be worth the complexity

The documents strongly advocate for a continuous latent purity variable
$\pi_r \in [0,1]$ as the ideal model, with the two-state mixture and
heuristic scoring as fallbacks. But the marginal value of continuous purity
over two-state is unclear for this application.

The calibration problem is inherently well-separated: most genomic regions
are either clearly gDNA-dominant (intergenic, unexpressed intronic) or
clearly RNA-contaminated (expressed exonic). The transition zone is narrow.
A two-state mixture with posterior state probabilities as weights would
capture the vast majority of practical benefit while being much simpler to
implement, debug, and validate.

**Recommendation:** Commit to the two-state mixture as the first
implementation target, not just as a fallback. Reserve the continuous model
for a future iteration if the two-state model proves insufficient.

### 2.6 Self-consistent EB iteration lacks practical specification

The alternating EB procedure for $\kappa_{\mathrm{sym}}$ is theoretically
sound but practically underspecified:

- How many iterations are needed? (Likely 2-3, but this should be stated.)
- What is the initialization of the provisional symmetry model?
- What convergence criterion should be used?
- What if the initial estimate is very wrong — is the procedure contractive?

The circularity is described as "mild" because splice and density channels
don't depend on $\kappa_{\mathrm{sym}}$. This is true, but the strand
channel may dominate in low-splice regions (intergenic), which are exactly
the regions that matter most for calibration.

**Recommendation:** Specify the initialization (e.g.,
$\kappa_{\mathrm{sym}} = 100$ as a weakly informative start), the
convergence criterion (e.g.,
$|\Delta\kappa| / \kappa < 0.01$), and the maximum iteration count.

### 2.7 The nRNA prior system is left untouched

The redesign focuses almost entirely on gDNA. But the EM_PARAMETER_AUDIT
reveals that nRNA initialization is equally problematic: `nrna_init` is a
noisy strand-excess estimate that gates component eligibility, and the
hierarchical Beta prior on nascent fraction is a complex three-level EB
system with its own fragilities.

The documents explicitly scope nRNA out, but the gDNA and nRNA systems
interact. If $\phi_\ell$ becomes free, it directly affects the nRNA
eligibility gate (which depends on strand excess). If the calibration stage
changes how gDNA density is estimated, it changes the intergenic density
input to the hybrid gDNA estimator, which feeds into `gdna_init`, which
gates the gDNA component, which changes how much signal is available for
nRNA.

**Recommendation:** At minimum, document the expected interactions between
the gDNA redesign and the nRNA prior system. Identify which nRNA parameters
will need to be re-tuned after the gDNA changes.

### 2.8 Computational cost is not discussed

The proposed calibration stage adds significant work before the main EM.
The current pipeline does a single BAM pass followed by direct EM. The
redesign inserts an entire calibration phase.

**Recommendation:** Design the calibration to piggyback on the existing
single BAM pass. The fragment buffer already contains `splice_type`,
`exon_strand`, `sj_strand`, and per-transcript overlap statistics. Region
assignment and evidence accumulation should be a post-scan aggregation step,
not a second pass.

---

## 3. Structural Suggestions

### 3.1 Write a bridge document

The biggest gap is between the EM_PARAMETER_AUDIT (what exists) and the
THEORETICAL_MODEL (what should exist). A bridge document should explicitly
map each current mechanism to its theoretical counterpart:

| Current mechanism | Theoretical counterpart | Gap |
|---|---|---|
| targeted excess penalty | Beta prior on $\phi_\ell$ | approximation gap |
| `gdna_init` eligibility gate | gDNA abundance prior | binary vs continuous |
| `nrna_init` eligibility gate | nRNA abundance prior | binary vs continuous |
| intergenic density anchor | calibration region purity model | single-point vs regional |
| OVR warm start + prior | Dirichlet simplex prior | dual role conflation |

### 3.2 Define concrete success metrics before building

The validation plan lists categories but no specific metrics. Before
building, define:

- What accuracy metric on synthetic data constitutes "does not regress"?
- What diagnostic on real data constitutes "improves interpretability"?
- What runtime overhead is acceptable?

---

## 4. Ideas for Further Consideration

### 4.1 FL estimation and purity estimation are the same problem

The current FL model has a known issue: it relies on intergenic fragments
for gDNA FL training, but modern annotations leave extremely sparse
intergenic space. The FL model needs "RNA-poor" regions to train on — which
is exactly what the purity model identifies.

This means FL estimation is not a *signal* for calibration in the first
pass. It is a *target* of calibration, alongside $\kappa_{\mathrm{sym}}$.
The bootstrap:

1. Use splice + density + strand (three channels) to compute initial purity
   weights — no FL needed
2. Estimate gDNA FL distribution from those purity-weighted regions (using
   both sense and antisense fragments from RNA-poor regions, not just
   intergenic)
3. Optionally: now FL is available as a fourth signal for refined purity
   weights
4. Re-estimate FL and $\kappa_{\mathrm{sym}}$ from updated weights
5. Iterate to self-consistency

This is a unification: the same purity-weighted calibration framework solves
both the $\kappa_{\mathrm{sym}}$ estimation problem AND the FL model training
problem. The current FL model's fragility (sparse intergenic space) and the
current gDNA density estimator's fragility (same sparse intergenic anchor)
share a common root cause that the calibration framework addresses.

### 4.2 Four-flag annotation partition

Instead of the complex annotation context enum proposed in the documents
(Section 6.1), use four boolean flags per genomic interval:

- `has_exon_pos` (yes/no)
- `has_exon_neg` (yes/no)
- `has_intron_pos` (yes/no)
- `has_intron_neg` (yes/no)

This gives 16 possible combinations. All-zero = intergenic. Various flag
combinations describe different overlapping-transcript scenarios. This is
simpler and more composable than the current IntervalType enum system, and
naturally handles antisense overlaps without an "ambiguous" catchall.

The existing `TranscriptIndex._gen_genomic_intervals()` tiling could be
refactored to produce four-flag intervals instead of typed intervals.

### 4.3 Exploit the BAM scanner's existing resolution context

The C++ fragment resolver (`resolve_context.h`) already computes per-fragment
`exon_bp`, `intron_bp`, `unambig_intron_bp`, and transcript overlap
statistics. These are stored in the fragment buffer. Region assignment and
evidence accumulation could be done entirely as a post-scan aggregation over
the buffer, with no changes to the scanner.

### 4.4 Consider whether $\kappa_{\mathrm{sym}}$ should vary by reference

The current EB hierarchy for gDNA density already has a per-reference level.
If gDNA strand balance varies systematically by chromosome (e.g., due to
replication timing or chromatin structure), a per-reference
$\kappa_{\mathrm{sym}}$ might be more appropriate than a single global value.

### 4.5 Asymmetric $\phi_\ell$ prior for known-bias scenarios

Some library preparations have known strand biases in gDNA (e.g., from
tagmentation orientation preferences). The symmetric Beta prior on
$\phi_\ell$ assumes no such bias. If the calibration data show a systematic
global shift away from 0.5, the prior could be extended to
$\operatorname{Beta}(a, b)$ with $a \ne b$.

### 4.6 Joint calibration of FL and $\kappa_{\mathrm{sym}}$

The gDNA fragment-length distribution and strand symmetry hyperparameter are
currently treated as independent calibration targets. But they share the same
purity-weighted calibration set. A joint estimation procedure that alternates
FL estimation and $\kappa_{\mathrm{sym}}$ estimation within the same EB loop
would be cleaner and potentially more stable.

---

## 5. The $\phi_\ell$ / nRNA Identifiability Problem

**This is confirmed as a real issue, not just a theoretical concern.**

The problem: at high-expression loci, if $\phi_\ell$ is free, the EM can
shift $\phi_\ell$ toward the dominant RNA strand, absorbing nRNA signal into
gDNA. This happens because:

- nRNA produces strand-asymmetric unspliced fragments (sense-dominated)
- gDNA with free $\phi_\ell$ can mimic this by shifting $\phi_\ell$ toward
  sense
- the EM sees this as a better local fit because it avoids the nRNA Beta
  prior penalty
- result: gDNA inflated, nRNA deflated, transcript quantification biased

The current targeted-excess penalty avoids this by construction: it protects
the symmetric baseline and only discounts the asymmetric excess. It never
*increases* the sense fraction of gDNA. This one-sided property is a feature.

### Implications for the redesign

A naive free $\phi_\ell$ under symmetric Beta prior does NOT have this
one-sided guarantee. The MAP update:

$$
\phi_{\ell}^{\mathrm{new}}
= \frac{R_{g+} + \kappa_{\mathrm{sym}}/2 - 1}
{R_{g+} + R_{g-} + \kappa_{\mathrm{sym}} - 2}
$$

will shift $\phi$ toward whichever strand has more gDNA-attributed mass. If
the EM has incorrectly attributed some sense nRNA to gDNA, $\phi_\ell$ will
chase the error.

### Possible solutions

1. **Constrained $\phi_\ell$**: Clamp $\phi_\ell$ to
   $[0.5 - \delta, 0.5 + \delta]$ where $\delta$ is informed by
   $\kappa_{\mathrm{sym}}$ from calibration.

2. **Asymmetric penalty structure**: Keep the theoretical Beta prior but
   implement the M-step update with the "protect symmetric baseline"
   structure from the current code.

3. **Joint $\phi_\ell$ / $\eta_n$ update**: Couple the gDNA strand fraction
   and nRNA nascent fraction updates, since they share a finite pool of
   unspliced sense-strand fragments.

4. **Evidence-gated $\phi_\ell$ freedom**: Allow $\phi_\ell$ to deviate from
   0.5 only when there is strong positive evidence for asymmetry that cannot
   be explained by RNA.

**Recommendation:** Solution 4 is the most principled. In a stranded library,
sense-strand excess at a locus is ambiguous between nRNA and gDNA asymmetry.
Antisense excess is much more likely to be genuine gDNA asymmetry. So in
T+N+2 terms:

- Allow $\theta_{g_{neg}} > \theta_{g_{pos}}$ freely (antisense-heavy gDNA)
- Constrain $\theta_{g_{pos}} > \theta_{g_{neg}}$ more strongly (sense-heavy could be
  nRNA leakage)

This asymmetric constraint is biologically motivated and would be a genuine
theoretical improvement over the current system.

---

## 6. Revised Recommendations Given User Priorities

Given that the goal is both accuracy AND theoretical cleanliness, and that
the $\phi_\ell$/nRNA interaction is a confirmed real problem, and that
runtime is secondary to correctness:

### Phase 0: Bridge analysis (before any code)

- Formally derive the relationship between the current targeted-excess
  penalty and the proposed Beta-prior MAP update
- Quantify the approximation gap on synthetic data
- Document the one-sided property and why it matters for nRNA

### Phase 1: Soft-weighted calibration regions (three-channel bootstrap)

- Build on existing `TranscriptIndex` interval tiling
- Implement two-state mixture purity model (not continuous — save that)
- Use three evidence channels initially: splice, density, strand
- FL is NOT a signal here — it is a calibration target (see Section 4.1)

### Phase 2: Joint EB calibration of $\kappa_{\mathrm{sym}}$ AND gDNA FL

- Estimate $\kappa_{\mathrm{sym}}$ from weighted Beta-Binomial objective
- Estimate gDNA FL distribution from purity-weighted regions using BOTH
  strands (not just intergenic fragments)
- This directly fixes the known FL model fragility (sparse intergenic space)
- Optionally iterate: FL becomes a fourth purity signal, update weights,
  re-estimate, converge

### Phase 3: T+N+2 EM refactor with asymmetric strand constraint

- Adopt two explicit gDNA strand components ($g_{pos}$, $g_{neg}$)
- Pre-compute deterministic strand likelihoods (no $\phi_\ell$ in E-step)
- M-step: couple $g_{pos}$ and $g_{neg}$ via Beta prior on their ratio
- Implement asymmetric constraint: allow $\theta_{g_{neg}} > \theta_{g_{pos}}$ freely
  (antisense excess = genuine gDNA asymmetry), constrain
  $\theta_{g_{pos}} > \theta_{g_{neg}}$ more tightly (sense excess = ambiguous with
  nRNA leakage)
- Feed calibrated $\kappa_{\mathrm{sym}}$ from Phase 2 into the Beta prior
- Output: collapse $\theta_{g_{pos}} + \theta_{g_{neg}}$ to total gDNA for public
  interface

### Phase 4: Validation and ablation

- Ablation studies: old vs new, each component individually
- Synthetic data: vary gDNA fraction, strand symmetry, nRNA abundance
- Real data: compare $\phi$ distributions, $\kappa_{\mathrm{sym}}$ estimates,
  gDNA FL
- Specific test: high-expression loci where nRNA/gDNA identifiability is a
  known problem

---

## 7. Summary Assessment

The theoretical framework is rigorous and the implementation plan is
well-sequenced. The confirmed $\phi_\ell$/nRNA identifiability issue is the
most important design constraint and should drive the strand-symmetry model
design. The main conclusions are:

1. **T+N+2 is better than T+N+1** for the existing solver architecture
   (avoids $\phi_\ell$-in-likelihood complication, plays better with SQUAREM,
   deterministic pre-computable strand likelihoods)
2. **Naive free $\phi_\ell$ will make the identifiability problem worse**
   unless the update has an asymmetric constraint (allow antisense excess
   freely, constrain sense excess tightly)
3. **FL estimation and purity estimation are coupled** — solve both with
   the same calibration framework
4. **The continuous purity model is premature** — start with two-state
5. **Four-flag annotation partition** (`has_exon_pos/neg`,
   `has_intron_pos/neg`) is simpler and more composable than the proposed
   context enum

The documents provide an excellent theoretical foundation. The main gaps are:

- **FL and purity estimation are coupled** — the calibration framework
  should solve both simultaneously, not treat FL as a separate upstream step
- **The $\phi_\ell$/nRNA interaction needs an asymmetric constraint**
  informed by biology (antisense excess is diagnostic, sense excess is
  ambiguous)
- **The bridge between current code and proposed theory** needs to be made
  explicit before committing to the full redesign
- **Two-state mixture is the right starting purity model**, not continuous
