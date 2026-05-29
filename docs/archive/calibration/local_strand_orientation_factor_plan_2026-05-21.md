# Local Strand-Orientation Factor For gDNA Leakage Control

Date: 2026-05-21

Status: concept and implementation plan

Primary target failure: true gDNA fragments in exon-rich regions being assigned to
mRNA or nRNA because a locally symmetric gDNA pileup can be split into two
strand-specific RNA explanations.

## Executive Recommendation

The outside reviewer's proposal is worth developing, but only in its genuinely
local, aggregated form. A new per-fragment strand term by itself would not solve
the problem. Rigel already has per-fragment strand likelihoods in the EM score:
for gDNA the strand factor is symmetric, while strand-compatible RNA states get
the strand-specific advantage. The failure occurs because the current objective
lets plus-orientation fragments choose one RNA state and minus-orientation
fragments choose another RNA state independently.

The promising idea is therefore not "add strand again." It is to ask, for each
local genomic tile, whether the tile's aggregate orientation pattern is better
explained by one symmetric gDNA source or by multiple coordinated RNA sources.
That is a model-selection question over a group of fragments. Integrating over
local RNA mixture proportions gives an Occam penalty for the RNA explanation:
two RNA states can fit a 50/50 pileup, but only by tuning two local abundances to
match what one gDNA source explains naturally.

Recommended path:

1. Build a diagnostic-only local Bayes factor over tiles and validate that it
   separates true-gDNA hotspots from true-RNA tiles in the VCaP mixture.
2. Add a conservative static EM regularizer by distributing each tile-level log
   Bayes factor across that tile's eligible unspliced gDNA candidates. This can
   be implemented with much less native EM surgery than a full hierarchical
   model.
3. Only if the diagnostic and static regularizer work, consider a full native
   variational/hierarchical tile model.

The feature is computationally modest in a static form: one extra per-unit
orientation array, one tile-count pass, and O(number of EM units) score offsets.
A fully dynamic tile-coupled EM is much harder because the native solver's
equivalence-class compression assumes independent rows with fixed likelihoods.

## Failure Mode Being Targeted

Current locus EM scores an ambiguous unspliced unit approximately as:

```text
score_ik = log_lik_ik + E[log theta_k] - log L_k
```

where `k` is an RNA component or the locus-level gDNA component. For gDNA, the
component effective length `L_g` can be tens of millions of bases in a mega-locus.
For RNA, component effective lengths are usually kilobases to tens of kilobases.

This creates the FLG2-like pathology:

```text
local truth:      gDNA + and gDNA - reads, approximately symmetric
current EM fit:   FLG2 mRNA(-) explains one orientation
                  antisense/nRNA(+) explains the other orientation
```

Fragment by fragment, each RNA assignment can look plausible. Tile by tile, the
fit is suspicious: two independent RNA components must conspire to reproduce a
symmetric local gDNA pattern.

The desired new signal is the aggregate local question:

```text
Does this tile look like one symmetric DNA source,
or like one/two locally active strand-specific RNA sources?
```

## Current Rigel Integration Points

The relevant current plumbing is:

- `src/rigel/scored_fragments.py` defines `ScoredFragments` and
  `LocusPartition`. `ScoredFragments` currently has `genomic_midpoint`, but it
  does not carry an explicit unit-level genomic orientation.
- `src/rigel/scan.py` receives native `StreamingScorer.finish()` arrays and
  constructs `ScoredFragments`. This is the natural place to add unit-level
  orientation and, eventually, tile IDs.
- `src/rigel/locus_partition.py` scatters per-unit arrays into per-locus
  partitions. It currently scatters `gdna_log_liks`, `is_spliced`, `frag_ids`,
  `frag_class`, and `splice_type`, but not midpoint or orientation.
- `src/rigel/estimator.py` calls the native `_batch_locus_em_partitioned`
  function with a 9-tuple per locus. The native ABI currently has no tile-level
  covariates.
- `src/rigel/native/em_solver.cpp` builds equivalence classes from rows with
  identical candidate sets. The E-step assumes fixed per-row/per-candidate
  likelihoods. Dynamic tile factors that depend on current posteriors would
  break or complicate this compression unless the factor is precomputed and
  folded into `log_liks`/`gdna_log_liks` before equivalence classes are built.

The simplest production shape is therefore a precomputed static tile offset:
compute tile-level evidence before EM, then add a bounded per-unit offset to the
gDNA log likelihood, and optionally to RNA candidate log likelihoods. A full
hierarchical model is possible, but it is a larger native solver project.

## Why A Pure Per-Fragment Factor Is Not Enough

Let `o_i` be the observed orientation of fragment `i`, with values `+` or `-`.
Let `q_k(o_i)` be the orientation likelihood under component `k`. A simple
per-fragment mixture objective is:

```text
L_ind = product_i sum_k theta_k b_ik q_k(o_i)
```

where `b_ik` contains non-orientation likelihood terms. This is basically the
current architecture: the EM compares independent rows, each row can place
responsibility on whichever local component best explains that row.

If a tile has 500 plus-orientation fragments and 500 minus-orientation fragments,
and if both a plus-compatible RNA component and a minus-compatible RNA component
exist, the independent mixture can fit the tile by assigning the plus group to
one RNA state and the minus group to the other. There is no penalty for the fact
that the two RNA abundances had to align locally. The per-fragment likelihood is
blind to that coincidence.

The reviewer suggestion becomes powerful only when the tile is treated as a
single aggregated observation with a shared local latent source structure.

## Tile-Level Generative Model

For tile `t`, define:

```text
n_t+ = number of plus-orientation ambiguous unspliced fragments
n_t- = number of minus-orientation ambiguous unspliced fragments
n_t  = n_t+ + n_t-
```

Define coarse orientation/source classes:

```text
G   = gDNA, symmetric orientation
R+  = plus-oriented RNA source
R-  = minus-oriented RNA source
N+  = plus-oriented nRNA source
N-  = minus-oriented nRNA source
```

The exact class names can be adjusted to Rigel's strand conventions. The key is
that gDNA is one symmetric local source, while RNA explanations are oriented and
may require multiple local sources.

For a class `h`, let:

```text
q_h = P(o = + | h)
```

Typical values are:

```text
q_G  = 0.5
q_R+ = 1 - s_err
q_R- = s_err
q_N+ = 1 - s_err
q_N- = s_err
```

where `s_err` is the strand error rate inferred from the library strand model.

Within a tile, local class proportions are:

```text
phi_t = (phi_tG, phi_tR+, phi_tR-, phi_tN+, phi_tN-)
phi_t ~ Dirichlet(alpha)
```

Given `phi_t`, the orientation count likelihood is:

```text
p(n_t+, n_t- | phi_t)
  = choose(n_t, n_t+) *
    (sum_h phi_th q_h) ^ n_t+ *
    (sum_h phi_th (1 - q_h)) ^ n_t-
```

The aggregated marginal likelihood is obtained by integrating out `phi_t`:

```text
p(n_t+, n_t- | alpha, q)
  = choose(n_t, n_t+) *
    integral [
      (sum_h phi_h q_h) ^ n_t+ *
      (sum_h phi_h (1 - q_h)) ^ n_t-
    ] Dirichlet(phi | alpha) d phi
```

This integration is where the useful penalty appears. A model with more local
source proportions can fit more patterns, but it pays for that flexibility after
marginalization.

## Exact Expansion Of The Marginal Likelihood

The integral above can be expanded as a finite sum. Let `a_h+` be the number of
plus observations allocated to class `h`, and `a_h-` the number of minus
observations allocated to class `h`. Define:

```text
sum_h a_h+ = n_t+
sum_h a_h- = n_t-
a_h = a_h+ + a_h-
alpha_0 = sum_h alpha_h
```

Then:

```text
p(n_t+, n_t- | alpha, q)
  = choose(n_t, n_t+) *
    sum_{a+} sum_{a-}
      [ n_t+! / product_h a_h+! ]
      [ n_t-! / product_h a_h-! ]
      [ B(alpha + a) / B(alpha) ]
      product_h q_h ^ a_h+ (1 - q_h) ^ a_h-
```

where `B(.)` is the multivariate beta function:

```text
B(alpha) = product_h Gamma(alpha_h) / Gamma(alpha_0)
```

This exact expression is useful for derivation and for small class sets. In
production we probably do not need the full five-class exact sum; closed-form
special cases or a Laplace approximation are enough.

## Core Bayes Factor: One Symmetric Source vs Two RNA Sources

The minimal hypothesis test compares:

```text
H_G:    one symmetric gDNA source, q_G = 0.5
H_R2:   two oriented RNA sources, q_+ = 1 - s_err, q_- = s_err
```

For `H_G`:

```text
p(n+, n- | H_G) = choose(n, n+) * 0.5 ^ n
```

For a perfect-strand two-RNA model with `s_err = 0` and
`phi = (phi_+, phi_-) ~ Dirichlet(a, a)`:

```text
p(n+, n- | H_R2)
  = choose(n, n+) * B(a + n+, a + n-) / B(a, a)
```

With the uniform prior `a = 1`:

```text
p(n+, n- | H_R2) = 1 / (n + 1)
```

For a balanced tile with `n+ = n- = n / 2`:

```text
p(n+, n- | H_G) approx sqrt(2 / (pi n))
p(n+, n- | H_R2) = 1 / (n + 1)
BF_G_R2 approx sqrt(2 / pi) * sqrt(n)
log BF_G_R2 approx 0.5 log n + constant
```

Examples:

```text
n = 100 balanced:   p_G approx 0.080, p_R2 approx 0.0099, log BF approx 2.1
n = 1000 balanced:  p_G approx 0.025, p_R2 approx 0.0010, log BF approx 3.2
```

This is not an enormous hammer, but it is real Bayesian evidence. With a sparse
RNA prior `a < 1`, the penalty against fortuitous 50/50 two-RNA explanations is
stronger. With a diffuse large `a`, the penalty weakens.

This derivation also explains why the current EM fails. If `phi_+` and `phi_-`
are optimized without integration or prior penalty, the two-RNA model can set
`phi_+ = n+ / n` and `phi_- = n- / n`. For a balanced tile, its maximum
orientation likelihood equals the gDNA orientation likelihood. The distinction
comes from the marginal likelihood, not from the maximum likelihood.

## BIC / Occam Interpretation

At a balanced tile, the best two-RNA orientation model and the one-gDNA model
can have nearly the same maximum likelihood. The difference is parameter count:

```text
H_G:  no local orientation mixture parameter
H_R2: one local mixture parameter phi_+
```

The integrated likelihood approximately subtracts:

```text
(d / 2) log n
```

for `d` fitted local parameters. This is the Occam factor. It is exactly the
kind of penalty Rigel is missing: the current EM gets to use two local RNA
states without paying for the fact that their local proportions were tuned to a
symmetric DNA-like pattern.

## Incorporating Existing Non-Orientation Evidence

A pure orientation-count test is too blunt. It must be conditioned on RNA anchor
evidence and ordinary fragment likelihoods.

Important auxiliary evidence:

- Spliced annotated reads in or near the tile. Strong same-orientation spliced
  support should make the RNA hypothesis less surprising.
- Implicit-splice support and exon junction context. These are RNA-specific and
  should reduce a gDNA guard.
- Intronic/nascent coverage beyond the exon. True nRNA should usually have
  broader local support than a pileup restricted to capture-rich exonic bait.
- Capture/exposure density. A symmetric unspliced tile in a high-gDNA-exposure
  region should be easier to call gDNA than the same count pattern in a region
  with little gDNA opportunity.
- Fragment class. The target population is unspliced, genic, ambiguous units.
  Spliced transcript-only units should not be penalized by this factor.

A practical tile Bayes factor should therefore be shaped as:

```text
log BF_t = log p(orientation counts | H_G)
         - log p(orientation counts | H_R)
         + log prior_odds_t(H_G / H_R)
```

where `prior_odds_t` is allowed to depend on anchored RNA evidence, regional
gDNA exposure, and calibration estimates.

## Proposed Static Regularizer

The safest first production implementation is a bounded static regularizer,
not a fully dynamic hierarchical EM.

For each tile `t`:

1. Count eligible ambiguous unspliced EM units by genomic orientation.
2. Compute `log BF_t` for symmetric gDNA vs oriented RNA split.
3. Gate the factor by minimum count, RNA anchor support, and evidence quality.
4. Convert the tile factor to a per-unit score offset.

Because the evidence is a tile-level factor, it must not be naively added once
per fragment at full strength. If a tile has `n_t` eligible units, distribute the
tile factor:

```text
delta_t = eta * clip(log BF_t, -B_max, B_max) / max(n_t, 1)
```

Then, for eligible units in that tile:

```text
gdna_log_lik_i <- gdna_log_lik_i + delta_t
```

This means that if all eligible units in the tile choose gDNA, the total added
support is approximately the tile-level log Bayes factor, not `n_t` times that
factor.

This gDNA-only positive guard is intentionally conservative. It does not require
rewriting RNA candidate likelihoods and it directly targets the observed false
negative for gDNA. A later version can add class-specific RNA penalties or
offsets if diagnostics show that gDNA-only offsets are insufficient.

Recommended default gates for a first experiment:

```text
eligible units:       unspliced, finite gDNA likelihood, genic ambiguous
minimum tile count:   30 or 50
tile size:            1 kb or 2 kb, plus optional 5 kb sensitivity run
positive cap B_max:   5 to 10 total log units per tile
negative cap:         0 initially, so the factor cannot hurt gDNA recall
eta:                  0.25, 0.5, 1.0 sweep
RNA anchor guard:     suppress or shrink if strong spliced RNA support exists
```

The positive-only first version is deliberately asymmetric. The current product
problem is gDNA being leaked into RNA. We should first test whether local
symmetric evidence can rescue gDNA without increasing RNA->gDNA too much.

## More Principled Variational Form

The static regularizer is an approximation. The cleaner generative model adds
tile-local class proportions to the EM.

Let `c(k)` map each component to a coarse class such as `G`, `R+`, `R-`, `N+`,
or `N-`. Let `phi_t` be the local tile class proportions. A variational E-step
would score component `k` for unit `i` in tile `t` as:

```text
score_ik = base_log_lik_ik
         + E[log theta_k]
         - log L_k
         + E[log phi_t,c(k)]
```

The tile-local posterior would update as:

```text
alpha_t,h' = alpha_h + sum_{i in tile t} sum_{k: c(k)=h} r_ik
E[log phi_t,h] = digamma(alpha_t,h') - digamma(sum_h alpha_t,h')
```

If the tile model also has a latent hypothesis `H_t` such as `gDNA-symmetric` vs
`RNA-oriented`, then:

```text
p(H_t | data) proportional p(H_t) p(data_t | H_t)
```

and `E[log phi_t,h]` is averaged over those hypotheses.

This is theoretically cleaner because the tile factor updates with the current
posterior assignments. It is also harder in Rigel because the native EM solver
builds equivalence classes with fixed row likelihoods. If `E[log phi_t,h]`
changes every iteration and differs by tile, equivalence classes either need to
include tile-specific offsets in their keys or the E-step must operate on rows
instead of compressed equivalence classes. That could substantially increase
runtime in mega-loci.

## Avoiding Double Counting Strand Likelihood

Current `log_liks` already include strand likelihood terms. A full generative
implementation should eventually separate:

```text
base_log_lik_ik = non-strand terms
orientation_log_lik_ik = log q_c(k)(o_i)
```

The first static implementation does not need this separation if the tile factor
is treated as an aggregate prior/regularizer rather than another copy of the
per-fragment orientation likelihood. But the documentation and code should be
clear: the offset is not a second strand model. It is a tile-level model-choice
adjustment derived from the aggregate orientation counts.

Long-term, the clean architecture is:

```text
fragment likelihood = non-strand alignment/FL/overlap likelihood
                    * per-fragment strand likelihood
                    * tile-level source-structure factor
```

The tile factor should be calibrated so that the total contribution per tile is
bounded and interpretable.

## Data Requirements

Production implementation needs these per EM unit fields:

```text
unit_ref_id          int32, or enough locus metadata to recover reference
unit_midpoint        int64, already present globally as genomic_midpoint
unit_orientation     int8  (+1, -1, 0 unknown)
unit_is_spliced      uint8, already present
unit_gdna_eligible   finite gdna_log_lik and unspliced
```

For a gDNA-only static offset, these are sufficient.

For RNA class-specific offsets, we also need per candidate or per component:

```text
component_class      G/R+/R-/N+/N-
candidate_orientation_class for each candidate entry, or a map from transcript
                     strand and nRNA status to coarse class
```

Current arrays probably contain partial information:

- `locus_count_cols` and `count_cols` encode sense/antisense relative to a
  transcript candidate.
- `index.t_to_strand_arr` gives transcript genomic strand.
- Native resolver/scanner already tracks `exon_strand`, `sj_strand`, and
  `ambig_strand` for buffered fragments.

However, the most robust production design is to explicitly propagate a
unit-level genomic orientation from the native scanner/scorer into
`ScoredFragments`. Relying on candidate-relative count columns risks confusing
observed genomic orientation with compatibility relative to whichever candidate
was selected as the locus representative.

## Tile Definition

Start with fixed genomic bins within each reference:

```text
tile_id = (ref_id, floor(midpoint / bin_size))
```

Candidate bin sizes:

```text
1 kb:   sensitive to local capture/exon pileups, but noisier
2 kb:   likely default starting point
5 kb:   smoother; useful for FLG2-like broad hotspots
10 kb:  diagnostic only, may blur distinct expression features
```

Use a minimum eligible count threshold so low-count tiles do not inject noisy
priors. For fragments without valid midpoint or orientation, do not apply the
factor.

Potential later refinement: define tiles by exonic/capture intervals rather than
fixed bins. That might align better with hybrid-capture biology, but fixed bins
are easier to validate and less likely to accidentally encode transcript/gene
assumptions.

## Implementation Plan

### Phase 0: Offline Diagnostic Bayes Factor

Goal: prove the local orientation factor has signal before touching production
quantification.

Tasks:

1. Add a diagnostic script under `scripts/debug/` that scans an annotated BAM or
   saved EM-unit table and builds fixed genomic tiles.
2. Count orientation by truth source and by Rigel predicted pool for VCaP.
3. Compute simple closed-form Bayes factors:
   - `H_G` symmetric gDNA.
   - `H_R1` one oriented RNA source.
   - `H_R2` two oriented RNA sources with Dirichlet concentration `a`.
4. Report tile-level ROC-like summaries:
   - log BF vs true gDNA fraction.
   - log BF in gDNA->RNA false-positive tiles.
   - FLG2 and top hotspot case studies.
5. Sweep tile size, minimum count, and `a`.

Estimated difficulty: low to moderate. This can be done without native code if
orientation can be recovered from BAM flags and library conventions in the
diagnostic script. For production parity, native unit orientation still needs to
be added later.

Exit criteria:

```text
log BF_t should be strongly positive in known true-gDNA symmetric hotspots,
including FLG2, and should not be broadly positive in tiles with strong spliced
RNA truth/support.
```

### Phase 1: Static gDNA Offset Regularizer

Goal: test whether tile-level evidence reduces VCaP gDNA->RNA leakage in full
quantification.

Tasks:

1. Add config fields under calibration or EM config:
   - `local_orientation_factor_enabled`
   - `local_orientation_tile_bp`
   - `local_orientation_min_count`
   - `local_orientation_dirichlet_alpha`
   - `local_orientation_strength_eta`
   - `local_orientation_logbf_cap`
   - `local_orientation_anchor_shrinkage`
2. Propagate `unit_orientation` from native scanner/scorer into
   `ScoredFragments`.
3. Compute tile counts before `partition_and_free()` while global midpoint and
   orientation arrays are available.
4. Compute `delta_t = eta * clipped_logBF_t / n_t`.
5. Add `delta_t` to `em_data.gdna_log_liks` for eligible units only.
6. Record summary diagnostics in `summary.json`:
   - number of tiles evaluated.
   - number of tiles gated out.
   - distribution of log BF and applied delta.
   - total positive/negative score mass applied.
7. Run focused tests and VCaP quantification.

Estimated difficulty: moderate. It touches native scan output, Python data
containers, partition-time data flow, config/CLI plumbing, and tests. It does
not require modifying the native EM solver if the final product is just a static
offset to `gdna_log_liks`.

Primary risk: if the offset is too strong or insufficiently gated, real
unspliced RNA in bidirectional/antisense regions can be pushed into gDNA.

### Phase 2: Candidate-Class Offsets

Goal: distinguish "boost gDNA" from "penalize unsupported RNA split".

Tasks:

1. Map every component/candidate to a coarse class: annotated mRNA plus/minus,
   synthetic/annotated nRNA plus/minus, and gDNA.
2. Compute tile-specific class offsets from `p(H_t | data_t)`.
3. Apply offsets to candidate `log_liks` and gDNA `gdna_log_liks` before native
   EM.
4. Keep total tile influence bounded by distributing offsets across eligible
   units/candidates.

Estimated difficulty: moderate to high. Mutating candidate `log_liks` before
partition is feasible, but getting class assignment exactly right for synthetic
nRNA states and candidate-relative strand columns requires care.

Expected benefit: better specificity than a gDNA-only boost. For example, a
tile with symmetric orientations but strong spliced FLG2 mRNA support might keep
the supported mRNA class while penalizing an unsupported antisense/nRNA split.

### Phase 3: Native Variational Tile Model

Goal: implement the full hierarchical local source model.

Tasks:

1. Extend the native EM ABI to accept per-unit tile IDs and per-component class
   IDs.
2. Maintain tile-class variational parameters during EM.
3. Update tile factors between E-steps.
4. Rework equivalence-class grouping so tile offsets are represented correctly.
5. Add profiling to measure EC fragmentation and runtime costs.

Estimated difficulty: high. This is a native EM redesign, not a calibration
patch. It is probably a 2-4 week effort including tests, benchmarks, and
performance cleanup.

Expected benefit: most principled model, but only justified if Phases 0-2 show
clear empirical value.

## Computational Cost

### Static Diagnostic / Regularizer

For `N` EM units and `C` candidate entries:

```text
tile counting:       O(N)
gDNA offset apply:   O(N)
RNA offset apply:    O(C), only if candidate-class offsets are enabled
memory:              unit_orientation int8[N], optional tile_id int32[N]
```

For millions of EM units, the memory overhead is small compared with existing
CSR candidate arrays. A rough budget is:

```text
orientation: 1 byte per EM unit
tile_id:     4 bytes per EM unit, optional if computed on the fly
delta:       no persistent array needed if applied in place
```

Expected runtime overhead for a gDNA-only static regularizer should be low,
probably below 5-10% of quantification runtime, and much less than the native EM
runtime in mega-loci.

### Exact Marginal Likelihood Computation

The minimal two-source formulas are closed form and O(1) per tile.

For more than two classes, exact summation over all allocations can become too
expensive for high-count tiles. Use one of:

```text
closed-form special cases for H_G, H_R1, H_R2
Laplace approximation for larger class sets
small fixed quadrature over phi for two-source mixtures with strand error
```

There is no need to enumerate fragment allocations in production.

### Dynamic Native Model

A full variational model can be substantially more expensive because tile
offsets change every iteration. If equivalence classes must be split by tile and
class-offset state, the number of equivalence classes can grow. In mega-loci,
that could increase both memory and E-step time. This is why a static offset is
the right first implementation.

## Expected Benefits

The feature directly targets the observed FLG2-style failure:

```text
many local unspliced ambiguous fragments
approximately 50/50 genomic orientation
little true RNA source by flowcell truth
current EM splits the pileup into mRNA plus nRNA
```

In such tiles, the gDNA hypothesis gets a tile-level evidence boost because one
symmetric source explains the aggregate pattern without tuning two RNA sources.

Expected improvements if the hypothesis is right:

- lower gDNA->mRNA in symmetric exonic hotspots.
- lower gDNA->nRNA where antisense/nRNA is used as the second half of a gDNA
  pileup.
- better behavior in mega-loci where scalar gDNA denominator changes have not
  solved the problem.
- clearer diagnostics: false RNA calls can be stratified by whether local
  orientation evidence favored gDNA.

## Risks And Failure Modes

1. Real bidirectional RNA can look symmetric.
   A promoter, antisense locus, or overlapping genes can produce real plus and
   minus RNA. The guard must shrink when spliced or otherwise RNA-specific
   anchors support both orientations.

2. nRNA can be unspliced and local.
   A true nascent signal may lack splice junctions in a short tile. The factor
   should not declare all unspliced symmetry to be gDNA.

3. Strand model misspecification can bias the test.
   If the library strand model has local or protocol-specific deviations, the
   assumed `q_R+`, `q_R-`, and `q_G` values can be wrong.

4. Tile size matters.
   Small tiles are noisy. Large tiles can merge separate biological processes.
   Validation must sweep tile size.

5. Double counting is easy.
   Because strand likelihood already exists per fragment, the tile factor must
   be documented and implemented as aggregate model evidence, with total
   influence capped per tile.

6. The Bayes factor may be too weak alone.
   The Occam penalty grows roughly like `0.5 log n` for the simple balanced
   two-RNA case. It may need to be combined with RNA anchor priors and regional
   gDNA exposure to overcome large RNA effective-length advantages.

## Validation Plan

### Synthetic Tests

Create synthetic loci with known truth:

1. symmetric gDNA-only exon tile with plus/minus orientations.
2. one-orientation mRNA tile with strong spliced support.
3. true bidirectional RNA tile with spliced support on both strands.
4. true nRNA tile with broad intronic support.
5. FLG2-like ambiguous tile with large gDNA denominator and two RNA-compatible
   states.

Expected outcomes:

```text
case 1: gDNA boost
case 2: no gDNA boost or very small boost
case 3: no false gDNA takeover
case 4: no false gDNA takeover when nRNA support extends beyond exon-only pileup
case 5: reduced split into mRNA+nRNA
```

### VCaP Real Mixture

Use flowcell truth as before.

Primary metrics:

```text
gDNA -> RNA total
gDNA -> mRNA
gDNA -> nRNA
RNA -> gDNA
gDNA recall
RNA recall
```

Required stratifications:

```text
by splice type: unspliced vs spliced
by fragment class: ambig_same_strand, ambig_opp_strand, multimapper
by tile log BF decile
by top hotspot windows, including FLG2
by locus size and gDNA effective length
```

Success criterion for Phase 1 should be strict:

```text
Reduce gDNA->RNA materially versus the current best kappa-units/exposure
baseline, without increasing RNA->gDNA enough to erase the gain.
```

A useful first target would be at least a 20-30% reduction in gDNA->RNA false
calls in high-log-BF tiles, with less than a 5-10% relative increase in RNA->gDNA
overall. If the effect only moves errors from mRNA to nRNA or increases RNA loss,
the model is not ready.

## Configuration Defaults For First Quant Experiment

Suggested experimental config:

```yaml
calibration:
  local_orientation_factor_enabled: true
  local_orientation_tile_bp: 2000
  local_orientation_min_count: 50
  local_orientation_dirichlet_alpha: 0.5
  local_orientation_strength_eta: 0.5
  local_orientation_logbf_cap: 8.0
  local_orientation_positive_only: true
  local_orientation_anchor_shrinkage: true
```

Suggested sweep:

```text
tile_bp:      1000, 2000, 5000
min_count:    30, 50, 100
alpha:         0.1, 0.5, 1.0
eta:           0.25, 0.5, 1.0
cap:           5, 8, 10
```

Do not expose this as a default production feature until it wins on VCaP and on
synthetic counterexamples.

## Implementation Difficulty Summary

```text
Phase 0 diagnostic:                  1-2 days
Phase 1 static gDNA offset:           3-5 days
Phase 2 candidate-class offsets:      1-2 weeks
Phase 3 full native variational EM:   2-4 weeks
```

The critical engineering unknown is unit-level orientation propagation. If the
native scanner can emit a clean `unit_orientation` alongside `genomic_midpoint`,
Phase 1 is straightforward. If orientation must be reconstructed indirectly from
candidate-relative count columns, the implementation is riskier and should be
kept diagnostic-only until validated.

## Bottom Line

The concept has a sound theoretical basis if it is implemented as a local
aggregated marginal likelihood or hierarchical prior. It is not a hack in that
form: it is a Bayesian model-selection penalty for explaining one symmetric DNA
source with multiple locally tuned RNA sources.

It will not solve every gDNA->RNA error. It specifically addresses the local
symmetric pileup failure mode seen in FLG2 and related exon-rich hotspots. The
right next step is diagnostic: measure whether the tile Bayes factor is high in
the false-RNA gDNA hotspots and low in genuine RNA regions. If that separation
exists, a bounded static gDNA offset is the lowest-risk implementation path.