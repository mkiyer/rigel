# Calibration System: Implementation vs Theory Survey

This document surveys the current calibration implementation against the
theoretical model defined in `THEORETICAL_MODEL.md`,
`THEORETICAL_EM_OBJECTIVE.md`, and `FIRST_PRINCIPLES_IMPLEMENTATION_PLAN.md`.
It identifies alignment, divergence, and missing pieces, then proposes
concrete integration steps.

## 1. Executive Summary

The calibration subsystem has achieved Phases 0–4 of the implementation
plan: region partition, region counting, regional evidence model, two-state
purity EM, and partial empirical-Bayes nuisance learning.  Stress testing
shows dramatic improvement over the histogram-based predecessor (gDNA F1:
0.946 → 0.987; SS=0.50 gDNA F1: 0.737 → 0.964).

However, the calibration system currently operates as a **standalone
diagnostic** — it is not wired into the production pipeline.  The main locus
EM still uses the old strand-only gDNA density estimation and a separate FL
model trained from unique mapper strand weighting.  Bridging this gap is the
single most important remaining integration task.

The theoretical model calls for `T + N + 2` per-locus components (two
strand-specific gDNA components) with a Beta prior on the strand fraction
φ_ℓ and empirical-Bayes κ_sym.  The current native EM uses `T + N + 1`
(single collapsed gDNA) with log(0.5) strand likelihood baked into scoring.
This collapsed form is the correct marginal of the two-strand model under
symmetry — but it does not support per-locus strand-fraction estimation or
overdispersion handling.

## 2. What Exists: Implementation Inventory

### 2.1 Region Partition (Workstream A — Complete)

**Theory target**: Non-overlapping genomic region partition with four
boolean annotation flags: `exon_pos`, `exon_neg`, `tx_pos`, `tx_neg`.

**Implementation**:
- `index.py::build_region_table()` — Boundary-sweep algorithm producing
  atomic half-open bins, merged when adjacent bins share identical flags.
  Four stored booleans match the theory exactly.
- Persisted as `regions.feather`, loaded into cgranges interval tree at
  index load time.
- Region classes (intergenic, intronic, exonic, ambiguous) are derived
  from the four flags, not stored separately.

**Assessment**: Fully aligned with theory.  No changes needed.

### 2.2 Fragment-to-Region Counting (Workstream B — Complete)

**Theory target**: Per-region evidence table with strand-resolved counts
and fragment-length observations.

**Implementation**:
- `region_evidence.py::count_region_evidence()` — Second BAM pass (via
  pysam) producing:
  - `region_counts`: (n_regions × 4) with `n_unspliced_pos`,
    `n_unspliced_neg`, `n_spliced_pos`, `n_spliced_neg`
  - `fl_table`: Fragment-length observations for unspliced single-region
    fragments (columns: `region_id`, `frag_len`)
- Fractional counting when fragments span multiple regions.
- Chimeric fragments excluded.

**Divergence from theory**: The implementation plan calls for evidence
extraction from the existing fragment buffer ("post-scan aggregation over
the existing fragment buffer ... not a second BAM pass").  The current
`count_region_evidence()` performs an independent second BAM pass via
pysam, which is:
1. Slower (full BAM re-read)
2. Not integrated with the C++ scan pass
3. May produce slightly different fragment assembly than the native scanner

**Assessment**: Functionally correct but architecturally misaligned.
Integration should either (a) port region counting into the C++ scanner
as a co-task of the existing single pass, or (b) operate on the
`FragmentBuffer` output rather than re-reading BAM.

### 2.3 Regional Evidence Model (Workstream C — Complete)

**Theory target**: Three convergent evidence channels — splice, density,
strand asymmetry.

**Implementation** (all in `calibration.py`):

| Channel | Model | Formula |
|---------|-------|---------|
| Splice | Hard constraint | `n_spliced > 0 → γ = 0` |
| Density | Gaussian on log(n/L + ε) | `LLR = log N(x|μ_G,σ²_G) − log N(x|μ_R,σ²_R)` |
| Strand | Binomial | `LLR = k·log(0.5/SS) + (n−k)·log(0.5/(1−SS))` |
| Frag length | Shape-normalized histogram | Shared Dirichlet(1/M) prior on proportions |

**Assessment**: All three channels are principled and parameter-free (no
hardcoded thresholds controlling behavior).  The Binomial strand model
naturally degrades to zero when SS = 0.5.  The Gaussian density model
avoids the histogram bin/smoothing artifacts of the predecessor.  The
shape-normalized FL model eliminates the systematic weight-dependent bias.

### 2.4 Two-State Purity Model (Workstream D — Complete)

**Theory target**: Per-region posterior weights
w_r = P(gDNA-dominant | data).

**Implementation**:
- `calibrate_gdna()` — Full EM with:
  - Two-phase seed initialization (expressed from spliced; gDNA from
    density eCDF + optional strand filter)
  - Bayesian π_soft (mixing proportion conditioned on unspliced regions
    only, not globally diluted by spliced regions)
  - Convergence on |Δπ_soft| < 10⁻⁴
  - Per-iteration FL model rebuilding
- Outputs `GDNACalibration` dataclass with:
  - `region_posteriors`: γ_r ∈ [0,1] per region
  - `gdna_density_global`, `gdna_density_per_ref`
  - `kappa_sym` (post-hoc Beta-Binomial MLE)
  - `gdna_fl_model` (posterior-weighted)
  - `mixing_proportion`, `expressed_density`

**Assessment**: Matches the "Workstream D: First-Pass Purity Model"
target.  The soft posterior weights are exactly the theory's w_r.

### 2.5 Empirical-Bayes Nuisance Learning (Workstream E — Partial)

**Theory target**: Joint weighted estimation of gDNA FL distribution and
κ_sym from purity-weighted regional evidence, with alternating refinement.

**Implementation**:
- **gDNA FL model**: `build_gdna_fl_model()` builds a posterior-weighted
  histogram from `fl_table` and region posteriors.  Updated every EM
  iteration.  ✅
- **κ_sym estimation**: `estimate_kappa_sym()` — Beta-Binomial MLE via
  golden-section maximization on γ-weighted regions.  Computed post-hoc
  after convergence.  ⚠️ Not fed back into the EM.
- **Self-consistent alternation**: The FL model alternates (rebuilt each
  iteration), but κ_sym does not.  The theory calls for joint iterative
  refinement of both.

**Assessment**: FL estimation is self-consistent within the EM.  κ_sym is
a one-shot post-hoc estimate — correct in value but not self-consistent
(not used to refine strand evidence during EM iterations).

### 2.6 Main Locus EM Integration (Workstream F — Not Started)

**Theory target**: Cleaner `T + N + 2` parameterization with calibrated
nuisance parameters feeding into the locus EM.

**Current state**: The calibration system is entirely disconnected from
the production pipeline:

| Calibration Output | How Pipeline Gets Same Info |
|---|---|
| `gdna_density_global` | `compute_gdna_density_hybrid()` in locus.py (strand correction) |
| `gdna_density_per_ref` | `_compute_ref_gdna_densities()` in locus.py |
| `gdna_fl_model` | `frag_length_models.gdna_model` via strand-weighted mixing in `frag_length_model.py` |
| `region_posteriors` | Not used; pipeline routes fragments directly to loci |
| `kappa_sym` | Not used; main EM has no strand symmetry prior |

The main EM uses a single collapsed gDNA component with log(0.5) baked
into per-fragment scoring — no per-locus φ_ℓ, no κ_sym prior.

## 3. Theory–Implementation Alignment Matrix

| Theoretical Concept | Status | Location |
|---|---|---|
| Region partition with 4 boolean flags | ✅ Complete | `index.py::build_region_table` |
| Fragment-to-region counting | ✅ Functional (not integrated) | `region_evidence.py` |
| Splice hard constraint | ✅ Complete | `calibration.py::_e_step` |
| Gaussian density channel | ✅ Complete | `calibration.py::_compute_density_llr_gaussian` |
| Binomial strand channel | ✅ Complete | `calibration.py::_compute_strand_llr_binomial` |
| Shape-normalized FL channel | ✅ Complete | `calibration.py::_compute_fl_llr` |
| Two-state purity EM | ✅ Complete | `calibration.py::calibrate_gdna` |
| Soft calibration weights w_r | ✅ Complete | `GDNACalibration.region_posteriors` |
| gDNA FL from purity weights | ✅ Complete | `calibration.py::build_gdna_fl_model` |
| κ_sym estimation (Beta-Binomial MLE) | ✅ Post-hoc | `calibration.py::estimate_kappa_sym` |
| β-Binomial strand LLR in calibration | ❌ Missing | Use Binomial; no overdispersion |
| Self-consistent EB alternation (κ_sym ↔ weights) | ❌ Missing | κ_sym not fed back |
| Two gDNA components per locus (g_pos, g_neg) | ❌ Missing | Single collapsed gDNA |
| Beta(κ/2, κ/2) prior on φ_ℓ in locus EM | ❌ Missing | log(0.5) in scoring only |
| Calibrated κ_sym → locus EM | ❌ Missing | Not connected |
| Calibrated gDNA FL → locus scoring | ❌ Missing | Pipeline uses separate FL model |
| Calibrated gDNA density → locus prior | ❌ Missing | Pipeline uses strand-only EB |
| Region counting from fragment buffer | ❌ Architectural gap | Uses 2nd BAM pass |
| Asymmetric gDNA constraint (sense-heavy suspect) | ❌ Not addressed | Theory §9.2 |

## 4. Key Divergences and Their Significance

### 4.1 Binomial vs Beta-Binomial Strand Model

**Current**: The calibration EM uses a simple Binomial strand LLR.  For
region r with n unspliced fragments and k sense-strand:

$$
\text{LLR}_r = k \log\frac{0.5}{\text{SS}} + (n-k)\log\frac{0.5}{1-\text{SS}}
$$

**Theory**: The Beta-Binomial marginal model accounts for overdispersion
in strand fractions across regions:

$$
K_r \mid N_r, \kappa_\text{sym} \sim \text{BetaBinomial}(N_r, \kappa_\text{sym}/2, \kappa_\text{sym}/2)
$$

**Significance**: The Binomial model assumes all gDNA regions share an
exact 50/50 strand split, with deviations due only to sampling noise.  In
reality, local genomic features (replication timing, chromatin state,
capture bias) can create real regional strand imbalance.  When κ_sym is
moderate (say 20–50), the Beta-Binomial marginalization properly relaxes
the strand signal for high-count regions where the Binomial would
over-confidently classify them based on small strand deviations.

**Risk of not fixing**: At high sequencing depth, the Binomial LLR grows
without bound for any region with even slight true strand asymmetry.  This
could create false-positive RNA classifications in deeply-sequenced gDNA
regions with genuine strand imbalance.

**Solution**: Replace the Binomial strand LLR with a Beta-Binomial LLR
parameterized by κ_sym.  The infrastructure already exists:
`estimate_kappa_sym()` computes κ_sym via MLE.  The fix requires:

1. Add `_compute_strand_llr_betabinom()`:
   $$
   \text{LLR}_r = \log\text{BetaBinom}(k|n,\kappa/2,\kappa/2) - \log\text{Binom}(k|n,\text{SS})
   $$
2. Initialize κ_sym from gDNA seed regions before the EM loop.
3. Optionally update κ_sym each iteration (self-consistent EB).

When κ_sym → ∞, this degenerates to the current Binomial model.  When
κ_sym is finite, it naturally tempers the strand LLR for regions whose
strand balance is within the expected overdispersion range.

### 4.2 Single gDNA Component vs Two-Strand gDNA Pair

**Current native EM**: One collapsed gDNA component per locus.  Strand
handled by adding log(0.5) per fragment at scoring time.  No per-locus
strand fraction φ_ℓ.  No strand symmetry prior.

**Theory**: Two explicit components g_pos and g_neg per locus with
deterministic strand likelihoods (p(s=pos|g_pos) = 1, p(s=neg|g_neg) = 1)
and a symmetric Beta prior on φ_ℓ = θ_{g_pos}/(θ_{g_pos} + θ_{g_neg}).

**Mathematical equivalence**: The theory document (§8, §9) establishes
that the collapsed single-gDNA model with p(s|g) = 0.5 is the **correct
marginal** of the two-strand model under a symmetric prior.  The two
formulations produce identical E-step posteriors when φ_ℓ is exactly 0.5.

**Where they differ**: The two-strand model allows φ_ℓ ≠ 0.5 at each
locus, regularized by Beta(κ/2, κ/2).  This matters when:
- A locus has genuinely asymmetric gDNA (e.g., from replication strand
  bias or capture asymmetry)
- The M-step can adjust φ_ℓ per locus to better fit the data while
  remaining penalized toward symmetry

The collapsed model forces φ_ℓ = 0.5 everywhere — it cannot adapt to
locus-level strand asymmetry.

**Practical significance**: For most loci with moderate gDNA, the
collapsed model is adequate.  The two-strand model primarily helps at:
- Loci with many gDNA fragments and genuine strand imbalance
- Loci where gDNA and nRNA overlap, creating strand-based identifiability
  issues (the sense-heavy gDNA excess vs. nRNA leakage problem from
  theory §18.2 and EM objective §9.2)

**Implementation cost**: Moderate.  Requires:
- Changing the native EM component layout from `T + N + 1` to `T + N + 2`
- Splitting gDNA candidates by observed strand
- Adding the Beta(κ/2, κ/2) prior term to the M-step
- Summing g_pos + g_neg at reporting time

### 4.3 Calibration Not Connected to Pipeline

This is the most significant practical gap.  The calibration system
produces high-quality outputs that the pipeline ignores:

**gDNA FL model**: The pipeline currently builds the gDNA FL model from
strand-weighted mixing in `frag_length_model.py::mix_models()`.  This
uses a heuristic: weight fragments by their probability of being
antisense (proxy for gDNA).  The calibration approach is more principled:
weight FL observations by region-level posterior P(gDNA-dominant).

**gDNA density priors**: The pipeline uses
`compute_gdna_density_hybrid()` — a strand-correction formula that
works only when 2·SS − 1 > 0.2.  The calibration provides a
posterior-weighted density estimate that works at any SS.

**κ_sym**: Not used anywhere in the main EM.  The native solver has
no strand symmetry prior on gDNA.

### 4.4 Second BAM Pass for Region Counting

The theory plan explicitly calls for region counting from the existing
fragment buffer, not a second BAM pass.  The current `count_region_evidence()`
re-reads the entire BAM through pysam.  This is:
- Architecturally wrong (violates single-pass constraint)
- Performance penalty (doubles I/O for large BAMs)
- Potential for inconsistency (different fragment assembly than C++ scanner)

The fix is to either:
1. Port region counting into the C++ `_bam_impl` scanner as a co-task
2. Operate on the `FragmentBuffer` output post-scan
3. Use the `ScoredFragments` data with region lookups

Option (2) is likely most practical for the first integration.

## 5. What Is Missing Entirely

### 5.1 Asymmetric gDNA Strand Constraint

The theory (§18.2, EM objective §9.2) identifies a real identifiability
risk: at loci with both gDNA and nRNA, sense-heavy gDNA excess is
confounded with nRNA leakage.  A purely symmetric Beta prior does not
address this.

**Proposed approach**: After the symmetric M-step update, apply a soft
asymmetric penalty that discounts sense-heavy gDNA excess relative to
antisense-heavy excess.  This can be formulated as a modified Beta prior:

$$
\phi_\ell \sim \text{Beta}(\alpha_{pos}, \alpha_{neg})
$$

where α_pos ≤ α_neg when the gene is on the positive strand (sense-heavy
gDNA is suspect).  The asymmetry can be controlled by a single
interpretable parameter derived from the nRNA evidence at the locus.

Alternatively, the existing nRNA initialization (strand-corrected intronic
evidence in `locus.py::compute_nrna_init()`) already separates nRNA from
gDNA using strand.  If nRNA initialization is accurate, the nRNA component
absorbs sense-heavy unspliced signal, and the residual gDNA should be
approximately symmetric.  In this view, the asymmetric constraint is a
safety net for when nRNA init is inaccurate, not a primary mechanism.

### 5.2 Self-Consistent κ_sym Estimation

Currently κ_sym is computed once post-hoc.  The theory calls for
alternating refinement:

1. Initialize provisional purity weights
2. Estimate gDNA FL distribution → update FL LLR
3. Estimate κ_sym → update strand LLR
4. Recompute purity weights
5. Iterate until stable

The FL channel already does this (rebuilt each iteration).  Adding κ_sym
to the loop requires:
- Computing κ_sym from current γ-weighted regions each iteration
- Using κ_sym in the Beta-Binomial strand LLR (not the Binomial)

This creates a fully self-consistent calibration.

### 5.3 Continuous Purity Model

The theory's long-run target (Workstream D) is a continuous latent purity
variable π_r ∈ [0,1] rather than a hard two-state mixture.  The current
implementation's soft posteriors γ_r are already continuous, but they arise
from a two-state mixture and are bounded by the mixing proportion π.

For the current stress test performance, the two-state model appears
sufficient.  Upgrading to a continuous purity model can be deferred.

## 6. Integration Plan

### Phase I: Connect Calibration to Pipeline (Highest Priority)

**Goal**: Wire calibration outputs into the production pipeline so that
locus EM benefits from calibrated gDNA parameters.

**Step 1: Region Counting from Fragment Buffer**

Replace the second BAM pass with buffer-based counting.  After
`scan_and_buffer()` produces the fragment buffer, extract region
evidence from the buffer's stored fragment coordinates.

Implementation approach:
- Add a `count_regions_from_buffer()` function that iterates the
  `FragmentBuffer`, looks up each fragment's genomic footprint against
  the region cgranges, and accumulates the same four count columns
  plus FL observations.
- This runs between `scan_and_buffer()` and `quant_from_buffer()`.

**Step 2: Run Calibration in Pipeline**

Insert `calibrate_gdna()` into `pipeline.py::quant_from_buffer()`:
```
scan_and_buffer() → count_regions_from_buffer() → calibrate_gdna()
    → [calibrated params] → build_loci() → locus EM
```

**Step 3: Feed Calibrated gDNA FL into Scoring**

Replace the strand-weighted FL mixing in `frag_length_model.py` with the
calibration's posterior-weighted gDNA FL model:
- `FragmentScorer` receives `GDNACalibration.gdna_fl_model` instead of
  `frag_length_models.gdna_model`
- The calibrated model already has the correct shape because it was built
  from purity-weighted regions

**Step 4: Feed Calibrated gDNA Density into EB Priors**

Replace `compute_gdna_density_hybrid()` with calibration's
`gdna_density_global` and `gdna_density_per_ref`:
- `compute_eb_gdna_priors()` uses calibration densities instead of
  strand-only estimation
- This immediately improves gDNA priors at low SS (where the strand
  correction formula breaks down)

**Step 5: Feed κ_sym into Native EM (requires C++ change)**

Pass `GDNACalibration.kappa_sym` through to the native EM solver.  In
the collapsed single-gDNA model, κ_sym can be used to modulate the gDNA
strand penalty:
- Current: fixed log(0.5) per fragment
- Improved: adjust the effective strand log-likelihood based on κ_sym
  and the locus-level observed strand balance

This step is a bridge until the full two-strand model is implemented.

### Phase II: Beta-Binomial Strand Model in Calibration

**Goal**: Replace Binomial strand LLR with Beta-Binomial to handle
overdispersed strand fractions.

**Step 1: Implement `_compute_strand_llr_betabinom()`**

For each region r:
$$
\text{LLR}_r = \log\text{BetaBinom}(k_r | n_r, \kappa/2, \kappa/2) - \log\text{Binom}(k_r | n_r, \text{SS})
$$

Use `_vec_lgamma` (already available) for the Beta-Binomial log-PMF.

**Step 2: Self-consistent κ_sym**

Move `estimate_kappa_sym()` inside the EM loop:
- After each M-step, re-estimate κ_sym from γ-weighted gDNA regions
- Use updated κ_sym in next E-step's Beta-Binomial LLR

**Step 3: Initialize κ_sym from seed regions**

Before the EM loop, compute an initial κ_sym from the gDNA seed
partition (the Phase 2 gDNA seeds from `_seed_initial_partition()`).

### Phase III: Two-Strand gDNA Model in Native EM

**Goal**: Implement the theoretically preferred `T + N + 2` component
layout with explicit g_pos and g_neg.

**Step 1: Component Layout Change**

Update `locus.py::build_locus_em_data()`:
- Component layout: `[0, n_t) mRNA + [n_t, n_t+n_nrna) nRNA +
  [n_t+n_nrna] g_pos + [n_t+n_nrna+1] g_neg`
- For each unspliced fragment, create TWO gDNA candidates:
  - g_pos candidate with deterministic strand likelihood (pos → 1.0,
    neg → ε)
  - g_neg candidate with deterministic strand likelihood (neg → 1.0,
    pos → ε)
- Remove the log(0.5) from per-fragment gDNA scoring (it is absorbed
  by the prior on φ_ℓ)

**Step 2: Beta Prior on φ_ℓ**

In the native EM M-step, add the Beta(κ_sym/2, κ_sym/2) prior term:
- After computing raw responsibility totals R_{g_pos} and R_{g_neg},
  apply:
  $$
  \theta_{g_{pos}}^{new} \propto R_{g_{pos}} + \kappa_\text{sym}/2 - 1
  $$
  $$
  \theta_{g_{neg}}^{new} \propto R_{g_{neg}} + \kappa_\text{sym}/2 - 1
  $$
  (with standard MAP normalization over the full simplex)

**Step 3: Public gDNA Reporting**

The `AbundanceEstimator` sums g_pos + g_neg to produce the collapsed
public gDNA abundance.  No change to output formats.

**Step 4: Asymmetric Constraint (Safety Net)**

For loci where the gene strand is known and nRNA is present, add a soft
penalty on the sense-dominant gDNA component:
$$
\alpha_{g_\text{sense}} = \kappa_\text{sym}/2 - \delta
$$
$$
\alpha_{g_\text{anti}} = \kappa_\text{sym}/2
$$

where δ > 0 is small and derived from the nRNA evidence strength at
the locus.  This makes sense-heavy gDNA slightly harder to sustain than
antisense-heavy gDNA, reducing the identifiability risk.

## 7. Parameter Audit: Remaining Thresholds

The current system has the following tunable parameters.  Each is
evaluated for whether it is principled or a hack.

### Calibration Parameters (calibration.py)

| Parameter | Value | Assessment |
|---|---|---|
| `_EPS` | 1e-12 | Numerical safety; principled |
| `GDNA_INIT_DENSITY_PERCENTILE` | 0.10 | Initialization heuristic; acceptable |
| `GDNA_INIT_MIN_REGIONS` | 100 | Initialization floor; acceptable |
| `max_iterations` | 50 | EM budget; acceptable |
| `convergence_tol` | 1e-4 | EM convergence; principled |
| `strand_specificity > 0.55` check in `_seed_initial_partition` | Threshold | Mild hack for seed supplementation; acceptable for init |

### Locus EM Parameters (config.py)

| Parameter | Value | Assessment |
|---|---|---|
| `prior_alpha` | 0.01 | Dirichlet; principled (sparse MAP) |
| `prior_gamma` | 1.0 | OVR coverage weighting; principled |
| `nrna_sparsity_alpha` | 0.9 | nRNA Dirichlet < 1; principled |
| `STRAND_DENOM_MIN` | 0.2 | Min 2SS−1 for strand correction; principled |

### Scoring Parameters (scoring.py)

| Parameter | Value | Assessment |
|---|---|---|
| `log(0.5)` strand penalty | -0.693 | Correct collapsed marginal |
| `log(0.01)` overhang penalty | -4.605 | Structural; principled |
| `log(0.1)` mismatch penalty | -2.303 | Structural; principled |
| `TAIL_DECAY_LP = log(0.99)` | -0.01 | FL model tail; principled |

**Assessment**: No remaining "hacks" that violate the theoretical model.
The initialization thresholds in `_seed_initial_partition()` are standard
EM initialization practice and do not affect the converged solution.

## 8. Stress Test Evidence

The calibration system (Gaussian density + Binomial strand + shape-
normalized FL) was validated against a 230-scenario stress test varying:
- Strand specificity: {0.50, 0.60, 0.70, 0.80, 0.90, 0.95}
- gDNA contamination: {0%, 10%, 30%, 60%}
- FL overlap: {none, partial, identical}
- nRNA presence: {none, low, high}

### v3 (histogram hacks) → v5 (principled models) Comparison

| Metric | v3 | v5 | Δ |
|---|---|---|---|
| Overall gDNA F1 | 0.946 | **0.987** | +0.041 |
| Overall RNA F1 | 0.980 | **0.991** | +0.011 |
| SS=0.50, gDNA F1 | 0.737 | **0.964** | +0.227 |
| SS=0.50, FL-identical, gDNA F1 | 0.287 | **0.902** | +0.615 |
| SS≥0.60, gDNA F1 (all) | ≥0.964 | **≥0.990** | — |
| Win/Tie/Loss | — | 26/203/1 | 0 regressions >0.05 |

## 9. Risk Assessment

### Low Risk
- **Region partition**: Complete, well-tested, no changes needed.
- **Calibration EM convergence**: Robust across 230 scenarios.
- **gDNA FL model**: Shape-normalized approach eliminates systematic bias.

### Medium Risk
- **Binomial → Beta-Binomial**: The Binomial works well in current stress
  tests.  The Beta-Binomial is theoretically better but may have minimal
  practical impact unless overdispersion is large.  Low regression risk.
- **Pipeline integration**: Mechanical wiring changes.  Main risk is
  ensuring the calibrated gDNA FL model matches the scoring system's
  expected format.

### Higher Risk
- **Two-strand gDNA model**: Requires C++ changes to the native EM.
  Risk of introducing bugs in the deterministic EM sort order, SQUAREM
  acceleration, or post-EM pruning logic.  Needs careful testing.
- **nRNA/gDNA identifiability**: The asymmetric constraint is not yet
  specified precisely.  Getting the penalty magnitude wrong could either
  (a) not help or (b) over-penalize legitimate sense-heavy gDNA.

## 10. Recommended Execution Order

1. **Phase I.1–I.4**: Connect calibration to pipeline (no C++ changes)
2. **Phase II**: Beta-Binomial strand model in calibration (Python only)
3. **Phase I.5**: Bridge κ_sym into existing collapsed gDNA model
4. **Stress test**: Validate end-to-end with integrated calibration
5. **Phase III**: Two-strand gDNA model (C++ changes)
6. **Ablation**: Measure contribution of each component

This order maximizes the ratio of benefit to risk: the biggest gain
(connecting calibration to pipeline) requires no C++ changes, while the
theoretically cleanest change (two-strand model) is deferred until the
integration is validated.

## 11. Summary

The calibration subsystem implements the theoretical model faithfully for
Workstreams A–E.  The primary gap is Workstream F: integration into the
production pipeline.  The implementation diverges from theory in two
specific places — Binomial vs Beta-Binomial strand model, and single vs
two-strand gDNA components — both of which have clear, mathematically
grounded remedies.  No remaining hardcoded thresholds or hacks control
the statistical behavior of the system.
