# Calibration v3: Aggregate-First Empirical EM

## Motivation

The current calibration module imposes parametric assumptions (Normal
density, LogNormal density, Beta-Binomial strand κ) that keep breaking:

1. **Sparse regions**: tiny regions (10–100 bp) with 0–5 fragments make
   per-region variance estimates collapse.
2. **Parametric assumptions without evidence**: we don't know the true
   shape of the gDNA or expressed density distribution.  Imposing Normal
   vs LogNormal kept failing in edge cases.
3. **RNA expression varies dramatically**: a single global mean λ_E
   cannot represent the expressed density distribution.  Gene expression
   spans 4+ orders of magnitude.  Mean and variance are heavily
   influenced by outliers and are difficult to estimate robustly.

The fix: **pool fragments by class, build global *empirical*
distributions (histograms of log-density), then score individual regions
against those empirical models.**  No parametric density assumptions at
all.

### Guiding principles

- **Work in log-space.**  Sequencing is inherently a log-space problem
  (PCR amplification, expression dynamic range).  Log-density is the
  natural scale.
- **Prefer empirical distributions.**  Build histograms / eCDFs from the
  data.  No parametric forms unless we have strong empirical evidence for
  them.
- **Exclude zero-count regions from estimation.**  Regions with
  `n_total = 0` carry no information about instrument fragment
  distributions.  Exclude from M-step aggregates and π estimation.
  Their posterior defaults to the prior π (see "Zero-Count Regions"
  section).
- **Use global aggregates for stability.**  Pool across regions within
  each class to build well-powered empirical distributions, even when
  individual regions are sparse.
- **Freeze what can be frozen.**  The expressed-class FL model comes
  from spliced fragments (gold standard) and is frozen before EM.  Bin
  edges for density histograms are fixed after initialization.  This
  guarantees monotonic EM convergence.

---

## Design Overview

### Two-class mixture model

Every region *with at least one observed fragment* belongs to one of two
classes:

- **G** = "not expressed" (gDNA only).  Uniform background.
- **E** = "expressed" (RNA + gDNA).

Regions with `n_total = 0` are **excluded from estimation** — not
included in M-step aggregates or π computation.  Their posterior
defaults to the prior: `γ_r = π` (the current mixing proportion).
This is the statistically correct uninformative assignment: with zero
evidence, the posterior equals the prior.

The EM iterates:

1. **M-step (aggregate)**: pool all fragments from G-weighted and
   E-weighted regions (n_total > 0 only) to build *global* empirical
   distributions.
2. **E-step (per-region)**: score each region (n_total > 0) against the
   global models to update γ_r = P(G | data_r).

### Three evidence channels

For each eligible region *r*:

```
LLR_r = LLR_density(r) + LLR_strand(r) + LLR_fl(r)
```

Then: `γ_r = σ(log(π/(1−π)) + LLR_r)` for soft regions.

Hard constraint: `n_spliced > 0 → γ = 0` (definitively expressed).
All other eligible regions are classified softly.

---

## Detailed Channel Design

### Channel 1: Density — Empirical histograms of log-density

**Log-density definition**:

```
log_d_r = log((n_total_r + 1) / L_r)
```

A pseudocount is **mandatory**.  Standard RNA-seq tools (DESeq2,
edgeR) apply pseudocounts before log transforms for the same reason:
without it, the derivative of log at counts near zero causes massive
artificial separation between $n=1$ and $n=2$ that is pure Poisson
sampling noise.

**Implementation detail**: the naïve `log((n+1)/L)` has a subtle
problem — regions with zero reads but different lengths get different
log-densities (`(0+1)/100` vs `(0+1)/10000`), artificially spreading
them across the histogram.  The correct approach adds a global
epsilon *after* the length-normalised density:

```
log_d_r = log(n_total_r / L_r + ε)
```

where `ε = 1 / median(L)` across all eligible regions.  This ensures
all zero-count regions pile up at exactly the same baseline log-density,
creating a sharp, clean left edge for the gDNA histogram.  Regions
with even 1 fragment are already well above ε for typical region
lengths, so the pseudocount has negligible effect on non-zero regions.

Only regions with `n_total > 0` AND `L_r > 0` participate in the EM
(zero-count regions are excluded from estimation and assigned γ = π).

**Model**: two empirical histograms of log-density — one for each class.

- **H_G**: histogram of `log_d_r` values, weighted by `γ_r`.
  Captures the log-density distribution of not-expressed (gDNA) regions.
- **H_E**: histogram of `log_d_r` values, weighted by `(1 − γ_r)`.
  Captures the log-density distribution of expressed regions.

Both use the same infrastructure as strand histograms:
`build_histogram` + Laplace smoothing + `_log_density_from_histogram`.

**Why this is better than Poisson or LogNormal**:

The Poisson model `n ~ Poisson(λ × L)` assumes a single constant rate
per class.  For gDNA this is reasonable — gDNA density is roughly
uniform genome-wide.  But for the expressed class, expression levels
vary by 4+ orders of magnitude (a few TPM to thousands of TPM).  A
single λ_E cannot capture this.  A Poisson LLR with λ_E = global mean
would heavily penalise both very-low-expression and very-high-expression
regions (both are far from the mean).  This is wrong — a region with
density 100× the mean is *extremely* expressed, not gDNA.

The empirical histogram avoids this entirely.  It captures the actual
shape of the expressed density distribution — heavy right tail,
multimodal peaks, whatever the data shows — without any parametric
assumptions.

**LLR**:

```
LLR_density(r) = log P(log_d_r | H_G) − log P(log_d_r | H_E)
```

This is just a histogram bin lookup for each class, exactly like the
strand channel.  Reuses `_log_density_from_histogram`.

**M-step**: build H_G and H_E from γ-weighted log-density values.

```python
log_d = np.log((n_total[eligible] + 1) / region_length[eligible])
H_G = build_histogram(log_d, weights=γ[eligible])
H_E = build_histogram(log_d, weights=(1 − γ)[eligible])
```

The bin range is **auto-ranged at initialization, then frozen** for
all subsequent EM iterations.  Compute `[min(log_d) − 0.5,
max(log_d) + 0.5]` from all eligible regions during Phase 3, establish
the bin edges (n_bins ≈ 50), and never change them.

*Rationale*: EM guarantees monotonic convergence only if the parameter
space (here, the histogram support) is fixed.  If bin edges shift
between iterations, likelihoods from iteration $k$ and $k+1$ are not
comparable and convergence checks may oscillate or fail.

**NumPy implementation**:
```python
# Phase 3: compute and freeze bin edges
hist_min = np.min(log_d[eligible]) - 0.5
hist_max = np.max(log_d[eligible]) + 0.5
density_bin_edges = np.linspace(hist_min, hist_max, num=n_bins + 1)

# M-step: use frozen edges for both classes
H_G_raw, _ = np.histogram(log_d[eligible], bins=density_bin_edges,
                           weights=γ[eligible])
H_E_raw, _ = np.histogram(log_d[eligible], bins=density_bin_edges,
                           weights=(1 - γ)[eligible])
```

Using `np.histogram` with explicitly provided `bins` guarantees that
H_G and H_E have identical support and align perfectly for the
element-wise LLR subtraction in the E-step.

**Global scalar summaries** (for output / diagnostics, NOT for the LLR):

```python
λ_G = Σ(γ × n_total) / Σ(γ × L)         # aggregate gDNA density
λ_E = Σ((1−γ) × n_total) / Σ((1−γ) × L) # aggregate expressed density
```

These are reported in `GDNACalibration` for downstream use but the
per-region classification score uses the empirical histograms, not
these scalars.

**What about gDNA?**  For gDNA, a single Poisson rate might actually be
adequate (gDNA density is roughly constant).  But using an empirical
histogram hurts nothing — if gDNA really is a tight unimodal peak
around `log(λ_G)`, the histogram will reflect that.  And if there's
any variation (e.g., GC bias, mappability effects), the histogram
captures it automatically.  So we use the same histogram approach for
both classes for simplicity and consistency.

### Channel 2: Strand (sense fraction)

Unchanged from v3 plan.

**Observation**: region *r* has `k_r` sense-strand fragments out of
`m_r` unspliced fragments, with gene strand `s_r ∈ {+1, −1, 0}`.

**Sense fraction**: gene-strand-normalized metric.
```
sf_r = 1 − (n_pos_r / m_r)  if s_r = +1    (R1-antisense convention)
sf_r = n_pos_r / m_r         if s_r = −1
sf_r = NaN                   if s_r = 0 or m_r < 2
```

Regions with `sf = NaN` get `LLR_strand = 0` (no strand signal).

**Model**: empirical histograms (already implemented and working).
- gDNA histogram: γ-weighted sense fractions, forced symmetric about 0.5.
- RNA histogram: (1−γ)-weighted sense fractions.
- Both use Laplace smoothing.

**LLR**:
```
LLR_strand(r) = log P(sf_r | gDNA_hist) − log P(sf_r | RNA_hist)
```

**First iteration bootstrap**: binomial LLR with `p_gDNA = 0.5` vs
`p_RNA` from gene strand + SS.

### Channel 3: Fragment length

**Observation**: region *r* has fragment lengths `{fl_1, ..., fl_k_r}`.

**Model**: two `FragmentLengthModel` histograms.
- `FL_E`: **Frozen.**  Built from spliced fragments *before* EM begins.
  Spliced reads are an unambiguous gold standard for the RNA fragment
  length distribution (both mature and nascent RNA produce spliced
  fragments with the true RNA FL profile).  This model is never updated
  during EM.
- `FL_G`: γ-weighted fragment lengths from all eligible regions.
  Updated each M-step.

**Why FL_E is frozen from spliced reads**: "expressed" regions contain
a mixture of RNA and gDNA fragments.  If we built FL_E from (1−γ)-
weighted fragments of unspliced reads in expressed regions, it would be
contaminated by gDNA fragment lengths in proportion to the gDNA
fraction.  Spliced reads are *pure* RNA by definition, so their FL
distribution is the correct reference.

**Diagnostic FL_E tracking** (optional): for monitoring purposes, we
can also build a (1−γ)-weighted FL histogram from unspliced fragments
in expressed regions.  This "impure" expressed FL model is NOT used
for scoring — it's purely diagnostic.  The difference between the
frozen (spliced) FL_E and the diagnostic (unspliced) FL_E reflects
the gDNA contamination level in expressed regions.

**LLR**:
```
LLR_fl(r) = Σ_{i=1}^{k_r} [log FL_G(fl_i) − log FL_E(fl_i)]
```

**Note on LLR validity**: even though FL_E is pure RNA and FL_G is
built from "not expressed" regions (which are essentially pure gDNA),
the LLR optimises for the *difference* in fragment length between
classes, which is theoretically valid for classification.

Enabled from iteration 2 onward (FL_G needs initial γ weights).
FL_E is available from the start (frozen from spliced reads).

---

## EM Algorithm

### Inputs

- `region_counts` (DataFrame): per-region `n_unspliced_pos`,
  `n_unspliced_neg`, `n_spliced_pos`, `n_spliced_neg`.
- `fl_table` (DataFrame): per-fragment `region_id`, `frag_len`.
- `region_df` (DataFrame): region metadata with `tx_pos`, `tx_neg`,
  `length`.
- `strand_specificity` (float): library SS ∈ [0.5, 1.0].

### Preprocessing

1. **Compute summary stats** (`compute_region_stats`): n_total, n_pos,
   n_neg, n_spliced, n_unspliced, region_length, gene_strand, etc.

2. **Compute log-density** for all regions:
   ```
   log_d = log((n_total + 1) / region_length)     # for L > 0
   ```

3. **Compute sense fractions** (`compute_sense_fraction`): per-region sf.

4. **Identify eligible regions**: `eligible = (n_total > 0) & (L > 0)`.
   Only these participate in the EM.  Zero-count regions get `γ = π`
   (prior — see "Zero-Count Regions" below).

### Zero-Count Regions

Regions with `n_total = 0` provide zero information about the fragment
properties of either class.  They are:
- **Excluded from M-step** aggregates (no fragments to pool).
- **Excluded from π** computation (no evidence to update the prior).
- **Assigned `γ = π`** (the current mixing proportion).  This is the
  Bayesian default: with no observed data, the posterior equals the
  prior, $P(G \mid \emptyset) = P(G) = \pi$.

We do NOT use NaN.  NaN propagates silently through numpy/pandas/C++
and invites bugs.  The prior π is semantically correct and numerically
safe.

### Initialization

The seed partition determines initial empirical distributions.
Bad seeds → bad histograms → EM stuck in wrong mode.  We use a
**two-phase** initialization that combines splice evidence, strand
information, and density.

#### Phase 1: Expressed seed (high confidence)

```python
expressed_seed = n_spliced > 0   # γ = 0 (definitively expressed)
```

This is unambiguous.  From these regions, build:
- **Expressed log-density eCDF**: `F_E(x) = P(log_d ≤ x | expressed)`.
- **Expressed strand eCDF**: distribution of sense fractions.

These are well-powered because most datasets have thousands of
regions with spliced reads.

#### Phase 2: gDNA seed (bootstrapped from density + strand)

We need to identify "not expressed" regions among the unspliced-only
set.  The key insight: **not-expressed regions should have the lowest
density**, because they contain only gDNA background, while expressed
regions have gDNA + RNA.

Using the expressed eCDF from Phase 1:

```python
# For each unspliced-only region, compute P(expressed | density)
# using the expressed density eCDF as reference
p_expressed = F_E(log_d_r)   # quantile in the expressed distribution

# Regions whose density falls well below the expressed distribution
# are unlikely to be expressed → strong gDNA candidates
gdna_seed = (n_spliced == 0) & (p_expressed < τ_density)
```

where `τ_density` is a percentile threshold — regions whose log-density
is below the `τ_density`-th percentile of the expressed distribution.

**Choosing τ_density**: we want high-confidence gDNA seeds, so a
conservative threshold.  `τ_density = 0.10` means: "if this region's
density is lower than 90% of expressed regions, it's probably gDNA."
This is conservative because:
- Lowly-expressed genes do exist at the bottom of the expressed
  distribution, but they are a minority.
- Because EM is self-correcting, a false *negative* (excluding a true
  gDNA region from the seed) is much safer than a false *positive*
  (poisoning the initial gDNA histogram with lowly-expressed RNA).

**Minimum seed size**: with real data we have hundreds of thousands of
regions, so the 10th percentile yields thousands of gDNA seeds.  But
for small datasets or unusual libraries, we add a minimum count
guarantee:

```python
GDNA_INIT_DENSITY_PERCENTILE = 0.10  # default
GDNA_INIT_MIN_REGIONS = 100          # minimum gDNA seeds

# Initial candidate set: unspliced-only below percentile
candidates = (n_spliced == 0) & (p_expressed < τ_density)
n_candidates = candidates.sum()

# If not enough, expand upward until we reach MIN_REGIONS
if n_candidates < GDNA_INIT_MIN_REGIONS:
    # Sort unspliced-only regions by p_expressed, take bottom N
    unspliced_only = (n_spliced == 0) & eligible
    sorted_idx = np.argsort(p_expressed[unspliced_only])
    candidates = sorted_idx[:GDNA_INIT_MIN_REGIONS]
```

Both thresholds are tunable parameters in `calibrate_gdna`.

**Supplementary strand filter** (when SS > 0.55): we can additionally
include symmetric-strand regions as gDNA seeds:

```python
# Strand-symmetric unspliced regions on known gene-strand
strand_gdna = (n_spliced == 0) & (abs(sf - 0.5) < 0.1) & (gene_strand != 0)
gdna_seed = gdna_seed | strand_gdna
```

This helps with cases where density alone is ambiguous (e.g., a low-
expression gene has similar density to gDNA, but its strand bias
reveals it as expressed).

**Combining**: the initial γ vector is:

```python
γ[n_total == 0]       = π_init  # no evidence → defaults to prior
γ[expressed_seed]     = 0.0     # definitively expressed
γ[gdna_seed]          = 1.0     # probable gDNA
γ[remaining eligible] = 0.5     # ambiguous
```

#### Phase 3: Build initial empirical distributions

From the Phase 2 seeds:

```python
H_G_init = build_histogram(log_d[gdna_seed], weights=1)
H_E_init = build_histogram(log_d[expressed_seed], weights=1)
strand_G_init = build_strand_histogram(sf[gdna_seed], weights=1, symmetric=True)
strand_E_init = build_strand_histogram(sf[expressed_seed], weights=1)
π_init = gdna_seed.sum() / eligible.sum()
```

**Fallback if gDNA seeds are empty** (e.g., unstranded library with
no low-density unspliced regions):
- Use the bottom 20th percentile of all eligible regions' log-density as
  a proxy gDNA distribution.
- Set `π_init = 0.3` (weak prior toward some gDNA).

**Fallback if expressed seeds are empty** (e.g., no spliced regions
observed — unusual but possible):
- Use the top 50th percentile of all eligible regions' log-density.
- Set `π_init = 0.5`.

#### Pristine sample detection

Some samples are essentially free of gDNA contamination.  In these
cases the gDNA seed will be very weak: few candidates, and those
candidates will have density not much lower than the expressed class.

After initialization, we check:

```python
# If initial π is very small AND the gDNA density histogram
# overlaps heavily with the expressed histogram → likely pristine
if π_init < 0.02 or gdna_seed.sum() < GDNA_INIT_MIN_REGIONS // 2:
    # Mark as pristine: run EM but expect π → ~0
    # The EM will naturally converge to π ≈ 0 and all soft regions
    # will get γ ≈ 0 (expressed).  No special intervention needed,
    # but we flag it in the output for downstream awareness.
    diagnostics['pristine_sample'] = True
```

The key insight: the EM is self-correcting.  If there truly is no gDNA,
the M-step density histogram for G will converge toward a degenerate
distribution and π will shrink toward zero.  We don't need a separate
code path — just let the EM run and report the result.  The diagnostic
flag alerts downstream code that gDNA priors may be uninformative.

### Iteration

```
for iteration in 1..max_iterations:
    # E-step: score each eligible region
    for each eligible region r:
        if n_spliced_r > 0:  γ_r = 0;  continue   # hard expressed

        LLR = LLR_density(r) + LLR_strand(r) + LLR_fl(r)
        γ_r = σ(log(π/(1−π)) + LLR)

    # M-step: re-estimate globals from soft assignments, eligible only
    π = Σ γ_r / n_eligible
    H_G = build_histogram(log_d[eligible], weights=γ[eligible])
    H_E = build_histogram(log_d[eligible], weights=(1−γ)[eligible])
    strand_G = build_strand_histogram(sf, weights=γ, symmetric=True)
    strand_E = build_strand_histogram(sf, weights=1−γ)
    FL_G = build_fl_model(fl_table, weights=γ[region_id])
    # FL_E is FROZEN from spliced reads — not updated

    # Scalar summaries (for output, not used in LLR)
    λ_G = Σ(γ × n) / Σ(γ × L)
    λ_E = Σ((1−γ) × n) / Σ((1−γ) × L)

    # Convergence check
    if max(Δπ, Δλ_G, Δλ_E) < tol:  break
```

### Outputs (GDNACalibration — interface preserved)

| Field | Source |
|-------|--------|
| `region_posteriors` | final γ array (π for n_total=0 regions) |
| `gdna_density_global` | final λ_G (scalar aggregate) |
| `gdna_density_per_ref` | per-chromosome λ_G |
| `kappa_sym` | post-hoc `estimate_kappa_sym(stats, final_γ)` |
| `gdna_fl_model` | final FL_G histogram |
| `mixing_proportion` | final π |
| `expressed_density` | final λ_E (scalar aggregate) |
| `n_iterations`, `converged` | loop metadata |
| `region_stats`, `iteration_history` | diagnostics |

**Note on region_posteriors**: previously, n_total=0 regions got γ=1
(assumed gDNA).  Now they get `γ = π` (the final mixing proportion).
This is a clean numeric value — no NaN propagation — and is the
Bayesian correct default for zero-evidence regions.

---

## Per-Region Likelihood: Why Empirical Histograms

The key modelling question is: how do we score a region with log-density
`log_d_r` under each class?

**Options considered**:

| Model | Pros | Cons |
|-------|------|------|
| Poisson(λ×L) | Exact for constant rate. No variance param. | Single rate per class. Penalises both low and high expression equally relative to mean. Cannot represent 4-OOM expression range. |
| Normal on log-density | Simple. Two params. | Still parametric. Assumes symmetric log-density. |
| LogNormal on density | Handles heavy tail. | Two params. Underflow at extremes. |
| **Empirical histogram on log-density** | **Non-parametric. Captures actual shape. No distributional assumptions. Handles multimodality, skew, heavy tails.** | **Need enough data to populate bins (solved by pooling across regions).** |

**The empirical histogram wins** because:

1. **RNA expression is not a single rate.**  It spans 4+ orders of
   magnitude.  No single-parameter model can represent this.  A Poisson
   with λ_E = global mean would give a huge LLR for a highly-expressed
   gene (far from the mean), incorrectly suggesting it's NOT consistent
   with the expressed class.  The histogram correctly shows that some
   expressed regions do have very high density.

2. **gDNA density might not be perfectly constant either.**  GC bias,
   mappability variation, and repeat content cause gDNA density to vary.
   The histogram captures this without needing to model it.

3. **Consistent with strand and FL channels.**  All three channels now
   use the same architecture: empirical histograms + LLR from bin
   lookup.  One set of building blocks.

4. **Working in log-space is natural.**  Sequencing involves PCR
   amplification; expression dynamic range is inherently logarithmic.
   Log-density is the right scale for comparing gDNA background to
   expressed regions.

5. **Pooling across regions provides power.**  Even if individual
   regions are sparse, the global histograms aggregate thousands of
   regions.  With ~50 bins and thousands of observations, each bin
   is well-populated.

---

## What Changes vs Current Code

### Removed

- **`_compute_density_llr`** (Normal vs LogNormal): replaced by
  empirical histogram lookup.
- **`_estimate_lognormal_params`**, **`_init_lognormal_from_spliced`**:
  no LogNormal anywhere.
- **`_init_gdna_density_from_strand`**: replaced by histogram-based
  init from seed partition.
- **`_estimate_density`**: replaced by M-step aggregate.
- **`log_mu_e`, `log_sigma_e`** parameters: gone.  Replaced by density
  histograms.

### Kept (unchanged)

- **`compute_region_stats`**: still computes per-region summaries.
- **`compute_sense_fraction`**: sense fraction computation.
- **`build_strand_histogram`**, **`_log_density_from_histogram`**,
  **`_compute_empirical_strand_llr`**, **`_compute_strand_llr_binomial`**:
  entire strand histogram infrastructure.  Also reused for density
  histograms (same `build_histogram` + `_log_density_from_histogram`
  pattern).
- **`build_gdna_fl_model`**: builds weighted FL histogram.  Generalise
  or add parallel `build_fl_model`.
- **`_compute_fl_llr`**: per-fragment FL LLR.
- **`FragmentLengthModel`**: upstream class with `log_likelihood`,
  `finalize`, weighted histogram.  Used for both FL_G and frozen FL_E.
- **`estimate_kappa_sym`**, **`_golden_section_max`**,
  **`_beta_binom_loglik`**, **`_vec_lgamma`**: post-hoc κ.
- **`GDNACalibration`** dataclass: same fields, same interface.

### New

- **`build_density_histogram(log_d, weights, n_bins)`**: builds an
  empirical histogram on a data-dependent domain (not fixed [0,1] like
  strand).  Returns `(bin_edges, density)`.  Internal helper, same
  pattern as `build_strand_histogram` but with auto-ranging.
- **`_compute_density_llr_empirical(log_d, H_G, H_E)`**: histogram
  bin lookup for density channel.  Thin wrapper around
  `_log_density_from_histogram`.
- **`_seed_initial_partition(stats, log_d, sense_frac, ss)`**: the
  two-phase initialization (splice seed → density eCDF → gDNA seed).

### Modified

- **`_e_step`**: takes `(stats, π, log_d, sense_frac, H_G, H_E,
  strand_G, strand_E, FL_G, FL_E)`.  All three LLRs are histogram
  lookups.
- **`_m_step`**: estimates `(π, λ_G, λ_E, H_G, H_E, strand_G,
  strand_E)` from γ-weighted aggregates.
- **`calibrate_gdna`**: new initialization + simplified EM loop.

---

## Convergence

Primary convergence check: `|Δπ| < tol`.

The mixing proportion π is the most globally robust indicator that
mass has stopped shifting between the two components.  We also track
λ_G and λ_E for diagnostics, but **do not gate convergence on λ_E**.

*Why not λ_E*: because λ_E (arithmetic mean of expressed density) is
hypersensitive to outliers.  A single highly-expressed region changing
its γ from 0.99 to 0.999 can shift λ_E enough to prevent convergence,
even when all histograms are completely stable.  π is immune to this.

Additionally, track the total observed-data log-likelihood:

```
LL = Σ_r log[ π × P(x_r | G) + (1−π) × P(x_r | E) ]
```

If LL does not increase between iterations (within floating-point
tolerance), the EM has converged regardless of what π is doing.

**Convergence criterion**:
```python
converged = (abs(Δπ) < tol) or (ΔLL < tol_ll)
```

where `tol = 1e-4` (default) and `tol_ll = 1e-6`.

Expected: 5–15 iterations.

---

## Implementation Steps

1. **Add `build_density_histogram`**: auto-ranging histogram builder
   for log-density values.
2. **Add `_compute_density_llr_empirical`**: thin wrapper.
3. **Add `_seed_initial_partition`**: two-phase init.
4. **Rewrite `_e_step`**: all three channels use histogram LLR.
5. **Rewrite `_m_step`**: aggregates + histogram construction.
6. **Rewrite `calibrate_gdna`**: new init + simplified EM loop.
7. **Remove**: `_compute_density_llr`, `_estimate_lognormal_params`,
   `_init_lognormal_from_spliced`, `_init_gdna_density_from_strand`,
   `_estimate_density`.
8. **Update tests**: rewrite `test_calibration.py` for new API.
9. **Run stress tests**: 64-scenario sweep + analysis.

---

## What Stays the Same

- **The public interface to `calibrate_gdna`**: same arguments, same
  `GDNACalibration` output.
- **Downstream pipeline**: fully decoupled (pipeline.py does not use
  calibration.py).
- **Strand histogram infrastructure**: reused as-is and extended to
  density.
- **FL model infrastructure**: passes through.
- **Region partition and counting**: upstream, unchanged.
- **κ estimation**: post-hoc diagnostic, untouched.

---

## Resolved Design Decisions

1. **Log-density pseudocount**: **Use `log(n/L + ε)` where
   `ε = 1/median(L)`.**  Adds a global constant after length
   normalisation so all zero-count regions share the same baseline
   log-density.  Avoids the naïve `log((n+1)/L)` pitfall where regions
   of different lengths get artificially different log-densities at
   n=0.  Same stabilisation rationale as DESeq2/edgeR pseudocounts.

2. **Density histogram bin range**: **Auto-range at initialization,
   then freeze.**  Compute `[min(log_d) − 0.5, max(log_d) + 0.5]`
   from all eligible regions during Phase 3 and lock the bin edges for
   the entire EM run.  EM guarantees monotonic convergence only when the
   parameter space (including histogram support) is fixed.

3. **gDNA seed threshold**: **Dual-parameter: `GDNA_INIT_DENSITY_
   PERCENTILE = 0.10` and `GDNA_INIT_MIN_REGIONS = 100`.**  If the
   percentile-based filter yields fewer than `MIN_REGIONS` seeds, expand
   upward until the minimum is met.  Both are tunable.  Conservative
   default is correct: false negatives (missing true gDNA) are safe
   because EM self-corrects; false positives (RNA in gDNA seed) poison
   the initial gDNA histogram.

4. **Zero-count region posteriors**: **Output `γ = π`, not NaN.**  NaN
   propagates silently through numpy/pandas/C++ and invites bugs.  The
   Bayesian correct answer for zero evidence is `P(G|∅) = P(G) = π`.

5. **FL model for expressed class**: **Freeze FL_E from spliced
   fragments.**  Spliced reads are the gold standard — pure RNA with no
   gDNA contamination.  The frozen model is computed upstream before
   calibration and passed in.  FL_G is built from γ-weighted fragments
   during EM.  Optionally track a diagnostic (1−γ)-weighted unspliced FL
   for monitoring, but it is NOT used for scoring.

---

## Laplace Smoothing

All empirical histograms use Laplace smoothing before normalising to
probabilities.  This prevents infinite LLRs when a region falls into a
bin that has zero observed weight in one class.

```python
# After building raw histogram counts:
hist_smoothed = hist_raw + α
prob = hist_smoothed / hist_smoothed.sum()
```

**Choice of α**: a value of `1.0` (full Laplace) or `0.5` (Jeffreys)
added to every bin before normalisation.  This ensures that even bins
with zero weight get a small but non-zero probability.

*Example*: a highly expressed region falls into a density bin where
`H_G` has exactly 0 observations.  Without smoothing,
`log P(log_d | H_G) = −∞` and the LLR would be `−∞` regardless of
strand or FL evidence.  With α = 1.0 and 50 bins, the smoothed floor
is `1/(total_weight + 50)`, giving a finite (but large) LLR penalty.

The strand histograms already use this pattern (`_STRAND_SMOOTH`).
The density histograms will use the same approach with a configurable
`_DENSITY_SMOOTH` constant.
