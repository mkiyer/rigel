# gDNA Calibration Redesign: Two-Component Mixture Model

## Problem Statement

The current calibration module (`calibration.py`) suffers from parameter soup:
`MIN_COVERAGE`, `SEED_FRACTION`, `MIN_SEED_REGIONS`, `FL_WEIGHT_THRESHOLD`,
`KAPPA_MIN`, `KAPPA_MAX`, `KAPPA_FALLBACK`, `KAPPA_MIN_OBS` — 8 hard-coded
domain-specific constants. These thresholds are fragile, non-portable across
libraries, and cannot adapt to the data.

The fundamental insight is that this is a **mixture model deconvolution**
problem. Regions are either *not expressed* (gDNA only) or *expressed*
(gDNA + RNA). Four convergent signals discriminate between them:

1. **Splice evidence** — spliced fragments are definitively RNA
2. **Strand balance** — gDNA is symmetric; RNA is strand-biased
3. **Coverage density** — gDNA is low, uniform; expression adds excess
4. **Fragment length** — gDNA and RNA may have different FL distributions

The algorithm should *learn from the data* which regions are which, and
estimate gDNA parameters from the deconvolution — no magic numbers.

## Biological Model

Each genomic region belongs to one of two classes:

- **Class 0 (not expressed)**: fragments come entirely from gDNA
  contamination. Coverage ~ Poisson(λ_g × L), strand ratio ~ Beta(κ/2, κ/2)
  centered at 0.5, splice rate = 0.

- **Class 1 (expressed)**: fragments come from gDNA + RNA.
  Coverage > gDNA level, strand ratio biased toward SS (for stranded
  libraries), may have spliced reads.

## Algorithm: EM Two-Component Region Deconvolution

### Phase 0: Feature Computation (no parameters)

For each region r:
- n_total, n_unspliced, n_spliced, n_pos, n_neg, length
- density = n_total / length
- strand_ratio = n_pos / n_unspliced (NaN if n_unspliced = 0)
- splice_rate = n_spliced / n_total (0 if n_total = 0)
- gene_strand, chromosome

### Phase 1: Hard Constraints

- **Spliced**: if n_spliced > 0, γ_r = 0 (definitively expressed).
  Spliced reads are indisputable RNA evidence.
- **Zero coverage**: if n_total = 0, γ_r = 1 (not expressed by default;
  gDNA simply missed this region or it's very short).

**Important**: a region with n_spliced = 0 is NOT automatically Class 0.
Single-exon genes and nascent RNA (nRNA) produce zero spliced reads but
are highly expressed.  Such regions enter Phase 3 and are classified
softly via density, strand, and FL signals.

### Phase 2: Initialization

The key insight (from the user): even in *expressed* regions, the strand
formula can estimate gDNA density. For a region with S sense and A
antisense unspliced reads and strand specificity SS:

    G = 2(A × SS − S × (1 − SS)) / (2SS − 1)

Applied globally (aggregating all unspliced reads), this gives an initial
λ_g estimate that's robust even when most regions are expressed.

**Initialization steps:**
1. λ_g = global strand-based density estimate (or median density of
   unspliced-only regions as fallback for unstranded)
2. μ_e = mean density of regions with n_spliced > 0 (or 10× λ_g as fallback)
3. κ = MoM from all unspliced-only regions' strand ratios around 0.5
4. π = 0.5 (uninformative; the EM will learn the mixing proportion)

If κ estimation fails (too few regions, variance ≥ 0.25), set κ so that
the strand LLR contributes ≈ 0 for all regions. This honestly disables
the strand signal when the data can't support it.

### Phase 3: EM Iteration

**E-step** — compute per-region posterior γ_r = P(not expressed | data_r):

For each unspliced-only region (n_spliced = 0, n_total > 0):

    log_odds = log(π / (1 − π))

    # Strand signal
    + LLR_strand(n_pos, n_unspliced, κ, SS, gene_strand)

    # Density signal
    + LLR_density(n_total, λ_g, μ_e, length)

    # Fragment length signal (after iteration 0)
    + LLR_fl(fragments_in_region, gdna_fl, rna_fl)

    γ_r = sigmoid(log_odds)

**M-step** — update all parameters from γ-weighted data:

    π = Σγ / N_eligible   (eligible = regions with n_total > 0)

    λ_g = Σ(γ × n_total) / Σ(γ × length)

    μ_e, σ_e = weighted MLE from (1-γ)-weighted log-densities:
        μ_e = Σ((1-γ) × log(d)) / Σ(1-γ)
        σ_e² = Σ((1-γ) × (log(d) − μ_e)²) / Σ(1-γ)

    κ = γ-weighted Beta-Binomial MLE (golden-section search on [0.01, 10000]):
        maximise Σ γ_r × log BetaBinom(k_r | n_r, κ/2, κ/2)
        Accurate at low per-region counts because the generative model
        inherently accounts for binomial sampling variance (no MoM).

    Update gDNA FL model from γ-weighted fragment lengths

**Convergence**: max(|Δπ|, |Δλ_g|/max(λ_g, ε), |Δκ|/max(κ, 1),
                     |Δμ_e|/max(|μ_e|, 1), |Δσ_e|/max(σ_e, ε)) < tol

### Phase 4: Post-processing

1. **Per-chromosome density**: λ_g^c from γ-weighted regions on each ref
2. **Empirical strand distribution**: γ-weighted strand ratios
3. **Final gDNA FL model**: γ-weighted fragment length histogram

## Likelihood Components

### Strand LLR (unchanged from current, well-validated)

    gDNA model: p ~ N(0.5, 0.25/(κ+1) + 0.25/n)
    RNA model:  p ~ N(p_rna, p_rna(1−p_rna)/n)

    where p_rna = SS for + genes, 1−SS for − genes, 0.5 for ambiguous

    LLR = loglik_gDNA − loglik_RNA

When SS ≈ 0.5: both models agree → LLR ≈ 0 (strand signal vanishes).
When gene_strand = 0: p_rna = 0.5 → LLR ≈ 0 (no directional info).

### Density LLR: Normal(gDNA) vs LogNormal(expressed)

Both components model density d = n_total / L:

    gDNA:      d ~ N(λ_g, λ_g / L)       [CLT of Poisson in density space]
    Expressed: d ~ LogNormal(μ_e, σ_e²)   [overdispersed expression]

The Normal approximation for gDNA is tight (var = λ_g / L), so regions
near the gDNA rate get high gDNA likelihood.  The LogNormal for expressed
regions handles extreme expression safely — a gene at 100× the mean
density gets finite log-likelihood instead of the numerical underflow
a two-Poisson model would produce.

Both λ_g and the LogNormal parameters (μ_e, σ_e) are estimated from
data via the M-step.  The M-step computes μ_e and σ_e² from the
(1−γ)-weighted log-densities of regions.  The output `expressed_density`
is derived as the LogNormal mean: exp(μ_e + σ_e² / 2).

### Fragment Length LLR (iterative, from iteration 1+)

    LLR_fl = Σ [log P(f | FL_gDNA) − log P(f | FL_RNA)]

where the sum is over all unspliced fragments in the region, and FL_gDNA
is the γ-weighted gDNA FL model from the previous iteration.

FL_RNA comes from the existing spliced-read FL model (already trained by
the pipeline). The gDNA FL model is built from the calibration's own
γ-weighted histogram.

## Outputs

```python
@dataclass(frozen=True)
class GDNACalibration:
    # Per-region posteriors
    region_posteriors: np.ndarray     # γ_r ∈ [0, 1], P(not expressed)

    # Global gDNA density
    gdna_density_global: float        # λ_g (frags/bp)

    # Per-reference density
    gdna_density_per_ref: dict[str, float]

    # Strand symmetry
    kappa_sym: float                  # Beta(κ/2, κ/2) concentration

    # Fragment length
    gdna_fl_model: FragmentLengthModel

    # Mixing proportion
    mixing_proportion: float          # π = P(not expressed)

    # Expressed density (for density LLR)
    expressed_density: float          # μ_e (frags/bp)

    # Convergence
    n_iterations: int
    converged: bool

    # Diagnostics
    region_stats: dict | None = None
    iteration_history: list[dict] | None = None
```

## Eliminated Parameters

| Old Parameter       | Replacement                                          |
|---------------------|------------------------------------------------------|
| MIN_COVERAGE        | No threshold. Regions with n_total=0 → γ=1 trivially |
| SEED_FRACTION       | EM learns π (mixing proportion) from data            |
| MIN_SEED_REGIONS    | No threshold. EM works with any number of regions    |
| FL_WEIGHT_THRESHOLD | Posteriors γ naturally weight FL contributions        |
| KAPPA_MIN           | κ ≥ 0 from MoM (no arbitrary lower bound)            |
| KAPPA_MAX           | No upper bound; large κ means low overdispersion     |
| KAPPA_FALLBACK      | If κ unestimable, disable strand signal (LLR=0)     |
| KAPPA_MIN_OBS       | No threshold; MoM gives best estimate from available data |

**Remaining parameters** (standard EM, not domain-specific):
- `convergence_tol = 1e-3` — relative change threshold for convergence
- `max_iterations = 50` — safety bound for EM iterations

## Edge Cases

- **No gDNA**: λ_g → 0, π → 0, all regions classified expressed. Correct.
- **All gDNA**: π → 1, all regions not expressed. Correct.
- **Unstranded (SS ≈ 0.5)**: strand LLR ≈ 0, relies on density + FL. Honest.
- **Very few regions**: EM still works; posteriors reflect uncertainty.
- **No spliced reads**: expressed model initialized from density alone.
- **All regions zero coverage**: bailout → no gDNA detected.

## Connection to Downstream EM

The calibration outputs feed into `compute_eb_gdna_priors()` in locus.py:

1. `gdna_density_global` → replaces the strand-formula global estimate
2. `gdna_density_per_ref` → replaces per-ref strand-formula estimates
3. `kappa_sym` → informs strand symmetry model in fragment scoring
4. `gdna_fl_model` → provides FL discrimination in fragment scoring
5. `region_posteriors` → diagnostic; can inform per-locus priors

The 3-level shrinkage (global → ref → locus) continues to work, but now
with calibration-supplied hyperparameters learned from the data.
