# gDNA Calibration via Iterative Seed-Extend

## Problem

Rigel needs two global gDNA nuisance parameters:

1. **κ_sym**: Beta concentration for strand symmetry (overdispersion of
   gDNA strand ratios around 0.5)
2. **gDNA fragment-length distribution**: separate from the RNA FL model

The chicken-and-egg: identifying gDNA-dominated regions requires κ_sym,
but estimating κ_sym requires knowing which regions are gDNA-dominated.

## Key Biological Insight

In a typical RNA-seq experiment, >50% of annotated genes are not expressed.
All regions (exonic + intronic) of these unexpressed genes contain purely
gDNA fragments. This is an abundant, high-quality signal — far superior to
intergenic regions, which are often unmappable/repetitive and absent in
hybrid-capture data.

## Available Signals

| Signal | Reliability | Use Case |
|--------|------------|----------|
| **Spliced counts** | Definitive: any spliced read = RNA | All data types |
| **Strand ratio** | Powerful when SS > ~0.65: gDNA→0.5, RNA→SS | Stranded data |
| **Coverage density** | Weakest, but sometimes the only signal | Unstranded data |

## Algorithm

### Phase 1: Per-Region Summary Statistics

From `region_counts` (4 columns) + `region_df` (annotation flags):

- `n_unspliced = unspliced_pos + unspliced_neg`
- `n_spliced = spliced_pos + spliced_neg`
- `strand_ratio = unspliced_pos / n_unspliced`
- `splice_rate = n_spliced / (n_unspliced + n_spliced)`
- `gene_strand`: +1 (sense gene only), -1 (antisense only), 0 (both/neither)
- `density = n_total / region_length`

### Phase 2: Seed Selection

Seed = regions where we are highly confident gDNA dominates.

**Common criteria** (both data types):
- Zero spliced fragments (`n_spliced == 0`)
- Sufficient coverage (`n_unspliced >= min_coverage`)

**Stranded data** (SS > 0.65): additionally require strand symmetry:
- `z = |2 * n_pos - n_unspliced| / sqrt(n_unspliced)`
- Require `z < z_threshold` (e.g., 1.5 for two-sided p ≈ 0.13)

**Unstranded data** (SS ≤ 0.65): additionally require low coverage:
- Coverage density (`n_total / length`) in bottom percentile (e.g., 50th)
  among qualifying zero-spliced regions

### Phase 3: Initial Parameter Estimation from Seed

**κ_sym** via Method-of-Moments (forced mean = 0.5):

```
observed_var = Σ(w_r * (p_r - 0.5)²) / Σ(w_r)
sampling_var = 0.25 * Σ(w_r / n_r) / Σ(w_r)    # binomial noise
true_var     = max(observed_var - sampling_var, ε)
κ_sym        = 0.25 / true_var - 1
```

where `p_r = n_pos_r / (n_pos_r + n_neg_r)` and `w_r` = seed weight.

**gDNA FL**: weighted `FragmentLengthModel` from FL observations in seed
regions.

**gDNA density** (for unstranded mode): weighted mean density across seed
regions.

### Phase 4: Score All Regions → Per-Region gDNA Weights

**Strand-based score** (stranded data, normal approximation):

For each region with strand_ratio `p`, unspliced count `n`:

```
var_gdna = 0.25/(κ+1) + 0.25/n       # overdispersion + sampling
var_rna  = p_rna*(1-p_rna)/n          # RNA sampling only
ll_gdna  = -0.5 * ((p-0.5)²/var_gdna + ln(var_gdna))
ll_rna   = -0.5 * ((p-p_rna)²/var_rna + ln(var_rna))
```

where `p_rna = SS` for + strand gene, `1-SS` for − strand gene, `0.5`
for ambiguous.

Strand weight = sigmoid(ll_gdna − ll_rna).

**Density-based score** (unstranded or supplementary):

```
expected_gdna = gdna_density * region_length
z_density = (n_unspliced - expected_gdna) / sqrt(max(expected_gdna, 1))
density_weight = sigmoid(-z_density)
```

**Splice adjustment**: `w_r *= (1 - splice_rate)`

**Combined weight**:
- Stranded: `w_r = strand_weight * (1 - splice_rate)`
- Unstranded: `w_r = density_weight * (1 - splice_rate)`

### Phase 5: Re-estimate and Iterate

1. Re-estimate κ_sym from weighted strand ratios (Phase 3 formula)
2. Re-estimate gDNA FL from weighted FL observations
3. Re-estimate gDNA density from weighted regions
4. Re-score all regions (Phase 4)
5. Check convergence: `|κ_new - κ_old| / max(κ_old, 1) < tol`
6. Repeat until converged or max iterations

## Outputs

```python
@dataclass(frozen=True)
class GDNACalibration:
    kappa_sym: float                   # global κ
    gdna_fl_model: FragmentLengthModel # global gDNA FL
    region_weights: np.ndarray         # (n_regions,) gDNA weights ∈ [0,1]
    n_iterations: int
    converged: bool
```

## API

```python
def calibrate_gdna(
    region_counts: pd.DataFrame,   # from count_region_evidence()
    fl_table: pd.DataFrame,        # from count_region_evidence()
    region_df: pd.DataFrame,       # from index.region_df
    strand_specificity: float,     # from strand model
    **kwargs,                      # tuning parameters
) -> GDNACalibration:
```

## Design Decisions

1. **Replaces** existing `_compute_intergenic_density()` with per-region
   computation.
2. **κ_sym is global** (library-prep property, not chromosome-specific).
3. **gDNA FL is global** (library-prep property).
4. **No scipy dependency** — uses normal approximation for scoring and
   MoM for κ estimation (numpy-only).
5. **Strand-specificity-informed RNA alternative**: uses the known SS
   from spliced fragment training to define `p_rna` per gene strand.
6. **Must support unstranded data** (SS ≈ 0.5) via density signal.
