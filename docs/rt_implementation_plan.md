# Implementation Plan: Making R_t (Unspliced-to-Spliced Ratio) Influence gDNA Estimation

## 1. Diagnosis: Why the Previous Implementation Was Inert

The geometric splicing expectation was injected into `gdna_init`, but `gdna_init`
is a **boolean gate** — the magnitude is discarded at two independent locations:

**Python** (`locus.py:388-390`):
```python
if gdna_init == 0.0:
    prior[gdna_idx] = 0.0
```

**C++** (`em_solver.cpp:1646-1648`):
```cpp
if (gdna_init == 0.0) {
    sub.prior[sub.gdna_idx] = 0.0;
}
```

In both cases, `gdna_init > 0` → `prior[gdna] = EM_PRIOR_EPSILON` (a tiny constant).
The actual prior mass comes from the **OVR (One Virtual Read) computation**
(`em_solver.cpp:808-862`), which overwrites the initial prior entirely.

## 2. How gDNA Prior Currently Works (Full Signal Chain)

```
gdna_init (from EB shrinkage)
    │
    ▼
Boolean gate: gdna_init > 0 → eligible[gdna] = 1.0
                gdna_init = 0 → eligible[gdna] = 0.0
    │
    ▼
OVR prior: coverage_totals[gdna] = Σ (share of each unspliced frag → gDNA)
           prior[gdna] = alpha_flat + gamma × coverage_totals[gdna] / n_ambiguous
    │
    ▼
E-step: posterior ∝ exp(log_lik[frag,gdna]) × theta[gdna] / eff_len[gdna]
    │
    ▼
M-step: theta_new[gdna] = unambig[gdna] + em_totals[gdna] + prior[gdna]
        + strand symmetry penalty (TEP)
```

Key: The OVR prior for gDNA scales with `coverage_totals[gdna]` — the sum of
coverage shares from ALL unspliced fragments that have gDNA as a candidate.
This includes **phantom mRNA unspliced fragments**, inflating the gDNA prior.

## 3. Per-Fragment Strand Discrimination (Already in E-step)

Strand evidence is already baked into per-fragment log-likelihoods:

| Fragment type     | mRNA log_lik         | gDNA log_lik  | mRNA advantage |
|-------------------|----------------------|---------------|----------------|
| Sense, SS=0.95    | log(0.95) + log_fl   | log(0.50) + … | 1.9×           |
| Antisense, SS=0.95| log(0.05) + log_fl   | log(0.50) + … | 0.1× (gDNA 10× better) |
| Sense, SS=0.80    | log(0.80) + log_fl   | log(0.50) + … | 1.6×           |
| Sense, SS=0.65    | log(0.65) + log_fl   | log(0.50) + … | 1.3×           |

At low SS, the per-fragment discrimination is weak. The M-step strand symmetry
penalty (TEP) provides a second line of defense at the aggregate level.

## 4. Why Splice-Class Probability Is a NO-OP

Adding `log(R_t / (1 + R_t))` to mRNA log-likelihoods for unspliced fragments
seems like it should help. But the per-fragment effective length correction
(BiasProfile) already accounts for this:

```
P(frag | transcript t, unspliced) × P(unspliced | t)
    = (1/L_total) × (R_t × L_spliced / L_total)
    = R_t × L_spliced / L_total²
```

When combined with the effective length normalization in `log_weights`, the
R_t factor cancels out. The E-step already implicitly handles the splice-class
probability through the geometric relationship between L_total and L_spliced.

## 5. Viable Intervention Points

### Option A: Promote gdna_init to a Magnitude-Sensitive Prior ★ RECOMMENDED

**Concept**: Replace the boolean gate with a scaled prior, analogous to how nRNA
uses informative Beta priors (`nrna_frac_alpha/beta`).

**Current behavior**:
```
gdna_init > 0  →  prior[gdna] = alpha_flat + OVR_contribution
gdna_init == 0 →  prior[gdna] = 0 (disabled)
```

**Proposed behavior**:
```
prior[gdna] = gdna_init_scaled + OVR_contribution
```

where `gdna_init_scaled` encodes the EB-estimated gDNA count, adjusted by R_t.

**Where R_t enters**: In `compute_eb_gdna_priors()` (`locus.py:688-857`), the
strand+density hybrid estimator computes gDNA density. Currently, the strand-based
estimate treats ALL unspliced fragments as potential gDNA signal. R_t predicts
how many are phantom mRNA:

```
observed_unspliced = true_gDNA + phantom_mRNA
phantom_mRNA ≈ R_t × spliced_count
adjusted_unspliced = max(observed_unspliced - R_t × spliced_count, 0)
```

This adjustment flows through the 3-tier EB shrinkage (global → ref → locus)
and produces a more accurate `gdna_init` magnitude.

**C++ changes required**:
1. `em_solver.cpp:1643-1649`: Use `gdna_init` as actual prior mass:
   ```cpp
   // Replace boolean gate with scaled prior
   sub.prior[sub.gdna_idx] = std::max(gdna_init, 0.0);
   // If gdna_init == 0, eligible will be 0 (same as before)
   ```
2. `compute_ovr_prior_and_warm_start()`: Add `gdna_init` to the OVR prior
   instead of `alpha_flat` for the gDNA component:
   ```cpp
   prior_out[gdna_idx] = (eligible[gdna_idx] > 0.0)
       ? (gdna_init_prior + ovr_gdna)
       : 0.0;
   ```

**Advantages**:
- EB machinery already computes meaningful magnitudes (just currently wasted)
- nRNA already uses magnitude-sensitive priors — consistent architecture
- R_t integrates naturally into the EB density estimation
- No new per-fragment arrays or EC bloat

**Risks**:
- OVR and gdna_init interact — calibration needed to avoid double-counting
  (the strand signal contributes to both)
- gdna_init from EB shrinkage may be noisy at low-coverage loci
- Need to define the right scaling (raw count? pseudo-count? fraction?)

### Option B: Scale OVR Prior for gDNA Using R_t

**Concept**: Pass per-locus R_t into the C++ EM. In `compute_ovr_prior_and_warm_start()`,
scale down gDNA's coverage contribution by the expected phantom fraction.

**Implementation**:
```cpp
// After computing coverage_totals[gdna_idx]:
double phantom_fraction = R_t_locus / (1.0 + R_t_locus);
coverage_totals[gdna_idx] *= (1.0 - phantom_fraction);
// = coverage_totals[gdna_idx] / (1.0 + R_t_locus)
```

**Where R_t enters**: New field `gdna_ovr_scale` in `LocusEMInput`, computed in
Python as `1.0 / (1.0 + R_t_locus)`. Passed through to C++.

**Advantages**:
- Minimal C++ change (one multiplication)
- Directly addresses the root cause (OVR inflation from phantom mRNA)
- R_t naturally vanishes for single-exon genes (R_t → ∞, scale → 0)

**Risks**:
- R_t quality degrades at high gDNA contamination (fragment length model is
  contaminated by gDNA fragments — R_t drops from 3.6 to 1.5 at extreme gDNA)
- Single-exon genes have undefined R_t (no spliced fragments to calibrate)
- Doesn't address the deeper issue of gdna_init being wasted information

### Option C: Ablation-First (No Code Change)

**Concept**: Before adding complexity, determine whether the M-step strand
symmetry penalty is even necessary. Run the benchmark suite with
`strand_symmetry_kappa = 2.0` (disables TEP entirely) and compare.

**Rationale**: The E-step already has per-fragment strand discrimination.
The TEP is a second layer. If the E-step alone is sufficient, we can remove
the TEP and simplify the model before adding R_t.

**Implementation**: Config change only — `EMConfig.strand_symmetry_kappa = 2.0`.

## 6. Recommended Approach: Option A with R_t-Adjusted EB

### Step-by-step plan

#### Phase 1: R_t Computation (Python only)

1. **Compute per-transcript R_t** in `estimator.py`:
   - Use `frag_length_model.compute_all_transcript_eff_lens()` with **spliced-only**
     lengths to get `L_eff_spliced[t]`
   - R_t = `(L_total[t] - L_eff_spliced[t]) / L_eff_spliced[t]`
   - Important: Use the spliced-only fragment length model (from spliced fragments
     only) if available, to avoid gDNA contamination of the length distribution

2. **Aggregate R_t to locus level**:
   - Weighted average: `R_locus = Σ(R_t × theta_t) / Σ(theta_t)` where theta_t
     is from unambiguous counts or a first-pass estimate
   - Or simpler: `R_locus = Σ(L_total - L_spliced) / Σ(L_spliced)` across
     transcripts at the locus

#### Phase 2: Adjust EB gDNA Density Estimation (Python)

3. **Modify `compute_gdna_density_hybrid()`** (`locus.py:538-601`):
   - Currently computes gDNA density from total sense + antisense unspliced
   - Subtract predicted phantom mRNA:
     ```python
     # Before: uses raw sense/anti counts
     # After: subtract predicted mRNA contribution
     predicted_phantom = R_locus * spliced_count
     adjusted_anti = max(anti - predicted_phantom * (1 - SS), 0)
     adjusted_sense = max(sense - predicted_phantom * SS, 0)
     ```
   - Or equivalently, adjust the total unspliced count:
     ```python
     adjusted_unspliced = max(total_unspliced - R_locus * spliced_count, 0)
     ```

4. **This flows through existing EB shrinkage** — no changes needed to the
   3-tier hierarchy (global → ref → locus). The shrinkage automatically
   handles noisy per-locus estimates.

#### Phase 3: Promote gdna_init to Magnitude (C++)

5. **Modify em_solver.cpp** to use `gdna_init` as prior mass:
   - Replace boolean gate (`gdna_init == 0.0`) with:
     ```cpp
     sub.prior[sub.gdna_idx] = gdna_init;  // magnitude, not gate
     // eligible remains: (prior > 0) ? 1.0 : 0.0
     ```
   - In `compute_ovr_prior_and_warm_start()`, use `gdna_init` instead of
     `alpha_flat` for the gDNA component:
     ```cpp
     double base_prior = (i == gdna_idx) ? gdna_init : alpha_flat;
     prior_out[i] = (eligible[i] > 0.0) ? (base_prior + ovr_i) : 0.0;
     ```

6. **Calibrate gdna_init scale**: The EB shrinkage produces
   `gdna_init = shrunk_density × locus_bp`. This is in units of expected
   fragment count. The OVR `alpha_flat` is typically 1e-10 (EM_PRIOR_EPSILON).
   We need gdna_init to be on the right scale relative to OVR contributions.
   Options:
   - Use raw: `gdna_init` as-is (units: expected fragment count)
   - Scale to pseudo-count: `gdna_init_scaled = gdna_init × pseudo_weight / median_gdna_init`
   - Use as fraction: `gdna_init / total_fragments × scaling_factor`

   The right calibration depends on empirical testing.

#### Phase 4: Benchmark

7. **Run benchmark suite** comparing:
   - Baseline (current code, boolean gate)
   - Option A (magnitude prior, no R_t adjustment)
   - Option A + R_t adjustment
   - Across SS values (0.65, 0.80, 0.90, 0.95, 0.99) and gDNA levels

## 7. Key Caution: R_t Contamination

The global fragment length model mixes mRNA and gDNA fragments. At high gDNA:
- R_t drops from ~3.6 (clean) to ~1.5 (extreme gDNA = 1000 frags/bp)
- This causes R_t to **underpredict** phantom mRNA → insufficient gDNA correction

**Mitigation**: Compute R_t using a **spliced-only fragment length model**:
- Build FLD from only `is_spliced == True` fragments
- These are guaranteed mRNA (gDNA cannot produce spliced fragments)
- Use this clean FLD for `compute_all_transcript_eff_lens()`

Check whether `FragLengthModel` can be built from a subset of fragments
(it uses a histogram, so filtering by `is_spliced` before building should work).

## 8. Single-Exon Gene Handling

For single-exon genes: `L_spliced = L_total`, so `R_t = 0`. This means:
- No phantom mRNA adjustment (correct — single-exon genes produce no
  distinguishable unspliced pattern)
- All unspliced fragments remain in the gDNA pool (correct behavior)
- gdna_init magnitude from EB shrinkage still provides useful prior information

## 9. Files to Modify

| File | Change |
|------|--------|
| `src/rigel/estimator.py` | Add `compute_per_transcript_Rt()` method |
| `src/rigel/locus.py` | Modify `compute_gdna_density_hybrid()` to accept R_t adjustment; modify `compute_eb_gdna_priors()` to pass R_t |
| `src/rigel/locus.py` | Modify `build_locus_em_data()` — gdna_init already passed through, just ensure magnitude preserved |
| `src/rigel/native/em_solver.cpp` | Replace boolean gate with magnitude-sensitive prior; modify OVR to use gdna_init as base |
| `src/rigel/pipeline.py` | Wire up R_t computation; optionally build spliced-only FLD |
| `src/rigel/frag_length_model.py` | Possibly add method to build from spliced-only subset |

## 10. What NOT to Change

- **E-step log_liks**: Strand evidence is already per-fragment. Adding splice-class
  probability is a NO-OP (already implicit in L_total normalization).
- **Effective lengths**: Currently all 1.0 because BiasProfile handles per-fragment
  correction. Changing component-level eff_len would double-count.
- **M-step TEP**: Keep the strand symmetry penalty for now. Ablation testing
  (Option C) can be done separately.
