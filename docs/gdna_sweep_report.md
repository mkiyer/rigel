# gDNA Siphoning Fix: Simulation Sweep Report

## Executive Summary

Three simulation sweeps (197 total runs) evaluated the Phase 1 (density-based gDNA
initialization) and Phase 2 (soft Beta symmetry penalty) changes. **Phase 1 is the
primary improvement.** Phase 2 (the symmetry penalty) improves mRNA accuracy at
moderate strength (kappa=3) but worsens nRNA estimation by transferring suppressed
gDNA mass into nRNA. **The current default kappa=6.0 is too aggressive and should
be reduced to kappa=2.0 (disabled) or kappa=3.0 (mild).**

---

## 1. Have We Improved?

**Yes — Phase 1 (density-based gDNA initialization) is the core improvement.**

With kappa=2.0 (penalty disabled, Phase 1 only), the pipeline already achieves
reasonable gDNA/mRNA separation without the pathological siphoning that occurred
when gDNA was initialized proportional to mRNA expression. The symmetry penalty
(Phase 2) provides marginal additional mRNA accuracy at the cost of worse nRNA/gDNA
component estimates.

## 2. Quantify Improvement

### Grid Search — Single Hard Scenario

(TA1=100, TA4=3, NTA=100, gDNA=100, SS=1.0)

| kappa | pseudo | mRNA err | gDNA diff | nRNA diff | Composite |
|-------|--------|----------|-----------|-----------|-----------|
| 2.0   | any    | 7.6%     | 116       | 80        | **231**   |
| 3.0   | 50     | 7.5%     | 45        | 80        | **160**   |
| 3.0   | 1      | 7.4%     | 62        | 97        | 194       |
| 6.0   | 10     | 0.6%     | 296       | 293       | 591       |
| 6.0   | 50     | 1.3%     | 266       | 272       | 545       |
| 15.0  | 50     | 11.6%    | 508       | 454       | 1016      |
| 20.0+ | any    | >100%    | >1000     | >500      | >1900     |

On this single scenario, kappa=3/pseudo=50 achieves the best composite score (160
vs 231 baseline). kappa=6 minimizes mRNA error (0.6%) but at 2.5× worse composite.

### Robustness Sweep — 4 gDNA levels × 4 SS × 2 NTA

(All conditions with gDNA > 0, 3 kappa values)

| kappa | Mean mRNA err | Max mRNA err | Mean gDNA diff | Mean nRNA diff | Composite |
|-------|---------------|--------------|----------------|----------------|-----------|
| 2.0 (disabled) | 18.6% | 39.1% | 142 | **85** | **285** |
| 3.0 (mild)     | **14.8%** | **32.7%** | **125** | 144 | 315 |
| 6.0 (default)  | 16.1% | 62.9% | 451 | 420 | 911 |

**Key observation:** kappa=3 improves mRNA error by 20% relative (18.6%→14.8%) and
actually improves gDNA separation (142→125). However, it worsens nRNA estimation
(85→144), dragging composite from 285→315. kappa=6 is categorically worse on every
aggregate metric except mRNA-only accuracy in narrow conditions.

### By gDNA Level

| gDNA level | k=2 mRNA err | k=3 mRNA err | k=2 composite | k=3 composite | k=6 composite |
|------------|-------------|-------------|---------------|---------------|---------------|
| 50         | 8.3%        | 6.7%        | 310           | **267**       | 697           |
| 100        | 14.0%       | 12.0%       | **259**       | 421           | 1103          |
| 200        | 17.4%       | 12.7%       | **244**       | 346           | 834           |
| 500        | 34.6%       | 28.0%       | 328           | **225**       | 1009          |

kappa=3 wins at gDNA=50 and gDNA=500; kappa=2 wins at gDNA=100 and gDNA=200. Neither
is universally better. kappa=6 loses everywhere.

### By Strand Specificity

| SS   | k=2 mRNA err | k=3 mRNA err | k=2 composite | k=3 composite | k=6 composite |
|------|-------------|-------------|---------------|---------------|---------------|
| 0.80 | 22.3%       | 18.1%       | 502           | **447**       | 1339          |
| 0.90 | 19.3%       | 15.5%       | **263**       | 400           | 1085          |
| 0.95 | 18.1%       | 14.7%       | **172**       | 280           | 770           |
| 1.00 | 14.7%       | 11.0%       | 203           | **133**       | 449           |

kappa=3 wins at SS=0.8 (low specificity) and SS=1.0 (perfect specificity). kappa=2
wins at SS=0.9 and SS=0.95. kappa=3 always has better mRNA error.

## 3. Hyperparameter Recommendations

### strand_symmetry_kappa

**Recommendation: Change default from 6.0 to 3.0**, with `pseudo=50.0`.

Rationale:
- kappa=6.0 is unambiguously too aggressive: 3.2× worse composite than disabled,
  62.9% max mRNA error (vs 39.1% disabled), catastrophic gDNA underestimation
- kappa=3.0 improves mRNA accuracy vs disabled (14.8% vs 18.6% mean, 32.7% vs
  39.1% max) at moderate cost to nRNA estimation
- kappa=2.0 (disabled) is the safest choice if nRNA accuracy matters as much as
  mRNA accuracy
- pseudo=50.0 provides meaningful regularization at kappa=3 but has minimal effect
  at kappa=2 (which is invariant to pseudo)

**Alternative: Default to kappa=2.0** if the conservative "do no harm" principle is
preferred. Phase 1 alone (density-based init) provides the primary gDNA siphoning
fix without the nRNA overestimation side effect.

### strand_symmetry_pseudo

At kappa=3, higher pseudo is modestly better (composite 194→160 as pseudo goes
1→50). At kappa=2.0, pseudo is irrelevant. **Recommend pseudo=50.0** if kappa>2.

## 4. Remaining Issues and Opportunities

### Fundamental Limitation of the Symmetry Penalty

The penalty suppresses ALL gDNA when strand asymmetry is detected, including
legitimate symmetric gDNA. The freed mass transfers to nRNA, causing systematic
nRNA overestimation. This is why higher kappa always increases nRNA error.

**Root cause:** When gDNA siphons mRNA reads, those absorbed sense-only reads bias
the gDNA strand ratio asymmetric (triggering the penalty — correct). But the penalty
then suppresses the entire gDNA component, not just the siphoned portion.

### High gDNA Contamination (500+ reads)

At gDNA=500 (5× the dominant transcript), even the best configuration shows 28%
mRNA error. This is fundamentally challenging because gDNA fragments that overlap
exons are nearly indistinguishable from mRNA fragments without strand information.

### Low Strand Specificity

SS=0.8 produces 18-22% mRNA error regardless of penalty settings. Low strand
specificity makes gDNA/mRNA separation inherently harder because there's less
information in the antisense signal to distinguish gDNA from mRNA.

### Potential Future Improvements

1. **Targeted penalty**: Instead of penalizing all gDNA, penalize only the
   asymmetric excess. This would require decomposing gDNA into symmetric and
   asymmetric components during the M-step.
2. **Fragment-level strand features**: Use per-fragment strand evidence rather
   than aggregate strand ratios for finer-grained classification.
3. **Iterative gDNA refinement**: After initial EM convergence, re-estimate gDNA
   density from the gDNA-classified fragments and re-run, allowing the model to
   better separate gDNA from mRNA.

---

## Sweep Configurations

| Sweep | Config | Runs | Purpose |
|-------|--------|------|---------|
| Baseline | `gdna_fix_baseline_comparison.yaml` | 32 | 4 patterns × 2 NTA × 2 gDNA × 2 kappa |
| Grid | `gdna_kappa_pseudo_grid.yaml` | 45 | 9 kappa × 5 pseudo on hardest scenario |
| Robustness | `gdna_robustness_sweep.yaml` | 120 | 5 gDNA × 4 SS × 2 NTA × 3 kappa |

Results in `sweep_results/{baseline_comparison,kappa_pseudo_grid,robustness}/`.
