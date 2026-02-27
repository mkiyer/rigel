# Adaptive Prior and Sparsity-Inducing EM: Analysis and Plan

## Date: 2026-02-25

## Context

After implementing VBEM as an opt-in mode alongside MAP-EM, a comprehensive
24-condition benchmark (FGFR2 + BRCA1 × 3 gDNA × 2 nRNA × 2 SS) shows:

- **VBEM cuts FPs by 51%** (355 → 174) across ALL 24 conditions
- **Spearman improves +0.048** (0.666 → 0.714)
- **Dropout increases modestly** (20 → 30, +10 across 24 conditions)
- Two structural FPs (ENST00000369058.7, ENST00000346997.6) persist in
  both modes — these share enough exon structure that no EM variant can
  distinguish them from expressed transcripts purely from fragment evidence

The question now is: can we do better, and what are the right levers?

---

## 1. The Prior Anatomy (Current State)

Each locus EM has `2*n_t + 1` components. The prior for component `k` is:

$$\alpha_k = \underbrace{\alpha_{\text{base}}}_{\text{flat Dirichlet}} + \underbrace{\frac{w_k}{N_{\text{ambig}}}}_{\text{One Virtual Read (OVR)}}$$

where:
- $\alpha_{\text{base}} = 0.5$ (the `em_pseudocount` pipeline default),
  applied per component
- $w_k$ = sum of coverage-weighted shares across all ambiguous units for
  component $k$ (proportional to how much trapezoidal area this component
  covers across all fragments)
- $N_{\text{ambig}}$ = number of ambiguous units in the locus
- OVR sums to exactly 1.0 across the locus (parameter-free)

Components with zeroed prior (nRNA without intronic evidence, gDNA without
unspliced fragments) also get zero OVR — they cannot receive warm-start
weight. This prevents phantom self-sustaining cycles.

### What the OVR does

The OVR serves two distinct purposes:

1. **SQUAREM stabilization**: Without it, the SQUAREM extrapolation can
   overshoot in flat-likelihood regions (many shared exons, balanced
   abundances) and oscillate. The OVR acts as a gentle anchor.

2. **Coverage-weighted differentiation**: The OVR encodes our geometric
   prior belief that a transcript covering more exonic area "near" a
   fragment is a stronger candidate. This especially matters at transcript
   starts/ends where tiny unique exons provide the only differentiating
   signal.

### The FP mechanism

In MAP-EM, the M-step always adds the prior:

$$\theta_k^{(t+1)} = \frac{\text{unique}_k + \text{em}_k^{(t)} + \alpha_k}{\sum_j (\cdots)}$$

A zero-truth transcript with shared exons accumulates a non-trivial
$\text{em}_k$ from posterior mass in the E-step (because its likelihood
for shared-exon fragments is identical to expressed transcripts). Combined
with $\alpha_k > 0$, this creates a permanent positive floor — the classic
"zombie transcript" false positive.

VBEM replaces $\log(\theta_k)$ with $\psi(\alpha_k) - \psi(\sum \alpha)$
in the E-step. Since $\psi(x) \to -\infty$ as $x \to 0^+$, low-alpha
components get exponentially suppressed. This is why VBEM cuts FPs so
dramatically. But it cannot eliminate FPs from structurally indistinguishable
transcripts.

---

## 2. Design Space: Three Independent Levers

### Lever A: Prior Scale Factor

Multiply the OVR by a constant $\gamma$:

$$\alpha_k = \alpha_{\text{base}} + \gamma \cdot \frac{w_k}{N_{\text{ambig}}}$$

- $\gamma = 1.0$ (current): 1 total virtual read. Parameter-free.
- $\gamma < 1.0$: Weaker geometric prior. EM converges faster, less
  anchoring. FPs may decrease (lower floor) but SQUAREM may become
  less stable; tiny-exon differentiation weakened.
- $\gamma > 1.0$: Stronger geometric prior. Better tiny-exon
  differentiation but higher FP floor. Essentially says "I trust the
  coverage geometry more than the fragment likelihood."
- $\gamma = 0$: No OVR. Prior reverts to flat Dirichlet. We tried this
  during the failed VBEM-as-default attempt — it broke gDNA separation
  in negative control scenarios (9 test failures).

**Key insight**: The OVR is additive with $\alpha_{\text{base}}$.
With $\alpha_{\text{base}} = 0.5$ and typical FGFR2 locus having 41
mRNA + 41 nRNA + 1 gDNA = 83 components, the total flat prior is
$83 \times 0.5 = 41.5$ while the OVR adds only 1.0. **The flat base
dominates.** Tuning $\gamma$ matters less than tuning $\alpha_{\text{base}}$.

### Lever B: Base Pseudocount ($\alpha_{\text{base}}$)

This is the `em_pseudocount` parameter (currently 0.5 per component).

- $\alpha_{\text{base}} = 0.5$: Current. Each component starts with
  half a virtual read just from the flat prior.
- $\alpha_{\text{base}} = 0.01$: Much weaker. In MAP-EM, the FP floor
  drops proportionally. In VBEM, $\psi(0.01) \approx -100.6$ provides
  extreme suppression of evidence-free components.
- $\alpha_{\text{base}} = 0.0$ (improper prior): Components with no
  evidence go to exactly zero. But creates numerical issues and violates
  Dirichlet regularity.

**Interaction with VBEM**: In VBEM the base pseudocount has a qualitatively
different effect than in MAP-EM. The digamma function is:

| $\alpha$ | $\psi(\alpha)$ | Effective log-weight |
|-----------|----------------|---------------------|
| 0.01      | −100.6         | Extreme suppression |
| 0.1       | −10.4          | Strong suppression  |
| 0.5       | −1.96          | Mild suppression    |
| 1.0       | −0.58          | Minimal             |
| 10.0      | +2.25          | Neutral             |

In MAP-EM, reducing $\alpha_{\text{base}}$ gives a proportional reduction
in the FP floor. In VBEM, it gives a **superlinear** reduction because
the digamma curve is concave for small arguments. This means:

> **VBEM is much more sensitive to $\alpha_{\text{base}}$ than MAP-EM.**
> A value that works well for MAP-EM may over-suppress in VBEM, causing
> dropout. Conversely, a value conservative enough for VBEM may leave
> MAP-EM with too many FPs.

### Lever C: Post-EM Pruning

After EM convergence, identify and zero out "zombie" transcripts:

**Candidate criteria** for pruning transcript $t$:
- `unique_counts[t] == 0` (no uniquely-mapping fragments)
- Final count below threshold $\tau$

**Threshold options**:

1. **Absolute threshold** ($\tau = c$): Simple but fragile. A threshold
   of 5 might be right for 50K fragments but wrong for 500K. Doesn't
   scale.

2. **Relative to prior** ($\tau = m \cdot \alpha_k$): "Did the data move
   this component significantly beyond its prior?" If the converged count
   is less than $m$ times the prior, the evidence is insufficient. This
   is more principled — it asks whether the posterior differs meaningfully
   from the prior. A natural $m$ might be 2–5.

3. **Evidence ratio** ($\tau = \text{em}_k / (\text{em}_k + \alpha_k)$):
   What fraction of the final count came from data vs prior? If < 0.5,
   the prior dominates and we're hallucinating.

4. **VBEM posterior credible interval**: For VBEM, the converged $\alpha$
   values parameterize a Dirichlet posterior. We can compute the 95%
   upper bound of the marginal Beta($\alpha_k$, $\sum_{j \neq k} \alpha_j$)
   distribution. If this upper bound is below some small threshold (e.g.
   $10^{-4}$), prune. This is the most Bayesian approach but only works
   for VBEM.

**After pruning**: Zero out the pruned component's prior, re-run EM
(or a few final iterations) to redistribute the freed mass. This is
fast because the fixed-point is nearby.

---

## 3. MAP-EM vs VBEM: Different Optimal Configurations

Given the analysis above, I expect the optimal configurations to differ:

### MAP-EM + Pruning
- **Keep $\alpha_{\text{base}} = 0.5$**: MAP-EM needs a reasonable prior
  to stabilize the EM and feed the coverage-weighted warm start.
- **Keep OVR $\gamma = 1.0$**: The OVR is small relative to the flat
  prior; its value is in SQUAREM stability, not FP control.
- **Add evidence-ratio pruning**: After convergence, prune components
  where `(em_k + unique_k) / (em_k + unique_k + alpha_k) < 0.3` and
  `unique_k == 0`. Re-run a few EM iterations.
- **Expected outcome**: Catches the medium FPs that MAP-EM currently
  leaves at 5–30 counts. Structural FPs (100+) won't be caught because
  they have genuine EM evidence from shared exons.

### VBEM (no explicit pruning needed)
- **Lower $\alpha_{\text{base}}$ to ~0.1**: VBEM's digamma naturally
  handles sparsity. A lower base lets it suppress more aggressively.
  But not too low — we need the OVR anchor for SQUAREM and to avoid
  gDNA separation issues.
- **Keep OVR $\gamma = 1.0$**: Critical for SQUAREM stability and
  coverage-weighted differentiation at transcript boundaries.
- **No explicit pruning**: VBEM already achieves 51% FP reduction
  without pruning. Adding pruning on top is possible but less necessary.
- **Expected outcome**: Already demonstrated. Further gains from lower
  $\alpha_{\text{base}}$ would reduce the remaining ~174 FPs further.

### The Coverage Weight Question

> "We want coverage weights to help us with tiny exon transcript
> starts/ends."

The coverage weight enters the system twice:
1. **Warm start** (`theta_init`): Distributes initial mass. This is
   important for convergence speed and for nudging the EM toward the
   geometrically-correct solution early.
2. **OVR prior** (`prior += coverage_totals / n_ambig`): Persistent
   bias toward coverage-plausible components.

For MAP-EM, both contribute. For VBEM, the warm start matters less
(VBEM converges to the same fixed point regardless of initialization
for well-identified models) but the OVR prior still shapes the
posterior by adding persistent virtual counts.

**Critical question**: Can we keep the warm start but decouple OVR
influence? If we set $\gamma = 0$ (no OVR) but keep the
coverage-weighted warm start, we'd let VBEM's natural sparsity handle
FPs while still getting fast convergence. The risk is SQUAREM stability.

The safest approach is to keep OVR at $\gamma = 1.0$ but lower
$\alpha_{\text{base}}$. The OVR adds only 1 virtual read total across
~83 components — it's already tiny. The coverage signal it adds to the
prior is meaningful (components with good coverage geometry get ~0.02
extra, while poor-coverage components get ~0.001). Lowering the flat
base from 0.5 to 0.1 makes the OVR **relatively** more important,
which is exactly what we want: the coverage geometry becomes the
dominant prior signal rather than the flat floor.

---

## 4. Parameter Tuning: Why We Need Large Benchmarks

### Overfitting risk

We are currently tuning on 2 regions (82 transcripts), 24 conditions.
This is good for understanding mechanism but dangerous for parameter
selection:
- FGFR2 has one enormous structural FP that dominates MAE
- BRCA1 has a different isoform overlap pattern
- Neither is representative of the transcriptome-wide distribution

### What we need

1. **Diverse complexity**: Regions ranging from 1 transcript (trivial)
   to 40+ (FGFR2/BRCA1-level). The transcriptome median is ~3
   transcripts per gene; we need the easy cases to ensure no regression.

2. **Scale**: 50+ regions, 100K+ fragments. This gives statistical
   power to detect small regressions (a parameter change that helps
   FGFR2 but hurts 30 simpler loci is a net negative).

3. **Expression diversity**: Uniform random abundance (current) plus
   realistic skewed abundance (few high, many low, some zero) to test
   dropout sensitivity.

4. **Cross-validation**: Tune on a training set of regions, validate
   on a held-out set. Two-fold CV on 50 regions = 25 train / 25 test.

### Proposed benchmark design

| Parameter | Value |
|-----------|-------|
| Regions | 50 randomly selected non-overlapping loci |
| Complexity | Mix of 1–5 tx (30 loci), 5–15 tx (15 loci), 15+ tx (5 loci) |
| Fragments | 200K per region (10M total) |
| gDNA rates | 0.0, 0.2, 1.0 |
| nRNA rates | 0.0, 0.3 |
| SS | 1.0, 0.9 |
| Abundance | random (seed 101) + skewed (seed 202) |
| EM configs | MAP + MAP+prune + VBEM, with grid over $\alpha_{\text{base}}$ ∈ {0.01, 0.1, 0.5} |

Total: 50 regions × 12 conditions × 2 abundance modes × 3 EM × 3 prior = ~10,800 runs.
This is computationally heavy but parallelizable.

### A faster incremental approach

Before the full grid, a targeted sweep on fewer dimensions:

**Phase A: Prior scale sweep (small)**
- Current 2 regions (FGFR2 + BRCA1), 4 conditions (pristine + gDNA=r20 + nrna=r30 + both)
- $\alpha_{\text{base}}$ ∈ {0.01, 0.05, 0.1, 0.25, 0.5, 1.0}
- Both MAP and VBEM
- 2 × 4 × 6 × 2 = 96 runs (fast, ~10 min)
- **Goal**: Identify the sensitivity curve for both modes

**Phase B: Pruning threshold sweep**
- Same 2 regions, 4 conditions
- MAP-EM with $\alpha_{\text{base}} = 0.5$ (current default)
- Evidence-ratio thresholds: {0.1, 0.2, 0.3, 0.5}
- With and without re-EM after pruning
- 2 × 4 × 4 × 2 = 64 runs
- **Goal**: Find pruning threshold that catches medium FPs without dropout

**Phase C: Validation on 10-region benchmark**
- Existing 10-region pristine benchmark config (FGFR2, EGFR, FHB, BRCA1,
  HBB, HOXA, GAPDH, ELANE, BCR, HES4)
- Best configs from Phase A + B
- Full condition grid
- **Goal**: Validate that improvements generalize beyond FGFR2/BRCA1

**Phase D: Large-scale validation**
- 50 random regions
- Best 2–3 configs from Phase C
- Full condition grid
- **Goal**: Final parameter selection with statistical confidence

---

## 5. Concrete Implementation Plan

### Step 1: Parameterize `em_pseudocount` for sweep (already done)

The pipeline already accepts `em_pseudocount` as a parameter, making it
trivially sweepable via benchmark YAML `hulkrna_configs`.

### Step 2: Implement post-EM pruning (new)

Add an optional pruning pass in `run_locus_em()`:

```python
def _prune_zombies(theta, alpha, unique_totals, prior,
                   evidence_threshold=0.3):
    """Zero out components with insufficient evidence."""
    total_count = alpha  # unique + em + prior
    evidence_ratio = (total_count - prior) / np.maximum(total_count, 1e-300)
    no_unique = unique_totals == 0
    prune_mask = no_unique & (evidence_ratio < evidence_threshold)
    # Never prune gDNA component
    prune_mask[-1] = False
    return prune_mask
```

After EM convergence and before returning:
1. Compute `prune_mask`
2. Zero the prior for pruned components
3. Run 10 more EM iterations to redistribute
4. Return the re-converged theta and alpha

Gate behind a new `em_prune_threshold` parameter (default `None` = disabled).

### Step 3: Build sweep infrastructure

Create a benchmark YAML generator that produces configs for the Phase A
sweep (prior scale × EM mode) and Phase B sweep (pruning threshold).

### Step 4: Run Phase A + B sweeps

Analyze results to select candidate configurations.

### Step 5: Validate on 10 regions, then 50 regions

Confirm improvements generalize. Select final defaults.

---

## 6. Summary of Intuitions

1. **$\alpha_{\text{base}}$ is the dominant lever**, not OVR scale. The flat
   prior (0.5 per component × 83 components = 41.5 total) dwarfs the
   OVR (1.0 total). Reducing $\alpha_{\text{base}}$ has the largest
   single effect on FP reduction.

2. **MAP-EM and VBEM want different $\alpha_{\text{base}}$ values.**
   MAP-EM benefits from moderate values (0.1–0.5) with explicit pruning.
   VBEM benefits from small values (0.01–0.1) and handles sparsity
   natively via digamma.

3. **Pruning is most useful for MAP-EM.** VBEM's natural sparsity makes
   explicit pruning largely redundant, though it could still help with
   the structural FPs where EM evidence is genuinely present.

4. **Coverage-weighted OVR should be preserved.** It's small in absolute
   magnitude but provides SQUAREM stability and boundary differentiation.
   Lowering $\alpha_{\text{base}}$ makes it relatively more influential,
   which is desirable.

5. **Tuning on 2 regions is insufficient for parameter selection.** We can
   identify promising regions of the parameter space with the current
   micro-benchmark, but final values must be validated on diverse loci
   to avoid overfitting.

6. **The incremental phased approach (A→B→C→D) reduces risk.** Each phase
   is a few hours of compute at most, and each provides a clear
   go/no-go signal before investing in the next phase.
