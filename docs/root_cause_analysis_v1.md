# Root Cause Analysis: gDNA Absorbing nRNA in the T+N+2 Model

*Date: March 13, 2026*

## Executive Summary

The 1944-run baseline sweep with `locus_config_simple.yaml` reveals that when nRNA
and gDNA coexist, the EM solver systematically overestimates gDNA (+107.9% mean
|error|) and underestimates nRNA (−65.3% mean signed error).  nRNA is almost
completely absorbed by gDNA in many scenarios.  Pure-mRNA and pure-nRNA scenarios
(without gDNA) produce near-perfect results, confirming the issue is specifically
in the nRNA/gDNA competition.

This document traces the root causes through the code and proposes hypotheses for
fixes.

---

## Issue 1: `nrna_init` is Dead Code in the Production Path

### Evidence

`compute_nrna_init()` in `locus.py` computes a per-nRNA warm-start value using
the strand-excess formula:

```
nRNA_int = (sense_int − anti_int) / (2SS − 1)
```

This is the correct analytical estimator for nRNA fragment count.  The result is
stored as `estimator.nrna_init`, passed through `estimator.run_batch_locus_em()`,
and reaches the C++ `batch_locus_em()` as `nrna_init_arr`.  In C++:

```cpp
const double* nri = nrna_init_arr.data();    // line 1717
```

This pointer is then passed to `extract_locus_sub_problem()` as `all_nrna_init`:

```cpp
extract_locus_sub_problem(
    sub, t_arr, n_t, u_arr, n_u, gdna_init,
    ...,
    unambig_row_sums.data(), nri,  // ← passed but never read
    ...);
```

**The pointer `all_nrna_init` is declared as a parameter at line 1206 but is NEVER
dereferenced anywhere in the function body.**  The computed nRNA warm-start evidence
is completely discarded.

### Impact

The EM solver's nRNA component receives no informed initialization.  Instead, it
relies entirely on the OVR coverage-weighted warm start (Issue 3), which distributes
mass uniformly among candidates and has no strand awareness.

---

## Issue 2: `gdna_init` is a Binary Gate, Not an Initializer

### Evidence

`compute_eb_gdna_priors()` in `locus.py` implements a sophisticated 3-level
hierarchical Empirical Bayes estimator (global → reference → locus) to compute
`gdna_init = shrunk_density × exonic_bp`.  Despite the careful estimation, the
value is used only as a binary test:

```cpp
// em_solver.cpp, extract_locus_sub_problem, line 1444
if (gdna_init == 0.0) {
    sub.prior[sub.gdna_pos_idx] = 0.0;
    sub.prior[sub.gdna_neg_idx] = 0.0;
}
```

When `gdna_init == 0`, gDNA components are eliminated permanently:
- `prior[g_pos] = prior[g_neg] = 0` → `eligible[g_pos] = eligible[g_neg] = 0`
- The OVR warm start gives them zero weight, zero OVR prior, zero theta_init.
- In the M-step, `alpha = unambig + em_totals + prior = 0 + 0 + 0` → dead forever.

When `gdna_init > 0`, the actual value is never used.  The component gets
`prior = EM_PRIOR_EPSILON = 1e-10` (immediately overwritten by OVR), and the
initialization comes from the OVR coverage share.

The name `gdna_init` is misleading — it is a **gate**, not an initializer.

### Impact

The carefully computed Empirical Bayes gDNA count estimate (which could provide a
good starting point for gDNA theta) is reduced to a binary signal.  The EM starts
from the uninformed OVR warm start instead.

---

## Issue 3: OVR Warm Start Ignores Strand Information

### Evidence

The OVR warm start (`compute_ovr_prior_and_warm_start()` in em_solver.cpp) drives
the actual EM initialization.  For each fragment, it distributes coverage weight
among candidate components:

```cpp
for (int j = 0; j < k; ++j) {
    double w = wt[i * k + j] * eligible[cidx[j]];
    row_sum += w;
}
// ...
for (int j = 0; j < k; ++j) {
    double share = (wt[i * k + j] * eligible[cidx[j]]) * inv_row_sum;
    theta_init_out[cidx[j]] += share;
}
```

It uses `wt_flat` (coverage weights), **not** `ll_flat` (log-likelihoods).
Coverage weights are based on fragment position within a transcript (edge
correction), not on strand or fragment-length compatibility.

For a sense, intronic fragment with candidates {nRNA, g_pos}, both have
`coverage_wt ≈ 1.0`.  The warm start gives each **50% share**, regardless of
strand specificity.

**Worse: for an antisense fragment at SS=1.0,** nRNA has log-likelihood = −∞
(impossible under the model), but the warm start still gives nRNA a coverage share
because `coverage_wt > 0`.  The warm start gives nRNA credit for fragments the
model says nRNA cannot produce.

### Impact

For a test case with 10,000 nRNA fragments + 720 gDNA fragments in the locus:
- ~10,000 sense intronic fragments → nRNA gets 5,000, g_pos gets 5,000
- ~360 antisense fragments → nRNA gets ~180, g_neg gets ~180

Initial g_pos theta ≈ 5000/10720 ≈ 0.47 (true: 0.034).
Initial nRNA theta ≈ 5000/10720 ≈ 0.47 (true: 0.93).

The initialization is catastrophically wrong for gDNA — 14× too high.

---

## Issue 4: The EM Cannot Escape the Bad Initialization

### The Identifiability Problem

For sense, unspliced fragments with identical FL distributions, nRNA and gDNA
have nearly identical likelihoods:

| Component | Log-likelihood (sense frag) |
|-----------|---------------------------|
| nRNA      | log(SS) + FL + bias_correction |
| g_pos     | 0.0 + FL + bias_correction |

The strand term difference is `log(SS) − 0 = log(SS)`.  At SS=0.95, this is
only −0.05 nats per fragment.  At SS=0.5, it is −0.693 (equal to gDNA).  At
SS=1.0, it is 0 (identical).

### The EM Dynamics

At the MAP-EM fixed point (ignoring coupling and prior):

```
r_new = theta_gpos_new / theta_nrna_new
      = (em_gpos / em_nrna)
      = theta_gpos / (theta_nrna × SS)
      = r / SS
```

This means:
- **SS = 1.0:** `r_new = r` — any ratio is a fixed point.  The EM stays at
  whatever the initialization gives (50/50).
- **SS < 1.0:** `r_new = r/SS > r` — gDNA ratio **grows** each iteration.  The
  EM systematically shifts mass from nRNA to gDNA.  At SS=0.95, gDNA grows by
  5.3% per iteration.

In both cases, the EM cannot escape the bad initialization.  The strand signal
is too weak per-fragment to overcome the 50/50 starting point, and at SS < 1.0
it actually makes things worse by giving gDNA a slight per-fragment advantage.

### The Coupling is Too Weak

The Beta(κ/2, κ/2) strand symmetry coupling:

```cpp
double phi = (R_pos + kappa/2 - 1) / (R_g + kappa - 2);
```

With default κ = 6.0 and R_pos ≈ 5000, R_neg ≈ 360:

```
phi = (5000 + 2) / (5360 + 4) = 0.933
```

The coupling barely moves phi from 0.933 toward 0.5.  Even κ = 10,000 only gives
phi = 0.655.  The coupling formula is a Bayesian update — when the "data" (R_pos)
is thousands of fragments, a prior with κ = 6 has negligible influence.

---

## Issue 5: Simulation Zero-Weight Fallback (separate)

When all transcript weights are 0, `compute_rna_split()` returns `(n_rna, 0)` and
`_compute_probs()` distributes mRNA fragments uniformly.  This makes the
zero-baseline tests invalid as they contain phantom RNA fragments.  This is a sim
bug, not a Rigel quantification bug.

---

## Issue 6: Unexpressed Transcript Leakage (TD, TE)

TD (single-exon, +strand, always weight=0) shows up to 3332 observed fragments.
TE (2-exon, −strand, always weight=0) shows up to 3434 observed.

For TD (single-exon), the nRNA prior is correctly zeroed (single-exon check).
But the mRNA component for TD remains eligible.  Any fragment landing in the
TD region (from gDNA or nRNA overshoot) has mRNA_TD as a candidate.  With the
same FL distributions, mRNA_TD is a valid candidate and absorbs fragments.

TE is on the opposite strand.  Antisense fragments from gDNA in the TE region
are candidates for mRNA_TE.  At low SS, this leakage is substantial.

---

## Hypotheses and Proposed Fixes

### Hypothesis A: Use Analytical Estimates as EM Initialization (Primary Fix)

**Problem:** The OVR warm start is strand-blind and creates ~50/50 splits.

**Fix:** Use `nrna_init` and `gdna_init` as actual initialization values, not just
gates.  Specifically, in `extract_locus_sub_problem()`:

1. Actually read `all_nrna_init[gnrna]` and place it into the nRNA component's
   `unambig_totals` or a new `init_totals` vector.
2. Use `gdna_init` numerically to set initial g_pos and g_neg mass (split 50/50).
3. These analytical estimates feed into `theta_init` (along with the OVR share),
   giving the EM a much better starting point.

**Expected outcome:** The EM starts near the true solution and converges correctly.
The analytical estimates are derived from strand information (which the OVR ignores),
providing the missing signal.

**Risk:** If the analytical estimates are wrong (e.g., at low SS where strand signal
is weak), the EM starts from a bad point.  Mitigate by blending analytical + OVR.

### Hypothesis B: Replace OVR Warm Start with Likelihood-Aware Initialization

**Problem:** The OVR uses `wt_flat` (coverage weights) not `ll_flat` (log-likelihoods).

**Fix:** In `compute_ovr_prior_and_warm_start()`, use log-likelihoods to weight the
coverage shares.  For each fragment, compute the posterior under uniform theta and
use that as the share:

```cpp
// Instead of: share = wt * eligible / row_sum
// Use: share = exp(ll + some_constant) * eligible / row_sum
```

This would correctly give nRNA near-zero share for antisense fragments at SS=1.0,
and proportional shares based on actual compatibility.

**Expected outcome:** The warm start respects strand information.  g_pos only gets
credit from fragments it can plausibly explain.

**Risk:** May need careful normalization.  The log-likelihoods can span wide ranges.

### Hypothesis C: Strengthen the Strand Symmetry Coupling

**Problem:** κ = 6.0 is negligible when R_pos is thousands of fragments.

**Fix:** Use much larger κ (e.g., 10^6), or replace the Beta coupling with a hard
constraint (`theta_gpos = theta_gneg` exactly, i.e., `phi = 0.5` always).

**Expected outcome:** g_pos is forced to match g_neg.  Since g_neg correctly captures
antisense gDNA, g_pos is calibrated to the correct level.

**Risk:** Hard coupling assumes perfect strand symmetry of gDNA.  In real data,
locus-specific strand biases might exist.  A moderate κ (100–1000) might be better.

**Limitation:** Even with perfect coupling, the EM has two fixed points (correct
and degenerate).  The initialization determines which one is reached.  So coupling
alone may not suffice — it must be combined with better initialization (A or B).

### Hypothesis D: Combined Fix (Recommended)

Combine A + B + C:

1. **Use analytical estimates for initialization** (fix the dead `nrna_init`, use
   `gdna_init` value).
2. **Make OVR likelihood-aware** so the coverage share respects strand compatibility.
3. **Increase coupling strength** to κ ≈ 100–1000 (or hard coupling) as a safety net.

This is the most robust approach: the analytical estimates and likelihood-aware OVR
provide a good starting point, the coupling provides a constraint, and the EM refines.

### Hypothesis E: Single gDNA Component (Model Simplification)

**Problem:** The T+N+2 model splits gDNA into g_pos/g_neg, creating the coupling
requirement.

**Fix:** Use a single gDNA component with an explicit strand probability of 0.5 in
its log-likelihood.  Remove the architectural strand routing entirely.  Every fragment
has a single gDNA candidate with `LL = log(0.5) + FL_gdna(footprint)`.

**Expected outcome:** No coupling needed.  The strand probability is baked into the
likelihood.  For sense fragments: nRNA gets log(SS) + FL, gDNA gets log(0.5) + FL.
At SS > 0.5, nRNA wins for sense; gDNA wins for antisense.

**Risk:** Loses the explicit strand decomposition that the T+N+2 model was designed
for.  May reduce accuracy of strand-specific gDNA diagnostics.

**Note:** This is essentially reverting to a single-gDNA model but with correct
strand accounting.  It was rejected previously in favor of T+N+2 for diagnostic
reasons, but might be reconsidered if the T+N+2 coupling proves unworkable.

---

## Supplementary: The Math (Why gDNA Absorbs nRNA)

### E-step for a sense unspliced fragment

```
P(nRNA | f) ∝ θ_nrna × exp(log(SS) + FL)
P(g_pos | f) ∝ θ_gpos × exp(0 + FL)
```

Likelihood ratio: `P(g_pos) / P(nRNA) = (θ_gpos / θ_nrna) × (1/SS)`

### M-step (MAP-EM, no coupling)

```
θ_gpos_new ∝ Σ P(g_pos | f_i)   for sense unspliced fragments
θ_nrna_new ∝ Σ P(nRNA | f_i)    for all nRNA-compatible fragments
```

At the fixed point:

```
r = θ_gpos / θ_nrna
r_new = r / SS
```

- SS = 1.0: r_new = r → any ratio is stable (non-identifiable)
- SS = 0.95: r_new = 1.053 × r → gDNA grows 5.3% per iteration
- SS = 0.75: r_new = 1.333 × r → gDNA grows 33% per iteration
- SS = 0.50: r_new = 2.0 × r → gDNA doubles each iteration

The T+N+2 model, by removing log(0.5) from gDNA to avoid double-counting the
strand routing, creates an intrinsic bias: gDNA has a 1/SS advantage per sense
fragment that the weak coupling cannot offset.

### With hard coupling (θ_gpos = θ_gneg)

```
g_pos = (em_gpos + N_anti) / 2    [averaged with antisense count]
```

Solving for the fixed point at SS=1.0 gives two solutions:
- g = 360 → correct (total gDNA = 720)
- g = 5360 → degenerate (almost everything is gDNA)

The OVR warm start puts the initial g near 5000, which is closer to the degenerate
solution.  The EM converges to the wrong fixed point.

This proves that **even hard coupling is insufficient without a better initialization**.
The fix MUST address initialization (Hypotheses A or B) in addition to coupling.
