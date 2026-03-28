# VBEM Convergence Issue on Mega-Loci

**Status:** Deferred — requires simulation with ground truth to validate changes
**Date:** March 28, 2026

## Problem

VBEM fails to converge on the large BAM's mega-locus (315K transcripts, 5.8M components), hitting the max iteration limit (1000 EM steps = 333 SQUAREM iterations). MAP-EM converges in 123 SQUAREM iterations on the same locus.

This makes VBEM **2.23× slower** on the large mega-locus (634.6s vs 285.1s SQUAREM time), despite being 18% faster per SQUAREM iteration.

## Root Cause

The convergence criterion is the L1 norm of normalized theta change:

```
delta = Σᵢ |theta_new[i] - theta_old[i]|
```

where `theta[i] = alpha[i] / Σ alpha` and `convergence_delta = 1e-6`.

With 5,777,023 components, reaching `delta < 1e-6` requires each component to change by < 1.7×10⁻¹³ on average — near machine epsilon (2.2×10⁻¹⁶).

VBEM's Dirichlet posterior distributes mass more smoothly than MAP-EM's point estimates. Small components that MAP-EM drives to exactly zero persist at tiny positive values in VBEM, creating persistent oscillations that prevent the L1 norm from dropping below the threshold.

## Evidence

### Large BAM Mega-Locus

| Metric | VBEM | MAP-EM |
|--------|------|--------|
| Components (n_units) | 5,777,023 | 5,777,023 |
| SQUAREM iterations | **333 (MAX)** | 123 |
| SQUAREM time | 634.6s | 285.1s |
| Per-iteration cost | 1.91s | 2.32s |
| Converged? | **No** | Yes |

### Small BAM Mega-Locus (179K components)

| Metric | VBEM | MAP-EM |
|--------|------|--------|
| Components (n_units) | 179,704 | 179,704 |
| SQUAREM iterations | 74 | 121 |
| Converged? | Yes | Yes |

VBEM converges fine on the smaller mega-locus — the problem is specific to very high-dimensional problems.

## Proposed Solutions (to evaluate with ground truth)

### Option A: Scale convergence delta by dimensionality

Use `convergence_delta * sqrt(n_components)` or `convergence_delta * n_components`. This accounts for the fact that the L1 norm naturally scales with the number of components.

**Code location:** `em_solver.cpp`, `run_squarem()`, lines ~929-948 (VBEM branch) and ~1028-1036 (MAP branch).

### Option B: Use L∞ (max) norm instead of L1

Replace `Σ |Δθᵢ|` with `max |Δθᵢ|`. This is dimensionality-invariant and checks that no single component is still changing significantly.

### Option C: Use relative ELBO change

Track the ELBO (evidence lower bound) and stop when relative change `|ELBO_new - ELBO_old| / |ELBO_old| < delta`. More principled for variational inference but requires computing the ELBO (additional cost per iteration).

### Option D: Increase max_iterations for mega-loci only

Detect mega-loci and allow more iterations. Not recommended — worsens performance without guaranteed convergence.

## Validation Plan

Before implementing any fix:
1. Run synthetic simulations with known ground truth abundances
2. Compare VBEM vs MAP-EM accuracy (relative error, correlation)  
3. Implement convergence fix
4. Verify accuracy is maintained or improved
5. Verify performance improvement on real BAMs
