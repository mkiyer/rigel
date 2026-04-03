# VBEM Convergence Fix: Variational Zero-Forcing

**Status:** Ready for implementation  
**Date:** April 2, 2026

## Problem

VBEM fails to converge on mega-loci because SQUAREM's quadratic extrapolation is incompatible with the digamma function's infinite gradient at the $\alpha \to 0$ boundary.

### Mechanism

In each SQUAREM iteration for VBEM:

1. **Extrapolation** pushes ghost-component $\alpha_i$ negative (line 933–934 of `em_solver.cpp`)
2. **Clamping** to `EM_LOG_EPSILON = 1e-300` (line 936)
3. **Stabilization E-step** computes `digamma(1e-300) ≈ -1e300`, which yields `exp(-1e300) = 0` in the posterior, so `em_totals[i] = 0`
4. **M-step** produces `alpha_new = 0 + 0 + prior ≈ 3e-6`
5. **Next iteration**: `state0 = 3e-6`, SQUAREM extrapolates, overshoots, clamped to `1e-300` again

Each micro-oscillation contributes $|3 \times 10^{-6} / \Sigma - 10^{-300} / \Sigma| \approx 3 \times 10^{-12}$ to the L1 norm. Summed over millions of ghost components, this creates a noise floor of $\sim 10^{-6}$ — exactly at the convergence threshold.

MAP-EM doesn't have this problem: `log(theta)` has gentle curvature near zero, and clamping to `0.0` with renormalization is a clean projection that stays at exactly zero.

### Evidence

| BAM | Components | VBEM iters | MAP iters | VBEM converged? |
|-----|-----------|------------|-----------|-----------------|
| Small (MEG-01 sub) | 179,704 | 74 | 121 | Yes |
| Large (MEG-01 full) | 5,777,023 | 333 (MAX) | 123 | **No** |
| Benchmark (dirty sim) | 216,179 | ~similar | ~similar | Marginal |

The problem scales with dimensionality: more ghost components → higher L1 noise floor.

## Solution: Variational Zero-Forcing

Two changes to `run_squarem()`, VBEM branch only.

### Change 1: Clamp extrapolated alpha to its prior (not to machine epsilon)

**Current** (line 935–936):
```cpp
if (state_extrap[i] < EM_LOG_EPSILON)
    state_extrap[i] = EM_LOG_EPSILON;
```

**Proposed**:
```cpp
if (state_extrap[i] < prior[i])
    state_extrap[i] = prior[i];
```

**Rationale**: The prior $\alpha_0$ is the variational analogue of "no data." When SQUAREM pushes a component below its prior, the data is actively saying this component should not exist. Clamping at the prior rather than at machine epsilon means:

- `digamma(prior)` is well-behaved (prior ≈ $10^{-10}$ → digamma ≈ $-10^{10}$, not $-10^{300}$)
- SQUAREM's quadratic model is accurate near this boundary
- The component sits at its prior — exactly what a Bayesian with no data should do

**Edge cases**:

- `prior[i] = 0` (ineligible components): `max(alpha, 0)` clamps to 0. The log-weight computation already guards with `digamma(max(alpha[i], EM_LOG_EPSILON))` on line 737. No issue.
- `prior[i] = EM_LOG_EPSILON` (eligible but zero coverage): These components have no EC entries, so `em_totals[i] = 0` always. They sit stably at `EM_LOG_EPSILON` regardless. Change 1 is neutral; Change 2 catches them.

### Change 2: Zero-force components with no net evidence

After the stabilization step, check if each component's posterior alpha exceeds its prior. If not, the data provides zero or negative net evidence — collapse it to the prior.

**Rationale**: The VBEM M-step (line 756) produces:

$$\alpha_{\text{new},i} = \text{unambig}_i + \text{em\_totals}_i + \text{prior}_i$$

So `alpha_new[i] <= prior[i] * (1 + tol)` is equivalent to:

$$\text{unambig}_i + \text{em\_totals}_i \leq \text{prior}_i \cdot \text{tol}$$

The data has contributed negligible evidence for this component. In variational inference, the KL divergence between the posterior and prior is minimized when they are equal — collapsing to the prior is the principled action.

**Irreversibility is a feature**: Once collapsed, `digamma(prior)` yields an extremely negative log-weight (e.g., $-3.3 \times 10^5$ for prior $= 3 \times 10^{-6}$), `fast_exp` returns 0, `em_totals` stays 0, alpha returns to prior. Self-reinforcing fixed point. Components cannot "wake up," but they shouldn't — any component with genuine support has accumulated evidence far exceeding $\text{prior} \times \text{tol}$ within the first few iterations.

**The tolerance uses the prior itself as the scale**: No arbitrary absolute threshold like $10^{-8}$. The relative tolerance `VBEM_ZERO_FORCE_REL_TOL` controls how much evidence (relative to the prior) counts as "no evidence."

### Implementation

Both changes fold into the existing convergence check loop (no extra pass):

```cpp
// --- New constant (with existing constants, line ~43) ---
static constexpr double VBEM_ZERO_FORCE_REL_TOL = 1e-6;
```

```cpp
// --- run_squarem(), VBEM branch ---

// SQUAREM extrapolation (lines 928-937):
// ... existing r_vec, v_vec, step computation ...
for (size_t i = 0; i < nc; ++i) {
    state_extrap[i] = state0[i] + 2.0 * step * r_vec[i]
                    + step * step * v_vec[i];
    // Change 1: clamp to prior, not machine epsilon
    if (state_extrap[i] < prior[i])
        state_extrap[i] = prior[i];
}

// Stabilisation step (unchanged)
vbem_step(state_extrap.data(), ...);

// Change 2: zero-forcing + convergence check (replaces lines 945-958)
double sum_old = 0.0, sum_new = 0.0;
for (size_t i = 0; i < nc; ++i) {
    // Zero-force: if posterior <= prior, data provides no evidence
    if (state_new[i] <= prior[i] * (1.0 + VBEM_ZERO_FORCE_REL_TOL)) {
        state_new[i] = prior[i];
    }
    sum_old += state0[i];
    sum_new += state_new[i];
}
double delta = 0.0;
if (sum_old > 0.0 && sum_new > 0.0) {
    double inv_old = 1.0 / sum_old;
    double inv_new = 1.0 / sum_new;
    for (size_t i = 0; i < nc; ++i) {
        delta += std::abs(state_new[i] * inv_new
                        - state0[i] * inv_old);
    }
}
```

### Constants

All constants in `em_solver.cpp` after this change:

| Constant | Value | Purpose |
|----------|-------|---------|
| `EM_LOG_EPSILON` | `1e-300` | Floor for `digamma()` input guard and ineligible component alpha |
| `SQUAREM_BUDGET_DIVISOR` | `3` | Each SQUAREM iter costs 3 EM steps |
| `ESTEP_TASK_WORK_TARGET` | `4096` | Work units per parallel E-step task |
| **`VBEM_ZERO_FORCE_REL_TOL`** | **`1e-6`** | **Relative tolerance for zero-forcing: collapse to prior if `alpha <= prior * (1 + tol)`** |

`VBEM_ZERO_FORCE_REL_TOL = 1e-6` means: for a typical OVR prior component with $\alpha_0 \approx 3 \times 10^{-6}$, the absolute evidence threshold is $3 \times 10^{-12}$ fragments — far below any biologically meaningful signal. For the gDNA component with $\alpha_0 \approx 0.5$, the threshold is $5 \times 10^{-7}$ fragments. Both are safely negligible.

## Files Modified

Only **one file**: `src/rigel/native/em_solver.cpp`

Changes:
1. Add `VBEM_ZERO_FORCE_REL_TOL` constant (line ~46)
2. Change alpha clamping floor in SQUAREM extrapolation (line ~935)
3. Add zero-forcing to convergence check loop (lines ~945-958)

## Why This Works

### Oscillation eliminated

With `state_extrap[i] = max(state_extrap[i], prior[i])`:
- Stabilization step sees `alpha = prior`, computes `digamma(prior) ≈ -3.3 \times 10^5`, E-step returns 0 for that component
- M-step: `alpha_new = 0 + 0 + prior = prior`
- Zero-forcing: `prior <= prior * (1 + 1e-6)` → collapsed to `prior`
- Next iteration: `state0[i] = prior`, `state1[i] = prior`, `r[i] = 0`, `v[i] = 0`
- Component drops out of SQUAREM acceleration entirely. L1 contribution = 0.

### Convergence restores

With millions of ghost components contributing exactly 0 to the L1 norm, only the $O(10^3)$ components with genuine data support determine convergence. The L1 norm over these active components behaves identically to MAP-EM (smooth, monotonically decreasing).

### No accuracy loss

Components forced to their prior had posterior ≤ prior — the data says they shouldn't exist. Setting them to exactly their prior is the Bayesian-optimal action. The only change from current behavior: ghost components that previously oscillated at $\sim 10^{-300}$–$10^{-6}$ now sit stably at $\sim 10^{-10}$. Both are indistinguishable from zero at any biologically meaningful scale.

## Expected Performance Impact

### Convergence improvement (primary goal)

The mega-locus should converge in ~100–150 SQUAREM iterations (comparable to MAP-EM) instead of hitting the 333-iteration cap.

### Per-iteration cost (unchanged)

The zero-forcing adds a single comparison per component to the convergence check loop — negligible. The E-step cost is unchanged: ghost components still appear in the EC data structure but their log-weights are so negative that `fast_exp` returns 0. The SIMD early-zero skip catches them per-lane but not necessarily per-vector (ghost components are not contiguous in the component array).

### Wall-clock improvement

On the MEG-01 mega-locus (previous profiling):
- Current: 333 iters × 1.91s/iter = 636s
- Expected: ~120 iters × 1.91s/iter = ~229s
- Savings: ~407s (64% reduction in mega-locus EM time)

## Alternatives Rejected

| Alternative | Problem |
|-------------|---------|
| Scale L1 by $\sqrt{N}$ | No principled scaling. $\sqrt{5.8M} = 2408$; $\delta$ becomes $2.4 \times 10^{-3}$ — may trigger premature termination |
| Scale L1 by $N$ | $\delta$ becomes 5.8 — converges on iteration 1 |
| L∞ norm | Masks symptom (summation noise) but doesn't fix cause (boundary oscillation). Still wastes iterations bouncing |
| Relative L∞ | Requires choosing $\epsilon$ floor (new parameter). Doesn't fix boundary oscillation |
| ELBO convergence | `sum(lgamma(alpha))` over 5.8M components per iteration — prohibitively expensive |
| MAP warm-start → VBEM polish | Doubles code complexity. Awkward transition: MAP zeros become VBEM positive values. New parameter (how many polish iterations?) |
| SQUAREM damping | Scalar step-length can't be both aggressive (real components) and gentle (boundary). Reduces acceleration when it's needed most |

## Validation Plan

### Step 1: Compile and test
```bash
pip install --no-build-isolation -e .
pytest tests/ -v --ignore=tests/test_calibration.py
```

If golden outputs drift (expected: small changes in low-abundance VBEM components):
```bash
pytest tests/ --update-golden
```

### Step 2: Benchmark
Re-run MAP vs VBEM benchmark on both conditions:
```bash
python -m scripts.benchmarking run -c scripts/benchmarking/configs/map_vs_vbem.yaml --force
python -m scripts.benchmarking analyze -c scripts/benchmarking/configs/map_vs_vbem.yaml \
    -o /scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_benchmarks/results/map_vs_vbem_v2
```

### Step 3: Verify

| Metric | Criterion |
|--------|-----------|
| VBEM Pearson R (clean) | ≥ 0.9928 (no regression from baseline) |
| VBEM Pearson R (dirty) | ≥ 0.9929 (no regression from baseline) |
| VBEM WARE (clean) | ≤ 0.156 (no regression) |
| VBEM WARE (dirty) | ≤ 0.158 (no regression) |
| VBEM wall time (dirty) | Reduced vs baseline 282s |
| Mega-locus convergence | VBEM converges (does not hit iteration cap) |

### Step 4: Real BAM validation (if available)

Run on MEG-01 BAM to verify mega-locus (159K transcripts, 5.8M components) converges and wall time improves.

## Future Work (Not in This Change)

### Active set elimination (Phase 3)

Once zero-forcing establishes which components are genuine, physically remove collapsed components from the EC data structure mid-EM. This would reduce E-step cost proportionally to the active fraction (likely <5% of total components). More invasive — requires rebuilding or masking the CSR structure. Deferred.

### Per-locus convergence diagnostics

Add iteration counts and convergence status to `loci.feather` output (requires changes to `estimator.py` and the C++ profiling structs). Useful for monitoring but not required for this fix.
