# Tripartite Prior Implementation Plan

## Summary

Decouple prior policy (Python) from EM engine (C++). Three component-type-specific
priors replace the uniform OVR treatment:

| Component | Prior α | Source |
|-----------|---------|--------|
| **mRNA** | `alpha_flat + OVR_share` | Unchanged — data-driven OVR |
| **nRNA** | `nrna_sparsity_alpha` (default 0.9, < 1) | Dirichlet sparsity prior |
| **gDNA** | `1.0 + gdna_prior_scale × gdna_init` | Spatial EB fragment count |

Default EM mode switches from MAP-EM to VBEM for SQUAREM stability with α < 1.

---

## Phase 1: Configuration Updates (Python — `config.py`)

1. Add `nrna_sparsity_alpha: float = 0.9` — must be < 1.0 for sparsity
2. Add `gdna_prior_scale: float = 1.0` — tunes EB prior strength
3. Change default: `mode: str = "vbem"` — enables digamma-based sparsity

## Phase 2: Prior Array Construction (Python — `locus.py`)

Build the complete `prior[T+N+1]` array explicitly before passing to C++:

- **mRNA [0, T)**: Set to `alpha_flat` (C++ OVR adds coverage on top)
- **nRNA [T, T+N)**: If eligible → `nrna_sparsity_alpha`; if dead → `0.0`
- **gDNA [T+N]**: `1.0 + gdna_prior_scale × gdna_init` (flat base + EB anchoring)

The Python prior now carries actual α values, not just a binary mask.

## Phase 3: Engine Modification (C++ — `em_solver.cpp`)

Modify `compute_ovr_prior_and_warm_start()`:

1. Accept `n_t` parameter to know how many mRNA components there are
2. **OVR restricted to mRNA only**: coverage sharing loop only adds to indices `[0, n_t)`
3. **nRNA + gDNA priors passed through**: the pre-computed α values from Python are
   preserved as-is in `prior_out`
4. Warm-start initialization unchanged (intronic evidence still seeds θ₀)

## Phase 4: Post-EM Cleanup (C++ — `em_solver.cpp`)

After VBEM convergence:
- Truncate: `if (theta[i] < 1e-5) theta[i] = 0.0`
- Re-normalize remaining θ to sum to 1.0

This snaps VBEM's asymptotically-small components to exact zero.

---

## Key Design Decisions

- **VBEM over MAP-EM**: digamma(α<1) creates smooth sparsity penalty;
  MAP-EM's max(0, N+α−1) creates discontinuity that destabilizes SQUAREM
- **OVR kept for mRNA only**: multi-mapper resolution still needs data-driven prior
- **gDNA prior = 1.0 + scale × gdna_init**: the +1.0 ensures a flat (non-sparsifying)
  base even when EB estimate is zero; scale parameter allows tuning
- **No changes to EM math**: VBEM M-step and E-step are unchanged; only prior
  construction and post-convergence cleanup are modified

## Validation Plan

Re-run `synthetic_sim_sweep.py` with `nrna_sparsity_alpha ∈ {0.5, 0.7, 0.9, 0.95}`:
- Target nRNA FP rate < 10% (currently 83%)
- mRNA accuracy unaffected
- gDNA accuracy improved at low SS
