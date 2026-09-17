# ψ native — the port's step (iii), the design (2026-09-17)

The frame: `ISSUES: performance-memory-bounded-solve` ③ (iii); the rulings ψ keeps are `DESIGN.md` §6b.15.3
(chunk-exact read-out), §6b.15.6 (one solver, float64), §6b.15.11 (the θ nodes follow the peak), §6b.15.13 (the
tilt atom) and §6b.15.14 (the ½-quantile read on λ); the derivations `EQUATIONS.md` §9e–§9f. The owner's
one-path ruling (2026-09-17) shapes the port: the C++ is the only implementation, the Python it replaces is
deleted in the same commit, and a small floating-point tolerance between the two is accepted for the speed.

## Where the time is

MO_3021's first sweep replayed (cProfile, after step (ii)): 12.4 s, of which ψ is 4.95 s — the cube's assembly
2.7 s (`_psi` 1.3, the strand term 0.67, the θ window 0.38, the cube rows 0.31), the read-out 1.8 s (`exp` and
the sums 0.7, the log-sum-exp 0.54, the quantile 0.35), the dispatcher 0.34. On VCaP ψ is 46 + 57 s of the
190 s sweep. The cost is numpy's passes over `(rows, K, K_t + 2)` temporaries — some fifteen of them per tile —
not the arithmetic: one AMBIG slot is 101 × 26 = 2,626 cells, one exp each; a single-strand slot 101.

## The shape

ONE kernel, per slot, in one pass over the slot's own cube (21 KB, in cache), for every slot the dispatcher
selects; nothing is tiled and nothing is shared between slots, so the read-out is chunk-exact by construction
and the kernel is trivially parallel over slots (step (iv)'s ψ half needs no block threading).

```
psi_solve(slots, u_pos, u_neg, allow_pos, allow_neg, fg_ref, fpos_ref, fneg_ref, kappa, od_g, od_r, lam,
          gdna_logprior | None, lam_logprior | None, cube rows (packed, see below), n_tilt,
          out: f_g, f_pos, f_neg, var_gdna)
  for each slot:
    var        = the frozen strand variance at the reference composition          (transfer_rows.h)
    columns    = 1 (single strand: τ = ±1)  |  n_tilt nodes across the θ window + the two atoms
    ψ[k, t]    = strand(u_pos, n, p(f_g, τ)) + ½ log f_g + ½ log(1 − f_g) + gdna prior + λ-factor row
                 + the delivered row read at (f_g, τ) + log weight;  the ruled-out atom = −inf
    post       = exp(ψ − max);  Z, the θ-marginal per λ, Σ post·f_pos, Σ post·f_neg     (one pass)
    f_g        = the continuous ½-quantile of the θ-marginal on λ's midpoint edges, through σ
    Var log f_g= E[log² f_g] − E[log f_g]²  over the θ-marginal, clipped at 0
    w_pos      = Σ post·f_pos / Σ post·(f_pos + f_neg)   (½ when there is no RNA mass)
    compose    = (1 − f_g)·w_pos on the admissible strands;   a slot with n = 0 reports zeros
```

The dispatcher stays in Python and keeps its name and signature (`_solve_regions_logodds_all`, the entry point
of `sweep._psi`, `region_init` and `region_geometry`): the reference defaults, the `signal` mask, the
`{slot: CubeRow}` delivery packed into parallel arrays (slot, profile per strand with a presence bit, `u`,
total, opportunity, ρ_ref), one native call, a `RegionDeconv` back. `_block_rows` stays (the landscape and the
capture efficiency tile on it). `CubeRow` stays as the record the message layer delivers; its `at` moves into
the kernel.

## One implementation, tested through its own entry points

The module `_psi_impl` (`src/rigel/native/psi_kernel.cpp`, on `transfer_rows.h`'s shared pieces — the
sigmoid, the strand variance and term, `interp`) binds four things, all the SAME code the solve runs:

| binding | for |
|---|---|
| `psi_solve` | production — the dispatcher's one call |
| `psi_cube` | the gates that read ψ itself: the cube `(m, K, C)`, the tilt grid, the two strand-fraction grids, the log-weights — `test_vertex_reference.py`'s reference-measure, θ-quadrature, tilt-atom and cube-row gates |
| `posterior_median` | the ½-quantile read-out on a given posterior — the read-out gates and the replay's budget self-test |
| `compose` | the admissibility gates |

Deleted from `simplex_logodds.py`: `_psi`, `_solve_logodds`, `_tilt_window`, `_posterior_median_fg`,
`_row_moment`, `_lse`, `_compose`, `_gdna_arm`, `_rna_arm`, `_mixture_strand_loglik`, `strand_row_logodds`
(unused since step (ii)'s cleanup), `_log_fg`, `_log1m_fg`, the two class masks, `CubeRow.at`. The docstrings'
derivations that are not yet in `EQUATIONS.md` (the quantile read on λ, the closure by parametrisation) move
there under the move rule; the C++ carries the formulas and cites the gates.

The gates that read the internals are rewritten to the bindings: `test_vertex_reference.py` (the cube through
`psi_cube`, the read-out through `posterior_median`, `compose`; its reference integrals keep their own
five-line numpy strand term — a test's oracle, not a production path), `test_sweep_replay_tolerance.py` (the
same), `test_strand_likelihood_reference.py` (the single-strand cube against the layer-4 two-component
reference plus the Jeffreys arms), `test_sweep.py` (`log_expit` for `_log_fg`), `test_transfer_rna_lanes.py`'s
bracket theorem (the cube through `psi_cube`). The derived-node-count gate passes `n_tilt = 60` explicitly
instead of monkeypatching a module global.

## The protocol, as for steps (i) and (ii)

1. A fresh capture of the committed cleanup tree first (`sweeps_MO_3021_step9`; step8 predates step (ii) and
   no longer replays bit-identical); step8 deleted for storage.
2. The kernel in-tree, built, NOT wired: a scratchpad harness wraps `_solve_regions_logodds_all` to run the
   Python and the native on every ψ call of sweep 0 — 852 dispatcher calls, both slot classes, the self-solve
   and the final solve — and compares the four outputs against the replay's derived budget, per slot.
3. Wire (the dispatcher's body), delete the Python, rewrite the gates, the suite, `sweep_replay.py replay
   --tolerance` on step9 for calls 0–3, the three references re-frozen with the reason, two interleaved
   timing pairs on VCaP at 8 threads against a worktree of the cleanup commit, the docs, one commit.

## Not in this step

Threads (step (iv): ψ over slots inside the kernel, the builders and the pass over blocks); the remaining
Python duplicates of the transfer (the per-hop pass kernel and the row constructors the unit gates recompute
with) — the convergence after this step, on the same pattern: bind the C++ pieces, rewrite the gates, delete.
