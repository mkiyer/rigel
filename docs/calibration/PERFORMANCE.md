# Calibration solver — performance, measured at GENOME SCALE

**Rewritten 2026-07-27.** The first version of this file was profiled on a 3,397-node synthetic and ranked
the targets from that. **Every one of those rankings was wrong**, and §1 is the reason. Re-measure before and
after every change; do not trust these numbers once the code moves.

---

## 1. ⭐ THE LESSON: the 10 Mb synthetic ranks the hotspots BACKWARDS

Profiling the same `refit=1` calibration at two scales — the 3,397-node synthetic, and a real cfRNA sample
against the human index v7 (752,654 regions ≈ 1.5 M chain nodes):

| cluster | 3.4 k nodes | 1.5 M nodes | |
|---|---|---|---|
| the relay's per-node scalar path (`_relay` cum) | 34 % | **81 %** | ranked 1st either way, but it is not "a third of the time", it is *the* calibration |
| the per-node ψ solves (`_solve_nodes_logodds_all` cum) | 9 % | **29 %** | the toy said "minor" |
| the NPMLE hyperprior EM (`npmle._lse`) | **28 %** | **0.7 %** | the toy's #1 hotspot is noise at scale |

The NPMLE inversion is the instructive one and it is structural, not incidental: **its cost is in collapsed
*cells*, not nodes**, and the cell count saturates — 902 `_lse` calls at 3.4 k nodes and 902 at 1.5 M. A
sub-agent produced a careful, correct, entirely wasted analysis of how to speed it up.

**So: profile on `cfrna:LBX0190` before ranking anything.** `scratchpad/perf_scale.py` runs both scales off
cached payloads (`~/Downloads/rigel_runs/cfrna/_calib_cache/*.pkl` — real samples, human index) and takes
~70 s per genome-scale run.

## 2. Where it stands

Genome scale, `refit=1`, `OMP_NUM_THREADS=1`, `cfrna:LBX0190`:

| | wall |
|---|---|
| 2026-07-26 (`c3c1b5ab`) | **266 s** |
| + scalar-native peel (`619465c8`) | 110 s |
| + cache-blocked ψ solves, scalar-native relay (`6d483004`) | 73 s |
| + the AMBIG log hoist | **67 s** |

**≈ 4×.** The 10 Mb synthetic moved 0.603 → 0.35 s (`refit=0`) and 1.128 → 0.62 s (`refit=1`).

Current profile (67 s; `tottime` self / `cumtime`):

```
  _solve_nodes_logodds_all  cum 49        ─ the per-node ψ solves          ~62 %
      _solve_ambig_logodds  self 12 / cum 25    the 2-D (λ,θ) cube
      _lse                  self 10             both solvers
      _mixture_strand_loglik self 5
      _local_loglik_logodds self 5 / cum 8      the 1-D single-strand ψ
      _regrid_global        self 2              (was 9.7 s before blocking)
  _relay                    self 11 / cum 18   ─ the sequential scan        ~22 %
  npmle (_lse + _em_weights) 1.4               ─                            ~2 %
```

## 3. What made it faster, and the two ideas behind it

**(a) A sequential scan should not call numpy.** `_relay` is Gauss-Seidel — it reads state an earlier
iteration wrote — so it cannot be vectorised, and it was calling array primitives once per node with a
*scalar*. A numpy op on a 0-d array costs 0.5–0.7 µs against ~0.02 µs for the float expression, and because
every `np.where` evaluates both arms, a node with no level still paid for four `log_ndtr` calls and a
`polygamma`. The fix is `*_scalar` twins of the shared primitives (`residual_level_scalar`,
`peel_continue_share_scalar`, `graft_frame_logvar_scalar`, `bp_solver._peel_share_scalar`) plus Python-list
operands (the `*_l` block in `node_sweep`). `_peel_share` cum went 181 s → 3 s.

**(b) The ψ solves were memory-bound, not FLOP-bound.** The 1-D path solved all 357,739 single-strand nodes
in one shot at K=256 — a **699 MB temporary** per intermediate, ~10 live, streamed from DRAM. Blocking to a
1 MB working set (`_SOLVE_BLOCK_BYTES`) runs the same arithmetic out of cache. Bit-identical because every
node solves independently and every reduction is *within a row*.

## 4. Bit-identity — how the gate is actually met

`bash scratchpad/verify_clean.sh ARM` (ruff + pytest + the 32-condition benchmark at both refits) is the
landing gate, but it takes ~22 min. The **inner loop** is `scratchpad/perf_identity.py`, which compares the
FULL solved per-node state (f_g / f_pos / f_neg, all three variances, the fused modes and precisions, and
every float array on `CalibrationResult`) over 8 conditions × both refits in ~8 s — strictly stronger than
the aggregate TSV, and 150× faster.

Two traps that cost real time here, both now guarded in `verify_clean.sh`:

* **an empty arm scored 32/32.** The comparison looped over the *new* arm's conditions, so a bench that died
  (OOM, under a concurrent genome-scale job) produced no rows and no differences. A gate that cannot fail is
  not a gate.
* **the `clean0` baseline was stale.** Unmodified HEAD no longer reproduced it — including the `mass` column,
  which comes from the accumulator payload and no solver change can reach; `_selfsolve_cache` had been
  regenerated underneath it. It reported 0/32 for a change that was in fact bit-identical. **Re-record the
  baseline from a `git stash` of HEAD, in the same session, at both refits, and if HEAD-vs-baseline is not
  32/32 then the baseline is what is broken.**

Three float facts this work established, each verified rather than assumed:

* `scipy.special.polygamma(1, k)` **is** `zeta(2, k)`, bit-for-bit (its `(−1)^(n+1)` and `Γ(n+1)` prefactors
  are exactly 1.0 at n=1) — and `polygamma` also evaluates `psi(k)` and discards it.
* `math.exp/log/sqrt` are bit-identical to their numpy counterparts here, but `np.maximum`/`np.minimum`
  **propagate nan** where a bare `if` swallows it, and `np.clip` **keeps x** inside the interval (preserving
  −0.0) rather than composing max-then-min. Hence `_fmax`/`_fmin`/`_fclip01` in `enrichment_frame`.
* numpy's **float32 `log` is monotone** on [0,1] — verified over all 1,065,353,217 values there — which is
  what licenses `max(log a, log b)` in place of `log(max(a, b))` on the AMBIG cube. Pinned by
  `test_float32_log_is_monotone_so_the_ambig_cube_may_hoist_it`.

## 5. ⛔ Measured and rejected — do not re-do

| idea | verdict |
|---|---|
| merge `_relay` and `_transport` | **15.7× slower** on the relay's per-op arithmetic. The twins' three deliberate edge-case differences are documented at both definitions; the relay is now scalar-native throughout, which makes a merge worse still |
| vectorise `_relay` | impossible — a sequential Gauss-Seidel scan (mechanically verified) |
| optimize the **NPMLE EM** | 0.7 % at genome scale (§1). A full analysis exists and found ~2–3 % *of that*. Not a target |
| in-place `psi` accumulation in `_local_loglik_logodds`; in-place `exp` buffer in `_lse` | **74.62 s vs 74.69 s — no effect.** Once blocking makes the arrays cache-sized the allocator reuses the same buffers, so removing allocations buys nothing. Reverted |
| block sizes below 256 KB | sharply worse (136 s at 64 KB) — per-call overhead takes over. The curve is flat 1–16 MB |

## 6. Where the remaining time is, if anyone wants it

The ψ solves are ~62 % and are now genuinely arithmetic-bound on the (m,K,K_t) cube; the honest next lever
there is **fewer cells**, not faster cells — `K = K_t = 60` on 12.6 % of nodes carries ~3× the cell count of
the entire 1-D path. That is a modelling question (grid resolution), not a cleanup one, so it is out of
scope for this phase.

`_relay` is ~22 % and is now ~11 s of straight-line Python over ~6 M edge-iterations. Nothing left short of
moving the scan to C++ — a real option, and a much larger change than this phase was scoped for.

## 7. ⏸ THIS PHASE IS PAUSED (owner, 2026-07-27)

Further optimization waits until the full pipeline is finished. What is left, and what it costs, is
§6. When it resumes, start by re-reading §1 and re-profiling at genome scale — these numbers will
have moved.

## 8. The rules for this phase

* **Every change must be verified BIT-IDENTICAL** at both refits (§4). A performance change that moves a
  number is a **bug**, not a trade.
* **Measure, change, re-measure — at genome scale.** Two ideas in this file measured flat or backwards, and
  the whole original target ranking was an artifact of the toy.
* Keep the code **clean, concise, readable**. A 5 % win that costs legibility is not a win.
* Floating-point **association order** changes results. Reordering a sum, fusing a reduction, or swapping
  `a*b/c` for `a*(b/c)` breaks bit-identity. If it does, that is a signal to stop and think, not to relax
  the gate.
