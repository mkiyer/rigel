# Calibration solver — performance baseline and optimization targets

**Measured 2026-07-26 at HEAD `922abfcb`.** This is the starting point for the optimization phase. Re-measure
before and after every change; do not trust these numbers once the code moves.

---

## 1. The baseline

One condition of the standard suite (`gdna300_ss_0.50_nrna_present_capture_off`), **3,397 chain nodes**,
single-threaded (`OMP_NUM_THREADS=1`):

| | wall clock |
|---|---|
| `calib_refit_iters=0` (pass-0 only) | **0.50 s** |
| `calib_refit_iters=1` (production: pass-0 + fit + re-solve) | **0.90 s** |

**⚠ Scale is the whole problem.** 3,397 nodes is a 10 Mb synthetic. A human index is ~300× that in sequence,
so a linear extrapolation puts a production calibration at **~4–5 minutes** — and the NPMLE fit may scale
worse than linearly in cells. **Nobody has profiled this at genome scale; do that first (§4).**

## 2. Where the time actually goes — SELF time, `refit=1`

```
   ncalls  tottime            function
      902    0.206  (23%)     npmle._lse                          ─┐
        2    0.048  ( 5%)     npmle._em_weights                    ├── CLUSTER B: the hyperprior EM
                                                                  ─┘
     3640    0.124  (14%)     enrichment_frame.residual_level     ─┐
     3640    0.110  (12%)     bp_solver._peel_share                │
        4    0.059  ( 7%)     bp_solver._relay                     ├── CLUSTER A: the relay's
    10920    0.046  ( 5%)     scipy polygamma  (+0.013 zeta)       │   per-node SCALAR path
     7280    0.015  ( 2%)     enrichment_frame.peel_continue_share │
    95088    0.009  ( 1%)     bp_solver._damp                      │
    99219    0.014  ( 2%)     numpy.asarray                       ─┘
        6    0.082  ( 9%)     simplex_logodds._solve_ambig_logodds     (the 2-D (λ,τ) cube)
       24    0.045  ( 5%)     simplex_logodds._lse
```

**Two clusters of roughly equal size dominate: A ≈ 34 %, B ≈ 28 %.**

## 3. The targets, in priority order

### ⭐ A. The relay's per-node scalar path — ~34 %, and the shape of the fix is known

`_relay` is a sequential Gauss-Seidel scan (verified: every one of `rg pg pp pn mg mp mn tau` is read at the
source index and written at the destination index, so a later node reads what an earlier one wrote). It
therefore cannot be vectorised, and its arithmetic is plain Python floats — **which is fast**
(0.083 µs/op measured).

**The problem is that it calls numpy anyway.** `_peel_share(i, …)` is invoked **once per node per direction**
(3,640 calls) with a *scalar* index, and inside it runs `residual_level` (which calls
`scipy.special.polygamma`, 10,920 calls), `peel_continue_share`, and a stack of `np.where`/`np.asarray` on
0-d arrays. A numpy op on a scalar measured **15.7× slower** than the equivalent Python float expression
(1.307 µs vs 0.083 µs).

So the win is **not** to vectorise the relay (impossible) and **not** to merge it with `_transport`
(measured, rejected — see §5). It is to give the relay's peel a **scalar-native path**: plain floats,
`math.*` instead of `np.*`, and a scalar trigamma instead of `scipy.special.polygamma` on a 0-d array.
`residual_level` and `peel_continue_share` are small closed forms; a `_scalar` twin of each, exercised by the
relay and asserted equal to the array version in a unit test, is the obvious shape.

⚠ `polygamma(1, k)` is `zeta(2, k)`; scipy's generic `polygamma` wrapper is doing dispatch work 10,920 times
for a single float. A direct scalar trigamma (or a small cache, since `k ≥ 1` and the values repeat) removes
most of that 0.059 s.

### ⭐ B. The NPMLE hyperprior fit — ~28 % of `refit=1`, and it is one function

`npmle._lse` is **902 calls / 0.206 s self** — the single largest self-time entry in the whole run. It is
called from `_em_weights`, i.e. the EM loop of `DensityNPMLE.fit`. Look for: recomputation that could hoist
out of the EM iteration, a `logsumexp` that could use the array-wide `scipy.special.logsumexp` (or an
in-place stable form) instead of a per-row Python-level reduction, and whether the EM's convergence
tolerance/iteration cap is doing more passes than the result needs. **Check the iteration count first** —
0.206 s over 902 calls is 0.23 ms per call, so the cost may be *count*, not per-call work.

### C. `_solve_ambig_logodds` — 9 %, 6 calls

The 2-D (λ,τ) cube. Already batched (`_AMBIG_BATCH`) to bound peak memory. Check whether the batch size is
right for cache, and whether the cube is being materialised where a marginal would do.

### D. Allocation churn

**99,219 `numpy.asarray` calls** and 95,088 `_damp` calls for 3,397 nodes. Most `asarray` calls are on
inputs that are already float64 arrays — the codebase defensively re-wraps. Cheap to remove where the type is
guaranteed by construction, but measure: it is only ~2 % here, and it may matter more at scale.

## 4. Do this FIRST — profile at genome scale

Every number above is from a 3,397-node chain. Before optimizing, establish how each cluster **scales**:
build or find a substantially larger index and re-profile. A hotspot that is 23 % at 3.4 k nodes and 60 % at
1 M nodes is a different problem from one that stays flat. In particular the NPMLE's cost is in *cells*, not
nodes, and may not track the chain length at all.

## 5. ⛔ Already measured and rejected — do not re-do

| idea | verdict |
|---|---|
| merge `_relay` and `_transport` into one polymorphic routine | **15.7× slower** on the relay's per-op arithmetic (0.083 → 1.307 µs) = **+5.5 s per calibration** at 74,494 nodes, to save ~60 lines. The twin relationship and its three deliberate edge-case differences are documented at both definitions |
| vectorise `_relay` | impossible — it is a sequential Gauss-Seidel scan (mechanically verified) |

## 6. The rules for this phase

* **Every change must be verified BIT-IDENTICAL** on the 32-condition benchmark at **both** refits:
  `bash scratchpad/verify_clean.sh ARMNAME` (runs ruff + pytest + the A/B and prints `32/32 IDENTICAL`).
  A performance change that moves a number is a **bug**, not a trade — unless it is a deliberate,
  owner-approved algorithmic change, in which case it needs a full A/B on both axes (see
  `SESSION_2026_07_26_HANDOFF_15.md` §4).
* **Measure, change, re-measure.** Never land an optimization on reasoning alone — this codebase has already
  produced two "obvious" wins that measured backwards.
* Keep the code **clean, concise, readable**. A 5 % win that costs legibility is not a win; the owner has
  been explicit that the code must reach Phase 2 *"efficient, clear, concise, readable, maintainable, not
  over-engineered, as simple as possible."*
* Floating-point **association order** changes results. Reordering a sum, switching to a fused reduction, or
  swapping `a*b/c` for `a*(b/c)` can break bit-identity. If it does, that is a signal to stop and think, not
  to relax the gate.
