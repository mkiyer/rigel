# Threads — the port's step (iv), the design (2026-09-17; put to the owner BEFORE building)

The frame: `ISSUES: performance-memory-bounded-solve` ③ (iv). The rulings this plan stands on: `DESIGN.md`
§6b.15.1–§6b.15.2 (the locus block is the unit of the solve, and the only information that crosses a block
boundary is the policy's library), §6b.15.3 (ψ's read-out is chunk-exact: every slot is solved on its own cube),
§6b.15.9 (threads were the wrong tool for the PYTHON sweep and the executor waited for the port — this plan is
what the port has made possible), and the owner's order of 2026-09-17: finish the block in native before threads,
and design threads on paper first. Nothing here is built.

## Where the time is (this tree, after the tables' allocations)

MO_3021 replayed (`sweeps_MO_3021_step11`), blocks of 5,000 slots — 426 blocks over 2,087,476 slots, 357,738 of
them solved per ψ call-set and 51,734 of those AMBIG. The four native calls are timed by wrappers (a nanobind call
is invisible to cProfile), the rest by cProfile; the machine quiet (`s10/kernel_times_quiet.log`,
`s10/psi_sweep_compare.log`, `s10/grid_per_sweep.log` in the synced scratchpad).

The first sweep and the refit sweeps are DIFFERENT PROBLEMS: the first sweep solves on K = 101 cells at the
window L = 10; the three refit sweeps on the landscape's own bracket, K = 202 at L = 20.09 (`calibrate` derives
the bracket from the fitted prior's support, `landscape.required_logodds_window`). Every per-cell cost — ψ, the
factory rows, the gDNA arm, the cache key's bytes — is therefore 2× on a refit sweep for the same slots (measured:
ψ's per-call time ratio 2.02 at the median over the 852 matched calls, the slot sets identical). In a production
run sweep 1 also MISSES the message cache and runs the whole layer at K = 202; sweeps 2–3 are served.

| stage, per sweep | first sweep (K = 101) | refit sweep (K = 202, cache-served) | scales with |
|---|---|---|---|
| ψ native, the self-solve + the final solve (`psi_solve`, 852 calls) | 1.57 s | 3.19 s | K × (n_tilt + 2) per AMBIG slot, K per single-strand slot |
| the builders native (`transfer_prepare`) | 1.61 s | served | K, the nodes |
| the pass native (`transfer_pass`, both directions) | 0.99 s | served | the hops, K |
| the policy's solve native (`transfer_solve`) | 0.12 s | served | |
| the message cache's key (`MessageCache.key`: blake2b over every array of the block context — the factory rows are n × K × 8 bytes of it, 1.7–3.4 GB a sweep) | — | 2.88 s | K |
| the gDNA arm (`landscape.logprior`: `np.interp` over an (m, K) matrix per block) | — | 1.51 s | K |
| the factory rows (`calibrate.__getitem__` → `density_lambda_factor` → `_log_negbinom`; `gammaln` is 85 % of it) | 0.46 s | 0.90 s | K, the introns |
| the cache's dense rows (`MessageCache.message`) | — | 0.13 s | |
| the glue: `solve_chain`'s loop, `_psi`'s row add (`factory + lam_rows`, an (n, K) add per block, 0.28 s), `build_region_init`'s strand evidence and factor precision, `block_slice`, the checks, the write-back | ~1.9 s | ~1.3 s | the slots |
| **the sweep** | **6.7 s** (native 64 %) | **10.0 s** (native 32 %) | |

On the deep library (VCaP, 18.6 M fragments, 8 threads, the solve commit's report `perf/port_solve_2026-09-17/`)
the four sweeps are 135 s of 264 s: ψ 53 s (the self-solves 20, the final solves 33), the pass 37 s, the builders
18 s, the solve 2 s, the rest the Python above. The allocation step moved these by about a second
(`perf/block_alloc_2026-09-17/`).

## What is parallel, and what it costs

### (iv-a) ψ over slots, inside `psi_solve` — first, and alone

Every slot is solved on its own cube with nothing shared but the read-only `Grid` and the delivered `Rows`
(`psi_kernel.cpp`: `solve_slot` reads `SlotInputs` and writes the slot's four outputs; no reduction crosses
slots). So a partition of the dispatcher's slot list over threads, each thread with its own `Scratch`, is
BIT-IDENTICAL to the serial loop for ANY thread count — the same arithmetic per slot in the same order within the
slot. That is the gate: the replay of the step11 capture at 1, 2 and 8 threads (BIT-IDENTICAL on all four
sweeps), the suite, the three identity references; then two interleaved timing pairs on VCaP.

* The partition is by COST, not by count: an AMBIG slot is n_tilt + 2 = 26 columns against a single-strand
  slot's one, and the AMBIG slots are 14 % of the solved slots but ~82 % of the cells. The EM's E-step already
  partitions its tasks by cumulative cost (`em_solver.cpp::parallel_estep`); the same greedy cut over the slot
  list, cost 26 or 1 per slot, balances the threads.
* The pool is the EM's `EStepThreadPool` (`thread_pool.h`: persistent workers, `run_parallel(fn(tid))`, the
  caller as worker 0). ONE pool for the module's lifetime — a static, built at the requested size and rebuilt
  when the size changes — never one per call: 852 ψ calls a sweep × 4 sweeps × 7 spawns at 20–50 µs is a second
  a run for nothing.
* The GIL: the EM releases it around its solve (`nb::gil_scoped_release`); ψ does the same. Nothing Python
  runs concurrently today, so this is hygiene until (iv-b).
* The thread budget: calibration has no thread knob. The EM's is `EmConfig.n_threads` (0 → every core;
  `estimator.py` passes it to the native call). PROPOSED: `CalibrationConfig.n_threads` with the same
  semantics, fed by the pipeline from the one `--threads` budget the scan and the EM already share. It is a
  resource budget, not a tunable of the answer (the answer is bit-identical at every value), so it is not a
  magic number — but it is a new config field, which is the owner's call.
* Expected: ψ 1.57 → ~0.25 s on the first sweep and 3.19 → ~0.45 s on a refit sweep at 8 threads (the machine
  has 12 performance cores of 16, an M3 Max; the AMBIG chunks bound the balance; the per-slot working set, ≤ 170 KB at
  K = 202, sits in cache, so the kernel is compute-bound on `exp`). On VCaP: ψ 53 s → ~8 s of 264 s, the run
  ~0.83×. Memory: a `Scratch` per thread, ~50 KB.

### (iv-b) the blocks in parallel — after the Python per block has shrunk

The blocks are independent (§6b.15.1–§6b.15.2), so any order and any interleaving give the same answer; each
block's write-back is a disjoint slice of `out`; `AssertionCounts.absorb` is a sum (order-free); the cache's
`put` needs a lock; the diagnostic capture (instruments only) stays serial. Two forms:

* PYTHON THREADS over blocks (`concurrent.futures.ThreadPoolExecutor`), the four kernels releasing the GIL.
  Bounded by the serial Python (Amdahl): on the first sweep 2.4 s of 6.7 s is Python, so the sweep is at best
  2.4 + 4.3/T = 2.9 s at T = 8 (2.3×); on a refit sweep 6.8 s of 10.0 s is Python, so at best 7.2 s (1.4×) —
  and once (iv-a) has taken ψ to ~0.4 s, a refit sweep is ~7 s of Python that threads over blocks cannot
  touch. So (iv-b) is worth building only after the Python per block has been made small.
* ONE NATIVE CALL PER BLOCK — `solve_block` in C++: the self-solve ψ, the builders, the two passes, the solve,
  the final ψ, the write-back — and a C++ pool over blocks. The Python per block is then the slicing and the
  loop (~1 ms). This is the block-in-native end-state the owner named; it subsumes the items below and is the
  larger step.

The Python per block that must shrink first — each EXACT (bit-identical, or inside the derived budget where
stated), each its own commit, in the order of its size:

1. **The cache key on the factory's INPUTS, not its rows** (2.88 → ~0.05 s a refit sweep). The rows are a
   pure function of the background's parameters, the per-slot count and eff (already per slot on the factory)
   and the grid; digesting those keeps the key content-keyed — every input the layer reads is digested, and a
   changed background, count, eff or grid still misses — at 1/K of the bytes. The context field stays the rows
   (the layer reads them); only what the KEY hashes for that field changes. Owner: this changes what
   "content-keyed" means for one field.
2. **The gDNA arm inside ψ** (1.51 s a refit sweep, and 4 MB a block of traffic). The kernel already takes a
   per-slot prior row; give it the landscape's curve (`log_rho`, `logP`) and the per-slot (mass, eff) instead
   of the (m, K) matrix and let it interpolate at its own K cells. `np.interp` and a C++ linear interpolation
   agree to the ulp only when written identically; hold it to the replay's tolerance budget.
3. **`_psi`'s row add** (0.28 s): the kernel takes the factory row and the delivered row apart (it sums per cell
   already); no (n, K) add in Python.
4. **The factory rows** (0.46–0.90 s): lgamma is 85 % of `_log_negbinom`, so a serial port saves ~0.1 s a sweep
   and moves numbers (libm against cephes) — NOT worth its own port and protocol; inside one native call per
   block, threaded over rows, it costs ~0.1 s at 8 threads. Part of the block port, not a step.
5. **`build_region_init`'s glue** (`strand_evidence`, `density_factor_precision`, the where's; ~0.3 s): part of
   the block port.

### What is not proposed

* Threads inside a pass: a hop reads the row the previous hop wrote, in chain order — sequential by
  construction. The builders are per node but 1.6 s a first sweep only; parallel over blocks or not at all.
* OpenMP: the tree links none; the EM's pool is the idiom.
* The underflow cut in ψ — skipping the cells whose log-weight lies below the truncation constant under the
  cube's maximum, a plausible 2× on ψ at a tolerance — not before (iv-a) is measured; a separate A/B if at all.
* An arena for the per-block tables: priced at ≤ 0.3 s of a run's first sweep (the first-touch faults the
  allocation step left; the refit sweeps allocate none of these tables under the cache) — not built.

## The order proposed

1. (iv-a) ψ over slots — one commit; the replay at 1/2/8 threads BIT-IDENTICAL, the suite, the references;
   two interleaved pairs on VCaP at `--threads 8`.
2. The cache key on the factory's inputs (exact; the largest Python item of a refit sweep).
3. The gDNA arm inside ψ (tolerance-gated); `_psi`'s row add with it.
4. The block in one native call with a C++ pool over blocks — or Python threads over blocks with GIL-releasing
   kernels — judged on the Python floor measured after 2–3.

## Open for the owner

1. The thread budget's home: a `CalibrationConfig.n_threads` mirroring `EmConfig.n_threads`, fed by the
   pipeline's `--threads`?
2. The cache key digesting the factory's inputs instead of its rows.
3. After (iv-a): the block in one native call (the end-state, the larger step), or Python threads over blocks
   (sooner, bounded by the Python)?
