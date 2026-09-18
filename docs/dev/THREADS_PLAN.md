# Threads — the port's step (iv), the design (2026-09-17; the owner's answers 2026-09-18: `CalibrationConfig.n_threads` fed by `--threads` — YES; the cache key on the factory's inputs — YES; the block in ONE native call over Python threads — PREFERRED). (iv-a) LANDED 2026-09-18

The frame: `ISSUES: performance-memory-bounded-solve` ③ (iv). The rulings this plan stands on: `DESIGN.md`
§6b.15.1–§6b.15.2 (the locus block is the unit of the solve, and the only information that crosses a block
boundary is the policy's library), §6b.15.3 (ψ's read-out is chunk-exact: every slot is solved on its own cube),
§6b.15.9 (threads were the wrong tool for the PYTHON sweep and the executor waited for the port — this plan is
what the port has made possible), and the owner's order of 2026-09-17: finish the block in native before threads,
and design threads on paper first. Nothing here is built.

## Where the time is (this tree, on the deep library — the optimisation target from 2026-09-17)

VCaP (18.6 M fragments) at 8 threads, two back-to-back runs (`perf/vcap_baseline_2026-09-17/run{1,2}.json`): wall
263.6 and 265.5 s, the stages drifting 0.98–1.07 between the runs. Calibrate is 154 s, the sweep 134 s of it; 110 s are
outside calibration (the scan 35, the second pass 22, the fragment-length models 8.3 twice, quant 37, the index load
6.5). The sweep's stages over the four sweeps: ψ 21.6 (the self-solves) + 32.7 (the final solves), the pass 36.5, the
builders 17.2, the solve 2.0, the Python between the kernels ~24.

Per sweep (`sweeps_VCaP_step13` replayed; the four native calls timed by wrappers, the transfer kernels counted by the
census — a scratch copy with counters swapped in for the three calls, `s12/census_run.py`; the machine quiet;
`s12/kernel_times_vcap.log`, `s12/census_vcap.log` in the synced scratchpad). The chain is the annotation's: 2,087,476
slots in 426 blocks of 5,000, 357,738 solved per ψ call-set, 51,734 of them AMBIG. The first sweep and the refit sweeps are
DIFFERENT PROBLEMS: the first solves on K = 101 cells at L = 10; the three refit sweeps on the landscape's own bracket,
K = 233 here (`landscape.required_logodds_window`), so every per-cell cost is 2.3× there for the same slots. In
production sweep 1 also MISSES the message cache and runs the whole layer at K = 233; sweeps 2–3 are served.

| stage, per sweep | first sweep (K = 101) | refit sweep (K = 233, cache-served) | scales with |
|---|---|---|---|
| ψ native, the self-solve + the final solve (`psi_solve`, 852 calls) | 6.7 s | 15.5 s | K × (n_tilt + 2) per AMBIG slot, K per single-strand slot |
| the pass native (`transfer_pass`, both directions): 4.1 M hops, 2.1 M through a rule, 1.55 M compositions written, 1.39 M levels emitted; 1.43 M BLURS at a mean of 56 taps over K cells = 8.1 G multiply-adds | 10.4 s | served (~24 s at the miss) | the blurs × taps × K |
| the builders native (`transfer_prepare`): the RNA lanes 2.7, the gDNA lane 1.2, the claims 0.3, the splice faces 0.3, the terminus and alternative-splice rules 0.2 | 4.7 s | served (~11 s at the miss) | K, the nodes |
| the policy's solve native (`transfer_solve`) | 0.5 s | served | |
| the message cache's key (blake2b over the block context, the factory rows most of its bytes) | — | 3.2 s | K |
| the gDNA arm (`landscape.logprior`: `np.interp` over (m, K)) | — | 2.3 s | K |
| the factory rows (`_log_negbinom`; lgamma 85 %) | 0.5 s | 1.1 s | K, the introns |
| the glue | ~1.8 s | ~1.2 s | the slots |
| **the sweep** | **24.7 s** (native 91 %) | **23.3 s** (native 67 %) | |

The kernels are the sweep on the deep library: 91 % of a first sweep and, once ψ is threaded, the whole of what is
left there. MO_3021 (the sparse capture library the port was developed on) had a different shape — a first sweep 64 %
native, 0.26 M blurs, 0.08 M mapped hops — and is retired as a target.

## What is parallel, and what it costs

### (iv-a) ψ over slots, inside `psi_solve` — LANDED 2026-09-18 (267 → 214 s and 261 → 213 s at 8 threads (0.80 / 0.82), ψ 57.4 → 8.3 s, the sweep 135 → 88 s); the record is `DESIGN.md` §6b.15.5

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
* Expected: ψ 6.7 → ~1 s on the first sweep and 15.5 → ~2.3 s on a refit sweep at 8 threads (the machine has 12
  performance cores of 16, an M3 Max; the AMBIG chunks bound the balance; the per-slot working set, ≤ 200 KB at
  K = 233, sits in cache, so the kernel is compute-bound on `exp`). On VCaP: ψ 58 s → ~8 s of 264 s, the run
  ~0.81×. Memory: a `Scratch` per thread, ~50 KB.

### (iv-b) the blocks in parallel — after the Python per block has shrunk

The blocks are independent (§6b.15.1–§6b.15.2), so any order and any interleaving give the same answer; each
block's write-back is a disjoint slice of `out`; `AssertionCounts.absorb` is a sum (order-free); the cache's
`put` needs a lock; the diagnostic capture (instruments only) stays serial. Two forms:

* PYTHON THREADS over blocks (`concurrent.futures.ThreadPoolExecutor`), the four kernels releasing the GIL.
  Bounded by the serial Python (Amdahl). On the deep library the first sweep — and the miss sweep, which is the
  same work at K = 233 — is 91 % native: 2.3 s of 24.7 s is Python, so the first sweep is at best 2.3 + 22.4/T =
  5.1 s at T = 8 (4.8×) and the miss sweep ~14 s from ~59 — the two layer-running sweeps are where this form pays.
  A cache-served refit sweep is 67 % native: 7.7 s of 23.3 s is Python, so at best 9.6 s (2.4×), and once (iv-a)
  has taken ψ to ~2 s it is ~8 s of Python that threads over blocks cannot touch — the per-block items below
  are what shrinks it.
* ONE NATIVE CALL PER BLOCK — `solve_block` in C++: the self-solve ψ, the builders, the two passes, the solve,
  the final ψ, the write-back — and a C++ pool over blocks. The Python per block is then the slicing and the
  loop (~1 ms). This is the block-in-native end-state the owner named; it subsumes the items below and is the
  larger step.

The Python per block that must shrink first — each EXACT (bit-identical, or inside the derived budget where
stated), each its own commit, in the order of its size:

1. **The cache key on the factory's INPUTS, not its rows** — LANDED 2026-09-18 (3.2 s → 0.27 s a refit sweep on VCaP). The rows are a
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

### (v) The blur's reduction — a summation-order change, priced by the census

The pass is the blur: 1.43 M calls a first sweep at a mean of 56 taps over K = 101 cells (8.1 G multiply-adds of
the pass's 10.4 s; ~2.3× at the miss sweep's K = 233). `blur_row`'s inner loop is a floating-point reduction over
the taps, which the compiler does not reassociate, so it runs scalar at about a nanosecond a tap. An explicit
four-lane partial sum changes the summation order and nothing else — a port-class change, held to the replay's
tolerance budget, worth ~−20 s a run on its own and multiplying with a pass threaded over blocks. The constant
edge pads are prefix sums, a smaller saving. A recursive (IIR) Gaussian would change the blur's VALUE, not its
summation order — an accuracy question, not proposed.

### What is not proposed

* Threads inside a pass: a hop reads the row the previous hop wrote, in chain order — sequential by
  construction. The builders are per node but 1.6 s a first sweep only; parallel over blocks or not at all.
* OpenMP: the tree links none; the EM's pool is the idiom.
* The underflow cut in ψ — skipping the cells whose log-weight lies below the truncation constant under the
  cube's maximum, a plausible 2× on ψ at a tolerance — not before (iv-a) is measured; a separate A/B if at all.
* An arena for the per-block tables: priced at ≤ 0.3 s of a run's first sweep (the first-touch faults the
  allocation step left; the refit sweeps allocate none of these tables under the cache) — not built.

## The order proposed

1. (iv-a) ψ over slots — one commit; the replay at 1/2/8 threads BIT-IDENTICAL on `sweeps_VCaP_step13`, the suite,
   the references; two interleaved pairs on VCaP at `--threads 8`. (58 s of 264.)
2. (v) the blur's four-lane reduction — one commit, tolerance-gated on the replay; the census re-run to show the
   taps unchanged and the time moved. (~36 s of 264, the pass.)
3. The cache key on the factory's inputs (exact; the largest Python item of a refit sweep).
4. The gDNA arm inside ψ (tolerance-gated); `_psi`'s row add with it.
5. The block in one native call with a C++ pool over blocks — or Python threads over blocks with GIL-releasing
   kernels, which on the deep library's layer-running sweeps already pays — judged on the Python floor measured
   after 3–4.

## Open for the owner

1. The thread budget's home: a `CalibrationConfig.n_threads` mirroring `EmConfig.n_threads`, fed by the
   pipeline's `--threads`?
2. The cache key digesting the factory's inputs instead of its rows.
3. After (iv-a): the block in one native call (the end-state, the larger step), or Python threads over blocks
   (sooner, bounded by the Python)?
