# The block in one native call — the design (2026-09-18), and the three commits that land it — ALL THREE LANDED 2026-09-18 (snapshots 18–20 awaiting the go)

The owner's decision of 2026-09-18: the block in ONE native call, a C++ pool over blocks, preferred over Python
threads over blocks. The rulings it stands on: `DESIGN.md` §6b.15.1–§6b.15.3 (the locus block is the unit of the
solve, the only information that crosses a block boundary is the policy's LIBRARY, ψ's read-out is chunk-exact),
§6b.15.9 (threads were the wrong tool for the Python sweep; the parallelism waits for the port — this is that port's
last step), the one-path ruling of 2026-09-17 (one production code path; the Python a validated kernel replaces is
deleted; a small floating-point tolerance is accepted for speed), and the working rules: one mechanism per commit,
each gated as `THREADS_PLAN.md`'s protocol says. This file is the design put on paper before anything is built; the
record of what landed is `DESIGN.md` §6b.15.5 and `ISSUES: performance-memory-bounded-solve`. As built: VCaP at 8 threads 198 → 144 s and 194 → 138 s (0.73 / 0.71), calibrate 85.6 → 33.7 s, the four sweeps 68.8 → 16.8 s and 66.3 → 15.9 s (0.24).

## 1. Where the time is, and what a native block buys

VCaP (18.6 M fragments) at 8 threads after the splice-out hoist (`perf/splice_out_2026-09-18/pair1_post.json`): the run
195 s, calibrate 89 s, the sweep 73 s. Inside the sweep's block loop (68 s over the four sweeps): the builders 17.1 s,
the pass 29.2 s, the policy's solve 1.9 s, the two ψ solves 7.7 s (threaded over slots since step 15), the backbone's
checks 0.2 s — and ~11 s of Python between the kernels: the block's slicing and context, the four allocations, the
cache's key, the gDNA arm's interpolation (`np.interp` over ``(n, K)``, 2.3 s a refit sweep), the factory rows (lgamma,
0.5–1.1 s a sweep), the row add, the write-back. Every kernel runs one block at a time on one thread, so 48 s of the 68
are serial native work that the blocks' independence divides: 426 blocks, none sharing a byte of output, the same
arithmetic per block whatever thread takes it. Expected after this design lands: the layer ~48 → ~7 s, the Python
~11 → ~1 s, ψ's wall unchanged (already threaded), the sweep ~73 → ~17 s, the run ~195 → ~140 s at 8 threads.

## 2. The target — what the tree looks like after

**The call.** `native.solve_blocks(...)`: ONE call per sweep. It takes the chain's arrays whole (the observations and
the geometry `ChainView` already holds, the incoming belief, the terminal bits), the block table (`start, stop, end`
per block from `region_chain.locus_blocks`), the solve's scalars, the policy (its name, its strand model, its library),
the priors as their INPUTS (the landscape's curve; the factory's background, intron mask, counts and opportunities),
the cache's served deliveries per block, and the thread budget; it writes the belief arrays and `has_composition` in
place and returns the assertion counts, the deliveries of the blocks whose layer ran (for the cache), and — asked for —
the diagnostics capture. Inside it, a pool of threads pulls the blocks one at a time; each block runs the whole pipeline
on its own thread and its own arena:

1. the prior rows: the gDNA arm interpolated from the landscape's curve at ``log f_g + log M − log E`` per slot per
   cell; the factory rows ``log NegBinom(f_g·C; ρ_bg·E, α_eff)`` at the intron slots — both built in the arena, never
   crossing to Python;
2. the SELF-SOLVE ψ (every slot with a live strand and a fragment, its own cube, the incoming belief the variance
   freeze) → ``fg_loc``;
3. ``tau_lam`` = the single-strand strand evidence at ``fg_loc`` + the factory rows' curvature → ``has_own_composition``
   (``> 0``) and the own-evidence predicate (``> 1e-9``, the one constant, read by the instruments' predicate);
4. the LAYER, unless the policy is silent or the cache served the block: the builders into the arena's tables (the
   claims, the faces, the three lanes), the forward pass, the backward pass (the two received tables in the arena), the
   solve → the fused λ rows and the cube delivery; the backbone's checks on what was delivered;
5. the FINAL ψ: the arm, the factory row and the delivered row added per cell in the kernel (no ``(n, K)`` add), the
   cube rows at the AMBIG slots;
6. the write-back on the owned slots (a solvable slot takes its solve, a locked or empty one keeps the incoming belief),
   ``has_composition`` = own evidence | structurally certain | a composition held from either side;
7. the counts of the four assertions on the owned slots — ``population_at_most_three``, ``population_reaches_three``,
   ``lam_rows_finite`` (where rows were delivered), ``cube_rows_finite`` (where a cube was), ``writeback_only_solvable``
   — as integers per block; the backbone sums them and RAISES exactly as `AssertionCounts.note` does today.

**Bit-identity by construction.** The kernels are the ones that run today (`transfer_kernel`'s builders, pass and solve;
`psi_kernel`'s `solve_slot`), called on the same per-block numbers in the same order; a block's outputs are a disjoint
slice of the chain's; nothing is reduced across blocks but the integer counts (a sum). So the call at any thread count is
bit-identical to today's serial loop — PROVIDED the two prior constructions the kernel absorbs are exact, which is what
commits 18 and 19 below establish first.

**Threads and memory.** The pool is the module's one `EStepThreadPool` (persistent, rebuilt on a budget change, the GIL
released around the call; `resolve_threads(n_threads, n_blocks)`). One block at a time by an atomic counter — no chunk
constant — balances the M3's asymmetric cores; a block's ψ runs serially on its thread (the same total work as the
per-slot pool, the same bits). Each thread owns an arena: the ``(n, K)`` tables of one block — the claims, the face row
store, the three lanes' own levels, the flux store, the two received tables (a composition and three level profiles
each), the arm, the factory rows, the fused rows — about 16 matrices of 5,000 × 233 doubles at the refit bracket,
~150 MB, resized per block and freed when the call returns: ~1.2 GB at 8 threads during a layer-running sweep, under
quant's 11 GB peak. The flux levels are kept in a store indexed per ``(node, side)`` like the face rows (a ``(n, 2, K)``
matrix per lane was 6 mostly-empty matrices).

**The Python that remains** (layer 6): `sweep.solve_chain` — the structure and the chain view once, the library, the
blocks, the cache's keys and lookups, the one call, the counts, the cache's puts, the capture (~150 lines against 770);
`messages/` — the two policies as a NAME, a strand model and a library reduction (`TransferPolicy.library` is the only
cross-block computation and stays numpy: three ratios of sums over structurally selected slots and one boolean);
`message_cache.MessageCache` — keyed on the block's slices of the very arrays the kernel reads; `blocks.SweepCapture` —
the diagnostics record and its gather; `region_init` — the strand protocol decision and the own-evidence predicate;
`density_deconv` — the background fits; `landscape` — the curve's fit; `simplex_logodds` — ψ's dispatcher for the
pre-sweep solve and the gates' read-outs, and `CubeRows`.

**What is deleted** (every reader accounted for by the census of 2026-09-18): `messages/faces.py` (`Faces`, `RowTable`,
`side_of`), `messages/lanes.py` (`LevelLane`); in `messages/__init__.py` `Received`, `Levels`, `PsiMessage`,
`BlockContext`, `Prepared`; in `messages/transfer.py` `_prepare`, `_Chain`, `_PreparedTransfer`, `MARGINAL_NODES`
(the nine equal-probability nodes move to the kernel's constants — the quadrature resolution of the splice-out marginal,
like ``n_tilt``), `read_column`; in `sweep.py` `_solve_block`, `_message_layer`, `_psi`, `_pass`, `_check_message`,
`_write_back`, `_block_diagnostics`, `_gdna_logprior`, `_Sweep`, `_RowsOfArray`, `_factory_of`; in `blocks.py`
`block_slice`; in `region_init.py` `build_region_init`, `strand_evidence`, `RegionInit`; in `density_deconv.py`
`density_lambda_factor`, `_log_negbinom`, `density_factor_precision`; in `landscape.py` `DensityLandscape.logprior`;
in `calibrate.py` `FactoryRows.__getitem__` / `__array__` (the factory keeps its inputs and its digest). The
`Prepared` protocol goes with them: a message policy is no longer a Python object the backbone calls per hop — it is a
name the kernel switches on (``silent`` runs no layer), a strand model, and a library. `_layers.py` drops the two
modules.

**The native module.** ONE module for the solve's kernels, `_solve_impl`, built from `native/solve_kernel.cpp` — the
block pipeline, the pool, every binding — which includes `transfer_kernel.h` (the builders, the pass, the solve, on
pointer views that the arena and the gates' numpy arrays both satisfy; today's `.cpp` with its bindings removed),
`psi_kernel.h` (ψ's `Grid`, `SlotInputs`, `solve_slot`, `Scratch`; today's `.cpp` likewise), `transfer_rows.h` and
`thread_pool.h`. `native.py` exports `solve_blocks`, `psi_solve`, the ψ read-outs for the gates, the transfer kernels
for the gates and `transfer_rows` as before. `_transfer_impl` and `_psi_impl` are gone; `CMakeLists.txt` builds one
target.

**The cache.** Kept (the owner's switch, `DESIGN.md` §6b.15.4). The key is a blake2b over the block's slices of the
arrays the kernel reads — every array field of the `ChainView`, the incoming belief's ``f_g`` — the grid's two scalars
and the strand liveness, the factory's digest for the block (`FactoryRows.digest`; rows given as one array digest by
content), the library and the policy's name and strand. ``has_own_composition`` leaves the key because it is a function
of what the key already holds (a count above zero, the protocol's liveness, the single-strand bits, the factory rows'
curvature); the incoming belief stays because the own strand claims read it. An entry is the block's delivery as the
kernel returns it: the delivered rows' local slots and the rows (sparse by construction — the solve knows which rows it
wrote), the cube table, the owned slots' held-composition bits. A served block hands the entry back to the kernel, which
skips the layer, counts the served rows' finiteness itself, and runs ψ. PRICED FOR THE OWNER after the timing: what the
cache saves at 8 threads (two served sweeps' layers) against what it costs (the memory of the entries at K = 233, the
digest per block, `message_cache.py`, the factory's digest, five gates).

**The gates.** One implementation, read through bindings that allocate and return fresh arrays: `transfer_prepare`
(one block's context in → the tables out, as a dict of arrays), `transfer_pass` (the tables and a received table in,
one pass or one hop written in place), `transfer_solve` (the two received tables in → the rows and the cube out), and
in `transfer_rows` the absorbed constructions — `factory_rows`, `factor_precision`, `strand_evidence`, `gdna_arm`,
`lgamma` — each held to its analytic properties and to scipy at a tolerance where scipy is the oracle. The containers
the gates read the tables through — `RowTable`, `Faces`, `LevelLane`, `Received`, `Levels` — move VERBATIM from `src/`
to `tests/calibration/_transfer_harness.py`: they hold arrays and index them, nothing more, and no production reader is
left. The backbone gates that today drive fake Python policies (`_Echo`, `_Quiet`, `_Full`) are rewritten against the
one backbone: the pass order and sides, the terminal rule and the two states are read off the kernel's received tables
on hand-built chains; the perturbation gates on `AssertionCounts` stand; the shape refusals that were Python checks are
structural now and go. `test_region_init`'s `build_region_init` and `strand_evidence` gates, `test_density_deconv`'s
factor gates and `test_landscape`'s two `logprior` gates read the bindings.

**The instruments.** `SweepCapture` keeps every field an instrument reads (the census: `fg_loc`, `fg_strand`, `tau_lam`,
the observations, `mass_global`/`eff_global`, `policy_name`, `solve_grid`, `backbone_assertions`) and gains `tau_fac`
— the factory rows' precision, which `solvability_audit.py` recomputed from the captured rows through
`density_factor_precision`; `intron_prior` (the chain-wide rows) leaves the capture with its last reader. The two
received tables it publishes become dicts of the kernel's arrays. `profiler.py`'s per-block rows (`_solve_block`,
`build_region_init`, `prepare`, `_pass`, `solve`, `_check_message`) go; the call itself is a row.
`sweep_replay.py` needs nothing: `solve_chain`'s signature stands and the captured `TransferPolicy` / `FactoryRows` /
`DensityLandscape` / `MessageCache` unpickle into the same classes (the captured cache's entries, keyed the old way,
miss — as they did at step 16).

**`policy_prototype.py`** installs a Python class in place of `TransferPolicy` for the `transfer` arm. With the layer in
the kernel there is no Python hop for a prototype to implement: the instrument's `--module` mechanism cannot work, and an
arm that silently ran the shipped policy under a prototype's name is exactly what CLAUDE.md forbids. In commit 20 the
instrument REFUSES a `--module` arm with the reason and keeps scoring the two shipped policies; whether it is retired, or
re-pointed at prototypes built in C++ and swapped in as the census instrument is, is the owner's call (§6).

## 3. The three commits

Each its own snapshot in `commits/`, each gated, in this order, so that the big one is a pure restructure.

**18 — ψ takes its priors apart** (`THREADS_PLAN.md` items 2–3). `psi_solve` / `psi_cube` take the gDNA arm as the
landscape's CURVE ``(log_rho, logP)`` with the per-slot support ``(mass, eff)`` and interpolate at each cell themselves
(`transfer_rows::interp`, numpy's end rules); the λ-factor row and the delivered row arrive as two inputs and are added
per cell. `sweep._gdna_logprior`, `DensityLandscape.logprior` and `_psi`'s ``(n, K)`` add are deleted; the dispatcher's
``gdna_logprior`` becomes ``gdna_prior=(log_rho, logP), gdna_support=(mass, eff)`` (arrays, so layer 3 imports nothing
from layer 5), its ``lam_logprior`` gains ``row_logprior`` beside it. GATE: the replay on `sweeps_VCaP_step16` — expected
BIT-IDENTICAL (the kernel's `interp` is numpy's formula in numpy's order and its `sigmoid` is scipy's `expit` to the bit
on every solve grid, checked 2026-09-18); else `--tolerance` and the moves recorded — the suite (the two `logprior` gates
re-pointed at `transfer_rows.gdna_arm`), the three references. No timing pairs: an enabling step (~2 s a refit sweep).

**19 — the factory's log-gamma is the kernel's.** `_log_negbinom` reads `transfer_rows.lgamma` (libm) instead of
scipy's `gammaln` (cephes). Nothing else moves. GATE: the replay with `--tolerance` (the last-bit moves of lgamma, the
budget beside them), the suite (the negative-binomial gates against `scipy.stats.nbinom` at their existing tolerance),
the three references — which DIFFER at the ulp and are RE-FROZEN with this reason. This is the one commit in the series
that moves a number, and it is a two-line change, so the move is attributable.

**20 — the block in one native call.** Everything in §2. GATE: the replay on the four VCaP sweeps at 8 threads AND at
1 thread, BIT-IDENTICAL to the tree at 19 (re-captured there if 18 or 19 moved a number); the suite, re-derived; the
three references BIT-IDENTICAL against the re-frozen set; `preflight.py --full`; two interleaved timing pairs on VCaP at
8 threads against a worktree of 19 carrying its own modules (the argument shapes change), the peak RSS read beside the
wall. The perturbation half: break the arena's reuse between blocks (a stale row) and watch the replay differ; skip a
block in the pool and watch it; poison the arena's tables at allocation and watch the poison gate (kept, on the
bindings' allocate-and-return tables).

## 4. The order inside commit 20

1. `transfer_kernel.cpp` → `transfer_kernel.h`: the views (`Chain`, `RowsOut`, `FacesOut`, `LaneOut`, `FacesView`,
   `LaneView`, `LevelsView`, `Held`, `SolveLane`) on raw pointers with their sizes; `prepare_block`, `pass_block`,
   `solve_block` as functions of those views; the flux store indexed per face. `psi_kernel.cpp` → `psi_kernel.h`.
2. `solve_kernel.cpp`: the `Arena` (per thread), the block pipeline, `solve_blocks` with the pool, the bindings —
   `solve_blocks`, `psi_solve` and the read-outs, `transfer_prepare` / `transfer_pass` / `transfer_solve` in their
   allocate-and-return forms, the `rows` submodule with the absorbed constructions. `CMakeLists.txt`, `native.py`.
3. Python: `sweep.py` rewritten around the call; `messages/` reduced; `message_cache.py` re-keyed; `blocks.py`,
   `region_init.py`, `density_deconv.py`, `landscape.py`, `calibrate.py` trimmed; `_layers.py`.
4. The gates: the harness takes the containers; the transfer gates re-pointed; the backbone gates rewritten; the
   absorbed functions' gates on the bindings; the poison gate on the bindings' tables; a thread gate (`solve_chain` at
   1, 2 and every core bit-identical on the chunk gate's substrate).
5. The instruments and scripts: the capture's fields, `profiler.py`'s rows, `policy_prototype.py`'s refusal,
   `scripts/README.md`.
6. The gates in order: build; the suite (derive the count); the replay ×4 at 8 and 1 threads; the references; the
   perturbations; the timing pairs; `preflight --full`; the docs (`DESIGN.md` §6b.15.5, `ISSUES`, `CLAUDE.md`'s message
   layer row and baseline, `TESTING.md` §6, `THREADS_PLAN.md`, `NEXT_SESSION.md`); the snapshot.

## 5. Risks, and what answers each

* A last-bit difference in the arm's interpolation or the factory's arithmetic against numpy — caught in 18 and 19,
  where each is alone; 20 then owes bits. The strand evidence and the factor precision reproduce numpy's expression
  order, the precision's row sums numpy's PAIRWISE summation (its blocked eight-accumulator scheme), so `tau_lam` is
  the same bits and the predicate cannot flip.
* Returning per-block arrays from worker threads: the workers fill `std::vector`s; the calling thread wraps them as
  numpy arrays (a capsule owner, zero-copy) after the pool returns, under the GIL. Every input pointer is taken before
  the GIL is released; an exception in a worker is stored and rethrown after the join.
* The arena's memory at 8 threads (~1.2 GB) — read on the timing pairs' peak; under quant's peak by design.
* The captured `TransferPolicy` pickles: the class keeps its module and its one slot, so old captures replay.

## 6. Open for the owner (after 20 lands)

1. The message cache: keep or delete, at the price §2 records after the timing.
2. `policy_prototype.py`: retire, or re-point at C++ prototypes swapped in like the census instrument.
3. The layer-4 `strand_likelihood` executable reference (unchanged here) and the blur's loop interchange (priced at
   step 17, not taken) — both still the owner's.
