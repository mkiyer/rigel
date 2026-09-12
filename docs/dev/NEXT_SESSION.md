# NEXT SESSION — the C/C++ port of the block solve (2026-09-11, after the locus sweep landed)

A handoff, provisional like everything in this directory. The references are `CLAUDE.md`, `docs/ROADMAP.md`,
`docs/DESIGN.md` §6b.15 (the ruling and every measured fact behind it) and `docs/ISSUES.md`
(`performance-memory-bounded-solve`).

## What landed (uncommitted on `main`; the owner drives commits)

Three steps, each gated before the next, in this order:

1. **The chain declares its terminals; the backbone enforces them.** `region_chain.locus_blocks(chain,
   terminal, block_slots)` cuts the chain at every terminal (`is_region & g1_locked`) and reference start
   and merges loci up to a block size; `sweep._pass` never asks a kernel for a hop INTO a terminal; the
   gDNA lane no longer lists faces FROM one (65,852 dead faces on the human chain). Bit-identical on the
   four captured sweeps and the three frozen references of the previous session.
2. **ψ's read-out is chunk-exact** (`simplex_logodds`: a contiguous ψ, `_row_moment` per-row sums instead
   of BLAS; `density_deconv.density_factor_precision` likewise). The one step that moved a number:
   ≤ 3.1e-15 per slot per sweep, no amplification through four sweeps and three refits, TPM and effective
   lengths bit-identical on LBX0190, `calibration_vs_oracle.py` at the last ulp on all 16 ladder rows with
   `ruler_n_moved` unchanged. Accepted by the owner as identical to a tolerance. The references were then
   re-frozen from this tree: `~/Downloads/rigel_runs/arms/locus_identity_*.json` and the captures
   `~/Downloads/rigel_runs/perf/sweeps_MO_3021_step2` (the old `sweeps_MO_3021` and `cleanup_identity_*` /
   `perf_identity_LBX0190.json` now differ from the tree by that priced ≤ 1e-14).
3. **The sweep solves a locus block at a time.** `solve_chain` is an orchestrator: the library reduction
   once (`Policy.library(ChainView)` — observations and geometry, no beliefs), then `_solve_block` per
   block on `_slots(...)` views of every input (a chain view re-bases its links), the owned prefix written
   back, `_gather` re-keying the diagnostic capture to the chain. `TransferPolicy(strand=)` reads the
   factory rows from the context (`factory_rows`) and the strand deadband's verdict from
   `ChainView.strand_live` (`region_init.strand_discriminability`). `CalibrationConfig.sweep_block_slots`
   (CLI `--sweep-block-slots`) is the working-set knob. Bit-identical to step 2's captures and references
   for the whole chain and for eight block sizes on a real 2.09M-slot sweep.

Suite: 3,414 passed / 2 xfail / 3,416 collected (the ten gates are in existing files). `preflight --full`
was run at the end of the session — re-run it first thing.

## The numbers that decide what comes next

* One real sweep (876k-fragment library, 2.09M slots) is ~40 s: `prepare` ~17 s, own claims ~4 s, the two
  passes ~10 s, the final ψ ~4 s, the policy's solve ~2 s. All of it Python, one core.
* Threads on the locus-split passes: 0.94× / 0.83× at 8 threads (GIL-bound); ψ's grid solve 2.06×; the
  same passes in 8 forked processes: 6.16×. Threads are not the tool; the executor is the port's.
* The whole-chain sweep's own peak allocation is 9.9 GB; per-block figures are in the config docstring.
* The chain: 33,120 terminals, 33,018 loci, median 19 slots, largest 2,477 (0.12 %) — a serial floor of
  ~840× that only a compiled block solve can cash.

## The agreed order (owner, 2026-09-11) — and where it stands

⓪ re-measure the deep library end to end (`main` vs this tree, two back-to-back pairs at 8 threads,
`profiler.py --compare`) — DONE: peak 33.2 → 14.9 GB, wall 0.97–0.98, sweeps 0.96–0.97, the reports in
`~/Downloads/rigel_runs/perf/ab_locus_2026-09-11/`; the peak is now `build_region_geometry`'s transient
and the pre-sweep `init_beliefs` solve, not the sweep; ① DONE — the whole message layer is
refit-invariant given the grid (proved on the captures: sweeps 1–3 deliver identical rows, cube rows and
hearings), so `calibrate` shares a content-keyed `sweep.MessageMemo` across the refit sweeps; the context
now carries `own_live` (the liveness bits) instead of the self-solve object; on the deep library the run
reads 0.65 of its wall (refits 2–3 served in 38 s each instead of 176 s) for +2.8 GB of peak (the memo
holds 2.68 GB, cube rows as float32; reports in `~/Downloads/rigel_runs/perf/ab_memo_2026-09-11/`) —
whether it should be switchable is the owner's; ② DONE — `messages.transfer.Faces`, the rules as
``(n, 2)`` typed tables over (destination, side) with five kinds (FORWARD, TRANSPORT, SPLICE_OUT, EDGE,
LEVEL) and a row store, the lanes' faces and two-sided faces as ``(n, 2)`` bits, the junction flux as a
row table; `_neighbour_pairs`, the closures and the face sets are gone; `Faces.apply` is the one home of
the rule arithmetic and what the port's pass reads; ③ the port of `_solve_block` behind a derived
tolerance gate; ④ the factory rows per block; ⑤ the scan and the second pass. `ISSUES:
performance-memory-bounded-solve` carries the same list with the reasoning; `ROADMAP.md` rank 1 the order.

## What is next — the port (owner, 2026-09-11)

`_solve_block` is the unit to port: one block's own claims, the policy's claims and rules, two passes, one
solve, one write-back, on slices of the per-slot arrays plus the library scalars. Its inputs are exactly
`_view_fields` + the beliefs + `factory_rows` + `scalars` + `library`; its outputs six belief arrays,
`informed` and the assertion counts. Parallelise inside the port over blocks. Keep every gate: the replay
(`sweep_replay.py replay --block-slots`), `rename_identity.py --check` against the `locus_identity_*`
references, the suite, `preflight --full`. Still genome-wide and worth folding into the block: the
memoised intron-factory rows `(n_slots, K)` (1–2 GB per grid) — a block can build its own from the
substrate slice.

## Decisions recorded

* The refit count stays (owner). The scan's thread split is a separate decision
  (`ISSUES: scan-thread-split-starves-the-workers`).
* The thread count reuses the existing `--threads`; the block size is a performance tunable with a measured
  default (`CalibrationConfig.sweep_block_slots`).
* Do not re-derive: a per-row NumPy solve is not chunk-exact by default (F-ordered fancy-index results,
  BLAS at one row); the gate `test_the_psi_solve_is_chunk_exact_so_a_block_split_moves_no_number` holds it.
* When patching `calibrate` for an arm or a spy, take the module from
  `importlib.import_module("rigel.calibration.calibrate")` — `import rigel.calibration.calibrate as C`
  binds the package's re-exported FUNCTION, and a patch on it is an ablation that never ran (it cost one
  A/B pair this session; `TRAPS: could-the-arm-have-fired` is the check that caught it: the layer's
  stages read 1.00 between "arms").
