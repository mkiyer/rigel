# NEXT SESSION — start here (2026-09-17 night: the one-path convergence and the VCaP baseline; the threads design awaits the owner)

This file is only how to begin. The port's plan is `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` §F and §G; the LIVE
status with numbers is `ISSUES: performance-memory-bounded-solve`; the day's record is the closing paragraphs of
`DESIGN.md` §6b.15.5; ψ's design is `docs/dev/PSI_PORT_PLAN.md`; THE THREADS DESIGN, put to the owner before
building, is `docs/dev/THREADS_PLAN.md`.

## The order (owner, 2026-09-17: the release ships the contraction as it stands; performance until the tool is fast;
## finish the block in native BEFORE threads; design threads on paper first)

1. ~~The ruler's repair, the yield's floors~~ — DONE 2026-09-16/17.
2. ~~The port, steps (i)–(iii): the pass, the builders, ψ~~ — DONE 2026-09-17, with the one-path cleanup.
3. ~~The block in native, step 1: the policy's solve~~ — DONE 2026-09-17 (`28bd174c`).
4. ~~The block in native, step 2: the tables' allocations~~ — DONE 2026-09-17, PREPARED AS A SNAPSHOT for the owner's
   go (`commits/12_alloc/` in the synced scratchpad: the files, `MESSAGE.txt`; the tree at `28bd174c` + these edits IS
   the snapshot). The factory rows were priced and not ported alone (`ISSUES: performance-memory-bounded-solve`).
5. ~~Step (iv-a), ψ over slots~~ — DONE 2026-09-18 with the owner's answers (`CalibrationConfig.n_threads` fed by
   `--threads`; the cache key on the factory's inputs; the block in ONE native call preferred over Python threads);
   VCaP 267 → 214 s and 261 → 213 s at 8 threads (0.80 / 0.82), ψ 57.4 → 8.3 s, the sweep 135 → 88 s; snapshot `commits/15_psi_threads/`. As built: the EM's pool
   (`native/thread_pool.h`), a cost-balanced partition of the slot list (26 columns an AMBIG slot, 1 a single-strand
   one), a `Scratch` per thread, the GIL released; gate: the replay of `sweeps_MO_3021_step11` at 1, 2 and 8 threads
   BIT-IDENTICAL on all four sweeps, the suite, the three references; two interleaved pairs on VCaP at `--threads 8`.
   Then (v) the blur's four-lane reduction (the pass IS the blur on the deep library: 1.43 M blurs a first sweep,
   8.1 G multiply-adds; a summation-order change behind the replay's tolerance budget), the cache key on the
   factory's inputs, the gDNA arm inside ψ, and the block in one native call — in that order, each its own
   commit, each gated as the plan says.
6. ~~The one-path convergence of the pass~~ — DONE 2026-09-17, prepared as snapshot `commits/13_one_path/`
   (the per-hop Python kernel and `transfer_rows.py` deleted; the constructors bound for the gates).
7. ~~The optimisation target moved to VCaP~~ — DONE 2026-09-17 (owner): the four sweeps captured
   (`perf/sweeps_VCaP_step13`), the whole run profiled twice (`perf/vcap_baseline_2026-09-17/`), the sweeps dissected
   with the census and the kernel wrappers, the ranked opportunities in `ISSUES: performance-memory-bounded-solve`;
   the MO_3021 captures and the timing worktree removed; snapshot `commits/14_vcap_baseline/` (docs only).
8. §G: the stages outside calibration — 110 s of 264 on VCaP: the scan 35, the second pass 22, the two
   fragment-length fits 16.6, quant 37 (the locus EM 15, the capture effective lengths 8), the index load 6.5.

## Where the time is now (VCaP, the deep library; `s12/kernel_times_vcap.log`, `s12/census_vcap.log`, `perf/vcap_baseline_2026-09-17/`)

The whole run at 8 threads: 264 s (two back-to-back runs 263.6 / 265.5, the stages drifting 0.98–1.07). Calibrate 154 s
(the sweep 134: ψ 54, the pass 36.5, the builders 17.2, the solve 2.0, Python ~24; the landscape fit 4.2; ψ before the
sweep 3.8); outside calibration 110 s (the scan 35, the second pass 22, the fragment-length models 8.3 + 8.3, quant 37, the
index load 6.5).

| per sweep | first sweep (K = 101, the layer runs) | refit sweep (K = 233, cache-served) |
|---|---|---|
| ψ native (self-solve + final) | 6.7 s | 15.5 s |
| the pass native (both): 1.43 M blurs at a mean of 56 taps — 8.1 G multiply-adds | 10.4 s | — (~24 s at the miss) |
| the builders native (the RNA lanes 2.7, the gDNA lane 1.2) | 4.7 s | — (~11 s at the miss) |
| the policy's solve native | 0.5 s | — |
| the message cache's key (blake2b over the block context) | — | 3.2 s |
| the gDNA arm (`np.interp` over (m, K)) | — | 2.3 s |
| the factory rows (lgamma 85 %) | 0.5 s | 1.1 s |
| the glue | ~1.8 s | ~1.2 s |
| the sweep | 24.7 s (native 91 %) | 23.3 s (native 67 %) |

⛔ The refit sweeps run on the landscape's bracket — K = 233 on VCaP against the first sweep's K = 101 — so every per-cell
cost is 2.3× there for the same slots, and sweep 1 misses the cache in production and runs the whole layer at that K;
profile `replay --call 1` beside `--call 0`. ⛔ The census is a SCRATCH instrument (`s12/census/`: a patched copy of
`transfer_kernel.cpp` and `transfer_rows.h` built as `_transfer_census` and swapped in for the three native calls by
`s12/census_run.py`); rebuild it from the current kernel sources before trusting it after any kernel change.

## The one-path ruling (owner, 2026-09-17) — what is still duplicated

One production code path; once native code is validated the Python it replaces is deleted; a small floating-point
tolerance is accepted. Applied to the builders, ψ, the policy's solve and now the pass: the per-hop Python kernel
(`_PreparedTransfer.propagate`, `Faces.apply`, `LevelLane.emit` / `receive`) and `messages/transfer_rows.py` are
deleted; the row constructors and flag predicates the gates need are the native ones, bound as
`native.transfer_rows`; a rule or a lane hop is gated by driving ONE hop of the native pass
(`_transfer_harness._hop` / `_rule`). What remains duplicated: the layer-4 `strand_likelihood` "executable
reference" module (its gate reads the native strand term through `psi_cube`) — the owner's call.

## The protocol for every step (unchanged)

1. `python scripts/design/preflight.py` first; the suite's standing count is in `CLAUDE.md`.
2. The capture: `perf/sweeps_VCaP_step13` (9.0 GB; the deep library, VCaP, taken from this tree at 8 threads, ~4 min)
   replays BIT-IDENTICAL — use it; take a fresh one only after a commit that moves numbers, and delete the superseded
   one once its successor replays. The MO_3021 captures are gone: the optimisation target is VCaP (owner, 2026-09-17).
3. DERIVE from the Python → the C++ in-tree → a scratchpad harness that runs both on every call of a captured sweep
   against the replay's derived budget → wire, delete the Python, rewrite the gates → the suite →
   `sweep_replay.py replay [--tolerance]` for calls 0–3 → the references (re-frozen only with the reason) → two
   interleaved timing pairs at 8 threads on VCaP against a worktree of the pre-step commit carrying the SAME
   binaries (`git worktree add --detach /tmp/rigel_pre <commit>`, the installed `site-packages/rigel/*.so` copied in,
   `pre_site_28bd174c/sitecustomize.py` on `PYTHONPATH`; `s10/time_pairs.sh`) — the worktree was removed after the
   allocation step's timing and is recreated when needed.
4. Each step its own commit, prepared as a snapshot; the owner drives the commit and the push.

## ⛔ Traps met in this session (also in memory)

* A `sed` revert with an unescaped `*` or `[i]` in its PATTERN silently matches nothing: two perturbations stayed in
  the source, the "clean" rebuild was doubly broken, and its binary was copied to the timing worktree. Perturb and
  revert with exact-string Python edits and `git checkout -- FILE`; assert the anchor count; check `git diff` after
  every revert.
* The timing worktree must carry the binaries the working tree IMPORTS (`site-packages/rigel/*.so`), not the
  `build/` directory's — the scanner's differed there.
* A zero-fill's cost does not vanish when the fill does: the allocator returns freed pages and the first touch in
  the kernel pays the fault — read the kernels' times beside the allocation line, never the allocation line alone.
* The refit sweeps are a different problem from the first (K = 202 against 101): a mechanism found on `--call 0`
  can be half the story.
* A `np.empty` table needs a gate that a `np.zeros` table never needed: the poison gate (NaN in every unmasked cell,
  outputs bit-identical) — and any gate that compared whole matrices (`test_pass_kernel`'s leaves, the `RowTable`
  gate's fresh-table check) had to be rewritten to compare under the bits.

## Storage (2026-09-17 night)

`perf/sweeps_VCaP_step13` (9.0 GB; the refit sweeps' pickles are 2.9 GB each, their message cache) is THIS tree's
capture and stays until a step moves numbers. Removed tonight: the MO_3021 capture (3.6 GB), the timing worktree
`/tmp/rigel_pre` (recreate from the pre-step commit with the installed binaries when a step is timed), the retired
modules' binaries (`s9/pre_so`). Still large and regenerable: `prototypes/2026-09-16_ruler_repair/s5/scratch_test_reference*`
(22 GB) — the owner's.

## Where everything is

* The synced scratchpad: `~/Downloads/rigel_runs/prototypes/2026-09-17_port_ii/` — `commits/committed/` (7–11),
  `commits/12_alloc/`, `commits/13_one_path/` (with its `DELETED.txt`), `commits/14_vcap_baseline/` (the snapshots
  and messages awaiting the go, in order), `s8/`, `s9/`, `s10/` (the allocation step), `s11/` (the convergence: the
  edit scripts, the replay and identity logs), `s12/` (the VCaP baseline: the census instrument, its harness, the
  kernel timings, the profiles' logs), `pre_site_28bd174c/`.
* Captures and reports: `perf/sweeps_VCaP_step13`; `perf/vcap_baseline_2026-09-17/` (the two runs of this tree);
  the earlier `perf/block_alloc_2026-09-17/`, `perf/port_*_2026-09-17/`.
* The identity references: `~/Downloads/rigel_runs/arms/review_identity_*.json`, frozen 2026-09-17 on `80b459bf`'s
  numbers, BIT-IDENTICAL on this tree.

## Decisions on record (unchanged, carried)

* Float64 for the whole of ψ; ONE solver; ONE λ lattice at a dimensionless step with a stated guarantee; no θ
  lattice anywhere; the tilt's hypothesis space is {pure +, pure −, mixed}.
* The strand channel's liveness is a protocol decision on the spliced 2×2; gDNA enters nowhere.
* The landscape trains only where a solve locates a slot; the ruler reads its located enriched mode, and admits
  failure where there is no gDNA to read (owner, 2026-09-15).
* Real data is a test input, never a design input; the four cfRNA libraries are re-run with the regime printed.
* A few high-quality instruments, kept current; no suite gate polices instruments; the source cites no doc.
* The message cache's on/off switch, the refit count and the scan's thread split are the owner's. CI runs on
  demand only.
* The certifier's FIELD gate flake on λ ≈ 7 boundaries is DEFERRED (owner, 2026-09-14).
