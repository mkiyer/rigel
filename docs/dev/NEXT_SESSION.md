# NEXT SESSION — start here (2026-09-17, after the tables' allocations; the threads design awaits the owner)

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
5. **THIS ONE — step (iv), threads: read `docs/dev/THREADS_PLAN.md`, take the owner's answers to its three open
   questions, then build (iv-a) — ψ over slots inside `psi_solve` — as its own commit**: the EM's pool
   (`native/thread_pool.h`), a cost-balanced partition of the slot list (26 columns an AMBIG slot, 1 a single-strand
   one), a `Scratch` per thread, the GIL released; gate: the replay of `sweeps_MO_3021_step11` at 1, 2 and 8 threads
   BIT-IDENTICAL on all four sweeps, the suite, the three references; two interleaved pairs on VCaP at `--threads 8`.
   Then the cache key on the factory's inputs, the gDNA arm inside ψ, and the block in one native call — in that order,
   each its own commit, each gated as the plan says.
6. ~~The one-path convergence of the pass~~ — DONE 2026-09-17, prepared as snapshot `commits/13_one_path/`
   (the per-hop Python kernel and `transfer_rows.py` deleted; the constructors bound for the gates).
7. §G: the scan (35 s) and the second pass (22 s), the stages that scale with depth.

## Where the time is now (MO_3021 replayed, this tree; `s10/kernel_times_quiet.log`)

| per sweep | first sweep (K = 101) | refit sweep (K = 202, cache-served) |
|---|---|---|
| ψ native (self-solve + final) | 1.57 s | 3.19 s |
| the builders native | 1.61 s | — |
| the pass native (both) | 0.99 s | — |
| the policy's solve native | 0.12 s | — |
| the message cache's key (blake2b over the block context, the factory rows most of it) | — | 2.88 s |
| the gDNA arm (`landscape.logprior`, `np.interp` over (m, K)) | — | 1.51 s |
| the factory rows (`_log_negbinom`; lgamma 85 %) | 0.46 s | 0.90 s |
| the glue (the loop, `_psi`'s row add, `build_region_init`, `block_slice`, the checks, the write-back) | ~1.9 s | ~1.3 s |
| the sweep | 6.7 s (native 64 %) | 10.0 s (native 32 %) |

⛔ The refit sweeps run on the landscape's bracket — K = 202 at L = 20.09 against the first sweep's K = 101 at
L = 10 — so every per-cell cost is 2× there for the same slots; profile a REFIT sweep (`replay --call 1`), not only
the first. In production sweep 1 also misses the cache and runs the whole layer at K = 202. VCaP at 8 threads: the
sweep is 135 s of 264 s, ψ 53 s of it; the allocation step 277 → 275 s and 271 → 268 s at 8 threads (0.99 / 0.99), the builders' stage 19.2 → 17.9 s and 18.7 → 17.3 s, the pass 38.4 → 38.0 s and 37.7 → 36.8 s.

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
2. The capture: `perf/sweeps_MO_3021_step11` (3.6 GB, taken from `28bd174c`) replays BIT-IDENTICAL on this tree — use
   it; take a fresh one (`sweep_replay.py capture` on MO_3021 at 8 threads, ~5 min) only after a commit that moves
   numbers, and delete the superseded one once its successor replays.
3. DERIVE from the Python → the C++ in-tree → a scratchpad harness that runs both on every call of a captured sweep
   against the replay's derived budget → wire, delete the Python, rewrite the gates → the suite →
   `sweep_replay.py replay [--tolerance]` for calls 0–3 → the references (re-frozen only with the reason) → two
   interleaved timing pairs at 8 threads on VCaP against a worktree of the pre-step commit carrying the SAME
   binaries (`/tmp/rigel_pre` + `pre_site/sitecustomize.py`, `s10/time_pairs.sh`).
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

## Storage (2026-09-17)

`perf/sweeps_MO_3021_step11` (3.6 GB) is THIS tree's capture and stays until a step moves numbers. The worktree
`/tmp/rigel_pre` stands at `28bd174c` with this tree's binaries — `git worktree remove --force /tmp/rigel_pre`
when it is next moved. Still large and regenerable: `prototypes/2026-09-16_ruler_repair/s5/scratch_test_reference*`
(22 GB) — the owner's.

## Where everything is

* The synced scratchpad: `~/Downloads/rigel_runs/prototypes/2026-09-17_port_ii/` — `commits/committed/` (7–11),
  `commits/12_alloc/` (this step's snapshot and message), `s8/`, `s9/`, `s10/` (this session: the edit scripts, the
  perturbation logs, the kernel timings, the ψ ablations, the replay and identity logs), `pre_site_28bd174c/`.
* Captures and reports: `perf/sweeps_MO_3021_step11`; `perf/block_alloc_2026-09-17/` (this step's four profiler
  reports), the earlier `perf/port_*_2026-09-17/`.
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
