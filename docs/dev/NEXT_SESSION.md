# NEXT SESSION — start here (2026-09-18 evening: the block in one native call is built, bit-identical; snapshots 18–20 await the go; §G — the stages outside calibration — is next)

This file is only how to begin. The port's plan is `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` §F and §G; the LIVE
status with numbers is `ISSUES: performance-memory-bounded-solve`; the day's record is the closing paragraphs of
`DESIGN.md` §6b.15.5; THE BLOCK'S DESIGN, written before it was built, is `docs/dev/BLOCK_NATIVE_PLAN.md` (the
threads design it grew from is `docs/dev/THREADS_PLAN.md`).

## The order (owner, 2026-09-17/18: the release ships the contraction as it stands; performance until the tool is fast;
## the block in ONE native call, designed on paper first, preferred over Python threads over blocks)

1. ~~The ruler's repair, the yield's floors~~ — DONE 2026-09-16/17.
2. ~~The port, steps (i)–(iii): the pass, the builders, ψ; the block in native steps 1–2; the one-path convergence~~ —
   DONE 2026-09-17 (landed through `a1e5be03`).
3. ~~The optimisation target moved to VCaP~~ — DONE 2026-09-17 (`9e5131aa`).
4. ~~ψ threaded over slots; the cache key on the factory's inputs; the splice-out marginal's hoist~~ — DONE 2026-09-18,
   LANDED AND PUSHED as `81d53c8b`, `5d76d55b`, `4ace210e` (main == origin/main).
5. ~~The block in one native call~~ — BUILT 2026-09-18 as THREE SNAPSHOTS awaiting the go (`commits/18_psi_priors/`,
   `19_lgamma/`, `20_block_native/`; `commit_series.sh` replays them onto `4ace210e`):
   * 18 — ψ takes its priors apart: the gDNA arm interpolated per cell inside the kernel from the landscape's curve, the
     factory row and the delivered row two inputs; `landscape.logprior` deleted. EXACT.
   * 19 — the factory's log-gamma is libm's (for scipy's cephes `gammaln`): THE ONE NUMBER-MOVING COMMIT, `f_g` by at most
     5.7e-14 (3.7e-6 of the budget); the three identity references RE-FROZEN (the previous set in
     `arms/pre_lgamma_2026-09-18/`); the capture re-taken as `perf/sweeps_VCaP_step19`.
   * 20 — `native.solve_blocks`: ONE call per sweep, a pool of threads pulling the locus blocks, each block end to end
     (the prior rows, the self-solve ψ, the layer unless silent or served, the final ψ, the write-back, the counts) on
     its own arena; the Python per block deleted (`messages/faces.py`, `messages/lanes.py`, the backbone's per-block
     functions, `region_init`'s solve, `density_deconv`'s factor functions); one module `_solve_impl`
     (`solve_kernel.cpp` + `transfer_kernel.h` + `psi_kernel.h` + `transfer_rows.h` + `thread_pool.h`). BIT-IDENTICAL
     on the four VCaP sweeps at 8 threads and at 1, at another block size, on the served-cache path, on the three
     references; the suite 3,430 passed / 5 xfail / 3,435 collected; VCaP at 8 threads 198 → 144 s and 194 → 138 s (0.73 / 0.71), calibrate 85.6 → 33.7 s, the four sweeps 68.8 → 16.8 s and 66.3 → 15.9 s (0.24).
6. §G: the stages outside calibration — on VCaP the scan, the second pass, the two fragment-length fits, quant (the locus
   EM, the capture effective lengths), the index load; they scale with depth and are now the larger half of a run.
   Where the time is now (VCaP at 8 threads, commit 20's post runs, `perf/block_native_2026-09-18/`): the whole run 144 / 138 s. Calibrate 33.7 / 31.9 s: the four sweeps 16.8 / 15.9 (the kernel 15.7 / 14.9 — the first sweep ≈ 2.8 s, the refit that misses the cache ≈ 7.4 s, the two served refits ≈ 3 s each), the landscape fits 4.5, ψ before the sweep 0.65, and about 11 s of calibrate's OWN Python outside every probe (the chain, the statics, the beliefs' init and reset, the deconv) — the next thing to dissect inside calibration. Outside calibration 104 s: the scan 36 / 32, the second pass 22 (the fl models 8.4, scoring 12), the second fl fit 8.4, quant 36 (the locus EM 14.5, the capture effective lengths 5.7–8.2, scoring 6–7), the index load 6.6.

## The one-path ruling (owner, 2026-09-17) — applied to the block

One production code path; once native code is validated the Python it replaces is deleted; a small floating-point
tolerance is accepted. The whole block is that one path now: the gates read the ONE implementation through bindings
that allocate and return fresh tables (`native.transfer_prepare` / `transfer_pass` / `transfer_solve`,
`native.transfer_rows.*`), and the tables' containers (`RowTable`, `Faces`, `LevelLane`, `Received`, `Levels`) live in
`tests/calibration/_transfer_harness.py`, not in `src/`. What remains duplicated: the layer-4 `strand_likelihood`
"executable reference" module (its gate reads the native strand term through `psi_cube`) — the owner's call.

## Open for the owner (asked in the handoff report)

* THE MESSAGE CACHE — keep or delete. in production the last two refit sweeps are SERVED (the pre run's policy prepare ran 852 = 2 × 426 times over four sweeps, the first sweep and the first refit missing); with the block native a refit sweep that misses replays at 7.4 s at 8 threads and a served one at 3.0 s, so the cache saves about 9 s of a 140 s run (6 %) and costs 2.3 GB held through calibrate at K = 233 (the deliveries: the written rows with their slots, the cubes, the held bits), the keys (0.27 s a sweep) and the deliveries' round trip through Python — message_cache.py, the served list, the deliveries return and their gates; keep or delete is the owner's call. The kernel returns and takes back deliveries only so the cache can
  exist; deleting it removes `message_cache.py`, the keys, `served`, the deliveries' round trip and their gates.
* `policy_prototype.py` — its mechanism (a Python class installed in the backbone) cannot exist any more; it now scores
  the two SHIPPED policies per gene type / per node class / slot by slot and refuses any other arm. Keep as that, or
  retire (`policy_benchmark.py --by-class` covers the per-class view; the per-type table and `dissect` are its own).
  CLAUDE.md's working rule "prototyped outside `src/`" now means, for a message mechanism, a C++ builder in a WORKTREE
  scored against the tree without it — the rule's wording is the owner's.
* `strand_likelihood` (above). The blur's loop interchange (−1.5 s a first sweep, a summation-order change 1.9e-7 of the
  budget; `s15/edit_blur.py`) — priced, not taken.

## The protocol for every step (unchanged in form)

1. `python scripts/design/preflight.py` first; the suite's standing count is in `CLAUDE.md`.
2. The capture: `perf/sweeps_VCaP_step19` (9.0 GB; VCaP at 8 threads from the log-gamma tree) replays BIT-IDENTICAL on
   the block tree — `sweep_replay.py replay --dir … --call N --threads 8 [--block-slots N] [--tolerance]`; take a fresh
   one only after a commit that moves numbers, and delete the superseded one once its successor replays. ⛔ The replay
   carries NO message cache: calls 1–3 replay as cache misses (the layer runs at K = 233), which is why a refit sweep
   replays at ~7.4 s while a production refit sweep served from the cache runs at ~3 s.
3. DERIVE → the C++ in-tree → the replay on calls 0–3 (bit or `--tolerance`) → the references (re-frozen only with the
   reason) → the suite with its count re-derived → PERTURB the kernel and watch the replay fire, restoring from a saved
   copy → two interleaved timing pairs at 8 threads on VCaP against a worktree of the pre-step commit carrying its OWN
   modules (`s17/pre_worktree.sh` builds the pre tree's modules from its own source when an argument shape changed;
   `s17/time_pairs.sh` runs the pairs; `pre_site/sitecustomize.py` on `PYTHONPATH` points the worktree's import at itself).
4. Each step its own commit, prepared as a snapshot (`commits/snapshot.sh N_NAME file…` — pass the files literally,
   never through an unquoted variable; `DELETED.txt` beside `FILES.txt` for removals); the owner drives the commit.

## ⛔ Traps met in this session (also in memory)

* A CONSTANT TYPED FROM MEMORY: the nine splice-out marginal nodes hand-typed in the header were off at 1e-8 and moved
  `f_g` by 1e-7 on the replay. Print the oracle's `repr` and paste it; check `np.array_equal` against the oracle.
* FMA CONTRACTION: clang fuses `a*b + c` inside one expression (`-ffp-contract=on`), one rounding fewer than numpy's
  separate multiply and add — the factory rows differed at 1e-14 until every product became a named temporary. Do NOT
  set `-ffp-contract=off` for the module: numpy's own `np.interp` IS fused, so `interp`'s `fp + slope*(x − xp)` matches
  numpy only with contraction on. numpy's `sum(axis=1)` is a PAIRWISE sum (8 accumulators, blocks of 128, halving) —
  reproduced as `pairwise_sum` so the factor precision is exact.
* BRANCH ORDER: a served block must be checked BEFORE the silent policy's branch, and a silent block stores an EMPTY
  delivery when deliveries are wanted — otherwise a silent sweep's served-cache gate has nothing to serve.
* A PERTURBATION MUST TARGET WHAT THE PATH RELIES ON: removing the `held` clear proved nothing (the transfer path writes
  every held bit; the silent path never sets one); removing the own-mask clear or the received-bits clear moved
  600,000+ slots. Skipping a block in the pool moved 2,028 slots. The small-chain suite gates do not see a stale
  arena between blocks — only the replay does.
* The captured pickled policy carries the fields of the constructor it was pickled with; the replay rebuilds it through
  the CURRENT constructor (instrument-side), never a shim in `src/`.
* An identifier containing `capture` in the backbone trips the backbone gate — the kernel's diagnostics local is `diag`.

## Storage (2026-09-18 evening; 90 GB free)

`perf/sweeps_VCaP_step19` (9.0 GB) is THIS tree's capture (step16 deleted). The worktree `/tmp/rigel_pre` stands at
`4ace210e` + snapshots 18 and 19 with `_transfer_impl` and `_psi_impl` built from its own source (the timing baseline for
commit 20): `git worktree remove --force /tmp/rigel_pre` when it is next moved. Still large and regenerable:
`prototypes/2026-09-16_ruler_repair/s5/scratch_test_reference*` (22 GB) — the owner's.

## Where everything is

* The synced scratchpad: `~/Downloads/rigel_runs/prototypes/2026-09-17_port_ii/` — `commits/committed/` (7–17, landed),
  `commits/18_psi_priors/`, `19_lgamma/`, `20_block_native/` (awaiting the go, in order; `commit_series.sh` replays them
  onto `4ace210e`), `s16/` (18 and 19: the edit scripts, the pre copies of the deleted `.cpp` kernels, the tolerance
  logs), `s17/` (20: the header and kernel drafts, `perturb*.sh` and their logs, `identity_checks.sh`,
  `pre_worktree.sh`, `time_pairs.sh`, `finish_docs_20.py`), `pre_site/`.
* Captures and reports: `perf/sweeps_VCaP_step19`; `perf/block_native_2026-09-18/` (commit 20's two interleaved pairs);
  `perf/psi_threads_2026-09-18/`, `perf/splice_out_2026-09-18/`, `perf/vcap_baseline_2026-09-17/` (the earlier steps).
* The identity references: `~/Downloads/rigel_runs/arms/review_identity_*.json`, RE-FROZEN 2026-09-18 on the log-gamma
  tree's numbers (the reason in `DESIGN.md` §6b.15.5), BIT-IDENTICAL on the block tree; the previous set in
  `arms/pre_lgamma_2026-09-18/`.

## Decisions on record (unchanged, carried)

* Float64 for the whole of ψ; ONE solver; ONE λ lattice at a dimensionless step with a stated guarantee; no θ
  lattice anywhere; the tilt's hypothesis space is {pure +, pure −, mixed}.
* The strand channel's liveness is a protocol decision on the spliced 2×2; gDNA enters nowhere.
* The landscape trains only where a solve locates a slot; the ruler reads its located enriched mode, and admits
  failure where there is no gDNA to read (owner, 2026-09-15).
* Real data is a test input, never a design input; the four cfRNA libraries are re-run with the regime printed.
* A few high-quality instruments, kept current; no suite gate polices instruments; the source cites no doc.
* One production path (owner, 2026-09-17); `CalibrationConfig.n_threads` fed by `--threads`; the cache key on inputs;
  the block in ONE native call (owner, 2026-09-18). The message cache's existence, the refit count and the scan's
  thread split are the owner's. CI runs on demand only.
* The certifier's FIELD gate flake on λ ≈ 7 boundaries is DEFERRED (owner, 2026-09-14).
