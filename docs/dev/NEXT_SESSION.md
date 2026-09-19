# NEXT SESSION — start here (2026-09-18 night: the block in one native call is LANDED AND PUSHED (`7dad00d7`); the owner's four decisions are executed as snapshots 21–23 awaiting the go; §G — the stages outside calibration — is next)

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
5. ~~The block in one native call~~ — LANDED AND PUSHED 2026-09-18 as `dbc3d660`, `3590837c`, `7dad00d7` (main ==
   origin/main); the three commits:
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
6. ~~The owner's four decisions (2026-09-18 night)~~ — EXECUTED as THREE SNAPSHOTS awaiting the go (`commits/21_cache_deleted/`,
   `22_prototype_retired/`, `23_strand_reference/`; `commit_series.sh` replays them onto `7dad00d7`):
   * 21 — THE MESSAGE CACHE DELETED: every sweep runs the whole layer in the kernel; `message_cache.py`, the digests, the
     kernel's `served`/`deliveries` and `ServedBlock`, six gates gone; the replay capture rewritten without the pickled
     cache (9.0 → 2.7 GB). BIT-IDENTICAL on the four sweeps (the refits replay 7.4 → 7.0 s); suite 3,421 / 5 / 3,426.
     PRICE, two interleaved pairs against a worktree at `7dad00d7` with its own module: the run 141 → 148 s and
     138 → 145 s (1.06 / 1.05), the four sweeps 16.1 → 24.8 s (1.54), calibrate's peak RSS −1.3 / −1.6 GB
     (`perf/cache_deleted_2026-09-18/`).
   * 22 — `policy_prototype.py` RETIRED; the working rule in CLAUDE.md says how a mechanism is prototyped now (C++ in a
     worktree, both trees scored with `policy_benchmark.py --by-class`); suite 3,417 / 5 / 3,422.
   * 23 — `strand_likelihood.py` CONVERGED: `strand_loglik` lives in `tests/calibration/_psi_reference.py` beside the other
     oracles; layer 4 is production only; suite 3,414 / 5 / 3,419. The three identity references on the final tree:
     BIT-IDENTICAL on all three (capture off, capture on, the LBX0190 library); preflight `--full` green.
7. §G: the stages outside calibration — on VCaP the scan, the second pass, the two fragment-length fits, quant (the locus
   EM, the capture effective lengths), the index load; they scale with depth and are now the larger half of a run.
   Where the time is now (VCaP at 8 threads, commit 20's post runs, `perf/block_native_2026-09-18/`): the whole run 144 / 138 s. Calibrate 33.7 / 31.9 s: the four sweeps 16.8 / 15.9 (the kernel 15.7 / 14.9 — the first sweep ≈ 2.8 s, the refit that misses the cache ≈ 7.4 s, the two served refits ≈ 3 s each), the landscape fits 4.5, ψ before the sweep 0.65, and about 11 s of calibrate's OWN Python outside every probe (the chain, the statics, the beliefs' init and reset, the deconv) — the next thing to dissect inside calibration. Outside calibration 104 s: the scan 36 / 32, the second pass 22 (the fl models 8.4, scoring 12), the second fl fit 8.4, quant 36 (the locus EM 14.5, the capture effective lengths 5.7–8.2, scoring 6–7), the index load 6.6.

## The one-path ruling (owner, 2026-09-17) — applied to the block

One production code path; once native code is validated the Python it replaces is deleted; a small floating-point
tolerance is accepted. The whole block is that one path now: the gates read the ONE implementation through bindings
that allocate and return fresh tables (`native.transfer_prepare` / `transfer_pass` / `transfer_solve`,
`native.transfer_rows.*`), and the tables' containers (`RowTable`, `Faces`, `LevelLane`, `Received`, `Levels`) live in
`tests/calibration/_transfer_harness.py`, not in `src/`. Nothing is duplicated any more: the strand likelihood's
two-component reference is a gates' oracle (`tests/calibration/_psi_reference.py`), and the kernel's contract is arrays
in, arrays out, integer counts (under a capture: the cube rows and the received tables too).

## Decided by the owner (2026-09-18 night) and executed above

* THE MESSAGE CACHE — DELETED. Its price when it went: it saved about 9 s of a 140 s run (the two served refit sweeps at 3.0 s against 7.4 s) and cost 2.3 GB held through
  calibrate at K = 233, a hand-wired key and the deliveries' round trip through Python (the record: `DESIGN.md` §6b.15.5).
* `policy_prototype.py` — RETIRED. CLAUDE.md's working rule now says: prototyped OUTSIDE THE MAIN TREE — anything
  inside the block in C++ in a worktree, the two trees scored against each other on the same conditions; what is still
  Python, in Python the same way.
* `strand_likelihood` — CONVERGED (the readable two-component form lives with the gates' oracles).

## Still open for the owner

* The blur's loop interchange (−1.5 s a first sweep, a summation-order change 1.9e-7 of the budget; `s15/edit_blur.py`) —
  priced, not taken.
* Calibrate's own Python outside the sweep (~11 s of a 148 s run) and §G's stages: the next performance targets.

## The protocol for every step (unchanged in form)

1. `python scripts/design/preflight.py` first; the suite's standing count is in `CLAUDE.md`.
2. The capture: `perf/sweeps_VCaP_step19` (2.7 GB; VCaP at 8 threads from the log-gamma tree, the refit sweeps' inputs
   rewritten without the pickled cache) replays BIT-IDENTICAL on the cache-free tree — `sweep_replay.py replay --dir …
   --call N --threads 8 [--block-slots N] [--tolerance]`; take a fresh one only after a commit that moves numbers, and
   delete the superseded one once its successor replays. Every call runs the whole layer (the refits at K = 233, ~7 s).
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
  the CURRENT constructor (instrument-side), never a shim in `src/`. A captured kwarg whose CLASS is about to be deleted
  (the pickled `MessageCache`) must be stripped from the pickles BEFORE the class goes (`s18/strip_cache.py`).
* A SEQUENCE OF EDIT SCRIPTS UNDER `set -e` WITH `&&` CHAINS DOES NOT STOP when a script fails mid-chain: two anchor
  failures left the tree mixed and two snapshots captured partial states. Recover by stashing, re-applying the last
  verified snapshot and re-running each script; an anchor must never contain a placeholder that an earlier step fills.
* An identifier containing `capture` in the backbone trips the backbone gate — the kernel's diagnostics local is `diag`.

## Storage (2026-09-18 evening; 90 GB free)

`perf/sweeps_VCaP_step19` (2.7 GB) is THIS tree's capture (step16 deleted). The worktree `/tmp/rigel_pre` stands at
`7dad00d7` with `_solve_impl` built from its own source (the timing baseline for commit 21): `git worktree remove
--force /tmp/rigel_pre` when it is next moved. Still large and regenerable:
`prototypes/2026-09-16_ruler_repair/s5/scratch_test_reference*` (22 GB) — the owner's.

## Where everything is

* The synced scratchpad: `~/Downloads/rigel_runs/prototypes/2026-09-17_port_ii/` — `commits/committed/` (7–20, landed),
  `commits/21_cache_deleted/`, `22_prototype_retired/`, `23_strand_reference/` (awaiting the go, in order;
  `commit_series.sh` replays them onto `7dad00d7`), `s16/` (18 and 19), `s17/` (20: the header and kernel drafts,
  `perturb*.sh`, `identity_checks.sh`, `pre_worktree.sh`, `time_pairs.sh`, `finish_docs_20.py`), `s18/` (21–23:
  `strip_cache.py`, `commit22.py`, `commit23.py`, `snapshot_n.sh`, `pre_worktree.sh`, `time_pairs.sh`, the logs),
  `pre_site/`.
* Captures and reports: `perf/sweeps_VCaP_step19`; `perf/block_native_2026-09-18/` (commit 20's two interleaved pairs),
  `perf/cache_deleted_2026-09-18/` (commit 21's);
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
* One production path (owner, 2026-09-17); `CalibrationConfig.n_threads` fed by `--threads`; the block in ONE native
  call; the message cache deleted, `policy_prototype.py` retired, the strand reference converged (owner, 2026-09-18).
  The refit count and the scan's thread split are the owner's. CI runs on demand only.
* The certifier's FIELD gate flake on λ ≈ 7 boundaries is DEFERRED (owner, 2026-09-14).
