# NEXT SESSION — start here (2026-09-17, after the port's step (ii))

This file is only how to begin. The port's plan is `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` §F (the order of its
four steps) and §G (the scan and the second pass after it); the LIVE status with numbers is
`ISSUES: performance-memory-bounded-solve`; steps (i) and (ii)'s record is the closing paragraphs of `DESIGN.md`
§6b.15.5.

## The order (owner, 2026-09-17: the release ships the contraction as it stands; performance until the tool is fast)

1. ~~The ruler's repair~~ — DONE 2026-09-16 (four commits, `e51d91c6` … `64b79e9a`).
2. ~~The yield's floors~~ — DONE 2026-09-17 (`a79004b5`).
3. ~~The port, step (i): the directional pass native~~ — DONE 2026-09-17 (`553b3dfc`; VCaP 526 → 403 s).
4. ~~The port, step (ii): the builders~~ — DONE 2026-09-17, two commits prepared for the owner's go (the
   scratchpad's `commits/7_layout` and `commits/8_prepare_kernel`, snapshots with their messages): (ii-a) the
   layout (`RowTable`, `_pack` deleted, bit-identical) and (ii-b) `native.transfer_prepare` (the builders
   20.0 → 2.5 s on MO_3021's first sweep; VCaP 409 → 315 s and 403 → 314 s at 8 threads, the sweep 0.68×, the
   builders 0.18×, ψ and every untouched stage 0.93–1.03 — `perf/port_ii_2026-09-17/`).
5. **THIS ONE — step (iii): ψ native** (`simplex_logodds`), then
6. step (iv): threads over blocks (the blocks are independent given the library; `solve_chain`'s loop), then
7. §G: the scan (33 s) and the second pass (22 s), the stages that scale with depth.

## Step (iii), concretely

After step (ii) the sweep's time on VCaP is ψ's: the stage tree of step (ii)'s port run (`perf/port_ii_2026-09-17/
pair1_port.log`) reads 45.4 s for the self-solve's ψ (`build_region_init`, 1,704 calls) and 57.3 s for the final
solve's (`_solve_block`, 1,704 calls) of a 191 s sweep, beside 18 s of builders (one native call per block plus the
Python allocation and wrapping), 36.7 s of native passes and 6.9 s of the policy's solve. ψ is `simplex_logodds._solve_regions_logodds_all` → `_solve_logodds` → `_psi`: one solver for
both slot classes on the `(m, K, K_t + 2)` cube in float64 (`DESIGN.md` §6b.15.6) — the strand term, the two
arms (`_gdna_arm`, `_rna_arm`), the λ-factor rows, the delivered `lam_rows` and, at AMBIG slots, the `CubeRow`
delivery evaluated at ψ's own θ nodes (`_tilt_window`, `_TILT_NODES` derived), then the read-out: `f_g` the
posterior median over the θ-marginal (`_posterior_median_fg`, a continuous quantile), `Var(log f_g)` its grid
moment (`_row_moment`), the tilt share, `_compose`. Everything is numpy over tiles of rows (`_block_rows`). The
port is the same shape as steps (i) and (ii): DERIVE the C++ from the Python (the executable specification),
build it in-tree beside `pass_kernel.cpp` / `prepare_kernel.cpp` (a `_psi_impl` module; the row pieces of
`native/transfer_rows.h` are not ψ's — ψ needs `exp`/`log` over the cube, the log-sum-exp `_lse`, the quantile),
a two-kernel harness on real captured blocks BEFORE wiring (copy the pattern of the scratchpad's
`s8/prepare_ab.py`: wrap the production call, run both, compare every output field), wire, gate
(`tests/calibration/test_psi_kernel.py` on the pattern of `test_prepare_kernel.py`: every output within the
budget, the wiring by a spy), the replay's `--tolerance` on `sweeps_MO_3021_step8` for calls 0–3, the suite, the
three references re-frozen with the reason, interleaved timing pairs on VCaP. ⛔ Float64 for the whole of ψ,
ONE solver, the λ lattice a parameter (`sweep_logodds_step` through `calibrate.lattice_points`) — decisions on
record. ⛔ Read `tests/calibration/test_vertex_reference.py`, `test_simplex_logodds*.py` and `test_sweep.py`
first: the read-out's rules (the median not the mean, the tilt atom at τ = ±1, the reference measure) are what
a port must keep to the bit-budget, not to bits.

## The protocol for every step (unchanged since the plan)

1. `python scripts/design/preflight.py` first; the suite's standing count is in `CLAUDE.md`.
2. DERIVE from the Python → the C++ in-tree (`src/rigel/native/`, a `nanobind_add_module` in `CMakeLists.txt`, the
   import in `src/rigel/native.py`; `pip install --no-build-isolation -e ".[dev]"` rebuilds in 30 s) → a two-kernel
   harness on real captured blocks before wiring → wire → the gate file (equality on every bit and count, a
   budget on every profile, a spy for the wiring) → `python scripts/profiling/sweep_replay.py replay --dir
   ~/Downloads/rigel_runs/perf/sweeps_MO_3021_step8 --call C --tolerance` for C in 0–3 → the suite →
   `rename_identity.py --check` on the three references (they MOVE under a port: re-freeze with the reason, the old
   ones kept under the scratchpad) → interleaved timing pairs at 8 threads on VCaP.
3. The timing arms: a worktree of the PRE-step commit with this build's binaries — `git worktree add
   /tmp/rigel_pre <SHA>`, copy `build/cp312-abi3-macosx_26_0_arm64/_*_impl.abi3.so` into `<worktree>/src/rigel/`
   (copy to a new name, then `mv`), run it with `PYTHONPATH=<dir holding sitecustomize.py>` (the file that strips the
   editable install's redirecting finder and puts the worktree's `src` first: `pre_site/sitecustomize.py` in the
   scratchpad; edit its `PRE` path) — then, in one sitting with nothing else running, the scratchpad's
   `s8/time_pairs.sh` pattern: pre, port, pre, port; `profiler.py --compare A.json B.json`; the untouched stages must
   read ~1.00. The libraries: VCaP `~/Downloads/rigel_runs/cfrna/mctp_vcap_rna20m_dna05m/bam/star.srt.rmdup.collate.bam`
   (18.6 M fragments, the timing library), MO_3021 (the capture library); the index `~/Downloads/rigel_runs/refs/rigel_index`.
4. Each step its own commit, PREPARED (`commits/N_name/files/` by `snapshot.sh`, the message beside it) and committed
   on the owner's go: move the committed snapshots to `commits/committed/`, `git stash push -u`, run
   `commit_series.sh`, verify `git diff stash@{0}` shows only the stash's untracked files and `cmp` those against
   HEAD, drop the stash, push.

## ⛔ Traps met in this session (also in memory)

* A pure re-layout is gated by DIGESTS of every table from two trees (`s8/prepare_digest.py`: the pass's kwargs
  dict per block from the pre-step worktree and from the tree), never by reading the code twice.
* Identity checks and the suite measure the tree they STARTED on: finish every `src/` edit of a commit first, then
  run them; a check started before the last edit is wasted (and editing `src/` while one starts is the step-(i)
  trap). Docs may be edited while they run; they read no docs.
* The jargon and docs-boundary gates collect `.h` files too: a header under `src/rigel/native/` is +2, like a `.cpp`.
* The compiler contracts `a*b + c` into a fused multiply-add: a C++ port agrees with numpy at a few ulps of the raw
  terms, never to the bit, even where `exp`, `log` and `expit` match libm exactly — so a gate on real data is a
  budget (the toy's small rows passed the existing gates' 1e-10..1e-12 tolerances unchanged).
* A scripted edit that asserts several anchors must write per file — one that asserts after writing loses nothing,
  one that asserts BEFORE writing loses the whole edit silently; `grep` after every scripted edit.

## An owner decision to take

The Python builders (`transfer._prepare_reference`, with `_claims`, `_splice_faces`, `_edge_level`,
`_terminus_rules`, `_alternative_splice_site`, `lanes.gdna_lane`, `lanes.rna_lanes`) no longer run in production:
they are the executable specification `test_prepare_kernel.py` holds the native call to, on the layer-4
`strand_likelihood` precedent. The alternative is to delete them and hold the C++ to the unit gates of
`test_transfer_*.py` alone (independent recomputes from `transfer_rows`, which those gates already pass on the
native builders). The same question will arise for ψ.

## Where everything is

* The synced scratchpad: `~/Downloads/rigel_runs/prototypes/2026-09-17_port_ii/` — `commits/` (`7_layout`,
  `8_prepare_kernel`, `snapshot.sh`, `commit_series.sh`), `s8/` (`prepare_digest.py`, `prepare_ab.py`,
  `time_pairs.sh`, every gate's log, `identity_prev/` with the references before the re-freeze), `pre_site/`.
* Captures and baselines: `~/Downloads/rigel_runs/perf/sweeps_MO_3021_step8` (this tree before step (ii), 8a9c25c4),
  `perf/port_2026-09-17/` (step (i)'s pairs), `perf/port_ii_2026-09-17/` (step (ii)'s pairs).
* The identity references: `~/Downloads/rigel_runs/arms/review_identity_*.json`, frozen on the (ii-b) tree.
* The pre-step worktree `/tmp/rigel_pre` (8a9c25c4) — remove it with `git worktree remove /tmp/rigel_pre` when the
  timing is no longer needed; `git worktree prune` after.

## Decisions on record (unchanged, carried)

* Float64 for the whole of ψ; ONE solver; ONE λ lattice at a dimensionless step with a stated guarantee; no θ
  lattice anywhere; the tilt's hypothesis space is {pure +, pure −, mixed}.
* The strand channel's liveness is a protocol decision on the spliced 2×2; gDNA enters nowhere.
* The landscape trains only where a solve locates a slot; the ruler reads its located enriched mode, and admits
  failure where there is no gDNA to read (owner, 2026-09-15).
* Real data is a test input, never a design input; the four cfRNA libraries are re-run with the regime printed.
* A few high-quality instruments, kept current; no suite gate polices instruments; the source cites no doc.
* The message cache's on/off switch, the refit count and the scan's thread split are the owner's; parallelism
  waits for the port. CI runs on demand only.
* The certifier's FIELD gate flake on λ ≈ 7 boundaries is DEFERRED (owner, 2026-09-14).
