# NEXT SESSION — start here (2026-09-17, after the port's step (iii))

This file is only how to begin. The port's plan is `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` §F (the order of its
four steps) and §G (the scan and the second pass after it); the LIVE status with numbers is
`ISSUES: performance-memory-bounded-solve`; steps (i)–(iii)'s record is the closing paragraphs of `DESIGN.md`
§6b.15.5; ψ's design is `docs/dev/PSI_PORT_PLAN.md`.

## The order (owner, 2026-09-17: the release ships the contraction as it stands; performance until the tool is fast)

1. ~~The ruler's repair~~ — DONE 2026-09-16.
2. ~~The yield's floors~~ — DONE 2026-09-17.
3. ~~The port, step (i): the directional pass native~~ — DONE 2026-09-17 (`553b3dfc`; VCaP 526 → 403 s).
4. ~~The port, step (ii): the builders~~ — DONE 2026-09-17 (`78d4a075`, `584fd9e6`; VCaP 409 → 315 s), and the
   one-path cleanup that followed it (`5c49447c`: the Python builders deleted).
5. ~~The port, step (iii): ψ native~~ — DONE 2026-09-17 (the commit this file rides in; ψ 4.3 → 1.6 s on
   MO_3021's first sweep, the replay 8e-6 of its budget; VCaP 341 → 285 s and 324 → 269 s at 8 threads, the
   sweep 0.73×, ψ 109 → 55 s, untouched stages 0.97–1.03 — `perf/port_iii_2026-09-17/`).
6. **THIS ONE — what is left of the sweep, and step (iv)**: read "Where the time is now" before choosing.
7. §G: the scan (32 s) and the second pass (22 s), the stages that scale with depth.

## Where the time is now (MO_3021's first sweep replayed, cProfile, 8.7 s)

| stage | s | what it is |
|---|---|---|
| `TransferPolicy.prepare` | 2.2 | the native builders ~0.7 s inside a Python wrapper of 1.45 s: `_Chain`'s `asarray` copies, the `LevelLane` constructions, and the TABLE ALLOCATIONS — `numpy.zeros` is 1.3 s of the whole sweep over 40,897 calls: every `RowTable` (own, three own levels, two `(n, 2, K)` flux tables) and both `Received` tables (composition + three level profiles, per pass) zero-filled per block, ~64 MB a block, ~27 GB a sweep |
| ψ (`_solve_regions_logodds_all`) | 1.65 | the native solve (1.6 s); the AMBIG cube's 2,600 exponentials per slot are the floor |
| `_PreparedTransfer.solve` | 1.34 | PYTHON: `_ceilings` 0.47 (per single-strand node, `rna_row_of_level`, `intersect`), `_cube_rows`, the held-level rows through `profile_of_level` (32 k calls) |
| `run_pass` (×2) | 0.85 | the native passes |
| the factory rows (`calibrate.__getitem__` → `density_lambda_factor`, `_log_negbinom`) | 0.6 | the intron factory's per-block λ rows, numpy |
| `solve_chain` / `_solve_block` self, `block_slice`, `view_fields`, the checks | ~1.0 | the block plumbing |

The native kernels are 3.2 s of 8.7; the Python around them is the majority. Two ways forward, to put to the owner:

* **(a) finish the block in native**: the policy's `solve` (the two held tables → ψ's `lam_rows` and `cube_rows`;
  the same shape as the builders — one call per block writing rows), and the allocations (`RowTable` and
  `Received` tables as `np.empty` where every reader goes through a mask — CHECK every reader first: the
  diagnostics capture concatenates the whole `Received` tables, so an instrument may read an absent row; or
  allocate once per sweep at the largest block and slice). Each a bounded step with the same gates.
* **(b) step (iv), threads**: ψ over slots inside `psi_solve` is deterministic (nothing is shared between
  slots) and one `thread_pool.h` loop away — but ψ is 19 % of the sweep now, so it caps at ~15 %; threads over
  BLOCKS parallelise everything, and need either the whole block in native (a) or the native calls releasing
  the GIL (`nb::call_guard<nb::gil_scoped_release>`) with Python threads driving `_solve_block` — the Python
  parts then still serialise. (a) first makes (b) worth more.

## The one-path ruling (owner, 2026-09-17) — what is still duplicated

One production code path; once native code is validated the Python it replaces is deleted; a small
floating-point tolerance is accepted. Applied to the builders (`5c49447c`) and to ψ (this commit). STILL
DUPLICATED, to converge next on the same pattern (bind the C++ pieces for the unit gates, rewrite the per-hop
gates to drive `run_pass` on small tables, delete the Python): the per-hop pass kernel
(`_PreparedTransfer.propagate`, `Faces.apply`, `LevelLane.emit` / `receive`) and the row constructors of
`transfer_rows.py` that only it and the gates' recomputes read (`face_map_lambda`, `edge_level_row`,
`level_map_lambda`, `level_bound_row`, `poisson_level`, `level_of_profile`, `rna_level_of_profile`, `flux_level`,
`transport_row`, `splice_out_row`, `level_row`, `blur_row`, the flag helpers, `strand_bits`);
`test_pass_kernel.py`'s two-kernel gate goes with them. What the policy's `solve` still reads in Python
(`intersect`, `lower_side`, `profile_of_level`, `rna_row_of_level`, `hop_price`, `count_logvar`) converges when
`solve` does. The layer-4 `strand_likelihood` "executable reference" module is the same kind of duplicate
(its gate now reads the native strand term through `psi_cube`) — the owner's call.

## The protocol for every step (unchanged)

1. `python scripts/design/preflight.py` first; the suite's standing count is in `CLAUDE.md`.
2. A fresh capture of the tree before any `src/` edit (`sweep_replay.py capture` on MO_3021 at 8 threads, ~5 min,
   3.5 GB; delete the superseded one).
3. DERIVE from the Python → the C++ in-tree (a `nanobind_add_module` in `CMakeLists.txt`; the import in
   `src/rigel/native.py` only AFTER the build) → a scratchpad harness that wraps the production call and runs
   both on every call of a captured sweep, comparing every output against the replay's derived budget
   (`s9/psi_ab.py`, `s8/prepare_ab.py`) → wire, delete the Python, rewrite the gates to the bindings → the
   suite → `sweep_replay.py replay --tolerance` for calls 0–3 → the three references re-frozen with the reason →
   two interleaved timing pairs at 8 threads on VCaP against a worktree of the pre-step commit carrying this
   build's binaries (`/tmp/rigel_pre`, `pre_site/sitecustomize.py`, `s9/time_pairs.sh`).
4. Each step its own commit; the owner drives the push.

## ⛔ Traps met in this session (also in memory)

* `np.finfo(np.float64).eps` is `std::numeric_limits<double>::epsilon()` (2⁻⁵²), NOT half of it: the quadrature's
  truncation constant read T = 36.7 instead of 36.0 and widened every AMBIG window — invisible on a toy,
  1e-4 on deep slots. Derive constants from the same definition on both sides and print them.
* A rotation recurrence for `sin` leaves [−1, 1] by an ulp at a domain end; a share `(1 ∓ τ)/2` below zero has no
  logarithm. Clamp to the function's range whenever a recurrence replaces the function.
* The harness must pass the PRODUCTION inputs: the toy agreed while real blocks (priors, message rows,
  delivered levels, a belief of `f_g = 1`) did not; catch the first mismatching call and dump the slot.
* A Python script that edits several files must write per file, and one that asserts several anchors in one
  file loses the whole file's edits on the first miss — re-run it, never patch by hand.
* Doc gates read `docs/` and `docs/dev/` while the suite runs: edit docs before launching it or after it ends.

## Storage (2026-09-17)

Freed today: the superseded test-reference sets (35 GB), the previous session's scratchpad (29 GB, synced first),
the `step6`/`step7`/`step8` captures (10.5 GB). Kept: `perf/sweeps_MO_3021_step9` (3.5 GB, this tree before
step (iii); the next step captures `step10` and deletes it). Still large and regenerable:
`prototypes/2026-09-16_ruler_repair/s5/scratch_test_reference*` (six scratch renders, 22 GB) — the owner's.

## Where everything is

* The synced scratchpad: `~/Downloads/rigel_runs/prototypes/2026-09-17_port_ii/` — `commits/committed/`
  (7–10 with their messages), `s8/` (step (ii)'s harnesses and logs), `s9/` (`psi_ab.py`, the step (iii) logs,
  `identity_prev/`, `time_pairs.sh`, `refreeze.sh`), `pre_site/`.
* Captures and reports: `perf/sweeps_MO_3021_step9`; `perf/port_ii_2026-09-17/`, `perf/port_iii_2026-09-17/`.
* The identity references: `~/Downloads/rigel_runs/arms/review_identity_*.json`, frozen on this tree.
* The worktree `/tmp/rigel_pre` (5c49447c) — `git worktree remove /tmp/rigel_pre` when the timing is done.

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
