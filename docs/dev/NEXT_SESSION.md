# NEXT SESSION — start here (2026-09-17, after ψ and the policy's solve went native)

This file is only how to begin. The port's plan is `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` §F and §G; the
LIVE status with numbers is `ISSUES: performance-memory-bounded-solve`; the day's record is the closing paragraphs
of `DESIGN.md` §6b.15.5; ψ's design is `docs/dev/PSI_PORT_PLAN.md`.

## The order (owner, 2026-09-17: the release ships the contraction as it stands; performance until the tool is fast;
## finish the block in native BEFORE threads)

1. ~~The ruler's repair, the yield's floors~~ — DONE 2026-09-16/17.
2. ~~The port, steps (i)–(iii): the pass, the builders, ψ~~ — DONE 2026-09-17 (`553b3dfc`, `78d4a075`, `584fd9e6`,
   `80b459bf`), with the one-path cleanup `5c49447c`.
3. ~~The block in native, step 1: the policy's solve~~ — DONE 2026-09-17 (the commit this file rides in: bit-identical;
   the cube delivery a table; the transfer's kernels one module).
4. **THIS ONE — the block in native, step 2: the tables' allocations, then the factory rows** (below).
5. Step (iv), threads: ψ over slots inside `psi_solve` (deterministic), the blocks in parallel once a block is
   native calls with little Python between them.
6. §G: the scan (32 s) and the second pass (22 s), the stages that scale with depth.

## Where the time is now (MO_3021's first sweep replayed, cProfile, 7.2 s)

| stage | s | what it is |
|---|---|---|
| `TransferPolicy.prepare` | 2.1 | the native builders ~0.7 s inside a Python wrapper of 1.4 s — mostly the TABLE ALLOCATIONS: `numpy.zeros` is 1.25 s of the whole sweep over 41,363 calls (every `RowTable` — own, three own levels, two `(n, 2, K)` flux tables — and both `Received` tables per pass, ~64 MB a block, ~27 GB a sweep) |
| ψ (`_solve_regions_logodds_all`) | 1.6 | native; the AMBIG cube's 2,600 exponentials per slot are the floor (52 k AMBIG against 306 k single-strand slot-solves) |
| `run_pass` (×2) | 0.8 | native |
| the factory rows (`calibrate.__getitem__` → `density_lambda_factor` → `_log_negbinom`) | 0.6 | the intron factory's per-block λ rows: a NegBinom log-likelihood over (introns, K), numpy with `gammaln` |
| `solve_chain` / `_solve_block` self, `block_slice`, `view_fields`, `asarray`, the checks | ~0.9 | the block plumbing |
| the policy's `solve` | 0.2 | native (0.13) plus its wrapper |

Step 2, concretely:

* **The allocations.** Every reader of a `RowTable`'s matrix and of a `Received` table's profiles goes through a
  mask (the kernels; `solve` natively; the cache stores `lam_rows` sparsely by non-zero rows — CHECK: a `Received`
  composition row that is absent must stay zero for the cache's sparsity to hold, so `Received.composition` keeps
  its zeros or the cache keys on `has_composition`); the diagnostics capture concatenates whole tables but no
  instrument reads a row unmasked (grep `scripts/` for `.composition`, `.profile`: none). So: `RowTable` on
  `np.empty` (its contract already says the mask, never the matrix, says whether a row is a claim), the level
  tables' profiles likewise, and `Faces.rows` is already `np.empty`. Alternatively one arena per sweep sized at
  the largest block, sliced per block. Gate: bit-identical (the replay on a fresh capture, the references, the
  suite); the reader audit written into the commit.
* **The factory rows.** `density_lambda_factor` per block is a `(n_intron, K)` NegBinom log-likelihood — small in
  C++ (`lgamma` per cell) but a layer-5 concern; port it only if it still shows after the allocations.
* Then step (iv).

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
* A capture pickles the refit sweeps' MESSAGE CACHE: retiring a class the cache held (`CubeRow`) makes those
  sweeps' captures unloadable on the new tree; sweep 0 (no cache entries) still replays. Re-capture after a step
  that changes a delivered type, and gate that step on sweep 0, the references and the suite.
* A gate that called one step of a pipeline alone (`_ceilings` on a zero row) expects that step's output; when
  the step is folded into the whole (`solve`), the gate must expect the whole — here the held composition fused
  with the ceiling.

## Storage (2026-09-17)

Freed today: the superseded test-reference sets (35 GB), the previous session's scratchpad (29 GB, synced first),
the `step6`–`step9` captures (14 GB). Kept: `perf/sweeps_MO_3021_step11` (3.5 GB, THIS tree; the next step
re-captures and deletes it). Still large and regenerable:
`prototypes/2026-09-16_ruler_repair/s5/scratch_test_reference*` (six scratch renders, 22 GB) — the owner's.

## Where everything is

* The synced scratchpad: `~/Downloads/rigel_runs/prototypes/2026-09-17_port_ii/` — `commits/committed/`
  (7–11 with their messages), `s8/` (step (ii)'s harnesses and logs), `s9/` (`psi_ab.py`, `solve_ab.py`, the
  step (iii) and solve logs, `identity_prev/`, the timing and re-freeze scripts), `pre_site/`.
* Captures and reports: `perf/sweeps_MO_3021_step11`; `perf/port_ii_2026-09-17/`, `perf/port_iii_2026-09-17/`,
  `perf/port_solve_2026-09-17/`.
* The identity references: `~/Downloads/rigel_runs/arms/review_identity_*.json`, frozen on this tree.
* The worktree `/tmp/rigel_pre` (80b459bf, the ψ commit) — `git worktree remove /tmp/rigel_pre` when the timing is
  done; the two retired modules' binaries it needs are under `s9/pre_so/`.

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
