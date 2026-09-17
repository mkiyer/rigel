# NEXT SESSION — start here (2026-09-17, after the port's step (i))

This file is only how to begin. The port's plan is `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` §F (the order of its
four steps) and §G (the scan and the second pass after it); the LIVE status with numbers is
`ISSUES: performance-memory-bounded-solve`; step (i)'s record is the closing paragraph of `DESIGN.md` §6b.15.5.

## The order (owner, 2026-09-17: the release ships the contraction as it stands; performance until the tool is fast)

1. ~~The ruler's repair~~ — DONE 2026-09-16 (four commits, `e51d91c6` … `64b79e9a`).
2. ~~The yield's floors~~ — DONE 2026-09-17 (`a79004b5`).
3. ~~The port, step (i): the directional pass native~~ — DONE 2026-09-17 (the commit this file rides in:
   `native.transfer_pass`; VCaP 526 → 403 s and 519 → 410 s at 8 threads).
4. **THIS ONE — step (ii): `prepare`'s builders write the tables the native pass reads**, then
5. step (iii): ψ native (`simplex_logodds`; the two ψ solves are 46 + 58 s of the 397 s sweep on VCaP), then
6. step (iv): threads over blocks (the blocks are independent given the library; `solve_chain`'s loop), then
7. §G: the scan (33 s) and the second pass (22 s), the stages that scale with depth.

## Step (ii), concretely

`TransferPolicy.prepare` (`src/rigel/calibration/messages/transfer.py`) builds, per block, every node's own claim
(`_claims`: a Python list of `(K,)` rows or `None`), the composition rules (`Faces`, already typed tables plus a
row store), and the level lanes (`lanes.gdna_lane`, `lanes.rna_lanes`: `LevelLane` objects whose `own_level` and
`flux_witness` are Python lists). `_PreparedTransfer._pack` turns those lists into `(n, K)` matrices with masks once
per block — 1.0 s a sweep on MO_3021 against the native kernel's 1.4 s — and the builders themselves are 101 s of
VCaP's 397 s sweep (`profiler.py`: RNA lanes 43, gDNA lane 22, claims 16, splice faces 11, alternative splice site 4,
terminus rules 4). Two sub-steps, each its own commit:

* **(ii-a) the layout** — the builders write matrices and masks directly (`own` as `(n, K)` + mask,
  `LevelLane.own_level` likewise, `flux_witness` as `(n, 2)` + mask, `Faces.rows` a matrix), `_pack` deleted. A pure
  re-layout: gate it BIT-IDENTICAL — the two-kernel harness pattern (`s7/pass_ab.py`, old prepare against new on the
  captured blocks, every table equal), the replay on a FRESH capture of this tree (`sweep_replay.py capture … --out
  perf/sweeps_MO_3021_step8`, five minutes; `step7` was captured before step (i) and already differs at ulp), the
  suite, and the three `review_identity_*` references, which must read BIT-IDENTICAL.
* **(ii-b) the arithmetic** — the builders' per-node work in C++ beside `pass_kernel.cpp` (one call per block,
  `transfer_prepare`): the strand profiles (`simplex_logodds.strand_row_logodds`), `level_of_profile`,
  `poisson_level`, `flux_level`, the face rules' rows (`face_map_lambda`, `level_map_lambda`, `edge_level_row`,
  `level_bound_row`), the face and lane bits. A port that changes summation order moves at ulp: gate it to the replay's
  `--tolerance` budget, re-freeze the references with the reason, time it on interleaved pairs.

Read the Python first; it is the executable specification, and every function above is unit-gated in
`tests/calibration/test_transfer_*.py` and `test_sweep_backbone.py`.

## The protocol for every step (unchanged since the plan)

1. `python scripts/design/preflight.py` first; the suite's standing count is in `CLAUDE.md`.
2. DERIVE from the Python → the C++ in-tree (`src/rigel/native/`, a `nanobind_add_module` in `CMakeLists.txt`, the
   import in `src/rigel/native.py`; `pip install --no-build-isolation -e ".[dev]"` rebuilds in 30 s) → a two-kernel
   harness on real captured blocks before wiring (copy `s7/pass_ab.py`: it monkeypatches `sweep._pass`, runs both
   kernels on every block and compares every field) → wire → `tests/calibration/test_pass_kernel.py`'s pattern for
   the gate file (equality on every bit and count, 1e-12 on profiles, a spy for the wiring) → `python
   scripts/profiling/sweep_replay.py replay --dir ~/Downloads/rigel_runs/perf/sweeps_MO_3021_step8 --call C
   --tolerance` for C in 0–3 → the suite → `rename_identity.py --check` on the three references → interleaved
   timing pairs at 8 threads on VCaP.
3. The timing arms: a worktree of the PRE-step commit with this build's binaries — `git worktree add
   /tmp/rigel_pre <SHA>`, copy `build/cp312-abi3-macosx_26_0_arm64/_*_impl*.so` into `<worktree>/src/rigel/`
   (copy to a new name, then `mv`), run it with `PYTHONPATH=<dir holding sitecustomize.py>` (the file that strips the
   editable install's redirecting finder is `shipped_site/sitecustomize.py` in the scratchpad copy) — then, in one
   sitting, `profiler.py --bam <VCaP> --index <index> --threads 8 --out pairN_pre.json` from the worktree and the
   same from the working tree, twice; `profiler.py --compare A.json B.json`; the untouched stages must read ~1.00.
   The libraries: VCaP `~/Downloads/rigel_runs/cfrna/mctp_vcap_rna20m_dna05m/bam/star.srt.rmdup.collate.bam`
   (18.6 M fragments, the timing library), MO_3021 (the capture library); the index `~/Downloads/rigel_runs/refs/rigel_index`.
4. Each step its own commit, PREPARED (`commits/N_name/files/` by `snapshot.sh`, the message beside it) and committed
   on the owner's go: move the committed snapshots to `commits/committed/`, `git stash push -u`, run
   `commit_series.sh`, verify `git diff stash@{0}` shows only the stash's untracked files and `cmp` those against
   HEAD, drop the stash, push.

## ⛔ Traps met in this session (also in memory)

* Never edit `src/` while a background measurement is starting runs: an import added before its module was built
  killed the second baseline pair and the capture at import.
* `import rigel.calibration.calibrate as X` binds the FUNCTION the package re-exports, not the module; patch
  submodule globals through `importlib.import_module`.
* A scripted edit that asserts several anchors must write per file, or a failed later anchor silently loses the
  earlier edits — `git status` and `grep` after every scripted edit.
* The jargon gate reads C++: a one-letter-one-digit type alias reads as a numbered rule label; name aliases.
* Counting: a `native/*.cpp` is +2 collected (jargon, docs-boundary); a `tests/` file +2 plus its tests.

## Where everything is

* The synced scratchpad: `~/Downloads/rigel_runs/prototypes/2026-09-16_ruler_repair/` — `commits/` (the snapshots,
  `snapshot.sh`, `commit_series.sh`), `s7/` (`pass_ab.py`, `pass_census.py`), `s6/` (the lever census and the
  yield draws), `s5/` (the ruler's arms), `real/` (the real-library harnesses), `shipped_site/sitecustomize.py`.
* Captures and baselines: `~/Downloads/rigel_runs/perf/sweeps_MO_3021_step7` (before step (i)),
  `perf/baseline_2026-09-17/pair1_*` (the pre-port tree, 524 / 528 s), `perf/port_2026-09-17/` (the two
  interleaved pairs of step (i), with the stage tree each side).
* The identity references: `~/Downloads/rigel_runs/arms/review_identity_*.json`, frozen on this tree.

---

## ⛔ Read first — three things the repair uncovered

* **The capture-OFF identity digest moves by one ulp** after commit 3 (one transcript's posterior mean, 1.1e-16 on
  the ladder's `g05 ss.50 OFF` row): the old floor's `w·span + (1−w)·span` arithmetic against the exact span. Every
  count and every calibration array is bit-identical; the references are re-frozen. Do not chase it.
* **A per-transcript Python loop over the real index is quadratic.** `transcript_piece_lengths` ran LBX0588 and
  VCaP for 70 minutes until `BaseTaper.interval_sums` took one template length per interval; now 160 s and 667 s
  against the shipped tree's 136 s and 616 s. Time every new per-transcript pass on `rigel_index` (457,371
  transcripts), never only on the ladder. Its perturbation found a hole in a green gate — no template length between
  one and two fragment lengths had been parametrised — closed with 600 and 998.
* **The ladder's panel spans junctions** (`ISSUES: ruler-witness-geometry-on-transcript-panels`, amended): the
  probed class reads 35 % within ±0.1 nat under every gDNA ruler, the certified true counts included, where the
  test chromosome's benign panel reads 99 %. A correction from the probe design is observable in principle; it is
  not this release's.

## The test chromosome moved (commit 3)

The tiny-exon block (`TESTING.md` §0a: ten 40 bp exons at 1,040 bp pitch, and the same run between two 1 kb exons,
unprobed and probed; 273 genes on 7.930 Mb, the budget unchanged) moved the test chromosome's benchmark rows: the
40 bp pieces are a new stress for the message layer (stranded ON 41,856 → 54,362, the ss 0.70 ON rows 64,174 →
105,997), the same on the shipped and the landed tree, so it is the substrate and not the mechanism. `DESIGN.md` §7
carries the new standing numbers; the superseded scenario sets are under
`~/Downloads/rigel_runs/test_reference_superseded_2026-09-16/`. The benign panel's rule gained one clause (a union
piece shorter than one probe gets a single probe centred on it).

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

## The session scratchpad (persists; nothing in the tree cites it)

Copy `/private/tmp/claude-503/-Users-mkiyer-proj-rigel/d290397d-5368-4194-86e0-5e1ccf118452/scratchpad/` to
`~/Downloads/rigel_runs/prototypes/2026-09-16_ruler_repair/` before it is lost: `commits/` (the four snapshots with
`MESSAGE.txt`, `snapshot.sh`, `commit_series.sh`), `s5/` (the arms `member_arms.py`, `ruler_arms.py`, `converge.py`,
the falsification harnesses, `gates3/` with every instrument's output on the landed tree, `identity_prev/` with the
identity references before commits 1 and 3), `real/` (`real_ruler.py`, `multimapper_check.py`, the four
`landed_<lib>/` outputs with `ruler.npz` and the check's tables), `head_tree/` + `shipped_site/sitecustomize.py`
(how the shipped tree is run beside an editable install: the scikit-build redirecting finder must be stripped).
