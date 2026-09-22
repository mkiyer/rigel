# CLEANUP — the ledger (started 2026-09-22)

The owner's instruction: restore the tree to a clean, concise, production-ready state WITH THE SHIPPED
INFRASTRUCTURE — no new mechanism, no new feature — as if 0.8.0 shipped what runs today. Every step is a
deletion or a simplification proven a numeric no-op (`rename_identity.py --check` bit-identical on its three
references; the suite green at a re-derived count), or an owner decision listed here until it is made. This file
is the working ledger; when the cleanup is done it is deleted and `NEXT_SESSION.md` points at the two problems.

## Principles

* One kind of change per commit: retire a file; delete a knob; converge a duplicate. The owner drives commits.
* A deletion is proven by the suite AND the identity check (`scripts/design/rename_identity.py --check`, the
  three references frozen at `71037f6c`), never by the suite alone (`TRAPS: a-green-suite-hid-five-dead-instruments`).
* What a knob does at its default is what ships; a knob nobody has derived and no instrument arm uses is deleted
  with its CLI flag, its manual row and its tests.
* Docs follow the code: a permanent doc is not a history (`CLAUDE.md`'s own rule); numbers live where an
  instrument re-derives them.

## Done

| date | what | proof |
|---|---|---|
| 09-22 | The per-fragment-prior derivation line STOPPED by the owner: `docs/dev/PER_FRAGMENT_GDNA_PRIOR{,_REVIEW}.md`, the four `CALIBRATED_LIKELIHOOD_*.md` notes, `NASCENT_SIPHON.md`, `CAPTURE_WITHOUT_THE_PANEL.md` retired; `scripts/design/gdna_prior_reliability.py` retired with its row and pin; `TWO_PROBLEMS.md` restated in the owner's words | docs and an instrument: no number moves |
| 09-22 | `CLAUDE.md`'s standing-baseline paragraph — a 60-line history — condensed to the count and the rule | text |
| 09-22 | `pipeline._DEFAULT_MEAN_FRAG = 200.0`, the magic fallback for an empty RNA length pool (never executed by the suite), replaced by a refusal that names the condition; the falsification test in `tests/test_d7_transcript_eff_lengths.py` was written first and watched fail | suite; no number moves (the branch was unreachable) |
| 09-22 | `EMConfig.gdna_em_llr_bias` and `--gdna-em-llr-bias` deleted (config, CLI, estimator, `em_solver.cpp`, manual, two tests): an underived odds knob, default 0, used by no arm | suite green; `rename_identity.py --check` BIT-IDENTICAL on both ladder references (log `~/Downloads/rigel_runs/arms/`) |

The suite after the above: 0 failed / 3,454 passed / 2 xfail, 3,456 collected — derived from the table in `CLAUDE.md`
(from 3,468: −4 the instrument, −1 its pin, −8 the retired `docs/dev/` notes, −2 the knob's two tests, +1 the refusal's test,
+1 this ledger, +1 the owner's `docs/dev/my_notes.md`) and matched by `--collect-only`.

## The census (2026-09-22)

* `src/rigel`: 46,303 lines (Python + C++), `scripts/`: 13,239, `tests/`: 133 files / 76,040 lines — the tests
  are 1.6× the source. `module_census.py`: no upward import, no dead public surface in `calibration/`.
* Coverage of `src/rigel` by the suite (2026-09-22, `pytest --cov=rigel`): 89.4 % of 10,839 Python statements
  (C++ is coverage-blind). Whole functions the suite never executes, by file — each is read and either kept as a
  named coverage gap (a shipped path the suite does not drive) or deleted:
  - `sim/whole_genome.py` (284 missing): the `rigel sim` driver — `main`, `build_arg_parser`, `run_simulation`,
    `ensure_index`, the abundance loaders (`assign_random_abundances`, `assign_file_abundances`,
    `_detect_abundance_format`, `_load_abundance_map`), `_build_nrna_pairs`, `nrna_label_for_ratio`. A coverage
    GAP by the 2026-09-14 census (the panel is built through `scripts/sim/panel.py`, which drives the same engine).
  - `sim/wgs_engine.py` (107): the sharded writers `_write_all_parallel`, `_run_shard`, `_shard_by_count`, and
    `_introduce_errors_batch` (sequencing errors: the panel simulates none). GAP or dead — the owner's call: if
    the panel never simulates errors, the error model is speculative.
  - `sim/orchestrator.py` (68): `run_condition_grid`, the whole-genome grid. GAP (2026-09-14).
  - `splice_blacklist.py` (29): `load_splice_blacklist_from_zarr` — the zarr blacklist at index build. GAP.
  - `buffer.py`: `_weak_cleanup` (the finalizer) and the spill writer's error path — defensive, keep.
  - `pipeline.py:511`: the `_DEFAULT_MEAN_FRAG = 200.0` fallback is NEVER executed — production reaches
    `_setup_geometry_and_estimator` with a pmf-built model. Replaced by a refusal (below).
  - `frag_length_model._normalized_probs`: the `total == 0` uniform fallback, the same smell, unexecuted.
  - Partly executed, all error/validation branches: `splice_graph.validate_graph`, `index.load`'s refusals,
    `scan_payload.from_dict`, `strand_summary.__post_init__`, `report/model._verdicts`.
  - `config.py`, `estimator.py`, `cli.py` were edited during the run and are re-measured next session.
* Rot markers (`TODO`, `legacy`, `deprecated`, `compat`, `unused`, `SRD v2`): 76 lines, most of them the word
  "compatible" in its biological sense. Actionable: the `SRD v2` spec label in 14 comments (`constants.h`,
  `resolve_context.h`, `scoring.cpp`, `bam_scanner.cpp`, `splice.py`, `buffer.py`) — a revision name nobody
  can look up, to be rewritten as what the code does; `estimator.run_batch_locus_em_partitioned`'s two
  "wrapper-level test compatibility" defaults (`rna_prior_count`, `gdna_eff_len` optional in production for the
  tests' sake); `index.py`'s format-version history comment; `bam_scanner.cpp:1673` "per the plan's TODO".
* Branches: ten stale local branches beside `main` and `calibrated-likelihood`; one worktree
  (`/Users/mkiyer/proj/rigel-pin`, the prototype switches `RIGEL_GDNA_OPP` / `RIGEL_PIN_GDNA` / `RIGEL_EFF_FORM`,
  its own build). Prototype harnesses and results under `~/Downloads/rigel_runs/prototypes/`.

## Decisions for the owner (each blocks a deletion)

1. **The per-transcript prior lane** `rna_prior_weight` (pipeline → estimator → `em_solver.cpp`) and the
   `EMConfig.warm_start` modes `prior` / `uniform`: plumbed end to end, no producer in `src/`, used only by
   `quant_accuracy.py`'s `oracle_alloc` / `warm_uniform` arms and five tests. Under "no speculative code" both go;
   under "problem 1 is a prior per transcript" the lane is the mechanism that work would fill. Delete now and
   re-add designed, or keep? (Deleting also removes the allocation ceiling arm.)
2. **`EMConfig.mode = "map"`**: VBEM ships; MAP is an A/B partner with gates (`EQUATIONS.md` §9b.2). Keep as a
   measured alternative or delete?
3. **`CalibrationConfig.background_abundance = "measured_total"`**: an unruled alternative pair
   (`ISSUES: background-abundance-pair-unruled`). Rule it or delete it.
4. **`CalibrationResult.rna_region_eff_len` / `rna_boundary_eff_len`**: no consumer in `src/`; the result's own
   docstring argues they stay. Keep or delete?
5. **The instrument shelf** (19 files after today): the bloat ruling says few. Proposed retirements, each with the
   reason: `pass0_vs_oracle.py` (its headline is the wrong yardstick by its own row; its `C_input` arms hand the
   post-capture law to geometry — the ledger's "rotten but live"), `worst_objects.py` (built on it),
   `toy_harness.py` + `zero_controls.py`'s toy path (the harvest calibrates undrained), `module_census.py`,
   `rename_census.py` (one-off vocabulary tools), `calibration_walk.py` (a ladder no session has read since the
   sweep went native), `simulator_gates.py` (run once per panel build). Kept: `preflight`, `policy_benchmark`,
   `calibration_vs_oracle`, `calibration_oracle`, `ruler_vs_truth`, `quant_accuracy`, `solvability_audit`,
   `prior_vs_oracle`, `zero_controls` (ladder path), `rename_identity`, `build_scan_cache`, the four `sim/` rows.
6. **The worktree and the branches**: remove `rigel-pin` (its prototype is refused or superseded) and the ten
   stale branches? Reversible only through the reflog.
7. **`pipeline._DEFAULT_MEAN_FRAG = 200.0`**: the magic fallback for an empty RNA length pool. Replace with a
   refusal (a library with no RNA length observation cannot be quantified) — if coverage shows the branch is
   never executed by the suite, it is a refusal with no test to move.

## Stages ahead

1. **`src/` knobs and lanes** — the decisions above, then each deletion alone with the identity check.
2. **Dead and duplicated code** — from the coverage table: every function the suite never executes is read and
   either gated or deleted; the ledger's duplicates (`ISSUES: hygiene-ledger` (c)) converged into `_shared.py`.
3. **Tests** — the vacuous gate clauses named in `ISSUES: hygiene-ledger` (b) deleted; tests of deleted code
   removed with it; the count re-derived from the table, never adjusted.
4. **Instruments** — decision 5, then `CLAUDE.md`'s table and `preflight.py`'s import sweep follow the disk.
5. **Docs** — `ROADMAP.md` carries numbers against its own rule; `CLAUDE.md` is 380 lines; `DESIGN.md` 2,015.
   Each cut by the move rule: a number goes to the instrument that re-derives it, a history to git.
6. **Release readiness** — `CHANGELOG.md`'s `[Unreleased]` written from git since 0.7.1; `pyproject.toml`
   classifiers; `pip wheel` from a clean checkout; the CLI smoke run on the test chromosome; `preflight.py --full`;
   `PUBLISHING.md`'s procedure walked once without publishing.
