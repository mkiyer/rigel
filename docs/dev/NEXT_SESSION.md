# NEXT SESSION — the state after the cleanup session (2026-09-10/11)

A handoff, provisional like everything in this directory. The references are `CLAUDE.md` (scope,
instruments, baseline), `docs/ROADMAP.md` (the ranking) and `docs/ISSUES.md` (every open case).

## What stands

The code review, code cleanup and documentation cleanup is done: four phases, no mechanism change.
Both ladder identity conditions read bit-identical after every step that touched `src/`, and the goldens
did not move. The first three phases are commit `3f7c5591`; the fourth (the tests) is its own commit.

| | before | after |
|---|---|---|
| permanent docs + `CLAUDE.md` + `README.md` | 10,645 lines | 6,046 lines |
| `scripts/design/` | 68 files, 29,038 lines | 43 instruments + `_shared.py`, 19,463 lines |
| `src/rigel/` Python | 38,295 lines, 202 dates, 9 doc citations | 37,665 lines, none of either |
| `tests/` | 165 files, 53,867 lines | 129 files, 52,151 lines; the same 2,670 tests |

The suite: **0 failed / 3,404 passed / 2 xfail, 3,406 collected**. The identity references frozen at
the start of the session are `~/Downloads/rigel_runs/arms/cleanup_identity_<condition>.json` for
`gdna_g05_ss_0.50_nrna_mid_capture_off` and `gdna_g05_ss_0.99_nrna_mid_capture_on`.

Structural changes a reader will meet: `strand_deconv` is folded into `gdna_strand`; the scripts load
each other through `scripts/design/_shared.py`; the tests take shared builders from
`tests/_index_builder.py` and `tests/_em_harness.py` (never `from conftest import …`, which is
ambiguous once a sub-directory has its own conftest) and the transfer-policy fixture from
`tests/calibration/conftest.py`; eight groups of test files were merged, with every test id kept.

## Owner decisions this session left open

1. **Deleted numbers with no home.** The agents cut measurements out of docstrings (bandwidth sweeps,
   old panel readings, profiling records, the v8 no-merge partition's justification in `splice_graph.py`,
   the `sweep_logodds_window` A/B, the closure failure rate in `result.py`). Git carries every one; none
   is recorded in `DESIGN.md` or `ISSUES.md`. Say which, if any, should be.
2. **A duplicate gate.** `test_an_ambig_flank_cannot_seed` and `test_an_AMBIG_flank_cannot_seed` in
   `tests/calibration/test_gdna_strand_fit.py` pin one rule on one fixture; the second is strictly
   stronger. The merge kept both because no test was to be lost.
3. **The two xfail reason strings** narrate history (a date, "the relay policy of the day"). They are
   executable records and were not touched.
4. **Not done, by judgement:** `splice_graph.py` and `index.py` are each four concepts and were not
   split; `strand_balance` stays separate from `gdna_strand` (a location estimand beside a dispersion
   one, and part of the pinned public surface); `cli._fragment_length_report` stays in the CLI (the
   `report/` package is behind an optional extra); `sim/net_flow.analyze_net_flow` is kept as that
   module's documented entry point.

## Worth knowing

- Earlier bulk renames damaged ordinary English in prose ("region_bounds" for *cuts*, "boundaries up"
  for *lines up*, "an boundary", "Region dependency" for *Node*). The ones found are fixed; run
  `rename_census.py --sense` before any future rename.
- A test that matches an assertion message by regex breaks when the message is reworded; one did in
  this session (the SILENCE gate in `test_sweep_backbone.py`) and was re-pointed at the sentence, not
  widened.

## What is next

`ROADMAP.md`'s ranking is unchanged by the cleanup: rank 1 is the ruler at zero gDNA
(`ISSUES: g00-shrinkage-upstream-repair`), then the rest of the pre-EM setup, the intron's own solve on
unstranded capture-OFF, the vertex atom, and the message policy only where a row is above the bar.
