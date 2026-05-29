# PR 0 — Burn the residue *(Step 1: remove stale/legacy/dead code)*

**Parent plan:** [`../00_implementation_plan.md`](../00_implementation_plan.md) §7.
**Type:** Python-only deletion + docstring fixes. No behavior change.
**Build required:** no.

## Goal

Delete every remnant of v5 calibration that still ships
([parent §2.3](../00_implementation_plan.md)). After this PR the tree
contains *no* dead calibration-prior machinery — but `rigel quant` still
aborts at `calibrate()` (`NotImplementedError`), exactly as before. This
PR is a clean, reviewable deletion that de-risks every PR after it.

**Non-goals.** No new calibrator code; no reorganization (that's PR 1); no
golden regeneration (that's PR 8, gated on real output changes). The sim
calibration reports are **stubbed**, not rewritten (rewrite is PR 8).

> Line anchors below are from `4124276` and **will drift** — every task
> says what to `git grep` for. Match on symbol, not line number.

## Pre-work — capture the pre-burn baseline (Decision D4)

Before any deletion, snapshot the synthetic baseline PR 7 compares against
(the live `main` no longer runs end-to-end, so it must come from a
pre-burn SHA — candidate `fc96902 "new calib system implemented"`):

```bash
git worktree add /tmp/rigel-preburn fc96902
# in that worktree: build, then run the synthetic sweeps into the MAIN repo's scratch/:
#   pip install --no-build-isolation -e .   (in /tmp/rigel-preburn)
#   for cfg in scripts/benchmark/configs/locus_simple_*.yaml; do
#     python scripts/sim/locus_sweep.py -c "$cfg" -o "/Users/mkiyer/proj/rigel/scratch/preburn/$(basename "$cfg" .yaml)"
#   done
git worktree remove /tmp/rigel-preburn
```

Record the resulting path and the exact SHA used in
`scratch/preburn/PROVENANCE.txt`. If `fc96902` does not build/run cleanly,
walk back to the nearest "all tests passing" SHA (`7a1c700`) and note it.
**This pre-work is the one part of PR 0 that must not be skipped** — the
baseline window closes once we delete more.

## Deletions

### T1 — `pipeline.py`: delete the dead v5 EM-prior driver

Verify no live caller first:
```bash
git grep -n "_run_locus_em_partitioned\|_build_locus_meta\|_calibration_strand_summary" src/ tests/
```
Expected: definitions only, no call sites (the pipeline aborts at
`calibrate()` before reaching them; `quant_from_buffer` is a stub).

Delete:
- `_run_locus_em_partitioned` — the entire function (`def` at ~546 through
  just before `def quant_from_buffer` at ~959). This includes its nested
  `_build_locus_meta` (~588) and all 11 `prior_*` parameters (~560–570),
  the `prior_*` dict population (~651–661), and the `prior_*[lid]`
  call-throughs (~780–805, ~899–926).
- `_calibration_strand_summary` — the function at ~190–196 (stop before
  `_warn_if_calibration_strand_unidentifiable` at ~197, which is **KEEP**).

**KEEP** (PR 6 reuses these): `_setup_geometry_and_estimator`,
`_score_fragments`, `_assign_locus_ids`, `_populate_em_annotations`
(~406–543), `_warn_if_calibration_strand_unidentifiable` (~197),
`_strand_summary_identifiable`, `scan_and_buffer`,
`_wire_calibration_regions`.

> **Note for PR 6.** PR 6 will rebuild a lean locus-EM driver. The
> non-v5 partition→EM wiring inside the deleted `_run_locus_em_partitioned`
> (building `partition_tuples`, calling `estimator.run_batch_locus_em_partitioned`,
> collecting per-locus results) is salvageable from
> `git show 4124276:src/rigel/pipeline.py`. We delete the whole function
> now (it is unreachable dead code) and re-harvest the clean parts in PR 6,
> rather than half-delete a function in PR 0.

### T2 — `estimator.py::get_loci_df`: drop the 11 `prior_*` output columns

```bash
git grep -n '"prior_' src/rigel/estimator.py
```
Delete from `get_loci_df`:
- the 11 `"prior_*"` entries in the `cols` list (~844–854),
- the per-row `r.get("prior_*", ...)` extractions (~914–924),
- the corresponding row-dict population (~953–963).

Leave the EM core and all non-`prior_*` columns untouched. These columns
are already constant-zero dead weight (their only producer was T1's
`_build_locus_meta`).

### T3 — `tests/test_cli.py`: delete the 3 legacy CLI tests

Delete (flags removed from `cli.py`; fields removed from
`config.CalibrationConfig` — these tests can only fail):
- `test_cal_quality_thresholds_flow_to_config` (~312–320)
- `test_cal_fl_scoring_prior_ess_flows_to_config` (~322–329)
- `test_yaml_cal_fl_scoring_prior_ess` (~331–337)

### T4 — `tests/test_golden_output.py`: drop the legacy golden column

Remove `"prior_ess_final"` from `_LOCI_NUMERIC_COLS` (~197). (The other
`prior_*` columns are not individually listed there; confirm with
`grep -n prior_ tests/test_golden_output.py`.) Golden *files* are not
regenerated here — that is PR 8, after real output changes exist.

### T5 — `sim/`: stub the dead `summary.json` / `CalibrationResult` consumers

These read `CalibrationResult.to_summary_dict()` / `.fl_models` and old
`summary.json` calibration keys that no longer exist — they crash on
access today. **Stub, don't rewrite** (the v6-schema rewrite is PR 8):
- `sim/locus_sweep.py` (~893–926): replace the calibration-report block
  with a guard that emits `"calibration report unavailable (v6 calibrator
  pending — see docs/acc_caljointmodel/)"` and skips the `cal.*` access.
- `sim/analysis.py` (~400–421, ~1163–1439): wherever it does
  `cal.get("region_calibration"|"fl_models"|"density_evidence"|...)`,
  guard so a missing/None calibration block yields empty/`None` report
  fields instead of touching deleted attributes. (Most already use
  `.get(...)`; the task is to ensure none assume the old nested schema and
  none reference `fl_models`/`to_summary_dict`.)

### T6 — `tests/test_sim_analysis.py`: neutralize the old-schema fixture

The `_write_summary` fixture (~42–67) builds a synthetic `calibration`
block with `density_evidence`/`fl_models`. Delete or `pytest.mark.skip`
the tests that assert the old v5 `summary.json` calibration schema,
consistent with T5. (Keep any non-calibration sim-analysis tests.)

## Docstring / doc-reference fixes

- `scoring.py:139–143` — drop the `CalibrationResult.fl_models.rna_scoring`
  / `.gdna_scoring` lines; describe the actual `FragmentLengthModel` args
  the function takes.
- `locus.py:11–12` — the docstring points at `rigel.calibration.assemble_priors`,
  which does not exist yet. Point it at "`rigel.calibration.priors`
  (lands in PR 6)" or drop the line.
- `pipeline.py` — module docstring (~13–20, "assemble fused regional gDNA
  priors") and the stale `CalibrationScanPayload` type comments (~85, ~238)
  → `AccumulatorPayload`. Also fix the `quant_from_buffer` stub message
  (~982, "assemble_em_inputs replacement") to reference PR 6.
- `config.py:217` and `calibration/__init__.py:11` — `04_outputs.md`
  (nonexistent) → `04_interface_contract.md`.
- `cli.py` (~1162–1163) — drop the "unbiased v5 prior" wording (minor).

## CHANGELOG

```
2026-05-29: Completed the calibration-v5 burn-down — removed the dead
locus-prior EM driver (pipeline._run_locus_em_partitioned), the 11
prior_* loci-df columns (estimator.get_loci_df), the dead
_calibration_strand_summary helper, the 3 legacy --cal-* CLI tests, and
stubbed the sim summary.json consumers left behind by the Phase-A burn.
Replacement is the v6 calibrator (docs/acc_caljointmodel/).
```

## Acceptance gate

- [ ] `python -c "import rigel"` succeeds.
- [ ] `pytest --collect-only tests/` is clean (no import/collection errors).
- [ ] Substrate suite green: `pytest tests/native/ tests/test_accumulator_payload.py tests/test_scanner_accumulator_integration.py`.
- [ ] Zero hits (exclude the golden *data* fixtures, which still embed the
      old column schema and are regenerated in PR 8):
      `git grep -nE 'prior_n_local_gdna|prior_rna_share_v5|prior_ess_final|prior_flags|prior_unspliced_total|prior_locus_weight|prior_shrink_weight|_run_locus_em_partitioned|_build_locus_meta|_calibration_strand_summary|to_summary_dict|\.fl_models' -- src/ tests/ ':(exclude)tests/golden/*'`
- [ ] Zero hits: `git grep -n '04_outputs.md' src/`
- [ ] `git grep -n 'CalibrationScanPayload' src/` → zero hits.
- [ ] `sim/` imports cleanly: `python -c "import rigel.sim.locus_sweep, rigel.sim.analysis"`.
- [ ] `rigel quant` on a scenario BAM still scans and then aborts at the
      `calibrate()` `NotImplementedError` (unchanged behavior; nothing dead
      in between). Equivalent: the *only* end-to-end failures are the
      pre-existing `calibrate()` aborts — no new failure modes, no
      `AttributeError`/`NameError` from a dangling reference.
- [ ] `scratch/preburn/PROVENANCE.txt` exists (pre-work / D4).
- [ ] `ruff check src/ tests/` no new findings vs baseline (deletions
      shouldn't introduce lint).

## Rollback

Revert the PR. Pure deletion; the accumulator substrate is untouched, so
the substrate test suite stays green regardless.

## Review checklist (for critique before implementation)

- Is deleting `_run_locus_em_partitioned` wholesale (vs half-deleting the
  `prior_*` machinery and keeping the partition→EM skeleton) the right
  call? (Plan's position: delete wholesale, re-harvest clean parts in PR 6.)
- Is the `fc96902` baseline SHA acceptable for D4, or prefer `7a1c700`?
- Stub-vs-delete for the `sim/` calibration reports — acceptable to defer
  the real rewrite to PR 8?
- Any KEEP item above that you actually want gone now?

## Execution notes (as implemented — 2026-05-29)

Implemented on branch `claude/pr00-burn-residue`. All acceptance-gate
checks pass. Deviations from the plan above, with rationale:

1. **Baseline capture deferred (D4).** The snapshot is pinned to an
   immutable SHA (`fc96902`), so the window does not close as later PRs
   land; building the pre-burn binary now would clobber the dev editable
   install. Provenance + exact procedure recorded at
   `scratch/preburn/PROVENANCE.txt`; the capture runs as the first step
   of PR 7.
2. **T5 reduced to `sim/locus_sweep.py` only.** That file had the hard
   attribute access (`cal.to_summary_dict()`, `cal.fl_models.gdna`) the
   gate flags — stubbed (keeps the FL truth columns, blanks the
   calibration estimates). `sim/analysis.py` reads the old schema only via
   null-safe `.get(...)` (not flagged by the gate) and is fully deferred
   to the PR 8 v6-schema rewrite.
3. **T6 not needed.** `tests/test_sim_analysis.py` is self-consistent (it
   writes its own synthetic old-schema `summary.json` and asserts
   `analyze_calibration` parses it); `analysis.py` is unchanged, so those
   5 tests still pass. Rewritten alongside `analysis.py` in PR 8.
4. **Bonus cleanup.** Deleting `_run_locus_em_partitioned` orphaned
   `import os` (F401) — removed. Also fixed a PRE-EXISTING `F841`
   (`calibration = calibrate()` in `run_pipeline`, present on `main`) by
   dropping the dead assignment and refreshing its stale "Phase D/E"
   comment. Repo is now fully `ruff`-clean.

### Acceptance-gate results
- `import rigel`, `import rigel.sim.*` — OK.
- `pytest --collect-only` — clean (867 tests; −3 deleted CLI tests).
- Substrate suite — 31 passed.
- Touched test files — `test_cli.py` 33, `test_estimator.py` 54,
  `test_sim_analysis.py` 5 — all pass.
- Grep gate (excl. `tests/golden/*`), `04_outputs.md`,
  `CalibrationScanPayload` — zero hits.
- `ruff check src/ tests/` — All checks passed.
- Full suite — 619 passed; 233 failed + 15 errors, **all** tracing to the
  single `calibrate()` `NotImplementedError` abort (248 `E`-lines, 100%
  `NotImplementedError`; zero new failure modes). End-to-end green returns
  at PR 6; golden fixtures regenerate at PR 8.
