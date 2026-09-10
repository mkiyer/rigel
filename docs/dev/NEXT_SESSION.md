# NEXT SESSION — the tool re-measured under the landscape prior's two landings (2026-09-10, uncommitted on `main`); NEXT IS THE RULER AT ZERO gDNA (handoff)

⭐⭐⭐ **THE REFERENCES:** `CLAUDE.md` (the scope, the instrument table, the standing baseline),
`ROADMAP.md` (the state and the ranking, both re-audited 2026-09-10), `ISSUES.md` (every open case;
`the-landscape-training-population-arms` in CLOSED / REFUSED with every arm's number), `DESIGN.md` §7.1
(the prior's three rulings and their measurements). This file is the state.

## WHAT STANDS (2026-09-10; the working tree, NOT committed — the owner drives commits)

1. **Two mechanisms landed on the gDNA landscape prior, one gate file** (`DESIGN.md` §7.1;
   `tests/calibration/test_landscape_training_population.py`, 10 gates, nine perturbations watched fire):
   a node whose only evidence is a bound does not train it (`RegionBelief.informed`); its location-free
   kernels (count < 1) are placed by the previous refit's landscape (`landscape._estep_kernels`); its grid
   spans every region and boundary (`fit_landscape(domain=…)`, the owner's ruling). The ladder's zero
   controls went from ~150k invented fragments to a few hundred; in scope 0.95–1.00×.
2. **THE SUITE: 0 failed / 3,609 passed / 2 xfail, 3,611 collected** (`CLAUDE.md` carries the +16
   accounting; `docs/dev/PLAN_measured_prior.md` retired as superseded, git carries it).
   `preflight.py --full` 17/17. Goldens regenerated; the magnitude is in `DESIGN.md` §7.1.
3. **THE NEW INSTRUMENT:** `scripts/design/landscape_training_census.py` (row in `CLAUDE.md`'s table).
4. **THE RECORDS:** the prototype harness and every arm's output under
   `~/Downloads/rigel_runs/prototypes/2026-09-10_landscape/`; the four instrument re-runs under
   `~/Downloads/rigel_runs/arms/2026-09-10_post_landscape/` (`calibration_vs_oracle.txt`, `walk/`,
   `v_base.jsonl` + `v_free.jsonl` + `v_compare.txt`, `qa_base.jsonl` + `qa_oracle.jsonl` + `qa_report.txt`).

## THE STATE, RE-MEASURED (every number below is from those four runs; re-derive, never quote onward)

* **`calibration_vs_oracle.py`, the 0.8.0 metric** (mass-weighted |Δ gDNA share| per object, region /
  boundary): stranded OFF 0.0067 / 0.0103, stranded ON 0.0102 / 0.0150, unstranded OFF 0.0085 / 0.0122
  — unchanged from the merge day within a few ten-thousandths; the `g00` zero controls **0.0000 / 0.0000**
  (828 and 921 fragments of 18M and 26M; they were 0.0081 / 0.0077); the deferred stratum 0.0717 / 0.1185.
  ⭐ **The RULER** (③): factor P 0.150 against O 1.000 on the `g00` rows — 996M bp over 51,543 transcripts —
  and P/O 1.033 / 1.043 on the in-scope capture-OFF strata where the instrument's own contract says 1.000.
* **`calibration_walk.py`** (C = local solve; the refit alone; messages on top): at every zero row the
  prior alone reaches a few hundred fragments and messages add a few hundred of noise; in scope the prior
  does most of the work and messages still remove 13k–129k on the stranded capture-ON and `g98` rows; on
  the deferred stratum the refit alone is +2.6M / +6.3M and the messages −4.9M / −9.9M (unchanged).
* **`vertex_ceiling.py`** (`vertex_free` against `base`, the final answer): within 1 % on every stranded
  in-scope row; 2 % / 2 % / 7 % of the unstranded capture-OFF rows at `g05` / `g50` / `g98`; 24–39 % of the
  deferred stratum; the zero rows 325 → 13. Unchanged in scope from the merge day.
* **`quant_accuracy.py`, the thermometer** (misassigned fragments, transcript level, `base` → `oracle`
  prior): stranded OFF 1.03×, stranded ON 0.98×, unstranded OFF 1.02× — a perfect prior is worth nothing
  in scope; the deferred stratum 0.73×; the `g00` rows 9.64M under both arms, the largest of any stratum
  and untouched by the prior — the ruler.
* **`policy_benchmark.py --by-class`** (the final table, `ladder_final.txt` in the prototype record):
  unstranded OFF `g50` — the intron class 46 % of the row at 2.3 % of its own fragments; stranded ON `g50`
  — exon|exon boundaries 49 % and walled exons 19 %. Same shape as the merge day.

## WHAT IS NEXT — the owner's decision (2026-09-10): a CODE REVIEW, CODE CLEANUP AND DOCUMENTATION CLEANUP session comes first; its kickoff prompt is `docs/dev/CLEANUP_SESSION_PROMPT.md`. After it, `ROADMAP.md`'s ranking, re-audited

1. **THE RULER AT ZERO gDNA** (`ISSUES: g00-shrinkage-upstream-repair`, re-priced: the composition it reads
   is now right and the factor is still 0.15, because `capture_eff_length._global_reference_density`
   detects a reference from any five slots with positive mass). The fix is the detector. First step:
   settle whether the instrument's "exactly 1.000 off capture" contract is stale (both P and O read
   0.92–0.97 there), then derive the enrichment test, prototype outside `src/`, judge on ③ per stratum.
2. The rest of the pre-EM setup (rank 2; re-run `prior_vs_oracle.py` first).
3. The intron's own solve on unstranded capture-OFF (rank 3); the vertex atom (rank 4); the message
   policy above the bar (rank 5).

## THE LESSONS THIS SESSION PAID FOR

* **Re-measure the whole page, not the number you changed.** The composition metric said "done"; the
  ruler on the same page said the largest in-scope defect on the metric is now somewhere else entirely.
* **An open entry's prescription can be completed and its defect survive.** "Repair the composition, not
  the function" was followed to 828 fragments and the factor did not move; the entry is re-priced, not
  closed, and the repair moved to the detector.
* **A census in the estimator's own currency before a mechanism; a decomposition by contributor before a
  derivation.** Both mechanisms this session followed from a table naming who carried the mass.
* **A prototype's evidence class is not the src predicate until it is checked slot by slot** (1,476 own-flux
  ceilings; 137k against 111k). **Diff the goldens column by column** (the tiny toys, not the ladder,
  exposed the grid collapse and the guard). **A shuffle or a mirror control only discriminates where the
  population it acts on is mixed.** **zsh arrays are 1-indexed** — check every sharded table's row count.
