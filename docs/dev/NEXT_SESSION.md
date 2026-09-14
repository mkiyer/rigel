# NEXT SESSION — start here (2026-09-14, after the code-review and cleanup session)

This file is only how to begin. The port is `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` §F; the θ thread's
derivations are `docs/dev/THETA_QUADRATURE.md`. Everything is committed and pushed; `e1358c55` is the last
change to code or instruments, and the handoff commit follows it.

## The agreed order of the sessions (owner, 2026-09-14)

1. ~~Code review and cleanup~~ — DONE 2026-09-14 (what it left is below).
2. **THIS ONE — the test chromosome's new structures**, then both panels remeasured: the long single-exon
   transcript encompassing a multi-exon transcript on the opposite strand; the single-exon antisense gene
   wholly inside a sense exon (the tilt atom's accepted limit, so the panel carries its number); the shared
   exon of two spliced genes on opposite strands (the deep stress); a single-exon gene inside the opposite
   strand's intron (the majority AMBIG class); head-to-head genes with overlapping UTRs. Edit
   `scripts/sim/test_reference/test_chr.yaml`, render, rebuild and re-certify by `docs/TESTING.md` §0a
   (`panel.py status --config scripts/sim/configs/test_reference.yaml` names each stage); then
   `calibration_vs_oracle.py` and `policy_benchmark.py --by-class` on both panels become the standing
   numbers (`DESIGN.md` §7).
3. **The ruler at zero gDNA** (`ISSUES: g00-shrinkage-upstream-repair`, the largest in-scope number on the
   metric page, the modal real case) and **the flux price's witness in whole-strand units**
   (`ISSUES: flux-price-witness-units`, the one accuracy item inside the port's unit). One mechanism at a
   time, each its own commit.
4. **The performance re-baseline** (two back-to-back profiler pairs at 8 threads, identity references frozen
   fresh on that tree, `sweep_replay.py capture`) **and THE PORT** (plan §F).

## The prompt for the next session (paste as the first message)

> Start from `main` (clean, at the handoff commit after `e1358c55`, or later). Read `CLAUDE.md`, then
> `docs/dev/NEXT_SESSION.md`. Run `scripts/design/preflight.py --full` and the suite before touching
> anything; `CLAUDE.md`'s baseline line is the count to reproduce (3,374 passed / 3 xfail / 3,377 collected).
>
> THE STANDING RULINGS are unchanged: elegant, simple, efficient, clear, concise, maintainable code; converge
> and delete — no legacy, no compatibility shims, no speculative code; the suite's count re-derived, never
> adjusted; each item its own commit, committed on my go; never `ruff format scripts/`.
>
> THIS SESSION IS THE TEST CHROMOSOME'S NEW STRUCTURES, and nothing lands in `src/`: a defect a structure
> exposes is recorded for session 3. One structure at a time, each its own commit: mirror it from the ladder
> locus it stands for, as the YAML's block comments do; add it to `test_chr.yaml` with its renders; rebuild by
> `docs/TESTING.md` §0a; and confirm on certified truth that the substrate carries the case — the new slots
> by class, never a panel number. `rename_identity.py --check` against the three `review_identity_*`
> references stays BIT-IDENTICAL throughout: they are ladder conditions and a real library, so a change that
> moves them is a change to code. When every structure is in, re-simulate, re-cache and re-certify the test
> chromosome's configs, then run `calibration_vs_oracle.py` and `policy_benchmark.py --by-class` on BOTH
> panels — the new standing numbers, every in-scope stratum, both zero controls and the deferred stratum
> reported apart. Run `preflight.py --full` first and last. End with the handoff for session 3.

## Practical notes for the test chromosome

* The YAML is the one hand-edited file; `build_test_reference.py` renders the GTFs, the abundance TSV and the
  probe panels into the repo and the runs directory, `--check` says whether the checked-in renders match, and
  `tests/test_test_reference_renders.py` refuses a drifted render — commit renders with the YAML edit.
* Each block's comment records what it did to the chromosome length and the fragment budget, holding the gDNA
  density per base and the RNA depth per existing transcript near constant; a new block follows the same rule
  and says so (the both-stranded and sj+terminus blocks are the pattern). Strands stay balanced (the builder
  refuses a difference of more than one) and a multi-isoform structure follows the replication rule.
* `panel.py status` now reads the index the config names, so the recipe needs no `--index`. The seven
  test-chromosome configs share one reference, index and abundance file and differ only in their outdir and
  the one axis each names; TESTING §0a says which to rebuild.
* The metric's certification is `calibration_oracle.py` (its six named gates; `slot_truth.npz` per
  condition), which `panel.py cache` runs; the census of the new slots by class can be read off it and off
  `policy_benchmark.py --by-class`.

## What the cleanup session left

* **Commits** `95bedc8e` … `e1358c55` (34 on main, each gated): 21 instruments retired by the owner's ruling
  (a few high-quality instruments kept current; CLAUDE.md's table is the set), the certifier's real-data
  crash since W11 repaired, dead and legacy code removed from `src/`, the native layer (one accumulator
  finalizer, the EM's never-read `locus_enable_gdna`, dead helpers, the unused SIMD runtime detection), tests
  and instruments; every false or stale docstring the census found rewritten against the code; the source
  cites no doc; `DESIGN.md` §6b.15 split into thirteen numbered sub-rulings; three TRAPS retired.
* **No number moved.** Every step: `rename_identity.py --check` 3/3 and `sweep_replay.py replay` on
  `sweeps_MO_3021_step6` calls 0–3, BIT-IDENTICAL; the suite re-derived at each retirement.
* **The identity references are `review_identity_*`** (`~/Downloads/rigel_runs/arms/`), frozen on
  `84923136` — the previous handoff's `port_identity_*` described the L5 tree, not `84923136`, and are
  deleted. `sweeps_MO_3021_step6` still replays bit-identical; the landscape's location floor landed after it
  was captured, so re-capture in session 4 as planned.
* **Recorded, not cut:** `ISSUES: hygiene-ledger` (a) rotten-but-live instrument paths, (b) vacuous gate
  clauses, (c) instrument duplicates, (d) claims not re-derived on the current tree — among them `fl`'s
  crossing pools called "gDNA by structure", which unspliced RNA contradicts.
* **CI** runs on demand only (`.github/workflows/ci.yml`, owner, 2026-09-14), so a push is untested by GitHub.
  The last three automatic runs before the postponement had failed and were not diagnosed; run it once from
  the Actions tab before restoring the triggers.
* **Disk:** the stale identity references and sweep captures, `test_reference_STALE_193genes_2026-09-08` and
  the stale `rigel 0.1.0` editable-install files were deleted; kept are `arms/review_identity_*`, the two dated
  `arms/` result directories, `perf/{sweeps_MO_3021_step6, baseline_2026-09-14, ab_locus_2026-09-11,
  ab_memo_2026-09-11}`.

## The two pricing questions inside the port's unit — what the owner decides

Both live in `messages/transfer_rows.py` and `messages/lanes.py`, the port's first step, so a fix after the
port is written twice; the decision is fix-before or accept-and-record.

* **`ISSUES: flux-price-witness-units`.** A junction's certified flux is a lower bound on the exon's RNA
  level; its price compares the junction's route rate (whole-strand units) with the exon's own count on the
  read column, which holds only `(1 − κ')` of the strand's RNA plus half the gDNA, so every flux ceiling pays
  `log(1 − κ')²` nats² of spurious disagreement — 0.48 on unstranded data, an in-scope stratum. The naive fix
  was REFUSED; the owed form is the strand's RNA count in whole-strand units. Recommendation: derive and A/B
  it in session 3, before the port carves the price.
* **`ISSUES: the-lower-bound-noise-ratchet`** (the encompassing flank's xfail). A level is delivered at the
  sender's sampled density, so a neighbour whose sample ran high pins the flank above the truth (0.8 % of one
  ladder row). Recommendation: accept and record; the port carves the level rule as it is.

## Decisions on record

* Float64 for the whole of ψ; ONE solver; ONE λ lattice at a dimensionless step with a stated guarantee; no θ
  lattice anywhere (the tilt count is derived); the tilt's hypothesis space is {pure +, pure −, mixed}.
* The strand channel's liveness is a protocol decision on the spliced 2×2; gDNA enters nowhere.
* The landscape trains only where a solve locates a slot.
* A few high-quality instruments, kept current; no suite gate polices instruments (owner, 2026-09-14).
* The source cites no doc (owner, 2026-09-14).
* `EMConfig.warm_start`'s `prior` and `uniform` arms stay (the parked per-transcript prior lane's
  foundation); the `gdna_none_` condition names stay (the example config uses them).
* The message cache's on/off switch, the refit count and the scan's thread split are the owner's; parallelism
  waits for the port.

## The session scratchpads (persist across sessions; nothing in the tree cites them)

This session: `/private/tmp/claude-503/-Users-mkiyer-proj-rigel/60889fc2-ba0b-410a-ab2e-ebff583aed2d/
scratchpad/cleanup/` — `CENSUS.md` and `LEDGER.md` (the rulings, with greps), `agents/` (the instrument,
TRAPS and docstring-claim censuses), `cov/` (the coverage data), `gates/` (every item's verdict),
`messages/`, `cites/` (every edit batch). The previous session's
`…/4ec3e3a5-268d-4828-a6e7-91f9cd89647e/scratchpad/` holds L3, L5 and the location floor's material.
