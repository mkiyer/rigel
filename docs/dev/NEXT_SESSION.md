# NEXT SESSION — start here (2026-09-19, night: the capture physics is corrected and the ladder is rebuilt)

This file is only how to begin. The ranked view is `docs/ROADMAP.md`, the open problems are `docs/ISSUES.md`,
the rulings and the record are `docs/DESIGN.md` (§0b carries the scope and its 2026-09-19 amendment), what
"done" means is `docs/SUCCESS.md`, the lessons are `docs/TRAPS.md` cited by name, and the release procedure is
`docs/MANUAL.md` / `docs/PUBLISHING.md`.

## What changed today

THE SIMULATOR'S CAPTURE PHYSICS (owner ruling, 2026-09-19). A molecule hybridises through ONE contiguous part of
a probe: a transcript holding the junction a probe spans binds it whole, while gDNA, a nascent span and an isoform
without the junction hold its parts apart and bind the better one. That geometry is gDNA's whole disadvantage —
the `gdna_split_penalty` that bound a gDNA half-match at a fifth of the identical cDNA one is gone, and the loader
refuses a config that still carries the key. THE LADDER WAS DELETED AND REBUILT under it (23 GB, 16 conditions,
cached and certified, simulator gates 6/6); its capture-OFF rows are bit-identical to the retired ladder's, which
is the rebuild's own control, and the capture-ON identity reference was re-frozen (the capture-OFF one still
passes `--check` bit-identically).

THE DELIVERABLE, on the rebuilt ladder, is in `ROADMAP.md`'s claim: stranded × capture-ON now reads 6.8 / 4.5 /
7.8 / 109.8 % at `g00` / `g05` / `g50` / `g98`, where the retired ladder read 7.4 / 10.5 / 13.9 / 98.8 %. Most of
what ranked first this morning was the simulator's asymmetry (`TRAPS: prove-the-substrate`).

## What is ranked now (`ROADMAP.md`)

1. **The EM does not hold calibration's gDNA split** (`ISSUES: em-overturns-the-calibrated-gdna-split`) — the
   dominant in-scope residual and the whole of `g98`: the table reads 0.4793 against 0.50 at `g50 ss.99 ON` while
   calibration reads +0.9 %. Neither the prior, the ruler nor the gDNA length moves it. Its capture-OFF half is
   sized at the realistic nascent share first (`ISSUES: nascent-stress-sensitivity`).
2. **The capture ruler where no gDNA witnesses it** (`ISSUES: ruler-witness-geometry-on-transcript-panels`) — worth
   4.2 points of stranded capture-ON at `g00` and 1.1 at `g05`. It waits on an owner decision, because the only
   observable that sees isoform-specific capture is the probe design and Rigel reads no panel.
3. **The per-transcript allocation** (`ISSUES: per-transcript-prior-lane`) — true weights halve every stratum.
4. The pre-EM prior chain, now ranked by `prior_vs_oracle.py` rather than by the table.

## What is stale, and what it costs to fix

The physics changes only a panel whose probes span junctions. The test chromosome's benign panel re-simulates
BIT-IDENTICAL (checked), so its 30 conditions and the depth/fl/odg variants stand. Its junction-probed twin
(`scenarios_probes_junction`, 271 split probe blocks) and the two fl-gap side panels (they share the ladder's
panel) are STALE: re-simulating the twin is minutes, the fl-gap panels hours, and the fl-gap pair is already
stale on its nascent model (`ISSUES: flgap-panels-stale-nascent-model`).

## The five xfails are proven defects, each deferred to its thread

`ISSUES: two-sided-exon-row`; `ISSUES: antisense-prior-assembly-casualty`;
`ISSUES: the-lower-bound-noise-ratchet`; `ISSUES: nested-antisense-leak-under-the-sane-ruler` (two rungs).
Closing one means repairing the thing or asserting the invariant structurally, never widening a bound.

## What gates the release (`ROADMAP.md` item 8)

`PUBLISHING.md` is two commands and a wait; the STATE is what gates it. The deliverable measured and not
regressed per stratum; the zero controls at 0.000 and 1.000 (`zero_controls.py`, and the `g00` rung); the
suite at its standing count with an empty failure set; `preflight.py --full` green; the standing risks re-read
(`ISSUES: capture-degeneracy-standing-risk`, `ISSUES: flgap-panels-stale-nascent-model`); and the manual true
of what ships.

## Standing rulings carried (unchanged)

* Float64 for the whole of ψ; ONE solver; ONE λ lattice; no θ lattice; the tilt's hypothesis space is
  {pure +, pure −, mixed}. The strand channel's liveness is a protocol decision on the spliced 2×2.
* The landscape trains only where a solve locates a slot; the ruler reads its located enriched mode and admits
  failure where there is no gDNA to read.
* Real data is a test input, never a design input. A few high-quality instruments, kept current. The source
  cites no doc. One production path: once native code is validated the Python it replaces is deleted.
* Unstranded × capture-ON stays DEFERRED — reported on every benchmark, never a development target, never
  ranked on a pooled total. The fragment-length composition channel stays retired.
* The owner drives commits. The refit count and CI's on-demand trigger are the owner's.

## Where everything is

* The synced scratchpad: `~/Downloads/rigel_runs/prototypes/2026-09-17_port_ii/` — `commits/committed/` (every
  snapshot, all landed), `s16/`–`s18/` (the port and the block), `s19/`–`s21/` (the performance campaign, its
  attributions and the taper study that the parked thread resumes from).
* Captures and reports: `perf/sweeps_VCaP_step19` (the sweep replay's capture, bit-identical on this tree);
  `perf/plan_final_2026-09-19/` is the deep run's current before-and-after.
* THE REBUILT LADDER's results: `~/Downloads/rigel_runs/suite/ladder/arms/` (the four `quant_accuracy` arms plus
  `oracle_alloc_seed`, and `calibration_vs_oracle.json`) with every stage's log in
  `~/Downloads/rigel_runs/logs/ladder_*.log` and the chain that produced them in `ladder_rebuild_2026-09-19.sh`.
* THE BASELINE on the RETIRED ladder: `~/Downloads/rigel_runs/arms/2026-09-19_e2e_baseline/` (sampled assignment, the per-stratum
  `qa_report.txt`, `decomposition.txt`, `alloc.txt`). THE STRANDED × ON DISSECTION:
  `~/Downloads/rigel_runs/arms/2026-09-19_stranded_on/` — `all_scenarios.txt` (the fractional panel),
  `stranded_on_arms.txt` (every ruler and prior arm), `tables/` (per-transcript tables per arm), `testchr/` (the
  benign-vs-junction control), `confusion_*` (per-fragment truth against assignment at `g50`), and the scratch
  runner `dissect_run.py` / `dissect_analyze.py` / `confusion.py` it came from.
* The identity references: `~/Downloads/rigel_runs/arms/review_identity_*.json`, BIT-IDENTICAL on this tree.
  They are the gate for any change that must not move a number.
