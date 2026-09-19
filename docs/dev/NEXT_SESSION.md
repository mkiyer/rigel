# NEXT SESSION — start here (2026-09-19, late: stranded × capture-ON is root-caused; the repair waits on an owner decision)

This file is only how to begin. The ranked view is `docs/ROADMAP.md`, the open problems are `docs/ISSUES.md`,
the rulings and the record are `docs/DESIGN.md` (§0b carries the scope and its 2026-09-19 amendment), what
"done" means is `docs/SUCCESS.md`, the lessons are `docs/TRAPS.md` cited by name, and the release procedure is
`docs/MANUAL.md` / `docs/PUBLISHING.md`.

## Where the tool is

The end-to-end baseline is measured (`ISSUES: end-to-end-error-unattributed`, CLOSED) and the owner made
stranded × capture-ON the focus, with a target below 5 %. Every benchmark now runs under fractional assignment
(`quant_accuracy.py --set em.assignment_mode=fractional`, owner 2026-09-19; the report refuses to mix modes).

THE ROOT CAUSE (`ISSUES: ruler-witness-geometry-on-transcript-panels`): the ladder's panel places probes along
transcripts, so 24 % of them span a splice junction, and such a probe captures only the isoforms that hold the
junction. Capture becomes isoform-specific, and the EM splits a gene's shared fragments by the ratio of its
isoforms' capture-aware lengths. The shipped ruler reads capture from gDNA, which has no junctions, and at zero
gDNA it has nothing to read, so on this panel it costs more than it corrects (stranded ON `g05` 10.5 %, 7.8 % with
it switched off). The simulator's own length (`quant_accuracy.py --arm oracle_ruler`) takes the stratum to 2.9 /
3.6 / 8.4 % at `g00` / `g05` / `g50`. The test chromosome is the control: its benign panel reads capture-OFF
levels, its junction-probed twin 35–41 % at every gDNA level. What remains at `g50` is the EM's gDNA under-call
under capture (`ISSUES: em-overturns-the-calibrated-gdna-split`).

## The decision it waits on

The one observable that sees isoform-specific capture is the probe design, and Rigel reads no panel
(`DESIGN.md` §7.2). So before any build: do the panels Rigel will meet span junctions (exome-style panels on
genomic exons do not; transcript-designed ones do), and does Rigel take the probe design as an input? Both are the
owner's. Everything a candidate needs to be judged is in place: `oracle_ruler` is the ceiling, `ruler_vs_truth.py`
scores a ruler per transcript (read its within-gene spread, `TRAPS: judge-a-ruler-by-its-within-gene-spread`), and
the test chromosome's two panels are the controlled pair.

Then `ROADMAP.md` items 2–4: the per-transcript prior lane, the EM's gDNA split in both directions (the
capture-OFF half sized at the realistic nascent share first), the rest of the prior chain.

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
* THE BASELINE: `~/Downloads/rigel_runs/arms/2026-09-19_e2e_baseline/` (sampled assignment, the per-stratum
  `qa_report.txt`, `decomposition.txt`, `alloc.txt`). THE STRANDED × ON DISSECTION:
  `~/Downloads/rigel_runs/arms/2026-09-19_stranded_on/` — `all_scenarios.txt` (the fractional panel),
  `stranded_on_arms.txt` (every ruler and prior arm), `tables/` (per-transcript tables per arm), `testchr/` (the
  benign-vs-junction control), `confusion_*` (per-fragment truth against assignment at `g50`), and the scratch
  runner `dissect_run.py` / `dissect_analyze.py` / `confusion.py` it came from.
* The identity references: `~/Downloads/rigel_runs/arms/review_identity_*.json`, BIT-IDENTICAL on this tree.
  They are the gate for any change that must not move a number.
