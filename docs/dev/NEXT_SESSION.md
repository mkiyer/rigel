# NEXT SESSION — start here (2026-09-19, evening: the end-to-end baseline is measured; the EM owns the residual)

This file is only how to begin. The ranked view is `docs/ROADMAP.md`, the open problems are `docs/ISSUES.md`,
the rulings and the record are `docs/DESIGN.md` (§0b carries the scope and its 2026-09-19 amendment), what
"done" means is `docs/SUCCESS.md`, the lessons are `docs/TRAPS.md` cited by name, and the release procedure is
`docs/MANUAL.md` / `docs/PUBLISHING.md`.

## Where the tool is

The tree is `cffca248` plus a docs-only change: the baseline's numbers, the re-ranked roadmap, the issues
it closed and opened. No `src/` line moved in the measurement session. The suite stands at 3,416 passed / 5
xfail, and `preflight.py --full` is green.

THE BASELINE (`ISSUES: end-to-end-error-unattributed`, CLOSED, has every number): in scope the transcript
table misassigns 4.4 / 3.9 / 12.8 % of the RNA (unstranded OFF / stranded OFF / stranded ON), 150 / 78 / 519 times its
reseed floor, and a perfect prior recovers only 3–5 % of that — almost all at `g98`, and nothing measurable at
`g05`. What survives splits in two, and neither half is the prior chain's:

* the per-transcript allocation — truth as the weights removes about half of the transcript error in every
  stratum, the gDNA-free `g00` rows included (`ISSUES: per-transcript-prior-lane`);
* the capture-OFF gDNA over-call — calibration's library split is right (`g05 ss.50 OFF` 0.049
  against 0.05) and the table's is not (0.104), under a perfect prior and under true allocation weights alike
  (`ISSUES: em-overturns-the-calibrated-gdna-split`), measured only at the nascent stress share.

## The next build: the per-transcript prior lane (`ROADMAP.md` item 1)

`rna_prior_weight` is built end to end and `pipeline.py` omits it. A wiring gap plus a support decision; two
weightings are already refused with their numbers, and the next candidate is a sparsity mechanism
(`ISSUES: per-transcript-prior-lane`). DERIVE → DESIGN → PLAN → PROTOTYPE → A/B before `src/`, judged on the
deliverable per stratum above `--arm base_reseed`, with both zero controls. `oracle_alloc_seed` is the
capability proof, never the headroom: it hands over the true support.

Then `ROADMAP.md` item 2 — size the capture-OFF over-call at the realistic nascent share
(`ISSUES: nascent-stress-sensitivity`) before anything is built on it — and item 3, the rest of the prior
chain, which the oracle arms price near zero on the ladder's table.

## ⛔ Two instrument defects the baseline found — repair before the next oracle arm

* `ISSUES: oracle-cache-key-hashes-a-thread-count` — `109d8aac` moved `bgzf_threads`' default and the cache
  digest hashes it, so `quant_accuracy.py` refuses every oracle arm, and `prior_vs_oracle.py` /
  `pass0_vs_oracle.py` quietly re-scan and overwrite the certified oracle caches. The repair is in
  `scan_cache.py` (the owner's call) and is a no-op on every number. Until it lands, the oracle arms run through
  the key-only wrapper `qa_keyed.py` beside the baseline's data (below); do NOT run `prior_vs_oracle.py`
  against the shared caches.
* `ISSUES: oracle-ruler-arm-cannot-reach-the-ruler` — the arm swaps count arrays the ruler stopped reading at
  `c44fc306`, so it refuses every condition. Repair (an oracle efficiency) or retire.

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
* THE BASELINE: `~/Downloads/rigel_runs/arms/2026-09-19_e2e_baseline/` — one `qa_<arm>.jsonl` per arm, the
  per-stratum `qa_report.txt` (`quant_accuracy.py --report`), the floor-and-arm table `decomposition.txt`, the
  allocation arm `alloc.txt`, calibration's own split on the capture-OFF rows (`cvo*.json`), `tree.txt`, and
  `qa_keyed.py`, the wrapper the oracle arms ran through.
* The identity references: `~/Downloads/rigel_runs/arms/review_identity_*.json`, BIT-IDENTICAL on this tree.
  They are the gate for any change that must not move a number.
