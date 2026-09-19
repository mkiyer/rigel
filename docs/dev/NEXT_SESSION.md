# NEXT SESSION — start here (2026-09-19: 0.8.0 is a RELEASE OF THE TOOL; the first job is the end-to-end baseline)

This file is only how to begin. The ranked view is `docs/ROADMAP.md`, the open problems are `docs/ISSUES.md`,
the rulings and the record are `docs/DESIGN.md` (§0b carries the scope and its 2026-09-19 amendment), what
"done" means is `docs/SUCCESS.md`, the lessons are `docs/TRAPS.md` cited by name, and the release procedure is
`docs/MANUAL.md` / `docs/PUBLISHING.md`.

## Where the tool is

Everything through `aaa0d8a2` is landed and pushed; the tree is clean and `main == origin/main`. Calibration is
measured against an oracle per stratum and is no longer the thing in the way; two machine campaigns made the
sweep one native call and took a deep run from 145 s to 125 s, and **not one of those commits moved a number**,
so every accuracy measurement recorded before them still stands. The machine thread is parked and resumable
(`ISSUES: performance-memory-bounded-solve`).

What changed on 2026-09-19 is the FRAME: the release is the whole tool, so the transcript table is a
first-class number beside the calibration metric. The two answer different questions and neither stands in for
the other (`DESIGN.md` §0b, `SUCCESS.md`).

## The first job: measure the deliverable (`ISSUES: end-to-end-error-unattributed`)

⛔ A MEASUREMENT SESSION. No `src/` change belongs in it, because the point is a baseline nothing has been
tuned against, and because the order of everything after it depends on what it says.

1. `python scripts/design/preflight.py --full`, then `python -m pytest tests/ -q` against the standing count
   in `CLAUDE.md`. ANY failure is a regression.
2. `quant_accuracy.py --arm base` on the ladder, and `--arm base_reseed` IN THE SAME SESSION. The floor is
   not optional: the deliverable is not reproducible by default, and a delta below the floor is sampling
   noise, not a result (`TRAPS: re-record-the-baseline`).
3. The decomposition, same session, same conditions: `oracle` (all three prior fields at truth),
   `oracle_gdna`, `oracle_rna`, `oracle_efflen`, and `oracle_ruler` — the only arm that reaches the
   effective-length shrinkage, because it substitutes at the `calibrate` boundary while every other arm wraps
   `assemble_priors`.
4. Read it PER STRATUM, never pooled (`TRAPS: never-pool-the-strata`), with the deferred stratum reported and
   not ranked. What remains under `oracle` belongs to the EM and the assignment by construction; what
   `oracle` recovers belongs to the prior chain.
5. Write the numbers into `ROADMAP.md`'s state claim for the deliverable, close or re-rank
   `ISSUES: end-to-end-error-unattributed`, and only then choose between roadmap items 2 and 3.

## Then: the pre-EM setup (`ROADMAP.md` item 2), in this order

The bridge between a good calibration and the user's number. Ranked by what is already measured:

1. **`ISSUES: per-transcript-prior-lane`** — `rna_prior_weight` is built end to end and `pipeline.py` omits
   it, so the shipped EM carries NO per-transcript information; a perfect version of it roughly halves
   in-scope gene-level error. It is a wiring gap plus a support decision, and two weightings are already
   refused with their numbers.
2. **`ISSUES: capture-blind-gdna-divisor`** — +6.0 % on all six capture-ON rows.
3. **`ISSUES: eb-shrinkage-magic-ess`** — a magic ESS, inert on the ladder and dominant on the fl-gap arm.
4. **`ISSUES: prior-fidelity-vs-deliverable`** — the anti-correlation; the baseline's arms may answer it.
5. **`ISSUES: antisense-prior-assembly-casualty`** — the assembler's alpha = 0 rule, which owns an xfail.

Each judged on BOTH numbers: the deliverable above its floor, and `prior_vs_oracle.py` for the assembler's own
error, so a repair that moves the prior and not the user's number is visible as exactly that.

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
* The identity references: `~/Downloads/rigel_runs/arms/review_identity_*.json`, BIT-IDENTICAL on this tree.
  They are the gate for any change that must not move a number.
