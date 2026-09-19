# NEXT SESSION — start here (2026-09-19: the machine is fast enough; THE METHOD is the focus again)

This file is only how to begin. The ranked view is `docs/ROADMAP.md`, the open problems are `docs/ISSUES.md`,
the rulings and the record are `docs/DESIGN.md`, the lessons are `docs/TRAPS.md` cited by name, and how
performance is judged is `docs/SUCCESS.md`.

## Where the tool is

Everything through `b834df80` is landed and pushed; the tree is clean and `main == origin/main`. Two campaigns
finished back to back and are recorded in `DESIGN.md` §6b.15: THE PORT, which made the calibration sweep one
native call over a pool of threads, and THE WORK OUTSIDE IT, six phases that took a deep run from 145 s to
125 s. ⭐ **Not one of those commits moved a number.** Every accuracy measurement recorded before them still
stands, and the three frozen identity references are bit-identical on the current tree.

So the machine is no longer the thing in the way, and the owner's order of 2026-09-19 is THE METHOD: what the
tool answers, judged on 0.8.0's metric per stratum, end to end.

## Start here, in this order

1. `python scripts/design/preflight.py` — one command, one verdict, before anything else. `--full` adds every
   instrument's self-test in seconds.
2. `python -m pytest tests/ -q` — the standing baseline is in `CLAUDE.md` and ANY failure is a regression.
3. Re-read `docs/ROADMAP.md`'s ranked list. It is the method's list now: the strand tilt where AMBIG slots
   carry RNA on both strands, the pre-EM setup (`priors` / `result` / `derive`), the intron's own solve on
   unstranded capture-OFF, the vertex atom, and a message mechanism only where a row is above the bar.
4. Before proposing anything, take the CURRENT picture with the instruments that own it, in this order:
   `calibration_vs_oracle.py` (0.8.0's metric, per stratum, ~5–12 s a condition), `solvability_audit.py`
   (which objects are solvable, solved wrong, or confidently wrong), `policy_benchmark.py --panel ladder
   --by-class` (where a policy's remaining error sits, halves NEVER pooled). ⛔ The numbers in the docs are
   from before the two campaigns; they are still valid, because nothing moved, but re-record a baseline in
   the same session before quoting a delta (`TRAPS: re-record-the-baseline`).
5. The method loop is the debug loop: the panel → the worst IN-SCOPE scenario → its highest-error objects
   (`worst_objects.py`, `calibration_walk.py`) → the mechanism → a gated fix → the panel again.

## The five xfails are the method's own worklist

They are executable records of proven defects, each deferred to its thread by ruling, and closing one means
repairing the thing or asserting the invariant structurally — never widening a bound:
`ISSUES: two-sided-exon-row`; `ISSUES: antisense-prior-assembly-casualty`;
`ISSUES: the-lower-bound-noise-ratchet`; `ISSUES: nested-antisense-leak-under-the-sane-ruler` (two rungs).

## The other kind of work, parked and resumable

`ISSUES: performance-memory-bounded-solve` is the MACHINE thread: judged on seconds, bytes and bit-identity
rather than on the metric, and parked by the owner rather than finished. It carries what a deep run costs now,
what is ranked next with its measured price, and two candidates already researched and refused as targets (the
scan, which is at its floor for an eight-thread budget, and the index load, which is already mostly native).
The next item there needs no re-derivation: the short-template taper table has two rewrites derived and
measured, one bit-identical and one priced.

⛔ The protocol, when it resumes: a capture replays the sweep (`perf/sweeps_VCaP_step19`, 2.7 GB, bit-identical
on this tree), but for anything OUTSIDE the sweep the gate is `rename_identity.py --check` on the three frozen
references — and only the `--bam` one runs the scan and the second pass. A timing is read from two interleaved
pairs at 8 threads on VCaP against a worktree of the pre-step commit carrying its own modules
(`s20/pre_worktree_p4.sh` and `s20/time_pairs_p4.sh` are the current forms).

## Standing rulings carried (unchanged)

* Float64 for the whole of ψ; ONE solver; ONE λ lattice; no θ lattice; the tilt's hypothesis space is
  {pure +, pure −, mixed}. The strand channel's liveness is a protocol decision on the spliced 2×2.
* The landscape trains only where a solve locates a slot; the ruler reads its located enriched mode and admits
  failure where there is no gDNA to read (owner, 2026-09-15).
* Real data is a test input, never a design input. A few high-quality instruments, kept current. The source
  cites no doc. One production path: once native code is validated the Python it replaces is deleted.
* Unstranded × capture-ON stays DEFERRED — reported on every benchmark, never a development target, never
  ranked on a pooled total.
* The owner drives commits. The refit count and CI's on-demand trigger are the owner's.

## Where everything is

* The synced scratchpad: `~/Downloads/rigel_runs/prototypes/2026-09-17_port_ii/` — `commits/committed/` (every
  snapshot 7–25, all landed), `s16/`–`s18/` (the port and the block), `s19/` (the first attribution),
  `s20/` (the performance campaign: the scan-split and arena measurements, the worktree and timing scripts,
  every phase's identity log), `s21/` (the attribution of the landed tree and the taper study).
* Captures and reports: `perf/sweeps_VCaP_step19`; `perf/plan_final_2026-09-19/` is the campaign's own
  before-and-after, with `perf/phase3_2026-09-19/` and `perf/phase4_2026-09-19/` beside it.
* The identity references: `~/Downloads/rigel_runs/arms/review_identity_*.json`, frozen 2026-09-18 on the
  log-gamma tree's numbers, BIT-IDENTICAL on the current tree.
