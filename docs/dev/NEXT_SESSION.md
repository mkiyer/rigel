# NEXT SESSION — performance, after the first pass on the real libraries (2026-09-11)

A handoff, provisional like everything in this directory. The references are `CLAUDE.md`, `docs/ROADMAP.md`
and `docs/ISSUES.md` (`performance-memory-bounded-solve` carries the measured shape of the problem).

## The substrate

Four real libraries under `~/Downloads/rigel_runs/cfrna/<lib>/bam/star.srt.rmdup.collate.bam`, human index
`~/Downloads/rigel_runs/refs/rigel_index` (2.09M chain slots). Three cfRNA libraries are the smoke tests
(155k, 876k and 1.3M fragments); `mctp_vcap_rna20m_dna05m` (18.6M fragments, ~25% spiked DNA) is the stress
case. Libraries above 100M fragments come next.

## The tools

* `scripts/profiling/profiler.py` — the pipeline as a tree of named stages (50 probes, patched at run time),
  per-stage peak and held RSS, `--set`, `--scan-only`, `--compare`, `--self-test`.
* `scripts/profiling/sweep_replay.py` — capture every calibration sweep from one real run, replay one in
  isolation, compare bit for bit. Captures of the 876k library: `~/Downloads/rigel_runs/perf/sweeps_MO_3021`.
* `scripts/design/rename_identity.py --bam` — end-to-end identity on a real library; reference
  `~/Downloads/rigel_runs/arms/perf_identity_LBX0190.json`, plus the two ladder references
  `cleanup_identity_*.json`.
* Timing is only honest as back-to-back A/B pairs: stages nobody touched drift 25–40% between runs taken at
  different times (page cache, concurrent work). py-spy needs sudo on macOS.

## Where the time goes, and what landed

Back-to-back A/B at 8 threads, committed source against the working tree (all changes bit-identical: the
three identity references, and every captured sweep replayed):

| library | before | after | ratio |
|---|---|---|---|
| 876k fragments (two pairs) | 285.7 s, 285.0 s | 199.5 s, 201.1 s | 0.70 |
| 18.6M fragments | 1008.8 s | 835.7 s | 0.83 |

What moved: prior assembly 0.10 (a quadratic region-to-locus scan made linear), capture effective lengths
0.12 (the incidence vectorized; the KDE one grid row at a time, no 4 GB temporaries), mature walls 0.20 (one
shared `region_arrays.overlapping_region_runs`), the RNA lanes 0.61–0.75 and the directional passes 0.76
(face sets built vectorized, the per-slot kernel tightened, trigamma as `zeta(2, x)`, `blur_row`'s edge pad).

What is left on the deep library (835.7 s, 32.5 GB peak): the four sweeps 706 s — directional passes 320,
`prepare` 210 (RNA lanes 87, gDNA lane 50), ψ grid solves about 150 (inside own claims and the sweep) —
scan 32, second pass 22, locus EM 8. All of calibration runs on one core.

## Decisions for the owner

1. **The per-reference parallel sweep** (`ISSUES: performance-memory-bounded-solve`). The passes couple slots
   only within a reference and ψ is per slot, so a split by reference groups stays bit-identical if the three
   genome-wide scalars in `prepare` are computed once. The largest chromosome is about 9% of the slots, which
   caps the speed-up near 10×; running a bounded number of groups at a time bounds the memory.
2. **Threads for the ψ block loop.** Large NumPy array work that releases the GIL; bit-identical; about 150 s
   on the deep library. Needs a thread budget calibration does not have today.
3. **The scan's thread split.** `resolved_scan_threads` gives BGZF `min(4, total − 1)` threads, so a 2–4
   thread budget runs one scan worker. Measured on the deep library, scan seconds by (BGZF, workers):
   total 4 — (3,1) 113.5, (2,2) 59.5, (1,3) 42.4, (0,4) 33.7; total 8 — (4,4) 34.4, (2,6) 26.4, (1,7) 24.8;
   total 16 — (4,12) 19.1, (1,15) 19.5, (2,14) 17.4. Any new rule is a tunable, and the knob is a CLI flag.
4. **The refit count.** Each refit is a full sweep; on the 876k library they move 31k, 11k and 3.7k gDNA
   fragments of 502k, converging by a third per step. Cutting one changes answers.
5. **Caching index-derived geometry in the index** (mature walls, incidence, opportunity tables): now about
   3 s a run after the vectorization, so low priority.
