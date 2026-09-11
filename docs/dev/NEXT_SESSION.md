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

## What is next — parallel calibration (owner, 2026-09-11)

Calibration is the tool's one unfinished component and the target is at least a hundredfold, C++ not yet
spent. **The decomposition is the LOCUS, exactly as the EM already does it.** An intergenic region
terminates message passing: it is solved, fixed and measured as pure gDNA. Verified on the real human
chain, not assumed — of 1,206,202 composition faces and 4,621,302 lane faces, NONE delivers to an
intergenic node, while 65,852 gDNA-lane faces are SENT by one (its measured level is the pure-gDNA anchor
its neighbour reads). So a locus needs only its two flanking intergenic claims, and those are local.

    chain slots 2,087,476   intergenic 33,120   loci 32,927
    locus size: median 19, p99 459, largest 2,477 = 0.12% of the chain

The largest locus is 0.12% of the work, so the serial floor is about 800x and the working set becomes per
locus instead of genome-wide, which is the memory half of the same problem. Threads are wanted for the
message passes and for the grid solves alike.

What must stay identical, and what that costs: the three genome-wide scalars inside `prepare` (the gDNA
lane's reference density, each RNA lane's, and whether the strand split is live) are sums over the whole
chain and have to be computed once, over the same arrays in the same order, then handed to the workers.
Everything else in a sweep is per slot or per face. `sweep_replay.py` replays a captured real sweep and
compares bit for bit in about 40 s, so each step is provable before the pipeline is run at all.

## The other decisions

* **The scan's thread split** — `ISSUES: scan-thread-split-starves-the-workers`: a 2-4 thread budget runs
  one scan worker. The numbers are in the entry; any new rule is a tunable, so the rule is the owner's.
* **The refit count** stays as it is (owner, 2026-09-11).
* **Caching index-derived geometry in the index** (mature walls, incidence, opportunity tables): about 3 s
  a run after the vectorisation, so low priority.
