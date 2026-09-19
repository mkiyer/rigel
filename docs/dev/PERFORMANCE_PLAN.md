# The performance plan after the block went native (2026-09-18) — six phases, ranked

The block in one native call landed and the balance changed: calibration is now 28 % of a deep run and the
four sweeps inside it 17 %. Everything else — the scan, the second pass, the two fragment-length fits,
quant, the index load — is the other 72 %, and almost all of the reducible part of it is PYTHON doing
per-object work that numpy or the C++ already does per array. This file is the ranked plan to take it out,
phase by phase. It is a working note: when a phase lands, its record goes to `DESIGN.md` §6b.15.5 and
`ISSUES: performance-memory-bounded-solve`, and this file keeps only the plan.

## The baseline this plan is ranked against

VCaP (`mctp_vcap_rna20m_dna05m`, 18,568,456 fragments) at `--threads 8`, the tree at `e921869c`, two
interleaved pairs taken 2026-09-18 (`perf/cache_deleted_2026-09-18/pair{1,2}_post.json`): the run 148.4 s
and 145.1 s, peak RSS 10,632 MB and 10,326 MB. The stage table below is `pair2_post`, the faster of the two;
read every saving against the STAGE row, never against the wall, because the wall's own drift between runs
is 25–40 % of a few seconds.

| stage | s | what it is |
|---|---|---|
| calibrate | 40.6 | the four sweeps 24.7, the landscape fits 4.2, ψ before the sweep 0.6, its own Python 10.6 |
| quant | 34.2 | the locus EM 14.3, the capture effective lengths 6.0, scoring 5.9, the priors 3.7, the partition 3.0 |
| the scan | 33.0 | one stage, C++, `BamScanner::scan` |
| the second pass | 21.6 | scoring the held fragments 11.8, a fragment-length fit 8.3, the drain 1.3 |
| a second fragment-length fit | 8.0 | on the DRAINED tally, before calibration, by design |
| the index load | 6.5 | |

**The attribution instrument.** One `profiler.py --cprofile` run on the same library
(`s19/vcap_cprofile.prof`, 2026-09-18). It inflates the wall to 172 s, so its numbers are SHARES, never
timings, and the phases below quote them as attribution only. What it found, by call count:

| call | count per run | where |
|---|---|---|
| `np.searchsorted`, scalar | 2,708,183 | the second pass's `_exact_region_bound` (1,853,986) and `_sj_id` (926,993) |
| `np.mean`, on a list of ≤ 2 | 1,429,451 | `fl._realized_gdna_counts`, once per exon per fit |
| `poisson_lower_mean` | 808 | `gdna_density.one_sided_rate`, 4 fits × 202 bisection steps |
| `_cum_short` | 1,982 | `effective_length.interval_sums`, one table per distinct short length |
| `_region_locus_shares` | 2 | `priors.assemble_priors`, the same triples computed twice |

**What the gates are, for work OUTSIDE the sweep.** The sweep replay
(`sweep_replay.py`, `perf/sweeps_VCaP_step19`) gates the sweep and nothing else, so it does NOT cover any
phase here except phase 5. The gate that does is `rename_identity.py --check` on the three frozen
references, which compares the quant digest and every array's content end to end: two panel conditions,
which run from a SCAN CACHE, and `review_identity_LBX0190.json`, which runs `--bam` from the real library
and is therefore the only one that exercises the scan and the second pass. Every phase below names which
of the three can see it.

## The phases

Ranked by measured value against implementation risk. Each numbered item is ONE commit, prepared as a
snapshot for the owner's go. The protocol is unchanged: `preflight.py` first, a falsification test before
the fix, the suite with its count re-derived, perturb the fixed code and watch each gate fire, the three
identity references, and two INTERLEAVED pairs at 8 threads on VCaP against a worktree of the pre-step
commit (`s18/pre_worktree.sh`, `s18/time_pairs.sh`).

### Phase 0 — the record (no code)

The roadmap's performance claim says the sweeps dominate a deep run, which stopped being true today; its
ranked item 1 has grown into a paragraph of landing history, which the file's own no-history rule forbids;
and the preamble to its ranked list names a session order agreed 2026-09-14 that is finished.
`ISSUES: performance-memory-bounded-solve` carries a 264 s baseline from 2026-09-17 and a ranked
opportunity list whose first four items landed 2026-09-18.

1. **Re-rank the record.** `ROADMAP.md`: the performance claim becomes "calibration is a quarter of a deep
   run and the stages outside it are the rest", item 1 loses its history and names this plan's phases,
   the stale preamble goes. `ISSUES: performance-memory-bounded-solve`: the baseline table above replaces
   the old one, and the ranked opportunities become the phases below. No number is claimed that is not in
   this file. Suite: content-only, so the collected total moves by +1 for this file alone.

### Phase 1 — the exact micro-wins (≈ 4 s, three tiny commits, all bit-identical by construction)

The cheapest work in the plan, and the safest: each is a few lines, each provably cannot move a number.

1. **The bisection stops when the bracket is closed** (`gdna_density.one_sided_rate`, ≈ 1.5 s).
   `_BISECT_STEPS = 200` halves a float64 bracket 200 times. After about 60 halvings `lo` and `hi` are
   adjacent floats and `mid = ½(lo + hi)` rounds to one of them, so every later step reassigns the same
   value: steps 61–200 are arithmetically dead. Replace the fixed count with the convergence the bracket
   itself states, `hi - lo > eps * hi`, keeping a hard cap. The magic number dies with it, which is the
   real point: the rule becomes derivable. BIT-IDENTICAL, and the gate is that fact — a test that runs the
   old loop to 200 steps and the new one to convergence and asserts `==` on the returned rate over a grid
   of brackets, plus the three references.
2. **The Poisson identity's log-gamma is a table** (`gdna_density.poisson_lower_mean`, ≈ 1 s). `k =
   floor(lam)` is an INTEGER, so `gammaln(k + 1)` is `log(k!)` at an integer, evaluated per element over
   the whole support at every bisection step. Build the table once per call, `gammaln(arange(k_max + 2))`,
   and index it; above the table fall back to `gammaln`. Bit-identical because the table's entries are
   that same function's values at those same points. Gate: an equality test against the direct form on a
   dense `lam` grid including the boundaries, plus the references.
3. **The region↔locus overlap is computed once, as its docstring says** (`priors`, ≈ 1.7 s).
   `_region_locus_shares` says "THE region↔locus overlap, computed exactly once" and is then computed
   twice per run: once in `assemble_priors` and once inside `_boundary_locus_shares`, on the same three
   arguments. Pass the triples in. Bit-identical: the same function, the same arguments, one call fewer.
   Gate: the prior gates, a spy asserting one call per `assemble_priors`, the references.

### Phase 2 — the scan's thread split (≈ 8 s, one commit, one owner decision)

`BamScanConfig.resolved_scan_threads` gives BGZF decompression `min(4, total − 1)` threads and the scan
workers the remainder. Measured on this library 2026-09-11 (`ISSUES:
scan-thread-split-starves-the-workers`), scan seconds by (bgzf, workers): total 4 — (0,4) 33.7, (1,3) 42.4;
total 8 — (4,4) 34.4, (2,6) 26.4, (1,7) 24.8; total 16 — (4,12) 19.1, (2,14) 17.4. The shipped rule picks
the WORST cell at every budget it was measured at.

1. **Derive the split from the measured ratio.** Re-measure the table on the current tree with
   `profiler.py --scan-only` (a scan is 33 s, so a full sweep of five splits is minutes, not an hour), then
   set the rule from what it says. `bgzf = total // 8` reproduces the best measured cell at all three
   budgets, and its justification is a measurement, not a taste: one decompression thread feeds about
   seven workers on this file. ⛔ THE OWNER RULES: it is a tunable and `--scan-bgzf-threads` is a
   user-facing flag.
   ⚠ **The risk that decides this phase**: the split changes how many WORKERS deposit, and a deposit is a
   float accumulation with a largest-remainder apportionment. If the merge is not order-independent the
   scan's output moves. `review_identity_LBX0190.json` is the gate that sees it, and it runs the scan; if
   it DIFFERS, this phase is not a no-op and becomes a priced decision rather than a free win.

### Phase 3 — the second pass's boundary lookups (≈ 7 s, two commits)

`score_held_fragments` loops over every held fragment in Python and, per hypothesis and per intron, asks
two questions of the region-bound axis: "is this position a region bound" (`_exact_region_bound`) and
"is this intron an annotated sj" (`_sj_id`). Each asks numpy for ONE scalar `searchsorted`, so the run
makes 2.7 M of them, and the cProfile puts 15.7 s of its inflated time there. Both helpers say in their own
docstrings that they mirror C++ that already exists: `Accumulator::sj_edge_id` and
`Accumulator::exact_region_bound` (`native/calibration/accumulator.h`). The second pass already HOLDS that
accumulator per reference, for `length_under`. So this is the one-path ruling and the speed-up in one
change.

1. **Bind the two lookups, batched.** Add `sj_edge_ids(starts, ends, strands) -> ids` and
   `exact_region_bounds(positions) -> ids` to the bound `Accumulator`: arrays in, arrays out, the loop
   inside C++ where it costs nothing. One call per reference replaces a million calls per run.
   ⛔ First prove the ID SPACES agree: the Python returns an index into the payload's sj axis, the C++ into
   the table `set_sj` was given. A one-off assertion pass over a whole library, both forms on every
   fragment, is the falsification test, and it runs BEFORE anything is deleted.
2. **Restructure the scoring loop around a pre-pass, and delete the Python mirrors.** The deferred payload
   already holds the hypotheses' introns as a CSR, so every lookup the loop will need can be gathered and
   answered in one batched call per reference before the loop starts; the loop then indexes an array.
   `_sj_id` and `_exact_region_bound` are deleted, and `_distinguishing_boundaries` takes its two bounds
   from the pre-pass. Gate: the second-pass and drain gates unchanged, the LBX0190 reference, the suite.
   Expected BIT-IDENTICAL — no arithmetic changes, only where the search runs.

### Phase 4 — the fragment-length fits' two index-scale loops (≈ 11 s of 16.3, two or three commits)

`build_fl_models` runs twice per pipeline by design, once on pass one's payload and once on the drained
tally, and costs 8.3 s and 8.0 s. Inside, `_realized_gdna_counts` holds two Python loops over
index-scale objects.

1. **The adjacent-pair loop, vectorised** (≈ 4 s). It walks every adjacent region pair of every reference
   twice (the `for _ in range(2)` refresh), keeps the exon|non-exon pairs, and accumulates two
   class sums plus a per-exon list. All of it is array work: the pair table is `ref_region_offsets`
   differenced, the class masks are comparisons, the sums are `np.bincount`. BIT-EXACT BY CONSTRUCTION,
   and the reason is worth stating because it is what makes this phase safe: `bincount` accumulates in
   input order, so generating the pairs in the loop's own order reproduces its summation order exactly.
2. **The per-exon filter, vectorised** (≈ 5 s). The second loop visits EVERY exon region to compute a mean
   enrichment, a weight and an excess, then does array work only for the exons whose excess is positive.
   The mean is over AT MOST TWO entries, because an exon has at most two flanking boundaries, so a
   grouped `bincount` mean is bit-identical to `np.mean` on the list (pairwise summation and sequential
   summation agree below eight elements). Vectorise the filter, keep the surviving loop as it is.
   ⛔ MEASURE THE SURVIVOR COUNT FIRST: if few exons survive on a real library the phase ends here, and if
   many survive the remaining loop is item 3.
3. **The surviving accumulation, only if item 2 leaves it hot** (≈ 2 s, and the only item in the plan that
   may move a number). The accumulated shape depends on the exon only through its LENGTH, so the sum over
   exons is a sum over distinct lengths weighted by a cumulative count — exact in exact arithmetic,
   different in summation ORDER in float64. That makes it a tolerance commit like the log-gamma one: the
   replay budget does not cover this stage, so the price is read from the three references and from
   `calibration_vs_oracle.py` on both arms. Do it alone, or not at all.

### Phase 5 — memory (one commit, plus one measurement)

The run's peak is 10.3 GB and it sits in QUANT, not in calibration: the scan's fragment buffer holds
3.6 GB alive until quant reads it, by design, the index holds 1.4 GB, and quant's scored candidates add
2.2 GB. Calibrate's own peak is 7.9 GB.

1. **The sweep's arena, priced and sized.** Each `solve_blocks` call allocates, per thread, sixteen
   `(block slots × grid)` float64 arrays; at the shipped 5,000 slots, the refits' grid of 233 and eight
   threads that is 1.19 GB, allocated and first-touched and freed four times a run. Two levers, and they
   are independent: `CalibrationConfig.sweep_block_slots` scales it LINEARLY and moves no number (the
   block size is gated to the bit), and keeping the arenas alive across calls trades that peak for the
   page faults. Measured already: `block_slots = 1000` replays the first sweep at 3.01 s against 2.82 s,
   so a fivefold memory cut costs about 7 % of one sweep. ⛔ THE OWNER RULES the default, because it is a
   memory-for-speed trade with no right answer.
2. **Measure quant's 2.2 GB before proposing anything.** The scored-candidate arrays are the EM's input,
   so cutting them is a streaming redesign, not a tuning knob. No commit is planned here until the
   measurement says what the 2.2 GB actually is.

### Phase 6 — the re-baseline and the record

1. **Two interleaved pairs on the finished tree**, against a worktree of `e921869c`, plus one
   `--cprofile` run to re-attribute what is left. The record goes to `DESIGN.md` §6b.15.5 and
   `ISSUES: performance-memory-bounded-solve`; `CLAUDE.md`'s baseline line is re-derived; this file is
   deleted or reduced to the phases that did not happen.

## What is NOT in this plan, and why

* **The sweeps' kernels (24.7 s).** The blur's loop interchange is priced at −1.5 s on a first sweep for a
  summation-order change at 1.9e-6 of the replay's budget, and it remains the owner's call. Everything
  else inside the block is native, threaded and bit-identical; there is no cheap arithmetic left.
* **The locus EM (14.3 s) and the capture effective lengths (6.0 s).** Both already native. `_cum_short`
  inside the effective lengths shows 3.8 s of self time on 1,982 table builds and is the one candidate
  there, but the table is built once per distinct length by design, so making it cheaper is a derivation,
  not a clean-up. It is the first thing to look at AFTER these phases.
* **The index load (6.5 s) and the scan's C++ itself (33 s minus the split).** Unexamined. The split is
  phase 2; the scanner's own throughput is a separate study.
* **Anything that trades accuracy for speed.** No phase here may move a calibration number except
  phase 4 item 3, which is priced and isolated for exactly that reason.

## The expected arithmetic

| phase | expected | confidence |
|---|---|---|
| 1, the exact micro-wins | −4 s | attributed by cProfile, bit-identical by construction |
| 2, the scan split | −8 s | MEASURED on this library, at the risk named above |
| 3, the second-pass lookups | −7 s | attributed, no arithmetic changes |
| 4, the length fits | −11 s | attributed; items 1 and 2 bit-exact, item 3 priced |
| 5, memory | −1.2 GB of calibrate's peak | measured from the arena's own arithmetic |

A run of 145 s becomes about 115 s if every phase lands, with the peak unchanged and calibrate's peak
down by a fifth. No kernel is touched and no number moves, except where phase 4 item 3 says otherwise.
