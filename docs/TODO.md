# TODO — the project's deferred work, ranked

**This is the one list.** Add here rather than starting another file, and **delete an item when it
lands** rather than marking it done — `LEDGER.md` records finished work with its gates and its reasoning.
Every item states **why it is deferred**, because "we'll get to it" is how a list stops being read.

Live handoff: `IMPLEMENTATION_PLAN.md` §0. Finished work: `LEDGER.md`. The suite: `BENCHMARK_SUITE.md`.

---

## ⭐ THE CRITICAL PATH

| | item | why now / why not yet |
|---|---|---|
| **1** | ⭐ **S5 — in progress: S5.e next** | ⛔ **Everything is still behind this** — `calibrate()` does not run, so the suite cannot produce a number and the goldens cannot be adjudicated. ⭐ S5.0/a/b/c/d have landed; the live plan is **`S5_DESIGN_LOG.md` §2**, which supersedes `IMPLEMENTATION_PLAN.md` §4/§5 |
| **2** | **S6 — delete** | `ruff check` undefined-name failures are the authoritative list; goldens regenerated once |
| **3** | Close the suite's two open requirements | (c) non-Poisson counts and (f) the low-gDNA corner. `suite_resolves.py` fails on both today and names them |
| **4** | The stress chromosome | The toy-scenario half of `testing_plan.md`. Needs S5 for the seed to be verifiable |
| **5** | A new benchmark skill | Only once the suite can produce a number. The old one was deleted; writing one now would be speculative |

---

## 1. S5 — in progress. **`S5_DESIGN_LOG.md` §2 is the live plan**

⭐ S5 is not a rewiring job. It was stopped on 2026-07-30 and turned into a derivation first, because the
design could not say which observables the rewiring should consume. `IMPLEMENTATION_PLAN.md` §4's R1–R4
ranking and §5's S5 row are **superseded**.

**Landed:** S5.0 (the derivation) · S5.a (`length_sum` on every population) · S5.b (`fl.py` → the five
pure pools) · S5.c (`effective_length.py` → the one placements formula) · S5.d (one substrate, the
`N E N … E N` chain). All in `LEDGER.md`.

**Next: S5.e** — `build_node_geometry` rewritten and `bp_solver`'s per-face consumers collapsed. It
clears the ~25 geometry failures, which currently fail with a message naming their step. Then **S5.f**
(`calibrate`, `CalibrationResult`, `priors`, `pipeline`), where calibration runs end to end and its
numbers become the FIRST entry of a new baseline.

⛔ **A7 must be ruled before S5.e finishes:** does reach enter the divisors, and how? An unspliced edge
has no single reach — the RNA part is bounded by its transcript, the gDNA part is not. `S5_DESIGN_LOG.md`
§1 A7.

## 2. The suite's two unmet requirements

`scripts/design/suite_resolves.py` fails on exactly these and prints why. Everything else resolves.

### (f) A low-gDNA × strong-capture corner — **the cheap one**

Real libraries live at 1–10 % gDNA; the pilot grid is `none`/`100 %` by the owner's 8-condition spec, so
nothing sits in that band. **Fix:** add one gDNA rate (e.g. `0.05`) to `chr22_pilot.yaml` and re-run those
conditions. `CARRY_FORWARD.md` §1 fact 15/16: capture destroys the intron signal 75×, and nascent-vs-gDNA
is unidentifiable under capture — so this corner is where the design's hardest failure mode lives.

### (c) Non-Poisson counts — **needs a mechanism decision first**

The simulator draws multinomial at fixed abundance (`wgs_engine._accumulate_pool`), so overdispersion is
**0 by construction** (measured ω < 5e-5). Deferred by owner ruling 2026-07-30.

⛔ **The mechanism first chosen — per-transcript Gamma with REALIZED truth — cannot work.** Overdispersion
is only visible as excess count variance relative to a predicted rate, and there are two places that rate
can come from: **truth**, which records the realized (post-Gamma) abundance, so counts are Poisson given
it; and **replicates**, which must share the Gamma draw or the conditions stop being comparable. Both
absorb it. The injected `od` would be a knob that measurably does nothing.

⭐ **Recommended instead: sub-transcript clumping** — a Gamma field along each transcript with its total
preserved. Per-transcript truth stays exact (so the mean comparison is fair) while **node** counts go
non-Poisson, and nodes are what calibration reads. ⚠ It also needs **replicate conditions** (same
parameters, different `sim_seed`) so ω is estimable at all.

## 3. The stress chromosome — the toy half of `testing_plan.md`

A designed synthetic reference appended to the suite backbone, carrying the cases a real chromosome does
not concentrate: a density step, alternative TSS/TES strictly inside exons, single-stranded neighbourhoods
flanking ambiguous nodes, 1 bp nodes, overlapping/abutting introns, and **strand-coincident junction
pairs** (`CARRY_FORWARD.md` §3 trap 24 — GENCODE has zero of them, so only a synthetic test can find it).

⭐ **The toy geometry that worked before is worth reusing** (from the deleted `toy_inject.py`): a full
3-exon transcript with **intergenic ends** — `intergenic — exon1(TSS) — intron — exon2 — intron —
exon3(TES) — intergenic`. The intergenic ends seed the baseline gDNA level that propagates through the
messages, which is what makes it a complete, self-grounded problem rather than a fragment.

⚠ **Blocked on S5 for the seed**, not for the reference: `scan_cache` step 4 needs `calibrate` to run.

## 4. The scan cache's step 4 — the population-prior seed

Steps 1–3 landed (`LEDGER.md` B3). Step 4 seeds a toy from a genome-scale scan.

⭐ **The mechanism already exists and needs wiring, not inventing.** `InjectedCalibrationPriors`
(`calibrate.py:82`) already carries exactly the population quantities a toy cannot fit — `rna_sense_frac`,
`n_rna_obs`, `n_gdna_obs`, both strand overdispersions, the enrichment NPMLE, the intron background, the
aggregate background reference — and `calibrate` already stashes the fitted bundle in
`_debug["calibration_priors"]`. The recipe is three lines: run `calibrate` on the cached population scan
with `_debug=`, pull the bundle, pass it as `injected_priors=` on the toy.

⚠ `testing_plan.md`'s constraint becomes a checkable assertion: the toy must run under the **same**
strand specificity, fragment-length distributions, gDNA level and capture model as the population scan.

**Why deferred:** extracting the priors requires `calibrate` to run. Pinned by a strict xfail in
`tests/test_scan_cache.py`.

## 5. The soft 3-pool surplus does not exist

`BENCHMARKING.md` names it as the **primary** pool metric — because the hard-label metric is nearly blind
to a calibration-prior change (a real change can move the soft pools by tens of thousands of fragments
while the hard-label net is byte-identical). `rigel.sim.analysis` implements only the hard-label version.

**Also missing:** the absolute per-transcript error alongside the net (`BENCHMARKING.md` caveat 2 — net
cancels). **Why deferred:** both need `rigel quant` to run end to end, i.e. S5.

---

## Smaller, self-contained, and each with its reason

### ⛔ The simulator hangs on an impossible fragment-length truncation

`sampling.truncated_normal_frag_lengths:30` rejection-samples in a `while` loop with **no termination
guard**. `Normal(206, 0)` truncated to `[200, 200]` never yields a valid draw and the simulator spins
forever instead of raising. Reachable from any config whose `frag_min`/`frag_max` exclude the mean —
including the obvious way to reproduce the deleted suite's `frag_std: 0`. It cost a 10-minute timeout on
2026-07-30. **Fix:** bound the loop and raise, naming the interval and the distribution.

### ⚠ Capture sampling is still O(probes per run) per draw

`LEDGER.md` B1b: after the batching work, `_run_landscape` is 47.8 s over 23,788 calls with `scatter` at
511,711 calls — the `(group slot × piece slot)` loop, ~21.5 iterations per call because at large fragment
lengths probe neighbourhoods merge into big runs. **The remaining lever is a rejection sampler**: propose
from the separable per-probe trapezoid (closed form, no per-start array), accept with `max_g / Σ_g` at the
one drawn position. ⚠ **Owner has authorised the RNG reordering this needs.** Deferred because the
measured state is already adequate (full pilot 741 s / 4.09 GB) and it is its own arm with a
distributional gate rather than a content-identity one.

### ⚠ Performance regressed at S4, and it is measured

| | pre-rework | S4 | re-recorded 2026-07-30 |
|---|---|---|---|
| per-fragment deposit | ~357 ns | 410.0 ns | **367.8 ns** |
| fixed partition cost | 0.108 s | 0.348 s | **0.315 s** |

The per-fragment figure sits inside the old 350–400 band and scatter explains it. ⛔ **The fixed
`O(partition)` cost does not** — it is still ~3× pre-rework. Consistent with `build_result`'s AoS→SoA
transpose and the per-worker zeroing, both `O(partition)`. Paid once per scan: amortised on a deep BAM,
dominant on a shallow one. **Why deferred:** owner ruling — *"i'm not worried about the performance gate.
we will make it fast eventually."* The non-negotiable part holds: the deposit allocates nothing per
fragment.

### ⛔ THE TWO "SPLASH" POOLS ARE NOT gDNA — interrogate `deposit()`'s exon/intron path

⭐ **Measured on LBX0190 (S5.b), and the arithmetic is damning.** With the pure pools giving gDNA 79.4 bp
and RNA 227.3 bp, a splash pool's mean implies its composition directly:

| pool | n | mean L | p50 | p90 | **implied gDNA share** |
|---|---|---|---|---|---|
| `DNA_INTERGENIC_EXON` | 123 | 231.1 | 210 | 380 | ⛔ **−2.6 %** — i.e. 100 % RNA |
| `DNA_INTRON_EXON` | 2,092 | 144.1 | 117 | 255 | ⚠ **56.3 %** — so ~44 % RNA |

`DNA_INTERGENIC_EXON` is statistically indistinguishable from the pure RNA pool (231.1 against 227.3).
**It is RNA bleeding across a transcript edge**, not on-target gDNA, and on that reading it should be
excluded rather than reported. `DNA_INTRON_EXON`'s p90 of 255 bp sits above the RNA *median* of 204 — a
long tail no gDNA population explains.

⚠ **Both are named as pure-by-construction in `ACCUMULATOR_DESIGN.md` §8 and neither is.** They are QC
only today (S5.b keeps them out of `gdna_fl_mass`), so nothing is currently mis-fitted — but the design's
claim is wrong as stated, and if they were ever folded in they would bias the gDNA model long by ~29 %.

**The owner's hypothesis, and it is the place to start.** An unspliced fragment can overlap **+1 or +2
bases into an intron** because the aligner could not place those bases on the next exon with so little
overhang. Those are *spliced RNA* molecules mis-aligned as unspliced, and they land in the intron flank —
exactly `DNA_INTRON_EXON`. Two things to establish:

1. **Does the aligner soft-clip those bases instead?** If it clips rather than extends, the fragment never
   reaches the intron and the mechanism is something else. Check the CIGARs directly.
2. ⭐ **The existing anchor machinery does NOT cover this case.** `splice_blacklist.py` +
   `max_anchor_left`/`max_anchor_right` (`constants.h:171`) reject *artifactual junctions* by requiring a
   minimum anchor on each side. That inspects **spliced** reads. The fragments hypothesised here carry
   **no junction at all** — that is the whole point — so the anchor check never sees them. Any fix is new
   logic, not a tuning of the existing gate.

**Other candidates to rule out before blaming the aligner:** an off-by-one in which flank types the
crossing pool is keyed by (`_pool` / `_SPLASH_POOL`, which sorts the flank pair — so an intron↔exon and an
exon↔intron seam are one bucket); the `sole_line` gate admitting multi-line crossings; and retained-intron
isoforms, which are genuinely RNA in an intron and would need the composition to be modelled rather than
excluded.

**Why deferred:** it is a debugging project of its own, with its own evidence gathering, and nothing
currently consumes these pools. Owner-agreed 2026-07-30. ⚠ It blocks nothing, but it does mean
`ACCUMULATOR_DESIGN.md` §8's purity claim for the two splash pools should be marked **unverified** until
it is settled.

### ⚠ `POOL_EB_PRIOR_ESS = 1000` now dominates the gDNA length model

⭐ Measured on LBX0190 after S5.b: the pure gDNA pool is **4,467** fragments, so a prior of 1000 pulls the
fitted mean **79.4 → 100.3** — a **+26 %** shift toward the global FL, i.e. toward RNA. The constant was
chosen when `gdna_counts` was the old four-pool aggregate, several times larger; the pure pools are
smaller *by design*, so the same number shrinks much harder.

⚠ It is load-bearing, not cosmetic: `μ_g − μ_r` is the determinant of every 2×2 in
`ACCUMULATOR_DESIGN.md` §6, and §7.2 measures a 10 % length-model error at 0.010–0.026 of composition.
Shrinking the gDNA mean 26 % toward RNA shrinks that determinant directly.

**Why deferred:** it is a magic number and changing it is the owner's call. Options are not obvious —
a smaller ESS, an ESS scaled to the pool, or shrinking toward the *splash* pools rather than the global
FL (they are on-target gDNA, so they are a better anchor than an RNA-dominated global). Needs a decision
before a number is picked.

### `strandedness` is not on the payload

Design §5.2 lists it as a header field — measure it, declare it, assert it in QC. Derivable from
`strand_observations`. ⚠ The new schema **cannot express** the bug it was meant to guard against — no
channel is labelled *sense*, so no label can be wrong — which is why this is a diagnostic to add rather
than a hole to plug.

### ⛔ `detect_chimera` is blind to two real populations

It considers only blocks with a non-empty transcript set (`constants.h:339-343`) and needs ≥2 mutually
disjoint transcript components, so it returns `CHIMERA_NONE` for **same-reference-orientation mates**
(now counted as `dropped_strand_undefined`: 158 / 615 / 246 on the three cfRNA libraries) and for the
**multi-megabase spans** (95 % carry a supplementary record; intergenic blocks have empty transcript
sets). Vocabulary exists: `CHIMERA_CIS_STRAND_SAME` / `CHIMERA_CIS_STRAND_DIFF`. **Why deferred:**
widening the gate reclassifies currently-accepted fragments, which is a change to *what counts as a
fragment*, not to how one is tallied — its own arm with its own before/after. Owner-agreed 2026-07-29.

### Single-reference *cis* chimeras still deposit on the intergenic path

`cr.chimera_type` is set for them and the resolved path honours it, but the intergenic path has no chimera
check. S3's adapter gates only on multi-reference, which is what a `FragmentPath` can express.
⭐ Measured: the narrow and broad gates admit **identical** fragments on LBX0190 and MO_3021, so nothing
is lost by deferring.

### The side buffer: multimappers and `path_ambiguous`

Design §9 has the recipe — group by `(candidate set, column)`, apportion each group's integer count across
candidates in proportion to their densities, rounding by **largest remainder**. ⚠ The assignment must be
**integral**: a fractional `1/NH` split re-creates the non-integer observable this redesign exists to
delete. Quantified: `path_ambiguous` is **70–76 %** of implicit splices (0.04–1.03 % of everything
accepted); multimappers are **59.3 / 39.5 / 20.2 / 16.7 %** of intergenic fragments and ~53.6 % of
intronic. ⚠ The gate thins the **one pure gDNA pool** non-randomly toward repeats, so this is a bias in
the length model, not only a recovered-fragment count.

### `contradictory_sj_strand` has never fired on real data

Measured 0 on all three cfRNA libraries — expected for STAR, which writes a consistent `XS` per record,
but it means one branch of the three-way `sj_strand` rule is untested outside the spec matrix. A
deliberately mixed-tag BAM would settle it.

### `sj.feather` duplicates junction coordinates

One row per *(junction, transcript)*, feeding `sj_map` at load (`index.py:1351`), so it cannot merge into
`edges.feather` (one row per distinct junction). But its `(ref, start, end, strand)` columns duplicate
what `(src, dst, strand)` already determines, and **nothing enforces they agree**. **Available
simplification:** re-key it to `(junction_edge_id, t_index)`. Owner's idea.

### `align_strand` → `strand` tree-wide

101 occurrences, and ⚠ **seven are string keys** (`buffer.py:193,342,381` — a parquet column in the spill
path — plus `resolve.cpp:38`, `resolve_context.h:307,352`, `tests/test_buffer.py:45`). A string key
survives compilation and fails at **runtime**, on the buffer→EM path. **Why deferred:** that path is
untouched by the accumulator rework, so bundling it in would put an unrelated runtime-failure mode inside
a step whose gate is byte-identity.

### The goldens moved at the index change and were never validated

`gdna_em_count` fell 16–52 %. The new suite should adjudicate it once S5 lands; the accumulator sequence
should not.
