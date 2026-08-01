# TODO — the project's deferred work, ranked

**This is the one list.** Add here rather than starting another file, and **delete an item when it
lands** rather than marking it done — `LEDGER.md` records finished work with its gates and its reasoning.
Every item states **why it is deferred**, because "we'll get to it" is how a list stops being read.

Live handoff: `IMPLEMENTATION_PLAN.md` §0. Finished work: `LEDGER.md`. The suite: `BENCHMARK_SUITE.md`.

---

## ⭐ THE CRITICAL PATH

⛔ **READ `docs/FRAGMENT_LENGTH_AUDIT.md` §0 FIRST — the critical path moved there on 2026-07-31.**
S5 is finished. The live sequence is C2 → C3 → the P2 A/B → regenerate the goldens once. The ranks below
are everything *else*, and rank 1 has changed: the goldens must be regenerated **after C3**, not now.

⭐ The gate that blocked everything is passed. `calibrate()` runs (S5.f, 2026-07-30) and the FIRST
BASELINE exists in `LEDGER.md`. Items 3–5 below were all "needs S5" and are now **unblocked**.

⚠ **The rank is the priority; the § is where the detail lives.** They are not the same numbering.

| rank | item | § | why now / why not yet |
|---|---|---|---|
| **1** | ⚠ **Regenerate the goldens — but AFTER C3, not now** | — | ⛔ **Changed 2026-07-31.** 21 goldens are stale because P1 changed the EM prior's units, and C2/C3 will move them **again** (the FL models change ⇒ scoring and eff-lengths change). Regenerating now would bake in a number that is about to change. ⭐ Still true when the moment comes: **regenerate twice and diff** (§7) — the goldens run under the default sampling mode, so a flaky expectation baked in now is permanent. *(Original note: ⚠ the suite reads 22 failed and 21 of them are expected* — that is exactly the state a real regression hides in, and everything after this is an A/B needing "did anything else move?" to be a strong signal.)* |
| ~~2~~ | ⛔ **S5.g — A7's taper: MEASURED AND REFUTED** | §1 | The A/B is done (`LEDGER.md` S5.g-2): **≤ 0.0002** on the library gDNA fraction. The 11.0 % was bp-weighted geometry; mass-weighted the taper is 0.9596 and **89 % of edge mass is on lines it never touches**. ⚠ Would not hold for a 3′-biased/degraded library — screening test is the mass-weighted taper ratio |
| **3** | ⚠ **`EMConfig.seed` does not make `assignment_mode="sample"` reproducible** | §7 | ⭐ **No longer a blocker**: the default draws from the posterior *by design*, and running an A/B under `map`/`fractional` drops the noise 8 orders of magnitude. What remains is that a **seeded** run still does not reproduce, which is what `seed` is for |
| **4** | Close the suite's two open requirements | §2 | (c) non-Poisson counts and (f) the low-gDNA corner. `suite_resolves.py` fails on both today and names them. ⭐ Now unblocked |
| **5** | The stress chromosome + the scan cache's prior seed | §3, §4 | The toy-scenario half of `testing_plan.md`. ⭐ Now unblocked — the seed needs wiring, not inventing |
| **6** | ⚠ `rna_sense_frac` is mirrored | §6 | A **labelling** defect on a correct inference (proven: forcing the nominal truth is 166× worse). Nothing downstream is wrong today, but the reported κ says the opposite of what it means |
| **7** | A new benchmark skill | — | ⭐ Now unblocked: the suite can produce a number. Wants ranks 3 and 4 first, or it will be a skill that reports noise |

---

## 1. S5 — **`S5_DESIGN_LOG.md` §2 is the live plan**

⭐ S5 is not a rewiring job. It was stopped on 2026-07-30 and turned into a derivation first, because the
design could not say which observables the rewiring should consume. `IMPLEMENTATION_PLAN.md` §4's R1–R4
ranking and §5's S5 row are **superseded**.

**Landed:** S5.0 (the derivation) · S5.a (`length_sum`) · S5.b (`fl.py` → the five pure pools) · S5.c
(`effective_length.py` → the one placements formula) · S5.d (one substrate, the `N E N … E N` chain) ·
S5.e (the faces dissolve; A7 ruled) · ⭐ **S5.f — `calibrate()` runs, and the FIRST BASELINE exists.**
All in `LEDGER.md`.

**Next: S5.g — A7 proper.** The contiguous-edge RNA reach is deliberately `UNBOUNDED_REACH` today, so
the first baseline carries a **measured 11.0 %** genome-wide gDNA over-call. Turning the taper on is the
first change that gets a real A/B, which is the whole reason A7 was deferred past S5.f.

✅ **S5.g-1 landed** (`LEDGER.md`): `build_contiguous_edge_reach_arrays` gives the per-edge, per-strand
reach, gated by 6 tests and 6 perturbations, and it independently reproduces fact 6's magnitude (tapered
/ unbounded = 0.84–0.86 against the recorded 0.8904).

⛔ **S5.g-2 IS BLOCKED ON ONE OWNER DECISION.** Turning the taper on makes `NodeGeometry.eff_rna`
`[n_slots, 2]`. Two consumers take that without a decision — `node_init._rna` is already called once per
strand, and `node_total_density` already sums strands for junctions. **`bp_solver.py:334` does not**: its
`E_r` feeds the transfer variance `Var(log ρ_tot) = 1/n + [(1/E_g − 1/E_r)/B]²·Var(f_g)`, a
strand-agnostic scalar. The three honest candidates for reducing two divisors to one:

| candidate | for | against |
|---|---|---|
| belief-weighted `M·(f_pos+f_neg)/ρ_r_total` | the only one reproducing the true total RNA density | belief-dependent, evaluated at the INPUT belief — trap 11's family |
| max over strands | at the 40 % of lines with one live strand it IS that strand's reach | ⛔ **no monotone safety argument**: the term is `(1/E_g − 1/E_r)²`, which is ZERO at `E_r = E_g` and grows both sides. On the suite's pools (`E_g` 156, `E_r` unbounded 205) the zero sits INSIDE the tapered range, and at a deep terminus (`E_r` 50) the damping is **79×** the unbounded value |
| leave it `UNBOUNDED_REACH` | minimal change; it is damping, not a density | ⛔ two meanings on one name — trap 27, the defect S5.e removed |

⚠ **This is a new heuristic and the standing rule is to stop and discuss before adding one.** Nothing is
guessed; the arrays are in place so the decision is one edit away.

⚠ **A/B on the CALIBRATION numbers, or end to end under a DETERMINISTIC assignment mode** (§7).
Calibration is bit-identical run to run; the EM's default `assignment_mode="sample"` draws from the
posterior and moves ~0.5 %, while `map`/`fractional` move ~1e-10.

**Then, in order (`S5_DESIGN_LOG.md` §2 Phase 2):** S5.a2 (how `length_sum` enters the solve — it is
stored on every population and consumed by nobody) · A6 then A3 (`node_spanning`: the largest single win
found, 0.000 → 0.758 efficiency at a 25 bp node, and **56.7 %** of human nodes are shorter than one
fragment — but A6 must settle the spanning⊂crossing overlap first) · does `spliced_count` enter the level.

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

⭐ **UNBLOCKED (S5.f)**: `scan_cache` step 4 needed `calibrate` to run, and it does. The reference half
was never blocked.

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

⭐ **UNBLOCKED (S5.f)**: extracting the priors required `calibrate` to run, and it does. ⚠ The strict
xfail in `tests/test_scan_cache.py` should now XPASS once wired — check it rather than assuming, since a
strict xfail that silently keeps failing is indistinguishable from one nobody touched.

## 5. The soft 3-pool surplus does not exist

`BENCHMARKING.md` names it as the **primary** pool metric — because the hard-label metric is nearly blind
to a calibration-prior change (a real change can move the soft pools by tens of thousands of fragments
while the hard-label net is byte-identical). `rigel.sim.analysis` implements only the hard-label version.

**Also missing:** the absolute per-transcript error alongside the net (`BENCHMARKING.md` caveat 2 — net
cancels). ⭐ **UNBLOCKED (S5.f)**: both needed `rigel quant` to run end to end, and it does. ⚠ Build the
metric against item 7's noise floor — a soft-pool difference smaller than the EM's own run-to-run spread
is not a result.

## 6. ⚠ `rna_sense_frac` REPORTS THE MIRROR OF WHAT ITS NAME SAYS

⭐ **Measured against ground truth for the first time (S5.f).** A chr22 pilot library simulated at
`strand_specificity = 0.99` calibrates to **κ = 0.0101**; at 0.50 it reads 0.4990–0.5002, where a mirror
is invisible by construction. This is `CARRY_FORWARD.md` §0 **C4**'s open question — "sense fraction is
0.002–0.012 on all four real cfRNA libraries (nearly fully antisense) — possibly a read-orientation
convention bug" — and the answer is that it **is** a convention flip, not biology.

⭐ **BUT THE FLIP IS CONSISTENT, SO THE INFERENCE IS CORRECT.** This was the thing worth establishing,
and it is a measurement rather than an argument. On a **zero-gDNA** stranded condition (truth
`f_gdna = 0` exactly), injecting κ via `InjectedCalibrationPriors`:

| κ used | `f_gdna` on `gdna_none_ss_0.99_capture_off` | on `..._capture_on` |
|---|---|---|
| **0.0100 — the FITTED value** | **0.0030** | **0.0016** |
| 0.99 — the simulated "truth" | ⛔ **0.4992** | ⛔ **0.4822** |
| 0.50 — uninformative (control) | 0.0792 | 0.0100 |

Forcing the nominal truth is **166× worse** than the fitted value and 6× worse than carrying no strand
information at all. So κ and the per-node sense columns are mirrored *the same way*: the strand
likelihood scores a mirrored observation against a mirrored `p`, the mirror cancels, and the
deconvolution is right. **Only the exported scalar is mis-labelled.**

⚠ **What that costs today, and it is not nothing.** `rna_sense_frac` leaves calibration, reaches the QC
report and `cli.py`'s summary, and is the number `CARRY_FORWARD.md` §1 fact 17 quotes for the four real
libraries. Read as written, "κ = 0.002" says those libraries are almost purely antisense — a striking
claim about the biology. Read correctly, they are ordinary **highly stranded** libraries. Anyone
reasoning from the reported number is reasoning about a mirror.

**What to do, and the order matters.** (1) Find where the mirror lives — the scanner's genome-strand
convention or the simulator's read orientation — by checking one known-strand read's flags against
`sj_strand` and the accumulator column it lands in. (2) Fix the **naming or the orientation, not both**,
and gate it on the table above: after the fix, the fitted κ must be ≈ 0.99 *and* `f_gdna` must stay at
0.0030. A change that moves κ without preserving `f_gdna` has broken the inference to fix a label.
(3) Correct `CARRY_FORWARD.md` §1 fact 17 and §0 C4.

**Why deferred:** it is a labelling defect on a correct inference, so nothing downstream is wrong today.
⚠ It is nonetheless a **trap 27** in the making — one quantity, two meanings, and the prose disagreeing
with the code is exactly how this project lost months before.

## 7. ⚠ THE EM's END-TO-END VARIATION IS `assignment_mode="sample"` — BY DESIGN

⭐ **Owner, 2026-07-30: the EM draws each fragment's assignment from its posterior by default
(`EMConfig.assignment_mode = "sample"`), and that is deliberate, not a defect.** It can be set
deterministic. Measured, four runs of one scenario, same seed, same BAM, only the mode differing:

| `assignment_mode` | total assigned, 4 runs | spread |
|---|---|---|
| **`sample`** (the default) | 5440.7 / 5458.8 / 5465.3 / 5455.8 | ⚠ **~25, i.e. 0.5 %** |
| `map` | 6002.246879554673 / …736 / …669 / …744 | ~**1e-10** relative |
| `fractional` | 6276.853933515619 / …675 / …646 / …627 | ~**1e-11** relative |

⭐ **So an A/B has a switch: run it under `map` or `fractional` and the noise drops eight orders of
magnitude.** The ~1e-10 residual is float accumulation order, the same family as `CARRY_FORWARD.md` §1
fact 11, and is far below any effect these steps chase. **This is the practical answer for S5.a2, A6/A3
and the benchmark skill: A/B under a deterministic assignment mode, not the default.**

⚠ **The three modes give materially different totals** (5441 / 6002 / 6277 here) — they are different
estimators, not the same answer at different precision. An A/B must hold the mode fixed on both arms,
and a number quoted from one mode is not comparable to one from another.

### ⛔ AND THE ONE FAILING SCENARIO IS NOT FLAKINESS — it is mode-dependent LEAK

⚠ **Established at S6, and it changes the diagnosis.** `test_nrna_double_counting[g20_n0_s100]`'s
negative control (a silent transcript, truth 0, limit 25) does not fail *because* the answer wobbles —
every draw is over the limit:

| `assignment_mode` | `t_ctrl` counts | scenario suite (120 tests) |
|---|---|---|
| `sample` (default) | 29 / 29 / 32 / 29 / 33 | 1 failed |
| `fractional` | 30 | 1 failed |
| **`map`** | ≤ 25 | ⭐ **120 passed** |

⛔ **Do NOT switch the harness to `map` just because it goes green.** MAP assigns each fragment to its
single highest-posterior component, so it is the mode that most aggressively suppresses low-posterior
assignment — and a negative control is a **one-sided** metric (`CARRY_FORWARD.md` §3 trap 19: in a
library with none of X, any change that lowers X scores better). Going green under MAP would be fitting
the test to the mode. The real question is whether ~30 counts of leak onto a silent transcript is
acceptable, and that is a **modelling** call, not a testing one.

### ⛔ The separate defect: `EMConfig.seed` does not make `sample` reproducible

The measurement above passed `seed=12345` **explicitly** and `sample` still varied. The seed IS plumbed —
`estimator.py:187` builds `np.random.default_rng(em_config.seed)`, `:369` draws `rng_seed`, and
`em_solver.cpp:2165` mixes it per locus as `SplitMix64(rng_seed ^ li·φ)`. So a seeded run **should**
reproduce and does not. `tests/scenarios/conftest.py:92` already sets `EMConfig(seed=PIPELINE_SEED)` and
its negative control still returns 29 / 29 / 32 / 29 / 32.

**Candidates, in the order worth checking:** `rng_seed` is drawn once *per call* of the batch entry
point, so a varying number of calls or a varying batch split shifts every subsequent seed; or the
per-locus `li` is not stable across runs; or the ~1e-10 posterior jitter flips draws near a decision
boundary (this would explain a *small* rate, not 0.5 %). ⚠ Note `EMConfig.n_threads` defaults to `0` =
all cores, so a run that does not also set `OMP_NUM_THREADS` is parallel regardless.

**Why it matters even though the mode switch exists:** reproducibility is what `seed` is *for*, and a
seed that silently does not reproduce is worse than no seed — it invites exactly the "I ran it twice and
it agreed" reasoning that is not available here. **First step is characterisation, not a fix:** run one
scenario N times under `sample` with a fixed seed, record the per-transcript spread, and state the floor.

**Why deferred:** the mode switch removes the blocker for every downstream A/B, so this is now a
correctness-of-`seed` question rather than a gate on progress.

---

## Smaller, self-contained, and each with its reason

### ⚠ `build_scan_cache.py`'s skip logic does not cover the payload schema

It printed `cached / skip` for all eight pilot conditions that `read_scan_cache` then **refused**
(`payload_schema_digest None != '66b41ea0b645209d'`). `--force` is the workaround; the digest belongs in
the skip test, which currently checks only that the directory exists. Same family as
`CARRY_FORWARD.md` §3 trap 25 — a cache key that does not cover the artifact it caches. Cost so far: one
confusing failure and a 67 s re-scan. **Fix:** attempt the load in the skip branch and rebuild on refusal.

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

`gdna_em_count` fell 16–52 %. ⭐ S5 has landed, so the new suite can adjudicate it — but ⛔ **not by
regenerating**: S6 regenerating the goldens records whatever the tree now does, which is not the same as
validating it. Adjudication is a separate comparison against the suite's truth.
