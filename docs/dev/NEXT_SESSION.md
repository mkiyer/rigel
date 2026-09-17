# NEXT SESSION — start here (2026-09-16, after the ruler's repair)

This file is only how to begin. The port is `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` §F; the θ thread's
derivations are `docs/dev/THETA_QUADRATURE.md`; the ruler's plan, EXECUTED, is `docs/dev/EXPECTATION_RULER_PLAN.md`.

## The agreed order of the sessions (owner, 2026-09-14; session 4 inserted 2026-09-16)

1. ~~Code review and cleanup~~ — DONE 2026-09-14.
2. ~~The test chromosome's new structures, both panels remeasured~~ — DONE 2026-09-14.
3. ~~The ruler at zero gDNA and the flux price's witness~~ — DONE 2026-09-15 (`aeb465fa`, `f5b7471f`).
4. ~~The ruler's repair~~ — DONE 2026-09-16, four commits PREPARED for your go (below).
5. **THIS ONE — the performance re-baseline and THE PORT** (`ISSUES: performance-memory-bounded-solve`, plan §F):
   two back-to-back profiler pairs at 8 threads on the deep library, identity references frozen fresh on that tree
   (they are fresh on this one: re-frozen after commit 3 with the reason in its message), `sweep_replay.py capture`
   (`sweeps_MO_3021_step6` predates the location floor, the ruler, the witness and the repair), then the port of the
   block solve, gated bit-identical on every stage. Take `ISSUES: multimapper-blind-support` right after the port: it
   is the first ruler question on real libraries and its repair is in the opportunity model, not the solve.

## What session 4 left, and how to commit it (owner's go)

Four commits are prepared, not made: the working tree holds all four and each one's file state is snapshotted with
its message under the session scratchpad's `commits/1_members`, `2_truth_instrument`, `3_expectation_ruler`,
`4_multimapper_check` (`commits/commit_series.sh` replays them in order and leaves the tree clean; `cat
commits/*/MESSAGE.txt` and `git diff` are the review). Every gate is on record in the messages: the falsification
tests verified failing first, every perturbation fired, the suite re-derived file by file (3,414 / 5 xfail / 3,419),
the goldens unchanged to the bit, the three `review_identity_*` references re-frozen after commit 3 with the reason,
both panels' instruments run on the landed tree against `DESIGN.md` §7's standing numbers.

1. **The reference's members** (`ISSUES: the-ruler-reference-on-sparse-real-libraries` CLOSED): members are
   kernels with a location, widths read among them, located iff more than √n members at a median knn width within
   1 nat; the regime on the result (`gdna_reference_members`, in `summary.json`). Both panels unchanged to the
   fragment; LBX0190 and MO_3021 `None → no contraction` (their walls-only landscapes had read a false reference).
2. **`ruler_vs_truth.py`** promoted from the prototype, `--self-test` 20/20, its index row and `TESTING.md` §0c; the
   depth-ladder configs given a home under `scripts/sim/configs/test_reference_depth_*.yaml`.
3. **The expectation ruler on the per-base length** (`ISSUES: ruler-multimapper-floor-caps-the-correction` CLOSED;
   `DESIGN.md` §7.2, `EQUATIONS.md` §11): efficiencies are posterior means under the landscape from own counts and
   apportioned crossings, published on the result; the transcript length a taper-weighted sum over its own bases;
   the locus prior's gDNA length the count's own objects at their efficiencies; floor, junction objects and flank
   imputation deleted. Every open question of the plan tried both ways, the refused arms with their numbers in §7.2.
   The ladder thermometer's capture-ON misassignment 0.45× stranded / 0.68× unstranded, OFF strata identical.
4. **The multimapper check** (`ISSUES: multimapper-blind-support` OPENED, read-only): on both captured real
   libraries the factor falls monotonically with the transcript's multimapper share (VCaP median 0.231 → 0.002 from
   below 1 % to ≥ 50 %); the repair is a mappable support in the opportunity model.

## 2026-09-17 addendum — the yield's endpoint, and a fifth commit prepared

The owner ruled TPM stays on the plain fragment-length-marginal length (kept constant for now) and asked for the
yield's variance to be modelled. The session's `s6/` scratchpad holds the work: `lever_census.py` (the pipeline
twice, contracted vs plain transcript yields), `census_flips.py`, `yield_draws.py` (every piece efficiency drawn
from its posterior; the count as the expectation over draws) and `draws_census.py`. What was found is in
`DESIGN.md` §7.2 (the yield's two consumers, and its endpoint) and two new entries, `ISSUES:
yield-variance-beside-the-count` (later, with the per-transcript prior lane: the analytic yield sd beside the
count, gated against the draws' spread) and `ISSUES: capture-premise-untested-on-cdna` (watch). **Commit 5 is this
commit** (`commits/5_floors`): the 1 bp
floors deleted for the cannot-emit rule, five gates, four perturbations fired, the thermometer neutral within the
reseed floor, the capture-ON identity reference re-frozen with the reason. The isoform flips the census found
(VCaP 39 % of multi-isoform genes) are the model's answer under the capture premise and are stable under the yield's
posterior — not a defect to fix, a limit to disclose (`count_unambig` beside `count`). The owner's order after it (2026-09-17):
the release ships this contraction as it stands; PERFORMANCE next — the tool must be blazingly fast and is painfully
slow — which is item 5 below; the variance metric and the prior's allocation across transcripts come with the
per-transcript prior lane, after.

## ⛔ Read first — three things the repair uncovered

* **The capture-OFF identity digest moves by one ulp** after commit 3 (one transcript's posterior mean, 1.1e-16 on
  the ladder's `g05 ss.50 OFF` row): the old floor's `w·span + (1−w)·span` arithmetic against the exact span. Every
  count and every calibration array is bit-identical; the references are re-frozen. Do not chase it.
* **A per-transcript Python loop over the real index is quadratic.** `transcript_piece_lengths` ran LBX0588 and
  VCaP for 70 minutes until `BaseTaper.interval_sums` took one template length per interval; now 160 s and 667 s
  against the shipped tree's 136 s and 616 s. Time every new per-transcript pass on `rigel_index` (457,371
  transcripts), never only on the ladder. Its perturbation found a hole in a green gate — no template length between
  one and two fragment lengths had been parametrised — closed with 600 and 998.
* **The ladder's panel spans junctions** (`ISSUES: ruler-witness-geometry-on-transcript-panels`, amended): the
  probed class reads 35 % within ±0.1 nat under every gDNA ruler, the certified true counts included, where the
  test chromosome's benign panel reads 99 %. A correction from the probe design is observable in principle; it is
  not this release's.

## The test chromosome moved (commit 3)

The tiny-exon block (`TESTING.md` §0a: ten 40 bp exons at 1,040 bp pitch, and the same run between two 1 kb exons,
unprobed and probed; 273 genes on 7.930 Mb, the budget unchanged) moved the test chromosome's benchmark rows: the
40 bp pieces are a new stress for the message layer (stranded ON 41,856 → 54,362, the ss 0.70 ON rows 64,174 →
105,997), the same on the shipped and the landed tree, so it is the substrate and not the mechanism. `DESIGN.md` §7
carries the new standing numbers; the superseded scenario sets are under
`~/Downloads/rigel_runs/test_reference_superseded_2026-09-16/`. The benign panel's rule gained one clause (a union
piece shorter than one probe gets a single probe centred on it).

## Decisions on record (unchanged, carried)

* Float64 for the whole of ψ; ONE solver; ONE λ lattice at a dimensionless step with a stated guarantee; no θ
  lattice anywhere; the tilt's hypothesis space is {pure +, pure −, mixed}.
* The strand channel's liveness is a protocol decision on the spliced 2×2; gDNA enters nowhere.
* The landscape trains only where a solve locates a slot; the ruler reads its located enriched mode, and admits
  failure where there is no gDNA to read (owner, 2026-09-15).
* Real data is a test input, never a design input; the four cfRNA libraries are re-run with the regime printed.
* A few high-quality instruments, kept current; no suite gate polices instruments; the source cites no doc.
* The message cache's on/off switch, the refit count and the scan's thread split are the owner's; parallelism
  waits for the port. CI runs on demand only.
* The certifier's FIELD gate flake on λ ≈ 7 boundaries is DEFERRED (owner, 2026-09-14).

## The session scratchpad (persists; nothing in the tree cites it)

Copy `/private/tmp/claude-503/-Users-mkiyer-proj-rigel/d290397d-5368-4194-86e0-5e1ccf118452/scratchpad/` to
`~/Downloads/rigel_runs/prototypes/2026-09-16_ruler_repair/` before it is lost: `commits/` (the four snapshots with
`MESSAGE.txt`, `snapshot.sh`, `commit_series.sh`), `s5/` (the arms `member_arms.py`, `ruler_arms.py`, `converge.py`,
the falsification harnesses, `gates3/` with every instrument's output on the landed tree, `identity_prev/` with the
identity references before commits 1 and 3), `real/` (`real_ruler.py`, `multimapper_check.py`, the four
`landed_<lib>/` outputs with `ruler.npz` and the check's tables), `head_tree/` + `shipped_site/sitecustomize.py`
(how the shipped tree is run beside an editable install: the scikit-build redirecting finder must be stripped).
