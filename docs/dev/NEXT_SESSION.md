# NEXT SESSION — start here (2026-09-19, end of day: the EM's gDNA split under capture, now a nascent-vs-gDNA question)

This file is only how to begin. The ranked view is `docs/ROADMAP.md`, the open problems are `docs/ISSUES.md`,
the rulings and the record are `docs/DESIGN.md` (§0b carries the scope, its 2026-09-19 amendment and the
restoration ruling), what "done" means is `docs/SUCCESS.md`, the lessons are `docs/TRAPS.md` cited by name, and
the release procedure is `docs/MANUAL.md` / `docs/PUBLISHING.md`.

## Where the tool is

Everything through `61e82758` is landed and pushed. On top of it, **two commits are prepared as snapshots and
await the owner's go** — nascent RNA's share of the EM's RNA prior, restored, and the re-measurement it forced.
The suite is 3,432 passed / 2 xfail / 3,434 collected, `preflight.py --full` is green, the tree is otherwise clean.

The deliverable on the rebuilt ladder, under fractional assignment, transcript-level Σ|Δ| as a share of the true
annotated RNA at `g00` / `g05` / `g50` / `g98` (the pre-restoration reading in brackets):

| stratum | numbers |
|---|---|
| unstranded × OFF | 1.7 / 1.9 / 2.5 / 17.8 % [3.1 / 3.4 / 4.2 / 35.0] |
| stranded × OFF | 2.0 / 1.6 / 2.5 / 15.4 % [3.3 / 2.9 / 3.9 / 24.8] |
| stranded × ON | 6.5 / 3.5 / 5.2 / 33.0 % [6.8 / 4.5 / 7.8 / 109.8] |
| unstranded × ON (deferred) | 7.3 / 10.9 / 10.5 / 102.0 % [8.3 / 12.4 / 20.8 / 739.5] |

All 16 conditions improved and gene level fell 3.4–5.8× in scope. The reseed floor is 0–2,395 fragments in scope
(0.03 points; it was 0–416, not the 0–82 the ROADMAP used to claim). A perfect prior still recovers nothing in
scope. ⛔ The capture-OFF magnitude is read at the panel's 20.2 % nascent stress share; realistic is ~4.2 %.

## What landed (the two prepared snapshots)

`ISSUES: nascent-gets-no-rna-prior` is CLOSED. The EM's RNA pseudocount now goes to every RNA component in
proportion to the evidence it carries, none singled out for zero; the eligibility test, `component_is_synthetic`,
the `t_is_synthetic` lane and the `index` parameter it was the only reader of are gone (`EQUATIONS.md`
§9b–§9b.2, `DESIGN.md` §0b). THE WEIGHT is `w_i = raw[i]` — the shipped weights, every component admitted — and
not an equal share, because that keeps the absorbing state at zero evidence for free and leaves
`component_rna_prior_weight` free for the per-transcript lane. It also closed
`ISSUES: nested-antisense-leak-under-the-sane-ruler` (both rungs): 3 xfails → 0, 5 → 2 overall.

## ① THE FIRST JOB — the EM's gDNA split under capture (`ISSUES: em-overturns-the-calibrated-gdna-split`)

The entry SPLIT BY SIGN when it was re-measured, and the half that is left is a different question from the one
the old entry asked.

**The capture-OFF half largely closed**, and its cause was the allocation rule rather than the EM's arbitration:
`g50 ss.50 OFF` reads 0.5127 against 0.50 where it read 0.574. Nothing is owed there.

**The capture-ON half got WORSE**: `g50 ss.99 ON` reports 0.4465 against 0.50 (it reported 0.4793), `g98 ss.99 ON`
0.9189 against 0.98. The mass did not go back to the annotated transcripts — it went to the NASCENT channel, which
now over-calls 4.6× at `g50 ss.99 ON` (691,648 against 150,432 true) and 100× at `g98 ss.99 ON`. So the retired
α = 0 rule had been MASKING a nascent-vs-gDNA competition under capture by suppressing one of the two
competitors, and that competition is now the whole of this entry. ⚠ Unlike the capture-OFF half this is not only
a stress reading: capture-ON runs at a 2.6 % nascent fragment share against a realistic 4.2 %.

The standing hypothesis is unchanged and now sharper: the EM gives each locus ONE gDNA rate spread uniformly
along it, while under capture gDNA's density is an order of magnitude higher at probed exons — which is exactly
where the nascent entity also sits, since it spans the whole gene. Inside a probed exon the model under-predicts
gDNA and the surplus goes to the component with the most opportunity there. Calibration already publishes
per-piece efficiencies, so the test is whether a position-dependent gDNA weight closes it. Start with a per-locus
attribution on `g50 ss.99 ON` (the EM's gDNA against the certified truth, ranked by mass) and read `nrna_est`
beside it — the nascent channel is where the mass went, and `quant_accuracy.py`'s own comment calls it "a THIRD
false-positive channel with nothing to cancel it". The per-fragment instrument is `confusion.py` beside the
dissection data (below).

Then `ROADMAP.md` items ② (the capture ruler, which waits on an owner decision — now worth 5.2 points at `g00`,
1.8 at `g05` and 2.2 at `g50`) and ③ (the per-transcript lane).

## The 2 xfails are proven defects, each deferred to its thread

`ISSUES: two-sided-exon-row`; `ISSUES: the-lower-bound-noise-ratchet`. Closing one means repairing the thing or
asserting the invariant structurally, never widening a bound.

## What gates the release (`ROADMAP.md`'s last item)

`PUBLISHING.md` is two commands and a wait; the STATE is what gates it. ⚠ `zero_controls.py` EXITS 1 on this
tree and did so before this session too — byte-identical across the restoration, so it is a standing state and
not a regression: five controls off by more than 0.01 on a constant truth, worst `TA_single_exon · ZERO gDNA`
at 0.5931. The `g00` rung of the ladder is solved (library gDNA fraction 0.0002–0.0005 against 0). The rest:
the deliverable measured and not regressed per stratum; the suite at its standing count with an empty failure
set; `preflight.py --full` green; the standing risks re-read (`ISSUES: capture-degeneracy-standing-risk`,
`ISSUES: flgap-panels-stale-nascent-model`); and the manual true of what ships.

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

* THE PREPARED SNAPSHOTS: `~/Downloads/rigel_runs/prototypes/2026-09-19_nascent_prior/commits/` —
  `1_nascent_prior_restored` and `2_remeasured`, each with `FILES.txt`, `DELETED.txt` and `MESSAGE.txt`.
  `commit_series.sh` replays them in order.
* THE RE-MEASUREMENT: `~/Downloads/rigel_runs/suite/ladder/arms/` (the four `quant_accuracy` arms, re-run), with
  the PRE-RESTORATION baseline preserved beside it in `arms_baseline_alpha0_2026-09-19/`. Logs are
  `~/Downloads/rigel_runs/logs/ladder_*_nascentprior_2026-09-19.log`.
* THE CONTROLS, all confirmed identical across the change: `calibration_vs_oracle.py` (every metric column),
  `zero_controls.py` (byte-identical) and `policy_benchmark.py --panel ladder` (identical but its wall-clock
  line). They run no EM, which is why they are the right controls.
* Captures and reports: `perf/sweeps_VCaP_step19`; `perf/plan_final_2026-09-19/`.
* THE STRANDED × ON DISSECTION: `~/Downloads/rigel_runs/arms/2026-09-19_stranded_on/` — `confusion_*` (per-fragment
  truth against assignment at `g50`) and the scratch runner `dissect_run.py` / `dissect_analyze.py` /
  `confusion.py`. ⚠ Measured under the RETIRED allocation; re-run it before quoting a number from it.
* The identity references: `~/Downloads/rigel_runs/arms/review_identity_*.json`. ⚠ They are NOT bit-identical on
  this tree any more — the restoration moves the EM's counts by design. Re-freeze before the next rename.
