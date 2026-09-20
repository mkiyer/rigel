# NEXT SESSION — the nascent RNA siphon: trace it to a root cause (2026-09-19, end of day)

This file is only how to begin. The ranked view is `docs/ROADMAP.md`, the open problems are
`docs/ISSUES.md`, the rulings and the record are `docs/DESIGN.md`, what "done" means is
`docs/SUCCESS.md`, the lessons are `docs/TRAPS.md` cited by name, the release procedure is
`docs/MANUAL.md` / `docs/PUBLISHING.md`.

## Where the tool is

Everything through `d8d563a2` is landed and pushed; the tree is clean, the suite is
**3,440 passed / 0 failed / 2 xfail / 3,442 collected**, `preflight.py --full` is green. Today landed
the RNA prior's restoration (`ISSUES: nascent-gets-no-rna-prior`, CLOSED — every ladder condition
improved, three xfails closed), its re-measurement, the per-scenario release report
(`quant_accuracy.py --markdown`) and the `ladder-report` skill that rebuilds and republishes it.

The deliverable, transcript-level Σ|Δ| as a share of the true annotated RNA at `g00`/`g05`/`g50`/`g98`:
unstranded OFF 1.7 / 1.9 / 2.5 / 17.8 %, stranded OFF 2.0 / 1.6 / 2.5 / 15.4 %, stranded ON
6.5 / 3.5 / 5.2 / 33.0 %, deferred 7.3 / 10.9 / 10.5 / 102.0 %.

## THE WHOLE JOB — `ISSUES: nascent-siphons-gdna-under-capture`

Read that entry first; it carries every number below and the two arms that already constrain the answer.

**What is measured and not in dispute.** Under capture the SYNTHETIC nascent entities take fragments
from gDNA almost one for one, with the annotated pool barely moving: at `g50 ss.99 ON` nascent is
+541,216 and gDNA −534,656 (annotated −6,560). Capture flips the sign — off capture gDNA takes from
nascent instead. `g00` capture-ON is a separate sub-case (no gDNA to take; nascent loses 132,921 to
annotated, and the true ruler fixes it outright) and must not be pooled with the rest.

**What the existing arms already rule out.** A perfect `LocusPriors` removes ~0 of it in scope
(`g50 ss.99 ON` 541,216 → 541,762). The simulator's own capture-aware lengths remove 13–28 %. So
neither the per-locus prior's magnitude nor the transcript ruler alone is the mechanism, and the bulk
is unexplained — that is the session's subject.

### The owner's two hypotheses, and the state of each

**(A) The capture contraction is derived from gDNA, which is unspliced, and nascent entities are
unspliced too.** A nascent entity is a single-exon span over the whole gene — geometrically
indistinguishable from gDNA, which is exactly what the per-piece capture efficiencies are measured on.
An annotated transcript is spliced, and a probe lying over a junction captures a molecule gDNA can
never produce there, so gDNA cannot witness that probe's effect on the mature isoform. Measured on
`capture_truth_on.npz`: the TRUE capture factor's median is **555.7 for annotated transcripts and 41.0
for synthetic entities**, a 13.5× gap, while a nascent entity's probed fraction is a median 0.0246 of
its span against 0.271 for an annotated one. PARTLY SUPPORTED — the true ruler removes 13–28 %, and all
of the `g00` sub-case. Not yet the bulk. Its twin is
`ISSUES: ruler-witness-geometry-on-transcript-panels`, the same witness geometry pointed at the isoform
split; start by reading that entry, because whatever is derived here may close both.

**(B) The prior lets silent synthetics rise as zombies and steal fragments.** UNTESTED IN ITS STRONG
FORM. The `oracle` arm is a perfect PER-LOCUS prior and it moves nothing, but the per-TRANSCRIPT
allocation is a different lever and its arm (`oracle_alloc_seed`) on disk is STALE — it predates the
prior's restoration and was never re-run. ⭐ Re-running it is the cheapest decisive move in the whole
session; do it first and let it run while reading. ⚠ Note what the restoration did and did NOT change:
the ABSORBING STATE survives, so a component with zero evidence still cannot be revived by prior mass
— but a nascent entity spanning a gene always has some coverage-weighted warm start, and the prior no
longer helps it decay (the rate fell from `kappa/(1 + P/A)` to `kappa = w_N/w_T`). "Zombie" here means
"decays too slowly", not "revived from exactly zero" (`EQUATIONS.md` §9b).

### A third candidate the closed entry leaves on the table

The EM gives each locus ONE gDNA rate spread uniformly along it, while under capture gDNA's density is
an order of magnitude higher at probed exons — which is exactly where the nascent entity also sits,
since it spans the whole gene. Inside a probed exon the model under-predicts gDNA and the surplus goes
to whichever component has the most opportunity there. Calibration already publishes per-piece
efficiencies, so the test is whether a position-dependent gDNA weight closes it. This came from
`ISSUES: em-overturns-the-calibrated-gdna-split` (now CLOSED into the siphon entry) and survives intact.

### How to work it

The debug loop, and the worst IN-SCOPE scenario is `g50 ss.99 ON` or `g98 ss.99 ON` — never the
deferred `ss_0.50` rows. Per-locus attribution first (the EM's gDNA against the certified truth, ranked
by mass) with `nrna_est` read beside it, then `ruler_vs_truth.py` per class and kind on the nascent
entities specifically, then `confusion.py` for per-fragment truth against assignment. ⛔ Derive on
paper, prototype outside the main tree, A/B against what ships on the same conditions, one mechanism at
a time. No magic numbers.

⛔ `calibration_vs_oracle.py`, `zero_controls.py` and `policy_benchmark.py` run no EM: they are the
CONTROLS for anything done here and must come back identical. If one moves, the change leaked upstream.

## The 2 xfails are proven defects, each deferred to its thread

`ISSUES: two-sided-exon-row`; `ISSUES: the-lower-bound-noise-ratchet`.

## What gates the release

⚠ `zero_controls.py` EXITS 1 on this tree and did before today too — byte-identical across the
restoration, so it is a standing state and not a regression: five controls off by more than 0.01 on a
constant truth, worst `TA_single_exon · ZERO gDNA` at 0.5931. The `g00` ladder rung is solved (library
gDNA fraction 0.0002–0.0005 against 0). The rest: the deliverable measured and not regressed per
stratum; the suite at its standing count; `preflight.py --full` green; the standing risks re-read
(`ISSUES: capture-degeneracy-standing-risk`, `ISSUES: flgap-panels-stale-nascent-model`); the manual
true of what ships.

## Standing rulings carried (unchanged)

* Float64 for the whole of ψ; ONE solver; ONE λ lattice; no θ lattice; the tilt's hypothesis space is
  {pure +, pure −, mixed}. The strand channel's liveness is a protocol decision on the spliced 2×2.
* The landscape trains only where a solve locates a slot; the ruler reads its located enriched mode and
  admits failure where there is no gDNA to read.
* Real data is a test input, never a design input. A few high-quality instruments, kept current. The
  source cites no doc. One production path: once native code is validated the Python it replaces goes.
* Unstranded × capture-ON stays DEFERRED. The fragment-length composition channel stays retired.
* The owner drives commits. The refit count and CI's on-demand trigger are the owner's.

## Where everything is

* THE ARMS: `~/Downloads/rigel_runs/suite/ladder/arms/` (base, base_reseed, oracle, oracle_ruler — all
  re-scored today on the restored prior), with the PRE-RESTORATION set preserved beside them in
  `arms_baseline_alpha0_2026-09-19/`. ⛔ `qa_ladder_oracle_alloc_seed.jsonl` is STALE — it was not
  re-run and must not be put in a table with the others.
* THE REPORT: `~/Downloads/rigel_runs/reports/ladder_accuracy_2026-09-19.md` and the Artifact at
  `https://claude.ai/artifact/Cek2wmKtitgbDyfM5gNqyj`. Rebuild and republish both with the
  `ladder-report` skill (`.claude/skills/ladder-report/`), which updates that same URL in place.
  ⛔ The report is a RENDERING: re-score the arms first or it will show stale numbers and look current.
* THE TRUE CAPTURE FACTOR: `~/Downloads/rigel_runs/suite/ladder/oracle_cache/capture_truth_on.npz`
  (`factor`, `probed_frac`, `L_plain` on the index's own 15,669-row axis, synthetics included).
* The snapshots of today's four commits:
  `~/Downloads/rigel_runs/prototypes/2026-09-19_nascent_prior/commits/`.
* The identity references `~/Downloads/rigel_runs/arms/review_identity_*.json` are NOT bit-identical on
  this tree any more — the restoration moves the EM's counts by design. Re-freeze before a rename.
