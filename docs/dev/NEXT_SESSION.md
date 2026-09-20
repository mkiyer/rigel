# NEXT SESSION — the nascent siphon's ROOT CAUSE is found; the repair is the job (2026-09-19, end of day)

This file is only how to begin. The ranked view is `docs/ROADMAP.md`, the open problems are
`docs/ISSUES.md`, the rulings and the record are `docs/DESIGN.md`, what "done" means is
`docs/SUCCESS.md`, the lessons are `docs/TRAPS.md` cited by name, the release procedure is
`docs/MANUAL.md` / `docs/PUBLISHING.md`.

## Where the tool is

Everything through `7eab2b47` is landed and pushed. Today's work is PREPARED AS SNAPSHOTS AND AWAITS THE
OWNER'S GO (`~/Downloads/rigel_runs/prototypes/2026-09-19_siphon_root_cause/commits/`). The suite is
**3,449 passed / 0 failed / 2 xfail / 3,451 collected**. ⛔ NOTHING IN `src/` CHANGED TODAY — the only code
change is one instrument, `scripts/design/quant_accuracy.py`, so every deliverable number stands exactly as
`docs/ROADMAP.md` records it and the release report needs no rebuild.

## THE FINDING — `ISSUES: nascent-siphons-gdna-under-capture`

Read that entry; it carries every number. In one line:

**`theta_n = 0` is an UNSTABLE fixed point of the shadow-vs-gDNA contest.** The EM gives a whole MultiLocus
ONE gDNA component with ONE opportunity `L_g` over the entire connected component, while every synthetic
nascent entity carries only its own gene's span `L_n`. A connected component is a union of gene spans, so
`L_g > L_n` STRUCTURALLY and grows with the component's gene count. Near zero the shadow's density is
multiplied by `L_g/L_n > 1` every iteration, so a shadow holding NOTHING climbs off zero and settles where
only the strand channel and `gdna_prior` stop it. The threshold is exactly 1, derived not chosen
(`EQUATIONS.md` §9b), and it is pinned against the SHIPPED solver in `tests/test_estimator.py`.

⭐ **Capture does not reverse the arbitration.** The same channel leaks 257,002 fragments off capture and
452,854 on it. The library-level sign flip is arithmetic: capture raises `L_g/L_n` 2.7× (3.57 → 9.70,
mass-weighted) and collapses the TRUE nascent pool 6.7× (1,013,400 → 150,405), so the compensating
under-call on LIVE shadows that was masking the false positive disappears.

## The owner's question, answered (2026-09-19, second round)

**"If calibration can tell us there is no evidence for a shadow, can we give it zero prior?"** Calibration
CAN tell us, it is RIGHT, and the EM never receives it.

* `gdna_prior_count` is accurate to **±1 % on every in-scope condition**; `rna_prior_count` to ±1 % on
  three of four (`g98 ON` over-states 64.5 %). At the INTRON-ONLY regions — the only place nascent RNA can
  sit — calibration says 5,951 fragments against a true 92, while the EM's shadows hold 476,270.
* On the 20 loci carrying the most shadow false positives at `g98 ss.99 ON` the EM **discards 196,217
  gDNA fragments calibration had right** (161,476 at `g50 ss.99 ON`).
* ⛔ `rna_prior_weight` is PLUMBED end to end and **NOTHING IN `src/` FILLS IT**. The solver therefore
  always takes `w_i = raw[i]`, which echoes its own belief and cannot contradict it.

## What the next session should do — two priced candidates, and the choice is the owner's

Judged on `quant_accuracy.py` per stratum above `--arm base_reseed`, fractional, with
`calibration_vs_oracle.py` / `zero_controls.py` / `policy_benchmark.py` as the controls:

1. ⭐⭐ **A PER-GENE gDNA OPPORTUNITY.** It needs NO new information: the severity is set by `L_g/L_n`, and
   at the ratio a per-gene opportunity gives (1.25) the shipped solver leaves the shadow at **0.00 %**
   with the gDNA pseudocount it ALREADY receives, against 45–49 % at the panel's mass-weighted 9.7. It
   changes `LocusPriors`, so all three controls become live — the expensive option, and the clean one.
   Derive before prototyping.
2. ⭐ **THE MEASURED PER-TRANSCRIPT PRIOR, on the tested components only.** A component reaching an object
   no other component's structure reaches has an independently measurable mass; one whose opportunity is
   wholly shared has none. That test selects 97.4 % of shadow spans and **0 % of annotated transcripts**.
   Priced (diagnostic, the untested half circular): `g50 ss.99 ON` +541,216 → **+22,187**, 96 % of the
   siphon, for 5.21 → 5.85 % of transcript error; `g50 ss.99 OFF` unharmed (2.50 → 2.53 %). ⛔ It needs
   `raw[i]`, which lives in the kernel, so the shippable form is a small `em_solver.cpp` change, NOT a
   Python producer filling the static lane — filling that lane wholesale reallocates the entire RNA
   pseudocount (2,793,710 fragments at `g50 ss.99 ON`) and takes the transcript table to 53 %.
   ⛔ `g98` first needs `rna_prior_count`'s +64.5 % over-call fixed: its pseudocount (180,806) is nearly
   the whole true RNA (194,011) and the shadows were acting as its SINK, so removing them without fixing
   it just moves the error to the transcript table (33.0 → 65.7 %).

⛔ **NOT a stronger gDNA pseudocount.** Closing the channel that way needs 1× the data at ratio 2, 5–10×
at 6.2 and 50× at 20 — and a prior many times the data is not a prior.

⛔ **NOT A LENGTH KNOB.** Scaling the shadows' EM length by 2 removes 97 % of the siphon and returns
437,761 fragments to gDNA — and makes the transcript table WORSE (252,376 → 386,614 Σ|Δ|), because it kills
the live entities with the dead ones. That probe is a falsification instrument and is recorded as one.

⚠ **AND SIZE IT FIRST.** `ISSUES: nascent-stress-sensitivity` is now the gating question, not a footnote:
the whole defect is a competition between gDNA and a nascent pool the ladder runs at a 0.50 `on_fraction`
DEVELOPMENT STRESS. Re-simulate `g50 ss.99 ON` at the realistic 0.10 and re-read the siphon before paying
for a mechanism — a repair worth 541,216 fragments at stress may be worth a fifth of that in the expected
case (`DESIGN.md` §0b).

## What is unexplained, and should not be quietly dropped

* **The damping factor.** The bare two-component fixed point predicts 44–53 % of the shadow-exclusive pool
  at the measured geometry; the panel leaks 21.1 % (ON) and 11.6 % (OFF). The mature isoforms competing at
  exonic positions and the per-locus gDNA pseudocount are the two damping forces and neither is closed
  quantitatively.
* **The LIVE shadows' under-call off capture** (−323,623 at `g50 ss.99 OFF`), which is what masks the false
  positive there. Not the same defect and not yet its own entry.
* **A perfect allocation makes `g98` capture-ON markedly WORSE** (33.0 → 49.6 %, and 102.0 → 321.9 % on the
  deferred stratum). Both the broken and the repaired arm show it. Unexplained.

## The instrument defect found on the way — read before trusting ANY oracle arm

`quant_accuracy.truth_weights` read `observed_mrna_fragments`, identically 0 on all 6,919 SYNTHETIC rows
(the nascent truth is `observed_nrna_fragments`), so `--arm oracle_alloc*` handed every shadow a weight of
ZERO and was re-running the retired `alpha = 0` rule under an oracle's name — tracking the pre-restoration
baseline to within 15 % on every in-scope condition. Fixed, gated three ways, each watched to fire.
`TRAPS: an-oracle-column-that-omits-a-population`. ⛔ The arm's headline in
`ISSUES: per-transcript-prior-lane` was re-measured; the old numbers are retired.

## The 2 xfails are proven defects, each deferred to its thread

`ISSUES: two-sided-exon-row`; `ISSUES: the-lower-bound-noise-ratchet`.

## What gates the release

⚠ `zero_controls.py` EXITS 1 on this tree and did before today — reproduced byte-identically this session,
so it is a standing state and not a regression: five controls off by more than 0.01 on a constant truth,
worst `TA_single_exon · ZERO gDNA` at 0.5931. The rest: the deliverable measured and not regressed per
stratum; the suite at its standing count; `preflight.py --full` green; the standing risks re-read
(`ISSUES: capture-degeneracy-standing-risk`, `ISSUES: flgap-panels-stale-nascent-model`); the manual true of
what ships.

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

* THE ARMS: `~/Downloads/rigel_runs/suite/ladder/arms/`. `base`, `base_reseed`, `oracle`, `oracle_ruler`
  are UNCHANGED from 2026-09-19 and current. ⭐ `qa_ladder_oracle_alloc_seed.jsonl` was RE-SCORED today on
  the repaired instrument and is the only arm whose numbers moved; the two superseded versions are kept
  beside it as `STALE_*.bak` (pre-restoration) and `*_MATURE_COLUMN_ONLY_*.bak` (the defect, on the current
  tree) — ⛔ neither may go in a table with the others.
* THE REPORT: `~/Downloads/rigel_runs/reports/ladder_accuracy_2026-09-19.md` and the Artifact at
  `https://claude.ai/artifact/Cek2wmKtitgbDyfM5gNqyj`. ⭐ NOT REBUILT TODAY AND CORRECTLY SO: it renders
  `base` / `base_reseed` / `oracle` / `oracle_ruler`, none of which moved, so a rebuild would republish the
  same numbers. Rebuild with the `ladder-report` skill the moment a repair lands.
* THE INVESTIGATION, kept: `~/Downloads/rigel_runs/prototypes/2026-09-19_siphon_root_cause/investigation/`
  — the per-locus dissection (`dissect2.py` and eight `d2_*.npz`, one per in-scope condition), the analyses
  (`an1`–`an10`), the closed-form fixed point (`derive.py`), the solver toy (`shadow_toy.py`), the λ probe
  (`lambda_probe.py`, with its two arm files) and the corrected-allocation prototype (`alloc_fixed.py`).
  Every number in `ISSUES: nascent-siphons-gdna-under-capture` is re-derivable from these.
* The snapshots: `~/Downloads/rigel_runs/prototypes/2026-09-19_siphon_root_cause/commits/`.
* The identity references `~/Downloads/rigel_runs/arms/review_identity_*.json` are NOT bit-identical on
  this tree — the RNA prior's restoration moved the EM's counts by design. Re-freeze before a rename.
