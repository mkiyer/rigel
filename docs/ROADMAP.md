# ROADMAP — the short ranked view

⭐ **WHAT THIS FILE IS: the brief overview and the ordered next steps — nothing else.** Three rules keep
it short: the SUBSTANCE of every item lives in `docs/ISSUES.md` (the issue log — OPEN entries plus the
append-only CLOSED/REFUSED record); **the changelog is git** (this file records no history); and the
NUMBERS POLICY (owner, 2026-08-22) — no figure lives here, a claim names the instrument that re-derives
it. Judging: `SUCCESS.md`. Rulings: `DESIGN.md`. Lessons: `TRAPS.md` (cite by NAME).

## THE 0.8.0 FRAME — owner ruling 2026-08-14; the full text is `DESIGN.md` §0b

Version on disk **0.7.1**; target **0.8.0**, a CALIBRATION release. **The metric is the calibration
result scored against oracle calibration** (`calibration_vs_oracle.py`, `solvability_audit.py`,
`prior_vs_oracle.py`) — the transcript number is a thermometer, never the ranking.

| stratum | 0.8.0 |
|---|---|
| unstranded × capture-OFF · stranded × capture-OFF · stranded × capture-ON | ⭐ **IN SCOPE** |
| unstranded × capture-ON | ⛔ **DEFERRED** — reported on every benchmark, never a development target; it carries most of the panel's error, so ⛔ never rank on a pooled total (`TRAPS: never-pool-the-strata`) |

⛔ The fragment-length **composition channel** is RETIRED until after 0.8.0 — do not propose it.
⭐ The ladder gives gDNA and RNA EQUAL fragment lengths on purpose: a gap lets the EM split origins on
length alone and mask calibration bugs (`TRAPS: a-length-gap-bypasses-calibration`).

## WHERE THE TOOL IS — one line per claim; run the named instrument for a current number

- **Library gDNA fraction**: accurate on the three in-scope strata, structurally BLIND on the deferred
  one (at κ = ½ no channel reaches an AMBIG slot; the θ-independent-channel search is CLOSED) —
  `solvability_audit.py`, `relay_pool_ab.py`.
- **Transcript assignment**: a sixth to a fifth of RNA fragments misassigned even under a perfect prior —
  calibration and assignment are two problems in two files — `quant_accuracy.py` (the thermometer).
- **Stage A (the accumulator)**: DONE; the fragment ledger closes exactly — `calibration_oracle.py`.
- **Fragment lengths**: CLOSED, both halves (2026-08-31). gDNA: two estimands, deconvolved by the
  two-pool contrast, coupled by resolution (`calibration/fl.py`, `gdna_density.py`; gates in
  `test_fl_realized.py`, `test_gdna_density.py`). RNA: the law is sound as shipped —
  `ISSUES: the-rna-length-law-fix` (CLOSED). Watch item: `ISSUES: capture-degeneracy-standing-risk`.
- **gDNA strand overdispersion**: robust to the annotation (`EQUATIONS.md` §6a–§6c, `DESIGN.md` §3.3a);
  on real data read `clamped_at_ceiling` and `effective_seeds`, never the bare value.
- **The message layer**: `relay` ships (propagation ON); `silent` is the measured floor; `MessagePolicy`
  ≡ silent at rung 0. The bar: **win on unstranded, minimal harm on stranded, never pooled** —
  `policy_benchmark.py`.
- **ψ**: the composition closes structurally on every published object (`test_composition_closes.py`);
  the reference location is DELETED (owner, 2026-08-24; the surviving form is `DESIGN.md` §6b.1); the
  λ-bracket widening is built and ships OFF — `ISSUES: psi-lambda-bracket-unshipped`.
- **The prior assembler**: with perfect masses its own error is parts per thousand —
  `prior_vs_oracle.py`, `mass_prior_ab.py`.
- **Two open defects reach the 0.8.0 metric**: the `g00` effective-length shrinkage
  (`ISSUES: g00-shrinkage-upstream-repair`) and the never-passed per-transcript prior lane
  (`ISSUES: per-transcript-prior-lane`).
- **Panels**: the sparse-nascent 16-condition ladder + the 30-condition test chromosome (anchored twin
  block), both cached and certified — `panel.py status`; the fl-gap side panels still carry the retired
  nascent model — `ISSUES: flgap-panels-stale-nascent-model`. ⚠ The ladder's nascent level is a
  DEVELOPMENT STRESS, not real data (`DESIGN.md` §0b).
- **Oracle FIELD certification**: 16/16 stamped, but the uniformity gate is vacuous on capture-ON and
  zero-gDNA rows — read the stamp with its vacuity flag — `calibration_oracle.py`.
- **Attribution floor**: the deliverable is not reproducible by default; no `quant_accuracy` delta below
  a few thousand fragments is attributable — re-derive `--arm base_reseed` in the same session.
- **Reading rules**: rank per stratum; quote `mwae_all`/Σ|err| and the SHIPPED column, never
  `solv%`/pass-0 (`TRAPS: the-intermediate-is-not-the-deliverable`).

## ⭐⭐⭐ NEXT — the recommended order (audited 2026-08-31)

**The method stays the owner's dissection loop** (2026-08-20): run the panel → worst IN-SCOPE scenario →
rank its objects by error mass (`worst_objects.py`, `calibration_walk.py`) → find the mechanism → gated
fix → **add the offending transcripts to the test chromosome** → re-run → repeat.

✅ **Measurement integrity is DONE (2026-08-31)**: the frame is ruled and migrated through the
loaders (`ISSUES: instruments-calibrate-undrained-cache` — only wave-3 bank-readers remain, migrate
as touched), `slot_truth` is re-certified, and **the in-scope baselines are re-derived in the drained
frame with same-day noise floors** — five instruments, whole ladder; run any of them for a current
number. ⭐ The re-derived facts a ranking can lean on: the worst IN-SCOPE scenario is
`g98 ss.99 capture-ON` on the 0.8.0 metric (almost entirely UNDER-called gDNA at near-pure objects —
the shortfall pattern holds on every stratum) and `g50 ss.50 capture-OFF` on `solvability_audit`'s
Σ|err| — dissect both before choosing; the relay hurts the solvable set on 9/12 contaminated
conditions while winning the blind rows and all four zero controls.

1. ⭐ **The drain's certified-RNA leak** — `ISSUES: drain-contaminates-certified-rna` (owner,
   2026-08-31: avoid all sources of false positives). Derive a repair; do not tune one.
2. ⭐ **The dissection loop on the worst in-scope scenario** (named above), answering
   `ISSUES: message-value-for-blind-slots` on the way (no solver needed; the anchored twin block was
   built for it).
3. ⭐ **The calibration build thread**: `ISSUES: measured-prior-rung-4` under the
   `ISSUES: reference-prior-refuted-at-concept-level` constraint — the one open BUILD aimed directly at
   the 0.8.0 metric — with `ISSUES: landscape-trains-on-real-substrate` as its payoff check.

Then, in standing order: `ISSUES: per-transcript-prior-lane` (the largest in-scope end-to-end lever,
now with a measured target population) · `ISSUES: performance-memory-bounded-solve` (owner: mandatory
before 0.8.0) · the ruler pair `ISSUES: u-ruler-arm` + `ISSUES: g00-shrinkage-upstream-repair` · the
capture-ON length pair `ISSUES: capture-blind-gdna-divisor` + `ISSUES: eb-shrinkage-magic-ess` · the
message-vs-prior questions `ISSUES: refit-vs-message-arbitration`,
`ISSUES: prior-fidelity-vs-deliverable`.

**Later / parked** (each has its entry): `expand-the-gdna-spectrum` · `psi-lambda-bracket-unshipped` ·
`alt-splice-rung-unverified` · `transfer-variance-premise` · `nascent-stress-sensitivity` ·
`relay-od-r-discontinuity` · `f32-strand-tilt-at-half` · `hygiene-ledger` ·
`oracle-effective-length-diagnostic` · `flgap-panels-stale-nascent-model` · `rename-the-drain` · `the-cancelling-pair`
(refused twice) · `crossing-pool-contrast` (blocked) · `parked-capture-pilot-sign` ·
`pure-rna-mirror-asymmetry` · `capture-degeneracy-standing-risk`.

## ⛔ DELIBERATELY NOT NEXT

The length composition channel (retired until after 0.8.0) · anything whose only target is the DEFERRED
stratum · every mechanism in `ISSUES.md`'s **CLOSED / REFUSED** section — read it before proposing
anything, because each entry is a build that was measured and turned down, with the number that killed it.
