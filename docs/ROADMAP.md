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
  `solvability_audit.py`, `message_pool_ab.py`.
- **Transcript assignment**: a sixth to a fifth of RNA fragments misassigned even under a perfect prior —
  calibration and assignment are two problems in two files — `quant_accuracy.py` (the thermometer).
- **Stage A (the accumulator)**: DONE; the fragment ledger closes exactly — `calibration_oracle.py`.
- **Fragment lengths**: CLOSED, both halves (2026-08-31). gDNA: two estimands, deconvolved by the
  two-pool contrast, coupled by resolution (`calibration/fl.py`, `gdna_density.py`; gates in
  `test_fl_realized.py`, `test_gdna_density.py`). RNA: the law is sound as shipped —
  `ISSUES: the-rna-length-law-fix` (CLOSED). Watch item: `ISSUES: capture-degeneracy-standing-risk`.
- **gDNA strand overdispersion**: robust to the annotation (`EQUATIONS.md` §6a–§6c, `DESIGN.md` §3.3a);
  on real data read `clamped_at_ceiling` and `effective_seeds`, never the bare value.
- **The message layer**: `transfer` SHIPS (the default since 2026-09-09; propagation ON), on the
  two-phase backbone with the level lanes (`DESIGN.md` §6b.12–§6b.14); `silent` is the measured floor;
  the relay and its baggage are deleted (2026-09-09; git carries them). The bar: **win on unstranded,
  minimal harm on stranded, never pooled** — `policy_benchmark.py --panel ladder`. The retired relay
  kept three zero-gDNA rows, which are the landscape's
  (`ISSUES: gdna-landscape-trains-on-false-positives`) — rank 1 below.
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

## ⭐⭐⭐ NEXT — the recommended order (audited 2026-09-09, the day the message layer merged)

**The method stays the owner's dissection loop** (2026-08-20): run the panel → worst IN-SCOPE scenario →
rank its objects by error mass (`worst_objects.py`, `calibration_walk.py`) → find the mechanism → gated
fix → **add the offending transcripts to the test chromosome** → re-run → repeat.

⭐ **The message layer is FINISHED** (`transfer` ships; the rulings `DESIGN.md` §6b.4–§6b.14; the
cleanup on the one policy done and gated bit-identical). ⭐ The facts this ranking leans on, each with the
instrument that re-derives it: every in-scope stratum's per-object composition error is about one
percent (`calibration_vs_oracle.py`) and every in-scope row misplaces under 2.5 % of its fragments
(`policy_benchmark.py --panel ladder --by-class`), so messages are at diminishing returns in scope;
the remaining in-scope error sits on the intron's own solve (the largest class by mass, 1–2 % of its
own fragments), on exon|intron boundaries at their counting floor, and on near-pure objects at high gDNA
(the vertex atom, rank 4); the zero-gDNA rows' false positives are the landscape prior's
(`ISSUES: gdna-landscape-trains-on-false-positives`); and the two rankers still disagree on the worst
IN-SCOPE scenario because they weigh objects differently, so dissect BOTH before choosing one.

1. ⭐⭐⭐ **THE gDNA LANDSCAPE PRIOR (owner priority 3, now first).** It owns the largest remaining in-scope
   numbers — the zero rows' false positives and the walled and terminus classes' residuals — and the
   vertex ceiling (rank 4) says its training population is where the in-scope value sits. The training
   population per node class on the zero rows first, then the estimator at bound-only nodes
   (`ISSUES: gdna-landscape-trains-on-false-positives`, `ISSUES: measured-prior-rung-4`,
   `ISSUES: landscape-trains-on-real-substrate`, under the `ISSUES: reference-prior-refuted-at-concept-level`
   constraint).
2. ⭐⭐ **THE POST-CALIBRATION, PRE-EM SETUP (owner priority 4).** `priors.py` / `result.py` / `derive.py`
   against `prior_vs_oracle.py` and the ruler column: `ISSUES: prior-fidelity-vs-deliverable`,
   `ISSUES: eb-shrinkage-magic-ess`, `ISSUES: g00-shrinkage-upstream-repair`, `ISSUES: capture-blind-gdna-divisor`,
   `ISSUES: per-transcript-prior-lane`, `ISSUES: u-ruler-arm`. First step: re-run `prior_vs_oracle.py` so the
   assembler's own error is re-recorded under the shipped policy before anything moves.
3. ⭐ **IMPROVE THE MESSAGE POLICY (owner priority 2), only where a row is above the bar.** One prototype
   arm at a time through `policy_prototype.py --module`, judged on the three panels and the ladder, halves
   apart, pass zero beside the pipeline: `ISSUES: flux-price-witness-units` (the open defect the landed
   flux level shares), `ISSUES: two-sided-exon-row`, `ISSUES: flux-floor-dispersion`,
   `ISSUES: ambig-node-as-a-gdna-source`, `ISSUES: message-layer-open-cases` (the substrate: `div`, the
   antisense's nascent variant, `nest`).
4. ⭐ **The vertex atom, PRICED under the shipped policy (2026-09-09, `vertex_ceiling.py`, re-pointed to
   the two-phase solve and to the PARAMETER-vertex population)**: worth within 1 % on every stranded
   in-scope row, a few percent of the unstranded capture-OFF rows rising with gDNA, and a quarter to a
   third of the DEFERRED stratum. The in-scope value sits on silent genes' regions and nascent-free
   introns — rank 1's training population — and arrives through the refit prior, not at pass zero;
   `EQUATIONS.md` §9a.1/§9d.4 carry the derivation and the spike with no new constant. Re-derive with
   `vertex_ceiling.py --arm base|vertex_free --oracle-cache …` and `--compare`.

Then, in standing order: `ISSUES: performance-memory-bounded-solve` (owner: mandatory before 0.8.0) ·
the message-vs-prior question `ISSUES: refit-vs-message-arbitration` (the pre-EM entries above sit under
rank 2).

**Later / parked** (each has its entry): `expand-the-gdna-spectrum` · `psi-lambda-bracket-unshipped` ·
`transfer-variance-premise` · `nascent-stress-sensitivity` ·
`f32-strand-tilt-at-half` · `hygiene-ledger` ·
`oracle-effective-length-diagnostic` · `flgap-panels-stale-nascent-model` · `rename-the-drain` · `drain-contaminates-certified-rna` (parked 2026-09-01: the ceiling refused the in-solve correction; two recorded follow-ups) · `the-cancelling-pair`
(refused twice) · `crossing-pool-contrast` (blocked) · `parked-capture-pilot-sign` ·
`pure-rna-mirror-asymmetry` · `capture-degeneracy-standing-risk`.

## ⛔ DELIBERATELY NOT NEXT

The length composition channel (retired until after 0.8.0) · anything whose only target is the DEFERRED
stratum · every mechanism in `ISSUES.md`'s **CLOSED / REFUSED** section — read it before proposing
anything, because each entry is a build that was measured and turned down, with the number that killed it.
