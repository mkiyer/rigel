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

## ⭐⭐⭐ NEXT — the recommended order (audited 2026-09-01)

**The method stays the owner's dissection loop** (2026-08-20): run the panel → worst IN-SCOPE scenario →
rank its objects by error mass (`worst_objects.py`, `calibration_walk.py`) → find the mechanism → gated
fix → **add the offending transcripts to the test chromosome** → re-run → repeat.

⭐ **The measurement frame is settled and the baselines are current** — the ruling is `DESIGN.md`
§4.3 (only wave-3 bank-readers remain: `ISSUES: hygiene-ledger`), `slot_truth` is certified in that
frame, and every in-scope baseline was re-derived on the whole ladder with a same-session noise floor.
⭐ The facts a ranking can lean on, each with the instrument that re-derives it: the worst IN-SCOPE
scenario differs between the two rankers — `calibration_vs_oracle.py` and `solvability_audit.py`
disagree because they weigh objects differently, so dissect BOTH before choosing one; the error is
almost entirely gDNA UNDER-called at near-pure objects, on every stratum; and the shipped policy leaves
the zero-gDNA controls' walled and both-stranded pieces to the landscape prior.

⛔ **Two local repairs were run to a verdict on 2026-09-01 and BOTH are refused** — the arcsine
magnitude coordinate (`ISSUES.md` CLOSED / REFUSED carries the full record) and, priced in the same
session, a PERFECT vertex answer at every reachable object with the whole chain re-solved, which nets
roughly a wash. ⭐⭐ **That is what re-ranked the list below: a correct local answer does not survive
propagation, so no local repair — coordinate, prior, or atom — can pay off until the message layer
is fixed.** Re-derive with `calibration_walk.py` (rung E vs F) and `vertex_ceiling.py`.

1. ⭐⭐⭐ **CODE AND DOC CLEANUP ON THE ONE SHIPPED POLICY (owner priority 1, 2026-09-09).** `transfer.py`'s
   `prepare` as named message builders and one lane class; the relay's leftover plumbing (15 unread
   `StepContext` fields, `RegionInit`'s unread precisions, two orphaned geometry functions); the
   `relay`-named labels renamed after a census; `test_transfer_policy.py` split by message; the sandbox
   pruned by the MOVE RULE and `DESIGN.md` converged on the shipped state. Every step gated bit-identical
   by `rename_identity.py --check` on the two frozen references.
2. ⭐⭐⭐ **IMPROVE THE MESSAGE POLICY (owner priority 2).** The debug loop on the worst in-scope rows of the
   benchmark page, one prototype arm at a time:
   `ISSUES: flux-price-witness-units`, `ISSUES: two-sided-exon-row`, `ISSUES: flux-floor-dispersion`,
   `ISSUES: ambig-node-as-a-gdna-source`; then the remaining both-stranded structures (`div`, the
   antisense's nascent variant, `nest`).
3. ⭐⭐ **THE gDNA LANDSCAPE PRIOR (owner priority 3).** `ISSUES: gdna-landscape-trains-on-false-positives`
   owns the zero-gDNA rows; the training population and the estimator at bound-only nodes first
   (`ISSUES: measured-prior-rung-4`, `ISSUES: landscape-trains-on-real-substrate`).
4. ⭐⭐ **THE POST-CALIBRATION, PRE-EM SETUP (owner priority 4).** `priors.py` / `result.py` / `derive.py`
   against `prior_vs_oracle.py` and the ruler column: `ISSUES: prior-fidelity-vs-deliverable`,
   `ISSUES: eb-shrinkage-magic-ess`, `ISSUES: g00-shrinkage-upstream-repair`, `ISSUES: capture-blind-gdna-divisor`,
   `ISSUES: per-transcript-prior-lane`, `ISSUES: u-ruler-arm`.
5. ⭐ **The calibration build thread** (part of priority 3): `ISSUES: measured-prior-rung-4` under the
   `ISSUES: reference-prior-refuted-at-concept-level` constraint, with
   `ISSUES: landscape-trains-on-real-substrate` as its payoff check.
6. ⭐ **The vertex atom, PRICED under the shipped policy (2026-09-09, `vertex_ceiling.py`, re-pointed to
   the two-phase solve and to the PARAMETER-vertex population)**: worth within 1 % on every stranded
   in-scope row, a few percent of the unstranded capture-OFF rows rising with gDNA, and a quarter to a
   third of the DEFERRED stratum. The in-scope value sits on silent genes' regions and nascent-free
   introns — the landscape prior's training population (rank 3) — and arrives through the refit prior,
   not at pass zero; `EQUATIONS.md` §9a.1/§9d.4 carry the derivation and the spike with no new constant.
   Re-derive with `vertex_ceiling.py --arm base|vertex_free --oracle-cache …` and `--compare`.

Then, in standing order: `ISSUES: performance-memory-bounded-solve` (owner: mandatory before 0.8.0) ·
the message-vs-prior question `ISSUES: refit-vs-message-arbitration` (the pre-EM entries above sit under
priority 4).

**Later / parked** (each has its entry): `expand-the-gdna-spectrum` · `psi-lambda-bracket-unshipped` ·
`alt-splice-rung-unverified` · `transfer-variance-premise` · `nascent-stress-sensitivity` ·
`f32-strand-tilt-at-half` · `hygiene-ledger` ·
`oracle-effective-length-diagnostic` · `flgap-panels-stale-nascent-model` · `rename-the-drain` · `drain-contaminates-certified-rna` (parked 2026-09-01: the ceiling refused the in-solve correction; two recorded follow-ups) · `the-cancelling-pair`
(refused twice) · `crossing-pool-contrast` (blocked) · `parked-capture-pilot-sign` ·
`pure-rna-mirror-asymmetry` · `capture-degeneracy-standing-risk`.

## ⛔ DELIBERATELY NOT NEXT

The length composition channel (retired until after 0.8.0) · anything whose only target is the DEFERRED
stratum · every mechanism in `ISSUES.md`'s **CLOSED / REFUSED** section — read it before proposing
anything, because each entry is a build that was measured and turned down, with the number that killed it.
