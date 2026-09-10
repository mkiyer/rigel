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
  Re-derived 2026-09-10: in scope the perfect prior no longer improves the transcript number at all
  (it is worse or equal on all three strata); it is worth a quarter of the deferred stratum's; and the
  `g00` rows carry the largest transcript error of any stratum under BOTH arms — the ruler above.
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
  minimal harm on stranded, never pooled** — `policy_benchmark.py --panel ladder`. Since the prior's
  two landings the zero rows sit at a few hundred fragments under BOTH policies, so the "beats silence"
  count is read on the contaminated rows, where every row favours `transfer`; `calibration_walk.py` says
  messages still carry the stranded capture-ON rows and are essential on the deferred stratum.
- **The gDNA landscape prior**: DONE for 0.8.0 (`DESIGN.md` §7.1, 2026-09-10) — a node whose only
  evidence is a bound does not train it, its location-free kernels are placed by the previous refit,
  its grid spans every region and boundary. The zero controls are solved (`calibration_vs_oracle.py`
  reads 0.0000 on both axes; `landscape_training_census.py` re-derives the population); in scope the
  per-object composition error is unchanged at about one percent.
- **ψ**: the composition closes structurally on every published object (`test_composition_closes.py`);
  the reference location is DELETED (owner, 2026-08-24; the surviving form is `DESIGN.md` §6b.1); the
  λ-bracket widening is built and ships OFF — `ISSUES: psi-lambda-bracket-unshipped`.
- **The prior assembler**: with perfect masses its own error is parts per thousand —
  `prior_vs_oracle.py`, `mass_prior_ab.py`.
- **The largest number on the 0.8.0 metric page is now the RULER, not the composition**: the `g00`
  effective-length shrinkage fabricates a reference from the few hundred residual fragments and
  contracts every transcript to 0.15 of its length (`ISSUES: g00-shrinkage-upstream-repair`, re-priced
  2026-09-10 — the fix is the detector); the never-passed per-transcript prior lane
  (`ISSUES: per-transcript-prior-lane`) is the other.
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

## ⭐⭐⭐ NEXT — the recommended order (audited 2026-09-10, after the landscape prior's two landings)

**The method stays the owner's dissection loop** (2026-08-20): run the panel → worst IN-SCOPE scenario →
rank its objects by error mass (`worst_objects.py`, `calibration_walk.py`) → find the mechanism → gated
fix → **add the offending transcripts to the test chromosome** → re-run → repeat.

⭐ **The facts this ranking leans on, each with the instrument that re-derives it** (every one re-run on
2026-09-10 under the two landings; the outputs are under
`~/Downloads/rigel_runs/arms/2026-09-10_post_landscape/`): the zero controls are solved on the 0.8.0
metric (`calibration_vs_oracle.py`, both axes 0.0000); every in-scope stratum's per-object composition
error is about one percent and did not move; the ruler's factor at `g00` is 0.15 of the oracle's 1.00 and
is the largest in-scope number on the metric page; end to end, a perfect prior is worth nothing in scope
and a quarter of the deferred stratum (`quant_accuracy.py`); the vertex ceiling is within 1 % on every
stranded row, 2–7 % of the unstranded capture-OFF rows rising with gDNA, a quarter to a third of the
deferred stratum (`vertex_ceiling.py`); by class the in-scope residual sits on the intron's own solve
(unstranded OFF) and on exon|exon boundaries and walled exons (stranded ON) (`policy_benchmark.py
--by-class`).

1. ⭐⭐⭐ **THE RULER at zero gDNA — `ISSUES: g00-shrinkage-upstream-repair` (rank 2's territory, now
   first).** A gDNA-free library is the modal real case, the composition there is now right, and the
   effective length the EM divides by is still 0.15 of the truth because the reference-density detector
   accepts any five slots with positive mass. Derive what "this library has an enriched gDNA mode" is
   evidence of (the landscape's enrichment detector, a boolean), prototype it outside `src/` in
   `_global_reference_density`'s caller, judge on `calibration_vs_oracle.py` ③ (the ruler's `factor`
   and `moved`) per stratum with both zero controls, then `quant_accuracy.py`. Settle first whether the
   instrument's "exactly 1.000 off capture" contract is stale: both P and O read 0.92–0.97 there.
2. ⭐⭐ **THE REST OF THE PRE-EM SETUP (owner priority 4)** — `priors.py` / `result.py` / `derive.py`
   against `prior_vs_oracle.py` (re-run it first) and the ruler column: `ISSUES: prior-fidelity-vs-deliverable`,
   `ISSUES: eb-shrinkage-magic-ess`, `ISSUES: capture-blind-gdna-divisor`, `ISSUES: per-transcript-prior-lane`,
   `ISSUES: u-ruler-arm`.
3. ⭐ **THE INTRON'S OWN SOLVE on unstranded capture-OFF** — 46 % of the in-scope error there is the
   intron class at 2.3 % of its own fragments (`policy_benchmark.py --by-class`), the factory profile's
   resolution against the intergenic background (`density_deconv`); dissect with `worst_objects.py`.
4. ⭐ **THE VERTEX ATOM** — priced unchanged (above); its in-scope value is 2–7 % of the unstranded
   capture-OFF rows at `g50`/`g98`, on silent genes and nascent-free introns; a mechanism for it is the
   prior's reference (`ISSUES: reference-prior-refuted-at-concept-level` constrains the form) or the
   intron's own solve, not a message.
5. **THE MESSAGE POLICY, only where a row is above the bar** (owner priority 2): one prototype arm at a time
   through `policy_prototype.py --module`, halves apart, pass zero beside the pipeline:
   `ISSUES: flux-price-witness-units`, `ISSUES: two-sided-exon-row`, `ISSUES: flux-floor-dispersion`,
   `ISSUES: ambig-node-as-a-gdna-source`, `ISSUES: message-layer-open-cases`.

Then, in standing order: `ISSUES: performance-memory-bounded-solve` (owner: mandatory before 0.8.0) ·
`ISSUES: refit-vs-message-arbitration` (re-read under the E-step: the walk now says the prior does the
unstranded rows and the messages the stranded capture-ON ones).

**Later / parked** (each has its entry): `expand-the-gdna-spectrum` · `psi-lambda-bracket-unshipped` ·
`transfer-variance-premise` · `nascent-stress-sensitivity` · `f32-strand-tilt-at-half` · `hygiene-ledger` ·
`oracle-effective-length-diagnostic` · `flgap-panels-stale-nascent-model` · `rename-the-drain` ·
`drain-contaminates-certified-rna` (parked 2026-09-01: the ceiling refused the in-solve correction; two
recorded follow-ups) · `the-cancelling-pair` (refused twice) · `crossing-pool-contrast` (blocked) ·
`parked-capture-pilot-sign` · `pure-rna-mirror-asymmetry` · `capture-degeneracy-standing-risk`.

## ⛔ DELIBERATELY NOT NEXT

The length composition channel (retired until after 0.8.0) · anything whose only target is the DEFERRED
stratum · every mechanism in `ISSUES.md`'s **CLOSED / REFUSED** section — read it before proposing
anything, because each entry is a build that was measured and turned down, with the number that killed it.
