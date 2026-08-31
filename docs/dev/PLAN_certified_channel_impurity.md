# PLAN — the certified channel carries its own measured impurity

    ⚠ **A DEV DOC and a SANDBOX.** Nothing here is authoritative and nothing may cite into it.
    ⛔ **NOTHING HERE IS BUILT.** The DERIVE measurements were taken 2026-08-31 by prototypes in the
    session scratchpad (`leak_derive.py`); the defect record is `ISSUES: drain-contaminates-certified-rna`.

## 1. What the DERIVE established (all gated on byte-identity with production's own drain choices)

1. **The leak is EXACTLY posterior sampling, and nothing else.** Realized gDNA-in-spliced (A) matches
   the sum of the drain's own per-record spliced posteriors over true-gDNA records (B) within the
   draw's binomial noise on every condition measured: `g98 ss.99 ON` 499 vs 499.0 ± 15.7 (0.0σ);
   `g50 ss.99 OFF` 75 vs 79.6 ± 3.4 (1.4σ); `g05 ss.99 ON` 133 vs 135.9 ± 4.1 (0.7σ).
   ⛔ So "make the drain smarter" has no headroom at fixed information, and any deterministic rule
   (MAP, an odds floor) trades true splices for gDNA at the calibrated posterior's own rate.
2. **No structural channel for gDNA.** Records with NO genomic survivor (every hypothesis spliced —
   a forced leak if gDNA) exist in bulk but are 100 % mature on all three conditions (3,293 / 66,685 /
   157,916 records, gdna 0 in each): a gDNA fragment's genomic span is its length and survives the
   length limit on these panels.
3. **The regimes differ, and the information to say so exists at drain time.** At `g98 ss.99 ON`
   the leakers' own P(genomic) sits at median 0.494 (the model KNEW it was torn); at low-gDNA rows
   the few held gDNA records lean spliced (RNA-dominated ρ) but are rare (188–277 records).
4. **The harm is concentrated and qualitative**: at `g98 ss.99 ON`, 982 leaked boundary-spliced
   deposits land on 382 boundaries, and **137 of those carry ZERO true spliced RNA** — a false
   certified-RNA claim at an object with none (291 deposits; the p90 leaked object's certified count
   is 100 % gDNA). The sj axis is nearly clean (1 such object). At `g50 ss.99 OFF`: 36 objects.
5. **Two candidate repairs are REFUTED by these numbers.** A provenance split (drained records leave
   the certified channel) evicts 134,850 correct records to remove 75 contaminants at `g50 ss.99 OFF`
   — and those are the systematically-LONG fragments whose re-inclusion is what repaired the RNA
   length law (`ISSUES: the-rna-length-law-fix`, CLOSED). A posterior-odds floor is dominated by
   expected-value accounting because the posterior is calibrated (point 1) — and any fixed floor is a
   tuned constant besides.

## 2. The design, forced by the measurements

**The deposits do not change. The CERTAINTY CLAIM does.** At drain time, for every gap-spliced
deposit, accumulate onto the receiving certified objects (boundary-spliced bank and sj bank):

    E_genomic[obj]   += q_null(record)              the record's own P(genomic path)
    V_genomic[obj]   += q_null · (1 − q_null)       its variance — free from the same number

Both are already computed by `score_held_fragments`; the drain draws from exactly this posterior, so
the accumulation adds no model and no constant. `E_genomic` estimates the expected GENOMIC-TRUTH
contamination of each certified count — unbiased when the posterior is calibrated, which point 1
verifies in aggregate.

**The layering does the origin split.** `q_null` is a PATH statement (genomic vs spliced) — the
drain's own currency. What the certified-RNA consumers need is the gDNA part of it, and the genomic
path's composition is exactly what CALIBRATION owns: contamination_gdna[obj] ≈ E_genomic[obj] ×
f_g(local unspliced composition). So the drain records the path-level pair and calibration applies
its own composition — each layer speaks its own estimand, and the recurring one-name-two-populations
trap is dodged by construction. (At `g98 ss.99 ON` the genomic-truth contamination is 499 gdna +
1 nascent, so the f_g factor is nearly 1 exactly where the harm is.)

**Consumer form (to be decided in DESIGN, then A/B'd one at a time):**
* the certified counts read by `region_geometry`/the anchor become
  `count − E_genomic·f_g` with precision widened by `V_genomic` — a mean-and-variance correction,
  no threshold; or
* precision-only (widen, do not shift) — the conservative half, if the mean shift measures badly.

**Where the numbers live**: the drain is Python (`second_pass.drain`), so the two per-object float64
arrays need NO accumulator/ABI change — they can ride the drained payload the way `DrainQC` does.
⚠ The scan cache stores pass one, so nothing is cached; instruments get them through the drain they
already run (`DESIGN.md` §4.3).

## 3. Sequencing

0. ⭐ Prototype the accumulation outside `src/` (extend the session's `leak_derive.py` mapping
   records → receiving objects) and verify per-object: `E_genomic` vs the oracle's realized leak,
   condition by condition — the per-object calibration check that aggregate point 1 already passed.
1. Prototype the CONSUMER correction on the certified banks inside a `calibrate` injection arm and
   score on `calibration_vs_oracle.py` per stratum, both zero controls, plus the fl-gap sign arms —
   ⛔ the RNA-length-law lesson: the drain's repair of the spliced pool must SURVIVE (the corrected
   channel must not resurrect the −4 bp bias), so `fl_pool_purity`-style means are checked beside
   the metric.
2. One consumer at a time (`TRAPS: one-thing-varied`): the anchor's certified flux and ψ's certified
   terms are separate arms.
3. Only then `src/`: the drain-side accumulation + the consumer wiring, each with a falsification
   test first (inject a known q pattern, assert the arrays; perturb, watch it fire).

## 4. Open

* ⚠ The mirror loss (true spliced drawn genomic — certified-channel FALSE NEGATIVES) is measurable
  with the same machinery (`Σ (1−q_null)` over drawn-genomic) and unpriced; the owner's stated
  priority is the false-positive side.
* ⚠ Mature records whose TRUE path is the genomic gap (another isoform's intron in the gap) make the
  full genomic-truth denominator larger than gdna+nascent; the per-object check in step 0 should
  count them via the mrna partition's drained banks rather than assume.
* ⚠ Whether the correction should also reach the DEFERRED stratum's consumers unchanged — it should
  (it is not a stratum-specific mechanism), but the A/B reports it per stratum regardless.
