# HANDOFF — ψ's composition closes; the ranked list is `ROADMAP.md` §1 and rank 0 is half-done

    ⚠ **A DEV DOC. Nothing may cite it.** ⛔ It is a HANDOFF, not the state. The state is `ROADMAP.md` §0,
    the ranked list is `ROADMAP.md` §1, the rulings are `DESIGN.md` §6b.1 / §6c, the derivations are
    `EQUATIONS.md` §9c.1 / §9c.2, and the lessons are `TRAPS.md`. Everything settled has already MOVED.

    ⭐ Session of 2026-08-17. Branch `fragment-length-gold-standard`. Target: **0.8.0**.

## ⭐⭐⭐ START HERE — where the tool is

Two things landed and are COMMITTED. Suite **0 failed / 3,496 passed / 0 skipped / 11 xfail**, lint clean.

1. **ψ's reference has a MEAN set from the annotation** (`CalibrationConfig.structural_reference`, default
   **ON**). Panel per stratum 0.930 / 0.908 / 0.925, `g00` control bit-identical.
2. **ψ's composition CLOSES, structurally** — 100.00 % of published objects, against 74.7 % / 77.2 %
   before. Panel 0.999 / 0.995 / 0.996, control 0.995.

⛔ Both are *correctness* wins sitting inside the panel's 0.996–1.013 noise floor. **Neither moved the
0.8.0 metric materially.** The items that would are `ROADMAP.md` §1 ranks 0–4, and rank 0 is half-run.

## ⛔⛔ WHAT TO DO FIRST — `ROADMAP.md` §1 RANK 1, THE ZERO-gDNA FALSE POSITIVE

⛔ **AND THE RANKING WAS CORRECTED BY THE OWNER, WHICH IS WORTH RECORDING.** This doc first ranked the
post-calibration ruler items on top. The owner's stance is to move through the pipeline IN ORDER — scan →
accumulator → **calibration** → post-calibration → EM → output — and `CLAUDE.md`'s own 0.8.0 scope agrees:
*"THE FOCUS IS CALIBRATION, and the metric is the calibration result against ORACLE CALIBRATION — not the
end-to-end transcript number."* The ruler and the per-transcript prior are post-calibration; they move the
THERMOMETER. Tuning them against a calibration with a 17.68 % hole would be tuning against a broken input.

**MEASURED 2026-08-17, all four `g00` conditions through `calibrate`.** Truth is `f_g = 0.0000` at every
object, so every gDNA fragment is a false positive with nothing to cancel it::

    unstranded x capture OFF   1,898,257   17.68 %   <- IN SCOPE, and 63x the stranded rate
    unstranded x capture ON      350,623    3.50 %
    stranded   x capture OFF      29,692    0.28 %
    stranded   x capture ON       21,348    0.21 %

    REGION intergenic      0 (0.0 %)          REGION exon   1,164,526 (61.3 %, aggregate f_g 0.2541)
    REGION intron          0 (0.0 %)          BOUNDARY        733,731 (38.7 %, 0.1192)

⭐ It is LIVE in the shipped configuration — unlike the relay items, which sit behind
`message_propagation = False` — and it is perfectly localised: the two classes the structural reference
pins contribute NOTHING, and all of it sits where mature RNA can be.

⚠ The oracle-effective-length diagnostic (§1 rank 0) is still half-run and still cheap; it is rank 0 only
because it re-ranks two POST-calibration items and costs half an hour, not because it outranks this.

## ⚠ THE HONEST STATE OF THE TWO WINS

* The panel cannot resolve either at 0.5 %. Both are justified on THEORY plus a non-regression, not on an
  accuracy gain, and `ROADMAP.md` says so.
* `structural_reference`'s central claim — nascent RNA ≈ 0 — **cannot be falsified by this panel**, which
  holds `nrna = 0` on all 16 rows. The two mechanisms that CAN refute it are built and measured live
  (`τ_fac = 161.4` at an intron, unstranded). An `nrna > 0` rung is low priority per the owner but it is
  the only way to score that claim.

## ⭐⭐ WHAT THIS SESSION COST, SO IT IS NOT REPEATED

Three defects were introduced and caught, all by gates or by adversarial hunting, none by inspection:

* deriving the RNA total from `f_g` made it inherit `f_g`'s grid SNAP, where the retired posterior-mean
  read-out had been exact (`test_relay_mass_pin` caught it) — `TRAPS: deriving-one-coordinate-propagates-its-error`;
* the ½-quantile was first interpolated in `f_g`-space, biasing a concentrated posterior toward ½ by
  2.7e-03 — `TRAPS: interpolate-on-the-axis-where-the-lattice-is-uniform`;
* `_compose` did not clamp the tilt share, so an out-of-range share gave a NEGATIVE fraction that still
  summed to 1 and passed every closure gate.

⛔ And a reporting failure worth more than any of them: **21 test failures were read off the last eleven
lines and reported as "21 golden failures". Two were real.** `TRAPS: read-the-whole-failure-list`.

## ⚠ OPEN, MEASURED, NOT FIXED

* **`ROADMAP.md` §1 rank 6 — the f32 cube manufactures a strand tilt at κ = ½.** Falsification already
  written: put the strand term in float64 and `w_pos → 0.500000000` at every depth. Negligible at panel
  scale (≤ 3.4e-4), pre-existing, one dtype to repair.
* **A mass-only slot** (spliced mass, zero unspliced counts) is dispatched and publishes `(0,0,0)`.
  Pre-existing and ungated; found while bug-hunting, not investigated.
* **The relay's composition licence** (§1 rank 2) knows transcript termini but not `mrna_active` flipping.
  Confirmed, localised to `lam_channel`, dormant because `message_propagation` is OFF.

## ⛔ THE DISCIPLINE THAT ACTUALLY FOUND THINGS

Every real defect this session came from a gate firing or from an adversarial pass, and **three separate
claims I had written into permanent docs turned out to be wrong when re-measured** (`f_g` bit-identical;
the flat-regime residual; the closure bound 1.11e-16 vs 2.22e-16). ⭐ Re-derive before quoting, and
re-measure after any change to `_posterior_median_fg` — it is on the critical path for every slot.
