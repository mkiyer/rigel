# NEXT SESSION — the coordinate thread ran to a verdict (2026-09-01)

    ⚠ **A DEV DOC, and it is a HANDOFF.** It says where things stand, not what to build — the ranked
    list is `ROADMAP.md`. MOVE anything that settles into the permanent docs and DELETE this file.

## What happened this session

`ISSUES: arcsine-magnitude-coordinate` was executed end to end — step ① the dissection, step ② the
derivation on paper for BOTH solvers, step ③ the prototype outside `src/` with a per-stratum A/B.
**`src/` is untouched.** Suite green at **3,745** (baseline 3,743 **+2** for the two `docs/dev/` files
added — the derivation and the preserved prototype, jargon only each, confirmed by bracket-matched
collection; `CLAUDE.md`'s standing-baseline line is updated).

**The verdict, in one line: the attribution is CONFIRMED, the headline premise is REFUTED, and the
thread needs re-aiming at a MODEL term.** The full record is `ISSUES: arcsine-magnitude-coordinate`
(entries ①–⑥ + the re-aim) with the derivation in `docs/dev/DERIVATION_arcsine_magnitude.md`.
⛔ Do not re-derive any of it; do not re-run the refuted claims.

The three things a next session most needs:

1. **The under-call scales as `1/√n` across four decades of depth** (obs/pred 0.63…1.21 against
   `w/4 = 1/(2√n)` on 34,207 solved exact-vertex slots). That is a posterior WIDTH, not a grid
   artifact. The solver is answering honestly; the error is the price of a CONTINUOUS prior with no
   atom at the vertex, scored against a truth that sits exactly on it.
2. **The coordinate's whole leverage is ~2.5 % of the defect**, bounded by K-invariance: the two
   coordinates each converge to their own limit while the gap between them holds at 6.5e-04 (the λ
   bracket's truncated DOMAIN — refining K cannot recover it) against a 0.026 shortfall.
3. **λ is FINER at the vertices than uniform-φ**, and the median cannot reach a vertex in ANY
   coordinate (`t = 1 − 0.5/m ≤ ½` for the outer bin). The issue's premise was backwards.

## Where the prototype is

⭐ **`docs/dev/arcsine_proto.py` is PRESERVED IN THE REPO** — the mechanism itself: three module
patches (`_logodds_grid` → the midpoint uniform-φ lattice, the two arms → no written reference,
`_posterior_median_fg` → interpolate in φ) plus **23 perturbed self-test gates, all passing**. Run it
directly: `python docs/dev/arcsine_proto.py`. `install("noop")` reinstalls the ORIGINALS through the
same mechanism and is byte-identical to production on every condition and both metrics — that is the
harness's own falsification, and it held throughout.

The A/B drivers were left in the session scratchpad and are thin by comparison (an afternoon to
rewrite against the preserved prototype): `arcsine_ab.py` (walk metric vs certified `slot_truth`),
`arcsine_release_ab.py` (0.8.0's release metric via `calibration_vs_oracle.measure_condition`),
`arcsine_kinvariance.py` (the acceptance test), `band_dissect.py` (truth-banded dissection),
`neighbour_information.py` (the no-solver neighbour ceiling — see the parallel section).

## The re-aim — proposed AND priced in the same session; it is ALSO not a win

The obvious next mechanism is **an ATOM at the vertex** ("this slot holds NO RNA" as a point mass —
a point-null/spike-and-slab that no continuous reference expresses in any coordinate), because the
source already names the socket and its price: `simplex_logodds._rna_arm` — nothing fits `logP_r`,
and the unfitted reference is a fixed **3.107 nats** repulsion at `f_g = 0.999` (a 22:1 handicap);
objects with true `f_g ≥ 0.999` carry **49–83 %** of all calibration error and read **0.13–0.23
below** the vertex. The target population is **half the ladder's live mass** (50.0 %, no solver).

⛔⛔ **But it was priced before recommending it, and the ceiling is roughly a WASH.**
`vertex_ceiling.py`, three conditions, oracle truth pinned and the whole chain RE-SOLVED (`noop` pins
0 and is byte-identical — falsification passing; `vertex_free` pins ~70–79 k objects):
`vertex_free` gives region Σ|err| **+18,887** and boundary **−20,208** — they cancel. The looser
`vertex_all` is clearly worse (+94,490 / −18,237). Handing a vertex object the exact answer changes
what it BROADCASTS and the relay over-propagates it.

⭐⭐ **What that leaves, and it is the session's real handoff**: the near-vertex error is (a) an
honest posterior width that scales as `1/√n`, and (b) NOT removable by local certainty either. Both
the coordinate and the local-prior repair are priced and neither wins. **That points the question at
the MESSAGE LAYER** — which is also where `ISSUES: message-value-for-blind-slots` (below) now has its
first numbers. ⚠ Owed before anyone builds an atom: the whole-ladder `vertex_ceiling` run and a
`calibration_vs_oracle` equivalent (the ceiling above is that instrument's pass-0-flavoured metric).

## Cautions that still stand

* ⛔ Nothing measured before 2026-08-31 is comparable to current numbers (the drained-frame ruling,
  `DESIGN.md` §4.3).
* ⛔ The A/B (whole ladder, all 16, noop byte-identical on every one) **splits cleanly by stratum**:
  stranded × OFF WINS and more as gDNA rises (1.008 / 0.939 / 0.898); unstranded × OFF LOSES every
  row (1.107 / 1.071 / 1.014); stranded × ON mild loss; 3 of 4 zero controls win. ⭐ The coordinate
  helps where the strand LIKELIHOOD is strong and hurts where it is flat (there the posterior IS the
  reference, so the grid's tail resolution moves the answer directly). That is a real trade, and
  losing all three in-scope unstranded rows is what disqualifies it — never pool
  (`TRAPS: never-pool-the-strata`), the ladder total would have hidden it.
* ⚠ **Rung B is derived but NOT measured**: φ-native message DELIVERY (map a claim's precision by
  `(dx/dphi)²` — derivation §5). Rung A delivers claims pointwise in λ's coordinate, which is the
  honest suspect for the F-stage losses, since the C-stage (message-free) improves nearly everywhere.
  It cannot exceed the 2.5 % bound.
* ⚠ `_posterior_median_fg`'s interpolation coordinate is load-bearing and easy to get wrong — the
  histogram edges must be midpoints in the coordinate whose measure is UNIFORM, or a one-hot posterior
  comes back biased. Production documents this for λ; the prototype re-derives it for φ.

## Parallel, cheap, no solver — unchanged and unaffected by the verdict

`ISSUES: message-value-for-blind-slots` on the anchored twin block: is the information a blind
unstranded/AMBIG slot needs actually PRESENT in its neighbours? It decides the message policy's
ceiling. ⚠ The owner's "the policy is built once, on the settled coordinate" ordering was premised on
the coordinate landing; with the thread re-aimed, that ordering is the owner's to re-rule.
