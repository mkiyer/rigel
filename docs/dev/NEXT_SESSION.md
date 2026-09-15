# NEXT SESSION — start here (2026-09-14, after the ruler and the flux-witness session)

This file is only how to begin. The port is `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` §F; the θ thread's
derivations are `docs/dev/THETA_QUADRATURE.md`.

## The agreed order of the sessions (owner, 2026-09-14)

1. ~~Code review and cleanup~~ — DONE 2026-09-14.
2. ~~The test chromosome's new structures, both panels remeasured~~ — DONE 2026-09-14 (six commits after `1f885a7e`).
3. ~~The ruler at zero gDNA and the flux price's witness~~ — LANDED 2026-09-14, two commits PREPARED for your go
   (below).
4. **THIS ONE — the performance re-baseline and THE PORT** (`ISSUES: performance-memory-bounded-solve`,
   plan §F): two back-to-back profiler pairs at 8 threads on the deep library, identity references frozen
   fresh on that tree (they are fresh on this one: re-frozen after each mechanism, reasons in the commits),
   `sweep_replay.py capture` (`sweeps_MO_3021_step6` predates the location floor, the ruler and the witness),
   then the port of the block solve, gated bit-identical on every stage.

## What this session left, and how to commit it (owner's go)

The two commits are prepared, not made (the standing ruling: committed on your go): the working tree holds
both mechanisms and the session's docs, and each mechanism's file state is snapshotted with its message
under the session scratchpad's `commits/7_ruler` and `commits/8_flux` (the first six snapshots, already on
`main`, moved to `commits/committed/`). `commits/commit_series.sh` replays the two in order and leaves the
tree clean; `cat commits/*/MESSAGE.txt` and `git diff` are the review. Every gate is on record in the messages:
the falsification tests verified failing first, the perturbations fired, the suite re-derived, the three
`review_identity_*` references re-frozen after each mechanism with the reason stated and checked 3/3
BIT-IDENTICAL on this tree, both panels' instruments run on the landed tree.

### The two mechanisms

TWO MECHANISMS, one commit each, both derived, prototyped outside `src/` on the test chromosome, confirmed on
the ladder, then landed with a falsification test verified failing first and every gate watched firing:

* **The ruler at zero gDNA** (`DESIGN.md` §7.2, `EQUATIONS.md` §11; `ISSUES: g00-shrinkage-upstream-repair`
  and `u-ruler-arm` CLOSED). The reference the EM's effective length contracts against is the located
  enriched mode of the fitted gDNA landscape, on the result as `CalibrationResult.gdna_reference_density`;
  no enriched mode ⇒ no contraction, exactly. Both panels: the zero controls 0.154 / 0.141 → 1.000 with
  nothing moved (51,436 / 5,108 transcripts had moved); both in-scope capture-OFF strata P = O = 1.000
  (from P 0.957 / 0.970, O 0.923 / 0.926 on the ladder); stranded capture-ON P/O 1.013 / 1.015 against
  1.011; the solve untouched (`policy_benchmark.py` identical). The private kernel density and its two
  constants deleted; the U arm retired. 13 goldens moved — gDNA-free and capture-OFF toys whose old rulers
  had contracted transcripts 49× and 74× against false-positive specks; the antisense_contained golden's
  antisense count 0.0 → 107.4 against a truth near 104 (its 0 was the ruler, not the tilt atom).
* **The flux price's witness** (`DESIGN.md` §6b.13, `EQUATIONS.md` §12; `ISSUES: flux-price-witness-units`
  CLOSED). The exon's witness of the junction→exon price is its column count on the protocol's share of the
  opportunity, `(c_s, κ_read · a_r)`: whole-strand units on both sides, the pair's agreement priced as
  counting alone at every κ. Correct by derivation and neutral on both panels in scope (within 0.3 %), −10 %
  on the test chromosome's deferred stratum, the ten zero rows within a fragment, the golden's record exon 0.235 →
  0.008 gDNA. The cost the issue named (0.48 nats² on every unstranded flux level) was real in the price and
  worth ~0 on the metric: the flux ceilings do not bind in scope. Three forms refused with numbers.

## What the mechanisms uncovered — recorded, not fixed

* **`ISSUES: nested-antisense-leak-under-the-sane-ruler`** (EM-side, priority later with the per-transcript
  prior lane): with the EM's ruler honest, a gDNA-free library's strand-flipped intronic nascent fragments are
  assigned to a nested antisense transcript nothing witnesses — 24 of 2,000 at SS 0.9, 124 at SS 0.65 on the
  negative control of `tests/scenarios/test_antisense_intronic.py`; the old bounds held only because a
  fabricated reference had contracted the host's nascent entity 3.9×. Two strict xfails (5 in the suite now).
* **The flux ceilings do not bind in scope**: the 0.48 nats² the witness issue named cost nothing measurable
  on either panel's in-scope strata — the unstranded residual sits elsewhere (`policy_benchmark.py
  --by-class`: introns 45 %, `exon|intron` boundaries 14 % on the ladder's unstranded OFF).
* **The FIELD gate flake** on λ ≈ 7 boundaries (DEFERRED by the owner, 2026-09-14) and the certifier's
  one-panel-at-a-time cache rule are in `docs/dev/` of the previous handoff and `TESTING.md` §0a.

## Decisions on record (unchanged, carried)

* Float64 for the whole of ψ; ONE solver; ONE λ lattice at a dimensionless step with a stated guarantee; no θ
  lattice anywhere; the tilt's hypothesis space is {pure +, pure −, mixed}.
* The strand channel's liveness is a protocol decision on the spliced 2×2; gDNA enters nowhere.
* The landscape trains only where a solve locates a slot; the ruler reads its located enriched mode.
* A few high-quality instruments, kept current; no suite gate polices instruments; the source cites no doc.
* `EMConfig.warm_start`'s `prior` and `uniform` arms stay; the `gdna_none_` condition names stay.
* The message cache's on/off switch, the refit count and the scan's thread split are the owner's; parallelism
  waits for the port. CI runs on demand only.
* The certifier's FIELD gate flake on λ ≈ 7 boundaries is DEFERRED (owner, 2026-09-14); the shadow floor on the
  test chromosome's zero rows is the designed control.

## The session scratchpad (persists; nothing in the tree cites it)

`/private/tmp/claude-503/-Users-mkiyer-proj-rigel/23defdd7-4a51-4398-89bc-ac9b16e623dc/scratchpad/`:
`commits/7_ruler`, `commits/8_flux` (the two snapshots with `MESSAGE.txt`; `commits/commit_series.sh` replays every
`[0-9]_*` snapshot in order — the first six are already committed), `s3/` (the harnesses `ruler_arms.py`,
`flux_arms.py`, `landscape_probe.py`, `leak_diag.py`, `nested_diag.py`; the arms' outputs under `s3/cvo`, `s3/flux`;
the landed runs under `s3/final` (ruler) and `s3/final2` (both); the derivation notes `DERIVATION_ruler.md`,
`EQUATIONS_11.md`, `EQUATIONS_12.md`; the previous identity references in `s3/identity_prev`).
