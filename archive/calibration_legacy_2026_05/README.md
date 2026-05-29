# Calibration v5 — Archived 2026-05-29

Salvageable algorithms preserved during the Phase A burndown documented in
[`../../docs/acc_caljointmodel/00_implementation_plan.md`](../../docs/acc_caljointmodel/00_implementation_plan.md).

## Contents

- `src/rigel/calibration/fl.py` — public FL surface (build_fl_models, Quality,
  per-pool ESS constants). Algorithm reference for the downstream EM scorer's
  per-fragment FL handling — not for re-use in the new calibrator (FL is no
  longer a calibration channel; see
  [`../../docs/caljointmodel/00_overview.md`](../../docs/caljointmodel/00_overview.md) §1).
- `src/rigel/calibration/_fl_sources.py` — counts-extraction helpers
  feeding `fl.py`.

## Why archived (vs deleted)

These modules contain working FL-mixture EM and Empirical-Bayes Dirichlet
smoothing logic. The downstream EM scorer may want to reuse these patterns
when its own FL handling is revisited.

## Do not resurrect into the calibrator

The new calibrator
([`../../docs/caljointmodel/03_inference.md`](../../docs/caljointmodel/03_inference.md))
deliberately excludes FL as a likelihood signal. The decision is recorded
in [`../../docs/accumulator/audit_phase1.md`](../../docs/accumulator/audit_phase1.md)
decision #6. Re-introducing FL into calibration without re-opening that
decision would regress the substrate-memory savings and the calibrator's
goal-directedness.
