# `docs/acc_caljointmodel/` — Consolidated Accumulator + Calibrator Rebuild

This folder holds the **active** implementation plan for the joint
fractional-accumulator + calibration-v6 rebuild. It supersedes the
phase plans in [`../accumulator/`](../accumulator/) and
[`../caljointmodel/`](../caljointmodel/) but **does not** replace
the design specs in those folders.

## Read in this order

1. **[`00_implementation_plan.md`](00_implementation_plan.md)** —
   the seven-phase rollout (A: burn → B: scanner wiring → C:
   scaffold → D: implementation → E: integrate → F: validate →
   G: cleanup). Includes file-touch matrix and acceptance gates.

## Authoritative specs (not duplicated here)

| Topic | Doc |
|---|---|
| Fractional accumulator design | [`../accumulator/00_design.md`](../accumulator/00_design.md) |
| Legacy accumulator audit | [`../accumulator/audit_phase1.md`](../accumulator/audit_phase1.md) |
| Calibration goals / scope / principles | [`../caljointmodel/00_overview.md`](../caljointmodel/00_overview.md) |
| Generative model | [`../caljointmodel/01_generative_model.md`](../caljointmodel/01_generative_model.md) |
| Failure audit of legacy calibration | [`../caljointmodel/02_failure_audit.md`](../caljointmodel/02_failure_audit.md) |
| Inference algorithm | [`../caljointmodel/03_inference.md`](../caljointmodel/03_inference.md) |
| Public API + locus prior pseudocount formulas | [`../caljointmodel/04_interface_contract.md`](../caljointmodel/04_interface_contract.md) |
| Validation plan | [`../caljointmodel/06_validation_plan.md`](../caljointmodel/06_validation_plan.md) |

## What lives here later

| File | Created in phase | Purpose |
|---|---|---|
| `validation_report.md` | F | Results vs acceptance criteria per [`../caljointmodel/06_validation_plan.md`](../caljointmodel/06_validation_plan.md) §7 |
| `postmortem.md` | G | What worked, what didn't, lessons |
