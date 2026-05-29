# `docs/acc_caljointmodel/` — Consolidated Accumulator + Calibrator Rebuild

This folder holds the **active** implementation plan for the joint
fractional-accumulator + calibration-v6 rebuild. It supersedes the phase
plans in [`../accumulator/`](../accumulator/) and
[`../caljointmodel/`](../caljointmodel/) but **does not** replace the
design specs in those folders.

The accumulator substrate and the Phase-A burn skeleton are already
landed; the remaining work is organized as a **series of PRs**, each with
its own implementation doc that we critique and implement individually.

## Read in this order

1. **[`00_implementation_plan.md`](00_implementation_plan.md)** — the
   master plan (rev. 3): ground-truth tree state, resolved decisions
   D1–D6, the **Recovered-Code Map** (what to harvest from pre-burn git),
   and the **9-PR series** (PR 0 burn → … → PR 8 cleanup).
2. **[`prs/`](prs/)** — per-PR implementation docs, written one ahead.
   Start with **[`prs/PR00_burn_residue.md`](prs/PR00_burn_residue.md)**.

## Authoritative specs (not duplicated here)

| Topic | Doc |
|---|---|
| Fractional accumulator design | [`../accumulator/00_design.md`](../accumulator/00_design.md) |
| Calibration goals / scope / principles | [`../caljointmodel/00_overview.md`](../caljointmodel/00_overview.md) |
| Generative model | [`../caljointmodel/01_generative_model.md`](../caljointmodel/01_generative_model.md) |
| Failure audit of legacy calibration | [`../caljointmodel/02_failure_audit.md`](../caljointmodel/02_failure_audit.md) |
| Inference algorithm | [`../caljointmodel/03_inference.md`](../caljointmodel/03_inference.md) |
| Public API + locus prior pseudocount formulas | [`../caljointmodel/04_interface_contract.md`](../caljointmodel/04_interface_contract.md) |
| Validation plan | [`../caljointmodel/06_validation_plan.md`](../caljointmodel/06_validation_plan.md) |

> The master plan §1/§3/§4 flag where these specs diverge from the shipped
> code (substrate schema, boundary "half-split", `kappa_rna` source); those
> docs are reconciled in PR 8.

## What lives here later

| File | Created in | Purpose |
|---|---|---|
| `prs/PR01_…` … `prs/PR08_…` | each PR (one ahead) | Per-PR implementation plans |
| `validation_report.md` | PR 7 | Results vs acceptance criteria (val plan §7) |
| `postmortem.md` | PR 8 | What worked, what didn't, lessons |
