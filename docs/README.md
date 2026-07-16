# Rigel documentation

Rigel is a Bayesian RNA-seq quantifier that jointly models mRNA, nascent RNA (nRNA), and genomic
DNA (gDNA) contamination. This is the documentation map. The **code is the source of truth**; these
docs describe the current system (v0.7.0).

## Start here

| doc | what it is |
|---|---|
| [`MANUAL.md`](MANUAL.md) | user manual — install, the `index` / `quant` / `sim` / `export` CLI, options, outputs |
| [`METHODS.md`](METHODS.md) | the scientific method (the mRNA/nRNA/gDNA model + the three pipeline stages) |
| [`ARCHITECTURE.md`](ARCHITECTURE.md) | the code map and data flow (scan → calibrate → quant; Python modules + C++ extensions) |
| [`parameters.md`](parameters.md) | configuration & CLI parameter reference |
| [`BENCHMARKING.md`](BENCHMARKING.md) | the synthetic benchmark suite + net-fragment-flow evaluation |
| [`SIMULATOR.md`](SIMULATOR.md) | the synthetic read simulator (`rigel sim`) |
| [`PUBLISHING.md`](PUBLISHING.md) | building + publishing `rigel-rnaseq` to PyPI |
| [`TODO.md`](TODO.md) | open items |

## Calibration (`calibration/`)

The gDNA-vs-RNA deconvolution stage: a bipartite region↔boundary **belief-propagation sweep** (single
forward-backward pass) over the accumulator payload.

**Anchors** — read these first:

| doc | what it is |
|---|---|
| [`CALIBRATION_ARCHITECTURE.md`](calibration/CALIBRATION_ARCHITECTURE.md) | **★ AUTHORITATIVE** — the information & precision model: the count-zero-information principle and the three composition sources (strand likelihood, cross-node imputation, global gDNA prior). Read this first; if any other note disagrees, this wins. |
| [`calibration_roadmap.md`](calibration/calibration_roadmap.md) | **★ START HERE** — the living roadmap (partially stale — its Phase 1 reference work has landed, and a fresh agreed roadmap is planned). Where calibration stands, the derived Jeffreys reference (`c=½`), why "prior-free" does not exist, why capture breaks the prior, the traced failing nodes, the ordered phases. |

**Current initiative (message layer):**

| doc | what it is |
|---|---|
| [`mature_crossing_gate.md`](calibration/mature_crossing_gate.md) | ⚠️ **DISMANTLED (2026-07-16, commit `5e54fdc5`)** — landed, then removed for a pristine gate-free relay (`message_model_derivation.md` §5). Retained as the diagnosis/provenance of the exon→intron mature leak (its §1 defect is still authoritative); the gate itself is retired. Helped 7/7 cached conds (ALL mwae 0.169→0.193), so the dismantle regresses the mature-heavy suite until `σ²_transfer` + the nascent factory land. The `mrna_active_*` mask + D4 strand-aware splice routing stay. |
| [`nascent_rna_sourcing_regression.md`](calibration/nascent_rna_sourcing_regression.md) | **KNOWN/ACCEPTED regression** from the gate — unstranded AMBIG introns lose RNA support and over-call gDNA. Deferred to a dedicated "nascent factory" session (density-discrepancy sourcing). Owner-accepted for now. |
| [`dof_pie_relay_derivation.md`](calibration/dof_pie_relay_derivation.md) | **➡ CURRENT — the derivation (2026-07-16), third-party REVIEWED (verdict: proceed).** The rigorous, reviewer-facing theory for relaying the free coordinate `(λ[,θ])`: the coherent EP relay on the composition manifold, the two-level model (stored coordinates / transmitted per-component density), the fold-as-chain-rule, the count-zero-info count-term `1/M_src`, and the honest `σ²_transfer=0` liabilities. Validated by `scripts/debug/dof_pie_relay_check.py` (C1–C9). |
| [`dof_pie_relay_implementation_plan.md`](calibration/dof_pie_relay_implementation_plan.md) | **✅ S1–S4 LANDED** — the coherent `(λ,θ)` relay is in (`n_src>M` 424–925/cond → 0; relay pie = 1.000 everywhere). §7 has the live PROGRESS + the (A) accuracy A/B (flagship unstranded+capture leak −21…−25%; stranded regresses) + the (B) diagnostics. Two behavioural test failures shown to be artifacts. |
| [`message_model_derivation.md`](calibration/message_model_derivation.md) | **★ CURRENT — the formal message model + the session-close decisions.** A message is NOT a pie: it's per-component density factors (BP edge potential factorizes); node state = one DOF, message = free to speak about one component. The residual errors are **precision** (`σ²_transfer=0` overrules correct self-solves — the load-bearing gap) and the **mature gate** (a band-aid — now **✅ DISMANTLED**, commit `5e54fdc5`). Next: build `σ²_transfer` (the load-bearing fix) → the nascent factor (`RNA − mature`). |
| [`dof_pie_model_fix.md`](calibration/dof_pie_model_fix.md) | the **measurement + stepwise plan** (item 2) for the above. Its §1 defect is re-derived POST-GATE (62–70% of solvable nodes incoherent, 424–925/cond with a component fraction >1, MAX 600× under capture; `scripts/debug/pie_probe.py`). The case stands independent of item 1. |
| [`boundary_spliced_channel_design.md`](calibration/boundary_spliced_channel_design.md) | **priority #3 (deferred)** — "mature is a MEASUREMENT, nascent is an IMPUTATION": the spliced count (the only strand-free evidence) is delivered at the nascent guess's precision (176× haircut). Un-gagging it is what fixes the nascent regression. |

**Reference & prior:**

| doc | what it is |
|---|---|
| [`reference_prior_derivation.md`](calibration/reference_prior_derivation.md) | the **RESOLVED** `c=½` Jeffreys reference derivation (§10); the `(λ,θ)` DOF parameterisation the pie model uses; a landed prereq. |
| [`prior_ramp_and_bp_roadmap.md`](calibration/prior_ramp_and_bp_roadmap.md) | the `ψ_prior = kde + 0.5·λ` ramp analysis + the corrected bare `ψ = strand + logP_g + logP_r`; the strand Fisher information `∝ n(½−κ)²`. |
| [`npmle_roadmap.md`](calibration/npmle_roadmap.md) | the gDNA prior **live in code on this branch** (`GdnaRatePrior` in `gdna_rate_prior.py`, fit in `calibrate.py`). |
| [`calibration_prior_production_reference.md`](calibration/calibration_prior_production_reference.md) | ⚠ **currency question** — documents the M2-background + KDE prior "on `main`", but that code (`gdna_density_prior.py`) is **deleted on this branch** in favour of the NPMLE (`npmle_roadmap.md`). Scope-check against `main` before trusting. |

**Theory, tooling, and pending-decision docs:**

| doc | what it is |
|---|---|
| [`oracle_and_benchmarking.md`](calibration/oracle_and_benchmarking.md) | the one trusted debug/benchmark path — the accumulator-consistent oracle (`scripts/debug/oracle.py`) and the net-fragment-flow metric. |
| [`calibration_bp_theory.md`](calibration/calibration_bp_theory.md) | BP theory + sandbox proofs (precision currency, DOF, KKT self-defense). *Pending: keep as the standalone theory home, or fold into the new roadmap.* |
| [`node_prior_design.md`](calibration/node_prior_design.md) | the "is mature RNA present?" solver-prior classifier design. *Pending: still-implemented vs folded into `reference_prior_derivation.md`.* |
| [`calibration_selfsolve_status.md`](calibration/calibration_selfsolve_status.md) | the earlier self-solve initiative (its "boundaries gagged" thesis is refuted; §1.5 unified-solver plan still cited). *Pending: absorb §1.5/§1.3 into the new roadmap, then archive.* |
| [`archive/`](calibration/archive/) | superseded design/plan/diagnosis docs (incl. `splice_junction_absorption_fix`, `message_precision_and_dof`, `calibration_design`, `calibration_initialization`, `npmle_struggles`, `gdna_prior_zero_handling`, `mynotes`) — provenance, **not** current. |

## Accumulator (`accumulator/`)

The C++ fractional accumulator (the calibration sufficient statistics).

| doc | what it is |
|---|---|
| [`00_design.md`](accumulator/00_design.md) | **the spec** — per-region/boundary mass + FL pools; matched byte-for-byte by `tests/native/_accumulator_reference.py` |
| `01_implementation.md`, `01_junction_strand.md`, `accumulator_fragment_span_redesign.md`, `audit_phase1.md` | historical build/design/audit records |

## Simulator (`sim/`)

Design notes behind the `rigel sim` simulator (usage lives in [`SIMULATOR.md`](SIMULATOR.md)).

| doc | what it is |
|---|---|
| [`sim/architecture_redesign.md`](sim/architecture_redesign.md) | simulator architecture |
| [`sim/hybrid_capture_simulation.md`](sim/hybrid_capture_simulation.md) | hybrid-capture modeling |
| [`sim/suite_driver_unification_plan.md`](sim/suite_driver_unification_plan.md) | suite-driver design |

## Performance

| doc | what it is |
|---|---|
| [`performance/README.md`](performance/README.md) | memory & compute optimization — a **historical snapshot (v0.6.2)**: profiling infrastructure + the ranked roadmap |

---

**Note on history.** The `docs/` tree has been pruned repeatedly of superseded planning documents (the
long calibration-evolution journey now lives under [`calibration/archive/`](calibration/archive/)). That
history remains in git. The docs above describe the **current** system;
[`calibration/CALIBRATION_ARCHITECTURE.md`](calibration/CALIBRATION_ARCHITECTURE.md) is the authoritative
reference when code and any older note disagree.
