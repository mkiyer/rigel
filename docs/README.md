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

The gDNA-vs-RNA deconvolution stage: a bipartite region↔boundary **belief-propagation sweep** over the
accumulator payload. **Consolidated 2026-07-18 to a single master reference + a small set of anchors; the
per-phase roadmaps, message-layer WIP, and investigation notes are now under [`archive/`](calibration/archive/).**

**Read these — the current, converged set:**

| doc | what it is |
|---|---|
| [`CALIBRATION_MASTER.md`](calibration/CALIBRATION_MASTER.md) | **★ START HERE** — the single master reference: the goal, the pipeline in order, the two NPMLE roles (enrichment/transfer-variance vs the gDNA hyperprior), the AMBIG two-root problem, the two-solver reconciliation, what's implemented vs. what needs restructuring, and the phased plan P1–P5. |
| [`CALIBRATION_ARCHITECTURE.md`](calibration/CALIBRATION_ARCHITECTURE.md) | **★ AUTHORITATIVE INVARIANT** — the count-zero-information principle and the three composition sources. If any note disagrees about the *invariant*, this wins; `CALIBRATION_MASTER.md` is authoritative about the *design and plan*. |
| [`reference_prior_derivation.md`](calibration/reference_prior_derivation.md) | reviewer-facing — the **resolved** Beta(½,½) reference (§10) and the exact 1-D/2-D solver reconciliation (θ=arcsin(τ) cancels the tilt; class gap 0). |
| [`background_reference_derivation.md`](calibration/background_reference_derivation.md) | reviewer-facing — the aggregate background DNA level `ρ_bg` (the data-gated global hyperprior). |
| [`transfer_variance_formal_derivation.md`](calibration/transfer_variance_formal_derivation.md) | reviewer-facing — `σ²_transfer` message precision via projection onto the enrichment NPMLE (the NPMLE's Role A). |
| [`density_imputation_precision.md`](calibration/density_imputation_precision.md) | reviewer-facing (**derivation + proposal, for prototyping/review**) — the message precision for imputing a **density** across nodes that differ in total density, effective length, and composition: the 3-term variance (per-component count `1/ρE` + enrichment gap `σ²_transfer` + coverage non-uniformity). |
| [`npmle_crush_dissection.md`](calibration/npmle_crush_dissection.md) | **node-level dissection** of the Phase-2 refit's unstranded+capture `f_g` collapse: one captured exon (21,860 gDNA fragments, true f_g 0.902) traced init → pass-0 → the exact NPMLE `logprior` arithmetic → refit 0.0001. Root cause = an enrichment-blind depleted-mode prior; fix = the additive KDE (`gdna_kde_restore_plan.md`). |
| [`enrichment_sensitivity_worklog.md`](calibration/enrichment_sensitivity_worklog.md) | **★★ detailed work log & discoveries** — the three workflows (landscape-explore, simplify, enriched-sensitivity) and what each found; the KEY CONCEPT (pass-0's gDNA under-call is **directional & bracketed** `obs ≤ true ≤ total`); the elegant core (unified landscape + asymmetric projection); and **the pivot**: the under-call is a pass-0 *solver* problem on unstranded data (measured: enriched single-strand under-call −0.13…−0.30 unstranded vs −0.01 stranded), so §6 lays out the dissection plan + root-cause hypotheses. §8 = the **dissection results**: the identifiability floor (strand-flat + capture-starved crossings), single-exon vs multi-exon, and the **TSS/TES pure-gDNA message we corrupt into the wrong sign**. |
| [`message_system_derivation.md`](calibration/message_system_derivation.md) | **★★ THE message-system review + derivation + plan (2026-07-21)** — verifies the belief stores `(λ,θ)` + two precisions; states what the current message system does (composition/density hybrid + the τ-emission gate) with `bp_solver` line refs; documents the **verified τ-gag bug** (52% of spliced junctions silenced); derives the **per-component `(density, precision, gate)` message**, the equal-gate composition shift (recipient-side), and the **unequal-gate gDNA-only lower bound** (the intergenic↔exon case — the source-physics origin of the asymmetric-upward projection); reconciles the gDNA prior as the *designed* resolver; gives the **phased implementation plan** + the **open issues & questions**. Formal companion to [`intergenic_boundary_behavior.md`](calibration/intergenic_boundary_behavior.md) (the owner's mental model). |
| [`dna_prior_session_resume.md`](calibration/dna_prior_session_resume.md) | **★★ SESSION RESUME GUIDE — start here to resume the DNA-prior work.** Starter prompt, the current elegant best (unified landscape + asymmetric-upward projection, `enr_recovery` +0.25), the full arc/reframes, the pure-numpy exploration infrastructure (`scripts/debug/gdna_explore_lib.py` + cached substrate), the artifact index, the open next steps, the owner principles, and the abandoned production gates to clean up. |
| [`enriched_mode_sensitivity_hypotheses.md`](calibration/enriched_mode_sensitivity_hypotheses.md) | **★ hypothesis-generation + benchmarking plan** for making the DNA prior sensitively DETECT the enriched mode, measured by the **projection** endpoint (observed density → anchor μ* vs oracle). Projection function is validated (self-consistency ±0.06; a real enriched mode lets it pull under-called nodes up). Organizes the levers into families — projection read-out (mode/MAP/kernel), enriched-mode building (boundaries, two-population, gentler weighting, adaptive bandwidth, spliced-resolved), and evidence signals (enrichment/spliced/intron) to break the variance-identifiability — plus the sensitivity-first benchmarking plan and 4 workflow families. |
| [`gdna_hyperprior_plan.md`](calibration/gdna_hyperprior_plan.md) | **★★ THE plan (Role B, authoritative)** — the prior is a **third line of defense** (after strand, message-propagation): a **gentle anchor**, not a solver. The core is the **projection** = a *sampling-likelihood* onto the DNA landscape (`DensityNPMLE.project`): given a node's observed DNA density `d`, `μ* = Σ r_j μ_j`, `r_j ∝ w_j·N(d;μ_j,h)` — "which DNA level was this **drawn from**?" A far mode (however tall) has ≈0 responsibility; a nearby hill wins; valleys give wide `v*` ⇒ weak anchor. Enters ψ as a **gentle Gaussian pull toward μ\*** on the retained symmetric reference, `p = w_anchor/v*` (firm in a mode, weak in a valley). Fit additive Poisson-lognormal on the **confidently-solved** nodes (non-circular — the invariants: intergenic/intron/spliced). **Discards** the δ-pin, the Stage-0 floor, the enrichment log-shift, and `π_r` — all solved the wrong task (making the prior *decide* `f_g`). §3 is the exact implementation (reuse `project`, one `_gdna_arm` edit, both solvers, gated). |
| [`gdna_hyperprior_from_scratch.md`](calibration/gdna_hyperprior_from_scratch.md) | **★ THE build doc (Role B only)** — clean-slate, simplify-first rebuild of the gDNA composition hyperprior organized around one principle (**gravity**: pull to the nearest mode, harder when closer, able to pull UP). The 3-stage ladder (Stage 0 smooth Jeffreys reference floor over the full range → Stage 1 counting-precision bandwidth → Stage 2 enrichment-conditioning `d_g=ρ_g/e(x)`), the 8 open questions answered, the reviewer's 4 fixes assessed, what to STRIP, the two-pole acceptance test, and code touch-points. Role A (transfer-variance NPMLE) is explicitly OUT of scope / untouched. |
| [`gdna_crush_dissection_node1055.md`](calibration/gdna_crush_dissection_node1055.md) | **★ reviewer-ready — proves the crush is NOT a solver bug.** The 5-hop δ-pin code path (fit → project → build → apply → read-out, with the exact snippets + file:line) and the decisive swap-the-prior experiment: hand the solver an enriched-mode prior and node 1055 lands at f_g=0.863 (end-to-end, Part C). Root cause = the *fitted* prior has no enriched mode, from **circularity** (trained on pass-0's own under-call), **minority suppression** (competitive EM), and the **Jeffreys replacement** (worse-than-pass-0). Escape = enrichment-conditioning. Runs via `scripts/debug/crush_dissection.py`. |
| [`gdna_hyperprior_clean_slate.md`](calibration/gdna_hyperprior_clean_slate.md) | **★ START HERE for the gDNA hyperprior (Role B)** — clean-slate study: the enrichment/depletion requirement (the δ-pin), the released v0.7.1 KDE + the epsilon floor, KDE-vs-NPMLE, the two crush mechanisms (valley-collapse vs silenced-enriched-mode), the training subset, precision as **weight vs width** ("one count, spread by precision"), bandwidth options, the levers/experiments, alternatives, and a reviewer's code-module map. Includes a new-session starter prompt. |
| [`cliff_message_derivation.md`](calibration/cliff_message_derivation.md) | reviewer-facing — the cross-node message **mode** across the enrichment cliff (the eff-length-frame log-odds shift); §9 records the rejected "hybrid enrichment-corrected density" unifiers and why the exon floor is identifiability, not the mode. |
| [`oracle_and_benchmarking.md`](calibration/oracle_and_benchmarking.md) | the one trusted debug/benchmark path — the accumulator-consistent oracle (`scripts/debug/oracle.py`) and the net-fragment-flow metric. |
| [`archive/`](calibration/archive/) | superseded roadmaps/plans, the message-layer WIP series (relay derivation, mature-gate, spliced channel), and per-session investigation notes — **provenance, not current**. Durable findings are folded into `CALIBRATION_MASTER.md` §5–§9. |

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
