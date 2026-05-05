# M9 — Validation + Documentation (Implementation Plan)

**Status:** in progress.
**Parent plan:** [calibration_v6_plan.md](calibration_v6_plan.md) §M9.
**Predecessor:** M8a/M8b/M8c (commits `2fadb1b`, `f5454e6`, `5f4754f`).
**Cleanup commit:** `dca1788` (ruff clean, goldens regenerated, orphan v1 test deleted).

This plan covers the *final* milestone of the v6 calibration redesign.
It has three workstreams, executed roughly in parallel:

1. **Plan/code reconciliation** — close or document the gaps between
   the v6 plan as written and the code that actually shipped in M8.
2. **Validation matrix** — run rigel on the synthetic, hybrid-capture
   and Armis2 benchmarks called out in `calibration_v6_plan.md` §M9
   and confirm the design's claimed accuracy targets.
3. **User-facing documentation** — `docs/MANUAL.md`, `docs/METHODS.md`,
   `docs/parameters.md`, `CHANGELOG.md`, and a tagged release.

---

## 1. Plan/code reconciliation

The v6 plan was written before implementation; M7 and M8 deviated from
it in several reviewable but unrecorded ways.  M9.1 either closes the
gap or amends the plan to match what shipped.  This step is *cheap*
and *prerequisite* — running validation against an undocumented surface
makes the validation results un-auditable.

### 1.1 Inventory of deviations

| # | Plan §       | Plan called for                                                                | Shipped (HEAD = `dca1788`)                                                                            | Decision              |
|---|--------------|--------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|-----------------------|
| 1 | M7           | `compute_pool_fl_models` returning `PoolFLModels`                              | `build_fl_models` returning `FLModels`                                                                | **Document** — already justified in `m8_implementation_plan.md`. Update §M7 of v6 plan to reference the shipped names. |
| 2 | M8 (CLI)     | `--cal-quality-good`, `--cal-quality-weak` driving `pool_quality_thresholds`   | Not implemented; `build_fl_models` uses fixed internal thresholds                                      | **Choose**: (a) ship the flags, or (b) drop them from the plan.  Recommendation: ship the flags — the thresholds materially affect the `weak`/`unusable` quality flag that downstream consumers branch on, and exposing them costs ~20 LOC. |
| 3 | M8 (Config)  | `CalibrationConfig.pool_quality_thresholds=(5_000, 200)`                       | Not present                                                                                            | Tied to #2.  If we ship the flags, we also add the dataclass field. |
| 4 | M8 (summary) | `summary.json` keys: `kappa_diagnostics`, `boundary_flux_gdna_summary`, `global_densities` | Keys: `global_densities`, `fl_models`, `diagnostics`, `n_multi_loci`, `c_base`, `mean_pi_gdna`         | **Reconcile**: `boundary_flux_gdna_summary` is folded into `diagnostics`; `kappa_diagnostics` is folded into `global_densities` (per-region table).  Rather than rename, document the actual schema in `docs/parameters.md` and amend the plan. |
| 5 | M8 (deprec.) | v1 CLI flags become deprecation-warning shims for one release cycle             | v1 flags removed outright in M8c                                                                      | **Document** — rigel has no released v1 surface (pre-1.0); skipping the deprecation cycle is acceptable.  Note this in the CHANGELOG. |
| 6 | M8 (`__init__.py`) | exports include `PoolFLModels`, `compute_pool_fl_models`                  | exports include `FLModels`, `build_fl_models`, plus `extract_global_counts`, `extract_rna_counts`, `extract_gdna_counts` (added during M7) | **Document** — superset of plan; no breaking change.  Update plan's export list. |

### 1.2 Tasks

- **M9.1.a** — Decide on issue #2 (`pool_quality_thresholds` flags).
  Default recommendation: **ship them**.  Estimated effort: half a day,
  3 unit tests, one integration test.  Skipped only if the user
  explicitly wants the surface kept minimal.
- **M9.1.b** — Amend `docs/calibration/calibration_v6_plan.md` §M7,
  §M8, §M8 export list, and the migration table to match the shipped
  symbol names (`build_fl_models`, `FLModels`, the actual
  `to_summary_dict()` keys).  Single commit, no code changes.

**Exit gate (M9.1):** plan and code agree symbol-for-symbol; CHANGELOG
draft started.

---

## 2. Validation matrix

The four scenarios in `calibration_v6_plan.md` §M9 are the acceptance
gate.  We rerun all of them with HEAD and record the results in
`results/m9_validation/`.

### 2.1 Synthetic scenarios (fast, local)

Driver: `scripts/sim/synthetic_sim_sweep.py` (already migrated to the
v6 API).  Configs live in `scripts/benchmark/configs/`.

| Scenario                | Input                                         | Acceptance metric              | Threshold        |
|-------------------------|-----------------------------------------------|--------------------------------|------------------|
| Pristine RNA-seq        | `gdna_fraction=0.0`                            | `mean_pi_gdna`                 | < 0.05           |
| 50/50 dna20m            | `gdna_fraction=0.5`, dna20m simulator          | `mean_pi_gdna`                 | ∈ [0.45, 0.55]   |
| TA1 nRNA-siphon         | high NTA1/TA1 ratio (existing `nrna_*` configs)| residual mRNA relative error   | < 5%             |
| Hybrid-capture          | low gDNA + capture targets                     | gDNA flux on internal exons    | within ±15% of truth |

**Tasks:**

- **M9.2.a** — Author/refresh sweep configs:
  `scripts/benchmark/configs/m9_pristine.yaml`,
  `m9_dna20m_50_50.yaml`,
  `m9_ta1_siphon.yaml`,
  `m9_hybcap.yaml`.
- **M9.2.b** — Run all four with `synthetic_sim_sweep.py`; collect
  outputs into `results/m9_validation/<scenario>/`.
- **M9.2.c** — Add `scripts/benchmark/m9_validate.py` that loads the
  four output dirs, computes the table above, and prints
  PASS/FAIL per row.  Fail loudly if any threshold is violated.
- **M9.2.d** — Capture the table in
  `docs/calibration/m9_validation_results.md`.

### 2.2 Armis2 13-condition matrix (cluster)

Driver: `scripts/benchmarking/` (already migrated).

| Step                                | Tool                                                           |
|-------------------------------------|----------------------------------------------------------------|
| Pre-flight (build/index up to date) | `pip install --no-build-isolation -e .`                        |
| Status                              | `python -m scripts.benchmarking status -c .../default.yaml`    |
| Run                                 | `python -m scripts.benchmarking run -c .../default.yaml`       |
| Analyze                             | `python -m scripts.benchmarking analyze -c .../default.yaml -o results/m9_armis2` |

**Tasks:**

- **M9.2.e** — Add a `m9_v6` rigel config to
  `scripts/benchmarking/configs/default.yaml` (no flag changes from
  the default v6 surface; the named config exists purely to namespace
  the output directory under `runs/<cond>/rigel/m9_v6/`).
- **M9.2.f** — Run the matrix on Armis2; copy the analyze output back
  to `results/m9_armis2/` in the repo (or link via README if too large).
- **M9.2.g** — Compare to the latest pre-M8 baseline (last published
  benchmark report).  Acceptance: no regression > 5% on per-condition
  median |log2FC| against truth for mRNA, nRNA, or gDNA classes.
- **M9.2.h** — Append an Armis2 results section to
  `m9_validation_results.md`.

**Exit gate (M9.2):** all synthetic thresholds met, Armis2 matrix
no-regression, results document published.

---

## 3. Documentation

### 3.1 `docs/MANUAL.md` — calibration section (≈1 page)

Audience: a user who runs `rigel quant` and wants to understand
what calibration *does* and *which knobs to turn*.

Outline:

1. What calibration solves (joint mRNA + nRNA + gDNA, library-wide
   `π_gdna` plus per-multi-locus priors).
2. Where the inputs come from (in-place C++ accumulator over the BAM
   scan; no extra pass).
3. The three knobs that exist:
   - `--cal-prior-ess` (Dirichlet pseudo-count for FL pool).
   - `--cal-nrna-weight` (relative prior weight for nRNA components).
   - `--cal-c-base` (base concentration for `assemble_priors`).
   - (If M9.1.a ships) `--cal-quality-good`, `--cal-quality-weak`.
4. What the `summary.json` calibration keys mean — short reference.
5. When to suspect calibration is misfiring (link to FAQ / debugging).

### 3.2 `docs/METHODS.md`

A condensed restatement of `calibration_v6_plan.md` §2 (Locked design)
suitable for inclusion in a methods paper or supplement.  Includes:
the region partition (8-state mask), global gDNA density derivation
(κ), the per-locus gDNA mass model, the Dirichlet prior assembly,
and the EM coupling via `prior_weight_rna`.

### 3.3 `docs/parameters.md`

Add a calibration subsection mirroring the existing parameter table
style.  One entry each for `prior_ess`, `nrna_weight`, `c_base`
(and the two quality thresholds if shipped).  Document the
`summary.json` schema (section 1.1, issue #4 above).

### 3.4 `CHANGELOG.md`

Major-release entry.  Headline: "v6 calibration redesign — joint
mRNA/nRNA/gDNA quantification with in-place per-region accumulator."
Sub-bullets: the M-by-M progression, the deleted v1 surface, the
new `summary.json` keys, the new CLI flags, the goldens regeneration.
Note the deviations from §M8 §M9 (issues #2 and #5 above).

### 3.5 Release tag

`v0.X.0` — first release with the v6 surface.  Tag from HEAD after
M9.1, M9.2, and M9.3 land.

**Exit gate (M9.3):** all four files updated; release tagged; PyPI
upload smoke-tested in dry-run.

---

## 4. Sequencing

```
M9.1.a  (optional)  ──┐
M9.1.b              ──┼─► M9.2.a-d  (synthetic)  ──┐
                      │                            │
                      └─► M9.2.e-h  (Armis2)    ──┴─► M9.3 (docs) ──► release
```

M9.1.b and M9.2.{a,e} can run in parallel.  M9.2 must finish before
the CHANGELOG can claim "validated."

---

## 5. Definition of done

- All ruff checks pass (`ruff check src/ tests/ scripts/`).
- All tests pass (`pytest tests/ -q`).
- `docs/calibration/m9_validation_results.md` exists and shows PASS
  on every synthetic threshold and "no regression" on the Armis2
  matrix.
- `docs/MANUAL.md`, `docs/METHODS.md`, `docs/parameters.md`,
  `CHANGELOG.md` updated and cross-linked.
- `calibration_v6_plan.md` amended to match the shipped surface
  (or the surface amended to match the plan, per M9.1.a).
- Release tag pushed.

---

## 6. Out of scope (deferred)

The following items from `calibration_v6_plan.md` §8 ("Open questions
deferred to future plans") remain deferred and are *not* part of M9:

- nRNA siphon root-cause (still flagged in repo memory as a fundamental
  identifiability issue when NTA1/TA1 ≥ 5).
- gDNA overestimation at low contamination (relative error 1.5–2.0× at
  `gdna_fraction=0.3`).
- Mappability-aware regions (separate plan in `docs/mappability/`).

These will be tracked in their own future plans; M9 closes the v6
calibration redesign as designed, not as we'd ideally want it after
two more rounds of root-cause analysis.
