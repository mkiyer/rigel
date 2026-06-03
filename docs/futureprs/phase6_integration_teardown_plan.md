# Phase 6 — Integration & Teardown: Implementation Plan (for critique)

The acyclic core (Phases 1–4) is built, shipped, and proven stable (Phase 5a). Phase 6
makes it the *only* calibration path: swap the engine inside `calibrate()`, slim the result
to the single-pass reality, rewrite the one consumer, and **tear out the entire old cyclic
machinery** — modules, tests, debug scripts, config knobs, error types, and stale comments.

> **Guiding constraint (user):** the result must be a *clean, concise, well-organized,
> readable, maintainable* codebase. Be tenacious — leave **no** vestigial concept behind.

---

## 0. Ground truth (corrects the stale docs)

Two long-standing notes are **wrong** and must be fixed as part of this phase:

1. **The pipeline is already wired end-to-end.** `run_pipeline` → `calibrate()`
   (`pipeline.py:867`) → `quant_from_buffer` → `assemble_priors` (`pipeline.py:744`) →
   per-locus EM. The CLAUDE.md "run_pipeline stops after calibrate() until PR 6" and the
   `pipeline.py:826` "skeleton … raises NotImplementedError" comment are **both stale**.
   `calibrate.py` is the *fully implemented old EM loop*. ⇒ **Phase 6 swaps the engine; it
   does not wire from scratch.** Low risk.

2. **The consumer contract is stable and barely changes.** `assemble_priors` already projects
   per-region quantities to loci by genomic overlap (`_project_regions_to_loci`) and emits
   `LocusPriors(alpha_gdna_add, alpha_rna_add, gdna_eff_len)`. The acyclic outputs map onto it
   almost 1:1. ⇒ **No C++ / EM changes. No `LocusPriors` change. No estimator change.**

### 0.1 Region↔locus geometry (the node→locus question, settled)

Regions (calibration) and loci (quant) are **not coordinate-aligned** — a region can straddle
two loci. There is exactly one mapping primitive in the codebase: `_project_regions_to_loci`
(`priors.py:40`), overlap-weighted fractional projection. **We reuse it unchanged.** Decision
C ("region mass = contained + left-side + right-side") is honored by summing a region's three
acyclic nodes *before* projection — exactly what the old `assemble_priors` already does for
RNA mass (`mass_d_contained + mass_d_left + mass_d_right`).

---

## 1. Target architecture (end-to-end, after Phase 6)

```
scan_and_buffer → AccumulatorPayload
        │
        ▼   calibrate()  (single pass, no loop)
   substrate ─▶ node_gdna_density ─▶ decode_regions / decode_sides ─▶ derive
        │            (count clue)         (joint count×strand)      (ρ₀, ω, eff-len)
        ▼
   CalibrationResult   { per-node gDNA/RNA mass, gdna_exposure_len, ρ₀, κ_rna }
        │
        ▼   assemble_priors()  (reuse _project_regions_to_loci)
   LocusPriors  { alpha_gdna_add, alpha_rna_add, gdna_eff_len }   ← UNCHANGED type
        │
        ▼   (unchanged) per-locus EM  →  n_t + 1 components
```

Every box except the three middle ones already exists. `calibrate`, `CalibrationResult`,
`assemble_priors` are rewritten; everything downstream is untouched.

---

## 2. New `CalibrationResult` (result.py — rewrite)

The single-pass acyclic calibrator has **no EM loop**, so all loop/posterior artifacts vanish.

### 2.1 Schema

```python
@dataclass(frozen=True, slots=True)
class CalibrationResult:
    # Per-region deconvolved mass — gDNA / RNA across the region's 3 nodes
    # (contained + the two boundary sides). The consumer sums the three.
    mass_g_contained: np.ndarray   # float64[R]
    mass_d_contained: np.ndarray
    mass_g_left:  np.ndarray       # right side of the left boundary
    mass_d_left:  np.ndarray
    mass_g_right: np.ndarray       # left side of the right boundary
    mass_d_right: np.ndarray

    # Exposure-weighted gDNA length per region: Σ_node ω_node · L_node.
    # The gDNA component's effective-length contribution (per-read normalizer).
    # Carries the clamped exposure (ω:=1 where a node saw nothing), so it is NOT
    # simply mass_g/ρ₀ in degenerate regions — this is why it is stored, not derived.
    gdna_exposure_len: np.ndarray  # float64[R], ≥ 0

    # Library scalars
    rho_0:     float               # ≥ 0  global gDNA density (mass/bp); headline diagnostic
    kappa_rna: float               # in [0,1]  RNA sense fraction used by the strand clue

    # Provenance
    n_regions: int
    config:    CalibrationConfig
```

### 2.2 Field-by-field teardown (every old field accounted for)

| Old field | Fate | Why |
|---|---|---|
| `mass_g_contained` … `mass_d_right` (6) | **KEEP** | now sourced from `decode_regions`/`decode_sides`; conservation-auditable |
| `omega`, `log_omega_var` | **REPLACE → `gdna_exposure_len`** | no Gamma posterior; ω is a pure ratio. Consumer needs Σ ω·L, not ω + its variance |
| `exposure_dispersion` (φ) | **PURGE** | no NB count-dispersion / exposure prior in the acyclic model (0 external consumers) |
| `rho_d_bb`, `rho_r_bb` | **PURGE** | strand overdispersion is `0.0` in the joint decode; never read except by logging |
| `n_iterations`, `converged`, `mass_change_history` | **PURGE** | single pass: no iteration, no convergence diagnostic |
| `rho_0` | **KEEP** | the headline contamination density (now from `derive`); diagnostic/QC |
| `kappa_rna` | **KEEP** | provenance: which strand balance fed the decode |
| `n_regions`, `config` | **KEEP** | provenance |

18 fields → **11**. `__post_init__` invariants shrink to: masses finite & ≥0;
`gdna_exposure_len` finite & ≥0; **`rho_0 ≥ 0`** (was `> 0` — must permit the graceful
zero-gDNA case, decision F); `kappa_rna ∈ [0,1]`. The non-finite "divergence sentinel" check
(and thus `CalibrationConvergenceError`) is **deleted**.

---

## 3. New `calibrate()` (calibrate.py — rewrite, ~15 lines of body)

Same signature as today (so `pipeline.py:867` call site is unchanged):

```python
def calibrate(payload, region_arrays, strand_model, gdna_fl_pmf, config) -> CalibrationResult:
    substrate  = CalibrationSubstrate.from_payload(payload, region_arrays)
    region_eff = region_eff_length(region_arrays.region_size_bp, gdna_fl_pmf)
    bside      = boundary_side_eff_length(gdna_fl_pmf, region_arrays.region_size_bp)
    mu_fl      = boundary_eff_length(gdna_fl_pmf)

    node_density = node_gdna_density(substrate, region_arrays, region_eff, mu_fl)
    kappa_rna    = fit_strand_balance(strand_model).kappa_rna

    reg          = decode_regions(substrate, region_arrays, node_density, region_eff,
                                  kappa_rna=kappa_rna, confidence=config.confidence, n_grid=config.n_grid)
    left, right  = decode_sides(substrate, region_arrays, node_density, bside,
                                kappa_rna=kappa_rna, confidence=config.confidence, n_grid=config.n_grid)
    derived      = derive(reg, left, right, region_eff, bside, region_arrays.ref_id)

    return CalibrationResult(
        mass_g_contained=reg.gdna_mass,  mass_d_contained=reg.rna_mass,
        mass_g_left=left.gdna_mass,      mass_d_left=left.rna_mass,
        mass_g_right=right.gdna_mass,    mass_d_right=right.rna_mass,
        gdna_exposure_len=derived.gdna_exposure_len,
        rho_0=derived.rho_0, kappa_rna=kappa_rna,
        n_regions=region_arrays.n_regions, config=config,
    )
```

Deleted from `calibrate.py`: the entire outer loop, `initial_hyperparameters`, the seed
constants `_RHO_D_BB_INIT` / `_PI_G_PRIOR_INIT`, and all imports of `density`, `estep`,
`exposure`, `mstep`, `sweep`. This is the heart of "acyclic": the module no longer references
a single piece of the cyclic machinery.

---

## 4. `derive.py` — one additive extension

`derive()` already computes the side-existence-masked node lengths internally. Extend
`DerivedExposure` with the consumer-facing aggregate so the result is self-contained and the
ω·L logic lives in exactly one place:

```python
gdna_exposure_len = region_omega*region_eff + left_omega*left_eff + right_omega*right_eff
```

(`left_eff`/`right_eff` are the existing `np.where(side_exists, bside, 0.0)` arrays.) The
per-node ω arrays stay on `DerivedExposure` for diagnostics (`proto_stability.py`,
`proto_derive.py`); only `gdna_exposure_len` flows into the result.

---

## 5. New `assemble_priors()` (priors.py — rewrite, lighter)

```python
def assemble_priors(calibration, region_arrays, multi_loci, *, prior_weight) -> LocusPriors:
    g_region = calibration.mass_g_contained + calibration.mass_g_left + calibration.mass_g_right
    d_region = calibration.mass_d_contained + calibration.mass_d_left + calibration.mass_d_right
    proj = _project_regions_to_loci(
        region_arrays, multi_loci, len(multi_loci),
        {"g": g_region, "d": d_region, "e": calibration.gdna_exposure_len},
    )
    w = float(prior_weight)
    return LocusPriors(
        alpha_gdna_add=w * proj["g"],
        alpha_rna_add =w * proj["d"],
        gdna_eff_len  =np.maximum(proj["e"], _GDNA_EFF_LEN_FLOOR),
    )
```

`_project_regions_to_loci` and `LocusPriors` are **unchanged**. The old version multiplied
`ρ₀ · (ω·L_phys)`; the new one uses the **deconvolved gDNA mass directly** (`g_region`), which
equals `ρ₀·Σω·L` in the non-degenerate case but is *more correct* in degenerate regions (it is
the actually-observed gDNA). The `ω·L`-with-clamping path survives only where it belongs — in
`gdna_eff_len`. The `inverse-variance weight w_r` from the old interface contract §6.2 is
**dropped** (decision E: no second shrinkage).

---

## 6. Teardown manifest (tenacious — nothing left behind)

### 6.1 Delete source modules (5)
`density.py`, `estep.py`, `exposure.py`, `mstep.py`, `sweep.py`
→ 22 calibration modules become **17**.

### 6.2 Delete tests (5)
`tests/calibration/test_density.py`, `test_estep.py`, `test_exposure.py`, `test_mstep.py`,
`test_sweep.py`.

### 6.3 Rewrite (3 + tests)
`calibrate.py`, `result.py`, `priors.py`; and `test_calibrate.py`, `test_result_schema.py`,
`test_priors.py` to the new schemas (no legacy fields).

### 6.4 Purge *within* surviving files
- **`errors.py`**: delete `CalibrationConvergenceError` (keep `CalibrationSubstrateError`).
- **`config.py` `CalibrationConfig`**: delete `max_outer_iterations`, `mass_rel_tol`,
  `exposure_dispersion_floor` and their `__post_init__` validation (all old-loop knobs;
  `boundary_split_factor` was already removed). Replace with the *only* real acyclic knobs —
  `confidence: float = 0.0` and `n_grid: int = 200` (today these are `joint_decode` function
  defaults; promoting them removes hidden magic). **Open Q (§10.1).**
- **`__init__.py`**: drop `CalibrationConvergenceError` from the export surface; keep the lean
  set actually imported elsewhere (`calibrate`, `CalibrationResult`, `CalibrationConfig`,
  `assemble_priors`, `LocusPriors`, `CalibrationSubstrate`, `SubstrateView`, `StrandBalance`,
  `CalibrationSubstrateError`).
- **`pipeline.py`**: rewrite the calibration log line (`875–886`) to the new fields
  (`R, rho_0, kappa_rna`, + a one-glance gDNA/RNA mass split); delete the stale `816–820`
  FL-note remnant and the `826–827` "skeleton raises NotImplementedError" comment.

### 6.5 Delete vestigial debug scripts (old-loop tracers)
Known: `scripts/debug/trace_zero_gdna.py`, `trace_calibration_oscillation.py`,
`phi_floor_exploration.py` (all import `initial_hyperparameters` / probe `exposure_dispersion`
oscillation). **Action:** sweep `scripts/debug/` in 6c and delete every script that imports a
deleted symbol; keep the acyclic protos (`proto_density_model`, `proto_strand_decode`,
`proto_joint_decode`, `proto_derive`, `proto_stability`).

### 6.6 KEEP (shared infra + acyclic core, 17 modules)
`signature`, `substrate`, `region_arrays`, `regions`, `effective_length`, `fl`,
`density_model`, `strand_decode`, `joint_decode`, `derive`, `strand_balance`, `strand_summary`,
`errors`(trimmed), `calibrate`(rewritten), `result`(rewritten), `priors`(rewritten),
`__init__`(trimmed).

---

## 7. Goldens & regression

The per-locus prior changes numerically (acyclic deconvolution ≠ old cyclic), so EM outputs
shift. Plan:
1. Land 6a–6c; run the full suite to see what moves.
2. **Inspect** the golden diffs for sanity (gDNA/RNA split direction, exposure factors) — do
   **not** blind-regenerate.
3. `pytest tests/ --update-golden` (78 golden files: 20 scenarios × ~4).
4. Reconcile `tests/test_estimator.py` `em_exposure_factor` assertions (Agent: lines ~578–611).
Native/accumulator tests are **unaffected** (verified).

---

## 8. Staging (reviewable, green between steps)

- **6a — engine swap.** Extend `derive` (§4); rewrite `result.py` (§2) + `calibrate.py` (§3);
  rewrite `test_result_schema.py`, `test_calibrate.py`. After 6a, `calibrate` imports none of
  the old modules ⇒ they are dead. *Gate: calibration unit tests green; `calibrate` produces a
  valid new result on the fixtures.*
- **6b — consumer.** Rewrite `priors.py` (§5) + `test_priors.py`. *Gate: `assemble_priors`
  green; end-to-end `run_pipeline` executes on a scenario without error.*
- **6c — teardown.** Delete the 5 modules + 5 tests + vestigial debug scripts; purge
  `errors`/`config`/`__init__`/`pipeline` remnants (§6.4–6.5). *Gate: `ruff` clean, no dangling
  imports, full non-golden suite green, module count = 17.*
- **6d — validate + goldens (folds in Phase 5b).** §7 golden regen + §9 validation gate.

Each sub-step is a commit; ship the set together (or as it stabilizes) per the established
rhythm.

---

## 9. Validation gate (Phase 5b lives here)

End-to-end is the *only* place the deferred questions become answerable:

1. **Vindication, end-to-end.** Re-run the `nrna_dc g20` regime *through the EM*: RNA recovery
   stable across seeds/depth (Phase 5a showed the calibration estimate is stable; confirm it
   survives the prior→EM hop).
2. **Risk-C (the real new risk).** The acyclic ω can reach ~0 (pure-RNA / capture-off-target),
   which the **old Gamma posterior never allowed** (it floored ω at the prior mean 1). A near-0
   ω shrinks `gdna_exposure_len`; once it hits `_GDNA_EFF_LEN_FLOOR = 1.0` the gDNA component
   looks like a 1 bp target (high per-read likelihood) and may **over-attract** unspliced reads
   despite a near-0 gDNA prior. **Measure** on a capture scenario: RNA recovery at off-target
   loci with `enable_gdna=1`. If depressed → add the *minimal* fix (raise the floor, or shrink
   ω→1 mildly) — **only if benchmarks warrant** (decisions B + E: measure first, no
   double-shrink by default).
3. **Nascent-confound watch.** The accepted gDNA over-call biases `alpha_gdna_add` high;
   quantify the RNA-recovery cost end-to-end. If severe, it strengthens the case for the
   deferred 3-component fix (out of scope here).
4. **Full suite + benchmark recovery** sane across scenarios; goldens reviewed (§7).

---

## 10. Open questions — RESOLVED (locked 2026-06-03)

1. **Knobs:** promote BOTH. `confidence` → high-value user parameter. `n_grid` → advanced/
   technical (not expected to change), but **prove it during benchmarking** (and watch its
   performance implications).
2. **Store ω explicitly:** YES — store the per-node exposures in the result for calibration-
   module clarity and straightforward debugging/diagnostics/QC. (Schema §2.1 gains
   `omega_contained/left/right` alongside `gdna_exposure_len`.)
3. **`prior_weight`:** keep default `1.0`; validate across benchmarks later.
4. **Risk-C:** measure-first confirmed — quantify the effect of the **raw** exposure-factor
   fits end-to-end; add shrinkage only if benchmarks warrant.
5. **Cadence:** one bundled Phase 6.

### (original open questions, for reference)

1. **`CalibrationConfig` knobs.** Promote the two real acyclic knobs (`confidence=0.0`,
   `n_grid=200`) into config, or keep them as `joint_decode` defaults and let
   `CalibrationConfig` become essentially empty (kept only for provenance/forward-compat)?
   *Recommendation: promote both — removes hidden magic, keeps the no-magic-numbers rule honest.*
2. **Store ω explicitly?** I store only `gdna_exposure_len` (the consumer-facing Σ ω·L) and let
   debug scripts recompute per-node ω. Alternative: also store the 3 ω arrays for at-a-glance
   QC. *Recommendation: keep lean; ω is reconstructable + lives on `DerivedExposure`.*
3. **`prior_weight` scale.** Default `1.0`. The new prior uses deconvolved gDNA mass (≈ the old
   `ρ₀·ω·L`), so the scale is ~unchanged — but I'll confirm empirically in 6d and flag if it
   needs a new default.
4. **Risk-C policy.** Confirm measure-first: ship 6d with **no** shrinkage and add it only if
   the capture benchmark shows off-target gDNA over-attraction.
5. **Ship cadence.** One bundled Phase 6, or ship 6a–6b (engine + consumer) first, then 6c–6d
   (teardown + goldens) as a follow-up? *Recommendation: bundle — the teardown is what keeps
   the tree honest, and goldens must regen against the final code.*

---

## 11. Budget check (Q6: ≤25 modules, ≤~8 calibration constants)

- **Modules:** 22 → **17**. Well under 25, and net-removing (the goal).
- **Constants (calibration):** `_PI_EPS`, `_JEFFREYS` (joint_decode), `_TRAFFIC_PSEUDOCOUNT`
  (density_model), `_GDNA_EFF_LEN_FLOOR` (priors), `POOL_EB_PRIOR_ESS`, `N_FL_POOLS` (fl) — 6,
  plus the soon-promoted `confidence`/`n_grid` *config defaults* (not buried constants). The
  teardown deletes `_RHO_D_BB_INIT`, `_PI_G_PRIOR_INIT`, `_RHO_D_BB_FALLBACK`, `_PI_CLIP`,
  `exposure_dispersion_floor`, and the old config knobs. **Net: fewer constants, all justified.**
```
