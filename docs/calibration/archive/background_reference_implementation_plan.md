# Background-reference implementation plan (2026-07-17)

**Scope.** Transition the aggregate DNA-background reference from theory (`background_reference_derivation.md`,
reviewer-confirmed) to production, **simplifying the calibration codebase as it lands**. Phased, each phase
independently A/B-validated on the full 36-condition suite (24 original + 12 low-DNA/capture-strength) with the
hard controls, and shippable on its own. Discipline: monkeypatch-A/B before wiring; pin `bp_solver` /
`simplex_logodds` / `gdna_rate_prior` md5 before/after each phase (a concurrent editor has touched these).

**What the plan delivers, and what it honestly does not (yet).** Phases 0–3 make the prior **correct, safe, and
simpler**: the off-target background is measured as a scalar-with-variance, enters as a one-sided log-floor that
protects against the crush where DNA is present and vanishes where it is not, and the NPMLE fit stops
double-counting / mis-locating the sub-resolution mass. This fixes the off-target / low-density / zero-DNA
behavior with **no regressions**. It does **not** by itself lift the **enriched on-target** flagship nodes —
their DNA density lives on-target, entangled with RNA, so it needs the **iterative global-abundance** estimate
(Phase 4, the refit anchored by the background). The correct engineering order is: build the safe, simple
foundation first, then the harder iterative lift on top of it.

---

## Architecture: before → after

**Before.** `GdnaRatePrior.fit(total_density)` → one NPMLE mixture that (i) conflates DNA/RNA and (ii) by
Kiefer–Wolfowitz atomicity mis-locates the faint background as an atom at 0. `.logprior()` projects it onto `f_g`
with a **clamped left tail** (`np.interp(..., left=logP[0])`) — the level-blind asymmetry that crushes
strand-blind gDNA nodes. Retired σ²_imp machinery (`message_precision.py`,
`adjacent_disagreement_variance`, `_poisson_moment_var`) still sits in the tree.

**After.** A `BackgroundReference` scalar `(log ρ_bg, σ_bg)` measured from the pooled pure-DNA regions;
`GdnaRatePrior.fit(total_density, background=...)` with a **pinned background component** (formula 2, no
double-count); `.logprior()` with a **one-sided log-floor** at `ρ_bg` (data-gated, recedes as `ρ_bg → 0`)
instead of the clamped tail; the dead σ²_imp code removed. Net: fewer moving parts, one honest object,
level-aware behavior.

---

## Phase 0 — Cleanup first (health, zero behavior change)

Simplify before adding. All changes are dead-code removal or comment fixes; the full suite must stay green
(golden outputs unchanged).

- **Delete `calibration/message_precision.py`** (`adjacent_imputation_variance`) — the backfired honest-σ²_imp,
  retained but unused (`CLEANUP_LOG`). Drop any imports.
- **Delete `bp_solver.adjacent_disagreement_variance`, `_poisson_moment_var`, `_adjacent_edges`,
  `_adjacent_log_density_residuals`** — the retired total-density σ²_imp scalar, dead in the solve. If any debug
  harness still imports them (`transfer_variance_diag.py`, `phase0_projection_sigma2.py`), inline the two helpers
  there or retire the harness.
- **Fix stale comments/docstrings**: `calibrate.py:231-237` still says "σ²_transfer is now ZERO" (FALSE — the
  projection is on); `config.py:312` references the retired `message_precision.adjacent_imputation_variance`;
  the `NodeBelief.var_gdna` docstring says "Var(f_c)" (it is `Var(log f_c)`, `simplex_logodds.py:321`).
- **Retire the broken KDE-era debug scripts** (`npmle_fusion.py`, `npmle_variance.py`, `pass0_kde_*.py`,
  `landscape_real.py`) that import the deleted `gdna_density_prior`.
- **Grep-confirm** no remaining references before each deletion; run `pytest tests/ -q` + `ruff`.

**Gate:** 1079+ green, goldens byte-identical, ruff clean. Net LOC: **down**.

---

## Phase 1 — Measure the background reference `(log ρ_bg, σ_bg)`

New, small, pure-addition (no wiring yet).

- **New `calibration/background_reference.py`** — `BackgroundReference` (frozen dataclass: `log_rho_bg: float`,
  `sigma_bg: float`, `n_counts: int`, `eff_total: float`) + `measure_background(substrate, region_arrays,
  region_eff_len) -> BackgroundReference`. It pools the **count-observable** regions
  (`density_model.count_observable_masks` — intergenic + RNA-free introns, the machinery already exists) and the
  count-observable boundary sides, computing:
  ```
      Σg = Σ_{i∈B} count_i ,   ΣE = Σ_{i∈B} E_i^gdna
      ρ_bg = Σg / ΣE ,   σ_bg² = Var(log ρ_bg) ≈ 1/Σg        (Poisson; +½ continuity for Σg=0 → a finite
                                                              upper bound `1/ΣE`, σ_bg = ∞-ish → "no lower bound")
  ```
  `Σg = 0` (true/hidden zero) yields `ρ_bg = 0` with `σ_bg` reflecting the whole-genome upper bound `1/ΣE` — the
  honest "faint or none, below `1/ΣE`". No new constants: it is a pooled Poisson MLE.
- **Test** `tests/calibration/test_background_reference.py`: on the probe's 36 conditions the scalar reproduces
  `background_reference_probe.py`'s `ρ_bg` (gdna_none = 0 at every capture; monotone in DNA level;
  strand-independent). Determinism.

**Gate:** matches the validated probe values; no production path touched.

---

## Phase 2 — The one-sided log-floor in the gDNA arm (load-bearing safeguard)

Replace the clamped left tail with the data-driven, one-sided floor — the reviewer's §4 safeguard verbatim.
This is the V4-fix done **right** (level-gated, not a fixed barrier).

- **`GdnaRatePrior`** stores `log_rho_bg`, `sigma_bg` (passed through from `fit`, default `None` ⇒ floor off).
- **`GdnaRatePrior.logprior(fg_grid, mass, eff)`** adds, to the NPMLE projection term, the one-sided floor
  (all in natural-log-density; `logprior` already has `mass`, `eff`, hence `log ρ_g = log fg + log(M/E)`):
  ```
      floor(f_g) = −½ · ( max(0, log ρ_bg − log ρ_g) / σ_bg )²
  ```
  Zero for `ρ_g ≥ ρ_bg`; a soft Gaussian penalty (width `σ_bg`) below it. **Data-gated**: `ρ_bg → 0` ⇒
  `log ρ_bg → −∞` ⇒ the penalty region `[0, ρ_bg)` is empty ⇒ the floor vanishes — *dormant when DNA is
  absent*, so it cannot break gdna_none by construction. **`σ_bg`-softened**: a poorly-measured floor (strong
  capture, few counts) is wide/gentle; a well-measured one is sharp. This term also **replaces the clamped
  tail** (`left=logP[0]` → a genuinely decaying Jeffreys-slope left tail, `slope ½` in natural-log), since the
  floor now owns the low-end behavior.
- **Never a denominator / scale** — `ρ_bg` appears only inside the `max(0, log ρ_bg − log ρ_g)` bound; the only
  division is by the observed `M/E`. Satisfies the reviewer's stability constraint.
- **Wire in `calibrate.py`**: measure the background once (Phase 1), pass to `GdnaRatePrior.fit`. One extra call,
  one extra argument.

**A/B (monkeypatch `logprior` first):** all 36 conditions. **Gates:** (a) gdna_none stays ≈ 0 at every capture
(dormant floor); (b) low-DNA/strong-capture rows do **not** manufacture DNA (floor softens via `σ_bg`);
(c) off-target / depleted-node crush relieved where `ρ_bg` is substantial; (d) enriched flagship: expected
*partial* improvement only (floor `ρ_bg/ρ_tot` is small for enriched nodes — the honest limit; full lift is
Phase 4). Report the grid; wire only if the gates hold.

---

## Phase 3 — The pinned background component in the fit (formula 2)

Make the density estimate itself correct (no double-count, no atomicity) — improves the σ²_transfer projection
and the enrichment shape.

- **`GdnaRatePrior.fit(..., background=(log_rho_bg, sigma_bg))`** adds ONE component to the mixture whose
  **location is pinned** at `log ρ_bg` (kernel width `max(σ_bg, h)`), with a free weight; the free kernels keep
  their grid locations. The EM (`_em_weights`) fits `{w_bg, w_j}` unchanged in form — only the kernel matrix
  gains a fixed column. The sub-resolution mass flows to the pinned component (accurately located) instead of
  the atomic pile at 0. `background=None` ⇒ today's behavior (safe default).
- Because the low mass is now owned by the pinned component, the free grid's bottom can be **raised toward
  `1/E`** (the resolved range), a small simplification of the grid-support logic.
- **Test**: the pinned fit places a component at `ρ_bg` with the right weight on synthetic bimodal data; the
  faint-uniform-background case no longer collapses to an atom at 0 (contrast the naive fit).

**A/B:** all 36. **Gate:** σ²_transfer projection + solve hold or improve; no regressions vs Phase 2. (This phase
is mostly a fit-quality / robustness improvement; large mwae moves are not expected — verify it does not hurt.)

---

## Phase 4 — The enriched on-target lift: iterative global abundance (deferred, scoped)

The flagship's remaining gap. The enriched-on-target DNA level is not in the off-target floor; it is the
**library-wide DNA abundance**, measurable only from the on-target deconvolution → **iterative**. This is the
refit (`refit_loop_study_findings`), now anchored: (i) the background floor pins the low end and cannot be
crushed through; (ii) each pass re-estimates the global DNA fraction from the solved on-target DNA, lifting the
enriched nodes; (iii) `σ²_transfer` and the prior are refit on the *solved* belief with the background as the
fixed lower anchor. Scope it as the **next initiative** after Phases 0–3 land; it depends on the refit being
made convergent (cold-restart + confidence-weighting, per the study) and is the larger piece.

Parallel future: **`logP_r`** (the symmetric RNA arm) — replace the lopsided Jeffreys RNA reference once RNA
rates can be deconvolved, closing the ψ asymmetry at both ends. Also refit-gated.

---

## Phase 5 — Full validation + golden regeneration

- Run the complete battery: 36 conditions × {current, Phase-2, Phase-3}; the flagship, gdna_none, and
  low-DNA/strong-capture controls tabulated together.
- `pytest tests/ --update-golden` once the calibration outputs intentionally shift (Phase 2/3), review the diffs
  region-by-region (the golden scalars: `rna_sense_frac`, overdispersions, `gdna_density_global`, per-region
  masses).
- Update `CALIBRATION_ARCHITECTURE.md` + `calibration_prior_production_reference.md` to the shipped
  background-reference + pinned-component + one-sided-floor design (they still describe the retired 2-pass
  M2/KDE regime — a `CLEANUP_LOG` item).

---

## Code-health ledger (what gets simpler)

- **−1 module** (`message_precision.py`), **−4 functions** in `bp_solver` (retired σ²_imp), **−4 broken debug
  scripts** (Phase 0).
- **+1 small module** (`background_reference.py`, ~50 LOC) reusing existing masks — no new infrastructure.
- **`logprior` clamped-tail special-case → a principled decaying tail + one-sided floor** (removes the
  load-bearing bug, one clear mechanism).
- **One honest object**: `GdnaRatePrior` now carries its background explicitly (`log_rho_bg`, `sigma_bg`) rather
  than the level being implicit-and-wrong in a total-density fit.
- **No new magic numbers**: `ρ_bg`, `σ_bg` are a pooled Poisson MLE; the floor is `½·z²`; the Jeffreys tail slope
  `½` is derived. Nothing tuned.

## Sequencing & risk

Phases 0 → 1 → 2 → 3 land in order; each is independently green + A/B-gated + revertible. Phase 2 is the only
one that intentionally changes calibration output (the flagship/off-target behavior) and thus the goldens; it is
the one to A/B most carefully. Phase 4 is a separate initiative. Recommend: **land Phase 0 immediately** (pure
cleanup, de-risks the tree), then Phase 1, then A/B Phase 2 by monkeypatch before committing.
