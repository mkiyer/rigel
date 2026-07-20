# Calibration cleanup log — deferred tidy-ups (living)

Running list of cleanup owed as the NPMLE projection-variance work lands. Append; check off when done. Keep the
feature commits focused — cleanup rides in separate, reviewable passes.

## Dead / retired production code — ✅ DONE (Phase 0, 2026-07-17)
- [x] **`bp_solver.adjacent_disagreement_variance` + `_poisson_moment_var` + `_adjacent_edges` +
  `_adjacent_log_density_residuals`** — RELOCATED to `scripts/debug/_disagreement_variance.py`; deleted from
  `bp_solver` (+ the now-unused `BOUNDARY` import). The 6 debug harnesses repointed to the util.
- [x] **`message_precision.py`** — DELETED (was unimported in production).
- [x] Confirmed no `src/` references to the retired names remain (`grep -rn` CLEAN).

## Broken debug scripts (imported the DELETED `gdna_density_prior`) — ✅ DONE
- [x] `npmle_prototype.py`, `npmle_fusion.py`, `pass0_kde_prototype.py` — RETIRED (deleted; imported the deleted
  `gdna_density_prior`, already unrunnable). (No `npmle_variance.py` / `pass0_kde_landscape.py` /
  `pass0_kde_zero.py` / `landscape_real.py` present in the tree.)

## Naming / structure
- [ ] `GdnaRatePrior` is now more than a gDNA prior — it holds the enrichment-landscape NPMLE and serves BOTH
  the prior (`logprior`) and the message variance (`project`). Consider renaming to reflect the enrichment-
  landscape role. Churn: `calibrate`, `diagnostics`, `config`, `cli`, `tests`. Low priority (cosmetic).
- [x] The `σ²_transfer` formula — **RESOLVED to F1** by the formal derivation
  (`transfer_variance_formal_derivation.md`, 2026-07-17): message variance `= V_src + σ²_T`,
  `V_src = v_logfc + 1/M_src` (the BP source variance, includes the NPMLE prior via the solve),
  `σ²_T = var_proj(dst) + gap²` (F1, the conditional transport). **F2 rejected** (its `var_proj(src)` double-
  counts `V_src`); **F0 was my premature suggestion — withdrawn** (its `2h²` floor is a heuristic; F1's
  `var_proj(dst)` is the exact destination second-moment). No code change — the shipped form is correct.
- [x] **External statistical review — SIGNED OFF** (`transfer_variance_formal_derivation.md` §10, 2026-07-17):
  §7 adjudicated — using `O_d` in `σ²_T` is legitimate (Framing A = data-dependent CRF edge potential, the
  load-bearing justification; precondition "O_d touches only the variance, not the message mean" verified in
  `bp_solver._scan` — `s2t` enters `pr` only, `mo` is source-pinned). Framing B (EP switching transport)
  reproduces F1 **after** one correction (crossing reverts to the destination's *local* prior, not the global
  mixture). All three sub-questions resolved for the shipped form (exclude `var_proj(s)` ✓; identity-mean gag
  ✓; `h²`↔max-precision cap ✓). Kept-on-purpose deviation: `gap²` uses belief-free `μ_proj(s)`, not the
  posterior `m_s` — the refit seam.

## Reviewer implementation-risk determination (2026-07-17) — do NOT reintroduce retired knobs
The reviewer's Risks 1–2 target the retired discrete-stratification + Poisson-subtraction design, not the
shipped continuous NPMLE projection. Determination (full table in the formal doc §10.1):
- [x] **Risk 1 (step-functions):** N/A — `project()` is a softmax (smooth). **Adopted the guard only:** the
  refinement-based continuity regression test `test_projection_is_continuous_across_the_valley`. **Did NOT** add
  the prescribed sigmoid-between-strata (no strata exist).
- [x] **Risk 2 (Poisson over-subtraction → zero-clip):** N/A — `var_proj = between + h²` has a structural floor
  `h²>0`, cannot clip to 0. **Did NOT** add the prescribed Bayesian-shrinkage (2 new knobs for an impossible
  failure). The opposite-sign Poisson *double-count* stays deferred below (second-order).
- [ ] **Risk 3 (anchor representativeness) → the refit motivation:** for F1 this becomes "is the total-density
  NPMLE mode structure genuine enrichment geometry or low-count/RNA-smear noise?". **Build the coverage-
  stratified `var_proj` diagnostic** (checkpoint 1), and it is the strongest argument for the per-component
  **refit** (checkpoint 2 = the 5-pass convergence test = the next-phase examination). Nascent-factor hook
  (checkpoint 3) already present (message = independent per-component factors).

## Docs still describing the retired prior regime (npmle_roadmap §N4 tail)
- [ ] `CALIBRATION_ARCHITECTURE.md`, `calibration_prior_production_reference.md`, repo `CLAUDE.md` — still
  describe the 2-pass M2-background + KDE floor/global. Rewrite to the shipped pass-0 NPMLE + projection σ².

## Open engineering (from the design + reference derivation)
- [ ] **★ LOAD-BEARING — the DOMINANT flagship error (proven 2026-07-17, `flagship_prior_asymmetry_diagnosis.md`).**
  The `logP_g` left tail is CLAMPED (`gdna_rate_prior.logprior`: `np.interp(..., left=logP[0])`, a constant), so
  the gDNA prior lost the f_g→0 barrier that the Jeffreys `½log f_g` arm provided, while the RNA arm keeps its
  full `½log(1−f_g)` barrier at f_g→1. The resulting ASYMMETRIC ψ crushes strand-blind gDNA nodes to f_g≈0
  (the unstranded-capture leak): the "weak" prior nearly DOUBLES the strand-only error (0.329→0.622); restoring
  the `½log f_g` barrier moves node 1909's ψ argmax 0.002→0.700. Fix = a genuinely decaying left tail (so
  logP_g→−∞ as ρ_g→0); the symmetric endpoint is to WRITE `logP_r` and retire the lopsided Jeffreys RNA arm
  (`calibration_roadmap.md` "logP_r UNWRITTEN"; `reference_prior_derivation.md` §10.8). Validate at FULL-solve
  level + zero-gDNA/stranded controls before adopting. NOT the σ²_transfer, NOT strand deposit, NOT messages
  (which HELP 34:1).
- [ ] The Poisson double-count (Phase-2 open): `σ²_transfer` (from raw densities) mildly double-counts the
  source Poisson already in `1/M_src`. Second-order; revisit if it matters.
