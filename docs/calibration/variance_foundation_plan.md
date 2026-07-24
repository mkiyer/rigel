# Variance foundation — implementation plan (v3: build approach E, verified)

**Status:** v3, the plan of record. Builds the derivation in `variance_foundation_proposal.md` (approach E — the
Schur-marginal scalar `τ_λ`), which was independently verified (`scratchpad/verify_foundation.py`, all claims
PASS) and hardened against `variance_foundation_critique.md` (27 findings). This plan does NOT build message
propagation; it solidifies the foundation the message variance model will build on.

---

## 1. The model (verified)

A node's composition is 2-DOF; `λ = logit f_g` (gDNA level), `θ = arcsin τ` (tilt). The strand Beta-Binomial is
**rank-1** (depends only on `p = ½+(κ−½)(1−f_g)sinθ`). The honest **foundational composition precision is ONE
scalar** — the Schur-marginal gDNA-level precision:

```
    τ_λ = τ_density  +  c·a² · 1[single-strand]         a = ∂p/∂λ = −(κ−½)sinθ·f_g(1−f_g),  c = N_eff/(p(1−p))
```

* **single-strand (1-DOF)**: `θ` structurally locked ⇒ `τ_λ = τ_density + c·a²` (strand **pins** f_g).
* **AMBIG (2-DOF)**: `θ` free ⇒ **strand cancels out of f_g** ⇒ `τ_λ = τ_density` (0 from strand, for all N).

Verified: Schur = `I_density` to 1e-13; the collapse is EXACT (in `(λ,t)`, `t=f₊−f₋`, the likelihood is exactly
block-diagonal, cross-derivative 0.0); the honest AMBIG λ-precision is N-invariant. **`1/n` (sampling) is NOT
part of `τ_λ`** — it is added only to a transported density at message time.

**The bug this fixes:** `node_init.strand_evidence` today credits the single-strand term `c·a²` to **every**
node, including AMBIG (verified: returns `τ_λ=1.33` where the honest value is `0`). It is a **bounded, weak**
phantom (overdispersion caps the count power at `1/ω`; `τ(2N)/τ(N)=1.04`, NOT ∝N) — a real correctness defect on
AMBIG nodes, not the dominant error. The fix = **gate `c·a²` to single-strand nodes**.

## 2. Current state (census — CORRECTED against the critique)

* `belief.var_gdna` — **LIVE**: weights `_fit_gdna_hyperprior` (calibrate.py:174). Only ψ-variance ship-path reader.
* `rna_pos/neg_frac_var` — **LIVE**: `node_init._rna` uses `isfinite(var_loc)` as the RNA liveness gate (node_init.py:271).
* `strand_likelihood.strand_loglik` — production-dead but the **test oracle** for `_mixture_strand_loglik` (test_rna_strand.py). **Keep.**
* `lam_var`/`theta_var` — dead in production; `message_precision_mc.py` (successor) references `lam_var` — announce on removal.
* `var_pos`/`var_neg` — no ship reader; asserted in ~5 tests (re-point, don't silently drop).
* `struct_lock = locked & is_region` (region-only; G1 seams excluded — a seam that is certain becomes a phantom-gDNA emitter).
* Relay uses one `pg_own` for BOTH the fusion weight AND the transport seed.

## 3. Phases (each: green suite + benchmark A/B checkpoint)

### Phase 1 — the `τ_λ` correction (approach E core; the bug fix). *Priority.*
* In `node_init`: gate the strand λ-term to single-strand nodes — `tau_lam = tau_density + where(free_pos ^ free_neg, i_strand, 0)`. AMBIG τ_λ = density only. (Keep `own_composition_logvar`/`own_precision`; only the `tau_lam` assembly changes.) Rename/clarify: `tau_lam` is now the honest Schur-marginal λ-precision.
* Unit tests: an AMBIG node → strand contributes **0** to `τ_λ` at every N (the N-invariance guard); a single-strand node → `τ_λ` includes `c·a²` (unchanged); unstranded → 0 (deadband, unchanged).
* **Gate (NOT byte-identical on AMBIG-stranded):** the §4 falsifiable per-condition A/B (refit=0 AND refit≥1). Expect near-neutral (the phantom is weak/bounded; unstranded + single-strand unchanged), no per-condition regression.

### Phase 2 — the composition/sampling separation (1/n out of the fusion weight).
* Foundation exposes `τ_λ` (composition, no 1/n) + the per-component counts. The relay:
  - **fusion weight** `w_c = 1/Var(log f_c)` (composition only, with the `(1−f_g)²`/`f_g²` Jacobians);
  - **transport/emit seed** = `1/(Var(log f_c) + 1/n)` (= today's `pg_own`, unchanged) → `_damp` adds σ²_transfer;
  - **struct_lock = HARD OVERRIDE** in fusion (adopt own belief, skip `_fuse` — never an ∞ weight → the interior-anchor nan). Emit precision stays `n`.
* Provide a **verbatim byte-identical adapter first** (reuse `own_composition_logvar+own_precision`, NOT `1/(1/τ+1/n)`) to prove the extraction, THEN flip the fusion weight. Preserve the `_capture` diagnostic keys (the debug loop runs the A/B).
* Tests: the interior-anchor **no-nan** test; re-value `test_measured_intergenic_is_poisson_precision` (fusion weight vs emit); audit the `strict=True` xfail `test_mature_measurement_disagreement_silenced` (may flip to XPASS — re-derive, don't widen its bound). **Gate:** §4 A/B.

### Phase 3 — cleanup + honest deletions.
* Before deleting `rna_*_frac_var`: re-express the RNA liveness gate as `free_s & n>0 & rho>_EPS`; test an unresolved AMBIG strand emits **zero** own RNA. Re-point the `var_pos/var_neg` assertions.
* `lam_var/theta_var`: delete after Phase 2; announce the retirement in `message_precision_mc.py`. Keep `strand_likelihood.py` (test oracle).
* Decide `var_gdna`→hyperprior: keep as the documented weight (low-risk) or migrate to `τ_λ` (A/B'd). Recommend keep.
* Goldens regenerated (note they exercise `_scan` at HEAD — run a `RIGEL_UNIFIED` golden pass or state they don't cover the change). Update ROADMAP/CLAUDE.md/memory.

## 4. Validation — falsifiable, per-condition, both refit arms

* **Baseline:** capture per-scenario `mwae`/`corr`/region-boundary/`gdna_none` on the CURRENT tree (fresh), at refit=0 AND refit≥1 → a committed baseline file (extend `scripts/debug/pass0_oracle_bench.py` with a refit arm).
* **Noise floor `ε`:** derive from the suite's run-to-run spread (`OMP_NUM_THREADS=1`).
* **PASS iff:** aggregate `mwae` not up > `ε`; **every** {stranded,unstranded}×{cap-on,cap-off} not up > `ε` (name **unstranded-capON** — the known landmine); `gdna_none` phantom mass not up; same at **refit≥1**. Any per-condition regression past `ε` FAILS — "more correct" does not override a measured regression.

## 5. Files
`node_init.py` (Phase 1 gate; Phase 2 expose τ_λ + counts, drop the 1/n from composition), `bp_solver.py`
(Phase 2 fusion-weight vs transport-seed split, struct_lock hard-override; Phase 3 dead-code), `calibrate.py`
(var_gdna decision), `simplex_logodds.py` (Phase 3 drop dead outputs after the gate re-express),
`tests/calibration/test_node_init.py` (+ the AMBIG-zero, no-nan, re-valued tests), `pass0_oracle_bench.py`
(refit arm). Keep `strand_likelihood.py`.

## 6. Open items (carried)
* **Measured-spliced boundaries** self-solve to `τ_λ=0` on unstranded data (own RNA is strand-only). Decision: document-and-defer (ride the relay graft) OR promote to a Poisson `1/n_spliced` own-RNA belief. **State it.**
* **Pass-dependence:** `τ_λ` is evaluated at the prior-shifted mode; document (not "reference-free").
* **τ_θ / tilt precision:** deferred; re-enters `λ` as `r_θ` when a tilt message arrives (message-layer task).
