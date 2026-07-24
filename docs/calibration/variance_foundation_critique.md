# Variance-foundation plan — adversarial critique ledger

Raw findings from a 6-lens adversarial critique of `variance_foundation_plan.md` (v1 draft), 2026-07-24.
27 findings: **4 blocker, 16 major, 7 minor**. Each is grounded in file:line. The v2 plan (§9 ledger) maps
every finding to its resolution. Kept so the next session can verify the reasoning and critique further.

## BLOCKERS

**B1 — nan at the struct_lock anchor (numerics).** `own_composition_logvar` returns `Var=0` for struct_lock
nodes; today `own_precision(n,0)=n` (finite). Phase C removes `1/n` from the fusion weight ⇒ the locked-node
weight becomes `1/Var=1/0=∞` ⇒ `_fuse` computes `(∞·a+p·b)/∞ = nan` (bp_solver `_fuse`). Intergenic anchors are
**interior** chain nodes whose running belief propagates ⇒ the nan poisons all downstream nodes. The anchor is
the load-bearing gDNA source on unstranded data. **Fix:** struct_lock = HARD OVERRIDE in the fusion (adopt own
belief, skip `_fuse`); emit precision stays `n`; add an interior-anchor no-nan test.

**B2 — strand Fisher is RANK-1; a diagonal `(τ_λ,τ_θ)` cannot represent it (τ_θ derivation).** The strand BB
likelihood depends on `(λ,θ)` only through `p = ½+(κ−½)·t`, `t=(1−f_g)sinθ` (simplex.py:41-47). ⇒ Fisher =
rank-1 `c·∇p∇pᵀ`, with off-diagonal `I_λθ = −c(κ−½)²f_g(1−f_g)²sinθcosθ ≠ 0` and a null ridge. A diagonal state
either drops the covariance or, populating both `τ_λ` and `τ_θ` from the same scalar, **double-counts the
strand**. Also `I_strand_θ=0` at single-strand (cosθ=0), contradicting the `τ_θ=∞` gate (that ∞ is STRUCTURAL,
not Fisher); and the arcsin-θ orthogonality is a **prior** property, not a **likelihood** one. **Fix:** DEFER
τ_θ; resolve rank-1 first (single-precision-on-`t` + null ridge / full 2×2 / derive from the ψ joint).

**B3 — Phase C gate is unfalsifiable (test).** "expect a shift, ideally ≤ current" permits a worse aggregate to
pass and pre-frames a regression as "more correct." **Fix:** explicit two-sided per-condition gate on a
committed baseline with a noise-derived tolerance; any per-condition regression FAILS.

**B4 — validation is refit=0-only; the ship path (hyperprior re-solve) is never gated (test).** The change
alters `f_g` and `var_gdna`, both fed to `_fit_gdna_hyperprior` (calibrate.py:173-174) ⇒ the refit moves, yet no
phase runs refit>0. **Fix:** add a refit≥1 A/B arm to Phase C and E.

## MAJORS

* **M1 (numerics/layering/test ×3) — Phase B adapter `1/(1/τ+1/n)` ≠ `own_precision`.** Differs at ULP from
  `n/(n·v+1)`; drops the per-arm Jacobian (`Var(log f_g)=(1−f_g)²/τ_λ`, `Var(log f_r)=f_g²/τ_λ`, not `1/τ_λ`);
  ignores the single shared `n_node`, the lock/∞/nan-mask, and the `live`+`isfinite(var)` gates. **Fix:** adapter
  = verbatim reuse of `own_composition_logvar+own_precision`; commit the byte-identical test (not a scratchpad
  script).
* **M2 (numerics + downstream) — `rna_pos/neg_frac_var` are LIVE.** `node_init._rna` uses `isfinite(var_loc)` as
  the RNA liveness gate (node_init.py:271). "Dead downstream" is wrong; deleting them (Phase D) is a hard break
  and dropping the gate turns on phantom RNA. **Fix:** re-express the gate (`free_s & n>0 & rho>_EPS`) first.
* **M3 (downstream) — deleting `strand_likelihood.py` breaks tests.** `test_rna_strand.py` imports `strand_loglik`
  as the `_mixture_strand_loglik` no-regression oracle. **Fix:** keep-as-oracle or port the assertions.
* **M4 (layering) — "reference-free" is false; `NodeState` is pass-DEPENDENT.** `tau_lam` is evaluated at
  `fg_loc` (the prior-shifted mode) ⇒ differs pass-0 vs refit. **Fix:** drop the label; document; optionally
  freeze at `fg_ref`.
* **M5 (completeness) — measured-spliced deferral zeroes the boundary self-solve on unstranded.** A boundary's
  own RNA precision is only `tau_lam` (=0 unstranded) ⇒ precision-0 own belief, riding entirely on the relay.
  The concept doc lists measured-spliced as a foundation source. **Fix:** promote it (Poisson `1/n_spliced`) or
  explicitly document boundaries aren't self-solvable + why.
* **M6 (completeness) — the debug diagnostic contract is omitted.** ~30 debug scripts consume `_ni`/`_capture`
  keys; the Phase-C A/B is run WITH `pass0_node_dissect.py`. A silent rename breaks the validation tool. **Fix:**
  §8 contract + a Phase-B "dissect replays byte-identical" checkpoint.
* **M7 (completeness) — `free_pos/free_neg` + per-strand structural locks omitted from `NodeState`.** `struct_lock`
  covers only pure-gDNA certainty; `n_pos=0` does not encode strand forbiddance. **Fix:** state `free_pos/free_neg`
  stay in `NodeStatics`, applied unchanged; test a forbidden strand stays `f=0` at any count.
* **M8 (test) — no per-condition guardrails.** unstranded-capON (the known landmine) can regress while the
  aggregate + stranded arm improve. **Fix:** per-condition floors; name unstranded-capON.
* **M9 (test) — test accounting incomplete.** `@xfail(strict=True) test_mature_measurement_disagreement_silenced`
  may flip to XPASS (hard CI fail); ~5 tests assert `var_pos/var_neg`; `test_measured_intergenic_is_poisson_precision`
  asserts `prec_g==n` on the field Phase C redefines. **Fix:** enumerate per-phase test edits; re-derive the xfail
  (don't widen its bound); re-value the intergenic test on BOTH fusion weight and emit.
* **M10 (τ_θ) — I_strand_θ=0 at single-strand contradicts the ∞ gate** (structural, not Fisher) — deferred.
* **M11 (τ_θ) — MC target ambiguous:** conditional curvature `I_strand_θ` ≠ marginal `theta_var` (rank-1 coupling
  inflates the marginal) — deferred; validate the conditional at fixed λ.
* **M12 (τ_θ) — prior-orthogonality (arcsin θ) ≠ likelihood-orthogonality** — deferred; the diagonal claim must
  come from the likelihood (which is rank-1).
* **M13 (layering) — refit not A/B'd** (= B4).

## MINORS

* **m1 — "double-count" overstated.** `Var(log f)` and `1/n` are variances of different quantities; their sum is
  the honest density-estimate variance — a modelling choice, not a proven bug. **Fix:** reframe; let the A/B
  adjudicate; don't preordain.
* **m2 — `pg_own` serves two roles** (fusion weight + transport seed). **Fix:** two named arrays; write
  `1/Var(log f_c)` explicitly, not the bare `τ`.
* **m3 — deleting `lam_var` undercuts the successor.** `message_precision_mc.py:94-96` references
  `NodeDeconv.lam_var`. **Fix:** announce the retirement in Phase D / I7.
* **m4 — θ deadband misuses gDNA terms.** `σ²_d` includes `1/N_gdna+ω_g`; the tilt is RNA-internal. **Fix:**
  derive the θ deadband from RNA stats alone — deferred with τ_θ.
* **m5 — struct_lock region-only scoping not carried into `NodeState`.** Risk of re-locking G1 seams (phantom
  emitter). **Fix:** encode as a NodeState invariant + regression test; reconcile the concept doc line.
* **m6 — adapter omits the live/finite-var gates** (subset of M1). **Fix:** carry the full predicate.
* **m7 — Phase A MC under-specified; baseline is a remembered number on a "NOT A WORKING VERSION" checkpoint;
  goldens run `_scan` not the unified path.** **Fix:** concrete toy+seed+tolerance; pin the baseline on a
  fresh run; flip `RIGEL_UNIFIED` for a golden run or note goldens don't cover the change. *(Positive note: τ_θ
  reuses the existing `disc`/`N_eff` — no new model constant; only test tolerances need justifying.)*
