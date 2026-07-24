# Variance foundation — plan (v4: approach E landed; 1/n separation deferred to the message task)

**Status:** v4. The foundational composition-precision model is **settled and its bug fix is landed + validated**
(Phase 1). Two independent analyses converge on it: the in-house 5-approach derivation workflow
(`variance_foundation_proposal.md`, numerically validated + independently re-verified in
`scratchpad/verify_foundation.py`, all claims PASS) **and** an external critique — both select **approach E**
(the single Schur-marginal scalar `τ_λ`) and both recommend **Option B** (defer the `1/n` composition/sampling
separation into the message-variance task). This plan records the settled foundation, the architectural
invariants, and a precise handoff of the deferred work.

Prereq: `variance_foundation_proposal.md` (the derivation), `CALIBRATION_ARCHITECTURE.md` (count-zero-info).

---

## 1. The settled model (approach E — verified)

The strand Beta-Binomial is **rank-1** (depends on `(λ,θ)` only through `p = ½+(κ−½)(1−f_g)sinθ`). So the honest
foundational composition precision is **ONE scalar** — the Schur-marginal gDNA-level precision:

```
    τ_λ = τ_density  +  c·a² · 1[single-strand]        a = ∂p/∂λ = −(κ−½)sinθ·f_g(1−f_g),  c = N_eff/(p(1−p))
```
* **single-strand (1-DOF):** tilt structurally locked ⇒ `τ_λ = τ_density + c·a²` (strand **pins** f_g);
* **AMBIG (2-DOF):** tilt a free nuisance ⇒ strand **cancels** out of f_g (Schur ⇒ 0) ⇒ `τ_λ = τ_density`.

A diagonal `(τ_λ, τ_θ)` is **prohibited** (it double-counts the rank-1 strand or drops its covariance). The tilt
precision `τ_θ` is not part of the gDNA-vs-RNA foundation; it re-enters `λ` only when a genuine tilt message
arrives (message task). `1/n` (density sampling) is **not** in `τ_λ` — it is a transport quantity.

## 2. Architectural invariants (codified — endorsed by both critiques)

| principle | rule |
|---|---|
| **Count-zero-information** | A raw count `N` carries ZERO composition information. It enters ONLY as statistical **power** — the strand Fisher `τ_λ` or the sampling noise `1/n` — never a composition vote (`CALIBRATION_ARCHITECTURE.md` §0). |
| **Parameter space** | Composition variances live in **log-odds** `λ = logit f_g` (precision `τ_λ`), emitted as `Var(log f_c)` via the `(1−f_g)²`/`f_g²` Jacobians. **Simplex fractions are NEVER used for variance storage.** (The code already honours this — confirm, do not change.) |
| **Local precision state** | Pass-0 local composition precision is a **single Schur scalar `τ_λ`**. A 2-DOF diagonal `(τ_λ, τ_θ)` is prohibited (rank-1 strand Fisher). |
| **Composition ⟂ sampling** | Fusion weight = composition precision `τ_λ` (mapped to `1/Var(log f_c)`). Transport/emit precision = composition ⊕ sampling `1/(Var(log f_c) + 1/n)`. The `1/n` belongs only on the transported density. |
| **Anchor fusion** | `struct_lock` (composition-certain) nodes use a **hard override** in fusion (adopt own belief), never an `∞` weight in soft `_fuse` — an `∞` weight makes `(∞·a+p·b)/∞ = NaN` and cascades through interior anchors. |

## 3. DONE — Phase 1 (the foundation + the bug fix), committed `c6df8c50`

* `node_init.build_node_init` gates the strand λ-term to single-strand nodes (`free_pos ^ free_neg`); AMBIG gets
  `τ_λ = τ_density`. Retires a live (bounded) phantom where AMBIG nodes were credited the single-strand term.
* Verified (`verify_foundation.py`): rank-1 exact; Schur = `I_density` to 1e-13; `(λ,t)` block-diagonalization
  exact; the AMBIG-strand-only λ-precision is N-invariant; the production over-credit is REAL but **bounded**
  (overdispersion caps the count power at `1/ω` — `τ(2N)/τ(N)=1.04`, NOT ∝N; a weak phantom, corrected vs the
  workflow's Poisson-toy claim).
* Unit test: an AMBIG node gets **zero** f_g precision from the strand at any N (`test_ambig_stranded_strand_gives_zero_fg_precision`).
* **A/B (improves, no regression):** refit=0 aggregate 0.1497→0.1491; refit=1 (ship path) 0.0910→0.0889,
  stranded/capON −0.0102; only stranded scenarios change (as the theory predicts); unstranded neutral.

## 4. DEFERRED to the message-variance task (Option B) — with the guard spec'd

The `1/n` is a **message-transport** quantity (the sampling noise of a density being *sent*), so the
composition/sampling separation belongs to the message task, where the message precision is re-derived anyway.
Doing it standalone now risks a spurious regression on a fragile condition (unstranded-capON) before the
propagation rules are finalized. The message task inherits these, precisely:

1. **Split the two roles.** `node_init` exposes the **fusion weight** `w_c = 1/Var(log f_c)` (composition only)
   AND the **transport seed** `1/(Var(log f_c) + 1/n)` (= today's `own_precision`, `bp_solver.py:377,418` use one
   `pg_own` for both today). Land the split **byte-identically first** (verbatim reuse of `own_composition_logvar
   + own_precision`, NOT the ULP-different `1/(1/τ+1/n)`), THEN flip the fusion to `w_c`.
2. **`struct_lock` HARD OVERRIDE in `_fuse`** — REQUIRED before the flip (else the `∞→NaN`, §2). This is a
   semantic change (anchors stop soft-blending), so it rides with the flip under the same per-condition A/B.
3. **Fold `1/n` into the transfer function** (the source's sampling noise, per its composition/enrichment form).
4. **Gate:** the falsifiable per-condition A/B (refit=0 AND refit≥1), unstranded-capON named; STOP if any
   condition regresses past the noise floor `ε` — "more correct" never overrides a measured regression.

## 5. Corrected census (for the message task — do NOT act on the stale draft)

* `belief.var_gdna` — **LIVE** (weights `_fit_gdna_hyperprior`, calibrate.py:174). Keep.
* `rna_pos/neg_frac_var` — **LIVE**: `node_init._rna` uses `isfinite(var_loc)` as the RNA liveness gate (a
  numeric no-op today — the ψ solve always returns finite `np.maximum` variances — but a real code reference).
  Re-express the gate as `free_s & n>0 & rho>_EPS` before ever deleting the field.
* `strand_likelihood.py` — production-dead but the **test oracle** for `_mixture_strand_loglik`
  (`test_rna_strand.py`). **Keep** (do not delete).
* `lam_var`/`theta_var` — dead in production, but `message_precision_mc.py` (the successor's harness) references
  `lam_var`. **Keep for the message task** (do not delete now).
* `var_pos`/`var_neg` — no ship reader; asserted in ~5 tests. Re-point, don't silently drop.
* `struct_lock = locked & is_region` — **region-only** by design (a G1 boundary seam is NOT struct_lock — a
  certain seam becomes a phantom-gDNA emitter). Preserve.

## 6. Foundation status: COMPLETE

The variance FOUNDATION is settled and landed (§1–§3). The message-variance task (`variance_model_handoff.md`)
is next; it builds the message precision on this foundation and absorbs the deferred §4 work. Nothing further is
needed here beyond keeping the invariants (§2) intact.
