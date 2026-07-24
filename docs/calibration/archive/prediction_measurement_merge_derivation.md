# Item E — the PREDICTION ⊕ MEASUREMENT merge: derivation

**Status:** derived, **not implemented**. Unblocks R2 (`message_arithmetic_reconciliation.md` Part D).
**Scope:** pass-0 message propagation. The RNA message a splice-junction boundary sends to its exon flank.
**Date:** 2026-07-22.

---

## 1. The problem

A splice-junction boundary emits ONE RNA message, but its RNA has **two provenances**:

| component | what it is | provenance |
|---|---|---|
| ρ_m | mature RNA, from the **spliced** junction fragments | **MEASURED** — motif-stranded, unambiguously RNA. No composition solve. |
| ρ_n | nascent RNA, from the deconvolved **unspliced** crossing | **IMPUTED** — requires the gDNA/RNA composition solve. |

They are lumped into one message. The owner's question: *does the imprecise nascent get "free confidence" by
riding along with the precise spliced?*

**Today's code answers this wrongly.** It does:

```python
pr = _pred_precision(n_eff, v_logfr, s2t)   # PREDICTION channel
if S > 0: pr += S                            # MEASUREMENT channel — precisions ADDED
```

This was measured to inflate `gdna_none` phantom gDNA **+49 %** when the σ²_transfer damping was removed
(Part D). The damping had been silently holding it together.

---

## 2. The theoretical violation: replicate-measurement algebra on an additive decomposition

There are two distinct combination rules, and the code uses the wrong one:

* **Two independent observations of the SAME quantity θ** → precisions **add**: `τ = τ₁ + τ₂`.
  (The Bayesian/Gaussian product. This is what `pr += S` implements.)
* **Two components that SUM to the quantity**, `θ = θ₁ + θ₂` → **variances** add, and in log space they add
  **share-weighted**.

Mature and nascent are not two looks at the same number. They are two **addends** of the RNA density:

```
    ρ_r = ρ_m + ρ_n
```

Adding their precisions asserts that observing the mature makes us more certain about the *nascent*, which is
false — the mature tells us nothing whatsoever about how much nascent there is. **That is the violation.**

### 2.1 The correct combination (delta method)

With ρ_m ⟂ ρ_n:

```
    Var(ρ_r) = Var(ρ_m) + Var(ρ_n)
    Var(log ρ_r) = Var(ρ_r)/ρ_r²
                 = w_m²·Var(log ρ_m)  +  w_n²·Var(log ρ_n)
    where   w_m = ρ_m/ρ_r,   w_n = ρ_n/ρ_r,   w_m + w_n = 1
```

**Share-weighted, and the weights are SQUARED.**

### 2.2 This answers the "free confidence" question exactly

* **Nascent dominates** (w_n → 1): `Var(log ρ_r) → Var(log ρ_n)` — the message inherits the *imputation's*
  full uncertainty. **No free confidence.**
* **Mature dominates** (w_m → 1): `Var(log ρ_r) → Var(log ρ_m)` — the measurement's small variance. Correct:
  the total really is well determined, because the uncertain part is a negligible share.
* **Mixed**: quadratic in the shares. A minority component contributes *quadratically* little — which is right:
  a large *relative* error on a small addend is a small *absolute* error on the sum.

So the weak component does **not** piggyback. But note the converse, which is the part today's code gets
backwards: the strong component does not get *diluted* either. Precision-addition is wrong in **both**
directions.

---

## 3. Provenance enters through the component variances

```
    Var(log ρ_m) = 1/n_m                                   MEASURED: Poisson counting ONLY
    Var(log ρ_n) = 1/n_u  +  Var(log f_r^src)              IMPUTED: counting + the COMPOSITION solve
```

* `n_m` — the integer **spliced** flux (`spliced_n_*`, already plumbed).
* `n_u` — the integer **unspliced** crossing flux (`n_unspl_*`, plumbed by R1).
* `Var(log f_r^src)` — the source's reference-free composition uncertainty (the existing τ term).

The measurement's variance contains **no composition term at all** — that is the entire content of "spliced is
pure RNA". Its provenance advantage is real, and it is expressed here and nowhere else.

---

## 4. The transfer: BOTH components are still predictions

> *"Both of these are still imputed. They still need to jump across the boundary onto the region. We're
> predicting that the flux of fragments at a single point — the boundary — is predictive of the contained
> fragments over an entire region. We aren't directly measuring the exon."* — owner

This is the load-bearing correction to my earlier R2 reasoning. The spliced count is a **direct measurement of
the junction**, and a **prediction about the exon**. So a transfer variance applies to *both* components:

```
    Var(mode) = w_m²·(1/n_m)  +  w_n²·(1/n_u + Var(log f_r^src))  +  σ²_xfer
```

`σ²_xfer` is the **point-flux → region-containment** imputation variance: the RNA density at the junction is
not exactly the RNA density averaged over the exon (5′/3′ coverage gradients, partial probe coverage, local
heterogeneity). It is a genuine, irreducible model-error term.

**This resolves the R2 confusion.** Two different things were conflated under one name:

| term | meaning | verdict for the measurement channel |
|---|---|---|
| σ²_transfer (enrichment-cliff proxy) | `(μ_proj[dst] − μ_proj[src])²` — a stand-in for a mode that was wrong across cliffs | ❌ must NOT damp a count (R2's original claim — still correct) |
| **σ²_xfer (point→region)** | spatial heterogeneity of the boundary→exon imputation | ✅ **applies to BOTH components, including the measurement** |

R2 was right to remove the first and wrong to conclude the measurement is undamped. It is damped — by the
*correct* term.

---

## 5. The mode, and the structural guarantee that kills the phantom

The message mode is the **sum**, not the prediction alone:

```
    mode = log(ρ_m + ρ_n) + [eff-length frame terms]        with ρ_n floored at 0
```

This is what prevents the +49 % phantom regression. Today the mature absorption is subtracted from the merged
density and can drive the whole thing to ~0, so the message can assert *"confidently no RNA"*. Under this
derivation that is structurally impossible:

> **The measured mature is a LOWER BOUND on the RNA that the prediction cannot erase.**
> If absorption drives ρ_n → 0, then ρ_r = ρ_m > 0, w_m = 1, and the message says *"there is at least this
> much mature RNA here"* at the measurement's (transfer-limited) confidence — a **positive** assertion.

`f_g → 1` by RNA-starvation therefore cannot happen through this channel, which is exactly the failure mode the
`gdna_none` guard caught.

---

## 6. One message, not two factors

Tempting alternative: emit the measurement and the prediction as two separate λ-factors and let the EP fold
combine them. **This is wrong** — factors on the same log-fraction axis combine as a *product* (precisions
add), which is the §2 violation again. The components are additive in **density**, not multiplicative in
**fraction**. They must be summed *before* the fraction is formed. One message, share-weighted variance.

---

## 7. Approximation notes (honest limits)

1. **First-order delta method.** `w_m`, `w_n` are themselves functions of the uncertain ρ_n, so §2.1 is a
   first-order propagation. The leading term dominates whenever the shares are not both ~½ with huge relative
   errors; a second-order treatment is not warranted before measurement says otherwise.
2. **Independence.** ρ_m and ρ_n are treated as independent. They share the same underlying transcript
   abundance, so they are positively correlated in *truth* — but their *estimation errors* come from disjoint
   fragment sets (spliced vs unspliced), which is the relevant independence. Sound.
3. **σ²_xfer is not yet measured.** Its magnitude is an empirical question — see R5, which measures the
   residual disagreement between a boundary's imputation and the adjacent region's realized density. Until
   measured, E should land with σ²_xfer carried symbolically (initially the existing s2t) so the *algebra* fix
   is separable from the *magnitude* question.

---

## 8. Implementation sketch

In `_scan`, replacing the `pr += S` pattern in the `emit_p` / `emit_n` / RNA-total blocks:

```python
rho_m = SPs[lsrc] / _esp                       # measured mature density (per-face, spliced frame)
rho_n = max(rho_pred_nascent, 0.0)             # imputed nascent density, floored (never negative)
rho_r = rho_m + rho_n
if rho_r > _EPS:
    w_m, w_n = rho_m / rho_r, rho_n / rho_r
    v_m = 1.0 / (SPNs[lsrc] + 1.0)             # MEASURED: counting only (integer spliced flux)
    v_n = 1.0 / (smn + 1.0) + v_logfr          # IMPUTED: counting + composition solve
    var_msg = w_m * w_m * v_m + w_n * w_n * v_n + s2t_xfer
    pr = 1.0 / var_msg
    mo = _mode_shift(rho_r * erd, _den, comp_fl)   # the SUM sets the mode
```

Note `rho_n` is floored at 0 **before** the sum — the absorption may say "no nascent", never "negative RNA",
and it can never cancel the measured mature.

### Gates
1. `gdna_none` phantom guard **as a delta vs HEAD** (this is what caught R2 — an absolute number is
   unreadable). Expect ≤ baseline; the §5 lower-bound guarantee should *reduce* phantom.
2. `test_mature_measurement_disagreement_silenced` — this test should now pass *for the right reason*: a
   depleted junction lowers `w_m` and raises the variance, rather than being silenced by damping.
3. `msg_audit` direction + the mature-dilution identity test.
4. Then R2 becomes a one-line removal of the enrichment-cliff proxy from this channel, re-gated.

---

## 9. Summary

| question | answer |
|---|---|
| Does the weak component get free confidence? | **No** — share-weighted, squared. Uncertain nascent dominates the variance exactly when it dominates the density. |
| Is the spliced measurement "exact" for the exon? | **No.** It is exact for the *junction*; the exon imputation carries σ²_xfer like everything else. |
| What was actually broken? | Precision **addition** — replicate-measurement algebra applied to an additive decomposition. |
| Why did it manifest as phantom gDNA? | The merged mode could be driven to ~0 by absorption while carrying measurement confidence ⇒ confident "no RNA" ⇒ f_g→1. |
| What structurally prevents that? | The measured mature is a **lower bound on RNA** once the mode is the SUM with ρ_n floored at 0. |
