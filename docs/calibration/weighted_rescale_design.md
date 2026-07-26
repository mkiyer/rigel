# The weighted rescale — derivation, prototype, and implementation plan

**Status: DERIVED + PROTOTYPED 2026-07-26. Not implemented in the solver.** Owner-directed. Read
`PASS0_FINISH_PLAN.md` §P1b first (the measurement that motivates it), then this.
Prototype: `scratchpad/wp1_prototype.py`. Instrumentation: the inert `_capture["_pin"]` hook in `bp_solver`.

---

## 1. The problem, in one paragraph and no jargon

A node observed `M` fragments of unspliced mass. A message arrives claiming "your gDNA density is `ρ_g`, your
RNA density is `ρ_R`". Those claims imply a fragment count: `ρ_g·E_g + ρ_R·E_r`, where `E_g` and `E_r` are the
effective lengths over which a gDNA and an RNA fragment respectively can be observed at this node. That number
**must** equal `M` — it is an identity, not an approximation. It does not: measured **1.24–1.31×** on the
arms where pass-0 is confidently wrong, and 1.00–1.02× on the arms where it is accurate. The size of the
violation predicts the error. Today the solver restores the identity by multiplying **every** component by
the same factor `k = M/S`, which punishes a component that is right in order to accommodate one that is wrong.

## 2. ⭐ Two structural facts established before any derivation

**(a) The violation is created by the FUSE, not by the messages.** Each message is already rescaled inside
`_transport`, so each one individually satisfies the identity. The combine then fuses the two messages
**per component** — gDNA weighted by the gDNA precisions, RNA by the RNA precisions — and a per-component
weighted average of two claims that each sum to `M` does **not** sum to `M`. **The fused claim is never
rescaled at all.** This is exactly what `_pin_v`'s own docstring already prescribes ("applied at EVERY node
rather than only at the final combine") and it is the one place it is not applied.

**(b) A COMMON rescale cannot change a composition — so today's operator is inert on the answer's split.**
`f_g = ρ_g E_g / (ρ_g E_g + ρ_R E_r)` is invariant under `ρ_c → k·ρ_c`. Confirmed in the prototype: the
`fused_common` arm is **bit-identical to shipped in every condition**. `_pin_v` moves the per-component
*levels* (which do matter — `mo_g` and `mo_p` are separate ψ factors, and the relay's next-hop reframe reads
the densities) but it can never move the gDNA-vs-RNA split. **Only a differential rescale can.**

## 3. The derivation

Write the claim's error in log space, `ρ_c = ρ_c^true·exp(ε_c)`. Choose the correction `a_c`
(`ρ_c → ρ_c·exp(a_c)`) that is most likely under the error model, subject to the identity:

```
    minimise   aᵀ Σ⁻¹ a        subject to   Σ_c m_c·exp(a_c) = M ,      m_c = ρ_c·E_c
```

Linearise the constraint about `a = 0`: with `α_c = m_c/S` (the claimed mass shares, `Σα_c = 1`) and
`δ = log(M/S)`, it becomes `αᵀa = δ`. The Lagrangian stationarity condition `2Σ⁻¹a = λα` gives `a = ½λΣα`, and
imposing the constraint fixes `λ`:

```
    a  =  δ · Σα / (αᵀΣα)                                                                    (★)
```

**Only the DIRECTION `s = Σα` comes from the error model; the magnitude is then re-solved EXACTLY** along that
direction — because the constraint is an identity and deserves to hold exactly, not to first order:

```
    find μ such that   Σ_c m_c·exp(μ·s_c) = M
```

`g(μ)` is strictly increasing whenever every `s_c ≥ 0` (true for the Σ below), so this is a unique root
reachable by bracketed bisection — **no step control, no cap, no tuned constant.**

### 3.1 The error model

Every component of a message shares the reframe `r` and the source's own fragment count, so those are
**common-mode**; whatever is left in a component's variance is its **own**:

```
    Σ = σ_cm²·11ᵀ + diag(w) ,     σ_cm² = σ²_transfer + 1/n_src ,     w_c = max(0, v_c − σ_cm²)
    ⇒   s_c = σ_cm² + α_c·w_c
```

**Two exact limits, and they are the whole justification:**

* **`w → 0`** (the only error is the shared frame): `s_c` is constant across components, so every component
  moves by the same factor. **That is exactly today's `_pin_v`** — the shipped operator is the
  zero-independent-variance limit of this one, not a different thing. Nothing is being replaced; a special
  case is being generalised.
* **`σ_cm² → 0` and `w_R → 0`** (no frame error, and the RNA arm is a measurement rather than an imputation):
  `s_R = 0`, so the RNA claim does not move at all and the gDNA absorbs the entire residual.

Capture-OFF drives the second limit (`r ≈ 1 ⇒ σ_cm² ≈ 0`), and the graft makes the RNA arm a measurement
exactly there. Capture-ON drives the first (`Var(log r)` large). **The regime switch is not a rule; it is
what the variances already say.**

### 3.2 What is preserved

`_pin_v`'s partial-claim semantics are load-bearing (rescaling all three components blindly regresses
capture-OFF 3.6×, `derive_2_relay_pin.py`) and are preserved exactly: a component the message does not
**supply** (precision 0) is filled from the node's own density for the mass budget and **does not move**
(`s_c := 0` there). A message that supplies gDNA only still delivers `f_g < 1`.

### 3.3 Why the anti-correlation term is NOT in Σ

The honest source covariance has a rank-1 composition part: a source's `f_g` error moves its two components
in *opposite* directions, with `u = (f_R, −f_g)` and `Var(λ) = 1/τ_λ` — the codebase's existing
`own_composition_logvar`. That term is deliberately **omitted**, and the reason is a result rather than a
simplification:

```
    uᵀα  =  ∂log S/∂λ  =  f_g f_R (k_g − k_R) / (f_g k_g + f_R k_R) ,     k_c = E_c^dst / E_c^src
```

so a pure composition error changes the claimed total **only** through the eff-length mismatch `k_g ≠ k_R`,
and `uᵀα → 0` when the two nodes' eff-length ratios agree. Under a strict rank-1 model the conservation
violation is therefore nearly *uninformative* about composition and (★) degenerates to the common rescale.
The observed violations (1.24–1.31×) are far larger than the eff-length bound (×1.04 contained, ×1.50 at a
boundary crossing, already on record), so **the violation is not composition drift from a single source** —
it is component-specific error, which the rank-1 model says cannot happen. It can happen because the RNA arm
receives an **independent additive measurement** (the grafted spliced count) that the gDNA arm does not.
`diag(w)` is the model that admits that; rank-1 is the model that forbids it.

## 4. Prototype results

`scratchpad/wp1_prototype.py`, applied to the solver's real fused claim. Metric: mass-weighted
`|f_g − oracle|` of the **composition the node is handed**, single-strand exons.

| condition | claim/obs | shipped | fused_common | **fused_weighted** | fused_rna_exact |
|---|---|---|---|---|---|
| `gdna300 ss0.50 present capOFF` | 1.28× | 0.2365 | 0.2365 | **0.2174** | **0.0807** |
| `gdna300 ss0.50 none capOFF` | 1.31× | 0.2364 | 0.2364 | **0.2127** | **0.0643** |
| `gdna100 ss0.50 present capOFF` | 1.24× | 0.2093 | 0.2093 | **0.1657** | **0.0612** |
| `gdna300 ss0.99 present capOFF` | 1.00× | 0.0614 | 0.0614 | **0.0579** | 0.0584 |
| `gdna1 ss0.50 present capON` | 1.02× | 0.0218 | 0.0218 | **0.0123** | **0.0056** |
| `gdna300 ss0.50 present capON` | 1.14× | 0.1556 | 0.1556 | **0.1538** | 0.1992 |
| `gdna300 ss0.99 present capON` | 1.18× | 0.1310 | 0.1310 | **0.1023** | 0.2150 |
| `gdna300 ss0.50 none capON` | 1.18× | 0.2066 | 0.2066 | **0.1759** | 0.3137 |
| `gdna100 verystrong` | 1.01× | 0.2442 | 0.2442 | 0.2456 | 0.2587 |

* **`fused_common` is identical to shipped in all nine** — §2(b), measured.
* **`fused_weighted` improves 8 of 9 and hurts none** (verystrong +0.0014, noise). Safe.
* **`fused_rna_exact` is 2.9–3.9× better on capture-OFF and low-gDNA, and much worse on gDNA-rich
  capture-ON.** It is the `w_R = 0` limit, and it is only correct where the RNA arm really is a measurement.

### 4.1 ⭐ Why `rna_exact` cannot simply be adopted — the over-claim FLIPS with capture

Mean `log(claimed / true)` per component, over the same exons:

| condition | gDNA claim | RNA claim |
|---|---|---|
| `gdna300 ss0.50 present capOFF` | **+1.59** (4.9× too big) | −0.02 (exact) |
| `gdna300 ss0.50 none capOFF` | **+1.54** (4.7× too big) | −0.01 (exact) |
| `gdna300 ss0.99 present capON` | −0.02 (exact) | **+1.31** (3.7× too big) |
| `gdna300 ss0.50 none capON` | −0.19 | **+1.73** (5.6× too big) |
| `gdna100 verystrong` | −0.67 | +0.73 |

**Which component is wrong reverses with capture**, and that is the whole reason a fixed rule fails and an
apportionment is required. It also sets the sharpest test of the variance model: does `w` reproduce the flip?
On capture-ON it does (the model already rates the RNA arm as the less certain one). **On capture-OFF it does
not** — the model still rates RNA as less certain when RNA is in fact the exact one. That is why
`fused_weighted` captures roughly a quarter of the available prize rather than all of it, and it is the single
most important open question for the implementation.

## 5. Implementation plan

**Where.** A new pure law in `enrichment_frame` (`conservation_rescale`, arrays in / arrays out, closed-form
unit-tested there as every other law is) plus its call sites in `bp_solver`. `_pin_v` becomes a thin caller of
it, so there is one operator and one set of semantics.

**Sequencing — each step its own commit and its own A/B, so a failure localises:**

| # | change | isolates | expectation |
|---|---|---|---|
| 1 | add `conservation_rescale` + unit tests; **call it nowhere** | the law and its limits | bit-identical A/B |
| 2 | apply it at the **combine**, in its `w = 0` limit (common factor) | that pinning the fused claim at all is safe | **bit-identical on composition** by §2(b); may move the per-component level channels `mo_g`/`mo_p`, so this measures exactly that side-effect |
| 3 | switch to the **weighted** direction at the combine | the thesis | the prototype says 8 of 9 better, none hurt |
| 4 | replace `_pin_v`'s body with the same law (per message) | consistency; one operator | small |
| 5 | publish `claim/obs` and the Mahalanobis `δ²/(αᵀΣα)` as a **precision** input | §6 | separate decision, separate A/B |

**Gates:** `pytest tests/ -q`, `ruff`, `message_variance_mc.py` 0 failures, the `√2·σ_own` DL inequality
untouched, per-condition A/B at refit=0 **and** refit=1, `errQ1conf`/`z2` reported alongside mwae, goldens
last.

**Unit tests to write with the law (all closed-form, no solver):**
1. `w = 0` ⇒ exactly the common factor `k = M/S` (the shipped operator is a limit).
2. `w_R = 0, σ_cm² = 0` ⇒ the RNA claim is unchanged and gDNA absorbs the whole residual, exactly.
3. the identity holds **exactly** after the rescale, to 1e-12, for random inputs across 6 orders of magnitude.
4. an unsupplied component (precision 0) does not move and still enters the mass budget from the node's own
   density — the partial-claim semantics.
5. `S = M` already ⇒ `μ = 0` ⇒ no change at all (idempotence).
6. monotonicity: a larger violation produces a larger correction, in the same direction, for every component.
7. no nan for `S = 0`, `M = 0`, a single supplied component, or an infinite variance.

**⚠ Prototype caveat to carry into the A/B.** The prototype scores the *composition* the node is handed. The
solver additionally hands ψ two per-component **level** channels (`mo_g` with `cm_g`, `mo_p`/`mo_n` with
`cm_p`/`cm_n`) which are **not** scale-invariant and therefore *do* respond to the common rescale. So step 2
is not guaranteed inert in the solver even though it is inert in the prototype — which is exactly why it is
its own step.

## 6. The second, larger opportunity this exposed

`claim/obs` — the size of the conservation violation — **predicts the error with no oracle**: 1.00–1.02× ⇒
composition error 0.02–0.06; 1.24–1.31× ⇒ 0.21–0.24. In the derivation's own terms the natural statistic is
the Mahalanobis distance `z² = δ²/(αᵀΣα)`: *how surprised should I be that this claim misses my observed mass
by this much?*

Today `_pin_v` **erases** that signal — it forces the identity and discards the residual. But it is a
DerSimonian–Laird-shaped quantity measured against a **hard observable**, and unlike M7's `b̂²` (which needs
the node's own composition evidence and is therefore inert wherever `τ_own = 0`, i.e. on all unstranded data —
83.9 % of suite error mass) **the node's own mass is always known**. This would be the first message-damping
term that is never inert. It is deliberately *not* part of the rescale change: separate derivation, separate
prototype, separate A/B.
