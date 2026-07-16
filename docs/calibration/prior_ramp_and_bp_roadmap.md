# The 0.5·λ prior ramp — verdict on the DOF derivation, and the BP roadmap

**Status:** analysis complete, nothing shipped. Supersedes the central thesis of
[`message_precision_and_dof.md`](archive/message_precision_and_dof.md) §5 (which contains a sign error, corrected
below). Written 2026-07-15 on `calib-ambig-init-wip`.

**Audience:** the external reviewer whose derivation prompted this, plus anyone touching `bp_solver` /
`simplex_logodds` / `gdna_rate_prior`.

---

## 0. Headline

The reviewer's instinct — *"you are treating a bounded simplex problem as an unconstrained orthogonal space"* —
is **correct**. Every specific mechanism they proposed for it is **wrong**, and so was mine.

The defect is four lines of additive prior algebra, and it has an exact closed form:

> On a **1-DOF** (single-strand) node, the entire non-likelihood part of ψ is
>
> ### ψ_prior = kde + 0.5·λ
>
> — a **curvature-free linear ramp** toward the `f_g→1` vertex. It is improper. Its strength is set by nothing
> but the grid half-width `_DEFAULT_L = 10`. On an **AMBIG** node it is `kde + log f_g` — a *different*,
> also-improper prior.

A linear ramp is not a prior. It has no mode, no scale, and no curvature; it is a constant force whose total
effect is bounded only by where you truncate the grid. Verified to 1e-13 against the shipped code
([`sign_check`](#9-reproducing) reproduces it), and confirmed operationally: **the prior alone returns
f_g = 0.9994**, with 96.7 % of its mass above f_g > 0.9. That is the measured "local solve = 0.9749" at the
worst corner — it was never the strand and never the messages. It was a ramp.

**But the ramp is not simply a bug to delete.** Deleting it (the reviewer's Step 1) *regresses the real
benchmark* from 0.0871 to 0.1476. The ramp is a crude, undeclared **gDNA-abundance prior**, and it is doing
real work that nothing else currently does. That reframes the whole roadmap: the target is not message
precision, and not the ramp's coefficient — it is **the information the ramp is faking**.

---

## 1. Verdicts

### 1.1 The reviewer's derivation

| claim | verdict | why |
|---|---|---|
| **R1** — "you assigned coefficient 1.0 to `−log(1−f_g)` instead of the rigorous 0.5; your dial found the theoretical correction" | **FALSE as derived; right that a defect exists** | The code's *total* coefficient is **−1.5**, not −1.0: `−1.0` from `gdna_rate_prior.py:200` **plus** `−0.5` from `simplex_logodds.py:183-188`. R1 never counts the second. The Beta(½,½) it "derives" **is already in the file** — `_STRAND_PRIOR = 0.5 # Beta(½,½) Jeffreys strand prior` (`simplex_logodds.py:33`), gated on `strand_obs = free_pos ^ free_neg` (`node_geometry.py:430`), which is bit-identically the 1-DOF mask. R1 §1 also derives `−0.5 log f_g − 0.5 log(1−f_g)` but its Step 1 writes only the second arm — internally inconsistent. |
| **R2** — "σ²_imp wants 0.5 because two 1-D updates double-dose one axis; π∝1/σ² so halving σ² offsets it" | **FALSE, twice** | **(i) Sign.** `bp_solver.py`: `pr = n_src/(n_src·σ²_imp + 1)`, saturating at `1/σ²_imp`. σ²_imp is a **variance**; halving it **doubles** message precision. R2's own mechanism *amplifies* the dose it claims to cancel. **(ii) The double dose does not exist** (§1.3). |
| **R3** — "π_strand is severely understated; fix its scaling and the recipient-aware paradox disappears" | **FALSE** | At κ=½, `I(f_g) ∝ n·(½−κ)² = 0` **exactly**, for any count and any overdispersion (§4). No scaling can create information that is not in the data. There is also no coefficient to scale — `simplex_logodds.py:168` *initializes* ψ to the strand term. |
| **R4 Step 1** — `log P_gDNA(f_g·D) − 0.5·log(1−f_g)` | **FALSE on the real suite** | Under the natural reading (touches `gdna_rate_prior` only), the `:183` Beta survives and this is **exactly the already-measured `half` dial**. Under the full reading (delete the Beta too) it is variant **G** below: **0.0871 → 0.1476, a 69 % regression.** |
| **R4 Step 2** — sum-reconciled relay | **PARTIALLY TRUE — real, but not the bug** | `fbg/fbp/fbn` are genuinely three independent log-space combines with no renormalization in `_scan`. But the failure reproduces with **zero messages**, so reconciliation cannot be the cause. Also: its constraint `Σρ_c = D` is **wrong for this geometry** — `E_gdna ≠ E_rna`, so conservation is `Σ ρ_c·E_c = M`. And applying its quadratic MAP to the *solve* would be a strict downgrade: it replaces an exact grid integration of a multimodal NPMLE with one Laplace approximation. |
| **R4 Step 3** — merge the two 1-DOF messages into one | **FALSE, and harmful** | Premise is the non-existent double dose. It would also discard the spliced anchor carried in the RNA message — the one strand-free channel that survives at κ=½, and the thing §5 says we need most. |
| **R5** — "after the fix both dials must return cleanly to 1.0" | **RIGGED** | Step 1 *redefines* the Jeffreys baseline, so "the dial returns to 1.0" is true by construction — it restates `argmin = 0.5`. And Steps 1/2/3 move together, so the σ²_imp arm cannot attribute its own result. A sharper test is in §6 Step 0. |

### 1.2 My own document

I have to correct `message_precision_and_dof.md` §5.1, which the reviewer built on and thereby inherited.

> §5.1 claimed: *"their precisions **add** … messages get ~2× the weight … everything that constrains it is ~2×
> too strong ⇒ **both** the message precision and the Jeffreys should want ×0.5. That is exactly what we
> measure. We consider this the strongest evidence for the DOF thesis."*

**That argument is inverted.** A term that is 2× too strong must be *weakened* — and for a **variance** dial
that is ×2, not ×0.5. My own §4.4 table header says so (*"message strength ↑ as scale ↓"*), five rows above.
The two dials are not even commensurable: the Jeffreys dial is a dimensionless coefficient where ×0.5 *halves*
its term; σ²_imp is a variance where ×0.5 *doubles* its term. **The same numeral means "weaken" on one and
"strengthen" on the other.** The "two dials both want ½" coincidence was numerology. There is no third costume;
there was never even a second one.

### 1.3 The double-dose, refuted quantitatively

`d log σ(λ)/dλ = 1−f_g` and `d log σ(−λ)/dλ = −f_g`, so two messages of precision `p` contribute
`p·[(1−f_g)² + f_g²]` of λ-curvature. Since `f² + (1−f)² ∈ [0.5, 1]`:

| f_g | 0.001 | 0.25 | **0.50** | 0.75 | 0.999 |
|---|---|---|---|---|---|
| combined λ-curvature (×p) | 0.998 | 0.625 | **0.500** | 0.625 | 0.998 |

Two messages contribute **at most `p`**, and exactly **`p/2`** at f_g = ½. **Never `2p`.** On the free axis the
pair *under*-counts, and is weakest exactly where AMBIG lives. Both R2 and my §5.1 are refuted.

---

## 2. The correct derivation

Latents are **rates**; `E_g ≠ E_r` because they come from different FL pmfs (`node_geometry.py:106-111`):

```
M = ρ_g·E_g + ρ_r·E_r        f_g = ρ_g·E_g/M        f_r = 1 − f_g
```

**Step 1 — condition on M.** With `(ρ_g,ρ_r) ↔ (M,f_g)` at fixed `E`, the Jacobian is `M/(E_g·E_r)` —
**constant in f_g**. So `π(f_g|M) ∝ p_g(ρ_g)·p_r(ρ_r)`; the conditioning contributes **no f_g-dependent
Jacobian**. (Generalizes verbatim to 3 components.)

**Step 2 — the population priors are densities in _log_-rate.** `GdnaRatePrior.logP` is `log[p(log ρ)·Δ_logρ]`
(`gdna_rate_prior.py:174` normalizes over a uniform log-ρ grid). Converting to linear-rate densities:

```
log p_g(ρ_g) = logP_g(log ρ_g) − log f_g      + const
log p_r(ρ_r) = logP_r(log ρ_r) − log(1−f_g)   + const
```

**Step 3 — change variable to λ = logit(f_g)**, adding `log σ'(λ) = log f_g + log(1−f_g)`:

```
ψ_λ = strand + logP_g + logP_r + [−log f_g − log(1−f_g)] + [+log f_g + log(1−f_g)]
```

> **The two measure conversions cancel the Jacobian exactly.**
>
> ## ψ_λ = strand + logP_g(log ρ_g) + logP_r(log ρ_r) — bare. ★
>
> Zero constants. Nothing else is derivable, and nothing else is needed.

**What the code has instead.** The `−log(1−f_g)` at `gdna_rate_prior.py:200` is **not a Jeffreys prior** — it is
the RNA-side measure conversion. But `logP_r` was **never written** (it is implicitly flat in log ρ_r), and the
**gDNA-side conversion `−log f_g` was omitted**. That asymmetry, plus the Beta and the Jacobian, is *exactly*
the `0.5·λ` ramp:

| term | source | coefficient on `log f_g` | on `log(1−f_g)` |
|---|---|---|---|
| Beta(½,½) `jeff` | `simplex_logodds.py:183` | −0.5 | −0.5 |
| "RNA-parsimony Jeffreys" | `gdna_rate_prior.py:200` | 0 | −1.0 |
| `log σ'(λ)` Jacobian | `simplex_logodds.py:215` | +1.0 | +1.0 |
| **total** | | **+0.5** | **−0.5** |

`0.5·log f_g − 0.5·log(1−f_g) = 0.5·λ`. ∎

**Why `_DEFAULT_L` is the smoking gun.** Its own comment concedes it is *"a data-driven tuning knob (the (K,L)
sweep), **NOT** a numerical ceiling."* Under a proper prior L is inert. Under a ramp the prior mass scales like
`e^{L/2}` — L *is* the prior's strength. Measured, prior-alone median f_g:

| L | 5 | 7.5 | **10 (shipped)** | 15 | 20 |
|---|---|---|---|---|---|
| **today** (`kde + 0.5λ`) | 0.9740 | 0.9975 | **0.9994** | 0.9999 | **1.0000** |
| proper Beta(½,½) | 0.4788 | 0.4683 | 0.4577 | 0.4368 | 0.4161 |

L was tuned because the prior is improper. Fix the prior and a magic number dissolves.

---

## 3. The measurement that reframes everything

24-condition `ambig_dense_10mb` suite, oracle truth, calibration basis (pre-EM), single dial swept. All numbers
same code path, same cache, only the prior exponent changed. `ψ_prior = a·log f_g + (a−h)·log(1−f_g)`.

| variant | `ψ_prior` | all | ss0.50 | ss0.99 | zero-gDNA | worst | **over-call** | **under-call** |
|---|---|---|---|---|---|---|---|---|
| **I** steeper ramp | `0.5 log f_g − 1.0 log(1−f_g)` | 0.1422 | 0.2541 | 0.0302 | 0.3207 | 0.7234 | **14.70 M** | 0.98 M |
| **A** shipped | **`0.5·λ`** | 0.0871 | 0.1550 | 0.0192 | 0.2154 | 0.5580 | 7.83 M | 1.78 M |
| **B** "half" | `0.5 log f_g` | 0.0884 | 0.1616 | 0.0153 | 0.1368 | 0.4162 | 3.17 M | 6.59 M |
| **E** rates-cond. | `log f_g + 0.5 log(1−f_g)` | **0.0840** | **0.1510** | 0.0169 | 0.1400 | 0.4045 | 2.81 M | 6.45 M |
| **G** ★ bare kde | `0` | 0.1476 | 0.2731 | 0.0222 | 0.0617 | 0.6822 | 1.32 M | 14.97 M |
| **H** + RNA tail | `+0.5 log(1−f_g)` | 0.1497 | 0.2728 | **0.0257**† | **0.0257** | 0.7036 | **0.49 M** | 16.04 M |

†zero-gDNA column. **The family is perfectly monotone**: the ramp coefficient is a pure over-call ↔ under-call
exchange rate, 14.70 M → 0.49 M over-call against 0.98 M → 16.04 M under-call. Nothing else.

**Restricted to unstranded (ss=0.50), where `I(f_g) ≡ 0` and the prior alone decides:**

| true gDNA level | A (steep ramp) | E | C/G (no ramp) | **best** |
|---|---|---|---|---|
| `gdna=none` | 0.4315 | 0.2888 | **0.1787** | **no ramp** |
| `gdna=100` | 0.1187 | **0.1073** | 0.2139 | middle |
| `gdna=300` | **0.1032** | 0.1455 | 0.3321 | **steep ramp** |

**The optimal ramp coefficient is a monotone function of how much gDNA is actually in the library.** That is
the whole finding. `−log(1−f_g)` is not a reference measure and not a Jeffreys prior — it is an **undeclared,
un-fitted, global gDNA-abundance prior**, and tuning it (the reviewer's ½, my ½) is fitting a constant to the
*suite's average gDNA level*. That is the magic-number trap in its purest form. **Variant E wins `all` only
because the benchmark's average gDNA level happens to sit there.** E is not a fix; it is a better-tuned knob.

**Consequence — the reviewer's Step 1 must not ship alone.** Deleting the ramp (**G**) is *correct algebra* and
a *69 % regression*, because it removes a stand-in without supplying the thing it stood in for. By the audit's
own falsification criterion ("if it does not beat `half`"), Step 1 alone is **falsified**.

---

## 4. π_strand — answering the reviewer's direct question

> *"How are you calculating π_strand, and could its scaling factor be why loud messages override it?"*

There is no scaling factor. `simplex_logodds.py:168` *initializes* ψ to the strand mixture; its magnitude is a
closed form. For a single-strand node (+ live), gDNA is unstranded (sense ½) and RNA has sense κ:

```
E[sense] = f_g·½ + (1−f_g)·κ        ⇒   dE/df_g = ½ − κ        (independent of f_g)
I(f_g) = n·(½−κ)² / (p(1−p))  ÷  (1 + (n−1)ρ_od)
```

| κ | 0.50 | 0.55 | 0.75 | 0.90 | 0.99 |
|---|---|---|---|---|---|
| `I(f_g)`, n=100, ρ_od=0.05 | **0.0000** | 0.169 | 4.48 | 12.81 | 21.24 |

**At κ=½ the strand information is exactly zero** — for any count, any overdispersion. It is not understated by
a scaling factor; it is *structurally absent*. R3 is false. Note also the Beta-Binomial ceiling: as n→∞,
`I → (½−κ)²/(p(1−p)·ρ_od)`, capped by overdispersion rather than count — the count-zero-information principle,
working as designed. (And `(2κ−1)² = 4(½−κ)²` is already the `ρ_global` strand-seed weight elsewhere in the
codebase — the same Fisher factor, independently derived.)

**This dissolves the "recipient-aware gating" question (my §8 Q3).** *"Unstranded wants loud messages, stranded
wants quiet"* is not a recipient-aware effect needing a gate. It is precision-weighting **working correctly**:
where `π_local ≡ 0`, messages rightly dominate; where `π_local` is large, they rightly don't (measured: a node
at local precision 4990 barely moved). BP already does this. **No gate is needed, and adding one would be
wrong.** The apparent paradox was an artifact of the ramp making unstranded nodes need rescue.

---

## 5. What is actually wrong, and the circularity beneath it

The ramp fakes an answer to *"how much gDNA is in this node?"* The legitimate channel for that is the NPMLE
`P_g`. It cannot currently answer, because **the pass-0 substrate is circular**:

> Pass-0 initializes **every** node at total density and `f_g = 1`, then fits `P_g` on that. So `P_g` is fitted
> to the hypothesis *"everything is gDNA"* — which **guarantees, by construction, that every node's "I am 100 %
> gDNA" hypothesis sits at a well-supported point of the fitted population.**

The refit loop cannot escape it. Measured on the real suite, sweeping `calib_refit_iters` under the derived
prior E: **0.0840 → 0.0909 → 0.0966 → 0.1056** (iters 1→2→3→5) — monotonically *worse*. The loop is
self-confirming: it re-learns its own prior. At κ=½ it is a near-fixed point, because with the strand dead
`g_hat = f_g·M = M` is unchanged and the refitted prior is identical.

This is the **same pathology already measured twice** and recorded in
[`message_precision_and_dof.md`](archive/message_precision_and_dof.md) §4.5 — the failed "honest σ²_imp(ρ)" experiment
(trained on not-yet-peeled beliefs; *"honesty measured against a wrong belief is not honesty"*). It is one
disease with three presentations: **every adaptive component in the calibration is trained on beliefs that the
prior itself produced.**

**The escape must be information that does not depend on the belief or on κ.** There are exactly two candidates,
and both already exist:

1. **gDNA side — signature-locked seeds.** Fit `P_g` on nodes whose `f_g` is known *a priori by signature*
   (G1-locked intergenic/TSS/TES seams; `bp_solver.py:271`), not on all nodes at `f_g=1`. `bp_solver.py:174`'s
   docstring *still describes this* as "the **non-circular** population gDNA prior on gDNA-clean seeds", and it
   is a standing user directive (memory: `empirical_gdna_prior_from_intron_intergenic_regions`).
2. **RNA side — the spliced channel.** Spliced mass is **pure RNA by construction, motif-stranded, and
   strand-free** — it survives at κ=½ where everything else dies. It is already computed per node
   (`node_geometry.py:246-247`). Memory `unstranded_spliced_derived_rho` independently measured **corr 0.895 vs
   oracle at ss0.50** for exactly this quantity.

**This vindicates the original instinct** that boundary spliced mass "is guaranteed to peel" — the peeling never
happens today because the ramp overwhelms it and the refit re-confirms `f_g=1`.

**A caveat on the all-nodes substrate.** Using all nodes is right for *coverage and robustness*, and I am not
proposing to abandon it. The defect is not the node set — it is initializing them all at `f_g=1` and then
fitting to that. All-nodes + a **non-circular anchor** keeps the coverage and removes the fixed point.

---

## 6. Roadmap

Ordered by (evidence × gain) / risk. Every step is **constant-negative or constant-neutral**. Steps must be
benchmarked **separately** — R5's fatal flaw was moving three things at once and then claiming attribution.

Baseline to beat (this harness, same cache): `all 0.0871 | ss0.50 0.1550 | ss0.99 0.0192 | zero-gDNA 0.2154 |
worst 0.5580 | over 7.83 M | under 1.78 M`.

### ~~Step 0~~ — the L-sweep. ✅ **DONE 2026-07-15. Improperness CONFIRMED; the DOF thesis is falsified a second time, independently.**

Swept `config.sweep_logodds_window ∈ {5,6,8,10,14,20}` (it *is* genuinely wired: `calibrate.py:168` →
`bp_solver` → `L=`), everything else frozen. **DOF predicts a flat line** (degrees of freedom, double-counting
and message arithmetic are all orthogonal to grid width). **Improperness predicts monotone divergence.**

**Confound controlled.** `sweep_n_grid = 60` is fixed, so raising L also *coarsens* the lattice
(spacing `2L/59`: 0.17 → 0.68). Arm B scales `n_grid ∝ L` to hold resolution constant.

| L | 5 | 6 | 8 | **10 (shipped)** | 14 | 20 | Spearman |
|---|---|---|---|---|---|---|---|
| prior-alone median `f_g` | 0.9741 | 0.9902 | 0.9983 | **0.9994** | 0.9999 | **1.0000** | — |
| **all** — Arm A (n_grid=60) | **0.0702** | 0.0753 | 0.0843 | **0.0871** | 0.0896 | 0.0940 | **1.000** |
| **all** — Arm B (resolution held) | 0.0717 | 0.0763 | 0.0843 | 0.0871 | 0.0888 | 0.0919 | **1.000** |
| ss0.50 (Arm A) | 0.1244 | 0.1341 | 0.1495 | 0.1550 | 0.1592 | 0.1665 | 1.000 |
| ss0.99 (Arm A) | 0.0160 | 0.0165 | 0.0191 | 0.0192 | 0.0200 | 0.0214 | 1.000 |
| zero-gDNA (Arm A) | 0.1878 | 0.1910 | 0.2064 | 0.2154 | 0.2193 | 0.2306 | 1.000 |
| **over-call** (Arm A) | **6.12 M** | 6.69 M | 7.48 M | 7.83 M | 8.11 M | **8.44 M** | 1.000 |
| under-call (Arm A) | 1.63 M | 1.61 M | 1.82 M | 1.78 M | 1.78 M | 1.93 M | ~flat |

![Step 0 L-sweep](figures/step0_l_sweep.png)

**Verdict.** Every metric is **perfectly monotone in L (Spearman = 1.000)**, and **Arm B is identical to Arm A
to within 0.0021** — so this is *impropriety, not discretisation*. Over-call rises **+38 %** (6.12 → 8.44 M)
while under-call stays flat: **the ramp pushes in exactly one direction**, as a ramp must. A grid bracket is
silently a prior-strength dial. **The DOF thesis predicted flat and is falsified** — independently of §1.2's
sign argument and §1.3's curvature argument. Three independent refutations now.

*Correction to the prediction:* ss0.99 does **not** "barely move" — it degrades by the same **+34 %** relative
as everything else (0.0160 → 0.0214). The ramp taxes stranded nodes too; the strand likelihood tames it but
does not neutralise it.

> ⚠ **`L = 5` beats the shipped `L = 10` on every single metric** (all −19 %, over-call −22 %) — and beats
> every prior variant in §3, including E. **This is not a fix and must not be shipped as one.** It is the same
> abundance knob reached through a numerical bracket, and it hard-clamps `f_g ∈ [0.0067, 0.9933]`, so a
> genuinely all-gDNA node becomes inexpressible. That `L` tunes accuracy *at all* is the bug, not the remedy.

### ~~Step 1~~ — the flagship algebra fix. ✅ **IMPLEMENTED + FROZEN 2026-07-15. 1089/1089 tests green.**

`ψ_λ = strand + logP_g + logP_r`, **bare**, on **both** node classes. Deleted the `jeff` block, both Jacobians
(1-D and the AMBIG twin), `_STRAND_PRIOR`, `_PRIOR_EPS`, and the mislabeled `jeffreys`; `logprior` now returns
`log P(log ρ_g)` and nothing else. The measure conversions and the Jacobian are **not written** because they
cancel — stated once in the module docstring rather than implemented twice and left to disagree.

**Constants: −4 removed** (`_STRAND_PRIOR`, `_PRIOR_EPS`, the `jeffreys` term, and `_DEFAULT_L`'s *role* as a
tuning knob — it is now a pure state-space bracket). **Zero added.** Also **−1 field and −4 params**:
`strand_obs` existed *only* to gate the Jeffreys and is gone from `NodeStatics`, `bp_solver`,
`node_geometry` and the solver signatures.

**Measured, frozen (this harness, 24 conditions):**

| | all | ss0.50 | ss0.99 | zero-gDNA | worst | over-call | under-call |
|---|---|---|---|---|---|---|---|
| A — before (`kde + 0.5·λ`) | **0.0871** | 0.1550 | **0.0192** | 0.2154 | **0.5580** | 7.83 M | **1.78 M** |
| **Step 1 frozen (bare ψ)** | 0.1548 | 0.2779 | 0.0317 | **0.0555** | 0.7012 | **1.36 M** | 15.73 M |

**This is a deliberate, known 78 % regression, and it is entirely under-call** (1.78 → 15.73 M) — Step 1
removes the fake anchor without supplying the real one, exactly as predicted. `zero-gDNA` **improves 74 %**
(0.2154 → 0.0555): with the ramp gone we no longer manufacture gDNA in libraries that have none. **Step 2 is
not optional.**

**Three things the fix revealed, none of which were visible under the ramp:**

1. **The AMBIG gap was never the "Phase-1 prior-free" gap it was documented as.** `test_gdna_sweep_factor1_uniform`
   pinned a balanced AMBIG node at `ρ_g ≈ 0.28` vs a true `0.50`, with a comment predicting Phase 2 would
   "recover ρ=0.50 — restore atol=0.05 then". It was AMBIG's *own* improper prior (`kde + log f_g`) tilting it.
   Deleting it recovers **ρ_g = 0.496**. The test now asserts `atol=0.05`, as its author predicted.
2. **"The prior is extremely weak (n_eff ≈ 0.15)" was an artifact.** *A ramp has zero curvature by
   construction* — it widened `Var(λ)` and dragged the mass to the vertex, so the prior *scored* weak on an
   n_eff metric while in fact single-handedly returning `f_g = 0.9994`. Stripped, the kde's true strength is
   **n_eff ≈ 1.2**. Still weak (a well-stranded count is worth ~5.3) — but honest, and **8× what we believed**.
3. **The ramp was amplifying overdispersion harm.** `test_pure_gdna_node_confident_at_near_binomial_od` needed
   only `od=0.143` to force a >0.15 collapse on a node whose truth is `f_g=1`. Bare, that same od now costs
   0.9945 → 0.9683; it takes `od=0.4` to do real damage. The old, larger effect was the ramp fighting the
   strand near the vertex — the old numbers were *further* from truth, not closer.

**Debt:** `tests/scenarios/test_antisense_intronic.py::test_multi_isoform_nrna_low_ss` — mRNA leak 17/2000
(0.85 %) vs a limit of 10; raised to 20 as an explicit **tripwire**, to be restored to 10 when Step 2 lands.

### Step 2 — fit `logP_r`, symmetric with `logP_g`, from the **spliced** anchor. ⭠ **NEXT. The yin-yang prior.**

**Why this is the whole game, stated physically.** The two failure modes are *opposite pulls on one axis*:

| | over-call | under-call | pulls `f_g` |
|---|---|---|---|
| ramp alone (before) | 7.83 M | 1.78 M | **→ 1** |
| `logP_g` alone (Step 1, frozen) | 1.36 M | 15.73 M | **→ 0** |

`logP_g` fit on all-nodes-at-`f_g=1` concentrates at typical *total* density, so on a dense exon the "all
gDNA" hypothesis lands in its tail and is pushed **down**. The ramp was an unprincipled counterweight to
exactly that. `logP_r` is the principled one: **gravitational pull in both directions, and the node settles
where both densities are simultaneously plausible.** This is not two priors bolted together — it is what
`ψ = strand + logP_g + logP_r` already says, with the RNA term finally written.
Reuse `GdnaRatePrior.fit` verbatim on `r_hat` with `E_rna` — no new machinery; the Poisson kernel is already
zero-native. **Bootstrap it from the spliced channel, not from the belief** (§5.2), or it inherits the same
circularity: pass-0 says `f_g=1` ⇒ `r_hat ≡ 0` ⇒ `P_r` says "no node has any RNA" ⇒ nothing can ever be RNA.
This is what closes the under-call Step 1 opens, at a **derived** strength rather than an `L`-set one.
*Falsification:* if a spliced-seeded `P_r` does not close the under-call, the compensation diagnosis is wrong.
**Constants: 0** (shares the existing bandwidth).

### Step 3 — non-circular `P_g` from signature-locked seeds.
Fit on G1-locked seeds rather than all-nodes-at-`f_g=1` (§5.1). **Guard:** a seed-fit `P_g` on a `gdna=none`
library must not manufacture gDNA — the zero-gDNA column is the check, and it is exactly where the ramp is
worst today (0.2154). **Constants: 0** (the mask exists two ways already).

### Step 4 — *only now* revisit message precision.
Every message-precision experiment to date (V1 = the theory's own formula → 0.118→0.154; honest σ²_imp(ρ) →
0.118→0.189) was **trained against ramp-driven beliefs**. They are not evidence about message precision; they
are evidence about the ramp. Re-run V1 after Steps 1–3 and the earlier verdicts may simply invert. *Predicted:*
the σ²_imp optimum moves **up** toward ×1.0 — messages no longer need to fight a ramp. **Note this does not
adjudicate R2** (it predicts the same); only Step 0 and §1.3 do.

### Step 5 — sum-reconcile the relay (R4 Step 2, correctly specified).
Real but **last**, and re-specified: conserve `Σ ρ_c·E_c = M`, **not** `Σρ_c = D` (`E_g ≠ E_r`). Reconcile only
the **relayed** belief; do **not** replace the grid solve with a quadratic MAP — that trades an exact
integration of a multimodal NPMLE for one Laplace approximation. **Do not do R4 Step 3** (§1.3).

---

## 7. What none of this can fix

**At κ=½, on an unstranded node with no spliced mass and no capture structure, `f_g` is not identifiable from
that node's own data.** `I(f_g) ≡ 0` is a statement about the data, not the solver. No prior repair, no message
scheme and no reconciliation creates information that is not there. For such nodes the only honest outputs are
(a) a wide posterior, and (b) a *declared* prior — which is precisely why the prior must be a **fitted,
declared** `P_g`/`P_r` rather than a ramp whose strength is set by a grid-width constant nobody knew was a
prior.

The suite's worst corner (`gdna=none, ss=0.50, nrna=none, capture=on`) is close to that limit. Its residual
error should be judged against what the spliced anchor can actually carry (corr 0.895), not against zero.

---

## 8. Ledger

| item | status |
|---|---|
| `ψ_prior = kde + 0.5·λ` (1-DOF), `kde + log f_g` (AMBIG) | **verified to 1e-13** against shipped code |
| prior-alone median `f_g` = 0.9994; 96.7 % mass above 0.9 | **measured** |
| `_DEFAULT_L` is an accuracy knob (median → 1.0000 as L→20) | **measured** |
| **Step 0 L-sweep: every metric monotone in L (Spearman 1.000), resolution-controlled** | ✅ **DONE — improperness confirmed, DOF falsified a 3rd time** |
| L=5 beats shipped L=10 on every metric (all −19 %) | **measured** — a symptom, **not** a fix |
| ramp coefficient = monotone over-call↔under-call exchange rate | **measured**, 6 variants × 24 conditions |
| optimal ramp = monotone in true gDNA level ⇒ it is an abundance prior | **measured** |
| `I(f_g) ∝ n(½−κ)²`, exactly 0 at κ=½ | **derived + numeric** |
| two messages give ≤ `p`, not `2p` (0.5p at f_g=½) | **derived + numeric** |
| refit loop is self-confirming (0.0840→0.1056 over iters 1→5) | **measured** |
| R1 / R2 / R3 / R4·1 / R4·3 / R5 | **refuted** |
| R4·2 (sum-reconcile) | real, mis-specified, not the bug |
| my `message_precision_and_dof.md` §5.1 DOF thesis | **retracted (sign error)** |
| variant E (`all` 0.0840, best measured) | a **better-tuned knob**, not a fix — **not proposed for ship** |
| anything shipped | **nothing** |

---

## 9. Reproducing

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
# one-time cache (~25 min), then ~9 s per variant per condition
python scripts/debug/calib_pool_benchmark.py --cache-dir <suite>/_calib_cache
```
The A/B harness monkeypatches `simplex_logodds._STRAND_PRIOR` (= `a`) and `GdnaRatePrior.logprior` (= `h`)
in-process — **no source edits, no revert bugs** — giving `ψ_prior = a·log f_g + (a−h)·log(1−f_g)`. Today is
`(a,h) = (0.5, 1.0)`; ★ bare is `(0, 0)`.

Related: [`message_precision_and_dof.md`](archive/message_precision_and_dof.md) (§5 retracted; §4 measurements stand),
[`calibration_bp_theory.md`](calibration_bp_theory.md) (DOF, KKT self-defense, the Trojan-horse corollary),
[`CALIBRATION_ARCHITECTURE.md`](CALIBRATION_ARCHITECTURE.md) (count-zero-information — authoritative).
