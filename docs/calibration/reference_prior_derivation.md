# The per-node reference prior — derivation, open questions, and what is *chosen* vs *forced*

**Status:** **REVIEWED 2026-07-15 — the reference is RESOLVED. See §10 for the outcome and the adopted form.**
Written 2026-07-15, branch `calib-ambig-init-wip`.
**Audience:** a reviewer who does not know the codebase. Everything needed is stated here.
**Companion:** [`calibration_roadmap.md`](calibration_roadmap.md) §1 states an earlier version of this
derivation; this document supersedes its §3 "OPEN — the 2-DOF reference is NOT settled" item.

> **§1–§9 are the document as submitted for review** (the questions, and our reasoning up to that point).
> **§10 records the review's answers, our verification of them, and the resolved design.** Where §10
> contradicts §1–§9, §10 wins — most importantly, §5.3's `Dir(½,¼,¼)` proposal is **superseded** by the
> exact Berger–Bernardo prior, and §3's framing of the measure residual `R` as separately-shippable is
> **corrected** (it is non-integrable alone).

---

## 1. The problem in one paragraph

We must estimate, at each of ~10⁶ genomic nodes, the composition of that node's observed fragment mass into
**gDNA** and **RNA** (and, at nodes where both genomic strands may carry transcripts, into RNA₊ vs RNA₋).
The estimator integrates a per-node log-posterior `ψ` on a fixed numerical grid. We have discovered that the
integrand as currently written omits every component's prior/reference term, and that **omission is not
neutral**: the grid's coordinate supplies a reference by default. We need to know what reference we *should*
write. This document derives the measure-theoretic bookkeeping exactly (§3–§4), and then separates what that
bookkeeping *forces* from what remains a modelling *choice* (§5–§6).

---

## 2. Setup and notation

### 2.1 Components and rates

A node has `n` **components** indexed by `c`. Each component `c` has

* a latent **rate** (density) `ρ_c > 0` — physically, fragments per effective base,
* a known **effective length** `E_c > 0` (a fixed, precomputed constant per node and component),
* an expected count `ρ_c · E_c`.

Two node classes occur:

| class | components | `n` |
|---|---|---|
| **single-strand** ("1-DOF") | gDNA `g`, RNA on the one live strand `r` | 2 |
| **AMBIG** ("2-DOF") | gDNA `g`, RNA₊ `p`, RNA₋ `m` | 3 |

### 2.2 Composition coordinates

Let `M` be the node's total observed mass and define the **composition**

```
f_c  =  ρ_c E_c / Σ_c' ρ_c' E_c'          (so f_c > 0, Σ_c f_c = 1)
```

We reparameterize the simplex by a **log-odds** coordinate `λ` and (for `n=3`) a **tilt** `τ`:

```
f_g   = σ(λ)              σ(λ) = 1/(1+e^{-λ}),   λ ∈ ℝ
n=2:  f_r   = 1 − f_g
n=3:  f_p   = (1 − f_g)(1 + τ)/2,   f_m = (1 − f_g)(1 − τ)/2,   τ ∈ [−1, 1]
```

`λ` is used because `ρ_g` spans 5–6 decades and `λ` resolves both vertices (`f_g → 0` and `f_g → 1`), which a
uniform grid in `f_g` cannot (see §6.3 for a numerical demonstration).

The estimator integrates `exp(ψ)` on a **uniform grid in `λ`** over `[−L, L]` (and, for `n=3`, a uniform grid
in `τ` over `[−1,1]`). `L` is intended purely as a **state-space bracket** — the widest `f_g` representable —
and *not* as a tunable parameter. **This intent is a testable claim: a correctly-specified `ψ` must be
`L`-invariant.** It is the main falsification test we use throughout.

### 2.3 The likelihood

`ψ = (log-likelihood) + (log-prior)`. The log-likelihood contributes two things:

1. a **strand likelihood** — a two-component Beta-Binomial on the node's per-strand counts `(u₊, u₋)`,
   parameterized by a library-level sense fraction `κ ∈ [½, 1]`. gDNA is symmetric (mean ½); RNA is tilted to
   mean `κ`. **Its Fisher information for `f_g` is `I(f_g) ∝ n(½ − κ)²`** — *exactly zero when `κ = ½`
   (unstranded libraries), for any count, any overdispersion*. This is structural. §7 Q1 turns on it.
2. **imputation messages** from adjacent nodes (Gaussians in `log f_c`). Not relevant to this document — they
   are ordinary likelihood factors and do not affect the measure bookkeeping.

**We are not asking about the likelihood.** We are asking what prior/reference to write for each component.

### 2.4 What we have and don't have

* `logP_g` — a **fitted** population prior on the gDNA log-rate, `P_g(log ρ_g)`, estimated by a
  nonparametric MLE (a fixed-kernel Poisson mixture) over all nodes. **A density in log-rate.** We have this.
* `logP_r`, `logP_p`, `logP_m` — the analogous RNA-side priors. **We do not have these.** They are unwritten.

The question is what to write in place of the missing ones.

---

## 3. The measure identity (the part that is forced)

> **Proposition.** Let each component's prior be given as a density `P_c` in its **log**-rate, i.e.
> `p_c(ρ_c) dρ_c = P_c(log ρ_c) d log ρ_c`, so `p_c(ρ_c) = P_c(log ρ_c) / ρ_c`. Assume the components are a
> priori independent: `p(ρ_1..ρ_n) = ∏_c p_c(ρ_c)`. Then, as a density in the coordinates `(λ, [τ,] M)`,
>
> ```
>     log p(λ, [τ,] M)  =  Σ_c logP_c(log ρ_c)  +  R  +  (terms in M alone)
> ```
> where the **measure residual** `R` depends only on the coordinates and *not* on the priors:
> ```
>     n = 2 :   R  =  0
>     n = 3 :   R  =  − log( (1 − τ²)/4 )  −  log 2
> ```

### 3.1 Proof

Write `ρ_c = f_c · M / E_c`. Split the change of variables into two steps.

**Step 1 — rates → (composition chart, M).** Use the chart `(f_1, …, f_{n−1}, M)` (the last `f_n` is
determined). The Jacobian is

```
    ∂(ρ_1..ρ_n) / ∂(f_1..f_{n−1}, M)   =   M^{n−1} / ∏_c E_c
```

*Verification for `n=2`:* rows `(∂ρ_g/∂f_g, ∂ρ_g/∂M) = (M/E_g, f_g/E_g)` and
`(∂ρ_r/∂f_g, ∂ρ_r/∂M) = (−M/E_r, (1−f_g)/E_r)`; determinant
`= M(1−f_g)/(E_gE_r) + M f_g/(E_gE_r) = M/(E_gE_r)`. ✓
*Verification for `n=3`:* expanding the 3×3 determinant gives `M²/(E_gE_pE_m) · (f_g + f_p + f_m) = M²/(E_gE_pE_m)`. ✓
In both cases **composition-independent** — this is the key structural fact.

Now assemble. Using `p_c(ρ_c) = P_c(log ρ_c)/ρ_c` and `∏_c ρ_c^{−1} = M^{−n} ∏_c f_c^{−1} ∏_c E_c`:

```
 p(f, M)  =  [∏_c P_c(log ρ_c)] · M^{−n} ∏_c f_c^{−1} ∏_c E_c  ·  M^{n−1}/∏_c E_c
          =  [∏_c P_c(log ρ_c)] · M^{−1} · ∏_c f_c^{−1}
```
so
```
    log p(f, M)  =  Σ_c logP_c(log ρ_c)  −  Σ_c log f_c  −  log M                     (★)
```
Every `E_c` has cancelled, as has all but one power of `M`.

**Step 2 — composition chart → `(λ[, τ])`.** Let `J₂ = |∂(f-chart)/∂(λ[,τ])|`.

```
 n=2:  J₂ = dσ/dλ = f_g(1 − f_g)                              ⇒  log J₂ = log f_g + log(1−f_g)
 n=3:  chart (f_g, f_p);   ∂f_g/∂λ = f_g(1−f_g),  ∂f_g/∂τ = 0
                            ∂f_p/∂λ = −f_g(1−f_g)(1+τ)/2,  ∂f_p/∂τ = (1−f_g)/2
       ⇒ J₂ = f_g(1−f_g) · (1−f_g)/2 = f_g(1−f_g)²/2
       ⇒ log J₂ = log f_g + 2 log(1−f_g) − log 2
```

Adding `log J₂` to (★) and defining `R := log J₂ − Σ_c log f_c` gives the Proposition. Evaluate `R`:

```
 n=2:  Σ_c log f_c = log f_g + log(1−f_g) = log J₂       ⇒  R = 0                                    ∎
 n=3:  Σ_c log f_c = log f_g + log f_p + log f_m
                   = log f_g + 2 log(1−f_g) + log((1−τ²)/4)
       ⇒  R = [log f_g + 2log(1−f_g) − log 2] − [log f_g + 2log(1−f_g) + log((1−τ²)/4)]
          = − log((1−τ²)/4) − log 2                                                                  ∎
```

### 3.2 Two consequences we want checked

**(a) For `n = 2` the measure terms cancel exactly.** This is the codebase's current claim and it is
**correct**: with both components' priors supplied as log-rate densities, no Jacobian and no measure term
should be written. The cancellation requires *both* components' `−log ρ_c` conversions; it is not a property
of the gDNA arm alone.

**(b) For `n = 3` they do *not* cancel.** A residual `R = −log((1−τ²)/4)` survives. It is **independent of
every prior** — pure measure. It is **not currently written anywhere in the code**, and it diverges (though
integrably, `∫₋₁¹ (1−τ²)^{−1/2} dτ = π`) as `τ → ±1`. Its effect is to *reward* strand-dominant tilts.

We would like (b) confirmed, since it means the AMBIG path has been missing a derived, prior-independent term
regardless of how the reference question is settled.

> **⚠ Post-review correction — see §10.2.** (b) was confirmed, but this framing is misleading. `R` is **not
> separately shippable**: written without a taming RNA-side reference the `τ` density scales as `(1−τ²)^{−1}`,
> which is **non-integrable**. `R` is real, but it is only meaningful *combined* with a reference that leaves
> the net `τ` exponent `> −1`. Do not "just add `R`".

---

## 4. What the missing RNA prior currently *is*

The code writes no RNA term. By the Proposition with `logP_r ≡ const`:

```
    "write nothing for RNA"   ⟺   P_r(log ρ_r) = const   ⟺   p_r(ρ_r) ∝ 1/ρ_r    (the Haldane prior)
```

Flat in log-rate **is** Haldane in linear rate — these are the same statement. The codebase describes this as
"`logP_r ≡ 0`, deliberate and declared"; the roadmap describes it as "silently selects Haldane". **Both are
correct**; they are the same fact with opposite emphasis.

The consequence is not cosmetic. Haldane places unbounded prior mass at `ρ_r → 0`, i.e. at `f_g → 1`. With
both arms Haldane, `ψ` is uniform in `λ`, which in composition space is **Beta(0,0)** — improper at *both*
vertices. This is the reference currently in force wherever `logP_g` is uninformative. It is symmetric, so it
hides while nothing breaks the symmetry, and it detonates the moment anything one-sided arrives.

`logP_g` does not rescue the vertices: the fitted prior is projected by linear interpolation with **clamped
ends** (`np.interp(..., left=logP[0], right=logP[-1])`), and the fitted `logP` is measured to be *constant* at
both ends. A constant cannot cancel a `1/f_g` divergence. So `logP_g` supplies real curvature in the interior
and **nothing at either vertex**.

### 4.1 What the AMBIG grid currently implies (and what it does *not*)

It is tempting to call the current AMBIG reference "Dirichlet(0,0,0)". **It is not.** A uniform grid in
`(λ, τ)` has, by the Proposition,

```
    p(f_g, f_p)  ∝  1 / J₂  =  2 · f_g^{−1} (1 − f_g)^{−2}
```

which involves `f_p`, `f_m` only through `(1−f_g)` — not individually, as a Dirichlet would. Marginalizing
`f_p` over `[0, 1−f_g]` gives `f_g^{−1}(1−f_g)^{−1}`, i.e.

```
    current AMBIG reference  =  Beta(0,0) on f_g   ⊗   Uniform on τ
```

*Verified by direct Monte Carlo on the grid measure (§6.1).* This matters for §5: **the code today is already
marginally consistent across the two node classes** — a single-strand node and an AMBIG node imply the *same*
(improper) Beta(0,0) prior on `f_g`. The defect is properness alone, not consistency.

*(Empirically — and we flag this only as motivation, not evidence: the dominant residual error in the
estimator is precisely an unstranded gDNA **over-call**, RNA hallucinated as gDNA, i.e. the `f_g → 1`
direction. That is the signature an unbounded `f_g → 1` vertex would produce. We are not treating this
agreement as confirmation; see §6.4 on why our benchmark cannot currently adjudicate.)*

---

## 5. The reference choice (the part that is *not* forced)

By the Proposition, choosing a reference for an unfitted component means choosing `P_c`. The proposal on the
table is **Jeffreys for a Poisson rate**:

```
    g_c ~ Poisson(ρ_c E_c)   ⇒   I(ρ_c) = E_c/ρ_c ∝ 1/ρ_c   ⇒   p_c(ρ_c) ∝ ρ_c^{−1/2}
    as a density in log-rate:  P_c(log ρ_c) = ρ_c · p_c(ρ_c) ∝ ρ_c^{+1/2}
    ⇒   logP_c  =  +½ · log ρ_c  =  +½ · log f_c  + const
```

So each unfitted component contributes **`+½ · log f_c`**. Substituting into the Proposition:

### 5.1 The two node classes, all components unfitted

```
 n=2:  ψ_ref = ½(log f_g + log f_r) + 0
             = ½ log f_g + ½ log(1−f_g)                      ≡  Beta(½, ½) on (f_g, f_r)

 n=3:  ψ_ref = ½(log f_g + log f_p + log f_m) − log((1−τ²)/4)
             = ½ log f_g + log(1−f_g) − ½ log((1−τ²)/4)      ≡  Dirichlet(½, ½, ½)
```

The `n=3` line is what the roadmap proposed. **Note the sign on the `τ` term**: a naive "add `+½ log f_c` per
component and stop" would give `+½ log((1−τ²)/4)`; the residual `R` flips it. Getting `R` wrong flips whether
the reference rewards or penalizes a strand-dominant tilt.

### 5.2 The classes disagree about `f_g` — and this is a problem

Dirichlet marginals aggregate: `Dir(a₁,a₂,a₃)` has `x₁ ~ Beta(a₁, a₂+a₃)`. Hence

| reference | 1-DOF `f_g` prior | AMBIG `f_g` prior | class-consistent? |
|---|---|---|---|
| **current** (nothing written) | Beta(0,0) | Beta(0,0) | **YES** — but improper |
| **`Dir(½,½,½)`** (the roadmap's proposal) | Beta(½,½) | **Beta(½, 1)** | **NO** |
| `Dir(½,¼,¼)` (§5.3) | Beta(½,½) | Beta(½,½) | **YES** — and proper |

*(Verified by Monte Carlo, §6.1: median `f_g` = 0.2505 under `Dir(½,½,½)` vs 0.5004 under `Dir(½,¼,¼)` vs
0.5000 for Beta(½,½).)*

Under `Dir(½,½,½)`, **an AMBIG node is a priori half as gDNA as its single-strand neighbour** (median `f_g`
0.25 vs 0.50), purely because the annotation says a second strand *could* be transcribed there. The gDNA
question's reference leaks the annotation's strand multiplicity. Two adjacent nodes exchanging messages then
reason from different priors about the same physical quantity.

This is worth emphasizing because it **inverts the stated motivation**. A code comment claims `Dir(½,½,½)`
*"CLOSES the 1-DOF/AMBIG asymmetry: the two classes previously carried DIFFERENT improper priors, which
confounded every dial ever fitted on a mixed benchmark."* Per §4.1 the two classes **currently agree**
(Beta(0,0) both). So `Dir(½,½,½)` does not close an asymmetry — **it opens one that does not presently
exist**, trading a symmetric-but-improper reference for a proper-but-asymmetric one.

### 5.3 A marginal-consistent alternative: `Dirichlet(½, ¼, ¼)`

If we require the reference's **`f_g` marginal to be the same for both classes** — which we argue is the right
requirement, because calibration's *deliverable is the gDNA-vs-RNA split*, the strand tilt `τ` is a **nuisance
parameter**, and per §4.1 **the code already has this property and we should not lose it** — then we need
`a_p + a_m = ½`, giving `a_p = a_m = ¼`:

```
 n=3:  ψ_ref = ½ log f_g + ¼(log f_p + log f_m)·... ⇒  (algebra below)
             = ½ log f_g + ½ log(1−f_g) − ¾ log((1−τ²)/4)
                ╰──────────── identical to the n=2 reference ────────────╯   ╰─ a pure-τ factor ─╯
```

*Derivation:* `Σ_c logP_c = ½log f_g + ¼(log f_p + log f_m)`; add `R = −log((1−τ²)/4)`; substitute
`log f_p + log f_m = 2log(1−f_g) + log((1−τ²)/4)`:
`= ½log f_g + ½log(1−f_g) + ¼log((1−τ²)/4) − log((1−τ²)/4) = ½log f_g + ½log(1−f_g) − ¾log((1−τ²)/4)`. ∎

**The reference factorizes exactly into (the 1-DOF `f_g` reference) × (a pure `τ` term).** The class asymmetry
vanishes identically rather than approximately. Verified numerically to 7e-15 (§6.2) and by Monte Carlo (§6.1).

This is not an ad-hoc patch. It is the shape that **Berger–Bernardo reference-prior theory** prescribes when
there is a **parameter of interest** (`f_g`) and a **nuisance parameter** (`τ`): one constructs the prior
marginally-then-conditionally in the order of interest, rather than taking the joint Jeffreys prior. Jeffreys'
joint rule is well known to misbehave in multiparameter problems — Jeffreys himself recommended against it
there (the standard example being the location-scale model, where the joint rule gives `σ^{-2}` and the
recommended independent construction gives `σ^{-1}`).

**Q2 in §7 asks whether this reasoning is sound and whether `Dir(½,¼,¼)` is the right formalization of it.**

### 5.4 With `logP_g` fitted (the actual production situation)

The rule is *"`logP_c` if fitted, else the reference"*. `logP_g` **is** fitted, so the gDNA arm takes its NPMLE
prior and **must not also take `+½ log f_g`**. Applying the Proposition:

```
 n=2:  ψ = strand + logP_g(log ρ_g)  +  ½ log(1−f_g)
 n=3:  ψ = strand + logP_g(log ρ_g)  +  ½ log(1−f_g)  −  ¾ log((1−τ²)/4)        [under Dir(½,¼,¼)]
       ψ = strand + logP_g(log ρ_g)  +    log(1−f_g)  −  ½ log((1−τ²)/4)        [under Dir(½,½,½)]
```

Neither contains `+½ log f_g`. **We believe the previously-attempted implementation added the full
`Dir(½,½,½)` including `+½ log f_g` while `logP_g` was already supplied — double-counting the gDNA arm.** That
spurious term repels `f_g` from 0, which is exactly backwards on a zero-gDNA test, and we reproduced the
recorded failure value to 3 decimal places (§6.5). We record this as a diagnosis of a previously-open item,
not as an argument for any particular reference.

---

## 6. Numerical verification

Scripts: `ref_derivation_check.py`, `theory_tests.py` (session scratchpad; pure numpy/scipy, no codebase).

### 6.1 The transform is right (Monte Carlo)
Sample from the target composition prior, map to `(λ, τ)`, compare the empirical density to the closed form.

| claim | max rel. error |
|---|---|
| `n=2`, `½log f_g + ½log(1−f_g)` ≡ Beta(½,½) | **0.012** ✓ |
| `n=3`, `½log f_g + log(1−f_g) − ½log((1−τ²)/4)` ≡ Dir(½,½,½) | **0.095** ✓ |
| `n=3`, the naive "+½ per component, no `R`" | **1.018** ✗ |

The naive form is decisively wrong — confirming `R` is real and its sign matters.

Marginal medians of `f_g` (3M samples): `Dir(½,½,½) → 0.2505`; `Dir(½,¼,¼) → 0.5004`; `Beta(½,½) → 0.5000`.
Confirms §5.2/§5.3.

Sampling the **current** grid measure directly (uniform `λ` × uniform `τ`, 4M samples) confirms §4.1: the
`f_g` marginal quantiles (10/50/90) are `[1e-4, 0.502, 0.9999]` for the AMBIG grid and `[1e-4, 0.501, 0.9999]`
for the single-strand grid — **identical**, and `τ` comes out uniform. So the current reference is
`Beta(0,0) ⊗ Uniform(τ)`, class-consistent and improper — not `Dir(0,0,0)`.

### 6.2 The closed forms are exact
Comparing each claimed `(λ,τ)` closed form against the exact transform of its Dirichlet, over a 601×601 grid,
the difference is **constant** to `4.9e-15` (`Dir(½,½,½)`) and `7.1e-15` (`Dir(½,¼,¼)`).

### 6.3 Properness — and a constraint the reference imposes elsewhere

`L`-invariance of the posterior median of the *pure reference* (no likelihood, no messages):

| reference | L=6 | L=10 | L=20 | L=30 | spread |
|---|---|---|---|---|---|
| Haldane both (shipped) | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.0000 |
| Jeffreys `½(log f_g + log(1−f_g))` | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.0000 |
| **flat `logP_g` + `½log(1−f_g)`** | **0.0903** | **0.0133** | **0.0001** | **0.0000** | **0.0903** |

Row 1 is `L`-invariant only by **symmetry** — Beta(0,0) is improper at both vertices, but symmetrically so, so
the median cannot move. Break the symmetry (any one-sided message) and it walks; that is the pathology this
work exists to fix. **Row 1's flatness here is not evidence of properness.**

Row 3 is the important one. It says: **if `logP_g` is flat in some region, adding the RNA reference makes `ψ`
improper at the `f_g → 0` vertex.** And `logP_g` *is* flat below its fitted grid — the projection clamps
(`np.interp(..., left=logP[0])`). So:

> **The reference and the gDNA prior's left tail are coupled. Adding any term that bounds `f_g → 1` converts a
> `f_g → 1` blow-up into an `f_g → 0` blow-up unless `logP_g` is given a genuine (non-clamped) left tail.**
> They must land together. The roadmap independently flagged this clamp as "inert today, load-bearing the
> moment `logP_r` exists"; this makes the mechanism precise and shows it applies to the *reference* too, not
> only to a fitted `logP_r`.

### 6.4 Reparameterization — and why the `λ` grid is the right one

Median of the pure prior computed two ways:

| prior | via `λ`-grid | via uniform `f`-grid | analytic |
|---|---|---|---|
| Beta(½,½) | 0.50000 | 0.50000 | 0.50000 |
| Beta(1,1) | 0.50000 | 0.50000 | 0.50000 |
| **Beta(½,1)** | **0.25004** | **0.21279** | **0.25000** |

The `λ`-grid recovers the analytic median; a uniform `f`-grid does **not**, because it under-resolves the
integrable `f^{−1/2}` vertex singularity. This validates the `λ` parameterization for exactly the
Beta(½,·)-type references under discussion.

### 6.5 The previously-"undiagnosed" failure, reproduced

The roadmap records that implementing the reference left a zero-gDNA regression test at `f_g = 0.564` against
a `< 0.50` bound (previously 0.44), undiagnosed, and the attempt was reverted. Re-running that exact test
chain with the reference injected in-process:

| variant | AMBIG `f_g` | vs. recorded |
|---|---|---|
| shipped (nothing) | 0.4362 | matches "0.44" |
| **`Dir(½,½,½)` incl. `+½log f_g`** | **0.5638** | **matches "0.564"** |
| ~~`logP_g` fitted → RNA-only reference~~ | ~~0.4362~~ | **⚠ MIS-SPECIFIED — see below** |

The recorded number is reproduced. But **the third row was our error, and it retracts our "gDNA
double-count" diagnosis of this particular test.** That test chain runs `node_sweep` with **no `gdna_prior`**
— it is prior-free. Applying the "`logP_g` fitted" variant to it wrote *neither* a fitted `logP_g` *nor* the
gDNA reference, leaving the gDNA arm bare Haldane: an **improper** ψ (exactly the §6.3 row-3 pathology, in the
mirror direction). Its 0.4362 was an artifact of a missing arm, not a passing result.

**Correctly specified, the prior-free reference gives 0.5638 on this chain — and that is right.** Both arms
take their reference ⇒ Beta(½,½) ⇒ `f_g = 0` is forbidden ⇒ the estimate moves *away* from the zero-gDNA
truth. The double-count diagnosis still explains why *`Dir(½,½,½)` with a fitted `logP_g`* is wrong (§5.4
stands), but it does **not** explain this test — this test moved because a **proper reference forbids the
vertex**, precisely as §10.5 predicts. The test's `< 0.50` bound was calibrated against the improper ψ.

Two further reasons the test cannot adjudicate: at its `n_grid=40, L=10` the two grid points straddling `λ=0`
are **exactly 0.4362 and 0.5638** — the recorded "0.44 → 0.564" is *one grid cell*; and the AMBIG node's counts
are balanced under `κ=0.95`, so its strand likelihood genuinely says "gDNA". The truth (`f_g = 0`) is not
identifiable from this node's data at all — only from a prior or a message, and the chain has neither
(no intergenic seeds). **It is a knife-edge on a node with no information.**

### 6.6 What we did **not** learn from the benchmark

We A/B'd the references on the 24-condition oracle benchmark. Means: shipped 0.1317, `Dir(½,½,½)` 0.1317,
§5.4-derived 0.1396. Split by strandedness the derived reference is **best on unstranded** (0.2315 vs 0.2360)
and **worst on stranded** (0.0477 vs 0.0275).

**We are not drawing any conclusion from this.** A correct reference can measure worse while the rest of the
pipeline is mis-specified, and there are known, documented defects upstream and downstream of this term. We
report it only to note that a benchmark improvement is *not* available as a validation criterion for this
change, and to record an observation that may matter later: **the simulated truth lives on the vertices**
(`f_g = 1` exactly for gDNA-only nodes, `f_g = 0` exactly for RNA-only nodes), and every Beta(½,·)-type
reference is precisely a statement that vertices are impossible. Each half of the Jeffreys reference helps
where the truth is at the *opposite* vertex. That is the behavioural signature of an **over/under-call
exchange rate** rather than a neutral regularizer — the same character as an improper `+0.5·λ` "ramp" that
this project previously identified as a disguised gDNA-abundance prior and removed. **Q3 in §7.**

---

## 7. Questions for review

### Q1 — Is "Jeffreys for a Poisson rate" the Jeffreys prior of a model we do not have? *(most important)*

The reference in §5 is derived from `g_c ~ Poisson(ρ_c E_c)` — i.e. from a model in which **each component's
count is observed separately**. It is not. We observe only the **total** mass `M` and a **strand split**
`(u₊, u₋)`. Under the actual observation model, the Fisher information for the quantity of interest is

```
    I(f_g)  ∝  n(½ − κ)²   →   0   as κ → ½
```

i.e. **identically zero for unstranded libraries** — which is exactly the regime where the reference is
load-bearing (it is the regime where `ψ` reduces to the prior alone). The model's own Jeffreys prior for `f_g`
is therefore `|I|^{1/2} ≡ 0` — degenerate, not a prior at all.

This is not a limiting statement. The strand term in the code is **bit-flat** at `κ = ½`: the mixture's mean
is `p = ½f_g + κf_p + (1−κ)f_m ≡ ½` identically on the simplex when `κ=½`, and the term's variance is
(deliberately, per an architectural principle) frozen at a fixed reference composition, so it contributes a
per-node constant that cancels under normalization. Measured spread over the whole `λ` grid: `ptp = 0.0`
exactly in float64. **At `κ = ½` the likelihood is exactly constant in `f_g` and the posterior *is* the
reference, whatever we choose it to be.** Roughly half of the target data (unstranded libraries) is in this
regime.

So the "derivation" appears to import the Jeffreys prior of the **complete-data** model as a prior for the
**observed-data** model. Our concerns:

1. Is that a legitimate, named construction (it resembles an independence/product-of-Jeffreys prior), or a
   category error?
2. This project operates under an explicit architectural principle — *"a fragment count carries **zero**
   intrinsic gDNA/RNA information; the count enters only as the strand likelihood's precision"*. A reference
   derived from **counting each component separately** seems to smuggle in exactly the information the
   architecture forbids. Is that a real contradiction, or is a prior on rates innocuous here because it is not
   conditioning on unobserved counts?
3. If Jeffreys is unavailable because `I(f_g) ≡ 0`, **what is the principled reference for a parameter the
   experiment carries no information about?** Is the honest answer that no data-derived reference exists and
   the choice is frankly subjective (in which case we should declare it as such and choose on other grounds,
   e.g. §5.3's marginal consistency), or is there a construction we are missing?

### Q2 — Is the marginal-consistency argument sound, and is `Dir(½,¼,¼)` its right formalization?

§5.3 argues: `f_g` is the parameter of interest, `τ` is a nuisance parameter, therefore the reference should be
built marginally-then-conditionally (Berger–Bernardo), not as a joint Jeffreys — and that this uniquely
selects `a_p = a_m = ¼` by requiring the `f_g` marginal to match the 1-DOF class. Specifically:

* Is "the two node classes must imply the same marginal prior on `f_g`" a defensible **requirement**, or are we
  entitled to let a node's component count change its `f_g` prior (as `Dir(½,½,½)` does)?
* Is `Dir(½,¼,¼)` actually the Berger–Bernardo reference prior for this problem with `f_g` as the parameter of
  interest, or merely a prior that *happens* to have the right marginal? (We derived it by matching marginals,
  not by running the reference-prior algorithm. We would like the genuine construction if it differs.)
* Does the exact factorization in §5.3 — reference = (1-DOF `f_g` reference) × (pure `τ` term) — have a name?
  It seems too clean to be a coincidence.

### Q3 — Is a proper reference *supposed* to bias a vertex-valued truth?

§6.6: our ground truth is exactly on the simplex vertices. Any Beta(½,·) reference asserts vertices have zero
probability, and so is a guaranteed bias for those nodes. Is that simply the price of properness (and the right
response is "the reference makes `ψ` proper so that *measurements are trustworthy*; the *information* must come
from the fitted `logP_g`/`logP_r`, and one should not expect the reference alone to improve accuracy")? Or does
a truth that genuinely lives on the boundary indicate that a **boundary-admitting** reference (a mixture with
atoms at the vertices, i.e. a spike-and-slab on "this node has no RNA at all" / "no gDNA at all") is the correct
model, and Beta(½,·) is the wrong family?

The physical claim is real, not an artifact of simulation: an intergenic node genuinely has **no RNA** (`f_g = 1`
exactly), and a plasma cfRNA library genuinely has **no nascent RNA**. Structural zeros exist here.

### Q4 — Confirmations sought on the forced part

* §3: is the Proposition correct as stated, in particular `R = −log((1−τ²)/4) − log 2` for `n=3` and `R = 0`
  for `n=2`? (Verified numerically to ~1e-14, but we want the algebra checked.)
* §6.3: is our reading right that a clamped-flat `logP_g` tail plus any `f_g→1`-bounding term is improper at
  `f_g→0`, so the two changes are **not separable**?
* Is there a reason the `n=3` residual `R`, which rewards `τ → ±1` unboundedly (though integrably), is
  problematic for a *numerical* grid integration on `τ ∈ [−1,1]` with closed endpoints? We currently clip at
  half a grid step; is there a standard quadrature treatment for an integrable endpoint singularity here?

---

## 8. Summary table — forced vs. chosen

| statement | status |
|---|---|
| `ψ = strand + Σ_c logP_c + R`, priors given as log-rate densities | **forced** (§3) |
| `R = 0` for single-strand nodes; no Jacobian to write | **forced** (§3) |
| `R = −log((1−τ²)/4)` for AMBIG nodes — currently **missing from the code** | **forced** (§3) |
| writing nothing for RNA ⟺ Haldane ⟺ unbounded `f_g→1` | **forced** (§4) |
| a fitted `logP_c` replaces that component's reference (no double-count) | **forced** (§5.4) |
| the current AMBIG reference is `Beta(0,0) ⊗ Uniform(τ)`, **not** `Dir(0,0,0)` | **forced** (§4.1) |
| the two classes **currently agree** on `f_g` (Beta(0,0) both); `Dir(½,½,½)` would **break** that | **forced** (§4.1, §5.2) |
| at `κ = ½` the likelihood is exactly flat in `f_g`, so the posterior **is** the reference | **forced** (§7 Q1) |
| **the reference should be Jeffreys (`c = ½`) at all** | **CHOSEN** — Q1 |
| **`Dir(½,¼,¼)` over `Dir(½,½,½)` (marginal consistency)** | **CHOSEN** — Q2 |
| **Beta-family at all, vs. a vertex-admitting mixture** | **CHOSEN** — Q3 |

---

## 9. The shortest version *(as submitted for review)*

Three things we believe are **theorems** and want checked: `R = 0` for `n=2` but `R = −log((1−τ²)/4)` for
`n=3` (§3); writing nothing for a component **is** choosing Haldane for it (§4); and `Dir(½,½,½)` gives an
AMBIG node a *different* `f_g` prior (Beta(½,1)) than a single-strand node (Beta(½,½)), whereas `Dir(½,¼,¼)`
gives them the same one and factorizes exactly into (the 1-DOF reference) × (a pure `τ` term) (§5.2–§5.3).

One thing we believe is a **choice masquerading as a derivation** and most want adjudicated: calling
`+½·log f_c` "the Jeffreys reference" (§7 Q1). It is Jeffreys for a model in which each component's count is
observed. Ours observes a total and a strand split, and its Fisher information for the estimand is
*identically zero* on half our data. We would like to know whether that makes the construction unprincipled,
or whether it is a recognized and defensible choice that we should simply declare as subjective and defend on
other grounds (e.g. §5.3's marginal consistency).


---

## 10. REVIEW OUTCOME — the resolved design *(added 2026-07-15, post-review)*

Verification of every claim below: `scripts/debug/reference_prior_bb_check.py` (runs standalone).

### 10.1 What the review confirmed

| our claim | verdict |
|---|---|
| the Proposition; `R = 0` (`n=2`), `R = −log((1−τ²)/4) − log 2` (`n=3`) | ✅ **confirmed** — "algebraically flawless" |
| writing nothing for a component ⟺ Haldane | ✅ confirmed |
| the left-tail coupling (§6.3) | ✅ **confirmed, "deep and correct"** — with a sharper statement of the mechanism: a clamped-flat `logP_g` makes the integrand approach a *constant* as `λ → −∞`, so the integral diverges **linearly in `L`** at the left boundary |
| `Dir(½,½,½)` gives an AMBIG node a Beta(½,1) `f_g` marginal (the class asymmetry) | ✅ confirmed |
| the gDNA double-count diagnosis of the reverted attempt | ✅ (not disputed) |

Independently re-verified here: the multinomial Fisher entries `I_ff = 1/(f_g(1−f_g))`, `I_ττ = (1−f_g)/(1−τ²)`,
**`I_{f_g,τ} = 0` (exact orthogonality, numerically ~1e-11)**; and `√det I = f_g^{−½}(1−τ²)^{−½}`, which is
**exactly `Dir(½,½,½)`** pushed into `(f_g,τ)` coords (max rel. dev. 5.6e-16) — the classic multinomial-Jeffreys
result, reached by a third independent route.

### 10.2 Two corrections to *our* document

1. **`R` must never be shipped alone** (review, §1). Written without a taming RNA-side reference, the `τ`
   density scales as `(1−τ²)^{−1}`, which is **non-integrable** (numerically: the integral diverges, vs 3.16
   for the `−½` exponent). Our §3 framing — "a derived, prior-independent term the AMBIG path is missing" —
   is true but misleading: `R` is **not separable** from the reference choice. Any adopted reference must
   leave the net `τ` exponent `> −1`.
2. **`Dir(½,¼,¼)` is superseded.** It was the best *Dirichlet* approximation to the right answer. Since we
   integrate on a grid, we are not confined to the Dirichlet family (§10.3).

### 10.3 The adopted reference: the exact Berger–Bernardo prior

Because `f_g` and `τ` are **information-orthogonal**, the Berger–Bernardo construction with `f_g` as the
parameter of interest and `τ` as nuisance gives, exactly:

```
    p(τ | f_g)  ∝  I_ττ^{1/2}  ∝  (1 − τ²)^{−1/2}                 (normalized; ∫ = π)
    p(f_g)      ∝  I_ff^{1/2}  ∝  f_g^{−1/2} (1 − f_g)^{−1/2}     = Beta(½, ½)
    ⇒  p(f_g, τ)  ∝  f_g^{−1/2}(1−f_g)^{−1/2}(1−τ²)^{−1/2}
```

and on our grid coordinates:

```
    ψ_ref(BB)  =  ½·log f_g  +  ½·log(1 − f_g)  −  ½·log(1 − τ²)      + const
```

**Properties (all verified):**
* `f_g` marginal is **identically Beta(½,½)** — the *same* as the single-strand class. Marginal consistency
  holds by construction, not by matching.
* `τ` marginal is **identically Beta(½,½)** on `[−1,1]` — the standard Jeffreys prior for a binomial split.
* It factorizes **perfectly**: `Beta(½,½)_{f_g} ⊗ Beta(½,½)_τ`. No leftover fractional exponents.
* It is **outside the Dirichlet family**: pulled back, `p ∝ f_g^{−½}(1−f_g)^{−½} f_p^{−½} f_m^{−½}`, i.e.
  **`BB = Dir(½,½,½) × (1−f_g)^{−½}`**. That extra factor is exactly what repairs the `Beta(½,1)` marginal.
* `τ` exponent `−½ > −1` ⇒ integrable ⇒ §10.2's constraint satisfied.

> **One erratum in the review**, which does not affect its result. §3 of the review writes
> `p(λ,τ) = p(f_g,τ) · J₂` with `J₂ = f_g(1−f_g)²/2`. But `J₂` is the chart for a density in `(f_g, f_p)`;
> `p(f_g,τ)` is a density in `(f_g, τ)`, whose chart to `(λ,τ)` is just `∂f_g/∂λ = f_g(1−f_g)`. Using the
> stated `J₂` would give `½log f_g + (3/2)log(1−f_g) − ½log(1−τ²)`. **The reviewer's final `ψ_ref` is the one
> the correct chart produces** (verified: spread 1.8e-15 against `∂f_g/∂λ`, vs 6.0 against the stated `J₂`) —
> so this is a slip in an intermediate, not in the answer.

### 10.4 The gap BB leaves, and how we close it

BB is a **joint** construction over `(f_g, τ)`. It does **not** factorize into per-component rate priors — the
`(1−f_g)^{−½}` factor is precisely what breaks that. So the rule *"`logP_c` if fitted, else the reference"*
**has no BB slot**: the review derives the reference for the all-unfitted case, but production has `logP_g`
**fitted** (the NPMLE) and eventually wants `logP_r` fitted too.

The resolution — and it is forced, not chosen — comes from noticing that BB's `f_g` marginal, `Beta(½,½)`, is
*identical* to the **two-group** (gDNA vs RNA-**total**) Jeffreys reference, for which `R = 0` and the
per-component rule applies cleanly. So:

> **Treat the composition as a TWO-GROUP split on the `f_g` axis (gDNA vs RNA-total) — which is exactly what
> calibration models — and treat the tilt `τ` as a nuisance carrying BB's own conditional.**

```
    (2-group Jeffreys on f_g)  +  (BB's τ conditional)  ≡  ψ_ref(BB)      [verified EXACT, spread 0.0]
```

This yields the production rule. With `logP_g` fitted and `logP_r` not yet:

```
    1-DOF  :  ψ  =  strand  +  logP_g(log ρ_g)  +  ½·log(1 − f_g)
    AMBIG  :  ψ  =  strand  +  logP_g(log ρ_g)  +  ½·log(1 − f_g)  −  ½·log(1 − τ²)
```

and once `logP_r` is fitted, `+½·log(1−f_g)` is replaced by `logP_r(log ρ_r)`, `ρ_r = (1−f_g)·M/E_r`.

Why this is the right shape:
* **Identical on the `f_g` axis for both node classes** — marginal consistency by construction. The AMBIG path
  differs from the single-strand path *only* by a term in `τ` alone, which is what "τ is a nuisance" means.
* **Reduces to BB exactly** when nothing is fitted (so we inherit the review's endorsement).
* **Leaves exactly two slots** — `logP_g` and `logP_r` — for the two NPMLE priors the roadmap wants, and they
  sit on the **RNA-vs-gDNA axis the architecture says calibration models**, with the strand split factored out
  as the nuisance it is. The tilt never needs a fitted prior.
* **No `+½·log f_g`** anywhere — the gDNA arm takes its fitted prior once (the reverted attempt's bug, §6.5).

### 10.5 Q3 answered: the boundary bias is real, and the long-term fix is spike-and-slab

The review confirms: **yes**, any proper continuous density on `[0,1]` assigns zero mass to the exact
boundaries, so a Beta(½,·) reference *will* systematically bias genuinely-vertex-valued nodes into the
interior. This is the "tax" for interior regularization. Since our structural zeros are **physically real**
(an intergenic node has no RNA; plasma cfRNA has no nascent RNA), the mathematically correct long-term model
is a **spike-and-slab**:

```
    p(f_g)  =  π₀·δ(f_g)  +  π₁·δ(f_g − 1)  +  (1 − π₀ − π₁)·Beta(f_g; a, b)
```

On a grid solver this is cheap: the boundary grid points become **discrete states carrying finite prior mass**;
the interior points carry the continuous density. **Deferred**, but this is the principled destination, and it
predicts recovering exactly the accuracy that §6.6 shows a pure Beta(½,·) reference costs on vertex-valued
truth. *(Note the structural nodes — intergenic, TSS/TES — are already handled outside the solver by a
signature-based lock; the spike-and-slab matters for nodes whose vertex-ness is **inferred**, not annotated.)*

### 10.6 Q1 answered: not a category error — but the license is explicit

The review's answer: using the complete-data Jeffreys prior when only observed data is available is the
well-established **structural / independence Jeffreys prior**, not a category error. The architectural
count-zero-information principle governs how the **likelihood** may process counts; the **prior** may encode
the physical capacity of the generative system. At `κ = ½`, where `I(f_g) ≡ 0` and the posterior *is* the
prior, the structural Jeffreys prior is the principled objective regularizer.

**We accept this, and record the license explicitly:** the reference is a *declared, subjective-but-principled
choice*, justified by (a) the structural-Jeffreys framework, and (b) marginal consistency + orthogonality
(§10.3). It is **not** forced by the observed-data likelihood, which is exactly flat in `f_g` on half our data.
Note also that BB (§10.3) is built from the **same** multinomial (i.e. complete-data) information — so §10.6 is
what licenses §10.3, and the two stand or fall together. This is worth remembering the next time an
accuracy regression is attributed to the reference: **the reference is a choice we can revisit**, and §10.5 is
the named alternative.

### 10.7 IMPLEMENTED — `simplex_logodds.py`, and how it validates

Shipped as `_JEFFREYS_REF = 0.5` with two helpers, `_gdna_arm` / `_rna_arm`, called identically by the 1-D and
2-D paths; `_tilt_grid` now returns **θ ∈ [−π/2, π/2]** with `τ = sin θ`. **No tilt term and no Jacobian are
written anywhere** — both cancel identically. `global_logprior=None` now means "this arm takes its reference",
not "omit this arm": *a prior-free solve is not a reference-free solve.*

**Validation (theory only, per the agreed scope):**

| test | result |
|---|---|
| reference reproduces Beta(½,½) on `f_g` (κ=½, no prior, no messages) | 0.50000 — matches the analytic median |
| **1-DOF vs AMBIG class gap** | **0.00e+00 — closed identically** |
| L-invariance, symmetric ψ | spread **0.0** (both classes) |
| the reference is weak — strand evidence overrides | `f_g` spans 0.0093 … 0.9953 as counts tilt |

**The discriminating test** (roadmap §1.2's node-1479 setup: κ=½ ⇒ strand bit-flat, plus ONE one-sided
message, so ψ = reference + message only). A *symmetric* ψ is L-invariant whether proper or not — symmetry
hides the pathology — so the honest test needs an asymmetric ψ, and the right question is **"does an L→∞ limit
exist?"**, not "is the spread small?":

| L | 6 | 10 | 14 | 20 | 30 | 40 | 60 |
|---|---|---|---|---|---|---|---|
| **shipped** (`_JEFFREYS_REF=0`) | 0.6752 | 0.9286 | 0.9894 | 0.9995 | 1.0000 | 1.0000 | 1.0000 |
| **implemented** (`=0.5`) | 0.5120 | 0.5287 | 0.5315 | 0.5325 | 0.5300 | 0.5300 | 0.5300 |

The shipped ψ marches to the **vertex** (`f_g = 1.0000` — the posterior degenerating to a point mass, which is
what improper *means* here). The implemented ψ converges geometrically to an **interior limit, 0.5300**; the
residual over any finite window is the bracket truncating a *decaying* tail. At production `L = 10` the
truncation error is **0.0013**. Roadmap §1.2 recorded spreads of 0.525 (shipped) / 0.045 (predicted for c=½);
we measure 0.4915 / 0.0475 on the isolated node — its prediction, reproduced.

**Test status (not chased, per scope):** exactly one calibration test fails —
`test_gdna_sweep_zero_gdna_pin_and_monotone` at `f_g = 0.5638` vs its `< 0.50` bound. Per §6.5 this is the
**correct** behaviour of a proper reference on a node with no information, and the bound was calibrated against
the improper ψ. Its chain *is* L-invariant under the new reference. Leave it red until initialization lands,
then re-derive the bound (or retire the chain — its own comment concedes it "cannot occur in a real
annotation"). The 147 already-regenerated golden files will move again.

### 10.8 Open engineering items this creates

1. **The `logP_g` left tail must land with the reference** (§6.3, confirmed). Replace the clamped
   `np.interp(..., left=logP[0])` with a genuinely decaying left tail. **Not separable** from §10.4.
2. **`τ` endpoint quadrature.** The `(1−τ²)^{−½}` singularity is integrable but the closed-endpoint uniform
   grid cannot represent it. The review's recommendation is **Gauss–Jacobi quadrature** (the exact tool for a
   `(1−τ²)^{−a}` weight); our current half-grid-step clip is the engineering fix. Since the AMBIG `(λ,τ)` cube
   is already the sweep's dominant cost, Gauss–Jacobi may pay for itself by allowing a coarser `τ` grid.
3. **`logP_r`'s target is now unambiguous**: the **RNA-total** rate `ρ_r = (1−f_g)M/E_r`, on the 2-group axis —
   not a per-strand prior. The tilt is a nuisance and needs no fit.

---
