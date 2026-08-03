# What a node should accumulate — the derivation

    Verification: `scripts/design/node_density_derivation.py`   (T0–T6, each perturbed)
    Information:  `scripts/design/observable_efficiency.py`     (756 FL pairs x 3 shapes x 4 phi)
    Consequences: `docs/calibration/S5_DESIGN_LOG.md` §1

⭐ **Re-run the scripts. Do not quote this file.**

---

## 1. The observations that started this

The owner's, and each turns out to be exactly right:

1. **An edge is the extreme case.** A 0-bp line has no width, so a fragment either flows across it or does
   not, and `1/L` measures exactly the flux. **There is therefore no reason to accumulate `Σ L` at an
   edge — we already have what we need there.**
2. **Nodes are a spectrum, not a different object.** A 1 bp node is very nearly an edge: essentially
   nothing is contained, everything spans, and it measures flux. As the node grows the picture changes,
   and it changes qualitatively once `node length > fragment length`.
3. **The two limits are known.**
   * at node length 0 the density is `1/fragment_length` — 100 % spanning, 0 % contained;
   * at node length ≫ fragment length the density converges to `count / node_length`.
4. **So there should be a blended equation** — one rule that interpolates between them — and the
   accumulator should deposit *that*, not one fixed weight.

## 2. The problem

The accumulator deposits `round(2³²/L)` at a node and `round(2³²/(L−1))` at a line
(`docs/accumulator/DESIGN.md` §10.1). The design already concedes that the node case is not model-free
(§6: *"at a node the `1/w` does not cancel `(ℓ−w+1)₊`; it is a better-conditioned second moment and
nothing more"*), and treats that as a fact of life.

It is not. `1/(L−1)` is not model-free at an edge *because it is a reciprocal length*. It is model-free
because `L−1` **is the edge's opportunity** — the number of start positions from which a length-`L`
fragment crosses that line. That is a coincidence of the edge frame, and reading it as "reciprocal
length" is what carried the wrong weight over to nodes.

## 3. The questions

| | question | status |
|---|---|---|
| **Q1** | Is there a weight that is model-free at *any* node length, and is it unique? | ✅ **answered — yes, and yes** |
| **Q2** | Does it reproduce the two known limits? | ✅ **answered — exactly, both** |
| **Q3** | Do the contained and spanning populations combine into one statement? | ✅ **answered — exactly, to one missing length** |
| **Q4** | Does making the density model-free *help* the gDNA/RNA split? | ⛔ **answered — NO, and the reason is structural** |
| **Q5** | What should the accumulator therefore store at a node? | ⚠ open — §7 states the trade, the ruling is the owner's |
| **Q6** | Does any of this change the edge? | ✅ **answered — no. Observation 1 is correct** |

---

## 4. The derivation

### 4.1 Setup

A node is the half-open interval `[0, ℓ)`. Component `c ∈ {gDNA, RNA}` has uniform start density `ρ_c`
and length pmf `f_c`. A fragment is `[s, s+w)`. Each population has an **opportunity** `A(w)` — the
number of integer start positions producing that event. Verified by exact enumeration (T0):

```
contained   the whole fragment lies inside          A(w) = (ℓ − w + 1)₊
spanning    it covers the node and both its lines   A(w) = (w − ℓ − 1)₊
crossing    it crosses one 0-bp line                A(w) = (w − 1)₊
```

Deposit some weight `h(w)` per qualifying fragment. Because starts are Poisson,

```
    E[Σ h]  =  Σ_c ρ_c · Σ_w f_c(w) · A(w) · h(w)                              (★)
```

### 4.2 T1 — the reciprocal-opportunity rule, and it is unique

We want (★) to equal `K · Σ_c ρ_c` with `K` **not depending on `f_c`** — that is precisely what
"model-free" means. Then `Σ_w f(w)·A(w)h(w) = K` must hold for *every* pmf `f` supported where `A > 0`.
Taking `f` to be a point mass at each such `w` in turn forces

```
    A(w) · h(w) = K   for every w with A(w) > 0        ⟹        h(w) = 1 / A(w)
```

**The deposit weight is the reciprocal of the opportunity — necessarily, and up to the scale `K`.** The
shipped edge rule `1/(L−1)` is the `A(w) = w−1` instance. It was never about length.

Substituting back into (★):

```
    E[Σ 1/A]  =  Σ_c ρ_c · P_c(A > 0)                                          (T2)
```

Model-free **up to that population's own support truncation** — the probability the fragment could have
qualified at all. Verified over the whole FL grid at every node length:

| weight | violations of "constant `K`" | worst deviation |
|---|---|---|
| **`1/A` (derived)** | **0** | **0.000000** |
| `1/L` (ships) | 266 | 1582× |
| `1` (raw count) | 282 | 19950× |

⭐ The last two rows are the perturbation: a weight *without* the property must show many violations, and
they do. `1/A` shows none, at machine precision.

### 4.3 T2/T3 — the two populations are complementary

`A > 0` means `w ≤ ℓ` for contained and `w ≥ ℓ+2` for spanning. Those cover every length except one, so
depositing the reciprocal opportunity on **both** and adding:

```
    E[ Σ 1/(ℓ−L+1)  +  Σ 1/(L−ℓ−1) ]  =  Σ_c ρ_c · ( 1 − f_c(ℓ+1) )            (T3)
```

⭐ **This is observation 4's blended equation.** One rule, every node size, and the only defect is the pmf
mass at the *single* fragment length `w = ℓ+1` — a fragment that exactly overhangs the node by one base
qualifies for neither population. Verified to **4.2e-16**, i.e. exactly, and confirmed end-to-end through
a simulated Poisson fragment process (all |z| < 1.3).

⚠ The `ℓ+1` fragment is not lost from the system — it still crosses one of the node's two lines and is
counted there. It is missing only from this node's own pair. For a smooth pmf the deficit is ~0.4–0.8 %,
it is component-specific (so it does not cancel in the composition), and it is a **known function of the
fitted FL model**, hence correctable post hoc if it ever matters.

### 4.4 T4 — both limits fall out

**`ℓ = 0`.** Contained opportunity `(0−w+1)₊ = 0` for every real fragment; spanning opportunity becomes
`w−1`, so the deposit becomes `1/(L−1)`. ⭐ **The node rule at `ℓ = 0` IS the edge rule, identically** —
verified to `1e-12` on every FL in the grid, with the contained channel exactly 0. Observation 2 is a
theorem: an edge is the `ℓ → 0` limit of a node, not a separate kind of object.

**`ℓ ≫ E[L]`.** Spanning empties; contained gives `1/(ℓ−w+1) → 1/ℓ`, so `Σ h → count/ℓ`. Observation 3's
second limit, exactly.

### 4.5 T6 — what this buys: the effective length, done exactly, with no model

`count/ℓ` is biased by `(ℓ − μ + 1)/ℓ`, which is why `effective_length.region_eff_length` exists and why
it needs a **fitted** FL model. The reciprocal deposit applies the correction per fragment, using that
fragment's own length, at deposit time:

| node ℓ | FL mean | `count/ℓ` (× truth) | `Σ 1/A` (× truth) |
|---|---|---|---|
| 151 | 200 | 0.0551 | 0.2642 *(× the containment probability — see below)* |
| 400 | 200 | 0.5046 | 0.9908 |
| 1000 | 200 | 0.8015 | **1.0000** |
| 3000 | 200 | 0.9338 | **1.0000** |
| 3000 | 50 | 0.9838 | **1.0000** |

⭐ At a 1000 bp node with 200 bp fragments the naive density is **20 % low** and the derived one is exact
to six figures. The 151 bp row is not a failure: `Σ1/A` on the contained population alone estimates
`ρ·P(w ≤ ℓ)`, and only 26 % of 200 bp fragments fit in a 151 bp node. Adding the spanning channel
recovers the missing 74 % — that is T3.

---

## 5. ⛔ Q4: the model-free density does NOT help the split, and that is a theorem too

Model-free means the coefficient `K` in (★) is the same for gDNA and RNA. Two components with **identical
coefficients** contribute a row `(K, K)` to the estimating system. Its determinant against any other row
is proportional to `K·(a_g − a_r)`, so the model-free channel contributes **exactly zero** discrimination
by itself:

> **A channel cannot simultaneously be model-free and informative about composition.**
> Being model-free *is* the statement that it cannot tell the components apart.

Measured, as `Σ1/A` alone over both node populations: efficiency **0.113 median at a 151 bp node, 0.000 at
1000 bp** — it recovers the level and almost nothing else. Exactly as the theorem requires.

⚠ **This is the same structural fact as the edge blind spot in `docs/calibration/S5_DESIGN_LOG.md` §1 A1**, seen from the
other side. There the pair `(count, Σ1/(L−1))` has determinant `μ_g − μ_r` and dies at equal means; here
the flat channel has determinant 0 outright. **Discrimination requires a deliberately length-TILTED
channel.** The level and the split are two jobs and no single number does both.

## 6. The measured trade

Efficiency = `Var(full length histogram) / Var(stored set)` — the fraction of the available information
kept. Median / **minimum** over 756 FL pairs, gamma (the ranking is unchanged under lognormal and normal,
and at `φ_g ∈ {0.05, 0.2, 0.5, 0.8}`):

| stored per node population | 25 bp | 151 bp | 1000 bp | 3000 bp | edge |
|---|---|---|---|---|---|
| `count, Σ1/L` — **ships** | 0.758 / 0.002 | 0.832 / 0.078 | 0.197 / 0.001 | 0.182 / 0.000 | 0.324 / 0.000 |
| `count, Σ1/A` — derived | 0.347 / 0.000 | 0.684 / 0.000 | 0.347 / 0.000 | 0.563 / 0.000 | *(≡ ships)* |
| `count, Σ1/L, ΣL` | 0.836 / 0.376 | 0.953 / 0.188 | 0.722 / 0.064 | 0.692 / 0.047 | 0.686 / 0.078 |
| `count, Σ1/A, ΣL` | 0.664 / 0.123 | 0.918 / 0.084 | 0.749 / 0.012 | **0.804** / 0.020 | *(≡ above)* |
| **`count, Σ1/A, Σ1/L, ΣL`** | **0.867 / 0.409** | **0.960 / 0.230** | **0.782 / 0.074** | **0.832 / 0.064** | **0.686 / 0.078** |

⭐ **`Σ1/A` does not replace `Σ1/L`; it is a different quantity that answers a different question.**
Swapping one for the other *loses* composition information at short nodes (0.832 → 0.684 at 151 bp).
Adding it gains, most at long nodes (0.692 → 0.832 at 3000 bp).

### 6.1 What `Σ1/A` costs in variance, and where

Per-unit-`ρ` effective sample size `n_eff = mean²/variance` of each channel, contained population:

| node ℓ | FL | `count` | `Σ1/L` | `Σ1/A` |
|---|---|---|---|---|
| 151 | 200 | 0.416 | 0.383 | **0.148** ← noisy |
| 400 | 100 | 15.08 | 12.96 | **14.84** |
| 3000 | 200 | 140.1 | 120.4 | **140.0** ← free, and exact |
| 3000 | 100 | 145.1 | 124.5 | **145.1** |

And on the spanning side, FL 200: at `ℓ = 1`, `n_eff` 9.88 (count) vs **8.63** (`1/A`); at `ℓ = 151`,
2.79 vs **0.87**.

⭐ **The mechanism, and it is tidy.** `1/A` blows up when typical fragments sit near their population's
boundary (`A ≈ 1`) — contained fragments that nearly fill the node, spanning fragments that barely cover
it. That happens only when **`ℓ ≈ fragment length`**. Away from the crossover in *either* direction `1/A`
is as efficient as the raw count *and* exactly unbiased:

* `ℓ ≪ FL` — spanning dominates, `A = w−ℓ−1 ≈ w`, so `1/A ≈ 1/w`. **This is why the edge rule works.**
* `ℓ ≫ FL` — contained dominates, `A = ℓ−w+1 ≈ ℓ`, so `1/A ≈ 1/ℓ`. Strictly better than `1/L` here on
  *both* counts (`n_eff` 140 vs 120, and unbiased vs 20 % low).
* `ℓ ≈ FL` — the crossover. `1/A` is noisy; the sum over both populations is still exactly `ρ(1−f(ℓ+1))`.

---

## 7. Q5 — where this leaves the accumulator

**Settled by the derivation:**

* **The edge is untouched.** `1/(L−1)` is already the reciprocal opportunity, and observation 1 stands:
  there is no reason to accumulate `Σ L` at an edge *for the level*. (`Σ L` at an edge is still on the
  table for the **split** — that is `docs/calibration/S5_DESIGN_LOG.md` §1 A1's equal-means blind spot, a different
  question.)
* **`weight = L` at a node is simply the wrong argument.** It is neither the opportunity nor a
  deliberately chosen tilt. `docs/accumulator/DESIGN.md` §10.1's *"weight is `L` at a node and `L−1` at an
  edge"* should read *"weight is the population's own opportunity"*, which then covers both.
* **`docs/accumulator/DESIGN.md` §6's `density = E[placements/weight]` becomes identically 1** at every object
  when `weight = placements`. The design wrote the general form and then instantiated it inconsistently.

**The open ruling.** Storing `Σ1/A` buys an exactly model-free density at every object — which removes
the fitted FL model from the level path entirely, makes `docs/SESSION_HANDOFF.md` §2's *"density is the
frame-invariant currency"* true at nodes as well as edges, and reduces `effective_length.py` (old R4) to
the count channel alone. It does **not** buy composition information. The composition needs a tilted
channel, and the measurement says the tilt worth adding is `Σ L`.

### 7.1 What each channel actually buys — the answer to "why four?"

⭐ **Four is not needed. Three is, and the case for the fourth is architectural rather than statistical.**

| channel | what it buys | can it be dropped? |
|---|---|---|
| `count` | the statistical power. A Beta-Binomial needs an integer, and `Var(log ρ_c) = 1/(f_c·n)` (`docs/SESSION_HANDOFF.md` §2) | **no** |
| one tilted channel (`Σ1/L` *or* `ΣL`) | the gDNA/RNA split at all. `count` alone scores 0.596 median / 0.000 min | **no** |
| a *second* tilted channel | the **equal-means blind spot** — the one hard failure. Min at 151 bp `0.078 → 0.188`; at an edge `0.000 → 0.078` | only if that failure is accepted |
| `Σ1/A` | the **model-free level**. Median 151 bp `0.953 → 0.960` — statistically almost nothing | **yes**, on information grounds |

⚠ **And the architectural claim for `Σ1/A` is weaker than it first appears.** It makes the *level*
model-free while the *split* still needs the fitted FL models, so it does not remove the length model
from calibration — only from one of the two paths. Its one substantial information gain is at long nodes
(0.692 → 0.832 at 3000 bp), which is where the pure gDNA background is anchored.

### 7.2 The three candidate rulings, priced

| | store per population | model-free level? | eff. 151 / 3000 bp / edge | per worker | ×8 workers |
|---|---|---|---|---|---|
| **R-a** | `count, Σ1/L` | no | 0.832 / 0.182 / 0.324 | 109 MB | 0.85 GB |
| ⭐ **R-b** | `count, Σ1/L, ΣL` | no | **0.953 / 0.692 / 0.686** | 179 MB | 1.40 GB |
| **R-c** | `count, Σ1/A, Σ1/L, ΣL` | yes | **0.960 / 0.832 / 0.686** | 249 MB | 1.94 GB |

Against a measured 8.6 GB peak RSS, so R-b costs +6 % of peak and R-c +13 %.

⭐ **R-b introduces no new numerical machinery at all: `Σ L` is a sum of integers.** No fixed point, no
scale, no rounding rule, no overflow scheme. Every objection to the density encoding (§7.3) simply does
not apply to it. R-c's `Σ1/A` is the only channel that needs the fixed-point apparatus extended.

⚠ **R-c's `1/A` needs no new *code*, though.** `density_quantum(placements)` already takes a placement
count; the deposit passes `L` where the derivation says it should pass `(ℓ−L+1)` or `(L−ℓ−1)`, and the
node length is in the cut array the deposit already binary-searches.

### 7.3 ✅ The `A = 1` overflow — priced, and it is a non-issue

Flagged in an earlier draft of this file as needing re-pricing before R-c could be ruled. Done:

| fragments landing exactly on `A = 1` | `Σ` at scale 2³² | share of uint64 |
|---|---|---|
| 10⁶ | 4.3e15 | 0.02 % |
| 10⁷ | 4.3e16 | 0.23 % |
| **10⁸** | 4.3e17 | **2.3 %** |
| 10⁹ | 4.3e18 | 23 % |

`A = 1` requires a fragment to *exactly* fill a node (or exactly overhang it by two), so it is a rare
event, and even assuming every one of 10⁸ fragments did it the accumulator is at 2 % of capacity.
⛔ **A second, smaller scale for the `1/A` channel is therefore unnecessary**, and would be actively
harmful: two scales in one schema is precisely the "do not give two quantities one column name" failure
`docs/accumulator/DESIGN.md` §10.1 warns about.

### 7.4 ⚠ Why integers at all — the stated reason is against a strawman

`docs/accumulator/DESIGN.md` §10.1 justifies uint64 fixed point by comparing it to **float32** (36.6 MB
against 73.3 MB at human scale) and `docs/SESSION_HANDOFF.md` §1 fact 11's measured nondeterminism — 17/28 and
20/28 cells differing, max relative **3.7e-7** — is a **float32** number.

⛔ **float64 is the same 8 bytes as uint64, and was never compared.** By the same mechanism its
worst-case merge nondeterminism would be ~`n·eps` = 1e-8 at 10⁸ fragments and ~1e-12 in practice, i.e.
some 10⁵× smaller than the figure the decision was made on. **The precision argument as written does not
support the conclusion.**

The argument that does survive is a testing one, and it should be stated as such: **byte-identity is the
S3 gate.** The C++ is accepted only if it reproduces `tests/native/_accumulator_reference.py` exactly.
With floats that gate becomes "agrees within a tolerance", every downstream check inherits a tolerance,
and this project's most expensive bug (`docs/SESSION_HANDOFF.md` §3 trap 2 — an exact factor of 2 hidden for
months behind a `min()` clip and cancelling fixtures) is exactly the class a tolerance conceals. Integer
accumulation makes a single differing bit a real defect rather than noise.

⚠ That is a reason to keep integers, not a reason to believe float64 would be numerically inadequate.
The two should not be conflated again.

### 7.5 A coarse length histogram instead of moments — tested, not recommended

The natural challenge to storing four moments is that the same bytes buy a **length histogram**, which
is the sufficient statistic and needs no deposit-weight decision at all (the exact `1/A` can be applied
post hoc, per node, with that node's own `ℓ`). Measured at equal storage, fixed global geometric bins:

| | 25 bp | 151 bp | 1000 bp | 3000 bp | edge |
|---|---|---|---|---|---|
| 3 bins (12 B) vs **R-a** (12 B) | 0.552 / **0.758** | 0.710 / **0.832** | **0.630** / 0.197 | **0.629** / 0.182 | **0.564** / 0.324 |
| 5 bins (20 B) vs **R-b** (20 B) | 0.748 / **0.836** | 0.851 / **0.953** | **0.803** / 0.722 | **0.799** / 0.692 | **0.756** / 0.686 |
| 7 bins (28 B) vs **R-c** (28 B) | 0.862 / **0.867** | 0.914 / **0.960** | **0.892** / 0.782 | **0.892** / 0.832 | **0.868** / 0.686 |

⭐ It wins clearly at long nodes and at edges and loses at short ones, because the bins are **fixed
globally** while the opportunity `A` depends on each node's own `ℓ` — a 25 bp node needs resolution
near 25 bp that global bins do not have. Uniform bins are worse than geometric everywhere.

It is a real alternative with two genuine advantages — pure integer counts (no fixed point at all) and
no weight fixed at scan time — but at equal bytes it does not dominate, and it is a far larger change.
**Logged, not recommended.** Revisit if per-node relative bin edges (a fraction of `ℓ`) are ever worth
exploring, since that is what would fix its short-node weakness.

## 8. ⛔ Still open

1. **Fixed-point headroom for `1/A`.** `A` can be **1** (a fragment exactly filling a node), so the
   quantum is `2³²` — the largest value the scheme can emit. At 10⁸ such fragments the sum is 4.3e17,
   inside uint64's 1.8e19 but with far less headroom than the `L ∈ [20,2000]` case
   `docs/SESSION_HANDOFF.md` §1 fact 23 priced. **Must be re-priced before R-c is ruled.**
2. **Junction edges are not in the grid yet** as their own frame with their own reach geometry.
3. **The `ℓ+1` hole** is quantified (~0.4–0.8 %) but not decided: correct it from the fitted pmf, or
   document it as a known deficit.
4. ⚠ **None of this is measured on real data.** Every number here is analytic or simulated. The five
   pure FL pools are the first model-free length measurement available and have never been read off a
   real cfRNA scan.
