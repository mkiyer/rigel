<!-- title: Message Precision — Why It Overshoots, and the Count-Currency Fix -->
# Calibration message precision: why it overshoots, and how to damp it

**Status:** design note for review. **Thesis:** the imputation precision is computed in *absolute density*
space, where a near-empty source looks *certain*; the fix is to bound it by the source's **fragment
count** — the fair currency — which we can do as a smooth, derived variance floor.

---

## 1. How precision is computed today

Each per-node solve maximizes `log ψ = L_strand(tilt) + Σ −½·τ_f·(f − μ_f)²` over the simplex. The
strand term's precision is the count-based Beta-Binomial Fisher information `≈ N·(2κ−1)²`. Each
imputation **message** (`bp_solver._message`, source → destination) contributes a `(μ_f, τ_f)`:

```
μ_f  = (ρ_src·E_dst − spliced_dst) / M_dst                     # imputed fraction (identity density)
τ_ρ  = 1 / ( σ²_bio(μ) + ρ_src/E_src )                          # precision in DENSITY space
τ_f  = τ_ρ · (M_dst / E_dst)²                                   # Jacobian density → fraction
```

with `ρ = C/E` (density = count / length). The two pieces of `τ_ρ`:

- **`ρ_src/E_src`** is the Poisson sampling variance of the source density: `Var(ρ) = Var(C)/E² =
  C/E² = ρ/E`. So its inverse — the sampling *precision* — is `E²/C`, the **Fisher information of a
  Poisson rate**.
- **`σ²_bio(μ)`** is the learned biological dispersion (the monotone Poisson-offset `var~mean` — the
  excess of the cross-node residual over the sampling floor). **This is where overdispersion lives, and
  it is already learned empirically** from the data.

In the regime that bites (adjacent densities similar, `σ²_bio → 0`): `τ_ρ ≈ E_src²/C_src` and therefore
`τ_f ≈ (E_src²/C_src)·(M_dst/E_dst)²`.

## 2. Why it overshoots

Three compounding effects, all rooted in working in **absolute density** space:

**(a) The Poisson Fisher information `E²/C` diverges at low count.** A density estimated from a finite
count has precision `E²/C`. As `C → 0` this → ∞: a source with `0.1` fragments over a 1 kb region has
`E²/C = 10⁶/0.1 = 10⁷` "precision" — because a near-zero density has near-zero *absolute* variance. But
(your point) 0.1 fragments carries essentially **no information**. The absolute-variance view says "I'm
certain the density is ~0"; the truth is "I have barely sampled." The `c₀=1` pseudocount bounds this at
`E²/c₀`, but for a kb-scale region that is still `≈10⁶`.

**(b) The `(M/E)²` Jacobian compounds it.** For a capture exon `M_dst ≈ 16000`, `E_dst ≈ 1000` ⇒
`(M/E)² = 256`. So `τ_f ≈ (E²/C)·256`, and `10⁷ · 256 ≈ 10⁹` — exactly the spike we measure.

**(c) Density blends count and length, so the model cannot tell them apart.** Your example:

| count | length | density | fragments of information |
|---|---|---|---|
| 1 | 10 bp | 0.1 | **1** |
| 10 | 100 bp | 0.1 | **10** |
| 100 | 1000 bp | 0.1 | **100** |

Same density, precision differing 100×. A `var~mean` indexed by *density* sees all three as identical.
Precision is a function of **(count, length, overdispersion)** — never of density alone.

**Net effect.** Substituting `M_dst ≈ C_dst` and `E_src ≈ E_dst`, the message precision is

```
τ_f  ≈  C_dst² / C_src
```

— **quadratic in the destination's count, inverse in the source's count.** This is the opposite of a
fair currency: the node with more mass dominates, and a near-empty source can carry enormous precision
about a node it knows nothing about. (This is precisely the objection we raised against naive count
messages — and it turns out the *current* density formula has it too, hidden in the Jacobian.)

## 3. First principles: the count is the currency

The strand likelihood already speaks in counts (`N·(2κ−1)²`). For a message to compete *fairly* it must
speak the same currency: **a message built from `C_src` fragments can carry at most ~`C_src`
observations of precision about the destination's composition, and should compete *linearly* with the
destination's strand reads (`C_src` vs `N_dst`), not quadratically.**

The information ceiling, derived: a message imputing the gDNA fraction is backed by the source's
`N_src := ρ_src·E_src` gDNA fragments. The Fisher information of a *fraction* estimated from `N`
fragments is `≈ N/(f(1−f)) ≈ N` (binomial). So the message's fraction-precision cannot exceed `N_src`:

```
τ_f  ≤  N_src  =  ρ_src · E_src        (the source's component FRAGMENT COUNT — data, not a knob)
```

**Why low counts must be imprecise, done right (your pseudocount point).** In absolute density space the
Poisson floor `ρ/E → 0` at low count, so it *looks* precise. In **count / log / relative** space the
Poisson noise is the coefficient of variation `≈ 1/√C`, which → ∞ at low count — *correctly* imprecise.
This is exactly why such models are normally fit in **log space**: the variance is relative, so few
counts ⇒ high relative variance ⇒ low precision, automatically. Deviating to absolute space (to bolt on
the Poisson offset) is what reintroduced the "low count looks certain" pathology.

## 4. Proposals

### A. A count-derived variance floor (short-term — recommended; smooth, no cliff)

The cleanest statement of the ceiling is a **variance floor**: a message based on `N_src` fragments
cannot have less than `1/N_src` variance on the fraction it imputes. Add it to the message variance:

```
Var_f_eff = 1/τ_f + 1/N_src        ⇔        τ_f_eff = τ_f · N_src / (τ_f + N_src),   N_src = ρ_src·E_src
```

- **Smooth.** This is the harmonic combination (precisions in series / a variance floor) — `τ_f_eff → τ_f`
  when `τ_f ≪ N_src`, `→ N_src` when `τ_f ≫ N_src`. No `min()` cliff; the operator stays continuous.
- **First-principles.** The floor `1/N_src` *is* the count-limited minimum variance: you cannot know a
  fraction to better than `1/C` from `C` fragments. No magic constant — `N_src` is the measured count.
- **Fixes the overshoot exactly.** The `10⁹` message from a `0.1`-fragment source → `Var += 1/0.1 = 10`
  ⇒ `τ_f_eff ≈ 0.1` (precisely your intent: 0.1 fragments ⇒ ~0.1 observations). A message from a 1000-
  fragment source → `≈ 1000`, competing one-for-one with the strand. A well-supported message
  (`τ_f ≪ N_src`) is unchanged.
- **Drop-in.** One line in `_message`; the rest of the machine is untouched.

Validate on the dissection (does the `10⁹` spike vanish; do reg231 / reg242 hold) **and** the net-flow
benchmark (does it shrink the capture-off / unstranded overshoot the global mean + uncapped messages
produced).

### B. Fit the reliability in count / log space (medium-term)

Re-parameterize the `var~mean` from `σ²_bio(density)` to the residual reliability versus the **source
count** (or in log-density). Then "low count ⇒ low precision" is *learned* rather than floored, and
**overdispersion is the high-count plateau** — the biological excess that survives as `C → ∞`. This puts
the Poisson/biological decomposition in the space where it behaves (your log-space observation), and
removes the absolute-density "low count looks precise" pathology at the *source* as well as the message.

### C. Native count-space messages (long-term — Alt A's narrow win)

Express each message as `~N_src` pseudo-observations at the imputed composition, injected into the count
likelihood instead of a quadratic fraction penalty. Then the precision is bounded by `N_src` **by
construction**, there is no Jacobian, and it is the fairest possible currency. This is the end state the
adversarial review endorsed for the global + spliced terms; (A) is its smooth approximation on the
existing machine, and (C) is where it converges.

## 5. Recommendation

Ship **A** now: a count-derived variance floor `1/N_src` is the smooth, first-principles damping that
bounds every message at its source's information content — no cliff, no magic, and it makes the count
the currency exactly as the strand term already does. It directly removes the `10⁹` spikes and should
shrink the capture-off / unstranded overshoot, after which precision is genuinely stable and we commit.
Then pursue **B** (the `var~mean` in count/log space — the right home for overdispersion) converging to
**C** (native count-space messages). Each step is validated on the dissection + net-flow benchmark
before commit; the `priors.py` factor-1-under-uniform invariant gets a test before any golden regen.

### Open question for the long-term

Is the destination-side `(M/E)²` Jacobian still needed once the floor is in? Under (A) it remains (we
cap, not remove). Under (C) it disappears (counts never leave count space). The interim question is
whether `N_src` alone is the right ceiling or whether the *destination's* count should also enter the
floor (a message onto a 2-read node should also be gentle) — i.e. `1/N_src + 1/N_dst`. Worth testing in
the (A) prototype.
