# The DNA Prior — Sampling-Likelihood Projection & Gentle Anchor (plan)

Authoritative plan for Role B (the gDNA composition hyperprior). This revision fixes the **projection**, which
is the whole point: the prior is a **third line of defense** that gently anchors the solver, not a replacement
for it.

---

## 0. What the DNA prior is (and is not)

**Goal:** measure *how much DNA is in this sample* — and, under hybrid capture (or cancer CNV), how much
**conditional on enrichment**. Off-capture this is one mode, read straight off. On-capture it splits into
"how much DNA when off-target" vs "when on-target" — and enrichment is a **spectrum**, not a binary (a node may
partially overlap a probe), so the landscape can be any shape: multiple modes, peaks, valleys.

**Three lines of defense** (the prior is the third):
1. **strand** — absent at ambiguous / unstranded nodes.
2. **message propagation** — present, imperfect (it assumes uniform coverage, which capture breaks), but better
   than nothing. **This is what solves the node.**
3. **the DNA prior** — a **gentle anchor** on top of line 2. It does *not* set `f_g`; it tells the solver which
   DNA level the node was most likely *drawn from*, and lets the solver gravitate there.

**The invariants that break the circularity** (they do not change across iterations):
- **intergenic** reads = pure DNA (the off-target/background mode).
- **intronic** nodes = DNA relative to intergenic (+ nascent RNA).
- **spliced** fragments = direct evidence of RNA.
These anchor the landscape so the iterate converges instead of drifting into its own error.

---

## 1. The algorithm

### 1.1 The landscape `P(a)`
An **additive Poisson-lognormal** density (D1) over DNA log-density `a = log ρ_g`, fit **per library** on the
**confidently-solved** nodes (§1.4). It naturally carries the depleted (off-target) mode, the enriched
(on-target) mode(s), and the enrichment spectrum between them. Each solved node deposits one unit-mass kernel;
no EM competition (so a minority enriched mode survives); kernel width `h` (bandwidth).

### 1.2 The projection — "what DNA level was this drawn from?"
This is the crux, and it is **not** "evaluate the landscape along the ray and take the tallest reachable mode."
Given the node's **observed** DNA log-density `d = log(f_g^cur · M/E)` (from the current solve), ask: *across
the whole landscape, which DNA level most plausibly generated this observation?* That is a **sampling
likelihood** — a Gaussian in log-density — integrated against the landscape:

```
r_j   ∝  w_j · N(d ; μ_j , h)            responsibilities  (Σ_j r_j = 1)
μ*(d) =  Σ_j r_j · μ_j                    the projected DNA level  (the population source)
v*(d) =  Σ_j r_j · (μ_j − μ*)²  +  h²     the projection variance
```

Read the properties directly off `r_j`:
- **Distance dominates, not height.** A mode 17 logs away contributes `e^{−(17/h)²/2} ≈ 0` responsibility no
  matter how tall — you cannot draw an observation of 18 from a mode at 1. A **nearby hill wins** a far
  mountain.
- **Valleys are honest.** An observation between two modes splits its responsibility ⇒ `v*` is **large** ⇒ a
  weak anchor (§1.3). A node sitting in a mode gets `v* ≈ h²` ⇒ a firm anchor.
- **The spectrum is native.** Partial enrichment (a node between depleted and enriched) yields an intermediate
  `μ*` with wide `v*` — exactly "partially on-target, and I'm not sure."

This is **exactly `DensityNPMLE.project`** (already implemented for Role A) — evaluated on the **Role-B**
landscape at the node's **DNA** density (pass `mass' = f_g^cur · M`, `eff' = E`, so `mass'/eff' = ρ_g`).

### 1.3 The anchor — a gentle pull, not a verdict
`μ*` is **not** the answer; the node's own evidence may already be the truth. The prior enters as a **gentle
Gaussian anchor** toward `μ*`, added to the always-present symmetric reference:

```
ψ_gDNA(λ)  =  reference(½·log f_g)   −   ½ · p · ( log ρ_g(λ) − μ* )²
                                         p = w_anchor / v*
```

- The **reference is retained** (symmetric Beta(½,½) with the RNA arm) — it keeps the `f_g→0/→1` bounds, so the
  prior can never *delete* the barrier (the crush cause). The anchor rides on top.
- `p = w_anchor / v*` makes the pull **firm in a mode, weak in a valley** automatically; `w_anchor` is the
  overall third-line gentleness (small — the messages must still win). The solver gravitates the DNA level
  toward the nearest mode, which anchors the RNA level too.

### 1.4 Iteration & non-circular training
`d` uses the current solve, so this is iterative: `pass-0 → project → anchor → re-solve → re-fit → …`. It
converges because the **invariants** (§0) don't move: fit the landscape on the **confidently-solved** nodes
(low pass-0 variance and/or strand-/splice-/structurally-resolved — intergenic, spliced exons, resolved
introns), **excluding** the ambiguous, high-variance nodes whose `f_g` is still in play. So the landscape
encodes the *true* per-library DNA amount, not the over-call it is trying to fix.

---

## 2. Why this is the right task (contrast with what shipped)

| | **δ-pin** (current / what I built) | **sampling-likelihood anchor** (this plan) |
|---|---|---|
| question | "over all `f_g`, where is the landscape tallest along the ray?" | "which DNA level was *this observation* drawn from?" |
| far tall mode | **wins** (evaluate it at its own peak by choosing `f_g`) | **ignored** (≈0 responsibility from `d`) |
| nearby hill | loses to any taller mode | **wins** |
| valley node | snaps to the global mode | wide `v*` ⇒ weak anchor (honest) |
| role | tries to *set* `f_g` (a solver) | *anchors* the solver (third line) |
| crush | yes (depleted mode is tallest) | no (anchors to the node's *own* nearby mode) |

The δ-pin, the Stage-0 floor, and the enrichment-conditioning log-shift were all attempts to patch the wrong
projection. This plan **removes all three** and replaces them with the projection above.

---

## 3. Implementation (exactly)

**Keep:** additive Poisson-lognormal Role-B fit (D1). **Role A untouched.** **Remove:** the δ-pin
(`gdna_prior.logprior` along the ray), the Stage-0 `logaddexp` floor, the enrichment-conditioning eff-scaling.

### 3.1 The projection + anchor (bp_solver.node_sweep, replacing the `global_lp` block)
```python
# THE DNA-PRIOR ANCHOR (gdna_hyperprior_plan.md §1.2–1.3): project the node's CURRENT DNA density onto the
# Role-B landscape → the sampling-likelihood (μ*, v*), then a GENTLE Gaussian pull toward μ* on the f_g grid.
if gdna_prior is not None:
    d_obs_mass = np.maximum(belief.f_g, 0.0) * mass_global            # ρ_g = f_g^cur · M  (÷E inside project)
    mu_star, v_star = gdna_prior.project(d_obs_mass, eff_global)      # REUSE project on the Role-B landscape
    log_me = np.log(mass_global) - np.log(eff_global)                 # (m,) = log(M/E)
    log_rho = np.log(solve_grid)[None, :] + log_me[:, None]           # (m,K) = log ρ_g at each grid point
    p = w_anchor / np.maximum(v_star, _EPS)                           # firm in a mode, weak in a valley
    global_lp = -0.5 * p[:, None] * (log_rho - mu_star[:, None]) ** 2 # the anchor, (m,K)
else:
    global_lp = None
```
`_gdna_arm` then returns **reference + global_lp** (additive — the anchor rides on the retained symmetric
reference), instead of today's "replace the reference." This is one edit to `_gdna_arm` and covers **both**
solvers (1-DOF `_solve_nodes_logodds` and 2-DOF `_solve_ambig_logodds`, which both consume `global_lp`).

### 3.2 The `w_anchor` gentleness
`w_anchor` is the third-line strength (a pseudo-observation weight). Derive it from `n_eff` (≈0.15,
`gdna_background_floor_derivation.md`) rather than tuning; **flag for owner sign-off** before it lands. It is
the *only* new constant, and `v*` already carries the per-node adaptivity, so it is a single global scalar.

### 3.3 Non-circular training subset (calibrate._fit_gdna_hyperprior)
Fit on `sel = live & isr & confident`, where `confident` = low pass-0 `var_gdna` (self-calibrating) ∪ the
structurally-resolved invariants (intergenic `gonly`, intron, clean spliced exons). Exclude the high-variance
ambiguous nodes. `d_train = log(f_g^cur · M/E)`; additive, unit mass, bandwidth `h`.

### 3.4 Gate + harness
Config flag `gdna_prior_anchor: bool` (default off ⇒ byte-identical). A/B on
`scripts/debug/gdna_hyperprior_eval.py`. **Metric shift** (per your reframe): success is *gentle gravitation
toward the nearest mode with honest variance*, per library — not nailing a value. Report, per node type:
- node 1055 (gdna100): anchored **up** from the messages' 0.375 toward the enriched mode, **variance stays
  wide** (not collapsed) — a nudge, not a bulldoze.
- none_* (all ss): anchored **down** (that library's landscape has no enriched mode) — no over-call.
- AMBIG nodes (both-strand, the two-root case the prior is *for*): moved toward truth without collapsing `var`.
- ss0.99 stranded: unchanged.

### 3.5 Bandwidth / floor (open knobs, after the projection lands)
- **Bandwidth `h`**: start fixed (0.15); then adaptive (per-node counting precision — "one unit of mass, spread
  by precision"). `h` sets both the landscape smoothness and the projection's `N(d;μ,h)` width.
- **Floor**: the depleted background mode from the intergenic invariant (§0) is the natural low anchor; no
  separate ε-floor needed once the projection is an anchor (it can't crush).

---

## 4. Acceptance
- The anchor **gravitates** each node's DNA level toward the nearest landscape mode, **gently** (messages
  dominate; `var` reflects identifiability — wide at genuine ambiguity).
- node 1055 ↑ (real DNA library), none_* ↓ (no-DNA library), ss0.99 unchanged, AMBIG improved — **all from one
  mechanism**, no floor, no log-shift, no per-case tuning.
- Role-B path readable in one sitting: additive fit on solved nodes → `project` → gentle Gaussian anchor on the
  retained reference → read-out, both solvers. Role A untouched.

---

## 5. What we discard from the prior drafts
The δ-pin, the Stage-0 `logaddexp` floor, the enrichment-conditioning eff-scaling, and the RNA prior `π_r` — all
were solving a different problem (making the prior *decide* `f_g`). The prior's job is to **anchor**, and the
right anchor is the sampling-likelihood projection above.

## Appendix — evidence & prior docs
`gdna_crush_dissection_node1055.md` (the δ-pin is what crushes; swap-the-prior proof), `gdna_hyperprior_from_scratch.md`
(the gravity principle — now realized as §1.2), `scripts/debug/{proximity_projection_test,gdna_hyperprior_eval}.py`,
`docs/figures/*_fit_vs_oracle.png`.
