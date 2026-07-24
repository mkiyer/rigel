# Variance foundation — the PROPOSAL (derived + numerically validated)

**Status:** proposal for review, 2026-07-24. Supersedes the incomplete `variance_foundation_plan.md` (its τ_θ
diagonal state was refuted). Produced by a 5-approach derivation workflow, each candidate **numerically
validated** against reference-free grid posteriors + Monte-Carlo, then adjudicated. The load-bearing result is
**exact algebra, not an approximation**, and all five approaches independently corroborate it.

---

## 1. The problem (recap)

A node's composition `(f_g, f_+, f_−)` is 2-DOF; coordinates `λ = logit f_g` (gDNA level), `θ = arcsin τ`
(strand tilt). The strand Beta-Binomial likelihood depends on `(λ, θ)` **only** through the scalar
`p = ½ + (κ−½)(1−f_g) sin θ` — so the strand Fisher information is **rank-1**: `I_strand = c·∇p∇pᵀ`,
`c = N_eff/(p(1−p))`. A diagonal `(τ_λ, τ_θ)` state cannot represent a rank-1 constraint — it either drops the
covariance or **double-counts** the one scalar the strand actually measures. That refuted the drafted plan.

## 2. The resolution — the Schur-reduced SCALAR gDNA-level precision `τ_λ`

The honest **marginal** gDNA-level precision is the Schur complement of the composition Fisher
`J = I_strand + diag(I_density, 0)`:

```
    τ_λ  =  J_λλ − J_λθ²/J_θθ  =  (c·a² + I_d) − (c·a·b)²/(c·b²)  =  I_d           (θ a FREE nuisance)
```

with `a = ∂p/∂λ = −(κ−½) sinθ · f_g(1−f_g)`, `b = ∂p/∂θ = (κ−½)(1−f_g) cosθ`, `I_d` = the NB
density-deconvolution (intron-factory) curvature on `λ` (0 on a pass-0 exon with no gDNA prior). **The strand
cancels out of `f_g` exactly.** With a tilt prior/message of precision `r_θ` on `θ`:

```
    τ_λ  =  τ_d  +  c·a² · r_θ/(c·b² + r_θ)                                        (general)
```

**Pass-0 has exactly two limits, nothing between:**

| node | why | `r_θ` | `τ_λ` |
|---|---|---|---|
| **single-strand (1-DOF)** | θ is STRUCTURALLY locked (`f_−≡0`) — not a free nuisance | ∞ | `τ_d + c·a²` (strand **pins** f_g) |
| **AMBIG (2-DOF)** | θ is a free nuisance | 0 | `τ_d` (strand gives **ZERO** f_g precision, for all N) |

Equivalently (the exact proof, approach D): reparametrize the tilt to the strand's actual identified scalar
`t = f_+ − f_−`. In `(λ, t)` the atomic likelihood is **exactly block-diagonal** — `∂²logL/∂λ∂t ≡ 0`, no
linearization — so `τ_λ = I_λ + I_t·1[single-strand]` with `I_t = N_eff(κ−½)²/(p(1−p))` the strand precision on
`t`. The `(λ,θ)` "coupling" is a coordinate artifact (grid correlation **0.64 in `(λ,θ)` → 0.035 in `(λ,t)`**,
and **2e-17** for a density-only node). **This is why the Schur collapse is exact, not a Laplace
approximation.**

## 3. It exposes and fixes a LIVE bug (real, but bounded — corrected by the independent verification)

The current `node_init.strand_evidence` computes `I_strand = c·(κ−½)²[f_g(1−f_g)]²` (the single-strand λ-Fisher
`c·a²`) and adds it to **every** solvable node's `τ_λ` — i.e. it credits the **single-strand** λ-term to
**AMBIG** nodes too, where the honest contribution is **zero** (the strand cancels out of `f_g`, §2).

**Independent verification (`scratchpad/verify_foundation.py`, ALL claims PASS):** the bug is REAL — for a
stranded AMBIG node `strand_evidence` returns `τ_λ = 1.33` where the honest value is `0`. **Severity correction
(vs the workflow's Poisson toy):** the workflow claimed a ∝N / 2.8×–1000× over-credit, but that is the *no-
overdispersion* limit. Production caps the count power via the overdispersion `n_str = N/(1+(N−1)ω) → 1/ω`, so
the AMBIG phantom **SATURATES** (measured `τ(2N)/τ(N) = 1.04`) — a **bounded, weak** phantom (`~1/ω·disc·jac`),
not an unbounded one. The reference-free grid marginal λ-precision for AMBIG-strand-only is correspondingly
**N-invariant** (0.028→0.054 across N=50→3200) while the *uncapped* `I_strand` grows ∝N (2.4→152). So this is a
real correctness defect on the AMBIG nodes calibration exists to resolve — worth fixing for correctness and the
clean skeleton — but **not** the dominant error source (consistent with the healthy stranded arm, mwae 0.026).

The fix is **surgical and simplifying**: the production `own_composition_logvar` already maps one scalar
`τ_lam` to both fraction arms, so **the architecture is already this proposal (approach E)**. The only change is
to **gate the `c·a²` strand λ-term to single-strand nodes** (`free_pos XOR free_neg`); AMBIG gets `τ_λ = τ_d`.

## 4. Numerical validation (reference-free grid + Monte-Carlo; rigel env)

* `τ_λ` matches a reference-free brute-force grid posterior of the **atomic** (BB strand × NB density) likelihood
  to **<10% everywhere**: single-strand 16.34 vs 15.81 (1.05); AMBIG+density 8.75 vs 8.23 (1.09) **while current
  code = 23.04 (2.8× overconfident)**; AMBIG-unstranded → 0 (correct).
* Schur identity confirmed to **machine precision** (`marginal-λ = I_density`, 5.3319 vs 5.3319).
* Rank-1/singular structure confirmed: `det(J) = 2.7e-15` for AMBIG-no-density, top eigenvector = `∇p` to 6
  digits (the node is honestly unsolved along the ridge, solved along `p`).
* MC (2000–4000 draws) cleanly **separates belief precision (N-invariant, correct) from point-estimate sampling
  (∝1/N, irrelevant to the state)** — the count-zero-info principle, made numerical.

## 5. The proposed per-node state (minimal)

Nothing new at pass-0 beyond **one scalar**:

1. `f_g` (= `σ(λ)`) — the gDNA-level MODE (existing ψ solve; unchanged).
2. `θ` (tilt point) — to split RNA into ± (existing `NodeDeconv.theta_mean`; unchanged).
3. **`τ_λ`** — THE single foundational composition precision (Schur-marginal λ). The only genuinely new number;
   it **replaces** the mislabeled `tau_lam`. Emitted to the eventual message layer via the exact Jacobians
   `Var(log f_g) = (1−f_g)²/τ_λ`, `Var(log f_r) = f_g²/τ_λ`.
4. `struct_lock` + `free_pos`/`free_neg` — already in `NodeStatics`, applied unchanged. **A struct_lock node is a
   HARD OVERRIDE in fusion, never an ∞ weight** (an ∞ poisons the interior chain with nan).
5. per-component counts `n_g, n_pos, n_neg` — carried **separately** for the message layer's `1/n`; **never
   folded into `τ_λ`.**

`τ_θ` (tilt precision) is **deferred** — it is not part of the gDNA-vs-RNA composition foundation; it re-enters
`λ` only as `r_θ` when a genuine tilt message arrives (the message-layer task).

## 6. Why this is the right skeleton (foundation ⊕ messages, unified)

The **same Schur formula IS the message-combination rule**: `r_θ` is the incoming tilt-message precision (0 at
pass-0, >0 when a tilt message arrives), and `1/n` is added only to the **transported density** at message time.
So the foundation and the message layer share **one derivation**; the foundation is the `r_θ=0`/no-`1/n`
instance. This is exactly "message-propagation precision built on top of a solidified foundation." The full 2×2
information matrix (approach A) is the general fusion algebra of which this scalar-per-axis form is the
rank-reduced instance — kept as the **conceptual model**, not stored state.

## 7. What changes in code (surgical)

* `node_init.strand_evidence` → gate the strand `c·a²` λ-term to single-strand nodes (the bug fix); AMBIG τ_λ =
  density curvature only. Fold the `sinθ` factor correctly (currently dropped).
* `own_composition_logvar`/`own_precision` → keep (they already map the scalar to the arms); **remove the `1/n`
  from the composition precision** (it moves to the message layer — the separation).
* No 2×2 matrix, no covariance storage. The state SHRINKS.

## 8. Open items (carried, not blockers)

* **Pass-dependence:** `τ_λ` is evaluated at the prior-shifted mode `fg_loc`, so it differs pass-0 vs refit. It
  is not "reference-free"; document for the successor.
* **The single-strand-vs-AMBIG gate is now load-bearing** — add a regression test that an unresolved AMBIG
  strand emits **zero** own f_g precision at every N.
* **Measured-spliced boundaries** self-solve to `τ_λ=0` on unstranded data (own RNA precision is strand-only).
  Either promote the measured spliced count to a Poisson `1/n_spliced` own-RNA belief, or explicitly document
  that boundaries ride the relay graft. **State it, don't silently defer.**
* **Validation is NOT byte-identical** (single-strand unchanged; AMBIG changes by design — it fixes the bug).
  Needs the falsifiable per-condition A/B (refit=0 AND refit≥1; unstranded-capON called out).
* **Laplace mapping** `τ_λ → Var(log f_c)` is 2nd-order; near the `f_g→1` vertex it under-summarizes ~10–18%
  (shared by all approaches and the current code) — not a blocker.
* **verify-agent caveat:** the dedicated math-verifier sub-agent returned a degenerate stub; verification is
  instead covered by the decision agent's hand-derivation + the five approaches' mutual numerical corroboration
  + the by-hand Schur identity above. A fresh independent verification pass is cheap insurance before coding.

## 9. Recommendation

**Adopt approach E** (the Schur-reduced scalar `τ_λ`), derived by D's exact `(λ,t)` block-diagonalization, with
A's 2×2 info-form as the conceptual fusion algebra. It resolves the rank-1 blocker, **shrinks** the state to one
honest scalar, fixes a live 2.8×–1000× AMBIG over-credit, and unifies the foundation with the message layer.
Next: revise `variance_foundation_plan.md` to build E (the surgical strand-gate + the 1/n separation), then
implement with the falsifiable per-condition gate.
