<!-- title: Honest-precision imputation messages — audit + remaining design -->
# Honest-precision imputation (belief propagation) — audit + remaining design

**Agreed design:** a message is **(density, precision)** — a density `ρ` with its *accompanying* precision,
not a density alone. The precision must be **honest**: in count/effective-count units, so a message competes
**one-for-one** with the destination's own observed fragments (emergent deference), never steamrolling them.

## 1. Audit — what is implemented, and the gap (causally confirmed)
`_message` (bp_solver.py) **computes** the honest density pair internally:
- density `ρ_src` (the source's facing density, from `node_densities`);
- density-precision `τ_ρ = 1/(σ²_bio(μ) + ρ_src/E_src)` — the learned biological between-node spread + the
  source's Poisson sampling noise. This is count- and length-aware and **source-bounded** (the `ρ_src/E_src`
  term caps it by the source's own evidence). **This part is honest.**

But it does **not deliver** `(ρ, τ_ρ)`. It converts to a fraction prior and returns `(μ_f, τ_f)`:
```
μ_f = (ρ_src·E_dst − spliced_dst)/M_dst        # a FRACTION
τ_f = τ_ρ · (M_dst/E_dst)²                       # the DEST-MASS Jacobian
τ_f = 1/(1/τ_f + f(1−f)·(1/M_src + 1/M_dst))     # binomial floor — VANISHES at μ_f→0/1
```
and `_local_loglik` applies it as a Gaussian `−½·τ_f·(f_c−μ_f)²`. **Two dishonesty sources:**
1. **The `(M_dst/E_dst)²` Jacobian scales the precision by the DESTINATION mass², not the source's evidence.**
   The message's information is `N_src ≈` the source's fragment count; the delivered fraction-precision instead
   grows like `M_dst²`, so a high-mass destination makes the message look far more confident than `N_src`
   warrants.
2. **The source-evidence bound (binomial floor) vanishes at the simplex walls** (`f(1−f)→0`). Exactly where a
   message asserts "component ≈ 0/1" the cap disappears, so `τ_f = τ_ρ·(M/E)²` is uncapped.

**Causal confirmation** (node463, stranded zero-gDNA): a +RNA message `μ=0.00, τ=1e7` overrode a node with
upos=16336, uneg=149 (its own data says f_pos≈1) → `f_pos→0 ⇒ f_g→1` phantom. The message — worth a few source
fragments — beat 16k of the node's own evidence because its delivered precision (1e7) was a `(M/E)²` artifact,
not `N_src`. Capping the message precision at the dest count collapsed the stranded-0 phantom 27,933→1 (and the
unstranded-0 16,420→9), proving the over-delivered message precision is the cause. (Same `(M/E)²` Jacobian that
afflicts the global prior — see `global_prior_count_space_derivation.md` §6b.)

## 2. The honest representation — message = (density, EFFECTIVE-COUNT)
The accompanying precision must be the source's **effective count** `N_src` (in fragments), not a
fraction-precision. Then the message enters the destination's solve as `N_src` **pseudo-fragments** at the
imputed composition, competing one-for-one with the destination's `M_dst` real fragments.

```
ρ_src                                              the message density (unchanged)
N_src = ρ_src² / (σ²_bio(μ) + ρ_src/E_src)         the source's effective count = 1/CV²_total = ρ_src²·τ_ρ
```
`N_src` is the **count-currency** capacity (count_space_solver_design.md §3): a sparse or biologically-uncertain
source carries few pseudo-fragments; a deep, homogeneous source carries many — but never more than its own
evidence. Crucially it is **independent of `M_dst`** — no dest-mass amplification.

**Application (Jacobian-free):** the message becomes a Beta/Dirichlet pseudo-count on the destination component
`c`, with the imputed fraction `μ_c = clip(ρ_src·E_dst/M_dst, 0, 1)`:
```
log-term:  α_c·log(f_c) + (N_src − α_c)·log(1 − f_c),     α_c = N_src · μ_c
```
evaluated on the f lattice. The curvature is set by `N_src` (the source's evidence), **not** `(M_dst/E_dst)²`.
At node463: `N_src` (the neighbour's few fragments) ≪ `M_dst=16485` ⇒ the node's strand likelihood dominates ⇒
the message defers ⇒ no phantom. At a genuine clean crossing (capture), the source's real `N_src` is large ⇒ it
appropriately pins. **Deference is emergent**, exactly as for the count-space global (already built behind
`RIGEL_GLOBAL_ALLNODES` as the Gamma-on-count term — this is the same pattern for the messages).

## 3. Remaining implementation pieces
1. **`_message` returns `(μ_c, N_src)`** instead of `(μ_f, τ_f)`: keep `μ_c = clip(ρ_src·E_dst/M_dst, 0, 1)`;
   compute `N_src = ρ_src²/(σ²_bio + ρ_src/E_src)`; **delete the `(M_dst/E_dst)²` Jacobian and the
   wall-vanishing binomial floor** (the `N_src` bound replaces them — it is the honest source-evidence cap and
   does not vanish at the walls). The spliced subtraction stays in `μ_c`.
2. **`_local_loglik` applies the gDNA + per-strand-RNA messages as count-space pseudo-counts** (Beta/Dirichlet
   `α_c·log f_c + (N−α_c)·log(1−f_c)`), replacing the three Gaussian `−½·τ·(f−μ)²` imputation terms. Reuse the
   same machinery as the global Gamma-on-count term so all priors are one currency.
3. **Unify with the global** (already count-space behind the flag): messages + global all enter as
   pseudo-counts ⇒ they and the strand likelihood (the real counts) are all log-likelihoods of counts on the
   simplex, competing one-for-one. This is the count-space local solve; it removes the `(M/E)²` over-pinning
   **everywhere** — the unifying fix the global-only change pointed to.
4. **Per-strand bookkeeping:** the RNA messages are per strand (f_pos, f_neg); each carries its own `N_src,s`
   and `μ_s`. The `free_s` continuity gate (τ=0 today) becomes `N_src,s=0` (no pseudo-count) — identical effect.
5. **Relay door (future):** with messages as pseudo-counts, forwarding through a transparent node is "pass the
   incoming `N` along" — the door the relay effort needs (kept open, not built here).

## 4. Validation plan
- **Unit (count-currency):** a message with `N_src` pseudo-fragments must move a destination's f_c by no more
  than `N_src` real fragments would — assert the node463 case (`N_src ≪ M_dst` ⇒ f_pos stays ≈1).
- **The 3 scenarios:** stranded-0 and unstranded-0 phantoms → ~0 (the message no longer overrides the strand);
  capture preserved (real `N_src` at clean crossings still pins); + `factor-1-uniform`; + the full suite.
- **No-magic-numbers:** the only constant is the shared pseudocount `a = _MSG_PSEUDOCOUNT = 1` (already
  decided); `N_src` and `μ_c` are derived.

## 5. Status
The (density, precision) message is honest in its **density-precision computation** but dishonest in its
**fraction delivery** (the `(M/E)²` Jacobian + wall-vanishing floor). The remaining work is to deliver the
precision as the **source effective-count `N_src`** and apply messages as **count-space pseudo-counts** —
unifying with the count-space global already prototyped. This is the message layer of the count-space solve, and
it is the precise fix for the over-pinning root cause (causally confirmed in §1).
