<!-- title: Count-Space Solver — Design -->
# Count-space solver: design (pseudo-count messages)

**Status:** design proposal — the convergence target for the solver's *precision* architecture. Builds on
the brainstorm + reviewer blueprint in `count_space_solver_idea.md`; this doc adopts the core, **corrects
one real flaw**, and scopes precisely what the pivot does and does not fix.
**Thesis:** a message is not a Gaussian penalty on a fraction — **it is a bundle of pseudo-fragments** at
the imputed composition, entering the *same count likelihood* as the real strand reads. Counts are the
currency end to end. This is the native form of everything we have been approximating (pseudocounts,
source-count ceilings, variance floors).

---

## 1. Why fraction-space overshoots (the reviewer's framing, which is correct)

The `idea.md` reviewer pinned the root cause cleanly: `σ²_bio` is **depth-independent biology** (how much
true template density fluctuates between adjacent regions), but the message precision multiplies
`1/σ²_bio` by the `(M/E)²` Jacobian. So:

```
strand Fisher  I(t) ∝ M           (linear in mass/depth)
message prec   τ_f  ∝ M²          (quadratic — via the Jacobian)
```

Because `M²` always outgrows `M`, **any deep/enriched node lets its messages steamroll the strand**, no
matter how we patch the variance. That is the structural reason the floor is a band-aid: it caps a
quantity that is in the wrong units. The cure is to make messages scale **linearly** in counts, like the
data.

## 2. The pivot — messages as pseudo-counts

Replace, in `_local_loglik`:

```
   −½·τ_f·(f − μ_f)²        →        Σ_c  α_c · log f_c          (c ∈ {pos, neg, g})
```

a Dirichlet/multinomial term whose mode sits at the source composition. It is in the *same units* as
`L_strand` (log-likelihood of discrete counts), so the two compete one-for-one. A zero channel
(`α_c = 0`) contributes `0·log f_c = 0` and drops out — no boundary anti-pinning, no Jacobian, no floor.

## 3. `N_pseudo` — the pseudo-count capacity (CORRECTING the blueprint)

How many pseudo-fragments may a source send? The blueprint proposes `N_pseudo = ρ²/σ²_bio = 1/CV²_bio`.
**That omits the source's sampling noise** and is the one real flaw: it credits a 1-fragment source and a
1000-fragment source at the same density with identical capacity — reinstating the "low count looks
certain" pathology this whole effort exists to kill. Keep the Poisson term:

```
N_pseudo = ρ_src² / ( σ²_bio(ρ_src) + ρ_src/E_src )          ← the effective count (1/CV²_total)
α_c      = N_pseudo · f_c[src]                               (partition by source composition)
```

This is `1/CV²` of the source's *total* uncertainty (biological + sampling), and it interpolates correctly:

| regime | `N_pseudo →` | meaning |
|---|---|---|
| sampling-limited (low count) | `ρ²/(ρ/E) = ρ·E = C_src` | the source's actual fragment count |
| biology-limited (high count) | `ρ²/σ²_bio = 1/CV²_bio` | the blueprint's hard ceiling (~100) |

A 2-read flank sends ≤ ~2 pseudo-counts; a million-read exon sends ≤ `1/CV²_bio`, no matter how deep. Use
the source's **observed** total count for the Poisson term so `N_pseudo` is frozen (no live-sweep feedback
— the reviewer's cascading-collapse hazard).

## 4. Effective length needs NO message-side normalization (the elegant part)

A subtle, correct point from the `idea.md` reviewer: **do not normalize the pseudo-counts by `E_dst`.**
The destination integrates `α` alongside its own real reads `N_dst`, and `N_dst` already scales with
`E_dst` (a longer node catches more fragments). So effective length balances the vote *automatically*
through the real-count term:

- huge dest (`N_dst = 1000`) + message (`N_pseudo = 100`) ⇒ local data holds ~91% (correct: deep node
  needn't be bossed around);
- tiny dest (`N_dst = 100`) + same message ⇒ ~50/50 (correct: a Poisson-noisy node leans on neighbours).

The graph self-balances in counts. No Jacobian anywhere — that term simply ceases to exist.

## 5. The whole objective unifies into counts

- **Global gDNA prior** → a gDNA pseudo-count (the review's `M`-free "Alt-A win").
- **Spliced lower bound** → a *hard* `+`/`−` RNA pseudo-count: observed spliced reads literally *are*
  strand-known RNA fragments — inject them as `α` directly (replaces the soft clipped-Gaussian).
- **Jeffreys Beta(½,½)** → `α_c = ½`.

Everything is `log-likelihood of (real + pseudo) counts on the simplex`.

## 6. What it fixes — and, honestly, what it does NOT

**Fixes (the precision dimension — the source of the whack-a-mole):**
- No `(M/E)²` Jacobian ⇒ no `M²` overshoot, no `1e9` spikes, no floor/cap to tune.
- Information symmetry: a neighbour competes with local reads one-for-one, capped at `N_pseudo`.
- Boundaries: `α_c=0` drops out ⇒ no anti-pinning; the Dirichlet repels `f_c` from 0 *only* when the
  source actually has component `c` (correct imputation; the zero-gDNA pin holds iff all neighbours are 0).

**Does NOT fix (orthogonal — do not expect the pivot to solve these):**
- **reg242 / capture-blind global *mean*.** Pseudo-counts change the *currency*, not the *mean*. If the
  global baseline is measured on un-enriched regions it is still wrong. Needs the capture-aware
  `ρ_baseline` + local-enrichment propagation — a separate track.
- **The balanced-AMBIG `f_g` ridge (identifiability).** A balanced count is still equally gDNA or
  balanced-RNA; pseudo-counts resolve it by *propagation* exactly as the Gaussian messages do — no better.
  The genuine resolver is the node-intrinsic **gDNA-FL vs RNA-FL** signal (deferred for accumulator cost).

The pivot is the cure for *precision* and nothing more — which is exactly the dimension that has been
generating the per-condition trade-offs.

## 7. Implementation scope — a term swap, not a reparameterization

Unlike ALR (which re-parameterizes the simplex and was rejected for its `f_g→0` singularity), this keeps
the lattice, the directional sweep, the posterior-median readout, and the strand likelihood **unchanged**:

1. `simplex_sweep._local_loglik`: swap the Gaussian imputation penalty for `Σ_c α_c·log f_c`
   (pre-compute the static `log f_c` lattice, clamped at `_PRIOR_EPS`).
2. `bp_solver._message`: return `(f_c[src], N_pseudo)` instead of `(μ_f, τ_f)`; delete the Jacobian and the
   variance floor; `N_pseudo` per §3.
3. global / spliced / Jeffreys priors → pseudo-counts (§5).
4. `variance_model.py`: docstring only (it now scales pseudo-counts, not precisions).
5. `priors.py` factor-1-under-uniform invariant: add a test with `f_g` from the new solve **before**
   golden regen (the one place a calibration change reaches the EM).

## 8. Open questions

1. **Per-component vs total `N_pseudo`.** `N_pseudo` is a total-density property, but `σ²_bio` is fit per
   component today. Fit a total-density `var~mean`, or combine the component curves.
2. **Readout.** The posterior is now strand × Dirichlet; verify median-`f_g` stays well-behaved (wide
   posterior when unanchored — the honest answer at AMBIG).
3. **Convergence.** Confirm the count-accumulating sweep contracts (frozen `N_pseudo` is what guards it).
4. **Spliced as hard pseudo-counts** vs the existing spliced floor — replace or complement?

## 9. Recommendation & sequencing

The count-space pivot is the right precision architecture — your instinct, the review's endorsement of the
Alt-A "narrow win," and the count-currency derivation all converge, and it is *more tractable than feared*
(a term swap). Sequence:

1. **Interim (done):** the binomial variance floor — ship-grade precision control while the pivot is built.
2. **Build the pivot** (§3 corrected `N_pseudo` + §5 pseudo-count priors).
3. **In parallel, separate tracks:** the capture-aware global **mean** (reg242) and the **FL** ridge
   resolver — neither is solved by the pivot and both still matter afterward.
