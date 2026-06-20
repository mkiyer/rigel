<!-- title: Count-space local solve — implementation plan + deep evaluation -->
# Count-space local solve (honest-precision messages + global) — implementation plan + evaluation

> **VERDICT (adversarial readiness review, 6 lenses): NOT READY — SUPERSEDED. Do not execute.** Three
> decisive findings: (1) the **Dirichlet `Σα·log f` form is wrong** — a one-sided `α_c·log(f_c)` is a MONOTONE
> push to the f_c=1 vertex, not a pull toward μ_c (measured: a "μ_g=0.2 / gDNA-is-low" message → median
> f_g≈0.77, the opposite of its content). **Keep the two-sided `(density, precision)` Gaussian pull** and fix
> only its precision. (2) This plan **contradicts the user-LOCKED plans-of-record** (`count_space_relay_*` D1,
> `count_space_execution_plan` invariant #3) which adopt `(density, precision)` and explicitly reject the
> Dirichlet; the `(n+a)/E²` message precision is **already shipped** (`dfc3360f`), and the binomial floor this
> plan proposed to delete is **load-bearing** (capping it cost capture −6%, reverted). (3) The motivating
> regression is **stale**: node463/stranded-0 is already **0** on the production default path (`225d65ac`, the
> intrinsic-solve); the +28K is a `RIGEL_GLOBAL_ALLNODES`-ON (default-off) artifact. The one genuinely-open
> phantom (unstranded-0 at κ≈½) is an **identifiability floor** the prior form structurally cannot cross — it
> needs the **FL signal**, not a prior rewrite. **Remaining real honest-precision work = the σ²_bio raw-density
> fit** (OQ-B / execution-plan §1c). Sections below retained as the record of the rejected approach.

**Goal:** make every prior in the per-node solve a **count-space pseudo-count** so messages, the global, and
the strand likelihood are all log-likelihoods of counts on the simplex, competing one-for-one (emergent
deference, no `(M/E)²` Jacobian). This is the unifying fix for the over-pinning root cause
(`honest_precision_message_design.md`, `global_prior_count_space_derivation.md` §6b).

## 1. The unified form — DIRICHLET pseudo-counts (corrects the prior Beta note)
Every prior adds per-component pseudo-fragments to the node's count pool; the log-prior is
```
Σ_c  α_c · log(f_c) ,   c ∈ {pos, neg, g}
```
a Dirichlet on the simplex, in the **same units** as the strand multinomial likelihood. (NOT the Beta
`α·log f + (N−α)·log(1−f)` from the earlier doc — that double-counts the complement and mis-handles the 3-way
simplex.) Key property: a zero channel `α_c=0` contributes `0·log f_c = 0` and **drops out** — no anti-pinning,
no override (this is what fixes node463: a `μ=0` +RNA message adds `α_pos=0` ⇒ vanishes ⇒ the strand wins).

## 2. The message — (density, effective-count), applied as a pseudo-count
`_message` returns `(μ_c, N_src)`:
```
μ_c   = clip((ρ_src·E_dst − spliced_dst)/M_dst, 0, 1)      the imputed fraction (unchanged)
N_src = ρ_src² / (σ²_bio(μ_q) + ρ_src/E_src)               the source effective count (1/CV²_total)
```
applied as `α_c += N_src · μ_c` for the message's component `c`. **Delete** the `(M_dst/E_dst)²` Jacobian and
the wall-vanishing binomial floor — `N_src` is the honest, non-vanishing source-evidence bound. gDNA message →
`α_g`; +RNA → `α_pos`; −RNA → `α_neg`. `free_s`-gated channels contribute `N_src=0` (no pseudo-count).

## 3. The global — also a pseudo-count, with the zero-pin caveat (OPEN)
The exposure-pooled rate `ρ_global=(Σĝ+a)/(ΣE+b)` over self-solvable nodes (validated: zero→0 by intergenic-
length dilution; capture→genome-average). As a Dirichlet pseudo-count: `α_g += N_global·μ_global`,
`N_global = ρ_global²/σ²_n`, `μ_global = clip(ρ_global·E_j/M_j,0,1)`.

> **OPEN QUESTION A — the global cannot PIN f_g→0 in the pure Dirichlet form.** In zero-gDNA `μ_global→0` ⇒
> `α_g=0` ⇒ the global drops out; f_g→0 then relies on the **RNA messages** (α_pos/α_neg from RNA-bearing
> neighbours pushing the RNA fractions up, squeezing f_g down). For an **isolated, balanced-count AMBIG** node
> (no RNA neighbour, no strand tilt) nothing pins f_g→0 — this is the genuine κ→½ identifiability floor (FL is
> the only resolver), so dropping out may be the *honest* behaviour. BUT it differs from the Gamma-on-count
> global already prototyped (whose rate term `k/g_pseudo→∞` actively pins g→0). **Must measure** whether the
> pure-Dirichlet global regresses the isolated-balanced-AMBIG nodes vs the Gamma form; if so, keep a rate-style
> zero-pin term for the global only (a principled exception, not full unification).

## 4. Migration steps (each measured against the 3 scenarios + suite)
1. `_message` → `(μ_c, N_src)`; drop the Jacobian + binomial floor.
2. `_local_loglik` → replace the 3 Gaussian imputation terms + the global term with the Dirichlet
   `Σ_c α_c·log(f_c)`, accumulating α_c from gDNA msg + ±RNA msgs + global + Jeffreys (`α_c += ½`).
3. The strand BB likelihood, the spliced floor, the lattice, the median readout — **unchanged**.
4. Retire the Gaussian path + the `RIGEL_GLOBAL_ALLNODES` flag once green (or keep flag through validation).

## 5. DEEP EVALUATION — issues I can already see (before execution)
1. **(form) Dirichlet is add-only ⇒ a message cannot pull a component DOWN, only up.** A message saying "this
   node's +RNA is LOW" contributes `α_pos≈0` (no force); f_pos only falls if *other* components' pseudo-counts
   rise. Usually fine (the strand / other messages compensate), but a node with spuriously-high f_c and a
   confident "low" neighbour will not be corrected directly. Risk: residual over-calls where a component should
   be pulled down. Needs a scenario.
2. **(global unification) §3 OPEN-A** — the zero-pin tension; the global may need a non-Dirichlet term.
3. **(precision honesty) `N_src` is only as honest as `σ²_bio`.** If the var~mean under-estimates the
   between-node spread, `N_src` is too large ⇒ the same over-confidence returns (now in count units). The
   var~mean honesty (raw-density fit, OQ4 from the FB review) is a *prerequisite*, not yet done.
4. **(competition scale) the strand likelihood's gDNA-vs-RNA information is `M·(2κ−1)²`, not `M`.** Emergent
   deference is "N_src vs M·(2κ−1)²": correct in principle (the strand wins only where informative), but the
   Dirichlet curvature `α_c/f_c²` vs the BB curvature must actually balance — verify a resolved node isn't
   over- or under-ridden across κ.
5. **(simplex coupling) three independent per-component pseudo-counts on a 2-simplex.** `Σα_c log f_c` is a
   proper Dirichlet only if the α_c are jointly the concentration; here each α_c comes from a *different* source
   (gDNA genomic, ±RNA strand-continuous). This is fine as a product-of-independent-priors but the normalization
   interplay with the strand multinomial needs a check (does the lattice posterior stay proper / no edge spikes
   at `_PRIOR_EPS`?).
6. **(numerics) `log(f_c)` at the f_c→0 lattice edge** is clamped at `_PRIOR_EPS=1e-3`; with large α_c the
   `α_c·log(_PRIOR_EPS)` term is a large constant that could dominate / distort the median near the wall.
   Re-examine the edge handling under the new (potentially large) α_c.
7. **(relay door) OQ1 conflict.** We rejected the Dirichlet message form for *relay* (it conflates density &
   precision — can't forward "same density, lower confidence"). The per-hop GS sweep doesn't relay, so the
   Dirichlet is fine *now*; but piece-5 "keep the relay door open" is then NOT served by this form — forwarding
   would need the (density, precision) carrier, a different representation. The plan must state honestly that
   this closes the relay door it earlier promised to keep open (or reconcile).
8. **(equivalence/validation) no bit-stable fallback** — this changes the prior form for every node; goldens
   shift. Need the 3-scenario + suite gates and the node463 unit test, not goldens-unchanged.

## 6. Open questions (consolidated)
- **A.** Global zero-pin in the pure Dirichlet (§3) — measure; keep a rate term if needed.
- **B.** `σ²_bio` honesty (raw-density var~mean fit) — prerequisite; is it in scope here or a separate fix?
- **C.** Add-only Dirichlet can't pull a component down (§5.1) — acceptable, or do we need a complement term?
- **D.** Relay door (§5.7) — does this form abandon it? If relay matters, the carrier must stay (density,
  precision), contradicting the Dirichlet application. Decide the precedence (honest-now vs relay-later).
- **E.** The `a`,`b` Gamma-prior + `_MSG_PSEUDOCOUNT` ledger through no-magic-numbers.

## 7. Readiness
**Not yet ready to execute as a single drop-in.** The form is right (Dirichlet pseudo-counts) and the core
fix (node463 over-pin) is well-justified, but four substantive design questions are unresolved (A: global
zero-pin; C: add-only down-pull; D: relay-door precedence; B: σ²_bio honesty prerequisite). Recommend: resolve A
& D by decision, B by sequencing, then execute the migration §4 sub-phase by sub-phase with the per-step gates —
not a big-bang. The next section is an adversarial review to surface what this self-evaluation missed.
