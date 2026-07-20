# Message transfer variance in a belief-propagation deconvolution — a formal derivation for review

**Status:** formal derivation prepared for external statistical/mathematical review, 2026-07-17. **Self-contained
— the reviewer is not assumed to know the codebase.** The question: we pass Gaussian messages along a chain to
deconvolve a mixed signal; each message needs a variance. We have built (a) a population density model (an
NPMLE) and (b) a per-node "projection" giving each node a posterior mode and variance for its latent rate. We
want to know, rigorously, **how the message variance should be assembled from these pieces** — and in particular
whether our formula is correct, or double-counts, or omits a term. Companion (design + empirics):
`npmle_projection_variance_design.md`; the empirical validation is in §8.

---

## 1. The problem

At each of ~10⁶ genomic **nodes** we observe a total fragment **mass** `M` (a count) over a known **effective
length** `E`, plus a strand split `(u₊, u₋)` of that mass. Each node's mass is a mixture of components — for this
document, treat it as two: **gDNA** and **RNA** (the real problem has RNA₊/RNA₋; the two-component case carries
the whole argument). We want each node's composition, summarised by `f_g ∈ [0,1]` = the gDNA fraction, so that
the component **densities** are `ρ_g = f_g·M/E` and `ρ_r = (1−f_g)·M/E`.

The nodes form a **chain** (a forest of linear genomic paths). Adjacent nodes' true compositions are related:
genomic gDNA density is spatially smooth *except* where a hybridisation-capture probe enriches one region and not
its neighbour (an **enrichment crossing**), and RNA density varies with expression. We exploit this by **belief
propagation** (BP): a node with a confident composition informs its neighbour. Because the graph is a forest of
linear paths, one forward + one backward Gaussian pass is exact (up to the Gaussian approximation).

**The estimand of this document is the variance of the message** one node sends its neighbour. Too small and a
dense node steamrolls a correct neighbour; too large and messages are inert. The count must enter *only* as
statistical power, never as a vote on composition (a domain invariant we hold to: a fragment count alone carries
no gDNA-vs-RNA information).

---

## 2. Notation

Work in **log-density** `u_c = log ρ_c` for a component `c`. Per node `n`:

| symbol | meaning |
|---|---|
| `u_c^n` | the node's **true** (latent) log-density of component `c` |
| `M_n, E_n` | observed total mass (a count), known effective length |
| `O_n` | all of node `n`'s own data (`M_n, E_n, u₊, u₋`) |
| `P(u)` | the population prior on the log-density (the NPMLE, §4) |
| `V_c^n` | node `n`'s **posterior** variance on `u_c^n` given its own data + prior (its "belief width") |
| `b_n(u_c)` | node `n`'s belief (a Gaussian `N(m_c^n, V_c^n)`) |

We fit `P(u)` **belief-free** on the observed total density `M/E` (asserting all mass is one component) — it is
the population landscape of enrichment, prior to any deconvolution (justification: capture enriches all nucleic
acid, so the *total*-density landscape is a valid, assumption-free profile of the enrichment geometry).

---

## 3. The Gaussian BP message (textbook, stated for reference)

For a pairwise Markov chain with node beliefs and an **edge potential** `ψ_{ij}(u_c^i, u_c^j)` encoding how
adjacent true log-densities co-vary, the message from source `s` to destination `d` is

```
    m_{s→d}(u_c^d)  =  ∫ ψ_{sd}(u_c^s, u_c^d) · b_s(u_c^s) / m_{d→s}(u_c^s) · du_c^s .
```

Take the edge potential to be a soft "adjacent true log-densities are equal", i.e. a Gaussian on their
difference with spread `σ²_T` (the **transport variance** — the biological/structural spatial dispersion):

```
    ψ_{sd}(u_c^s, u_c^d)  =  N( u_c^d − u_c^s ; 0 , σ²_T ) .
```

With the (cavity) source belief `N(m_c^s, V_c^s)`, the convolution gives a Gaussian message

```
    m_{s→d}(u_c^d)  =  N( u_c^d ; m_c^s , V_c^s + σ²_T ) .                           (★)
```

**So the message variance is `V_c^s + σ²_T`: the source's own posterior variance plus the transport variance.**
This is the object we must specify. `V_c^s` is *knowledge* (shrinks with data); `σ²_T` is a *world* property
(the spatial dispersion), count-independent. The rest of the document determines each from our pieces.

---

## 4. The two pieces we built

**(A) The population prior `P(u)` (NPMLE).** A fixed-kernel mixture `P(u) = Σ_j w_j·N(μ_j, h²)` on a grid of
locations `μ_j`, weights `w_j` fit by EM on the belief-free total density, bandwidth `h`. Under capture it is
multi-modal (a depleted off-probe mode, enriched on-probe mode(s)); off capture, unimodal.

**(B) The per-node projection.** Given a node's observed log-density `d_n = log(M_n/E_n)`, define responsibilities
over the mixture components and read a posterior mode and variance for the node's latent log-rate:

```
    r_j^n ∝ w_j · N(d_n ; μ_j , h²) ,   Σ_j r_j^n = 1
    μ_proj(n)  = Σ_j r_j^n · μ_j                                   (the denoised rate: snaps toward a mode)
    var_proj(n) = Σ_j r_j^n · (μ_j − μ_proj(n))²  +  h²            (between-mode ambiguity + within-mode h²)
```

`var_proj(n)` is small (`≈ h²`) when the node's density sits unambiguously in one mode, and large in a valley
between modes (mode-assignment ambiguity). It is **count-free** (uses only `d_n` and the fixed mixture).

---

## 5. Assembling the message variance

We must fill in `V_c^s` and `σ²_T` in (★).

### 5.1 `V_c^s` — the source posterior variance (the term BP already prescribes)

`V_c^s = Var(u_c^s | O_s, prior)`. In the solver, the source's belief on `f_g` is produced by a per-node solve
that combines its **strand likelihood**, the **NPMLE prior**, and its **incoming messages**; its posterior
variance on `log f_c` is `v_logfc`. Converting composition→density, `u_c = log f_c + log(M/E)` and `M` is a
count with `Var(log M) ≈ 1/M`, so

```
    V_c^s  =  Var(log f_c^s)  +  1/M_s   =   v_logfc  +  1/M_s .                      (5.1)
```

This is exactly the term textbook BP wants — **the source's own posterior variance**, and it already
incorporates the NPMLE prior (through the solve). *(This answers a reviewer-style objection we received: "aren't
you failing to use the source's variance-projection that BP would use?" We use `V_c^s`, the full posterior
(5.1). The belief-free `var_proj(s)` is a **different, NPMLE-only** quantity — a redundant proxy for the part of
`V_c^s` the prior already contributes. Using it in addition would double-count the prior; see §5.3.)*

### 5.2 `σ²_T` — the transport variance (the crux)

`σ²_T = ` the spread of the destination's true rate around the message's claim. The message's mean is the
transported source rate; its variance about the destination's **true** rate is

```
    σ²_T  =  E[ ( u_c^d,true − (message mean) )² ]
          =  Var(u_c^d,true)  +  ( E[u_c^d,true] − (message mean) )² .                (5.2)
```

Model the destination's true rate by projecting its observation onto the population mixture:
`u_c^d,true` has mean `μ_proj(d)` and variance `var_proj(d)`. Take the message mean at the source's denoised
rate `μ_proj(s)` (identity transport in the enrichment landscape). Then (5.2) becomes

```
    σ²_T  =  var_proj(d)  +  ( μ_proj(d) − μ_proj(s) )²   ≡   F1 .                     (5.3)
```

**Two terms, both endpoints present:** the destination's own mode-ambiguity `var_proj(d)`, and the squared
**mode gap** `(μ_proj(d) − μ_proj(s))²` — this is where the **source** enters `σ²_T`. Same mode ⇒ gap ≈ 0 ⇒
`σ²_T ≈ var_proj(d) ≈ h²` (the floor); an enrichment crossing ⇒ gap = the mode separation ⇒ `σ²_T` large ⇒ the
message is gagged.

### 5.3 Why `var_proj(s)` is **absent** from `σ²_T` (the double-count, made precise)

It is tempting to symmetrise (5.3) to `var_proj(s) + var_proj(d) + gap²` — call it **F2**. Indeed F2 is the exact
**second moment of the difference of two independent mixture draws**:

```
    E[ (u_c^d − u_c^s)² ]  =  var_proj(s) + var_proj(d) + (μ_proj(d) − μ_proj(s))²      (verified: §8)
```

(with `var_proj = between-mode + h²`). But F2 is the variance of the difference when **both** true rates are
unknown. In the BP message (★), the source's uncertainty is **already carried by `V_c^s`** — the message
variance is `V_c^s + σ²_T`, and `σ²_T` is the *conditional* transport `Var(u_c^d | u_c^s)`, i.e. the spread of
the destination given the source's rate. Conditioning on `u_c^s` removes `var_proj(s)`:

```
    Var(u_c^d,true | u_c^s,true)  ≈  var_proj(d)  +  gap²   =   F1 .
```

**So `var_proj(s)` belongs to `V_c^s`, not to `σ²_T`. Putting it in both (F2) double-counts the source's own
uncertainty.** F1 is the conditional transport; the source term is `V_c^s`; the total message variance is
`V_c^s + F1 = v_logfc + 1/M_s + var_proj(d) + gap²`. **This is the shipped formula.**

### 5.4 Where each projected quantity lands — the answer to "where does the variance projection go?"

| quantity | role | appears in |
|---|---|---|
| `μ_proj(s), μ_proj(d)` (modes) | the transport **magnitude** (mode gap) | `σ²_T` via `gap²` |
| `var_proj(d)` (dst ambiguity) | the destination's mode-assignment uncertainty | `σ²_T` (F1) |
| `var_proj(s)` (src ambiguity) | the source's rate uncertainty — **already in `V_c^s`** via the solve | `V_c^s`, **not** `σ²_T` |
| `h²` (within-mode floor, part of `var_proj`) | the count-zero-information max-precision cap (`msg precision ≤ 1/h²`) | the floor of `σ²_T` |

So the variance projection **is** used — the destination's in the transport (F1), the source's inside the
source's own posterior `V_c^s`. Nothing is wasted, and nothing is double-counted.

---

## 6. The count-zero-information check

`σ²_T = var_proj(d) + gap²` is **count-free** — it is mixture geometry, no `M`. The only count entries in the
message variance `V_c^s + σ²_T` are (i) the strand Fisher information inside `v_logfc` (statistical power of the
source's *own* strand signal) and (ii) `1/M_s` (the source's Poisson sampling of its density scale). Both are
statistical power, neither a composition vote. A high count therefore buys no false confidence in a message:
`σ²_T` floors the precision at `1/h²` no matter how large `M` is. This is the domain invariant, satisfied by
construction.

---

## 7. The one genuinely open question for the reviewer

`σ²_T` (5.3) uses `μ_proj(d)` and `var_proj(d)`, which are functions of the **destination's** observed density
`O_d`. Textbook sum-product BP forbids the message `m_{s→d}` from depending on the destination's own data (it
would enter the destination's posterior twice — once through the destination's local factor, once through the
message). **Is our use of `O_d` in the transport variance legitimate?** We see two defensible framings and would
like the reviewer to confirm one, reject both, or supply the correct treatment.

**Framing A — a data-dependent (fixed) edge potential.** The observed total density `M/E` is a **fixed
covariate** of the edge (like the effective length), not the evolving belief. The transport variance is a
physical property of the edge — *are the two sides in different capture regimes?* — and the observed densities
reveal it. So `σ²_T(edge) = f(O_s, O_d)` is a fixed edge-potential parameter computed once, belief-free. The
destination's density then modulates *how much the source's evidence counts* (a reliability weight), which is
not the same as adding evidence about `u_c^d`. Is a covariate-dependent edge potential of this form a sound,
recognised construction, and does it preserve the guarantees we want (no double-counting of `O_d` as evidence)?

**Framing B — expectation propagation against a switching transport.** The true edge potential is a **mixture**:
"same regime" (`u_d ≈ u_s`, tight) with some probability, or "crossing" (`u_d ≈ a different mode`, shifted) —
because a probe boundary genuinely switches the regime. Under this switching transport the exact message is
**bimodal**; projecting it to a Gaussian by **expectation propagation** (EP) makes its moments depend on the
cavity, which includes `O_d`. Our `μ_proj(d)/var_proj(d)`-dependence would then be the EP moment-match of the
switching-transport message. (Our solver already does EP-style moment-matching for the fold.) Is F1 the correct
EP moment-match of a two-component switching transport, or an approximation to it — and if the latter, what is
the exact form?

**Sub-questions.**
1. Is the decomposition `message variance = V_c^s + σ²_T` with `V_c^s = v_logfc + 1/M_s` and
   `σ²_T = var_proj(d) + gap²` (F1) correct, given the double-count argument in §5.3? Specifically, is excluding
   `var_proj(s)` (F2 → F1) the right call, or is there a regime where the symmetric F2 is required?
2. Is the identity-mean transport (`message mean = μ_proj(s)`) appropriate, or should the mean carry a
   regime-shift (a non-zero conditional-mean term) that we are currently folding into the *variance* as `gap²`?
   (I.e. should a crossing down-weight the message via variance, as we do, or *relocate* its mean?)
3. The floor `h²` doubles as the max-precision cap. Is tying the regulariser (bandwidth) to the message-
   precision ceiling principled, or should they be separated?

---

## 8. Empirical validation (for context; the questions above are theoretical)

On synthetic data with a known ground-truth split (an oracle):

* **F2 second-moment identity** — `var_proj(s)+var_proj(d)+gap²` matches Monte-Carlo `E[(u_d−u_s)²]` to <0.3%
  (`scripts/debug/` numeric check). Confirms §5.3's algebra.
* **F1 reproduces the true transport variance** — the belief-free `σ²_T = var_proj(d)+gap²` matches the oracle's
  *raw* adjacent gDNA log-density disagreement variance under capture to ~15%, across strata
  (same-mode ≈ 0.3, crossing ≈ 10–25), and is bandwidth-robust in the crossing (only the floor scales with `h`).
* **F0/F1/F2 differ by <5%** on the suite (the `gap²` dominates), so the choice is a matter of *correctness*, not
  fit — which is why we want the derivation adjudicated rather than chosen by benchmark.
* **Wired into the real solver**, F1 helps every regime that has a fallback (off-capture; stranded capture,
  mwae 0.30→0.04) and regresses only the *unstranded-capture* corner, where the strand is uninformative and the
  weak prior leaves the node leaning on the (now-gagged) crossing messages — an information-starvation effect
  we expect a per-component prior refit to address, not a defect of `σ²_T`.

---

## 9. Summary — what we claim, and what we ask

**We claim (and ask you to check):** the message variance is `V_c^s + σ²_T`; the source term `V_c^s =
v_logfc + 1/M_s` is the source's posterior and *is* the BP source variance; the transport `σ²_T =
var_proj(d) + (μ_proj(d) − μ_proj(s))²` is the conditional transport variance, using the destination's mode
ambiguity and the squared enrichment-mode gap; `var_proj(s)` is deliberately excluded because it is already in
`V_c^s` (F2 would double-count). Everything is count-free except the two legitimate statistical-power channels.

**We ask:** is §7's use of the destination's observation in `σ²_T` justified (Framing A, Framing B, or neither),
and is F1 the right assembly of the pieces — or should a crossing relocate the message *mean* rather than only
inflate its *variance* (§7 sub-question 2)?

---

## 10. Reviewer adjudication (external statistical review, 2026-07-17) — resolved

An external statistical reviewer adjudicated §7. **Verdict: the use of `O_d` is legitimate and F1 is correct.**
Both framings were confirmed as dual representations of the same graphical reality; the sub-questions are all
answered in favour of the shipped form. Summary of the ruling and the two places we tightened it:

**Framing A (data-dependent edge potential) — confirmed, and it is the load-bearing justification.** In a
conditional random field `P(X|Y) ∝ ∏ᵢ φᵢ(xᵢ,Yᵢ) ∏ ψᵢⱼ(xᵢ,xⱼ,Y)`, the observations `Y` are fixed covariates, so
an edge potential `ψ_sd(u_s,u_d; O_s,O_d) = N(u_d−u_s; 0, σ²_T(O_s,O_d))` **may legally depend on the entire
observation set**, including `O_d`. Because `O_d` enters `ψ` *only* as a variance parameter, the sum-product
message `m_{s→d}(u_d)=∫ ψ · q^{∖d}(u_s) du_s` keeps its mean pinned at `m_s` — `O_d` acts purely as a **coupling
gate** (a reliability weight), never shifting the message toward the destination's own value. No state evidence
about `u_d` is transmitted; no double-count. **Precondition check (verified in code):** this argument requires
that `O_d` touch *only* the variance, never the mean. In `bp_solver._scan`, `s2t` enters only the message
precision `pr` ([`bp_solver.py:482`](../../src/rigel/calibration/bp_solver.py#L482)); the message mode `mo` is
computed from the *source's* belief alone ([`:481`](../../src/rigel/calibration/bp_solver.py#L481)). So our
implementation sits exactly on Framing A's precondition.

**Framing B (EP against a switching transport) — confirmed as an exact second-moment projection, with one
correction to the reviewer's write-up.** Model the transport as a mixture: with prob `1−π` "same regime"
(`N(u_d−u_s; 0, σ²_local)`) and with prob `π` a "crossing". Propagating the source cavity and moment-matching the
EP message *under the identity-mean constraint* `m_{s→d}=N(u_d; m_s, Σ)` gives
`Σ = E[(u_d − m_s)²] = Var(u_d) + (E[u_d] − m_s)²`. The reviewer wrote the crossing branch as a reversion to the
**global** prior `P(u_d)` but then read off its moments as `μ_proj(d), var_proj(d)` — those are the **local**
projection at the destination's density, not the global mixture's moments, so the two are not equal as written.
**The fix (which yields F1 exactly):** the crossing branch reverts to the destination's *own local* prior
`u_d ~ N(μ_proj(d), var_proj(d))` — "the source tells you nothing; fall back on what the destination's density
says about its mode" — not to the global mixture. Under this *local-reversion* switching model the EP moment-match
is `Σ = var_proj(d) + (μ_proj(d) − m_s)² ≡ F1` with `m_s ≈ μ_proj(s)`. So F1 **is** the EP projection of a
switching transport — but only of the local-reversion variant, and that variant is itself just Framing A wearing
a generative hat (its crossing branch *is* "use `O_d`"). Framing A is therefore the primitive justification;
Framing B (local-reversion) is a consistency check that reproduces F1.

**Sub-questions — all resolved for the shipped form.**
1. **Exclude `var_proj(s)` (F2 → F1): correct.** `V_c^s` already carries the source's posterior width (which
   incorporates the prior via the solve); adding `var_proj(s)` to `σ²_T` too would count the source prior twice.
   F1 treats the source rate as the conditioning given; `V_c^s` handles the source's own width. No regime needs F2.
2. **Identity-mean transport (gag by variance, do not relocate the mean): correct — and necessary.** Relocating
   the message mean to `μ_proj(d)` on a crossing would actively carry `O_d` back to the destination, creating the
   feedback loop Framing A avoids. Variance-inflation is the sound Gaussian-BP treatment of a crossing.
3. **Tying the floor `h²` to the max-precision cap (`≤ 1/h²`): principled.** `h` is the resolution limit of the
   empirical landscape; a message may never convey rate information finer than the population landscape resolves.

**One belief-free deviation we keep on purpose (noted, not an error).** The reviewer's EP derivation carries the
source's *posterior* mean `m_s`, then substitutes `m_s ≈ μ_proj(s)`. We ship the substitution — `gap² =
(μ_proj(d) − μ_proj(s))²` uses the **belief-free** source projection, not the running posterior — precisely to
keep `σ²_T` belief-free and fixed across passes (the circularity firewall). This is the exact seam where a
*refit* pass would differ: with the solved belief available, `gap²` could use `m_s`. Whether that is worth the
circularity exposure is a design question for the refit-loop examination, not a defect of pass-0.

### 10.1 Implementation-risk determination (which reviewer prescriptions to adopt)

The reviewer flagged three implementation risks. Two of the three target a **discrete regime-stratification +
Poisson-subtraction** design (`[0,1.6,25]` strata) that we **superseded** with the continuous NPMLE projection;
adopting their prescriptions would reintroduce the very knobs we removed. The determination:

| Reviewer risk | Applies to shipped F1? | Determination |
|---|---|---|
| **1. Step-function discontinuities → refit oscillation** | **No** — `project()` is a softmax over the fixed mixture ([`npmle.py:251-257`](../../src/rigel/calibration/npmle.py#L251-L257)); `mu_proj`/`var_proj` are smooth, no strata, no thresholds. | **Do not** add the prescribed sigmoid interpolation (there are no strata to interpolate). **Adopt** the *guard*: a refinement-based continuity regression test (`test_projection_is_continuous_across_the_valley`) so a future reversion to discrete bins fails loudly. |
| **2. Poisson over-subtraction → zero-clip → phantom telegraph** | **No** — `σ²_T` does **zero** Poisson subtraction; `var_proj = between-mode + h²` has a *structural* floor `h² > 0`, so it can never clip to 0. This is a design **strength** of F1 vs. the prototype. | **Do not** add the prescribed Bayesian-shrinkage estimator (two new knobs `ε_v, α` for a failure mode we cannot have). The *related* real item is the opposite-sign Poisson **double-count** (`σ²_T` from raw densities mildly re-counts the `1/M_src` Poisson) — second-order, stays on `CLEANUP_LOG`. |
| **3. Structural-anchor spatial-representation bias** | **Transmutes** — we do not fit `σ²_T` on a separate structural anchor; `var_proj` reads off the NPMLE fit on **all** nodes' *total* density. The live analog: is that total-density landscape's mode structure driven by genuine enrichment geometry, or contaminated by low-count noise / the RNA smear? | **Adopt** the diagnostic (stratify `var_proj` by node coverage depth; confirm the mode structure is structural, not count-driven). The deeper answer **is** the per-component **refit** (deconvolve total → gDNA/RNA NPMLEs so gDNA `var_proj` reflects gDNA enrichment, not the RNA smear) — i.e. Risk 3 is the strongest argument *for* examining the refit. |

**The three validation checkpoints the reviewer signed off on are all adopted** (they are diagnostics/tests, not
knobs): (1) the coverage-stratified `var_proj` diagnostic (Risk 3 analog); (2) the **multi-pass convergence test**
(5 refit iters; assert MWAE decreases or stabilises; watch for post-pass-2 oscillation) — this *is* the
refit-loop examination the next phase runs; (3) the **nascent-factor handoff hook** — the message is already a
product of independent per-component density factors (gDNA + RNA-total + per-strand), so the nascent factory
(`ρ_nascent = ρ_RNA − ρ_mature`, intron-baselined) folds in by updating the RNA-channel *observation* before the
`(λ,θ)` fold, with no change to `σ²_T`.

**Net:** no change to the shipped `σ²_T` assembly (F1 is theoretically vindicated). The actionable items are one
continuity guard test (done) and the refit-phase diagnostics/harness (checkpoints 1–2), which are the next step.
