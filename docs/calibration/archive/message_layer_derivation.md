# The pass-0 message layer — derivation from first principles

**Status:** derivation. **Date:** 2026-07-22. **Branch:** `calib-ambig-init-wip`.
**Supersedes, on every point where they disagree:** `message_arithmetic_reconciliation.md`,
`gdna_factory_and_seam_derivation.md`, `message_layer_implementation_plan.md` §0–§12.
Those documents contain the ideas this one formalizes; §8 lists exactly where they drifted.

**Method.** Every structural claim below is read off the code, not the prior docs
(`bp_solver.py`, `node_geometry.py`, `effective_length.py`, `simplex_logodds.py`, `node_chain.py`).
Every mathematical claim is verified numerically (`scratchpad/derive_check.py`; results inline).
Every empirical claim is measured on the 20 unstranded `ambig_dense_10mb` scenarios
(`scripts/debug/identifiability_census.py`; table in §9).

---

## 0. The one mechanism

> Every node holds a per-component **density** vector `ρ = (ρ₊, ρ₋, ρ_g)` and knows its own per-component
> effective lengths `E = (E_r, E_r, E_g)`. Its observation is the single linear constraint `ρ·E = M`.
> Every message is the same object: a neighbour's `ρ` with a log-variance, attenuated by a structural
> **transmission** `τ_c` and displaced by an **unknown shared log-scale** `δ_e`. A node's solve is that
> constraint, plus its strand likelihood, plus one Gaussian per incident edge on the log-density vector with
> covariance `diag(v) + V_e·11ᵀ`, plus its reference. **Identifiability is the rank of that system.**

There is no vehicle predicate, no peel, no graft, no honest clamp, no emission gate, no `τ`-evidence
bookkeeping. Those are all special cases or symptoms; §7 is the deletion ledger.

The two "message modes" the codebase has been oscillating between are **the two limits of one Gaussian**:

| | `V_e → ∞` (enrichment step unknown) | `V_e → 0` (enrichment step known) |
|---|---|---|
| what survives | the **contrasts** between components | the **absolute** level of each component |
| the old name | composition **SHIFT** | **DENSITY** imputation |
| rank an edge supplies | `\|C_e\| − 1` | `\|C_e\|` |
| a gDNA-only edge | supplies **nothing** | supplies 1 |

That last row is the whole of the capture problem, and §9 measures it.

---

## 1. Ground truth — what the code actually holds

Read from source, so the derivation cannot drift from the implementation.

**Topology** (`node_chain.py`). Per reference the chain is strictly `B R B R … R B`, `k` regions and `k+1`
boundaries, `left`/`right` giving the single neighbour of the other kind, `−1` at a terminal. **Therefore every
edge has exactly one boundary endpoint and exactly one region endpoint.** No per-edge storage is needed: an
edge's structural properties are the boundary endpoint's own arrays.

**Per-node geometry** (`node_geometry.build_node_geometry`, `effective_length.py`).

| node | mass | `E_g` | `E_r` | spliced |
|---|---|---|---|---|
| REGION | contained unspliced, same both faces | `E_f[max(0, L−ℓ+1)]` on the **gDNA** FL | same on the **RNA** FL | `≡ 0` |
| BOUNDARY | per-side crossing, `m_l + m_r` = total | `E_f[min(ℓ,R)]/2` per side | same on the RNA FL | one-sided, motif-stranded, routed to its exon flank |

Two facts follow that the current `_scan` does not use:

* **A boundary's two faces observe the same fragments.** The face split is mass *bookkeeping* (which bases of a
  crossing fragment lie in which flank), not a second sampling. So the two faces cannot sit at different
  capture efficiencies, and `ρ = (m_l+m_r)/(E_l+E_r)` is the single, unbiased, maximally-precise estimate of the
  crossing density. This is exactly `node_global_geometry`. **A node has ONE density per component, not two.**
* **`E_f[max(0,L−ℓ+1)]` and `E_f[min(ℓ,R)]/2` are two views of one physical quantity** — the opportunity for a
  molecule of the local density to be observed at that node. So `ρ` is comparable *between* nodes of different
  kinds, which is what makes a density the message currency at all.

**Solve** (`simplex_logodds.py`). Coordinates `λ = logit(f_g)` and `θ = arcsin(τ)`, `τ = (f₊−f₋)/f_r`.
`ψ = strand + gDNA arm + RNA arm + Gaussian message factors on log f_g / log f₊ / log f₋`. The arms are the
Jeffreys reference `+½·log f` unless a prior is fitted. Messages therefore enter as **quadratic forms in
log-fraction space** — which is precisely the interface §5 needs, unchanged.

**Message construction today** (`bp_solver._scan`). Per face, per component, the source computes a
*destination-frame log-fraction* mode and a precision; `_comb` then averages the two directions' **normalized
log-fractions**. That is the defect the user's requirement (1) names: the sender is doing destination-specific
arithmetic, and the receiver is averaging compositions rather than pooling densities.

---

## 2. The observation constraint and the DOF theorem

A node observes total unspliced mass `M` with integer flux `n`, split `(n₊, n₋)` by genome strand.
Its unknowns are the three component densities. By the definition of the effective length,

```
        ρ₊·E_r  +  ρ₋·E_r  +  ρ_g·E_g   =   M                    (the OBSERVATION constraint — hard)
```

This is exact and it is already encoded: `f_c ≡ ρ_c E_c / M` and `Σ_c f_c = 1`, i.e. **the simplex the solver
already works on IS this constraint**. A structurally-absent strand is not an unknown; it is pinned to 0.

> ### Theorem (DOF). A node's free degrees of freedom = its number of ACTIVE RNA strands.

| node class | active strands | unknowns | − constraint | **free DOF** |
|---|---|---|---|---|
| intergenic region; TSS/TES seam; any no-RNA-crossing boundary | 0 | 1 | 1 | **0** |
| single-strand region or boundary | 1 | 2 | 1 | **1** |
| AMBIG (both strands) | 2 | 3 | 1 | **2** |

The 0-DOF case is not "solvable"; it is **already solved by measurement**: `ρ_g = M/E_g`, exactly, with
counting variance `1/n` and no composition solve. That is the owner's factory, and it is a *theorem*, not a
policy — which is why a 0-DOF node must emit at counting precision and must never be reseeded from a solver
that skipped it.

### 2.1 The strand likelihood supplies rank 1 — never 2

With `p` the genome-plus rate of the unspliced mixture, `κ` the library sense fraction, `τ` the RNA tilt:

```
    p = ½f_g + κf₊ + (1−κ)f₋
      = ½ + ½·f_r·τ·(2κ−1)
```

*One* scalar equation, and it constrains only the **product** `f_r·τ`.

* single-strand: `τ = ±1` is fixed by structure ⇒ the equation determines `f_r` ⇒ **rank 1 solves 1 DOF**.
* AMBIG: `f_r` and `τ` are both free ⇒ **rank 1 against 2 DOF** — one short, always.
* unstranded (`κ = ½`): the equation is `p ≡ ½` ⇒ **rank 0**, identically, for any count.

This is the count-zero-information principle in its sharpest form, and it explains without any further
assumption why AMBIG is the hard class and why unstranded data needs an external source.

---

## 3. SENDING — one object, no gating

> **Requirement (1).** A node sends `{ρ₊, ρ₋, ρ_g}` and `{v₊, v₋, v_g}` and forgets. No gating, no filtering,
> no knowledge of the destination.

### 3.1 What is sent

```
    ρ_c = f_c · M / E_c                    (node-level, both faces pooled — §1)
    v_c = Var(log ρ_c)
```

Three numbers and three variances, computed once per node per sweep direction, sent to both neighbours
unchanged. **All per-face message construction disappears**, and with it `sf`/`df`, `EGs/EGd/ERs/ERd`,
`MSs/MSd`, and the `_dst_mat`/`_src_mat` directional selectors.

### 3.2 The variance — the decomposition is right, the substitution is wrong

Write `n_c = f_c·n` for the component's share of the flux. If the total flux is Poisson and the composition
multinomial, `n_c` is Poisson and

```
    Var(log ρ_c) = Var(log n_c) = 1/n_c = 1/(f_c·n)
```

and this factors **exactly**:

```
    1/(f_c·n)  ≡  (1−f_c)/(n·f_c)  +  1/n  =  Var(log f_c)│multinomial  +  1/n
```

*(verified: `f,n` = (0.5,100) → 0.02 = 0.01+0.01; (0.05,100) → 0.2 = 0.19+0.01; (0.9,40) → 0.0278 both sides.)*

So the code's shape `Var(log f_c) + 1/n` is **structurally correct** — this is a real result, not a
convention. What is wrong is *what is substituted for* `Var(log f_c)`.

### 3.3 The `∞` is an artifact of a bookkeeping proxy, not a fact about the physics

`_scan` sets `v_logfg = ∞` whenever the reference-free evidence `τ_λ = 0`, i.e. on every composition-vacuous
source — which on unstranded data is *most boundaries*. Measured consequence (recorded in
`message_arithmetic_reconciliation.md` §E.4 and reproduced here): **`prec_g = 0.000` on every boundary→exon
edge**, so every mode correction in the whole message thread is multiplied by zero on exactly the data it was
built to fix, and the only surviving channel is `pr += S` — which is why that unsound term is load-bearing.

The infinity is manufactured. `f_c` is a fraction on `[0,1]`; "no evidence" means the node's posterior *is its
reference*, and the reference has a finite log-variance:

```
    f ~ Beta(½,½)   ⇒   Var(log f) = ψ′(½) − ψ′(1) = π²/3 = 3.289868   ⇒   precision 0.3040
    f ~ Uniform     ⇒   Var(log f) = 1                                  ⇒   precision 1
```

*(both verified numerically.)* Neither is infinite. But neither may be **transmitted**, and that is the real
principle the `τ` machinery was groping for:

> **A message must carry the node's EVIDENCE, never its posterior.** Two nodes' compositions are coupled — the
> transport says they are nearly the same parameter — so applying an independent `Beta(½,½)` at every node
> along a vacuous run multiplies one prior `N` times. The reference belongs to the node; the message must be
> the likelihood.

### 3.4 The primitive: divide the reference out

This is standard EP, it needs no new machinery, and it is exact in the Gaussian coordinates the solver already
uses. Let the node's cavity posterior on `log ρ_c` (excluding the message arriving from the destination — the
forward/backward scans already give this) be `N(m_post, 1/p_post)`, and its own reference contribution be
`N(m_ref, 1/p_ref)`. The outgoing message is the quotient:

```
    p_msg = max(p_post − p_ref, 0)
    m_msg = (p_post·m_post − p_ref·m_ref) / p_msg          (undefined and unused when p_msg = 0)
```

Three properties fall out, all of them things the current code special-cases:

1. **A composition-vacuous source emits `p_msg = 0` automatically.** No `τ`, no evidence gate, no `∞`, no
   `_pred_precision` early return. "Send and forget" is satisfiable precisely because the precision self-zeroes.
2. **A 0-DOF node emits at counting precision.** Its posterior is a measurement, not a reference —
   `p_post ≫ p_ref` — so the factory propagates with no unlock and no `struct_lock` special case.
3. **No double counting of the reference along a chain**, which is what the `τ` accumulator, its Jacobian²
   frame conversion, and its self-bootstrap guard were all compensating for.

`p_msg` can hit the clamp at 0 when a message widened the posterior; that is the ordinary EP negative-variance
case and is handled by the clamp, not by a model change.

### 3.5 Additive routes — the sender sums before it speaks (item E, retained)

The one place a sender must do arithmetic is where **one component reaches the destination by two routes**. At
a splice junction the RNA that leaves toward the exon flank is the RNA that continued contiguously *plus* the
RNA that departed through the spliced channel and reappears there. These are two **addends** of one density,
so their variances combine share-weighted (delta method), and precisions must **not** be added:

```
    ρ_r = ρ_ν + ρ_μ ,      w_ν = ρ_ν/ρ_r ,  w_μ = ρ_μ/ρ_r
    Var(log ρ_r) = w_ν²·Var(log ρ_ν) + w_μ²·Var(log ρ_μ)

        Var(log ρ_ν) = 1/n_unspliced + Var(log f_r)     IMPUTED  — counting + the composition solve
        Var(log ρ_μ) = 1/n_spliced                       MEASURED — counting only, NO composition term
```

`prediction_measurement_merge_derivation.md` derives this and it survives intact. Two clarifications it
did not make:

* This is **sender-side**, and it is the *only* per-component arithmetic a sender does. It does not depend on
  the destination, so it does not violate "send and forget".
* Its structural payoff is that `ρ_μ > 0` floors `ρ_r`, so the message can never assert a **confident**
  "no RNA" — the `f_g → 1` phantom pathway. That guarantee is what `pr += S` currently fakes and what
  `σ²_transfer` was silently holding together (the measured +49 % `gdna_none` regression when R2 removed it).

---

## 4. THE EDGE — transmission, and the enrichment step

The sender knows nothing about the edge; the edge is where structure lives.

### 4.1 Transmission `τ_c ∈ [0,1]` — routing, not species

Adopting the owner's reframe (there is one RNA; the question is where it goes):

```
    τ_g     = 1              gDNA is genomic and always continues
    τ_{r,s} = 0              strand s not present on both endpoints (transcript end, strand discontinuity)
    τ_{r,s} = 1              contiguous — mid-region, and an exon↔exon annotation seam
    τ_{r,s} = measured       at a splice junction, and DIRECTIONAL:
                               toward the EXON flank  — continued + departed both arrive  (→ §3.5's sum)
                               toward the INTRON flank — only the continued part arrives
```

`τ` is a pure ratio, so it carries no enrichment factor and does not disturb §4.3.

**The exon→junction direction is the one `τ` is NOT known at the source**, because an exon's solve returns a
lumped `f_r` covering both routes; the simplex has no nascent/mature axis. Under the old framing this forced
the `c_b` "peel". It is not the sender's problem: the **junction measures its own departing density** `ρ_μ`,
so the receiver subtracts what it knows departed. **The peel is a receiver-side operation on the receiver's own
measurement, and `c_b` — which imported a *destination* measurement into a *source's* mode — disappears.**
That also resolves §10.1 of the implementation plan, where a TSS seam had `S_B = 0` and the peel silently
failed to fire: a seam has `τ_r = 0`, so there is nothing to peel in the first place.

### 4.2 `δ_e` — the shared latent log-enrichment step

Capture makes observation efficiency position-dependent. Write the observed density
`ρ̂_c(x) = ε_c(x)·ρ_c(x)`. For an edge `s → d`, define

```
    δ_e = log[ε(d)/ε(s)]           the edge's log-enrichment step
```

`δ_e` is **shared across components on that edge** whenever `ε` is component-agnostic (§6 treats the exception).
Then the transported claim is, per component,

```
    log ρ_c(d)  =  log ρ̂_c(s) + log τ_c(e) + δ_e + noise(v_c)
```

with `δ_e ~ N(m_e, V_e)` a per-edge latent. **This single line contains both message modes.**

### 4.3 Theorem — SHIFT and DENSITY are the `V_e → ∞` and `V_e → 0` limits

Substituting `log ρ_c(d) = log f_c + log M_d − log E_c^d`, each edge contributes a Gaussian on the vector
`u = (log f_c)_{c ∈ C_e}` with mean `a` and covariance `Σ = diag(v) + V_e·11ᵀ`. By Sherman–Morrison the
quadratic form is closed-form and cheap:

```
    Q(u) = (u−a)ᵀ Σ⁻¹ (u−a)
         = Σ_c r_c²/v_c  −  V_e·(Σ_c r_c/v_c)² / (1 + V_e·Σ_c 1/v_c) ,      r_c = u_c − a_c
```

*(verified against a brute-force inverse at |C| = 1, 2, 3 — agreement to 10 d.p.)* Its limits:

| limit | `Q` reduces to | meaning | verified |
|---|---|---|---|
| `\|C\|=1`, any `V` | `r²/(v+V)` | **`σ²_transfer` is exactly this** — the single-component special case | ✔ 0.1000 = 0.09/0.9 |
| `V → 0` | `Σ r_c²/v_c` | independent absolute claims = the **DENSITY** vehicle | ✔ |
| `V → ∞`, `\|C\|=2` | `(r₁−r₂)²/(v₁+v₂)` | pure contrast = the **composition SHIFT** | ✔ |
| `V → ∞`, `\|C\|=1` | **0** | a gDNA-only edge across an unknown cliff carries **nothing** | ✔ 9e−14 |

So `σ²_transfer` is not a legacy proxy to be neutralized (R4) nor the mechanism holding a floor together
(§E.8); it is **the degenerate case of the right object**, and the information the current solver is throwing
away lives entirely in the off-diagonal — the multi-component contrast it never forms.

The last row is the honest statement the "one-sided gDNA floor" and the "factory unlock" were both trying to
approximate: an intergenic→exon seam sends gDNA alone, so **under unknown capture it identifies nothing about
the exon**, and asserting a floor there is what pinned a vacuous exon at `f_g = 0.9616`.

### 4.4 Verification on the closed-form fixture

`tests/calibration/test_bp_solver.py::_mature_exon_chain` — `intron⁺ | exon⁺ | intron⁺`, `ρ_g = 0.5`,
`ρ_m = 1.0`, δ-FLs (gDNA 300, RNA 200), 2000 bp regions, **no capture** (so `δ = 0`, `V = 0`). Computed
against the real `effective_length` functions:

```
    E_g(contained) = 1701.0   E_r(contained) = 1801.0   E_g(crossing, per face) = 150.0
    exon md = 2651.5   TRUE f_g = 0.320762
```

| direction | receiver solves | result | truth |
|---|---|---|---|
| boundary → exon | `ρ_g·E_g + ρ_r·E_r = md` with `ρ_g = 0.5` | **f_g = 0.320762**, and `ρ_r` recovers **1.000000** | 0.320762, 1.0 |
| exon → boundary | same rule, boundary's own `md_B` | **f_g = 1.000000** | 1.000000 |

**Exact in both directions, from one rule.** No vehicle predicate, no peel, no graft. The docs' §11.4 table
(where DENSITY was exact both ways and SHIFT was catastrophically wrong) is reproduced as a *consequence* of
`V = 0`, not as a vehicle choice — and unlike a raw density claim it degrades gracefully as `V` grows.

Note also what the receiver got for free: **`ρ_r = 1.000000`, the exon's mature density, recovered from a
gDNA-only message plus the node's own observation.** The "shift asserts the exon's mature is absent" bug cannot
occur, because an uninformed component is a free parameter the hard constraint solves for — it is never
normalized to zero.

---

## 5. RECEIVING — the lazy solve

### 5.1 The receiver's objective

For node `d`, over its free coordinates `(λ, θ)`:

```
    ψ(λ,θ) =  strand log-likelihood                       (§2.1 — rank 1, or rank 0 when κ=½)
            + reference  ½·log f_g + ½·log(1−f_g)         (unchanged; makes ψ proper)
            + Σ_{e ∋ d}  −½·Q_e(u(λ,θ))                   (§4.3 — ONE rank-1-corrected quadratic per edge)
            + gDNA hyperprior                             (pass-2 only)
```

`Q_e` is a drop-in for the existing per-component Gaussian factors: same `log f_c` arguments, same grid, plus
one rank-1 correction term per edge. It is fully vectorizable over `(λ,θ)`.

**This subsumes the pooling question.** The two flanks are two edges with two *independent* `δ`s, so they are
not summed into one normalizer — which is correct, and is the honest answer to §9.3 of the factory doc
(pooling two flanks does **not** cancel enrichment, because they sit at different coverages). What pooling
buys is rank, not cancellation.

### 5.2 Identifiability — the same divide-out primitive

> **A node may solve iff its evidence precision matrix on its free coordinates is non-singular.**

Evidence precision = (posterior precision) − (reference precision) on `(λ, θ)` — the §3.4 quotient again, now
used as a verdict rather than as a message. This is one primitive with two consumers, and it needs no
threshold beyond "non-singular", because the reference precision is a *derived* scale:

```
    reference precision on log f  =  1/Var_{Beta(½,½)}(log f)  =  1/3.2899  =  0.3040
```

Any message precision above `0.304` is **stronger than the node's own prior** and is therefore driving the
answer. That is the yardstick §9 uses, and it is derived, not tuned.

This is materially different from the refuted §6B DOF gate, which kept the `f_g = 1` init — an all-gDNA
**lock** — for withheld nodes. A node that fails this test must be marked **precision 0**, not defaulted; the
handoff contract to the hyperprior is that pass-2 resolves marked nodes, and a marked node's `f_g` must never
be read as an estimate (this is exactly the metric trap that produced a false "refutation" once already).

### 5.3 The rank ledger

| supply | rank | condition |
|---|---|---|
| own strand likelihood | 1 | `κ ≠ ½` beyond the derived noise floor; **1 even for AMBIG** (§2.1) |
| incident edge `e` | `\|C_e\| − 1` | `δ_e` unknown (capture on, no measured gradient) |
| incident edge `e` | `\|C_e\|` | `δ_e` known (capture off; or a measured coverage gradient) |
| adjacent spliced channel | 1 | on that strand's RNA axis; own `δ` (§6.3) |
| gDNA hyperprior | 1 | on the gDNA axis — **pass-2** |

`|C_e|` = 1 (gDNA) + the number of RNA strands live on **both** endpoints.

### 5.4 The exception list — which nodes cannot self-solve

Derived, then measured (§9).

1. **RNA-isolated nodes under capture.** A node that admits RNA but has *no* incident edge transporting any
   RNA strand — topologically, **the exon of a single-exon transcript**, flanked by two TSS/TES seams. Both
   edges are gDNA-only ⇒ `|C_e| = 1` ⇒ rank 0 under unknown `δ`, rank 1 under known `δ`.
   **Sharp prediction: this class is identified with capture OFF and unidentified with capture ON, while no
   other class flips.** Measured: corr(f_g, oracle) **0.715 → 0.034**, message precision **0.71 → 0.05**.
   Confirmed, and the solver is already *honest* here (precision far below the 0.304 reference scale).
2. **AMBIG anywhere unstranded, unless a neighbour carries both strands and is itself anchored.** 2 DOF,
   strand supplies 0. Unsolvability propagates, so the criterion is reachability from an anchor.
3. **Probe-depleted junctions.** A centred probe depletes junction-crossing reads (~300× measured, trap #3),
   zeroing `S_B`. That removes the rank-1 spliced measurement and demotes a junction boundary to a seam.
   Measured: junction boundaries correlate 0.303 with capture off and **−0.034** with capture on.
4. **Loci with no intergenic flank** (chromosome termini): no 0-DOF anchor in the connected component.
5. **Copy-number change**, which makes an anchor *wrong* rather than absent (kills `ρ_g ≡ ρ_bg`; detectable as
   intergenic `ρ_g` dispersion).

> **"Anchored" is a PRECISION statement, not a boolean** (owner, 2026-07-22). gDNA flows from the intergenic
> flanks, so every node in a locus is technically anchored — the question is only *how weakly*. §6.5 makes
> that quantitative: an edge's precision is set by how well its `δ_e` is known, and `δ_e` is **measurable on
> most edge types**. So class 1 above is not "unanchored"; it is anchored only through edges whose `δ` must be
> taken from the population, and its measured precision (0.05, vs a 0.304 reference scale) is exactly that.

**Two candidate exceptions removed after owner correction, both wrong:**

* ~~Zero-count interior nodes break the chain.~~ **They must not, and the solver has a relay for it** — a node
  folds its incoming message onto its running belief and passes it on, so a node with no data of its own
  should transmit its neighbour's claim, not its own (vacuous) one. Many regions are *structurally* zero: a
  region shorter than a fragment has `E_f[max(0,L−ℓ+1)] ≈ 0` and can contain no fragment **by construction**,
  while its flanking boundaries carry full crossing counts. Measured on the suite: of **33,960** interior
  region nodes, **11,313 (33 %) carry zero mass**, and **5,013 of those have a live neighbour — so the relay
  dies there today**, because emission is gated on the source's own mass (`emit_g = sm > _EPS`). 3,565 are
  short enough that `E_g < 1`. **This is a live bug, not an identifiability class** — see §11 P1.
* ~~Unannotated transcription inside intergenic.~~ **Out of scope: Rigel models only annotated transcripts**
  (owner). Not a design consideration.

**And one class the derivation did not predict, which the measurement forced (§9): every boundary node.**

---

## 6. ENRICHMENT — the model and its one genuine asymmetry

### 6.1 The model

```
    ε_c(x)  =  π(x) · γ_c(x)
```

`π(x)` — positional probe coverage, shared by all nucleic acid at `x`.
`γ_c(x)` — a component-specific hybridization efficiency at `x`.

The owner's three statements map onto this exactly: probes enrich *all* nucleic acid (`π` dominates);
enrichment follows probe overlap (`π` is a geometric function of the probe layout); a probe crossing a splice
junction enriches RNA over DNA (`γ` becomes **position-dependent**, at junctions specifically).

### 6.2 What each vehicle is invariant to

For an edge with the transported claim of §4.2:

| | invariant to `π` | invariant to component-constant `γ` | invariant to position-dependent `γ` |
|---|---|---|---|
| contrast / SHIFT limit | ✅ cancels identically | ✅ same `γ_c` at both ends | ❌ |
| absolute / DENSITY limit | ❌ error is exactly `log[π(d)/π(s)]` | ✅ cancels against a same-`γ` reference | ❌ |

`bp_solver`'s in-code claim that "capture is nucleic-acid-agnostic" is **false but almost harmless**: the toy
simulator itself carries `gdna_split_penalty = 0.2`, yet a *constant* `γ_g` cancels in the contrast and cancels
against the intergenic reference in the factory. The comment should be corrected; the model should not.

### 6.3 The junction-spanning probe — DEFERRED, with one consequence worth stating now

*(Owner, 2026-07-22: defer modelling this.)* Deferred — but it is the reason a splice-junction boundary can
legitimately carry a **higher total density than its adjacent exon region**: a probe spanning the junction is
complementary to the spliced RNA across its full length and captures those fragments preferentially. So a
boundary denser than its neighbour is **not** evidence of an error, and nothing downstream should treat it as
such. Family C's measured median of −0.32/−0.61 under capture (§6.5) is the aggregate of this against the
opposing geometric effect, and it is modest.

The mechanism, recorded for when it is taken up:

A probe spanning a splice junction is complementary to the **spliced** RNA across its whole length, but to
genomic DNA (and to unspliced nascent RNA) only in half-length blocks separated by the intron. Therefore at a
junction

```
        γ_spliced  ≫  γ_unspliced                      (position-dependent by construction)
```

Two consequences, both of which contradict current documents:

1. **`c_b = log1p(S_B/D_B)` is NOT enrichment-invariant.** It is a ratio of the boundary's spliced density to
   its unspliced crossing density at the same position: `π` cancels, but `γ_spl/γ_uns` does **not**. The claim
   "enrichment-invariant, zero constants" holds only under `γ_spl = γ_uns`, which a junction probe is
   specifically designed to violate. (`c_b` is deleted anyway by §4.1 — but the *reason* recorded for it in
   `exon_boundary_mature_dilution_plan.md` and the roadmap is wrong and should not be reused.)
2. **The spliced measurement must carry its own `δ`, not share the unspliced edge's `δ`.** In §3.5's merge the
   two addends `ρ_ν` and `ρ_μ` reach the destination through *different* efficiency channels. They share a
   node and a route, not an `ε`. This is a refinement item E did not make.

### 6.4 Where `δ` is identified

`δ_e` is identified only where composition is known on **both** endpoints — because in general the observed
total-density gap decomposes as

```
    log ρ̂(d) − log ρ̂(s)  =  δ_e  +  log(composition/content ratio)
```

and the second term is what we are solving for. Using the raw gap as `δ_e` is circular. It is legitimate on a
**pure↔pure edge** (intergenic region ↔ intergenic-exon seam, both 0-DOF by §2), where the second term is
identically zero and the gap **is** `δ_e` exactly. That is the only place in the chain where an enrichment
ratio is data-identified without a solve — and it is RNA-invariant (measured: <6 % across nascent 50→200),
hence learnable from pure-gDNA structure alone, with no solve and no circularity.

### 6.5 `δ_e` IS MEASURABLE — on every edge whose two endpoints share a component set

*(Owner taxonomy, 2026-07-22. This replaces two earlier wrong positions of mine: a "population spread over
adjacent pure-gDNA pairs" estimator, which has no data; and its withdrawal to `V_e ≡ ∞`, which was
over-conservative. The truth is in between and it is measurable.)*

**The circularity breaks on one observation.** `δ_e` is confounded by composition (§6.4) — but the confound
**vanishes when the two endpoints carry the same active component set**, because then the composition ratio is
1 by construction. And component sets are given by the **annotation**, not by a solve. So on those edges

```
        δ_e  =  log ρ̂(dst) − log ρ̂(src)          from OBSERVED densities alone — belief-free, no solve
```

is unbiased. This is not circular and it is available before the sweep, exactly like every other pass-0
population input.

| # | edge | src set | dst set | equal? | `δ` from |
|---|---|---|---|---|---|
| **A** | intergenic region ↔ TSS/TES seam boundary | `{g}` | `{g}` | ✅ | unspliced density |
| **B** | intron region ↔ splice-junction boundary | `{g,ν}` | `{g,ν}` | ✅ | **unspliced only — exclude spliced** |
| **C** | exon region ↔ splice-junction boundary | `{g,ν,μ}` | `{g,ν,μ}` | ✅ | **TOTAL — boundary INCLUDES its spliced** |
| **D** | TSS/TES seam boundary ↔ exon region | `{g}` | `{g,ν,μ}` | ❌ | **not measurable** |

**A splice-junction boundary measures TWO ratios** — B against its intron flank and C against its exon flank —
because it knows its own spliced density. That is the structural payoff of the spliced channel, and it is why
a junction boundary is the *best*-connected node in the chain rather than the worst.

Generalizing D: **any boundary↔region pair with different active components is unmeasurable.** In this
topology that is exactly the TSS/TES-seam ↔ exon edge — i.e. exactly the RNA-isolated / single-exon case.

#### Measured (`scripts/debug/enrichment_ratio_census.py`, 20 unstranded scenarios)

```
capture                               family |   edges |   median      sd |              IQR
    OFF             A intgc-reg <-> seam-bnd |    2066 |   -0.090   0.370 | [ -0.27,  0.07]
    OFF    B intron-reg <-> junc-bnd (unspl) |    7546 |   -0.022   0.356 | [ -0.19,  0.13]
    OFF      C exon-reg <-> junc-bnd (TOTAL) |    6212 |    0.008   1.271 | [ -0.22,  0.19]
    OFF    D exon-reg <-> seam-bnd  [UNMEAS] |    1552 |   -1.116   1.916 | [ -2.88, -0.29]
    ON              A intgc-reg <-> seam-bnd |    1685 |    0.093   1.949 | [ -0.33,  0.58]
    ON     B intron-reg <-> junc-bnd (unspl) |    6435 |    0.553   2.457 | [ -0.06,  4.43]
    ON       C exon-reg <-> junc-bnd (TOTAL) |    5429 |   -0.321   1.766 | [ -1.08,  0.15]
    ON     D exon-reg <-> seam-bnd  [UNMEAS] |    1269 |   -1.320   2.921 | [ -6.05, -0.36]
   VSTR             A intgc-reg <-> seam-bnd |      91 |    7.199   2.555 | [  4.03,  8.81]
   VSTR    B intron-reg <-> junc-bnd (unspl) |     969 |    7.859   2.703 | [  5.80,  9.03]
   VSTR      C exon-reg <-> junc-bnd (TOTAL) |    1403 |   -0.614   1.878 | [ -1.60, -0.25]
   VSTR    D exon-reg <-> seam-bnd  [UNMEAS] |     121 |   -4.642   2.037 | [ -5.83, -2.48]
```

**E1 — the taxonomy is empirically correct.** At capture OFF (no enrichment, so a *measurable* family must
centre on 0) the three measurable families give medians **−0.090, −0.022, +0.008** while the unmeasurable one
sits at **−1.116** — a factor of 3.1. Measurable ≈ 0, unmeasurable ≠ 0, exactly as predicted.

**E2 — this independently validates the effective-length frames.** Three different pairings exercising three
different formulas (`region_eff_length` contained, `boundary_side_eff_length` crossing, and the spliced mass
increment) agree with each other to within 0.09 nats at capture off. That is a stronger check of
`effective_length.py` than any unit test currently in the suite, and it also validates **nascent continuity**
(assumption A6) via family B's −0.022.

**E3 — the residual −0.090 on family A is the seam's RNA contamination**, the effect that motivated
`I_struct`'s restriction to intergenic regions. It is ~9 %, i.e. small and bounded — worth knowing, not worth
a mechanism.

**E4 — the spread is the honest prior width, and it is measured, not chosen.** On measurable families
`sd(δ)` is ≈0.36 at capture off and ≈1.8–2.5 under capture. Family D's edges are drawn from the same
population of probe geometries, so that spread **is** `V₀`, the prior width for the one unmeasurable edge
type. No constant is introduced: `V_e` = counting noise where `δ` is measured, `V₀` where it is not — and
`V₀` is a per-library statistic fit once, belief-free, like `κ` and the overdispersions.

This is what makes the owner's §5.4 point precise: a single-exon exon *is* anchored, through a family-D edge,
at precision `1/V₀` — low, but not zero.

---

## 7. The consolidation ledger — what this deletes

| deleted | why it existed | why it goes |
|---|---|---|
| `use_shift` / `use_shift_g` / `edge_vehicle` | choosing between two modes | they are the two limits of one `V_e` (§4.3) |
| `c_b` (`mature_dilution`) and `_var_mat` | restoring mature inside a composition normalizer | the receiver peels with its own measurement (§4.1) |
| the mature **graft** `+SPs/_esp`, `absorb_p/absorb_n` | the mirror of the peel | subsumed by §3.5's route sum + §5.1's constraint |
| the honest clamp `n_eff = ρ·E when ρ ≤ 1/E` | stopping a saturated subtraction becoming a confident zero | `ρ_μ` floors `ρ_r` structurally (§3.5) |
| `pr += S` | the only surviving channel once `v_log = ∞` gagged the rest | replaced by the share-weighted merge; the gag is gone (§3.3) |
| `StrandEvidence` / `tau0_lam` / `tau_lam` / the Jacobian² relay conversion | keeping the reference out of messages | the EP quotient does it exactly (§3.4) |
| `struct_lock`, the `& is_region` restriction, the factory "unlock" | letting a certain node emit | a 0-DOF node's posterior is a measurement; nothing to unlock (§2) |
| `_pred_precision`'s `return 0` branch and the `∞` state | a vacuous source | `p_msg = 0` emerges (§3.4) |
| `emit_g` / `emit_p` / `emit_n` | emission gating | the sender does not gate; `τ_c` is the edge's (§4.1) |
| per-face message geometry (`sf`/`df`, `EGs/EGd`, `_dst_mat`/`_src_mat`) | facing-side construction | one density per node (§1, §3.1) |
| `σ²_transfer`'s separate identity | cliff damping | it is `Q_e` at `\|C\|=1` (§4.3) — kept, generalized, not neutralized |
| `_comb` averaging normalized log-fractions | combining two directions | edges are separate `Q_e` terms in one `ψ` (§5.1) |
| `comp_fl = 1/md` (the "one-fragment floor") | a resolution floor on the mode | **a fraction is a simplex coordinate, not a quotient** (§5.1) — and `md` is a fractional mass, so the floor exceeds 1 on 27–67 % of edges (§9.1 M5) |

Two things are **kept and promoted**: the Jeffreys reference (it is the identifiability yardstick, §5.2), and
`σ²_transfer`'s projection machinery (it becomes the estimator for `V_e`).

---

## 8. Where the current documents drifted

Recorded so they are not re-used as premises.

| document | claim | status |
|---|---|---|
| `message_arithmetic_reconciliation.md` A.1 | "the density mode should be retired, not tuned" | **wrong** — it is the `V→0` limit and is *exact* where `δ` is known (§4.4) |
| same, B.1 / roadmap R4 | `σ²_transfer` is an obsolete proxy to neutralize | **wrong** — it is the `\|C\|=1` case of the correct object (§4.3) |
| same, §E.4 | a vacuous source needs a chosen variance: population prior / uniform 1/12 / Popoviciu ¼ | **question dissolves** — the EP quotient gives precision 0 with no choice (§3.4). None of the three options is needed. |
| same, §E.8 / factory §7 | "the missing piece is a sound statement of the FLOOR" | **wrong framing** — no floor is needed; a gDNA-only edge under unknown `δ` supplies rank 0, and that is the whole content (§4.3) |
| `gdna_factory_and_seam_derivation.md` §3 | "the SHIFT is invalid into an exon, period" | **too strong** — it is invalid as a *normalizer*, but the contrast it encodes is exactly right; the defect was normalizing, not shifting (§4.4) |
| `message_layer_implementation_plan.md` §11.2 | "SHIFT iff `src.active == dst.active`" | **superseded** by its own §12 reframe *and* by §4.3: the vehicle is not a set question |
| same, §11.5 | Step 2 must ship together with Step 5 or it regresses under capture | **restated correctly**: the coupling is `V_e`, and the fix is to carry `V_e` honestly, not to sequence two changes together |
| roadmap §8.2 / `exon_boundary_mature_dilution_plan.md` | `c_b` is "enrichment-invariant, zero constants" | **wrong** — invariant to `π`, not to `γ_spl/γ_uns` at a junction probe (§6.3) |
| `bp_solver.py` `_scan` docstring | "capture is nucleic-acid-agnostic" | **false**, harmless for a constant `γ`, load-bearing at junctions (§6.2–6.3) |
| `bp_solver.py` §9 header comment | "do not re-attempt a composition-conserving mode on exons" | **the rejection was of an unnormalized shift**; the contrast form is untested and is the derived object |
| `bp_solver.py` mode-header comment | `comp_fl = 1/M_dst` is "the one-fragment resolution floor … NOT an arbitrary epsilon" | **the reasoning is right and the quantity is wrong** — `md` is a fractional MASS, not a count, so the floor exceeds 1 on 27–67 % of edges and asserts impossible fractions up to 1e9 (§9.1 M5) |
| `message_arithmetic_reconciliation.md` §E.5 defect 1 | exon→boundary emits `exp(mode_g) = 1.107` | **real, and understated ~2–4×** — real-suite median 2.202 / 1.414 / 4.368, on **100 %** of live edges, and the cause is the `_den = Mg` collapse, not the peel (§9.1 M6) |
| `effective_length.spliced_side_eff_length` docstring | half-triangle is "2–199× low" on continuing faces | **magnitude unverified**. For a δ-FL the exact ratio is `ℓ/R` for `R ≤ ℓ` and 1 above; measured 8.0× at `R=25, ℓ=200`, 2.0× at `R=100`, **1.000 for `R ≥ ℓ`**. 199× needs `R ≈ 1 bp`. Real, bounded, concentrated on very short flanks — and it becomes load-bearing the moment §3.5 un-gags the spliced channel. |

---

## 9. Measured — the identifiability census

`scripts/debug/identifiability_census.py`, 20 unstranded (`ss=0.50`) `ambig_dense_10mb` scenarios, **pass-0
only** (`calib_refit_iters=0`, no hyperprior). Correlation is computed **within** each scenario (mass-weighted)
then averaged, because a cross-scenario pooled correlation is confounded by the library gDNA level — every
node tracks that, including an uninformed one. `r(pr,err)` = corr(message precision, |error|): it must be
**negative** if the solver is honest.

```
capture                            class |  nodes  mass% |    corr |  |err| |   prec | r(pr,err)
------------------------------------------------------------------------------------------------
    OFF        0DOF pure-gDNA (bnd-seam) |   2066   0.2% |     n/a |  0.000 |  10.15 |     n/a
    OFF             0DOF pure-gDNA (reg) |   1296  30.0% |     n/a |  0.000 |   1.13 |     n/a
    OFF 1DOF single-strand (bnd-junction)|   5703   1.0% |   0.303 |  0.443 |   0.90 |    +0.093
    OFF    1DOF single-strand (bnd-seam) |   1587   0.8% |  -0.016 |  0.595 |   0.71 |    +0.223
    OFF         1DOF single-strand (reg) |   6896  57.5% |   0.581 |  0.236 |   4.35 |    -0.708
    OFF        2DOF AMBIG (bnd-junction) |   1096   0.5% |   0.177 |  0.562 |   0.77 |    +0.319
    OFF                 2DOF AMBIG (reg) |   1743   9.2% |   0.461 |  0.325 |   2.50 |    -0.691
    OFF               RNA-ISOLATED (reg) |    193   0.8% |   0.715 |  0.523 |   0.71 |    -0.616
    ON  1DOF single-strand (bnd-junction)|   4171   3.6% |  -0.034 |  0.418 |   0.62 |    +0.263
    ON     1DOF single-strand (bnd-seam) |   2216   1.6% |   0.064 |  0.493 |   0.74 |    +0.092
    ON          1DOF single-strand (reg) |   6454  69.0% |   0.455 |  0.227 |   2.56 |    -0.596
    ON         2DOF AMBIG (bnd-junction) |    783   1.7% |   0.110 |  0.424 |   0.37 |    +0.247
    ON                  2DOF AMBIG (reg) |   1666  18.7% |   0.543 |  0.238 |   2.77 |    -0.608
    ON                RNA-ISOLATED (reg) |    172   2.1% |   0.034 |  0.320 |   0.05 |    -0.060
   VSTR 1DOF single-strand (bnd-junction)|   1380   4.8% |  -0.024 |  0.521 |   0.22 |    +0.604
   VSTR         1DOF single-strand (reg) |   1993  71.4% |   0.240 |  0.278 |   2.75 |    -0.652
   VSTR        2DOF AMBIG (bnd-junction) |    317   2.2% |  -0.043 |  0.517 |   0.13 |    +0.665
   VSTR                 2DOF AMBIG (reg) |    607  18.0% |   0.368 |  0.247 |   3.65 |    -0.739
   VSTR               RNA-ISOLATED (reg) |     50   1.5% |  -0.235 |  0.379 |   0.00 |     n/a
```

**M1 — the sharp prediction holds.** RNA-ISOLATED corr **0.715 (OFF) → 0.034 (ON) → −0.235 (VSTR)**. No other
class flips like this: single-strand regions go 0.581→0.455, AMBIG regions 0.461→0.543. The collapse tracks
exactly what §4.3 says it must: `|C_e| = 1` supplies rank 1 when `δ` is known and rank 0 when it is not.
**This class is already honest** — its message precision collapses in step (0.71 → 0.05 → 0.00), far below the
derived 0.304 reference scale. It is the hyperprior's job, not a bug.

**M2 — `r(pr,err)` inverts sign between regions and boundaries.** Every region class is **−0.60 to −0.74**
(confident ⇒ correct: honest). Every boundary class is **+0.09 to +0.67** (confident ⇒ *wrong*). The inversion
is systematic across all three capture states and both DOF classes, so it is not a mass or expression confound
— that would push both kinds the same way.

**M3 — boundaries are the large exception class, and they are confidently wrong.** Boundary nodes carry
**7.0 % of mass at capture-on and 8.9 % at very-strong**, correlate ≈ 0 with the oracle, at `|err| ≈ 0.42–0.52`,
while emitting precision **0.62–0.90 — two to three times the 0.304 reference scale**, i.e. their messages
dominate the receiving node's own prior. By the §5.2 yardstick this is a **ship-criterion violation**
("never confident-wrong"), and it is far larger than the single-exon case that motivated the search.

**M4 — the spliced channel is what capture destroys.** Junction boundaries are the only boundary class with
any signal at capture-off (0.303 vs −0.016 for seams) and lose all of it at capture-on (−0.034). Consistent
with both §6.3 (`γ_spl ≠ γ_uns` at a junction probe) and trap #3 (a centred probe depletes junction reads
~300×). These two have opposite signs and are **not yet separated** — see §10.

### 9.1 The root cause of M2/M3 — measured, and it is a mass-vs-count frame error

`scripts/debug/exon_to_boundary_probe.py` reads the per-edge `_edge_modes` records the inert `_capture` hook
already emits, and reports what each message actually **asserts**. Live messages only (`prec_g > 0`; `amg` is
zero-initialised, so a non-emitting edge reads back as `exp(0) = 1` and must be excluded — the first run of
this probe was contaminated exactly that way).

```
capture      edge family |   edges |  median   >1.0   >0.95 | med prec  pr>ref |  med md   md<1
-----------------------------------------------------------------------------------------------
    OFF      bnd -> exon |    7440 |   1.000  27.0%   59.0% |    0.259   45.5% | 101.500  27.0%
    OFF       bnd -> reg |    5446 |   0.917   0.3%   33.7% |    0.386   54.6% |1854.000   0.3%
    OFF      exon -> bnd |    8126 |   2.202  71.2%  100.0% |    0.170   40.8% |  14.570  17.7%
    OFF   nonexon -> bnd |    7327 |   0.777  29.1%   31.1% |    0.682   69.3% |  14.535  12.9%
    ON       bnd -> exon |    5103 |   1.000  28.8%   60.0% |    0.055   25.5% |  23.000  28.8%
    ON        bnd -> reg |    3727 |   0.881   1.7%   24.2% |    0.049   32.0% | 151.000   1.7%
    ON       exon -> bnd |    6206 |   1.414  72.0%  100.0% |    0.104   33.4% |   2.915  30.4%
    ON    nonexon -> bnd |    5881 |   0.768  39.8%   41.0% |    0.200   39.7% |   2.145  33.4%
   VSTR      bnd -> exon |     422 |   0.896  40.0%   48.3% |    0.014    3.6% | 122.500  40.0%
   VSTR       bnd -> reg |     655 |   1.000  48.4%   66.3% |    0.011    0.0% |   1.000  48.4%
   VSTR      exon -> bnd |    1444 |   4.368  85.2%  100.0% |    0.046    9.0% |  12.933  39.5%
   VSTR   nonexon -> bnd |     943 |   1.0e9  67.4%   67.6% |    0.174   11.8% |   0.000  67.4%
```

**M5 — messages routinely assert gDNA "fractions" ABOVE 1, and the cause is `comp_fl = 1/md`.**
Compare the `>1.0` and `md<1` columns: **27.0/27.0, 0.3/0.3, 28.8/28.8, 1.7/1.7, 40.0/40.0, 48.4/48.4,
67.4/67.4** — an exact one-to-one identity on every boundary→region family. A message asserts an impossible
fraction **iff the destination face carries less than one fragment of mass**.

The mechanism is in two lines of `_scan`:

```python
    comp_fl = 1.0 / md                                   # "the one-fragment resolution floor"
    _mode_shift(mass_c, den, comp_fl) = log(max(mass_c/den, comp_fl))
```

`comp_fl` is documented as *"a node can never be more certain than its one-fragment opportunity"* — which is
sound **only if `md` is an integer COUNT**. It is a fractional **MASS**: the accumulator splits one fragment
across two boundary faces and across several nodes, so a boundary face's mass is routinely below 1 (median
**14.5** at capture-off but **2.1** at capture-on and **0.000** at very-strong). When `md < 1` the "floor"
exceeds 1 and the message asserts the destination is *more than entirely* gDNA; when the face is empty,
`md = _EPS = 1e−9` and it asserts **1e9**. The receiving factor `−½·p·(log f_g − m)²` is deliberately
un-clipped ("verify-don't-clip"), so this is a **maximal saturating pull to `f_g = 1`** — carried, on 40–69 %
of edges, at precision above the 0.304 reference scale.

This is the same MASS-vs-COUNT confusion R1 fixed in `_pred_precision`, left unfixed one expression away. The
integer flux it needs — `geometry.n_unspl_left/right` — is **already plumbed**.

**M6 — the `exon → boundary` edge additionally asserts `f_g = 1` by construction, on 100 % of live edges.**
Independent of `md`. For that edge the code sets `rho_pos = SPs/_esp − absorb_p`, and a REGION has
`spliced_* ≡ 0`, so `rho_pos ≤ 0` ⇒ `Mp = Mn = 0` ⇒ the shift normalizer `_den = Mg + Mp + Mn` collapses to
`Mg` ⇒ `mode_g = log(Mg/Mg) = 0` ⇒ `exp = 1.000`, and `+c_b` then pushes it above. Every splice-junction
boundary receives a "you are pure gDNA" assertion from its exon flank. `message_arithmetic_reconciliation.md`
§E.5 saw this on the toy at 1.107; on the real suite the **median is 2.202 / 1.414 / 4.368**.

M5 and M6 are present at capture **off** as well as on, so they are not an enrichment problem. Together they
account for M2 and M3 without any appeal to `δ`: boundaries are confidently wrong because they are being told,
loudly and impossibly, that they are pure gDNA.

**M7 — capping the asserted fraction at 1 is INERT. Measured, and it matters.** The obvious one-line,
zero-constant patch (`mode ← min(mode, log 1)`, since a fraction cannot exceed 1) was implemented behind a
temporary flag and A/B'd on both instruments:

| | census corr, `bnd-junction` OFF / ON | `r(pr,err)` ON | `gdna_none` |
|---|---|---|---|
| baseline | 0.303 / −0.034 | +0.263 | 3,821,731 |
| fraction capped at 1 | 0.303 / −0.035 | +0.262 | 3,833,935 (**+0.3 %**) |

Every class moves by ≤0.001 and the phantom guard gets marginally **worse**. The reason is instructive: the
receiving factor is `−½·p·(log f_g − m)²` and `log f_g ≤ 0` everywhere on the grid, so a mode at `m = +0.79`
and a mode at `m = 0` are maximized at the *same* grid edge. **The saturating pull to `f_g = 1` is caused by
`m ≈ 0` — "you are pure gDNA" — not by `m > 1`.**

So **M6 is the load-bearing defect and M5 is a symptom of the same frame confusion**, not a second lever. The
flag was removed; the finding is recorded in `_mode_shift`'s docstring. This is exactly the kind of patch that
would have been landed on plausibility and would have consumed a session — it was disproved in two runs.

> **All of this is UNREPRESENTABLE in the derived design.** There is no per-face `md` (§1: one pooled density
> per node), and a fraction is a **coordinate on the simplex** produced by the hard constraint (§5.1), never a
> quotient that can leave `[0,1]` or collapse to `Mg/Mg`. This is the strongest argument for the
> reformulation: it does not fix these defects, it removes the space in which they can be written — and M7
> shows that patching *within* that space does not work.

---

## 10. What is derived, what is measured, what is still open

**Derived and verified:** the DOF theorem (§2); strand rank 1 (§2.1); the `Var(log ρ) = Var(log f) + 1/n`
identity (§3.2, exact); the finiteness of the no-evidence variance (§3.3); the EP-quotient message (§3.4); the
shift/density unification and its four limits (§4.3, Sherman–Morrison verified); exactness on the closed-form
fixture in both directions (§4.4); the invariance table (§6.2).

**Measured:** the census (§9), M1–M4; the root cause of M2/M3 (§9.1), M5–M6.

**Open — needs derivation or measurement before implementation:**

1. ~~Why boundaries are confidently wrong.~~ **CLOSED by §9.1** — two concrete, capture-independent defects
   (`comp_fl = 1/md` on a fractional mass, and the `_den = Mg` collapse on exon→boundary), both confirmed by
   direct observation rather than ablation. The candidate `τ_r` explanation offered when this item was written
   is *not* the cause; the cause is a frame error and a degenerate normalizer. Both are removed by
   construction in the derived design, so **no separate patch should be written for them**.
2. ~~`V_e`'s estimator.~~ **WITHDRAWN — the estimator has no data (§6.5).** `V_e ≡ ∞` is fixed by fiat, is
   zero-parameter, and needs nothing measured. This removes the only blocking dependency in the plan.
3. **Should `exon → junction boundary` stay silent, or recover `τ_r` from the destination's own split?**
   §6.5 makes it silent (`τ_r` unknown at the source). But the *destination* can measure it:
   `τ_r = D_B/(D_B + S_B)` is a ratio at one node, so positional coverage `π` cancels — this is `c_b`
   re-derived as a transmission coefficient, applied by the receiver. It restores rank 1 on that edge, at the
   cost of the `γ_spl/γ_uns` caveat (§6.3) and of trap #3 (a depleted junction reads `S_B ≈ 0` ⇒ `τ_r = 1` ⇒
   over-transports RNA ⇒ *under*-calls gDNA). **Decide with data, after P0 measures the silent version.**
4. **`γ_spl/γ_uns` at junction probes (§6.3).** Sign and magnitude unmeasured; separable from trap #3's `π`
   depletion by a probe-layout sweep (`toy_inject._write_probes` already supports `full`/`centre`/`junction`).
   Only needed if item 3 goes the `τ_r`-recovery way.
5. **Does the EP quotient actually retire `τ` without regression?** §3.4 says the `τ` machinery is
   compensating for a double-counted reference. That is an argument, not a measurement; it must be A/B'd on the
   `gdna_none` guard **as a delta**, on stranded *and* unstranded arms. Note the contrast form is already
   partly self-protecting: a vacuous source's contrast variance is `2·π²/3 = 6.6` ⇒ precision **0.15**, i.e.
   *half* the destination's own reference precision (0.304), so one such message cannot dominate. A long
   vacuous chain still can, which is what the quotient fixes.
6. **The `spliced_side_eff_length` frame** (§8 last row) — real, bounded, and load-bearing once §3.5 lands.

**Not open:** the "vacuous-source variance" choice (§E.4's three options) and the "sound floor" (§E.8) are both
**dissolved**, not deferred. No constant is introduced anywhere in this document.

---

## 11. THE PASS-0 COMPLETION PLAN — implementation

Nothing is blocked: `δ_e` is measurable on families A/B/C and its population spread gives `V₀` for family D
(§6.5). The receive rule is **one expression**, now with a *finite* `V_e` (§4.3's general Sherman–Morrison
form, not the `V→∞` special case):

```python
def _edge_logfactor(r, v, V):
    """psi += -0.5 * Q.  r_c = log f_c(lam,theta) - a_c ;  a_c = log(rho_c^src * tau_c * E_c^dst).
    V = the edge's log-enrichment-step variance: counting noise where delta_e is MEASURED (families A/B/C),
    the measured population spread V_0 where it is not (family D).

        Q = sum_c r_c^2/v_c  -  V*(sum_c r_c/v_c)^2 / (1 + V*sum_c 1/v_c)

    All degenerate cases automatic (verified): |C|=1 -> r^2/(v+V); any v_c=inf -> that component drops out;
    V->0 -> independent absolute claims; V->inf, |C|=2 -> (r1-r2)^2/(v1+v2).  Note a_c carries no M_dst.
    """
```

Because `a_c` carries no `M_dst`, **the receiver never divides by its own mass** — `md`, `comp_fl` and the M5
fractional-mass defect are absent by construction, as is M6's `_den` collapse (there is no normalizer to
collapse).

### The two surgical fixes first — both measured, both small

| # | fix | what | measured basis |
|---|---|---|---|
| **P0** | **M6 — stop zeroing the exon's RNA on `exon → boundary`.** The branch sets `rho_pos = SPs/_esp − absorb_p`, and a REGION has `spliced_* ≡ 0`, so this is `≤ 0` ⇒ `Mp=Mn=0` ⇒ `_den` collapses to `Mg` ⇒ the edge asserts `f_g = 1`. Restore the exon's own density `fp_s·sm/_er` (same form the intron branch already uses) and keep `+c_b` as the total→unspliced frame conversion. **~2 lines.** | 100 % of live `exon→bnd` edges assert `f_g≥0.95`; median 2.20/1.41/4.37 (§9.1 M6). §6.5 family C shows this edge is component-set-EQUAL and therefore one of the *best* edges, not one to silence. |
| **P1** | **The relay must survive a zero-mass node.** Emission is gated on the source's own mass (`emit_g = sm > _EPS`), so a structurally-empty region stops the forward-backward pass dead. A node with no data must transmit its *cavity belief* (the relayed neighbour claim), not its own vacuous density. | **5,013 of 33,960** interior region nodes (15 %) are zero-mass **with a live neighbour**; 3,565 have `E_g < 1`, i.e. cannot contain a fragment by construction (§5.4). |

**Neither depends on the restructure**, both are testable on their own, and P0 was greenlit. Sequence: P0 →
measure → P1 → measure → decide how much of P2–P8 is urgent.

### P0 — EXECUTED, MEASURED ON THE FULL BATTERY, and **REJECTED** (default OFF, 2026-07-23)

Two lines: `rho_pos/rho_neg = fp_s·sm/_er` on the `exon → boundary` branch (the form the intron branch already
uses). `+c_b` kept. Behind `RIGEL_P0_EXON_RNA`, **default 0**.

> ⚠️ **A first report of P0 called it a −13.1 % win. That was wrong, and the error was methodological:** it
> leaned on the `gdna_none` phantom guard (a **zero-gDNA** instrument) plus a census restricted to the 20
> **unstranded** scenarios. The guard is **one-sided** — in a library with no gDNA, *any* change that lowers
> `f_g` scores better — and P0 lowers `f_g` systematically. The full 32-scenario oracle battery
> (`scripts/debug/pass0_oracle_bench.py`, both arms, per scenario) reverses the verdict.

**Full battery, pass-0 vs oracle, mass-weighted:**

| slice | mwae base → P0 | corr base → P0 | scenarios |
|---|---|---|---|
| **ALL 32** | 0.1396 → **0.1416** ⬆ worse | 0.671 → 0.658 | **16 better / 12 worse / 4 flat** |
| stranded `ss_0.99` (12) | 0.0209 → **0.0316** ⬆ **+51 % rel** | 0.988 → 0.975 | 2 better / 7 worse |
| unstranded `ss_0.50` (20) | 0.2238 → 0.2196 ⬇ | 0.502 → 0.489 | 14 better / 5 worse |
| capture off (14) | 0.1053 → **0.1172** ⬆ worse | 0.808 → 0.758 | 6 better / 7 worse |
| capture on (14) | 0.1420 → 0.1397 ⬇ | 0.668 → 0.674 | 6 better / 5 worse |
| capture verystrong (4) | 0.2964 → **0.2709** ⬇ best gain | 0.225 → 0.273 | 4 better / 0 worse |

**The organizing axis is gDNA LEVEL, not capture.** P0 biases `f_g` **downward**, so it helps wherever the
truth is low and hurts wherever the truth is high:

| gDNA level | direction | worst / best cell |
|---|---|---|
| `none`, `gdna1`, `gdna5` | **uniformly better** | `gdna5_ss_0.50_off` 0.3267 → **0.2098** (−0.117) |
| `gdna100`, `gdna300` | **uniformly worse** | `gdna300_ss_0.50_nrna_none_off` 0.1090 → **0.2039** (+0.095), corr 0.859 → **0.616** |
| stranded, any gDNA | **worse** | `gdna300_ss_0.99_nrna_none_off` 0.0067 → 0.0374 |

**Why:** the exon's RNA now transports, but `+c_b` corrects only the **gDNA** mode — the RNA mode is left
unpeeled — so the exon's within-exon **mature** arrives at the boundary as nascent. Over-called RNA ⇒
under-called gDNA. This is exactly what `test_mature_no_nascent_hallucination_in_introns` asserts, and it is
why that test reverted xpass → xfail under P0.

**Conclusion: the suppression P0 removes was a load-bearing patch for an unfixed mature-accounting bug.**
Removing it alone trades a phantom-gDNA bias for a phantom-RNA bias — visibly better on zero-gDNA instruments,
worse on the benchmark. **P0 must not land alone.** Its prerequisite is the coherent mature peel, which P0b
attempted and which failed on the `spliced_side_eff_length` frame — so that frame is now the blocker, and it
is *independently* confirmed as such by two separate experiments.

**Test states, both diagnostic rather than incidental:**

* `test_geomode_suppresses_exon_unspliced_rna_into_boundary` **FAILS** — it pins exactly the suppression P0
  removes, on the authority of the geo-mean crossing mode's "change 1". **That mode was itself retired** (R3,
  `rho_g_cross` removed as measurably inert). The test now pins an orphaned behaviour. It needs an owner
  decision, not a quiet deletion.
* `test_mature_no_nascent_hallucination_in_introns` reverts **xpass → xfail**, i.e. back to its documented
  state. Its docstring says the exon's ~95 %-mature payload leaks into the flanking introns as nascent
  "because the RNA-total factor does not yet subtract mature". The suppression had been *masking* that leak;
  P0 removes the mask. The leak is `c_b`'s job and `c_b` is applied to the **gDNA mode only** — the E.5
  defect-2 incoherence, now the top open item.

**P0b tried and REJECTED.** Making the peel coherent (subtract `absorb_*` in the RNA numerator so it flows
through the shared normalizer; drop `c_b`) is worse on every axis: guard 3,321,740 → 3,524,635, boundary
`r(pr,err)` −0.205 → +0.033, AMBIG-region corr unchanged. Probable cause: `absorb_*` is a density in the
**spliced** eff-length frame while `rho_pos` is in the **contained RNA** frame, and that spliced frame is
known-wrong wherever the exon continues past the flank (§8, `spliced_side_eff_length`). **Fix that frame
before retrying the numerator peel** — it just moved from "deferred housekeeping" to a blocker.

### P0c — the DIRECTIONAL mature peel (exon↔exon boundaries that are also splice junctions). LANDED

*(Owner, 2026-07-23.)* "Is a splice junction" and "mature crosses contiguously" are **independent** properties
of a boundary, and `mature_dilution` conflated them by summing **both** faces' spliced flux. Owner's example:

```
    TA+ (1000,2000)(10000,11000)      TB+ (1000,2000)(9000,11000)
```

The signature changes at 10000, so there is an exon↔exon boundary there — and it is *also* TA's splice
**acceptor** (2000→10000), while TB's exon runs straight through. Mature both **departs** (TA, into the
spliced channel) and **continues unspliced** (TB) at the same point. The peel must remove only the mature
terminating on the **source's** side, and the accumulator already routes spliced one-sided to its exon flank —
so the source-facing face's flux *is* that quantity. Summing peeled the full junction flux from both directions.

**Measured** (`scripts/debug/exon_exon_junction_census.py`, all 32 scenarios):

| boundary class | count | % of bnd | **% of bnd mass** |
|---|---|---|---|
| seam (no spliced, mature does not cross) | 820 | 48.3 % | 11.8 % |
| classic junction (spliced, mature does not cross) | 758 | 44.6 % | 50.9 % |
| exon↔exon (no spliced, mature crosses) | 43 | 2.5 % | 18.2 % |
| **exon↔exon + splice junction (BOTH)** | **78** | **4.6 %** | **19.1 %** |

Of the 78: **38 carry spliced on one face** (the zero-flux flank was being falsely peeled) and **40 on both**
(a donor one side, an acceptor the other — each direction peeled by the sum). So **37.3 % of boundary mass**
sits where c_b's "the crossing is mature-free" premise is false.

**Fix:** `mature_dilution` / `_var_mat` become per-FACE arrays; `_scan` indexes the face facing the source
(`df`). The summed-eff *frame* is deliberately left untouched, so exactly one thing changes. **Byte-identical
on a classic junction by construction** (spliced lands only on the exon flank, and `c_b` is applied only on an
exon-source edge, so the facing face already holds the whole flux). 313 tests pass; `gdna_none`
3,821,731 → 3,821,084.

**The correction is LARGE at those nodes but INERT on the benchmark — and the reason is diagnostic.** At
exon↔exon+junction boundaries `c_b` is a median **0.4986** nats and the directional error is a median
**0.5367** nats (p90 1.46, max 2.15). Yet all 32 scenarios are flat. Because with the exon's RNA suppressed
(P0 off), `_den` collapses to `Mg`, so the emitted "fraction" is `e^0.50 = 1.65` before the fix and `≈0.96`
after — **both saturate the pull to `f_g = 1`** (the M7 effect). Combined with P0 the battery is 0.1415, i.e.
indistinguishable from P0 alone (0.1416).

> **The three defects are entangled and cannot be validated one at a time.** `c_b`'s correctness is
> unobservable while the normalizer is degenerate; P0 fixes the normalizer but over-transports mature because
> the peel is incoherent; the coherent peel (P0b) fails on the spliced eff-length frame. This is the empirical
> case for the restructure over serial patching — each patch is individually correct and individually
> unmeasurable. **P0c is landed anyway: it is free correctness that removes a demonstrably false peel.**

### A1 — the SPLICED eff-length frame: DECIDED by brute force (2026-07-23)

The one blocking question was decidable without the solver, by enumerating the accumulator's own deposit rule
(`accumulator/00_design.md` §4.3: a slice deposits `(slice_len/ℓ)/n_cross`). Junction at 0, exon flank `[0,R)`,
`a` = the fragment's bases on this side:

```
    a ≤ R : near slice lies inside the flank      → END slice,      n_cross=1 → a/ℓ
    a > R : it overruns into the next region      → INTERIOR slice, n_cross=2 → R/(2ℓ)
    Σ_a  =  R²/(2ℓ) + (ℓ−R)·R/(2ℓ)  =  R/2  =  min(ℓ,R)/2          ← the exon CONTINUES
    Σ_a  =  R²/(2ℓ)                                                 ← the exon TERMINATES at R
```

| case | correct divisor | half-triangle error |
|---|---|---|
| exon **continues** past the flank | `E[min(ℓ,R)]/2` — i.e. **`boundary_side_eff_length`**, already implemented | **ℓ/R** (measured: 2.0× at R=100/ℓ=200, 8.0× at R=25, 12.0× at R=25/ℓ=300) |
| exon **terminates** at the flank's far edge | the half-triangle `E[min(ℓ,R)²/(2ℓ)]` | exact |

Brute force matches each to 4 d.p. (the 0.5 offset at `R ≥ ℓ` is the `a = 1…ℓ−1` discretization). So §8's
"2–199× low" was right in kind — 199× needs `R ≈ 1` — and the **fix needs no new function**: select per face
by whether the exon continues, which is the far boundary's `mrna_active_s`. **Zero constants.**

This is the frame `absorb_*` lives in, hence the direct explanation of P0b's failure.

### A2 — LANDED (2026-07-23). Plus a correction to A1's *predicate*

Implemented in `node_geometry.build_node_geometry`: `eff_spl_left/right` now select per face.

**A1's formula was right; its first predicate was wrong, and the accumulator reference caught it.**
I initially selected on `mrna_active_s` alone ("exon on both sides of the far boundary"). Reading the
reference crossing path (`tests/native/_accumulator_reference.py`) shows a slice is INTERIOR iff **another
slice follows it** — `crosses_right = i < n_slices−1` — which is a property of the *fragment's slice list*,
not of region geometry. For a spliced fragment the next block lands in a different region whether the far
boundary is an exon↔exon seam **or** an exon↔intron splice junction. So an internal exon continues in both
cases, and only a genuine transcript END terminates:

```
    continues_s(far)  =  mrna_active_s(far)   OR   far is a splice junction on strand s
```

The wrong predicate corrected only **66** faces; the right one corrects **~280** — i.e. **27 % of all
spliced-carrying faces**, at a median ratio of **3.1×** and a maximum of **199.6×**, which pins §8's
"2–199× low" exactly. It is self-limiting: the two divisors coincide once `R ≫ FL support`, so only short
flanks move — precisely where the geometry matters.

| gate | result |
|---|---|
| `pytest tests/calibration/ tests/native/` | **313 passed, 2 xpassed** |
| δ-census (the frame regression test) | family medians unchanged (C at capture-off 0.008, sd 1.270 vs 1.271) |
| `gdna_none` guard | 3,821,731 → **3,813,747** (−0.2 %) |
| full 32-scenario battery | mwae 0.1396 → **0.1395**, corr 0.671 → **0.672**, **32/32 flat** |

Neutral-to-marginally-positive everywhere, no regression in any slice — exactly what `effective_length.py`
predicted ("*inert today because the spliced channel is heavily down-weighted… becomes load-bearing once the
spliced-measurement channel is un-gagged*"). **Landed as free correctness.**

**But it did NOT rescue P0b.** Retried with the corrected frame, the coherent numerator peel is still
net-negative on the battery: mwae **0.1443** vs base 0.1396 (and vs P0's 0.1416), 15 better / 12 worse, with
the same signature — every low-gDNA scenario better, every high-gDNA and every stranded one worse. So **the
peel mechanism is not the blocker**: transporting the exon's RNA into the boundary at all is what biases
`f_g` downward, under either peel form. The remaining suspect is the *magnitude* of the peel — `S_B` is a
junction **flux**, and whether it correctly represents the exon's mature **content** is the A3 question.

### A3 — the `c_b` DENOMINATOR: decided exactly, LANDED (2026-07-23)

`c_b = log1p(S_B/D_B)` converts an exon's mature-**inclusive** composition into a boundary's mature-**free**
crossing composition. `S_B` is the boundary's spliced (mature) density; spliced mass lands on **one** face, so
it must be divided by **that face's** opportunity. Production divided by the **sum of both faces'** — counting
an opportunity that can never receive mass. Exact closed-form check (`scripts/debug/cb_denominator_check.py`,
`ρ_g=0.5, ρ_ν=0.3, ρ_μ=1.0` ⇒ `r_true = 1.25`):

| denominator | `r = S_B/D_B` | `c_b = log1p(r)` | |
|---|---|---|---|
| SUMMED-EFF (production) | 0.714 | 0.5390 | **under-peels 33 %** |
| **PER-FACE** + ideal `D_B` | **1.250** | **0.8109** | **exact** |
| PER-FACE + production `D_B` | 1.429 | 0.8873 | over by the `D_B` frame error alone |

**This closes the historical puzzle.** The old comment "*the per-side density sum over-corrects ~2×*" was a
real observation of a **different** bug: `eff_spl` was the half-triangle, which A1/A2 showed is up to 199× low
on a continuing face (1.88× at R=100). Summing the denominators cancelled part of that. **Two compensating
errors** — A2 fixed the frame, so the correct per-face denominator can now be used.

**Residual, deliberately not fixed here:** `D_B` divides a mass containing **nascent** (RNA FL) by the
**gDNA**-FL opportunity, so it is not `ρ_g+ρ_ν` — measured `D_B/(ρ_g+ρ_ν) = 0.875`, i.e. `c_b` still
over-peels ~14 %. Removing it requires the composition, which is what the solve produces — so it belongs in
the restructure, not in a patch.

**Measured on the full battery (all 32, mwae):**

| arm | mwae | vs base 0.1396 |
|---|---|---|
| A3 alone (P0 off) | **0.1397** | flat, 32/32 — safe to land |
| P0 **before** A3 | 0.1416 | +0.0020 |
| **P0 with A3** | **0.1405** | **+0.0009 — the gap more than halves** |
| P0b (numerator peel; A3 does not apply) | 0.1443 | unchanged, as expected |

> **The predicted direction held.** Under-peeling was a real part of P0's regression: too little mature was
> removed, so the exon's mature arrived at the boundary as nascent, RNA was over-called and `f_g` biased
> downward. Fixing the denominator closed **55 %** of P0's gap to baseline. The residual is the `D_B` mixed-FL
> frame error — which cannot be fixed without the composition. That is the cleanest possible statement of why
> the restructure is required rather than optional.

Gates: 313 passed / 2 xpassed; `gdna_none` 3,821,731 → **3,818,168**; battery flat with P0 off.

---

## 12. B1 — THE RELAY THROUGH DATA-FREE NODES. Design (2026-07-23)

**Design only. Not implemented.** Owner directive: a node with no fragments must not stop the
forward–backward pass. Many regions are *structurally* empty — a region shorter than a fragment has
`E_f[max(0, L−ℓ+1)] ≈ 0` and can contain nothing **by construction** — while its flanking boundaries carry
full crossing counts.

### 12.1 There are TWO defects, not one, and their order matters

**D1 — EMISSION.** `_scan` builds the message density from the source's *own* mass:

```python
    rho_g = fg_s * sm / _eg          #  sm = MSs[lsrc], the source's facing mass
    emit_g = sm > _EPS               #  ⇒ a data-free node sends NOTHING
```

The density a data-free node *should* report is its **neighbour's** — gDNA is genomically continuous and
nascent is continuous within a transcript. Its own mass is the wrong numerator.

**D2 — FALSE CERTAINTY.** `bp_solver.py:467,558-570`:

```python
    solvable = (fp | fn) & (statics.mass_unspliced > 0.0)
    locked   = ~solvable
    lam_loc  = np.where(locked, lam_locked, lam_loc)   # lam_locked = logit(f_g_init) = +L  (all-gDNA)
    lvar_loc = np.where(locked, 0.0, lvar_loc)         # variance 0  ⇒  CERTAIN
```

One predicate covers two different things:

* `~(fp|fn)` — **no admissible RNA strand** (intergenic, a no-RNA-crossing seam). Its unspliced mass is
  *necessarily* gDNA. Certainty is **correct**, and this is what the G1-emission fix was for.
* `(fp|fn) & mass==0` — **admits RNA, but has no observation**. `init_beliefs` itself sets
  `var_gdna = inf` ("no information"); the sweep then overrides that to `0` ("certain, pure gDNA").

> **D1 must not be fixed without D2.** The 790 falsely-certain nodes are currently *silent* only because they
> have no mass. Restore emission alone and every one of them starts broadcasting "I am certainly pure gDNA"
> at full confidence — a phantom-gDNA emitter on 23 % of the chain. **D2 is the prerequisite.**

### 12.2 Measured (`scripts/debug/relay_break_census.py`, all 32 scenarios, per scenario)

```
   nodes                                                3,397  (100.0%)
   live                                                 2,328  ( 68.5%)
   dead (mass == 0)                                     1,069  ( 31.5%)
     dead & STRUCTURAL gdna (correctly certain)           279  (  8.2%)
     dead & has RNA strands (FALSELY certain)  ← D2       790  ( 23.3%)
       ... REGION 463   BOUNDARY 327
   dead interior w/ live neighbour (relay BREAKS) ← D1    565  ( 16.6%)
     ... eff_g < 1 (cannot contain a fragment)            185  (  5.4%)
```

**Chain segmentation** — maximal runs of live nodes, split at every dead node; an *anchor* is a live
structural-gDNA intergenic node (the only measured pure-DNA reference):

```
   segments reaching an anchor :  104    holding 64.6% of mass
   segments with NO anchor     :  303    holding 35.4% of mass
   segment size: median 3, mean 5.7, p90 15, max 67
   singleton segments (1 isolated node): 38.8%
```

> **74 % of segments never reach a measured gDNA reference, and they hold 35 % of the mass.** This corrects
> the earlier claim (§5.4) that gDNA "flows from the intergenic flanks so every locus is anchored, weakly" —
> **it does not flow, because the chain is severed**, and 38.8 % of segments are a single node talking to
> nobody. B1 is therefore not a tidy-up; it is what makes the anchoring argument true in the first place.

### 12.3 The correct semantics — three node states, not two

| state | predicate | belief about own composition | emits |
|---|---|---|---|
| **STRUCTURAL gDNA** | `~(fp\|fn)` | `f_g = 1`, **certain** | its own measured `ρ_g` at counting precision |
| **OBSERVED** | `(fp\|fn) & mass>0` | solved | own `ρ_c` **fused with** what it heard |
| **DATA-FREE** | `(fp\|fn) & mass==0` | **unknown** (`var = ∞`, as `init_beliefs` already says) | **what it heard**, transported |

### 12.4 Mechanism — the zero-mass case is a LIMIT, not a branch

The relay must carry the **density**, not only the composition `λ`. Per component `c`, a running
`(m_c, p_c)` = (log-density mode, precision):

```
    own_c     = log(f_c · M_i / E_c^i)        precision  1/(Var(log f_c) + 1/n_i)     → 0 when M_i = 0
    running_c = precision-weighted fusion of (own_c, incoming_c)
    outgoing_c = (running_m , 1/(1/running_p + σ²_transfer))
    emit_c    ⟺  running_p > 0
```

* `M_i = 0` ⇒ own precision `= 0` ⇒ `running = incoming` **exactly**: pass-through, with no special case.
* `M_i` large ⇒ own observation dominates.
* A node with *little* mass blends honestly — the discontinuity at `mass = 0` disappears, which the current
  `sm > _EPS` gate creates artificially.
* The emission gate becomes a consequence rather than a rule.

**This is exactly Stage C's `send()` primitive**, arriving early for the degenerate case — so B1 is a down
payment on the restructure, not a throwaway patch.

**Note a property this buys:** under the current design `σ²_transfer` does **not** accumulate along a relay
(the `τ` mechanism re-derives precision from evidence at each hop). Under the running-density form the
variance accumulates naturally per hop, which is the correct behaviour for a claim travelling several nodes
and is what bounds long-range over-confidence.

### 12.5 Sequencing, risks, gates

| step | change | expected | risk |
|---|---|---|---|
| **B1a** | split the predicate: `struct_gdna = ~(fp\|fn)` keeps `lvar=0`; `no_data` gets the reference `lvar = π²` | **near-inert** — see §12.6 for what actually happened | low — and it is the safety prerequisite for B1b |
| **B1b** | running per-component density + precision-weighted fusion; retire `emit_* = sm > _EPS` | behavioural: 565 broken relays reconnect | medium |

**Gates.** B1a: byte-identical or near — if it moves anything materially, the false certainty was
load-bearing and that must be understood before B1b. B1b: `gdna_none` **as a delta** (the phantom risk is
real — claims now travel further); the full 32-scenario battery per condition, stranded **and** unstranded;
and a **re-run of the segmentation census**, where the fraction of mass reaching an anchor must rise from
64.6 % toward ~100 % (bounded by references with no intergenic node at all).

**Numerical care.** `_fold_lambda` does `var = max(var, _EPS)`; feeding `∞` must be checked — prefer an
explicit "no belief" sentinel over a literal `inf` in the fold.

**Out of scope for B1.** Marking data-free nodes precision-0 for the hyperprior handoff (that is P7).

**Output inertness — VERIFIED, not assumed.** A data-free node's own `f_g` never reaches the result:
the region projection gives `gdna_mass = f_g·mass_u = 0` trivially, and the boundary-side projection
(`chain_boundary_side_deconv`) applies a *boundary's* `f_g` to the *region's* side mass — which could in
principle be nonzero at a dead boundary. Checked directly: **0 of 6,792 region-sides** have a zero-mass
boundary paired with nonzero side mass. So B1 changes only what flows *through* these nodes, never what is
read *out of* them — which is what makes B1a's expected inertness a testable prediction rather than a hope.

### 12.6 B1a — LANDED (2026-07-23), and a design correction the gate caught

Implemented as designed: `struct_gdna = ~(fp|fn)` keeps the `λ=+L, σ²=0` lock; a data-free node keeps its
skipped-solve `λ = 0` and takes the **derived** reference variance

```
    Var(λ) = Var(logit f) under Beta(½,½) = ψ′(½) + ψ′(½) = 2·(π²/2) = π²      (verified: 9.86960440108936)
```

No tuned constant: it is the same `_JEFFREYS_REF` the solver already applies, read in the λ coordinate. And
`E[λ] = ψ(½) − ψ(½) = 0` is exactly what the skipped local solve already returns, so only the variance is
written.

| gate | result |
|---|---|
| `pytest tests/calibration/ tests/native/` | **313 passed, 2 xpassed** |
| `gdna_none` guard | **3,818,168 — bit-identical to pre-B1a** |
| full 32-scenario battery | 30/32 bit-identical; 2 differ by **max \|Δmwae\| = 6.0e−08** |

**The byte-identity prediction was FALSIFIED, for a reason worth keeping.** I argued these nodes "are silent
today, so removing their false certainty cannot change what they emit." Wrong: **`mass_unspliced == 0` does
not mean "no data"** — a *boundary* can have an empty unspliced crossing and still carry **spliced** mass, and
the spliced measurement channel bypasses the unspliced emission gate entirely
(`emit_p = … and (sm > _EPS or SPs[lsrc] > _EPS)`). Measured: **65 such spliced-only boundaries** in
`gdna100_ss_0.99_nrna_none_capture_on`, and **0** in `gdna300_ss_0.50_..._capture_off` — i.e. they exist in
exactly the two scenarios that moved, and nowhere else. The path is `lvar → the fold → μ_λ → the Jacobian² τ
conversion → v_logfp → the emitted RNA prediction precision`.

**This makes B1a *more* correct, not less.** A spliced-only boundary was previously seeded "certainly pure
gDNA" (`λ=+L, σ²=0`) while simultaneously emitting a spliced **RNA** measurement — a self-contradiction. It is
now honestly uninformed on `λ` and still emits its measurement.

**Carried into B1b:** the emission gate cannot simply become "`running_p > 0`". The spliced channel is a
*separate provenance* (§3.5) that is already independent of the unspliced gate and must stay so — a node with
no unspliced mass but real spliced counts is not information-free, it is uninformed *about λ only*.

### 12.7 The intergenic FACTORY and the seam — reconciled with B1b (owner, 2026-07-23)

**First, two variances that this document had been sloppy about.** At an intergenic region:

| quantity | value | variance |
|---|---|---|
| **composition** `f_g` | ≡ 1 — no admissible RNA strand | **0 — genuinely certain** |
| **density** `ρ_g = M/E_g` | a *measurement* | **`1/n` — counting precision** |

`σ²_λ = 0` is the *composition* lock; it is not the message precision. The code already separates them:
`_pred_precision(n, v_log=0, s2t) = n/(n·s2t + 1)` — the **counting** precision, damped by `σ²_transfer`.
So an intergenic region already sources gDNA at count precision, exactly as the owner's model requires.

**Measured** (`scripts/debug/intergenic_factory_probe.py`, all 32 scenarios):

```
capture                    edge family |   edges |  prec_g=0  med prec   >ref |  med src n
    OFF      intergenic-REGION -> seam |    4320 |     0.0%    8.2130  99.8% |    4114.50
    OFF SEAM -> exon  (the factory hop)|    3786 |     0.6%    0.4072  60.6% |      11.02
    ON       intergenic-REGION -> seam |    4276 |     0.0%    7.7402  66.0% |     388.00
    ON  SEAM -> exon  (the factory hop)|    3258 |    14.8%    0.1341  35.5% |       1.74
   VSTR      intergenic-REGION -> seam |     606 |     0.0%    0.1906  18.3% |       6.00
   VSTR SEAM -> exon  (the factory hop)|     159 |   100.0%    0.0000   0.0% |       6.80
```

**The factory's first hop works; its second hop is crushed.** The intergenic region speaks at precision
**8.21** off a count of **4,114**. One hop later the seam speaks at **0.41** (capture off), **0.13** (on) and
**exactly zero on 100 % of edges** (very strong) — off a count of **11 / 1.7 / 6.8**.

**Why: the seam re-derives its message from its OWN mass instead of relaying what it just heard.** That is the
same `rho_c = f_c · sm / E_c` construction as D1 — a node's message is built from `sm`, its own facing mass.
A boundary's crossing count is two to three orders of magnitude smaller than a region's contained count, so
the seam discards a 4,114-count measurement and substitutes an 11-count one.

> **This is D1 at a different mass scale.** Zero mass (§12.1) is the extreme case; low mass — the seam — is
> the common one. Both are the same defect, and B1b's fusion rule fixes both **with no special case**:
>
> ```
>     running_c = precision-weighted fusion( own_c , incoming_c )
> ```
>
> * **intergenic** (`n ≈ 4,114`): its own term dominates ⇒ it **SOURCES**, and overwrites whatever arrived —
>   which is precisely the owner's "sink": a message that reaches an intergenic region dies there.
> * **seam** (`n ≈ 11`): the incoming term dominates ⇒ it **RELAYS** the factory onward at near-full precision.
> * **exon** (large `n`, uncertain composition): the two genuinely blend.
>
> **The owner's source / sink / relay taxonomy is not three rules — it is one fusion, and the precisions
> decide which behaviour each node exhibits.** RNA is the mirror image: a splice junction sources and sinks
> RNA (its spliced channel is the measurement), and no RNA crosses an intergenic node (`free_s = False`), so
> the two anchor types are already structurally separated by the existing gates.

**Companion change: `struct_lock` must become structural.** `struct_lock = locked & is_region`
(`bp_solver.py:288`) excludes boundaries, so a seam is not composition-certain even though no RNA can cross
it — it falls through to the `τ` path and emits `v_logfg = ∞ ⇒ pr = 0`. Composition certainty is a property
of the **component set**, not the node kind. Dropping `& is_region` was tried and reverted (Part E.8) because
it pinned a vacuous exon at `f_g = 0.9616` — but §9.1 **M6** now explains that: the exon-facing message had a
collapsed normalizer (`_den = Mg`) asserting `f_g = 1`, and the unlock merely supplied the confidence to make
it stick. The mode defect is diagnosed; the unlock should be retried **with** the B1b fusion, not alone.

**Sequencing consequence.** B1b is now scoped to cover both: (a) data-free nodes pass through, (b) low-count
nodes relay rather than overwrite. `struct_lock` structural rides with it and is gated the same way.

### 12.8 B1b — IMPLEMENTED (2026-07-23). The first change this session that improves the benchmark.

Behind `RIGEL_B1B` (`1` = density relay; `2` = also `struct_lock` structural), **default off**. Flag-off is
bit-identical (guard 3,818,168, 313 pass).

**What it does.** A running per-node fused gDNA log-density, in the frame-free density currency:

```
    own(i)  = ( log(f_g^loc · M_node/E_g^node) , 1/(Var(log f_g^loc) + 1/n_node) )     node-level, both faces pooled
    run(i)  = precision-weighted fuse( own(i) , transported run(src) )
    message = ( run(src) , 1/(1/p_src + σ²_transfer) )        emit ⟺ p_src > 0
```

Replaces `rho_g = f_g·sm/E_g` (the source's own facing mass) and `emit_g = sm > _EPS`. Because the scan is
genomic, `run(src)` is already the cavity for the message to `i`, so nothing is fed back into its own origin.

**Results — all 32 scenarios, pass-0 vs oracle:**

| slice | mwae base → B1b | corr | scenarios |
|---|---|---|---|
| **ALL 32** | 0.1396 → **0.1361** ⬇ | 0.671 → **0.688** | **16 better / 6 worse / 10 flat** |
| unstranded `ss_0.50` | 0.2238 → **0.2175** | 0.502 → **0.528** | 15 better / 5 worse |
| **stranded `ss_0.99`** | 0.0209 → 0.0214 | 0.988 → 0.989 | **1/1, 10 flat — neutral** |
| capture off | 0.1053 → **0.1000** | 0.808 → **0.831** | 7 / 2 |
| capture on | 0.1420 → 0.1419 | 0.668 → 0.676 | 5 / 4 |
| capture verystrong | 0.2964 → **0.2850** | 0.225 → **0.253** | **4 / 0** |
| `gdna_none` guard | 3,818,168 → **3,704,635 (−3.0 %)** | | |

**The stranded arm is neutral** — this does *not* break what already works, unlike P0 (+51 % relative there).

**The factory gate (§12.7) passes.** Fraction of messages emitted at exactly zero precision:

| edge family | before | after |
|---|---|---|
| `SEAM → exon`, capture off / on / verystrong | 0.6 % / 14.8 % / **100 %** | **0 % / 0 % / 0 %** |
| `intron → junction`, off / on / verystrong | 29.0 % / 37.4 % / 85.4 % | 16.0 % / 13.8 % / 18.7 % |

Nothing is silenced any more. The residual zeros are segments that reach no anchor at all.

**`struct_lock` structural (mode 2) is a wash** — 0.1359 vs 0.1361, one more scenario worse. Exactly as §12.7
predicted: the fusion *subsumes* it, because a seam now relays the intergenic's density rather than
re-deriving a composition from its own 11 counts. **Recommend mode 1** — the smaller change, and it avoids
re-litigating the E.8 revert.

**One red test:** `test_mature_measurement_disagreement_silenced` — `|Δf_g| = 0.0612` vs a `0.05` tolerance
(0.134 depleted vs 0.073 normal). It measures sensitivity to a **probe-depleted junction** (trap #3), a known
one-sided bias that no precision term can represent. B1b makes a depleted junction slightly more
consequential because the relay now carries more weight. Owner decision: re-derive the bound or accept it.

**Hygiene:** a strict `-W error::RuntimeWarning` run caught `0·∞ = nan` in the masked own-precision
expression (`np.where` evaluates both branches) — the same hazard `_pred_precision` avoids by branching.
Fixed by substituting a finite value before the product; both arms are now warning-clean.

### 12.9 THE SINK MODEL — where messages must die (owner, 2026-07-23)

> *Messages should always be emitted. Each node decides for itself what to do with what it receives.*

Emission becomes **unconditional**; all gating moves to the **destination**, per component, and splits into
**two distinct operations** the current code does not distinguish:

* **ACCEPT** — integrate into this node's own belief.
* **RELAY** — include in what this node passes onward.

A sunk component is neither accepted nor relayed. This is **TSS/TES logic**: a boundary catches a message
emitted past the start or end of a transcript and kills it.

| node class | gDNA accept / relay | RNA₊ | RNA₋ | why |
|---|---|---|---|---|
| **intergenic REGION** (with mass) | **SINK / SINK** — sources its own | SINK | SINK | the **locus wall** |
| intergenic REGION (no mass) | relay / relay | SINK | SINK | nothing to defend — must stay transparent |
| **intergenic↔transcript SEAM** | accept / relay | **SINK** | **SINK** | TSS/TES — no RNA continues |
| single-strand node, strand *s* absent | accept / relay | per *s* | per *s* | **SINK** the unsupported strand |
| intron / exon / AMBIG | accept / relay | per `free_s` | per `free_s` | ordinary conductor |

**Two different justifications, and they should not be conflated.** RNA sinks are **structural**: the
component genuinely cannot continue (a transcript ends; a strand is absent) — already encoded by `free_s`, and
already enforced today via `emit_p = fp[lsrc] and fp[i]`. The **gDNA sink at an intergenic region is a
modelling CHOICE**, since gDNA *is* genomically continuous. Its justification is threefold: (a) `ρ_g`
uniformity is only local — CNV and mappability break it at range; (b) an intergenic region measures `ρ_g`
directly (`f_g ≡ 1`, no solve), so it has nothing to learn from inside the locus; (c) it bounds error
propagation to one locus. State it as a choice, not a theorem.

**Letting precision decide instead is NOT safe — measured.** Under B1b without a wall, the fraction of
intergenic-sourced gDNA messages that differ from the region's own measurement is **11.3 % at capture-off and
60.9 % at capture-on**. The magnitude is negligible today (max 0.002 nats) *only* because the intergenic count
(n≈4,114) out-shouts the incoming — and that count falls to **n≈6 under very strong capture** (§12.7), where
it would lose. **The wall was absent, merely out-shouted.** That is precisely the "messages run wild" risk.

**Implemented** (conditioned on having a measurement, so a dead intergenic region stays transparent and does
not re-sever the chain B1b just reconnected). **Behaviourally neutral today** — battery identical to B1b
without it (0.1361 / corr 0.688, 16 better / 6 worse) — i.e. it costs nothing now and removes the failure mode
that unconditional emission would expose.

**What remains before emission can be made unconditional:**

1. **RNA has no relay yet.** B1b relays gDNA only; RNA still uses the per-face construction gated at the
   source by `fp[lsrc] and fp[i]`. When RNA gains a running relay (Stage C), the same ACCEPT/RELAY split must
   be applied — and `fp[lsrc]` becomes unnecessary (a source without strand *s* has precision 0 for it
   automatically), while `fp[i]` must move to the destination as the ACCEPT gate.
2. **ACCEPT ≠ RELAY must become explicit in code.** Today one predicate does both. The intergenic wall is the
   first case where they differ (accept: no; relay: no; but it still *sources*), and the exon↔exon seam will
   be another.
3. **A sunk message must not reach the belief either** — currently `_fold_lambda` folds whatever arrives.

### 12.10 Where the error lives AFTER B1b — and the next plan

Census re-run with `RIGEL_B1B=1`. Two things changed that the aggregate mwae did not show:

**The ship criterion is now essentially met.** `r(pr,|err|)` — the honesty diagnostic — was **positive on
every boundary class** at baseline (+0.09 … +0.67, "confidently wrong"). It is now **negative on nearly every
class**, and boundary message precision has collapsed from **0.62–0.90 → 0.02–0.06**, far below the 0.304
reference scale. Boundaries still correlate ≈ 0 (they remain unidentifiable) but they now say so **honestly
and weakly** instead of loudly and wrongly. That was the stated ship blocker.

**The RNA-isolated class partially solved itself.** capture-on corr **0.034 → 0.382** — from coin-flip to
meaningful. The relay now carries the flanking seams' gDNA into single-exon exons, which is exactly the
mechanism §6.5 predicted was missing.

**Remaining error is in REGIONS**, which hold ~88 % of mass at capture-on (`|err|` 0.225 / 0.241, corr
0.465 / 0.546). Boundaries are ~7 % and honestly weak. So the next lever is **not** boundaries.

#### The plan

| # | step | why now | risk |
|---|---|---|---|
| **N1** | **Regression tests for the FRAME facts** — the A1/A2 spliced-eff selector vs the brute-forced accumulator rule; the A3 per-face `c_b` denominator vs its closed form; the P0c directional peel on the owner's exon↔exon+junction topology (`TA+(1000,2000)(10000,11000)`, `TB+(1000,2000)(9000,11000)`). | Pure geometry — **behaviour-independent, zero risk, landable immediately**. Every one of these was silently wrong for months *because nothing pinned it*, and two of them were held in place by compensating errors. | none |
| **N2** | **Land B1b as the default**: flip `RIGEL_B1B`, resolve `test_mature_measurement_disagreement_silenced` (0.0612 vs 0.05 — re-derive the bound or accept the probe-depletion sensitivity), regenerate goldens **last**. Add tests for the relay through a zero-mass node, the intergenic wall, and the B1a `struct_gdna` / `no_data` split. | Converts a **measured** win into shipped error reduction: mwae 0.1396→0.1361, corr 0.671→0.688, guard −3.0 %, stranded neutral, and it removes the confident-wrong ship blocker. | low — already A/B'd on all 32 |
| **N3** | **Cleanup**: retire `RIGEL_P0_EXON_RNA` (measured net-negative; keep the finding in comments, delete the code path); decide `component_model.py` — its `edge_vehicle` has no production consumer and is superseded, but it is one rename away from being the **ACCEPT/RELAY predicate** §12.9 needs; drop the inert `_capture` vehicle census. | Removes three dead scaffolds and one misleading flag before they accrete more. | none |
| **N4** | **RNA relay — the mirror of B1b.** Running per-strand density with the §12.9 sink model, `ACCEPT`/`RELAY` split explicit, `fp[lsrc]` retired (a source without strand *s* has precision 0 automatically) and `fp[i]` moved to the destination. | B1b relayed **gDNA only**; RNA is still rebuilt from each node's own mass — the same defect, untouched. Regions hold ~88 % of mass and are where the error now is. This is the next error-down lever. | medium — gate as B1b was |

**Ordering rationale.** N1 is free and it locks in a session's worth of findings that nothing currently
protects. N2 banks the measured win. N3 clears the ground. N4 is the next real reduction. Stage C's full
`send`/`transmit`/`receive` restructure follows N4, by which point three of its four primitives
(density currency, the sink model, the fusion) will already be in production and tested.

---

## 13. N4 — THE RNA RELAY. Design + implementation plan (2026-07-23)

**Design only.** B1b relays **gDNA**; RNA is still rebuilt at every node from that node's own facing mass —
the same defect, untouched. Regions hold ~88 % of mass and are where the error now lives (§12.10).

### 13.0 THERE IS ONE RNA. We model ROUTING, not species. (owner, 2026-07-23)

> *We split "mature" and "nascent" for the algorithm's convenience, but we do not need to preserve the
> distinction. It is the same RNA — it just goes one way or another. All we model is: how much RNA goes
> **across** the splice junction, and how much flows **contiguously**. And the contiguous unspliced RNA is
> the only RNA that competes with genomic DNA — that is the whole reason the separation exists.*

This is not a rewording; it removes a question the design never has to ask. **We never need to estimate "how
much mature is in this exon."** A region lumps its RNA because there is nothing to split there. The only
place a split exists is a junction, and the junction **measures both routes directly** — the spliced channel
and its own unspliced crossing.

**That reframes both failed attempts.** P0 and P0b were species-subtractions: *remove the exon's mature*.
That framing needs a quantity nothing measures. Under routing the question becomes *what share of the RNA
arriving here continues?* — which is a **local, measured ratio**. Same arithmetic; a premise that exists.

Notation below writes `ρ_R` for the one RNA density and splits it only where a junction splits it, by route:
`ρ_R = ρ_cont + ρ_dep` (continues contiguously / departs through the spliced channel). Where earlier sections
of this document say "nascent" and "mature", read "the contiguous route" and "the departing route".

### 13.1 ⚠ The central flaw: RNA is NOT continuous, so B1b cannot simply be copied

B1b works because `ρ_g ≈ ρ_bg` is **locally uniform**, so "the neighbour's density is my density" is a sound
transport. **That premise fails for RNA at a junction — not because RNA is two species, but because the
routing splits it:**

```
    ρ_R(exon)              = ρ_cont + ρ_dep      both routes are present in the exon's unspliced pool
    ρ_R(boundary crossing) = ρ_cont              only the contiguous route is an unspliced crossing
    ρ_R(intron)            = ρ_cont
```

The unspliced density **drops by `ρ_dep` across the junction**. A naive per-strand copy of B1b would carry
the departing share into the intron as if it continued — which is **exactly the P0 regression**, re-entering
through a new door. Every design decision below follows from this.

**Corollary — why `c_b`'s gDNA-inclusive denominator is CORRECT.** Under routing,
`f_g(bnd)/f_g(exon) = (ρ_g+ρ_R)/(ρ_g+ρ_cont) = (D_B+S_B)/D_B` exactly, with `D_B` the **total** unspliced
crossing (gDNA included) — because it converts a *composition*, whose denominator includes gDNA. So the
gDNA-inclusion is right; the separately-flagged **FL-frame** error in `D_B` (dividing an RNA-containing mass
by the gDNA-FL opportunity, `D_B/(ρ_g+ρ_cont)=0.875`, §A3) is a different and still-open defect. Do not
conflate them.

### 13.2 The two directions are NOT mirror images — and that is forced, not incidental

| direction | operation | who holds the information |
|---|---|---|
| **exon → junction boundary** | **ATTENUATE** — only the contiguous share crosses | the **exon cannot** decompose its own RNA (a region has no spliced channel; `spliced_* ≡ 0`). The **boundary** measures the split as its own `S_B/D_B`. |
| **junction boundary → exon** | **SUM** — contiguous crossing **plus** departed spliced, which reappears in the exon | the boundary measures **both** terms directly |

Two corrections this forces to earlier documents:

* **`τ ∈ [0,1]` is wrong** (`message_layer_implementation_plan.md` §12.3). Going boundary→exon the ratio
  `(D_B+S_B)/D_B` **exceeds 1**. It is not a transmission coefficient at all in that direction — it is an
  **additive merge of two measured provenances** (item E). Only the exon→boundary direction is an attenuation.
* **"Send and forget" survives, but only because the RECEIVER applies the split.** The source genuinely cannot
  compute its own outgoing RNA density at a junction. The sender sends `ρ_r`; the receiving boundary applies
  the attenuation it measured. This is the same principle as §12.9's sink model — **all edge logic lives at
  the destination** — and it is why the design is coherent rather than a special case.

### 13.3 What N4 SUBSUMES — the cleanup is most of the value

The spliced mass currently enters the solve by **three** separate paths, and any two of them together
double-count (E.3 measured exactly that: `c_b` + graft = a systematic gDNA under-call at every exon):

| mechanism | what it patches | fate under N4 |
|---|---|---|
| `mature_dilution` / `c_b`, `_var_mat` | converts the exon's mature-inclusive composition to the boundary's mature-free one **on the gDNA mode only** | **retired** — with RNA relaying at the correct attenuation, `f_g = M_g/(M_g+M_r)` has the mature removed *by construction*, and the gDNA and RNA λ-factors become coherent (E.5 defect 2: they currently sum to 0.75, not 1) |
| the graft `+SPs/_esp` in `rho_pos/rho_neg` | puts the departed mature back into the exon | **retired** — it *is* the boundary→exon SUM, now expressed once |
| `absorb_p` / `absorb_n` | removes the source exon's mature | **retired** — it *is* the exon→boundary attenuation |
| the honest clamp `n_eff = ρ·E if ρ ≤ 1/E else smn` | stops a saturated subtraction becoming a confident zero | **retired** — with `ρ_μ` floored ≥ 0 inside the SUM, `ρ_r ≥ ρ_μ > 0` structurally, so "confidently no RNA" is unreachable (item E §5) |
| `pr += S` (measurement precision bolted onto a prediction's mode) | the only surviving RNA channel once `v_log = ∞` gagged the rest | **retired** — replaced by the share-weighted merge; this is what unblocks **R2** |

So N4 is not an addition — it replaces five patches with one mechanism.

### 13.4 The mechanism

Per strand `s`, mirroring B1b but with a **transport** term the gDNA relay does not need:

```
    own_s(i)  = ( log(f_s^loc · M_node/E_r^node) , 1/(Var(log f_s^loc) + 1/n_node) )
    run_s(i)  = fuse( own_s(i) , transport_s(src→i) ∘ run_s(src) )

    transport_s(src→i):
        not continuous (free_s fails either end, or a TSS/TES)   →  SINK: contributes nothing (§12.9)
        continuous, boundary is NOT a junction on s              →  identity          (~48 % of boundaries)
        exon → junction boundary                                →  × D_B/(D_B+S_B)   ATTENUATE, receiver-applied
        junction boundary → exon                                →  + ρ_μ = S_B/E_spl SUM (item E, share-weighted var)
```

`Var(log f_s^loc)` is already computed and returned (`rna_pos_frac_var` / `rna_neg_frac_var`).
`E_r^node` is the RNA-FL twin of `eff_global` and needs adding alongside it.

### 13.5 Flaw hunt — what could go wrong, and the stop conditions

| # | risk | mitigation / test |
|---|---|---|
| **R1** | **N4 re-runs P0 and regresses the same way.** P0 transported the exon's RNA and biased `f_g` down. | This is the design's **falsifiable prediction**: P0b failed *specifically* because `absorb_*` sat in a broken spliced eff-length frame, and **A2 fixed that frame while A3 fixed the `c_b` denominator (closing 55 % of P0's gap)**. N4 should now succeed where P0 failed. **Stop condition:** if the high-gDNA and stranded arms regress the same way again, the cause is the residual **`D_B` mixed-FL frame error** (`D_B/(ρ_g+ρ_ν) = 0.875`, §A3) — which needs the composition and therefore Stage C. Do not iterate; stop and report. |
| **R2** | **Nascent is less continuous than gDNA** (5′→3′ co-transcriptional gradient), so the transport is a weaker assumption. | Already measured: δ-census **family B** (intron ↔ junction, unspliced only) sits at median **−0.022** at capture-off, i.e. continuity holds well. But the RNA transport variance must be fit **separately** from gDNA's, not shared. |
| **R3** | **Double counting.** Four spliced paths exist; N4 adds a fifth unless the other four are removed **in the same commit**. | §13.3 is a checklist, not a wish-list. Gate: assert the emitted `f_g + f_r ≈ 1` per edge (E.5 measured 0.75 — a direct regression test). |
| **R4** | **`S_B` is a junction FLUX, not the exon's mature CONTENT.** They coincide only if the spliced density is a correct per-base density — which A2 has now made true *on continuing faces* — and if the junction is not probe-depleted. | Trap #3: a centred probe depletes junction reads ~300×, so `S_B → 0`, attenuation → identity, and mature leaks. This is the **already-xfailed** `test_mature_measurement_disagreement_silenced`. N4 must not make it worse; ideally it fixes it "for the right reason" (item E §5). |
| **R5** | **Three λ-factors must stay coherent** (gDNA, RNA-total, per-strand). | Regression test that the per-edge emitted fractions sum to 1. |
| **R6** | The **intergenic wall** has no RNA analogue — RNA is already sunk there by `free_s`. | Verify, don't assume: assert zero RNA precision on every intergenic-incident edge. |

### 13.6 Staging

| step | scope | risk | gate |
|---|---|---|---|
| **N4a** | Per-strand relay with **identity transport only** — continuous edges where the boundary carries **no spliced channel**. Pure B1b symmetry, no junction logic, nothing retired yet. | low | battery + guard; the junction path is untouched so `c_b` etc. stay valid |
| **N4b** | **Junction transport**: exon→B attenuation and B→exon SUM, retiring `c_b`, `_var_mat`, the graft, `absorb_*` **together**. | **high — this is P0's territory** | full battery per condition, stranded **and** unstranded; the R3/R5 coherence tests; R1's stop condition |
| **N4c** | Retire the honest clamp and `pr += S` (the share-weighted merge = item E), unblocking **R2**. | medium | `gdna_none` delta; the xfailed depletion test should now pass *for the right reason* |

**⚠ Measured correction to my own staging.** A first draft claimed N4a covers "820 of 1,699 boundaries, 48 %".
That counted seams *including* those carrying no RNA at all. Measured properly:

```
   boundaries total                          1699
   ... carrying RNA (free_s on some strand)  1267
   ... of those WITHOUT a spliced channel     261   <- N4a scope, 21 % of RNA-carrying boundaries
   ... of those WITH a spliced channel       1006   <- N4b scope, 79 %
```

So **N4a is not where the value is** — the junctions are, and they are the risky step. N4a is therefore a
**scaffolding** step, and should be sold as such: it lands the per-strand relay machinery on the safe 21 %
and proves it neutral-to-positive, so that **N4b changes only the transport term** rather than changing the
machinery and the transport in one commit. Expect N4a to be near-inert on the benchmark (as P0c, A2 and A3
each were); that is the point, not a disappointment.

**R6 pre-verified:** intergenic edges with a live RNA strand on **both** ends: **0**. RNA is already sunk at
the intergenic wall by `free_s`, structurally — checked, not assumed.

### 13.6b N4b DE-RISKED IN ADVANCE — the routing subtraction now WINS (measured 2026-07-23)

Under §13.0 the exon→junction direction *is* P0b: send `ρ_R`, subtract the measured departing density. P0b was
last tested against the **wrong A2 predicate** (66 faces corrected, not ~280). Re-ran it as a temporary probe
on top of corrected-A2 + A3 + B1b, then reverted:

| | mwae (all 32) | vs its baseline |
|---|---|---|
| P0b **before** the frame fixes | 0.1443 | **+0.0047 WORSE** than base 0.1396 |
| **routing subtraction, after A2+A3+B1b** | **0.1332** | **−0.0029 BETTER** than 0.1361 |

**A swing of 0.0076 — the frame fixes flipped the sign.** That is the strongest confirmation yet that the
A1/A2/A3 chain was the real prerequisite, and it retires N4b's principal risk (R1) *before* any code is written.

Per-slice, and the trade is real:

| slice | mwae | corr | scenarios |
|---|---|---|---|
| unstranded | 0.2175 → **0.2093** | 0.528 → 0.522 | **15 better / 5 worse** |
| capture verystrong | 0.2850 → **0.2702** | 0.253 → 0.260 | **4 / 0** |
| capture off | 0.1000 → **0.0966** | 0.831 → 0.832 | 7 / 5 |
| **stranded** | 0.0214 → **0.0261** ⬆ **(+22 % rel)** | 0.989 → 0.984 | **3 / 6** |
| ALL 32 | 0.1361 → **0.1332** | 0.688 → 0.683 | 18 / 11 |

**R1's stop condition is PARTIALLY triggered**: stranded still regresses, and aggregate `corr` dips even as
mwae improves. This is the P0 signature *attenuated, not eliminated* — consistent with the residual `D_B`
FL-frame error (0.875, ~14 % over-peel) being the remaining cause, exactly as R1 predicted. So N4b should be
built, but with the expectation that it trades stranded accuracy for unstranded, and the decision on that
trade is the owner's, not mine. **Do not iterate past it** — the fix is Stage C.

### 13.6c N4a — IMPLEMENTED, MEASURED, and its PREMISE REFUTED (2026-07-23)

Landed behind `RIGEL_N4A`, **default OFF**. Battery: **0.1361 → 0.1372, 0 better / 4 worse / 28 flat.**
Mildly negative, not the neutral scaffolding I predicted — and the dissection says why.

> **⚠ `free_s` is NOT a continuity condition for RNA.** It says a strand is *present* on both flanks. It does
> **not** say it is the **same transcript**. gDNA is a genomic background and is continuous; **RNA is
> transcript-specific and jumps at every structural boundary** — a transcript ends, another begins, a
> different exon set starts. A "no-spliced boundary" is a pure *signature change*, which is precisely where
> that happens. My reasoning that it was therefore the *safe* subset was backwards.

**Measured** (`scripts/debug/stranded_regression_dissect.py`, 12 stranded scenarios) — the observed RNA-density
step across boundaries where a strand is continuous, `|log ρ_R(right) − log ρ_R(left)|`:

| boundary class | n | median | p75 | p90 |
|---|---|---|---|---|
| **no-spliced boundary (N4a scope)** | 2,353 | **2.314** (≈10×) | 5.591 | 7.531 (≈1800×) |
| splice junction (N4b scope) | 6,213 | **4.460** (≈86×) | 7.395 | 8.606 |

Compare gDNA's step across measurable families: **0.02–0.09** (§6.5 E1). RNA is two orders of magnitude less
continuous. At a junction a large step is *expected* — it is exactly the routing split `ρ_dep` that N4b's
transport measures. At a **no-spliced** boundary there is no routing split, so identity transport predicts
≈0 and reality is ≈10×. **N4a's transport was simply wrong**, and the battery reflects it.

**Consequence for N4b — it is strengthened, not weakened.** N4b never assumes continuity: its transport is
*measured* from the junction's own spliced/unspliced split. The routing probe (§13.6b) already showed that
measured transport is a net win (0.1361 → 0.1332). **The lesson is that RNA transport must always be measured,
never assumed** — so the "identity transport" case should be deleted from the design rather than kept for
boundaries that lack a spliced channel. Those boundaries have no measurement, hence no licence to relay RNA.

### 13.6d The stranded regression — dissected, and my hypothesis was WRONG

The natural hypothesis — *messages are drowning a pristine strand signal* — is **refuted**. Mass-weighted
`|err|` by how strongly the incoming message precision outweighs the node's own strand Fisher information
`tau0` (stranded scenarios, N4a arm, overall `|err|` = 0.0227):

| regime | n | mass % | \|err\| |
|---|---|---|---|
| msg ≪ strand (<0.1×) | 1,007 | 13.9 % | **0.0311** |
| msg < strand | 3,815 | 29.8 % | **0.0314** |
| msg > strand | 4,562 | 11.5 % | 0.0186 |
| msg ≫ strand (>10×) | 10,556 | 10.3 % | 0.0204 |
| no strand evidence (`tau0 = 0`) | 9,888 | 34.5 % | **0.0138** |

**Error is HIGHEST where strand dominates and LOWEST where messages dominate** — the opposite of the
interference story. Nodes whose answer is driven by messages are *more* accurate, not less. So the stranded
regression is not messages overpowering strand; the residual error in stranded data lives in the
strand-solved nodes themselves.

**Honest limit of this measurement:** it is the error *distribution within one arm*, not a *paired per-node
delta between arms*, so it rules out the proposed mechanism but does not yet positively attribute the
+0.0013 (N4a) / +0.0047 (routing) stranded deltas. **The next measurement is a paired per-node diff** — same
scenario, both arms, ranked by mass-weighted error increase — which will name the exact nodes. That is a
small script and should run before N4b lands.

### 13.6e THE PAIRED DIFF — the regression DIAGNOSED (2026-07-23)

`scripts/debug/paired_node_diff.py`, 12 stranded scenarios, base → routing arm (`Δmwae = +0.00470`).

**Damage concentration** — 10,404 nodes harmed, 3,023 helped; worst 1 % of nodes carry 46 %, worst 5 % carry
85 %. Diffuse enough to be systematic, concentrated enough to name a class:

| class | mass % | damage % | \|err\| base → route | mean Δf_g |
|---|---|---|---|---|
| **intron adj-junction [1-strand]** | 11.8 % | **33.3 %** | 0.0166 → 0.0359 | **−0.0254** |
| **intron adj-junction [2-strand]** | 1.5 % | **12.7 %** | 0.1219 → 0.1784 | **−0.1089** |
| boundary (junction) [1/2-strand] | 3.2 % | 12.8 % | 0.042/0.181 → 0.069/0.211 | −0.058 / −0.228 |
| intron (not adjacent) | 0.7 % | 5.3 % | | |
| **exon adj-junction [1-strand]** | **47.6 %** | **0.8 %** | 0.0135 → 0.0136 | −0.0002 |
| exon (all other) | 14.4 % | 0.7 % | | |
| intergenic | 19.2 % | **0.0 %** | 0.0000 → 0.0000 | +0.0000 |

**The damage is at INTRONS and JUNCTIONS, not exons** — exons hold 57 % of mass and take 1.5 % of the damage.
`Δf_g` is **uniformly negative**: RNA is over-transported *into* introns and junction crossings, displacing
gDNA. Worst individual nodes: oracle 0.974 → base 0.927 → **route 0.376**.

**Root cause — the peel is under-sized, measured against the ORACLE.** Required peel =
`ρ_R(exon) − ρ_R(boundary crossing)`; actual peel = `S_B/E_spl`:

```
   measured / required :  median 0.9139   p10 0.3056   p25 0.6193   p75 1.0268
   under-sized (<1)   : 68.2 %      by >2x : 18.0 %      by >10x : 4.9 %
```

**The arithmetic closes.** A median 8.6 % residual of the exon's RNA, transported into an intron whose RNA is
~86× smaller (§13.6c Q1), over-states the intron's RNA by ≈ `0.086 × 86 ≈ 7×` — which is exactly the observed
`f_g` collapse.

### 13.6f THE ACTUAL DEFECT: the peel has NO variance. This is fixable.

**It is not a modelling misspecification and not message-vs-strand contention.** Note the shape: the peel is
*right in the median* (0.914) and *catastrophically wrong in a tail* (p10 = 0.31, 4.9 % below 0.1). A bias
would shift the median; this is a **precision** failure — junctions with few spliced reads peel just as
confidently as junctions with thousands.

And the code confirms it: the subtraction `rho_pos = f_r·sm/E_r − absorb_p` is applied as a **point estimate**.
`_pred_precision(n_eff, v_logfp, s2t)` carries the source's composition variance and the count — **nothing for
the peel**. (`c_b` at least had `_var_mat`; the routing subtraction has no equivalent.)

**The fix is item E's delta method, applied to a difference instead of a sum.** With
`ρ_cont = ρ_R − ρ_dep`:

```
    Var(log ρ_cont) = [ρ_R/ρ_cont]²·Var(log ρ_R) + [ρ_dep/ρ_cont]²·Var(log ρ_dep) ,   Var(log ρ_dep)=1/n_spliced
```

The weights **blow up as `ρ_cont → 0`** — i.e. exactly when the peel is nearly complete, which is the intron
case. A near-total subtraction *must* yield a near-worthless residual, and the solver would then ignore it
instead of crushing the intron's `f_g`. **This single term converts the regression from a wrong answer into an
honestly weak one**, which is the ship criterion.

> **Answer to the owner's question.** The stranded regression is **not** messages overpowering strand
> (§13.6d refutes that directly — error is *lower* where messages dominate). It is a **missing variance on a
> difference of two noisy densities**, whose relative error diverges as the difference approaches zero. There
> is a real, permanent contention between messages and strand, but this is not it — this is a bug, and it has
> a derived fix with no free parameters.

### 13.6g NEXT STEP (concrete)

**N4b-var:** implement the share-weighted variance on the routing subtraction *before* flipping the transport.
Order matters — the transport without the variance is what produced this regression.

1. `Var(log ρ_dep) = 1/n_spliced` — `spliced_n_*` is already plumbed.
2. `Var(log ρ_R^src)` — the source's composition variance, already available.
3. Combine by the delta method above; floor `ρ_cont` at 0 and let the variance carry the uncertainty.
4. Gate: the paired diff re-run — the intron-adj-junction classes must lose their damage share, and the
   worst-node `f_g` collapses (0.97 → 0.38) must become *weak* rather than *wrong*.

Only then re-run the routing A/B. Expected: the unstranded gain (0.2175 → 0.2093) is retained while the
stranded regression (0.0214 → 0.0261) shrinks toward neutral.

### 13.6h N4b-var — IMPLEMENTED, and my §13.6f DIAGNOSIS IS REFUTED (2026-07-23)

Implemented the delta-method variance on the routing difference (`_routing_precision`), behind `RIGEL_N4B`,
**default OFF**. Flag-off bit-identical. Then measured — twice, because the first attempt missed a channel.

| arm | stranded mwae | all-32 mwae |
|---|---|---|
| default (no routing) | **0.0214** | 0.1361 |
| routing, **no** variance | 0.0261 | 0.1332 |
| routing + variance on the per-strand factors | 0.0260 | 0.1331 |
| routing + variance on the **RNA-total** factor too | **0.0257** | 0.1345 |

**The stranded regression is essentially unmoved by the variance** (+0.0047 → +0.0044), and fixing the last
channel made the *aggregate* worse by damping the unstranded gain as well.

**Why the fix failed — I conflated precision with bias.** The delta method says the residual is *uncertain*;
it cannot say the residual is *too big*. And the residual **is** too big: the peel under-removes by a median
**8.6 %** of the exon's RNA, and the exon's RNA is a median **86×** the intron's (§13.6c Q1), so the leftover
over-states the intron's RNA by

```
    0.086 × 86  ≈  7.4×          exactly the observed f_g collapse (0.97 → 0.38)
```

A variance term does not move a mean. This is the very error this document repeatedly warns about — *a mode
bias paid for as variance* (§4.3, §8) — and I made it.

### 13.6i THE REAL DIAGNOSIS: the exon→junction RNA transport is ILL-CONDITIONED

`ρ_cont = ρ_R − ρ_dep` is a **difference of two large, nearly-equal numbers whose result is ~1 % of either** —
textbook catastrophic cancellation. The accuracy the peel would need:

| residual within … of the intron's true RNA | peel must be accurate to |
|---|---|
| 100 % | **1.16 %** |
| 50 % | **0.58 %** |
| 20 % | **0.23 %** |

The spliced measurement delivers a median of 91.4 % accuracy with a p10 of 30.6 %. **It is three orders of
magnitude away from what this subtraction requires**, and no variance, frame fix, or precision term changes
that — the conditioning is set by the 86× density ratio, which is physical.

> **Answer to the owner's question, corrected.** The stranded regression is **not** a bug, and it is **not**
> messages overpowering strand (§13.6d refuted that). It is an **ill-conditioned subtraction**: transporting
> an exon's RNA across a junction requires resolving a 1 % residual of a quantity 86× larger than the target,
> from a measurement good to ~9 %. That is a **structural** limit of the exon→junction RNA edge, not a
> misspecification we can tune away.

### 13.6j CONSEQUENCE — the exon→junction RNA edge should stay SILENT

This inverts N4's premise, and the current default already does the right thing for the right reason:
RNA is suppressed on exon→boundary. What N4 should keep is the *other* direction and the *other* source:

* **boundary → exon**: the SUM (contiguous + departed) is *well*-conditioned — both terms are measured at the
  boundary and neither is a small difference. Keep, and this is where item E's merge belongs.
* **intron → junction**: also well-conditioned (§6.5 family B, median δ = −0.022). The intron's nascent should
  inform the boundary, **not** the exon's leftover.
* **exon → junction**: **structurally silent for RNA.** gDNA still transports (B1b, already landed and a win).

So an intron's nascent should be informed by *its own* data and its *other* flank — never by an exon's
residual. That is a design conclusion with a measured basis, and it retires N4b as originally scoped.

**Next step:** re-scope N4 to `boundary → exon` + `intron → junction` only, and verify the ill-conditioning
claim holds on unstranded data too (where the 86× ratio may differ). The unstranded gain the routing probe
showed (0.2175 → 0.2091) came from *somewhere* — if it is not the exon→junction edge, it is worth isolating,
because it is the largest unclaimed improvement on the board.

### 13.6k ⚠ §13.6i WAS WRONG — I compared the wrong pair. Correction (2026-07-23)

§13.6i claimed the exon→junction subtraction is "ill-conditioned, amplification ≈86×". **That is an
arithmetic error.** I took the 86× from the exon↔**intron** RNA-density ratio (§13.6c Q1), but the message
subtracts at the exon↔**boundary** pair. Measured condition numbers (`ρ_R/(ρ_R − ρ_dep)` and `1 + S_B/D_B`,
all 32 scenarios):

| formulation | n | median | p75 | p90 | p99 |
|---|---|---|---|---|---|
| **mine** `ρ_R(exon) − ρ_dep` (RNA-only) | 9,169 | **1.1** | 1.3 | 1.8 | 21.9 |
| **owner's** `total_B − spliced_B` (total) | 22,032 | **2.15** | 5.56 | 11.13 | 74.85 |

**Neither is ill-conditioned**, and mine is *better* at the median — the opposite of both my §13.6i claim and
the prediction I opened this section with. The boundary's RNA is ~9 % of the exon's, so the residual is ~91 %
of the subtrahend. **My "catastrophic cancellation" framing is withdrawn.**

**Where the 7× really comes from — the CHAIN, not the subtraction.** Re-deriving with the right pair: the peel
recovers 91.4 % of the required amount, so the residual at the *boundary* is ≈1.9× too large — significant but
benign. That boundary error then propagates to the **intron**, whose true RNA is ~1.2 % of the exon's, and
`0.08·ρ_R(exon) / 0.012·ρ_R(exon) ≈ 6.7×`. **The amplification is in the transport from boundary to intron,
where the scale drops by ~86×, not in the subtraction at all.**

### 13.6l THE OWNER'S DESIGN — I cannot refute it, and it targets the real failure

Owner's proposal: compute **two** enrichment ratios at a junction and scale each incoming message by its own —

```
    r1 = density(exon region) / density(boundary)     TOTAL frame,     boundary INCLUDES spliced
    r2 = density(boundary)    / density(intron)       UNSPLICED frame, spliced excluded
```

resolving the boundary's composition (needed because density depends on composition through the FL) from the
**intron's** composition when the boundary has none of its own.

**This is right, and it is better than what §6.5 proposed**, for a reason worth stating:

* `r1` compares the exon's `{g,ν,μ}` against the boundary's **total**, which is also `{g,ν,μ}` — a
  **component-set-matched** pair. `r2` compares two `{g,ν}` pools — also matched. These are exactly the
  **family C and family B** edges whose `δ` I already measured as *unbiased* (§6.5 E1: capture-off medians
  **+0.008** and **−0.022**, against an unmeasurable family D at −1.116). **The estimator exists and is
  already validated.**
* Borrowing the **intron's** composition for the boundary's unspliced pool is sound for the same reason: both
  are `{g,ν}`. It is not an approximation of convenience; it is the matched-set condition.
* And it addresses the failure §13.6k just localized — **transport across a scale change** — rather than the
  subtraction, which was never the problem.

**The one caveat**, which is a caution rather than a refutation: `r1` and `r2` are ratios of *solved* or
*partially solved* quantities, so they carry the composition uncertainty into the scale factor. That must
propagate (it is the `V_e` of §4.3, now *measured per edge* rather than assumed), or we repeat exactly the
mistake of §13.6h — treating an uncertain correction as exact.

**Concrete next step.** Implement the two-ratio scaling at junction boundaries, taking the boundary's
unspliced composition from the intron message when its own precision is 0, and carrying each ratio's variance
into the message precision. Gate on the paired diff: the `intron adj-junction` classes (46 % of the damage)
must lose their damage share, and stranded must not regress. The measurement infrastructure for all of this
already exists (`paired_node_diff.py`, `enrichment_ratio_census.py`).

### 13.7 Codebase cleanup to land alongside

* **Delete `component_model.py` + `test_component_model.py`.** `edge_vehicle` implements the set-equality
  rule that §12 **superseded**, it has had no production consumer since the census was removed, and the
  ACCEPT/RELAY predicate N4 needs is a *different, simpler* function. Keeping the file invites someone to
  reuse a retired rule. (−2 files toward the ≤25-module budget.)
* **Extract one `_relay_fuse(run, own, transport, s2t)` helper.** The gDNA and per-strand RNA relays are
  structurally identical apart from the transport term; writing the fusion three times is how the mature
  accounting got out of sync in the first place.
* **`_scan` is ~250 lines with more comment than code.** Once the five patches in §13.3 are gone, move the
  superseded history (P0/P0b, geo-mean, gate-equality, τ-gag) **out of the hot loop** into this document and
  leave one-line pointers. The loop should read as: densities → transport → fuse → modes → fold.
* **`_pred_precision`'s last consumer** is the RNA prediction channel; after N4c it can be retired with the
  `τ` machinery (`StrandEvidence`, `tau_lam`, the Jacobian² conversion) in favour of the §3.4 EP quotient.
* **Retire `RIGEL_B1B`** once N4a lands and the pre-B1b path is no longer a meaningful A/B baseline.

### Then the restructure

```python
def _edge_logfactor(r, v):
    """The V->inf CONTRAST form: an edge constrains only the DIFFERENCES between the components it
    transports, because their common (unknown) enrichment scale is marginalized out.

        Q = sum_c r_c^2/v_c  -  (sum_c r_c/v_c)^2 / (sum_c 1/v_c) ,     r_c = log f_c(lam,theta) - a_c
        a_c = log(rho_c^src * tau_c * E_c^dst)        <- NOTE: no M_dst. It cancels in every contrast.

    psi += -0.5 * Q.   One line, no parameter, and every degenerate case is automatic (all verified):
        |C| = 1               -> Q = 0    a gDNA-only edge is SILENT under unknown enrichment
        any v_c = inf         -> that component drops out; a contrast needs BOTH ends known
        |C| = 2               -> Q = (r1-r2)^2 / (v1+v2)          the composition SHIFT, correctly
        |C| = 3               -> the full 2-contrast quadratic     (matches brute inverse to 5 d.p.)
        r += const            -> Q unchanged                       scale-invariant by construction
    """
```

Because `a_c` carries no `M_dst`, **the receiver never divides by its own mass** — `md`, `comp_fl` and the
whole M5 fractional-mass defect are absent *by construction*, and so is `_den`'s M6 collapse (there is no
normalizer to collapse). Both measured bugs are deleted by the rule that replaces them, not patched.

| # | step | change | gate | risk |
|---|---|---|---|---|
| **P2** | **`delta_estimator`** — a new ~60-line module: classify each edge into families A–D from the signatures, measure `δ_e` on A/B/C from observed densities, and take `V₀` = the population spread for D. Fit **once, belief-free**, before the sweep — same contract as `κ` and the overdispersions. Replaces `σ²_transfer`'s NPMLE projection. | new, pure | family medians ≈ 0 at capture off (E1); no solve touched | low |
| **P3** | **`send()`** — one node-level `(ρ, v)` triple per node (`node_global_geometry`, both faces pooled), carrying the **cavity belief** so P1's relay is structural. Deletes the per-face frame (`sf/df`, `EGs/EGd/ERs/ERd`, `MSs/MSd`, `_dst_mat/_src_mat`). | refactor + frame change | goldens move; census must not regress | medium |
| **P4** | **`transmit()`** — structural `τ_c` per edge from `free_s` / `mrna_active_s`, one small pure function. | pure, testable | unit tests on the §4.1 table | low |
| **P5** | **`receive()`** — the `Q` form above replaces `use_shift`/`use_shift_g`/`_mode_shift`/`_mode_density`/`comp_fl`/`_den`/`c_b`/graft/absorb/`_comb`. **The single largest deletion in the plan.** | the main behavioural step | full census + `gdna_none` + per-condition A/B, stranded **and** unstranded | high |
| **P6** | **Precision** — EP quotient (§3.4) replaces `StrandEvidence`/`tau_lam`/the Jacobian² relay. Retires `_pred_precision`'s `∞`/`return 0` and the honest clamp. | behavioural | §10.5's A/B | medium |
| **P7** | **Identifiability marking** (§5.2) — a node whose evidence precision is singular is marked **precision-0**, not defaulted to `f_g=1`. Handoff contract to pass-2. | honesty, not accuracy | `r(pr,err)` negative in **every** class | low |
| **P8** | **Single-exon holdout** (roadmap item J) — exclude family-D-only nodes from hyperprior training. §9 M1 shows they are already honest, so this is the *whole* treatment they need. | one mask | hyperprior fit quality | low |

**Out of scope for shipping pass-0**: the coverage gradient `π(E)/π(B_L)`, `γ_spl/γ_uns` (§6.3, deferred by
owner), and the seam-mode derivations R6. **R4 is subsumed** — `σ²_transfer`'s `(μ_d−μ_s)²` term is a
variance stand-in for `δ_e`, so it is *replaced* by P2's measured `δ` and removed with P5, never A/B'd on its
own.

Standing gates, unchanged: `gdna_none` phantom guard **as a delta** (working tree = **3,821,731**);
`pytest tests/calibration/ tests/native/` (313 passed, 2 xpassed); the grounded toy with **injected**
population priors (never toy-fitted); per-condition benchmark A/B on stranded **and** unstranded arms;
goldens regenerated last. No new constants: the plan introduces zero.
