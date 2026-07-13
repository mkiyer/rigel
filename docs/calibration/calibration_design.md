# Calibration deconvolution — design & derivation

**Status:** design document (branch `calib-ambig-init-wip`). This is the implementation-facing design for
the calibration deconvolution: the model, the equations, and — grounded in the current code — what must be
built. It **supersedes** the two narrative notes it grew out of (`calibration_bp_theory.md`,
`calibration_initialization.md`), which remain only as history. It does **not** authorize code changes; the
theory must be complete and sandbox-proven first (§10).

**Scope of this revision.** Five issues are open, in priority order:

1. **Initialization** — the value *and precision* every node holds before any message. *(designed — §5)*
2. **Boundary self-solve** from spliced/unspliced — the load-bearing strand-free gDNA/RNA readout. *(designed — §6)*
3. **Imputation (transfer) variance** — requires an enriched-vs-depleted node model. *(placeholder — §7.2)*
4. **Message interdependence** — a joint constrained integration so a weak message cannot flip a confident
   component. *(node solve designed — §4; the sweep-relay fix is a placeholder — §9)*
5. **Count conservation** — counts are fixed; total density is **not** (density depends on fragment length,
   which depends on composition). The node solve conserves counts, not densities. *(designed — §1, §4)*

This document develops #1 and #2 to implementation depth and states #3–#5 precisely, leaving the not-yet-designed
pieces as explicit placeholders.

---

# Part I — The theory (target design)

## 1. Model, notation, and the three quantities

A single-pass scan deposits fragment mass into **nodes** arranged as a linear **chain** per reference
(`…region — boundary — region…`):

- **regions** — exon / intron / intergenic interiors (fragments *contained* within);
- **boundaries** — the junctions between regions (fragments *crossing* the point).

Calibration **deconvolves** each node into up to three **components**:

- **gDNA** (`g`) — genomic contamination; present *everywhere*, its abundance modulated by capture enrichment;
- **RNA+ / RNA−** (`p` / `n`) — sense / antisense RNA (mature + nascent), present *only where a transcript of
  that strand is active*.

Calibration separates **RNA vs gDNA only**; the per-locus EM splits mature vs nascent downstream.

### 1.1 What is observed, what is latent, what is delivered

The single most important distinction in this design — and the correction that reshapes the node solve (#5) — is
between three quantities that were previously conflated:

| quantity | symbol | role | conserved? |
|---|---|---|---|
| **count** | `N_c` | fragments attributable to component `c`. Observed total `N` is a **hard fact**; the split is latent. It is also the **deliverable** (the EM prior is a count-mass split). | **yes:** `Σ_c N_c = N` |
| **density** | `ρ_c = N_c / E_c` | the **currency of comparison** across nodes: it normalizes counts by an effective length so a long region and a short boundary are comparable, and it is where spatial continuity lives (the gDNA field, the RNA coverage field). | **no** |
| **composition** | `f_c = ρ_c / D`, `D=Σρ_c` | the emergent **output**. Never a message, never conserved. | — |

**Observed per node.** For a **region**: the per-strand *unspliced* contained counts `(n_+, n_−)`, total
`N = n_+ + n_−`, and the region's genomic length `L`. For a **boundary**: the crossing counts, split into
*spliced* `(n_{s+}, n_{s−})` (motif-stranded) and *unspliced* `n_u`, on each side (left/right).

**Latent per node.** The count split `(N_g, N_p, N_n)` with `Σ = N`. Equivalently the count fractions
`a_c = N_c / N` on the simplex.

**Delivered.** The gDNA-vs-RNA **count** split feeds `assemble_priors` as the two per-locus Dirichlet scalars.
We *compare* in density but *deliver* in counts.

### 1.2 The effective length is per component — this is why density is not conserved (#5)

A component's count and density are linked by an **FL-marginal effective length**:

```
ρ_c = N_c / E_c ,   E_c = Σ_ℓ p_c(ℓ) · g(ℓ ; geometry)
```

where `p_c(ℓ)` is component `c`'s fragment-length pmf and `g` is the geometric acceptance:

- **contained region** (length `L`):  `g = max(0, L − ℓ + 1)`  ⇒  `E_c = E_{p_c}[max(0, L−ℓ+1)]`;
- **boundary crossing**:  a fragment of length `ℓ` crosses a point iff its start lies in a window of width `ℓ`,
  so `g = min(ℓ, L)` per side and the crossing effective length is `≈ μ_{p_c}` (the FL mean) for an
  interior point.

**gDNA and RNA have different fragment-length distributions** (`p_g ≠ p_r`), so `E_g ≠ E_r`. Therefore

```
D = Σ_c ρ_c = N_g/E_g + (N_p + N_n)/E_r      is NOT fixed as the split varies,
```

even though `N` is. Shifting one fragment from gDNA to RNA changes `D` by `1/E_r − 1/E_g ≠ 0`. **The conserved
quantity is the count `N`, not the total density `D`.** (Validated: `scripts/debug/count_conservation.py` TEST 3
— with `E_g=80, E_r=250` the count split `a_g=0.24` differs from the density composition `f_g=0.50` by 0.26.)

The old design (`calibration_bp_theory.md` §5, and `bp_theory.py`) imposed `Σρ_c = D`. That is the **equal-FL
special case** (`E_g = E_r`), where the count- and density-simplices coincide (TEST 2, exact match). The general
constraint is count conservation.

---

## 2. Precision is the currency

Beliefs are Gaussians in a log coordinate with a **mean** and a **precision** `π = 1/Var`. Precision — not the
value — is what everything combines on; a value without its precision is meaningless.

Two poles, both used at initialization:

| precision | meaning | examples |
|---|---|---|
| `π = 0` (Var = ∞) | **no information** — a placeholder that *will move* | a blank exon; the unstranded strand tilt |
| `π = ∞` (Var = 0) | **locked structural fact** — unmovable, absolute defense | intergenic gDNA; an off / inactive RNA strand |

> **Code convention (load-bearing, note the inversion).** The implementation stores the state as a **variance**:
> `var = 0` ⇒ `π = ∞` (locked); `var = ∞` ⇒ `π = 0` (blank).

The components live on the count simplex `Σ N_c = N`, so a node's **degrees of freedom** = (#active components
− 1): intergenic **0** (locked), single-strand **1**, AMBIG **2**.

Beliefs combine by precision-weighted sums (the Gaussian product). All log coordinates below are natural log.

---

## 3. Honest precision (derivation)

A component density is a count over an effective length, `ρ̂ = n / E`. Its precision comes from the count:

- **Poisson leading order.** `n | ρ ~ Poisson(ρ·E)` ⇒ `Var(log ρ̂) = Var(n)/E[n]² = 1/n` ⇒
  `precision(log ρ) = n` — **the count**. Given a global FL model, per-fragment lengths are *ancillary* to `n`,
  so a per-fragment `1/ℓ` estimate is a red herring; the FL **mean** normalizes, the FL **spread** shapes `E`.
- **Overdispersion cap.** Real counts over-disperse (GC/mappability, clustering, FL-model error). With a
  relative-variance floor `φ`:  `Var(log ρ̂) = 1/n + φ`  ⇒  `precision = n/(1 + nφ)`, **saturating at `1/φ`**.
  A large count cannot claim more precision than the overdispersion ceiling — this is the honest cap a bare
  Poisson `π=n` would violate.
- **Composition precision (the load-bearing quantity).** The gDNA/RNA log-odds `λ = log(ρ_g/ρ_c)` combines two
  independent log-densities, so `Var(λ) = Var(log ρ_g) + Var(log ρ_c) = (1/n_g + φ) + (1/n_c + φ)`. **This is
  nonzero for unstranded data whenever both counts exist** — it comes from the two components' *counts*, not
  from strand. It is exactly why a boundary can self-solve its gDNA/RNA split with no strand information (§6).

**Where fragment length enters precision.** The FL mean sets `E` (hence the density value). A boundary crossing
has a short `E ≈ μ_FL`, so it collects **few counts** and its count precision is inherently modest — the true
sense of "a short boundary is a variance sink." The FL spread enters `E`'s shape (matters for short regions) and
`φ` (FL-model error). It does **not** independently inflate the count variance.

---

## 4. The node solve — a joint constrained MAP with count conservation (#4, #5)

Each node forms, per component, a single log-density target `(m_c, π_c)` by the precision-weighted product of
its local belief and all incoming messages. It then solves the split **jointly, under count conservation** —
never one component at a time.

### 4.1 The objective

```
minimize over  (N_g, N_p, N_n ≥ 0)   subject to   Σ_c N_c = N :

    F(N)  =  − ℓ_strand(n_+, n_− | N ; κ, over-dispersions)          [ intrinsic strand likelihood ]
             + Σ_c  π_c · ( log ρ_c − m_c )² ,   ρ_c = N_c / E_c      [ density beliefs + messages ]
```

Two clean layers meet here:

- **Count space** owns the observation and the constraint. The strand term `ℓ_strand` is the two-component
  gDNA/RNA Beta-Binomial mixture likelihood of the observed `(n_+, n_−)` given the split — the **only intrinsic
  gDNA/RNA signal** (each component predicts an expected `+`-fraction: gDNA `½`, RNA+ `κ`, RNA− `1−κ`; both
  overdispersions applied symmetrically, so κ=½ is uninformative). The constraint `Σ N_c = N` is count
  conservation.
- **Density space** owns the cross-node information: `m_c` are the per-component log-density targets from the
  prior and neighbour messages (§7); `π_c = 0` for an active-but-uninformed component.
- The **link** `ρ_c = N_c/E_c` is the only coupling, and it is where the per-component FL enters.

### 4.2 The KKT self-defense derivation (now in log-count space)

Reparametrize by `x_c = log N_c`. Then `log ρ_c = x_c − log E_c`, so the density term is
`Σ_c π_c (x_c − m̃_c)²` with the **FL-shifted target** `m̃_c = m_c + log E_c`. The constraint is
`Σ_c e^{x_c} = N`. Setting aside the smooth strand term, the Lagrangian is
`L = Σ_c π_c (x_c − m̃_c)² + μ(Σ_c e^{x_c} − N)`:

```
∂L/∂x_c = 2π_c (x_c − m̃_c) + μ e^{x_c} = 0
   ⇒   δ_c ≡ (x_c − m̃_c) = −(μ/2) · e^{x_c}/π_c = −(μ/2) · N_c / π_c .
```

`μ` is shared (it enforces the single count sum). So **each component's log-count deviation from its target is
`δ_c ∝ N_c/π_c`** — inversely proportional to its precision. This is **self-defense**, unchanged in form from the
density case but with the lever now the component **count** `N_c`:

- a **high-π (confident)** component barely deviates — it defends itself;
- a **low-π (weak)** component deviates freely and **absorbs the count-conservation constraint**;
- the `N_c` factor is the log-lever: a large confident component can supply/absorb some *count* at a small
  log-deviation — finite-precision protection is **proportional**, not absolute;
- **absolute** protection is the `π=∞` structural lock (`δ_c → 0`).

Validated: `scripts/debug/count_conservation.py` TEST 1 — confident gDNA (π=50) deviates `δ_g = +0.010`, weak
RNA+ (π=2) absorbs `δ_p = +0.555`, and the implied `μ` is identical from both components (−0.03187), confirming
the KKT.

### 4.3 Interdependence — solve the DoF, never a component (#4)

The count constraint makes the components a **pie chart**: you cannot move one without the others moving. The
inference is over the **degrees of freedom** (single-strand 1, AMBIG 2), solved **jointly**.

The failure to avoid is the **Trojan horse**: updating the addressed (weak) component independently, then making
a confident component the residual — which drags the confident component even though nothing informed it. The
joint solve is the fix: by §4.2 the confident component moves only in proportion to `(other π)/(its π)`, and in
the AMBIG case a message aimed at one weak component is absorbed by the *other* weak component, sparing the
confident one. (Demonstrated in `bp_dependency.py`: naive 30→25.6 vs joint 30→27.3.)

**Reconciliation of differing totals.** Messages set per-component *targets*; the node's own observed `N` is the
hard constraint; the KKT routes the residual count to the lowest-π components. A boundary with 10 crossing
fragments does not overwrite a region's 100 contained fragments — it informs the *split*; the region allocates
its own `N`. Composition is emergent, never transmitted.

### 4.4 Uniqueness

The 1-D (single-strand) and 2-D (AMBIG) solves are empirically unimodal across the regimes tested (excess,
deficit, two-confident conflict, extreme precision gap). The clean monotonicity proof requires
`2π_c + μ e^{x_c} > 0`, which can fail in the excess (`μ<0`) regime; unimodality still holds empirically there,
so the production solver should use **bisection on `μ` or a grid**, not a bare Newton, in that regime. *(This
was verified under density conservation; it must be re-verified under count conservation — §10 open item.)*

---

## 5. Initialization (#1 — primary)

Every node self-solves a local belief **before any message**. Init runs in three layers over the two precision
poles; the output per node is a count-split belief with an honest per-component precision.

### 5.1 Default — the blank node
Every node starts at `{f_p, f_n, f_g} = {0, 0, 1}` (100 % gDNA) at **precision 0**. The value is irrelevant;
π=0 means "no information, will move."

### 5.2 Structure locks (from the region signature, never from counts)

The 4-bit region signature → `strand_class ∈ {POS, NEG, NONE, AMBIG}` → the lock pattern:

- **Intergenic regions and intergenic↔exon boundaries (TSS / TES):** pure gDNA `{0,0,1}` at **π=∞**. They are
  (a) **RNA sinks** — an RNA message dies at them (RNA cannot exist intergenic); (b) **barriers** — they neither
  move nor relay composition across themselves, so one gene cannot leak into the next; (c) sources for the
  *population* gDNA background (the depleted mode of the prior), **not** a per-node message into a genic
  neighbour.
- **Single-strand nodes** (POS or NEG): the absent RNA strand's count is locked at `0` (π=∞) — a permanent sink.
  The remaining axis (the ON strand vs gDNA) is the 1-DoF quantity to solve.
- **AMBIG nodes:** all three components active; 2 DoF; no strand lock.

*(The exact signature-bit → class → lock mapping is grounded against the code in §8.)*

### 5.3 Strand-tilt self-solve (from the per-strand counts)

Each node self-solves an initial gDNA/RNA **tilt** from the two-component strand Beta-Binomial likelihood
(`ℓ_strand`, §4.1) of `(n_+, n_−)`, with a precision that **scales with strand-specificity**:

- unstranded (κ ≈ ½): the likelihood is flat ⇒ **π → 0** (the tilt carries no information);
- stranded (κ ≫ ½): the imbalance is informative ⇒ **high π**, scaling with the discriminability
  `w(κ) = (2κ−1)²` and the count.

For the flagship unstranded case this contributes nothing — correct, and precisely why the boundary self-solve
(§6), which is strand-free, is load-bearing.

### 5.4 The mandatory message-free self-solve on all nodes

Before any message, **every** node runs its local self-solve (§5.3–§6) to set its `(count-split, precision)`.
This must run on **all** nodes, including locked and empty ones — a locked all-gDNA seam must still *emit* its
confident gDNA belief (the "silent locked seam" bug, §8/§9). The local solve is **likelihood-only** here (no
gDNA hyperprior), so init establishes honest precisions the first message pass then combines.

---

## 6. Boundary self-solve from spliced vs unspliced (#2 — primary, load-bearing)

A boundary is the only node that directly observes two physically distinct fragment classes at one genomic
point, giving a **strand-free** gDNA/RNA readout.

### 6.1 The two classes

- **Spliced = mature RNA, with known strand.** A splice junction carries a single-stranded genomic motif
  (canonical GT/AG), so a spliced fragment's strand is **fixed by the motif**, not inferred. Placed directly into
  the correct RNA-strand density:
  ```
  ρ_p = n_{s+} / E_spl ,   ρ_n = n_{s−} / E_spl        (E_spl = crossing eff-length under the RNA FL)
  ```
- **Unspliced = gDNA (+ sparse nascent).** Nascent RNA crossing an unspliced junction is rare, so at init the
  unspliced crossing count is assigned to gDNA:
  ```
  ρ_g = n_u / E_uns                                     (E_uns = crossing eff-length under the gDNA FL)
  ```
  With stranded data the unspliced count can be strand-deconvolved into nascent vs gDNA (Beta-Binomial, its own
  precision); unstranded, it is all gDNA. *(The nascent split under `nrna_present` is a placeholder — §9.)*

Note `E_spl` (RNA FL) `≠ E_uns` (gDNA FL): the count→density map differs by class, so count conservation (§1.2)
applies at the boundary too.

### 6.2 The self-solved composition + honest precision

The boundary self-solves a full `{ρ_p, ρ_n, ρ_g}` with the honest density precision of §3, **with no prior and
no strand**. Worked example (`+` junction, `μ_RNA = 200`, `μ_gDNA = 100`, `φ = 0.02`):

```
n_{s+} = 100 spliced(+)  → ρ_p = 100/200 = 0.5
n_u    = 300 unspliced   → ρ_g = 300/100 = 3.0
λ = log(ρ_g / ρ_p) = log 6 = 1.79        ⇒   f_g = σ(1.79) = 0.857
Var(λ) = (1/300 + φ) + (1/100 + φ) = 0.053   ⇒   π(λ) ≈ 18.8    (NONZERO, no strand used)
```

The spliced:unspliced **ratio sets the direction** of the split; the counts set the honest precision. Boundaries
are the only nodes that localize the gDNA/RNA composition without a prior in unstranded data.

### 6.3 Per-side structure and the flux into regions

A boundary owns its spliced RNA **one-sided** (a junction faces one exon). Its per-side gDNA/RNA flux feeds the
adjacent region's deconvolution as the neighbour density message (§7). *(The exact per-side geometry and the
region-projection are grounded in §8.)*

### 6.4 Pass-0 enriched mode from self-solved boundaries

Because boundaries self-solve an honest, prior-free composition, the gDNA-density prior's **enriched mode** can
be built directly from them **before** the first message pass. Select the structurally clean gDNA-enrichment
readouts — **intron↔exon, single-strand, enriched (high unspliced density), RNA-free (low spliced density)** —
and take their adjacent exon densities as the enriched training sample. (Validated on the flagship: such
boundaries' adjacent exons are 99.7 % gDNA at density ≈ 26 vs the true plateau ≈ 30; `gdna_none`/capture-off
select **zero**, so it self-gates. This replaces the depleted-only KDE trained on already-crushed pass-1
densities.) This is a *consequence* of the self-solves, not a hand-set floor.

---

## 7. Message currency & propagation (summary; #3 is a placeholder)

### 7.1 Currency — per-component density fields (settled)

**Composition is not transmissible** (the active set changes along the chain, so composition jumps where a
strand switches on/off), and **absolute totals are not** (the receiver reallocates to its own `N`). The currency
is **per-component density**, each a field with its own continuity:

- **gDNA** — active everywhere; a smooth log-density field **modulo a multiplicative capture enrichment factor**;
  its message is `ρ_g`, transferable across any edge, carrying an **enrichment-transfer variance**.
- **RNA+/RNA−** — coverage fields active **only along their transcript**, exactly 0 elsewhere; the message is
  `ρ_s`, transferable only where strand `s` is active on **both** endpoints; a deactivation is a sink.

Message precision is honest, capped, and decays per hop:

```
π_msg(c, i→j) = 1 / ( 1/π_c^belief(i) + σ²_transfer(c, i→j) )
```

— never more precise than the source knows, nor than the edge reliably carries. A message relayed `k` hops has
`π = 1/(1/π_src + k·σ²_hop)`, so far nodes are rightly less sure. The chain is a forest of linear paths, so BP is
**exact in one forward + one backward pass**, each node emitting the message derived from its belief *excluding*
the message it received from that direction (true tree BP), followed by the joint solve (§4).

**Self-defense protects the confident, not the blank.** A blank node (π=0) adopts whatever dominates; the
flagship collapse was a blank enriched exon dominated by a **confident-but-wrong depleted prior**, not by
messages. Hence the first pass must be **prior-free** so the blank node is governed by honest messages + its own
count constraint. This makes message precision honesty non-negotiable: a stronger message is safe **only** if it
is genuinely more reliable.

### 7.2 The transfer variance — PLACEHOLDER (#3, requires an enriched/depleted model)

`σ²_transfer(g, i→j)` is the remaining quantitative piece and is **not yet designed**. What is established:

- The current production model uses a **single global** total-density disagreement variance `σ²_imp`. On real
  region↔region pairs it is ≈ the true gDNA transfer variance (not composition-inflated as an earlier synthetic
  suggested), because the gDNA density genuinely **jumps** at every exon↔intron enrichment boundary.
- The transfer variance is **bimodal by enrichment regime** (measured on the flagship, oracle-validated):
  within-regime (both-enriched `Var(Δlog ρ_g) ≈ 0.30`, both-depleted `≈ 0.25`) is **tiny**; across-regime
  (enriched↔depleted, the enrichment jump) is **huge** (`≈ 37`). A single global value is a bad average of a
  0.25-vs-37 bimodal.
- The regime is **observable** from the two nodes' total densities (no deconvolution), which **breaks the
  circularity**: classify each edge as within- or across-regime from the observed densities, and assign the
  stratified variance. Within-regime → strong messages (junction→exon); across-regime → weak (intron→exon).

**Open:** formalize the enriched-vs-depleted node model that assigns `σ²_transfer` per edge from observables
(the density-regime stratification above is the current best candidate); the RNA continuity leg; the
irreducible unobservable-probe floor. To be designed and sandbox-proven before implementation
(`scripts/debug/transfer_var_estimator.py` has the regime-stratification evidence).

---

# Part II — Status quo & what must be implemented

*(This part is grounded in the current calibration code — filled from the status-quo grounding pass. Pending.)*

## 8. Current calibration behavior (grounded, with `file:line`)

Grounded against `node_geometry.py`, `simplex_logodds.py`, `simplex.py`, `strand_likelihood.py`, `signature.py`,
`effective_length.py`, `substrate.py`, `strand_deconv.py`, `bp_solver.py`, `fl.py`.

### 8.1 State and the precision convention

Each chain node (region *and* boundary) carries one `NodeBelief` = the composition **fractions**
`(f_pos, f_neg, f_g)` over the node's **unspliced** mass, plus a per-component **variance**
`(var_pos, var_neg, var_gdna)` = `Var(log f_c)` (`node_geometry.py:199-217`). The variance *is* the precision
state, with the inversion of §2: `var = 0` ⇒ locked, `var = ∞` ⇒ no information. The unspliced mass `M` and all
effective lengths are **static** (`build_node_statics`, `build_node_geometry`, computed once); the sweep mutates
only the fractions and variances.

### 8.2 Initialization (`init_beliefs` → `_type_belief`)

Init runs in two steps (`node_geometry.py:421-468`):

1. **A likelihood-only strand solve on every node** — `_solve_nodes_logodds_all(...)` with **no global prior and
   no imputation message** (`node_geometry.py:448-462`): the bare per-node strand Beta-Binomial + a Jeffreys
   reference. This matches the §5.4 "likelihood-only self-solve" requirement.
2. **A signature-binary override** (`_type_belief`, `node_geometry.py:267-314`) keyed on the admissible-strand
   booleans `free_pos/free_neg` (a region's own ± transcript bits; a boundary's ± continuity across the seam,
   `node_geometry.py:336-341, 359-360`):

   | class | mask | belief | precision |
   |---|---|---|---|
   | default | (all) | `{0,0,1}` | `var_g=∞`; `var_p=∞ if free_pos else 0`; `var_n=∞ if free_neg else 0` (`:288-294`) |
   | **G1** (neither free) | `~free_pos & ~free_neg` | `{0,0,1}` | `var_g=0` **locked** (`:312-313`) — intergenic / TSS-TES sink |
   | **G2** (one free, `M_unspl>0`) | `free_pos ^ free_neg & M>0` | strand-deconv posterior `gdna_frac, rna_*_frac` | the deconv `*_var` (`:305-310`) |
   | **G3** (both free = AMBIG) | `free_pos & free_neg` | `{0,0,1}` | `var=∞` (no info) — resolved later by the sweep |

The strand solve carries the count only through the Beta-Binomial variance (Fisher information) — the
**count-zero-information** principle. `κ = rna_sense_frac`, and both overdispersions enter the variance
symmetrically (`simplex.py:33-45`), so κ=½ is uninformative.

**Init gaps (`node_geometry.py`):** (i) a single-strand node with `M_unspl == 0` (G2 but not `g2_active`) falls
silently into the *same* no-information hold as AMBIG rather than being a locked/solved single-strand node
(`:296-298`); (ii) the strand overdispersions default to `0.0` in the signature (`:428-429`) — a silent Binomial
fallback if a caller omits them; (iii) `logodds_window = 10.0` is a magic clamp (`:432`); (iv) the G2 branch
writes the forbidden strand's *fraction* from `deconv` but its *variance* from a different mask, relying on the
solver returning ≈0 on the forbidden axis.

### 8.3 The node solve — a pure-fraction log-odds MAP (no explicit conservation)

`_solve_nodes_logodds_all` reparametrizes the 2-simplex to the single gDNA-vs-RNA log-odds
`λ = logit(f_g) = log ρ_g − log ρ_rna` and grids `λ` on a fixed window `[−L, L]`, `L=10`, `f_g = σ(λ)`
(`simplex_logodds.py:74-78`). Single-strand is an exact 1-D solve (the live strand carries all `1−f_g`); AMBIG
adds an inner tilt `τ∈[−1,1]` with `f_pos=(1−f_g)(1+τ)/2`, `f_neg=(1−f_g)(1−τ)/2`, marginalized by logsumexp
(`:155-157, 313-318`). The local log-posterior ψ is the additive sum (`_local_loglik_logodds`, `:127-197`):

```
ψ(λ) =  strand-mixture BB loglik            [simplex.py:33-45; the only intrinsic gDNA/RNA signal]
      + strand_obs · Jeffreys Beta(½,½)     [:162-168]
      + global_logprior[grid]               [:170-171; pre-evaluated, count-space; built in bp_solver]
      − ½ · gdna_imp_prec · (log f_g − gdna_imp_mode)²          [:176-181; message on log-FRACTION]
      − ½ · Σ_s rna_imp_prec[s] · (log f_active/f_s − rna_imp_mode[s])²   [:182-191; on log-FRACTION]
      + ( log f_g + log(1−f_g) )             [:196; change-of-variable Jacobian]
```

Extraction: `f_g` = posterior **median** over the grid (`:81-90, 247`), `f_pos/f_neg` = posterior **means**
(`:249-250`); the message-currency precisions are `Var(log f_c)` moment-matched on the grid (`:256-265`).

**The constraint / conservation reality (issue #5).** The solve is over **pure fractions**, with the simplex
enforced *structurally* by the `σ(λ)`/`τ` parametrization — **not** by a fixed total density `D` and **not** by
an explicit count sum. Component mass is recovered post-hoc, `gdna_mass = f_g·M_unspl`,
`rna_mass = (1−f_g)·M_unspl + mass_spliced` (`:273-282`). Because all components share the one unspliced mass
`M`, the *unspliced count* `f_c·M` is conserved by construction; the spliced mass is added to RNA **outside** the
ψ objective. So: production already conserves the unspliced count in its projection — **my old sandbox's
`Σρ_c = D` (density conservation) was the erroneous frame, not production's.** What production lacks is the clean
joint per-component-**density** message integration of §4; its messages are instead log-**fraction** penalties
(below).

**Node-solve gaps:** the ψ has **no spliced-floor term** despite the module docstring claiming one
(`simplex_logodds.py:11-12, 144` are stale — do not cite as an existing term); `L=10`, `_PRIOR_EPS=1e-3`, and the
AMBIG `1/(n+1)` RNA-floor are fixed constants; the AMBIG `f_g` median is quantized on the coarse `n_grid=60` grid
(single-strand de-quantizes on a finer `n_grid_ss`); AMBIG uses the neutral measure `log f_g + log(1−f_g)`, not
the exact uniform-simplex Jacobian (a deliberate anti-gDNA-prior choice).

### 8.4 Message currency and per-component effective length (issue #5 substrate)

Effective lengths are geometric functionals of a FL pmf, computed **per component** (`node_geometry.py:108-116`):

```
region_eff_length(L, pmf)      = E_pmf[max(0, L−ℓ)]        = L·F(L) − S(L)        [effective_length.py:57-73]
boundary_side_eff_length(pmf,L)= E_pmf[min(ℓ, R)]          = S(R) + R·(1−F(R))    [:76-99]   (crossing density divisor)
spliced_side_eff_length(pmf,L) = E_pmf[min(ℓ,R)²/(2ℓ)]     = ½S(R) + ½R²(H_tot−H(R))  [:102-128] (one-sided; → μ/2)
```

gDNA uses `gdna_fl_pmf`, RNA (nascent + spliced) uses `rna_fl_pmf` — **separate EB-shrunk pools**
(`fl.py:152-196`; gDNA = intergenic+intronic, RNA = spliced-annotated; `POOL_EB_PRIOR_ESS=1000.0`), so
`E_gdna ≠ E_rna` is real and deliberate. The **message currency** is a per-component density that divides the
*same* face mass `M` by the *component-specific* `E` (`node_geometry.py:246-253`):

```
ρ_g = f_g·M/E_gdna ;   ρ_p = f_pos·M/E_rna + spliced_pos/E_spl ;   ρ_n = f_neg·M/E_rna + spliced_neg/E_spl
```

**This is issue #5 in the code:** one mass `M`, three divisors ⇒ `ρ_g + ρ_p + ρ_n` is **not** invariant under a
fraction shift, even though `f_c·M` (count) is. The FL asymmetry is baked into the currency, exactly as §1.2
derives.

*Caveat:* `node_densities` is documented as "the sweep's message currency" but is in fact consumed **only** by
the Phase-2 KDE substrate (`gdna_density_prior.py:136`); the sweep's `_scan` recomputes message densities inline
with a different construction (nascent + mature *measurement* − absorption, `bp_solver.py:835-879`). Two parallel
density formulas coexist.

### 8.5 Boundary handling (issue #2)

A boundary is solved on the chain **exactly like a region** — its `f_g` comes only from (a) the strand
Beta-Binomial over its **pooled unspliced** counts (`u_pos=max(left,right)`, `mass_unspl=left+right`,
`node_geometry.py:344-347`), (b) the global prior, (c) neighbour messages. **Its spliced mass never enters its
own ψ solve:** `node_sweep` passes `_zero_spl` into every per-node solve (`bp_solver.py:616, 632`), and there is
no spliced ψ term (§8.3). Spliced instead does two things only:

1. rides as a **fixed additive RNA term** at projection: `rna_mass = (1−f_g)·M_unspl + mass_spliced`
   (`bp_solver.py:1002, 1033`; all spliced counted as RNA with certainty, no gDNA share);
2. is **emitted** in the outbound RNA message to the exon flank as a mature *measurement*
   (`ρ = n_nasc/E_rna + n_mat/E_spl − ρ_mat_dst`, `bp_solver.py:839-874`), routed to the strand set by the
   **junction motif** `js` (`node_geometry.py:143-156`).

The boundary's solved `f_g` then splits each flanking region's side crossing mass into gDNA/RNA for the per-locus
prior (`chain_boundary_side_deconv`, `bp_solver.py:1009-1039`; note it drops the per-strand RNA split, setting
`rna_pos_frac=rna_neg_frac=0` on the sides). One-sided spliced eff-length uses the half-triangle `E_spl ≈ μ_RNA/2`
to avoid the ~2× mature-density understatement (`effective_length.py:102-128`).

**This is the core of issue #2:** a boundary that is obviously mature RNA (large spliced mass) but has
balanced/ambiguous unspliced counts **cannot use its own spliced evidence** to pull `f_g` down — the strand-free
gDNA/RNA readout of §6 does not exist in the solve.

### 8.6 The sweep (brief)

The sweep is a forward-backward pass over the chain (`bp_solver.node_sweep`): each node's running belief is
seeded from its local solve, gDNA and per-strand RNA density messages are emitted to neighbours at the
message precision `pr = n_src/(n_src·σ²_imp + 1)` (honest count ⊗ imputation), and the global prior + Phase-2 KDE
enter every per-node ψ. *(Grounded partially — from `bp_solver.py:835-879` and the existing docs; the dedicated
sweep-driver grounding reader did not complete, so §8.6 is the one area to re-verify directly before touching
sweep code.)* The theory-level issues it carries are the **Trojan-horse relay** (#4, §9) and the **not-prior-free
first pass** (§9).

## 9. The gap and the implementation plan (init + boundary first)

### 9.1 Spec → status → change

| # | design (Part I) | current status | change | where |
|---|---|---|---|---|
| A | **Boundary self-solve** from spliced+unspliced (§6) | **missing** — spliced never enters ψ (`_zero_spl`); `f_g` from unspliced only | add a spliced/unspliced density self-solve → boundary `{ρ_p,ρ_n,ρ_g}` + honest precision | `node_geometry` (statics), `simplex_logodds` (ψ term), `bp_solver:616,632` |
| B | **Count-conservation** joint solve with per-component `E` (§1.2, §4) | fraction-space solve conserves unspliced count via `f·M`; per-component `E` exists (§8.4); but messages are log-**fraction** penalties, spliced added outside ψ | make the message integration a joint per-component-**density** solve; carry `E_gdna≠E_rna` explicitly in the count↔density map; fold spliced into the conserved count | `simplex_logodds:127-197`, `bp_solver:835-879` |
| C | **Prior-free first pass** (§7.1) | **diverges** — global prior + depleted floor in *every* ψ; pins enriched exons | remove the pre-imposed floor; let depleted + enriched modes emerge from self-solves (§6.4); replace the floor's bleed-stopping role with honest messages | `bp_solver` global-prior path |
| D | **Locked seam emits confident gDNA** (§5.4) | silent — locked seam emits `π≈0`; local solve skips locked/empty | promote `emit_locked`; local solve must run on all nodes | `bp_solver` (`emit_locked`) |
| E | **Joint relay, no Trojan horse** (§4.3) | per-component running belief relayed without sum-reconciliation | sum-reconcile (joint) before relay; exclude the recipient's own message | `bp_solver._scan` |
| F | Strand-tilt precision → 0 at κ=½ (§5.3) | rises with κ but a reference floor keeps it finite | make it vanish cleanly, or document as intended | `gdna_strand.py`, `simplex_logodds` |

### 9.2 Ordered plan — initialization + boundary first (#1, #2)

The user's priority is #1 (init) and #2 (boundary), which are the same first move: **give boundaries the
strand-free spliced/unspliced self-solve (§6)**, since that self-solve *is* the load-bearing init for the
unstranded case. Concretely, in order:

1. **Boundary self-solve (change A) + honest precision.** Implement §6 in the init/statics path: a boundary's
   spliced (motif-stranded, `E_spl`) and unspliced (gDNA, `E_uns`) counts produce `{ρ_p, ρ_n, ρ_g}` with the
   density precision of §3 — no prior, no strand. Prove in the sandbox first (extend `bp_sandbox.py`), then wire.
2. **Count-conservation in the message integration (change B).** Recast the node's message reconciliation as the
   §4 joint solve over the count simplex with per-component `E` (the FL-shifted targets `m̃_c = m_c + log E_c`).
   Re-verify uniqueness under count conservation (§10 open item 1) before implementing.
3. **Locked-seam emission (change D)** and the **relay fix (change E)** — both are message-path correctness fixes
   that the boundary self-solve depends on to propagate.
4. **Prior-free first pass (change C)** — the highest-risk change; only after A/B/D/E are proven to carry the
   deconvolution, and only with the 24-condition bleed check (§10 open item 4).

Changes to #3 (transfer variance) and the full #4 message model follow, once §7.2 is designed.

## 10. Sandbox validation status

- **Self-defense, recovery, self-gating, relay/decay, Trojan-horse, AMBIG 2-DoF, uniqueness** — proven under the
  *density-conserved* solve (`bp_theory.py`, `bp_reconcile.py`, `bp_dependency.py`).
- **Count conservation (#5)** — self-defense KKT preserved (identical `μ`), reduces to the old solve at equal FL,
  and exposes the count-vs-density FL bias (`count_conservation.py` TESTs 1–3, this doc §1.2, §4.2).
- **Transfer-variance regime stratification (#3)** — evidence gathered (`transfer_var_estimator.py`); model not
  yet designed.

**Open before any code:**
1. Re-verify node-solve **uniqueness** under count conservation (§4.4) — was proven only for density conservation.
2. Design the **transfer-variance / enriched-vs-depleted** model (§7.2) and prove it in the sandbox.
3. The **nascent-vs-gDNA** unspliced split under `nrna_present` (§6.1).
4. Prove the **prior-free first pass** does not reintroduce RNA bleed across all 24 conditions before removing
   the floor.
