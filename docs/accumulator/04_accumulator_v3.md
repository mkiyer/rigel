> ⛔ **SUPERSEDED by `05_accumulator_v5.md` (2026-07-28) — the standalone implementation spec, which folds in every decision and correction below. Do not cite this draft.**

# Accumulator v3 — contained nodes, flux edges, and the length moment

**2026-07-28. Written from scratch.** Supersedes `02_redesign_derivation.md` and
`03_path_accumulator.md` (both bannered). The path/cell store is **withdrawn**: it is expensive and no
consumer needs it. This version keeps the architecture Rigel already has — nodes own contained fragments,
edges own crossing — and changes four things.

Foundation: `../index/00_splice_graph_design.md` (the graph; unchanged and still correct).
Prototypes: `scratchpad/acc_proto_d.py` (the central result), `acc_proto.py`, `acc_proto_c.py`.

---

## 1. THE FOUR CHANGES

1. **TSS/TES become partition cuts**, so a node is internally homogeneous in RNA density. Today a TSS
   inside another transcript's exon creates no boundary at all, and the region spans a place where the RNA
   density genuinely steps. That is a modelling error in the *partition*, not merely a missing label.
2. **Splice junctions become real edges** from donor node to acceptor node, carrying the **full** effective
   length `fl_mean(RNA)`. Today a spliced fragment is split across two boundary objects, each with a
   region-size-dependent partial divisor.
3. **Edges keep OWNING their crossing fragments — but the currency becomes an integer count, not
   fractional mass.** This is what dissolves the multi-region partitioning problem (§2). Boundary
   two-sidedness, fractional mass, and the `n_cross` halving all disappear.
   ⚠ **This is NOT "all mass moves to the regions"** — that is `dag_design.md`'s second draft, which the
   owner rejected and so does this document. An edge still has its own data, its own density
   (`flux / fl_mean`), its own belief, and its own place in the chain. Nothing is demoted to flow-tracking.
   The word *flux* is used here in the sense the codebase already uses it (`boundary_flux_left/right`, the
   integer crossing count that today sits **beside** the mass) — v3 keeps that number and drops the mass.
4. ⭐ **Every object stores a COUNT and a RECIPROCAL-LENGTH SUM `Σ 1/L`.** The owner's idea, and it is
   worth more than a density convenience: at an edge `Σ 1/L` is a **divisor-free, FL-model-free** estimate
   of total density, and with the count it is a **complete two-component deconvolution** — prior-free,
   strand-free, 1.6× more precise than the `Σ L` variant an earlier draft substituted (§4).

**Explicitly not doing** (owner decisions): no path/cell store, no loopy BP, no per-node FL histograms,
no anchor rule.

---

## 2. WHY FLUX, NOT MASS — the partitioning problem dissolves

### 2.1 The problem

A fragment crossing `A → B → C` crosses two boundaries. Its slice in `B` is interior; today it is split
`b/(2L)` to each flanking boundary side so that per-fragment mass conserves to 1. The owner's objection:
*"how do you partition the fragment's length across those edges? The current solution is not right."*

### 2.2 It is the wrong question

A boundary's job is to estimate the local density. The number of fragments **spanning a genomic point** `p`
is a clean observable: a fragment `[s, s+w)` spans `p` iff `s ∈ [p−w, p−1]`, which is `w` start positions,
so

```
    E[ flux(p) ]  =  ρ · Σ_w f(w) · w  =  ρ · fl_mean            — exactly, for any FL pmf
```

**Integer observable. Exact Poisson variance. Independent of both flanking region sizes.** Nothing is
partitioned because nothing is deposited — the fragment crossing three regions crosses two points, and
`+1` at each is the correct count of crossing events at each.

*(This formula is already in the codebase — `effective_length.boundary_eff_length` returns exactly
`fl_mean` and its docstring derives it. The solver does not use it; it uses per-side mass instead.)*

### 2.3 The encompassed middle node

Node `B` fully spanned by the fragment gets **no contained count**, and that is correct — no fragment was
contained in it. Its two flanking edges each get `+1` flux, with divisors `fl_mean` that do not depend on
`|B|`. **There is nothing to share, so the sharing question does not arise.** `B` relays, and the messages
it relays are now properly normalised.

### 2.3b ⭐ THE DEPOSIT RULE BY NODE COUNT — the question that keeps coming up

A fragment **crosses** an edge whenever it spans that point. Deposit `(+1, +1/L)` on **every edge it
crosses**. If it crosses none, it is contained.

| fragment overlaps | edges crossed | deposit |
|---|---|---|
| **1 node** | 0 | node: `count += 1`, `recip += 1/L` |
| **2 nodes** | 1 | that edge: `+1, +1/L` |
| **3 nodes** | 2 | **each** edge: `+1, +1/L` |
| **4 nodes** | 3 | **each** edge: `+1, +1/L` |
| **N nodes** | N−1 | **each** edge: `+1, +1/L` |

**There is no partitioning at any node count.** Worked: a 200 bp fragment over
`A(50) → B(60) → C(50) → D(40)` crosses `A|B`, `B|C`, `C|D`; each gets `+1, +1/200`. Today the same fragment
deposits six fractional numbers — `50/200`, `60/400`, `60/400`, `50/400`, `50/400`, `40/200` — whose values
depend on the region sizes. That asymmetry *is* the partitioning problem, and it exists only because a
mass is being conserved.

**Why this is not double-counting.** Each edge answers its own question — *"how many fragments span
me?"* — and each has its own expectation `ρ·fl_mean`. A fragment over four nodes genuinely spans three
points, and each of those points genuinely has one more fragment over it. There is no total to conserve
across edges.

**Why it is not the "+1 per junction" bias retracted in `02_redesign_derivation.md` §0.** That retraction
was correct *in its own frame*: the path scheme declared "one fragment = one count on one path" and then
added +1 per junction to the same accounting, breaking its own rule. Here each object is an independent
observable, and `E[flux(J)] = ρ_mature·fl_mean` holds per junction regardless of how many other junctions
the fragment crosses. **This is also not new — it is what already ships:** `00_design.md` §6 invariant 3,
*"Flux is the integer per-side crossing count and is not split by the 50/50 rule — only mass is."*
v3 keeps flux, drops mass, adds `Σ 1/L`.

**The real consequence is CORRELATION, not bias.** `flux(e₁)` and `flux(e₂)` share fragments when the
edges are closer than a fragment length — a variance question for BP, computable as `P(w > |B|)`, and
present identically today. §11.2.

⭐ **A simplification falls out.** Every intron endpoint is a cut in the v8 partition, so a spliced
fragment's blocks always land in different nodes — **a spliced fragment can never be contained.** Nodes
therefore need only **2 channels** (unspliced ±), not 4. Contiguous edges still need 4 (a spliced fragment
crosses them inside an exon); junction edges need 1–2, the strand being known. (Exception: a fragment on an
*unannotated* junction — Q4.)

### 2.4 What conservation is replaced by

`Σ mass = n_fragments` goes away; flux is a crossing-event count, not a conserved mass. The replacement
invariants are stronger because they are integer:

```
    Σ_nodes contained_count            == n_contained_fragments        (exact)
    Σ_edges flux                       == Σ_f (#edges fragment f crosses)  (exact, computable at deposit)
    Σ_nodes contained_recip            == Σ_{contained f} round(2^32/L_f)   (exact, integer)
```

### 2.5 What is lost — nothing the solver uses ⭐

*(An earlier draft claimed dropping per-side mass would collapse the reframe's two faces and needed its own
A/B. That was wrong on two counts — owner-corrected, then verified in the code.)*

**The reframe's two faces are not left-mass vs right-mass. They are spliced-in vs spliced-out:**

```python
    def _rho_faces(fgc):                                  # bp_solver.py
        ru, rw = node_total_density(chain, geometry, fgc)
        rs = rw - ru                                      # the one-sided spliced density
        return (ru, ru + where(accept_l, rs, 0.0), ru + where(accept_r, rs, 0.0))
```

`ru` is already **node-pooled** — `node_global_geometry` returns `mass_left + mass_right` over
`eff_left + eff_right` for a boundary. The two faces differ only by whether `rs` is added.

**And per-side unspliced mass is already dead in the solve.** Grepped: `MS[0]/MS[1]`, `EG[0]/EG[1]`,
`ER[0]/ER[1]` occur **only** inside `_capture.update(...)`, the debug-only diagnostic dict; line 323 sums
`ER[0] + ER[1]` for boundaries. The only genuinely per-face structure in the production path is the
**spliced** channel (`SP[k]`, `SN[k]`, `ESP[k]`, and the `accept_*` gate built from them) — which is exactly
what §3 turns into junction edges.

**The swap is exact in expectation, not approximate:**

```
    E[m_l + m_r] / (e_l + e_r)  =  ρ  =  E[flux] / fl_mean
```

— the same estimand. Flux counts the *same fragments* with unit weight instead of overlap weight, and the
overlap distribution given crossing carries no further information about `ρ`, so **flux is the sufficient
statistic and mass is a lossy, region-size-dependent, fractional version of it.** `Σlength` (§4) then
recovers the length information that mass was garbling, cleanly separated from the count.

**So `_rho_faces` survives intact and both of its inputs improve:** `ru` becomes an integer flux over a
region-free `fl_mean`; `rs` gets the junction edge's full `fl_mean(RNA)` instead of the 2.7×-off partial
divisor. The one real behaviour change nearby is `accept_l/accept_r` moving from observational
(`SP+SN > 0`) to structural (does a junction edge attach on this side) — that is **C3**, already scoped
with its own A/B.

⚠ One genuine consequence remains, and it is a *gain*: Experiment C measured that overlap-weighted
attribution is a **left-looking FL-weighted smooth** of `ρ`, biased by up to +0.57 in `f_g` within one
fragment length of a density step. So part of the long-standing *"the boundary is a slope, not a cliff"*
effect is an artifact of the attribution rule rather than physics, and dropping mass removes it.

---

## 3. JUNCTION EDGES AND THE FULL EFFECTIVE LENGTH

A junction edge runs **donor node → acceptor node**, directed by the annotated motif strand, and a spliced
fragment deposits `(+1, +L)` on it. It deposits **nothing** on the contiguous edges it splices over — that
is the "half edges inside boundary nodes" this replaces.

### 3.1 The divisor — a trapezoid, not simply `fl_mean` ⚠ (v3 first draft over-stated this)

A mature fragment spans `J` iff it has `a ≥ 1` exonic bases on the donor side and `w − a ≥ 1` on the
acceptor side, **and both must fit inside the transcript**. With `R_d`, `R_a` the transcript's *cumulative*
exonic runway either side of `J` (a mature fragment reaches back through earlier junctions, so these are
transcript-wide, not the flanking exon alone):

```
    E[flux(J)] = ρ_mature · E_J ,   E_J = E_f[ min(w−1, R_d, R_a, R_d + R_a − w + 1)_+ ]
```

Computed on a `N(200, 50)` RNA FL:

| `R_d` | `R_a` | `E_J` | vs `fl_mean` |
|---|---|---|---|
| 5000 | 5000 | 199.0 | **1.00** |
| 550 | 550 | 199.0 | **1.00** |
| 300 | 300 | 198.2 | 0.99 |
| 200 | 200 | 160.1 | 0.80 |
| 147 | 147 | 87.8 | 0.44 |
| 100 | 100 | 19.6 | 0.10 |
| 50 | 5000 | 50.0 | 0.25 |

⭐ **`fl_mean` is exact for a typical internal junction** — median human transcript is 1099 bp, so a
mid-transcript junction has `R ≈ 550` each side and `E_J = 199.0`. It **tapers** where a runway is short,
i.e. at junctions near a transcript **end** — the same `reach` mechanism as §8.1, and the graph's TSS/TES
typing is again what makes it computable. ⚠ `R_d`/`R_a` are per transcript, so the isoform-ambiguity caveat
of §8.1 applies identically.

### 3.2 What this actually buys — corrected

⛔ The first draft claimed *"today's divisor is 73.5 vs v3's 200, a 2.7× difference."* **That comparison was
invalid** — it set today's per-side *mass* divisor against v3's *count* divisor, which price different
observables. Both are nominally unbiased. The real gains are three:

| | today | v3 |
|---|---|---|
| objects one spliced fragment touches | **2** boundary sides, each a fraction — using both double-counts, using one halves the data | **1** junction edge, integer |
| divisor | `E[min(ℓ,R)]/2` **or** a half-triangle, chosen by ~40 lines reconstructing "does the mature transcript continue past this flank" from region signature bits — documented wrong by **2–199× on ~26 % of junction faces** when the branch is picked wrong | `E_J` above, **read from the graph**, no inference |
| donor↔acceptor connected | no — the two flanks are unrelated objects | yes — one edge, one strand, one count |

This deletes `spliced_side_eff_length`, its continues/terminates selector, and the one-sidedness routing —
the largest single block of inferential machinery in `node_geometry.py`.

---

## 4. ⭐ THE LENGTH MOMENT — `(count, Σ 1/L)` identifies the two components

**⚠ v3 first draft stored `Σ L`. That was my substitution for the owner's original `Σ(1/fl)`, made
silently, and it is WRONG — `Σ 1/L` is better on every axis measured. Corrected here.**

### 4.1 The identity, and why `1/L` is the right weight AT AN EDGE

At a genomic point the crossing opportunity for a length-`w` fragment is **exactly `w` start positions**,
so depositing `1/w` cancels it identically:

```
    E[ Σ 1/L ]  =  Σ_c Σ_w ρ_c · w · f_c(w) · (1/w)  =  Σ_c ρ_c · Σ_w f_c(w)  =  ρ_total
```

⭐ **A divisor-free, FL-MODEL-FREE estimate of total molecular density.** No `E_c`, no effective length,
no FL pmf — and it holds for a *mixture* of components with different FL distributions. This is the
property `Σ L` does not have, and it is why the owner's original formulation is the right one.

The count supplies the second equation, and the pair identifies the split:

```
    N       = ρ_g · E_g[w]  +  ρ_r · E_r[w]
    Σ 1/L   = ρ_g           +  ρ_r
```

— identified iff the two **mean** lengths differ. The design's second row is `[1, 1]`, which is what makes
it well conditioned. Equivalently `N / Σ(1/L)` is the abundance-weighted mean fragment length, and `f_g`
falls out of a linear mixing equation.

⚠ **This is an EDGE property.** At a node the opportunity is `max(0, B−w+1)`, which `1/w` does *not*
cancel, so there is no divisor-free total there — `Σ 1/L` is simply another moment (and still the better
one, §4.2).

### 4.2 Measured (`scratchpad/acc_proto_g.py`) — gDNA FL 50, `f_g = 0.15`

| | RNA 200 | RNA 150 | RNA 100 | RNA 70 |
|---|---|---|---|---|
| `(N, Σ L)` sd / cond | 0.028 / 1137 | 0.027 / 726 | 0.029 / 452 | 0.043 / 459 |
| **`(N, Σ 1/L)` sd / cond** | **0.017 / 283** | **0.019 / 250** | **0.023 / 250** | **0.041 / 371** |

**1.6× more precise and 4× better conditioned at the 4× FL contrast** — the regime that matters — and
never worse. Both are unbiased (|bias| ≤ 0.001).

`ρ_total` from `Σ 1/L` **alone**, no model: **1.0024 / 0.9989 / 1.0004 / 1.0002 ×** truth.

**Robustness — the fitted gDNA FL wrong by 10 %:**

| | `f_g = 0.15` | `f_g = 0.50` |
|---|---|---|
| `ρ_total` via `(N, Σ L)` model | 0.9915× | **0.9652×** |
| `ρ_total` via `Σ 1/L` alone | **1.0012×** | **1.0000×** |
| `f_g` bias, `(N, Σ L)` | −0.0069 | −0.0137 |
| `f_g` bias, `(N, Σ 1/L)` | +0.0053 | +0.0174 |

The **composition** needs the FL model either way and is similarly sensitive. The **total density** is
*exactly* model-free under `Σ 1/L` and off by 3.5 % under `Σ L`. That is the robustness win.

At nodes it is a wash to a modest win: cond 147/184 at `B`=300, **242/126** at 1000, **275/124** at 5000.

### 4.3 Node vs edge remain complementary

Unchanged from §11b: a short node is gDNA-selective (a 200 bp RNA fragment cannot fit in a 147 bp node),
an edge is 4× RNA-selective. Each class is precise on one component and nearly blind to the other; the
2×2 Fisher reports which, and the solver must weight by it.

## 5. WHAT TSS/TES BUY

1. **Node homogeneity.** The partition's premise is that annotation is constant across a node. An interior
   TSS/TES violates it. 59.5 % of human transcript termini currently fall strictly inside a node.
2. **Structural zeros.** At a seam no strand-`s` transcript spans, `E_rna(edge, s) = 0` — so the edge's
   flux is **pure gDNA**, measured against `fl_mean(gDNA)`. This is the structural fix for the reframe
   defect (`gdna_reframe_terminus.md`): the face is *declared* RNA-free rather than *observed* to be, and
   no ratio against the destination's own belief is ever formed.
3. ⭐ **New clean gDNA anchors.** Every such terminus seam is a pure-gDNA measurement site with a
   region-size-independent divisor. There are 383,187 terminus positions on the human index; today the only
   clean anchors are intergenic regions. The gDNA hyperprior's substrate could grow substantially — worth
   measuring (Q3).
4. **A continuing SHARE at mixed seams.** ⚠ Honest scope: the binary *"does any RNA cross"* is already
   answerable from the signature (exon/intron bits on both flanks), and where a second transcript spans a
   terminus the signature's "RNA crosses" is **correct**. What the graph adds that the signature cannot is
   the per-strand **count** of transcripts starting / ending / continuing — which turns the binary into
   M10's peel share `w = ρ_ν/(ρ_ν+ρ_μ)`, currently estimated from data. Do not overclaim item 4; items 1–3
   are the load-bearing ones.

---

## 6. DATA MODEL

```
NODE  (per region, 2 channels: unspliced +/−  — a spliced fragment can never be contained, §2.3b)
    contained_count[c]     uint32
    contained_recip[c]     uint64      Σ round(2^32 / L_f)  over contained fragments

EDGE  (contiguous: src=i, dst=i+1  |  junction: donor node -> acceptor node)
    flux[c]                uint32
    flux_recip[c]          uint64      Σ round(2^32 / L_f)  over crossing fragments
```

* `L_f` is the fragment's **path length** — exonic bases, introns excluded — the same quantity in both
  equations of §4.1, for spliced and unspliced fragments alike.
* `Σ 1/L` is accumulated in **fixed point** (`2^32` scale). Range check: `L ∈ [20, 2000]` ⇒ each term
  ≤ 2.1e8; at 1e8 fragments the sum is ≤ 2.1e16 against a `uint64` ceiling of 1.8e19 — 800× headroom, and
  a single term keeps ~21 significant bits at the longest `L`.
* Channel counts differ per object class and this is deliberate (§2.3b): **nodes 2** (a spliced fragment
  cannot be contained), **contiguous edges 4** (a spliced fragment does cross them inside an exon),
  **junction edges 1–2** (the strand is known from the annotated motif) — Q6.
* **No mass. No `left`/`right`. No floats.** Memory: node `2×(4+8) = 24 B`; contiguous edge `4×12 = 48 B`;
  junction edge `2×12 = 24 B` — at 992 K nodes + 992 K contiguous + 404 K junction edges that is
  **~81 MB**, flat and predictable, with no hash map and no data-dependent growth.
* ⭐ **Integer accumulation is associative, so the result is bit-exact at any thread count.** This *removes*
  the known ~2.6 % cross-process nondeterminism (`calibrate_cross_process_nondeterminism`) rather than
  extending it — a strict improvement on today's float32 masses.

---

## 7. DEPOSIT

```
deposit(spans, channel):
    L = Σ span lengths          # R(L) = round(2^32 / L), fixed point
    walk the spans across the node cut array:
        if all spans fall in ONE node n:
            contained_count[n][ch] += 1 ;  contained_recip[n][ch] += R(L) ;  return
        for each contiguous node interface the fragment crosses:
            flux[e][ch] += 1 ;  flux_recip[e][ch] += R(L)
        for each inter-span gap:
            e = junction_edge(donor, acceptor, strand)      # index lookup
            flux[e][ch] += 1 ;  flux_recip[e][ch] += R(L)
```

Cheaper than today's deposit: no slice list, no per-slice division, no `n_cross` arithmetic, no float
stores — one integer pair per crossed object.

---

## 8. EFFECTIVE LENGTHS — the complete table

`F_c` = CDF, `S_c` = component support (reference for gDNA; pre-mRNA span for nascent; transcript for
mature). `reach_c(s)` = bases of `S_c` remaining from `s`.

| object | component | first moment (÷ count) | second moment (÷ Σ 1/L) |
|---|---|---|---|
| node, contained | gDNA | `Σ_w f_g(w)·max(0,B−w+1)` | `Σ_w f_g(w)·max(0,B−w+1)/w` |
| node, contained | nascent / mature | same with `f_r`, restricted to `B ∩ S_c` | same |
| contiguous edge | gDNA | `fl_mean(gDNA)` | **`1`** (the §4.1 identity) |
| contiguous edge | nascent | `fl_mean(RNA)`, **0 if no strand-`s` transcript spans the seam** | **`1`**, likewise |
| contiguous edge | mature | `fl_mean(RNA)` if mature crosses contiguously, else 0 | **`1`**, likewise |
| junction edge | mature | **`E_J`** (§3.1) | **`1`** |

### 8.1 ⚠ "REACH" — what it means, in plain terms

**Reach is simply: how many bases are left for the molecule to occupy, starting here.** It exists because a
fragment must *fit inside the thing it came from*.

* **gDNA** comes from the chromosome. A fragment can start anywhere and run on — reach is effectively
  unbounded. Every start position works for every fragment length.
* **Mature RNA** comes from a transcript, which **ends** at the polyA site. A fragment starting 50 bases
  before that end can only be ≤ 50 bases long — a 200 bp mature fragment starting there **cannot exist**.

So near a transcript end the RNA *opportunity* shrinks while the gDNA opportunity does not. Divide both by
the same length and you will conclude there is less RNA than there is — i.e. **you over-call gDNA**. The
correction is to count only the starts that could actually have produced a molecule:

```
    E_c(node) = Σ_{s ∈ node ∩ support_c}  F_c( reach_c(s) )
```

`F_c` is the fragment-length CDF, so `F_c(reach)` is *the probability that a length drawn from component c
fits in the space remaining*. Far from an end, `reach ≫ FL` and `F_c = 1`, so `E_c` is just the length —
this only bites within one fragment length of an end.

Measured on the real human index: mean `E_M/|A|` over **all** exonic bp is **0.893** at FL 200 — a **10.7 %
systematic gDNA over-call** if ignored — rising to a **+0.36** `f_g` error in the last region before a TES
(Experiment B). **Required, not optional**, and TSS/TES in the graph is what makes `reach` computable.

⚠ **Q1: what is `reach` for NASCENT RNA?** Mature has a hard end (the polyA site). Nascent RNA is
*partially transcribed* — its 3′ end is wherever polymerase happened to be, which differs molecule to
molecule, so there is no single position to measure reach against. **Recommendation: no taper for nascent**
(`F_V ≡ 1` inside the gene body). Justification: there is no hard fitting constraint, only a soft 5′→3′
coverage gradient, and a gradient is an *abundance* effect that is already absorbed into the per-node `ρ_V`.
⚠ It is an approximation — the constraint is length-dependent so it does not absorb perfectly — and it must
be measured on the `nrna_present` scenarios, where the oracle knows the answer.

---

## 9. ⭐ CAPTURE EFFECTIVE-LENGTH SHRINKAGE — proved unchanged

**The question.** *"If boundaries don't have their own mass, how do we shrink the effective length under
hybrid capture? Show me the equations."* Answer: **the shrinkage machinery is untouched.** Here is why,
from the code (`capture_eff_length.transcript_capture_eff_lengths`).

### 9.1 What the shrinkage actually computes

```
                Σ_{n ∈ t}  S_n · min(ρ_n / ρ_ref, 1)
    factor_t =  ───────────────────────────────────── ,        eff_em_t = fl_t · factor_t
                        Σ_{n ∈ t}  S_n

    then shrunk toward 1 on sparse evidence:   factor ← w·factor + (1−w),   w = C_t/(C_t + 1)
```

over three node classes on the transcript's own node set, with `ρ_ref` the global enriched-mode gDNA
density (`_global_reference_density`, fitted on contained region nodes):

| node class | support `S_n` | density `ρ_n` |
|---|---|---|
| contained region `r` | `E_g[max(0, L_r − ℓ)]` (`gdna_region_eff_len`) | `mass_gdna_contained[r] / S_r` |
| **pooled seam** between `r`, `r+1` | `½·(gdna_boundary_len[r] + gdna_boundary_len[r+1])` | `(mass_gdna_right[r] + mass_gdna_left[r+1]) / S_s` |
| splice-junction seam | `½·(side_len[jl] + side_len[jr])` | **imputed** `½·(ρ_left + ρ_right)` |

### 9.2 The two facts that make it immune

**Fact 1 — `S_n` is pure geometry × FL. The accumulator never enters it.**
`gdna_region_eff_len` and `gdna_boundary_len` are `effective_length.region_eff_length` and
`boundary_side_eff_length`: functions of a region length and an FL pmf, nothing else. **They are computed
identically under v3, bit for bit.** Whether the accumulator stores mass or counts is irrelevant to them.

**Fact 2 — the seam is ALREADY a single pooled node with a single density.** `_pooled_seam_arrays` forms
`m_s = mass_gdna_right[r] + mass_gdna_left[r+1]` — the two sides of the same boundary, summed — over the
averaged support. **The shrinkage never uses the two sides separately.** So there is no per-side quantity
here to lose.

### 9.3 The substitution, and the proof it is exact

Only `ρ_n` touches the accumulator, and only through the ratio. Per side,
`E[per-face mass] = ρ · E[min(ℓ, R)]/2` (`boundary_side_eff_length`), so

```
    E[m_s]  =  ρ·E[min(ℓ,L_r)]/2  +  ρ·E[min(ℓ,L_{r+1})]/2  =  ρ · S_s
    ⇒  E[m_s] / S_s  =  ρ                                        ← today's seam density
    E[flux_s] = ρ · fl_mean         ⇒  E[flux_s] / fl_mean  =  ρ  ← v3's seam density
```

**Same estimand, exactly.** Substituting one for the other changes `factor_t` in *variance* only — and in
the favourable direction, because flux counts the same fragments at unit weight and is the sufficient
statistic (§2.5). Per class:

| | today | v3 | status |
|---|---|---|---|
| contained region `ρ_r` | `f_g · (integer contained count) / S_r` | identical | **unchanged** — a contained fragment already deposits `+1` |
| pooled seam `ρ_s` | `m_s / S_s` | `f_g · flux_s / fl_mean(gDNA)` | same expectation, lower variance |
| junction seam `ρ_j` | imputed from flanking exon densities | identical | **unchanged** — a junction carries no gDNA, so there is nothing to measure there and the imputation stands |
| `S_n` (all classes) | geometry × FL | identical | **unchanged** |
| evidence `C_t` for `w = C/(C+1)` | a *mass* used as a count | an actual integer count | **strictly better** — the rule was designed count-based |

### 9.4 Capture-off bit-identity survives

Under uniform gDNA every `ρ_n = ρ_ref`, so `min(ρ_n/ρ_ref, 1) = 1`, so `num = Σ S_n = span_full` and
`factor = 1`. **This holds for any choice of `S_n` whatsoever**, so the capture-off bit-identity property —
the one the module's docstring leans on — is preserved by construction, independent of everything above.

### 9.5 ⚠⚠ A LIVE FACTOR-2 FOUND WHILE PROVING THIS — in shipped code, independent of the rework

Working §9.3 end to end turned up a discrepancy. Three lines, all from the repo:

1. `boundary_side_eff_length`'s own docstring: *"the per-face expected mass per unit ρ is … = min(ℓ,R)/2
   exactly, for every R"* ⇒ **`E[per-face mass] = ρ · gdna_boundary_len(L_face)`**.
2. `_pooled_seam_arrays`: `seam_m = mass_gdna_right[r] + mass_gdna_left[r+1]` — the two faces of one
   boundary ⇒ `E[seam_m] = ρ·(gbl[r] + gbl[r+1])`.
3. `_pooled_seam_arrays`: `seam_S = 0.5·(gbl[r] + gbl[r+1])`, and `calibrate.py:226` sets
   `boundary_eff_len = boundary_side_eff_length(...)` — i.e. `gbl` is **already halved**.

```
    E[seam_m] / seam_S  =  ρ·(gbl_r + gbl_{r+1}) / [½·(gbl_r + gbl_{r+1})]  =  2ρ
```

**Measured** (`scratchpad/acc_seam_check.py`, the current deposit rule driven under uniform density, repo
effective-length functions, 395 seams):

| region length | 2000 | 500 | 200 |
|---|---|---|---|
| shipped `m_s/S_s` ÷ truth | **1.994** | **2.002** | **1.981** |
| summed `m_s/(S_l+S_r)` ÷ truth | 0.997 | 1.001 | 0.990 |
| v3 `flux/fl_mean` ÷ truth | 0.997 | 1.001 | 0.987 |

`_gdna_region_node_arrays`'s **own docstring gives the correct formula** —
`S_s = ½·(E[min(ℓ,L_r)] + E[min(ℓ,L_{r+1})])`, which since `gbl = E[min]/2` is the **SUM** of the two
`gdna_boundary_len` values. The prose beside it says *"the AVERAGE of the two flanking per-side density
lengths gdna_boundary_len"*, and the code follows the prose. Formula and prose disagree; the formula is
right. This is precisely the ½-confusion `boundary_side_eff_length`'s docstring warns about
(*"anyone deriving a message variance from the old wording lands 2× off"*).

**Why no test caught it.** Under uniform gDNA the clip saves it: `min(m_s/ρ_ref, S_s) = min(2·S_s, S_s) =
S_s`, so the contraction factor is 1 and the *"factor-1-under-uniform bedrock"* invariant — the test
designed to catch exactly this class of error — **passes anyway.**

**Impact.** Under capture, a seam whose true density lies in `(ρ_ref/2, ρ_ref)` reads `≥ ρ_ref`, clips, and
**contributes no contraction when it should contribute some**. Seams systematically under-contract, over a
band exactly one octave wide. Both consumers are affected — `transcript_capture_eff_lengths` and
`priors._gdna_region_node_arrays` share the model deliberately.

⚠ **This is a live defect in `main` today, not something the rework introduces, and it is not mine to
fix.** It is a one-line change (`0.5*(a+b)` → `(a+b)`) with golden churn and an A/B, and it interacts with
`ρ_ref`, so it is an owner call. Recorded here because §9 claims the shrinkage carries over unchanged —
**it does, including this.** v3 does not inherit the bug (`flux/fl_mean` reads 1.00×), so landing v3 would
silently *change* the seam densities by 2× unless this is fixed first or the arm is scoped to expect it.

### 9.6 What DOES need attention

Not the shrinkage — the **partition**. `S_r = E_g[max(0, L_r − ℓ)]` collapses toward 0 as `L_r` shrinks, and
the v8 partition makes the median human exon region 147 → 102 bp. Contained region nodes therefore carry
*less* of `span_full` and seams carry *more*. `factor_t` stays well-defined (seam supports do not collapse),
but the **weighting between the two classes shifts**, so `factor_t` moves. **That is a real, measurable
change and it belongs to the index arm (V1/V2), not to the accumulator arm.** It must be measured against
the gDNA benchmark suite, where the siphon this module fixes is visible.

---

## 10. WHAT IS DELETED

`boundary_mass_left/right` · `boundary_flux_left/right` (replaced by one flux per edge) · the `mass`/`n`
duality through `substrate.py` and `node_geometry.py` · `boundary_side_eff_length` and its ½ (a property of
the deposit rule) · `spliced_side_eff_length` + the continues/terminates selector (~40 lines) · the
one-sided spliced routing · `boundary_junction_strand` and its `0`-is-identity worker merge (becomes
index-time) · terminal boundary objects · `BoundarySubstrate`'s two-sided view · the `_EPS` floors that
exist because divisors collapse.

---

## 11. HONEST LIMITATIONS

1. **Crossing observables are windowed.** `E[flux(p)] = Σ_w f(w)·mean(ρ over [p−w, p−1])` — a left-looking
   average over one fragment length. Only contained counts are strictly local. This is physics, it applies
   to the current design equally, and BP is what reconciles it. It is why contained counts still matter
   even though they are sparse.
2. **Adjacent edges of a short node are correlated.** A fragment spanning a 50 bp node crosses both of its
   edges. Their fluxes share fragments, so treating them as independent factors overstates precision. The
   overlap is computable (`P(w > |B|)`) and must be priced, not ignored. **Q2.**
3. **Dropping per-side mass removes the reframe's two faces.** `_rho_faces` builds a left and a right
   `ρ_tot`; with one density per object they collapse to one. That is a simplification, but it is a
   *behaviour change in the solver*, not just the accumulator, and it needs its own A/B.
4. **The length moment is a count-based claim on a Poisson design.** §4.3.

---

## 11b. ⭐ VERY DIFFERENT FRAGMENT LENGTHS (gDNA 50 / RNA 200)

The concern: short gDNA fragments cross far fewer boundaries, so edges see almost no gDNA. **Correct — and
it is a SELECTION effect that `E_c` prices exactly, not a bias.**

`E[flux_c] = ρ_c·E_c[w]`, so at an edge RNA is over-represented per molecule by
`E_r[w]/E_g[w] = 200/50 = 4×`. With a cfRNA-like truth of `f_g = 0.150` **of molecules**, an edge sees
gDNA at **0.042 of its events**. Measured (`scratchpad/acc_proto_f.py`):

| object | events | gDNA share of events | `f_g` estimate | bias |
|---|---|---|---|---|
| edge, thin | 178 | 0.042 | 0.148 ± 0.079 | **−0.0017** |
| edge, typical | 1,775 | 0.042 | 0.148 ± 0.027 | **−0.0017** |
| edge, deep | 17,751 | 0.042 | 0.150 ± 0.008 | **+0.0001** |
| node L=147 | 180 | **0.818** | 0.156 ± 0.037 | +0.0065 |
| node L=300 | 1,238 | 0.304 | 0.150 ± 0.009 | −0.0001 |
| node L=1000 | 8,235 | 0.173 | 0.150 ± 0.005 | +0.0001 |

**Unbiased at every depth; the cost of sparsity is variance.** ⭐ And the two object classes are
**complementary — they measure different components**:

| component sensitivity `E_g/E_r` | edge | L=100 | L=147 | L=300 | L=1000 | L=5000 |
|---|---|---|---|---|---|---|
| | **0.25** | **115.7** | 25.5 | 2.5 | 1.19 | 1.03 |

A short node is *massively* gDNA-selective — a 200 bp RNA fragment cannot fit in a 147 bp node — while an
edge is 4× RNA-selective. **Each class is nearly blind to one component and precise on the other.**
`L = 100` shows the failure mode honestly: 95 % gDNA events ⇒ the RNA column collapses ⇒ `f_g` bias +0.205,
sd 0.372. That node is a good `ρ_g` measurement and says nothing about `ρ_r`, and the 2×2 Fisher reports it.

**Consequences — none for the deposit, three for the consumers:**

1. ⛔ **Normalisation must NOT happen at deposit.** Dividing by `L` at deposit would require knowing which
   component the fragment came from — which is the thing being estimated. The normalisation lives in `E_c`,
   which is per component. *This is why we deposit raw counts and lengths.*
2. ⚠ **The solver must carry per-component precision**, because a single object is nearly blind to one
   component. `NodeBelief` already carries `var_pos/var_neg/var_gdna` separately, so the shape exists;
   `node_init` would gain the 2×2 Fisher from the moment solve.
3. ⚠ **FL accuracy matters asymmetrically.** Where a component holds 4 % of an object's events, an error in
   the *other* component's `E_c` dominates its estimate. Experiment D measured a 10 % gDNA-FL error costing
   0.0096 in `f_g` at a 60/200 contrast and a balanced mixture; **re-measure at 50/200 with `f_g = 0.15`
   before wiring** (M3).
4. ⭐ **An opportunity v3 creates:** junction-edge fragments are pure RNA by construction and intergenic
   contained fragments are pure gDNA, so each supplies a clean, unbiased FL moment
   (`flux/flux_recip` = the abundance-weighted mean length) for its own component — better substrate than today's
   region-type × compartment pools.

---

## 12. OPEN QUESTIONS

* **Q1 — `reach` for nascent RNA** (§8.1 explains `reach`). Nascent has no fixed 3′ end, so there is no
  hard fitting constraint to measure against. Recommendation: **no taper** (`F_V ≡ 1` in the gene body),
  because the 5′→3′ nascent gradient is an abundance effect already carried by `ρ_V`, not an opportunity
  effect. Alternatives: taper at the TES like mature, or model the gradient. Changes `E_V` at every TES —
  measure on `nrna_present`.
* **Q2 — correlated flanking edges at short nodes.** How is the shared-fragment overlap priced in BP?
* **Q3 — how many new clean gDNA anchors do terminus seams supply**, and does the hyperprior's substrate
  grow enough to matter? (§5.3)
* **Q4 — unannotated junctions.** A novel junction has no edge. Hold out, or deposit on a synthetic
  `(donor, acceptor, UNANNOTATED)` edge? Needs a real-cfRNA frequency measurement first; the toy has none
  by construction.
* **Q5 — is the SPANNING count worth one more integer per node?** Fragments that fully span node `B` form
  a fourth class with effective length `Σ_w f(w)·max(0, w−B+1)`, which is *large* exactly where the
  contained class is empty. It is not derivable from the two flanking fluxes. Deferred — §4.2 suggests the
  flanking edges already cover short nodes, so measure before adding.
* **Q6 — junction edge channels**: 2 or 4?

---

## 13. PHASES

```
V0  this spec reviewed; Q1 and Q6 answered (they change the schema)
V1  INDEX — 00_splice_graph_design.md, X0-X6 unchanged and still correct
V2  ⭐ P1g off the index alone: C1 (reframe at terminus seams), C2 (omega_graft per
    class), C3 (structural accept_l/accept_r).  Pre-registered A/B of P1G_SCOPE §8.
    This is the first measurable win and it needs no accumulator change.
V3  ACCUMULATOR — §6 schema, §7 deposit, C++ + reference; contiguous edges only,
    dual-writing today's payload.        Gate: legacy payload BIT-IDENTICAL 32/32
V4  JUNCTION EDGES + the full fl_mean divisor (§3).   Gate: its own A/B; the
    spliced-channel divisor changes ~2.7x, so this is a real arm
V5  EFFECTIVE LENGTHS — §8 table incl. the admissible-reach taper
V6  drop per-side mass; collapse the reframe faces (§11.3).   Own A/B
V7  the LENGTH MOMENT wired into node_init (§4).     Gate: REAL cfRNA, overdispersion
V8  deletions (§9), then goldens LAST
```
