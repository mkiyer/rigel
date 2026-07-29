> ⛔ **SUPERSEDED by `../accumulator/05_accumulator_v5.md` (2026-07-28). The path/cell store this document marginalises is WITHDRAWN, and its §4 mass-allocation ladder is resolved there. Do not cite.**

# Paths → calibration: marginalization, mass allocation, and BP on a splice graph

**v1, 2026-07-28.** The third of three. `../index/00_splice_graph_design.md` builds the graph;
`../accumulator/03_path_accumulator.md` counts paths over it; **this document is how calibration consumes
them**, and it is the least settled of the three by design — the index and accumulator are mechanical, this
is a modelling change.

Prototypes: `scratchpad/acc_proto.py` (A, B), `scratchpad/acc_proto_c.py` (C).

---

## 1. THE PROBLEM

Paths nest. The owner's example: `A = [1,2,5,6,7]` (spliced), `B = [1,2]`, `C = [2]` (both unspliced).
`C ⊂ B ⊂ A` as node sets. The solver works per node and per edge. So we must turn a set of nested,
variable-length path counts into per-node and per-edge quantities, without counting a fragment twice and
without a length bias.

Two things make it tractable:

1. **Fragment length is handled by dividing by it.** A fragment of length `L` occupying node `n` with
   overlap `o_n` contributes `o_n/L` — the owner's `1/L` density, made exact by the overlap weight.
   §2 proves this is unbiased with an effective length of exactly `|n|`, for any FL pmf.
2. **A path's interior nodes are always FULLY covered.** The graph cuts at every exon endpoint, so a
   junction donor is always a node boundary. Only the two end nodes have variable overlap. §3 turns that
   into three accumulators per path.

---

## 2. THE TWO ATTRIBUTIONS, AND WHY THIS IS A CHOICE

A path count is **one** observation touching several unknowns. Turning it into per-node evidence needs an
attribution rule. There are exactly two coherent ones and they are not equivalent.

### 2.1 COVERAGE attribution (the owner's `1/L`)

Node `n` accumulates `V_n = Σ_f o_{f,n} / L_f` over **every fragment that touches it**.

*Unbiasedness.* At a genomic position `x`, fragments of length `w` covering `x` have `w` admissible starts,
so their rate is `ρ·w·f(w)` and each contributes `1/w`:

```
    E[ Σ_{f ∋ x} 1/L_f ]  =  Σ_w ρ · w · f(w) · (1/w)  =  ρ · Σ_w f(w)  =  ρ
```

Summing over `x ∈ n` gives `E[V_n] = ρ·|n|`. **Effective length `|n|`, exactly, for any FL pmf — so it
cancels from the composition.** This is the property the owner was reaching for and it is real.

*The catch.* The same calculation with a spatially varying `ρ` gives

```
    E[ Σ_{f ∋ x} 1/L_f ]  =  Σ_w f(w) · mean( ρ over [x−w+1, x] )
```

— a **left-looking, FL-weighted local average of `ρ`.** The coverage statistic estimates a *smoothed*
density, with a smoothing width equal to that component's mean fragment length. gDNA smooths over ~60 bp
and RNA over ~200 bp, so **the two components are smoothed differently and the composition is biased at
any sharp density step.** An exon boundary under hybrid capture is exactly such a step.

### 2.2 ANCHOR attribution

Node `A` accumulates the count of paths whose **first** node is `A`, divided by the admissible-start
effective length `E_c(A) = Σ_{s ∈ A ∩ S_c} F_c(reach_c(s))` (`02_redesign_derivation.md` §2.2).

Exact under a Poisson-start process — `E[N_p] = ρ(a(p))·E_c(p)` depends on the anchor and nothing else —
so it is unbiased **at every node, including at a step**. It uses each fragment once, so it is noisier.

### 2.3 The measurement (Experiment C)

gDNA `ρ=0.02` FL 60 everywhere; RNA `ρ=0.05` FL 200 on a contiguous block only — a sharp step at each end.
True `f_g` = 0.286 inside the block, 1.000 outside.

| | anchor `f_g` bias | coverage `f_g` bias | anchor CV | coverage CV |
|---|---|---|---|---|
| **100 bp regions**, block interior | **+0.0000** | +0.0000 | 0.344 | **0.251** |
| **100 bp regions**, at the step | **+0.0000** | **+0.155** | — | — |
| **25 bp regions**, block interior | **+0.0000** | +0.0000 | 0.79 | **0.29** |
| **25 bp regions**, first region inside the block | **+0.000** | **+0.565** (0.851 vs 0.286) | — | — |

The coverage bias decays over ~10 regions at 25 bp — i.e. **one RNA fragment length**, exactly as §2.1
predicts.

> **Verdict.** Coverage is **1.4–2.7× more precise** and **biased by up to +0.57 in `f_g` within one RNA
> fragment length of a density step**. Anchor is exactly unbiased everywhere and correspondingly noisier.
> Short regions sit at exon boundaries; exon boundaries under capture are density steps. **The coverage
> attribution's bias lands precisely where calibration is hardest.**

⭐ **This identifies a defect that is already shipping.** Today's boundary mass *is* the coverage
attribution (`o_n/L` split across boundary sides). The long-standing note that *"the boundary is a slope,
not a cliff — three enrichment ratios, not one"* (`density_composition_reconciliation.md` §4) is this
smoothing, observed from the other side. It has been treated as a physical effect to be modelled; §2.1
says a large part of it is an **artifact of the attribution rule**.

### 2.4 ⭐ A second, independent argument against the coverage attribution

`V_n = Σ_f o_{f,n}/L_f` is **fractional**. So is every per-strand split of it. But the strand likelihood is a
**Beta-Binomial on integer counts** (`strand_likelihood.py`), and its statistical power is `n_obs` — which
is why `node_geometry` already carries `n_unspl_*` *and* `mass_*` as parallel arrays and its own comment
says `1/mass` is "wrong in two independent ways".

> **Adopting the coverage attribution re-creates the mass/flux duality that this whole redesign exists to
> delete** — for the strand channel, which is the only intrinsic gDNA/RNA information source there is.
> The anchored attribution keeps every observable an integer, end to end.

That was not obvious when §2.1 was written and it is a stronger argument than the bias measurement, because
it is structural rather than regime-dependent.

⚠ Not yet measured: how much of the *real* enrichment slope is physical (probe capture genuinely retains
molecules overlapping a probe) versus attribution blur. Both exist. **M-cal1.**

---

## 3. THE EXACT MARGINALIZATION

Per path `p`, per channel, the accumulator stores `N_p`, `D_p = Σ 1/L_f`, `F_p = Σ o_first/L_f`
(`03_path_accumulator.md` §3). Let `S_int(p) = Σ_{interior nodes} |n|`. Then, **exactly**:

```
    COVERAGE   interior node n :  V_n += |n| · D_p
               first node      :  V_n += F_p
               last node       :  V_n += N_p − S_int(p)·D_p − F_p
    ANCHOR     node a(p)       :  n_A += N_p
    EDGE FLUX  every edge on p :  φ_e += N_p          (integer)
```

The first/last identity follows from `o_first + o_last = L_f − S_int(p)`, summed after dividing by `L_f`.
Self-check, per fragment: `Σ_n o_n/L = 1` — accumulator test A8.

**Per-node effective lengths** (`02_redesign_derivation.md` §2.2, §6): `|n ∩ S_c|` in the interior of the
component's support, tapering within one FL of a support end.
⚠ **The taper applies to BOTH attributions.** Experiment B measured it for the anchor (`f_g` over-called
by +0.36 in the last region before a TES); the coverage statistic has the same cause and a different shape
— at the transcript's final base only fragments *ending exactly there* can cover it, one start per length,
so its effective length collapses by `E[1/w]`. Neither attribution escapes it; the graph's TSS/TES typing
is what makes it computable.

---

## 4. ⭐ THE MASS-ALLOCATION DECISION — nodes vs edges

The owner's framing: today edges **own** mass so that fragments are not double-counted on regions; the
alternative is a **flux** model where edges only track flow. Here is what the code actually does and what
each option costs.

### 4.1 What boundary nodes do today (code audit)

Boundary nodes are **first-class variables**, not bookkeeping. Four distinct roles:

| # | role | where |
|---|---|---|
| B1 | **own composition belief** — solved from their own crossing strand counts `u_pos/u_neg = max(left,right)` | `node_geometry._boundary_strand_stats`, `init_beliefs` |
| B2 | **the reframe faces** — `ρ_tot` per side feeds `r = ρ_tot(dst)/ρ_tot(src)` | `bp_solver._rho_faces`, `MS = (mass_left, mass_right)` |
| B3 | **the mature-RNA measurement** — one-sided motif-stranded spliced mass, grafted into exons / peeled out | `NodeGeometry.spliced_*`, the graft/peel relay |
| B4 | **the hyperprior's density basis** — boundary mass over summed per-side eff-length | `node_global_geometry` |

They are also **half the chain**: `node_chain` interleaves `B0 R0 B1 R1 … Bk`, so removing them halves the
node count and changes every index in `bp_solver`.

### 4.2 The three options

**Option N — nodes own everything (coverage marginal), edges are pure flux.**
Every node gets `V_n`; edges carry only integer crossing counts used as flow/structure.
*Cost:* B1–B4 all need re-homing. B3 in particular has no home — the mature measurement is intrinsically a
*seam* quantity. And §2.3 says the coverage marginal is biased at exactly the seams. **Not recommended.**

**Option E — keep the bipartite structure; edges keep beliefs (status quo shape).**
Marginalize paths onto today's region/boundary objects; the solver is unchanged.
*Cost:* a path crossing several edges contributes to each, so the "consume a partition, never one of its
own marginals" rule (`02_redesign_derivation.md` §2.4) is violated — exactly as today. Measured scope:
**12.97 % of human exonic crossers cross more than one interface.**
*Benefit:* it is the only option that ships the index and accumulator without touching the solver, and it
is the **required dual-write gate** regardless (bit-identical 32/32).

**Option A ⭐ — anchored cells: the edge becomes a LABEL on an anchored partition.**
For node `A`, the paths anchored at `A` partition by where they end: `(A,A)` contained, `(A,A+1)`,
`(A, Z>A+1)`, `(A, …junction…, Z)`. This is a genuine partition of `n_A` (Claim 4 holds), each cell has an
exact effective length, and §3.1 of `02_redesign_derivation.md` shows the cells are **identified by FL
contrast** (Fisher cond 2.5 at gDNA 60 / RNA 200, best at short regions).
*Cost:* paths are left-anchored, so the two-sidedness of a boundary disappears — a real loss, since
`left`/`right` faces currently carry the flanking densities separately (B2). The mature measurement (B3)
becomes the count on junction-bearing cells, which is a natural home but a different one.
*Benefit:* no double counting, unbiased at steps, and short nodes stop being information-free.

### 4.3 Recommendation

⭐ **Amended.** Two findings since this section was drafted push the decision **earlier**, not later:
§2.4 (coverage re-creates the mass/flux duality for the strand likelihood) and
`03_path_accumulator.md` §5 (Option A needs `N_p` alone — a **4× smaller** accumulator schema, and the
schema must be frozen before the C++ is written). Both point the same way: **Option A**.

⭐ **And the node-count arithmetic favours it too.** Today's solver chain is `R + B` = 752 K + 753 K =
**1.5 M nodes**. Option A drops boundary nodes: **992 K** — *fewer than today*, on a 32 % finer partition.
Option E keeps them: 992 K + 992 K = **1.98 M**, a 32 % increase on a solve that is already 67 s at genome
scale. So Option E is not the cheap option at run time either.

**Ship Option E as the migration gate** (it is the only bit-identical arm and it is worth having), **but
take the Option A decision before the accumulator schema is frozen**, not after. Concretely:

1. Dual-write today's payload from paths → bit-identical 32/32. Nothing about the model changes.
2. Add the anchored-cell view alongside, unconsumed.
3. One A/B arm per role: B1 (the boundary's own strand belief) → the anchored cells; B4 (the hyperprior
   basis) → the anchored node density; then B2 and B3, which are the hard ones and should go last.

**Do not attempt Option N.** It puts the coverage attribution's blur exactly where §2.3 measured its worst
bias, and it deletes B3 with nowhere to put it.

---

## 5. ⚠⚠ BP ON A SPLICE GRAPH — the loopiness problem

The owner's plan is to traverse the graph with BFS so sources are processed before sinks. Two points, one
easy and one hard.

**The easy one: no BFS is needed.** Every edge satisfies `src < dst` in genomic order
(`00_splice_graph_design.md` I8), because a splice junction's donor always precedes its acceptor
genomically on both strands. **Genomic order IS a topological order.** A left-to-right sweep is already
correct; the existing `chain.order` needs no change in principle.

**The hard one: the graph is no longer a tree, and BP stops being exact.**

Today's chain is a linear path per reference — a tree — so one forward + one backward pass is **exact**.
That exactness is stated in `node_sweep`'s docstring and is load-bearing: it is why there is no outer
fixed-point loop and no convergence check.

Adding junction edges closes a loop for **every** junction. The cyclomatic number is
`E − V + C = 1.396 M − 0.992 M + ~200 ≈ 404,000` independent cycles on the human graph — one per junction,
as expected, because a junction plus the contiguous path it bypasses is a cycle. Consequences:

* BP becomes **loopy**: not exact, not guaranteed to converge, and it **double-counts evidence around
  cycles**, which manifests as over-confidence. That is the failure mode this project has fought twice
  already (`_far`'s 100 % structural data reuse; `_pin_v`'s self-confirming message), and the trust metric
  `z2` is the instrument that catches it.
* The governing principle — *"pass-0 must be WEAK and correctable"* — argues strongly against introducing
  404 K over-confidence loops to buy message paths.

**Options, in order of increasing risk:**

| | approach | exactness | change |
|---|---|---|---|
| **L1** ⭐ | **junction edges are FACTORS, not message channels.** Keep the contiguous chain as the BP tree; a junction's count is an observation attached to its two endpoint nodes. | **exact FB preserved** | small — this is close to what the solver does today with spliced mass |
| L2 | spanning tree = the contiguous chain; junctions handled as a correction pass after convergence | approximate, bounded | medium |
| L3 | full loopy BP with damping + convergence monitoring | not exact, may not converge | large |

**Recommend L1.** It gets TSS/TES, the path counts, the structural effective lengths and the terminus fix
without giving up exact inference. L3 is the destination the owner described, but it should be entered
**after** the mode defects are fixed and with `z2` as the gate — not as part of the same change.

⚠ **This is the single largest open modelling question in the rework**, and it is not answerable from a
design document. It needs a prototype against the real solver on real payloads.

---

## 6. FILE-BY-FILE IMPACT

| file | change |
|---|---|
| `substrate.py` | `CalibrationSubstrate` / `BoundarySubstrate` become projections of the path table; `_make_view`'s `(counts, mass)` pair collapses to one integer count |
| `node_chain.py` | unchanged under Option E; under Option A the boundary interleave goes away |
| `node_geometry.py` | `NodeGeometry`'s 20 per-face arrays shrink; `spliced_side_eff_length`'s continues/terminates reconstruction (~40 lines) is **replaced by the graph's edge flags** |
| `effective_length.py` | rewritten to the admissible-reach formulas; `boundary_side_eff_length`'s ½ (a property of the deposit rule) **deleted** |
| `node_init.py` | the four init sources keep their shape; the MEASURED source gains the FL-contrast cell likelihood (§2 of `02_redesign_derivation.md` §3) |
| `bp_solver.py` | `accept_l/accept_r` (observational) → the graph's `rna_crosses_contiguously` (structural, C3); `_rho_faces` per-component under the terminus typing (C1); `graft_premise_logvar` refit per structural class (C2) |
| `simplex_logodds.py` | unchanged |
| `priors.py`, `derive.py` | unchanged (they read the deconvolved masses) |
| `density_model.py`, `gdna_strand.py` | unchanged inputs, finer partition |

---

## 7. MEASUREMENTS

```
M-cal1  how much of the measured enrichment "slope" at a boundary is PHYSICAL
        capture vs the coverage attribution's blur?           -> §2.3, and it re-opens
                                                                 density_composition_reconciliation §4
M-cal2  anchored-cell vs coverage marginal on REAL cfRNA payloads: bias against the
        oracle, variance, and z2                              -> §4.3, the Option A/E decision
M-cal3  does the FL-contrast cell likelihood help AMBIG nodes on real data, where
        the strand likelihood is inert and od_r saturates?    -> the largest potential prize
M-cal4  loopiness: on the real graph, how much do junction messages change the
        answer, and does z2 degrade?                          -> §5, the L1/L3 decision
M-cal5  terminus taper: recompute omega_graft per structural class once the graph
        carries TSS/TES                                       -> P1g C2, already specified
```

⚠ All of these on real cfRNA. The synthetic suite has no alternative TSS/TES inside exons, no novel
junctions, an exon-region size 4.6× too generous, and is Poisson by construction.

---

## 8. WHAT IS SETTLED AND WHAT IS NOT

**Settled.** The `o_n/L` density identity and its `|n|` effective length (§2.1). The exact three-counter
marginalization (§3). The admissible-reach taper applying to both attributions (§3). Genomic order being a
topological order (§5). The measured bias/variance trade between attributions (§2.3).

**Not settled, and needing prototypes against the real solver, not more derivation.**
The mass-allocation option (§4.3). The loopy-BP question (§5).

### 8.1 ⚠⚠ THE BIGGEST RISK TO THE MOST EXCITING RESULT

The FL-contrast identifiability of §3.1 of `02_redesign_derivation.md` (Fisher cond 2.5, `f_g` bias +0.001)
is a **count-based** claim, and this project has been burned by count-based claims on real data three times
(P1e overstated 26–300× by the toy; `ω_graft` 15× apart across four real samples; the strand overdispersion
saturating its `Beta(2,2)` ceiling on 4/4 real libraries where the toy fits 0.0008–0.0017). Three specific
ways it could fail:

1. ⛔ **Overdispersion.** The design matrix is Poisson. Real counts are not. An overdispersed count fed to a
   Poisson deconvolution is **over-confident**, which is the direction the governing principle forbids
   (*"pass-0 must be WEAK and correctable"*). The synthetic suite is Poisson **by construction** and
   therefore **cannot** detect this. It must be tested on real cfRNA against an oracle, or not trusted.
2. ⛔ **Circularity.** The prototype used *exactly known* `f_g` and `f_r`. In production they are **fitted**,
   from the same library — and the gDNA FL pool is itself selected using region types. If the FL fit leans
   on the same fragments the deconvolution then explains, the identifiability is partly self-confirming.
   That is the `_pin_v` failure pattern in a new place. **Audit the FL fit's independence before trusting
   the contrast.**
3. **Regime.** The prototype's contrast was 3.3× (60 vs 200). Where the two FLs are close the design goes
   singular *smoothly* — so the term must degrade to inert on its own, with no gate and no constant. Verify
   that it does, rather than assuming it.

**Recommendation: treat the FL-contrast likelihood as a Phase-2 candidate, not part of the rework.** The
rework's value does not depend on it. Landing it inside the same change set would make an over-confidence
regression impossible to attribute.
