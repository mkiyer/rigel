# Phase C — the iterative propagation deconvolution: implementation plan

**Status:** dry-run plan, rev 3 (2026-06-12). Third phase of
[`deconvolution_roadmap.md`](deconvolution_roadmap.md), merging former Phase B. **Rev 3 generalizes the
per-node likelihood model (rev 2) into a system solve:** the per-node prototype proved a single region
cannot be solved in isolation (a symmetric AMBIG node is under-determined without its neighbours' RNA
level), so the deconvolution is a **graph of interdependent nodes solved by propagation from seeds to
convergence.** This subsumes Phase A's mature imputation, the gDNA density imputation (`density_model`), the
strand module, and the splice fraction into one solver. Builds on `mature_density.py` (`83b3610`) and the
Step-0 record ([`phaseB_mature_subtraction_plan.md`](phaseB_mature_subtraction_plan.md)).

## 1. The reframe — a system, not a per-node solve

The rev-2 prototype (`scripts/debug/phaseC_node_posterior.py`) showed: strand-tilt + mature-lower-bound +
gDNA-lower-bound **under-determine `f_g`** for a symmetric AMBIG node — the MAP collapses to the depleted
gDNA lower bound (region 224: 0.37 vs oracle 0.98). The missing information is the **per-strand RNA level**,
which lives in the node's **neighbours** (where a gene is single-strand and the strand resolves its RNA).
So:

1. a node's solution depends on its neighbours; neighbours depend on theirs; up to a **seed**;
2. solving a node **propagates** new information to its neighbours, which may then become solvable (or
   change);
3. **iteration to convergence is a requirement, not an enhancement** — it is how the seed information
   reaches the interior AMBIG runs.

**The key revelation (yours):** a **boundary**'s *unspliced* crossing counts also deconvolve into all three
`{RNA+, RNA−, gDNA}` — not just gDNA. So boundaries are full 3-term nodes too, and information flows both
ways: region→boundary→region.

## 2. The graph

Nodes = **regions** (contained view) **and boundaries** (the seam between adjacent regions). Per reference
they form a 1-D chain: `… region[i] — boundary[i|i+1] — region[i+1] …`. Edges connect each region to its two
bounding boundaries (and each boundary to its two flanking regions). Each node carries a 3-term partition of
its **unspliced** mass into `(f₊, f₋, f_g)` plus a deterministic **mature** term from its **spliced** counts.

| node | unspliced = | spliced = | seed? |
|---|---|---|---|
| intergenic region (`NONE`) | gDNA | — | **yes** → `f_g=1` |
| single-strand intron region | gDNA + nascent_s | — | **yes** (strand resolves) |
| single-strand exon region | gDNA + nascent_s + mature_s | mature_s (contained) | **yes** (strand resolves) |
| exon↔intron / exon↔intergenic **boundary** | gDNA + nascent (crossing) | mature_s (the junction) | **yes** (clean split) |
| **AMBIG** region/boundary | gDNA + nascent₊ + nascent₋ | mature₊ + mature₋ | **no** → needs neighbours |

## 3. What each node emits and consumes — the messages

A **solved** node emits per-bp **densities** (mass ÷ effective length):
`ρ_gdna`, `ρ_nasc₊`, `ρ_nasc₋` (gDNA + per-strand nascent). Mature is **local** (from the node's own /
bounding spliced crossings, Phase A geometry) — it is exon-specific, not carried.

- **nascent density is carried across the gene body** (introns + exons): it is the species whose density is
  ~uniform along a gene, so a single-strand neighbour's `ρ_nasc_s` is the right estimate for an AMBIG node
  on the same strand's gene. Carried along the chain by strand, **gated to the strand's own exon/intron run**
  (never across an unrelated gene).
- **gDNA density** is carried locally (adjacent nodes), but see §5 — for an AMBIG region the gDNA is taken
  as the **residual** of its own count, not the carried (depleted) density.

## 4. The per-node solve (the likelihood combination, rev 2, now with neighbour priors)

For a node with own counts `(U_pos, U_neg)`, own/bounding mature `(M₊, M₋)`, and neighbour-carried densities
`(ρ̄_nasc₊, ρ̄_nasc₋, ρ̄_gdna)`, find the MAP partition by combining:

- **own strand likelihood** `BetaBinom(U_pos | U, ½+(κ−½)·t, od)`, `t = f₊−f₋` — pins the **tilt** (sharp if
  κ≠½, flat if κ≈½). Reuses the existing BB posterior core.
- **mature** `M₊, M₋` — the spliced-derived contained mature, a component of `RNA_s = nascent_s + mature_s`.
- **neighbour nascent priors** `nascent_s ≈ ρ̄_nasc_s · eff_len` — the carried per-strand nascent level
  (two-sided, with the carry uncertainty). **This is what pins the flat direction.**
- **gDNA**: `f_g = 1 − f₊ − f₋` with `RNA_s = nascent_s + mature_s` — a *residual*, plus a soft consistency
  with `ρ̄_gdna` where available.

A **seed** has no (or weak) neighbour priors and is solved from its own strand alone (single-strand /
intergenic / clean-junction → the partition is determined). An **AMBIG** node consumes the neighbour nascent
priors to resolve the gDNA-vs-balanced-RNA confound, then refines the tilt with its own strand.

## 5. Why this fixes the capture depletion (the crucial property)

For an AMBIG region: `gDNA = U − RNA₊ − RNA₋`, where `RNA_s = nascent_s (carried, NOT capture-depleted) +
mature_s (junctions, NOT depleted)`. So the gDNA is the **residual of the region's own capture-enriched
count** minus non-depleted RNA estimates → it sits at the region's **true enriched level**, not the depleted
boundary level. This is exactly the depletion that defeated the boundary-gDNA-crossing approach (rev-2
prototype region 224: the carried nascent ≈ 0 ⇒ `RNA ≈ mature = 218` ⇒ `f_g ≈ 0.97`, recovering the oracle).

## 6. The iteration (the main loop)

```
build the node graph (regions + boundaries, adjacency, observed counts, eff-lengths, class).
solve all SEED nodes from own data; mark solved; emit densities.
repeat (forward then backward sweep along each reference chain):
    for each node:
        gather current densities from neighbours
        (re)solve the node's partition (§4); update its densities
until max |Δ(f₊,f₋,f_g)| across nodes < tol  (or max_iters reached).
read off each region's f_g → the calibration prior (derive/priors unchanged).
```

**Convergence is guaranteed** in the sense the user noted: every reference chain is anchored by seeds
(intergenic / single-strand / clean-junction nodes are always present and self-solvable), and an AMBIG run
between two seeds is filled by propagation from both ends — the bidirectional sweep is the existing
`run_fill` generalized to the 3-term solve. The `tol` is a numerical convergence tolerance (not a model
constant); a damping factor is available if a long balanced-AMBIG run oscillates (§9).

## 7. Implementation — what it builds, consolidates, retires

- **New `propagation.py`**: the node graph (regions+boundaries), the per-node solve (§4, reusing the BB
  posterior + `mature_density` geometry + the eff-lengths), and the seed→sweep→converge loop (§6).
- **Consolidates**: `density_model` (gDNA imputation → the gDNA message + residual), `mature_density`
  (Phase A → the local mature term), `strand_deconv`'s per-node blend (→ the §4 solve), `run_fill`
  (→ the sweep). The 3-term {gDNA, nascent₊, nascent₋} *is* the generative model (`deconvolution_*`).
- **Retires**: `region_splice_gdna_frac` (the splice fraction), the heuristic `w=I/(I+I₀)` blend, and the
  `ρ` weight from rev 2 — all subsumed by the likelihoods + the propagation.
- **Unchanged**: `derive`/`priors`/the EM hand-off (still read `f_g`·mass per region); the BB posterior
  core `strand_loglik`; `splice_junction_eligibility` (mature anchors).

## 8. Validation (prototype-first, the Phase-A/Step-0 discipline)

1. **Prototype the solver standalone** (`scripts/debug/phaseC_propagation.py`): build the graph for locus
   21, run the seed→sweep→converge loop, and check the AMBIG regions (226/231/236) recover to oracle and the
   single-strand seeds match their current values. **Before any production wiring.**
2. **Toy AMBIG sweep** — gDNA 0→high × nascent {none, both, one-strand} × κ {0.99, 0.5}: confirm seeds
   solve, propagation fills the AMBIG run, convergence in a few sweeps; calibrate the nascent-carry
   uncertainty here.
3. **`test_ambig_scenario`** — both cases pass.
4. **No regression** on single-strand / simple / no-gDNA / unstranded conditions (seeds reduce to today's
   strand result); golden refresh + review.
5. **Flagship + full suite** — the 8.7% leak drops (the AMBIG + the eff-len-coupled part both improve as the
   per-region gDNA gets the right enriched level); confirm the eff-len contraction follows.

## 9. Open questions / risks

- **Nascent-carry geometry** — which neighbour, same-strand same-gene, how far to carry, the carry
  uncertainty growing with distance. The toy sweep is the bed; gate to the strand's exon/intron run.
- **Boundary 3-term solve** — confirm the boundary's unspliced crossing semantics (the per-side flux) feed
  the same solve; the boundary's gDNA is the *depleted* one (it's the exon edge) — used as a message but the
  region's residual (§5) is the authority for the region's gDNA.
- **Convergence / damping** — does a long balanced-AMBIG run converge monotonically? Add damping only if a
  test oscillates; prefer none.
- **Performance** — iteration over the node graph; bounded sweeps per reference chain. Vectorize the
  per-node solve across nodes per sweep. Deferred to the perf phase, but design the sweep to be vectorizable.
- **Scope change**: rev 3 **lifts the "one feed-forward pass" lock** — iteration is now required, but it is
  bounded (a few bidirectional sweeps per reference chain to convergence), not an open-ended EM.

## 10. Phase D (unchanged)

Validate at scale: golden refresh; full gDNA suite + flagship; confirm the leak drop and no regression;
update the authoritative theory doc + CLAUDE.md.
