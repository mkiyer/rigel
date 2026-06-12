# The propagation deconvolution — implementation plan (`propagation.py`)

**Status:** implementation spec, 2026-06-12. The production design for the iterative propagation solver,
validated by two prototypes: the region-chain core (`scripts/debug/phaseC_propagation.py`, `1810dc1` —
seeds → propagate per-strand nascent → AMBIG resolves) and the boundary physics
(`scripts/debug/phaseC_boundary_propagation.py`, `244fd51` — unspliced bidirectional / spliced one-sided /
flux enrichment-local). Supersedes the rev-2 per-node cascade in
[`phaseC_ambig_resolution_plan.md`](phaseC_ambig_resolution_plan.md) and the roadmap §2/§6. Theory:
[`deconvolution_generative_model.md`](deconvolution_generative_model.md),
[`deconvolution_implementation.md`](deconvolution_implementation.md). **§11 lists the remaining issues — the
focus before/while coding.**

## 1. The object — a node graph, solved by propagation from seeds

Per reference, a 1-D chain alternating **region** and **boundary** nodes:
`r₀ — b₀ — r₁ — b₁ — … — r_{n−1}` (n regions, n−1 internal boundaries). Every node carries a 3-term
partition of its **unspliced** mass into `(f₊, f₋, f_g)` = (RNA+, RNA−, gDNA fractions, summing to 1) plus a
deterministic **mature** term from its spliced counts. The calibration reads each **region**'s `f_g`; the
boundary partitions feed the per-locus prior's flux transport (as today).

`strand_class` constrains the simplex (one solve for all node types): `NONE`→`f₊=f₋=0`; `POS`→`f₋=0`;
`NEG`→`f₊=0`; `AMBIG`→full. The solve is driven from **seed** nodes (self-solvable) and propagates to the
rest, iterating to convergence.

## 2. Node observations (from the substrate)

| node | unspliced | spliced | eff-lengths |
|---|---|---|---|
| **region** `r_i` | contained `n_unspliced_pos/neg` | contained `n_spliced_*` (mature, usually ~0 for single-exon-piece regions) | `E_g=region_eff_length(gdna)`, `E_rna=region_eff_length(rna)` |
| **boundary** `b_i` (between `r_i,r_{i+1}`) | crossing flux, **both sides** equal (`sub.right[i]` ≈ `sub.left[i+1]`) — bidirectional, mature-free | crossing flux, **one side** (the exon side) — mature, one strand | `fl_mean_gdna`, `fl_mean_rna` (crossing) |

Confirmed physics (`244fd51`): boundary unspliced is bidirectional (deposits both sides), boundary spliced
is one-sided (the junction read's body is on the exon side). A region's **mature** comes from its bounding
**boundaries'** spliced flux (the Phase-A `E_mu/E_ms` geometry: `mature_s = Σ_bounding S_s · E_mu/(n·fl_rna)`).

## 3. The per-node solve (the likelihood combination)

Reparam `r = f₊+f₋ = 1−f_g`, tilt `t = f₊−f₋`. The node's own strand sees only the tilt:
`p_pos = ½ + (κ−½)·t`. Inputs: own `(U_pos,U_neg)`; own mature `(M₊,M₋)`; **propagated** per-strand nascent
densities `(ρ̄_nasc₊, ρ̄_nasc₋)` and gDNA density `ρ̄_gdna` from neighbours (precision-weighted §6). Then:

```
RNA_s   = ρ̄_nasc_s · E_rna + M_s            # nascent (propagated, enrichment-matched) + mature (local)
f_s     = RNA_s / U
f_g     = 1 − f₊ − f₋                         # gDNA = the region's OWN enriched count minus RNA
```

refined by the strand likelihood `BetaBinom(U_pos | U, ½+(κ−½)t, od)` and a soft consistency with `ρ̄_gdna`.
- **seed** (single-strand / clean-junction / intergenic): no propagated nascent needed — the strand alone
  pins `f_g` (and emits the node's own `ρ_nasc_s = (RNA_s−M_s)/E_rna`). Reduces to today's strand module.
- **AMBIG**: the propagated `ρ̄_nasc_s` pins the otherwise-flat gDNA-vs-balanced-RNA direction; `f_g` is the
  residual; the strand tilt checks `f₊−f₋` for consistency.

Each solved node **emits**: `ρ_gdna = f_g·U/E_g`, `ρ_nasc₊`, `ρ_nasc₋` (the messages), each with a
**precision** (from its strand info `N·(2κ−1)²` and counts).

## 4. Seeds (guaranteed self-solvable — the propagation anchors)

| seed | unspliced is | how solved |
|---|---|---|
| intergenic region (`NONE`) | gDNA only | `f_g=1` (highest confidence) |
| single-strand intron region | gDNA + nascent_s | strand splits (no mature) |
| single-strand exon region | gDNA + nascent_s + mature_s | strand splits gDNA vs RNA_s; mature from junctions |
| clean-junction boundary (exon↔intron, exon↔intergenic) | gDNA + nascent_s (one strand) | strand splits (mature-free) |

`splice_junction_eligibility` classifies clean-junction boundaries. Every reference chain has seeds at its
intergenic/single-strand stretches; an AMBIG run is bracketed by seeds and filled from both ends.

## 5. Boundary solve (physics-specific)

- **Unspliced** (both sides) → deconvolve into `{gDNA, nascent₊, nascent₋}` by the same §3 solve (it is
  **mature-free**, so simpler — 2 components per strand). Emits `ρ_gdna`, `ρ_nasc_s` to **both** flanking
  regions (bidirectional).
- **Spliced** (one side) → `mature_s` on the exon-side strand → emitted to the **one** exon-side region
  (unidirectional), feeding that region's `M_s`.
- A clean-junction boundary is a seed (single strand). An AMBIG boundary (both strands cross) is solved by
  propagation, like an AMBIG region.

## 6. Propagation — precision-weighted, enrichment-matched

Messages flow region↔boundary↔region along the chain. A node's propagated `ρ̄_nasc_s` is the
**precision-weighted combination** of its neighbours' emitted `ρ_nasc_s` (weight = the emitter's precision,
`N·(2κ−1)²` — derived, replacing the prototype's count threshold). Two properties the prototypes proved
necessary:
- **enrichment-matched**: nascent is capture-enriched (high on exons, depleted in introns), so a region
  draws its nascent from its **bounding boundaries** (local, at its own enrichment) rather than a distant
  same-strand seed. The boundary nodes make this automatic — they sit at the region's edges.
- **precision-weighted**: tiny-mass nodes (sliver regions, low-flux boundaries) carry near-zero precision
  and cannot poison the propagation (the prototype's failure mode).

## 7. The solve loop

```
build the node graph (regions + boundaries, adjacency, observations, eff-lengths, strand_class).
solve all SEED nodes (§4) from own data; emit (ρ_gdna, ρ_nasc±) with precision.
repeat (forward then backward sweep per reference chain):
    for each node:
        gather neighbours' current emitted densities (precision-weighted, §6)
        (re)solve the node (§3 / §5); update its emitted densities + precision
until max |Δf_g| over all nodes < tol   (numerical tolerance; not a model constant)
read each region's f_g  →  CalibrationResult (mass split); boundaries → flux transport.
```

Convergence is anchored by the seeds (every chain has them) and the bidirectional sweep (= `run_fill`
generalized to the 3-term solve). §11 lists the convergence questions to settle.

## 8. Output & hand-off (unchanged downstream)

Per region: `mass_gdna_contained = f_g·U`, `mass_rna_contained = (1−f_g)·U + spliced`. Boundary partitions
→ `mass_gdna_left/right` for `priors._transport_boundary_flux`. `CalibrationResult` schema, `derive`,
`priors`, the EM hand-off — **all unchanged**. (The IPR eff-len in `priors` benefits automatically: a
correct per-region `f_g` un-amplifies the squared coupling the user identified.)

## 9. Consolidation map (issue #3 — code health falls out)

| current module | becomes |
|---|---|
| `density_model` (gDNA density + run-fill imputation) | the gDNA component of the node solve + the `ρ_gdna` message; run-fill → the propagation sweep |
| `mature_density` (Phase A) | the boundary→region mature term (`E_mu/E_ms` geometry kept) |
| `strand_deconv` per-node blend (`w=I/(I+I₀)`) | the §3 likelihood solve; `strand_posterior_gdna_frac` → the strand likelihood; **`w` blend retired** |
| `run_fill` | the propagation sweep |
| `splice_junction.region_splice_gdna_frac` (the fraction) | **retired** (subsumed) |
| `splice_junction_eligibility` | **kept** (clean-junction → seed boundaries + mature anchors) |
| `strand_loglik`, the BB posterior core | **kept** |
| the rev-2 `ρ` weight | **retired** (the propagated nascent does its job, principled) |

New module `propagation.py` (the graph + solve + loop). `calibrate.py` calls it in place of the
`node_gdna_density → cleaned → splice → deconv` block. `derive`/`priors`/`result` untouched.

## 10. Code organization (`propagation.py`)

```
NodeGraph.from_substrate(substrate, region_arrays, fl_models, config)
    -> nodes (regions+boundaries), adjacency, observations, eff-lengths, strand_class, seed mask
solve_node(node, neighbour_messages, kappa, ods, eff) -> (f_pos, f_neg, f_g), emitted densities + precision
propagate(graph, kappa, ods, ...) -> per-region (f_pos, f_neg, f_g)   # the seed→sweep→converge loop
   # vectorize the per-node solve across nodes per sweep (design for issue #4)
```

Keep single-purpose: the graph (structure), the solve (one node's likelihood), the loop (propagation).
Reuse the BB posterior + the eff-length functions. Aim for the single-strand path to **byte-reduce** to the
current strand result (the golden-stability guard).

## 10.5. The AMBIG solve, refined by the #1 resolution (validated `nrna_rnd`)

The user's #1 insight sharpens the AMBIG solve and is **validated** (`scripts/debug/phaseC_propagation_v2.py`,
`nrna_rnd` locus 21: AMBIG mean error **0.056**, region 224 the v1 failure 0.34→0.96):

> An AMBIG node has ≤2 strands; **one is always seedable** — a gene with single-strand stretches seeds
> directly, and a *nested* gene's TSS/TES boundary exposes the OTHER gene's nascent as a 2-component
> (strand-solvable) node. So **propagate the seedable strand's nascent, subtract it + both matures, and the
> residual = `gDNA + the nested strand's nascent` → the node's OWN strand resolves it** (gDNA symmetric, the
> remaining nascent tilted). Propagation supplies one strand; the node's strand finishes the other.

This means the AMBIG solve **uses the node's own strand** (not `w=0`) after subtracting the propagated
other-strand RNA — the strand does the final gDNA-vs-(one-strand-)nascent split it is good at. **Enrichment
matters** (validated): the propagated nascent must be enrichment-class-matched — intron regions get the
off-target (≈0, capture-depleted) nascent, exon regions the on-target nascent; carrying an exon seed's
nascent to an intron region was the v1 failure.

## 11. REMAINING ISSUES — the focus for next

1. **The gDNA-vs-balanced-nascent confound — RESOLVED (§10.5).** The user's insight (≤2 strands, one always
   seedable, node-strand resolves the residual) dissolves it; validated on `nrna_rnd` (0.056). The genuine
   residual is now narrower: a node where *neither* strand is seedable (a pathological all-AMBIG reference,
   or both strands fully nested) → the global gDNA prior is the deferred ultimate fallback (issue #6).
2. **Precision / uncertainty model — the now-primary issue.** The simple plan (issue #2 decision): the
   strand weight `w = N·(2κ−1)²` governs; count/propagation gets `(1−w)`; a confident node strand is not
   overridden by an uncertain propagated nascent. The prototype's nrna_none regression (0.042→0.117) is
   exactly this: a no-nascent exon region's confident own-strand says "symmetric ⇒ gDNA," but a small
   *spurious* seed-nascent (from the seed's strand-fit error) was over-subtracted. The fix is to weight the
   propagated nascent by its precision (the seed's strand info, decayed with distance) so a confident node
   strand wins. Derive it; no magic number. **This is the main thing to get right.**
3. **Enrichment-matching via boundaries — confirm it's automatic.** Does drawing nascent from bounding
   boundaries fully replace explicit exon/intron class gating, or are there transitions (e.g. an AMBIG run
   spanning an exon→intron flip) where a region still pulls cross-enrichment nascent? Test on the toy sweep.
4. **AMBIG boundary deconvolution + the depleted boundary gDNA.** A boundary's own gDNA crossing is
   capture-depleted (it's the exon edge). For the region we use the residual (un-depleted); for the
   **boundary**'s own `f_g` (which feeds the flux transport) we must decide: depleted crossing, or
   reconciled with the flanking regions' residual? Get the boundary gDNA mass right for the transport.
5. **Convergence & damping.** Does a long balanced-AMBIG run converge monotonically? Need: a convergence
   proof/argument or an empirical bound on sweeps; a damping factor only if oscillation appears; a
   max-iters backstop. Define `tol`.
6. **Seed sufficiency / ordering.** Confirm every reference chain is seed-anchored; handle a chain with no
   seed (all-AMBIG reference — pathological but define behaviour). Sweep order (forward/backward) vs a
   priority queue (solve highest-precision first).
7. **Mass conservation** across the graph (Σ gDNA + Σ RNA = Σ fragments) — an invariant test.
8. **Single-strand / simple-condition no-regression** — the solve must reduce to the current strand result
   there (golden byte-stability where expected); the validation guard for the big rewrite.
9. **Performance (issue #4).** Per-node solve × nodes × sweeps. Vectorize the solve across nodes; bound the
   sweeps. Profile on the real dataset once correct.

## 12. Validation plan (prototype-first, then production)

1. Extend the prototypes into a full graph solve on locus 21 + the toy AMBIG sweep (gDNA 0→high × nascent
   {none, both, one-strand} × κ {0.99, 0.5}); settle issues #1–#5 there.
2. `test_ambig_scenario` + new propagation unit tests (seed solve, boundary solve, convergence, mass
   conservation).
3. Golden refresh + review; **no regression** on single-strand/simple/unstranded/no-gDNA conditions.
4. Full gDNA suite + flagship: the 8.7% leak drops; the eff-len contraction follows; confirm per-condition.
5. Then issue #3 (cleanup is largely done by the consolidation) and issue #4 (performance).
