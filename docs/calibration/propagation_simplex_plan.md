# The simplex propagation deconvolution — execution plan

**Status:** execution spec, 2026-06-12. Redesigns the propagation solver after the flagship dissection
(`ff1a55a`) traced the catastrophic over-subtraction to its root: the count-cascade subtracts independently
estimated counts (`gDNA = U − mature − nascent`) that need not fit inside `U`, and the nascent estimate
(`RNA − mature`) is corrupted whenever `mature_density` misses a region's mature (single-exon / no-junction
— Step-0's 29% gap), so the mature is relabelled nascent and over-subtracted. **Two fixes, together:**

1. **Solve in the simplex (the pie), not by count subtraction.** Each node's total unspliced mass `U` is
   one pie cut into `(f₊, f₋, f_g)` summing to 1. Evidence *pushes* the point on the simplex; the answer is
   always a normalized partition. Over-subtraction is **structurally impossible**; an under-constrained
   slice stays **wide/uncertain** rather than forcing a wrong value elsewhere.
2. **Take the nascent from a mature-free source.** A boundary's **unspliced** crossing is physically
   `gDNA + nascent` (mature spans junctions only *spliced*), so the boundary nascent needs **no mature
   subtraction** — the corruption cannot occur.

Supersedes the count-cascade in `propagation.py` (`propagate`/`_solve_ambig`) and its count-floor band-aid.
Builds on the rev-2 likelihood model ([`phaseC_ambig_resolution_plan.md`](phaseC_ambig_resolution_plan.md)
§2) and the generative theory ([`deconvolution_generative_model.md`](deconvolution_generative_model.md)).

## 1. The principle — evidence as precision-weighted pushes on the pie

Precision lives in the **evidence** (counts → densities → likelihoods, carrying Poisson + propagation
uncertainty); the **pie** is the normalized result of combining them. Counts give precision; the simplex
gives safety. The per-node negative-log-posterior over `(f₊, f₋, f_g)` on the 2-simplex is a sum of
precision-weighted terms; the answer is its MAP (or posterior mean + width). Reparameterize by
`r = f₊+f₋ = 1−f_g` (RNA fraction) and `t = f₊−f₋` (tilt); `f_s = (r±t)/2`.

## 2. The five evidence terms (each carries its own precision)

| term | from | constrains | precision | form |
|---|---|---|---|---|
| **strand** | the node's own `U_pos/U_neg` | the **tilt** `t` (`p_pos = ½ + (κ−½)·t`) | `N·(2κ−1)²` (BB curvature) | BetaBinomial(`U_pos`\|`U`, `p_pos`, od) — sharp if κ≠½ & deep; flat at κ≈½ |
| **mature₊ / mature₋** | bounding **spliced** junction flux `S_s` (Phase-A `E_mu/E_ms`) | **lower bound** `f_s·U ≥ mature_s` | `1/S_s` (flux Poisson) | soft one-sided (RNA can't be below its mature); single-exon mature invisible ⇒ no push (fine) |
| **nascent₊ / nascent₋** | bounding **boundary unspliced** crossing, deconvolved (**mature-free**), propagated | **two-sided** estimate of the nascent part of `f_s` | crossing-flux Poisson **×** propagation-distance decay | Gaussian on `f_s` around `(mature_s+nascent_s)/U` |
| **gDNA lower bound** | the bounding boundary's deconvolved **gDNA** (capture-depleted) | **lower bound** `f_g·U ≥ gDNA_boundary` | boundary-flux Poisson | soft one-sided (subsumes `density_model`) |
| **prior** | weak Dirichlet | regularizes; keeps unconstrained slices sane | fixed weak | `Dir(½,½,½)`-like |

**Why this kills the bug class:** the mature & nascent now *push* `f_s` (bounded in [0,1], summing into the
pie) instead of being *subtracted* from `U`. A spurious "lots of nascent" is one soft push the strand and
the gDNA-lower-bound push back against — it can never zero the gDNA. And the nascent comes mature-free, so
the spurious value rarely arises in the first place. **Single-exon-transcript AMBIG regions** (no mature, no
nascent, no usable strand) are simply **left unspecified above the gDNA floor** — the honest "unknown" you
described, not a manufactured value.

## 3. The mature-free boundary node (the nascent source)

A **boundary** is itself a node, solved by the *same* simplex machinery with **no mature term** (its
unspliced crossing is `gDNA + nascent` by physics; its spliced crossing is the mature it *exports* to the
exon-side region, not part of its own pie). So a boundary's pie is `(nascent₊, nascent₋, gDNA)`:
- **clean-junction / single-strand boundary** → a **seed**: its own strand resolves `gDNA` vs `nascent_s`.
- **AMBIG boundary** → solved by propagation, like an AMBIG region.

A solved boundary **emits** per-bp densities `ρ_nasc₊, ρ_nasc₋, ρ_gDNA` (each with precision) to its two
flanking regions — the **local, enrichment-matched** nascent (term 3) and gDNA-lower-bound (term 4) the
region needs. The boundary's spliced flux feeds the exon-side region's **mature** (term 2).

## 4. The propagation (seeds → converge), all in the simplex

```
build the node graph (regions + boundaries, adjacency, observations, eff-lengths, strand_class).
SEED: intergenic regions (f_g=1); single-strand regions & clean-junction boundaries (own strand) -> emit
      (ρ_nasc±, ρ_gDNA) with precision.
WORKLIST: queue nodes whose neighbour emissions changed.
repeat (pop a node):
    gather neighbours' current (ρ_nasc±, ρ_gDNA) [precision-weighted, enrichment-matched]
    re-solve its pie (§2 MAP); if its emitted densities moved > tol, queue its neighbours.
until the worklist drains (seeds are fixed anchors -> termination; oscillation guarded by tol/optional damping).
read each region's f_g -> CalibrationResult; boundary pies -> flux transport (priors).
```

Precision **decays with propagation distance** (a carried nascent is less certain far from its seed), so a
confident node-strand always out-weighs a distant propagated nascent — your "strand governs" requirement,
emergent from the precisions, not a hand-set `w`.

## 5. The solver (tractable, elegant)

`p_pos` depends only on `t`, so **profile**: the strand gives a posterior over `t`; the mature/nascent give
`f₊,f₋` estimates (hence `t` and `r`); the gDNA-lower-bound caps `r`. Reduce to a 1-D MAP in `f_g`
(equivalently a small 2-D grid over `(r,t)` on the valid simplex), evaluated vectorized across nodes per
sweep. Reuses the BB posterior core. For a single-strand node (`f₋≡0`) it collapses to today's 1-D strand
posterior — the **no-regression guard** (verify byte-reduction when the nascent/mature terms are flat).

## 6. Code structure & consolidation

`propagation.py` rewritten around the simplex:
- `solve_node(obs, evidence, prior) -> (f₊,f₋,f_g, precision)` — the §2 MAP (one node; regions and
  boundaries share it, boundaries drop the mature term).
- `node_graph(...)` — regions + boundaries, adjacency, observations, eff-lengths, seed mask.
- `propagate(...)` — the §4 worklist loop → `(regions, left, right)` NodeDeconv.
- **Consolidates & retires**: `density_model` (→ the boundary gDNA term), `mature_density` kept (→ the
  mature term), the count-cascade `_solve_ambig` + the count-floor (→ deleted), `strand_deconv` per-node
  blend (→ the strand term), `run_fill` (→ the worklist). `splice_junction_eligibility` kept.
`derive`/`priors`/the EM hand-off unchanged.

## 7. OPEN ISSUES / UNCERTAINTIES (review together)

1. **The likelihood shapes.** The one-sided forms (mature & gDNA lower bounds) — clipped-Gaussian vs a
   soft-plus/log-barrier vs a Poisson tail — and the two-sided nascent (Gaussian on the fraction vs on the
   count). These are *shape* choices; **no magic constants** (every width is a derived flux precision), but
   the shape affects how hard a lower bound bites. Calibrate on the toy sweep. *Which shapes feel right?*
2. **Boundary→region geometry.** A region has two boundaries that may disagree (different enrichment /
   neighbours). Precision-weighted mean of the two? Nearest? And the boundary nascent is read at the
   *boundary's* enrichment — is that the region's? (For an exon region, both its boundaries are exon-edge —
   matched. For a region split mid-exon, both edges are exonic — matched. The mismatch case to check:
   an exon↔intron boundary feeding an exon region vs an intron region.)
3. **Single-exon mature invisibility.** The pie leaves `f_s` unspecified above the gDNA floor (honest).
   Accept, or attempt single-exon mature another way (hard — no junction; maybe the region's own contained
   spliced, or a length/coverage model)? *I lean: accept the "unknown," let the gDNA floor + the EM handle
   it.*
4. **Precision decay with distance.** The functional form of how a propagated nascent's precision falls
   with hops/bp from its seed. Inverse-variance accumulation along the chain? A fixed per-hop factor?
   Needs a principled choice (or derive from the run-fill variance).
5. **Convergence / damping / `tol`.** Does the worklist converge monotonically on long balanced-AMBIG runs?
   Seeds are fixed anchors (should terminate). Add damping only if a test oscillates. Define `tol`.
6. **The boundary's own gDNA for the transport.** It's capture-depleted (the exon edge). Use as-is (a
   lower bound), or reconcile with the flanking regions' residual? (Deferred decision from the prior plan.)
7. **Solver cost (perf, issue #4).** The per-node MAP × nodes × sweeps. Profile-1-D keeps it ~current cost;
   vectorize across nodes. Deferred to the perf phase but design for it.
8. **κ≈½ (unstranded).** The strand term goes flat (correct); the pie rests on mature (motif-based) +
   nascent (crossing-based) + gDNA — all strand-independent — with the nascent-vs-gDNA confound left wide.
   Confirm graceful degradation on the ss=0.5 conditions.
9. **No-regression + mass conservation.** Single-strand/simple/no-gDNA must reduce to today's strand
   result (golden); the simplex conserves by construction (Σ=1). Invariant tests.

## 8. Validation plan

1. **Per-node solve in isolation** on locus 21 (both nrna) + the toy AMBIG sweep (gDNA 0→high × nascent
   {none,both,one} × κ {0.99,0.5}): the pie recovers the oracle; **no over-subtraction** anywhere; settle
   the §7 shape/geometry choices here.
2. **Boundary nodes**: their mature-free pies recover the oracle nascent/gDNA at the boundaries; feed the
   regions correctly.
3. `test_ambig_scenario` + new simplex unit tests (single-strand reduces to strand; AMBIG resolves;
   mass conservation; the 4030271-type catastrophe is gone *without* the floor).
4. **Flagship end-to-end** (`--use-propagation`): the gross leak drops further than the +floor's −4%, and
   region 4030271 reaches ~0.94 (the cure, not the floor); single-strand conditions unchanged.
5. **Full gDNA suite**; flip the default; consolidate/retire modules; perf.

## 9. Migration

The current count-cascade `propagation.py` stays the opt-in path until the simplex version validates
end-to-end; then swap (still behind `use_propagation`), retire the count-floor, then flip the default.
The dissection tooling (`propagation_flagship_dissect.py`) is the regression bed.
