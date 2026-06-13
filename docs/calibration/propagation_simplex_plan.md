# The simplex propagation deconvolution — execution plan

**Status:** execution spec, rev-2 2026-06-13. Redesigns the propagation solver after the flagship
dissection (`ff1a55a`) and a design review (`docs/calibration/my_notes.md`) that **simplified the model**:

> **Calibration models only RNA vs gDNA — not mature vs nascent.** The goal of calibration is to solve
> for the gDNA fraction at each node; it does **not** need to separate nascent from mature RNA — the
> per-locus EM already does that well downstream. So every node's unspliced mass is one pie cut into
> three slices `(f_rna₊, f_rna₋, f_g)` summing to 1. There is **no nascent slice and no mature slice**.

This deletes the layer that produced the catastrophic bug. The old count-cascade estimated
`gDNA = U − mature − nascent` and was corrupted whenever `mature_density` missed a region's mature
(single-exon / no-junction); the mature was relabelled nascent and over-subtracted to `f_g=0`. Under the
simplified model **we never subtract mature** — mature stops being a tracked quantity and becomes just
*one kind of evidence* (spliced fragments = a sided lower bound on RNA). Over-subtraction is structurally
impossible and the corruption cannot arise.

Two principles, together:

1. **Solve in the simplex (the pie), not by count subtraction.** Evidence *pushes* the point on the
   simplex; the answer is always a normalized partition. An under-constrained slice stays **wide/uncertain**
   rather than forcing a wrong value elsewhere. Part of the pie may be **unknown** (not zero, not guessed).
2. **RNA must be earned; gDNA is the default.** The RNA prior is **zero** (RNA abundance is unpredictable —
   neighbouring isoforms vary >10000×, so we cannot assume any baseline). gDNA has a **weak nonzero prior at
   every node**. With no evidence, a node is gDNA; RNA only appears where evidence (strand, spliced flux, or
   propagated count) reduces RNA's uncertainty enough to overcome the gDNA prior.

Supersedes the count-cascade in `propagation.py` (`propagate`/`_solve_ambig`) and its count-floor band-aid.

## 1. The pie and its reparameterization

Per node, total unspliced mass `U`, pie `(f_rna₊, f_rna₋, f_g)`, `Σ=1`. Reparameterize:
`r = f_rna₊ + f_rna₋` (total RNA fraction = `1 − f_g`) and `t = f_rna₊ − f_rna₋` (strand tilt);
`f_rna_s = (r ± t)/2`. The strand sees **only** `t` (`p_pos = ½ + (κ−½)·t`); the count sees **only** `r`
(RNA-vs-gDNA magnitude); they are orthogonal, which is what lets them combine cleanly.

## 2. Initialization from the node signature

The signature → 3-term **fraction** and **uncertainty** vectors. An absent strand is a fact (its slice is
exactly 0, fully certain); an active strand is unknown (∞ uncertainty); gDNA always starts at its weak prior.

| node signature | `(f_rna₊, f_rna₋, f_g)` init | `(unc₊, unc₋, unc_g)` init |
|---|---|---|
| intergenic / intergenic↔exon (NONE) | `(0, 0, 1)` | `(0, 0, 0)` — fully specified |
| single-strand POS (e.g. `in+`, `E+`) | `(nan, 0, prior)` | `(∞, 0, prior_unc)` |
| single-strand NEG | `(0, nan, prior)` | `(0, ∞, prior_unc)` |
| AMBIG (both strands) | `(nan, nan, prior)` | `(∞, ∞, prior_unc)` |

`prior` = the weak gDNA prior (§5). A NONE node is a **seed** with zero work. A single-strand node with
usable strand data is a **seed** (strand deconvolution solves it ab initio). AMBIG and thin/unstranded nodes
start unsolved and are filled by propagation — or left honestly "unknown" above the gDNA floor.

## 3. The evidence terms (each carries its own precision)

| term | from | constrains | precision | form |
|---|---|---|---|---|
| **strand** | the node's own `U_pos/U_neg` | the **tilt** `t` | `N·(2κ−1)²` (BB curvature) | BetaBinomial(`U_pos`\|`U`, ½+(κ−½)t, od) — sharp if κ≠½ & deep; flat at κ≈½ |
| **spliced RNA (sided)** | bounding **spliced** junction flux on the **exon side** | **lower bound** `f_rna_s·U ≥ spliced_s` | `1/S_s` (flux Poisson) | soft one-sided; single-exon ⇒ no junction ⇒ no push (fine); spliced is SIDED → feeds only the exon-side region |
| **count (propagated)** | the **unspliced** count field imputed from neighbours across the boundary (bidirectional) | the RNA-vs-gDNA magnitude `r` (equivalently `f_g`) | NB / loess `var~mean` (Poisson floor to start, §4) | gDNA-density imputation; ignores the sided spliced; decays with propagation distance |
| **gDNA prior** | weak, nonzero, derived (§5) | keeps `f_g` from collapsing with no evidence; graceful at κ≈½ | fixed weak | `Beta(½,½)`-like; ~0.5 count strength |
| **RNA prior** | — | **zero** (RNA earned, never assumed) | — | no push |

**Why this is safe:** strand and spliced flux *push* the RNA slices up (bounded in the pie); the count and
gDNA prior anchor `f_g`. No term subtracts a count from `U`, so nothing can drive a slice negative. A node
with no evidence rests at the gDNA prior; the EM still allocates RNA downstream (a zero calibration prior on
RNA means "no evidence," not "forbidden").

## 4. Modelling uncertainty (the one real open build)

Precision is how the terms compete fairly for the pie. Three sources, three models:

- **Strand** — `N·(2κ−1)²` (BB curvature). **Have it**, validated.
- **Count** — Poisson **overestimates** precision (a thin AMBIG node's imputation would wrongly overpower a
  weak strand). The fix: a nonparametric **`variance ~ mean`** model. Observations = the imputation residual
  between a region and its two anchoring boundaries — extend to **(left-boundary, region, right-boundary)
  triplets** for richer residuals. **Build plan:** ship with a Poisson floor (proves the frame), then swap
  in the loess/NB precision as a clean drop-in. **Not yet built.**
- **gDNA** — a `variance ~ mean` model is wanted but harder: under capture only a subset of nodes is fully
  specified at init, and fitting needs deconvolution first → compounded uncertainty. Second step, after count.

## 5. The gDNA prior (the one sanctioned starting constant)

RNA prior = 0; gDNA prior = weak nonzero, **eventually derived** (straightforward for non-capture, nontrivial
under capture). **Start:** a `Beta(½,½)`-like Jeffreys prior with ~0.5-count strength as a tiebreaker.
Per the no-magic-numbers rule this is logged as a **placeholder pending derivation**, not a silent heuristic;
deriving the per-node gDNA expectation (the count-density global) is tracked as follow-up.

## 6. The propagation (seeds → converge), all in the simplex

```
build the node graph (regions + boundaries, adjacency, observations, eff-lengths, signature).
INIT every node from its signature (§2): absent strands -> 0/0; active -> nan/inf; gDNA -> weak prior.
SEED (ab-initio-solvable = strand/signature resolves it): NONE nodes (f_g=1); single-strand nodes whose
     strand deconvolution reduces RNA uncertainty below the gDNA prior. Seeds emit (count/gDNA density) + precision.
WORKLIST: queue nodes a seed/changed-neighbour can inform.
repeat (pop a node):
    gather neighbour emissions [inverse-variance accumulated -> order-independent], disregard SIDED spliced.
    re-solve its pie (§3 MAP). Track the IMPROVEMENT in uncertainty.
    propagate the uncertainty-improvement to neighbours, DAMPED on each merge
      (UNLESS the neighbour's prior uncertainty was inf -> the signal passes intact).
    a signal that no longer reduces uncertainty DIES (this is the convergence guarantee).
until the worklist drains (seeds are fixed anchors; tol on the uncertainty improvement).
read each region's f_g -> CalibrationResult; boundary pies -> flux transport (priors, §8).
```

Two design rules make this elegant and correct:

- **Order-independence:** merge evidence by **inverse-variance accumulation** (accumulate Σprecision and the
  precision-weighted mean), not sequential blending — so the worklist fixed point does not depend on visit
  order. This is the "ideally non-deterministic" property: the system of equations has one solution.
- **Solvability is non-binary:** "solvable" = "is there information that reduces uncertainty." A high-precision
  node (deep, strand-specific) barely moves under propagation — your "strand governs" requirement, emergent
  from the precisions rather than a hand-set blend weight `w`.

## 7. The solver (tractable)

`p_pos` depends only on `t`; the count/gDNA depend only on `r`. **Profile:** strand → posterior over `t`;
spliced lower-bounds + count → `r`; reduce to a 1-D MAP in `f_g` (or a small 2-D grid over `(r,t)` on the
valid simplex), vectorized across nodes per sweep. A single-strand node (`f_rna₋≡0`) collapses to today's
1-D strand posterior — the **no-regression guard**. Boundaries use the **same** solver (their pie is also
`(f_rna₊, f_rna₋, f_g)`; their spliced flux is exported to the exon-side region, not part of their own pie).

## 8. Boundary transport (deferred fold-in)

The existing boundary-flux transport (`priors._transport_boundary_flux`) — a post-deconvolution mass shift
that helps hybrid capture and is a no-op off-capture — **stays as-is** for now. Whether it can become part of
signal propagation (rather than a final step) is an open design question to revisit once the simplex solver
is validated. No change this phase.

## 9. Code structure & consolidation

`propagation.py` rewritten around the simplex:
- `solve_node(obs, evidence, prior) -> (f_rna₊, f_rna₋, f_g, precision)` — the §3 MAP (one node; regions and
  boundaries share it).
- `node_graph(...)` — regions + boundaries, adjacency, observations, eff-lengths, signature → init vectors.
- `propagate(...)` — the §6 worklist loop → `(regions, left, right)` NodeDeconv.
- **Consolidates & retires:** the count-cascade `_solve_ambig` + the count-floor (→ deleted); the per-node
  strand-blend `w` in `strand_deconv` (→ the strand term + emergent precision); `run_fill` (→ the worklist);
  `density_model` becomes the **count-evidence / imputation** provider. **`mature_density` is removed as a
  pie quantity** — its only surviving role (spliced flux as a sided RNA lower bound) reads directly from the
  substrate's spliced channels, so the separate per-strand mature-density imputation is no longer needed.
  `splice_junction_eligibility` kept. `derive`/`priors`/the EM hand-off unchanged.

## 10. OPEN ISSUES / UNCERTAINTIES (post-review status)

Resolved in review (`my_notes.md`): the model simplification (RNA vs gDNA only, no mature/nascent slices);
RNA prior = 0 / gDNA weak prior; init-from-signature vectors; single-exon = honest "unknown"; the
uncertainty-improvement propagation + damping rule; solvability = uncertainty-reduction; κ≈½ via the gDNA
prior; no-regression + mass conservation. **Remaining:**

1. **The count uncertainty model (§4).** The real build: loess/NB `variance~mean` on region↔boundary (triplet)
   residuals, replacing the Poisson floor. The gDNA `var~mean` (capture-compounded) is a second step.
2. **The one-sided shape** for the spliced RNA lower bound and the gDNA prior — clipped-Gaussian vs soft-plus
   vs Poisson-tail. Shape only (no magic constants); calibrate on the toy sweep.
3. **Boundary→region geometry.** Two boundaries may disagree → inverse-variance mean (§6). Confirm the
   enrichment match on the exon↔intron-boundary-feeding-an-exon-vs-intron-region case.
4. **The gDNA-prior derivation (§5).** Replace the ~0.5-count placeholder with the derived per-node gDNA
   expectation (the count-density global), capture-aware.
5. **Boundary transport fold-in (§8).** Keep post-hoc for now; revisit.
6. **Solver cost / parallelism.** Profile-1-D, vectorized across nodes; design for the parallel
   system-of-equations the user envisions. Deferred to the perf phase.

## 11. Validation plan

1. **Per-node solve in isolation** on locus 21 (both nrna) + the toy AMBIG sweep (gDNA 0→high × nascent
   {none,both,one} × κ {0.99,0.5}): the pie recovers the oracle; **no over-subtraction**; settle the §10
   shape/geometry choices here.
2. **Boundary nodes** recover the oracle RNA/gDNA at the boundaries; feed the regions correctly.
3. `test_ambig_scenario` + new simplex unit tests (single-strand reduces to strand; AMBIG resolves; mass
   conservation; the 4030271-type catastrophe is gone *without* any floor).
4. **Flagship end-to-end** (`--use-propagation`): the gross leak drops below the +floor's −4%, region 4030271
   reaches ~0.94 (the cure, not the floor); single-strand conditions unchanged.
5. **Full gDNA suite**; flip the default; consolidate/retire modules; perf.

## 12. Migration

The current count-cascade `propagation.py` stays the opt-in path until the simplex version validates
end-to-end; then swap (still behind `use_propagation`), retire the count-floor, then flip the default.
The dissection tooling (`propagation_flagship_dissect.py`) is the regression bed.
