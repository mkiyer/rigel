# Propagation = message passing on a chain — theory, problem statement, roadmap, increment-3 dry-run

**Status:** theory + execution doc, 2026-06-13. Companion to
[`propagation_simplex_plan.md`](propagation_simplex_plan.md) (rev-2: calibration models RNA-vs-gDNA only;
each node is a pie `(f_rna₊, f_rna₋, f_g)`). This doc answers the design questions in
`my_notes.md` about **how signals propagate between nodes**, derives and verifies the theory, reconciles
the two mental models (loglik-residual vs inverse-variance), states the problem cleanly for literature
search, and gives the phased roadmap + the increment-3 dry-run implementation plan with its open issues.

---

## 1. The headline result

> **A locus is a 1-D chain of nodes (a tree). On a tree, belief propagation is EXACT and converges in
> exactly two sweeps — forward (left→right) then backward (right→left) — independent of the schedule.**

This is precisely the algorithm `my_notes.md` derived from first principles ("iterate left to right, then
right to left … every node talks to every other node … the answer must not depend on the direction").
It is the **forward–backward algorithm** (= the Kalman/RTS smoother when messages are Gaussian; = the
sum–product algorithm on a tree). The user's instinct is correct and it is a theorem, not a heuristic:
two sweeps suffice, and the fixed point is unique (order-independent).

So **yes** to every question posed:
- *All nodes emit signals?* Yes — every node sends a message to each neighbour.
- *L→R then R→L sufficient for solvency?* Yes — exactly sufficient and exact on a chain.
- *Order-independent?* Yes — tree BP has a unique fixed point regardless of message schedule.

---

## 2. The graphical model (what the nodes and signals actually are)

Latent state per node `i`: the pie `θ_i = (f_rna₊, f_rna₋, f_g)_i` on the 2-simplex. Two kinds of factor:

**(a) Local evidence `ψ_i(θ_i)`** — the node's OWN observations, a log-likelihood over the simplex:
- the strand term (the 3-component mixture loglik, `solve_node`) on its `(u₊, u₋)`;
- the sided spliced flux (the RNA lower bound);
- the weak gDNA prior.
This is exactly `solve_node` **minus the count term**.

**(b) Pairwise coupling `φ_{ij}(θ_i, θ_j)`** — what *links* neighbours. The key modelling commitment:

> **Only gDNA propagates. RNA does not.** gDNA contamination density is a smooth property of the library
> (it varies slowly along the genome); RNA abundance is local and unpredictable (neighbouring isoforms
> differ >10000×, `my_notes.md`). So the coupling ties the **gDNA density** `ρ_g = f_g·U/eff_len` of
> adjacent nodes (`ρ_g,i ≈ ρ_g,j`) with a tolerance = the **count `var~mean`** imputation variance
> (increment 4). The RNA slices are **uncoupled** — each node earns its RNA locally.

This is the precise meaning of "the signal that propagates": **the gDNA-density estimate**. The "count
term" in `solve_node` (the Gaussian pull of `f_g` toward `count_gdna_frac`) is therefore **not local
evidence — it is the incoming message** from the neighbours via this coupling. Local = strand + spliced +
prior; propagated = gDNA density. That clean split is the whole architecture.

**The graph is a chain.** The calibration partition is sorted in genomic order; nodes (regions and
boundaries, interleaved) are adjacent only to their genomic neighbours within a reference. Even an AMBIG
node (overlapping opposite-strand transcripts) is a *single* node on the chain carrying both strands
internally — it creates no loop. A chain is a tree ⇒ §1 applies.

**Geometry → factors → abstract message passing.** `my_notes.md`: "abstract away transcriptome geometry
… replace it with 3-term tuples of log-likelihoods." **Yes, with one nuance:** geometry (signature,
eff-lengths, spliced flux, adjacency, capture enrichment) is *compiled* into the factors `ψ_i` and
`φ_{ij}`; once compiled, inference is pure message passing over loglik functions on the simplex — fully
geometry-agnostic. The "3-term tuple" is the message representation (a loglik over `(f₊,f₋,f_g)`, or its
sufficient statistics).

---

## 3. The message algebra — and reconciling the two mental models

**Messages add in log space.** A message `m_{i→j}(θ_j)` is a loglik over `θ_j` summarising "everything on
`i`'s side." The two operations:

```
message i→j :  m_{i→j}(θ_j) = log ∫ exp[ ψ_i(θ_i) + m_{left→i}(θ_i) + φ_{ij}(θ_i,θ_j) ] dθ_i
belief at i :  b_i(θ_i)     = ψ_i(θ_i) + m_{left→i}(θ_i) + m_{right→i}(θ_i)   (+ normalize)
```

The outgoing message to `j` includes everything `i` knows **except what `j` told `i`** (`m_{left→i}` only,
when sending right). Equivalently, from the belief: `m_{i→j} = b_i − m_{j→i}`.

**This IS the user's "residual," made precise.** `my_notes.md` proposed `propagate = input_loglik −
current_node_loglik` and "only propagate the productive aspect." The correct subtracted quantity is **the
message previously received from the target** (`m_{j→i}`), not the node's own loglik. That subtraction is
exactly what prevents double-counting and is what makes the result order-independent. So the instinct
("subtract to avoid sending back what improves nothing") is right; the precise bookkeeping is the BP
message-exclusion rule.

**"If it doesn't productively change the solution, stop."** Correct as a *convergence* rule. On a tree,
after the backward sweep no message changes any belief — the system has settled, and we stop (two sweeps,
automatic). The "productive change" threshold is the general stopping rule (and the one you'd need on a
loopy graph); on our chain it is satisfied exactly after two sweeps.

**"A signal dies into a seed."** Right, with one correction. A seed (intergenic `f_g=1`; or a deep,
strand-resolved node) has a **sharp** local evidence `ψ_i` (high precision), **not** `loglik = 0`. Its
outgoing message is dominated by `ψ_i`; an incoming message from the far side barely perturbs it, so the
far-side "news" cannot pass through — the seed is an information **sink/source** that resets propagation.
What goes to zero is the *uncertainty* of the seed (posterior variance → 0) and the *information gain* a
distant message provides past it — not the loglik. (Conflating "uncertainty → 0" with "loglik → 0" is the
one place the mental model needs adjusting.)

### 3.1 Likelihood vs inverse-variance — the same thing at two resolutions

This is the stated theoretical gap. The resolution:

- **Log-likelihood functions are the general, exact currency.** Combining independent evidence = adding
  logliks. The grid `solve_node` already does this exactly (it evaluates the full loglik on a simplex
  lattice). No approximation.
- **Inverse-variance accumulation is the Gaussian special case.** If each message is *approximated* by a
  Gaussian in the simplex coordinates `(r, t)` — keep only a mean + precision — then "add logliks" becomes
  "add precisions, precision-weight the means." That is inverse-variance weighting. It is **not a rival
  theory**; it is a cheap parametric shadow of the loglik message (Gaussian belief propagation / a
  moment-matched EP site / a Laplace approximation of each message).
- **Precision decays with distance because the coupling adds variance.** Each hop through `φ_{ij}` convolves
  the message with the coupling's imputation variance `σ²_couple` (the count `var~mean`), so a propagated
  gDNA-density precision falls as `1/(1/prec₀ + Σ_hops σ²_couple)`. The decay is **derived** from the count
  model, not a hand-set per-hop factor. (Answers simplex-plan §10.4.)

**Practical consequence — you don't have to choose.** Use the **grid message** (reuse increment-1's
lattice) for a first correct, constraint-respecting, inequality-factor-exact implementation; swap to
**Gaussian/inverse-variance** messages later only if perf demands (it's a clean drop-in with the same
algebra). The grid is the exact loglik; inverse-variance is its O(1) approximation.

---

## 4. The problem statement (self-contained, for literature search)

> **Inference target.** A genome is partitioned into nodes `i = 1..n` in linear (genomic) order. Each node
> has a latent composition `θ_i = (f₊, f₋, f_g)_i` on the 2-simplex (RNA-plus / RNA-minus / gDNA fractions
> of the node's `Uᵢ` observed unspliced fragments). Observations per node: a strand split `(u₊, u₋)_i`
> (Beta-Binomial mixture: gDNA at sense ½, RNA₊ at κ, RNA₋ at 1−κ, fitted overdispersions), a sided
> spliced-junction flux (a lower bound on the corresponding RNA slice), and a weak gDNA prior.
>
> **Coupling.** Adjacent nodes are tied only through the **gDNA density** `ρ_g,i = f_g,i·Uᵢ/Lᵢ`
> (effective length `Lᵢ`): `ρ_g,i ≈ ρ_g,i+1` with variance from an empirical `variance ~ mean` model. RNA
> slices are uncoupled. Under hybrid capture the gDNA density is **discontinuous at exon edges** (capture
> enriches exons), so the coupling must be enrichment-aware (capture-aware mean shift across an
> exon↔intron interface).
>
> **Graph.** Linear chain (a tree); seeds = intergenic nodes (`f_g=1`) and informative single-strand
> nodes. Goal: each node's posterior marginal over `θ_i` (especially `f_g` and its uncertainty).

**This is:** Bayesian inference on a **linear-chain factor graph / 1-D Markov random field** with
**compositional (simplex-constrained) latent states**, **one-sided (inequality) likelihood factors**, and
a **smoothness coupling on a derived quantity** (density, not the latent fraction directly). Exact
inference = the **forward–backward / sum–product algorithm**; tractable approximations = **Gaussian belief
propagation**, the **Kalman/RTS smoother**, **expectation propagation** (for the non-Gaussian one-sided and
simplex factors), or a **grid/discretized message** (exact up to lattice resolution).

**Search terms:** *sum-product algorithm on a chain*, *forward-backward algorithm*, *belief propagation
tree exactness*, *Gaussian belief propagation*, *Kalman / RTS smoother*, *expectation propagation factor
graph*, *message passing with truncated / inequality factors*, *state-space model compositional /
Dirichlet state*, *hidden Markov model continuous state smoothing*.

---

## 5. Phased roadmap (the whole line of work)

| phase | deliverable | status |
|---|---|---|
| **0** | rev-2 simplex plan; this theory doc | ✅ `d3bd8fe`, this commit |
| **1** | `solve_node` — per-node simplex MAP (local evidence + count + prior), isolated tests | ✅ `12a287b` |
| **2** | `init_from_signature` + chain adjacency (the graph skeleton) | ✅ `2494d2e` |
| **3** | **`propagate`** — build the node chain, compile local factors `ψ_i` + coupling `φ_{ij}`, run forward→backward (grid messages), emit per-region/per-side `f_g`. Replaces the count-cascade `_solve_ambig`. | ⏳ dry-run §6 below |
| **4** | the **count `var~mean`** coupling precision (loess/NB on region↔boundary triplets), capture-aware coupling mean shift — replaces the Poisson floor | ⏳ |
| **5** | flagship end-to-end behind `--use-propagation`; 4030271 → ~0.94 *without* a floor; single-strand unchanged | ⏳ |
| **6** | full gDNA suite; flip the default; consolidate/retire (`density_model`→coupling, count-cascade `_solve_ambig`+floor deleted, `run_fill`→the sweep, `strand_deconv` blend→local factor); revisit boundary transport | ⏳ |
| **7** | perf — Gaussian-message swap if needed; parallelize across independent chains (loci/references) | ⏳ |

**Are we closer?** Substantially. Phase 0–2 turned a heuristic-with-band-aids into a **well-posed graphical
model with an exact, order-independent inference algorithm** that matches the user's own derivation. The
remaining unknowns are now narrow and concrete: the coupling variance (phase 4) and the capture-aware
coupling mean shift (the one genuine modelling risk, §6 issue A). The architecture is settled.

---

## 6. Increment-3 dry-run — `propagate` implementation plan + open issues

**Goal.** Replace the opt-in count-cascade with the forward–backward message passing, keeping the
`(regions, left, right)` `NodeDeconv` output interface so `calibrate` is unchanged.

**Sketch.**
```
1. linearize the chain: per reference run, interleave [region r, boundary(r,r+1), region r+1, ...].
2. compile local factors ψ_i (grid loglik over the lattice):  strand (3-comp mixture) + sided spliced
   lower bound + gDNA prior.  (region nodes: full;  boundary nodes: strand + prior, spliced is the
   exon-side region's.)   Seeds: intergenic -> delta at f_g=1.
3. compile coupling φ_{ij}: gDNA-density continuity ρ_g,i≈ρ_g,j, variance σ²_couple (Poisson floor now,
   count var~mean in phase 4), capture-aware mean shift across exon↔intron interfaces.
4. FORWARD sweep i=1..n: m_{i→i+1} = collapse[ ψ_i + m_{i-1→i} + φ ]  (grid integral over θ_i).
5. BACKWARD sweep i=n..1: m_{i→i-1} = collapse[ ψ_i + m_{i+1→i} + φ ].
6. belief_i = ψ_i + m_{i-1→i} + m_{i+1→i}; read posterior-mean (f₊,f₋,f_g) (the solve_node combine).
7. regions -> CalibrationResult f_g; boundary nodes -> flux transport (priors).
```

**Open issues / potential problems (flagged for review):**

- **A — couple density, not fraction; state in fractions (RESOLVED in review).** The design review settled
  the make-or-break question. **State = fraction** (the pie): a propagated value can only *redistribute* a
  node's observed `U`, never inject phantom mass — the capping / no-injection safety. **Coupling = gDNA
  density** `ρ_g = f_g·U/L`, NOT the fraction. Reason: `f_g = ρ_g/(ρ_g+ρ_RNA)` is contaminated by the local
  RNA expression (the 10000×-variable, non-propagating quantity), so it is node-specific; the gDNA *density*
  is the library-wide contamination property (piecewise-uniform by enrichment class), so it is the
  transferable quantity. Decisive case: two adjacent on-target exons with the *same* `ρ_g` but very
  different expression have very different `f_g` (e.g. 0.05 vs 0.60) — propagating the *fraction* would call
  the low-expression AMBIG exon 5% gDNA (**exactly the gDNA→RNA leak we are fixing**, and capping does not
  save a high-count node); propagating the *density* recovers ~0.60. So capping bounds blast-radius (thin
  nodes), density-coupling removes the bias (high-count nodes) — complementary.
  The capture discontinuity (exon enriched, intron/intergenic baseline) is handled by **enrichment-matching
  the density**: exon-class and off-target-class form two interleaved chains (each still a tree → still
  exact in two sweeps), OR a single genomic chain coupling density in **off-target-equivalent units** (a
  per-node capture factor divided out — lean: single chain, reuses increment-2 adjacency). This is where the
  `density_model` capture logic ports into `φ_{ij}`.
- **B — boundary node = one or two?** The accumulator gives two side views (region `r`'s right, region
  `r+1`'s left) of the *same* crossing fragments. Model the boundary as **one** chain node (its observation
  = the shared crossing counts), coupled to each flank with that flank's eff-length? Or keep two nodes? One
  node is cleaner and avoids double-counting the crossing; confirm against the accumulator semantics.
- **C — spliced flux ownership.** A junction's spliced flux is a lower bound on the *exon-side region's*
  RNA, not the boundary's pie. Make sure the sided spliced factor attaches to the correct region node and
  is not double-applied to the boundary.
- **D — one-sided & simplex factors break Gaussianity.** The spliced lower bound (inequality) and the
  simplex constraint are non-Gaussian → the **grid message is the right first choice** (it represents them
  exactly); a Gaussian-message version would need EP-style moment matching. Defer Gaussian to phase 7.
- **E — parallel sweeps + parallel chains (confirmed).** The forward and backward sweeps are **mutually
  independent** (neither consumes the other's messages; the belief combines them only at the end), so the
  two directions run fully in parallel. Independent chains (loci/references) parallelize too. Only the scan
  *within* one direction is sequential (message `i` needs message `i−1`). Cost O(n · P) per sweep
  (P = lattice points), linear in nodes — cheap.
- **F — circularity of the count `var~mean` fit (iteration is sound).** The coupling variance is fitted from
  the deconvolution, which depends on the coupling. Start **one pass**: fit `σ²_couple` from the seeds
  (intergenic + single-strand) *before* the sweep, propagate once. The outer loop (fit `var~mean` from the
  deconvolution, re-propagate) is **theoretically sound** — iteratively-reweighted estimation / a poor-man's
  EM, typically settling in 1–2 outer passes — so it is a principled extension, not a hack, if one pass
  under-fits.
- **G — boundary transport subsumed?** Boundary nodes now carry gDNA mass directly; the post-hoc
  `priors._transport_boundary_flux` may become "read the boundary node's gDNA" rather than a separate mass
  shift. Keep transport as-is for phase 3 (no-regression); revisit in phase 6.
- **H — seed delta vs sharp factor.** Pin intergenic `f_g=1` as a near-delta local factor (numerically a
  very sharp grid spike), so the §3 "signal dies into a seed" behaviour emerges from the math rather than a
  special case in the loop.

**No-regression guards.** A single-strand chain with no AMBIG must reproduce today's strand result
(messages on `f_g` from a confident strand node dominate); mass is conserved per node by the simplex
(Σ=1); `gdna_none` conditions must stay false-positive-free (the gDNA prior + zero RNA evidence ⇒ low RNA).
