# Calibration — the bipartite region/boundary node graph (PLAN v6)

**Status:** authoritative architecture plan (2026-06-16). **Supersedes** `CALIBRATION_PLAN_v5.md` and the
`CALIBRATION_PHASE_BC_IMPLEMENTATION_PLAN.md` dry-run. v5 unified the per-node pie + reliability layer but left
two things ambiguous that the user has now made **firm decisions** on (session 2026-06-16):

1. **Boundaries are first-class nodes, not "future work."** The calibration graph is a **linear bipartite
   chain** of interleaving **Region** and **Boundary** nodes. Both own fragment mass *and* effective length.
   (v5 §4 wanted this; v5 §8.4 called it premature — §8.4 is now overruled.)
2. **`I₀` (`gdna_strand_info_scale`) is deleted.** No buried, un-interpretable magic constant ships. Deference
   of count→strand must be **emergent from honest precisions** (§6). A temporary accuracy regression during
   recalibration is **expected and accepted**.
3. **Spliced mass is owned by Boundary nodes only** (one-sided, stranded by the splice motif). Measured: on
   `quick_3to1_5mb/ss0.99`, **100.0% of spliced mass is on boundary objects, 0.0% in `region_contained`**
   (0 of 1221 regions carry any), boundary spliced **97.8% one-sided**. The user's model is empirically exact.

Where this and v5 disagree, this wins. The v5 *mechanisms* survive verbatim — the 2-simplex pie, the
grid sum-product, the per-component `var~mean` reliability, compose-by-precision. v6 changes the **graph** the
mechanisms run on (region-only chain → region↔boundary bipartite chain) and **retires `I₀`**.

---

## 0. TL;DR — the whole architecture in one picture

```
reference chain (genomic order):

  B0 ── R0 ── B1 ── R1 ── B2 ── R2 ── … ── R_{k-1} ── Bk
  └terminal┘      └seam┘     └seam┘            └terminal┘

  Region node R   : a genomic interval [start,end). Owns CONTAINED unspliced mass → pie {f₊,f₋,f_g}.
                    Effective length = E_g[max(0,L−ℓ)] (the contained support). Owns ~0 spliced.
  Boundary node B : a genomic point with a LEFT side (in R left) and a RIGHT side (in R right). Owns
                    CROSSING unspliced mass (two-sided) → pie {f₊,f₋,f_g}, AND CROSSING spliced mass
                    (one-sided, stranded by the GT-AG motif) = a FIXED mature-RNA floor.
                    Effective length = the per-side density length(s) E_g[min(ℓ,L_side)].
```

> **Every node — region OR boundary — carries a pie `(f₊,f₋,f_g)` over its UNSPLICED mass and is solved by
> the SAME 3-tier ψ.** A boundary additionally carries a one-sided spliced floor (mature RNA). **The
> communication between nodes is the PIE ITSELF (the node's state) — NOT transported mass, NOT a moved
> density.** A source node shares its current pie `(f₊,f₋,f_g)` (together with its owned mass + effective
> length — its geometry); the destination **integrates that pie into its own pie**, weighted by an honest,
> measured reliability. Nothing is moved between nodes. The reliability of each shared component is the
> `var~mean` machine (§7); deference among the signals is **emergent from honest precisions** — there are
> **no magic coupling/blend constants** (both `I₀` and `q_rna` are deleted, §6/§9). The deconvolution is the
> exact grid sum-product over the bipartite chain, iterated to a fixed point over the all-gDNA bootstrap.

**What this buys.** Boundary nodes are where the hybrid-capture signal is richest: an exon|intron seam is
**count-observable** (its crossing unspliced is gDNA+nascent, mature-free) and carries the **mature floor**
(its one-sided spliced). So boundary nodes are the anchors that resolve the non-observable exon region nodes
between them — the imputation backbone of the whole capture model. The region-only chain threw this away
(it solved regions and treated boundaries as mere edges / a one-shot `deconv_sides` side-calc).

---

## 1. Empirical grounding (measured this session, not assumed)

`scripts/debug/measure_spliced_ownership.py` on `quick_3to1_5mb/gdna_gdna300_ss_0.99_nrna_rnd_capture_on`
(4.18M deposited fragment-mass):

| quantity | region_contained | boundary objects |
|---|---|---|
| **spliced** (c2/c3) | **0.0 (0.000%)** | 286,621 (**100.000%**) |
| unspliced (c0/c1) | 2,503,050 (64.2%) | 1,392,803 (35.8%) |

- **0 of 1221 regions carry any contained spliced.** Reason (verified against the accumulator spec,
  `tests/native/_accumulator_reference.py:203-211`): a fragment is `region_contained` **iff all its slices
  fall in one region**; a multi-exon spliced read has its exon blocks in two *different* exon regions (the
  intron between them is its own region), so it is always a **crossing** whose spliced mass lands one-sided on
  the flanking exon|intron seam boundaries.
- **Boundary spliced is 97.8% one-sided** (596 boundaries carry spliced; the ~2.2% two-sidedness is the T5
  contiguous-spliced straddle of an internal same-signature seam + boundaries that are donor for one junction
  and acceptor for another — both legitimate, not a bug).
- **Corollary (a real latent defect this fixes):** the region-only sweep's spliced floor reads
  `substrate.contained.n_spliced_sense` ([`simplex_sweep.py:215`](../../src/rigel/calibration/simplex_sweep.py#L215))
  — which is **≈0 on realistic annotations**. The spliced floor is effectively dead today; the mature signal
  lives on boundaries the sweep never solves.
- **No simulation bug (Q-B checked):** the user flagged that *contained* spliced mass in **simulation** data
  would indicate a bug (sim has ground truth, all junctions annotated → every spliced read should be a
  crossing). The measurement found **exactly 0.0% contained spliced (0 of 1221 regions)** — the sim is clean,
  no bug. (In *real* data, contained spliced *can* legitimately occur from **un-annotated** junctions not in
  the reference GTF — not modeled by the region partition. For calibration, **contained spliced is ignored**
  either way: it is null in sim and out-of-model in real data. Rigel distinguishes **ANNOTATED**, **IMPLICIT**
  (intron-skip implied but not directly captured by the alignment), and **ARTIFACT** junctions — an alignment
  axis upstream of the calibration node model, which only consumes the boundary-owned spliced.)

---

## 2. The node model (what each node owns)

| | **Region node** | **Boundary node** |
|---|---|---|
| **geometry** | interval `[start,end)`, length `L` | point; two sides, in the left region (`L_left`) and right region (`L_right`) |
| **count** | R nodes (`region_contained`) | `B_obj_total` nodes (`boundary_*`); `k+1` per ref (2 terminal + `k−1` internal) |
| **unspliced mass** | contained `c0/c1` | crossing `c0/c1`, **two-sided** (a contiguous crossing straddles both flanks) |
| **spliced mass** | **≈0** (single-region spliced only) | crossing `c2/c3`, **one-sided** (mature on the exon flank), **stranded by motif** |
| **pie `{f₊,f₋,f_g}`** | over its contained unspliced mass | over its (pooled two-side) crossing unspliced mass |
| **eff length** | `E_g[max(0,L−ℓ)]` (contained support) | per-side density length `E_g[min(ℓ,L_side)]`; pooled support `½(E[min(ℓ,L_left)]+E[min(ℓ,L_right)])` |
| **count-observable?** | iff no exon bit (intergenic/intron) | iff no shared exon across its two sides (exon↔intron, exon↔intergenic) |
| **spliced floor** | from its ≈0 contained spliced | from its one-sided crossing spliced (the real mature anchor) |

**The data already exists.** Boundary nodes are built **directly** from the payload's per-boundary arrays
`boundary_mass_{left,right}` / `boundary_flux_{left,right}` ([`scan_payload.py:28-34`](../../src/rigel/scan_payload.py#L28)),
indexed by boundary; `region_arrays.boundary_region_indices` already maps each boundary → `(left_region,
right_region)` ([`region_arrays.py:156`](../../src/rigel/calibration/region_arrays.py#L156)). **No C++/accumulator
change is required for the linear bipartite graph** (see §3 — this is the load-bearing distinction).

**A boundary is ONE belief unit** (user Q5). Its two sides are mass-accounting (where the crossing landed), not
two beliefs: the boundary solves **one** pie `f_g_b` and applies it to both sides' masses. This both (a) matches
the physics (one physical junction) and (b) replaces the old `deconv_sides`, which solved the two sides as
*independent* nodes.

---

## 3. CRITICAL reconciliation — linear genomic-seam bipartite (THIS PLAN) vs transcript-adjacency SJ graph (DEFERRED)

This is the single most important thing to get right, because two different "bipartite" graphs are in play and
conflating them introduces a fatal, expensive misunderstanding.

| | **(A) Linear genomic-seam bipartite — v6, BUILD NOW** | **(B) Transcript-adjacency SJ graph — DEFERRED** |
|---|---|---|
| boundary = | a **genomic seam** between two genomically-adjacent regions | a **splice junction** connecting two transcript-adjacent exons across a skipped intron |
| topology | **linear** chain `R–B–R–B–…` (genomic order) | **non-linear** (cassette exons → cycles) |
| adjacency | from region partition order (already in the payload) | from intron coords `(start,end,strand)` — **thrown away before `deposit`** |
| C++ change | **NONE** (boundary objects already exist, spliced already one-sided) | **YES** — `deposit` signature must carry cut-intron coords; new SJ store; payload reshape; byte-for-byte reference rewrite |
| BP | exact two-sweep on the chain (a chain is a tree) | loopy-BP / junction-tree (cycles) |
| status | **this document** | `splice_junction_node_architecture.md` (verified-fatal-as-naively-designed; deferred) |

The user's Q1 ("a **linear** bipartite graph with interleaving region↔boundary↔region") is unambiguously **(A)**.
Q3's "a single boundary node can represent one and only one splice junction" is true **at the facet level**: a
genomic-seam boundary at an exon|intron transition owns **one facet** (donor *or* acceptor) of a junction's
one-sided spliced mass. The *full* junction (exon→exon across the intron) touches **two** seam boundaries
(donor seam + acceptor seam) and is the **transcript-adjacency edge** of graph (B) — which we are **not**
building now. **v6 builds (A); the transcript-adjacency edge that would let exon-1 borrow mature odds from
exon-2 across a long intron remains the deferred SJ-node feature.**

> **Guard against the trap:** do NOT key boundary nodes by intron coordinates, do NOT add transcript-adjacency
> edges, do NOT touch the accumulator. v6 is pure Python over the existing payload. If a future need for (B)
> arises, it is the separate, greenlit-gated initiative in `splice_junction_node_architecture.md`.

---

## 4. Spliced ownership, the one-sided floor, and spliced effective length

**Ownership (Q3/Q5).** In the calibration model the **boundary node owns its one-sided crossing spliced mass**
as a *fixed, stranded mature-RNA source*. Region nodes own ≈0 spliced (the rare single-region case, §1 —
handled as a region floor for completeness, but it is empirically null). Conservation is preserved exactly: the
accumulator already deposits each spliced fragment's mass fractionally across the boundary sides it crosses,
summing to 1 per fragment ([`_accumulator_reference.py:254-263`](../../tests/native/_accumulator_reference.py#L254)).
**No re-attribution, no double-count, no C++ change** — v6 simply *reads* the spliced from the boundary node
(where it already is) instead of from a (null) region channel.

**The one-sided floor in ψ.** At a boundary node, the same-strand spliced crossing mass on strand `s` is a
**hard lower bound** on that strand's RNA: it was directly observed as mature. The boundary's unspliced pie sets
the gDNA/RNA± split of its *unspliced* crossings *subject to that floor* (v5 §2, now correctly homed on the
boundary). The floor is **one-sided** — it exists only on the exon-facing side; the intron-facing side has
spliced=0, which is correct (no mature RNA crosses into an intron).

**Spliced effective length (Q3 — "must reflect the 1-sided nature").** The spliced floor is a *density* on the
exon flank: `ρ_spliced,s = spliced_s / E_rna[min(ℓ, L_exon_side)]` — divided by the **one** exon-side RNA
density length, NOT a two-sided sum and NOT the count length `rna_fl_mean = E_rna[ℓ]`. This is the RNA twin of
the gDNA per-side density length and the symmetric fix already shipped for the unspliced node-pair
([`variance_model.py:443-458`](../../src/rigel/calibration/variance_model.py#L443)). The floor is converted to
an unspliced-mass-fraction equivalent for the pie exactly as the RNA prior target is (§6, §8).

---

## 5. Naming clarity — "own" demystified (Q4)

The cryptic `own` in [`density_model.py:321`](../../src/rigel/calibration/density_model.py#L321)
(`own = region_count_observable & (region_eff_len > eps)`) means **"this node measures its own density
directly, from its own counts, because it is count-observable"** — as opposed to **"imputed"** (it has no
direct gDNA observation, so it borrows a density from its neighbours). A count-observable node (an intergenic/
intron region, or an exon↔intron boundary) is **gDNA by construction**, so its own crossing/contained count
*is* a clean gDNA-density measurement. A non-observable node (an exon region) is contaminated by mature RNA, so
it cannot measure its own gDNA — it is imputed from observable neighbours.

**v6 renames** for interpretability (no behavior change): `own` → **`directly_observed`** (or `is_direct`);
`imputed` stays; the DIRECT/IMPUTATION `var~mean` split (`fit_direct_varmean` / `fit_pair_imputation_varmean`)
already uses the right words — align the variable names to them. The principle: **a node is either
*directly-observed* (count-observable; measures its own component density) or *imputed* (borrows from chain
neighbours via the bipartite imputation).** This is the same DIRECT-vs-IMPUTATION distinction across regions
*and* boundaries now.

---

## 6. The unified per-node ψ — principled deference, NO `I₀` (Q2)

Every node (region or boundary) is solved by the **same 3-tier ψ** on the 2-simplex lattice:

```
ψ(f₊,f₋,f_g) =  L_strand(tilt | u₊, U, κ, od_g, od_r)         # LIKELIHOOD — capture-invariant DIRECTION
             +  [ one-sided spliced floor on f_s ]            # HARD — boundary's mature anchor (gDNA-free)
             −  ½·τ_count·(f_g − μ_g)²                        # PRIOR — gDNA magnitude (imputed pie, var~mean)
             −  ½·τ_rna₊·(f₊ − μ₊)² − ½·τ_rna₋·(f₋ − μ₋)²      # PRIOR — RNA magnitude (imputed pie, var~mean)
             −  ½·τ_global·(f_g − μ_global)²                  # PRIOR — global gDNA hyperprior (MAD spread)
```

### 6.1 Why `I₀` was there, and why it can go

Today the count precision is artificially suppressed by `(1−w_strand)`, `w_strand = I/(I+I₀)`
([`simplex_sweep.py:234`](../../src/rigel/calibration/simplex_sweep.py#L234)) — an *extra* hand-tuned deference
on top of the honest precisions, to stop the count's `~N` precision over-riding the solve at AMBIG. It is
exactly the brittle, un-interpretable constant the user will not ship.

**The principled replacement is plain Bayes — no weight.** The three signals carry their *honest* precisions and
compose additively in log-space; the posterior is dominated by whichever is sharpest:

- **`L_strand` is a LIKELIHOOD**, not a prior. Its precision for the gDNA-vs-RNA split is the **intrinsic
  curvature** of the Beta-Binomial loglik w.r.t. `f_g` via the tilt — `I_strand ∝ N·(2κ−1)²`. This is **already
  computed honestly** by `_mixture_strand_loglik` ([`simplex.py:61`](../../src/rigel/calibration/simplex.py#L61));
  no knob sets it. At a confident single-strand node it is sharp and **dominates** any diffuse prior. At κ=½
  (unstranded) or AMBIG (no defined sense) it is **flat** (`I_strand→0`) and yields to the priors. *This is the
  "strand is the reliable, unbiased priority" principle, realized by the math, not a constant.*
- **`τ_count`** (gDNA magnitude) is honest **iff the `var~mean` reports high variance where the count is
  unreliable** — which it does by construction at non-observable/AMBIG nodes (the boundary→region imputation
  error is large there → low `τ_count`). *This is the "count is important but unreliable, defers to strand"
  principle: it defers automatically because its own reliability layer humbles it where it can't see.*
- **`τ_global`** (the global gDNA hyperprior) is the MAD between-node spread — the population width an
  unanchored node faces. *This is the "global hyperprior must be present to stabilize" principle.* It governs
  only nodes the likelihood and the local priors both leave undetermined.

**So deference is emergent:** sharp likelihood ⇒ strand wins; flat likelihood + reliable count ⇒ count wins;
flat likelihood + unreliable count ⇒ global wins. No `I₀`.

### 6.2 Why this is *safe enough to commit to* (and the recalibration it demands)

The risk the dry-run flagged (R2): removing `(1−w_strand)` lets `τ_count` over-ride propagation at AMBIG, the
historical complex-locus win. v6 has **two new defenses the old region-only chain lacked**, which is why the
user's "we'll recalibrate" is well-founded:

1. **The RNA prior (`τ_rna±`, §8) gives AMBIG nodes a magnitude signal the count previously monopolized.** At an
   AMBIG node the strand can't split ±, but its flanking boundaries' one-sided spliced + stranded unspliced RNA
   *can* — so the RNA prior, not the raw count, carries the AMBIG magnitude. The count no longer has to be
   suppressed by hand because it is no longer the only magnitude voice.
2. **Boundary nodes are now solved**, so an AMBIG region sits between *solved* boundary anchors whose
   imputation `τ` is honest — the var~mean humbling is computed on real bracketed neighbour estimates, not the
   frozen one-shot `deconv_sides` values.

**The recalibration (expected, accepted regression):** the var~mean must be *trusted* to humble the count
correctly. The single must-do is to **verify, on the complex AMBIG battery + net-flow, that with `I₀` removed
the var~mean `τ_count` at AMBIG is low enough that strand+RNA-prior+propagation govern** — and if a residual
over-call remains, the fix is to *improve the reliability layer* (the honest lever), never to re-introduce a
suppression constant. `gdna_strand_info_scale` is **deleted from `config.py`** at the end of Phase C.

> Note `cleaned_gdna_count`'s upstream `I₀` ([`strand_deconv.py:517`](../../src/rigel/calibration/strand_deconv.py#L517))
> is a *different* use (it cleans the boundary count *before* the count module). It is retired together with
> `deconv_sides` when the boundary node solve replaces that whole pre-pass (Phase C, §12) — the boundary node's
> own ψ does the strand/count combine by precision, so no separate cleaning step survives.

---

## 7. The reliability layer (the `var~mean` fits) on the bipartite graph

Unchanged in spirit from v5 §3; extended to the bipartite graph:

- **gDNA** (`fit_pair_imputation_varmean`, LIVE) and **RNA per-strand** (`fit_pair_imputation_rna_varmean`,
  INERT → activated). One node-pair point per `(observable source → imputed dest)` adjacency,
  `raw_var = (d_dest − d_source)²` in density space, dof=1, Bernoulli-clamped. DIRECT (own-count Poisson) for
  directly-observed nodes.
- **NOW BIDIRECTIONAL (the v5 §3 deferral lifts).** Because boundaries are *solved nodes* that query `τ` at
  `μ = boundary_density`, the `region→boundary` direction is now a real queried axis — so the builder must emit
  **both** directions (`boundary→region` and `region→boundary`). v5 §3 deferred this "until boundaries become
  solved nodes in Phase C" — that is now. Add a `bidirectional=True` mode (the points are symmetric pairs; the
  GCV/IRLS-corruption concern v5 raised applied only when one direction's axis was never queried — both are now).
- **Bracketed** (fit on the previous pass's estimate) — the two-prior convergence guard (§11).

---

## 8. The message is the PIE (the node's state) — no transport (the mechanism, user-clarified 2026-06-16)

> **A source node shares its current pie `(f₊,f₋,f_g)` with a destination node; the destination integrates
> that pie into its own pie. No mass is transferred; no density is moved.** The communicated thing is the
> node's *state*. (The user's explicit correction: we are *not* "communicating density and transferring mass"
> — we are "communicating the current state of the node.")

**Why the pie alone is not enough — and what "integrate" means.** A node also *owns its mass and effective
length* (§2). The shared state is therefore the pie **plus that geometry**, and the destination integrates by
recovering the comparable, physically-smooth quantity from it:

- **gDNA — compare in DENSITY, because the gDNA fraction is NOT directly comparable across node types.** An
  intron's `f_g≈1` and an adjacent exon's low `f_g` differ for a *biological* reason (the exon has RNA), so
  naively pulling the exon's `f_g` toward the intron's is wrong. The smooth, comparable quantity is the gDNA
  **density** `ρ_g = f_g·M/E_g`, recovered from the shared pie + the neighbour's owned mass/eff-len. gDNA is
  genomically smooth (modulo spatially-smooth capture enrichment), so a neighbour's `ρ_g` is an honest
  predictor of the destination's `ρ_g` → the destination's `f_g` prior `μ_g = ρ̂_g·E_g/M_u`. **This is not
  transport** — the neighbour's mass stays the neighbour's; only its *state* (pie + geometry) is read to form a
  prior. Reliability `τ_count` = the gDNA `var~mean` (the measured density-residual; §7).
- **RNA — continue along same-strand exons only.** RNA is spiky/isoform-driven (the user's foundational
  correction — RNA does **not** chain globally). The same-strand RNA continues across an exon↔exon-or-boundary
  step; the destination's `μ₊/μ₋` prior is the neighbour's same-strand RNA density (incl. the boundary's
  one-sided spliced) carried over, `μ_s = ρ̂_s·E_rna/M_u`. Reliability `τ_rna_s` = the RNA `var~mean`. Where no
  same-strand RNA neighbour exists, `τ_rna_s = 0` (prior off). No-anchor fallback = 0, never a global RNA density.
- **spliced floor:** the boundary's own one-sided spliced (§4), a hard lower bound on its own pie's RNA — local,
  not communicated.

**Coherence is structural:** one node, one pie on the 2-simplex lattice, `f₊+f₋+f_g=1` by construction. The
destination's *posterior* pie is its own ψ (strand likelihood + spliced floor + global foundation) combined
with the neighbour-derived component priors — always a normalized point on the simplex. No coherence guard.

> **Density is the *unit of comparison/reliability*, not a *transferred quantity*.** Fitting the `var~mean` on
> density residuals (the proven factor-1-under-uniform machine) measures *how predictive* a neighbour's pie is;
> it moves nothing. This is the precise sense in which we "communicate state, not mass."

---

## 9. Propagation — ONE reliability-weighted mechanism; `q_rna` DELETED (Q-D, firm)

`q_rna = 0.25` is **deleted**, for the same reason as `I₀`: an un-interpretable, brittle-in-production magic
constant (Q-D, firm). There is **one** cross-node mechanism — the pie-communication of §8, weighted by the
honest `var~mean` reliability — and it must be **parameter-free**. Two principled realizations; **prototype and
measure, the scalar 0.25 is gone either way:**

- **(i) Iterated one-hop priors (simplest; removes the `(P,P)` edge machinery entirely).** Each pass, neighbour
  pies set the per-node component priors `(μ, τ)` in ψ (§8); solve every node from its ψ; refit; iterate.
  **Transitive continuity emerges from iteration** — a node's pie informs its neighbour next pass, that
  neighbour's the next-next pass. No within-pass log-odds edge, no `q_rna`, no `(P,P)` transition matrices. The
  grid sum-product collapses to independent per-node lattice solves with neighbour-derived priors. This is the
  most direct reading of "share the pie, integrate it."
- **(ii) Derived edge coupling (keeps within-pass propagation, parameter-free).** If a within-pass coupling
  proves to add value (Phase-1 measured the *old* region-only `q_rna` helped +1–5%), keep the grid sum-product
  edge but set its coupling variance **from the RNA `var~mean`** (the per-edge measured log-odds-difference
  variance), not a scalar. `q_rna` becomes the same reliability machine that gives `τ_rna` — one source, no
  constant.

**Recommendation:** prototype **(i)** first (it deletes the most machinery and is the cleanest realization of
the user's framing), then A/B **(ii)** — does a derived within-pass edge beat pure iterated one-hop on the
bipartite graph? Whatever wins, `q_rna=0.25` does not survive. *(The Phase-1 +1–5% was measured on the
region-only chain without boundary nodes or an RNA prior; both new mechanisms may already supply that
continuity, making (i) sufficient — must be re-measured.)*

**Tree-exactness** (if (ii) is kept): the bipartite chain is still a chain (a tree) → forward+backward exact,
the v5 guarantee holds. (Cycles arise only in the deferred transcript-adjacency graph (B), §3.)

---

## 10. The result schema + `assemble_priors` — first-class per-node fields (Q-E, REQUIRED)

The output contract **must change** (user, firm): stop emitting per-region `mass_{gdna,rna}_{contained,left,right}`
and instead make **every node — region AND boundary — a first-class result element** carrying its pie and
geometry. The eff-len IPR code already did this promotion for the gDNA node set
([`priors._gdna_region_node_arrays`](../../src/rigel/calibration/priors.py#L103) builds region + pooled-seam
nodes with mass + effective support) — **it is the template**.

**New `CalibrationResult` (per-node, node set = regions ∪ boundaries):** each node emits
- `f_pos, f_neg, f_g` — its pie (the fractions, the canonical state);
- `mass_gdna = f_g·M_u`, `mass_rna_pos = f₊·M_u`, `mass_rna_neg = f₋·M_u` — the pie applied to its unspliced
  mass (the strand split is now carried, not collapsed to one `mass_rna`);
- `mass_spliced` — boundary nodes only (the one-sided mature floor), region nodes ≈0;
- `eff_len` — its effective support (region: `E_g[max(0,L−ℓ)]`; boundary: the per-side density length(s)).

Mass and effective length are **aggregable**, so densities are computable at any granularity (node, locus,
library). The layout (a single unified node array with a `is_boundary` flag, or parallel region/boundary
arrays) is a Phase-C detail; the boundary↔region adjacency (`region_boundary_indices` /
`boundary_region_indices`) keys the aggregation.

**`assemble_priors` rewrite:** aggregate per-locus over the **node set** (regions + boundaries) by genomic
overlap — `gdna_prior_count = Σ mass_gdna`, `rna_prior_count = Σ (mass_rna_pos + mass_rna_neg)` (spliced still
withheld — guaranteed-RNA in the EM), `gdna_eff_len` = the IPR over the node masses/supports (the existing
Laplace-smoothed, factor-1-under-uniform form, now reading the first-class boundary nodes directly instead of
re-deriving the pooled seam from `mass_gdna_{left,right}`). This **subsumes** the current `_gdna_region_node_arrays`
seam-pooling: boundaries *are* the seam nodes now. The bedrock invariant (factor-1-under-uniform-gDNA) and the
contained-evidence multimap-blindness shrinkage are preserved verbatim — they already operate on per-node
masses/supports.

*Consequence:* `priors.py`, `capture_eff_length.py` (shares `_gdna_region_node_arrays`), `result.py`, and the
`pipeline.py` QC totals all touch the new schema — this is the largest blast radius of v6 and lands in Phase C.
The intrinsic invariants (`result.py::__post_init__`: non-negative, finite, shapes) carry over per-node.

---

## 11. The two-prior convergence must-pass (carried from v5 §7)

Two iterating priors (gDNA + RNA) on one simplex is a coupling the region-only 2a path lacked. Bracketing +
frozen κ/overdispersions should make it contract, but **assert monotone convergence** on a capture-on+nascent
toy (the `calibration_feature_dev_toy_scenarios` recipe — small genome, multi-exon strand training, an
overlapping opposite-strand AMBIG node with splice carry-over, gDNA swept 0→high + nascent) before trusting.
Damp by under-relaxing the RNA estimate fed to the next fit (`ρ̂_next = (1−η)ρ̂_prev + η·ρ̂_new`) **only if it
oscillates** — documented as a numerical-stability constant, not a model tunable.

---

## 12. Implementation sequence (prototype first, then optimize — Q6)

Each phase is independently committable/revertable; we commit to `origin/main` directly (no branches). A
temporary accuracy regression during the `I₀`/`q_rna` recalibration is **expected and accepted** (user) — to be
recovered by improving the reliability layer, never by re-adding a suppression constant.

- **Phase A — boundary node substrate + RNA fit activation + payload cache (plumbing, low risk).**
  - New `BoundaryArrays` / boundary-node view: build boundary nodes directly from `boundary_mass/flux_{left,right}`
    indexed by boundary; pool the two sides for the belief, keep per-side mass for accounting. Reuse
    `boundary_region_indices` for adjacency, `_side_strand_orientation` for sense.
  - Activate `fit_pair_imputation_rna_varmean` in the pass loop (compute `τ_rna±`, log; not yet consumed).
  - **Payload cache utility** (Q6): persist the scanned payload so calibration re-runs without re-scanning the
    500 MB BAM — the foundation for thousands of fast iterations.
  - *Validate:* boundary nodes reconstruct the old `substrate.left/right` exactly (a re-indexing identity);
    suite green; RNA fit real/monotone in-loop. **No solver change yet.**

- **Phase B — the RNA prior + the one-sided spliced floor in ψ (region chain first).**
  - `node_rna_density` (RNA analog of `node_gdna_density`): per-node per-strand imputed RNA density for the
    apply set. `μ±/τ_rna±` (§8). One-sided eff-len-consistent spliced floor (§4, Q-A confirmed).
  - `_local_loglik` += two RNA pulls (gated on `μ± is not None`, toy back-compat).
  - *Validate:* zero-RNA bit-identical; coherence/conservation; complex battery ≤ 4502; net-flow + mature FP;
    two-prior convergence (§11). Goldens regenerated post-gate.

- **Phase C — the unified bipartite solver + the schema change + delete `I₀` and `q_rna`.**
  - **Interleave boundary nodes** into the per-reference chain; solve every node (region + boundary) by the same
    ψ (their spliced floor + direct-Poisson count prior where observable + node-class prior).
  - **Delete `q_rna` (Q-D):** realize the cross-node communication as **(i) iterated one-hop priors** (no `(P,P)`
    edge, no `q_rna`); A/B against **(ii) a `var~mean`-derived edge coupling** only if (i) under-performs (§9).
  - **Delete `I₀` (Q2):** full `τ_count`, no `(1−w_strand)`; the **`I₀`-removal A/B** (§6.2) confirms the
    var~mean humbles the count at AMBIG. `compose_by_precision` is the transitional flag, then the old path is
    deleted.
  - **Schema change (§10, Q-E, REQUIRED):** promote boundary nodes to first-class result fields; rewrite
    `CalibrationResult` + `assemble_priors` (+ `capture_eff_length`, `pipeline` QC) to aggregate per-node pies
    (`mass_gdna`, `mass_rna_pos/neg`, `mass_spliced`, `eff_len`); stop emitting `{contained,left,right}`.
  - Un-freeze the var~mean predictor (boundaries iterate). Enable bidirectional node-pair (§7).
  - **Delete `gdna_strand_info_scale` + `q_rna` from config.py; delete `deconv_sides`/`_deconv_per_node`/`cleaned_gdna_count`.**
  - *Validate:* battery + net-flow + FP + convergence + eff-len bedrock; accept a traced/recovered recalibration
    regression; the schema-conservation (`Σ node mass = total`) byte-level check.

- **Phase D — performance (after a correct prototype).** §13.

- **Phase E — (separately) the deferred transcript-adjacency SJ graph (B), §3.** Gated on
  `splice_junction_node_architecture.md`'s greenlight conditions (needs the C++ deposit-signature change). Not
  required for v6.

---

## 13. Performance plan (Q6 — measure; thousands of runs)

Performance is a first-class requirement (the tool runs thousands of times during development). **Prototype the
bipartite solver in Python/numpy first, then optimize.** Concrete plan:

- **Instrument now:** a `--profile-calibration` path that reports per-phase wall time (substrate, strand fits,
  the two `var~mean` refits/pass, the grid sum-product/pass, #passes). v5 §8.6 already flags the refit as the
  bottleneck and the `_sweep` matrix materialization as a fixed OOM wall (edge-matrix caching, decoupled `O(P)`
  short-circuit — already in `_sweep_chain`).
- **Expected hotspots after the bipartite change:** (a) **2× node count** (regions + boundaries) in the chain →
  more lattice evaluations; (b) **two var~mean fits/pass** (gDNA + RNA), now bidirectional → ~2× the fit points.
- **Optimization levers (deferred to Phase D, in ROI order):** fit caching across passes (the curve changes
  slowly); vectorize the per-node ψ over the whole chain (already mostly array-wise); a C++ kernel for
  `_local_loglik` / the grid sum-product only if profiling shows it dominates at genome scale. The var~mean SCAM
  fit is ~10³–10⁴ points, sub-second — likely *not* the genome-scale bottleneck; the per-node ψ over millions of
  nodes is. **Measure before porting.**
- **Benchmark harness:** keep `quick_3to1_5mb` as the fast iteration suite (scan once, cache the payload to
  disk, re-run calibration from the cached payload to avoid re-scanning the 500 MB BAM each iteration). A
  payload-cache utility is worth building in Phase A (it accelerates every subsequent measurement).

---

## 14. Magic-number ledger (v6 state)

| constant | verdict |
|---|---|
| **`I₀` `gdna_strand_info_scale`** | **DELETED** (Q2, firm) — all uses retired: the sweep `(1−w_strand)` (compose-by-precision, §6), `deconv_sides` (replaced by the boundary node ψ), `cleaned_gdna_count` (retired with `deconv_sides`). |
| **`q_rna` `0.25`** | **DELETED** (Q-D, firm) — the cross-node communication is the `var~mean`-weighted pie message; realized as iterated one-hop priors (§9 (i)) or a `var~mean`-derived edge (§9 (ii)). No scalar coupling survives. |
| **`σ²_global`** | **KEEP MAD** — the principled between-node population spread (the global gDNA hyperprior the user requires, §6; v4-measured, `DIRECT.predict` reverted). |
| **Jensen `Δ_k`, Bernoulli clamp, MAD `1.4826`, SCAM `k=18`/GCV, lattice `K=60`, `_PRIOR_EPS`** | canonical / numerical-resolution — kept, documented. |
| **under-relaxation `η`** | introduced **only if** the two-prior solve oscillates (§11); numerical, not a model tunable. |
| **the coherence guard** | **RETIRED** — structural (one node, one pie). |

**After v6 the calibration carries no un-interpretable model tunable.** Every remaining constant is canonical
(Jeffreys ½, MAD 1.4826), a numerical resolution (lattice K, GCV grid, EPS), or a measured reliability (the
`var~mean`). The two brittle hyperparameters the user flagged — `I₀` and `q_rna` — are both gone.

---

## 15. Resolved this session + the remaining live questions

**Resolved (user, 2026-06-16) — now firm in the plan:**
- **Q1 → C2 firm:** linear bipartite region↔boundary chain; boundaries first-class (§0–§2).
- **Q2 → `I₀` deleted, deference emergent** (§6); recalibration regression accepted.
- **Q3/Q5 → spliced boundary-owned, one-sided, motif-stranded** (§4); measured 100%/0% (§1).
- **Q-A → one-sided spliced eff-len** from the side the spliced mass is allocated to (§4). ✔
- **Q-B → no sim bug** (0% contained spliced measured); contained spliced ignored for calibration; ANNOTATED/
  IMPLICIT/ARTIFACT is an upstream alignment axis (§1). ✔
- **Q-C → direct count = Poisson now**, overdispersed-NB later (deferred) (§6). ✔
- **Q-D → `q_rna` deleted** (§9); communication is the `var~mean`-weighted pie message. ✔
- **Q-E → first-class per-node result fields**, `assemble_priors` rewrite (§10), modeled on the IPR code. ✔
- **Mechanism → communicate the PIE (state), no transport** (§8). ✔

**Remaining live design decisions (settle during execution, flagged):**
- **§9 (i) vs (ii):** prototype iterated one-hop priors (no edge, no `q_rna`) first; A/B a `var~mean`-derived
  edge only if it under-performs. *(Recommend (i) first.)*
- **§10 schema layout:** unified node array with `is_boundary` flag vs parallel region/boundary arrays.
  *(Recommend whichever `assemble_priors`' per-locus overlap aggregation reads most cleanly — likely parallel
  arrays keyed by the existing region/boundary offsets.)*
- **§6.2 / §11:** the size of the accepted transitional regression from removing `I₀`+`q_rna` together — measure
  on the complex battery + net-flow and confirm it is recovered by the reliability layer, not papered over.

---

## 16. Validation gates (every build phase)

- **Suite** green; goldens regenerated only after the accuracy gate.
- **Coherence** (`f₊+f₋+f_g=1` structural; one-sided spliced floor respected); **conservation**
  `gdna+rna = total` per node and `Σ all spliced = n_spliced_fragments` (the byte-for-byte accumulator gate is
  untouched — no C++ change).
- **Complex battery** TOTAL ≤ 4502 (the 2a baseline) — **except** a *traced, recovered* recalibration
  regression in Phase C (§6.2, §11) is acceptable.
- **Net-flow** non-regressing on capture-on ss0.99/ss0.5 + zero-DNA (identical re-simmed BAMs) — **mature-RNA
  FP rate must not rise** (the RNA prior must not manufacture RNA).
- **Monotone convergence** of the two-prior simplex (§11).
- **eff-len bedrock**: factor-1-under-uniform-gDNA preserved (`priors.py`).

---

## 17. Doc status

- **Authoritative:** this plan (v6).
- **Superseded:** `CALIBRATION_PLAN_v5.md` (boundaries-as-future-work / `I₀`-survives framings), the
  `CALIBRATION_PHASE_BC_IMPLEMENTATION_PLAN.md` dry-run (folded in here).
- **Deferred (NOT this plan):** `splice_junction_node_architecture.md` — the transcript-adjacency SJ graph
  (graph (B), §3), requiring the C++ deposit-signature change. Unchanged status.
- **Keep:** `per_node_deconv_hierarchy_design.md` §3 (per-node math); `effective_length_redesign_plan.md` §8
  (the eff-len IPR + pooled-seam supports — the boundary-node eff-len uses exactly these).
- Update `docs/README.md` to point at v6.
