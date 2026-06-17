# Calibration — the unified node solver with full-pie message passing (PLAN v5)

**Status:** authoritative architecture plan (2026-06-16). **Supersedes** `CALIBRATION_PLAN_v4.md` and
`phase2b_rna_imputation_plan.md`. v4 shipped the gDNA-side var~mean redesign (Phase 2a) and scoped the RNA
imputation + `q_rna` as separate bolt-ons; v5 unifies them: **one imputation that passes a node's whole pie,
with a per-component reliability layer**, which is the same thing as the unified node solver — so they are
built as **one step** (user decision, 2026-06-16). Where this and v4 disagree, this wins. Shipped code
(`main@8d9df6bd`, Phase 0 + 2a) is the foundation and is unchanged by this document until the build phases.

---

## 0. TL;DR — the whole architecture collapses to one idea

> **Every genomic node — region OR boundary — carries a pie `(f₊, f₋, f_g)` over its UNSPLICED mass.**
> **Imputation passes the full pie between adjacent nodes (bidirectional, node-agnostic), projected onto the
> destination's observed counts.** A **per-component `var~mean` reliability layer** (gDNA + per-strand RNA, both
> fit on all pairs from the current node state, bracketed) says how much to trust each axis of an imputed pie.
> **Spliced fragments are a fixed, direct RNA floor** at the nodes that own them; the pie determines only the
> *unspliced* split, subject to that floor. The exact grid sum-product BP then solves all nodes jointly.

This single mechanism subsumes four things v4 treated separately: the gDNA imputation (2a — the `f_g` slice of
the pie), the RNA imputation (the `f±` slices), the boundary-side solve (`deconv_sides` — boundaries are just
nodes), and the `I₀` strand/count blend (emergent from the per-component precisions). It also **retires the
v4/2b §3 coherence guard**: one deconvolution → one pie → the components partition and sum to 1 *by
construction*, not by assertion.

---

## 1. Foundation (shipped, unchanged until the build)

- **The exact grid sum-product engine** on the 2-simplex K-lattice, per-reference chains, edge-matrix caching,
  decoupled-edge `O(P)` short-circuit, posterior **median** readout + variance (`simplex_sweep.py`).
- **The `var~mean` machine** (`MonotoneVarMean`: monotone SCAM P-spline, GCV-λ, Jensen df-offset, power-law
  fallback) and the **gDNA imputation** (2a): boundary unspliced gDNA → region gDNA, all-pairs, current-density
  axis, Bernoulli-clamped → `τ_count`. **This is the `f_g` slice of the v5 pie.**
- **`σ²_global` = the robust MAD between-node density spread** (the foundation-prior width; the principled
  population spread, not a single node's sampling variance — v4 measured `DIRECT.predict(ρ_global)` ~5.5× too
  tight, +11.3% battery, reverted).
- **Frozen, capture-invariant:** κ (`rna_sense_frac`) and the two strand Beta-Binomial overdispersions. Only the
  *magnitude/composition* quantities iterate.

---

## 2. The unified pie imputation (the core mechanism)

**Every node is a pie over its UNSPLICED mass.** The simplex `(f₊, f₋, f_g)` is the gDNA/RNA± split of a node's
*unspliced* fragments. **Spliced fragments are separate, fixed, and pure mature RNA** (gDNA cannot splice; the
splice motif gives the strand) — they sit *on top of* the unspliced pie as a known RNA quantity, not inside it.

**Imputation = project a source node's pie onto a destination node's observed counts**, in **density space**:

1. The source node's current pie `(f₊, f₋, f_g)` + its RNA evidence per strand `s` = its **same-strand
   (spliced + unspliced) RNA density** `ρ̂_rna,s = (spliced_s + unspliced_RNA_s) / eff_rna` and its gDNA density
   `ρ̂_g = unspliced_gDNA / eff_g`.
2. Project onto the destination: the imputed destination composition is the source's *densities* carried over
   (gDNA is smooth; RNA continues along the transcript), then **multiplied onto the destination's own observed
   counts** — so the fraction is what propagates and the **smear is capped by the destination's evidence** (the
   reason the latent is a pie, not a count).
3. **Bidirectional and node-agnostic:** every adjacent `(node, node)` pair imputes both ways — region→boundary
   *and* boundary→region — exactly as the BP sweep already runs forward + backward.

**The spliced floor (load-bearing — vindicated by the user's `T1+/T2+` example).** At a node that owns spliced
mass, the spliced count is a **fixed lower bound** on that strand's RNA (it was directly observed as mature). The
projected pie determines the **unspliced** split *subject to that floor*: if the projection implies less RNA than
the spliced floor, the spliced (direct observation) dominates. Worked example: region `(1000,2000)` pie
`{f₊=0.9, f₋=0, f_g=0.1}` projected onto boundary `2000` (200 counts, 160 spliced fixed) ⇒ 180 RNA₊ total ⇒
spliced 160 (fixed) + unspliced RNA₊ 20 + gDNA 20. The 180 ≥ 160 consistency is the floor; the imputation only
ever sets the *unspliced* 40-count split.

**FL/density conversion (load-bearing — flagged by the user).** Spliced and unspliced cross with **different RNA
effective lengths**; the projection must be in density space with the correct FL marginals so the spliced and
unspliced add coherently. Reuse the `splice_junction.py` FL-consistency machinery (`region_eff_len_rna`,
`rna_fl_mean`, the `ρ_RNA` density form) — do not impute on raw counts.

**Coherence is now structural** (retires the v4/2b §3 guard): because one imputation carries the *whole* pie from
*one* deconvolution, `f₊ + f₋ + f_g = 1` holds by construction and the gDNA and RNA priors can never encode
contradictory boundary splits. No assertion needed.

---

## 3. The per-component reliability layer (the `var~mean` fits)

> **CORRECTION (2026-06-16, user): the imputation reliability is a NODE-PAIR prediction error, not a triplet
> disagreement.** A pair is **exactly two adjacent nodes — one region, one boundary**. For each adjacency
> `(observable boundary → imputed region)`, `raw_var = (d_region − d_boundary)²` in **density space**
> (eff-len-normalized, log-log; the boundary density *is* the single-predictor of the region), **one point per
> adjacency at `mean = d_region`** — the *queried* axis (the sweep only reads `τ_count` at region density). It is
> **NOT bidirectional for the fit** (the boundary-dest direction sits at an axis never queried until boundaries
> become solved nodes in Phase C; emitting it would also corrupt GCV/IRLS via duplicated residuals). A region
> with one observable flank still contributes a point (the densification the both-sides triplet missed). **One
> builder, symmetric** for gDNA and RNA — they differ only in the data. **Two regression vectors to A/B-isolate
> vs the shipped triplet** (re-opening 2a): (i) the full residual `(d_region−d_boundary)²` is ~4× the triplet's
> variance-of-the-mean `¼(d_L−d_R)²` ⇒ `τ_count` ~4× smaller ⇒ the count defers more to strand+global (correct
> for a single-predictor imputation error, but a real level change); (ii) the node-pair compares a
> region-*contained* density against a boundary-*crossing* density, which the triplet's crossing-vs-crossing
> form self-cancelled — sound under the factor-1-under-uniform eff-len bedrock but must be **measured** (the
> net-flow gate; hybrid fallback = keep the crossing-vs-crossing numerator emitted at `μ=d_region` if it
> regresses). The gDNA boundary density uses the **deconv'd** `gdna_mass/fl_mean` (RNA-removed), not raw
> crossings. This **supersedes the shipped 2a `fit_imputation_varmean_current`**, which wrongly used the
> boundary↔region↔boundary **triplet** disagreement `¼(d_L−d_R)²` (it required *both* flanking sides — the very
> "both-sides" model we abandoned, and the cause of the RNA sparsity). `DIRECT` (own-count Poisson) stays for
> observable nodes; the node-pair model is the **imputation** reliability for imputed nodes. Correcting 2a
> re-opens its validation gate (battery ≤ 4502, net-flow, goldens) and lets the triplet machinery
> (`varmean_points`'s both-sides logic, the observability-split) be **deleted** (the cleanup).

The pie message says *what* to impute; the reliability layer says *how much to trust each component*. Both
components use the **same node-pair machine** (`MonotoneVarMean`, current-density axis, Jensen `dof=1` on the
node-pair residual, Bernoulli clamp, all adjacent pairs both directions, **bracketed** — fit on the previous
pass's estimate, never the pass being solved), differing only in the data:

| component | predictor (source, density) | target (destination) | data scope |
|---|---|---|---|
| **gDNA** (2a, shipped) | boundary unspliced gDNA | region/node gDNA | all pairs with an observable gDNA crossing |
| **RNA, per strand `s`** (new) | source `(spliced_s + unspliced_RNA_s)` | dest unspliced RNA on `s` | **all pairs**, both directions, **excluding pairs with 0 RNA on `s`** |

- **Per-strand = two independent observations per pair**, both pooled into one RNA reliability curve (the
  reliability is a property of the imputation geometry, not the strand identity). Corrected per-strand form
  (same-strand evidence on both sides): `log(dest RNA_s) ~ log(impute(src spliced_s + src unspliced_RNA_s))`.
- **Iterates with the solve:** RNA mass is a current estimate, refined each pass; the curve improves as the
  deconvolution sharpens — identical bootstrap to gDNA. (At the all-gDNA init, RNA ≈ 0 ⇒ the RNA reliability is
  inert ⇒ the prior degrades gracefully to the spliced floor; it activates from pass 1+.)
- **Data-richness confirmed** (`phase2b_activation_rate.py`): the RNA imputation activates on ~52% of
  single-strand exon nodes on flagship-like genomes (52 pairs, real spline fit) and 15.8% on the AMBIG-heavy
  battery (20 pairs, just clears the 18-point spline-vs-power-law threshold). 3–13× smaller than the gDNA set
  but never starved; the `MonotoneVarMean` power-law fallback handles the thin-data tail.

These reliabilities weight the corresponding components of the unified pie message as **per-component precision-
weighted priors** in the node solve (§4). The gDNA prior pulls `f_g`; the per-strand RNA priors pull `f±`; the
spliced floor is the hard constraint underneath.

---

## 4. The unified node solver (the one-step build)

Regions and boundaries become **one per-reference chain of nodes** (`region–boundary–region–…`). Each node —
region or boundary — is solved by the **same 3-tier `ψ`** on the 2-simplex lattice:

```
ψ(f₊,f₋,f_g) =  L_strand(tilt)                              # capture-invariant DIRECTION (frozen κ, overdispersions)
             +  [ spliced floor: f_s ≥ observed spliced density on s ]   # hard, direct, gDNA-free (where owned)
             −  ½·τ_count ·(f_g − μ_g )²                     # gDNA component of the imputed pie (2a reliability)
             −  ½·τ_rna₊ ·(f₊  − μ₊ )²  −  ½·τ_rna₋·(f₋ − μ₋)²  # RNA components of the imputed pie (new reliability)
             −  ½·τ_node ·(f_g − μ_node)²                    # foundation: Jeffreys (single-strand) / global (AMBIG)
```

where `(μ_g, μ₊, μ₋)` is the **one imputed pie** from the neighbours (§2) and `(τ_count, τ_rna±)` are its
per-component reliabilities (§3). The bidirectional grid sum-product solves the whole chain exactly; readout is
the posterior median pie + per-node variances.

**This one solver merges and retires:**
- the separate boundary-side solve `deconv_sides` → boundaries are chain nodes;
- the `I₀` strand/count blend → the strand likelihood and the count/RNA priors compose by their *honest
  precisions* (a sharp strand likelihood dominates a wide prior; deference is emergent — no hand weight);
- the v4/2b §3 coherence guard → structural (§2);
- the ad-hoc spliced-lower-bound *weight* → the spliced floor stays as a hard constraint, but the RNA *magnitude*
  is now the derived RNA reliability, not a raw-count weight.

**Boundaries are one node** (the two physical sides are mass-accounting / edge-specific evidence, not two
beliefs); the simplex is over the boundary's unspliced crossings, with its spliced as the fixed RNA floor.

---

## 5. Sequenced after the unified solver (NOT subsumed)

- **`q_rna` — the region↔region RNA propagation coupling** (still `0.25`). The one-hop pie imputation (§2)
  carries *local magnitude*; `q_rna` carries *transitive continuity* across multiple nodes (the canonical
  innermost-AMBIG-core case, v4 §2.1). They are **complementary, not redundant** — the local imputation cannot
  reach a node whose flanking boundaries are themselves unobservable. Phase-1 A/B measured propagation helps
  (+1.1–5.0%), so `q_rna` is **sequenced after** the unified solver, with that A/B as its live arbiter.
- **The transitive gDNA chaining** (the smooth-field edge coupling — the gDNA analog of `q_rna`): likewise after.

---

## 6. Magic-number ledger (v5 state)

| constant | verdict |
|---|---|
| **`I₀` `gdna_strand_info_scale`** | **DIES** — the strand/count/RNA priors compose by honest precisions; deference emergent. (Survives only in the upstream `cleaned_gdna_count` until measured redundant.) |
| **`σ²_global`** | **KEEP MAD** — the principled between-node population spread (v4-measured; `DIRECT.predict` reverted). |
| **`q_rna` `0.25`** | **DEFERRED** (§5) — derived later as the region↔region coupling reliability; not subsumed. |
| **Jensen `Δ_k`, Bernoulli clamp, MAD `1.4826`, SCAM `k=18`/GCV grid, lattice `K=60`, `_PRIOR_EPS`** | canonical / numerical-resolution — kept, documented (per v4). |
| **the §3 coherence guard** | **RETIRED** — structural under the unified pie. |

---

## 7. Implementation sequence (one architectural step, tractable sub-phases)

The unification is one architectural step; it builds in sub-phases for bisectability, each with execution-time
code-grounding (as Phases 0/2a/2b had) *before* building. The RNA reliability fit is independent and can land
first; the imputation restructure and the node-chain unification are the substance.

- **Phase A — the RNA `var~mean` reliability fit** (straightforward, reusable; fits from current node state, no
  solver change). All-pairs, per-strand-pooled, exclude 0-RNA, Jensen, current-density axis. *Validate standalone:
  the fit is real (>18 pts), spans the exon μ-range, no extrapolation.* This is the only piece that is "just a
  fit"; it does not change the solve yet.
- **Phase B — the unified pie imputation** (the mechanism restructure). Generalize the gDNA-only
  `node_gdna_density` into a **full-pie imputation** that carries `(μ_g, μ₊, μ₋)` from one deconvolution, in
  density space with the FL conversion, respecting the spliced floor. Wire the per-component reliabilities as the
  precisions. *Validate: coherence (pie sums to 1 structurally); the gDNA slice reproduces 2a; the RNA slice
  matches Phase A.*
- **Phase C — the unified node solver** (boundaries as chain nodes; one 3-tier `ψ`; retire `deconv_sides`, `I₀`,
  the §3 guard). *Validate: complex battery ≤ 4502; net-flow non-regressing (before/after on identical BAMs,
  esp. mature-RNA FP rate); conservation `gdna+rna+spliced = total` per node; the factor-1-under-uniform eff-len
  bedrock preserved.*
- **Phase D — `q_rna` + transitive gDNA chaining** (§5), gated on the Phase-1 A/B and a measured run-interior need.

**The new must-pass gate (from the 2b critique):** the unified solver now has **two iterating priors (gDNA +
RNA) on the same simplex** — a coupling 2a did not have (gDNA↑ → f±↓ → next pass RNA predictor sees less → f±↓
more → frees f_g). Bracketing kills the within-pass loop and frozen κ/overdispersions should make it contract,
but **assert monotone convergence on a capture-on + nascent toy** before trusting; damp by under-relaxing the
RNA estimate fed to the next fit if it oscillates.

---

## 8. Load-bearing details & residual open issues (user-reviewed 2026-06-16)

1. **Spliced floor consistency (hard) — confirmed.** Projected RNA on a strand must be `≥` the node's fixed
   spliced density on that strand; the direct spliced observation dominates any imputation that implies less.
   Hard constraint in `ψ`.
2. **FL/density conversion — RESOLVED + PROVEN (2026-06-16).** A real eff-len bug was found and fixed: the
   node-pair side **density** divided the side **MASS** by `fl_mean = E[ℓ]` (the *count* effective length) instead
   of `boundary_side_eff_length = E[min(ℓ, R_side)]` (the per-side *density* length — the divisor for fractional
   crossing mass, the same `gdna_boundary_len` `priors._gdna_region_node_arrays` uses for the seam). Under uniform
   gDNA this under-read short flanks by **up to ~52%** (`mass/E[ℓ] = ρ·E[min(ℓ,R)]/E[ℓ] < ρ`), a systematic
   cross-type offset vs the correct region density `contained/E[max(0,L−ℓ)] = ρ` — inflating the residual and
   over-humbling the count. Fixed: side mass ÷ `boundary_side_eff_length` (gDNA + RNA, symmetric). **Proven by
   exact ground-truth tests** (`test_node_pair_imputation_correctness.py`): uniform gDNA ⇒ fixed `d_side = ρ` for
   all flank lengths, residual exactly 0 (uniform = nothing to learn); capture step ⇒ residual exactly
   `(ρ_hi−ρ_lo)²` (genuine structure, not a normalization artifact); RNA factor-1 recovery. The eff-len module's
   own docstrings document this `count` (`fl_mean`) vs `density` (`E[min(ℓ,R)]`) distinction; the triplet dodged
   it by differencing two same-type crossings (which also blinded it to real structure). NOTE: this fix had ~0
   battery effect (4591 vs 4590) — confirming the +88-vs-triplet is the triplet's **unjustified overconfidence on
   predictable sim data** (user-dismissed), NOT a node-pair defect; the node-pair's correctness stands
   independent of the battery. The real ship gate is **net-flow** (capture-on, where the count is genuinely less
   accurate and the node-pair's honest humility should help).
3. **Multi-prior simplex convergence (MEDIUM) — expected to converge.** The gDNA↔RNA cross-pass coupling: RNA and
   gDNA imputations may disagree locally, but with bracketing + frozen κ/overdispersions the method should
   contract to a fixed point. Prove it via the monotone-convergence test before trusting; damp (under-relax the
   estimate fed to the next fit) if it oscillates.
4. **Boundaries are INDEPENDENT NODES — the two-sided model persists but is being proven unnecessary.** The key
   development (shipped, the eff-len IPR redesign) is that each boundary is an **independent owner of mass AND
   effective length** — a truly independent node, not an attribute of its flanking regions. The two-sided
   internal representation **remains the current architecture and suffices for this implementation**; collapsing
   it to a single-node boundary is **future work**, gradually justified as we prove the two sides are not
   load-bearing for the node abstraction. (The bipartite SJ-node split remains separately DEFERRED — EM-side.)
5. **Train on current estimates (confirmed).** The RNA (and gDNA) reliability fits train on the current node
   estimates, refined over iterations (bracketed). At AMBIG sources the per-strand RNA is the current estimate —
   the same acceptance gDNA makes for imputed nodes; the reliability absorbs the extra uncertainty as wider var.
6. **Cost (eventual optimization, agreed).** Two `var~mean` fits per pass (~doubles the refit, already the
   wall-time bottleneck vs the sweep); absorbed by the eventual C++ kernel / fit caching. Profile both the refit
   and `_local_loglik` at genome scale.
7. **TSS/TES boundaries are FLUX ANCHORS, not part of the transcript (architectural principle).** A transcript is
   a linear set of nodes (regions + boundaries) bounded by a TSS and a TES; **the anchoring TSS/TES boundaries are
   NOT part of the transcript — they measure flux into/out of it.** Fragments crossing a TSS/TES are NOT
   compatible with that transcript (they belong to other transcripts' exons, nascent RNA from overlapping
   introns, or gDNA). Consequences: **(a)** the transcript's mature-RNA continuity — and therefore the RNA
   imputation edges and the deferred RNA propagation — lives *within* the transcript's node span and must **not
   cross its own TSS/TES anchors**; **(b)** the nascent-smoothness "break" at a TSS/TES is *structural* (the
   transcript ends), not merely high-variance to learn; **(c)** TSS/TES boundaries are intergenic↔exon
   transitions and do **not** own splice-junction spliced mass, so they are already naturally excluded as
   RNA-imputation sources (which require same-strand spliced) — Phase A's reliability fit is therefore largely
   unaffected; the principle binds the imputation edges (Phase B) and the deferred propagation (Phase D).

---

## 9. Validation gates (every build phase)

- **Suite** green; goldens regenerated only after the accuracy gate (the v4/2a discipline).
- **Coherence** (`f₊+f₋+f_g=1` structural; spliced floor respected); **conservation** `gdna+rna+spliced=total`.
- **Complex battery** TOTAL ≤ 4502 (the 2a baseline); anchor-starved AMBIG degrades gracefully (no vertex-push).
- **Net-flow** non-regressing on capture-on ss0.99/ss0.5 + zero-DNA (before/after on identical re-simmed BAMs) —
  **mature-RNA false-positive rate must not rise** (the RNA prior must not manufacture RNA).
- **Monotone convergence** of the two-prior simplex (the §7 must-pass).
- **eff-len bedrock**: factor-1-under-uniform-gDNA preserved (`priors.py`).

---

## 10. Doc status

- **Authoritative:** this plan (v5).
- **Superseded by v5:** `CALIBRATION_PLAN_v4.md` (the separate-imputation framing), `phase2b_rna_imputation_plan.md`
  (the RNA bolt-on + §3 guard — now structural). Keep both for the design rationale / shipped-2a record.
- **Shipped foundation (unchanged):** Phase 0 (`main@5c4627a0`) + Phase 2a (`main@8d9df6bd`).
- **Keep:** `per_node_deconv_hierarchy_design.md` §3 (the per-node math); `effective_length_redesign_plan.md` §8
  (the eff-len IPR); `splice_junction_node_architecture.md` (the deferred bipartite SJ — EM-side).
- Update `docs/README.md` to point at v5.
