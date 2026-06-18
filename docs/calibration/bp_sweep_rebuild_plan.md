# Calibration solver rebuild — identity-density belief propagation, bidirectional sweep (v3)

**Status:** v3 (2026-06-18), after two 7-lens adversarial critiques (v1: 19 BLOCKERs; v2: 19 BLOCKERs). v3
drops v2's flawed enrichment correction (tautological + blow-up), keeps the shipped factor-1 eff-lengths, uses
per-component fragment-length supports, and describes the sweep as the new machinery it actually is.

**Authority chain (first wins):** `CALIBRATION_ARCHITECTURE.md` → `calibration_nodes.md` →
`CALIBRATION_PLAN_v6.md` → `rna_imputation_transcript_structure.md` → `imputation_variance_model.md`.

**Resolved decisions (user, 2026-06-18):** (1) **identity-density** message, **defer** explicit capture
enrichment; (2) pooled-seam density divisor stays the shipped **AVERAGED** `½(E[min(ℓ,L_l)]+E[min(ℓ,L_r)])`
(the "sum of region sizes" is the *mass ceiling*, a separate statistic, not the density divisor).

---

## 0. The model in one paragraph

Calibration partitions each node's **unspliced** fragments into `{RNA+, RNA−, gDNA}`. Every node (region OR
boundary, on the linear bipartite chain) carries a pie `(f₊, f₋, f_g)` + variances, and exposes three
**densities** `(ρ₊, ρ₋, ρ_g)`. Nodes inform each other by a **pairwise, bidirectional sweep**: a source's
density is the destination's prior mean — the **plain identity** `μ_dst = ρ_src` (`ARCHITECTURE §1.2`) —
precision-weighted into the destination's own evidence. **gDNA and RNA use the same machinery**; the boundary
additionally owns the **fixed, single-strand, one-sided spliced** RNA, added to its motif strand's RNA density
(a floor on the boundary's own RNA, never a magnitude pushed onto a region). **No fraction transfer, no
enrichment factor in the message, no two-boundaries-onto-one-region blend.**

---

## 1. Node densities — per-component fragment length (the v2 fix)

Each component uses **its own FL** (gDNA crossings are gDNA-FL; RNA — nascent + spliced — is RNA-FL):

- **Region r:** `ρ_g = f_g·M_u,r / E_gdna_contained(r)`, `ρ_s = f_s·M_u,r / E_rna_contained(r)`, where
  `E_x_contained(r) = region_eff_length(L_r, x_fl) = E_x[max(0, L_r−ℓ)]`.
- **Boundary b, per SIDE σ** (the side inside the flanking region on that side):
  `ρ_g(σ) = f_g·M_u,σ / E_gdna_cross(σ)`,
  `ρ_s(σ) = (f_s·M_u,σ + spliced_s·[σ is the spliced side]) / E_rna_cross(σ)`,
  where `E_x_cross(σ) = boundary_side_eff_length(x_fl, L_σ) = E_x[min(ℓ, L_σ)]`, `M_u,σ` is that side's
  unspliced crossing mass, and `spliced_s` is the boundary's spliced mass on motif-strand s (`_spliced_floor`,
  sense+antisense summed; 0 on the non-exon side and on regions).

The boundary solves ONE pooled pie (the crossing composition is one belief), but presents the **side facing the
destination** (re-indexing identity) — resolving v1's pooled-pie-wrong-side. `f·M/E` is the definitional
density↔fraction jacobian (`ARCHITECTURE §2`), not a count vote.

---

## 2. Effective lengths — ONE owner (`effective_length.py`); shipped factor-1 kept

| quantity | formula | FL | use |
|---|---|---|---|
| `region_eff_length(L, fl)` | `E[max(0,L−ℓ)]` | per-component | region contained support |
| `boundary_side_eff_length(fl, L_σ)` | `E[min(ℓ,L_σ)]` | per-component | a side's crossing support; the spliced same-side support (RNA-FL, = v6 §4) |
| pooled-seam density support | **`½(E[min(ℓ,L_l)]+E[min(ℓ,L_r)])`** (AVERAGED — shipped, factor-1-proven, `main@3f068374`) | gDNA-FL | priors EM gDNA eff-len (UNCHANGED) |
| `fl_mean(fl)` | `E[ℓ]` | per-component | unconstrained |

**Consolidation, no behaviour change to the shipped seam:** route every eff-len call site through
`effective_length.py`; delete the ad-hoc re-derivations (my `calibrate.py` `boundary_support`, etc.). **Do NOT**
switch the pooled-seam to a sum-constrained divisor (v2's error — breaks the shipped factor-1). The
"sum-of-region-sizes" is the **mass ceiling** (`min(ℓ, L_l+L_r)` on total deposited overlap) — record it as a
distinct statistic only if a consumer needs it (none today).

---

## 3. The message — identity density, factor-1 shown

Source → one adjacent destination, per component, in density space. **Mean = plain identity** `μ_ρ,dst = ρ_src`;
destination converts to a fraction prior via its own geometry `μ_f = ρ_src · E_dst / M_u,dst` and integrates
precision-weighted into its ψ.

**Factor-1 proof (the algebra the committed fraction-transfer code demanded, `imputation.py:137`).** Under a
uniform component-c density `ρ_c` (capture aside — §6): a region's contained mass `= ρ_c·E_c[max(0,L−ℓ)]` ⇒
`ρ_c(region) = ρ_c`; a boundary side's crossing mass `= ρ_c·E_c[min(ℓ,L_σ)]` ⇒ `ρ_c(side) = ρ_c`. **Both equal
`ρ_c`** because each uses its own component-FL support — so `μ_dst = ρ_src` is unbiased. The committed code's
"density is an arithmetic error" was a *single-FL, mismatched-support* density (it divided the crossing by
`M_crossing/S` and multiplied by the region's `E/M` with a different FL); per-component, per-node supports make
the density factor-1, which the dimensionless fraction also is — but the density keeps the per-component
magnitude the fraction collapses (needed across exon↔intron biology). The spliced adds to its strand's RNA
crossing density (`spliced_s/E_rna[min(ℓ,L_exon_side)] = ρ_mat` under uniform mature) — one-sided, RNA-FL.

- **Boundary → Region:** source = the boundary's density on the side facing the region (§1). The intron-side
  has no spliced → a clean nascent/gDNA crossing density (the intron-phantom rescue). The spliced rides on the
  exon side (the boundary's owned RNA).
- **Region → Boundary:** source = the region's `ρ_s`/`ρ_g`. The boundary's spliced is FIXED; the region informs
  the boundary's **unspliced** pie only.
- **Transcript-structure gate** (`rna_imputation_transcript_structure.md`): RNA on strand s flows only where the
  strand-s bit is continuous across the edge; TSS/TES blocks RNA; gDNA flows genomically.

**Capture/enrichment is NOT in the message (v2 fix).** It is handled by each node observing its own local
density + the global prior, and within a smooth-capture region the identity density already preserves the pie
(adjacent `e` cancels). The explicit cross-capture-edge correction is a **deferred phase** (P5) — it is an open
problem (`calibration_nodes.md` flags it), and a per-node `e = ρ_g/ρ_global` correction is tautological for the
gDNA channel (it cancels the gDNA message) and blows up at RNA-only nodes.

---

## 4. The sweep — BUILD the genomic-order directional traversal (honest: it is new)

The current code (`deconv_regions_sweep`, `deconv_boundaries_sweep`) is an **all-at-once vectorized** per-node ψ
solve — there is no genomic traversal. v3 builds the unified region∪boundary chain in genomic order (per
reference; `−1` terminals are sinks) and the directional sweep:

```
init nodes (§6); build global hyperprior + var~mean models
repeat (outer iterations to convergence):
  snapshot := all node densities                              # FROZEN: previous-iteration final state
  fit gDNA + RNA var~mean on the snapshot's adjacent residuals  # §7 — frozen ⇒ no σ²→0 collapse
  L→R: each node (genomic order) ← message from its LEFT neighbour's CURRENT (Gauss-Seidel) density   # §3
  R→L: each node (reverse)       ← message from its RIGHT neighbour's CURRENT density
  recompute global hyperprior
converge when max over nodes |Δpie|_∞ between outer iterations < δ
```
- **var~mean reliability** is fit on the **frozen previous-iteration** densities (the genuine cross-node
  disagreement) — never on post-integration densities. The **message MEAN** uses the live Gauss-Seidel neighbour
  (so a strong anchor propagates along the chain in one sub-sweep — the point of a sweep). The two are coherent:
  σ²_bio(μ) is a population-level reliability curve applied per-node by its density level; both converge to the
  fixed point.
- **Directional ⇒ no recirculation:** L→R uses only left neighbours (each already carries the left chain), R→L
  only right. (v1's both-neighbour-every-pass recirculation is gone.)
- Keep `_solve_nodes`' per-node ψ core; replace the all-at-once **driver** with this per-node loop. Drop any
  "tree-exact 2-pass" claim (false for the BB ψ) — it is a Gauss-Seidel fixed-point iteration; **must-pass =
  monotone convergence** on the capture+nascent toy (`CALIBRATION_PLAN_v6 §11`); re-tune `sweep_max_passes` for
  `outer iteration = {fit, L→R, R→L, refit}`.

---

## 5. Precision — two terms (the authority form), AMBIG is the live regime

`τ_message = 1 / ( σ²_bio(μ_dst) + ρ_src/L_src )` — exactly the two terms of `imputation_variance_model.md §4`
(the learned biological/structural dispersion + the predictor's Poisson sampling noise). **Drop v2's third
"source posterior variance" term** (the source's sampling uncertainty already lives in `ρ_src/L_src`; the extra
term contradicts §4). Jacobian-convert to fraction precision `τ_f = τ_ρ·(M/E)²` for the ψ pull.

**The real lever is honest `σ²_bio`, not a scale trick.** The over-trust the memory measured (τ×0.05→4939) is
because a boundary is a *predictor* of a region (its nascent crossing ≠ the region's contained density), so the
genuine prediction-error must be LARGE there — and the var~mean, fit on the true cross-node residuals (frozen
snapshot, §7), *is* that prediction error. A high-precision single-strand destination is dominated by its own
strand by construction; an **AMBIG destination has a flat strand**, so the message governs — **D6 is a P3 GATE**
(the dry-run, §8, prints `τ_message` vs the strand Fisher per traced node; complex-loci is the measure).

---

## 6. Initialization

- **G1 (0-D seeds):** intergenic regions + intergenic↔exon boundaries → `{0,0,1}`, var 0 (sinks).
- **G2 (1-D single-strand):** intron / single-strand exon / single-strand exon↔intron boundary → strand
  deconvolution alone. Exon regions flagged **floor-anchored** (their RNA is confirmed by the boundary's spliced
  floor via the message, not imputed-magnitude). Boundary G2 = `_solve_nodes` with the de-dup one-side count +
  summed spliced floor.
- **G3 (AMBIG):** `{0,0,1}` at MAX variance.
- **Model init:** global hyperprior `ρ_global` from the **signature-known / strand-resolved gDNA nodes only**
  (carry forward `calibrate.py`'s current mechanism: MAD spread, jacobian, 1-pseudo-obs floor — avoids the
  inflated baseline). The per-node enrichment factor `e_i` is a **global/prior** concept retained for the
  deferred capture phase (P5) — it does **not** enter the message (§3). var~mean undefined until the first
  snapshot.

---

## 7. var~mean fit

Per outer iteration, over all adjacent (source→dest) density pairs, both directions, on the **frozen snapshot**:
Poisson-offset fit (`MonotoneVarMean.fit_offset`) of `σ²_bio(μ)` = excess over the computed Poisson floor. Two
models — gDNA (gDNA-FL density residual) + RNA (RNA-FL density residual, same-strand edges incl. the
spliced-bearing boundary edges). Residual axis = density (consistent with the §3 density message). Reuse
`variance_model.pair_imputation_points`/`MonotoneVarMean`.

---

## 8. Dry-run (the message is a PRIOR; the strand dominates where informative)

- **r7 (gB exon2, single − strand, strong ss):** B→R density prior ≈ ρ_RNA−; the region's own strand (217
  reads) is high-precision → f_g≈0 from the strand, the prior consistent. Print τ_message vs Fisher. Gate: f_g≈0.
- **antisense intron (in+, ss=0.65):** flat strand → message governs; the exon↔intron boundary's intron-facing
  side = clean nascent crossing density (no spliced) → ρ_RNA+>0 → f_g↓. Gate: phantom→0.
- **reg-82 (`ex+|ex−` AMBIG):** each flanking boundary's facing-side density imputed per side in its sweep
  direction (no cherry-pick, no >1 pie); f_g set by precision-weighted messages + global. Gate: complex-loci ≤5024
  and improving.
- **factor-1:** on a uniform-gDNA AND a uniform-RNA toy, assert each node's density ≈ true ρ (the bedrock).

---

## 9. Critique-blocker mapping (v1 + v2)

- mature/rule-(c)/double-count → §0/§3: density (not magnitude); spliced owned + floor; regions unspliced-only.
- `E_spliced`/read_len/v6 §4 → §1/§2: spliced = `boundary_side_eff_length(rna_fl, L_exon_side)`; per-component FL.
- fraction-vs-density contradiction w/ committed code → §3: factor-1 algebra shown (per-node, per-component
  supports make density unbiased; the committed "arithmetic error" was a single-FL mismatch).
- enrichment tautology/blow-up/identity-mean conflict (v2) → §3/§6: enrichment OUT of the message; identity
  mean; capture deferred (P5).
- pooled-seam averaged↔sum (v2) → §2: keep shipped AVERAGED; sum = mass ceiling, not the divisor.
- per-component FL (v2) → §1/§2.
- directional sweep doesn't exist / ψ recirculation / var~mean collapse → §4: BUILD the traversal; frozen-fit +
  live-mean; directional (no recirculation).
- τ third term (v2) → §5: two-term τ per imputation_variance_model §4.
- precision-scale deferred-but-claimed-safe → §5: honest σ²_bio is the lever; AMBIG is a P3 gate.
- P0 revert target = HEAD (v2) → §10: corrected.
- per-node state vs ARCHITECTURE §6 stateless / v6 §9(i) → §10-note: this rebuild IS Step 3 (the first-class
  bipartite node); record the supersede of §8's deferral.

---

## 10. Phasing + open items

**P0** — `git checkout 60043391 -- src/rigel/calibration/{calibrate,imputation,boundary_nodes}.py` (HEAD is
`60043391`; the v1 drift is the UNCOMMITTED working-tree edits; post-revert = gDNA-imputation-from-deconv_sides,
NO RNA imputation, NO boundary co-solve). Strip the TEMP trace. Build the unified region∪boundary genomic chain
+ `NodeBelief` densities; **consolidate `effective_length.py`** (no shipped-seam change). *Gate: suite green.*
**P1** — init §6. *Gate: zero-gDNA introns init ~0 via strand; intergenic locked.*
**P2** — gDNA-only directional sweep (§3 gDNA identity density, §4, §5). *Gate: monotone convergence; zero-gDNA
clean; complex-loci ≤5024; net-flow non-regressing; factor-1 on uniform-gDNA toy.*
**P3** — RNA density message + spliced (§3 RNA). *Gate: intron phantom→0; r7=0; AMBIG via precision (τ-vs-Fisher
printed); complex-loci improves; factor-1 on uniform-RNA toy.*
**P4** — per-node schema + priors/derive rewire; retire `deconv_sides`/I0/`gdna_deconv_quantile` in dependency
order (`ARCHITECTURE §6.4`). Golden update LAST.
**P5 (deferred)** — explicit capture-enrichment correction across capture edges (the open problem); only if a
measured capture regression warrants it.

**Open items:** (D-schema) per-node first-class (v6 §10 REQUIRED) vs a transitional region-keyed projection for
P3 — specify the crossover if transitional. (D-passes) outer-iteration δ + `sweep_max_passes` re-tune. (D-floor-
anchor) the exon-region floor mechanics in G2 (how the boundary spliced floor reaches the region's f_s).
