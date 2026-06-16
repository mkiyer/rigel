# Calibration — the inverted-model unified-solver redesign (PLAN v4, EXECUTION-READY)

> **⚠️ SUPERSEDED by `CALIBRATION_PLAN_v5.md` (2026-06-16, same day)** for the imputation architecture: v5 unifies
> the per-component imputations (gDNA 2a + the RNA bolt-on) into **one full-pie message + a per-component
> reliability layer** = the unified node solver (one step), retiring the §3 coherence guard as structural. v4's
> **shipped** content (Phase 0 teardown, Phase 2a gDNA var~mean, the σ²_global=MAD reversal, the ledger) remains
> the factual record; its *forward* plan (separate 2b, the §3 guard) is superseded by v5.

**Status:** SUPERSEDED execution plan (2026-06-16). Supersedes `CALIBRATION_PLAN_v3.md`. v3 formalized the
completion of the *region-only* iterative bootstrap; v4 incorporates two corrections from design review that
change the spine of the algorithm — the **gDNA-chains / RNA-local inversion** and the **unified node solver** —
and re-scopes the deferrals into an ordered, measurement-gated implementation plan. Every `file:line` below was
cross-checked against the live code (`main@3f068374`). Companion docs kept: `per_node_deconv_hierarchy_design.md`
(per-node math), `effective_length_redesign_plan.md` §8 (eff-len IPR, SHIPPED), `variance_model.py` docstring
(the SCAM fitter). Where this and any older note disagree, this wins.

Production bar (the acceptance test for every decision): **ROBUST** (no brittle tunable — derived, canonical, or
a documented numerical-resolution knob), **ELEGANT** (behaviour emerges from honest precisions), **EFFICIENT**
(genome-scale), **ACCURATE** (optimal given identifiability).

---

## Amendment A1 — Phase-1-driven refinement (2026-06-16)

Phase-1 measurement + design review refine the var~mean and the propagation model (updates §7, §2, the ledger,
and Phase 2). The §0/§2 "gDNA chains / RNA local" framing below is **superseded** by points (1)-(2) here.

**(1) The var~mean is the node-PAIR imputation-error model, fit on ALL adjacent node pairs — NOT
boundary↔region↔boundary triplets restricted to count-observable nodes.** Phase-1 confirmed the current
triplet/observable-only fit flat-extrapolates at **71% of captured-exon nodes**: count-observable ⇔ no-exon-bit
⇔ intron/intergenic ⇔ the depleted low-μ regime, so the fit never sees the enriched exon μ. The fix: for every
sequential adjacent pair `(source, dest)`, fit the imputation residual `dest_estimate − impute(source_estimate)`
as a `var~mean` curve over **all** nodes, using the iteration's **current** gDNA/RNA estimate at every node
(captured exons included) as the "actual." Every node contributes its current-μ residual ⇒ the fit spans the full
μ range ⇒ it **interpolates** at exons ⇒ **no extrapolation, and the confidence(μ) guard is DISSOLVED** (it was a
band-aid for the wrong model — retracted, not deferred). DIRECT (own-count Poisson) stays for observable nodes;
the **IMPUTATION** model becomes this all-pairs version. Built for BOTH gDNA and RNA.

**(2) RNA propagation STAYS, with a DERIVED coupling reliability — `q_rna=0.25` dies as a *derived* quantity, not
by removing RNA chaining.** Phase-1 A/B showed RNA-odds propagation helps (+1.1–5.0% on the complex battery). RNA
carries cross-node correlation and should propagate, but the coupling strength must be derived: the **RNA
pair-model `var~mean`** supplies the per-μ coupling reliability on the (annotation-gated) same-strand-exon edges.
The hard "gDNA chains / RNA local" rule is SUBSUMED by: **each component propagates on its structurally-gated edges
with its own derived reliability** — gDNA's high where genomically smooth; RNA's reflecting expression coherence
(low ≈ local where isoforms/capture break continuity). Behaviour emerges from the learned reliabilities, not a
hard rule or a constant. *Caveat:* the battery is the favorable case (capture-off, single-abundance transcripts);
a multi-isoform/capture-on stress locus is still needed to test where RNA coupling should decouple.

**(3) Phase 2 is reframed** around building this pair-model `var~mean` (gDNA + RNA) as the foundational step — it
delivers `τ_count` (no extrapolation), the derived `q_rna`, and `σ²_global`, with the Jensen df-offset + Bernoulli
clamp applied to it. The all-pairs fit needs an estimate at every node, so it composes naturally with the unified
node solver (Phase 5); it can be prototyped on the current region + boundary-side estimates.

**A1.1 — Phase-2a result (2026-06-16, SHIPPED to working tree, suite-pending-goldens):** the IMPUTATION-axis fix
(fit & query on the region's CURRENT density), the Jensen df-offset, and the Bernoulli clamp landed and
**eliminated the extrapolation** (0% of captured-exon nodes extrapolate, was 71%) while **improving** the complex
battery (4595 → 4502). `sweep_max_passes` 4→6; `1e-3` stop promoted to `CalibrationConfig.sweep_convergence_delta`.
**`σ²_global = DIRECT.predict(ρ_global)` was tried and REVERTED** — A/B-measured ~5.5× too tight, +11.3% battery
(see ledger). Kept the MAD between-node spread. Downstream: 4 golden scenarios shift ~4% (the intended
prior change); golden regen gated on a net-flow accuracy confirmation (release-suite BAMs need a rebuild).

---

## 0. TL;DR — the two corrections and the elegance

Calibration deconvolves each node's unspliced mass into the 2-simplex pie `(f₊, f₋, f_g)` = sense-RNA /
antisense-RNA / gDNA (RNA-vs-gDNA **only**; the EM separates mature/nascent). `f_g` feeds the per-locus EM prior;
the deconvolved gDNA mass also drives the effective-length IPR contraction — so **gDNA accuracy is doubly
load-bearing**.

**Correction 1 — gDNA CHAINS, RNA stays LOCAL** (inverts v3). gDNA is the only genomically transitively-continuous
quantity (uniform modulo CNV/capture), so it is what propagates along the genomic chain. RNA is spiky,
isoform-driven, capture-variable — it does **not** chain (isoform structure and arbitrary probe placement break any
spatial smoothness; v3's RNA-odds propagation rested on a rejected premise). A node's RNA is `(1−f_g)` split
**locally** by the strand tilt + local spliced evidence. **AMBIG resolution = the smooth gDNA field pins `f_g`; the
local tilt splits the rest** — propagation flows inward from solvable outer anchors (intergenic `f_g=1`, introns) to
the innermost ambiguous exon↔exon boundary.

**Correction 2 — the unified node solver.** Regions and boundaries are the *same* node, solved by the *same* 3-tier
posterior. A boundary **owns** its spliced fragments; a region **imputes** spliced one-hop from its boundaries —
that is the *only* difference (the user's review-point 8). This eliminates the separate `deconv_sides`, the `I₀`
strand/count blend, and the side-freezing. All nodes iterate; only the capture-invariant direction quantities
(κ, the two strand overdispersions) stay frozen.

**The elegance, made precise:** one `MonotoneVarMean` machine, refit each pass on the running gDNA estimate, emits
**every** precision. The two magic constants die: **`I₀` (`gdna_strand_info_scale=10.0`)** because the count enters
at its honest prediction-error precision (no `(1−w)` down-weight — strand deference *emerges* from the Fisher info
`I=N(2κ−1)²` overwhelming a wide count Gaussian); **`q_rna` (`0.25`)** because RNA does not chain. The enrichment-
ratio circularity (`local/global → 1.0`) is dissolved by a **ratio-free** neighbour `var~mean`, not per-node ratios.

---

## 1. First principles (carried from v3, re-grounded)

1. **Strand = unbiased DIRECTION, priority.** Fisher info `I = N·(2κ−1)²`; Beta-Binomial overdispersion; capture-invariant.
2. **Count = local MAGNITUDE, fallback.** Down-weight by reliability, never bias-correct; the honest prediction-error
   reliability *converts* unmodelable capture bias into measured variance.
3. **Global baseline = FOUNDATION** for truly-unanchored nodes; never overtakes strand or count.
4. **The pie caps capture smear:** a *fraction* projected onto a node's observed counts is bounded by local evidence.
5. **Behaviour emerges from honest precisions, not hand weights** — the reason `I₀`/`q_rna` can be deleted.
6. **Capture-awareness is partitioned, never duplicated:** DIRECTION → Tier-1 strand; LOCAL MAGNITUDE → Tier-2 count
   (crossing flux) + the eff-len IPR; BASELINE → Tier-3 global (structure-free).
7. **Continuity is the propagation gate:** gDNA is genomically smooth → it chains (edge-preserving, decoupling at
   capture/CNV jumps); RNA has no spatial continuity → it stays local. (This is the corrected statement of v3's
   principle 7.)

---

## 2. Architecture (final)

**SETUP (once, frozen).** κ (`rna_sense_frac`); the two strand Beta-Binomial overdispersions (gDNA mean-½, RNA
mean-κ); `strand_deconvolve` of the boundary crossings; the count crossing-flux imputation `region_count_frac`
(Tier-2 capture-awareness). Fit once; never iterate (re-fitting would feed the deconvolution's own gDNA weight back
into the strand model — a feedback loop with no upside).

**INIT — all-gDNA.** `gdna_c = u_tot` (`calibrate.py:238`): the conservative worst case; makes `ρ_global` a
deliberate over-estimate and gives the `var~mean` support across the full depleted→enriched range from pass 0.

**THE ITERATIVE LOOP** (≤ `sweep_max_passes`; stop on `sweep_convergence_delta`). Each pass re-fits the
*magnitude/baseline* quantities on the current gDNA estimate: `ρ_global` (count-observable nodes only,
`calibrate.py:247-249`) and the gDNA `var~mean` (DIRECT + IMPUTATION, Jensen-corrected); derives the precisions
(`τ_count` per node; `σ²_global = DIRECT.predict(ρ_global)·jac²`); then solves every node. Frozen across passes: κ,
the strand overdispersions.

**THE UNIFIED PER-NODE 3-TIER SOLVE** (exact grid sum-product on the 2-simplex K-lattice; the engine
`_sweep_chain`, `simplex_sweep.py:143-188`, is unchanged). Every node — region OR boundary — carries the same `ψ`:

```
ψ = strand_mixture(u₊, u₋ | f, κ, od_g, od_r)   # Tier 1 — capture-invariant DIRECTION (both kinds)
  + sided_spliced_lower_bound                    # boundary OWNS its spliced; region IMPUTES one-hop
  − ½·τ_count·(f_g − μ_count)²                    # Tier 2 — MAGNITUDE, RAW precision (no (1−w); I₀ dead)
  − ½·τ_node ·(f_g − μ_node )²                    # Tier 3 — Jeffreys (single-strand) / global Gaussian (AMBIG/NONE)
```

gDNA continuity flows along the chain (interim: one-hop boundary→region imputation + the global tier; deferred:
the transitive edge-preserving coupling, §Phase 7). The local strand tilt `p = ½f_g + κf₊ + (1−κ)f₋` splits
`(1−f_g)` into `(f₊, f₋)`. Readout: posterior **median** `f_g` (avoids the overdispersion skew) + posterior
**variance** (`_fg_var`, surfaced).

**OUTPUTS.** Per-node `{gdna_mass, rna_mass, gdna_frac, gdna_frac_var}`. The surfaced variance reaches
`CalibrationResult` and (gated) precision-weights the per-locus Dirichlet split *direction*. Boundaries become
**freestanding prior nodes** (own mass, own effective length, own genomic span).

### 2.1 The worked example (your 3-transcript locus, verified against `regions.py`/`signature.py`)

`T1+` exons (1000,2000),(10000,13000),(15000,16000); `T2+` single exon (1000,16000); `T3−` exons
(5000,6000),(11000,12000),(20000,21000). The region partition (13 regions incl. the two intergenic flanks):

```
 [..1000) NONE  ┐ anchor (f_g=1)                                      observable anchors ┌ NEG  NEG  NONE
   intergenic   │                                                                        │ intron exon  intergenic
  ─B(obs)─[1000,2000)─[2000,5000)─[5000,6000)─…(AMBIG run)…─[15000,16000)─B(obs)─[16000,20000)─B─[20000,21000)─B─[..)
              POS        POS       AMBIG ───────────────────────── AMBIG          NEG          NEG
              exon+      exon+     (T1+ × T3−, core exon↔exon @ [11000,12000))    intron−      exon−
              └──────────────── ALL exon+ (T2 single exon spans 1000–16000) ──────────────┘
```

**The hard structure your example creates:** because `T2+` is one exon spanning 1000–16000, **the entire gene body
[1000,16000] is exonic ⇒ no interior region or boundary is count-observable.** The only gDNA anchors are the
*flanks*: the intergenic↔exon boundary at 1000, the exon↔intron boundary at 16000, the intron region [16000,20000],
and the boundaries at 20000/21000. So the whole POS+AMBIG span is solvable *only* by gDNA flowing **inward from the
two ends** — exactly your point-5 mental model. At the innermost exon↔exon core [11000,12000]: the gDNA field
(interim: the global `ρ_global`) pins `f_g`; the local strand tilt splits `(1−f_g)` into `f₊` (T1) and `f₋` (T3).
This is the canonical case the transitive chain (Phase 7) targets and the global fallback covers in the interim.

---

## 3. Phased implementation plan

**Ordering constraint:** Phases 0–4 are executable **without** the deferred transitive chaining. Phase 5 (the
unified node solver) is the core. Phase 7 (transitive gDNA chaining) is gated on a measured insufficiency of the
interim. The "when to remove RNA-odds" question is resolved by the blocking Phase-1 A/B.

### Phase 0 — Truth-in-docs + dead-code teardown (zero behaviour change)

- Fix the "acyclic / single feed-forward / no EM loop" docstrings → "iterative bootstrap": `config.py:226-232`,
  `calibrate.py:1,4,83`, `result.py:40`(+2), `substrate.py:1`, `__init__.py:1`, `priors.py:188`, `derive.py:1`.
  Fix `calibrate.py:4` axes (sense/antisense/gDNA, **not** mature/nascent).
- Delete `count_trust_beta` (`simplex_sweep.py:48,89-90,193,260`).
- Delete the dead LOESS path in `density_model.py` (`_LOESS_SPAN/_MIN_FIT/_BISQUARE_C:97-99`, `_loess`,
  `_node_disagreement`, `density_variance_curve`, `_count_fraction_variance`, the `need_count_variance` branch,
  `NodeDensity.count_gdna_frac_var/count_rel_var`). *Gated:* production keeps `gdna_deconv_quantile=½` (the value is
  multiplied by `Φ⁻¹(½)=0`); if the QC quantile knob is retained, substitute a minimal own-count `1/N` variance.
- Delete `simplex.solve_node` + `SimplexSolution`/`SignatureInit`/`init_from_signature`/`region_adjacency`
  (`simplex.py:51-114,154-247`). **Keep** `_simplex_lattice` + `_mixture_strand_loglik`.
- Delete `rna_variance.py` + `mature_density.py` (no production importer; the q_rna source is removed in Phase 2).
- Delete dead tests (`test_simplex.py`, `test_rna_variance.py`, `test_mature_density.py`); prune
  `test_spliced_boundary_onesidedness.py` to the `splice_junction` refs only.
- Delete the broken `scripts/debug/dissect_efflen.py` (imports deleted `_transport_boundary_flux`); `git clean` orphaned `.pyc`.
- **Gate:** full suite green (1122 baseline); zero production-path change; no golden regen.

### Phase 1 — Instrument before building (BLOCKING measurement)

Settles the decisions that govern Phases 2/5/6/7.
- `complex_loci_benchmark.py`: add an A/B flag — current sweep with RNA-odds ON (`q_rna=0.25`) vs OFF (all edges
  decoupled). **This is the single arbiter of "when to remove RNA-odds."** (The vindicating TSVs are stale post-teardown — re-measure.)
- Instrumentation script: distributions of `τ_count`, the density→fraction Jacobian² (the v3 `geom2`), `var_d`, and
  the realized count-prior influence `τ_count·jac²`; per-pass mass-delta + prediction-error trace; count + μ-coverage
  of (a) all-three-observable triplets and (b) adjacent count-observable region pairs (the Phase-7 `s_edge` support);
  per-pass `var~mean` refit + `_local_loglik` wall-time and peak per-chain edge memory; whether `ρ_global` sits inside
  the DIRECT fit range `[x_lo,x_hi]`.
- **Decisions emitted:** (1) RNA-odds ON vs OFF on the complex battery; (2) does `τ_count` survive the Jacobian;
  (3) does the loop converge in ≤4 passes without oscillation; (4) do triplets/pairs span the exon μ-range;
  (5) is the refit sub-second and `ρ_global` in-range.

### Phase 2 — Cheap unconditional correctness (the robustness derives)

- **σ²_global** (`calibrate.py:268-274`): delete the MAD block; `sig2_glob = max(direct.predict([ρ_global])[0]·jac², 1e-12)`
  — the zero-DNA pin. *Gated* on Phase-1 (5) DIRECT support near `ρ_global`.
- **Jensen df-offset** (`variance_model.py`): add `kcount` to `VarMeanPoints` and apply
  `Δ_k = log((k−1)/2) − ψ((k−1)/2)` (positive; `ψ`=digamma) to `log(raw_var)` inside `MonotoneVarMean.fit`, before
  both the spline and the power-law fallback. Removes the verified +3.56×/+1.78× over-confidence at k=2/3.
- **Bernoulli clamp** (`calibrate.py:267`): `sig2_frac = min(var_d·jac², count_frac·(1−count_frac))`.
- **Promote the stop**: `CalibrationConfig.sweep_convergence_delta = 1e-3` (golden-bit-compatible); replace the inline
  literal at `calibrate.py:291`; log per-pass `mean|Δf_g|` + the prediction-error trace.
- **q_rna removal** (default = remove; deferred to Phase 7 only if Phase-1 (1) says ON helps): delete `q_rna`,
  `_log_odds` (`simplex_sweep.py:39-44`), the same-strand-exon edge build (`:264-281`). All edges become `Q=∞` ⇒ the
  `O(P)` decoupled short-circuit (`:170-188`) fires everywhere — a pure perf win.
- **Gate:** zero-DNA phantom tol TIGHTENS; gDNA-release net-flow ≥ baseline; complex battery TOTAL ≤ baseline
  (anchor-starved AMBIG degrades toward `ρ_global`, not vertex-push); full suite green; unit-test `Δ_k` at k=2/3/∞.

### Phase 3 — Surface per-node uncertainty (unconditional)

- `result.py`: add `gdna_frac_var_{contained,left,right}` (float64[R]); register in `__post_init__`. (Verified:
  `deconv_sides` returns `NodeDeconv` which already carries `gdna_frac_var` — no struct extension needed.)
- `calibrate.py`: capture `regions.gdna_frac_var` + the side variances; wire into the `CalibrationResult` build.
- **Gate:** schema test (finite, non-negative); goldens regen for field additions; no behaviour change.

### Phase 4 — Boundaries as freestanding prior nodes (unconditional correctness)

- `priors.py`: add `_gdna_nodes_table(...)` → `R + n_internal_seams` rows with
  `(node_start, node_end, node_ref, mass, partic=m²/S, support=S, rna, contained_ev, fgvar)`. Region rows: span
  `[start,end)`, mass `mass_gdna_contained[r]`, support `gdna_region_eff_len[r]`. Seam rows (internal-seam mask
  `ref_id[:-1]==ref_id[1:]`): span `[end[r],start[r+1])`, mass `mass_gdna_right[r]+mass_gdna_left[r+1]`, support
  `½(gdna_boundary_len[r]+gdna_boundary_len[r+1])`, `contained_ev=0`.
- Generalize `_project_regions_to_loci` → `_project_nodes_to_loci`; add a **point-overlap rule** for 0-width seams
  (the common case): assign to each locus block with `start≤p<end` (the area-overlap rule yields 0 share for
  0-width spans and would *drop* the seam mass — load-bearing for conservation).
- `assemble_priors` calls the node-table builder + node projector; per-locus sums are byte-identical whenever every
  node of a region chain falls in one locus (only locus-straddling seams move). The Laplace IPR closed form is unchanged.
- **Keep** `_gdna_region_node_arrays` as-is for `capture_eff_length` (transcript footprint is single-locus; its
  region-keyed seam fold is correct there). `capture_eff_length.py` is **unchanged**.
- **Gate:** `test_priors`/`test_capture_eff_length` green at rtol=1e-9; NEW: a seam straddling two loci splits by its
  own span; `Σ gdna_prior_count == Σ deconvolved gDNA mass`; 0-width point-overlap unit test.

### Phase 5 — The unified node solver (CORE; executable on the interim, no transitive chaining)

- `simplex_sweep.py`: `deconv_regions_sweep` → `deconv_nodes_sweep`. Add `_build_node_chain(...)`: per ref, interleave
  region + internal-boundary indices in genomic order (2k−1 nodes; terminal boundaries dropped). Pool a boundary's
  inputs: `u_{pos,neg} = substrate.right.n_unspliced_{pos,neg}[r] + substrate.left.n_unspliced_{pos,neg}[r+1]` (the
  two physical sides of one crossing = **one belief** — verified `substrate.py:105-114`); spliced = pooled crossing
  spliced, oriented via `_side_strand_orientation`.
- Build ONE `psi` over all chain nodes via `_local_loglik`. Region spliced lower bound = the **max** over its two
  flanking boundaries' owned spliced bound (the spliced count carries its own honest Poisson precision; no RNA
  reliability built). Boundary spliced = its own pooled crossing spliced.
- **Delete the `(1−w_strand)` count down-weight** (`simplex_sweep.py:235-240`): pass `count_precision` **raw**. `I₀`
  for the solve dies here; strand deference emerges from `I=N(2κ−1)²` overwhelming the wide count Gaussian.
- Edges are region↔boundary; with q_rna removed all decouple → `_sweep_chain` runs on the longer node list via the
  `O(P)` short-circuit. Return region AND boundary beliefs.
- `calibrate.py`: **delete** the `deconv_sides` call (`:218-229`) and `gdna_left/right`. The `var~mean` IMPUTATION
  reference becomes the **previous pass's** boundary-node `f_g·mass` (bracketed, never in-loop — preserves convergence).
- `result.py`: add `mass_gdna_boundary`, `mass_rna_boundary`, `gdna_frac_var_boundary` (keyed to left-flank r, 0 at
  terminal/cross-ref); drop `mass_gdna_left/right` + `mass_rna_left/right` (subsumed). `priors._gdna_nodes_table` reads
  `mass_gdna_boundary[r]` directly.
- `strand_deconv.py`: delete `deconv_sides`, `_deconv_per_node`, `_compute_side`/`_SideQuantities`/`_left_right_neighbors`.
  **Keep** `strand_deconvolve`, `cleaned_gdna_count`, `boundary_side_seeds`, `_side_strand_orientation`.
- `config.py`: `gdna_strand_info_scale` is no longer read by the sweep; retained ONLY for `cleaned_gdna_count`
  (a separate upstream mechanism) until measured redundant. Re-document.
- **Gate:** complex battery AMBIG net-flow ≤ baseline (the exact failure the `(1−w)` blend prevented); per-pass
  trace converges without oscillation; node conservation `gdna+rna=total` with the pooled node; the factor-1-under-
  uniform bedrock (`priors.py:120-131`) preserved; full suite green; goldens regen for the schema change.
- **Phase-1-dependent:** gated on (2) `τ_count` survives the Jacobian and (3) the loop converges. If raw `τ_count`
  over-trusts at AMBIG, sequence the Phase-6 IMPUTATION re-cut + Jensen first.

### Phase 6 — IMPUTATION re-cut + diagnostics (gated)

- `variance_model.py`: add `imputation_pairs(...)` over (observable-boundary-side, count-observable-region)
  adjacencies — predictor `x = cleaned_crossing/fl_mean` (**frozen**, capture-invariant), held-out truth
  `raw_var = (gdna_c[r]/eff_len[r] − x)²`. Re-cut `fit_imputation_varmean` to consume these (the user's
  predicted-vs-actual model) instead of `~region_observable`. Keep `fit_direct_varmean` own-count, Jensen-corrected.
- Diagnostics feature: `CalibrationConfig.emit_diagnostics` / `--calibration-diagnostics` → `varmean_gdna`,
  `per_region`/`per_boundary`, `global_scalars` dataframes (feather+TSV). No matplotlib in `rigel quant`.
- **Gate:** ss0.50 + capture-on net-flow improve; complex battery non-regressing. If triplets/pairs are scarce or
  low-μ-only (Phase-1 (4)), keep `raw_var~mean` + Jensen and record the decision.

### Phase 7 — Transitive gDNA chaining + prior confidence blend (LATER, gated)

The principled replacement for the interim one-hop+global fallback — only if the interim is measured insufficient.
**"When to remove RNA-odds":** removed in Phase 2 if Phase-1 (1) shows OFF ≥ ON; otherwise kept as a documented
crutch through Phases 2–6 and removed here, when the gDNA chain matches-or-beats it on the same battery (the
couplings are mutually exclusive per edge; `_sweep_chain` is coupling-agnostic ⇒ a localized matrix-builder swap).
- `simplex_sweep.py`: `_gdna_logdensity(f_g, mass, eff_len)`; `_edge_loggdna(ld_i, ld_j, s_edge)` = a **redescending
  Cauchy** potential `−log(1 + (Δlog-gDNA-density / s_edge)²)` (edge-preserving — smooths within capture-uniform
  stretches, decouples at jumps; the saturated edge still triggers the `O(P)` short-circuit).
- `variance_model.py`: `fit_neighbour_varmean` over genomic-adjacent BOTH-observable region pairs:
  `(mean=½(d_i+d_{i+1}), raw_var=¼(d_i−d_{i+1})²)`. **Ratio-free** — genuine neighbour disagreement with strictly
  fewer DOF than the per-node densities, which is *precisely* why it does not reintroduce the `local/global → 1.0`
  tautology. `s_edge = sqrt(var_nbr(μ_edge))/μ_edge`.
- **Prior confidence blend** (`CalibrationConfig.priors_confidence_blend`, default False = bit-identical): per-locus
  `c` = inverse-variance, mass-weighted pooling of `node_fgvar` (floor `var_ref = (1/K)²/12`, the lattice-cell
  variance); `a_g' = c·a_g + (1−c)·a_g_global`, `a_r'` likewise, blending toward the **global split** (not zero),
  preserving depth `a_g+a_r` exactly (direction-only); discount prior-pinned (AMBIG/NONE) nodes.
- **Gate:** complex battery AMBIG `|f_g−oracle|` ≤ the RNA-odds baseline; net-flow leak ≤ the one-hop+global baseline
  (no capture-edge smear); zero-DNA no phantom; `O(P)` short-circuit preserved; confidence blend flipped only if
  ss0.50/AMBIG net-flow improves and zero-DNA does not balloon; assert `a_g'+a_r' == a_g+a_r`.
- **Phase-1-dependent:** gated on (4) neighbour-pair μ-coverage AND a measured run-interior failure of the interim.
  If one-hop+global already matches the RNA-odds baseline, the chain is **deferred entirely** (honest minimal ship).

### Perf PR (separate, byte-for-byte gated) — raise K without cost (your point 3)

`_sweep_chain` is already cached + `O(P)`-short-circuited; the cost is `_local_loglik` materializing ~9 dense `(R,P)`
arrays for the strand mixture (~1.5 GB / ~2 s at R=100k, K=60, ×4 passes). Port `_local_loglik` to a fused nanobind
kernel (`src/rigel/native/calibration/local_loglik.cpp`) that emits `psi` + median + variance without the
intermediates (~5–15× wall, >5× memory headroom; lets K rise to ~100). With q_rna removed all edges decouple, so the
readout is a per-node normalize the kernel can emit directly. Mirror a Python reference byte-for-byte
(`tests/native/_local_loglik_reference.py`). This is perf-only, not a correctness phase.

---

## 4. Definitive magic-number ledger

| Constant | Current | Verdict | Derived replacement |
|---|---|---|---|
| **σ²_global** | `1.4826·MAD²·jac²` (`calibrate.py:269-274`) | **KEEP-MAD** (Phase-2a REVERTED the DERIVE) | `direct.predict(ρ_global)` was measured ~5.5× too tight (it's one node's *sampling* variance, not the between-node *spread* an unanchored node faces) ⇒ over-pinned anchorless AMBIG ⇒ **+11.3% complex battery** (A/B-isolated). The MAD IS the principled population spread (1.4826 = normal-consistency, not a tunable); → 0 at zero-DNA = the pin. |
| **`I₀` `gdna_strand_info_scale=10.0`** | `config.py:288`; sweep down-weight `simplex_sweep.py:236`; side blend; side-clean | **DELETE** (sweep+side, P5) | **DIES** — count enters at raw `τ_count`; deference emerges from `I=N(2κ−1)²`. Survives only in `cleaned_gdna_count` until measured redundant |
| **`q_rna` `0.25`** | `simplex_sweep.py:194,264-281` | **DERIVE** (Amendment A1; P1 A/B: RNA-odds helps) | **DIES as a constant, not as a mechanism** — replaced by the per-μ RNA pair-model imputation-error reliability on same-strand-exon edges; RNA propagation STAYS |
| **Jensen df-offset** | absent (`variance_model.py:282`) | **DERIVE-NOW** (P2) | `Δ_k = log((k−1)/2) − ψ((k−1)/2)`; +3.56×@k=2, +1.78×@k=3 |
| **`sig2_frac` upper bound** | uncapped (`calibrate.py:267`) | **DERIVE-NOW** (P2) | Bernoulli `count_frac·(1−count_frac)` |
| **convergence stop** | inline `1e-3` (`calibrate.py:291`) | **PROMOTE** (P2) | `sweep_convergence_delta` + prediction-error trace |
| **per-node enrichment factor** | proposed `local/global` | **REJECT** | circular identity; dissolved by ratio-free neighbour `var~mean` (P7) |
| **`s_edge` (gDNA chain scale)** | new (P7) | **DERIVE** | `sqrt(var_nbr(μ))/μ` from `fit_neighbour_varmean` |
| **`var_ref` (confidence floor)** | new (P7) | **DERIVE** | `(1/K)²/12` — min representable posterior variance |
| **`count_trust_beta`** | `simplex_sweep.py:48,89-90,193` | **DELETE** (P0) | dead at default |
| **dead LOESS** `_LOESS_*` | `density_model.py:97-99` | **DELETE** (P0) | unreachable at quantile ½ |
| **`τ_global` floor `1.0`** | `calibrate.py:276` | **KEEP-canonical** | 1-pseudo-observation foundation; document |
| **Tukey `4.685`, MAD `1.4826`, Beta(2,2) od=0.2** | various | **KEEP-canonical** | 95%-eff Tukey; normal-consistent; hard model bound |
| **`_PRIOR_EPS=1e-3`** | `simplex_sweep.py:35` | **KEEP-numerical-knob** | Jeffreys edge floor; document |
| **`sweep_n_grid K=60`** | `config.py:296` | **KEEP-numerical-knob** | lattice resolution; raise via the C++ kernel; document K=20/60/100 sweep |
| **SCAM `k=18`, GCV λ-grid** | `variance_model.py` | **KEEP-numerical-knob** | basis upper bound; GCV-λ smooths |
| **`n_grid=200`** | `config.py:254` | **KEEP-numerical-knob** | 1-D side posterior in `cleaned_gdna_count` |
| **`POOL_EB_PRIOR_ESS=1000`** | `fl.py:59` | **FLAG** (separate) | over-shrinks tiny panels; scale or config |

**Confirmed: `I₀` and `q_rna` both die.**

---

## 5. Cross-cutting resolutions (the subtle interactions)

1. **Sides freeze (v3) vs all-nodes-iterate (v4).** The boundary node's *belief* iterates; the `var~mean` IMPUTATION
   *reference* for a pass uses the **previous pass's** boundary-node `f_g·mass` (bracketed, never in-loop — the
   stability v3 §12 required). Sound because κ + overdispersions are frozen. Phase-1 trace must confirm no oscillation.
2. **Fold-boundaries vs transitive chaining — orthogonal & sequenced.** Phase 5 makes boundaries *nodes* for the
   local solve + prior; all edges decoupled. Phase 7 adds the *coupling* on the already-unified `R—B—R—B` chain (the
   boundary's local crossing — gDNA+nascent, mature-free — is the sharp anchor competing with the smooth field).
3. **Prediction-error needs a stable reference.** DIRECT (own-count) needs none; IMPUTATION's *predictor* is the
   frozen cleaned-crossing density (capture-invariant), only the *destination* iterates.
4. **Removing `deconv_sides` does not lose the `var~mean` training set** — its boundary anchor re-points to the
   (bracketed) boundary-node estimate; the cleaned crossing density for Phase-6 pairs comes from `cleaned_gdna_count`
   (survives, upstream).
5. **No strand double-count at the unified boundary node.** Its Tier-2 `μ_count` uses the **RAW** pooled crossing
   density (strand enters once via Tier-1); the count-module *region* imputation keeps using the cleaned crossings
   (a separate consumer).
6. **`cleaned_gdna_count`'s `I₀` is a distinct mechanism** (count cleaning, not the solve blend); it keeps
   `info_scale` until a measured-redundancy check. The field stays in config but is no longer read by the sweep.
7. **`SideDeconv` variance worry is moot** — `deconv_sides` returns `NodeDeconv`, which carries `gdna_frac_var`.
8. **`capture_eff_length` is unchanged** — `_gdna_region_node_arrays` is kept for the transcript footprint;
   `_gdna_nodes_table` is additive, used only by `assemble_priors`.

---

## 6. Residual measurement-gated open questions

1. **RESOLVED (Phase 1): the count lever SURVIVES the Jacobian.** `geom2` spans ~6 orders under capture but `var_d`
   shrinks in lockstep ⇒ `sig2_frac` stays O(1) ⇒ count-inert = 0% everywhere; count is a live lever on 84–87% of
   nodes under capture. Jensen + the IMPUTATION re-cut ARE visible to net-flow.
2. **RESOLVED (Phase 1 → Amendment A1): the OBSERVABLE-only fit does NOT span the exon μ-range (71% extrapolation).**
   Hence the switch to the all-pairs node-pair var~mean (fit on every node's current estimate), which DOES span it.
3. **RESOLVED (Phase 1 A/B): RNA-odds HELPS (+1.1–5.0%).** Keep it; derive `q_rna` (Amendment A1) rather than remove.
   Open sub-question: build a multi-isoform/capture-on stress locus to test where RNA coupling should decouple.
4. **PARTLY RESOLVED (Phase 1): the var~mean refit (90–246 ms) DOMINATES the sweep (25–33 ms) at R=157.** The perf
   PR must profile BOTH the refit and `_local_loglik` at genome scale, not just `_local_loglik`. Also bump
   `sweep_max_passes` (capon converged only on the last allowed pass).
5. **Does the chain smooth a true capture/CNV edge?** Phase-7 toy: a sharp captured exon flanked by depleted introns;
   assert the chain is a no-op on the IPR under uniform gDNA (factor-1) and `s_edge` widens where neighbours disagree.
6. **Robust-loss form (P7): Cauchy vs Huber** — Cauchy decouples harder (redescending), Huber is convex; choose by
   the capture-edge toy. Sub-question: `s_edge = sqrt(var_nbr)/μ` vs fitting `var(Δlog d) ~ mean log d` directly.
7. **Is the per-locus prior a live accuracy lever?** The user is certain (gDNA mass drives both the split *and* the
   eff-len IPR). The confidence blend ships behind a flag; flip on a measured ss0.50/AMBIG gain.

---

## 7. Doc status

- **Authoritative:** this plan; `effective_length_redesign_plan.md` §8; `per_node_deconv_hierarchy_design.md` §3
  (per-node math — update its RNA-odds-propagation framing to the inverted model).
- **Superseded by v4:** `CALIBRATION_PLAN_v3.md` (the RNA-chains framing is overturned), `CALIBRATION_PLAN_v2.md`,
  `CALIBRATION_PLAN.md`, `iterative_bootstrap_design.md`, `simplex_ambig_count_fallback_design.md`,
  `deconvolution_implementation.md`, `stage2_wiring_dryrun.md`.
- **Scratch, assimilated:** `calibration_nodes.md`, `my_notes.md`.
- **Keep:** `splice_junction_node_architecture.md` (deferred bipartite SJ — EM-side), `nascent_efflen_investigation.md`,
  `em_gdna_strand_likelihood_fix.md` (operative conclusion §9.4 only; an EM-side front, not calibration completion).
- Rewrite the `docs/README.md` calibration section to point at v4.
