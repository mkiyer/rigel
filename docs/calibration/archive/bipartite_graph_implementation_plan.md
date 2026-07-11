# Bipartite region↔boundary↔region graph — implementation plan

**Status:** the authoritative plan for promoting boundaries to **first-class solved nodes** (the R↔B↔R linear
seam graph). Primary authority: `CALIBRATION_ARCHITECTURE.md` (the §0 count-zero-information invariant, the
precision model, signature-binary init, the solve flow). Mechanism: `CALIBRATION_PLAN_v6.md` (the 3 firm
decisions), `splice_junction_node_architecture.md` (linear seam vs the deferred transcript-SJ graph),
`rna_imputation_transcript_structure.md` (mature anchored / nascent imputed), `imputation_variance_model.md`
(the Poisson-offset var~mean, the density→fraction jacobian, identity mean).

This plan was produced by a map → draft → adversarial-critique pass; the BLOCKER/MAJOR findings are folded in
below as **resolved decisions** (§4) with the few genuinely-open items in §7.

---

## 0. What already exists, and the scope guard

**The data layer is DONE — no C++.** `substrate.BoundarySubstrate.from_payload` (substrate.py:174-189, built,
tested, currently *unconsumed*) gives each boundary node `b ∈ [0, b_obj_total)`: its two sides
(`left`/`right` `SubstrateView`s = `boundary_*_left[b]`/`boundary_*_right[b]`), its flanking regions
(`left_region`/`right_region`, **−1** at a reference terminal), and the one-sided motif-stranded spliced.
The chain adjacency is `region_arrays.region_boundary_indices` / `boundary_region_indices` (exact inverses).
The re-indexing identity (tested): `BoundarySubstrate.left[b] == CalibrationSubstrate.right[left_region[b]]`,
`BoundarySubstrate.right[b] == CalibrationSubstrate.left[right_region[b]]`.

**The solver is already node-agnostic.** `simplex._mixture_strand_loglik`, `simplex_sweep._local_loglik`, and the
reductions `_fg_median`/`_fg_var`/`_axis_mean` read only per-node arrays + the lattice (no `substrate`/
`region_arrays` access). `NodeDeconv` already carries `f_pos/f_neg/f_g` + masses. So a boundary can be fed
through the *same* solve once three pieces exist (the boundary node-class helper, the boundary RNA side eff-len,
and the region→boundary var~mean direction).

**SCOPE GUARD (firewall the deferred C++ graph).** This builds ONLY the linear genomic-seam graph: a boundary
node is keyed **solely by its seam index `b`** (`boundary_region_indices`). It MUST NOT be keyed by intron
coordinates, MUST NOT add transcript-adjacency edges across skipped introns, and MUST NOT touch the accumulator
deposit signature — those are the **deferred** transcript-SJ graph (`splice_junction_node_architecture.md`, needs
C++). `splice_junction.splice_junction_eligibility(sig_l, sig_r)` (reads only the two flank signatures) is the
sanctioned motif-strand classifier.

---

## 1. The boundary node

A boundary node `B` is **one belief unit** (Decision Q5: its two sides are mass-accounting, not two beliefs).

- **Pie** `{f_pos, f_neg, f_gdna}` over its **crossing UNSPLICED** mass, the same 2-simplex latent as a region.
  One `f_g_b` applies to *both* sides' masses.
- **Mass (the pie split base — CONSERVED, summed):** `mass_unspl_b = left.mass_unspliced[b] +
  right.mass_unspliced[b]` (the fair partitioned mass; §4.1).
- **Trial count `N` for the strand BB — ONE SIDE, not the sum (§4.1, verified):** a contiguous unspliced
  crossing credits *both* sides for the *same* fragment, but `flux_left[unspliced] == flux_right[unspliced]`
  exactly (measured 1209/1209). So `N_b = max(flux_left, flux_right)[0:2]` = the distinct crossing-fragment count
  (one side; `max` guards terminals); per-strand `u_pos_b = max(fl,fr)[0]`, `u_neg_b = max(fl,fr)[1]`. Count and
  mass are DIFFERENT quantities. `_mixture_strand_loglik` is used **verbatim** (count enters as its honest power;
  no closed-form `N(2κ−1)²` substitution — that would drop the overdispersion that makes deference emergent).
- **Spliced mature anchor (ONE-sided, ONE strand, a FIXED floor — §4.16):** a boundary is exactly one splice
  junction = one genomic position = one strand, so its spliced is on exactly ONE side (the exon flank for that
  junction's motif strand; the other side is 0). Take the spliced from the **appropriate side per strand** via
  `splice_junction_eligibility(sig_l, sig_r)` (which returns the exon side + the motif strand) — **never**
  `left + right` (summing conflates the two distinct junctions of a rare opposite-strand seam). It enters
  `_local_loglik` through the EXISTING sided-spliced lower-bound term (a hard floor `f_motif·N_b ≥ spliced` on
  that motif strand), which is ~dead on regions today and becomes LIVE on boundaries. **§0 status (§4.6):** the
  spliced is a *directly observed* mature component, not a count vote — a fourth legitimate input.
- **Eff-lengths (two, both exist):** (a) the pooled crossing DENSITY support
  `S_b = ½(boundary_side_eff_length(L_left) + boundary_side_eff_length(L_right))` (`= E_g[min(ℓ,L_side)]`
  averaged; the SAME `0.5·(gdna_boundary_len[r]+gdna_boundary_len[r+1])` priors.py:167 already uses — must be
  **one shared function**, §4.7) — used for the f_g density + the global/imputation jacobian; (b) the strand BB
  uses COUNTS, no eff-len. A `−1` flank contributes 0 mass and 0 support (the present-flank length stands alone).
- **Strand class / allow mask (the one new derivation, §4.5):** a boundary has no own `strand_class`; derive
  `strand_obs_b` from `strand_deconv._side_strand_orientation` (the `{POS,POS}`/`{POS,NONE}`-wildcard consistency
  rule, a keeper). The **gDNA** `allow_pos_b`/`allow_neg_b` = flank-OR (`has_s(left) | has_s(right)`); a `{POS,NEG}`
  seam is strand-UNOBSERVABLE yet admits both strands (AMBIG-like). The **RNA imputation** allow mask is STRICTER
  (the TSS/TES black hole, §4.5): RNA crosses a boundary only where the strand-s transcript bit is on **both**
  flanks.

---

## 2. The chain + the region↔boundary conversion

**Chain (existing index arithmetic, no C++):** region `r`'s RIGHT boundary == region `r+1`'s LEFT boundary == a
seam slot `b`. `R_total + (#nonempty refs) == B_obj_total`. `−1` flanks are masked exactly like the region chain's
`same_ref_left_right` (§4.4 — never let Python index `region[-1]`).

**Conversion (the message is the PIE + geometry; NO mass moved — compare in DENSITY space, each node by its own
eff-len):** a source shares `ρ̂ = f_src · M_src / eff_src`; the dest forms the prior mean
`μ = ρ̂ · eff_dst / M_dst` (**identity** in density, slope 1 by the factor-1-under-uniform eff-len — a slope ≠ 1
is an eff-len bug, `imputation_variance_model.md` §3). gDNA region density uses `region_eff_length` (`E_g[max(0,
L−ℓ)]`); boundary crossing density uses the averaged `S_b`. The precision is `1/(σ²_bio(μ) + ρ_src/L_src)`
(Poisson-offset var~mean + the predictor sampling floor) via the `(eff/mass)²` jacobian.

**RNA needs a NEW eff-len (§4.5):** add `boundary_side_eff_length(rna_fl_pmf, region_size_bp)` (mirrors the gDNA
one at calibrate.py:94). The boundary **unspliced-RNA** pie density uses the **averaged two-flank** RNA side
length; the **spliced mature** density uses the **single exon-flank** RNA side length (NOT averaged).

---

## 3. The unified solve

Factor the per-node-array assembly + reduction out of `deconv_regions_sweep` into a shared
`_solve_nodes(u_pos, u_neg, n_trial, spliced_pos, spliced_neg, allow_pos, allow_neg, strand_obs, mass_unspl,
mass_spliced, kappa, od_g, od_r, n_grid, mu_global, global_tau, gdna_imp_*, rna_imp_*) → NodeDeconv` core (it just
calls `_local_loglik` + the reductions — already node-agnostic). Two thin wrappers build the arrays:

- `deconv_regions_sweep` (region): arrays from `substrate.contained` + `region_arrays.strand_class` (today's
  derivations). **Regions keep their own spliced floor LIVE** (reads `contained.n_spliced_sense`, §4.2).
- `deconv_boundaries_sweep` (boundary, NEW): arrays from `BoundarySubstrate` — pooled `u_pos/u_neg/mass`, the
  de-duplicated trial count `N_b` (§4.1), the one-sided motif-stranded spliced floor (`splice_junction_
  eligibility` selects the strand), `allow_*`/`strand_obs` from the boundary node-class helper (§4.5), the
  boundary global baseline `μ_global_b = clip(ρ_global · S_b / M_pooled_b)` (§4.8). Replaces `deconv_sides`.

The 3-tier ψ (strand mixture + node-class prior [Jeffreys / global] + imputation pulls) + the spliced floor apply
**uniformly** to region and boundary nodes. **No I0**, no `g_count` composition vote — deference is emergent (the
BB likelihood's intrinsic curvature dominates a diffuse prior when sharp, is flat at κ=½/AMBIG).

---

## 4. Resolved design decisions (the dry-run BLOCKER/MAJOR fixes)

**4.1 (BLOCKER) Count vs mass at a boundary are DIFFERENT quantities (verified in the accumulator).** A
contiguous unspliced crossing credits BOTH sides' flux for the SAME fragment (accumulator.cpp:203,208), so
`flux_left + flux_right` double-counts. MEASURED on the sim payload: `flux_left[unspliced] ==
flux_right[unspliced]` **exactly** at every internal boundary (1209/1209) — each contiguous crossing increments
each side once. So the three boundary quantities are distinct and NOT interchangeable:
- **trial count `N` (BB statistical power) = ONE side** `= max(flux_left, flux_right)[0:2]` — the distinct
  crossing-fragment count (both sides provably equal for unspliced; `max` guards terminals where one side is the
  off-edge zero). NOT the sum (doubles the Fisher info), NOT `½(sum)` (a band-aid). One observation per fragment
  per boundary — §0-clean. The per-strand tilt counts `u_pos = max(fl,fr)[0]`, `u_neg = max(fl,fr)[1]`.
- **pie mass (the gdna/rna split base) = `mass_left + mass_right`** — the CONSERVED, partitioned, *fair* value
  (a fragment's deposits across the whole graph sum to 1; mass < one-side-flux globally because multi-boundary
  "encompassing" crossings split mass but are a full observation at each boundary).
- **spliced floor (per strand) = the spliced on the exon-flank SIDE for that junction's motif strand** (via
  `splice_junction_eligibility`, §4.16) — NOT `flux_left + flux_right`. A boundary is ONE junction = one strand,
  so its spliced sits on exactly one side; the rare both-sided seam (137/596) is two DISTINCT junctions on
  opposite strands, each feeding its OWN strand's floor from its OWN side. Summing would conflate the strands.
Keep `_mixture_strand_loglik` verbatim (count enters as its honest power; no closed-form substitution). *Gate:* a
pure-gDNA boundary with N crossing fragments has the SAME `f_g` posterior variance as a region with N contained
fragments (no √2 sharpening); mass conservation Σ(mass over the graph) = the fragment count.

**4.2 (RESOLVED — user ruling) Contained spliced is IGNORED; there is NO region spliced floor.** Decision: we do
NOT model unannotated within-region transcripts. Contained spliced is exactly 0 in simulated data (MEASURED:
0/1221 on the 5Mb suite, 0/53 on a fresh small sim — every sim junction is annotated, hence a region boundary,
so a spliced fragment always crosses a boundary). So spliced is **100%-boundary-owned by construction**, the
region pie is **pure unspliced** (no spliced floor term on regions), and `mass_rna_spliced` is built from the
**boundary** spliced only (§4.3). The critic's "keep the region floor for real-data robustness" is overruled by
the no-unannotated-modeling decision.

**4.3 (BLOCKER) `mass_rna_spliced` + the prior withhold from the BOUNDARY spliced.** Boundary `b`'s one-sided
spliced projects onto its **exon-flank region's** `mass_rna_spliced` (re-indexing identity); the unspliced
crossing `(1−f_g)·M` projects onto `mass_rna_left/right`. The `assemble_priors` subtraction
(`rna_region = contained+left+right − spliced`, priors.py:250-256) is **unchanged**. Contained spliced is ignored
(§4.2 — 0 in sim). *Gate:* `Σ mass_rna_spliced == Σ boundary spliced`; `rna_prior_count` never negative-clamps.

**4.4 (BLOCKER) `boundary_side_seeds` re-homed off the deleted helpers, kept PRE-solve.** The gDNA-strand-
overdispersion fit runs before the solve and must stay non-circular (its seed weight is the strand mean ½). Before
deleting `_compute_side`/`_SideQuantities`/`_left_right_neighbors`, re-express the seed extraction over
`BoundarySubstrate` using the keepers (`_side_strand_orientation` + `density_model.count_observable_masks` + raw
crossing counts). *Gate:* `gdna_strand_overdispersion` unchanged after the deletion (regression test).

**4.5 (MAJOR) RNA black-hole allow mask, separate from gDNA strand-observability.** The gDNA strand-observability
uses `_side_strand_orientation` (NONE wildcard — keep). The **RNA** imputation across a boundary requires the
strand-s bit on **both** flanks; an intergenic flank BLOCKS RNA (rule (a) TSS/TES black hole). Do not let the
gDNA flank-OR admit RNA into an intergenic region.

**4.6 (MAJOR) §1.4 annotation — the spliced floor is a directly-observed mature component.** Add to
`CALIBRATION_ARCHITECTURE.md` §1: the spliced is a guaranteed-RNA *observation* (the fragment physically skipped
an intron; no gDNA hypothesis), so it is a hard lower bound whose confidence scales with the observed spliced
count — observation, not a composition vote, distinct from the two count→precision channels. The SOLVE floor is
count-space (`spliced/N_b`); the DENSITY conversion uses the one exon-flank RNA length (§2).

**4.7 (MAJOR) ONE shared pooled-seam support function.** The boundary solve's `S_b` and
`priors._gdna_region_node_arrays`' seam support must be the *identical* formula
(`0.5·(gdna_boundary_len[r]+gdna_boundary_len[r+1])`). Factor it into one helper called by both. *Gate:* the
boundary node's support == priors' seam support for every internal boundary; factor-1-under-uniform holds.

**4.8 (MAJOR) Boundary global baseline geometry + FL-consistency.** `μ_global_b = clip(ρ_global · S_b / M_pooled_b,
0, 1)`. Confirm `ρ_global` (fit over contained `E_g[max(0,L−ℓ)]`) is density-consistent applied to a crossing node
(support `E_g[min(ℓ,L)]`); if not, apply the eff-len-ratio correction (`splice_junction.density_frac_to_count_frac`
idiom). *Gate:* factor-1-under-uniform for the boundary global baseline.

**4.9 (MAJOR) The region→boundary var~mean direction is BUILT, not toggled.** `pair_imputation_points` emits only
boundary-side→region today (variance_model.py:410-411 defers the reverse). Add the boundary-as-dest point: `mean
= boundary_density[b]`, `raw_var = (boundary_density − region_density[flank])²` per present flank, `offset =
ρ_b/S_b + ρ_region/L_region`. Pool both directions into ONE curve (more data). **Until built, the boundary
imputation prior is a no-op (τ=0)** → the boundary falls to strand + global (a correct degenerate, not a corrupt
prior).

**4.10 (RESOLVED — user deferred to me) Minimal-churn region-keyed projection (schema unchanged).** Chosen for
the smaller blast radius + zero fixture/`__post_init__`/package-surface churn, and because the pooled-seam node
IS the boundary node so the projection is exact (not a lossy compromise); first-class per-boundary `float64[B]`
fields remain a clean future option if perf/clarity ever warrants. Keep the per-region `CalibrationResult` schema;
project the boundary solve back via the re-indexing identity — written EXPLICITLY (easy to invert):
`mass_gdna_right[left_region[b]] = f_g_b · left_side.mass_unspliced`;
`mass_gdna_left[right_region[b]] = f_g_b · right_side.mass_unspliced` (RNA analogs, spliced onto the exon-flank
region). Guard `−1`. This keeps `result.py`/`_gdna_region_node_arrays`/`capture_eff_length`/all fixtures
structurally intact (only values shift). `_gdna_region_node_arrays` MUST keep returning **region-keyed** arrays
(`capture_eff_length` addresses them by `_exon_region_incidence` — re-keying by boundary silently breaks the
transcript contraction). *Gate:* `test_calibrate` per-side conservation totals exact; `test_priors`/
`test_capture_eff_length` factor-1 green.

**4.11 (RESOLVED — user ruling) `gdna_deconv_quantile` — DELETE.** `_fg_median` stays at q=0.5 (the current
default — a no-op, so no golden churn from the quantile's removal). Delete `config.gdna_deconv_quantile` + its
`>0` validator + the quantile tests (`test_calibrate.py:94-112`, the `test_strand_deconv` FP-quantile tests)
atomically with the `deconv_sides` retirement (Phase 4). No re-plumbing into the simplex solve.

**4.12 (MAJOR) Complete the deletion inventory.** Retire: `deconv_sides`, `_deconv_per_node`, `_compute_side`,
`_SideQuantities`, `_left_right_neighbors`, `cleaned_gdna_count`, and `strand_deconvolve` + `StrandSplit` (their
only role is feeding `cleaned_gdna_count`); `config.gdna_strand_info_scale` (I0) + its validator. **Trace before
deleting** `_grid_posterior_quantile` / `strand_posterior_gdna_frac` — keep whichever the surviving solve/quantile
needs. Rewrite/delete the orphaned tests (`test_strand_deconvolve.py`, the `_deconv_per_node`/FP-quantile tests,
`test_calibrate.py` quantile + hardcoded mass-value tests).

**4.13 `derive.gdna_density_global`** sums over the **B boundary nodes** once (pooled gdna over `S_b`) + the region
contained nodes — avoiding the double-count of a crossing's two halves. QC-only; confirm nothing asserts its exact
value (cli/pipeline/sim only log it).

**4.14 Package surface unchanged.** Keep `BoundarySubstrate`/`deconv_boundaries_sweep`/the boundary helpers
INTERNAL (not re-exported); `test_package_surface.py` set-equality stays green. The stable external interface is
the 3 per-locus prior scalars → C++ EM, which does not change.

**4.15 Geometry deferred (no premature lock-in).** Ship the per-component density-space imputation (identity mean,
`(eff/mass)²` jacobian, Poisson-offset var~mean). Record the per-component-vs-log-ratio precision question as a
POST-Phase-C measurement gated on the complex-loci battery — only revisit if `f_g`/`f±` are shown to fight under
the simplex. Do NOT introduce a log-ratio reparameterization in this build.

**4.16 (RESOLVED — user ruling) The boundary spliced is ONE-sided / ONE-strand; use the appropriate side per
strand, never `left + right`.** A boundary is exactly one splice junction: a single genomic position is one
junction on one strand (+ or −, not both). So the boundary's spliced mass sits on exactly ONE side (the exon
flank for that junction's motif strand; the other side is structurally 0). The boundary's spliced floor is built
**per strand** from the spliced on the side `splice_junction_eligibility(sig_l, sig_r)` identifies for that motif
strand — `left.n_spliced[exon_side]` **or** `right.n_spliced[exon_side]`, **not** their sum. Summing is invalid:
for the rare opposite-strand seam (measured 137/596 boundaries with spliced on both sides — two DISTINCT junctions
sharing a genomic position), `left + right` would mix the `+` junction's mature with the `−` junction's; instead
the `+` floor reads its side and the `−` floor reads the other. For a single-junction boundary only one strand's
floor is set. (This is consistent with the gДНК `N_b` being one-sided too, §4.1 — both reflect "a boundary credits
one side for one event/junction"; the unspliced both-sided equality is a contiguous-crossing artifact, the spliced
one-sidedness is the junction's physical strand.) *Gate:* extend `test_spliced_boundary_onesidedness` to assert
the per-strand floor reads only its own exon-side, and an opposite-strand seam feeds the two strand floors from
the two sides separately (never summed).

---

## 5. The phases (resequenced — each independently testable, with a gate)

> **Resequencing (dry-run finding):** do NOT commit the uncommitted region→region RNA imputation — it is
> architecturally superseded (it skips the boundary that owns the nascent, so it cannot fix the intron
> false-gDNA, and committing+goldening it double-churns). Set it aside as scratch. The first committed checkpoint
> is the pure refactor.

**Phase 1 — extract `_solve_nodes` (pure refactor, bit-identical).** Factor the node-agnostic core out of
`deconv_regions_sweep`; the region wrapper builds its arrays and calls it. NO behavior change.
*Gate:* `test_simplex_sweep` (the 3 per-node behaviours) green UNCHANGED; full calibration suite bit-identical;
golden unchanged. *Risk: low* (the cut line is clean — `_local_loglik` is already node-agnostic). **First commit.**

**Phase 2 — the boundary node-class helper + `deconv_boundaries_sweep`, solved in ISOLATION.** Build
`boundary_node_class(left_region_strand, right_region_strand) → (allow_pos, allow_neg, strand_obs)` (reuse
`_side_strand_orientation`); `deconv_boundaries_sweep` over `BoundarySubstrate` with the de-duplicated `N_b`
(§4.1), the one-sided motif-stranded spliced floor, the averaged `S_b`, `μ_global_b` (§4.8). Add
`boundary_side_eff_length(rna_fl_pmf, …)`. The boundary imputation prior is a **no-op** here (§4.9 not yet built).
NOT yet consumed by `calibrate`'s loop or `priors`.
*Gate:* NEW boundary unit tests (the `test_simplex_sweep` 3-behaviour template on a boundary: intergenic-flank →
gDNA vertex; `{POS,POS}` resolves; `{POS,NEG}` admits both + global breaks degeneracy; the one-sided spliced floor
`f_motif·N_b ≥ spliced`); the `N_b` de-dup test (§4.1); the re-indexing identity (`test_boundary_substrate`) green.
*Risk: medium* (the flank-OR allow mask + the motif-strand floor placement are the new derivations).

**Phase 3 — wire boundaries into the per-pass chain + route imputation region→boundary→region.** Replace
`deconv_sides` (once, pre-loop) with `deconv_boundaries_sweep` INSIDE the per-pass loop (co-converge). Build the
region→boundary var~mean direction (§4.9). Rewire `gdna_imputation_prior` predictors from the frozen deconv_sides
sides to the SOLVED boundary nodes; **REPLACE** (not layer) the region→region RNA adjacency with
region←boundary←region (the RNA black-hole mask §4.5). Per pass: (a) re-fit `ρ_global` + `τ_global` + the
bidirectional var~mean; (b) solve regions from flanking boundary nodes; (c) solve boundaries from flanking region
nodes; (d) carry both; convergence on max(region, boundary) `f_g` deltas.
*Gate:* MONOTONE convergence on a capture-on+nascent toy (under-relaxation η only if it oscillates, documented);
`rna_imputation_diagnostic` intron false-gDNA DROPS (the structural fix — nascent now crosses exon↔intron via the
boundary; the mature is absorbed at the boundary); `test_ambig_scenario` green; net-flow non-regressing, mature-FP
not rising. *Risk: high* (the central coupling; doubled node count — re-verify `sweep_max_passes`).

**Phase 4 — retire `deconv_sides`/I0; rewire `priors`/`derive` to boundary nodes (minimal-churn).** Delete the
§4.12 set + I0 + validator. Re-home `boundary_side_seeds` (§4.4). Project boundary masses back region-keyed (§4.10,
explicit crossover). `derive.gdna_density_global` over B nodes (§4.13). Re-plumb the quantile (§4.11). Keep
`mass_rna_spliced` + the withhold correct (§4.3).
*Gate:* `test_priors` factor-1 + mass conservation green (preserved by the projection); `test_capture_eff_length`
factor-1 green; the seed-fit regression test (§4.4); rewrite/delete the orphaned strand_deconv tests; golden
`--update-golden` ONLY AFTER the calibration-benchmark + complex_loci gate confirms the shift is an improvement.
*Risk: medium* (deletions strand tests — expected, rewrite; the I0 regression is pre-accepted).

**Phase 5 — validation, tolerances, cleanup, docs.** Run the full validation surface; tighten the tolerances this
build is predicted to improve (the weak-SS phantom; nrna corners). Refresh stale docstrings (`sweep_max_passes`;
`simplex_sweep` "Boundary sides keep deconv_sides"; the `calibrate` preamble). Mark v6 Phase C / ARCHITECTURE
Step 3 done.
*Gate:* `complex_loci_benchmark` TOTAL ≤ the pre-bipartite baseline (record it in Phase 1) except a traced/recovered
transient; net-flow non-regressing; full suite green; the schema mass-conservation invariant (Σ node mass = total).

---

## 6. Invariants that must hold at every phase

- **Zero-gDNA pin:** `rna_imputation_diagnostic` introns AND exons → `f_g ≈ 0` (the strand alone).
- **Mass conservation:** per-node `gdna + rna = total`; the projection preserves per-side totals exactly; a
  spliced fragment is counted in exactly one view.
- **Factor-1-under-uniform:** the boundary `S_b` == priors' seam support (one shared function); the density
  conversions slope-1.
- **§0:** count enters only as statistical power — the strand BB (`N_b`, de-duplicated) + the imputation Poisson
  floor; the spliced is a directly-observed mature component (§4.6); NO count vote, NO I0, NO magic constant.
- **The stable external interface** (the 3 per-locus prior scalars → C++ EM) does not change; the build is internal
  to the calibration package.

---

## 7. Decisions (was open questions — all resolved)

1. **`N_b` count form** — RESOLVED (§4.1, verified): use **one side** (`max(flux_left, flux_right)[0:2]`), the
   distinct crossing count, since `flux_left == flux_right` exactly for unspliced. NOT the sum, NOT `½(sum)`.
2. **Schema layout** — RESOLVED (user deferred to me, §4.10): **minimal-churn region-keyed projection** (smaller
   blast radius; the pooled-seam node IS the boundary node so the projection is exact).
3. **`q_rna` replacement** — RESOLVED (user): **(i) iterated one-hop priors** (the current shape). A/B the
   within-pass edge later only if the complex battery shows the one-hop chain under-propagates.
4. **`gdna_deconv_quantile`** — RESOLVED (user): **DELETE** (§4.11). `_fg_median` stays q=0.5; drop the config
   field + validator + the quantile tests atomically with the `deconv_sides` retirement.
5. **Region spliced floor** — RESOLVED (user): **drop it** (§4.2). Contained spliced is ignored (0 in sim; we do
   not model unannotated within-region transcripts). Spliced is 100%-boundary-owned; regions are pure unspliced.
6. **Boundary spliced side** — RESOLVED (user, §4.16): a boundary is exactly ONE splice junction (one genomic
   position = one strand), so its spliced lives on exactly ONE side. Use the spliced on the **appropriate side**
   (the exon flank for that junction's motif strand, via `splice_junction_eligibility`) **per strand** — do NOT
   sum left+right (summing would conflate the two distinct junctions of a rare opposite-strand seam).

---

## 8. New code (deletions outnumber additions — net health positive)

- `simplex_sweep._solve_nodes` (extracted core) + `deconv_boundaries_sweep` (wrapper).
- `boundary_node_class(...)` helper (small; next to `_side_strand_orientation`, or a tiny `boundary_nodes.py` if
  the array-assembly grows).
- One shared pooled-seam support function (§4.7), called by the boundary solve + `priors`.
- `pair_imputation_points` gains the region→boundary point emission (§4.9, in-place).
- One new eff-len call: `boundary_side_eff_length(rna_fl_pmf, region_size_bp)`.
- **Deleted:** `deconv_sides`/`_deconv_per_node`/`_compute_side`/`_SideQuantities`/`_left_right_neighbors`/
  `cleaned_gdna_count`/`strand_deconvolve`/`StrandSplit` + I0. One node model, one solve core, one ψ, no I0
  blend, no separate count-cleaning pre-pass.
