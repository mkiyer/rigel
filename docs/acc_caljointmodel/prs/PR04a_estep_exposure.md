# PR 4a — Calibrator core: density seed + E-step (G2/G3) + exposure (G4)

**Parent plan:** [`../00_implementation_plan.md`](../00_implementation_plan.md) §4 D1/D7, §7 (PR 4).
**Spec:** [`../../caljointmodel/03_inference.md`](../../caljointmodel/03_inference.md) §3–§8, [`../../caljointmodel/01_generative_model.md`](../../caljointmodel/01_generative_model.md) §4/§7.
**Type:** Python-only (no index rebuild, no C++). **Build required:** no.
**Status:** **IN IMPLEMENTATION.** Decisions III.1–III.7 resolved by the user
(see [`../../TODO.md`](../../TODO.md)); this is the active sub-PR. The AMBIG
boundary-sweep is split out into **PR 4b** (drafted next, with pseudocode for
review before coding — III.3).

PR 4a turns the placeholder `calibrate()` into a real **single E-step pass**:
per region (and boundary-projected), split observed mass into gDNA (`M_g`) vs
RNA (`M_d`) — **G2** contained, **G3** boundary — and recover the per-region
gDNA exposure `ω` — **G4**, seeded by a principled global gDNA density `ρ_0`.
**Hyperparameters are fixed** at their initial values; fitting them and
iterating to convergence is **PR 5** (M-step + outer loop). PR 4a produces a
one-iteration result (`n_iterations = 1`, `converged = False`).

---

# Part I — Theory & design

## I.1 Global gDNA density `ρ_0` (III.2 — the Q2 resolution)

`ρ_0` is the library-wide gDNA fragment density (fragments per bp of accessible
gDNA). It is the **count channel's gDNA-only rate** and the **exposure
denominator** (`β_post = 1/φ + ρ_0·L_eff`). Earlier code seeded it by
**excluding the top 1 % density** (`top_t_fraction` — the clip-top-X% heuristic
the burn removed). PR 4a replaces that with a **signature-based Gamma posterior
with no percentile clip**.

**Seed set (III.2, user decision).** Seed `ρ_0` from the regions whose unspliced
mass is gDNA-dominated *by annotation* — **intergenic (NONE)** *and* **intronic**
regions, **together with the boundary mass crossing into them**. Exonic regions
are excluded.

> **Why not intergenic alone (the correction to the earlier draft).** Seeding
> from intergenic regions only **fails on hybrid-capture libraries**, where
> intergenic gDNA is depleted relative to the captured (genic) gDNA.
> Boundary-crossing fragments and intronic regions near captured exons carry the
> representative gDNA signal, so the seed must include **intergenic *and*
> intronic regions and their boundaries**.

The estimate is the `ρ_0` M-step (doc 03 §5.1) **restricted to the seed regions
with `ω ≡ 1`**, regularized by a unit-strength Gamma prior:

```
seed   = regions with no exon bit (intergenic ∪ intronic)
M_g_seed = Σ_seed (contained + left + right) unspliced mass     # D1 aggregation
ρ_0   = (α₀ + M_g_seed) / (β₀ + Σ_seed L_eff),    α₀ = β₀ = 1
```

**No top-fraction clip** (III.2): "the algorithm may work just fine without it;
we look at benchmarks and decide." With no non-exonic region anywhere (degenerate
— effectively never for a real genome) it falls back to the doc 03 §7
half-and-half density over all regions.

> **Known approximation (transparent, not a cliff).** Intronic unspliced mass
> also contains **nascent RNA**, so the seed is upper-biased where introns are
> actively transcribed. This is an explicit, benchmarked approximation; the
> planned remedy is the **expressed/unexpressed latent** (learns which genes are
> unexpressed → their introns/exons become clean gDNA seeds, even in AMBIG
> territory). That latent is **out of scope for PR 4a** (future PR).

## I.2 E-step (doc 03 §3) — identical on all three views

Per region, combine two log-Bayes-factors into a gDNA mixing proportion `π_g`,
then soft-allocate mass. Run the identical code on each substrate view
(`contained`, `left`, `right`).

- **Count (NB), doc 03 §3.1.** gDNA-only mean `μ_g = ω·ρ_0·L_eff` vs the
  mixture `μ_g + μ_d`: `LLR_count = logpmf_NB(n_u | μ_g) − logpmf_NB(n_u | μ_g+μ_d)`
  (`n = 1/φ`). **Silent when `μ_d = 0`** (the two hypotheses coincide) — which
  is exactly the single pass (no RNA mass estimated yet); it activates once the
  PR 5 outer loop feeds `M_d` back. See §I.5.
- **Strand (BB vs BB), doc 03 §3.2.** gDNA `BB(n_u, κ_d=0.5, ρ_d_bb)` vs RNA
  `BB(n_u, κ_rna, ρ_r_bb)` (PR 3's `κ_rna`/`ρ_r_bb`). The sense count `k_sense`
  is oriented by `ts_class`: POS → `n_unspliced_pos`, NEG → `n_unspliced_neg`,
  NONE → `pos` (arbitrary; gDNA is unstranded so neutral). The library
  direction is inside `κ_rna` (= `p_r1_sense`, measured with the same
  align-matches-transcript convention), so **no separate flip**.
- **D7 — AMBIG skips the strand term.** A region with transcripts on both
  strands has no valid sense split, so `LLR_strand = 0`; `π_g` comes from the
  count channel + prior. Intergenic (NONE) keeps the strand term (neutral).

Combine and allocate (doc 03 §3.3–3.6):

```
π_g       = expit(logit(π_g_prior) + LLR_count + LLR_strand)   # π_g_prior = 0.5
M_g_unspl = π_g · mass_unspliced                               # D1: π from flux, mass allocated
M_g       = M_g_unspl + ε_s · mass_spliced                     # spliced = deterministic RNA
M_d       = (mass_unspliced − M_g_unspl) + (1 − ε_s) · mass_spliced
```

**D1 (power from flux, density from mass).** The log-Bayes-factors are driven by
the integer **flux** (`n_unspliced`, `k_sense`); that `π_g` allocates the
fractional **mass**. For the contained view flux == mass; for the boundary views
the crossing flux drives `π_g` and the fractional per-side mass is split.

## I.3 Exposure posterior (G4, doc 03 §4) — D1 side-attribution (III.4)

Aggregate per-region gDNA mass with the **D1 side-attribution — no ½**:

```
M_g_tot = M_g_contained + M_g_left + M_g_right     # each boundary side already
                                                   # attributed to its one region
α_post  = 1/φ + M_g_tot
β_post  = 1/φ + ρ_0·L_eff
ω        = α_post / β_post          log_omega_var = 1/α_post
```

Closed form, `O(R)`. Each boundary side carries its own region's exposure
(III.4): the left/right sides of a boundary belong to different regions and fold
into those regions' `ω`. (Docs 01 §7 / 03 §4 still show the superseded `½`; a
Phase-8 doc reconciliation, tracked separately.)

## I.4 AMBIG regions in PR 4a (deviation — see §III)

AMBIG regions flow through the **same unified E-step** (strand term skipped via
`ts_class`); they are **not** special-cased with a separate fallback. Their gDNA
recovery is the count channel (density; live in PR 5) + the boundary sweep
(PR 4b). On the single pass — where the count channel is silent (§I.5) — AMBIG
(and every strand-uninformative region) gets `π_g = 0.5`. This is interim and
explicit; PR 4b/PR 5 complete the D7 recovery.

## I.5 What the single pass is (and is not)

The single pass uses `ω = 1`, `π_g_prior = 0.5`, and `μ_d = 0` (no prior RNA
mass — doc 03 §3.1). With `μ_d = 0` the **count log-Bayes-factor is identically
0**, so the pass is **strand-driven + spliced-deterministic**. The count channel
(and thus the intergenic→gDNA and paralog-rescue behaviour of doc 01 §10) only
engages once PR 5's outer loop feeds `M_d` back. **PR 4a builds and unit-tests
the full count channel** (with a non-zero `μ_d`), but the single-pass
*integration* invokes it dormant. This is faithful EM iteration 1, surfaced
explicitly so the one-iteration result is not mistaken for the converged
deconvolution.

---

# Part II — Implementation

```
src/rigel/calibration/
  density.py    # NEW: ρ_0 Gamma-posterior seed (intergenic+intronic+boundaries, no clip)
  estep.py      # NEW: count LLR + strand LLR + π_g + soft mass over one view (G2/G3)
  exposure.py   # NEW: closed-form Gamma exposure (G4), D1 aggregation (no ½)
  calibrate.py  # EDIT: seed ρ_0; one E-step pass over 3 views; aggregate exposure; iter-1 result
```

`density.py` / `estep.py` / `exposure.py` are internal (not re-exported from the
package surface — `__init__.__all__` is unchanged). `initial_hyperparameters`
drops `ρ_0` (density.py owns it) and returns `(φ, ρ_d_bb, ε_s)`.

**Constants (Q6, all in doc 03 §8 or unit-strength priors):**
- `density.py`: `_RHO_PRIOR_ALPHA = _RHO_PRIOR_BETA = 1.0` (unit-strength Gamma
  prior — one pseudo-fragment over one pseudo-base; floors `ρ_0 > 0`).
- `estep.py`: `_PI_CLIP = 1e-6` (doc 03 §8 — keeps `π_g ∈ (ε, 1−ε)`; machine
  guard, not a tunable). NB/BB use scipy; `φ`, `κ_rna`, `ρ_r_bb`, `ρ_d_bb`,
  `ε_s` are the already-set hyperparameters.

No clip-top-X%, no `confidence_floor`, no outlier cliff.

## Tests

- `test_density.py`: seed mask = non-exonic (intergenic+intronic, excludes exon);
  Gamma posterior math incl. boundary mass; fallback when no seed; `ρ_0 > 0`.
- `test_estep.py`: strand favours RNA → `π_g→0`; strand neutral (`κ_rna=0.5`) →
  `π_g≈0.5`; AMBIG skips strand; **count channel with non-zero `μ_d`**
  (gDNA-consistent count → `LLR>0`; excess count → `LLR<0`); mass conservation
  (`M_g+M_d == mass_unspliced+mass_spliced`); flux-drives-π / mass-allocated.
- `test_exposure.py`: closed form vs hand calc; D1 aggregation (no ½);
  `log_omega_var = 1/α_post`.
- `test_calibrate_single_pass.py` (rewrite of the placeholder test): invariants —
  per-view conservation, `0 ≤ M_g ≤ mass`, `ω > 0`, `n_iterations=1`,
  `converged=False`, `ρ_0` from the density seed, `κ_rna` from the StrandModel.

## Acceptance gate

- New tests pass; PR 1–3 suites green; `ruff` clean.
- `calibrate()` returns a schema-valid result; full-suite failure mode stays the
  post-calibrate `NotImplementedError` (PR 6).

---

# Part III — Decisions (resolved) & deviations from the original draft

**Resolved (user, TODO.md):** III.1 split 4a/4b ✓ · III.2 seed = intergenic +
intronic + boundaries, **no clip** ✓ · III.4 exposure D1 no-½, per-side
exposure ✓ · III.6/III.7 single pass, `π_g_prior=0.5`, `n_iterations=1`,
`converged=False` ✓ · III.5 magic numbers — satisfied ✓.

**Deviations from the original combined PR04 draft (noted per the user's "note
any deviations"):**

1. **No AMBIG `min(n_u, ρ_0·L_eff)` fallback in 4a.** The draft proposed a
   special global-density fallback for AMBIG regions. Removed: AMBIG flows
   through the unified E-step (strand skipped), and its density recovery comes
   from the count channel (PR 5) + sweep (PR 4b). This **removes a heuristic
   `min()` cliff** (aligns with Q6) and is more faithful to doc 03. Interim
   single-pass effect: AMBIG gets `π_g=0.5` (§I.4).
2. **`ρ_0` seed includes boundary mass and intronic regions** (III.2), not just
   intergenic contained mass — and uses the D1-aggregated mass (no clip),
   matching the doc 03 §5.1 `ρ_0` M-step restricted to seed regions.
3. **Count channel is dormant on the single pass** (`μ_d=0` ⇒ `LLR_count=0`).
   Not a code change — a faithful property of EM iteration 1 (doc 03 §3.1),
   surfaced explicitly (§I.5). The channel is fully built and unit-tested.

## PR 4b (next, drafted for review before coding — III.3)

AMBIG boundary-sweep imputation: propagate gDNA-density evidence from
strand-decodable neighbours into AMBIG territory. The recovered
`boundary_sweep`/`boundary_model` mechanism, **scrubbed of `confidence_floor`**
and replaced by a smooth, unit-pseudocount, evidence-proportional transfer
weight. The algorithm + pseudocode will be drafted in `PR04b_ambig_sweep.md`
for the user's review **before** any 4b code (III.3).

## Rollback

Revert the sub-PR. `calibrate()` reverts to the PR 3 placeholder (no-gDNA /
unit-exposure + the real strand model). No on-disk artifacts change.
