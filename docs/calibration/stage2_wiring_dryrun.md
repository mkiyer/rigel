# Stage 2 wiring — implementation dry-run ((2) builders + (3) wire the real precisions)

A walk-through of the code path before writing it, to pin the **units** and **data flow** (where the bugs
hide). Scope: a **single pass** = "Pass 0" of the iterative loop (`CALIBRATION_PLAN_v2.md` §2). The iteration
control flow (§8.4) wraps this unchanged. Confidence guard is **removed** (all means seen; the var~mean
learns high variance at RNA-rich means directly).

## The objects in play

- Per node we have, in `calibrate`: `region_eff_len[i]` (gDNA region eff-len `L_i`), the boundary eff-len
  `fl_mean`, `substrate.{contained,left,right}` unspliced counts, `node_density` (the count module:
  `count_gdna_frac`, `density`, `global_density`=ρ_global, `region_count_observable`).
- The var~mean fitter `MonotoneVarMean` models **variance of a gDNA *density* estimate vs the density
  mean** (it is fit on density-space points). Its `predict(μ)` returns `var_density(μ)`.

## Units — the one thing to get right

The count prior in the sweep is a Gaussian on the **fraction** `f_g`:
`−½·τ_count·(f_g − count_gdna_frac)²`. The var~mean gives a **density** variance. Convert once, per node:

```
f_g = ρ · L / mass                          (density → fraction; ρ = gDNA density, mass = unspliced mass)
⇒  Var(f_g) = Var(ρ) · (L / mass)²          (linear map)
⇒  σ²_frac(i) = var_density(μ_i) · (L_i / mass_i)²
⇒  τ_count(i) = 1 / σ²_frac(i) = mass_i² / ( var_density(μ_i) · L_i² )
```

where `μ_i` = the node's gDNA **density** at which we read the var~mean. **In Pass 0 (all-gDNA init)
`μ_i = total_unspliced_i / L_i`** (the `varmean_points` mean). The same identity gives the **global prior**
variance: evaluate the DIRECT curve at ρ_global and convert with the node's own geometry:

```
σ²_global(i) = var_density_DIRECT(ρ_global) · (L_i / mass_i)²
τ_global(i)  = 1 / σ²_global(i)
μ_global(i)  = clip(ρ_global · L_i / mass_i, 0, 1)            (unchanged from today)
```

So `τ_global` becomes **per-node** (was the flat scalar 1.0), tight where ρ_global is well-determined
(low-variance regime ⇒ zero-DNA pin) and wide where not — exactly the intended foundation.

## (2) The builders — which curve for which node

`varmean_points(substrate, region_arrays, region_eff_len, fl_mean)` already yields per-node
`(mean=μ_i, raw_var, region_observable)` from the boundary↔region↔boundary triplets (the all-gDNA total at
Pass 0). Split:
- **DIRECT** = `fit_direct_varmean(points)` — over `region_observable` nodes (the region anchored its own
  estimate). Used for: count precision at count-observable nodes, and `σ²_global`.
- **IMPUTATION** = `fit_imputation_varmean(points)` — over the imputed nodes (boundaries only). Used for:
  count precision at imputed/AMBIG nodes (honestly high at enriched exons ⇒ count yields there).

Per-node `var_density(μ_i)` = `DIRECT.predict(μ_i)` if `region_count_observable[i]` else
`IMPUTATION.predict(μ_i)`.

## (3) Wiring — what changes in `calibrate` (sweep branch) and `simplex_sweep`

`calibrate` (the `use_propagation` branch), replacing the β placeholder and the flat global τ:

```python
pts   = varmean_points(substrate, region_arrays, region_eff_len, fl_mean)
direct = fit_direct_varmean(pts)
imp    = fit_imputation_varmean(pts)
mu_i   = total_unspliced / region_eff_len                       # Pass-0 density (all-gDNA)
var_d  = np.where(node_density.region_count_observable, direct.predict(mu_i), imp.predict(mu_i))
geom2  = (region_eff_len / np.maximum(mass_unspl, eps))**2
tau_count  = mass_unspl**2 / np.maximum(var_d * region_eff_len**2, eps)     # = 1/(var_d·geom2)
sig2_glob  = direct.predict(np.full_like(mu_i, rho_global)) * geom2
tau_global = 1.0 / np.maximum(sig2_glob, eps)
# clip both precisions to a sane band so a degenerate var (→0 or →∞) can't dominate/vanish
tau_count  = np.clip(tau_count,  TAU_LO, TAU_HI)
tau_global = np.clip(tau_global, TAU_LO, TAU_HI)
```

Pass `count_precision=tau_count` (already a param) and a new **per-node** `global_tau=tau_global` to
`deconv_regions_sweep`; in `_local_loglik` the global term becomes
`−½·global_tau[:,None]·(f_g−μ_global)²` (currently a scalar 1.0). Drop the `count_trust_beta` use on this
path (keep the field deprecated). The `(1−w)` count down-weighting **stays** (it makes the count yield to
the strand+propagation; orthogonal to the precision magnitude).

## Edge cases & guards (the dry-run's payoff)

1. **Degenerate var → ∞ precision.** `var_d→0` (a near-constant fit region) would make `τ_count→∞` and pin
   `f_g` to a possibly-wrong count. **Clip `τ` to `[TAU_LO, TAU_HI]`** — derive the band, don't magic it:
   `TAU_HI` = the precision of a single fragment at this node (`mass_i`, the Poisson cap — can't be more
   certain than your own count); `TAU_LO` = the global τ (a node is never *less* certain than the
   foundation). So the clip is data-derived, not a tunable.
2. **mass_i = 0** (empty node) → `f_g=0` already handled upstream (the `active` mask); skip.
3. **The `(1−w)` interaction.** Final count term precision = `(1−w_i)·τ_count(i)` — `w` yields to the
   strand, `τ_count` is the count's intrinsic reliability. Both apply (they answer different questions:
   "should the count speak here" vs "how precise is it"). Verify ss0.5 (`w=0` ⇒ full `τ_count`) and complex
   AMBIG (`w` high ⇒ count yields) both behave.
4. **σ²_global at zero-DNA.** ρ_global≈0 ⇒ DIRECT.predict(≈0) = the low-mean variance (small, the depleted
   seeds agree) ⇒ `σ²_global` small ⇒ `τ_global` high ⇒ tight pin at `μ_global≈0` ⇒ **no phantom**. This is
   the zero-DNA fix; assert it on `gdna_none`.
5. **DIRECT vs IMPUTATION at the same μ.** They're fit on different node subsets (observable depleted vs
   imputed enriched) so they live in different μ ranges; each is read only in its own regime — no
   cross-regime comparison. (Last turn's "ratio" confusion was exactly this; not a concern when each is
   used per-node-class.)

## Validation (single pass = Pass 0)

- Scoreboard fusion vs sweep (now var~mean precision instead of β): flagship ss0.99 / ss0.50 / zero-DNA /
  capture-off. Targets: ss0.50 ≤ fusion (the var~mean count precision should beat the flat β=10), zero-DNA
  phantom ≈ fusion (the per-node σ²_global pin), complex battery 2D/1D ≤ 0.95 preserved.
- Per-node `gdna_frac_var` populated from the lattice belief (the posterior `f_g` variance) — the per-node
  confidence (replaces the dropped guard), feeding diagnostics.
- Full suite green; `gdna+rna=total`.

## NOT in this increment

The all-gDNA **iteration** loop (re-fit var~mean on the solved estimates → re-solve → converge) is §8.4 —
this wires Pass 0. The boundary→region **prediction-error** convergence signal + cross-iteration stability
land with the loop.
