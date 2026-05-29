# Locoregional gDNA Density Redesign

**Date:** 2026-05-07
**Status:** Reviewed plan for locoregional density diagnostics. The
local-evidence prior-strength section is superseded by
[`bayesian_prior_plan.md`](../bayesian_prior/bayesian_prior_plan.md).
**Supersedes (in spirit):** the per-locus assembly path in
[`src/rigel/calibration/locus_prior.py`](../../src/rigel/calibration/locus_prior.py),
the EXON–INTRON branch in
[`src/rigel/calibration/density_global.py`](../../src/rigel/calibration/density_global.py),
and the C++ boundary accumulator in
[`src/rigel/native/calibration/accumulator.cpp`](../../src/rigel/native/calibration/accumulator.cpp).
Companion to [`boundary_flux_improvement.md`](./boundary_flux_improvement.md).

---

## 1. Motivation

The companion critique enumerates six modeling/robustness gaps in the current
SRD-v6 locoregional gDNA prior. They share a single root cause: the per-locus
assembly path mixes three responsibilities (FL-aware geometry, evidence
accounting, prior strength) inside one function and pairs **whole-region**
numerators with **clipped-locus** denominators. The redesign extracts FL-aware
exposure into its own primitive and decomposes the gDNA mass estimate into
explicit additive terms. The original idea of scaling Dirichlet ESS with
observed local evidence is retained only as a historical diagnostic path; it
should not drive the EM prior after the Bayesian prior redesign.

### Reviewer correction: keep units separate

The original draft correctly identified the boundary-flux problem, but it
mixed two different quantities:

- **raw boundary event count**: an observed fragment count; this can be added
   to `n_gdna_total` and compared to `n_obs`.
- **boundary density numerator / exposure-normalized contribution**: a rate
   estimator; this must be multiplied by an exposure before it becomes a
   fragment count.

Therefore the implementation MUST NOT replace `u_left`/`u_right` with only
`Σ 1/(ℓ−1)` floats. We need raw event counts for observed boundary mass,
evidence strength, diagnostics, and NB-style overdispersion. If we later want
Horvitz-Thompson (HT) localization, it should be stored as an additional
weighted channel, not as the only boundary-flux state.

## 2. Per-issue analysis & resolution

### Issue 1 — Local numerator/denominator geometry mismatch
*Symptom.* `_shrink_one_type` and `_shrink_exon_intron` sum
`count_col[rids]` (whole-region counts) but pair them with `leff` computed
from clipped spans (`max(start, locus.start) … min(end, locus.end)`). A
locus that touches only 5 kb of a 100 kb intron gets a numerator scaled to
100 kb against a denominator scaled to 5 kb — densities can inflate by a
large factor.

*Resolution.* Prorate per-region counts by the FL-aware exposure ratio:
```
eff_full_r,  eff_clip_r  = contained_exposure(start_r, end_r, clip_lo, clip_hi, fl)
ratio_r                  = eff_clip_r / eff_full_r       (0 if eff_full_r == 0)
n_loco_r                 = full_count_r × ratio_r
```
Boundary flux uses an analogous rule: an eligible side contributes to a
locus-local estimate when its boundary position lies in the locus window. If
we later add fractional side-window proration, it must be applied identically
to both the raw boundary-event numerator and the boundary-crossing exposure.
Numerator and denominator must always describe the same geometric object.

### Issue 2 — Wrong boundary denominator
*Symptom.* `density_global._density_exon_intron` uses
`(1_L + 1_R) × gdna_fl.mean` as the per-side exposure. For a fragment of
length ℓ, the number of start positions that **strictly** cross a boundary
is `ℓ − 1`, not `ℓ`. The current code is conceptually off by a +1 bp shift
and slightly density-deflating.

*Resolution.* Use the exact boundary-crossing exposure as the default
estimator:
```
B_cross = Σ_ℓ h(ℓ) × max(ℓ − 1, 0)
ρ̂_boundary = N_crossing / (n_sides × B_cross)
```
This keeps the numerator as an event count and fixes the `mean_FL` versus
`mean_FL - 1` geometry without changing the scanner schema.

Do not use an accumulator like `boundary_eff_len += (ℓᵢ - 1)` over only the
observed crossing fragments as the denominator. That quantity is conditional
on already having crossed the boundary and is therefore length-biased. For a
fixed fragment length it would give `N / (N × (ℓ - 1)) = 1 / (ℓ - 1)`, which
is independent of the true local gDNA density. It can be a useful diagnostic
for the observed crossing-fragment length mix, but it is not the physical
opportunity space. The opportunity space is either `n_sides × Σ h(ℓ)(ℓ - 1)`
from the gDNA FL model, or the HT density channel below.

The per-fragment Horvitz-Thompson estimator remains theoretically useful:
```
ρ̂_boundary_HT = (1 / n_sides) × Σᵢ 1 / (ℓᵢ − 1)
```
but it is a density estimator, not an observed-count field. If implemented,
store it as `u_left_ht`/`u_right_ht` alongside raw integer
`u_left`/`u_right`, and keep the raw counts as the source of
`n_gdna_boundary_observed`.

Variance / filtering note: if a `MIN_FL_FOR_BOUNDARY` threshold is used, the
denominator must also be truncated:
```
B_cross(min_fl) = Σ_{ℓ >= min_fl} h(ℓ) × max(ℓ − 1, 0)
```
For HT, filtering requires either division by `P(ℓ >= min_fl)` or projection
onto long-fragment-only exposure. Otherwise the estimator is biased low. The
recommended first implementation is therefore **no HT and no C++ length
filter**: just replace `mean_FL` with `B_cross` in Python.

### Issue 3 — Boundary mass conflated with imputed exon-contained mass
*Symptom.* Today the locus-prior multiplies ρ_boundary by the **full**
exonic L_eff in the locus. That predicts only the exon-contained gDNA mass.
The observed boundary fragments are never added back, but they *do* enter
`n_obs` (the EM-unit denominator of `π_gdna`). Result: `π_gdna` is biased
**downward** in loci where boundary fragments are a non-trivial share of
the EM units.

*Resolution.* Decompose `n_gdna` into four explicit additive terms, each
with its own numerator and exposure:
```
n_gdna_intron      = ρ_intron_loco       × L_intron_contained_in_locus
n_gdna_intergenic  = ρ_intergenic_loco   × L_intergenic_in_core_locus
n_gdna_boundary    = N_boundary_observed       ← raw observed boundary events
n_gdna_exon_only   = ρ_boundary_loco     × L_exon_contained
n_gdna_total       = sum of the above
```
`L_exon_contained` is the FL-weighted contained exposure for exon-only
fragments. Do **not** subtract boundary-crossing exposure from it: contained
exon-only fragments and boundary-crossing fragments are disjoint start-position
events. By construction, `n_gdna_total` and `n_obs` now describe the same
population of EM units, so `π_gdna = n_gdna_total / n_obs` is arithmetically
consistent.

### Issue 4 — Boundary anchor not RNA-free in nRNA-rich loci
*Symptom.* EXON–INTRON fragments can come from nascent RNA, not only from
gDNA. In nRNA-active loci the boundary anchor slightly overstates gDNA.

*Resolution (deferred, with telemetry).* Today, EM has explicit nRNA and
gDNA components and posteriors will arbitrate ambiguous boundary
fragments. We do not yet have evidence that the anchor bias materially
degrades π̂. Action this round:

- Add a `nrna_active` boolean to `LocusGdnaEstimate` (true when the
  multi-locus contains a synthetic-nRNA component).
- Surface it in the prior diagnostics so we can stratify benchmark
  results by it.

If post-deployment data shows a measurable bias, follow up by
down-weighting ρ_boundary in nRNA-active loci (or by iterating the
calibration with EM posterior gDNA expectations — a true two-stage SRD).

### Issue 5 — Locoregional intergenic branch is empty by construction
*Symptom.* `Locus` intervals are transcript-defined and almost never
contain intergenic regions, so the intergenic locoregional branch
collapses to ρ_global for nearly every locus. Local capture or
contamination patterns are invisible to the prior.

*Resolution.* For estimating **local background density**, expand the
intergenic query window by `INTERGENIC_FLANK_BP = 5000` (clamped to
`[0, ref_length]`). Intron and exon branches keep using the unflanked locus
interval.

Important: flanking intergenic fragments usually do **not** enter the EM unit
set for the transcript locus, so their predicted count must not be added to
`n_gdna_total` unless the core locus itself overlaps intergenic territory.
Use the flank to inform `ρ_intergenic_loco` and `n_evidence`; use only the
core-locus intergenic exposure for the mass term in `π_gdna`.

### Issue 6 — Fixed Dirichlet ESS erases evidence strength
*Symptom.* `α_gdna = c_base × π̂` makes a locus with 2 boundary events and
a locus with 200 boundary events produce identical prior strengths. Strong
local evidence is wasted.

*Supersession note.* This issue is real for the legacy `π̂ -> ESS` prior
family, but the recommended Bayesian redesign does **not** fix it by making
local evidence a stronger prior. The cleaner fix is to stop using in-locus
local evidence as prior mass at all. Local `n_evidence`, `c_loco`, and
`π̂_gdna` may remain useful diagnostics or warm-start inputs, but the EM prior
should use independent expected gDNA counts as described in
[`bayesian_prior_plan.md`](../bayesian_prior/bayesian_prior_plan.md).

*Former resolution, now deprecated for the EM objective.* Scale Dirichlet ESS
with observed evidence:
```
N_evidence = N_intron_loco + N_intergenic_loco + N_boundary_loco
c_extra    = min(C_LOCO_CAP - c_base,
                 BETA_EVIDENCE × N_evidence,
                 C_OBS_FRAC × max(n_obs, 1))
c_loco     = c_base + c_extra
α_gdna     = c_loco × π̂_gdna
α_rna      = c_loco × (1 - π̂_gdna)
```
Defaults (subject to the Phase 5 sweep): `c_base = 10.0` (unchanged for
zero-evidence loci → bit-identical legacy fallback), `BETA_EVIDENCE = 1.0`,
`C_LOCO_CAP = 1000.0`, `C_OBS_FRAC = 1.0`. The `C_OBS_FRAC` cap prevents
flanking or high-count evidence from producing a prior that is much stronger
than the locus's own EM likelihood. `c_loco` is recorded in
`LocusGdnaEstimate` for diagnostics.

---

## 3. Architecture

### 3.1 New module: `_exposure.py`

Single owner of FL-aware geometry. Contains:

| Function | Returns | Replaces |
|---|---|---|
| `contained_exposure_clipped(starts, ends, clip_lo, clip_hi, fl)` | `(eff_full[R], eff_clip[R])` | `l_eff_contained` (moved here, augmented to return both) |
| `boundary_crossing_exposure(fl, min_fl=2)` | `Σ h(ℓ) max(ℓ−1,0)` | `gdna_fl.mean` in boundary denominator |
| `boundary_side_in_window(boundary_pos, clip_lo, clip_hi)` | `bool` or `0/1` | (new — local side eligibility) |

The first implementation should keep per-fragment HT weighting out of the
C++ accumulator. If an HT diagnostic is later added, `_exposure.py` still
owns the denominator/truncation convention so count-based and HT estimates
remain comparable.

### 3.2 New schema: `LocusGdnaEstimate`

| Field | Meaning |
|---|---|
| `locus` | Source `Locus` object. |
| `n_obs` | EM units in the locus denominator. |
| `n_gdna_intergenic` | Core-locus intergenic gDNA mass; usually zero for transcript-defined loci. |
| `n_gdna_intron` | Intron-contained gDNA mass. |
| `n_gdna_boundary_observed` | Raw/prorated observed boundary event mass. |
| `n_gdna_exon_only` | Exon-contained gDNA mass imputed from boundary density. |
| `n_gdna` | Exact sum of the four `n_gdna_*` mass fields. |
| `pi_gdna` | `n_gdna / n_obs`, clipped to `[0, 1]`. |
| `rho_loco` | `(intergenic, intron, boundary)` local densities. |
| `leff_loco` | Exposure terms used to compute the mass fields. |
| `n_eligible_boundaries` | Eligible local exon-intron boundary sides. |
| `n_boundary_events` | Raw/prorated boundary event-count diagnostic. |
| `n_evidence` | Evidence count used to scale `c_loco`. |
| `c_loco` | Per-locus Dirichlet ESS. |
| `nrna_active` | True when the multi-locus contains a synthetic-nRNA component. |
| `fallback_flags` | Bit flags explaining global fallback behavior. |

The old compatibility aggregate `n_gdna_exon_intron` has been retired. Use the
explicit split fields `n_gdna_boundary_observed` and `n_gdna_exon_only`.

### 3.3 Boundary accumulator

Default implementation: keep [`accumulator.h`](../../src/rigel/native/calibration/accumulator.h)
and [`accumulator.cpp`](../../src/rigel/native/calibration/accumulator.cpp)
unchanged. `u_left`/`u_right` remain raw integer event counters.

Optional later HT diagnostic: add new `std::vector<double>` arrays
`u_left_ht`/`u_right_ht`. For each crossing fragment with `fl > 1`, add
`1 / (fl - 1)` to the matching HT side array while still incrementing the raw
integer side counter.
Do not drop raw event counts. `MIN_FL_FOR_BOUNDARY`, if introduced, must be
exported through the native module or summary machinery and must also be used
by the Python exposure denominator.

### 3.4 Density formula site

[`density_global.py:_density_exon_intron`](../../src/rigel/calibration/density_global.py)
becomes (essentially):
```python
n_cross = (u_left * bf_left + u_right * bf_right).sum()      # raw count
n_sides = (bf_left + bf_right).sum()
B_cross = boundary_crossing_exposure(gdna_fl)
denom   = n_sides * B_cross
rho     = n_cross / denom if denom > 0 else 0.0
```
No scalar `gdna_fl_mean` factor remains, but `gdna_fl` still matters through
the full FL-PMF exposure.

### 3.5 Per-locus assembly

`estimate_locus_gdna` should read top-to-bottom as this algorithm:

1. Query region IDs for the core transcript locus and for the
intergenic-flanked window.
2. Estimate intergenic density/evidence from the flanked window, but compute
`n_gdna_intergenic` from core-locus intergenic exposure only.
3. Estimate intron density and mass from core-locus intron exposure with
count proration.
4. Estimate boundary density from raw boundary event counts and boundary-side
crossing exposure, then compute both `n_gdna_boundary_observed` and
`n_gdna_exon_only`.
5. Compute `n_gdna_total = n_gdna_intergenic + n_gdna_intron +
n_gdna_boundary_observed + n_gdna_exon_only`.
6. Compute `n_evidence = n_intergenic_evidence + n_intron_evidence +
n_boundary_events`.
7. Compute `c_extra = min(c_loco_cap - c_base, beta_evidence * n_evidence,
c_obs_frac * max(n_obs, 1))` and `c_loco = c_base + c_extra`.
8. Return `LocusGdnaEstimate` and derive `alpha_gdna` / `alpha_rna` from
`c_loco` and `pi_gdna`.

Steps 7 and 8 are superseded for the production EM prior by the Bayesian prior
plan. Keep `n_evidence` and any legacy `c_loco` value as diagnostics only.

A single helper `_density_term_prorated` covers both INTRON and INTERGENIC.
`_boundary_term_prorated` is the only branch with both raw boundary-event
mass and exon-only imputed exposure.

### 3.6 Test plan

| Test file | Coverage |
|---|---|
| `tests/test_exposure.py` (new) | `contained_exposure_clipped` invariants: clip = full → ratio 1; clip ⊂ full → ratio = clipped_eff / full_eff; clip = ∅ → 0. `boundary_crossing_exposure` and `boundary_side_in_window` edge cases. |
| `tests/test_density_global.py` (update) | Boundary density: raw count divided by `n_sides × Σ h(ℓ)(ℓ−1)`. |
| `tests/test_per_locus_gdna_mass.py` (update) | New decomposition fields validated to sum to `n_gdna`; flanking intergenic evidence does not add mass outside the EM locus. |
| `tests/test_calibration_result.py` / golden outputs | Per-locus diagnostic columns using the split boundary/exon-only fields. |
| `tests/test_pipeline_integration_v6.py` | Update the invariant from `alpha_gdna + alpha_rna == c_base` to `== c_loco` / recorded ESS. |
| `tests/test_em_prior_weight.py` | Should pass unchanged (no API change to `build_prior_weight_rna`). |
| Synthetic-sim suite | Re-run after Phase 2, Phase 3, Phase 5; compare gDNA recovery and π̂_gdna calibration vs current baseline. |

## 4. Phasing

| Phase | Scope | Touches | Risk |
|---|---|---|---|
| **1** | `_exposure.py` module + tests | New file only | None |
| **2** | Boundary denominator fix (Issue 2), count-based MLE | density_global, tests | Low; no C++ change required |
| **3** | Locoregional proration (Issue 1) + decomposition (Issue 3) | locus_prior.py, schema | Largest schema change; gates on goldens |
| **4** | Intergenic flank as density/evidence only (Issue 5) | locus_prior.py, config/summary | Medium; must not add flanking mass to `π_gdna` |
| **5** | Legacy evidence-scaled ESS diagnostics only (Issue 6) | locus_prior.py, config/summary | Superseded as an EM prior by Bayesian prior plan |
| **6** | nRNA telemetry (Issue 4) — no model change | locus_prior.py field, summary.json | None |
| **7** | Optional HT diagnostic/channel | accumulator.{h,cpp}, scan_payload, density_global | Defer unless benchmarks show MLE is insufficient |

Each phase is one commit. Goldens regenerate at phases 2, 3, 5. Synthetic
suite re-runs after each of those.

## 5. Legacy defaults (superseded for EM prior strength)

The following constants describe the former local-evidence ESS proposal. They
should not be exposed as production prior-strength knobs after the Bayesian
prior redesign.

| Constant | Value | Rationale |
|---|---|---|
| `MIN_FL_FOR_BOUNDARY` | `2` initially | Avoid truncation bias. Raise only if denominator is truncated consistently. |
| `INTERGENIC_FLANK_BP` | `5000` | Captures ~1 mean fragment-length of intergenic context per side without dragging neighbouring loci. |
| `BETA_EVIDENCE` | `1.0` | One real evidence fragment ≈ one prior pseudocount. |
| `C_LOCO_CAP` | `1000.0` | Saturates at ~100× the legacy `c_base`; loci with thousands of fragments still let likelihoods dominate. |
| `C_OBS_FRAC` | `1.0` | Caps extra prior ESS by observed EM units in the locus. |
| `c_base` | `10.0` (unchanged) | Backward-compatible floor for zero-evidence loci. |

A single sweep — `BETA_EVIDENCE ∈ {0.25, 0.5, 1.0, 2.0}` against the synthetic
suite — fixes Phase 5 defaults empirically.

## 6. Risks & mitigations

- **Goldens regenerate multiple times.** Each phase is a separate commit so the
  diff is interpretable and bisectable.
- **Boundary units drift.** Raw counts, weighted HT sums, and densities have
   different units. Preserve raw counts as the mass/evidence quantity and keep
   HT weighted sums optional.
- **Schema/API churn.** The per-locus diagnostic schema changed to the split
   `n_gdna_boundary_observed` and `n_gdna_exon_only` fields; the old
   `n_gdna_exon_intron` aggregate has been removed.
- **C++ ABI change if HT is added.** Optional HT arrays require native
   payload/schema updates and recompilation. The default denominator fix avoids
   this risk.
- **Phase 5 ESS aggression.** A small locus with 1–2 evidence fragments
   could swing prior strength. The `+ c_base` floor, absolute cap, and
   observed-unit cap protect both ends; the sweep validates the multiplier.
- **HT variance for short fragments.** Deferred. If HT is implemented, length
   filtering must be mirrored in the denominator or the estimator is biased.

## 7. What this redesign delivers

1. One module owns FL-aware exposure; `locus_prior.py` reads as
   physics, not as flag-juggling.
2. `LocusGdnaEstimate` becomes self-explanatory — every field is a
   physical quantity, and the four `n_gdna_*` fields literally sum to
   `n_gdna`.
3. Numerator/denominator geometry is consistent at every density site
   (Issue 1 closed).
4. Boundary density estimator uses the correct FL-PMF crossing exposure
   rather than scalar `mean_FL` (Issue 2 closed without a native schema
   change).
5. Observed boundary mass is no longer dropped on the floor (Issue 3
   closed).
6. Local intergenic evidence is recoverable via the flank window
   (Issue 5 closed).
7. Local evidence strength is observable without being double-counted as EM
   prior mass (Issue 6 superseded by the Bayesian prior plan).
8. nRNA cross-anchor bias (Issue 4) is observable and ready for a
   follow-up if data warrants it.

The default implementation is mostly Python-side: exposure primitives,
boundary denominator replacement, per-locus decomposition, diagnostics, and
config/summary plumbing. Native C++ changes are deferred to the optional HT
diagnostic phase, so the first implementation can improve correctness without
an ABI/schema break in the scanner payload.
