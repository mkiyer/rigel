# Bayesian EM Prior Redesign Plan (v3)

**Date:** 2026-05-08
**Status:** Implementation-ready synthesis of v1 (engineering specifics) and
v2 (asymmetric-prior theory).
**Scope:** Replace the current `pi_gdna → c_base · pi_gdna` Dirichlet prior
with an asymmetric pseudocount prior derived from independent background
evidence. Decouple gDNA component eligibility from prior mass in the native
batch EM. Remove the `c_base`, `BETA_EVIDENCE`, `C_LOCO_CAP`, `C_OBS_FRAC`
heuristics.

---

## 0. Summary

```text
α_gdna  =  η_g       (expected gDNA count from external background evidence)
α_rna   =  0         (no informative RNA prior; let likelihood decide)
enable_gdna  =  finite gDNA likelihood candidates exist for the locus
                  (boolean, no longer derived from α_gdna > 0)
```

Everything else in this document is supporting machinery, migration choreography,
or guards against silent breakage.

---

## 1. Theoretical foundation

### 1.1 Asymmetric prior

The mixture's components are biologically asymmetric and must be treated
asymmetrically by the prior.

| Component | Nature | Predictable from background? | Appropriate prior |
|-----------|--------|-------------------------------|-------------------|
| **gDNA**  | Physical artefact: incomplete DNase digestion, capture off-target. Subject to genome-wide and regional density. | **Yes.** Background regions (intergenic, intron, flank) are genuine i.i.d.-ish samples of the same contamination process. | Informative pseudocount $\eta_g$ = expected count from background. |
| **RNA (mRNA + nRNA)** | Regulated biological signal. Locus-specific. Whether a transcript is on/off cannot be inferred from genome-wide rates. | **No.** | Uninformative ($\alpha_{rna} = 0$). Likelihood owns it. |

This is *not* a heuristic choice; it is the structurally correct prior
given the data-generating processes. A symmetric prior on RNA would mean
"a priori, every transcript is expected to be expressed at some non-zero
rate," which is biologically false.

### 1.2 Strict data partitioning

Partition each multi-locus's evidence into:

| Set | Definition | Role |
|-----|------------|------|
| $A_\ell$ | Fragments routed to the EM locus. | Likelihood. |
| $C_\ell$ | External evidence: global calibration counts, optionally a flanking window that excludes the EM locus. | Prior. |

The prior must be a function only of $C_\ell$. The current locoregional code
violates this by feeding in-locus boundary observations into both the prior
and the likelihood — that is the bug this redesign fixes.

This partitioning is a constraint of the current scalar Dirichlet-EM solver,
not a claim that boundary flux is biologically uninformative. In-locus
boundary-crossing fragments are legitimate evidence for local gDNA density,
and therefore for exon-contained gDNA, **if** the model can apply that
evidence only to the target exonic pool or remove/fix the boundary observations
before EM. The current solver has one locus-level `alpha_gdna`; using
in-locus boundary observations to set it would push on the same boundary units
that also appear in the likelihood. v3 therefore uses external evidence only.
A later stratified prior or full hierarchical model can use in-locus boundary
flux for exon-only fragments without double-counting.

### 1.3 Solver semantics

Rigel's native EM treats `alpha_gdna` / `alpha_rna` as **physical pseudocount
increments over the mode-specific Dirichlet baseline**:

| Mode  | Native baseline | Python-side $\alpha$ semantics |
|-------|-----------------|--------------------------------|
| MAP-EM | 0.0 (so MAP weight is $\alpha - 1$) | extra pseudocount on top |
| VBEM  | 0.5 (Jeffreys) | extra pseudocount on top |

Therefore $\eta_g$ is passed as `alpha_gdna` directly. No translation, no
rescaling, no `c_base` multiplier.

---

## 2. The prior formula

For each core locus $\ell$, using the existing FL-aware exposure pipeline
(`contained_exposure_clipped`, `boundary_crossing_exposure`,
`boundary_side_in_window`):

$$
\eta_g(\ell) \;=\; \rho^{\text{ig}} L^{\text{ig}}_\ell
            \;+\; \rho^{\text{in}} L^{\text{in}}_\ell
            \;+\; \rho^{\text{b}} \bigl( s_\ell \cdot B_{\text{cross}}
                                 \;+\; L^{\text{ex}}_\ell \bigr)
$$

where:

- $\rho^{\text{ig}}, \rho^{\text{in}}, \rho^{\text{b}}$ are global gDNA
  densities (intergenic, intron, exon-intron-boundary), produced by the
  existing `density_global` module.
- $L^{\text{ig}}_\ell, L^{\text{in}}_\ell, L^{\text{ex}}_\ell$ are
  FL-PMF-weighted contained effective lengths for each region type inside
  the locus's core (unflanked) interval.
- $s_\ell$ = number of eligible boundary sides inside the locus.
- $B_{\text{cross}} = \mathrm{boundary\_crossing\_exposure}(\text{gdna\_fl})$.

**Critical distinction from current code:** the boundary term uses
$\rho^{\text{b}} \cdot s_\ell \cdot B_{\text{cross}}$ — the *expected*
boundary-crossing mass — not the *observed* in-locus boundary count
$u_L + u_R$. This is what makes $\eta_g$ a function of $C_\ell$ alone.

**Implementation shortcut:** $\eta_g$ is mathematically equivalent to the
current `estimate_locus_gdna` evaluated at $\kappa = \infty$ (full shrinkage
to global), with the boundary-observed-count term replaced by its global
expectation. Reuse the existing geometry pipeline; do not build a parallel
exposure path.

### 2.1 Leave-one-locus-out guard

For mega-loci that contribute non-negligibly to global exposure, the global
densities are not strictly external. Add a per-branch guard:

```
if L_core_density[branch] / L_global_density[branch] > 0.01:
    use ρ_global_loo := (N_global − N_core) / max(L_global_density − L_core_density, ε)
else:
    use ρ_global
```

In practice this only fires for a handful of large gene clusters and
pericentromeric regions in human; it is a cheap correctness guarantee.

Implementation detail: exact LOO needs both global density-estimation totals
and core-locus density-estimation numerators/exposures. `GlobalDensityTable`
already stores global `N_global` and `L_global_density` per branch via
`n_fragments` and `eff_length_bp`, but `N_core` requires the calibration
payload. This is allowed only as a **negative subtraction** from the global
calibration pool; it must not be added as local evidence.

Keep the pure projection helper separate:

```python
expected_gdna_count_global(locus, region_index, region_arrays,
                           global_densities, gdna_fl) -> ExpectedGdnaPriorParts
```

Then layer exact LOO on top only when needed:

```python
expected_gdna_count_global_loo(..., payload_arrays, global_totals, threshold=0.01)
```

The LOO trigger uses density-estimation exposure for the relevant branch:
contained exposure for intergenic/intron, and boundary-crossing exposure
`s_l * B_cross` for the boundary density. The exon-contained projection
shares the boundary density; it does not define the LOO trigger denominator.

### 2.2 Why no softening factor

A reviewer will ask: "Should we not soften with $c \cdot \eta_g$ for some
$c < 1$?" The answer is no. $\eta_g$ already encodes the prior sample size
in the Bayesian sense ("as if I had observed $\eta_g$ external gDNA
fragments here"). For a small locus where $\eta_g \approx n_{obs}$, the
prior contributing roughly half the total information is **correct
propagation of independent evidence**, not over-regularization. If $\eta_g$
is empirically too strong, fix the density estimator or move to
leave-one-locus-out / flank-based sources — do not reintroduce a tuning knob.

---

## 3. Native ABI change (Phase 0 blocker)

The current batch EM overloads `alpha_gdna > 0` to mean **both** "give the
gDNA component prior mass" **and** "make the gDNA component eligible." This
must be split before $\alpha_{rna} = 0$ ships, otherwise we silently disable
gDNA wherever the prior happens to be zero.

### 3.1 New native signature

```cpp
batch_locus_em_partitioned(
    ...,
    const double*  locus_alpha_gdna,    // prior pseudocount, may be 0
    const double*  locus_alpha_rna,     // prior pseudocount, may be 0
    const uint8_t* locus_enable_gdna,   // NEW: 0/1 component eligibility
    ...);
```

### 3.2 Semantics

| Field | Old meaning | New meaning |
|-------|-------------|-------------|
| `alpha_gdna` | Pseudocount **and** eligibility switch **and** warm-start override target. | Pseudocount only. |
| `enable_gdna` | (implicit, derived from `alpha_gdna > 0`). | Explicit 0/1 flag. Set by Python from "any finite gDNA likelihood candidate exists for this locus." |

Use `uint8_t` / `np.uint8` rather than C++ `bool` storage in the native ABI;
this matches the existing native conventions for boolean per-unit arrays and
avoids packed-bool surprises.

### 3.3 Warm-start rule

When `alpha_rna == 0`, do **not** override the coverage-derived warm start
with `theta_gdna := alpha_gdna`. The coverage warm start is the least biased
initialization in the no-RNA-prior regime.

Concretely, `compute_ovr_prior_and_warm_start` should keep using
coverage-derived `theta_init_out` unless both `alpha_gdna > 0` and
`alpha_rna > 0` make the existing gDNA/RNA prior-ratio override meaningful.
`enable_gdna` controls eligibility before this function; it is not itself a
warm-start strength.

### 3.4 Reviewable change procedure

1. Add a characterization test that pins the **current** behavior
   (`alpha_gdna = 0` ⇒ gDNA disabled).
2. In the same PR that ships the new ABI, flip that test to assert the new
   behavior. The diff makes the semantic change explicit.

---

## 4. Out of scope (explicit)

Three things are deliberately *not* fixed by this redesign. Calling them
out prevents scope creep and bad reviews:

1. **nRNA-vs-gDNA arbitration.** `prior_weight_rna` becomes objective-neutral
   under $\alpha_{rna} = 0$. This is intentional. nRNA siphon is a separate
   modeling problem (likelihood / model evidence / spliced-vs-unspliced
   features), and shoehorning a "likelihood weight" hack to preserve nRNA
   suppression would replace one heuristic with another. It will be
   addressed in a follow-up project after this foundation is stable.
2. **Resolved-data augmentation** (deterministic gDNA fragments removed from
   EM, added back to outputs). Theoretically sound but a much larger
   accounting change. Not in scope.
3. **Hierarchical Gamma-Poisson contamination model.** A more exact prior
   form, but incompatible with the current Dirichlet-EM solver. Out of scope.

---

## 5. Migration phases

Each phase is independently testable and bisectable. Goldens are regenerated
only after Phase 3.

### Phase 0 — Native eligibility decoupling

- Add `locus_enable_gdna: uint8[n_loci]` to `batch_locus_em_partitioned`
  bindings and the C++ implementation.
- Replace the extractor's gDNA eligibility check with the new flag. Keep
  `alpha_gdna > 0` only where it truly means "positive gDNA prior mass."
- Update warm-start logic per §3.3.
- Recompile (`pip install --no-build-isolation -e .`).
- Tests:
  - Pin old behavior, then flip assertion in same PR (§3.4).
  - New: `alpha_gdna = 0, enable_gdna = True, finite gDNA scores` →
    EM emits non-trivial gDNA assignments.
  - New: `enable_gdna = False` → no gDNA mass regardless of `alpha_gdna`.

### Phase 0.5 — `prior_weight_rna` policy lock

- Document `prior_weight_rna` as an "RNA-prior allocation weight." Add a
  module docstring stating it is objective-neutral when `alpha_rna = 0` and
  must not be ported into the likelihood as part of this redesign.
- Add a unit test asserting that varying `prior_weight_rna` with
  `alpha_rna = 0` produces bit-identical EM output.
- No code-path changes.

### Phase 1 — Schema + diagnostics

- New fields on `MultiLocusPrior`:
  - `gdna_prior_count: float`  (canonical, = $\eta_g$)
  - `rna_prior_count: float`   (canonical, = 0 in this redesign)
- New array on `PriorTable`:
  - `gdna_prior_count: float64[n_ml]`, `rna_prior_count: float64[n_ml]`
  - `enable_gdna: uint8/bool[n_ml]`
- Demote `c_loco` array to a deprecated alias of `gdna_prior_count +
  rna_prior_count` for one release; remove `c_base_value` from new summaries.
- Keep `pi_gdna` as a *diagnostic* on `LocusGdnaEstimate`, computed from the
  legacy locoregional pathway, but never consumed by prior assembly.
- Introduce `pi_gdna_local` as the explicit alias for downstream code that
  needs to migrate.

### Phase 2 — Global-only prior helper

- New helper `expected_gdna_count_global(locus, region_index, region_arrays,
  global_densities, gdna_fl, *, ref_length=None) → ExpectedGdnaPriorParts`
  implementing §2. Reuses `_exposure.py` primitives. Pure function; no
  `payload_arrays` dependency (proves at the type level that local counts
  cannot leak into the ordinary global-only prior).
- `ExpectedGdnaPriorParts` should carry at least:
  - `total`
  - `intergenic_contained`
  - `intron_contained`
  - `boundary_crossing_expected`
  - `exon_contained_expected`
  - `density_exposure_intergenic`
  - `density_exposure_intron`
  - `density_exposure_boundary`
- Add an optional leave-one-locus-out wrapper (§2.1). Implementation:
  pre-compute global $N$ and $L$ totals once from `GlobalDensityTable`; when
  the threshold trips, compute the core $N$ terms from `payload_arrays` only
  for subtraction from the global calibration pool.
- Add `enable_gdna_for_multilocus(ml, em_data) -> bool`:

  ```python
  units = ml.unit_indices
  return bool(np.any((~em_data.is_spliced[units]) & np.isfinite(em_data.gdna_log_liks[units])))
  ```

  This must run before `partition_and_free`, while `em_data` still owns the
  global per-unit arrays.
- Wire into `assemble_priors`:

  ```python
  for ml in multi_loci:
      eta_g = sum(
          expected_gdna_count_global(loc, ...) for loc in ml.loci
      )
      gdna_prior_count[idx] = eta_g
      rna_prior_count[idx]  = 0.0
      enable_gdna[idx]      = enable_gdna_for_multilocus(ml, em_data)
      alpha_gdna[idx]       = eta_g
      alpha_rna[idx]        = 0.0
  ```

- Delete `compute_c_loco`, `BETA_EVIDENCE_DEFAULT`, `C_LOCO_CAP_DEFAULT`,
  `C_OBS_FRAC_DEFAULT` and their `__all__` entries. Mark `C_BASE_DEFAULT`
  deprecated in docs and CLI behavior; do **not** emit a warning merely on
  import, because many tests and downstream modules import constants.

### Phase 3 — Output schema + tests + goldens

- `loci_df.gdna_prior` redefined for one release as
  `gdna_prior_count / max(n_em_fragments, 1)` (a *rate*, not a fraction).
  Add new explicit columns `gdna_prior_count`, `rna_prior_count`.
  Document the redefinition in `CHANGELOG.md` and `docs/MANUAL.md`.
- Update tests:
  - Drop `alpha_gdna + alpha_rna == c_loco` invariant in
    `test_pipeline_integration_v6.py`. Replace with
    `alpha_rna == 0 and alpha_gdna >= 0`.
  - Drop `c_base` references in `test_calibrate_orchestrator.py`,
    `test_calibration_result.py`, `test_assemble_priors.py`.
  - Add: two loci with the same legacy `pi_gdna_local` but different
    geometric exposures get different `gdna_prior_count`.
  - Add: changing in-locus `u_left`/`u_right` does not change
    `gdna_prior_count` for the pure global helper, or for loci below the
    LOO trigger, when global densities and exposures are fixed (the whole
    point of the ordinary global-only path).
  - Add: $n_{obs} = 0, \eta_g > 0$ → no fragments are assigned (zero data,
    zero E-step responsibility); `theta_gdna` may be non-zero in
    parameter space but no read is generated.
- **Acceptance gate** (before regenerating goldens): run a synthetic sweep
  on the Phase 5 regression scenarios (gDNA=0 + nRNA-heavy). Global-only
  prior must recover those scenarios. If it does not, stop and root-cause —
  the theory predicts it should, so failure indicates a deeper issue.
- Regenerate goldens.

### Phase 4 — Independent-flank prior (production destination)

`global_only` is the correct *foundation*, but it assumes uniform
contamination across the genome. For capture protocols and CNV-rich regions
this is wrong. This phase is **not optional for production**; it is the
shipping default for capture-protocol users.

- Add `gdna_prior_source: Literal["global_only", "independent_flank"]` to
  `PipelineConfig` (default still `global_only` until benchmarked).
- For `independent_flank`:
  - Define $D_{\text{flank}} = [\ell_{start} - W, \ell_{start}) \cup
    (\ell_{end}, \ell_{end} + W]$, clamped to ref length, $W$ default 5000.
  - Subtract overlap with neighboring loci from $D_{\text{flank}}$.
  - Estimate $\rho_\ell^{\text{flank}}$ per branch from $D_{\text{flank}}$
    counts and exposures.
  - Empirical-Bayes shrink $\rho_\ell^{\text{flank}}$ toward $\rho^{\text{global}}$
    using existing `shrink_to_loco`.
  - Project shrunken density onto core exposure as in §2.
- Open question: should flank evidence include only intergenic, or also
  intron and external boundary? For capture protocols intergenic alone may
  be depleted. Resolve via benchmark.
- Benchmark sweep: legacy vs `global_only` vs `independent_flank` on the
  full synthetic suite; promote whichever wins on capture-protocol
  scenarios to default.

---

## 6. Validation strategy

### Unit invariants

| Invariant | Test |
|-----------|------|
| `gdna_prior_count` independent of in-locus boundary observations in ordinary `global_only`. | Vary `u_left`, `u_right` below the LOO trigger; expect identical `gdna_prior_count`. |
| `gdna_prior_count` scales with exposure, not fraction. | Two loci, same `pi_gdna_local`, different sizes → different counts. |
| `alpha_rna == 0` always under `global_only`. | Trivial. |
| `alpha_gdna = 0, enable_gdna = True` keeps gDNA component active. | Native batch test. |
| `prior_weight_rna` is objective-neutral when `alpha_rna == 0`. | Sweep weights, expect bit-identical θ. |
| Zero-data loci do not produce read assignments from prior. | $n_{obs} = 0, \eta_g > 0$ → empty assignment matrix. |
| `loci_df.gdna_prior` is a nonnegative rate not a fraction. | Type/finite/nonnegative test; do not require `<= 1`. |
| LOO guard fires above 1% exposure share. | Construct synthetic mega-locus. |

### Synthetic regression scenarios (block goldens until passing)

- **Phase 5 regression set:** gDNA=0 + nRNA={30,70} + ss={0.9, 1.0}. Total
  RNA recovery error must drop back below 5% (currently 5–27% under Phase 5).
- **Boundary-rich loci:** boundary fragments must move via likelihood, not
  prior. Verify by varying `u_left`/`u_right` and observing no change in
  prior, but expected change in posterior.
- **gDNA contamination sweep** (`gdna_none/low/high` × `ss_{0.5,0.9,1.0}`):
  global-only prior should provide stable baseline regularization.
- **Same legacy $\pi_{gdna}$, different locus size:** prior counts diverge.

### Full benchmark

Compare on the Armis2 unified-benchmark conditions:

- mRNA / nRNA / gDNA relative error per condition,
- nRNA siphon magnitude (NTA1/TA1 ratio sweep),
- gDNA recovery at low contamination (`gdna_low`),
- posterior calibration: $\eta_g$ vs realized assigned gDNA count
  (Pearson correlation, slope, intercept),
- MAP vs VBEM sensitivity.

---

## 7. Code touch points

| File | Change |
|------|--------|
| `src/rigel/native/em_solver.cpp` | Add `locus_enable_gdna`; decouple eligibility from `alpha_gdna > 0`; fix warm-start when `alpha_rna == 0`. |
| `src/rigel/native.py`, nanobind bindings | Pass new array. |
| `src/rigel/calibration/locus_prior.py` | New `expected_gdna_count_global`; rewrite `assemble_priors`; delete `compute_c_loco`; deprecate `C_BASE_DEFAULT`. |
| `src/rigel/calibration/_result.py` | Add `gdna_prior_count`, `rna_prior_count` columns; deprecate `c_base` in summaries. |
| `src/rigel/calibration/density_global.py` | Add LOO subtraction helper. |
| `src/rigel/config.py` | Add `gdna_prior_source` (default `global_only`); deprecate `c_base`. |
| `src/rigel/estimator.py` | Redefine `loci_df.gdna_prior` as rate; add new columns. |
| `src/rigel/cli.py` | `--cal-c-base` accepted but warns and is no-op. |
| `tests/test_pipeline_integration_v6.py` | Drop `c_loco` invariant. |
| `tests/test_assemble_priors.py` | Rewrite for new schema. |
| `tests/test_calibration_result.py` | Drop `c_base`/`c_loco` from `PriorTable` constructions. |
| `tests/test_calibrate_orchestrator.py` | Same. |
| `tests/test_em_impl.py` | Add native eligibility tests. |
| `tests/test_per_locus_gdna_mass.py` | Split: legacy local-mass tests stay (diagnostic); add new prior-count tests. |
| `tests/golden/` | Regenerate after Phase 3 acceptance gate passes. |
| `docs/MANUAL.md`, `docs/METHODS.md`, `CHANGELOG.md` | Document new prior, deprecate `c_base`, document `gdna_prior` semantic migration. |

---

## 8. Pre-implementation cleanup

Before starting the Python prior rewrite, **revert or neutralize Phase 5**
(the `compute_c_loco` evidence scaling) so the diff for this redesign is
clean. Phase 0 can land before or after this cleanup because it changes only
native eligibility semantics.
Keep Phases 1–4 of the locoregional plan: the FL-aware geometry
(`_exposure.py`, `B_cross`, intergenic flank) is correct and reused by §2.

Concrete revert (one-line equivalent):

```python
# In assemble_priors, restore:
alpha_gdna[idx] = c_base * ml_prior.pi_gdna
alpha_rna[idx]  = c_base * (1.0 - ml_prior.pi_gdna)
c_loco_arr[idx] = c_base
```

This restores the Phase 4-style baseline and a clean canvas for the redesign.
Run the focused calibration/prior tests after the neutralization before
claiming the baseline is green.

---

## 9. Open questions (track but do not block)

1. **LOO threshold.** 1% is a heuristic. Sweep on real mega-loci and tune.
2. **Boundary mass in $\eta_g$.** Included per §2 (expected, not observed).
   Some readers will argue exonic gDNA is always observable as boundary
   crossings only and so the exon-contained term double-counts. The current
   formulation treats them as additive contributions to *expected* mass and
   relies on the exposure pipeline to be FL-correct. Verify against
   simulation.
3. **`independent_flank` evidence categories** (§5 Phase 4): intergenic only,
   or all external? Defer to benchmark.
4. **RNA prior, ever?** Default no. Only revisit if independent transcript
   abundance evidence (e.g. matched orthogonal assay) is available.
5. **`c_loco` deprecation timing.** Alias for one release vs immediate
   removal. Prefer alias to soften downstream impact.

---

## 10. Recommendation

Proceed with v3. The asymmetric pseudocount prior is the structurally
correct treatment of the gDNA-vs-RNA distinction; everything else in this
plan is migration discipline. The prior formula is two lines of code on top
of geometry that already exists. The expensive parts are the native ABI
change (Phase 0) and the discipline to *not* fold nRNA suppression into the
likelihood as a side effect.

Acceptance order:

1. Phase 0 lands and is bisectable on its own.
2. Phase 0.5 freezes the `prior_weight_rna` semantics so nothing regresses.
3. Phases 1–3 land together; goldens regenerate only after the Phase 3
   acceptance gate passes.
4. Phase 4 lands as a follow-up driven by the Armis2 benchmark sweep on
   capture-protocol simulations.
