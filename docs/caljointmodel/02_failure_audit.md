# 02 — Failure Audit

This document justifies burning the legacy calibration system. The
audit is unchanged in substance from previous drafts — the legacy code
carries 91 hand-tuned numeric constants, ~25 of which are qualitative
cliffs, and produces the paralog phantom-RNA cascade described below.
The new design replaces all of it with the per-region sufficient
statistic model in [`01_generative_model.md`](01_generative_model.md).

## 1. The 91 magic numbers

Audit of `src/rigel/calibration/` (excluding KEEP files): **91 numeric
literals** acting as decision thresholds, prior strengths, damping
coefficients, fallback amplitudes, or schedule controls.

### Tier 1: Qualitative cliffs (~25)

Each flips a downstream quantity between qualitative regimes. Examples:

| Constant | Where | Cliff behavior |
|---|---|---|
| `rna_lower > 0.1` | `latent_states.py::build_logbf_expression` | Flips `p_unexpressed` between $e^{-5}$ and $\sim 0.5$. |
| `damping = 0.5` | `background_model.py` | Two-pass schedule with hand-tuned step. |
| `pass_index < 2` | `pipeline.py` | Different code path for first two passes. |
| `confidence_floor = 0.6` | `boundary_sweep.py` | Boundary evidence used or discarded. |
| `top_t_fraction = 0.01` | `bootstrap.py` | Seed mask trims top 1%. |
| `no_support_v_obs = 1e6` | `latent_states.py` | Artificial precision injection when evidence absent. |
| `mu_floor = 0.05` | `background_model.py` | gDNA mean clamped — phantom mass in empty regions. |
| `kappa_d_default = 0.5` | `_fl_mixture.py` | Hard fallback when fewer than `min_observations` gDNA fragments. |
| `min_observations = 50` | multiple | Below this, library defaults; above, regional estimates. |
| `fusion_weight = 0.7` | `latent_states.py` | Weight on strand vs density in fused soft label. |
| ... (15 more) ... | | |

Property: small data perturbations (one read added or removed) can
flip these gates, causing discontinuous downstream changes. This is
the architectural source of the paralog regression.

### Tier 2: Numerical guards (~30)

Floors, ceilings, and tolerances that should be set relative to
machine epsilon but in legacy code are hard-coded heuristics: `1e-3`
jitter values, `1e-6` Newton tolerances, `1e-12` log-clip floors.
Individually harmless; collectively fragile to input scale.

### Tier 3: Schedules and budgets (~36)

Iteration caps, batch sizes, EM pass counts, warmup schedules with
no Bayesian interpretation: `max_em_iters=15`, `warmup_passes=2`,
`chunk_size=1024`. Replaced in the new design by one outer iteration
cap and one mass-change convergence tolerance.

### Aggregate

| Tier | Count | New design |
|---|---|---|
| 1 (cliffs) | ~25 | **0** |
| 2 (guards) | ~30 | 4 documented `np.finfo` floors |
| 3 (schedules) | ~36 | 1 outer cap + 1 tolerance |
| **Total** | **91** | **≤ 8** |

## 2. Paralog read-by-read trace

### 2.1 Scenario

`tests/scenarios_aligned/test_paralogs.py`: t1 and t2 share an
internal multimap-only region with unique flanks. Simulated 500 RNA
fragments per transcript; expected output ~`(500, 500)`.

### 2.2 What the calibrator sees

- **No unique RNA in the multimap region.** Multimappers excluded by `bam_scanner.cpp:1347` (NH>1 filter) and `bam_scanner.cpp:1648`.
- **Pair-anchored gDNA straddlers do reach calibration.** One mate uniquely mapped in a flank, partner multi-mapping in the paralog → fragment recorded as NH=1. These are gDNA contamination by construction.
- Straddler strand split is noisy, e.g. 126 sense / 26 antisense in the failing seed.

### 2.3 What the legacy does

1. `background_model.py` estimates regional gDNA mean ~150 expected fragments.
2. `compartment_strand_deconv.fuse_density_and_strand` interprets 126/26 as RNA evidence (`κ_eff ≈ 0.83 > κ_d ≈ 0.5` → "sense excess must be RNA") → `rna_mass = 26`.
3. `latent_states.build_logbf_expression` sees `rna_lower = 26/150 ≈ 0.17 > 0.1` (Tier-1 cliff), flips `is_expressed_bf` high → `p_unexpressed = e^{-5}`.
4. `compute_locus_priors_from_partitions` reads `p_unexpressed = 0.007` and `rna_mass = 26` → 26 RNA pseudocount on the paralog region.
5. Locus EM converges to `(4, 638)` instead of `(500, 500)`.

### 2.4 Why this cannot be patched

Pushing the cliff `rna_lower > 0.1` up moves the cliff without removing
it. Damping fusion hides this case but breaks others. Setting
`rna_mass = 0` for unstranded fragments breaks stranded libraries. Each
fix introduces a new constant; no way out.

### 2.5 What the new system does

The paralog region is now rescued by three weaker but architecturally
clean channels firing together (FL is no longer a calibration channel
— see [`00_overview.md`](00_overview.md) §1 and
[`01_generative_model.md`](01_generative_model.md) §2):

1. **Count channel (NB).** The unspliced count $n_r^{\text{u}} \approx 152$
   is consistent with the gDNA-only mean
   $\mu_r^{(g)} = \hat{\omega}_r \rho_0 L_r^{\text{eff}} \approx 150$
   under the Gamma exposure posterior pooled across surrounding
   intergenic regions. Adding an RNA component $\mu_r^{(d)} > 0$
   pushes the NB mean off its mode for negligible likelihood gain.
   The count channel votes gDNA.
2. **Strand channel (BB).** A 126/26 split under
   $\mathrm{BetaBinom}(N=152, \kappa_d = 0.5, \rho_d^{\text{BB}})$
   with the library-fit $\rho_d^{\text{BB}}$ has finite tail mass at
   this asymmetry; the BB log-evidence ratio between the gDNA and
   RNA strand hypotheses is essentially zero. The legacy Binomial
   misread the same data as decisive RNA signal because it lacked
   the overdispersion term.
3. **Spliced channel.** Zero spliced fragments on the region;
   $n_r^{\text{s}} = 0$ contributes no RNA mass via the deterministic
   pathway.

Combined soft allocation: $M_r^{(g, \text{cont})} \approx 145$,
$M_r^{(d, \text{cont})} \approx 5$. The locus prior gets a near-symmetric
small RNA pseudocount across paralog transcripts. The EM splits
multimapper mass near-symmetrically. No cliff anywhere on this
pathway.

> **This rescue is empirically weaker than the prior
> (count + FL + strand) design.** The FL channel was the most
> decisive piece of evidence on paralog regions (the gDNA-characteristic
> short straddler-fragment distribution is a near-perfect
> discriminator). Without it, rescue depends on all three channels
> above firing together. Synthetic and scenario benchmarks must
> verify this rescue holds on production data before the new system
> replaces the legacy path — see
> [`06_validation_plan.md`](06_validation_plan.md).

## 3. Three architectural anti-patterns

### 3.1 Hard staging

Legacy pipeline has explicit stages (density estimation, boundary
scoring, fusion, deconvolution, posterior smoothing) with hand-off
through serialized intermediate objects (`BackgroundDensity`,
`BoundarySoftLabels`, `FusedSoftLabels`, `LatentStates`). Once a stage
finalizes a wrong answer, no later stage can recover.

**Replacement.** Single E-step that updates all per-region masses and
exposures simultaneously. Single M-step that updates all library
hyperparameters. No intermediate handoff objects.

### 3.2 Discriminative rules dressed as Bayes factors

Several legacy functions name their output `log_bf` but compute it as
a hand-tuned monotonic function of a feature, not as a ratio of
likelihoods under a stated model. Example: `build_logbf_expression`
computes `log_bf = sigmoid_steepness * (rna_lower - threshold)`.

**Replacement.** Every Bayes factor is derived from the per-region
likelihoods in [`01_generative_model.md`](01_generative_model.md) §5.

### 3.3 Asymmetric evidence treatment

Legacy code treats evidence *for* expression with one set of rules
(strand-fusion path) and evidence *against* with another (background
density floor). This is what allows pair-anchored gDNA straddlers to
be reinterpreted as RNA: the gDNA evidence is processed by
`background_model.py` and quietly floored, while the strand-based
"RNA evidence" goes through fusion without ever being asked to clear
a likelihood test against the gDNA hypothesis.

**Replacement.** The per-region likelihood is symmetric:
$\mathcal{L}_r(\pi_r^{(g)}, \omega_r)$ contains both hypotheses
under the same factorization. No hypothesis-specific code paths.

## 4. Burn list per file

KEEP / DELETE / ARCHIVE classifications. Definitive list confirmed by
the Phase 1 archive script via `git mv`.

| File | Classification |
|---|---|
| `__init__.py` | REWRITE (Phase 3) |
| `_arrays.py` | KEEP |
| `_categorize.py` | DELETE |
| `_fl_empirical_bayes.py` | ARCHIVE |
| `_fl_mixture.py` | ARCHIVE |
| `_result.py` | DELETE |
| `_simple.py` | DELETE (G1 algorithm lifted into new module) |
| `accumulator_view.py` | DELETE |
| `adaptive_prior.py` | DELETE |
| `background_model.py` | DELETE |
| `boundary_sweep.py` | DELETE |
| `boundary_model.py` | DELETE |
| `bootstrap.py` | DELETE |
| `compartment_strand_deconv.py` | DELETE |
| `coverage_weight.py` | DELETE |
| `density_model.py` | DELETE |
| `density_observation.py` | DELETE |
| `density_global.py` | DELETE |
| `exposure.py` | DELETE (replaced by `01_generative_model.md` §4.2) |
| `fl.py` | ARCHIVE |
| `fractional_evidence.py` | KEEP |
| `fusion.py` | DELETE |
| `latent_states.py` | DELETE |
| `locus_partition.py` | KEEP if used by `locus.py`; else DELETE |
| `prior.py` | DELETE |
| `region_count_ledger.py` | KEEP |
| `region_partition.py` | KEEP if used by `locus.py`; else DELETE |
| `regions.py` | KEEP |
| `scan_payload.py` | KEEP |
| `strand_deconv.py` | DELETE |

## 5. Tests to delete

~28 test files target deleted modules. Deleted outright, not adapted:

```
tests/test_adaptive_prior.py
tests/test_background_model.py
tests/test_bayesian_prior_acceptance.py
tests/test_boundary_model.py
tests/test_boundary_sweep.py
tests/test_calibrate.py
tests/test_calibration_accumulator.py
tests/test_calibration_integration.py
tests/test_calibration_iteration.py
tests/test_calibration_prior.py
tests/test_calibration_result.py
tests/test_compartment_strand_deconv.py
tests/test_coverage_weight.py
tests/test_density_global.py
tests/test_density_model.py
tests/test_density_observation.py
tests/test_exposure.py
tests/test_fl_eff_len_cache.py
tests/test_fl_models.py
tests/test_latent_states.py
tests/test_locus_partition.py
tests/test_region_persist.py
tests/test_region_unspliced_mass.py
tests/test_strand_*               (review case-by-case)
```

KEEP: `tests/test_region_count_ledger.py`, `tests/test_region_index_native.py` (substrate tests).

Replacement test suite in `tests/calibration/` per
[`06_validation_plan.md`](06_validation_plan.md).

## 6. Acceptance criteria for "the burn is done"

1. `git grep -nE 'background_model|fusion|latent_states|boundary_sweep|strand_deconv|adaptive_prior|p_unexpressed|fit_status|fused_soft_label' src/` returns zero hits.
2. `wc -l src/rigel/calibration/*.py` ≤ 600 lines post-rewrite.
3. Magic-number audit: ≤ 8 literals justified by comments citing this doc or the model doc.
