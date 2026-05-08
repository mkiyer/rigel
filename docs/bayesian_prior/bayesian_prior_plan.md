# Bayesian EM Prior Redesign Plan

**Date:** 2026-05-07
**Status:** Verified implementation plan
**Scope:** Replace the current `pi_gdna -> c_base * pi_gdna` prior assembly
with a cleaner Bayesian prior derived from independent evidence counts.

## 1. Executive summary

The current prior path is over-parameterized because it compresses physical
gDNA count estimates into a fraction and then tries to recover prior strength
with tuning constants. After the latest locoregional changes, the concrete
path is now:

```text
n_gdna_estimate -> pi_gdna = n_gdna_estimate / n_obs
c_loco          -> c_base + clipped(beta_evidence * local_evidence)
alpha_gdna      -> c_loco * pi_gdna
alpha_rna       -> c_loco * (1 - pi_gdna)
```

That fraction bottleneck loses the scale of the evidence. A locus with 2
expected gDNA fragments and a locus with 2000 expected gDNA fragments can have
the same `pi_gdna`, so extra knobs are needed to reconstruct strength.

The deeper problem is that much of `n_gdna_estimate` is currently derived from
the same in-locus fragments that EM is about to process. Shrinkage toward a
global estimate reduces variance, but it does not make the local observations
independent. Using those same observations to build a prior and then feeding
them into EM double-counts evidence.

The recommended redesign is:

```text
prior pseudocounts = independently predicted gDNA count in the EM locus
```

In-locus local boundary/intron evidence should remain in the EM likelihood, or
be removed from EM and fixed as resolved data in a separate mode. It should not
also become prior mass.

This removes `c_base`, `BETA_EVIDENCE`, `C_LOCO_CAP`, and `C_OBS_FRAC` from
the gDNA prior derivation. The current locoregional mass estimate remains
valuable as a diagnostic, but its local evidence count should no longer set
the EM prior objective.

## 1.1 Verification addendum

The plan has been checked against the current codebase after the latest
locoregional changes.

Current verified state:

- `_exposure.py` now owns FL-aware containment and boundary-crossing geometry.
- `density_global.py` uses `boundary_crossing_exposure(gdna_fl)` for the
  EXON-INTRON denominator.
- `locus_prior.py` now records split local mass fields, `n_evidence`, and
  active `c_loco` evidence scaling.
- `PriorTable` now carries a `c_loco` array and tests assert
  `alpha_gdna + alpha_rna == c_loco`.
- The focused calibration/prior tests currently pass against that local-ESS
  implementation.

New implementation constraint discovered during verification:

`alpha_gdna` currently has two meanings in the native batch EM path. It is
both a prior-count increment and a gDNA component enablement / warm-start
switch. Specifically, the batch extractor disables the gDNA component when
`alpha_gdna <= 0`, and the warm-start override uses `alpha_gdna` when it is
positive. Therefore Phase 0 must decouple component eligibility and warm-start
behavior from prior mass before a true `alpha_rna = 0` global-only prior can
be considered clean.

Second verification correction: the current `prior_weight_rna` mechanism is a
prior-allocation weight, not a likelihood term. Under the new default
`alpha_rna = 0`, it intentionally becomes inactive. Do **not** move this
temporary nRNA-prior deweighting into the likelihood as part of the gDNA prior
redesign; that would replace one heuristic prior force with a heuristic score
force. Nascent-RNA-versus-gDNA resolution should be revisited after the
gDNA-versus-total-RNA foundation is made theoretically clean.

## 2. Design principle

For each multi-locus, separate the information sources:

| Symbol | Meaning | May contribute to prior? |
|---|---|---|
| `D_core` | Fragments routed to this EM multi-locus | No, unless removed/fixed before EM |
| `D_global` | Calibration evidence outside the specific locus, or effectively genome-wide | Yes |
| `D_flank` | Flanking evidence not routed to this EM multi-locus | Yes, if the window excludes `D_core` |
| `D_local` | In-locus boundary/intron evidence also present in EM | No; likelihood/diagnostic only |

The EM posterior should be:

```text
p(theta | D_core, D_external) proportional to p(D_core | theta) p(theta | D_external)
```

The prior `p(theta | D_external)` must be built from evidence external to the
EM data. If a fragment is in `D_core`, its role is the likelihood term unless
we explicitly choose a fixed-assignment/data-augmentation model.

## 3. Solver semantics: MAP-EM versus VBEM

The mathematical Dirichlet prior is usually written as:

```text
theta ~ Dirichlet(alpha)
```

For MAP estimation, the log prior contributes:

```text
sum_k (alpha_k - 1) log(theta_k)
```

For VBEM, the posterior Dirichlet parameters are approximately:

```text
alpha_post_k = alpha_prior_k + expected_count_k
```

Rigel's native solver currently treats the Python `alpha_gdna` and
`alpha_rna` arrays as **physical pseudocount increments**, not as raw
mathematical Dirichlet alpha values:

| Mode | Native baseline | Python value should mean |
|---|---:|---|
| MAP-EM | `0.0` | extra pseudocount increment, equivalent to `alpha - 1` |
| VBEM | `0.5` | extra pseudocount increment on top of Jeffreys baseline |

Therefore the implementation should pass an external expected count `eta_g`
as the current `alpha_gdna` value. It should **not** pass `1 + eta_g` to the
current MAP path. A later cleanup should rename these public fields to
`gdna_prior_count` and `rna_prior_count`, while preserving `alpha_*` aliases
for native ABI compatibility.

### 3.1 Native eligibility and warm-start semantics

The current batch EM path also uses `alpha_gdna` to decide whether the gDNA
component is eligible. That is not the desired long-term semantics. A zero
prior count should mean "no external gDNA pseudocount," not "remove gDNA from
the likelihood."

Preferred implementation:

1. Keep gDNA eligibility based on model/candidate availability, not on
    `alpha_gdna > 0`.
2. Let `alpha_gdna = 0` pass through as zero objective prior mass, with only
  the existing numerical floor used internally where needed for logs.
3. Do not let the gDNA warm-start override replace coverage-derived warm
  starts when `alpha_rna == 0`. With no RNA prior budget, the coverage
  warm start is the least biased initialization.

Concrete native ABI change:

```text
batch_locus_em_partitioned(...,
    locus_alpha_gdna: float64[n_loci],
    locus_alpha_rna: float64[n_loci],
    locus_enable_gdna: bool[n_loci],
    ...)
```

`locus_enable_gdna[i]` controls whether the gDNA component is eligible.
`locus_alpha_gdna[i]` controls only prior pseudocount mass. Existing callers
can initially pass `locus_enable_gdna = isfinite_any_gdna_likelihood_for_locus`
or an equivalent scorer/model-derived flag. A compatibility overload can be
kept briefly if needed, but the implementation should not infer enablement
from `alpha_gdna > 0`.

Fallback implementation, if we defer the native cleanup for one iteration:
pass a tiny objective-neutral enablement floor to native while storing the
diagnostic `gdna_prior_count` separately. This is less clean because the
native `alpha_gdna` array would no longer be an exact copy of the diagnostic
prior count, so the preferred native cleanup should be done first.

## 4. Recommended default: independent background prior

Use only external/background evidence to predict the number of gDNA fragments
that should appear in the EM locus.

For one locus, compute an independent prior count:

```text
eta_g = E[N_gdna_core | D_external]
```

Then pass:

```text
alpha_gdna = eta_g
alpha_rna  = eta_r
```

where `eta_r = 0` by default unless we have independent RNA evidence. Do not
set `eta_r = n_obs - eta_g`; that would reuse the EM data as prior evidence.
The RNA abundances are the primary quantities being inferred, so they should
usually be learned from the EM likelihood with only the mode-specific baseline
or a deliberately separate external RNA prior.

### 4.1 Global-only prior source

The simplest theoretically sound source is the global calibration density.
For each core locus, use global densities and core-locus exposures to compute:

```text
eta_g = rho_global_intergenic * L_core_intergenic
      + rho_global_intron     * L_core_intron
        + rho_global_boundary   * B_core_boundary
        + rho_global_boundary   * L_core_exon_contained
```

where:

```text
B_core_boundary = n_eligible_boundary_sides_core * boundary_crossing_exposure(gdna_fl)
```

This is independent in the practical sense that no in-locus boundary or intron
counts from the same EM units are used to set `eta_g`. If exact independence is
required, the global density can be made leave-one-locus-out by subtracting the
core locus's calibration counts and exposures from the global totals. In real
genome-scale data the difference should usually be negligible, but the
leave-one-locus-out form is the clean theoretical endpoint.

The helper that computes `eta_g` should project global densities onto the core
EM locus exposure only. It must not consult local `u_left`, `u_right`, intron
counts, intergenic flank counts, `n_evidence`, `pi_gdna`, or `c_loco`.

Implementation shortcut: this is the same exposure pipeline as the current
locoregional estimator with full shrinkage to global density, except that the
boundary-observed local count is replaced by its global expectation. In other
words:

```text
contained intergenic term = rho_global_intergenic * L_core_intergenic
contained intron term     = rho_global_intron     * L_core_intron
boundary crossing term    = rho_global_boundary   * B_core_boundary
exon-contained term       = rho_global_boundary   * L_core_exon_contained
```

This lets the implementation reuse `contained_exposure_clipped`,
`boundary_side_in_window`, and `boundary_crossing_exposure` instead of creating
a parallel geometry path. Do **not** reuse the current diagnostic `n_gdna`
field verbatim, because that field deliberately adds observed local boundary
events; the prior count must use expected boundary-crossing mass instead.

Leave-one-locus-out guard: if a core locus contributes a non-negligible share
of the global exposure for any density branch, compute that branch's global
density after subtracting the core locus contribution. A practical trigger is
`core_exposure / global_exposure > 0.01` for the relevant branch. Below that,
ordinary global density is acceptable.

### 4.2 Independent flanking prior source

If global density is too blunt for capture or regional contamination effects,
use a flanking window that excludes the core EM locus:

```text
D_flank = [locus_start - W, locus_start) union (locus_end, locus_end + W]
```

Estimate a local external density from `D_flank`, shrink that density toward
the global density, and project it onto the core exposure:

```text
rho_external = shrink(n_flank, L_flank, rho_global, kappa)
eta_g        = rho_external * L_core_gdna_opportunity
```

This keeps the locoregional idea while preserving data separation. The flank
window is not magic prior strength; it defines which observations are external
to the EM locus. The first implementation should be `global_only` because it
is the smallest theoretically clean foundation. It is appropriate as the
default for approximately uniform-contamination libraries, but it is not the
final production answer for capture protocols. Independent-flank priors are a
required follow-up before claiming capture-robust locoregional gDNA modeling.

### 4.3 What happens to in-locus local evidence?

In-locus boundary and intron fragments should still be measured and reported,
but not used as prior mass in the default path.

Their roles become:

1. EM likelihood evidence.
2. Diagnostics: `n_gdna_local_estimate`, `n_boundary_events`, `n_evidence`.
3. Optional warm-start information, if needed, but not objective-function
prior mass.

If boundary reads are truly gDNA-like, their likelihood should assign them to
gDNA or nRNA naturally. If they are ambiguous between gDNA and nRNA, that is
real model uncertainty; forcing them through a prior is not the right fix.

## 5. Alternative path: resolved-data augmentation

A second theoretically clean design is to remove unambiguous evidence from EM
before inference:

1. Identify fragments that are deterministic gDNA under the model.
2. Assign them directly to gDNA output counts.
3. Remove them from the EM unit set.
4. Run EM only on ambiguous fragments with a weak prior.

This avoids double-counting because resolved fragments are data, not prior.
However, this should not be the first implementation:

- EXON-INTRON fragments are mature-mRNA-free but not necessarily RNA-free;
  nascent RNA can generate them.
- Removing/fixing units changes output accounting, posterior assignment, and
  diagnostics more deeply than changing prior construction.
- It requires a clear definition of which fragments are truly deterministic.

Keep this as a later option for fragments that are genuinely outside the RNA
model, such as deterministic intergenic gDNA observations.

## 6. Proposed implementation shape

### 6.1 Replace the fraction bottleneck

Current assembly point:

```python
c_loco_ml = compute_c_loco(n_evidence_ml, ml_prior.n_obs, c_base=c_base)
alpha_gdna[idx] = c_loco_ml * ml_prior.pi_gdna
alpha_rna[idx] = c_loco_ml * (1.0 - ml_prior.pi_gdna)
```

New assembly point:

```python
alpha_gdna[idx] = ml_prior.gdna_prior_count
alpha_rna[idx] = ml_prior.rna_prior_count
```

`pi_gdna` remains a diagnostic, not the quantity that creates prior strength.
`c_loco` should no longer be the source of `alpha_gdna + alpha_rna`.

### 6.2 Split local estimate from prior evidence

The current `LocusGdnaEstimate` should be split conceptually into two parts:

| Field group | Meaning | Used for prior? |
|---|---|---|
| `n_gdna_local_*` | Same-locus estimate using local boundary/intron evidence | No |
| `n_gdna_prior_*` | Independent external expected gDNA mass | Yes |
| `pi_gdna_local` | Local diagnostic fraction | No |
| `gdna_prior_count` | External pseudocount increment | Yes |

For compatibility, existing columns can remain for one release, but new
diagnostic names should make the data provenance explicit.

The current `n_gdna` and `pi_gdna` fields should be treated as local diagnostic
estimates during migration. If they stay under their current names for one
release, add explicit `n_gdna_local` / `pi_gdna_local` aliases and make the
new prior fields canonical.

### 6.3 Default prior formula

For Phase 1, use global-only independent prior counts:

```text
gdna_prior_count = n_gdna_global_expected_core
rna_prior_count  = 0.0
```

If we later add an independent RNA prior, it must come from external RNA
evidence, not from `n_obs - gdna_prior_count`.

For native execution, `alpha_gdna` and `alpha_rna` should remain aliases of
these prior counts after the Phase 0 eligibility cleanup. Before that cleanup,
do not claim the implementation is fully Bayesian: any positive floor used
only to keep gDNA enabled must be recorded separately from the diagnostic
prior count.

### 6.4 Component allocation

Rigel currently passes one scalar `alpha_rna` budget and lets the native solver
distribute it across RNA components using coverage totals and
`prior_weight_rna`.

Under the recommended default, `alpha_rna = 0`, so this RNA allocation path is
effectively inactive. That is acceptable: transcript and nRNA abundances should
be inferred from the locus likelihood unless there is independent expression
evidence.

Longer term, if we want RNA priors, add a separate external RNA prior source
and distribute it transcript-centrically, never by gene-level summaries.

### 6.4.1 nRNA prior-weight policy

The current `prior_weight_rna` name is easy to misread. It is not a general
nRNA suppression model; it is only a multiplier applied when distributing a
positive `alpha_rna` prior budget across RNA components. Under the recommended
default `alpha_rna = 0`, this path should become a no-op.

Do not migrate `prior_weight_rna` into the E-step likelihood during this
redesign. Likelihood weights should represent evidence models, not temporary
regularization preferences. For now, the theoretically clean behavior is:

```text
gDNA prior: external expected gDNA count
RNA prior:  none
nRNA handling: likelihood/model evidence only
```

After this gDNA-versus-total-RNA foundation is stable, nascent RNA versus gDNA
resolution can be improved with explicit model terms, calibration diagnostics,
or independent RNA evidence. That work should not block this prior cleanup.

### 6.5 Output diagnostics

The current locus output column `gdna_prior` is computed as:

```text
alpha_gdna / (alpha_gdna + alpha_rna)
```

That interpretation breaks under the recommended default because
`alpha_rna = 0` does not mean "the prior says this locus is 100% gDNA." It
means "there is no independent RNA prior." Therefore implementation must
replace or deprecate `gdna_prior` in `loci_df` and golden outputs.

Recommended diagnostics:

| Field | Meaning |
|---|---|
| `gdna_prior_count` | External gDNA pseudocount passed to EM. |
| `rna_prior_count` | External RNA pseudocount, zero by default. |
| `gdna_prior_rate` | Optional diagnostic `gdna_prior_count / max(n_em_fragments, 1)`. |
| `pi_gdna_local` | Same-locus diagnostic estimate from locoregional evidence. |

For one release, keep the existing `gdna_prior` column but redefine it as a
rate diagnostic:

```text
gdna_prior = gdna_prior_count / max(n_em_fragments, 1)
```

Add explicit `gdna_prior_count` and `rna_prior_count` columns at the same time
so downstream analysis can migrate without guessing. Document the semantic
change in the manual and changelog; do not silently preserve the old
`alpha_gdna / (alpha_gdna + alpha_rna)` interpretation. Keep raw
`alpha_gdna` / `alpha_rna` only as solver/debug fields while the native ABI
still uses those names.

## 7. Configuration changes

Remove or deprecate these gDNA-prior-strength knobs:

```text
c_base
BETA_EVIDENCE
C_LOCO_CAP
C_OBS_FRAC
```

Replace them with provenance choices, not strength heuristics:

| Config | Default | Meaning |
|---|---|---|
| `gdna_prior_source` | `global_only` | Source for independent gDNA prior evidence. |
| `gdna_prior_flank_bp` | unused for `global_only` | Window size only for `independent_flank`. |
| `rna_prior_source` | `none` | External RNA prior source, if ever added. |

The only default needed for the first implementation is
`gdna_prior_source = global_only`. Even that can be internal at first rather
than exposed as a user-facing knob.

The current CLI and manual expose `--cal-c-base`. During migration, keep the
flag accepted for compatibility but mark it deprecated and no-op for
`gdna_prior_source = global_only`. Remove or rename user-facing documentation
that describes it as the primary prior strength control.

## 8. Migration phases

### Phase 0: Native semantics guardrail

- Decouple gDNA component eligibility from `alpha_gdna > 0` in the batch EM
  path.
- Add a separate `locus_enable_gdna: bool[n_loci]` native input.
- Preserve gDNA likelihood candidates whenever the scorer emits finite gDNA
  likelihoods and the locus model includes a gDNA component.
- Adjust the gDNA warm-start override so `alpha_rna = 0` keeps the
  coverage-derived warm start instead of forcing `theta_gdna = alpha_gdna`.
- Add a characterization test for the current behavior
  (`alpha_gdna = 0` disables gDNA), then update it in the same change that
  ships Phase 0 so the semantic change is reviewable.
- Add native tests proving `alpha_gdna = 0` with `enable_gdna = true` means
  zero prior mass, not gDNA component removal.
- Recompile after the C++ change.

### Phase 0.5: nRNA prior-weight cleanup

- Rename/document `prior_weight_rna` as an RNA-prior allocation weight, not a
  likelihood suppression term.
- Verify `alpha_rna = 0` makes `prior_weight_rna` objective-neutral.
- Do not move nRNA deweighting into the E-step likelihood in this redesign.
- Keep nRNA-versus-gDNA improvements as a follow-up modeling project after the
  gDNA-versus-total-RNA prior is corrected.

### Phase 1: Diagnostic and schema cleanup

- Add `gdna_prior_count` and `rna_prior_count` to `MultiLocusPrior` and
  `PriorTable` diagnostics.
- Add `prior_count_total` as the canonical total-prior diagnostic.
- Keep `alpha_gdna` and `alpha_rna` arrays for native compatibility, but
  document them as prior-count increments.
- Keep `pi_gdna` as a diagnostic only.
- Deprecate `c_base_value` and `c_loco` in `PriorTable` and summaries. If a
  transition alias is needed, make `c_loco` equal `prior_count_total`; do not
  keep local-evidence ESS under that name.
- Replace `loci_df.gdna_prior` with explicit prior-count diagnostics, or keep
  it as a deprecated compatibility column with a clearly documented new
  meaning.

### Phase 2: Global-only independent prior

- Add a helper that computes expected core gDNA mass using only global
  densities and core-locus exposure.
- Implement it by reusing the existing FL-aware exposure pipeline with the
  local-count contribution zeroed out / fully shrunk to global density. For
  boundary evidence, use expected `rho_global * B_core_boundary`, not observed
  local boundary events.
- Add the leave-one-locus-out guard for loci contributing more than 1% of a
  branch's global exposure.
- Do not use local `u_left`, `u_right`, intron counts, or intergenic flank
  counts in the prior path.
- Set `alpha_gdna = gdna_prior_count` and `alpha_rna = 0.0`.
- Retain local gDNA estimates in diagnostics for comparison.

### Phase 3: Tests and golden update

- Update tests that assert `alpha_gdna + alpha_rna == c_loco` or
  `alpha_gdna + alpha_rna >= c_base`.
- Add tests showing two loci with the same `pi_gdna` but different expected
  independent gDNA counts receive different prior counts.
- Add tests showing local boundary count changes do not alter the global-only
  prior when global density and exposure are fixed.
- Add MAP/VBEM tests documenting that Python prior values are increments over
  the native mode baseline.
- Add batch-native tests for the new gDNA eligibility behavior with
  `alpha_gdna = 0` and finite gDNA likelihoods.
- Add output-schema tests showing `gdna_prior` is no longer interpreted as
  `alpha_gdna / (alpha_gdna + alpha_rna)` under `alpha_rna = 0`.
- Before regenerating goldens, run a legacy-vs-new synthetic sanity check on
  the gDNA=0 plus nRNA-heavy scenarios. If global-only does not recover the
  Phase 5 regression, stop and root-cause before updating golden files.
- Regenerate golden outputs after intentional prior changes.

### Phase 4: Independent flank source for capture protocols

- Add `gdna_prior_source = independent_flank` behind a config flag.
- Compute density from flanking windows that exclude the EM core locus.
- Shrink flanking density to global density using the existing EB machinery.
- Project the external density onto core exposure to produce `gdna_prior_count`.
- Benchmark against `global_only` before changing defaults.
- Treat this phase as required before advertising capture-robust locoregional
  priors. `global_only` is the temporary default and the correct foundation;
  independent flank is the production path for non-uniform contamination.

### Phase 5: Optional resolved-data mode

- Define deterministic gDNA observations.
- Remove or fix those units before EM.
- Add their mass back to final gDNA outputs.
- Keep the EM prior weak/independent.

This phase is larger and should not block the cleaner prior derivation.

## 9. Validation strategy

Run three levels of validation.

### Unit invariants

- `gdna_prior_count` is nonnegative and derived without in-locus local counts
  in `global_only` mode.
- `pi_gdna` no longer determines `alpha_gdna` by multiplication with a fixed
  ESS.
- `alpha_rna` is zero unless an independent RNA prior source exists.
- `alpha_gdna = 0` does not disable gDNA likelihood eligibility after the
  native guardrail phase.
- MAP and VBEM mode baselines remain native-owned.
- `loci_df.gdna_prior` is removed, deprecated, or redefined so it cannot be
  misread as a mixture fraction when `alpha_rna = 0`.
- `prior_weight_rna` has no effect when `alpha_rna = 0` and is not applied to
  likelihood scores.
- Zero-data loci do not create observed fragment assignments from prior mass.
  If a native empty subproblem is evaluated with `gdna_prior_count > 0`, it may
  have a prior-only theta, but no reads should be assigned from that theta.

### Synthetic scenarios

- Boundary-rich loci: boundary reads should move gDNA/nRNA through likelihood,
  not prior force.
- Same `pi_gdna`, different locus size: prior counts scale with expected
  independent count, not fraction alone.
- nRNA-rich boundary loci: removing local boundary prior should avoid pushing
  all boundary evidence into gDNA before EM can arbitrate nRNA.
- Global gDNA contamination sweeps: global-only prior should provide stable
  baseline regularization without local double-counting.
- gDNA=0 plus nRNA-heavy regression scenarios: removing local ESS should
  recover the Phase 5 failures before any golden regeneration.

### Full benchmark comparison

Compare legacy, global-only, and optional independent-flank priors on:

- transcript mRNA error,
- nRNA siphon behavior,
- gDNA recovery at low contamination,
- posterior calibration of `gdna_prior_count` versus assigned gDNA count,
- sensitivity to MAP versus VBEM.

## 10. Expected code touch points

| File | Change |
|---|---|
| `src/rigel/calibration/locus_prior.py` | Add independent prior-count helpers; stop using `c_base * pi_gdna`. |
| `src/rigel/calibration/_result.py` | Add prior-count diagnostics; deprecate `c_base` summary. |
| `src/rigel/config.py` | Remove/deprecate `calibration.c_base`; optionally add `gdna_prior_source`. |
| `src/rigel/pipeline.py` | Stop threading `c_base` as prior strength; pass prior-source config if exposed. |
| `src/rigel/estimator.py` | Replace/redefine `loci_df.gdna_prior`; add prior-count diagnostics. |
| `src/rigel/cli.py` | Deprecate `--cal-c-base`; add provenance flags only if exposed. |
| `src/rigel/native/em_solver.cpp` | Required Phase 0 cleanup: decouple gDNA eligibility/warm-start from prior mass. |
| `src/rigel/native` bindings | Add/pass `locus_enable_gdna` separately from prior-count arrays. |
| `tests/test_pipeline_integration_v6.py` | Replace fixed ESS invariant. |
| `tests/test_calibrate_orchestrator.py` | Update `PriorTable.empty` / summary expectations. |
| `tests/test_em_impl.py` or native batch tests | Pin zero-prior gDNA eligibility semantics. |
| `tests/golden/` | Regenerate intentional prior/output changes. |
| `docs/MANUAL.md`, `docs/METHODS.md`, `CHANGELOG.md` | Remove or deprecate `c_base` as a prior-strength control; document `gdna_prior` semantic migration. |

## 11. Open questions

1. Should global density be exact leave-one-locus-out? The clean answer is
   yes, but the practical benefit may be tiny compared with genome-scale
   calibration evidence.
2. Should `global_only` include expected boundary-crossing gDNA mass as prior
   count? The principled answer is yes, as long as it is expected from global
   density rather than observed local boundary events.
3. Should RNA receive any prior budget? Default no. Add it only when we have
   independent transcript-level RNA evidence.
4. Should local evidence be used as warm start? Possibly, but only as an
   optimization aid. It must not change the MAP/VBEM objective.
5. Is a Dirichlet prior the best long-term form for gDNA? A hierarchical
   Gamma-Poisson contamination model may be more exact, but independent
   expected counts are the minimal compatible improvement for the current EM
   solver.
6. Should the transition keep `c_loco` as a deprecated alias of
  `prior_count_total`, or remove it from `PriorTable` immediately? Prefer the
  alias for one release if it lowers golden/schema churn, but do not keep the
  old local-evidence ESS under that name.
7. Should `independent_flank` use only intergenic flanks, or also external
  intron and boundary evidence outside the EM core? For capture protocols,
  intergenic evidence alone may be depleted, so the production flank design
  should consider all external gDNA-informative categories.

## 12. Recommendation

Proceed with the global-only independent prior first, but only after the
native eligibility guardrail is in place.

This gives Rigel a simple, defensible Bayesian prior:

```text
alpha_gdna = expected gDNA count from independent background evidence
alpha_rna  = 0 unless independent RNA evidence exists
```

The EM likelihood then owns the in-locus local reads. This removes the
fraction bottleneck, avoids double-counting local evidence, and deletes the
heuristic prior-strength constants rather than trying to tune them into
stability.

Do not add a softening multiplier such as `c * eta_g` in the first
implementation. The Bayesian interpretation is that `eta_g` is already the
external expected count, i.e. the prior sample size. If a small locus has
`n_obs = 5` and external evidence predicts `eta_g = 5`, then the prior being
roughly half of the total information is the correct propagation of external
evidence. If that is too strong in practice, fix the uncertainty model or the
external density estimate; do not reintroduce an unprincipled ESS knob.

Before starting implementation, revert or neutralize the Phase 5 local
evidence-scaled ESS behavior while preserving the Phase 3/4 geometry and
diagnostics. That gives the redesign a clean baseline: correct locoregional
measurement, no local evidence as prior force.
