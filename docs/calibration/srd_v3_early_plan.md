# SRD v3 — Early Plan / Roadmap

**Date:** 2026-04-28
**Status:** Roadmap (no code change yet). Phase 1 implementation plan
lives in `srd_v3_phase1_strand_relabel.md`.
**Audience:** Future contributors; future-Claude pickup point.

---

## 1. Why v3

SRD v2 shipped a working 1-D fragment-length mixture for gDNA
calibration. It is sound for what it is, but it has one structural
limitation that everything else inherits:

> **The 1-D FL mixture cannot separate nascent RNA from gDNA when
> their fragment-length distributions overlap** (which they always do
> in practice). `srd_v2_phase2plus_handoff.md` §7a documented this as
> an "accepted identifiability limit". It is not a true identifiability
> limit — it is the limit of the chosen statistic. Adding a strand
> axis breaks the degeneracy whenever library strand specificity is
> well separated from 0.5.

SRD v3 is the natural next step: keep the v2 architecture (single-pass
C++ scanner → columnar buffer → Python calibration), but augment the
calibration model with strand information **expressed in the
transcript frame of reference** so that nRNA and gDNA become
distinguishable in the INTRONIC and EXON_INCOMPATIBLE buckets.

## 2. Summary of carryover items from v1 / v2

### 2.1 Closed (no action needed)

| Item | Disposition | Reference |
|---|---|---|
| v2 Phase 5 cleanup (intergenic accumulator removal) | Done | Phase 6 handoff §1 |
| Q1 version revert 0.6.0 → 0.4.1 | Done | Phase 6 handoff §1.1 |
| Q2 intergenic FL into global histogram | Done | Phase 6 handoff §1.1 |
| Q5 nRNA-not-in-cgranges architectural fix | Done | `nrna_derive_plan.md` |
| Q3 chimera leak plan (Option B) | **REJECTED** | `srd_v2_chimera_leak_plan.md` (banner) |
| OOR drop policy (drop FL>max from training) | Done — Phase 7 | `_simple.py:155-181`, `_result.py` |
| `n_pool` denominator tightened to in-range | Done — Phase 7 | `_simple.py` |

### 2.2 Open / carried into v3 scope

| Item | Note |
|---|---|
| nRNA siphon (v2 §7a) | **Primary motivation for v3.** Joint (FL × strand) mixture resolves it. |
| SS-gated antisense-exonic pool inclusion (v2 final-plan §Pillar 4) | Subsumed by v3: with proper transcript-frame strand labels, the EXON-incompat / contained antisense subset becomes one of the natural data sources. |
| INTRONIC POS/NEG strand asymmetry diagnostic | Replaced wholesale: POS/NEG in genomic frame is meaningless; SENSE/ANTISENSE in transcript frame is informative. |
| Phase 4 dna40m + dna80m sweep completion | Carryover ops task; do once v3 Phase 1 lands so the new diagnostic is captured. |
| FL mixture iteration cap (>1000 at dna00m) | Performance tuning; not a v3 prereq. Re-evaluate after the joint mixture changes the convergence landscape. |
| `pi_pool` wrong-direction monotonicity in VCaP mixture sweep | Unresolved in v2; suspected to evaporate once nRNA component is split out. Re-test after v3 Phase 2. |
| Per-block preservation for SPLICE_ARTIFACT | Deferred indefinitely — held out of calibration; downstream doesn't need it. |
| Mappability-corrected effective length / cross-region density | Out of v3 scope. Separate effort. |

### 2.3 Documentation debt

| Doc | What needs updating after v3 phases land |
|---|---|
| `docs/MANUAL.md`, `docs/METHODS.md`, `CHANGELOG.md` | Describe v3 calibration; cite this roadmap and the phase plans. |
| `srd_v2_phase2plus_handoff.md` §7a | Mark resolved with link to v3 Phase 2. |
| `srd_v2_results.md` "Long-standing artifact" section on INTRONIC=0 | Re-explain in v3 frame: INTRONIC=0 was never about antisense — it was about the strict zero-bp definition. |

## 3. Theoretical foundation for v3

### 3.1 The two strand operations

Two distinct operations on strand information are conflated in v2 code.
v3 separates them:

**Operation A — Frame transformation (deterministic).**
Given a fragment overlapping a transcript on strand `T` and a read
aligned to genomic strand `R`, define:
- `R == T` → **sense** to the transcript
- `R != T` → **antisense** to the transcript

This is a coordinate change. No probabilities. v2 skips it
entirely — `_categorize._strand_label()` reports "which transcript
strand has overlap", which is meaningless: it depends on which
genomic strand the local gene happens to live on, not on what the
fragment is doing relative to that gene.

**Operation B — Probabilistic interpretation (stochastic).**
Once a fragment is labeled in transcript frame, the strand-specificity
parameter `SS` learned by `StrandModel` from spliced unique mappers
gives source-conditional likelihoods:

| Source | P(sense \| source) | P(antisense \| source) |
|---|---|---|
| mRNA | SS | 1 − SS |
| nRNA | SS | 1 − SS |
| gDNA | 0.5 | 0.5 |

mRNA and nRNA share the same SS form (both arise from the same
underlying transcription orientation); gDNA is symmetric by
construction (sonicated dsDNA ends up on either genomic strand
50/50). When SS is well separated from 0.5, the strand axis has
direct discriminative power between gDNA and (m/n)RNA components.

### 3.2 Identifiability of the v3 mixture

In any pool bucket where multiple sources contribute, label every
unique-mapper UNSPLICED fragment with `(L, S)` where `L ∈ [0, max_size]`
is FL and `S ∈ {sense, antisense}`. The likelihood for a 3-component
mixture (gDNA + mRNA + nRNA) becomes:

$$
p(L, S \mid \pi_g, \pi_m, \pi_n) =
\pi_g \cdot p_{\text{gDNA-FL}}(L) \cdot 0.5
+ \pi_m \cdot p_{\text{mRNA-FL}}(L) \cdot p_S(\text{SS})
+ \pi_n \cdot p_{\text{nRNA-FL}}(L) \cdot p_S(\text{SS})
$$

with $p_S(\text{SS}) = \text{SS}$ if $S = $ sense and $1 - \text{SS}$
otherwise. Per pool bucket we collapse this further depending on
which sources are physically possible (see §4).

The mRNA and nRNA components share the same strand likelihood; they
remain distinguishable only via FL (mRNA after rRNA depletion vs.
nRNA which carries an additional sonication path). In the INTRONIC
bucket, mRNA cannot contribute by definition, leaving a clean
**2-component (gDNA, nRNA) joint (FL, S) mixture** — the simplest
case that v3 needs to solve to retire the §7a limitation.

### 3.3 What this is NOT

- v3 is **not** a strand-aware re-derivation of the per-locus EM. The
  per-locus EM is already strand-aware via `is_same_strand`; nothing
  changes there.
- v3 is **not** a refactor of `StrandModel` itself. SS estimation
  from spliced unique mappers is left intact.
- v3 does **not** add per-locus or per-region SS variation.
  Library-wide SS is the only strand specificity v3 consumes.
- v3 is **not** a chimera-detection refresh. The drop-by-max-FL policy
  already in place stays.

## 4. v3 thrusts and pool semantics

For each pool bucket, list the sources that can physically contribute
and the resulting mixture model:

| Bucket | Possible sources | Joint model |
|---|---|---|
| INTERGENIC | gDNA only | 1-D FL on gDNA-FL only (the cleanest gDNA-FL anchor) |
| INTRONIC | gDNA + nRNA | 2-D (FL, S) mixture; nRNA gets SS, gDNA gets 0.5 |
| EXON_INCOMPATIBLE | gDNA + nRNA + (rare unannotated mRNA) | 2-D (FL, S); treat as gDNA + (m/n)RNA where m and n share strand likelihood |
| EXON_CONTAINED | mRNA + (rare nRNA single-exon overlap) + (rare gDNA) | Held out (mature-mRNA dominated); diagnostic only |
| (chimeras / OOR-FL) | all sources | Held out (Phase 7 OOR policy) |

### 4.1 v3 phases

**Phase 1 — Strand relabeling.** Replace genomic-frame POS/NEG strand
labels with transcript-frame SENSE/ANTISENSE/AMBIG. No model change yet;
purely a data-pipeline relabeling that makes the strand axis
interpretable. `category_counts` shape changes from `(N_CATEGORIES, 4)`
to `(N_CATEGORIES, 3)`. `n_pool_intronic_strand_*` diagnostics rename
correspondingly. Detailed plan: `srd_v3_phase1_strand_relabel.md`.

**Phase 2 — Joint (FL × strand) mixture for INTRONIC.** Extend
`fit_fl_mixture` to a 2-D Bernoulli-FL mixture parameterised by SS
(treated as fixed input from the strand model). Initial scope:
INTRONIC bucket only; gDNA + nRNA two-component case. Output: a
new `pi_intronic_nrna` field on `CalibrationResult`, and an
nRNA-FL model that can later replace the per-locus prior.

**Phase 3 — Generalize to the full pool.** Re-cast the pool mixture
as a per-bucket weighted sum where each bucket contributes only the
sources it can host. Mixing weights `π_g`, `π_m`, `π_n` are shared
across buckets (library-wide); the per-bucket FL distributions
(`p_gDNA-FL`, `p_mRNA-FL`, `p_nRNA-FL`) come out of the joint fit.
Output replaces `π_pool` with a 3-vector of source-fractions.

**Phase 4 — Validate / golden-update / docs.** Same harness as v2
Phase 4; plus a new strand-identifiability diagnostic:

- For each library, plot observed (sense, antisense) fractions in
  INTRONIC vs. SS. A library at SS=0.95 with INTRONIC sense fraction
  ≈ 0.9 indicates ~80% nRNA in INTRONIC, computable as
  `(observed_sense_frac − 0.5) / (SS − 0.5)`. This is the unit-test
  for Phase 1 even before Phase 2 ships.

**Phase 5 — Cleanup.** Remove the v2 §7a "accepted limitation" text;
remove the `pi_pool` shorthand in favor of the source-fraction vector;
bump version.

### 4.2 Acceptance criteria preview

Phase 1 (strand relabeling):
- All existing categorization tests pass with the new enum.
- Phase 4 sweep on VCaP-mixture libraries: INTRONIC sense-fraction
  scales monotonically with library SS (cross-library). For a single
  library, INTRONIC sense-fraction lies between 0.5 (pure gDNA) and SS
  (pure nRNA). If observed > SS or < 0.5, diagnostic surfaces it.
- No quant-level regression vs. Phase 6 baseline (`pi_pool`,
  `gdna_fraction` headline numbers stable to within ±1pp; v3 Phase 1
  is purely a bookkeeping change, not a model change).

Phase 2 (joint mixture):
- On VCaP-mixture libraries, decomposed `π_intronic_nrna` is non-zero
  for unspiked dna00m (where nRNA contamination is the only INTRONIC
  source) and grows monotonically with the nRNA spike axis (when one
  is added).
- gDNA-FL recovered from the joint fit matches v2's gDNA-FL within
  ±5% bin-by-bin on the INTERGENIC-only fit (cross-check).
- Headline `gdna_fraction` on dna80m moves from v2's ~0.53 toward the
  nominal ~0.80 by reattributing the previously nRNA-polluted INTRONIC
  mass.

## 5. Risks & mitigations

| Risk | Mitigation |
|---|---|
| SS too close to 0.5 (unstranded library) → strand axis non-discriminating | Detect at fit time; fall back to v2 1-D mixture with a `quality="weak_strand"` flag. |
| `StrandModel` SS misestimated due to limited spliced fragments | Already surfaced via `n_strand_trained` telemetry; require a minimum count before v3 model engages. |
| AMBIG (both-strand) fragments large in dense regions | Drop from per-strand histograms; report count as diagnostic. |
| 2-D mixture EM convergence harder than 1-D | Start from v2 1-D fit as initialization; add max-iter cap and SQUAREM acceleration as needed. |
| Schema break for downstream consumers | Bump version; document `summary.json` key changes in `CHANGELOG.md`. |

## 6. Out of scope for v3 (preserved from v2)

- Cross-region density signals.
- Mappability-corrected effective length.
- Region partition resurrection.
- Per-block preservation for SPLICE_ARTIFACT.
- StrandModel internals refactor.
- Per-locus / per-region SS estimation.
- New chimera detection heuristics (OOR drop policy stands).

## 7. References

- `srd_v2_phase6_handoff.md` — v2 final state, Phase 7 OOR policy
  closure note.
- `srd_v2_phase2plus_handoff.md` §7a — original nRNA siphon write-up
  that motivates v3.
- `srd_v2_results.md` — v2 baseline numbers; Phase 4 comparison
  reference for v3.
- `srd_v3_phase1_strand_relabel.md` — Phase 1 implementation plan
  (this is what gets implemented first).
