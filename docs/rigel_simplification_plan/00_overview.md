# Rigel EM Simplification Plan

## Motivation

The 1,024-run complex locus ablation study
([complex_locus_ablation_studies_v1.md](../complex_locus_ablation_studies_v1.md))
revealed that Rigel's EM initialization and prior framework is over-engineered.
The system has 27 EMConfig parameters, 5 distinct prior mechanisms, and
multiple layers of empirical Bayes shrinkage that actively hurt accuracy —
especially for moderate nRNA signals and low strand specificity.

This plan breaks the simplification into phases, ordered by **expected accuracy
impact** and **implementation simplicity**. Each phase follows the same
workflow:

1. Publish a detailed implementation plan
2. Implement the changes
3. Re-run the 1,024-run complex locus ablation
4. Compare results against the baseline and prior phases

The baseline is the current system (pre-simplification ablation results in
`complex_locus_ablation_studies_v1.md`).

---

## Phases

### Phase 1: Decouple nRNA from mRNA — Eliminate the η Parameterization

**Root causes addressed**: RC1 (nRNA killed at SS ≤ 0.6), RC2 (3-tier nRNA EB
over-shrinks moderate nRNA), RC3 (nRNA exonic mass leaks to mRNA because η
can't reach 1.0)

**Core insight**: The nRNA fraction (η) parameterization that links each nRNA
component to its "child" mRNA transcripts via a Beta prior is unnecessary
complexity. nRNA spans already compete with their overlapping mRNA transcripts
for the same fragments through the EM likelihood. The η constraint tries to
reduce degrees of freedom but:

- **Biologically unjustified**: η varies enormously across genes (from ~0 for
  stable mRNAs with slow transcription to ~1 for rapidly transcribed,
  unstable transcripts). A shared prior across genes is fundamentally wrong.
- **Actively harmful**: The Beta MAP prevents η from reaching 0 or 1, creating
  false positives when only nRNA is present (+864% mRNA error) and suppressing
  nRNA when moderate.
- **Over-parameterized**: 10 EMConfig fields, 3-tier hierarchical shrinkage
  with MoM κ estimation — all to bias a parameter the EM can estimate on its
  own.

**What changes**: Remove the hierarchical M-step from the C++ EM solver. Remove
the 3-tier nRNA EB prior computation. Fix nRNA component gating so nRNA is not
killed at low SS. nRNA becomes a regular independent EM component.

**Parameters removed**: 9 (all `nrna_frac_*` fields)

**Difficulty**: Moderate — touches C++ EM solver, Python prior code, config,
CLI, and tests. But the changes are predominantly *deletion*.

**Detailed plan**: [phase1_decouple_nrna.md](phase1_decouple_nrna.md)

---

### Phase 2: Replace OVR Prior with Flat Prior + Informed Warm Start

**Root causes addressed**: RC3 (remaining — OVR seeds mRNA with gDNA mass),
RC5 (single-exon gene false positives)

**Core insight**: The OVR coverage-weighted prior distributes ambiguous
fragment mass to all eligible components proportional to coverage. This seeds
unexpressed genes with nonzero mass, creating false positives — especially for
single-exon genes that can't be disambiguated from gDNA by splice junctions.

**What changes**: Remove the `prior_gamma` OVR term. Use flat α only. Redesign
warm start to use unambiguous counts (splice-informed) rather than total
coverage.

**Parameters removed**: 1 (`prior_gamma`)

**Difficulty**: Moderate — localized to `compute_ovr_prior_and_warm_start()` in
em_solver.cpp.

---

### Phase 3: Scale Strand Symmetry Penalty by Strand Specificity

**Root causes addressed**: RC1 (synergy with Phase 1), RC6 (gDNA mass
misallocation)

**Core insight**: The strand symmetry penalty assumes balanced strand coverage
implies gDNA. At SS=0.5, *all* RNA coverage is balanced — the penalty wrongly
boosts gDNA. The penalty should scale with available strand information.

**What changes**: $\kappa_{eff} = \kappa \cdot (2 \cdot SS - 1)^2$. Fold
`strand_symmetry_pseudo` into the κ scaling.

**Parameters removed**: 1 (`strand_symmetry_pseudo`)

**Difficulty**: Low — single formula change in the C++ M-step.

---

### Phase 4: Simplify gDNA Initialization

**Root causes addressed**: RC6 (gDNA EB over-regularizes)

**Core insight**: The 3-tier gDNA EB (global → reference → locus) adds
complexity for marginal benefit. A direct per-locus strand-based estimate
floored at the intergenic density would be simpler and possibly more accurate.

**What changes**: Replace `compute_eb_gdna_priors()` with a direct per-locus
estimate. Remove MoM κ estimation for reference/locus levels.

**Parameters removed**: 4 (`gdna_kappa_ref`, `gdna_kappa_locus`,
`gdna_mom_min_evidence_ref`, `gdna_mom_min_evidence_locus`)

**Difficulty**: Moderate — touches locus.py gDNA initialization and the EB
machinery.

---

### Phase 5: Revisit nRNA-mRNA Exonic Ambiguity

**Root causes addressed**: RC3 (fundamental identifiability)

**Core insight**: After Phases 1-4, the exonic ambiguity between nRNA and mRNA
may be adequately handled by the EM's natural competition. If not, explore
splice-junction-informed disambiguation: spliced reads are evidence for mRNA,
unspliced exonic reads are ambiguous between mRNA and nRNA.

**Depends on**: Phase 1-4 results — this phase may be unnecessary if the EM
handles the competition well once the η constraint is removed.

**Parameters removed**: 0

**Difficulty**: High — may require changes to how fragment evidence is
attributed.

---

### Phase 6: Convergent Gene nRNA Handling

**Root causes addressed**: RC4 (cross-strand nRNA misattribution)

**Core insight**: When nRNAs from different strands overlap, intronic reads
from one strand's nRNA can be misattributed to the other strand's mRNA.

**Depends on**: Phases 1-4 — the severity may change once priors are
simplified.

**Parameters removed**: 0

**Difficulty**: High — may require strand-aware nRNA component indexing.

---

## Parameter Reduction Summary

| Phase | Params Removed | Running Total |
|-------|---------------|---------------|
| Baseline | — | 27 |
| Phase 1 (decouple nRNA) | 9 | 18 |
| Phase 2 (OVR prior) | 1 | 17 |
| Phase 3 (strand penalty) | 1 | 16 |
| Phase 4 (gDNA EB) | 4 | 12 |
| Phase 5 (exonic ambiguity) | 0 | 12 |
| Phase 6 (convergent genes) | 0 | 12 |
| **Final** | **15** | **12** |

---

## Validation Protocol

After each phase:

1. Re-run the 1,024-run complex locus ablation sweep
   (`scripts/complex_locus_config.yaml`)
2. Run the analysis script (`scripts/analyze_complex_locus.py`)
3. Compare key metrics against baseline:
   - Total mRNA accuracy by nRNA level × SS × gDNA fraction
   - nRNA accuracy (signed error)
   - gDNA accuracy (signed error)
   - False positive counts and magnitudes
   - Pattern 15 stress test (nRNA-only)
   - Convergent gene accuracy
4. Document results in a phase-specific ablation report
5. Run the full test suite (`pytest tests/`) to ensure no regressions
