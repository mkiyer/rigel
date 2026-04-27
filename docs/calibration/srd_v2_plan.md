# SRD v2 plan — scanner-integrated calibration & region-machinery removal

**Status**: Proposal, written 2026-04-24 in response to a parallel-Claude
plan that proposed integrating SRD v1 with scanner augmentation.
**Predecessor**: SRD v1 (`docs/calibration/srd_v1_implementation.md`,
`docs/calibration/srd_v1_results.md`) is **shipped and working** —
8/8 VCaP mixture libraries converge with `gdna_fl_quality=good` and
the per-locus γ → α prior (Pass 4) is in place.

This document **does not subsume v1**. It scopes a follow-up that
(a) removes the dead region-partition machinery flagged in v1 Phase 4b,
(b) opportunistically picks up the *one* genuinely good architectural
idea from the parallel plan, and (c) leaves alone the parts of v1 that
are doing real work and have benchmark validation behind them.

---

## 1. Assessment of the incoming v2 proposal

The parallel plan proposes four headline changes:

| # | Proposal | Verdict |
|---|----------|---------|
| 1 | Augment scanner to emit collapsed `exonic_bp` / `intronic_bp` / `intergenic_bp` / `exon_strand_flags` per fragment | **Mixed.** Genuinely cleaner data model, but C++ change + buffer schema migration + ~13 B/fragment memory cost on top of an already-working categorize. Ship as a careful refactor *after* region removal, not bundled with it. |
| 2 | Delete the region partition (`build_region_table`, `RegionAccumulator`, `regions.feather`, `region_cr`, `mappable_effective_length` per region) | **Yes.** Already on v1's Phase 4b carryover backlog. Independent of (1). |
| 3 | Replace categorize logic to read collapsed-buffer fields | **Conditional on (1).** No win without (1) shipped. |
| 4 | Pass 4: per-fragment γ → per-locus α | **Already shipped in v1.** See `src/rigel/locus.py::compute_locus_priors_from_partitions`. |

The strongest claim in the plan — *"once per-fragment metadata is in
place, the region partition disappears"* — is actually **false** in the
current codebase, because the only thing the region machinery still
feeds is `mappable_effective_length` (and stale v5 calibration that is
already gone). Region removal does not require scanner augmentation.

The plan also mis-states the bp-sum invariant: it claims
`exonic_bp + intronic_bp + intergenic_bp == genomic_footprint`, but
`genomic_footprint` is the genomic span (with introns covered by `N`
operations), while the sum of base classifications is bounded by
`read_length` (the sum of M/D/=/X exon-block spans, not query length).
Any C++ implementation must use `read_length` as the conservation
target — see `.github/copilot-instructions.md`.

## 2. What is already done (do not re-do)

From SRD v1, currently in main:

- **Pass 0** — vectorized categorization keyed off per-candidate `exon_bp`
  with the `read_length - max(exon_bp) <= tol` rule
  (`src/rigel/calibration/_categorize.py`). Regression-tested with the
  exact dna80m bug scenarios that motivated v2's proposal #3.
- **Pass 1** — `Pool = EXON_INCOMPATIBLE ∪ INTRONIC ∪ INTERGENIC` with the
  C++ scanner's `frag_length_models.intergenic` accumulator folded in
  (`_simple.py`). The intergenic-FL fix does part of what v2 proposal #1
  was meant to enable.
- **Pass 2** — 1-D mixture EM (`_fl_mixture.py`) with `max_iter=1000` after
  Phase 5 showed the mid-π regime needs >200 iters.
- **Pass 3** — Empirical-Bayes Dirichlet smoothing (`_fl_empirical_bayes.py`).
- **Pass 4** — Per-fragment γ → per-locus α
  (`compute_locus_priors_from_partitions` in `locus.py`).
- Buffer-schema-stable; no C++ changes; 962/962 tests green; goldens
  refreshed.

## 3. What this plan delivers

Two independent, sequenced workstreams. **Each is shippable on its own.**

### Workstream A — Region machinery removal (Phase 4b carryover)

Pure deletion. No new functionality, no behavioral change, no C++ schema
churn beyond removing dead fields.

**Targets:**
- `src/rigel/index.py::build_region_table` and the `regions.feather`
  artifact in the index format.
- `TranscriptIndex.region_df`, `region_cr`, lazy-load helpers.
- `mappable_effective_length` per-region column on the index. (Per-transcript
  effective length stays; it is computed in scoring from the FL model.)
- `src/rigel/native/resolve_context.h::RegionAccumulator` and the
  `fl_table` / `region_id` paths (resolve_context.h ~lines 1225–1430).
- `region_counts` / `region_evidence` plumbing in `pipeline.py`,
  `buffer.py`, `scan.py` (final remnants).
- `tests/test_regions.py`, region portions of `tests/test_mappability.py`,
  `tests/test_calibration_integration.py` (anything still referencing the
  region API).

**Deliverable**: a single PR-sized commit that deletes ~600–800 LOC
(per the v2 plan's estimate). Recompile required. Index format bumps
minor version; older indexes still load (the field is just absent now).

**Risk**: low. The audit step verifies nothing in the live pipeline reads
region tables. SRD v1 already proved this — calibration runs without
touching them.

### Workstream B — Scanner buffer augmentation (optional refactor)

Done **only after** A ships and only if the categorize-Phase-0 cost or
clarity becomes a profiling/maintenance issue. Today categorize is <100 ms
per million fragments and is tested against the production failure modes;
there is no acute pressure.

**If/when we do it, the right scope is:**

1. Add **only** the strand-presence bit-pair `exon_strand_flags` (1 byte)
   and a single per-fragment `read_length_classified_exonic_bp` (4 bytes).
   Skip the redundant `intronic_bp` / `intergenic_bp` triple — the
   invariant `read_length = exonic + intronic + intergenic` makes one of
   the three derivable, and categorize uses only `exonic_bp` plus the
   strand flags. Net cost: **5 bytes/fragment** (vs the v2 plan's 13).
2. Compute the scalar `exonic_bp_collapsed` in the scanner inner loop
   from the **union of exon overlaps across all candidate transcripts**.
   This is what `np.maximum.reduceat(exon_bp, ...)` is doing in Python
   today — the scalar value is identical because every byte covered by
   at least one candidate's exon counts once. Move the reduction into
   C++.
3. Replace the per-candidate `intron_bp` *semantic footgun* by renaming
   it to `t_within_span_nonexon_bp` (it is per-transcript and that's
   what it actually measures). Keep it for the per-locus EM, which uses
   `intron_bp` via the locus partition; do not delete.
4. New invariant test: `exonic_bp + (read_length - exonic_bp) == read_length`
   trivially; the real test is that `exonic_bp ≤ read_length` and that
   the scalar matches `np.maximum.reduceat(exon_bp_per_candidate, ...)`
   on a chunk of representative real data.

**What we do NOT do**:
- We do not introduce `intronic_bp` and `intergenic_bp` as separate
  collapsed buffer fields. They are not used anywhere downstream and
  add 8 B/fragment for no benefit.
- We do not delete the per-candidate `exon_bp` / `intron_bp` arrays.
  The per-locus EM uses them for transcript-level scoring.
- We do not rewrite the categorize module wholesale. The existing
  `np.maximum.reduceat` call becomes `chunk.exonic_bp_collapsed[mask]` —
  a one-line change. Tests stay; semantics are identical because the
  scalar value is identical.

### Workstream C — Documentation and CHANGELOG (the original Phase 6)

- `CLAUDE.md`, `README.md`, `docs/MANUAL.md`, `docs/METHODS.md` updated
  to remove regional-density language, document SRD v1, document the new
  index format (no `regions.feather`).
- `CHANGELOG.md` entry covering both A and (if shipped) B.
- Bump version (minor: index format change, no breaking API change for
  end users).

## 4. Sequencing & exit criteria

| Step | Action | Exit criterion |
|------|--------|----------------|
| C.1  | Update `CHANGELOG.md`, `CLAUDE.md`, `README.md`, `docs/MANUAL.md`, `docs/METHODS.md` for SRD v1 as currently shipped. Bump version. | Docs match code; one commit. |
| A.1  | Audit: confirm zero live-path consumers of `region_df` / `region_cr` / `mappable_effective_length` (per-region). | Audit doc with grep evidence. |
| A.2  | Delete Python region surface (`index.py`, `pipeline.py`, `mappability.py`). | `pytest tests/ -v` green; goldens unchanged. |
| A.3  | Delete C++ region surface (`resolve_context.h::RegionAccumulator`, `fl_table`, `region_id` paths in `bam_scanner.cpp`). Recompile. | `pytest tests/ -v` green; one VCaP smoke run completes with identical numbers. |
| A.4  | Delete dead tests (`test_regions.py`, region parts of `test_mappability.py`, `test_calibration_integration.py`). | LOC drop measured; all tests green. |
| —    | **Decision gate**: re-run VCaP mixture; verify no regression. Optionally proceed to B. | Numbers match `srd_v1_results.md` ±0.5pp. |
| B.1  | (Optional) Scanner emits scalar `exonic_bp_collapsed` and `exon_strand_flags`. | New invariant test green; categorize switches to scalar field; behavior identical. |
| B.2  | (Optional) Rename `intron_bp` → `t_within_span_nonexon_bp` everywhere. | Tests green; docstring updated. |

A and C are the actual deliverables. B is documented as available work
for a future cycle.

## 5. What we are explicitly *not* doing (and why)

### 5.1 Not bundling region removal with scanner augmentation

The parallel plan's central architectural pitch is "scanner augmentation
**enables** region removal." That's not true today: region removal is
already enabled (v1 calibration doesn't use regions). Bundling them
forces a high-risk C++ buffer-schema change into a low-risk deletion PR.

### 5.2 Not re-implementing Pass 0 from collapsed fields

The parallel plan's Pass 0 spec walks the same logic as the existing
`categorize_chunk`, just with collapsed inputs:

```python
# Parallel-plan formulation
intronic_mask = (exonic_bp == 0) & (intronic_bp > 0)
incompat_mask = (exonic_bp > 0) & ((intronic_bp > tol) | (intergenic_bp > tol))
```

vs current

```python
# Already shipped, semantically equivalent
intronic_mask = has_cands_unspliced & (max_exon == 0)
incompat_mask = has_cands_unspliced & (max_exon > 0) & ((rl - max_exon) > tol)
```

The only difference is whether `max(exon_bp_per_candidate)` is computed
in C++ (parallel plan) or in numpy (current). Workstream B captures
this if we decide it matters; otherwise it is a tied boxing match.

### 5.3 Not changing per-locus γ aggregation

Pass 4 (`compute_locus_priors_from_partitions`) is shipped, has scenario-
test coverage, and matches the parallel plan's design exactly. No work.

### 5.4 Not adding the strand-flag bitfield to address antisense overlap

The current `categorize_chunk` already handles the antisense-overlap
case via `chunk.ambig_strand` (set by the C++ scanner) plus the
per-candidate `t_to_strand[t_indices]` consensus. The bit-pair
formulation is mildly cleaner but observably equivalent. It is part of
Workstream B if we ship B at all.

## 6. Open questions

1. **Index format compatibility.** Removing `regions.feather` is a
   breaking change for users with existing indexes. A clean path:
   bump index format minor version; load returns `region_df = None`
   for old indexes (tolerated); writers omit it. No re-indexing forced.
2. **Effective-length computation.** Currently `mappable_effective_length`
   on the *region* table is being deleted. Per-*transcript* effective
   length is already computed by `frag_length_models.rna_model.compute_all_transcript_eff_lens`
   (`pipeline.py:325`) and is unaffected.
3. **Workstream B trigger.** What benchmark or profiling evidence would
   make B worth doing? Suggested triggers: (a) categorize cost > 5% of
   wall time on a 100M-fragment library, (b) a third bug like the
   `intron_bp` semantic footgun shows up.

## 7. Honest impressions

The parallel plan is *technically sound* and the bit-pair strand flags
are genuinely elegant. But it's optimising for the wrong axis: it
proposes a high-touch rewrite (C++ scanner change, buffer migration,
schema bump, full pipeline rewire) for benefits that are
**already captured** by the shipped v1 (the bug fix, the intergenic-FL
fold-in, Pass 4, the convergence cap bump).

The one piece of v2 that still has uncaptured value — collapsing the
per-candidate `exon_bp` reduction into the scanner — is small and easy
to ship later as a single targeted refactor *after* the dead region
machinery is gone. Doing it now bundles two unrelated risks.

**Recommended action**: ship Workstreams C (docs) and A (region removal).
File Workstream B as a follow-up, gated on profiling evidence. Do not
re-do Pass 0–4; v1 already implements what v2 wanted.

---

## 8. Future optimizations and algorithmic improvements

A roadmap of opportunities ordered by **expected payoff per unit of
risk**. Each item lists what it costs, what it buys, and what the
gating evidence should be before we commit. Items below the line are
research/scoping work, not engineering tasks.

### 8.1 SRD calibration — quick wins (Workstream B+ extensions)

#### 8.1.1 Drop per-bin smoothing in the FL mixture EM

`src/rigel/calibration/_fl_mixture.py` currently adds `smoothing=1.0` to
every gDNA-FL bin per iteration (~1000 bins). On a `n_pool ≈ 30k` library
(dna00m) that is ~3% of the total mass injected as pseudo-counts each
step. SRD v1 results §5.3 attributes the high-π underestimate
(`dna40m`/`dna80m`: π plateaus at 0.68 vs nominal ~0.85) to this floor.

**Cost**: 5 LOC change + sweep test on the 8-library VCaP ladder.
**Expected payoff**: +5–10pp recovered gDNA fraction in the 60–80%
contamination regime; π_pool tracks nominal across the full ladder.
**Risk**: dna00m corner becomes noisier (need to verify the pure-RNA
library still reports < 5% gDNA).
**Gating evidence**: re-run `analyze_srd_vcap_mixture.py` with smoothing
in {1.0, 0.1, 0.01, 0.001}; pick the value that minimises max-error
across the ladder.

#### 8.1.2 Switch convergence rule to log-likelihood per observation

Today's `|Δπ| < 1e-4` rule needs `max_iter=1000` to handle the long
shallow ridge in the mid-π regime (Phase 5 dna10m converged at iter 235).
A per-observation logL rule (`|ΔlogL| / n_pool < 1e-7`) is well-scaled
and would converge in tens of iterations on big libraries.

**Cost**: 10 LOC + tighten the unit tests in `test_fl_mixture.py`.
**Expected payoff**: 5–10× speedup of Pass 2 on large libraries; remove
the `max_iter` knob entirely.
**Risk**: low — logL is monotone-non-decreasing in EM, so this is
strictly more principled.
**Gating evidence**: profile Pass 2 on dna80m before and after.

#### 8.1.3 Replace 1-D mixture with a 2-component generalized linear FL model

The current `pool ∝ π·gDNA(L) + (1-π)·RNA(L)` treats RNA-FL as fixed.
On gDNA-rich libraries (dna40m+) the *pool's* RNA component is dominated
by intronic mRNA and short-fragment artefacts, NOT spliced fragments —
so anchoring to `spliced_FL` underestimates the real RNA contribution
to the pool and inflates π. A model

```
pool(L) ∝ π · gDNA(L) + ρ · RNA_short(L) + (1-π-ρ) · RNA_spliced(L)
```

with `RNA_short(L)` learned from the unspliced-exonic-sense category
would fit the pool tail more faithfully.

**Cost**: ~80 LOC; new EM loop; a new diagnostic field.
**Expected payoff**: tightens the high-π underestimate from ~10pp to
~3pp on synthetic data; resolves SRD v1's main residual systematic.
**Risk**: identifiability of `RNA_short` is non-trivial without
strand-balanced controls; may need an SS-conditional gate.
**Gating evidence**: synthetic sweep (in `scripts/synthetic_sim_sweep.py`)
showing two-component RNA recovers π within ±2pp at gdna_fraction ∈ {0.7, 0.8}.

#### 8.1.4 Stranded library bonus: per-strand pool decomposition

For strongly-stranded libraries (`strand_specificity > 0.95`) the
antisense-exonic category becomes a near-pure gDNA channel — currently
diagnostic only, never folded into the pool. Adding it is a "free" 30%
boost to pool size for the libraries that matter most.

**Cost**: 15 LOC in `_simple.py` Pass 1; a new SS gate (single threshold,
documented).
**Expected payoff**: tighter `pi_pool` estimates on stranded libraries;
particularly helpful at low contamination where pool is smallest.
**Risk**: at imperfect SS, antisense pool is RNA-contaminated; this is
exactly why the parallel plan kept it out. The gate must be sharp.
**Gating evidence**: synthetic sweep at SS ∈ {0.7, 0.85, 0.95, 0.99};
include a regression test for the unstranded case.

### 8.2 Per-locus EM — identifiability fixes

#### 8.2.1 Address the nRNA siphon at high gDNA contamination

SRD v1 results §5.4 documents nRNA fraction growing 2.6% → 10.2% across
the gDNA ladder. The root cause: intronic gDNA fragments and intronic
nRNA fragments have the same coordinate signature; only the strand
model can distinguish them at finite SS. Three mitigations, in order of
ambition:

1. **Strand-conditional nRNA prior**. Today the per-locus α prior treats
   nRNA and gDNA symmetrically. Adding `α_nRNA[ℓ] ∝ (1 − π_global) ·
   antisense_share[ℓ]` would suppress nRNA when antisense fragments are
   present (a gDNA signature).
2. **Per-locus γ broken into γ_gdna + γ_nrna** using a 3-component
   mixture in the per-fragment posterior. Requires a per-fragment
   nRNA-FL likelihood, which we currently take as `RNA_FL` (assumption
   noted in v1 docs).
3. **Locus-level coverage profile prior**. Real nRNA has a 5'→3' bias;
   gDNA is uniform. A simple coverage-profile statistic per locus
   (already computable from `chunk.genomic_start`) would discriminate.

**Cost**: (1) ~30 LOC, (2) ~150 LOC, (3) ~250 LOC including new C++ path.
**Expected payoff**: (1) ~30% reduction in nRNA siphon at dna80m;
(2) ~70% reduction; (3) near-elimination if the bias signal is real.
**Risk**: increasing — (1) is principled, (3) requires empirical
validation that VCaP's nRNA actually has the assumed coverage shape.

#### 8.2.2 Locus partitioning for mega-loci

The per-locus EM currently treats every connected component of the
multimapper graph as one locus. On real human data with high-NH reads
(repeat regions, gene families) this produces mega-loci with thousands
of transcripts and millions of fragments. The EM solver is robust
(`.github/copilot-instructions.md`) but slow.

**Idea**: spectral or modularity-based subdivision of the multimapper
graph; run EM on each block; reconcile boundary fragments via a small
outer EM. Standard technique in metagenomic binning.

**Cost**: ~400 LOC; needs new C++ kernel for the spectral step.
**Expected payoff**: 3–10× speedup on libraries with lots of repeats;
better convergence in the per-block EMs.
**Risk**: subdivision can break true biological multimappers;
correctness gate is preserving the equivalence-class likelihood.
**Gating evidence**: profile per-locus EM time on the largest VCaP
locus; if < 5% of total wall, defer.

### 8.3 Scanner & pipeline architecture

#### 8.3.1 Workstream B itself (scanner emits collapsed `exonic_bp`)

Already scoped above. Gating evidence: categorize > 5% wall time.
Currently it is far below.

#### 8.3.2 Memory-mapped FragmentBuffer

For multi-100M-fragment libraries the buffer can spill to disk
(`spill_dir` in BamScanConfig). The current spill path uses pickle;
a memory-mapped feather format would let downstream EM read columns
without copying.

**Cost**: 200 LOC across `buffer.py` and `scan.py`.
**Expected payoff**: 30–50% reduction in peak RSS on 200M+ fragment
libraries.
**Risk**: low; format change but fully internal.

#### 8.3.3 Streaming SRD calibration

Currently SRD reads the entire FragmentBuffer in chunks but accumulates
into per-library histograms before Pass 2 EM. Pass 2 could be invoked
periodically with monotonically-growing pool histograms; the EM is so
cheap that streaming convergence checks are essentially free.

**Cost**: 50 LOC.
**Expected payoff**: enables early-termination ("calibrated, stop
scanning") for QC mode.
**Risk**: trivial.

### 8.4 Scoring & EM kernel

#### 8.4.1 Mappability-corrected effective length

The `mappable_effective_length` per-transcript field has been on the
roadmap forever but never wired into scoring. Today scoring uses
`compute_all_transcript_eff_lens` which is a pure FL convolution against
transcript length; it ignores mappability. For repeat-rich transcripts
this overstates expected fragment count.

**Cost**: ~150 LOC + a new index artefact (per-transcript mappability
profile).
**Expected payoff**: corrects systematic over-quantification of repeat
transcripts; resolves a class of known false-positive isoforms.
**Risk**: needs validation that the mappability track (uniqueome /
crg36) actually matches the aligner being used.

#### 8.4.2 Replace VBEM with online stochastic VBEM

VBEM iterates the full equivalence-class graph each step. For large
loci this is the dominant cost. Online VBEM (Hoffman et al. 2010)
processes mini-batches of equivalence classes per step; converges in
~5× fewer wall-clock seconds with the same final ELBO.

**Cost**: ~200 LOC C++ change in `_em_impl`.
**Expected payoff**: 3–5× speedup on mega-loci; same final estimates.
**Risk**: well-trodden territory; salmon and kallisto both do variants
of this.

#### 8.4.3 Per-fragment SQUAREM acceleration in mixture EM

The 1-D mixture EM in Pass 2 is a textbook target for SQUAREM; the
existing per-locus EM already uses it (`_em_impl`). One-line change.

**Cost**: 20 LOC.
**Expected payoff**: convergence in ~30 iters instead of ~250 on
mid-π libraries; obviates the need for the `max_iter=1000` cap.
**Risk**: trivial.

### 8.5 Observability & QC

#### 8.5.1 Per-fragment classification audit trail

Today SRD emits aggregate `category_counts`. A per-locus drill-down
(N_EXONIC / N_INCOMPATIBLE / N_INTRONIC / N_INTERGENIC / N_SPLICED)
would let users debug surprising gDNA estimates at specific loci.

**Cost**: 80 LOC; new `loci_quality.feather` artefact.
**Expected payoff**: massively improves debuggability; enables loci-level
dashboards.
**Risk**: minor — adds an artefact to the output directory.

#### 8.5.2 Calibration sanity checks at end of run

Surface a one-line summary in `summary.json`:
- π_pool vs locus-aggregated γ̄ (sanity: should be within 5pp).
- Pool size vs n_unique_unspliced_fragments (should be 5–30%).
- gDNA-FL mean vs RNA-FL mean (sanity: should be smaller for gDNA).

**Cost**: 30 LOC; expand `CalibrationResult.to_summary_dict`.
**Expected payoff**: catch silent calibration failures before they
reach the analyst.
**Risk**: trivial.

### 8.6 Out-of-scope research questions

These are *research* problems, not engineering tasks. Listed for
completeness; each requires a literature review and at least one
prototype before scoping.

- **Joint multi-sample calibration** for batch-effect correction in
  high-throughput sequencing centres.
- **Bayesian model selection** between SRD v1's 2-component pool
  mixture and a non-parametric mixture (DPM).
- **Allele-specific quantification** under SRD's gDNA-aware
  framework (currently allele assignment is post-hoc).
- **Single-cell SRD**: per-cell gDNA contamination estimation in droplet
  RNA-seq, where pool size per cell is tiny and EB across cells is
  essential.
- **Long-read SRD**: the `read_length - max(exon_bp) <= tol` rule
  generalises trivially to long reads; the mixture model probably
  needs reparameterisation against length-dependent error rates.

### 8.7 Suggested next-cycle priorities (PI's view)

If forced to pick three for the next cycle:

1. **§8.1.1** (drop per-bin smoothing) — fixes the headline residual
   from SRD v1 results §5.3, ~5 LOC, no risk.
2. **§8.4.3** (SQUAREM in mixture EM) — eliminates the convergence
   issues in §8.1.2 with less work; one line.
3. **§8.5.1 + §8.5.2** (observability) — pays for itself the first time
   a benchmark surprises us.

Everything else can wait until §8.7.1 raises a new question.
