# SRD v1 Implementation Plan (Consolidated)

**Status:** APPROVED — phases below are the contract for execution.
**Supersedes:** `simplification_overhaul_v1.md`, `simplification_plan_other_claude.md`,
`srd_v1_plan.md`. Phase 0 audit findings live in `srd_v1_audit.md`.

This document merges the user's plan and the agent's revised plan. Substance
is identical between the two; this version locks in module names, the
continuous cluster-γ formula, the dual-channel Pool A, and the
`partition.py` → `locus_partition.py` rename, and aligns with the Phase 0
audit findings.

---

## 1 — Context and goals

Rigel's calibration step (`src/rigel/calibration/`) currently runs a v5
three-channel mixture EM (density + strand + FL) over annotation regions.
It fails on hybrid-capture data because the density channel assumes a
single global gDNA rate `λ_G` while real capture libraries exhibit 100–140×
spatial variation. The fix is not more parameters; it is **changing what we
classify**: stop classifying *regions* with density, classify individual
*fragments* with cheap geometric and strand signals.

Calibration's only load-bearing output to downstream is the gDNA
fragment-length model, `gdna_fl_model`, consumed by the C++ scorer. SRD
preserves that contract and discards everything else.

### Design goals

- **Library-agnostic.** Stranded and unstranded both supported. **No
  magic-number SS threshold** anywhere.
- **Blazing fast.** Closed-form per-fragment classification; one well-posed
  1-D mixture fit. No region-level EM, no quadrature, no dispersion.
- **No cross-region density.** Immune to hybrid-capture spatial bias.
- **First-principles pools.** Exon fit, splice compatibility, and
  strand-balance within overlap clusters. No annotation heuristics.
- **Graceful degradation.** Empirical-Bayes shrinkage to a global FL
  handles QC-fail libraries, low SS, and small pools without crashing or
  branching.
- **Delete, don't deprecate.** Pre-1.0 single-repo project. Dead code goes
  to `_deprecated/`; no zombie compat fields on `CalibrationResult`.

### Headline Phase 0 findings (locked)

1. **C++ requires zero source changes.** The C++ scorer reads only
   `gdna_fl_model`. The C++ EM solver reads only per-locus `α_gdna`/`α_rna`
   arrays — agnostic to how they are computed. With `compute_locus_priors`
   deleted and zeros passed for `α`, the EM warm-start is uninformative.
2. **`partition.py` stays but is renamed** to `locus_partition.py`. It is
   the per-locus CSR scatter for the per-locus EM, *not* region
   partitioning. Tests in `tests/test_partition.py` keep their assertions;
   the file is renamed to `tests/test_locus_partition.py`.
3. **`compute_locus_priors` is the sole downstream consumer** of
   `region_e_gdna`/`region_n_total`. Deleting it cleanly severs the entire
   region data flow.
4. **`FragmentLengthModels.global_model` already exists** and
   `FragmentLengthModel.finalize(prior_counts=..., prior_ess=...)` already
   implements EB shrinkage. `_fl_empirical_bayes.py` is a thin wrapper.
5. **`region_id` is NOT in `FragmentBuffer`.** Region IDs only live in the
   optional `fl_table` returned by the scanner. Removing region machinery
   has zero buffer-schema impact.

---

## 2 — SRD design

### 2.1 Signals

| Signal | Strength | Library | Used in |
|---|---|---|---|
| Spliced N-junction CIGAR | Hard — 100% RNA | All | SPLICED category |
| Exon-fit geometry (with bp tolerance) | Strong — excludes mature mRNA | All | EXON_INCOMPATIBLE category |
| Antisense vs tx_strand | Strong at high SS, weaker at low | Stranded only | Pool A (antisense channel) |
| Strand-balanced overlap cluster | Strong — per-cluster γ from imbalance | Stranded only | Pool A (cluster channel) |
| Fragment length (after training) | Soft — posterior weight | All | Pool B (1-D mixture) |

**Explicitly dropped:** cross-region density, Poisson/NegBin rate
estimation, LogNormal expression priors, Beta-Binomial dispersion, region
Dirichlet warm-starts, capture-class composite Poisson.

### 2.2 Pass 0 — Fragment categorization (single buffer walk)

Walk `FragmentBuffer` once. Classify each fragment using geometric tests
against the annotation index. **Reuse the C++ scanner's existing
per-fragment fields** (`splice_type`, `t_indices`, `ambig_strand`,
`exon_strand`, `sj_strand`); the scanner already determines tx_strand
ambiguity (mixed-strand exon overlap → AMBIG). The Python categorizer adds
only the **exon-fit geometry test** with overhang tolerance.

Categories (uint8 enum):

1. **`SPLICED`** — any N-CIGAR. Read directly from `splice_type`.
2. **`UNSPLICED_SENSE_EXONIC`** — fits inside one annotated exon
   (within `exon_fit_tolerance_bp`); strand matches a unique tx_strand.
3. **`UNSPLICED_ANTISENSE_EXONIC`** — fits inside one annotated exon;
   strand opposite to a unique tx_strand.
4. **`UNSPLICED_EXONIC_AMBIG`** — fits inside an exon, but tx_strand is
   undefined: unstranded library, mixed-strand exon overlap, or unknown.
5. **`EXON_INCOMPATIBLE`** — unspliced; does not fit any annotated exon
   within tolerance (overhangs an exon-intron boundary, or longer than any
   containing exon). Cannot be mature mRNA.
6. **`INTRONIC`** — unspliced; no annotated exon overlap; ≥1 transcript
   span overlap.
7. **`INTERGENIC`** — unspliced; no transcript overlap at all.

Equivalences for FL training: `EXON_INCOMPATIBLE ∪ INTRONIC ∪ INTERGENIC`
form the "non-mRNA pool" (gDNA + nRNA mix) used by Pool B.

Implementation notes:
- Numpy first, exhaustive synthetic tests, profile-driven C++ port deferred.
- `exon_fit_tolerance_bp` config knob, default 5 bp, applied symmetrically
  to both exon ends.
- All strand assignments use the existing R2-flip convention from the BAM
  scanner. No new strand semantics.

### 2.3 Pass 1 — Overlap-cluster identification (unspliced fragments only)

Form connected components of overlapping unspliced fragments (any category
except `SPLICED`):
- Sort unspliced fragments by `g0`.
- Single linear sweep: fragment $j$ joins the current cluster iff
  $j.g_0 \le \max_{k \in \text{cluster}} k.g_1$, else start a new cluster.
- Each fragment receives a `cluster_id`. O(n log n) with sort, O(n) sweep.

Per cluster, record `N_total`, `S` (sense, stranded only), `AS`
(antisense, stranded only), and contained fragment indices.

This is **not** a region partition. Clusters are fragment-data-driven, not
annotation-driven, and produce no persistent index.

Optional safeguard against deep-locus megamergers: cap cluster span (e.g.
≤ 10 kb). Defer until profiling proves it matters.

### 2.4 Pass 2A — Strand-based gDNA pool (Fit A; stranded libraries)

For each cluster with `N ≥ min_cluster_N` (default 5) and a defined
tx_strand:

- Observed sense fraction: $\text{sf} = S / (S + AS)$.
- Model: pure-gDNA (γ=1) gives $\text{sf}=0.5$; pure-RNA (γ=0) gives
  $\text{sf}=\text{SS}$. Linear interpolation:
  $$\text{sf} = \gamma \cdot 0.5 + (1-\gamma) \cdot \text{SS}$$
- Closed-form estimator:
  $$\gamma_{\text{cluster}} = \frac{\text{SS} - \text{sf}}{\text{SS} - 0.5}$$
  clipped to $[0, 1]$. **Defined only when SS > 0.5**; for unstranded
  libraries Pool A's cluster channel is empty by construction.
- Per-cluster confidence weight from binomial variance:
  $\mathrm{Var}(\text{sf}) = \text{sf}(1-\text{sf})/N$.

**Pool A is the union of two channels:**

- **Cluster channel** — fragments from clusters with $\gamma_{\text{cluster}} \ge \gamma_{\min}$
  (default 0.8). Per-fragment weight = $\gamma_{\text{cluster}}$.
- **Antisense channel** — all `UNSPLICED_ANTISENSE_EXONIC` fragments
  with confident tx_strand. Per-fragment weight = $\max(0, 2 \cdot \text{SS} - 1)$.
  Catches gDNA in clusters too small for the cluster channel.

If a fragment qualifies for both channels, take the max of the two weights.

**Fit A produces `gDNA_FL_strand`** as a weighted histogram over Pool A
fragment lengths, smoothed via the EB prior (Pass 3).

When SS ≈ 0.5, Pool A's cluster channel is undefined and the antisense
channel collapses; Fit A produces nothing and Fit B carries everything.

### 2.5 Pass 2B — Exon-incompatibility 1-D mixture (Fit B; all libraries)

Pool B = `EXON_INCOMPATIBLE ∪ INTRONIC ∪ INTERGENIC` fragments. This pool
is gDNA + nRNA.

Under the assumption $\text{nRNA\_FL} \approx \text{RNA\_FL}$ (both reflect
post-library-prep fragmentation of long precursor molecules), fit the 1-D
mixture with `RNA_FL` fixed (from SPLICED fragments):

$$
\text{pool\_FL}(l) = \pi \cdot \text{gDNA\_FL}(l) + (1-\pi) \cdot \text{RNA\_FL}(l)
$$

- `gDNA_FL(l)` is free as a Dirichlet-smoothed histogram.
- `π` is the pool-level gDNA fraction (scalar).
- EM converges in <30 iterations; well-posed because `RNA_FL` is fixed.

**Fit B produces `gDNA_FL_geom`** and `π_pool_b`.

If Pool B is too small (default `< min_pool_b = 300`), Fit B is skipped.

### 2.6 Pass 3 — Empirical-Bayes fusion and global-FL fallback

**Global FL prior.** Histogram over *all* fragments (unfiltered). Used as
the Dirichlet prior for all FL models with effective sample size
`fl_prior_ess` (default 500). Reuses the existing
`FragmentLengthModel.finalize(prior_counts=..., prior_ess=...)` primitive.

**`RNA_FL` training.** Histogram SPLICED fragment lengths, smoothed by the
global-FL prior. Near-zero SPLICED count (QC-fail) → collapses to global FL.

**`gDNA_FL` fusion.**

- Only Fit A converged → `gdna_fl_model = gDNA_FL_strand`.
- Only Fit B converged → `gdna_fl_model = gDNA_FL_geom`.
- Both converged → compute KL divergence between them.
  - If $\mathrm{KL} < \text{kl\_tolerance}$ (default 0.1): fuse via weighted
    average using each fit's effective pool size as weight.
  - Else: log warning; prefer Fit A (strand-based purity is
    higher-confidence when it works).
- Neither converged → fall back to the global FL.

**`gdna_fl_model_quality` enum:**
- `"good"` — pool ≥ minimum, converged, both fits agree if both ran.
- `"weak"` — converged but small pool, or only one fit ran.
- `"fallback"` — no reliable signal; using global FL.

**Graceful degradation:**
- **gDNA ≈ 0**: Pool A empty, Fit B returns $\pi \to 0$. `gdna_fl_model`
  falls back to global FL with `quality = "fallback"`. Per-locus FL LLR ≈ 0
  everywhere → EM relies on splice + strand alone. Correct behavior.
- **RNA ≈ 0** (QC-fail): SPLICED ≈ 0 → `RNA_FL` collapses to global → so
  does `gDNA_FL`. They become indistinguishable. Emit loud warning;
  `quality = "fallback"`. Tool completes; downstream EM should be treated
  skeptically.

### 2.7 New `CalibrationResult` schema

Fields kept (or newly added):

- `gdna_fl_model: FragmentLengthModel`
- `rna_fl_model: FragmentLengthModel` — hoisted for provenance
- `global_fl_model: FragmentLengthModel` — the EB anchor
- `gdna_fl_model_quality: Literal["good", "weak", "fallback"]`
- `strand_specificity: float`
- Category counts: `n_spliced`, `n_sense_exonic`, `n_antisense_exonic`,
  `n_exonic_ambig`, `n_exon_incompatible`, `n_intronic`, `n_intergenic`
- Pool diagnostics: `n_pool_a`, `n_pool_b`, `pi_pool_b`,
  `fit_a_converged`, `fit_b_converged`, `fit_ab_kl_divergence`
- Cluster diagnostics: `n_clusters_total`, `n_clusters_gdna_dominant`

Fields **deleted** (no zombie compat): `region_e_gdna`, `region_n_total`,
`lambda_gdna`, `region_gamma`, `region_gamma_strand`, `mu_R`, `sigma_R`,
`mixing_pi`, `mixing_pi_soft`, `kappa_G`, `kappa_R`, `em_n_iter`,
`em_converged`, `n_eligible`, `n_soft`, `n_spliced_hard`, `lam_G_on`,
`lam_G_off`, `capture_class_mode`.

### 2.8 Why this works where v5 failed

| v5 failure | SRD behavior |
|---|---|
| Single global `λ_G` misspecified under hybrid capture | No density channel. Spatial variation irrelevant. |
| Poisson / LogNormal / Beta-Binomial / κ dispersion fitting | None. Histograms + closed-form per-fragment / per-cluster. |
| 50-iter fused-LLR EM with quadrature | Closed-form everywhere except the 1-D mixture (<30 iters, well-posed). |
| Brittle hyperparameters | Three knobs with robust defaults: `exon_fit_tolerance_bp`, `min_cluster_N`, `gamma_min`. |
| Regional γ from one global number | Per-cluster γ from observable strand imbalance; per-fragment γ from cluster + antisense channels. |
| No viable signal for unstranded libraries | EXON_INCOMPATIBLE pool + 1-D mixture. |
| Stranded "antisense purity ≈ 0.8" inflated by unannotated antisense RNA | Cluster-balance check filters clusters that don't match gDNA expectation. |
| Regional prior load-bearing but misfit | No regional prior at all. `compute_locus_priors` deleted. |

---

## 3 — Phased execution plan

Each phase ends with a green test run and a commit. Phase numbering matches
prior plans for cross-reference.

### Phase 0 — Audit (✅ COMPLETE)

Deliverable: [srd_v1_audit.md](srd_v1_audit.md). Findings locked in §1.

### Phase 1 — Park v5 in legacy branch

**Goal:** preserve a rollback path; main is not yet modified.

1. Create branch `legacy/calibration-v5-em` from current `main` head; push.
2. Record branch SHA in [srd_v1_audit.md](srd_v1_audit.md) §exit-summary.
3. Switch back to `main`. No file changes on main yet.

**Exit criteria:** branch pushed; SHA recorded.

### Phase 2 — Build SRD core (parallel to v5)

**Goal:** new code lives next to v5, fully unit-tested, zero pipeline impact.

**New files** under `src/rigel/calibration/`:

- `_categorize.py` — Pass 0. Per-fragment categorization. Reuses C++
  scanner outputs (`splice_type`, `ambig_strand`, `exon_strand`, `t_indices`).
  Adds the exon-fit geometry test with `exon_fit_tolerance_bp`.
- `_cluster.py` — Pass 1. Single-pass overlap-cluster builder. Returns
  `cluster_id[n_frags]` and per-cluster `(N_total, S, AS)` arrays.
- `_fit_a_strand.py` — Pass 2A. Closed-form γ per cluster
  ($\gamma = (\text{SS}-\text{sf})/(\text{SS}-0.5)$ clipped) + antisense
  channel. Builds weighted `gDNA_FL_strand` histogram.
- `_fit_b_mixture.py` — Pass 2B. 1-D mixture EM with `RNA_FL` fixed.
  Returns `(π_pool_b, gDNA_FL_geom, converged)`.
- `_fl_empirical_bayes.py` — thin wrapper around
  `FragmentLengthModel.finalize(prior_counts=..., prior_ess=...)`. Builds
  the global FL anchor and applies EB shrinkage to `RNA_FL` and `gDNA_FL`.
- `_simple.py` — orchestrator `calibrate_gdna(buffer, index, models, config)`.
  Runs Passes 0→4; constructs new `CalibrationResult`.

**New tests** under `tests/`:

- `test_categorize.py` — exhaustive synthetic exon layouts: per category;
  fragment fitting one exon (SENSE / ANTISENSE / AMBIG); fragment longer
  than its containing exon (EXON_INCOMPATIBLE); spanning exon-intron
  boundary (EXON_INCOMPATIBLE); within-tolerance overhang (valid SENSE);
  beyond-tolerance overhang (EXON_INCOMPATIBLE); INTRONIC; INTERGENIC;
  SPLICED bypass; opposite-strand overlapping transcripts → AMBIG.
- `test_cluster.py` — singleton cluster; chained overlaps; disjoint
  clusters; zero-length fragments; exact-boundary touches.
- `test_fit_a_strand.py` — injected strand imbalance; γ recovery within
  ±0.05; SS=0.5 produces no Pool A; antisense-channel weighting
  $(2 \cdot \text{SS} - 1)$ verified.
- `test_fit_b_mixture.py` — known-mixture recovery with KL tolerance and
  $\pi$ within ±0.02; degenerate cases ($\pi=0$, $\pi=1$, identical FLs).
- `test_fl_empirical_bayes.py` — small pool collapses to global FL; large
  pool approaches raw posterior; ESS behavior.

**Decisions baked in:**
- No SS-threshold branch.
- `gdna_fl_model_quality` emitted.
- `exon_fit_tolerance_bp` default 5.
- `min_cluster_N` default 5; `gamma_min` default 0.8; `fl_prior_ess`
  default 500; `kl_tolerance` default 0.1; `min_pool_b` default 300.

**Exit criteria:** new tests green; v5 still runs unchanged in main pipeline.

### Phase 3 — Wire SRD; replace `CalibrationResult`; gut deprecated knobs

**Goal:** SRD becomes the default; v5 code remains importable but unused
(deletion in Phase 4).

**Modified files:**

- `src/rigel/calibration/_result.py` — replace contents per §2.7. Drop
  every v5 field. Keep class name `CalibrationResult`.
- `src/rigel/calibration/__init__.py` — export `calibrate_gdna` from
  `_simple`. Drop v5 exports.
- `src/rigel/config.py` — remove dead v5 knobs (κ bracket, eligibility
  floor, capture-class flags, μ_R/σ_R seeds). Add a `CalibrationConfig`
  dataclass:
  - `exon_fit_tolerance_bp: int = 5`
  - `min_cluster_N: int = 5`
  - `gamma_min: float = 0.8`
  - `fl_prior_ess: float = 500.0`
  - `kl_tolerance: float = 0.1`
  - `min_pool_b: int = 300`
- `src/rigel/pipeline.py` — replace the `calibrate_gdna(...)` call site
  per the new signature. Remove `region_counts`/`fl_table` plumbing in
  `scan_and_buffer`. Remove `_compute_priors` helper. In
  `quant_from_buffer`, pass zero `α_gdna`/`α_rna` arrays to the C++ EM.
- `src/rigel/estimator.py` — remove any reads of `region_e_gdna` /
  `lambda_gdna`.
- `src/rigel/locus.py` — **delete** `compute_locus_priors`. (Body moved to
  `src/rigel/_deprecated/compute_locus_priors.py` for git history; not
  imported anywhere.)
- `src/rigel/native/*` — **no source changes** per Phase 0 audit. Recompile
  in Phase 4 only after dead-code removal.

**Tests:**

- `tests/test_calibration_simple.py` — six scenarios:
  1. **Stranded pure-RNA** (SS=0.95, zero gDNA): Pool A empty,
     `pi_pool_b ≈ 0`, `gdna_fl_model_quality = "fallback"`.
  2. **Stranded 20% gDNA** (SS=0.95, distinct gDNA FL): Fit A and Fit B
     both converge; KL within tolerance; fused gDNA_FL matches injection;
     `pi_pool_b ≈ 0.20 ± 0.02`.
  3. **Unstranded 20% gDNA** (SS=0.5): Fit A skipped; Fit B recovers
     $\pi \approx 0.20 \pm 0.03$ and gDNA_FL.
  4. **Hybrid-capture-like** (stranded, 100× density variation, uniform
     10% gDNA): SRD produces consistent cluster-level γ; no density
     artifact. v5 fails this.
  5. **QC-fail near-empty RNA**: graceful fallback, loud warning,
     `quality = "fallback"`.
  6. **Schema contract test**: `CalibrationResult` fields match new
     schema; no v5 fields present.
- Delete: `tests/test_calibration.py`, `test_calibration_integration.py`,
  `test_calibration_stress.py`, `test_capture_class.py`,
  `test_locus_priors.py`, `test_calibrated_density.py`.
- Salvage `_make_region_df` / `_make_region_counts` helpers from old
  `test_calibration.py` to a new `tests/_helpers.py` if useful.
- Rename `tests/test_partition.py` → `tests/test_locus_partition.py`
  (after Phase 3.5 below).

**Phase 3.5 — Rename `partition.py` → `locus_partition.py`:**

The Phase 0 audit confirmed `partition.py` is the per-locus CSR scatter
that builds `LocusPartition` objects, not region partitioning. Rename to
reflect actual role. This is a focused rename done as its own commit.

- `git mv src/rigel/partition.py src/rigel/locus_partition.py`
- `git mv tests/test_partition.py tests/test_locus_partition.py`
- Update imports: `from .partition import partition_and_free` →
  `from .locus_partition import partition_and_free` (one site in
  `pipeline.py`).
- Update docstring to clarify scope ("Per-locus CSR scatter for the
  per-locus EM. Unrelated to region/annotation partitioning.").

**Goldens:** regenerate `tests/golden/` with `pytest --update-golden` after
eyeballing diffs. Splice-anchored loci should not regress;
hybrid-capture loci should improve.

**Exit criteria:** `pytest tests/ -v` green; smoke `rigel quant` runs on a
stranded BAM and an unstranded BAM; eyeball `summary.json` for sanity.

### Phase 4 — Bulk deletion

**Goal:** the death phase. Done last so main is never half-broken.

**Move to `src/rigel/_deprecated/calibration_v5/`:**
- `src/rigel/calibration/_em.py`
- `src/rigel/calibration/_stats.py`
- `src/rigel/calibration/_fl_model.py` (old; superseded by `_fl_empirical_bayes.py`)
- `src/rigel/calibration/_calibrate.py`
- `src/rigel/calibration/_annotate.py`

**Move to `src/rigel/_deprecated/`:**
- The body of the deleted `compute_locus_priors` (already done in Phase 3,
  reaffirmed here).

**Delete from `src/rigel/index.py`:**
- `build_region_table` function (lines 525–641).
- `region_df` property (line 713).
- `region_cr` property (lines 725–738).
- `uniform_region_exposures` if unique-consumer (review line 891–916).
- The region build invocation in `TranscriptIndex` constructor.

**Delete from `src/rigel/native/bam_scanner.cpp`:**
- The `fl_table` / `region_id` accumulation paths.
- **KEEP** all per-fragment exon resolution, strand classification, and
  splice-type detection — these are SRD primitives.
- **Recompile** with `pip install --no-build-isolation -e .` and run
  `pytest tests/ -v`.

**Archive to `docs/deprecated/calibration_v5/`:**
- `docs/calibration/strand_channel_theory_and_redesign.md`
- `docs/calibration/density_nb_model_plan.md`
- `docs/calibration/gdna_locus_prior_deep_diagnostic.md`
- `docs/calibration/capture_class_density_plan.md`
- `docs/calibration/channel_fusion_hybrid_capture.md`
- `docs/calibration/effective_length_audit.md` (review for mappability content)
- `docs/calibration/gdna_implementation_log.md`
- `docs/calibration/gdna_length_model_options.md`
- `docs/calibration/root_cause_locus_em_nrna_siphon.md`
- `docs/calibration/SESSION_HANDOFF_2026-04-23.md`
- `docs/calibration/simplification_overhaul_v1.md`
- `docs/calibration/simplification_plan_other_claude.md`

Keep on `docs/calibration/`: `srd_v1_plan.md` (this doc), `srd_v1_audit.md`,
plus the new `srd_v1_results.md` from Phase 5.

**Archive scripts:**
- `scripts/calibration/stress_mini/` (uses v5 API).
- `scripts/benchmark/compare_vcap_optionb.py` — either update to read new
  SRD diagnostics or archive.

**Verify:**
- `pytest tests/ -v` green.
- `ruff check src/ tests/` clean.
- `grep -r "from rigel._deprecated" src/` returns nothing.
- `src/rigel/calibration/` LOC dropped from ~2100 to ~600–800.

**Exit criteria:** all of the above; commit pushed.

### Phase 5 — Integration validation

**Goal:** confirm SRD beats v5 on realistic data and produces no
regressions on splice-anchored loci.

1. Run all 13 Armis2 conditions:
   `python -m scripts.benchmarking run -c scripts/benchmarking/configs/default.yaml`
2. Analyze vs v5 baseline:
   `python -m scripts.benchmarking analyze -c scripts/benchmarking/configs/default.yaml -o results/srd_v1_report`
3. **Pass criteria:**
   - Stranded pure-RNA: relative-error change ≤ +5%.
   - Stranded + gDNA high-contamination: improvement vs v5.
   - Unstranded + gDNA: SRD produces non-trivial gDNA separation.
   - Hybrid-capture-like synthetic: SRD recovers near-uniform gDNA;
     `quant.gdna_total / strand_inferred_gdna < 1.1×`.
4. **Profile:**
   - Pass 0 (categorize) < 1 s / M fragments.
   - Pass 1 (cluster) < 200 ms / M fragments.
   - Fit B (mixture) < 100 ms.
   - Total calibration < 25% of v5 wallclock.

**Deliverable:** `docs/calibration/srd_v1_results.md` with before/after
tables and profile traces.

### Phase 6 — Docs & changelog

- Update `CLAUDE.md`, `README.md`, `docs/MANUAL.md`, `docs/METHODS.md` to
  describe SRD; remove v5 references.
- Bump version. Write `CHANGELOG.md` entry for the breaking changes:
  - `CalibrationResult` schema change (v5 fields removed).
  - `compute_locus_priors` removed.
  - `partition.py` renamed to `locus_partition.py`.
  - v5 config knobs removed; new `CalibrationConfig` knobs added.
- Archive v5 docs in `docs/deprecated/`.

### Phase 7 (deferred / separate effort) — Mappability-corrected effective length

Not part of this overhaul. Tracked as a follow-up: incorporate mappability
into the per-transcript effective length used by `_scoring_impl`'s
log-likelihood. Independent of calibration.

---

## 4 — Risk register

| Risk | Mitigation |
|---|---|
| Pass 0 misclassifies overlapping-gene fragments | Exhaustive `test_categorize.py`; reuse C++-side mixed-strand detection |
| Overhang tolerance too tight or too loose | Unit-tested at multiple values; default 5 bp; config knob |
| Cluster filter too strict on sparse loci | EB shrinkage to global FL; `gdna_fl_model_quality` exposes it |
| Cluster filter too loose → false-positive RNA inclusion | $\gamma_{\min}$ default 0.8; verified on stranded ground-truth data |
| Single-linkage creates mega-islands at deep loci | Optional cluster-span cap; defer until profiling proves it matters |
| Fits A and B disagree on real data | Log KL; prefer Fit A; surface as `fit_ab_kl_divergence` diagnostic |
| Golden test churn hides real regressions | Eyeball every golden diff before regenerating |
| Removing region partition breaks something not found in Phase 0 | Phase 4 done as separate commit so it can be reverted alone |
| C++ smoke recompile reveals a hidden region dependency | Recompile and run full suite at end of Phase 4; revert if broken |

---

## 5 — Out of scope (deliberate)

- Any density-based signal (including on `NONEXONIC`). Revisit only if
  diagnostics justify.
- Changes to the C++ per-locus EM scoring math.
- Mappability-corrected effective length (Phase 7 follow-up).
- Strand-model or RNA-FL training changes.
- Any protocol-detection (capture vs non-capture) heuristic — SRD is
  immune by construction.
- Backward-compat shims for v5 fields. No zombies.
