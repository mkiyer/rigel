# SRD v1 — Revised Plan (post-amendments)

This document captures the revised SRD overhaul plan after incorporating the
ten user amendments to `simplification_overhaul_v1.md`. Phase 0 audit notes
are written separately to `srd_v1_audit.md`.

## Revised commentary on the new ideas

### Additional Point 4 — strand-balanced unspliced clusters (keystone)

The reasoning: **gDNA is intrinsically unstranded.** A genuine gDNA-dominated
locus must show sense/antisense ≈ 50/50 in unspliced reads (within Binomial
noise on `n`). Conversely, any cluster with a strong sense skew is
RNA-contaminated, regardless of nominal strand specificity.

This converts the impurity-at-low-SS problem into a *cluster-level QC filter*:

- Build "unspliced coverage islands" by single-linkage on overlapping
  unspliced fragments (cheap; one sweep).
- For each island compute `n`, `S`, `AS`. Accept iff
  `Binomial(S | n, 0.5)` is *not* rejected at α (e.g. two-sided p > 0.01),
  with a minimum-`n` floor (≥ 10).
- All fragments inside accepted islands → contribute to `gDNA_FL`. Strand
  label is irrelevant once the cluster is accepted.
- Library-agnostic (works at SS=0.5, 0.7, 1.0). Subsumes antisense-exonic
  pool naturally: a 100% antisense exon island fails balance; a balanced
  exon-incompatible island passes.

**This becomes the primary gDNA pool** in v1; antisense-exonic is kept as a
diagnostic secondary pool.

### Empirical-Bayes FL shrinkage (Point 10)

Train a single **global FL** model from all unique-mapper fragments. Let
`gDNA_FL_raw` be the histogram from the accepted-cluster pool. Final
`gDNA_FL` is the posterior of a Dirichlet–Multinomial with prior ESS
`α₀ ≈ 500` centered on `global_FL`:

$$
\text{gDNA\_FL}[l] = \frac{\text{counts}[l] + \alpha_0 \cdot \text{global\_FL}[l]}{\sum_l \text{counts}[l] + \alpha_0}
$$

Same shrinkage applied to `RNA_FL` (from spliced reads). Pool size << α₀ →
collapses to `global_FL` (kills FL channel safely). Pool size >> α₀ →
posterior ≈ raw histogram. One parameter, clear semantics. Replaces both
"fallback on small pools" and "QC-fail RNA ≈ 0" graceful-degradation
requirements.

### Fused stranded/unstranded — concretely

Run both paths every time:
- `gdna_fl_strand`: built from antisense-exonic pool weighted by `(2·SS−1)`.
- `gdna_fl_cluster`: built from balanced-cluster pool (keystone).

Default consumed by downstream EM = `gdna_fl_cluster`. Both shipped in
`CalibrationResult` for diagnostics, plus KL divergence as a confidence
score. Once they agree on real data we may drop the strand variant in v2.

### EXON_INCOMPATIBLE tolerance

Single config knob `exon_overhang_tolerance_bp: int = 5` (absorbs soft-clip
residue, indels, aligner edge cases). Fragment fits iff
`g0 ≥ exon_start − tol ∧ g1 ≤ exon_end + tol` for some exon.

### Mappability (Point 9)

Out of scope here. Moves to scorer's effective-length term in a separate
effort. Calibration loses *all* mappability dependencies — clean break.

### Reuse of region-partitioning code

Per-fragment exon-fit / strand-overlap classification logic lives inside
`partition.py` / `_annotate.py` / C++ region resolution. We **port the
inner per-fragment logic** into the new categorizer; we **delete the
region-aggregation scaffolding**.

---

## Revised phased implementation plan

### Phase 0 — Audit (no code changes)

Goals:
1. Kill-list: every `CalibrationResult` field's consumers.
2. Port-list: per-fragment classification primitives in
   `partition.py`/`_annotate.py`/C++ region resolution that we'll re-use.
3. Deletion inventory: full list of v5 + region-partitioning files.
4. C++ scorer audit: confirm `_scoring_impl` consumes only `gdna_fl_model`
   (and not `lambda_gdna` or per-region quantities).
5. Global FL anchor status: confirm pipeline trains a `global_FL` usable
   for EB shrinkage; if only spliced-conditional models exist, schedule
   training code in Phase 2.

**Deliverable:** `docs/calibration/srd_v1_audit.md`.

### Phase 1 — Park v5 in `legacy/calibration-v5-em` branch

Tag and push. No deletion on main yet.

### Phase 2 — Build SRD core (parallel to v5)

New files:
- `src/rigel/calibration/_categorize.py` — Pass 0 categorizer.
  Categories: `SPLICED`, `UNSPLICED_SENSE_EXONIC`,
  `UNSPLICED_ANTISENSE_EXONIC`, `UNSPLICED_EXONIC_AMBIG`,
  `EXON_INCOMPATIBLE`, `INTRONIC`, `INTERGENIC`. Mixed-strand exon
  overlap → AMBIG. Numpy first, ports per-fragment logic from partition
  framework.
- `src/rigel/calibration/_clusters.py` — keystone module. Single-linkage
  unspliced coverage islands; per-island `n`, `S`, `AS`, Binomial(0.5)
  two-sided p-value. Returns per-fragment `in_balanced_cluster` boolean.
- `src/rigel/calibration/_fl_eb.py` — EB shrinkage
  `shrink_fl(counts, global_fl, alpha0) -> FragmentLengthModel`.
- `src/rigel/calibration/_fl_mixture.py` — 1-D mixture EM (cross-validation
  diagnostic; not load-bearing in v1).
- `src/rigel/calibration/_simple.py` — orchestrator `calibrate_gdna(...)`:
  1. Pass 0 categorize.
  2. Pass 1a: `gdna_fl_cluster` from balanced-cluster pool.
  3. Pass 1b: `gdna_fl_strand` from antisense-exonic pool (diagnostic).
  4. EB shrink both, anchored by `global_fl`.
  5. KL(`cluster` ‖ `strand`) for confidence.

New tests: `test_categorize.py`, `test_clusters.py`, `test_fl_eb.py`,
`test_fl_mixture.py`.

### Phase 3 — Wire SRD; replace `CalibrationResult`; gut deprecated knobs

Modified:
- `_result.py` — new fields: `gdna_fl_model` (= cluster), `gdna_fl_strand`,
  `kl_strand_vs_cluster`, category counts, `n_balanced_clusters`,
  `n_fragments_in_balanced_clusters`, `gdna_fl_quality` enum
  `{good, weak, fallback_global}`. Drop *every* v5 field.
- `__init__.py` — export `calibrate_gdna` from `_simple` only.
- `config.py` — remove all v5 knobs. Add `CalibrationConfig`:
  `exon_overhang_tolerance_bp: int = 5`,
  `cluster_min_fragments: int = 10`,
  `cluster_strand_balance_alpha: float = 0.01`,
  `fl_prior_ess: float = 500.0`.
- `pipeline.py` — call new `calibrate_gdna`; thread `global_fl`; remove
  `region_counts`/`fl_table` plumbing.
- `locus.py::compute_locus_priors` → moved to `_deprecated/`.
- `estimator.py` — drop prior consumption.

Tests: replace `test_calibration*.py` with `test_calibration_simple.py`
covering six original scenarios plus SS=0.7 cluster vs strand agreement,
KL threshold check, QC-fail RNA-near-zero fallback. Regenerate goldens.

### Phase 4 — Bulk deletion

Move to `src/rigel/_deprecated/`:
- `calibration/_em.py`, `_stats.py`, `_fl_model.py` (old), `_calibrate.py`,
  `_annotate.py`
- `partition.py`
- `compute_locus_priors` and helpers from `locus.py`

Delete from main:
- Region machinery in `index.py` and `pipeline.py`.
- Region-related fields in `FragmentBuffer` if no other consumer.
- Any C++ region-aggregation code (per Phase 0). Per-fragment exon-fit
  primitive stays.
- Dead docs to `docs/deprecated/`.
- Dead scripts (`scripts/calibration/stress_mini/`).

### Phase 5 — Integration validation

Run all 13 Armis2 conditions; analyze; report cluster-vs-strand KL per
condition. Pass criteria:
- Stranded pure-RNA: relative-error change ≤ +5%.
- Stranded with gDNA: improvement vs v5.
- Hybrid-capture / VCaP titration: SRD recovers near-uniform gDNA;
  `quant.gdna_total / strand_inferred_gdna < 1.1×`.
- Profile: categorize + cluster < 1 s / M fragments; total calibration
  < 25 % of v5 wallclock.

### Phase 6 — Docs & changelog

Hard-breaking release. CHANGELOG documents removed config knobs and
deleted result fields. Update `CLAUDE.md`, `README.md`, `docs/MANUAL.md`,
`docs/METHODS.md`. Archive v5 docs.

### Phase 7 (deferred) — Mappability-corrected effective length

Separate follow-up effort in the scorer's effective-length term. Not part
of this overhaul.

## Risk register

| Risk | Mitigation |
|---|---|
| Cluster-balance filter too strict on sparse loci → tiny gDNA pool | EB shrinkage to `global_fl`; `gdna_fl_quality` flag exposes it |
| Cluster-balance filter too loose → false-positive RNA inclusion | α=0.01 + n≥10 floor; verified on known-stranded data |
| Single-linkage creates mega-islands at deep loci | Cap island span or split on coverage troughs; defer until profiling proves it's a problem |
| Ported per-fragment logic drags in coupling | Phase 0 identifies clean cut points; re-implement rather than wrap |
| C++ scoring still consumes a deprecated field | Phase 0 explicit audit + C++ smoke compile after Phase 4 deletion |
