# SRD v1 — Phase 0 Audit

Read-only audit of the Rigel codebase establishing the kill-list, port-list,
deletion inventory, C++ dependency status, FL anchor status, and test
migration plan for the SRD v1 calibration overhaul.

This document is the contract for Phases 1–4. Any deviation in later phases
must update this audit first.

## Headline findings

1. **C++ requires NO changes.** The C++ EM solver and scorer receive
   `gdna_fl_model` (a `FragmentLengthModel` object) and per-locus `α_gdna` /
   `α_rna` arrays. They are agnostic to how these are computed. Replacing the
   v5 calibrator with SRD only requires swapping the Python prior-computation
   path. Verified by reading `pipeline.py` lines 555–580, 730–850; `estimator.py`
   lines 740, 820–900; `scoring.py` `FragmentScorer` constructor.

2. **`partition.py` STAYS.** Despite the name, this module is the per-locus
   CSR array scatter (`scatter_candidates_*`, `scatter_units_*`). It is
   unrelated to region partitioning and is required by the per-locus EM. All
   tests in `tests/test_partition.py` remain valid.

3. **`region_id` is NOT in FragmentBuffer.** Buffer schema has only per-fragment
   transcript candidate set + strand/splice fields. Region IDs only exist in
   the optional `fl_table` produced by the scanner for v5 calibration. Removing
   region machinery has zero impact on buffer consumers downstream.

4. **Global FL anchor is ALREADY trained.** `FragmentLengthModels.global_model`
   exists, is populated during the BAM scan, and is *already* used as the
   Dirichlet prior when finalizing `gdna_model` (pipeline.py:940–950 via
   `finalize(prior_counts=..., prior_ess=...)`). Phase 2 does not need new
   training code — only a new way to populate the gDNA histogram.

5. **`compute_locus_priors` is the sole downstream consumer of
   `region_e_gdna` / `region_n_total`.** No other module reads these fields.
   Deleting `compute_locus_priors` removes the entire region→prior data flow
   in one step.

---

## Deliverable 1 — `CalibrationResult` field kill-list

Source: [src/rigel/calibration/_result.py](src/rigel/calibration/_result.py).
Consumers grepped across `src/`, `tests/`, `scripts/`.

| Field | Read by | Action |
|---|---|---|
| `region_e_gdna` | `_calibrate.py:81`, `locus.py:175` (`compute_locus_priors`), `pipeline.py:434` (via `_compute_priors`) | **DELETE** (with `compute_locus_priors`) |
| `region_n_total` | `_calibrate.py:90`, `locus.py:176` | **DELETE** |
| `gdna_fl_model` | `pipeline.py:942` (passed to scorer), `_calibrate.py:90`, `_result.py:60` | **KEEP** — primary product |
| `lambda_gdna` | `_result.py:50`, `pipeline.py:928` (logging only), `_calibrate.py:123` | **DELETE** — no functional consumer |
| `strand_specificity` | `pipeline.py:931`, `_result.py:51`, `_calibrate.py:124` | **KEEP** — used by strand model wiring |
| `region_gamma` | written only | **DELETE** |
| `region_gamma_strand` | written only | **DELETE** |
| `mu_R` | written only | **DELETE** |
| `sigma_R` | written only | **DELETE** |
| `mixing_pi` | written only; `scripts/benchmark/compare_vcap_optionb.py:75` (script) | **DELETE** + update script |
| `mixing_pi_soft` | written only | **DELETE** |
| `kappa_G` | written only | **DELETE** |
| `kappa_R` | written only | **DELETE** |
| `em_n_iter` | written only | **DELETE** |
| `em_converged` | written only | **DELETE** |
| `n_eligible` | written only | **DELETE** |
| `n_soft` | written only | **DELETE** |
| `n_spliced_hard` | written only | **DELETE** |
| `lam_G_on` / `lam_G_off` | written only | **DELETE** |
| `capture_class_mode` | written only | **DELETE** |

**New SRD fields to add** (per the plan):
`gdna_fl_strand`, `kl_strand_vs_cluster`, category counts, `n_balanced_clusters`,
`n_fragments_in_balanced_clusters`, `gdna_fl_quality` enum.

---

## Deliverable 2 — Port-list (per-fragment classification primitives)

The new `_categorize.py` does NOT need to re-implement transcript resolution.
The C++ BAM scanner already produces every primitive we need; we only need
to add the **exon-fit geometry test** in Python and consume existing fields.

### Reusable primitives already produced by C++ scanner / FragmentBuffer

Source: [src/rigel/buffer.py](src/rigel/buffer.py), populated by
[src/rigel/native/bam_scanner.cpp](src/rigel/native/bam_scanner.cpp).

| Primitive | Buffer field | Where set | Reuse |
|---|---|---|---|
| Splice status (`SPLICED_ANNOT` / `SPLICED_UNANNOT` / `UNSPLICED`) | `splice_type` | C++ scanner CIGAR analysis | Direct read; SPLICED → fixed category |
| Transcript candidate set | `t_indices` (CSR) | C++ scanner via `_resolve_impl` | Used to look up tx_strand |
| Mixed-strand exon overlap | `ambig_strand` | C++ scanner | **Direct port** — drives `UNSPLICED_EXONIC_AMBIG` |
| Exon strand of overlap | `exon_strand` | C++ scanner | Direct read |
| SJ strand | `sj_strand` | C++ scanner (annotated SJ lookup) | Direct read |
| Fragment genomic span `[g0, g1]` | derived from buffer alignment columns | C++ scanner | Used for exon-fit test |
| Strand qualification compound check | `is_strand_qualified` property | `buffer.py:110-120` | Reference logic for sense/antisense determination |

### What we MUST add in Python (`_categorize.py`)

1. **Exon-fit test with overhang tolerance**
   ```python
   def fits_in_exon(g0: int, g1: int, exon_starts: np.ndarray,
                    exon_ends: np.ndarray, tol: int = 5) -> int:
       """Returns index of containing exon, or -1."""
   ```
   Vectorized via cgranges interval lookup against exon table. The exon
   interval index already exists in `TranscriptIndex` (used for resolution).

2. **Strand-of-tx lookup with mixed-strand handling**
   - Already done by C++ via `ambig_strand`. Python just consumes.

3. **Intronic vs intergenic split**
   - "Intronic" = no exon overlap but ≥1 transcript span overlap.
   - "Intergenic" = no transcript overlap at all.
   - Transcript-span overlap can be derived from `t_indices` length (0 → intergenic).

### Code that should NOT be ported (deprecated logic patterns)

| File | Function | Why not |
|---|---|---|
| `_stats.py` | `compute_region_stats` | Operates on per-region aggregates; SRD is per-fragment |
| `_stats.py` | `compute_sense_fraction` | Per-region; subsumed by per-fragment `exon_strand` |
| `_annotate.py` | `annotate_capture_class` | Capture-class is deleted entirely |
| `_em.py` | `run_em` and helpers | Entire mixture EM is replaced |

---

## Deliverable 3 — Deletion inventory

### A. Move to `src/rigel/_deprecated/calibration_v5/`

| Path | Reason |
|---|---|
| `src/rigel/calibration/_em.py` | v5 region EM mixture (largest single deletion) |
| `src/rigel/calibration/_stats.py` | Region stats |
| `src/rigel/calibration/_fl_model.py` | v5 region-weighted gDNA FL builder (replaced by `_fl_eb.py`) |
| `src/rigel/calibration/_calibrate.py` | v5 orchestrator (replaced by `_simple.py`) |
| `src/rigel/calibration/_annotate.py` | Capture-class annotation |

### B. Move to `src/rigel/_deprecated/`

| Path | Reason |
|---|---|
| `src/rigel/locus.py::compute_locus_priors` (and `_compute_priors` wrapper in `pipeline.py`) | Sole consumer of `region_e_gdna`/`region_n_total`. Move as standalone module `_deprecated/locus_priors.py`. |

### C. Delete from `src/rigel/index.py`

| Symbol | Lines | Action |
|---|---|---|
| `build_region_table` | 525–641 | Delete (move to deprecated module if we want history) |
| `region_df` property | 713 | Delete |
| `region_cr` property | 725–738 | Delete |
| `uniform_region_exposures` (if only consumer is build flow) | 891–916 area | Delete |
| Region build invocation in `TranscriptIndex` constructor | 891–916 | Delete |

### D. Delete from `src/rigel/pipeline.py`

| Symbol | Action |
|---|---|
| `region_counts` / `fl_table` build path in `scan_and_buffer` (lines 246, 298–304, 338) | Delete |
| `_compute_priors` wrapper (line 430-area) | Delete |
| Calibration call site (line 893–915) | Replace with new `calibrate_gdna(buffer, index, models, config)` |
| `quant_from_buffer` priors plumbing | Pass zeros (or skip) when `calibration_prior_strength == 0.0` |

### E. C++ — region-related code

| File | Action |
|---|---|
| `src/rigel/native/bam_scanner.cpp` — `fl_table` / `region_id` emission paths | Delete the per-region FL table accumulation; **keep** all per-fragment exon resolution and strand classification (still needed). |
| Per-fragment exon overlap, strand classification, splice-type detection | **KEEP** (these are SRD primitives) |

A C++ smoke recompile is required after Phase 4 to confirm no header dependency was missed.

### F. Tests — full kill / migrate map

| Test file | Action | Notes |
|---|---|---|
| `tests/test_calibration.py` | **DELETE** | Pure v5 region EM tests. Salvage the `_make_region_df` / `_make_region_counts` synthetic-data helpers into a test fixture module if useful for SRD tests. |
| `tests/test_calibration_integration.py` | **DELETE** | v5 end-to-end. |
| `tests/test_calibration_stress.py` | **DELETE** | v5 EM stress. |
| `tests/test_capture_class.py` | **DELETE** | Capture-class is gone. |
| `tests/test_locus_priors.py` | **DELETE** | `compute_locus_priors` is gone. |
| `tests/test_calibrated_density.py` | **DELETE** | v5 prior plumbing. |
| `tests/test_partition.py` | **KEEP** | Per-locus CSR scatter; unrelated to regions. |
| `tests/test_golden_output.py` | **KEEP, regenerate goldens** | Eyeball diffs before `--update-golden`. |
| All other `tests/test_*.py` | **KEEP** | Not calibration-coupled. |

**New test files for Phase 2/3:**
- `tests/test_categorize.py` — exhaustive exon-fit / strand category tests
- `tests/test_clusters.py` — balanced-cluster Binomial QC
- `tests/test_fl_eb.py` — EB shrinkage behavior
- `tests/test_fl_mixture.py` — 1-D mixture diagnostic
- `tests/test_calibration_simple.py` — end-to-end SRD scenarios

### G. Docs — archive to `docs/deprecated/calibration_v5/`

- `docs/calibration/strand_channel_theory_and_redesign.md`
- `docs/calibration/density_nb_model_plan.md`
- `docs/calibration/gdna_locus_prior_deep_diagnostic.md`
- `docs/calibration/capture_class_density_plan.md`
- `docs/calibration/channel_fusion_hybrid_capture.md`
- `docs/calibration/effective_length_audit.md` (review — may have content relevant to mappability/effective-length follow-up)
- `docs/calibration/gdna_implementation_log.md`
- `docs/calibration/gdna_length_model_options.md`
- `docs/calibration/root_cause_locus_em_nrna_siphon.md`
- `docs/calibration/SESSION_HANDOFF_2026-04-23.md`

### H. Scripts — archive

- `scripts/calibration/stress_mini/` — uses v5 calibrator API.
- `scripts/benchmark/compare_vcap_optionb.py` — uses `mixing_pi` field; either update to read new SRD diagnostics or archive.

---

## Deliverable 4 — C++ dependency audit (CRITICAL)

### What flows from calibration → C++ EM solver

Read from `pipeline.py` lines 555–580 (FragmentScorer construction), 730–850
(`_run_locus_em_partitioned`), and `estimator.py` lines 740, 820–900
(`run_batch_locus_em_partitioned`):

```
FragmentScorer(
    log_p_sense=...,             # from StrandModels (NOT calibration)
    log_p_antisense=...,
    r1_antisense=...,
    overhang_log_penalty=...,    # from FragmentScoringConfig
    mismatch_log_penalty=...,
    splice_type_penalties=...,
    gdna_fl_model=                # ← FROM CALIBRATION
        frag_length_models.gdna_model,
)

run_batch_locus_em_partitioned(
    partition_tuples,             # CSR (offsets, indices, log_liks, ...)
    locus_t_lists,
    alpha_gdna,                   # ← per-locus prior, from compute_locus_priors
    alpha_rna,                    # ← per-locus prior, from compute_locus_priors
    em_iterations,
    em_convergence_delta,
    ...
)
```

### Verdict per `CalibrationResult` field

| Field | Read by C++ EM? | Read by C++ scorer? | Notes |
|---|---|---|---|
| `gdna_fl_model` | indirectly (via `gdna_log_liks` already pre-computed by scorer in Python) | **YES** | The scorer reads this `FragmentLengthModel` object |
| `lambda_gdna` | NO | NO | Only used for logging |
| `region_e_gdna` / `region_n_total` | NO | NO | Only used by `compute_locus_priors` to produce `α` arrays |
| `region_gamma*`, `mu_R`, `sigma_R`, `mixing_pi*`, `kappa_*`, `lam_G_on/off`, `capture_class_mode` | NO | NO | Diagnostic only |

### Conclusion

> **The C++ layer is calibration-state-agnostic.** It receives only:
> 1. A `FragmentLengthModel` object for gDNA likelihoods (kept in SRD).
> 2. Per-locus `α_gdna`, `α_rna` float arrays (in SRD: zeros by default since
>    `calibration_prior_strength = 0.0`).
>
> Replacing v5 with SRD requires **zero C++ source changes.** The only C++
> work in Phase 4 is removing the dead `fl_table`/`region_id` emission code
> in `bam_scanner.cpp` and recompiling.

---

## Deliverable 5 — Global FL anchor status

`FragmentLengthModels` (definition trail: `frag_length_model.py` plus
`pipeline.py:110-150`):

| Attribute | Trained from | SRD use |
|---|---|---|
| `exonic_spliced` | unique-mapper RNA, spliced | Source of `RNA_FL` |
| `exonic_unspliced` | unique-mapper RNA, unspliced exonic | Diagnostic |
| `intergenic` | unique-mapper intergenic fragments | Diagnostic |
| `global_model` | aggregate of all three | **EB shrinkage anchor for both `RNA_FL` and `gDNA_FL`** |
| `gdna_model` | populated by calibration | **OUTPUT: SRD writes here** |

Methods to reuse:
- `observe(length, splice_type)` / `observe_intergenic(length)` — already wired
- `finalize(prior_counts=..., prior_ess=...)` — **this is the EB shrinkage
  primitive** the new `_fl_eb.py` will call
- `log_likelihood(length)`, `mean`, `var`, `total_weight`

**Phase 2 implication:** No new training code is needed. `_fl_eb.py` becomes a
thin wrapper that constructs a `FragmentLengthModel`, observes the
balanced-cluster pool length histogram, and calls `finalize(prior_counts=
global_model.counts, prior_ess=cfg.fl_prior_ess)`. Same pattern for the
diagnostic strand-pool variant.

---

## Deliverable 6 — Region-related fields in FragmentBuffer

Grep of [src/rigel/buffer.py](src/rigel/buffer.py) for `region`: **zero matches.**

The buffer schema has no `region_id` column. Region IDs exist only in the
optional `fl_table: pd.DataFrame` returned alongside the buffer from
`scan_and_buffer` (pipeline.py:298–338) and consumed only by `_calibrate.py`.

Removing `fl_table` from the scanner's return tuple cleanly excises all
region propagation from the runtime.

---

## Deliverable 7 — Test deletion / migration table

Final classification (consolidated; see Deliverable 3F):

| File | Status |
|---|---|
| `tests/test_calibration.py` | DELETE; salvage `_make_region_df` / `_make_region_counts` to `tests/_helpers.py` if useful |
| `tests/test_calibration_integration.py` | DELETE |
| `tests/test_calibration_stress.py` | DELETE |
| `tests/test_capture_class.py` | DELETE |
| `tests/test_locus_priors.py` | DELETE |
| `tests/test_calibrated_density.py` | DELETE |
| `tests/test_partition.py` | KEEP (per-locus EM scatter, not region partitioning) |
| `tests/test_golden_output.py` | KEEP; regenerate goldens in Phase 3 after eyeballing diffs |
| `tests/test_categorize.py` | NEW (Phase 2) |
| `tests/test_clusters.py` | NEW (Phase 2) |
| `tests/test_fl_eb.py` | NEW (Phase 2) |
| `tests/test_fl_mixture.py` | NEW (Phase 2) |
| `tests/test_calibration_simple.py` | NEW (Phase 3) |

---

## Phase 0 exit summary

- C++ requires zero source changes (only dead-code removal in `bam_scanner.cpp`).
- `partition.py` stays; misnamed but unrelated to regions.
- `global_model` already exists as the EB anchor — Phase 2 needs no new
  training code.
- `compute_locus_priors` is the *sole* downstream consumer of region state;
  removing it cleanly severs the entire region data flow.
- ~3500–4000 LOC scheduled for deletion (calibration v5 + region machinery
  + dead tests).
- New SRD calibration module budget: ~600–800 LOC across
  `_categorize.py`, `_clusters.py`, `_fl_eb.py`, `_fl_mixture.py`, `_simple.py`,
  plus a leaner `_result.py`.

**Ready to proceed to Phase 1 (park v5 in legacy branch).**
