# Source Review — `src/rigel`

**Date:** 2026-04-20  
**Scope:** source-only review of `src/rigel/`  
**Goal:** identify the highest-value code quality issues blocking a more elegant, concise, readable, and durable codebase.

## Overall assessment

The strongest parts of this codebase are its clear top-level intent, solid configuration layer, and a purposeful columnar/native-performance architecture. The biggest threats to long-term quality are not superficial style problems. They are:

1. contract drift between public APIs and real runtime requirements,
2. semantic drift in the transcript/nRNA data model,
3. dict-shaped internal interfaces that are already starting to rot,
4. legacy paths that survived architecture shifts and now add cognitive weight.

The prioritized findings below focus on those problems.

## Prioritized findings

### 1. Enforce index manifest/version compatibility during load

**Severity:** High  
**Status:** New finding

**Evidence:**
- [../../src/rigel/index.py#L67](../../src/rigel/index.py#L67)
- [../../src/rigel/index.py#L78](../../src/rigel/index.py#L78)
- [../../src/rigel/index.py#L907](../../src/rigel/index.py#L907)
- [../../src/rigel/index.py#L948](../../src/rigel/index.py#L948)

The index builder writes a manifest with `format_version`, and the module exposes `load_manifest()`, but `TranscriptIndex.load()` does not validate that manifest before reading the persisted feather files. That weakens the on-disk contract and makes schema drift more likely to fail as partial compatibility or late runtime errors instead of a clean refusal.

**Why it matters:** persistent formats are part of the public API. Silent or delayed incompatibility is expensive to debug and will age badly as the index evolves.

**Improvement direction:** make manifest validation the first step of `TranscriptIndex.load()`, reject unsupported or missing versions unless a migration path is explicitly supported, and centralize required-column checks in one place.

### 2. Synthetic nRNA rows still encode single-gene ownership in the canonical state space

**Severity:** High  
**Status:** Already tracked in transcript-centric cleanup

**Evidence:**
- [../../src/rigel/index.py#L173](../../src/rigel/index.py#L173)
- [../../src/rigel/index.py#L322](../../src/rigel/index.py#L322)
- [../../src/rigel/index.py#L323](../../src/rigel/index.py#L323)
- [../../src/rigel/index.py#L324](../../src/rigel/index.py#L324)
- [../../src/rigel/index.py#L925](../../src/rigel/index.py#L925)

`create_nrna_transcripts()` still assigns synthetic implied nRNA rows the `g_id`, `g_name`, and `g_type` of one representative transcript. The derived gene table then groups the whole transcript table by `g_index`. That bakes a convenience projection into the core index, even though the project’s stated design is transcript-centric and overlapping loci make one-to-one gene ownership unreliable.

**Why it matters:** this is exactly the kind of hidden semantic shortcut that spreads wrong assumptions into reporting, locus metadata, and future model changes.

**Improvement direction:** keep synthetic implied nRNA rows gene-neutral in core tables and derive gene-oriented summaries after quantification from contributor metadata.

### 3. `is_nrna` currently conflates structure, annotation role, and modeled pool membership

**Severity:** High  
**Status:** Already tracked in transcript-centric cleanup

**Evidence:**
- [../../src/rigel/transcript.py#L172](../../src/rigel/transcript.py#L172)
- [../../src/rigel/scoring.py#L210](../../src/rigel/scoring.py#L210)
- [../../src/rigel/annotate.py#L99](../../src/rigel/annotate.py#L99)
- [../../src/rigel/estimator.py#L366](../../src/rigel/estimator.py#L366)
- [../../src/rigel/estimator.py#L539](../../src/rigel/estimator.py#L539)

Single-exon transcripts are marked `is_nrna` during GTF ingest, annotation flags treat those rows as nRNA-like outputs, but scoring explicitly says only synthetic rows belong to the nRNA pool. One flag is carrying multiple meanings depending on subsystem.

**Why it matters:** this is a semantic footgun. Any future change in reporting, routing, or scoring can easily interpret the flag the wrong way and produce a biologically inconsistent result.

**Improvement direction:** split this into explicit orthogonal metadata such as `is_single_exon`, `is_nascent_equiv`, and `is_synthetic_nrna`, then let each subsystem consume the flag it actually means.

### 4. `quant_from_buffer()` advertises optional calibration but dereferences it unconditionally

**Severity:** High  
**Status:** Already tracked in transcript-centric cleanup

**Evidence:**
- [../../src/rigel/pipeline.py#L666](../../src/rigel/pipeline.py#L666)
- [../../src/rigel/pipeline.py#L725](../../src/rigel/pipeline.py#L725)

`quant_from_buffer()` still exposes `calibration: "CalibrationResult" = None`, but the implementation immediately accesses `calibration.gdna_fl_model`. The signature and docstring say one thing; the runtime contract says another.

**Why it matters:** this is a public API trap. It guarantees confusing failures for direct callers and makes the code harder to trust.

**Improvement direction:** make calibration mandatory in the signature or fail immediately with a clear domain-specific error before any dereference.

### 5. The legacy per-locus EM extraction path is still present but appears unused

**Severity:** Medium-High  
**Status:** Partially tracked; dead-path aspect is new

**Evidence:**
- [../../src/rigel/locus.py#L126](../../src/rigel/locus.py#L126)
- [../../src/rigel/pipeline.py#L378](../../src/rigel/pipeline.py#L378)

`build_locus_em_data()` is a large, nontrivial implementation, but the live pipeline now goes through `_score_fragments()` and the partitioned native EM path. A search in `src/rigel` shows the function definition, but no production call site.

**Why it matters:** dead complex code is worse than simple dead code. It creates false architectural branches, increases review burden, and makes future refactors riskier because readers must first determine which path is real.

**Improvement direction:** remove the legacy path if the partitioned architecture is canonical, or move it into a clearly labeled compatibility/experimental module.

### 6. Internal locus-result metadata already shows schema drift

**Severity:** Medium  
**Status:** New finding

**Evidence:**
- [../../src/rigel/pipeline.py#L515](../../src/rigel/pipeline.py#L515)
- [../../src/rigel/pipeline.py#L525](../../src/rigel/pipeline.py#L525)
- [../../src/rigel/estimator.py#L688](../../src/rigel/estimator.py#L688)

The pipeline writes per-locus metadata as plain dictionaries with keys like `alpha_gdna` and `alpha_rna`, while `get_loci_df()` still documents and reads `gdna_init`. That is already a drifted internal schema.

**Why it matters:** dict-shaped internal protocols do not fail fast when they drift. They quietly produce wrong defaults, missing output, and long-tail maintenance bugs.

**Improvement direction:** replace the raw locus-result dicts with a dedicated dataclass or named tuple and normalize terminology across prior, init, and alpha.

### 7. Annotation-flag translation is duplicated in two different execution paths

**Severity:** Medium  
**Status:** New finding

**Evidence:**
- [../../src/rigel/scan.py#L189](../../src/rigel/scan.py#L189)
- [../../src/rigel/pipeline.py#L438](../../src/rigel/pipeline.py#L438)
- [../../src/rigel/pipeline.py#L477](../../src/rigel/pipeline.py#L477)
- [../../src/rigel/annotate.py#L92](../../src/rigel/annotate.py#L92)

The logic mapping winner transcript state to `AF_TRANSCRIPT`, `AF_NRNA_RESOLVED`, and `AF_SYNTH_RESOLVED` appears once in deterministic scan-time annotation and again in EM postprocessing.

**Why it matters:** duplicated domain policy is a slow-burn maintenance bug. Any future semantic change to what qualifies as nRNA or synthetic can easily update one path and miss the other.

**Improvement direction:** move winner-to-flag translation into one shared helper beside the flag definitions in `annotate.py` and call it from both paths.

### 8. Reporting code in `AbundanceEstimator` is too tightly coupled to raw `t_df` schema

**Severity:** Medium  
**Status:** Partially tracked in existing cleanup notes

**Evidence:**
- [../../src/rigel/estimator.py#L372](../../src/rigel/estimator.py#L372)
- [../../src/rigel/estimator.py#L378](../../src/rigel/estimator.py#L378)
- [../../src/rigel/estimator.py#L563](../../src/rigel/estimator.py#L563)
- [../../src/rigel/estimator.py#L566](../../src/rigel/estimator.py#L566)

The reporting methods repeatedly reach into `index.t_df`, branch on optional columns, duplicate TPM math, and rebuild nRNA-child relationships with ad hoc scans. The cleanup plan already flags some local duplication, but the larger issue is that estimator output depends directly on low-level dataframe schema details.

**Why it matters:** every index-schema or output-rule change fans out into multiple reporting methods. That is the opposite of stable, award-winner code.

**Improvement direction:** introduce a stable reporting/access layer from `TranscriptIndex`, centralize TPM computation, and precompute reusable nRNA-parent relationships instead of reconstructing them inside output methods.

### 9. `create_nrna_transcripts()` has a misleading return contract

**Severity:** Medium-Low  
**Status:** New finding

**Evidence:**
- [../../src/rigel/index.py#L173](../../src/rigel/index.py#L173)
- [../../src/rigel/index.py#L200](../../src/rigel/index.py#L200)
- [../../src/rigel/index.py#L337](../../src/rigel/index.py#L337)

The type signature says the function returns two values, the docstring documents four, and the implementation returns four.

**Why it matters:** this is concrete refactor drift in a core builder helper. It makes callers harder to reason about and signals that tuple-position APIs are carrying too much meaning.

**Improvement direction:** fix the signature immediately and strongly consider returning a named result object instead of a positional tuple.

### 10. `FragmentRouter` no longer buys enough abstraction to justify its indirection

**Severity:** Low  
**Status:** New finding

**Evidence:**
- [../../src/rigel/scan.py#L35](../../src/rigel/scan.py#L35)
- [../../src/rigel/pipeline.py#L400](../../src/rigel/pipeline.py#L400)

`FragmentRouter` still presents itself as the routing/scoring abstraction, but the live implementation is essentially a thin wrapper around the native scorer path with one call site. It also stores `strand_models`, which the class does not use in the visible path.

**Why it matters:** unnecessary wrappers make hot-path code harder to follow without adding real flexibility.

**Improvement direction:** either reduce the class to a narrowly named native scan adapter or flatten it into a small helper around the streaming scorer.

## Areas that look strong

- [../../src/rigel/config.py](../../src/rigel/config.py) has a clear, practical configuration shape.
- [../../src/rigel/buffer.py](../../src/rigel/buffer.py) is a strong example of purposeful memory-oriented design.
- [../../src/rigel/scored_fragments.py](../../src/rigel/scored_fragments.py) and [../../src/rigel/stats.py](../../src/rigel/stats.py) keep their responsibilities narrow.
- [../../src/rigel/mappability.py](../../src/rigel/mappability.py) reads as a focused orchestration module instead of a dumping ground.

## Recommended execution order

1. Fix API and persistence contracts first: findings 1, 4, and 6.
2. Finish the transcript-centric data-model cleanup next: findings 2 and 3.
3. Remove or quarantine dead/duplicated execution paths: findings 5 and 7.
4. Simplify reporting and helper abstractions last: findings 8, 9, and 10.