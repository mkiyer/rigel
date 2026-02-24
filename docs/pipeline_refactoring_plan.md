# Pipeline Refactoring Plan

## Problem Statement

`pipeline.py` is a 2,182-line monolith that tangles together six
logically independent concerns: BAM scanning, candidate scoring,
CSR data construction, locus graph building, EM initialization
(gDNA priors, nRNA priors), and top-level orchestration.  Key pain
points:

| Problem | Evidence |
|---------|----------|
| **Monster function** | `_scan_and_build_em_data` is 620 lines with 5 nested closures that capture 15+ locals |
| **Triplicated scoring** | mRNA WTA scoring is inlined 3 times (unique path, MM mRNA pool, MM nRNA pool) |
| **Scattered constants** | 17 module-level `_UPPER_CASE` values spread across 1,600 lines, some only used by one function |
| **Dead code** | `_score_candidate`, `_score_nrna_candidate`, `_UnionFind` — never called |
| **Leaky internals** | Tests import `_`-prefixed private functions (`_score_gdna_candidate`, `_compute_nrna_init`, etc.) |
| **Re-export indirection** | `annotate.py` imports `pair_multimapper_reads` from `pipeline.py` — it lives in `resolution.py` |
| **C/C++ barrier** | Nested closures with captured state are impossible to translate to C; the hot-path scoring must be plain functions with explicit arguments |

The goal is to make the code clean, testable, self-documenting, and
trivially translatable to C/C++ module-by-module.

---

## Target Module Layout

After refactoring, the `src/hulkrna/` package gains three new modules
and `pipeline.py` shrinks to a thin orchestrator (~400 lines):

```
src/hulkrna/
├── pipeline.py          ~400 lines  (scan_and_buffer, count_from_buffer, run_pipeline)
├── scoring.py           ~250 lines  (ScoringContext + pure scoring functions)
├── locus.py             ~350 lines  (locus graph, EB gDNA priors, nRNA init)
├── scan.py              ~500 lines  (CSR builder — the "hot loop")
├── estimator.py                     (unchanged)
├── annotate.py                      (unchanged, fix import)
├── buffer.py                        (unchanged)
├── ... (others unchanged)
```

---

## Module Responsibilities

### 1. `scoring.py` — Pure Scoring Functions (~250 lines)

**Principle**: Every scoring function is a **free function with
explicit arguments** — no closures, no captured state, no side
effects.  Perfect C translation targets.

```python
# --- Constants (all scoring-related, in one place) ---
LOG_SAFE_FLOOR = 1e-10
LOG_HALF = math.log(0.5)
DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT = 0.01
DEFAULT_OVERHANG_ALPHA = 0.01
DEFAULT_MISMATCH_ALPHA = 0.1

GDNA_SPLICE_PENALTIES = {
    SpliceType.UNSPLICED: 1.0,
    SpliceType.SPLICED_UNANNOT: DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT,
    SpliceType.SPLICED_ANNOT: 0.0,
}

def overhang_alpha_to_log_penalty(alpha: float) -> float: ...


@dataclass(frozen=True)
class ScoringContext:
    """Pre-computed scoring parameters, built once per pipeline run.

    Replaces the 15+ captured locals in the old closure-based design.
    Every scoring function receives this as its first argument.
    """
    log_p_sense: float
    log_p_antisense: float
    anti_flag: bool              # True when protocol is RF
    log_ig_p: float              # intergenic strand log-prob
    ig_p: float                  # intergenic strand raw prob
    overhang_log_penalty: float
    mismatch_log_penalty: float
    gdna_splice_penalties: dict[int, float]
    # Fragment-length LUT (pre-finalized)
    fl_log_prob: np.ndarray | None
    fl_max_size: int
    fl_tail_base: float
    # Precomputed arrays from index (borrowed, never copied)
    t_strand_arr: np.ndarray     # int8[n_transcripts]
    g_strand_arr: np.ndarray     # int8[n_genes]
    t_to_g: np.ndarray           # int32[n_transcripts]
    nrna_base: int               # offset for nRNA indices in CSR

    @staticmethod
    def from_models(
        strand_models, frag_length_models, index, counter,
        overhang_log_penalty=None, mismatch_log_penalty=None,
        gdna_splice_penalties=None,
    ) -> "ScoringContext": ...


# --- Pure functions (C-translatable) ---

def frag_len_log_lik(ctx: ScoringContext, flen: int) -> float: ...

def score_mrna_candidate(
    ctx: ScoringContext,
    exon_strand: int, t_strand: int, splice_type: int,
    frag_len_ll: float, overhang_bp: int, nm: int,
) -> tuple[float, int]: ...
    """Returns (log_lik, count_col_idx)."""

def score_nrna_candidate(
    ctx: ScoringContext,
    exon_strand: int, gene_strand: int,
    frag_len_ll: float, overhang_bp: int, nm: int,
) -> float: ...

def score_gdna_candidate(
    ctx: ScoringContext,
    exon_strand: int, splice_type: int,
    genomic_footprint: int,
) -> float: ...
```

**What moves here from `pipeline.py`:**
- `_LOG_SAFE_FLOOR`, `_LOG_HALF`, all `DEFAULT_*` and `_DEFAULT_*` constants
- `_STRAND_POS/NEG`, `_SPLICE_*` — keep as int aliases of the enums
- `_GDNA_SPLICE_PENALTIES`
- `overhang_alpha_to_log_penalty()`
- `_score_candidate()`, `_score_gdna_candidate()`, `_score_nrna_candidate()`
  (rewritten as pure functions taking `ScoringContext`)
- New `ScoringContext` dataclass

**What dies:**
- `_score_candidate` (dead code) → replaced by `score_mrna_candidate`
- `_score_nrna_candidate` (dead code) → replaced by `score_nrna_candidate`

**Tests to update:**
- `test_gdna.py`: `from hulkrna.scoring import score_gdna_candidate, GDNA_SPLICE_PENALTIES`
- `test_gdna.py`: `from hulkrna.locus import compute_gdna_rate_from_strand, compute_nrna_init`

---

### 2. `locus.py` — Locus Graph & EM Initialization (~350 lines)

**Principle**: Everything between "buffer scan complete" and
"per-locus EM loop" — the graph-theoretic partitioning, Empirical
Bayes priors, and nRNA initialization.

```python
# From pipeline.py:
class UnionFind: ...           # (currently dead, keep for potential future use or delete)

def build_loci(em_data, index) -> list[Locus]: ...

def build_locus_em_data(
    locus, em_data, counter, index, mean_frag, gdna_init,
) -> LocusEMInput: ...

def compute_nrna_init(
    transcript_intronic_sense, transcript_intronic_antisense,
    transcript_spans, exonic_lengths, mean_frag, strand_models,
) -> np.ndarray: ...

def compute_gdna_rate_from_strand(sense, antisense, ss) -> float: ...

def compute_eb_gdna_priors(
    loci, em_data, counter, index, strand_models,
    *, k_locus=20.0, k_chrom=50.0,
) -> list[float]: ...

# Constants
EB_K_LOCUS = 20.0
EB_K_CHROM = 50.0
```

**What moves here from `pipeline.py`:**
- `_UnionFind` → `UnionFind` (or delete if truly dead)
- `_build_loci` → `build_loci`
- `_build_locus_em_data` → `build_locus_em_data`
- `_compute_nrna_init` → `compute_nrna_init`
- `_compute_gdna_rate_from_strand` → `compute_gdna_rate_from_strand`
- `_compute_eb_gdna_priors` → `compute_eb_gdna_priors`
- `_EB_K_LOCUS`, `_EB_K_CHROM`
- The scipy imports (`coo_matrix`, `connected_components`)

**Drop the leading underscores**: these are now public symbols of a
focused module.

---

### 3. `scan.py` — CSR Builder / Buffer Scanner (~500 lines)

**Principle**: The "hot loop" that turns a `FragmentBuffer` into
`ScanData`.  One class, no closures, explicit state.

```python
class EmDataBuilder:
    """Builds the global CSR (ScanData) from a FragmentBuffer.

    Replaces the 620-line _scan_and_build_em_data function and its
    5 nested closures.  State that was captured via closure is now
    explicit instance attributes set in __init__.

    Every method is a candidate for Cython/C acceleration.
    """

    def __init__(
        self,
        ctx: ScoringContext,
        counter: AbundanceEstimator,
        stats: PipelineStats,
        index: HulkIndex,
        annotations: AnnotationTable | None,
    ): ...

    # --- Candidate scoring (delegate to scoring.py) ---
    def add_mrna_candidates(self, bf) -> tuple[float, int]: ...
    def add_nrna_candidates(self, bf) -> None: ...
    def finalize_unit(self, bf) -> None: ...
    def flush_mm_group(self) -> None: ...

    # --- Annotation helpers ---
    def _annotate_deterministic(self, bf, chunk, i) -> None: ...
    def _annotate_chimeric(self, chunk, i) -> None: ...

    # --- Main driver ---
    def scan(self, buffer: FragmentBuffer, log_every: int) -> ScanData:
        """Single pass over buffer.  Returns packed ScanData."""
        ...
```

**Why a class?**
- The 15+ accumulation lists (offsets, t_indices, log_liks, etc.) and
  the multimapper group state (`mm_pending`, `mm_fid`) become
  instance attributes — explicit, visible, inspectable.
- Methods replace nested closures with the same logic but no closure
  capture.
- The class is straightforward to convert to a C extension: `__init__`
  maps to a C struct init, methods map to C functions taking
  `self`/struct pointer.
- Each method is independently testable.

**What moves here from `pipeline.py`:**
- The entire body of `_scan_and_build_em_data` (lines 548–1167)
- The 5 nested functions become methods

The triplicated WTA scoring in `_flush_mm_group` is de-duplicated
by extracting a shared `_score_wta_pool` method that both the
single-fragment and multimapper paths call.

---

### 4. `pipeline.py` — Thin Orchestrator (~400 lines)

After extraction, `pipeline.py` retains only:

```python
# --- Imports ---
from .scoring import ScoringContext, overhang_alpha_to_log_penalty, ...
from .scan import EmDataBuilder
from .locus import build_loci, build_locus_em_data, ...
from .annotate import AnnotationTable, write_annotated_bam

# --- PipelineResult ---
@dataclass
class PipelineResult: ...

# --- scan_and_buffer (BAM scan, ~240 lines, unchanged) ---
def scan_and_buffer(...): ...

# --- count_from_buffer (orchestrates scan → loci → EM, ~200 lines) ---
def count_from_buffer(...): ...

# --- run_pipeline (top-level entry, ~150 lines) ---
def run_pipeline(...): ...
```

`count_from_buffer` now reads like pseudocode:

```python
def count_from_buffer(buffer, index, strand_models, ...):
    ctx = ScoringContext.from_models(strand_models, frag_length_models, index, counter, ...)
    builder = EmDataBuilder(ctx, counter, stats, index, annotations)
    em_data = builder.scan(buffer, log_every)
    counter.nrna_init = compute_nrna_init(...)
    loci = build_loci(em_data, index)
    gdna_inits = compute_eb_gdna_priors(loci, em_data, counter, index, strand_models)
    for locus, gdna_init in zip(loci, gdna_inits):
        lem = build_locus_em_data(locus, em_data, counter, index, mean_frag, gdna_init)
        theta, alpha = counter.run_locus_em(lem, ...)
        gdna_count = counter.assign_locus_ambiguous(lem, theta, ...)
        ...
    return counter
```

---

## Cleanup Checklist

### Dead Code Removal
- [ ] Delete `_score_candidate` (L430–468) — fully inlined, never called
- [ ] Delete `_score_nrna_candidate` (L511–541) — fully inlined, never called
- [ ] Delete `_UnionFind` (L1176–1201) — replaced by scipy

### Fix Leaky Imports
- [ ] `annotate.py` imports `pair_multimapper_reads` from `pipeline` →
      change to import from `resolution`
- [ ] Stop exporting `_`-prefixed symbols; give public names after
      moving to their own modules

### Test Updates
- [ ] `test_gdna.py`: update imports to `scoring` and `locus`
- [ ] `test_pipeline_routing.py`: update import to `scan.EmDataBuilder`
- [ ] Script imports (`analyze_em_bottleneck.py`, etc.): update paths

### Constants Consolidation
All scoring/penalty constants move to `scoring.py`.  Locus EB
constants move to `locus.py`.  `pipeline.py` has zero module-level
constants.

### Deprecation of Test-Only Imports

After refactoring, the public API of each new module replaces the
private `_`-prefixed functions that tests were importing.  No test
should ever need to import a `_`-prefixed symbol.

---

## Execution Order

This is a mechanical refactoring — no behavioral changes, no test
modifications beyond import paths.  Execute in this order to maintain
a green test suite at every commit:

| Step | Action | Risk | Tests Stay Green? |
|------|--------|------|:-----------------:|
| **0** | Delete dead code (`_score_candidate`, `_score_nrna_candidate`, `_UnionFind`) | None — dead code | ✅ |
| **1** | Extract `scoring.py` (constants + `ScoringContext` + pure functions) | Low — functions are stateless | ✅ (re-export from pipeline) |
| **2** | Extract `locus.py` | Low — functions are self-contained | ✅ (re-export from pipeline) |
| **3** | Extract `scan.py` (`EmDataBuilder` class) | **Medium** — largest change, replaces closures | ✅ |
| **4** | Remove re-exports from `pipeline.py` | Low — update importers | ✅ |
| **5** | Fix `annotate.py` import of `pair_multimapper_reads` | Trivial | ✅ |
| **6** | Update test imports | Trivial | ✅ |
| **7** | Update script imports | Trivial | ✅ |
| **8** | De-duplicate WTA scoring in `EmDataBuilder` (`_score_wta_pool`) | Low — internal to new class | ✅ |

Each step is one commit.  Run `pytest tests/ -x -q` after every step.

---

## C/C++ Translation Pathway

After refactoring, the C migration path is:

1. **`scoring.py`** → C first (pure functions, no Python objects needed).
   `ScoringContext` maps to a C struct. `score_mrna_candidate` etc.
   become `static inline` functions in a header.

2. **`scan.py`** → Cython wrapper around a C `EmDataBuilder` struct.
   The hot loop (`scan()`) iterates buffer chunks in C, calling scoring
   functions via function pointers or inlined calls.

3. **`locus.py`** → `build_locus_em_data` is the next bottleneck after
   EM.  Its inner loops over units/transcripts translate directly to C.

4. **`estimator.py`** → The SQUAREM EM core (`run_locus_em`) is already
   numpy-vectorized and fast; C migration optional.

The key enabler is that **no function captures closure state** after
the refactoring.  Every function takes explicit arguments and returns
explicit results — pure input→output, no side-channel state.

---

## Line Count Projection

| Module | Before | After |
|--------|-------:|------:|
| `pipeline.py` | 2,182 | ~400 |
| `scoring.py` | — | ~250 |
| `locus.py` | — | ~350 |
| `scan.py` | — | ~500 |
| **Total** | 2,182 | ~1,500 |

Net reduction of ~680 lines (31%) from dead code removal,
de-duplication of WTA scoring, and elimination of excessive function
parameter documentation that repeats higher-level docstrings.
