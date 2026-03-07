# Python ↔ C++ Redundancy Cleanup Plan

## Background

The codebase is transitioning from pure Python to C/C++ for performance-critical
paths.  The hot path is now fully C++: **BAM scanning → Fragment construction →
Resolution → Buffering** (all in `bam_scanner.cpp` via `BamScanner.scan()`).
Python orchestrates cold-path work: EM estimation, output, annotation.

The original Python modules that *were* the hot path still exist essentially
intact, creating parallel redundant implementations.  This plan identifies and
eliminates that redundancy.

**Design principle:** No Python fallback code.  C++ is required.  Python modules
serve as thin wrappers, orchestration, and output formatting only.

---

## Module-by-Module Assessment

### 1. `bam.py` — BAM parsing (annotation-only)

| Function | C++ replacement | Still needed? |
|---|---|---|
| `parse_bam_file()` | `BamScanner.scan()` (`bam_scanner.cpp`) | **Yes** — annotation 2nd pass (pysam tag writing) |
| `parse_read()` | `parse_cigar()` / `build_fragment()` (`bam_scanner.cpp`) | **Yes** — via `fragment.py` → annotation |
| `detect_sj_strand_tag()` | `detect_sj_strand_tag()` (`bam_scanner.cpp`) | **No** — `pipeline.py` uses C++ version exclusively |
| `_group_records_by_hit()` | `process_read_group()` (`bam_scanner.cpp`) | **Yes** — via `parse_bam_file` → annotation |

### 2. `categories.py` → renamed `splice.py`

Defines `SpliceType`, `SpliceStrandCol`, and derived column sets.  Despite the
vague name, it is purely splice-junction classification.

| Symbol | C++ equivalent | Redundant? |
|---|---|---|
| `SpliceType` (0/1/2) | `SPLICE_*` in `constants.h` | Values duplicated (acceptable — C++ can't import Python) |
| `SpliceStrandCol` (0–5) | Not in C++ | No — Python-only count array indexing |
| `ANTISENSE_COLS`, `SPLICED_COLS` | Not in C++ | No — Python-only |

### 3. `fragment.py` — Fragment dataclass (annotation-only)

| Symbol | C++ equivalent | Still needed? |
|---|---|---|
| `Fragment` dataclass | `AssembledFragment` (`bam_scanner.cpp`) | **Only by** `annotate.py` |
| `Fragment.from_reads()` | `build_fragment()` (`bam_scanner.cpp`) | **Only by** `annotate.py` |

### 4. `types.py` — Core types and enums

| Python type | C++ equivalent | Verdict |
|---|---|---|
| `Strand` (IntEnum) | `STRAND_*` in `constants.h` | **Keep** — has methods (`from_str`, `opposite`, etc.) |
| `IntervalType` (IntEnum) | `ITYPE_*` in `constants.h` | **Keep** — index build |
| `GenomicInterval` (NamedTuple) | `ExonBlock`/`IntronBlock` structs | **Keep** — annotation/index |
| `MergeOutcome` (IntEnum) | `MC_*` in `constants.h` | **Keep** — buffer output logic |
| `MergeResult` (dataclass) | `MergeResult` struct in `constants.h` | **Redundant** — never instantiated in production |
| `EMPTY_MERGE` | — | **Dead** with `MergeResult` |
| `ChimeraType` (IntEnum) | `CHIMERA_*` in `constants.h` | **Keep** — buffer/output |
| `Interval` (NamedTuple) | — | **Keep** — index build |
| `AnnotatedInterval` (NamedTuple) | — | **Keep** — index build |

### 5. `resolution.py` — Redundant Python implementations

| Function | C++ equivalent | Production callers | Verdict |
|---|---|---|---|
| `merge_sets_with_criteria()` | `constants.h` impl | **None** (tests only) | **Redundant** |
| `compute_frag_lengths()` | `resolve_context.h` impl | **None** (tests only) | **Redundant** |
| `resolve_fragment()` | Via `FragmentResolver` C++ | `annotate.py` | **Keep** (thin wrapper) |
| `ResolvedFragment` dataclass | `ResolvedFragment` (C++ nanobind) | **None** in production | **Redundant** |
| `pair_multimapper_reads()` | Not in C++ | `annotate.py` | **Keep** |

### 6. `scoring.py` re-exports

`scoring.py` defines `STRAND_POS = int(Strand.POS)`, `SPLICE_ANNOT = int(...)`,
etc.  These were for hot-path performance (avoid IntEnum overhead).  Now that the
hot path is C++, these re-exports add confusion.

---

## Redundancy Map

```
Python                           C++ (constants.h / resolve_context.h / bam_scanner.cpp)
─────────────────────────────    ──────────────────────────────────────────────────
types.py::Strand                 constants.h::STRAND_*                    REDUNDANT values
types.py::IntervalType           constants.h::ITYPE_*                    REDUNDANT values
types.py::MergeOutcome          constants.h::MC_*                       REDUNDANT values
types.py::MergeResult            constants.h::MergeResult struct         REDUNDANT type
types.py::ChimeraType            constants.h::CHIMERA_*                  REDUNDANT values
categories.py::SpliceType        constants.h::SPLICE_*                   REDUNDANT values
scoring.py::STRAND_POS/NEG       constants.h::STRAND_POS/NEG            TRIPLE redundancy
scoring.py::SPLICE_*/ANNOT       constants.h::SPLICE_*                   TRIPLE redundancy
resolution.py::merge_sets_*      constants.h::merge_sets_with_criteria   REDUNDANT impl
resolution.py::compute_frag_*    resolve_context.h::compute_frag_lengths REDUNDANT impl
resolution.py::ResolvedFragment  resolve.cpp::ResolvedFragment             REDUNDANT type
fragment.py::Fragment            bam_scanner.cpp::AssembledFragment          REDUNDANT type
bam.py::parse_read               bam_scanner.cpp::parse_cigar+build_frag REDUNDANT impl
bam.py::detect_sj_strand_tag     bam_scanner.cpp::detect_sj_strand_tag  REDUNDANT impl
bam.py::_group_records_by_hit    bam_scanner.cpp::process_read_group    REDUNDANT impl
```

---

## Tier 1 — Safe cleanup, no behavioral change ✅ DONE

1. **Rename `categories.py` → `splice.py`** — name matches content.
2. **Remove `detect_sj_strand_tag()` from `bam.py`** — dead code, pipeline uses
   C++ version exclusively.
3. **Remove `MergeResult` / `EMPTY_MERGE` from `types.py`** — never instantiated
   in production; only imported by `resolution.py`.
4. **Remove `merge_sets_with_criteria()` and `compute_frag_lengths()` from
   `resolution.py`** — redundant Python implementations.  Rewrite tests to
   exercise C++ via `FragmentResolver`.
5. **Remove `ResolvedFragment` from `resolution.py`** — redundant with C++
   `ResolvedFragment`.  Rewrite test mocks.
6. **Remove `scoring.py` re-exports** (`STRAND_POS`, `SPLICE_ANNOT`, etc.) —
   have consumers import directly from `types.py` / `splice.py`.

## Tier 2 — Moderate effort, annotation path refactor ✅ DONE

7. Port the annotation BAM pass to C++ (htslib write), eliminating `bam.py`,
   `fragment.py`, and `pair_multimapper_reads()` in `resolution.py`.

   **Completed:**
   - `BamAnnotationWriter` C++ class added to `bam_scanner.cpp`
   - `annotate.py` rewritten to use C++ writer (pysam dependency removed)
   - `bam.py` deleted (all BAM parsing now in C++)
   - `fragment.py` deleted (replaced by `make_fragment()` in `resolution.py`)
   - `pair_multimapper_reads()` removed from `resolution.py` (C++ handles pairing)
   - `test_bam.py` and `test_fragment.py` deleted (tested dead code)

## Tier 3 — Long-term, as C++ transition progresses ✅ DONE

8. Expose C++ enum constants via nanobind (`_resolve_impl.STRAND_POS`) to
   eliminate the "must match exactly" contract between `types.py`/`splice.py`
   and `constants.h`.

   **Completed:** All constants from `constants.h` are now importable from
   `rigel._resolve_impl`: `STRAND_*`, `SPLICE_*`, `MC_*`, `CHIMERA_*`,
   `ITYPE_*`.

9. Move `SpliceStrandCol` column indexing to C++ if the EM/quantification path moves
   there. **Deferred** — EM remains in Python; `SpliceStrandCol` stays in
   `splice.py`.
