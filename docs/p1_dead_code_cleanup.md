# P1 Dead Code Cleanup — Multimapper C++ Migration

*Created: 2026-03-08*

## Summary

P1 moved the entire multimapper (MM) scoring path from Python to C++
inside `fused_score_buffer` in `scoring.cpp`.  The Python-side
scoring/emit/MM pipeline and the per-fragment C++ wrappers are now dead
code.  This document inventories every item to remove.

---

## 1. `src/rigel/scan.py` — `FragmentRouter`

The production path is now `scan()` → `_scan_native()` →
`nc.fused_score_buffer()`.  All Python-side scoring, emit, and MM
accumulation code is dead.

### Dead imports (top of file)

| Import | Used only by |
|--------|-------------|
| `LOG_SAFE_FLOOR` (top-level, L25) | `_gdna_log_lik` — `_scan_native` re-imports locally |
| `LOG_HALF` (L26) | `_score_wta_mrna`, `_score_wta_nrna`, `_gdna_log_lik` |
| `frag_len_log_lik` (L27) | `_score_wta_nrna`, `_gdna_log_lik` |
| `genomic_to_transcript_pos_bisect` (L28) | `_score_wta_mrna` |
| `compute_fragment_weight` (L29) | `_score_wta_mrna`, `_score_wta_nrna` |
| `_TAIL_DECAY_LP` (L33) | `_score_wta_mrna` |

### Dead state in `__init__`

| Item | Purpose |
|------|---------|
| `from array import array as _array` + 15 `array.array` accumulators (L86–103) | CSR accumulation for Python merge path |
| `self.mm_pending` / `self.mm_fid` (L106–107) | MM group Python state |

### Dead methods

| Method | Lines | Purpose |
|--------|-------|---------|
| `_score_wta_mrna` | L116–241 | Python mRNA scoring dispatcher |
| `_score_wta_nrna` | L247–365 | Python nRNA scoring dispatcher |
| `_emit_mrna` | L371–397 | CSR append for mRNA winners |
| `_emit_nrna` | L399–416 | CSR append for nRNA winners |
| `_gdna_log_lik` | L422–434 | Per-fragment gDNA score |
| `_finalize_unit_metadata` | L436–447 | Per-unit metadata append |
| `_add_single_fragment` | L453–502 | Per-fragment score+emit loop |
| `_flush_mm_group` | L510–641 | MM group merge/emit in Python |

**Total: ~530 lines to remove from scan.py**

---

## 2. `src/rigel/native/scoring.cpp`

Three standalone methods that returned Python objects per-fragment
(called from the now-dead Python methods above).

| Method | Lines | Purpose |
|--------|-------|---------|
| `score_wta_mrna` | L305–443 | Per-fragment mRNA scoring → `nb::dict` |
| `score_wta_nrna` | L449–544 | Per-fragment nRNA scoring → `nb::dict` |
| `score_emit_fragment` | L556–778 | Fused score+emit → `nb::tuple` of bytes |

Plus their nanobind `.def()` bindings (~20 lines).

**Total: ~450 lines to remove from scoring.cpp**

---

## 3. `src/rigel/scoring.py` — Test-only functions

These functions are no longer called from production code but are still
imported by test modules.  **Keep for now; mark as test-only.**

| Function | Lines | Still used by |
|----------|-------|--------------|
| `frag_len_log_lik` | L253–263 | `test_frag_length_model.py` |
| `genomic_to_transcript_pos` | L271–322 | `test_coverage_weight.py` |
| `genomic_to_transcript_pos_bisect` | L361–419 | `test_coverage_weight.py` |

---

## 4. `src/rigel/buffer.py` — Test-only code

| Item | Still used by |
|------|--------------|
| `BufferedFragment` dataclass | `test_buffer.py`, `test_pipeline_routing.py` |
| `_FinalizedChunk.__getitem__` | `test_buffer.py`, `test_pipeline_routing.py` |
| `FragmentBuffer.__iter__` | `test_buffer.py` |

**Keep for now — production never calls these but tests do.**

---

## Cleanup procedure

1. **Validate** — Profile on complex minimap2 BAM to confirm no regressions
2. **Remove scan.py dead code** — imports, `__init__` accumulators, 8 methods
3. **Remove scoring.cpp dead code** — 3 methods + 3 bindings
4. **Update docstrings** — `_scan_native` no longer mentions Python MM path
5. **Build + test** — 545+ tests must pass
6. **Re-profile** — Confirm identical output and no performance regression
