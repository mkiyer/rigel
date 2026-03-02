

  

Research codebase for BAM parsing C++ conversion

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna)

Read [](file:///Users/mkiyer/proj/hulkrna/src)

Read [](file:///Users/mkiyer/proj/hulkrna/CMakeLists.txt)

Read [](file:///Users/mkiyer/proj/hulkrna/pyproject.toml)

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/native)

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/cgranges)

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/fragment.py)

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/types.py#1-1), lines 1 to 200

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/types.py#200-200), lines 200 to 234

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/bam.py#1-1), lines 1 to 200

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/categories.py)

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/resolution.py#1-1), lines 1 to 200

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/bam.py#200-200), lines 200 to 471

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/resolution.py#200-200), lines 200 to 500

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/buffer.py#1-1), lines 1 to 300

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/stats.py)

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/resolution.py#500-500), lines 500 to 901

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/scan.py#1-1), lines 1 to 200

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/index.py#1-1), lines 1 to 200

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/config.py)

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/index.py#200-200), lines 200 to 500

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/pipeline.py#1-1), lines 1 to 200

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/native/resolve.cpp#1-1), lines 1 to 200

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/native/scoring.cpp#1-1), lines 1 to 200

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/native/resolve.cpp#200-200), lines 200 to 400

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/pipeline.py#200-200), lines 200 to 637

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/index.py#500-500), lines 500 to 886

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/strand_model.py#1-1), lines 1 to 200

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/strand_model.py#200-200), lines 200 to 601

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/frag_length_model.py#1-1), lines 1 to 200

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/native/resolve.cpp#400-400), lines 400 to 700

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/cgranges/cgranges_bind.cpp)

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/native/resolve.cpp#700-700), lines 700 to 1000

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/buffer.py#300-300), lines 300 to 778

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/frag_length_model.py#200-200), lines 200 to 492

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/native/resolve.cpp#1000-1000), lines 1000 to 1340

Read [](file:///Users/mkiyer/Library/Application%20Support/Code/User/workspaceStorage/b0333c28ca532cc65598cf64915f0588/GitHub.copilot-chat/chat-session-resources/c3a97ad7-dcb8-4947-ab51-920e8db0698c/toolu_bdrk_01KUVHTS4T8KkhuZ1xoa8n3E__vscode-1772478769228/content.txt)

Now let me gather additional detail on the build system and the existing C++ patterns:  

Read [](file:///Users/mkiyer/proj/hulkrna/CMakeLists.txt)

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/fragment.py#1-1), lines 1 to 160

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/fragment.py#160-160), lines 160 to 167

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/native/resolve.cpp#1-1), lines 1 to 100

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/native/resolve.cpp#100-100), lines 100 to 200

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/native/resolve.cpp#1050-1050), lines 1050 to 1340

Now let me look at the resolution module's multimapper pairing and the buffer module:  

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/resolution.py#715-715), lines 715 to 901

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/stats.py#1-1), lines 1 to 65

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/stats.py#65-65), lines 65 to 161

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/types.py#1-1), lines 1 to 50

Now let me look at the config module to understand ScanConfig and how the conda environment/htslib is set up:  

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/config.py#1-1), lines 1 to 120

Read [](file:///Users/mkiyer/proj/hulkrna/mamba_env.yaml)

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/config.py#120-120), lines 120 to 187

Read [](file:///Users/mkiyer/proj/hulkrna/pyproject.toml#1-1), lines 1 to 50

Now let me check the existing strand/frag-length model training interfaces and understand what data flows back to Python: 

Research model training and buffer interfaces

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/strand_model.py#1-1), lines 1 to 150

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/frag_length_model.py#1-1), lines 1 to 150

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/buffer.py#580-580), lines 580 to 680

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/resolution.py#1-1), lines 1 to 100

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/pipeline.py#140-140), lines 140 to 250

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/pipeline.py#250-250), lines 250 to 300

Searched for text `pair_multimapper_reads` (`**/src/hulkrna/resolution.py`), 2 results

Searched for text `def resolve_fragment` (`**/src/hulkrna/resolution.py`), 1 result

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/resolution.py#748-748), lines 748 to 901

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/resolution.py#396-396), lines 396 to 500

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/strand_model.py#150-150), lines 150 to 601

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/pipeline.py#1-1), lines 1 to 140

Searched for text `intronic` (`**/src/hulkrna/pipeline.py`), 6 results

Searched for text `class FragmentLengthModels` (`**/src/hulkrna/frag_length_model.py`), 1 result

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/frag_length_model.py#367-367), lines 367 to 492

Searched for regex `class ResolvedFragment|class ResolvedResult` (`**/src/hulkrna/resolution.py`), 1 result

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/resolution.py#295-295), lines 295 to 395

Let me check one more important area — the `detect_sj_strand_tag` and the annotate second-pass dependency on pysam:  

Searched for text `pysam` (`**/src/hulkrna/annotate.py`), 4 results

Searched for regex `import pysam|from pysam` (`**/src/hulkrna/**/*.py`), 7 results

Searched for text `pysam` (`**/src/hulkrna/index.py`), 2 results

Good — pysam in index.py is only for FASTA reading during index building, not the hot path. Now I have everything needed to write the plan. Let me also quickly verify the native directory structure:  

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/native)

Read [](file:///Users/mkiyer/proj/hulkrna/src/hulkrna/cgranges)

Read [](file:///Users/mkiyer/proj/hulkrna/pyproject.toml#50-50), lines 50 to 69

Now I have a complete picture. Let me write the plan document.I don't have file editing tools available. I'll present the complete plan as a document for you to save. Here's the comprehensive plan:

---

## Plan: C++ htslib BAM Parsing — Phase F

**TL;DR**: Replace the entire Python/pysam BAM reading → Fragment construction → hit grouping pipeline with a single C++ extension (`_bam_impl`) using htslib directly. This eliminates 21.6M Python-level pysam attribute accesses, the `groupby` lambda, and Python-side Fragment construction. The C++ scanner will stream BAM records, filter/group by query name, parse CIGAR, build `ExonBlock`/`IntronBlock` vectors, group into hits, and route directly into the existing C++ `ResolveContext._resolve_core()` and `NativeAccumulator`. Model training observations will be batched and returned to Python. Projected savings: ~35s from a 133s baseline (profiled on 1.7 GB BAM, 21.6M records, 88.8% duplicates).

The key architectural decision is where to draw the C++/Python boundary. The most performant design keeps the entire inner loop (BAM read → resolve → buffer) in C++ and returns only aggregated model-training observations and statistics to Python, avoiding per-fragment Python round-trips.

**Steps**

1. **Add htslib as a build dependency**
   - Add `htslib` to mamba_env.yaml (conda-forge provides `htslib` headers and shared library)
   - Update CMakeLists.txt to `find_library(hts)` and `find_path(htslib/hts.h)` from the conda prefix (`$ENV{CONDA_PREFIX}/lib`, `$ENV{CONDA_PREFIX}/include`)
   - Update pyproject.toml cmake args if needed to pass the conda prefix include/lib paths

2. **Create the new C++ extension source file `src/hulkrna/native/bam_scanner.cpp`**
   - New nanobind module: `_bam_impl`
   - Include `<htslib/sam.h>`, `<htslib/hts.h>`, nanobind headers, and the existing cgranges headers
   - Mirror the constants already defined in native/resolve.cpp (Strand, SpliceType, ChimeraType, MergeCriteria values). Extract these into a shared header `src/hulkrna/native/constants.h` to avoid duplication.

3. **Create shared header `src/hulkrna/native/constants.h`**
   - Extract all `static constexpr` enum-mirror constants from resolve.cpp into this header
   - Also extract the `ExonBlock`, `IntronBlock`, `SJKey`, `SJKeyHash`, `CoreResult`, `MergeResult`, and `ChimeraResult` structs from resolve.cpp
   - Update resolve.cpp to `#include "constants.h"` instead of inline definitions
   - This shared header enables `bam_scanner.cpp` to call `ResolveContext::_resolve_core()` directly

4. **Extract `ResolveContext` into a shared header `src/hulkrna/native/resolve_context.h`**
   - Move the `ResolveContext` class declaration (with `_resolve_core()`, overlap index, SJ map, metadata) from resolve.cpp into a header
   - Keep resolve.cpp as the nanobind binding file that `#include`s the header
   - `bam_scanner.cpp` will `#include "resolve_context.h"` to call `_resolve_core()` directly in C++ without crossing the Python boundary
   - Similarly, move `ResolvedResult`, `NativeAccumulator` declarations into headers so `bam_scanner.cpp` can produce and accumulate resolved results directly

5. **Implement `BamScanner` C++ class in `bam_scanner.cpp`**
   - **Constructor**: accepts BAM file path (string), `ResolveContext*` (pointer to the existing C++ resolve context from `index._resolve_ctx`), and config flags (`skip_duplicates`, `include_multimap`, SJ strand tag names)
   - **`open()`**: open BAM via `hts_open()`, read header via `sam_hdr_read()`, build `ref_name → ref_id` mapping consistent with the `ResolveContext`'s `ref_to_id_` map
   - **`scan()`**: the core C++ scanning loop (replaces `parse_bam_file` + `Fragment.from_reads` + `resolve_fragment` + `buffer.append`):

     ```
     while (sam_read1(bam, hdr, record) >= 0):
       // Group by query_name (sequential — BAM is name-sorted)
       // Filter: skip QC-fail, unmapped; count duplicates
       // Parse CIGAR → ExonBlock[] + IntronBlock[]
       // R2 strand flip
       // Merge overlapping exon blocks
       // Group into hits (HI-tag or primary/secondary logic)
       // For each hit:
       //   Call ResolveContext::_resolve_core() directly
       //   Accumulate into NativeAccumulator
       //   Collect model-training observations
     ```

6. **Implement CIGAR parsing in C++** (replaces bam.py)
   - For each BAM record: `bam_get_cigar(b)`, iterate `n_cigar` ops
   - `BAM_CREF_SKIP` → intron (SJ), advance `pos` by `oplen`
   - `BAM_CMATCH | BAM_CEQUAL | BAM_CDIFF | BAM_CDEL` → advance reference within exon
   - Build `ExonBlock` vector (ref_id, start, end, strand) and `IntronBlock` vector (ref_id, start, end, sj_strand)
   - SJ strand: read `XS` tag via `bam_aux_get(b, "XS")` → `bam_aux2A()`, convert `'+'/'−'` to strand int. For `ts` tag: flip if `bam_is_rev(b)`. Check tags in priority order matching Python logic in bam.py

7. **Implement query-name grouping in C++** (replaces `groupby(bam_iter, key=lambda r: r.query_name)`)
   - Track `prev_qname` string; when `bam_get_qname(b)` differs, process the completed group
   - Use `strcmp` on the raw qname pointer — no Python string creation
   - This eliminates the 3.5s lambda overhead (21.6M calls)

8. **Implement BAM record filtering in C++** (replaces the inner loop of bam.py)
   - All flag checks via bitwise ops on `b->core.flag`: `BAM_FQCFAIL`, `BAM_FUNMAP`, `BAM_FDUP`, `BAM_FPAIRED`, `BAM_FSECONDARY`, `BAM_FSUPPLEMENTARY`, `BAM_FREAD1`, `BAM_FREAD2`, `BAM_FPROPER_PAIR`, `BAM_FMUNMAP`
   - NH tag: `bam_aux_get(b, "NH")` → `bam_aux2i()`
   - NM tag: `bam_aux_get(b, "NM")` → `bam_aux2i()` (sum across records for Fragment `nm`)
   - 88.8% of records are duplicates — short-circuit these immediately after the flag check, before any CIGAR parsing

9. **Implement hit grouping in C++** (replaces bam.py)
   - `HI` tag path: group by `bam_aux_get(b, "HI")` → `bam_aux2i()`
   - No-HI path: separate primary+supplementary from secondaries; same logic as Python but operating on C++ record structs
   - For the no-HI path, each "parsed record" needs: flag bits, ref_id, ref_start, CIGAR-derived exon/intron blocks, sj_strand, NM value. Store these in a lightweight `ParsedRecord` struct (no pysam `AlignedSegment`)

10. **Implement Fragment merging in C++** (replaces fragment.py)
    - For each hit (r1_records, r2_records): merge exon blocks by (ref_id, strand) with overlap merging, collect unique introns, sum NM tags
    - R2 strand flip: `strand ^ 3` (POS↔NEG, using the bitwise OR semantics where POS=1, NEG=2, and `opposite()` = `3 - strand` for POS/NEG)
    - Produce sorted `ExonBlock[]` and `IntronBlock[]` vectors, plus `genomic_footprint` = last exon end − first exon start
    - These feed directly into `_resolve_core()` — no Python `Fragment` object ever created

11. **Implement multimapper secondary pairing in C++** (replaces resolution.py)
    - For reads without HI tags: independently resolve each secondary R1 and R2 location via `_resolve_core()`
    - Perform the 4-tier pairing strategy (STRICT transcript-set intersection → FALLBACK same-ref proximity → CROSS-PAIR → SINGLETONS) in C++
    - This eliminates per-secondary round-trips through Python `Fragment.from_reads()` + `resolve_fragment()`
    - Uses `ref_id` and `ref_start` from the parsed records (already available in C++ structs)

12. **Implement model training observation collection**
    - Model training (strand + fragment length) currently happens in Python during the scan loop ([pipeline.py lines 207-258](src/hulkrna/pipeline.py#L207-L258))
    - **Design choice**: Collect training observations in C++ as batched arrays, return to Python for model updates
    - Define a `TrainingObservations` struct:
      - `strand_exonic_spliced`: pairs of `(exon_strand, sj_strand)` for qualified fragments
      - `strand_exonic`: pairs of `(exon_strand, gene_strand)` for unique-gene fragments
      - `strand_intergenic`: pairs of `(exon_strand, Strand::POS)` for intergenic unique-mappers
      - `frag_length_obs`: tuples of `(frag_length, splice_type_or_sentinel)` for unambiguous unique-mapper fragments
      - `frag_length_intergenic`: intergenic fragment lengths
    - The C++ scanner populates these vectors during the scan loop, accessing `t_to_g_arr` and `g_to_strand_arr` from the `ResolveContext` metadata
    - After `scan()` completes, return these as numpy arrays to Python, where `scan_and_buffer` feeds them into `StrandModels.observe()` and `FragmentLengthModels.observe()` in a batch loop

13. **Implement statistics collection in C++**
    - All BAM-level stats (`total`, `qc_fail`, `unmapped`, `duplicate`, `secondary`, `supplementary`, `n_read_names`, `unique`, `multimapping`, `proper_pair`, `improper_pair`, `mate_unmapped`) as C++ `int64_t` counters
    - Resolution-level stats (`n_fragments`, `n_chimeric*`, `n_intergenic*`, `n_with_exon`, `n_with_annotated_sj`, `n_with_unannotated_sj`, `n_unique_gene`, `n_multi_gene`, `n_multimapper_groups`, `n_multimapper_alignments`) as C++ counters
    - Model training stats (`n_strand_trained`, `n_strand_skipped_*`, `n_frag_length_unambiguous`, `n_frag_length_ambiguous`, `n_frag_length_intergenic`) as C++ counters
    - Return all as a Python dict at scan completion; `scan_and_buffer` merges into `PipelineStats`

14. **Implement SJ strand tag auto-detection in C++** (replaces bam.py)
    - Open BAM independently via htslib, scan up to 1000 spliced reads for `XS`/`ts` tags
    - Return detected tag names to Python
    - This is a one-time scan, not performance-critical, but avoids a separate pysam BAM open

15. **Add `_bam_impl` to the CMake build**
    - New `nanobind_add_module(_bam_impl NOMINSIZE STABLE_ABI ...)` target in CMakeLists.txt
    - Sources: `src/hulkrna/native/bam_scanner.cpp`, cgranges.c
    - Link against: `hts` (htslib shared library)
    - Include dirs: cgranges, native, conda prefix `include`
    - Install destination: `hulkrna`

16. **Update Python `scan_and_buffer()` to use the C++ scanner**
    - In pipeline.py: try to import `_bam_impl`; if available, use `BamScanner` instead of the Python path
    - Pass `index._resolve_ctx` to the scanner constructor
    - Call `scanner.scan()` → returns `(native_accumulator, training_obs, stats_dict)`
    - Feed `training_obs` into `StrandModels` / `FragmentLengthModels` via batch loops
    - Merge `stats_dict` into `PipelineStats`
    - Convert `native_accumulator` → `FragmentBuffer` via `finalize()`
    - Keep the Python fallback path intact for environments without htslib

17. **Update Python bam.py — preserve for fallback**
    - Keep bam.py intact as the Python fallback path
    - Add a module-level `_HAS_NATIVE_BAM` flag controlled by import availability
    - `detect_sj_strand_tag()` can optionally dispatch to C++ but the pysam path works fine (called once, not hot)

18. **Update `run_pipeline()` to pass BAM path to C++ scanner**
    - Currently pipeline.py opens the BAM with `pysam.AlignmentFile` and passes the iterator to `scan_and_buffer`
    - With the C++ scanner: pass the BAM **file path** (string) to the scanner, which opens via htslib internally
    - Keep the pysam open for the annotated BAM second pass (if enabled) — this is not performance-critical

19. **Handle reference name mapping consistency**
    - The `ResolveContext` builds a `ref_to_id_` map from the index at construction time ([resolve.cpp `build_overlap_index`](src/hulkrna/native/resolve.cpp))
    - The BAM header has its own `tid` mapping (integer reference IDs from htslib)
    - `BamScanner::open()` must build a `bam_tid → resolve_ref_id` translation table by matching `sam_hdr_tid2name(hdr, tid)` against `ResolveContext::ref_to_id_`
    - Contigs in the BAM but not in the index get `resolve_ref_id = -1` (fragments on these contigs are intergenic)

20. **Testing strategy**
    - Add `tests/test_bam_native.py` with parity tests: run both Python and C++ paths on the same BAM, assert identical `PipelineStats`, `FragmentBuffer` contents, and model training observations
    - Use existing test BAMs from scenarios_aligned
    - Add a specific test for the 88.8% duplicate short-circuit
    - Test multimapper pairing parity (HI-tag and no-HI-tag paths)
    - Test edge cases: unmapped mates, chimeric fragments, supplementary records, single-end detection error

**Verification**
- Run `pytest tests/` — all existing tests pass (Python fallback still works)
- Run `pytest tests/test_bam_native.py` — C++/Python parity on all test scenarios
- Profile on the 1.7 GB BAM: compare wall time, `PipelineStats` exact match, and final transcript counts within floating-point tolerance
- Expected wall time reduction: ~133s → ~95-100s (Phase F alone), stacking with prior Phase G savings

**Decisions**
- **Single C++ scan loop (not separate phases)**: The entire BAM→resolve→buffer pipeline runs in one C++ `scan()` call, avoiding per-fragment Python round-trips. Model training observations are batched and returned, not streamed.
- **Shared headers for resolve types**: Extract `ResolveContext`, `ResolvedResult`, `NativeAccumulator`, and constants into headers so `bam_scanner.cpp` links directly without code duplication.
- **Keep pysam as runtime dependency**: The annotated BAM second pass, index building, and simulation modules still use pysam. Removing pysam entirely is a separate future effort.
- **htslib from conda**: Use the conda-forge `htslib` package rather than vendoring. The conda environment already provides samtools (which depends on htslib).
- **Model training stays in Python**: Strand and fragment-length models have complex Python logic (Laplace smoothing, exponential tail decay, per-category routing). Moving them to C++ has low ROI since they're called O(unique_fragments) not O(BAM_records). Batch the observations instead.
- **Annotation table support**: The C++ scanner must track `frag_id` for each buffered fragment (already done by `NativeAccumulator`). The `AnnotationTable` is populated later during `count_from_buffer`, so no changes needed in the scanner.

---

Please enable file editing tools so I can save this plan to docs/bam_cpp_conversion_plan.md, or you can copy the above into that file. Let me know if you'd like any refinements before we begin implementation.