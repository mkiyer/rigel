## Scope read

I read all of: `src/rigel/buffer.py` (1–541), `scan.py` (1–186), `scan_payload.py` (1–912), `scan_cache.py` (1–596), `splice.py` (1–93), `splice_blacklist.py` (1–177), `annotate.py` (1–459), `stats.py` (1–112), `pipeline.py` (1–1135).

To check claims I also read parts of other files:
- native: `scoring.cpp` 200–279, 480–668, 673–1158; `bam_scanner.cpp` 99–124, 541–557, 860–1000, 1310–1325, 1440–1515, 1625–1660, 1755–2080, 2183–2325, 2372–2821; `resolve_context.h` 200–306, 795–830, 1230–1340; `calibration/accumulator.h` 214–226; `em_solver.cpp` 1410–1450, 1600–1680
- Python: `config.py` 160–270, 425–455, 780–835; `index.py` 425–440, 545–560, 900–967, 1030–1050, 1570–1610; `strand_model.py`; `second_pass.py` 15–40, 760–860; `frag_length_model.py` 60–130
- scripts and tests: `build_scan_cache.py`; `panel.py` 205–245; the relevant parts of `test_scan_cache.py`, `test_buffer.py`, `test_annotate.py`, `test_d7_transcript_eff_lengths.py`
- docs: the ISSUES entries listed under "Already recorded", the two TRAPS rules named there, and TESTING 100–112

Nothing was built, run or edited. "Confirmed" below means I traced the code path by reading it.

## Findings

**1. Gated fragments are stamped "resolved intergenic" in the annotated BAM.** `scan.py:140-171`, `pipeline.py:623-683`; the fallback is at `bam_scanner.cpp:2663-2667`:
```
                    stamp_and_write_hit(
                        hit_raws, out, hdr,
                        ".", ".", ".",
                        -1, -1,
                        intergenic_rna ? AF_RNA_INTERGENIC : AF_GDNA_INTERGENIC, 1.0f,
```
- **Category / severity:** defect, medium.
- **What is wrong:** the annotation table only gets rows for chimeric hits, deterministic fragments and EM units. A fragment the scorer gates gets no row (`scoring.cpp:748-751` on the deterministic path, 905-907 for EM units, 663-664 for multimapper groups). The writer treats a missing row as intergenic. So a genic fragment that the quant tables never counted (it appears only in `n_gated_out`) is written as ZF 0x25, or 0x21 if spliced, with ZW 1.0.
- **Evidence:** at the shipped alphas (0.1) this happens only when every candidate has `exon_bp ≤ 0`. At `--mismatch-alpha 0` or `--overhang-alpha 0`, which `config.py:116-121` documents as hard gates, it happens to every mismatched or overhanging fragment. A SPLICED_ANNOT read with NM > 0 then reads "RNA that no annotated transcript holds". `test_annotate.py:830-865` covers only a true intergenic sentinel.
- **Change:** have `StreamingScorer.finish` return the gated frag_ids (as it does for chimeric ones) and stamp them unresolved with a reason. Amend ZF invariant 1 (`annotate.py:119`).
- **Moves numbers:** no (BAM tags only). **Confidence:** confirmed.

**2. The scan-cache key misses two inputs of the tally: the splice blacklist and the BAM.** `scan_cache.py:35-45`, `433-491`; `index.py:937-967`:
```
THE KEY NEEDS FOUR PARTS
------------------------
* ``graph_hash`` — the region partition plus the sj CSR. The payload already carries it.
```
- **Category / severity:** contract, medium.
- **What is wrong (blacklist):** the blacklist strips junctions during the scan (`bam_scanner.cpp:541-557`). It also holds SPLICE_ARTIFACT fragments out of the accumulator (`:1642`). Both change the cached payload and strand table. But `graph_hash` hashes only the partition and the sj CSR, and `reach_digest` only the reaches. So a cache built against an index with no blacklist loads clean against one built with `--alignable-zarr`. The same holds for a different `--splice-blacklist-min-count`, or a manifest that no longer records `sources.alignable_zarr` (`index.py:1586-1604`).
- **What is wrong (BAM):** `manifest["bam"]` is written (`scan_cache.py:424`) and nothing reads it. `build_scan_cache.py:98-100` skips an existing directory unless `--force`, so a re-simulated BAM keeps its old tally.
- **Evidence:** grep finds no comparison of `provenance["bam"]` anywhere. TESTING:102-107 claims "an instrument fed a stale cache refuses it". The general rule is `TRAPS: a-hash-that-misses-its-artifact`; neither case is recorded.
- **Change:** add a blacklist digest (the rows as loaded, or "none") and a BAM identity (size + mtime, or a header + EOF digest) to the manifest, and refuse on mismatch. Have `build_scan_cache` check them before it skips a directory.
- **Moves numbers:** no. **Confidence:** confirmed.

**3. Only `scan_and_buffer` sets the resolver's transcript strands, and the comment says they are diagnostic.** `pipeline.py:281-283`:
```
    # Per-transcript strands drive the all-exonic diagnostic strand model.
    resolve_ctx = index.resolver
    resolve_ctx.set_transcript_strands(index.t_to_strand_arr.tolist())
```
- **Category / severity:** contract and stale comment, medium.
- **What is wrong:** `t_strand_arr_` also sets every fragment's `ambig_strand` (`resolve_context.h:1322-1336`). That drives FRAG_AMBIG_OPP_STRAND, the same/ambig-strand counts and strand training. It also sets the inferred strand of every gap hypothesis (`:808-811`). Nothing sets it at index load: an unconfigured resolver gives `ambig_strand` 0 everywhere. The annotation writer silently depends on this call, and on the max-fragment-length and anchor-tolerance setters next to it.
- **Evidence:** `pipeline.py:283` is the only caller in src, tests or scripts.
- **Change:** install the strands when the index builds the resolver, fix the comment, and pass the scan's resolver settings to the writer explicitly.
- **Moves numbers:** no. **Confidence:** confirmed.

**4. ZF documents transcript kinds as pools named "mRNA" and "nRNA".** `annotate.py:108-110` (and 120, 13-17, 52-59):
```
#   bit 1 (0x02)  is_mrna                 pool = mRNA
#   bit 2 (0x04)  is_gdna                 pool = gDNA (unspliced intergenic + EM-gDNA)
#   bit 3 (0x08)  is_nrna                 pool = nRNA (annotated or synthetic)
```
- **Category / severity:** naming, medium.
- **What is wrong:** invariant 2's "exactly one of {is_mrna, is_gdna, is_nrna}" is a population split that Axiom 0 forbids. The bits only record which kind of transcript won.
- **Evidence:** `is_nrna` is also true on an *annotated* single-exon transcript that covers a gene span (`index.py:550-556`; `TRAPS: nrna-does-not-mean-synthetic`). Fragments of that annotated isoform read "pool = nRNA" and fall outside a user's `ZF & 0x02` filter. `pipeline.py:1056` ("the gene's whole mature output") uses the same wording.
- **Change:** describe each bit as the winner's kind and keep the bit values, since the schema is user-facing.
- **Moves numbers:** no. **Confidence:** confirmed.

**5. `DeferredFragments.hypothesis_t` / `hypothesis_t_offsets` are dead.** `scan_payload.py:357-358`:
```
    hypothesis_t_offsets: np.ndarray  # int64[n_hypotheses + 1]
    hypothesis_t: np.ndarray  # int64[n_supporting] — which transcripts imply this path
```
- **Category / severity:** dead, low.
- **What is wrong:** the production scan fills them (`bam_scanner.cpp:1693-1699` → `accumulator.cpp:59-60`). They are validated (`scan_payload.py:453`, `474-479`), counted in the schema digest and written to every cache `.npz`. No production code reads them. `accumulator.h:220-221` says "the second pass weights hypotheses by their abundance", but `second_pass.py:20-23` says the score "needs no transcript abundance".
- **Evidence:** in src the only hits are in `scan_payload.py`. The only readers are two tests that call `supporting_t_of`: `test_gap_hypothesis_arbitration.py:559` and `test_gap_introns_are_cut.py:313`.
- **Change:** delete them in hygiene-ledger (b)'s re-cache, and fix the `accumulator.h` comment.
- **Moves numbers:** no (it re-caches every panel). **Confidence:** confirmed.

**6. `"auto"` means two things, and its resolution is written six times.** `pipeline.py:974-976` against `295`; `annotate.py:355-357`, `397-400`:
```
    if scan.sj_strand_tag == "auto":
        detected_spec = _native_detect_sj_tag(bam_path)
        scan = _replace(scan, sj_strand_tag=detected_spec)
```
- **Category / severity:** duplicate and contract, low.
- **What is wrong:** only `run_pipeline` resolves "auto". `scan_and_buffer` passes it straight through, and `parse_sj_tag_spec`'s fallthrough (`bam_scanner.cpp:105`) reads it as XS,ts. Five instruments copy the resolution: `build_scan_cache.py:41`, `prior_vs_oracle.py:359`, `_oracle_arms.py:149`, `quant_accuracy.py:126`, `profiler.py:345`. Tests that pass "auto" to `scan_and_buffer` scan as XS,ts: `test_gap_introns_are_cut.py:194`, `:405`, and `test_scanner_accumulator_integration.py:568`.
- **Also:** `write_annotated_bam` defaults `include_multimap=False` (pass 1's default is True) and `sj_strand_tag="auto"`, next to its own "Must match the value used in Pass 1". It also restates `_sj_tag_to_spec` inline.
- **Change:** resolve "auto" inside `scan_and_buffer`, delete the copies, and give the writer the resolved `BamScanConfig` with no defaults.
- **Moves numbers:** no on the shipped path; unknown for direct callers. **Confidence:** confirmed.

**7. The region ledger identity is never checked at the door, and two banks are coerced rather than asserted.** `scan_payload.py:782` and `789`:
```
        pool_lengths = np.ascontiguousarray(cal["pool_lengths"], dtype=np.int64)
        deposited_lengths = np.ascontiguousarray(cal["deposited_lengths"], dtype=np.uint32)
```
- **Category / severity:** contract, low.
- **What is wrong:** line 99 and `BANK_AXES` call `deposited == sum(region_start_count)` *the* ledger identity. `from_scan_result` checks only the length-histogram identity (799-804). The ledger identity is asserted only in tests (`test_accumulator_payload.py:402`, `test_scan_order_independence.py:533`). The two histograms are cast to a dtype, while `_bank` (859-862) argues that a widened array must be refused. Today the native types match, so the cast is a no-op that would hide a future change.
- **Change:** add the ledger check and assert both dtypes.
- **Moves numbers:** no. **Confidence:** confirmed.

**8. The zero-regions guard only postpones a crash.** `pipeline.py:305-307` (and the annotation at 248):
```
    # a genome whose references are all zero-length has no regions and nothing to deposit.
    if len(index.regions_df) > 0:
        _wire_calibration_regions(
```
- **Category / severity:** dead and stale comment, low.
- **What is wrong:** when the guard skips `set_regions`, the native result's calibration is None (`bam_scanner.cpp:2321-2322`). `from_scan_result` (722-725) then raises "set_regions was not called" after the full scan. Only an all-zero-length genome reaches this (`splice_graph.py:650-655`). The `AccumulatorPayload | None` return annotation is also false.
- **Change:** drop the guard and fix the annotation.
- **Moves numbers:** no. **Confidence:** confirmed.

**9. Three comments describe the side buffer as holding only gap-ambiguous fragments.** `scan_payload.py:600-602` (also 305; `pipeline.py:415-417`):
```
    #: The side buffer — fragments whose gap has more than one surviving explanation, held WHOLE for
    #: the second pass. NOT a loss: ``deposited + deferred + dropped_* == offered``, and this bank holds
    #: the fragments ``qc.deferred_undetermined_gap`` counts.
```
- **Category / severity:** stale comment, low.
- **What is wrong:** the bank also holds multi-placement records, and lines 813-820 require it to hold `deferred_undetermined_gap + deferred_multiple_placements`. `pipeline.py`'s "systematically the long ones" no longer describes the held population.
- **Moves numbers:** no. **Confidence:** confirmed.

**10. Literal numbers restate a memory budget and a config setting.** `buffer.py:291`, `463`; `splice_blacklist.py:59`, `126`:
```
    def __init__(self, spill_dir: Path | None, capacity: int = 2):
            self._spill_writer = _SpillWriter(self._spill_dir, capacity=2)
```
- **Category / severity:** bare constant, low.
- **What is wrong:** the spill capacity sets how many chunks stay in memory beyond `buffer_size_bytes`: a submitted chunk leaves `_memory_bytes` at line 449 but stays resident until written (350-354). It is written twice and is not in `CONSTANTS.resources`. Separately, the blacklist's `min_count=2` restates `IndexConfig.splice_blacklist_min_count` twice; production always passes the config value.
- **Change:** add a `CONSTANTS.resources` entry for the capacity, and drop the `min_count` defaults.
- **Moves numbers:** no. **Confidence:** confirmed.

**11. The FRAG_* class codes are restated across the native boundary with no parity test.** `buffer.py:39-43` against `native/constants.h:59-63`.
- **Category / severity:** duplicate, low.
- **What is wrong:** Python computes the class array; `scoring.cpp` routes on it (701-742) and the writer labels it. The SPLICE_* codes are exported (`resolve.cpp:127-131`) and pinned (`test_splice.py:40-42`). FRAG_* is neither, which contradicts hygiene-ledger (i)'s "exported, never restated".
- **Change:** export FRAG_* from native and pin it in a test.
- **Moves numbers:** no. **Confidence:** confirmed.

**12. The strand-model half of the cache is listed by hand and is outside the schema digest.** `scan_cache.py:328-330`:
```
_SJ_COLUMNS = ("ref_id", "start", "end", "motif_strand", "n_sense", "n_antisense")
_STRAND_SUBMODELS = ("exonic_spliced", "exonic")
_COUNT_FIELDS = ("pos_pos", "pos_neg", "neg_pos", "neg_neg")
```
- **Category / severity:** test gap, low.
- **What is wrong:** every `StrandModel` / `StrandModels` field has a default, so a new field would be dropped on write and silently defaulted on load. The round-trip test (`test_scan_cache.py:164-178`) restates the same lists, so it would stay green.
- **Change:** derive the lists from `dataclasses.fields` and add the strand schema to the digest.
- **Moves numbers:** no. **Confidence:** confirmed.

**13. "Second pass" and "scan" each name two different things.** `pipeline.py:943`, `1118-1121` against `413`, `998`; `scan.py:1-8`, `68`.
- **Category / severity:** naming, low.
- **What is wrong:** `before_second_pass` runs after the drain, which `pipeline.py:413` and `second_pass.py` call "the second pass". The hook actually fires just before the annotated-BAM re-read. Likewise `FragmentRouter.scan` in `rigel.scan` scores the buffer, while "scan" everywhere else in this slice means the BAM scan.
- **Change:** rename the hook `before_annotated_bam` and give the scoring module a scoring name.
- **Moves numbers:** no. **Confidence:** confirmed.

**14. `splice.py`'s docstring describes old behaviour.** `splice.py:4-6`, `15-16`:
```
Internal count arrays use 10 columns: 5 splice types × 2 strands
(sense/antisense).  Strand is an internal signal used for gDNA
estimation and is not exposed in user-facing output.
```
- **Category / severity:** stale comment, low.
- **What is wrong:** the antisense columns are only summed (`estimator.py:456`; `em_solver.cpp:1667-1669`), not used to estimate gDNA. The `(n_exon_bp, n_intron_bp)` it names do not exist anywhere.
- **Moves numbers:** no. **Confidence:** confirmed.

There is no alignment-repair (out-of-scope) code in this slice. The blacklist removes junction entries, and the scanner holds the fragment out rather than rebuilding its span (`bam_scanner.cpp:1635-1642`).

## Contracts

Assumptions marked ⚠ are not enforced anywhere.

- **`_FinalizedChunk.from_raw`:** ⚠ assumes the native arrays are C-contiguous, have the declared dtypes, and that `t_offsets` has `size+1` entries starting at 0.
- **`fragment_classes`:** ⚠ assumes `1 ≤ num_hits ≤ 65535`. `num_hits` is cast to uint16 without a check (`resolve_context.h:252`). A value of 0, or a wrap to 1 from more than 65,535 pairs, classifies a multimapper hit as unique (this is code due for rewrite). Intergenic rows are classed FRAG_UNAMBIG.
- **`FragmentBuffer`:** assumes one consumer, chunks consumed once and in order. Resident memory can exceed `max_memory_bytes` by the chunks still being spilled.
- **`FragmentRouter.scan`:**
  - ⚠ `log_every ≥ 1` is not validated; 0 raises ZeroDivisionError after the scan and calibration have run.
  - It hard-codes the native rule that deterministic fragments are SPLICED_ANNOT (`scan.py:170`).
  - It writes one row per chimeric hit. A multimapper with chimeric hits relies on the writer's last-row-wins and on EM rows being appended later.
- **`scan_and_buffer`:** assumes a name-grouped BAM (only `SO:coordinate` is refused) and a sj tag already resolved by the caller. It mutates the shared `index.resolver`.
- **`from_scan_result`:** ⚠ trusts `ref_region_bound_offsets` and `ref_sj_offsets` to be monotone (it checks only their length and end).
- **`with_drain`:** ⚠ assumes the `DrainQC` describes this payload. That holds for its only caller, `second_pass.drain`.
- **`read_scan_cache`:** ⚠ assumes the caller chose the right BAM's directory and an index with the same blacklist. Without `scan_config` the scan settings go unchecked (documented).
- **`write_annotated_bam`:**
  - ⚠ assumes the same BAM, the same filters and the same grouping as pass 1, so that frag_ids stay in step.
  - ⚠ assumes the resolver is still configured by the scan.
- **`_setup_geometry_and_estimator`:** its `n_observations ≤ 0` refusal can never fire on the production path (already pinned in `test_d7_transcript_eff_lengths.py:132-155`).
- **`aggregate_splice_blacklist`:** ⚠ assumes the store's columns exist (KeyError otherwise) and that anchors fit in int32.

## Already recorded

- **hygiene-ledger:**
  - (b) `region_end_count` / `region_span_count`; `sj_strand` / `merge_criteria` dropped at `from_raw`
  - (c) the resolver's `intron_bp`
  - (e) the buffer finalizer and writer-error tests; the one-worker exception gate; the ZS set-membership check; the report test's blacklist keys; the blacklist-loads test
  - (f) the `--no-mappability` name (`splice_blacklist.py:156`)
  - (g) the AF bits; `quant_command`'s annotated-BAM pass
  - (h) `_apply_scan_stats`' `.get(key, 0)`; the writer re-implementing `parse_bam_record` and multimap detection; FragmentScorer's copied fields
  - (i) native constants not gathered in `constants.h`
- **latent-defects:** the scan-cache items (`deposit_digest`'s fixture, `_payload_from_parts` skipping checks, in-place writes); the finalizer deadlock; `--scan-buffer-size` not bounding the native queues; the writer's `hts_close`; the census counting placements (`stats.py:90-97`); `finalize` dropping the leading 0 of `t_offsets`; StreamingScorer's missing `keep_alive`; antisense at p = ½.
- **input-mistakes-pass-silently:** a BAM with no `@HD SO`; an unpaired record.
- **sj-strand-tag-chosen-from-the-first-reads:** the XS_TS fallback for an unknown spec.
- **an-unrecorded-splice-blacklist-is-dropped-silently**
- **splicing-artifacts:** the blacklist's untested rules; the counters in three units.
- **the-em-prior-share-excludes-drawn-multimappers** (`n_calibrated`) and **the-pseudocount-strength-is-not-derived** (`em_pseudocounts`).
- **performance-memory-bounded-solve ⑥:** `add` called per fragment (`scan.py`'s loops, and with them the scalar `winner_flag` twin of `winner_flags`); `build_region_partition_arrays` and `build_sj_arrays` built repeatedly.
- **oracle-cache-key-hashes-a-thread-count** (CLOSED: `_WORK_ONLY_FIELDS`).
- **TRAPS:** `a-hash-that-misses-its-artifact` (the rule behind finding 2); `nrna-does-not-mean-synthetic` (behind finding 4).
- **CLAUDE.md / `build_scan_cache` docstring:** the cache key does not hash `resolve.cpp`.
