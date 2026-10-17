## Scope read

- `src/rigel/native/bam_scanner.cpp`: lines 1–3473. The helpers in 1–1160 were read in full (parse/build/group/pair, stats, WorkerState). The slice itself, 1160–3473, was read line by line.
- Native code, partial reads: `resolve_context.h` 55–262, 1020–1060, 1200–1285, 1376–1395; `constants.h` 1–260; `thread_queue.h` (all); `calibration/accumulator.cpp` 1010–1240 (the set ctor, `set_sj`, `offer`); `accumulator.h` 97–146 plus a grep for the scratch types; `scoring.cpp` 525–545, 595–745, 1060–1100, plus a grep for the gates.
- Python, partial reads: `pipeline.py` 95–376, 486–530, 600–700, 960–1000, 1110–1145; `config.py` 160–280; `annotate.py` 1–120, 340–460; `scan.py` 95–200; `buffer.py` 84–170; `stats.py` (all); `strand_model.py` 197–230, 370–460; `cli.py` 330–600; `report/model.py` 140–260.
- Tests: `test_scanner_accumulator_integration.py` 1–560; `test_annotate.py` 820–900; `test_buffer.py` 735–775; `test_golden_output.py` 110–150; the accountability helpers in `scenarios*/conftest.py`.
- Docs: ISSUES entries hygiene-ledger, latent-defects, splicing-artifacts, sj-strand-tag-chosen-from-the-first-reads, performance-memory-bounded-solve, multimapper-intergenic-alignments and the-calibration-count-is-blind-to-multimappers; DESIGN 890–930; MANUAL 64–80; a grep of TRAPS.

## Findings

**1. The intergenic append for unique mappers does nothing, except mislabel inter-chromosomal intergenic pairs as chimeric** — `bam_scanner.cpp:1881–1897`
```
// unique-mapper zero-candidate fragments are
// appended to the buffer so that calibration's
// geometric categorization can label them INTERGENIC.
...
if (is_unique_mapper) {
    ResolvedFragment ig_result = ResolvedFragment::from_core(cr);
```
- **Category:** dead + defect + stale-comment. **Severity:** medium.
- **What is wrong:** calibration no longer reads the buffer: it receives intergenic hits through `collect_placement` (line 1902). The scorer skips every row with no candidate (`scoring.cpp:738`, `if (n_cand <= 0) continue;`). So these rows do nothing, except when the pair's blocks sit on two references:
  - `_resolve_core` sets CHIMERA_TRANS from the block references alone (`resolve_context.h:1044–1045`).
  - `fragment_classes` then marks the row FRAG_CHIMERIC (`buffer.py:163`).
  - The scorer's chimeric branch runs before any candidate check (`scoring.cpp:701–707`).
  - `scan.py:141–152` adds an AF_CHIMERIC annotation row.
  - Pass 2 then stamps ZF=0x40 on a fragment the scanner counted as `n_intergenic_unspliced` and that summary.json books as gDNA.
  - The comment's "so they don't increment stat_gated either" is also wrong for such a row: it increments `stat_chim`.
- **Evidence:**
  - The fixture `multi_reference_bam` (`test_scanner_accumulator_integration.py:181`) is exactly this pair. Run it with `annotated_bam_path` and the "trans" records should read 0x40 by this path trace (not run); the single-reference sentinel in `test_annotate.py:853` reads 0x25.
  - Nothing in `calibration/`, `second_pass.py` or `scan_payload.py` touches FragmentBuffer.
- **Proposed change:** delete the append and its comment, and add that falsification test. The repair recorded under `ISSUES: multimapper-intergenic-alignments` would buffer intergenic alignments again, but with a consumer.
- **Moves a number:** no for quantification and summary.json; yes for those records' ZF and the logged `buffer.total_fragments`.
- **Confidence:** confirmed by tracing the path; the test was not run.

**2. `proper_pair` counts hits, but summary.json and the report present it as a share of alignment records** — `bam_scanner.cpp:1797–1812`
```
// Track hit-level mate stats
for (const auto& [r1s, r2s] : hit_result.hits) {
    if (r1s.empty() || r2s.empty()) {
        stats.mate_unmapped++;
    } else {
        bool found_proper = false;
```
- **Category:** defect. **Severity:** medium.
- **What is wrong:** the unit is a hit of `group_records_by_hit`, taken before secondary pairing. That is neither a record, a read name nor a fragment.
  - `cli.py:549–555` files it under "Record-level counts", and `report/model.py:179` divides it by `total_reads`. A perfectly paired unique library therefore reads "proper pairs 50 %".
  - STAR's secondary HI groups hold no primary R1, so every secondary locus counts as improper.
  - Without HI, a multimapper's primary pair goes to `sec_r*_locs` and is counted nowhere.
- **Evidence:** a BAM with one unique pair gives total=2, proper_pair=1. A pair with STAR NH:3 gives proper 1, improper 2.
- **Proposed change:** count primary records carrying `BAM_FPROPER_PAIR` in the reader, so the number matches its label. Alternatively count per read name and move it to the read-group block. Redefine or delete `improper_pair` and `mate_unmapped` in the same unit.
- **Moves a number:** no for quantification; yes for summary.json `proper_pairs` and the report row.
- **Confidence:** confirmed.

**3. Pass 2 closes a read-name group where pass 1 does not, and nothing checks that the two stay in step** — `bam_scanner.cpp:2714–2726` against `1461–1487`
```
if (filtered) {
    // A filtered record that opens a new read name: write the held group first, ...
    if (!light_group.empty() && current_qname != bam_get_qname(b)) {
        process_and_write();
```
- **Category:** contract. **Severity:** medium.
- **What is wrong:** pass 1 drops a filtered record before it compares names, so `[A/1, X(unmapped), A/2]` is one group. Pass 2 flushes at X and opens A again, so its frag_id runs one ahead and every later read name is stamped with another fragment's annotation, silently. Only `SO:coordinate` is refused (`pipeline.py:204`), so an uncollated `SO:unsorted` BAM gets through.
- **Evidence:** that BAM followed by a genic pair B: A/2 carries B's tags and B falls to the intergenic default. `write_annotated_bam` only logs `n_read_groups`, and `pipeline.py:1125` discards the summary.
- **Proposed change:** Python raises when `summary["n_read_groups"] != stats.n_read_names`. The two are equal by construction on collated input, since both count every group including dropped multimappers.
- **Moves a number:** no.
- **Confidence:** confirmed by reading; how often real input does this is unknown and probably rare.

**4. The writer assumes "no row means intergenic" and "one row per frag_id"** — `bam_scanner.cpp:2415–2420`, `2656–2670`
```
std::unordered_map<int64_t, int64_t> frag_to_row;
...
    frag_to_row[fid_ptr[i]] = i;
```
- **Category:** contract. **Severity:** low.
- **What is wrong:**
  - (a) A fragment the scorer gates out (`++stat_gated`, so no row) is stamped AF_GDNA_INTERGENIC or AF_RNA_INTERGENIC, although transcripts hold it and summary.json books it as `n_gated_out`.
  - (b) A multimapper with any chimeric hit — including `pair_multimapper_reads`' cross pairs, which lines 1913–1918 append — gets an AF_CHIMERIC row (`scan.py:141`) and an EM row (`pipeline.py:673`). `frag_to_row` keeps the last one, which is correct only because `add_batch` runs after the router. If such a multimapper is gated out, it is stamped chimeric while the scanner counts it as not chimeric.
- **Evidence:** (a) a fixture run at overhang alpha 0. (b) duplicate frag_ids in `annotations.frag_ids[:n]` for any multimapper with a cross-reference cross-pair.
- **Proposed change:** `write()` refuses duplicate frag_ids, and Python passes the reason a fragment has no row instead of the writer inferring it.
- **Moves a number:** no; only the tags of gated-out fragments.
- **Confidence:** the mechanism is confirmed; frequency is small at the default alphas.

**5. About 20 scan counters are checked by no test** — `bam_scanner.cpp:339–375` (`merge_from`), `2041–2080` (export)
- **Category:** test-gap. **Severity:** low–medium.
- **What is wrong:** no test checks these against a real scan: `proper_pair`, `improper_pair`, `mate_unmapped`, `unique`, `multimapping`, `secondary`, `supplementary`, `duplicate`, `qc_fail`, `unmapped`, `n_with_exon`, `n_with_unannotated_sj`, `n_same_strand`, `n_ambig_strand`, `n_chimeric_trans`/`_cis_strand_same`/`_cis_strand_diff`, `n_multimapper_groups`/`_alignments`, `n_strand_skipped_*` and `n_sj_observed`. `total` is only asserted `> 0` (`test_bam_tag_parsing.py:754`). Each counter passes through three hand-written lists and a `.get(key, 0)`, so a dropped increment reads 0. Finding 2 survived this way.
- **Evidence:** an `rg` of tests/ for each name.
- **Proposed change:** one hand-built BAM with one record of each class, every counter asserted exactly, at `n_workers` 1 and 4.
- **Moves a number:** no.
- **Confidence:** confirmed that nothing reads them; that the suite stays green under such a mutation is suspected (not run).

**6. Defaults restated and never used; two contradict the config** — `bam_scanner.cpp:1305–1308`, `3312–3313`, `3330–3333`, `3418–3419`, `2869`, `3470`, `99–105`
```
int64_t chunk_size = 1000000,
int n_workers = 1,
int n_decomp_threads = 2,
int qname_batch_size = 512)
```
- **Category:** dead / duplicate / bare-constant. **Severity:** low.
- **What is wrong:** every caller (`pipeline.py:296`, `328`; `annotate.py:403`) passes every argument, so these defaults never run.
  - They restate BamScanConfig's defaults, and `include_multimap = false` contradicts `BamScanConfig.include_multimap = True`.
  - `write_annotated_bam` repeats the False and defaults `sj_strand_tag="auto"`. `parse_sj_tag_spec` has no "auto" case and silently returns XS_TS.
  - The one live value, `max_spliced_reads = 1000` in the tag detection, is written twice and lives outside `constants.h`.
- **Proposed change:** make the arguments required, throw on an unknown spec, and move the 1000 into `constants.h` (or delete it with the per-record repair).
- **Moves a number:** no.
- **Confidence:** confirmed by grepping every caller.

**7. A dead WorkerState member and stale comments**
- **Category:** dead + stale-comment. **Severity:** low.
- **What is wrong:**
  - `bam_scanner.cpp:424`, `DepositScratch deposit_scratch;`, is dead. The scan path uses `ws.offer_scratch` (line 1766), whose `.deposit` member (`accumulator.h:647`) is the scratch actually in use. The member appears only in its declaration and in two comments, and the comments at 415–423 and 2957–2958 ("The scan path uses `WorkerState::deposit_scratch`") are false.
  - Lines 2136–2139 and 2942–2945 describe "a 48 B `Region` interleaving four channel pairs". `Region` is now 8 B (`accumulator.h:97–100`, enforced by a static_assert).
  - Line 13 says "Returns everything to Python in one call"; the results actually stream through `chunk_callback`.
  - Line 2574 contains a literal `—`.
- **Proposed change:** delete the member and rewrite the comments.
- **Moves a number:** no.
- **Confidence:** confirmed.

**8. The exonic strand diagnostic ships per-fragment labels to compute four counts** — `bam_scanner.cpp:1969–1982`, `455–466`
```
strand_obs.exonic_obs.push_back(static_cast<int8_t>(result.align_strand));
strand_obs.exonic_truth.push_back(static_cast<int8_t>(t_strand));
```
- **Category:** duplicate / naming. **Severity:** low.
- **What is wrong:**
  - `StrandModel.from_labels` (`strand_model.py:197`) reduces the two vectors to a 2×2.
  - `merge_strand_obs` never frees the workers' copies, so about 4 B per unique genic fragment is held at `build_result`.
  - "truth" is the annotated strand, not truth. "exonic" includes intron-only fragments, and so does `n_with_exon`: the resolver gives transcript-only overlaps candidates (`resolve_context.h:1240`).
- **Proposed change:** four counters per worker, as the sj table already does; rename.
- **Moves a number:** no.
- **Confidence:** the mechanism is confirmed; the memory figure is derived, not measured.

**Rule checks with no finding:** this slice has no population vocabulary, no doc citations, no Greek letters in identifiers, and no alignment-repair code outside the pairing search the rewrite replaces. The ZF "mRNA/nRNA" pool names live in `annotate.py` and `constants.h`, outside this slice.

## Contracts

- **`BamScanner(...)` / `BamAnnotationWriter(...)`:** any sj spec string is accepted, and an unknown one becomes XS_TS. Only BamScanConfig validates the spec.
- **`set_regions`:** a single call is enforced. That `n_refs` and the reference order match the resolver's ids is unchecked; a single Python source builds both.
- **`set_sj`:** the ordering (after `set_regions`, before `scan`) and the CSR lengths are enforced.
- **`scan`:**
  - It is single-use, and nothing enforces that. A second call adds into `stats_` and `acc_set_`, exports the already moved-out exonic vectors, and `build_tid_mapping`'s `resize` keeps stale ids.
  - Paired-end input is enforced (it throws).
  - Collated input is only partly enforced: Python refuses `SO:coordinate` and nothing else.
  - `set_transcript_strands` must be called first; nothing enforces it, and the exonic diagnostic is silently empty otherwise.
- **Reader thread:** the filters run before grouping, and frag_id is the group index including dropped multimappers; pass 2 relies on this.
- **`process_qname_group_threaded`:**
  - "Unique" is decided twice. `is_multimap` (NH>1 or a secondary flag) drives the gate and the `unique`/`multimapping` counters. `num_hits == 1` drives the buffer class, strand training and the intergenic append.
  - They disagree if a read name has several HI groups with NH 1 and no 0x100 flag (for example, a record without HI is grouped as HI 0). This is suspected to be reachable; the rewrite should keep one predicate.
  - `num_hits` is unbounded here and `FragmentAccumulator::append` narrows it to uint16. A wrap to 1 would make a multimapper's rows look like unique fragments.
- **`collect_placement`:**
  - It receives `cr` after `from_core` has moved out `t_inds`, `frag_lengths`, `t_exon_bp` and `t_intron_bp`, so it may read only the copied fields. Nothing enforces this.
  - Its consecutive intron de-duplication needs `f.introns` sorted, which only `build_fragment`'s `std::set` guarantees.
  - The `f.exons.empty()` (line 1644) and `!any` (1726) branches are unreachable.
- **`offer_placements`:** at least one placement (guarded), and arena pointers are taken only after the arenas stop growing.
- **`build_result`:** per-reference pool sizes must be equal; enforced (it throws).
- **`BamAnnotationWriter::write`:**
  - It must see the same BAM, grouping and settings as pass 1; unchecked (finding 3).
  - The nine arrays must all have length `ann_size`; unchecked, and a mismatch reads out of bounds.
  - It assumes one row per frag_id and that no row means intergenic (finding 4).
- **`detect_sj_strand_tag_native`:** a BAM with fewer than 1,000 spliced reads is read to the end of the file without decompression threads.

## Already recorded

- **latent-defects:**
  - The census counts placements; its identity gate omits `deferred_multiple_placements` and runs on NH 1 only. The per-fragment claims in the comments at 314–336 and 1625–1628 and in `stats.py` belong to that repair, and so do `n_deposit_not_offered` and the artifact hold-out, which are also per placement.
  - `read_groups` never counts a fully unmapped pair.
  - `mapped_reads` counts QC-failed unmapped records as mapped.
  - The writer's `hts_close` is unchecked and its handles leak on a throw.
  - The native queues are unbounded by `--scan-buffer-size`.
  - An outward-facing (RF) pair deposits as an ordinary molecule.
  - `_resolve_core`'s non-empty-exons precondition is stated only in a comment.
- **sj-strand-tag-chosen-from-the-first-reads:** the detection from the first 1,000 spliced reads, and `parse_sj_tag_spec`'s XS_TS fallback. Finding 6 adds the `write_annotated_bam` default that still reaches the fallback.
- **multimapper-intergenic-alignments:** a multimapper's intergenic hits are never buffered. Its citation `bam_scanner.cpp:1808-1845` is stale; the append is now at 1892–1897.
- **splicing-artifacts:** the counters are in three units. Its statement "the census counts uniquely mapped fragments" contradicts latent-defects.
- **hygiene-ledger:**
  - (b) `region_span_count` and `region_end_count` are exported and read by nothing; `sj_strand` and `merge_criteria` are dropped.
  - (e) The worker-exception gate runs one worker; the ZS check is set membership only.
  - (g) The AF bits in `constants.h` are unused in C++.
  - (h) `write` re-implements `parse_bam_record`, multimap detection is written twice, and `_apply_scan_stats` hand-lists its keys.
  - (i) The scanner's batch sizes and C++ literals are unaudited.
- **performance-memory-bounded-solve:** `gap_introns` is copied per hit, and `--annotated-bam` adds rows one fragment at a time.
