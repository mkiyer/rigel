## Scope read

- `src/rigel/native/bam_scanner.cpp`: 1–1199 in full. To follow calls: 1300–2029 (scan, the worker, `collect_placement`), 2286–2320 (`build_result` qc and calibration), 2370–2808 (`BamAnnotationWriter`), 2860–2920 (tag auto-detect), 3090–3314 (the Accumulator and AccumulatorSet bindings).
- `resolve_context.h`: 40–349, 417–444, 470–500, 590–634, 985–1018, 1015–1389.
- `constants.h`: 53–56, 94–235, 256–265, 336–410.
- `scoring.cpp`: 480–560, 690–900.
- `calibration/accumulator.cpp`: 560–585, 905–925, 1130–1160. `accumulator.h`: 269–300, 402–421.
- Python: `pipeline.py` 98–103, 185–238, 965–995, 1115–1135; `annotate.py` 349–420; `cli.py` 336–364; `buffer.py` 135–165; `strand_model.py` 197–227, 395–440; `config.py` 240–255; `index.py` 1620–1645; `second_pass.py` 530–575, 780–830, 875–905.
- Tests: `test_bam_tag_parsing.py` 572–632 and 848–930, `test_annotate.py` 600–700, `calibration/test_accumulator_span_unbiased.py` 180–262, `test_scanner_accumulator_integration.py` 375–465 and 560–580. I also grepped tests/ and `src/rigel/sim/` for CIGAR fixtures.
- Docs: ISSUES.md 20–284, 571–635, 767–881, 1084–1128, 2600–2612, 3144–3184; TRAPS.md 825–925.

## Findings

**1. A BAM without NM is scored as if no read had a mismatch, with no warning**
- **Where:** `bam_scanner.cpp:636–638` (the writer's copy is at 2756–2758).
  ```
  rec.nm = 0;
  uint8_t* nm_aux = bam_aux_get(b, "NM");
  if (nm_aux) rec.nm = bam_aux2i(nm_aux);
  ```
- **Category / severity:** contract · medium.
- **What is wrong:** STAR's default `--outSAMattributes Standard` writes NH HI AS nM and no NM. On such a BAM every record reads NM = 0, so two things switch off without a word:
  - the per-hit mismatch weight on multimappers (`scoring.cpp:499–500`, 521–522);
  - the hard gate at `mismatch_alpha = 0` (`scoring.cpp:748`).
  
  For a unique mapper the term is shared by every candidate and by the gDNA term, so it cancels there.
- **Evidence:** nothing counts or reports records without NM. A grep of `pipeline.py`, `cli.py` and the scanner finds no check. The MANUAL's BAM requirements and its STAR recipe never mention NM. The aligned scenarios use minimap2, which writes NM, so no test sees a BAM without it.
- **Proposed change:** add a counter `n_records_without_nm` to `BamScanStats`, export it to summary.json, and warn when it is above 0 (this needs no threshold). Add `NM` to the MANUAL's requirements and to the STAR line (`--outSAMattributes NH HI AS nM NM`).
- **Moves a number:** no.
- **Confidence:** confirmed for the code path. The STAR default comes from STAR's own manual.

**2. Pairing secondaries without HI accepts every mate combination inside one gene, and the cross-pair comment is false** (code due for rewrite: defects only)
- **Where:** `bam_scanner.cpp:938–951` and `985–991`.
  ```
  if (has_intersection(r1_resolved[i].t_inds,
                       r2_resolved[j].t_inds)) {
      paired.push_back({*r1_resolved[i].reads, *r2_resolved[j].reads});
  ...
  // combinations (e.g. R1 chr12 × R2 chr12) are valid PE hits
  ```
- **Category / severity:** defect · medium.
- **What is wrong, part (a): the STRICT step pairs by transcript overlap, and any two locations in one gene overlap.** Each mate is resolved alone, and its `t_inds` include:
  - every transcript whose body covers the location (`resolve_context.h:1250–1266`);
  - the synthetic nRNA span, which is added to every block (`resolve_context.h:1133–1172`; `index.py:1633–1644` sets it for every multi-exon transcript).
  
  So any two mate locations inside one gene intersect, and STRICT emits all n1 × n2 combinations, not only the true pairs. Chimera detection does not catch the extra ones. `detect_chimera` runs before the nRNA spans are added and ignores blocks with no exon hit (`constants.h:341–347`). An intronic cross-combination is therefore non-chimeric at any distance, and it is buffered and offered as a placement.
- **What is wrong, part (b): the CROSS-PAIR step can only make trans chimeras, not what its comment says.** The FALLBACK greedy pass looks at every unpaired same-reference pair with `ref_id >= 0`, and paired status only grows. So no unpaired same-reference pair is left after it. CROSS-PAIR therefore yields only:
  - pairs on two different references (trans chimeras, which are never deposited and are skipped by the scorer, `scoring.cpp:701`);
  - pairs of mates on references outside the index (see finding 3).
  
  The writer's comment at 2557–2561 ("overlapping hits for EM scoring") repeats the same wrong claim.
- **Evidence:** the call path `resolve_location` → `_resolve_core`, read above. An input that reaches it: a BAM without HI (minimap2) where a gDNA fragment's mates each have two or more locations in one gene's introns, for example intronic Alus.
- **Proposed change:** the rewrite's bounded R1 × R2 search must check genomic geometry, never `t_inds`. Its pairing fixture (A, B / A′, B′) should put both copies inside one gene, where the current rule over-pairs.
- **Moves a number:** yes, on input without HI that has multimappers.
- **Confidence:** confirmed by reading. How often it happens is unmeasured.

**3. Mates on references outside the index are all filed under one reference, −1**
- **Where:** `bam_scanner.cpp:690–693`, with line 658 (`rec.ref_id = mapped_ref_id;`, which is −1 for an unindexed reference).
  ```
  auto make_key = [](int32_t ref_id, int32_t strand) -> int64_t {
      return (static_cast<int64_t>(ref_id) << 32) |
             (static_cast<uint32_t>(strand));
  };
  ```
- **Category / severity:** defect · low.
- **What is wrong:** `_check_bam_header` explicitly allows references beyond the index (`pipeline.py:230–237`), so −1 reaches `build_fragment`. Two problems follow:
  - `(-1) << 32` is undefined behaviour in C++17, the standard nanobind builds with. In practice it decodes back to −1.
  - Mates on two different unindexed contigs share key −1. Their blocks are merged as if they sat on one reference: the trans pair is not flagged as a chimera, and its footprint spans two coordinate systems.
- **Evidence:** a pair with R1 on decoy A and R2 on decoy B. The resolver sees one reference, so the pair is not CHIMERA_TRANS, and it goes into the buffer as an intergenic row with chimera_type NONE.
- **Proposed change:** keep the BAM tid for unindexed references (key on the tid, or give each unindexed tid its own id).
- **Moves a number:** no quantification number moves; both cases already count as intergenic. The buffer row's class and the annotated BAM's label (gDNA-intergenic against chimeric) do change.
- **Confidence:** confirmed.

**4. No test fixture exercises the CIGAR ops D, I, =, X or H**
- **Where:** `bam_scanner.cpp:65–68`.
  ```
  return op == CIG_MATCH || op == CIG_DEL ||
         op == CIG_EQUAL || op == CIG_DIFF;
  ```
- **Category / severity:** test-gap · low.
- **What is wrong:** the code is correct against the SAM spec, but a regression here would shift every block and junction coordinate.
- **Evidence:**
  - A grep of tests/ for D/I/=/X/H CIGAR strings, or for pysam's `CDEL`/`CINS`/`CEQUAL`, finds only the `"10S"` and `"5S40N5S"` fixtures.
  - The oracle BAM writer emits only M and N (`sim/bam.py:80–92`).
  - The simulator injects substitutions only (`sim/reads.py:28`), so minimap2's alignments of its reads carry almost no indels.
  - So dropping `CIG_DEL`, or letting an insertion advance the position, most likely leaves the suite green. Not run.
- **Proposed change:** a small scan fixture whose records carry D, I, S, H, =/X before and after an N, asserting the exon blocks, the junction and both anchors.
- **Moves a number:** no.
- **Confidence:** confirmed that no fixture exists; suspected that the suite stays green.

**5. `MarshalledFragment` reads its hypotheses twice and silently empties a one-shot iterable**
- **Where:** `bam_scanner.cpp:1039–1042`.
  ```
  // ⛔ Reserved up front: the GapHypothesis spans below point INTO these vectors, so a reallocation
  // while filling them would dangle every pointer already handed out.
  std::size_t n_implied = 0, n_supporting = 0;
  for (nb::handle h : hypotheses) {
  ```
- **Category / severity:** contract and stale-comment · low.
- **What is wrong:** the spans hold offsets while the vectors fill, and are turned into pointers only afterwards (1060–1066). So the reservation pass guards against a hazard the offsets already remove, and the comment's "point INTO" is untrue. The pass also iterates `hypotheses` twice. Given a generator, the second pass sees nothing:
  - `deposit` then answers as if the fragment were unspliced only (`accumulator.cpp:571–576`);
  - `length_under` returns 0 (`accumulator.cpp:915`).
  
  Separately, `offered` points into the struct's own vectors. A move keeps those buffers, but a copy would leave dangling pointers, and nothing deletes the copy constructor.
- **Evidence:** today's callers pass lists (`second_pass.py:546`) and tuples (the tests), so it is latent.
- **Proposed change:** drop the pre-pass, or take `nb::list(hypotheses)` once; delete the copy constructor and copy assignment.
- **Moves a number:** no.
- **Confidence:** confirmed.

**6. The exonic strand diagnostic stores two labels per fragment where a 2×2 count would do**
- **Where:** `bam_scanner.cpp:249–250` and `455–465`.
  ```
  std::vector<int8_t> exonic_obs;
  std::vector<int8_t> exonic_truth;
  ```
- **Category / severity:** duplicate (a redundant representation) · low.
- **What is wrong:** `StrandModel.from_labels` (`strand_model.py:197–208`) only reduces these arrays to a 2×2. The spliced model already dropped per-fragment labels for exactly this reason (the comment at 1950–1951). Memory grows with depth: about 2 bytes per qualifying fragment in each worker, plus the merged copy. `merge_strand_obs` clears the source's sj table but not these vectors, and the workers live until `scan()` returns.
- **Proposed change:** accumulate the four counts in C++ and export them.
- **Moves a number:** no; the counts are exact.
- **Confidence:** confirmed.

**7. The nine qc keys are hand-written three times**
- **Where:** the helpers' claim at `bam_scanner.cpp:1083–1086` ("key strings live here and nowhere else"), against the qc dicts at 2298–2308, 3102–3112 and 3289–3299.
- **Category / severity:** duplicate · low.
- **What is wrong:** a field added to `DepositCounters` can reach one export and miss another. All three copies agree today.
- **Proposed change:** one `qc_dict(const DepositCounters&)` helper beside `gap_census_dict`.
- **Moves a number:** no.
- **Confidence:** confirmed.

**8. `ParsedAlignment::sj_strand` is written and never read, and the comment on `n_sj_blacklisted` names the wrong consumer**
- **Where:** `bam_scanner.cpp:192–200`, with 653–656 and the writer's 2772–2775.
  ```
  // blacklist.  Populated by filter_blacklisted_sjs in Pass 1 (for
  // scanner stats parity) and Pass 2 (for the ZB:i annotated-BAM tag).
  int32_t sj_strand;      // from XS/ts tag, STRAND_NONE if absent
  ```
- **Category / severity:** dead and stale-comment · low.
- **What is wrong:**
  - `sj_strand` is only copied into each `SJCigarEntry::strand`, and that is the value `build_fragment` reads.
  - `n_sj_blacklisted`'s pass-1 consumer is the SPLICE_ARTIFACT promotion (1859–1866 → `resolve_context.h:1288–1292`). The stat it names is accumulated through `parse_bam_record`'s pointer argument instead.
- **Evidence:** grep finds no read of `rec.sj_strand` or `->sj_strand` on a `ParsedAlignment`.
- **Proposed change:** make `sj_strand` a local variable; reword the comment.
- **Moves a number:** no.
- **Confidence:** confirmed.

**9. Two comments describe older behaviour**
- **Where and what is wrong:**
  - Line 13, "8. Returns everything to Python in one call": the scan streams chunks through `chunk_callback` (1303–1305, 1551–1556), and the header's step list leaves out the calibration offer.
  - Line 331, the identity `… == qc.deposited + Σ qc.dropped_* + n_deposit_not_offered`, omits `deferred_undetermined_gap`. The gate adds it (`test_scanner_accumulator_integration.py:412`), and `accumulator.h:412–414` states the full form. `stats.py:94` and the gate's docstring carry the same stale form.
- **Category / severity:** stale-comment · low.
- **Proposed change:** correct the text.
- **Moves a number:** no.
- **Confidence:** confirmed.

## Contracts

- **`parse_cigar`:** assumes positions fit in int32 and that every N op has reference-consuming ops on both sides. A leading or trailing N gives a junction with a 0-bp anchor, which certifies RNA unless the blacklist rejects it. Unchecked.
- **`read_sj_strand`:** assumes a present XS or ts carries `+`/`-`. A present but unusable value (`XS:A:.`, an `XS:i` score) returns STRAND_NONE (0, which is ≥ 0), so in XS_TS or TS_XS mode it blocks the second tag; only a missing tag falls through. It also assumes ts is relative to the read (it is flipped on BAM_FREVERSE). Unchecked.
- **`parse_bam_record`:**
  - NM is assumed present: an absent tag reads as 0 (finding 1).
  - NH and HI are assumed integer-typed: `bam_aux2i` returns 0 for other types, so an `NH:Z` record reads as unique.
  - References are assumed to be in the index: an unknown one maps to −1, which the header check warns about.
- **`filter_blacklisted_sjs`:** assumes blacklist coordinates are 0-based, half-open intron coordinates on the resolver's ref ids. That convention is pinned only by fixtures written in the same convention.
- **`build_fragment`:**
  - assumes an FR library with exactly one of READ1/READ2 per record;
  - assumes the mates agree wherever they overlap (a mate aligned straight through the other's intron yields an exon block that contains an observed intron);
  - assumes `ref_id ≥ 0` (finding 3);
  - assumes supplementary pieces lie in line with their primary;
  - `genomic_footprint` mixes references on a trans chimera (unread, since the scorer skips chimeras).
  
  None of this is checked.
- **`group_records_by_hit`:** a record carrying neither READ1 nor READ2 is dropped on the HI path but filed as R2 on the path without HI. A missing HI joins group 0. This is in the code due for rewrite.
- **`pair_multimapper_reads`:** every worker shares one non-const `FragmentResolver&`. Thread safety rests on `_resolve_core(…, scratch)` writing only to the scratch, which no `const` enforces. FALLBACK has no distance bound.
- **The hit loop (consumer of the slice):** a hit with no exon blocks (every record's CIGAR is `*`) is skipped with no counter, so it appears only in `n_fragments` and reads as "genic" (`cli.py:567`). `num_hits` is narrowed to uint16 unchecked (`resolve_context.h:252`); only "> 1" is read, so a count of 65,536 or more cross pairs could read as unique, which is unreachable at aligner defaults.
- **`BamScanner::scan`:** single use. `stats_` and `strand_obs_` accumulate, and `build_tid_mapping`'s `resize` keeps stale entries; nothing refuses a second call (`pipeline.py` makes one call per scanner).
- **`gap_census_dict` / `deferred_dict`:** their keys must equal the reference's `GapResolution` values and the `deferred_arrays` keys. Enforced: the parity gates compare them.

## Already recorded

- **ISSUES: hygiene-ledger**
  - (h): `BamAnnotationWriter::write` re-implements `parse_bam_record`, and multimap detection is written twice.
  - (i): the scanner's batch sizes are kept outside `constants.h`, and inline literals are unaudited (`reserve_chunk`'s 3/2).
  - (e): the one-worker exception gate, and the ZS set-membership check.
- **ISSUES: splicing-artifacts:**
  - the blacklist's `<=`, the OR of the two anchors and the inert stored 0;
  - anchors cumulative across earlier introns (`parse_cigar`'s back-fill);
  - rejection per record before the mates are joined, so one surviving junction shields the rejected ones;
  - a rejected fragment's blocks stay split at the rejected N;
  - mates disagreeing on an acceptor;
  - the counters in three units (`n_sj_observed` / `n_sj_blacklisted` over every record);
  - Rigel reads NM only.
- **ISSUES: sj-strand-tag-chosen-from-the-first-reads:** `parse_sj_tag_spec`'s fallback to XS_TS. Two additions:
  - The native fallback is still reachable from `write_annotated_bam`'s default `sj_strand_tag="auto"` (`annotate.py:357`); its only caller passes the resolved tag.
  - The members' `= SJTagMode::XS_ONLY` defaults (1136, 2375) are always overwritten by the constructor.
- **ISSUES: latent-defects:**
  - a motif tag that contradicts the annotated strand;
  - the census counting placements, and its gate omitting `deferred_multiple_placements`;
  - an outward-facing (RF) pair depositing as an ordinary molecule;
  - `_resolve_core`'s non-empty-exons precondition.
- **ISSUES: input-mistakes-pass-silently:** the unpaired record's RuntimeError; a coordinate-sorted BAM without `@HD SO`.
- **ISSUES: multimapper-intergenic-alignments:** the EM never sees a multimapper's intergenic hits (comment at 1899–1901).
- **ISSUES: performance-memory-bounded-solve:** `gap_introns` copied per hit.
- **TRAPS: a-transcript-predicate-must-not-silently-drop-a-molecule:** the same family as finding 2.
- **TRAPS: credit-exactly-one-sj:** the `sj_strand_table` rule.
