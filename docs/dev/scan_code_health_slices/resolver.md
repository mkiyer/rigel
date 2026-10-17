## Scope read

- **In depth, all lines:** `src/rigel/native/resolve.cpp` 1–132, `src/rigel/native/resolve_context.h` 1–1437, `src/rigel/native/constants.h` 1–417.
- **Context, read in part:**
  - `bam_scanner.cpp` 375–1030 (AssembledFragment, parse_cigar, the blacklist filter, build_fragment, grouping, pairing), 1596–2030 (the per-fragment loop and collect_placement), 2088–2125, 2540–2650.
  - `index.py` 655–680, 760–830, 1060–1110, 1336–1425, 1520–1650.
  - `pipeline.py` 270–300; `config.py` 160–256; `buffer.py` 30–170; `scoring.cpp` 680–720.
  - `calibration/accumulator.h` 205–240; `accumulator.cpp` 40–195; `second_pass.py` 540–610 and 880–897; `splice_blacklist.py` 85–115.
  - Tests: `tests/_resolution_reference.py` (all), `test_implicit_splice.py` (all), `test_resolution.py` 640–763, `test_transcript_space_fl.py` 1–200, `test_buffer.py` (the strand asserts), `test_types.py`, `test_splice.py`.
  - Docs: ISSUES (hygiene-ledger, latent-defects, splicing-artifacts, performance, binary-cuts, the closed tx_frag_length entry), DESIGN §3.1b-ii, §3.1d and §3.4, and TRAPS grepped for this slice.
  - `git show 750cc8ee`.

## Findings

**1. The tolerance K is documented as a label rule, but it shapes every fragment's gap hypotheses.**
`resolve_context.h:677-679` (the same wording is at `resolve.cpp:119-123` and `resolve_context.h:492`):
```
    /// Set the splicing-anchor tolerance K (bp) used for one-sided slack
    /// in the SPLICED_IMPLICIT per-intron whole-containment discriminant.
    /// Negative K is clamped to 0; default 0 (strict containment).
```
- **Category:** stale-comment. **Severity:** medium.
- **What is wrong:** No per-intron discriminant exists any more; commit 750cc8ee replaced it with whole-path enumeration. K is read at lines 797, 807 and 813 for every non-chimeric fragment, SPLICED_ANNOT included, and by `transcript_contradicts_blocks`. The implied introns K selects reach the accumulator, so K moves calibration's L and its deferral decisions. `config.py:212-222` ("Used only by implicit-splice resolution … It tunes the resolver alone") and `pipeline.py:285-288` mislead in the same way.
- **Evidence:** the call path `enumerate_gap_hypotheses` → `cr.gap_*` → `collect_placement` (`bam_scanner.cpp:1689-1700`) → `offer`. The tests U1 and U3 run K=3 on spliced fragments and assert the hypotheses.
- **Also:** the negative-K clamp cannot be reached from the pipeline (`BamScanConfig` raises at `config.py:243-247`), while `set_max_fragment_length` throws on a negative value. The resolver's default is 0 and the config's default is 3, and no reason is given for 3.
- **Proposed change:** reword the four comments to say what K governs. Make the setter throw. State why the default is 3.
- **Moves numbers:** no. **Confidence:** confirmed.

**2. The index's resolver gets its strands and run settings only when `scan_and_buffer` runs, so the unit tests run a different resolver.**
`resolve_context.h:1325-1327` (with the same silent fallback at 809-811):
```
        for (int32_t t : cr.t_inds) {
            if (t >= 0 && t < static_cast<int32_t>(t_strand_arr_.size())) {
                int32_t s = t_strand_arr_[t];
```
- **Category:** test-gap and contract. **Severity:** medium.
- **What is wrong:** `TranscriptIndex.load` builds the resolver with no strand table, K=0 and max_fragment_length=0. Only `pipeline.py:282-292` installs these, by changing the shared `index.resolver` in place. So every test that calls `resolve_fragment` runs without strands, which means `ambig_strand` is always 0 and every hypothesis strand is NONE. With max_fragment_length 0, every pair with disjoint transcript sets is a chimera, unless the test sets the limit itself.
- **Evidence:**
  - `rg set_transcript_strands` finds only `pipeline.py:283`; nothing in `tests/` or `index.py`, and `_index_builder` sets nothing on the resolver.
  - The assertions at `test_buffer.py:144, 362, 375, 758` all expect `ambig_strand == 0` or `is_same_strand`, and they still pass if lines 1322-1336 are replaced by `cr.ambig_strand = 0`.
  - No test asserts `ambig_strand == 1` or a nonzero `n_strand_skipped_ambig_strand`.
- **Proposed change:** install the strands in `TranscriptIndex.load`, since they are index data. Require a strand array of length n_transcripts instead of falling back silently. Pass K and the length limit from the scan config instead of setting them on the shared object. Add a resolver test on an antisense overlap that asserts `ambig_strand == 1`, then break the branch and watch it fire.
- **Moves numbers:** no.
- **Confidence:** confirmed for the unit-test path. Suspected that no end-to-end golden would catch a broken branch: the fragment class it flips is read only by counters, the ZC label and strand-training admission.

**3. Out of scope: K's purpose is to correct where the aligner placed a read end.**
`resolve_context.h:693-696`:
```
    /// ⚠ Abutting is fine: an intron may START where a block ends. And the overlap must exceed the
    /// SPLICING ANCHOR TOLERANCE K before it counts -- the same slack the gap match uses, because it is
    /// the same physical thing: an aligner routinely runs a read a base or two past a donor. Without
    /// this the predicate rejects exactly the fragments K exists to accept.
```
- **Category:** out-of-scope. **Severity:** medium.
- **What is wrong:** An annotated intron that overlaps up to K aligned bases is cut whole, at lines 846-857, so those sequenced bases leave L, and the transcript is not treated as contradicted. In effect this re-aligns the splice point to the annotation.
- **Evidence:** `test_c_protrusion_within_K`: blocks (149,202) and (299,350) accept intron (200,299) at K=3.
- **Proposed change:** none yet, per the ruling. It is for the owner's scope decision.
- **Moves numbers:** yes, if removed. **Confidence:** confirmed for the mechanism; the out-of-scope classification is the owner's call.

**4. `detect_chimera` measures the span over exon-touching blocks only.**
`constants.h:405-410`:
```
    if (chimera_type == CHIMERA_CIS_STRAND_SAME) {
        int32_t lo = exons[item_idx[0]].start;
        int32_t hi = exons[item_idx[0]].end;
        for (int i = 1; i < n; i++) {
            lo = std::min(lo, exons[item_idx[i]].start);
            hi = std::max(hi, exons[item_idx[i]].end);
```
- **Category:** contract. **Severity:** low.
- **What is wrong:** `item_idx` holds only blocks whose exon transcript set is non-empty (lines 341-344). An outer block that lies in an intron or in intergenic space is therefore left out of the "outermost start to outermost end" that the comment at 402-404 and DESIGN §3.1b-ii promise. With one or fewer exon-touching components the check never runs, so nothing at this stage bounds a CHIMERA_NONE fragment's span. Such a fragment gets gap hypotheses over any distance and can be promoted to SPLICE_IMPLICIT. Only downstream limits catch it: L ≤ max_length in the accumulator and the scorer's tail decay.
- **Evidence:**
  - Example input: A=[100,200) in tA's exon, B=[400,480) in tB's exon, C=[5000,5050) intergenic (for instance a minimap2 supplementary piece, which `group_records_by_hit` groups with the primary), all on one strand, limit 1000. The span over A and B is 380, so the result is NONE; the true outermost span is 4950.
  - The Python twin computes the same subset (`_resolution_reference.py:179`), so the parity test cannot see it.
- **Proposed change:** the owner should pick the rule: either measure over all blocks, or document the subset and the unbounded span. Either way, compute the extent once in `_resolve_core`. Today it is computed separately for genomic_start (1374-1376), in `compute_frag_lengths` (968-973), in the caller's footprint and in `collect_placement`.
- **Moves numbers:** yes, if the span is taken over all blocks. **Confidence:** confirmed for the code; suspected as to which rule is intended.

**5. Fragment length: the single-block shortcut is a second rule.**
`resolve_context.h:960-965`:
```
        // Single exon block: FL is trivially block length (same for all candidates)
        if (exons.size() == 1) {
            int32_t fl = exons[0].end - exons[0].start;
            if (fl > 0) std::fill(frag_lengths.begin(), frag_lengths.end(), fl);
            return;
        }
```
- **Category:** duplicate. **Severity:** low.
- **What is wrong:** DESIGN §3.4 and `tx_frag_length` define the length as the span less every intron of the candidate that lies strictly inside it. The shortcut skips that subtraction, so it disagrees whenever the block reads through an intron. The same molecule then gets a different length depending on whether its mates overlapped.
- **Evidence:** a candidate with exons [0,100) and [180,300): one block [50,250) gives 200; blocks [50,120) and [130,250) give 120. The single-block test covers only a block inside one exon, and the enumeration gate covers two-mate fragments only.
- **Proposed change:** delete the shortcut. `tx_frag_length` already returns the span when no intron lies inside.
- **Moves numbers:** yes, for candidates that a merged block reads through, which are common on unspliced fragments over short introns. **Confidence:** confirmed.

**6. The per-hypothesis supporting-transcript lists have no production reader.**
`constants.h:236-237`:
```
    std::vector<int32_t>     gap_supporting;          // flat candidate transcript indices
    std::vector<int32_t>     gap_supporting_offsets;  // n_hypotheses + 1, always starts at 0
```
- **Category:** dead (speculative). **Severity:** low.
- **What is wrong:** These lists are built with an O(n) in-place insert (lines 888-891), copied into the scanner's arenas, and stored, merged and exported by the deferred bank as `hypothesis_t`. The only production use is the last tie-break in `compare_placements` (`accumulator.cpp:153-157`). The comment at `accumulator.h:220-221`, "the second pass weights hypotheses by their abundance", describes code that does not exist.
- **Evidence:**
  - `second_pass.py:896` hands over an empty `supporting_t_inds` ("neither does the scorer").
  - The second pass weights hypotheses by junction and boundary densities (`second_pass.py:563-585`).
  - `supporting_t_of` is called only from tests.
- **Proposed change:** either record the abundance weighting as an open item, or delete the lists end to end and correct the comment.
- **Moves numbers:** no expected, since the tie-break only orders records that are equal on every other key; the payload schema and the caches change. **Confidence:** confirmed.

**7. An unreachable branch in the junction merge.**
`resolve_context.h:1227-1232`:
```
            } else {
                std::vector<std::vector<int32_t>> all_sets;
                all_sets.insert(all_sets.end(),
                                exon_t_sets.begin(), exon_t_sets.end());
                all_sets.insert(all_sets.end(),
                                sj_t_sets_vec.begin(), sj_t_sets_vec.end());
```
- **Category:** dead. **Severity:** low.
- **What is wrong:** `sj_t_sets_vec` keeps only non-empty lookups (lines 1200-1202 pop the empty ones), and `has_annotated_sj` guarantees at least one. `merge_sets` over non-empty sets returns either an intersection or a non-empty union (`constants.h:290-323`), so `sj_merge.is_empty()` can never be true. For the same reason the junction merge criterion is only ever INTERSECTION or UNION.
- **Proposed change:** delete the branch.
- **Moves numbers:** no. **Confidence:** confirmed.

**8. Three comments describe old or wrong behaviour.**
- (a) `resolve.cpp:99`, "Resolve a Fragment object to its compatible transcript set." `t_inds` is not the compatible set: `resolve_context.h:689-691` says so, and `test_h2` keeps a contradicted transcript in it.
- (b) `resolve_context.h:84-85`, "…and whether the candidate transcripts disagreed about them." The disagreement flag (`out_ambiguous` / `implicit_ambiguous`) was deleted in 750cc8ee, the same commit that wrote this sentence.
- (c) `resolve_context.h:1138` and 1147-1150, the `n != ti` guard that supposedly filters annotated nascent-equivs. It cannot fire: `index.py:1639-1643` maps every non-synthetic parent to −1, and synthetics never enter the overlap index (`index.py:670-671`).
- **Category:** stale-comment. **Severity:** low.
- **Proposed change:** fix the three comments, drop the guard, and check the real contract in `set_nrna_parent_index` (each parent is an in-range, single-exon synthetic span).
- **Moves numbers:** no. **Confidence:** confirmed.

**9. `ResolverScratch`'s hand-written move operations miss a member.**
`resolve_context.h:368-370`, with the member declared at 339:
```
    ResolverScratch(ResolverScratch&& o) noexcept
        : t_exon_bp(std::move(o.t_exon_bp)),
          t_transcript_bp(std::move(o.t_transcript_bp)),
```
- **Category:** defect (latent). **Severity:** low.
- **What is wrong:** Both the move constructor and the move assignment (388-408) list members one by one and omit `hypothesis_path`. It is harmless today because that buffer is cleared before every use. But any member added later will be dropped silently on move, and moves happen: `set_metadata` move-assigns the scratch, and WorkerState is moved.
- **Proposed change:** hold `buf` in an RAII holder (a `unique_ptr` with `free`) and `= default` both moves.
- **Moves numbers:** no. **Confidence:** confirmed.

**10. Vocabulary: "nascent" used as an origin beside gDNA.**
`resolve_context.h:821-823`:
```
        // ⭐ The UNSPLICED hypothesis. Present when some compatible transcript implies nothing in the
        // gaps (a retained-intron isoform), and ALWAYS when no annotated sj was sequenced -- then
        // the molecule may be gDNA or nascent, and the gap is real template.
```
- **Category:** naming. **Severity:** low.
- **Also:** `constants.h:232` ("synthetic nascent-shadow candidate") and `resolve_context.h:1148` ("annotated nascent-equivs").
- **Proposed change:** "gDNA or RNA that has not spliced there". The `nrna_*` identifiers follow a repo-wide schema and are out of this slice to rename.
- **Moves numbers:** no. **Confidence:** confirmed.

**11. Undocumented magic numbers.**
`constants.h:131`, `h ^= std::hash<int32_t>()(k.start)  + 0x9e3779b9 + (h << 6) + (h >> 2);`. The constant appears five times across two copies of the same hash-combine. There is also a `-999` sentinel at `resolve_context.h:1323`.
- **Category:** bare-constant. **Severity:** low.
- **Proposed change:** one `hash_combine` helper with a one-line reason for the value (2³²/φ), and a bool flag instead of the sentinel.
- **Moves numbers:** no. These maps are only looked up, and the strand table is sorted on export (`bam_scanner.cpp:2088-2101`). **Confidence:** confirmed.

**Constants audit:** every constant in `constants.h` is used, except four flag bits that the ledger already records and the merge-criteria codes, which feed only a column the buffer drops (also recorded). `TAIL_DECAY_LP` has no stated reason for its value (recorded).

## Contracts

Assumptions marked **(unenforced)** are not checked anywhere.

- **The setup calls (`set_ref_names`, the `build_*` calls, `set_metadata`, `set_nrna_parent_index`):** assume the order used at `index.py:1538-1644` **(unenforced)**. A second `build_overlap_index` leaks `cr_`. Calling `_resolve_core` before `set_metadata(n)` indexes the scratch out of bounds.
- **`set_metadata(n)`:** n must exceed every transcript index in the interval sets, the junction map and the parent array **(unenforced)**.
- **`set_nrna_parent_index`:** each parent must be an in-range single-exon transcript; the code reads `exon_starts_[exon_offsets_[n]]` without a check **(unenforced)**.
- **`set_transcript_strands`:** optional; a short array silently means "no strand" (finding 2) **(unenforced)**.
- **`_resolve_core`:**
  - `exons` must be non-empty (recorded).
  - `cr` must be fresh, except that `n_sj_blacklisted` is an input field inside the output struct **(unenforced)**. In a reused `cr`, `t_inds` (never assigned on the intergenic branch), `t_exon_bp` / `t_intron_bp` (appended) and `merge_criteria` would carry over, even though lines 1180-1183 reset the junction key "for reuse".
  - Introns must be sorted for the first-annotated-junction key: `build_fragment` guarantees this, `resolve_fragment` does not **(unenforced)**.
  - A candidate missing from every overlap set gets exon_bp = read_length and intron_bp = 0 (lines 1316-1318); only a CIGAR with an N and no M reaches it.
- **`detect_chimera`:** per-block sets sorted (guaranteed); one reference (the caller gates on it).
- **`enumerate_gap_hypotheses` / `tx_frag_length`:** each transcript's exons sorted and disjoint; `t_idx` in range **(unenforced; the index guarantees both)**.
- **`ResolvedFragment::from_core`:** moves out `t_inds`, `exon_bp`, `intron_bp` and `frag_lengths`, so callers may read only the scalars and `gap_*` afterwards; `bam_scanner.cpp:1909→1988` relies on this.
- **`FragmentAccumulator::append`:** `exon_bp` and `read_length` are range-checked, but `num_hits` and `nm` are narrowed to uint16 unchecked; `num_hits` wraps at 65,536 hits.
- **`build_sj_blacklist_map`:** a duplicate key silently overwrites the earlier entry. Unreachable today because `splice_blacklist.py:98-104` aggregates by key, though the comment credits alignable with that.
- **`build_overlap_index`:** accepts any interval type but only EXON and TRANSCRIPT are read; the INTERGENIC rows that `index.py:771-773` emits are returned by queries and dropped.

## Already recorded

- `ISSUES: hygiene-ledger`:
  - (b) the accumulator's `sj_strand` and `merge_criteria` columns, and with them the merge-criteria codes.
  - (c) the resolver's `intron_bp`.
  - (g) the four unused flag bits, two of which also carry mature/nascent names.
  - (i) the native constants are unaudited (`DIRTY_INDICES_RESERVE`; no reason given for `TAIL_DECAY_LP`).
- `ISSUES: latent-defects`: `finalize` drops the leading 0 of `t_offsets_`; `_resolve_core`'s non-empty-exons precondition; an RF pair passes the strand test; a motif tag that contradicts the annotation (`sj_lookup_into` has no POS/NEG fallback).
- `ISSUES: performance-memory-bounded-solve`: `from_core` copies `gap_introns` per hit.
- `ISSUES: binary-cuts-on-continuous-quantities`: `TAIL_DECAY`; `max_frag_length` in three roles.
- `ISSUES: splicing-artifacts`: a rejected junction is left as a gap; artifact fragments' blocks stay split; the untested `<=` and OR rules; the disagreements between calibration and the EM.
- `ISSUES: a-start-in-an-intron-was-measured-from-the-wrong-exon` (closed).
- `TRAPS: one-reference-hides-refid-bugs`; `TRAPS: a-transcript-predicate-must-not-silently-drop-a-molecule`.
- The known scanner bug: `pair_multimapper_reads` resolving single mates, and `num_hits = max(NH, hits)`.
