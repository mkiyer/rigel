## Scope read

- In full: `src/rigel/native/calibration/accumulator.cpp` (1–1222), `accumulator.h` (1–718), `src/rigel/second_pass.py` (1–897), `tests/native/_accumulator_reference.py` (1–1395), `src/rigel/scan_payload.py` (1–912).
- For context: `bam_scanner.cpp` 400–460, 1015–1130, 1595–1780, 2128–2327, 2925–3306; `resolve_context.h` 735–910, 1170–1240; `pipeline.py` 379–486, 965–995; `scan_cache.py` 205–270, 363–522; `strand_model.py` 276–300, 415–455; `sim/wgs_engine.py` 725–750.
- Tests: `test_accumulator_drain.py` (20–300, 560–720 and the test list), `test_second_pass_placements.py`, `test_second_pass_strand.py`, `test_accumulator_native_parity.py` 660–850, `test_gap_hypothesis_arbitration.py` (1–60, 200–240), `test_multimapper_offer.py` 195–210 and 365–395, `test_constants.py` 95–130.
- Docs: `ISSUES.md` (`hygiene-ledger`, `latent-defects`, `splicing-artifacts`, and grep hits), and a grep of `TRAPS.md`.

## Findings

**1. Supporting-transcript ids are stored, sorted on and exported, but nothing reads them.**
- Where: `accumulator.h:220-221`
  ```
  const std::int32_t* supporting_t;   // candidate transcripts implying this path -- the second pass
  std::size_t         n_supporting;   //   weights hypotheses by their abundance; pass one never reads
  ```
- Category: dead, stale-comment. Severity: medium.
- What is wrong: `hypothesis_t` and `hypothesis_t_offsets` are copied into every held record. They are part of the canonical key (`accumulator.cpp:152-158`) and are exported and validated. The second pass never reads them. Its own docstring says "the score needs no transcript abundance" (`second_pass.py:21-24`), and `_Span.supporting_t_inds` says "neither does the scorer". Three comments describe a consumer that does not exist: the one above, the spec's `GapHypothesis.supporting_t_inds`, and the spec's `DepositOutcome.DEFERRED` ("…and the transcript abundances").
- Evidence: `rg "supporting_t_of|hypothesis_t\b|hypothesis_t_offsets" src/ scripts/` finds only the producer, the exporter (`bam_scanner.cpp:1124-1125`) and the validator (`scan_payload.py:357,453,474,502`). The only readers are two tests (`test_gap_introns_are_cut.py:313`, `test_gap_hypothesis_arbitration.py:559`).
- Proposed change: delete them from the hypothesis, the bank, the key, the export and the spec, riding the re-cache of `hygiene-ledger` (b). At minimum, correct the three comments.
- Moves a number: yes. Records that differ in supporting transcripts and also in a later key field reorder, so the draw stream shifts. No deposit changes.
- Confidence: confirmed.

**2. `Accumulator::merge_from` checks the histogram width only after merging.**
- Where: `accumulator.cpp:969-971` and `994-995`
  ```
  for (std::size_t i = 0; i < deposited_lengths_.size(); ++i) {
      deposited_lengths_[i] += other.deposited_lengths_[i];
  ```
  ```
  if (pool_lengths_.size() != other.pool_lengths_.size()) {
      throw std::invalid_argument(
  ```
- Category: defect. Severity: medium (a latent crash path).
- What is wrong: the `max_length` check runs after the census, the bank, the regions, `deposited_lengths`, the boundaries and the sj have already been merged. If `other` has a smaller `max_length`, the loop reads past the end of `other.deposited_lengths_`. On any mismatch, `this` is left half-merged. The comment "Throws rather than skipping" promises a guard that the ordering defeats.
- Evidence: production workers share one `max_length` (`make_accumulator_set`), so no production path reaches it. The bound `Accumulator.merge_from` (`bam_scanner.cpp:3213`) accepts any pair. The only test call (`test_accumulator_native_parity.py:768`) merges accumulators of the same width.
- Proposed change: check `max_length_` first, with the other preconditions. Add a falsification test: a mismatch throws and leaves `this` unchanged.
- Moves a number: no.
- Confidence: confirmed by reading.

**3. The second pass restates the deposit's motif-strand rule two other ways.**
- Category: duplicate. Severity: low-medium.
- (a) An AMBIGUOUS motif reaches the sj lookup.
  - Where: `second_pass.py:278-280`
    ```
    observed = np.asarray(d.sj_strand, dtype=np.int64)[placement_of]
    implied = np.asarray(d.hypothesis_sj_strand, dtype=np.int64)[hyp_of]
    motif = np.where(observed != int(Strand.NONE), observed, implied).astype(np.int32)
    ```
  - What is wrong: the motif can be AMBIGUOUS when the mates' tags disagree, or when the supporting transcripts imply one path on both strands (`resolve_context.h:887`). `sj_edge_ids` matches a non-definite motif at any strand, so the spliced reading is scored on its sj's `sj_inv_length_sum`. The deposit trusts no splice on AMBIGUOUS (`accumulator.cpp:680-686`), and `placement_abundance` skips it (`second_pass.py:429`). The C++ comment at `accumulator.cpp:427`, "AMBIGUOUS never reaches here", is false. `test_the_batched_sj_lookup_IS_the_scalar_rule…` tests only POS, NEG and NONE.
- (b) The strand factor decides "this placement is RNA" from the motif tag.
  - Where: `second_pass.py:598-610`
    ```
    if observed_motif != int(Strand.NONE):
        strand[slot], _ = strand_terms(
    ```
  - What is wrong: the deposit and `A_p` decide it from an observed intron resolving to an annotated sj. `cr.sj_strand` combines the tags of annotated introns only (`resolve_context.h:1179-1192`). So an annotated observed splice arrives as NONE when the tag is missing or contradicts the annotation. Its ∅ hypothesis then gets gDNA's ½, while the deposit books it as spliced.
  - When it matters: on untagged libraries it is masked, because `p_r1_sense = ½` (`strand_model.py:284-295`). It is live on tagged stranded libraries, for the motif/annotation disagreement recorded in `latent-defects`. The ladder cannot see it: the simulator writes XS on every spliced read (`wgs_engine.py:740-742`), and `test_second_pass_strand.py` uses XS only.
- Proposed change: derive each placement's annotated-sj use and strand once, by the deposit's rule (nothing on AMBIGUOUS), and use it for rho, `A_p` and the strand factor. Fix the C++ comment. Add AMBIGUOUS cases to the lookup gate.
- Moves a number: yes, on rare held placements only; none on the ladder.
- Confidence: confirmed (call paths read).

**4. Three QC counters are read by nothing, and the production drain diverges from the spec's.**
- Where: `scan_payload.py:696-701`, inside `with_drain`
  ```
  deposited=self.qc.deposited + drain.deposited,
  dropped_too_long=self.qc.dropped_too_long + drain.dropped_too_long,
  dropped_empty=self.qc.dropped_empty + drain.dropped_empty,
  dropped_strand_undefined=(
      self.qc.dropped_strand_undefined + drain.dropped_strand_undefined
  ),
  ```
- Category: dead. Severity: low.
- What is wrong:
  - `unannotated_introns`, `contradictory_sj_strand` and `introns_absorbed` are counted, exported, required by `ScanQC.from_dict` and cached, but nothing in `src/` or `scripts/` reads them, and no output carries `payload.qc`.
  - The production drain drops their drain-time increments. The spec's `Accumulator.drain` re-enters `deposit` on the same tally and counts them. So `drain`'s claim that byte-identity with the spec "holds for free" is false for the qc.
  - `introns_absorbed` also counts introns that were clipped away (`accumulator.cpp:488-491`), not only ones "merged away".
- Evidence: `rg` for the three names in `src/` and `scripts/` finds only the producers and the validator; `rg "payload\.qc\." src/` finds only `with_drain`.
- Proposed change: delete the three, or report them and add the drain's increments. Add a gate comparing the production drain's qc with the spec's.
- Moves a number: no.
- Confidence: confirmed.

**5. The deposit merges the mates' conflicting introns into one: an alignment repair.**
- Where: `accumulator.cpp:456-458` (`normalise_introns`)
  ```
  /// ⚠ Overlapping AND ABUTTING introns both merge. Overlapping ones are contradictory observations of one
  /// molecule -- a real BAM produces them when the mates disagree about an acceptor. Abutting ones imply a
  /// zero-length exon, which is physically impossible, so no single molecule can legitimately use both.
  ```
- Category: out-of-scope. Severity: low.
- What is wrong: when the mates disagree, the deposit cuts the union of their introns, a path neither mate aligned, and computes L, the crossings and the sj lookup on it. That repairs contradictory alignments. `second_pass._covered_segments` (378-384) repeats it.
- Evidence: the stage disagreement this causes is already recorded (`splicing-artifacts`); only the scope classification is new.
- Proposed change: none yet, per the owner's ruling. Revisit when the pairing rewrite decides whether such a template is offered at all.
- Moves a number: no.
- Confidence: confirmed.

**6. `placement_abundance` re-derives the deposit's geometry in Python.**
- Where: `second_pass.py:447-450` (also 374-392 and 454-461)
  ```
  for seg_start, seg_end in segments:
      first = int(np.searchsorted(bounds, seg_start, side="right"))
      last = int(np.searchsorted(bounds, seg_end, side="left"))
      crossed += [boundary_base + b - 1 for b in range(first, last)]
  ```
- Category: duplicate. Severity: low.
- What is wrong: the path's segments, the boundaries it crosses and its region are Python copies of `normalise_introns`/`build_segments`, the crossing loop and `region_of_pos`. `region_of_pos` is already bound (`bam_scanner.cpp:3131`). The module's own rule says a scorer with its own definition of these is a second definition.
- Evidence: I read both sides and they are equivalent today (same sort, clip and merge; same searchsorted sides). The risk is drift.
- Proposed change: one kernel call that returns a placement's objects (sj ids, crossed boundaries, contained region).
- Moves a number: no.
- Confidence: confirmed.

**7. `chose_genomic` counts multimapper placement choices.**
- Where: `second_pass.py:824-825`
  ```
  introns = [tuple(pair) for pair in deferred.hypothesis_introns_of(h).tolist()]
  chose_genomic += not introns
  ```
- Category: naming. Severity: low.
- What is wrong: "genomic" here means "no implied intron". A multimapper whose drawn placement holds a single ∅ hypothesis, often spliced at its observed introns, therefore counts as "chose the genomic path" in `DrainQC` and in the `[SP2]` log (`pipeline.py:474-482`). That is not the gap composition the docstring describes (`scan_payload.py:192-196`).
- Proposed change: count only gap-held records, here and in the spec.
- Moves a number: no.
- Confidence: confirmed.

**8. The draw's clamp can pick a hypothesis scored zero.**
- Where: `second_pass.py:767` and `771`
  ```
  picked = np.searchsorted(cumulative, base + draws, side="right")
  return np.clip(picked, starts, ends - 1) - starts
  ```
- Category: defect. Severity: low.
- What is wrong: the clamp keeps a choice inside its own run, but not on a positive score. When a run ends in zero-score hypotheses, an overshooting uniform is clamped onto that zero.
- Evidence: a numpy check (no rigel import). Take 10⁶ records scored (0.3, 0.7), then one scored (1.0, 0.0), and draw 1 − 2⁻⁵³: the record picks local 1, whose score is 0. The rate is about half an ulp of the cumulative sum per record, under 10⁻³ expected events per library. `test_the_draw_never_leaks_into_the_NEXT_records_run` passes either way.
- Proposed change: clamp to the run's last positive-score slot.
- Moves a number: no, in practice.
- Confidence: confirmed (verified numerically).

**9. Vocabulary: "nascent" used as a population.**
- Where: `accumulator.h:210-211`
  ```
  /// the gap is real template, i.e. the molecule is gDNA -- or nascent RNA, which is the same unspliced
  /// span. That is why the accumulator needs no separate "could this be gDNA?" flag, and why the nascent
  ```
- Category: naming. Severity: low.
- What is wrong: the spec's twin already says "RNA that has not spliced there" and "a spanning shadow transcript".
- Proposed change: use the spec's words.
- Moves a number: no.
- Confidence: confirmed.

**10. Stale comments.**
- Category: stale-comment. Severity: low. Moves a number: no. Confidence: confirmed.
- In the slice:
  - `accumulator.h:278-295`: the `DeferredFragments` doc block now sits above `OfferedStorage` and says "Two nested variable-length levels … two offset arrays". There are three levels and five offset arrays.
  - `accumulator.h:296-297`: `OfferedStorage` claims "the parity surface", but only `offer` calls `offered_of`.
  - `ScoredHypothesis::absorbed` (`accumulator.h:384`) is written but never read.
  - `second_pass.py:3-5`, the score line at 25, and `HypothesisTerms` ("The three factors") omit the multimapper records and `A_p`.
  - `test_accumulator_drain.py:108-110` says qc is "excluded here by name". The assertion `differing <= {"gap_resolution"}` actually requires qc to agree, and it does.
- Adjacent, in `bam_scanner.cpp`'s bindings for these classes:
  - 2136-2139 and 2942-2945 describe a 48 B `Region`; `accumulator.h:100` asserts 8 B.
  - 3212 says "Every channel is an integer now"; five banks are float64.

## Contracts

- **`Accumulator(...)`**: bounds strictly increasing, `region_types` sized to the regions, `max_length ≥ 1`, `ref_id ≥ 0`. All checked.
- **`Accumulator::set_sj`**: only the sizes are checked. Nothing checks that offsets are monotone, that `boundary_right` is in range, or the (donor, acceptor, strand) slot order that the id-is-its-rank rule relies on.
- **`deposit`**: unchecked assumptions:
  - the introns are on this reference (`IntronBlock::ref_id` is never read);
  - implied and observed introns are disjoint (if not, they are silently merged);
  - each offer is one whole template. Single-mate fragments from the pairing bug pass unseen.
  - An empty hypothesis set is substituted with the unspliced one, while the spec raises IndexError on an explicit `()`.
- **`arbitrate` and `deposit` when every hypothesis is over the limit**: if two or more remain they are held, not rejected. `deferred_undetermined_gap` and the gap census then include molecules that no choice can deposit. This is by design (`test_every_pair_over_the_limit…`); the comment "kTooLong counts them, as it did before" (`accumulator.cpp:614-616`) is true only at the drain.
- **`AccumulatorSet::offer`**: one call per template, `n ≥ 1` checked. An out-of-range ref throws when `n == 1` but is silently counted EMPTY when `n > 1`. The identical-placement de-duplication absorbs duplicate pairings from the pairing bug; revisit it after the rewrite.
- **`AccumulatorSet(...)`**: nothing checks that offsets are monotone, and an oversized `region_types` is silently accepted.
- **`length_under`**: an out-of-range hypothesis index returns 0 instead of throwing.
- **`canonicalise`**: ties are identical records only because every stored field is in the key. That holds today.
- **`score_held_fragments` and `drain`**:
  - `region_types` and `sj` must come from the payload's index. `_Accumulators` checks only each reference's starting sj slot, not the count, and `_gather_delta` writes `block.shape[0]` rows without checking the reference's extent. In practice the guarantee is `read_scan_cache`'s graph_hash check.
  - They treat an observed motif of NONE as "not certified RNA" (Finding 3b).
- **`strand_terms`**: reads ∅ as pure gDNA (π_g = 1), while rho and f model ∅ as all unspliced traffic. The docstring's argument that unspliced RNA "cancels" holds only at π_g = 1.
- **`choose_hypotheses`**: assumes each run's scores are normalised. It guarantees the choice stays in its run, not that it lands on a positive score (Finding 8).

## Already recorded

- `hygiene-ledger` (b): `region_span_count` and `region_end_count` are dead; `accumulator.h`'s "only the channels something READS" is false.
- `hygiene-ledger` (g): the header's "reciprocal-opportunity sum"; the spec and the C++ word the mass rule differently.
- `latent-defects`:
  - `lift_choices` keys records on scalars only;
  - the leading intron (`build_segments`, `left_sj`);
  - a motif tag contradicting the annotated strand;
  - the splice census counts placements;
  - `set_sj` narrows to int32 unchecked;
  - `deposit_digest`'s `hypotheses=()` fixture;
  - `_payload_from_parts` skips `from_dict`'s checks.
- `splicing-artifacts`:
  - an all-unannotated-intron fragment is banked as UNSPLICED (`spliced = !sj_ids.empty()`);
  - mates disagreeing on an acceptor are merged into an unannotated intron;
  - two annotated junctions on opposite motif strands;
  - the stages apply the length rule differently.
- `binary-cuts-on-continuous-quantities`: the zero-density veto on held records; `max_frag_length = 1000` in three roles.
- `drain-contaminates-certified-rna`.
- `the-em-prior-share-excludes-drawn-multimappers`.
- TRAPS: `draining-breaks-the-oracle`, `an-equal-length-panel-defeats-the-lift`, `a-cancellation-is-conditional-on-its-support`, `a-purity-filter-is-a-length-filter`.
- The three pass-2 items you said were already planned: the length law chosen from implied introns, the empty contested set giving density 0, and NM never read.
