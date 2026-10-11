# The multimapper second pass — design (2026-10-10, sandbox, not authoritative)

**Status: BUILT in the worktree `~/proj/rigel-multimapper` (branch `multimapper-second-pass`), uncommitted; see §13.** The owner's framing (2026-10-10): the accumulator deposits uniquely aligned
fragments for calibration and discards multimapping ones; it already buffers the fragments whose unsequenced
gap admits several fragment lengths and assigns them in a second pass; multimapping fragments are to be
handled the same way — buffered, then assigned to the probabilistically best placement by the abundance of
uniquely aligned fragments — reusing that infrastructure. The added complexity: a multimapping fragment is ONE
fragment with several disparate alignments, each of which may be compatible with several isoforms or with
gDNA; every alignment and every assignment traces to that one fragment. The problem it fixes is
`ISSUES: the-calibration-count-is-blind-to-multimappers`; the design below is written to be built in the
order DERIVE → spec → native prototype in a worktree → A/B → `src/`, each stage with its named gate.

Sections: §1 the problem and its size · §2 the machinery as it stands, precisely · §3 the rule for a
multimapper · §4 the record · §5 pass one · §6 the second pass · §7 what the EM sees · §8 reproducibility,
cost, the cache · §9 gates · §10 stages and the A/B · §11 the decisions that are the owner's · §12 what to
expect.

## 1. The problem and its size

A fragment whose read name carries `NH > 1`, or secondary alignments, is a multimapper. Pass one never
deposits it into the calibration accumulator (`bam_scanner.cpp`: `deposit_to_accumulator` is reached only
under `is_unique_mapper`, lines 1812 and 1902; a multimapper's hits go to the EM's buffer alone, line
1909). So every region whose fragments all multimap reads as EMPTY to calibration, however dense its gDNA,
and everything downstream reads the emptiness as fact: the composition solve, the landscape (one more zero
anchor), the capture reader (an empty object beside a dense one is depleted), the locus prior (a gDNA count
of 0 and, under the reader, a yield of 0, which the EM reads as "cannot emit"). The paralog gate shows all of
it at once: two identical 500 bp exons at gDNA abundance 100, calibration count 0 / 0 on 155 bp of support
beside intergenic gDNA at 0.23–0.31 fragments/bp; the EM calls 154 / 0 (shipped) and 123 / 79 (with the
reader) against a truth of 75 / 71, totals 154 and 202 against 146.

The populations, from the production scans (`prototypes/2026-10-10_honest_reader/prod_*_main.log`) and the
footprint (`mm_footprint_*.json`; a "region with multimappers" is one holding at least one multimapper
alignment start):

| library | fragments | multimapper molecules | multimapper alignments | mean alignments per molecule | regions holding multimappers and no unique fragment (exon regions) | regions where multimappers outnumber unique fragments, holding this share of all multimapper alignments | gap-held fragments today, drained in |
|---|---|---|---|---|---|---|---|
| LBX0190 (plasma) | 155,352 | 11,720 | 104,140 | 8.9 | 5,488 (762) | 6,504 — 95 % | 3,216 in 0.9 s |
| MO_3021 (plasma) | 875,670 | 70,856 | 583,407 | 8.2 | 23,443 (7,018) | 35,781 — 94 % | 8,760 in 1.2 s |
| VCaP mix (deep) | 18,568,456 | 136,188 | 544,065 | 4.0 | 10,099 (5,305) | 13,792 — 72 % | 355,180 in 13.1 s |

Multimappers are 7.5 % of a plasma library's molecules and they land, with their eight placements each, in
the few thousand regions where calibration counts almost nothing. No benchmark panel saw this: every panel
is an oracle BAM (`NH:i:1` on every record), and the only substrates with real multimappers are the aligned
scenarios (`tests/scenarios_aligned`, minimap2) and the real libraries.

## 2. The machinery as it stands, precisely

**The scanner's per-fragment loop** (`bam_scanner.cpp` 1552–1950). One read name → `records` → `nh` (the
NH tag, 1 if absent) → `is_multimap = nh > 1 || has_secondary` → `group_records_by_hit` (by the HI tag when
present, sorted by HI; else by pairing secondaries, `pair_multimapper_reads` 850–995: transcript-set
intersection first, same-reference closest distance next, then every remaining R1 × R2 combination, then
singletons) → `all_hits`, `num_hits = max(nh, all_hits.size())`, `is_unique_mapper = num_hits == 1`. Then
per hit: `build_fragment` → `AssembledFragment` (blocks with `ref_id`, introns, `nm`) → `_resolve_core` →
`RawResolveResult cr` (`t_inds`, `splice_type`, `align_strand`, `sj_strand`, `chimera_type`, and the gap
hypotheses: `gap_intron_offsets / gap_introns / gap_sj_strand / gap_supporting`). A hit with empty `t_inds`
is intergenic: appended to the EM buffer and deposited only when `is_unique_mapper`. A chimeric hit is
appended to the EM buffer as chimeric and never deposited. A resolved non-chimeric hit: strand-model
observations and the deposit only when `is_unique_mapper`; appended to the EM buffer always; `n_buffered_mm`
counts a multimapper's appended hits.

**The deposit adapter** (`deposit_to_accumulator`, 1565–1691) turns one hit into an `OfferedFragment`
(`accumulator.h` 237): the extent on ONE reference (`[start, end)`, leftmost block start to rightmost block
end, the mate gap included; a hit with blocks on two references is not offered), the observed CIGAR-N
introns on that reference (de-duplicated, cut under every hypothesis), `align_strand` and the OBSERVED
`sj_strand` from `cr`, and the gap hypotheses re-presented as `GapHypothesis` spans (implied introns, the
implied strand, the supporting transcripts). It calls `ws.acc_set->at(ref_id).deposit(offered, scratch)`.

**The accumulator's deposit** (`calibration/accumulator.cpp` 468–720; the executable specification
`tests/native/_accumulator_reference.py`, `Accumulator.deposit` 712–830): the strand check first (NONE or
AMBIGUOUS → `dropped_strand_undefined`, and it wins over the deferral because that population is not
recoverable by the second pass), the clip to the reference (`L` is the clipped length), then ARBITRATION:
every hypothesis's `L` via `hypothesis_length` (the observed introns unioned with the implied ones,
normalised, clipped), the one filter `L ≤ max_length` — unless it would empty the set, in which case all
stand — then: more than one survivor → the WHOLE offered fragment (every hypothesis, filtered ones included,
"what was OFFERED") is appended to this reference's `DeferredFragments` bank with the clipped extent, the
umbrella census is recorded, `deferred_undetermined_gap` increments, outcome `kDeferred`; exactly one →
the write-out: the sj the path uses, `region_start_count` / `region_end_count` at the path's first and last
covered base, the crossings per contiguous segment (two channels on the spliced bank, four on the
unspliced: `unspliced_count`, `unspliced_inv_length_sum += 1/(L−1)`), the strict spans, the conserved mass
per slice over boundaries and sj together, the sj `count` and `inv_length_sum`, `contained_count` when the
whole path lies in one region, the length pool, `deposited_lengths`, `deposited`.

**The bank** (`DeferredFragments`, `accumulator.h` 292; `scan_payload.DeferredFragments` 290): one record
per held fragment — `ref, start, end, align_strand, sj_strand` (`_DEFERRED_RECORD_FIELDS`), the observed
introns as a CSR, the hypotheses as a CSR (`hypothesis_sj_strand`, implied introns as a CSR, supporting
transcripts as a CSR). Every array `int64`. Per worker, per reference; merged across workers by CSR
concatenation (`merge_from`); exported in ONE canonical order — sorted on the record's own content, the
specification's key `(ref, start, end, align_strand, sj_strand, observed_introns, hypotheses)` with
Python's prefix rule (`canonicalise`, 117–175), per reference, and the references concatenated in reference
order, which is already canonical since every record carries its reference (`bam_scanner.cpp` 2104–2112).
At the door the payload refuses a bank whose size is not `qc.deferred_undetermined_gap`, a record with
fewer than two hypotheses, and a gap census whose three `deferred_*` do not sum to that counter
(`scan_payload.py` 431–445, 757–775).

**The second pass** (`pipeline._drain_side_buffer` 349–440, `second_pass.py`). Between the scan and
calibration, from pass one's state alone. `score_held_fragments`: per record, per hypothesis, `L` from
`Accumulator.length_under` (the same C++ that will compute it at drain time); `rho`: a spliced path takes
the bottleneck (minimum) of the `sj_inv_length_sum` of the annotated sj it uses (0 if unannotated), the
genomic path the bottleneck of `boundary_unspliced_inv_length_sum` over the DISTINGUISHING boundaries — the
boundaries at region bounds inside the contested introns' ranges, endpoints included (`_distinguishing_boundaries`)
— scored over exactly what the competing paths jump so that ∅ is not penalised for touching more objects
than a path that jumps them; `f(L)`: `rna_pmf[L]` for a spliced path, `global_pmf[L]` (the unconditional
anchor) for the genomic one; the strand term only when no motif was observed (`strand_terms`: `p` / `1−p`
for a spliced path against the implied strand, `½` for the genomic path, because DNA is double-stranded).
`combine_factors`: the product normalised within the record, applied strongest-evidence-first with an
all-zero factor dropped as uninformative and a partial zero kept decisive; `undecided` when the lead is
tied. `choose_hypotheses`: ONE multinomial draw per record, vectorised over the canonical order from one
seeded stream (`second_pass_seed`), never keyed on content. `drain`: each record re-enters
`Accumulator.deposit` on ITS reference with the chosen hypothesis alone, so arbitration is degenerate and
the ordinary rules decide; the delta is added to pass one's arrays; the bank is emptied, its counter and the
census's `deferred_*` zeroed; `DrainQC(offered, deposited, dropped_*, chose_genomic, chose_spliced,
census_before)` with `deposited + dropped == offered` refused otherwise. `lift_choices` replays the whole
library's choices inside origin partitions for the oracle, keyed on `_DEFERRED_RECORD_FIELDS`.

**What the EM sees** (`scoring.cpp` 205–230, 577): every non-chimeric hit of a multimapper as a candidate,
with the gDNA candidate per hit by the one rule (`gdna_can_explain`) and pool-separated pruning
(`flush_mm_group`); the locus priors from the calibration result (`priors.assemble_priors`).

## 3. The rule for a multimapper, derived

A unique mapper is a fragment with one placement and one or more gap hypotheses at it. A multimapper is a
fragment with P ≥ 2 placements, each with its own gap hypotheses. The candidate set is the union

    C = { (p, h) : p a placement of the fragment, h ∈ H_p },

and the molecule took exactly one element of it. Pass one holds the fragment WHOLE with all of C (as the
gap bank holds every hypothesis, filtered ones included); the second pass scores every (p, h) with the
existing three factors and ONE draw over C; the drain deposits the chosen (p*, h*) at placement p*'s
reference with hypothesis h* alone. The unique mapper is the P = 1 case of the same rule, and its behaviour
is unchanged bit for bit (§5, §6).

**The score.** `score(p, h) = A_p · g_p(h) · f(L_{p,h}) · s(p, h)`, normalised over C:

- `A_p` — the owner's term, the abundance of uniquely aligned fragments at the placement: the pass-one
  unique traffic at the objects the placement's own path deposits on, read through the deposit rule's
  own geometry. Take the placement's path with NO implied intron (its observed introns only — the
  hypothesis ∅ at p). If that path is CONTAINED in one region `r`, `A_p = contained_count[r] /
  contained_eff_length(ℓ_r, global_pmf)` (`effective_length.contained_eff_length`: fragments per start
  opportunity, the same frame calibration's density uses; both strand columns summed). If it crosses
  boundaries and uses no annotated sj, `A_p` is the bottleneck of `boundary_unspliced_inv_length_sum` over
  the boundaries it crosses. If it uses annotated sj (an observed spliced placement), the bottleneck of
  `sj_inv_length_sum` over them — the spliced bank carries counts but no `inv_length_sum`, and this is
  exactly today's rule for a spliced hypothesis. The bottleneck, not a mean, for the reason the existing
  rule gives: a molecule that took this placement was present at every object on it.
- `g_p(h)` — today's within-placement gap score, unchanged: `rho(h)` as §2 defines it (sj densities for a
  spliced path, the contested boundaries for ∅), normalised over `H_p`; 1 when `H_p` is the unspliced path
  alone. For a P = 1 record `A_p` is a common factor and is NOT applied, so the record's scores are today's
  bit for bit — a design constraint, since the gap-held population is gated by
  `tests/test_second_pass_scoring.py` and must not move under a multimapper change (one mechanism at a time).
- `f(L_{p,h})` and `s(p, h)` — per placement and hypothesis, exactly as today: `L` from `length_under` on
  that placement's accumulator; the pmf by whether `h` splices; the strand term from that placement's
  `align_strand` and `h`'s implied strand when no motif was observed there.

The product goes through `combine_factors` unchanged, so its two rules carry over and are the whole of the
degenerate-case policy — no new constant, floor or tie-break:

- every placement with `A_p = 0` (the paralog pair: neither identical exon holds a unique fragment) →
  the factor is flat zero among the survivors, uninformative, DROPPED; the draw falls to `f` and `s`, which
  tie across identical placements → 1/P each, and the record counts in `n_undecided`;
- some placements with `A_p = 0` beside others with `A_p > 0` → the zero is decisive and stays hard: the
  placement with no unique evidence is never chosen. This is the owner's rule read literally, and it has a
  named consequence (§11, decision 1).

**What is derived, what is chosen.** The union C, the one draw, the drain through `deposit`, `f` and `s`
per placement, and `combine_factors`' zero rules are the existing design applied to the enlarged candidate
set. The one new quantity is `A_p`, and it is read off the deposit rule's own objects and densities with no
constant; the only choice in it is the bottleneck over a path's objects, the choice the existing rule
already made and argued.

## 4. The record

One bank for both populations: a record is a fragment, a fragment holds PLACEMENTS, a placement holds
hypotheses, a hypothesis holds introns. The gap-held unique mapper is a record with one placement, and
every existing field keeps its name and meaning at the placement level:

```
placement_offsets        int64[n_fragments + 1]           NEW: fragments → placements
ref, start, end          int64[n_placements]              per placement (today: per fragment)
align_strand, sj_strand  int64[n_placements]              per placement
observed_intron_offsets  int64[n_placements + 1]          per placement; observed_introns int64[2·n_observed]
hypothesis_offsets       int64[n_placements + 1]          placements → hypotheses (today: fragments → hypotheses)
hypothesis_sj_strand, hypothesis_intron_offsets, hypothesis_introns, hypothesis_t_offsets, hypothesis_t   unchanged
```

`n_fragments = len(placement_offsets) − 1`, `n_placements = len(ref)`. A fragment's hypotheses are one
contiguous run of the flat hypothesis axis — `[hypothesis_offsets[placement_offsets[i]],
hypothesis_offsets[placement_offsets[i+1]])` — so `choices[i]` stays a LOCAL index into the record's run
and identifies (p, h) uniquely; `choose_hypotheses` is unchanged except that its runs are records' rather
than placements' (its cumulative-sum arithmetic is over the record's run). Within a record the placements
are stored in THEIR canonical order (each placement's own key, below), so a record's content does not
depend on the aligner's hit order or the HI numbering.

**The canonical order** extends the specification's key by one level: a record compares as the sequence of
its placements, each placement as today's record key `(ref, start, end, align_strand, sj_strand,
observed_introns, hypotheses)`, with the prefix rule at both levels. Two records that tie are identical
records, as before. A P = 1 record's key is today's key, so the order among gap-held records is today's.
Because a multi-placement record spans references, the bank is no longer "per reference, concatenated in
reference order": the export merges every per-reference gap bank and the set-level multimapper bank (§5)
and canonicalises the union with ONE sort (`canonicalise` is idempotent and is already called at export).
`_DEFERRED_RECORD_FIELDS`, which `lift_choices` keys records on, becomes the per-placement tuple sequence
(the key the sort uses), defined once in `scan_payload` beside the arrays.

**The door** (`DeferredFragments.from_dict`): the nested CSRs re-derived as today; a record must hold at
least two (p, h) in total (a P = 1 record at least two hypotheses, as today); the bank's size must equal
`qc.deferred_undetermined_gap + qc.deferred_multiple_placements` (§5); the gap census's three `deferred_*`
must still sum to `deferred_undetermined_gap` alone.

## 5. Pass one

**The accumulator** (`accumulator.cpp`, the specification alongside). `Accumulator::deposit` is factored
into the three steps it already performs — the strand check and clip, `arbitrate(offered, start, end,
scratch) → survivors` (the lengths, the `L ≤ max_length` filter, the all-stand rule), and
`deposit_survivor(...)` (the write-out) — with `deposit` itself unchanged in behaviour and the existing
parity battery pinning it. `AccumulatorSet` gains one entry point and one bank:

```
DepositOutcome AccumulatorSet::offer(const Placement* placements, std::size_t n, OfferScratch&)
    // Placement = { ref_id, OfferedFragment }  — every non-chimeric hit of ONE fragment
```

whose rule is: per placement, the strand check (an undefined strand excludes THAT placement) and the clip
to its reference (an empty clip excludes it); if no placement remains → `dropped_strand_undefined` if every
exclusion was the strand, else `dropped_empty` — ONE count for the fragment; identical placements
(the same ref, clipped extent, strands, observed introns and hypotheses, which the pairing combinatorics can
produce) collapse to one, as two transcripts implying one path are one hypothesis; `arbitrate` on each
remaining placement; the survivor set is the UNION over placements of `(p, h)` with `L ≤ max_length`, and
if that union is empty every (p, h) stands (the ordinary `kTooLong` rejection counts the chosen one at the
drain, as today); then:

- exactly one placement remains → `at(ref).deposit(offered)` — THE existing path, bit for bit: a unique
  mapper, or a multimapper every other hit of which was chimeric or excluded, deposits or is gap-deferred
  exactly as today and counts under today's counters (`deferred_undetermined_gap`, the gap census);
- several placements and exactly one survivor → `at(ref_{p*}).deposit` with that single hypothesis (the
  fragment is a multimapper by NH, but only one placement and path is a molecule this library contains);
- several placements and two or more survivors → the whole fragment (every placement, every hypothesis)
  is appended to the set-level bank, `++deferred_multiple_placements`, outcome `kDeferred`. The gap census
  is NOT recorded for these records: its classes answer "how was the gap resolved" at one placement, and a
  multi-placement record has no single answer; the record carries its placements, so a census is derivable
  from the bank and needs no second ledger.

The new counter `deferred_multiple_placements` joins `DepositCounters` / `ScanQC` and the identity becomes
`deposited + deferred_undetermined_gap + deferred_multiple_placements + dropped_* == offered`, where
`offered` counts FRAGMENTS (one `offer` per read name), not hits. The set-level bank merges across workers
exactly as the per-reference banks do (`AccumulatorSet::merge_from`) and is sorted into the one bank at
export (§4).

**The scanner** (`bam_scanner.cpp`, the hit loop). The adapter `deposit_to_accumulator(frag, cr)` becomes
`collect_placement(frag, cr)`: it builds the `OfferedFragment` exactly as today but appends it, with its
backing storage (observed introns, hypothesis spans), to a per-worker `placements` arena cleared per
fragment, instead of depositing; the intergenic-hit branch and the resolved non-chimeric branch call it for
EVERY fragment, unique or not (the `is_unique_mapper` condition around the deposit goes); chimeric hits are
not collected (a fragment with blocks on two references deposits nothing, today's rule). After the hit
loop, one call: `ws.acc_set->offer(placements)`. For a unique mapper this is one placement and the set
delegates to `deposit`, so pass one's tally on `NH = 1` data is unchanged to the byte (gate §9.6). The EM
buffer appends are untouched: a multimapper's intergenic hits are still not appended (the EM never sees
them today and this design does not change what the EM sees, §7), but they ARE placements for
calibration — gDNA sits anywhere, and an intergenic placement is where a repeat's gDNA often is.

**The splice census, strand-model observations and `n_intergenic`** stay keyed on `is_unique_mapper` /
`any_hit_resolved` as today; this design moves only the deposit.

## 6. The second pass

`score_held_fragments` iterates records; within a record it iterates placements and, within each, the
hypotheses, with ONE accumulator per reference as today (`_Accumulators`). Per placement it computes `A_p`
(§3) once — `region_of_pos` for the ∅-path's first and last covered base, the crossed-boundary range by the
same `searchsorted` the distinguishing-boundary rule uses, `sj_edge_ids` for observed introns — and then
today's per-hypothesis loop verbatim: `L` via `length_under` on that placement's accumulator, `rho(h)`,
`f(L)`, `s`. The record's scores are `combine_factors(A · rho, f, s)` over the record's run, with `A` the
placement's `A_p` broadcast over its hypotheses and omitted (≡ 1) when the record has one placement. The
`HypothesisTerms` gain `placement_abundance` (per hypothesis, its placement's `A_p`) so a regression can be
attributed to the new term alone, as the three existing terms are kept apart for that reason.

`choose_hypotheses`: unchanged in form; its runs are records (§4).

`drain`: for record `i` the local choice maps to (p*, h*); the replay is `accumulators[ref_{p*}].deposit(
start_{p*}, end_{p*}, observed_introns_{p*}, align_strand_{p*}, sj_strand_{p*}, hypotheses=(h*,))` —
the same call as today on a placement's own reference. `chose_genomic` / `chose_spliced` count by `h*`.
`DrainQC` gains `offered_multimapper` (pass one's `deferred_multiple_placements`), so `offered ==
deferred_undetermined_gap + deferred_multiple_placements` is checkable, and keeps `deposited + dropped ==
offered`. After the drain both counters are 0 and the bank is empty.

`lift_choices`: the key is the record's placement-tuple sequence (§4); on oracle substrates (no
multimappers) nothing changes.

## 7. What the EM sees

Nothing new. The EM buffer still receives every non-chimeric hit of a multimapper with `num_hits`, and the
EM allocates it among its candidates by abundance and the locus priors. What changes is the locus prior:
calibration now counts a multimapper's fragment at the placement the second pass drew, so a repeat exon's
locus has a gDNA count, a density and a capture weight where it had zero. The second pass's placement and
the EM's allocation of the same fragment are not reconciled: they answer different questions (the field
calibration reads, against the final count) and run at different stages, exactly as the gap-held
population's drawn path and the EM's allocation are not reconciled today.

## 8. Reproducibility, cost, the cache

- **Order.** The bank's content-based canonical order covers both record kinds; one draw stream in that
  order; the drain per record. `test_the_DRAINED_payload_is_byte_identical_at_1_2_4_8_WORKERS` extends to
  a fixture with multimappers (gate §9.6).
- **Cost.** The second pass is a Python loop over records at ~37 µs per gap-held record (VCaP: 355,180 in
  13.1 s). A multimapper record carries ~8 placements on a plasma library and ~4 on VCaP; at ~40 µs per
  placement: LBX0190 11,720 × 8 → ~4 s, MO_3021 70,856 × 8 → ~23 s, VCaP 136,188 × 4 → ~22 s. Acceptable
  for the landing; a native scorer is a later speed item and is not part of this fix. Memory: ~300 bytes per
  record → ~40 MB on VCaP.
- **The cache.** The payload's schema digest is part of the scan-cache key (`scan_cache.py`), so every
  cache written before this change is refused at the door by construction; `build_scan_cache.py --force`
  stays the documented step, since the key hashes the deposit rule's field list and not `resolve.cpp`.
- **Hits absent from the BAM.** `NH` may exceed the alignments present (`--outSAMmultNmax`); the draw is
  over the present placements, the only possible rule; `num_hits = max(nh, hits)` is untouched.

## 9. Gates, each written and failing first; each deliberate defect must fire one

| gate | what it holds |
|---|---|
| `spec-offer` | the specification's `offer` on brute-force cases: P = 1 is byte-identical to `deposit` (the whole existing battery replayed through `offer`); identical placements collapse to one; a strand-undefined placement excludes itself and only itself, every placement undefined → one `dropped_strand_undefined`; every (p, h) over the limit → all stand, deferred; one survivor across several placements → deposited at its placement; two or more → deferred whole, every placement and hypothesis retained |
| `spec-canonical` | the flattened bank round-trips, is independent of deposit order and of the hits' order within a fragment, and a P = 1 record's key is today's |
| `native-parity` | `AccumulatorSet::offer` against the specification case by case, and the ten-thousand-random-fragments battery extended with random multimappers (P ∈ 1..4, mixed references, mixed strands) — arrays, counters and the bank byte-identical |
| `score-placements` | on a two-contig hand-built fixture with known traffic: unique traffic 3 : 1 between two placements scores 3 : 1 and the draw follows across seeds; both zero → ½ each and `n_undecided` counts it; one zero beside one positive → the positive always; a P = 1 record's scores are bit-identical to today's (every existing scorer test unchanged); the processed-pseudogene geometry (a spliced placement whose genomic path is long beside a contiguous placement) is decided by `f` as the length law says |
| `drain-identity` | draining (p*, h*) equals offering that placement with that hypothesis alone, bit-identically, on each reference; the bank empties; `offered_multimapper` reported; `deposited + dropped == offered` |
| `workers-identity` | the drained payload byte-identical at 1, 2, 4 and 8 workers with multimappers in the fixture; and on an `NH = 1` library the WHOLE change is a numeric no-op — `rename_identity.py --check` against a frozen reference of the ladder's scan caches and a frozen VCaP payload restricted to its unique fragments |
| `the paralog gate` | `TestParalogMultimapping.test_gdna_sweep[gdna_100]` turns green: the collapse branch is deleted as the test itself instructs and the even-split assertion covers it; the total approaches 146 |
| `aligned-scenarios` | every class in `tests/scenarios_aligned/test_multimap_counting.py` scored against the simulator's truth before and after (`run_benchmark`), per scenario, the table in the landing note; none worse beyond its own noise, the paralog and pseudogene classes better |
| `real-libraries` | LBX0190, MO_3021, VCaP through the production path, pinned, before and after: the footprint re-measured on the DRAINED payload (regions with multimappers and no unique fragment must now carry deposits), VCaP's gDNA fraction against 0.2518, the capture weights at those regions, the second pass's time and the share of records in each `A_p` class (all zero / mixed / all positive) |

## 10. Stages and the A/B

1. **The specification first** (`tests/native/_accumulator_reference.py`): `DeferredFragment` with
   placements, `Partition`-level `offer`, the canonical key, `deferred_arrays` with `placement_offsets`,
   `drain` through the placement; gates `spec-offer`, `spec-canonical`, `drain-identity` written against the
   current spec and failing (it has no `offer`). The derivation is the spec; nothing native until it holds.
2. **Native, in a worktree** (`cpp_prototype_isolation` recipe): `Accumulator::arbitrate` factored out,
   `AccumulatorSet::offer` and its bank, the export's union sort, the scanner's `collect_placement` and the
   single `offer` per fragment, the counter; gates `native-parity`, `workers-identity`.
3. **The second pass** (Python): the record level, `A_p`, the terms, the drain's placement, the QC, the
   lift's key; gate `score-placements`.
4. **The A/B** against main on the same conditions: the aligned scenarios against truth; the three real
   libraries through the production path; the `NH = 1` no-op proof. One mechanism: no other change rides
   along, the reader's rules untouched.
5. **Landing**: the suite re-derived, `build_scan_cache.py --force`, `preflight.py --full`; docs by the move
   rule — `DESIGN.md` (the accumulator's multimapper rule beside the deposit rule), `EQUATIONS.md` §1 (the
   offer and `A_p`), `ISSUES.md` (the entry closed with the paralog number and the real-library footprint
   after), `TESTING.md`, CLAUDE.md's baseline and its paralog-gate sentence. The commit is the owner's.

## 11. The decisions that are the owner's

1. **A placement with no unique evidence beside one with some.** The literal rule — proportional to the
   abundance of uniquely aligned fragments — makes the zero decisive: a silent processed pseudogene with
   no unique fragment never receives a multimapper while its expressed parent does, so the pseudogene's
   gDNA is attributed to the parent and the pseudogene region stays empty in calibration. (A perfect
   paralog pair, both without unique fragments, splits evenly; the paralog gate is that case.) Recommended:
   build the literal rule, and let the `real-libraries` gate report the share of multimapper records in the
   mixed class before deciding whether a gDNA-uniform component is wanted — a new rule, A/B'd on its own.
2. **A multimapper's intergenic placements** are placements for calibration (gDNA sits anywhere) although
   the EM never sees them. Recommended: yes, as §5 states.
3. **One bank or two.** One bank with a placement level (§4), so there is one order, one draw stream, one
   drain and one door; the alternative (a separate multimapper bank beside the gap bank) duplicates all
   four. Recommended: one.
4. **The gap census for multi-placement records**: not recorded (§5). Recommended as stated.

## 12. What to expect

On the paralog gate: each exon's 36 true gDNA fragments reach calibration at about half each of the 72
shared molecules, the density reads the intergenic level, the reader reads weights near the intergenic
regions', the locus priors hold the gDNA, and the EM's total falls from 202 toward 146 with the split within
the even-split assertion. On the real libraries: the regions that today hold multimappers and no unique
fragment (5,488 / 23,443 / 10,099) receive deposits in proportion to the draws; the plasma libraries' gDNA
fractions, which rose 3–15 % under the reader with no truth to judge them, move again and VCaP's is read
against 0.2518. On every oracle panel and every `NH = 1` fragment: nothing moves, by construction, and that
is the no-op gate.

## 13. Built (2026-10-11), as designed, with one finding

The owner's decisions (2026-10-10): the hard zero as the rule for now; intergenic alignments are
placements and feed the intergenic background; one bank with a placement level ("two banks" meant a
separate multimapper bank beside the gap bank, which would have duplicated the order, the draw stream,
the drain and the door); no gap census for multi-placement records.

**What was built, in the order of §10.** (1) The specification: `Placement`, `DeferredFragment.placements`,
`Accumulator.offer` with `_admit` / `_arbitrate` / `_deposit_survivor` factored out of `deposit` (unchanged
bit for bit), `_locate`, the three-level flattening with `placement_offsets`, the record-level canonical key,
the placement-aware drain with `offered_multimapper`; 21 gates in `tests/native/test_multimapper_offer.py`.
(2) Native: `DeferredFragments` gains `placement_offsets`, `append_placement` / `close_record`,
`copy_placement`, `offered_of`, `compare_placements`, `clear` and a record-level `canonicalise`;
`Accumulator::admit` / `arbitrate` / `deposit_survivor` / `restore_gap_census`; `AccumulatorSet::offer` with
the set-level bank and counters, merged across workers; the scanner collects every non-chimeric hit as a
placement into per-worker arenas (the resolver's result dies with the hit, so implied introns and
supporting transcripts are copied) and offers the fragment once after its last hit; the export merges every
reference's gap bank with the set's bank and sorts the union; `rigel.native.AccumulatorSet` binds the set
for the parity gate. `tests/native/test_multimapper_native_parity.py`: the one-placement battery, the
multi-placement cases, the interleaved union order and 2,000 random fragments (one to three placements,
two references, random strands, observed introns and hypothesis sets) — byte-identical throughout.
(3) Python: the payload's placement level and door (`deferred_multiple_placements`, the per-record
"at least two pairs" check, the multi-record count against the counter), `second_pass.placement_abundance`
exactly as §3 defines `A_p`, `combine_factors(..., placement_abundance)` as its own factor between the
strand term and `rho`, the record-level draw, the placement-aware drain, the lift's per-record key,
`DrainQC.offered_multimapper`, the pipeline's log line; `tests/test_second_pass_placements.py` (3:1 unique
contained depth scores exactly 0.75 / 0.25 and the draw follows within four binomial standard deviations; no
unique evidence → ½ / ½ and undecided; one side only → a hard zero and nothing reaches the silent copy; the
drain deposits exactly where it drew); the two-contig drain fixture carries four `NH = 2` fragments, so the
worker-identity gate covers multi-placement records.

**The suite under the build: 0 failed / 2,835 passed** (2,805 + 30 new gates); `ruff` clean; the docs gate
passes. Docs by the move rule: `DESIGN.md` §4.2, `EQUATIONS.md` §1.8 and §10, `TESTING.md`.

**The paralog receipt** (`paralog_dissect_site3.log`, gDNA 100): calibration now counts the two identical
exons at 0.24 and 0.28 fragments/bp against 0.24 beside them (from 0 / 0), the reader's weights there are
1.0 and 1.2 (from 1/2000 of the intergenic regions'), and the EM's total is 160 against a truth of 146 (154
while the priors were zero, 202 with the reader reading the zero as depletion). Every aligned multimapper
scenario passes: 49 / 49.

**The one finding, flagged to the owner.** The EM's split of two SEQUENCE-IDENTICAL templates is its own
degeneracy — bimodal, even or a vertex — and which one was decided by the accidental symmetry of two EMPTY
priors. With both loci now holding the drawn multimappers, their priors differ by the draw's noise and the
degenerate EM leaves the knife-edge for a vertex at every gDNA level and strand level (0 / 387 at gDNA 20,
0 / 551 spliced, 0 / 367 at ss 0.65 and 0.9), where the suite had asserted an even split within 15–30 %.
The totals are right in every case. The assertions were rewritten to what the EM can be held to on
identical templates — the TOTAL within a band, and the SHAPE: even within the tolerance or a clean vertex,
never between (`assert_identical_paralogs`) — with the old `gdna == 100` collapse branch generalised rather
than deleted. This is the judgment call in the landing: the even split was never a stable answer, and
`TestDistinguishableParalogs` remains the arm where a split is identifiable.

**The NH = 1 no-op, proven.** `rename_identity.py --freeze` on main (6d08f790) over the ladder's
`gdna_g00_ss_0.99_nrna_mid_capture_off` (an oracle BAM, `NH:i:1` throughout), `--check` under the build:
BIT-IDENTICAL — the quant digest (8,750 rows) and every payload and calibration array unchanged
(`identity_main.json`, `identity_check.log`).

**The real libraries, pinned production path, main (6d08f790) against the build** (`prod_*_main2.json`,
`prod_*_site3.json`):

| library | wall s | gDNA fraction | mRNA | nRNA | gDNA (EM) | second pass |
|---|---|---|---|---|---|---|
| LBX0190 | 22.1 → 26.4 | 0.0976 → 0.0997 | 123,407 → 123,476 | 8,958 → 8,580 | 11,953 → 12,262 | 3,216 held in 0.9 s → 13,703 held (10,371 on their placements) in 3.8 s |
| MO_3021 | 70.1 → 83.5 | 0.1631 → 0.1653 | 624,525 → 624,445 | 74,246 → 72,518 | 83,128 → 84,936 | 8,760 in 1.1 s → 77,452 (68,182 on their placements) in 15.2 s |
| VCaP mix (truth 0.2518) | 531.5 → 468.3 (both under a running suite; not a timing) | 0.2401 → 0.2528 | 13,334,271 → 13,323,232 | 746,774 → 523,324 | 4,365,879 → 4,600,368 | 355,180 in 13.2 s → 481,880 (126,154 on their placements) in 27.1 s |

Of LBX0190's 11,720 multimapper molecules 10,371 are held on their placements (the rest resolve to one
admissible placement, or are excluded whole); of MO_3021's 70,856, 68,182; of VCaP's 136,188, 126,154. The
second pass's cost is the projected one (+3, +14, +14 s at ~0.1–0.4 ms per held multimapper). **VCaP's gDNA
fraction moves from 0.2401 to 0.2528 against its truth of 0.2518**: the 234,000 fragments the EM gains as gDNA
are the 223,000 it loses as nascent RNA — the multimapped gDNA at repeat loci that an empty prior had let it
call RNA. The plasma libraries' gDNA fractions rise 2 % relative, with no truth to judge them.

The ISSUES entry is CLOSED with these numbers; CLAUDE.md's baseline is re-derived in the worktree (0 failed /
2,835 passed). The commit and the landing are the owner's.
