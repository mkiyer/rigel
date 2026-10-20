# One length rule: design (roadmap item 11, rule B)

2026-10-19 · code base `1fc2f609` · for the owner's review; rule B chosen by the owner the same day. Every claim
about today's code cites a line re-read for this document. Approved and landed 2026-10-20; the measurement is the
roadmap's item 11 (`PAIRED_END_ROADMAP.md`).

## 1. The change in brief

The EM scores each candidate transcript by the fragment length the molecule would have on it. Today that length
comes from a rule of its own, which cuts every intron of the candidate between the fragment's two ends, even one the
reads sequenced through. Calibration measures the same molecule differently.

After this change, a candidate's length is **the length of the path it implies, measured as calibration measures a
deposit**: the fragment's span less its observed junctions and less the candidate's own introns that lie in an
unsequenced gap, within K of the gap's ends. A candidate never cuts an intron the reads cover.

So one molecule has one length on a candidate, whatever its read length or mate overlap, and the EM and calibration
share one definition. The single-block branch goes: a fragment with no gap is its span on every candidate, which is
what that branch already gave.

## 2. Today

**The EM's length** (`compute_frag_lengths`, `resolve_context.h:1148-1177`), computed for every non-chimeric
fragment (`:1509-1511`):

- **one block** (the mates overlap or abut): the block's length, for every candidate (`:1159-1163`);
- **two or more blocks:** `tx_frag_length` (`:1111-1141`), the span less every intron of the candidate lying strictly
  inside it. It reads only the two ends: "a read-through of a whole intron is not modelled" (its docstring).

**Calibration's length** (`hypothesis_length`, `accumulator.cpp:485-505`): the bases of the span outside the observed
junctions and the hypothesis's introns, by one segment arithmetic (`normalise_introns`, `build_segments`,
`:432-481`). Its comment: "ONE definition of L ... Deriving it any other way invites two formulas for one quantity".

**Calibration's hypotheses** come from the resolver's gap search (`enumerate_gap_hypotheses`, `resolve_context.h:946-1032`):

- the gaps are the holes between consecutive blocks that are not an observed junction;
- a candidate an intron of which covers an aligned block by more than K implies nothing
  (`transcript_contradicts_blocks`, `:702-724`);
- any other candidate implies its introns that overlap a gap and lie within K of the gap's ends
  (`collect_transcript_introns_in_gap`, `:1039-1065`).

A fragment whose hypotheses disagree on length is held for the second pass (DESIGN §4).

**The validator** applies the same K rule to its witness path (`validate_pair`).

**The EM reads the length twice:** the fragment-length likelihood (`scoring.cpp:547, 818`), and the coverage weight,
which places the fragment on the transcript at `[start, start + length)` (`scoring.cpp:243-262`).

## 3. The rule

For a non-chimeric fragment with aligned blocks on one reference, observed junctions O, span [s, e) (leftmost block
start to rightmost block end, the extent the deposit uses), and a candidate t:

- **G_t**: t's introns that overlap an unsequenced gap and lie within K of its ends. This is
  `collect_transcript_introns_in_gap` over the fragment's gaps, as the gap search computes it.
- **L_t**: the bases of [s, e) outside O ∪ G_t, by the deposit's segment arithmetic.

**On a candidate the reads support,** L_t is exactly the length of the hypothesis calibration deposits for t. On those
candidates today's rule already agrees, except for an intron of 1 to K bases inside a block, which it cuts and this
rule keeps (§5).

**On a candidate the reads contradict** (an intron covered by more than K aligned bases), the covered intron is
retained: the reads show its bases. Its gap introns are still spliced, as t splices them, so L_t is t's molecule with
the covered intron retained. Calibration assigns such a candidate no hypothesis at all. The EM still scores it, with an
overhang penalty for the covered bases; this gives it a length a molecule could have.

**Every other case** follows from the rule:

- a fragment with no gap (one block, or blocks joined only by observed junctions) has L_t = the span less O on every
  candidate;
- the nascent span and single-exon transcripts have no introns, so L = the span less O;
- a chimeric fragment still has no length (unchanged).

The item's example (T2's 60-bp intron at 1200–1260; a 170-bp molecule from 1150 to 1320):

| read as | T2 today | T2 under the rule |
|---|---:|---:|
| 2×100, mates overlapping, the intron inside the block | 170 | 170 |
| 2×75, each mate 25 and 15 bases into the intron | 110 | **170** |
| mates stopping at the exon edges, the whole intron in the gap | 110 | 110 (and 170 on the nascent span, held) |

## 4. One definition, in code

- `normalise_introns` and `build_segments` move from `accumulator.cpp` into a small native header that the accumulator
  and the resolver both include. L then has one formula in C++, as it has one in the specification
  (`tests/native/_accumulator_reference.py`'s `_hypothesis_length`).
- The resolver computes each candidate's G_t in the gap search's candidate loop, which already collects them for the
  candidates the reads support. A contradicted candidate's are collected for its length only, and join no
  hypothesis. Each candidate's length is written from the same G_t its hypothesis uses, so the two cannot drift.
- `compute_frag_lengths`, its single-block branch and `tx_frag_length` are deleted. The coverage weight's docstring,
  which names `tx_frag_length` for the way the start is projected, is reworded; the projection itself is unchanged.

## 5. What it changes, measured on g50 κ0.99 capture-OFF (10 M templates)

The counts are unique, unspliced, inward primary pairs (`$RIGEL_RUNS/arms/2026-10-19_length_rule/`):

| pairs | mRNA | gDNA | nascent |
|---|---:|---:|---:|
| mates overlapping, an annotated intron inside the block (unchanged by the rule) | 2,231 | 807 | 420 |
| mates apart, an annotated intron in the span covered by more than K (changed) | 17,815 | 4,739 | 2,195 |
| RNA whose own transcript's intron is overrun by more than K (the rule's risk) | 0 | — | — |

- **Who moves.** For mRNA the covered intron is another transcript's, since a molecule's bases lie on its own
  transcript's exons. STAR overruns a molecule's own intron by 1 to 3 bases in 6,349 pairs, all within K, and never
  further.
- **Expected effect.**
  - Today a gDNA read through an intron gives the transcript it contradicts a shortened, plausible length. That
    stops, so less gDNA should be counted as RNA.
  - Isoforms whose intron an mRNA read covers lose the same advantage.
  - The candidates that move all carry an overhang penalty already, so the totals should move little.
- **The coverage weight.** For a contradicted candidate, the fragment's interval on the transcript grows by the
  retained intron's length (clipped to the transcript).
- **Not changed.** A deletion that stands for a micro-intron is read as covered sequence, under both rules. STAR
  writes gaps shorter than 21 bp as deletions, and annotated introns that short number 4 of 45,609 in the panel's
  index (all 1–3 bp, inside K) and 72 of 1.41 M in the human index. An intron of 1–3 bp inside a block, cut today, is
  kept.

## 6. Tests and breaks

**Failing on today's code:**

1. **Row 2 of the example** (`resolve_fragment`, K 3): mates 25 and 15 bases into T2's intron give T2 170, not 110.
2. **The enumeration, re-derived** (`test_transcript_space_fl.py`'s `TestLengthRuleByEnumeration`). Its base-by-base
   reference becomes the path rule: every base of the span counts except the observed junctions and the candidate's
   introns that lie in a gap within K. It runs at K 0 and K 3. Its mates are long enough to run into and through
   introns, the case the endpoint reference excluded by design (`_MATE = 20`, `:799`).
3. **One definition**: for every candidate of random fragments on the test index, the resolver's length equals
   `Accumulator.length_under` for the hypothesis that candidate implies (`bam_scanner.cpp:3040`).

**Passing today, kept as guards:** row 1 (one block, 170); row 3 (the intron in the gap: T2 110, the nascent span
170).

**Changed:** `test_internal_intronic_overhang_not_inflating` (mates 5 and 4 bases into the intron) is restated at K 3
as a covered intron, kept; its overruns within K (3 and 2 bases) as a cut intron.

**Breaks, each seen to fire:**

- cut covered introns again (the endpoint rule) → 1, 2;
- measure L by a second formula, the span less the sum of the introns → 3;
- give a single block anything but the span → the row-1 guard;
- drop K from the gap cut → 2's K 3 arm.

## 7. Measurement

- **The STAR panel against main**, pinned and fractional, per stratum, transcript and gene |error|, with the census by
  true origin, and read beside every capture efficiency at 1, since the reader amplifies every scan change.
- **The oracle ladder:** simulated alignments cover other isoforms' introns too, so a small move is expected; its size
  is recorded.
- **The pinned scan time**, since each candidate's path is collected where `tx_frag_length` took two binary searches.

## 8. At landing

- `DESIGN.md` §3.4: the EM's per-candidate length is the deposit's definition on the candidate's path, replacing "the
  span less every intron of the transcript lying strictly inside it".
- A CHANGELOG line.
- Roadmap item 11 ticked.
