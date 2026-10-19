# Geometric read pairing: design

2026-10-19 · code base `9cb85fae` · revision 2. The owner approved revision 1 and asked why the search does not pair by
shared transcripts. It should, for spliced paths: §5 adopts the owner's rule and withdraws revision 1's span bound.
Status: step 1 (§5) landed 2026-10-19; its numbers are in `ISSUES: the-scanner-builds-fragments-from-one-mate`.
Steps 2–4 landed 2026-10-19: the owner took the search for every input, HI not read (`DESIGN.md` §1.1; the four arms in
CLOSED `ISSUES: star-hi-pairs-against-rigels-own-pairing`). The arms here were named `hi`, `hi-collapsed`, `search`,
`legacy`. Step 5, `SPLICED_UNANNOT` deleted, is next.
Every claim about today's code cites a line re-read for this document. Roadmap Phase A items 3, 5 (in part), 6 and 8;
plan §6–§7.1.

## 1. The change in brief

Rigel pairs a template's mates in two ways today. With HI tags (STAR), each HI group is one hit. Without them
(minimap2), a heuristic pairs the mates' alignments by shared transcripts. This change replaces the heuristic with
**the search**: every R1 alignment is tried with every R2 alignment on the same reference with the opposite strand,
and the hit verdict that already judges STAR's pairs keeps the ones that can be one molecule.

That verdict follows the owner's rule (2026-10-19):

- an **unspliced** pair needs its mates within `max_frag_length` of each other;
- a **spliced** pair, an implicit splice included, needs **one transcript compatible with both mates**:
  - the transcript overlaps both mates;
  - it holds every observed junction;
  - no aligned base lies more than K bases inside one of its introns.

The same search then answers the open HI question (`ISSUES: star-hi-pairs-against-rigels-own-pairing`): run it on the
STAR panel with HI ignored and compare it end to end with HI pairing. The owner's rule (2026-10-19): if HI-agnostic
pairing gives similar or better results, it becomes the one approach for every aligner.

## 2. Pairing, taught by one template

A template is one read pair, R1 and R2. An aligner may place each mate in several places. Say R1 aligns at X1 and X2,
and R2 at Y1 and Y2. There are four ways to make a molecule: X1–Y1, X1–Y2, X2–Y1, X2–Y2.

- **STAR** reports pair alignments: HI 1 = (X1, Y1), HI 2 = (X2, Y2). It never offers X1–Y2 or X2–Y1, even when one of
  them is the true molecule. Each HI group is a whole pair, so STAR repeats a mate's record in every pair that uses it:
  a template whose R1 has one place and R2 two writes R1's record twice.
- **minimap2** reports each mate's places on its own and leaves the pairing to the reader.
- **The search** tries all four combinations and keeps those the owner's rule accepts (§1, §5). This is the pair
  validator that already judges STAR's pairs (`resolve_context.h:799-886`), with §5's change. The annotated-path rule
  then applies to each kept pair as it does today (`bam_scanner.cpp:999-1009`).

Two details the rule needs beyond its three conditions, both already in the validator:

- **The mates must face each other.** The forward mate lies upstream of the reverse mate, each 5′ soft clip walked
  back along the genome first, since a clip can hide the true start.
- **The length holds on the spliced path too.** Mates on exon 2 and exon 30 of one transcript are compatible with it,
  but kilobases apart along it. So the transcript must also give a path of 1 to `max_frag_length` bases between the
  mates, its introns there cut whole by the K rule. On the STAR panel, 183 k long pairs share a transcript and have no
  such path; none is correct (§5).

Measured at the template level on 47.9 M STAR templates (2026-10-18, per-base truth):

- the validator rejects no correct pair;
- wherever a correct pair exists among a template's places, the search finds it;
- STAR's HI groups miss it in 10,108 of 3.22 M multi-place templates (0.31 %);
- the search adds 475 k combinations STAR did not report, 2 % of them correct.

What is still unknown is whether those extra combinations help or hurt the counts, end to end. That is §7.

## 3. Today's pairing without HI, and why shared transcripts failed there

Without HI (`bam_scanner.cpp:898-949`):

- each secondary record is a location of its own;
- a secondary whose mate is unmapped becomes a one-mate hit;
- the primaries, with any supplementary record, are one more location each.

Then `pair_multimapper_reads` (`bam_scanner.cpp:1010-1147`):

1. resolves every location as a fragment of one mate (`_resolve_core`, `:1033-1046`);
2. pairs every R1 and R2 location whose transcript sets intersect (`:1067-1082`);
3. pairs the locations left over, nearest first on one reference (`:1084-1112`);
4. pairs every location still left with every other (`:1114-1136`);
5. makes one-mate hits of whatever remains (`:1138-1144`).

Sharing a transcript is the right test for a spliced path (§5). Step 2 fails because it uses that test as the
pairing itself:

- **The nascent span counts as a shared transcript.** The resolver adds a gene's synthetic nascent span to the
  candidates of every fragment that touches the gene (`resolve_context.h:1279-1328`). So any two places in one gene
  "share", however far apart they are, and step 2 pairs them.
- **A matched location is never tried again.** A location matched in step 2 is left out of steps 3–5. So when its
  true partner lies outside the gene, they are never tried together.
  - Example: R1 in a gene's last exon, its true R2 just past the annotated 3′ end, and a second R2 place in an Alu
    inside the same gene's intron. Step 2 pairs R1 with the Alu. The validator rejects that pair (it faces outward),
    the true R2 becomes a one-mate hit, and the template is withheld.
- **Unspliced pairs need no transcript.** Intergenic mates share none and fall through to step 3, which is greedy:
  one partner each.

Two more defects, apart from transcripts:

- **Every record is a location.** A STAR BAM read without HI repeats a mate's record in each pair alignment. Each
  repeat pairs again, so the template gets duplicate hits and its `num_hits` rises.
- **It resolves a single mate before any verdict.** This is the last place where a single mate is resolved, and the
  last thing keeping the `SPLICED_UNANNOT` label alive (roadmap item 10). It is also the slow part: the full resolver
  runs on every location before any pair is formed.

## 4. The search

**Mate locations** (plan §6, its first rule). Each mapped, non-supplementary record with at least one aligned base is
an alignment of its mate. Records of the same alignment are merged into one location; the others are its aliases.
The identity is:

- the mate;
- the BAM reference;
- the strand;
- the leftmost aligned position;
- the CIGAR, with M, = and X read as one operation (adjacent runs merged), and soft and hard clips kept.

So STAR's repeats of one record merge. Two different alignments of one mate at one place stay apart, for example
`50M3D50M` and `50M3N50M` (a deletion and an intron), or a spliced and a contiguous reading: those are the competing
readings the scan keeps. The representative is the primary record when an alias is primary, else the first in record
order.

To support this, two fields join `ParsedAlignment`: the canonical CIGAR, and the BAM's reference id. Today only the
index's reference id is kept (`bam_scanner.cpp:725,751`), and every reference outside the index shares −1. The
annotation pass parses records with its own copy of the parser (`bam_scanner.cpp:2919-2958`); it calls
`parse_bam_record` instead, so the two new fields have one source.

**Supplementary records** stay with the same mate's primary location, as today without HI (`bam_scanner.cpp:903-908`).
Bundle item 13 decides their fate.

**Candidates.** Every R1 location × every R2 location on the same BAM reference with opposite strands. Each candidate
goes, unchanged, through:

1. the hit verdict: missing mate, reference, strand, outward, length (with §5's witness);
2. the annotated-path rule (with §5's holder).

All usable candidates are the template's hits, and `num_hits` is their number. There is no cap and no early stop
(plan §7.1: a nearest-first stop misses a farther spliced partner).

**No candidate.**

- If one mate has no location, each location of the other becomes a one-mate hit, rejected as `missing_mate`, as today.
- If both have locations but none on one reference with opposite strands, the first location of each, in record
  order, is offered as the template's one hit. The verdict rejects it for reference or strand. So a unique pair is
  counted exactly as today, and a multimapper's reason once, not once per combination.

**Order.** Each R1 location in record order, with each R2 location in record order, so a run is deterministic.

**The mate counters.** Today `proper_pair`, `improper_pair` and `mate_unmapped` count over hits
(`bam_scanner.cpp:1936-1951`). Under the search they are counted once per template, from its primary records.

**The annotated BAM.** An alias record is written with its location's first hit. Today the end-of-group sweep would
write it as an intergenic pass-through (`bam_scanner.cpp:2847-2869`).

## 5. The witness: one transcript shared by both mates (owner, 2026-10-19)

The validator has two of the owner's three conditions: every observed junction held by one transcript, and the K rule
on aligned bases inside its introns. It lacks the first, that the transcript overlaps both mates.

- **An implicit splice.** For an unspliced pair longer than `max_frag_length` along the genome, the witness is any
  transcript whose span overlaps the molecule anywhere (`resolve_context.h:731-745`, called at `:835-837`). So a gene
  lying wholly between the two mates can shorten their path: a path that leaves that transcript at both ends.
- **A spliced pair.** One transcript must hold every junction (`resolve_context.h:753-784`), but it need not overlap
  the other mate.

**The change.** A witness, and a spliced pair's junction holder, must overlap both mates. A transcript overlaps a
mate when its span, from its first exon's start to its last exon's end, overlaps the mate's aligned span. "Overlapping
the boundaries", read as overlap, not containment: a mate that runs past an annotated end still overlaps, and its
extra bases are the overhang the EM already prices. A spliced pair whose holders all miss one mate is `unsupported`:
no one transcript explains the whole hit.

**Measured** with HI ignored, every pair of distinct locations on one reference with opposite strands, the validator's
own path rule applied under both witness rules, and each pair they judge differently scored against per-base truth
(`$RIGEL_RUNS/arms/2026-10-19_geometric_pairing/witness_rule.py`; five STAR conditions, 50 M templates):

| | pairs | share |
|---|---:|---:|
| all pairs | 76,349,702 | 100 % |
| facing outward: rejected, no index lookup | 10,640,182 | 14 % |
| 1..`max_frag_length` along the genome: accepted, no witness needed | 52,087,144 | 68 % |
| longer along the genome: a witness needed | 13,622,376 | 18 % |

| the long pairs | today's witness | the owner's |
|---|---:|---:|
| candidate transcripts walked | 1,970,709,473 | 27,781,708 |
| pairs with any candidate | 7,158,206 | 2,994,689 |
| pairs accepted | 2,812,090 | 2,811,226 |

- **Rejected by the owner's rule alone:** 864 long pairs, 0 correct. Also 1,643 spliced pairs short along the
  genome, whose holders all miss the other mate, 0 correct.
- **Sharing is not enough.** 183 k long pairs share a transcript and have no path of valid length along it, or an
  aligned base inside one of its introns. None is correct: every correct pair passes the validator (§2).
- **Cost.** The owner's rule walks 71× fewer transcripts. Each location's transcripts are looked up once per template;
  a pair intersects two small sets.

This replaces revision 1's span bound (`max_frag_length` + the reference's longest intron total). That bound cut only
the pairs more than about a megabase apart; the 10.3 M long pairs within a megabase would still have queried every
transcript between their mates.

**The work, whole.** On g50 ss.50 capture-OFF the search judges 18.3 M pairs in 10 M templates; the scan already
judges 11.1 M STAR pairs there today. Over the five conditions, 82 % of pairs are decided with no index lookup:
facing outward, or within `max_frag_length` along the genome. Today's pairing without HI runs the full resolver on
every location first (§3). If an aligner ever reports hundreds of places per mate, the loop
can become two joins with the same verdicts: by position for the unspliced test, by transcript for the spliced one.
Built only if the scan's time says so.

## 6. Two steps, one mechanism each

§5 changes the validator, which also judges STAR's HI pairs. So it lands first, on its own, measured against main on
the STAR panel with HI as shipped. On HI input it moves at most the few pairs §5 counts. The search follows on the
tree that has it.

## 7. The decision: four arms

One worktree build. An environment switch in the prototype only, never in `src/`, picks the pairing. Each arm runs on
the five STAR conditions, pinned and fractional, scored against truth.

| arm | pairing | what it is |
|---|---|---|
| `hi` | HI groups; the search without HI | ships for STAR input; must be bit-identical to the tree after §6 |
| `hi-collapsed` | HI groups, duplicate groups collapsed; the search without HI | STAR's pairs, deduplicated |
| `search` | the search, HI ignored | the HI-agnostic candidate |
| `legacy` | today's pairing without HI, HI ignored | what ships for input without HI |

The comparisons, each varying one thing:

- **`hi-collapsed` against `search` is the owner's question.** Both arms hold the same pairs except the combinations
  STAR did not report. Every valid HI pair is a search candidate, and the search's locations collapse STAR's duplicates
  as `hi-collapsed` does. So the
  difference is exactly the extra combinations: 475 k on the panel, 2 % correct, and among them the correct pair in
  10,108 templates where HI missed it.
- **`legacy` against `search` is the repair against what ships**, on input without HI.
  - ⚠ `legacy` on a STAR BAM also carries the duplicate-record defect (§3), which minimap2 input does not.
- **`hi` against `hi-collapsed`** re-measures the held duplicate collapse on today's tree. In g00 it moved SIK1B by
  27 k fragments (`ISSUES: the-em-stops-before-it-converges`). Keeping `hi-collapsed` separate keeps that paralog noise
  out of the pairing answer.

**What is read**, per stratum and never pooled:

- transcript and gene |error|;
- the census by true origin (true gDNA counted as RNA; RNA lost);
- the hit and withheld counters;
- scan seconds, from back-to-back pinned runs;
- the oracle ladder, identical to main in every arm: its templates are unique and carry no HI, and the search forms
  the one pair today's code forms.

Wherever `hi-collapsed` and `search` differ, the difference is dissected by gene. Exact paralogs are listed apart. A difference sitting
there is re-run at a 20× EM budget, to tell the stopping point from the converged answer. VCaP, a deep real library,
gives the cost check: the search's scan time and counters against `hi`'s.

**"Similar or better"** is the owner's call on those numbers. I will not set a bar in advance.

## 8. What lands afterwards

The search replaces `pair_multimapper_reads` in either outcome; it is needed for input without HI.

- **If the owner takes the search for every input:** `src/` pairs every template by the search. The HI grouping, its
  counters (`n_nh_above_reported`) and the mixed-template question (plan §7.1) are deleted. HI is no longer read.
- **If HI wins:** `src/` keeps HI groups for templates that carry HI (`hi` or `hi-collapsed`, as the owner rules on the duplicate
  collapse) and the search for templates without it.

Then, in a separate step, `SPLICED_UNANNOT` is deleted. Once no single mate is resolved before its verdict, every
resolved hit has its junctions held by one transcript (roadmap item 10).

## 9. Tests and breaks

New gates in `tests/test_hit_verdict.py`, without HI unless stated.

For §5, both fail on today's code:

1. **A gene between the mates is no witness.** Two unspliced mates farther apart than `max_frag_length`, with a
   transcript lying wholly between them whose introns would shorten the path. Today: a fragment. Now: withheld as
   too long.
2. **A spliced pair needs a holder over both mates.** F spliced on T's junction, R unspliced past T's annotated end,
   short along the genome. Today: usable. Now: withheld as `unsupported`.

The existing witness gates keep passing: the tolerance case (tK overlaps both mates) and the hybrid case (T3 holds
the junction and overlaps both mates).

For the search, the first three fail on today's code:

3. **The partner outside the gene** (§3's example): the true pair is found, and the template is a fragment, not
   withheld.
4. **A repeated record is one location:** a mate's record repeated twice, as STAR writes it, gives one location, so
   the fragment has two hits, not four.
5. **A location matched once is still tried with its other partners:** R1 in a gene, R2 there and again just past the
   gene's end, both valid: two hits. Today the gene pair matches by transcripts and the other R2 becomes a one-mate
   hit.
6. **A unique pair without HI is unchanged**, and **a mate with no location** gives `missing_mate`, as today.
7. **Two references only** (an index holding both): withheld as unpaired, counted once as `reference`. With the test
   index's one reference, the other mate is off the index and passes through unjudged (§10).
8. **An alias in the annotated BAM** carries its location's labels, not an intergenic pass-through.

Existing test whose premise changes: `test_a_lone_secondary_without_hi_is_no_hit`. Its lone R1 secondary is now an R1
location. The search tries it with R2 and rejects it as outward, so the primary pair still stays calibration's one
placement. The assertion moves from `missing_mate` to `outward`.

Breaks, each seen to fire:

- query the molecule's span for the witness again → gate 1;
- drop the both-mates check on holders → gate 2;
- restore transcript-first pairing → gates 3 and 5;
- drop the alias merge → gate 4;
- drop the alias write, so an alias falls to the end-of-group sweep → gate 8.

Hit order has no test, on purpose: it is not a correctness claim. It changes floating-point sums, so it is read in the
A/B numbers.

## 10. Out of scope

- References outside the index: the identity keeps them apart, and the verdict still passes them through unjudged
  (`bam_scanner.cpp:985`; bundle item 12).
- Supplementary chains and hard-clipped 5′ ends (item 13).
- Clip walks along a transcript path (item 9).
- NM conflicts between aliases (item 12): an alias is assumed to carry its location's NM.
- The duplicate-collapse ruling, except as the `hi-collapsed` arm measures it.

## 11. The steps

1. §5 alone, in a worktree: gates 1–2 and their breaks; against main on the STAR panel (HI as shipped) and the
   ladder; land.
2. The search, in a worktree: locations, the alias writes, the switch; gates 3–8 and their breaks.
3. The four arms, the ladder, VCaP; a receipt in the arm; the numbers to the owner for the HI ruling.
4. `src/` per the ruling; suite and preflight; commit, push, record.
5. `SPLICED_UNANNOT` deleted.
