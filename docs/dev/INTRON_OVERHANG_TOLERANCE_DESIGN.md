# The intron-overhang tolerance K, use by use: design (roadmap item 14)

2026-10-20 · code base `5482a10c` · for the owner's review. Every claim about Rigel's code cites a line re-read for this
document. Every number is measured; the receipts are in `$RIGEL_RUNS/arms/2026-10-20_k_uses/`. Approved 2026-10-21
(owner: keep both sites; the outer end's own arm run before the issue closes; one hold test only if it proves simpler
and better code); the arm ran and steps 2 and 3 landed the same day, their measurements the roadmap's item 14
(`PAIRED_END_ROADMAP.md`).

## 1. The recommendation in brief

K (`--intron-overhang-tolerance`, default 3) lets an aligned block's end lie up to K bases inside an annotated intron.
The code reads it in four tests. Measured, it acts at two places in a fragment, and on one rare annotation quirk:

1. **The edge of a gap the reads did not sequence.** A mate runs a base or two into the intron that spans the gap.
   - This is implicit splicing, the owner's scope.
   - **Keep.** It gives the molecule's own transcript its true length in 17.6–35.2 k fragments per 10 M templates.
     K = 0 never gives a truer one.
   - It also keeps 12.5–25.4 k pairs per 10 M that the pair rule would otherwise withhold, 99.0–99.3 % of them true
     molecules at their true length.
2. **The fragment's outer end.** A mate starts or ends 1–3 bases inside an intron of its transcript.
   - This is the one use outside implicit splicing.
   - **Keep: it is justified.** It is the same aligner habit: running a read's end across a splice site rather than
     splicing a short anchor.
   - Refused, it would withhold 5.1–10.2 k true pairs per 10 M, 99.3 % of them at exactly their true length.
   - It would also drop the molecule's own transcript from calibration's hypotheses in 8.4–16.7 k fragments.
   - Its own measured arm (K = 0 at the outer end alone) confirmed it (§5, step 1).

The annotation quirk: the real annotation has 46 "introns" of 1–3 bases (26 of them a single base), which a read simply
sequences through. The hold slack keeps such a transcript, rightly. This happens 1–440 times per real library and never
on the panel; on cfRNA-C all 34 cases are a 1-base "intron" inside a read's block. K acts nowhere else. It never acts
beside an observed junction, which is matched exactly.

Two fixes, neither a change to the rule:

- **Three descriptions omit the pair rule's use:** the config docstring, the CLI help and the MANUAL row.
- **The tests:** the outer end has none, and nothing pins that K never matches an observed junction.

The ISSUES entry also asked about the cut "moving the splice point to the annotation". It is right. The bases a read
runs into the intron are the next exon's, so cutting the intron whole gives the true length (§4).

## 2. What K is (the owner's rulings)

- **2026-10-16.** K "belongs to implicitly spliced fragments": mates aligned to different exons with a whole intron
  between them, where up to K aligned bases may sit inside the intron. The owner keeps it there, "but we must be
  judicious about deciding which cases apply it" (`ISSUES: the-anchor-tolerance-applies-beyond-implicit-splicing`).
- **2026-10-17, the wave-1 plan §7.5.**
  - A use is classified by the gap it tests, not by the hit's splice label.
  - Each use is tested at 0, K and K + 1, once with an implied gap and once beside an observed splice.
  - No use may shift an observed junction, delete alignment evidence or a candidate, or forgive two observed skips
    that conflict.
  - K is renamed `intron_overhang_tolerance`.
- **2026-10-19, the pair rule.** A witness transcript overlaps both mates, holds every observed junction, "and no
  aligned base sits more than K inside one of its introns".

## 3. Where K is read

| test | where | what it decides | what it reaches |
|---|---|---|---|
| **W-cut** | `validate_pair` (`resolve_context.h:906`) → `collect_transcript_introns_in_gap` (`:1055-1081`) | a witness's intron in an unaligned stretch of the molecule is cut whole when it lies within the stretch widened by K | the pair's length along the witness, so usable or withheld `length` |
| **W-hold** | `validate_pair` (`:896-902`, K at `:901`) | a witness offers no path when one of its introns in the molecule, not observed, overlaps an aligned block by more than K | the same |
| **G-cut** | `enumerate_gap_hypotheses` (`:1033`) → `collect_transcript_introns_in_gap` | a candidate's intron in an unsequenced gap is cut whole within the gap widened by K; gaps equal to an observed junction are not gaps (`:998-1009`) | the candidate's path: its length for the EM (`:1034`), its hypothesis for calibration (`:1040`), and the SPLICE_IMPLICIT label (`:1454`) |
| **G-hold** | `enumerate_gap_hypotheses` (`:1035`) → `transcript_contradicts_blocks` (`:688-710`) | a candidate one of whose introns overlaps an aligned block by more than K joins no hypothesis; it keeps its length | calibration's hypotheses, and so the label |

- The enumeration runs for every fragment with an unsequenced gap, whatever its splice label (`:968-1010`). So a
  spliced hit's mate gap is searched too: an implicit splice beside an observed one.
- K is read nowhere else:
  - not by junction matching, which is exact (the annotated-path rule);
  - not by scoring, whose overhang penalty prices a candidate's out-of-exon bases by its own alphas;
  - not by calibration's accumulator, which arbitrates the hypotheses it is given;
  - not by the splice blacklist's anchor test.

## 4. Where it acts, measured

**Method** (`k_uses.py`):

1. Each template's primary pair goes through a twin of the current pair rule and annotated-path rule, once with each
   slack removed.
2. The kept pairs with a gap are resolved by the native resolver at K = 3 and at K = 0.
3. Every difference is attributed to the tests whose slack it exercised. No difference is unexplained: 0 twin
   mismatches in 50 M panel templates and 21.3 M real ones.
4. Truth comes from the read name: the molecule's transcript and genomic blocks. A use is right when it gives the
   molecule's own transcript its true length.

**On the STAR panel** (recommended parameters), per 10 M templates, for g00 / g05 / g50 κ0.50 / g50 κ0.99 OFF / ON:

| | g00 | g05 | g50 κ0.50 | g50 κ0.99 OFF | g50 κ0.99 ON |
|---|---:|---:|---:|---:|---:|
| kept pairs with an unsequenced gap | 5.31 M | 5.27 M | 4.84 M | 4.84 M | 6.38 M |
| resolutions K changes | 58.0 k | 55.3 k | 28.9 k | 29.0 k | 41.8 k |
| **the gap's edge** (G-cut, with G-hold at the same intron) | 41.4 k | 39.7 k | 20.6 k | 20.7 k | 29.5 k |
| — mRNA: K gives its transcript the true length | 35.2 k | 34.1 k | 17.6 k | 17.6 k | 25.1 k |
| — mRNA: short by exactly the fragment's outer soft clips | 0.55 k | 0.52 k | 0.28 k | 0.29 k | 0.36 k |
| — mRNA: the slack serves only another isoform | 4.0 k | 3.5 k | 1.8 k | 2.0 k | 2.2 k |
| — gDNA: a spliced hypothesis added (beside the unspliced one) | 0 | 4 | 22 | 24 | 405 |
| — K = 0 gives the true length where K does not | 0 | 0 | 0 | 0 | 0 |
| **the fragment's outer end** (G-hold) | 17.6 k | 16.5 k | 8.8 k | 8.8 k | 13.0 k |
| — mRNA: its transcript keeps its hypothesis at its length | 15.9 k | 15.0 k | 8.0 k | 8.0 k | 11.6 k |
| — the label turns implicit (the transcript's exact gap intron, which K = 0 vetoes) | 5.8 k | 5.5 k | 3.0 k | 2.9 k | 4.8 k |
| — mRNA: another isoform kept | 0.8 k | 0.7 k | 0.4 k | 0.4 k | 0.7 k |
| — gDNA: a spliced hypothesis added | 0 | 0 | 8 | 11 | 132 |
| **pairs kept only by the slack** (the pair rule) | | | | | |
| — the gap's edge, W-cut / W-hold | 23.5 / 25.4 k | 22.5 / 24.3 k | 11.6 / 12.5 k | 11.6 / 12.6 k | 16.6 / 18.1 k |
| — the outer end, W-hold | 10.2 k | 9.5 k | 5.1 k | 5.2 k | 7.7 k |

Of the pairs the slack keeps, 99.0–99.3 % in every condition are mRNA whose own transcript gives their exact length,
at the gap's edge and at the outer end alike; in the rest the transcript is no witness and another one is.

**How far a true molecule's reads run into its own transcript's introns, whatever K** (mRNA, kept pairs; g00, g50
κ0.50 within 0.4 points of it):

| site | cases | 1 bp | 2 bp | 3 bp | 4 bp | 5 bp or more |
|---|---:|---:|---:|---:|---:|---:|
| a gap's edge | 37.8 k | 71.5 % | 28.0 % | 0.1 % | 0.0 % | 0.3 % |
| the fragment's outer end | 82.8 k | 57.1 % | 28.2 % | 9.9 % | 2.4 % | 2.4 % |
| anywhere else in the fragment | 1.5 k | 0 | 0 | 29.0 % | 10.3 % | 60.7 % |

- K = 3 covers 99.6 % of the gap edge's overruns and 95.2 % of the outer end's.
- At a gap's edge a read almost never runs 3 bases in. STAR at its default splices an annotated junction from a 3-base
  anchor (`--alignSJDBoverhangMin 3`, which the panel leaves unset), and extends or clips a shorter one.
- Overruns anywhere else are other misalignments, mostly 10 or more bases. K does not forgive them, and should not.

**On the real libraries** (no truth):

| | cfRNA-A | cfRNA-B | cfRNA-C | `vcap_mix` |
|---|---:|---:|---:|---:|
| kept pairs with an unsequenced gap | 21.9 k | 9.0 k | 87.6 k | 8.83 M |
| resolutions K changes | 571 | 273 | 2,417 | 136.0 k |
| the gap's edge (G-cut) | 408 | 145 | 1,834 | 118.7 k |
| the outer end (G-hold) | 179 | 135 | 656 | 19.6 k |
| pairs kept only by the slack: gap's edge (W-cut / W-hold), outer end | 187 / 210, 116 | 64 / 76, 95 | 595 / 681, 369 | 55.0 / 57.3 k, 9.8 k |
| a 1-base annotated "intron" a read sequences through (occurrences) | 2 | 1 | 34 | 440 |

**Why cutting the intron whole is right.** A mate that runs d bases past a donor read the next exon's first d bases,
which happen to match the intron's first d bases. The molecule is the left mate up to the donor plus everything from
the acceptor. Its length is the extent less the whole intron, which is what K's cut gives. It is exact in every
case where it changes the transcript's length: 98.5 % give the true length, and in the rest the shortfall equals the
fragment's outer soft clips, which no rule here sees (checked case by case on 2 M templates of g00).

## 5. The plan

1. **The outer end's own measured arm.** This is the issue's requirement for a use outside implicit splicing.
   - A prototype applies the hold slack at a gap's edge only: K for an intron the path cuts, 0 for any other.
   - It is run against main on the STAR panel, pinned and fractional, per stratum, with every capture efficiency at 1
     beside the shipped run, and the census.
   - Expected from §4: true pairs withheld and true transcripts' hypotheses lost.
   - Its numbers close the issue as justified. Nothing lands from the arm itself.
   - **Run 2026-10-21** (`$RIGEL_RUNS/arms/2026-10-21_k_outer_end/`). With every capture efficiency at 1, transcript /
     gene error +1.19 / +13.35 % unstranded capture-OFF, +1.51 / +20.12 % stranded capture-OFF, −0.16 / −0.65 % stranded
     capture-ON; the census: 4.8–9.6 k more fragments withheld per condition. Shipped: −0.94 / +3.57 %, −1.05 / +4.26 %,
     −0.63 / +2.58 %, the transcript gain the capture reader's (`ISSUES: the-reader-reads-false-gdna-as-capture`). The
     oracle ladder: all 16 conditions move, within 0.37 % in every stratum.
2. **The record, one commit, no change to the rule.**
   - The config docstring, the CLI help and the MANUAL row say that K also decides, along a witness transcript, whether
     two mates are one molecule, and that it acts at a gap's edge and a fragment's outer end.
   - The tests of §6.
   - `ISSUES: the-anchor-tolerance-applies-beyond-implicit-splicing` closed with §4's numbers and the arm's.
3. **Optional, its own step: one hold test.**
   - `validate_pair`'s inline test (`:896-902`) re-implements `transcript_contradicts_blocks` (`:688-710`). It skips
     observed introns and introns outside the molecule, and tests each mate's blocks rather than the merged ones.
   - One function for both callers is a simplification.
   - It is not a numeric no-op where two overlapping mates each overlap an intron by K or less but together by more, so
     it needs its identity check and, if anything moves, its own A/B.
   - **Landed 2026-10-21.** The witness passes its mates' blocks unmerged, as its inline test read them, so the case of
     two overlapping mates does not arise. The inline test's two skips go, both unreachable: an intron outside the
     molecule overlaps no aligned block; and an observed junction can hold aligned bases only of the other mate, whose
     span then overlaps this one's, so the junction covers what lies between them and the pair, at most its reads' bases
     and clips along the genome, passes before any witness is asked. Ten lines of code fewer. Bit-identical on oracle
     g00, STAR g50 ss.50 OFF, cfRNA-A and cfRNA-C (frozen on `5482a10c`); the shared test at K + 1 and at K − 1 fires
     both callers' outer-end cases (`breaks_hold3/`); the pinned scan of STAR g50 ss.50 OFF 38.8 → 38.6 s (A B A B).

## 6. Tests and breaks

**Today's gates** (`tests/test_implicit_splice.py` a–m and U1–U4; `test_hit_verdict.py`
`test_the_tolerance_lets_a_mate_run_k_bases_into_an_intron`; `test_transcript_space_fl.py`
`test_an_overrun_within_K_is_cut`) cover the gap's edge: in the resolver at 0, K and K + 1, and in the pair rule at K
and K − 1 bases of slack.

**New** (`test_hit_verdict.py`, 14 cases):

- **The outer end, in the resolver:** a mate starting K bases inside its transcript's intron leaves the transcript a
  hypothesis (and its gap intron found); at K + 1 it is contradicted; at 0 no slack is used.
- **The outer end, in the pair rule:** a long pair whose mate's outer end runs K bases into a witness's intron is one
  molecule along the witness; at K + 1 it is withheld `length`.
- **A gap's edge beside an observed splice,** in each: a spliced mate runs 0, K and K + 1 bases into the intron that
  spans the mate gap. Within K the intron is cut whole (the witness's path; the candidate's length), at K + 1 it is not.
- **Beside an observed junction:** a read whose junction lies 1 base from the transcript's intron is a novel junction
  at any K. K never matches a junction.

**Breaks, each seen to fire** (`$RIGEL_RUNS/arms/2026-10-21_k_outer_end/breaks_record/`):

| break | fires |
|---|---|
| the outer end's hold slack at 0 (the arm) | both outer-end cases at K |
| the witness's hold slack at K + 1 | its outer-end case at K + 1 |
| the witness's hold slack at K − 1 | its outer-end and beside-a-splice cases at K, and the older pair-rule gate |
| the resolver's hold slack at K + 1, at K − 1 | its outer-end case at K + 1, at K |
| the cut window at K + 1 | the resolver's beside-a-splice case at K + 1, and three length-rule enumerations |
| the cut window at K − 1 | both beside-a-splice cases at K, and six older gates |
| junction matching widened by K | both junction cases |

The witness's beside-a-splice case at K + 1 fires under no single break, by construction. At a gap's edge the hold test
and the cut window ask one question, whether the block runs at most K bases into the intron, and the witness holds
before it cuts, so widening either one alone changes no verdict: the other still refuses. The resolver records a
candidate's length before its hold test, so its cut window is seen alone.

## 7. Measurement

- **The arm (step 1):** the STAR panel against main, pinned and fractional, per stratum, with every capture
  efficiency at 1 beside the shipped run, and the census by true origin.
- **The oracle ladder:** a gDNA read placed at truth can end a few bases inside an intron by chance, so the arm runs on
  it too.
- **Step 2** changes no number: the suite, preflight, and rename_identity bit-identical on oracle g00 and STAR g50 ss.50
  OFF.
