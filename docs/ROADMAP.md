# ROADMAP — the state and the next steps

The tool's state in one paragraph, the claims that still hold with the instrument that re-derives each, and the ordered
next steps. Every open problem's substance is in its GitHub issue (https://github.com/mkiyer/rigel/issues); how
performance is judged is `SUCCESS.md`; the changelog is git.

## State

`src/`'s numbers are 1f97ff96's, the tree with the junction cap (2026-11-13, `DESIGN.md` §9.2) and the EM's gDNA share
over one tally (2026-11-14, §9.3), with the capture weight landed on 2026-11-15 (§7.2): the capture reading is raised to
the share of the library's reads whose own gDNA would show. That tree is the baseline every change is now read against.
It reproduces the prototype's rows exactly on the 94 rows of its screen: the 56 of `SUCCESS.md`, the test chromosome's
depth panels and the one-probe panel, pinned and fractional. Against 1f97ff96, transcript / gene error:
- At 5 % gDNA and above the weight is 0.98–1.00 under capture and 0.87–1.00 off it, and no in-scope whole-genome row
  rises by more than about 1 %: the largest rise is the simulator's stranded capture ON 5 %, genes +1.0 % (+518). Off
  capture a weight below 1 removes some of the shipped reader's false capture (STAR stranded capture OFF 5 %: genes
  −7.4 %). The deferred stratum's simulator g98 rises +1.0 % / +1.2 %.
- The hole between 0 and 1 % gDNA under capture closes on genes. The low-gDNA panel's stranded capture ON 0.1 % reads
  1,049,271 / 228,305 against the no-capture line's 625,202 / 232,081 (before 4,662,866 / 438,698): transcripts stay
  1.7× the line, which the release ruling reports and does not require. At 1 % the reading stands (weight 0.97).
- STAR's zero-gDNA rows fall: stranded capture OFF, the switch row, 1.90 M / 415 k → 0.78 M / 120 k (the line 0.20 M /
  27 k); stranded capture ON 1.19 M / 386 k → 0.67 M / 250 k; unstranded capture OFF 342 k / 147 k → 256 k / 55 k. The
  switch row stays about four times the line, because calibration's false gDNA counts as visible there
  (issue #1).
- The one-probe panel keeps its probe from 5 % gDNA up: weight 0.96–1.00, genes within 0.2 %.
- The test chromosome's transcripts move by up to 5 % with its genes and pools flat, which is
  issue #34. At 1/100 depth the weight falls even at 5 % gDNA (0.44 under capture),
  and the genes there improve (−38 %).
- On three real captured cfRNA libraries the weight is 0.25, 0.95 and 0.32. On the 93 %-gDNA one it calls some DNA RNA
  (issue #2).

The scan phase landed on 2026-11-18: the splice blacklist is retired (`DESIGN.md` §0b), a leading intron keeps every
later junction's ends, the splice census counts fragments, and a multimapper's candidates are pruned within each
place. Pinned and fractional, against the tree before it:
- On the simulator's alignments nothing moves: the four ladder identity references and the test chromosome's 90 rows
  are bit-identical.
- On the STAR rung only the pruning acts. Pools move by under 100 fragments on every row but stranded capture ON, and
  transcript error by −0.2 % to +1.0 % (it rises on 12 of 16 rows). At stranded capture ON the kept candidates link 19
  to 40 more loci into the one locus that holds most of the row (issue #78), where the EM stops at its cap
  (issue #34). What holds at a tighter stop is gDNA −435 (g05) and −693 (g50) fragments, mostly into nascent RNA.
- On the VCaP mix, whose BAM was aligned at STAR's defaults, the retirement moves 2,398 more DNA fragments into
  annotated RNA (issue #7); the pruning moves under 250 fragments of either library.

The intron abundance model landed on 2026-11-19 (`DESIGN.md` §5.4, `EQUATIONS.md` §7.1–7.2): every intron's count is
read against the gDNA law of the library's own intergenic regions, with a departure prior a stranded library measures
from its introns' strands. Pinned and fractional, against the tree before it:
- On unstranded capture OFF transfer now beats silent on every row. Calibration's Σ|gDNA error| on the ladder falls
  from 618,686 to 337,648, and the simulator's gDNA net from −325,530 to −77,872 fragments. Bar 1's worst in-scope row
  goes from −3.03 % to −1.76 % of the library's gDNA.
- The stranded strata hold, transcripts and genes within 0.7 %, and bar 5 holds and fails on the same rows.
- On the VCaP mix the introns' strands read the departure prior as 0.995, so the model steps aside there: every pool
  is within 79 fragments of the tree before it.
- STAR's unstranded capture OFF transcripts rise from 657,739 to 711,513, all of it g05's giant locus at the EM's cap
  (issue #34, issue #78); every other row of that stratum improves.
- Sparse unstranded libraries book about 1 % of their gDNA as nascent RNA (issue #82). The deferred stratum's annotated
  error rises from 430,814 to 520,047 on the simulator's alignments (issue #4).

All five identity references are re-frozen on this tree (`TESTING.md` §6).

The capture reader's first fixes landed the same day. The junction cap counts a transcript once per region where the
mask splits an exon (issue #61). An efficiency of exactly 1 gives RNA's own length, and only capture's excess over 1 is
carried between length tables (issue #65). A landscape trained on anchors alone reads one basin, gated by a test
(issue #6, which no longer reproduced). Off capture the identity references are bit-identical; under capture
transcript error moves by at most 0.35 % and gene error by 0.9 %, and bar 5 holds and fails on the same rows. The
capture-ON references and the VCaP mix's are re-frozen on that tree.

The validation of ec996bba is closed. On the 16-condition ladder, on the simulator's
alignments and on STAR's, HEAD beats the tagged 0.7.1 in the three in-scope strata: Brier 0.23–0.88 of 0.7.1's on
identical reads, and transcript error 0.10–0.47 of 0.7.1's. It also states eight failures, each an open entry
but those marked closed:
- unstranded sparse capture, where most true gDNA is called RNA in every build, and the deferred stratum (unstranded ×
  capture ON) (issue #4);
- capture outside the annotation's partition: with exact probe copies, gDNA called RNA rises ×2.7–2.8 at stranded
  capture ON, mostly away from the copies (issue #3);
- on STAR at g00, the capture tail, whose few uncertain gDNA bookings swing ×2 across draw seeds, and the switch row,
  where under 0.1 % of templates move transcript error ×9.5, worse than 0.7.1
  (issue #1). With the capture weight both beat 0.7.1 on the pinned seed (the
  switch row 0.78 M transcript errors against 1.38 M), and the switch row stays four times the no-capture line;
- two gene regressions against 0.7.1, both closed: STAR g98 κ0.99 OFF with the share over one tally (37,304 → 7,595,
  against 22,898), and the STAR g00 tail with the capture weight (386 k → 250 k against 328 k, one seed);
- VCaP exome DNA called RNA ×1.48 against 0.7.1 on identical reads (issue #5);
- annotation omission: unannotated unspliced RNA booked gDNA, and retained genes damaged ×1.4–3.0
  (issue #12).

Against 0.7.1 re-scored on the same reads (2026-11-16), HEAD's transcript and gene errors are lower on every in-scope
rung of the ladder, on both inputs. At 5 % gDNA and above, four stranded capture-ON rows call more of the library
annotated than 0.7.1 does: the simulator's g05, g50 and g98, and STAR's g50 (issue #76). The release holds for them
(`DESIGN.md` §0b).

No preview ships (the owner's ruling, 2026-11-04): 0.8.0 waits for every `release-blocker` issue, in the order under
Next. The validation's receipts are in `$RIGEL_RUNS/prototypes/2026-11-02_preview_validation/`, and the re-score's in
`$RIGEL_RUNS/prototypes/2026-11-16_v071_rescore/`.

## Claims that hold

| claim | re-derived by |
|---|---|
| HEAD beats 0.7.1 on transcripts and genes on the ladder's three in-scope strata, on both inputs | `quant_accuracy.py` per stratum, pinned and fractional, `--bam` for the STAR rung, beside the tag |
| HEAD's source probabilities beat 0.7.1's by the proper score on identical reads, and most of the exchange is compatible with approximately calibrated source uncertainty at the resolution tested | the validation's `paired.py`, outside the tree (issue #36) |
| Scarce gDNA at 1/10 and 1/100 depth shows no blow-up, and dense probe panels hold | `quant_accuracy.py` on the test chromosome's depth and probe panels |
| Calibration's library gDNA is near truth on the in-scope strata; its per-locus count errs both ways, mostly its own error rather than the pooled share | `calibration_vs_oracle.py`, `prior_vs_oracle.py` |
| Under capture the capture reader carries most of the accuracy, and at STAR g00 most of the error | `quant_accuracy.py --bam` beside `--arm no_capture`, `ruler_vs_truth.py` |
| The transfer policy wins on unstranded × capture ON and the contaminated stranded rows, and loses to silent on unstranded × capture OFF | `policy_benchmark.py --panel ladder` |
| A refactor or speed-up is a numeric no-op | `rename_identity.py --check` against the frozen reference |
| Every instrument runs and regenerates | `preflight.py --full` |

## Next

The `0.8.0` milestone (https://github.com/mkiyer/rigel/milestone/1) is what the release fixes or states; its
`release-blocker` label is what it fixes. In order:

1. **The capture reader, bounded at two sessions, with the owner:** issue #1 against a pass bar the owner sets before
   any candidate runs. If no candidate passes, issue #1 ships stated.
2. **The annotated pool against 0.7.1:** issue #9 and issue #76, against `SUCCESS.md`'s bar 5. Dissect the largest
   gap first, the simulator's stranded capture-ON g50, and bring the dissection to the owner before fixing.
3. **The release:** one ladder re-score; a MANUAL limitation for issue #2, issue #3, issue #4 with issue #5, issue #7,
   issue #10, issue #12, issue #13, issue #23, issue #46, issue #78 and issue #82, and for issue #1 if it falls back;
   the CHANGELOG; then issue #11.

After that, the open issues outside the milestone, in issue order.
