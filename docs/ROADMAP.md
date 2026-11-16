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

The identity references are re-frozen on this tree (`TESTING.md` §6).

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

No preview ships (the owner's ruling, 2026-11-04): a pre-release waits for the capture reader and a review of the open
issues. The validation's receipts are in `$RIGEL_RUNS/prototypes/2026-11-02_preview_validation/`.

## Claims that hold

| claim | re-derived by |
|---|---|
| HEAD beats 0.7.1 on the ladder's three in-scope strata, on both inputs | `quant_accuracy.py` per stratum, pinned and fractional, `--bam` for the STAR rung, beside the tag |
| HEAD's source probabilities beat 0.7.1's by the proper score on identical reads, and most of the exchange is compatible with approximately calibrated source uncertainty at the resolution tested | the validation's `paired.py`, outside the tree (issue #36) |
| Scarce gDNA at 1/10 and 1/100 depth shows no blow-up, and dense probe panels hold | `quant_accuracy.py` on the test chromosome's depth and probe panels |
| Calibration's library gDNA is near truth on the in-scope strata; its per-locus count errs both ways, mostly its own error rather than the pooled share | `calibration_vs_oracle.py`, `prior_vs_oracle.py` |
| Under capture the capture reader carries most of the accuracy, and at STAR g00 most of the error | `quant_accuracy.py --bam` beside `--arm no_capture`, `ruler_vs_truth.py` |
| The transfer policy wins on unstranded × capture ON and the contaminated stranded rows, and loses to silent on unstranded × capture OFF | `policy_benchmark.py --panel ladder` |
| A refactor or speed-up is a numeric no-op | `rename_identity.py --check` against the frozen reference |
| Every instrument runs and regenerates | `preflight.py --full` |

## Next

1. **The capture-reader session, developed with the owner.** It starts from what the capture weight leaves,
   issue #1 and issue #2, with issue #3, issue #4, issue #5 and issue #6 read beside it. Its A/Bs need the paired
   proper score in the tree (issue #36). The trade-off between strata is the owner's to set.
2. **A review of the open issues**, each read for a problem a pre-release cannot carry.
3. **The pre-release**, with its supported uses and its stated failures: issue #4, issue #3, issue #1, issue #2,
   issue #10, issue #5, issue #12, as the capture reader leaves them; before any tag, issue #11.

After that, the rest of the `0.8.0` milestone (https://github.com/mkiyer/rigel/milestone/1), in issue order.
