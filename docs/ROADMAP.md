# ROADMAP — the state and the next steps

The tool's state in one paragraph, the claims that still hold with the instrument that re-derives each, and the ordered
next steps. Every open problem's substance is in `ISSUES.md`; how performance is judged is `SUCCESS.md`; the changelog
is git.

## State

`src/`'s numbers are 81fc67b6's with the junction cap landed on 2026-11-13 (`DESIGN.md` §9.2), and that tree is the
baseline every change is now read against. (On the ladder; on the VCaP mix the tree before the cap already sat 30
fragments of 13.2 M from 81fc67b6, `ISSUES: instrument-ledger`, and the cap moves nothing there.) On the screen (`SUCCESS.md`; 56 rows, pinned and fractional) it is identical
to its predecessor on 40 rows, better on 7, within 1 % on 6 and worse on none, and three rows trade one axis for the
other (`ISSUES: the-junction-cap-reads-bookings`). Only stranded capture-ON rows with gDNA to anchor a junction move.
On the simulator's alignments, transcript / gene error is −2.5 % / −6.6 % at g98, +3.7 % / −20.0 % at g50 and +1.6 % /
−6.9 % at g05, and the annotated pool's excess falls 10 %, 28 % and 15 %; on the test chromosome's junction-probe
design it is −19 % / −36 %, −9 % / −46 % and −5 % / −41 %. The STAR rung moves by at most 1.2 %. The identity references
are re-frozen on this tree (`TESTING.md` §6). The cap is an accuracy gain under capture, not the capture reader's
repair: the hole between 0 and 1 % gDNA stands (`ISSUES: the-reader-reads-false-gdna-as-capture`).

The validation of ec996bba is closed. On the 16-condition ladder, on the simulator's
alignments and on STAR's, HEAD beats the tagged 0.7.1 in the three in-scope strata: Brier 0.23–0.88 of 0.7.1's on
identical reads, and transcript error 0.10–0.47 of 0.7.1's. It also states eight failures, each an open entry:
- unstranded sparse capture, where most true gDNA is called RNA in every build, and the deferred stratum (unstranded ×
  capture ON) (`ISSUES: unstranded-capture-calls-gdna-rna`);
- capture outside the annotation's partition: with exact probe copies, gDNA called RNA rises ×2.7–2.8 at stranded
  capture ON, mostly away from the copies (`ISSUES: capture-outside-the-annotation-partition`);
- on STAR at g00, the capture tail, whose few uncertain gDNA bookings swing ×2 across draw seeds, and the switch row,
  where under 0.1 % of templates move transcript error ×9.5, worse than 0.7.1
  (`ISSUES: the-reader-reads-false-gdna-as-capture`);
- gene regressions against 0.7.1: STAR g98 κ0.99 OFF, 37,304 against 22,898
  (`ISSUES: the-em-prior-share-excludes-drawn-multimappers`), and the g00 tail, 386 k against 328 k;
- VCaP exome DNA called RNA ×1.48 against 0.7.1 on identical reads (`ISSUES: vcap-exome-dna-called-rna`);
- annotation omission: unannotated unspliced RNA booked gDNA, and retained genes damaged ×1.4–3.0
  (`ISSUES: annotation-omission`).

No preview ships (the owner's ruling, 2026-11-04): a pre-release waits for the capture reader and a review of the open
issues. The validation's receipts are in `$RIGEL_RUNS/prototypes/2026-11-02_preview_validation/`.

## Claims that hold

| claim | re-derived by |
|---|---|
| HEAD beats 0.7.1 on the ladder's three in-scope strata, on both inputs | `quant_accuracy.py` per stratum, pinned and fractional, `--bam` for the STAR rung, beside the tag |
| HEAD's source probabilities beat 0.7.1's by the proper score on identical reads, and most of the exchange is compatible with approximately calibrated source uncertainty at the resolution tested | the validation's `paired.py`, outside the tree (`ISSUES: instrument-ledger`) |
| Scarce gDNA at 1/10 and 1/100 depth shows no blow-up, and dense probe panels hold | `quant_accuracy.py` on the test chromosome's depth and probe panels |
| Calibration's library gDNA is near truth on the in-scope strata; its per-locus count errs both ways, mostly its own error rather than the pooled share | `calibration_vs_oracle.py`, `prior_vs_oracle.py` |
| Under capture the capture reader carries most of the accuracy, and at STAR g00 most of the error | `quant_accuracy.py --bam` beside `--arm no_capture`, `ruler_vs_truth.py` |
| The transfer policy wins on unstranded × capture ON and the contaminated stranded rows, and loses to silent on unstranded × capture OFF | `policy_benchmark.py --panel ladder` |
| A refactor or speed-up is a numeric no-op | `rename_identity.py --check` against the frozen reference |
| Every instrument runs and regenerates | `preflight.py --full` |

## Next

1. **The capture-reader session, developed with the owner.** It starts from
   `ISSUES: the-reader-reads-false-gdna-as-capture`, with `ISSUES: capture-outside-the-annotation-partition`,
   `ISSUES: unstranded-capture-calls-gdna-rna`, `ISSUES: vcap-exome-dna-called-rna` and
   `ISSUES: an-anchor-only-landscape-reads-as-capture` read beside it. Its A/Bs need the paired proper score in the
   tree (`ISSUES: instrument-ledger`). The trade-off between strata is the owner's to set.
2. **A review of the open issues**, each entry read for a problem a pre-release cannot carry.
3. **The pre-release**, with its supported uses and its stated failures: `ISSUES: unstranded-capture-calls-gdna-rna`,
   `ISSUES: capture-outside-the-annotation-partition`, `ISSUES: the-reader-reads-false-gdna-as-capture`,
   `ISSUES: the-em-prior-share-excludes-drawn-multimappers`, `ISSUES: vcap-exome-dna-called-rna`,
   `ISSUES: annotation-omission`, as the capture reader leaves them; before any tag, `ISSUES: the-release-builds`.

After that, `ISSUES.md`'s next tier, in its order.
