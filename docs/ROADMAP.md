# ROADMAP — the state and the next steps

The tool's state in one paragraph, the claims that still hold with the instrument that re-derives each, and the ordered
next steps. Every open problem's substance is in `ISSUES.md`; how performance is judged is `SUCCESS.md`; the changelog
is git.

## State

`src/`'s numbers are 0645adea's, the tree with the junction cap (2026-11-13, `DESIGN.md` §9.2), with the EM's gDNA
share read over one tally, landed on 2026-11-14 (`DESIGN.md` §9.3). That tree is the baseline every change is now read
against. On the screen (`SUCCESS.md`; 56 rows, pinned and fractional) it is identical to 0645adea on 1 row, within 1 %
on 40, better on 13 and worse on 2.
- The simulator's alignments hold no multimapper, so the two tallies were already nearly one: its 28 rows (the
  ladder's 16, the low-gDNA panel's 4, the fl-gap panels' 8) and the test chromosome's 12 move by at most 1.5 %, and
  no pool by more than 304 fragments.
- The STAR rung's contaminated in-scope rows fall on both axes. Transcript / gene error at g98, g50 and g05: stranded
  capture ON −41.1 % / −48.3 %, −8.4 % / −18.5 %, −2.9 % / −8.9 %; stranded capture OFF −53.1 % / −79.6 %,
  −6.3 % / −29.7 %, −1.2 % / −6.1 %; unstranded capture OFF +0.3 % / +0.4 %, −15.8 % / −40.1 %, −16.1 % / −46.5 %.
- The STAR rung's pools now carry the simulator's signs and nearly its sizes, so the errors the share's excess had
  cancelled show on both inputs. The annotated pool under stranded capture ON is an excess of +30.1 k, +35.0 k and
  +12.7 k at g98, g50 and g05 (before: −60.7 k, −7.7 k, −0.1 k; the simulator's alignments: +25.7 k, +53.1 k,
  +27.6 k), which is the capture reader's. gDNA called nascent RNA at g98 unstranded capture OFF is +302 k (before
  +70 k; the simulator's +292 k), which is `ISSUES: the-pooled-q-in-the-gdna-count`.
- Two rows rise. STAR g98 of the deferred stratum: +13.9 % / +15.4 %, the annotated excess +161.7 k → +191.0 k
  (`ISSUES: unstranded-capture-calls-gdna-rna`). STAR g00 unstranded capture OFF: +3.0 % / +7.5 % with better pools,
  all of it one duplicated locus whose RNA the clipped share had called gDNA and the EM now puts on the wrong copy
  (`ISSUES: gdna-in-duplicated-blocks-is-called-rna`).

The junction cap is an accuracy gain under capture, not the capture reader's repair: the hole between 0 and 1 % gDNA
stands (`ISSUES: the-reader-reads-false-gdna-as-capture`), and three of its rows trade one axis for the other
(`ISSUES: the-junction-cap-reads-bookings`). The identity references are re-frozen on this tree (`TESTING.md` §6).

The validation of ec996bba is closed. On the 16-condition ladder, on the simulator's
alignments and on STAR's, HEAD beats the tagged 0.7.1 in the three in-scope strata: Brier 0.23–0.88 of 0.7.1's on
identical reads, and transcript error 0.10–0.47 of 0.7.1's. It also states eight failures, each an open entry
but the one marked closed:
- unstranded sparse capture, where most true gDNA is called RNA in every build, and the deferred stratum (unstranded ×
  capture ON) (`ISSUES: unstranded-capture-calls-gdna-rna`);
- capture outside the annotation's partition: with exact probe copies, gDNA called RNA rises ×2.7–2.8 at stranded
  capture ON, mostly away from the copies (`ISSUES: capture-outside-the-annotation-partition`);
- on STAR at g00, the capture tail, whose few uncertain gDNA bookings swing ×2 across draw seeds, and the switch row,
  where under 0.1 % of templates move transcript error ×9.5, worse than 0.7.1
  (`ISSUES: the-reader-reads-false-gdna-as-capture`);
- a gene regression against 0.7.1 at the STAR g00 tail, 386 k against 328 k (the other one, STAR g98 κ0.99 OFF at
  37,304 against 22,898, closed with the share over one tally: 7,595);
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
   `ISSUES: gdna-in-duplicated-blocks-is-called-rna`, `ISSUES: vcap-exome-dna-called-rna`,
   `ISSUES: annotation-omission`, as the capture reader leaves them; before any tag, `ISSUES: the-release-builds`.

After that, `ISSUES.md`'s next tier, in its order.
