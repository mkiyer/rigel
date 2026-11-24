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
  from 618,686 to 337,648, and the simulator's gDNA net from −325,530 to −77,872 fragments. The worst in-scope library
  net goes from −3.03 % to −1.76 % of the library's gDNA.
- The stranded strata hold, transcripts and genes within 0.7 %, and the annotated library net holds and fails on the
  same rows.
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
transcript error moves by at most 0.35 % and gene error by 0.9 %, and the annotated library net holds and fails on
the same rows. The
capture-ON references and the VCaP mix's are re-frozen on that tree.

The capture reader's clamp went on 2026-11-20 (`DESIGN.md` §7.2), on the owner's go with its cost stated: an object
reads its own level, not its enriched basin's peak. Pinned and fractional, against the tree before it:
- Under stranded capture the annotated library net is below 0.7.1's on every row at 5 % gDNA and above but the
  simulator's g05 (+27,644 → +18,647, against 18,082). The simulator's g50 goes
  from +53,207 to +30,002, its g98 from +25,398 to +18,575 and STAR's g50 from +34,964 to +22,297.
- Transcripts on that stratum rise 4.6 % on the simulator and 9.7 % on STAR, 24 % and 29 % at 5 % gDNA, where
  fragments move between the isoforms of genes with many (issue #83); genes fall 5.0 % on the simulator and
  hold on STAR. Bars 3 and 4 fail against the clamped reader there.
- Off capture the simulator's rows are unchanged to the digit. STAR's rows that read a false capture level rise by
  up to 1.5 % (issue #1); its zero-gDNA switch row reads 2.2 times the no-capture line on transcripts and genes.
- On the VCaP mix, DNA called annotated RNA falls by 4,990 fragments and RNA-library fragments called gDNA rise by
  8,899.

The capture-ON references and the VCaP mix's are re-frozen on that tree.

A locked object reads its strands, landed on 2026-11-21 (`DESIGN.md` §5, `EQUATIONS.md` §7.3; issue #1), after the
simulator's probes came to bind by their own GC and the capture panels were regenerated. Pinned and fractional,
against the tree before it:
- Every unstranded row, STAR at 98 % gDNA and every ladder and fl-gap row on the simulator's own alignments are
  bit-identical.
- On STAR's stranded rows transcript error falls on all six, most at no gDNA (431,143 → 287,469 off capture), and
  gene error rises on one, 5 % under capture, by 421 fragments, in duplicated gene families multimappers join into
  one EM locus.
- On the VCaP mix the DNA-origin share holds: 0.257245 against 0.257226, the truth 0.2518. Its identity reference is
  re-frozen on this tree.

The cuts' prices landed on 2026-11-22 (`EQUATIONS.md` §11, the best part; issue #76). An RNA placement that leaves a
cut's exons across a junction is priced at that junction's crossing level, read from the gDNA at its two splice sites
and capped by the junction cap, where gDNA's own placement there ran into the intron. Calibration does not move.
Pinned and fractional, per object, against the tree before it:
- Stranded capture ON, the annotated over-call per locus: simulator g05 25,142 → 20,972, g50 34,802 → 21,903, g98
  flat; STAR g05 11,104 → 8,849, g50 27,385 → 16,119, g98 29,672 → 25,796. The mass leaves annotated RNA mostly for
  nascent RNA, which holds no junction: simulator g50 annotated −18,577, nascent +15,019, gDNA +3,558.
- On the 24 in-scope rows at 5 % gDNA and above, transcripts fall on 8, rise on 3 and hold within 0.1 % on 13. The
  ruler's gDNA / isoform gap closes (g05 1.138 → 1.104, g50 1.049 → 1.007).
- It costs where the capture reading is false or thin and passes that level on: STAR unstranded capture OFF 5 %
  transcripts +9.8 % (that stratum +4.6 %, bar 4), STAR stranded capture OFF at no gDNA genes +4.7 %, 0.1 % gDNA under
  capture transcripts +7.7 % and genes +2.2 %, and a probe across the junction alone (the junction-probe design)
  transcripts +3.4 %.
- On the VCaP mix the DNA-origin share is 0.257370, and 6.3 % of the transcript mass moves between isoforms (genes
  0.30 %), which no truth there can score.

The index's one rule landed on 2026-11-23 (`DESIGN.md` §2.2): nascent entities whose start and end sites each lie
within `--end-tolerance` are one entity, transitively, and a single-exon transcript is its own nascent RNA. The suite
index's spans fall from 6,919 to 6,788. On the ladder's 16 conditions, on the same reads, pinned and fractional, the
in-scope transcript error moves −4,117 fragments and genes −50, every row within the EM's near-tie noise. The E-step's
row merge, a speed-up held to the tolerance (issue #84), landed beside it. The design phase on fused loci and the
EM's cap then closed for 0.8.0 (the owner, 2026-11-24): VBEM ships (`DESIGN.md` §9.6), MAP with removal by evidence is
issue #86, and annotated aliases are issue #87.

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

Every error is read per object (the owner, 2026-11-21; `SUCCESS.md`). Against 0.7.1 re-scored on the same reads
(2026-11-21), at 5 % gDNA and above, on every in-scope row of the ladder on both inputs and of both fl-gap panels,
HEAD's error is below 0.7.1's on transcripts (at most 0.72 of it) and genes (0.55), so bar 5 holds, and on each pool
per locus (annotated 0.64, gDNA 0.74) but nascent at STAR's stranded capture-OFF 5 %, 9,711 against 9,709. The
library nets that missed were 0.7.1's cancellation: at STAR's stranded capture-ON 50 %, its −1,380 is +16,265 over and
−17,645 under across loci. Below 5 %, stranded capture ON at 0.1 % gDNA is worse than 0.7.1, with about 120,000
nascent fragments read as annotated (genes 267,694 against 192,554), on the no-capture line too.

Calibration per object, Σ|P − O| less chance as a share of the library's gDNA: at 5 % gDNA 4.5 % unstranded capture
OFF, 4.5 % stranded capture OFF and 9.5 % stranded capture ON, which fails bar 1; at 50 % and 98 %, and on both fl-gap
panels, at most 3.0 %.

No preview ships (the owner's ruling, 2026-11-04): 0.8.0 waits for every `release-blocker` issue, in the order under
Next. The validation's receipts are in `$RIGEL_RUNS/prototypes/2026-11-02_preview_validation/`, the re-score's in
`$RIGEL_RUNS/prototypes/2026-11-21_regen/v071/`, and the per-object scores in `$RIGEL_RUNS/prototypes/2026-11-21_objects/`
and `$RIGEL_RUNS/prototypes/2026-11-21_bars/`.

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
`release-blocker` label is what it waits for. The last sessions have no pass bar (`DESIGN.md` §0b): each improves the
tool's current state, simply and robustly across many kinds of data, and what is left ships stated. In order, one
`src/` writer at a time:

1. **Performance** (issue #84): compute and memory, every change a proven numeric no-op.
2. **Code review** (issue #85): clean, concise, readable and maintainable, with no hard-coded constants.
3. **The release:** one ladder re-score, after rebuilding the suite and VCaP indexes with the one rule (their caches
   with them; the panels' reads were drawn on the old spans, so a truth-fed arm needs the span crosswalk or regenerated
   panels); a MANUAL limitation for issue #1, issue #2, issue #3, issue #4 with issue #5,
   issue #7, issue #9, issue #10, issue #12, issue #13, issue #23, issue #46, issue #78, issue #82 and issue #83, and
   for the cuts' prices' costs; the CHANGELOG; then issue #11.

CI runs alongside from now (issue #11). On its first run since September (fe0e3183) the native build compiled on Linux
and macOS. Two tests stopped the runs: one asserts a native kernel equals NumPy to the bit, which x86-64 misses by
1e-15, and one needs `minimap2`, which the CI environment lacks.

After that, the open issues outside the milestone, in issue order.
