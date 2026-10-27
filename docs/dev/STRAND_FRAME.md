# The strand frame — gDNA's length table in RNA's capture frame

Sandbox document for a session devoted to it (paused by the owner 2026-10-27; the open record is
`ISSUES: the-strand-frame-under-the-shared-law`). Nothing is in `src/`. The round-by-round notebook of 2026-10-25/26
(rounds one to three, SF3–SF12) is archived at `$RIGEL_RUNS/archive/docs_dev_2026-10-27/STRAND_FRAME_rounds1-3.md`.

## 1. The problem

The E-step scores an unspliced fragment's length by the ratio of the two origins' tables, `T_R(l) / T_G(l)`. The capture
weight of a placement is common to both origins and cancels from that ratio, so the ratio must be the origins' length
laws read in ONE frame (`DESIGN.md` §3.6b). The tool measures the two in different frames: RNA's law is its spliced
census, captured wherever RNA is; gDNA's is read off target, uncaptured. The shared length law keeps one law unless the
two differ by more than capture can make, and where they truly differ under capture it keeps the cross-frame
difference. So at a captured gene gDNA's table is its off-target law against an on-target RNA census, and the E-step
reads the frame as chemistry.

**The strand frame.** At an exon annotated on one strand only, an unspliced fragment contained in the exon is RNA or
gDNA at the same placements, under the same capture weight. Strand separates them: RNA reads sense at the library's
strand specificity `p`, gDNA reads either strand evenly. Fitting the split per exon × length reads gDNA's law in the
exon frame, the frame RNA's census is in.

**Why it matters now.** RNA's own-count length law (2026-10-26/27) removed a pseudo-count that had been paying back
gDNA's off-target frame. Without the frame, gDNA-long stranded capture-ON g05–g50 regress: the worst row is gDNA-long
g50 ss.99 capture-ON at +23.5 % genes, 0.31 points of true RNA. Only a frame correction recovers it. With the frame
beside the new law, that row reads −44.2 % genes.

## 2. The current construction (SF12a)

1. **The cells.** One-strand-exon × length counts of unspliced, exon-contained fragments, R1 on the exon's strand or
   against it, label-free (`cells_observed.py`), oriented by production's fitted R1 sense.
2. **The penalized fit.** `π = ½ + (p − ½) σ(a_e + d_b)`, with `p` production's strand specificity from spliced reads
   (`StrandModels.strand_specificity`). The contrast `d` is a centred first-difference walk inside the likelihood, its
   step at the Laplace evidence maximum. The exon odds `a_e` are free (profiled).
3. **gDNA's exon-frame law.** Each cell's gDNA share `1 − u` is spread over each exon's contained placements by a
   fixed point (`frame_fit.gon_from_fit`).
4. **The table.** `blend_tilt`: the inverse-variance blend of the exon-frame law and the off-target law times a
   random-walk tilt, which holds its edge where the cells are blind. It is applied where the shared law keeps two laws
   (`w > 0`), on the geometric path at `w`.

At ss 0.99 production's strandedness changes nothing that matters. gDNA's exon-frame law is within H² 7·10⁻⁶–6.5·10⁻⁴
of what the labels imply, means within 0.3 bp, and the fit's gDNA count at the one-strand exons is within 0.1 % of the
labels'.

## 3. The tradeoffs, measured (the four-arm A/B, 2026-10-27)

Frozen trees, pinned and fractional, against the tree before RNA's own-count law. Transcripts are read outside the
test chromosome's twelve terminus-cluster genes, whose isoform split the EM's start decides. Both changes:

| where | genes | transcripts | \|gDNA\| |
|---|---:|---:|---:|
| test gDNA-long, stranded capture-ON | −23.7 % | −3.0 % | −35.0 % |
| test RNA-long, stranded capture-ON | −10.2 % | −17.4 % | −6.5 % |
| suite RNA-short, stranded capture-ON | −17.4 % | −5.6 % | −92.0 % |
| suite RNA-long, stranded capture-ON | +1.9 % | −0.3 % | +90.9 % |
| the ladder (equal lengths), stranded | +0.1 to +0.4 % | −0.1 to +0.1 % | within ±0.2 % |
| test gDNA-long, stranded capture-OFF | −2.4 % | **+11.4 %** | +11.0 % |
| test RNA-long, stranded capture-OFF | −4.0 % | **+11.8 %** | −12.7 % |
| a 25 bp gap either way, stranded | within 0.4 % | **+0.8 to +3.4 %** | — |

**Where it fails, by mechanism.**
- **Weak strandedness with a scarce origin.** The free exon odds invent the minority origin at exons where it is
  absent: the positive part of a noisy fraction. The invented counts carry the majority's length law, and the blend
  trusts them, because their variance is read as Poisson counts' and the exon odds' uncertainty is left out.
  - gDNA-long g05 ss0.70 capture-OFF: the fit's gDNA law reads 116 bp against the labels' 259 (H² 0.70). gDNA's table
    goes 259 → 118 bp and transcripts cost +85 %.
  - RNA-long g05 ss0.70 is the mirror: 256 bp against 113 (H² 0.80), +64 %.
  - Under capture the same rows read H² 0.02–0.04, since capture enriches gDNA at the probed exons. Even so,
    gDNA-long g05 ss0.70 capture-ON costs +18 % transcripts while gaining 14 % on genes.
- **No gDNA (g00).** With no gDNA the efficiencies give `κ = 0`, the shared law keeps two laws, and the frame builds
  gDNA's table from RNA read as gDNA. Transcripts cost +9 % on the equal panel's g00 capture-ON rows and up to +25 %
  on RNA-long's.
- **Off capture.** The shared law's bound is 0 there, so `w = 1` on every resolved difference and the frame acts
  where it should be the identity.
- **A small real gap.** The frame's own error is 8–10 bp on the suite reference (exon frame against census frame), a
  third of a 25 bp gap.

**Why a gate is not the fix.** Off capture `κ` is noise around zero (3·10⁻⁴ at g00), not a clean switch. A gate also
leaves the weak-strand row under capture with its cost. The defect is the decomposition, and a gate only hides it.

## 4. Refused, with the number that killed each

- The fit's law used directly, unsmoothed: suite RNA-short +2,147 %.
- The census frame `normalize(r e^{−d})`, three forms: RNA-long g98 +81.6 / +172.1 %, suite RNA-short +146.3 %. Read only
  where both origins are measured: RNA-long g98 ss0.70 +1,349 %. Where an origin is absent, the contrast is identified
  only one-sidedly.
- Fixed-strength probability mixing of the two tables: RNA-long g98 +22.3 % at s = 0.5; gDNA's table alone +55.1 % at
  λ = 0.25; DNA-long g50 +42.5 % at s = 0.92. Where one origin is absent it floors the minority's share or caps it.
- SF7, the blend on the unsmoothed fit: harms g00 (ladder stranded capture-ON g00 +9.6 %).
- SF12b, SF12a's table times a census-frame factor: holds at the null, fails wherever one origin is scarce (test RNA-long
  ss.70 +468 %, suite RNA-short +677 %).
- The exon odds under a normal prior (2026-10-27, `integrated_fit.py`): worse. The invented RNA at g98 rises from 1.08–1.17×
  the labels' to 1.17–1.29×, because a normal prior on the log-odds cannot hold an exon with no minority and pulls
  exactly those toward the library's level. gDNA's law moves no nearer the labels'.

## 5. The principle the construction must meet

From the three review rounds: **a contrast the cells do not resolve is not used, and the two tables are then one law,
exactly; a contrast they resolve is used where it is measured; and where only one origin is measured, that origin's own
law carries the table.** Neither half shrinks toward equality by evidence weight. Where one origin is absent,
uncertainty in the log-ratio is certainty that the minority is absent.

## 6. The next step: an unbiased decomposition

gDNA's per-length count at the one-strand exons has an unbiased moment estimator:

    Ĝ_l = N_l − (S_l − N_l / 2) / (p − ½),     Var Ĝ_l = Σ_e N_{e,l} π_{e,l} (1 − π_{e,l}) / (p − ½)²  ≤  N_l / 4(p − ½)²,

pooled over exons, with `S_l` the sense count. Nothing is clipped at a cell, so nothing invents a minority, and a
minority the cells cannot see carries the variance that says so. gDNA's table is `g e^t`, `t` a first-difference walk
fitted to `Ĝ_l` at those variances. That is `g` itself where the cells say nothing (g00, weak strand at a scarce
origin, and off capture, where the exon frame is the off-target frame in expectation), and the exon frame where they
resolve it, with no gate. Open within it:

1. The placement de-tilt. A contained placement's count depends on each exon's length, so the per-length counts need
   each exon's gDNA weight. The per-exon moment estimates are unbiased but noisy, and clipping them is the bias again.
   Derive whether the de-tilt can read the pooled form (each exon's support weighted by its own unbiased count), or
   needs a presence prior: an atom at no minority, the shape a normal prior lacks.
2. Antisense RNA at one-strand exons reads as gDNA and pulls `τ` toward 1 (`ISSUES: the-strand-frame-under-the-shared-law`,
   exception 8). Can the cells carry a third, antisense-RNA term identified by its length law being RNA's?
3. Whether a library-level decision for exact equality is still needed once the variance is honest.
4. `p`'s own uncertainty, and a length-dependent strandedness (measured not to fool the penalized fit at the panel's
   strength, `scenarios_mech_lenp_equal`).

Then the four arms again. Read first the weak-strand g05 rows, the g00 rows and the small-gap panels, then every panel
per stratum. After that come the native banks: per one-strand exon × length bin × orientation counts of unspliced,
exon-contained fragments at the deposit, under 100 MB at int32 on GRCh38.

## 7. Instruments and receipts

- `$RIGEL_RUNS/prototypes/2026-10-25_shared_law/`: `penalized_fit.py` (the fit, its self-checks), `frame_fit.py`
  (`gon_from_fit`), `tilt_smooth.py` (`blend_tilt`), `sf12_arm.py`, `pen_inputs.py`, `sframe/` (the 15 rows' cells at
  the labels' `p`).
- `$RIGEL_RUNS/prototypes/2026-10-27_four_arm/`: the frozen trees (`site_head`, `site_wt`, `use_proto.py`), `four_arm.py`
  (one arm: a tree, optionally the frame), `cells_observed.py`, `frame_observed.py` (cells at production's `p`),
  `score_four.py` and `summ_four.py` (per stratum, transcripts outside the clusters), `cells/`, `frames/`, `runs/`.
- `$RIGEL_RUNS/prototypes/2026-10-27_strand_frame/`: `integrated_fit.py`, `compare_fits.py`, `gon_vs_oracle.py` (a fit's
  gDNA law against the labels'), and the small-gap configs (`test_reference_fl_smallgap_gdna25.yaml`, `_rna25`; the
  panels are simulated under `$RIGEL_RUNS/test_reference/`).
- Panels: the gap panels (test `scenarios_fl_gdna_long`, `scenarios_fl_rna_long`; suite `flgap_rna_long`,
  `flgap_rna_short`). The small-gap panels, simulated 2026-10-27, are test-chromosome prototypes and not part of the
  repo. The null and certification panels are `scenarios`, `scenarios_fl_equal200`, `scenarios_depth_d10`, `_d100`,
  `scenarios_mech_*_equal` and `seed43/scenarios`.
- Rules for the A/B: frozen trees, one mechanism per arm, pinned (`scan.total_threads=1`), fractional, against truth by
  read name; transcripts read beside genes, outside the clusters too; per stratum, never pooled.

## 8. Starting the session

> We are resuming the strand frame for Rigel 0.8.0. Read `docs/dev/STRAND_FRAME.md` first, then
> `ISSUES: the-strand-frame-under-the-shared-law` and `DESIGN.md` §3.6–§3.6b. The construction recovers stranded capture-ON
> gaps, but fails where one origin is scarce at weak strandedness, at g00 and at a small gap. Its decomposition
> invents the minority and reads the invented counts as precise. Derive the unbiased per-length moment decomposition
> of §6, settle the placement de-tilt (open question 1) on paper, and prototype it outside the tree on the cached cells
> (`$RIGEL_RUNS/prototypes/2026-10-27_four_arm/cells/`). Check it against `gon_vs_oracle.py` on the failing rows first.
> Then run the four arms on frozen trees and bring me the per-stratum table before anything enters `src/`.
