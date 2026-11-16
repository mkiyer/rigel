# ISSUES — the open problems

One entry per open problem in the current source. An entry is a `### kebab-name` heading, then a `priority:` line
(now / next / later / parked) naming its `stage:` and `kind:`, then what is wrong, the evidence, where it sits in the
code and what is known about a fix. The stage is where in the tool's run the problem first acts, in run order: `index`,
`scan`, `accumulator`, `calibration`, `pre-EM`, `EM`, `post-EM`; work outside the run is `panels`, `instruments` or
`release`, and a ledger is `several`. Within a tier the entries follow that order. Cite an entry as `ISSUES: name`. The
plan is `ROADMAP.md`, rulings are `DESIGN.md`, lessons are `TRAPS.md`. Numbers are pinned, fractional and scored by
read name unless an entry says otherwise; "the STAR rung" is the 16-condition ladder aligned by STAR at the recommended
settings.

---

## OPEN — now: the capture reader

### the-reader-reads-false-gdna-as-capture
`priority: now — the largest error on aligned input · stage: calibration · kind: defect`
`capture_efficiency._read` turns each object's deconvolved gDNA mass, ψ's posterior-median booking, into a capture
efficiency by reading it as a Poisson count on its support against the landscape's bulk, one grid cell. Calibration's
false gDNA on a few objects then reads as enormous capture, which every EM length carries. Where the bulk sits at the
resolution wall (STAR g00, 5.0×10⁻⁷ per bp) a booking ψ prefers over none by e^4.6 reads as e^250; an object holding under
one expected gDNA fragment cannot tell depletion from capture. Since 2026-11-15 the reading is raised to the capture
weight (`DESIGN.md` §7.2), which fades it where the library's reads sit on objects holding no gDNA; the false readings
remain, at a lower weight. Transcript / gene error, the shipped reader of 1f97ff96 → the weight landed, beside the
no-capture line (one seed):
- STAR g00 κ0.99 OFF, the switch row: 1,897,401 / 414,763 → 780,582 / 120,075 at weight 0.187, against the line's
  200,330 / 26,837 and 0.7.1's 1.38 M transcripts and 213 k genes. Calibration books gDNA of its own on objects holding
  about a fifth of the read mass in a library that has none, so the weight counts it as visible, and the reader keeps
  9.3× of its 149,367× saturation. Before the weight, the 229 regions and 581 boundaries read above 10× carried the
  gap: setting only them to 1 closed 99 %.
- STAR g00 κ0.99 ON, the capture tail: 1,188,241 / 385,519 → 672,234 / 250,244 at weight 0.052, against the line's
  652,274 / 239,710 and 0.7.1's 328 k genes. Its few uncertain gDNA bookings swung RNA called gDNA ×2 across draw seeds
  under the shipped reader; the other seeds are not re-run.
- STAR g00 unstranded OFF: 342,493 / 146,647 → 255,694 / 55,325 at weight 0.259, against the line's 204,248 / 29,026.
- Under capture the reader is most of the accuracy (every efficiency at 1: stranded ON transcripts ×6.7, genes ×14), and it
  amplifies every upstream change, so a scan or scoring change is read beside `quant_accuracy.py --arm no_capture`.
- **The hole between 0 and 1 % gDNA under capture**, the low-gDNA side panel (simulator's alignments, stranded ON): at
  0.1 % 4,662,866 / 438,698 → 1,049,271 / 228,305 at weight 0.129, against the line's 625,202 / 232,081, so genes are
  on the line and transcripts 1.7× it (the release ruling reports transcripts there and does not require them). At 1 %
  the weight is 0.97 and the reading stands (513,558 / 70,556 → 509,915 / 69,957; the line 648,791 / 259,961).
- **False depletion off capture**, the same panel at capture OFF: 0.1 % 373,318 / 26,528 → 183,280 / 12,238 at weight
  0.280, and 1 % 137,013 / 10,926 → 131,453 / 11,113 at 0.597, against the line's 140,017 / 9,567 and 126,849 /
  11,188. The reader reads below 1 on purpose: its bulk is the largest basin, which under a captured majority is the
  captured level (`tests/calibration/test_capture_efficiency.py`), and no panel holds that case. The rest of this error
  is the floor's, which still needs a ruling.

The 2026-11-11 audit also finds enriched basins on STAR capture-OFF at g50, g05 and g00 under both strandednesses.
Removing the reader's cap therefore changes all six rows despite absent capture; it is an exact no-op only at g98
on STAR and at every tested level on the simulator's alignments. Of those six STAR rows only g00 stranded was run
end to end: transcript error rises 0.2 % (1,896,916 → 1,899,995). The other five were not run because the capture-ON
A/B already rejected the candidate (`$RIGEL_RUNS/prototypes/2026-11-10_length_gap/`, section 7).

Known: candidates centred on ψ's booking keep the tail (widening the count by ψ's Var(log f_g) is refused); candidates
reading only the object's own evidence lose the deferred stratum. Library weights that lose a stratum or a panel: from
ψ's variance, or from the count the strand solve vouches for (unstranded capture ON fades), and from each object's count
against the landscape's spread (one probe fades: `$RIGEL_RUNS/prototypes/2026-11-15_capture_fade/`). The landed weight
cannot tell false gDNA from a small panel: both are a small set of objects where gDNA shows. Open: the bulk's
uncertainty propagated jointly with the object's; a regularised estimator that shrinks an object toward uniform capture
from its own evidence, with no count threshold or clamp. Not candidates, each screened end to end on 2026-11-10 to
11-13: a reader fitted on the located objects alone (best at g98 ON, and at zero gDNA off capture it reads capture, +446
% / +736 %), and a reader in which a gene's objects share a prior probability of enrichment (its zero-gDNA limit is a
library switch that does not hold on other probe designs). Receipts: `$RIGEL_RUNS/prototypes/2026-11-10_length_gap/`,
sections 10 and 15–18. Instruments: `quant_accuracy.py --arm no_capture`, `--bam`, `ruler_vs_truth.py`; the low-gDNA
side panel (`scripts/sim/configs/gdna_ladder_lowg.yaml`).

### off-target-rna-dilutes-the-capture-weight
`priority: now — DNA called RNA at high gDNA on a captured library · stage: calibration · kind: defect`
The capture weight is a mean over the library's one-strand reads (`DESIGN.md` §7.2), so RNA on objects whose own gDNA is
too thin to show lowers it for the whole library, a readable panel included. In a deep library the off-target exons hold
gDNA of their own and the weight stays near 1; where the off-target gDNA is thin it does not.
- cfRNA-B (real, stranded, captured; 1.39 M unspliced fragments; 93 % gDNA in the EM; 17,306 regions read above 2×,
  its panel read consistently at the saturated level): weight 0.946. Off-target RNA on highly expressed exons (ACTB,
  HBB, TMSB4X, EEF1A1: 5.4 % of the unspliced mass, its booked gDNA about 0.01–0.1 per object) fades the saturation
  from 4,239× to 2,702×. No truth; read by the strand of each gene's reads (reverse-stranded: a gene's RNA sits on one
  strand, gDNA on both). Of the 1,000 largest gene moves against the shipped reader, 470 of the net +487 fragments go
  to 816 genes whose reads are strand-balanced and already over-called (9,328 → 9,798, against 7,350 by strand), and 12
  genes go from 0 to 5 or more (PAX1 0 → 18).
- On the simulator the cost is small: the one-probe test chromosome pooled with the whole-genome capture-OFF library at
  5 % gives 0.884, and the probed gene at a forced 0.87 reads +70 fragments at 5 % gDNA and nothing at 25–98 %.

It is DNA called annotated RNA at high gDNA, the error the release ranks first. Not run: a weight from whether CAPTURED
gDNA would show on each object (the saturated level × the bulk × its support), which would keep cfRNA-B but trusts the
saturated level that the hole and STAR g00 inflate. Receipts: `$RIGEL_RUNS/prototypes/2026-11-15_capture_fade/`
(README, the external review, points 1 and 7).

### capture-outside-the-annotation-partition
`priority: now — a stated failure of the preview · stage: calibration · kind: defect`
Capture is read per object of the annotation's partition, and a probe also binds sequence the annotation does not
separate: an exact copy of its target in an intron, a paralog or intergenic space. In the counterfactual that gives every
exact probe copy the probes' binding opportunity, stranded capture-ON gDNA called RNA rises ×2.7–2.8 (simulator 797 k →
2.17 M, STAR 637 k → 1.76 M) and MCB ×26–76, while 0.7.1 stays flat on the simulator's alignments. The error lands mostly
away from the copies (g50 κ0.99: 473 k → 1.34 M beyond 500 bp of any changed site, 14 k → 84 k at them), so the copies
shift calibration's reading genome-wide, by a mechanism not established. At unannotated copy sites the region factor
reaches only 4.6 where 253 k gDNA fragments are newly captured. A copy whose fragments all multimap also reads the
uncaptured bulk (`ISSUES: the-draw-reads-unique-traffic`). Instrument: the copies panel and `copies.py` in
`$RIGEL_RUNS/prototypes/2026-11-02_preview_validation/`.

### unstranded-capture-calls-gdna-rna
`priority: now — two stated failures of the preview, one of them the deferred stratum · stage: calibration · kind: defect`
On an unstranded captured library the strand channel is dead and capture depletes the intergenic space that measures the
gDNA density, so a gene's unspliced coverage fits RNA, mostly the synthetic nascent spans, as well as gDNA.
- Sparse probes at κ 0.50 on the test chromosome: most true gDNA is called RNA in every build, confidently (g98: 919,552 of
  1,146,600 under HEAD, 1,026,110 under 0.7.1; MCB 805,566 of a 827,074 Brier). The capture reading reads about 1 there,
  so the probe-enriched gDNA goes to RNA. Dense probes at κ 0.50, and κ 0.70 or 0.99, do not fail.
- The ladder's deferred stratum (unstranded × capture ON) on STAR: gDNA called RNA 3.30 M against 0.7.1's 2.97 M, a one-way
  under-forecast (MCB 458 k of a 1.55 M Brier). Until 2026-11-14 the EM's share was read 4–7 % too high on the STAR rung
  and hid part of this. With the share over one tally the rung reads as the simulator's alignments do: the gDNA pool's
  deficit is 2.08 M, 0.95 M and 96 k at g98, g50 and g05 (before: 1.67 M, 0.72 M, 51 k; the simulator's alignments:
  1.96 M, 0.92 M, 86 k), and g98's transcript / gene error rose 13.9 % / 15.4 %.

The one discriminating signal identified is the spliced fraction: mature RNA at a multi-exon transcript makes
junction-crossing fragments at a rate its geometry and length law set, and gDNA read as RNA carries almost none. No
likelihood term reads it. This stratum is read, never pooled, beside every capture-reader change.

### vcap-exome-dna-called-rna
`priority: now — a stated failure of the preview, on real data · stage: calibration · kind: defect`
On the VCaP mix (an exome-captured DNA library and a stranded transcriptome library, truth by library of origin) HEAD
calls 447,371 of the exome DNA's fragments RNA against 0.7.1's 301,566 on identical read IDs (×1.48; Brier on the DNA
×1.21), about half of it booked to synthetic nascent spans. HEAD's overall Brier edge (×0.95) rests on the
transcriptome library, whose origin is only an upper bound on RNA. The DNA library is unstranded and captured, so this
is `ISSUES: unstranded-capture-calls-gdna-rna`'s mechanism inside a mixed library; an older measurement also put
calibration's deficit on the mix's no-RNA objects (P/O 0.922 at full depth). Not dissected by stage. Instruments:
`paired.py`, `score_real.py` in `$RIGEL_RUNS/prototypes/2026-11-02_preview_validation/`.

### an-anchor-only-landscape-reads-as-capture
`priority: now — a library with no gDNA reports capture · stage: calibration · kind: defect`
With no gDNA evidence the landscape trains on its zero-count anchors alone. That density has one maximum, but `_render`'s
column-normalised Gaussian (`landscape.py`) makes ripples at both grid edges and `census` reads the top one as an enriched
basin: an empty or one-fragment BAM reports `gdna_capture_saturation` 9.65 beside `gdna_density_global` 0 and lengthens
every `em_effective_length` by about 6 %, where an uncaptured library should read 1. Repair the render's edges or the
census and keep the anchor-only fit. Test first: anchors alone give one basin, saturation 1 and zero-count efficiencies
of exactly 1.

## OPEN — next: the preview's other stated failures, and the release

### splicing-artifacts
`priority: next — the mechanism is built; its validation on real artifact-bearing data waits for the cluster · stage: scan · kind: defect`
A gDNA read the aligner splices reads as certified RNA. At the recommended STAR settings the aligner also reports the
contiguous reading and a multimapper's hits are priced by their sequence at the library's own mismatch rate, so on the
STAR rung gDNA called RNA through a spliced hit is under 0.11 % of all gDNA called RNA. Unvalidated: what the rung cannot
hold (STAR itself, retrocopies, real errors, qualities and clips), to be read on the VCaP mix at the recommended settings.
The index's blacklist (`splice_blacklist.py`, applied by the resolver) errs both ways:
- on the VCaP RNA library it rejected 598,010 records where a pure-gDNA rate predicts about 20,000, and a rejected genuine
  read becomes an unspliced, gDNA-eligible fragment whose footprint spans the intron;
- the `alignable` catalogue is 60.5 % half-blind (one side's anchor 0), its anchors stop at half the longest read, it
  lists 20.7 % of annotated junctions, and it was simulated at STAR's production settings, not the recommended ones;
- `aggregate_splice_blacklist` applies `min_count` (2, underived) per read-length row before summing, so a junction seen
  once at each of several read lengths never enters;
- untested rules: the anchor's `<=`, the OR of the two anchors, cumulative anchors across introns, a junction rejected on
  one mate and kept on the other, and one kept junction shielding a rejected one.

The spliced fragments train κ and RNA's length law, so artifacts bias both; κ already fits only the genuine junctions.

### the-junction-cap-reads-bookings
`priority: next — the landed cap's stated limits · stage: pre-EM · kind: limitation`
The junction cap (`capture_eff_length.fit_junction_cap`, landed 2026-11-13) is one library level, fitted from each
junction's spliced mass against the unspliced RNA calibration booked in its own exons (`EQUATIONS.md` §11). What it
assumes without checking:
- **The anchors are bookings, nascent RNA included.** With the true annotated RNA at the anchors the fitted cap is
  within 2 % of the best cap against truth; with the bookings it is 6 % low on the ladder at g98 and 5 % high on the
  exon-only test panel. Taking nascent RNA off by the junction's intron is unstable at g98.
- **It reads opportunity through gDNA's length law**, which is off on aligned input
  (`ISSUES: the-gdna-length-contrast-is-ill-conditioned-on-aligned-input`). At STAR g50 stranded ON the counts then ask
  for a lower cap, and a fit allowed to lower it falls to the best object beside every junction (transcripts +12.6 %
  against the saturated level's). So the fit only raises the cap, and stays one-sided until that law is right.
- **One level per library.** The cap bounds a junction's price and does not set it: on the sparse probe design
  junction levels span 62–791× the bulk, and the cap moves 358× → 830× on six junctions with no visible effect. With
  weak evidence the likelihood is flat and the first maximum is kept.
- **The trades.** Simulator stranded ON, transcripts against genes: +3.7 % for −20.0 % at g50 (+7,213 fragments for
  −17,132) and +1.6 % for −6.9 % at g05; on the exon-only design at g98, genes +2.7 % for transcripts −3.0 %.
- An exon between two junctions anchors both, so the joint likelihood counts its RNA twice; a junction's multiplier
  is the unweighted mean over its transcripts; the maximum is found on a 600-level grid.
- It does nothing at 0–1 % gDNA or on unstranded libraries (the counts raise nothing there), and it leaves the cuts
  inside exons, half the length's within-gene scatter on the ladder, with no candidate.
Receipts: `$RIGEL_RUNS/prototypes/2026-11-10_length_gap/`, sections 12, 13 and 18, and
`$RIGEL_RUNS/prototypes/2026-11-13_landing/`. Instruments: `quant_accuracy.py` on the screen, `ruler_vs_truth.py --scale`.

### the-capture-likelihood-leans-toward-rna
`priority: next — under capture the EM calls gDNA RNA even when calibration is right · stage: EM · kind: defect`
On the ladder's g98 stranded capture-ON simulator row, the true-share prior with shipped lengths leaves +7.8 % net
RNA; with each component's exact capture yield it leaves −0.1 %. Without a count prior, exact yields reduce net RNA
from +70.0 % to +8.0 %. The prior masks much of the length error. This closes aggregate source imbalance in that
cell, not every transcript error or the gross source-confusion floor. The length discrepancy includes the object's
capture readout, junction pricing, and applying gDNA efficiencies to RNA opportunities. The simulator-yield
intervention is not a validated STAR oracle: it worsens STAR, for reasons not yet established.

Removing the reader's enriched-basin cap must not land. It improves g98 simulator transcript/gene error by
20.4 %/20.9 %, but worsens STAR by 3.9 %/2.0 %, transcript error at g50 by 9–11 %, and transcript error at g05 by
about 39 % on both inputs. Uncapping improves class-mean lengths while increasing within-gene length-error spread;
class-wide corrections can instead transfer annotated RNA into synthetic spans. The cap's flattening helps the tested
isoform splits, but this does not establish clipping to a basin peak as the generally correct regularizer. Receipts:
`$RIGEL_RUNS/prototypes/2026-11-10_length_gap/`, sections 1 and 5–7; the exact-yield comparison is in
`$RIGEL_RUNS/prototypes/2026-11-10_error_map/`.

The older unequal-length finding remains distinct: on the test chromosome's gDNA-long g50 stranded capture-ON row,
the oracle gDNA count and gDNA table in the capture frame still leave 5.7 k excess RNA fragments, 0.5–1.7 % per captured
gene. A candidate mechanism is pricing each object at one length-independent efficiency. The equal-length ladder
does not settle that case. Instruments: `ruler_vs_truth.py`, `quant_accuracy.py --arm oracle`.

### gdna-in-duplicated-blocks-is-called-rna
`priority: next — the STAR rung's next error now that the share's two sides are one tally · stage: EM · kind: defect`
On an aligner's input, gDNA reported at several places is called RNA where the reference repeats itself. STAR g98
stranded capture ON, per read with the share over one tally (2026-11-12, an arm on HEAD 5341ea06): the 610,692 gDNA
read pairs at several places are called nascent 20,974 and annotated 11,261 fragments (5.3 %), against 2,932 and 5,746
on the simulator's alignments (1.4 %), where each sits at its true place. 85 % of the nascent fragments are in five
1-Mb bins of chr21 (5–8 Mb, 34 Mb, 44 Mb), and three quarters are read pairs with exactly two places. The same blocks
carry two thirds of the STAR rung's extra transcript error at g50 and g05: about 300 transcripts holding 3 % of the
true fragments, every read a multimapper, the split between copies wrong (one gene −16.6 k of 36.3 k at g05, a
neighbouring copy +8.8 k on 14.6 k). Until 2026-11-14 a share clipped at one hid the gDNA half of this. Located, not
dissected; part of the transcript half compares against the simulator's alignment, which knows the copy a read came
from. Receipts: `$RIGEL_RUNS/prototypes/2026-11-12_star_gap/` (sections 4.1, 4.2).

The clipped share hid part of the split too, and the one in-scope row that rose when the share was repaired is this
entry (STAR g00 unstranded capture OFF, transcript / gene error +3.0 % / +7.5 %). One locus of two chr21 copies (6.08 Mb
and LINC00319 at 43.4 Mb; 59 k units, 98 % multimappers) was handed a share of one, and 18.2 k of its RNA fragments
were called gDNA in a library with none. With the share over one tally 3.1 k are. Of the 15.1 k returned, 12.6 k go to
the copy already over-called (true 25.1 k; 30.4 k before, 43.0 k after) while LINC00319 stays 21 k short of 35.9 k: the
pools improve by 15 k and the gene error rises by 10.4 k, all of it in that locus. Receipts:
`$RIGEL_RUNS/prototypes/2026-11-14_share_frame_landing/` (`why_row.py`).

### the-release-builds
`priority: next — the owner's actions before any tag · stage: release · kind: release`
The build is portable by default, the smoke tests load every native extension, and every wheel and sdist carries
`THIRD_PARTY_NOTICES.txt`. Left, the owner's:
- dispatch `ci.yml`, and `publish.yml` with `publish: no`, on the release commit: the native code, restructured into
  `scan/`, `resolve/`, `calibration/`, `testing/` and the `em_*.h` headers, has been built only with clang on macOS
  arm64, and CI has not passed since it became dispatch-only;
- the bioconda recipe: the auto-bump bot copies only the version and sha256, so the first PR carries `conda/build.sh` and
  `meta.yaml`'s `additional-platforms` (`linux-aarch64`, `osx-arm64`) by hand;
- before any cluster quant, check that the production index's `manifest.json` records `sources.alignable_zarr`, or rebuild
  with `--alignable-zarr`: `TranscriptIndex.load` applies `splice_blacklist.feather` only then, and otherwise warns and
  runs with splice-artifact detection off (`summary.json`'s `sj_blacklist_loaded: false`).

### annotation-omission
`priority: next — a stated failure of the preview, which assumes a complete annotation · stage: several · kind: limit`
Rigel admits RNA only where the annotation does, so a missing transcript's RNA leaves the tables and damages what stays.
On the STAR rung with annotation withheld from the quantifier:
- unannotated unspliced RNA is booked gDNA by the population rule: 56 % of a withheld lncRNA's RNA (0.7.1: 81 %), most of
  the rest other genes' nascent spans;
- unannotated splices are withheld by the scan's annotated-path rule and go unscored: 26 %;
- retained genes are damaged ×1.4–3.0 and retained transcripts ×1.4–4.9 with lncRNA withheld, and withholding processed
  pseudogenes (1 % of the RNA) raises the g00 rows' transcript error ×1.7–1.9;
- under isoform omission about a third of the missing isoforms' RNA leaves the gene table;
- on real captured libraries the capture weight releases such RNA to a transcript on the other strand: on cfRNA-A and
  cfRNA-C a single-exon lncRNA gains 687 and 530 fragments from reads 99.8 % antisense to it. Calibration books that RNA
  gDNA, and the shipped reader read its region as each library's most extreme capture (2.6 M×, 37 k×), which kept it in
  gDNA; at the weight the region reads far less. The weight also takes 50–302 fragments off two lncRNAs of the same kind.

The same rule books real unannotated transcription as gDNA in calibration (all of the test chromosome's `test_blank`
shadow RNA), biasing any pooled structural rate. Two routes are open: admitting RNA at unannotated slots needs a ruling,
since that RNA has no transcript to go to; and a template whose mates fit one gene along a path no transcript holds could
be counted to the gene alone, though that would re-admit splice artifacts, which concentrate on unannotated junctions.
Instruments: `omission.py`, `omission_table.py` in `$RIGEL_RUNS/prototypes/2026-11-02_preview_validation/`.

## OPEN — later

### overlapping-synthetic-shadows
`priority: later · stage: index · kind: design`
`index.create_nrna_transcripts` clusters TSS and TES within a tolerance per strand and makes one synthetic nascent span
per merged cluster, so a gene with several isoform clusters carries near-coincident twins: 6,366 of the ladder index's
6,919 spans overlap a same-strand twin, 5,879 at ≥ 90 %. The EM cannot tell twins apart and puts a family's mass on an
arbitrary one, so a per-span nascent count is not identifiable (the MANUAL says so) and `quant.feather`'s
`count_nascent` can read about 0 while a sibling twin holds the family's mass. The pool rows and the transcript table do
not see it. The choice: merge in the index (one span per gene per strand, which changes the EM's components and the
panels' nascent truth ids) or in the scoring.

### clips-are-walked-along-the-genome-only
`priority: later — rare on STAR input at the recommended settings · stage: scan · kind: defect`
The pair rule walks each mate's 5′ soft clip back along the genome (`FragmentResolver::validate_pair`). A clip that hides
a junction walks into the intron instead, so an inward molecule reads outward and its read pair is withheld: 0–11 per
10 M templates on the STAR rung, 131–522 per cfRNA library, 4,078 on the VCaP mix. Walking the clip along the witness
transcript's path recovers them. No length downstream reads a clip: the deposit and the scorer measure a candidate on its
aligned blocks. It matters for an aligner that clips short anchors (minimap2).

### input-mistakes-pass-silently
`priority: later · stage: scan · kind: defect`
- A BAM without `@HD SO` that is in fact coordinate-sorted still passes. Refusing it needs the count of read pairs whose
  mate the flag says is mapped and the record run lacks (`summary.json`'s `read_pairs` holds a looser count), and a
  threshold on it is the owner's call.
- An unpaired record stops the scan with the scanner's own RuntimeError: a clear message, but a traceback.

### the-tail-rule-waits-for-truth
`priority: later — needs a simulation with adapter read-through · stage: scan · kind: question`
Every chimeric alignment is withheld, including one whose foreign piece covers only a read's 3′ part past its mate's 5′
end. Dropping that piece and keeping the pair would recover 4.9 % of cfRNA-A's templates, which suggests adapter
read-through; on the STAR rung it accepts 700 gDNA templates of g50 κ0.50 OFF, every one misplaced, because the rung's
chimeras are STAR's false mate merges. The simulator writes no read-through and trims no 5′ end, so no data with truth
holds the mechanism; it is decided by `quant_accuracy.py` on a simulation that does.

### drain-contaminates-certified-rna
`priority: later · stage: accumulator · kind: defect`
The second pass's drain draws a spliced hypothesis for some true-gDNA fragments, so the certified-RNA banks hold gDNA
deposits (233 records at g50 ss.99 OFF and 1,482 at g98 ss.99 ON, 1.9 % of that channel, on an earlier draw). The leak is
exact posterior sampling, and a no-leak counterfactual moved calibration's metric by at most −0.60 %, so the harm is the
certainty claim and no in-solve correction is licensed. The oracle records the leak per bank rather than refusing the
frame (`scripts/design/_oracle.py`). To build: `DrainQC` (`scan_payload.py`) records the expected leak, Σ q_null.

### the-draw-reads-unique-traffic
`priority: later · stage: accumulator · kind: gap`
The second pass's placement term (`second_pass.placement_abundance`) reads pass one's traffic, which counts only the
fragments pass one deposited, so a copy's term is its rate times its own mappability; gDNA's floor bounds it below at
the pass-one bulk. A copy that capture enriched but whose fragments all multimap reads the uncaptured bulk. On the STAR
rung's held multimappers the floor removes most of the per-region gDNA error off capture (×0.07 of the draw without it)
but little under capture (×0.77). The exact rate is the count over the copy's own unique opportunity, which needs the
reference's mappability, or rates that count the held fragments, iterated to a fixed point. Instrument:
`$RIGEL_RUNS/prototypes/2026-10-29_draw_census/truth_regions.py`.

### the-gdna-length-contrast-is-ill-conditioned-on-aligned-input
`priority: later — what is left needs the STAR rung's own alignable store · stage: calibration · kind: defect`
`fl.deconvolved_gdna_counts` reads gDNA's length law as a two-pool contrast, g = [(1 − a1)·f0 − (1 − a0)·f1] / (a0 − a1),
which amplifies any difference between the pools by 1/(a0 − a1). With the region lengths over the observable genome the
law is near the simulator's, but on the STAR rung the intergenic purity still reads 0.79 (off capture) and 0.62 (on)
where the oracle clips at 1: the pericentromeric DNA the aligner places only half of (chr21:10.6–13.0 Mb and
chr22:12.6–15.3 Mb hold 69 k of the rung's 71.6 k unplaced gDNA pairs) reads emptier than it is, because the GRCh38 store
masks only 29 kb beyond N there. The resolution weight also reads the purities' scalar error, never the per-length noise
the division amplifies. Instrument: the STAR rung with the oracle's gDNA law injected.

What swings the contrast on the STAR rung is the second pass's deposits (2026-11-12, `pipeline.library_fl_models` on
the cached STAR scans, on pass one's tally and on the drained one). The drain adds the held multimappers to the DNA
length pools, two thirds again of the intergenic pool and a quarter to a third of the intron pool, 9–15 bp longer than
the unique fragments there, and the fitted law's mean moves from 217.3 to 248.4 bp at g50 stranded capture ON and from
222.0 to 210.4 at g98 (the simulator's law is 217.1). On the simulator's alignments the drain puts nothing into those
pools and the law does not move. With the true law injected, the annotated pool of STAR g50 stranded capture ON moved
+33 k fragments. Untested: the law read on pass one's tally in the drained one's place. Receipts:
`$RIGEL_RUNS/prototypes/2026-11-12_star_gap/` (section 4.3).

### eb-shrinkage-magic-ess
`priority: later · stage: calibration · kind: defect`
gDNA's length law is shrunk toward the all-fragment pmf, mostly RNA's where gDNA is a minority, at a fixed ESS
(`CONSTANTS.fragment_length.pool_prior_ess` = 1000, `fl._smooth_eb`); RNA's law already reads its own counts
(`fl.own_count_law`). On the test chromosome's fl-gap panels at g05 under capture, gDNA-long's law reads 246.6–248.2 bp
against a uniform-frame truth of 261 and RNA-long's 125.9–127.9 against 114. The contrast is itself noisy there (H²
0.05–0.22), so gDNA's own-count law must read the two contained pools as counts at their purity weights. The
equal-length ladder cannot see it.

### the-gdna-prior-enters-psi-twice
`priority: later — its fix fails alone; it needs a per-object partner · stage: calibration · kind: defect`
ψ's grid (`psi_kernel.h`) holds the Jeffreys reference ½·log f + ½·log(1 − f) and adds the fitted landscape's arm on top,
so wherever a landscape is fitted gDNA's ½·log f enters twice: up to 0.85 nats on log f_g, and an unstranded
no-information slot moves from f_g 0.50 to 0.75. Replacing the half with the arm removes the ladder's zero-gDNA false
gDNA (−97 %) but costs every capture-ON stratum (stranded × ON object error ×1.04–1.24 across panels, the deferred
stratum ×1.4–3.8): the extra half props enriched exons up against one landscape fitted across every class. A repair
needs a continuous per-object partner, and the log-odds bracket must cover the tail. `simplex_logodds.py`'s module
docstring says the arm replaces the half; the code adds it.

### psi-reads-kappa-where-the-strand-channel-is-dead
`priority: later · stage: calibration · kind: defect`
Where the protocol's Bayes factor calls the strand channel dead, the discriminability gates only `strand_evidence`
(`solve_kernel.cpp`); ψ's grid still reads the estimated κ̂. At od = 0 κ̂'s sampling error gives the strand term a slope
n(½ − κ̂), so real overdispersion reads as composition (SD(κ̂) ≈ 0.08 at about 40 spliced reads). On a test-chromosome
panel with planted gDNA strand overdispersion 0.05, ss 0.50 g98 (κ̂ = 0.515), it costs +1,178 in scope at capture OFF and +29,106 on the deferred capture-ON row. Candidate: ψ reads κ = ½
where the channel is dead, a pinned A/B on the unstranded rows (`policy_benchmark.py --by-class`).

### the-gdna-landscape-collapses-at-low-depth
`priority: later — re-measured on HEAD before it is ranked · stage: calibration · kind: defect`
On sparse real libraries calibration's gDNA has fallen with depth: cfRNA-B (stranded, captured) 0.11 / 0.45 / 0.83 gDNA
per deposited fragment at 1 % / 10 % / full depth, the VCaP mix P/O 0.125 / 0.52 / 0.92, on an older tree; cfRNA-B's EM
followed. The refits' landscape carried 63–86 % of it: a slot trains as located only when it holds a gDNA fragment, so at
low depth almost every kernel is a zero-count anchor and the refit converges to one mode at the grid floor, against which
the reader's efficiencies then rise. Shallow objects also lose gDNA to the plug-in, which conditions an object's gDNA
fraction on its own fragments (−21 / −11 / −2 % at 1 / 10 / 50 fragments; the leave-one-out posterior is exact). On the
test chromosome's 1/10 and 1/100 depth panels the preview found no blow-up. A repair keeps the mixture's implied total
equal to the pooled count, makes a fit at depth f the full fit shifted by log f where the data resolve it, and keeps the
zero controls. Instruments: `calibration_vs_oracle.py --suite` on the depth family; the VCaP mix at 1 % and 10 %.

### the-bulk-can-be-calibrations-nascent-rna
`priority: later · stage: calibration · kind: defect`
On a small, weakly stranded library with heavy nascent RNA, calibration books the genes' nascent RNA as gDNA, and that
population can outweigh the true intergenic gDNA and become the bulk every object is read against. On the
`combo_extreme` golden (κ 0.65, no capture, 13 training objects) the genic deconvolved gDNA, 0.29 per bp, is the bulk,
the true intergenic density, 0.06, reads 0.2, and transcripts are 36.8 % wrong; given the true gDNA deposits the reader
is right, so the origin is calibration. How often it happens is not established.

### the-transfer-price-arithmetic
`priority: later — each a one-line kernel change, prototyped in C++ and A/B'd · stage: calibration · kind: defect`
Three width-only inconsistencies in the transfer kernel's prices (`transfer_kernel.h`), each moving few fragments:
- an empty piece writes its junctions' witness opportunity in whole-strand units and the first full node past it
  re-prices it in `lane_receive` with no `κ_read`, so at κ = ½ the hop pays log(κ_read)² = 0.48 nats² of false
  disagreement; write it as `a_sum / κ_read`. The empty-piece gate in `test_transfer_rna_lanes.py` runs at κ = 0.99,
  where the forms coincide, so test first at κ = ½;
- the level rule's width divides by a mode built on `n_u + n_s` while its value reads `n_u / a_g`, inflating `f_pred` by
  1 + n_s/n_u (at most 5.6 fragments of 1.1 M); build the mode on the density;
- the gDNA lane prices a total over one component's opportunity (`hop_price(n_u, a_g, n_u, a_g)`), so equal component
  densities differ by the composition term: negligible at equal lengths, up to about 6 nats² on the fl-gDNA-long panel
  (124.7 fragments of 1.15 M). Either rule that only values must be densities, or A/B a composition-free total.

### the-transfer-policy-loses-to-silence-on-unstranded-capture-off
`priority: later · stage: calibration · kind: defect`
The transfer policy ships because it must win on unstranded libraries, yet on the ladder's unstranded × capture-OFF rows
it loses to the silent policy (whole-library |gDNA − truth|: g05 ×1.02, g50 ×1.07, g98 ×1.31), while it wins unstranded ×
capture ON about 4× and every contaminated stranded row (×0.61–0.97). The likely source is the removed intron factory's
local composition information. There an intron region's own solve has no strand channel, so its gDNA fraction is
whatever the landscape and the messages deliver, and which of the two should decide it is the open question. The
transcript-level cost is unmeasured. Instrument: `policy_benchmark.py --panel ladder --by-class`.

### the-strand-frame-under-the-shared-law
`priority: later — its prototype fails at weak strandedness · stage: pre-EM · kind: gap`
The shared length law (`calibration/shared_length_law.py`) reads no strand. Where gDNA and RNA lengths truly differ under
capture it keeps only the difference capture cannot make, so on the fl-gap side panels the stranded gDNA-long capture-ON
rows regress (genes +4.4 / +16.0 / +23.5 % at g05 / g25 / g50 ss0.99), where a strand frame reads −17.5 / −43.8 / −50.2 %.
The frame splits one-strand exons by strand, where capture cancels; its prototype invents the minority origin at weak
strandedness (gDNA-long g05 ss0.70 off capture: gDNA's law 116 bp against 259, transcripts +85 %). Open: an unbiased
moment estimator of gDNA's per-length count at one-strand exons, `N_l − (S_l − N_l/2)/(p − ½)`, with its exact binomial
variance, then the A/B again, the small-gap panels first. The equal-length ladder cannot see it. Prototypes:
`$RIGEL_RUNS/prototypes/2026-10-27_four_arm/`.

### the-pooled-q-in-the-gdna-count
`priority: later — calibration's own count is the larger half, a separate modelling project · stage: pre-EM · kind: defect`
`priors.assemble_priors` converts each boundary's gDNA mass to fragments by `q = boundary_mass_per_crossing`, pooled over
gDNA and RNA, so where RNA dominates a boundary the locus's gDNA count is over-stated even under a perfect calibration.
That is the small half: across the ladder O − S reads +11.3 k off capture and +47.0 k on (+0.13 % / +0.31 % of the
stratum's gDNA). Calibration's own error P − O is larger and of the other sign (unstranded capture-OFF −365.6 k, g98 P/O
0.945). It arises in the solve under the fitted landscape: at g98 ss.50 OFF in introns, as a floor of about −4.7 % that
does not shrink with size; at g50 ss.50 ON partly in the refit training on its own under-called captured exons. A
plausible, unestablished mechanism: ψ's composition arm reads the landscape with no term for the slot's own count noise.
gDNA's own share in the count alone improves capture-ON pools (stranded |gDNA| −26 %) and worsens capture OFF (+2 to
+7 %), so it lands only with the count repaired. On the STAR rung the EM's share was read 4–7 % too high until
2026-11-14, which hid calibration's under-count: at g98 unstranded capture OFF gDNA called nascent RNA read +70 k then
and reads +302 k now, the simulator's +292 k. Instrument: `prior_vs_oracle.py` (P, O, S).

### pruning-uses-three-references
`priority: later · stage: pre-EM · kind: defect`
`scoring.cpp` prunes annotated candidates against the best annotated one and synthetic spans against the best span
(`emit_pruned`), and the locus gDNA term against the best kept RNA candidate (`gdna_competes`). A span and the gDNA term
read the same footprint and under equal length laws differ by log(2·p_sense), yet the span survives any distance below the
annotated reading while gDNA is cut at 9.21 nats, so an implicit unit with a long gap can go to a transcript or a span but
never to gDNA. Candidate: one reference, the unit's best candidate of any kind, A/B'd per stratum. The routing tests set
`is_synthetic` False everywhere, so pool separation is untested. The intergenic background is never pruned.

### the-pseudocount-strength-is-not-derived
`priority: later · stage: pre-EM · kind: question`
The EM's two pseudocounts carry one pseudo-fragment per unit with a gDNA candidate (`pipeline.em_pseudocounts`), and
nothing derives that strength. Under MAP it cancels at the neutral point; under the shipped VBEM it does not, and the
grouped update lets it move the split within RNA (a 6-fragment isoform reads 4.0 / 5.0 / 5.8 at `P_R` = 0 / S / 10S).
The exact VB update for the same prior holds 4.0 at every strength with no constant (a toy in
`$RIGEL_RUNS/prototypes/2026-09-24_spliced_rule/`). A prior is needed: with none the EM's gDNA error is −731.8 k at g98
ss.99 ON.

### per-transcript-prior-lane
`priority: later — no candidate mechanism · stage: pre-EM · kind: gap`
The EM has no per-transcript RNA prior: `rna_prior_weight` is an argument of
`AbundanceEstimator.run_batch_locus_em_partitioned`, passed to `em_solver.cpp`; `pipeline.py` never passes it, and the
solver falls back to the EM's own belief. Handed the true weights
(`quant_accuracy.py --arm oracle_alloc_seed`, an oracle intervention), stranded capture-ON transcript error at g50 falls
5.50 → 2.07 % of the true RNA. Most of that sits in the proportions among mosaic isoforms that share every object; a
coverage-derived weight is far worse than the EM's likelihood (5.21 → 53.37 %). The lane and `warm_start='prior'` serve
the oracle arm.

### scoring-penalties-are-underived-constants
`priority: later · stage: pre-EM · kind: question`
The overhang price, log(`overhang_alpha`) = log 0.1 per base, RNA only (`config.py`), is underived, and it ties with the
pruning bound exactly: four overhang bases cost −log 10⁻⁴, so a candidate four bases behind its pool's best sits on the
bound and `-ffast-math` rounding decides its side (a scorer refactor moved 9,705 such candidates in and 8,906 out of 92 M
on one STAR-rung library). The splice blacklist's minimum count of 2 is underived too.

### ruler-witness-geometry-on-transcript-panels
`priority: later · stage: pre-EM · kind: defect`
Only a transcript holding the junction a probe spans binds that probe whole, so junction-spanning capture is
isoform-specific. Two probe panels with identical genomic probe parts have identical exact gDNA efficiencies, yet
their RNA yields differ by 7.1 % when the parts belong to junction-spanning baits in one panel and separate baits in
the other. Thus gDNA alone does not identify junction-specific RNA capture. The owner's real panel includes
junction-spanning baits (confirmed 2026-11-10); capture inference must nevertheless remain probe-free (`DESIGN` §7.2).

Neither junction-cap choice is universal. On the junction-only test panel, removing both reader and junction caps
puts every class's length/yield ratio within 2 % of one and improves gene error 33–49 %. On the panel without junction
baits, it makes isoforms 7–8 % too long and increases gene error 41–366 %. A common class correction also fails;
within-gene length differences matter. The observation that capping flattens useful RNA readings describes these
tested layouts, not a rule that every partially or sparsely probed transcript has one capture level.

The remaining data-derived route requires RNA-specific evidence. A feasibility check using true annotated-RNA
counts tracks junction/exon capture ratios across the three designs but remains 10–20 % low. It uses selected exon
regions feeding one junction and at least 50 true fragments on each side; it has not been demonstrated using
observable inputs without that selection. The maximal-reach junction opportunity is a suspected source of the low
reading. Spliced counts certify RNA but do not alone separate capture from abundance or isoform mixture. No
per-junction RNA-witness estimator has been accepted; the junction cap landed on 2026-11-13 is one library level read
from every junction's own footprint (`ISSUES: the-junction-cap-reads-bookings`). Receipts: `$RIGEL_RUNS/prototypes/2026-11-10_length_gap/`, sections 2 and 5–7.
Instruments: `ruler_vs_truth.py`, `quant_accuracy.py --arm oracle_ruler` (the older unanchored arm is described in
`ISSUES: instrument-ledger`).

### the-em-stops-before-it-converges
`priority: later · stage: EM · kind: defect`
The locus EM (VBEM with SQUAREM, 1,000 iterations, δ 1e-6, `em_solver.cpp`) runs out of budget on degenerate loci —
near-identical isoforms, antisense overlaps, capture clusters, exact paralogs — and its answer there is where it stopped.
The accelerator is rounding-sensitive: a common length factor of 1 + 1e-12 moves 40,860 fragments over the test panel's
30 conditions, under 2e-8 with ordinary updates. On aligned input a 1e-5 upstream perturbation can move tens of thousands
of fragments between exact paralogs (SIK1 and SIK1B), part of the g00 rows' instability. Converged, the ladder's
capture-ON rows move only −0.13 to +0.12 points. The repair is explicit regularisation of the poorly determined
directions with a stable solver run to convergence, decided on the ladder.

The cap is SQUAREM's iterations // 3 = 333. At ladder g05 stranded × ON on the simulator it stops the 34–36 largest
loci, a quarter of the EM units, and running them to convergence (cap × 100, same δ) adds 0.6–2.5 k fragments of
transcript error under each of three gDNA priors. The cap also widens differences between priors. Two priors 0.02 % apart
differ by 2.4 k fragments of transcript error at the cap and 0.9 k converged. At LARGE1 one stopped on a flat ridge, 1.6 k
fragments of isoform error for 0.3 nats; uncapped, it reaches the other's split. At CLTCL1 the converged split moves
1.1 k fragments with a 3-pseudocount change of the gDNA prior, and δ six decades tighter moves neither solution. Whether
those are two basins is unproved: that needs both starts under one prior and the VB objective compared. So a transcript
difference of about 2 k between two prior mechanisms at g05 is mostly where capped solves stopped
(`$RIGEL_RUNS/prototypes/2026-11-09_g05_isoforms/`).

### the-x86-e-step-is-scalar
`priority: later — with the performance session · stage: EM · kind: build`
The E-step's exponential runs NEON on AArch64 and scalar on every other machine, so an x86-64 build (a Linux cluster's)
takes the slowest path. The AVX2 and AVX-512 kernels were deleted because no build compiled them: the portable build that
ships targets baseline x86-64, and the native build that would compile them was never run on x86-64. They are in f871b669
(`native/fast_exp.h`, `native/simd_detect.h`, and the E-step loop, then in `native/em_solver.cpp` and now in
`native/em_estep.h`). Reviving them needs an x86-64
machine and dispatch at load time, since one wheel must run on every x86-64: each kernel built for its own ISA and chosen
by the CPU it finds. Each must pass `test_fast_exp.py`, state how far it moves the EM's numbers against the scalar path
(an FMA polynomial does, as NEON's does on AArch64), and be timed in back-to-back pairs.

### instrument-ledger
`priority: later — (a) and (b) before the capture-reader session's A/Bs · stage: instruments · kind: instrument`
What the instruments do wrong or leave out, each fix its own commit with its `--self-test` or a gate:
(a) The paired proper score the preview was judged by (Brier and its calibration component on identical read IDs) lives
outside the tree; no in-tree instrument reports it.
(b) `quant_accuracy.py --arm oracle_ruler` is wrong under capture (g98 ss.99 ON 475.87 % of the transcript truth): it
scales the simulator's factor so the fully probed class reads 1, while the shipped efficiencies read 1 at the uncaptured
bulk. Re-anchor on that convention and prove it on the capture-OFF identity and one capture-ON row.
(c) `calibration_vs_oracle.py` has no arm carrying the true capture efficiencies and no gDNA-only object row.
(d) The gdna-field-uniformity gate (`calibration_oracle.py`) z-scores every slot as an independent Poisson draw, but the
boundaries of a run of short regions share one draw, so one high draw fails a row; seven capture-OFF test-chromosome
rows certify composition only. Z-score such a run once, on its total.
(e) `slot_truth.npz` counts boundary gDNA per crossing, so under capture its sums exceed the library's true gDNA, while
`policy_benchmark.py` calls its Σ|Δ| "in fragments".
(f) No instrument output records the code revision, the `--set` overrides or the prototype arm.
(g) `build_test_reference.py` puts `src` first on `sys.path` unconditionally; the ladder-report skill hard-codes the four
ladder rungs; `panel.py --index` does not reach `simulate`.
(h) A refactor's no-op proof has not included a real library. The code pass (5341ea06) was proven on three ladder
conditions; on the VCaP mix its transcript table differs from 81fc67b6's on 8 transcripts (30 of 13.2 M fragments,
exchanged between isoforms of four genes; one zero-count transcript changes locus), with every payload and calibration
array identical. Which change did it is not established. The references of 2026-11-13 hold the VCaP mix
(`TESTING.md` §6): check it at every stage.

### debug-capture-memory-is-unbounded
`priority: later · stage: instruments · kind: instrument`
`calibrate(diagnostics=…)` fills one `SweepDiagnostics`, overwritten by each sweep, with no size or region bound; on the
whole human genome a record per sweep reached 25 GB in 20 s, and the one record is still unbounded. `ruler_vs_truth.py` always passes `diagnostics` and
does not refuse a real index. The fix: `diagnostics` refuses an index whose manifest does not mark a panel, and the
record is bounded to named regions. Until then one real-genome job runs at a time.

### latent-defects
`priority: later · stage: several · kind: defect`
Correctness defects found by review, each its own commit. MOVES A NUMBER:
- *index:* a `gene_id` with transcripts on two references takes its reference and strand from one and its span across both.
- *scan:* the splice census counts placements, so a multimapper's secondary placements enter `splice.*` (an NH:2 fixture:
  20,500 for 20,000 fragments); a library aligned without a motif tag loses the annotated junction strand the resolver
  already looks up.
- *accumulator:* a first intron starting at the fragment start makes `build_segments` emit no leading segment and
  `left_sj` mislabel every junction end; fix the executable specification first, then re-cache.
- *calibration:* "has its own composition" is τ > 0 in the transfer kernel and τ > `OWN_EPS` (1e-9) in `solve_kernel.cpp`'s
  landscape-training predicate, so a slot with 0 < τ ≤ 1e-9 sends a claim but trains as composition-free.
- *pre-EM:* the coverage weight's plateau is min(f, L/2) (`scoring.cpp`) where enumeration gives min(f, L − f) (L = 260,
  f = 250: 1.93 against 1.02), and its gate restates the formula; it acts through the EM's warm start.
- *EM:* `assign_posteriors` weighs `log(theta + 1e-300)`, so a dead component keeps a posterior of about 1e-300: a
  fractional run writes counts of about 1e-299, and sample mode can hand a fragment to a dead component. Give it −inf.
- *post-EM:* sample mode's `posterior_mean` is Σp²/Σp, and NaN for a transcript with no EM unit, which the MANUAL does
  not say.
- *instruments:* `lift_choices` keys records without the introns and hypotheses the native key holds, so tied records
  that differ in introns pool and oracle truth can replay the wrong one.

LATENT, no current path reaches it: `--scan-buffer-size` does not bound the native queues; `splice_graph._ref_slices`
drops exons of unsorted input silently; a junction efficiency clipped at 0 can leave an `em_effective_length` of about
1e-304; a capture entry with neither `enabled:` nor probes yields an empty `CaptureConfig`; the scan cache
(`scripts/design/_scan_cache.py`) writes in place, not temp-and-rename, and `_payload_from_parts` skips `from_dict`'s
checks; an arm injecting `rna_sense_frac` without `n_rna_obs` kills the strand channel; `strand_variance`
(`transfer_rows.h`) carries `(n·f)²` where `n(n−1)·f²` is right, which matters only at od > 0.

### performance-memory-bounded-solve
`priority: later — re-timed before it is ranked · stage: several · kind: build`
A deep run must be fast enough to iterate on, and memory-bounded. The baseline (VCaP, 18.6 M fragments, eight threads,
on an older tree): 125 s — calibrate 39.5, quant 33.7, the scan 28.4, the second pass 13.0, the index load 6.4; peak RSS
10.3 GB, in quant. Ranked by a cProfile share: the capture effective lengths (quant's 5.2 s); the second pass's factor
combiner, 355,180 calls on 2–4-element arrays; the landscape fits, about 5 s; the sweeps' 24.0 s, where the blur's loop
interchange saves 1.5 s at the price of a summation-order change (the owner's call). The scan and the index load are at
their floor. Numeric no-ops left, each proven with `rename_identity.py --bam`: `build_region_partition_arrays` runs 3+
times per run and `build_sj_arrays` twice; the EM clears `local_map` per locus and allocates per E-step; the solve builds whole cubes
for slots it discards.
Instrument: `profiling/profiler.py`.

### hygiene-ledger
`priority: later · stage: several · kind: hygiene`
What the reviews left, each its own commit, proven a numeric no-op unless marked:
- Recomputed outputs: `loci.feather` rebuilds the annotated count as max(rna_total − count_nascent, 0) + count_unambig
  behind a clamp (`estimator.py`), and the nascent table re-codes the TPM denominator; removing them moves ulps.
  `gdna_opportunity.contained_opportunity` and `sj_opportunity._ramp_sum` compute the same Σ(ℓ − w + 1)₊ with different
  float arithmetic; one home moves ulps in the gDNA length law.
- Constants: the simulator's default seed, `CONSTANTS.simulator.scenario_seed`, is restated as 42 in `ReadSimConfig`,
  `MutableGenome` and `SimulationParams`.

## OPEN — parked

### the-edge-density-floor-under-capture
`priority: parked — an approximation; the measured alternative was not selected · stage: calibration · kind: limit`
The transfer layer's EDGE row reads a boundary's structurally pure-DNA count as a one-sided lower bound on the adjacent
exon's gDNA density, and under capture that ordering fails: an exon [0, 2000) with a probe on [0, 100) and fragments of
200 puts the edge's density at 25.9× the exon's average for every nondecreasing binding function. Withdrawing the row
trades useful imputation for less bias (the deferred test row g50 ss0.50 ON region error 3.69 → 43.47 %; a weak-strand
OFF fixture's intron error 61.4 → 6.6 %).

### the-atom-at-an-unwitnessed-both-strand-slot
`priority: parked — a limit of the information · stage: calibration · kind: limit`
On stranded data at a both-strand slot whose minor strand has no witness (a single-exon gene wholly inside another gene's
exon), the tilt atom admits "all of this slot's RNA is on one strand", and the strand split cannot tell "pure with gDNA"
from "both strands without", so it pushes toward false gDNA: the golden `antisense_contained` reads 200.8 false gDNA
fragments in a gDNA-free locus. On the ladder, where the landscape trains on about 15 k regions, the same slots cost a few
fragments each; the landscape's location floor is the lever. On unstranded data the atom is inert.

### two-sided-exon-row
`priority: parked — waits for an enrichment witness · stage: calibration · kind: defect`
On unstranded data an exon's held row is the intron's composition through the face map, whose upper side is the map's
plateau above its ceiling — a lower bound on gDNA — so at pass zero an unstranded licensed exon reads about 9× its true
gDNA, and forwarding it compounds the bias: 4–6 % of the ladder's unstranded error, 20–32 % on the test chromosome, where
transfer is worse (2,647 → 4,056 at g50 ss.50 OFF). Every cap, two-sided row or wall tried fails where junction probes
enrich the flux more than the crossing. A two-sided level needs a witness that this library's gDNA is enriched (silent
genes' exons against the intergenic density, the landscape's job); derive that witness first. `policy_benchmark.py
--by-class`.

### the-lower-bound-noise-ratchet
`priority: parked — waits for the same enrichment witness · stage: calibration · kind: defect`
A level from an RNA-rich node's own strand profile has a mode that is noise around zero; its lower side bounds its
neighbours and the tightest noisy neighbour wins (ladder g05 ss.99 OFF 44,714 → 45,076; test chromosome g00 ss.99 ON
64 → 216, most of it at `capcluster_ab`'s inner termini). A two-sided own-profile level keeps fewer fragments overall, so
lower-only stays. The gDNA lane's edge level does the same: a shallow single-strand flank (truth 0.530) reads 0.596 under
its intergenic neighbour's Poisson level. The cure is the witness `ISSUES: two-sided-exon-row` waits for. Instrument: the
ladder's g00 rows, `calibration_vs_oracle.py`.

### message-layer-open-cases
`priority: parked · stage: calibration · kind: question`
Residual questions in the message layer, none a hole: (a) an exon with both faces speaking sums the arrivals, and the
price where both carry the same intron's claim is unmeasured; (b) the chain of termini, an empty outside piece (median
12 bp) whose far face is another terminus, is half the ladder's terminus-boundary error; (c) whether a hop's transfer
variance (`hop_price`, `transfer_rows.h`) prices a ratio built on a handful of counts: TRANSPORT blurs by
`trigamma(n_u+½) + trigamma(n_s+½)`, about 4.93 at a zero certified count, and SPLICE_OUT's rate sd is `trigamma(ΣJ+½)`
where the delta method gives ΣJ/A² over (ΣJ/A)²; (d) whether the flux level's floor prices its own dispersion.
`policy_benchmark.py --by-class`.

### capture-premise-untested-on-cdna
`priority: parked — no library in hand can test it · stage: calibration · kind: risk`
The contracted length reads the panel from gDNA and applies the same efficiencies to cDNA, assuming off-target cDNA is
depleted as off-target gDNA is; the simulator satisfies that by construction, and cross-hybridisation or nonspecific
binding of cDNA could make it far milder. Its exposure is the isoform split: on VCaP 4,677 of 12,039 multi-isoform genes
with ≥ 20 fragments changed dominant isoform between the contracted and the plain length (an older length rule). The
experiment that settles it
is one sample sequenced both ways, panel and whole transcriptome. `count_unambig` beside `count` tells a user which
isoform assignments rest on shared fragments alone.

### capture-degeneracy-standing-risk
`priority: parked — watch · stage: calibration · kind: risk`
gDNA's two-pool length contrast survives capture by a degeneracy: the shared-contaminant assumption is false under
capture, and it is safe only because the intergenic pool is depleted, not impure, so `a_0` clips to 1 and the algebra
collapses to `g = f_0` (`fl.deconvolved_gdna_counts`). A panel that put RNA into intergenic space would break it silently.
The one-sided off-target rate (`gdna_density.one_sided_rate`, fitted on intergenic and intron counts) leans on the same
purity: RNA spread over intergenic space lifts it (a toy fixed point 202 bp against 217), and intron regions beside probes
are capture-enriched because probes overhang exon edges (g05 ON on the test chromosome: pure intron seeds at a median 1.89× the
off-target rate, 99th percentile 44×). Neither is measured on real data.

### binary-cuts-on-continuous-quantities
`priority: parked — tracking; each raised when its mechanism is next touched · stage: several · kind: tracking`
The rule is no yes/no gate on a continuous quantity, and these still cut: `calib_refit_iters = 3` stops unconverged (an L1
change of 85–34 k fragments at a fourth refit); a kernel with count < 1 gets no location; pruning on a k·ln10 lattice; the
zero-density veto, which hits 27–43 % of held records; one unit merges two loci; `TAIL_DECAY`; `max_frag_length = 1000`,
in three roles; the splice blacklist's count ≥ 2 and anchor ≤ the longest seen. An inventory ranking 55 such cuts by harm
is `$RIGEL_RUNS/prototypes/2026-09-24_cut_inventory/INVENTORY.md`.

---

## CLOSED / REFUSED

Every closed and refused entry, with the number that closed or killed it, was archived at the 2026-11-03 reset:
`$RIGEL_RUNS/archive/cleanup_2026-11-03/docs/ISSUES.md`. Read it before re-proposing a mechanism.
