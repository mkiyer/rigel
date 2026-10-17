# Review brief: the capture weights' noise and the unstranded × capture-ON stratum (2026-10-11)

> **ERRATUM (2026-10-11, after the review it was written for — `CAPTURE_IMPLEMENTATION_REVIEW.md` in this
> directory).** Four statements below were refuted by measurement and stand corrected there: (1) §1's claim
> that calibration uses fragment-length *composition* evidence — the direct channel is unbuilt (no calibration
> source reads the inverse-length banks; zeroing them changes every calibration array by exactly zero); the
> "composition row" in the reader's likelihood is the solve's own gDNA:RNA composition message, not a length law;
> (2) §4.3's diagnosis of *shallow Poisson noise* — the extreme weights are rare mode switches created by the
> delivered composition and DNA-level factors at objects of any depth (696 of 70,176 log weights above 1 on the
> flat stranded field, depths to 3,769, median 6); (3) §4.3's open asymmetry — the unstranded rows are flat
> because the reader reads essentially nothing there (maximum log weight 0.00005 against 13.1 stranded), not
> because they are robust; (4) §5.3's exon-to-intron contrast as an enrichment witness — enriched gDNA
> (700, 1) and flat gDNA plus exon RNA (1, 1, 699) have identical likelihoods under a dead strand channel and
> equal length laws. The width-based local shrinkage of §4.4 was built as the review's candidates 1–3 and
> refused: it repairs the flat rows and breaks the captured gap rows. The review's own next step is shared
> capture-field uncertainty inside the allocation, derived there and not yet validated.

*Written for an external reviewer who knows RNA-seq quantification and Bayesian inference but not this
codebase. Sandbox document (`docs/dev/`): nothing here is authoritative; the rulings and derivations it
summarises live in the nine permanent docs under `docs/`. State described: `main` at `0fdb9dce`. Every number
below was measured in the session that wrote this and is reproducible with the commands in Appendix A.*

**What we ask of you.** Two accuracy problems are open, each with a measured attribution and, for the first,
two repairs that were derived, prototyped and refuted by measurement. We want a critique of the diagnosis,
of the proposed lines of investigation, and of anything we are not seeing — not a rewrite of the tool. The
owner's constraints in §6 are fixed; proposals that violate them will not be taken, however well they would
score.

## 1. The tool, in one page

**The problem.** RNA-seq libraries contaminated with genomic DNA — cell-free RNA from plasma above all, where
gDNA can be a quarter or more of the fragments and hybrid-capture panels enrich both populations — must be
quantified at the transcript level with the gDNA removed, not merely flagged.

**The axiom.** There are three populations and no fourth: gDNA, RNA on the plus strand, RNA on the minus
strand. "Nascent" is not a population: RNA inside an intron is RNA that has not spliced at that position. A
spliced fragment is certified RNA and needs no deconvolution; an unspliced one is the whole problem. The
annotation index manufactures one synthetic "nascent" entity per multi-exon gene so that unspliced RNA has
somewhere to go in the transcript table; the simulator plants nascent RNA at a 20 % share, a stress level,
where real libraries run near 4 %.

**Capture is a spectrum.** There is no capture-ON/OFF status for a library and no detector of one: every
object carries its own continuous capture weight, and "no enrichment" is the limit in which the weights are
all equal. This is an owner ruling with evidence behind it (plasma panels enrich scarce transcripts that stay a
tiny fraction of the library, so any library-level statistic reads the uncaptured background).

**The pipeline.**

1. *The scan.* A single pass in C++ over a name-sorted BAM resolves every fragment to the transcripts it is
   compatible with and deposits it on *objects*: regions (contiguous genomic pieces at exon/intron
   resolution) and boundaries (the cuts between them), with a column per strand. Fragments that are ambiguous
   between placements or multiply mapped are held and assigned in a second pass by the abundance of the
   unambiguous ones.
2. *Calibration.* Per object, a solve infers the gDNA density and the RNA amounts from the strand columns
   (gDNA is strand-symmetric; RNA follows the protocol's measured strand specificity), from the fragment-length
   composition evidence (gDNA and RNA have different length laws, learned from the library), and from messages
   passed between neighbouring objects along a gene (an exon's gDNA level is imputed from its neighbours). A
   *landscape* prior — a density over log gDNA density across all objects, fitted on the solve's own
   deconvolved gDNA and refitted three times — regularises every object. The outputs are per-object
   gDNA/RNA pies and the library's gDNA fraction. Since 2026-10-09 the solve carries no separate
   intron-background constraint ("the factory"); its messages are observational (the strand likelihood of the
   observed columns), and in ambiguous cases it allocates conservatively, gDNA over RNA.
3. *The capture reader.* For every object the reader evaluates the likelihood of its two strand columns as a
   function of its gDNA density, with the RNA amount integrated out under the factors the solve delivered
   (the composition row, the held neighbour levels), and takes the posterior **mode** under the landscape:

       log L(x) = log ∫ Pois(u; ρEg/2 + q·r) · Pois(v; ρEg/2 + (1−q)·r) · C(log(ρEg/r)) · H₊ · H₋ · r^(−1/2) dr,  ρ = eˣ
       x* = argmax_x [log L(x) + log P(x)],   w = exp(x* − median x*)

   `u, v` the strand columns, `Eg` the object's gDNA opportunity (its length-like exposure), `q` the strand
   probability for RNA, `C` the composition evidence, `H` the neighbour levels, `P` the landscape. The weight
   `w` is relative to the typical read object (the median), and the downstream EM is invariant to that unit.
4. *The EM.* Per locus, a variational EM over the annotated transcripts, the synthetic nascent entity and a
   gDNA component, with the gDNA prior from calibration. Each component's effective length is contracted by
   capture: the sum over the objects its fragments deposit on of its conserved share there times that object's
   weight, `L_T = Σ_pieces S_T(p)·w_p + Σ_cuts M_T(k)·w_k`. Benchmarks read the fractional assignment; production
   draws a sampled one.

**What the weights do.** On a captured library a probed exon reads a weight of hundreds relative to the
typical object, so a transcript whose exons are probed gets a longer effective length and the EM does not
mistake enrichment for expression. On a capture-OFF library the true weights are all one level.

## 2. How performance is judged

- *The ladder.* Sixteen simulated libraries of 10 M fragments: gDNA at 0, 5, 50 and 98 % × strand
  specificity 0.99 ("stranded") or 0.50 ("unstranded") × capture ON or OFF. gDNA and RNA are given **equal**
  fragment-length laws on purpose, so the length channel cannot separate them and calibration's other
  channels are what is tested. Three strata are in scope for the release; **unstranded × capture-ON is
  deferred** (reported, never a development target) because calibration is structurally blind there: no strand
  channel, and capture enriches both populations alike.
- *The gap panels.* Two side panels at 50 % gDNA with a fragment-length gap in each direction (RNA longer,
  RNA shorter); four conditions each.
- *The test chromosome.* A small hand-designed chromosome for seconds-long development runs.
- *Real data.* The VCaP mix (an exome DNA library and a transcriptome library of the same cell line, mixed;
  read names carry the truth: gDNA fraction 0.2518), and two plasma cfRNA libraries without truth.
- *Metrics.* Per condition, the three pools (gDNA, nascent, annotated) against the simulator's realised
  per-fragment truth; the transcript-level Σ|estimate − truth| over every transcript as a share of annotated
  truth; the same at gene level; the library's gDNA fraction. Everything is read **per stratum, never pooled**.
- *The owner's guide (2026-10-08).* Roughly 5 % gDNA/RNA separation error, about 5 % transcript error and
  under 2 % gene error would be very good; robustness and generalisability across input distributions,
  especially fragment-length distributions, are the objective, not aggregate simulated accuracy.

## 3. Where the tool stands (main `0fdb9dce`)

The ladder, transcript Σ|Δ| as a share of annotated truth with the gene level in parentheses, the arm scored
on 2026-09-30 against today's:

| stratum | 2026-09-30 | today |
|---|---|---|
| stranded × capture OFF (in scope) | 1.99 % (0.21) | **7.51 % (0.53)** |
| stranded × capture ON (in scope) | 6.27 % (1.62) | **9.71 % (1.87)** |
| unstranded × capture OFF (in scope) | 2.01 % (0.24) | 2.07 % (0.24) |
| unstranded × capture ON (deferred) | 10.20 % (3.50) | **40.42 % (15.40)** |

Per condition, the rows that moved (transcript %):

| condition | 2026-09-30 | today |
|---|---|---|
| g00 stranded OFF | 1.94 | 13.61 |
| g05 stranded OFF | 1.57 | 2.71 |
| g50 stranded OFF | 2.36 | 4.12 |
| g00 stranded ON | 6.46 | 15.10 |
| g05 unstranded ON (deferred) | 10.72 | 48.92 |
| g50 unstranded ON (deferred) | 10.73 | 65.16 |
| g98 unstranded ON (deferred) | 119.24 | 674.97 |

Every other row is flat. Between the two arms four things landed: the junction price's clip (2026-10-05), the
count foundation without the intron factory (2026-10-09), the detector-free capture reader (2026-10-10) and the
multimapper second pass (2026-10-11). Real data moved the right way: VCaP's gDNA fraction 0.2395 → 0.2528
against a truth of 0.2518 (the multimapper pass, +0.0127); the plasma libraries read 0.100 and 0.165. On the
gap panel with shorter RNA, the stranded × OFF row's false-capture catastrophe went from 42.3 % to 5.0 %
transcripts when the reference-based reader was replaced by the detector-free one. The pipeline runs the
20 M-fragment VCaP mix in 253 s with the scan pinned to one thread.

## 4. Issue 1 — the transcript (isoform) regression on the stranded strata

### 4.1 Attribution

The worst in-scope row, g00 stranded OFF (zero gDNA, capture off), was run as a pinned pair with everything
identical except the capture weights forced to one:

| g00 stranded OFF | transcripts | genes | mass assigned to silent isoforms | calibration's gDNA estimate |
|---|---|---|---|---|
| as shipped | 13.61 % | 0.83 % | 251,216 | 1,126 fragments |
| weights forced to 1 | 1.93 % | 0.14 % | 17,058 | 1,112 |
| the 2026-09-30 arm | 1.94 % | 0.14 % | 17,129 | 125 |

The weights own the whole of the row's regression; calibration does not move. The gene level holds, so the
damage is **isoform allocation inside genes**: with the weights all truly one, the published ones are each
object's own Poisson noise around that level, largest where objects are shallow; two isoforms that share most
objects but differ in a few shallow ones get random effective-length ratios, and the EM moves mass to silent
isoforms whose contracted length happens to look favourable. The price was known qualitatively when the
reader landed ("the price of no detector") but its size on the ladder was not: the terminal report excludes
the zero-gDNA rows, which is where it bites.

### 4.2 Two repairs derived, built and refuted

Both replace the plug-in mode by a posterior average of the same grid posterior the reader already holds
(no detector, no reference, no constant), published in the same median unit. Both were A/B'd pinned and
fractional, both arms from one build, on the test chromosome, the ladder, both gap panels and VCaP.

| readout | g00 stranded OFF | ladder stranded OFF | stranded ON | unstranded OFF | RNA-long stranded ON | RNA-short stranded ON | VCaP gDNA fraction (0.2518) |
|---|---|---|---|---|---|---|---|
| the mode (shipped) | 13.61 % | 7.52 % | 9.71 % | 2.07 % | 6.09 % | 14.24 % | 0.2528 |
| E[eˣ], the expectation of the density | 29.94 | 15.57 | 14.11 | 11.26 | 7.57 | 15.39 | 0.2489 |
| E[x], the mean of log density | 9.50 | 7.86 | 13.00 | 4.16 | 12.42 | 23.83 | 0.2412 |

The weights' quantiles on VCaP say why. Under E[eˣ] the 10th/90th percentiles went from 0.37/425 to
0.009/6.6: under a posterior as broad as the prior, the expectation of eˣ is the prior's high tail, and every
shallow object was lifted above the deep ones. Under E[x] the 10th percentile went to 0.09: the mean of a
one-sided posterior sits below its mode (0.14–0.79 nats at a zero-read object, measured earlier), so shallow
objects lost enrichment they have and every capture-ON stratum paid. The one row that improved under E[x] is
the one whose noise is symmetric.

### 4.3 Our diagnosis

The mode is the location-faithful readout of a one-sided posterior and should stay. The defect on a flat field
is the **variance of the modes across shallow objects**, and a posterior average moves every shallow object's
location (up or down) rather than reducing that variance. A repair must reduce the variance without moving the
location: shrink each object's mode toward a **local** reference by its own posterior width, where "local"
means its neighbours along the gene, never the library. That is the direction the owner has already approved
twice: the local capture footprint (a probe enriches adjacent regions and boundaries; 2026-10-07) and a proper
local enrichment prior, uniform in the retained background fraction `1/C`, density `1/C²` for `C ≥ 1`
(2026-10-09, an isolated prototype exists, not integrated).

An observation we cannot yet explain and want the reviewer's eye on: the same noise does not hurt the
**unstranded** capture-OFF rows (g00 unstranded OFF 1.59 %, against 13.61 % stranded). In an unstranded
library the strand columns carry no RNA-strand information and every RNA object is read through the share
integral; in a stranded one the reader resolves the share. Either the modes are noisier under a resolved
share (a narrower likelihood in `r` leaves the density's posterior wider on the low side?) or the EM's isoform
allocation is more sensitive when the strand is known. We have not measured which.

### 4.4 Lines of investigation we propose

1. **Characterise the noise before repairing it.** On g00 stranded OFF and g00 unstranded OFF, per object:
   the mode's deviation from the true level against the object's depth (reads), its posterior standard
   deviation, and its class (exon, intron, boundary). The question is whether the mode's variance tracks the
   posterior width (then a width-based shrinkage is derived, with the shrinkage scale the landscape's own
   spread and no new constant) and why the stranded rows differ.
2. **Dissect the damage by isoform structure.** Which transcripts gain silent mass: those differing from an
   expressed isoform by a single shallow object? By how much does their contracted-length ratio move? This says
   whether the repair belongs at the object (the weight) or at the component (the length rule's sensitivity to
   one object).
3. **The local shrinkage, derived.** Each object's mode pulled toward its neighbours' modes with weights from
   the posterior widths — a proper local prior, so that a deep enriched object keeps its weight, a shallow one
   next to it inherits the footprint, and a shallow one in a flat field inherits the flat level. The falsifiers:
   bit-identical where every object is deep; the single-probe case (one enriched gene in a flat field) must
   survive; the gap panels and both capture-ON strata read beside the OFF rows.
4. **Refused in advance:** a depth cut below which the weight is published as one (a binary cut on a
   continuous quantity, and a detector in disguise); any library-level reference level.

## 5. Issue 2 — the deferred stratum, unstranded × capture-ON

### 5.1 What happens

At the gDNA rungs, captured gDNA at exons is booked as nascent RNA and the gDNA fraction is under-called;
the zero-gDNA row is unchanged.

| g50 unstranded ON | 2026-09-30 | today |
|---|---|---|
| transcripts / genes | 10.7 % / 4.2 % | 65.2 % / 24.9 % |
| gDNA fraction (truth 0.50) | 0.442 | 0.343 |
| nascent estimate (truth 150 k) | 644 k | 1,052 k |

### 5.2 Attribution and the evidence of hope

This is the count foundation's accepted tradeoff, not the reader's. When the intron-background constraint
was removed (2026-10-09), the review of that change found that with the strand channel dead and the
constraint gone, every delivered composition row is a one-sided gene-edge floor, and that withdrawing those
edge rows "repairs the fixture but kills the deferred stratum's only enrichment witness on unstranded data".
The owner accepted the tradeoff knowingly, with the deferred stratum authorised for investigation "in depth
for robustness".

The stratum is not hopeless. Before that change, with a fragment-length gap between gDNA and RNA (the gap
panel, 2026-10-03 arm), unstranded × ON scored 8.2 % transcripts and 1.0 % genes with the gDNA fraction exact
(0.499). Today the same row is 34 % (RNA longer) and 68 % (RNA shorter). So two witnesses exist in principle —
the withdrawn edge rows and fragment length — and the ladder, by its equal-length design, can show neither;
the gap panels are the instrument for this stratum.

### 5.3 Lines of investigation we propose

1. **A replacement witness for enrichment on unstranded data,** in order of evidence: (a) the exon-to-intron
   boundary contrast — a measured finding of 2026-10-04 read about 700× on probed objects against 1.0 off
   capture from raw totals — a local, derived contrast that needs no strand; (b) the refused arm "a one-sided
   delivered row is a bound for admission", whose refusal predates the factory's removal and which the review
   itself said could be reopened after integration; (c) fragment length, which already carries the stratum
   wherever a gap exists.
2. **Dissect g50 unstranded ON by object class** (the per-class instrument) and read calibration against the
   oracle calibration there, to see which class books the gDNA as nascent and whether the landscape still
   forms its enriched mode on captured unstranded data.
3. **Measure on the gap panels first,** since the ladder cannot show a win on this stratum by construction.

## 6. Constraints that are fixed

- No capture detector, library-level capture status, or reference density: every weight is per object and
  continuous.
- No tuning to the simulated panels; constants must be derived and are refused otherwise.
- Robustness and generality across input distributions over aggregate simulated accuracy; conservative
  allocation under uncertainty (gDNA over RNA).
- Every score read per stratum; one mechanism per A/B; a numeric no-op proven before any speed change.
- The deferred stratum stays deferred for the release, and is investigated for robustness.

## 7. Questions for the reviewer

1. Is the diagnosis of Issue 1 (variance of the modes across shallow objects, not their location) sound, and
   is width-based local shrinkage the right family? What would you measure first?
2. Is there a formulation of the EM's effective length that carries the weights' uncertainty into the isoform
   allocation without moving each object's location — and without a constant?
3. The stranded/unstranded asymmetry at g00 (§4.3): what mechanism would you suspect?
4. For the deferred stratum: is a local exon-to-intron contrast a legitimate enrichment witness under the
   no-detector constraint, or does it smuggle a reference in? Which of the three witnesses in §5.3 would you
   pursue, and what would falsify it?
5. Anything in §1 that you believe is wrong in principle for cfRNA, independent of these two issues.

## Appendix A — instruments and data

- `scripts/design/quant_accuracy.py` scores a panel's arm; `scripts/sim/panel.py score --config
  scripts/sim/configs/<panel>.yaml --jobs 8 --arms base` runs it with the benchmark protocol (fractional
  assignment, the scan pinned). Panels: `gdna_ladder`, `flgap_rna_long`, `flgap_rna_short`, `test_reference`.
- `scripts/design/policy_benchmark.py --by-class` dissects calibration's error by object class;
  `scripts/design/calibration_vs_oracle.py` scores calibration against an oracle calibration.
- The arms: `$RIGEL_RUNS/suite/<panel>/arms/qa_<panel>_base.jsonl` (today's), with the preserved
  earlier arms beside them (`arms_2026-09-30_before_reader_and_multimapper/`,
  `arms_2026-10-03_before_reader_and_multimapper/`). The two refuted readouts' arms:
  `$RIGEL_RUNS/prototypes/2026-10-11_weight_expectation/{ab,ab_exp_density}/`. The attribution:
  `$RIGEL_RUNS/prototypes/2026-10-11_reader_speed/attrib/`.
- The ladder report page (rendered from today's arm): https://claude.ai/artifact/Cek2wmKtitgbDyfM5gNqyj.
- The permanent record: `DESIGN.md` §7.2 (the reader and its rulings), `EQUATIONS.md` §11 (the reader's
  likelihood and the conserved-share length), `ISSUES: the-capture-weights-are-noisy-on-a-capture-off-library`
  (Issue 1 with the attribution and the two refusals), `ISSUES: a-pure-gdna-library-reads-as-nascent-rna` and
  `DESIGN.md` §0b (the deferred stratum's rulings), `docs/TESTING.md` (the panels).
