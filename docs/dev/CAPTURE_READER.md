# The capture reader — false gDNA read as capture

Sandbox document for the session devoted to the reader. The owner called it a major cause of poor accuracy and
misprediction, said hybrid capture performs poorly, and gave it a dedicated effort (2026-10-28). It is taken now that
the second pass's draw is done (2026-10-29). The open record is `ISSUES: the-reader-reads-false-gdna-as-capture`; the
ruling is `DESIGN.md` §7.2, the derivation `EQUATIONS.md` §11. Nothing here is in `src/`. Every number is from
81fc67b6 unless dated.

## 1. What the reader is for, and what it does

Under hybrid capture a fragment is retained at a rate set by its overlap with the probes, whatever its origin. gDNA is
one template at a uniform rate before capture, so its density afterwards measures the capture yield at each object;
RNA's density is confounded with expression. So the tool reads the yield off gDNA (the ruler) and gives every EM
component — the locus gDNA component, every nascent span, every transcript — its conserved share of each object it
deposits on at that object's yield (the one shared rule). Off capture the right answer is 1 everywhere, and any other
value is error.

`capture_efficiency.capture_efficiencies`, per region (its contained gDNA count) and per boundary (its crossing count):

- **Prior**: the gDNA landscape the final sweep was solved under, read as Efron's two groups. The largest basin
  (`landscape.census`, `split_basins`) collapses to its peak `x_b`, the BULK, carrying its mass `π₀`; elsewhere the
  prior is the landscape's own density.
- **Likelihood**: the object's deconvolved gDNA mass `k = f_g·M` on its gDNA support `S`, as a Poisson count. `f_g` is
  ψ's posterior MEDIAN of the slot's λ-marginal (`psi_kernel.h`, `solve_slot`), and `M` the slot's unspliced count.
- **Readout**: `log c = E[min(x, s(x))] − x_b`, each enriched basin saturating at its own peak `s(x)`. An object with
  `S = 0` reads 1, and a one-basin landscape reads 1 everywhere.

Five consumers read the result (`pipeline.py`):

| consumer | what it does with the efficiencies |
|---|---|
| `capture_eff_length.transcript_capture_eff_lengths` | every transcript's and nascent span's EM length: Σ share·c over its pieces and cuts; a junction at `c_lo + c_hi − ½(c_intron,lo + c_intron,hi)`, its two exon\|intron boundaries less the intron pieces, clipped to `[0, max(S, every object it touches)]` |
| `priors.assemble_priors` | the locus gDNA component's length, `Σ S_r c_r + Σ M_e c_e` |
| `capture_eff_length.place_efficiencies` → `place_log_weights` | a multimapper's place weight `Σ m_o c_o`; an intergenic place reads 1 |
| `capture_eff_length.rna_captured_fraction` | the shared length law's captured fraction |
| `priors.intergenic_gdna_rate` | the bulk density itself, as the intergenic background's rate |

## 2. Today's measurement

`reader_census.py` runs one pinned scan, production's drain and calibration, then the EM twice from the replayed buffer,
fractional, scored by read name. The two arms are shipped and `c1`: every efficiency and the saturation at 1, with the
bulk kept. Transcript and gene error are Σ|Δ| against truth; `×` is shipped / c1, so below 1 the reader helps.

| condition | STAR tx × | STAR genes × | simulator tx × | simulator genes × | STAR saturation | STAR regions / boundaries >10× | STAR regions / boundaries <0.9× |
|---|---|---|---|---|---|---|---|
| g00 κ0.99 OFF | 9.44 | 14.95 | 1.00 | 1.00 | 149,367 | 229 / 581 | 0 / 0 |
| g05 κ0.99 OFF | 1.04 | 1.17 | 1.00 | 1.00 | 595 | 9 / 11 | 44 / 0 |
| g50 κ0.99 OFF | 1.00 | 1.00 | 0.98 | 1.00 | 32 | 4 / 6 | 67 / 1 |
| g98 κ0.99 OFF | 1.00 | 1.00 | 1.00 | 1.00 | 1 | 0 / 0 | 91 / 0 |
| g00 κ0.50 OFF | 1.53 | 3.21 | 1.00 | 1.00 | 1,743,403 | 64 / 75 | 0 / 0 |
| g05 κ0.50 OFF | 1.16 | 1.09 | 1.14 | 1.02 | 13 | 24 / 20 | 793 / 3,061 |
| g50 κ0.50 OFF | 1.08 | 1.08 | 1.00 | 1.00 | 11 | 4 / 5 | 67 / 0 |
| g98 κ0.50 OFF | 1.00 | 1.00 | 1.00 | 1.00 | 1 | 0 / 0 | 80 / 0 |
| g00 κ0.99 ON | 1.82 | 1.61 | 1.00 | 1.00 | 12,729 | 59 / 118 | 0 / 0 |
| g05 κ0.99 ON | 0.48 | 0.26 | 0.38 | 0.14 | 807 | 6,006 / 17,379 | 7,802 / 16,981 |
| g50 κ0.99 ON | 0.14 | 0.07 | 0.12 | 0.06 | 974 | 6,782 / 17,848 | 300 / 0 |
| g98 κ0.99 ON | 0.04 | 0.03 | 0.03 | 0.02 | 835 | 6,710 / 17,953 | 126 / 0 |
| g00 κ0.50 ON (deferred) | 1.04 | 1.07 | 1.00 | 1.00 | 12,472 | 22 / 28 | 0 / 0 |
| g05 κ0.50 ON (deferred) | 0.56 | 0.32 | 0.50 | 0.25 | 765 | 4,862 / 14,520 | 13,204 / 20,027 |
| g50 κ0.50 ON (deferred) | 0.17 | 0.12 | 0.15 | 0.10 | 755 | 5,962 / 15,309 | 8,224 / 10,851 |
| g98 κ0.50 ON (deferred) | 0.04 | 0.04 | 0.04 | 0.04 | 742 | 5,786 / 16,277 | 5,988 / 7,706 |

| stratum, Σ g05–g98 | STAR transcripts, shipped → c1 | STAR tx × | STAR genes × | simulator tx × | simulator genes × |
|---|---|---|---|---|---|
| stranded OFF | 397,299 → 390,030 | 1.02 | 1.06 | 0.99 | 1.00 |
| unstranded OFF | 475,804 → 425,464 | 1.12 | 1.08 | 1.08 | 1.01 |
| stranded ON | 754,768 → 5,033,090 | 0.15 | 0.07 | 0.12 | 0.05 |
| unstranded ON (deferred) | 1,240,323 → 8,201,754 | 0.15 | 0.09 | 0.14 | 0.07 |

What it says:

- **Under capture the reader is most of the tool's accuracy.** With every efficiency at 1, stranded capture-ON is 6.7
  times worse on transcripts and 14 times on genes.
- **Off capture, its cost now sits in the zero-gDNA controls and at weak strandedness.** On the STAR rung g00 κ0.99 is
  ×9.4 on transcripts and ×15 on genes, g00 κ0.50 ×1.5 and ×3.2, and the unstranded stratum ×1.12 and ×1.08. The
  stranded stratum is ×1.02 and ×1.06. On 2026-10-19, g50 κ0.99 capture-OFF read −34 % at efficiency 1; today it reads
  ×1.00. The changes since then (the scan's repairs, the observable genome, the multimapper places, the draw) removed
  it; which of them did is not separated.
- **The g00 capture-ON controls are hurt too** (κ0.99 ×1.8): with no gDNA there is no ruler.
- **On the simulator's alignments the enrichment side reads 1, but the depletion side does not.** g05 κ0.50 OFF costs
  ×1.14 on transcripts there, with 519 regions and 3,341 boundaries below 0.9× — the same as on STAR (§3, c′). Off
  capture, no object reads above 1 on the simulator's alignments in any condition.
- **Under capture at g05**, thousands of objects read below 0.9× (κ0.99: 7,802 regions and 16,981 boundaries). The
  reader still wins there by a wide margin; why they read below 1 is unexamined.
- **On a real library**, the VCaP transcriptome-only library (no capture, no spike, no truth) read 3,300 objects above
  e× the bulk and a saturation of 2.31e7× (2026-10-14). That is unphysical.

## 3. The mechanism, dissected

`reader_dump.py` dumps three STAR-rung conditions, which `reader_dissect.py` reads offline. "Located" means ψ's
`Var(log f_g) ≤ 1 nat²` (`CONSTANTS.landscape.located_var`).

| | g00 κ0.99 OFF (false enrichment) | g05 κ0.50 OFF (both sides) | g50 κ0.99 ON (true capture) |
|---|---|---|---|
| bulk, fragments per bp | 5.0e-7 = 1/2.0 Mb | 6.4e-3 = 1/155 bp | 1.8e-3 = 1/544 bp |
| landscape training: located + zero-count anchors | 120 + 5,186 | 7,959 + 1,383 | 13,828 + 3,044 |
| `π₀`; mass above the bulk; saturation | 0.994; 0.006; 149,367× | 0.989; 0.002; 13.2× | 0.662; 0.334; 974× |
| regions / boundaries above 10× | 229 / 581 | 24 / 20 | 6,782 / 17,848 |
| of them located | 64 / 101 | 16 / 20 | 6,461 / 17,460 |
| regions / boundaries below 0.9× | 0 / 0 | 793 / 3,061 | 300 / 0 |

**(a) At zero gDNA the bulk is a resolution wall, not a density.** The landscape trains on located compositions and
on zero-count anchors (`calibrate._fit_gdna_hyperprior`, replayed). At g00 the anchors' floor is the largest basin.
The bulk, 5.0e-7 per bp, sits 0.15 nat above the grid's bottom, 1/max(S): the resolution wall of the reference's
longest region, 2.3 Mb of support. So the saturation is a ratio against the reference's geometry, and every false-gDNA
object is read against it.

**(b) The false spectrum is a few dozen located false-gDNA kernels.** At g00, 39 of the 120 located kernels sit in four
basins 11–15 nat above the bulk, which hold 0.6 % of the mass: 5 intergenic regions (1,452 gDNA fragments), 28
single-strand exons and 6 single-strand introns. True capture at g50 κ0.99 ON is a different population: 5,448 located
exons with 1.9 M gDNA fragments in one basin holding 33 % of the mass, read against a bulk that is a density.

**(c) Most false capture is ψ's decision read as data.** At a slot both strands admit (AMBIG, `region_geometry`'s
G3), strand cannot tell gDNA from balanced sense and antisense RNA. The composition is then the prior's and the
messages'. At g00, 6,121 of 6,434 such slots resolve to RNA (`f_g < 0.05`); 68 resolve to `f_g` 0.6–0.99 and hold
15,750 of their class's 22,430 false gDNA fragments. Above 10×, AMBIG objects are 131 regions and 245 boundaries, with
a median Var(log f_g) of 14–17 nat², and the AMBIG introns alone carry 13,747 false gDNA fragments. The reader then
takes `f_g·M` as a sharp Poisson count, with √k precision, against the wall.

**(c′) The readings below 1 are empty objects against a floor basin.** At g05 κ0.50, the 793 regions and 3,061
boundaries below 0.9× hold no fragment at all. They are short (median support 91 and 215 starts) and expect 0.6 and 1.4
bulk fragments. The landscape keeps a basin at the grid's bottom, the zero-count anchors' resolution wall, 9.5 nat below
the bulk with 0.9 % of the mass. It takes 1.6–3.5 % of an empty object's posterior, and the log readout turns that into
0.7–0.9×. A junction priced from two such boundaries reads about 0.4, and 416 of g05 κ0.50's transcripts have EM
lengths under 0.9× their plain length. The simulator's alignments read the same (519 regions, 3,341 boundaries, 364
transcripts), so this half is not the aligner's. The 2026-10-19 picture — RNA-swamped exons, 11,055 below 0.9× at
g05 — is gone from today's tree.

The median is also discontinuous. Where a two-rooted posterior's mass crosses ½, the point jumps from one root to the
other. That is the untested candidate for why the reader amplifies every scan change.

**(d) Some false capture is located.** At g00, 64 of the 229 regions and 101 of the 581 boundaries above 10× are
located. They are single-strand exons and introns whose antisense reads ψ books as gDNA at κ0.99, and G1 objects
(intergenic regions and boundaries no RNA crosses), whose every unspliced fragment is gDNA by Axiom 0. On the
simulator's alignments nothing reads above 1 off capture, so these come from the aligner. On 2026-10-19 the
intergenic ones were STAR's unique placements of RNA in intergenic repeats, 231–435 fragments per region; where the
single-strand ones come from is not named. They are also the seeds of the false spectrum in (b).

**How it reaches the table.** A falsely captured exon|intron boundary raises the junction prices beside it, which
inflates the mature transcripts' EM lengths. A falsely captured intron subtracts from the junction price (½ c_intron)
and can clip it to 0. On 2026-10-19 the fragments moved from annotated transcripts to the less-inflated nascent spans
(annotated −386 k, nascent +366 k). Today at g00 κ0.99, 1,191 of 8,750 transcripts have EM lengths above 2× their plain
length, and the 99th percentile is 115,478×.

**Which objects carry the error.** `reader_ablate.py` shares one scan and calibration, then sets one population's
efficiencies to 1 per arm, leaving the saturation as read. Share of the shipped − c1 gap closed, transcripts / genes:

| arm | objects moved, regions / boundaries | g00 κ0.99 OFF | g05 κ0.50 OFF |
|---|---|---|---|
| only the readings above 10× | 229 / 581 | **99 % / 100 %** | |
| only the readings in (1, 10×] | 22,773 / 34,459 | 0 % / 0 % | |
| every reading above 1 | | | 76 % / 98 % |
| every reading below 1 | | | 23 % / −3 % |
| every boundary | 0 / 35,040 | 50 % / 30 % | |
| every region | 23,002 / 0 | 4 % / −2 % | |
| AMBIG objects | 3,760 / 4,671 | 33 % / 31 % | −2 % / −1 % |
| single-strand, not located | 10,372 / 18,922 | 40 % / 33 % | 5 % / 24 % |
| single-strand, located | 8,622 / 8,828 | 1 % / 17 % | 129 % / 262 % |
| G1 objects | 248 / 2,619 | 0 % / 1 % | 2 % / 11 % |
| the widened reader (§4) | 15,952 / 22,848 | 24 % / 26 % | |

The objects-moved column is g00's. At g05, every reading above 1 is 73 regions and 1,126 boundaries, and every reading
below 1 is 20,144 and 33,914. The gaps are, shipped → c1, transcripts 1,897,003 → 201,046 and genes 414,360 → 27,710
at g00, and 281,156 → 243,400 and 83,459 → 76,224 at g05.

- **At g00 the whole error is about 810 objects read at 10× to 4e6×.** The prior's small lift on tens of thousands
  of short objects — 1.01–1.1× on 11,000 regions with a median support of 7 — costs nothing, and neither do the G1
  objects' own readings.
- **Most of those 810 are unlocated** (AMBIG and single-strand), and half their transcript damage enters through
  the boundaries: the junction prices and contiguous cuts. 391 of the 581 boundaries above 10× are contiguous
  exon|exon boundaries, where an exon is split at another isoform's edge, so one isoform is inflated against its
  neighbour; 147 are exon|intron and 43 exon|intergenic.
- **At g05 κ0.50 on the STAR rung the enrichment side carries the cost** (76 % / 98 %), and the objects ψ calls
  located carry it. With strand dead at κ0.50, "located" there means the messages or the prior pinned the slot, not
  its own data. Setting them to 1 beats c1. The depletion side costs 23 % of the transcript gap and none of the
  genes'; on the simulator's alignments, which have no enriched basin, the same row's ×1.14 is all depletion.
- **The widened reader removes a quarter** of g00's gap end to end: transcripts 1,897,003 → 1,498,328 against c1's
  201,046.

## 4. The plan of 2026-10-19 does not hold

The 2026-10-19 plan read each object's gDNA "through ψ's Var(log f_g) instead of as a point count". That is the reader
refused on 2026-10-13 (`ISSUES: the-located-mode-capture-reader`, "the log readout widened by the
solve's Var(log f_g)"): the count's Poisson curve convolved in log ρ with N(0, Var(log f_g)), under the same prior, as
`tglh` and, with the shipped saturation, `tglhb`.

- **2026-10-13**, the test chromosome on the simulator's alignments, `tglb` (the shipped form) → `tglhb`:
  - the zero-gDNA capture-OFF rows gained on transcripts: flrl g00 κ0.99 7.36 → 4.99 %, flgl g00 κ0.70
    9.16 → 6.37 %, odg05 g00 κ0.70 10.83 → 7.94 %;
  - scarce gDNA under capture lost: lowg g001 κ0.99 ON 21.00 → 37.99 %, d10 g01 ON 30.32 → 39.30 %, with genes flat
    (0.27 → 0.12 %, 0.91 → 0.98 %). The loss is isoform allocation.
- **Today on the STAR rung**, offline: it leaves most of the false capture.
  - At g00 κ0.99 the objects above 10× go 229 → 166 regions and 581 → 413 boundaries. The AMBIG introns' false gDNA
    above 10× goes 13,747 → 13,295 fragments.
  - At g05 κ0.50, below 0.9× goes 793 → 792 regions and 3,061 → 3,061 boundaries.
  - At g50 κ0.99 ON nothing moves.
  - End to end at g00 it closes 24 % of the transcript gap and 26 % of the gene gap (§3's ablation).

Why: the convolution is centred on ψ's median. A decision 11–15 nat from the bulk stays several standard deviations
away even at 15 nat². The readings below 1 are empty objects, whose zero count is exact, so there is nothing to widen.
Var(log f_g) says how unsure ψ was; it does not undo where ψ's prior put the median. The ISSUES entry and ROADMAP step 7
record the withdrawal (2026-10-29).

## 5. The crux

The target is narrow. At g00 the whole error is about 810 objects read above 10×, most of them unlocated. At g05 κ0.50
it is a few hundred objects whose "location" is the prior's or the messages'. But per object, the false capture at g00
and the true capture at scarce gDNA have the same evidence: an object whose own data cannot locate its gDNA, under a
population whose enriched spectrum is thin.

- The shipped reader reads ψ's prior-and-message decision as data. That keeps sparse true capture — the widening is
  the only difference between `tglb` and `tglhb`, so all of lowg's advantage lives in objects ψ located loosely — and
  it is also what makes the false capture.
- A reader of the population prior alone erases both.

What separates the two cases is not in the object's own count. It is whether the spectrum's located kernels are true
(5,448 exons and 1.9 M fragments) or false (39 kernels, 0.6 %), and whether the object's own neighbourhood carries
capture. Two smaller questions sit beside it, both about walls read as levels. At zero gDNA the bulk is the deepest
wall, so a located false object reads 1e5–1e6×. And the floor basin is made of walls, so an empty object reads
0.7–0.9× (the 2026-10-13 review's lesson: a zero-count kernel's centre is a resolution wall, not a location).

## 6. What a fix must respect

- **Capture is a spectrum** (owner, 2026-10-05): no detector, library-level test or on/off status. Every weight is per
  object or per placement and continuous; "no enrichment" is the limit where the weights are 1.
- **Capture has a local footprint** (owner, 2026-10-07). Enrichment at an object is reason to expect related enrichment
  at adjacent regions and boundaries, junctions included, at the existing object resolution. The coupling is
  approximate and local: not equality of averages, no fixed flank ordering, and no requirement that neighbours show
  capture first. "No particular coupling strength or new prior has yet been selected."
- **Capture is local** for a junction: priced from the objects within one fragment of it. **Junctions are never
  pooled.** **Capture stays outside the EM.** **The one shared rule** for every component (`DESIGN.md` §7.2, rulings
  1–5).
- **An intergenic place is uncaptured** (owner, 2026-10-27). Today it is applied to the multimapper weights and the
  intergenic rate only.
- **No constant, floor or threshold without discussion**, and no magic number.
- **gDNA called RNA is the priority bleed**; RNA called gDNA is accepted for now (owner, 2026-10-27/28).
- **Real data is a test input**: the VCaP transcriptome-only library is a test, never a design input.

## 7. Refused — do not rebuild

The full records are in `ISSUES.md`; read them before re-proposing a reader.

- `ISSUES: the-located-mode-capture-reader`: a mode located and switched off when absent, plus ten readers of the same
  two-group prior. Among them are the linear mean (test κ0.99 g00 OFF 6.99 → 32.61 %), the null ratio with Var
  widening, the log readout with Var widening (§4), one saturation level, NPMLE, typed base rates or weights, the nested
  null and the smoothed-bulk null.
- `ISSUES: the-one-sided-null`: it erases a well-resolved ratio where capture is the majority.
- `ISSUES: the-likelihood-fitted-null-location`: even the oracle's density deepens the sparse losses.
- `ISSUES: the-unclipped-capture-efficiency`.
- `ISSUES: calibration-detects-capture-on-a-capture-off-library`'s refusals: the `eff ≥ 1` training guard, the
  evidence-fitted reference, the posterior-median and quasi-score rulers.
- `ISSUES: an-unconditional-class-conditional-landscape`.
- The junction family: `ISSUES: pooling-junctions`, `ISSUES: the-spliced-read-junction-price`,
  `ISSUES: the-junction-price-at-its-neighbours-mean`, `ISSUES: a-junction-price-read-from-a-fitted-capture-map`.

## 8. Directions, each derived on paper before any code

1. **Read ψ's evidence once.** The reader's likelihood becomes ψ's own evidence for the object's gDNA density: the
   slot's λ-posterior with ψ's landscape term divided out (the strand term and the messages), mapped to
   `x = log f + log M − log E`, with the two-group prior applied once. That removes the double prior EQUATIONS §11
   admits.
   - **Predicts**: an object whose data cannot locate its gDNA reads the population prior, which is about 1 at g00 and
     g05 off capture, where `π₀ ≥ 0.99`; a located object is unchanged, as at g50 ON.
   - **Needs** ψ's λ-posterior per slot: a Python prototype can export it on the ladder (70 k slots × the lattice); the
     shipped form would be a native read-out of the final sweep (a C++ worktree).
   - **Open in the derivation**: whether the messages belong in the likelihood (they are neighbours' evidence, and a
     gDNA level crosses a capture boundary unscaled); the reference term; and how the doubled ½ log f of
     `ISSUES: the-gdna-prior-enters-psi-twice` divides out.
   - **Risk**: it can repeat the refusal at scarce gDNA under capture, so measure lowg and d10 first.
2. **Local co-enrichment for an object that cannot read its own gDNA** — the 2026-10-07 ruling's open coupling. An
   unlocated object's prior would come from its neighbours' located readings (its boundaries and adjacent regions,
   within a fragment) instead of the library's whole population.
   - At g00 the AMBIG introns' single-strand neighbours read no gDNA. At scarce gDNA a probed exon's boundaries carry
     its probe.
   - This is the candidate for keeping sparse true capture under direction 1. The coupling must be derived, not chosen,
     and a junction is still priced from its own neighbours, never pooled.
3. **The located false gDNA, upstream.** Name g00's located false objects on the STAR rung by read name: intergenic
   repeats, antisense placements at single-strand exons, or something else. Their repair belongs in the scan or the
   multimapper handling, in pipeline order. They also seed the false spectrum, so directions 1 and 2 are measured with
   and without them, using an oracle arm that drops them from the landscape's training.
4. **ψ's doubled half.** Part 1 of `ISSUES: the-gdna-prior-enters-psi-twice` removed the ladder's zero-gDNA false gDNA
   on the simulator's alignments (−97 %) and was never run on the STAR rung. One arm: does it move g00's AMBIG false
   gDNA? It was not shipped for its capture-ON cost, which is itself read through this reader.
5. **Walls are not levels.** The grid's bottom is 1/max(S), and the zero-count anchors' kernels pile there. A basin
   made of walls is the floor the empty objects read against, and at zero gDNA it is the bulk itself. Derive what the
   two-group prior should do with a basin whose members locate nothing. The one-sided null is refused, so the answer
   is not "everything below the bulk is the null".

## 9. The session's order

1. Read this document, `ISSUES: the-reader-reads-false-gdna-as-capture` and the entries it cites, `DESIGN.md`
   §7.1–7.2 and `EQUATIONS.md` §11. In the source: `capture_efficiency.py`, `landscape.census` / `split_basins`,
   `calibrate._fit_gdna_hyperprior` and `psi_kernel.h`'s `solve_slot`.
2. Dissect, with nothing in `src/`:
   - the ablation on the other rows the reader costs on the STAR rung (g00 κ0.50 OFF, g50 κ0.50 OFF, the g00
     capture-ON controls) and on one capture-ON row;
   - at g05 κ0.50, which objects ψ called located and what pinned them;
   - g00's located false objects, named by read name;
   - shipped against widened at lowg g001 ON and d10 g01 ON: which objects carry the difference, and what evidence ψ
     has at each (strand counts, Var, messages, neighbours);
   - the median-flip hypothesis: between two trees one small scan change apart, count the AMBIG slots whose median
     changes root.
3. Derive the reader on paper, where EQUATIONS §11's reader paragraph sits, with hand-built falsification tests:
   - a balanced-strand AMBIG slot under a landscape with a false enriched basin reads 1;
   - a located captured exon reads its level;
   - an empty short object reads 1 against a landscape whose floor basin is made of walls;
   - an object ψ calls located only through its prior or messages reads no more than its own data allow;
   - a sparse captured exon beside located captured boundaries keeps its capture.

   Break each gate on purpose and watch it fire.
4. Prototype outside `src/`:
   - `ruler_vs_truth.py --module` for the lengths with no EM, reading the probed and unprobed classes apart;
   - `reader_census.py` extended with the arm, end to end;
   - a C++ worktree for anything inside ψ.
5. A/B per stratum, pinned and fractional, by read name, with c1 beside, on:
   - the STAR rung and the simulator's alignments (16 + 16);
   - the scarce-gDNA capture-ON rows (lowg, d10, d100), the refusal's own panels;
   - the capture panels (single, sparse, junction, the binding panels) and the fl-gap panels;
   - then VCaP transcriptome-only, one real-genome job at a time.
6. **The bar**:
   - every capture-OFF STAR condition near its c1 row, the g00 controls included;
   - capture-ON no worse than shipped on the ladder and the capture panels;
   - lowg and d10 capture-ON no worse than shipped;
   - gDNA called RNA not up;
   - VCaP transcriptome-only's objects above e× far below 3,300, with its saturation near 1.

   Bring the per-stratum table to the owner before anything enters `src/`.

## 10. Instruments and receipts

- `$RIGEL_RUNS/prototypes/2026-10-29_reader_census/`:
  - `reader_census.py` (shipped against c1, one condition), `run_all.sh` (32 rows, four at a time), `report.py`
    (`report.md`), `runs/`;
  - `reader_dump.py` and `reader_dissect.py`, whose dumps and readings are in `dump/`;
  - `reader_ablate.py`, whose output is in `ablate/`.
- `$RIGEL_RUNS/prototypes/2026-10-13_modefree/`: every 2026-10-13 reader offline (`readers.py`; `honest_curves` is the
  widening), the arms' receipts per panel (`out/mf_*`), the refusal tally (`review/refusals.txt`) and the real
  libraries (`prod_real.py`, `out/real/`).
- `$RIGEL_RUNS/arms/2026-10-19_annotated_path/dissect/`: `calib_dump.py` and `run_cond.py`. ⛔ They read
  `star_panel/<condition>/recommended.bam`, the old STAR panel; the rung is now
  `suite/ladder/<condition>/rustar_recommended.bam`.
- `scripts/design/`:
  - `ruler_vs_truth.py` (`--module`: a file of `ARMS`, each handed calibrate's landscape, chain and belief);
  - `quant_accuracy.py`, which has no arm with every efficiency at 1, so the census is that arm;
  - `calibration_vs_oracle.py`, which cannot see the efficiencies.
- The panels:
  - the ladder: the STAR rung and the simulator's alignments;
  - the test chromosome's depth ladder (`scenarios_depth_d10`, `_d100`, `_full_lowg`): simulated before the 2026-10-21
    re-simulation, on the simulator's alignments only, with scan caches from 2026-09-16, so read their BAMs directly;
  - the capture panels (`scenarios_probes_*` and the binding panels);
  - the fl-gap panels.
- The A/B baseline is 81fc67b6. The draw's kept `t_counts` are in `$RIGEL_RUNS/prototypes/2026-10-29_draw/runs/`.

## 11. Starting the session

> We are starting the capture reader for Rigel 0.8.0. Read `docs/dev/CAPTURE_READER.md` first, then
> `ISSUES: the-reader-reads-false-gdna-as-capture` and the refusals it cites, `DESIGN.md` §7.1–7.2 and `EQUATIONS.md`
> §11. Under capture the reader is most of the tool's accuracy. Off capture on the STAR rung it costs the g00 controls
> ×9.4 on transcripts and the unstranded stratum ×1.12. At g00 the whole cost is about 810 objects read above 10×, most
> of them ψ's prior-driven decisions taken as data, against a bulk that is a resolution wall. The 2026-10-19 plan
> (widen the count by ψ's Var(log f_g)) is the reader refused on 2026-10-13 and is withdrawn: it erases sparse true
> capture and closes only a quarter of g00's gap.
>
> Start with the dissection in §9 step 2, with nothing in `src/`: the ablation on the other costly rows, g00's located
> false gDNA named by read name, and which objects carry the shipped reader's advantage at scarce gDNA under capture.
> Then derive a reader on paper. Direction 1 (ψ's evidence read once) is the first candidate; direction 2 (local
> co-enrichment) is for whatever direction 1 erases. Prototype it outside the tree, and bring me the per-stratum table,
> with c1 beside it, before anything enters `src/`.
