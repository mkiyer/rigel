# External review: the capture reader, measured

*2026-10-10, revising the 2026-10-09 version after the owner's four objections: rule on benchmarks, not on a
ruling; the gDNA background is not a trustworthy reference under efficient capture; prove the transfer
limitation across inputs and probe placements; deletions only for code truly unneeded, and everything
committed and pushed first. Sandbox review; nothing here is authoritative. Nothing in `src/`, `tests/`,
`tests/golden/` or the installed native changed; nothing was committed, pushed or deleted. The prototype
reader and every benchmark script live in a session scratch directory (Appendix B) and read only cached
panels and the real BAMs. The machine was shared with a foreign workload for most of the run (load 20–170 on
16 cores), so no wall time below is a speed claim. Every run named below finished; nothing is pending.*

## Verdict, revised

**The owner's premise is confirmed by the data and decides the frame.** Anchoring capture on the abundant
population, not on the depleted background, is what works. The first version's background-relative readout
(`slab_bg`) is the worst arm on every captured stratum, and its failure is the count-dependent shrinkage that
version predicted for itself. The same honest observation likelihood read under the gDNA landscape the count
solver already fits (`land_med`, the posterior median of density: no background, no reference density, no
located-mode test, no clip, no `None`) matches the shipped reader where the shipped reader works and does two
things the shipped reader cannot:

- it reads a single probe, the owner's cfRNA case: the shipped reader forms no located mode and reads the one
  100× enriched gene as uncaptured on every row (partial class −4.59 nats); `land_med` reads it within ±0.1 on
  the stranded rows with every other object left at 1; end to end the gene error at g50 falls from 3.9 to 1.7
  and from 2.6 to 1.4 %;
- it survives capture strength 100 per base, where the off-target pool holds 25–75 reads in the whole 5.9 Mb
  chromosome: probed transcripts within 0.1 on 99 %, where the background-relative arm drops to 79–88 %.

**Its in-scope cost is real, measured, and has one mechanism, now located.** End to end on the test chromosome
(transcripts % / genes %, means over the stratum's non-zero rows): stranded × ON 22.0 / 10.1 → 24.4 / 10.6;
stranded × OFF 14.9 / 4.0 → 16.3 / 4.4; unstranded × OFF 16.2 / 4.6 → 17.2 / 4.8. The loss sits at g05–g25
(transcripts +2…+10 points on captured rows, +0.3…+5.5 on capture-OFF rows, genes +0.0…+1.2); g50 is a small
win on most rows; g98 is a wash. The mechanism is not the likelihood and not the landscape's location: at
every object with data the reader is right (slots with ≥ 5 true gDNA reads: −0.02…+0.10 nats), and the fitted
landscape's depleted mode sits within 0.1 nat of the true background on every normal-depth row. It is the
posterior **median** at objects with less than one expected read: their posterior is the prior's broad,
left-skewed depleted part, whose median lies below its mode by 0.15 (g50) to 0.36 (g05) nats, 0.39 under
strength-100 capture and 0.79 at a tenth of the depth. The shipped reader inherits the same breadth with the
opposite sign (its clipped posterior mean over-reads the same objects, +0.12…+0.37). The cheapest repair
needs no background: read the posterior **mode**, which at a data-poor object is the population's typical
density and at a data-rich one its own. It passes its gate at both levels (§1.4): at
calibration level data-poor objects move from −0.14…−0.39 to within ±0.06 of the true background on every
normal-depth row, captured objects stay where the median put them, the g05 unprobed class recovers from
−0.36 / −0.21 to −0.07 / +0.10 (the shipped reader: +0.12 / +0.37) and the probed class returns to the
shipped reader's 0.99; end to end the four captured g05–g25 stranded rows go from 6.27 / 0.79 (shipped)
and 10.17 / 1.25 (median) to 6.94 / 0.87, the single probe's genes from 2.12 to 1.22 with transcripts
9.22 → 9.81, strength 100 to within 0.4 / 0.3 of the shipped reader, and the junction-probed rows with the
touched cap from 33.0 / 1.6 to 25.0 / 2.2. What remains is the capture-OFF transcript cost, now
+0.4…+3.1 points on the g05–g25 rows with genes equal (6.88 / 1.00 → 8.73 / 1.01 over four stranded rows):
the price of reading each object's own density where the shipped reader's `None` reads exactly 1.

**Two rule defects were found, one fixed.** Capping the junction at the larger of its two exon pieces, which
the first version proposed, is wrong wherever probes do not tile exon interiors: on the junction-probed
layout it reads every fully probed transcript as uncaptured (partial class +5.2 nats). Capping at the largest
weight among every object a junction fragment touches (the simulator's own law: a fragment is captured by its
best single probe part) repairs it (−0.2 / +0.3, better than the shipped reader's +1.2 / +3.2; end to end
over the layout's eight stranded captured rows 39.2 / 7.5 → 30.2 / 8.1). On the
sparse-probe layout no honest arm matches the shipped reader: a partly covered exon's density sits between
the landscape's modes and the median snaps to the lower one (−0.6 nats at the object); the arithmetic mean
interpolates there but manufactures false capture on zero-DNA objects (10 % above 90× on the same rows) and
is rejected. Open.

**The deferred stratum stays deferred for a structural reason, now measured:** on unstranded captured data
every honest reader is far worse than the shipped one (30 / 18 → 45 / 21 %). Local evidence cannot see
capture where the strand channel is dead; the count solver's population feedback can, and the shipped reader
borrows it by reading the inferred count.

**Recommendation.** Do not land `land_med`: the in-scope transcript table gets worse, and that table is
what the release ships on. Land the mode readout (`land_map`) with the touched-objects junction cap, under
the gates of §7, because on the test chromosome's in-scope strata it is within a point or two of the shipped
reader on transcripts and equal on genes while removing the detector, reading a single probe, surviving
strength-100 capture and contracting the plasma library the shipped reader leaves uncontracted (LBX0190:
half a percent of objects read enriched, gDNA fraction 8.5 → 10.0 %, measured with the median readout); and
declare its one measured cost, the capture-OFF transcript points at low gDNA, in the release notes. One
reservation stands against that recommendation, measured where it counts: at genome scale, and most on the
RNA-short gap panel (gDNA fragments longer than RNA's, the plasma regime), both honest readouts trail the
shipped census on the probed class at calibration level (0.64–0.67 within 0.1 against 0.96 on that panel's
stranded captured row; 0.70–0.83 against 0.96–0.98 on the ladder's), and the deficit is a lower tail: a tenth
of the probed transcripts read 15–30 % too short (§1.3, §1.4). The per-object scoring shows the honest readouts right at every object with gDNA reads and wrong at
exons shorter than a gDNA fragment, which carry no DNA evidence of their own and receive the neighbour-imputed
DNA level 5 % of the time; but the cheap consumer rule that prices such pieces at their boundaries was run and
does not close the tail, so the tail's objects are not yet named. On the ladder end to end, the release's own
metric, that spread costs +0.3 transcript points on the captured rows with genes better by 0.8, and +1.1 on
the capture-OFF rows with genes equal (§1.4): no row moves against the shipped reader by more than the test
chromosome's cost, so the recommendation stands, with the tail's dissection as the first item after landing. The alternative is
(b): ship 0.8.0 with the current reader and its two known failures declared (LBX0190 and MO_3021 quantified
uncontracted; a single probe read as uncaptured), keeping the honest reader as the line after 0.8.0. The
owner chooses between a declared cost on synthetic low-gDNA capture-OFF rows and a declared blindness on
real plasma libraries; the numbers for both are in §1. The third failure on the first
version's list, the RNA-short false reference on the genome-scale panel, is not reproduced by the current
tree at calibration level: on that row the shipped reader now forms no reference and reads every weight as
exactly 1 (§1.3); its end-to-end number was not re-run here. The numerical contract of the first version stands (one fixed lattice, no
adaptive integrator), with the two lattice rules its own step-halving check caught (§6). The background
estimator stays deleted; the `1/C²` slab and the arithmetic mean are refuted and not to be revived.

## 1. Which reader — the data (owner's point 1)

### 1.1 What was built and how it was read

One evaluator, one consumer, four readouts of the same likelihood, run as prototype arms of the existing
truth instrument and of an in-process pipeline harness; nothing in `src/` changes.

- **Evaluator.** For every region and boundary, the observation likelihood over gDNA density of the brief's
  model: the two Poisson strand columns with the RNA amount integrated under its `r^(−1/2)` reference, the
  delivered composition row at `log(ρ·Eg/r)`, the held DNA level at `ρ`, the held RNA levels at their own
  amounts, the two pure-strand atoms and the arcsine continuum under the witness bits. Fixed lattices only
  (§6). The typed factors are the production message tables, rebuilt from the final refit's sweep inputs by
  the transfer gates' own harness; the rebuilt delivered rows equal production's bit for bit on every one of
  the 260 conditions run (`rows_max_diff = 0`).
- **Readouts.** `slab_bg`: the brief's Gamma intergenic background (Jeffreys location; moments shape combined
  with the mean's posterior width, as the deleted estimator did), the `1/C²` slab at equal odds, the counterarm
  log score clipped at 1; absolute, background-relative. `land_med`, `land_geo`, `land_mean`: the posterior
  median, geometric mean and arithmetic mean of density under the fitted gDNA landscape, relative weights (the
  EM is invariant to a common factor: `test_the_split_is_invariant_to_a_common_thinning_of_every_yield`).
  `land_map` (§1.4): the posterior mode.
- **Consumer.** The shipped conservation sum over pieces and cuts, with the junction capped either at the
  larger exon piece (`pieces`, the first version's proposal) or at the largest weight among the two pieces and
  the two boundaries a junction fragment touches (`cap4`).
- **Controls.** `shipped` (production) and `oracle_gdna` (production fed the true gDNA counts).
- **Instruments.** Calibration level: the scoring of `ruler_vs_truth.py` (each transcript's effective length
  against the simulator's own capture-aware yield, anchored on the probed class, transcripts with ≥ 20 RNA
  fragments: the probed class's share within 0.1 nat, the partial and unprobed classes' median log error, the
  within-gene sd). Per object: each weight against the simulator's gDNA yield at that very object, so the
  reader's error is read apart from the DNA-to-RNA transfer; and in absolute units against the library's true
  off-target density, so a bias at data-poor objects is read apart from a bias at data-rich ones. End to end:
  one scan and one calibration per condition, every arm's EM on the same replayed buffer, scan pinned,
  fractional assignment, seed 0; the control arm reproduces the published gene numbers exactly and the
  transcript numbers within the EM's known last-bit sensitivity (g50 ss0.99 ON: 10.33 / 1.20 against
  10.17 / 1.20). Real libraries: the same harness on the cfRNA BAMs, pools and runtime only.
- **Panels.** All cached, nothing re-simulated: the test chromosome's benign layout (30 conditions), the
  junction-probed and sparse-probed twins (ON rows), three length-gap arms (RNA shorter, RNA longer, equal at
  200 bp; 30 each), two depths (1/10 and 1/100 of the benign 1.17 M fragments, with 0.1 % and 1 % gDNA rows
  added), three capture strengths (0.1, 1 and 100 per base against the benign 10), the single-probe panel,
  three zero-control seeds; the suite's two genome-scale gap panels and the ladder's in-scope rows at
  calibration level, and the ladder's six stranded rows end to end. Every number is read per stratum, never pooled.

### 1.2 Benign layout, equal lengths: where the shipped reader is right, so is `land_med`; where it is not, neither is

Calibration level, stranded × ON, g25–g98 (6 rows): probed within 0.1, `shipped` 0.99–1.00, `land_med`
0.97–1.00; partial-class median error `shipped` +0.00…+0.07, `land_med` −0.01…+0.01; within-gene sd 0.003
against 0.003–0.005. At g05 `land_med` trails: probed within 0.94–0.96, unprobed class −0.21 / −0.36 against
the shipped +0.37 / +0.12 (both wrong, opposite signs: §1.4). `slab_bg` on the same rows: probed within
0.30–0.97, partial −0.04…−0.55, within-gene sd 0.004–0.033. `land_geo` reads 10 % of zero-DNA objects above
7× on captured rows (q90 +2.0): a bimodal posterior's geometric mean lands between the modes.

End to end, per condition (transcripts % / genes %; `shipped` → `land_med`):

| stranded × ON | shipped | land_med | | stranded × OFF | shipped | land_med |
|---|---|---|---|---|---|---|
| g05 ss0.70 | 6.11 / 0.77 | 16.37 / 1.94 | | g05 ss0.70 | 6.24 / 0.98 | 8.49 / 1.07 |
| g05 ss0.99 | 4.95 / 0.35 | 7.39 / 0.40 | | g05 ss0.99 | 6.42 / 0.89 | 11.90 / 0.90 |
| g25 ss0.70 | 6.73 / 1.26 | 10.27 / 1.85 | | g25 ss0.70 | 7.57 / 1.15 | 8.80 / 1.26 |
| g25 ss0.99 | 6.93 / 0.77 | 6.71 / 0.81 | | g25 ss0.99 | 7.27 / 0.98 | 8.68 / 0.99 |
| g50 ss0.70 | 9.88 / 2.24 | 9.68 / 2.33 | | g50 ss0.70 | 10.84 / 1.40 | 11.15 / 1.62 |
| g50 ss0.99 | 10.33 / 1.20 | 9.38 / 1.36 | | g50 ss0.99 | 9.28 / 1.16 | 8.92 / 1.23 |
| g98 ss0.70 | 71.94 / 45.54 | 76.03 / 47.21 | | g98 ss0.70 | 36.59 / 13.47 | 37.25 / 15.00 |
| g98 ss0.99 | 58.90 / 28.26 | 59.35 / 28.94 | | g98 ss0.99 | 34.63 / 11.62 | 35.17 / 12.94 |
| mean (8) | 21.97 / 10.05 | 24.40 / 10.60 | | mean (8) | 14.86 / 3.96 | 16.29 / 4.38 |
| zero controls (2) | 7.95 / 1.42 | 7.87 / 1.40 | | zero controls (2) | 6.67 / 0.87 | 6.76 / 0.87 |

Unstranded × OFF (in scope): g05 7.91 / 1.09 → 9.13 / 1.23; g25 8.07 / 1.16 → 9.67 / 1.27; g50 10.50 / 1.54 →
11.28 / 1.66; g98 38.13 / 14.73 → 38.71 / 15.13; mean 16.15 / 4.63 → 17.20 / 4.82. Unstranded × ON (deferred):
30.38 / 17.64 → 45.21 / 20.70 (`land_geo` 56 / 30, `slab_bg` 80 / 51). The g98 rows are where every arm, the
oracle included, fails on 2 % RNA, and they dominate every stratum mean; the reader question lives at g05–g50.

Reading: on the capture-OFF rows the shipped reader's `None` gives every object exactly 1, and `land_med`
gives data-rich objects their true density but data-poor ones 0.15–0.18 nats less (§1.4), a within-gene
contrast between a gene's long and short exons that re-splits near-tie isoforms; genes do not move. On the
captured rows the same contrast adds to the transfer error of the unprobed class. The mode readout removes
the contrast (§1.4) and leaves the per-object spread of any honest reading, which costs +0.4…+3.1 transcript
points on the g05–g25 capture-OFF rows with genes equal.

### 1.3 Robustness axes

| axis (panel) | shipped | land_med | reading |
|---|---|---|---|
| capture strength 0.1 / 1 / 10 / 100 per base; stranded rows | probed within 0.1: 0.91–1.00; partial −0.00…+0.03 | 0.86–0.99; partial −0.08…+0.00 | tracks the shipped reader across a 1000× range; end to end 10.30 / 1.25 → 11.13 / 1.59 (0.1), 10.61 / 1.83 → 10.73 / 1.95 (1), 10.60 / 3.31 → 10.81 / 3.67 (100) |
| capture strength 100: the background | reference located from thousands of kernels | off-target pool 25–75 reads in 5.9 Mb; probed within 0.99 | the landscape frame needs no background; but data-poor objects read −0.39 (§1.4) |
| single probe, 6 conditions | partial class −4.59 on every row (no located mode: every weight 1) | −0.10…+0.02 on stranded rows; −2.1 unstranded; zero-DNA objects q90 ≤ 0.04, max ≤ 0.54 | locality and no-detector demonstrated; end to end genes 3.87 → 1.68 and 2.55 → 1.42 at g50, transcripts +2.6…+6.0 at g05 (the data-poor mechanism) |
| zero controls, g00 × 3 seeds, OFF and ON | 1 everywhere | q90 +0.01, max ≤ 1.2 on captured-library zero rows | no false contraction; end to end within ±0.6 transcripts, genes identical |
| depth 1/10 and 1/100 (12 conditions each, stranded, with 0.1 % and 1 % gDNA rows) | 1/10: probed within 0.92 / 0.99 / 0.99 at g05 / g25 / g50; 1/100: 0.08 / 0.87 / 0.89 | 1/10: 0.51 / 0.83 / 0.96; 1/100: 0.18 / 0.54 / 0.39 | the honest reader degrades faster with sparsity: at 1/10 depth and g05 1,999 of 2,103 uncaptured objects have no read and read −0.79 (§1.4); the shipped census pools the whole library. End to end (12 rows each, transcripts / genes): 1/10 depth 17.2 / 1.2 → 19.9 / 1.3; 1/100 depth 30.4 / 2.5 → 32.0 / 2.6; the loss sits on the captured g01–g05 rows (1/10: g01 26.9 → 38.4, g05 15.0 → 21.1) |
| length gaps, RNA shorter / RNA longer / equal | §3 | §3 | the transfer error is of the same size in both readers |
| probe layouts, junction-probed and sparse | §3 | §3 | the cap defect, fixed; the sparse-layout snap, open |
| genome-scale suite (70,176 slots per condition): RNA-short gap panel, stranded × OFF g50 | every weight exactly 1 (within-gene sd 0.000): no reference forms on the current tree, so the 42 % false-reference failure recorded on 2026-10-03 is not reproduced at calibration level | unprobed class −0.01, within-gene sd 0.028, zero-DNA q90 0.00 | both read the row clean |
| genome-scale suite: RNA-short gap panel, stranded × ON g50 | probed within 0.96; partial +0.24; unprobed +2.98 (the oracle fed the true counts: 0.84 / +0.20 / +1.50) | median 0.67 / +0.20 / −0.41; **mode 0.64 / +0.21 / +1.76** | **the one place the honest readouts trail the shipped census in scope at calibration level**, and the per-object scoring locates it (§1.4): both readouts are right at every object with gDNA reads or with ordinary opportunity and none (no-read uncaptured objects −0.02, captured objects with ≥ 5 reads −0.00 under the mode), and wrong at the 3,342 captured exons that are shorter than a gDNA fragment: a fraction of a base of gDNA opportunity, 14 bases of RNA opportunity, four reads, all RNA. There the reading is the prior (median +0.54, mode +0.62, 10th–90th percentiles −6…+1.8), because the object's own DNA evidence is nil by construction and the message layer delivers a DNA level to 5 % of them (the composition row is present at 81 % but cannot place a density when the DNA opportunity is a fraction of a base). A coin-flip weight at a few percent of a transcript's share moves its length by up to 30 %, which is where the probed class goes. The shipped reader's clipped posterior mean cannot exceed its reference and lands between the modes there. The unprobed class reads nearer the oracle under both honest readouts (+1.76 against the shipped +2.98): the shipped reader's unprobed class is 20× too long here and the oracle's 4.5×, a transfer error of the junction rule at this length gap. The ladder's rows follow below |
| genome-scale suite: RNA-long gap panel, stranded × ON g50 | 0.98 / +0.00 / +1.02 (oracle 0.95 / −0.15 / +0.69) | median 0.82 / −0.19 / +0.61; mode with the touched cap 0.87 / −0.17 / +0.63 | the same pattern: the honest readouts trail on the probed class and are nearer the oracle on the unprobed one; the capture-OFF twin reads clean for both (every shipped weight 1; median sd 0.031) |
| **the ladder** (genome scale, equal lengths), stranded × ON g05 / g50 / g98 | probed within 0.96 / 0.98 / 0.97; partial +0.08 / +0.08 / +0.07; unprobed +1.24 / +0.24 / —; within-gene sd 0.040 / 0.035 / 0.022 (the oracle: 0.97 / 0.98 / 0.97) | mode with the touched cap: 0.70 / 0.75 / 0.83; partial +0.04 / +0.01 / +0.01; unprobed +0.98 / +0.08 / —; sd 0.065 / 0.042 / 0.032 (the median: 0.59 / 0.73 / 0.70) | at genome scale the honest readouts' probed class is spread 1.5–2× wider than the shipped census's on every captured row, at equal fragment lengths too; the partial class is nearer zero. The capture-OFF rows read clean for every arm (mode sd 0.015–0.022, zero-DNA q90 ≤ 0.02). End to end these rows are the release metric: §1.4 |
| genome-scale suite: RNA-long gap panel, unstranded × ON g50 (deferred) | probed within 0.10, within-gene sd 0.41 | 0.80 / +0.07 / +7.02, sd 0.034 | both fail the deferred stratum at genome scale in different places |
| real plasma library LBX0190 (147,000 fragments, 2,087,476 slots, no truth) | reference `None` from 0 kernels: every weight 1, gDNA fraction 0.0850, pools mRNA 126,512 / nascent 7,699 / gDNA 10,109 / intergenic 2,357 | `land_med`: gDNA fraction 0.0997, pools 123,618 / 8,433 / 12,270 / 2,357; weights q50 / q90 / q99 = 0.00 / 0.01 / 0.01 (relative log), 0.49 % of objects above 2 nats, max +13.5 (`slab_bg`: 0.0979; `land_geo`: 0.0998) | the honest reader contracts the library the shipped reader cannot: half a percent of objects read enriched 7× or more with the bulk at 1, the EM's gDNA pool rises 21 %, mRNA falls 2.3 %, nascent rises 9.5 %. Whether that is right is unknowable here (no truth); it is the direction a captured plasma library should move. The NumPy prototype needed 4.8 h at four workers and 37 GB peak for this one library (§5); the mode readout was not yet built when this run started and was not repeated |

### 1.4 The one mechanism behind the in-scope cost, located

Three measurements, in order, on the benign stranded rows (log-density error against the simulator's own
off-target density `ρ0`, medians):

| row | uncaptured objects, all | uncaptured, no read | captured objects | objects with ≥ 5 true gDNA reads |
|---|---|---|---|---|
| g05 ss0.70 ON | −0.36 | — | +0.04 | +0.10 |
| g05 ss0.99 ON | −0.30 | −0.30 (n = 1,592) | +0.01 | +0.04 |
| g25 ss0.99 ON | −0.11 | — | −0.00 | −0.01 |
| g50 ss0.99 ON | −0.08 | −0.14 (n = 935) | −0.02 | −0.01 |
| g98 ss0.99 ON | −0.10 | — | −0.02 | −0.02 |
| g05 ss0.99 OFF | −0.05 | −0.15 (n = 603) | — | −0.08 |
| strength 100, g50 ss0.99 ON | −0.39 | −0.39 (n = 1,567) | −0.03 | +0.00 |
| depth 1/10, g05 ss0.99 ON | −0.79 | −0.79 (n = 1,999) | −0.03 | −0.14 |

1. **Wherever an object has data the reader is right**, including under strength-100 capture; the whole error
   is at objects with no read, and it grows as the background thins: 0.14 (g50) → 0.30 (g05) → 0.39 (strength
   100) → 0.79 (a tenth of the depth).
2. **The landscape's depleted mode is not the cause.** Staging the same conditions and reading the fitted prior
   directly: the depleted mode sits at +0.08 (g05 ON), +0.02 (g50 ON), −0.03 (g05 OFF) and +0.01 (strength 100)
   nats of `ρ0`. Only at a tenth of the depth does it collapse to the grid floor (−3.0 nats, 10 % of the mass).
   My first-version hypothesis that the depleted anchor sits below the background is refuted.
3. **The cause is the breadth of the smoothed prior under a median.** A zero-read object's posterior is the
   prior times `exp(−ρ·Eg)`; the exponential annihilates the enriched mode and leaves the depleted part, which
   is rendered from kernels of the training slots' own posterior widths (about one nat) and is broad and
   left-skewed: its ±0.5 nat window holds 24 % (g05), 48 % (g50) and 22 % (strength 100) of the mass. A
   synthetic zero-read object at the median opportunity, read under the landscape alone, has posterior median
   −0.30 (g05), −0.17 (g50), −0.39 (strength 100) and −0.79 (depth 1/10): exactly the measured values. The
   shipped reader reads the same breadth through a clipped arithmetic mean and lands on the other side
   (unprobed class +0.37 at g05 ss0.99).

The repair this implies needs no background and no new population: read the posterior **mode**. At a
data-poor object it is the depleted mode, within 0.1 nat of `ρ0` on every normal-depth row; at a data-rich
object it is the object's own likelihood peak, which the median already reads right. Falsification test
(calibration level, cached, a minute per condition): `land_map` must read no-read uncaptured objects within
±0.1 of `ρ0` where `land_med` reads −0.14…−0.39, leave captured objects within ±0.05 of their `land_med`
reading, leave the g25–g98 probed-within and partial classes unchanged, and not manufacture capture at
objects whose true factor is 1 (read per truth class, not from the realized-zero census, which also holds
short captured objects that happen to have no read). Result, the same lattice argmax refined by the parabola
through its neighbours:

| row | no-read uncaptured objects: med → map | captured objects: med → map | unprobed class: shipped / med / map | probed within: shipped / med / map | uncaptured exons q90 / max: med → map |
|---|---|---|---|---|---|
| g05 ss0.70 ON | −0.36 → +0.00 (n = 1,559) | +0.04 → +0.09 | +0.12 / −0.36 / −0.07 | 0.99 / 0.96 / 0.99 | −0.20 / +3.3 → +0.07 / +7.1 |
| g05 ss0.99 ON | −0.30 → +0.04 (1,592) | +0.01 → +0.05 | +0.37 / −0.21 / +0.10 | 0.99 / 0.94 / 0.96 | +0.24 / +6.6 → +0.47 / +7.0 |
| g25 ss0.70 ON | −0.17 → −0.02 (1,142) | −0.00 → +0.04 | +0.04 / −0.04 / +0.01 | 0.99 / 0.98 / 0.99 | −0.01 / +0.2 → +0.08 / +7.1 |
| g50 ss0.99 ON | −0.14 → −0.01 (935) | −0.02 → +0.02 | +0.04 / −0.01 / +0.03 | 0.99 / 0.99 / 0.99 | +0.09 / +0.7 → +0.12 / +7.0 |
| g05 ss0.99 OFF | −0.15 → −0.06 (603) | — | 0.00 / +0.01 / 0.00 | — | +0.21 / +1.0 → +0.13 / +1.0 |
| strength 100, g50 ss0.99 ON | −0.39 → −0.00 (1,567) | −0.03 → +0.00 | — | 0.99 / 0.99 / 0.99 | −0.14 / +9.3 → +0.04 / +9.3 |
| sparse layout, g50 ss0.99 ON | −0.15 → −0.04 (413) | −0.04 → −0.00 | +0.13 / +5.93 / +6.02 (touched cap: −0.94 both) | 1.00 / 1.00 / 1.00 | +0.08 / +5.8 → +0.07 / +5.9 |
| single probe, g05 ss0.99 ON | −0.15 → −0.05 (586) | −0.54 → −0.56 (the one gene, 2 objects) | partial: −4.59 / −0.10 / −0.07 | — | +0.20 / +0.9 → +0.12 / +0.9 |
| depth 1/10, g05 ss0.99 ON | −0.79 → −3.03 (1,999; the prior's floor) | −0.03 → +0.06 | — | 0.92 / 0.51 / 0.66 | −0.44 / +6.2 → −0.03 / +7.0 |

The gate passes on every normal-depth row: the systematic bias at data-poor objects is gone (−0.00…−0.06 on
every row, strength 100 included), captured objects stay within 0.05 of the median's reading, the g25–g98
classes are unchanged, the within-gene sd is the same or smaller (0.003–0.008), the step-halving movement
is ≤ 0.13 (the mode does not flip where the median flipped by 0.6–4.4), and the median error of every
uncaptured class is within ±0.05 of zero. Per object the centring artefact disappears with the bias:
captured exons read +0.09 / −0.03 / +0.07 / +0.01 / +0.01 relative on the five captured rows where the median
read +0.41 / +0.28 / +0.17 / +0.07 / +0.36. The price is the upper tail: the 90th percentile of uncaptured
introns is +0.17…+0.50 (the median's +0.22…+0.79) and of uncaptured exons +0.05…+0.42 (the median's
+0.11…+0.54), and the few RNA-rich exons whose DNA content is unidentified jump all the way to the enriched
mode (maxima +7) instead of part of the way. Whether that tail costs more than the bias it removes was the
end-to-end gate, arms `shipped`, `land_med`, `land_map`, `land_map_cap4` on one scan and one calibration per
condition (transcripts % / genes %):

| rows | shipped | land_med | land_map | land_map_cap4 |
|---|---|---|---|---|
| benign, stranded × ON, g05–g25 (4) | 6.27 / 0.79 | 10.17 / 1.25 | 6.94 / 0.87 | 6.95 / 0.87 |
| benign, stranded × OFF, g05–g25 (4) | 6.88 / 1.00 | 9.47 / 1.06 | 8.73 / 1.01 | 8.74 / 1.01 |
| benign, unstranded × OFF, g05–g25 (2) | 7.99 / 1.13 | 9.40 / 1.25 | 8.40 / 1.19 | 8.41 / 1.19 |
| benign, unstranded × ON, g05–g25 (2, deferred) | 10.89 / 2.32 | 27.73 / 4.82 | 12.33 / 2.62 | 12.33 / 2.62 |
| single probe, stranded × ON (4) | 9.22 / 2.12 | 10.94 / 1.28 | 9.81 / 1.22 | 9.81 / 1.22 |
| strength 100, stranded × ON (2) | 10.39 / 3.31 | 10.97 / 3.67 | 10.79 / 3.61 | 10.76 / 3.59 |
| strength 100, unstranded × ON (1, deferred) | 14.83 / 6.57 | 30.99 / 10.48 | 15.00 / 7.04 | 15.01 / 7.05 |
| junction-probed, stranded × ON, g25 and g50 ss0.99 | 32.98 / 1.64 | — | — | 25.01 / 2.18 |
| junction-probed, stranded × OFF, g25 and g50 ss0.99 | 8.28 / 1.07 | — | — | 9.00 / 1.10 |
| **the ladder** (genome scale, the release metric), stranded × ON g05 ss0.99 | 3.03 / 0.58 | 11.02 / 2.49 | 5.62 / 0.81 | 5.68 / 0.80 |
| the ladder, stranded × ON g50 ss0.99 | 4.48 / 2.14 | 7.03 / 2.25 | 6.64 / 2.28 | 5.97 / 1.79 |
| the ladder, stranded × ON g98 ss0.99 | 28.54 / 20.52 | 26.12 / 19.24 | 26.12 / 18.99 | 25.41 / 18.25 |
| the ladder, stranded × ON, mean (3) | 12.02 / 7.75 | 14.72 / 7.99 | 12.79 / 7.36 | 12.35 / 6.95 |
| the ladder, stranded × OFF g05 / g50 / g98 ss0.99 | 1.57 / 0.16; 2.33 / 0.28; 15.17 / 4.45 | 4.28 / 0.34; 4.48 / 0.38; 16.11 / 4.50 | 2.59 / 0.23; 4.10 / 0.36; 15.52 / 4.40 | 2.71 / 0.23; 4.12 / 0.36; 15.50 / 4.41 |
| the ladder, stranded × OFF, mean (3) | 6.35 / 1.63 | 8.29 / 1.74 | 7.41 / 1.66 | 7.44 / 1.67 |

On the release metric the mode with the touched cap costs +0.3 transcript points on the captured rows with
genes better by 0.8 (the g05 row +2.7 / +0.2, g50 +1.5 / −0.35, g98 −3.1 / −2.3) and +1.1 transcript points
on the capture-OFF rows with genes equal (g05 +1.1, g50 +1.8, g98 +0.3). The genome-scale probed-class spread
of §1.3 shows up as the g05 captured row and the OFF rows; nothing moves by more than the test chromosome's
cost, so the reservation of the verdict is not triggered.

Row by row on the captured benign rows: g05 ss0.70 6.56 → 6.08, g05 ss0.99 4.93 → 6.91, g25 ss0.70
6.73 → 7.50, g25 ss0.99 6.84 → 7.28, genes within +0.3. The mode passes: the low-gDNA loss is gone on the
captured rows, the deferred stratum is nearly recovered, the single probe's genes nearly halve, and strength
100 is a wash. The capture-OFF rows keep a cost of +0.4…+3.1 transcript points (g05 ss0.99 OFF 6.42 → 9.51
is the largest) with genes equal. On that row the mode's reading is unbiased in every class and every bin:
exons −0.01…+0.03 whether they carry 4 or 4,000 RNA reads and whether 0 or 10 gDNA reads, introns −0.05,
boundaries −0.01. What remains is the spread of an honest per-object reading: the exons' 10th–90th
percentiles span −0.07…+0.16 nats (the median's −0.11…+0.23), the introns' −0.18…+0.09, and the EM's
isoform split is sensitive to a within-gene spread of that size (`TRAPS: judge-a-ruler-by-its-within-gene-spread`)
where the shipped reader's `None` gives exactly 1. That is the price of having no detector; it is declared,
not hidden (§1.2), and it shrinks with gDNA depth (+0.4…+1.4 at g25 against +1.3…+3.1 at g05). The mode does nothing for the deferred stratum at g50 or for the
sparse-layout snap (§3).

**The genome-scale deficit is a different object class, and it is the cfRNA regime.** On the suite's
RNA-short captured row (gDNA fragments longer than RNA's, as in plasma), 3,342 captured exons are shorter
than a gDNA fragment: their gDNA opportunity is a fraction of a base, so no gDNA read can ever sit in them,
their four reads are RNA, and the composition row cannot place a density at a vanishing opportunity. Both
honest readouts return the prior there (median +0.54 / mode +0.62 against the truth, 10th–90th percentiles
−6…+1.8), and the message layer delivers a neighbour-imputed DNA level to 5 % of them (`zero_opportunity_probe.py`,
Appendix B), against 47 % of captured boundaries. Every object with opportunity reads right (−0.02 / −0.00).
This is not a readout question and no summary of the posterior fixes it: an object with no DNA evidence of
its own must take its capture from objects that are witnessed. The cheapest consumer-side rule was run and
is **refuted**: pricing every piece in which fewer than one gDNA fragment fits (16,554 of the row's 35,135
regions) at the larger of its two boundaries' weights, the shipped consumer's own intron fallback extended,
leaves the probed class where it was (0.67 → 0.66 under the mode, 0.67 → 0.62 under the median), moves the
unprobed class from +1.75 to −0.49, and changes nothing on any test-chromosome row (`short_piece_rule.py`).
The probed class's own distribution says where the deficit is: relative to its median, the shipped reader's
probed transcripts sit within [−0.06, +0.04] (10th–90th percentiles; 97 % within 0.1), the oracle's within
[−0.11, +0.07] (85 %), the mode's with the touched cap within [−0.16, +0.10] (68 %) and the median's within
[−0.31, +0.10] (67 %): a lower tail of probed transcripts read 15–30 % too short, which the boundary rule does
not touch, so the lost length is at witnessed objects reading low or at pieces whose boundaries read low with
them. Dissecting that tail transcript by transcript (which objects carry the lost share, and what they
received) is the next step and was not done here; the honest reader's claim stays confined to the test
chromosome's scale until it is; what it costs on the release metric is the ladder table above: +0.3
transcript points on the captured rows with genes better by 0.8, +1.1 on the capture-OFF rows with genes
equal. At a tenth of the depth
the prior itself collapses to its grid floor and the mode follows it there (−3.0 at every data-poor object;
the probed class 0.66 against the median's 0.51 and the shipped census's 0.92); no readout fixes a prior, and
that row stays a declared limit of a local reader.

## 2. The background objection (owner's point 2)

Agreed, and measured. The background-relative readout needs `(μ, α)` from the intergenic pool; at capture
strength 100 that pool holds 25–75 reads in the whole test chromosome, and on the plasma libraries it is the
same picture (LBX0190: 2,357 intergenic fragments of 147,000). The slab arm is the worst arm on every captured
stratum (benign stranded × ON 33.5 / 14.3 against 22.0 / 10.1 end to end; within-gene sd 0.010–0.033 against
0.003–0.005), and the per-object scoring shows why: its weights are unbiased at the object (captured exons
−0.01, data-poor objects 0.00) but dispersed by each object's own count (captured boundaries' 10th percentile
−1.0; q10–q90 −0.30…+0.15 against +0.10…+0.29 for `land_med`), the count-dependent shrinkage the first
version's Q2 predicted. Even where its absolute scale is right (strength 100: −0.03 at the object) it loses at
the transcript level (probed within 0.88).

The honest likelihood does not need a background. Under the landscape the count solver already fits it
anchors on the abundant population, which is the owner's premise exactly, and keeps every property the brief
wanted: the object's observations and licensed local factors in place of an inferred count read as an
observation; no detector, no reference density, no `None`, no clip; a weight of 1 only where the evidence says
so. Against the shipped reader the change is the likelihood and the summary, nothing else. The background
estimator stays deleted. The one thing the background did well, placing data-poor objects at a narrow,
correctly located level, the posterior mode does without it (§1.4).

Two limits of the landscape frame are measured, not hidden. The median's discontinuity: a bimodal posterior
with nearly balanced mass flips between modes under a tiny mass change (step-halving moves 0.6–4.4 in log
weight at isolated objects; the geometric mean's largest is 0.3). And the deferred stratum: with the strand
channel dead the local likelihood cannot separate enriched DNA from RNA at a both-strand object, so every
honest readout under-reads capture there by 10–30 points of transcripts. Both-strand objects in stranded
libraries are read right (benign ON: median error +0.0 at the object).

## 3. The transfer limitation, measured (owner's point 3)

Calibration level, stranded × ON (8 rows each), partial / unprobed class median error:

| arm | RNA shorter than gDNA (`fl_gdna_long`) | RNA longer (`fl_rna_long`) | equal, 200 bp (`fl_equal200`) | benign (equal, 206 bp) |
|---|---|---|---|---|
| shipped | +0.04 / −0.21 | +0.02 / +0.40 | +0.04 / +0.07 | +0.03 / +0.11 |
| land_med | +0.00 / −0.42 | −0.02 / +0.20 | +0.00 / −0.12 | −0.01 / −0.10 |
| slab_bg | −0.17 / −0.29 | −0.24 / +0.33 | −0.20 / +0.02 | — |

The unprobed class carries the transfer error, its sign follows the gap direction, and it is of the same size
in both readers (0.2–0.4 nats): RNA is priced at gDNA's capture at every object, which is right only when the
two components overlap probes alike. On RNA-longer gaps `land_med` halves the shipped error; on RNA-shorter
gaps it doubles it, because the data-poor under-read of §1.4 adds in the same direction. Per object the reader's
own error is unchanged across the arms (captured exons +0.20 / +0.13 / +0.17 relative on the three arms, the
same centring artefact of §1.4), so the arm-to-arm movement is the transfer, not the reader. End to end
(stranded × ON means): RNA shorter 24.5 / 12.0 → 27.4 / 12.9; RNA longer 21.0 / 11.1 → 24.0 / 11.7; capture-OFF
rows of both arms: genes improve (4.27 → 3.93, 4.69 → 4.45), transcripts cost +1…+2 points. The RNA-shorter
capture-OFF rows, where the genome-scale panel makes the shipped reader form a false reference (42 %
transcripts), read clean on the test chromosome for both readers (no mode forms at this size); the
genome-scale capture-OFF row reads clean for both as well (every shipped weight exactly 1; §1.3), and the
genome-scale captured row is where the honest readouts' probed-class spread appears (§1.3, §1.4).

Worst cases by probe placement (stranded × ON, 8 rows; partial / unprobed / within-gene sd):

| layout | shipped | land_med, pieces cap | land_med, touched cap | per object: captured exons, land_med |
|---|---|---|---|---|
| benign (probes tile exon interiors) | +0.03 / +0.11 / 0.004 | −0.01 / −0.10 / 0.004 | −0.01 / −0.10 / 0.005 (no junction-spanning probe: the caps agree) | +0.17 [+0.10, +0.29], 94 % within 0.5 nat |
| junction-probed (probes span junctions) | +1.22 / +3.21 / 0.034 | +5.23 / +6.59 / 0.176 | −0.22 / +0.30 / 0.051 | +0.08 [−0.25, +0.27], 88 % |
| sparse (one probe inside an exon) | −0.05 / +0.10 / 0.005 | +5.65 / +5.89 / 0.008 | −1.21 / −0.99 / 0.007 | −0.60 [−0.79, +0.16], 88 % |
| single probe | −4.59 / +0.00 / 0.000 | −0.04 / +0.00 / 0.007 | the same (the probe sits inside an exon); mode readout −0.07 | — |

Reading: the junction-probed failure was the cap, not the weights (the objects read right; the fully probed
transcripts were capped at their uncaptured pieces); the touched cap repairs it and beats the shipped reader
there. End to end over the eight stranded captured rows the pieces cap costs 39.2 / 7.5 → 52.8 / 20.1 and the
touched cap (median readout) gives 30.2 / 8.1: transcripts better by 9 points on average (g25 ss0.70
32.5 → 19.4, g50 ss0.70 34.9 → 23.5), genes within +0.5 (g05–g25 +0.4…+0.7, g50 ss0.70 3.1 → 4.2, g98
ss0.70 33.2 → 30.9); with the mode readout the two ss0.99 rows read 23.2 / 1.3 (g25) and 26.8 / 3.1 (g50).
On the deferred stratum the touched cap halves the error (141 / 104 → 82 / 41). This layout's transcript
floor is structural: every arm, the shipped reader included, reads 35–40 % on its zero-gDNA rows, where
every weight is 1 by construction. The
sparse failure is the weights: a partly covered exon sits between the landscape's modes and the median snaps
down (−0.60 at the object); with the touched cap the anchors then read too high relative to the rest. The
arithmetic mean interpolates there (partial +1.28 → −1.12 with the cap) but invents capture at zero-DNA
objects (q90 +4.49 on the same rows, +0.94 on the g00 control) and is rejected. The shipped reader's clip at
the located mode behaves like a global cap and reads both layouts tolerably; a global cap in relative units
needs a library-wide enrichment level, which is the detector again. Open item, not a transfer question: both
caps are the junction rule's layout dependence, which no gDNA witness removes
(`ISSUES: ruler-witness-geometry-on-transcript-panels`).

## 4. Deletions and the commit (owner's point 4)

Nothing was deleted or committed by this review. Conditional on a detector-free reader landing:

- delete: the located-mode census (`landscape.py` lines 394–498: `_census`, `split_basins`, `LandscapeMode`,
  `LocatedMode`, `located_enriched_mode`), `capture_efficiency.py`, the `gdna_reference_density` /
  `gdna_reference_members` fields and their validators, the CLI summary and report fields, the clip to `[0, 1]`
  in `_cut_efficiencies` and the `None` return in `transcript_capture_eff_lengths`, the `np.minimum(eff_len,
  span)` clamp in `assemble_priors`, `test_capture_efficiency.py`, the six mode tests in `test_landscape.py`,
  the two reference-schema tests;
- delete (archive): the research integrator family and the per-run import hooks; the representation campaign
  stays parked with its three red tests;
- keep deleted: `density_deconv` (the background estimator has no consumer: the slab is refuted);
- if the shipped reader is kept for 0.8.0 (recommendation (b)), nothing above is deleted; the located-mode
  census is then the production reader and stays.

The order the owner set stands: commit and push the current working tree (the count foundation, the reviewed
goldens, the docs) before any deletion. The tree was last certified at 2,814 passes by the landing receipts;
it must be re-certified (suite, ruff, `preflight.py --full`) immediately before that commit, which this review
did not run while the benchmarks held the machine.

## 5. What this benchmark cannot show

- Real libraries carry no truth beyond the VCaP mixture's gDNA fraction; the LBX0190 row reports pools,
  runtime and the weight census only, for the median readout; the other three libraries were not run (one
  whole-genome process at a time, and this one took the session's remaining hours).
- The sparse-probe layout and the single-probe panel are the only transcript-designed panels; their
  "unprobed" classes contain transcripts enriched by proximity, which the anchored score reads as error.
- A posterior median or mode under a flat-tailed prior will keep producing isolated large weights at RNA-rich
  exons whose DNA content is unidentified within one percent of their reads (maxima of +7 on g00 capture-OFF
  rows); end to end they cost nothing measurable, but they are not zero.
- Every end-to-end number comes from one harness and one seed; within a run every arm's EM reads the same
  scan and calibration, so within-row contrasts are exact. Across runs the control arm reproduced itself to
  ≤ 0.02 points of transcripts on six of the seven g05–g25 rows run twice and moved 0.45 on one (g05 ss0.70
  ON, 6.11 against 6.56): the pinned pair's own variation is a few tenths of a point on a sensitive row, and
  no claim below rests on a difference of that size.
- The reader prototype is NumPy and holds every slot's factors in memory: 30–90 s per test-chromosome
  condition at ten workers, 3–17 min per genome-scale suite row, and 4.8 h with 37 GB peak for the
  2.1 M-slot plasma library at four workers on the shared machine. A production reader is a native kernel
  inside the block solve, streaming with the blocks (§6); the prototype's cost says nothing about it, but it
  does say that no Python per-slot pass over a whole genome is shippable.

## 6. The numerical contract, and what the step-halving check caught

The first version's claim stands: the three integrals of the proper-prior readout, and the landscape readouts,
are weighted sums of `exp(log L)` on fixed lattices (the count solver's 0.2-nat step, one node per standard
deviation around the own term's mode and around the full integrand's coarse maximum, the same in the strand
share), with no adaptive integrator. The six frozen objects reproduce the archived adaptive reference to
within 9e-4 in log weight and the two both-strand objects that timed out complete in under a second. Every
condition re-evaluates every 40th slot and six deep both-strand slots with all steps halved; the check caught
two rules the first version did not state and would have shipped wrong:

1. the background integral's lattice must follow the Gamma's own scale on both sides of the mean (one node
   per `1/√α`), not the count's step: at shallow objects half the background's mass was integrated with two
   nodes (log-weight errors of 0.2–0.3 at g05);
2. the strand-share lattice must be clustered at the window step, or the per-density windows multiply
   without bound (a 2,135-read both-strand object cost 1.6 s instead of 0.5).

With both rules every single-strand slot converges to ≤ 0.03 under halving on every panel; the remaining
large halving changes (0.6–4.4) are median flips at isolated bimodal objects (§2), not lattice error.

## 7. Gates, stopping criteria, and the decisions that remain the owner's

Gates for landing a detector-free reader (each written and failing first; deliberate defects fire named
gates): the six frozen objects and the analytic limits (flat, pure DNA); the lattice rules of §6; bit identity
of the rebuilt tables with production's delivered rows on the cached conditions; the §1.4 data-poor gate; the
junction cap's three layouts; per-stratum end to end against this benchmark's `shipped` rows with the deferred
stratum reported; real libraries serially; then the deletions and the move rule.

Stop and report rather than add machinery if: the native kernel disagrees with the prototype's frozen
objects after one round of window-rule debugging; the sparse layout needs a library-wide level to read right
(that is the detector); any capture-OFF gene error rises on the full panels; a real library's calibration
time doubles.

The owner's decisions, each with its number beside it:

1. **Which reader ships 0.8.0.** (a) the mode readout with the touched cap: on the test chromosome it costs
   +0.7 / +0.1 on the captured g05–g25 rows, +0.4…+3.1 / 0.0 on the capture-OFF rows, and +0.4 / +0.3 under
   strength-100 capture; it reads the single probe (genes 2.1 → 1.2), the junction-probed layout
   (33.0 → 25.0) and the plasma library the shipped reader cannot (LBX0190, gDNA fraction 8.5 → 10.0 %); at
   genome scale it trails the shipped census on the probed class at calibration level (0.70–0.83 against
   0.96–0.98 within 0.1 on the ladder's captured rows), which on the ladder end to end costs +0.3 transcript
   points on the captured rows with genes better by 0.8 and +1.1 on the capture-OFF rows with genes equal. Or
   (b) the current reader with its two failures declared. My recommendation is (a); the trade is the owner's.
2. **Whether the depth-1/10 and deferred-stratum results are acceptable declared limits** of any local reader
   (they are structural: a collapsed prior and a dead strand channel; the mode recovers the deferred stratum at
   g05–g25 but not at g50).
3. **The commit and push of the current tree** before any deletion, as the owner specified.

## Appendix A — the tables

Rendered from the result files by `render_tables.py`; strata never pooled; a `—` is an empty class (no transcript
of that class with ≥ 20 RNA fragments). Transcript-level columns: the probed class's share within 0.1 nat, the partial and
unprobed classes' median log error, the within-gene sd; the census columns read realized-zero objects (which include short
captured objects with no read), so the per-truth-class scoring of §1.4 is the false-capture measure. Zero controls, the
deferred stratum, per-condition rows and the mode readout's per-object tables are in the archive named in Appendix B.

### A.1 Calibration level, every panel (`results/`): `shipped`, `land_med`, `slab_bg`; in-scope strata, zero controls and the deferred stratum omitted (they are in the archive)

**benign** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x OFF | 8 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| stranded x OFF | 8 | land_med | — | — | +0.00 | 0.005 | +0.00 / +0.45 | 0.02 |
| stranded x OFF | 8 | slab_bg | — | — | +0.00 | 0.002 | +0.00 / +1.29 | 0.00 |
| stranded x ON | 8 | shipped | 0.99 | +0.03 | +0.11 | 0.004 | — / — | — |
| stranded x ON | 8 | land_med | 0.98 | -0.01 | -0.10 | 0.004 | +0.18 / +7.23 | 0.71 |
| stranded x ON | 8 | slab_bg | 0.79 | -0.21 | +0.05 | 0.010 | +0.01 / +7.57 | 0.07 |
| unstranded x OFF | 4 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| unstranded x OFF | 4 | land_med | — | — | -0.00 | 0.002 | +0.00 / +0.15 | 0.00 |
| unstranded x OFF | 4 | slab_bg | — | — | -0.00 | 0.001 | +0.00 / +0.44 | 0.01 |

**probes_junction** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x ON | 8 | shipped | 1.00 | +1.22 | +3.21 | 0.034 | — / — | — |
| stranded x ON | 8 | land_med | 1.00 | +5.23 | +6.59 | 0.176 | +0.06 / +3.89 | 0.13 |
| stranded x ON | 8 | slab_bg | 1.00 | +5.25 | +6.67 | 0.190 | +0.00 / +1.64 | 0.04 |

**probes_sparse** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x ON | 8 | shipped | 1.00 | -0.05 | +0.10 | 0.005 | — / — | — |
| stranded x ON | 8 | land_med | 1.00 | +5.68 | +5.89 | 0.007 | +0.09 / +5.84 | 0.59 |
| stranded x ON | 8 | slab_bg | 1.00 | +5.30 | +5.99 | 0.013 | +0.00 / +2.59 | 0.02 |

**fl_equal200** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x OFF | 8 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| stranded x OFF | 8 | land_med | — | — | +0.00 | 0.005 | +0.01 / +0.79 | 0.01 |
| stranded x OFF | 8 | slab_bg | — | — | +0.00 | 0.002 | +0.00 / +1.97 | 0.02 |
| stranded x ON | 8 | shipped | 0.99 | +0.04 | +0.07 | 0.003 | — / — | — |
| stranded x ON | 8 | land_med | 0.98 | +0.00 | -0.12 | 0.004 | +0.19 / +7.29 | 3.08 |
| stranded x ON | 8 | slab_bg | 0.83 | -0.20 | +0.02 | 0.010 | +0.00 / +7.73 | 0.06 |
| unstranded x OFF | 4 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| unstranded x OFF | 4 | land_med | — | — | -0.00 | 0.002 | +0.00 / +0.12 | 0.00 |
| unstranded x OFF | 4 | slab_bg | — | — | -0.00 | 0.002 | +0.00 / +0.33 | 0.01 |

**fl_gdna_long** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x OFF | 8 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| stranded x OFF | 8 | land_med | — | — | +0.00 | 0.006 | +0.01 / +0.81 | 0.03 |
| stranded x OFF | 8 | slab_bg | — | — | +0.00 | 0.004 | +0.00 / +1.96 | 0.01 |
| stranded x ON | 8 | shipped | 0.99 | +0.04 | -0.21 | 0.003 | — / — | — |
| stranded x ON | 8 | land_med | 0.97 | +0.00 | -0.42 | 0.005 | +0.19 / +7.44 | 0.65 |
| stranded x ON | 8 | slab_bg | 0.79 | -0.17 | -0.29 | 0.012 | +0.01 / +7.96 | 0.02 |
| unstranded x OFF | 4 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| unstranded x OFF | 4 | land_med | — | — | -0.00 | 0.002 | +0.00 / +0.42 | 0.00 |
| unstranded x OFF | 4 | slab_bg | — | — | -0.00 | 0.003 | +0.01 / +0.95 | 0.01 |

**fl_rna_long** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x OFF | 8 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| stranded x OFF | 8 | land_med | — | — | -0.00 | 0.005 | +0.00 / +0.42 | 0.02 |
| stranded x OFF | 8 | slab_bg | — | — | +0.00 | 0.001 | +0.00 / +1.12 | 0.00 |
| stranded x ON | 8 | shipped | 0.99 | +0.02 | +0.40 | 0.003 | — / — | — |
| stranded x ON | 8 | land_med | 0.98 | -0.02 | +0.20 | 0.004 | +0.17 / +6.76 | 0.92 |
| stranded x ON | 8 | slab_bg | 0.79 | -0.24 | +0.33 | 0.010 | +0.01 / +4.77 | 0.17 |
| unstranded x OFF | 4 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| unstranded x OFF | 4 | land_med | — | — | -0.00 | 0.001 | +0.00 / +0.08 | 0.00 |
| unstranded x OFF | 4 | slab_bg | — | — | -0.00 | 0.001 | +0.00 / +0.40 | 0.01 |

**depth_d100** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x OFF | 5 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| stranded x OFF | 5 | land_med | — | — | +0.00 | 0.003 | +0.07 / +2.73 | 0.01 |
| stranded x OFF | 5 | slab_bg | — | — | +0.00 | 0.002 | +0.00 / +2.92 | 0.01 |
| stranded x ON | 5 | shipped | 0.77 | +0.52 | — | 0.022 | — / — | — |
| stranded x ON | 5 | land_med | 0.59 | +0.55 | — | 0.027 | +0.29 / +7.76 | 0.01 |
| stranded x ON | 5 | slab_bg | 0.45 | +0.64 | — | 0.040 | +0.04 / +8.46 | 0.01 |

**depth_d10** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x OFF | 5 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| stranded x OFF | 5 | land_med | — | — | +0.00 | 0.007 | +0.09 / +2.24 | 0.01 |
| stranded x OFF | 5 | slab_bg | — | — | +0.00 | 0.002 | +0.02 / +1.80 | 0.02 |
| stranded x ON | 5 | shipped | 0.82 | +0.33 | — | 0.011 | — / — | — |
| stranded x ON | 5 | land_med | 0.65 | -0.41 | — | 0.024 | +0.28 / +7.74 | 0.23 |
| stranded x ON | 5 | slab_bg | 0.48 | +0.21 | — | 0.035 | +0.01 / +7.59 | 0.03 |

**binding_0p1** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x ON | 2 | shipped | 0.99 | +0.02 | +0.03 | 0.003 | — / — | — |
| stranded x ON | 2 | land_med | 0.98 | -0.04 | +0.00 | 0.004 | +0.01 / +2.50 | 0.06 |
| stranded x ON | 2 | slab_bg | 0.83 | -0.05 | -0.02 | 0.006 | +0.00 / +2.80 | 0.03 |

**binding_1** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x ON | 2 | shipped | 0.99 | +0.01 | +0.01 | 0.003 | — / — | — |
| stranded x ON | 2 | land_med | 0.99 | -0.02 | -0.00 | 0.005 | +0.03 / +4.66 | 0.02 |
| stranded x ON | 2 | slab_bg | 0.90 | -0.29 | -0.02 | 0.006 | +0.00 / +3.24 | 0.02 |

**binding_100** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x ON | 2 | shipped | 0.99 | +0.01 | — | 0.003 | — / — | — |
| stranded x ON | 2 | land_med | 0.99 | -0.00 | — | 0.004 | +0.60 / +9.50 | 0.11 |
| stranded x ON | 2 | slab_bg | 0.88 | -0.31 | — | 0.007 | +0.00 / +1.01 | 0.03 |

**single_probe** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x ON | 4 | shipped | — | -4.59 | +0.00 | 0.000 | — / — | — |
| stranded x ON | 4 | land_med | — | -0.04 | +0.00 | 0.007 | +0.01 / +0.53 | 0.01 |
| stranded x ON | 4 | slab_bg | — | -0.14 | +0.00 | 0.003 | +0.00 / +1.20 | 0.01 |

**g00_seed43** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|

**flgap_rna_short** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x OFF | 1 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| stranded x OFF | 1 | land_med | — | — | -0.01 | 0.028 | +0.06 / +1.27 | 0.01 |
| stranded x OFF | 1 | slab_bg | — | — | +0.00 | 0.017 | +0.31 / +4.13 | 0.04 |
| stranded x ON | 1 | shipped | 0.96 | +0.24 | +2.98 | 0.080 | — / — | — |
| stranded x ON | 1 | land_med | 0.67 | +0.20 | -0.41 | 0.092 | +6.51 / +8.60 | 4.13 |
| stranded x ON | 1 | slab_bg | 0.36 | +0.18 | -0.46 | 0.132 | +0.56 / +9.91 | 0.05 |
| unstranded x OFF | 1 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| unstranded x OFF | 1 | land_med | — | — | -0.01 | 0.009 | +0.00 / +0.27 | 0.01 |
| unstranded x OFF | 1 | slab_bg | — | — | +0.00 | 0.001 | +0.00 / +2.07 | 0.01 |

**flgap_rna_long** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x OFF | 1 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| stranded x OFF | 1 | land_med | — | — | -0.01 | 0.031 | +0.03 / +0.25 | 0.02 |
| stranded x OFF | 1 | slab_bg | — | — | -0.00 | 0.007 | +0.00 / +1.38 | 0.07 |
| stranded x ON | 1 | shipped | 0.98 | +0.00 | +1.02 | 0.036 | — / — | — |
| stranded x ON | 1 | land_med | 0.82 | -0.19 | +0.61 | 0.075 | +0.11 / +6.22 | 0.24 |
| stranded x ON | 1 | slab_bg | 0.62 | -0.24 | +0.70 | 0.110 | +0.00 / +5.56 | 0.35 |
| unstranded x OFF | 1 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| unstranded x OFF | 1 | land_med | — | — | -0.00 | 0.008 | +0.00 / +0.10 | 0.01 |
| unstranded x OFF | 1 | slab_bg | — | — | +0.00 | 0.001 | +0.00 / +0.86 | 0.01 |

**ladder** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x OFF | 3 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| stranded x OFF | 3 | land_med | — | — | -0.01 | 0.033 | +0.05 / +1.42 | 0.03 |
| stranded x OFF | 3 | slab_bg | — | — | -0.00 | 0.011 | +0.19 / +3.68 | 0.22 |
| stranded x ON | 3 | shipped | 0.97 | +0.08 | +0.74 | 0.032 | — / — | — |
| stranded x ON | 3 | land_med | 0.67 | +0.03 | -0.21 | 0.058 | +0.06 / +8.51 | 0.40 |
| stranded x ON | 3 | slab_bg | 0.39 | +0.01 | +0.19 | 0.099 | +0.33 / +9.79 | 0.24 |
| unstranded x OFF | 2 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| unstranded x OFF | 2 | land_med | — | — | -0.00 | 0.011 | +0.01 / +0.26 | 0.01 |
| unstranded x OFF | 2 | slab_bg | — | — | +0.00 | 0.000 | +0.00 / +2.03 | 0.01 |

### A.2 Calibration level, the mode readout (`results_map/`): the decisive test-chromosome conditions and the suite's rows

**benign** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x OFF | 1 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| stranded x OFF | 1 | land_med | — | — | +0.01 | 0.012 | +0.03 / +0.44 | 0.01 |
| stranded x OFF | 1 | land_map | — | — | +0.00 | 0.008 | +0.01 / +0.21 | 0.01 |
| stranded x OFF | 1 | land_map_cap4 | — | — | +0.00 | 0.008 | +0.01 / +0.21 | 0.01 |
| stranded x ON | 4 | shipped | 0.99 | +0.05 | +0.14 | 0.003 | — / — | — |
| stranded x ON | 4 | land_med | 0.97 | -0.01 | -0.16 | 0.005 | +0.26 / +7.23 | 0.09 |
| stranded x ON | 4 | land_map | 0.98 | +0.03 | +0.02 | 0.004 | +6.91 / +7.06 | 0.12 |
| stranded x ON | 4 | land_map_cap4 | 0.98 | +0.03 | +0.02 | 0.005 | +6.91 / +7.06 | 0.12 |

**probes_sparse** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x ON | 1 | shipped | 1.00 | +0.02 | +0.13 | 0.003 | — / — | — |
| stranded x ON | 1 | land_med | 1.00 | +5.91 | +5.93 | 0.003 | +0.07 / +5.82 | 0.02 |
| stranded x ON | 1 | land_map | 1.00 | +6.00 | +6.02 | 0.003 | -0.01 / +5.85 | 0.02 |
| stranded x ON | 1 | land_map_cap4 | 1.00 | -0.96 | -0.94 | 0.003 | -0.01 / +5.85 | 0.02 |

**depth_d10** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x ON | 1 | shipped | 0.92 | +0.15 | — | 0.009 | — / — | — |
| stranded x ON | 1 | land_med | 0.51 | +0.05 | — | 0.017 | +0.66 / +7.74 | 0.17 |
| stranded x ON | 1 | land_map | 0.66 | +0.12 | — | 0.014 | +10.04 / +10.22 | 0.09 |
| stranded x ON | 1 | land_map_cap4 | 0.66 | +0.12 | — | 0.016 | +10.04 / +10.22 | 0.09 |

**binding_100** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x ON | 1 | shipped | 0.99 | +0.00 | — | 0.003 | — / — | — |
| stranded x ON | 1 | land_med | 0.99 | -0.00 | — | 0.005 | +0.50 / +9.23 | 0.11 |
| stranded x ON | 1 | land_map | 0.99 | -0.00 | — | 0.005 | +9.12 / +9.12 | 0.13 |
| stranded x ON | 1 | land_map_cap4 | 0.99 | -0.00 | — | 0.005 | +9.12 / +9.12 | 0.13 |

**single_probe** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x ON | 1 | shipped | — | -4.59 | +0.00 | 0.000 | — / — | — |
| stranded x ON | 1 | land_med | — | -0.10 | +0.00 | 0.011 | +0.02 / +0.53 | 0.01 |
| stranded x ON | 1 | land_map | — | -0.07 | +0.00 | 0.006 | +0.01 / +0.59 | 0.01 |
| stranded x ON | 1 | land_map_cap4 | — | -0.07 | +0.00 | 0.006 | +0.01 / +0.59 | 0.01 |

**flgap_rna_short** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x ON | 1 | shipped | 0.96 | +0.24 | +2.98 | 0.080 | — / — | — |
| stranded x ON | 1 | land_med | 0.67 | +0.20 | -0.41 | 0.092 | +6.51 / +8.60 | 4.13 |
| stranded x ON | 1 | land_map | 0.64 | +0.21 | +1.76 | 0.087 | +6.91 / +8.80 | 0.05 |
| stranded x ON | 1 | land_map_cap4 | 0.67 | +0.21 | +1.75 | 0.083 | +6.91 / +8.80 | 0.05 |

**flgap_rna_long** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x ON | 1 | shipped | 0.98 | +0.00 | +1.02 | 0.036 | — / — | — |
| stranded x ON | 1 | land_med | 0.82 | -0.19 | +0.61 | 0.075 | +0.11 / +6.22 | 0.24 |
| stranded x ON | 1 | land_map | 0.85 | -0.19 | +0.64 | 0.072 | +0.00 / +6.59 | 0.17 |
| stranded x ON | 1 | land_map_cap4 | 0.87 | -0.17 | +0.63 | 0.064 | +0.00 / +6.59 | 0.17 |

**ladder** (calibration level: transcript lengths against the simulator's yield; means over the stratum's conditions)

| stratum | n | arm | probed within 0.1 | partial med | unprobed med | within-gene sd | false capture: zero-DNA objects q90 / max log w | step halving max |
|---|---|---|---|---|---|---|---|---|
| stranded x OFF | 3 | shipped | — | — | +0.00 | 0.000 | — / — | — |
| stranded x OFF | 3 | land_med | — | — | -0.01 | 0.033 | +0.05 / +1.42 | 0.03 |
| stranded x OFF | 3 | land_map | — | — | -0.00 | 0.018 | +0.01 / +1.53 | 0.14 |
| stranded x OFF | 3 | land_map_cap4 | — | — | -0.00 | 0.020 | +0.01 / +1.53 | 0.14 |
| stranded x ON | 3 | shipped | 0.97 | +0.08 | +0.74 | 0.032 | — / — | — |
| stranded x ON | 3 | land_med | 0.67 | +0.03 | -0.21 | 0.058 | +0.06 / +8.51 | 0.40 |
| stranded x ON | 3 | land_map | 0.74 | +0.02 | +0.54 | 0.052 | +0.51 / +7.49 | 0.13 |
| stranded x ON | 3 | land_map_cap4 | 0.76 | +0.02 | +0.53 | 0.046 | +0.51 / +7.49 | 0.13 |

### A.3 End to end, every panel (`e2e/`): stratum means

**benign** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 | land_mean_cap4 | land_geo | slab_bg |
|---|---|---|---|---|---|---|---|
| stranded x OFF | 8 | 14.86 / 3.96 | 16.29 / 4.38 | — | — | 16.99 / 4.48 | 16.82 / 4.18 |
| stranded x ON | 8 | 21.97 / 10.05 | 24.40 / 10.60 | — | — | 27.09 / 12.09 | 33.50 / 14.28 |
| unstranded x OFF | 4 | 16.15 / 4.63 | 17.20 / 4.82 | — | — | 17.64 / 5.16 | 17.89 / 5.01 |

**probes_junction** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 | land_mean_cap4 | land_geo | slab_bg |
|---|---|---|---|---|---|---|---|
| stranded x ON | 8 | 39.20 / 7.53 | 52.83 / 20.09 | — | — | 53.23 / 20.57 | 53.74 / 20.57 |

**probes_sparse** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 | land_mean_cap4 | land_geo | slab_bg |
|---|---|---|---|---|---|---|---|
| stranded x ON | 8 | 36.17 / 12.26 | 38.27 / 13.40 | — | — | 40.82 / 15.93 | 67.32 / 40.58 |

**fl_gdna_long** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 | land_mean_cap4 | land_geo | slab_bg |
|---|---|---|---|---|---|---|---|
| stranded x OFF | 8 | 15.37 / 4.27 | 17.36 / 3.93 | — | — | 18.18 / 3.94 | 19.34 / 4.27 |
| stranded x ON | 8 | 24.54 / 11.99 | 27.36 / 12.89 | — | — | 29.29 / 13.98 | 35.80 / 16.57 |
| unstranded x OFF | 4 | 17.16 / 5.30 | 17.33 / 4.96 | — | — | 17.59 / 4.95 | 18.79 / 5.33 |

**fl_rna_long** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 | land_mean_cap4 | land_geo | slab_bg |
|---|---|---|---|---|---|---|---|
| stranded x OFF | 8 | 15.31 / 4.69 | 16.21 / 4.45 | — | — | 16.72 / 4.46 | 17.34 / 4.71 |
| stranded x ON | 8 | 20.99 / 11.11 | 24.00 / 11.70 | — | — | 24.87 / 12.17 | 30.80 / 14.56 |
| unstranded x OFF | 4 | 17.78 / 7.22 | 17.56 / 7.00 | — | — | 18.11 / 7.00 | 18.34 / 7.50 |

**depth_d100** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 | land_mean_cap4 | land_geo | slab_bg |
|---|---|---|---|---|---|---|---|
| stranded x OFF | 5 | 28.75 / 2.31 | 30.20 / 2.45 | — | — | 30.03 / 2.46 | 29.98 / 2.40 |
| stranded x ON | 5 | 34.31 / 3.20 | 36.42 / 3.43 | — | — | 37.96 / 3.48 | 39.37 / 4.38 |

**depth_d10** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 | land_mean_cap4 | land_geo | slab_bg |
|---|---|---|---|---|---|---|---|
| stranded x OFF | 5 | 16.93 / 1.30 | 18.88 / 1.35 | — | — | 19.88 / 1.35 | 18.89 / 1.33 |
| stranded x ON | 5 | 18.24 / 1.22 | 22.98 / 1.39 | — | — | 24.19 / 1.45 | 26.49 / 2.50 |

**binding_0p1** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 | land_mean_cap4 | land_geo | slab_bg |
|---|---|---|---|---|---|---|---|
| stranded x ON | 2 | 10.30 / 1.25 | 11.13 / 1.59 | — | — | 12.33 / 1.73 | 17.27 / 2.36 |

**binding_1** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 | land_mean_cap4 | land_geo | slab_bg |
|---|---|---|---|---|---|---|---|
| stranded x ON | 2 | 10.61 / 1.83 | 10.73 / 1.95 | — | — | 11.14 / 2.32 | 15.06 / 3.41 |

**binding_100** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 | land_mean_cap4 | land_geo | slab_bg |
|---|---|---|---|---|---|---|---|
| stranded x ON | 2 | 10.60 / 3.31 | 10.81 / 3.67 | — | — | 10.66 / 3.65 | 15.98 / 5.47 |

**single_probe** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 | land_mean_cap4 | land_geo | slab_bg |
|---|---|---|---|---|---|---|---|
| stranded x ON | 4 | 9.22 / 2.12 | 10.94 / 1.28 | — | — | 11.79 / 1.33 | 11.62 / 1.27 |

**g00_seed43** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 | land_mean_cap4 | land_geo | slab_bg |
|---|---|---|---|---|---|---|---|

### A.4 End to end: the mode readout's gate (`e2e_map/`), the touched cap (`e2e_pass3/`), the ladder (`e2e_suite/`)

**benign** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_map | land_map_cap4 |
|---|---|---|---|---|---|
| stranded x OFF | 4 | 6.88 / 1.00 | 9.47 / 1.06 | 8.73 / 1.01 | 8.74 / 1.01 |
| stranded x ON | 4 | 6.27 / 0.79 | 10.17 / 1.25 | 6.94 / 0.87 | 6.95 / 0.87 |
| unstranded x OFF | 2 | 7.99 / 1.13 | 9.40 / 1.25 | 8.40 / 1.19 | 8.41 / 1.19 |

**probes_junction** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_map | land_map_cap4 |
|---|---|---|---|---|---|
| stranded x OFF | 2 | 8.28 / 1.07 | — | — | 9.00 / 1.10 |
| stranded x ON | 2 | 32.98 / 1.64 | — | — | 25.01 / 2.18 |

**binding_100** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_map | land_map_cap4 |
|---|---|---|---|---|---|
| stranded x ON | 2 | 10.39 / 3.31 | 10.97 / 3.67 | 10.79 / 3.61 | 10.76 / 3.59 |

**single_probe** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_map | land_map_cap4 |
|---|---|---|---|---|---|
| stranded x ON | 4 | 9.22 / 2.12 | 10.94 / 1.28 | 9.81 / 1.22 | 9.81 / 1.22 |

**benign** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 |
|---|---|---|---|---|
| stranded x ON | 8 | 22.02 / 10.05 | 24.37 / 10.60 | 24.39 / 10.60 |

**probes_junction** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 |
|---|---|---|---|---|
| stranded x ON | 8 | 39.20 / 7.53 | 52.83 / 20.09 | 30.17 / 8.05 |

**probes_sparse** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 |
|---|---|---|---|---|
| stranded x ON | 8 | 36.17 / 12.26 | 38.27 / 13.40 | 38.26 / 13.47 |

**fl_gdna_long** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 |
|---|---|---|---|---|
| stranded x ON | 8 | 24.54 / 11.99 | 27.36 / 12.89 | 27.38 / 12.90 |

**fl_rna_long** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 |
|---|---|---|---|---|
| stranded x ON | 8 | 20.99 / 11.11 | 24.00 / 11.70 | 23.77 / 11.46 |

**binding_100** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 |
|---|---|---|---|---|
| stranded x ON | 2 | 10.39 / 3.31 | 10.97 / 3.67 | 10.94 / 3.65 |

**single_probe** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_med_cap4 |
|---|---|---|---|---|
| stranded x ON | 4 | 9.22 / 2.12 | 10.94 / 1.28 | 10.89 / 1.28 |

**ladder** (end to end, pinned scan, fractional assignment, seed 0: transcripts % / genes %, means over the stratum's conditions)

| stratum | n | base | land_med | land_map | land_map_cap4 |
|---|---|---|---|---|---|
| stranded x OFF | 3 | 6.35 / 1.63 | 8.29 / 1.74 | 7.41 / 1.66 | 7.44 / 1.67 |
| stranded x ON | 3 | 12.02 / 7.75 | 14.72 / 7.99 | 12.79 / 7.36 | 12.35 / 6.95 |

### A.5 Per object against the simulator's gDNA yield (`results/`, median and slab readouts; captured rows only)

**benign** (per object, against the simulator's own gDNA yield at the object; means over the stratum's conditions)

| stratum | n | readout | captured exons: median log error [q10, q90] / share within 0.5 nat | captured boundaries: median [q10, q90] | uncaptured exons: median / q90 (false capture) | uncaptured introns: q90 |
|---|---|---|---|---|---|---|
| stranded x ON | 8 | med | +0.17 [+0.10, +0.29] / 0.94 | +0.16 [+0.04, +0.25] | +0.00 / +0.20 | +0.35 |
| stranded x ON | 8 | slab | -0.01 [-0.30, +0.15] / 0.92 | -0.08 [-1.00, +0.05] | +0.00 / +0.02 | +0.00 |

**probes_junction** (per object, against the simulator's own gDNA yield at the object; means over the stratum's conditions)

| stratum | n | readout | captured exons: median log error [q10, q90] / share within 0.5 nat | captured boundaries: median [q10, q90] | uncaptured exons: median / q90 (false capture) | uncaptured introns: q90 |
|---|---|---|---|---|---|---|
| stranded x ON | 8 | med | +0.08 [-0.25, +0.27] / 0.88 | +0.00 [-0.75, +0.08] | +0.01 / +0.14 | +0.13 |
| stranded x ON | 8 | slab | -0.06 [-0.36, +0.15] / 0.91 | -0.05 [-0.80, +0.03] | +0.00 / +0.05 | +0.00 |

**probes_sparse** (per object, against the simulator's own gDNA yield at the object; means over the stratum's conditions)

| stratum | n | readout | captured exons: median log error [q10, q90] / share within 0.5 nat | captured boundaries: median [q10, q90] | uncaptured exons: median / q90 (false capture) | uncaptured introns: q90 |
|---|---|---|---|---|---|---|
| stranded x ON | 8 | med | -0.60 [-0.79, +0.16] / 0.88 | -0.22 [-1.06, +0.15] | +0.01 / +0.14 | +0.14 |
| stranded x ON | 8 | slab | -0.74 [-1.56, +0.02] / 0.80 | -0.40 [-1.46, +0.06] | +0.00 / +0.05 | +0.00 |

**fl_gdna_long** (per object, against the simulator's own gDNA yield at the object; means over the stratum's conditions)

| stratum | n | readout | captured exons: median log error [q10, q90] / share within 0.5 nat | captured boundaries: median [q10, q90] | uncaptured exons: median / q90 (false capture) | uncaptured introns: q90 |
|---|---|---|---|---|---|---|
| stranded x ON | 8 | med | +0.20 [+0.13, +0.38] / 0.84 | +0.18 [+0.08, +0.26] | +0.01 / +0.21 | +0.35 |
| stranded x ON | 8 | slab | +0.01 [-0.33, +0.19] / 0.91 | -0.06 [-0.98, +0.07] | +0.00 / +0.04 | +0.00 |

**fl_rna_long** (per object, against the simulator's own gDNA yield at the object; means over the stratum's conditions)

| stratum | n | readout | captured exons: median log error [q10, q90] / share within 0.5 nat | captured boundaries: median [q10, q90] | uncaptured exons: median / q90 (false capture) | uncaptured introns: q90 |
|---|---|---|---|---|---|---|
| stranded x ON | 8 | med | +0.13 [+0.07, +0.20] / 0.97 | +0.19 [-0.07, +0.30] | -0.02 / +0.17 | +0.28 |
| stranded x ON | 8 | slab | -0.04 [-0.39, +0.10] / 0.93 | -0.13 [-1.07, +0.06] | +0.00 / +0.02 | +0.00 |

**depth_d100** (per object, against the simulator's own gDNA yield at the object; means over the stratum's conditions)

| stratum | n | readout | captured exons: median log error [q10, q90] / share within 0.5 nat | captured boundaries: median [q10, q90] | uncaptured exons: median / q90 (false capture) | uncaptured introns: q90 |
|---|---|---|---|---|---|---|
| stranded x ON | 5 | med | -3.85 [-4.23, -2.37] / 0.19 | -3.86 [-5.51, -2.45] | -0.05 / +0.09 | +0.44 |
| stranded x ON | 5 | slab | -4.35 [-6.58, -2.88] / 0.21 | -4.68 [-6.83, -3.33] | +0.00 / +0.00 | +0.13 |

**depth_d10** (per object, against the simulator's own gDNA yield at the object; means over the stratum's conditions)

| stratum | n | readout | captured exons: median log error [q10, q90] / share within 0.5 nat | captured boundaries: median [q10, q90] | uncaptured exons: median / q90 (false capture) | uncaptured introns: q90 |
|---|---|---|---|---|---|---|
| stranded x ON | 5 | med | -1.92 [-2.52, -0.79] / 0.26 | -2.31 [-2.54, -0.84] | -0.04 / +0.21 | +0.71 |
| stranded x ON | 5 | slab | -2.78 [-3.99, -1.39] / 0.48 | -2.92 [-4.23, -1.67] | +0.00 / +0.00 | +0.05 |

**binding_0p1** (per object, against the simulator's own gDNA yield at the object; means over the stratum's conditions)

| stratum | n | readout | captured exons: median log error [q10, q90] / share within 0.5 nat | captured boundaries: median [q10, q90] | uncaptured exons: median / q90 (false capture) | uncaptured introns: q90 |
|---|---|---|---|---|---|---|
| stranded x ON | 2 | med | +0.04 [-0.04, +0.15] / 0.97 | +0.06 [-0.05, +0.15] | +0.02 / +0.10 | +0.04 |
| stranded x ON | 2 | slab | +0.00 [-0.18, +0.17] / 0.93 | -0.04 [-0.20, +0.09] | +0.00 / +0.01 | +0.00 |

**binding_1** (per object, against the simulator's own gDNA yield at the object; means over the stratum's conditions)

| stratum | n | readout | captured exons: median log error [q10, q90] / share within 0.5 nat | captured boundaries: median [q10, q90] | uncaptured exons: median / q90 (false capture) | uncaptured introns: q90 |
|---|---|---|---|---|---|---|
| stranded x ON | 2 | med | +0.05 [-0.02, +0.19] / 0.98 | +0.05 [-0.05, +0.13] | +0.02 / +0.11 | +0.08 |
| stranded x ON | 2 | slab | +0.01 [-0.09, +0.13] / 0.95 | -0.02 [-0.12, +0.08] | +0.00 / +0.00 | +0.00 |

**binding_100** (per object, against the simulator's own gDNA yield at the object; means over the stratum's conditions)

| stratum | n | readout | captured exons: median log error [q10, q90] / share within 0.5 nat | captured boundaries: median [q10, q90] | uncaptured exons: median / q90 (false capture) | uncaptured introns: q90 |
|---|---|---|---|---|---|---|
| stranded x ON | 2 | med | +0.37 [+0.32, +0.51] / 0.88 | +0.34 [+0.26, +0.42] | -0.02 / +0.17 | +0.74 |
| stranded x ON | 2 | slab | -0.03 [-0.14, +0.09] / 0.94 | -0.11 [-0.20, -0.02] | +0.00 / +0.00 | +0.00 |

**single_probe** (per object, against the simulator's own gDNA yield at the object; means over the stratum's conditions)

| stratum | n | readout | captured exons: median log error [q10, q90] / share within 0.5 nat | captured boundaries: median [q10, q90] | uncaptured exons: median / q90 (false capture) | uncaptured introns: q90 |
|---|---|---|---|---|---|---|
| stranded x ON | 4 | med | — [—, —] / — | — [—, —] | +0.03 / +0.16 | +0.05 |
| stranded x ON | 4 | slab | — [—, —] / — | — [—, —] | +0.00 / +0.11 | +0.00 |

### A.6 The real library (`real/`)

**lbx0190** (real library, in process; the pools and the gDNA fraction, no truth except the VCaP mixture)

| arm | gDNA fraction | mRNA / nascent / gDNA (EM) / intergenic | reader s | EM s | peak GB | weights: q50 / q90 / q99 / max (log, relative) |
|---|---|---|---|---|---|---|
| base | 0.0850 | 126512 / 7699 / 10109 / 2357 | 0 | 2 | 3.7 | — |
| slab_bg | 0.0979 | 124098 / 8219 / 12003 / 2357 | 17284 | 5 | 37.2 | +0.00 / +0.00 / +0.00 / +13.43 |
| land_med | 0.0997 | 123618 / 8433 / 12270 / 2357 | 0 | 5 | 37.2 | +0.00 / +0.01 / +0.01 / +13.47 |
| land_geo | 0.0998 | 123665 / 8378 / 12277 / 2357 | 0 | 5 | 37.2 | +0.00 / +0.01 / +0.01 / +13.28 |

## Appendix B — the prototype, the receipts, and how to re-run it

Archived at `.cache/rigel_runs/2026-10-10_capture_reader_benchmark/` (git-ignored, beside the earlier campaign receipts); the
session scratch directory it was copied from is not in the repository:

| file | role |
|---|---|
| `reader_proto.py` | the evaluator (`log_L`, fixed lattices), the readouts (`slab`, `geo`, `med`, `mean`, `map`), the consumers (`lengths_from_weights`, caps `pieces` / `touched`), the arms, the identity cache |
| `bench.py` | calibration-level driver: `python bench.py OUT_DIR LABEL[:filter] …`; writes `LABEL__COND.json` (headline, classes, within-gene, census, halving, timing) and `.npz` (per-slot arrays) |
| `object_truth.py`, `object_truth_all.py` | per-object truth from the sampler's capture law; `LABEL__COND__objects.json` and `__truth.npz` |
| `landscape_probe.py` | the fitted prior's modes against the true background; the synthetic zero-read object (§1.4) |
| `e2e.py`, `e2e_real.py` | in-process end to end: `python e2e.py OUT.jsonl PANEL_DIR INDEX_DIR arm,arm,… COND …` |
| `render_tables.py`, `aggregate_e2e.py` | the tables |
| `results/`, `results_cap4/`, `results_pass3/`, `results_map/`, `e2e/`, `e2e_pass3/`, `real/` | the result files; `sweep*.log`, `e2e/*.log`, `objects_*.log` the receipts |
| `review_v1_2026-10-09.md` | the first version of this review, superseded |
| `short_piece_rule.py`, `zero_opportunity_probe.py` | the genome-scale dissection of §1.4: the refuted boundary-pricing rule and the delivered-factor census |
| `tables_*.md`, `appendix_a.md` | every rendered table |

Environment: the `rigel` conda env, `OMP_NUM_THREADS=1`, `READER_WORKERS=4` on a shared machine; every panel
from `~/Downloads/rigel_runs/test_reference/` and `~/Downloads/rigel_runs/suite/`, every real BAM from
`~/Downloads/rigel_runs/cfrna/`. Repository state at the end: only this file differs from the owner's working
tree as it was at the start (verified with `git status`).
