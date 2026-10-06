# The RNA-short defect: what the fix requires, and the implementation plan

*Sandbox document (`docs/dev/`): provisional, not authoritative, cited by nothing outside the sandbox. Written
2026-10-05 for an external reviewer who has not seen this project. Every number was measured on 2026-10-05 unless a
date says otherwise, pinned (the BAM scan on one thread, fractional EM assignment) and scored against per-fragment
simulation truth. The prototypes and their data are in `~/Downloads/rigel_runs/prototypes/2026-10-05_spectrum_ruler/`;
`SPECTRUM_RULER.md` beside this file holds the full ledger of the ruler prototype in §3.4. What we ask of you is §7.*

## 1. Rigel, in the terms this document uses

Rigel quantifies RNA-seq transcripts from libraries contaminated with genomic DNA (gDNA), in three stages.

1. **The scan** deposits every fragment on genomic *objects*: *regions* (an exon, an intron, an intergenic stretch) and
   *boundaries* (the line between two adjacent regions). A fragment wholly inside a region deposits there; one that
   crosses a boundary deposits on the boundary.
2. **Calibration** solves every object for its composition: gDNA, RNA on the + strand, RNA on the − strand. The
   per-object solver is called **ψ**. Its evidence is the read strand (in a stranded library RNA has one orientation,
   gDNA both), a population prior over gDNA density called **the landscape** (fitted from the objects' own solves and
   fed back in a refit loop), and messages between neighbouring objects. Its output per object is a gDNA count `k` on
   the object's **gDNA opportunity** `S`: the number of positions where a gDNA fragment can start inside it, which
   depends on gDNA's fragment-length distribution. A 100 bp exon has almost no opportunity for 250 bp fragments.
3. **A per-locus EM** assigns RNA to transcripts and gDNA to a gDNA component, dividing counts by effective lengths.

**Capture.** Hybrid-capture libraries enrich probed exons, about 1,000-fold in the simulations. Rigel gets no probe
file; it learns capture from gDNA, which is uniform before capture. The EM's effective lengths are contracted by
per-object capture weights; the component that computes them is called **the ruler**. The EM's per-component weight is
`log θ − log L`, so a factor common to every length in a locus cancels: only relative weights matter.

**The shipped ruler.** Each object's weight is the posterior mean of `min(ρ / ρ_ref, 1)` under the landscape, where
`ρ_ref` is the landscape's "located enriched mode". When no such mode is found, `ρ_ref` is `None` and every weight is
exactly 1.

**Panels.** A 30-condition test chromosome (seconds per condition); a 16-condition genome-scale *ladder* whose gDNA and
RNA share one fragment-length law by design; two *gap arms* on the same genome, RNA-short (RNA 78 bp, gDNA 250 bp) and
RNA-long (the reverse); four real libraries (three plasma cfRNA, one VCaP mix whose reads carry their origin). Results
are read per stratum: stranded or unstranded, capture on or off. Unstranded × capture-on is reported but deferred.

**Truth.** An oracle splits each simulated BAM by read origin and runs the production accumulator on each part. That
gives certified per-object gDNA counts (`slot_truth.npz`), against which any object's calibration can be checked.

## 2. The defect and the constraints

On the RNA-short arm, the stranded × capture-off library reads transcript error **42.30 %** (genes 3.28 %) on the
current tree against **6.31 %** (0.44 %) on the 0.7.1 release.

The owner's constraints:

1. Achieve or exceed 0.7.1 on every stratum of every panel, robust to a fragment-length gap in either direction.
2. **Capture is a spectrum.** No capture detector, no on/off status, no "reference or `None`", no library-level test of
   capture. Every capture weight is per object and continuous. A plasma panel enriches scarce transcripts that stay a
   tiny fraction of the library, so any library-level statistic reads it as uncaptured.
3. No probe file; capture is learned.
4. No unexplained constants: every fixed number is derived and documented.
5. One mechanism per A/B. Prototype outside the tree first; changes to the native (C++) solver in a separate worktree.
6. Already refused, with measurements: a training guard admitting only objects with at least one gDNA position (a
   threshold proxy that leaves over-calls at one to three positions); an unconditional exon/non-exon split of the
   landscape (loses every capture-off stratum); a reference read from a strand-only evidence fit (reads no capture on
   unstranded captured libraries).

## 3. What was measured

### 3.1 The damage enters through the ruler, from calibration's per-object errors

Calibration's library-level gDNA on the defect row is right (5.01 M against a truth of 5.00 M). Setting every capture
weight to 1 reads the row at **3.53 %** (genes 0.23 %). The shipped ruler loses it because about 400 short exons, each
with under one gDNA start position, carry a fraction of a fragment of strand noise each. Divided by that opportunity,
their density reads about 200 times the library's. They train a false mode in the landscape (1.9 % of its mass above
30 times the off-target density), the mode is read as the capture reference, and every length contracts unevenly.

### 3.2 Calibration's per-object counts are wrong under a gap, in both directions

Against the oracle, counting objects whose gDNA count is off by more than threefold and more than five Poisson standard
deviations (a reporting criterion for this diagnosis, not a model constant):

| class | capture-off row | objects | signature |
|---|---|---|---|
| short-exon over-call | RNA-short, stranded | 173 regions, +2.4k fragments | median opportunity 1.4 positions; solver density median 101× off-target; ψ variance small (0.17 nats²); fed by the landscape's false mode |
| intron under-call | RNA-short, stranded | 156 regions, −16k fragments | 151 are single-strand introns; median opportunity 1,524; ψ books 0.1 % of their fragments as gDNA, truth about 40 %; confident (variance 0.11) |
| neighbourhood over-call | RNA-long, unstranded | 1,431 boundaries and 136 regions, +86k fragments | ψ books 78 % of a boundary's fragments as gDNA at 2–20× off-target, truth near 1×; confident (variance 0.009); the landscape is unimodal |

The neighbourhood over-call is coherent: whole RNA-rich neighbourhoods are booked as gDNA together, an exon and both
its boundaries. One exon receives 2,071 gDNA fragments of its 3,037; its truth is 47.

### 3.3 The shipped ruler's `None` hides all three

On an uncaptured library the shipped ruler finds no enriched mode and ignores every per-object count. That is why the
current tree scores well on every capture-off row except the RNA-short one, where the short-exon over-call builds a
false mode. The spectrum ruling forbids the `None`, so the next ruler must face the per-object counts directly.

### 3.4 A ruler with no detector, prototyped

Five steps, no constant chosen:

1. Each object's level is its posterior median under a nonparametric maximum-likelihood population of gDNA densities,
   fitted on the landscape's grid with each object's Poisson likelihood. Each object votes by its gDNA opportunity, so
   an object gDNA cannot sit in carries no weight.
2. A floor at the off-target density, the intergenic regions' pooled count over their pooled opportunity, because
   capture only enriches.
3. Coherence: an exonic region takes the opportunity-weighted median of its own level and its two boundaries'; a
   non-exonic region keeps its own, held below its boundaries' larger level where its own opportunity is the smaller;
   every boundary lies between its two flanks.
4. Weights are relative: no reference, no clip, no `None`.
5. A junction is priced by the shipped sum by conservation of bases, capped at the larger of its two pieces' levels.

| transcripts / genes, % | 0.7.1 | current tree | prototype |
|---|---|---|---|
| RNA-short · stranded · off (the defect) | 6.31 / 0.44 | 42.30 / 3.28 | **3.55 / 0.23** |
| RNA-short · stranded · on | 20.33 / 3.52 | 19.97 / 1.53 | 23.74 / 1.69 |
| RNA-long · unstranded · off | 5.29 / 0.59 | 2.13 / 0.22 | **11.71 / 0.91** |
| RNA-long · stranded · on | 7.38 / 1.62 | 4.48 / 0.68 | 7.09 / 0.74 |
| test chromosome · stranded · on | 8.58 / 2.45 | 7.19 / 0.96 | 7.11 / 1.18 |
| test chromosome · stranded · off | 8.20 | 7.39 / 1.08 | 7.52 / 1.08 |
| test chromosome, junction-spanning probes · stranded · on | 40.01 / 1.46 | 32.12 / 1.15 | 37.36 / 0.82 |
| test chromosome, one probe centred per exon · stranded · on | 74.76 / 63.06 | 24.39 / 3.15 | 22.36 / 2.85 |

It fixes the defect exactly and, on the sparsest plasma library, finds the enriched levels the shipped mode reader
misses. It cannot ship: the neighbourhood over-call is coherent, so no rule inside the ruler can tell it from capture,
and the RNA-long unstranded row regresses from 2.13 to 11.71 %.

### 3.5 Refuted on the way

Each was run end to end or scored against the oracle, and fails as stated:

- Relative posterior means under the shipped landscape: an exposure bias at low depth (zero-gDNA rows 11–14 % against
  5.7–7.0 %).
- A quasi-likelihood using ψ's own variance as the count's noise: that variance is sharpened by the prior where the
  count is false (the short-exon over-call) and far too wide where messages pinned an accurate count. RNA-short
  stranded · on reads 56 %.
- The arithmetic and geometric posterior means: each follows a far population level held by a handful of objects.
- Pooling a region's count with its boundaries' crossings: biased at every capture edge.
- A contamination model, each count either right or uniform on zero-to-total: it calls 85 % of captured objects wrong.
- A global junction cap at the captured population's median level: worse than the local cap.

## 4. What the fix requires

- **Honest counts. Calibration's per-object gDNA counts must be right under a fragment-length gap, in both directions and every
  library type.** Gate: the §3.2 census near zero on both gap arms, with the ladder and the test chromosome no worse.
  This is the root cause. No detector-free ruler can hide coherent errors, and with honest counts the shipped `None` is no
  longer what keeps the capture-off rows right.
- **A ruler with no detector**, built on those counts: §3.4's design, with its open items (§5, Phase 3).
- **The release bar.** At or above 0.7.1 on every stratum of every panel; within noise of the current tree on
  capture-off strata.
- **Out of scope here.** The RNA-short stranded · on row is mostly the coupling between capture and fragment length.
  A 78 bp RNA fragment overlaps less of a probe than a 250 bp gDNA fragment, so no gDNA-read level is RNA's level. That
  needs a per-placement capture model, planned separately.

## 5. The implementation plan

Every phase changes one mechanism, has a falsification test that fails on the unfixed code and fires when the fix is
broken, and passes the census and the panels before the next begins.

**Phase 0: the measurement loop, cheap by construction.** The 2026-10-05 session ran end-to-end arms with every stage
on one thread, one condition at a time, from the BAM. That cost 18 minutes per ladder condition per arm, against 6.4
minutes on the standard protocol (scan pinned, every other stage on all cores, conditions sharded across processes).
The loop from here:

| step | cost | decides |
|---|---|---|
| per-object census from a calibration dump against `slot_truth.npz` | about 2 min per genome-scale condition, seconds per test condition | whether a calibration change removes the three error classes |
| offline ruler screen on the same dumps | seconds | whether a ruler variant prices objects right |
| test chromosome end to end, sharded | minutes | every stratum, three probe layouts |
| gap arms and ladder end to end, standard protocol, sharded | tens of minutes | the release strata |
| real libraries, one whole-genome process at a time | about 15 min each | plasma and VCaP sanity |

Nothing reaches the expensive rows until the cheap ones pass.

**Phase 1: diagnose the intron under-call and the neighbourhood over-call to their mechanism.** The short-exon
over-call's mechanism is known (§3.1). The other two are confident wrong answers, so evidence from outside the object
must be pinning them. Two hypotheses lead:

- **For the intron under-call:** an RNA level bound delivered by the message layer (its "level lanes" are lower
  bounds on RNA) over-states the nascent RNA in a stranded intron and squeezes out gDNA the strand data support.
- **For the neighbourhood over-call:** with no strand channel, a count-frame conversion somewhere in the message rows
  or the intron background uses the wrong length law's opportunity. At these boundaries RNA's crossing opportunity is about 3.2 times gDNA's
  (247 against 77 positions). Alternatively, ψ adds the gDNA prior twice: the reference's gDNA term and the fitted
  landscape, already an open defect, tilting gDNA up by half a nat per nat of log density wherever no data constrain it.

Decisive measurement: for each such object, ψ's posterior with each received message row removed in turn, and with the
duplicated prior removed, read against the oracle count. The row whose removal restores the truth is the mechanism.
This needs the solver's per-slot debug capture on one gap condition; on the suite genome that peaked at 7.2 GB, inside
budget (it must never run on the whole human genome, where it reached 25 GB).

**Phase 2: fix calibration, one mechanism at a time.**

- **2a. The landscape's training votes by its information about density.** Each training object's vote is
  `min(1, ρ_off · S)`: the gDNA fragments it expects at the off-target density, capped at one object. A slot that
  cannot hold one fragment at the off-target density carries a fraction of a vote, continuously, with no admission
  threshold. This is Python (the landscape fit); it targets the short-exon over-call. Its risk is the deferred stratum: unstranded captured
  exons take their composition from the landscape, and an enriched mode with less mass could drain them, as two earlier
  arms did. The test chromosome's unstranded · on rows and the ladder's deferred rows are its watch.
- **2b. The fixes for the other two classes**, as Phase 1 finds them. If a message row converts with the wrong
  opportunity, the fix is to convert in the density frame. If it is the duplicated prior, it is the existing prototype that lets the landscape
  replace the reference's gDNA term; that prototype needs a partner keeping probed exons' enriched prior under
  capture, which is an open owner ruling. Native-kernel changes are prototyped in a worktree.
- **Gate for each:** the census on both gap arms, then every panel per stratum.

**Phase 3: the ruler.** Land §3.4's design once the counts are honest, replacing the located-mode reader, the reference, the
clip and the `None` path. Three items are open:

- **The boundary clamp.** Where a probe sits at an exon's end, the boundary there is legitimately more captured than the
  exon's average level, and the clamp under-prices it. It costs the junction-spanning layout: partly probed isoforms
  spread 4 times more within a gene. Its only purpose was one false boundary on a zero-gDNA row, which honest counts
  should remove.
- **The exon assumption.** Probes sit on exons, so an intron's contained fragments are never captured. Without it, short
  introns between probed exons read as captured.
- **Speed.** The population fit is an EM over grid weights; 99 % of the medians are final at 200 iterations, and the
  prototype runs 3,000 unaccelerated ones (2–8 minutes per genome-scale condition, 13–15 on a real library). An
  accelerated solver and a convergence test on the medians are required before it lands.

**Phase 4: the release A/B.** Every panel and stratum against 0.7.1 and the current tree, the zero-gDNA controls, and
the real libraries, VCaP against its read-name truth (gDNA fraction 0.2518).

## 6. Risks

- Phase 1 may find the two unexplained classes are not one mechanism each, which would mean more phases.
- 2a may trade the short-exon over-call for the deferred stratum. Its unconditional relatives failed that way.
- With honest counts, the detector-free ruler may still trail the current tree by 0.05–0.3 points on capture-off strata, since
  its weights are estimated rather than switched off. Constraint 2 makes that the price; it still beats 0.7.1 there.
- The coupling remains: the RNA-short stranded · on row stays near 0.7.1's 20 % until the per-placement model exists.

## 7. Questions for the reviewer

1. Are honest per-object counts the right place to fix this, or is there a detector-free ruler that tolerates coherent
   per-object errors like the neighbourhood over-call, which we have missed?
2. Is `min(1, ρ_off · S)` a sound training vote? We considered a vote proportional to opportunity (which shrinks the
   enriched mode with the probed share of the genome) and a deconvolved population in place of the kernel sum.
3. For the two unexplained classes, is removing one message row at a time the most discriminating measurement, given that the solver
   iterates to a fixed point with the landscape?
4. The ruler snaps each object to a population level through the posterior median. Is that compatible with "capture
   is a spectrum", given that the levels and their shares are learned and continuous while each object takes one?
5. Can the boundary clamp and the exon assumption be replaced by one rule derived from the probe-binding physics
   (a fragment binds its best single contiguous probe part)?
6. Is anything in §3's evidence weaker than we treat it?
