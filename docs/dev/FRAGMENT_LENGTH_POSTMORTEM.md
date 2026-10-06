# The fragment-length campaign: what was built, what was measured, why nothing landed, and the simple path

*Sandbox document (`docs/dev/`): provisional, not authoritative, cited by nothing outside the sandbox. Written
2026-10-04 against `main` at `bdfd8709` plus the uncommitted working tree of 2026-10-03, for external reviewers who have
not seen this project. Every number below was measured in the 2026-10-03 session unless a date says otherwise; the
tables and the harness are in `~/Downloads/rigel_runs/prototypes/2026-10-03_fl_arms/` (`VERDICT.md`). Open questions for
reviewers are collected in §7.*

## 0. What Rigel is, in the terms this document uses

Rigel quantifies RNA-seq transcripts from libraries contaminated with genomic DNA (gDNA). It has three stages:

1. **A scan** tallies fragments onto genomic *objects*: *regions* (an exon, an intron, an intergenic stretch) and
   *boundaries* (the line between two regions, holding the fragments that cross it). Each object's tally is a *slot*.
2. **Calibration** deconvolves every slot into three populations, gDNA, RNA on the + strand and RNA on the − strand.
   Its channels are the *strand* of the reads (a stranded library's RNA reads have one orientation, gDNA's both) and
   the *density* of fragments per unit of *opportunity*. The opportunity of an object for gDNA, `E_g`, is the number of
   positions where a gDNA fragment can sit inside it, which depends on gDNA's fragment-length distribution; a 100 bp
   exon has almost no opportunity for 250 bp fragments. A *landscape* is a population prior over gDNA density
   `ρ = count / E_g`, fitted on the slots and fed back into every slot's solve as the composition prior (the solve is
   called ψ). Messages pass information between neighbouring slots.
3. **A per-locus EM** assigns the RNA to transcripts. It divides each transcript's count by an *effective length*.

**Capture.** Many real libraries are hybrid-capture enriched: probes pull down fragments overlapping targeted exons,
RNA and gDNA alike, by roughly a thousandfold on the test panels. A captured transcript's usable length is only its
probed part, so the EM's effective lengths must be *contracted*. Rigel has no probe input: it **learns** capture from
the data. The learned landscape of gDNA density is read for an *enriched mode* above its main (depleted) basin; the
mode's density is the *capture reference* `ρ_ref`, the fully captured level; each object's *efficiency* is
`E[min(ρ_o/ρ_ref, 1)]` under the landscape; the *ruler* contracts each transcript's length by the efficiencies of
the objects it spans (a junction's efficiency is inferred from the gDNA objects beside it by conservation of bases).
A library with no enriched mode has no reference, every efficiency is 1 and nothing contracts.

**Strata.** Every number is read per stratum: stranded or unstranded library × capture ON or OFF. Three strata are
in scope for the 0.8.0 release; unstranded × capture-ON is *deferred* (reported, never a target). The benchmark panels
are a 16-condition *ladder* (equal gDNA and RNA fragment lengths, by design), a small *test chromosome* (30 conditions)
and two *fragment-length gap arms* on the same genome as the ladder: RNA 250 bp against gDNA 75 bp, and the reverse.
All are simulated with per-fragment truth; the ladder's equal lengths exist so that the EM cannot separate the two
origins by length alone and hide calibration bugs.

## 1. The defect

On the RNA-short gap arm (RNA 78 bp, gDNA 250 bp realised), the **stranded × capture-OFF** library loses its isoforms:
transcript-level error 42.3 % of true RNA fragments against 3.8 % on the unstranded twin of the same library, 80k true
mRNA fragments booked as nascent RNA. The chain, every step re-measured on 2026-10-03:

1. In a stranded library the strand deconvolution credits a short exon (about 100 bp) with a fraction of a gDNA
   fragment: its RNA's few wrong-strand reads, plus the reference prior's own mass.
2. That exon has almost no gDNA opportunity: the slots that make the false mode have a median `E_g` of 0.16 bp of
   admissible starts, 91 % below one position, none at ten.
3. So its gDNA *density* is about 10 per bp, 200× the library's true 0.05.
4. These slots pass every landscape training rule. The rules ask whether the slot's *composition* is located
   (`Var(log f_g) ≤ 1 nat²`), never whether the slot carries information about *density*.
5. They form a located enriched mode: 10.07 per bp from 393 members.
6. The ruler reads that mode as the capture reference, so on a capture-OFF library the regions' efficiencies read a
   median of 0.0054 (10th–90th percentile 0.0048–0.020).
7. Every EM length contracts by the efficiencies of the objects a transcript spans, unevenly across isoforms; the most
   contracted isoforms absorb the shared fragments.

Why the shape: only *stranded* (unstranded exons have no own composition channel and never train), only *RNA shorter
than gDNA* (on the mirror arm gDNA is the short component and short exons do hold its fragments), only *at depth* (at
10 % depth too few slots are located to form a mode). Calibration's own composition metric on the failing row is
unremarkable (region error 1.1 % of mass): the defect is downstream of the composition, in what the ruler reads.

## 2. The root cause, stated once

Calibration reasons in the **count frame**: a slot's gDNA is a count, with a precision in `log f_g`. The landscape and
the ruler need the **density frame**: `log ρ = log f_g + log M − log E_g`. The two frames agree whenever gDNA's and
RNA's opportunities are the same number (`E_g = E_r`), which the equal-length ladder guarantees by construction
(mean `|log(E_g/E_r)|` 0.004) and the gap arms break (up to 6 nats). Under a gap, a *fraction of a fragment* of count
noise, which is harmless in the count frame, is divided by a sub-fragment opportunity and becomes an arbitrary density.

The deconvolved gDNA count of every stranded exon has a noise floor of order one fragment that does not shrink with the
exon's opportunity: the Poisson noise of its wrong-strand reads, and the mass a Beta(½, ½) reference prior puts away
from zero. Over `E_g = 1,000 bp` that floor is 0.001 per bp and invisible; over `E_g = 0.16 bp` it is 6 per bp. The
"located" test (`Var(log f_g) ≤ 1`) is satisfied by such a slot: one fragment known within a factor e *is* a located
composition, and `M/E_g` is a known constant, so the density is "located" too, at the noise floor. No rule in the
training path or the mode reader asks the question that matters: *how much information does this object carry about
the gDNA density?* For a Poisson count that information is the expected count `ρ·E_g`, which for these slots is about
0.01 at any plausible density. The earlier guard, "train only where `E_g ≥ 1` position", was a threshold proxy for this
quantity; the owner refused it as a band-aid (2026-10-02), rightly: over-calls persist at one to three positions.

Everything below follows from this one statement. A mechanism that uses each object's own **likelihood** (the strand
likelihood, with its proper Poisson expectation) weights these slots correctly, because their data are explained by
RNA plus strand error with no gDNA. A mechanism that uses a deconvolved count as if it were a measurement of density
cannot.

## 3. The decisions taken, in order, and what each one measured

**The panels (step 1).** The two gap arms and the junction-probed test panel had been simulated under an older capture
physics and, for the gap arms, an older nascent-RNA model; every number on them was void. They were deleted and
re-simulated on the current simulator with the ladder's nascent block verbatim, so the fragment-length gap is their only
difference from the ladder. Realised lengths: RNA-long 78.58 / 247.63 bp, RNA-short 249.59 / 78.43 bp (capture OFF).
Certified: 4/4 and 4/4 (composition and field); the junction panel 28 field + 2 composition-only rows, the same two rows
as the main test panel.

**The baselines (step 2).** Re-recorded pinned (one scan thread, fractional EM assignment) on the new panels. The
defect reproduces exactly as traced (§1). The RNA-short capture-ON rows read 19.8 % stranded and 18.8 % unstranded,
RNA-long 6.2 % and 8.2 %. The ruler's own instrument shows the second, separate fragment-length defect on the gap arms
under capture: the ruler prices RNA's capture at gDNA's efficiency, so partially probed transcripts read a median log
error of +0.43 to +0.51 on the RNA-short arm and −0.15 to −0.24 on the RNA-long arm, with the sign following the gap,
even when fed the oracle gDNA.

**The simulator deletion (step 3).** The retired `fragment_share` nascent mode and its inert config blocks were
deleted; the suite reads 2,650 passed (2,652 − 2 tests). No calibration or EM source moved.

**Part 1: the prior enters ψ once (prototyped in C++, 2026-10-02; measured end to end 2026-10-03).** The kernel adds the
reference prior's gDNA half *and* the fitted landscape, tilting gDNA by +½ per nat wherever a landscape is fitted. The
fix (the landscape *replaces* the gDNA half) passes its analytic gate to 1e-15 where the shipped kernel fails, cuts
false gDNA on the zero-gDNA controls by 87–96 %, and loses every capture-ON stratum: ladder stranded × ON transcripts
6.14 → 8.98 %, the deferred stratum 12 → 75 %. The extra tilt was propping enriched exons up against a landscape that
one population fit serves for every region class (under capture about 75 % of exon slots are enriched against about
0 % of the others). **Verdict: correct, unlandable alone.**

**Part 2 as ψ's prior (stage 0, 2026-10-02).** A population fit from each region's own likelihood (one shared RNA
amount per region; gDNA-only regions and one-RNA-strand regions). Right on the synthetic falsifiers, right about the
false mode (it reads none), not ready as ψ's prior: its atom aliases on the solver's lattice, the RNA nuisance has no
clean treatment at weak strand specificity (profiled it collapses, conditioned it reads RNA-free regions 13 % low), and
every variant fails at ss 0.70. **Verdict: not as the prior, for 0.8.0.**

**Part 2 as the ruler's reference source only (Q3, 2026-10-03; the owner's scope).** Only `calibrate`'s reference
reading was swapped: the evidence fit, rendered as a landscape, read through the *shipped* mode reader; the sweep, the
shipped landscape in ψ and under the efficiencies, the ruler and the EM untouched (verified: the composition metric is
bit-identical to shipped). Every panel, every stratum:

| stratum | shipped → evidence reference, transcript error % |
|---|---|
| RNA-short, stranded × OFF (the defect) | 42.30 → **3.53** (genes 3.28 → 0.23; mRNA error 79.7k → 0.5k) |
| every other capture-OFF row, every g00 row | identical |
| ladder, stranded × ON | 6.14 → 6.00 |
| RNA-long, stranded × ON | 6.17 → 5.77 |
| test chromosome, stranded × ON | 13.53 → 14.58 |
| ladder, unstranded × ON (deferred) | 12.23 → **52.27** (genes 3.7 → 45.7; gDNA pool −7.7M of 15.3M) |
| test chromosome, unstranded × ON | 18.05 → 44.05 |

The repair is exactly the "no reference" arm's figure: the fit reads no enriched mode on a capture-OFF library. The
collapse: on an **unstranded** library the fit has no evidence class beyond the RNA-free regions (no one-strand class
exists without a strand channel), finds no mode, and declares *no capture where capture is real*; the ruler contracts
nothing and the EM books captured gDNA as RNA. Where the strand exists the evidence reference sits 12–19 % below the
shipped one on the test chromosome and 51 % above it on the RNA-long arm, with small mixed end-to-end effects. Adding
Part 2's step-two class (RNA-admitting regions without strand information, `Poisson(N | X + Y)`) makes it worse: a
reference three orders of magnitude low on one captured unstranded row, and an invented one on a capture-OFF row.
**Verdict: does not land as the one reference source.** It is right wherever it has a likelihood.

**Q4: Part 1 with a class-conditional landscape, four arms (2026-10-03).** The hypothesis from Part 1's verdict: one
landscape serves every region class, so fit ψ's landscape per class (exon slots on the exon-class training set, every
other slot on the rest; two sweeps on one lattice merged by class, exact because the messages never read the prior),
reading the reference and the efficiencies off the shipped solve so that only ψ's prior moves.

| panel · stratum | shipped | Part 1 | class split | both |
|---|---|---|---|---|
| test chromosome, calibration object error, stranded × ON | 1.000 | 1.071 | **0.848** | 0.846 |
| test chromosome, same, unstranded × OFF | 1.000 | 0.990 | **3.719** | 2.470 |
| test chromosome, same, ss 0.70 × OFF | 1.000 | 0.994 | 1.733 | 1.495 |
| ladder, transcripts, stranded × ON | 6.14 % | 8.98 | 6.26 | 6.34 |
| ladder, transcripts, unstranded × ON (deferred) | 12.23 % | 75.02 | 10.55 | 11.50 |
| ladder, g00 stranded × OFF, false gDNA in the EM's pool | 127 | 16 | **15,043** | 16 |
| RNA-short, stranded × OFF | 42.30 % | 35.90 | 45.49 | 37.28 |

The split finds the prize (stranded × ON calibration −15 %, the deferred stratum's transcripts) and fails in scope:
off capture the two classes share one true density, so the split only adds noise; the exon class's landscape is trained
on few located exon slots and, by construction, on no zero-count anchor (anchors are non-exon regions), so exon slots
read a poorly located prior wherever the strand channel is weak, and on the ladder a zero-gDNA control reads 15k false
gDNA fragments. Paired with Part 1 the two failures cancel and the pair ends flat to slightly worse than shipped in
scope. A side finding: a class-wise belief fed to the pooled landscape fit makes it bimodal and the shipped reader then
invents a capture reference on a capture-OFF unstranded row (0.165 per bp from 533 members). **Verdict: none of the four
lands.** Recorded as `ISSUES: the-gdna-prior-enters-psi-twice` (open) and
`ISSUES: an-unconditional-class-conditional-landscape` (refused).

**Stage 3, the ceiling (2026-10-02/03).** With the simulator's *exact* capture-weighted yields as the EM's lengths
(`B_match`), capture-ON transcript error falls 74 % on the ladder (stranded × ON 874,880 → 229,502 misassigned
fragments; genes −64 %) and 39–47 % on the test panels. One efficiency per transcript transferred from gDNA's length
law (`B_transfer1`) equals it at equal lengths and fails under a length gap. This is a ceiling, not a mechanism: it
needs the probe geometry.

## 4. What is structurally missing

1. **Capture geometry.** Everything the ruler needs is a geometric fact about where the probes are. Rigel infers it
   as the mode of a mixture over objects, which is fragile exactly where objects carry little information: short exons,
   long gDNA fragments, unstranded libraries. The simulator knows the geometry; the tool does not.
2. **Information weighting in the density frame.** The training and reference rules test the composition's precision,
   not the object's information about density (`ρ·E_g`). The eff ≥ 1 guard was a threshold proxy; the evidence fit is
   the principled form, and it is right wherever a likelihood exists.
3. **A composition channel for unstranded RNA-admitting objects.** None exists (the plan's §8). In an unstranded
   captured library the evidence for capture lives in the RNA-free neighbours of probed exons, the intron flanks and
   the exon|intron boundaries where enriched gDNA crosses, which the Part 2 prototype excluded (a crossing count books
   one fragment once per boundary). Even with them, a both-stranded exon's own level is imputed, not measured.
4. **The capture × length coupling.** Under capture with unequal lengths, RNA's capture must be priced from gDNA's
   *length-resolved* mass on RNA's own placements, which Rigel does not collect; the shipped ruler prices it at gDNA's
   scalar efficiency (`ISSUES: the-scorer-reads-a-census-length-law`). With geometry it is a computed integral
   (exactly what `B_match` evaluates).

The defect of §1 is a symptom of 1 and 2 together; the Q3 collapse of 1 and 3; the stage-3 gap of 1 and 4.

## 5. The simple path: make capture an input and delete the learned reference

**Proposal.** `rigel quant` takes an optional probe BED. With it, capture is geometry:

- every object's capture weight is its overlap with a contiguous probe part (the simulator's one-contiguous-part
  physics, which the owner ruled on 2026-09-19), scaled by **one scalar learned from the data**, the on/off-target
  enrichment ratio, measured as the gDNA density of probed RNA-free objects against unprobed ones. That is a two-group
  ratio over known members, not a mode found in a mixture;
- every transcript's capture-weighted effective length and the locus gDNA component's are computed exactly from the
  geometry and the two fragment-length laws (the `B_match` integral), so the capture × length coupling is solved by
  construction, and a junction's efficiency is geometry, not an inference from its neighbours;
- the landscape's enriched mode, the capture efficiencies as posterior means and the junction price are **deleted**
  (`located_enriched_mode`, `capture_efficiencies`, the junction rule in `capture_eff_length`; about the size of the
  reference machinery, with its grid-point and within-gene-noise issues);
- without a BED, Rigel assumes no capture and contracts nothing, and says so. A learned, strand-only reference (the
  evidence fit of §3) could remain behind a flag for stranded libraries with an unknown panel, or be dropped.

**What it fixes.** The false-reference defect becomes impossible (there is no reference to invent). The 74 % ceiling of
stage 3 becomes reachable in the part the geometry explains. Part 1 can land, because the landscape no longer has to
carry an unknown enriched class: a probed exon's gDNA prior is the off-target density times its known weight, and the
landscape returns to one population. The within-gene junction-price noise disappears. The method gets smaller.

**What it costs and what is unproven.** Users without a probe BED lose capture-aware lengths; in practice every
hybrid-capture kit publishes its BED and the owner's own captured samples come from known panels. The simulator's
binding model (weight linear in overlap bases, one contiguous part, a single off-target weight) is an assumption on real
data; the learned ratio absorbs the scale, but the *shape* (linearity in overlap, no probe-to-probe variation) is
untested there. Probe-to-probe efficiency variation would show as residual dispersion in probed objects' gDNA density;
whether it matters is measurable on the VCaP mix, whose exome-DNA half carries per-fragment truth and whose panel is
known.

**How to decide it, cheaply, before building.** Three measurements, all outside the tree:
1. `B_geometry`: `B_match` with the simulator's constants replaced by the one learned ratio, on the ladder and both gap
   arms; it should sit near `B_match` where the binding model is the simulator's (it is), which measures what the
   learned scalar costs.
2. The same arm with a *wrong* BED (the sparse-probe or junction-probe design of the test chromosome against a library
   simulated with the standard one) to price the failure mode a user's mislabelled panel would cause.
3. On the VCaP mix with its exome BED: probed against unprobed gDNA density on RNA-free objects, its dispersion, and
   the end-to-end gDNA fraction against the read-name truth.

**The alternative without new input**, if capture must stay learned: declare capture only when the enriched mode is
corroborated by RNA-free objects adjacent to probed-looking exons (the gDNA witness), weighted by their information.
This is a new rule with its own thresholds, it keeps the mixture reader, and it does not touch the capture × length
coupling. It is the path of least change and most residual risk.

## 6. What was NOT the problem

- Not the strand model: the tree before the strand changes gives the same failing row.
- Not the EM's start or its strand term: isolated arms moved the row by under one point.
- Not the simulator's truth: read-name scoring reproduces the tables.
- Not the composition: calibration's object error on the failing row is ordinary. The defect is one scalar, the
  reference, read off a population statement that the objects behind it could not support.

## 7. Questions for reviewers

1. Is the density-frame diagnosis of §2 complete? In particular: is "the object's information about density is its
   expected gDNA count" the right statement, and does it predict every row in §3?
2. Should a reference be read from a mixture mode at all, in any library type, or only from a two-group comparison
   over objects whose membership is known or strongly witnessed?
3. The evidence fit is right wherever a likelihood exists and has none for unstranded RNA-admitting objects. Is there a
   principled likelihood for an unstranded exon's gDNA that does not route through a neighbour, or is §4.3's statement
   final?
4. Does the probe-BED route's binding model (linear in overlap bases, one contiguous part, one learned off-target ratio)
   have known failure modes on real hybrid-capture data that the three measurements of §5 would not catch?
5. Is there a reason to keep the learned reference as a fallback once the BED route exists, given that its only
   trustworthy regime is stranded libraries, or should a library with no BED simply be quantified as uncaptured?

## Appendix: where everything is

- Tables and verdicts: `~/Downloads/rigel_runs/prototypes/2026-10-03_fl_arms/VERDICT.md`; the data in `q3/`, `q4/`,
  `step2/`; the harness `harness/sitecustomize.py` (`RIGEL_ARM=evidence_ref | evidence_ref_blind | classwise`,
  `RIGEL_BUILD=part1`), inert unless armed, one stderr line per process when it fires.
- Part 1: the worktree `~/proj/rigel-part1`, `prototypes/2026-10-02_part1/` (gates, class diagnosis, the patch).
- Part 2: `prototypes/2026-10-02_fl_stage0/` (the fit, the falsifier bench, the verdict).
- Stage 3: `prototypes/2026-10-02_stage3/` (`B_match`, `B_transfer1`, 35 checks with 28 deliberate breaks).
- The campaign's state: `prototypes/HANDOFF_2026-10-03.md`; the plan: `docs/dev/FRAGMENT_LENGTH_PLAN.md`; the five
  reviews that shaped it: `docs/dev/FRAGMENT_LENGTH_REVIEWS.md`.
- The record in the permanent docs: `ISSUES: calibration-detects-capture-on-a-capture-off-library`,
  `ISSUES: the-scorer-reads-a-census-length-law`, `ISSUES: the-gdna-prior-enters-psi-twice`,
  `ISSUES: an-unconditional-class-conditional-landscape`; `docs/TESTING.md` for the panels.
