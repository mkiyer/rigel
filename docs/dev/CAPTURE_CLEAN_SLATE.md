# Learning hybrid capture from the data: a clean-slate design, revised after review

*Sandbox document (`docs/dev/`): provisional, not authoritative, cited by nothing outside the sandbox. First written
2026-10-04 at the owner's direction ("start over, clean slate; the method must learn capture; it needs some gDNA to
learn from; no probe panel"); REVISED 2026-10-05 after an external adversarial review of that version, whose
implementation plan is `CAPTURE_IMPLEMENTATION_PLAN.md` beside this file. The measurements this design rests on are in
`~/Downloads/rigel_runs/prototypes/2026-10-04_clean_slate/` and `2026-10-04_capture_compromise/`; the campaign record
is `FRAGMENT_LENGTH_POSTMORTEM.md`. Nothing here is built.*

**The bar (owner, 2026-10-05).** No probe BED now; capture is learned. The tool must achieve or exceed the shipped
0.7.1 release end to end, per stratum, and it must be robust to a gap between the gDNA and RNA fragment-length laws in
either direction, which is the defect that opened this work.

## 0. What this revision takes from the review, and what it does not

The review found the first version's physics and architecture sound and its estimator unidentified in three places,
and it was right on each. Taken, with the argument in the section named:

1. **Conditioning on exact length removes `ρ_off`, the nascent level and both length laws.** Conditioning each
   boundary's crossing count on the total of boundary and interior placements of the SAME exact length makes the
   uncaptured probability pure geometry; the gap-arm deviations the first version had to explain (0.79 and 1.40, §5)
   vanish by construction. (The review used this as a library-level test; the owner's ruling of 2026-10-05 abandons
   every such test, §6.3, and keeps the conditioning as a device for fitting the kernel's shape.)
2. **Unspliced RNA at a probed boundary is captured too** (the census reads it at 700–1,200× its interior), so the
   crossings alone cannot say whether gDNA is there to learn from; the kernel's likelihood carries the strand evidence
   as a nuisance model for that reason, never as a gate (§6.3).
3. **A probed fraction is not geometry.** Two exons with the same probed fraction and the probe at different
   positions have different RNA yields under a long RNA law (the review's fixture: a 1,000 bp template, gDNA 100 bp,
   RNA 250 bp, one 120 bp part at `[110, 230)` against `[400, 520)`: RNA overlap integrals 20,460 against 30,000, the
   same to every gDNA summary). The map is a one-part INTERVAL per exon context with a retained uncertainty set, read
   from the length- and offset-resolved crossings at the exon's two boundaries, which the first version pooled (§6.4).
4. **The yield is a sum over placements, not a share times an average weight.** The first version's §6.3 step 5 was
   the incidence-average shortcut, which the review's covariance fixture breaks (three 2-base pieces, fragment length
   4, a one-base probe at 0: yield `3 + β`, the shortcut `3 + 0.75β`). Every yield is the direct per-placement
   integral `B_match` already evaluates, under each component's SOURCE law, over the component's routed support (§6.5).
5. **The one new observation is a length-resolved unspliced tally.** The payload keeps per-object counts summed over
   length and five genome-pooled length histograms; nothing in it can condition on length per object or read a probe
   edge from a crossing's offset. Items 1, 3 and 4 all need the same bank (§6.2), built spec-first as every
   accumulator change is.
6. **The census was mis-labelled and the ratio hid its zeros.** Re-run 2026-10-05 (`/private/tmp/rigel_capture_review_audit.py`,
   read-only): 2,175 of the suite's 17,620 projected probe blocks overlap, so the probed fraction summed to above 1 on
   206 exons (max 5.3) and the classes are 799 / 1,744 / 3,166 by union coverage, not 834 / 1,709 / 3,166; and of the
   834 "probed" exons 135 had a positive boundary count over a ZERO intron count, so the medians in §5 are over the
   699 finite ratios. The conditional form of §6.3 has no division and keeps every such exon.
7. **The staging.** The review's phase order — baselines and the corrected census, the evidence bank, the detector
   with quantification untouched, a decision-only end-to-end arm, then geometry, laws, the field, denominators,
   numerators, composition, release — isolates one mechanism per step better than the first version's A/B/C, and the
   decision-only arm is a complete deliverable on its own (§8).

Not taken from the first review, and where the second review met it: the twelve-field risk protocol became a small
decision contract with every unchosen value `null` (adopted: the campaign's `protocol.json`); the mandatory Monte Carlo
p-value became an exact small-case distribution first and a certified scalable tail only where needed (adopted); the
familywise budget is required only where a report promises simultaneous control (adopted). The decision's LIKELIHOOD
is the review's and is settled by derivation (§6.3); the RULE on it is the owner's (§10), and the strand gate's Bayes
factor at equal odds, which the first revision offered as the project's precedent, is not a fixed-size guarantee
(§6.3). The statuses the reviews proposed are abandoned with the decision (§6.3).
Held-out confidence sets by Markov's inequality, cross-fitted priors and certified global optimisation over the
retained set are the review's third level of uncertainty work for the field stage, built only if the finite one-part
enumeration and its profile envelopes prove inadequate on their falsifiers. The review's plan is otherwise the source
of §8 and §9.

**The second review (2026-10-05)** read this revision and answered §10's four questions; the plan adopts its answers
pending the owner's word (§10). It also corrected five statements here, each now fixed in place: the overlap of a
crossing with a probe part is the full interval intersection, not `min(x, q) − p` (§6.4); the junction strand table
already counts each fragment once and lacks only layout stratification (§6.2); a nonlinear hotspot the data cannot
distinguish from capture is a stated limit, not an automatic `unresolved` (§6.3); Phases 0–3 do change rows beyond the
defect's, and which ones is listed and priced (§8); and the strand gate's Bayes factor is a model comparison under its
own prior, not a fixed-size guarantee (§6.3). **Phase 0 began the same day** (§8; the campaign directory
`~/Downloads/rigel_runs/prototypes/2026-10-05_capture_learning/`): the tagged 0.7.1 is built, indexed and validated on
one condition through a neutral adapter that reproduces the instrument to the printed digit, and the pooled-length
diagnostic on 17 cached conditions reads the geometric null within 1 % on every capture-OFF row at genome scale, both
gap arms included, and 10–20× above it on every captured row.

## 1. The question, and the constraints that shape the answer

A hybrid-capture library is enriched for molecules overlapping a probe, by about a thousandfold on the test panels, and
the probes' positions are not an input. Rigel must learn, from the library itself, how much each piece of the genome
was enriched, and then quantify every component (the locus gDNA, every synthetic span, every annotated transcript)
against its CAPTURED opportunity rather than its plain length. Two constraints are given:

- **Capture is learned from gDNA.** gDNA is the one template whose pre-capture abundance is known: uniform across the
  genome. Its post-capture placement density is the capture kernel itself. RNA cannot play that role, because a
  probed exon's RNA count confounds its abundance with its enrichment. A library with no gDNA learns no kernel and is
  quantified at identity weights, and says so (§6.3). That is accepted.
- **No probe panel.** The design treats the panel as something the data must reveal; a panel file could later supply
  the same intervals directly, and nothing here depends on it.

The campaign before this document found one defect and refused four mechanisms (`FRAGMENT_LENGTH_POSTMORTEM.md`). The
defect: on a stranded capture-OFF library whose RNA fragments are shorter than its gDNA fragments, the shipped tree
invents a capture reference from 393 short exons whose gDNA "density" is a fraction of a fragment of strand noise over a
sub-fragment opportunity, and loses 42 % of its transcripts. The design below makes that impossible by construction:
no density is ever read off an object's own sub-fragment opportunity, and no mode is found in anything.

## 2. What 0.7.1 does, what the current tree does, and what the swap cost

**0.7.1 (released 2026-07-12).** The fully captured level is the rightmost significant peak of a MASS-WEIGHTED kernel
density over the per-region gDNA densities (`_global_reference_density`: bandwidth 0.4 in log density, prominence 0.05
of the tallest peak, at least five nodes with mass, the mode snapped to a node density). Each node's efficiency is
`min(ρ_node / ρ_ref, 1)`; every EM component's length is contracted by the efficiencies of the nodes its fragments sit
on, seams included.

**The current tree (since 2026-09-14).** The reference is the located enriched mode of the fitted gDNA landscape, the
basin above the depleted one holding the most LOCATED members (kernels that counted at least one fragment), a mode only
if its members' nearest-neighbour widths resolve it; efficiencies are posterior means under the landscape; every
component's length is one conserved-share rule with a junction priced from its neighbours (`EQUATIONS.md` §11). The
swap was made because the 0.7.1 rule read a mode from strand-noise specks on the zero-gDNA controls (factors 0.51 /
0.12 / 0.62 on `g00` rows) and 0.92 instead of 1.000 on capture-OFF rows.

**What the swap cost.** The located-member census weights a slot by MEMBERSHIP, so a short exon carrying one noise
fragment over 0.16 bp of opportunity counts as much as a probed exon carrying two thousand. Measured on the existing
caches (`mass_fraction.py`, the last refit's training set; "bulk" is the mass-weighted median density):

| condition | gDNA mass within ×3 of the bulk level | located kernels there | mass within ×3 of the shipped reference | located kernels there |
|---|---|---|---|---|
| RNA-short gap arm, stranded × OFF (the defect) | 99.9 % at 0.051/bp | 95.9 % | **0.024 %** at 10.07/bp | 1.9 % (393) |
| ladder g50 stranded × ON | 92.4 % at 1.63/bp | 39.8 % | 92.4 % at 1.56/bp | 39.9 % |
| ladder g50 unstranded × ON | 85.3 % at 1.47/bp | 25.7 % | 85.3 % at 1.50/bp | 25.6 % |
| ladder g05 stranded × ON | 92.8 % at 0.162/bp | 60.4 % | 92.8 % at 0.159/bp | 60.5 % |
| ladder g98 stranded × ON | 92.3 % at 3.20/bp | 35.7 % | 92.3 % at 3.22/bp | 35.7 % |
| test chromosome g50 unstranded × ON | 93.9 % at 1.81/bp | 31.3 % | 93.9 % at 1.80/bp | 31.3 % |
| every capture-OFF row measured | 100 % at one level | 99 % | none | — |

This table is a DIAGNOSIS of the member census, not an estimator: on the defect row the false mode holds a
four-hundredth of a percent of the gDNA mass, and on every captured row the mass sits where the shipped reader also
reads. It does not make "the mass-weighted median" a universal information estimator, and agreement with the shipped
reference is not truth (the review's closing caution; this revision reads no mode at all).

**What has not been measured: 0.7.1 itself on the current panels.** Its mass weighting predicts it would not read the
defect row's false mode (0.024 % of the mass, below a 0.05 prominence); its 8 % capture-OFF contraction is known. Both
are predictions until Phase 0 runs the tagged release on every panel and the real libraries (§8). The owner's bar is
0.7.1's table, so that table comes first.

MEASURED 2026-10-05 on the defect row (RNA-short, stranded × OFF, pinned, receipted): 0.7.1 reads no enriched mode
(`capture.enriched: false`, peak-to-peak fold 1.0 over 20,557 nodes) and transcripts 6.31 % / genes 0.44 % against the
current tree's 42.30 % / 3.28 % with its 10.07/bp reference from 393 members. The prediction held; the rest of 0.7.1's
table followed the same day (the simulated half; the real libraries after it):

| panel | stratum | 0.7.1 tx % / genes % | current tx % / genes % | ahead |
|---|---|---|---|---|
| ladder | stranded × OFF | 3.95 / 0.44 | 2.02 / 0.25 | current |
| ladder | stranded × ON | 9.91 / 2.52 | 6.14 / 1.12 | current |
| ladder | unstranded × OFF | 10.29 / 1.49 | 2.29 / 0.31 | current |
| ladder | unstranded × ON (deferred) | 33.21 / 14.25 | 12.23 / 3.69 | current |
| ladder | g00 false gDNA (str OFF / str ON / unstr OFF / unstr ON) | 41,473 / 25,843 / 401,708 / 21,285 | 127 / 183 / 126 / 335 | current |
| RNA-short gap arm | stranded × OFF (the defect) | 6.31 / 0.44 | **42.30 / 3.28** | **0.7.1** |
| RNA-short gap arm | stranded × ON | 20.33 / 3.52 | 19.82 / 1.52 | current |
| RNA-short gap arm | unstranded × OFF | 6.00 / 0.61 | 3.78 / 0.27 | current |
| RNA-short gap arm | unstranded × ON (deferred) | 29.76 / 6.19 | 18.84 / 2.15 | current |
| RNA-long gap arm | stranded × OFF | 5.23 / 0.46 | 2.26 / 0.19 | current |
| RNA-long gap arm | stranded × ON | 7.38 / 1.62 | 6.17 / 0.62 | current |
| RNA-long gap arm | unstranded × OFF | 5.29 / 0.59 | 2.13 / 0.22 | current |
| RNA-long gap arm | unstranded × ON (deferred) | 28.99 / 5.84 | 8.24 / 0.98 | current |
| test chromosome | stranded × OFF | 8.20 / 1.20 | 7.39 / 1.08 | current |
| test chromosome | stranded × ON | 8.58 / 2.45 | **13.89 / 4.01** | **0.7.1** |
| test chromosome | ss 0.70 × OFF | 9.44 / 1.54 | 8.61 / 1.25 | current |
| test chromosome | ss 0.70 × ON | 10.28 / 3.32 | **15.09 / 5.75** | **0.7.1** |
| test chromosome | unstranded × OFF | 9.25 / 1.76 | 8.78 / 1.33 | current |
| test chromosome | unstranded × ON (deferred) | 18.94 / 10.31 | 18.15 / 7.80 | current |

Read it this way. The current tree is behind the release in exactly two places: the RNA-short defect row, where the
no-reference figure (3.53 %) would beat 0.7.1's 6.31 %, so the decision-only deliverable clears the bar there; and the
test chromosome's captured stranded rows (stranded × ON 13.89 against 8.58, ss 0.70 × ON 15.09 against 10.28), a loss
the ladder does not show (6.14 against 9.91 the other way) and that the decision-only arm, which keeps the shipped
reference on captured rows, cannot touch. Both trees detect capture on those rows, so that loss sits downstream of
detection, in what each tree does with the reference, and it is the one place where the bar as stated (every required
stratum of every panel) is not met by Phases 0–3 alone. The release's own weaknesses are elsewhere: it books 128k to
586k unspliced-RNA fragments as gDNA on the capture-OFF rows of the genome-scale panels and reads 21k to 402k false gDNA
fragments on the ladder's zero controls, which its transcript table does not see; the current tree's zero controls are
two to three orders of magnitude cleaner. The receipted current-tree ladder rows reproduce the 2026-10-03 pinned record
to the digit, so the CLI path and the instrument agree at genome scale as they did on the pilot.

## 3. The criteria, stated once

Three habits of the last month made good rules look like failures and bad rules look like fixes: optimising the
zero-gDNA controls (truth exactly zero; nothing can be "right" there, and the fix for them opened a real regime to a
real defect); reading sub-point transcript moves as verdicts (the 2026-10-04 BED prototype was stopped on a +0.59-point
transcript move at `g00` with gene error and false gDNA unchanged to the digit, which is the EM's start-dependence, not
a falsifier); and requiring every arm to beat shipped on every stratum at once. The criteria to keep:

- **Failure direction first.** A mechanism may be imprecise; it may not invent capture on a capture-OFF library in
  any stratum, and it may not drop capture on a captured one where the evidence is there (the Q3 collapse, deferred
  stratum 12 → 52 %, is a required fallback test though its stratum is deferred).
- **Capture is never switched off.** There is no library-level decision; a weight reads 1 only where the local evidence
  says so (owner, 2026-10-05).
- **The capture-ON ruler is read against the simulator's own yield per probed class** (`ruler_vs_truth.py --scale`),
  class means on one scale and the within-gene spread beside them; genes and pools beside the transcript table; strata
  never pooled; `g00` reported, never optimised.
- **Every mechanism runs on the real libraries as a test input** (the three locally inventoried plasma libraries and the VCaP mix; extend the manifest if another plasma input is supplied) before
  it is called done, and the real data never chooses a rule.

## 4. The physics, as a generative model

A molecule covering genomic bases `B` is captured with weight

    w(B) = 1 + β · max_parts |B ∩ P_k|,

where the `P_k` are the CONTIGUOUS parts of the probe set in the molecule's own coordinate space and `β` the binding
per probed base relative to the off-target weight (`β = 10` on the panels). A molecule hybridises through one
contiguous stretch, so it binds its best single part, never the sum or the union: a probe spanning a junction is one
part in a transcript that holds the junction and two parts in gDNA, in a nascent span and in an isoform without it
(the owner's physics, 2026-09-19; `sim/capture/sampler.py`). The library after capture is one draw from all molecules
with probability proportional to `w`; both origins obey the same rule at the same bases. Three facts carry the design:

1. **gDNA is uniform before capture**, so its placement density after capture, resolved by start and length, is
   `ρ_off · w(z)`: gDNA's own length-resolved placements ARE the capture kernel, up to one scale.
2. **Spliced RNA never crosses an exon|intron boundary contiguously.** The fragments crossing such a boundary
   contiguously are gDNA and RNA that has not spliced there, whatever the library's strandedness. This is
   annotation-conditional: a boundary is eligible only where no annotated spliced path can account for the continuous
   observation (an exon|exon boundary is not eligible), and a retained intron or an unannotated isoform crossing it
   is a nuisance to measure, not to assume away (§7).
3. **Capture is local.** A probe enriches every molecule overlapping it, including the gDNA and unspliced RNA that
   cross the exon's edges into the intron, and it does not touch the intron's interior beyond the probe's reach. A
   probe tiling an exon edge binds the intron bases it overhangs, as real panels do (`ISSUES: intron-seeds-near-probes-are-capture-enriched`:
   near-probe intron seeds read 1.9× the off-target level pooled on `odg05 g05 ON`), so a fitted interval is not
   confined to the exon.

Facts 2 and 3 give an observation of capture that needs no deconvolution: the placements crossing an exon's boundaries
against the placements inside its flanking introns, at the same fragment length. The interior is the local
uncaptured level of exactly the origins that cross; the boundary adds the probe's weight; and a boundary's crossing
opportunity is about one fragment length whatever the exon's size, so the sub-fragment noise floor of the defect
cannot reach it.

What is NOT identifiable without the panel: the capture of a base no gDNA fragment covers (zero gDNA; a probe farther
from every boundary than gDNA's reach); whether two genomic parts join across a junction in transcript space (the
isoform-specific capture of `ISSUES: ruler-witness-geometry-on-transcript-panels`); and, where the RNA law has mass
at lengths the gDNA law has none, the kernel's value at those lengths except through a geometric model (§6.4). Each is
reported as a limit or an interval, never repaired.

## 5. The second measurement: the boundary witness, corrected

For every exon with an intron on both sides, the density of its two boundaries' crossings (count over crossing
opportunity) over the density inside its two flanking introns (contained count over contained opportunity), from the
certified truth and from the RAW totals an unstranded library sees, by the exon's probed class
(`witness_census.py`; median and 10th–90th percentile over the exons with a finite ratio):

| panel · condition | probed, truth gDNA | probed, raw totals | unprobed, raw totals | n probed / unprobed |
|---|---|---|---|---|
| test chromosome g50 ss 0.99 ON | 696 [37–1,040] | 701 [37–1,040] | 0 [0–3.5] | 83 / 117 (of 203) |
| test chromosome g05 ss 0.50 ON | 420 [33–1,630] | 450 | 0 [0–1.7] | |
| test chromosome g98 ss 0.50 ON | 752 | 752 | 0.8 [0–2.5] | |
| ladder g50 ss 0.99 ON | 645 [330–1,180] | 636 | 0 [0–3.7] | 834 / 3,166 (of 5,709) |
| ladder g05 ss 0.50 ON | 400 [132–1,040] | 491 | 0 [0–2.0] | |
| RNA-short gap arm g50 ss 0.99 ON | 682 [354–1,320] | 654 | 0 [0–3.6] | |
| RNA-long gap arm g50 ss 0.99 ON | 363 [247–557] | 384 | 0 [0–3.3] | |
| **every capture-OFF row, every panel, truth gDNA** | **0.97–1.05** [0.54–1.57] | | 0.99–1.00 | |
| capture-OFF, raw totals: ladder and test chromosome | 0.97–1.05 | | 1.00 | |
| capture-OFF, raw totals: RNA-short arm (the defect row) | 0.79 [0.45–1.22] | | 0.85 | |
| capture-OFF, raw totals: RNA-long arm | 1.40 [0.72–2.75] | | 1.29 | |

Corrections from the 2026-10-05 audit (§0 item 6): the suite's classes by union coverage are 799 / 1,744 / 3,166; the
probed-class medians above are over 699 of 834 exons (135 had a positive boundary count over a zero intron count, an
infinite ratio), the unprobed over 2,691 of 3,166 (244 positive over zero, 231 zero over zero). The test chromosome's
203 exons are unaffected (single-block probes). Partially probed exons read in between, 190–470. Read it this way:

- Under capture the witness separates probed from unprobed exons by two to three orders of magnitude in every
  stratum and at every gDNA level from 5 % to 98 %, and the raw totals track the truth: no strand channel and no
  deconvolution is needed. Unspliced RNA crossing a probed boundary is captured too (its own ratio reads 700–1,200),
  which is why the raw total is the right observable: the witness measures the probe, not the origin. It is also why
  the witness alone cannot say whether gDNA is present (§6.3).
- Off capture the gDNA witness reads 1.0 on every panel, including the defect row, whose exons have a median gDNA
  opportunity of 0.4 bp. The raw totals deviate from 1 only on the two gap arms (0.79 and 1.40), through the unspliced
  RNA part: its crossing opportunity scales with RNA's law and the ratio above is in gDNA's units. A raw ratio read
  against 1 is therefore not a decision; a ratio is also the wrong statistic, since it discards every exon whose
  intron holds no fragment. Conditioning on exact length (§6.3) removes both faults.
- **The shadow.** On the test chromosome the deliberately unannotated transcription on `test_blank` lands in the
  "unprobed, RNA-free" class: at `g00 ss 0.50 OFF` 8,748 of its 8,748 fragments are RNA (the BED prototype's
  `witness_census.jsonl`), 18–21 at `g00 ON`, 7–12 RNA deposits among about 8,700 counts at `g50 ON`. An
  annotation-conditional "pure gDNA" class is exactly that, and the design treats it so (§6.3 admission, §7).

## 6. The design

### 6.1 Four principles

1. **Capture is a spectrum (owner, 2026-10-05).** Every capture quantity is per placement and continuous; "no
   enrichment" is the limit where every estimated weight is 1 because the local evidence says so. No library-level
   decision, no detector, no status gates a correction.
2. **One new observation, and the right one.** A length-resolved tally of determinate unspliced placements, built
   spec-first; the kernel and the yields read it, and nothing else in the payload changes.
3. **Capture enters per placement.** Every component's yield is the integral of the learned kernel over its own
   placements under its own source law; a share times an average is never formed.
4. **Saturation is physics.** A fragment binds its best single probe part, so capture saturates at the probe's length;
   no rule may add the capture of two objects as if capture were additive over a fragment's bases (§10: the test
   chromosome's regression is exactly that error).

### 6.2 The one new observation: the length-resolved unspliced tally

The accumulator gains one bank, `unspliced_placements`: for every deposited fragment whose one surviving path is the
unspliced one (determinate, non-chimeric, a single genome strand, length within the limit, exactly as `deposit`
already decides), the row `(ref_id, start, end, align_strand, read_layout, acceptance)` with an integer multiplicity (the layout and
acceptance ids stay in the key until a single-class restriction is proved and recorded: a key without them cannot
carry layout-specific exposures), canonically sorted and
reduced so the export is identical at any worker count (the deferred queue's pattern). `end − start` is the molecular
length. What is excluded is counted by reason (explicit splice, deferred path, strand undefined, too long, empty), and
two exact identities sit beside the existing ledger invariants, taken at the one hook every fragment passes
(`Accumulator::deposit`, where `deposited + deferred + dropped_* == offered` already holds): every offered fragment
lands in exactly one of retained or an excluded reason, `retained + Σ excluded == offered`, and on the accepted-unspliced
domain `retained + excluded_accepted_unspliced == deposited unspliced`; explicit splices, deferred paths and rejected
fragments are never added to that second side. The
bank must not change when the second pass drains a deferred gap: a drained path is a model choice, not raw evidence,
and adding it later means marginalising its hypotheses. Size: at most one row per unspliced fragment, spilled in sorted chunks as
the fragment buffer is, with the budget derived from the finalised dtypes and measured (array bytes, the sort and merge
peak, the spill volume, the cache size), never assumed; whether the production bank keeps only placements within a
fragment length of an annotated gene span (the witness set) or every placement is decided by the observability census
of Phase 1, not here.

Beside it, one small bank the kernel's strand likelihood needs: `spliced_admission`, the certified-spliced fragments
counted ONCE each (not once per junction) as sense or opposite to their transcript's strand, by the layout strata the
RNA strand-error model conditions on. The existing junction strand table already counts each qualified fragment once, at its
leftmost annotated junction, and lacks only the layout stratification, which is what the new bank adds (its raw
marginal serves a single-stratum first test where eligibility and error rate are common); the fitted strand model is
not an observation.

The chain: the Python specification (`tests/native/_accumulator_reference.py`) first, then the native accumulator
byte-identical to it, the payload schema, the spill and reload, the cache digest (`payload_schema_digest` and
`deposit_digest` move; incompatible caches rescan into separate destinations, never over a baseline), and the parity,
order-independence and conservation gates extended to the two banks. The region, boundary and sj banks and every
existing field are unchanged, and that is a gate.

**Cells** (the kernel's observation geometry; under the ruling they feed the field fit, never a test). A pure geometry function maps any hypothetical placement to exactly one of `boundary` (it crosses at least
one eligible exon|intron boundary; ownership by the smallest `(ref, coordinate)` among those it crosses, a partition
rule that never looks at a count), `interior` (wholly inside an eligible intron control), or `excluded(reason)`. The
start axis is partitioned into STRUCTURAL WINDOWS: maximal runs over which the set of continuous RNA templates able to
contain the whole placement is the same (strand and finite reach included, so overlapping genes and gene edges split
windows). A cell is `(window, exact length, align strand, layout class)`; it holds the boundary count, the interior
count and the two EXPOSURES, the exact numbers of admissible starts of that length in that window that the ownership
rule files under the boundary and under the interior, enumerated by a slow exact enumerator first and by an interval
sweep second, equal to the last start. Zero-count cells are enumerated too. A cell with one of the two exposures zero
carries no contrast and says so; no pseudocount, no minimum count, no minimum opportunity. The existing boundary
divisor is not reused: at contiguous boundaries it takes RNA's reach as unbounded by ruling, and a window's exposure
is finite by construction.

### 6.3 No decision (owner ruling, 2026-10-05)

There is no capture decision. A library-level test of capture — the conditional boundary statistic as a gate, its
affine envelope, the size-α or Bayes-factor rule on it, a gDNA admission test, the four statuses and the veto arm
built on them — was the 2026-10-05 plan's first deliverable and is ABANDONED by the owner's ruling: capture is a
spectrum, never on or off. In plasma cfRNA the panel enriches extremely scarce transcripts that remain a tiny
fraction of the library's fragments, so any library-level statistic is dominated by the uncaptured background and
would read "no capture" while the probed loci are enriched a thousandfold. The shipped tree's own reader is such a
decision (`located_enriched_mode` returns a mode or `None`, and `None` switches every correction off); the Phase 0
sweep measured its cost: the two sparse plasma libraries are quantified as uncaptured where 0.7.1 contracts them.

What survives from the section: the length-conditioned crossing geometry (§6.2's cells) is the kernel's
observation; conditioning on a cell's total, under which the off-target level, the unspliced RNA level and both
length laws cancel, is a way to fit the kernel's shape, not a test; the strand evidence enters the kernel's
likelihood as a nuisance model, not as an admission gate. Where gDNA is absent the kernel's likelihood is flat in
the binding amplitude and the weights stay at 1 by estimation; that is the accepted limit, and it is a property of
the fit, not a branch.

### 6.4 The kernel: from gDNA's own placements, geometry where gDNA cannot reach

**Where gDNA is identifiable, no geometry is recovered.** By fact 1 the gDNA placements of a contiguous support `A`
sample `ρ_off · w(z)` resolved by start and length, so the RNA yield on `A` is the gDNA placements re-weighted by the
law ratio, `Σ_{gDNA z ∈ A} src_rna(len z) / src_gdna(len z) / ρ_off`, exact in expectation wherever both laws have
support and acceptance is common (the review's §5.2). It needs to know which placements are gDNA: the truth in the
diagnostic, the raw strand likelihood on one-RNA-strand supports in a stranded library, and nothing in an unstranded
RNA-admitting exon. And it fails exactly on the gap arms, where RNA's law has mass at lengths gDNA's has none (RNA 250
against gDNA 78 bp): the ratio is unidentified there. So the direct transfer is the CONTROL the geometric model is
checked against where both exist, never the production form.

**The production form is one part per exon context, one amplitude per library, both falsifiable.** For each exon
context (the exon, its two eligible boundaries, its intron controls within the window) the model is zero or one
contiguous probe part `[p, q)` with the global `β`, the part free to extend past the exon edge into the intron (fact 3)
and to end outside the observable support (a censored end, not a clipped one). The observations are the
`(length, offset)`-resolved crossings at the two boundaries, which read the part's distance from each edge directly (a
crossing of length `w` with `x` bases inside the exon overlaps `[p, q)` by `max(0, min(x, q) − max(x − w, p))`, the full interval intersection, since a part may extend past
the fragment's start), and, in a
stranded library, the exon's own opposite-strand contained placements. The fit enumerates every integer endpoint pair
plus the null for the context, evaluates each with the physical operator, fits `β` and the RNA nuisance by the raw
count likelihood (zero-count placements in the normaliser; never only the observed starts), and RETAINS the set of
endpoint pairs the data do not distinguish rather than a point estimate. Every requested yield is minimised and
maximised over that set; a wide interval is reported, and an infinite one is `unresolved`. One `β` shared across the
library is the physics (one binding chemistry) and is what lets a sparse library's contexts borrow strength; it is a
model choice tested by the no-full-coverage-anchor and the binding-variation fixtures (§9) and by the per-context
amplitudes' dispersion on the VCaP exome half, whose panel is known. The first version's "`f_e = β_e / β` at the top
of the mass-weighted distribution" is gone: it needed a fully covered anchor and read a fraction where the data read a
position.

**What the kernel cannot see.** Two genomic parts beside a junction may be one part in transcript space or two; no
gDNA placement distinguishes them, and the RNA reads at the junction are refused as a witness (`ISSUES:
a-junction-price-from-its-genes-own-rna`, `ISSUES: the-spliced-read-junction-price`). A junction-spanning RNA placement is
therefore priced as an INTERVAL between the two topologies, its extrema established by the operator (the shipped
additive price is not necessarily either end), a joined part contiguous on the proposed RNA path and consistent across
every isoform it touches; where that
interval is material the component's yield is `unresolved`, which is what the within-gene residual the shipped rule
cannot remove (within-gene sd 0.065 / 0.056 with noise-free counts, truth / sum 0.847 where both sides are probed)
actually is. A probe beyond gDNA's reach from every boundary is invisible to the crossings; in a stranded library the exon's own
contained opposite-strand placements may still see it, and the result names the resolving channel or reports the
coverage unsupported.
At zero gDNA there is no kernel.

### 6.5 The yields

Every EM component's expected count per unit abundance is the direct sum over its admitted placements,

    Y_T  =  Σ_z  src_T(len z) · D_T(z) · w(z),

`src_T` the component's SOURCE length law, `D_T(z)` its routing or deposit share at `z` (one for a transcript's own
placements; the locus's deposit function for the gDNA component), `w` the learned kernel. The conserved-share frame
(`EQUATIONS.md` §11) is this sum at `w ≡ 1` and is kept as the identity check at `β = 0`; it is not the formula under
capture. The gDNA component integrates over the locus's ROUTED support (stage 3's `locus_structure` premises: block
edges on region bounds, nothing shared between loci, no escapable boundary; where a premise fails the union and the
conserved domains differ by 1.3–50 % in the fixture and the row is a sampling-yield diagnostic until the routed domain
is used). Incidence integrals (a molecule counted at several boundaries) stay diagnostics.

**Source laws, not census laws.** The uniform-frame gDNA law exists (`FLModels.gdna_pmf`); an RNA SOURCE law does
not: the RNA pmf is the spliced census, already capture-selected, and integrating it against the kernel selects twice
(`ISSUES: the-scorer-reads-a-census-length-law`). The two-length fixture is mandatory: source `(½, ½)` and weighted
opportunities `(1, 3)` give an observed `(¼, ¾)`, and treating the census as the source yields 2.5 against the correct
2.0. The first yield arm therefore runs with the oracle source laws as a control (stage 3's `laws.py`) beside the
measured laws, so that law error and kernel error are read apart; estimating the RNA source law is its own phase
(§8, Phase 5) and its first form is one shared RNA law fitted jointly with the abundance nuisance from spliced
observations whose selection is enumerated, never from a pooled histogram (which identifies only `src · Σ_t
abundance · opportunity`).

**Numerators.** With capture in the denominators only, the E-step is exact where every candidate placement of a
fragment carries the same weight (a contained exon fragment) and wrong where they differ: an implicit splice against
its unspliced reading, with distinct molecular lengths and overlaps, reads a posterior of 0.870 where the physical one
is 0.904 (stage 3's exhibit). The per-fragment likelihood carries the candidate's log kernel before pruning,
competition and alternative-hit merging, in the ordinary and the multimapper paths alike, and duplicate encodings of
one path are merged while distinct placements are summed. Measured last, as its own stage, with the merger's own
change priced first at capture OFF.

### 6.6 What the EM receives, and what changes nowhere else

The EM receives per-component yields in place of per-component contracted lengths, the same slots in `LocusPriors`
and `effective_lengths_em`, every enabled component of a locus on one kernel, one law contract and one scale (a locus
that cannot be represented coherently stays diagnostic; no mixing of a new yield for some transcripts with old
reference-normalised lengths for their competitors). Calibration's composition (the sweep, the messages, the landscape
as ψ's prior) is untouched until its own stage. Two later, separate arms the design makes possible and does not
depend on:

- **The map in ψ** (Part 1's missing partner, `ISSUES: the-gdna-prior-enters-psi-twice`): fit the landscape on the
  residual coordinate `log ρ_o − log(captured / plain gDNA opportunity of o)` under the learned kernel and shift it
  into each object's solve through a dedicated prior-offset array; off capture every offset is zero and one
  population is fitted, so the refused unconditional class split is not reinstated. First with the oracle map, then
  with the learned one on an independent molecule split, so no observation trains the prior it is then read against.
- **The gDNA count's conserved share** (`ISSUES: the-pooled-q-in-the-gdna-count`): the count converts a boundary's
  gDNA mass by the mixture's `boundary_mass_per_crossing`; gDNA's own expected share under its law and the kernel is
  the appropriate one, audited against per-origin conserved truth first, scored apart from the yields.

### 6.7 What is deleted, and when

After their replacements are complete and scored, never before: `landscape.located_enriched_mode`, `split_basins`,
the basin census and the nearest-neighbour width test; `capture_efficiency.capture_efficiencies`; the junction price
in `capture_eff_length._cut_efficiencies`; `CalibrationResult.gdna_reference_density` and `_members` as the ruler's
inputs, replaced by the amplitude and the retained map. The landscape stays as ψ's composition prior; every landscape
consumer is enumerated before any deletion.

## 7. Failure directions and limits, stated before any run

- **Zero gDNA:** the kernel's likelihood is flat in the amplitude and the weights stay at 1 by estimation. The
  nascent-stress panel's unspliced RNA would witness the probes at `g00 ON`, but real libraries have little or no
  unspliced RNA and the owner's ruling stands: nothing is learned there.
- **Sparse gDNA** (the plasma libraries: 3,434 and 65,874 gDNA fragments as the shipped calibration counts them): each
  probed locus is fitted from its own crossings, and the retained sets are wide and borrow through the shared `β`.
  Failure direction: under-resolution within a gene, reported as wide yields, never a library switched off.
- **Unannotated transcription** booked as gDNA raises the off-target level in the kernel's scale (the shadow; VCaP's
  intergenic +12 %). Failure direction: a mild under-contraction uniform across genes.
- **Nascent gradients and hotspots.** The kernel's background is fitted per window; a nonlinear hotspot at a boundary
  is a falsifier (§9) whose failure is an inflated local weight, measured, not a claim of protection.
- **Real probe physics.** Binding saturates at the probe length (a fragment binds its best single part), is not
  exactly linear in overlap and varies probe to probe; the retained sets widen and the per-context amplitudes
  disperse, measurable on the VCaP exome half against its known panel.
- **Junction topology:** whether the parts beside a junction are one probe or two is not seen by gDNA, and the right
  price differs by up to 65 % between them (§10); the kernel prices both and the choice between them is the owner's.
- **A length gap in either direction:** the kernel's extrapolation across lengths is the one-part model's, checked
  against the direct transfer where the laws' supports overlap and unchecked where they do not; `ruler_vs_truth.py
  --scale` per probed class on both gap arms is the read-out.

## 8. The plan: one mechanism per phase, each with its falsifier and its yardstick

Every phase runs outside the tree first (a prototype directory with a manifest naming the commit, the dirty diff, the
caches and the source hashes), pinned (`--set scan.total_threads=1 --set em.assignment_mode=fractional`), with a
same-session shipped baseline, one real-genome job at a time, and fitting inputs physically separate from the
evaluator's truth (a test proves no oracle reaches the estimator through configuration, cache metadata or read names).
Native work goes in a worktree. Order on every expensive run: the test chromosome, the RNA-short defect row, the other
gap rows, the ladder, the real libraries. The decision-only first deliverable (the old Phases 2–3) is abandoned by the
owner's ruling (§6.3); what replaces it as the 0.8.0 candidate is the owner's call (§10).

| phase | the one mechanism | its falsifier | its yardstick |
|---|---|---|---|
| **0 Audit and baselines** | (a) two arms — `release_tag` (0.7.1 at its tag, its own build, environment and format-7 index) and `current` — the tag run on every panel and the real libraries and scored through one neutral adapter outside the tree, validated on today's outputs against the instrument (DONE 2026-10-05: the full sweep, 54 simulated conditions and four real libraries per arm, §2); (b) the census re-published with union coverage and the zero-denominator counts; (c) the zero-cost length-resolved diagnostic from the existing pools (`kDnaIntronExon` against `kDnaIntronic` per length, genome-pooled: a shape, not a decision; RUN 2026-10-05 on 17 conditions: the pooled boundary count over its constant-intensity expectation 0.98–1.00 on every capture-OFF row of every panel, the RNA-short defect row 0.99 and the RNA-long arm 0.99 included, 10–20 on every captured row and 11.9 on the gDNA-free captured row); (d) fact 2 on real data: the share of the VCaP transcriptome half's unspliced fragments crossing eligible boundaries; (e) shipped calibration, ruler, priors and end-to-end re-recorded; `B_match`'s laws and support verdict recorded | 0.7.1's RNA-short stranded × OFF row: does its mass-weighted KDE invent the mode (predicted not) | the bar, per stratum, for everything after |
| **1 Evidence** | the two banks: the Python specification on synthetic offered events first (schema, ownership, observability and the tiny-case statistics; stop for repair if the fixture window loses its contrast through conditioning), the native bank in a worktree second, the cache digest third (today's `deposit_digest` hashes top-level arrays only, so a nested bank needs canonical recursive hashing or explicit exports); the cell enumerator and the ownership rule; the observability census (molecules and exposures surviving each conditioning dimension; one-sided cells broken out) | native ≠ spec; a deposit or drain changed; a worker-count difference; a fixture window with both exposures losing its contrast through the real ledger | every existing payload field unchanged; `Σ multiplicity + excluded == deposited unspliced` |
| **4 Geometry** | the physical operator (best single part, candidate molecular paths, source-law injection, incidence against conserved integrals), the independent molecule-centric enumerator, the EM's routed support audit per locus, log-kernel and log-yield APIs with exact zero support | the covariance fixture (`3 + β`, not `3 + 0.75β`); two 40-base overlaps contribute 40, not 80; template and reference ends; shared loci | every per-placement and per-object identity exact; tolerances from error analysis, never chosen |
| **5 Laws** | the RNA source law in isolation: the oracle-law control first, then one shared law fitted jointly with the abundance nuisance from enumerated spliced observations; the census-as-source negative control kept | the two-length fixture; unsupported tails stay unsupported; no law borrowed from the simulator in the fit | predicted per-class observed length histograms; uncertainty covering repeated draws |
| **6 Field** | the direct transfer with truth origins, then with the raw strand likelihood; the one-part enumeration with retained sets; yield intervals; the two-part family only if the one-part falsifiers demand it | observationally equivalent fields keep their different yields; alternating-exon targeting borrows nothing; the no-anchor and binding-variation fixtures; the invisible central probe is reported, not imputed | held-out count predictions calibrated; supported yields against the simulator's; unidentified loci listed |
| **7 Yields** | denominators only, composition and count priors frozen: the law-only control, the oracle-map matched yields and the learned-map matched yields, read back off the estimator | frozen-input hashes; a learned/oracle gap attributable to the field or the laws, never to support or scale | `B_match` as the qualified ceiling (ladder stranded × ON −74 % transcripts, −64 % genes); `ruler_vs_truth.py --scale` class means and within-gene spread; the pinned table per stratum, both gap arms beside the ladder |
| **8 Likelihood** | capture in the numerators before pruning, competition and merging, both scorer paths; the merger's own change measured first at capture OFF | the implicit-splice fixture reproduces the physical posterior; duplicate encodings change no probability; unique placements reduce to the single-path formula | the same kernel, law and support hashes reach the scorer and the denominator builder |
| **9 Composition** | the residual-coordinate prior with the oracle map, then learned on an independent split; the prior-accounting fix beside it; the gDNA count's own share; each scored alone, then together | OFF offsets identically zero; no new `g00` pseudo-mass; the refused class split not reinstated | `calibration_vs_oracle.py` and `prior_vs_oracle.py` per stratum |
| **10 Release** | the smallest complete passing scope named as such; port into the layers (`_layers.py` registered; types below consumers; `sweep.py` untouched); §6.7's deletions after every consumer is enumerated; the amplitude and the map's support in results, summary and report; the suite's count re-derived, `preflight.py --full`, fresh pinned panels and the real libraries | no undeclared approximation, no oracle dependency, no mixed scales in a locus | the permanent docs receive what settled, by the move rule; the owner commits |

## 9. The falsifiers, named

Each is a synthetic or panel case with its required read-out; each protection is paired with a deliberately broken
implementation whose failure is recorded. Merged from the first version and the review.

| case | construction | required read-out |
|---|---|---|
| defect row | RNA-short arm, stranded × OFF | weights 1 by estimation, no false reference; transcripts at the no-reference 3.5 %; genes 0.23 % |
| reversed laws | equal laws, RNA 75 / gDNA 250, and the reverse; zero and positive unspliced RNA; capture OFF | the fitted amplitude near 0 independent of the gap; the RNA-long OFF row (raw ratio 1.40) fits no enrichment |
| nonlinear hotspot | an RNA peak at a boundary, high depth | the local weight's inflation measured; no claim of protection |
| shadow and retained intron | the `test_blank` shadow; an unannotated retained intron; readthrough | the off-target level's inflation exposed where unannotated; eligibility moves where annotated |
| finite reach | short introns, gene termini, overlapping same- and opposite-strand genes, single-exon genes | exact window signatures and exposures or explicit exclusion; never unbounded RNA reach |
| repeated boundary | one molecule across several short regions | one event; incidences diagnostic |
| zero denominator | positive over zero and zero over zero cells | retained; no pseudocount, no infinite ratio read as precision |
| junction topology | the benign (separate parts) and the junction-probed (one part across the junction) test panels | each priced right under its own topology; the regression of §10 does not recur on either |
| thread and spill identity | the same scan at 1 / 2 / 4 / 8 workers and every spill size | the two banks and every old field identical |
| coverage union | duplicated and overlapping probe blocks in the evaluator | coverage within `[0, 1]`; part identities kept apart |
| same summary, different yield | the review's `[110, 230)` against `[400, 520)` | summaries agree; RNA yields 20,460 against 30,000; the positional evidence distinguishes them or the yield is an interval |
| part competition | two separate 40-base overlaps against one 80-base part | the best-part response, never a sum or a union |
| conserved covariance | three 2-base pieces, fragment 4, a one-base probe at 0 | `3 + β`; the share-times-average shortcut fails |
| shared and neighbour loci | shared objects, short outside neighbours, multi-block loci, reference ends | each count's routed domain matches its denominator |
| source against census | the two-length selection example | selection applied once; 2.0, not 2.5 |
| structural zero | a transcript shorter than every RNA length; zero counterfactual gDNA support | exact zero or unidentified; no rescue by a pmf floor |
| implicit splice | candidate weights 701 against 501 at distinct lengths | capture in the candidate odds; the posterior matches independent enumeration |
| duplicate paths | distinct alternatives plus duplicate encodings of one | alternatives summed, duplicates not |
| alternating exons | only alternate exons of a gene targeted | no all-exons borrowing; uncertainty follows local evidence |
| no full-coverage anchor | every target partial | no false global amplitude read from the top of a distribution |
| invisible central probe | a probe farther from every boundary than gDNA reaches | unsupported coverage reported |
| binding variation, GC | variable amplitudes, saturation, GC-dependent acceptance, both gaps | errors attributed; no shape assumption read as a fraction |
| sparse plasma | thinnings to the plasma libraries' usable gDNA, and a panel whose captured transcripts are a small minority of the library | per-locus weights recovered where the locus's own gDNA reaches; no library-level reading switches capture off |
| prior and count isolation | the truth map with the current against the residual prior; mixture against origin shares | each change's calibration and pool effect alone |
| oracle isolation | the fit with truth arrays, probe files and simulator parameters inaccessible | identical estimator output |

Every row is read by strand specificity and capture, `g00` apart, the ss 0.70 rows apart, per class as well as in
total: a favourable total cannot hide a class the model is unidentified for.

## 10. The owner's ruling, the regression's root cause, and what is open (2026-10-05, evening)

**Capture is a spectrum (owner).** The detector, its rule, the admission gate, the statuses and the decision-only
deliverable are abandoned (§6.3). Phase 0's measurements stand; the evidence bank, the geometry, the laws, the field,
the yields and the likelihood stay, now as the main line.

**The test chromosome's regression against 0.7.1 is the junction price**, measured end to end on six probe layouts
and recorded with every number in `ISSUES: the-junction-sum-over-prices-separately-probed-exons`. In one line: the
sum `c_lo + c_hi − ½(…)` assumes capture adds over a fragment's bases, a fragment binds its best single probe part,
and on exons tiled to their edges the sum prices junctions 65 % high; clipping it at 1 alone takes stranded × ON from
13.89 to 7.19 % (0.7.1: 8.58 %) and is neutral or better on every other layout except genes on transcript-coordinate
panels. It is principle 4 of §6.1 applied, and it changes the plan in one place: the junction is priced from gDNA's
offset-resolved crossings at its two boundaries — `max(e_A(x), e_B(w − x))` under separate parts, `e_A + e_B` under
one part spanning the junction — which needs §6.2's bank and no probe geometry, and leaves the topology as the one
choice gDNA cannot make (§7).

**The plan, re-ordered under the ruling** (§8's table stands for the phases it keeps):

1. **Phase 0 Audit** — done (§2's tables).
2. **The interim junction price** — one line in `_cut_efficiencies` (the sum clipped at 1), if the owner takes it;
   it reverses the 2026-09-23 junction ruling and `ISSUES: a-junction-price-clipped-at-one`.
3. **Phase 1 Evidence** — the length-and-offset-resolved bank, Python specification first, native second.
4. **The exact junction price** from the offset profiles, A/B'd on the six layouts (benign, sparse and
   junction-probed test panels, the ladder, both gap arms) and the real libraries.
5. **Geometry → Laws → Field → Yields → Likelihood → Composition → Release**, the kernel for every object; the
   reference, its `None` and the clip retire with it.

**For the owner:**

1. Take the interim clip now, or go straight to the exact price behind the bank?
2. When both sides of a junction are captured within a fragment of it, which topology is the default: separate
   probes on each exon (genomic and exome designs) or one probe across the junction (transcript designs)? The
   answer is a fact about the panels your libraries use.
3. What is the 0.8.0 candidate now that the veto is gone: the clip plus the existing reader, or the bank and the
   exact price, with the RNA-short defect repaired when the kernel replaces the reference?

## Appendix: the record

`~/Downloads/rigel_runs/prototypes/2026-10-04_clean_slate/`: `mass_fraction.py` (§2's table), `witness_census.py` and
`suite_probes_bed.py` (§5's table and the projected suite panel), the full outputs per class and row.
`2026-10-04_capture_compromise/`: the BED ceiling prototype (`bed_geometry.py`, brute-force checked; the `g00`
stop; `witness_census.jsonl` with the shadow finding). `2026-10-02_stage3/`: `oracle.py`'s per-placement integrals,
`laws.py`'s oracle source laws, 35 checks with 28 deliberate breaks, `B_match`'s tables. `2026-10-03_fl_arms/VERDICT.md`:
Q3 and Q4. `/private/tmp/rigel_capture_review_audit.py`: the review's census audit (copy it into the Phase 0 artifact
set before relying on it). The 0.7.1 rule: `git show v0.7.1:src/rigel/calibration/capture_eff_length.py`,
`_global_reference_density`. The review's plan, from which §8 and §9 are drawn: `CAPTURE_IMPLEMENTATION_PLAN.md`.
