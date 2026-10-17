# Rigel 0.8.0 — the release plan (2026-10-12)

> **Status, same day.** Stages 1–4 have been executed; the measurements are in `RNA_RELEASE_0_8_0_RESULTS.md`.
> Stage 1 named the held DNA level as the one factor needed at both-strand objects, so the candidate is the split
> keeping only that factor ("split·DNA"). Stage 3's ladder g00 stranded-OFF reading is 3.81 % against the tag's
> 3.70 %, and three failure classes outside that bar are measured with their mechanisms. The owner's decision on
> the listed tradeoffs comes before stage 5.
> One statement below was wrong: unstranded data is not flat under the shipped reader on captured rows. On the
> test chromosome's unstranded × capture-ON rows, main's neighbour factors at single-strand objects do carry
> capture, and the split loses them (results §4).

Sandbox, not authoritative. A plan only: no implementation is authorized by this file, and the owner drives
commits. It replaces the next-step sections of `RNA_RELEASE_REASSESSMENT_VERDICT.md`, whose head lists that
document's corrections. Priorities, in order: robustness across inputs (strandedness, capture strength,
fragment-length laws, depth, probe layouts, real libraries), competitive accuracy, then simplicity.

The new evidence behind this plan is a read-only, unit-aligned screen of archived weights
(`.cache/rigel_runs/2026-10-12_independent_review/split_screen.py`, outputs `split_screen*.out`). It ran no
calibration and no EM; it measures lengths against the simulator's yield, not allocation.

## 1. Settled by the review exchange

- **The shared capture field leaves the release path.** Its archive stays. There is no further sampler,
  transport or joint-optimizer work.
- **Two comparators, two questions.** Tagged 0.7.1 answers whether 0.8.0 is a worthwhile successor to the
  released product. Main answers what a given change improved and what it damaged. No candidate has to
  dominate every improvement accumulated on main, but every loss against main is reported.
- **Fixed constraints.** No detector, no library reference, no new empirical constant, no panel tuning.
- **Tradeoffs are explicit.** A small shortfall against the tag may be accepted when it buys robustness or
  simplicity. It is listed with its size and cause, and never hidden in an average.
- **No automatic follow-on mechanisms.** Each further change is chosen from a demonstrated cause. A
  deterministic likelihood correction comes before any new weighting.

## 2. The defect, by observation model

At an object the reader sees two strand counts. Their expected values are `(G/2 + κR₊ + (1−κ)R₋,
G/2 + (1−κ)R₊ + κR₋)`, with `G` the gDNA and `R±` the RNA on each strand.

- **Single-strand object** (one RNA strand admitted): two equations, two unknowns. It is invertible for every
  `κ ≠ ½`, with information about `G` scaling as `(2κ−1)²`.
- **Both-strand object**: two equations, three unknowns. `G` is exactly confounded with equal sense and
  antisense RNA for **every** `κ`, since `(½, ½) = ½[(κ, 1−κ) + (1−κ, κ)]`. Own evidence bounds `G` from
  above and cannot locate it.
- **Unstranded data** (`κ ≈ ½`): every object is singular. The shipped reader is flat there on capture-OFF rows
  (maximum log weight 5e-5). On captured rows its neighbour factors carry capture, which the split removes at
  single-strand objects (results §4).

| Check | Objects | Shipped | Own evidence | Split (§3), offline |
|---|---|---|---|---|
| Flat-field false capture (ladder g00 str OFF) | single-strand, shallow, gDNA-visible | 696 weights above e; within-gene SD 0.963; tx 13.61 % | 53; 0.403; tx 2.79 % | 102 (49 both-strand); 0.405 |
| Captured both-strand exons (test g98/g50 str ON) | both-strand, long, well covered | antisense factors 84 / 84 / 21 | 0.95–0.98 at g98; at g50 0.88, 0.92 and one at 136; g98 tx 72.21 %, genes 40.83 % | 85 / 85 / 21 (restored) |
| Short captured exons, RNA-short str ON | single-strand, too short to contain gDNA | within-gene SD 0.095; tx 14.24 % | 0.137; tx 20.43 % (tag 20.33 %) | 0.138 (not recovered) |
| Captured rows generally (ladder g05/g50/g98 str ON; RNA-long str ON) | mixed | 0.094 / 0.049 / 0.039 / 0.108 | 0.083 / 0.064 / 0.040 / 0.075 | 0.089 / 0.055 / 0.040 / 0.095 |
| Weak strand (test κ = 0.70 ON) | single-strand | 7.88 % | 10.34 % (tag 10.28 %) | expected ≈ own evidence |

The "Split" column comes from merging the two archived arms with units aligned. Each arm publishes weights
relative to its own median object, and at ladder g05 stranded ON those units differ by 4.88 nat. Unaligned, the
merge reads a false 0.341.

## 3. The candidate: admit neighbour information exactly where the own model is singular

**Rule.** The capture reader uses each object's own strand evidence everywhere. The delivered composition row and
held levels stay in the likelihood only at objects admitting both RNA strands (`statics.free_pos &
statics.free_neg`), where the own model has a structural null direction.

**Why this needs no threshold.** The criterion is annotation eligibility, which already exists. At single-strand
objects the conditioning `(2κ−1)²` is carried continuously by the own likelihood. At unstranded objects the reader
is flat with or without neighbour rows. Nothing is detected, no cutoff is chosen, and no constant is added.

**What it fixes, and what it does not.**
- It removes the delivered lower bounds where they manufacture capture: single-strand objects on flat fields.
- It keeps them where they carry the only information: both-strand captured exons.
- It does not recover short captured exons under a length gap. Pre-registered as a tradeoff: RNA-short stranded
  ON at about tag parity.
- It probably does not recover weak-strand captured rows. Pre-registered the same way.
- It keeps 49 both-strand false weights at ladder g00 stranded OFF. None is among the 93 most extreme, and the
  length screen moves only from 0.403 to 0.405. It must still face every zero control.

**Simplicity.** One branch in the reader block of `solve_kernel.cpp`, plus deletion of the single-strand
RNA-level assembly that no longer has a consumer. If the ablation in stage 3 shows one factor suffices at
both-strand objects, the other factors' reader plumbing is deleted too.

## 4. The release standard — the owner's first decision, before stage 2

1. **Product bar.** Per in-scope stratum (ladder, both gap panels, test chromosome), transcripts and genes:
   candidate ≥ tag, or a shortfall the owner accepts by name.
2. **No hidden condition.** Every condition, in-scope and robustness panels alike, where the candidate is worse
   than both tag and main is listed with its mechanism and an owner decision. The ten largest per-gene overcalls
   of every condition are read beside the averages; that is how the antisense failure surfaced.
3. **Zero controls.** Calibration gDNA at g00 stays at or below the tag's. Capture weights at g00 read flat on the
   ruler: within-gene SD at or below own evidence's.
4. **Calibration.** The six count arrays are identical to main's in every condition, because the reader sits
   downstream of the counts. `P − O` is therefore unchanged and reported.
5. **Real libraries.** VCaP origin pools reported against 0.2518 and against main. A near-gDNA-free real library
   (§5, stage 4) reads near-zero gDNA, with its weight distribution reported. All three plasma libraries run
   cleanly.
6. **Engineering.** Full suite plus the new gates green, `preflight.py --full` green, the golden diff read before
   updating, and VCaP runtime not worse than main's in an idle back-to-back pair.

Numeric tolerances in items 3, 5 and 6 are the owner's to set. None is proposed here as a gate.

## 5. Stages, budgets and stop rules

**Stage 1 — finish the diagnosis (read-only, about 30 min).**
- Done: the unit-aligned split screen above.
- Remaining: attribute the both-strand rescue to its factor (composition row vs held DNA level vs held RNA
  levels). Use the implementation review's factor adapter (`diagnose.py`, `curves.py`) on the test chromosome's
  both-strand captured slots and the ladder's 49 both-strand g00 false weights.
- **Stop** if no factor separates them; then the split keeps all three.

**Stage 2 — build outside main (half a day of engineering, about 5 min build).**
- Build the split arm in a worktree from the review's prototype tree, where `own_field_slot` already exists. If
  stage 1 names one factor, add at most one variant that keeps only that factor at both-strand objects.
- Gates:
  - Calibration arrays identical to main's.
  - Single-strand weights identical to the archived own-evidence arm.
  - Both-strand weights identical to main's.
  - Perturbing delivered rows moves no single-strand weight, and the gate must fire on main.
  - A flat-field fixture (zero-gDNA stranded chain, one false neighbour level) yields weight 1.
  - A captured both-strand fixture keeps its neighbours' level, and the gate must fire under own evidence
    everywhere.
- **Stop** if the identities are not exact.

**Stage 3 — the decisive screen (about 60–75 min of compute; arms: split, the optional variant, and main and own
evidence where not archived).**
- Test chromosome, all specificities, fast: base, junction, sparse, single-probe, binding strength 0.1/1/100,
  g00 seeds 43–45, gDNA overdispersion 0.05, fragment-length rna-long / gdna-long / equal-200, depth d10 / d100
  / low-g. About 260 conditions, roughly 13–15 min per arm.
- Suite: ladder g00 stranded OFF and ON, g05 and g50 stranded ON, RNA-short and RNA-long stranded ON.
- Readouts per condition: transcripts, genes, pools, silent-transcript mass, top-ten gene overcalls.
- **Go** if all of these hold:
  - The antisense overcalls are gone and test g98 stranded ON is at or below main's.
  - Ladder g00 stranded OFF is at or below the tag.
  - No condition is worse than both tag and main without an identified mechanism.
- Use the g00 replicate seeds to say what counts as noise.
- **Stop and diagnose** otherwise. Never append a mechanism by default.

**Stage 4 — full validation of the survivor (about 90–120 min).**
- Ladder (16 conditions) and both gap panels (8); main and own evidence are archived.
- The tag, once, on the robustness panels the October 5 sweep lacks.
- The ladder's low-nascent condition.
- VCaP origin pools.
- A real near-gDNA-free stranded library: the VCaP mix filtered to the transcriptome instrument's reads
  (the RNA library's instrument prefix), name order preserved, about 14 M fragments.
- The three plasma libraries.
- Calibration-array hashes.
- `ruler_vs_truth.py --scale` on capture-ON rows.
- Output: the release table against tag and main, per stratum **and** per condition, then the owner's decision on the listed tradeoffs.

**Stage 5 — production port (about one day).**
- Worktree from the pushed main.
- Make the reader branch, the deletion and the Python reference specification change together.
- Add the stage-2 fixtures as permanent gates.
- Run the suite, read the goldens, run `preflight.py --full`, and an idle VCaP timing pair.
- **Stop** if any number differs from stage 4's receipt.

**Stage 6 — release (about half a day).** By the move rule:
- DESIGN §7.2 takes the reader ruling.
- ISSUES closes `the-capture-weights-are-noisy-on-a-capture-off-library` with its killing numbers, records the
  shared field as REFUSED with its numbers, and opens one entry per pre-registered limitation.
- ROADMAP, then the manual's limitations section.
- The version bump, `PUBLISHING.md`'s procedure, and a manual CI dispatch (CI is dispatch-only).
- Superseded campaign notes in this directory are deleted once their findings have moved. The owner commits
  and tags.

Total: about three working days and four to five hours of compute, if no stop rule fires.

## 6. Not on the 0.8.0 path

Each item below needs a demonstrated cause before it is opened; none blocks the release.
1. The gDNA kernel/normaliser length-frame mismatch: a deterministic correction, tested alone under VB with the
   field and membership held.
2. Short captured exons under length gaps.
3. Weak-strand single-strand objects.
4. Zero-efficiency junctions in the additive junction rule.
5. A non-circular admission of both-strand objects to the landscape (below).
6. A locus-conditional capture-contrast reader, for the residual own-strand false weights.
7. Strand overdispersion (od = 0 by ruling).
8. The deferred stratum's matched-footprint length channel.
9. An optional probe/target BED. It is the only route to identifying capture without gDNA, because at zero gDNA
   the isoform/capture ridge is exact.

## Note — the landscape and both-strand objects

The gDNA landscape is count-native. It is fit first on pass-0 (the prior-free solve), then refit on each
re-solve, three fits in all (`calib_refit_iters = 3`). Each fit sums Poisson kernels of a region's deconvolved
gDNA mass at population (kNN) resolution, weighted by reliability; zero-count intergenic and intronic anchors
carry full weight. The last fit is the composition prior of the final solve and the prior of the capture reader.

Both-strand regions are **excluded from training** (`_fit_gdna_hyperprior`: `sel = expressed & ~(fp & fn)`), as
are boundaries and slots without a located composition. They are **inside the grid's domain**, so the prior is
read at them. The exclusion is principled: their gDNA is the null direction above, so their solved value is the
prior's and the messages' answer, and training on it is circular. It is also measured: admitting them, even
non-circularly, was refused on most conditions and regressed the zero-gDNA guard.

Their one-sided information is real: an upper bound near twice the smaller strand count. A fit could take it
without circularity as a one-sided kernel, as it already takes anchors and location-free kernels; that is item 5.
It would not fix the antisense failure, because a population prior is not spatial. Only the identified gDNA of a
both-strand exon's single-strand neighbours says that this exon is captured, which is why the split keeps
neighbour information there.
