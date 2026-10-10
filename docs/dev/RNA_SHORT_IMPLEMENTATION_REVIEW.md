# RNA-short repair: implementation record and request for external review

**Prepared 2026-10-06. This is a review snapshot in the development sandbox, not a new ruling.**
The owner is requesting external review to improve the next implementation plan before execution.
No new likelihood-profile, prior-partner or capture-opportunity implementation is authorized by this document.

**Review incorporated, 2026-10-06.** Read the [independent review](RNA_SHORT_INDEPENDENT_REVIEW.md)
and [revised plan and response](RNA_SHORT_FIX_PLAN.md#10-ruler-checkpoint-and-proposed-revision--2026-10-06).
The law-frame wording and g98 attribution have been corrected. The discrepancy-center change needs
a separate A/B, and the existing tests do not validate the retained uncertainty model. The proposed
order is now population-first, with no separate population fitted inside the ruler. The review's
finite-count profile remains a model proposal; it is not exactly today's conditional ψ.

The implemented change repairs calibration's composition-message arithmetic under unequal RNA and gDNA
fragment lengths. It is present in the working tree, with tests and documentation, but **not committed or
pushed**. The detector-free ruler and the broader capture-aware opportunity model remain unimplemented.

Review these two questions separately: whether the count-map repair is correct and sufficiently tested;
and whether the proposed next design is sound enough to turn into an implementation plan. A successful
count repair does not establish the correctness of that proposal.

## 1. Review map and exact scope

Start with this document's implementation sections, then read the
[big-picture proposal below](#7-big-picture-proposal-for-capture-aware-opportunities) and
[RNA_SHORT_FIX_PLAN.md §10](RNA_SHORT_FIX_PLAN.md#10-ruler-checkpoint-and-proposed-revision--2026-10-06).
That section is the proposal to **replace point-count inputs with per-object likelihood profiles and
develop the partner that applies the density prior once**. It is a proposed contract, not a completed
statistical derivation or a working candidate.

The implementation baseline is `main` at `a7b03103dda26941aae7cfa0ce2c65eebb99c133`, initially clean.
In the tables below, **current** always means that baseline, not the now-modified working tree.
The release comparator is the receipted 0.7.1 release. Review the working-tree diff; `git show HEAD`
alone will not show these changes. The new test file is untracked and must be included in the review.

| File or artifact | Exact change |
|---|---|
| [transfer_rows.h](../../src/rigel/native/transfer_rows.h), `splice_out_row` | Reverse splice profile now receives both RNA opportunities and the measured central splice rate; gDNA opportunity is no longer substituted for RNA's. |
| [transfer_kernel.h](../../src/rigel/native/transfer_kernel.h) | Corrected licensed splice, intron–boundary, alternative-splice and outside-terminus composition maps; added the shared-density coordinate map; expanded face tables; corrected the alternative-splice discrepancy center. |
| [solve_kernel.cpp](../../src/rigel/native/solve_kernel.cpp) | Wired three additional face parameters through native arenas, preparation, views, table import/export and the row binding. |
| [messages/transfer.py](../../src/rigel/calibration/messages/transfer.py) | Documentation only: describes FORWARD as preserving density mixture in the recipient's count frame and updates terminus wording. |
| [test_splice_out_opportunities.py](../../tests/calibration/test_splice_out_opportunities.py) | New independent analytic density-mixture gates: 48 parameterized cases. |
| [test_transfer_faces.py](../../tests/calibration/test_transfer_faces.py) | Updated existing expected-profile calculations, composed intron-to-exon expectations, and native row calls to the corrected contract. |
| [_transfer_harness.py](../../tests/calibration/_transfer_harness.py) | Added the three new arrays to empty face tables used by tests. |
| [tests/golden/](../../tests/golden/) | Regenerated 42 artifacts for six scenarios after inspecting their numerical changes and checking selected toy counts against BAM-origin truth. |
| [DESIGN.md](../DESIGN.md), [EQUATIONS.md](../EQUATIONS.md), [ISSUES.md](../ISSUES.md) | Updated message semantics, added **Component opportunities in a composition profile**, and recorded the measurements and remaining failures under the move rule. |
| [CLAUDE.md](../../CLAUDE.md) | Updated the re-derived suite baseline from 2,651 to 2,699. |
| [RNA_SHORT_FIX_PLAN.md](RNA_SHORT_FIX_PLAN.md) | Added independent review, corrected the plan's assumptions/order, recorded the execution checkpoint and proposed the next evidence contract. |

The six golden scenarios are `combo_extreme`, `combo_moderate`, `gdna_heavy`, `gdna_light`,
`nrna_heavy_ss90` and `nrna_moderate_ss90`. Each has transcript, gene and locus tables in Feather and
TSV, plus scalar JSON. The new review brief itself is a later documentation-only addition.

## 2. What was measured before changing the maps

The original RNA-short regression is a stranded capture-OFF library with RNA about 78 bp and gDNA
about 250 bp: transcript error 42.30%, versus 6.31% in 0.7.1. Library totals concealed bad local
gDNA counts, which fed a false enriched landscape mode and then the capture ruler.

Fresh calibrations of validated cached scans reproduced the original census on both capture-OFF
strata of both gap panels. Calibration and oracle per-object total counts agreed within `7.3e-12`.
The diagnostic flag was `abs(k-g) > 5 sqrt(max(g,1))` together with `k/max(g,1) > 3` or `< 1/3`.
It is a reporting cutoff, not a calibrated confidence interval for an inferred count.

Complete-locus replays, with the landscape, opportunities, strand model and lattice held fixed,
reproduced production outputs bit-for-bit. Removing the implicated composition messages restored
much of the truth at inspected objects; removing the gDNA or RNA level lanes did not. This localized
the damaging channels before the maps were tested analytically. It did not prove that every residual
count error has the same cause.

The permanent evidence is in
[ISSUES: calibration-detects-capture-on-a-capture-off-library](../ISSUES.md#calibration-detects-capture-on-a-capture-off-library).
The original plan's proposed landscape weight `min(1, rho_off * S)` was withdrawn as an implementation
instruction: it was not sufficiently derived. No landscape-training change was needed to remove the
false mode in the map experiment.

## 3. Exact arithmetic implemented

### Shared notation and invariant

Let a boundary contain `U` unspliced fragments, with gDNA count fraction `f`. Let `Eg_b, Er_b` be
the boundary's gDNA and RNA opportunities, and `Eg_x, Er_x` the corresponding flank opportunities.
Let `s` be the certified RNA rate added on that flank. The map uses

```text
g = U f / Eg_b
u = U (1-f) / Er_b
lambda_x(f) = log(g Eg_x) - log((u+s) Er_x).
```

Here `lambda` is count log-odds. All four opportunities use the existing **capture-blind geometric
operator**, evaluated on each component's supplied length law. **The laws are in different frames:**
gDNA's is uniform-frame; RNA's spliced law is junction-de-tilted but capture-selected. The earlier
wording wrongly implied that both inputs were capture-blind. No opportunity calculator or
fragment-length estimator was changed; the repair makes additional messages read these inputs.
It does not replace the divisors with capture-aware opportunities or solve their frame mismatch.

The forward direction transports the boundary's likelihood profile through this map. The reverse
direction evaluates the flank's likelihood profile at `lambda_x(f)`, retaining the existing marginal
over splice-rate uncertainty. These are likelihood profiles, not probability densities being changed
between coordinates; no Jacobian was added. The familiar `f_b = f_x (U+S)/U` is only the appropriate
equal-opportunity special case.

### Licensed intron–exon splice face

The forward boundary-to-exon map already read both component opportunities and the route rate.
The reverse exon-to-boundary path did not: it used `Eg_b` and `Eg_x` for RNA as well, and used
`n_s/Eg_b` as its central rate. `splice_out_row` now takes `a_r_b`, `a_r_e` and `splice_rate` explicitly.
Its RNA terms are `U(1-f)/Er_b` and `(u+s) Er_x`; the central `s` is the same route sum as the
forward map, `sum_J J/A_J`.

The existing uncertainty scale `sqrt(trigamma(n_s+0.5) + trigamma(U+0.5))`, marginal nodes, row
normalization and interpolation remain. The rate, rather than `n_s/Eg_b`, is multiplied by each
node's exponential perturbation. Nonpositive component opportunities produce no informative reverse
row. Existing strand and face licences remain in force.

### Intron–boundary shared-density face

An intron and its licensed boundary share component density odds, not generally count odds:

```text
lambda_dst - lambda_src = log(Eg_dst/Eg_src) - log(Er_dst/Er_src).
```

`FacesOut.shared_density` constructs the inverse coordinate map needed to read a source profile at
the destination grid. Equal opportunity ratios use the original identity path. Otherwise FORWARD
now sums the source's own and held composition profiles, interpolates that sum at the shifted
coordinates with the existing endpoint behavior, and normalizes it. No new blur is applied.
Both directions are corrected; the boundary-to-intron strand licence remains. A face with a
nonpositive required component opportunity is not installed.

### Alternative splice site, exon on both sides

Both flank maps now use each component's own opportunities. The continuing flank reads the
contiguous spliced crossing rate `S/Er_b`; the exonic flank reads that rate plus the appropriate
route sum. Both forward and reverse directions receive those same central quantities.

The existing excess-disagreement width also changes its center. It now compares
`lambda_x(f_b, measured s)` with the observed flank log-odds. Previously it compared unmapped
count odds using an equal-opportunity correction. The width remains
`max(0, disagreement^2 - (v_b + v_x + v_ratio))`, with the existing variance terms unchanged.
This is a change to the center used by that rule, not a new derivation of its uncertainty model.
It also acts at equal component opportunities, so it is a second mechanism bundled with the map
repair. Its contribution has not been isolated by an A/B.

### Outside flank of an exon–exon terminus

The outside composition pair now uses `Eg_b, Er_b, Eg_outside, Er_outside`. Its RNA rate is
`S/Er_b`, plus the route rate when the junction's exonic side is the outside flank. Its count
argument for the existing counting width includes the corresponding splice flux.
The inside-terminus **level rule is unchanged**.

### Native plumbing and implementation boundaries

Face parameter arrays expand from six to nine: `a_r_b`, `a_r_x` and `splice_rate` join the existing
`n_u`, `n_s`, `a_b`, `a_x`, `width` and `var`. Both production `solve_one_block` and the
`transfer_prepare`/`transfer_solve` instrument path carry them. The `splice_out_row` native binding
requires the new arguments; test callers were updated. No compatibility branch or new tunable was added.

There are no changes here to the scanner/deposit rule, ψ's prior, landscape training, RNA level lanes,
the capture-efficiency estimator, the EM algorithm or fragment-length estimation. The existing
reference reader, clipping and `None` switch remain in production. The junction cap at one predates
this diff and is not part of this implementation.

## 4. Tests, perturbations and their limits

The new tests construct counts directly as density times admissible contained/crossing placements,
using `(gDNA length, RNA length)` pairs `(250,78)`, `(78,250)` and `(150,150)`. They include zero and
positive certified splice rate where applicable and test the delivered profile's mode against the
known recipient mixture. They do not generate expected values by calling the production face map.

| New test family | Final cases | Failing cases observed before its repair |
|---|---:|---:|
| Reverse licensed splice | 6 | 4 |
| Intron–boundary, both directions | 6 | 4 |
| Alternative splice, four directions | 24 | 18 |
| Outside terminus, both directions | 12 | 8 |

All 48 pass after the cumulative repair. Perturbation scripts then replaced the corrected RNA
opportunities or central rate with the old quantities, or restored the intron identity in prepared
native tables; the corresponding gates failed again. These were deliberate table mutations on the
tested native path, not six separately rebuilt production binaries. The 15 existing face tests also
pass after updating their expected arithmetic and compositions.

Review the limits of these gates. The new analytic fixtures use deep counts, a sharp source profile,
positive opportunities and a positive-strand donor/TSS configuration. They strongly test the map
centers; they are not a complete proof of shallow-count coverage, zero-opportunity behavior, every
orientation, or the uncertainty of a sum of differently normalized junction rates. Existing integration
tests and panels add coverage, but those mathematical questions deserve explicit review.
The fixtures set junction opportunity equal to the boundary's RNA opportunity and exercise only
one route side. Thus a builder substituting `S/Er_b` for the route sum can pass these gates.
Prepared-table perturbations do not establish coverage of those builder choices. The revised
plan calls for unequal opportunities, mirrored builders and shallow-count uncertainty tests.

Validation already performed on the implemented code:

- Re-derived baseline: 2,651 tests; after adding 48 cases: **2,699 passed, zero skipped** in both the
  isolated native worktree and the rebuilt main working tree. One existing quadrature warning remains.
- `ruff check src/ tests/ scripts/`, `ruff format --check src/ tests/` and
  `python scripts/design/preflight.py --full` pass. Final documentation gates also pass.
- A selected fresh cached-scan calibration on the installed main build reproduces every per-object
  gDNA count of the A/B prototype bit-for-bit.
- Golden changes were read before regeneration. Maximum transcript-count change was 0.175 fragments.
  The extreme mixed toy's true locus gDNA is 227: estimate 192.631 → 204.506 improves. The moderate
  toy's truth is 189: 189.138 → 192.609 worsens. These are accuracy changes, not a numerical no-op.

## 5. Measured effects and regressions

Each map family was introduced cumulatively in an isolated native worktree and screened on fresh
calibrations of all four gap-panel capture-OFF rows. Intermediate repairs that still failed the
count screen did not receive expensive transcript A/B runs. The final candidate was compared with
fresh `a7b03103` runs and the receipted 0.7.1 release on both gap panels, the equal-length ladder,
and all three test-chromosome probe layouts. All strata and zero controls were read separately.
The alternative-splice stage included the separate discrepancy-center change described above;
these receipts cannot attribute its effect independently of the maps.

The flagged RNA-short stranded region under-calls fell from 156, missing 16,115 gDNA fragments, to
zero. RNA-long unstranded flagged boundary over-calls fell from 1,431, with 55,105 excess fragments,
to 10, with 117 excess. Residual errors remain. Exact before/after censuses, object-class errors and
the complete six-panel tables belong to the linked ISSUES entry, not to a new permanent ledger here.

Selected release-facing results, **transcripts / genes (%)**:

| Panel and stratum | 0.7.1 | Current `a7b03103` | Implemented count repair |
|---|---:|---:|---:|
| RNA-short · stranded OFF | 6.31 / 0.44 | 42.30 / 3.28 | 3.53 / 0.23 |
| RNA-short · stranded ON | 20.33 / 3.52 | 19.97 / 1.53 | 13.33 / 1.40 |
| RNA-short · unstranded OFF | 6.00 / 0.61 | 3.78 / 0.27 | 3.79 / 0.27 |
| RNA-short · unstranded ON | 29.76 / 6.19 | 19.41 / 2.34 | 17.57 / 2.23 |
| RNA-long · stranded OFF | 5.23 / 0.46 | 2.26 / 0.19 | 2.27 / 0.19 |
| RNA-long · stranded ON | 7.38 / 1.62 | 4.48 / 0.68 | 4.23 / 0.66 |
| RNA-long · unstranded OFF | 5.29 / 0.59 | 2.13 / 0.22 | 2.16 / 0.22 |
| RNA-long · unstranded ON | 28.99 / 5.84 | 6.86 / 0.95 | 6.55 / 0.94 |

Every transcript stratum in the six-panel comparison beats 0.7.1. That statement does not mean
every condition improves, every metric passes, or the release is ready. Important counterexamples:

- Test chromosome `g98 / ss0.70 / capture ON`: 0.7.1 **55.98 / 39.60**, current **46.92 / 19.51**,
  repaired **69.50 / 42.44**. Crossed frozen-input replays under both native implementations reproduce
  the inspected objects' losses according to the refitted landscape alone. This localizes the effect;
  it does not excuse the regression or prove a remedy. The archived strand-only replay does recover
  near-truth counts. Local evidence is being overridden by a deeper landscape valley; it is not
  absent. The landscape's broad band masses agree to three decimals. A perturbation sweep is needed
  to quantify ψ's interpolated-median sensitivity.
- Ladder unstranded capture-ON: **33.21 / 14.25**, **10.42 / 4.33**, **12.40 / 4.47**, respectively.
- Ladder unstranded capture-ON zero-gDNA control: gene error remains 3.218%, versus 3.183% in 0.7.1.

The RNA-short captured improvement to 13.33% also refutes treating about 20% transcript error as an
established capture-by-length floor. Real-library release validation has not been run for this change.
Small OFF transcript losses have near-identical genes and pools, consistent with EM amplification;
this does not prove that each loss is noise. Keep them reported until an affected-condition check
separates numerical sensitivity from systematic isoform error.

## 6. Ruler experiments that were not implemented in production

Three ruler variants were screened outside `src/` on the repaired counts, each on all 30 standard
test-chromosome conditions, including zero controls. They retain the earlier prototype's population
fit, off-target floor, spatial treatment and local junction rule unless stated otherwise.

| Prototype | Change and decisive failure |
|---|---|
| `pwfx_own_med` | Previous opportunity-weighted population and posterior-median ruler. It jumps 100-fold between density atoms 1 and 100 when a posterior tie is crossed, for count perturbations down to ±1e-9. It also fails unstranded OFF: 9.41% transcript error against 9.25% in 0.7.1. |
| `pwfx_own_geo` | Replaces only the posterior readout with `exp(E[log rho])`. Its conditional-readout continuity gates pass, including deliberate restoration of the failing median. Unstranded OFF still reads 9.41%. This is not a proof of continuity of the whole population-fitting pipeline. |
| Geometric readout plus inferred-count uncertainty | Changes only the working count likelihood, using delta-method variance `d=k² Var(log f_g)` and quasi-score `(k-mu)/(mu+d)`. Unstranded OFF improves to 8.81%, but ss0.70 capture-ON worsens to 18.73%, versus 10.28% in 0.7.1 and 6.91% on the current tree. |

The quasi-score's algebraic tests pass and fail under deliberate removal of its uncertainty term.
That does not validate the approximation: its capture-ON A/B fails. All three arms stopped before
gap, ladder or real-library runs. None is a candidate to land as written.

Targeted oracle substitutions changed only the counts read by the geometric ruler, leaving the
calibration output used for EM priors unchanged. On g25 / g50 unstranded OFF, transcript error
9.26 / 11.29% becomes 8.08 / 9.77% with all oracle counts. Correcting only intergenic/gene-edge
counts yields 9.13 / 10.72%; correcting only the other objects yields 8.11 / 9.86%. Thus the
remaining failure is not solely a contaminated off-target floor.

One implicated exon estimates 95.803 gDNA fragments against 55 true. Its reported log-fraction
variance is 0.341: the delta-method deconvolution standard deviation is about 56, while treating
the estimate as an observed Poisson count gives about 10. The ruler makes it 1.878 times background.
This motivates retaining more of the evidence, but **does not prove that likelihood profiles will
solve the problem or that no other count repair is possible**.

The experiment harness was extended to arm every sharded child through environment variables.
The population fitter's likelihood exponentials were cached; on the checked fixture the fit and
geometric weights agree with the original implementation to floating-point precision. That is a
prototype optimization, not a validated whole-genome solver. The old prototype's stride subsampling,
iteration cap and spatial assumptions are not approved production choices.

## 7. Big-picture proposal for capture-aware opportunities

**Recommendation for review, not an implemented fix:** use one data-learned description of capture
to calculate component-specific opportunities, while preserving the uncertainty in the evidence
used to learn that description. A shared scalar measured on gDNA fragments cannot generally stand
for RNA fragments with a different length law or different reach. Fixing its input counts alone
does not fix that transfer of capture between components.

**Revised scope after review:** the placement integral below is a specification, not the next
estimator to implement. First price a per-object scalar's transfer with the defined `B_match` and
`B_transfer1` oracle arms. If it is adequate at the release bar, a scalar under an explicit
within-object uniformity assumption may suffice. Otherwise identify the additional statistics
needed before proposing a richer estimator. Their common scale, support and numerator assumptions
must match the existing arm definitions.

There are two distinct gaps in the current system:

1. **Evidence:** the ruler treats a calibrated gDNA point estimate as another observed count. The
   fitted prior and messages helped produce that estimate; a Poisson likelihood around it can
   invent precision. A mean and variance need not preserve a multimodal profile either.
2. **Opportunity:** the ruler applies gDNA-derived object efficiencies to RNA. Capture can depend on
   fragment length and placement, including reach into a nearby captured area. The score's length
   law, the calibration opportunity and the EM normalizer must use compatible definitions.

The proposed profile interface addresses the first gap. It is not by itself a solution to the second.
The broader design must connect them without counting the same capture effect or prior twice.

### A common measurement model, with each component's own geometry

A useful target equation, which still needs an identifiable estimator, is

```text
E_cap[c,i] = sum_w P0_c(w) sum_{a in A_c(w)} D_i(a) C_theta(a).
```

`P0_c` is the explicitly defined pre-selection length law; `A_c(w)` is the component's allowed
fragment placements at that length, including its template/reach; `C_theta(a)` is capture at a
placement, inferred from data; and `D_i(a)` is the deposit operator for the quantity being predicted.
The notation does not add a calibration population: those remain gDNA, RNA+ and RNA−.

`D_i` must distinguish **count incidence** from **conserved mass**. A fragment can cross several
boundaries and increment each incidence count, while its conserved mass sums to one. Calibration
count likelihoods and EM effective lengths therefore cannot silently interchange these operators.
Use the same learned capture description with the appropriate operator and each component's own
length law, rather than substituting one component's averaged efficiency for another's.

This equation is a specification of what needs to be integrated, not a claim that `C_theta` is
recoverable from today's object totals. Different placement fields can produce the same gDNA
contained/crossing averages and different RNA yields. Review must determine whether current
statistics suffice, whether additional length-resolved statistics are necessary, and what reduced
model is defensible where the detailed field is unidentifiable. No probe BED is available to inference.

The scale of density and capture also needs an explicit convention: rescaling one and inversely
rescaling the other can leave expected counts unchanged. Normalization must be common to the
components and compatible with their likelihoods. It must not become a located reference or a
library-level capture switch.

### A consistent treatment of fragment-length selection

Specify whether every consumed length law is pre-capture or already selected by capture. Applying
capture weighting to a selected law again can count selection twice. Estimating an uncaptured law
from a selected pool is itself an inference problem, not a change of variable that can be assumed
correct. A single global multiplicative response can miss local reach at capture edges.

Do not infer one origin's law by multiplying the other origin's law by a ratio with unsupported
tails. Do not introduce a new length-discrimination channel merely because the two estimated means
differ. Any discrimination attributable solely to a length-law gap must vanish continuously as
that gap closes and must account for uncertainty when an origin is scarce or absent.

The existing message fix uses a capture-blind operator with a capture-selected RNA input law.
Replacing its opportunities with `E_cap` would require re-deriving both sides of its maps and its
RNA route rates in the same frame. A one-sided replacement is not the proposed implementation.

### Binding limits on this proposal

Capture stays continuous and local; scarce captured transcripts may remain a tiny fraction of the
library. No detector, on/off status, library-level test, reference-or-`None` gate, probe BED,
unconditional exon/non-exon landscape split, or strand-only sole capture reference is proposed.
No new junction-topology model is authorized; retain the existing junction-price cap at one.
In particular, the abstract placement equation is not permission to invent probe locations across
splice junctions. The admissible approximation there remains a review question.

For historical context, [FRAGMENT_LENGTH_REVIEW.md §9](FRAGMENT_LENGTH_REVIEW.md#9-toward-the-real-fix)
discusses length-resolved opportunities. It predates the current rulings and measurements. Its BED
option, background-mixture prescription and recommended ordering are **not** the current plan.
[FRAGMENT_LENGTH_POSTMORTEM.md](FRAGMENT_LENGTH_POSTMORTEM.md) also contains superseded BED/detector
proposals; it should not be read as an instruction to implement them. The current proposal is the
constrained direction in this section together with the evidence contract below.

## 8. The likelihood-profile and single-prior proposal

The concrete starting point is
[RNA_SHORT_FIX_PLAN.md §10](RNA_SHORT_FIX_PLAN.md#10-ruler-checkpoint-and-proposed-revision--2026-10-06).
Its central proposed interface is an object evidence profile `L_i(rho)` that retains shape and
multimodality. After defining the base measure, an object posterior would have the form

```text
p[i,j] proportional to L_i(rho_j) pi[j].
```

Here `pi[j]` must have an explicit interpretation as prior mass, including any required grid
quadrature weights. The fitted gDNA-density prior enters once. A point count can remain a reporting
output; it is not re-observed as Poisson data by the ruler.

The first proposed consumer is now the **calibration landscape**, which also turns inferred counts
into Poisson kernels. Investigate one population for ψ and the downstream ruler, rather than a
second population inside the ruler. Separate changes to training weights, the population estimator,
the evidence profiles and the prior application. A combined successful arm would not identify its
effective mechanism. The independent review's finite-count expression is not exactly equivalent
to the current conditional solver; the revised plan gives a one-count counterexample and the
required distinction between a density per log-ρ and a density per λ.

This is a proposed direction for the open partner in
[ISSUES: the-gdna-prior-enters-psi-twice](../ISSUES.md#the-gdna-prior-enters-psi-twice).
The known ψ defect adds the fitted landscape to the gDNA half of the reference measure instead of
replacing that half. Its earlier isolated correction helped zero controls and harmed captured rows;
it was not included in the current count-map repair. A profile interface is **not yet a demonstrated
partner** for that correction. Review should reject that claim if the necessary derivation fails.

The implementation plan needs answers to these specific problems before coding:

| Required derivation or decision | Why it matters |
|---|---|
| Define exactly which factors belong to `L_i` and its reference measure. | Removing a stored landscape array from a posterior does not automatically recover a likelihood. Account for the gDNA reference half, RNA marginalization, factory terms and coordinate measure. |
| Separate explicit prior use from hidden prior influence. | Message widths and other inputs can depend on prior-conditioned beliefs. First test cancellation with the other state frozen, then determine what must change when messages are rebuilt. |
| Account for reused fragments and neighbour evidence. | Adjacent objects' profiles can be correlated. Multiplying them or counting them as independent training observations can duplicate evidence even if the explicit prior appears once. |
| Derive the population fit and its weighting. | The previous opportunity vote is not automatically justified for uncertain, overlapping object evidence. Keep one population across annotation classes; establish how a genuinely enriched minority survives. |
| Define the connection from density profiles to capture-aware opportunities. | Profiles alone do not identify a placement-dependent capture field or translate gDNA capture to a different RNA law. State the identifiable quantity and any approximation. |
| Define empty, uninformative, zero-gDNA and unsupported-length limits. | These limits must carry appropriate uncertainty and avoid manufacturing enrichment, ratios from unsupported tails, or a library gate. |
| Specify convergence and bounded storage. | A stable bulk statistic does not bound error at rare captured objects. Use bounded blocks; never retain whole-human debug profiles. |

The known **Jeffreys-shaped landscape reproduces the reference-only ψ** identity is a necessary
prior-replacement test, not sufficient proof of the full model. Tests must also distinguish a
prior-free evidence profile from a posterior with prior influence hidden in its inputs. No source
change implementing this proposal has begun.

## 9. Review deliverable and proposed sequence after feedback

Please return an independent critique, with each finding tied to a function, equation or measured
counterexample. Separate implementation defects, evidence gaps, modelling choices and new owner
rulings. In particular, assess:

1. The direction and inverse of each count-frame map, the RNA route-rate units, and the behavior at
   zero opportunity or zero flux. Identify any untested orientations or compositions.
2. Whether retaining the existing counting widths and alternative-splice variance terms is
   justified after correcting the central maps. Suggest a falsifying example if it is not.
3. Whether the g98 landscape regression or another measured loss changes the verdict on the current
   repair. Stratum averages must not hide individual failures.
4. Whether per-object likelihood profiles are the right next abstraction, and whether they can
   provide a partner for the single-prior repair without a class split or a capture gate.
5. What capture-aware opportunity is actually identifiable from the available observations, and
   the smallest additional measurement or explicit assumption needed if it is not.
6. Which steps should be removed, reordered or separated into distinct falsification experiments.

The [revised sequence](RNA_SHORT_FIX_PLAN.md#10-ruler-checkpoint-and-proposed-revision--2026-10-06)
now starts with the law-frame diagnostic, the isolated discrepancy center and missing map/uncertainty
gates. Population-weighting and estimator comparisons precede profile implementation and the
prior-once partner; the ruler follows that contract. Scalar transfer is priced before a richer
capture model. None of these proposed model experiments was executed while incorporating the review.

There is no need to run new whole-genome benchmarks to review this record. For any future experiment,
check `uptime` and state expected wall time first. Pin only the scan, use fractional EM, retain
default calibration/EM thread budgets and shard independent conditions. Find a cheaper discriminator
if a step would exceed about 30 minutes. The current task ends with this review packet and awaits
the owner's external feedback.

## 10. Evidence package and reproduction references

Local archive: `~/Downloads/rigel_runs/prototypes/2026-10-06_rna_short_counts/`.
The archive is local evidence, not a public download; a reviewer outside this workspace needs the
owner to provide it or the relevant receipts. The formulas, file map and selected results above are
included so the request can be understood without access to that directory.

| Archive item | Purpose |
|---|---|
| `README.txt`, `checkpoint.json`, `receipt_manifest.json` | Run conventions, baseline identity, final status and hashes of the saved evidence. |
| `working_tree.patch` | Frozen uncommitted implementation, tests, goldens and documentation from the preceding implementation session; includes the new test. It predates this review brief. |
| `census.py`, `replay.py`, `messages.py`, `cross_prior.py` | Count census, exact locus replay, channel ablations and crossed-landscape diagnosis. |
| `current_*.jsonl`, `terminus_*.jsonl`, `*_comparison.txt` | Fresh baseline and final count-map receipts, all six panels, with stratum tables and individual conditions. `terminus` denotes the final cumulative map arm, not a terminus-only intervention. |
| `dependencies/release_tag.jsonl` | The receipted 0.7.1 comparison rows. |
| `*_gate_before.log`, `*_gate_after.log`, `perturb_*.log` | Falsification and deliberate perturbation records at each cumulative stage. |
| `golden_review/`, `landed_suite.log`, `landed_preflight.log`, `landed_identity.log` | Golden inspection and validation of the rebuilt working tree. |
| `ruler_{med,geo,uncertain}_test.jsonl`, `ruler_oracle_counts*_test.jsonl` | Failed ruler screens and targeted ruler-only oracle-count substitutions. |
| `uncertain_count_ruler.py`, `test_count_information.py`, `test_ruler_continuity.py`, `fast_population.py`, `harness/` | Non-production ruler experiments, gates and inherited instrumentation hook. |

The archived patch SHA-256 is
`584abf69bdd34166541ab09503a1ae7eaafbef9918ce1900bdbd8ba9dbee37c2`.
The main environment was rebuilt after the native repair. A new run on main therefore measures the
repaired code; it cannot stand in for `a7b03103`. Use the saved baseline receipts or an isolated
checkout of that commit for a new baseline. The archive's README describes the prototype wheel hook.

Permanent sources of settled facts remain [DESIGN.md](../DESIGN.md), [EQUATIONS.md](../EQUATIONS.md)
and the named [ISSUES.md](../ISSUES.md) entries. This snapshot and the proposed plan are provisional
material for external criticism, not approval to continue implementation.
