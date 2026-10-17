# Independent review: paired fragments, intergenic gDNA and splice artifacts

2026-10-15 · Packet: `INTERGENIC_GDNA_DESIGN.md`, Parts A and B · source: `41899ae4`.

Scope: the evidence bundle, current scanner/calibration/scorer, both prototype diffs, and the subsequent splice
derivation and implementation proposals. This is a design review, not authorization to merge either prototype.
The owner's latest requirements govern: complete paired hypotheses, mappable opportunity, design before code,
length correctness blocking release, and STAR validation on the cluster. No production source was changed.

## 1. Verdict

**Do not build arm 2 as specified. Preserve the competing-placement idea and repair its foundation first.**

- **Blocker:** the scanner admits single-mate records as whole fragments; a native reproducer deposits 100-bp lengths from molecules with valid 300-bp pairs.
- **Blocker:** the inferred length distribution amplifies different pool contamination; the bundle's raw global 100-bp fraction is 2.71%, versus 76.33% in inferred gDNA.
- **Blocker:** opportunity omits mapping/acceptance; removing all-N truth reads alone does not fix abundance.
- **Blocker:** arm 2 mixes an absolute external count density with normalized local EM weights, without a derived exposure or count update.
- **Blocker:** summation, capture, candidate deduplication and opportunity are not yet one normalized model.
- **Blocker for the full splice proposal:** its base-only calibration hold rule and universal sample-allele interpretation are not justified.

## 2. Findings

**1 — blocker · §§2, 6–9: a BAM record is being promoted to a fragment.** *Code-derived and measured.*
`bam_scanner.cpp:785` accepts incomplete HI groups; the no-HI branch emits mate-unmapped secondaries and accepts
a primary when either mate exists. `pair_multimapper_reads`, at line 1011, explicitly emits unmatched mates.
`AssembledFragment`, at line 382, retains merged blocks but no mate-completeness certificate. Its outer span then
becomes a fragment length, and `collect_placement`, at line 1621, admits it to calibration.

A six-case native scan/drain experiment used 64 molecules per case. Complete pairs deposited 64 lengths of 300 bp.
A missing primary mate deposited 64 lengths of 100 bp. A valid primary pair plus a lone intergenic secondary
deposited **33 lengths of 100 and 31 of 300**: this occurred with HI grouping, the mate-unmapped secondary branch,
and the unmatched-mate fallback. Supplying both remote mates restored **64 lengths of 300**. The baseline already
has this bug; the intergenic prototype broadens its downstream exposure. See the executable
[`paired_fragment_audit.py`](paired_fragment_review_evidence/paired_fragment_audit.py) and
[`paired_fragment_results.json`](paired_fragment_review_evidence/paired_fragment_results.json).

**Test:** the missing-mate case must produce no ordinary abundance hypothesis or length deposit; adding a lone
secondary to a valid pair must leave its accepted candidates, length histogram and abundance unchanged unless mate
rescue constructs another valid pair. Include both HI/no-HI routes. Break each route deliberately after the fix.

**2 — blocker · §8 length diagnostic: the contrast magnifies the defect.** *Measured from the supplied TSVs.*
At 100 bp the aligned intergenic pool has **8.86%**, the intronic pool **2.01%**, and all deposited fragments
**2.71%**. Inferred gDNA has **76.33%**, mean **99.19 bp**; the oracle has **0.242%**, mean **216.69 bp**.
Thus 76% is not a measured fraction of single-mate molecules. The inversion
`g = ((1-a1) f0 - (1-a0) f1)/(a0-a1)` assigns between-pool differences to gDNA, then clips negative bins and
renormalizes. **The full pinned rerun** gives the contrast's internal density **0.0007978**, not final calibration's
0.0521; its inferred DNA shares are **0.024326** intergenic and **0.010909** intronic. The de-tilted 100-bp bin is
0.088354 versus 0.019654. Substitution into the formula produces **5.0843** at that one bin, with total positive
mass **6.6608** and negative mass **−5.6608**. Clipping/normalizing yields approximately **0.7633**. These are not
probability distributions before projection. Under the assumed two-common-shape model, even the single-bin
difference cannot exceed `abs(a0-a1)=0.013417`; the observed difference is **0.068700**.

The second-pass drain raises the intergenic 100-bp count from **8,218 to 147,420**. The fixture proves incomplete
alternatives can cause such deposits; it does not assign provenance to every one of these library-wide additions.
Differential malformed-placement contamination violates the model, and real selection/opportunity differences can
also violate it after pairing is repaired. Positive opportunity assigned to unmappable sequence is a specific
hypothesis for the internal density failure, not yet an isolated causal measurement. A correct final bulk density
does not validate this separate fit. See the [audit script](paired_fragment_review_evidence/length_contrast_audit.py)
and [complete measured result](paired_fragment_review_evidence/length_contrast.json).
**Test:** trace raw and de-tilted pools, purities, signed contrast, negative mass and final normalization. Repeat after
complete-pair repair; vary read length and genuine insert length independently. Keep real overlapping 100-bp pairs.
Removing bin 100, imposing a read-length floor, smoothing harder, or using a survival likelihood for the lone mate
would conceal the structural error. Quantifying its contribution to the full library requires the repair A/B.

**3 — fix before building · §§2, 9: annotation-first greedy pairing is another source of wrong candidates.**
*Code-derived.* The scanner first pairs mates sharing transcripts, then greedily consumes unmatched locations by
distance, then cross-pairs remaining records. A location used in one pair can disappear from later legitimate
alternatives. A shared transcript or closest coordinate is not a complete molecule validation. Supplementary
records are also merged into primary collections without preserving alternative query chains.
**Test:** exhaustive small-reference paired placement enumeration; permute record order, remove redundant HI tags,
duplicate records, and vary annotation while holding molecular paths fixed. Ordinary DNA candidate geometry must
not change with annotation. Require orientation, compatible query coverage, both mates, and a valid molecular span.

**4 — blocker · §§4, 9: the denominator must describe accepted molecules.** *Derived and code-verified.*
`gdna_opportunity_from_index` constructs opportunities from region and reference lengths, without a mappability
input. The owner is right: zero-observability sequence cannot contribute ordinary counting opportunity. Neither
non-N base count nor uniquely mappable base count is sufficient. A repeated interval can produce a valid pair;
two mappable ends can span an unsequenced N gap. Zero observed coverage is not evidence of unmappability.
The subsequent `SPLICE_ARTIFACTS_PLAN_REVIEW.md` §3.2 proposes requiring the entire insert to be N-free; that is
too restrictive. A 300-bp fragment with unique 100-bp ends and an unsequenced 100-bp N gap has a valid observable
pair, although no N-free segment can contain its entire insert. Its proposed segment formula would give zero.
**Test:** add an unmappable decoy contig to a miniature reference. Existing-fragment density must stay invariant.
Add an observable identical copy instead: its distinct physical starts must enter opportunity, while each observed
pair still contributes total mass one. Enumerate starts, lengths and orientations to verify both cases exactly.

**5 — blocker · §§4, 9: fixed external density does not yet have compatible units.** *Derived from the solver.*
`em_solver.cpp:762` uses `log(theta)-log(effective_length)` with normalized `theta`; VBEM uses the corresponding
normalized digamma weights. Calibration's `rho` is fragments per opportunity. Adding `rho*c*f*e` directly beside
`theta/ell*f*e` mixes count intensity and probability density. Leaving the external alternative out of the M-step
while keeping the old local normalization is not a derived EM update. **Test:** duplicate an entire library and its
calibration evidence. With unchanged proportions and controlled prior strength, responsibilities must be invariant.
A finite reference enumeration must integrate each component density to one and conserve local plus external mass.

**6 — blocker · §§4, 8–9: sum is correct only over distinct valid origins with the matching normalization.**
*Derived; measured intervention.* Mean dilution is real, but the packet incorrectly says sum and mean agree when
one placement dominates. For two terms `1` and `epsilon`, sum tends to 1 and mean to 1/2. Summing malformed,
duplicate, clipped or inconsistently normalized alternatives is not a repair. RNA currently keeps a best hit per
transcript while DNA averages hits; both need the same definition of distinct latent origins, though their physical
support differs. **Test:** duplicate an alignment representation: no change. Add a distinct near-zero-likelihood
origin: no material change. Add a real identical genomic copy: reproduce enumerated probabilities and opportunity.
Repeat at gDNA truth zero. The packet's additional **12,288 false gDNA calls** in that condition remain unexplained.

**7 — fix before building · §§4, 9: capture cancellation and plug-in accuracy are overstated.** *Derived.*
Capture cancels only when the competing hypotheses share the same capture probability. A genomic fragment and a
spliced transcript fragment can have different probe overlap. The present scorer also reads a capture-selected
length law, so multiplying by a new placement capture factor needs a conditional-density derivation to avoid
selecting twice. A fixed density error generally perturbs a responsibility at first order (§4 below).
**Test:** enumerate a small probe design with equal and unequal candidate capture probabilities, equal-length
controls and both length-gap directions. Perturb external density above and below truth, including no-unique-support
repeats; compare plug-in responsibilities and uncertainty with a joint reference calculation.

**8 — blocker for the splice design · implementation §§1, 3 and derivation §7: a sequence tie is not an artifact
decision, and sequence preference is not a safe hold filter.** *Derived.* Both a true splice and a contiguous origin
can explain the same bases. Strand, length, opportunity and abundance can reverse a small base preference. The old
review already supplies such a partial-screen counterexample. Deposit-then-move is bookkeeping, not automatically
independent evidence: the focal row can already have changed junction support, length models or strand parameters.
**Test:** compare every stage against exhaustive full-factor scoring, including wrong-strand, unequal-length and
all-held junctions. Remove the focal contribution before fitting/scoring, or specify a coherent joint fit. Merely
subtracting its deposit after scoring must fail a deliberately self-supported fixture.

**9 — blocker for the splice design · derivation §9: a sample-specific DNA allele is not inferred merely because
reference search failed.** *Derived.* With unstranded data a novel junction count can have expectation
`mu_R + alpha*rho*W`. Many RNA/allele combinations produce the same observation. Reference-spanning coverage helps
only under additional assumptions about RNA contamination, copy number, capture and opportunity. `rho*W` bounds an
expected contribution, not the realized count; low-expression genuine RNA need not exceed it. Annotation site
density is not the distribution of *reported* DNA breakpoints, and annotation reuse does not identify novel RNA's
class distribution or a Beta-spike prior. **Test:** genuine withheld low-depth junctions versus matched DNA deletion
junctions, with identical bases, orientation and coverage in the unstranded arm. An asserted data-driven confident
separation must fail. Require explicit allele geometry, opportunity and an identifiable prior/fit before claiming
this source solved. Failure to find a reference placement is otherwise an unresolved coverage limitation.

## 3. The intergenic candidate: ranked models

| Rank | Model | Assumptions, failure mode and cost |
|---|---|---|
| 1 | **One locus gDNA component, with a normalized spatial law over its complete candidate support** | Preserves settled per-locus EM. Use calibration to shape relative regional rates, not to lend an arbitrary gene rate to a remote copy. Recompute opportunity and prior for the same observation population. Moderate implementation; routing/overlap normalization must be derived first. This is a corrected extension of arm 1, not arm 1 unchanged. |
| 2 | **Fixed regional density** | Sensible when independently informative calibration determines a region's intensity. Arm 2 needs an explicit conditional mixture or Poisson exposure/count model first. Repeats with no independent support have uncertain rates; a global fallback is an assumption, not a measurement. Moderate cost, substantial model-validation burden. |
| 3 | **One genome-wide gDNA component** | Statistically clean common intensity with a calibrated spatial/capture field; one scalar alone assumes homogeneity. Global coupling departs from the settled solver and can be expensive. Useful as a tiny-reference correctness oracle. |
| 4 | **A component per intergenic region, jointly coupled** | Can represent heterogeneous rates, but many repeats provide no information separating region rates. Requires regularization and globally consistent coupling, with substantial state and solver changes. Not justified merely to repair pairing. |
| 5 | **One component per chromosome** | Arbitrary biological boundary; does not fix local capture, copy number or cross-chromosome repeats. Adds coupling without a compelling observation model. |

For rank 1, conditioning matters. Let a deterministic rule route each accepted fragment to one connected transcript
component, with an additional intergenic-only destination. Its DNA opportunity is the expected accepted starts
routed there, including remote starts. Do not independently add the same remote region's entire opportunity and
prior to every linked locus. With shared ambiguous repeats, the route probabilities or equivalent global accounting
must partition the exposure. If this cannot be represented efficiently, the cost of the settled per-locus ruling
must be brought back explicitly; it cannot be hidden in a new fixed term.

## 4. A normalized derivation

Let `h` specify a whole physical fragment: source/path, start, molecular length and orientation. Define
`u_k(h)=f0_k(L) s_k(o) c_k(h)`, where `f0` is a length proposal before start/end selection and capture, and `c` is
a relative capture weight for that whole molecule (its common scale is absorbed into exposure). Let `A` be the
stated mapping and complete-pair acceptance procedure, and `a_k(h)=P(A | h,k)` its probability. Then

```
E_k = sum_h u_k(h) a_k(h)
p(y | k,A) = sum_h u_k(h) p(y,A | h,k) / E_k
P(k,h | y,A) proportional to theta_k u_k(h) p(y,A | h,k) / E_k
```

Here `theta` is an accepted-fragment mixture share. This retains DNA, RNA+ and RNA− as the populations. Transcript,
copy and allele paths are latent explanations within them. Distinct origins are **summed**; redundant SAM records
are not distinct origins. For fixed observed accepted `y`, do not multiply an unconditional sequence emission by
`a(h)` again. Equivalently use `a(h)*p(y | h,A)`—one convention throughout. Density is expected accepted count over
the corresponding effective opportunity. Both numerator and denominator include multimappable origins.

The length kernel's frame matters even without capture. A model that draws a length and then chooses uniformly
among `T(L)` valid starts uses a `1/T(L)` factor. A realized length PMF already selected by available starts must
not be substituted for the preselection kernel without that conversion. Verify boundary effects on short references.

This formulation does not require every genomic position to be uniquely identifiable. It does require a specified
mapping/acceptance model, measured or approximated with uncertainty. A fragment-weighted acceptance probability is
more general than a binary genome mask. On tiny references enumerate it; at genome scale validate a documented
approximation against the actual reporting and mate-rescue procedure. A non-N mask is only a lower-level diagnostic.

**Capture.** A common `c(h)` cancels from a particular posterior ratio; unequal placement capture does not. Effective
lengths still contain capture even in the common-factor case. If starting instead with the observed length law,
factor the density as `p(L | A,capture)*p(start,orientation | L,A,capture)`; a naive additional `c(start)` generally
counts part of capture twice. Capture can depend on length and the full path, not merely its first base. The packet's
unique-fragment numeric-no-op claim therefore needs explicit equal-capture conditions.

**Sum versus mean.** For equal per-start intensity, every distinct feasible origin contributes probability. A
read-specific `1/H` is not an opportunity normalization. A mean can be valid in another parameterization—for
example, conditional on a fixed equivalence class with an explicitly uniform origin prior—but its class weight and
opportunity must then reflect that definition. The shipped count of reported hits does not establish this model.
Correct summation is not a theorem that every finite-sample abundance error decreases. The OFF regression can
reflect invalid support, wrong emissions, opportunity/prior mismatch or an old bias accidentally compensating another.
These are testable possibilities, not established explanations.

**Plug-in error.** With competing RNA weight `r`, external DNA weight `d`, and error `delta`,
`z(d)=d/(d+r)` changes by `r*delta/(d+r)^2 + O(delta^2)`. This is first order. Only a specified expectation with
mean-zero error, or an explicitly orthogonal estimator, can cancel that term. Unique support can make error small;
it does not change its order. No unique support can leave an entire family of regional rates observationally equal.

**Booking and prior.** Sum RNA plus all local/remote DNA responsibilities to one per accepted fragment. Split the DNA
responsibility among its physical origins before adding regional tracks; do not also book its complete weight at
the transcript locus. A calibration hard draw remains one integer deposit, separately from EM's soft accounting.
Changing accepted observations or DNA support requires re-deriving the calibrated prior's expected count and
denominator. An unchanged prior is justified only when the update is a true reparameterization with those quantities
unchanged. Fragment duplication, locus splitting/merging, and shared-copy fixtures must verify conservation.

## 5. Capture-ON diagnosis

**The mean-dilution mechanism is demonstrated; the complete causal explanation is not.** The packet's genic-origin
DNA class moves from **−6,416** additional DNA assignments under the mean to **+41,361** under the sum. The sum also
adds **33,562** false DNA assignments among the two reported RNA multimapper classes. That establishes sensitivity
to aggregation; it does not isolate capture efficiency. Approximately **40.1%** of the NM-comparison categories are
labelled one-mate, so those measurements compare incomplete explanations as well as valid pairs. NM also ignores
soft-clipped bases. The statement that missing capture is the only remaining cause is therefore a conjecture.

The settling experiment is staged: complete-pair support first, then consistent whole-pair emissions and opportunity;
freeze those and compare sum/mean; then sum with/without the derived conditional capture factor. Stratify true
origin, primary/secondary status, fully paired versus rescued, capture overlap and sequence likelihood ratio.
Use oracle capture weights as a diagnostic control, not the release algorithm. Read both origin errors and transcript/
gene tables per condition. Do not bundle all changes into a single favorable ON number.

## 6. Aligner dependence and the full splice design

The minimap2 manual specifies that `splice:sr` uses independent mate chaining and disables secondaries by default;
the panel deliberately enables secondaries. These are search/reporting choices, not whole-fragment posterior
probabilities. STAR's different output density does not excuse the scanner bug. The SAM proper-pair bit is defined
by the aligner; it is useful metadata, not a substitute for validating both mates. Conversely, requiring that bit
would discard valid reconstructions. See the [minimap2 manual](https://lh3.github.io/minimap2/minimap2.html) and
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

Use raw read bases, qualities, CIGAR and reference to evaluate every complete candidate on the same observed bases.
NM alone cannot locate discriminating mismatches; it includes indels and excludes clipping. Missing NM is missing
evidence, not zero error. MD helps reconstruct aligned differences but does not provide an alternative reference's
sequence or missing clipped bases; see [SAM tags](https://samtools.github.io/hts-specs/SAMtags.pdf).
Do not multiply the resulting emission by an aligner score or MAPQ as independent evidence. An independent symmetric
substitution model gives the familiar likelihood ratio `epsilon/[3(1-epsilon)]` per discriminating base, but estimating
epsilon from NM gives alignment divergence, not pure sequencing error. Indels, clips and shared alleles need stated
models and controls. High-NM or short clipped secondaries must not get rewarded for explaining fewer observations.

The revised full specification is `SPLICE_ARTIFACTS_IMPLEMENTATION.md`: whole-pair validation shared by reported,
local and remote candidates; independent handling of mixed genuine/artifact junction paths; explicit uncertainty
for sample-specific structural alleles; and consistent calibration/training support. Both mates need to fit the
same molecule, **not both independently cross a splice junction**. Discrepant reported mates can reflect missing
search/reporting, repeats or alignment errors as well as a chimera. An ordinary candidate must be validated; a
chimeric explanation needs its own coherent path and justified model, never a free orphan fallback.

## 7. Validation gates

Before any implementation, finish the candidate schema, paired acceptance rule, normalized likelihood, observable
opportunity, calibrated-prior accounting and training-selection design. Compare a finite enumerated reference with
each proposed equation. No arbitrary threshold, unexplained density conversion or source-C identifiability claim
may remain hidden in a formula. The owner has made this design gate the first priority.

Before source integration, use an isolated worktree and one-mechanism A/Bs. The minimum fixture set covers every
orphan route; overlapping true short inserts; clipping/indels; duplicate and reordered alignments; absent/inconsistent
HI; complete remote copies; overlapping-mate consistency; mixed genuine/artifact multi-junction reads; unannotated
genuine RNA; and genuine DNA structural alleles. Require exact candidate/count conservation and brute-force agreement.
Every regression gate must first fail on baseline and then fire when its repair is deliberately broken. The packet's
unrun second break gate remains incomplete.

Rebuild invalid scan caches. On the same read identities, compare oracle and aligned input for raw/pool/inferred
length laws, densities, priors, both classification errors and transcript/gene abundance. Separate all simulated
molecules, accepted paired molecules, and mapping losses; pair-conditioned truth is a scoring view, not permission to
discard opportunity correction. Run gDNA-zero controls, equal-length controls, both length-gap directions, and all
three in-scope strata; report the deferred unstranded capture-ON stratum separately. Pin scan threads to one and EM
assignment to fractional. The existing length-law evidence used 16 scan threads; the audit repeats it pinned.

Finally run STAR on SLURM with recorded version/settings, the full aligned panel and independent artifact-containing
evidence; VCaP remains a test set. Retire the blacklist only through a paired comparison of genuine RNA loss and DNA
leakage. Stop release for any incomplete ordinary hypothesis, unaccounted mass, opportunity normalization failure,
unexplained length-law collapse, or an unexplained in-scope error regression. A clean STAR result alone cannot waive
the confirmed scanner defect. No new numerical success margin is invented here: use the standing release criteria
and disclose unresolved ambiguity instead of promising perfect origin identification.
