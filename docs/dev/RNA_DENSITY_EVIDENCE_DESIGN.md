# Honest DNA-density evidence: design for the next checkpoint

*2026-10-07. Provisional design in the development sandbox. The owner approved the evidence and
count-versus-opportunity checkpoint, not a new production capture prior or reader landing.
Execution order is [RNA_SHORT_FIX_PLAN.md](RNA_SHORT_FIX_PLAN.md).*

## The change in ordinary language

An object is a genomic piece or boundary on which Rigel counts fragments. Today calibration
reduces the object's evidence to an estimated DNA count. The prototype capture reader treats
that estimate like a count of positively identified DNA fragments.

The replacement keeps a curve: **how well would these observations be explained if the DNA
density were this value?** A broad curve means several explanations remain possible. A narrow
curve means the data discriminate between them. No variance cutoff decides whether the object
is allowed to speak.

This is an information-preserving interface. It is not a new detector, a promise of perfect
identification, or permission to fit another population model. The count posterior and the
capture posterior may ask different questions of the same evidence, with their priors explicit.

```text
Observed counts + component opportunities + licensed local observations
                              |
                    DNA-density evidence curve
                       /                  \
            DNA count prior          local capture prior
                 |                         |
         count inference             capture correction
```

The count prior can remain the existing library-wide landscape. It must not be smuggled into
the evidence curve and then counted again. A reader using that count posterior directly would
still inherit the population prior; calling it local would not make it so.

## What “distant RNA information” means

There are several distinct paths in the current implementation.

| Path | What it does now | Intended treatment |
|---|---|---|
| Largest total-read density anywhere in the library | Sets the capture slab's upper endpoint; a highly expressed unrelated RNA object changes everyone else's prior | Remove from the intended reader. It has no local biological justification. |
| Landscape learned from inferred DNA on many regions | Supplies the count solver's population prior; RNA/DNA deconvolution errors elsewhere can affect this prior | Retain its declared count-inference role while investigating the profile-based partner. Keep it out of capture's local evidence. |
| Belief-frozen strand widths and evidence precision | The incoming count belief affects a local row's width or admission, indirectly bringing the landscape into messages | Make the primitive depend on observations and the declared observation model. Test rebuilt messages, not just an exported posterior with a term subtracted. |
| Library RNA/DNA reference densities in message lanes | Label internal log-rate coordinates; the RNA coordinate is computed from exon reads | Prove changes of coordinate cancel. If they do not, repair the coordinates/numerics; do not interpret them as local abundance evidence. |
| Library-level lane availability, including `split_live` | Some message capabilities depend on whether a counted single-strand exon exists somewhere | Audit explicitly. A new capture profile must not gain local RNA evidence just because an unrelated exon is expressed. Changing an existing count-channel rule requires a separate ruling and A/B. |
| Strand fidelity and fragment-length laws | Estimate shared properties of the library from reads, including RNA elsewhere | Keep as explicit technical inputs for this checkpoint; hold them fixed in locality tests. Their adequacy is an assumption, particularly for the capture-selected RNA length law. |
| Intergenic DNA background | Provides a shared background location and spread | Retain as the admitted global biological baseline. Unannotated RNA can contaminate it; it is not guaranteed pure by a name or annotation. |
| Connected local observations | Junctions and adjacent objects constrain possible RNA/DNA explanations through the existing geometry | Necessary, especially without strand information. Genomic proximity alone is not a license to share RNA abundance. |

A locality guarantee must therefore say what is held fixed. The target is:

> Given the same local observations, opportunities, technical calibration, background and
> declared capture prior, changing unrelated RNA abundance cannot change this object's capture
> evidence or correction.

Count estimates may still change through the explicitly retained population prior. A genuine
change to a shared technical estimate or the background is also a different experiment. There
is no claim that every output of the tool is independent of every other locus.

“Local” means connected by the existing fragment/strand/splice relationships, not within a new
chosen number of base pairs. It does not mean borrowing a nearby gene's expression level.

## Counts and opportunities answer different questions

- A **realized count** is how many observed fragments originated from DNA.
- A **rate** predicts the distribution of such counts over repeated sampling.
- An **opportunity** says how many placements a component's fragment-length distribution permits.

For DNA density `rho` and DNA opportunity `Eg`, the expected DNA count is `rho Eg`.
Its realized count is random. An unknown origin count is not an observed Poisson count.

RNA uses its own opportunity `Er` wherever RNA rates are compared or transported between objects.
No repair substitutes DNA's length law for RNA's except as an explicitly labelled oracle diagnostic.
A rate posterior is allowed above `observed total / Eg`: a Poisson expectation is not bounded by
its realization. The actual latent DNA count remains between zero and the observed total.

This separation is essential to the gap problem. A tiny DNA opportunity does not turn a noisy
fraction of an inferred fragment into precise evidence of a huge DNA density.

## The observation model and the primitive

Keep the existing three populations: DNA, RNA on the positive strand, RNA on the negative strand.
For an object's two observed columns `u+, u−`, define:

```text
a = rho Eg                         expected DNA fragments
r+, r−                             expected RNA fragments admitted at this object
mean+ = a/2 + kappa r+ + (1−kappa) r−
mean− = a/2 + (1−kappa) r+ + kappa r−
u+ ~ Poisson(mean+)
u− ~ Poisson(mean−)
```

These are the substrate's integer incidence counts, not its conserved fractional mass.
Separate objects can still count the same fragment; this does not establish independence
between their observations.

An unadmitted RNA strand has zero RNA amount. A channel with zero RNA opportunity also has
zero RNA expectation; do not integrate a nonexistent component. RNA rate comparisons at positive
opportunity use `r_s / Er`, never `r_s / Eg`.
This is the existing Poisson/strand observation assumption made explicit. It is not a claim that
mapping errors, duplicate molecules or protocol heterogeneity never violate it on real data.

**Recommended nuisance treatment:** retain the existing RNA reference while summing over the
unknown RNA amount. For one admitted RNA strand this is proportional to `r^(−1/2) dr`.
For both strands, retain the existing total-RNA reference and the existing pure-positive,
pure-negative and mixed-tilt hypothesis measure. Do not silently introduce independent RNA
reference priors that change the total-RNA measure.

The resulting function is a **marginal evidence curve**, not a normalized probability distribution
over DNA density. “Free of the DNA prior” does not mean free of all modelling assumptions: the RNA
reference remains explicit. Nuisance integration and maximization are different models; the old
own-read prototype's profiled maximum is an ablation, not silently the same formula.

In an unstranded object, the own reads constrain the total amount but cannot identify its split
into RNA and DNA. Any preference along that ambiguity comes from the stated RNA reference,
not from a newly discovered strand witness. This is a reason to retain justified local
observations and test reference sensitivity, not to promise identification from two columns alone.

For a local observation factor `H`, the intended calculation is schematically

```text
L(rho) = integral
           Poisson(u+; mean+) Poisson(u−; mean−)
           H(rho, RNA amounts, tilt)
           [RNA reference measure].
```

With no additional factor, the one-strand calculation below is exact under this observation model.
For neighbours, `H` is a specification to derive and verify, not a license to multiply today's
already-fused rows and call the result an exact likelihood.

### An independent small-count reference

Let `M=u+ + u−`, `k` be the unknown number of realized DNA fragments, and suppose the admitted RNA
strand has positive-column probability `q` (`kappa` or `1−kappa`).
Conditional on `k`, the positive column is a sum of two binomials:

```text
DNA-positive ~ Binomial(k, 1/2)
RNA-positive ~ Binomial(M−k, q)
S(k) = P(DNA-positive + RNA-positive = u+).
```

Then, up to a constant independent of `rho`,

```text
L(rho) = sum over k=0..M of
           Poisson(k; rho Eg)
           Gamma(M−k+1/2) / Gamma(M−k+1)
           S(k).
```

The Gamma ratio comes from integrating the RNA amount with its existing reference. No DNA prior
is present in this expression. The two-dimensional latent-allocation sum and direct integration
of the original two Poisson columns agree in 140 checks to relative error below `9e−15`.

**A correction to the earlier review's exactness claim:** evaluate the conditional convolution
`S(k)`, not the usual strand-mixture likelihood at `f=k/M`. For `M=2, k=1, u+=1, kappa=.99`,
the conditional split probability is 0.5; the plug-in binomial gives 0.37995.
The distinction matters at small counts. At large counts any faster approximation must be
checked against the correct reference, not against another implementation of the plug-in.

Limits and invariants:

- No admitted RNA: the density-dependent factor is `Poisson(M; rho Eg)`.
- No observed reads: `L(rho)` is proportional to `exp(−rho Eg)`. Zero observations are evidence;
  they are not automatically a flat likelihood.
- `Eg=0`: this object's own observations do not depend on `rho`; for RNA-explainable data the
  own curve is constant. Positive counts impossible under every admitted component are an
  opportunity/model inconsistency, not a reason to invent a density.
- The Poisson tail above `M/Eg` remains. There is no count-derived prior ceiling.
- Multiplying one object's likelihood by a positive constant changes no posterior or population
  fit. Normalizing a curve for numerical stability is allowed; averaging such curves as if each
  were an observed density distribution is not the population likelihood below.

The independent reference and illustration are in
`.cache/rigel_runs/2026-10-07_density_evidence_design/`
(`check_derivation.py`, `derivation.json`, `density_evidence.png`).
This calculation proves the primitive, not the neighbour model or release accuracy.

## Neighbour evidence: the part that must be certified

Own reads alone remain insufficient and are not proposed as the release reader.
Unstranded and both-strand objects need the existing local relationships. The implementation
must distinguish three ingredients without building a new inference framework:

1. Raw observations and their likelihood factors.
2. Background-derived constraints, especially the intron factory.
3. The target object's DNA density prior.

Use separate arguments/arrays and a diagnostic provenance ledger. Do not pass their sum as
“evidence.” Do not reconstruct it by subtracting a landscape term from a posterior.

Audit these existing sites first:

- `native/transfer_rows.h`: the strand likelihood and its frozen variance.
- `native/transfer_kernel.h`: own claims, licensed faces, RNA level rows and their rate coordinates.
- `native/solve_kernel.cpp`: factory construction, evidence precision, received rows and final assembly.
- `calibration/messages/transfer.py`: shared coordinates and lane-availability inputs.

A local RNA factor belongs **inside** the nuisance integration when it depends on the RNA amount.
A row about an intensity must not be reinterpreted as a row about a realized count through
`k/M` without a derivation. Spliced fragments certify RNA; they do not measure all unspliced RNA,
and retained RNA is not a fourth population.

The required tests rebuild the local evidence while changing only a distant expression tally,
the incoming posterior belief, a coordinate origin, or a density-prior curve. Holding already
contaminated rows frozen proves only the arithmetic downstream of those rows.

Shared fragments and repeated messages are another concern: two boundary counts need not be
independent. Keep the existing distinction between bounds and evidence, and test that supplying
the same physical observation twice does not sharpen a curve. Do not solve a failure with a new
arbitrary variance multiplier.

**The factory is an explicit unresolved assembly issue.** Its background-predictive DNA-count
row is not a new independent DNA observation. It must not be trained into the landscape as if
it were data and then multiplied beside that landscape at the same object. The first audit
will show which neighbouring claims depend on it. Do not automatically delete it, substitute
a count posterior as an RNA message, or choose separate exon/intron landscapes: those changes
either lose established information or revisit refused mechanisms. Bring a concrete
single-prior assembly for the affected face back for the owner's ruling before a count-changing
implementation that requires it. **2026-10-07 update:** the owner authorized the assembly
specified in [the checkpoint](RNA_DENSITY_EVIDENCE_CHECKPOINT.md), including factory removal.
This resolves permission to investigate, not the statistical adequacy of that replacement.

## The proposed partner for applying the DNA prior once

The landscape currently trains on median-derived point counts with Poisson kernels and
additional reliability/admission rules. The candidate partner is to fit the **same one
population prior** from evidence curves instead.

On a density grid with atom masses `pi_j`:

```text
fit pi by maximizing sum_i log(sum_j pi_j L_i(rho_j))
subject to pi_j >= 0 and sum_j pi_j = 1

posterior_i(j) proportional to pi_j L_i(rho_j).
```

A flat curve contributes no information about the shape of `pi`. A broad curve expresses its
own uncertainty. There is no need to turn its fitted median into a precise pseudo-observation
or to invent a second prior for the capture reader.

This is a candidate statistical fit, not an established cure. Initial contrasts keep the
training objects, grid and existing weights fixed to isolate the kernel change. Unit object
weighting and removal of posterior-based admission/weighting form a subsequent explicit
contrast and require amending the relevant training ruling before landing. Preserve the
structural training exclusions initially; do not silently add boundaries or both-strand regions.
Where neighbouring profiles reuse data, this is a composite likelihood, not a product of
independent observations; its influence and uncertainty must be tested accordingly.

The numerical fitting reference is now independently tested. The larger-panel continuation
also makes an extra control necessary: changing the averaging estimator to this likelihood
fit must be measured with its input curves held fixed, before replacing those inputs.
[RNA_LANDSCAPE_PROFILE_CHECKPOINT.md](RNA_LANDSCAPE_PROFILE_CHECKPOINT.md) specifies that
order and the narrowly scoped posterior-variance admission question. The reference has not
been connected to calibration, and it does not certify the missing local factors.

After the partner passes, the count consumer applies the chosen DNA prior once. The old DNA
Jeffreys half must not remain as an additional tilt when a fitted DNA prior replaces it.
The RNA reference stays. Keep the existing count-summary policy in the first comparison;
a change from medians to means is a separate mechanism, not part of “prior once.”
A density marginal and a realized-count marginal are different outputs; do not report
`rho Eg` as if it were the observed allocation.

The prior-once change does not land alone. If the profile fit cannot preserve supported capture
while fixing false density evidence, reject the partner rather than add class-specific offsets.

## Capture: the information boundary, not another new fit

The intended capture reader consumes the local evidence plus the intergenic background and an
explicit capture prior. It does not consume the count posterior, a library RNA-density maximum,
or a second learned population of capture levels.

There is a remaining prior decision. The current normalized log-uniform slab requires an upper
endpoint; making that endpoint numerical rather than statistical does not fix its normalization.
The current `1/n` mixture weight also depends on annotation size. Neither choice is ratified
by the owner's approval of this checkpoint.

Do not adopt the previous review's `1/C²` candidate automatically. A proper unbounded tail could
remove the endpoint, but it is still a new prior assumption and may suppress thinly witnessed
capture. The checkpoint can be completed without choosing that tail: compare likelihoods,
counts and opportunities first; keep the old reader only as a frozen diagnostic consumer.
A new reader still has its owner checkpoint.

Nor does a correct DNA density identify RNA capture perfectly under a length gap. DNA and RNA
average different placements. Expected-yield comparisons must separately test the inherited
DNA-to-RNA transfer. No probe BED enters inference; simulator probe information is an oracle only.

## Implementation size and acceptance

Start with one reference function evaluating log evidence at requested densities, plus the
existing local factor construction and one prior-application site. Work in log space, include
the zero-rate endpoint analytically, and normalize only after a prior is applied.

The latent-count sum is a small-case oracle, not an instruction to enumerate every possible
origin count for every whole-genome object. Use the analytic pure-DNA path; derive and validate
a bounded-memory numerical evaluation for mixed objects. Numerical refinement controls an
error tolerance, not a biological count threshold. Stream through the existing blocks instead
of storing an object-by-density cube for the human genome. Profile before a panel campaign.

Do not add an inference framework, a detector, new per-class tuning constants, a distance window,
a second population fit, or retained production switches for rejected experiments.

The checkpoint succeeds when:

- raw evidence and prior-derived information have an explicit, tested boundary;
- the core and admitted neighbour constructions match independent references;
- the count/opportunity contrasts locate the important residual errors by object class;
- the population-partner experiment has a falsifiable specification and passes before prior-once
  is paired with it;
- the release reader's still-open prior and transfer decisions are stated rather than hidden.

There is no guarantee that local observations can distinguish capture from copy number,
mappability, unannotated transcription or RNA where DNA is scarcely observed. The design must
remain uncertain in those cases. Simulated accuracy is evidence about a candidate, not a reason
to add an otherwise unjustified assumption.
