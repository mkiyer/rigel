# Paired-fragment splice-artifact design

2026-10-15 · Revised after the owner's paired-end and mappable-genome requirements.

**Status: design under review; not ready for production implementation.** This replaces the earlier build sequence
and base-only calibration hold rule. It does not claim that sample-specific structural alleles are identifiable
from every RNA-seq library. `INTERGENIC_GDNA_REVIEW.md` contains the code trace, reproduced failure, statistical
review and build gates. `SPLICE_ARTIFACTS_REVIEW.md` remains the review of the original proposal.

The hold-rule revision is a review recommendation supported by the earlier full-likelihood counterexample;
it is not a claim that a new permanent ruling or implementation has already been adopted.

## 1. The unit is one physical fragment

A read pair is one observation. An ordinary hypothesis must place **both mates on one physically consistent
molecular path**, with compatible orientation and sequence. A reported alignment for one mate is a candidate
location for that mate; it is not an abundance or calibration hypothesis for the fragment.

There are three populations: DNA, RNA+ and RNA−. Local stitching, remote copies and sample structural alleles are
ways a population can produce an observation, not additional populations.

Both mates must support the same molecular explanation. **They need not both contain a CIGAR splice.** A valid RNA
fragment can have one junction-crossing mate and one entirely exonic mate. A reference-DNA explanation must place
both mates on a contiguous genomic molecule; neither may retain a biological splice. RNA can also be contiguous.
If a mate-gap RNA path has an implicit splice, its molecule length is computed on that RNA path; the DNA path has
its own genomic length.

An unmatched secondary is incomplete evidence. It can result from the aligner's search or reporting limits, a repeat,
or alignment error. It does not prove a chimeric molecule. If a valid primary pair exists, keep it and admit a remote
alternative only after recovering and validating both mates there. If no ordinary complete hypothesis survives,
record a rejected/unresolved fragment with a reason; do not deposit a read-length fragment, invent its mate, or score
a single-end survival likelihood. This is bookkeeping about missing inference, not a fourth biological population.

A real chimera is possible. A structural-allele or chimeric interpretation requires a complete path connecting the
observations and a justified probability model; unrelated placements cannot acquire an ordinary genomic insert
length by taking their coordinate extent. The existing chimera exclusion remains until such a model is justified.
Any proposed chimera prevalence parameter must be stated and supported, not supplied as a new hidden constant.

## 2. One candidate representation, before either stage

Each candidate records:

- physical source/path and component compatibility;
- explicit R1 and R2 query-to-path mappings, including how every retained sequenced base is explained;
- the alignment chains used for each mate, with secondary and supplementary provenance;
- mate orientation and outer molecule endpoints on that path;
- path-dependent molecular length, observed junctions, implied mate-gap junctions and population eligibility;
- sequence emission or its sufficient statistics, capture/acceptance context and origin coordinates;
- validation status and a canonical identity independent of SAM record order, HI numbering or duplicate records.

Mate completeness remains available downstream. Both calibration and EM assert that every ordinary admitted
candidate is complete and valid. The representation must make an orphan impossible to mistake for a short insert.
One legitimate overlapping 100-bp molecule remains a 100-bp molecule; read length is not a rejection threshold.

Treat HI and reciprocal mate coordinates as association evidence, with checks appropriate to the aligner's output.
They do not replace path validation. Do not require the proper-pair bit as a universal gate: it reflects the aligner's
criteria. Recover the original mate sequences from the query group when a secondary omits SEQ; missing information
must be explicit. Missing NM must not silently become perfect sequence agreement.

Construct coherent per-mate alignment chains before pairing. A supplementary record is a segment of a chain, not
another independent molecule and not automatically part of every primary hypothesis. Overlapping query segments
cannot be counted twice. Overlapping mates must agree with one common underlying path/allele, allowing the stated
sequencing-error process.

Enumerate feasible pairs without consuming a mate location greedily: the same mate alignment can participate in
multiple alternative complete hypotheses. A shared annotated transcript helps establish RNA compatibility; it must
not suppress a legitimate DNA pair elsewhere. Deduplicate alternative representations of the same origin, while
retaining distinct genomic starts/copies. Bound expensive searches only with a documented completeness guarantee
or a measured truncation limitation. Search limits do not become molecule likelihoods.

## 3. Generate the competing paths

| Candidate source | Construction and validity requirement | What it can establish |
|---|---|---|
| Reported alignments | Assemble and validate both mates, including independently reported secondary locations | Whole-fragment alternatives already present in the BAM |
| Local contiguity | Reconstruct each mate against contiguous reference paths anchored at its aligned blocks; pair and evaluate every base | A short-anchor stitch, including mismatch-tolerant local alternatives |
| Remote reference copies | Use reference copy mappings to propose locations, then place and score both actual reads there | A retrocopy/paralog explanation when the full fragment fits |
| Sample structural alleles | An explicit allele path, both mate mappings, reference/alternate opportunity and an identified rate model | A DNA molecule crossing a sample-specific deletion/rearrangement; currently an open design problem |

For a local candidate, moving a short block may require moving the other mate or changing the implied insert. Never
leave a contradictory mate attached simply because the short overhang matches. For multi-junction reads, a DNA
candidate must explain the whole fragment without RNA splicing. RNA alternatives can retain a genuine junction while
repairing another; independent per-junction penalties do not establish that combined path's existence.

A transcript-to-genome copy map is a proposal index. Covering a transcript interval does not by itself certify the
read pair: mate orientation, sequence, indels, bounds and molecular length still require evaluation. Copies of novel
unannotated junction paths require another proposal/search route. An annotated-only index must expose that coverage
gap. Retain any valid complete remote alternative reported by the aligner, including one without an RNA candidate.

The RNA side must remain honest too. A new or repaired junction path cannot be attributed to an annotated transcript
that does not contain it. The representation/reporting of unannotated RNA paths must be specified before the model
uses DNA as an alternative; omitted annotation must not turn genuine RNA into apparent DNA by construction.

For sample alleles, `expected junction DNA = alpha * density * opportunity` is a possible model term, not a proof
of origin. It requires allele-specific geometry, copy-number/capture assumptions, reference-allele competition and
an explicit prior or identifiable fit for alpha. On unstranded libraries the observed junction mean can be
`RNA abundance + alpha * DNA density * opportunity`; these contributions can be indistinguishable. Expected
capacity is not a hard count bound. Absence of a reference hit is not proof of an allele, and a low-expression true
RNA junction need not exceed expected DNA capacity.

Known splice-end class, motif and distance can be informative only through a stated, selection-conditioned model.
Uniform random breakpoint positions and annotation reuse frequencies do not automatically describe the aligner's
reported junctions. Do not manufacture a prior-free structural-allele discriminator from them. This source remains
release-relevant until its model and withheld-annotation controls are resolved; failure to resolve it limits the
claimed scope rather than justifying confident certification.

Sample-specific retrocopy insertions also belong here: they can mimic annotated junctions at several positions
without any contiguous reference copy. An unannotated-only allele mechanism misses them. A shared insertion path
and its rate need joint evidence; multiplying independent per-junction allele explanations would count it repeatedly.

## 4. Sequence, length, capture and opportunity form one likelihood

Use the normalized generative model in `INTERGENIC_GDNA_REVIEW.md` §4. For each complete physical origin, combine
one sequence emission, one path-dependent length term, one orientation term and one correctly conditioned capture
term. Sum over distinct compatible origins. Do not divide by the number of reported alignments. The same rule applies
to multiple physical placements within a transcript; duplicate alignment representations add no probability.

Evaluate the same observed bases for every candidate. A clipped alternative must explain its clipped sequence through
a specified alignment/error/adapter model; dropping inconvenient bases gives it an unfair advantage. If those bases
are unavailable, record incomplete evidence or integrate over them under an explicit model, not a zero penalty.
Ordinary adapter trimming must be distinguished from genomic evidence, and candidate molecule endpoints must follow
the chosen model. A smaller aligned extent is not automatically a shorter original molecule.

The symmetric substitution likelihood is a useful starting specification, not a license to equate NM with error
probability. At discriminating positions it yields the ratio `epsilon/[3(1-epsilon)]` per additional mismatch.
At an exact sequence tie the ratio is one: the splice is **ambiguous**, not disproved. Use bases and reference, with
qualities where supported; calibrate or state divergence/error assumptions. Indels and shared sample variants require
their own consistent treatment. Do not multiply an aligner score/MAPQ by the sequence likelihood as independent
evidence. Retain no sequence-only calibration hold cutoff: other factors can reverse the sequence preference.

For DNA, define effective opportunity from possible **accepted paired fragments**, integrated over length and
orientation, and capture where applicable. A fragment can be observable through unique or multiple placements;
physical genomic starts are distinct, whereas several records for the same origin are not. Completely unobservable
starts contribute zero. Zero coverage must never be used as a mappability mask.

This is a fragment-level quantity, not simply genome size minus Ns. Paired ends can map across an unsequenced N gap;
non-N repetitive or low-complexity sequence can still be unobservable under the reporting procedure. A reference
mask or acceptance approximation must match the actual candidate generation, mate rescue and acceptance rules.
Validate it by exact enumeration on small genomes and a documented sampling/approximation audit at large scale.
Simulation-generated all-N molecules can remain in a mapping-loss ledger, but they cannot inflate count opportunity.

Choose one length/capture frame. With a pre-capture length law, integrate capture and acceptance in opportunity and
include their joint observation likelihood consistently. With a capture-selected length law, use the conditional
position/orientation law given length and acceptance. Multiplying the current observed length law by an extra capture
factor without the corresponding normalization is not valid. Cancellation for a unique pair requires equal capture
probability across its competing physical explanations.

## 5. Calibration, EM and learning

The stages share the validated candidate space and geometry. Calibration retains the owner's hard seeded draw;
the EM retains soft responsibilities and the per-locus solver. They need not have identical fitted weights at
different stages, but neither may invent a mate or certify an origin that the shared representation cannot explain.

Hold an ambiguous fragment until all required likelihood factors are available, or use a mathematically valid bound
under the existing pruning rule. Do not hold only when the alternative wins the base comparison. A bound must state
which abundance/strand/length ranges it covers; an unbounded remaining factor makes an early exclusion unsafe.

Provisional deposits may be an implementation optimization only after proving equivalence to the held calculation.
Subtracting a count later does not erase its earlier effect on junction support or fitted nuisance parameters.
Specify independent support, exact exclusion of the focal evidence, or a coherent joint estimation scheme. An entirely
held junction has unavailable independent support, not evidence that its RNA path has probability zero.

Learn length and strand from admissible full-fragment observations with the correct selection model. If high-confidence
training is selected by ambiguity, length or geometry, its opportunity changes. Soft RNA training weights are distinct
from calibration's integer deposit ledger. Introduce any refit or joint-fit change explicitly; do not assume that
adding weights fixes selection bias or preserves the current one-pass fit/draw contract.

The length-law audit shows why this matters. The internal contrast density is about 0.000798 while final bulk density
is about 0.0521. Inferred pool purities imply a signed 100-bp probability above five before clipping. Pair validity,
acceptance opportunity, purity estimation and the common-contaminant assumption all need auditing. Clipping a
misspecified inversion into a probability distribution does not establish a valid estimator. After removing invalid
hypotheses, inspect each of these assumptions rather than introducing an arbitrary smoothing strength or bin filter.

Every admitted fragment contributes one total observation. A calibration draw deposits once at its chosen origin.
EM responsibilities sum to one across RNA and all DNA origins; regional DNA tracks receive the split origin mass.
Remote DNA is not also booked as a complete fragment at the transcript locus. Recompute the calibrated prior and
opportunity when the accepted population or support changes. Calibration's absolute density cannot be inserted beside
normalized local mixture weights without an explicit conversion and valid count update.

The preferred design preserves one DNA component per connected locus with a calibrated, normalized spatial law.
Its effective opportunity conditions on the same routing rule as the observed fragments, including remote copies.
Shared remote regions cannot donate their entire prior independently to several loci. Whether this conditioning can
be computed tractably is an open derivation gate; fixed external density and a globally coupled solver remain explicit
alternatives, not silent changes to settled machinery.

## 6. Falsification fixtures

Each fixture gets an exhaustive reference answer before an optimized implementation. The gate must fail on the
defective baseline, pass after repair, and fail when the corresponding repair is deliberately broken.

| Fixture | Required behavior |
|---|---|
| Valid pair plus lone remote secondary, every scanner branch | Original pair unchanged unless a complete remote pair is reconstructed |
| Only one mapped mate | No ordinary abundance hypothesis or fragment-length observation |
| Real overlapping short pair, several read lengths | Preserve its true length; no special rejection of the read-length bin |
| One spliced mate and one exonic mate | Admit their common RNA path |
| Local stitch with a contradictory mate | Reject the invalid DNA path; test legitimate mate relocation separately |
| Complete remote copy pair versus one-end-only copy | Admit only complete compatible origins |
| Genuine and erroneous junctions on the same fragment | DNA and partial-RNA candidates match exhaustive molecular paths |
| Clips, indels, missing NM/SEQ, supplementary chains | Explicit evidence treatment; no gain from silently omitting bases |
| Duplicate/reordered records, HI/no-HI, alternative mate reuse | Same physical candidate set and same probabilities when information is equivalent |
| Unmappable decoy and observable duplicate copy | Correct exposure invariance or increase, with one observed count per molecule |
| Different capture overlap and length laws | Agreement with enumerated probabilities; cancellation only where derived |
| Entirely held and self-supported junctions | No spurious zero support or circular certainty |
| Withheld novel RNA versus indistinguishable DNA allele | No claim to identify information absent from the observations |

## 7. Decisions that must be settled before coding

The owner has settled the complete-pair requirement, mappable denominator, release-blocking length correction,
design-first order and later cluster STAR arm. No local STAR installation or source implementation is part of this
review. The remaining technical design work is concrete:

1. Specify candidate normalization and locus-routing opportunity, including shared remote origins and prior mass.
2. Specify a validated paired mappability approximation and its uncertainty for production-sized references.
3. Settle the read-error/indel/clipping model and training-selection correction; expose any additional assumptions.
4. Establish what fixes the contrast's internal density and shape assumptions after complete-pair acceptance.
5. Specify representation of genuine unannotated RNA and the identifiable scope of sample structural alleles.
6. Verify the completeness/cost contract for the recorded choices of stored reference sequence and index-time
   `mappy`; preserve the explicitly accepted unannotated-copy search gap as a reported limitation, not a solved class.

These are derivations and validation questions, not invitations to add default constants. The next phase is to
close them on paper and with finite reference calculations, then prototype one approved mechanism at a time outside
the main tree. The execution gates are in `RELEASE_0_8_0_FIX_PLAN.md`.
