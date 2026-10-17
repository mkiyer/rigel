# Review of the staged mate-pairing implementation

2026-10-15. Review of `<worktree>`, based on commit
`41899ae4b8209893f7b5ae0cf893304e6796dabf`, and the supplied implementation map and evidence bundle.
This is a review and an execution recommendation, not a source repair.

The shared assembler, complete-pair enforcement, exhaustive pairing and separate stage comparisons are useful
progress. The observed-path clip walk in implementation-map section 8 is necessary and agrees with the endpoint
model already reviewed. I would implement that walk after correcting location normalization and making the
path-to-target contract explicit. I would not build section 8 verbatim or treat stages 2–3 as fully validated yet.

I read the native changes, Python integration, reference and gates, transition scripts and outputs, mutation
summary, full-suite log and oracle identity logs. The current changes to the four native files match the saved
`stage3.patch`. I ran the supplied stage-3b wheel through the pairing and scanner/accumulator gates: **36 passed,
five failed**, with exactly the reported future clip/chain failures. I inspected, rather than reran, the full-suite
result of 2,982 passed/five failed and the mutation sweep. I also ran nine additional synthetic probes and replayed
both supplied disappeared-from-buffer subsets through the native assembler. No production source, worktree source,
compiled extension, benchmark input or installed environment was changed.

The new probes use native module SHA-256
`7560ae249204bcca428aa2f1aae3743f3a9ea10d0c7c4301dc66625d5a26ea25` from `site_s3b`.
Their [script](paired_fragment_review_evidence/pairing_stage3_audit.py) and
[results](paired_fragment_review_evidence/pairing_stage3_audit_results.json) are reproducible.
These synthetic alignments test representation and geometry, not sequence-alignment accuracy.

**1. Fix record-to-location normalization before implementing clips. High priority; reproduced.**

[`collect_mate_locations`](../../../rigel-pairing/src/rigel/native/bam_scanner.cpp:888) currently puts all records
for one mate in an HI group into one location. Deduplication runs afterward. Consequently:

- Repeating the identical R1 record within one HI group changes the fragment's NM from **2 to 4**, with no
  duplicate conflict reported. The repeat must be an alias, not additional sequence evidence.
- Two distinct nonsupplementary R1 records in one HI group are concatenated into one mate location. My probe
  uses starts 1100 and 1110; the assembler reports one placement containing both R1 alignments. An HI association
  cannot certify that these records are segments of one alignment.
- The location key omits clipping and query coordinates. At the same position, `5S95M` and `95M5S` collapse to one
  location. Under the approved endpoint model they imply different molecule starts and lengths, 300 versus 295 bp
  in this fixture. Stage 4 cannot recover a hypothesis already discarded during deduplication.

Parse each nonsupplementary alignment into its own location first; deduplicate equivalent representations; attach
HI/primary associations and raw record aliases as provenance. Keep supplementary-chain metadata separate. If one
HI/mate association contains conflicting standalone alignments, retain them as separate reported locations or
diagnose the malformed association; never merge their blocks into a synthetic mate. Missing HI is not HI zero.

The key must preserve every property needed by validation and scoring: reference/orientation, query-to-reference
geometry, raw/retained junctions, terminal soft and hard clips, and endpoint/chain status. Separate matched query
coverage from deletion spans. The same outer blocks alone do not establish equivalence. Homopolymer-shifted indels
can still be deduplicated when their normalized representation and model evidence are equivalent; raw CIGAR-string
equality is not the proposed replacement.

The least-`(NM, motif strand)` representative is also a policy, not merely deterministic bookkeeping. A probe with
equal geometry and opposing XS tags selects the positive-strand report and becomes annotated-splice evidence;
the conflict counter does not undo that evidence. Inspect the 110 real conflicts before freezing this behavior.
For genuinely equivalent alignments, reconcile evidence using a justified normalization/recomputation when
available. Preserve unresolved strand conflict as ambiguity, not an integer-order preference. If differing
sequence alignments remain, preserve their emission alternatives under an explicit rule without counting repeated
reports as extra physical placements. Do not infer that a lower NM is automatically the correct report.

**2. The new splice-tolerance rule is not a certificate of one molecular path. High priority; reproduced.**

[`splices_contradict`](../../../rigel-pairing/src/rigel/native/bam_scanner.cpp:1009) permits an overlap of at most
three bases between a mate's aligned block and the other's skipped interval. This accepts more than small
uncertainty at a donor boundary. My forward mate is `50M3N50M` at 1000; the reverse mate is `100M` at 1050. It reads
through the entire three-base skip. Stage 3 accepts the pair and deposits **147 bp**, removing three bases that
the other mate explicitly aligns. A strict mate-consistency control rejects it.

Reusing an existing numerical tolerance does not establish the new interpretation. A literal common path cannot
both retain and remove those bases. The independent reference duplicates this permissive rule, so agreement with
that reference cannot prove physical consistency.

Recommendation: make the structural mate-versus-mate witness strict. Keep the existing annotation/emission
tolerance separate. If short boundary errors are to be rescued in this fix, represent an explicit repaired
alignment hypothesis and its evidence; do not silently erase aligned bases and call the result exact. That rescue
would need its own derivation and A/B. This recommendation adds no new cutoff.

There is a related completeness hole: `100D100S` has no query base aligned to a reference base, but is accepted as
a mate and deposits a 300-bp fragment with its partner. `parse_cigar` puts D in its reference-consuming blocks, and
the new `any_block` guard therefore passes it. Require positive M/=/X query-to-reference support for each location;
do not equate reference advance with aligned sequence. This is an adversarial parser case, not a measured real-data
prevalence. SAM distinguishes these CIGAR operations in its
[specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

**3. Carry a validated path per interpretation; do not repair unsupported transcript assignments by copying a
short length. High priority; partly acknowledged stage-5 work, with a stage-4 design dependency.**

The existing [`shortest_path_length`](../../../rigel-pairing/src/rigel/native/bam_scanner.cpp:1047) answers whether
some witness fits. It does not establish which transcript/source interpretations fit or even which target the
witness supports. The synthetic intergene splice is 300 bp on its observed path, but reaches the buffer with four
11,200-bp candidates: both annotated transcripts and their two synthetic RNA entities. This reproduces the supplied
3,544-fragment concern in a tiny example.

Changing all four lengths to 300 is not a correction: neither neighboring annotated transcript contains that
observed intergene splice. Physical plausibility, an annotated transcript explanation, and RNA-versus-gDNA evidence
are different questions. Disabling the old gene-based chimera heuristic after geometric validation is reasonable;
it does not certify any RNA target returned by the resolver's fallback unions.

A canonical path record should carry the retained segments/cuts, inferred endpoints, length, supporting targets,
and source eligibility. Observed alignment evidence remains separate. Acceptance is existential over surviving
paths; each scored interpretation must point to its own surviving path. A minimum is a diagnostic, not a substitute
for this set. Invalid/nonpositive paths must be filtered before taking a minimum.

The owner subsequently clarified the scope: annotated transcript estimation remains the task, a separate
gene-centric estimator is deferred, and unsupported spliced placements are reported separately. The operational
contract is now in [Spliced placements without a supported transcript](MATE_PAIRING_FIX.md#spliced-placements-without-a-supported-transcript).
This resolves the reporting/abundance scope raised during this review; it is not yet a production behavior.

Stage 5 also needs explicit cases for a valid gDNA path with zero RNA candidates, and for valid alternatives
elsewhere in a multimapper. The scorer currently conditions some gDNA handling on `n_cand > 0`, and treats
nonpositive lengths as neutral. Neither is a safe way to represent a path deliberately removed by the bound.
Use explicit support and ensure invalid gDNA cannot reappear through a neutral score or zero-candidate fallback.
This does not settle the separately parked mean-versus-sum multimapper likelihood.

**4. Approve the clip walk in principle, but revise its path and chain contracts before coding it.**

Walking the observed path across a skip, as in section 8's 1600–2000 example, is correct. The endpoint at 1595 is
on that path; 1995 is in its removed intron. Per-hypothesis endpoints must reach occupancy, deferred storage,
second-pass length evaluation and start-dependent coverage scoring. Those changes implement the already reviewed
model rather than introduce another biological assumption.

The unresolved details are consequential:

- Define the complete candidate path first and walk that same path at both ends. An observed-plus-gap witness and
  a transcript's full exon chain are not automatically the same path. Preserve the observed skips and only add
  compatible, supported gap/extension continuations; do not union incompatible splice alternatives.
- `tx_frag_length(u, v, t)` is not a general path-length function. It subtracts t's internal introns, ignores novel
  observed skips, and can subtract an intron covered by a read. Section 8 must not describe that value as the
  validated path length. Compute length from the path's retained segments once and consume it consistently.
- Section 4's overhang exception needs an explicit boundary. Keeping an existing error-emission penalty is a
  separate choice from asserting that its candidate has a physically supported path. Extending along t's exons
  when an anchor is inside them, but contiguously otherwise, does not by itself define one coherent target path.
  Keep the legacy approximation visibly separate in diagnostic arms; do not let it certify a target or bypass
  containment, exhaustion or the bound in the final implementation.
- Chain detection cannot rely only on SA presence. SA is an optional tag; supplementary records may establish an
  unresolved chain without it. Inspect chain metadata before discarding records, associate it conservatively,
  and withhold locations needing reconstruction. Preserve another independently represented valid alternative.
  Add a no-SA supplementary-chain test. See the [SAM tag definitions](https://samtools.github.io/hts-specs/SAMtags.pdf).
- Store endpoint coordinates with explicit conventions: increasing genomic/path direction versus transcript
  strand direction. Test reverse-strand transcripts, not just swapping R1 and R2. Make deferred-bank comparison,
  compression, append/merge, serialization and rehydration preserve the new fields.

The input contract remains the one already supplied by the owner: remove non-template 5′ sequence and preserve
template 5′ bases. Do not reopen that policy. A STAR command line alone is not evidence about earlier preprocessing;
record the benchmark assumption rather than claim its history was verified.

**5. Strengthen the reference gate to test the promised object. High priority before accepting stages 4–5.**

The gate is independent of the native interval arithmetic, which is valuable. Its advertised coverage is broader
than its assertions. The native seam exports the shortest length, not every path. The test compares that minimum
and lengths for returned targets; it does not require the exact expected target set. Synthetic RNA entities are
skipped, as are overhang candidates. Its `transcript_length` deliberately reproduces the old EM rule rather than
the new path contract. See
[`test_the_placements_are_exactly_the_specified_pairs`](../../../rigel-pairing/tests/test_paired_placements.py:632).

Expose and compare complete canonical path sets, endpoints, lengths, target/source support and record membership.
Use a tiny base-by-base reference that independently checks that each claimed aligned base lies on the chosen
path, with any explicitly modeled alignment error identified separately. Test missing as well as extra targets.

Add mutation gates that remove a nonshortest valid path, retain an out-of-bound candidate beside an in-bound one,
replace an observed-path length with `tx_frag_length`, exchange two equal-length paths' starts, drop a new endpoint
field during deferred-bank compression, and revive an invalid gDNA path via a sentinel. Each must fail even when
the placement count and shortest length stay unchanged. Add duplicate-inside-HI, clip-distinct-location, conflicting
tag, deletion-only mate, and supplementary-without-SA cases. Do not weaken the tests merely to preserve identity.

**6. Correct the transition interpretation and pin the numerical comparison. Measured.**

The supplied `transition_dump.py` explicitly says that absence from the buffer can mean an intergenic multimapper,
not only rejection. `unpaired_breakdown.py` then treats the absent IDs as unpaired and diagnoses only their primary
pair. I replayed those extracted BAMs through the native assembler:

| supplied disappeared-from-buffer subset | subset size | actually unpaired | still placed, all alternatives with zero transcript candidates |
|---|---:|---:|---:|
| minimap2 | 12,071 | **12,055** | **16** |
| VCaP | 46,875 | **46,865** | **10** |

The full-run ledger's net increase in `n_unpaired` is **12,055** for minimap2 and **47,096** for VCaP. The VCaP
disappeared-buffer subset therefore misses additional unpaired transitions as well as including ten placements.
Use an explicit outcome per eligible fragment, before buffer routing, for the final transition census.
The [replay script](paired_fragment_review_evidence/pairing_stage3_fate_audit.py),
[minimap2 result](paired_fragment_review_evidence/pairing_stage3_fates_mm2.json), and
[VCaP result](paired_fragment_review_evidence/pairing_stage3_fates_vcap.json) preserve this check.

The primary-pair reason breakdown is diagnostic, not the assembler's exact final reason distribution. For example,
the claimed 2,669 VCaP clip rescues include 68 with SA, which the proposed chain policy withholds. A successful
contiguous containment calculation does not establish survival of every other gate or absence of a valid alternative.
Likewise, a simulated true molecule of at most 500 bp does not alone prove every rejected alignment is impossible:
unmodeled paths, clipping and implementation errors must be ruled out using truth coordinates and rejected
candidates. The inspected examples establish errors in those examples. Outward orientation alone does not establish
a circular-RNA back-splice.

The supplied transition runs use **four scan threads** (`run_transitions.sh` and the stats JSONs). Their identity
joins can still count deterministic placement changes, but the post-drain length histogram is not the pinned
one-thread A/B required for a numerical attribution. Recheck the length-law changes with one scan thread and
`OMP_NUM_THREADS=1`, then use fractional EM for abundance comparisons. The large removal of singleton deposits
remains compelling; the smaller subsequent histogram difference needs that controlled comparison.

**Recommended execution order.**

1. Freeze this review's failing microfixtures. Repair normalization/aliases and require actual aligned query
   support. Settle strict structural splice consistency versus an explicit error-repair extension; the
   recommendation here is strict consistency for this fix. Inspect the 110 real duplicate conflicts.
2. Make path witnesses and target/source projections explicit. Implement the owner's clarified
   no-supported-transcript disposition. Extract shared geometry below its consumers, so calibration and scoring
   consume one result rather than re-derive different lengths.
3. Implement endpoint storage and round-trip gates, then the conditional soft-clip walks and conservative chain
   handling in their own comparisons. Clip-free, unchanged-support fixtures should remain numerically identical;
   newly corrected defects get separate arms and explanations.
4. Implement per-interpretation bounds and target filtering. Test a short RNA/long genomic span, two RNA paths on
   opposite sides of the bound, gDNA-only support, a novel-only path, and survival of a different valid placement.
   Carry the validated start and length into both unique and multimapper scoring branches.
5. Complete the per-fragment outcome ledger, annotation parity, cache/signature invalidation and capacity checks.
   Include clear missing-mate, unresolved-endpoint, unresolved-chain, incompatible-geometry and length reasons;
   counters with overlapping meanings must not be presented as a partition. Do not silently redefine
   `mate_unmapped`. Check all new parallel arrays and integer narrowing before promotion.
6. Run the final one-thread, fractional comparisons on the required strata, including oracle and aligned controls.
   Read pool occupancy, the complete length PMF, the signed length contrast, final calibration, both origin errors,
   and transcript/gene tables. A reduced 100-bp spike is not proof that the separate density/contrast defect is fixed.
   Re-measure the already parked problems on this accepted-pair implementation before the full splice-artifact
   likelihood work resumes. Promote only after the full reference, mutation, persistence and A/B gates pass.

The immediate work does not require a new clip threshold, aligner installation or source-emission heuristic. The
owner has resolved the novel-path reporting/abundance scope through the linked contract. The other findings and
implementation gates remain open, including explicit calibration eligibility and mixed supported/unsupported
alternative accounting; separate reporting alone does not solve their likelihoods.
