# Mate pairing: reviewed design and execution plan

2026-10-15 · Reviewed against `41899ae4b8209893f7b5ae0cf893304e6796dabf`.
Sandbox design; no production implementation or merge approval is implied. This replaces the earlier proposal
in this file, preserved in `paired_fragment_review_evidence/mate_pairing_proposal_reviewed.txt`.

## Verdict and decisions

**Repair the paired-placement foundation before resuming the splice-artifact prototype, but do not implement the
original proposal unchanged.** Requiring two mates, sharing assembly between scanning and annotation, filtering
multimappers after assembly, and excluding unmatched mates from ordinary inference are correct directions.
The proposed geometry, hit counting, supplementary handling and validation promises need the changes below.

The owner settled four decisions during the initial review, then reopened the 5′-clip part after the input census.
The other decisions stand:

- Apply the same molecule-compatibility checks to aligner-associated and constructed pairs. Measure their effects
  separately; do not ship the asserted-pair exemption.
- Stop with an actionable input diagnostic when withheld alignments are established. Do not manufacture candidates
  or report ordinary quantification from a known truncated reported set.
- Make `scan.max_frag_length` an explicit configurable acceptance bound on each candidate's **molecular** length.
  Measure this semantic change separately. The current default remains 1,000 bp.
- Withhold candidate locations requiring supplementary-chain reconstruction, report the reason, and retain other
  fully defined pairs. Ordinary 3′ clipping remains supported. Supplementary reconstruction remains separate work.

**Revised 5′-clip recommendation under the owner's stated tool scope:** replace blanket withholding of terminal
soft clips with the explicit templated endpoint model below, measured in its own stage. Users must remove
non-template 5′ sequence and should preserve the 5′ insert bases. The length model assumes those template bases
were not trimmed away without retained offsets. It is exact conditional on its assumptions, not an exact endpoint
measurement from CIGAR. Do not introduce a small-clip cutoff based on the observed 1–5-bp majority.

This work repairs representation and acceptance. It does not change the two-pool length estimator, capture model,
mean versus sum, intergenic EM routing, alignment rescue, mappability opportunity or the full splice-artifact
likelihood. Those remain in `PARKED_WORK.md`. Correct pairing is necessary, but does not by itself establish
accurate abundance or clear the release-blocking length defect.

## Critique of the original proposal

| Finding | Consequence | Revision |
|---|---|---|
| Asserted pairs bypass validation | A primary flag or matching HI tag can still admit a physically incompatible pair. | One validator; aligner association supplies provenance only. |
| Shared transcript bypasses the length limit | A positive resolver length can exceed 1,000 bp; sharing an exon proves little about molecule size. | Verify each path and bound its actual molecular length. |
| Facing predicate is too weak | Forward `[1100,1200)` and reverse `[1050,1150)` pass `forward.start <= reverse.end`, although both protrude beyond their implied 5′-endpoint interval. | Check both aligned paths against both endpoints; retain legitimate overlap. |
| Merging blocks loses contradictory evidence | One mate can skip an intron while the other aligns continuously through it. | Validate separate mate CIGAR paths before merging. |
| Withheld NH entries count as accepted placements | Their mates, geometry and lengths cannot be validated. | Separate reporting status from accepted count; stop on established truncation. |
| Supplements remain attached to primaries | Different query chains can be combined and NM/query sequence counted twice. | Withhold chain-dependent locations under the approved bounded policy. |
| Shared helper is claimed to guarantee writer parity | Producer filtering, IDs, config, fallback tags and first-hit tagging can still disagree. | Share the whole assembly contract and record membership; test parity. |
| The 1,000-bp rule is described as existing behavior | Config calls it a histogram maximum; scoring extrapolates beyond it. | Implement the newly approved support consistently, with its own comparison. |
| Approximate census becomes an exact prediction | The script uses greedy matching and a 2,000-bp start-distance rule without the proposed path/orientation checks. | Recompute exact transitions; withdraw 501k unique conversions and the 318-fragment output bound. |
| Cache rebuilding is optional | Existing cache keys exercise deposits, not BAM assembly. | Mandatory, tested rejection of stale assembly caches. |

New native probes in [`pair_geometry_audit.py`](paired_fragment_review_evidence/pair_geometry_audit.py), with
[`results`](paired_fragment_review_evidence/pair_geometry_results.json), confirm several gaps. The current scanner
buffers a 4,100-bp pair sharing one exon; accepts an outward primary pair and deposits 400 bp; deposits 150 bp from
the protruding overlap above; and deposits **125 bp** from an incompatible spliced/continuous pair whose second
mate alone aligns **150 bases**. Genuine paired 100-bp and 50-bp inserts pass the controls. These are synthetic
representation probes, not alignment-search or accuracy benchmarks.

Key code sites are `native/bam_scanner.cpp` (parsing, `build_fragment`, grouping, scanner and writer),
`native/constants.h::detect_chimera`, `native/resolve_context.h::tx_frag_length`, `native/scoring.cpp` (both length
likelihoods), and `scan_cache.py`. In particular, `detect_chimera` skips blocks with empty transcript sets and can
return early for connected transcript sets; `tx_frag_length` subtracts internal transcript introns without proving
that neither read aligned through them. Neither is a whole-molecule certificate.

The [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf) distinguishes primary, secondary and
supplementary records; proper pairing is an aligner judgment. The
[SAM tag specification](https://samtools.github.io/hts-specs/SAMtags.pdf) defines NH in terms of reported query
alignments, not validated molecules. STAR's output loop supports alignments with only one mapped mate and can
write fewer alignments than its declared count; see the
[STAR output implementation](https://raw.githubusercontent.com/alexdobin/STAR/2.7.11b/source/ReadAlign_outputAlignments.cpp)
and [output-cap parameters](https://raw.githubusercontent.com/alexdobin/STAR/2.7.11b/source/parametersDefault).

## Molecular contract and representation

A placement is a pair of reported mate locations with at least one compatible ordinary molecular path. It places
**both mates**. Transcript/gap alternatives within that pair are not additional genomic placements. Only one mate
need cross a junction; the other must be compatible with that same molecule.

An unmatched record is neither a fragment-length observation nor evidence of a biological chimera. A rejected
cross-combination cannot make an otherwise valid fragment chimeric. This fix has no single-mate likelihood and
does not infer a missing mate from TLEN, PNEXT or the mate-unmapped flag.

Return a shared assembly result rather than anonymous pairs of record vectors:

```
FragmentAssembly
    fragment identity and stable scan ID
    raw-record identities and filter outcomes
    normalized R1 and R2 locations
    reporting status and aligner-associated pair provenance
    accepted placements in canonical order
    raw-record -> accepted-placement membership
    terminal fragment outcome and rejection reasons

Placement
    one R1 location and one R2 location
    canonical identity and raw-record aliases
    verified molecular paths and their lengths
    eligible transcript/gap hypotheses and genomic support
```

Each mate location retains reference, orientation, full CIGAR/query intervals, splice skips, separate soft/hard clip
counts, endpoint basis (aligned, conditionally inferred, or unresolved) and scoring evidence. The merged
`AssembledFragment` becomes a downstream projection of a validated pair, not its
certificate. The production boundary must reject an incomplete placement in release builds; a debug assertion
alone is insufficient. No fallback can manufacture a placement from an unused record.

### Normalize and enumerate the reported locations

1. Group raw name-collated records before per-record eligibility filters. Scanner and writer use the same producer
   contract. Retain raw records for output and reporting checks. Require one template per QNAME; reject conflicting
   read groups under a QNAME instead of pairing different templates. For this paired-end input contract, each
   eligible record must identify exactly one of R1/R2 and be flagged paired; ambiguous mate-role flags are input
   errors. One record cannot supply both mates. Reconcile available query identity/length evidence within each mate
   without requiring secondary records to repeat SEQ.
2. Apply mapping/QC/duplicate filters once, with reasons. Filtered or unmapped records cannot supply a mapped mate.
   A secondary flagged mate-unmapped may participate if an actual compatible other-mate location was reported;
   metadata does not override the observed geometry.
3. Deduplicate locations by canonical query-to-reference path and relevant evidence, ignoring redundant primary,
   secondary and HI labels. Retain original records as output aliases. Duplicate representation cannot increase NM,
   multiplicity or probability. Conflicting scoring evidence on otherwise duplicate representations gets an input
   diagnostic, not an order-dependent choice of score.
4. Enumerate R1 × R2 without consuming locations, for HI and no-HI input alike. HI records the aligner's reported
   associations; it does not forbid another physically compatible combination among reported locations. Validate
   every combination and deduplicate accepted pairs.
5. Sort accepted pairs by canonical geometry/path identity. Record order, HI numbering and primary designation of
   equivalent records must not change this order. Keep the existing second-pass random-number algorithm.

Additional cross-HI combinations can therefore be accepted on STAR input. Even direct candidate changes are not
bounded by the orphan census. Check this enumeration against an independent exhaustive oracle.

### Validate the molecule before merging mates

Require a common reference and opposite genomic orientations. For the ordinary inward-facing library, let `a` and
`b` be the complete molecule's forward and reverse endpoint boundaries under the candidate path hypothesis, in
half-open coordinates. For unclipped reads they come directly from the aligned 5′ anchors; for supported terminal
soft clips, infer them on the candidate path as specified below **before** testing orientation/containment.
Require `a < b` and both aligned paths inside `[a,b)`. For unclipped contiguous alignments this includes
`forward.start <= reverse.start` and `forward.end <= reverse.end`. Complete overlap is allowed; no minimum insert
of one read length is imposed. Either R1 or R2 can be forward. RNA strandedness is a different check.

Validate splice consistency using the separate CIGARs. A base aligned by one mate cannot lie in an interval removed
from the common path. Overlapping nonidentical splice skips cannot both describe that linear path. Distinguish
aligned bases, insertions, deletions and reference skips; test their projection explicitly. Do not silently turn
a deletion into a splice or invent an allele to reconcile contradictory reported paths. Preserve raw splice
geometry when a later classifier blacklists a junction; removing its label must not conceal a mate contradiction.

Construct witnesses from the reference path, mutually consistent observed splice paths, and compatible
annotation-supported paths through **unsequenced** gaps. An unannotated observed splice does not need a shared
annotated transcript to form a valid paired RNA placement. Conversely, annotation overlap or a positive resolver
length is insufficient. Do not infer an intron across a mate's aligned bases or combine mutually exclusive
transcript paths. Deduplicate equivalent gap paths while preserving their supporting transcript identities.

For an unclipped monotone reference-coordinate path, or after endpoints have been projected back to genomic
coordinates with the correct candidate-path continuation:

```
L(path) = b - a - length of the union of introns actually removed by that path
```

Count a removed interval once even if both mates observe it. A genomic path uses `L = b-a`. RNA can have a short
molecular length across a long genomic intron. Verify projection and length by enumerating retained reference bases
on tiny examples, independently of `tx_frag_length`. This retains reference-path length conventions; it does not
newly infer sample-specific indel alleles.

Require `1 <= L(path) <= scan.max_frag_length` for every retained path. Keep a pair if at least one supported path
survives. Remove invalid transcript/gap interpretations before scoring or deposition: one short transcript cannot
admit all of that pair's longer interpretations. A 250-bp RNA path spanning 100 kb can survive while its 100-kb
gDNA interpretation cannot. Unknown/nonpositive lengths must not become neutral likelihoods. Existing source/splice
eligibility can further restrict support; this fix never turns a spliced record into a newly aligned contiguous
alternative. Do not use `detect_chimera` as the validator.

Keep pair validation distinct from the existing source-emission policy. In particular, an otherwise compatible
pair classified `ARTIFACT` currently has a legacy gDNA interpretation. Preserve that classifier/emission rule in
this experiment, subject to the new genomic-length bound; do not silently remove or redesign it while fixing
pairing. It is not a newly discovered contiguous alignment, and this plan does not certify its biological
interpretation. A cross-mate path contradiction still rejects the reported pair before that classifier runs.
Include compatible artifact-classified controls so the pairing A/B does not conceal a second splice-model change.

### Revised 5′-soft-clip policy: a conditional endpoint model

The first-3M-fragment census found 126,556 VCaP primary pairs with a 5′ clip (4.219%), including a reported 4.7%
of RNA-library fragments, and 44,672 minimap2 panel pairs (1.489%). This prices a substantial selection effect.
It is not an exact fragment-loss estimate: alternative placements, eligibility filters and overlap with other
rejection reasons remain to be evaluated. The census checks whether a primary mate has any splice, not its clip's
distance to a junction. Its length summary uses the maximum clip of the mates and reads only the outermost CIGAR
operation; refine this for adjacent hard/soft operations. First-3M counts do not establish whole-library proportions.

The owner confirmed the intended tool scope: instruct users to remove non-template 5′ sequence and discourage
trimming the insert's own 5′ bases. Put this in the input manual and provenance. With unrecorded template trimming,
the original insert endpoints cannot be recovered from this BAM; do not silently call the remaining interval the
original insert. If retained offsets are supported later, apply them once, distinctly from alignment soft clips.
This policy does not require restoration of 3′ adapter read-through to molecule length.

Replace blanket exclusion with a model for ordinary self-contained terminal soft clips. Under the model, a read's
input start is the insert endpoint, and each clipped query base consumes one template base on a specified candidate
path. Substitution errors can occur; unmodelled sequencing insertions/deletions, retained non-template prefixes,
unrecorded 5′ template trimming and incomplete chains violate that simple correspondence. Such inputs need their
own endpoint model/provenance, not an assertion that all clips have known template lengths. SAM soft clipping
consumes query but not reference; it does not identify the reason or reference extent of a clipped prefix.
See the [SAM CIGAR definitions](https://samtools.github.io/hts-specs/SAMv1.pdf).

For each candidate path `P`, use a path coordinate increasing from the forward-aligned end toward the
reverse-aligned end, regardless of transcript strand. Let `u` be the forward mate's aligned 5′ boundary, `v` the
reverse mate's aligned 5′ boundary, and `cF`, `cR` their templated 5′ soft-clip lengths. Then:

```
start_P = u - cF
end_P   = v + cR
length_P = end_P - start_P = (v - u) + cF + cR
```

`v-u` is a signed anchor separation on that path, not the union span of aligned blocks, not TLEN, and not an
absolute distance. The two aligned anchors may even cross after heavy clipping; the completed endpoints and both
aligned paths still have to pass containment. Unclipped outward pairs remain invalid.

An overlapping 100-bp insert illustrates the distinction: forward alignment `[105,200)` with a 5-base prefix clip,
reverse alignment `[100,195)` with a 5-base prefix clip. The aligned union spans 100 bp. Adding ten would give 110;
the endpoint calculation gives `(195-105)+5+5 = 100`. Even a clip on only one overlapping mate can be double-counted
by the union-span shortcut. Never implement this policy as `genomic_footprint += clipped_bases`.

Infer endpoints along the candidate path, checking its supported boundaries and requiring both aligned mate paths
inside the resulting molecule. Extension on a transcript may cross an annotated intron; subtracting a clip from a
genomic coordinate is not equivalent. A path that runs out of supported sequence cannot be rescued by a plausible
scalar length. Enumerate supported alternative continuations if needed; unresolved continuation remains unresolved.
Keep genomic and transcript endpoint hypotheses separate and apply the hard length bound to each completed path.

Store two distinct kinds of geometry:

- **Observed alignment:** original aligned blocks, query coverage, NM and observed junctions. A clipped base does
  not become matched coverage or junction support. No new reference alignment or CIGAR is invented.
- **Hypothesized molecule:** path-specific inferred endpoints and length. Calibration fragment occupancy,
  deferred-path handling, fragment-length scoring, start-dependent capture/coverage weights and opportunity must
  use this geometry consistently, while preserving the observed/implied distinction.

Updating length alone while leaving `coverage_weight(genomic_start, flen)` at the old start is inconsistent.
Use path-specific molecule coordinates where required, rather than overloading the existing observed-start field.
Do not add a junction crossed only by inferred extension to observed-junction evidence. The same endpoint model
must reach scan, calibration, scoring and annotation without changing the raw BAM alignment.

This is a length/geometry repair, not a clip-emission repair. NM excludes soft clips under the
[SAM tag specification](https://samtools.github.io/hts-specs/SAMtags.pdf), and the current scorer's NM term does not
score their bases. Retaining them in molecule length does not establish sequence support for that placement.
Preserve their metadata for the full splice-artifact likelihood, which must score or marginalize omitted query
sequence when comparing differently clipped alternatives. Do not add an arbitrary per-clip penalty or treat each
soft-clipped base as a confirmed mismatch. Source-emission redesign remains a separately measured mechanism.

For a more general endpoint model, the consumed template distances `dF,dR` are latent and the likelihood sums over
them and their compatible paths. This bounded model chooses `dF=cF`, `dR=cR`; it does not prove those equalities
from the data. Record that assumption in the input/analysis contract. A protocol with extra 5′ sequence should
remove/document it before this model is applied; unknown preprocessing cannot justify a claim of exact insert length.

The independent [`five_prime_clip_audit.py`](paired_fragment_review_evidence/five_prime_clip_audit.py) enumerates
98,692 physical molecule/read/clip combinations under the model. The anchor calculation recovers all lengths; the
outer-span-plus-clips shortcut fails in 9,075. Additional counterexamples cover a transcript boundary, extension
across an intron, non-template prefixes and sequencing indels in the clipped prefix. These verify conditional
geometry, not the proportion of real clips satisfying the model or a production implementation.

### Remaining unresolved endpoints and chains

- Support self-contained mate locations with aligned endpoints or the conditional soft-clip endpoints above,
  including ordinary 3′ clipping. Parse both CIGAR ends correctly under reverse alignment. Keep terminal S and H
  lengths distinct, including adjacent H/S operations; do not read only the first or last operation.
- Keep hard-clipped or otherwise unresolvable endpoints withheld unless original query coordinates and appropriate
  provenance have been recovered. The soft-clip revision does not authorize supplementary-chain reconstruction.
- Never concatenate supplementary records onto primaries by QNAME/HI. Withhold a location requiring chain assembly
  as `unsupported_alignment_chain`. Another self-contained location may still form a valid pair.
- Use SA/query intervals to establish chain dependence. Supplementary coverage outside a primary's aligned query
  portion prevents certifying that clipped primary as self-contained. If association is ambiguous, withhold the
  affected incomplete location rather than guessing. An unrelated fully represented alternative remains usable.
- Count rejected locations even when the fragment is placed elsewhere. No supplementary segment supplies a
  standalone mate location or contributes its NM twice.

These molecules may exist; their unresolved representations are outside this fix's supported subset. Price the
remaining loss after supported soft clips and all alternative pairs have been considered. Supplementary
reconstruction, new mate alignment and the broader endpoint/error model remain deferred.

## Counts, reporting completeness and filtering

`num_hits` counts distinct accepted paired placements. Never add `NH - observed_HI`. Maintain reporting status
separately: established incompleteness, unknown, or no incompleteness detected under a verified aligner contract.

Check reporting on raw groups **before Rigel filters**, so dropping duplicates, QC failures or incomplete HI groups
cannot mimic aligner withholding. For supported STAR HI/NH conventions, reconcile distinct HI groups with consistent
declared counts, including one-mate groups. A missing declared group indicates missing output; it does not certify
that the withheld alternative was a complete pair.

On established truncation, stop with the QNAME, aligner/profile, observed versus declared information, and a remedy
(for STAR, regenerate without the output cap). Do not leave outputs advertised as completed quantification.
Inconsistent tags get an input diagnostic. Generic NH without verified paired-reporting semantics is not an exact
paired count; record that completeness is unknown. Absence of truncation evidence cannot prove an aligner searched
every origin. This policy concerns missing declared output, not a claim of exhaustive alignment search.

| Reported input | Outcome |
|---|---|
| Valid pair plus incompatible orphan secondary | One placement; orphan unused |
| Three complete HI groups | Count distinct validated combinations, not automatically three |
| Verified STAR NH 3, one HI group present | Stop for incomplete reporting; no phantom placements |
| NH 2, one complete and one one-mate HI group | One placement if the extra mate has no compatible partner |
| Pairs at two distant copies, no cross-pair support | Two placements |
| Multiple transcript paths for one genomic pair | One placement with several path interpretations |

Apply `--no-include-multimap` after assembly: retain exactly one accepted placement and exclude more than one.
Secondary flags alone never decide. Check reporting errors regardless of this option.

## Consumer, accounting and annotated-BAM contract

Scanner and writer share grouping, eligibility, config, normalized locations, validation, ordering and ID allocation.
Assign IDs consistently before outcome-dependent skips. Test excluded and unpaired groups interleaved before a
known annotated fragment. Pass the configured length limit to both consumers explicitly.

Only accepted placements and valid paths enter ordinary inference. Apply path support to initial calibration offers,
deferred draws, RNA lengths and genomic lengths. Absent support stays absent after smoothing, likelihood floors and
tail extrapolation. Finite-support length laws and opportunity consumers use the same bound; check normalization on
finite examples. Do not change the contrast estimator or smoothing to force agreement with the old oracle.

Existing intergenic routing and artifact hold rules remain separate downstream policies. This work does not claim
that all accepted intergenic alternatives reach today's EM or that the current multimapper likelihood is correct.

### Spliced placements without a supported transcript

Owner scope clarification, 2026-10-15: Rigel estimates annotated transcript abundances and derives ordinary gene
abundance from those estimates. It does not assemble or estimate unannotated transcripts. A future estimator could
use fragment paths compatible with one gene's isoforms collectively, but that gene-centric mechanism is deferred.

A complete, physically compatible spliced placement may therefore be **placed, no supported transcript**. Retain
its alignment/path evidence and report it separately. Do not force it onto a nearby annotated transcript, a
synthetic unspliced RNA entity, or gDNA merely because its transcript candidate set is empty. An inferred 300-bp
observed path does not make an incompatible annotated target valid.

This applies to within-gene unmodeled splice paths, paths joining genes, and truly intergenic spliced placements.
Report their genomic/gene context as metadata where available; gene association is not an abundance assignment.
Do not add these diagnostic fragment counts to the ordinary gene table, change its transcript-derived totals, or
mix them into its abundance normalization. Reporting status is separate from the terminal placed/unpaired outcome
and from biological source: these are not another population.

Apply target support per placement. Preserve independently supported transcript or gDNA alternatives, and report
the unsupported alternatives as such. A fragment with no modeled target on any accepted placement remains
unassigned and is counted once. Mixed supported/unsupported cases need explicit accounting and a regression
fixture: omitting an unmodeled path does not supply its likelihood or establish posterior certainty over all
possible origins. The full splice-artifact likelihood remains separate work.

Use an observed-splice diagnostic, not a claim of certain RNA origin: the reported alignment can itself be an
artifact. This reporting decision does not authorize novel-only fragments as RNA-only calibration anchors or alter
the separately reviewed artifact eligibility rule. Calibration use must have an explicit evidence-based rule;
an empty transcript set is never sufficient evidence for gDNA. Annotated-BAM fallback must preserve the unassigned
status and must not invent a transcript/gene assignment or a certain source posterior.

### Fragment outcomes and annotation

| Terminal fragment outcome | Meaning | Ordinary calibration / EM |
|---|---|---|
| Placed | At least one accepted pair; passes multimap option | Existing routing, with validated path support |
| Unpaired | Eligible mapped evidence but no accepted ordinary pair or retained reported-chimera outcome | Neither |
| Reported chimeric | No ordinary pair; explicit complete aligner-associated chimeric explanation | Diagnostic only |
| Excluded | No eligible mapped evidence, or valid multimapper excluded by user option | Neither |

A missing mate, cutoff failure, unsupported chain or incompatible arbitrary cross-product is not proof of a
biological chimera. Retain operational chimera reporting for explicit complete associations, but do not promote
new geometry failures into confident biological labels. Rejection reasons distinguish missing eligible mate,
incompatible geometry/path, out-of-range support, unknown endpoint and unsupported chain.

The input ledger is:

```
templates present in BAM = placed + unpaired + reported chimeric + excluded
placed = genic + intergenic
excluded = no eligible mapped evidence + multimap disabled
```

Split no-eligible evidence into wholly unmapped versus QC/duplicate-only where possible. Genic/intergenic status
uses the accepted set. The simulator's molecules absent from the BAM belong in a separate truth ledger.

Add explicit fragment and alignment counters, update consumers, and document changed denominators. `unpaired` is
a fragment outcome; leftover records of a placed fragment are `unused_alignment`. Count a raw record once even if
it belongs to several placements. Keep raw mate-unmapped flags distinct from fragments missing an eligible mate;
do not silently change the units of the old `mate_unmapped` counter.

Annotated BAM requirements:

- Write every raw record exactly once, including unmapped/filtered records. Preserve aligner tags and flags.
- Compute tags from all accepted memberships. First-containing-hit logic must not hide a shared record's support
  for the assigned component in a later pair.
- Use `ZC:unpaired` for an unpaired fragment and `ZC:unused_alignment` for unused alternatives of a placed fragment.
  Excluded records have an explicit nonassigned outcome.
- Nonassigned records have `ZF=0`, `ZW=0`, `ZH=0`, `ZN=0` and cleared assignment IDs/labels. Clear stale assignment
  tags on reannotation, including write-through records; recompute per-record diagnostics separately.
- Define `ZH` as support for the assigned component through an accepted placement, not a unique physical winner
  when the inference selected no such winner. Carry that support into annotation; winning transcript ID alone is
  insufficient. A record shared by multiple pairs is still emitted once.
- The no-EM-result fallback may classify only an accepted placement under existing intergenic routing. It never
  stamps an unused record or an unsupported spliced placement as certain intergenic gDNA or certain RNA. The
  no-supported-transcript status above remains unassigned.

## Falsification and invariant tests

Build an independent tiny-reference oracle projecting individual CIGARs onto explicit molecules and enumerating
R1 × R2 pairs. It must not use the production resolver to choose expected answers. Compare accepted sets and lengths,
not merely final count tables.

| Family | Required cases |
|---|---|
| Original defect | Previous six fixtures; missing primary mate, incomplete HI group, mate-unmapped secondary and unmatched secondary; complete primary and remote-pair controls |
| Identity/enumeration | One mate reused in several pairs; nearby/distant copies; duplicates; primary/secondary relabeling; permutations; equivalent HI/no-HI representations |
| Endpoints/orientation | Both R1 orientations; same/opposite strands; cross-reference; inward/outward; valid overlap and protrusion; genuine 50-bp clipped and 100-bp pairs |
| Path consistency | Spliced/exonic pair; intron versus mate read-through; conflicting skips; repeated observed junction counted once; insertion/deletion projection; negative transcript; unannotated observed splice; compatible legacy-artifact controls |
| Length support | Exact bound and bound+1; long shared single exon; long intron/short RNA; RNA-only acceptance with gDNA absent; invalid isoform beside valid one; nondefault bound in scan and writer |
| Clipped endpoints | Both 5′ orientations; one/both clipped overlapping mates; crossed aligned anchors with a valid completed molecule; transcript-junction continuation; exhausted path; bound and bound+1 after extension; no invented observed coverage/junctions; clip-free numerical identity |
| Unresolved locations | Separate S/H and adjacent H/S; non-template-prefix/trim/indel counterexamples to unconditional exactness; 3′ controls; missing/ambiguous SA; overlapping query chains; complete alternative beside unsupported chain; no double NM |
| Reporting/filtering | Complete and truncated STAR reporting; inconsistent/mixed/missing tags; NH absent; one filtered mate; duplicate flags; multimap option after assembly |
| Accounting/annotation | One terminal template outcome; each record once; shared record in losing/winning pairs; no orphan gDNA; stale-tag removal; ID synchronization through exclusions |
| Determinism/capacity | Stable order under permutations; pinned reproducibility; existing parallel-order gates; high multiplicity; checked narrowing; no silent cap/overflow |
| Cache | Old assembly cache rejected even with unchanged payload/deposit rules; repaired-build round trip accepted |

Record baseline failures first; valid controls should pass. After repair, deliberately restore each singleton route,
old hit counting, weak facing, shared-transcript bypass, asserted-pair bypass, contradictory-path acceptance, stale
cache acceptance and BAM fallback separately. Each matching gate must fail. Save mutation/failure evidence and revert
each mutation before continuing. Do not import sandbox helpers into production tests.

## Execution sequence

### Freeze baseline and measure exact scope

Use an isolated worktree and native build from the reviewed SHA. Record source/diff, compiled extension identity,
dependencies, input/index hashes, aligner provenance, config and seeds. Do not reuse either intergenic prototype's
wheel. Re-record the same-input baseline.

First census the raw input conventions: flags, clipping, supplementary records, HI/NH consistency, reported pair
associations and record multiplicity. Do not label approximate mate matches as final accepted placements. After the
reference tests and native assembler exist, obtain the exact transition census from their validated output before
the full A/B: unique, multimapper, unpaired, chimeric and excluded, by rejection reason and provenance. Include
new cross-HI pairs and beyond-bound paths. That census replaces the approximate 501k/318 predictions.

**Gate:** all inputs satisfy the supported reporting contract, counts have explicit units/denominators, and the
standing decisions plus the reopened clip model's assumptions are frozen in fixtures. Record the library/5′
preprocessing provenance or explicitly declared analysis assumption; do not infer it from CIGAR S. A potentially truncated input is investigated before it
can supply a completed quantification comparison.

### Reference tests, then staged native implementation

Create the falsification tests and independent oracle first; save failures on baseline. Implement the validated
representation and shared scanner/writer assembly in the worktree. Keep temporary experimental snapshots isolating:

1. Complete-pair enforcement and accepted-placement counting.
2. Exhaustive construction and physical/path checks on constructed pairs.
3. Identical physical/path checks on aligner-associated pairs.
4. Supported 5′ soft clips through path-specific inferred endpoints, occupancy and lengths; compare with blanket
   withholding on the same input and report the assumption and exactly which fragments/path supports return.
5. New hard acceptance support on all retained interpretations and consumers.

These snapshots are diagnostic arms, not shipped modes. Final code has one validator and no singleton/legacy
asserted-pair fallback. Run small controls and the worst affected aligned condition between snapshots; measure
asserted-pair changes on VCaP separately. Do not rerun the full ladder for every intermediate snapshot without need.
No length-estimator, capture or mean/sum change enters an arm.

**Gate:** independent enumeration agrees exactly on sets and lengths, every required mutation gate fires, scan/write
parity and conservation hold, and every changed small-fixture result has a candidate/support explanation. Save the
exact whole-input candidate/fate transition census now; use it to identify which benchmark identity conditions hold.

### Integrate reporting, persistence and capacity

Update stats, CLI summaries, report model, annotation documentation, simulation/accounting consumers and their tests
together. Audit `num_hits` uint16 storage: check before narrowing and fail clearly if unrepresentable, or widen it
consistently through native/Python/serialization. Never silently wrap or cap candidate counts.

A default aligner setting does not bound arbitrary BAM input. Enumerate lazily; reject incompatible references and
orientations early. Prove any coordinate/path indexing optimization equal to exhaustive enumeration. Test high
multiplicity and measure pairs examined/accepted, memory and time. Add no arbitrary candidate cap.

Extend the cache behavior signature to cover native assembly and relevant config, using deterministic assembly
probes consistent with the existing behavioral-digest approach. Missing/mismatched signatures must fail loading with
a rebuild instruction. `--force` is an operational rebuild, not the correctness mechanism. Audit persisted buffer
artifacts too. Merely adding a payload field is not a proof that upstream behavior is covered.

**Gate:** stale caches rejected, repaired-cache round trip passes, consumers agree, no overflow/order dependence.
Run the full suite and `scripts/design/preflight.py --full`; record actual collected count and complete failure set.
Read golden diffs before regeneration. Name unrelated pre-existing failures separately and resolve them before
claiming a clean release gate.

### Final same-input validation

Pin `scan.total_threads=1` and `OMP_NUM_THREADS=1`; use `em.assignment_mode=fractional` and fixed seeds in both arms.
Build caches under each implementation's own signature. Compare final prototype against the freshly recorded
baseline, retaining the staged contrasts for attribution.

- **Oracle ladder and test chromosome:** preflight actual records against all new conditions, including paths,
  clipping and length. Where accepted sets, order and support are unchanged, require bit-identical numerical
  payloads/tables. Extra diagnostics/provenance need not make whole files byte-identical. Separately explain any
  control subset outside the approved support. Having no secondaries alone is not an identity proof.
  The sampled oracle includes 47,097 unspliced primary pairs with genomic spans above 1,000 bp (1.570% of 3M).
  Validate the reported intron-spanning RNA explanation per path. Surviving RNA support does not imply unchanged
  quantification when long gDNA support disappears. The clip-only stage must be a numerical no-op for clip-free
  input; the hard-bound stage is not required to preserve the old oracle output on this affected subset.
- **All five minimap2 conditions:** measure exact fates, paired provenance of remaining 100-bp observations, pre/post
  drain pools, length laws, calibration versus oracle, both RNA/gDNA confusion directions, and transcript/gene error
  by truth origin and stratum. Include gDNA-zero. Keep all-simulated, BAM-observed and accepted-fragment truth apart.
- **Length audit:** rerun the existing contrast audit; retain internal density, purities, de-tilted pools, signed
  contrast, negative mass and final normalization. Compare oracle truth matched to the accepted population. Bin 100
  need not vanish or match the old 0.24% oracle fraction. Vary read length independently of true insert length.
- **VCaP:** inspect exact candidate/fate changes, then all fitted models and outputs. Small direct changes can
  propagate globally through calibration; deferred-bank changes can also change later random draws. The 318 count
  is no output bound. Use frozen-model replay only as a diagnostic to separate direct effects from refitting.
  Profile back-to-back baseline/prototype runs.
  For the clip stage, separate soft/hard clips, source library, splice class, each mate's clip length, alternatives,
  residual chain dependence, path-boundary failures and maximum-length crossings. Report actual accepted-fragment
  changes, not only clipped-primary prevalence. Measure clip-to-junction distances before claiming junction-local
  enrichment; avoid treating a library label as molecular RNA/gDNA truth.
- **Strata and length controls:** judge unstranded/OFF, stranded/OFF and stranded/ON separately; report deferred
  unstranded/ON. Include both length-gap directions beside equal-length controls when the accepted length population
  changes, under the existing repository rules for length mechanisms.
- **STAR:** use existing STAR BAMs locally. Build/run the controlled panel STAR arm on SLURM after local tests and
  implementation; no Mac STAR installation is part of this task. VCaP is not a replacement for that comparison.

**Gate:** zero incomplete ordinary placements or invalid retained paths; exact template/record accounting; valid
controls retained; numerical identity where its preconditions hold; every material aligned-accuracy regression
explained before promotion. Use existing success criteria, justifying any new tolerance before examining outcomes.
Do not tune acceptance to stabilize real-data tables or declare structural correctness solely from a headline gain.

### Review, promote and resume parked designs

Deliver patch, manifest, independent-oracle results, mutation evidence, transition tables, timing/memory, length
audit and per-stratum comparisons for owner review. Only after that gate promote to the main source tree. Move
settled definitions to permanent documentation, update release notes, and let the owner drive the commit.

Pairing is complete when its invariant and integration pass these gates. **The release remains blocked** until
the repaired-input length law and controlled STAR/cluster validation are acceptable, along with the other release
requirements. A remaining pathological contrast resumes its parked design; it is not solved by declaring orphan
placements gone.

Resume mappability opportunity and the full splice-artifact design on this accepted-pair foundation. Their
denominators must reflect pairing, length and filtering, including observable repeated regions and two mapped ends
spanning an unsequenced gap. Removing all-N truth reads alone does not repair the abundance denominator.

## Readiness and residual uncertainty

The non-clip decisions stand. The conditional endpoint model above replaces the reviewer's earlier
blanket-withholding recommendation, and the owner has supplied the intended input scope: remove non-template
5′ sequence and preserve the insert's own 5′ bases. No further protocol-scope decision is needed to design the
scoped prototype. Verify or record the assumption for benchmark inputs; STAR local alignment alone does not
establish their preprocessing history. Do not advertise unconditional exactness for sequencing indels or inputs
outside that contract. If reporting truncation is established, regenerate the input rather than changing hit
arithmetic.

The unknowns are empirical: loss from unresolved endpoints/chains, added valid combinations, performance at high
multiplicity, and whether the length contrast recovers. They have explicit measurements above. A review cannot
guarantee a not-yet-written implementation; confidence in promotion requires the enumeration, break tests and
same-input comparisons. Estimate effort after the census/native prototype, not as a promised few hundred lines
and a few hours.
