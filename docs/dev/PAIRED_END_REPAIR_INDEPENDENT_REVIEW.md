# Independent review of the scan-stage paired-end repair

2026-10-16 · Packet dated 2026-10-16 · Baseline `41899ae4b8209893f7b5ae0cf893304e6796dabf`

## 1. Verdict

**Do not implement the complete design as written.** Preserve complete pairs, exhaustive pairing, overlap eligibility and the existential length witness.

- **Blocker 1:** tied alternatives encounter best-hit RNA merging and mean gDNA aggregation; native probes demonstrate lost RNA contributions and gDNA dilution (§§7.2, 7.4).
- **Blocker 2:** reconciliation confuses sequence preference with identification, and does not define a once-counted, common sequence likelihood (§§7.2, 7.8).
- **Blocker 3:** the proposed alternatives lack a valid calibration/training contract; holding all junction support can force its RNA reading to zero (§§7.2–7.8).
- **Blocker 4:** the claimed exact clip-containment walk rejects a valid 100-bp TA2 molecule when neither alignment observes its junction (§5).

Resolve these contracts before approval; implement the revised scanner components in separate arms. The regression diagnosis identifies major mechanisms, but overstates exclusivity and experimental isolation.

## 2. Findings

The packet was read first and completely. All **125 manifest files** match the archive hashes; all corresponding workspace files match, at the stated base commit. I ran the supplied baseline counterexamples and additional finite/scorer probes; results are in [baseline output](paired_end_repair_review_evidence/baseline_counterexamples.txt) and [review audit](paired_end_repair_review_evidence/review_audit.json), reproducible with [the review-only script](paired_end_repair_review_evidence/review_audit.py). Recorded library measurements were inspected, not rerun. No production code changed.

All recommendations retain gDNA, RNA+ and RNA− as the only populations, the per-locus EM with calibrated gDNA prior, and calibration's hard seeded draw. Reporting dispositions are not populations.

**A/B contract for every recommendation below:** vary that mechanism alone against its parent and the baseline, with `scan.total_threads=1`, `OMP_NUM_THREADS=1`, fractional EM and recorded calibration seeds. Freeze the declared gDNA laws and ruler in both arms, assert the actual arrays consumed, and score against truth. Report both artifact RNA mass and genuine RNA loss, plus transcript/gene/calibration errors, separately in the three required strata; report the deferred fourth. The tests below add mechanism-specific falsifiers.

**1 — blocker · §§7.2, 7.4, 11.4: a sequence tie does not validate the multimapper path.** *Derived; measured on the baseline scorer.* [scoring.cpp](../../src/rigel/native/scoring.cpp), `score_mm_alignment`, chooses the lowest-overhang hit per transcript, then the greatest likelihood; `flush_mm_group` averages gDNA-compatible **genic** hits. The audit gives RNA log-likelihood **0 → 0** for two equal contributions to the same transcript, where distinct equally weighted origins require a **ln 2** numerator gain. Adding a gDNA origin with relative likelihood **10⁻⁶** changes the existing term by **0.500000506**, rather than the sum's **1.000001**. A splice hit itself does not enter the gDNA denominator.

Define physical-origin identity and the normalized placement sum for both RNA and gDNA before activating alternatives, including intergenic support, capture and the corresponding opportunity/prior. A bare mean-to-sum edit is insufficient. **Test:** enumerate finite origins, duplicate only a report, add a negligible origin, and add an observable copy. Match enumerated probabilities and conserved mass; change neither probability nor opportunity for a duplicate representation. Repeat with two readings sharing one transcript and unequal overhangs.

**2 — blocker · §7.2 rules 1–3, §§7.4, 7.8: preference, ambiguity and cancellation are conflated.** *Derived.* `Δm ≤ 0` includes improvements; only **Δm = 0** gives factor one. Even a tie for the moved bases does not prove S is the molecule: C may also be valid. Likewise, one surviving local reading's factor cancels only if it is common to **all** competing explanations of that template, including remote placements. The rule does not establish that condition. If corrected NM also changes, multiplying by the proposed factor counts the same mismatch difference again.

Keep the owner's within-K correction as construction of a complete S reading; do not mistake it for proof excluding C. Specify both anchoring directions, joint multi-junction edits and a common sequence reference across all placements. **Test:** use S/C ties, `Δm = ±1`, a valid remote pair, and an alternative reachable from two different reported parents. The same final reading must receive the same score regardless of its generating parent or report order. Changing one discriminating base must change odds by exactly the derived factor once.

**3 — blocker · §§7.2–7.8: adding the alternatives changes calibration before EM.** *Derived from current code.* [second_pass.py](../../src/rigel/second_pass.py), `placement_abundance`, reads annotated-junction support for S and contiguous traffic for C. When every S-supporting template is held, S can have zero independent support while C has positive support; `combine_factors` makes that a hard zero. A tie factor cannot restore S. Also, `score_held_fragments` chooses its length PMF from **implied** introns, whereas observed introns are separate. [accumulator.cpp](../../src/rigel/native/calibration/accumulator.cpp), `deposit_survivor`, treats a path as spliced only when an annotated junction resolves; a wholly novel observed splice therefore enters different channels from its EM interpretation.

Specify shared candidate support, eligibility, training selection and unavailable-support semantics before these additions. Preserve the hard seeded draw. **Test:** all-held versus independently supported junctions, annotated versus annotation-withheld RNA, motif present/absent, and mixed genuine/artifact junctions. Trace factors, pools and the single deposit through both passes. Holding the focal evidence must not manufacture decisive evidence against its own RNA path.

**4 — blocker · §5.2, §8: the observed clip walk is only conditionally exact.** *Derived; independently enumerated.* On TA2, take forward **`50S50M@2000`** and reverse **`50M50S@1150`**, with no observed N. The packet's continuous completion gives **a = 1950, b = 1250**, rejecting containment. Yet both 100-base reads fit the same 100-base molecule `[1150,1200) + [2000,2050)`, with **a = 1150, b = 2050**. Neither alignment supplies the missing junction; the reverse mate does not determine the entire walk to the forward anchor.

Permit coherent annotation-supported clip continuations in the acceptance witness, without imposing transcript endpoints or filtering candidates. Where continuation is unavailable, report unresolved geometry rather than claiming exact impossibility. **Test:** enumerate both endpoints on the toy path, mirror strands and swap mate identities; accept this case and the packet's observed-skip example. Deliberately restrict walks to observed skips: this new gate must fail.

**5 — fix before implementing · §§4–5 identity: “identical geometry” is underspecified.** *Derived; historical measurements.* Blocks joined across D omit the query mapping needed to compare bases. `5S95M` and `95M5S` at one start have identical blocks but different endpoints. Homopolymer-shifted insertions may be equivalent, but require an explicit equivalence rule. Opposing XS/ts evidence, missing NM, and different retained/raw skips cannot be resolved by selecting the smaller NM. The supplied stage-3 review measured duplicated R1 evidence changing **NM 2 → 4**, and reported **110** real duplicate conflicts; these are prototype findings, not evidence that the new design fixes them.

Retain query-to-reference CIGAR mapping, clip/query offsets, evidence validity and aliases. Recover bases before filtering their carrier record; orient and slice recovered sequence in query coordinates. **Test:** duplicate within/across HI, normalized indels, clip-distinct records, conflicting tags, missing NM/SEQ, hard-clipped carriers and supplementary-without-SA. Permutations must preserve canonical hypotheses; independent valid alternatives must survive withholding a chain.

**6 — fix before implementing · §6: the ledger is necessary but not sufficient.** *Derived from scanner code.* A “placed” template may have no EM row, no calibration deposit, or only intergenic/unsupported-splice support. Those dispositions need explicit counts. Withheld locations also need a state distinct from absent locations, or “missing mate” consumes the “unresolved location” reason. Define reason selection over **all combinations**, not the first failing pair, and reconcile §7.2's template-wide “unpaired” wording with §5's surviving-placement rule.

The claimed old raw-flag unit is wrong: `bam_scanner.cpp` increments `mate_unmapped` for incomplete **hit groups**. Preserve that meaning or explicitly rename it and add a raw-flag count. **Test:** a mixed BAM with all-filtered groups, one geometry failure beside a too-long consistent pair, an unresolved pair beside a valid alternative, empty-target placements and unused aliases. Reconcile raw templates, records, accepted placements, deposits, routed units and reasons; every record appears once with stale assignment tags cleared.

**7 — fix before implementing · §§2, 10: the diagnostic claims exceed their controls.** *Measured and derived.* The archive contains **no joint law-and-ruler injection arm**. Capture-OFF baseline ruler log10-factor SD is **0.33279** unstranded and **0.42150** stranded, contradicting “the ruler is flat here.” The wrapper also saves its **own, pre-substitution** ruler; that file is not necessarily the ruler applied. “Law fixed” changes only the uniform-frame calibration law.

**Test:** record/hash the arrays actually consumed, then cross calibration-law replacement with ruler replacement on a fixed build and repeat across builds. Separately fix the EM's realized law and second-pass inputs. If the purported isolated effect disappears, revise the attribution; details follow in §3.

**8 — fix before implementing · §7.3: define gene membership at the junction.** *Derived.* A gene's outer span is not proof that it supplies both splice flanks; nested/antisense genes expose the distinction. Conversely, one unannotated flank does not establish two different genes. Specify flank-associated gene sets, splice-strand compatibility, the no-common-gene condition, and treatment of annotation disagreement. Apply the owner's rule to unannotated inter-gene junctions at any span; do not extend it to ordinary unspliced cross-gene pairs.

**Test:** same-gene novel junction, two nearby genes, nested/antisense overlap, a read-through gene annotating both flanks, one unassigned flank, and a valid nonchimeric alternative. Require exact dispositions. Measure true RNA excluded by policy separately from false classifications of ordinary RNA.

**9 — consider later · §5 enumeration: optimize only after completeness is demonstrated.** *Derived.* Cartesian pairing costs `nR1 × nR2`; joint junction alternatives can grow further. **Test:** compare any reference/orientation/interval indexing optimization against exhaustive tiny-reference enumeration and back-to-back real-data timing. No hit cap or mismatch cutoff without a proven bound. Count-width overflow and persistence checks remain immediate correctness requirements.

## 3. The diagnosis (packet §2)

**The main mechanisms are supported; “the whole regression” is too strong.** Under capture ON, stage 1's own transcript error rises **564,547 → 801,015**, while gene error falls **342,837 → 341,628**. The reciprocal ruler swaps are persuasive interventions: stage 1 with baseline ruler scores **488,421 / 275,533**; baseline with stage-1 ruler scores **907,999 / 442,235**. Thus the ruler can reverse the build ranking. This demonstrates mediation, not that all other effects vanish or that their contributions add independently.

The contrast failure is concrete. From the recorded exact values, stage 1 uses **63.820 f0 − 62.820 f1**: its signed mean is **739.77 bp**, negative mass **2.47886**, and final fitted mean **317.80 bp**. Stage 3 uses **52.168 f0 − 51.168 f1**, with signed mean **668.03**, negative mass **2.03638**, final mean **317.17**. The 318-bp result therefore includes projection, blending and smoothing; it is not simply the raw linear inversion. Near-purity on the oracle makes the contrast effectively one-pool, despite its also-small purity separation. Aligned input loses that protection. The common-contaminant assumption and observability remain possible failures alongside conditioning.

The stage-3 ON internal density **0.000240** versus oracle **0.001820** supports a separate internal-fit defect. Agreement of final calibration densities (**0.0673 / 0.0687**) does not validate that internal estimate.

The law intervention strongly supports the calibration route: stage 1 improves **801,015 → 577,866**, and the centred-log ruler discrepancy improves **0.101 → 0.045**. But `stage_run.py` replaces only `gdna_pmf` on the `run_pipeline` call. It leaves the second-pass fit and EM's `gdna_realized_pmf` unchanged; the RNA law also remains build-specific. Calibration density, prior and ruler may all respond. This is a valid **single-input intervention**, not a fixed complete length model or an isolated ruler mechanism. Ruler substitution is likewise valid, but replaces full transcript effective lengths, including scale and RNA-length effects, while leaving gDNA's effective length/prior alone.

At zero gDNA, flattening the ruler removes **93.79%** of baseline transcript error and **92.90%** of stage-3 error. False capture is therefore the dominant demonstrated mediator here. It does not establish whether the root is false calibrated DNA, normalization near zero bulk, or both. The flat-ruler residual increases by **1,256 transcripts / 6,566 genes** in absolute-error units. Attributing that entire difference to **31,848 unpaired templates** requires a matched transition decomposition; altered candidates, fitted laws and seeded draws can affect surviving templates too.

The separating experiments are specific:

- **Law versus ruler versus prior:** within-build law × ruler interventions, then fixed realized-EM-law and fixed-prior comparisons; inspect isoform shape and RNA:gDNA scale separately. These are diagnostic injections, never release defaults.
- **Internal density cause:** hold accepted placements fixed; enumerate paired observation opportunity, then compare the fitted rate and truth-conditioned pool shapes. An N-only mask is insufficient: two observable ends can span an unsequenced N gap, and non-N repeats can be unobservable.
- **Unpaired loss versus changed inference:** compare common accepted template IDs with frozen nuisance inputs, and account separately for each lost template's truth and independent alternative. Repeat calibration seeds to distinguish sampling variation.

The OFF law-fixed improvements are real recorded results (**310,430 → 266,983** unstranded; **309,207 → 249,401** stranded), but neither proves performance with both downstream defects held fixed. Nor do four conditions establish all ladder rungs.

## 4. The components

**Locations and pairing (§§4–5).** Same-reference, opposite orientation and inward containment are correct for the declared linear FR model; either mate may be forward. Proper-pair, HI, TLEN and mate-coordinate assertions are not physical validators. Require valid R1/R2 identity, positive M/=/X support and reference/query bounds. Keep the input's templated-5′-clip assumption explicit: SAM clipping alone does not establish it. Reverse-aligned SEQ is stored reverse-complemented; secondary recovery must respect that and hard-clip offsets. [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

STAR truncation detection is useful only for verified tag/output semantics, before filters and deduplication. Check per-mate missing groups and inconsistent NH/HI as well as the template union. Missing tags mean unknown completeness. The appropriate STAR output remedy is `--outSAMmultNmax -1`, subject to its mapping filter; it cannot recover hits never found or discarded upstream. [STAR parameter definitions](https://raw.githubusercontent.com/alexdobin/STAR/2.7.11b/source/parametersDefault).

The **existential length witness is sound as an acceptance screen**. Its additional introns must lie in genuinely unobserved gaps, come from one coherent path and never erase aligned bases through K slack. Its transcript is a witness, not an eligibility filter. TA1's 550-bp witness must not replace TA2's 350-bp length or the synthetic span's 1,150-bp emission. Keep these explicitly approximate baseline lengths; a witness is not a proof that each scoring interpretation satisfies a hard physical bound.

The bit-identity promise needs a narrower domain: unchanged complete standalone record pairs, normalized evidence, splice tags, hit multiplicity and consumer route. Duplicate NM removal, chain withholding, splice-strand union changes and geometry ordering are deliberate exceptions. Identical per-hit emissions also do not imply identical post-draw calibration if ordering or the candidate universe changes.

**Unpaired templates (§6).** Retain exclusion from all fragment likelihoods and length pools. Add the disposition audit in finding 6, using explicit scanner outcomes rather than absence from the buffer. Preserve genomic/locus context when RNA targets are empty; an unsupported observed splice must not become gDNA by fallback. Report rates against both all raw templates and eligible templates, with units. Write unpaired, excluded, chimera and unused records consistently in annotation's second traversal.

**Reconciliation (§7.2).** Re-place bases; do not delete them. Recompute the union footprint/overhang, but preserve both mates' sequenced observations for sequence likelihoods. The packet's two corrections remove **4.605** and **6.908 nats** of overhang penalty—**100× and 1,000×** RNA-weight gains. Thus “already certified RNA, so costs nothing” is false for abundance, competing remote placements, pruning and training, even when a hard label stays unchanged. Geometry and sequence are distinct evidence; correcting geometry is not itself double-counting. Retaining the old overhang penalty, separately rewarding the same repaired matches, or updating NM and multiplying the same mismatch correction would be.

K is an approved local-repair bound, not a biologically derived maximum. The census's k=6–8 classes contain **437/466 = 93.78%** mRNA-origin templates OFF and **441/473 = 93.23%** ON. That supports testing longer corrections, not selecting a new cutoff. Origin does not verify the reported junction or repair; VCaP's **2,165 RNA-library / 5 DNA-library** observations are provenance, not 99.8% demonstrated repair accuracy. Recount after full pair validation, with all disagreement types retained rather than only the first class.

Construct S/C jointly across both mates and every junction, considering both local anchoring directions. Preserve partial-RNA alternatives; gDNA must explain every apparent splice contiguously. Resolve ambiguous shifted junctions without an arbitrary annotated winner when multiple paths remain sequence-equivalent. Revalidate containment after edits. The interim unequal-reading exclusion is defensible only as an explicitly measured incomplete development stage, not a release solution; it must not erase an independent accepted pair silently.

**Read-throughs and alternatives (§§7.3–7.5).** Any-span exclusion follows the owner's ruling once finding 8 defines the predicate. The **3,544 long-span cases** price prevalence, not the incremental short-span rule: baseline already rejects the long toy example. Report true read-through loss, ordinary RNA loss and artifact leakage separately. Ties require no *relative sequence factor*, but still require findings 1 and 3. The earlier source-A **2.5–3.3%** improvement must be remeasured on complete pairs. Carrying intergenic placements can be an inert schema step only if existing mean denominators and routing truly stay unchanged.

**Emission factor (§7.8).** For the same n observed A/C/G/T bases under independent symmetric substitutions, `e(m) = (1−ε)^(n−m)(ε/3)^m`, hence `e(m+Δm)/e(m) = [ε/(3(1−ε))]^Δm`. This derivation permits negative Δm. It requires equal explained observations and fixed indel/clipping treatment, not arbitrary differently clipped alignments. Use one common comparison convention; a relative factor anchored separately to each reported parent is not invariant.

If NM is updated from m to m+Δm and the new factor is also applied, the relative multiplier becomes `(α_mm β)^Δm`, not `β^Δm`, where `β=ε/[3(1−ε)]`. Specify which term owns those changed bases while preserving untouched emissions.

The proposed ε estimates **alignment-selected substitution divergence**, not sequencing error. Define the sampled complete pairs, count each sequenced observation once rather than once per report, and handle missing/inconsistent NM, ambiguity bases, empty training and ε=0 without an invented floor. SAM NM includes inserted/deleted bases and counts ambiguous symbols as differences. [SAM tags specification](https://samtools.github.io/hts-specs/SAMtags.pdf). Shared variants and overlapping mates limit the independent-error interpretation. Carry the same sequence evidence into both stages, with each stage's declared other factors; adding the multiplier alone does not settle these contracts.

## 5. Order and scope (§§10–11)

First close the four blocking derivations and canonical data/disposition schema, then prototype locations, pairing and unpaired reporting independently. Move the already-approved **reference store before active reconciliation**: doing an unchecked repair and later adding its required check creates avoidable, confounded intermediate behavior. Test its reference identity, offsets, N/ambiguity policy, masking, reverse complements, contig boundaries and missing-store failure.

Read-through classification can proceed in its own arm. Resolve placement aggregation/opportunity and calibration support **before enabling even tied alternatives**; these downstream dependencies cannot wait behind §7.4. Store alternatives inertly if necessary. Keep unequal alternatives inactive until the sequence contract is complete. Length-estimator and false-capture repairs can follow the core scanner under properly frozen diagnostics; genuine final abundance benefits require subsequent own-law/own-ruler measurements.

**Implementation gate:** revised equations and independent finite path/origin enumeration; explicit exceptions to identity; complete candidate and disposition schemas; tests that expose each claimed defect; intervention hooks proven to reach their named consumers. Do not reuse the withdrawn target-containment expectations.

**Release gate:** each isolated repair passes its failing/broken-again tests, full required suite and cache/persistence checks; then the complete system passes truth-scored aligned and oracle ladders, zero-gDNA controls, both length-gap directions, STAR and minimap2, VCaP transitions and independent artifact-containing validation. Read both origin errors per stratum, calibration/prior/ruler diagnostics, counts, runtime and memory. Retain the blacklist until its matched two-error comparison justifies removal. Observable opportunity, intergenic origin accounting, clip lengths and the reference-copy/allele scope remain release obligations under the owner's “fix everything” ruling. Indistinguishable sample alleles require an honest scope decision, not an invented discriminator.

§10's absolute “no stratum gets worse” gate is not yet met by its recorded zero-gDNA arm; a correct policy exclusion can also have a nonzero truth cost. Quantify such costs separately. Any unavoidable regression needs explicit acceptance; do not silently loosen a tolerance or hide it behind another mechanism.

## 6. Tests (§9)

These values were recomputed by interval arithmetic and agree with the executed baseline resolver:

| Case | Independent calculation | Expected exon_bp / overhang / L |
|---|---|---|
| One-base overlap | coverage 200; overlap 1; outer span `1001−700` | TA1, TA2, synthetic A: **1 / 199 / 301** |
| Both outer boundaries | TB1 overlap `50+50`; coverage 200; span 250 | **100 / 100 / 250** |
| One overhanging mate | TA1 overlap `50+100`; `L=500−300`; TA2's 800-bp intron is not wholly inside the span | TA1 **150 / 50 / 200**; TA2 **50 / 150 / 500**; synthetic A **200 / 0 / 500** |
| Long-intron witness | span 1,150; subtract TA1's `300+300`, or TA2's 800 | TA1 **200 / 0 / 550**; TA2 **200 / 0 / 350**; synthetic A **200 / 0 / 1,150** |
| Acceptor-side correction | old union `50+100=150`; corrected union `50+98=148`; span `948−800` | TA2 **148 / 2 / 148 → 148 / 0 / 148** |
| Donor-side correction | corrected union `97+50=147`; span `447−300` | TA1 **147 / 3 / 147 → 147 / 0 / 147** |

The short insert also gives **60/0/60**; antisense overlap gives **150/50/200 for all five candidates**. The counterexample script passes premerged blocks to the resolver: it verifies emissions, **not** record parsing, pairing, sequence restoration or calibration.

Missing gates include crossed clipped anchors (finding 4), R1-forward/R2-forward and minus-strand mirrors, adjacent H/S, indels crossing repair boundaries, incompatible overlapping skips, multiple annotated shifted alternatives, unknown bases, missing tags and chains without SA. Include exact K/K+1 and length-bound edges, raw blacklisted skips, an RNA witness whose gDNA span is long, and a gap intron that would erase aligned bases. Test multi-junction S/C combinations, same-transcript alternatives, remote competitors, all-held support, canonical draw ordering, buffer spill/reload and hit-count overflow.

Add deliberate breaks for asymmetric RNA/gDNA aggregation, double-applied NM, parent-dependent Δm, fabricated junction support, loss of raw skips, dropped sequence/evidence fields during serialization, and omitted filtered-template output. Distinguish **regression invariants** from **baseline defect demonstrations**: one-base overlap and existing long read-through rejection should already pass baseline. Only the added short-span read-through fixture can establish that removing the new rule breaks that repair.

## 7. The four open decisions (§12)

**1. Approve the design?** No, pending findings 1–4 and the specified contracts/tests. Approve the complete-pair architecture and overlap-preserving witness in principle; do not approve a production implementation or the supplied stage-3 prototype by inheritance.

**2. Read-throughs at any span?** Yes. A length threshold has no role in whether a known inter-gene junction is a read-through under the owner's policy. Define flank membership, preserve nonchimeric alternatives, and measure the incremental short-span cost separately from already-rejected long cases.

**3. Extend beyond K with the reference?** Yes for proposing and evaluating complete corrected readings. Do not turn `Δm ≤ 0` into automatic exclusive S selection, or tune a new k cutoff to the observed bump. The approximately 93–94% mRNA-origin bump motivates a falsification cohort; it does not identify the remaining templates or prove corrections.

**4. Include the emission factor now as straightforward?** Not as a standalone multiplier. Close the common-base/NM convention, aggregation, calibration support and training contracts now, because the alternatives depend on them. Implement that reviewed dependency in separate arms before activating alternatives; deferring it while shipping ties through the existing scorer would leave a known defect active.
