# Parked work: every open thread, set aside to be resumed (2026-10-15)

For the owner. On 2026-10-15 the owner paused the open threads so the single-mate bug can be fixed first
(`MATE_PAIRING_FIX.md`). This file lists each thread, where it stands, where its materials are and what resumes it.
Sandbox: when a thread settles, its result moves to its one home (`ISSUES`, `DESIGN`, `TRAPS`, `EQUATIONS`).

**Almost every thread waits on the pairing fix.** Existing aligned-input measurements used the defective assembler
and may be affected, so each thread below resumes by re-measuring on the fixed tree before deciding anything.

## 1. The release

| item | state | resumes with |
|---|---|---|
| The release sequence (owner, 2026-10-14) | Steps 1–5 done. Step 6 (release issues) is open: the EM thread fix landed and was pushed (f02851e5, 7805ede7). | Unchanged: step 7 follows once the blockers below are closed. |
| Local commits, not pushed | e7dc6f6f (the splice-artifact ISSUES and ROADMAP text) and 41899ae4 (the second pass's strand factor for a sequenced motif) | The owner's push. Before it, re-derive the `CLAUDE.md` suite baseline: 2,959 collected since 41899ae4, against the 2,958 recorded. |
| The owner's release actions | CI dispatched by hand; the bioconda recipe PR | `ISSUES: the-release-builds` |
| Release-blocking (owner) | Splice artifacts (2026-10-14); the length defects on aligned input (2026-10-15) | §3 and §5 |

## 2. The active item: complete, compatible paired placements

`MATE_PAIRING_FIX.md` now contains the reviewed execution plan. The owner approved common validation for asserted
and constructed pairs, rejection of established reporting truncation, and a configurable hard molecular-length
acceptance bound. After the census, the 5′-clip policy was reopened: the reviewer now recommends retaining terminal
soft clips through candidate-path endpoint inference, under the owner's stated scope of removing non-template
5′ sequence and preserving insert bases. Hard-clipped/unresolved endpoints and supplementary-chain reconstruction
remain separate. These changes are measured in separate prototype comparisons; no production implementation has
been made by the review.

Evidence: the six-case reproduction, length audit and new geometry probes (`paired_fragment_review_evidence/`).
The earlier scope counts (`$RIGEL_RUNS/prototypes/2026-10-15_mate_pairing/`) are approximate; an exact
candidate/fate census replaces predictions of how many multimappers become unique or how narrowly VCaP changes.
`five_prime_clip_audit.py` adds an independent conditional-geometry check: adding clips to the outer aligned span
is wrong for overlapping mates; use each mate's aligned 5′ anchor on the candidate path and keep inferred molecule
extent distinct from observed alignment evidence.

## 3. The length law on aligned input (release-blocking)

- **What is known.** On the minimap2 panel the gDNA length law put 76.3 % of its mass at 100 bp, the read length;
  the oracle BAM of the same reads gives 216.7 bp. VCaP (STAR) did not show that spike; this does not establish
  immunity to the pairing defects or the revised acceptance rule. The reviewer traced the mechanism
  (`INTERGENIC_GDNA_REVIEW.md`, finding 2):
  - calibration's second pass raises the intergenic 100-bp deposits from 8,218 to 147,420;
  - the length module's internal gDNA density is 0.000798, against final calibration's 0.0521;
  - its two-pool contrast gives a signed 100-bp bin of 5.08, with negative mass of −5.66;
  - clipping and renormalizing produce the 76 %.
- **Open:**
  1. how much the pairing fix alone repairs;
  2. why the internal density is 65 times too low (unobservable opportunity is a hypothesis, not yet isolated);
  3. the contrast's common-shape assumption, and its clipping of negative mass.
  No bin removal, read-length floor or extra smoothing: the reviewer and the owner refused those workarounds.
- **Materials:** `paired_fragment_review_evidence/length_contrast_audit.py`, `length_contrast.json` and
  `length_contrast_arrays.npz`; the oracle and aligned runs in `prototypes/2026-10-15_mm_intergenic/diag_fl/`.
- **Resumes:** rerun the reviewer's audit on the fixed tree, then decide whether the estimator needs its own design.

## 4. Observable opportunity (owner's ruling: abundance per mappable base)

- **What is known.** The index's region lengths include sequence no read can come from. N alone is 36 % of the panel
  genome's intergenic length and 11.4 % of GRCh38's (150 of 1,318 Mb); the genic share is 0.05 %. On the oracle
  ladder the simulator places gDNA reads in the N gaps (0 secondaries, 0 unmapped records in its BAMs), so the
  defect never showed. Aligned, those reads are lost: 0.99 M per capture-OFF `g50` condition.
- **The reviewer's correction to my proposal.** Requiring the whole insert to be N-free is wrong: two mapped ends can
  span an unsequenced gap. The right quantity is the opportunity for an accepted paired fragment. It includes
  multimappable origins and excludes unobservable ones, and is approximated with stated uncertainty at genome scale.
  Zero coverage is never a mask.
- **Open:** the paired-acceptance opportunity's definition and its approximation; a mapping-loss ledger for the
  simulator's all-N reads; what the oracle ladder does with them.
- **Materials:** `prototypes/2026-10-15_mate_pairing/n_in_regions.py` and its outputs; `INTERGENIC_GDNA_REVIEW.md`
  finding 4 and §4.

## 5. Splice artifacts (release-blocking)

- **State.** The design has been reviewed twice. The current specification is the reviewer's revision of
  `SPLICE_ARTIFACTS_IMPLEMENTATION.md`:
  - whole-pair candidates;
  - local, remote and allele sources;
  - one normalized likelihood;
  - no base-only hold rule;
  - allele identifiability stated as a limit.

  Its §7 lists six design questions that must be settled before any code. My step-back review
  (`SPLICE_ARTIFACTS_PLAN_REVIEW.md`) is superseded where the two disagree. It differs on the opportunity rule and
  the hold rule, and the reviewer's position governs.
- **Prototypes, not approved to merge:** source A step 1, the local tie as a second EM hit (worktree
  `<worktree>`; `prototypes/2026-10-15_splice_a/`).
  - On VCaP it moved 2,044 of 2,363 DNA tie artifacts to gDNA.
  - Its isoform redistribution side effect is unjudged.
  - It places both mates, but it must be re-measured on the fixed tree.
- **Owner decisions pending:**
  - the hold rule: the reviewer recommends removing the base-only rule I proposed and you adopted;
  - source C's scope, including non-reference retrocopies;
  - an external review of the complete design.
- **Materials:** `SPLICE_ARTIFACTS_DESIGN.md` (packet 1 and its §12 verdicts), `SPLICE_ARTIFACTS_REVIEW.md`
  (review 1), `SPLICE_ARTIFACTS_DERIVATION.md` (with the reviewer's status note), `SPLICE_ARTIFACTS_IMPLEMENTATION.md`
  (the current specification) and `SPLICE_ARTIFACTS_PLAN_REVIEW.md`.
- **Resumes:** after the pairing fix and the length-law re-measurement, with the §7 design questions, on paper first.

## 6. Multimappers: intergenic placements and the gDNA term

- **What is known:**
  - the shipped mean over gDNA placements dilutes gDNA's term when a poorer placement joins;
  - the sum alone over-states gDNA under capture;
  - every measurement was made with one-mate placements in the input.
- **Open (from the second review):**
  - units: calibration's absolute density cannot sit beside the EM's normalized weights without a derived conversion;
  - the sum is valid only over distinct valid origins with matching opportunity;
  - capture needs a conditional derivation;
  - the 12,288 false gDNA at `g00` under the sum are unexplained.

  The reviewer ranks first one gDNA component per locus, with a normalized spatial law over its complete support and
  a routing-conditioned opportunity.
- **Prototypes, not approved:** arm 1 (`<worktree>`) and the sum arm (`<worktree>`). Both
  fed one-mate intergenic placements to the EM. The break test of arm 1's second gate was never run.
- **Materials:** `INTERGENIC_GDNA_DESIGN.md` (with the reviewer's status note), `INTERGENIC_GDNA_REVIEW.md`, and the
  run directories `prototypes/2026-10-15_mm_intergenic/` and `prototypes/2026-10-15_mm_sum/`.
- **Resumes:** after the pairing fix, as part of the splice design's normalized likelihood (§5).

## 7. The aligned panel and validation

- **The local minimap2 panel:** five ladder conditions in `prototypes/2026-10-15_splice_a/aligned/` (13 GB, including
  the BAMs). It is the only local artifact-containing truth. Its baseline must be re-recorded on the fixed tree.
- **The truth ledgers** the reviewer asks for: all simulated molecules, molecules accepted as complete pairs, and
  mapping losses.
- **The `g00` isoform gap:** transcript error is 33 times the oracle's (gene error 65 times) with nearly right pools.
  It is dissected after the fix. The approximate census flagged 346 k of its 1.05 M "multimappers" as affected by
  one-mate hits; that is not a prediction that all become unique under the reviewed validator.
- **STAR (owner):** on the SLURM cluster after implementation and testing, using the drafted aligned ladder
  (`$RIGEL_RUNS/aligned_ladder/`). No local installation.

## 8. Housekeeping

- **`tests/test_docs_boundary.py::test_no_numbered_rule_labels` fails** on two older review files
  (`RNA_RELEASE_0_8_0_CAPTURE_REFERENCE_INDEPENDENT_REVIEW.md`, `RNA_RELEASE_0_8_0_MODEFREE_INDEPENDENT_REVIEW.md`).
  It no longer fails on mine: `SPLICE_ARTIFACTS_PLAN_REVIEW.md` now uses names.
- **Worktrees from this campaign:** `rigel-splice-a`, `rigel-mm-intergenic` and `rigel-mm-sum`,
  all at 41899ae4 plus their changes. Older ones are kept: `rigel-part1`, `rigel-psi1`, `rigel-sched`,
  `rigel-solverdiag`, `rigel-v071`, `rigel-release` and those under `/private/tmp`. None is deleted without the owner.
- **0.8.1 and later:** unchanged, in `ROADMAP.md`.

## 9. Gene-centric estimation from unannotated isoform paths

- **Owner scope, 2026-10-15:** a paired path might be supported by the combined isoforms of one gene while matching
  no annotated transcript. Using that evidence for gene abundance would require a separate estimator; it is not
  part of the current pairing repair or transcript quantifier.
- **Current handling:** follow [the unsupported-transcript contract](MATE_PAIRING_FIX.md#spliced-placements-without-a-supported-transcript).
  Gene-context diagnostics must not become additional ordinary gene-abundance counts.
- **Resumes:** only with an explicitly scoped future design for its likelihood, opportunity, normalization and
  relationship to transcript-derived gene abundance. No implementation is authorized by recording this direction.
