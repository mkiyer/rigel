# Rigel 0.8.0: design and validation plan

Revised 2026-10-15 after the owner's review. This supersedes the implementation-first order and the proposed
single-mate survival likelihood. Sandbox; no production source change or commit is authorized by this plan itself.

The owner has decided: the length defect blocks release; each ordinary hypothesis must place both mates; abundance
uses mappable genomic opportunity, including observable multimapping sequence; establish the full splice-artifact
design before implementing it; run STAR on SLURM after local implementation/testing, without a Mac installation.

## Evidence and immediate conclusion

`INTERGENIC_GDNA_REVIEW.md` records the independent review. The current scanner emits incomplete primary, HI-grouped
and secondary hypotheses. A native six-case reproducer confirms that adding a lone secondary to a valid 300-bp pair
can put the fragment into calibration at 100 bp. Complete paired remote alternatives preserve the 300-bp length.

The pinned full-library diagnostic reproduces the collapse. The length contrast's internal density is 0.000798,
distinct from final calibration's 0.0521. Its signed 100-bp mass is 5.084 before nonnegative projection; clipping and
renormalization produce approximately 76.3%. Rejecting malformed hypotheses is necessary, but has not yet been
measured to repair the entire shape. Mappability opportunity and contrast assumptions are part of the investigation.

The original intergenic sum/mean experiments establish mean dilution, but are confounded by invalid candidates and
the length-law failure. Neither prototype is approved to merge. The reported capture-OFF false-DNA increase, gDNA-zero
transcript error and the unrun second break test remain open.

## Phase 0 — finish the full design (M; current work)

0.1 Establish the complete-pair candidate contract in `SPLICE_ARTIFACTS_IMPLEMENTATION.md`: per-mate chains, whole
molecule geometry, local/remote paths, true short overlaps, supplementary alignments, rescue and rejection reasons.

0.2 Derive one normalized likelihood and opportunity over the same accepted paired observation population. Resolve
capture/length conditioning, mappability, distinct origins versus duplicate records, locus routing, external density
units, prior mass and regional booking. Preserve the three populations and hard calibration draw.

0.3 Resolve the calibration hold/training design. Remove the base-only hold rule. No provisional deposit may become
its own support through counts or fitted parameters. Audit the internal length-density estimator and common-shape
assumption; no bin removal or ad hoc smoothing workaround.

0.4 Define the treatment of genuine unannotated RNA, reference-copy search gaps and sample structural alleles. Admit
identifiability limits explicitly. Do not infer an allele because no reference placement was found or turn expected
DNA capacity into a hard count ceiling.

**Gate:** all candidate, likelihood, denominator, training and conservation definitions agree on finite enumerated
references. Remaining assumptions and unsupported classes are explicit. The full design, not merely the first bug
patch, must be reviewable before source implementation. Read-only diagnostics and falsification calculations continue.

## Phase 1 — establish the instruments and failing fixtures (M)

1.1 Preserve the supplied baseline and evidence. Build the paired-fragment fixtures listed in the design, including
all existing orphan branches and the annotation-first/greedy pairing counterexamples. Verify failures on the current
tree. Distinguish geometry, probability normalization and count-ledger tests.

1.2 Define three truth ledgers: all simulated molecules, molecules accepted as complete pairs, and mapping/acceptance
losses. Evaluate accepted fragments against their original source identities. All-N molecules remain visible as losses;
do not let them enter observable opportunity or use scoring exclusion as a substitute for a denominator correction.

1.3 Create a shared per-fragment diagnostic census: candidate completeness, observed/proposed placement provenance,
true origin, length before/after drain, pool, contrast coefficients, capture and final responsibilities. Check the
length contrast's internal density separately from final calibration density.

**Gate:** hand-enumerated counts/opportunities and the simulation ledger reconcile; every planned mechanism has a
specific falsification test. Local minimap2/oracle input suffices for these instrument gates; STAR remains later.

## Phase 2 — prototype complete-pair acceptance and observable length calibration (M–L)

2.1 In an isolated worktree, repair candidate assembly and validation. No ordinary orphan deposits or abundance
hypotheses. Preserve valid overlapping and rescued pairs. Repair greedy pairing/chain handling with separately
attributable tests. Expire scan caches whose candidate/deposit semantics changed.

2.2 A/B this mechanism on identical read identities, then separately introduce the approved opportunity correction.
Re-estimate component length laws with the correct selected population. Decompose remaining contrast errors; if its
common-shape/purity model fails, derive and prototype the replacement separately.

2.3 Re-measure both origin errors, calibration, prior and abundance tables. The gDNA-zero aligned condition is a
required control; dissect its residual transcript/gene error by complete candidate class after pairing is repaired.

**Gate:** invalid hypotheses are impossible, length/selection mathematics agrees with truth under controlled
generation, and no unexplained spike follows read length when true inserts are held fixed. The oracle mean alone is
not the target: aligned acceptance can legitimately change the observed length distribution. Break every repair and
verify its gate fires, including the previously outstanding prototype gate.

## Phase 3 — prototype the complete splice likelihood (L)

3.1 Implement the agreed whole-pair base/error model and local alternatives, with mixed genuine/artifact RNA paths.

3.2 Introduce validated reported remote pairs and the copy-proposal index. Compare origin summation and its matching
opportunity/prior accounting; isolate the capture factor in its derived conditional frame. Any fixed external rate
requires the normalized model settled in Phase 0, not the original arm 2 shortcut.

3.3 Implement the agreed calibration/training contract. Implement structural-allele treatment only after its geometry,
prior and identifiability gates pass. Test withheld annotation and unstranded ambiguities explicitly.

**Gate:** one-mechanism A/Bs pass the finite-reference and two-error tests; no hidden truncation or unsupported
sequence is certified as a complete alternative. No capture-OFF or gDNA-zero cost is waived because capture-ON wins.
Both false DNA and false RNA and their transcript/gene consequences are explained per stratum.

## Phase 4 — local integration and regression validation (M)

Run the oracle ladder, both fragment-length-gap panels, the existing minimap2 panel and the VCaP source mix. A
mechanism that cannot apply must be a numeric no-op under unchanged support/opportunity; do not demand identity where
the denominator or accepted population deliberately changes. Compare blacklisting with the replacement on both
artifact leakage and genuine RNA loss.

Pin `scan.total_threads=1`, `OMP_NUM_THREADS=1` and fractional EM assignment for A/Bs. Score each in-scope stratum
separately: unstranded capture-OFF, stranded capture-OFF and stranded capture-ON. Report deferred unstranded capture-ON
separately. Re-record baseline from the current source; do not mix unrelated changes into a comparison.

**Gate:** meaningful fixture and required regression checks pass, break tests fire, all affected instruments run, and
no unexplained in-scope regression remains. Existing unrelated documentation-gate failures must be reported honestly;
they are not evidence that a new mechanism passed or failed. No blanket golden update.

## Phase 5 — STAR and independent aligned validation on SLURM (L)

Run the same simulated reads through STAR with recorded version and settings, then the complete aligned panel and
independent artifact-containing libraries. Include alternate reporting/settings that stress candidate completeness,
clipping and annotation dependence. Compare with minimap2 on accepted molecule truth and exposure as well as total
mapping losses; a lower secondary-record count is not a correctness proof.

**Gate:** complete-pair/length/opportunity invariants survive both aligners; both source errors and transcript/gene
accuracy meet the standing release criteria in each required stratum. A clean STAR arm does not waive a confirmed
scanner defect or unvalidated source-specific mechanism. Retire the blacklist only if the paired comparison supports it.

## Phase 6 — source integration and release (S–M)

Only after the derived, designed, prototyped and measured mechanisms pass: converge the chosen implementation into
production source, rerun required checks and instruments, reconcile permanent documentation, and prepare release
artifacts. Preserve decision provenance and update any index/cache manifest whose semantics changed. The owner drives
commits and publication. Do not silently promote provisional design claims into permanent rulings.

**Gate:** no open release blocker, complete validation record, reproducible build/release checks and the owner's
release actions. This plan intentionally supplies no invented numerical tolerance or promise of perfect inference
where two complete physical explanations are observationally identical.
