# Independent review: spliced-fragment artifacts

2026-10-14 · Design review of `SPLICE_ARTIFACTS_DESIGN.md` · tree `7805ede7`.

Scope: packet, bundle and current source; all 17 compared bundle files matched. Review experiments used the native EM and local derived junction table; no realignment or end-to-end A/B was performed. **RNA-source** denotes provenance, not verified genuine junctions. The three populations, solver and hard-draw ruling are preserved.

## 1. Verdict

**Build the competing-placement architecture after revision; do not build §5 as specified.**

- **Blocker finding 1:** `r` conflates an emission probability, a sequence likelihood ratio and an observed mixture share; their product has no stated generative interpretation.
- **Blocker finding 2:** clipping changes the hypothesis behind the sequence penalty, and independent junction edits do not establish a physically possible whole-fragment gDNA origin.
- **Blocker finding 3:** pruning on `r` alone is not the existing likelihood rule; a concrete counterexample changes candidate eligibility.

Singletons remain prior-dependent; long-anchor novel artifacts need another mechanism. Decide the blacklist’s fate through a matched two-error A/B. Shipping requires independent artifact-containing validation.

## 2. Findings

**1 — blocker · §§4.1, 5.1, 6: define what the weight conditions on.** *Derived.* A class share `a_c = P(gDNA | reported splice, class c)` already reflects gDNA/RNA prevalence. It is neither `P(reported splice | gDNA)` nor a Bayes factor. Multiplying it by a sequence ratio and then by locus abundance generally reuses prevalence; class and `h` also need not be independent. A strand-conditioned fragment posterior would additionally reuse that fragment’s strand. Fitting nuisance parameters on the same data is not, by itself, double-counting; a coherent joint likelihood or empirical-Bayes fit can do that legitimately. **Test:** generate mixtures with fixed component emissions and varying gDNA abundance. Fit each proposed weight independently; posterior odds must follow Bayes’ rule and recover mixture proportions without a second prevalence effect. Include an unstranded case: class information must disappear when its only evidence is strand.

**2 — blocker · §§5.1–5.2, 9: local ambiguity is not a whole-fragment origin.** *Derived.* The `h` penalty assumes every moved base came from contiguous reference. Deleting those bases instead makes it a different emission model. A terminal `h=0` junction does not erase another well-supported junction or a contradictory mate. Alternative origins are summed, not multiplied as if each mechanism must occur. Multiplication by a class share cannot rescue a very small local-contiguity likelihood. *Measured:* §4.3’s harmful subset contains **2,116 unannotated units of 2,944**; most have long anchors and `h≥2`. **Test:** construct paired fragments with terminal/interior artifacts, genuine-plus-artifact junctions and remote unspliced origins. Compare candidates with exhaustive whole-read placement enumeration. Any gDNA candidate requiring unexplained bases or a retained biological splice fails; so does claiming coverage of remote artifacts using only local edits.

**3 — blocker · §5.3: the early pruning rule changes the admitted set.** *Derived, numerically checked.* With the packet’s illustrative `ε=10⁻³`, `κ=1.26×10⁻⁴`, neutral class evidence and equal length terms, an `h=2` sequence ratio is **1.113×10⁻⁷**. On the wrong strand it becomes **4.418×10⁻⁴** after multiplying by `0.5/κ`, exceeding the existing **10⁻⁴** pruning bound. The early screen drops it. Fragment-length and overhang differences can likewise reverse the screen. **Test:** compare full scoring with screened scoring over the supported likelihood terms; the admitted sets must agree. A fast rejection needs a proven upper bound on the remaining terms. The packet’s “about 2%” is a workload hypothesis, not a bound.

**4 — fix before building · §§3, 5.3–5.4: the stage contract is incomplete.** *Derived from code; existing checks passed.* Contrary to §3, calibration defines spliced as `!sj_ids.empty()`: wholly unannotated observed introns use unspliced channels while their gaps remain cut. The EM excludes gDNA for them. Moreover, the second pass skips its strand term whenever an observed motif exists; a spliced flag alone will not restore competition. See [accumulator.cpp](../../src/rigel/native/calibration/accumulator.cpp), `deposit_survivor`, and [second_pass.py](../../src/rigel/second_pass.py), `score_held_fragments`. **Test:** replay annotated, unannotated, motif-present/absent, contradictory-mate and entirely-held junction fixtures through both stages. Require the declared geometry, likelihood factors, training eligibility and conservation identities; an absent pass-one junction must not become decisive evidence against RNA. Two existing unannotated specification tests and the native parity battery passed during review, confirming the baseline correction.

**5 — fix before building · §10.8: the calibrated prior is not a classification guarantee.** *Derived and measured on a synthetic native-solver input.* For one eligible fragment, one RNA candidate and one gDNA candidate, all prior mass favoring gDNA still leaves **0.94737 RNA under MAP** and **0.94694 under VBEM’s reported assignment** when the effective likelihood ratio is `B=0.05` (§3 below). A right-strand sequence tie can reach that ratio through effective lengths alone. **Test:** cross calibrated share, effective-length ratio, strand, fragment-length ratio and independent RNA support; compare MAP/VBEM and the existing coverage/uniform starts. Reject a claimed singleton remedy if it merely changes certification while leaving the false transcript abundance. Exact observational ties cannot require perfect source-label recovery.

**6 — fix before building · §§4.3, 5.1, 6: alignment validity and molecular origin are different latent facts.** *Derived.* RNA can be misaligned, and real novel RNA may lack a compatible annotated transcript. Existing overlap candidates retaining the novel intron are not necessarily an adequate RNA alternative. Thus a wrong-strand excess is not automatically gDNA: with a reversed genuine fraction `b`, `q=κ+(1/2−κ)a+(1−2κ)b`. One pooled rate cannot identify both `a` and `b`; at `κ=1/2` it identifies neither. The packet’s 28% inference is model-dependent, not origin truth. **Test:** remove known transcripts/junctions from the quantifier’s annotation while retaining their true RNA reads, and include reversed and misaligned RNA controls. Failure means genuine omitted isoforms flow into gDNA or unrelated annotated transcripts. Any unresolved transcript identity must remain distinguishable from molecular-origin uncertainty within gDNA/RNA+/RNA−.

**7 — fix before building · §5.1: NM does not determine the proposed sequence ratio.** *Derived.* The approximation assumes the read matches the spliced reference at every discriminating base. The same `h` and total NM can instead describe a read matching the contiguous reference, reversing the evidence. **Test:** synthesize equal-NM reads with the mismatch moved onto/off the discriminating positions, plus insertions/deletions and clips. Compare against explicit base likelihoods. A single fixed `(ε/3)^h` must fail this test; exact equality at `h=0` should survive it where both complete placements are valid.

**8 — fix before building · §4 tables: repair the measurement contract before ranking methods.** *Measured/code-derived.* The “library-wide” screen uses the NPMLE’s **junction-weighted mean, 0.00355**, as a fragment prior. The local table has **7,120 wrong-strand units / 14,932,781**; the simple two-class pooled estimate at the packet’s κ is **0.0007018**, a different estimand, not a corrected truth estimate. The NPMLE also reached its 500-step cap. Separately, **32,019 of 82,001 RNA-source `h=0` units** came from multi-junction records, so many local ties are not whole-fragment ties ([local table audit](\<scratch>/rigel_splice_review_screen_audit.json)). **Test:** refit the library baseline on the same fragment-weighted likelihood, check numerical convergence, then score complete fragments. Separate missing XS, multimappers and multi-junctions. A changed ranking refutes the original screening comparison; unchanged ranking strengthens it. DNA-majority at one junction must not be treated as proof that its gene is unexpressed.

## 3. Improvements, ranked by expected benefit against risk and effort

**1. Specify component likelihoods and candidate support first — highest benefit, moderate effort.** For fragment observations `x`, retain the settled abundance update and write the competing weights as

\[
w_D=(\theta_D/\ell_D)\,f_D(L_D)\,s_D(o)\,e_D(x),\qquad
w_t=(\theta_t/\ell_t)\,f_R(L_t)\,s_t(o)\,e_t(x).
\]

Here `e` holds sequence/alignment evidence not already in the other factors. An additive `log r` is valid when it implements the appropriate **emission likelihood ratio**, with the RNA reference likelihood accounted for; such a ratio can exceed one. Distinct physical origins contribute additive terms normalized by their genomic opportunities; duplicate alignment records add no probability. Keep three populations; local stitching and remote mapping are observation mechanisms, not extra populations. A contiguous geometry can also be an admissible RNA+/RNA− reading: misalignment alone must not exclude it. Count any existing NM penalty once: retain a common baseline and a correctly derived relative correction, or replace the corresponding base term consistently.

With matching conditioning, `odds(a_c)=odds(π)×P(c|D)/P(c|R)`. Opportunity, ascertainment and correlation with sequence features must match; dividing by an arbitrary global gDNA share is invalid. Use strand once; junction pooling requires a coherent parameter fit. **A/B:** local emission, then conditional class evidence separately; apply findings 1/3 and read both errors plus transcript tables.

**2. Explain all sequenced bases and both mates — high benefit, moderate/high effort.** The largest block is a useful seed, not a uniquely justified origin. Enumerate plausible contiguous placements anchored on either side and resolve both mates together; deduplicate identical physical hypotheses. For `[s−a,s)|N|[e,e+b)`, the right-anchored contiguous interval is `[e−a,e+b)` and the left-anchored interval `[s−a,s+b)`. Score the relocated bases at either `h=0` or `h>0`; reference disagreement should change likelihood, not delete molecule length. Keep original soft-clipped sequence distinct from proposed clipping. Unknown clipped endpoints require marginalization or an explicit missing-sequence model, not a shortened fragment asserted as truth.

For multiple junctions, compare the complete read sequence under complete candidate paths. Multiply base emissions once per sequenced observation, not overlapping junction penalties; shared biological variation across overlapping mates needs consistent treatment. A gDNA hypothesis must explain every apparent splice without biological splicing. Remote origin is a hypothesis, not the only explanation: a genomic deletion allele can also resemble a splice and belongs in the controls. Preserve partial-artifact RNA paths when another junction is real. Carry coordinates, observed/implied introns, strand, length and emission evidence together into calibration and EM. **A/B:** change geometry alone, then candidate enumeration; finding 2’s exhaustive fixtures precede paired-fragment measurements.

**3. Use a defensible sequence likelihood — high benefit at `h>0`, moderate effort.** Under independent symmetric substitutions, when all `h` discriminating bases match the RNA reference,

\[
B_{seq}=\left[\frac{\varepsilon}{3(1-\varepsilon)}\right]^h.
\]

A base matching the contiguous reference contributes the reciprocal; a third base contributes one under that model. At `h=0`, the compared sequences are identical, so base evidence cancels even if the read has errors. This does not make mate geometry, length or origin priors equal.

Estimate **aligned substitution divergence** as `Σ(NM−inserted bases−deleted bases)/Σ(M,=,X bases)` over the declared unspliced sample. Variation, misalignment and selection prevent interpreting this as pure sequencing error; missing NM is missing information. MD localizes mismatches but supplies reference bases, not observed alternatives, and omits skipped introns. Exact evaluation needs SEQ/alternative reference, with QUAL when available. [SAM optional-field specification](https://samtools.github.io/hts-specs/SAMtags.pdf).

Retain mismatch positions, not only prefix counts, and provide reference access for novel junctions. **A/B:** NM approximation versus direct base likelihoods, then QUAL separately; inspect disagreements, allele-bearing reads and scan cost.

The aligner’s output is a function of sequence and settings, not independent molecular evidence. Do not exponentiate its score into an uncalibrated probability. STAR 2.7.11b adds `sjdbScore` to annotated stitches; its negative genomic-length coefficient penalizes larger genomic span. Test changed settings and missing competing alignments. [STAR stitching source](https://raw.githubusercontent.com/alexdobin/STAR/2.7.11b/source/stitchAlignToTranscript.cpp), [parameter definitions](https://raw.githubusercontent.com/alexdobin/STAR/2.7.11b/source/parametersDefault).

**4. Finish the single-drain contract before training changes — high benefit, moderate risk.** Keep the hard seeded draw and downstream soft allocation separate. Share candidate geometry and observation evidence, not the realized assignment: conditional on frozen weights, a draw has expectation `p` and variance `p(1−p)`, whereas conditioning EM on that draw discards alternatives. Equal emissions do not imply equal posteriors because the stages use different priors and currently different genomic length laws.

Define held-junction support from independent pass-one evidence, with absence caused by holding distinguished from observed zero support; do not add a floor. Update the observed-motif strand bypass as well as spliced flags, deposit paths, denominators and training eligibility. Fit from pass one, score once, draw/deposit once, then fit final models; do not feed final EM decisions upstream. Training weights must mean probability of RNA under a valid RNA geometry, not merely `1−P(gDNA)`. Keep integer deposits and weighted parameter-fit statistics distinct, preserving opportunity corrections and testing selection bias from the held set. **A/B:** scorer/deposit contract, then training reweighting separately; fixed seeds for comparisons and repeated seeds to measure the draw’s downstream variance.

**5. Quantify singleton limits rather than promising recovery — low effort, essential risk control.** Let `B=(L_D/ℓ_D)/(L_R/ℓ_R)` include all fragment likelihood factors and effective lengths. For one eligible unit, calibration share `q` gives the actual pseudocounts `P_D=q`, `P_R=1−q`. The MAP fixed point for gDNA responsibility `z` obeys

\[
z=\frac{B(z+q)}{B(z+q)+2-z-q}.
\]

Thus `B=1` gives `z=q`; for `q=1`, the interior solution is `z=B/(1−B)` when `B<1/2`. A pure-gDNA **pseudocount** is not a hard constraint. VBEM instead uses `z=logistic(log B+ψ(z+q)−ψ(2−z−q))` internally. The shipped final assignment uses normalized alpha with MAP-style weights, so its reported responsibility differs from that internal `z`.

| Synthetic input | MAP reported gDNA | VBEM reported gDNA |
|---|---:|---:|
| `q=0.5, B=1` | 0.50000 | 0.50000 |
| `q=0.9, B=1` | 0.90000 | 0.94999 |
| `q=1, B=0.05` | 0.05263 | 0.05306 |

These are review measurements, agreeing from coverage/uniform starts; they are counterexamples, not tuning values. Independent density, mate/sequence evidence and other RNA fragments can change the odds. Without distinguishing evidence, neither an extra constant nor another EM iteration creates identifiability. **A/B:** finding 5’s matrix before and after each mechanism, including the prior induced by the actual calibration draw. Reproduction: [review-only audit script](\<scratch>/rigel_splice_review_audit.py), [results](\<scratch>/rigel_splice_review_audit.json).

## 4. Unannotated junctions

**Ranking: 3 > 5 > 4 > 2 > 1.** These are complementary operations, not five interchangeable complete models.

| Rank / owner strategy | Recommendation |
|---|---|
| 1 / **3: both readings** | Best architecture when both are physical, scored hypotheses and the RNA side can represent genuine novel splicing. |
| 2 / **5: known-end classes** | Useful first grouping for a correctly normalized observation model, combined with 3. Match reference, strand and low/high endpoint role. The screen’s coordinate-only membership is insufficient. |
| 3 / **4: penalize RNA** | Equivalent to adding a relative gDNA term only with a valid likelihood ratio and a competing candidate. A penalty common to all RNA candidates alone cancels on normalization. |
| 4 / **2: longest-block unspliced only** | Retains some abundance evidence but misrepresents genuine novel geometry. It need not lose every genuine read: unspliced candidates can still be RNA. That does not justify the placement. |
| 5 / **1: discard all** | Removes 2,253 known DNA artifacts at the cost of 114,519 RNA-source units; genuine loss within the latter is unknown. Throws away recoverable evidence. |

These classes describe annotation, not truth; record annotation completeness and the aligner’s annotation. Near-zero estimated shares do not prove impossibility: report uncertainty without arbitrary smoothing. Consider motif and continuously modeled intron length/distance as conditional emission features, never VCaP-tuned bins. Site distance tests shifted-junction mechanisms. Depth and orientation support need sampling opportunity: R1/R2 from one molecule are not independent confirmations, nor must stranded RNA have opposite molecular orientations. Test each feature separately using annotation withholding, depth thinning and independent libraries.

**Unstranded:** do not assign the VCaP shares, zero, or an invented fallback share. Omit an unidentified class factor—likelihood ratio one means no additional class evidence, not 50% artifacts—and let physically valid sequence/length/density evidence compete. Report the resulting class share as model-inferred, not strand-measured. A stronger class term requires identifiable non-strand information, such as exposure-normalized matched-aligner simulations or independent origins. If all relevant emissions are identical, retain uncertainty; the class share is not estimable from that library alone. Missing remote origins remain an explicit coverage gap.

## 5. The blacklist and `alignable`

**Target: retire the blanket cut; retain simulation evidence with a derived likelihood.** Flags/counts/maxima lack exposure denominators and reporting probabilities. Retain origins, read/fragment lengths, settings and primary/secondary reporting behavior, separating local ambiguity from remote contributions without counting an origin twice.

The **1,155 annotated `h≥1` catches** are a substantial residual—**16.4% of 7,032 annotated single-junction DNA records** in that old run. They are not dispensable because another class improves. Conversely, **598,010 RNA-source rejections versus about 20,000 extrapolated from DNA** argues strongly against retaining the blanket cut. The roughly 30-fold figure is neither an exact genuine-error count nor the cost/benefit ratio of the 1,155-record subclass: different read lengths, capture and record definitions prevent that inference.

Do not keep a long-anchor-only hard cut without measuring its own genuine-read loss. Rebuild comparable current-tree arms: no blacklist, current production blacklist, corrected candidate model, and candidate model with simulation evidence; test any residual hard-cut arm separately. Read both errors and transcript tables in precisely that residual class. `alignable` can become optional without deleting a compact validated splice/origin product; its 26 GB packaging is a separate deployment choice. If removal exposes unacceptable residual false expression, that stops release pending evidence, rather than justifying an unmeasured cut.

## 6. Validation and shipping gates

Freeze mechanisms before testing. Re-record matched baselines; pin `scan.total_threads=1`, `OMP_NUM_THREADS=1`, and `em.assignment_mode=fractional`. Rebuild scan caches after candidate/deposit changes. Compare each addition against its parent and original baseline. Separate the three release strata; report deferred unstranded/capture-ON separately.

For each library, report both **DNA-source RNA posterior mass** and **confirmed genuine-RNA mass lost**, with RNA-source loss as a separate imperfect proxy. Also report hard certification, zero-truth transcript expression, transcript/gene errors, calibration/prior errors, κ, length-law changes and runtime/memory. Break down known ends, `h`, strand, singleton support, junction multiplicity and mapping multiplicity. VCaP source labels cannot supply absolute per-transcript abundance truth; read its table shifts and DNA-origin leakage alongside the truth-scored synthetic tables.

| Minimum evidence | What stops shipping |
|---|---|
| **Mechanism fixtures:** exact ties, discriminating bases/indels, both mates, interior/multiple junctions, allele differences, remote origins, missing tags, all-held support, omitted true annotations and finding 5’s singleton matrix. | Invalid molecule geometry; duplicate/lost fragments; screened/full candidate disagreement; unexplained stage differences; unacknowledged prior-driven false expression. |
| **Artifact-free oracle ladder and both fragment-length-gap panels**, plus existing tests and calibration instruments. | Any existing gate fails or a reproducible in-scope transcript/calibration regression lacks an accepted explanation/tradeoff. These panels test harm, not artifact sensitivity. |
| **VCaP matched A/B**, including both source libraries separately and the production-blacklist comparison on the current tree. | Artifact RNA mass does not fall; apparent improvement only relabels fragments; confirmed genuine loss increases; or gains conceal deterioration in the harmful/long-anchor classes. RNA-source provenance alone cannot certify specificity. |
| **Independent artifact-containing evidence before release:** an aligned truth panel covering the three release strata and relevant mechanisms, or independent truth mixes covering missing strata plus mechanistic aligned controls. | No independent evidence; failed unstranded/capture validation; benefit disappears on a new library; or genuine-novel/remote-origin controls fail. An aligned panel with no artifacts does not satisfy this gate. |
| **Aligner-setting controls on known-origin reads:** vary annotation bonus, stitching and reporting settings one at a time; exercise another aligner before claiming general support. | Learned rates silently transfer across incompatible settings, discarded alternatives become systematic false expression, or absent optional tags change the interpretation without a defined fallback. |

Repeat draws to separate bias from sampling variance; judge aggregate mass/abundance, not perfect labels on ties. Use existing acceptance criteria. No genuine-loss budget has been authorized: a new tradeoff needs an explicit decision. Local evidence supports prototyping, not release.

## 7. Answers to §10

**1. Rule/statistics.** Competing gDNA hypotheses are sound with valid emissions and strand counted once. A class posterior share is not an emission factor. Coherent parameter estimation does not itself double-count evidence.

**2. Tie, ε, MD, aligner score.** The corrected ratio requires matching-base and symmetric-error assumptions. NM estimates aligned divergence; MD localizes mismatches; SEQ/alternative reference resolve their direction. STAR scores govern reporting, not molecular probabilities.

**3. Novel classes.** Start with strand-aware known-end classes and annotation-withholding controls. Add features separately. On unstranded data, leave unidentified class evidence neutral or obtain independent emissions; no strand-derived fallback exists.

**4. Geometry.** Largest-block anchoring proposes an origin. Score both directions, both mates and complete multi-junction paths; retain or model clipped bases. Separate local and remote mechanisms.

**5. Calibration.** Keep separate hard and soft decisions with shared candidates/emissions. Repair the motif bypass and withheld-support prior; measure seed variance. Avoid upstream feedback from EM.

**6. Blacklist.** Retire the blanket cut after matched validation. Preserve normalized simulation evidence. The long-anchor cut needs its own two-sided error measurement before becoming a fallback.

**7. Shipping.** Require mechanism fixtures, oracle-panel harm checks, VCaP’s paired errors and independent artifact-containing validation across supported strata. Cluster unavailability does not waive the last requirement.

**8. Isolated read.** MAP follows the calibrated share at a true likelihood tie; effective lengths can defeat all-gDNA pseudocounts. VBEM differs. Independent evidence can resolve ambiguity; no rule identifies the true source of indistinguishable observations.
