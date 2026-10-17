# Rigel 0.8.0: reassess the model before adding more machinery

**External scientific and implementation review, 2026-10-12. Sandbox, not a release decision.**

We need the smallest defensible change that makes Rigel a reliable transcript quantifier in the presence of genomic DNA and capture bias. We have identified regressions and corrected several mathematical and numerical mistakes, but the latest shared-field implementation has not earned its complexity through better results. Please reconsider the direction, including a narrower repair or selective rollback. Continuing the shared-field project is an option to justify, not the premise of this review.

This packet separates measurements, mathematical results, and hypotheses. The earlier campaign notes preserve failed approaches and superseded recommendations; their latest “next step” is not an instruction to you. In particular, the proposal to build a better joint-field optimizer is now subject to this review.

## 1. What the release must deliver

Rigel scans a name-sorted BAM, calibrates gDNA/RNA composition, and allocates fragments among compatible transcripts and a locus gDNA component. Biological populations are **gDNA, RNA+, and RNA−**. Synthetic unspliced “nascent” transcripts provide annotation coverage and a robustness test; they are not a fourth biological population.

Three strata are required: stranded capture-OFF, stranded capture-ON, and unstranded capture-OFF. Unstranded capture-ON remains deferred, but must be reported. Fragment-length differences in **both directions** are in scope. ON/OFF labels describe simulator conditions; production capture is a continuous local quantity.

Two independent primary measurements matter: calibration against origin-split oracle calibration, and delivered transcript counts against truth. Transcript error below is `100 × Σ|estimated count − true count| / Σ annotated true count`; gene and origin-pool errors accompany it. Calibration's six count arrays do not measure the capture-contracted lengths. Those need `ruler_vs_truth.py`, including within-gene length-ratio distortion, followed by allocation truth. A raw length/yield scale difference alone does not prove a likelihood error: numerator and denominator units must be examined together.

Fixed constraints: **no capture detector, no library-level reference density, no new empirical constant, no panel tuning**. Use local continuous evidence. Preserve conservative gDNA allocation under ambiguity. Do not silently relax these constraints; if an objective is impossible under them, specify the missing information and the exact conflict. Numerical iteration budgets and diagnostic tolerances are computational controls, not permission to fit biological constants.

For new comparisons: `scan.total_threads=1`, `OMP_NUM_THREADS=1`, `em.assignment_mode=fractional`; other stages retain their normal threading. Read every stratum and rung, including zero controls, separately. Derive before prototyping outside main; budget before execution; keep the machine idle for timing claims. The owner makes commits. Full requirements and instruments are in [SUCCESS](../SUCCESS.md), [DESIGN](../DESIGN.md), and [CLAUDE.md](../../CLAUDE.md).

## 2. Three different baselines

| State | Identity | What can be claimed |
|---|---|---|
| Released 0.7.1 | Tag `v0.7.1`, commit `ff26d8e59b7e0290b9a01a920714eb006173aec5` | An actual release, with a receipted October 5 comparison. |
| Current main | `624dfade4afce37d795890ea8ab537086ea20a93` | Runtime `src/` matches `0fdb9dce`; the later commit fixes the ruler truth instrument. Standing suite: 2,837 passed. |
| Shared-field prototype | Isolated source and adapter in the [implementation archive](../../.cache/rigel_runs/2026-10-11_shared_field_implementation/README.md) | Research code; no production install or inference change. No release candidate selected. |

Main still declares version 0.7.1; use the commit identities above to distinguish it from the released tag.

**0.7.1 was useful, but was not uniformly more accurate.** Its October 5 sweep used the release wheel, separate compatible indexes, fractional assignment, seed 0, and CLI `--threads 1`. Selected historical transcript errors:

| Condition | Tagged 0.7.1 | October 5 comparison tree |
|---|---:|---:|
| Ladder g00 stranded OFF | 3.70% | 1.93% |
| Ladder g00 unstranded OFF | 20.09% | 1.60% |
| RNA-short gap, stranded OFF | 6.31% | 42.30% |
| RNA-long gap, stranded ON | 7.38% | 6.17% |

The tag also assigned 41,473 and 401,708 false gDNA fragments on those two ladder g00 OFF rows. It undercounted synthetic RNA on gap OFF rows by roughly 127,000–276,000 fragments. These are reasons to audit a rollback, not reject it without testing. The 42.30% defect is historical; today's corresponding RNA-short stranded-OFF error is 4.9906%.

The historical `current_cli` column is **not today's main**. Its ladder summaries exclude g00 from contaminated-row aggregates; newer tables sometimes include it. Do not compare those averages directly. September 30 arms are also not the released tag. The [handoff archive](../../.cache/rigel_runs/2026-10-12_release_reassessment/README.md) includes the historical summary, manifest, and receipts so a reviewer can recover a matched comparison instead of relying on release nostalgia.

## 3. Current defects and what actually owns them

### A. Stranded transcript regression: the capture-weight path

The decisive **ladder** g00 stranded-OFF ablation holds the pipeline fixed and forces capture weights to one:

| Readout | Shipped | All weights one |
|---|---:|---:|
| Transcript error | 13.61% | 1.93% |
| Gene error | 0.83% | 0.14% |
| Mass on silent isoforms | 251,216 | 17,058 |
| End-to-end gDNA fragments, true zero | 1,126 | 1,112 |

The last row is an **EM pool**, not the primary calibration result. The latter reports 1,593.742 gDNA fragments in the review's cached score. The ablation attributes the transcript regression to the weight/length path; it does not authorize flat production weights.

The original explanation, “independent shallow Poisson noise,” was inadequate. On 70,176 ladder objects, 696 published log weights exceed 1, 352 exceed 5, and 93 exceed 10 despite flat capture truth. Their depths reach 3,769. Delivered composition and DNA-level factors cause many posterior mode switches. The unstranded reader is nearly flat because it is largely uninformative there: maximum log weight about 0.00005 versus 13.1 stranded. Its apparent robustness is not evidence of accurate field recovery. At g00 there is no gDNA density to estimate; constant capture truth does not supply a nonzero gDNA-density reference.

Two independent posterior readouts and three local/own-likelihood readers failed broader tests. Examples: `E[exp(x)]` changes ladder stranded-OFF transcript error 7.52→15.57%; `E[x]` changes stranded-ON 9.71→13.00%. The own-likelihood-only reader improves ladder averages but worsens RNA-short stranded-ON 14.2427→20.4344%. Both locally coupled readers also damage captured gap rows. Smoother weights are not sufficient evidence of a better quantifier. See [the expectation A/B tables](RNA_CAPTURE_WEIGHT_EXPECTATION_PLAN.md) and [implementation review](CAPTURE_IMPLEMENTATION_REVIEW.md).

### B. Deferred captured unstranded data: a missing discriminator

After removal of the intron-background constraint, captured exonic gDNA is often allocated to synthetic unspliced RNA. Ladder g50 unstranded-ON transcript/gene errors moved approximately 10.7/4.2→65.2/24.9%; inferred gDNA fraction 0.442→0.343 at true 0.50. Current deferred transcript errors are 34.1078% with longer RNA and 68.0552% with shorter RNA. An earlier longer-RNA arm achieved about 8.2% with gDNA fraction 0.499; this motivates causal examination of the removed constraint, not automatic restoration.

Two corrections are essential. First, **the current calibration has no direct fragment-length composition likelihood**: zeroing the inverse-length banks changes every calibration array by exactly zero. Length affects opportunities and downstream scoring, which is different. Second, an exon/intron total contrast cannot distinguish captured DNA `(700, 1)` from flat DNA plus exonic RNA `(1, 1, 699)` when strand and matched length laws provide no discriminator. Bulk observed length differences are not proof of local composition information. These findings correct the [original brief](CAPTURE_WEIGHTS_AND_DEFERRED_STRATUM_REVIEW_BRIEF.md).

Keep useful advances in view: the multimapper second pass brought VCaP's gDNA fraction to 0.252753, near its approximate read-name origin fraction 0.2518. That reference uses aligned primary first mates, not exactly the filtered pipeline denominator, and provides no transcript truth. Plasma libraries have no truth. Neither can certify isoform accuracy.

## 4. What the proposed shared field means

Current production estimates per-object capture weights, contracts each component's conserved-share length, and runs per-locus VB allocation with calibration-derived grouped priors. The proposal instead integrates uncertainty in one field shared by the locus's competing components.

For fragment `f`, component `k`, candidate path `h`, and deposit object `o`, let `d_fkho` be the conserved path share and `phi_fkh` its base weight. Define:

```text
B_fko = Σ_h phi_fkh d_fkho
K_fk(c) = Σ_o B_fko c_o
A_ko = Σ_(all admissible placements f) B_fko
L_k(c) = Σ_o A_ko c_o
p(f | a,c) = Σ_k a_k K_fk(c) / L_k(c),       Σ_k a_k = 1.
```

The gDNA kernel and opportunity both include the same existing fractional object-to-locus membership; retained partial shares are not renormalized. One fragment contributes one mixture likelihood factor. Several deposit objects are not independent observations. The field numerator cancels between components **only when their candidate field responses agree**. The submitted design's universal cancellation and one-object-per-fragment representation were wrong.

This normalization identity is established for the finite geometric model. Actual alignment errors, candidate pruning, and compatibility rules remain conditional scoring approximations. The prototype uses the same **observed** gDNA length law in kernel and opportunity, and likewise the observed RNA law. Algebraic consistency under this empirical base measure does not establish a physically correct capture-selection model.

With captured proportions `a` fixed, the intended objective integrates the product of fragment mixtures over the field prior, then adds the existing grouped MAP penalty `P_g log(a_g) + P_R log(Σ_R a_k)`. In log-field coordinates `x=log(c)`, sampled integer component/object allocations give the spatial conditional:

```text
log p(x | z,a,D) = log π_x(x) + Σ_j m_j log c_j(x)
                  − Σ_k N_spatial,k log L_k(c) + constant.
```

Junctions can be dependent functions of free coordinates. Missing/out-of-template positions have a separately marginalized channel and support restrictions. Integer statistics drive field transitions; conditional responsibilities are averaged for fractional output. “Collapsed” refers to integrating auxiliary variables out of the sampler, not to a closed-form solution of the field or abundance problem.

A separate rate model holds positional rates `lambda` fixed, with `W=Σ_k lambda_k L_k(c)` and `p(f|lambda,c)=Σ_k lambda_k K_fk(c)/W`. Its captured proportions depend on the field. At a known positive field these are reparameterizations; after field integration they are different models. The rate experiment moves the captured-origin prior inside the field integral. It also fails the first truth screen.

The reader's posterior and delivered factors are excluded from the new field likelihood. **Evidence reuse has not disappeared:** the landscape prior, length laws, calibration pseudocounts, and placement decisions are still estimated from the same library. The reviewer should decide whether this empirical-Bayes construction is justified and whether it duplicates information consequentially.

## 5. What the implementation has earned—and has not

### Mathematical corrections

| Claim tested | Result | Implication |
|---|---|---|
| Hybrid integrated responsibilities plus mean-field field update share an objective | On noiseless 50/50 toy: exact 50.101/49.899%, consistent mean-field 50.045/49.955%, hybrid 40.614/59.386%; objective decreases on 12/16 fixtures | Reject the hybrid, not every variational method. |
| Second-order log-length expansion handles broad uncertainty | Responsibility errors 0.0360 for lengths `{1,9}`, 0.1550 for `{1,30}`; negative probability for `{1,1000}` | This expansion cannot be the production shortcut. |
| Isoform structure identifies arbitrary local capture | Exact equal-likelihood ridge; 5/6 selected test pairs and 1/6 selected real-annotation pairs locally confounded | Identification is conditional. These selected, fully observed-path audits are not prevalence estimates. |
| Junction capture is bounded by adjacent gDNA-object averages | Simulator fixture: true junction average 44.5593 versus cap 27.7797 | The cap is not a universal physical identity. Free junctions are not thereby validated. |
| Correct conditionals imply useful posterior sampling | Gamma-augmented chains miss an exact count by 17.78–41.05 fragments; collapsed chains recover it within 0.0655 | A finite success does not establish large-locus mixing. |

The [follow-up review](SHARED_CAPTURE_FIELD_REVIEW.md) and its [executable references](../../.cache/rigel_runs/2026-10-11_shared_field_followup/README.md) contain derivations, enumeration, and deliberately broken alternatives. The native work adds correct path representation, continuous bounded slice updates, stable nonnegative length sums, and exact equal-kernel factorization. Path and fixed-field identities agree at roughly `10^-15`. These are correctness results for the implemented assumptions, not proof of biological adequacy.

### End-to-end screens

The following are **test-chromosome g00**, not the ladder's 13.61% row. Each field arm uses 1,000 MCEM updates ×16 sweeps, burn 256, 2,048 retained draws, seed 0; these are finite-budget screens, not certified optima.

| Stratum | Shipped | Shared field, dependent junction | Rate model | Free junction |
|---|---:|---:|---:|---:|
| Stranded OFF | 7.6610% | 10.9547% | 11.4295% | 20.8163% |
| Stranded ON | 8.2243% | 14.2810% | 15.4507% | 40.5560% |
| Unstranded OFF | 5.6250% | 14.9057% | 15.0008% | 25.5974% |
| Unstranded ON, deferred | 9.3964% | 14.6496% | 15.5310% | 43.0029% |

The dependent field also worsens unstranded-OFF gene error 1.0137→4.9256%. Free-junction stranded-ON synthetic RNA is 280,656 against 9,476 truth. All six calibration arrays remain unchanged across these four-condition comparisons. The new allocation has not fixed calibration.

The **g50 test controls** expose a major attribution problem:

| Stratum | Shipped VB | Shipped geometry, MAP | New geometry, flat field, MAP | New geometry, fixed reader field, MAP | Integrated field, shorter fit |
|---|---:|---:|---:|---:|---:|
| Stranded OFF | 9.2363% | 11.2056% | 6.2933% | 6.8669% | 7.2854% |
| Stranded ON | 9.8912% | 10.3054% | 32.1087% | 6.5782% | 7.5200% |
| Unstranded OFF | 10.5474% | 11.7934% | 6.9352% | 7.3864% | 9.3934% |
| Unstranded ON, deferred | 15.3084% | 17.4215% | 55.1909% | 11.9844% | 15.0050% |

“New geometry” includes candidate numerators, locus-membership handling, and the matched observed gDNA opportunity measure. This is **not a complete factorial isolation** of those changes. The integrated arm uses 256×4 sweeps and 1,024 retained draws. A fixed reader field beats that integrated arm on every transcript row. Therefore its apparent improvement over shipped cannot be credited to integration. A smaller numerator/normalizer repair remains a hypothesis worth separating.

### Unresolved optimization and cost

Deep g00 runs still have approximately 6–11-fragment fixed-point residuals. Increasing ordinary fits from 1,000 to 10,000 updates changes selected locus allocations by over 1,000 fragments. A frozen-field objective with safeguarded acceleration passes small references, but two starts at cached locus 123 differ by 1,965 fragments; their transcript truth errors are 144.7% and 89.3%. Posterior bridge comparisons lack adequate overlap. Small within-chain errors or residuals do not certify the solution.

The three complete g00 field screens cost about 114, 120, and 108 seconds across four conditions, with aggregate peak RSS around 4.5–4.6 GiB. An earlier four-worker gap attempt reached 55,587,088 KiB and was stopped after a numerical gate failed. These are cost receipts, not idle-machine speed comparisons. Exact factorization made a cached large locus cheaper; it did not validate whole-pipeline cost.

**There is no completed new-field ladder/both-gap accuracy matrix or VCaP score.** Earlier reader campaigns do not fill that gap. The full 2,837-test suite was not rerun against this research implementation. Public sampled assignment, read annotations, production configuration, and replacement of the reader's outputs outside EM loci are unfinished. Prototype annotation requests are explicitly refused.

## 6. Questions that should decide the direction

| Open question / hypothesis | Cheapest discriminating evidence |
|---|---|
| Are conserved bookkeeping shares a sufficient physical capture model? The simulator selects by best contiguous probe overlap, not necessarily a linear object field. | Enumerate placements for tiny exons/junctions under the actual simulator. Test representability of **both** candidate probabilities and integrated opportunities using one field. A normalization identity alone does not answer this. |
| Is observed length being conditioned on correctly, or is capture selection applied again to an already selected length law? | Derive the likelihood conditional on observed fragment length and compare it with the current joint kernel on an enumerable locus, for both gap directions. Do not invent an unseen precursor law to close the derivation. |
| Can zero/low-gDNA data identify the field separately from isoform expression and missing annotation? | Compare known-field, known-abundance, and joint inference on the same finite data; exhibit indistinguishable alternatives and prior sensitivity. More posterior computation cannot create absent information. |
| Is the smallest useful repair in candidate numerators/normalizers rather than field integration? | Split the fixed-reader g50 control into single changes; hold abundance update, length law, support, and reader fixed wherever possible. Price the g00 ladder defect as well as the test chromosome. |
| Are calibration messages or empirical priors generating the harmful modes through information reuse? | Trace one affected locus from raw counts through delivered factors to lengths and allocations. Remove one evidence path at a time; validate captured and sparse single-probe cases. Existing own-only failures rule out a blanket factor deletion. |
| Can length supply the missing local composition observation on gap data? | Derive and test a matched-footprint conditional likelihood from available observations. State what information is missing before requesting new banks. Equal-law controls must contain no invented discriminator. |
| Would a smaller, better identified junction relation suffice? | Compare exact probe-layout truth with a derived constrained relation. Neither the refuted hard cap nor an independent free coordinate is a validated answer. |

The competing release paths are: a narrowly attributed repair of current main; selective restoration of an older mechanism with its defects measured; a genuinely smaller conditional capture/composition model; or joint-field inference only if simpler models are demonstrably insufficient. “Defer unstranded capture-ON” is already a scope decision and cannot excuse the stranded regressions. More sampler engineering is justified only after showing the target model deserves to be optimized.

## 7. Code and evidence to inspect

Read the code in this order, following one fragment and one locus through it:

| Contract | Production source | Review concern |
|---|---|---|
| Scan and held placements | [accumulator.cpp](../../src/rigel/native/calibration/accumulator.cpp), [second_pass.py](../../src/rigel/second_pass.py) | Conservation, evidence ownership, choices made before EM. |
| Composition and reader | [calibrate.py](../../src/rigel/calibration/calibrate.py) (`_solve`, `capture_weights`), [honest_reader.h](../../src/rigel/native/honest_reader.h), [solve_kernel.cpp](../../src/rigel/native/solve_kernel.cpp) | Which observations/factors infer composition and capture; unsupported mode certainty. |
| Lengths and priors | [capture_eff_length.py](../../src/rigel/calibration/capture_eff_length.py), [priors.py](../../src/rigel/calibration/priors.py), [fl.py](../../src/rigel/calibration/fl.py) | Candidate measure versus opportunity measure; junction relation; observed law; grouped pseudocounts. |
| Allocation | [scoring.cpp](../../src/rigel/native/scoring.cpp), [pipeline.py](../../src/rigel/pipeline.py), [em_solver.cpp](../../src/rigel/native/em_solver.cpp) | Numerators, pruning, VB/MAP difference, initialization, output count/TPM contracts. |
| Capture truth | [sampler.py](../../src/rigel/sim/capture/sampler.py), [ruler_vs_truth.py](../../scripts/design/ruler_vs_truth.py) | Best-part probe physics and what the oracle actually measures. |

The new kernel is [shared_capture.cpp](../../.cache/rigel_runs/2026-10-11_shared_field_implementation/source/src/rigel/native/shared_capture.cpp) (`deposits`, `pack_locus`, `infer_compact`). Its [Python adapter](../../.cache/rigel_runs/2026-10-11_shared_field_implementation/site/shared_pipeline.py) replaces quantification through monkeypatching; the archive's `source/` and `tracked_changes.patch` retain scorer/Python/build changes. **The patch alone omits the new, originally untracked C++ file.**

Implementation debt is concrete: a roughly 40 KB dense native file combines packing, several statistical models, sampling, and optimization; environment flags select research branches; the adapter uses module globals and implicit hooks; old scripts contain machine-specific paths. These are acceptable experiment scaffolds, not an elegant production design. The proposed production boundary must state what disappears, not just where these branches will move. Exact reference cases and provenance should survive simplification; failed mechanisms need not survive in runtime code.

The source/evidence handoff contains the current and tagged source snapshots, the prototype, mathematical references, selected JSONL arms and receipts, historical summaries, and the documents cited here. See its [inventory and limitations](../../.cache/rigel_runs/2026-10-12_release_reassessment/README.md). Full BAMs, indexes, cached loci, and native binaries are excluded. Existing binary artifacts are platform-specific; some earlier binary/adapter combinations lack complete standalone snapshots. Do not claim a replay from unavailable inputs. Experimental receipts identify the model; a legacy JSONL `arm=base` value does not.

## 8. Requested review and bounded next step

Please return:

1. A verdict on whether our estimand, observation model, and evidence ownership are sound, with equations or a counterexample for each fundamental objection.
2. Findings ranked as model error, non-identifiability, approximation, implementation defect, or missing measurement. Cite the relevant code and state what would refute each finding.
3. A **keep / delete / replace** decision for each major mechanism, including whether to abandon the shared-field route. Explain what accuracy each retained complexity buys.
4. One preferred, minimal candidate with derivation, smallest code boundary, falsifiers, and an explicitly budgeted validation sequence; explain why the simpler alternatives fail.

Before new large runs, reuse the archived tagged comparison and make the exact release/main/prototype estimands comparable. Then settle the finite observation-model question and isolate one change on cached loci/test data. Historical test receipts suggest seconds for a fixed-field control and roughly two minutes for a four-row field screen; these are planning inputs, not promised convergence budgets. A surviving mechanism must face every in-scope ladder condition, both gap panels, controls, independent calibration/ruler scores, and VCaP before selection. Budget those from a representative measured pilot; no new genome-scale field timing estimate is established.

**Our present recommendation is to withhold promotion and reassess the statistical contract before extending the optimizer.** We refuse to choose a start by truth, tune to panels, hide g00 in aggregates, revive flat weights from an attribution experiment, or call passing numerical gates a release. A concise model with explicit limits and better measured behavior is the desired outcome; a larger inference framework is not itself progress.
