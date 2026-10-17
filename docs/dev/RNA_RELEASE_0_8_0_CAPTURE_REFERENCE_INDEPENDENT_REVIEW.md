# Independent review: Rigel 0.8.0 capture reference

2026-10-13. Reviewed the [complete packet](RNA_RELEASE_0_8_0_CAPTURE_REFERENCE_REVIEW.md), all four figures, both readers, the one-sided diff, analysis scripts and archived results. Verified all **229 manifest hashes**, archive contents and 48 arm receipts; independently reproduced the **265-condition / 118-stratum** comparisons and integrated/release JSONL identity excluding timings. Reran the two small stress scenarios and numerical counterexamples. This was not a rerun of the full matrix or real BAMs. [Audit and reproduction records](../../.cache/capture_reference_independent_review/README.md).

## 1. Verdict

Choose the **integrated two-sided, geometric reader** as the 0.8.0 candidate; retain its bulk unit for continuity.
Do **not** restore the captured-mode gate or ship the present one-sided null.
Treat largest-basin identity as a working assumption, not an identified uncaptured population.
Hold release sign-off on the unresolved **in-scope stress regression** and the missing causal sparse-bulk experiment.
Run the spike-location oracle below before choosing a replacement baseline estimator.
If that arm explains the sparse losses, test a likelihood-fitted null location with the rest of the reader frozen.
Do not add a library switch, gene pooling, a new threshold, or a calibration redesign to resolve this review.

## 2. Findings

**1. Fix before release — separate normalization from inference.** Packet §§2, 4.3, 7.6.

*Derived:* multiplying every component's effective length in a locus by the same positive constant cancels in the EM. Changing the denominator alone therefore cannot improve quantification. But `x_b` also determines the null's likelihood and collapsed density; these change object ratios. Likewise, changing a saturation density changes a cap, not merely its unit. The old/new comparison changes the prior, readout, gate and saturation as well as the denominator. The new junction ceiling is only equivalent to the old ceiling when the selected physical saturation levels agree and the touched-object term does not exceed it.

**Test:** run an output-rescaling control, including junction saturation, gDNA, both RNA orientations and zero-support defaults: allocations must remain unchanged. Separately move only the null spike's physical location, keeping basin membership, mass, alternative and physical caps fixed. Only the latter arm can attribute losses to the estimated null location.

**2. Fix before release — the largest basin does not identify the biological null.** Packet §3.2, §5.2b; Figure 3c.

*Derived:* the census measures reliability-weighted, rendered region mass, not unprobed genomic area or fragment count. An abundance argument about the genome does not establish a majority in this census. Even a correctly measured, captured majority may legitimately be the largest basin. *Measured:* the stress landscape's biased level has **56%**, versus **44%** at the lower level. Calling that higher level the unit is harmless in principle; collapsing it and attributing all density differences to shared capture is consequential. Efron's empirical-null construction also relies on identifying assumptions; the two-groups label supplies none for this application. [Efron, 2004, §§1–3](https://www.stat.cmu.edu/~jiashun/Teaching/F08STAT756/Lectures/Efron.pdf).

**Test:** hold two physical rates, supports and object likelihoods fixed; vary only population mass through equal basin weights. Inspect ratio continuity and length ratios, then repeat with genuinely captured objects and with RNA leakage. A discontinuity or loss of a resolved contrast refutes universal null identification. These two biological causes can have identical reader inputs; density alone cannot always distinguish them.

**3. Blocker for adopting the one-sided variant — it erases resolved capture when capture dominates.** Packet §6.

*Derived and reproduced:* with rates **0.001 and 0.1**, strong object evidence, and upper-level mass **60%**, the integrated reader returns **0.01 and 1**, preserving the **100×** ratio. The exact one-sided prior transformation returns **1 and 1**. At upper mass **40%**, both readers return **1 and 100**. Thus the proposed “safety” can discard an arbitrarily well-measured capture contrast solely because population weights change. This survives infinite sequencing depth; it is not sparse-count uncertainty.

**Test:** [the saved counterexample](../../.cache/capture_reference_independent_review/reader_checks.json) already confirms the failure. Any replacement must preserve the ratio across this mass reversal, then pass an end-to-end captured-majority simulation. Fixing the OFF stress case by suppressing this ON signal is not a general repair.

**4. Fix before release — pooling needs a likelihood and defensible membership.** Packet §§5.2a, 7.1; Figure 3a.

*Measured:* `d10 g05 ON` has **140** uncaptured region fragments: peak/truth **0.0464**, pooled/truth **0.8336**. Conversely, `lowg g01 ON` has **267**: peak/truth **0.9292**, pooled/truth **0.00443**. *Derived:* the script selects members using the same observed `k/S` it averages. Truncating a Poisson count below a rate cutoff biases its selected mean downward; including zeros does not remove that selection bias. Under known membership and correct Poisson exposures, 140 total counts imply roughly **8.5%** relative standard error, not a 22-fold ambiguity. Deconvolution and membership error can be much larger.

**Test:** independently vary location estimation and membership: census peak; current hard-selected pool; truth-geometry membership with inferred counts; same membership with true counts. Report count totals and support totals separately. Recovery only after membership replacement implicates selection; additional recovery with true counts implicates calibration. “Sparse” is not itself a diagnosis.

**5. Blocker to unconditional release approval — the OFF stress loss is real, but its remedy is unproved.** Packet §5.2b and §6.

*Measured, independently rerun:* pinned scan, fractional assignment, EM seed 0 and `OMP_NUM_THREADS=1` reproduce **12.4438 → 36.7715%** error; one-sided restores **12.4438%**. There are only **54** true transcript fragments, so the difference is **13.14** additional absolute-count-error units. The moderate companion improves **3.5412 → 2.4862%** under both mode-free variants. This establishes a specific failure, not its prevalence. The stress is capture-OFF and in scope even when calibration reads its strand channel as unstranded.

**Test:** cross the fixed reader with current versus oracle gDNA deposits, then repeat this mechanism across simulation seeds and genome/annotation sizes. If oracle deposits remove the loss, the origin is upstream; if the loss remains, inspect reader and consumer pricing. A release requires either a validated local correction or an explicit acceptance of this measured limitation. A changed golden file alone is not resolution.

**6. Fix before release — distinguish zero counts from no information.** Packet §§4.1–4.3, 5.1, 5.2c; Figures 2 and 4.

*Derived:* at `k=0, S>0`, the likelihood is `exp(−Sρ)`, not flat. Figure 4 stratifies by **true deposited count**, not necessarily the reader's inferred count. Consequently its **×190** versus **up to ×7** comparison does not isolate a prior-only error. The release also explicitly returns **1 at S=0**, unlike its positive-support flat-likelihood limit. *Reproduced:* the two-level fixture gives **1 at S=0**, **6.3096 as S→0+**, and **1 at S=1000**, all at zero count.

**Test:** stratify independently by true capture geometry, inferred count and support; report errors weighted by the object's contribution to component length. Sweep support through zero on RNA-supported pieces. A material allocation jump warrants a fallback-policy experiment; an unweighted zero-count histogram cannot establish that cost.

## 3. The reference question

**Largest basin:** defensible as an empirical approximation when the uncaptured population dominates the actual training measure, is resolved, and nuisance/deconvolution heterogeneity does not generate a larger competing population. The **92** saved positive-truth calibrations with at least **1,000** uncaptured region fragments have peak/truth **0.9545–1.1061**, which supports that regime. It does not establish a universal rule or a production threshold. There are **123** rows overall, **99** with positive truth; these are selected saved calibrations. The fragment tally covers contained unprobed regions, not every uncaptured fragment in the library.

**Sparse null:** for known uncaptured membership `U`, maximize the joint count likelihood:

\[
\ell(\rho_0)=K_U\log\rho_0-E_U\rho_0+\text{constant},\qquad
\widehat\rho_0=K_U/E_U,
\]

where `K_U=Σk_i` and `E_U=ΣS_i`. Zero counts contribute exposure and belong in the calculation. This uses information more efficiently than taking the peak of a sum of separately rendered kernels. With no events, report one-sided uncertainty; the grid floor is not an estimate supported by the data.

Without known membership, the first candidate after a successful oracle is a **reader-only fit of the null location**: freeze current basin masks, mass and alternative, and maximize the mixture's marginal Poisson likelihood over the training objects, using the existing training weights. Posterior null responsibilities give the rate update `Σw_i q_i k_i / Σw_i q_i S_i`. Use soft responsibilities from that likelihood, not `k_i/S_i` cutoffs. Joint inference is not inherently circular; treating selected members as certain is the problem. This candidate adds no empirical constant, but assumes the chosen null component is the right one. It does not fix the stress's biological ambiguity. Fractional deconvolved counts also make this a working likelihood; assess uncertainty by simulation before treating its intervals as calibrated.

Do not add zero inflation merely because counts are small: ordinary Poisson sampling already generates zeros. An extra structural-zero component needs independent predictive evidence and is a later model change. Probe geometry is useful for **oracle diagnostics** here; production use needs an explicit change to the owner's data-only ruling. Gene pooling remains inappropriate for sparse probes.

**Neither end need be located as a mode.** For independently sampled counts from two known groups,

\[
K_C\sim\operatorname{Pois}(\lambda R E_C),\quad
K_U\sim\operatorname{Pois}(\lambda E_U)
\quad\Longrightarrow\quad
K_C\mid K_C+K_U=N\sim\operatorname{Binom}
\left(N,\frac{R E_C}{E_U+R E_C}\right).
\]

Conditioning removes the nuisance scale `λ`; it does **not** remove uncertainty about `R`. For positive large counts, `Var(log R̂)≈1/K_C+1/K_U`. With `K_U=0<K_C`, the likelihood has no finite upper bound on the ratio without additional assumptions. For a spectrum, the corresponding conditional likelihood is multinomial with probabilities proportional to `S_i r_i`; regularization and uncertain membership still matter. A later hierarchical reader can carry a joint posterior over relative rates without locating either mode. It cannot manufacture absent information or distinguish capture from identical gDNA-only nuisance signals. Any new prior must be specified as a model assumption, never fitted to benchmark success.

**Per-library unit:** choose one for numerical stability if useful, then transport it consistently through every consumer. A well-measured captured unit does not rescue a poorly measured uncaptured ratio. Keep the existing unit in 0.8.0; switching between the old and new inferential rules by library would add model selection, not merely change units. Null/alternative groups here remain a capture prior within gDNA, not additional biological populations.

## 4. The one-sided null

**Refuse the present formulation.** The physical inequality is capture yield relative to the **true uncaptured baseline** ≥1. It says nothing about yield relative to an arbitrary largest basin. Values below that basin can encode real uncaptured objects. Moreover, §6 changes the prior and posterior membership; it is not just clipping final answers at 1.

The measured trade is mixed. These are per-stratum error percentages; deltas are percentage points:

| Stratum | Integrated | One-sided | Delta |
|---|---:|---:|---:|
| OFF extreme stress | 36.77 | 12.44 | −24.33 |
| lowg κ.99 OFF, contaminated | 7.70 | 6.96 | −0.74 |
| fle200 κ.70 ON | 9.29 | 8.77 | −0.53 |
| bind100 κ.99 ON | 7.85 | 8.72 | +0.88 |
| flgl κ.99 / κ.70 ON | 9.05 / 10.49 | 9.37 / 10.85 | +0.32 / +0.35 |
| Ladder κ.99 ON | 3.87 | 3.90 | +0.02¹ |
| RNA-short / RNA-long gap κ.99 ON | 13.50 / 5.17 | 13.50 / 5.17 | 0.00 |

¹ Delta calculated before rounding. At the packet's descriptive ±0.3-point band, in-scope better/same/worse versus foundation is **7/72/17 → 6/73/17**. This band is not an acceptance threshold or a significance test. Deferred κ.50 ON results remain reported in the [audit](../../.cache/capture_reference_independent_review/audit.json).

The bind100 cost's isoform mechanism is still **conjectural**: changing 1,080 regions to 1 does not identify which length ratios caused the loss. Replay those object changes through fixed consumers, inspect the affected transcript/gDNA lengths and then allocations. A future physical one-sided model would need an identified baseline and a way to separate shared capture from gDNA-only nuisance; that is not this patch. There is no one-sided real-library arm in the supplied bundle.

## 5. Objects without evidence

For 0.8.0, keep the **geometric prior prediction for positive-support objects with effectively flat likelihood**, plus the existing explicit zero-support convention pending finding 6's test. A known-unprobed object should read the uncaptured level. An unknown object with no informative observation cannot be assigned that knowledge merely because the reviewer knows its simulation truth.

If exchangeability with the chosen prior holds, a prior prediction is the coherent Bayesian answer. For a fixed target and posterior, geometric prediction minimizes squared log error; arithmetic prediction minimizes squared linear error. Neither loss is automatically the downstream transcript-allocation loss. The reduction from ×190 to at most about ×7 accompanies the geometric candidate; it does not isolate its readout from the other changes. A zero-count override would also suppress captured objects with the same weak evidence. No support cutoff or new pseudo-count is justified by these plots.

## 6. Experiments to run first

**Common protocol for every candidate:** freeze scan/deconvolution inputs where the arm permits; scan pinned, fractional assignment, `OMP_NUM_THREADS=1`, EM seed 0; one mechanism at a time. Use paired simulation seeds for uncertainty, never to select a favorable result. Read transcript error and component-length ratios per stratum, separating zero-gDNA rows. Promote any candidate to **ladder and both gap panels**, all **13 robustness panels**, and the **three plasma libraries plus VCaP/matched RNA-only**. Report unstranded ON separately. Real libraries test behavior, conservation, runtime and available truth; they do not choose constants or certify accuracy without truth.

1. **Oracle bulk location, with a normalization control.** Start with the seven sparse positive-truth calibrations, especially `d10 g05`, `d10 g01`, `lowg g001` and the apparently correct-peak `lowg g01`. A0 only rescales all outputs and ceilings: error must not move. A1 inserts the oracle density into the **null spike's likelihood and readout**, consistently rebasing saturation and fallbacks; retain original masks, prior masses, slab, absolute cap locations and calibration. Use the simulator's expected uncaptured yield on matching supports where available; otherwise identify `Σtrue counts/Σtrue support` as a noisy oracle and assess it across seeds. Recovery prices the spike-location mechanism. Non-recovery refutes its sufficiency, not every possible baseline defect. Add oracle membership separately. For the old reader, separate cap-location override from gate bypass: a reference supplied behind a closed gate cannot price the sparse-capture failure.

2. **Stress causality and captured-majority preservation.** Use finding 5's deposit oracle and replicated stress simulations, alongside finding 3's mass-reversal test. These decide whether a local reader repair can address the regression while preserving real capture. Do not turn the oracle into a production switch. If only upstream correction works, document that boundary and obtain an explicit release-risk decision; do not redesign calibration in this release review.

3. **Sparse-null estimator A/B, conditional on experiment 1.** Hold membership and all consumers fixed; compare peak location with the marginal-likelihood location fit in §3. Separately test membership uncertainty. Recovery of location and transcript errors across sparse and well-supported strata supports adoption; a better rate without better allocations does not. Include OFF controls and grid-range sensitivity: estimates tracking an arbitrary lower grid edge fail this test.

4. **Weak-evidence and saturation attribution.** Apply finding 6's length-weighted analysis, then perturb basin partitioning under simulated resampling on the split-level `odg05` and RNA-long cases. Hold the junction rule fixed. If losses follow a cap split rather than baseline placement, a reference change will not repair them. Any eventual fallback or cap change requires its own full A/B; no spatial smoothing or gene pooling is inferred from this diagnostic.

5. **Direct-ratio prototype, later.** First test the conditional likelihood on known Poisson counts with fixed truth ratios while varying exposure and group mass; check uncertainty coverage and zero-event limits. Then add uncertain membership and deconvolved inputs separately. It earns a release comparison only after passing those identification checks. Shared uncertainty must be propagated jointly; independent per-object error bars would miss a common uncertain baseline.

## 7. Answers to packet §8

**1. Null identification.** Largest mass is a useful conditional approximation, not a biological identifier. Its success requires dominance in the training measure and tolerable nuisance variation. Replace the unconditional claim now; benchmark membership alternatives before changing the rule.

**2. Sparse null.** Use a joint exposure-aware count likelihood, with zeros included and membership uncertainty retained. The present hard-selected pool is biased. Geometry-defined oracle membership separates that problem from calibration error; zero inflation is not justified by sparsity alone.

**3. Unit and ratio.** Either unit is acceptable if transported consistently. Conditional ratio inference eliminates an absolute scale but retains the information limit from scarce uncaptured events. Per-library normalization is harmless; per-library model switching is not equivalent.

**4. Below-bulk levels.** Keep them in the current candidate. Below the largest basin is not necessarily below the true baseline. The one-sided variant trades a stress rescue for measured ON losses and a derived captured-majority failure.

**5. No evidence.** Use the stated prior prediction for an unknown object, distinguishing weak positive support from the explicit zero-support fallback. A universally correct read of 1 requires uncaptured membership information the reader may not possess.

**6. Stress.** The originating ambiguity is upstream; the reader amplifies it through unequal component lengths. “Lowest well-populated level” introduces a new selection rule and excludes unresolved cases without resolving causality. Establish the deposit-oracle result first; preserve the three-population model.

**7. Release.** Prefer the integrated two-sided reader, conditionally. Require the oracle-bulk and stress disposition before sign-off; reject the current one-sided patch and old capture gate as universal remedies. None of the supplied evidence establishes a safe, unconditional replacement baseline estimator.
