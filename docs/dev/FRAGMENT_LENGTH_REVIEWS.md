# The fragment-length reviews

*Sandbox record (`docs/dev/`): the external reviews of the fragment-length work, verbatim, oldest first. The plan
they shaped is [`FRAGMENT_LENGTH_PLAN.md`](FRAGMENT_LENGTH_PLAN.md); the problem statement the first one reviewed is
[`FRAGMENT_LENGTH_REVIEW.md`](FRAGMENT_LENGTH_REVIEW.md). The plan versions the later reviews read are kept outside the
repo, in `~/Downloads/rigel_runs/prototypes/2026-10-02_fl_review/`, as `plan_round1.md` to `plan_round3.md`.*

## The first review, verbatim

It reviewed `FRAGMENT_LENGTH_REVIEW.md`.

> **I recommend testing an opportunity-aware background model first, using the existing deposits.** That is the smallest plausible repair for false capture. Length-dependent capture then needs its own correction; a softer global length likelihood cannot fix both problems.
>
> The review’s measurements are compelling, but two mathematical claims need correction:
>
> - The count fraction is a valid coordinate when opportunities differ. The identity is
>   \(\operatorname{logit}f_g=\log(\rho_g/\rho_r)+\log(E_g/E_r)\).
>   Making the reference symmetric in density instead would change the prior.
> - With count and opportunity fixed,
>   \(\operatorname{Var}(\log\rho_g)=\operatorname{Var}(\log f_g)\).
>   Changing the variance test’s name would do nothing. The missing uncertainty concerns the opportunity estimate and the observation model.
>
> The August history also shows that “make length fuzzier” is insufficient by itself. Nearly flat likelihoods became endpoint estimates with nonzero precision. The current [factor-precision calculation](/Users/mkiyer/proj/rigel/src/rigel/native/solve_kernel.cpp:196) still has that near-flat discontinuity.
>
> Here is my step-0 plan, before building anything.
>
> 1. **Use the existing full-depth failure as the development substrate.**
>    Its cached calibration loop is documented at roughly 5–12 seconds per condition. I verified that the failing condition’s scan and gDNA-oracle payload digests match the current tree. This avoids thinning away the false mode.
>
>    Next, create an isolated test-chromosome arm with the actual failing distributions: RNA **75±20**, gDNA **250±60**. The existing arm uses **100±50 / 250±150**; its much broader gDNA short tail makes it a different test. Require the new arm to reproduce the intermediate failure—small opportunity, exaggerated density, false enriched mode—before calling it a reproducer. Preserve the current nascent model.
>
> 2. **Re-record baselines from `bdfd8709`.**
>    Start with the failing stranded/OFF condition, its unstranded and capture-ON twins, and the capture-reference-disabled diagnostic. Then cover the equal-length ladder, both gap directions, equal-length controls, and g00/g05/g98.
>
>    Record calibration and prior errors against certified truth; gDNA error by opportunity; landscape membership and capture reference; transcript, gene and pool errors. Read `ruler_vs_truth` by probe class beside every capture-ON arm. Pin scan threads, fractional assignment and `OMP_NUM_THREADS=1`. The quoted 42.1% and 3.5% remain historical results until re-recorded.
>
> 3. **Derive and test the background model without changing fragment accumulation.**
>    For an exon, let expected counts be \(x=\rho_gE_g\), \(y=\rho_rE_r\), with \(N\sim\mathrm{Poisson}(x+y)\). The external background constrains \(\rho_g\); RNA remains locally free. Thus a short exon with \(E_g=0.09\) and background density 0.05 expects **0.0045 gDNA fragments**, continuously rather than through a one-position cutoff.
>
>    The existing [factory row](/Users/mkiyer/proj/rigel/src/rigel/native/solve_kernel.cpp:179) evaluates a continuous extension of \(\mathrm{NB}(Nf_g)\). That is an approximation worth testing, not an exact observed-count likelihood.
>
>    A compact mathematical comparator is available. With background density \(\rho_g\sim\mathrm{Gamma}(a,b)\), RNA’s Jeffreys rate reference, and \(\beta=b/E_g\), integrating total intensity gives the following density on the existing logit grid:
>
>    \[
>    p(\lambda\mid N,\text{strand})
>    \propto L_{\mathrm{strand}}(f_g)
>    \frac{f_g^a(1-f_g)^{1/2}}
>         {(1+\beta f_g)^{N+a+1/2}}.
>    \]
>
>    This replaces the corresponding reference/background factors. Its parameters come from the measured background and its uncertainty. Compare it with the existing factory extension in separate experiments.
>
> 4. **Make the falsification gates address the actual risks.**
>    Before prototyping, establish failing tests for low-opportunity over-attribution and false capture. Require improvement across opportunity bins, preservation of true capture, flat ladder performance, and no increase in g00 false gDNA.
>
>    Two additional gates are essential:
>
>    - A background-driven exon must not acquire apparent independent evidence merely by training the landscape that subsequently reinforces it. The [current evidence and training path](/Users/mkiyer/proj/rigel/src/rigel/calibration/calibrate.py:248) makes a blanket factory extension vulnerable to this.
>    - The capture slab must be justified separately. A lower endpoint near background does **not** collapse a slab whose upper endpoint is the RNA-rich exon’s total density. The proposed OFF-collapse argument is incomplete.
>
>    Prototype C++ outside the main tree, vary one mechanism per A/B, and deliberately break each successful fix to verify its gates. Background evidence may remain at equal lengths; only added **length discrimination** must vanish as the local laws coincide.
>
> 5. **Preserve the right observations if length evidence is still needed.**
>    My preferred additional payload is **length-bin × genome-strand counts** for contained regions and unspliced boundaries, keeping certified-spliced evidence separate. For length-dependent conserved attribution, boundaries also need conserved mass in matching bins. A crossing count and a conserved fragment share answer different questions.
>
>    Store observations, not fitted origins or likelihood scores. Use the resolved molecular length already defined by the accumulator. Derive each bin’s opportunity from the exact deposit population; fit strand and length jointly or conditionally. Multiple boundary crossings must not become independent fragments.
>
>    “Soft” should mean that uncertain laws produce weak evidence throughout the solve—including messages and landscape admission. It must cover shared model error under capture, not merely finite histogram sampling. Missing gDNA training evidence must remain unmeasured.
>
> For capture, the required quantity is component-specific, length- and placement-weighted opportunity. A single global correction cannot generally recover probe reach. I would establish that mechanism with truth-fed capture opportunities before choosing additional production storage or a response model.
>
> No files were changed, prototypes built, or benchmarks rerun. This is the requested pre-build plan.

## The second review, verbatim

It reviewed `plan_round1.md`, and its section numbers refer to that version.

> **The Gamma algebra is sound. The claimed echo guarantee is false on the shipped grid.** I checked the cited code, re-derived the formulas, and ran small calculations through the existing native slot solver. Nothing was edited or built.
>
> 1. **Agree on the Gamma derivation; unsure about the landscape plug-in’s accuracy.**
>    The Jacobians and exponent \(N+a+\tfrac12\) are correct, including the proposed \(a=\tfrac12\) term. For a log-density prior \(q\), my derivation gives the exact expression
>    \[
>    p(\lambda\mid N)\propto L(f)\sqrt{1-f}\,
>    \mathbb E_{T\sim\Gamma(N+1/2,1)}
>    [q(\log(fT/E_g))].
>    \]
>    Replacing that expectation by its value at \(T=N\) is an approximation; the integration kernel’s mode is \(N-\tfrac12\). The smallest check is numerical quadrature versus the plug-in for sparse counts and narrow or multimodal landscapes.
>
> 2. **Disagree with the guarantee; \(a=\tfrac12\) is a defensible experimental choice, not a derived biological dispersion.**
>    It imposes density CV \(=\sqrt2\), regardless of measured heterogeneity. More decisively, using the [shipped grid](/Users/mkiyer/proj/rigel/src/rigel/config.py:288), reference plus proposed term **alone**, at \(\mu=0.0045\), produces native `Var(log f)` **0.709 at \(N=30\)** and **0.106 at \(N=300\)**—both “located.” Widening the bracket from 10 to 40 gives **4.968 and 4.938**: truncation manufactured the precision. The inverse-variance addition argument also fails for these non-Gaussian distributions; even without truncation, it supplies no universal echo bound.
>
> 3. **Unsure about admissibility under the ruling; disagree with §4.7’s justification.**
>    This is a located density prior conditioned on total count, mathematically different from the refuted fixed-share prior in [EQUATIONS §9c](/Users/mkiyer/proj/rigel/docs/EQUATIONS.md:787). However, its precision does **not** keep growing with \(N\): at fixed \(\mu\), `Var(log f)` approaches \(\pi^2/2\). If sequencing depth increases, \(\mu\) increases too, so the claimed \(1/N\) share decrease need not occur. A tiny depth-scaling experiment should distinguish increasing RNA alone from increasing the whole library; accepting this as a population-prior change still requires reconciling the ruling explicitly.
>
> 4. **Agree: pass 0 only for the first experiment.**
>    Multiplying the fitted population prior by the same background prior again would assert an additional constraint. But removing the term after pass 0 does not remove its influence: the [refit loop](/Users/mkiyer/proj/rigel/src/rigel/calibration/calibrate.py:448) learns from its answers. Record training membership and density modes at every refit; recovery after three refits is an experimental question.
>
> 5. **Agree with the squeeze estimate only as a local approximation.**
>    I derive \(-c/(2P+c)\) as one Newton step when \(N\gg c\). Under the same Gaussian log-count approximation, solving the displacement fully gives \(-W(c/(2P))\): at \(c=100,P=25\), **−0.853**, versus **−0.667** in the table. I would try the pooled-exon centre first as the cheapest diagnostic A/B, with the evaluated exon excluded from its estimate. It is an opportunity-weighted, potentially biased mean; it neither equals background automatically OFF nor protects a rare highly enriched population.
>
> 6. **Agree that the declared landscape semantics imply an extra prior tilt; disagree that the proposed reconciliation is automatically correct.**
>    [The kernel](/Users/mkiyer/proj/rigel/src/rigel/native/psi_kernel.h:182) unquestionably adds both halves, while [the landscape](/Users/mkiyer/proj/rigel/src/rigel/calibration/landscape.py:48) declares itself a population density in log-rate. Under that interpretation, the extra \(f^{1/2}\) reweights \(q\) by \(\rho_g^{1/2}\); it is not a necessary Jacobian. Continuing slope \(+\tfrac12\) below the grid repairs the lower tail, but continuing it toward infinite density does **not** define a proper population prior—even though an individual slot’s plug-in density is bounded by \(N/E_g\). Test tail completion and bracket/grid-extension invariance together; dropping the half while retaining constant tails would confound the comparison with impropriety.
>
> 7. **Agree to a separate exact-form intron experiment; disagree that prior curvature is independent own evidence.**
>    Preserve the existing background shape when comparing integration against `NB(Nf)`; otherwise two mechanisms change. The [current kernel](/Users/mkiyer/proj/rigel/src/rigel/native/solve_kernel.cpp:534) counts factory-row precision as own evidence, so the provenance concern is real. External background measurements can inform introns, but reproducing their constraint at many introns does not create independent measurements of the background. Start with a small intron chain and vary background uncertainty and repeated use separately.
>
> 8. **Disagree with calling the proposed subtraction derived information.**
>    It fixes continuity, but inverse-variance subtraction and clipping remain a heuristic. An informative bimodal row can have greater variance than the flat distribution and receive zero precision; the result also depends on the bracket. Prefer information derived from the observation model. The smallest discriminator is a collection of flat, weakly tilted, concentrated and bimodal rows evaluated across bracket widths—not another genome run.
>
> 9. **Agree with the original `B_match`; disagree with the ruler-only causal test described here.**
>    Feeding simulator yields into an instrument scored against those yields mostly verifies plumbing. The [specified experiment](/Users/mkiyer/proj/rigel/docs/ISSUES.md:289) also requires `B_transfer1`, matched scorer laws, the gDNA denominator on the same scale, and downstream quantification controls. For production, I favour optional probe geometry plus a measured capture response: it represents reach explicitly. The response and behaviour without probe information remain to be established.
>
> 10. **Agree with the rejection bars, with additional numerical and provenance gates.**
>     Reject the proposed implementation if its apparent success comes from clipping posteriors against the grid, suppressing genuine capture, increasing g00 contamination, or teaching the landscape its own imposed background. Also reject a transcript-only win that leaves the low-opportunity calibration error intact. The “term alone never locates” gate already fails under the existing numerical settings; that must be resolved before panel A/Bs.
>
> Additional errors or unsupported claims in the document:
>
> - **§4.4:** the Gamma statement is asymptotic. For finite \(N\) and very small \(\mu\), \(f/(2\mu)\) approaches a beta-prime distribution with shapes \(\tfrac12,N+\tfrac12\); its log-variance is \(\psi_1(\tfrac12)+\psi_1(N+\tfrac12)\).
> - **§4.3 and Appendix B:** the contained-opportunity formula applies to regions. [Boundaries use crossing opportunity](/Users/mkiyer/proj/rigel/src/rigel/calibration/region_geometry.py:193). Small global TV error also does not bound relative error in a tiny opportunity.
> - **§8 numbers:** the review’s table gives oracle-gDNA unprobed RNA-short error **+1.62**, not **+2.96**; the latter is shipped. RNA-long unprobed error is **+2.13**, so “within ±0.1” needs restriction to the partially probed classes. [Source table](/Users/mkiyer/proj/rigel/docs/dev/FRAGMENT_LENGTH_REVIEW.md:264).
> - **§7 citation:** the supposed Python reference [calls the same native cube](/Users/mkiyer/proj/rigel/tests/calibration/_psi_reference.py:236); it is not independent mathematical corroboration.
> - **§2 and §9:** the ladder and gap arm share the annotation, so “the ladder’s exons are long” is misleading. Uniform Poisson simulation also does not guarantee the finite-sample [dispersion estimator](/Users/mkiyer/proj/rigel/src/rigel/calibration/density_deconv.py:106) returns infinity.
>
> The main misreading of my previous review is treating the comparator as an endorsement of fixed shape \(a=\tfrac12\). I proposed comparing models whose uncertainty was justified; I did not establish that shape, the echo guarantee, or prior-only causation.
>
> **I would first resolve the numerical counterexample**, then compare exact integration and the existing factory on small known-truth examples. Add simulated density heterogeneity and opportunity error, and track background-induced landscape admission. Those checks can decide whether this is a viable first mechanism before expensive runs.

## The third review, verbatim

It reviewed `plan_round2.md`, and its section numbers refer to that version. It answered: "Do you agree that this plan
is ready for implementation? Is there a simpler or more elegant design? Considering the entire tool workflow, is there a
way to solve this definitively?"

> **Not yet. The revised plan is ready for focused experiments, but I would not approve it as an implementation specification.** Mechanism two has a sound derivation. Mechanism one is an adaptive exposure guard whose claimed guarantee remains unproved—and, in its general form, is false.
>
> I checked the selection, landscape estimator, solver, accumulator and downstream scoring code. The numerical counterexamples below are my own derivations.
>
> 1. **The expected-count cut does not remove selection bias.** The current selector requires a located composition ([calibrate.py:248](/Users/mkiyer/proj/rigel/src/rigel/calibration/calibrate.py:248)). Adding `ρ_ref E_g ≥ 1` leaves that outcome-dependent selection in place.
>
>    Even in the simpler case `K ~ Poisson(μ)`, retaining `K ≥ 1` gives
>    \[
>    E[K\mid K\ge1]=\frac{\mu}{1-e^{-\mu}}.
>    \]
>    At the proposed boundary `μ = 1`, the selected mean is **1.582**, not 1. There is no transition to unbiased measurement at one expected fragment.
>
>    More importantly, a deconvolved gDNA count is not an observed Poisson count. With 1,000 RNA fragments and 1% strand error, strand noise alone gives the unconstrained gDNA count estimator an SD of approximately **6.4 fragments**. The plan’s “one or two fragments” argument does not bound that noise. If the new rule instead *replaces* the location requirement, it needs a separate argument preventing uninformative, prior-determined estimates from training.
>
> 2. **I would remove the pooled-located-exon fallback from the proposed fix.** A ratio of sums addresses unstable small denominators; it does not correct selection of the numerator.
>
>    For example, take equally exposed objects with true density `0.1`, exposure `1`, and retain positive counts. Their pooled density converges to **1.051—10.5 times the truth**. Using that estimate to decide admission can admit the very population that inflated it. The assertion that this fallback cannot put noise above the enriched mode is unsupported.
>
>    The census is useful, but surviving membership alone cannot certify preservation of capture. It must also measure reference displacement and truth error through the refits.
>
> 3. **I agree with replacing the extra gDNA Jeffreys half, with two implementation qualifications.** Under the declared interpretation of `q` as a population density in log rate, adding the full reference introduces an additional `√ρ_g` tilt ([psi_kernel.h:182](/Users/mkiyer/proj/rigel/src/rigel/native/psi_kernel.h:182)).
>
>    The proposed lower-tail continuation fixes the integrability problem **only if evaluated using the unclipped log fraction**. The current evaluator clips `f` at `10⁻¹²` ([psi_kernel.h:65](/Users/mkiyer/proj/rigel/src/rigel/native/psi_kernel.h:65)); continuing a slope after that clipping eventually produces a flat tail again.
>
>    The upper-grid argument works for the current `T = N` evaluation, provided every consumer is covered. It does **not** extend to exact integration over `T`, which has unbounded support. That version needs an explicit upper-tail definition. I would prototype this correction independently before adopting a new training cut.
>
> The simpler design principle is: **retain uncertainty as evidence instead of turning an uncertain estimate into a new observation.** Currently, calibration passes `f_g M` to the landscape, which treats it as a Poisson count ([calibrate.py:255](/Users/mkiyer/proj/rigel/src/rigel/calibration/calibrate.py:255), [landscape.py:138](/Users/mkiyer/proj/rigel/src/rigel/calibration/landscape.py:138)). That is the deeper statistical seam.
>
> If training needs redesign, the principled target is a population fit using the observation likelihood:
> \[
> \sum_o \log\!\int L_o(D_o\mid z)\,q(z)\,dz,
> \qquad z=\log\rho_g,
> \]
> with RNA nuisance quantities handled explicitly. A likelihood constant in `z` contributes no information about `q`. This requires preserving genuinely independent evidence, rather than recycling posterior beliefs or counting shared messages repeatedly.
>
> This is **a candidate to prove**, not an endorsement of resurrecting the old estimator. The refused `row_kernel` and `estep_all` experiments matter: averaging normalized likelihoods is different from this objective, and the counted-minority failure remains a mandatory falsifier ([ISSUES.md:2563](/Users/mkiyer/proj/rigel/docs/ISSUES.md:2563)).
>
> For the whole workflow, I see a more promising route than estimating an uncertain global gDNA length distribution and giving it another calibration role: **measure the relative RNA:gDNA length law within objects where capture cancels.** This is already identified in [ISSUES.md:278](/Users/mkiyer/proj/rigel/docs/ISSUES.md:278).
>
> For a one-RNA-strand object, when both origins share the same placement and capture factor at length `b`,
> \[
> \operatorname{logit}q_{ob}=\alpha_o-\log r_b,
> \]
> where `q` is the gDNA fraction and `r_b` is the RNA:gDNA length-law ratio. Strand counts observe that mixture; the common capture factor cancels. This avoids requiring an accurately reconstructed uncaptured global law for that particular measurement. The cancellation must be established for each placement class; it cannot simply be assumed across splice junctions.
>
> That suggests a focused payload extension:
>
> - Preserve **joint length × strand counts** per object, with certified spliced evidence distinguished.
> - Preserve the corresponding **conserved-mass sums** for boundary attribution.
> - Derive opportunities for precisely those deposits.
> - Supply fragment scores and their effective-length integrals from the same component-specific observation model.
>
> The current region payload retains only two strand totals ([accumulator.h:97](/Users/mkiyer/proj/rigel/src/rigel/native/calibration/accumulator.h:97)). Global length pools cannot recover the lost association. Uncertainty in the relative length law should be propagated through the model; that gives the requested softness without inventing a likelihood-temperature constant.
>
> **What I would do first:** build a small truth-controlled falsifier varying RNA depth at fixed gDNA opportunity, alongside genuinely enriched short exons. Use it to judge the cut and any replacement training method. Separately prototype the prior correction, and run the specified `B_match`/`B_transfer1` experiment before committing to a capture implementation. Keep the background prior deferred.
>
> A definitive repair can make the measurement model coherent and eliminate these artificial sources of information. It cannot guarantee identification where the observations do not distinguish origins or where unknown capture and abundance remain confounded. The design should preserve that uncertainty rather than manufacture a located answer.

## The fourth review, verbatim

It reviewed `plan_round3.md`, and its section numbers refer to that version.

> This is a substantially better research plan. I agree with starting Stages 0, 1 and 3, with the corrections below. Part 2 still has an unresolved statistical assumption; Part 3 still has an identifiability gap. Neither should yet be treated as an implementation specification.
> 1. Own-evidence likelihood: partly agree.
>    The no-RNA and one-RNA-strand expressions are correct under the stated Poisson model and RNA reference. But integrating against Y⁻½ does not produce reference-free evidence, and the total-only expression is not merely an upper bound. For example, my derivation gives
>    \[
>    \int_0^\infty \mathrm{Pois}(1\mid X+Y)Y^{-1/2}\,dY
>    =\sqrt{\pi}e^{-X}(X+\tfrac12),
>    \]
>    which peaks at X = ½, despite there being no strand information. I verified this by quadrature. The RNA reference itself favours a positive gDNA intensity here. This is the first new falsifier to add.
>    Also, the both-stranded row discards useful strand evidence: a 90:10 split constrains gDNA even when both RNA strands are allowed. Integrate the RNA strand allocation, consistent with the solver’s existing treatment, rather than collapsing immediately to N ([simplex_logodds.py (line 37)](/Users/mkiyer/proj/rigel/src/rigel/calibration/simplex_logodds.py:37)).
>    Excluding messages, the factory row and overlapping boundary observations from this first population fit is sound. But call the result a likelihood marginalized under an explicit RNA reference, and test that reference’s influence.
> 2. Full fit versus frozen variant: use the full fit as the candidate; reject the frozen variant as specified.
>    The full fit has a defined objective. The frozen variant reintroduces selection and uses a normalized likelihood as a population contribution—the same conceptual problem that motivated this redesign.
>    There is a sharper mathematical issue: for RNA-admitting regions, generally L(z) → L₀ > 0 as z → −∞. Consequently, normalizing that likelihood over log density requires an artificial lower boundary. Its variance, and therefore “located” status, depends on that boundary.
>    Population 3 is necessary, but vary minority size, exposure and RNA depth. Require recovery only where the observations can distinguish the minority. Many individually weak regions can jointly locate a population; requiring every capture-reference member to have a located individual likelihood would discard that collective information. The proposed “mostly in the basin” and inherited √n rules therefore need their own justification.
> 3. Fit once before the sweep: agree conditionally.
>    With the likelihood inputs and fitted prior fixed, this simplification is sound. The current code resets the beliefs before every refit sweep; repeating the same prior and grid repeats the same calculation ([calibrate.py (line 448)](/Users/mkiyer/proj/rigel/src/rigel/calibration/calibrate.py:448)).
>    What remains unproved is whether the reduced evidence produces an adequate prior. Losing message-derived training information affects unstranded capture-OFF too, not only the deferred stratum. Include that comparison explicitly.
>    The smallest check is to compare one sweep against repeated, fully reset sweeps using an identical fixed prior. Then assess the new prior’s accuracy separately. Removing redundant sweeps should not be bundled with an unmeasured loss of evidence.
> 4. Part 1 implementation: mostly agree, but amend the domain and convergence gates.
>    The unclipped axis and lower-tail slope are correct. However, “no consumer reads above the grid” is not currently guaranteed: the grid domain excludes opportunities at or below 10⁻⁹, while the native arm evaluates opportunities down to its 10⁻¹² floor. A positive-count consumer with E_g = 10⁻¹⁰ can therefore read above the grid ([calibrate.py (line 240)](/Users/mkiyer/proj/rigel/src/rigel/calibration/calibrate.py:240), [psi_kernel.h (line 72)](/Users/mkiyer/proj/rigel/src/rigel/native/psi_kernel.h:72)). Make the domain agree with the actual consumers.
>    Also inspect [required_logodds_window (line 72)](/Users/mkiyer/proj/rigel/src/rigel/calibration/landscape.py:72): it clips the fraction and only ensures that the bracket reaches the landscape floor, not that it includes the completed tail.
>    Replace unconditional “invariance at 10, 20, 40” with convergence after controlling omitted tail mass and moments. Even reference-only Jeffreys gives Var(log f) ≈ 2.804 at L = 10, versus 3.290 at L = 40, by continuous quadrature. Properness alone does not make a short bracket accurate.
> 5. Stage 3: agree. Part 3: promising, but the production chain is incomplete.
>    Run the specified matched-table/integral experiment unchanged ([ISSUES.md (line 289)](/Users/mkiyer/proj/rigel/docs/ISSUES.md:289)). A partial improvement would establish a contributing mechanism; failure to close the entire gap would not prove that the mechanism lies elsewhere.
>    Within-object cancellation is valid for shared placements and a shared capture response. But relative odds do not identify both absolute length laws and capture. Even the ratio equation has the invariance
>    r_b → a r_b, α_o → α_o + log a. An independently established RNA law can fix that scale on shared support; a capture-biased spliced histogram cannot automatically supply that law. Stage 3 cannot determine how to estimate the response without a BED.
>    A simplification: prototype with sparse, exact integer lengths. Capture cancellation need not survive pooling lengths into a bin when r varies within it. Introduce compression only after measuring its error. Preserve conserved mass jointly with length and strand wherever attribution uses both.
> 6. The bench is missing several decisive cases.
>    First, add the analytic cases above: N = 0/1, κ = ½, both RNA strands, zero opportunity, zero RNA, and completely uninformative likelihoods. An unidentified fit must not publish its initialization as a learned population.
>    Second, test pooled zero-count resolution. For zero-count gDNA-only regions,
>    \[
>    L(\rho)=e^{-\rho\sum_o E_o}.
>    \]
>    The informative scale is 1/ΣE, not 1/max(E). Reusing the current landscape grid therefore excludes densities the pooled observations distinguish ([landscape.py (line 118)](/Users/mkiyer/proj/rigel/src/rigel/calibration/landscape.py:118)). Specify zero-rate handling and lower-support convergence before fitting q.
>    Third, use repeated realizations, held-out loci, exposure correlated with enrichment, uncertain strand specificity, and oracle-prior controls. These separate sampling error, model misspecification, self-influence and estimator failure. For Part 3, add disjoint length support and capture variation within proposed bins.
> 7. Rejection criteria: agree with the intent; make them more precise.
>    Reject a candidate that learns reference-induced composition from blind regions, creates locations through numerical boundaries, erases a detectable enriched minority, increases in-scope truth error, or succeeds only by disabling genuine capture.
>    Reject the frozen workaround if it is needed to conceal those failures. Do not require recovery of an observationally invisible minority.
>    Two wording changes matter: “displaces the capture reference” must mean away from truth, and “any capture-ON row” must respect the deferred-stratum ruling. Also, the guard remains a comparison arm; failed candidates do not establish it as an acceptable production fallback.
> I would simplify the next revision by dropping the frozen variant and SQUAREM initially, retaining ordinary EM with a checked objective, and putting the RNA-reference and zero-support falsifiers ahead of the large benchmark work. Start the likelihood prototype with structurally gDNA-only and one-RNA-strand stranded regions; add other classes only after their contribution is derived. Structural eligibility avoids the rejected outcome threshold and keeps the first experiment interpretable.

## The fifth review, verbatim

It reviewed the fourth version of the plan. Its `plan:38` is that version's domain-floor bullet in Part 1, and its
`plan:71` that version's zero-support rule in Part 2.

> Yes—go ahead with Stages 0, 1 and 3, and the bench-only Python prototype. The plan is improving, and its scope is now appropriate. The prior arithmetic is a correction; the training redesign remains a candidate.
> There is one substantive qualification to profiling, plus two small implementation clarifications.
> 1. Profiling gives the intended local curves, but profiling inside the population mixture can still manufacture population structure.
>    The count-only profile is indeed flat for X ≤ N. For the stranded case, the stationary equation is quadratic in Y; the log-likelihood itself is not quadratic.
>    However, the proposed population objective uses
>    \[
>    \sum_z q_z\,\max_Y L_o(z,Y),
>    \]
>    whereas profiling a region’s RNA amount in the population model gives
>    \[
>    \max_Y\sum_z q_z\,L_o(z,Y).
>    \]
>    These generally differ. The first allows a separately optimized RNA amount for every candidate gDNA density. Consequently, this is more than ignoring uncertainty at tiny counts.
>    I checked a concrete null case by enumerating the count outcomes: pure RNA, expected RNA count 10, strand specificity 0.99, and true gDNA zero. Under the proposed objective, assigning 0.1% population weight to an alternative with expected gDNA count 1 improves the expected objective by approximately 0.0000208 nat per region. That improvement persists as the population grows. The bias is small in this example; its practical effect needs the bench with structural anchors.
>    Add this cheap population-null check before the larger prototype. Use the shared-Y expression above as a comparator: it requires one scalar nuisance per region and no RNA quadrature. Neither approach should be declared adequate before testing.
>    Also, profiling removes the RNA integral from training. It does not settle the solver’s separate T ≈ N approximation.
> 2. One shared opportunity floor needs one shared consumer domain.
>    At [plan:38 (line 38)](/Users/mkiyer/proj/rigel/docs/dev/FRAGMENT_LENGTH_PLAN.md:38), specify that the grid is constructed from the same floored opportunities and masses for every actual arm consumer. Merely making the two epsilon constants equal does not suffice: a domain filter can still exclude a slot that the evaluator processes after flooring.
>    The decisive gate remains: enumerate actual consumer coordinates and verify that none exceed the represented upper domain.
> 3. 1/ΣE is a resolution scale, not a justified hard lower boundary.
>    For pooled zeros, the likelihood at that rate is e⁻¹, while its supremum is at zero. Substantial likelihood remains below the proposed floor.
>    At [plan:71 (line 71)](/Users/mkiyer/proj/rigel/docs/dev/FRAGMENT_LENGTH_PLAN.md:71), check convergence of predictions and consumed quantities as the lower support extends. Checking the fitted mass at the floor is insufficient: it can remain 100% while the implied rate keeps changing. This does not require introducing a zero atom.
> The capture-reference readout and unstranded capture-OFF adequacy are correctly left open. I would measure both before adding junction evidence or choosing a new membership rule.
> I agree with postponing the broader test programme. The profiling null check above belongs in the first tier because it takes seconds and directly tests the proposed objective. Start there, alongside the prior-arithmetic prototype and the specified capture experiment.
