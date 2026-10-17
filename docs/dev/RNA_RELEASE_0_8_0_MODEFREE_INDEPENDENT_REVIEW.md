# Independent review: the mode-free capture reader

2026-10-13 · Foundation `3e7e0fbc`; frozen `site_tgl3`, arm `tglb`.

## 1. Verdict on option (b)

**Conditional GO as a release candidate; not approval of the literal integration plan.**
The sparse-target improvements and reported accuracy tradeoffs reproduce.
The important correction is that §6.1 does **not** establish flank flattening as the cause of the ladder loss.
**Blocker B1:** do not replace the tiled implementation with the unbatched §4.2 listing: one full matrix alone requires 4.37 GB; several coexist.
The frozen implementation already avoids B1. Preserve its batching while adopting the closed-form spike likelihood.
Correct the evidence and isolate the junction-cap change before integration; complete the planned integrated release gates before shipping.
The known statistical regressions merit explicit acceptance, but no specified accuracy threshold makes them automatic release blockers.
**Option (c): no-go for 0.8.0; reconsider after a corrected diagnosis.**

## 2. Findings

**1. Fix before integration — the principal flank diagnosis labels unprobed regions as capture flanks.**
Location: packet §6.1; `analysis/review_flank.py:18–30`; foundation `native/calibration/accumulator.cpp:737`; `sim/capture/sampler.py:12`.

*Measured independently:* **all 244** ladder g05 regions labelled “flank” have **zero overlap with any genomic probe part**. They comprise 191 intronic, 39 intergenic and 14 exonic regions; their median true gDNA count is **3**. The original −0.46/−1.00 median log errors reproduce. On test g50, **38/71** labelled flanks also lack probe overlap. Among the other 33, the medians are **−0.008/−0.105**, substantially smaller than the pooled −0.57/−0.74. Positive controls confirm probe overlap for **3,142/3,142** ladder and **363/363** test regions in the captured class. [Geometry audit](../../.cache/modefree_independent_review/flank_geometry.json)

*Derived:* a fragment wholly contained in an unprobed region cannot bind an external probe under this simulator. Proximity can enrich **crossing** fragments, which belong to boundaries. Moreover, selecting realized density `d ≥ 2ρ_b` forces the log error of an exactly bulk-valued prediction to be at most `−log 2`. The median is affected, not just the RMS. Thus this table cannot identify a systematic capture under-read.

**Test:** classify objects by expected capture yield from placement enumeration, independently of realized counts; score regions and boundaries separately, then `ruler_vs_truth.py --scale` and transcript allocations. The diagnosis survives only if genuinely enriched shoulder objects are under-read and repairing those objects improves the deciding strata. The current ladder classification already fails its geometric premise.

**2. Blocker B1 — preserve bounded working memory.**
Location: packet §9 step 1; `reference_reader.py:38–49` versus frozen `capture_efficiency.py:224–230`.

*Derived:* `2,100,000 × 260 × 8 = 4,368,000,000` bytes per dense array. The listing retains `lam`, `ll`, `a_alt` and `w_alt`, plus temporaries: over **17 GB** for those four arrays alone. The measured production path uses `_block_rows`; a literal transplant is a substantial memory regression.

**Test:** retain that existing allocation policy; compare several chunk sizes against the reference on saved dumps, then measure peak RSS at genome scale. Confirm unchanged efficiencies and bounded temporary memory. The reported 4.3-second timing does not establish either property of an unbatched integration.

**3. Fix before integration — “unit-free” does not mean behavior-preserving.**
Location: packet §§2.3, 4.3, 6.1; frozen `capture_eff_length.py:191–201`.

*Derived and reproduced:* with exon and boundary efficiencies 0.4 and intronic efficiencies 0.1, the old junction rule returns **0.7**, while the touched-object maximum returns **0.4**, in identical units. Scaling an old ceiling of 1 by a common factor is a unit conversion; replacing it with a local maximum is an additional model change. Equal calibration counts do not isolate the reader from this consumer change. Likewise, type-weighted linear `tgtla` versus `tglb` does not bound the log-readout contribution by 0.09 point.

**Test:** first run a cap-only arm on foundation efficiencies; then compare reader arms with identical consumers. For a readout-only diagnostic, hold the prior, spike, saturation and consumers fixed. Measure component-length scale and within-gene spread beside end-to-end error. Attribution is supported only if these controls explain the observed losses. This is a diagnostic, not a recommendation to ship a refused linear reader.

**4. Fix before integration — correct the claimed invariants and evidence scope.**
Location: packet §§3.2, 5.1, 8.

*Measured:* “every contaminated capture-OFF row matches” is false: lowg stranded-OFF is **6.97→7.51**, d10 **16.56→16.97**. Each gap panel contains **four g50 conditions**, not the ladder’s 16-condition grid; this matches the panel configurations. *Derived:* one-basin identity is valid, but capture-OFF need not imply one basin. The lone-level property bounds influence for a flat likelihood; it does not prevent influence, nor guarantee basin mass `1/n` under reliability weighting and smoothing. Log-density information scales with expected count `Sρ`, not `kS`: increasing support at fixed count does not make log-density precise.

**Test:** replace universal claims with their actual conditions; generate the coverage inventory and OFF results from JSONL. Test the one-basin algebra separately from simulated OFF behavior. Use actual basin mass in lone-object influence tests, including unequal reliability weights.

**5. Consider later — false gDNA is an upstream error amplified by this reader’s likelihood.**
Location: packet §§3.3, 6.2; foundation `simplex_logodds.py:310–314`.

*Measured:* the four highlighted g00 OFF losses are **+2.32 to +2.82 points**. *Derived:* published `f_g` is a **posterior median**, while `Var(log f_g)` is a posterior grid moment. A fractional `k=Mf_g` is therefore not an observed Poisson count; the likelihood is a working approximation. Gamma continuation accommodates fractional values algebraically but does not propagate deconvolution uncertainty. The solve already uses the landscape, so the reader’s posterior is not an independent second measurement. Its probabilities need empirical calibration before becoming footprint seeds.

**Test:** compare false micro-boundaries and genuine scarce-gDNA captured objects at matched `(k,S,f_g,Var(log f_g))`, using certified origins only for evaluation. If their summaries overlap, no deterministic reader of those summaries can perfectly distinguish them. Any proposed uncertainty treatment must improve all four g00 controls while retaining lowg g001 and d10 g001/g01 capture signal. Blanket variance widening has already failed that test; do not reinstate it.

**6. Consider later — saturation inherits discrete census failures.**
Location: packet §6.3; frozen `landscape.py:410–449`, `capture_efficiency.py:205–210`.

*Measured:* odg05 κ0.70 ON loses **3.11 points** in aggregate; its g25 condition is **15.47→23.84**. *Derived:* basin splitting changes the saturation targets, and a largest-basin swap changes which distribution is collapsed, not merely the unit. Endpoints are also omitted whenever an interior maximum exists. A synthetic density `[.70,.10,.05,.02,.03,.08,.02]` on `x=0,…,6` produces one basin at `x=5`, suppressing a real floor-versus-enriched contrast. This is a reproduced boundary case, not a measured cause of the benchmark losses.

**Test:** perturb inputs by one fragment and refine/shift the numerical grid diagnostically; track basin assignments, efficiencies and isoform ratios on odg05, sparse and junction panels. Exercise endpoint mixtures and tied basin masses explicitly. Do not merge basins or introduce a smoothing threshold based on benchmark improvement.

**7. Fix before integration of option (c) — its seed probability is wrong.**
Location: packet §8 step 1; `reference_reader.py:48–50`.

*Derived:* `w_null` is unnormalized. Even `1−w_null/den` includes below-bulk alternatives. Use

`P_o⁺ = Σ[j in above, outside bulk] w_alt[o,j] / den[o]`.

*Measured independently:* on ladder g05 boundaries, raw `1−w_null` differs from the normalized non-null probability by as much as **0.936**; below-bulk states create a further distinction. No footprint was implemented, so this is not a defect in option (b)’s efficiencies.

**Test:** hand-computed spike/upper/lower mixtures must reproduce all three posterior masses. Include a strongly supported depleted object: its enriched probability must approach zero, not one. Test predictive calibration before using these values as neighbour evidence.

## 3. Improvements to option (b), ranked

All proposed numerical changes require **one mechanism per pinned A/B**: `OMP_NUM_THREADS=1`, `--set scan.total_threads=1`, `--set em.assignment_mode=fractional`, seed 0; ladder, both gap directions, every robustness panel and all five real inputs. Read each stratum, gene errors and pools separately; report the deferred stratum. Diagnostic constants are not new production tunables.

| Priority | Improvement and expected benefit | Risk, effort and deciding test |
|---|---|---|
| 1 | Preserve batching, finite-value validation and correct units; avoids operational failure without changing inference. | Low; dump identity, memory measurement and consumer scaling tests. |
| 2 | Correct the truth classification and separate cap/readout effects. Highest value for avoiding a misdirected statistical change. | Low–moderate; Findings 1 and 3. |
| 3 | Investigate uncertainty-aware reading of micro-boundaries; the measured g00 losses exceed the ladder’s aggregate gap. | High statistical risk; Finding 5’s matched-summary and rare-gDNA tests must precede a candidate. |
| 4 | Revisit the bulk shoulder only if expected-yield diagnostics establish a real missing level. | Moderate–high; geometry-defined classes and within-gene error, then the full A/B. |
| 5 | Explore saturation stability without a new global cap. Keep the measured per-level rule for this candidate. | Moderate–high; Finding 6. Neither “no saturation” nor one common level survived the existing matrix. |
| 6 | Retain the log readout and peak bulk unit pending an isolated benefit. | Low benefit established so far; pooled bulk was flat, and the 0.09-point attribution is confounded. |

The collapsed null is a defensible **regularizer**, not an identified distribution of all unprobed objects. Efron’s two-groups framework addresses null specification; it does not supply Rigel’s largest-basin assumption. That assumption remains this model’s responsibility. [Efron, *Microarrays, Empirical Bayes and the Two-Groups Model*](https://arxiv.org/abs/0808.0572)

There is nevertheless a genuine structural limitation: where the null removes an intermediate true density, more counts cannot restore absent prior support. For a two-point approximation with levels `1,R`, the large-support decision changes at `r=(R−1)/log R`; lower intermediate levels converge to 1. This motivates testing real shoulders, but does not validate the packet’s selected “flanks.”

Linear use of efficiencies does not make their arithmetic posterior mean automatically optimal: EM subsequently divides by lengths and normalizes competing components. For saturated log efficiency `Z`, `log E(exp Z)−E(Z)≥0`, but that identity alone gives no transcript-error attribution. Any exact arithmetic mean retains the rare-outlier contribution `μ(R−1)`; obtaining logarithmic robustness requires changing the loss, the posterior evidence or the model. No safe, otherwise identical linear replacement is established here.

## 4. Option (c): the probe-footprint mechanism

**The proposed causal claim is presently unsupported.** Probe effects are spatially local at the fragment level, but the packet confuses a fragment’s footprint with an object’s genomic neighbourhood. An adjacent long intron’s contained fragments cannot reach an exon probe. Conversely, one partly probed long region can contain a capture gradient without any enriched neighbouring object. Annotation reaches locate template opportunities; they do not locate probes inside an object.

If corrected diagnostics justify a prototype **after 0.8.0**, derive locality from placement geometry. Let `d_o(u,l)` be the deposit indicator for the reader’s object, `f_g(l)` the existing gDNA length law, and `A` the latent set of probe-active sequence parts. A defensible target for the locality weight is

`q_o = [Σ_l f_g(l) Σ_u d_o(u,l) Pr(fragment(u,l) overlaps A | neighbour evidence)] / S_o`.

This averages a **union event**, preventing duplicate credit from several neighbours. It gives zero to contained fragments that cannot reach any active part. It uses the whole length law, including its tails, not “one fragment length” chosen as a cutoff. Use the RNA law separately when validating transcript/junction consequences; there is no justified average of the two laws.

The difficult ingredient is the posterior over `A`. Corrected `P_n⁺` gives probability of an elevated **object-average density**, not probe position, affinity or overlap geometry. Known probe support or an explicitly justified latent-position model is necessary. `max P_n`, an independence product, or uniform probe position within a region introduces an assumption; none follows from fragment length alone. Shared crossing fragments and calibration messages also correlate neighbours. Validate using disjoint fragment evidence or paired simulations before allowing evidence to travel out and back to its source.

A conservative, mass-preserving **candidate** shoulder allocation is derivable without a tunable amplitude. Set `H={x_b<x≤bulk.hi}` and `h=Σ_H p(x)`, then

`G_o = (π₀−q_o h)δ_xb + q_o p(x)1_H + p(x)1_outside_bulk`.

It restores at most the existing shoulder mass, rather than reallocating the entire spike to a normalized shoulder. It remains a modelling choice: the rendered shoulder contains measurement noise, and this formula does not identify which part is physical capture. If `h=0` or `q_o=0`, it is exactly option (b).

**Gates:** physical placement overlap and continuous evidence only—no capture detector, probability threshold, gene pooling or library switch. Test isolated probes, long unprobed introns, mixed-length laws, junction probes, adjacent false micro-boundaries and annotation subdivisions. First require better expected-yield and within-gene errors; then run the pinned full matrix. Permuting neighbour locations while preserving their marginal evidence must remove a truly spatial advantage. Reject a gain driven by fitting the selected realized-count “flanks.”

Restricted shoulder opening differs from the refused readers that coupled every object’s level, so their exact failure need not recur. Their risks—correlated evidence, noise propagation and junction interactions—remain. Development cost is substantial, gain currently unquantified. The **0.47/1.23/1.12** foundation gaps are neither measured recoverable gains nor mathematical ceilings. **No-go as sketched; a small diagnostic investigation after release is justified, implementation is not yet justified.**

## 5. Integration risks and test strategy

**Units and surfaces.** Efficiencies are relative yields, not probabilities; `em_effective_length` may exceed the plain length. TPM’s length divisor remains plain, but TPM values can still change through revised counts. Rename metadata consistently; training count must not masquerade as captured membership. Update `report/assets/report.js:395` as well as `report/model.py`: the JavaScript currently says an enriched mode was located whenever reference density exists. That would misdescribe ordinary one-basin OFF libraries.

**Consumers.** Remove the locus span ceiling and use a homogeneous junction rule, while preserving Finding 3’s distinction between invariance and changed physics. Scale every participating region and boundary efficiency by a common positive factor: all transcript, synthetic-RNA and gDNA lengths must scale together; fractional allocations and count pseudocounts must remain invariant. Search source, native code, serialization, reports and instruments for old bounds and reference-based bypasses.

**Edges.** No landscape returns ones; one basin returns ones even with counts inconsistent with its peak. Empty object arrays must work. Below-bulk alternatives legitimately give efficiencies below 1. Reject non-finite inputs and all-zero/non-normalizable priors; skip zero-mass basins without NaNs. Preserve and test valley-bin ownership: the census masses use half-open intervals while the reader’s masks are inclusive. Exercise endpoint and tied-peak cases rather than silently altering the census during integration.

For `S_g=0`, returning 1 is a fallback, not a measured efficiency. With unequal length laws, `S_RNA` can be positive; such a piece can still affect a transcript or the junction cap. Test a region between the supports of the two fragment-length laws. Apply unit changes to that fallback too; otherwise “only ratios matter” does not protect the result.

Keep Python. Rewrite the named tests around algebraic properties, edge cases and meaningful consumer outcomes; verify each new gate by mutation. Inspect truth-scored golden differences before regeneration. Finally require the integrated full matrix, real-library checks, complete suite and `preflight.py --full`; port the multimapper pass separately and retest the combined candidate because it changes the reader’s inputs.

## 6. Evidence check

**Verified here:** all **632** manifest hashes and the ZIP CRC/SHA; all six frozen source files; the 16 mode-free result receipts and seed/assignment records; the **14 better / 83 within ±0.3 / 21 worse** tally; switch-off identity on **108 rows** excluding runtime; and reference-reader equality on seven saved dumps, **91,284 object evaluations**, maximum log difference **4.53×10⁻¹⁴**. The packet’s 70,176 is the ladder dump’s size, not their combined size. Its validator also reproduces all 33 quoted values it checks. [Audit results](../../.cache/modefree_independent_review/audit.json)

The §1 numerical gains and costs are supported. Separate the tally by release scope: **96 in-scope strata: 10 better, 70 within ±0.3, 16 worse**; **22 deferred: 4/13/5**. These are descriptive counts of correlated conditions, not significance tests or 118 independent replications.

Ratio-of-sums weights conditions by true RNA count. Ladder stranded-ON g98 contributes only **1.36%** of the stratum denominator; its error is **28.54→29.86**, beside the stratum’s **3.87→4.34**. Retain per-condition and probed/unprobed results; odg05 g25’s **+8.37 points** also exceeds its aggregated loss materially.

Realized counts are appropriate for the stated transcript-count target. They are insufficient to label latent capture efficiency, particularly after selecting on those same counts (Finding 1). The simulator’s best-contiguous-probe-overlap rule supports geometry checks, but agreement under it does not establish robustness to sequence-dependent affinity, mapping bias or nonuniform genomic templates. Use independent simulations for generality; real libraries remain validation inputs, never tuning targets.

The real-library enrichment column uses incompatible raw units for the two readers; the packet acknowledges this. Rescale to a common bulk unit before comparing counts of enriched objects. The two plasma gate failures are demonstrated behavior, not demonstrated accuracy gains. VCaP recovery **91.9% versus 91.8%** gives no superiority claim, and raw primary-pair truth is not the pipeline-filtered denominator. `main` changes other components and is not a reader-only control. Frozen builds and switch-off identity are strong provenance checks; they do not replace the isolated consumer A/B or the final combined release run.

I reran the audits and targeted geometry/numerical checks, **not** the full BAM-to-quantification matrix or production test suite.

## 7. Answers to packet §10

**1. Null.** The spike is useful noise suppression, but largest rendered mass does not identify the biological null or separate flanks from noise. Keep it for this candidate. Diagnose expected yields before reopening support; the packet’s ladder flank evidence is invalid.

**2. Readout.** Retain the geometric readout. Linear exposure aggregation does not specify the optimal decision loss for downstream normalized allocations. There is no demonstrated safe arithmetic replacement, and the claimed 0.1-point contribution needs a readout-only experiment.

**3. Saturation.** A partition-free smooth shrinkage rule is possible, but a physical saturation density is not identified by these summaries. Probe-overlap physics could supply a bound; an arbitrary smooth cap would add assumptions or constants. Test census stability first; do not restore the refused common cap or no-cap variants.

**4. Uncertainty.** A median and one log-variance do not determine tail direction or the probability of genuine enrichment. A bounded, asymmetric posterior could help if available, but cannot be uniquely reconstructed from these outputs. Reusing it would also require accounting for the prior already applied. Do not redesign the solve or promise a moments-only cure for 0.8.0.

**5. Remaining losses.** Missing prior support and scarce-count indistinguishability are structural given the current inputs. Census sensitivity and consumer pricing are potentially reader/consumer repairs. Micro-boundary amplification may be reduced, but cannot be guaranteed removable while preserving real capture at identical summaries. The deciding-panel loss still needs causal attribution.

**6. Footprint.** No-go before release. Derive placement overlap, correct the enriched posterior and conserve shoulder mass as in §4 only after establishing a real spatial residual. Region adjacency alone is insufficient; the existing “ceiling” does not price the mechanism.

**7. Integration.** Preserve tiling, test homogeneous consumers including zero-gDNA-support/RNA-supported pieces, update the JavaScript report and metadata, validate priors and finite arrays, and exercise basin boundaries. Keep statistical changes separate from the port and from the multimapper pass.

**8. Evidence.** The numerical results support a sparse-target benefit with known costs. They do not support the flank mechanism, its promised recovery, universal OFF identity or plasma accuracy improvement. Correct those interpretations while preserving the measured tables and all per-condition losses.

Reproduction: [audit.py](../../.cache/modefree_independent_review/audit.py) and [check_flank_geometry.py](../../.cache/modefree_independent_review/check_flank_geometry.py), run with the `rigel` environment’s Python and `OMP_NUM_THREADS=1`. They read the frozen bundle/live artifacts and write only review diagnostics. Production source and the supplied packet were not changed.
