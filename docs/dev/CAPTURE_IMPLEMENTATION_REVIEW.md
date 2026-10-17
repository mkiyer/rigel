# Capture implementation review — 2026-10-11 (sandbox)

The owner clarified that this session supplies the independent critique and pursues a
shippable 0.8.0. No external text is missing. Main baseline: `0fdb9dce`. No commits.

## Review outcome

**No inference candidate is ready for 0.8.0.** The best ladder improvement fails
RNA-short stranded capture-ON, 14.242734→20.434400% transcript error. The two
coupled readers also fail both captured gap strata. Main inference remains unchanged;
the corrected ruler instrument and this review are the only new working-tree changes.
No commits. The deferred stratum remains deferred; its failure is not substituted
for the three release targets.

### Verdicts on the five review points

1. **The weight path owns the regression, but “unbiased shallow noise” is not a
   sufficient diagnosis.** Among 70,176 ladder objects, 696 log weights exceed 1,
   352 exceed 5, and 93 exceed 10 on the flat stranded field. Their depth reaches
   3,769 (median 6). Received composition/DNA factors cause many mode switches;
   removing RNA-level factors changes none of the 101 sampled modes by more than
   one log unit. In that sample, mode versus posterior-SD rank correlation is
   .694592 and mode versus depth is .506565; this does not calibrate a shrinkage
   scale. Those 101 were deliberately selected, not a population sample.
   Own-only g00 transcript error falls 13.609200→2.789140%, silent mass
   251215.742→32061.638; coupling own likelihoods reaches 1.992229% and 16879.440,
   but loses captured truth. The top ten baseline silent isoforms carry 138426.279
   fragments; only 3354.086 are on isoforms one exact exon-set element from an
   expressed sibling. A one-shallow-exon explanation is too narrow.
2. **The asymmetry is already present upstream of EM.** Maximum log weight is
   13.102353 stranded versus .0000476 unstranded on the ladder. Only 49/696 elevated
   objects are both-strand, none of the 93 most extreme. The test chromosome's
   3,518 genic slots are essentially flat and cannot decide the tail. Reader-only
   A/Bs preserve the six count arrays yet move isoform errors strongly. This
   establishes an upstream cause, not equal EM sensitivity under an unrun field
   transplant.
3. **Neither failed marginal average establishes that every mode is faithful;
   the tested local-prior repair is refused.** Both coupled versions fail the
   captured gap tests; removing all delivered factors without coupling still fails
   RNA-short by +6.191666 percentage points. The same two marginal length laws
   give expected responsibilities .800000, .720166 or .640333 under different
   correlations. Carrying shared uncertainty into allocation remains a distinct,
   untested mechanism, derived below. Deep-object bit identity is not promised by
   adding a proper prior; the derived finite-curvature displacement bound is d/k.
4. **A local exon/intron contrast is allowed, but is not an identified capture
   witness by itself.** Enriched DNA (700,1) and flat DNA (1,1) plus exon RNA 699
   give identical total/strand/equal-length likelihoods: max difference 0. The
   claimed direct length-composition channel is unbuilt: zeroing both inverse-length
   banks changes every calibration array by 0. Primary deferred RNA-long gDNA is
   1.865093 M against 5 M truth; its walled-exon signed gDNA error is −1.161884 M
   object incidences. Incidences are not additive library counts. Matched length
   evidence may identify composition when the laws differ; bulk means do not
   establish that evidence. Reopening a one-sided bound as a point observation is
   unsupported. No new deferred-stratum estimator was selected.
5. **The truth instrument contained a real bug; raw scale claims need a unit
   audit.** The ruler used RNA's generating PMF for gDNA on both gap panels. It now
   reads `manifest.gdna` with the simulator's own defaults. Equal-law ladder output
   remains exactly unchanged. There are 36 passing self-test gates, both loader
   mutations fail, the PMF-consumer mutation fails, 17 repository boundary/index
   tests pass, and lint is clean. Corrected gDNA/isoform raw scale ratios 2.248805
   and .510912 do not alone prove an allocation defect: a finite model with raw
   ratio 1.5 has exactly correct K/L probabilities (max difference 0). Changing
   L alone to align those scales breaks normalization.

### Complete deciding A/B

Transcript absolute error / true fragments, percent. “Coupled delivered” adds
only the local field to the shipped reader. “Coupled own” removes delivered
reader factors from that field. “Own only” removes the coupling as the final
factorial ablation. Every pair is pinned and fractional; all count inference is
held fixed. No empirical coefficient was chosen or tuned.

| Panel / stratum | Shipped | Coupled delivered | Coupled own | Own only |
|---|---:|---:|---:|---:|
| Ladder stranded OFF | 7.515343 | 3.697799 | 2.159653 | 2.598793 |
| Ladder stranded ON | 9.711430 | 10.806300 | 9.071331 | 6.465506 |
| Ladder unstranded OFF | 2.068994 | 2.071055 | 2.034908 | 2.037749 |
| Ladder deferred ON | 40.415036 | 48.509202 | 35.352037 | 34.347237 |
| RNA-long stranded OFF | 3.517133 | 2.847446 | 2.473139 | 3.465777 |
| RNA-long stranded ON | 6.093690 | 11.508120 | 11.468872 | 6.541159 |
| RNA-long unstranded OFF | 2.547063 | 2.506615 | 2.250067 | 2.212726 |
| RNA-long deferred ON | 34.107782 | 28.740063 | 14.660790 | 14.752888 |
| RNA-short stranded OFF | 4.990618 | 4.858462 | 4.297629 | 4.558286 |
| RNA-short stranded ON | 14.242734 | 17.462804 | 18.360167 | 20.434400 |
| RNA-short unstranded OFF | 3.951550 | 3.869995 | 3.891142 | 3.937907 |
| RNA-short deferred ON | 68.055249 | 65.837002 | 22.776427 | 63.708295 |

The ladder's captured contaminated-only transcript error is 6.047457% shipped,
11.752182% coupled delivered, 10.836850% coupled own and 6.365602% own only.
Thus even the own-only aggregate improvement includes a tradeoff hidden by g00.
For own only, captured g05/g50/g98 read 5.694003→6.039520,
5.944245→6.225588 and 25.413098→25.351186%. The .70-specificity test control
worsens 8.042554→10.007251%; test g98 stranded ON worsens
52.686796→72.213579%. Single-probe stranded ON improves
11.462933→8.670760%, but cannot override those failures.

At ladder g50 stranded ON, within-gene log(L/Y) SD is .049226 shipped,
.344738 coupled delivered, .235997 coupled own and .063994 own only.
RNA-short own-only silent mass rises 128607.837→186888.081; genes rise
1.498370→1.740974%. The deliverable, not the smoother-looking field, decides.

VCaP: the repeated baseline transcript files and pools are byte-identical.
Baseline gDNA fraction .252753061 becomes .241557054 with coupled delivered
and .238243947 with coupled own. The approximate aligned-primary origin fraction
is .2518; it is not the filtered pipeline denominator and provides no transcript
truth. Own only gives .253232601, nascent RNA 517025.204 versus 523324.296
shipped; it does not improve agreement with that approximate reference. All five
new VCaP executions produced complete outputs; reported peak process RSS stayed
below 11.1 GB. The final own-only supervisor failed on a sandbox-denied `ps`, so
its child exit status was unobserved and the live RSS guard was inactive. The
254,319-row transcript file has unique IDs, finite nonnegative counts, the exact
reported RNA sum, conserved filtered total and a final log matching the JSON; an
authorized final process inspection found no workers. Its receipt records this
recovery explicitly. The final candidate-only campaign consumed about 27.5 minutes
against its announced 25–30-minute budget (accuracy, not a timing comparison).

Primary calibration was independently scored on all 24 ladder/gap conditions.
Ladder region/boundary absolute errors as percent of observed mass are
.591179/.632989 stranded OFF, .757625/.982237 stranded ON,
1.409653/1.014926 unstranded OFF, and 6.600975/10.952600 deferred ON.
All reader variants leave this primary result unchanged. At g00 stranded OFF,
primary calibration reports 1593.742 gDNA fragments; 1126.426 is the baseline
end-to-end pool, not that calibration number. The zero control is retained.

The full evidence lives in
[the review archive](../../.cache/rigel_runs/2026-10-11_implementation_review/README.md):
270 condition comparisons across six panels, all rung/gene/pool readouts,
239 saved paired six-array equality checks, and 31 earlier conditions whose array
hook was absent (explicitly missing, with separate cached target checks).
Standing suite count is 2,837; it was not rerun or presented as validation of
these isolated candidates. Elapsed times are budget receipts under overlapping
accuracy work, not speed measurements. No performance claim is made.

### Refusals

- Do not ship any of these three reader candidates, the two previously refuted
  independent averages, or flat weights justified by the attribution ablation.
  Captured and gap truth must survive the repair.
- No depth cut, clipped fallback, library capture detector/reference density,
  empirical coefficient search or tuning to a favorable panel.
- Do not call an exon/intron total contrast a DNA-composition observation, or
  restore the removed intron constraint merely because an older arm scored better.
- Do not repair a raw cross-population scale ratio by changing normalizers alone.
- Do not accept a pooled/g00-dominated score, a deferred-stratum gain, passing unit
  gates, or approximate VCaP pool truth as a release decision.

The next mechanism and its mathematical, evidence-ownership and full-panel
falsifiers are derived at the end of this note. It is not a release claim.
Earlier sections below preserve the order of derivation and measurement.

## Claims to test

1. **Attribution is established; the noise model is not.** Forcing weights to one
   identifies their downstream path, but does not establish independent Poisson
   noise, unbiased modes, or a relationship between posterior width and error.
   At g00 the true gDNA density is zero: a relative gDNA-density truth is undefined.
   The capture-OFF truth is instead a constant capture field and the simulator's
   transcript yields. Measure mode differences by object class, opportunity, RNA
   depth, strand ambiguity, and delivered factors; compare within-gene length ratios
   and false-positive transcript mass. Test chromosome first, ladder to decide.
2. **The stranded asymmetry may be upstream of the EM.** Compare both g00 OFF
   conditions on identical annotation. Different fields implicate the reader before
   EM; the flat-weight attribution and count-preserving reader A/B then test the
   downstream error path. This does not separately identify the EM's relative
   sensitivity to a common transplanted field.
3. **The two failed marginal averages do not validate the mode or local shrinkage.**
   The marginal expected length is not generally the likelihood marginalized over
   a shared uncertain field. Any replacement must state its loss/model and retain
   correlations of shared objects. Inspect the actual likelihood and integration
   before deriving a candidate; avoid a new whole uncertainty framework without a
   measured need.
4. **The deferred failure needs an observed information-path audit.** Quantify
   signed gDNA/RNA error by class, composition evidence and fitted length laws on
   the gap panels. Equal configured lengths are a control, not a proof that all
   information is absent after capture. A local exon/intron contrast can also encode
   real RNA abundance or splicing, so it is not automatically a capture witness.
5. **Release validation must include zero controls and recoverable artifacts.**
   Archived ladder/gap A/B tables reproduce their reported figures. The prior
   campaign's test-chromosome logs survive but its driver copied the wrong result
   path; the per-arm test JSONL files are absent from those archives. Fresh arms
   must write to isolated explicit paths, record provenance and include every rung.

## Execution discipline

No detector, library reference, new empirical constant, or panel tuning. Derive a
mechanism here before prototyping outside main. Freeze a current baseline in the
same build/session; pin scan threads to one, use fractional EM, and leave calibration
and EM at their normal thread settings. Screen cached observations before end-to-end
runs. Every candidate gets all three in-scope strata and both gap panels, with the
deferred rows separately reported, before selection. VCaP is a truth-bearing test,
not a fitting input. No production source change before that evidence.

Initial budget announced: test g00 OFF pair ~1 minute; corresponding ladder cached
diagnostics ~2–3 minutes, no EM. System storage services currently consume about two
cores; these are accuracy measurements, not performance timings. Do not terminate
unrelated processes or call their loaded runtimes speed measurements.

## First measurements

Cached g00 OFF calibration on main: test chromosome 3,518 slots, ladder 70,176
slots. On the ladder, the stranded row has 696 objects with log weight above 1,
352 above 5, and 93 above 10; the maximum is 13.102353. The unstranded maximum
is 0.000048. Those are descriptive bins, not proposed thresholds. Only 49 of the
696 are both-strand objects, and none of the 93 most extreme are both-strand.
Their depths span zero to 3,769 reads (median 6): a depth cutoff is unjustified.
The test chromosome's ordinary genic objects stay essentially flat in both rows;
it cannot decide this particular defect. Its designed intergenic shadow is separate.

An adapter reconstructing the native reader's delivered factors reproduces every
sampled mode exactly. Removing composition or DNA-level factors separately removes
many extreme jumps; removing RNA-level factors does not. Own strand evidence also
creates a few jumps. For example slot 35579, counts (800,21), has log weight 13.102
and posterior SD 4.899; removing every delivered factor still gives 13.126. Slot
30200, counts (207,7), has log weight 12.708 and SD 4.368; removing only its DNA
factor returns it to the ordinary level. Thus the simple independent-noise account
is incomplete: one-sided messages and multimodal posterior decisions matter.

Receipts and executable diagnostics: `<scratch>/rigel-review-20261011/`.

## Candidate derivation, before prototype: a local capture-field prior

Keep the entire shipped per-object log posterior `ell_i(x)` and all count inference
fixed. Couple adjacent owned objects within each existing locus block; no link
crosses blocks. This first experiment uses genomic adjacency only, not an invented
distance window or inferred probe position.

The approved local enrichment prior has `p(C)=C^-2`, `C>=1`. Under `d=log(C)`,
`p(d)=exp(-d)`, `d>=0`. Neither adjacent object is a privileged reference. Giving
the two directions equal odds yields the normalized symmetric difference prior

    p(x_i-x_j) = exp(-abs(x_i-x_j))/2.

The factor 1/2 normalizes the two directions, and the unit exponent follows from
the approved prior's change of variables; neither is a new fitted constant. The
field model, conditional on the same evidence summaries that ship, is

    p(x | summaries) proportional to
        product_i exp(ell_i(x_i)) * product_adjacent exp(-abs(x_i-x_j))/2.

This is an explicit approximate field model: adjacent deposit/message summaries
are correlated, so multiplying them is not a claim of independent raw fragments.
That limitation is a falsifier, especially in dense chains and overlapping genes.
The experiment changes only the local prior, not the summaries, the count solve,
the median unit, the junction rule, or the EM. It neither reads a library reference
nor decides whether capture exists.

Compute exact marginal posteriors of this chain on the existing density lattice
by sum-product; publish each marginal mode using the existing parabolic refinement.
For a neighbour message and uniform grid step h,

    left[j] = logaddexp(a[j], left[j-1]-h)
    right[j] = logaddexp(a[j], right[j+1]-h)
    message[j] = logaddexp(left[j], right[j+1]-h)

with absent terms omitted. This sums `exp(a[k]-abs(x[j]-x[k]))` without counting
the centre twice. The two spatial passes cost O(objects * density-grid-points),
with memory bounded by the existing locus block. A common additive shift in every
log density leaves the prior and all published ratios invariant.

### Falsifiers and gates

- The chain marginals must agree with brute-force enumeration of a tiny joint
  posterior. Deleting a coupling edge must make that gate fail.
- One object must return the shipped full-grid mode; translating the grid must
  translate every mode by the same amount. Baseline switch OFF must reproduce main.
- A sharply measured isolated enriched object must retain enrichment with flat
  neighbours. Deep-object *bit identity* is not mathematically promised by adding
  a proper prior: movement vanishes with uncertainty, rather than being exactly
  zero at finite depth. Claiming both is an invalid falsifier.
- On flat truth, reduce within-gene length-ratio error and silent-transcript mass;
  do not merely shrink the plotted object spread. Fresh pinned/fractional A/B on
  every ladder rung and both gap panels decides. Report all deferred rows.
- Single-probe and junction-probe test layouts, zero controls, and VCaP truth must
  expose erased real enrichment, amplification of correlated messages, or a
  worsened gDNA/RNA split. A failure refuses this candidate, not the data.

## Channel audit and prototype gates

The brief's statement that calibration already uses fragment-length composition
evidence is incorrect. The direct channel is unbuilt. No calibration source reads
`boundary_unspliced_inv_length_sum` or `sj_inv_length_sum`. At test g50 unstranded ON,
setting both banks to zero removes inverse sums 1276.991346 and 432.600086 yet changes
every result array by exactly zero. Learned fragment laws still enter geometry, and
the EM scores individual fragment lengths; these are distinct information paths.
`length_channel.json` records the perturbation.

The four local-field gates fail before implementation, then pass; 11 existing reader
gates also pass. Removing coupling fires the joint-enumeration and propagation gates.
On the current ladder target the prototype switch OFF reproduces main's cached
counts and weights exactly. The switch ON keeps the count arrays identical and
reduces within-gene log-length SD over annotated multi-exon isoforms from 0.968027
to 0.237380; objects with log weight >5 fall from 352 to 54. This is a diagnostic,
not an end-to-end victory.

Fresh test-chromosome fractional/pinned A/B (five rungs per stratum): transcript
error stranded OFF 8.82235→8.49405%, stranded ON 8.18443→7.66183%, unstranded OFF
8.07834→8.17289%, deferred ON 12.33493→12.77846%. Gene error stranded ON
0.82682→1.06024%; deferred ON 3.28562→5.57649%. The .70-specificity stress rows
are separate: OFF 8.44313→8.27820%, ON 8.04255→9.36993%. The candidate is mixed,
not selected. Accuracy-run elapsed times 35.8/39.8 seconds are budget receipts,
not idle-machine speed measurements.

Fresh ladder g00 stranded OFF A/B: transcripts 13.609200→4.844192%, genes
0.832732→0.193961%, silent-isoform mass 251215.742→63194.910; estimated gDNA
1126.426→1148.866 fragments, true zero. Count inference is unchanged; the small
end-to-end gDNA movement is downstream allocation. The flat-weight attribution
ceiling remains 1.93%, so this candidate recovers much, not all, of the regression.
The pair took 192.1/193.7 seconds under background OS load; no speed claim.

The equal-length control does not mean unmatched bulk populations have identical
observed length distributions. Origin-split oracle histograms at ladder g50 unstranded give
mean observed lengths OFF: gDNA 216.697, annotated RNA 208.256, synthetic-span RNA
216.418 bp; ON: 240.772, 225.483, 238.806 bp. Geometry and capture condition the
shared generating law differently. These bulk means do **not** establish a local
composition witness: geometry and capture cancel only inside a matched footprint.
The owner's earlier matched exon control was 219.00 versus 218.98 bp; equal laws
carry no composition information there. RNA-long captured bulk means are 81.693/250.550/257.789 bp, and RNA-short
258.776/81.538/81.712 bp, respectively. Receipt: origin_length_means.json.

A useful finite-depth bound follows from the chosen kernel. A neighbour message is
`M(x)=sum_y a(y) exp(-abs(x-y))`, so its log derivative is a weighted average of
`-sign(x-y)` and has magnitude at most one. With d adjacent factors, the added
log-posterior slope is bounded by d. If the object's original log posterior is
strongly concave with curvature at least k between its old and new modes, then
`abs(new_mode-old_mode) <= d/k` (up to lattice/refinement error). Thus strongly
measured objects move vanishingly little without a depth gate, while broad or
multimodal objects are allowed to move. The condition on curvature matters; it
is not a universal guarantee for a multimodal posterior.

## Uncertainty in effective lengths: what the two failed averages do not test

For abundance masses `a_t`, conditional EM responsibility has the form
`r_t(W) = (a_t k_ft / L_t(W)) / sum_s(a_s k_fs / L_s(W))` when the fragment's
capture multiplier is common to compatible components. Integrating this ratio is
not the ratio evaluated at an integrated weight or length. Furthermore, the
components share objects, so their length errors are correlated. Independent
posterior averages discard this information. A fully joint likelihood would
integrate the shared field once across fragments, not re-use its posterior as
independent evidence for every fragment. The current prototype only tests a local
field prior upstream of the same plug-in EM; it is not that joint likelihood.

A finite enumeration (0.04 s, no data fit) makes the limitation exact. Take abundance
masses (0.8, 0.2), equal fragment compatibilities, and each component's length
marginal equally likely 1 or 9. Both marginal means are 5 in all three models.
Expected responsibility for the first component is 0.800000 with common lengths,
0.720166 with independent lengths, and 0.640333 with opposed lengths. Identical
marginal uncertainty does not identify the allocation. Thus neither an independent
mode nor an independent mean can promise to carry uncertainty into isoform
allocation while preserving every location. `uncertainty_identity.py/json` is the
executable counterexample, not a proposed production estimator.

## Probe-layout robustness and truth-domain audit

The single-125-bp-probe test (two rungs per strand level) gives stranded ON
transcript error 11.462933→10.623224%, genes 1.078929→1.072989%; deferred ON
transcripts 9.178974→9.164645%, genes 1.743071→1.754455%. The local signal survives.
The junction-only panel is a counterexample to uniform robustness: stranded ON
transcripts 28.207586→29.124832%, genes 1.053550→1.229744%; deferred ON transcripts
44.458006→44.203812%, genes 4.426931→5.126374%. Its deferred g98 ON row worsens
188.955293→247.917963%. OFF rows reproduce the test chromosome because those BAMs
share the OFF construction. These are not additional independent OFF replicates.
The single-probe pair consumed 25.3/27.4 s; junction 263.9/125.2 s, under changing
load. The latter's initial three-minute joint budget was revised to nine minutes
when its baseline finished. No speed conclusion is drawn.

VCaP's cited 0.2518 is the aligned primary first-in-pair read-name fraction
(4,708,907 DNA; 13,994,657 RNA), not exactly the production filtered population.
The earlier memory `capture_design_revised_2026_10_05.md` already documents the
18.70 M versus ~18.5 M domain mismatch. It is an approximate origin-pool check;
no transcript-level truth exists here. A movement closer to 0.2518 cannot by itself
establish an exact matched-population accuracy gain.

## Instrument correction, derived before changing it

`ruler_vs_truth.load_sim` constructs only the RNA generating length PMF from
`manifest.simulation`; `gdna_space_yield` then uses that same PMF. The gap manifests
have a separate `manifest.gdna` law (RNA-long: gDNA mean 75, RNA mean 250; RNA-short:
gDNA mean 250, RNA mean 75). Thus the current gDNA-space truth on these panels is
wrong. RNA truth and the equal-law ladder are unaffected.

The simulator yield for a gDNA component is
`Y_g = sum_w p_g(w) sum_{s: fragment(s,w) overlaps footprint} capture_g(s,w)`;
RNA yields use `p_r(w)` in their own transcript coordinates. The correction is to
carry both generating laws and use the gDNA law only in `gdna_space_yield`, with
per-field fallback to `GDNASimConfig` defaults, matching the simulator's configuration (not the RNA parameters).
No measured fragment histogram replaces a generating PMF. Gates: an uncaptured
interior footprint has yield `length + E_g[w] - 1`; captured yields must agree with
direct start/width enumeration on a small reference; changing RNA widths alone
must not change gDNA yield, and changing gDNA widths must change it. Equal-law
results must remain bit-identical. Correct outside main first, then re-score the
saved A/B lengths against the corrected truth; this changes the instrument, not
the inference mechanism.

## Follow-up derivation: local field with each likelihood admitted once

The first candidate fails in scope: ladder stranded OFF 7.515343→3.697799%, but
stranded ON 9.711430→10.806300% (genes 1.874753→2.569200%); RNA-long stranded ON
6.093690→11.508120%. At ladder g50 stranded ON the within-gene log(L/Y) SD rises
0.049226→0.344738. This is not a release candidate.

One identifiable defect of that model is reuse of evidence. Its node potential
`ell_i` already includes messages `H_i(D_neighbours)` and a landscape fitted to
the observed sample. The new local-field message then sends those same neighbours'
data back along the chain. Exact sum-product on these composite node potentials
does not make the underlying data likelihood exact. The g00 factor ablations
showed that delivered composition/DNA factors create many of the rare jumps.

The next single experiment is a **reader with local messages that exclude the
recipient's incoming message**, using only own observational likelihoods as node
potentials. Keep the fitted landscape and count solve fixed. Let

    phi_i(x_i) = P(x_i) * integral Pois(u_i; exp(x_i) Eg_i/2 + q_i r)
                                Pois(v_i; exp(x_i) Eg_i/2 + (1-q_i) r)
                                r^(-1/2) dr.

Use the existing both-strand share integral when both RNA strands are admissible,
and no RNA integral at a structurally RNA-inadmissible object. Do not multiply
received composition, DNA-level or RNA-level factors into this reader potential.
The already-derived local enrichment prior stays `exp(-abs(x_i-x_j))/2`. Its
message is then

    m_(i→j)(x_j) = sum_(x_i) phi_i(x_i) exp(-abs(x_i-x_j))/2
                            product_(k in neighbours(i) except j) m_(k→i)(x_i).

Each own likelihood occurs once on a tree. The count foundation's messages are
unchanged; this is a reader experiment. It still has a limitation: adjacent object
counts can contain the same fragment, and the fitted landscape is empirical Bayes.
Removing delivered-message reuse does not resolve those remaining correlations.

Falsifiers, before the prototype: the old switch must reproduce candidate 1;
count arrays must remain identical; mathematical chain gates still apply;
removing delivered factors must make the reader invariant to perturbing those
factors, while own strand counts still affect it. On truth, recovery on the flat
row must coexist with preserved captured within-gene ratios, all three in-scope
strata, both gap directions, the single/junction probe controls and VCaP. If the
captured regression persists, removal of delivered-message reuse is insufficient;
there is no coefficient to tune. The original candidate is already refused, and
this new hypothesis is not selected in advance.

## Probability-unit audit, derived before the finite check

The generating-law correction above is a factual instrument fix. Interpreting its
raw cross-population scale ratios is a separate question. A normalized component
likelihood is `K_c(f)/L_c`. Replacing both numerator and denominator by
`K_c/z_c` and `L_c/z_c` leaves it exactly unchanged, even when `z_c` differs by
population. Consequently equality of `L_c/Y_c` across populations is necessary
only after establishing that the compared numerators share the same unit.

Construct a periodic footprint with the same number W of starts for every width,
two distinct observed-origin generating PMFs `f_g`, `f_r`, and a common capture
weight `C(w)`. Oracle raw kernels and yields are `f_c(w) C(w)` and
`Y_c = W Z_c`, where `Z_c = sum_w f_c(w) C(w)`. The observed PMF is
`h_c(w)=f_c(w) C(w)/Z_c`. A model using `h_c(w)` and `L_c=W` gives exactly the
oracle conditional fragment probability, but `L_g/Y_g` and `L_r/Y_r` differ when
`Z_g != Z_r`. This is a finite probability identity, not an estimator using an
unobservable generating law. It tests the interpretation of simulator truth.

The planned check enumerates both origins and widths, verifies normalized
probabilities exactly, and recovers a known mixed-origin fraction. It must also
show that changing the normalizer alone to force the raw class scales to agree
breaks normalization. This would qualify the raw scale readout, not excuse any
measured transcript regression or within-gene ratio damage.

## Final ablation derivation: own observational mode, without a field prior

Candidate 2 improves ladder stratum-average transcript error (stranded OFF
7.515343→2.159653%, stranded ON 9.711430→9.071331%, unstranded OFF
2.068994→2.034908%), but the captured average hides damage: g05/g50/g98 stranded
ON worsen 5.694003→11.501174, 5.944245→8.314588, and 25.413098→42.335223%.
Its g00 stranded OFF improvement is 13.609200→1.992229%, near the flat-weight
attribution ceiling. The three contaminated captured rows, and genes
1.874753→2.533589%, refuse acceptance by the average alone.

Both coupled candidates share a prior across annotation neighbours. The minimal
remaining ablation is to remove that coupling and retain only the own likelihood
`phi_i` derived above:

    x_i* = argmax_x log phi_i(x),    w_i = exp(x_i* - median_j x_j*).

This is the shipped posterior-mode estimator with its delivered composition and
level factors removed, not a new shrinkage family. There is no new prior, spatial
scale, statistic, detector, reference, or constant. Count inference remains intact.
Its rationale is specifically the measured false modes created by delivered
factors; it does not claim that the remaining own strand modes are noiseless.
The prior calibration still uses all evidence; only the reader omits these factors.

Falsifiers: base/previous arms unchanged with their switches; all six count arrays
unchanged; standalone native own-mode curves match the integrated reader's output;
the flat row's residual own-strand false modes must be priced, not hidden; all
three in-scope strata and both gap directions must hold or improve, alongside gene,
pool, single/junction-probe and VCaP readings. Deferred rows remain reported.
No improvement there is assumed. This completes the factorial attribution of
received factors versus the new coupling; there is no coefficient search.

Baseline reuse is explicit: the same native wheel is used, and all 15 comparable
ladder baseline conditions across the two builds have byte-identical six count
arrays, weight arrays, and transcript counts (maximum difference 0). The existing
fresh baseline rows will be copied with their source receipts and hashes; candidate
runs will have new receipts. The diagnostic wrapper changes only the final
reader options when its new arm is explicitly enabled.

## Local contrast: an identifiability gate, before any estimator

With opportunities absorbed into the rates, two unstranded unspliced objects
observe totals `T_E = G_E + R_E` and `T_I = G_I + R_I`. The local ratio `T_E/T_I`
is an observation, and computing it introduces no library reference or detector.
It is not, by itself, an identified gDNA enrichment ratio `G_E/G_I`.

For example, `(G_E,G_I,R_E,R_I)=(700,1,0,0)` and `(1,1,699,0)` give exactly
the same Poisson total likelihood and 700-fold local contrast. With dead strand
information and equal length laws their marked-fragment likelihoods also agree.
One explanation is enriched gDNA; the other is flat gDNA plus exon RNA. These
numbers construct a counterexample, not an estimator constant. The contrast may
be used jointly with a measured composition witness, or as its legitimately
derived one-sided bound. Treating it as a direct capture observation assumes the
answer. Removing a previous constraint does not invalidate this counterexample.

The finite check will enumerate total-count log likelihoods and two-bin fragment
marks for both explanations. Equal marks must give exactly zero distinguishing
information. Distinct *matched* origin laws may distinguish them; unmatched bulk
histograms are not those laws. This is a sub-second algebra check, not a new
prototype or a recovery claim for the deferred stratum.

## The single mechanism I would take next: shared uncertainty in allocation

This is a next-step derivation, **not a selected or implemented release fix**.
The coupled readers reduce the flat-field failure but distort captured isoform
ratios; deleting delivered factors removes real information. I would stop changing
the published marginal location and test whether allocation must integrate the
shared field. This follows from the finite covariance counterexample above, not
from another desired panel score. It does not add the refused adjacency prior.

Let `c` be one latent capture field for a locus, `a_t` the component's observed
fragment fraction, and `K_t(f;c)` its fragment kernel. Define the normalizer from
the **same** kernel and observational frame:

    L_t(c) = sum_over_admissible_fragments K_t(f;c)
    p(f | a,c) = sum_t a_t K_t(f;c)/L_t(c)
    p(D | a) = integral p(c) product_(f in D) p(f | a,c) dc.

The latent field is shared across the whole locus; it is not integrated afresh
for each fragment. Writing `q(c) = p(c | D,a_old)`, the allocation update is

    r_ft(c) = [a_t K_t(f;c)/L_t(c)] / sum_s[a_s K_s(f;c)/L_s(c)]
    expected_count_t = sum_(f in D) integral q(c) r_ft(c) dc.

Existing component priors would enter the same abundance objective they enter
today. There is no new temperature, reference density, threshold, or shrinkage
coefficient. A common multiplication of every kernel and normalizer cancels.
When capture at an observed fragment is common to its compatible components, its
numerator multiplier cancels from `r_ft(c)`, but correlated normalizers remain.
Replacing them with separate means or modes does not evaluate this integral.

The essential implementation constraint is ownership of evidence. The current
reader posterior already consumes observations in `D`; multiplying it as a prior
by the same fragment likelihood would count them twice. A valid small prototype
must enumerate each primitive observation once, use the existing prior/legitimate
licensed relationships, and derive the normalizers in their actual units. The
calibration result may initialize numerical integration; it cannot silently become
a second likelihood. Empirical-Bayes training on the same locus must be exposed
and tested, or excluded for the finite identity. A product of current marginal
posteriors is not an acceptable substitute for the shared posterior.

Start with exact enumeration of two overlapping isoforms and two atomic capture
objects, plus a gDNA component, before adding any production machinery. The first
gate is feasibility, not accuracy on the ladder. Stop if the required joint
information cannot be obtained from the current summary representation without
reusing observations. Do not turn that representation failure into an arbitrary
independence assumption or a fourth local smoother.

Falsifiers:

- The enumerated `K/L` must sum to one for every field state and every population;
  forcing raw `L/Y` agreement is not this gate.
- A point-mass field must recover the shipped conditional allocation. A common
  uncertain multiplier must cancel. The common/independent/opposed example must
  reproduce .800000/.720166/.640333 with identical length marginals.
- Grouping identical observed fragments must agree with explicit repetition;
  depositing one fragment on two objects must not create two observations.
- The exact small model must preserve an isolated captured exon and a junction
  probe, as well as a flat field, without choosing a spatial range or cut.
- If the small model passes, the next isolated A/B must improve/hold **each**
  in-scope stratum, both gap directions, captured contaminated rungs, genes,
  silent mass and pool truth. Ruler ratios must be read in matching units, all
  zero and deferred controls retained, and VCaP checked in its correct truth
  domain. Only then would a runtime/memory budget for a scalable implementation
  be justified. A win on g00 alone is already known to be insufficient.

This is a narrower question than redesigning all calibration: does carrying the
correct shared uncertainty into allocation remove the observed rare-mode damage?
Its cost and data representation are unresolved; no release improvement is claimed.
