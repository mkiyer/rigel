# The shared capture field inside the allocation — design and implementation plan (2026-10-11, sandbox, for external review)

> **Follow-up review, 2026-10-11:** this submitted plan is preserved as reviewed, but its inference and
> representation are **not approved for a literal production port**. See `SHARED_CAPTURE_FIELD_REVIEW.md`
> for the six decisions, measurements and a corrected executable reference. The integrated E-step and
> sufficient-statistic mean-field update below optimize different objectives; the second-order shortcut
> fails on broad posteriors; conserved deposits require candidate-path share vectors; and the current
> scorer/geometry length-law split does not establish the claimed normalized kernel. A junction cap is
> also not a universal truth constraint. Collapsed posterior sampling repairs the inference update on
> enumerable models, but the pipeline representation, common kernel measure and full accuracy A/B remain
> unfinished. No production inference change or release-readiness claim follows from the toy results.

*Status: DERIVED and screened on an exactly enumerable locus; nothing built in the tree. This is the mechanism the
2026-10-11 implementation review (`CAPTURE_IMPLEMENTATION_REVIEW.md`, this directory) named as its single next step,
carried from its one-paragraph derivation to a design the tree can hold and a plan with falsifiers, budgets and
stages. It is submitted for external review before any prototype touches production structures. Nothing here is
authoritative; the rulings it rests on are cited from the permanent docs. State described: `main` at `624dfade`.*

**What we ask of the reviewer.** Break the derivation (§2–§3), the evidence-ownership argument (§2.5), the
approximation and its gates (§3.3, §6), the claim of what the mechanism can and cannot identify (§2.4), and the cost
(§5). The constraints of §9 are fixed.

## 0. The mechanism in one page

Today the EM contracts every component's effective length by capture weights read one object at a time by a
separate stage (the reader: each object's posterior mode of its own gDNA density under the landscape, with the
RNA amount integrated out and the solve's delivered factors multiplied in), then treats those weights as known.
The review measured what that costs: on a flat field the delivered factors create rare false modes at objects of
any depth (696 of 70,176 log weights above 1 on the g00 stranded ladder row), and every per-object repair that
moved those modes — two posterior averages, a local prior, own-only evidence — fixed the flat rows and broke the
captured ones. It also showed, by a finite identity, why no per-object average can carry the uncertainty: the
responsibilities depend on how the components' length uncertainties are *correlated* through the objects they
share (identical marginals give 0.800, 0.720 or 0.640).

The mechanism: treat the capture field of a locus as **one latent variable shared by every fragment of the locus**,
give it the landscape as its prior, and integrate it inside the allocation. Each fragment then enters the
likelihood exactly once, through the same factor the EM already uses, with its rate under every component scaled
by the field at the object it deposits on and normalised by that component's field-contracted length. The reader
disappears as a separate likelihood: the field is inferred jointly with the abundances from the same fragments,
with the gDNA component's uniform-before-capture structure and the isoforms' shared exons as the two sources of
identification. On an exactly enumerable locus (§7) the mechanism passes the review's six falsifiers and does one
thing no gDNA-based reader can: at a captured exon in a library with no gDNA it reads the capture from the isoform
structure itself (allocation error 0.63 → 0.05 against an oracle's 0.04).

The claim it must earn, per stratum and before anything lands: the flat rows return toward the attribution
ceiling (g00 stranded OFF 13.6 → ~2 %) **and** every captured row, both gap panels, the probe layouts, the zero
controls and VCaP hold. The honest cost is a Monte-Carlo E-step, roughly S times today's E-step, S between 8 and
16, with a second-order expansion as the candidate to bring it near today's cost if it passes its own gate.

## 1. The facts this rests on (all measured, 2026-10-11)

- The capture weights own the stranded strata's regression: g00 stranded OFF transcripts 13.61 % → 1.93 % with
  every weight forced to one, calibration unchanged. The damage is isoform allocation through the contracted
  lengths (the gene level holds).
- The false weights are rare mode switches created by delivered composition and DNA-level factors, at any depth;
  removing those factors repairs the flat rows and breaks the captured gap rows, so they carry real information on
  captured libraries.
- A gDNA-based reader cannot see capture where there is no gDNA (the g00 capture-ON rows: 15.1 % stranded) nor on
  unstranded data (maximum log weight 0.00005), because its only handle is the gDNA's uniformity before capture.
- Responsibilities under uncertain lengths depend on the correlation of the components' lengths through shared
  objects; independent per-object averages cannot represent it (`uncertainty_identity` in the review archive).
- The conserved-share length rule (`EQUATIONS.md` §11; `calibration/capture_eff_length.py`): a component's EM
  length is its plain fragment-length-marginal length times the share-weighted mean of the weights over its pieces
  and cuts; a junction cut's weight is derived from the objects beside it and capped at the most captured object
  it touches.
- The EM (`native/em_solver.cpp`, `estimator.py`): per locus, units (fragments or equivalence classes) carry a
  candidate component set and a per-candidate log-likelihood from the scorer (fragment length, strand, position)
  plus the gDNA component's; the E-step weights each candidate by `log θ_k − log L_k` and normalises; VBEM with
  the pseudocounts of `DESIGN.md` §3.1c; SQUAREM; the count-first draw of §3.1f. A unit carries no object index.

## 2. The generative model

### 2.1 Objects, components, field

A locus has objects `o`: its regions (pieces) and its cuts (boundaries between consecutive pieces, junctions
among them) — exactly the objects of the conserved-share rule (`TranscriptObjects`: pieces with contained shares,
cuts with cut shares). Components `k`: the annotated transcripts, the synthetic nascent entities, and the locus's
gDNA component `G`. The field `c = (c_o)_o`, `c_o > 0`, is the capture multiplier at object `o` relative to the
library's typical read object (the median unit of today's weights, `DESIGN.md` §7.2); only ratios of `c` matter
to the RNA components, and a common factor cancels everywhere (gate 5).

### 2.2 Lengths as functions of the field

Every component's contracted length is today's rule with the field in place of the plug-in weights:

    L_k(c) = fl_k · Σ_o s̃_k(o) · c_o,        Σ_o s̃_k(o) = 1,

`fl_k` the component's plain length and `s̃_k(o)` its normalised conserved shares over its pieces and cuts (for
`G`: its opportunities over the locus's objects, as `assemble_priors` forms `gdna_eff_len`). A cut that is a
junction is an object of its own with its own `c`; §4.3 says what becomes of the cap rule.

### 2.3 Fragments

A fragment `f` deposits on one object `o(f)`: the region that contains it, or the cut it crosses (the deposit rule
of `DESIGN.md` §3). Its rate under component `k` is

    rate_k(f; a, c) = a_k · φ_k(f) · c_{o(f)} / L_k(c),

`a_k` the component's abundance (its share of the locus's fragments, the EM's θ), `φ_k(f)` the scorer's existing
per-fragment likelihood under `k` (the unit's `log_liks` for transcripts, `gdna_log_liks` for `G`: fragment length,
strand with the protocol's κ for RNA and ½ for gDNA, position). The normalisation is the conserved-share identity:
summing `φ_k(f) · c_{o(f)}` over the fragments `k` can emit gives `L_k(c)` (gate 1). The likelihood of the data
given the abundances and the field is one factor per fragment,

    p(D | a, c) = Π_f  c_{o(f)} · Σ_k a_k φ_k(f) / L_k(c),

and the marginal the allocation maximises is

    p(D | a) = ∫ p(c) · p(D | a, c) dc,        p(c) = Π_o p_o(c_o).

This is the review's form, with the kernel and normaliser in the same units (both are the scorer's φ times the
field, summed over the same admissible fragments).

### 2.4 What identifies the field

The numerator `c_{o(f)}` is common to every component compatible with `f` (the deposit object is genomic), so it
cancels from every responsibility: the field enters the allocation **only through the components' lengths**, and
the components' lengths are correlated exactly through the objects they share. The field itself is identified by
three things and nothing else:

1. **The gDNA's uniformity before capture.** `G`'s share of object `o` is its opportunity; a count at `o` above
   `a_G · s̃_G(o) / L_G(c)` that the strand split calls gDNA is capture. This is the reader's handle, now inside
   the joint model.
2. **The isoform structure.** The objects two isoforms share fix their abundances; an object one of them owns
   with a count above that abundance's share is capture. This is what the enumerable locus shows at zero gDNA
   (§7), and what no per-object reader can do.
3. **The prior**, where 1 and 2 are silent (a single-exon gene with no gDNA): the posterior is the prior and the
   weight is the typical object's — the correct statement of ignorance, with no false mode.

Two blindnesses are unchanged by design and must be said: on unstranded data the strand split is dead and `G` is
told from RNA by length and position only; and a junction with no spliced fragments and no gDNA has no evidence.

### 2.5 Evidence ownership

Each fragment contributes exactly one factor. Nothing from the reader enters: its posterior consumed the same
fragments through the strand columns, and multiplying it in as a prior would count them twice. The solve's
delivered factors (composition rows, held levels) do not enter either; the neighbour information they carry
travels instead through the shared parameters — one abundance per isoform across all its exons, one `a_G` and one
`L_G(c)` per locus. Whether that recovers what the delivered factors contribute on captured rows (the review's
RNA-short stranded ON finding) is a stage-2 measurement, not an assumption. The landscape prior is empirical
Bayes on the same library; §8 says how that is exposed.

## 3. Inference

### 3.1 The bound

With `q(c) = Π_o q_o(c_o)` on a grid per object, maximise

    F(a, q) = E_q[ log p(D | a, c) ] − KL( q ‖ p ).

The EM bound on the log-sum over components, with responsibilities `r_fk`, makes the field's part of `F` a
function of two sets of sufficient statistics only — every object's fragment count and every component's total:

    E_q[ log p(D | a, c) ]  ≥  Σ_o n_o E_q[log c_o]  −  Σ_k N_k E_q[log L_k(c)]  +  Σ_f Σ_k r_fk [ log a_k + log φ_k(f) ]  −  Σ_f Σ_k r_fk log r_fk,

`n_o` the number of fragments deposited on `o` (the E-step's units by multiplicity plus the deterministic spliced
fragments, which skip the E-step and enter every M-step, `DESIGN.md` §3.1c), `N_k = Σ_f r_fk` plus the
deterministic counts. The field never needs the fragment list: it needs `n_o`, `N_k`, the shares `s̃_k(o)` and the
prior.

### 3.2 The three updates

- **Allocation (E-step).** `r̄_fk = E_q[ r_fk(c) ]` with `r_fk(c) ∝ a_k φ_k(f) / L_k(c)` — the responsibility
  integrated over the shared field, which is the review's requirement and the thing no plug-in does.
- **Abundances (M-step).** Unchanged: the grouped update with the pseudocounts of `DESIGN.md` §3.1c, VBEM's
  Dirichlet update, SQUAREM on the abundance state as today (`DESIGN.md` §3.1e; the accelerator never decides).
- **Field.** Coordinate ascent on each `q_o` over its grid `{c_g}`:

      log q_o(g) = log p_o(c_g) + n_o log c_g − Σ_{k ∋ o} N_k E_{q_{−o}}[ log L_k(c) | c_o = c_g ] + const.

  The expectation over the other objects is over the same field draws the E-step uses (§3.3), so the field step
  costs `O(objects × grid × S × components per object)` — negligible.

The schedule per locus: E-step, M-step (with SQUAREM), field step, repeat to the EM's convergence rule; then the
count-first draw of `DESIGN.md` §3.1f from `r̄`. The published capture weight per object is the marginal mode of
`q_o` (the location-faithful readout the review kept), relative to its median over the library's read objects;
the prototype also reports `E_q[log c_o]`.

### 3.3 Evaluating the integrated responsibility

`r̄_fk` depends on the joint distribution of `(log L_k(c))_k` under `q`, which has no closed form.

- **Scheme A, the reference: Monte Carlo with common random numbers.** Draw `S` fields per locus from `q`
  (per object independent on its grid), compute every component's `log L_k(c^(s))`, and average the softmax over
  `s`. The same draws serve every fragment of the locus and the field step; the seed is the EM's own
  (`SplitMix64`, `DESIGN.md` §3.1f), so a run is reproducible. Cost: `S` evaluations of the E-step's inner loop.
- **Scheme B, the candidate: a second-order expansion.** Around `μ_k = E_q[log L_k]`, with the mean-field
  covariance `Σ_km = Cov(log L_k, log L_m)`, which is exact in closed form through the shared objects
  (`Cov(L_k, L_m) = fl_k fl_m Σ_{o ∈ k∩m} s̃_k(o) s̃_m(o) Var_q(c_o)`), `r̄ ≈ r(μ) + ½ tr(H Σ)` with `H` the
  softmax's Hessian in the log-lengths — `O(candidates²)` per unit, candidates being two to five. It carries the
  correlation the review insists on and costs about twice today's E-step. It is an approximation, and it ships
  only if gate 3b holds: agreement with scheme A on the enumerable locus and on the test chromosome's loci to the
  gate's tolerance, including where `q` spans orders of magnitude.
- **Adaptive S.** A locus whose every `q_o` is a spike (its entropy below the grid's resolution) needs `S = 1`,
  which is the plug-in exactly; the cost concentrates where the uncertainty is. Measured, not assumed.

### 3.4 The prior

`p_o(log c) = P_land(log c + x_med)`: the fitted landscape (`DESIGN.md` §7.1) shifted so that the typical read
object's density is `c = 1` — the same unit the weights are published in today. No new density, no constant. Two
derived alternatives are kept for the prototype's sensitivity gate: the landscape's own mode as the unit, and the
owner-approved local enrichment prior (`p(C) = C^−2, C ≥ 1`, `DESIGN.md` §7.2) on the background's mode. The
landscape is fitted on the same library's deconvolved gDNA; the enumerable model excludes that (its prior is
given), and stage 2 measures the sensitivity by perturbing the prior's location and width.

## 4. What it replaces, what stays, what it may delete

### 4.1 Stays

The scan and the deposit rule; calibration entire (the count foundation, the landscape, the strand model, the
messages); the scorer's `φ`; the EM's components, pseudocounts and odds (§3.1c), the whole-count draw (§3.1f), the
accelerator (§3.1e); the conserved-share definition of a length; the published capture surfaces' meaning (a
weight per object relative to the typical read object).

### 4.2 Replaced

The plug-in: `transcript_capture_eff_lengths` evaluated once at the reader's modes becomes `L_k(c)` evaluated
under `q` inside the EM. The reader as a separate likelihood stage (`native/honest_reader.h`, the reader loop in
`solve_kernel.cpp`, `calibrate.capture_weights`) is then redundant for the EM. The design proposes deleting it
once the field's weights are published from the EM and the A/B holds — one path, no legacy — and recording the
reader's rulings as superseded in `DESIGN.md` §7.2. If the reviewer or the owner prefers the reader kept as a
diagnostic surface, that is a one-line decision; the design does not depend on it.

### 4.3 The junction cap

Today a junction cut's weight is derived from the objects beside it and capped at the most captured object it
touches, because a junction has no gDNA to read and the simulator's law captures a fragment by its best single
probe part (`DESIGN.md` §7.2, 2026-10-05). In the field model a junction cut is an object with its own `c`, its
own prior and its own evidence — the spliced fragments at it against the isoforms' abundances (identification 2
of §2.4). The cap becomes a falsifier: on the junction-probed layout the field must not read a junction above the
objects it touches by more than the pair's own variation; if it does, the cap is re-derived as a constraint on the
cut's prior support, not as a rule applied after the fact. This is flagged for the reviewer (§8, question 4).

### 4.4 Units gain an object

A unit (fragment or equivalence class) needs its deposit object to form `n_o`. The placement is in the fragment
buffer by `frag_id`, and the object bounds are calibration's chain; the second pass already performs this lookup
(`second_pass._resolve_intron_lookups`, `placement_of_hypothesis`). The pipeline attaches an `int32 object` per
unit at scoring, and the deterministic spliced fragments' cut objects likewise. Equivalence classes merge units
with the same candidate set *and* the same object (a small increase in class count, measured).

## 5. Representation and cost

Per locus: the objects of its components (tens to a few hundred); the sparse share matrix `s̃_k(o)` from
`TranscriptObjects` and `G`'s opportunities; `n_o`; the field grid per object (the landscape's support at the count
lattice's 0.2-nat step: about a hundred nodes — a derived lattice, ψ's own); `S` field draws (`objects × S`
doubles). The native E-step kernel (`em_step_kernel_range`) takes `S` weight vectors instead of one and loops over
`s` inside the row; the field step is a per-locus pass. Thread identity holds as today (per-locus determinism, the
seed per locus).

Cost on VCaP (22,787 loci, 17.3 M ambiguous fragments; the EM today about 27 s of wall on 16 threads, measured
2026-10-11): scheme A at `S = 16` is about 7 min of wall; at `S = 8` about 3.5 min; scheme B about 1 min; adaptive
`S` lands between, depending on how many loci carry broad field posteriors — on a capture-OFF library many
shallow objects do, so the saving is not assumed. These are estimates from the kernel's cost model, to be replaced
by measurement at stage 3; a performance claim comes only from a back-to-back pair on an idle machine.

## 6. Gates — each written before the code it gates, and shown to fail under a deliberate defect

| gate | holds | the defect that must fire it |
|---|---|---|
| 1 · normalisation | `Σ_f φ_k(f) c_{o(f)} / L_k(c) = 1` for every field state and component, on the enumerable locus and on real `TranscriptObjects` | a share omitted from `L_k` |
| 2 · point-mass identity | with every `q_o` a point mass at today's reader mode (`S = 1`), the transcript counts equal the shipped pipeline's (the Python prototype to 1e-6 relative; the native port bit-identically) | the field read at the wrong object |
| 3a · the correlation identity | 0.800000 / 0.720166 / 0.640333 with identical marginals (the review's finite check) | independent averaging of the lengths |
| 3b · scheme B against scheme A | responsibilities and transcript counts within the gate's tolerance on the enumerable locus and the test chromosome, including loci where `q` spans orders of magnitude; else B is refused and A ships | the covariance term dropped |
| 4 · one observation per fragment | `Σ_o n_o` equals the locus's fragment count; a spanning fragment counted once on its cut; grouping identical units equals repetition | a fragment deposited on both of two objects |
| 5 · invariances | a common multiplier of every `c` leaves every count unchanged; a translation of the grid translates every published weight by the same amount; a zero-width prior reproduces today's plug-in | the unit dropped from the prior shift |
| 6 · the approximation's price | mean-field + scheme A against exact enumeration on the enumerable locus: responsibilities within tolerance, the zero-gDNA captured exon still recovered | a cross-object correlation planted in truth |
| 7 · identification where it is silent | a single-exon gene with no gDNA reads the prior's typical weight; a deep uncaptured object stays at 1; an isolated captured exon keeps its weight with flat neighbours; the single-probe layout survives | the prior replaced by a flat one |
| 8 · the deliverable, per stratum | pinned, fractional, both arms from one build: every in-scope stratum, both gap panels, the single-probe and junction-probed layouts, the zero-gDNA and `.70`-specificity controls, VCaP against its origin reference and its byte-identical baseline; the flat row not worse than the weights-forced-to-one arm by more than the pair's variation; no captured contaminated row worse; the deferred rows reported, never pooled | — |

## 7. The enumerable locus (done, `prototypes/2026-10-11_shared_field/shared_field_model.py`, 5 s)

Objects A (exon shared by T1 and T2), AB (the boundary between A and B, a spanning fragment deposited once), B
(exon of T2 only), I (intron); components T1 = {A}, T2 = {A, AB, B}, G over all; the field on {1, 3, 10, 30} per
object (256 states) with a spectrum prior (0.85, 0.10, 0.04, 0.01); fragments by (object, strand) with s = 0.99 for
RNA and ½ for gDNA. Gates 1, 2, 3a, 4, 5 and the one-observation gate pass to machine precision. Experiments (60
replicates, 300 fragments per unit rate; the plug-in reads each object's field from its own two columns by the
posterior mode with the RNA amount integrated out — the shipped reader's logic without delivered factors):

| truth | oracle | plug-in | shared field |
|---|---|---|---|
| flat field, zero gDNA | 0.082 | 0.082 | 0.084 |
| flat field with gDNA | 0.101 | 0.101 | 0.107 |
| one exon captured 10×, with gDNA | 0.094 | 0.110 | 0.101 |
| one exon captured 10×, zero gDNA | 0.043 | 0.626 | 0.050 |

(Σ|Δ abundance| over the three components.) On a flat field with own-only evidence the per-object modes are all
correct, so the production flat-field defect — the delivered factors — is not in this model and the model cannot
show that repair; what it shows is the identification by isoform structure (last row) and that integrating the
field costs nothing where the plug-in is right.

## 8. Risks, stated

1. **Mean-field loses cross-object correlation.** Adjacent objects under one probe are positively correlated in
   truth (the footprint); `q = Π q_o` cannot share strength between them except through the components. Gate 6
   prices it on the enumerable locus with a planted correlation; a structured `q` over adjacent pairs is the
   fallback, not a prior.
2. **The landscape is empirical Bayes on the same library.** Stage 2 perturbs the prior's location and width and
   reports the counts' sensitivity; the enumerable model is free of it.
3. **Cost.** §5. Scheme B or adaptive S, or the honest answer that the mechanism costs minutes on a deep library.
4. **The deferred stratum and unstranded data** are not helped: the strand split is dead there. Said in §2.4.
5. **The junction cap** (§4.3).
6. **Every number moves.** The goldens, the ladder report and every receipt are re-recorded; the A/B is the
   argument, not the goldens' update.
7. **The draw.** `r̄` is a mixture of responsibilities; the count-first draw uses it as today's posterior. The
   whole-count gates (`tests/test_estimator.py`) hold unchanged.

## 9. Constraints that are fixed

No capture detector, library capture status or reference density; no new constant (the grid is the count
lattice's, the prior the landscape's, `S` is a measured budget, not a tunable — the prototype reports the
convergence of the counts in `S`); no tuning to the panels; every score per stratum, one mechanism per A/B;
numeric no-ops proven before speed work; the owner drives commits.

## 10. Stages and budgets

0. **This document; external review.**
1. **The enumerable locus, extended** (hours, no production code): mean-field + scheme A against exact
   enumeration (gate 6); scheme B against A (gate 3b); a planted cross-object correlation; the single-transcript
   locus and the zero-gDNA single-exon gene (gate 7); the junction object; the prior's sensitivity.
2. **The Python prototype on the test chromosome** (days, outside `src/`): per-locus structures built from the
   pipeline's own objects — `ScoredFragments` with the attached deposit object, `TranscriptObjects`,
   `assemble_priors`' gDNA shares, the landscape from `CalibrationResult` — the joint EM in NumPy at the
   equivalence-class level with scheme A (`S` = 8 and 16, the counts' convergence in `S` reported) and scheme B;
   gates 1, 2, 4, 5, 7; then the A/B on the test chromosome's strata, the single-probe and junction-probed layouts
   and the `.70`-specificity controls (seconds per condition). Budget announced before each run.
3. **The ladder and the gap panels through the prototype** (hours): 16 + 8 conditions, scheme A at the `S` stage 2
   selected; the deliverable gate 8 read per stratum; VCaP through the prototype if the panels hold (one
   whole-genome process, budgeted).
4. **The native port** (days): `em_solver.cpp` — the unit's object array, `S` weight vectors in the E-step kernel,
   the field step per locus, the published weights from the field; thread identity; the Python prototype as the
   port's executable specification (bit-identity of the field step, 1e-9 on the counts for scheme A with the same
   draws); the reader stage deleted or kept as a diagnostic (§4.2); the suite, `preflight.py --full`.
5. **The landing**, by the move rule: `DESIGN.md` §3 and §7.2 (the field replaces the reader's rulings as the
   capture mechanism), `EQUATIONS.md` §11 (the joint likelihood, the bound, the field step, scheme A/B),
   `ISSUES: the-capture-weights-are-noisy-on-a-capture-off-library` closed with the killing numbers or kept with
   the residual; the goldens re-recorded with the diff read; the ladder report republished. The commit is the
   owner's.

## 11. Questions for the reviewer

1. Is the bound of §3.1 the right objective — in particular, is the field step's use of the EM bound (the
   log-sum replaced by its responsibility-weighted expectation) acceptable, or must the field step use the exact
   `log Σ_k`? The latter costs a fragment loop per grid node per object.
2. Scheme B: is a second-order expansion in the log-lengths acceptable for production if gate 3b holds, or should
   scheme A ship with adaptive `S`?
3. The prior: the landscape relative to the median unit (§3.4), or the approved `C^−2` enrichment prior on the
   background's mode? The first is what ships today; the second is what you approved on 2026-10-09.
4. The junction cut as a free field object (§4.3): acceptable, or keep the cap as a hard constraint from the start?
5. Deleting the reader stage once the field publishes the weights (§4.2): premature, or the right convergence?
6. Anything in §2.4's identification claim you believe is wrong — above all whether the isoform-structure
   identification survives real annotation (overlapping genes, many isoforms sharing most exons, synthetic
   nascent spans over introns that also carry gDNA).

## Appendix — notation and the production structures

`o` object (region, boundary, junction cut); `k` component (transcript, synthetic nascent entity, `G`); `f`
fragment; `o(f)` its deposit object; `φ_k(f)` the scorer's per-fragment likelihood (`ScoredFragments.log_liks`,
`gdna_log_liks`); `a_k` the abundance (the EM's θ, VBEM's α); `s̃_k(o)` the normalised conserved share
(`TranscriptObjects.contained`, `cut_share`; `assemble_priors` for `G`); `fl_k` the plain length
(`TranscriptGeometry.effective_lengths`); `c_o` the field; `L_k(c)` the contracted length; `n_o`, `N_k` the field's
sufficient statistics; `q_o` the variational marginal on the grid; `S` the field draws; `P_land` the landscape
(`DensityLandscape`); `x_med` the median read object's log density (`calibrate.capture_weights`'s unit).
