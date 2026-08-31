# ISSUES — the issue log

⭐ **WHAT THIS FILE IS.** One entry per problem, question, decision or risk, in two sections: **OPEN**
(what could still change the tool) and **CLOSED / REFUSED** (a permanent, append-only record of what was
measured and turned down — the reason a refused mechanism is not rebuilt). ⭐ `ROADMAP.md` is the short
ranked view and points here; this file holds each issue's substance so the roadmap never grows.
⛔ **The changelog is git** — nothing here records what was done, only what is open and what was refused.

⭐ **An issue is keyed by its NAME** — a kebab-case heading, greppable as one string — never by a number
(`tests/test_no_jargon_labels.py` and the rule in `CLAUDE.md`). Cite as `ISSUES: <name>`.
Each open entry carries `priority` (now / next / later / parked), a `kind`, and its stamp.

⭐ **THE NUMBERS POLICY APPLIES** (owner, 2026-08-22): an OPEN entry states direction + magnitude class
and names the instrument that re-derives the figure; a precise number stays only when a ranking turns on
it. ⚠ CLOSED entries are the deliberate exception — a refusal keeps its stamped measurements forever,
because a graveyard row without its number is an invitation to rebuild.

---

## OPEN

### rename-the-drain
`priority: later · kind: decision · stamped: 2026-08-31 (owner: "might consider")`
"Drain" carries no intuition — the concept is: pass one BUFFERS any fragment whose unsequenced mate
gap admits more than one explanation (an annotated intron may lie in the gap), and the second pass
performs DATA-DRIVEN ASSIGNMENT of each buffered fragment against the whole library's densities and
length models. The buffer is already well named (`payload.deferred`, `DeferredFragments`); the verb
family is not (`drain`, `DrainQC`, `_drain_side_buffer`, `payload.drain`, `lift_choices`' docs).
⚠ ~530 sites across src/scripts/tests/docs — `arm`-rename scale, so it needs its own
`rename_census.py` sense registration and a staged pass with `rename_identity.py --freeze/--check`,
never a tail-end sweep. ⚠ Candidate verbs and their collisions: *resolve* (collides with
`resolve.cpp`, fragment construction), *assign* (collides with the EM's fragment→transcript
assignment), *settle*, *place*, *adjudicate* — the owner picks. The QC dataclass and the
`payload.drain is not None` frame test rename with it.

### drain-contaminates-certified-rna
`priority: next (owner, 2026-08-31: "extremely concerned … avoid all sources of false positives") · kind: defect · stamped: 2026-08-31`
**The second-pass drain deposits some TRUE-gDNA fragments into the certified-RNA banks** — the
whole-library drain draws a spliced hypothesis for a held gDNA fragment whose mate gap admits an
annotated intron, so production's tally violates "gDNA cannot splice" as a statement about DEPOSITS.
⭐ Proven independent of the lift: on `flgap_rna_short g50 ss.99 OFF` the lift ambiguity is exactly 0
and the leak is still 15 records; on `flgap_rna_long` ambiguity is 5,805 and the leak is 0 (a long-RNA
library never mistakes a gDNA length for a spliced one). Measured on the ladder: 233 records at
`g50 ss.99 OFF` (~1e-4 of the certified channel) but **1,482 at `g98 ss.99 ON` — 1.9 % of that
condition's whole certified-RNA channel**, an IN-SCOPE stratum. ⚠ The certified-flux anchor and every
"spliced ⇒ certainly RNA" consumer treat this channel as exact; at high gDNA under capture it is ~98 %
pure, not 100 %. Invisible to every undrained-frame instrument, which is how it went unmeasured
(the frame ruling is `DESIGN.md` §4.3 — landed, so the leak is now visible: `calibration_vs_oracle` reports it per row and every certified
`slot_truth.npz` carries the drained-frame report verdict).
⭐⭐ **DERIVED 2026-08-31, and the measurements FORCE the design** (the derivation and numbers:
`docs/dev` carries the plan; prototypes gated on byte-identity with production's drain choices).
① The leak is EXACTLY posterior sampling — realized gDNA-in-spliced matches Σ P(spliced|record) over
true-gDNA records within the draw's own noise on every condition measured (0.0σ / 1.4σ / 0.7σ) — so
"make the drain smarter" has no headroom, and a posterior-odds floor is dominated (the posterior is
calibrated; a floor trades true splices for gDNA at its own rate, and is a tuned constant besides).
② A provenance split is refused by arithmetic: at `g50 ss.99 OFF` it evicts 134,850 correct drained
records from the certified channel to remove 75 contaminants, resurrecting the −4 bp spliced-pool
bias the drain exists to repair. ③ No structural gDNA channel: records with no genomic survivor are
100 % mature. ④ The harm is concentrated and qualitative: at `g98 ss.99 ON`, 137 boundaries carry
certified RNA where truth has ZERO spliced RNA. ⑤ The leakers' own P(genomic) is known at drain time
(median 0.494 at `g98 ss.99 ON`) — the information to price the false certainty already exists.
⭐ **The design**: deposits unchanged; the drain accumulates per certified object
`E_genomic = Σ q_null` and `V_genomic = Σ q(1−q)` (its own posterior, no new model, no constant), and
CALIBRATION — which owns composition — converts path-contamination to gDNA-contamination with its own
local `f_g`. Consumers correct the certified count's mean and/or precision. Score per stratum on
`calibration_vs_oracle.py`, both zero controls, both fl-gap sign arms, AND the spliced-pool means
(the drain's length-law repair must survive). Sequencing and open points: the plan doc in `docs/dev`.

### measured-prior-rung-4
`priority: now · kind: build · stamped: 2026-08-26`
**The measured prior — ψ's proper prior FITTED from composition-free observables before any solve**
(`DESIGN.md` §3.1a-i–ii, `EQUATIONS.md` §2.3b). Rungs ① MEASURE (`total_abundance`), ② INTEGRATE
(`fit_intron_background`, behind `CalibrationConfig.background_abundance`) and ③ FIT
(`AbundanceLandscape`, censused by `abundance_landscape_census.py`, A/B'd by `landscape_head_to_head.py`)
are done. **Rung ④ consumes `rho_0`, the per-class enrichment responsibility `w`, and the enrichment
DETECTOR — NOT `span_R`**, which is grid-fragile (`TRAPS: a-mode-count-is-not-a-well-posed-quantity`).
It lands behind its own `CalibrationConfig` flag (`composition_reference`, default bit-identical); the
door is ψ's location term widening from one scalar to `(m_lo, m_hi, w)`.
⭐ **The measured requirements spec** (each a measurement): ① the location's range must span BOTH lattice
ends; ② `rho_bg` must be ON-RATE for the destination's enrichment class — the intergenic pool carries
exactly zero probe bases and is the UNPROBED rate, the single named cause of every capture-ON failure;
③ exact at `g00` (structurally true today — keep it as a free falsification); ④ strength stays one
pseudo-fragment; ⑤ reference and message interact and must be priced jointly; ⑥ both zero controls on
every arm, pure-gDNA/RNA-bearing never pooled. Rung ⑤ prices per stratum against a SHUFFLE control
(`TRAPS: attribution-must-survive-a-shuffle`), split by destination-had-own-evidence.
⛔ **Read `ISSUES: reference-prior-refuted-at-concept-level` first** — the location-tilt concept is
refuted; rung ④'s reference must be a reparameterisation or a density-channel extension, not a better
tilt. ⛔ Constraints from the 2026-08-22 audit: ⓐ the intron-inclusive anchor pool is CONTAMINATED on the
sparse panel (`object_composition.PURE_GDNA_STRATA` includes `R intron`; shipped `fit_intron_background`
pools intergenic-only and is clean — the defect is the INSTRUMENT's pool); ⓑ the enrichment detector
survives as a BOOLEAN, never a calibrated LEVEL (`TRAPS: a-total-density-ratio`). ⛔ Refused cheap
alternatives (measured): the library-wide truth mean, both single-mode per-object locations, and ④a —
widening the shipped estimator's selector to `solvable_exon`, BUILT, priced on all 16, DELETED 2026-08-24
(wins all 8 capture-OFF rungs and all four `g00` rows exactly, loses all 6 contaminated capture-ON rungs
— the cleanest statement of why a data-derived location fails).

### reference-prior-refuted-at-concept-level
`priority: now · kind: design-constraint · stamped: 2026-08-24 (owner + external review)`
**The reference prior's location tilt is refuted at the concept level — do not repair it.** The root
cause is a COORDINATE mismatch: on the λ = logit(f_g) axis the data's information
`I ∝ N_eff·disc·[f_g(1−f_g)]²` vanishes at the vertices and is identically zero at κ = ½, while the
tilt holds fixed nats there — the prior-to-data ratio diverges exactly where truth lives. The information
does NOT vanish in `f_g`-space (`I_{f_g} = N(½−κ)²/(p(1−p))` is bounded away from zero), so
`[f_g(1−f_g)]²` is purely the Jacobian. Measured: the tilt is overturned at N = 3 fragments where the
strand channel is alive and NEVER at any depth at κ = ½ (0.7471 from N = 10 to N = 10⁶); 82–95 % of
scored pass-0 error at unstranded conditions sits on slots the tilt decides.
⛔ Measured-refused repairs: a different constant (0.75 optimal at none of 16; `σ(L)` 3.3×/10.8× worse
at the zero controls; the required location changes SIGN across the panel); a data-derived location
(circular; its exon form built, priced, DELETED); re-weighting the tilt (helps only zero controls, up to
2.8× worse elsewhere); an information-weighted tilt `w = 1 − exp(−γI)` (introduces γ; vanishes at true
vertices). ⭐ **The candidates worth prototyping**: a REPARAMETERISATION away from λ (judged on
`L`-invariance and the Berger–Bernardo cancellation BEFORE any panel number) and extending the density
CHANNEL (`density_deconv.density_lambda_factor`, already a likelihood, ships at ss-intron REGIONs).

### message-value-for-blind-slots
`priority: now · kind: question · stamped: 2026-08-27`
**The decisive message-policy question, answerable with no solver**: is the information a blind
unstranded or AMBIG slot needs actually PRESENT in its neighbours? The anchored twin block was designed
for exactly this. The measurement that survived the 2026-08-27 tear-down: propagation is net-harmful
wherever the local solve HAS evidence, and its value concentrates where the solve is BLIND — so the
follow-up is whether a policy that speaks ONLY into destinations with no own evidence beats one that
speaks everywhere. Score split by destination-had-own-composition-evidence on `policy_benchmark.py`,
never pooled, substrate named on every claim. ⭐ The laws a policy must obey are NAMED TRAPS rules
(`zero-the-precision-with-the-value`, `an-imputation-must-cost-something-every-hop`,
`off-grid-message-mode`, single-source-may-only-reduce — the foundation spec enforces the last at
runtime). ⚠ Re-baseline first: all policy numbers predate the strand-estimator and fragment-length work.

### per-transcript-prior-lane
`priority: next · kind: build · stamped: 2026-08-31`
**The largest in-scope end-to-end lever measured, and the weighting function IS the work.**
`rna_prior_weight` is built end to end; the production call site in `pipeline.py` omits it, so the
shipped EM carries zero per-transcript information. A perfect per-transcript prior roughly halves
in-scope gene-level error and cuts false-positive mass by orders of magnitude; it does NOT rescue the
deferred stratum, so the prize is in scope. ⛔ Two weighting functions are built and refused
(`ISSUES: refused-transcript-weights`, `ISSUES: refused-soft-min-path-weighting`, CLOSED;
`TRAPS: an-upper-bound-is-not-an-estimate`); what they settled: the target is TOTAL abundance and **the
support problem is the whole problem** — no expressed transcript is wrongly zeroed, thousands of silent
ones are not. The next candidate must be a SPARSITY mechanism, scored by `quant_accuracy.py`'s existing
arms. ⭐ **New ammunition (2026-08-31)**: the rna-pmf dose–response (CLOSED, `rna-length-law-fix`)
measured a lower bound on what EM assignment at short-exon genes is worth — a crude prior nudge moves up
to ~7 % of transcript error per condition, concentrated on expressed multi-exon transcripts with median
exon ≤150 bp, capture-independent. That is this issue's target population.

### performance-memory-bounded-solve
`priority: next · kind: build · stamped: 2026-08-17 (owner: mandatory before 0.8.0)`
The grid solve in memory-bounded parallel chunks, with an advanced CLI flag spanning one object at a
time through many in parallel; compute-versus-memory is the dial. ⛔ Optimise on HIGH-DEPTH REAL RNA-seq,
not cfRNA (this REVERSES the older "profile on cfRNA" instruction). ⚠ Calibration is depth-independent
and scales with the INDEX while the EM scales with the DATA — the panel profile and the genome-scale one
are inverted and both correct (`TRAPS: toys-rank-hotspots-backwards`).

### u-ruler-arm
`priority: next · kind: measurement · stamped: 2026-08-2x`
**Build the `U` ruler arm — the uniform-gDNA null ruler priced end to end.** The oracle ruler is NOT the
ceiling: a perfect-composition ruler is a LOSS by ~2× on the two capture-OFF in-scope strata, because
the correct factor there is exactly 1.000 and the oracle ruler sits FURTHER from it than the shipped one
— `oracle_ruler` is an A/B between two wrong rulers. The null `U` (the oracle's gDNA total laid down at
exactly uniform density) reads essentially 1.000 with no fitting; `calibration_vs_oracle.py` already
carries the column. Admissible at capture-OFF only; read `ruler_n_moved`, never the aggregate — the
aggregate barely moves while most transcripts are redistributed.

### g00-shrinkage-upstream-repair
`priority: next · kind: defect · stamped: 2026-08-2x`
At the zero-gDNA control the shipped effective-length shrinkage contracts the large majority of
transcripts where the correct factor is exactly 1.000, capture-OFF included, with `rho_ref` fabricated
from false-positive gDNA. ⛔ It is a SYMPTOM: feed the SHIPPED function correct composition arrays and it
returns the correct factor at `g00` and on the deferred stratum alike — **repair the COMPOSITION it
reads, not the function**. `priors.py` imports `_global_reference_density` from `capture_eff_length.py`,
so one repair serves both consumers. `calibration_vs_oracle.py` is the only instrument whose patch point
is upstream of the ruler, which is why nothing else has ever priced this.

### capture-blind-gdna-divisor
`priority: next · kind: defect · stamped: 2026-08-31`
**The gDNA opportunity divisor is blind to the probe panel** — with `eb-shrinkage-magic-ess` it owns the
already-priced −5.90 % capture-ON length ceiling. `gdna_opportunity_from_index` is computed from the
INDEX alone, so under capture it removes only ~6 bp of a ~30 bp length selection. This is what
`fl_anchor_gap.py`'s `G-gdna` control has reported as "impossible, therefore a bug" (+6.0 % on all six
capture-ON rows) since 2026-08-17: gDNA has no introns to miss, so a moving control is a divisor error.
`capture_eff_length` already models the probe panel and the divisor should use it. It is also what makes
the crossing pools' weights unestimable under capture (estimated `a_2`/`a_3` 0.002/0.002 vs truth
0.980/1.000), which blocks `ISSUES: crossing-pool-contrast`.

### eb-shrinkage-magic-ess
`priority: next · kind: defect · stamped: 2026-08-31`
`POOL_EB_PRIOR_ESS = 1000.0` shrinks the gDNA pmf toward `global_pmf` — a mixture that is mostly RNA
whenever gDNA is a minority — at a magic ESS whose own comment says to revisit it. Measured INERT on the
ladder (0.01 bp; pools 1.5–3.8 M) and DOMINANT on the fl-gap arm at `g05` capture-ON where the gDNA
pools collapse (`ship−pool` −23.7 of −31.7 bp). Same shape as the deleted `Beta(14,14)`: a
constant-weighted pull toward a contaminated anchor. ⭐ The derived replacement is the strand fit's own
move — reconcile the pools against EACH OTHER by their precision (`EQUATIONS.md` §6c), never toward a
mixture.

### refit-vs-message-arbitration
`priority: next · kind: design · stamped: 2026-08-2x`
The unstranded × capture-OFF exon cell is a REFIT-versus-MESSAGE arbitration question, smaller than
recorded ("+21 %" was one rung of four; the stratum is a small net win). Not a reference problem — a
`solvable_exon` never received a reference location. What solves those exons is the REFITTED gDNA prior;
the message is the accurate voice there and is displaced by the refit — two imputations at one slot with
nothing arbitrating them. Re-derive with `pass0_claimed_ab.py` (its STAGE reading prints
`calib_refit_iters` beside the shipped count). The arbitration belongs with the landscape prior
(`ISSUES: measured-prior-rung-4`), where the refit's own status is decided.

### prior-fidelity-vs-deliverable
`priority: next · kind: question · stamped: 2026-08-2x`
Why is prior fidelity anti-correlated with deliverable quality? Leading answer: it is not the prior, it
is the messages — at the worst slots the self-solve with the fitted prior is nearly correct and the
message layer destroys it, with the certified-RNA channel alone recovering truth at most of them. ⛔ That
dissection ran at a retired rung; confirm on a second stratum before closing. ⛔ Exclude first: every
prior-injection arm left the effective-length shrinkage unsubstituted, so part of "a better prior does
not show up" may be the RULER (`ISSUES: u-ruler-arm`, `ISSUES: g00-shrinkage-upstream-repair`).
`TRAPS: the-intermediate-is-not-the-deliverable`.

### landscape-trains-on-real-substrate
`priority: next · kind: measurement · stamped: 2026-08-2x`
Confirm `_fit_gdna_hyperprior` then trains on data rather than on the prior: it trains on
`belief.f_g·mass` over expressed REGIONs INCLUDING exons, so too strong a reference means it learns the
prior back. Measure the no-evidence share of training mass before and after with
`composition_evidence_census.py`; do not accept a change that moves only it. Nothing to build — this is
`measured-prior-rung-4`'s payoff check.

### expand-the-gdna-spectrum
`priority: later · kind: decision · stamped: 2026-08-2x`
The owner wants the gDNA spectrum filled (1, 5, 10, 25 % up past 90) without multiplying into dozens of
benchmarks. Levels are informative where behaviour CHANGES — each new level must be justified by a
measured transition and should cross a REDUCED set of the other axes until an interaction is shown. The
cost per condition is real (simulate, scan cache, oracle cache, certification); pay it when the current
scenarios are exhausted, which the dissection loop establishes. ⚠ See
`ISSUES: flgap-panels-stale-nascent-model` before comparing across panels.

### flgap-panels-stale-nascent-model
`priority: later · kind: decision · stamped: 2026-08-22`
The two fl-gap side panels were NOT regenerated in the sparse-nascent rebuild and still carry the
RETIRED uniform nascent model. Each is internally consistent with its own config, so both remain valid
on their own terms — but a claim spanning the ladder and a side panel varies TWO things and is not
comparable. Re-simulating them is the owner's call, not a prerequisite anyone should assert.

### psi-lambda-bracket-unshipped
`priority: later · kind: decision · stamped: 2026-08-2x`
ψ's λ bracket was too narrow to express its own prior; `DensityLandscape.required_logodds_window`
derives the correct bracket with no chosen constant (predicted matches measured on every stratum).
Built, gated, priced: nearly every in-scope condition improves, one dense capture-ON rung regresses
marginally, `g05` improves on both strand settings. **Ships OFF** pending two unpriced costs: memory at
genome scale (a small multiple on `sweep_n_grid`) and the end-to-end thermometer. Re-derive with
`ladder_arm_ab.py` / `arm_sweep.py`.

### alt-splice-rung-unverified
`priority: later · kind: question · stamped: 2026-08-2x`
Do we solve ALTERNATIVE SPLICING correctly? The `alt_splice` toy rung exists and is unverified
(`toy_harness.py --list`). Cheap, and the only structure where several splice junctions share a
BOUNDARY — in scope on all three shipping strata.

### transfer-variance-premise
`priority: later · kind: question · stamped: 2026-08-2x`
Does the message transfer variance correctly price a ratio built on a handful of counts? — PARTLY
answered by `transfer_variance_audit.py`: the shipped `transfer_logvar` is a counting term plus a
composition term, so every term shrinks as either slot deepens and a deeply-counted transport arrives
essentially undamped. The obvious substitute is refuted in both directions (the landscape's per-slot
posterior is TIGHTER than counting variance; its population spread over-states a fitted premise ~10×).
What survives is a per-hop-type premise; the honest seed is that under capture the posterior LOOSENS
(mode-membership ambiguity). `EQUATIONS.md` §3.5d.

### nascent-stress-sensitivity
`priority: later · kind: question · stamped: 2026-08-22`
Does any in-scope verdict depend on the nascent STRESS level? The ladder runs `on_fraction 0.50`
(development stress); realistic is ~0.10 (`DESIGN.md` §0b). No new panel needed: re-simulate the single
worst in-scope scenario at the realistic level and check whether any RANK moves. A verdict that only
holds at stress is a robustness finding and must be labelled as one.

### relay-od-r-discontinuity
`priority: later · kind: defect · stamped: 2026-08-30`
The relay is discontinuous in `od_r` at ~1e−5 (`g98 ss0.50 capture-OFF`: error 217,531 at
`od_r ≤ 1e−7`, 212,581 at 1e−5 — a threshold in the relay/anchor path, not a response;
`TRAPS: a-constant-parked-a-value-off-a-knife-edge`). In the ANCHOR, not the strand estimator.
⭐ BOUNDED: a 1e−5 nudge across all 30 test conditions moves 1/30 rows more than 0.5 % per policy (worst
1.65 %), so do not believe a single-row policy difference below ~2 %, and a relay comparison crossing
the edge is not attributable.

### f32-strand-tilt-at-half
`priority: later · kind: defect · stamped: 2026-08-2x`
At κ = ½ the strand mean is ½ identically, but the AMBIG cube evaluates the sum in float32 and departs
τ-dependently — a manufactured tilt growing linearly in depth. Falsified decisively: the same cube with
only the strand term in float64 returns `w_pos → 0.500000000` at every depth. Negligible at panel scale
and pre-existing; bites only a very deep library with a dominant locus. ⛔ Repair the strand term inside
`_solve_ambig_logodds`, not the cube (its f32 storage is an authorised memory choice);
`TRAPS: panel-before-src`.

### hygiene-ledger
`priority: later · kind: hygiene · stamped: 2026-08-31`
Each its own commit (`TRAPS: one-thing-varied`), none moving the 0.8.0 metric:
**the wave-3 frame migration** — the standalone bank-readers still on pass one
(`structural_claims_audit`, `held_flux_census`, `gdna_pool_census`, `landscape_head_to_head`,
`abundance_landscape_census`, `transport_dispersion`, `fl_pool_purity`, `certified_q_census`,
`anchor_opportunity_census`, `calibration_truth_ab`), migrate as touched per the frame ruling
(`DESIGN.md` §4.3; a claim spanning frames must say so); the index's duplicate map
(an ALIAS MAP `dropped_t_id → kept_t_id`, not re-admission); `mass_*_boundary` → `count_*_boundary`
(crossing INCIDENCES); restore the moment tests deleted with the length channel; the ledger of dead
surface and stale sibling references; the stale comment at `pipeline.py:416` (the drain's fl models are
pass one's, not "the SAME pool the calibrator reads" — line 382 has it right). ⚠ The duplicate map needs
an index rebuild but not a panel re-scan — verify with `rescan_panels.py`; `reach` is covered by no
other hash.

### oracle-effective-length-diagnostic
`priority: later · kind: measurement · stamped: 2026-08-2x`
Started, not finished (~half an hour); re-ranks the two ruler issues. ⛔ Needs a stashed pre-closure arm
or it measures nothing — one arm is not an A/B. Demoted because it ranks post-calibration items while
calibration's own path ranks higher.

### the-cancelling-pair
`priority: parked (refused twice) · kind: design · stamped: 2026-08-26`
`struct_lock` rescoped to `g1_locked ∧ REGION` AND the `intergenic|exon` boundary claiming its
RNA-contaminated crossing mass as gDNA — neither half has an honest price alone
(`TRAPS: a-cancelling-defect-pair`); five xfails go green iff the pair lands; priceable only with
`--messages on` (`TRAPS: an-ablation-that-never-ran`). ⛔ RE-PRICED 2026-08-26 with the measured intron
reference as replacement load (`ladder_arm_ab --arm stage1_pair{,_onesided}`) — **STILL REFUSED**:
marginal to the reference alone it worsens two of three in-scope strata; wins confined to `g00`.
**The analysis, kept whole:** the certified-RNA channel is a LOWER BOUND delivered two-sided; making it
one-sided (`−½·p·max(0, mo − log f)²`, no new constant) is the only mechanism the zero-gDNA control has
ever endorsed on every row, yet the panel reads a small regression — the two-sided term was doing two
jobs (the bound, and by its upward side a de-facto gDNA LEVEL channel), and removing the accident
exposes the real gap on exactly the stratum with no working level channel. ⛔⛔ **The level channel is
structurally disconnected — a theorem**: the chain is strictly REGION·BOUNDARY·REGION·…, so a
one-slot-step channel is BIPARTITE; the only licensed originators of a gDNA level are structurally
pure-gDNA REGIONs, so no REGION can ever receive one. Measured at `g00`: thousands of BOUNDARY
receivers, exactly zero of every REGION class. Two repairs already refused: a BOUNDARY originating a
level (patches a symptom); a two-step level (dominated by intergenic anchors — hands every exon the
off-probe floor). New fact: gene-edge BOUNDARIES carry EXACTLY ZERO RNA fragments on 32/32 conditions
(`structural_claims_audit.py`), so the mis-scoped mask's load is zero-count slots emitting certainty.
Run any revival against `--messages on` at `g05 ss0.50 capture_on`
(`TRAPS: all-small-singly-large-jointly`).

### crossing-pool-contrast
`priority: parked (blocked) · kind: question · stamped: 2026-08-31`
Should the gDNA length model run a second contrast on the CROSSING pools? Mirrors the contained pair
(mature RNA cannot cross an intron|exon boundary). Measured with ORACLE weights on both fl-gap sign
arms: under capture the crossing route BEATS the shipped contained one by a lot (TV 0.076 vs 0.136,
+2.8 bp vs −9.8; other arm 0.078 vs 0.182, +1.6 vs −23.9); OFF capture much worse (pool 3 starved,
29–630 fragments). ⛔ Blocked twice: ① the weight estimator does not transfer under capture
(`ISSUES: capture-blind-gdna-divisor`); ② the substrate cannot test the premise — every shadow
transcript lives on `test_blank` with NO intergenic|exon boundaries, so pool 3 measures exactly 1.0000
pure and the contrast silently reduces to pool 3 alone (`TRAPS: purity-is-a-property-of-the-annotation`).
⭐ Reshaped 2026-08-31: per-pool enrichment over the off-target rate is ~1.0 for both CONTAINED pools
even under capture and 294–338× for both CROSSING pools, and the two crossing enrichments are EQUAL
(ratio 1.002/1.034), so the weight RATIO `a_2/a_3` survives capture from the index alone — only the
common LEVEL (one scalar) is missing. **Would answer it**: a shadow transcript overlapping an annotated
gene edge on `test_chr`; the missing scalar via a one-sided fit within the probed-boundary stratum or a
§6c-style reconciliation. ⛔ Non-negativity identification is already refuted for the contained pair.

### capture-degeneracy-standing-risk
`priority: parked (watch) · kind: risk · stamped: 2026-08-31`
The gDNA two-pool contrast survives capture by a DEGENERACY, not its premise: under capture the
shared-contaminant assumption is false (TV 0.95 vs 0.06–0.14 off), and it is safe only because the
intergenic pool is depleted-not-impure, so `a_0` clips to 1 and the algebra collapses to `g = f_0`
(verified to 3e-17). **A probe panel that put RNA back into intergenic space would break it silently.**
`_deconvolved_gdna_counts` carries the derivation; nothing on the current panels can fire this.

### pure-rna-mirror-asymmetry
`priority: parked (tracked, not blocking) · kind: defect · stamped: 2026-08-2x`
Two exact per-fragment mirrors of a PURE-RNA library deconvolve differently in `mass_gdna_region` by a
few percent, neither boundary-only nor monotone in strandedness — a few percent of the zero-gDNA
false-positive channel. An R1-sense library is now simulable, so it is measurable whenever wanted.

### parked-capture-pilot-sign
`priority: parked · kind: question · stamped: 2026-08-13`
The two capture-ON pilot rows that disagreed about the SIGN of every length correction. Both panels were
deleted 2026-08-13; the correction concerned is inside the retired length channel, so it cannot change a
0.8.0 decision. ⛔ If revisited: do not average the two rows — find which one is lying. The fl-gap
panels are not a drop-in replacement (retired nascent model).

---

## CLOSED / REFUSED — ⛔ do not rebuild these; append-only

⛔⛔ **These entries point FORWARD — *do not rebuild* — so they are never deleted, and every row keeps
its stamped measurement exactly as recorded** (the deliberate exception to the numbers policy: a
graveyard row without its number is an invitation).
⚠ Where a mechanism's only target was **unstranded × capture-ON**, that target is DEFERRED, so the row
is moot as a 0.8.0 candidate on top of being refused; what is never moot is the `g00` zero-control
column. ⚠ **PANEL STAMP**: a row measured on "all 36 conditions" or quoting `g01`/`g10`/`g25`/`g75`/
`g90` predates the ladder retired 2026-08-13; the verdict stands as a record — re-opening one means
re-running it on the current panel. ⚠ "the RNA fragment-length model" row below is the accumulator's FL
*geometry* (ships in 0.8.0); the length-channel retirement is of a CALIBRATION COMPOSITION channel.

| | closed by | verdict |
|---|---|---|
| **the gDNA scale rule** · **the mass pin** · **TSS/TES as the population licence** | landed 2026-08-04 | ✅ `EQUATIONS.md` §3.5/§3.5b/§3.5c, gates in `test_gdna_scale_rule.py`, `test_relay_mass_rescale.py`, `test_terminus_population_licence.py`. ⚠ The ceiling says the mass pin cost the panel **nothing** (+0.0002 to delete it outright); it landed on the derivation and on being free |
| **face (I) of the `intron\|exon` BOUNDARY** | re-solve ceiling + panel arm | ⛔ **DO NOT BUILD.** The derivation (`EQUATIONS.md` §3.6) is re-verified and is not what failed: handing both BOUNDARIES the ORACLE truth and re-solving is worth **−0.000** off capture, and the ladder prototype is **negative** (mwae 0.0413 → 0.0426, confidently-wrong +10.7 %). TRAPS: panel-before-src |
| **a LEVEL transfer from the intron** | toy + panel | ⛔ **REFUTED**, +0.207 on capture-ON × unstranded — capture inverts which side is well-counted (TRAPS: capture-inverts-the-counted-side) |
| **the RNA fragment-length model** | `length_ceiling.py`, one pmf at a time | ⛔ **−0.02 %** at pass-0, **+0.21 % (worse)** over all objects. Root cause exact (`pi(w)` scores sj *crossing*, the pool requires the splice to be *seen*). ⭐ Its value is the BOUND: the whole fragment-length-model cluster costs ≤0.43 % of the shipped solve. TRAPS: price-the-halves-separately |
| ⛔⛔ **the-rna-length-law-fix — return the contrast's `r_hat`, model sj observability, or improve `rna_pmf` at all** | step-0 re-measure + a mean-shift dose–response + `calibration_vs_oracle.py`, sparse ladder, reseed floor re-recorded same session, 2026-08-31 | ⛔ **REFUTED — `rna_pmf` is SOUND AS SHIPPED and no fix is licensed.** ① The −4 bp sj-observability diagnosis was measured on the **UNDRAINED** pool against undrained deposits (reproduced to the last digit: 218.35→204.08 vs 208.18); production fits the **DRAINED** payload, where the residual vs mature truth is **−0.24/−0.11 bp** off capture and **−1.06 bp** under it — the drain re-includes the gap-hidden spliced fragments (empirical selection ≈1.00 off / ~0.97 tail on). One name over two populations, a fourth instance. Dividing by `pi·o(w)` on the real pool OVERCORRECTS **+17 bp**. ② The consumer split re-measured (`g05 ss.99 ON`, floor 1,545): scorer +273 / drain −18 / EM-geometry +586 / **`calibrate(rna_fl_pmf)` −17,754** — but the truth pmf is NOT that lever's optimum: interp-shifting the SHIPPED shape +1.31 bp wins **−43,648**, +10 bp wins **−171,342** (and **−42,858** of 542 k at `g50 ss.99 OFF`), monotone far past truth while library `gdna_frac_est` and gene fp_mass degrade — so the transcript win is **COMPENSATION**, dissected to within-gene isoform reallocation at short-median-exon (≤150 bp) multi-exon genes: the standing EM assignment error (`ISSUES: per-transcript-prior-lane`) nudged through the per-locus prior, and a caution against ranking any calibration input on the transcript thermometer. ③ On the 0.8.0 metric the true RNA pmf moves misplaced mass **±0.3–2.9 %, mixed sign** — nothing. ⭐ What SURVIVES: `r_hat` (the contrast's discarded contaminant law) is accurate off capture (**−0.04/+0.19 bp** vs nascent truth) and broken under it (−11/−20 bp) — a fact worth keeping, licensing nothing; and the mature-vs-all-RNA estimand mismatch is real but ≤1.7 bp and unpriceable above the floor |
| **TRAPS: pure-and-length-censored's κ residue, as an ACCURACY fix** | κ injected at exactly ½, all 36 conditions | ⛔ **−0.2 %** unstranded, worse on the shipped solve. ⭐ But the *general* defect — a boolean licence flipped by a small residue — is **the-capture-level-residual**, and the destruction control taught TRAPS: honesty-metrics-reward-ignorance |
| **a nascent-bearing ladder condition** | toy, 36 conditions × 7 rungs | ⚠ **−5 %**, and the wrong way on one stratum. Keep it as a harness arm (`--nrna 60`); it no longer justifies re-simulating the panel |
| **the gDNA prior's BIMODAL CAPACITY, and "give the prior more signal"** | a read of `gdna_landscape.py` + the production refit on real conditions | ⛔ **BOTH BRANCHES CLOSED.** The prior already renders the landscape correctly — **2.98 decades** of mode separation at `g75 ss0.99 capture_ON`, 30× more enriched mass ON than OFF, a single pile at the wall at `g00`. And a prior fitted from ORACLE truth is the same prior (0.04 dec). Not capacity, not signal, not location. ⭐ Why an evidence-free object cannot reach the vertex at all — and why that is the value of missing information rather than headroom — is `EQUATIONS.md` §9a |
| **the Jeffreys MEAN density location** | `--arm eta`, the `g00` zero control | ⛔ **REFUTED at +96,299 %.** It cannot say ZERO (`region_init.rho_g` is an exact 0 at 60,544/70,176 slots — the statement earning the −98 % at `g00`), and the TRAPS: a-ratio-cannot-carry-zero benefit it was credited with belongs to that fix, not to this arm. ⭐ If revisited the derived form is the Gamma **MODE** `max(a−½,0)/E`, which is exactly 0 at a zero count |
| ⛔⛔ **refused-transcript-weights — a SPECIFIC per-transcript allocation RULE: soft-min over exclusive objects with a per-object Jeffreys half** | built end to end and A/B'd on all 36 conditions, seed pinned | ⛔ **REFUSED — worse on EVERY stratum and on the zero control**, transcript Σ\|err\| **57.5 M → 81.6 M (1.42×)**, and a length-proportional variant **2.10×**. ⭐ **The MECHANISM is not what failed** — the gDNA:RNA split moved **+0.2 %**, exactly as the conservation identity requires, so the A/B priced the ALLOCATION alone. Three defects, each measured: exclusivity hard-zeroed **38.7 %** of transcripts; estimating a density on a tiny exclusive region and extrapolating over the whole transcript amplified variance up to **6,534×** (44.6 % of weighted transcripts had their density from <200 bp); and a per-object `+½` revived the silent half of the annotation (`frac_expressed: 0.5`), taking false-positive mass **18.6 M → 41.6 M**. ⛔ The rule and its config flag were DELETED; the LANE it rode on was kept (`ISSUES: per-transcript-prior-lane`). ⚠ The Jeffreys-mean half is its own row above — see `TRAPS: a-trap-names-the-defect-not-the-repair` |
| ⛔⛔ **refused-soft-min-path-weighting — the owner's own theorem, built faithfully and REFUSED** | 12 arms (4 modes × 3 multipliers) on `g00 ss0.99 capture_off`, 3 re-run on the blind stratum `g50 ss0.50 capture_on`, base re-recorded in the same session | ⛔ **REFUSED. Worse than `base` at TRANSCRIPT level on every rung of every arm and on both strata** — 1.317–1.604× at `g00`, 1.262–1.331× on the blind stratum — with transcript false-positive mass 1.76–2.20× worse. ⭐ **The one encouraging number does not survive:** at `g00` the same weights took GENE error to **0.395–0.527×** (against `oracle_alloc`'s 0.128×), and on the blind stratum that collapses to **1.006–1.041×**, i.e. nothing. ⭐⭐ **The MECHANISM is `TRAPS: an-upper-bound-is-not-an-estimate`, structural rather than tuning:** the theorem bounds a transcript by the thinnest object on its path, but **3,644 of 4,839 silent transcripts (75.3 %) share an object with an expressed one** and inherit its bound. The zero-weight SET was byte-identical across all twelve arms — a bound is zero only when every object is, a property of the data. ⭐ Retreating to GENE granularity did NOT rescue it (1.340× vs 1.317×), so the damage was never the within-gene split. ⚠ **Two things ESTABLISHED to keep:** the dial is monotone in the theorem's favour on both axes and both strata (`min` < `harmonic` < `geometric` < `arithmetic`), so the pooled `Σmass/Σopportunity` control is the WORST rung and the soft min does real work; and **0.0 % of expressed transcripts were ever zeroed**. `scripts/design/transcript_weights.py`, `tests/calibration/test_transcript_weights.py` (31 gates) |
| **a threshold anywhere in the licence family** | TRAPS: a-threshold-on-a-fitted-residue implemented and refuted one | ⛔ τ is continuous across the region, so any floor is a tuned constant (TRAPS: a-threshold-on-a-fitted-residue, TRAPS: a-licence-with-no-floor, TRAPS: a-multiplication-gated-by-a-trace — refused three times) |
| **simulator captures a pre-mRNA through every probe its genomic blocks span** | NASCENT SCOPE RULING, 2026-08-22 | ✅ **WON'T-FIX** — nascent × capture fidelity is out of scope (`DESIGN.md` §0b); the sparse rebuild shrank the residual it explains (capture depletes nascent ~an order of magnitude), and no verdict depends on it (`TRAPS: the-panel-enriches-nascent-by-its-own-probes`). Re-open only if a real library forces it |

### CLOSED / REFUSED — the reference-mean family (2026-08-15/16)

⭐ Kept separate because the table above is rules for resolving DOUBT at an evidence-free slot; these are
attempts to give ψ's reference a MEAN — plus two out of SHIPPING it — each built, measured on the panel
with both zero controls, and refused. The form that survived is `DESIGN.md` §6b.1.

⛔⛔ **giving `τ_λ` the location term's curvature** ("the asymmetry with the intron factory"). Built,
measured, refused on three counts, and the motivating reasoning was wrong at every step: the 3,227× fall
in `τ_λ` at a pinned slot is **~98 % the `[f(1−f)]²` Jacobian** (nothing was lost); the contribution is
a **boolean gate flip** releasing the full COUNT precision, not a ¾-unit increment (τ = 0.029 and
τ = 1e6 both return 850.44 of a 850.50 ceiling); and it carries no count, so it credits **data-free**
slots (`n = 0` ⇒ `prec_g` 0 → 0.2026) — the very population the structural reference's safety argument
rests on being empty. Measured: bit-identical on the deliverable on all 32 panel rows, moving only
`has_own_composition_evidence`. `TRAPS: a-priors-curvature-is-not-the-datas-information`.

⛔ **softening the prior to a per-object one-pseudo-fragment floor** (`m_i = E[g]_i/(E[g]_i+1)`). Worse
on **every** stratum (0.609 / 1.045 / 0.580 / 1.000 against 0.381 / 0.659 / 0.363 / 0.800): on a
structurally pure-gDNA object the truth IS `f_g = 1`, and a soft floor pulls it off the vertex. What
replaced it introduces no constant — the lattice's own top point `σ(L)` (`EQUATIONS.md` §9c.1).

| # | mechanism | why it was refused |
|---|---|---|
| 1 | **a fitted RNA density `logP_r`**, the mirror of the gDNA landscape | ⛔ the only non-circular form (fit from the solver's own belief) reads **0.988 / 0.997 / 1.037** — nothing, then worse: feeding ψ a density fitted from ψ's own belief tells it what it already believes. ⭐ And the ORACLE version's gain did not survive a shuffle — a shape that is wrong on purpose BEAT the true one at `g98` (0.786 vs 0.854), so the attribution was never established (`TRAPS: attribution-must-survive-a-shuffle`) |
| 2 | **a library-wide Beta mean, `a = f_lib`** | ⛔ `g05` regresses **1.43×**; `f_lib` is calibration's own output so the loop has positive feedback with both vertices attracting; and moving `a`/`b` sets the TAILS as well as the location, so `b = 0.03` leaves **57 %** of the prior outside `L = 10` |
| 3 | **the OBJECT-weighted mean instead of `f_lib`** | ⛔ **the two split by STRAND and neither wins everywhere**: object-weighted 0.584 / 0.452 on the two stranded strata but **5.570×** on unstranded × capture-OFF. ⚠ The sweep that motivated it was `ss_0.99` on three of its four rows — `TRAPS: never-pool-the-strata`, met on a sweep rather than a panel |
| 4 | **a stratified ASSERTION** — pure-gDNA strata claim `f_g = 1`, reweighted by stratum size | ⛔ reads 1.000 at every condition **including `g00`**: an assertion cannot see a library with no gDNA in it. ⭐ What replaced it is a per-object DENSITY, which needs no reweighting at all — the strata select the training set, not the answer |
| 5 | **a pooled RNA density from sj flux** | ⛔ RNA spans six decades with no genomic autocorrelation, so a pooled flux is not a population parameter (owner). It scored well only by sitting on the mass-weighted centre — `TRAPS: a-mean-hits-the-mass-weighted-centre-by-luck`. ⭐ Replaced by RNA-as-residual, which predicts no RNA at all |

⚠ **One more that is a CAUTION rather than a refusal:** `f_g ≤ 1 − S/M` as an assumption-free bound from
certified RNA. `boundary_spliced` is a SEPARATE bank from `boundary_unspliced`, not a subset, so the bound
is simply false — the truth violated it by **302**. The correct statement is the identity
`ρ_r·E_r = unspliced_RNA + S`, i.e. S SUBTRACTS.

### the-truncation-free-region-bank — BUILT, PRICED, REFUSED (2026-08-20). Do not rebuild it as a drop-in.

The shipped REGION bank reads `rho·P(w<=ell)` (`TRAPS: a-cancellation-is-conditional-on-its-support`),
and the truncation-free repair — `region_start_count / ell`, the STARTS-IN relation, no fragment length
in the weight — was built behind `CalibrationConfig.region_abundance_bank`, gated by an absolute fill
gate fired three ways, proven inert outside the currency arm (16/16), and priced on all 16 conditions.
⛔ **REFUSED: it improves the currency arm's pass-0 exon solve (no-evidence ≤0.05 mass coverage
45.8 % → 66.1 % at `g98 ss0.99 ON`) and regresses its DELIVERABLE — the four zero controls 2.18×
(each 2.0–3.3×), the deferred stratum 1.84×, stranded × ON 1.20×; the only stratum win is 0.843× on
stranded × OFF.** ⭐ The named suspect (one leg measured): under the shipped bank **58.6 % of live exon
hops carry a REGION bank of exactly 0** and `enrichment_ratio` returns its 1.0 default — an ACCIDENTAL
MUTE the currency arm's standing partly rested on; the live bank un-mutes them and the policy's
machinery converts a better input into a worse answer
(`TRAPS: the-intermediate-is-not-the-deliverable`, measured inside one arm). ⭐ What survives: the
truncation algebra (`EQUATIONS.md` §2), the fill gate, the two-names rename, and the wall-exposure
numbers (FLUSH 3.7–9.1 % / BINDING 7.4–22.6 % of exonic starts by template population, spliced
coordinates, MAX collapse; the full table stays in the dev sandbox until a consumer lands).
⛔ The knob was DELETED after pricing (converge-and-delete); the full implementation is one commit
before the deletion (`a2b81b34`). Re-opening this requires a policy whose DELIVERABLE improves under a
better level channel — none exists today.

### the-doubt-graveyard — ELEVEN MECHANISMS PRICED, ELEVEN REFUSED. Do not rebuild these.

⭐ Promoted from a working doc when it was deleted (2026-08-07). Every row is a real build that was
measured and refused. ⛔ **`g00` is the owner-required ZERO-gDNA control: its truth is exactly 0, so
every fragment there is a false positive with nothing to cancel it** — which is why a mechanism can look
good on its target and still be inadmissible.

⛔⛔ **THE TABLE IS UNCHANGED BY THE 0.8.0 SCOPE AND IS NOT TO BE EDITED — it is eleven measurements.**
⚠ How to read it now: the **`its target`** column is, for most rows, a win on **unstranded × capture-ON**,
which is the **DEFERRED** stratum — so those wins are moot as 0.8.0 arguments and were never admissible
anyway. ⭐ The **`g00`** column is the one that decided each row, it is in scope for all four strata, and
it is why the pattern paragraph below is the real result. ⚠ `zc_struct_lock_g1` is the one row still live
(`ISSUES: the-cancelling-pair`) — as **half of a pair**, never alone, exactly as its row says. ⚠ Its `g00`
column was RE-PRICED 2026-08-18 on the licence + SPLICE IN-fix relay: **0.71×** on top of the SPLICE IN
fix, with the four `nrna_none` zero controls 1.8–15× worse — the sign on the total flipped, the verdict
(half of a pair) did not; the table row is the 2026-08-11 measurement and stays.

| candidate | what it did | `g00` | its target | why it died |
|---|---|---|---|---|
| `zc_jeffreys_mean` | `ρ_g = ½/E_g` at zero mass | ⛔ +7,269 % | −13.9 % | moves the mode UP |
| `zc_logmean` | `ρ_g = e^{ψ₀(½)}/E_g` | ⛔ +6,264 % | −11.3 % | moves the mode UP |
| `zc_anchor_mute` | no `prec_g` at empty locked slots | ⛔ +5,554 % | −7.7 % | kills the zero-gDNA win |
| `zc_struct_lock_g1` | scope `struct_lock` to `g1_locked ∧ REGION` | ⛔ +3,207 % | −1.2 % | ⭐ the MIS-SCOPED mask is load-bearing |
| `zc_reference_var` | `Var(f_g) = ⅛` where `τ = 0` | ✅ +0.0 % | −0.3 % | ⭐ passes the control and is INERT |
| `zc_discrepancy` | `+½ log D` shift, `(log D)²/12` | ⛔ +982 % | panel +4.5 % | moves the mode UP |
| `zc_disc_var` | the variance alone, mode untouched | ⛔ +255 % | panel +0.9 % | damping cannot bite |
| `zc_ref_prior` | own belief = ψ's reference, `τ + 1/π²` | ⛔ +3,792 % | −14.9 % | moves the mode UP |
| `zc_ref_prior_damp` | the two above, PAIRED | ⛔ +3,809 % | −15.5 % | ditto |
| the `eta` rebuild | a clean frame-free re-derivation | ⛔ unbounded | +85–103 % | see `DESIGN.md` §6.1 |
| the mean-location as a structural floor | the same idea in the LEVEL channel | ⛔ +96,299 % | — | it cannot say ZERO |
| `struct_lock = g1_locked ∧ REGION` **re-priced** | the standing strict xfail's own fix, on top of the SPLICE IN precision repair | ⭐ **0.71×** | in-scope **+2.5 / +2.1 / +0.5 %**, deferred +17 % | ⛔ **2026-08-18**, and the sign on `g00` FLIPPED from the 2026-08-11 row above once the relay was fixed underneath it — but the verdict did not: the four `nrna_none` zero controls go **1.8–15× WORSE**. The empty exons' "gDNA = 0 @ 0.2026" is LOAD-BEARING at AMBIG `exon\|exon` boundaries in an RNA-only library (an RNA+ level claim with no RNA− claim drifts ψ to 0.38 without it). ⭐ Still half of a pair |
| the mass rescale refuses a ZERO-MASS source (`pinM`) | `may_share_composition` additionally requires `M[src] > 0` | ⭐⭐ **0.31×** (134 k — better than the pre-licence relay's 154 k) | unstr × OFF 0.999, str × OFF 1.004, ⛔ **str × ON 1.032, six of six worse** | ⛔ **2026-08-18.** Under capture the empty slots between probe-covered stretches are the CONDUITS a relayed composition travels through (`TRAPS: the-divergence-was-a-barrier`), so the rescale at those hops is load-bearing. ⭐ The sharp predicate is NEITHER this nor the row above: refuse a source's OWN zero-count artefact, KEEP a relayed composition passing through an empty slot |

⭐⭐ **THE PATTERN, AND IT IS THE REAL RESULT: every one of the eleven was a rule for how to resolve DOUBT,
and at `g00` the doubt must resolve to NO gDNA.** A rule that lifts an evidence-free slot off zero is
inadmissible there however well it scores elsewhere. ⛔ The only candidate the control has ever ENDORSED is
the one-sided certified-RNA bound (−81.9 %, 8/8), and that one is panel-negative alone — it is half of
`TRAPS: a-cancelling-defect-pair`; see `ISSUES: the-cancelling-pair`.

⚠ Four more `zc_*` arms exist as decomposition REVERTS used to attribute the 39 % win, not as proposals:
`zc_own_count`, `zc_live_count`, `zc_total_n` (inert) and `zc_transfer`, which reproduces the pre-fix tree.
