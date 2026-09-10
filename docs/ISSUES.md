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


### splice-out-premise-bias-uncorrected
`priority: later · kind: decision · stamped: 2026-09-02`
**The exon → intron|exon boundary message assumes spliced and unspliced fragments at one face share
capture affinity, and the assumption is MEASURED to fail as a BIAS under capture — recorded, not
corrected.** The two-witness estimator (`messages/transfer_rows.splice_out_row`'s premise: the exon's
and the boundary's own strand solves imply the face's spliced-to-unspliced ratio, the face measured
one) reads `log a ≈ 0` off capture, `+0.28` under benign capture (the enrichment shoulder — the
contiguous ratio `DESIGN.md` §6b.2 measured at 0.775, seen from the spliced side), `+0.78` on
junction probes, each with negligible spread. A fitted widening therefore does nothing (measured
identical to a few fragments). A bias correction would subtract a fitted level — the kind of
cross-locale fudge the rebuild exists to avoid — and is the owner's call; re-derive with the
estimator before deciding (the session harness `item1_proto.py`'s `PREMISE` arm carries it; promote
it if the decision is taken up). Cost of leaving it: the message over-claims gDNA at nascent-bearing
probed boundaries by ~`a` in the ratio, bounded by saturation at pure gDNA.

### gdna-landscape-trains-on-false-positives
`priority: now · kind: defect+question · stamped: 2026-09-02 (owner: "figure out the implications")`
**At zero-gDNA conditions the fitted gDNA landscape is trained ENTIRELY on pass-1 false
positives, and the refit loop then entrenches them** — the prior that exists to rescue blind
slots is taught by the blind slots. Measured (test chromosome, `landscape_poison_study` in the
session scratchpad; re-derive by spying `_fit_gdna_hyperprior`'s inputs against `slot_truth`):
at `g00 ss.50` **100 % of the training gDNA-mass** sits on slots whose certified gDNA is zero.
Two mechanism facts: the training gate (`fp ^ fn`) is an ANNOTATION test that cannot see κ = ½
(unstranded exons train at full membership), and `_reliability` consumes the SOLVED posterior's
`var_gdna` — small at blind slots — rather than the own-composition variance
(`own_composition_logvar`, ∞ at τ = 0) that would catch them; measured mean weight 0.90–0.93.
⛔ **The naive fix is REFUTED (2026-09-02): excluding κ-dead exons from training destroys the
bootstrap** — the blind slots' training values carry the MESSAGE ROWS' corrections after sweep
1, and re-fitting on them is how local true evidence generalizes population-wide (whole-library
`g50 ss.50 ON` regressed 2,691 → 56,422 under the exclusion; the stranded control was
byte-identical). So the open question is a robustness mechanism that discounts PURE-ECHO
training mass without starving message-corrected mass — candidates to derive, not assume:
sweep-aware weighting (sweep-1's blind values are echo; later sweeps' are message-informed),
own-OR-DELIVERED-evidence weighting, or capping the landscape's claim strength by its
evidence-bearing mass. ⚠ Downstream symptom already priced: `silent`'s residual `g00 ss.50 OFF`
whole-library error (~1.25 M on the ladder) is largely this loop; the transfer policy's rows
break it wherever they deliver. ⚠ A second, separate defect found by the same study: a
consumer reading `log_rho[-1] − log_rho[0]` as "the landscape's span" reads the GRID (built
from `mass/eff` — what is expressible), not the fitted mass's support.
⭐ **A second exposure (item 2 of the message rungs, 2026-09-02).** On the ladder at `g05 ss.99 ON`
the boundary → intron message improves 5,144 near-empty introns (7,300 fragments in all) by 26 % at
their own slots, and the refit prior then moves the thin boundary class — which receives nothing —
by +1,104 (1.9 %), for a whole-library 1.005× against silence. The two were told apart only by
scoring the message at its DESTINATIONS (node-local) beside the whole-library number; the reversed
row moved the same boundary class by a different amount in the same direction, so the boundary
movement is the prior's response to changed training beliefs, not the message's sign. Same regime
as rung 1's recorded +169. `DESIGN.md` §6b.5 carries the measurement.
⭐⭐ **A third exposure, and the re-ranking (item 5, 2026-09-02; owner: "an important priority").** The
by-class census (`policy_benchmark.py --by-class`) puts 62 % of the stranded capture-ON error and 85 %
of the zero controls' at exon|exon boundaries and the exons they wall — and item 5, the message built
for them, moved their error by ±1 % on the contaminated rows (they sit at the resolution of their own
evidence plus the prior, 1.8 % of mass) and −4.6 % at the `g00 ss.50 OFF` zero control, whose error IS
this prior's false-positive training (silent 327k → the messages 97k → nothing further). The
messages have done what messages can at these slots; the prior is the lever. ⛔ PARKED behind the
message policy's completion (owner ruling 2026-09-02, `DESIGN.md` §0c.0e): first in the parked list,
taken up the moment the checklist reads ✅.

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

### rename-row-and-face
`priority: later · kind: decision · stamped: 2026-09-09 (moved from the owner's sandbox note of 2026-08)`
Two more terms the owner flagged as obscure beside `drain` (`ISSUES: rename-the-drain`). **"row"**:
the calibration solver's word for a slot's max-normalised log-profile over the solve grid — a
message's composition claim, ψ's `lam_rows`, the intron factory's per-slot factor. The owner: "This
is obscure. It is not a biological term. It typically refers to rows of a table. We need to understand
and then rename this term." ⚠ It IS a row of a ``(n_slots, K)`` array, which is where the word came
from; a candidate is *profile*, which the transfer policy's docstrings already use for the same
object. **"face"**: the transfer policy's word for one DIRECTED side of a boundary — the pair
``(source, destination)`` a rule is keyed by, distinct from the boundary itself (one boundary carries
two faces). The owner: "a 'face' is a synonym for boundary, but we should try to keep the terminology
standardized." ⚠ Both are `rename_census.py --sense` passes with `rename_identity.py --freeze/--check`,
like the drain; the owner picks the words.

### drain-contaminates-certified-rna
`priority: later (PARKED by owner, 2026-09-01 — diminishing returns; the ceiling refused the in-solve correction) · kind: defect · stamped: 2026-08-31`
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
⭐⭐ **DERIVED 2026-08-31, and the measurements FORCE the design** (every prototype was gated on
byte-identity with production's own drain choices).
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
⑥ **The per-record reliability** (oracle truth walk, 100 % matched): the posterior is essentially
PERFECT at the extremes (~75 % of records) and UNDER-calls genomic in the torn middle — realized
null-truth 0.70 at q 0.50, same shape both regimes — so `Σ q_null` recovers ~70 % of the true
contamination; the mature-genomic-gap term is real (405 of 542 at `g50 ss.99 OFF`) and is RNA.
⑦ ⛔ **The NO-LEAK COUNTERFACTUAL CEILING REFUSES any in-solve correction**: rebuilding BOTH worlds
with the leaked choices flipped to genomic (residual leak exactly 0) moves the release metric by
**−0.60 %/−0.05 %** at the worst condition (`g98 ss.99 ON`) and ≤0.10 % mixed-sign elsewhere — all
under the ~2 % attribution floor. The leak's harm is the CERTAINTY CLAIM, not calibration accuracy.
⭐ **What remains to build** — both are recorded here, in full, and nowhere else:
① honest labeling — the drain records its own expected certified-channel impurity
(`Σ q_null`, variance free) into `DrainQC`/the QC report; no solve behaviour changes; and
② the one true false-positive REDUCER — repair the middle-bin posterior bias (suspect: the genomic
side's `_bottleneck` min over several noisy boundary densities vs one sj flux; a derived pooled
aggregation, no constant), which cuts actual leak draws with no trade-off. Each its own
DERIVE → prototype → A/B; the reliability curve is the re-check instrument.

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
`priority: now (its "reparameterisation away from λ" candidate IS `ISSUES: arcsine-magnitude-coordinate`) · kind: design-constraint · stamped: 2026-08-24 (owner + external review)`
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

⭐⭐ **FIRST MEASUREMENT, 2026-09-01 — no solver, certified `slot_truth` only, both substrates.** For
each slot, answer with its two chain neighbours' mass-weighted TRUE `f_g` (an ORACLE message: no
estimation error, no transport loss, no precision bug) and score it against the best CONSTANT answer;
`skill = 1 − err_neighbour/err_constant`, so **skill ≤ 0 means no policy built on neighbour transport
can win there, however well engineered**. Read the in-gene column (`R exon`, `R intron`,
`B exon|exon`, `B exon|intron` — where an unstranded library has no channel):

* **The ladder: 7/12 defined rows positive, median +0.295.** Strongly positive at mid-contamination
  (`g50` +0.53…+0.69) and at `g98` capture-OFF (+0.52); **NEGATIVE at every `g05` row**
  (−0.44…−0.93); ≈0 at `g98` capture-ON (−0.001, +0.068).
* **⛔ The anchored twin block: only 6/24 defined rows positive, median −1.635**, and catastrophically
  negative at low gDNA (`g05` −12…−28, `g25` −1.9…−4.6). ⚠ **The two substrates DISAGREE in sign**
  (`TRAPS: a-toy-and-a-panel-can-disagree-in-rank`) — the ladder is the shipping judgement.
* The pattern on both: neighbour information exists at MID-to-HIGH contamination and is actively
  MISLEADING at low contamination, where "almost all RNA" is already an excellent global answer that
  a noisy neighbour can only degrade. The four `g00` rows are UNDEFINED (constant truth ⇒ the best
  constant is exact), not zero.

⚠ **Read the baseline honestly**: the comparator is the best constant *computed from truth*, which is
an oracle-informed and therefore STRONG baseline — a real policy does not know it either. So a
negative row says "the neighbour is a worse predictor than a good global prior", not "messages are
worse than the shipped local solve". The decision-relevant contrast against `SilentPolicy` still
needs `policy_benchmark.py`. What the measurement DOES settle is the issue's literal question — the
information is present at `g50`/`g98`-OFF and absent-to-negative at `g05` and at `g98` capture-ON.

### the-lower-bound-noise-ratchet

**Priority: with the enrichment witness.** The level lane's residue at the LOW-gDNA end. A level made
from a node's own strand profile at an RNA-rich node is a measurement whose mode is noise around zero
(±0.003 in share at 2,000 crossings); its LOWER SIDE becomes a soft bound at its neighbours, and among
nodes of one density the tighter of the noisy neighbours wins — a ratchet of noise, not of evidence.
Measured (2026-09-05): the ladder's `g05 ss.99 OFF` 44,714 → 45,076 (+0.8 % after the intersection form
recovered two thirds of the product form's +2.4 %), and the rebuilt test chromosome's zero-gDNA control
`g00 ss.99 ON` 64 → 216, 187 of it at `capcluster_ab`'s nine inner terminus boundaries at share
0.007–0.011 with no gDNA at all. A TWO-SIDED own-profile level would average the noise away, and
measured −0.9 % / +1.5 % on `g50 ss.99 ON` / `g98 ss.50 OFF` against lower-only's −4.7 % / −0.5 % — so
the lower-only law keeps more fragments and this residue. What would remove it without giving up the
law: the same enrichment witness `two-sided-exon-row` waits for (where the library's gDNA is not
enriched, a level between nodes of one gene is two-sided). Until then, report it on every zero control.

### flux-floor-dispersion
`priority: with the transport-dispersion decomposition · kind: question · stamped: 2026-09-08 (re-stamped after the owner's ruling)`

The certified flux at an exon's junction is that strand's RNA level at the exon, an ESTIMATE of the
exon's abundance priced by the node pair (the junction's count at its rate against the exon's own
count of that strand per RNA opportunity, `count_price`; `DESIGN.md` §6b.13). The route rate scatters
around the exon's true body density BEYOND counting (stage 0, 2026-09-04: median −3 %, 5–9 % at depth;
nine readings on the block 0–40 % over), and the pair's price sees that scatter only where the exon's own
strand count disagrees with the rate — at a gDNA-rich exon it does (gDNA's half inflates the count and
the price widens, the conservative direction), at a pure-RNA exon it does not, and a lucky over-read
there is a sharp floor a few points too high. The floor is lower-sided (the two-sided estimate refused
at the cliff, §6b.13), so only over-reads cost. The instrument is `transport_dispersion.py`; the relay's
answer was a POOLED left-tail centre fit, refused with the relay's pooling.

### flux-price-witness-units
`priority: next · kind: problem · stamped: 2026-09-09`

The flux level's price (`DESIGN.md` §6b.13) compares the junction's route rate — in WHOLE-STRAND units — with
the exon's COLUMN count per RNA opportunity (`count_price(c_j, c_j / r_j, cnt[x, col_read], a_r[x])`). A
column holds ``(1 − κ)`` of the strand's RNA (plus gDNA's half and the other strand's leak), so at
κ = 0.31 every flux floor pays a systematic ``log(1 − κ)² ≈ 0.14`` nats² that is no disagreement, and on
unstranded data (κ = ½, the column is half of everything) ``0.48`` nats² — a 0.7-nat blur on every flux
floor of the half where a policy must WIN. The record: the golden scenario `strand_ss65_multi_iso`'s nested
exon (t1's exon inside t2's intron, no gDNA, 96 fragments) reads 0.152 gDNA through the pipeline (0.32 at
pass zero) under `transfer` against the relay's 0.000: its ceiling from the junction's flux reads
``f_g ≤ 0.38`` where the flux itself says ≤ 0 — 0.42 nats² of price, of which the κ term is a third and the
rest the rate's 12 % over-read of the exon's total. ⛔ The obvious repair is NOT a clean win: a prototype
with the exon's witness in the strand's units (the column split's asymmetry where the channel is live,
the total where it is dead; `fluxw_proto.py`, 2026-09-09) read WORSE at pass zero on the test chromosome
(`g05 ss.50 OFF` +12.6 %, `g05 ss.99 ON` +20 %; the full pipeline within 0.1 %), because the split's
asymmetry is ``R_s − R_s'`` and reads ZERO at an equal-abundance overlap exon — the both-stranded exons the
lanes were built for lost their flux floor to a counting-on-nothing price — and because a sharper flux
floor exposes the rate's own over-read (`flux-floor-dispersion`). What the price needs is a witness that
is the strand's RNA count at single-strand exons (the asymmetry, or the total less nothing) and a bounded
one at both-stranded exons (the asymmetry is a floor on ``R_s``, the column an upper bound), charged only
where the rate falls outside the bound — designed and A/B'd as its own step, halves apart. ⚠ The same
asymmetry witness on the HOP price (landed 2026-09-09, `_RnaLane.witness`) has the same blind spot: at
a node where both strands are lit it under-reads the weaker strand, so a dim claim about that strand
can arrive sharper than it should; on the ladder no such case moved a row (the dark host intron beside a
lit antisense exon is the case that occurs, and there the claim is true).

### ambig-node-as-a-gdna-source
`priority: after phase 2 · kind: decision, measured once · stamped: 2026-09-08`

The owner ruled that a determined gDNA level may propagate (2026-09-08). Built as: a both-stranded
node's own strand counts on the ``(λ, θ)`` cube plus the RNA levels it holds from its far side and its
own flux estimate, marginalised over the tilt, read as its gDNA level through its total, emitted like any
measured level (lower side across a face, priced by the totals' disagreement). Gated right on a
hand-built node (the mode at the bracketed density; nothing emitted without an RNA profile; no echo).
⛔ REFUSED AS BUILT on the chain: +2.6 % on `g05 ss.70 OFF` on all three test panels (11,247 → 11,534),
up to +13.7 % on the sparse panel's `g05 ss.99 ON`, +4…+38 % on the weak-κ zero controls; on the LADDER
`g05 ss.50 OFF` 1.745×, `g05 ss.99 OFF` 1.619×, `g05 ss.99 ON` 1.465× through the pipeline. At κ = 0.7 and
low gDNA the bracket has little leverage, so the emitted level's mode is noise, and a noisy level travels
as a floor — `the-lower-bound-noise-ratchet` from a new source. The target it was for (the walled host
exon of a `span` locus under host-only capture, 0.586 against 0.645 through the pipeline) gains ~100
fragments; the cost elsewhere is larger. Re-open with a gate on the emitted level's own width (emit
only where the bracket is tighter than the node's counting) once phase 2's ceiling is in place; until
then the walled overlap exon keeps its lower side from the landscape prior.

### message-layer-open-cases
`priority: next · kind: question · stamped: 2026-09-09 (moved from the sandbox tracker when it was deleted)`
The completion contract's checklist reads ✅ on every row (the ten messages, the level lane, the RNA
level lanes, the ceiling, sj+terminus — `DESIGN.md` §6b.4–§6b.14); what the tracker still carried as
open, each a candidate for the debug loop of `ROADMAP.md` rank 2 and none a hole: **(a) the exon solve
with every face speaking** — an exon's two faces' arrivals are summed at the solve (`_fuse`); whether
the two-witness sum is priced right where both faces carry the SAME intron's claim through two maps is
unmeasured (`ISSUES: transfer-variance-premise` is the neighbouring question); **(b) the factory on a
region carrying BOTH exon and intron bits** — the intron factory runs only where no exon bit is set
(a mixed region is an exon by ruling); whether the density-against-background measurement is valid
there, and under capture, is open "when first needed"; **(c) the chain of termini** — an empty outside
piece (median 12 bp on the ladder) whose far face is another terminus, half the ladder's
terminus-boundary error by the 2026-09-02 census: the level lane now crosses the empty piece
lower-sided; the short-range gDNA-level continuity under one probe footprint as an UPPER bound was
refused with every other upper side (`ISSUES: the-edge-upper-side`,
`ISSUES: levels-always-travel-for-the-gdna-lane`) and the prior serves these slots to ~1.8 % of mass;
**(d) the substrate** — `nest` (a region walled by two termini; measured 2026-09-02: the refit prior
already serves it, 0.666 vs 0.630), `div` and the antisense's nascent variant, one YAML block each
(`docs/TESTING.md` §0a). ⚠ The rung-4 prototype's refusals are `DESIGN.md` §6b.9's record; the two
session census instruments it used (an exon|exon boundary census by flag class, a directed reachability
census, a terminus pair-gap measurement on certified truth) are gone and re-derivable from
`slot_truth.npz` plus the chain in an afternoon — promote one only when a case above needs it.

### two-sided-exon-row
`priority: now · kind: problem · stamped: 2026-09-04`
**What an exon needs before the scan can forward anything on unstranded data, and where it can come
from.** Rung 2's transported row is the intron's composition pulled through the s = 1 face map: its
upper side is the intron factory's own upper tail (soft at mid-gDNA, a cliff at a zero control), so an
exon's held row is a LOWER bound on gDNA wherever the intron cannot distinguish 0.9 from 1.0, and a pass
that forwards it compounds the bias (`policy_prototype.py --by-class` under a forward-backward
prototype: the in-scope unstranded row +10 % when today's rows are forwarded, the stranded half within
4 %). The flux LEVEL as the two-sided pin is REFUTED (`ISSUES: the-certified-flux-row-as-a-level`).
What remains: (a) the crossing composition itself sharpened on the gDNA side — the intron's own solve
(`ROADMAP.md` rank 3, the vertex atom) is the source, because the composition paradigm is invariant to
a common enrichment step and its only residual is the TAPER — the ratio of the spliced fragments'
enrichment to the crossing fragments' at a probed exon's edge (measured +0.25…+0.3 nats under exon
probes, ≈0 under junction probes, 0 off capture; `flux_stage0.py`, session scratchpad; the tool never
sees the probe panel, so the taper is learned or bounded, never read); (b) a bounded marginal over the
taper in the face map in place of the flat top. Neither is built. ⭐ READ AT PASS ZERO FIRST (the
owner's method, `DESIGN.md` §6b.12): in the first pass today's messages leave unstranded exons at
silence (licensed exons 131k silent / 145k transfer on the in-scope test-chromosome row), so every
unstranded exon number in the full pipeline is the prior's work; the refuted flux profile is the only
candidate that solves them first-pass (→ 6k) where its transport is 1, so the actual problem is a
PER-LIBRARY transport the exons can learn (one probe design per library is the reason they share it).
Judge at `R exon (licensed)` and the walled classes, halves apart, pass zero beside the full pipeline.
⭐ STEP F's MEASUREMENT (2026-09-04): the one-sidedness is the face map's PLATEAU above its ceiling (an
unbounded enrichment tolerance), not the crossing count's noise (correlation with the exon's error
−0.09); at pass zero an unstranded licensed exon reads ~9× its true gDNA (+2,000 % at `g25 ss.50 OFF`).
A cap priced by the pair's discrepancies fixes the first pass (−20 to −24 %) and is REFUSED on the
junction-probed panel's weakly stranded rows (`g25 ss.70 ON` 1.122×): under junction probes the flux is
enriched more than the crossing beside it, the ceiling sits below the truth, and no total or weak strand
mode can witness it. The plateau is honest; the first-pass remedy is the solve's (honest ignorance at a
node with a one-sided profile and no own evidence) and the landscape's training population — not a
message. The remaining message-side candidate is a sharper INTRON profile (its own solve).
⛔ FORMS ALREADY REFUSED on the test chromosome (2026-09-03, whole-library, licensed exons in
brackets; the scan's prototype): a TWO-SIDED POISSON row (the crossing count's likelihood under the
exon's hypothesised share, marginalised over the intron row) fixes capture-OFF (`g50 ss.50 OFF`
11,242 [4,071] → 8,945 [1,874]) and is catastrophic capture-ON (`g50 ss.99 ON` 6,367 → 15,558):
under capture the exon interior's gDNA exceeds what its tapered edge crossing implies, and the plateau
was tolerating exactly that; an ABUNDANCE-BOUNDED row (the crossing fixes the level up to an
enrichment step bounded by the exon's total against what the crossing and flux predict) keeps the OFF
wins and still fails at `g50`/`g05` ON (12,333 / 3,090 against 6,367 / 2,293) — a bound on the TOTAL
is blind to the gDNA's enrichment where RNA dominates the total, which is most probed exons; a FLUX
CAP (the route rate × the exon's RNA opportunity as the exon's RNA) claims 19 % gDNA at the zero
control's exons (16,819 → 121,263) — the route rate under-states the contained RNA systematically.
⭐ THE LEVEL LANE'S MEASUREMENT (2026-09-05, `DESIGN.md` §6b.12): with every node reached, a TWO-SIDED
level lane reads −26 % at pass zero on `g50 ss.50 OFF` (the minimum total density of an exon complex
bounds its gDNA from above under capture-OFF) and +33 % on `g50 ss.99 ON` (the same bound is false under
capture); the landed lane is LOWER-ONLY and reads nothing at pass zero on unstranded rows. The gate on the
prize is ONE library-level fact — is this library's gDNA enriched — and its gDNA-specific witness on
unstranded data is the exons of SILENT genes (no spliced fragment at any junction at this depth): their
density against the intergenic density is the enrichment spectrum, at pass zero, with no strand channel.
Where it is flat, the lane may be two-sided. ⭐ OWNER RULING 2026-09-06: the enrichment witness IS the gDNA
landscape prior, and nodes whose only evidence is a bound do not train it; the message rules finish first.
⭐ THE CLIFF INSIDE AN EXON (2026-09-05, the terminus-cluster block under sparse probes): ten isoform ends
126–147 bp into one exon put ten overlapping probes on the cluster and none on the exon's rest; a
lower bound from the ten-times-probed piece over-claims at the boundaries beside it (`g50 ss.99 ON`
10,100 → 12,834, all at the three probed cluster genes). The lane's lower-only law assumes the recipient
is at least as enriched as the source; a probe edge inside an exon complex breaks it, as junction probes
broke step F's cap. A wall-reading strand witness recovers 57 % at +0.8 % elsewhere (not landed; the
plan §5h). The enrichment witness this issue waits for is the cure for this too.
⭐ THE WALL, RE-DERIVED AND RE-REFUSED AT THE SHIP PROTOCOL (2026-09-09): with `transfer` the shipped
default, the toy harness's gate `test_the_harness_REPRODUCES_the_intron_composition_dependence` (an
unstranded two-exon transcript at 60 % gDNA beside a pure-gDNA intron, capture OFF) reads the exon at
|Δf_g| 0.848 dry against 0.098 wet — the relay passed it because its reframe carried a two-sided
Gaussian. The mechanism is step F's: `transport_row` extends the intron's row FLAT above the face map's
ceiling, so a channel-free exon holds two floors (the edge's level and rung 2's plateau) and no ceiling,
and sits at ψ's measured intron reference, which on a dry chromosome is pure gDNA. A WALL above the
ceiling priced by the junction–exon pair's disagreement (`count_price`, the flux level's own price)
closes the gate (0.123 / 0.126) and wins the LADDER's target rows at pass zero (`g05 ss.50 OFF` 0.921×,
`g50 ss.50 OFF` 0.847×) and is REFUSED where step F's cap was: through the pipeline the stranded half
0/6 (`g50 ss.99 ON` 1.142×, `g98 ss.99 OFF` 1.079×, `g98 ss.99 ON` 1.057×), the deferred rows 1.15–1.43×,
the junction and sparse panels' capture-ON stranded rows 1.4–2.5×. The junction–exon pair does not see
the cliff the level crossed (intron → junction), and at `g98 OFF` the flux's own scatter puts the wall
below the truth on exons whose gDNA is the row. The gate is a strict xfail citing this entry; the plateau
stays honest; the remedy is the enrichment witness above, not a wall.

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
`priority: later (dissolves into `ISSUES: arcsine-magnitude-coordinate` — a bounded coordinate has no bracket) · kind: decision · stamped: 2026-08-2x`
ψ's λ bracket was too narrow to express its own prior; `DensityLandscape.required_logodds_window`
derives the correct bracket with no chosen constant (predicted matches measured on every stratum).
Built, gated, priced: nearly every in-scope condition improves, one dense capture-ON rung regresses
marginally, `g05` improves on both strand settings. **Ships OFF** pending two unpriced costs: memory at
genome scale (a small multiple on `sweep_n_grid`) and the end-to-end thermometer. ⚠ The arm harness
that priced it (`ladder_arm_ab.py`) retired with the relay (2026-09-09); re-derive as a
`policy_prototype.py --module` arm on the window, judged by `calibration_vs_oracle.py`.

### alt-splice-rung-unverified
`priority: later · kind: question · stamped: 2026-08-2x`
Do we solve ALTERNATIVE SPLICING correctly? The `alt_splice` toy rung exists and is unverified
(`toy_harness.py --list`). Cheap, and the only structure where several splice junctions share a
BOUNDARY — in scope on all three shipping strata.

### transfer-variance-premise
`priority: later · kind: question · stamped: 2026-08-2x`
Does the message transfer variance correctly price a ratio built on a handful of counts? — PARTLY
answered under the retired relay (its audit instrument retired with it, 2026-09-09): its transfer
variance was a counting term plus a composition term, so every term shrank as either slot deepened and
a deeply-counted transport arrived essentially undamped. The transfer policy prices every hop by both
witnesses' counting plus the pair's own disagreement beyond it (`transfer_rows.hop_price`), which is
the per-hop premise this entry asked for; whether it is right where a pair agrees by coincidence is the
open half. The obvious substitute is refuted in both directions (the landscape's per-slot
posterior is TIGHTER than counting variance; its population spread over-states a fitted premise ~10×).
What survives is a per-hop-type premise; the honest seed is that under capture the posterior LOOSENS
(mode-membership ambiguity). `EQUATIONS.md` §3.5d.

### nascent-stress-sensitivity
`priority: later · kind: question · stamped: 2026-08-22`
Does any in-scope verdict depend on the nascent STRESS level? The ladder runs `on_fraction 0.50`
(development stress); realistic is ~0.10 (`DESIGN.md` §0b). No new panel needed: re-simulate the single
worst in-scope scenario at the realistic level and check whether any RANK moves. A verdict that only
holds at stress is a robustness finding and must be labelled as one.

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
reference as replacement load (the relay-era arm harness's `stage1_pair{,_onesided}` arms) — **STILL REFUSED**:
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

### the-empty-flux-source-at-the-junctions-counting-alone — the sharper price of the empty-piece flux source, PROTOTYPED, A/B'd on the ladder, REFUSED (2026-09-09); the source itself LANDED at the counting price. Do not rebuild the sharper price without the witness-units repair.

`flux-source-skipped-at-an-empty-exon-piece` (CLOSED by landing, 2026-09-09): the transfer policy now
builds a junction's flux level at an EMPTY exon piece too (`transfer._rna_lanes`), priced by `hop_price`
on the piece's zero count — both counts' counting, the rule every hop pays — and the piece emits it with
the flux's own witness (the pooled spliced count on the pooled route opportunity), so the next full node
prices the hop as a full exon prices its flux. ⭐ THE SUBSTRATE: on the ladder over half of all
junction-adjacent exon pieces are empty (2,100–2,200 per strand; a piece shorter than a fragment has no
contained opportunity) but only ~70 per strand have any lane face — a face carrying none of the strand's
bits, the OTHER strand's boundary — and those reach the AMBIG exon|exon boundaries of the overlapping
loci; the test chromosome has none (its junction-adjacent pieces are full), so only the ladder can judge
it. LANDED FORM, the ladder through the pipeline: every in-scope row within 0.5 % (worst `g98 ss.99 ON`
1.0048×, 175,321 → 176,162), the stranded zero controls 0.977× / 0.958× (14,658 → 14,324; 16,043 →
15,365), the unstranded 0.997× / 0.994×; pass zero 6/8 on both halves, worst 1.0051×. ⛔ REFUSED: the
SHARPER price — the junction's counting ALONE at a piece with no RNA opportunity ("nothing to witness
with") — wins the zero controls more (0.944× / 0.970×) and costs the stranded capture-ON rows 1.0179×
(`g98 ss.99 ON`, +3,144 at the AMBIG exon|exon boundaries) and 1.0091× (`g50 ss.99 ON`): a sharp RNA
floor from a probed junction over-reading the exon body reaches a gDNA-rich node whose recipient price
reads the column count (`flux-price-witness-units`), the landed flux level's own open defect on a new
source. The zero count's counting is what keeps it inside the bar.

### relay-od-r-discontinuity — a defect of the RETIRED relay's anchor path, CLOSED with the relay (2026-09-09). Do not rebuild the anchor to look for it.

The relay was discontinuous in `od_r` at ~1e−5 (`g98 ss0.50 capture-OFF`: error 217,531 at `od_r ≤ 1e−7`,
212,581 at 1e−5 — a threshold in the relay/anchor path, not a response; `TRAPS: a-constant-parked-a-value-off-a-knife-edge`).
Bounded at the time: a 1e−5 nudge across all 30 test conditions moved 1/30 rows more than 0.5 % per policy
(worst 1.65 %). The relay, its anchor and the path that carried the threshold were deleted on 2026-09-09;
the transfer policy has no such constant (its hop prices are counting terms and measured disagreements).
The lesson survives as the trap. ⚠ The bound's other half stands for any policy: do not believe a
single-row policy difference below ~2 % without the noise floor re-recorded in the same session.

### levels-always-travel-for-the-gdna-lane — DERIVED (phase 0's finding 1), PROTOTYPED, gated, A/B'd on three panels and the ladder, REFUSED BY THE BAR pending the upper side (2026-09-08). Do not rebuild the gDNA lane on every face before single-strand recipients read the RNA ceiling.

The gDNA lane emitting on EVERY directed face (a composition alongside where a map exists; the solve
reading the composition from a side that sent both) reaches the AMBIG stretches the landed lane starves
(9 of 9 nodes on `span_ab` against 1 of 9) and wins every capture-ON ladder row (unstranded × ON −15 to
−17 %, stranded × ON −1 to −3 %; `capspan` on `g50 ss.99 ON` 957 → 641). It loses the ladder's
`g05 ss.99 OFF` +7.1 % (106,559 → 113,250 at pass zero; 45,076 → 48,295 full) and `g05 ss.50 OFF` +7.8 %
(50,948 → 54,926), the test chromosome's `g00 ss.70 OFF` +7.6 % (9,340 → 10,052) and `g50 ss.50 OFF`
+3.6 %, all one mechanism: a one-sided floor at a node with no channel of its own for the gDNA share is
a TILT, not a floor — the hop price's discrepancy term (an intron against an exon, log r ≈ 4–8) blurs a
step into a slope across the grid, and a monotone likelihood on a flat local posterior moves the median
up the line. On the ladder the floors land on AMBIG walled exons and exon|exon boundaries at near-zero
true gDNA (slot 37345: 10,188 fragments, truth 0, 3,673 → 3,934); on the test chromosome's weak-κ zero
row on one exon's strand profile reading noise as gDNA (f_g ≈ 0.08 at 2,631 fragments), made a floor for
its whole gene through the FORWARD faces an empty boundary could not cross. Stacked with the RNA lanes
delivered at AMBIG nodes it still loses those rows (`g05 ss.50 OFF` 1.124×, `g05 ss.99 OFF` 1.031×).
The RNA lanes alone, on the landed gDNA lane, win every non-zero ladder row of both halves — landed.
⚠ THE FIRST REFUSAL (2026-09-05, the terminus-cluster block): completing the reach by re-reading what a
node holds across KINDS — a held level read as a composition for a composition rule, a held composition
as a level for a lane face, rules into empty recipients dropped — reads 0 % unreached in the census
and worse rows: test `g50 ss.50 OFF` 13,009 → 13,471 (+3.6 %); ladder `g50 ss.50 OFF` 144,288 → 148,516
(+2.9 %), `g05 ss.99 OFF` 45,076 → 46,532 (+3.2 %), `g98 ss.50 OFF` +1.0 %, `g50 ss.99 ON` −1.5 %, g00
identical. Dropping the rules into empties ALONE reads identical on the test chromosome and
+0.2…+5.4 % on the ladder: a composition already crosses a dark exon through the composition rules on
both of its faces (the maps read the boundaries' numbers, not the empty's), worth 5 % on
`g98 ss.50 OFF`. So the law holds for the lane too — what a node holds as a composition is never
re-issued as a level; the information that would reach a cluster's inner boundaries is a lower bound at
an RNA-rich node with weak own evidence, and that is what harms.
Re-judged with phase 2's ceiling in place (2026-09-08, `lat2`): WORSE — `g25 ss.50 OFF` 1.59×,
`g05 ss.70 ON` 1.18×, the junction panel's `g25 ss.70 ON` 1.69× through the pipeline, the weak-κ zero
control 1.8–3.4×. The ceiling reaches only faces without a composition, so it opposes none of the floors
the wider gDNA lane adds at licensed exons, and the two together train the prior on sharper false
brackets. The gDNA lane keeps its landed faces. What would re-open this is a ceiling that reaches the
same exons the floors reach without counting the flux twice — the two-sided exon row
(`two-sided-exon-row`), not this step.

### the-edge-upper-side — DERIVED, PROTOTYPED, A/B'd on three panels and the ladder, REFUSED (2026-09-04). Do not rebuild an upper side on the edge's level, nor a zero-count claim.

Rule 5 (the intergenic|exon edge → exon) is a LEVEL: the exon has at least the edge's gDNA density, at
the count's Poisson width. Three upper sides were tried because a two-sided level wins the zero controls
(the relay's lead there is exactly this: with its anchor off it reads 8,623 at the in-scope zero control
and 39 at the capture-ON one, against the transfer policy's 14,501 and 22,299): (i) undampened — a
first-pass disaster under capture on stranded rows (`g98 ss.99 ON` +75 %, `g98 ss.70 ON` +360 %);
(ii) dampened both sides in log level by the pair's discrepancies — every main-panel row equal or better,
REFUTED on the sparse-probe panel (`g98 ss.99 ON` 36,645 vs 6,981: a dark edge beside a probed exon is
darkness, not absence, and the centre sits far below the truth); (iii) the split form, dampened above
only, a zero count vacuous — safe on all three probe panels (worst 1.03×) and REFUSED BY THE LADDER: no
in-scope win, stranded 0/8, the deferred unstranded capture-ON rows 1.065–1.146×. The physics is
one-sided (capture only enriches the interior over its edge) and no local witness prices the size of
that enrichment on unstranded data, so the upper side cannot be honest. The one-sided level stands; the
zero controls are the landscape's (`gdna-landscape-trains-on-false-positives`).

### the-abundance-discrepancy-map — LANDED as item 6 (2026-09-02), found MISSPECIFIED on the owner's review and REPLACED by the level rule (2026-09-04). Do not rebuild a fitted step or a hypothesis mix.

Item 6 carried the boundary's own strand row into the inside exon through a map with a step `s` between
the "enrichment" and "new RNA" hypotheses, the step's spread FITTED across the served pairs' two
witnesses and a cap `s ≤ r`. Refused on four grounds: (1) a pooled premise — the per-pair discrepancy
rule (`the-pooled-hop-step`, 2026-09-03) refused exactly this for item 7 and was never swept back;
(2) it chose between hypotheses by a fitted prior where the owner ruled that nothing can be assumed —
value kept, precision dampened; (3) the composition currency where the ruling is a level; (4) only
exon|exon termini served. Measured on the pass-form policy: removing it improved the in-scope
unstranded row (`g50 ss.50 OFF` 12,328 → 11,044) and `g25 ss.50 ON` (28,414 → 19,854) and was within
0–4.8 % on stranded rows. The replacement (the owner's design, 2026-09-04, `DESIGN.md` §6b.12): the
level rule from the boundary's MEASUREMENT only, shape-preserving through the level-kept map, per-pair
widths, the crossing total's upper bound without a claim. ⚠ Two forms were tried and refuted on the
way: (i) a Gaussian summary of the boundary's outgoing profile as the level — on unstranded data that
profile is the one-sided curve it holds from the outside exon, and summarising a plateau into (mean,
variance) invents a value the sender never claimed (`g50 ss.50 OFF` 12,328 → 13,911); (ii) the
shape-preserving form made from what the boundary HOLDS — helps unstranded rows but forwards the
imputation (11,503, still behind no-rule at 11,044) and costs 3–5 % on stranded capture-ON rows. The
law that settled it: a level is made from the sender's measurement (its own claim and its total),
never from what it holds. The dampening earns its place: removing it costs 47 % on `g50 ss.99 ON`.
Stage 0 on certified truth (test chromosome and four ladder rows): across terminus faces the inside
region's gDNA density matches the boundary's off capture (median log ratio +0.06 to +0.20, the spread
counting on 13–60 crossings) and shows the taper under capture where the inside exon is probed (+0.44
to +0.89 at exon|intron termini); the totals' discrepancy bounds the gDNA error in 75–100 % of pairs.

### the-certified-flux-row-as-a-level — DERIVED, STAGE 0 ON TRUTH, PROTOTYPED, A/B'd on three probe panels and the ladder, REFUTED by probe placement (2026-09-04). Do not rebuild a flux LEVEL into an exon without a per-face transport.

The candidate: at a CERTIFIED interface the exon's spliced-in RNA is predicted by the face's
route-summed flux rate × the exon's RNA opportunity, the contiguous RNA by the crossing's composition
row × the opportunity ratio, and the exon's gDNA is the remainder of its OWN count — a two-sided
NegBinomial row over the exon's composition, replacing rung 2's one-sided transported row at certified
faces (the relay's anchor re-derived as a composition-currency message, `DESIGN.md` §6b.3). Stage 0 on
the exon-probed test chromosome was clean: the flux-implied RNA is unbiased for the exon's spliced-in
RNA (median log ratio −0.03, mad 0.04–0.09 at depth ≥ 50, on and off capture, probed or not; the
overhang hypothesis refuted; the two faces of one exon agree within counting, median excess variance
0.000, so the local discrepancy rule is a no-op). On that panel it WON: unstranded 9/10 vs the shipped
transfer, ≤ silent 10/10 (`g50 ss.50 OFF` 11,242 → 8,876 vs silent 9,349, licensed exons 4,071 →
1,803); stranded 16/20, worst +1.2 %. **Refuted where the probes move**: on the junction-probed panel
the spliced fragments are enriched e^3.0…e^3.6 (20–35×) MORE than the exon's contained fragments (stage
0, probed exons, mad 0.04), so the row over-claims RNA — stranded capture-ON rows 1.5–4.1× the shipped
error (`g50 ss.99 ON` 6,502 → 26,497), unstranded `g98 ss.50 ON` 1.48×; on the sparse-probed panel
worst 17.4× (`g25 ss.70 ON`) and 5.0× (`g05 ss.50 ON`); on the full ladder unstranded 4/8 with the zero controls at 3.3× (`g00 ss.50 OFF`, 299,380 →
995,880) and 6.9× (`g00 ss.50 ON`) and `g05 ss.50 ON` at 4.1×, stranded within 2 % except the two
`g00 ss.99` rows (1.2–1.3×); at `g00 ss.50 OFF` the mechanism being tiny exons (3 contained fragments, median) beside
a large flux where a count-based row has no resolution, whose false modes the refit prior then trains
on (`ISSUES: gdna-landscape-trains-on-false-positives`; walled exons 106k → 371k, terminus boundaries
93k → 296k) — ⚠ a PRIOR-MEDIATED number: on the test chromosome's zero control the same row is slightly
BETTER than transfer at pass zero (266,954 vs 278,108) and worse only through the pipeline. The relay's own transport-centre estimator (`rna_anchor.left_fit_center_spread`) REFUSES
on 9 of 10 conditions tried and misfires where it accepts (`g50 ss.99 ON` 6,006 → 32,235). Joining the
row beside rung 2 (the relay's placement) double-counts the flux and the intron row and fails the same
way. A hard s ≥ 1 truncation of the crossing marginal is refuted on its own (5,712 vs 1,803 at the
in-scope unstranded row). This is `DESIGN.md` §6b.9's founding refusal re-measured: a level carried
between locales under capture is refuted by probe placement alone, and the tool never sees the probe
panel. What survives is recorded in `ISSUES: two-sided-exon-row`.

### the-pooled-hop-step — DERIVED, PROTOTYPED, A/B'd, landed for a day, REFUSED by the owner (2026-09-03). Do not rebuild a per-library premise.

Item 7's messages (`DESIGN.md` §6b.8) first carried a POOLED premise: per hop kind, the precision-
weighted mean of the sighted pairs' disagreement (the boundary's own strand mode against the flank's
mapped to it), applied to every pair as a shift with its standard error as width — justified by a
named mechanism (capture's taper of gDNA at a probed exon's edge) and by measurement (the un-premised
form's +1.9 % at ladder `g50 ss.99 ON` turned into −0.7 %). Refused on three grounds: (1) the other
pairs' behaviour is not known to predict a given pair's — truth showed the offset locus-dependent
(−0.07…+0.49 nats across the test chromosome, none on the ladder's C hop, length-dependent on its E
hop); (2) the pooled shift carries ~1/n of the recipient's own mode back into the message it receives;
(3) on the ladder the owner's per-pair discrepancy rule alone stands within 0.25 % of it on every row and
ahead on four of six stranded rows (local vs pooled against the pre-item-7 policy: `g05 ON` −1.34 vs
−1.38 %, `g50 ON` −0.44 vs −0.67 %, `g98 ON` −2.70 vs −2.49 %, `g98 OFF` −0.76 vs −0.75 %, `g05 OFF`
+0.12 vs +0.17 %, `g50 OFF` +0.06 vs +0.08 %; unstranded identical). Where the pooled shift wins is
the test chromosome's low-gDNA capture-ON rows and the junction-probed panel (junction `g05` +8…+9 %
local vs +3.5…+5.6 % pooled; benign `g50 ss.99 ON` 3,933 vs 3,873; sparse 4,118 vs 4,011): a
systematic offset of the licence below each pair's counting, which only pooling can see and which the
owner declines to extrapolate. ⚠ A first "local" measurement reported the pooled form's E hop by
mistake (a leftover branch in the prototype) — corrected the same day; log what an arm APPLIES before
reporting it. Owner: start with the simplest elegant form; a global model of disagreement (the
production relay's projection of each node onto the total-abundance landscape) is an option for after
the tool works end to end, not before.

### the-flux-factor-hop-premise — DERIVED, PROTOTYPED, A/B'd, REFUSED (2026-09-02). Do not rebuild it as a fitted parameter.

The E hop of the alternative splice site (`DESIGN.md` §6b.8) reads the boundary through `U/(U+S_b+F)`;
where a panel captures the junction's spliced fragments differently from the unspliced crossing, the
measured flux `F` is off by a factor, and the premise "fit `1/delta` on the flux from the sighted
pairs" was built two ways. ALONE (the flux factor as the E hop's only premise): fitted 2.07 ± 0.26 on
the benign `g50 ss.99 ON` row and HARMED more than the un-premised arm (boundaries 231 → 280 node-
locally; the sparse panel 394 → 562) — this simulator's probes are tiled in transcript space, so there
is no junction depletion there, and what the pairs see is E's own edge taper, an odds step. JOINTLY
with the step (a profile over `k` with the step in closed form): at n ≈ 14 the two are degenerate — the
profile picked the step on the junction panel (k = 1.14 ± 0.30, E flanks 294 → 210) and the flux
factor on the benign row (k = 2.28 ± 0.20) and on the sparse panel (k = 2.59 ± 0.13 with χ² = 47.7 for
14 pairs: a tight, wrong commitment; boundaries 382 → 526). The mechanism is real at junction-probed
panels (the junction panel's E flanks carry +106 fragments under the step-only form, recorded as an
owner decision) but is not identifiable from the sighted pairs beside the step; one mechanism per arm.

### arcsine-magnitude-coordinate — DERIVED, PROTOTYPED, A/B'd, REFUSED (2026-09-01). Do not rebuild it.

⭐ Replacing the simplex solver's magnitude axis `λ = logit(f_g)` with the arcsine coordinate
`f_g = sin²φ`, so that a vertex becomes a finite interior endpoint. Steps ①–③ all ran (dissection →
derivation → prototype outside `src/` → per-stratum A/B); `src/` was never touched. The prototype and
the full derivation are in commit `48be3d0a`, one commit before their deletion.

⛔ **THE HEADLINE PREMISE IS FALSE, and that is the durable finding.** (a) A posterior MEDIAN cannot
reach a vertex in ANY coordinate: for the outer bin holding mass `m` the CDF crossing is
`t = 1 − 0.5/m ≤ ½`, so the read-out never passes the outer bin's midpoint. (b) In `f`-space the λ
grid is **FINER at the vertices than uniform-φ** — `σ(±L)` already sits `4.54e-5` from the vertex; at
K=60 the top `f`-gap is `1.83e-5` (λ) vs `1.37e-3` (φ), and the one-hot read-out is `0.999954602` (λ)
vs `0.999828662` (φ). `logit` CONCENTRATES resolution at the vertices; uniform-φ spends it evenly, so
the change is a resolution REALLOCATION (~3× finer mid-simplex, ~75× coarser at the extremes).

⛔ **THE LEVERAGE IS BOUNDED AT ~2.5 % OF THE DEFECT** (K-invariance, `g98 ss.99 ON`, rung C). Each
coordinate converges to its OWN limit (mass-wtd `|K60 − K240|`: base 1.47e-04, arcsine 1.31e-04) while
the gap BETWEEN them holds at ~6.5e-04 and grows slightly with K — the λ bracket's truncated DOMAIN,
which refining a grid cannot recover — against a 0.026 shortfall in the vertex band.

⛔ **THE A/B SPLITS BY STRATUM AND LOSES WHERE IT MATTERS** (whole ladder, release metric Σ|Δ|,
`arcsine/base`, noop byte-identical on all 16): stranded × OFF **1.008 / 0.939 / 0.898** (wins, and
more as gDNA rises); stranded × ON 0.998 / 1.010 / 1.037; **unstranded × OFF 1.107 / 1.071 / 1.014 —
loses EVERY row**; deferred 1.024 / 1.007 / 1.002; `g00` controls 0.824 / 0.925 / 0.978 / 1.037.
⭐ The mechanism the split names: the coordinate helps where the strand LIKELIHOOD is strong and hurts
where it is flat — there the posterior essentially IS the reference, so the grid's tail resolution
moves the answer directly. Losing all three IN-SCOPE unstranded rows is what refused it; a pooled
total would have hidden it (`TRAPS: never-pool-the-strata`).

⭐⭐ **WHAT SURVIVES AND IS WORTH KNOWING.** The derivation is CORRECT and was confirmed numerically:
under `f = sin²φ` the measure conversion `log(dλ/dφ) = −½·log f − ½·log(1−f) + log 4` is EXACTLY minus
the two written Jeffreys reference halves, so under φ **neither a Jacobian NOR a reference is
written**; uniform-φ weights reproduce `Beta(½,½)` bin mass to `3.3e-15`; the two solvers agree to
`2e-8`; the tilt's fact-3 cancellation survives (the two Jacobians are independent). ⚠ **Rung B —
φ-native message DELIVERY, mapping a claim's precision by `(dx/dφ)²` — was derived but NEVER
MEASURED**, and is the honest candidate for the shipped-stage spread (the message-free local solve
behaves: C median ≈ 0.997). It cannot exceed the 2.5 % bound. ⚠ Also never run on the arcsine arm:
`zero_controls.py`, the 30-condition test chromosome, and a same-session reseed floor.

⭐ **AND THE OBVIOUS RE-AIM WAS PRICED IN THE SAME SESSION AND IS ALSO NOT A WIN.** The near-vertex
under-call is an HONEST posterior width — on the STRANDED condition measured it tracks
`w/4 = 1/(2·√n)` across four decades of depth (obs/pred 0.63…1.21), reproducing `EQUATIONS.md` §9a's
recorded log-log slope of −0.5221. ⛔ §9a also records that on the real UNSTRANDED fit the shortfall
is depth-INDEPENDENT, so **never quote an `n^(−1/2)` shrinkage on an unstranded stratum.** §9a/§9a.1
already carry the theorem and its spike-and-slab exception, and §9d.4 already carries a fully derived
atom with no new constant. **What this session ADDS is the pin-and-re-solve**: `vertex_ceiling.py`
(3 conditions; `noop` pins 0 and is byte-identical, `vertex_free` pins 69,850–78,820) gives
`vertex_free` region Σ|err| **+18,887** and boundary **−20,208** — **a WASH** — and `vertex_all`
+94,490 / −18,237, clearly worse. ⭐⭐ **Supplying the vertex information PERFECTLY nets ≈0, because
handing an object the exact answer changes what it BROADCASTS and the relay over-propagates it.**
⛔ That is a stronger statement than §9a's "value of missing information", and it is why the next
thread is the MESSAGE LAYER and not the prior: if a perfect local answer cannot survive propagation,
no local improvement can pay off. ⚠ Three conditions, on `vertex_ceiling`'s pass-0-flavoured metric —
re-price the atom on the whole ladder AFTER the message layer is repaired, not before.

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

### the-message-policy-campaign — SIX MECHANISMS BUILT, MEASURED AND REFUTED; the campaign closed 2026-08-27. Do not re-run these without new evidence.

⭐ Moved from the sandbox record when it was deleted (2026-09-09). A campaign to replace the relay
with a derived message policy ran for several sessions and did not reach the bar; the code every row
describes was DELETED on 2026-08-27 and git carries it. ⛔ **The bar it missed is the one the transfer
policy was later judged by: NOT to beat `SilentPolicy`** — on strand-specific data a sighted exon's own
solve is excellent and a message can mostly only disturb it; the goal is to perform on UNSTRANDED data
while doing minimal harm on stranded data. The campaign repeatedly optimised a pooled total, which
hides a sign flip between the two halves. ⭐ The one measurement that survived everything: propagation
is net-harmful wherever the local solve HAS its own evidence, and its value concentrates where that
solve is BLIND — the reason the two halves are judged against different bars.

| mechanism | what happened |
|---|---|
| **A composition-transporting policy** (`CurrencyPolicy`) | Best zero controls ever measured, but lost every in-scope contaminated stratum to silence. Deleted. |
| **The gDNA-continuity rule** (an unsupplied source's gDNA level crosses unscaled) | Built THREE ways — a static per-slot licence (a value RATCHET, gDNA densities to 3.9e+32, from breaking the knob's telescoping cancellation), a running-state licence (killed the ratchet, still lost capture-ON), and a fuse-based pure-gDNA re-anchor lattice (too weak — a fuse negotiates where the relay's mass rescale overwrites). Halved one row, lost others. Only safe beside a scan-time mass rescale. ⭐ The transfer policy's LEVEL LANE is this rule rebuilt as a one-sided profile with a priced hop (`DESIGN.md` §6b.12) — a different mechanism, judged separately. |
| **The premise's exon-end scoping** | The dispersion decomposition proves intron-end hops carry no COMPOSITION cost, yet scoping the charge to exon-end hops REGRESSED the panel: freeing intron chains before a measured LEVEL charge exists releases un-priced level drift. Restoring the pooled charge recovered `g98 ss.50 ON` 4.37 M → 2.44 M. |
| **A class-keyed method-of-moments fit on the observed log-ratio** (as the runtime law for the transport variance) | Tracks truth at intron and plain classes, REFUTED at sj classes: the route-summed flux cancels the visible step exactly where the true error is largest (0.08 observed vs 4.4 true). |
| **A totals-form pair fit** (for the same variance) | Refuted by construction — the knob consumes the totals, so transported totals agree by the (1−w) algebra and carry almost no information. |
| **`FanOutPolicy`** | Measured dominated once the certified-flux anchor gave its destinations own evidence. Deleted 2026-08-24. |

**Derivations the campaign left, recorded so they are not re-derived.** (i) The conservation identity
is a COUNT identity, ``Σ_c ρ_c·E_c = M`` — each component's density weighted by its OWN opportunity,
summing to the slot's observed unspliced count; summing raw densities against the reciprocal-opportunity
total is a DIFFERENT identity, exact at boundaries and WRONG at regions, where that total reads
``ρ·P(w ≤ ℓ)`` — measured on the ladder, the median exon's P is 0.452 and the 5th percentile 0.044, so it
was up to 23× too low at half of all exon slots. (ii) Two variance kinds: a reframe's cost multiplies
every lane of a message identically (a LEVEL statement) while each component also carries its own (a
COMPOSITION statement); spending the shared part as per-component variance converts a common-mode
level error into a composition error — the pathology where a near-zero gDNA claim eats an unstranded
slot's unexplained RNA mass. (iii) The transport dispersion decomposed against certified truth (four
`g50` corners, noise subtracted): intron↔boundary hops are FREE; all structural error lives on
exon↔boundary hops and is predominantly COMMON-MODE, exploding under capture, from two derivable sources
— the truncation frame term ``log(1/P(w ≤ ℓ))`` (pure geometry) and the capture step (common-mode
because probes bind gDNA and RNA alike). `transport_dispersion.py` is the instrument;
`ISSUES: flux-floor-dispersion` is the open thread.

**Method lessons** (the durable ones have named homes in `TRAPS.md`): attribute before iterating — three
laws landed together, the panel regressed, and only an attribution factorial (each law removed alone, all
16 conditions) named the carrier; a cancelling defect pair reads as success
(`TRAPS: a-cancelling-defect-pair`) — fixing the conservation identity made low-gDNA WORSE because the
truncated total had been suppressing claims at short exons and hiding phantom gDNA; and layering fixes on
fixes is how a policy becomes unmaintainable — the campaign ended with three conservation operators behind
flags, two of which were provably limits of the third.

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
