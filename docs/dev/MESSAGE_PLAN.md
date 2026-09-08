# THE MESSAGE PLAN — the ten rules, the missing rules, and the repairs (owner review, 2026-09-04)

    ⚠ A DEV DOC and a PLAN. It answers the owner's review of the landed policy and lays out the work
    to an intact message design for every face — the bar for the flip. Rulings live in `DESIGN.md`
    §6b.12; the case-by-case state in `MESSAGE_RUNGS.md`; this file says what is wrong, what is
    missing, and in what order it gets fixed. MOVE anything that settles.

## 0. The owner's rulings this plan is built on (2026-09-04)

* A message carries up to five optional lanes: the gDNA-vs-RNA composition PROFILE, a TILT profile
  (the RNA+ vs RNA− degree of freedom, off at single-stranded nodes), and a LEVEL claim
  `(log rate, log-variance)` each for gDNA, RNA+ and RNA−.
* Across a face where composition cannot cross, a LEVEL is what travels. Its VALUE is the source's
  level, kept; its PRECISION is dampened by the degree of abundance discrepancy between source and
  destination. ⛔ No hypothesis is chosen for the discrepancy (capture, new transcription, noise — we
  cannot know), so nothing is rescaled and nothing is fitted across pairs.
* Decisions on the recorded residues are deferred until the contract is fulfilled (no missing rules).
* The message architecture is finished end to end, becomes the default, the older policies retire,
  and only then the landscape prior's training population is taken up.

## 1. One word, defined once: "flank"

The code says *flank* for THE REGION ON EITHER SIDE OF A BOUNDARY (`outside_flank`, `junction_flanks`
in `transfer_rows.py`). Example, rule 9: an alternative 5′ splice site inside an exon splits it into
two exon pieces, A on the left and B on the right, with the junction boundary between them; isoform 1
splices out at that boundary, isoform 2 continues into B. "Flank → boundary" is "exon piece A (or B)
→ the junction boundary". This plan writes "the region on each side of the boundary" and keeps the
code's names only in backticks.

## 2. The ten rules today — status, and what is wrong with each

| # | face | currency | status | the defect |
|---|---|---|---|---|
| 1 | intron → intron\|exon boundary | composition, FORWARD | sound | — |
| 2 | intron\|exon boundary → intron (one shared strand) | composition, FORWARD | sound | — |
| 3 | boundary → exon at a licensed face | composition, the splice-in face map | **to refine** | ONE-SIDED: the wall ("at least this much gDNA") is set by the face's crossing count (12–25 fragments, ±30–50 % noise) — a one-way ratchet; the plateau is the intron's own soft high side; the flat top is an implicit enrichment allowance (§5.F) |
| 4 | exon → boundary at that face | composition, the map read backwards | **to refine** | no discrepancy widening: under capture the pair's two witnesses disagree by the taper (measured 1.3× benign, 2.2× junction-probed) and the rule charges nothing for it (§5.G) |
| 5 | intergenic\|exon edge → exon | LEVEL, consumed as a one-sided profile bound | **misspecified** | the profile over an UNBOUNDED enrichment nuisance makes a ZERO-count edge vacuous: at a zero-gDNA library every edge says nothing, where the relay's level transfer says "gDNA density is zero" — the relay's zero-control lead (ladder `g00 ss.50 OFF`: 58,840 vs 233,420) is largely this (§5.E) |
| 6 | region outside a terminus → the terminus boundary | composition, splice-out with the spliced crossing | sound | — |
| 7 | terminus boundary → the region outside | composition, the face map | sound | — |
| 8 | terminus boundary → the region INSIDE | "abundance-discrepancy" map on the composition axis | **misspecified** | (i) a POOLED step spread across served pairs, refused for rules 9–10 the next day and never swept back; (ii) mixes the two hypotheses through that fitted prior instead of keeping the value and dampening; (iii) composition currency where the ruling is LEVEL; (iv) registered only where both sides are exons — the exon\|intron terminus's inside exon (4,130 ladder boundaries) gets nothing (§5.A) |
| 9 | exon piece → alternative splice site boundary | composition, splice-out with the pair's own discrepancy width | sound | the low-gDNA capture-ON residues are deferred decisions |
| 10 | alternative splice site boundary → exon piece | composition, the face map with the pair's discrepancy width | sound | as 9 |

**How rule 5 was handled before.** The relay sent the edge's gDNA density into the exon UNSCALED (its
`_may_share` licence fails for an edge, which has no RNA component, so the reframe ratio was 1 — a
density transfer), two-sided, with its precision damped by the two slots' total-density log-variance.
That is a level transfer with counting-based damping. The transfer policy replaced it with the profile
likelihood `sup_{s ≥ 1} Pois(n_b; c/s)`, which is exactly right about the sign (capture can only enrich
the interior relative to its edge) and exactly wrong at zero: any enrichment explains zero counts, so
the bound is vacuous where the evidence is strongest. Yes, this is a problem.

**How rule 8 slipped.** Item 6 landed on 2026-09-02 with the pooled spread; the per-pair ruling that
refused pooling came with item 7 on 2026-09-03 and was applied forward only. A ruling must trigger a
sweep of every earlier rule carrying the same premise; this plan is that sweep.

## 3. The solve with mixed message kinds — one procedure

At solve time a node may hold two composition messages, a composition and a level, or two levels, plus
its own evidence and the hyperprior. The procedure is one: every LEVEL lane is converted AT THE NODE
into a profile over the node's own composition through the node's own total (an observation the
contract allows) — a two-sided level becomes a Gaussian in log-level read along the node's λ grid, a
lower-bound level a one-sided profile — and then every profile adds: the two held messages, the node's
own evidence, the prior. The tilt lane adds on the θ grid the same way. Nothing else is needed for
(a), (b) or (c), and a lane that is `None` costs nothing.

## 4. THE LEVEL RULE — the owner's design, written as a rule

For a directed face `s → i` where composition cannot cross:

* **What crosses.** The gDNA level always (gDNA is genomically continuous). A strand's RNA level only
  where that strand's population continues across the face — as a LOWER bound where new transcription
  can add on the far side (the inside of a terminus: `RNA_known = the unspliced RNA crossing + the
  spliced-in flux`), two-sided where the population is unchanged. At a strand-change face no RNA
  crosses.
* **The value.** The source's level per population, in counts per base of that population's
  opportunity, kept. At a boundary the gDNA level is its composition profile times its crossing
  density (the profile's median and width, delta method, in log-level); at a region likewise on its
  contained density.
* **The precision.** Counting, plus the dampening the owner ruled: the excess of the totals'
  discrepancy over its counting variance, PER PAIR — `max(0, (log r)² − var_count(log r))`, with
  `r = (destination total / its opportunity) / (source total / its opportunity)` — added to the level's
  log-variance. The same shape as rules 9–10's per-pair rule; nothing pooled, no hypothesis, no shift.
  At `r = 1` the level transfers with counting width alone.
* **The conversion.** At a node that has a total, the level becomes a profile over that node's λ grid
  (§3). At an EMPTY node (no total: an empty chain piece) the rule FORWARDS the level lane unchanged —
  the discrepancy cannot be measured there and is measured at the next node that has a total.
* **What it serves.** The inside region of every terminus (exon|exon and exon|intron), strand-change
  faces, termini both ways, the empty chain pieces between termini, and — with the composition rule —
  the junction-plus-terminus faces and the AMBIG region's two degrees of freedom.

## 5. The designs, case by case

**A. Rule 8 re-specified as the level rule (the inside of every terminus).** Faces: exon|exon termini
(7,973 ladder boundaries) and exon|intron termini (4,130). Lanes: `level_gdna` two-sided, dampened by
the pair's discrepancy; `level_rna_<strand>` as the lower bound `RNA_known`. Stage 0 on certified
truth: the inside region's gDNA density against the boundary's (1.0 off capture; the taper under
capture — that is what the dampening prices) and the RNA lower bound's slack (the new-transcription
share). Falsification: the orientation reversed (the outside region served as if inside) must fail;
the dampening removed must fail on the sparse-probe panel where the discrepancy is capture.
⛔ Replaces `abundance_row`, `abundance_map` and the pooled `v_step`.

**B. Strand-change faces (1,138) and termini both ways (48).** `level_gdna` only, dampened. Stage 0:
gDNA density continuity across those faces on truth. The both-stranded locus has no substrate yet;
`test_chr.yaml` grows one structure for it.

**C. The empty chain pieces (5,588 terminus boundaries; 3,792 reached only through another
terminus).** The forwarding case of §4: the empty piece forwards the level lane; the next boundary
with a total converts and dampens. Substrate: `nest` (a region walled by two termini) is the
structure that falsifies it.

**D. Junction plus terminus on one face (1,968 + 386).** The licence is the derivation owed: the
population that crosses the boundary is the outside region's (unspliced) plus what enters spliced
through the junction; the inside region receives the level rule (A) with `RNA_known` including the
flux; the junction's other side receives rules 9–10's composition where it shares the population.
⚠ Not a simulator artefact: the counts are the ladder's, carved from the real annotation (alternative
last exons whose end is another isoform's donor, and the like). Substrate: `instart`.

**E. Rule 5 as a level.** `level_gdna` = the edge's gDNA density with its counting variance (the Gamma
posterior of the count, which at a zero count says "below about one fragment per opportunity", never
"nothing"). Two forms to A/B, because the owner's dampening and the one-sided allowance pull apart
here: (i) one-sided — the exon's gDNA is at least the edge's level (capture can only enrich the
interior), consumed as a one-sided profile that is NOT vacuous at zero; (ii) two-sided at the edge's
level with the per-pair discrepancy dampening. ⚠ At an expressed gene the totals' discrepancy is RNA
and always large, so form (ii) is weakest exactly where form (i) is informative; the zero controls
and the capture-ON rows decide, halves apart, pass zero beside the full pipeline.

**F. Rule 3 made two-sided at the face.** ⚠ Stage 0 (`ratchet_stage0.py`, 2026-09-04) corrected the
mechanism before anything was built: the crossing count's draw does NOT predict the licensed exon's
signed error (correlation −0.09 to +0.2), so the ratchet below is minor. The mechanism is the PLATEAU:
`transport_row` clamps the preimage above the face map's ceiling, leaving the profile flat to f = 1 as
an unbounded enrichment tolerance, and the posterior median of that plateau in log-odds sits near 1 —
at pass zero an unstranded licensed exon is over-estimated ~9× on the in-scope row (+2,000 % at
`g25 ss.50 OFF`), and the landscape repairs it while training on it. The design that follows steps A
and E: the tolerance above the ceiling is priced by the pair's discrepancy — a log-level Gaussian fall
of variance counting + the totals' disagreement beyond counting — never unbounded (`cap3_proto.py`).
**Step F's verdict (2026-09-04): the cap is REFUSED, by the same law that refused the edge's upper side.**
With the strand witness the junction panel's `g25 ss.99 ON` came back to 1.022×, but its weakly
stranded `g25 ss.70 ON` reads 1.122× — at κ = 0.7 the exon's own strand mode cannot witness the
disagreement beyond counting, its own solve is too weak to resist, and the cap pulls the exon toward a
ceiling that under junction probes sits below the truth (the flux is enriched more than the crossing
beside it, which no total can see). A cap wide enough to be safe there is no cap at all on unstranded
data, which is where its value was (first pass −20 to −24 %). So the plateau stays: it is the honest
statement of an enrichment ambiguity with no local witness. What that leaves: the first-pass estimate of
an unstranded exon under a one-sided profile is ~all gDNA, and the remedy belongs to the SOLVE and the
LANDSCAPE, not the message — a node with a one-sided profile and no own evidence is honest ignorance
(`solvability_audit.py`'s rule), and the landscape's training population must exclude it (the parked
issue, the owner's "which regions and boundaries train the landscape"). Recorded in
`ISSUES: two-sided-exon-row`.

The earlier plan for the wall's position stands as a refinement: the wall's position comes from the INTRON's gDNA level —
its density (hundreds of fragments) times its composition — carried by the opportunity ratio, not
from the 12–25-fragment crossing count; the crossing count becomes the pair's discrepancy WITNESS
(its excess over counting dampens), which is also where the capture taper shows. That removes the
ratchet. The plateau above the wall is the intron's own high side and is the intron's own solve to
sharpen (parked), unless a bounded taper marginal replaces the flat top; both are A/B'd after A–E.

**G. Rule 4's premise — MEASURED NEUTRAL, not landed (2026-09-04).** The per-pair width on rule 4
(`rule4_proto.py`: the boundary's own strand mode against the exon's mapped through the splice-out law,
the excess over counting widening the pair's message) fires at 2–65 faces per condition and moves the
eight stranded rows by −0.2…+1.0 %; applied both ways −2.2…+2.7 %. The premise bias under capture is a
systematic OFFSET (1.3–2.2×) that a per-pair width cannot see — the same conclusion as the refused pooled
shift. Nothing to land; recorded.

**G (as first planned). Rule 4's premise.** The per-pair discrepancy rule of rules 9–10 applied to rule 4 (and 3): where
the exon's mapped composition and the boundary's own strand composition disagree beyond counting,
widen this pair's messages by the excess; on unstranded data there is no second witness and the rule
stays at counting width.

**H. The AMBIG node.** Its two degrees of freedom are imputed one at a time by single-stranded
neighbours' RNA levels (`level_rna_pos`, `level_rna_neg`) and by the tilt profile where a neighbour is
itself both-stranded. A ruling first (rung 5, last by the owner's order), then the substrate.

## 5b. STEP A — what was derived, measured and landed (2026-09-04)

**Stage 0 on certified truth** (`level_stage0.py`, test chromosome and four ladder rows): across terminus
faces the inside region's gDNA density matches the boundary's off capture (median log ratio +0.06 to
+0.20; the spread is counting on 13–60 crossing fragments) and shows the taper under capture where the
inside exon is probed (+0.44 to +0.89 at exon|intron termini); the totals' discrepancy bounds the gDNA
error in 75–100 % of pairs; the RNA lower bound holds (violations 3–12 %) with a median slack of
0.2–1.2 nats (new transcription).

**Three forms tried, one survived.** (i) A Gaussian summary of the boundary's outgoing profile as the
level lane — REFUTED: on unstranded data the boundary's outgoing profile is the one-sided curve it holds
from the outside exon, and summarising a plateau into (mean, variance) invents a value the sender never
claimed (`g50 ss.50 OFF` 12,328 → 13,911). (ii) The shape-preserving form — the profile read through the
level-kept map, blurred by counting and the discrepancy — helps unstranded rows but still forwards the
imputation (11,503, still behind no-rule at 11,044) and costs 3–5 % on stranded capture-ON rows.
(iii) THE LAW THAT SETTLED IT: **a level is made from the sender's MEASUREMENT (its own claim and its
total), never from what it holds** — an imputation is not re-issued as a level. The dampening earns its
place: removing it costs 47 % on `g50 ss.99 ON`.

⚠ **A prototype bug, caught by the per-slot identity of the landed form against it.** The prototype
popped rules by geometry at every terminus boundary and thereby popped rung 3 at every GENE EDGE (an
intergenic|exon face carries a terminus flag), so every `level8_*` arm above was measured WITHOUT the
edge bound, and its "removing item 6 helps 10 %" was partly "removing rung 3 helps". The clean
attribution is source against source (`level_off.py`: the landed policy with only the level rule
removed): the level rule itself is INERT on the in-scope unstranded row (12,030 vs 12,025 — only the
upper bound can fire there), a small harm at `g50 ss.99 ON` (+1.2 %), small wins at `g98 ss.99 ON`
(−0.7 %) and `g25 ss.50 ON` (−1 %), and a large win at the capture-ON zero control (26,594 → 22,299,
the upper bound capping phantom gDNA). Against the committed policy (item 6 in place) the landed policy
wins the main panel's unstranded half 9/10 (worst 1.009×) and the stranded half 11/20 (worst 1.013×) —
most of that gain is item 6's removal. ⛔ Lesson for the method: a rule popped "by face" can be a
neighbour's rule; the per-slot identity gate against the source is not optional.

**Landed** (`messages/transfer.py`, `transfer_rows.level_map_lambda / level_row / level_bound_row`):
a rule now receives the sender's own claim and what it holds AS TWO ARGUMENTS and decides — a
composition rule composes and maps, the level rule reads the measurement only; item 6's `abundance_map`,
`abundance_row` and pooled `v_step` are deleted; both terminus kinds are served. Gates: the level map's
arithmetic (the level kept, monotone, the bound vacuous below the total), the rule at every served pair
against an independent recompute, the held imputation never crossing, the upper bound without a claim,
the rule vanishing with the terminus bits, and NO POOLING (another pair's counts leave this pair's
message unchanged). Standings of the landed policy: §5c below.

## 5c. Standings of the landed policy after step A — against the committed policy (item 6 in place)

`compare_logs.py`, the two halves apart; "old" is the committed pass-form policy, "new" the landed
level rule. Main test panel: unstranded new wins 9/10 (worst 1.009× at `g98 ss.50 ON`; `g00 ss.50 ON`
26,594 → 22,299, `g50 ss.50 OFF` 12,328 → 12,030), stranded 11/20 (worst 1.013×). Junction-probed
panel: unstranded 8/10 (worst 1.003×), stranded 12/20 (worst 1.007×). Sparse-probed panel:
unstranded 9/10 (worst 1.036×), stranded 11/20 (worst 1.022×). Suite: 3,757 passed / 8 xfail / 0
failed. The ladder (`landed_level_ladder.out`): stranded new wins 6/8, worst 1.001× — the minimal-harm
bar met; unstranded 4/8 — the zero controls 233,420 → 195,268 (0.84×) and 305,626 → 230,800 (0.76×),
the OFF rows within 0.2 %, and the three DEFERRED unstranded capture-ON rows 8 % worse (`g05` 1.084×,
`g50` 1.081×, `g98` 1.079×): item 6's forwarded lower bound was helping there, and the measurement law
(a level from the sender's own claim only) removes it; on unstranded data the terminus face now carries
only the crossing total's upper bound. Reported, not a target (`CLAUDE.md`'s scope ruling); the
two-sided exon profile (step 3) is what would give those rows a first-pass claim to forward.

## 5d. STEP E — rule 5 as a level: derived, measured, and the ladder's verdict (2026-09-04)

**Stage 0 on truth** (`edge_stage0.py`): off capture the exon's gDNA density equals its edge's (median log
ratio +0.04 to +0.05, counting spread on ~12 edge fragments); under capture a probed exon's interior is
1.25× its edge on the test chromosome (q10 +0.16, q90 +0.33) and 2.3× on the ladder — the taper,
one-sided as argued; zero-count edges are 30–47 % of gene edges at `g05` and capture-ON (63–91 % at
unprobed exons) and 100 % at `g00`. **The relay's zero-control lead is this level, not its anchor**:
with the anchor off the relay reads 8,623 at the in-scope zero control (unchanged) and 39 at the
capture-ON one.

**Forms measured** (`edge_proto.py`; `pass0_score.py` for the first pass), each against the committed
policy: (i) the two-sided Poisson level, undampened — zero controls 14,501 → 8,640 and 22,299 → 2,223,
the in-scope unstranded row 12,030 → 9,627 (first pass 171,260 → 77,851), and a first-pass disaster on
stranded capture-ON rows (`g98 ss.99 ON` +75 %, `g98 ss.70 ON` +360 %); (ii) both sides dampened in
log level by the pair's discrepancies — equal or better on every main-panel row, and REFUTED on the
sparse-probe panel (`g98 ss.99 ON` 36,645 vs 6,981): under sparse capture a probed interior is captured
many-fold more than its dark edge, the level's centre sits far below the truth, and a node with no
evidence of its own settles at the centre — a variance cannot fix a bias; (iii) the split form, exact
Poisson below and a dampened log-level Gaussian above, a zero count vacuous — safe on all three probe
panels (worst 1.03×), wins on the main panel's unstranded rows (`g05/g25/g50 ss.50 ON` 0.82–0.86×),
and REFUSED BY THE LADDER: in scope within 0.5 % (no win), stranded 0/8 wins (worst 1.004×), the
deferred unstranded capture-ON rows 1.065–1.146× — the upper side, even dampened, pulls an unstranded
probed exon toward a centre 2.3× below the truth.

**What stands.** The one-sided level: the exact Poisson below the edge's level ("the exon has at least
the edge's gDNA density", at counting width), NOTHING above it, a zero count vacuous — rung 3's formula,
now derived as the level rule's lower side rather than a profile over a nuisance, with the reasons the
upper side cannot exist recorded: no local witness prices a probed interior's enrichment over its edge,
and a dark edge under capture is not an empty one. `edge_level_row` carries that derivation; the
policy's rule reads the count alone. Byte-identical to the committed step-A standings (checked on three
conditions). The zero controls' win that "zero means zero" would buy is the relay's capture-blind lever;
those rows are the landscape's to win (its parked issue). Lesson for the method: the exon-probed test
chromosome passed two forms the sparse panel and the ladder refused — all three panels AND the ladder,
every time.

## 5e. Standings after step E

Identical to §5c (the committed step-A policy) by construction. Suite: 3,757 passed / 8 xfail / 0 failed.

### 5h. THE LEVEL LANE — the census that re-ranked the plan, and the one mechanism that closes the holes (2026-09-04)

**The census (`reach_census.py`, `why_unreached.py`, `empty_census.py`, ladder, pass zero).** On the
ladder — not the test chromosome — the landed policy leaves HALF the pass-zero error at nodes no message
reaches at all:

| ladder row | whole-library \|err\| | unreached (SILENCE both sides) | one side | both sides |
|---|---|---|---|---|
| g50 ss.50 OFF (in scope) | 1,394,354 | **716,190 (51.4 %)** at 28,436 nodes | 398,348 (28.6 %) | 270,128 (19.4 %) |
| g50 ss.99 ON (in scope) | 331,677 | **197,026 (59.4 %)** at 30,743 nodes | 59,668 (18.0 %) | 74,933 (22.6 %) |

On the test chromosome the same census reads 0.2 % / 4.4 % unreached: its genes are apart and its
exons are whole, so it does not exercise the hole. **The mechanism is the EMPTY REGION.** 13,478 of the
ladder's 70,176 nodes have no total — 12,486 of them exon pieces (52 % of all exon nodes): 7,123 with
ZERO gDNA opportunity (a piece shorter than a fragment, cut by closely spaced termini and splice sites
inside one exon complex) and 5,363 dark (opportunity, no fragment). Every empty run has length ONE (a
region between two boundaries), and 5,403 of them sit between two exon|exon boundaries. Every rule's
licence asks the flank for a total or an opportunity (`a_g[o] > 0`, `n_u[i] > 0` — `probe_item5.py`),
so item 5 and the level rule both STOP at 2,697 same-strand terminus exon|exon boundaries, and the
boundaries on both sides of an empty piece hear nothing: `B exon|exon [TES+/TSS−/…]`, same-strand
flanks, orientation resolved, NO-RULE both faces — 24 % of the in-scope unstranded row's pass-zero error
by themselves, sitting at honest ignorance (½ of their count) with a true gDNA share near 0.1. The
stranded capture-ON row's unreached mass is the AMBIG (+−) exon complex — the both-stranded locus
(step H) — whose own channel is dead and whose only witnesses are the single-strand nodes at its ends.

**The design: a level is an ABSOLUTE quantity, so it needs no map and no knowledge of its recipient.**
`Message.level_gdna` becomes a PROFILE over `u = log(rho / rho_ref)` on the solve grid (the same `K`
and window as `lam`; `rho_ref` = the library's structurally pure gDNA density, Σ count / Σ opportunity
over intergenic regions — a coordinate choice, not a constant), carried with the last full node's
`(n, a)` and its class. Two coordinate changes, both at a node with a total `(n, a)`:

* **own composition → level:** `u(lam) = log(sigma(lam) · n / (a · rho_ref))`; above the node's total
  the level falls as the total's Poisson tail (`edge_level_row`'s form — the total bounds the level).
* **level → composition at the recipient:** the same map read backwards, widened by the HOP's price
  (the owner's rule 8: both totals' counting plus the abundance discrepancy beyond counting, per hop,
  nothing pooled), and made ONE-SIDED by the class pair — a level from an unprobed class (intergenic,
  intron, intron|exon boundary) into the exon class is a LOWER bound (step E's law: the interior may be
  enriched over its edge, never depleted); exon → unprobed an UPPER bound; the same class two-sided.

Three laws then close every hole at once: **(1) the level rule is the DEFAULT of every directed face
without a composition rule** (strand-change faces, termini both ways, the AMBIG complex, the inside of
termini, the edge — rule 5 and rule 8 become two instances of the one conversion, `edge_level_row` and
`level_map_lambda`/`level_row`/`level_bound_row` retire); **(2) an EMPTY node is transparent** — it holds
levels only and forwards them unchanged (a few base pairs of the same gDNA density); **(3) a full node
emits the PRODUCT of its own level and the priced level it holds** — forward-backward's rule — so the
intergenic region's whole count reaches the exon complex through the gene edge, hop by hop, each hop
priced at reception. A face with a composition rule sends composition only (the map already carries
the level), so no witness is counted twice. On unstranded data the honest content that reaches an exon
complex at pass zero is a lower bound from the nearest structurally pure gDNA and each node's own total
as an upper bound — a one-sided profile — and the point estimate a one-sided profile yields (the plateau's
median) is the ESTIMATOR's question, owed with the landscape's training population, not the message's:
the census is judged pass zero and full, halves apart.

**The first cut before the lane (`level_default_proto.py`): step A's level rule registered at every
boundary → region face without a rule** found ZERO such faces on the test chromosome and 3,585–4,907 on
the ladder's rows, moving them by −0.5…+1.5 % — because the faces that matter lead into EMPTY regions,
which that rule's licence (a recipient total) refuses. That is what forced the lane.

**v1 measured (ladder, full pipeline, `level_lane_proto.py` / `level_lane_arms.py`, 2026-09-05).**

| arm | g50 ss.99 ON (stranded, in scope) | g98 ss.50 OFF (unstranded, in scope) |
|---|---|---|
| `transfer` (landed) | 232,392 | 126,391 |
| lane v1: two-sided by class, products, every total's bound emitted | **310,211 (+33 %)** | 132,397 (+4.8 %) |
| no products | 265,720 | 129,519 |
| no bound claims | 230,239 | 128,321 |
| **every level LOWER-ONLY** (products on) | **221,545 (−4.7 %)** | **125,746 (−0.5 %)** |
| lower-only, no products, no bounds | 229,070 | 126,036 |

The harm sits exactly at the classes the lane newly reaches (`R exon (walled)` 43,190 → 61,136,
`B exon|exon [term]` 65,575 → 93,080, `[sj]` 48,303 → 78,199 on the ON row) and every bit of it is the
UPPER side: a total's bound from a low-total piece of an exon complex, multiplied along the chain, says
"gDNA is at most this" to a probed neighbour — darkness under capture read as absence, step E's refuted
upper side at the scale of a whole complex. **THE LAW, measured: a level that crosses a face says "at
least this much gDNA" and nothing more.** Lower-only levels win on the stranded capture-ON row (the
AMBIG complexes take lower bounds from the single-strand nodes at their ends) and hold everywhere else.

**The prize the law forgoes, priced.** On the in-scope UNSTRANDED row `g50 ss.50 OFF` at PASS ZERO the
two-sided lane reads 1,394,354 → **1,037,431 (−26 %)** — `B exon|exon [term]` 419,460 → 232,586, walled
exons 428,178 → 332,190 — and the lower-only lane reads 1,393,350 (nothing): under capture-OFF the
minimum total density of an exon complex truly bounds its gDNA from above, and a lower bound at an
RNA-rich node is dampened to nothing by rule 8's discrepancy price (the totals disagree by the RNA).
Full pipeline both read 144.3–145.1k against 144,571: the landscape repairs the first pass either way,
while training on it. So the first-pass gain on unstranded data is real, large, and gated on ONE fact
the message layer cannot see: whether this library's gDNA is enriched. That is not a message question —
it is the ENRICHMENT WITNESS (`ISSUES: two-sided-exon-row`), and there is a gDNA-specific one on
unstranded data: the exons of SILENT genes (no spliced fragment at any of their junctions at this depth)
carry pure gDNA, so their density against the intergenic density is the library's gDNA enrichment
spectrum, measurable at pass zero with no strand channel. Where that spectrum is flat, every level is
two-sided; where it is not, lower-only — a library-level fact, learned, never a gate (the owner's
ruling on capture). Parked behind the architecture with the landscape's training population.

**v2 refused; the landing form is v1 lower-only (2026-09-05).** Re-expressing rules 5 and 8 through the
lane (`level_lane_v2.py`) broke the zero-gDNA controls on the test chromosome — `g00 ss.50 ON` 22,299 →
26,594 (1.193×), `g00 ss.70 ON` 1.147× — because rule 8's total bound and its two-sided own-profile level
are what pull a terminus's inside to zero when the truth is zero, and they measured safe at that one
face; the culprit was rule 8's replacement, not the edge's (`lane_v2_edge` reads the same). A middle
form — an own measurement keeps both sides for its FIRST hop, forwarded levels lower-only (`lane_b`) —
held the g00 controls but harmed `g98 ss.70 ON` 1.023× with the rules kept and 1.060× with them replaced,
and on the ladder read 226,666 / 127,063 / 17,763 against lower-only's 221,545 / 125,746 / 17,777 on
`g50 ss.99 ON` / `g98 ss.50 OFF` / `g00 ss.99 ON`. The nearest-witness form (no products) measured equal
on the panels and 229,070 on the ON row. The intergenic REGION's whole count as a level made no
difference on any panel and was dropped (unannotated transcription contaminates it; the gene edge's
crossing is the structural source). So: lower-only, products, empties transparent, no intergenic
levels, rules 5 and 8 as landed. Landed in `src/` 2026-09-05, byte-identical to the prototype on 29/30
test conditions; the prototype's sixteen-row ladder: unstranded 7/8 at or below the landed policy
(worst 1.001×, best 0.972× at `g98 ss.50 ON`), stranded 6/8 (worst 1.024× at `g05 ss.99 OFF`, best
0.913× at `g98 ss.99 ON`), the four g00 rows identical.

**The residue the landing left, dissected (`lane_delta.py`, 2026-09-05): the RATCHET.** The landed lane's
one in-scope harm is `g05 ss.99 OFF` 44,714 → 45,798 (+1,084 fragments, 509 slots moved). 766 of them sit
at terminus exon|exon boundaries and the worst are a CHAIN of nine consecutive terminus boundaries in one
highly expressed complex — 4,200 crossings each, ONE true gDNA fragment each — every one moved 2 → 29.
The mechanism: each full node emits own × held, so nine soft one-sided claims (each node's own strand
mode is noise around zero, ±0.003 in share; its lower side is a half-nat rise) multiply along the chain
into a hard bound at the noisiest node's mode: the product of CENSORED likelihoods ratchets upward. A
two-sided product would converge to the truth; a lower-only product converges to the maximum of the
upward noise. Under the lane's own law a level is a BOUND, and two bounds on one density combine by
INTERSECTION — the pointwise minimum of two non-decreasing log-profiles, the tighter wins at each
density, nothing sharpens — not by product. `lane_isect.py` measured that form against the products form:
`g05 ss.99 OFF` 45,798 → 45,076 (two thirds of the harm recovered; 44,714 before the lane), `g50 ss.99 ON`
220,674 → 220,404, `g98 ss.50 OFF` 125,797 → 125,753, `g00 ss.99 ON` identical, and every row of the test
and sparse panels identical (no chain of full nodes exists there — which is also why the toy's gates
could not see the ratchet, so the landed form carries a hand-built three-node gate). LANDED: `emit` is the
intersection (`transfer_rows.intersect`, `lower_side`), and at the solve the two sides' levels intersect
before the constraint joins the composition evidence. The residual +0.8 % on `g05 ss.99 OFF` is the
tightest noisy witness of the chain — one node's 1-sd over-claim, no longer nine multiplied.

**THE LANDED FORM ON THE LADDER (`policy_benchmark.py --panel ladder --by-class`, 2026-09-05; the
pre-lane column is commit `2315b5af`'s policy).**

| ladder row | silent | relay | transfer before the lane | transfer WITH the lane | vs before |
|---|---|---|---|---|---|
| g00 ss.50 OFF | 1,254,145 | 58,840 | 195,268 | 195,268 | 1.000× |
| g00 ss.50 ON | 454,560 | 152,534 | 230,800 | 230,800 | 1.000× |
| g05 ss.50 OFF | 50,435 | 75,048 | 51,003 | 50,948 | 0.999× |
| g05 ss.50 ON | 518,535 | 349,481 | 225,825 | 223,499 | 0.990× |
| g50 ss.50 OFF | 155,660 | 189,133 | 144,571 | 144,288 | 0.998× |
| g50 ss.50 ON | 6,141,095 | 1,801,961 | 1,264,540 | 1,235,011 | 0.977× |
| g98 ss.50 OFF | 165,259 | 219,158 | 126,391 | 125,753 | 0.995× |
| g98 ss.50 ON | 12,030,888 | 2,433,908 | 2,214,300 | 2,152,433 | 0.972× |
| **unstranded** | | | **8/8 at or below the pre-lane policy**, worst 1.000× | below silence 7/8 (`g05 ss.50 OFF` 1.01×, as before) | |
| g00 ss.99 OFF | 30,606 | 7,421 | 17,814 | 17,814 | 1.000× |
| g00 ss.99 ON | 20,787 | 17,980 | 17,777 | 17,777 | 1.000× |
| g05 ss.99 OFF | 44,519 | 51,047 | 44,714 | 45,076 | 1.008× |
| g05 ss.99 ON | 85,294 | 97,957 | 81,407 | 80,416 | 0.988× |
| g50 ss.99 OFF | 125,634 | 155,496 | 116,580 | 116,165 | 0.996× |
| g50 ss.99 ON | 260,629 | 409,168 | 232,392 | 220,404 | 0.948× |
| g98 ss.99 OFF | 126,467 | 208,306 | 95,467 | 93,430 | 0.979× |
| g98 ss.99 ON | 298,597 | 456,838 | 221,971 | 201,578 | 0.908× |
| **stranded** | | | **7/8 at or below the pre-lane policy**, worst 1.008× (`g05 ss.99 OFF`, the tightest noisy witness) | below silence 7/8 | |

Pass zero on the two key rows (`pass0_score.py`): `g50 ss.50 OFF` 1,394,354 → 1,393,639 (a lower bound at
an RNA-rich node is priced to nothing — the unstranded first pass is the enrichment witness's, above);
`g50 ss.99 ON` 331,677 → **290,822 (−12 %)**, walled exons 90,343 → 77,167 at pass zero and 43,190 → 40,121
through the pipeline. Panels: the test chromosome, junction-probed and sparse-probed panels read
identity to 1.003× on every row but the sparse panel's `g00 ss.70 ON` (3,135 → 3,300: one gene type,
one node pair, a neighbour's own 3-sd strand error carried as a lower bound). Suite 3,764 passed / 8 xfail.

**THE TERMINUS-CLUSTER BLOCK'S FIRST READING (2026-09-05, the rebuilt test chromosome).** The block
mirrors MIR99AHG's ten transcript ends 126–147 bp into a shared last exon (twelve genes: `cluster` /
`capcluster` × ab / ba / eq × blocks 2 and 4); the index cuts the exon into nine pieces of 1–10 bp,
every one empty, on both strands. What it showed at once:

* **The census now sees the hole on the test chromosome**: 4.0 % of the pass-zero error on
  `g50 ss.50 OFF` (8,006 fragments at 124 terminus exon|exon boundaries) and 1.5 % on `g50 ss.99 ON`
  sit at nodes no message reaches — the cluster's INNER boundaries. Item 5's composition reaches the
  two outer ones from the outside exons and stops: a node that holds a COMPOSITION has nothing to send
  across a lane face, and its inside face into the first empty piece has no composition rule.
* **Completing the reach measured harmful and is REFUSED** (`lane_convert.py`: at a full node a held
  level is read as a composition for a composition rule and a held composition as a level for a lane
  face; rules into empty recipients dropped). The census reads 0 % unreached with it, and the rows read:
  test `g50 ss.50 OFF` 13,009 → 13,471 (+3.6 %), `g50 ss.99 ON` 7,690 → 7,637; ladder `g50 ss.50 OFF`
  144,288 → 148,516 (+2.9 %), `g50 ss.99 ON` 220,404 → 217,085 (−1.5 %), `g98 ss.50 OFF` +1.0 %,
  `g05 ss.99 OFF` 45,076 → 46,532 (+3.2 %), g00 identical. Dropping the rules into empties ALONE
  (`lane_dropdead.py`) reads identical on the test chromosome and +0.2…+5.4 % on the ladder: a
  composition already crosses a dark exon through the composition rules on both of its faces (rung 2 in,
  item 1 out — the maps read the boundaries' numbers, not the empty's), and that path is worth 5 % on
  `g98 ss.50 OFF`. So: THE LAW HOLDS FOR THE LANE TOO — what a node holds as a composition is never
  re-issued as a level (step A's law, measured). The information that WOULD reach the inner boundaries is
  a lower bound at an RNA-rich node with weak own evidence, and that is what harms: the reach is
  complete by construction, the message the lower-only law permits is not worth sending there.
* **The zero control's residue has a name now**: `g00 ss.99 ON` reads 64 → 216 on the rebuilt panel, 187
  of it at `capcluster_ab`'s nine inner boundaries, each at share 0.007–0.011 with 2,000 crossings and no
  gDNA at all — each boundary's OWN strand mode is noise around zero (±0.003), its lower side becomes a
  soft bound at its neighbours, and the tighter of two noisy neighbours wins: the NOISE RATCHET of
  lower-only levels among RNA-rich nodes of one density (`ISSUES: the-lower-bound-noise-ratchet`). A
  two-sided own-profile level would average the noise away — `lane_nb` measured −0.9 % on `g50 ss.99 ON`
  and +1.5 % on `g98 ss.50 OFF` against lower-only's −4.7 % / −0.5 % — so lower-only keeps more
  fragments and this residue.
* Standings on the rebuilt benign panel (new reads, every number moved; `policy_benchmark.py --panel
  test`): unstranded 16/20 rows below silence (worst 1.23× at `g50 ss.50 OFF`, the recorded one-sided
  first-pass weakness), stranded 7/10 (worst 3.38× at `g00 ss.99 ON`: 64 → 216, above). Cluster types at
  `g50 ss.99 ON`: `capcluster_eq` 699 → 584, `capcluster_ba` 280 → 404, unprobed clusters identical.

* **The SPARSE-probed panel shows the lane's limit, and it is the owner's call.** With one 125-bp probe
  centred on every annotated exon, the ten isoforms' 126–147-bp last exons put TEN overlapping probes on
  the cluster and none on the rest of the exon: an enrichment cliff INSIDE one exon complex. On
  `g50 ss.99 ON` (stranded, in scope) the row reads silent 10,100 / transfer 12,834 (1.27×), all of it at
  the three probed cluster genes (`capcluster_ab` 530 → 1,783, `_ba` 658 → 1,395, `_eq` 834 → 1,771;
  the lane switched off reads 10,025): the − gene's inner boundaries, true share 0.85, are pushed to
  0.92–0.95 by a lower bound from the ten-times-probed piece beside them — the interior is NOT enriched
  over that source, the bound is false, and the totals cannot show it (`lane_witness.py`: the strand
  witness carried as a scalar sees only the agreeing neighbour and recovers 1 %). Reading the arriving
  bound's WALL against the recipient's own strand mode (`lane_witness2.py`: the owner's rule, a widening
  only) recovers 57 % of it — 12,834 → 11,285 — at +0.8 % on the ladder's `g98 ss.99 ON` (201,578 →
  203,178) and +0.8 % on the sparse panel's, everything else identical: 1,550 fragments back on one
  adversarial row against 1,600 lost on one in-scope ladder row. NOT landed; recorded. The same physics
  refused step F under junction probes: no local witness prices capture's enrichment of one piece over
  its neighbour. The owner's domain call: whether probe designs that target isoform-specific ends are a
  case to protect, at that price, before the enrichment witness exists.

### 5i. THE ABLATION — what each level mechanism is worth, by node class (2026-09-05)

`level_ablation.py` switches one mechanism off at a time on the landed policy: the lane, rule 5 (the
edge's count reaches nothing by any path), rule 8 (nothing crosses a terminus's inside face), and all
three (composition rules only). Read on the ladder with `--by-class`, full pipeline; the delta is
"removed − landed", so a positive number is what the mechanism EARNS at that class.

| ladder row | lane | rule 5 (edge) | rule 8 (terminus inside) | all three |
|---|---|---|---|---|
| g50 ss.50 OFF (unstranded, in scope) | +283 (terminus boundaries +260) | **−661** (edge exons −569, licensed exons −470; introns +173, intron boundaries +218) | +80 | −297 |
| g98 ss.50 OFF (unstranded, in scope) | +638 (terminus boundaries +605) | +2,092 (edge exons +695, introns +498, intron boundaries +450) | −2 | +2,742 |
| g00 ss.50 OFF (zero control) | 0 | 0 | **+45,873** (walled exons +18,457, terminus boundaries +14,087, sj boundaries +7,091) | +45,873 |
| g50 ss.99 OFF (stranded, in scope) | +415 | +331 | −49 | +710 |
| g50 ss.99 ON (stranded, in scope) | **+11,988** (terminus boundaries +5,163, walled exons +3,069, sj boundaries +2,676) | +36 | −743 (walled −974) | +11,527 |
| g98 ss.99 ON (stranded, in scope) | **+20,393** (terminus boundaries +13,281, sj boundaries +3,318, walled +1,972) | +300 | +3,671 (intron boundaries +1,331, sj +1,257) | +24,929 |
| g05 ss.99 OFF (stranded, in scope) | −362 (walled −241: the noise ratchet) | −50 | −54 | −381 |

Pass zero on the two key rows: the levels together take `g50 ss.50 OFF` from 1,484,499 to 1,393,639
(−6 %, walled exons 480,727 → 428,294 by rule 8's pass-zero bound, licensed exons 214,762 → 196,309),
while the edge's own destinations read 4 % WORSE with it (edge exons 73,524 → 76,754: the plateau above
a lower bound); the lane takes `g50 ss.99 ON` from 331,784 to 290,822 (−12 %). On the test chromosome
(`ablation_test_pass0.out`) the signs on the unstranded capture-OFF rows flip — the edge costs 8 % of
`g50 ss.50 OFF` through the pipeline, rule 8 costs 16 % of `g00 ss.50 OFF` through the pipeline while
helping at pass zero — the recorded toy-versus-panel disagreement, and the landscape retrain is what
flips it. The ladder decides.

What it says. THE LANE is the message layer's largest single win: 5 % and 9 % of the two stranded
capture-ON rows, 12 % at pass zero, entirely at the classes the census named — the AMBIG and probed
exon complexes' walled exons and exon|exon boundaries, which held no message before; it earns nothing
on unstranded rows (a lower bound at an RNA-rich node is priced to nothing) and costs the noise ratchet
at low-gDNA stranded rows and the sparse cliff. RULE 8 is what holds the zero-gDNA control (19 % of
that row: its UPPER side pulls a terminus's inside to zero when the truth is zero — which is why it
cannot become lower-only) and is ± elsewhere. RULE 5 helps at high gDNA, where the lower bound sits
near the truth and the plateau above it is short, and HURTS at mid-gDNA unstranded rows, where the
plateau above the wall is what the estimate reads — the one-sided profile's point estimate, the same
mechanism as rung 2's licensed-exon weakness (`R exon (licensed)`: silent 7,480 / landed 10,196 on
`g50 ss.50 OFF`, 9,784 with no levels at all — mostly the composition face map's plateau).

**Where the wall is.** Not in the plumbing: every face has a rule, every node is reached or refuses a
message that measured harmful. It is informational. On unstranded data no local measurement says
whether a neighbour's gDNA density equals this node's (capture-OFF: uniform) or exceeds it (capture-ON:
a probe edge), so every message must be one-sided, and a one-sided profile's point estimate lands
above its wall. Three things, none a message: the estimator at nodes whose only evidence is a bound;
the ENRICHMENT WITNESS (silent-gene exons against the intergenic density on unstranded data; exon
against intron gDNA densities on stranded) that licenses two-sided levels where the library is not
enriched; and the landscape's training population, which today learns from the plateau.

**Order as it was planned.** v1 (`level_lane_proto.py`): the lane on the faces that have no rule today, the landed rules
untouched — the increment from reaching the unreached, on the ladder's four key rows, then the panels.
v2: rules 5 and 8 re-expressed through the lane (byte-identical or better), the old level rows deleted.
Then the RNA levels (`level_rna_pos/neg`) for the AMBIG complex by the same two conversions (step H).

## 6. The order of work (re-ranked 2026-09-06 by the AMBIG census; owner rulings 2026-09-06)

⭐ **Owner rulings, 2026-09-06.** (1) THE ENRICHMENT WITNESS IS THE gDNA LANDSCAPE PRIOR — the second
pass is where a level becomes two-sided where the library is not enriched; it is not a message rule.
(2) Nodes that cannot be solved — whose only evidence is a bound — DO NOT TRAIN the landscape. (3) The
message propagation rules are finished FIRST; the landscape (the estimator at bound-only nodes, the
training population) comes after.

**What remains of the rules, by measured mass** (`ambig_census.py`, the landed policy, full pipeline):

| case | ladder nodes | share of the remaining error | status |
|---|---|---|---|
| **H. THE AMBIG COMPLEX** — both strands admitted: the RNA levels per strand and the tilt | 9,912 (14 % of nodes) | **38 % of `g50 ss.99 ON`, 50 % of `g98 ss.99 ON`**, 16.5 % of `g50 ss.99 OFF`, 13.5 % of `g50 ss.50 OFF` — at the overlapping loci's exon\|exon boundaries and walled exons | gDNA reached by the lane (lower bound); the RNA split unimputed |
| **D. sj+terminus** — one boundary carrying a junction and a terminus | 386 | 0.5–0.7 % of any row | the lane serves the level; the composition part refused by `outside_flank` |
| strand-change faces, termini both ways, the empty piece, the walled classes, the dark exon between licensed faces | — | — | ✅ served (the lane; the composition rules across a dark exon) |

1. **H. THE AMBIG COMPLEX** — the last lane and the largest remaining case. ⭐ DESIGNED (`AMBIG_DESIGN.md`,
   owner-approved 2026-09-07/08; phase 0 built; phase 1 approved as written there). ⭐ THE DESIGN: `docs/dev/AMBIG_DESIGN.md` (2026-09-06). `level_rna_pos` /
   `level_rna_neg` as ABSOLUTE profiles over the log RNA density per strand (the gDNA lane's two
   conversions, with the strand's own opportunity `a_r`); a single-strand node's own strand profile
   gives its live strand's RNA density (`(1 − f) · total / a_r`), which continues across a face where
   that strand's population is unchanged; an AMBIG recipient converts the two arriving RNA levels and
   its gDNA level through its own total into a constraint on its TILT (its second degree of freedom),
   which its own strand mode alone cannot fix (the gDNA share cancels from the strand mean). The
   substrate: a BOTH-STRANDED block on the test chromosome mirroring a real ladder locus, as the
   cluster did (owner authors; the ladder's overlapping loci are the deciding measurement). Judged at
   its destinations, halves apart, pass zero beside the pipeline.
2. **D. sj+terminus** — small and self-contained: `outside_flank` resolves the orientation from the
   terminus flag with a junction present (the junction does not change which flank is inside); item
   5's outside map gains the junction's leaving flux as item 7's does; rule 8's inside level is
   unchanged. `instart` on the test chromosome is its structure.
3. **The ship protocol** (`MESSAGE_RUNGS.md`'s checklist): the default flip, `silent` / `relay`
   retired. The four g00 rows where the relay still leads (its two-sided reading of a dark edge) are
   the LANDSCAPE's to win under ruling (2): trained on intergenic and intron nodes at g00 it reads
   zero, and the second pass pulls every exon there.
4. **Then the landscape** (the owner's points 1 and 3): the estimator at bound-only nodes and the
   training population — the −26 % first pass on unstranded data, the cluster's inner boundaries,
   the noise ratchet, the sparse cliff all wait for it.

## 7. On "43–45 % of the error is the intron class"

Those shares are FRAGMENTS of misplaced gDNA (the by-class table sums |estimate − truth| in fragments),
so the intron's 65,242 at ladder `g50 ss.50 OFF` is absolute error mass, not a percentage of small
truths. Relative to each class's OWN mass the picture is different and more useful: exon|intron
boundaries 8.3 %, exon|intron termini 7.8 %, alternative splice sites 2.5 %, introns 2.3 %, terminus
boundaries 1.8 %, licensed exons 1.6 %, walled exons 1.0 % (capture-ON stranded: introns 6.4 % of a
small mass, every boundary class 3.0–3.4 %, exons 0.8–1.6 %). The intron's absolute error is largest
because its mass is largest; per fragment it solves as well as the exons. Whether absolute or relative
error is the one that matters depends on the deliverable: the library gDNA fraction sums fragments,
a per-transcript count cares about its own exons and boundaries. Both are reported from here on.
