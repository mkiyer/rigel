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

**F. Rule 3 made two-sided at the face.** The wall's position comes from the INTRON's gDNA level —
its density (hundreds of fragments) times its composition — carried by the opportunity ratio, not
from the 12–25-fragment crossing count; the crossing count becomes the pair's discrepancy WITNESS
(its excess over counting dampens), which is also where the capture taper shows. That removes the
ratchet. The plateau above the wall is the intron's own high side and is the intron's own solve to
sharpen (parked), unless a bounded taper marginal replaces the flat top; both are A/B'd after A–E.

**G. Rule 4's premise.** The per-pair discrepancy rule of rules 9–10 applied to rule 4 (and 3): where
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

## 6. THE ORDER, and the method every step keeps

1. **A** — rule 8 re-specified as the level rule, both terminus kinds. Fixes a misspecified rule and
   opens the level family with the case that has the most substrate.
2. **E** — rule 5 as a level, the two forms A/B'd. The zero controls.
3. **F + G** — rules 3 and 4 refined. The ratchet and the premise.
4. **B** — strand-change faces and termini both ways.
5. **C** — the empty chain pieces (grow `nest`).
6. **D** — junction plus terminus (grow `instart`).
7. **H** — the AMBIG ruling and the both-stranded locus.
8. The SHIP LIST: the 0.8.0-metric pricing under `transfer` (`calibration_vs_oracle.py`,
   `solvability_audit.py`), the flip, the obsolescence pass (`relay.py`, `variance.py`,
   `rna_anchor.py`, the instruments and tests that name the relay), `preflight --full`, goldens, docs.

Every step: stage 0 on certified truth → the simplest LOCAL form, no constant, nothing pooled →
prototype through `policy_prototype.py` (`--by-class`, all three probe panels, `pass0_score.py`'s
pass-zero beside the full pipeline, halves apart) → the ladder → `src/` with a fail-first gate per law
and every perturbation watched firing → the per-slot identity of the landed form against the
prototype. One rule per step; a ruling made mid-way triggers a sweep of every earlier rule carrying
the same premise.

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
