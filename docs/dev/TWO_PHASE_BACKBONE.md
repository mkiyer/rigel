# THE TWO-PHASE BACKBONE — an audit of the message layer against the owner's architecture, and the proposal (2026-09-04)

    ⚠ A DEV DOC — a working argument for one owner decision, written the day the owner asked
    "what are we publishing, and to whom?" Nothing here is settled. What settles MOVES: a ruling to
    `DESIGN.md` §6b, a lesson to `TRAPS.md`, an open problem to `ISSUES.md`.

## 0. The words, in plain terms (added after the owner asked, 2026-09-04)

* **A "row"** is my shorthand for what a message in the transfer policy actually contains: a CURVE
  over the possible gDNA shares of the destination — for every candidate share from 0 to 1, how well
  the sender's evidence supports it (a log-likelihood). It is called a row only because the code stores
  one such curve per slot as one row of an array (`PsiMessage.lam_rows`). The honest name is a
  **composition profile**, and the rename should say so.
* **Two-sided vs one-sided.** A two-sided profile has a peak: it penalises BOTH more gDNA and less
  gDNA than its best value, so it says "this much". A one-sided profile rises to a plateau and stays
  flat: it says "at least this much" and nothing about how much more. Today's message from a face into
  an exon is one-sided on the gDNA side (§4), because it is the intron's own profile pushed through the
  face's arithmetic, and the intron's density-against-background evidence cannot tell 90 % gDNA from
  100 %. When a pass forwards a one-sided profile through another face, each hop can only add another
  "at least", and the exon's answer drifts toward more gDNA — the measured harm on the unstranded rows.
* **The flux level.** The spliced fragments at a junction, divided by the junction's opportunity, give
  an ABSOLUTE rate of mature RNA (molecules per base). Using that rate to predict how many of the exon's
  own fragments are RNA — and calling the rest gDNA — is using a LEVEL. A composition (a fraction) is
  invariant to a common enrichment; a level is not. Under capture the junction and the exon interior are
  enriched by different amounts depending on where the probes sit, which the tool cannot see, so the
  level's prediction is off by that unknown factor (30× on the junction-probed panel).
* **The face map** is the arithmetic that converts the composition measured AT the boundary (the
  crossing fragments, spliced and unspliced) into the composition expected INSIDE the exon. It assumes
  every fragment at the face is enriched by the same amount. **The taper** is where that assumption
  bends: at a probed exon's edge a crossing fragment lies half in the intron and overlaps the probe
  partially, while a spliced fragment overlaps it fully, so the crossing is under-enriched relative to
  the spliced flux by a modest factor (measured 1.3× under exon probes, 1.0 under junction probes and
  off capture). **A bounded marginal over the taper** means: instead of assuming that factor is exactly
  1, let the face map average over the range of factors fragment geometry allows, so the exon's profile
  becomes honestly wider on the gDNA side instead of flat-topped.

## 1. The architecture, as the owner stated it — and what it commits us to

Two phases, and nothing straddles them.

**Phase 1 — `propagate`.** A forward pass, then a backward pass. At each hop a node RECEIVES a
message, reconciles it with its OWN claim, and passes the result on to the next node. Beliefs do
not change. When both passes end, every node holds TWO messages — one from each neighbour — and
the two nodes at a chain's ends hold one. A hop that carries nothing still arrives, as an
explicitly uninformative message (`DESIGN.md` §6b.11).

**Phase 2 — `solve`.** Every node solves once, from four things: its own evidence, the message
from its left, the message from its right, and the gDNA hyperprior.

This is sum-product belief propagation on a chain, and on a chain it is EXACT — no evidence is ever
counted twice — provided two laws hold at every hop:

* **What a node sends is its OWN CLAIM composed with what it holds from its FAR side** — never its
  belief (a belief already contains the prior and the neighbours, so sending it would echo). The
  message `s → i` is `map_{s→i}( own_s ⊗ held[s ← far side of s] )`.
* **No echo.** What `s` sends toward `i` never includes what `s` holds FROM `i`. The two passes
  give this for free: the forward pass only ever composes left-side arrivals, the backward pass
  right-side ones.

Everything else in the policy — the licence, the map, the width, the discrepancy rule — is the
RECIPIENT's decision at the moment of receiving: **STOP** (composition cannot cross: deliver
Silence), **FORWARD** (the shared-population identity), or **MODIFY** (the face's map with its
counting width and the pair's own discrepancy). That is the owner's ruling made structural: the
three decisions ARE the map registry, keyed by directed face.

## 2. The audit — what the code does today, piece by piece

| piece | what it does | verdict against §1 |
|---|---|---|
| `messages/foundation.py` (ratified 2026-08-26) | The two timepoints exactly (`PropagationModel.propagate(own, incoming, hop)`, `SolveModel.solve(own, forward, backward)`), four laws the skeleton enforces, and a `Message` of five Gaussian `Claim(abundance, precision, measured)` LANES — per-population LEVELS | The SKELETON is the owner's architecture. The MESSAGE TYPE is the pre-rebuild currency: the rebuild ruled (2026-09-01) that every message is a COMPOSITION on the λ grid, because a level cannot cross a capture cliff. `MessagePolicy` runs on it and is byte-identical to silence — a scaffold nobody builds on |
| `messages/transfer.py` `prepare()` | Builds all ten messages BY HAND at prepare time: for each face it maps the source's own row and adds it into `rows[destination]`. Two of the ten are hand-built TWO-hop messages — rung 2 (intron → boundary → exon) and item 5's composed transport (the exon's arrivals → the terminus boundary) — with a hand-kept no-echo (`rows5` snapshot) and a `consumed` set so nothing is counted twice | **The propagation phase collapsed into `prepare`.** Correct for one hop; it is why every second hop needs bookkeeping |
| `transfer.py` `_Ledger` + `scan()` | A forwarding pass layered on top: arrivals with NAMES, `maps`, `consumed`, `depth`, the integer `hops` budget; `publish()` returns `sent[s] = held[i]` so the backbone's gather-at-source yields the destination's row | Scaffolding that compensates for the row above. Under proper forward-backward the ledger, the names, the consumed set and the budget all disappear: what a node holds from its far side IS what it forwards |
| `sweep.solve_chain` | `relay = policy.prepare(ctx)`; `_scan` ×2 calling `step(s, i)`; `publish()` gathered at the SOURCE (`_at_source`) into `NeighbourState`; `relay.deliver(left, right)` | The LOOP is right — two passes, chain order, then one solve. The VOCABULARY is wrong: `relay` names both the protocol and the frozen policy; `publish` says nothing about who receives; `deliver` is phase 2's first half under a phase-1 word. The gather-at-source exists for the relay policy (which publishes belief-valued state and transports it at solve time), not for the contract |
| `NeighbourState` (source-indexed) | The structural protection against `TRAPS: a-message-from-the-destinations-belief` | Right law, wrong place. The protection belongs on the KERNEL's signature: a recipient's decision receives the source's outgoing message and the hop, never a belief of the destination |
| `silent.py` | `scan → None`, `deliver → silent` | Fine; it becomes `receive → Silence`, `solve → silent` |
| `relay.py` (frozen, shipped) | Scalar `step(s, i)` fusing into a running state at `i`; `publish` the state; `deliver` re-transports the SOURCE's state (the "vectorised twin") | Must keep running byte-identically until the flip. It does, under the proposal — see §3 |

**The one-line diagnosis.** The transfer policy has the right messages and the wrong PHASE
STRUCTURE: it computes messages where the architecture says it should compute CLAIMS, and then needs
a ledger to do what a pass does by construction. The foundation spec has the right phase structure
and the wrong CURRENCY. Neither is the owner's architecture; each holds half of it.

## 3. The proposal — one skeleton, the owner's names

```
prepared  = policy.prepare(ctx)                 # per sweep: every node's OWN CLAIM + the face rules
from_left = backbone.pass(prepared.propagate(backward=False))   # every node holds what its LEFT sent
from_right= backbone.pass(prepared.propagate(backward=True))    # ... and its RIGHT
evidence  = prepared.solve(from_left, from_right)               # phase 2, the policy's half -> PsiMessage
belief    = psi(own, evidence, prior)                           # phase 2, the backbone's half
```

**The kernel.** `propagate(backward)` returns `receive(source, destination) -> Message`. The
backbone owns the pass: in chain order, for every node `i` with a neighbour `s` on that side,
`held[i] = receive(s, i)`. Inside `receive` the policy composes what `s` sends — `own_s ⊗ held[s]`
(what `s` holds from ITS far side, which the same pass wrote one step earlier) — and applies the
recipient's decision for the face `s → i`. The backbone hands the kernel two INDICES and nothing
else; the policy reads `held[s]` from its own pass state. A kernel that wanted the destination's
belief has nowhere to get it.

**The laws the backbone enforces** (each a gate, each watched firing against a broken build):

1. after each pass, every node with a neighbour on that side holds a `Message` (never `None`); a
   node without a neighbour holds `NoNeighbour`, which is not the same thing as `Silence`;
2. both passes run in chain order, and the forward pass reads only `left`, the backward only `right`;
3. the kernel is called with indices only; `held` is written by the backbone, never by the policy
   (a policy cannot reach past its hop);
4. `solve` receives the two held arrays indexed AT THE RECIPIENT and returns a `PsiMessage` whose
   channels pass the existing domain assertions.

**The message.** One type, two LANES by what they can cross:

```
Message(composition: row over lam | None,     # scale-free: the rebuild's currency
        level: (log gDNA rate, log-variance) | None)   # counts/bp: crosses where composition cannot
Silence = Message(None, None)
```

The composition lane is today's row. The LEVEL lane is the checklist's owed "level message": gDNA is
genomically continuous, so its rate crosses ANY face (a terminus, a strand change, an empty chain
piece), priced by the opportunity ratio and the enrichment premise; rung 3's edge bound is already its
first instance (an edge's gDNA count converted at the exon, one-sided, one hop). A level is converted
into a composition row only AT A RECIPIENT, using the recipient's own count — an observation, which
the contract allows at either end of a hop.

**The recipient's three decisions, per lane, per directed face**, are one registry:
`face[(s, i)] = STOP | FORWARD | MODIFY(map, width, discrepancy)`. Today's policy already holds every
map; what it lacks is the pass that composes them.

**Compatibility with the frozen policies, without an adapter.** `SilentPolicy`: `receive` returns
`Silence`, `solve` returns silence. `RelayPolicy`: its scalar `step(s, i)` runs unchanged inside
`receive` (updating its own running state at `i`) and RETURNS the source's state tuple as the
message — which is its semantics ("what you hold from `s` is `s`'s published state") — so its `solve`
is today's `deliver` on the stacked tuples, byte-identical, and `_at_source` leaves the backbone.
Gated by the goldens and `arm_identity.py`.

## 4. What the forward-backward form CHANGES in what is delivered — the honest list

At every face the sender composes `own_s ⊗ held_far[s]`. Against today's ten messages:

| change | today | under forward-backward |
|---|---|---|
| the boundary's own strand row toward the EXON | not delivered (own_b goes to the intron only, item 2) | travels through the face map with the intron's row — a new witness of the same crossing composition |
| forwarding through an exon | only into a terminus boundary (item 5's composed transport) | through every face the exon has: the scan's measured hops, in both directions |
| forwarding through a boundary into an intron | never | the exon's splice-out row reaches the intron (the "first real hop" gate already measures this) |
| rung 2 and the composed transport | hand-built two-hop messages with a `consumed` set | two ordinary hops; no bookkeeping |
| chain ends and refused faces | no message | `NoNeighbour` and `Silence`, explicit |

**The zero point is therefore stated per family, not as one byte-identity**: with forwarding OFF the
pass reproduces every one-hop message exactly (gate), and rung 2 and the composed transport are
reproduced as the composition of their two hops (gate); the three rows above are the CASE under
study, measured `--by-class` at their destinations, halves apart.

**Why forwarding harmed the in-scope unstranded row (`DESIGN.md` §6b.10) — restated in these
terms.** An exon's held row from a licensed face is a LOWER BOUND (flat above the face map's
ceiling); an unstranded exon has no own claim, so it forwards the bound unchanged into its other
face, where the next map's blur cannot make it two-sided. A pass compounds whatever shape it is
given. The prerequisite is a two-sided row at exons — §5.

## 5. The certified-flux message into exons — the derivation, in the composition currency

⛔ **REFUTED on 2026-09-04 by the adversarial probe panels and the ladder's zero control — see §6d and
`ISSUES: the-certified-flux-row-as-a-level`. This section stays as the derivation record: the row is a
LEVEL crossing a face, and its transport is set by probe placement the tool cannot see.**

**What an exon's RNA is made of.** Fragments contained in exon `e` that are RNA came in by exactly
two routes: CONTIGUOUSLY across a face (the same molecules the face's unspliced crossing measures,
carried by the opportunity ratio `e_r(e) / a_r(b)`), or SPLICED IN through the face's junctions
(measured at the face as the route-summed flux rate `r_b`, certified RNA, opportunity `e_r(e)`).
The exon's gDNA is whatever remains of its own count `C_e`.

**Today's face map takes the other road.** `face_map_lambda` computes the exon's gDNA from the
CROSSING's gDNA level scaled by the opportunity ratio, and the exon's RNA from the crossing's RNA plus
the flux. Under capture the exon interior is more enriched than its edge crossing by an unknown step
`s ≥ 1`, so the row's top is left FLAT to tolerate it — which is exactly the one-sidedness §4 names.
Every candidate the thread record measured (`COMPOSITION_TRANSFER_STAGE01.md`, THE SCAN) sits on
this road: `s` free (today, flat top), `s` bounded by the totals (the abundance-bounded form; fails
capture-ON where RNA dominates the total), `s = 1` (the two-sided Poisson form; catastrophic
capture-ON).

**The road that has no `s` on it.** Write the exon's RNA as a prediction and let the exon's own count
supply the total:

    RNA_e(lam_u)  =  n_u (1 - sigma(lam_u)) / a_r(b) * e_r(e)  +  r_b * e_r(e)
    C_e ~ Poisson( gDNA_e + RNA_e ),   gDNA_e = f_e * C_e

so for a hypothesised exon composition `f_e` the implied RNA count is `(1 - f_e) C_e`, and the row is
its likelihood under the prediction: the Poisson–Gamma (NegBinomial) marginal over the flux rate's
posterior `Gamma(flux_b + 1/2, A_b)`, averaged over the crossing's composition row (`own_b ⊗ own_i`,
what `b` holds and sends). The enrichment step never appears — it only ever scaled the gDNA, and the
gDNA is now the remainder. The row is TWO-SIDED because `C_e` is a two-sided measurement and the
prediction is a two-sided likelihood. This is the relay's anchor (`rna_anchor._quadrature_rows`:
route sum, NB marginal, the unspliced-crossing RNA as the "nascent" nodes) re-derived as the face's
MODIFY rule, with the crossing's own row in place of the intron's excess-over-background posterior,
and it REPLACES rung 2's transported row rather than joining it — the two share the flux and the
intron count, so fusing them would count both twice.

**What it uses, and from which end.** `C_e`, `e_r(e)` — the destination's OBSERVATION and CONSTANT
(allowed). `n_u`, `a_r(b)`, `r_b`, `flux_b`, the crossing's row — the source's. No belief of `e`.

**The width, local.** Counting: the NB size `flux_b + 1/2` and the crossing's `trigamma(n_u + 1/2)`.
The premise: at an exon with TWO certified faces the two flank rates are two witnesses of one
quantity; where they disagree beyond counting, widen this exon's rows by the excess —
`max(0, d² − (trigamma(flux_L+½) + trigamma(flux_R+½)))`, `d = log(r_L / r_R)` — nothing pooled
(the relay pooled it library-wide as `route_pair_log_variance`; the owner's rule is per pair).
A single-face exon carries counting width alone; the pooled centre fit and MAD the relay adds are
what the local form deliberately does not have, and stage 0 says whether they are needed.

**Why the "rough form" claimed 19 % gDNA at the zero control.** It took `r_b e_r(e)` as the exon's
WHOLE RNA — dropping the contiguous route. On the panel the contiguous share of RNA is ~20 %
(the sparse-nascent design, `DESIGN.md` §0b), so a zero-gDNA exon reads ~20 % gDNA. Stage 0 (§6)
tests this attribution directly, against the alternative the handoff suggested (an overhang mismatch
in the sj opportunity): the context's route rate is `Σ flux_J / crossing_eff_length(reach)` — the
same opportunity the anchor's route table uses.

## 6. Stage 0 — the flux-implied RNA against certified truth (`flux_stage0.py`, session scratchpad)

Seven conditions, exons with at least one certified interface, `log(truth / prediction)` where the
prediction is the face's route-summed rate × the exon's RNA opportunity (the context's `route_rate`,
the same crossing opportunity the anchor's route table uses):

| condition | exons (1 face / 2) | (1) log(spliced-in RNA / pred), depth ≥ 50 | (4) two-face excess var beyond counting (median / mean) |
|---|---|---|---|
| g00 ss.50 OFF | 354 / 44 | −0.031 (mad 0.074) | 0.000 / 0.007 |
| g50 ss.50 OFF | 356 / 43 | −0.038 (mad 0.074) | 0.000 / 0.065 |
| g98 ss.50 OFF | 276 / 33 | −0.045 at depth 10–50 (mad 0.25; 174 of 276 below depth 10) | 0.000 / 0.150 |
| g00 ss.50 ON, probed | 178 / 22 | −0.005 (mad 0.050) | 0.000 / 0.030 |
| g05 ss.99 ON, probed | 178 / 22 | −0.021 (mad 0.057) | 0.000 / 0.021 |
| g50 ss.99 ON, probed | 178 / 22 | −0.020 (mad 0.067) | 0.000 / 0.001 |
| g98 ss.99 ON | 157 / 19 | +0.029 (mad 0.086) | 0.000 / 0.275 |

**What it settles.** (a) At a CERTIFIED interface the flux-implied RNA is an unbiased predictor of the
exon's spliced-in RNA — a 2–4 % over-prediction, inside counting — on every condition, probed or not,
on and off capture: the rate TRANSFERS under capture (a mature molecule's probe continues in transcript
space). (b) The overhang hypothesis is refuted: the route rate is not under-reading. (c) The "19 % at
the zero control" of the rough form is NOT the contiguous route on this substrate (1.2 % of exon RNA
here; the ladder's is 20 %): a rough-form replay at certified interfaces claims 1.8 % (OFF) / 2.1 %
(ON) where truth is 0. The refuted arm was delivered at LICENSED faces, not CERTIFIED interfaces —
exons with routes the index never sees (the shadow transcripts) read their unseen RNA as gDNA. The
message is gated on the structural claim, as the anchor is. (d) The two certified faces of one exon
agree within counting everywhere (median excess variance 0.000; the share of pairs beyond 2σ is
0.02–0.14), so the LOCAL discrepancy rule is a no-op on this substrate and the pooled centre fit and
MAD the relay carries are not needed — the simplest local form is the measured form. (e) Off the
main road: the crossing's unspliced RNA carried by the opportunity ratio predicts the exon's contiguous
RNA to +0.25…+0.3 nats LOW under capture at probed exons (the edge taper the handoff's last lesson
names) and within noise off capture — the ladder's 20 % contiguous share will see this; recorded, not
built.

## 6b. The prototype, three arms of one mechanism (`flux_proto.py`, through `policy_prototype.py --by-class`)

`flux_join` adds the flux row beside rung 2's row (the relay's placement); `flux_replace` puts it in
rung 2's place at certified interfaces (rung 2 stays at licensed-but-uncertified faces);
`flux_trunc` adds the s ≥ 1 law as a hard truncation of the crossing marginal. Whole-library |error|
in fragments, silent / transfer / join / replace:

| condition | silent | transfer | join | replace |
|---|---|---|---|---|
| g00 ss.50 OFF (zero control, unstranded) | 43,412 | 16,819 | 15,863 | 19,024 |
| g50 ss.50 OFF (in scope, unstranded) | 9,349 | 11,242 | 9,960 | **8,876** |
| g98 ss.50 OFF | 6,620 | 6,363 | 5,846 | **5,539** |
| g50 ss.99 OFF | 8,706 | 8,381 | 8,355 | **8,231** |
| g05 ss.99 ON | 2,371 | 2,293 | 2,313 | **2,277** |
| g50 ss.99 ON | 7,186 | 6,367 | 6,038 | **6,006** |
| g98 ss.99 ON | 8,340 | 6,854 | 4,492 | **4,466** |
| g00 ss.50 ON (deferred) | 40,315 | 26,860 | 25,382 | 25,572 |

At the destinations (`R exon (licensed)`): g50 ss.50 OFF 1,667 silent / 4,071 transfer / 1,803
replace — the in-scope unstranded row goes below silence for the first time; g98 ss.99 ON 4,674 /
4,519 / 2,132. `flux_trunc` is REFUTED (5,712 at g50 ss.50 OFF): a hard cut on a crossing count of
9–25 fragments is a bias, not a law.

**The one loss, understood.** At the unstranded zero control replace is worse than transfer at licensed
exons (6,615 → 8,672) and join is better (6,039). Rung 2's row is legitimately SHARPER there: the
intron's thousands of fragments certify the crossing has no gDNA, and through the equality map
"no gDNA at the crossing" is "no gDNA in the exon" for ANY enrichment (s · 0 = 0), while the flux
row's width is the flux's own counting. The two rows carry different information — the RNA pin
(two-sided, from the flux and the exon's count) and the crossing's gDNA carried by a step s ≥ 1 —
and the exact single message is their JOINT likelihood over (λ_u, s) with a premise on s, which is
item 6's abundance-discrepancy structure applied to this face. That is a second mechanism and is the
NEXT case; join's number is its upper bound at the zero control. Joining as-is double-counts the flux
and the intron row (the relay did exactly this).

## 6c. The forward-backward skeleton, prototyped (`bp_proto.py`)

The transfer policy re-expressed as §3: every node's own claim, a rule per directed face (absent =
STOP, identity = FORWARD, a map = MODIFY), the backbone's pass written as it will be in `sweep.py`
(`held[i] = receive(s, i)` in chain order, both directions), `solve` fusing the two held messages.
The edge bound is the LEVEL lane's own claim (the edge's count; the exon's rule converts it), so it
forwards like everything else. `bp_today` restricts composition to the two faces today forwards
across and sends parts apart, as today's `prepare` does; **it reproduces the shipped transfer per slot
to the last bit on 1,768 of 1,782 slots** (unstranded and stranded conditions). The 14 that differ are
all exon|exon terminus boundaries in the walled block's CHAINS of termini (Σ|Δ| 0.2 and 7.9
fragments): today's `prepare` forwards a snapshot taken before items 6–7 delivered, a pass forwards
what the node holds — today's construction ORDER was an implicit rule, and the pass replaces it.


## 6d. THE REFUTATION — the row is a LEVEL, and probe placement decides its transport (measured after §6b)

Three probe worlds and the ladder, `policy_prototype.py --all`, the halves apart (`halves.py`):

| panel | half | `flux_replace` vs transfer | ≤ silent | worst |
|---|---|---|---|---|
| test (exon-probed) | unstranded 10 | wins 9 | 10/10 | 1.13× (`g00 ss.50 OFF`) |
| test | stranded 20 | wins 16 | 19/20 | 1.012× |
| test_junction | unstranded 10 | wins 6 | 10/10 | **1.48×** (`g98 ss.50 ON`) |
| test_junction | stranded 20 | wins 7 | 12/18 | **4.08×** (`g50 ss.99 ON`: 6,502 → 26,497) |
| test_sparse | unstranded 10 | wins 5 | 8/10 | **5.04×** (`g05 ss.50 ON`) |
| test_sparse | stranded 20 | wins 8 | 11/20 | **17.4×** (`g25 ss.70 ON`) |
| ladder `g00 ss.50 OFF` | — | 299,380 → **995,880** | — | walled exons 106k → 371k, terminus boundaries 93k → 296k |
| ladder (16) | unstranded 8 | wins 4 | 5/8 | **6.86×** (`g00 ss.50 ON` 342,666 → 2,350,977); `g05 ss.50 ON` 4.09×; `g50`/`g98 ss.50 ON` 0.86–0.88× |
| ladder | stranded 8 | wins 4 | 7/8 | 1.32× (`g00 ss.99 OFF`), 1.23× (`g00 ss.99 ON`); every other row within 2 % |

**Stage 0 on the junction panel names the mechanism**: at probed exons the flux over-predicts the exon's
spliced-in RNA by e^3.0…e^3.6 (mad 0.04 at two-face exons) — the spliced fragments sit on the junction
probe, the exon's contained fragments do not. The rate does NOT transfer; it transferred on the
exon-probed panel because a mature molecule's probe continues in transcript space there. The relay's
transport-centre fit (`rna_anchor.left_fit_center_spread`) REFUSES on 9 of 10 conditions tried and
misfires where it accepts (`g50 ss.99 ON` 6,006 → 32,235): the transport is bimodal by probe state
(−3.5 probed, −0.2 unprobed) and no single number fits it. On the ladder's zero control the mechanism
is different and worse: tiny exons (median 3 contained fragments) beside a large flux, where a
count-based row has no resolution and its false modes train the refit prior
(`ISSUES: gdna-landscape-trains-on-false-positives`).

**The conclusion, and it is §6b.9's founding refusal re-measured.** A LEVEL carried between locales
under capture is refuted by probe placement alone, and the tool never sees the probe panel. The
composition paradigm is invariant to a COMMON enrichment step; its only residual at an intron|exon face
is the TAPER — the spliced fragments' enrichment against the crossing fragments' at the probe edge
(+0.25…+0.3 nats under exon probes, ≈0 under junction probes, 0 off capture; stage 0's check (3)) —
which the flat top tolerates crudely. So the two-sided exon row cannot come from the flux level; it can
only come from (a) a sharper crossing composition — the intron's own solve, `ROADMAP.md` rank 3 — and
(b) a bounded marginal over the taper in place of the flat top. Recorded as
`ISSUES: the-certified-flux-row-as-a-level` (refused) and `ISSUES: two-sided-exon-row` (open).

## 6e. The formal forward-backward form, measured with today's rows (`bp`)

| panel | half | `bp` vs transfer | ≤ silent | worst |
|---|---|---|---|---|
| test | unstranded 10 | wins 6 | 8/10 | 1.53× (`g05 ss.50 ON`) |
| test | stranded 20 | wins 13 | 15/20 | 1.037× |
| test_junction | unstranded 10 | wins 6 | 8/10 | 1.13× |
| test_junction | stranded 20 | wins 14 | 15/20 | 1.037× |
| test_sparse | unstranded 10 | wins 5 | 7/10 | 1.32× |
| test_sparse | stranded 20 | wins 9 | 14/20 | 1.08× |
| **ladder (16)** | **unstranded 8** | **wins 7** (the eighth a tie, `g05 ss.50 OFF` 1.000×) | 7/8 | zero control 299,380 → 233,420 (0.78×); `g50 ss.50 ON` 0.81×, `g98 ss.50 ON` 0.82×, `g98 ss.50 OFF` 0.91× |
| **ladder** | **stranded 8** | **wins 7** | 7/8 | 1.006× (`g05 ss.99 OFF`); `g98 ss.99 ON` 0.84×, `g98 ss.99 OFF` 0.91× |

The stranded half meets the minimal-harm bar on every panel; the unstranded half is mixed, and the
in-scope row (`g50 ss.50 OFF`, 11,242 → 12,328) worsens exactly as §6b.10 predicted: one-sided rows
compound when forwarded. The zero control improves (16,819 → 15,352). On the LADDER — the shipping judgement — the formal form with today's rows WINS BOTH HALVES: 7/8 and
7/8 against the shipped transfer, below silence 7/8 and 7/8, the walled exons at the zero control
106,083 → 80,480 (the class the scan exists for) and licensed exons at `g00 ss.50 ON` 44,842 → 39,575.
The test chromosome and the ladder disagree in rank on the unstranded half
(`TRAPS: a-toy-and-a-panel-can-disagree-in-rank`); the ladder decides. So the FORM is right, measured,
and ladder-positive WITHOUT the two-sided row; the two-sided row remains what would let the unstranded
in-scope test-chromosome row stop worsening, and is the open problem above. Full tables:
`bp_ladder.out` in the session scratchpad.

## 6f. PASS ZERO beside the full pipeline — the owner's method applied (`pass0_score.py`, 2026-09-04)

`calib_refit_iters = 0` (the first-pass solve, which is what trains the landscape) beside the shipped
pipeline; whole-library |error| and the licensed-exon class, silent / transfer / `bp` / `flux_replace`:

| condition | frame | silent | transfer | `bp` | `flux_replace` | licensed exons (silent / transfer / bp / flux) |
|---|---|---|---|---|---|---|
| g00 ss.50 OFF | pass 0 | 297,218 | 278,108 | 277,240 | 266,954 | 256,664 / 238,616 / 237,761 / 227,462 |
| | full | 43,412 | 16,819 | 15,352 | 19,024 | 29,323 / 6,615 / 5,495 / 8,672 |
| g50 ss.50 OFF (in scope) | pass 0 | 156,001 | 169,133 | 171,788 | **29,830** | 130,867 / 145,432 / 146,541 / **6,128** |
| | full | 9,349 | 11,242 | 12,328 | 8,876 | 1,667 / 4,071 / 4,983 / 1,803 |
| g98 ss.50 OFF | pass 0 | 22,094 | 12,411 | 11,308 | 7,801 | 10,603 / 6,976 / 7,033 / 2,367 |
| g05 ss.50 ON (deferred) | pass 0 | 294,995 | 305,021 | 318,980 | **46,215** | 257,358 / 267,122 / 280,927 / **8,315** |
| g50 ss.99 OFF | pass 0 | 10,433 | 9,305 | 9,264 | 9,133 | 2,744 / 2,373 / 2,379 / 2,201 |
| g50 ss.99 ON | pass 0 | 9,096 | 7,408 | 6,691 | 6,631 | 4,912 / 4,873 / 4,520 / 4,096 |
| g98 ss.99 ON | pass 0 | 9,697 | 7,478 | 5,054 | 4,883 | 5,081 / 4,842 / 3,394 / 2,247 |
| g05 ss.99 ON | pass 0 | 4,217 | 3,418 | 2,710 | 3,326 | 2,641 / 2,201 / 1,839 / 2,110 |

**What pass zero says that the full pipeline hid.** (1) On unstranded rows today's messages do not solve
the exons AT ALL in the first pass — transfer and the formal form sit at or above silence at licensed
exons (131k → 145k → 147k on the in-scope row); every unstranded exon number in the full pipeline is the
PRIOR's work, and the messages' effect there is the harm the prior sees. (2) The two-sided flux
profile is the only thing measured so far that solves unstranded exons in the first pass — 20–30×
lower error where the transport is 1 — which is exactly the population the owner wants the landscape
trained on; its refutation (§6d) is the transport under other probe designs, not its shape. (3) On the
stranded half the formal form's wins are first-pass wins (`g98 ss.99 ON` 7,478 → 5,054), so they are
real and not the prior's. (4) The test-chromosome zero-control loss of the flux row was PRIOR-MEDIATED:
at pass zero it is slightly better than transfer there (266,954 vs 278,108). The ladder's pass-zero rows
(`pass0_ladder.out`): the formal form `bp` beats the shipped transfer at pass zero on all six rows tried
(`g00 ss.50 OFF` 2,711,072 → 2,532,590; `g50 ss.50 OFF` 1,537,644 → 1,494,766; `g98 ss.50 ON` 4,463,243
→ 3,947,880; `g05 ss.50 ON` 3,238,391 → 3,208,461; `g50 ss.99 ON` 337,883 → 331,882; `g98 ss.99 OFF`
131,086 → 118,436), so its ladder win (§6e) is a first-pass win and not the prior's.

**So the exon question, restated under the owner's method:** the first pass needs a two-sided exon
profile on unstranded data, and the only candidate that delivers one carries a transport factor set by
probe placement. Learning that factor per library — from the exons themselves, with a reason every
exon shares (one probe design per library) — is the open problem, `ISSUES: two-sided-exon-row`, and
it is taken up after the architecture works end to end (the owner's order).

## 7. The decisions

**Ruled by the owner on 2026-09-04 (`DESIGN.md` §6b.12):** `re-found-or-keep-two-homes` — RETIRE the
foundation's Gaussian-lane `Message` and the `message` policy now, re-found the skeleton as §3 designs
it; `the-names` — accepted (a wider rename campaign waits behind an accurate, working policy; the owner
keeps further items in personal notes); and the METHOD: every policy idea is judged at PASS ZERO and
with the prior, apart, because the landscape prior is trained on the first pass and returns its errors.
⚠ §6d's ladder zero-control number was read through the full pipeline only; its pass-zero reading is
`pass0_score.py` (session scratchpad) and is recorded in §6f.

**Already answered by standing rulings, so not re-asked:** the currency is the composition row (the
rebuild, §6b.9); two phases and the recipient decides (foundation + §6b.11); every node holds two
messages (§6b.11); the simplest LOCAL form of the discrepancy rule, nothing pooled
(`ISSUES: the-pooled-hop-step`); the relay stays frozen and byte-identical until the flip.

**Open, and the owner's:**

* **`re-found-or-keep-two-homes`.** Retire `foundation.py`'s Gaussian-lane `Message` and
  `MessagePolicy` (the `message` policy name) NOW, re-founding the skeleton on the row currency in
  `messages/__init__.py`; or keep them until the flip. Recommendation: retire now — they are a
  scaffold byte-identical to silence, and two homes for one architecture is the divergence the MOVE
  RULE exists to prevent. The four laws they enforce survive as the backbone's gates.
* **`the-names`.** `prepare / propagate / receive / solve`; `Message`, `Silence`, `NoNeighbour`;
  `Prepared` for the per-sweep object (retiring `Relay` the protocol, which collides with
  `RelayPolicy`); `PsiMessage` stays (it is what ψ receives) unless the owner prefers `Evidence`.
* **`skeleton-first`.** Land the skeleton first — byte-identical for silent and relay, per-family
  identical for transfer (the walled block's terminus chains excepted, where the pass is the rule),
  the four laws gated — then the cases. Recommendation: skeleton first; nothing the skeleton
  delivers changes until a case is ruled.
* **`the-two-sided-rows-source`** — re-posed in plain words after the owner asked (§0). On
  unstranded data an exon's answer comes from its faces. What its faces send today says "at least this
  much gDNA" (one-sided), and forwarding that through more faces makes the exon drift toward gDNA. The
  fix that used the junction's absolute RNA rate is refuted (§6d). Three ways remain, and the question
  is which one, and when: (a) sharpen the INTRON's own answer — the face only relays the intron's
  composition, so a two-sided intron profile gives a two-sided exon profile (the intron's own solve,
  parked at `ROADMAP.md` rank 3); (b) make the face map honest about the taper (§0) so the exon's
  profile is wide rather than flat; (c) neither yet — the formal form with today's messages already wins
  both halves of the ladder (§6e), so land it and take the exon profile as its own later case.
  Recommendation, given the owner's order (the architecture end to end, then the prior): (c) now,
  (a) when the intron's own solve is taken up.

## 8. LANDED (2026-09-04, after the owner's rulings)

**Stage A — the backbone and the protocol.** `messages/__init__.py` carries `Message` (composition
lane, level lane), `SILENCE`, `NO_NEIGHBOUR`, the `Prepared` protocol (`propagate(backward) → receive`,
`solve(from_left, from_right)`) and `Policy`; `sweep._pass` runs each pass in chain order and REFUSES a
kernel that leaves a real hop unspoken; the relay's `propagate` returns its scalar step as `receive`
and its `solve` rebuilds the source-indexed state it always read (its `NeighbourState` now lives in
`relay.py`) and publishes `fwd_*`/`bwd_*` itself; `silent.py` is four lines; `foundation.py`,
`policy.py`, their two gate files and the `message` policy name are gone. Gates: the four laws in
`test_sweep_backbone.py` (every node holds a message from each neighbour; chain order and one side per
pass; a mute kernel REFUSED; a policy that sends nothing leaves SILENCE; the solve receives the held
lists). Byte-identity: the suite's goldens unchanged; `policy_benchmark.py` standings of silent, relay
and transfer identical between the committed tree and the new one on three conditions.

**Stage B — the transfer policy as claims and rules.** `messages/transfer.py` rebuilt as §3: `prepare`
states every node's own claim and a rule per directed face (the ten messages as rules; the edge's count
as the level lane's claim), `propagate` composes what a node holds from its far side with its own claim
and applies the recipient's rule, `solve` adds the two held profiles. Bit-identical to `bp_proto.py`'s
formal arm on three conditions. Gates (`test_transfer_policy.py`): the backbone's rows against an
independent RECURSIVE reference of the passes; the no-echo perturbation (a node's claim replaced by a
spike never returns to it while it travels elsewhere); a claim-and-rule gate per family reusing the
existing independent recomputations; the constructor and predicate gates unchanged. The scan-seam gates
and the one-hop "beside nothing else" gates retired with the seam.

**Standings of the landed policy (`policy_benchmark.py --panel test`, silent / relay / transfer):** the
unstranded half — the zero controls 43,412 / 8,623 / 15,352 (OFF) and 40,315 / 6,372 / 26,594 (ON), the
in-scope `g50 ss.50 OFF` 9,349 / 12,813 / 12,328, `g98 ss.50 OFF` 6,620 / 8,817 / 6,006, and the deferred
capture-ON rows 0.05–0.18× silence; the stranded half — transfer at or below silence on 17 of 20 rows
(worst 1.10× at `g50 ss.70 OFF`), the relay above silence on every row. The ladder re-derivation (`landed_ladder.out`, silent / relay / transfer) reproduces the prototype's §6e
numbers exactly: the landed transfer beats the shipped relay on 13 of 16 rows and is at or below silence
on 14 of 16 (the two above are `g05 ss.50 OFF` 1.01× and `g05 ss.99 OFF` 1.00×); the relay leads only the
three `g00` rows its anchor pins (58,840 vs 233,420 at `g00 ss.50 OFF`; 152,534 vs 305,626 ON; 7,421 vs
17,851 at `g00 ss.99 OFF`) and ties `g00 ss.99 ON`. Where the relay is worst the landed policy is best:
`g98 ss.99 ON` 456,838 relay vs 224,402 transfer vs 298,597 silent.

## 9. THE TWO-SIDED EXON PROFILE PROBLEM, SHOWN ON THREE REAL EXONS (2026-09-04, for the owner)

Three `clean` exons of the test chromosome at `g50 ss.50 OFF` (unstranded, in scope), each with a
pure-gDNA intron beside it. The rows are log-likelihood in nats (0 = best supported, −30 = ruled out),
read at seven candidate gDNA shares of the exon; the truth is in the left column.

| exon (true gDNA share) | profile | 0.05 | 0.20 | 0.40 | 0.60 | 0.80 | 0.95 | 0.995 |
|---|---|---|---|---|---|---|---|---|
| 126 (**0.10**), 401 fragments | the intron's own (391 fragments) | −272 | −160 | −78 | −33 | −8.5 | −0.8 | 0 |
| | mapped into the exon through its face | −31 | −0.3 | 0 | 0 | 0 | 0 | 0 |
| 130 (**0.07**), 399 fragments | mapped into the exon | −10 | −0.1 | 0 | 0 | 0 | 0 | 0 |
| 186 (**0.04**), 1,051 fragments | mapped into the exon | −2.1 | 0 | 0 | 0 | 0 | 0 | 0 |

**What each row says.** The intron's own profile is a wall on the low side (a 10 % gDNA intron is ruled
out by 160 nats) and a gentle slope on the high side: 80 % costs 8.5 nats, 95 % costs 0.8, 100 % costs
nothing — density against the background cannot tell 95 % from 100 %. The face map carries that shape
into the exon and adds the face's counting: the result says "the exon has AT LEAST about 20 % gDNA" and
is flat from 40 % to 100 %. That is a one-sided profile. It never says "this much".

**Why the flat top is not merely uninformative but HARMFUL.** The position of the wall — the "at least"
— is set by the face's crossing count, which is 12–25 fragments here: Poisson noise of ±30–50 %. Exon
126's crossing drew 25 fragments where its intron's density predicts 12 (the crossing opportunity is
unbiased: pooled ratio 0.99–1.005 over 10,600 ladder faces on certified truth), so the wall landed at
about 18 % while the truth is 10 %, and the message charges the truth about 8 nats. When the crossing
draws LOW the wall lands below the truth and the profile is silent about everything above it — no
harm, no help. So the counting noise acts as a one-way ratchet: half the faces push their exon above
the truth and the other half say nothing. Under forwarding the ratchet compounds: an exon with no claim
of its own (every exon on unstranded data) passes the "at least" on through its next face, where it
becomes another "at least". That is the whole of the unstranded harm on the probe panels, and the
reason the prior — which trains on these exons — is fed a gDNA excess (`ISSUES:
gdna-landscape-trains-on-false-positives`).

**A two-sided profile** would put a peak at the value the evidence supports and charge both sides. It
cannot come from the crossing count alone (it is a count of the crossing, not of the exon); the flux
level was refuted (§6d); what remains is a sharper INTRON profile on its high side — the intron's own
solve — or an honest bounded marginal in the face map in place of the flat top.

**Is it the largest residual? No.** By node class on the ladder under the landed policy: off capture,
the INTRON class carries 43–45 % of the in-scope error (65k of 145k at `g50 ss.50 OFF`), identical under
every policy — the intron's own solve, not a message problem; exon|intron boundaries 12–14 %; terminus
boundaries 12–15 %; walled exons 8–10 %; alternative splice sites 8–9 %; the licensed exons this
section is about 6–7 %. On capture-ON stranded rows the terminus boundaries (28 %), alternative splice
sites (21 %) and walled exons (18 %) lead, licensed exons 10 %. The two-sided profile matters for two
other reasons: it is what a pass forwards INTO the walled and terminus classes, and at pass zero it is
the only thing that would give the landscape an unstranded exon population to train on (§6f).
