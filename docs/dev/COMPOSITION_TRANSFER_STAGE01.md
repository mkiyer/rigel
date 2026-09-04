# THE COMPOSITION-TRANSFER THREAD — the working record of the `transfer` policy (2026-09-01 → )

    ⚠ A DEV DOC — working results of the owner-approved plan in `message_notes.md`. Nothing here
    is settled; when a verdict lands it MOVES to its permanent home and this file shrinks.

**The mechanism under test (owner design, approved 2026-09-01):** one node type at a time,
starting with intron|exon BOUNDARIES; solve them by COMPOSITION TRANSFER from the adjacent
intron REGION (the shared-population pair: mature RNA cannot cross an intron|exon boundary, so
both objects' unspliced population set is {gDNA, unspliced RNA}); the exon-side message is OFF;
the open question is the transfer's HONEST PRECISION. Strictly one mechanism, zero `src/`
changes — the prototype runs through the shipped `message_policy = "message"` seam by patching
`rigel.calibration.calibrate.MessagePolicy` in-process, so every phase and refit runs it.

Substrate: the test chromosome, all 30 conditions, **100** intron|exon pairs. ⛔ Node-local
numbers below are NOT the 0.8.0 metric and NOT the panel — whole-library scoring and the ladder
confirmation are still owed before any shipping claim (`TRAPS: a-toy-and-a-panel-can-disagree-in-rank`).

## Rungs 1–3 and the rung-4 excursion — SHIPPED / SET ASIDE; the record moved (2026-09-03)

Stages 0–5 of rung 1 (the certified pair gap, the epsilon ladder, the two precision laws, the blur
constant deleted, the promotion), rung 2's derivation and design (the route reformulation, the
anchor's cross-locale assumption refused, the face-composed transfer, the ingredient width, the three
probe worlds, the ladder confirmation), rung 3 (the sign-certified lower bound, the mono block, the
ceiling refused) and the rung-4 excursion (three mechanisms prototyped together, the composed
transport neutral, the inside bound refuted) now live in `DESIGN.md` §6b.9, with the systemic issue in
`ISSUES.md` (`gdna-landscape-trains-on-false-positives`). Git holds the full working text
(commit `08b03c46` and before).

# THE RESET (owner, 2026-09-02, late)

Rungs 1 and 2 were unfinished (the exon-side messages had been nullified) and come first, one message
at a time; rung 4 grows one structure per step (`altstart` only); `MESSAGE_RUNGS.md` is the tracker.

## THE BASELINE ON THE TRIMMED SUBSTRATE (55 genes: twin + mono + `altstart`; 480 k fragments; 2026-09-02)

`policy_benchmark.py --panel test --policies silent transfer`, the shipped rungs 1–3 against silence —
the number every item in `MESSAGE_RUNGS.md` is measured from. Unstranded: transfer beats silence on
**14/20**, worst **1.16×** (`g00 ss.50 ON`, the deferred zero control, 64,143 vs 55,290); stranded
**8/10**, worst **1.00×**. In-scope unstranded OFF rows: `g05` 1.00×, `g25` **1.06×** (12,738 vs
12,042 — the one in-scope row where the shipped policy costs), `g50` 0.99×, `g98` 0.92×. Blind
capture-ON rows: `g05` 0.08×, `g25` 0.76×, `g50` 0.05×, `g98` 0.07×. Full table in the session's
`baseline_trimmed_silent_transfer.txt`; every earlier test-chromosome number is on another substrate.

## Item 1 — SHIPPED (2026-09-02); the record moved

The exon → intron|exon boundary message: derivation, the certified-truth check, the prototype's
tables, the transfer variance (counting + premise, the two-witness estimator, the marginal width, the
geometric-opportunity form), the three-panel and ladder tables and the landing all MOVED to their
homes — the ruling and measurements to `DESIGN.md` §6b.4, the premise-bias decision to
`ISSUES: splice-out-premise-bias-uncorrected`, the harness double-count lesson to
`TRAPS: a-harness-on-the-parent-class-dies-when-the-parent-gains-the-mechanism`, the gates to
`tests/calibration/test_transfer_policy.py`. The session harness (`item1_proto.py`, scratchpad) was
promoted as `scripts/design/policy_prototype.py`.

## Item 2 — SHIPPED (2026-09-02); the record moved

The intron|exon boundary → intron message: the derivation (one shared population, the s = 0 map the
identity, the boundary's own strand row verbatim, the derived deadband at the boundary), stage 0 on
certified truth (zero excess variance over counting on all four substrates; the plug-in Poisson null
that exposed the digamma estimator's vertex artifact; the ladder's mate-gap mature crossings at
0.3–0.7 % of boundaries), the information budget, the node-local and whole-library tables on the three
test panels and the ladder, the reversed-row falsification and the prior-mediated `g05 ss.99 ON`
residue all MOVED to `DESIGN.md` §6b.5; the prior's second exposure to
`ISSUES: gdna-landscape-trains-on-false-positives`; the gates to `tests/calibration/test_transfer_policy.py`.
Two session instruments (the node-local scorer, the pair-gap census with its null) are described in
`NEXT_SESSION.md` for promotion. The prototype (`item2_proto.py`, scratchpad) subclassed the shipped
policy with `strand` passed THROUGH so item 1 stayed on in both arms — one thing varied — and is
retired with the landing.

# ITEM 5 — SHIPPED (2026-09-02); the record moved

The exon|exon terminus boundary from the outside flank: the orientation table, the certified
directed licence on both substrates, the census of the ladder's terminus boundaries by what lies
beyond the outside flank (half are chains of termini a dozen bases apart), the REFUTATION of the
verbatim own-row exchange and its diagnosis by component (the mature depletion of the unspliced
crossing), the corrected licence with the spliced crossing, the three messages, the node-local and
whole-library tables on all four substrates, the structural no-echo law and the two stage-0 lessons
all MOVED to `DESIGN.md` §6b.6; the gates to `tests/calibration/test_transfer_policy.py`; the chain-
of-termini derivation owed to the tracker's item 6b. The prototypes (`item5_proto.py` and its
refuted verbatim form `item5_proto_v0_verbatim.py`, scratchpad) are retired with the landing.

# ITEM 6 — SHIPPED (2026-09-02); the record moved

The owner's abundance-discrepancy rule into the inside flank of a terminus: the certified bracket
check, the three marginal forms (uniform, profile, fitted step) and their agreement, the falling-totals
sign error every prototype shared and its corrected support, the isolation of the ingredients (the
forwarded arrivals — held for the scan with a per-hop premise; the reverse direction — held), the
self-fitted hop premise, and the node-local, three-panel and ladder tables all MOVED to `DESIGN.md`
§6b.7; the gates to `tests/calibration/test_transfer_policy.py`. The prototype (`item6_proto.py`,
scratchpad, eleven arms) is retired with the landing.

## Item 7 — SHIPPED (2026-09-02); the record moved

The alternative splice site: the two licences (C through §6b.6's map with `S_b`, E through §6b.4's map
with `S_b + F`), the stage-0 certification on both substrates and by flank length, the component
dissection that named the capture taper, THE HOP PREMISE (the fitted step with its error, the excess,
the per-pair discrepancy), the refused flux-factor forms and every measurement now live in
`DESIGN.md` §6b.8 and `ISSUES.md` CLOSED/REFUSED (`the-flux-factor-hop-premise`); the substrate's
`altss` block and the REPLICATION RULE in `TESTING.md` §0a's YAML header. Instruments left in the
session scratchpad: `item7_stage0.py` / `item7_stage0_len.py` (the licences on certified truth, by type
and by flank length), `item7_components.py` (the by-component dissection), `item7_nodelocal.py` (the
node-local scorer at alt-ss objects), `item7_proto.py` (every arm tried: own rows, width-only fit,
step forms, flux forms, joint fit, per-pair discrepancy).

## The pooled step questioned and REFUSED (owner, 2026-09-03); the record moved

The owner's critique of item 7's per-library step, the purely local form measured against it
(on the ladder within 0.25 % of the pooled form on every row and ahead on four of six stranded rows;
behind on the test chromosome's low-gDNA capture-ON rows), the refusal and its grounds, and the lesson (simplest elegant form first; end-to-end before subtler accuracy work) live in
`DESIGN.md` §6b.8 and `ISSUES.md` CLOSED/REFUSED (`the-pooled-hop-step`). The local rule is what ships.

## ⛔ A MEASUREMENT ERROR OF MINE, corrected (2026-09-03)

The "purely local" prototype arm whose numbers I reported to the owner (node-locally 177 at the g50 ON
boundaries; the ladder within +0.45 % of the pooled form) was NOT purely local: a leftover branch in
the prototype re-set the E hop's width and shift to the POOLED values after my local block, so it was a
hybrid — a local C hop and a pooled E hop. Found by logging the widths each implementation applies
(the landed local form: 96 non-zero widths summing to 57.8; the prototype's "local" arm: 184 summing
to 10.8, the recurring 0.0509 being the pooled E width). The lesson: before reporting an arm, log
what it APPLIES, not what its flags say. The landed code is truly local on both hops and the corrected
prototype arm reproduces it exactly (3,933 / 229 / 121 / 310 on benign g50 ss.99 ON).

**The true local form against the pre-item-7 policy (full 30-condition sweeps, whole-library):**
benign — capture-OFF stranded rows within ±0.2 %, unstranded identical on 9/10; capture-ON stranded
g05 +1.7 / +2.1 %, g25 −0.8 / +1.5 %, g50 −1.0 / +2.0 %, g98 −4.1 / −3.2 %. Junction panel — ON: g05
+9.0 / +7.7 %, g25 +6.1 / +3.1 %, g50 0.0 / −0.6 %, g98 +0.5 / −2.2 %. Sparse panel — ON: every row a
win (−1.6 … −10.4 %) except g50 ss.99 (+2.4 %). Beside the pooled form: better where a pair's
disagreement is large enough to see (junction g50 ss.99 ON 3,029 vs 3,149), worse where a systematic
offset hides below each pair's counting (junction g05 / g25; benign g50 ss.99 ON 3,933 vs 3,873; sparse
g50 ss.99 ON 4,118 vs 4,011) — exactly the case the per-pair rule gives up by design. Node-locally at
the 30 alt-ss boundaries (g50 ss.99 ON): 126 pre-item-7, 249 un-premised, 229 local, 171 pooled.
The corrected ladder A/B (`proto7L_ladder_*`: pre7 / local / pooled) is the standing.

**Corrected ladder A/B (pre-item-7 → local / pooled):** g05 OFF +0.12 / +0.17 %, g05 ON −1.34 / −1.38 %,
g50 OFF +0.06 / +0.08 %, g50 ON −0.44 / −0.67 %, g98 OFF −0.76 / −0.75 %, g98 ON −2.70 / −2.49 %;
unstranded identical. The local rule ships.

# THE SCAN — the machinery built, the zero point held, and what it exposed (2026-09-03)

**The machinery — LANDED INERT in `src/` on 2026-09-03 (`DESIGN.md` §6b.10 carries the laws and the
zero point); what follows is the prototype it was measured with (`scan_proto.py`, session scratchpad,
which also holds the three candidate exon rows below).** The landed policy's `prepare` restructured
so every one-hop delivery is recorded with its SOURCE, its KIND and the adjacent MAP that made it
(`arrivals`, `maps`, `consumed`); `scan(backward)` returns the backbone's `(step, publish)`: the forward
pass carries each slot's LEFT-side arrivals (those today's policy has not already carried across that
hop — rung 1's row across the boundary into the exon, and the outside exon's arrivals into the terminus
boundary, are marked consumed) plus the state it received, through the adjacent map, into the next
slot; the backward pass mirrors it; the published state is what a slot SENDS (the backbone gathers at
the source); `deliver` adds the two forwarded terms to today's rows. Maps are registered by GEOMETRY,
not inside the deadband guards, so an unstranded library's chains do not break at the first exon.
`HOPS = 0` is BYTE-IDENTICAL to the landed policy (checked exactly on four conditions, the unstranded
in-scope row and the zero control among them) — the owner's one-hop fallback, structural.

**Scan v0 (today's rows forwarded, map widths only), whole-library, test chromosome:**

| row | transfer | 1 hop | 2 hops | 6 hops |
|---|---|---|---|---|
| g50 ss.50 OFF (in scope, unstranded) | 11,242 | 11,115 | 11,909 | 12,727 |
| g98 ss.50 OFF | 6,363 | 6,361 | 6,117 | 6,149 |
| g50 ss.99 ON | 6,367 | 6,275 | 6,102 | 6,123 |
| g05 ss.99 ON | 2,293 | 2,271 | 2,282 | 2,171 |
| g98 ss.99 ON | 6,854 | 6,482 | 6,083 | 5,196 |
| g00 ss.50 OFF (zero control) | 16,819 | 16,819 | 16,490 | 15,447 |

The walled exons: g98 ss.99 ON 135 → 115 at six hops, the zero control 60 → 44; on the in-scope
unstranded row 125 → 225 — worse with every hop, and the licensed exons 4,071 → 5,316. The forwarded
rows compound a bias that was already in them.

**What the bias is — rung 2's row is a LOWER BOUND.** The first reading of the rebuilt chromosome (the
in-scope unstranded row 1.04× → 1.20× against silence, at licensed exons of every type) was traced to
the substrate (the parked panels read the same slots near truth with the current code), then to the
row's SHAPE (`rung2_shape.py`, real factory rows): `transport_row` reads the intron row at the map's
preimage and takes the flat limit above the map's ceiling, so the transported row penalises exon shares
BELOW the ceiling and is flat above it — "at least this much gDNA", never a peak. Two faces fuse to the
HIGHER of two noisy ceilings (a clean exon of true share 0.324: ceilings 0.41 and 0.26, plateau from
0.5 up, the exon lands at 0.52); the ceiling's noise is the crossing count's (9–25 fragments), so the
bias is largest at low depth, and the walled block added many low-depth licensed exons.

**Four forms of the row, node-locally on six rows (whole-library; licensed exons in brackets):**

| row | transfer | two-sided Poisson | abundance-bounded | flux cap (route rate) |
|---|---|---|---|---|
| g50 ss.50 OFF | 11,242 [4,071] | **8,945** [1,874] | **8,936** [1,868] | 11,027 |
| g98 ss.50 OFF | 6,363 [2,529] | 5,604 [1,711] | 5,556 [1,666] | 5,712 |
| g50 ss.99 ON | 6,367 [4,216] | 15,558 [13,086] | 12,333 [9,970] | 10,548 |
| g05 ss.99 ON | 2,293 [1,426] | 3,188 [2,247] | 3,090 [2,160] | 3,297 |
| g98 ss.99 ON | 6,854 [4,519] | 7,483 [5,113] | 5,660 [3,331] | 4,573 |
| g00 ss.50 OFF | 16,819 [6,615] | 15,506 [5,757] | 15,508 [5,758] | 121,263 [96,147] |

(1) The TWO-SIDED Poisson form (the likelihood of the crossing count under the exon's hypothesised
share, marginalised over the intron row's composition) fixes capture-OFF — the in-scope unstranded row
goes BELOW silence — and is catastrophic capture-ON: under capture the exon interior's gDNA exceeds
what its tapered edge crossing implies, and the flat top was tolerating exactly that. (2) The
ABUNDANCE-BOUNDED form (the owner's discrepancy rule on this face: the crossing count fixes the level up
to an enrichment step bounded by the exon's measured total against what the crossing and flux predict;
`transport_row_bounded`) keeps the OFF wins, lets the scan's hops help the walled exons on unstranded
rows for the first time (125 → 113, 115 → 92, the zero control 60 → 36), and still fails at g50/g05
ON: a bound on the TOTAL abundance is blind to the gDNA's enrichment where RNA dominates the total,
which is most probed exons. (3) What measures the exon's gDNA LEVEL two-sidedly on and off capture is
the exon's own count minus the RNA its certified flux implies — THE RELAY'S ANCHOR, which `DESIGN.md`
§6b.3 ruled a message and which the transfer policy never received (it stayed in the relay behind
`certified_flux`; the relay's lead at every zero control is that anchor). (4) The quick flux-cap arm
(`route_rate × the exon's RNA opportunity` as the exon's RNA) is REFUTED at the zero control: it claims
19 % gDNA at exons that have none — the route rate under-states the contained RNA systematically (the
sj opportunity's overhang requirement is the likely cause), which is precisely what the anchor's
estimator (`rna_anchor`: the sj opportunity, the route sum, the NB marginal) exists to get right.

**The verdict and the next case.** The scan machinery is sound and its zero point holds; forwarding
helps stranded capture-ON rows and the zero control already, and cannot help the unstranded rows it
exists for until the exon's rows are two-sided. That needs the certified-flux message into exons —
the anchor, face-local, as a transfer message, with the anchor's own estimator of the exon's RNA — which
also turns rung 2's lower bounds into peaks and explains the walled block's first reading. Order: land
the scan machinery at HOPS = 0 (no behaviour change; the seam and its bookkeeping in `src/`, gated
byte-identical), then the certified-flux message as its own case (derive from `rna_anchor`'s estimator;
stage 0 on truth: the flux-implied RNA against the exon's true RNA, by depth and probe state), then the
hops with the rows two-sided, priced per hop by the discrepancy rule.
