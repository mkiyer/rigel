# NEXT SESSION — THE SCAN AS FORMAL FORWARD-BACKWARD BELIEF PROPAGATION (owner ruling 2026-09-04)

    ⚠ A DEV DOC, and a HANDOFF. It says where things stand and how to start, not what is settled —
    rulings are `DESIGN.md`, the ranked list `ROADMAP.md`, the open problems `ISSUES.md`, the
    case-by-case state `MESSAGE_RUNGS.md`'s COMPLETION CHECKLIST. MOVE anything that settles.

## THE FRAME — the owner's ruling, which is the definition of done

⭐⭐⭐ **The sender does not decide; it sends. The RECIPIENT decides — forward, modify, or stop —
during the propagation phase. And the propagation is FORMAL FORWARD-BACKWARD: when it ends, EVERY node
holds TWO messages, one from each neighbour, except the two nodes at a chain's ends, which hold one.**
(`DESIGN.md` §6b.11; the completion contract it serves is §0c.0e.) ⛔ A hop that carries nothing must
still ARRIVE as an explicitly uninformative message — a node with no message from a side is a node the
policy never spoke to, and the solve cannot tell that from "nothing to say".

⛔ This thread does not switch away until every row of the checklist reads ✅. Other priorities are in
the tracker's PARKED-PRIORITIES LOG. ⛔ Every case starts with the SIMPLEST LOCAL form of the owner's
discrepancy rule; pooling and global models come after the tool works end to end
(`ISSUES: the-pooled-hop-step`).

## WHERE THE POLICY STANDS (branch `message-layer`, all committed)

Ten one-hop messages ship (`messages/transfer.py`, rows in `messages/transfer_rows.py`; every ruling in
`DESIGN.md` §6b.4–§6b.9), each priced by the owner's per-pair discrepancy rule where two witnesses
exist. **The scan seam is landed and INERT** (§6b.10): every delivery is recorded with its source, a
name and the adjacent map that carried it; the prepared object answers the backbone's `scan(backward)`
with its `(step, publish)` kernel and `deliver` fuses what each side forwarded; the budget `hops` is 0,
which builds no ledger and relays nothing, so the shipped answer is exactly the one-hop answer (gated
per slot against the previous policy on seven conditions, three gates, seven perturbations watched).

**Standings, re-derived on the committed tree** (`policy_benchmark.py --panel test|ladder --policies
silent relay transfer`; the seam is inert, so these are the ten one-hop messages' numbers):
ladder — beats silence 15/16, the relay 12/16; the four the relay still wins are the zero-gDNA controls
and `g98 ss.50 ON`, where its anchor pins the exons the transfer policy leaves to the landscape prior.
Test chromosome — 30 conditions, both halves apart. The published page:
https://claude.ai/code/artifact/f91ce211-1d5c-4d4f-842c-4f2d10161b3c

## WHAT IS NOT IMPLEMENTED — the whole list, in the owner's order

`MESSAGE_RUNGS.md`'s COMPLETION CHECKLIST is the authority (it carries slot counts, the IN/OUT of each
case and the measured state). In short, and in the order to take them:

1. ⭐⭐⭐ **THE SCAN, in the ruled form** (`multi-hop-scan`): a message per side at every node, the
   recipient's three decisions explicit, each hop priced. Its PREREQUISITE is (2). Its first cases are
   already measured un-premised: the forwarded terminus arrivals (−2 % at the in-scope zero control,
   +5 % at inside exons on `g50 ss.50 OFF`), the reverse direction (+0.8 %), and the chains of termini
   (`terminus-outside-empty`: 5,588 ladder boundaries, 3,792 of them reached only from another
   terminus). Subsumes `item 3`, the exon solve with every face speaking.
2. ⭐⭐⭐ **THE CERTIFIED-FLUX MESSAGE INTO EXONS** — the relay's anchor, ruled a message in §6b.3 and
   never delivered by this policy. Without it rung 2's transported row is a one-sided LOWER bound (flat
   above the face map's ceiling), two faces fuse to the higher of two noisy ceilings, and forwarding
   compounds that bias — which is why the scan's first measurement worsened the in-scope unstranded row.
   ⛔ A rough form (route rate × the exon's RNA opportunity) is REFUTED at the zero control: it claims
   19 % gDNA where there is none. Use `rna_anchor`'s own estimator (the sj opportunity with its overhang
   requirement, the route sum, the NB marginal). Two other forms are measured and recorded in the thread
   record: a two-sided Poisson row fixes capture-OFF and breaks capture-ON; the owner's abundance bound
   fixes OFF and still breaks ON where RNA dominates the total.
3. **THE LEVEL MESSAGE** — gDNA continuity where composition cannot cross: `strand-change-face`
   (1,138 faces), `termini-both-ways` (48), `intron-exon-face-terminus`'s inside exon (4,130), the empty
   chains. Two-sided, counting width, the abundance-discrepancy premise; the one-sided profile was
   refuted at exon|exon termini on sparse probes.
4. **`sj-plus-terminus`** (1,968 boundaries): the two landed maps composed.
5. **`ambig-tilt`** (4,671) and the **both-stranded locus** — rung 5, last by owner ruling; the tilt's
   ruling is whether it carries any message at all.
6. **Substrate still owed for two structures**: `nest` (a region walled by two termini) and `instart`
   (an alternative first exon inside an intron) — the isoform block grows ONE structure per step.
7. **`item 4`**: the factory on a region carrying both exon and intron bits — ☐ when first needed.
8. Then **the SHIP LIST**: the 0.8.0-metric pricing under `transfer`, the default flip, the obsolescence
   pass over relay/variance/rna_anchor/foundation, `preflight --full`, goldens.

## HOW TO START THE SCAN'S FORM

* **Read first**: `DESIGN.md` §6b.10 (the seam: the ledger, the two passes, the four gated laws) and
  §6b.11 (this ruling); then `messages/transfer.py`'s `_Ledger` and `_PreparedTransfer`.
* **The backbone's contract**: `sweep.solve_chain` calls `relay.scan(backward=False)` and then
  `(backward=True)`; each returns `(step, publish)`, `step(source, destination)` runs over the chain
  order (the forward pass reading each slot's LOW neighbour, the backward its HIGH), `publish()` hands
  back arrays the backbone gathers AT THE SOURCE, and `deliver(left, right)` receives them as two
  `NeighbourState`s. ⛔ `NeighbourState` is source-indexed by construction — a kernel may never read a
  destination's belief (`TRAPS: a-message-from-the-destinations-belief`).
* **The gates to extend**: `tests/calibration/test_transfer_policy.py`'s three scan gates — the inert
  budget, the kernel on a hand-built ledger, and the first real hop into an intron. The formal form
  needs one more: EVERY node holds a message from each side (the chain's ends one), which is a direct,
  falsifiable statement of the ruling.
* **The instruments**: `policy_prototype.py --by-class` (the node-local view: the error at each node
  class, which is how a message is judged at its destinations), `--panel test|test_junction|test_sparse`
  with `--all`, then `policy_benchmark.py --panel ladder`, halves apart. In the session scratchpad:
  `scan_proto.py` (the seam as a prototype plus three candidate exon rows), `walled_census.py`,
  `rung2_shape.py` (which showed the lower-bound shape), `scan_landing_identity.py` (per-slot identity
  against a committed policy — the pattern for judging any inert change).
* **The substrate**: `test_chr.yaml`, 133 genes at 720 k fragments per condition, four blocks — twin,
  mono, isoform (`altstart`, `altss`), and the WALLED block (`chain`, `tssalt`, `tandem`, `altlast`:
  72 walled pieces, the scan's stress test, designed from the ladder's walled-exon census).
  `TESTING.md` §0a has the rebuild recipe; ⚠ kill a rebuild's shard workers before relaunching.

## OWNER DECISIONS OUTSTANDING (recorded, not re-litigated)

Item 2's `g05 ss.99 ON` residue; item 5's +46 at `g05 OFF`; item 6's +1.2 % inside exons at `g50 ON`;
item 7's low-gDNA capture-ON rows on the benign and junction panels (`g05`, `g25`: +1.5…+9 % — an offset
of the licence below each pair's counting that only pooling would see, and pooling is refused) and
`g50 ss.99 ON` benign/sparse (+2.0 / +2.4 %). ⚠ The pre-walled panels are parked at
`~/Downloads/rigel_runs/test_reference/pre_walled/` (14 GB) and are deletable once the new panels are
trusted.

## THE LESSONS THIS THREAD PAID FOR (do not re-learn)

* **Judge a message at its DESTINATIONS** (`--by-class`, node-locally) beside the whole-library number,
  which carries the refit prior's response.
* **A width fit alone is blind to a consistent bias** under wide counting; and a message to a node whose
  belief the prior has already sharpened is disturbed by any residual bias.
* **Two premises the pairs cannot separate must not both be fitted** (a step and a flux factor at n ≈ 14
  are degenerate — `ISSUES: the-flux-factor-hop-premise`).
* **Log what an arm APPLIES before reporting it** — a leftover branch once made a "local" arm a hybrid.
* **Compare src-vs-src across a landing** (`TRAPS: a-harness-on-the-parent-class-dies-when-the-parent-gains-the-mechanism`):
  rebuild the previous policy from its own source text and compare per slot.
* **The simulator's capture tapers gDNA at a probed exon's edge** while a mature molecule's probe
  continues in transcript space — the mechanism behind item 7's step and rung 2's capture-ON behaviour.
