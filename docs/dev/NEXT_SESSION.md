# NEXT SESSION — FINISH THE MESSAGE POLICY: THE SCAN (multi-hop through the derived maps, each hop priced by the owner's rule), then sj+terminus (owner ruling 2026-09-02, `DESIGN.md` §0c.0e)

    ⚠ A DEV DOC, and it is a HANDOFF. It says where things stand and how to start, not what is
    settled — rulings are `DESIGN.md`, the ranked list is `ROADMAP.md`, the open problems are
    `ISSUES.md`, the message-rung order is `MESSAGE_RUNGS.md`. MOVE anything that settles.

## Where the thread stands (2026-09-02, branch `message-layer`, after items 2, 5, 6 and 7 landed)

* The `transfer` policy ships TEN messages (`messages/transfer.py`, rows in `messages/transfer_rows.py`):
  intron → boundary, boundary → exon, edge → exon, exon → boundary (item 1), boundary → intron
  (item 2), item 5's three at exon|exon TERMINUS boundaries (through the spliced-crossing map), item 6's
  abundance-discrepancy message into the inside flank, and item 7's at the ALTERNATIVE SPLICE SITE
  (both flanks, both directions, each carrying the fitted HOP PREMISE). Rulings: `DESIGN.md` §6b.2–§6b.8. Ladder, both bars, never pooled: unstranded
  ≤ the pre-item-5 policy on 6/8 (0.985× at the in-scope zero control), stranded ≤ on 7/8 and below
  silence on 8/8.
* ⭐ **Where the remaining error sits is now measured** (`policy_benchmark.py --panel ladder --policies
  silent transfer --by-class`, the tracker's SHIP LIST): exon|exon boundaries and walled exons hold
  62 % of the stranded capture-ON error and 85 % of the zero controls' — and item 5 showed why a
  message barely moves them: on the contaminated rows they sit at the resolution of their own
  evidence plus the prior (1.8 % of mass), and at the zero controls the error IS the prior's
  false-positive training. Introns off capture (42–46 %) are the vertex-atom thread. So the lever is
  now `ISSUES: gdna-landscape-trains-on-false-positives` (owner: "an important priority",
  2026-09-02), and rung 4's remaining structures are priced against a prior that already serves them.
* The owner's rulings of 2026-09-02, late: there is NO "one-hop paradigm" — multi-hop propagation with
  per-hop dampening is the direction (the composed transport is its first instance; the per-hop
  premise is fitted, never a constant); and the landscape prior's poisoning is a priority.
* The substrate: `scripts/sim/test_reference/test_chr.yaml` is the one hand-edited file; all seven
  panels certified 30/30; `TESTING.md` §0a is the recipe.
* Instruments left in the session scratchpad, worth promoting when next needed: a NODE-LOCAL scorer
  (|err| at a message's destinations beside every other class — the instrument that separated item
  2's and item 5's effects from the refit prior's response; `policy_prototype.py` is its natural
  home), a certified-truth pair-gap census by COMPONENT (the instrument that found the spliced
  crossing — test the MEAN gap in f-space by component; excess variance and a plug-in null cannot see
  a bias), and the terminus-boundary decomposition by what lies beyond the outside flank.

## START HERE: THE SCAN — multi-hop composition through every case already derived

⛔ The owner's ruling (2026-09-02, late): this thread does not switch away until every case in
`MESSAGE_RUNGS.md`'s COMPLETION CHECKLIST reads ✅. Other priorities are LOGGED in the tracker's
parked list, not taken up. ⛔ Start every case with the SIMPLEST LOCAL form of the owner's rule; pooling
and global models are for after the tool works end to end (`ISSUES: the-pooled-hop-step`).

**Why the scan is next.** Ten one-hop messages ship, each from a node's OWN evidence (a factory row or
a strand row). On unstranded data the strand rows are dead, so a walled exon — the largest error class
on the ladder (11,350 slots, 2.0 M mass; 35 % of the in-scope zero control's error, where the relay's
anchor still beats the transfer policy 5×) — receives nothing from the policy and is left to the
landscape prior. What can reach it is the composition its NEIGHBOUR received: the intron row that
arrived at the next exon through a licensed face, carried one hop further through the terminus or
alt-ss boundary between them. That is the scan: the policy's forward and backward step kernels on the
backbone's two directional passes, each hop through the map its boundary case already has, each hop
priced by the owner's discrepancy rule and nothing pooled.

**The backbone's protocol (`sweep.py` `_scan`, `messages/__init__.py`).** `Relay.scan(backward)`
returns `None` (relay nothing — what `transfer` does today) or a pair `(step, publish)`: the backbone
calls `step(s, i)` for every slot `i` in chain order with `s` its neighbour of the other kind
(`i−1` forward, `i+1` backward; a `−1` reference terminal is skipped, so nothing crosses a reference),
then `publish()` returns the state the pass produced. `deliver(left, right)` then receives two
`NeighbourState`s — each pass's published arrays gathered AT THE SOURCE slot with a `valid` mask — and
may fuse them into the slot's `PsiMessage.lam_rows`. ⛔ `NeighbourState` is indexed at the source by
construction (`TRAPS: a-message-from-the-destinations-belief`): a kernel carries what the SOURCE holds,
never a destination belief. The relay's own kernel (`relay.py` `scan`) is the worked example of the
shape; its content (a scalar level with `_damp`) is what the rebuild replaced.

**The design to derive, one piece at a time.** (1) The state a pass carries per slot: a composition
ROW in λ (the currency of every message) plus what is needed to price the next hop — the row's own
width is inside its shape; the per-pair discrepancy needs the two witnesses at the hop. (2) The step
kernel per boundary case, each reusing the landed map: intron|exon licensed face → `face_map_lambda`
+ `transport_row` (rung 2's form); exon|exon terminus → §6b.6's spliced-crossing map (`splice_out_row`
into the boundary, the face map with the spliced density out of it), the inside flank through §6b.7's
abundance map; alt-ss → §6b.8's two maps; a strand change or a refused face → no composition crosses
(the level message, phase A, is a separate case). (3) The pricing at each hop where the destination
has no strand witness (unstranded data, the purpose): the owner's ABUNDANCE-DISCREPANCY rule — the
ratio of the two objects' total abundances is measured on every library, and `abundance_row` already
turns it into width with the two hypotheses; where both strand witnesses exist, the per-pair rule of
item 7. Nothing pooled. (4) No echo: a forwarded row never returns to the slot it came from (item 5's
snapshot rule made this structural for one hop; the scan needs it for every hop — the backbone's
direction split gives it for free if the forward state is built only from forward arrivals).

**The first cases, already measured un-premised** (`DESIGN.md` §6b.6, §6b.9): forwarding a terminus
boundary's arrivals one hop further (−2 % at the in-scope `g00 ss.50 OFF` control, +5 % at inside exons
on `g50 ss.50 OFF` — over-claiming without a premise), the reverse direction into the boundary
(+0.8 %), and the chains of termini (`DESIGN.md` §6b.6's census: a hop of ≤ 20 bases carries
composition within counting; half the terminus-boundary error is empty outside pieces whose far face
is another terminus). Stage 0 on certified truth first: the composition transported k hops through the
landed maps against the truth at the destination, by hop count and boundary case, with the counting
width beside it — the instrument that says how fast composition degrades along a chain and which
premise each case needs.

**Acceptance, unchanged:** improves or stable with minimal harm on the test chromosome (all three probe
panels, `policy_prototype.py --all --by-class`), the ladder (halves apart) — node-locally at the
destinations beside whole-library; every gate fail-first, every perturbation watched; faithfulness of
the landed code against the prototype's RECORDED numbers on a parent WITHOUT the mechanism.

**After the scan, in the tracker's order:** sj+terminus (item 10, the two maps composed), the level
message where composition cannot cross (phase A: strand-change faces, termini both ways, the empty
chains), the AMBIG tilt and the both-stranded locus, then the ship protocol.

⭐ OWNER DECISIONS OUTSTANDING (each landed with the residue recorded, none re-litigated): item 2's
`g05 ss.99 ON` residue; item 5's +46 at `g05 OFF`; item 6's +1.2 % inside exons at `g50 ON`; item 7's
low-gDNA capture-ON rows on the benign and junction panels (`g05`, `g25`: +1.5…+9 % — an offset of the
licence below each pair's counting that only pooling would see, and pooling is refused) and
`g50 ss.99 ON` on the benign and sparse panels (+2.0 / +2.4 %).

## THE LESSONS FROM THE RUNG-4 EXCURSION (2026-09-02), carried so they are not re-learned

The next work is rung 4 — the exon|exon boundary types, starting with the TERMINUS
case on `altstart` (tracker items 5–9), before item 3 (the exon solve with every face). This session
began that rung and pulled back; nothing it learned is to be re-derived. The record is the thread
record's RUNG 4 sections (the census, the certified directed licence, the derivation of the three
mechanisms, the prototype tables, the ladder) and, for the rules that became rulings, `DESIGN.md` §6b.4.
The lessons, so they are carried and not lost:

1. **The census re-framed the hole.** Of the ladder's 12,811 exon|exon boundaries, 7,604 are internal
   TERMINI (a TSS/TES of one isoform inside another's exon; 950 k crossing mass) and 4,838 alternative
   splice sites; 10,259 exons (1.05 M) are unreachable by any chain of licensed faces, walled by
   termini. The terminus is the mechanism to build; the splice-site class is licensed already.
2. **The directed licence is CERTIFIED on ladder truth.** At a terminus boundary the flank the
   terminating transcripts do NOT cover (the outside) shares the crossing's composition to within
   counting noise on every condition, capture on or off; the inside flank does not. Direction from the
   flag alone: TSS+ / TES− bodies extend right (outside = left flank); TES+ / TSS− extend left; mixed
   flags → no side. Verified on both strands of the built `altstart` structures.
3. **Composition into the boundary from the outside flank is the existing currency** (the s = 0 face
   map, an opportunity shift). What was NOT settled: where the outside EXON's composition comes from
   on unstranded data (a composed transport of one intron row through a chain of licensed faces was
   prototyped: join steps only, leave steps refused; NEUTRAL on the ladder, −6 % at the g00 zero
   control; splice-out steps would reach 1,741 more exons). Now that item 1 exists, an outside exon
   on STRANDED data has its own row to send — the same message as item 1 without the subtraction.
4. **The inside flank: a gDNA lower bound was derived (the edge rung's profile generalised — the source
   row's running maximum through the level map) and REFUTED for shipping**: inert where a witness or
   the refit prior exists, +8 % harm on sparse probes because the enrichment sign is not certified at
   an exon|exon terminus. Do not rebuild it without a new measurement; what may cross into the inside
   is an open question, not a bound.
5. **The refit prior already reaches walled exons** (a `nest` wall read 0.666 vs truth 0.630 under the
   shipped policy where silence reads 0.001): the 1 M walled mass did not turn into whole-library
   error, so a terminus mechanism is judged on the BOUNDARIES it solves and the in-scope rows it must
   not harm, not on a large expected win.
6. **Two rules of arithmetic proven this session apply to every message to come**: both components
   convert counts to densities on ONE, capture-blind opportunity (a capture-aware opportunity on one
   component alone re-introduces a level across locales), and a width is a marginal over the measured
   ingredient, not a uniform blur in log-odds (noise in a ratio moves the log-odds by `σ/(1−f)`).
7. **Substrate**: one structure per step. `altstart` is present (an exon|exon terminus, both strands);
   `altss`, `nest`, `instart` are one YAML edit each (their geometry is in the thread record's RUNG 4
   substrate section and the YAML header), rebuilt by `panel.py` in ~40 minutes for all seven panels.

## Standing cautions

* One message per item; a moved number must have ONE cause; adversarial probe panels in every loop.
* A claim below its own evidence is silence, never a near-zero row.
* Accepted errors, do not re-litigate: the zero-gDNA edge residual; the deferred `g05 ss.50 ON` row;
  the splice-out premise bias (`ISSUES: splice-out-premise-bias-uncorrected`, the owner's call).
* `ISSUES: gdna-landscape-trains-on-false-positives` is the exposed systemic issue; the owner ranks it.
* The owner's `docs/dev/rename.md` (`row`, `drain`) is a rename campaign of its own — `rename_census.py
  --sense` first, never a tail-end sweep.
