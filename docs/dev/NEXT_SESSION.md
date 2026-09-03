# NEXT SESSION — FINISH THE MESSAGE POLICY: the scan (multi-hop with per-hop premises), then sj+terminus (owner ruling 2026-09-02, `DESIGN.md` §0c.0e)

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

## START HERE: the completion checklist — the scan, then sj+terminus

⛔ The owner's ruling (2026-09-02, late): this thread does not switch away until every case in
`MESSAGE_RUNGS.md`'s COMPLETION CHECKLIST reads ✅ — no nullified message, no skipped boundary, multi-hop
through chains with per-hop dampening, gDNA always conveyed where composition cannot cross, every node
solved from two honest messages. Other priorities are LOGGED in the tracker's parked list, not taken up.

Items 5, 6 and 7 landed (`DESIGN.md` §6b.6–§6b.8). Item 7 is the template for every hop the scan will
take: a licence certified on truth, then THE OWNER'S DISCREPANCY RULE PER PAIR — where a pair's two
witnesses (the boundary's own strand mode and the flank's mapped to it) disagree beyond counting, that
pair's messages are widened by the excess; nothing pooled, no mode shifted. ⛔ The lesson the owner
paid for on 2026-09-03: a per-library POOLED step was built first and refused as over-engineering —
"how do we know the behaviour of other node pairs is predictive of a global pattern?" Start with the
simplest elegant form, finish the tool end to end, and only then take up subtler accuracy work (a
global model of disagreement, such as the production relay's projection of a pair onto the
total-abundance landscape, is that kind of later work). Two more lessons: a second-moment fit alone is
blind to a consistent bias under wide counting, and two premises the pairs cannot separate must not both
be fitted.

⭐ OWNER DECISIONS OUTSTANDING (each landed with the residue recorded, none re-litigated): item 2's
`g05 ss.99 ON` residue; item 5's +46 at `g05 OFF`; item 6's +1.2 % inside exons at `g50 ON`; item 7's
low-gDNA capture-ON rows on the benign and junction panels (`g05`, `g25`: +1.5…+9 % against the
pre-item-7 policy — a systematic offset of the licence that hides below each pair's counting; a shift
would recover it and shifts are refused, `ISSUES: the-pooled-hop-step`) and `g50 ss.99 ON` on the benign
and sparse panels (+2.0 / +2.4 %); the ladder itself reads four wins and two harms under 0.12 %; and
the commit of the local form.

Next: THE SCAN — the policy's forward/backward step kernels on the backbone's two scans, one per
boundary case, each hop priced by the item-7 rule (the pair's own disagreement beyond counting); its first
cases are the two HELD pieces (forwarding a boundary's arrivals one hop further, and the reverse
direction into the boundary — both measured, both over-claiming without a premise) and the chains of
termini (`DESIGN.md` §6b.6: a hop of ≤ 20 bases carries composition within counting). What to fit on an
UNSTRANDED library, where no strand modes exist, is the open derivation. Then sj+terminus (the two maps
composed), the rule at strand-change faces and termini both ways, the tilt ruling and the both-stranded
locus, then the ship protocol.

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
