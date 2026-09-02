# NEXT SESSION — item 2 of the message rungs: the intron|exon BOUNDARY → INTRON message

    ⚠ A DEV DOC, and it is a HANDOFF. It says where things stand and how to start, not what is
    settled — rulings are `DESIGN.md`, the ranked list is `ROADMAP.md`, the open problems are
    `ISSUES.md`, the message-rung order is `MESSAGE_RUNGS.md`. MOVE anything that settles.

## Where the thread stands (2026-09-02, branch `message-layer`, HEAD after the cleanup commit)

* The `transfer` policy ships four messages (`messages/transfer.py`, rows in `messages/transfer_rows.py`):
  intron → boundary, boundary → exon, edge → exon, and — item 1, landed 2026-09-02 — exon → boundary.
  Rulings: `DESIGN.md` §6b.2–§6b.4. Ladder, both bars, never pooled: unstranded 7/8 wins, stranded
  every capture-ON row won by item 1 (0.987–0.995×), capture-OFF within 33 fragments.
* The tracker `MESSAGE_RUNGS.md` carries the order. Rung 1 is COMPLETE; rung 2 is partial (items 2
  and 3); rung 4 holds one structure (`altstart`) and resumes after items 2–3; rung 5 (both-stranded)
  is last.
* The substrate: `scripts/sim/test_reference/test_chr.yaml` is the one hand-edited file (55 genes, both
  strands balanced; twin + mono + `altstart`); all seven panels certified 30/30; the recipe is
  `TESTING.md` §0a and `panel.py cache` now does the g00 pre-warm, `_main` copy and certification itself.
* Baselines to measure from: suite `CLAUDE.md`; `policy_benchmark.py --panel test --policies silent
  transfer` (unstranded 14/20 worst 1.16× at the deferred zero control; stranded 10/10 after item 1).

## START HERE: item 2 — boundary → intron

The question, as the owner framed it: an intron has its own composition measurement (the density
factory against the intergenic background); can its two boundaries add anything, and what? Both share
the intron's population exactly (mature RNA crosses neither), so the transfer is the opportunity shift
only — no splice-out, no premise. What a boundary holds that the intron does not is its OWN strand
evidence at the face and, under capture, depth (a probed exon's faces carry hundreds of crossings; the
unprobed intron ten or twenty). So the honest outcomes are "the intron takes the boundaries' strand
rows, opportunity-shifted, where the deadband declares them live" or "nothing" — both measured, not
assumed. Derive first, teach every piece, check on certified truth (the pair gap of stage 0 already says
the composition is shared), prototype through `scripts/design/policy_prototype.py` (subclass the shipped
policy — ⚠ `TRAPS: a-harness-on-the-parent-class-dies-when-the-parent-gains-the-mechanism`), falsify
by reversing the boundary's claim, ladder both halves apart, then `src/` with fail-first gates.

## THE ORDER AFTER ITEM 2 (owner, end of 2026-09-02): the exon|exon boundaries, terminus structures first

After item 2 ships, the next work is rung 4 — the exon|exon boundary types, starting with the TERMINUS
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
