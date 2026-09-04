# NEXT SESSION — THE TWO-PHASE BACKBONE IS LANDED; NEXT IS THE TWO-SIDED EXON PROFILE (handoff, 2026-09-04)

    ⚠ A DEV DOC, and a HANDOFF. It says where things stand and how to start, not what is settled —
    rulings are `DESIGN.md`, the ranked list `ROADMAP.md`, the open problems `ISSUES.md`, the
    case-by-case state `MESSAGE_RUNGS.md`'s COMPLETION CHECKLIST. MOVE anything that settles.

## WHAT HAPPENED ON 2026-09-04

The owner asked for the message layer's ARCHITECTURE to be solidified — two phases, `propagate` and
`solve`, the naming debt paid — and ruled the same day: retire the foundation scaffold, re-found the
skeleton as designed, judge every idea at pass zero and with the prior apart, finish the architecture
end to end before the landscape prior. **Both stages LANDED** (`DESIGN.md` §6b.12; the record is
`TWO_PHASE_BACKBONE.md` §8): the two-phase backbone with the relay byte-identical and the foundation
scaffold gone; and `messages/transfer.py` rebuilt as claims + rules + the passes, bit-identical to the
prototype's formal arm, which wins both halves of the ladder against the one-hop policy (7/8 + 7/8, at
pass zero and through the pipeline). The prerequisite the previous handoff named — the certified-flux
row into exons — was derived, stage-0'd, prototyped and REFUTED as a level by probe placement
(`ISSUES: the-certified-flux-row-as-a-level`); what remains of it is `ISSUES: two-sided-exon-row`.

## WHERE THE POLICY STANDS (uncommitted in the working tree — the owner drives commits)

`policy_benchmark.py --panel test --policies silent relay transfer` (30 conditions): the stranded half
at or below silence on 17/20 (worst 1.10×), the relay above silence on every stranded row; the
unstranded half — the deferred capture-ON rows 0.05–0.18× silence, the in-scope `g50 ss.50 OFF`
12,328 vs silent 9,349 (the one-sided exon profile compounding, the open problem), the zero controls
0.35× / 0.66× silence where the relay's anchor still leads (8,623 / 6,372). The ladder table is
`landed_ladder.out` (session scratchpad; §6e of the note by construction).

## THE NEXT CASES, in the owner's order

1. **`ISSUES: two-sided-exon-row`** — what makes an exon's profile two-sided on unstranded data, now
   that forwarding is structural: (a) the intron's own solve (`ROADMAP.md` rank 3), or (b) a bounded
   taper marginal in the face map. ⛔ Not the flux level. Judge at `R exon (licensed)` and the walled
   classes, pass zero beside the full pipeline (`pass0_score.py`), all three probe panels, the ladder.
2. `sj-plus-terminus` (1,968 ladder boundaries): the two maps composed — a RULE at that face.
3. THE LEVEL MESSAGE where composition cannot cross (strand-change faces, termini both ways, the empty
   chains): the `Message.level` lane exists and carries one rule (the edge bound); the rest are owed.
4. `ambig-tilt` and the both-stranded locus (rung 5, last by ruling).
5. The SHIP LIST: the 0.8.0-metric pricing under `transfer`, the default flip, the obsolescence pass
   (`relay.py`, `variance.py`, `rna_anchor.py`), `preflight --full`, goldens.

## HOW TO START

* Read `TWO_PHASE_BACKBONE.md` §0 (the words), §3 (the skeleton as landed), §6d–§6f, §8.
* The prototypes and logs are in the 2026-09-04 session scratchpad
  (`/private/tmp/claude-503/-Users-mkiyer-proj-rigel/d7adfcc1-…/scratchpad/`): `bp_proto.py`,
  `flux_proto.py`, `flux_stage0.py`, `bp_identity.py`, `pass0_score.py`, `halves.py` and every `.out`.
  ⚠ Scratchpads are per session: `flux_stage0.py`, `bp_identity.py` and `pass0_score.py` are worth
  promoting into `scripts/design/` (each adds +4 collected cases and needs a docstring the gates accept).
* Every A/B: `policy_prototype.py --panel test|test_junction|test_sparse --all --by-class`, the halves
  apart, pass zero beside the full pipeline, then the ladder. ⛔ Never quote one panel.
* The suite baseline after the landing is in `CLAUDE.md`; re-derive, never adjust.

## OWNER DECISIONS OUTSTANDING FROM BEFORE (recorded, not re-litigated)

Item 2's `g05 ss.99 ON` residue; item 5's +46 at `g05 OFF`; item 6's +1.2 % inside exons at `g50 ON`;
item 7's low-gDNA capture-ON rows on the benign and junction panels and `g50 ss.99 ON` benign/sparse.
The pre-walled panels are parked at `~/Downloads/rigel_runs/test_reference/pre_walled/` (14 GB).

## THE LESSONS THIS SESSION PAID FOR

* **A level carried between locales is refuted by probe placement alone, and the tool never sees the
  probe panel** — `DESIGN.md` §6b.9's founding refusal, re-measured at 4–17× on the adversarial panels
  after a clean stage 0 and a 9/10 win on the exon-probed one. Run all three panels before believing
  any mechanism that reads a rate across a face.
* **A count-based row has no resolution at a three-fragment exon, and the refit prior trains on its
  false modes** (`ISSUES: gdna-landscape-trains-on-false-positives`): the ladder's zero control turned
  a +13 % test-chromosome loss into 3.3×.
* **Today's `prepare` ORDER was an implicit propagation rule** (the walled block's terminus chains):
  a pass replaces it, and the per-slot identity gate is the instrument that finds such rules.
* **The relay's pooled transport-centre fit refuses on small substrates and misfires where it accepts**;
  it is not a drop-in for the transfer policy.
