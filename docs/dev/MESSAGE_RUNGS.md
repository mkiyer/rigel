# THE MESSAGE RUNGS — the tracked list of messages and boundary cases still to finish (owner, 2026-09-02)

    ⚠ A DEV DOC and a TRACKER. It says what is done, what is next and in which order; the
    derivations and measurements live in `COMPOSITION_TRANSFER_STAGE01.md`, the substrate in
    `test_chr.yaml`'s header. Move a finished item's verdict to its permanent home and mark it here.

**The paradigm (owner rulings 2026-09-01/02):** one node type, one message, one boundary case at a
time; DERIVE → DESIGN → PLAN → PROTOTYPE (outside `src/`) → A/B on the test chromosome against
per-object truth, with the perturbation watched firing → confirm on the ladder, two halves apart,
adversarial probe panels in the loop → only then `src/`. A moved number must have ONE cause, so the
substrate grows ONE structure per step and a step is a whole session or more, validation included.
⛔ Nothing below is small; do not bundle two items to save a rebuild.

## Where the rungs stand

| rung | substrate | what it is | state |
|---|---|---|---|
| 1 | twin block | multi-exonic, single-isoform, single-stranded: the intron\|exon BOUNDARY | ✅ COMPLETE 2026-09-02 — intron → boundary (rung 1) and exon → boundary (item 1) both ship |
| 2 | twin block | the same substrate: the EXON region | ⚠ PARTIAL — the boundary → exon face map (rung 2) and the boundary → intron row (item 2, 2026-09-02) ship; item 3 (the exon solve with every face) is owed |
| 3 | mono block | single-exon transcripts: the intergenic\|exon EDGE | ✅ the sign-certified lower bound ships; the ceiling REFUSED (accepted error) |
| 4 | isoform block | multi-exonic, MULTI-isoform, single-stranded: exon\|exon boundaries | ⏳ substrate holds ONE structure (`altstart`); item 5 (the terminus boundary from the outside flank, spliced crossing included) SHIPS 2026-09-02; the chain-of-termini bound, `altss`, `nest`, `instart` remain |
| 5 | — | strand-change faces, both-stranded loci, the AMBIG tilt channel | ⛔ LAST, owner ruling |

## ⭐⭐⭐ THE COMPLETION CHECKLIST (owner ruling, 2026-09-02, `DESIGN.md` §0c.0e) — every case, its state

The policy is finished when every row below reads ✅. A row is one node type or one boundary case, with
the messages it must receive (IN) and send (OUT); "level" means the gDNA-level message (owed: the one
derivation the list hinges on); "composition" the face-map family (built). Ladder slot counts are one
condition's, from `policy_benchmark.py --by-class` and the terminus census (`DESIGN.md` §6b.6).

| row | node / boundary case | ladder slots | IN | OUT | state |
|---|---|---|---|---|---|
| intergenic-region | intergenic region | 1,312 | nothing (structural pure gDNA) | the edge bound (rung 3) | ✅ |
| gene-edge | intergenic\|exon edge | 2,620 | nothing (a count claim, locked) | the lower bound into the exon (rung 3) | ✅ (the ceiling REFUSED, accepted) |
| intron-single-strand | intron, single-strand | 9,805 | both faces' own strand rows (item 2) | the factory row (rung 1) | ✅ |
| intron-ambig | intron, AMBIG | (in intron-single-strand) | rung 1 (matching AMBIG sides) | the factory row | ✅ for λ; the tilt is ambig-tilt |
| intron-exon-face-sj | intron\|exon face, sj only | 15,480 | intron row (rung 1) + exon row (item 1) | to the intron (item 2), to the exon (rung 2) | ✅ |
| intron-exon-face-terminus | intron\|exon face WITH a terminus (the intron is the outside) | 4,130 | rung 1 (intron row) | item 2 to the intron; **level** into the inside exon | ⏳ level owed |
| exon-with-licensed-face | exon with a licensed intron face | 11,418 | rung 2 per face, rung 3, item 5's boundary row | item 1 to its faces | ✅; item 3 (the two-face sum re-priced) owed |
| exon-edge-only | exon, edge only (single-exon) | 1,250 | rung 3 | nothing | ✅ |
| exon-inside-or-walled | exon INSIDE a terminus, and the WALLED exon (no licensed face) | 11,350 | THE ABUNDANCE-DISCREPANCY message (`DESIGN.md` §6b.7): the boundary's own row through the map, the step's spread fitted | its own row back (owed) | ✅ item 6 SHIPPED 2026-09-02 for the boundary's own row (stranded data; ladder ≤ shipped 5/6, below silence 16/16); ⏳ owed to the scan phase with per-hop premises: the forwarded arrivals one hop further (the −2 % zero-control win, the +5 % unstranded-OFF harm) and the reverse direction (+0.8 %); walled exons beyond one hop likewise |
| exon-mixed-bits | exon with mixed exon+intron bits | (in exon-with-licensed-face/exon-inside-or-walled) | as its faces license | as its faces license | ✅ by ruling: the factory runs only where no exon bit is set (a mixed region is an exon) |
| terminus-outside-populated | exon\|exon TERMINUS, one direction, no sj, outside populated | 2,260 | outside exon's row + the composed transport (item 5) | its row to the outside exon (item 5); **level** into the inside | ✅ outside; ⏳ level inside |
| terminus-outside-empty | exon\|exon TERMINUS whose outside piece is EMPTY | 5,588 | the composite through the empty piece: from a licensed intron\|exon far face (716), an alt-ss (491), a gene edge (589); from ANOTHER terminus (3,792): **level** only | level onward | ⏳ the composite (multi-hop through an empty piece) + level |
| alt-splice-site | exon\|exon ALTERNATIVE SPLICE SITE (sj, no terminus) | 4,838 | the intron-side flank C shares the full crossing (item 5's map with `S_b`); the exon-of-both flank E holds the crossing plus T's mature RNA, MEASURED as the face's route flux `F` (item 1's map with `S_b + F`); every message carries THE HOP PREMISE — a step fitted per flank kind from the sighted pairs with its error, and the owner's discrepancy rule per pair | both flanks, both directions | ✅ item 7, LANDED 2026-09-02 (`DESIGN.md` §6b.8): ladder stranded 4 wins / 2 within 0.2 %, unstranded identical; sparse panel 8/8 ON; residues at `g05` capture-ON (+2.3…+5.6 %, mostly the refit prior's response) and the junction panel's `g50 ss.99 ON` (+3.4 %, the flux capture asymmetry at junction-probed panels) are OWNER DECISIONS |
| sj-plus-terminus | exon\|exon with sj AND terminus | 1,968 | the two maps composed | both | ⏳ item 10, after terminus-outside-populated/alt-splice-site |
| termini-both-ways | boundary with termini BOTH ways | 48 | **level** only | level | ⏳ with the level message |
| ambig-tilt | AMBIG boundary / the tilt channel | 4,671 | rungs 1–3 for λ | as its class | ⏳ rung 5: rule whether the tilt has any message (likely local-only, measured) |
| strand-change-face | strand-change face | 249 + 889 | **level** only (composition refused: membership changes) | level | ⏳ with the level message |
| multi-hop-scan | MULTI-HOP through any chain of the above | — | each hop's map + width + fitted premise; forward and backward; every node fused from two arrivals | — | ⏳ the transfer policy's `scan` (the backbone's two directional scans; one step kernel per boundary case) |

**The phases, in order.** ⭐ (A) is DONE as the owner's abundance-discrepancy rule (item 6, `DESIGN.md` §6b.7): the level message IS the composition map with the step between enrichment and new RNA, fitted. What (A) still owes: the same rule at strand-change faces and at termini both ways (the rule needs no orientation there; measure). (A, as first written) DERIVE THE LEVEL MESSAGE — gDNA continuity through an unmeasured
population change, two-sided, with counting width and a FITTED enrichment-step premise (stage 0: the
gDNA-density log-ratio across such faces on certified truth, its spread beyond counting, by capture,
range and probe state); prototype on `altstart`'s inside piece (exon-inside-or-walled/terminus-outside-populated), then intron-exon-face-terminus, termini-both-ways, strand-change-face. (B) alt-splice-site on
`altss` ✅ (its HOP PREMISE — a fitted step with its error plus the per-pair discrepancy — is the scan's per-hop
dampening template); sj-plus-terminus after it. (C) THE SCAN — the policy's forward/backward step kernels, one per boundary
case, which subsume terminus-outside-empty's composite, multi-hop-scan and item 3, with per-hop premises fitted and multi-hop
dampening measured on chains (a chain structure in the YAML). (D) ambig-tilt, and the both-stranded locus.
(E) The ship protocol (the SHIP LIST below).

## ⛔ THE PARKED-PRIORITIES LOG (owner, 2026-09-02: "keep a log, come back after the policy is done")

* ⭐⭐ `ISSUES: gdna-landscape-trains-on-false-positives` — the lever for the zero controls and the
  terminus/walled classes once the messages have done what they can (three exposures recorded there).
* The vertex atom / the intron's own solve off capture (`ROADMAP.md` rank 3; 42–46 % of in-scope
  off-capture error).
* `ISSUES: splice-out-premise-bias-uncorrected`; the item-2 `g05 ss.99 ON` prior-mediated residue.
* The fl-gap side panels' regeneration (`ISSUES: flgap-panels-stale-nascent-model`) — the only
  substrate that can price the per-component opportunity shift both directions of a hop carry.
* Promoting the session instruments (node-local scorer, the by-component pair-gap census).

## The ordered list (history of the items, in the order they were taken)

Finishing rungs 1 and 2 (owner: high priority, first):

| # | message / case | what must be derived | state |
|---|---|---|---|
| 1 | **exon → intron\|exon boundary** | DERIVED 2026-09-02 (owner's rescale/subtract/rescale ⇒ `f_b = f_E·(U_b+S_b)/U_b`, the enrichment ratio cancels; = the shipped splice-in map read backwards; checked unbiased on certified truth). Owed: the honest width (delta method through the map, diverging at `U_b → 0`), the imputation-cost premise priced on the adversarial panels, prototype + falsifiers, A/B | ⏳ PROTOTYPED 2026-09-02: wins every stranded/part-stranded capture-ON row (0.83–0.99×), exact silence elsewhere, both falsifiers fire; ONE benign cost (`g25 ss.99 ON` 1.014×, nascent-bearing probed boundaries) and sparse-panel harm (1.124×) = the premise cost, measured. Transfer variance DERIVED fresh: counting (1/S+1/U) + premise (log a, fitted by the two-witness estimator in log-ρ units — reads 0 off capture; under capture a BIAS a≈1.3 benign / 2.2 junction with ~no spread, recorded not corrected); width applied as a MARGINAL over log ρ; opportunities must be capture-blind for both components (geometric form). `formb_g` wins every stranded capture-ON row on all three panels. LADDER 2026-09-02: unstranded byte-identical, stranded capture-ON 0.987–0.995×, capture-OFF +8…+33 fragments, reversal fires. ✅ SHIPPED 2026-09-02 (`messages/transfer.py`, `simplex_logodds.strand_row_logodds`; ruling `DESIGN.md` §6b.4; 3 new gates, 3 perturbations watched firing) |
| 2 | **boundary → intron** | whether an intron with its own factory takes anything from its boundaries; possibly "nothing", measured | ✅ SHIPPED 2026-09-02: the boundary's OWN strand row VERBATIM (one shared population — no splice-out, no premise; the s = 0 map is the identity under the one-opportunity rule), deadband-gated at the boundary, AMBIG faces refused, terminus faces served (the intron is always the outside flank). Stage 0 certified the shared composition on all four substrates (zero excess variance over counting; the mean gap reproduced by a plug-in null). Node-local at the receiving introns −26…−58 % under capture, ±2 % off; the reversed row 1.7–30× worse everywhere. Ladder: unstranded byte-identical 8/8, stranded ≤ silent 7/8 (the `g05 ss.99 ON` 1.005× residue is prior-mediated — owner's call). Ruling `DESIGN.md` §6b.5; 3 gates fail-first, 4 perturbations watched |
| 3 | **the exon solve with every face speaking** | both faces' incoming rows plus the exon's own evidence, once item 1 exists on both sides; the two-witness sum re-priced | ☐ |
| 4 | the factory on a region carrying BOTH exon and intron bits (raised 2026-09-02: today the factory runs only where no exon bit is set) | whether the density-against-background measurement is valid there and under capture | ☐ when first needed |

⭐ ORDER (owner, end of 2026-09-02): item 2 ✅, then rung 4's terminus case (items 5–9), then item 3.

Rung 4, one structure per step (each step = one YAML change, one licence or message):

| # | boundary case | structure | what must be derived | state |
|---|---|---|---|---|
| 5 | **exon\|exon boundary WITH a terminus, solved from the OUTSIDE flank** | `altstart` (present) | the orientation table (TSS+/TES− body right, TES+/TSS− body left; mixed → no side), the opportunity shift, the solve; verified per object, then with orientation reversed | ✅ SHIPPED 2026-09-02 — ⚠ NOT verbatim: the licence counts the boundary's SPLICED crossing (`f_b = f_O (U_b+S_b)/U_b`, item 1's map; the verbatim form harmed +6…+10 % on the ladder and was refuted on certified truth). Three messages (own rows both ways through the map; the COMPOSED TRANSPORT of rungs 2–3's rows into the boundary, the no-echo law structural). Ladder: unstranded ≤ shipped 6/8 (0.985× at the g00 zero control), stranded ≤ shipped 7/8 (0.974× g98 ON), below silence everywhere; small at the destinations because they already hold own evidence + the prior. Ruling `DESIGN.md` §6b.6 |
| 6 | the region INSIDE the terminus | `altstart` | what, if anything, crosses into it (a level bound was derived and REFUTED on the sparse panel — record, do not re-litigate without a new measurement) | ☐ |
| 6b | **the CHAIN of termini** — an empty outside piece (median 12 bp) whose far face is another terminus: HALF the ladder's terminus-boundary error (`DESIGN.md` §6b.6's census) | needs a structure (two TSSs a dozen bases apart inside an exon) | composition reaches none but the outermost; derive the short-range gDNA-LEVEL continuity under one probe footprint as an UPPER bound on the inner boundary's gDNA share; price it against the prior, which already serves these slots to ~1.8 % of mass | ☐ after the landscape prior |
| 7 | exon\|exon boundary at an alternative SPLICE SITE | add `altss` | the crossing loses the transcript that splices out (join vs leave direction); the composed transport was prototyped and measured neutral on the ladder | ☐ |
| 8 | a region WALLED by two termini | add `nest` | what reaches it; measured 2026-09-02: the refit prior already serves it (0.666 vs 0.630) | ☐ |
| 9 | an alternative FIRST exon inside an intron (`exon\|intron[term]`) | add `instart` | the intron is always the outside flank; the inside exon's other face | ☐ |
| 10 | terminus faces with sj+term flags, chains of termini | (ladder only so far) | after 5–9 | ☐ |

Rung 5 (last): strand-change faces (246 exon|exon + 888 exon|intron on the ladder), the both-stranded
locus, the AMBIG tilt channel.

## THE SHIP LIST (2026-09-02, after item 2) — what stands between `transfer` and the default

Re-derive every number here with `policy_benchmark.py --panel ladder --policies silent transfer
--by-class` (and `--panel test`); the shares are of each row's own |err| under `transfer`.

**Where the shipped policy's remaining error sits, by node class (ladder):**

| row | exon\|exon (terminus + alt-ss) | R exon WALLED | B exon\|intron (+ terminus) | R exon licensed / edge-only | R intron |
|---|---|---|---|---|---|
| g50 ss.99 OFF (in scope) | 21 % | 7 % | 20 % | 7 % | 45 % |
| g98 ss.99 OFF (in scope) | 23 % | 9 % | 18 % | 7 % | 43 % |
| g05 ss.99 ON (in scope) | 39 % | 18 % | 30 % | 11 % | 3 % |
| g98 ss.99 ON (in scope) | 47 % | 14 % | 26 % | 10 % | 3 % |
| g50 ss.50 OFF (in scope) | 20 % | 8 % | 18 % | 9 % | 46 % |
| g00 ss.50 OFF (in-scope zero control) | 51 % | 35 % | 1 % | 13 % | 0 % |
| g50 ss.50 ON (deferred) | 44 % | 16 % | 29 % | 11 % | 1 % |

Reading it: (i) **exon|exon boundaries receive NO message today** — the policy's face loop requires an
intron flank — and with the exons they wall they hold 62 % of the stranded capture-ON error and 85 %
of the zero-control error; the 5× gap to the relay at `g00 ss.50 OFF` (303,826 vs 58,840) sits there.
Rung 4 is the hole, as the owner ordered. (ii) **Introns off capture (42–46 %) are NOT a message
hole**: the intron's factory row is one-sided and its strand row vanishes at the pure-gDNA vertex, so
the prior decides there — the vertex-atom / prior thread (`ROADMAP.md` rank 3, after the message
layer). The only message that could reach them is the exon → boundary → intron TWO-HOP path (item 1's
row carried on by item 2's identity) — the COMPOSED TRANSPORT of rung 4 (ii), which the owner
confirmed (2026-09-02) is the direction: multi-hop propagation with each hop's own dampening, never a
sudden switch of paradigm; ⚠ the rebuild's "one hop" is its current STRUCTURE (`scan` returns None),
not a ruling. (iii) exon|intron boundaries (18–30 %) are served by rung 1 + item 1; their residue is their
own thin evidence and the prior. (iv) On the test chromosome the off-capture rows are the designed
shadow-transcription floor (82–100 % at `R intergenic`, identical in every arm) and the on-capture rows
sit at the exons' own strand resolution: the twin block cannot indict rung 4, so the isoform block must
grow (`altss`, `nest`, `instart`), one structure per step.

**What ships today, per node type and message — the accounting:**

| node type | own solve | messages IN | messages OUT | the hole |
|---|---|---|---|---|
| intergenic region | structural (pure gDNA) | none | rung 3 | none in scope |
| intergenic\|exon edge | count claim | none | rung 3's lower bound → exon | the ceiling REFUSED (accepted error) |
| intron region | factory + strand | **item 2** (both faces' own strand rows) | rung 1 (the factory row) | item 4 — the factory runs only where no exon bit is set (mixed exon+intron regions) |
| intron\|exon boundary | strand | rung 1 + item 1 | item 2; rung 2 (through the face map) | terminus faces: served boundary-side (the intron is the outside flank), refused exon-side (the inside; the bound REFUTED) |
| exon region | strand | rung 2 (each licensed face) + rung 3 (edge) | item 1 (to licensed faces) | item 3 (the solve with every face, re-priced); the WALLED exons (rung 4) |
| exon\|exon boundary | strand | **NONE** | **NONE** | rung 4, items 5–10: terminus 7,973 slots (0.99 M mass OFF / 2.2 M ON), alt-ss 4,838 (0.46 M) |
| strand-change faces | strand | none | none | rung 5: 249 exon\|exon + 889 exon\|intron slots, 0.3 % of mass |
| AMBIG slots (both strands admitted) | λ from messages only (the strand row is Schur-cancelled) | rungs 1–3 for λ; NO tilt channel | as their class | rung 5: 5,241 regions + 4,671 boundaries, 6–13 % of mass; the tilt is local-only |

**The order (owner, unchanged): items 5–9 → item 3 → item 4 when first needed → rung 5 → the flip.**
Each rung-4 item is one YAML structure (`panel.py` rebuilds all seven panels in ~40 min), one licence or
message, prototyped through `policy_prototype.py`, scored NODE-LOCALLY at the slots it serves and
whole-library, falsified by reversing its orientation, laddered both halves apart.

**The ship protocol, once the holes are closed or ruled accepted:**

1. Price the rebuild on the 0.8.0 metric — `calibration_vs_oracle.py` and `solvability_audit.py`
   under `transfer`, per stratum (never yet done; the policy benchmark is the fragment-count view).
2. Flip `CalibrationConfig.message_policy` to `"transfer"`; `rna_anchor` (live iff relay) goes with it.
3. The obsolescence pass, converge-and-delete, one owner ruling per module: `relay.py` (1,254 lines),
   `variance.py` (727, its toolbox), `rna_anchor.py` (the certified-flux stream — the one relay
   capability the rebuild never re-derived; rung 1's stage 1 measured it wanted at `g98 ss.99 ON`
   node-locally), `policy.py` + `foundation.py` (the rung-0 skeleton). 15 `scripts/design/`
   instruments and 8 test files name the relay (`grep -rl RelayPolicy`); each is re-pointed or deleted.
4. `preflight.py --full` (a default flip is one of its three triggers), `--update-golden` with the diff
   read and its magnitude recorded first, `docs/MANUAL.md`, and CLAUDE.md's message-layer section
   rewritten to the shipped state.

**Decisions owed to the owner** (each recorded where it lives): the prior-mediated `g05 ss.99 ON`
1.005× residue (`DESIGN.md` §6b.5); `ISSUES: splice-out-premise-bias-uncorrected`; the zero-gDNA edge
residual (accepted); whether the exon → intron two-hop is a paradigm exception; and
`ISSUES: gdna-landscape-trains-on-false-positives`, exposed again by item 2.

## Recorded and set aside (do not rebuild)

* The rung-4 prototype (`rung4_proto.py`, session scratchpad 2026-09-02): composed join-only transport
  NEUTRAL on the ladder (worst 1.006×/1.008×, −6 % at `g00 ss.50 OFF`); the inside bound REFUTED
  (in-scope +1.3 %, sparse probes +8 %); flips fire. Its pieces re-enter only as items 5–8 earn them.
* The rung-4 census instruments (an exon|exon boundary census by flag class, a directed
  reachability census, a terminus pair-gap measurement on certified truth) lived in a session scratchpad
  and are gone; their findings are in the thread record's RUNG 4 sections and are re-derivable from
  `slot_truth.npz` + the chain in an afternoon. The generic prototype harness survived as
  `scripts/design/policy_prototype.py`.
* `panel.py cache` now pre-warms the g00 rows, copies `_main` and certifies (2026-09-02).
