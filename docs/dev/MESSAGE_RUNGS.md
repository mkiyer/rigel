# THE MESSAGE RUNGS — the tracked list of messages and boundary cases still to finish (owner, 2026-09-02)

    ⚠ A DEV DOC and a TRACKER. It says what is done, what is next and in which order; the rulings
    and measurements of everything landed live in `DESIGN.md` §6b.4–§6b.9, the open working notes in
    `COMPOSITION_TRANSFER_STAGE01.md`, the substrate in `test_chr.yaml`'s header. Move a finished
    item's verdict to its permanent home and mark it here.

**The paradigm (owner rulings 2026-09-01/02):** one node type, one message, one boundary case at a
time; DERIVE → DESIGN → PLAN → PROTOTYPE (outside `src/`) → A/B on the test chromosome against
per-object truth, with the perturbation watched firing → confirm on the ladder, two halves apart,
adversarial probe panels in the loop → only then `src/`. A moved number must have ONE cause, so the
substrate grows ONE structure per step and a step is a whole session or more, validation included.
⛔ Nothing below is small; do not bundle two items to save a rebuild.

## Where the rungs stand

| rung | substrate | what it is | state |
|---|---|---|---|
| 1 | twin block | multi-exonic, single-isoform, single-stranded: the intron\|exon BOUNDARY | ✅ COMPLETE — intron → boundary (rung 1), exon → boundary (item 1), boundary → intron (item 2) ship (`DESIGN.md` §6b.9, §6b.4, §6b.5) |
| 2 | twin block | the same substrate: the EXON region | ⚠ PARTIAL — the boundary → exon face map ships (§6b.9); item 3 (the exon solve with every face speaking) is owed and is subsumed by THE SCAN |
| 3 | mono block | single-exon transcripts: the intergenic\|exon EDGE | ✅ the sign-certified lower bound ships; the ceiling REFUSED (accepted error; §6b.9) |
| 4 | isoform block + WALLED block | multi-exonic, MULTI-isoform, single-stranded: exon\|exon boundaries, and exon pieces with NO licensed face | ⚠ PARTIAL — items 5, 6 (`altstart`; §6b.6–§6b.7) and 7 (`altss`; §6b.8) ship, every message priced by the owner's per-pair discrepancy rule; the WALLED BLOCK (`chain` · `tssalt` · `tandem` · `altlast`, 48 genes, owner-approved 2026-09-03) is the scan's substrate; `nest`, `instart`, the chain of termini beyond it and sj+terminus remain |
| — | every block | THE SCAN — multi-hop through every case above | ✅ LANDED 2026-09-04 as the two-phase backbone's passes (`DESIGN.md` §6b.12); the seam of §6b.10 is superseded. Owed: the two-sided exon profile (`ISSUES: two-sided-exon-row`), which is what would make the unstranded probe-panel rows win |
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
| alt-splice-site | exon\|exon ALTERNATIVE SPLICE SITE (sj, no terminus) | 4,838 | the intron-side flank C shares the full crossing (item 5's map with `S_b`); the exon-of-both flank E holds the crossing plus T's mature RNA, MEASURED as the face's route flux `F` (item 1's map with `S_b + F`); every message carries the owner's DISCREPANCY RULE PER PAIR (widened by its own pair's disagreement beyond counting; nothing pooled — the per-library step REFUSED 2026-09-03) | both flanks, both directions | ✅ item 7, LANDED 2026-09-02 (`DESIGN.md` §6b.8): ladder stranded 4 wins / 2 within 0.12 %, unstranded identical; sparse panel 7/8 ON; residues at the low-gDNA capture-ON rows of the benign and junction panels (`g05`, `g25`: +1.5…+9 %, an offset below each pair's counting that only pooling would see) and `g50 ss.99 ON` benign/sparse (+2.0 / +2.4 %) are OWNER DECISIONS |
| sj-plus-terminus | exon\|exon with sj AND terminus | 1,968 | the two maps composed | both | ⏳ item 10, after terminus-outside-populated/alt-splice-site |
| termini-both-ways | boundary with termini BOTH ways | 48 | **level** only | level | ⏳ with the level message |
| ambig-tilt | AMBIG boundary / the tilt channel | 4,671 | rungs 1–3 for λ | as its class | ⏳ rung 5: rule whether the tilt has any message (likely local-only, measured) |
| strand-change-face | strand-change face | 249 + 889 | **level** only (composition refused: membership changes) | level | ⏳ with the level message |
| multi-hop-scan | FORMAL FORWARD-BACKWARD through every case above | — | ⭐ the owner's ruling 2026-09-04 (`DESIGN.md` §6b.11): the sender just sends; the RECIPIENT decides to FORWARD, MODIFY or STOP; when the two passes end, EVERY node holds a message from each neighbour (the chain's two end nodes hold one) — a hop that carries nothing still ARRIVES, explicitly uninformative | — | ✅ LANDED 2026-09-04 (`DESIGN.md` §6b.12): the two-phase backbone (`sweep._pass`, `messages/__init__.py`'s `Message`/`SILENCE`/`NO_NEIGHBOUR`/`Prepared`), the foundation scaffold retired, the relay byte-identical, and `messages/transfer.py` rebuilt as claims + rules + the passes (bit-identical to the prototype's formal arm; gated on an independent recursive reference of the passes, the no-echo perturbation, and a claim-and-rule gate per family). Measured: wins BOTH halves of the ladder against the one-hop policy 7/8 + 7/8 (pass zero and full); probe panels stranded within 4–8 %, unstranded mixed — the forwarded exon profile is one-sided, `ISSUES: two-sided-exon-row`. ⛔ The prerequisite as first proposed (the certified-flux row) is REFUTED as a level (`ISSUES: the-certified-flux-row-as-a-level`) |

## ⭐⭐⭐ THE TEN RULES, BY CURRENCY (2026-09-04) — what each directed face does today

A rule is what the RECIPIENT does with what arrives; every rule below is a function of the sender's
claim and the two nodes' observations, never a belief. "Composition" = the gDNA-vs-RNA profile crosses
by a map; "level" = a gDNA abundance crosses and is converted at the recipient through its own total.

| # | face (source → destination) | rule | currency | ingredients | owner ruling |
|---|---|---|---|---|---|
| 1 | intron → intron\|exon boundary | FORWARD (identity) | composition | the intron's factory profile; one shared unspliced population | rung 1, `DESIGN.md` §6b.9 |
| 2 | intron\|exon boundary → intron | FORWARD, if the pair shares ONE strand | composition | the boundary's strand profile | item 2, §6b.5 |
| 3 | boundary → exon at a LICENSED face | the splice-in face map + the face's counting width | composition (the flux enters the MAP as a measured spliced density) | crossing count, both opportunities, the route rate | rung 2, §6b.9 — ⚠ one-sided (`ISSUES: two-sided-exon-row`) |
| 4 | exon → boundary at that face | the splice-in map read backwards, marginalised over the face's spliced/unspliced ratio | composition | the exon's strand profile | item 1, §6b.4 — the premise bias under capture recorded, not corrected |
| 5 | intergenic\|exon edge → exon | THE EDGE'S LEVEL, one-sided: the exon has at least the edge's gDNA density, at the count's Poisson width; nothing above; a zero count vacuous (darkness under capture is not absence) | **level** | the edge's count, the opportunity ratio, the exon's total | `MESSAGE_PLAN.md` step E, 2026-09-04: every upper side REFUSED by the sparse panel and the ladder (`ISSUES: the-edge-upper-side`); rung 3's formula stands, re-derived |
| 6 | outside exon → exon\|exon TERMINUS boundary | splice-out with S = the boundary's SPLICED crossing | composition | the outside exon's claim (and what it holds from beyond) | item 5, §6b.6 |
| 7 | terminus boundary → outside exon | the face map with the spliced density | composition | the boundary's strand profile | item 5, §6b.6 |
| 8 | terminus boundary → the region INSIDE (exon\|exon AND exon\|intron) | THE LEVEL RULE: the boundary's OWN strand profile through the level-kept map (its share × its crossing density × the inside's opportunity / the inside's own total), shape preserved, blurred by both totals' counting and the pair's own discrepancies (totals; strand modes where live); with no own claim, the crossing total's upper bound | **level** (a measurement only; what is held never crosses) | the boundary's own claim, both totals, the spliced crossing | `MESSAGE_PLAN.md` step A, landed 2026-09-04; item 6's map, cap and pooled spread deleted (`DESIGN.md` §6b.7 superseded) |
| 9 | flank → ALTERNATIVE SPLICE SITE boundary | splice-out with S_b (+ the route flux F on the exon-of-both flank), blurred by the pair's own disagreement | composition | the flank's strand profile | item 7, §6b.8 |
| 10 | alt-ss boundary → flank | the face map with the matching spliced density, blurred by the pair's disagreement | composition | the boundary's strand profile | item 7, §6b.8 |

**The missing rules (STOP by omission), and what exists for each:** STRAND-CHANGE faces (1,138) and termini BOTH WAYS (48) — rule 8's map is the candidate, the
totals' comparability at those faces is the derivation owed; the EMPTY chain pieces (5,588 terminus
boundaries, 3,792 reached only through another terminus) — no implementation: the node has no total to
convert with, so the gDNA level must be FORWARDED (the `level_gdna` lane's first real use) and converted
at the next node that has one; SJ + TERMINUS faces (2,354) — the constructors exist, the licence (which
flank shares which population when a junction and a terminus share a face) is the derivation owed;
the AMBIG tilt and the both-stranded locus — no substrate, no rule, a ruling first.

## ⭐⭐⭐ THE SHIP AUDIT (2026-09-04, after the two-phase backbone landed) — what stands between `transfer` and the flip

The bar (owner, 2026-09-04): an INTACT design for every node type and face before the policy ships —
no case may be a STOP by omission — and then the standings. Under the landed policy a face with no rule
is a hole; a node with no claim of its own on unstranded data is the other kind of hole. Ladder error
shares are the landed policy's (`policy_benchmark.py --panel ladder --by-class`).

| # | case | kind of hole | ladder slots | error share (in scope) | what it needs | state |
|---|---|---|---|---|---|---|
| 1 | **THE LEVEL RULE** — the inside of every terminus (✅ landed 2026-09-04 as rule 8, both terminus kinds), then strand-change faces, termini both ways, the EMPTY chain pieces between termini | STOP by omission (no rule) at the remaining faces | 1,138 + 48 + 5,588 boundaries (3,792 reached only through another terminus) | terminus boundaries 12–15 % OFF / 28 % ON; walled exons 8–10 % / 18 % | the landed rule's laws (a level from the measurement only, shape preserved, per-pair widths, the total's upper bound) extended to the faces where composition cannot cross; the empty piece needs the level LANE forwarded | ⏳ the inside of termini ✅; the rest owed |
| 2 | **sj + terminus** faces | STOP by omission (`outside_flank`/`junction_flanks` return none) | 1,968 + 386 | inside the terminus / alt-ss shares | item 5's and item 7's maps COMPOSED into one rule | ⏳ |
| 3 | **the two-sided exon profile** (unstranded data) | no CLAIM at exons; the forwarded profile is one-sided | 11,418 licensed + 11,350 walled exons | licensed 6–10 %; feeds walled/terminus via forwarding; the prior's training population at pass zero | the intron's own high-side sharpness, or a bounded taper marginal in the face map; ⛔ not the flux level | ⏳ `ISSUES: two-sided-exon-row` |
| 4 | **AMBIG tilt** and the both-stranded locus | a ruling, not a hole (rung 5) | 4,671 + 5,241 | 6–13 % of mass | rule whether any message exists; likely local-only | ⏳ last by owner ruling |
| 5 | **substrate**: `nest`, `instart` | nothing to falsify 1–2 on | — | — | one YAML structure per step | ⏳ |
| 6 | the INTRON's own solve | not a message | 9,805 | **43–45 % off capture**, identical under every policy | `ROADMAP.md` rank 3 (parked: after the policy ships) | parked |
| 7 | THE SHIP LIST | — | — | — | the 0.8.0-metric pricing under `transfer` (`calibration_vs_oracle.py`, `solvability_audit.py`), the flip, the obsolescence pass (`relay.py` 1,254 lines, `variance.py` 727, `rna_anchor.py` 576, 15 instruments and 8 test files naming the relay), `preflight --full`, goldens, `MANUAL.md`, `CLAUDE.md` | after 1–4 |

**THE ORDER (owner review, 2026-09-04 — `MESSAGE_PLAN.md` carries the designs).** (A) rule 8
re-specified as THE LEVEL RULE for the inside of every terminus (it is misspecified: a pooled spread, a
hypothesis mix, the wrong currency, half its faces missing) → (E) rule 5 as a level, one-sided vs
two-sided A/B'd (a zero-count edge is vacuous today: the relay's zero-control lead) → (F, G) rules 3 and
4 refined (the wall from the intron's level, the per-pair discrepancy) → (B) strand-change faces and
termini both ways → (C) the empty chain pieces → (D) junction plus terminus → (H) the AMBIG ruling →
the ship list. Every case: the simplest LOCAL form, nothing pooled; stage 0 on certified truth;
prototype through `policy_prototype.py`; judged at its DESTINATIONS, halves apart, PASS ZERO beside the
full pipeline, all three probe panels, then the ladder; then `src/` with fail-first gates and every
perturbation watched. ⛔ Residue decisions are deferred until no rule is missing.

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
| 6 | the region INSIDE the terminus | `altstart` | what, if anything, crosses into it (a level bound was derived and REFUTED on the sparse panel — record, do not re-litigate without a new measurement) | ✅ SHIPPED 2026-09-02 — THE OWNER'S ABUNDANCE-DISCREPANCY RULE: the boundary's own row travels `f_X = f_c·s/r`, the step `s` between enrichment and new RNA, never above `f_c`, its spread fitted from the served pairs' two witnesses; ladder stranded ≤ shipped 5/6, unstranded byte-identical; +1.2 % inside exons at `g50 ON` is an owner decision. Ruling `DESIGN.md` §6b.7 |
| 6b | **the CHAIN of termini** — an empty outside piece (median 12 bp) whose far face is another terminus: HALF the ladder's terminus-boundary error (`DESIGN.md` §6b.6's census) | needs a structure (two TSSs a dozen bases apart inside an exon) | composition reaches none but the outermost; derive the short-range gDNA-LEVEL continuity under one probe footprint as an UPPER bound on the inner boundary's gDNA share; price it against the prior, which already serves these slots to ~1.8 % of mass | ☐ after the landscape prior |
| 7 | exon\|exon boundary at an alternative SPLICE SITE | `altss` (added under the REPLICATION RULE) | the crossing loses the transcript that splices out (join vs leave direction) | ✅ SHIPPED 2026-09-02 — both licences certified on the ladder at every flank length; every message carries the owner's DISCREPANCY RULE PER PAIR, nothing pooled (the per-library step landed for a day and was REFUSED as over-engineering, `ISSUES: the-pooled-hop-step`); ladder stranded 4 wins / 2 within 0.12 %, unstranded identical; the low-gDNA capture-ON rows of the benign and junction panels are owner decisions. Ruling `DESIGN.md` §6b.8; the refused pooled step and flux-factor forms in `ISSUES.md` |
| 8 | a region WALLED by two termini | add `nest` | what reaches it; measured 2026-09-02: the refit prior already serves it (0.666 vs 0.630) | ☐ |
| 9 | an alternative FIRST exon inside an intron (`exon\|intron[term]`) | add `instart` | the intron is always the outside flank; the inside exon's other face | ☐ |
| 10 | terminus faces with sj+term flags, chains of termini | (ladder only so far) | after 5–9 | ☐ |

Rung 5 (last): strand-change faces (246 exon|exon + 888 exon|intron on the ladder), the both-stranded
locus, the AMBIG tilt channel.

## THE SHIP LIST (re-stamped 2026-09-03, after item 7) — what stands between `transfer` and the default

Re-derive every number here with `policy_benchmark.py --panel ladder --policies silent relay transfer
--by-class` (and `--panel test`); the shares are of each row's own |err| under `transfer`.

**Where the shipped policy's remaining error sits, by node class (ladder, 2026-09-03):**

| row | exon\|exon (terminus + alt-ss) | R exon WALLED | B exon\|intron (+ terminus) | R exon licensed / edge-only | R intron |
|---|---|---|---|---|---|
| g50 ss.50 OFF (in scope, unstranded) | 20 % | 8 % | 18 % | 9 % | 46 % |
| g98 ss.50 OFF (in scope, unstranded) | 23 % | 9 % | 17 % | 8 % | 44 % |
| g50 ss.99 OFF (in scope) | 21 % | 7 % | 20 % | 7 % | 45 % |
| g98 ss.99 OFF (in scope) | 23 % | 8 % | 19 % | 7 % | 43 % |
| g05 ss.99 ON (in scope) | 40 % | 18 % | 28 % | 11 % | 3 % |
| g98 ss.99 ON (in scope) | 47 % | 14 % | 27 % | 10 % | 3 % |

⭐ Read with the g00 census (`DESIGN.md` §6b.9's zero-control note): at the in-scope zero control the
transfer policy's error sits at walled exons (35 %), exon|exon terminus boundaries (31 %) and alt-ss
boundaries (19 %) — the objects the relay's anchor pins and the transfer policy leaves to the prior.
⚠ The own-row messages of items 5–7 serve the exon|exon boundaries but barely move them on the
contaminated rows (g98 ss.99 ON: terminus boundaries 79,582 silent → 72,812; alt-ss 55,906 → 51,053):
they sit at the resolution of their own evidence plus the prior. The scan is what changes the walled
exons' and the zero controls' standing; the boundaries' remaining share is the prior's question
(`ISSUES: gdna-landscape-trains-on-false-positives`, parked behind the policy).
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
