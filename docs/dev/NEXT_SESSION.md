# NEXT SESSION — THE SHIP PROTOCOL IS DONE: `transfer` IS THE DEFAULT AND THE RELAY IS DELETED (2026-09-09); THE REFACTORING AND THE OPEN QUESTIONS ARE NEXT (handoff)

⭐⭐⭐ **THE REFERENCES:** `docs/dev/AMBIG_DESIGN.md` §4e (today's finding, the price of the two-sided hop,
every form with its numbers), `DESIGN.md` §6b.13 (the ruling of record), `docs/dev/MESSAGE_RUNGS.md` (the ship
protocol's checklist). This file is the state.

## WHAT STANDS (2026-09-09, end of session; the working tree is UNCOMMITTED — the owner drives commits)

1. **The base is committed**: `d68cd0bf` on `message-layer` (phases 1–2, the sj+terminus rule and block).
   Everything below is in the working tree.
2. **THE DEFAULT IS FLIPPED**: `CalibrationConfig.message_policy = "transfer"`. `rna_anchor` stays `True`
   but is inert (live iff the policy is relay). Judged on the ladder (`policy_benchmark.py --panel ladder
   --policies silent relay transfer`, 09:24 today):

   | | vs silence | worst row | vs relay |
   |---|---|---|---|
   | unstranded (must win) | 7/8 | 1.01× (`g05 ss.50 OFF`, +305 on 50,435) | wins 6/8 (the relay keeps `g00 ss.50 OFF/ON`) |
   | stranded (minimal harm) | 7/8 | 1.00× (`g05 ss.99 OFF`, +174 on 44,519) | wins 7/8 (the relay keeps `g00 ss.99 OFF`) |

   and on the 0.8.0 metric (`calibration_vs_oracle.py --message-policy <policy>`, the flag added today;
   mass-weighted |Δ gDNA| per object, region axis / boundary axis):

   | stratum | silent | relay | **transfer** |
   |---|---|---|---|
   | stranded × capture OFF | 0.0072 / 0.0138 | 0.0118 / 0.0154 | **0.0067 / 0.0104** |
   | stranded × capture ON | 0.0124 / 0.0218 | 0.0229 / 0.0295 | **0.0104 / 0.0149** |
   | unstranded × capture OFF | 0.0091 / 0.0169 | 0.0136 / 0.0181 | **0.0088 / 0.0124** |
   | unstranded × capture ON (deferred) | 0.5723 / 0.4818 | 0.1456 / 0.1146 | **0.0742** / 0.1166 |
   | g00 zero control | 0.0509 / 0.0320 | **0.0046 / 0.0059** | 0.0081 / 0.0077 |

   The relay's remaining lead is the zero-gDNA control, where the transfer policy leaves walled exons and
   both-stranded pieces to the landscape prior (`ISSUES: gdna-landscape-trains-on-false-positives`).
3. **LANDED TODAY, each gated (perturbations watched firing) and A/B'd on three panels and the ladder,
   halves apart, pass zero beside the pipeline:**
   * Two lane fixes: the gDNA lane's coordinate falls back to the library-wide density when the intergenic
     count is zero (a zero-gDNA library had no lanes at all: a both-stranded overlap read 0.95 gDNA against
     0); the RNA lanes' witness count and coordinate are read on the column κ points to (`read_column`).
   * ⭐ **THE SPLIT WITNESS** (`_RnaLane.witness`, `Level.rna_count` / `rna_count_var`): every RNA hop pays the
     pair's price — both column counts' counting plus the disagreement between the two nodes' estimates of
     the strand's abundance read from the column split's asymmetry, carried with the level across empty
     nodes; the column count where the library's strand channel is dead (`tau_lam` at single-strand exons,
     the deadband's verdict). The counting-only exemption on two-sided faces (landed 2026-09-08 on a
     measurement made with the broken witness column) carried a lit intron's upper side across a 170-fold
     probe cliff (the ladder's `g05 ss.99 ON`: a 93 % RNA junction read as 86 % gDNA). Ladder, full pipeline,
     against the exemption: stranded `g00 ss.99 ON` 0.767× (now below silence), `g05 ss.99 ON` 0.955×,
     `g50 ss.99 ON` 0.996×, `g98 ss.99 ON` 1.016×, capture-OFF rows within 0.1 %; unstranded `g00 ss.50 ON`
     0.706×, `g05 ss.50 OFF` 0.988×, in-scope rows within 1.2 %, deferred rows up to 1.037×; the three test
     panels within 0.1 % through the pipeline (the unstranded half is `pair2`'s form there, measured).
   * The strict xfail on `test_antisense_intronic::test_strand_sweep_with_nrna` removed (repaired by the RNA
     lanes); the toy harness's intron-independence gate made a strict xfail citing
     `ISSUES: two-sided-exon-row` (the wall that closes it is refused, below); the two anchor tests that
     assumed the relay default now name the relay; `calibration_vs_oracle.py --message-policy`;
     `config.py`'s propagation docstring cut to the shipped state; the goldens regenerated.
   **Suite: 0 failed / 3,783 passed / 7 xfail, 3,790 collected** (from 3,781 / 8 xfail: +1 gate in
   `test_transfer_policy.py`, the two-case antisense xfail now green, the toy gate xfailed; run on the final
   tree after the goldens' regeneration, and written into `CLAUDE.md`'s baseline line).
4. **REFUSED today, each with its numbers** (`docs/dev/AMBIG_DESIGN.md` §4e; `ISSUES.md`): the pair's price on
   the column counts (`pair2`: junction panel `g25 ss.99 ON` 1.240×, ladder `g98 ss.99 ON` 1.079×); the lower
   side everywhere and two-sided only at own junction bits (both ~1.21× on `g98 ss.99 ON`); a dilation of the
   whole profile by the totals' ratio (1.138×); a dark node's witness at its own noise (1.062×); THE WALL above
   rung 2's face-map ceiling (step F's cap re-derived: wins unstranded capture-OFF at pass zero 8–15 %, loses
   the stranded half 0/6 through the pipeline, the probe panels' capture-ON rows 1.4–2.5×); the flux price's
   witness in the strand's units (`ISSUES: flux-price-witness-units`: worse at pass zero, the split's
   asymmetry reads zero at an equal-abundance overlap exon).
5. **The goldens moved, recorded before `--update-golden`:** 18 of 21 scenarios by less than one fragment of
   gDNA; `antisense_contained_ss90` 0.04 → 51.9 false gDNA fragments of 1,000 (one both-stranded region with
   no witness of its second strand: the transfer policy equals silence there, the relay's reframe imputed the
   flanks' composition across the strand change); `antisense_overlap_ss90` 0.01 → 11.6; `strand_ss65_multi_iso`
   0.03 → 14.3 (a nested exon at κ = 0.31 whose flux ceiling reads f_g ≤ 0.38 where the flux says ≤ 0 —
   `ISSUES: flux-price-witness-units`); `combo_extreme` 321 → 264 gDNA with 120 → 181 nascent.
6. **`preflight.py --full` read 17/18 self-tests green** on the flipped tree; the one failure was
   `ladder_arm_ab.py`, a relay-era instrument whose arms rebind the relay's globals: pinned to
   `message_policy="relay"` (its docstring says so; it retires or is re-pointed in stage D), and its six
   flux-delivery arms (`msgfree_*`, `msgscale_*`) re-pointed from the retired `_PreparedRelay.deliver` to
   the two-phase `solve` — they had been dead since the backbone landed on 2026-09-04, which no `--full`
   run had been asked since. Its self-test rerun after the fix reads "every arm reaches the solver"
   (12:44), so `preflight.py --full` is 18/18 on this tree.
7. **The prototypes and A/B logs**: `~/Downloads/rigel_runs/prototypes/2026-09-09_ship/` (`two_sided_proto.py`
   with every arm, `wall_proto.py`, `fluxw_proto.py`, `cube_dump.py` — the cube delivery around one slot under
   two arms — `arm_vs_arm.py`, `pass0b.py`, `halves_pass0.py`, the `*_pass0_*.out` logs). ⛔ They subclass
   `TransferPolicy`; compare `src` against `src` from here.

## THE RETIREMENT, DONE THE SAME DAY (stages B–D; the owner's "proceed with the retirement stages")

8. **Deleted from `src/`**: `messages/relay.py` (1,303 lines), `messages/variance.py` (727; `count_logvar`
   moved to `messages/transfer_rows.py`, THE one home of the counting term, which `hop_price`,
   `count_price`, the transfer policy and `region_init` now read), `rna_anchor.py` (576), `RelaySwitches`,
   `config.rna_anchor`, the Gaussian channels on `PsiMessage` (gdna/rna/lam/theta mode + precision,
   `rna_one_sided`), their plumbing in `simplex_logodds.py` (`ONE_SIDED_RNA`, `_rna_residual`, every
   ``*_imp_*`` argument) and `sweep.py` (the two coordinate assertions and the share sum;
   `_KNOWN_VIOLATIONS` is empty). `PsiMessage` is two row channels, `lam_rows` and `cube_rows`. The
   backbone now stamps `policy_name` into the capture (the witness an instrument's "the arm ran"
   assertion reads) and publishes `order`/`n_slot`/`left`/`right` at the top level of the capture (the
   relay's `_uni_static` bank is gone).
9. **Deleted tests (8 files, 130 cases)**: `test_enrichment_frame`, `test_gdna_scale_rule`,
   `test_relay_mass_rescale`, `test_rna_anchor`, `test_splice_flux_reframe`,
   `test_terminus_population_licence`, `test_lambda_message`, `test_message_sidedness`. Re-pointed:
   `test_sweep.py` runs the sweep's shape under `SilentPolicy` (five relay-operator gates cut; the two
   Gaussian-message tests re-expressed as λ-rows), `test_sweep_backbone.py` (the channel and switch gates
   cut, a λ-rows shape/finiteness gate added), `test_region_init.py` (the relay's corner-variance xfail pair
   cut), `test_vertex_reference.py` (its claims delivered as λ-rows in each coordinate — every gate of the file holds).
10. **Deleted instruments (9)**: `ladder_arm_ab`, `arm_score`, `arm_sweep` (the relay-era arm harness and
    its readers), `reframe_walk`, `transfer_variance_audit`, `toy_ceiling`, `toy_dissect`,
    `toy_trace_error`, `psi_channel_ablation` (every one read the relay's per-slot banks or ablated its
    channels). Re-pointed: `backbone_parity.py` (arms `transfer` / `silent` / `module:<file.py>:<arm>` — the
    per-slot view of one prototype mechanism), `certified_rna_audit.py` (check (c) reads whether the
    shipped policy READS the bank — a flux level or a composition map — through a spy on `prepare`; it is
    green on the `tes_readthrough` rung and found `ISSUES: flux-source-skipped-at-an-empty-exon-piece`),
    `zero_controls.py` (the relay-derived column dropped), `calibration_walk.py` and `relay_pool_ab.py`
    (assert the policy off the capture's `policy_name`), `pass0_claimed_ab.py` (silent / transfer),
    `toy_harness.py` (the relay-bank helpers replaced by a policy stamp), `transport_dispersion.py`
    (keeps its own fifteen-line route table), `policy_benchmark.py` / `policy_prototype.py` (no relay
    arm, no `rna_anchor`), `benchmark_report.py`, `module_census.py`.
11. **Gated and measured**: `rename_identity.py --check` BIT-IDENTICAL on both frozen references
    (`g05 ss.99 ON`, `g05 ss.50 OFF`, at `~/Downloads/rigel_runs/arms/retire_identity_*.json`) after stage B,
    after stage C and on the final tree; the suite **0 failed / 3,595 passed / 4 xfail, 3,599 collected**
    (from 3,790 by −191, every case accounted in `CLAUDE.md`'s baseline line); `preflight.py --full`
    **16/16** self-tests, everything present. `DESIGN.md` §6b.3 is marked the record of a retired
    integration; `ISSUES: relay-od-r-discontinuity` closed with the relay; `EQUATIONS.md` §3.5 and §9d.1
    now say the transfer policy's map and lane are the built mechanism.
12. ⚠ **Naming debt the retirement leaves**: `relay_pool_ab.py` keeps its file name (its arms are OFF/ON of
    the shipped policy; a rename needs `rename_census.py --sense relay` first); `docs/dev/` still names the
    relay throughout (sandbox, historical); `EQUATIONS.md` §3.5–§3.6 derive the retired reframe and its
    licence — kept as the derivation the face map and the level lane satisfy, an owner call whether to
    move them under a record heading.

## WHAT IS NEXT — the owner's order (`ROADMAP.md` rank 1)

1. **The refactoring the retirement made possible**: `transfer.py`'s `prepare` as named rule builders
   (rung 1 / item 2 / rung 2 / item 1; rule 5; item 5 and the level rule; item 7; the gDNA lane; the RNA
   lanes) and ONE lane class for the gDNA and RNA lanes; `Message.tilt` (unused by construction); the
   split of `test_transfer_policy.py` (42 gates) by rung. Byte-identity gated by `rename_identity.py
   --check` on the two frozen references after each step.
2. **The open message-side questions, each its own step** — `ISSUES: flux-source-skipped-at-an-empty-exon-piece`
   (small; prototype on the terminus-cluster and sj+terminus blocks), `ISSUES: flux-price-witness-units`
   (the unstranded half's floors; needs a bounded witness at both-stranded exons),
   `ISSUES: two-sided-exon-row`, `ISSUES: flux-floor-dispersion`, `ISSUES: ambig-node-as-a-gdna-source`.
3. **The landscape's training population** (`ISSUES: gdna-landscape-trains-on-false-positives`) — the
   zero-gDNA rows the retired relay led on are its.
4. **The remaining both-stranded structures** (`div`, the antisense's nascent variant, `nest`).

## THE LESSONS THIS SESSION PAID FOR

* **A price chosen while a witness was broken is not a measurement.** The counting-only exemption rested on
  a run where every hop was priced on the wrong column; when the column was fixed, the exemption carried a
  false claim sharp. Re-measure every price that was chosen under a defect once the defect is repaired.
* **Two cliffs that look alike at the sender are told apart at the recipient.** The dark host intron beside
  a probed antisense exon and the lit intron beside a probed antisense exon send the same claim across the
  same jump in gDNA; only the recipient's own split says whether the strand's RNA is there.
* **A liveness test must be the channel's own.** `tau_lam` at an intron carries the factory's precision; the
  strand channel's verdict is `tau_lam` at single-strand exons. The prototype and the first landing both
  gated on the wrong one and every unstranded row ran on noise nets until the check was made.
* **The relay's remaining wins are one mechanism**: a two-sided reframe imputing the flanks' composition
  across a strand change, right whenever the truth is zero gDNA and unsound in general. The two golden
  scenarios that moved most are that mechanism, not two new defects.
* **Delete behind a frozen identity, module by module, suite after each stage.** Two references, one per
  half, made three deletion stages provable; the suite's collected count closed to the case on the table's
  rules. Nine instruments that read a retired policy's private banks were the real weight of the baggage,
  not the policy file.
