# NEXT SESSION — PRIORITY 1 (the cleanup) IS DONE, THE EMPTY-PIECE FLUX SOURCE IS LANDED AND `vertex_ceiling.py` RE-POINTED, ALL UNCOMMITTED; PRIORITIES 2–4 STAND (handoff, 2026-09-09)

⭐⭐⭐ **THE REFERENCES:** `DESIGN.md` §6b.4–§6b.14 (every ruling of the shipped policy, each with its
measurement), `ISSUES.md` (every open case and every refusal with its number), `CLAUDE.md`'s message-layer
section and instrument table. This file is the state.

## WHAT STANDS (2026-09-09; the working tree on `message-layer`, on top of `593ba28a`)

1. **UNCOMMITTED**: the cleanup session's five stages are in the working tree (55 paths: `git status`).
   The owner drives commits. Every stage was gated before the next began: the suite green, ruff clean,
   and `rename_identity.py --check` BIT-IDENTICAL on both frozen references
   (`~/Downloads/rigel_runs/arms/retire_identity_*.json`: `g05 ss.99 ON`, `g05 ss.50 OFF`) after each of
   the five stages and on the final tree. `preflight.py --full` on the final tree: everything present,
   **16/16** instrument self-tests. `module_census.py`: no upward import, no dead public surface.
2. **THE SUITE: 0 failed / 3,592 passed / 2 xfail, 3,594 collected** (`CLAUDE.md`'s baseline line carries the
   accounting from 3,599: −7 test functions retired with the relay's own-precision arithmetic, +14 for the
   seven `tests/calibration/` files the transfer gate file became, −9 for the nine `docs/dev/` files moved
   out or deleted, −3 with `RegionInit.struct_lock`). The 2 xfails: the toy harness's intron-independence
   gate and the antisense t2 prior-assembly casualty.
2b. **THE OWNER'S THREE DECISIONS (2026-09-09, end of the session), applied and gated bit-identical:**
   `docs/dev/message_notes.md` deleted; the plan's rung/item labels removed from `DESIGN.md` §6b.4–§6b.9
   (titles and bodies name the messages: the intron's forward, the splice-in map, the splice-out map, the
   edge's bound, the terminus's outside map, the level rule, the alternative splice site); and
   `RegionInit.struct_lock` DELETED with its xfail pair — nothing on any path read it (the sweep and the
   policy read `f_*` and `tau_lam`; every instrument computes its own certainty mask from
   `region_geometry.g1_locked`, the one predicate), so `strand_evidence` returns `I_strand` alone and
   `RegionInit` is four fields. ⚠ FOUND on the way: **`vertex_ceiling.py`'s pin is INERT under the transfer
   policy** — it rewrites `RegionInit.f_*`, which reach only the diagnostics capture (`fg_loc`); the sweep
   solves from the incoming belief and ψ, and the policy reads `tau_lam` only. Its comparator would report
   the vertex arm "byte-identical, did not fire". Before `ROADMAP.md` rank 6 (re-price the vertex atom) the
   instrument must pin where the solve reads: the incoming `belief` handed to `solve_chain`, or ψ's own
   evidence rows.
3. **Nothing any message computes changed** — that is what the identity gate proves per stage. Every
   refactor below is bit-identical on both halves; anything that would have moved a number is listed
   under "FOUND, NOT CHANGED" for priority 2.

### Stage 1 — `messages/transfer.py` (904 → 949 lines, one structure)
* `prepare` is a table of contents: one named BUILDER per shipped message — `_claims` (every node's own
  claim), `_splice_faces` (the intron|exon face: forward both ways, the splice-in map, the splice-out map),
  `_edge_level`, `_terminus_rules` (the outside map and THE LEVEL RULE, the sj+terminus placement),
  `_alternative_splice_site`, `_gdna_lane`, `_rna_lanes` — each a function of `_Chain`, the context
  unpacked once (with `strand_profile` and `strand_mode`, the two strand witnesses the builders share).
  The module docstring lists the messages by builder name.
* `_Lane` and `_RnaLane` are ONE class, `_LevelLane(population, ...)`: the lanes differ by the WITNESS
  the hop price reads (`count`/`a`, plus `other` where the strand channel is live) and by `two_sided`
  (empty on the gDNA lane) — parameters, not classes. `emit(s, x, far)`, `receive(level, s, x)`,
  `row(profile, x)` (a held level as the node's composition row: `profile_of_level` for gDNA,
  `rna_row_of_level` for a strand). `_PreparedTransfer.lanes` is a dict keyed `"gdna"`/`"pos"`/`"neg"`;
  `_SolveSite` is what the two RNA deliveries read at a node. The propagate kernel loops over the lanes.
* `Message.tilt` is deleted (four lanes; the tilt has no lane — the two RNA levels constrain it in the
  cube). `transfer_rows.count_price` folded into `hop_price` (one formula; the zero-count guard kept —
  bit-identical on the gDNA lane, whose counts are always positive); `priced_level` deleted (the lane's
  `receive` spells it: the lower side unless the face is two-sided, then the hop's blur).
* The plan's rung/item numbers are gone from the source and the gates: a message is cited by its name.

### Stage 2 — the relay's plumbing
* `StepContext` is the 24 fields the policy and the backbone read (from 41): deleted `mass`,
  `inv_abundance`, `inv_sj_lo/hi`, `eff_gdna` (the per-face one; `eff_gdna_global` RENAMED `eff_gdna`),
  `eff_sj`, `route_count_lo/hi`, `left/right_interface_certified`, `ss_intron_boundary`, `geometry`,
  `order`, `left_list`/`right_list` (the backbone's own lists, now locals), `solve_grid`, `capture`. Its
  docstring's "measured debt" paragraph (the relay's reframe reading `belief_fg` at both ends) is gone:
  the shipped policy reads `belief_fg` once, at `prepare`, source-side by construction.
* `RegionInit` is `(f_g, f_pos, f_neg, struct_lock, tau_lam)`: `rho_*`/`prec_*` (the relay's own
  densities and precisions), `own_precision`, `own_composition_logvar` and `region_init`'s
  `region_gdna_geometry`/`count_logvar` imports deleted; the location-term refusal comment shrunk to the
  ruling (its numbers are `ISSUES.md`/`TRAPS.md`'s). `has_own_composition_evidence` stays (the
  instruments' predicate).
* `region_geometry.terminus_flank_gain` (the population test in genomic terms — `transfer_rows.outside_flank`
  is its one home now), `region_total_density` and the `_rate` helper only it used, the `splice_graph`
  import they needed; `gdna_strand.binomial_scale`; `structural_claims.interface_masks` (dead once the
  context dropped its masks; `structural_claims` is now an instrument-facing entry point).
* `sweep.py`: the `structural_claims` import and the interface-mask locals gone; the assertion
  `flux_rows_finite` renamed `lam_rows_finite`; the docstrings say the two-phase shape (the self-solve,
  two passes, one solve, one write-back, four assertions); the vocabulary firewall now also bans `level`
  and `lane` from the backbone's identifiers.
* Tests re-expressed on the slimmed surface: `test_zero_count_is_a_measurement.py` (the counting term's
  gates on `count_logvar`, and the two shipped consumers — `hop_price` at a zero count, `poisson_level`
  at a zero-count gene edge), `test_region_init.py` (the sources gated on `f_*`/`struct_lock`/`tau_lam`),
  `test_solvability_denominator.py` (the one-home gate on the predicate, agreeing with the policy's `tau > 0`
  on every published value), `test_region_geometry.py` (the three `region_total_density` gates deleted).
  `vertex_ceiling.py` rebuilds the five-field `RegionInit`.

### Stage 3 — vocabulary
* `rename_census.py --sense relay` (the token added to its list; the per-site dump is the census):
  `relay_pool_ab.py` → **`message_pool_ab.py`** (`git mv`; `benchmark_report.py`, `CLAUDE.md`, `ROADMAP.md`
  re-pointed; self-tests 11/11 and 10/10); `pass0_vs_oracle.py`'s solver class `relay_only` →
  **`message_only`** (its test and `SUCCESS.md` with it); `solvability_audit.py`'s `relay_delta` →
  **`message_delta`** (column `msg Δ`); the instruments' and `worst_objects.py`'s prose say "the messages".
  What remains of the word in `src/` is one-line history ("retired with the relay, 2026-09-09").
* `region_chain.py`'s "relay across a reference boundary" (a verb), `simplex_logodds.py`'s citation of a
  deleted test, `config.py`'s factory line — reworded.

### Stage 4 — the gates
* `tests/calibration/test_transfer_policy.py` (42 gates, 2,554 lines) is SEVEN files by message, each
  naming its subject: `test_transfer_policy.py` (the protocol, the silence identity, the passes against
  the recursive reference, the no-echo law, own claims, the completion contract),
  `test_transfer_splice_faces.py`, `test_transfer_edge_and_terminus.py`, `test_transfer_alt_splice.py`,
  `test_transfer_gdna_lane.py`, `test_transfer_rna_lanes.py` (the lanes and the cube),
  `test_transfer_ceilings.py`; the captured sweep and the shared builders are `_transfer_harness.py`
  (`capture_sweep_inputs`, wrapped by each file's module-scoped `sweep_inputs` fixture — a second
  `conftest.py` collides with the root one, which other tests import by name). All 42 pass.
* The relay mentions in test docstrings are rewritten as the shipped state (`test_ambig_scenario`,
  `test_toy_harness`, `test_antisense_intronic`'s xfail reason — "measured under the relay policy of the
  day", the numbers kept — `test_region_chain`, `test_pass0_vs_oracle`); the sandbox-doc names
  (`MESSAGE_PLAN.md step A/E`) and the rung/item numbers are gone from every gate.

### Stage 5 — the docs, by the MOVE RULE
* **`docs/dev/` 3,128 → 907 lines.** DELETED after their settled content was placed (each an ISSUES
  or TRAPS edit in the same change): `HONEST_PRECISION.md` → `ISSUES: the-message-policy-campaign`
  (six mechanisms refuted with numbers, the three derivations, the method lessons; `CLAUDE.md`
  re-pointed); `rename.md` (the owner's note) → `ISSUES: rename-row-and-face` (the owner's words kept,
  beside `rename-the-drain`); `BASELINE_2026-09-01.md` (a pre-ship snapshot under the relay, superseded
  by the ship-day page; git carries it); `MESSAGE_RUNGS.md` → `ISSUES: message-layer-open-cases` (the
  exon solve with every face speaking, the factory on mixed-bit regions, the chain of termini, the
  substrate `nest`/`div`/antisense-nascent) and its parked log's items, every one already an ISSUES
  entry or a ROADMAP rank; `COMPOSITION_TRANSFER_STAGE01.md` → the scan's three refused forms into
  `ISSUES: two-sided-exon-row`, the hybrid-arm lesson as form ④ of `TRAPS: an-ablation-that-never-ran`;
  `MESSAGE_PLAN.md` → step A's two refuted forms and its stage 0 into `ISSUES: the-abundance-discrepancy-map`,
  the 2026-09-05 reach refusal into `ISSUES: levels-always-travel-for-the-gdna-lane`, the popped-by-face
  lesson as form ⑤ of the same trap; `TWO_PHASE_BACKBONE.md` and `AMBIG_DESIGN.md` (their rulings and
  refusals were already `DESIGN.md` §6b.12/§6b.13 and ISSUES; the two citations of them in `DESIGN.md`
  and `CLAUDE.md` re-pointed). `PLAN_measured_prior.md` (priority 3's record) keeps its plan; its
  duplicated Part B is one paragraph and its re-scan incidents live in `TESTING.md` §2 (`RIGEL_SCRATCH`
  must be exported; a widened bank is `--expect-changed`). ⚠ **`message_notes.md` is the OWNER's own
  design notes and was left untouched** — deleting an owner-authored note is the owner's call.
* **`DESIGN.md`** converged on the shipped state where it described deleted mechanisms, each ruling's
  numbers kept: §0c.1 (the mechanism is the transfer policy's splice-in face map), §0c.2 (the mute of
  2026-08 as a record, with the ruling that a number measured with a known defect live is a record of
  the defect), §6 (the module list) and §6.1 (the layout table: `silent`, `transfer`, `transfer_rows`,
  `messages/__init__`; the two-phase interface; the belief-field paragraph), §6b retitled "the message
  layer and ψ's reference" with its first paragraph a record, §6b.2 stamped a record of the retired
  anchor factor (§6b.3, §6b.7, §6b.10 were already stamped), §6c MOVED after §6b.14 (it sat between
  §6b.12 and §6b.13), and the four stale citations (`terminus_flank_gain` ×2, `test_splice_flux_reframe`,
  the mass pin's module) re-pointed. ⚠ NOT done: the §6b.4–§6b.9 section titles still carry the plan's
  "rung N / item N" labels (the record's own names; other docs cite the sections by number), and the
  sections' bodies were not rewritten — an owner call whether a title should name the builder instead.
* `EQUATIONS.md`: the six stale code citations in §3.5–§3.6 and §9c re-pointed or marked as the relay's;
  the §3.5–§3.6 derivations themselves stay under the record stamp the retirement gave them (the owner
  call recorded in the previous handoff stands). `TRAPS.md`: the two stale citations fixed; two forms
  added to `an-ablation-that-never-ran` (above). `ISSUES.md`: three deleted-instrument citations
  re-pointed (`ladder_arm_ab`, `transfer_variance_audit`). `SUCCESS.md`: three "the relay" sentences
  now say the messages. `TESTING.md` §0b's toy family and `SUCCESS.md`'s run order name only existing
  instruments (checked against the disk).

## THE TWO ITEMS AFTER THE CLEANUP (the owner's direction, 2026-09-09, late)

1. **`ISSUES: flux-source-skipped-at-an-empty-exon-piece` — LANDED** (`transfer._rna_lanes`, `_LevelLane.emit`,
   `flux_witness`; the gate `test_an_empty_exon_piece_beside_a_lit_junction_is_a_flux_source` in
   `test_transfer_rna_lanes.py`, written first and watched failing, then three perturbations watched
   firing: empties-only-forward, the empty source stamped with its own zero count, the junction's counting
   alone). ⭐ The census first: on the ladder over half of all junction-adjacent exon pieces are EMPTY
   (2,100–2,200 per strand) but only ~70 per strand have any lane face (the other strand's boundary), and
   those reach the AMBIG exon|exon boundaries; the test chromosome has NONE, so the ladder judged it.
   ⛔ The issue's own example (the sj+terminus piece between TB's junction and TA's TES) cannot be
   reached at all: a strand's RNA level stops at that strand's own terminus by the face rule. A/B on all
   16 rows, both frames (`flux_source_ab.py`, session scratchpad): the landed form — the empty piece's
   level priced by `hop_price` on its zero count, emitted with the flux's witness — every in-scope row
   within 0.5 % through the pipeline (worst `g98 ss.99 ON` 1.0048×), the stranded zero controls 0.977× /
   0.958×, pass zero 6/8 both halves; the sharper price (the junction's counting alone at a no-opportunity
   piece) REFUSED at 1.0179× / 1.0091× on the stranded capture-ON rows
   (`ISSUES: the-empty-flux-source-at-the-junctions-counting-alone`). The landed source reproduces the
   judged prototype to the fragment (176,162 on `g98 ss.99 ON`, 15,365 on `g00 ss.99 ON`). No golden
   moved. ⭐ NEW IDENTITY REFERENCES frozen on the landed tree:
   `~/Downloads/rigel_runs/arms/fluxsrc_identity_{gdna_g05_ss_0.99_nrna_mid_capture_on,gdna_g05_ss_0.50_nrna_mid_capture_off}.json`
   — every future cleanup gates against THESE (the `retire_identity_*` pair is the pre-landing tree).
2. **`vertex_ceiling.py` RE-POINTED and its population corrected.** The pin now enters where the two-phase
   solve reads: the pinned node's OWN CLAIM (`transfer._claims`, patched) and a delta row at its truth
   delivered into ψ through the policy's `solve` (`_PinnedPolicy`, installed as `calibrate.TransferPolicy`);
   the classification runs in the `sweep.build_region_init` wrapper. `--self-test` 27/27, `noop`
   byte-identical to `base` on the ladder row measured. ⛔ FOUND ON THE WAY: the instrument's REALIZED-vertex
   population (truth exactly 0 or 1 per object) priced CHANCE — on `g50 ss.50 OFF` most of its 22,201
   evidence-free vertex slots were boundaries of 6–19 crossings that all happened to be gDNA; pinned as
   certain and propagated they drove RNA-rich exon neighbours from 0.46 to 0.64 (truth 0.01) and the region
   axis 23 % WORSE (Σ|err| 728 k → 894 k) while the pinned slots went to zero error (`vertex_pin_diag.py`,
   session scratchpad). The population is now the PARAMETER vertex (`_parameter_vertex`, gated with
   perturbations): every region of a silent gene no expressed gene overlaps, every intron of a gene with no
   nascent fragment, every counted slot of a zero-gDNA row; it needs `--oracle-cache`. On that row
   `vertex_free` then IMPROVES every column, by 0.4 % of the region axis (728,378 → 725,752).
   ⭐ **THE CEILING, RE-MEASURED ON THE WHOLE LADDER UNDER THE SHIPPED POLICY (2026-09-09; `base` vs
   `vertex_free`, 16 rows, ~1 h; `ROADMAP.md` rank 6's number).** `noop` byte-identical. Pins per sweep and
   the FINAL solve's Σ|err| ratio (vertex_free / base), region axis, boundary axis:

   | row | pins/sweep | region | boundary |
   |---|---|---|---|
   | `g00` (all four) | 8,900–16,000 (every counted object) | 0.000–0.004 | 0.000 |
   | `g05 ss.99 OFF / ON` | 11 / 18 | 1.000 / 0.999 | 1.000 / 0.999 |
   | `g50 ss.99 OFF / ON` | 44 / 48 | 1.000 / 0.999 | 1.000 / 0.999 |
   | `g98 ss.99 OFF / ON` | 104 / 236 | 0.999 / 0.991 | 0.998 / 0.988 |
   | `g05 ss.50 OFF` | 448 | 0.983 | 0.995 |
   | `g50 ss.50 OFF` | 1,104 | 0.978 | 0.986 |
   | `g98 ss.50 OFF` | 2,111 | 0.931 | 0.956 |
   | `g05 / g50 / g98 ss.50 ON` (deferred) | 445 / 993 / 2,631 | 0.648 / 0.742 / 0.765 | 0.883 / 0.920 / 0.910 |

   Reading: the vertex information is worth ≤ 1 % on every stranded in-scope row (the strand term already
   reaches the vertex), 2–7 % of the unstranded capture-OFF rows rising with gDNA (silent genes' regions
   and nascent-free introns, the landscape prior's training population — priority 3's territory, not a
   message's), and 25–35 % of the DEFERRED stratum; the zero rows are total by construction. Pass zero
   reads the same way (unstranded OFF 0.996–0.999 at pass zero, so the gain arrives through the refit
   prior). ⚠ Every number in the instrument's 2026-08-05 record is now triply historical (the retired
   36-row ladder, the relay, a population that priced chance).

## MEASURED ON THE CLEANED TREE (2026-09-09; both instruments, the ladder; re-derive, never quote)

`calibration_vs_oracle.py` (the 0.8.0 metric, mass-weighted |Δ gDNA share| per object, region / boundary
axis; the library gDNA fraction's |P − O| beside it): stranded × OFF **0.0067 / 0.0104** (0.003);
stranded × ON **0.0104 / 0.0149** (0.004); unstranded × OFF **0.0088 / 0.0124** (0.003); the `g00` zero
controls 0.0081 / 0.0077 (0.008); the DEFERRED unstranded × ON 0.0742 / 0.1166 (0.090). `policy_benchmark.py
--panel ladder --by-class` (689 s): transfer beats silence 7/8 on each half (worst 1.01× / 1.00×); the
misplaced gDNA as a share of every non-intergenic fragment is **0.6–2.4 % on every in-scope contaminated
row** (silent 0.6–3.2 %) and 1.9 % on `g00 ss.50 OFF` (150 k false fragments; silent 15.6 %). By class,
the share of the CLASS's own fragments misplaced: the intron 1–2.3 % off capture (43–45 % of the absolute
error because its mass is largest; 4–19 % on capture-ON rows where it holds almost no fragments);
exon|intron boundaries **4–8 %** off capture (15–17 fragments per slot: about one fragment each, the
counting floor; the class the messages moved most, 30 k → 20 k at `g50 ss.50 OFF`); every exon and
exon|exon class ≤ 2.5 % except the `g98` capture-OFF rows, **4–9 %** across all of them — near-pure gDNA
under-called, the vertex atom's signature (`ROADMAP.md` rank 6), not a message's.

## FOUND, NOT CHANGED — candidates for priority 2 (each would move a number, so each is a proposal)

* **An edge boundary's own strand profile is overwritten by the edge level's marker.** `_claims` gives a
  single-strand boundary with a live channel its strand profile; `_edge_level` then sets `own[b] = zeros`
  at every intergenic|exon edge as the level claim's marker, so the RNA lanes never read an edge
  boundary's strand profile as an RNA source (the gDNA lane is unaffected: a gene edge's level is its
  Poisson count). Preserved exactly; whether the edge's strand profile should be an RNA source is a
  measurement.
* **Two liveness predicates.** The policy's test of a node's strand channel is `tau_lam > 0`
  (`_Chain`), the instruments' one home is `has_own_composition_evidence` = `tau_lam > 1e-9`. They agree
  on every value observed (τ is exactly 0 or a Fisher information far above the guard) and the one-home
  gate now asserts that agreement on a vector; a ladder census of `0 < tau_lam ≤ 1e-9` would close it for good.
* **The `struct_lock` xfail pair's rationale is the relay's.** Nothing on the solve path reads
  `struct_lock` any more (the policy and the backbone read `f_*` and `tau_lam`); the mask reaches only
  the instruments that classify objects by it. The measured objection to scoping it to
  `g1_locked ∧ REGION` (the zero-control +3,207 %) was the relay's; the repair is now a change to an
  instrument-facing mask and an owner call (`tests/calibration/test_region_init.py`'s note says so).
* `PsiMessage.lam_rows`' docstring still narrates the retired relay's certified-flux stream as the
  channel's origin; harmless history, one paragraph.
* Nine rung/item labels remain in `DESIGN.md` outside §6b.4–§6b.9 (§6b.10–§6b.13's bodies) and in
  `ISSUES.md`'s older entries; the owner's decision covered §6b.4–§6b.9.

## WHAT IS NEXT — the owner's priorities 2–4 (unchanged from the previous handoff, re-pointed)

### 2. Improving the message policy (`ROADMAP.md` rank 2)
The debug loop on the worst IN-SCOPE rows of the published page (`g05 ss.50 OFF` at 1.01×, the by-class
tables' walled exons and exon|exon boundaries), one `policy_prototype.py --module` arm at a time, A/B'd
on the three panels and the ladder, halves apart, pass zero beside the pipeline (the empty-piece flux
source is landed, above): `ISSUES: flux-price-witness-units`,
`ISSUES: two-sided-exon-row` (its refused forms are now listed in the entry), `ISSUES: flux-floor-dispersion`,
`ISSUES: ambig-node-as-a-gdna-source`, then `ISSUES: message-layer-open-cases` (the substrate:
`nest`, `div`, the antisense's nascent variant). ⛔ Compare `src` against `src`: the prototypes under
`~/Downloads/rigel_runs/prototypes/2026-09-09_ship/` subclass the pre-cleanup `TransferPolicy` and its
`_Lane`/`_RnaLane` names, and will not import against the folded class without a one-line rename.

### 3. The gDNA landscape prior (`ROADMAP.md` rank 3)
`ISSUES: gdna-landscape-trains-on-false-positives` (the owner's 2026-09-06 ruling in its entry: the
enrichment witness IS the landscape prior; bound-only nodes do not train it). First step: census the
training population per node class on the zero rows, then the estimator at bound-only nodes;
`ISSUES: measured-prior-rung-4`, `ISSUES: landscape-trains-on-real-substrate`, `docs/dev/PLAN_measured_prior.md`.

### 4. The post-calibration, pre-EM setup (`ROADMAP.md` rank 4)
`priors.py` / `result.py` / `derive.py` against `prior_vs_oracle.py` and `calibration_vs_oracle.py`'s
ruler column; first step: re-run `prior_vs_oracle.py` on the committed tree so the assembler's own error
is re-recorded under the shipped policy before anything moves.

## THE LESSONS THIS SESSION PAID FOR

* **Two references, one per half, make a refactor provable one stage at a time.** Five stages, ten
  identity checks, zero bits moved — and the one place a fold could have moved a bit (`count_price` into
  `hop_price`) was proved safe by the gDNA lane's counts being positive, not by the suite.
* **The census counts `src/`; the readers are in `tests/` and `scripts/`.** Every deletion had a test or an
  instrument reading it (`ctx.order` in eight gates, `RegionInit` in `vertex_ceiling.py`,
  `region_total_density` in three gates), and the suite named each one.
* **A second `conftest.py` shadows the root one.** Other tests import `conftest` by module name, so a
  package-level fixture file broke ten unrelated modules at collection; a plain builder in a helper module,
  wrapped by a fixture per file, costs a few lines and nothing else.
* **Delete a thread record only after grepping its numbers in the permanent docs.** Three of the eight had
  refusals whose numbers lived nowhere else; the MOVE RULE is what turned those into ISSUES entries instead
  of losses.
