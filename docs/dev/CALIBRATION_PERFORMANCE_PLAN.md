# Calibration performance — where we are and the plan (2026-09-11)

A working document, provisional like everything in this directory. The permanent homes are `CLAUDE.md`,
`docs/ROADMAP.md` (rank 1 carries the order without numbers), `docs/ISSUES.md`
(`performance-memory-bounded-solve` carries the numbers), `docs/DESIGN.md` §6b.15 (the rulings) and
`docs/TESTING.md` (how to run each gate). This file holds what none of them should: the whole picture in
one place, the reasoning behind each next step, and the implementation notes a session needs on day one.

**Two owner rulings frame everything below.** Cleanup comes first: the code must be simple and elegant
despite what it does, and a port lands only on a Python design and implementation we are satisfied
with. And names must be intuitive: no new word for a thing that already has one, and none the owner has
to look up (see §5).

## 1. Where we are — commit `38755cb1`

Calibration is the tool's one unfinished component. On the 18.6M-fragment library (`mctp_vcap_rna20m_dna05m`,
8 threads, every number a back-to-back pair with untouched stages at 1.00):

| tree | wall | peak RSS | the four sweeps |
|---|---|---|---|
| `main` before the session (`710ed209`) | 892–917 s | 32.6–33.2 GB | 747–769 s |
| + the locus blocks | 875–892 s (0.97) | **14.9–15.0 GB** | 727–738 s |
| + the shared message layer | **515 s** (0.65 of its off-arm) | 17.8 GB (2.68 GB shared) | ~390 s |
| + the face tables | a no-op refactor, not re-measured at scale | | |

Net: wall ≈ 900 → 515 s (1.75×), `calibrate` ≈ 800 → 424 s (1.9×), peak 33 → 17.8 GB (−46 %), the
per-sweep working set 9.87 → 0.19 GB. Where the 515 s goes: scan 31, second pass 22, `calibrate` 424 —
sweep 0 with its full message layer 136 s, refit 1 with its layer 176 s, refits 2–3 served from the
shared layer 38 s each, the init ψ 27, the landscape fits 7 — quant 23. What is left inside calibration
is two message layers (~250 s: prepare ~95, the passes ~145, the policy's solve ~15) and four ψ solves
(~120 s); the hundredfold on calibration is theirs and needs the port. Only within-pair ratios are quoted:
the same tree drifted 875 → 795 s between sittings. Reports: `~/Downloads/rigel_runs/perf/ab_locus_2026-09-11/`
(main vs the blocks) and `ab_memo_2026-09-11/` (shared layer off vs on).

What landed, each gated before the next (`DESIGN.md` §6b.15 has the rulings and every measured fact):

1. **A terminal receives nothing.** A region admitting no RNA strand (`is_region & g1_locked`) is
   structurally pure gDNA, solved before any message exists; the backbone never asks a kernel for a hop
   into one (`sweep._pass`), and no lane lists a face from one. The chain breaks at every terminal into
   33,018 loci on the human chain (median 19 slots, the largest 0.12 %). `region_chain.locus_blocks`
   merges loci to `CalibrationConfig.sweep_block_slots = 5000` and `sweep.solve_chain` solves each block on
   its own slice of every input, reading one terminal beyond itself (`_solve_block`; `blocks.block_slice`, `blocks.gather`).
2. **The library is the only cross-block information.** `Policy.library(ChainView)` runs once over the
   whole chain on a view with no belief field; `prepare(ctx, library)` sees one block. The factory rows
   travel on the context (`factory_rows`); the context carries the own-composition bits
   (`has_own_composition`), not the self-solve object.
3. **ψ's read-out is chunk-exact.** The one step that moved a number (≤ 3.1e-15 per slot, no
   amplification, TPM bit-identical, the oracle metric at the last ulp; owner: identical to a tolerance).
   Every block size gives the same bits, so the block size is a working-set knob.
4. **The refit sweeps share their message layer** (`message_cache.MessageCache`). Everything the layer reads is on
   the context, the library and the grid — never the prior — and the belief is reset before every sweep,
   so sweeps 1–3 deliver identical messages; the content-keyed cache serves refits 2–3 their rows
   (38 s instead of 176 s each) for 2.68 GB held.
5. **The rules are typed tables** (`messages.faces.Faces`): every directed face is one of a node's
   two sides, so a rule is a KIND (FORWARD / TRANSPORT / SPLICE_OUT / EDGE / LEVEL) and its parameters
   at ``(destination, side)``, with a row store for the ``(K,)`` maps; `Faces.apply` is the one home of
   the rule arithmetic; a face carries one rule (a second write is refused — the documented precedence
   had no instance). The lanes hold their faces as ``(n, 2)`` bits and their junction flux as a row table.

Threads were measured and refuted for the Python passes (0.83–0.94× at 8 threads; forked processes
6.16×); the owner deferred parallelism to the port.

## 2. The rules of the game (unchanged, and load-bearing)

* **Every step is proven before it is believed.** A pure restructure is bit-identical: `sweep_replay.py
  replay --dir ~/Downloads/rigel_runs/perf/sweeps_MO_3021_step5 --call 0..3` (and `--block-slots N|none`,
  `--tolerance` for the report of what moved) and `rename_identity.py --check --reference
  ~/Downloads/rigel_runs/arms/port_identity_*.json` (two ladder conditions and `--bam` LBX0190). A change
  that may move numbers is judged on the oracle metric per stratum, the panel, and profiler pairs
  (owner, 2026-09-12: minuscule changes are not a concern; elegance is the bar). The suite (`CLAUDE.md`'s baseline line), `ruff`, `preflight.py --full`.
  A timing is read only from back-to-back pairs (`profiler.py --compare`); patch `calibrate` for an arm
  through `importlib.import_module("rigel.calibration.calibrate")` — `import … as` binds the re-exported
  function and the arm never runs (it cost one pair this session).
* **A falsification test first, verified failing; then break the fixed code and watch each gate fire.**
  This session's gates each had a perturbation cycle; two cycles that did not fire exposed a design fact
  (a dead precedence) rather than a gate weakness — read a silent cycle that way first.
* **One mechanism at a time; no magic numbers** (a knob is derived or measured and its ladder recorded in
  its docstring, as `sweep_block_slots`'s is).
* **The owner drives commits.**

## 3. The ordered plan

Each step: what, why, the design, the gate. Steps A–E are Python and precede the port by ruling.

### A. Names and cleanup (first, and the owner's paramount priority)

Bit-identical throughout; the gates are the replay, the references and the suite.

1. **The naming pass** — DONE 2026-09-11 (§5 records the ruling and the names as applied).
2. DONE 2026-09-12. **Split `sweep.py`** (~1,000 lines, three concerns) into the backbone (`_pass`, `_check_message`,
   `solve_chain`, `_solve_block`), the block plumbing (`blocks.py`: `_block_slice`, `_view_fields`,
   `_gather`) and the cache (`message_cache.py`); `messages/transfer.py` (1,200 lines): the lanes
   (`LevelLane`, `gdna_lane`, `rna_lanes`) into `messages/lanes.py`, and the face table with its three
   helpers (`side_of`, `norm`, `fuse`) into `messages/faces.py`, because the gDNA lane reads the table's
   kind while the policy imports the lanes. Layer 6 throughout; `_layers.py` lists them; `test_layering.py`
   and `module_census.py` confirm. Sizes: `sweep` 777, `blocks` 120, `message_cache` 126, `transfer` 680,
   `faces` 237, `lanes` 346.
3. DONE 2026-09-12 (`MessageCache.key(ctx, library, policy)`; a skipped field fires the perturbation gate).
   **The cache key iterates `BlockContext`'s dataclass fields** rather than a hand-built dict, so a new field
   cannot be left out of the digest (correctness, not tidiness; the six perturbation channels in
   `test_sweep_backbone` gate it).
4. **FOLDED INTO D (owner, 2026-09-12).** The single-strand tiling in `_solve_regions_logodds_all` is
   redundant only at the per-block callers; `init_beliefs` calls the same dispatcher on the WHOLE chain
   (the pre-sweep solve D names as the run's high-water mark), where dropping it would materialise the
   regridded ``(n_ss, 513)`` priors for the chain at several GB apiece, and the tiling's own comment calls
   it a cache-locality knob too. A row bound the block caller passes would be a knob. When D makes the init
   solve per block, the tiling goes for every caller at once.
5. DONE 2026-09-12: `build_region_init`'s dead `chain` parameter; `n_slot` / `spliced_slot` are properties
   of `ChainView`, which carries `spliced_count` per strand as it carries `unspliced_count`. Still open,
   lower priority: the 25-key `_capture` dict as a typed record — DONE 2026-09-13 as W7 (`blocks.SweepCapture`).

### B. The received messages as tables — DONE 2026-09-12 (bit-identical on every gate)

As built: `Received` / `Levels` in `messages/__init__.py` with one column beyond the sketch below —
`Levels.has_witness`, where today's RNA witness was `None` on the gDNA lane and a dead strand channel; the
kernel is `prepared.propagate(received, backward=…)` returning `receive(source, destination)`, which
writes the destination's row; `sweep._pass(seq, nbr, prepared, n_grid, …)` allocates and returns the
table; a lane's `emit(s, x, levels) -> bool` writes row `x`, `receive(levels, s, x)` re-prices it in
place; `blocks.gather` concatenates the tables (`Received.take` / `concat`). The design as agreed:

**What a node holds after a pass.** Today, after each directional pass, every node holds one Python object
from the neighbour on that side — `Message(composition, level_gdna, level_rna_pos, level_rna_neg)`, each
level a `Level(profile, n, a, rna_count, rna_count_var)` — or `SILENCE` (the neighbour spoke and had nothing
to say) or `NO_NEIGHBOUR` (an open end); the backbone keeps two lists (`from_left`, `from_right`), the policy
keeps its own copy of the far-side state during a pass (`_PreparedTransfer.held`), and `solve`, the informed
predicate, the diagnostics capture and the census read the objects node by node.

**The design.** One table type replaces the objects, in words the interface already uses:

```
Received                      # what every node of a block received from ONE side, after one pass
  has_neighbour  (n,)  bool   # the side exists — False is today's NO_NEIGHBOUR
  composition    (n,K) f64    # the composition profile that arrived, where has_composition
  has_composition (n,) bool
  level_gdna, level_rna_pos, level_rna_neg : Levels
Levels                        # one population's level lane, received
  present        (n,)  bool
  profile        (n,K) f64
  count, opportunity, rna_count, rna_count_var  (n,) f64   # the Level's n, a, rna_count, rna_count_var
```

Silence needs no word: a node with a neighbour and nothing present received SILENCE. `from_left` and
`from_right` are two `Received` tables per block. The kernel `receive(source, destination)` writes row
``destination`` of the table for the pass it runs in (the backbone still owns the table and the order: it
hands the kernel the table and calls it in chain order, so "the backbone writes held" keeps its meaning);
the far-side state the policy consulted from its own copy is the same table's row ``source`` — one copy,
not two. `solve` reads the two tables and becomes array code over nodes (`_fuse` of the compositions,
the gDNA level rows through `profile_of_level`, the ceilings, the cube rows). The informed predicate reads
`has_composition` from both tables. `_check_message` is unchanged (it checks the delivered `PsiMessage`).
The capture stores the two tables; `landscape_training_census.held_evidence` reads them. The message
cache stores what it stores today — the delivered `PsiMessage` and the composition bits.

**What does not change:** the rulings of `DESIGN.md` §6b.11–12 — two phases, the recipient decides, every
node ends the passes with one thing from each neighbour it has, silence is received and an open side is
not — and the four backbone assertions. `Message` and `Level` retire; `SILENCE` and `NO_NEIGHBOUR` survive
as the two states the table expresses.

**Order of work.** (i) `Received`/`Levels` in `messages/__init__.py` with their docstrings; (ii) the
backbone's `_pass` creates the table and drives the kernel; (iii) the transfer policy's `receive` writes
rows and reads its far side from the table; `solve`, `_ceilings`, `_cube_rows` read tables (vectorise where
it stays legible); (iv) the informed predicate, the capture, `_gather`, the census; (v) the ~30 test sites
that build or inspect `Message` objects (`test_sweep_backbone`'s echo policies, the harness's
`_drive_the_backbone`, the RNA-lane and policy tests) re-expressed against rows; (vi) the docs' wording.
**Gate:** the replay bit-identical on all four captures and `--block-slots 5000`, the three references,
the suite, the block-invariance and cache gates; break-the-code cycles on `has_neighbour` (an open side
read as silence) and `has_composition` (a level read as a composition).

### THE PRE-PORT WORKLIST (owner, 2026-09-12) — the code becomes pristine in Python before any port

Ruled 2026-09-12: bit-identity is no longer the bar for the remaining Python work — the owner is not
concerned with minuscule changes; the bar is elegant, simple, efficient, clear, concise, maintainable code,
judged on the oracle metric (`calibration_vs_oracle.py`, per stratum, both zero controls), the panel
(`policy_benchmark.py --panel test|ladder`), the suite, and timing on back-to-back pairs. Every item below
precedes F (the port). Status is kept HERE; tick an item by writing DONE and the date beside it.

| # | item | what | judged by |
|---|---|---|---|
| W1 | **C, shrunk to the report** — DONE 2026-09-12 (`EQUATIONS.md` §9d holds the derivation; the report reads zero moved on the current tree) | `sweep_replay.py replay --tolerance` prints per output array the slots moved, max abs Δ, max relative Δ, beside the bit verdict; the derived budget (½·ε·T·c_κ·N on fractions, L̃² of it on log-variances; §C) printed as a sanity bound; no frozen-array companions, no verdict on the transcript table | the report on the current tree reads zero moved; `--self-test` |
| W2 | **One ψ solver, one precision** — DONE 2026-09-12 (`_solve_logodds`; the unread strand log-variances deleted with it; metric and panels identical to the printed precision; timing pair in DESIGN §6b.15) | a single-strand slot is the cube with a tilt grid of one cell, so `_solve_regions_logodds` and `_solve_ambig_logodds` become one solver with `K_t` a parameter, in float64 throughout (the float32 cube was a memory choice the tiling made moot); `f32-strand-tilt-at-half` closes with it; the cache stores what the solve produces | the oracle metric and the panel, W1's report on the replay captures, timing pairs, the vertex-reference and strand-reference gates |
| W3 | **`calibrate.calibrate` as named stages** — `calibrate` DONE 2026-09-12 (141 lines of orchestration over `_fit_strand`, `_IntronFactory`, `_background_pair`, `_abundance_landscape`, `_policy`, `_Solve` + `_init_belief` / `_sweep` / `_solve`, `_result`, `_log_summary`; bit-identical on the three references); `_solve_block` (109 lines over `_psi`, `_composition_arms`, `_message_layer`, `_write_back`, `_block_diagnostics`) and `solve_chain` (136, over the `_Structure` and `_Sweep` records) DONE 2026-09-12, every step bit-identical on the three references, the replay and the suite | the 600-line function becomes the walk's rungs (init → strand → local → messages → refits → shipped) as functions with one job each; `_solve_block` (319) and `solve_chain` (158) likewise | `rename_identity.py --check` where a step is a pure restructure; the metric otherwise |
| W4 | **Memory (D)** — DONE 2026-09-12: measured first; the crossing divisor in closed form (the 9 GB sj matrix gone), the landscape's kernels tiled (the 3 GB peak gone), the factory rows per block (`FactoryRows`); peak 19.2 → 11.4 GB, wall 0.99, metric untouched (3.4e-15); the ψ tiling stays, since the whole-chain `init_beliefs` solve still runs on it and is tiled to 1 MiB anyway — the plateau's remaining owners are the pre-calibration floor and the cache (E) | `build_region_geometry`'s 14.5 GB transient; `init_beliefs` per block, which frees the ψ tiling for every caller (A.4); the factory rows built per block | `profiler.py` peak and held per stage, pairs |
| W5 | **The grid study** — DONE 2026-09-13 (`DESIGN.md` §6b.15: one λ lattice, `sweep_logodds_step` 0.2, `sweep_n_tilt` 60 explicit; the second grid, `_regrid_global` and `_scaled_grid` deleted; refusals in `ISSUES: the-second-lambda-grid-and-its-regrid`; the θ quadrature filed) | understand ψ's grid before designing it: accuracy vs K per stratum and slot class, the read-out's quantisation, the bracket coupling, the regrid's cost, time and memory vs K; then the simple design — §6 | the oracle metric per stratum with both zero controls, the panel, profiler pairs |
| W6 | **The tunables census** — DONE 2026-09-13: `CalibrationConfig` had 11 fields, not 49, every one read; four were unearned and are gone — `message_propagation` (folded into `message_policy`, the one read site was `not propagation or policy == "silent"`), the `intron_factory` switch (an off branch no instrument or test set), `gdna_prior_strength` with its `strength` plumbing (a temperature nothing moved; exact Bayes only), the `abundance_landscape` switch (a QC-only fit nothing turned off); the oracle instrument's two bespoke flags folded into `--set`; the dead pre-sweep AMBIG cube in `init_beliefs` removed on the way; `background_abundance` kept as an unruled decision (`ISSUES: background-abundance-pair-unruled`); every kept knob carries its derivation or its measured ladder; every step byte-identical on the default rows of both substrates | `module_census.py`, the suite |
| W7 | **The capture as a typed record** — DONE 2026-09-13: `blocks.SweepCapture`, 25 typed fields (the 38 keys' 13 unread ones gone, and with the message-free variance the extra ψ solve that fed it); the request an empty record the sweep fills; `SweepCapture.gather` one typed concatenation per field kind; eight instruments and six test files read fields; the parity instrument compares fields | every capture reader byte-identical before/after on one condition per substrate; the suite; the default rows' identity |
| W8 | **Vocabulary rulings** — RULED 2026-09-13 (owner): no rename; `drain`, `row` and `face` stay (the two entries closed as rulings; §7 keeps the census). The `hygiene-ledger` items stay open as a later item of their own | — |
| W9 | **The two xfails** — RULED 2026-09-13 (owner): both DEFERRED, neither worth the pre-port thread's effort. The antisense prior-assembly casualty belongs to the prior-assembly session, which plans to change the very rule the xfail records (synthetic nascent at Dirichlet alpha = 0; today the leak is 72 against the test's 50) — `ISSUES: antisense-prior-assembly-casualty`. The two-sided exon row belongs to the calibration-accuracy thread: the toy reads a factor of 88 today (0.762 beside a pure-gDNA intron, 0.0086 beside a nascent-bearing one), but on the ladder its class is 4–6 % of the unstranded error with the transfer policy at parity there, while the introns carry 38–46 % — `ISSUES: two-sided-exon-row`, whose first step is the enrichment witness's derivation. Both xfails stay as executable records | — |
| W10 | **The hygiene ledger's pure cleanups** (§8) — DONE 2026-09-13, nine commits `0f3ca3e3`…`bdc1089c`: the comment; the rename EXTENDED to all five deconvolved arrays (`count_<population>_<axis>`); 33 moment tests restored; four bank-readers RETIRED after an instrument census (owner's ruling) and two migrated, `fl_pool_purity` now pricing the model production builds (true gDNA mean exactly at capture-OFF). The index alias map stays in the ledger | fresh identity references (`hygiene_identity_*`, bit-identical throughout); the suite 3,452 / 2 xfail / 3,454; each migrated instrument's recorded before/after in its commit |
| W11 | **The coverage census** (§8) — MEASURED, RULED and LANDED 2026-09-13 (owner: all tiers, converge on the production surface; thirteen commits `97947ede`…`d4934c5c`, the record in `ISSUES: hygiene-ledger`; suite 3,419 / 2 xfail / 3,421; identity references and oracle rows identical throughout). The measurement: the suite (3,452 tests) and all ten `--self-test`s under `coverage` (data in the session scratchpad, not the tree): 13,487 statements in `src/rigel`, 2,191 never executed (84 %), 28 files fully covered; 63 whole functions (1,567 statements, 72 % of the misses) — `sim/locus_sweep.py` entirely (466, one thin wrapper, no manual names it), `sim/net_flow.py` 80 % (no entry point but its test), `sim/suite.py main` + `simulate_suite.py`, the CLI command bodies, the simulator's sharded writers; the calibration package ~107 lines, nearly all guards and degenerate-input branches, plus a dead local helper (`simplex_logodds._sp`), four dead one-line properties, the never-fed `CompositionPriors.rna` socket and the never-set `boundary_rna_reach` taper arm. The owner's decisions are owed on the tiers (the session handoff lists them) | the census itself; anything removed is byte-identical on the default rows |
| W12 | **The θ quadrature** — DONE 2026-09-13 (step 1: `DESIGN.md` §6b.15, `EQUATIONS.md` §9e, `ISSUES: theta-quadrature-at-zero-gdna` CLOSED with its numbers; the nodes follow the strand term's peak, `K_t` = 24 derived; the mechanism corrected — one deep interior-tilt slot, a comb, never strand purity; `ISSUES: strand-marginal-volume-factor` opened). Step 2 DONE the same day: the lanes deliver a row's ingredients (`simplex_logodds.CubeRow`) and ψ evaluates them at its nodes; `sweep_n_tilt`, `_tilt_grid` and the row interpolation deleted — no tilt count exists | the one open design item inside ψ, settled in Python before ψ is ported | the marginal against adaptive quadrature, the oracle metric per stratum with both zero controls, the shared-exon deep stress, the AMBIG both-strand census |
| W13 | **The port's prerequisites** (§8) — DONE 2026-09-13 on `23a431a9`: `sweeps_MO_3021_step5` (four calls, each replays BIT-IDENTICAL: 30.0 / 12.7 / 12.7 / 12.7 s) and `port_identity_{gdna_g05_ss_0.50_nrna_mid_capture_off, gdna_g05_ss_0.99_nrna_mid_capture_on, LBX0190}.json`, each `--check`ed on the tree that froze it | re-capture the deep library's sweeps on the current tree (`sweep_replay.py capture`), since `sweeps_MO_3021_step4` predates the lattice; fresh identity references for the port thread | the captures replay bit-identical on the tree that made them; `rename_identity.py --check` |


### THE LANES WORKLIST (owner, 2026-09-13) — precedes the port; each item its own commit, stressed on the encompassing locus

The encompassing-transcript audit (a single-exon TB− over a two-exon TA+) found the message layer switched
off, and RNA levels not passed, by couplings a whole chromosome never shows. The owner's ruling: whenever a
strand is active and live (nonzero fragments) an RNA level is passed, alongside composition; real data will
stress every one of these; everything is fixed, stressed and proven before the port.

| # | item | what | judged by |
|---|---|---|---|
| L1 | **The lanes exist whenever their own coordinate does** — DONE 2026-09-13 (local commit) | `prepare` builds the layer without factory rows; the RNA lanes without the gDNA lane; the rung-0 identity gate retired | three gates; byte-identical on both panels |
| L2 | **One RNA coordinate, so a junction's flux is always a source** — DONE 2026-09-13 (local commit) | a level is absolute and the coordinate only an origin; the library's RNA coordinate is the pooled single-strand exon density of both strands, falling back to all exons; `CubeRow` carries one `rho_ref`; a strand with no single-strand exon still builds its flux levels | a gate on the audit's context (a + junction at a both-stranded exon builds and delivers a + level); the panels |
| L3 | **The strand channel's gate is a protocol decision** — DONE 2026-09-14 (`DESIGN.md` §6b.15, `EQUATIONS.md` §5.2b, the issue CLOSED with its numbers; the one-line form had been REFUSED 2026-09-13, `g00 ss.50 OFF` 499 → 21,484) | the channel is live iff the Bayes factor on the spliced 2×2 reads the protocol as strand-preserving (a free κ under the fit's own prior against κ = ½ exactly; closed form, no constant); `disc = 4(κ̂−½)²` where live; gDNA enters nowhere (`n_gdna_obs` deleted throughout); a gDNA-free stranded library keeps its strand channel and its strand-derived levels; the belief-read RNA level REFUSED as the relay | the gates (the structural invariant, the ladder's own scalars, the Occam scaling, four perturbations); both panels with both zero controls; the encompassing locus on a gDNA-free donor; the deep stress on the `g00` donor |
| L4 | **The encompassing locus in the toy ladder and the suite** — DONE 2026-09-13 (local commit; the exon∩exon solve gate xfails until L5) | TA+ (1,000–2,000, 10,000–11,000) inside TB− (0–20,000), four regimes: TA+ ≫ TB−, TA+ ≪ TB−, TA+ ≈ TB−, both ≈ 0 (all gDNA); gates: the − level delivered at every AMBIG slot when TB− is expressed, the + flux level at TA+'s exons when TA+ is, the exon∩exon solves within tolerance, all-gDNA reads gDNA | the suite (a stranded donor fixture) |
| L5 | **The witnessed atom** | as approved: tests, `src/`, docs; judged on the metric per stratum, both zero controls, the census bands, the shared-exon and encompassing stresses | — |

Then the deep-library baseline, fresh references, the port.

Not on the list, and why: the arcsine coordinate (REFUSED with its numbers; logit is finer at the vertex);
the vertex atom (a prior-family change whose whole ceiling is ≤ 1 % on stranded in-scope rows — see the
plan's owner note of 2026-09-12 — parked with that number).

### C. The tolerance gate — shrunk to W1 (owner, 2026-09-12); the derivation stands

A language port cannot be bit-identical (libm and summation order), so the replay needs a second verdict
beside BIT-IDENTICAL: `sweep_replay.py replay --tolerance`, reporting per output array the slots moved and
the max |Δ| and max relative Δ (the scratch comparator of this session did exactly this: per sweep
≤ 3.1e-15, no amplification). The budget must be DERIVED and written down — the empirical scale is Step 2's
≤ 1e-14 through four sweeps; argue it from the read-out's conditioning (a max-normalised log-profile's
exp/log through K = 138 cells) rather than choose it — and the same verdict is added to
`rename_identity.py --check` (digests give no tolerance; it needs the arrays or a fresh dump).

### D. Memory, now that the sweep is not the peak

The run's high-water mark is `build_region_geometry`'s transient (14.5 GB) and the pre-sweep `init_beliefs`
solve, not the sweeps (11.5 GB while they run). Measure both with `profiler.py` (peak and held per stage) and
fix what is found; then the factory rows built per block from the substrate slice (the last genome-wide
``(n_slots, K)`` arrays, 1–2 GB per grid, and most of the cache key's hashing). Expected peak near 8–10 GB.

### E. The shared layer's on/off switch — the owner's ruling

The cache trades ~280 s for 2.68 GB held. If a memory-constrained run should be able to turn it off, that
is a `CalibrationConfig` tunable with no effect on results; it was not added unasked.

### F. The port (③) — only when the PRE-PORT WORKLIST is done and the design is agreed

The unit is `sweep._solve_block`: one block's own claims, the policy's claims and rules (`Faces`, the lanes'
tables, the row stores), the two passes writing `Received` tables, the solve, the write-back. Order inside
the port: (i) the passes and `transfer_rows` (the highest ratio of Python overhead to arithmetic), (ii)
`prepare`'s builders, (iii) ψ (SIMD exp/log, the AMBIG cube), (iv) threads over blocks — the parallelism
the owner deferred to the port. Every step behind the tolerance gate of C, the references, the suite. Keep
the λ lattice a parameter (`sweep_logodds_step`, its point count following the bracket through
`calibrate.lattice_points`: an accuracy ruling, not a performance one). "Satisfied with the Python" means: A–D landed, one representation for the received
messages, the layering clean, the gates and the tolerance instrument in place, the docs current.

### G. Then the scan and the second pass (⑤)

31 s and 22 s at 8 threads on the deep library; the stages that scale with depth and the whole problem at
100M+ fragments. The scan's thread split is the owner's decision (`ISSUES: scan-thread-split-starves-the-workers`).

## 6. W5 — the grid study (the design as written before the study; RULED 2026-09-13, `DESIGN.md` §6b.15)

**What is there.** ψ's λ grid is ``λ ∈ [−L, L]`` with ``K`` points (`_logodds_grid`); ``f_g = σ(λ)`` is
read out as the posterior MEDIAN by a continuous quantile over the grid's histogram, interpolated on λ
(`_posterior_median_fg`, transform-invariant). Two grids: the coarse ``sweep_n_grid = 60`` — the AMBIG
cube's λ axis, the tilt axis ``K_t`` (``sweep_n_tilt``, default ``= n_grid``), every message row
(`lam_rows`, `cube_rows`), the factory rows, the composition arms — and the fine
``sweep_n_grid_single_strand = 256`` for the single-strand read-out, the arms and rows regridded onto it
per tile by `_regrid_global` (linear interpolation in ``f``). Both scale with the bracket: when the
landscape prior's support exceeds ``L``, `_scaled_grid` widens ``L`` and grows ``K`` to hold ``dλ``
fixed (138 / 513 on the deep library at ``L ≈ 23``). The fine grid's recorded reason (the config
docstring) is de-quantising ``f_g`` at high-mass single-strand slots; its recorded cost reason is the
cube's memory, which the tiling has since bounded to 1 MiB tiles.

**What to measure, in this order — each a config value, nothing in `src/` until a ruling.**

1. *The instrument.* `calibration_vs_oracle.py --set SECTION.FIELD=VALUE` (as `profiler.py` has), so a
   grid arm is `--set calibration.sweep_n_grid=K --set calibration.sweep_n_grid_single_strand=K_ss`.
   `policy_benchmark.py` likewise if the panel is wanted per arm.
2. *Accuracy vs K, one grid.* ``K = K_ss`` over a ladder of values (say 30 · 60 · 100 · 138 · 200 ·
   300 · 513), the metric per stratum with both zero controls; and today's pair (60 / 256) as the
   reference row. Read where the curve plateaus per stratum, and whether the deferred stratum differs.
3. *Where the quantisation sits.* Per slot class (`--by-class` thinking): the single-strand slots by
   mass decile — the median read-out's error against K should fall as the posterior narrows relative
   to ``dλ``; find the mass above which today's coarse grid is too coarse, if such a mass exists on the
   panel and on the deep library's slot distribution.
4. *The bracket coupling.* At a fixed ``K`` the bracket widening changes ``dλ`` unless ``K`` scales; hold
   ``dλ`` fixed as `_scaled_grid` does and ask whether the answer depends on ``dλ`` or on ``K``.
5. *The regrid.* Today's two-grid path against a one-grid path at ``K = 513`` for single-strand slots
   alone: the interpolation's own error, isolated.
6. *Cost.* Time and memory vs K on a profiler pair: the message rows, the cache (≈ 4 GB at K = 138)
   and the cube all scale with K; ``K_t`` separately.

**Candidate elegant designs, to be chosen by the measurements, not before:** one grid at a K derived
from the read-out's own resolution requirement (a rule, not a number); a read-out whose accuracy is
insensitive to K, making the grid a pure cost knob; ``K_t`` decoupled from ``K`` if the tilt's resolution
requirement is different. A two-grid design is the answer only if the study shows a single grid cannot
serve both at acceptable cost — and then it is a ruling with its numbers, not an inheritance.

## 7. W8 — the vocabulary rulings and the hygiene ledger (prepared 2026-09-13; RULED the same day: no rename, the terms stay; the ledger stays open)

**What W8 is.** Two issue entries hold three word rulings that are the owner's (`ISSUES: rename-the-drain`,
`ISSUES: rename-row-and-face`), and `ISSUES: hygiene-ledger` holds five small items, each its own commit and
none moving the metric. A rename is a pure restructure and is proven bit-identical; the hygiene items are
judged by the suite and the instruments they touch.

**Two facts found while preparing it.** `rename_census.py --sense` knows nine tokens (acceptor, cut, donor,
edge, junction, line, node, relay, seam) and none of W8's three, so the census must be taught `drain`, `row`
and `face` before any site is ruled on. And every identity reference on disk (`memory_identity_*`,
`locus_identity_*`, …) predates the one-lattice landing, so W8 freezes its own first.

**The opening steps, mechanical, in order.**

0. Freeze fresh identity references from the current tree — `rename_identity.py --freeze` on the two ladder
   conditions the earlier references used (`gdna_g05_ss_0.50_nrna_mid_capture_off`,
   `gdna_g05_ss_0.99_nrna_mid_capture_on`) and `--bam` LBX0190 — named for this thread (`vocab_identity_*`).
   Every rename stage below ends with `rename_identity.py --check --stage <name>` against them.
1. Teach the census the three tokens: `TOKENS` in `scripts/design/rename_census.py` gains `"drain": None`,
   `"row": None`, `"face": None` (each carries two senses, below), and `EXEMPT` gains the non-vocabulary
   senses the dump shows (`interface`, `surface`; the instruments' table rows; numpy rows). Gate: the dump for
   each token classifies every site, and `test_scripts_index.py` still passes.

**The three rulings — the owner picks the words; the census counts and the candidates, so the ruling is a
sentence at the start of the session.** Counts are words containing the token (src / scripts / tests / docs).

* `drain` — 141 / 146 / 298 / 50. What it names: the second pass's operation on the deferred fragments (pass
  one buffers a fragment whose mate gap admits more than one explanation; the second pass decides it), and by
  extension the FRAME after it (`DRAINED` / `undrained`, the frame the truth is certified in). Identifiers:
  `drain`, `drained`, `DrainQC`, `undrained`, `drain_seed`, `_drain_side_buffer`, `with_drain`,
  `lift_drain_parts`. Candidates: **`settle`** (the deferred fragments are settled; the settled frame; unused
  anywhere in the tree — the recommendation), `decide` (the chooser is `choose_hypotheses`), `place`;
  `resolve` collides with `resolve.cpp`'s fragment construction and `assign` with the EM's assignment.
* `row` — 1,236 / 752 / 915 / 335. Two senses entangled by the word's generality: a slot's max-normalised
  log-profile over the solve grid (`lam_rows`, `cube_rows`, `blur_row`, `transport_row`, `level_row`,
  `Faces.rows`, `transfer_rows.py`) against a numpy or table row everywhere else (`_row_moment`,
  `_block_rows`, the instruments' printed rows). Candidate: **`profile`**, the word the transfer policy's
  docstrings already use for the thing — but `Levels.profile` already names a level lane's profile, so the
  ruling must say whether a delivered row and a lane's profile are one word. The largest of the three and the
  one the census must classify site by site first.
* `face` — 303 / 26 / 237 / 125. One directed side of a boundary, the `(destination, side)` pair a rule is
  keyed by (`Faces`, `FaceRule`, `face_is_licensed`, `face_map_lambda`, `_splice_faces`, `M_face`); the
  collisions are `interface` and `surface` in prose, not the concept. Recommendation: **keep `face`** and close
  the entry with the ruling — it is a geometric word for a geometric thing, and SIDE is already taken by
  left / right.

**The hygiene ledger, in order of risk, each its own commit.**

* the stale comment in `pipeline.py` above the second pass's `build_fl_models` (the text "The SAME
  de-tilted RNA pool the calibrator will read" — the second pass's fl models are pass one's; the lift's
  docstring above it has it right). A three-line edit; the suite.
* `mass_*_boundary` → `count_*_boundary` on `CalibrationResult` (crossing incidences, 23 files): a rename
  gated by the identity check (content, never names), the suite and a golden regeneration that must show only
  column names moving.
* the wave-3 frame migration: all six bank-readers still read pass one's `cache.payload` and none calls
  `calibration_inputs` — `structural_claims_audit`, `gdna_pool_census`, `abundance_landscape_census`,
  `transport_dispersion`, `fl_pool_purity`, `calibration_truth_ab`. Each moves to the drained frame the truth
  is certified in (`calibration_inputs(cache, index)["payload"]`, as the oracle instrument does), one at a time;
  their numbers WILL move with the frame, so the gate is each self-test plus a before/after recorded in the
  commit. `calibration_truth_ab` prints undrained against drained on purpose: read it before migrating it.
* the moment tests deleted with the length channel: find them with `git log --diff-filter=D -- tests` and
  `git log -S"moment"`, restore those that gate something still live.
* the index's duplicate map as an alias map `dropped_t_id → kept_t_id`: an index rebuild, no panel re-scan,
  verified with `rescan_panels.py` (`reach` is covered by no other hash). The largest item; last.

**What closes W8.** The two rename entries closed with their rulings (a kept word is a ruling too), the ledger's
items ticked in place, every identity check bit-identical, the suite at its count (renames add no file), the
docs following each rename in the same commit (the move rule), `preflight --full`.

## 8. The remaining steps before the port (agreed 2026-09-13)

Plan §F's "satisfied with the Python" is met: memory, one representation for the received messages, the
layering, the gates and the tolerance instrument, the docs. Four items remain, in this order, and then F.

**W10 — the hygiene ledger's pure cleanups**, each its own commit, none moving the metric:
1. DONE 2026-09-13. the stale comment in `pipeline.py` above the second pass's `build_fl_models` (the text "The SAME
   de-tilted RNA pool the calibrator will read": the second pass's fl models are pass one's; the lift's
   docstring above it has it right) — three lines, the suite;
2. DONE 2026-09-13, EXTENDED to the two region arrays: four instruments build the field name from the axis (`f"…_{axis}"`), so the prefix must be one word — `count_<population>_<axis>` throughout, the form `count_rna_sj` already had. `mass_*_boundary` → `count_*_boundary` on `CalibrationResult` (crossing incidences, 23 files) — a rename
   gated by `rename_identity.py --check` (content, never names), the suite, and a golden regeneration that
   must show only column names moving;
3. the wave-3 frame migration, RULED 2026-09-13 (owner, after an instrument census: 44 files, 19,532 lines;
   none of the six is loaded by a test or named by the roadmap): four of the six bank-readers are RETIRED
   rather than migrated, each its own commit, gated by the suite (−4 collected per file) and preflight —
   `structural_claims_audit` DONE 2026-09-13; `gdna_pool_census` DONE 2026-09-13; `calibration_truth_ab` DONE 2026-09-13;
   `abundance_landscape_census` DONE 2026-09-13 — and two migrate to the drained frame the truth is certified in
   (`calibration_inputs(cache, index)["payload"]` plus the partition lift where partitions are read):
   `fl_pool_purity` DONE 2026-09-13 (an open defect at priority next names it; it priced a length model production
   never builds — after the migration the shipped model reads the true gDNA mean exactly at capture-OFF on both
   fl-gap arms, and −2.9 / −8.5 bp at capture-ON) and `transport_dispersion` DONE 2026-09-13 (the decomposition the owner asked for; an open question waits on
   it; in the drained frame the deep pairs' excess disagreement falls ~10 % and the pooled centre halves), each
   with a recorded before/after in its commit;
4. DONE 2026-09-13 (33 cases into `test_effective_length.py`; both perturbations fired). the moment tests deleted with the length channel: `git log --diff-filter=D -- tests` and
   `git log -S"moment"`, restore those that gate something still live.
The index alias map (`dropped_t_id → kept_t_id`, an index rebuild verified with `rescan_panels.py`) stays in
`ISSUES: hygiene-ledger` as a later item.

**W11 — the coverage census.** DONE 2026-09-13. The suite and every instrument's `--self-test` ran under
coverage; the never-executed lines of `src/` were reviewed one by one and ruled dead, rotten-but-live, or a
coverage gap; every dead item went in its own commit, byte-identical on the default rows of both substrates
(`arms/2026-09-13_preport/oracle_*_default.json`) and on the three identity references. The record is in
`ISSUES: hygiene-ledger`; the coverage gaps kept are listed there too.

**W12 — the θ quadrature at zero gDNA** (`ISSUES: theta-quadrature-at-zero-gdna`). The one open design item
inside ψ, and the port will carve ψ into C++, so it is settled in Python first. What is known: K_t 30 breaks the
zero control 6.3× on one row (`g00 ss.99 capture ON`, with the message layer off as well), 120 and 240 equal
60; a Chebyshev mesh clustered at the ends is refused (9.9× at 30 nodes); the per-slot bias is small
(≤ 4 fragments at 50k, prior-free) and the refits amplify it; the peak-width analysis does not predict the
pipeline, so the bias is a λ-DEPENDENCE of the quadrature's error, not its size. Begin with that derivation;
the design target is a quadrature whose accuracy does not depend on node count, which is also the cube's
cost lever (K × K_t). DERIVE → PROTOTYPE outside `src/` (patch `_tilt_grid` and `_psi` in both their bindings,
as `w5/theta_mesh.py` did) → A/B on both panels with both zero controls → then `src/`.

**W13 — the port's prerequisites.** DONE 2026-09-13: `sweeps_MO_3021_step5` captured on `23a431a9` and every
call replays bit-identical; `port_identity_*` frozen and checked (two ladder conditions and LBX0190). The
earlier captures and references describe trees before the θ quadrature. Then F.

## 4. Two things not to do

* Micro-optimise the Python passes (a silent-hop early exit would halve them, bit-identically): the port
  deletes it.
* Bake the refit grid into anything.

## 5. Names — the rule and the proposals

**The rule.** A name says what the thing is in the words the interface already uses; no new noun for a
thing that has one; nothing the owner has to look up. The established words: TERMINAL (a region that
receives nothing), LOCUS and BLOCK (`LocusBlock`), the LIBRARY (the whole-chain reductions), a FACE and its
SIDE (left = 0, right = 1), a rule's KIND, a node's OWN CLAIM and OWN LEVEL, a LANE, a COMPOSITION and a
LEVEL, `from_left` / `from_right`, `SILENCE` / `NO_NEIGHBOUR`, the FACTORY ROWS, `solvable`.

**The naming pass, RULED and APPLIED (owner, 2026-09-11; bit-identical on every gate):**
`MessageMemo` → `MessageCache` (`message_memo=` → `message_cache=`); `StepContext` → `BlockContext`;
`StepContext.own_live` → `BlockContext.has_own_composition` — the bit is `tau_lam > 0`, the node's OWN
composition evidence (strand term OR factory row), not a strand channel, so the plan's `strand_live` was
wrong in substance; `ChainView.strand_live` and the transfer policy's private `_Library.split_live` keep
their names (three liveness bits, one renamed removes the collision); `RegionBelief.informed` →
`has_composition` (a vocabulary ruling: the same word carries the same fact from a received table to the
prior's training population); `_slots` → `_block_slice` (`_block_view` refused: "view" already means
`ChainView` in the same file); `_view_fields` and `_gather` keep their names (`_context_fields` refused:
the function builds `ChainView`'s fields and serves the library's whole-chain view too). The split's file
names: `sweep.py` (backbone), `blocks.py`, `message_cache.py`, `messages/lanes.py`.
