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

* **Every step is proven before it is believed.** Bit-identity: `sweep_replay.py replay --dir
  ~/Downloads/rigel_runs/perf/sweeps_MO_3021_step2 --call 0..3` (and `--block-slots N|none`) and
  `rename_identity.py --check --reference ~/Downloads/rigel_runs/arms/locus_identity_*.json` (two ladder
  conditions and `--bam` LBX0190). The suite (`CLAUDE.md`'s baseline line), `ruff`, `preflight.py --full`.
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
   lower priority: the 25-key `_capture` dict as a typed record (moderate churn).

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
| W5 | **One grid?** | with one solver and tiled memory, does one λ grid serve both classes, deleting the per-tile regrid of the priors? The refit grid (`sweep_n_grid_single_strand`) is an accuracy ruling, so this is an A/B | the oracle metric per stratum |
| W6 | **The tunables census** | `CalibrationConfig`'s 49 fields: live / derived / dead, each derivation named; anything unearned removed | `module_census.py`, the suite |
| W7 | **The capture as a typed record** | the 25-key `_capture` dict (75 keyword lines in `_solve_block`) | the instruments that read it (`landscape_training_census.py`, `backbone_parity.py`, the walk) |
| W8 | **Vocabulary rulings** | `ISSUES: rename-the-drain`, `ISSUES: rename-row-and-face`; the `hygiene-ledger` items | `rename_census.py --sense`, `rename_identity.py --check` |
| W9 | **The two xfails** | `ISSUES: two-sided-exon-row` (priority now; the fix is the landscape's enrichment witness) and the antisense prior-assembly casualty | each xfail's own test turning green structurally |

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
the refit sweeps' grid a parameter (`n_grid_ss = 513` on the refits is an accuracy ruling, not a
performance one). "Satisfied with the Python" means: A–D landed, one representation for the received
messages, the layering clean, the gates and the tolerance instrument in place, the docs current.

### G. Then the scan and the second pass (⑤)

31 s and 22 s at 8 threads on the deep library; the stages that scale with depth and the whole problem at
100M+ fragments. The scan's thread split is the owner's decision (`ISSUES: scan-thread-split-starves-the-workers`).

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
