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
   its own slice of every input, reading one terminal beyond itself (`_solve_block`, `_slots`, `_gather`).
2. **The library is the only cross-block information.** `Policy.library(ChainView)` runs once over the
   whole chain on a view with no belief field; `prepare(ctx, library)` sees one block. The factory rows
   travel on the context (`factory_rows`); the context carries the strand channel's liveness bits
   (`own_live`), not the self-solve object.
3. **ψ's read-out is chunk-exact.** The one step that moved a number (≤ 3.1e-15 per slot, no
   amplification, TPM bit-identical, the oracle metric at the last ulp; owner: identical to a tolerance).
   Every block size gives the same bits, so the block size is a working-set knob.
4. **The refit sweeps share their message layer** (`sweep.MessageMemo`). Everything the layer reads is on
   the context, the library and the grid — never the prior — and the belief is reset before every sweep,
   so sweeps 1–3 deliver identical messages; the content-keyed memo serves refits 2–3 their rows
   (38 s instead of 176 s each) for 2.68 GB held.
5. **The rules are typed tables** (`messages.transfer.Faces`): every directed face is one of a node's
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

1. **The naming pass** (§5 lists the proposals — present them to the owner at the start of the session,
   apply unless overruled): `MessageMemo` → `MessageCache`; `ChainView.strand_live` (the library's
   deadband is open) and `StepContext.own_live` (this node's channel is live) are two near-identical names
   for two scopes — `strand_channel_open` and `strand_live`; `StepContext` → `BlockContext`; `_slots` →
   `_block_view`; `_view_fields` → `_context_fields`; and a question for the owner: `RegionBelief.informed`.
2. **Split `sweep.py`** (~1,000 lines, three concerns) into the backbone (`_pass`, `_check_message`,
   `solve_chain`, `_solve_block`), the block plumbing (`blocks.py`: `_block_view`, `_context_fields`,
   `_gather`) and the cache (`message_cache.py`); `messages/transfer.py` (1,200 lines): the lanes
   (`_LevelLane`, `_gdna_lane`, `_rna_lanes`) into `messages/lanes.py`. Layer 6 throughout; `_layers.py`
   lists them; `test_layering.py` and `module_census.py` confirm.
3. **The cache key iterates `ChainView`'s dataclass fields** rather than a hand-built dict, so a new field
   cannot be left out of the digest (correctness, not tidiness; the six perturbation channels in
   `test_sweep_backbone` gate it).
4. **Drop the redundant single-strand tiling** in `_solve_regions_logodds_all` inside a block
   (`_block_rows` stays for the AMBIG cube only): recovers the ~7 s/sweep the final ψ lost at 5,000-slot
   blocks and removes a layer of nesting. Chunk-exact, so bit-identical; the tiling gate in `test_sweep`
   holds it.
5. Small things: `build_region_init`'s dead `chain` parameter; `n_slot` / `spliced_slot` carried beside the
   arrays they derive from (properties on the context); later, the 25-key `_capture` dict as a typed record
   (moderate churn, lower priority).

### B. The received messages as tables — the contract, end to end (owner: agreed)

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

### C. The tolerance gate — before any compiled code exists

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

### F. The port (③) — only when A–D are done and the design is agreed

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

**Proposals for the naming pass** (A.1): `MessageMemo` → `MessageCache` (and `message_memo=` →
`message_cache=`); `ChainView.strand_live` → `strand_channel_open`; `StepContext.own_live` → `strand_live`;
`StepContext` → `BlockContext`; `_slots` → `_block_view`; `_view_fields` → `_context_fields`. **For the
owner:** `RegionBelief.informed` — "this slot's answer rests on a composition (its own evidence, structural
certainty, or one received), not only on a bound", the landscape prior's training population
(`DESIGN.md` §7.1) — is established in the code and docs; if it is to change, `has_composition` is the
plain candidate, and it is a vocabulary ruling (`DESIGN.md` §0).
