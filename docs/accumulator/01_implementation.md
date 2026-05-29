# Fractional Accumulator — Phased Implementation Plan

Companion to [`00_design.md`](00_design.md). This document is the
step-by-step playbook. Each phase has explicit file lists,
acceptance gates, and rollback plans. **No phase begins until the
previous phase's acceptance gate is green and is reviewed.**

---

## Cross-cutting conventions

- **Build between every phase.** `pip install --no-build-isolation -e .`
  must succeed before any test command runs.
- **Tests scoped per phase.** Each phase adds its own test file
  under `tests/native/` or `tests/` and is mergeable as a stand-alone
  commit. Never delete a legacy test until the parity gate (Phase 5)
  signs off.
- **No premature deletions.** Legacy code (`RegionCountLedger`
  boundary fields, [`boundaries.py`](../../src/rigel/calibration/boundaries.py))
  remains compiled-and-callable through Phase 3. Cutover
  (Phase 4) replaces the legacy accumulator in place; downstream
  consumers migrate in the same phase. Dead-code deletion happens
  in Phase 5.
- **FL is not accumulated.** Per audit_phase1.md decision #6, the
  the substrate does not collect per-region or per-boundary FL
  histograms. FL stays in the downstream EM scorer
  ([`scoring.py`](../../src/rigel/scoring.py),
  [`frag_length_model.py`](../../src/rigel/frag_length_model.py))
  where it can be applied per-fragment without binning.
- **Naming.** Old type: `RegionCountLedger`. New type:
  `Accumulator`. Old Python adapter: `RegionCountsPayload`
  (whatever it's called today). New: `AccumulatorPayload`. The
  legacy types are removed entirely during Phase 4 cutover.

---

## Phase 1 — Audit & spec lock

**Goal.** Document the *exact* current per-fragment accounting
semantics in `src/rigel/native/accumulator.cpp` and
`bam_scanner.cpp`. Fix any per-block double-counting bugs found
**inside the legacy implementation** so the parity gate in Phase 5
compares two correct accumulators. Lock the accumulator spec by writing the
reference test corpus.

### 1.1 Code review tasks (read-only)

Files to read top-to-bottom and annotate with a written memo
(`docs/accumulator/audit_phase1.md`):

1. `src/rigel/native/accumulator.h`
2. `src/rigel/native/accumulator.cpp`
3. `src/rigel/native/bam_scanner.cpp` — focus on the per-fragment
   loop that calls accumulator deposit helpers.
4. `src/rigel/native/region_count_ledger.{h,cpp}` (or equivalent)
5. `src/rigel/calibration/boundaries.py`
6. Tests that touch the accumulator:
   `tests/test_region_count_ledger.py`,
   `tests/test_region_partition.py`,
   `tests/test_region_persist.py`,
   `tests/test_boundary_model.py`,
   `tests/test_boundary_sweep.py`.

### 1.2 Memo contents (`audit_phase1.md`)

Required sections:

1. **Per-fragment vs per-block.** Quote each line in
   `accumulator.cpp` where a deposit occurs. State whether it is
   inside a per-block loop or a per-fragment loop, and whether the
   normalization divides by per-block $\ell_k$ or per-fragment $L$.
2. **Boundary-mass-mirroring inventory.** Today's
   `RegionCountLedger` stores `boundary_left_*` on region $r$ and
   `boundary_right_*` on region $r-1$ for the same boundary. Confirm
   that the two slots are never both incremented independently.
3. **Fully-spans-region path.** Document the current logic for
   fragments where one block straddles a region. Today this triggers
   the `w/2` split rule (see [`00_design.md`](00_design.md) §2.3).
4. **Spliced-flag attribution.** State exactly when a fragment is
   marked spliced vs unspliced in the boundary deposit. Compare to
   the spec (per-junction-event spliced flag).
5. **Strand attribution.** State how strand $s \in \{+, -\}$ is set
   per fragment.
6. **Found bugs.** Enumerate any discrepancies between the
   documented intent and the code. Each bug gets a ticket-style
   entry: file:line, observed behavior, expected behavior, fix
   proposal.

### 1.3 Fix any bugs found in legacy code

If §1.2.6 surfaces bugs, fix them **inside the legacy
`accumulator.cpp`** with their own commits and tests in this phase
(not the new accumulator). The point is to make the parity gate in Phase 5 meaningful.

### 1.4 Lock the accumulator spec via reference test corpus

Create `tests/native/test_accumulator_spec.py` with **expected
outputs only** (the native implementation does not exist yet; tests are
xfail-marked until Phase 2 lands). Each test is a tiny synthetic
fragment list paired with the expected `AccumulatorPayload`
contents, computed by hand or by a Python reference implementation
in the test file itself.

Tests required (from [`00_design.md`](00_design.md) §9):

| ID | Name |
|---|---|
| T1 | `test_contained_single_block` |
| T2 | `test_contained_multi_block_spliced` |
| T3 | `test_two_block_adjacent_regions` |
| T4 | `test_two_block_non_adjacent_regions` (user-verified §4.5.1 walkthrough) |
| T5 | `test_three_block_all_adjacent` (user §4.5.2 walkthrough) |
| T6 | `test_fully_spans_region` |
| T7 | `test_mass_conservation_random` |
| T8 | `test_flux_dedup_adjacent_regions` |
| T10 | `test_strand_attribution_negative_strand` |
| T11 | `test_spliced_flag_per_junction` |

Plus a pure-Python reference implementation in
`tests/native/_accumulator_reference.py` (~150 LoC) that any
human can read and confirm matches the spec. Phase 2's native
implementation will be tested against this Python reference.

### 1.5 Acceptance gate

- [ ] `docs/accumulator/audit_phase1.md` exists with all 6 required sections.
- [ ] All bugs found in §1.2.6 either fixed (with tests) or filed as deferred (with justification).
- [ ] `tests/native/test_accumulator_spec.py` exists with all 11 tests, currently `pytest.mark.xfail(reason="accumulator not yet implemented")`.
- [ ] `tests/native/_accumulator_reference.py` is a working Python implementation.
- [ ] `pytest tests/ -v` passes (legacy tests still green; new spec tests xfail).

### 1.6 Rollback

This phase is read-only / docs / fix-legacy-bugs. Rollback = revert
the (small) legacy-bug-fix commits if they regress something.

---

## Phase 2 — Native rewrite (`accumulator.{h,cpp}`)

**Goal.** Land the C++ implementation. Compile alongside legacy
accumulator (both are live in the binary). No call sites switched
yet.

### 2.1 Files added

1. `src/rigel/native/accumulator.h` — full struct + class API per
   [`00_design.md`](00_design.md) §3, §7.1.
2. `src/rigel/native/accumulator.cpp` — implementation.

### 2.2 Files modified

- `src/rigel/native/CMakeLists.txt` — add the new `.cpp` to
  `_bam_impl` target sources. Bind module unchanged in this phase.

### 2.3 Implementation outline (`accumulator.cpp`)

```cpp
void Accumulator::deposit(const FragmentEvent& fe) {
  // fe carries: ref_id, blocks[K] (start, end), strand_sign,
  //             per_block_region_id[K], per_junction_spliced[K-1].

  // 1. Compute L = sum(block lengths).
  // 2. Classify contained vs crossing:
  //      crossing = any (per_block_region_id[k] != per_block_region_id[0])
  //                 OR any block straddles a boundary (rare; §4.5.4).

  if (!crossing) {
    int r = per_block_region_id[0];
    bool spliced_any = std::any_of(...);
    // Increment contained channel.
    if (spliced_any) {
      if (strand_pos) ++regions_[r].contained_spl_pos;
      else            ++regions_[r].contained_spl_neg;
    } else {
      if (strand_pos) ++regions_[r].contained_unspl_pos;
      else            ++regions_[r].contained_unspl_neg;
    }
    return;
  }

  // Crossing path.
  // For each junction (k, k+1):
  std::array<int, MAX_JUNCTIONS> touched_boundaries{};
  int n_touched = 0;
  for (int k = 0; k + 1 < K; ++k) {
    int Rk   = per_block_region_id[k];
    int Rkp1 = per_block_region_id[k+1];
    int b_out = right_boundary_of(Rk);        // = ref_boundary_offset + local_region_index_k + 1
    int b_in  = left_boundary_of(Rkp1);       // = ref_boundary_offset + local_region_index_kp1
    bool spliced = per_junction_spliced[k];
    float mass_k   = float(block_len[k])   / float(L);
    float mass_kp1 = float(block_len[k+1]) / float(L);

    // block-on-left (B_k is left of b_out) -> mass_left
    deposit_left (b_out, strand_pos, spliced, mass_k);
    // block-on-right (B_{k+1} is right of b_in) -> mass_right
    deposit_right(b_in,  strand_pos, spliced, mass_kp1);

    // Flux dedup: each touched boundary gets one flux per fragment-junction.
    flux_increment_dedup(touched_boundaries, n_touched, b_out, strand_pos, spliced);
    flux_increment_dedup(touched_boundaries, n_touched, b_in,  strand_pos, spliced);
  }
}
```

`flux_increment_dedup` checks whether the boundary index is already
in `touched_boundaries[0..n_touched]`; if not, appends and
increments. The "channel" granularity (spl_pos/spl_neg/etc.) is
finer than the dedup key, but per [`00_design.md`](00_design.md)
§4.4 the dedup is at boundary level — a single junction event over
adjacent regions is one observation, one flux increment.

### 2.4 The fully-spans-region case (§4.5.4)

Block $B$ straddles boundaries $B_a, B_b, \ldots$. Decompose $B$
into per-region intersections $\ell^{(R)}$ for each region $R$ the
block overlaps. Treat the resulting sequence of $(R, \ell^{(R)})$
slices **as if they were adjacent blocks with no inter-block gap**
and apply the §4.3 rule. Mark each implicit boundary-crossing as
unspliced (it's a contiguous CIGAR M run, not an N).

Implementation: a helper `expand_block_to_region_slices(block,
out)` runs before the main per-block loop. After expansion, the
deposit loop is identical to the multi-block case.

### 2.5 Tests added

- `tests/native/test_accumulator_spec.py` — remove `xfail` markers.
- `tests/native/test_accumulator_against_python_reference.py` —
  property tests: 1000 random synthetic fragment streams, assert
  the native `Accumulator` output matches the Python reference
  byte-for-byte (after the float32→bit equality dance).

### 2.6 Acceptance gate

- [ ] `pip install --no-build-isolation -e .` succeeds.
- [ ] All 11 spec tests pass (no `xfail`).
- [ ] Python-reference property tests pass on 1000 random fragments.
- [ ] Legacy tests still pass (we haven't switched callers yet).

### 2.7 Rollback

Drop the new source files from CMakeLists.txt. Nothing else is
touched.

---

## Phase 3 — nanobind binding + Python adapter

**Goal.** Expose `Accumulator` to Python via the existing
`_bam_impl` nanobind module. Add `AccumulatorPayload` dataclass.
Still no call sites switched.

### 3.1 Files added

1. `src/rigel/scan_payload.py` — augmented (or new module if missing)
   with the `AccumulatorPayload` dataclass from
   [`00_design.md`](00_design.md) §7.3.

### 3.2 Files modified

- `src/rigel/native/bam_scanner.cpp` — nanobind binding block at the
  bottom: expose `Accumulator::regions`, `boundaries`, and offset
  arrays as numpy arrays (zero-copy with `nb::ndarray`).
- `src/rigel/native.py` — re-export the new bindings.

### 3.3 Numpy view contract

Each of the 16 boundary mass/flux fields and 4 region contained
fields is exposed as a separate 1-D numpy array. Lifetime is tied
to the `Accumulator` C++ object via nanobind's owner-reference
machinery.

### 3.4 Tests added

`tests/test_accumulator_payload.py`:

- Construct a tiny synthetic accumulator from Python (using a thin
  Python factory exposed by the binding).
- Confirm field dtypes (`uint32`, `float32`), shapes, and contiguity.
- Confirm zero-copy: writes to the numpy view mutate the C++ state
  (or, if we choose to expose read-only views, writes raise).
- Confirm garbage collection: the numpy arrays keep the underlying
  C++ object alive until the last reference is dropped.

### 3.5 Acceptance gate

- [ ] Build succeeds.
- [ ] `pytest tests/test_accumulator_payload.py -v` passes.
- [ ] No regression in existing tests.

### 3.6 Rollback

Revert the nanobind binding block in `bam_scanner.cpp` and remove
`AccumulatorPayload` from `scan_payload.py`.

---

## Phase 4 — In-place cutover

**Goal.** Replace the legacy `RegionCountLedger` boundary path with
the new `Accumulator`. This is a single atomic cutover: the legacy
deposit call is removed in the same commit (or PR) that migrates
every downstream consumer. Before starting, snapshot the
pre-cutover golden outputs for the Phase 5 regression gate.

### 4.1 Snapshot pre-cutover golden (read-only)

Before any code change:

1. Build the **current `main`** binary (legacy accumulator still
   live).
2. Run rigel on 3 representative conditions:
   - `gdna_none_ss_1.00_nrna_none`
   - `gdna_high_ss_0.50_nrna_rand`
   - `gdna_low_ss_0.90_nrna_none`
3. Save the resulting transcript-level `quant.feather` files to
   `runs/<condition>/rigel/precutover_golden/`. These are the
   reference for Phase 5.

### 4.2 Consumer audit (read-only)

Run `grep_search` for these symbols and enumerate every callsite
that needs migration:

1. `BoundaryTable` — primary consumer in `src/rigel/calibration/`.
2. `build_boundary_table`, `validate_boundary_table` — from
   `src/rigel/calibration/boundaries.py`.
3. `RegionCountLedger` boundary-slot reads (anything reading
   `boundary_left_*` or `boundary_right_*` off the ledger).
4. `compute_locus_priors_from_partitions` — confirm whether it reads
   from boundary table or directly from ledger.
5. Test fixtures that construct synthetic `BoundaryTable`s — they
   get rebuilt against synthetic `AccumulatorPayload`s using the
   Python reference impl from Phase 1.

### 4.3 Files modified

- `src/rigel/native/bam_scanner.cpp` — replace the legacy
  `RegionCountLedger` deposit with `accumulator_.deposit(fe)`. No
  feature flag, no dual-write — the legacy call is removed.
- `src/rigel/scan.py` — replace the legacy payload surface with
  `AccumulatorPayload`.
- `src/rigel/calibration/_simple.py`, `_categorize.py`,
  `src/rigel/locus.py`, `src/rigel/priors.py` — switch reads to
  `AccumulatorPayload`. Exact set to confirm via §4.2 grep.
- Test fixtures (`tests/test_boundary_model.py`,
  `tests/test_boundary_sweep.py`, `tests/test_region_count_ledger.py`,
  etc.) — update to the new substrate.

### 4.4 Strategy

A single PR with these commits, in order:

1. Replace deposit call in `bam_scanner.cpp`. Build must still
   succeed; the legacy `BoundaryTable` consumers will fail at this
   point until commit 2 lands.
2. Migrate each consumer in sequence, smallest-blast-radius first.
3. After every commit, run `pytest tests/ -v`; partial failures are
   expected mid-PR, but the final commit must restore green.

### 4.5 Acceptance gate

- [ ] Pre-cutover golden snapshots saved to
      `runs/<condition>/rigel/precutover_golden/`.
- [ ] Build succeeds.
- [ ] No remaining read of `RegionCountLedger.boundary_left_*` or
      `.boundary_right_*` fields.
- [ ] Full `pytest tests/ -v` passes.

### 4.6 Rollback

Revert the cutover PR. The pre-cutover snapshot path remains.

---

## Phase 5 — Regression gate against pre-cutover golden

**Goal.** Re-run rigel on the 3 representative conditions with the
post-cutover binary and confirm transcript-level outputs match the
Phase 4.1 snapshot. This is the single critical-path correctness
gate of the project.

### 5.1 Tasks

1. Re-run rigel on the same 3 conditions used in §4.1. Output goes
   to `runs/<condition>/rigel/accumulator_cutover/`.
2. Diff transcript-level `quant.feather` against
   `runs/<condition>/rigel/precutover_golden/`.
3. Allowed numerical noise: ~1e-6 relative tolerance on
   `mrna_abundance`. Larger discrepancies must be root-caused.

### 5.2 Diagnostics (`scripts/debug/accumulator_parity.py`)

Optional helper for root-causing divergences. Loads an
`AccumulatorPayload` and the saved legacy `RegionCountLedger`-style
projection (kept as a one-off serialization captured in §4.1) and
reports per-boundary, per-region, per-channel diffs.

### 5.3 Acceptance gate

- [ ] All 3 conditions match the pre-cutover golden within 1e-6
      relative tolerance on `mrna_abundance`.
- [ ] If any diff exceeds threshold: revert Phase 4 and root-cause
      first. Most likely culprits: implicit-splice classification
      mismatch ([`00_design.md`](00_design.md) §4.3 spliced flag), or
      fully-spans-region attribution (§4.5.4).

### 5.4 Rollback

Revert Phase 4 PR.

---

## Phase 6 — Dead-code deletion

**Goal.** Remove the legacy boundary path and any now-orphaned
helpers. The contained counts may stay if `RegionCountLedger` still
serves other purposes; otherwise, delete the whole ledger.

### 6.1 Files deleted

- `src/rigel/calibration/boundaries.py`
- `tests/test_boundary_model.py`, `tests/test_boundary_sweep.py`
  (or rewrite as tests of `AccumulatorPayload` if their coverage is
  still valuable).
- `RegionCountLedger::boundary_left_*` / `::boundary_right_*` fields
  in `src/rigel/native/region_count_ledger.h`.

### 6.2 Files modified

- `src/rigel/native/region_count_ledger.cpp` — strip boundary
  accounting helpers.
- `tests/test_region_count_ledger.py` — strip boundary tests.

### 6.3 Documentation updates

- `CLAUDE.md` — update the "Python Module Roles" table: remove
  `boundaries.py` row, update `scan_payload.py` description.
- `.github/copilot-instructions.md` — same.
- `docs/caljointmodel/04_interface_contract.md` §2 — remove the
  "Status (2026-05-29)" pending banner; mark substrate as live.
- `CHANGELOG.md` — entry for the rewrite.

### 6.4 Acceptance gate

- [ ] Full `pytest tests/ -v` passes.
- [ ] `ruff check src/ tests/` clean.
- [ ] Phase 5 regression gate re-run still passes.
- [ ] No grep hits for `BoundaryTable`, `boundaries.py`,
      `boundary_left_`, `boundary_right_`.

### 6.5 Rollback

Revert the deletion commit if any test surface lights up.

---

## Phase 7 — Calibration v6 substrate hookup

**Goal.** Re-open [`../caljointmodel/`](../caljointmodel/) and wire
the calibrator to consume `AccumulatorPayload` directly per docs
01–03. This phase is *out of scope* for the fractional accumulator
plan itself but is the immediate downstream consumer; flag it here
so the rewrite doesn't drift away from its purpose.

### 7.1 Trigger

Phase 6 acceptance gate green.

### 7.2 First task

Read [`../caljointmodel/05_implementation_plan.md`](../caljointmodel/05_implementation_plan.md)
Phase 3 ("substrate adapter") and re-scope: the substrate is now
the live `AccumulatorPayload`, not an extension. Update that doc to
mark it as completed-by-accumulator-rewrite.

---

## Risk register summary

| Phase | Top risk | Mitigation |
|---|---|---|
| 1 | Audit memo skipped | Required artifact at gate |
| 2 | Fully-spans-region case is subtle | Pre-decomposition helper + dedicated test T6 |
| 3 | Nanobind lifetime bugs | Explicit owner-reference + GC test |
| 4 | Mid-PR build/test breakage during cutover | Single PR, ordered commits, final commit restores green |
| 5 | Numerical divergence vs pre-cutover golden | Diagnostics script + per-boundary report; revert Phase 4 if structural |
| 6 | Deletion regresses unrelated tests | Revert commit; tests are rebuildable |
| 7 | Out of scope | Documented as future work |

---

## File-touch matrix

| File | P1 | P2 | P3 | P4 | P5 | P6 | P7 |
|---|---|---|---|---|---|---|---|
| `accumulator.{h,cpp}` (legacy) | read+fix | **replace** | — | — | — | delete residual boundary fields | — |
| `bam_scanner.cpp` | read | edit (CMake) | edit (binding) | edit (replace deposit call) | — | — | — |
| `scan.py` | — | — | — | edit | — | — | — |
| `scan_payload.py` | — | — | **add** | edit | — | — | — |
| `native.py` | — | — | edit | — | — | — | — |
| `calibration/boundaries.py` | read | — | — | — | — | **delete** | — |
| `calibration/_simple.py` | read | — | — | edit | — | — | — |
| `calibration/_categorize.py` | read | — | — | edit | — | — | — |
| `locus.py` | read | — | — | edit if needed | — | — | — |
| `priors.py` | read | — | — | edit if needed | — | — | — |
| `tests/native/test_accumulator_spec.py` | **add (xfail)** | unxfail | — | — | — | — | — |
| `tests/native/_accumulator_reference.py` | **add** | — | — | — | — | — | — |
| `tests/test_accumulator_payload.py` | — | — | **add** | — | — | — | — |
| `tests/test_boundary_model.py` | read | — | — | rewrite for new payload | — | delete | — |
| `tests/test_boundary_sweep.py` | read | — | — | rewrite for new payload | — | delete | — |
| `scripts/debug/accumulator_parity.py` | — | — | — | — | optional **add** | — | — |
| `docs/accumulator/audit_phase1.md` | **add** | — | — | — | — | — | — |
| `docs/caljointmodel/04_interface_contract.md` | — | — | — | — | — | edit (drop banner) | — |
| `CLAUDE.md` | — | — | — | — | — | edit | — |
| `.github/copilot-instructions.md` | — | — | — | — | — | edit | — |
| `CHANGELOG.md` | — | — | — | — | — | edit | — |

---

## Sign-off

This plan covers Phase 1 through Phase 7. Approval to begin Phase
1 = approval of the plan. Each phase gate is its own checkpoint.
