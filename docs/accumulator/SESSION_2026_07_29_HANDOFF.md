# Accumulator v5 — SESSION HANDOFF, 2026-07-29

> # ⛔ SUPERSEDED BY `SESSION_2026_07_30_HANDOFF.md` — START THERE.
> This file is kept only for its §1 (the design in one page) and §4 (things not to re-derive), both
> still accurate. Its §6 roadmap is stale: W1b, W1b-clean, W1c and the W1 review have all landed.

> ## ⛔ SUPERSEDED IN PART — W1b AND W1b-clean HAVE LANDED. Read this box, then §0.
>
> **Owner directive, 2026-07-29: no backwards compatibility, no legacy retention, converge and
> delete.** Two consequences that change the roadmap in §6:
>
> * ✅ **W1b LANDED**, and its premise was wrong: the 32-condition bench **cannot see a partition
>   change**. `ambig_dense_10mb`'s v8 node partition is byte-identical to its v7 region partition
>   (1,698 == 1,698, row-for-row) — the suite has no alternative TSS/TES. The gate inverted to
>   *bit-identical*, and it held 32/32 at both refits. Goldens did **not** move. `partition_hash`
>   did **not** change there, so no re-scan. **THE PARTITION EFFECT was measured on
>   `quick_3to1_5mb`** (+25.0 % nodes): **+0.0069 at pass-0, −0.0014 at refit=3**, splitting on
>   whether the library has gDNA at all. Ledger W1b.
> * ✅ **W1b-clean LANDED — the v7 partition is GONE.** `calibration/regions.py` deleted whole,
>   `regions.feather`/`boundaries.feather` out of the index, `index.region_df`/`.boundary_df` gone,
>   P3 replaced by the stronger **I3b**. Bit-identical 32/32. src nets −434 lines. Ledger W1b-clean.
> * ⛔ **SHADOW MODE (W4) IS CANCELLED.** W3 is unchanged; W4 replaces the C++ accumulator outright,
>   gated on the Python reference and the oracle bench.
>
> * ✅ **W1c LANDED** — the graph's 8 structural flags reach `NodeStatics` (raw `uint16`, not
>   pre-derived predicates). Bit-identical 32/32. ⛔ **Its own test found a real bug in the W1a
>   builder:** plan F1's flags filter `~is_synthetic & ~is_nrna` deleted the TSS/TES of **26,475 real
>   single-exon transcripts** (52,104 terminus positions) — on a non-synthetic row `is_nrna` means
>   "single-exon", not "manufactured span". Now ONE filter, `~is_synthetic`. That revealed a missing
>   invariant, now shipped as **I13** (the flags ARE the events, both directions, with a bookended-exon
>   exemption). Ledger W1c.
>
> * ✅ **W1-REVIEW LANDED** — 5 adversarial lenses, 64 agents, 42 confirmed findings, all fixed,
>   bit-identical. ⛔ The one that mattered: **I13 was supplying FALSE assurance** — it called the
>   builder's own event emitter, so deleting the NEG-strand TSS/TES swap passed the entire suite AND
>   `validate_graph`. It now re-derives the events by a different algorithm and catches it. ⭐ Also:
>   `validate_graph` at load **1.91 s → 0.23 s**, build **4.6 s → 3.1 s**, and the human termini
>   figure is **53.4 %, not 59.5 %** (the old number was computed under the buggy filter). Ledger
>   W1-REVIEW.
>
> **Baseline in force: r0 `0.079005` / r3 `0.046675` on the v8 partition** (numerically identical to
> the v7 one on this suite — the objects are the same). Suite **1280 pass**, ruff + format clean.
> **Next: W2a/b/c** (P1g off the index alone), then W3.

**Read this first, then `06_implementation_plan.md`, then `accumulator_ledger.md`.**
Everything below was measured or run in the session that wrote it. Nothing is quoted from an older doc
without being re-verified.

---

## 0. STATE IN TEN LINES

| | |
|---|---|
| branch / HEAD | `main` @ `3c293038`. **Nothing committed — the owner drives commits.** |
| test suite | **1281 pass**, ruff clean on `src/ tests/ scripts/`, `ruff format` parity, warning-clean under `-W error::RuntimeWarning` |
| baseline (32-cond `ambig_dense_10mb`, mass-weighted mwae) | **r0 `0.079005` / r3 `0.046675`** — re-recorded from this tree and reproduced **32/32 exactly** |
| phases done | **W0** (instruments + the D6 arm) · **W0b** (coarsening map) · **W1a** (v8 graph, gate passed) |
| next phase | **W1b** — wire the scanner to the v8 node cut array. The first arm expected to MOVE numbers. |
| C++ rebuilt? | yes (`merge_from` binding). Re-run `pip install --no-build-isolation -e .` after any `src/rigel/native/` edit. |
| indexes rebuilt? | **no** — and none needed yet. The region builder is untouched so `partition_hash` is unchanged and every cache is valid. **W1b changes that.** |
| goldens | regenerated ONCE (at W0.5, for D6). Two more expected: **W1b** and **W7**. |

---

## 1. THE DESIGN IN ONE PAGE — read this before touching code

Rigel deconvolves each genomic object's fragments into gDNA vs RNA. The accumulator is what turns a BAM
into per-object observations. **v5 changes four things and deletes a great deal.**

```
    NODE = a genomic interval.  Owns fragments CONTAINED in it.
    EDGE = a transition.        Owns fragments CROSSING it.
             CONTIGUOUS(i, i+1)                 — a genomic point
             JUNCTION(donor, acceptor, strand)  — an annotated intron

    deposit(fragment):
        crosses 0 edges  ->  node:  count += 1,  recip += 1/L
        crosses k edges  ->  EACH:  count += 1,  recip += 1/L

    estimator at object n, component c:   rho_c = (its share of the count) / E_c(n)
```

1. **TSS/TES become partition cuts.** 59.5 % of human transcript termini currently fall *inside* a merged
   region, so today's partition cannot see them. v8 removes the merge — **six lines** (`regions.py:128-133`)
   — and nothing else about the cut set.
2. **Splice junctions become first-class directed EDGES** with their own exact divisor.
3. **Objects carry an integer COUNT, not fractional mass.** Boundary two-sidedness, fractional mass and the
   `n_cross` halving all disappear. Integer addition is associative ⇒ bit-exact at any thread count.
4. **Every object stores `(count, Σ 1/L)`** — two moments, a complete 2-component deconvolution.

**Why `Σ1/L` exists, stated correctly** (this was sharpened this session, see §4): away from a terminus the
crossing opportunity is ∝ `w`, so the observed sample is **length-biased**; `1/L` cancels exactly that bias,
giving `E[Σ1/L] = ρ_total` with no divisor and no FL model. It does **not** correct the *placement* loss near
a terminus — that is what `reach` is for. Two different corrections.

**Junction edges are FACTORS, not BP variables.** Every junction closes a cycle (~404 K on human) and loopy
BP is settled as out of scope. The message chain stays linear — `node ↔ contiguous-edge ↔ node` — and a
junction's `(flux, Σ1/L, E_J)` enters as a local measurement at its incident nodes. That is what
"graft/peel move to junction edges" means.

⛔ **Settled, do not re-derive:** no path/cell store · no loopy BP · no anchor rule · no per-node FL
histograms · `Σ(1/L)` not `ΣL` · no fractional mass · the deposit rule needs no partitioning at any node
count · D1 no taper for nascent · D2 maximal reach.

---

## 2. WHAT TO READ, IN ORDER

| # | doc | why |
|---|---|---|
| 1 | **this file** | state, roadmap, gates |
| 2 | `06_implementation_plan.md` | **F1–F16** are measured amendments to the spec; **§3** is the phase plan; **§5** records the closed decisions Q1–Q5 |
| 3 | `accumulator_ledger.md` | every arm's measured delta, appended as it landed. **W1a-Q** holds two owner challenges settled by simulation |
| 4 | `05_accumulator_v5.md` | THE SPEC. Standalone. ⚠ Read it *through* the plan's F-findings — four of its statements are amended |
| 5 | `../index/00_splice_graph_design.md` | the v8 graph. ⚠ Its **P2** gate is unsatisfiable; use **P2′** (plan F2) |
| 6 | `../calibration/P1G_SCOPE.md` §7/§8 | the W2 contract (C1/C2/C3) |
| 7 | `../calibration/ROADMAP.md` | calibration state of play; §6 item 8 is P1g |

⛔ **Do NOT read as guidance** (all bannered SUPERSEDED): `02_redesign_derivation.md`,
`03_path_accumulator.md`, `04_accumulator_v3.md`, `../calibration/PATH_MARGINALIZATION.md`,
`../calibration/archive/**`. `dag_design.md` is the owner's concept draft — context, not spec.

---

## 3. WHAT LANDED

### W0 — instruments, and one arm

* **W0.1** baseline re-recorded, **32/32 exact** at both refits ⇒ the stored reference is valid.
* **W0.2** `scripts/debug/z2.py` — THE `z2` denominator. `Var(log f_g)` was being divided into a linear
  squared error in **4 of 5 tools**; only `pass0_error_table.py` had the 2026-07-28 fix. All five now import
  `lin_var`. Any pre-fix `z2` from the other four is not comparable to a post-fix one.
* **W0.3** `pass0_oracle_bench.py` rows carry `partition` / `mass_kind` / `refit`; `--report` **refuses**
  cross-partition diffs without `--coarsen`; the writer refuses a schema-mismatched append.
* **W0.4** `TranscriptIndex.partition_hash` — blake2b-8 over the `(boundary_positions, ref_pos_offsets,
  region_types)` triple the scanner actually receives. **Computed on demand, never stored** (a hash beside
  the feathers can go stale against them). One edit in `selfsolve_diag._scan_and_truth` namespaces both the
  cache and the work dir, covering all 55 call sites. `calib_cache.load(path, index=...)` raises on mismatch.
* **W0.5** ⭐ **D6 landed** — the pooled-seam factor-2, both sites. See §5.
* **W0.6** ⛔ **CANCELLED — claim refuted by its own gate.** See §5.
* **W0.7** `Accumulator::merge_from` bound + a 6-test A9 **control**. Integer channels are bit-identical at
  any worker count; float mass channels are **not** (17/28 and 20/28 cells differ, max rel 3.7e-7) — the
  documented ~2.6 % nondeterminism, now pinned by a test. ⚠ **That float test must be INVERTED to assert
  exact equality when v5's fixed-point `uint64` recip lands. That inversion IS the §6/A9 deliverable.**

### W0b — the coarsening map (`scripts/debug/coarsen.py`)

Validated on the **real** v8→v7 relationship (1,043,881 → 752,654, all four invariants, 0.7 s). Two measured
limits, both load-bearing:

* ⚠ **27.9 % of v8 boundary slots have NO v7 counterpart** — excluded, and the excluded share is reported.
* ⭐ **The Jensen gap is 0.739×** — a coarse mwae reads 26 % better than the fine one with *no solver
  change*. This is why the standing gate scores on the **fine** population (Q4).
* Eff-lengths / densities / variances are **REFUSED**, not approximated (superadditive: `ΣE(children)/E(whole)`
  = 0.765).

### W1a — the v8 splice graph, built and validated, NOT wired — ⭐ GATE PASSED

`src/rigel/calibration/splice_graph.py` + `tests/calibration/test_splice_graph.py` (60 tests).

| gate | result |
|---|---|
| ⭐ bit-identical 32/32, both refits | `+0.000000` — nothing reads the graph, nothing moved |
| ⭐ **P3** merge v8 back == `regions.feather` | **EXACT, 752,654 rows, human scale** |
| **P2′** v7 interfaces ⊆ v8 cuts | ✅ all 286 refs |
| G1–G18 · I1–I12 (validators proven to FIRE) · P1–P5 · T-D1 byte-identical rebuilds | ✅ |
| §8 budgets | build **4.4 s** / builder allocates **0.12 GB** (v7's allocates 0.74 GB); validate **4.6 s** |

Census: **1,043,881 nodes** · median 151 bp · 15,687 of 1 bp · 1,043,595 contiguous + **404,168** junction
edges.

---

## 4. THE SEVEN THINGS A NEW SESSION MUST NOT RE-DERIVE

1. ⭐ **The v8 census is 1,043,881 nodes, not 992,068** (plan F1). The graph doc counted with
   `~is_synthetic & ~is_nrna`; the builder uses `~is_synthetic` only, and 26,475 `is_nrna &
   ~is_synthetic` rows DO cut the partition.
   ⛔ **F1's "TWO filters" is STRUCK — there is ONE filter, `~is_synthetic`** (ledger W1c). F1 gave
   the flags and reaches `~is_synthetic & ~is_nrna`; measured, all 26,475 of those rows are
   **single-exon real transcripts** (`n_exons == 1`, none a `RIGEL_NRNA_*` row), and every
   manufactured span is already `is_synthetic`. The extra clause deleted **52,104 real terminus
   positions**. Invariant **I13** now makes that class of bug impossible to reintroduce silently.
2. ⭐ **Graph-doc P2 is unsatisfiable** (F2). The merge *deletes* cuts. Use **P2′** (subset) + **P3** (merge
   reproduces v7). Both green.
3. ⭐ **`Σ1/L` is divisor-free ONLY where reach does not bind** (F5b). Simulated: `Σ1/L ÷ ρ_M` = 0.992
   mid-transcript, **0.103** at 20 bp from a TES. `Σ1/L` corrects the LENGTH BIAS; `reach` corrects the
   PLACEMENT LOSS. Neither substitutes for the other. **This matters for W5e's moment solve.**
4. ⭐ **`reach` is `reach_lo`/`reach_hi`, PER STRAND, on BOTH edge kinds** (Q1, C1–C3).
   Not `donor`/`acceptor`: `src < dst` is genomic order, so `src` is genomically left whatever the strand,
   and a NEG junction's biological donor is on the RIGHT. Measured: contiguous edges carry MORE divisor error
   than junctions (mean 0.750 vs 0.886 of `fl_mean`).
   ⚠ **reach is INERT today** — the spliced channel is down-weighted, `f_g` moves ≤1e-4. **Claim nothing
   until W5c**, then A/B it. If it loses, delete four columns.
5. ⭐ **True `f_g` genuinely rises toward a terminus (0.244 → 0.774).** gDNA dominates there from the counts
   alone. `f_g` (composition) and `ρ_M` (density) are different quantities; only `ρ_M` is broken. Relevant to
   P1g/C1, whose framing has been that the solver *wrongly* calls gDNA at termini.
6. ⭐ **Strand-coincident junctions are biologically impossible** (`GT..AG` reverse-complements to `CT..AC`);
   0 in GENCODE. The edge sort key is `(src, kind, dst, strand)` — **strand is required for a total order**,
   and `kind` stays ahead of `dst` for the CSR grouping contract. `validate_graph` **warns**; it does not
   raise, because G18 builds this input to prove the guard works.
7. ⭐ **`region_types` and the FL pools SURVIVE** (F4). §11.4's "intergenic-contained = pure gDNA" pool would
   **destroy the gDNA FL estimate under hybrid capture**, where the signal comes from exon↔intron "splash"
   reads — intergenic is not captured. Owner domain fact, plan §6.

---

## 5. THE TWO ARMS THAT CHANGED BEHAVIOUR, AND WHY ONE DIDN'T

**W0.5 — D6, the pooled-seam factor-2. LANDED.**
`gdna_boundary_len` is already `E[min(ℓ,L)]/2`, so summing two faces and dividing by their *average* read
**2ρ** (reproduced: 1.994 / 2.002 / 1.981). Two sites — `capture_eff_length.py:214` and `:318` — and they
are **one definition used inconsistently**, not two bugs: site 2's mass is imputed as `ρ_avg·s_j` so its
density is unbiased, but halving `s_j` halves the junction seam's WEIGHT in `span_full`, and fixing only
site 1 re-creates the nascent<mature inversion.
⚠ **The fixtures were why it hid**: all four stored the UN-halved `E[min(ℓ,L)]` and deposited `ρ·bl/2` —
same face mass, doubled length — exactly cancelling the code's spurious ½. Repaired, not relaxed.
Result: **bit-identical 32/32 on the bench** (it is an EM-side consumer downstream of `f_g`), goldens moved
~6e-5. ⭐ With the fix a seam's support is `gbl_r + gbl_{r+1}` → `fl_mean`, i.e. the corrected v7 arithmetic
lands exactly on v5 §10.3.

**W0.6 — the "dead" gDNA seed weight. CANCELLED, claim REFUTED.**
Pre-registered gate: assert `weight == 1.0` on 4 cfRNA payloads + 32 synthetic conditions **before**
deleting ~200 lines. Measured: **23–66 % of real-cfRNA seeds differ from 1.0**, 100 % on `gdna_none`,
1,067,456 seeds total. The reasoning held only for contained-region seeds; boundary-side and
imputed-density seeds are a different path and are the majority. **The channel is LIVE**, and
`boundary_side_eff_length` keeps a real non-geometry consumer until W5. `src/` untouched.

> **This is the methodology working.** Both arms were pre-registered with a falsification test. One passed
> and landed; one failed and was dropped for the cost of a measurement. Keep doing this.

---

## 6. THE ROADMAP FROM HERE

### ✅ W1b + W1b-clean — DONE. See the banner at the top of this file and ledger W1b / W1b-clean.

<details><summary>original plan, kept for the reasoning (its premise about the bench was wrong)</summary>

### W1b — WIRE THE v8 PARTITION

Flip `build_region_partition_arrays` to the node cut array so the scanner sees the v8 partition.

```
1. build_region_partition_arrays(index) -> use index.nodes_df when present, else region_df.
   Same 3-tuple, same ref_names iteration order. region_types = coarse_type_array(node signature).
2. RegionArrays.from_region_df -> from_nodes_df (keep the old name as a thin alias; 100 call sites).
3. INDEX_FORMAT_VERSION 7 -> 8, with the history note. The loader now REQUIRES nodes/edges.
4. Rebuild all 8 indexes. partition_hash CHANGES here ⇒ every _selfsolve_cache and calib_cache
   is invalidated by construction. That is the design working; re-scan.
5. pipeline.py:281 — replace the silent `getattr(index, "region_df", None)` skip with an explicit
   v7/v8 dispatch that RAISES. Today a v8 index would disable calibration with NO error. (plan F16)
6. CalibrationResult._check_region_array — relax the float64-only dtype gate before any integer
   count arm reaches it (plan F16).
```

**Gate — this arm is EXPECTED TO MOVE, and that is the point.**
* Record the delta in the ledger as **THE PARTITION EFFECT**. It is one of the four entangled effects and
  this is the only arm that isolates it.
* ⭐ **Record the NEW baseline here (both refits).** It becomes the reference for W2 onward *and* for W4's
  legacy-vs-legacy gate (plan F7). The pre-W1b baseline retires.
* **Regenerate goldens (2 of 3).**
* ⭐ Rely on this: **the FL pmfs do NOT move at W1b** — all 118,195 splitting v7 regions are `coarse_type
  EXON`, and the gDNA pool is intergenic+intronic. So W1b isolates the geometry change cleanly.
* ⚠ `node_chain.build_node_chain` hard-asserts `(b1−b0) == k+1` with the message "rebuild the index" — it
  fires first on a v8 partition and points at the wrong thing. Fix the message or gate it.
* ⚠ `test_capture_eff_length.py`'s module docstring locks a regression that depends on the merge being
  present; its scenario must be rewritten to keep that regression alive under v8.

</details>

### ✅ W1c — DONE. See the banner and ledger W1c.

<details><summary>original plan</summary>

### W1c — thread the structural flags to `calibrate()`

`index.edges_df` → `build_boundary_flags_array(index) -> uint16[B]` aligned to
`payload.ref_boundary_offsets` (⚠ the flag↔boundary index spaces are **off by one per reference**: v8 has
`k−1` contiguous edges for `k` nodes; the accumulator has `k+1` boundary slots — assert the offset and the
two zeroed terminals) → a new keyword on `calibrate()` → new bool fields on `NodeStatics`.
**Gate: bit-identical — nothing reads them yet.**

</details>

### ⭐ W2 — take the P1g win off the INDEX ALONE. THREE arms, order INVERTED vs P1G_SCOPE §10 *(NEXT)*

C1/C2/C3 are welded through one 10-line function, `_rho_faces` (bp_solver.py:537-546), so they cannot be one
arm (plan F9). Baseline re-recorded between each.

* **W2a — C3.** ⚠ **NOT a 4-line substitution** (plan F10). The doc's
  `rna_crosses_contiguously = exon_s(src) AND exon_s(dst) AND NOT terminus` is nearly the **complement** of
  what `accept_*` means at a junction: only **13.7 %** of junction seams have `exon_s` on both flanks, so a
  literal swap drops the mature term at ~86 % of them and fails P1G_SCOPE's own falsification test. The
  correct analogue is **per-face**: `DONOR_s or ACCEPTOR_s` on the flank carrying `exon_s`. Also needs
  `node_total_density` to return the two per-face spliced densities separately.
* **W2b — C1.** Terminus gate on the **gDNA leg only** of the reframe. ⛔ Do NOT divide `r` by `f_g(dst)` —
  that is the pin bug rebuilt. Falsification: **junction-only edges must be UNMOVED**.
* **W2c — C2.** `ω_graft` per structural class. ⚠ Watch the fit-population bias: one-seam exons never enter
  the fit and are the worse-behaved half, and terminus-flanked exons are disproportionately one-seam.

### W3 — Python reference accumulator + tests, written to target and FAILING

A1–A14, E1–E4, S1–S2. Lift `acc_proto_e.py` whole as S1 and the A3/A4 fixture (its
`A=500 B=50 C=10 D=1 E=1000` geometry already covers the 1 bp node and the encompassed middle);
`acc_proto_g.py` as S2; `acc_proto_d.py`'s enumerators as the E1–E4 oracle. Run it on a REAL cfRNA payload
through a shim before any C++ exists. Build `pass0_oracle_bench`'s v5 object mapping here.

⚠ **Freeze three contract items in the header BEFORE W4**, do not discover them: the exact rounding of
`R = round(2^32/L)` (`llround` vs banker's differ at ties, and byte-identity is undefined without it);
whether `L` is pre- or post-clip (observable in A14); and the per-object-class channel API.

### W4 — C++ accumulator + SHADOW MODE

`set_graph(...)` as a **second** method, not an overload (the one-shot guard at bam_scanner.cpp:1118-1122
throws today). ⚠ **Budget peak RSS here, not W6**: ~94 MB/copy × 9 workers for the legacy arm plus v5's
81 MB/copy ⇒ **~1.6 GB peak in shadow mode**.
⚠ **Gate restated** (plan F7): legacy-vs-legacy on ONE FIXED v8 cut array across the deposit rewrite, exact
on integer channels, at `total_threads=1` (or tolerance-bounded on floats, stated). The original "legacy
bit-identical 32/32" is unsatisfiable — the strict `!=` version gate makes a v7 index unloadable after W1b.
⚠ **Invert the W0.7 float determinism test here.**

### W5 — consumers, one observable per arm

`W5a` contained counts · `W5b` contiguous-edge flux · `W5c` junction edges + `E_J` **(and the reach A/B —
this is where reach earns its keep or gets deleted)** · `W5d` the §4 effective lengths — ⚠ structural-zero
masks BEFORE removing `_EPS` (plan F6: node divisors collapse to exactly 0.0 on 12.4 % of v8 nodes; §8's
"delete the floors" is right for edges and wrong for nodes) · `W5e` the `(count, Σ1/L)` moment into
`node_init`, **gated on real cfRNA**.

⚠ Highest-risk mechanical change in the whole rework: the `_relay`/`_transport` and
`_peel_share`/`_peel_share_scalar` twin pairs are deliberate, explicitly forbidden from merging (a measured
15.7×/op), and the face index threads through six functions. **Every face-removal edit must be hand-mirrored
into both arms with nothing enforcing it.** Re-measure on `cfrna:LBX0190`, never the toy.

### W6 — optimization (budgeted) · W7 — deletions, then goldens (3 of 3)

⚠ W6: the object count is 1.04 M nodes + 1.04 M contiguous + 404 K junction ≈ **2.49 M vs today's 1.5 M
chain nodes**. v5 §14's "may get FASTER" holds only if junction edges are factors, not ψ solves.
W7: run `ruff check src/ tests/ scripts/` and treat undefined-name failures as the authoritative delete
list. ⚠ `scratchpad/` is NOT linted — enumerate its ~40 gate scripts by hand.
Deletion census: **≈3,800 lines out, ≈1,600 in.** Plan §4 has it by area, with a **⛔ DO NOT DELETE** list of
seven items that were each proposed and are each wrong.

---

## 7. STANDING METHODOLOGY — non-negotiable

* **No magic numbers.** Pause and discuss before ANY new constant or heuristic.
* **One thing varied per arm**, with a **pre-registered falsification test**. W0.6 is the proof this pays.
* **Re-record the baseline from the current tree in the same session**, both refits. If HEAD-vs-baseline is
  not 32/32, the BASELINE is broken, not your change.
* **Append to `accumulator_ledger.md` at EVERY arm.** Four entangled effects (partition, deposit rule,
  junction divisor, seam factor-2) are only attributable if each delta is recorded as it lands.
* **Held-fixed `z2` must not regress**, scored on the **fine** population (Q4).
* **ruff + `ruff format` + full suite on every arm.** Goldens at W1b and W7 only.
* ⚠ **Profile and gate on REAL cfRNA, never the 10 Mb toy**: it overstates contained efficiency 4.6×,
  understates multi-crossing 3.8×, has no alternative TSS/TES, and is Poisson by construction.
* **The owner drives commits.** Do not commit unless asked.

---

## 8. ENVIRONMENT AND TOOLS

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
# ALWAYS OMP_NUM_THREADS=1.  After any src/rigel/native/ edit:
pip install --no-build-isolation -e .
```

| tool | use |
|---|---|
| `scripts/debug/pass0_oracle_bench.py --arm NAME` | THE A/B. `P0_REFIT=0\|3`, `P0_BENCH_OUT=...`. `--report --vs A --new B` |
| `scripts/debug/z2.py` | THE `z2` denominator — `lin_var`, `z2`. Import it; never re-implement |
| `scripts/debug/calib_cache.py` | cfRNA fast loop. `load(path, index=...)` verifies `partition_hash` |
| `scripts/debug/graph_human_gate.py` | the human-scale graph gate (I1–I12 incl. I3b, census, §8 budgets) |
| `scratchpad/reach_worth_it.py` | the reach / `Σ1/L` simulation |
| suites | `~/Downloads/rigel_runs/ambig_dense_10mb` (32 cond, oracle) · `.../cfrna/_calib_cache` (4 real) · `.../refs/rigel_index_v7` (human) |

⚠ **`/tmp/rigel_selfsolve` and `_selfsolve_cache` are now namespaced by `partition_hash`** — two sessions no
longer corrupt each other, and a stale payload cannot silently load.

---

## 9. OPEN QUESTIONS FOR THE OWNER

1. **The third golden regeneration.** Decision Q3 was "regenerate twice (W1b, W7)", taken when D6 sat at
   W5f. D6 moved to W0.5 and moves 21 goldens by ~6e-5, so they were regenerated there too. Revert to
   carrying them red, or accept three?
2. **`partition_hash` is computed, not stored in `manifest.json`** as the plan's W0.4 said. Deliberate — a
   stored hash can go stale against the feathers beside it — but it is a deviation.
3. **`reach` earns its keep at W5c or gets deleted.** Currently inert. Four columns, ~17 MB.
4. **Should nRNA spans stay in the event set?** They must for P3 to hold. If the owner would rather exclude
   them, P3 becomes "reproduces v7 modulo the nRNA cuts" and the partition changes — a separate decision.
