# Accumulator v5 — the arm ledger

**Append one row per arm, as it lands. Never retroactively.**

The rework moves four entangled effects — the **v8 partition**, the **deposit rule**, the **junction
divisor**, and the **seam factor-2** — and each one's delta is only attributable if it is recorded at the
moment it lands, against a baseline recorded from the same tree in the same session. This file is that
record. Plan: `06_implementation_plan.md`. Spec: `05_accumulator_v5.md`.

Standing gates on every arm: ruff · full suite · one thing varied · a pre-registered falsification test ·
held-fixed `z2` must not regress (scored on the **fine** population against the W1b baseline, per plan Q4) ·
goldens regenerated only at **W1b** and **W7** (plan Q3).

---

## W0.1 — BASELINE (2026-07-29)

Recorded from the working tree at `3c293038` with `src/` and `tests/` clean, in one session, both refits.
Tool: `scripts/debug/pass0_oracle_bench.py --arm head`, suite `ambig_dense_10mb`, 32 conditions,
`OMP_NUM_THREADS=1`.

| | mass-weighted mwae | recorded reference | delta |
|---|---|---|---|
| refit = 0 (pass-0) | **0.079005** | 0.079005 | `-0.000000` |
| refit = 3 (production) | **0.046675** | 0.046675 | `+0.000000` |

**32/32 exact on both refits** ⇒ the stored reference is valid and HEAD reproduces it. Every W0–W7 arm has a
comparator. Provenance: `partition=v7`, `mass_kind=mass_frac`.

Artifacts: `/tmp/w0_baseline_r0.tsv`, `/tmp/w0_baseline_r3.tsv`.

⚠ This baseline is **retired at W1b**, when the partition changes. W1b records its successor, which becomes
the reference for W2 onward and for W4's legacy-vs-legacy gate (plan F7).

---

## W0.2 — z2 denominator unified (instrument, not an arm)

`Var(log f_g)` was being divided into a linear squared error in four of five tools — the mixed-scale bug
fixed in `pass0_error_table.py` on 2026-07-28 and nowhere else. `scripts/debug/z2.py` is now the single
definition (`lin_var`, `z2`); the five consumers import it.

**No production code touched.** Effect is on the *reported* statistic only: raw log-space denominators
under-state `z2` (the suite total read 0.046 when it is ≈1), so any pre-2026-07-28 `z2` from
`ablate_report.py` / `v_p4_paired_z2.py` / `ho_subz2.py` / `v_ho_subcheck.py` is not comparable to a
post-fix one.

---

## W0.5 — ARM: D6, the pooled-seam factor-2 (v5 §10.4) — PRE-REGISTRATION

**One thing varied:** the seam SUPPORT definition in `capture_eff_length.py`. Nothing else. Pure v7
arithmetic — independent of the index, the accumulator and the solver, which is why it lands first.

**The defect.** `gdna_boundary_len` is already `E[min(ℓ,R)]/2` (`calibrate.py:226` →
`effective_length.boundary_side_eff_length`). `_pooled_seam_arrays` SUMS the two faces' mass but divides by
their AVERAGE, `0.5·(gbl_r + gbl_{r+1}) = (E[min_r]+E[min_{r+1}])/4`, against an expectation of
`ρ·(E[min_r]+E[min_{r+1}])/2` ⇒ **the pooled-seam density reads 2ρ.** The correct divisor is the SUM,
which is what `_gdna_region_node_arrays`'s own docstring already says (`S_s = ½·(E[min_r]+E[min_{r+1}])`,
and since `gbl = E[min]/2` that IS `gbl_r + gbl_{r+1}`). The prose beside the code says "AVERAGE" and the
code followed the prose.

⭐ Reproduced this session (`scratchpad/acc_seam_check.py`, uniform truth, 395 seams):

| region length | 2000 | 500 | 200 | 100 | 50 |
|---|---|---|---|---|---|
| shipped (averaged) ÷ truth | **1.994** | **2.002** | **1.981** | 1.972 | 1.847 |
| summed ÷ truth | 0.997 | 1.001 | 0.990 | 0.986 | 0.923 |

**TWO sites, and they are not two independent bugs — they are one definition used inconsistently.**

* **Site 1** `capture_eff_length.py:214` — `seam_S = 0.5·(side_len[r] + side_len[r+1])`. Here `seam_m` is a
  MEASURED mass, so halving the divisor genuinely doubles the density. **A real bias.**
* **Site 2** `capture_eff_length.py:318` — `s_j = 0.5·(side_len[jl] + side_len[jr])`, the splice-junction
  seam. ⚠ Here the mass is IMPUTED as `m_j = ρ_avg·s_j`, so `m_j/s_j = ρ_avg` **regardless of `s_j`'s
  scale** — the junction seam's *density* is unbiased either way. What the halving corrupts is its
  **WEIGHT** in `span_full`. Fixing only site 1 therefore doubles every genomic seam's weight while leaving
  junction seams at half — and a spliced mRNA (few interior genomic seams, many junction seams) loses
  footprint against its nascent parent, producing the physically impossible inversion. Confirmed: patching
  site 1 alone fails `test_no_nascent_mature_inversion_under_capture` with *"INVERSION: nascent eff_em
  310.17 < mature eff_em 412.48"*.
  ⚠ The `0.5·(rho_l + rho_r)` on the same line is a genuine AVERAGE OF DENSITIES (the junction's imputed
  density is the mean of its two flanks) and **must not be touched**.

**Predictions, pre-registered:**

1. **The pooled-seam density moves by exactly 2×.** Expected and must NOT be attributed to any accumulator
   or partition change later.
2. **`factor_t` moves only where the clip does not bind.** Under uniform gDNA `min(2·S, S) = S` ⇒ factor 1
   either way, so **capture-OFF stays bit-identical**. Under capture, seams whose true density lies in
   `(ρ_ref/2, ρ_ref)` currently clip to no contraction and should now contribute some. So: capture-OFF
   flat, capture-ON/verystrong move.
3. **The suite metric may not move at all**, because `pass0_oracle_bench` scores calibration `f_g`, and
   `capture_eff_length` is an EM-side consumer downstream of it. A flat bench is a PASS, not a null result.
4. ⭐ **FALSIFICATION TEST.** The reason no test caught this is that the clip rescues the uniform case, so
   the "factor-1-under-uniform bedrock" invariant passes with the bug present. A new test must construct a
   seam whose true density sits strictly inside `(ρ_ref/2, ρ_ref)` with `ρ_ref` genuinely detected, and
   assert the contraction is non-trivial. **If that test passes against UNPATCHED code, the fix is not the
   fix and this arm is wrong.**

**Gates:** ruff · full suite · `acc_seam_check.py` promoted to a unit test asserting 1.99–2.00× pre-fix and
0.98–1.01× post-fix · the three fixture builders repaired rather than the assertions relaxed.

### W0.5 — RESULT (2026-07-29): landed, all four predictions held

**Falsification tests written FIRST and verified failing against unpatched `3c293038`:**

| test | unpatched | patched |
|---|---|---|
| `test_pooled_seam_density_recovers_rho` | **2.000×** ρ (0.074 vs 0.037) | 1.000× to 1e-12 ✅ |
| `test_seam_inside_the_clip_band_still_contracts` | reads **1.4–1.7×** ρ_ref ⇒ clips ⇒ no contraction | < 1 ⇒ contracts ✅ |

**Code changed — two lines, one definition:**
`capture_eff_length.py` `_pooled_seam_arrays`: `0.5·(side_len[r] + side_len[r+1])` → `side_len[r] + side_len[r+1]`;
and the junction seam `s_j`: same. The `0.5·(rho_l + rho_r)` average-of-densities on the same line is
untouched, as pre-registered.

**Fixtures repaired, not assertions relaxed.** All four builders (`test_capture_eff_length._uniform_field_cal`
and `._field_cal`; `test_priors._calibration`, `._blind_seam_cal`, `._global_bimodal_cal`, and two inline
fixtures) stored `gdna_boundary_len` as the UN-halved `E[min(ℓ,L)]` while depositing `ρ·bl/2` per face — the
same face mass, but a **doubled length**, which exactly cancelled the code's spurious ½. That is why the bug
was invisible to 29 tests. They now store the production semantic `E[min(ℓ,L)]/2` and deposit `ρ·gbl` per
face. ⭐ With that, a seam's support is `gbl_r + gbl_{r+1}` → `fl_mean` for long flanks — **the corrected v7
arithmetic converges exactly on v5 §10.3's `ρ = E[flux]/fl_mean`**, which is independent evidence the fix is
the right one.

**Predictions vs outcome:**

| # | prediction | outcome |
|---|---|---|
| 1 | pooled-seam density moves exactly 2× | ✅ 2.000× to 1e-12 |
| 2 | capture-OFF bit-identical (the clip binds), capture-ON moves | ✅ consistent with 3 |
| 3 | the **calibration** bench does not move — this is an EM-side consumer downstream of `f_g` | ✅ **bit-identical 32/32 at BOTH refits.** r0 `0.079005 → 0.079005`, r3 `0.046675 → 0.046675`, `delta +0.000000`, 0 better / 0 worse / 32 flat, and 32/32 conditions equal to the last recorded digit |
| 4 | falsification: the new tests must fail unpatched | ✅ both did |

**Suite:** 1215 pass, ruff clean on `src/ tests/ scripts/`.

⚠ **GOLDENS MOVED, and this is a deviation from decision Q3 that needs owner awareness.** Q3 chose to
regenerate at W1b and W7 only — decided when D6 sat at W5f. D6 now lands at W0.5, and it moves 21 golden
scenarios by **≈ 6e-5 relative** (e.g. `wide_intron/transcript.count` 311.9371 → 311.9554) because the EM
effective lengths genuinely change. I regenerated them here rather than carry 21 red goldens into W1a, whose
entire gate is *"bit-identical 32/32 **and** full suite green"*. **That makes three regenerations: W0.5, W1b,
W7.** 147 golden files, 218 insertions / 218 deletions.

---

## W0.6 — ⛔ CANCELLED: the gDNA seed weight is **NOT** dead. Claim REFUTED by its own gate.

**Claim** (review critic, plan §4): `count_gdna_frac` evaluates to exactly 1.0 on every fitted seed except
two narrow masks, so the gDNA strand-fit seed-weight channel is arithmetically dead and ~200 lines across
`density_model.py` / `gdna_strand.py` / `strand_deconv.py` can be deleted — which would also have removed
`boundary_side_eff_length` from the strand-fit call graph.

**Gate** (pre-registered): assert `weight[valid] == 1.0` on the four cached cfRNA payloads and all 32
synthetic conditions BEFORE deleting anything. Harness: `calibrate()`'s prefix only (substrate →
eff-lengths → κ → `node_gdna_density` → `fit_gdna_strand_from_substrate`), never the BP sweep.

**Measured — the gate FAILS:**

| population | n seeds | `weight ≠ 1.0` | min |
|---|---|---|---|
| cfRNA LBX0190 | 450,521 | **297,565 (66.05 %)** | 0.000 |
| cfRNA LBX0588 | 514,110 | **240,468 (46.77 %)** | 0.000 |
| cfRNA MO_3021 | 864,467 | **276,651 (32.00 %)** | 0.000 |
| cfRNA vcap | 1,053,265 | **242,903 (23.06 %)** | 0.000 |
| synthetic `gdna_none` × 4 | 726 each | **726 (100 %)** | 0.000 |
| synthetic, remaining 28 | 1.4–3.2 K | 2.22 – 44.07 % | 0.000 |
| **total** | | **1,067,456 seeds with weight ≠ 1.0** | |

**Why the reasoning broke.** It is sound only on the *contained-region* seed path: there
`density[own] = contained_gdna[own]/region_eff_len[own]` and `contained_gdna == contained_mass` (the
substrate passes `region_contained` as both the count and the mass, `substrate.py:102`), so
`count_gdna_frac[own] ≡ 1`. But that is not the seed population. **Boundary-side seeds**
(`strand_deconv.boundary_side_seeds`) and every region outside `own` — which takes an *imputed* density
from `runfill_bidirectional` or the global baseline — go through a different path and land anywhere in
`[0, 1]`. On a zero-gDNA library the weight is 0 everywhere, which is the channel working correctly.

**Verdict: the channel is LIVE. No deletion. `boundary_side_eff_length` keeps a real non-geometry
consumer**, so it cannot be retired until that consumer moves in W5. Plan §4's "≈200 lines" and the
"removes `boundary_side_eff_length`'s last non-geometry consumer" note are both struck.

*Nothing was changed in `src/`. Cost: one measurement. This is the gate paying for itself.*

---

## W1a-Q — TWO OWNER CHALLENGES, SETTLED BY SIMULATION (2026-07-29)

### Q-A — "does `Σ1/L` already compensate for the terminus taper?" → **NO, it is ~10× short**

**Owner's argument.** Each fragment deposits `1/L` at deposit time — a measured single length, not a pmf.
Near a terminus only short fragments fit, so each deposits a larger `1/L`, and the density self-corrects.
If so, `reach` is over-engineering.

**Measured** (`scratchpad/reach_worth_it.py`, true physics: mature fragments must lie entirely inside the
transcript; gDNA unconstrained; both densities uniform per bp):

| d(TES) | n obs/rep | mean OBSERVED L | mean 1/L | `Σ1/L` ÷ ρ_M | ρ̂_M/ρ_M no reach | with reach |
|---|---|---|---|---|---|---|
| 1500 | 6.00 | **214.2** | 0.00496 | **0.992** | 1.039 | 1.045 |
| 200 | 5.36 | 207.4 | 0.00511 | 0.914 | 0.920 | 1.025 |
| 100 | 2.96 | 202.0 | 0.00531 | 0.524 | 0.497 | 0.998 |
| 50 | 1.46 | 201.4 | 0.00534 | 0.260 | 0.237 | 0.950 |
| 20 | 0.59 | **202.2** | 0.00528 | **0.103** | 0.098 | 0.983 |

**The compensation is real but is worth +6.5 % against a −90 % count loss.** Mean observed length moves
214.2 → 202.2 (−5.6 %) while the crossing count falls 6.00 → 0.59 (−90 %).

⭐ **Why the fragments do NOT get shorter.** At a point with `R_a` bases downstream, a length-`w` fragment
has `trapezoid(w) = min(w−1, R_d, R_a, R_d+R_a−w+1)` placements — and for `R_a = 50`, any `w > 51` gives
**exactly 50, independent of `w`**. A 200 bp fragment still crosses a point 50 bp from the TES; it just sits
asymmetrically. The terminus removes **placements**, not **lengths**.

⭐ **And that reveals what `1/L` is actually for.** Away from a terminus the opportunity is ∝ `w`, so the
crossing sample is **length-biased**: predicted mean `E[w²]/E[w] = 212.5`, measured **214.2**. `1/L` cancels
exactly that, giving `Σ1/L = 0.992 ρ_M`. Near a terminus the opportunity saturates at `R_a` and stops
depending on `w`, the length bias **disappears** (observed mean returns to the plain FL mean 200 — measured
202.2, converging down to it and never below), and `1/L` has nothing left to cancel, so it over-divides.

> **`Σ1/L` handles the LENGTH BIAS; `reach` handles the PLACEMENT LOSS. They are different corrections and
> neither substitutes for the other.**

⚠ **SPEC AMENDMENT — v5 §3.3's "divisor-free" claim for `Σ1/L` is true only where reach does not bind.**
Where it does, the correct divisor is `E_f[trapezoid(w)/w]`, not 1. v5 leans on that identity for `ρ_total`,
so this stands regardless of what is decided about reach.

**Where the error lives** (real human annotation, `E_reach/fl_mean`):

| edge class | live (mature crosses) | mean | off > 10 % | off > 50 % |
|---|---|---|---|---|
| JUNCTION | 404,168 | 0.886 | 23.8 % | 9.6 % |
| **CONTIGUOUS** | 391,796 | **0.750** | 36.7 % | **26.0 %** |

Contiguous edges carry MORE divisor error than junctions — a contiguous edge where mature crosses is often
inside a terminal exon, whereas a junction has exons on both sides by construction. This is the measured
justification for decision Q1 (reaches on both edge kinds).

⭐ **The owner's physical intuition is CONFIRMED, and it is free.** True `f_g` rises **0.244 → 0.774**
approaching a TES: gDNA genuinely dominates there because mature fragments cannot fit. That is in the counts
and no estimator arranges it. **`f_g` (composition) and `ρ_M` (density) are different quantities — only
`ρ_M` is broken, and reach fixes only `ρ_M`.**

⚠ **HONEST COUNTERWEIGHT, and it scopes the work:** `effective_length.py`'s own docstring records that
today's mature-divisor error (2–199× on ~26 % of junction faces) moves `f_g` by **≤ 1e-4**, because the
spliced channel is heavily down-weighted. Reach corrects the same class of quantity, so **it is inert until
v5 makes junction edges a primary observable (W5c)**. Decision: keep the columns (index-time, ~17 MB,
already built and tested), claim nothing until W5c, and A/B it there against the oracle. If it loses, four
columns are deleted and nothing else changes.

### Q-B — strand-coincident junctions: the owner is right on the biology; the guard STAYS, and now warns

Splice motifs are **non-palindromic** — a `GT..AG` intron reverse-complements to `CT..AC` — so the same
interval cannot be a valid intron on both strands. **Biologically impossible; measured ZERO in GENCODE.**
A strand-coincident junction therefore means the ANNOTATION is wrong (a simulator, or a hand-filled strand
column), not that a rare case occurred.

* The sort key keeps `strand` — `(src, kind, dst, strand)`. (The owner suggested `(src, dst, strand, kind)`;
  same guard, but `kind` is kept ahead of `dst` so a node's out-edges stay grouped by kind, which is the
  design doc's CSR contract.)
* `validate_graph` now **warns** (`RuntimeWarning`) with the count and the reason.
* ⛔ **Not a hard failure**, deliberately: G18 constructs exactly this input to prove the guard works, and
  raising would make the guard untestable. Two tests: the warning fires and the graph is still correct
  (two distinct edges); and a normal annotation warns not at all. Suite is warning-clean under
  `-W error::RuntimeWarning`.

---

## W1a — THE v8 SPLICE GRAPH, built and validated, NOT wired — ⭐ GATE PASSED

`src/rigel/calibration/splice_graph.py` (+ `tests/calibration/test_splice_graph.py`, 58 tests).

**⭐ GATE: bit-identical 32/32 at BOTH refits.** r0 `0.079005 → 0.079005`, r3 `0.046675 → 0.046675`,
`delta +0.000000`. Nothing reads the graph, so nothing may move — and nothing did.
Suite **1279 pass** (1221 + 58), ruff + format clean.

| gate | result |
|---|---|
| **P3** merge v8's equal-signature nodes == `regions.feather` | ⭐ **EXACT, 752,654 rows, at human scale** |
| **P2′** every v7 interface is a v8 cut | ✅ all 286 references |
| G1–G18 toy matrix · I1–I12 (validators proven to FIRE) · P1–P5 | ✅ 58 tests |
| **T-D1** byte-identical rebuilds of `nodes.feather` / `edges.feather` | ✅ |
| §8 build budget < 10 s / < 1 GB | ✅ **4.5 s**, builder allocates **0.12 GB** |
| §8 `validate_graph` budget < 5 s (incl. the I11 walk) | ✅ **4.7 s** |
| §3.4 census | 1,043,881 nodes · median 151 bp · 15,687 of 1 bp · 1,043,595 contiguous + **404,168** junction edges |

⭐ **The v8 builder is 6× leaner than the v7 one it replaces**: +0.12 GB vs +0.74 GB head-to-head on the
same fixture, while producing 38.7 % more nodes *plus* 1.45 M edges with flags and reaches. The dict-per-row
sweep was the memory cost, and it is gone.

### ⭐ Two real bugs the matrix caught, both invisible on real data

1. **The design doc's edge sort key `(src, kind, dst)` is NOT a total order.** Two strand-coincident
   junctions (G18) share all three and differ only in `strand`, so the order is ambiguous and the
   uniqueness check reads them as a duplicate. **GENCODE contains zero such junctions**, so no amount of
   real-data testing would ever have found it. Fixed by sorting on `(src, kind, dst, strand)` — a
   refinement of the documented order, so the "out-edges are contiguous ⇒ CSR is one searchsorted"
   contract is unchanged.
2. `validate_graph` indexed an empty junction key array when a reference had introns but the graph had no
   junction edges — i.e. exactly the corruption I11 exists to report. Guarded.

### Deviations from the plan, and why

* ⭐ **`INDEX_FORMAT_VERSION` stays at 7** (plan C4). The loader gate is a strict `!=`, so bumping here
  would make all 8 existing indexes unloadable and W1a's own bit-identity gate unrunnable until every one
  was rebuilt. `nodes.feather`/`edges.feather` are written **additively** and loaded **optionally**
  (`index.nodes_df is None` on an older index, and it stays fully usable — asserted by a test). W1b bumps,
  when the loader starts *requiring* the graph. This is graph-doc §9's own recommendation (b).
* **`reach_lo`/`reach_hi`, not `reach_donor`/`reach_acceptor`** (plan C1). `src < dst` is genomic order, so
  `src` is genomically left whatever the strand — but a NEG-strand junction's biological donor is on the
  RIGHT. The divisor is symmetric so no number changes, but the spec's names mislabel ~half of 404,168
  junctions.
* **Reach is PER STRAND — 4 int32 columns, not 2** (plan C2). The mature-crossing gate is per strand, and
  at a coincident junction the two strands have genuinely different reaches (pinned by a test where a
  strand-agnostic maximum would over-state the `+` junction's downstream reach **10-fold**).
* **Reaches on CONTIGUOUS edges too**, per decision Q1 — that is where the taper near a TES bites.
* `build_index_artifacts` now returns 5 frames, not 3 (one test caller updated).
* `tests/conftest.py::build_test_index` gained an optional `refs=` for G17; every existing caller is
  untouched.

⚠ **Not yet done, deliberately**: the CSR accessors (`ref_node_offsets`, `out_edge_offsets`, the reverse
CSR) are derived-at-load per graph doc §4 and have no consumer until W1c/W2 — building them now would be
untested plumbing. No index has been rebuilt: since the region builder is untouched, `partition_hash` is
unchanged and every cache stays valid.

---

## W0b — THE COARSENING MAP (trust instrument, not an arm) — `scripts/debug/coarsen.py`

Without a v8→v7 map the standing gate *"held-fixed z2 must not regress"* and every mwae A/B are
**undefined** from W1 onward, because the object set itself moves. Built before W1, as the plan requires.

**Validated on the REAL v8→v7 relationship, not just a synthetic one.** The true v8 cut set was built from
the annotation and mapped back onto the shipped `regions.feather`, at full human scale, in **0.7 s**:

| | |
|---|---|
| v8 nodes → v7 regions | **1,043,881 → 752,654** (1.387×) |
| **I-a** every v7 interface IS a v8 cut (**P2′**) | ✅ asserted, all 286 references |
| **I-b** every v8 node lies inside exactly one v7 region | ✅ |
| **I-c** surjective — every v7 region covered | ✅ |
| **I-d** v8 nodes tile each v7 region exactly (bp) | ✅ |
| v8 nodes per v7 region | mean 1.39, median 1, **max 115** |
| v7 regions left UNSPLIT | **634,459 (84.3 %)** |

⭐ This is an independent second confirmation of plan finding **F2**: v7 interfaces are a strict SUBSET of
v8 cuts, so the graph doc's P2 (*"the cut sets are EQUAL"*) is unsatisfiable and P2′ is the right gate.

**Self-test** (`--self-test`, run on both the toy and the human index) proves four things on a synthetic
midpoint refinement: the invariants; an EXACT round-trip when sub-nodes are homogeneous; that a
non-refinement is REJECTED rather than silently mapped; and that the refused kinds actually refuse.

### ⛔ The honest scope — what this instrument CANNOT do

* **Effective lengths, densities and variances are REFUSED**, not approximated. Eff-lengths are
  superadditive (`Σ E(children)/E(whole)` = 0.765 measured, 0.092 for 305 → 145+160), and precision-pooling
  a variance assumes independence, which BP violates by construction. `coarsen(kind=...)` raises with the
  reason rather than returning a plausible wrong number.
* ⚠ **27.9 % of v8 boundary slots (291,227) have NO v7 counterpart** — they are the cuts the v7 merge
  deleted. Their crossing fragments were never a v7 object, so folding them into the coarse region would
  invent mass v7 never had. They are EXCLUDED, and `dropped_interface_frac()` reports the excluded share so
  the exclusion is never silent. **Better than a quarter of the v8 edge population cannot participate in a
  coarse comparison at all** — a real limit on this instrument, stated rather than buried.
* ⭐ **The JENSEN GAP is real and large.** Measured on the human index with heterogeneous sub-nodes: mwae
  fine 0.23135 vs coarse 0.17086 — **the coarse score is 0.739× the fine one purely from averaging**, with
  no change in the solver whatsoever. Comparing a fine mwae to a coarse one would read as a 26 %
  improvement that does not exist. This is exactly why the standing gate is scored on the fine population
  (decision **Q4**) and why `pass0_oracle_bench --report` refuses cross-partition diffs (W0.3).

---

## W0.4 — cache provenance: `TranscriptIndex.partition_hash` (instrument, not an arm)

A cached scan payload is keyed to the exact partition it was scanned against, and **nothing checked
that**. `_selfsolve_cache` was keyed by condition name alone; `calib_cache` stored `index_dir` and never
read it back — and it pickles `region_arrays` *alongside* the payload, so the two stayed mutually
consistent while both went stale against the index on disk. Every W5e / W6 measurement after an index
rebuild would have run v7 geometry while the operator believed it was on v8.

`TranscriptIndex.partition_hash` — a blake2b-8 over the `(boundary_positions, ref_pos_offsets,
region_types)` triple `BamScanner.set_regions` actually receives, plus the reference lengths. Two indexes
hash equal **iff** a scan against them yields an identically-keyed payload, which is exactly the
reusability condition.

⭐ **Computed on demand, never stored** — a deviation from the plan's "index_hash into manifest.json", and
a deliberate one: a derived value written beside the feathers can go stale against them, which is the
failure mode this whole item exists to prevent. It also means **every existing index gains a valid cache
key with no rebuild**, and the key keeps working unchanged when the v8 graph replaces the region partition
behind `build_region_partition_arrays`. Cost: 39 ms on the human index, 0 ms on the toy.

Wired at the one place that matters — `selfsolve_diag._scan_and_truth` namespaces **both** the cache dir
and the split-BAM work dir, which covers all 55 call sites across `scripts/debug` without touching one of
them. `calib_cache.load(path, index=...)` now raises on mismatch or on an unrecorded hash.

Verified: the 32 existing `ambig_dense_10mb` payloads migrated into `_selfsolve_cache/4006135f7b855240/`,
and the bench re-ran in **12 s** (cache hit, not a rescan) reproducing `0.079005` exactly. The four cfRNA
`_calib_cache` pickles report `<unrecorded>` and will now refuse to load against an index — by design.

---

## W0.7 — `Accumulator::merge_from` bound + the A9 determinism CONTROL (instrument, not an arm)

`merge_from` was defined in C++ (`accumulator.cpp:234`) and called only from `BamScanner::scan` — it was
**never bound to Python**, so v5's starred test A9 ("N workers in any order → bit-identical") had zero
infrastructure and the v5 §6 claim had no *control*: nothing recorded what today's accumulator does.

Bound through nanobind + the `_accumulator.py` façade, plus `TestWorkerMergeDeterminism` (6 tests):

| | result |
|---|---|
| K workers + merge == one accumulator | ✅ exact on integer channels |
| integer channels at `n_workers ∈ {1,2,4,8}`, shuffled shardings | ✅ **bit-identical** — integer addition is associative |
| float mass channels across shardings | ⚠ **NOT bit-identical**: 17/28 and 20/28 cells differ, max rel **3.7e-7** |

⭐ That last row is the documented ~2.6 % cross-process nondeterminism **caught at its source**, and it is
now pinned by a test. `mass_left`/`mass_right` accumulate `float32 += share` over a data-dependent worker
partition, and float addition is not associative.

⚠ **The float test must be INVERTED to assert exact equality when v5's fixed-point `uint64` recip lands.**
That inversion *is* the §6 / A9 deliverable — if it cannot be inverted, the claim that v5 removes the
nondeterminism is false.

---

## W0.3 — bench provenance (instrument, not an arm)

`pass0_oracle_bench.py` rows now carry `partition` (v7/v8), `mass_kind` and `refit`. `--report` **refuses**
to diff arms recorded on different partitions unless `--coarsen` is passed, and the writer refuses to append
a different column set to an existing TSV. Rationale: plan F8 — the v8 partition refines v7 and effective
lengths are superadditive (`Σ E(children)/E(whole)` = 0.765), so a raw per-object mwae diff across the
partition change is not a measurement of anything.

---

## W1b — ARM: WIRE THE SCANNER TO THE v8 PARTITION (2026-07-29)

**One thing varied:** which partition `build_region_partition_arrays` hands to
`BamScanner.set_regions` — `index.region_df` (v7, merged signatures) → `index.nodes_df` (the v8
splice graph). Everything downstream follows through one accessor.

### ⛔ THE PREMISE WAS WRONG, AND MEASURING IT FIRST IS WHY THIS ARM IS CLEAN

The plan and HANDOFF called W1b *"the first arm expected to MOVE numbers"*, whose *"whole point is to
isolate THE PARTITION EFFECT"*. **On the 32-condition bench there is no partition effect to
isolate.** Measured before writing any code, by building the v8 graph on every index and comparing
the actual `set_regions` triple:

| index | v7 regions | v8 nodes | Δ | scanner arrays identical |
|---|---|---|---|---|
| **`ambig_dense_10mb` (THE bench)** | 1,698 | **1,698** | **+0.0 %** | ⭐ **True** |
| `gdna_benchmark_5mb` / `gdna_shortfl` / `nascent_5mb` | 1,082 | 1,226 | +13.3 % | False |
| `efflen_binding_sweep` | 441 | 541 | +22.7 % | False |
| `quick_3to1_5mb` | 1,221 | 1,526 | +25.0 % | False |
| human (`refs/rigel_index_v7`, `rigel_runs/rigel_index`) | 752,654 | 1,043,881 | +38.7 % | False |

Stronger than the table: on `ambig_dense_10mb` the v8 `nodes_df` is **row-for-row equal** to the v7
`region_df` on `(ref_name, start, end, length, signature)`. The merge never fires there — the suite
has no alternative TSS/TES, which is plan §15's own warning about the toy, now measured on exactly
this axis. ⭐ **The gate therefore inverts from "expected to move" to "must not move", which is the
stronger gate**, and the partition effect has to be measured on a suite that actually splits.

### PRE-REGISTERED PREDICTIONS vs OUTCOME

| # | prediction | outcome |
|---|---|---|
| P1 | 32-cond bench **bit-identical**, both refits | ✅ **32/32 exact on all 7 columns.** r0 `0.079005 → 0.079005`, r3 `0.046675 → 0.046675`, `+0.000000` |
| P2 | `partition_hash` unchanged ⇒ cache hit, no re-scan | ✅ `4006135f7b855240` → `4006135f7b855240`; the bench re-ran off cache |
| P3 | goldens do not move (AST census: 0 merges in all 20 scenarios) | ✅ **21/21 pass unmodified — NO regeneration.** Decision Q3's "regenerate at W1b" does not apply |
| P4 | a v7 index becomes unloadable, with an actionable message | ✅ and it fired for real mid-session, killing a v7 bench pass (see below) |
| P5 | geometry aligns 1:1 with the payload on every rebuilt index | ✅ no `CalibrationSubstrateError` anywhere |

⭐ **FALSIFICATION TEST, written first and verified failing against `3c293038`**
(`tests/calibration/test_partition_wiring.py`, 7 tests). Since the bench cannot see this change, the
arm needs its own proof it did anything. The case is **G4 — an alternative TSS strictly interior to
another transcript's exon**, which both flanks make `exon_pos`, so the v7 merge deletes the cut.
**6 of 7 failed unpatched, each for the right reason**; the 7th asserts the pre-state and passes by
design. The sharpest one: two annotations whose **v7 partitions are identical** and whose v8
partitions differ by one node hashed to the *same* `partition_hash` (`645a0dd8aa560236`) — a payload
scanned against one would have loaded silently against the other. That is the stale-cache failure
`partition_hash` exists to prevent, and pre-W1b it was live.

### ⭐ THE SCAN-COST ACCEPTANCE NUMBER (graph doc §8) — REAL cfRNA, +38.7 % partition

§8 makes this a required measurement, not an afterthought, and says to take it on the real
annotation. Exact A/B — same BAM, same process, arms alternating, **only the cut array differs**
(`scratchpad/w1b_scan_cost.py` patches the one function `_wire_calibration_regions` consumes):

| sample | v7 | v8 | Δ |
|---|---|---|---|
| cfRNA LBX0190 (37 MB) | 0.38 s | 0.42 s | **+9.7 %** |
| cfRNA MO_3021 (228 MB) | 1.38 s | 1.47 s | **+6.9 %** |

Cut array 6.0 MB → 8.4 MB. §8 asked for *"no measurable regression"*: it **is** measurable, and it
is ~7–10 %, not the cache cliff §8 feared. Honest input to D5/W6, not a pass.

### WHAT CHANGED IN `src/`

* `build_region_partition_arrays(index)` reads `index.nodes_df`. ⚠ **Deviation from the handoff's
  step 1 ("else `region_df`"): there is NO fallback.** After the version bump the graph is always
  present on a loadable index, so the fallback is reachable only when something is wrong, and it
  would calibrate on the merged partition while the caller believed otherwise.
* ⭐ **`RegionArrays.from_index(index)`** — the single accessor, reading the same frame. ⚠
  **Deviation from step 2** (rename `from_region_df` → `from_nodes_df` + alias): a rename touches 229
  sites and does not close the real hazard, which is *which frame is passed*. `from_index` does.
* `INDEX_FORMAT_VERSION` 7 → **8**, history note; the loader now **requires** nodes/edges.
* `pipeline.py:281`'s silent `getattr(index, "region_df", None)` skip → explicit dispatch: **missing
  graph raises**, empty graph (zero-length genome) skips exactly as before.
* `CalibrationResult._check_region_array` admits integer dtypes (inert here; W5a needs it). float32
  is still refused — the gate's real job is precision discipline, and integers are exact.
* `node_chain.build_node_chain`'s "rebuild the index" assert message now names the real condition:
  both offset arrays come from ONE payload, so a mismatch is not a stale index and rebuilding cannot
  fix it.
* ⭐ `substrate.PARTITION_MISMATCH_HINT` — both alignment guards now name the likely cause
  (`from_region_df(index.region_df)` instead of `from_index`). This is why the **~125 unlinted
  `scratchpad/` scripts were left unmigrated**: one message turns them from silent traps into
  self-diagnosing failures, which is worth more than 125 mechanical edits. All **69
  `scripts/debug/`** call sites WERE migrated.

### INDEXES — MIGRATED ADDITIVELY, NOT REBUILT (`scripts/debug/migrate_index_v8.py`)

Re-running `rigel index` re-derives transcripts, intervals, splice junctions and the splice-artifact
blacklist; any of those coming out differently would land inside this arm indistinguishably from the
partition change. The v8 graph is a pure function of `(transcripts, ref_lengths)` and **both are
already on disk**, so the four v7 indexes were migrated in place: build the graph, write two files,
bump the manifest, touch **nothing else**. ⭐ **Gated on P3 per index** — re-applying the merge must
reproduce that index's own `regions.feather` exactly — which simultaneously proves the reconstruction
from feathers is faithful and that the graph differs from v7 by the merge alone. Exact on all four.

The remaining four (`efflen_binding_sweep`, `gdna_benchmark_5mb`, `gdna_shortfl_5mb`,
`nascent_benchmark_5mb`) sat at **format 5** and were already unloadable before this arm; they were
fully rebuilt. Their node counts came out at exactly the pre-flight census (541 / 1,226 / 1,226 /
1,226) — an independent check on the census. **8/8 indexes at v8.**

### ⚠ TWO OPERATIONAL TRAPS THIS ARM CREATED, BOTH HIT IN THIS SESSION

1. ⭐ **The version gate is now a one-way door mid-session.** The `quick_3to1_5mb` v7 baseline ran
   refit=0 successfully, then refit=3 started a **fresh process**, picked up the bumped constant, and
   correctly refused the v7 index — leaving half a baseline. **Record every refit of a pre-change
   baseline BEFORE touching `INDEX_FORMAT_VERSION`.** Recovered without reverting the tree
   (`scratchpad/w1b_v7_from_cache.py`: score from the warm payload cache + `regions.feather`, which
   needs no index load); its refit=0 **control reproduced the recorded rows 16/16 exactly**, which is
   what makes its refit=3 rows a usable baseline.
2. **`--report` refuses `head`(v7) vs `w1b`(v8)** even though the object sets are provably identical
   here. That is W0.3's guard working as designed — it cannot know the partitions coincide. The
   comparison was done directly against the TSV. `pass0_oracle_bench` also gained a `--suite`
   argument and a **`suite` provenance column** (never waivable, not even with `--coarsen`: two
   suites are different genomes and no aggregation makes them comparable).

### TESTS

**1290 pass** (1281 + 7 partition-wiring + 2 net), ruff + `ruff format` clean, warning-clean.
Three tests were inverted or rescoped rather than deleted, each with its reason in place:

* `test_splice_graph.py::test_graph_is_optional_at_load` → **`..._is_REQUIRED_at_load`**. W1a's whole
  point was that the graph was optional; W1b's is that it is not.
* `test_region_index_alignment.py::test_neighbour_differs_per_ref` → scoped to
  **`..._holds_for_the_RETIRED_v7_partition`**, cross-referenced to the live statement
  (`test_adjacent_nodes_may_share_a_signature`). Neighbour-differs was v7's defining property and v8
  drops exactly it.
* `test_result_schema.py::test_rejects_non_float64_array` → split into accepts-integers /
  still-rejects-float32 / still-rejects-negative.
* `test_capture_eff_length.py`'s three `misaligned_index` tests **stay on the v7 frame deliberately**
  — v8 can no longer *produce* an interior exon edge from an index, so switching them to `from_index`
  would leave them passing while testing nothing. Added
  `test_incidence_is_correct_on_the_v8_partition_too` for the live path. ⚠ Both halves are on the
  W1b-clean delete/rewrite list below.

### ⛔ OWNER DIRECTIVE, 2026-07-29 — CONVERGE AND DELETE; NO LEGACY

> *"We don't need to maintain any backwards compatibility whatsoever… What I worry about is by trying
> to keep the old version intact for the sake of comparing, we run the risk of having a harder time
> cleaning it up and disambiguating the old legacy behavior from the new."*

Two decisions taken, and they amend the plan:

1. ⭐ **W7's partition deletion moves FORWARD to the next arm (W1b-clean).** After W1b the v7
   partition is already dead at runtime — nothing reads `regions.feather` or `boundaries.feather`
   during a scan, calibration or quant. Carrying it through W2–W6 buys only P3-against-a-shipped-
   artifact, and P3 survives as an in-test comparison (both partitions are pure functions of the same
   transcripts). Gated on **bit-identity**: nothing reads it, so nothing may move.
2. ⭐ **SHADOW MODE IS CANCELLED (plan W4).** Emitting both payloads from one scan — ~1.6 GB peak
   RSS, two payload classes, `scan_payload_v5.py` beside `scan_payload.py`, every consumer written
   twice — existed purely to give each W5 arm a legacy A/B. That is precisely the machinery the
   directive rules out. **W3 is unchanged** (Python reference + byte-for-byte tests, written first);
   W4 then replaces the C++ accumulator outright, gated on the reference and the oracle bench. The
   reference-oracle gate is the stronger of the two anyway — it checks against truth, not against the
   previous implementation.

### ⭐ W1b — THE PARTITION EFFECT, MEASURED (`quick_3to1_5mb`, 16 oracle conditions, +25.0 % nodes)

The bench cannot show it, so it was taken on the suite that splits most and has an oracle. v7
baseline recorded **before** the version bump; v8 after; identical scored mass (37,047,841) both
sides — the same fragments over more objects.

| | v7 (1,221 regions) | v8 (1,526 nodes) | Δ | better / worse / flat |
|---|---|---|---|---|
| **refit = 0 (pass-0)** | 0.0629 | **0.0698** | **+0.0069** (+11 %) | 3 / 8 / 5 |
| **refit = 3 (production)** | 0.0410 | **0.0395** | **−0.0014** (−3.4 %) | 5 / 3 / 8 |

⭐ **The Jensen confound is ~ZERO on this suite, so the sign is interpretable.** W0b measured a
0.739× gap on the human index — a coarse mwae reading 26 % better than the fine one from averaging
alone — which would otherwise make any v7-vs-v8 diff meaningless. Measured here by scoring **the same
v8 solve** at both granularities (region nodes, mass-weighted pooling — the only pooling the oracle
permits): **0.998× at refit=0, 0.999× at refit=3.** The gap is a property of sub-node heterogeneity,
and this suite's split nodes are nearly homogeneous. ⚠ Do not carry 0.998 to another suite.

⭐ **The effect is structural, not noise, and it splits on whether there IS gDNA:**

| stratum | refit=3 behaviour |
|---|---|
| `gdna300` (real gDNA present) | **v8 WINS**, 5/8 better — biggest at unstranded capture-off (**−0.0124**, **−0.0152**) |
| `gdna_none` (ZERO gDNA — truth is `f_g ≡ 0` everywhere) | **v8 LOSES**, worst at unstranded (**+0.0306**, **+0.0836**) |

**Reading, stated plainly.** A finer partition splits the same fragments over more, smaller objects,
so each carries less evidence. Where gDNA exists the extra resolution finds it and more than pays for
the thinner counts. Where there is *none*, every non-zero `f_g` is pure error and there is nothing to
win — the thin nodes just drift off zero. That is v5 §11.2's honest failure case (a short node is a
good `ρ_g` measurement and says nothing about `ρ_r`) appearing at the partition level, and §12.3's
"density cannot be resolved below the fragment length".

⚠ **Pass-0 pays and the refit recovers it** (+0.0069 → −0.0014): the imputation/hyperprior fills in
what the thinner nodes cannot measure. **The prior-free first pass is where the finer partition
costs**, which is a live input to the P1g/hyperprior work — not a reason to reverse the partition,
since the terminus visibility it buys is the entire point of v8.

---

## W1b-clean — ARM: DELETE THE v7 PARTITION (2026-07-29) — ⭐ BIT-IDENTICAL

**Owner directive** (above): converge, retire the old implementation, keep no legacy. After W1b the
v7 partition was dead at runtime — nothing read `regions.feather` or `boundaries.feather` during a
scan, calibration or quant. Carrying it to W7 bought only *P3-against-a-shipped-artifact*.

**⭐ GATE: bit-identical 32/32 at BOTH refits.** r0 `0.079005`, r3 `0.046675`, exact on all 7 columns
against the W1b arm. Nothing read it, so nothing could move — and nothing did. Suite **1280 pass**,
ruff + `ruff format` clean.

**Deleted:** `calibration/regions.py` **whole** (548 lines: the merge builder, the boundary builder,
both loaders, both validators) · `regions.feather` / `boundaries.feather` / their TSVs from the index
schema, the builder and the loader · `coarsen_nodes_to_regions` · `TranscriptIndex.region_df` /
`.boundary_df` · `scripts/debug/coarsen.py` · `scripts/debug/migrate_index_v8.py` (its migration is
done and its input no longer exists) · `tests/calibration/test_regions.py`,
`test_boundary_partition.py`. **src/ nets −434 lines.**

**Moved rather than duplicated:** `build_region_partition_arrays` → `splice_graph.build_node_partition_arrays`
(it is the graph→scanner adapter, and `regions.py` was its last tenant) · `load_ref_lengths` →
`index.py` (it reads an index feather) · `RegionArrays.from_region_df` → `.from_frame` (it takes a
partition frame; the v7 name was the last thing calling it a region table).

### ⭐ P3 DIED; I3b REPLACES IT, AND IS STRONGER

P3 — *merging v8's equal-signature neighbours reproduces `regions.feather`* — was the only
independent check on the v8 **signature** computation, and it cannot outlive its comparator.
Deleting it bare would have left `validate_graph`'s I3 checking nothing but the value RANGE.

So **I3b is now implemented**: `validate_graph` recomputes every node's signature from its
**midpoint**, by direct interval containment, and asserts equality with the stored value. Both halves
consume ONE definition of which intervals set which bit (`_signature_intervals`, shared with the
builder) and differ in the *evaluation* — cumulative difference arrays over the node index versus
midpoint containment — so the two can only agree by both being right. It needs no v7, it runs on
every graph including the human one, and it is checked at build **and** at load.

⚠ **Proven to FIRE** (`test_I3b_FIRES_on_a_corrupted_signature`): a hand-corrupted signature raises
`I3: … has signature 2, recomputed 0`. A validator that cannot fail is worth nothing.

### THE FOUR v7-PINNED TESTS, RESOLVED

* `test_capture_eff_length.py`'s three `misaligned_index` tests guarded a real bug — the incidence
  off-by-one on a partition whose regions contain interior exon edges — using the merged partition as
  the fixture. No index emits that geometry any more. **The coarse partition is now built by hand
  (`_coarsened`), in the test**, which is what it should always have been: the guarded property
  belongs to the FUNCTION, not to a shipped artifact. The fixture asserts it really does coarsen, so
  it cannot silently degenerate.
* `test_region_index_alignment.py::test_neighbour_differs_…` deleted — its subject is gone, and the
  live statement (`test_adjacent_nodes_may_share_a_signature`) already exists.

### ⚠ WHAT IS STILL OUTSTANDING

* **~125 unlinted `scratchpad/` scripts** call `index.region_df` and now break with an
  `AttributeError`. They are finished one-offs from closed investigations; deleting them is a
  separate decision for the owner (plan W7 already flags "enumerate its ~40 gate scripts by hand").
  `substrate.PARTITION_MISMATCH_HINT` covers the subtler failure mode for the ones worth keeping.
* `scripts/debug/**` (69 call sites) WAS migrated to `RegionArrays.from_index`.

### ⚠ ONE MEASURED BUDGET BREACH, REPORTED NOT BURIED

`validate_graph` with transcripts now costs **5.6 s** on the human annotation against graph doc §8's
**< 5 s** budget — I3b's midpoint recomputation added ~0.9 s. Two facts that scope it: it is a
**build-time** cost only (the load-time call passes no transcripts, so I3b/I4/I11 are skipped, and
load is unaffected), and it buys the only independent check on the signature that exists now that P3
is gone. Left as-is deliberately; the obvious recovery if it ever matters is to hoist the per-(ref ×
bit) `np.sort` out of the I3b loop. **Build total is 4.6 s + 5.6 s = 10.2 s, once per index.**

---

## W1c — ARM: the structural flags reach the solver (2026-07-29) — ⭐ BIT-IDENTICAL

**One thing varied:** `NodeStatics` gains the splice graph's 8 structural bits at every boundary
node. **Nothing reads them** — W2 will — so the gate is bit-identity, and it held: **32/32 exact at
both refits**, r0 `0.079005`, r3 `0.046675`. Suite **1289 pass**, ruff + format clean.

`splice_graph.build_boundary_flags_array(index) -> uint16[B]` → a `boundary_flags` keyword on
`calibrate()` → `NodeStatics.boundary_flags`, plus `is_terminus` / `is_splice_site`.

⚠ **The two index spaces are off by one per reference, in opposite directions**: `k` nodes give
`k − 1` contiguous edges but `k + 1` accumulator boundary slots. Slot 0 and slot `k` are reference
terminals — genuinely flag-less, not padding (a terminal is not a transition and never carried a
deposit either). A mistake here shifts every flag by one seam, is invisible in aggregate, and **would
pass this arm's own bit-identity gate**, so the tests check the mapping by **genomic coordinate**
against a really-scanned payload rather than by re-deriving the arithmetic.

⚠ **Deviation from the handoff's "new bool fields on `NodeStatics`": the RAW `uint16` is carried, not
pre-derived predicates.** Each W2 consumer wants a different combination, and plan F10 already
measured that P1G_SCOPE's specified predicate was nearly the **complement** of what it was meant to
replace. Carrying the bits makes a wrong predicate a one-line fix in its consumer instead of another
plumbing arm.

### ⛔⛔ A REAL BUG, FOUND BY THIS ARM'S OWN TEST — plan F1's second filter was WRONG

Writing the flags fixture surfaced that a real single-exon transcript's TSS and TES were **missing**.

Plan F1 specified **two** transcript filters: the event set `~is_synthetic`, and the flags/reaches
`~is_synthetic & ~is_nrna`, reasoning *"an nRNA span's ends are not real transcript termini"*. The
reasoning is right. **The predicate is not.** Measured on the human annotation:

| | |
|---|---|
| rows `is_nrna & ~is_synthetic` | **26,475** |
| ...of which `n_exons == 1` | **26,475 (all of them)** |
| ...of which are `RIGEL_NRNA_*` manufactured spans | **0** |
| rows `is_synthetic` | 203,052, **every one `is_nrna`** |

⭐ On a **non-synthetic** row, `is_nrna` does not mean "manufactured span" — it means **this real
transcript is single-exon, so its mature and nascent forms coincide**. Every manufactured span is
already `is_synthetic`, so `~is_synthetic` alone does exactly what F1's reasoning requires. The extra
clause deleted the termini of 26,475 real transcripts — **52,104 distinct terminus positions** — i.e.
precisely the visibility the entire v8 partition exists to buy.

**Fix: ONE filter, `~is_synthetic`.** F1's "TWO filters, two purposes" is struck. It is also a
simplification: `ev` and `ma` were two full `_Exons` builds of the same population.

**Effect on the artifact:** nodes UNCHANGED on all 8 indexes (so `partition_hash` is unchanged and
every cached payload stayed valid); flagged contiguous edges **992,068 → 1,043,595** at human scale,
**+51,527**. Bit-identical on the bench, because nothing reads flags yet.

### ⭐ AND IT REVEALED A MISSING INVARIANT — I13

1,043,595 is **every** contiguous edge. That is not a coincidence, it is forced: every cut is an exon
endpoint (I4), and every exon endpoint of a real transcript is either a terminus or a splice site.
The old filter left **95.1 %** — and nothing noticed, because the flags had no consumer, so no number
moved. **A whole class of defect was undetectable.**

**I13, now in `validate_graph`: each structural bit is set at EXACTLY the positions that generate
it — both directions.** No event without its flag, no flag without its event. It subsumes the design
doc's I10, needs no v7, and is **proven to fire both ways** (`I13: … is missing flag bit 0x0001` /
`… wrongly carries flag bit 0x0001`).

⚠ **One exemption, and it is real:** an interior interface may carry NO flag when adjacent exons of
one transcript are **bookended** (a zero-length intron) — an exon endpoint that is neither a terminus
nor a splice site. Measured **zero** times in GENCODE; constructed by G14. This is why I13 is stated
as *"the flags are the events"* and not as the tempting *"every interface carries a flag"*, which is
what I first wrote and which G14 immediately refuted.

---

## W1-REVIEW — ADVERSARIAL REVIEW OF W1b + W1b-clean + W1c (2026-07-29) — ⭐ BIT-IDENTICAL

Five independent review lenses (index-space · wiring · elegance · tests · doc-accuracy) over the whole
change set, each finding refuted-by-default by a separate verifier that had to reproduce the defect
before it counted. **64 agents, 59 findings raised, 42 survived.** Fixes below; the bench is
**32/32 exact at both refits** after all of them, r0 `0.079005` / r3 `0.046675`. Suite **1280 pass**.

### ⛔ THE ONE THAT MATTERED — I13 WAS SUPPLYING FALSE ASSURANCE

**`_terminus_and_splice_events` holds the only strand-dependent line in the flag path** (a NEG
transcript's 5′ end is its genomically HIGH edge). Nothing pinned it: every flag fixture was `+`-only,
and `is_terminus` ORs TSS|TES so it cannot expose a swap. Worse — **I13 called the builder's own event
emitter**, so it validated the searchsorted-vs-isin *placement* and never the event *definition*.

⭐ **Verified by the reviewer, not argued:** delete the minus-strand swap and the full **1,289-test
suite still passes and `validate_graph` still accepts the graph.** A control mutation
(donor/acceptor swap) fails 2 tests, so the suite discriminates — it was specifically blind to this.

That is exactly the failure mode I13 was added to end, on a different axis, and it is worse than no
check because it reads as one. **Two fixes, both landed:**

* ⭐ **`_events_independently`** — I13 now re-derives all eight event sets by a **different
  algorithm**: per-transcript `min(start)`/`max(end)` via `np.minimum.at`, against the builder's
  cumulative-exonic-offset arithmetic (`before == 0`). The I3b discipline, now applied to the flags.
  Confirmed: the swap mutation is **REJECTED** — `I13: position 500 wrongly carries flag bit 0x0002`.
* Three tests: `test_NEG_strand_termini_are_the_other_way_round` (with its POS mirror, so it cannot
  pass by being wrong both ways), `test_I13_FIRES_on_a_SWAPPED_neg_strand_terminus`, and —

**⚠ A SECOND HOLE in I13, also found:** the builder's emitter appends an entry only for a class that
*has* events, so **a bit whose class is empty on a reference was never compared**. A spurious
`TSS_NEG` on a POS-only reference passed silently. `_events_independently` returns all eight, empty
position set included. Pinned by `test_I13_FIRES_on_a_bit_whose_event_class_is_EMPTY`.

### ⭐ PERFORMANCE — THREE MEASURED WINS, ALL FROM THE ELEGANCE LENS

| | before | after |
|---|---|---|
| `validate_graph` at **load** (human) | **1.91 s** | **0.23 s** (8.3×) |
| `build_splice_graph` (human) | 4.6 s | **3.1 s** |
| `validate_graph` **with transcripts** | 5.6 s | **2.8 s** — ⭐ back inside §8's 5 s budget, *while now also running I13* |

* **I7 was doing the exact O(n_refs × n_nodes) string scan the I1 comment 45 lines above says it
  removed** — 286 full-array object-dtype comparisons over 1.04 M rows, **1.68 s of a 1.91 s** load
  validation, computing a number `slices` already held. The comment warning about the trap did not
  stop the trap being re-laid three blocks later.
* **`_group` and `_ref_slices` were two implementations of "reference name → rows"**, the first a
  Python loop over every exon row (1.0 s of the build). One name, one implementation.

### SIMPLIFICATIONS LANDED

* ⭐ **`_is_real` — ONE transcript filter, module-level.** It was duplicated verbatim in the builder
  and the validator, *in the module whose headline defect was a filter divergence*: I13 was
  validating the flags against its own private copy of the predicate, so the two could drift and
  still agree with each other.
* ⭐ **Ten names for four objects deleted.** `ev = ma = …`, `ev_i = ma_i = mi = …` were aliases left
  from the retired two-filter design; ruff's F841 flagged them the moment the first alias was
  touched. Now `ex` / `intr` / `by_ref` / `intr_by_ref`, one name each, threaded through every helper.
* **`pipeline._check_region_payload_alignment` deleted** — a verbatim duplicate of
  `CalibrationSubstrate._check_alignment`, called microseconds before it. The `payload is None`
  guard moved into the survivor. Two copies of an invariant is one too many.
* **`RegionArrays.order` deleted** — a stored field with zero readers anywhere in the tree.
* **Three "index predates version 8" guards deleted** — the loader cannot produce that state.

### CORRECTED MEASURED CLAIMS (the doc-accuracy lens, all re-measured here)

| claim | was | ⭐ measured |
|---|---|---|
| human termini invisible to a merged partition | 59.5 % | **53.4 %** (232,451 / 435,291) |
| positions that are BOTH terminus and splice site | "the majority" | **0.99 %** (10,337 of 1,043,595) |

⚠ **59.5 % was itself computed under the buggy `~is_synthetic & ~is_nrna` filter** — the headline
justification for the whole partition was still reading the annotation the way the bug did. Under the
shipped filter it is 53.4 %. Also measured, and it validates I13's exemption: **zero** contiguous
edges carry neither bit, i.e. GENCODE contains no bookended exons.

### OTHER CONFIRMED DEFECTS FIXED

* **`INDEX_FORMAT_VERSION` 8's history entry was false in four places** — it named a deleted function
  (`build_region_partition_arrays`), and claimed `regions.feather`/`boundaries.feather` and the P3
  gate were still live. It is the block a loader author reads.
* **The load-time comment claimed I3b runs at load. It does not** — I3b/I4/I11/I13 are all gated on
  `transcripts`, which `load()` never passes. Both the comment and `validate_graph`'s docstring now
  say exactly which invariants run where, and that a drifted `signature`/`flags` column loads clean.
* **`test_a_pre_graph_index_is_refused_at_load` never reached the guard it claimed to test** — it also
  rewrote `format_version` to 7, and the version gate fires ~80 lines earlier. Deleted;
  `test_graph_is_REQUIRED_at_load` does it properly.
* ⭐ **`calib_cache` no longer caches `boundary_flags`.** It did, with a comment asserting they
  "cannot drift from the payload independently of `partition_hash`" — **false**: that hash covers
  `nodes.feather` only, and today's flag fix rewrote every `edges.feather` while leaving every
  `nodes.feather` byte-identical, so a stale-flags cache would have verified CLEAN and fed every W2
  A/B the pre-fix flags. They are derived from the index now. `partition_hash`'s docstring says so.
* **The D6 fix left "the AVERAGE" prose at five sites** in `priors.py` / `capture_eff_length.py` —
  the exact prose D6's own analysis blames for the bug. `gbl = E[min]/2`, so `gbl_r + gbl_{r+1}` IS
  `½·(E[min_r]+E[min_{r+1}])`: one quantity, two forms, and naming it after the wrong one is how the
  ½ got applied twice. All now state the SUM form first.
* `signature.py` — the lowest layer everyone imports — still defined a region as "a **maximal**
  interval over which this signature is constant", the exact property v8 abolished, and cross-
  referenced the deleted module.
* Dead no-op in a new test; a duplicate test; stale `P2′`/`P3` section headers.

### ⚠ NOT ACTED ON, DELIBERATELY

* **`CalibrationResult`'s widened dtype gate guards a producer that does not exist yet** (confirmed,
  low). Kept: the owner named it explicitly as a W1b item, and W5a is two arms away. Recorded so it is
  a decision rather than an oversight.
* **17 findings were REFUTED** by their verifiers, including a claimed `grp.index` hazard in
  `build_boundary_flags_array` and a claimed `partition_hash` collision — the refutations are in the
  workflow journal.
