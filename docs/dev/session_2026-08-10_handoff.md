# Session handoff — 2026-08-10

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When an item settles, MOVE it to its permanent
    home and delete it here in the same edit.

    Everything this session did, in the order it happened, with the state a successor needs. Commits
    `f470a570 … 69a85be2` plus the segfault fix. ⭐ **Read §0 and §6 first.**

---

## 0. ⭐⭐⭐ WHERE THINGS STAND

**Suite: 0 failed / 3,224 passed / 2 skipped / 9 xfailed.** `ruff check src/ tests/ scripts/` clean.
`src/rigel/` builds; the C++ was rebuilt after every native edit.

⛔⛔ **EVERY SCAN CACHE IS INVALID AND NOTHING HAS BEEN RE-SCANNED.** `payload_schema_digest` moved from
`19ee4ba867ff0441` to `cc12ee0a113d8c76` when `node_contained_inv_length_sum` was renamed. **No design
instrument under `scripts/design/` will run until the panels are rebuilt** — that is phase 5 of the plan
in `counts_densities_paradigm.md` §9 and it is the first irreversible step. Everything before it is a
code revert.

| flag | state |
|---|---|
| `CalibrationConfig.message_propagation` | `False` (owner, unchanged) |
| `CalibrationConfig.length_likelihood` | ⛔ **GONE** — the channel was deleted |

---

## 1. THE LENGTH CHANNEL — MEASURED TO A VERDICT AND DELETED

`length_likelihood` is gone (`b7ed7a0b`). The verdict, on the DRAINED arm: **its answer is not a function
of the pmf gap at all.** Interpolating `gdna_pmf` toward `rna_pmf` over twelve orders of magnitude, at a
gap of ~1e-9 bp it reported **0.72 / 0.59 / 0.72** on libraries whose truths were **0.00 / 0.00 / 0.57**,
and *growing* the gap made every one better. So the fragment-length MODELS were exonerated and every
"fix the pmf" repair is refused.

**Mechanism**: a Gaussian log-likelihood is asymptotically LINEAR in the composition, so its argmax is a
SIGN saturated at a grid endpoint, while participation (an exact `!=` guard) and precision (`1/Var` of a
near-uniform posterior) are THRESHOLDS that fire the instant the gap is nonzero. Amplitude fades with the
gap; influence does not. `TRAPS: a-linear-likelihood-emits-a-sign`,
`TRAPS: amplitude-fades-influence-does-not`.

⛔ **Do not rebuild it.** The instruments that priced it were deleted with it; they are in `f470a570`.

---

## 2. THE PRIOR AUDIT — the density→count reconciliation is CORRECT

`f_g` is a **count share** everywhere, never a density share. `node_gdna_geometry` returns `eff_gdna`, so
`rho_c = f_c·M/E_c` and `Σ_c rho_c·E_c = M` — which IS the statement that the fractions apportion the
COUNT. The one density→fraction conversion (`density_model.py`) multiplies by `eff_gdna`, the gDNA's own
effective length. The fragment-length difference is carried by `E_g` vs `E_r` and never by the fraction.

**Three currencies, not two** (`counts_densities_paradigm.md`): INCIDENCE (`count`, not conserved),
FRAGMENT (`mass`, conserved), DENSITY. ⭐ The symmetry is complementary — at a NODE incidence = fragment
is free and density needs a model; at an EDGE the reverse.

### 2.1 ⛔ TWO DEFECTS FOUND, NEITHER FIXED

1. **`pipeline.py:993-998` sums INCIDENCES and calls them fragments.** Applying `q = mass/count` to the
   edge terms moves the reported library `f_gdna` **0.3781 → 0.4214** against a truth of 0.5085 — **+4.3 pp
   toward truth**, on all four conditions tested. ⭐ This is the highest-value cheap fix outstanding.
2. **`mass_gdna_edge` / `mass_rna_edge` are incidences under a `mass_` name** while `mass_gdna_node` is a
   real fragment count. Defect 1 is the direct consequence of defect 2.

⚠ And **there is no `sj_mass`**, so the junction axis cannot be converted to fragments at all — a fully
conserved library fragment count is not computable today.

### 2.2 MEASURED AND CLOSED

`q` is measured pooled and applied per component. Real, but **placement-driven, not length-driven** (the
equal-length null carries the LARGEST error), and bounded at `Δphi` ≤ **0.6 pp** on the total prior —
below calibration's own noise floor. **Recorded, not fixed.** Instrument: `edge_q_population.py`.
`TRAPS: a-pooled-conversion-applied-per-component`.

---

## 3. THE NODE DEPOSIT IS NOW `1/OPPORTUNITY` (`69a85be2`)

`node_contained_inv_opportunity_sum` deposits `1/(ell − w + 1)₊`, so the node channel is a **density**:
`E[Σ] = rho` for any length distribution. Previously `1/L`, which read **25.67** density units for short
fragments and **1.60** for long ones at the same true density — a 16× swing driven by the length set alone.

* Falsification written first and **verified failing**; then the fixed code was broken three ways
  (off-by-one, revert to `1/L`, use `1/ell`) and **3/3 gates fired**.
* C++ byte-identical to `_accumulator_reference.py` across every named case, every arbitration outcome
  and 10,000 random fragments.
* ⭐ **The rename IS the cache-key fix.** `payload_schema_digest` hashes field NAMES, so a deposit-RULE
  change moves it not at all — stale caches would have been silently accepted
  (`TRAPS: a-hash-that-misses-its-artifact`, third occurrence). The bank had to be renamed anyway, and
  the rename invalidates every cache. Verified: a real ladder cache is now REJECTED.

⛔ **Nothing consumes the new channel yet.** Wiring `density_model` to it is phase 6 and is subject to
`panel-before-src`; it is worth ~1.1 % on the density surface (§8.2 of the paradigm doc).

---

## 4. ⭐⭐⭐ THE SEGFAULT — DIAGNOSED AND FIXED

**Symptom.** `rigel._bam_impl.Accumulator(...).deposit(start=…, end=…, hypotheses=())` →
`EXC_BAD_ACCESS (code=1, address=0x0)`. Deterministic.

**Root cause — `accumulator.cpp`, the arbitration block.** With `n_hypotheses == 0`:

1. the scoring loop never runs, so `survivors` is empty;
2. the `remove_if` leaves it empty; the "would empty the set" re-fill iterates `n_offered = 0`;
3. `survivors.size() > 1` is false, so the deferral is skipped;
4. **`const GapHypothesis& chosen = fragment.hypotheses[survivors.front().index];`** —
   `survivors.front()` on an **empty vector**, indexing a **null** `hypotheses` (the binding sets
   `offered.hypotheses = spans.data()`, and `.data()` on an empty vector is `nullptr`). 💥

**Why it survived.** The executable specification defaults `hypotheses` to `UNSPLICED_ONLY` — *"the
degenerate case is the general case, not a branch"* — and the scanner always offers the genomic path, so
**nothing in production or in the suite could reach it.** The native binding, by contrast, has **no
default**, so a caller writing the natural `hypotheses=()` hits it immediately. It was found while trying
to build a deposit-behaviour digest.

**The fix.** The accumulator now treats an empty set as the unspliced-only set, which is not a fallback —
it is what the specification says the default *is*:

```cpp
static const GapHypothesis kUnsplicedOnly{nullptr, 0, STRAND_NONE, nullptr, 0};
OfferedFragment offered_or_default = fragment;
if (fragment.n_hypotheses == 0) { … n_hypotheses = 1; }
const OfferedFragment& arbitrated = offered_or_default;
```

**The gate.** `test_an_EMPTY_hypothesis_set_is_the_UNSPLICED_ONLY_set_and_does_not_CRASH` asserts not
"does not crash" but that an empty set **deposits identically to `UNSPLICED_ONLY`**, with full
reference/native parity. ⛔ Perturbed: disabling the guard **segfaults the test process (exit 139)**;
restoring gives 10 passed.

### 4.1 ⚠ TWO METHOD ERRORS I MADE FINDING IT — worth not repeating

1. ⛔⛔ **`python … | tail` reports TAIL's exit code, not python's.** Three "successful" isolation runs
   were segfaults whose `OK` line simply never printed. **Always capture the exit code without a pipe**
   (`cmd > file 2>&1; echo $?`) when probing for crashes.
2. ⚠ **Unflushed stdout hides where a crash happened.** The original probe printed nothing at all, which
   read as "died at construction"; it had actually died at the end, with the loop output still buffered.
   Use `flush=True`.

### 4.2 ⛔ STILL OWED — the deposit-behaviour digest

The rename closed *this* cache hole, but **a future deposit-rule change that keeps every field name would
slip through again**. The intended repair is a behavioural digest (run the accumulator over a fixed tiny
synthetic case at key time and hash the banks). ⭐ **It is now possible** — the segfault that blocked it is
fixed and `Accumulator` is directly constructible:

```python
a = Accumulator(cuts=…, node_types=…, max_length=…, ref=0)
a.set_junctions(np.zeros(n_cuts+1, np.int32), np.zeros(0, np.int32), np.zeros(0, np.int8))
a.deposit(start=…, end=…, hypotheses=())      # now safe
```

⚠ It needs no version number and is deterministic *because* the channels are integers.

---

## 5. THE FIXED-POINT QUESTION — SETTLED, INTEGERS STAY

Worst case (a fragment exactly filling its node, `A = 1`, quantum `2³²`, **every** one of 1e8 fragments on
that one object): `4.295e17` against a `uint64` ceiling of `1.845e19` — **42.9× headroom**, on a scenario
already physically impossible. ⛔ And the integers are not about precision: **integer addition is
associative**, which is what makes the tally byte-identical across worker counts on a multi-threaded
scanner. Going float would end the C++↔reference byte-identity gate. The owner raised switching to float;
the arithmetic means it is not needed.

---

## 6. ⭐⭐ WHAT TO DO NEXT, RANKED

1. ⭐⭐⭐ **`pipeline.py`'s `q` conversion** (§2.1 defect 1) — arithmetic, +4.3 pp of the deliverable, and
   the spliced half needs `q_spliced = edge_spliced_mass / edge_spliced_count`. ⛔ Needs a re-scanned
   panel to verify, so it follows phase 5.
2. **Phase 5 — re-scan the panels.** Nothing measurable happens until this is done. ⚠ Owner has not yet
   chosen full (36 ladder + 8 flgap) vs the working subset (g00/g50/g98 × capture, plus flgap).
3. **The deposit-behaviour digest** (§4.2) — now unblocked.
4. **`mass_*_edge` → `count_*_edge`** rename (§2.1 defect 2).
5. **Re-price message propagation** — the only remaining lever for unstranded × capture-ON, one config
   flag, already measured to help exactly that stratum. `ROADMAP.md` §0 says re-price, never inherit.
6. **Moment tests** — `contained_moments`/`crossing_moments`/`build_slot_moments` moved to
   `effective_length.py` **without their tests**, which went with the deleted channel. Untested geometry
   in a live layer-2 module.
7. `sj_mass` decision; the two dead banks (`node_contained_length_sum`, `edge_unspliced_length_sum`);
   `second_pass.py:443-445`'s comment; `module_census.py` re-run.

⛔ **`edge_spliced_mass` is NOT dead — keep it.** It is the spliced half of the fragment conversion and
the reason it looks unused is that the conversion is not being done (item 1).

---

## 7. THE DOCS THIS SESSION LEFT

| doc | what is in it |
|---|---|
| `counts_densities_paradigm.md` | ⭐ the three currencies, the conversion matrix, the site audit, the ruling, the oracle-free seam gate (§7), the model-free node density (§8), the implementation plan (§9) |
| `cleanup_and_prior_audit.md` | the density→count audit, the bank ledger, the remaining-work table |
| this file | the session record |

⭐ Permanent homes already updated: `TRAPS.md` (three new named rules), `ROADMAP.md` §0/§1.4,
`DESIGN.md` §3.1c, `CLAUDE.md`'s instrument table.
