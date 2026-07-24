# `_scan` / message-system refactor — survey + meticulous plan (Phase A structural)

**Status (2026-07-21):** the survey-first plan for restructuring `bp_solver.node_sweep._scan` (the per-edge
message construction) into a clean, intuitive, bug-resistant form — **behavior-preserving** (golden byte-identical)
— after the τ-gag fix landed (`243bd5ef`). Follows [`message_system_implementation_plan.md`](message_system_implementation_plan.md)
Phase A. The behavioral derivation is [`message_system_derivation.md`](message_system_derivation.md).

---

## 1. Survey — what `_scan` is today, and why it's hard to read

`node_sweep` (~600 lines) contains, in one function: the geometry/statics unpacking, the τ-seed computation
(I_strand deadband + struct_lock), the **`_scan` inner closure** (~200 lines, the per-edge loop, run twice —
forward + backward), `_comb`, and the final solve + the `_capture` diagnostic dump.

**`_scan` interleaves EIGHT distinct concerns in one loop body** (lines ~502-697), per edge:

| # | concern | current form |
|---|---|---|
| 1 | edge geometry | `md, egd, erd, sm` from `MSd/EGd/ERd/MSs` face arrays |
| 2 | transfer variance `s2t` | `var_proj[i] + (mu_proj[i]−mu_proj[lsrc])²` |
| 3 | source composition | `fg_s, fr_s, fp_s, fn_s` from the running `(μ_λ, μ_θ)` belief |
| 4 | **evidence / identifiability** | `lock_s, lam_ev, ev_lam, th_ev, ev_th` → `v_logfg/fr/fp/fn` |
| 5 | **emission gates** | `emit_g/emit_p/emit_n`, `has_comp_ev` (the τ-gag fix lives here) |
| 6 | source densities | `rho_g, rho_pos, rho_neg` (with mature ±absorb) |
| 7 | **message mode** | `use_shift` → the cliff shift (`Mg/Mp/Mn`) OR the density mode |
| 8 | **message precision** | prediction (τ-gated) + spliced measurement, per component |
| — | fold + τ cavity | `_fold_lambda` + the Jacobian-converted `tau_lam[i] +=` |

**Why it's bug-prone (the τ-gag was exactly this):** concerns 4/5/8 are the identifiability logic, but they're
scattered across ~40 interleaved lines with the geometry and mode math. The spliced-measurement-vs-gate bug hid
precisely because "when does a message emit" (5) and "how confident is it" (8) were not one readable unit. There
are **six parallel accumulator arrays** (`amg/apg/amp/app/amn/apn`) hand-indexed per component. The
**prediction-vs-measurement dichotomy** (the heart of the fix) is implicit in a `has_comp_ev` ternary.

**Two hard constraints the refactor must respect:**
1. **Performance.** `_scan` is a **sequential** loop (each fold updates the running belief the next edge reads),
   so it *cannot* be vectorized, and it deliberately uses Python-scalar locals (the `order_list` int-list, inlined
   arithmetic) because per-edge numpy indexing / function-call overhead is the measured cost. A naive
   "extract every concern into a function called per edge" would regress the sweep's runtime.
2. **Golden byte-identity.** The arithmetic must be *bit-for-bit* preserved through the refactor. Any reordering
   of float ops that changes rounding fails the golden. So extraction must move whole expressions unchanged.

**Also in scope (cleanup):** the intermingled **default-off δ-pin experimental gates** committed pending removal
— `config.gdna_prior_enrichment_condition`, `simplex_logodds._STAGE0_FLOOR` (`RIGEL_STAGE0_FLOOR`),
`calibrate` `RIGEL_ADD_GONLY`, and the `enrichment_condition` plumbing threaded through `node_sweep`/`calibrate`.
These are dead weight that obscure the real logic; removing them is byte-identical (all default-off).

---

## 2. Target design — three readable units, one hot loop

The refactor's spine is to make the **identifiability model explicit** and to separate the **prediction** and
**measurement** channels structurally, *without* turning the hot loop into a function-call storm.

### 2.1 The `Evidence` unit (extract — runs ONCE per sweep, not per edge)
The τ-seed block (I_strand deadband + struct_lock, `node_sweep` lines ~414-452) becomes a small pure helper
returning a named structure:

```python
@dataclass(frozen=True)
class StrandEvidence:          # computed once, before the forward/backward scans
    tau0_lam: np.ndarray       # I_strand (differential-κ deadband) — the composition Fisher seed
    tau0_th:  np.ndarray       # the tilt Fisher seed
    struct_lock: np.ndarray    # I_struct: signature composition-certainty (intergenic REGIONS today)

def _compile_strand_evidence(statics, fg_loc, kappa, od_g, od_r, n_gdna_obs, n_rna_obs, chain) -> StrandEvidence
```

This is a *pure, unit-testable* function (no loop, no perf concern) — and the natural home for the later
`I_spliced` seed (Phase B remainder) and the `I_struct`-at-TSS/TES extension (Phase D). **It gives the
identifiability model one name and one place**, which is the single biggest bug-resistance win (the τ-gag was a
hole in exactly this).

### 2.2 The message-mode functions (extract — pure, called per edge but cheap)
The two mode formulas become named module-level pure functions (they already are self-contained expressions):

```python
def _shift_mode(rho_c, E_c_dst, denom, comp_floor) -> float      # cliff-invariant composition shift
def _density_mode(rho_c, E_c_dst, md) -> float                    # observed-total density anchor
```

Named, docstring'd, individually testable against `cliff_message_derivation.md` — and the call sites read
`mo = _shift_mode(...) if use_shift else _density_mode(...)` instead of two inline `math.log(max(...))` blobs.
(Micro-cost: a function call per component per edge; measure — if it regresses, keep them as module functions but
inline the two callers via a local alias. Expected negligible vs `_fold_lambda`'s grid work.)

### 2.3 The precision block — make prediction vs measurement STRUCTURAL (the bug-resistance core)
The heart of the fix, currently a `has_comp_ev` ternary, becomes an explicit, commented two-channel computation
kept **inline** (hot path) but sectioned and named so the dichotomy is unmissable:

```python
# CHANNEL 1 — the deconvolution PREDICTION: trusted only with composition evidence (τ>0 or a lock).
pred_pr_p = n_eff / (n_eff * (v_logfp + s2t) + 1.0) if has_comp_ev else 0.0
# CHANNEL 2 — the spliced MEASUREMENT: independent evidence, always credited (a splice read IS pure RNA).
meas_pr_p = SPs[lsrc] / (1.0 + SPs[lsrc] * s2t) if SPs[lsrc] > _EPS else 0.0
pr_p = pred_pr_p + meas_pr_p
```

Same arithmetic as today (byte-identical) — but the two channels are *named and adjacent*, so "a measurement must
never be gated by the prediction's evidence" is self-evident in the code, not a comment. A **regression test**
already pins it (`test_tau_gag_fix_*`); this makes the *code* say what the test checks.

### 2.4 Keep — the loop, `_comb`, the fold, the f32 AMBIG cube, the `_capture` dump
The sequential loop stays a loop (perf). `_fold_lambda`, `_comb`, the final `_solve_nodes_logodds_all`, and the
diagnostic `_capture` are already reasonable units — untouched except for consuming the extracted helpers.

---

## 3. Incremental, test-gated steps (each byte-identical unless noted)

Every step gates on: `pytest tests/calibration` (249) + the **golden** suite (byte-identical) + lint. Behavioral
steps additionally gate on the `gdna_none` FP guard + the `pass0_bench.py` mwae table.

- **A0 — remove the δ-pin gates (byte-identical, high-value cleanup). ✅ DONE (`7a51351c`).** Deleted
  `gdna_prior_enrichment_condition` (config), `_STAGE0_FLOOR` (simplex_logodds), `RIGEL_ADD_GONLY` (calibrate),
  and the `enrichment_condition` plumbing (node_sweep/calibrate) + the now-unused `import os`. 0 golden files
  changed; 1132 tests; lint clean.
- **A1 — extract `_compile_strand_evidence` → `StrandEvidence`. ✅ DONE (`7db05195`).** The τ-seed block moved
  verbatim into a frozen dataclass + pure function; `node_sweep` calls it. Unit tests: the deadband kills
  I_strand unstranded + gates a gDNA-free library; `struct_lock` = locked REGIONS only. 0 golden files changed;
  1132→1134 tests. This is the home for the §6B solvability model (Phase C) + I_spliced (B) + I_struct@TSS/TES (D).
- **A2 — extract the mode functions `_shift_mode` / `_density_mode`.** Replace the two inline blocks; unit-test
  each against a hand-worked cliff example. *Golden byte-identical.* Measure sweep runtime (expect flat).
- **A3 — section + name the precision block (prediction vs measurement).** Introduce `pred_pr_*` / `meas_pr_*`
  locals (2.3) for all three RNA channels + the RNA-total factor. *Golden byte-identical* (pure renaming of the
  existing sum). The `test_tau_gag_fix_*` tests already guard it.
- **A4 — a short module docstring / section banner** mapping the eight concerns (the §1 table) to their code
  regions, so the next reader has the map this survey gave us.

**Explicitly NOT in Phase A** (they are behavioral — separate phases, each measured): `I_spliced` as a τ seed
(rest of B), the `ev_lam=∞` quirk fix (B), skip-unidentifiable-solves (C), `I_struct` at TSS/TES + the
unequal-gate lower bound (D). Phase A only makes the code *ready* for them by giving the identifiability model a
home.

---

## 4. Bug-resistance measures (the "protects itself" goal)

- **One name for identifiability** (`StrandEvidence` + `_compile_strand_evidence`) — the τ-gag was a hole here;
  a named unit with its own tests makes a future hole visible.
- **Prediction vs measurement is structural** (§2.3) — the exact dichotomy the bug violated is now spelled in
  the code, not a comment, and pinned by falsifier tests.
- **Mode functions are pure + tested** — the cliff arithmetic can't silently drift.
- **The §1 concern-map banner** — the loop body's eight jobs are documented in-place.
- **No new magic numbers; no env flags** (A0 removes the last ones). Every constant stays canonical/derived.
- **Golden byte-identity is the gate on every step** — the refactor cannot change behavior by construction.

---

## 5. Risks & mitigations

| risk | mitigation |
|---|---|
| golden drift from float-op reordering | move whole expressions verbatim; gate every step on the golden |
| perf regression from per-edge calls | only §2.2 adds per-edge calls (cheap, pure); measure sweep time at A2; fall back to inlined aliases if needed |
| δ-pin removal breaks a caller | A0 is default-off everywhere; run the full suite + a `calibrate` smoke on one scenario |
| scope creep into behavioral change | Phase A is byte-identical ONLY; behavioral work is B/C/D, explicitly fenced |

---

## 6. Recommendation

Do **A0 → A1 → A2 → A3 → A4** in order, each its own byte-identical commit gated on golden + tests. A0 (δ-pin
removal) delivers immediate cleanliness; A1 (the `Evidence` unit) delivers the structural bug-resistance the τ-gag
motivated. Then resume the behavioral phases (B remainder → C → D) on a codebase that now has a home for the
identifiability model.
