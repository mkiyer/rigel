# Graft precision — implementation plan (item 1)

**Status:** PLAN, for review before any code change. **Date:** 2026-07-24.
**Derivation:** `graft_message_precision_design.md` (T1/T2/T5, MC-validated). **Simplification (owner):**
Poisson counts, `ω = 0` — precision knowingly optimistic at high counts on real data.
**Primary goal (owner):** land the fix *and* leave the codebase **simpler, more intuitive, more
maintainable** than before. This plan is written to reduce net complexity, not add a patch.

---

## 0. The health problem this fix must not worsen — and should help

`bp_solver.py` is **1870 lines** with **two** solve paths (`_scan` legacy ~470 lines + `_unified_solve`
~300-line nested closure) and a thicket of inline lambdas (`_rho_faces`, `_pin`, `_pin_v`, `_relay`,
`_transport`, `_fuse_v`, `_damp`). The message *math* is smeared through that closure. Meanwhile the codebase
already has the RIGHT pattern next door: `enrichment_frame.py` and `gdna_intron_factory.py` are **pure,
frame-agnostic arithmetic modules** (arrays in, arrays out, no chain/substrate/solver state) pinned by
closed-form unit tests that "cannot drift with solver behaviour."

**T1/T2/T5 are pure arithmetic.** They belong in that layer, not inlined. Putting them there makes the
closure SMALLER and the math FINDABLE and TESTABLE — the opposite of a patch. That is the design principle of
this task: **every line I add to a pure module lets me delete more than one line of inline closure.**

## 1. What changes, in one sentence

Move the message-precision math into the pure arithmetic layer as **named, unit-tested functions**, and make
`_transport` *call* them — replacing the inline `_dv` (σ²_transfer damping) + `pr += S` (measurement
addition) block. The **mode path is untouched** (the density relay + reframe already give a correct mode:
measured exon msg f_g 0.682 vs truth 0.677); only how confident each factor is changes.

## 2. Why this is the right first fix (the mechanism, measured)

Undecided boundaries (34 % of boundary mass, local f_g ≈ 0.5) currently emit RNA at **2.6× the precision of
decided ones** and poison the exons — most of the unified solver's 0.1280-vs-0.0949 gap. Under T2 an undecided
boundary has a wide `Var(λ)`, so `Var(log k)` is large, so **its message correctly goes quiet.** Mechanistically
targeted at the measured defect, not a general tune.

---

## 3. The refactor — where the math lands

### 3.1 Home: extend `enrichment_frame.py` (do NOT create a new file)

`enrichment_frame.py` already owns the message arithmetic: `reframe_density`, `density_mode_logfrac`,
`k_from_belief`, `f_g_from_k`, `composition_logvar`. The precision laws are the same family. Adding them there
keeps **all** message math in ONE pure module a reader can hold in their head — more intuitive than scattering
a second file. Two new pure functions, each mapping 1:1 to a validated theorem:

```python
def graft_var_log_k(w_mu, f_g_src, var_lam_src, n_b, n_s):
    """T2 — Var(log k) for a boundary→exon message, in the source's logit(f_g) frame.
    Share-weighted (w_mu) so a measurement's confidence does not inflate the prediction (§3.1).
    Poisson counts (ω=0, owner). Degenerate: w_mu=0 ⇒ returns Var(λ) exactly."""

def message_factor_precision(var_log_k, f_g_dst):
    """T5 — the destination Jacobian: (prec on log f_g, prec on log f_R) at the DESTINATION's fraction.
    Returns precisions (1/Var), 0 where Var is non-finite."""
```

* **1:1 with the MC harness.** `message_precision_mc.py`'s T2/T5 blocks become the docstring examples and the
  unit tests — the functions are *defined by* the validated algebra.
* **`composition_logvar` (Var log ρ_tot) vs `graft_var_log_k` (Var log k) are DIFFERENT quantities** (a sum of
  densities vs a ratio) — co-located but clearly named and documented so no one conflates them.
* **Naming debt noted:** once the peel (item 2) + transfer (item 3) precision also live here, the file is the
  message-arithmetic module, not just "enrichment." Rename `enrichment_frame.py → message_arithmetic.py` in
  Phase 2 (a pure rename, deferred so this task stays reviewable). Flagged, not done now.

### 3.2 Call-site: `_transport` gets SMALLER

Replace the precision block (currently `_dv(pg/pp/pn)` + the two `pr += S` graft lines, ~lines 1627–1634)
with: gather the six inputs (§3.4) → `graft_var_log_k(...)` → `message_factor_precision(...)` → assign
`tpg/tpp/tpn`. The `_damp`/`s2t`/`_dv` machinery on the precision path is **deleted** from `_transport` (T1:
the reframe is enrichment-exact, so σ²_transfer here is a pure double-count — design §4). Net: fewer inline
lines, and the arithmetic is behind a name.

### 3.3 No new flag

The whole unified path is already gated (`RIGEL_UNIFIED`, default off) and unshipped, so there is nothing to
protect. **A/B via git**, not a runtime toggle — adding `_UNIFIED_GRAFT_PREC` would be exactly the flag-cruft
we are trying to shed. The old precision is simply replaced; if it regresses, revert. (This also means one
fewer flag than my first draft proposed.)

### 3.4 Inputs and where each comes from (all already in scope — no new plumbing)

| quantity | source | notes |
|---|---|---|
| `Var(λ^s)` | `lvar_loc[s]` (local solve, line 686) | wide for undecided → quiet msg. **Guard inf/nan.** |
| `f_g^s` | `rg[s]·E_g[s]/M[s]`, clipped | source's relayed implied gDNA fraction |
| `f_g^i` | `tg·E_g[i]/M[i]`, clipped | THIS message's implied dst fraction — T5 uses the DST |
| `ρ_μ(s)` | `gp + gn` (graft add, `spl_*_f[sf][s]`) | 0 off-junction ⇒ w_μ = 0 |
| `ρ_ν(s)` | `rp[s] + rn[s]` (relayed continuing RNA) | src frame, consistent with ρ_μ |
| `n_b`, `n_s` | `_n_node[s]` (line 869); `SPN[sf][s]+SNN[sf][s]` | integer counts; 0 ⇒ drop 1/n_s |

## 4. The two limits this must satisfy (assertable, not assumed)

* **Non-junction edge** (`w_μ = 0`): `Var(log k) = Var(λ^s)` ⇒ `prec_g = 1/[(1−f_g^i)²·Var(λ^s)]` — the pure
  **destination-Jacobian composition precision**. So this is NOT graft-only: it corrects EVERY edge's
  composition precision (the old code used the source Jacobian, wrong by 45–132 % per T5's MC). The graft is
  just where `w_μ ≠ 0`. This is why the change *simplifies* — one precision law for all edges, replacing three
  ad-hoc terms.
* **Structural lock** (intergenic, `Var(λ) → 0`): `prec_g → ∞` ⇒ the anchor emits confidently. **Verify
  `lvar_loc ≈ 0` at locked nodes first** (§8 Q3); floor if not.

## 5. Simplifications, stated (each measurable, each a known limit)

1. **`Var(λ^s)` from the LOCAL solve**, not the relayed belief — conservative, well-defined, and exactly what
   is wide for an undecided boundary. Revisit only if measured to matter.
2. **RNA-total → per-strand:** T5 gives `Var(log f_R)`; each active strand receives it (the ± split rides the
   separate tilt axis). No per-strand `k` derived.
3. **`ω = 0`** (owner). Optimistic at high counts; the synthetic battery cannot test it anyway
   (`synthetic_suite_is_poisson_omega_zero`). Land a `CalibrationConfig.omega_spliced` hook, default 0.0, so
   the constant has a home the day real data is available — **no magic number buried in code.**
4. **`f_g^i` via ÷M density mode** ≈ the k-transport `f_g(dst)` — exact for a full composition, first-order for
   a partial message; the clip bounds it.

## 6. Scope fence — what this task does NOT touch (so it stays reviewable)

* **`_relay`'s precision** — it accumulates *context* (fusion weights for the MODE), a separate concern for
  item 3. The measured defect is at the COMBINE, which is what this fixes.
* **`_scan`** and its flags — deleted wholesale in Phase 2, not pruned piecemeal now.
* **The mode / reframe / pin / fuse** — untouched; validated already.
* Everything outside `_unified_solve` ⇒ flag-off (`RIGEL_UNIFIED=0`) byte-identical **by construction**.

## 7. Validation (in order; each gates the next)

1. **Unit:** the two new pure functions, tested against `message_precision_mc.py`'s T1/T2/T5 numbers.
2. **Invariant:** `unified_message_audit.py` — Σ_c f_c ≈ 1 unchanged (mode untouched).
3. **The measured mechanism:** the undecided-boundary probe — decided boundaries must now emit at HIGHER
   precision than undecided ones (the 2.6×-wrong ratio should flip to > 1 the right way).
4. **Worst-scenario dissection** (gdna300_ss_0.50_nrna_none_capture_on): the exon `S+P+R` column should stop
   being dragged below `S+P+G`.
5. **Battery A/B** (RIGEL_UNIFIED=1, vs current-unified AND vs legacy 0.0949), **per condition, stranded AND
   unstranded** (R4).
6. **Flag-off byte-identity** (guard 3,704,635) + `pytest tests/calibration tests/native` + ruff.

## 8. Success + stop condition

* **Success:** unified unstranded exon mwae drops materially, no stranded regression, Σf_c holds, aggregate
  moves toward/past legacy 0.0949 — AND `_transport` is shorter and the math is in a tested pure module.
* **Stop / report if:** the change does NOT close the exon gap → the defect is the *mode* or the *relay
  context*, not this precision; re-diagnose rather than tune. **No magic numbers to rescue a negative result.**

## 9. Bigger picture — the consolidation this STARTS (not finishes)

This is step 1 of making `_unified_solve` an **orchestrator over named pure operations** instead of a 300-line
closure. The end state: `enrichment_frame.py` (→ `message_arithmetic.py`) owns reframe + mode + all precision
laws (graft now, peel item 2, transfer item 3); `_unified_solve` reads as "relay → transport → fuse → solve"
with the arithmetic behind names. Each of items 1–3 moves one more block of math out of the closure. When all
three land, `_scan` is deletable (Phase 2) and `bp_solver.py` roughly halves. **This task must leave that
trajectory strictly monotonic — more math in the pure layer, less in the closure, every time.**

---

## 10. Open questions / issues for the owner

1. **Module home & rename.** Extend `enrichment_frame.py` now, rename to `message_arithmetic.py` in Phase 2 —
   agree? Or rename now (cleaner target name, but touches imports mid-stream)?
2. **`composition_logvar` coexistence.** It (`Var log ρ_tot`, for the reframe *ratio*) and `graft_var_log_k`
   (`Var log k`, for the message *factor*) are genuinely different quantities that both live in the same file.
   Keep both with sharp docstrings, or is there a unification I'm missing? I believe they're irreducibly
   distinct (sum-of-densities vs ratio) — confirm you're comfortable with two similarly-named variance
   functions side by side.
3. **`lvar_loc` at locked nodes.** T5 needs `Var(λ) ≈ 0` at intergenic anchors for them to emit confidently.
   I'll MEASURE this before wiring; if it's not ~0, the honest fix is a structural floor (locked ⇒ Var 0),
   NOT a magic epsilon. Flagging that a floor may be needed.
4. **`ω_spliced` config hook.** Add `CalibrationConfig.omega_spliced = 0.0` now (with the caveat docstring) so
   the value is never a bare literal in the math — or defer the field until real-data ω is wanted? I lean add
   now (no-magic-numbers hygiene).
5. **Relay/combine precision consistency.** I am deliberately fixing only the combine precision (§6). This
   leaves the relay using the old σ²_transfer damping to weight the MODE fusion. If measurement (validation 3)
   shows the relay's mode is now inconsistent with the combine's precision, item 3 pulls forward. Acceptable
   staging, or do you want the relay precision reconsidered in the same pass?
