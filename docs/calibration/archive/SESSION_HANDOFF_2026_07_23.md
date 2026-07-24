# Pass-0 session handoff — the enrichment-ratio solver (2026-07-23)

**Read this first. It is self-contained.** Branch `calib-ambig-init-wip` — do **NOT** push to main.
**Nothing is committed.** The working tree holds the whole session.

---

## 0. STARTER PROMPT (paste after a context reset)

> Continuing Rigel calibration **pass-0** on branch `calib-ambig-init-wip` (do NOT push to main). Read
> `docs/calibration/SESSION_HANDOFF_2026_07_23.md` in full, then
> `docs/calibration/enrichment_ratio_generalization.md` (the framework),
> `docs/calibration/junction_enrichment_scaling.md` (the junction instance + its ⚠ correction), and
> `docs/calibration/message_layer_derivation.md` §12–§13 (the relay, the sink model, the measured diagnoses).
>
> **The task: implement the generalized enrichment-ratio solver**, phase E0 first (see §5 of the handoff).
> Every incoming message is scaled into the destination node's frame by a *measured* enrichment ratio derived
> from total densities, then fused and solved.
>
> **Constraints:** design → implementation plan → execute, one step at a time, pause and re-evaluate.
> Elegance > accuracy. **No magic numbers.** Every behavioural change gated by the `gdna_none` phantom guard
> **as a delta** and the **32-scenario** oracle battery **per condition, stranded AND unstranded**. Run
> everything in the activated `rigel` conda env with `OMP_NUM_THREADS=1`. Owner drives commits and sequencing.
>
> **Standing lesson from the last session: measure before concluding.** Four derivations were stated
> confidently and then refuted by measurement (see §4). Treat any derivation as a hypothesis until a number
> confirms it.

---

## 1. Where we are — the numbers that matter

| instrument | session start | **now** |
|---|---|---|
| 32-scenario oracle battery, mwae | 0.1396 | **0.1361** |
| … correlation | 0.671 | **0.688** |
| `gdna_none` phantom guard | 3,821,731 | **3,704,635** |
| suites (`tests/calibration/ tests/native/`) | 313 pass, 2 xpass | **326 pass, 1 xfail, 2 xpass** |
| goldens | 5 red (pre-existing) | 5 red (unchanged, pre-existing) |

**The pass-0 ship blocker is cleared.** `corr(message precision, |error|)` on boundaries flipped from
**+0.09…+0.67** (confidently wrong) to **negative on nearly every class**, and boundary message precision fell
from 0.62–0.90 to **0.02–0.06** — honestly weak instead of loudly wrong.

### Landed (default ON)
* **P0c** directional mature peel — `mature_dilution`/`_var_mat` are now per-FACE, indexed by the face that
  faces the source. Fixes exon↔exon boundaries that are *also* splice junctions (78/scenario, 19 % of boundary mass).
* **A2** spliced effective length — per-face selector: `boundary_side_eff_length` where the transcript
  continues past the flank, half-triangle where it terminates. Corrects ~280 faces (27 % of spliced-carrying),
  median 3.1×, max 199.6×.
* **A3** `c_b` denominator — per-FACE, not summed-both-faces. Exact against the closed form (r = 1.250 vs
  truth 1.250; the old form gave 0.714).
* **B1a** predicate split — `struct_gdna` (certain) vs `no_data` (unknown, seeded at the derived reference
  variance `Var(λ)=π²`), instead of one `locked` flag asserting "certainly pure gDNA" on 23 % of the chain.
* **B1b** the **gDNA density relay** (`RIGEL_B1B=1`, default) — a running per-node fused
  `(log ρ_g, precision)`; a data-free node passes it through, a low-count node blends. Plus the **intergenic
  wall** (a locus wall: sources its own measured gDNA, sinks what arrives).

### Behind flags, default OFF (measured, kept for A/B)
* `RIGEL_N4A=1` — per-strand RNA relay at identity transport. **Premise refuted** (§4).
* `RIGEL_N4B=1` — routing subtraction + delta-method variance. **Net positive overall but regresses
  stranded** (§4).

---

## 2. The design we converged on — the generalized solver

> For the node being solved, with a left and a right incoming message:
> 1. compute the **total density** of the left source, of the node itself, and of the right source;
> 2. form `r_L = ρ_tot(left)/ρ_tot(node)` and `r_R = ρ_tot(node)/ρ_tot(right)`;
> 3. **scale** each incoming density into the node's frame;
> 4. fuse with the node's belief and **solve**.

Full derivation in `enrichment_ratio_generalization.md`. Two results make it work:

**The bounding lemma.** `ρ_tot = M·[f_g/E_g + (1−f_g)/E_r]` is a convex combination, so
`M/max(E_g,E_r) ≤ ρ_tot ≤ M/min(E_g,E_r)` — **a totally wrong composition still pins total density to within
1.04–1.5×** (measured on the real FL models). Against cliffs of 10²–10³× that is negligible. This is why the
"100 % gDNA at zero precision" fallback is safe. **Exception: short regions (`L ≲ fl_mean`) reach 4×+ and must
be excluded as enrichment references.**

**The r₂/r₁ asymmetry.** In `r₂` (intron↔boundary) the `(k+1)` **cancels exactly** (verified to 5.8e−16 over
20 000 draws), leaving `r₂ = (M_B/M_I)·(k·E_g(I)+E_r(I))/(k·E_g(B)+E_r(B))` — dominated by *observed masses*.
`r₁` (exon↔boundary) does **not** cancel and depends on the exon's own composition — genuinely circular.
**Hence the solve is step-wise:** resolve the boundary from the intron (or the fallback) first, then `r₁`
becomes a measurement rather than an assumption.

**Transport `k = ρ_g/ρ_R`, never `f_g`** — `k` is enrichment-free; `f_g` depends on the node's own
crossing-vs-contained effective lengths.

---

## 3. ⚠ THE OPEN PREREQUISITE (from adversarial review, 44 defects raised → 1 survived)

`junction_enrichment_scaling.md` §3.3 claimed the enrichment estimator was "already validated on all 32
scenarios." **False on two counts, both verified in-repo:**

1. **Wrong estimator.** `enrichment_ratio_census.py` measures `ρ̂` via `node_global_geometry`, which returns
   `eff = eff_gdna_*` — i.e. `M/E_g`, **belief-free and composition-free**. The framework's `ρ_tot` is
   **composition-blended and taken from a solve**. They coincide only at `f_g = 1`.
2. **Wrong scenario set.** The census filters `"0.50" in d.name` — **all 20 conditions unstranded**. The 12
   `ss_0.99` were never censused, and stranded is exactly the arm that regressed in both prior attempts.

What *is* validated is the **taxonomy** (which edges are component-matched, and in which frame). **Phase E0
below closes this gap and is a hard prerequisite.**

---

## 4. Four derivations that measurement REFUTED — do not re-derive them

| claim | refuted by |
|---|---|
| "Silence the exon→boundary edge; it supplies rank 0" | P0 measured **net-negative on the benchmark** (0.1416 vs 0.1396) despite −13.1 % on the `gdna_none` guard. **The guard is one-sided** — in a zero-gDNA library any change lowering `f_g` scores better. Always use the 32-scenario battery. |
| "N4a's no-spliced boundaries are the safe subset for an RNA relay" | **Backwards.** `free_s` says a strand is *present* on both flanks, not that it is the *same transcript*. Measured RNA-density step across those boundaries: median **2.314** (≈10×), p90 7.53. gDNA's is 0.02–0.09. |
| "The peel has no variance; add it and the intron collapse becomes honestly weak" | Implemented; **stranded unmoved** (+0.0047 → +0.0044). I conflated precision with bias — a variance term does not move a mean. |
| "The exon→junction subtraction is ill-conditioned (~86×)" | **Arithmetic error** — I used the exon↔*intron* ratio where the message subtracts at exon↔*boundary*. Measured condition number: **1.1 median**. The 86× amplification is in the **transport across a scale change**, not the subtraction. |

**The pattern: every one was a derivation stated past a number.** Measure first.

---

## 5. THE PHASED IMPLEMENTATION PLAN

Each phase is independently gated. Phases E0–E1 are non-behavioural.

### E0 — close the validation gap — ✅ **DONE, GATE PASSED** (2026-07-23)

`scripts/debug/enrichment_ratio_census.py` rewritten: **all 32 conditions**, stranded and unstranded reported
separately, and **both estimators** side by side (`gDNA-frame M/E_g` vs the framework's composition-blended
`ρ_tot = M·[f_g/E_g + (1−f_g)/E_r]`).

**Capture-off medians under the BLENDED estimator (the gate):**

| | A intgc↔seam | B intron↔junc (unspl) | C exon↔junc (TOTAL) | D exon↔seam *(control)* |
|---|---|---|---|---|
| **unstranded** | −0.090 ✅ | −0.122 ✅ | −0.191 ✅ | **−2.684** (far from 0 ✅) |
| **STRANDED** | −0.103 ✅ | −0.078 ✅ | −0.174 ✅ | **−2.418** (far from 0 ✅) |

**All six measurable cells pass; the unmeasurable control is cleanly separated in both arms.
The framework's premise holds — E1 may proceed.**

Three findings worth carrying forward:

1. **Stranded is now measured, and it is CLEANER than unstranded** — sd 0.236 / 0.403 / 1.514 (A/B/C) versus
   the unstranded arm's much wider spreads. The arm that regressed in every prior attempt is the arm where
   the enrichment estimator is *best* behaved. That is encouraging for E2.
2. **Family A is bit-identical between the two estimators** — as it must be, since intergenic regions and
   seams are structurally `f_g ≡ 1` so the blend degenerates to `M/E_g`. A free internal consistency check.
3. ⚠ **The blended estimator is slightly WORSE-centred than the gDNA-frame one** on B and C (e.g. stranded C:
   +0.050 → −0.174; B: −0.018 → −0.078), because it inherits the solve's composition error. Still inside the
   bounding lemma's tolerance and inside the ±0.25 gate — but it means the blend is **not free**.
   **Carry `gDNA-frame M/E_g` as a viable fallback and A/B the two at E2**; the framework wants the blend on
   theoretical grounds, but today the simpler estimator is marginally better centred.

### E1 — `enrichment_frame.py` (pure module, no solver coupling)
```python
total_density(mass, f_g, E_g, E_r)            -> rho_tot
k_from_belief(f_g, E_g, E_r)                  -> rho_g/rho_R        # enrichment-free
boundary_unspliced_from_k(k, mass, E_g, E_r)  -> (rho_g, rho_R)     # the step-wise resolve
enrichment_ratio(rho_tot_src, rho_tot_dst)    -> (r, var_log_r)
is_valid_reference(...)                       -> bool               # short-region + §5 exclusions
```
**Gate:** unit tests against closed forms (mirror `tests/calibration/test_message_frames.py`); production
byte-identical. **Note the module budget: 31 files against a ~25 target — see §7.**

### E2 — the junction instance, behind a flag
Wire the step-wise solve at splice-junction boundaries: resolve `B`'s unspliced composition from the intron's
`k` (or the 100 %-gDNA-at-zero-precision fallback), form `ρ_tot(B)`, then `r₂` then `r₁`, scale both messages,
fuse, solve. **Refuse** the fallback where §5 of the generalization says it is invalid (AMBIG, exon↔exon,
`S_B=0`, retained-intron) — emit nothing rather than a wrong frame.

**Gate:** `gdna_none` delta; 32-scenario battery per condition; **`paired_node_diff.py` on stranded — the
`intron adj-junction` classes (46 % of the damage) must lose their damage share and stranded must not regress.**

### E3 — variance propagation
`r₁`, `r₂` are ratios of solved quantities. Carry `var_log_r` into the message precision. **Without this we
repeat the §13.6h failure** — an uncertain correction applied as exact.

### E4 — generalize beyond the junction
Per-case total-density recipes, in order: **AMBIG boundaries** (per-strand `k` transfer — *not* covered by the
junction derivation), **exon↔exon-that-is-also-a-junction**, then TSS/TES seams (which should remain the
honest exception — family D is genuinely unmeasurable).

### E5 — cleanup (see §7)

---

## 6. HOW TO RUN EVERYTHING

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
cd /Users/mkiyer/proj/rigel

# HARD GATE — always a DELTA, never the absolute number
OMP_NUM_THREADS=1 python scripts/debug/gdna_none_guard.py                    # current: 3,704,635

# THE benchmark: 32 scenarios vs oracle, per condition (arms are labels; join with --report)
OMP_NUM_THREADS=1 python scripts/debug/pass0_oracle_bench.py --arm mylabel
OMP_NUM_THREADS=1 python scripts/debug/pass0_oracle_bench.py --report --vs base --new mylabel
#   results accumulate in /tmp/pass0_oracle_bench.tsv (delete to reset)

# WHICH NODES regressed, paired between two arms — the diagnostic that finally worked
OMP_NUM_THREADS=1 python scripts/debug/paired_node_diff.py --arm base
RIGEL_N4B=1 OMP_NUM_THREADS=1 python scripts/debug/paired_node_diff.py --arm route
OMP_NUM_THREADS=1 python scripts/debug/paired_node_diff.py --report

# supporting censuses
OMP_NUM_THREADS=1 python scripts/debug/enrichment_ratio_census.py       # the delta taxonomy (E0 edits this)
OMP_NUM_THREADS=1 python scripts/debug/identifiability_census.py        # DOF classes + honesty diagnostic
OMP_NUM_THREADS=1 python scripts/debug/relay_break_census.py            # chain segmentation / dead nodes
OMP_NUM_THREADS=1 python scripts/debug/intergenic_factory_probe.py      # does the factory reach the exon
OMP_NUM_THREADS=1 python scripts/debug/exon_exon_junction_census.py     # the boundary taxonomy
OMP_NUM_THREADS=1 python scripts/debug/stranded_regression_dissect.py   # strand-vs-message interrogation
OMP_NUM_THREADS=1 python scripts/debug/cb_denominator_check.py          # closed-form c_b check

# suites
OMP_NUM_THREADS=1 python -m pytest tests/calibration/ tests/native/ -q   # 326 pass, 1 xfail, 2 xpass
OMP_NUM_THREADS=1 python -m pytest tests/calibration/ -q -W error::RuntimeWarning   # nan/inf hygiene
ruff check src/ tests/ scripts/
```

**Cached substrate:** `~/Downloads/rigel_runs/ambig_dense_10mb` — 32 conditions, `_selfsolve_cache/` so only
`calibrate` re-runs. 20 are `ss_0.50` (unstranded), 12 are `ss_0.99` (stranded).

---

## 7. CODE SURVEY — state and cleanup owed

* `bp_solver.py` is **1,444 lines** and more comment than code in `_scan`. Superseded history (P0/P0b,
  geo-mean, gate-equality, τ-gag) should move **out of the hot loop** into the docs with one-line pointers.
  The loop should read: densities → transport → fuse → modes → fold.
* **31 calibration modules against a ~25 target.** `component_model.py` + `test_component_model.py` are
  **dead**: `edge_vehicle` implements the set-equality rule §12 superseded and has had no production consumer
  since the census was removed. **Delete both** (−2), which pays for `enrichment_frame.py`.
* **Three env flags in production** (`RIGEL_B1B`, `RIGEL_N4A`, `RIGEL_N4B`). `RIGEL_B1B` can retire once the
  pre-B1b path stops being a useful baseline; `N4A` should be deleted outright (premise refuted).
* Extract one `_relay_fuse(run, own, transport, s2t)` helper — the gDNA and per-strand RNA relays are
  identical apart from transport, and writing the fusion repeatedly is how the mature accounting desynchronised.
* **Known-red, deliberate:** 5 golden tests (pre-existing on this branch, verified by stashing) and
  `test_mature_measurement_disagreement_silenced` (xfail, `strict=True`) — the latter documents a real
  coupling: `_boundary_spliced_mass_increment` folds the spliced density into σ²_transfer, so probe depletion
  of the spliced channel perturbs the **gDNA** relay. Different components should not be coupled; the fix
  belongs with E3/Stage C.

---

## 8. TRAPS — each cost real time

1. **The `gdna_none` guard is ONE-SIDED.** In a zero-gDNA library any change lowering `f_g` scores better.
   It is a *hard gate against phantom*, never evidence of accuracy. **Use the 32-scenario battery to decide.**
2. **Always A/B as a DELTA**, never an absolute number.
3. **Per-condition, always.** Aggregates hid a +51 %-relative stranded regression twice.
4. **`mass_unspliced == 0` does not mean "no data"** — a boundary can have zero unspliced crossing and real
   **spliced** mass, and the spliced channel bypasses the unspliced emission gate.
5. **Never fit calibration priors on a toy** (κ, the NPMLE, ρ_bg need genome scale). Use
   `scripts/debug/toy_inject.py`'s injection.
6. **`np.where` evaluates both branches** — `0·∞ = nan`. Substitute a finite value *before* the product.
   Caught by `-W error::RuntimeWarning`; run it.
7. **Do not use a regex to edit `bp_solver.py`.** A greedy `.*?` deleted ~120 lines of `_scan` this session.
   Use targeted `Edit` calls. (Recovered from `git HEAD` + re-applied; verified bit-identical.)
