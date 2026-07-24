# Pass-0 message arithmetic — SESSION HANDOFF (2026-07-22)

Read this first. It is self-contained: state, what to do next, how to run everything, and the traps.

---

## 0. STARTER PROMPT (paste this after a context clear)

> We're continuing Rigel's calibration **pass-0** solver on branch `calib-ambig-init-wip` (do **NOT** push to
> main). Read `docs/calibration/pass0_message_arithmetic_handoff.md` in full first, then
> `docs/calibration/message_arithmetic_reconciliation.md` (the MODE/PRECISION reconciliation and the R1–R6
> stage list) and `docs/calibration/prediction_measurement_merge_derivation.md` (item **E**, derived, not yet
> implemented).
>
> **Context:** pass-0's message *mode* arithmetic was structurally wrong across the hybrid-capture enrichment
> cliff — the density mode divides a depleted source density by the enriched destination's observed total, so
> the capture enrichment `e(dst)` never cancels (~600× gDNA under-imputation measured at boundary→exon). That
> is the root cause of "unstranded + capture-ON has never worked". Commit `a48e5970` landed the
> enrichment-invariant mode (composition shift + the mature-dilution term `c_b`), the MASS→COUNT precision fix
> (R1), the retirement of the geo-mean hack (R3), and a toolkit that injects **population-scale priors** into a
> small **full-transcript toy** (a tiny toy cannot fit κ / the NPMLE / ρ_bg, so toy-fit priors are worthless).
>
> **Immediate task: implement item E** — the prediction⊕measurement merge — per §8 of
> `prediction_measurement_merge_derivation.md`. The current code adds precisions (`pr += S`), which is
> replicate-measurement algebra applied to an *additive decomposition*; the fix is a share-weighted
> (delta-method) log-variance combination whose mode is the **sum** with the nascent floored at 0. Landing E
> unblocks **R2**, then **R4** (the σ²_transfer cliff term — the big lever, needs a per-condition A/B on
> stranded AND unstranded arms).
>
> **Constraints:** mantra design→derive→plan→execute, one step at a time, pause and re-evaluate. Elegance >
> accuracy. **No magic numbers** — pause and discuss before adding any constant. Every behavioural change is
> gated by the `gdna_none` phantom guard **as a delta against HEAD** (an absolute number is unreadable — this
> is what caught the +49 % R2 regression). Develop on the grounded toy, validate on the suite. Run everything
> in the activated `rigel` conda env with `OMP_NUM_THREADS=1`. Owner drives commits and sequencing.

---

## 1. Where we are

**Committed:** `a48e5970` on `calib-ambig-init-wip`.

Landed:
* **MODE — enrichment-invariant across the cliff.** intron↔IE-boundary = composition invariance (shift);
  exon↔IE-boundary = shift `± c_b`, `c_b = log1p(S_B/D_B)`. Zero free constants; `e` cancels identically.
* **R1 MASS→COUNT.** `n_unspl_left/right` plumbed through `NodeGeometry`; `_pred_precision` uses the integer
  flux (`Var(log ρ)=1/n`, not `1/mass` — mass is fractional *and* split across nodes by the accumulator).
* **R3.** `rho_g_cross` (the unweighted pre-scan geo-mean) retired — measured **inert**.
* **Toolkit.** `calibrate(injected_priors=…)` + `scripts/debug/toy_inject.py` + the full-transcript toy.

**Deliberately red (gated, not overlooked):**
* **5 golden tests** — behavioural change. Goldens are regenerated **LAST**, only after the per-condition
  benchmark A/B. Do not regenerate them casually.
* **`test_mature_measurement_disagreement_silenced`** — 0.053 vs 0.05 tol. It *quantifies* the probe-depletion
  sensitivity. **Item E is expected to fix this for the right reason** (see §3).

---

## 2. THE NEXT TASK — implement item E

Everything is derived in `prediction_measurement_merge_derivation.md`. The essentials:

**The bug:** `pr += S` adds precisions. Precision-addition is for two observations of the *same* quantity.
Mature and nascent are two **addends** of one quantity (`ρ_r = ρ_m + ρ_n`), so their **variances** combine,
share-weighted in log space:

```
Var(log ρ_r) = w_m²·Var(log ρ_m) + w_n²·Var(log ρ_n) + σ²_xfer ,   w_m = ρ_m/ρ_r, w_n = ρ_n/ρ_r
Var(log ρ_m) = 1/n_m                        MEASURED — Poisson counting ONLY (no composition term)
Var(log ρ_n) = 1/n_u + Var(log f_r^src)     IMPUTED  — counting + the composition solve
```

**Two things to get right:**
1. **The mode is the SUM**, `log(ρ_m + ρ_n)`, with `ρ_n` **floored at 0**. This makes the measured mature a
   *lower bound on RNA that the prediction cannot erase* — which is what structurally prevents the confident
   "no RNA" ⇒ `f_g→1` ⇒ phantom-gDNA pathology.
2. **σ²_xfer applies to BOTH components.** The spliced count is a direct measurement *of the junction* and a
   *prediction about the exon* — a point flux imputed onto a whole region. Nothing here is a direct
   measurement of the exon. Land E with σ²_xfer carried symbolically (start with the existing `s2t`) so the
   **algebra** fix is separable from the **magnitude** question, which R5 measures.

**Do NOT** emit the measurement and prediction as two separate λ-factors — factors on the same log-fraction
axis multiply (precisions add), reintroducing the same error. One message, share-weighted variance.

Implementation sketch + gates: §8 of the derivation doc.

---

## 3. After E — the stage list (roadmap `pass0_roadmap.md` §8)

| # | item | note |
|---|---|---|
| **E** | prediction⊕measurement merge | **NOW** — unblocks R2 |
| R2 | remove the enrichment-cliff proxy from the measurement channel | one-line, once E lands; re-gate |
| **R4** | neutralize the σ²_transfer cliff term `(μ_dst−μ_src)²` | **the big lever.** Behind a switch; **per-condition A/B on stranded AND unstranded** — it currently *protects* stranded data, so it can regress there while fixing unstranded. That asymmetry is the main risk in the whole reconciliation. |
| R5 | measure the residual disagreement after the corrected transfer | decides whether ANY damping survives, and gives σ²_xfer its magnitude |
| R6 | seam derivations: TSS/TES (intergenic↔exon) and exon–exon AMBIG | **completing R6 is what lets the DENSITY MODE be retired entirely** |
| A′ | solve-gate as an **honest precision-0 state** | *Not* the refuted DOF gate, which kept the `f_g=1` init (an all-gDNA lock). The intended paradigm is: always emit; a node that cannot solve stays at precision 0 and the gDNA hyperprior resolves it later. Needs derivation + the hyperprior handoff contract. |
| L | nascent≫mature over-call | **re-measure** — the mode fixes may have addressed it. Do not carry it as open on stale evidence. |

**Ship criterion (unchanged):** pass-0 is HONEST — never confident-wrong, weak where unidentifiable,
`gdna_none` guard held. The gDNA hyperprior owns the residual identifiability floor. Pass-0 does **not** need a
big benchmark number.

---

## 4. How to run everything

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
cd /Users/mkiyer/proj/rigel

# HARD GATE — always A/B against HEAD, never read the absolute number
OMP_NUM_THREADS=1 python scripts/debug/gdna_none_guard.py            # current
git stash push -- src/rigel/calibration/bp_solver.py src/rigel/calibration/node_geometry.py
OMP_NUM_THREADS=1 python scripts/debug/gdna_none_guard.py            # baseline
git stash pop
#   reference: baseline(a48e5970^) 3,741,537 → landed 3,766,743 (+0.7 %); R2 spiked it to 5,586,288 (+49 %)

# per-edge mode/precision trace on the grounded toy (population priors injected, pass-0 only)
OMP_NUM_THREADS=1 python scripts/debug/msg_audit.py

# the mature-dilution identity vs oracle (premise + identity + frame)
OMP_NUM_THREADS=1 python scripts/debug/mature_dilution_check.py                 # one condition + raw dump
OMP_NUM_THREADS=1 python -u scripts/debug/mature_dilution_check.py --grid --fl-sweep > /tmp/md.txt

# suites
OMP_NUM_THREADS=1 python -m pytest tests/calibration/ tests/native/ -q
OMP_NUM_THREADS=1 python -m pytest tests/test_golden_output.py -q     # 5 red BY DESIGN
```

**Toolkit:** `scripts/debug/toy_inject.py` — `extract_priors(cond)` pulls population priors from the cached
`ambig_dense_10mb` scenario (34 conditions in `_selfsolve_cache/`, holds `payload` + `strand_model` + FL, i.e.
exactly `calibrate`'s inputs); `build_toy(...)` generates the full-transcript toy at ambig-matched FL
(RNA=200/gDNA=100) and capture. **`calibrate(injected_priors=…)` is additive — default `None` is
byte-identical.**

---

## 5. Traps — every one of these cost real time

1. **NEVER fit calibration priors on a toy.** κ, the enrichment NPMLE, the intergenic ρ_bg and the strand
   overdispersions need genome-scale (or many-gene) data. A single-gene toy fit κ=0.27 when the truth was 0.50
   and produced a *degenerate flat* NPMLE — and a whole day of "findings" derived from it were worthless.
   Always inject.
2. **The phantom guard is a DELTA.** Absolute numbers are unreadable across scenario sets. Always A/B vs HEAD.
3. **Probe placement is load-bearing.** Junction reads live at exon *edges*; centred probes deplete them ~300×
   (331 → 1 spliced alignments), which silently zeroes `S_B` and the mature correction. `toy_inject._write_probes`
   defaults to full-exon for this reason; probe layout is a sweep dimension (`full` / `centre` / `junction`).
   A probe *over* the junction enriches the boundary *above* the exon — the rule survives (e cancels) but
   under-corrects ~19 %.
4. **`tail`/pipes buffer Python output** — background greps see nothing until exit. Use `python -u > file`.
5. **MASS vs COUNT.** Mass is the numerator for a **density**; the integer count is what a **Poisson variance**
   needs. The accumulator also splits one fragment across several nodes, so mass is not conserved per node.
6. **σ²_transfer is the only damping mechanism left** in production (the old total-density disagreement
   estimator was relocated to `scripts/debug/`). Nothing holds the line if R4 removes it — replace
   deliberately, do not delete casually.

---

## 6. The one-paragraph orientation

Pass-0 deconvolves each node into gDNA vs RNA from three sources only: the strand likelihood (intrinsic), the
cross-node messages (imputation), and the global gDNA prior (out of scope — pass-2). Under hybrid capture the
message *mode* was wrong across the enrichment cliff, so ρ_bg could never enter a transcript body and every
interior node fell to the no-evidence ~0.5; σ²_transfer had been damping those messages precisely *because*
they were wrong, and strand-specific data was strong enough to mask the whole thing. The mode is now
enrichment-invariant. What remains is to make the **precision** honest — E (this merge), then R2, then R4
(retiring the cliff damping that is now cancelling the fix) — and then finish the seam modes (R6) so the
density mode can be retired entirely.
