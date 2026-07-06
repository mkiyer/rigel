# Robust (redescending) gDNA floor prior + message precision — design & test plan

**Status:** design + study (env-gated), authored 2026-07-04 on branch `calib-gdna-accuracy`. Supersedes the
two-sided quadratic pass-1 global prior. Companion to `dispersion_aware_message_precision.md` (the message
side) and `exon_gdna_peak_smoothing_diagnosis.md` / `effective_length_state_and_roadmap.md`.

Project rules that bind every choice: **no magic numbers / new tuned constants**; **messages communicate
densities, not fractions**; principled/derivable, not band-aids. The scale everywhere here is a quantity the
solver already estimates from the data (the node's own variance / the background density) — no new constant.

---

## 1. The diagnosis (proven, not hypothesized)

The pass-1 global gDNA prior is `glob(f_g) = −½·N·r²`, a Gaussian in `r = log f_g − target`,
`target = log(ρ_global·E/M)`. Its **pull** (restoring force) is `−N·r` — **linear in the distance and
unbounded**. It is a Hooke spring: it pulls *hardest* on the nodes *farthest* from the floor, which are exactly
the capture-enriched exons we must not drag. Measured on the flagship (`demo_floor_influence.py`), mean |pull|
by distance-from-floor `|r|`:

| \|r\| bin | n | current (spring) | redescending `r/max(v,r²)` |
|---|---|---|---|
| 0.0–0.5 (at floor) | 516 | 0.15 | 1.65 |
| 4–8 (enriched) | 228 | 7.20 | 0.14 |
| 8+ | 15 | 8.32 | 0.12 |

The spring's influence is **inverted** from what we want: weak at the floor, brutal far away. Node 1085 (enriched
AMBIG exon, r≈7) is dragged with force 7.0. It is also **two-sided** (nodes below the floor are pulled up), which
is the correct *regularizer* role — so the fix is to change the **shape**, not to make it one-sided (a one-sided
floor was considered and rejected: it cannot restore the beneficial downward regularization).

The remedy is a **redescending influence** (robust M-estimator): strong regularization near the floor, pull
vanishing for large residuals. This is the same principle the message precision already uses
(`dispersion_aware_message_precision.md`). The correct scale is **the node's own measurement noise** `v`, so a
sparse/uncertain node is regularized over a wide range while a precise/enriched node releases as soon as it
clears its own noise — this is the "sparse counts shrink to the floor; a precise far observation is trusted"
behaviour, magic-number-free.

---

## 2. `max` vs Cauchy — the two redescending shapes

Both take a squared residual and a base variance and return an inflated variance (→ precision = 1/var):

```
max    (method-of-moments):   v_eff = max(base, r²)      pull = r / max(base, r²)
cauchy (Student-t, ν=1):      v_eff = base + r²          pull = r / (base + r²)
```

`max` is what the message code uses today; it was **not** chosen over Cauchy deliberately — it fell out of the
MoM derivation `σ²_edge = max(resid² − expected, 0)` (a variance estimate can't be negative → `max(·,0)`), and
Cauchy was never in the alternatives table. `max` buys a **dead-zone** (within the noise, full precision) that
keeps the smooth relay undamped; Cauchy is **smooth, proper, and simpler** (drops the `vloc` credit and the
clamp) at the cost of the dead-zone.

| | `max` `1/max(v,r²)` | Cauchy `1/(v+r²)` |
|---|---|---|
| smoothness | kink at r=√v | C∞ (no kink) |
| agreeing input r≈0 | full 1/v (dead-zone) | ≈1/v |
| mild disagreement | full trust to r=√v | attenuates immediately |
| far tail | ~1/r² | ~1/r² (same) |
| proper prior | no (flat tails) | yes (Student-t) |

**Where the difference bites:**
- **Floor prior** — the residual is the *grid variable* (`r = log f_g − target`), so the shape *is* the ψ
  landscape; a kink perturbs the posterior median. **Cauchy's smoothness is a real advantage here.**
- **Messages** — the precision is a *per-edge scalar* (resid fixed at the disagreement), so the kink is harmless
  to grid smoothness. The live question is whether dropping the **dead-zone** damps the long-range relay of
  agreeing messages. **Test empirically across all scenarios.**

---

## 3. The implementation (cleaner + more concise)

One shared primitive expresses the redescending inflation for both uses:

```python
def _robust_var(resid_sq, base):        # base within noise; ~resid² beyond → redescending
    return base + resid_sq              # Cauchy (smooth, proper)
    # max variant: np.maximum(base, resid_sq)   # method-of-moments (dead-zone)
```

### 3.1 Floor prior → `_floor_logprior`
Replace `_global_logprior` (two-sided spring + pooled-ρ_global + `floor_mask` override + pseudocount cap) with:

```python
target = log(clip(ρ_floor · E / M, eps, 1))                 # background-implied f_g, ALL nodes
r      = log(f_g_grid) − target
glob   = −½ · log1p(r² / v_node)         # Cauchy  (or −½·min(r²/v_node, 1) for max)
```
- `ρ_floor` = the genomic background gDNA density (from intergenic+intron regions, `_floor_estimate`) — one
  density for every node; the pooled `ρ_global` and the `floor_mask` special-case disappear.
- `v_node` = the **strand-only** solve variance `Var(log f_g)` (the node's own evidence uncertainty; count-aware,
  magic-free). Obtained from a message-free, global-free local solve.

### 3.2 Message precision → shared primitive in `_scan`
Swap the three `emit_g/p/n` blocks' `max(resid²−(base+vloc),0)` for `_robust_var((mo−lf_loc)², base_var)`
(Cauchy drops the near-cosmetic `vloc` credit + the clamp).

### 3.3 Deletions (the concision win — verified: used only by the old floor)
`_gdna_seed_estimate`, `_fit_seed_varmean`, `ρ_global`/`gdna_vm`/`var_mean`, and the `MonotoneVarMean` usage for
the global spread all become dead once the floor uses only `ρ_floor`. ~50 lines removed; one background density,
one shape, one scale (the node's own variance) everywhere.

---

## 4. Honest prediction

The redescending floor **will** fix capon (enriched exons at r≈7–9 are released — pull ~50× weaker) and should
**partially** recover the capture-off-unstranded regression — not because a floor can tell gDNA from RNA (it
can't), but because enriched-gDNA exons sit far from the background (large r → released) while moderate RNA
over-calls (r≈2–3) stay partly regularized. If capoff recovery is insufficient, that discrimination belongs to
the messages/KDE, not the floor. **The benchmark decides.**

---

## 5. Test plan (env toggles, commit the winner)

Toggles (default = current behaviour, so unset ≡ baseline): `RIGEL_FLOOR_SHAPE ∈ {baseline, off, cauchy, max}`,
`RIGEL_MSG_SHAPE ∈ {baseline, cauchy, max}`.

1. **Floor shape — affected subset first.** All the benchmark action is in the **8 unstranded (ss0.50)
   conditions** (every ss0.99 was inert, |Δ|<720). Run `baseline / off / cauchy / max` on those 8; confirm the
   capon win is kept and measure capoff recovery. Cheap, high signal.
2. **Message shape — all 16.** `cauchy` vs `max`(baseline) across the full suite; watch relay-dominated
   conditions for decay from losing the dead-zone.
3. Winner(s) → `pytest tests/ --update-golden` → `pytest tests/` green → commit + push.

Primary metric: net gDNA↔RNA fragment flow (`net_flow_per_condition.tsv`); secondary: abundance Spearman/MARD/FP.
Study scripts: `scripts/debug/{demo_floor_influence,compare_pass1_global,compare_ab_results}.py`.
