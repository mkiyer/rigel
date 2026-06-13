# Calibration simplex deconvolution — current roadmap (authoritative)

**Status:** roadmap, 2026-06-13. The forward plan for the calibration gDNA/RNA deconvolution, consolidating
the simplex pie, the count-trust `β`, and the propagation sweep. Supersedes the sequencing in the older
propagation docs; those remain for theory/history.

## 1. Where we are (phase 1, at production parity, opt-in)

```
scan → calibrate → quant/EM
        │
        └─ per node: strand balance κ; gDNA/RNA strand Beta-Binomial overdispersions;
           node_gdna_density (count clue, spatially imputed: boundary anchors + run-fill);
           splice-junction gDNA-FRACTION upgrade (region_count_frac, count-mean-bias-corrected);
           ── DECONVOLUTION ──
           use_propagation=False (current default): deconv_regions/deconv_sides
               g = w·g_strand + (1−w)·g_count,  w = I/(I+I₀),  I=N·(2κ−1)², I₀=10
           use_propagation=True (phase 1, simplex): deconv_regions_simplex
               POS/NEG: solve_node(own strand [N·(2κ−1)²] + count at fixed trust β) on the pie
               AMBIG/intergenic: count-only (no valid sense — production's w=0 rule)
               boundary sides: deconv_sides (production, unchanged)
```

The simplex path is the **same fusion** as the old combine, re-expressed as a per-node log-likelihood MAP
on the 2-simplex `(f_rna₊, f_rna₋, f_g)`: the strand precision `I=N·(2κ−1)²` vanishes at κ=½, and `β`
(=`count_trust_beta`, the successor to the hard-coded `I₀`) is the explicit count trust. It adds the
no-over-subtraction safety and makes the count trust a first-class, calibratable knob. **Flagship A/B =
parity** (pool gDNA within ~0.7% across ss×capture; the catastrophic flat-average gut is gone).

## 2. What IS and ISN'T propagating (be precise)

| signal | spatially propagated in phase 1? | how |
|---|---|---|
| **count density** (gDNA magnitude) | **yes** | `density_model`: boundary-anchored imputation + bidirectional run-fill across regions (upstream of the deconvolution) |
| **strand** (direction) | **no** | each node uses only its own `u_pos/u_neg` |
| **the simplex RTS sweep** | **NO — built but not wired** | `propagate_simplex`/`_rts_smooth` exist + pass order-independence tests, but `calibrate` uses the per-node `deconv_regions_simplex` |

**Why the RTS sweep is off:** wiring it to propagate the **count density** *smeared* (averaged unrelated
RNA-rich and gDNA-rich exons → regression). The lesson: the sweep must propagate the **strand-derived**
density (unbiased) into the **AMBIG** nodes that lack strand — not the count. That is **phase 2b**.

## 3. The phases (count-trust `β` track + the propagation track)

**Count-trust `β` (per-node misspecification penalty; the integration lever):**
- **Phase 1 — single hard-coded `β`** ✅ (parity; `β` = I₀ successor, documented placeholder).
- **Phase 2 — 2-level `β` by count-observability.** `β_high` on count-observable nodes (introns/intergenic,
  measured MAE 0.005 — trust the count); `β_low` on imputed-across-capture exons (MAE 0.47 — distrust). The
  first expected win over the old fixed `I₀`, targeting the capture-on exon leak.
- **Phase 3 — continuous `β(observability, var~mean, capture-class-crossing)`.**
- **Phase 4 — derived/calibrated `β`** on the benchmark (fused MAE < either signal alone); retires the
  `I₀`/`β` magic number.

**Propagation (carry strand quality to count-only nodes):**
- **Phase 2b — strand→AMBIG propagation.** Turn on the RTS sweep to carry single-strand neighbours'
  strand-derived density (bias −0.15) into AMBIG exons, where a low `β` lets it govern the biased count
  (−0.38). Needs the **exclude-self** (BP message) rule so it does not smear or double-count a node's own
  strand. This is the genuine payoff of the bidirectional sweep.

The two tracks compose: `β` decides *how much* to trust the local count vs the (propagated) strand; the
sweep *delivers* the strand quality to nodes that have none.

## 4. Productionization plan for phase 1 (this effort)

Goal: ship the simplex+β path as the production deconvolution at parity, with a clean codebase, before
resuming phases 2/2b.

1. **Validate across all 20 scenarios** (`gdna_benchmark_5mb`): net `gdna→rna` flow, pool fraction,
   `gdna_none` FP — ON vs the OFF baseline. Require **no glaring regression** (parity within tolerance).
2. **All unit + golden tests green** with the simplex path.
3. **Wire to production**: flip `use_propagation` default → `True` (or make the simplex path
   unconditional). Temporary commit/push.
4. **Teardown the prior framework** (after parity confirmed): retire `deconv_regions`/`deconv_sides`'s
   region/side **combine** that the simplex replaces, the opt-in count-cascade `propagation.py`, and any
   now-dead helpers (`run_fill` if unused, the `strand_deconv` per-node blend). Keep what the simplex still
   uses (the strand posterior, the count clue, `region_splice_gdna_frac`, `deconv_sides` for boundary
   sides — until the sides are also moved to the simplex).
5. **Full code cleanup**: docstrings, dead-code removal, CLAUDE.md/docs index, the config (`I₀` vs `β`).
6. **Resume phases 2 → 2b → 3 → 4** from the clean production base.

**Caveat to decide before teardown:** phase-1 is *parity*, not an improvement, and the **sweep is off**. If
the production foundation should already include the strand→AMBIG propagation (phase 2b), do that *before*
the teardown; otherwise ship per-node phase-1 and add 2b next. The teardown is irreversible — confirm the
scope first.

## 5. Risks / watch-items

- The small consistent `β`-combine −0.005 pool-fraction vs the old `w` blend (the pie trusts the count
  marginally more) — within parity; revisit if phase 2 doesn't absorb it.
- The tiny `gdna_none` FP (+0.004) — a zero-gDNA library should stay at 0; confirm it's negligible on the
  full suite.
- Boundary sides still use the production `deconv_sides` — a hybrid until they move to the simplex; ensure
  the flux transport is consistent.
