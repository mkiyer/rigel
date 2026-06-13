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

## 2b. Phase-1 net-flow finding (2026-06-13): a small regression, NOT parity

The 16-condition net-flow A/B (the *primary* metric) shows phase-1 (per-node, no sweep) leaks **+0.5–2 pt
more** gDNA→RNA than production across every gdna300 condition (ss0.99 cap-on 8.7→9.6 / 7.6→9.3; ss0.5
cap-on 20.3→21.2 / 17.8→18.5; cap-off similar), 6/16 worse, 0 better, +8.7 pt total. `gdna_none` stays
FP-free. The earlier "parity" call was from the *pool gDNA fraction* (within 0.7%) — but on net flow that
0.7% is a uniform extra leak. **So phase-1 is not yet shippable; the productionization is gated on reaching
≥ parity on net flow.** The per-node β combine is *meant* to be equivalent to the old `w=I/(I+I₀)` blend, so
this gap is a discrepancy to **close** (diagnose `solve_node` MAP vs the linear blend; the 3-component
strand mixture vs the 2-component posterior; `β`=10 vs `I₀`=10), not an inherent cost.

## 3. The plan (revised priorities — the sweep is the foundation)

**The integration principle (one method, not two signals).** Each node resolves locally by **one elegant
log-likelihood integration** of strand + count on the pie — the likelihood-form analogue of the old
`w=I/(I+I₀)` that worked well (κ=0.99 ⇒ w≈0.96, strand dominates; κ=½ ⇒ w→0, count is all that's left). The
**count term carries a penalty `β`** (the `I₀` successor) so the strand governs where available and the
count is a *fallback tier* that surfaces only when the strand is silent. We do **not** split propagation
into separate strand- and count-signals; there is one integrated belief per node.

**Why the sweep should work now (the key correction).** The first sweep attempt smeared **because it
predated the `β` penalty** (and the full commitment to fractions). The smear was biased count propagating
onto good-strand nodes. With `β` penalizing the count, a propagated count signal arriving at a
confident-strand node has **low likelihood → it cannot tarnish the strand**. The penalty is precisely what
makes propagation safe. So the sweep is no longer expected to smear.

**Priority 1 — the propagation sweep (with `β`).** Turn the bidirectional RTS sweep back on, now that the
count is penalized. Each node resolves locally (strand + count, count penalized), then propagates its
resolved gDNA-density belief; neighbours integrate the received signal **with the count component
penalized** (the penalty is for *propagation/imputation* — a signal-communication event — distinct from a
node using its own count). A good-strand node barely moves; an AMBIG node (no strand) is governed by the
propagated neighbour density over its own biased count. This is the foundational payoff and the reason
parity-or-better is expected. Open: the **exclude-self** (BP message) rule; the imputation penalty
magnitude (may exceed the local `β`).

**Priority 2 — configure `β`.** Tune the count penalty (and, if warranted, a separate, larger imputation
penalty for propagated count) so strand cleanly dominates at κ=0.99 and count takes over at κ=½ — and so
the per-node integration is ≥ the old `w=I/(I+I₀)` (closing §2b).

**Later — `β` sophistication.** 2-level by count-observability (introns MAE 0.005 vs imputed exons 0.47) →
continuous (`observability, var~mean, capture-class`) → derived/calibrated; retires the magic number.

## 4. Productionization — gated on ≥ parity (net flow)

The teardown ships only once the simplex path is **≥ production on net flow** (it is currently behind by
§2b). Sequencing:

1. **Close the per-node gap** (§2b): make `solve_node`(strand + count·`β`) reproduce the old
   `w=I/(I+I₀)` result to ≥ parity on the 16-condition net flow. The two are meant to be equivalent.
2. **Turn on the sweep** (priority 1): propagate the resolved gDNA-density belief; the `β`/imputation
   penalty prevents biased-count smear. Validate it does **not** regress (and ideally improves AMBIG /
   capture-on leak — the expected payoff).
3. **Confirm ≥ parity across all scenarios** (regenerate `gdna_benchmark_5mb`'s 20 conditions — its BAMs
   were cleaned up — + the `quick_3to1_5mb` 16). All unit + golden tests green.
4. **Then** wire to production (flip the `use_propagation` default / make it unconditional), commit/push,
   **teardown** the prior combine (`deconv_regions`; the count-cascade `propagation.py`; dead helpers —
   keeping the strand posterior, the count clue, `region_splice_gdna_frac`, and `deconv_sides` until the
   boundary sides also move to the simplex), and **full code cleanup** (docstrings, dead code, CLAUDE.md /
   docs index, `I₀` vs `β` config).

The teardown is irreversible (recoverable only via git); we do it once, from a validated ≥-parity base —
not on the current small regression.

## 5. Risks / watch-items

- The small consistent `β`-combine −0.005 pool-fraction vs the old `w` blend (the pie trusts the count
  marginally more) — within parity; revisit if phase 2 doesn't absorb it.
- The tiny `gdna_none` FP (+0.004) — a zero-gDNA library should stay at 0; confirm it's negligible on the
  full suite.
- Boundary sides still use the production `deconv_sides` — a hybrid until they move to the simplex; ensure
  the flux transport is consistent.
