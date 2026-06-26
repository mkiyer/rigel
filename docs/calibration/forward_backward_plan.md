<!-- title: Forward-backward calibration sweep — detailed implementation plan -->
# Forward-backward sweep + the symmetric-refit gating fix — implementation plan

**Status:** designed, not started. Phase 1 (nested loop, GS inner, symmetric var~mean refit) is shipped in the
working tree (`bp_solver.node_sweep`, `config.sweep_max_outer/…`). This plan replaces the **GS inner kernel**
with **forward-backward (FB)** and fixes the symmetric-refit regression the nested loop exposed.

## 0. Why FB (three wins, one change)

Measured facts that motivate this (all from `scripts/debug/gs_convergence_probe.py` + the timing probe):

- **GS does not converge with the entangled refit** (limit cycle); even at fixed var~mean GS needs ~10-12 inner
  passes and leaves a few knife-edge flip-floppers (a node caught between a gDNA-leaning left and an RNA-leaning
  right, because GS updates each node twice per pass from ONE neighbour each, never both at once).
- **GS is slow and unvectorizable.** Calibration is **~128 s** for a 2443-node chain at `max_outer=4`
  (~1 ms/node/pass), because the sweep is a Python loop calling `_local_loglik` on size-1 slices. GS's in-place
  sequential update *forbids* batching (node i needs node i−1's just-written value).
- **GS kills relays** (the high-priority RELAY problem): a low-/zero-count intermediate node transmits its weak
  current belief and erases the upstream signal.

FB fixes all three with one change:
1. **Convergence:** on a chain (our graph is a forest of linear paths), FB is near-exact in one forward + one
   backward pass — ~2-3 outer iterations instead of ~12 inner GS passes.
2. **Speed:** the expensive per-node lattice solve becomes **batchable**. `_local_loglik` is ALREADY vectorized
   over nodes (`(m,)` inputs → `(m,P)` ψ). FB needs only **2 batched lattice solves per outer iteration** (the
   local moments + the final combine), plus two cheap moment-space scans. Expect 10-100× on the sweep.
3. **Relay:** the forward message accumulates the upstream belief, so a thin node *passes it through* instead of
   killing it.

## 1. The graph, components, and gating (recap)

- **Chain** = a forest of linear paths (`chain.left/right`, −1 at reference ends), genomic order
  R-B-R-B-…-R (region/boundary interleave). `chain.order` is the genomic node order.
- **Components** `c ∈ {g, p, n}` (gDNA, RNA-pos, RNA-neg). Each propagates along edges with its OWN gate:
  - `g`: **every** edge (gDNA is genomically continuous).
  - `p`: edges where `free_pos` on BOTH endpoints (a continuous + strand run); RESETS at a + discontinuity.
  - `n`: edges where `free_neg` on BOTH endpoints.
- **Per-node state** (`NodeBelief`): fraction `f_c` + per-component posterior variance `Var(f_c)` (the precision
  state). `prec_c = 1/Var(f_c)` (0 = unsolved, ∞ = locked/forbidden strand).
- **Per-face geometry** (`NodeGeometry`): `mass_{left,right}`, `eff_{gdna,rna}_{left,right}`,
  `spliced_{pos,neg}_{left,right}`. Region presents the same contained geometry both ways; a boundary presents
  its per-side geometry. `node_densities` turns `(f_c, face)` into the message density `ρ = f_c·M_face/E_face`.

## 2. The message (unchanged — reuse `bp_solver._message`)

From source `s` to destination `d` on one edge, per component `c` (gDNA: `spliced=0`):

```
ρ_s       = f_c[s] · M_s^face / E_s^face                      # source facing density
mode      = (ρ_s · E_d^face − spliced_d) / M_d^face           # imputed fraction at d
Var_own_s = (M_s^face / E_s^face)² · Var(f_c[s])              # source's own density uncertainty
prec      = (M_d^face / E_d^face)² / (Var_own_s + σ²_bio(μ) + pois)
```

A message is the pair `(mode, prec)`. `prec=0` ⇒ no-op (no facing flank, or gate off). This is exactly today's
`_message`; FB changes only WHERE the source belief comes from (the forward/backward belief, not the in-place
neighbour).

## 3. The algorithm (one outer iteration)

Two cheap moment-space scans bracket two batched lattice solves.

```
# (A) LOCAL MOMENTS — one batched solve, no messages.  ψ_local = strand mixture + global prior + Jeffreys
#     + spliced floor.  Gives (f_local_c, prec_local_c) per node per component.
psi_local = _local_loglik(ALL nodes, …, gdna_imp=None, rna_imp=None)      # (m, P)
f_local_c, var_local_c = marginals(psi_local)                            # batched _axis_mean/_fg_median/_fg_var

# (B) FORWARD SCAN — genomic order L→R, O(m) cheap arithmetic, per component, gated.
#     α_i = the message INTO node i from its left context.
for i in genomic order:
    L = left[i]
    if L < 0:  α_i,c = (·, prec=0)               # path start: no left message
    else:
        # L's FORWARD belief = local ⊕ α_L  (Gaussian product; left context only, NOT β)
        prec_fwd = prec_local_c[L] + prec(α_L)_c
        f_fwd_c  = (prec_local_c[L]·f_local_c[L] + prec(α_L)_c·mode(α_L)_c) / prec_fwd
        Var_fwd_c = 1/prec_fwd
        # project L(right face) → i(left face) via _message, per component, with the edge gate
        α_i,c = _message(ρ = f_fwd_c·M_L^right/E_L^right, …, Var_own = (M/E)²·Var_fwd_c)   if gate_c(L,i) else prec0

# (C) BACKWARD SCAN — symmetric, genomic order R→L, gives β_i.

# (D) COMBINE + FINAL SOLVE — one batched lattice solve.
#     combined message per node per component = α_i ⊗ β_i (Gaussian product of the two (mode,prec)).
mode_comb_c = (prec(α)_c·mode(α)_c + prec(β)_c·mode(β)_c) / (prec(α)_c + prec(β)_c)
prec_comb_c =  prec(α)_c + prec(β)_c
psi = _local_loglik(ALL nodes, …, gdna_imp=(mode_comb_g,prec_comb_g), rna_imp=(…p, …n))   # (m, P)
f_g = _fg_median(psi); f_pos/f_neg = _axis_mean(psi); var_* = _fg_var(psi)                # batched
```

**Relay** falls out of (B): at a thin node `L` (`prec_local≈0`), the forward belief ≈ `α_L`, so the upstream
message passes through `L` instead of being overwritten by `L`'s weak local evidence.

**Outer loop** (the var~mean fixed point): run (A)-(D), then refit the reliability curves on the final belief
(see §5), recompute the global log-prior, and repeat until the belief stops moving (`sweep_outer_convergence_delta`,
~2-4 iters). Per outer iter: **2 batched lattice solves** (A, D) + **2 O(m) scans** (B, C).

## 4. Accuracy notes (what's exact, what's approximate)

- On a true linear-Gaussian chain, F+B+combine is EXACT. Ours is moment-matched (the message summarizes a
  non-Gaussian lattice belief as `(mode, prec)` per component) — so it's **Expectation-Propagation on a chain**,
  not bit-exact. The FINAL belief (D) uses the full lattice ψ with the combined messages, so only the message
  *propagation* is approximated, not the final per-node posterior.
- If the EP approximation proves loose at bimodal AMBIG nodes, iterate (A)-(D) 2-3× within an outer iter
  (recompute local moments against the combined message — full EP). Start WITHOUT this; add only if measured.
- The scans (B,C) are sequential but O(1)/node. A Python loop is fine at 2443 nodes; for genome scale, vectorize
  per-path as a wavefront or accept a ~seconds scan (the lattice solves dominate either way).

## 5. The symmetric-refit gating fix (REQUIRED — independent of FB)

The nested loop exposed a regression: the **symmetric gDNA var~mean refit on ALL converged nodes is circular**.
Ablation (`scripts/debug/nascent_phantom_ablation.py`, `gdna_none_ss_0.50_nrna_rnd`): the zero-gDNA+nascent
phantom is 98% on EXON nodes (nascent unspliced mass mis-split to gDNA); refitting gDNA on the converged belief
learns that phantom (`ρ_g,global` 0.031→0.037) and amplifies it +18% over the outer iterations; freezing gDNA
recovers ~77% of the growth. Yet the SAME ungated refit HELPS gDNA-present (flagship −11%, capture-off −31%) by
learning the real enriched-exon regime. A recall/precision tradeoff rooted in nascent-vs-gDNA identifiability.

**Resolution — gate to SINGLE-STRAND nodes (`strand_obs`), NOT seed-only and NOT all-nodes.** Seed-only
(intergenic/intron/crossings) was rejected: it discards the single most valuable real-data gDNA-on-exon signal.
In real data MOST genes are unexpressed in any sample → many **exons carry ONLY gDNA** (no RNA), and that is
where hybrid-capture enrichment is strongest and where the intron-exon boundary crossings *underestimate*.
The fit must be allowed onto exons. The compromise that stays non-circular: fit on nodes that are
**definitively strand-solvable** — single-strand nodes (`statics.strand_obs = free_pos ^ free_neg`, already
computed per region AND boundary). These include single-strand exons (an unexpressed one → pure gDNA, the
training gold) but EXCLUDE the hard AMBIG nodes (where gDNA-vs-RNA is not strand-separable and the f_g is the
deconv's own guess — the circular part). So:

- `fit_gdna_varmean`: `live = strand_obs` (was `live = all`). Symmetric criterion: "fit where the component is
  definitively strand-solvable." Single-strand exons train the enriched gDNA-on-exon regime directly.
- `fit_rna_varmean`: the symmetric reading gates each strand to single-strand-of-that-strand
  (`free_pos & strand_obs` / `free_neg & strand_obs`), i.e. drop AMBIG from RNA training too. Measure vs keeping
  the current `free_s` (AMBIG included) — `varmean_policy_compare.py` policies `ss` vs `ss_g`.
- `ρ_global`: no RNA analog (the absolute-level anchor, most circular) — keep anchored (seed-based) for now.

The sim DOES contain unexpressed pure-gDNA exons (168 exon regions, 140 single-strand, 251K gDNA mass on the
gDNA-present flagship), so it CAN validate the recall side — real data just has proportionally far more, so the
sim *understates* the single-strand policy's benefit, it doesn't miss it. **Note the spliced-imputation ~2×
bug** (`spliced_fragment_length_deposit_handoff` / re-audit: one-sided crossing ÷ full eff-len): it
under-subtracts mature from exon unspliced mass → inflates the exon gDNA phantom. Fixing it (separate session)
should remove much of the residual zero-gDNA+nascent phantom independent of the gating. **This policy is
PROVISIONAL** — revisit after the spliced-imputation fix and after FB.

This fix is orthogonal to FB (it's about the refit, not the inference) and lands FIRST (Phase 1), so the
accuracy policy is settled before the algorithm swap.

## 6. Phasing

1. **Refit gating fix (§5)** — gate `fit_gdna_varmean` to gDNA-clean edges; keep `ρ_global` anchored. Re-run the
   nested loop (still GS inner) and confirm the zero-gDNA+nascent regression is gone while gDNA-present holds.
   Small diff, isolates the accuracy fix from the algorithm swap.
2. **FB kernel** — implement §3 behind the existing `node_sweep` signature (same inputs/outputs: it returns
   `(NodeBelief, outer_deltas)`). Replace the inner GS pass loop with local-moments + F-scan + B-scan +
   combine + batched solve. Keep the outer loop + refit (now gated). Drop `sweep_max_passes` semantics to "EP
   inner iterations" (default 1-2).
3. **Vectorize the marginals** — `_fg_median/_axis_mean/_fg_var` over the full `(m,P)` belief in one shot (they
   already are; just call once). Confirm the 10-100× speedup on the timing probe.
4. **Tune** — `sweep_max_outer` (aggregate is stable by iter 1-2; per-node by ~3), EP inner iters.

## 7. Validation

- **Unit:** existing `tests/calibration/` (186) must stay green; add an FB-vs-GS agreement test on a small chain
  (they should converge to the same fixed point at fixed var~mean, within lattice resolution).
- **Convergence:** `scripts/debug/gs_convergence_probe.py` / `outer_convergence_probe.py` — FB outer should
  settle in ~2-3 with far fewer total lattice solves.
- **Speed:** the timing probe — target calibration back to ≲ the old single-loop (~16 s on the 2443-node chain),
  ideally faster.
- **Accuracy:** `scripts/sim/evaluate_suite.py` on `quick_3to1_5mb` + `scripts/debug/precision_benchmark_report.py
  <suite> before_nested` (and `before_precision`). Must: keep the gDNA-present gains (flagship/capture-off), and
  FIX the zero-gDNA+nascent regression (the §5 gating). `scripts/debug/nascent_phantom_ablation.py` to confirm
  the phantom no longer grows with outer iterations.

## 8. Risks / open questions

- **EP looseness at AMBIG bimodal nodes** — mitigate with 2-3 EP inner iters if measured; start at 1.
- **The relay benefit needs a real-data test** — the sim has short genes (short paths); FB's relay advantage
  shows most on long low-count runs (real introns, mappability gaps). Note as a real-data validation item.
- **Genome-scale scan cost** — the O(m) Python scans (B,C) may need per-path vectorization at genome scale;
  the lattice solves are already batched.
- **`ρ_global` / `ê(z)` refit policy** (§5) — anchored vs refit-on-clean; decide empirically with the gating fix.

## 9. Code map

- `bp_solver.node_sweep` — the kernel to replace (currently nested GS).
- `bp_solver._message` — reuse verbatim (the message form).
- `simplex_sweep._local_loglik` — the batched ψ (already `(m,P)`); call once for (A) and once for (D).
- `simplex_sweep._fg_median/_axis_mean/_fg_var` — batched marginals.
- `bp_solver.node_densities` — `(f_c, face) → ρ` message density.
- `bp_solver._edge_varmean` / `fit_rna_varmean` / `fit_gdna_varmean` — the refit; §5 gates the gDNA one.
- `bp_solver._gdna_seed_estimate` — the seed cleanliness mask to reuse as the gDNA-clean edge gate.
- `scripts/debug/{gs_convergence_probe,outer_convergence_probe,outer_meandelta_probe,nascent_phantom_ablation}.py`
  — the convergence + regression instrumentation.
