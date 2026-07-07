# Mature-RNA (spliced) message channel — dissection & design

**Status:** design (2026-07-06). Empirically grounded on the AMBIG-dense 10 Mb benchmark
(`~/Downloads/rigel_runs/ambig_dense_10mb`, condition `gdna300 / ss0.99 / capture-on`). Successor to the
region-236 diagnosis (`region236_ambig_tau_nullspace_intron_ceiling` memory) and the conflation map. Companion
to `dispersion_aware_message_precision.md`.

---

## 1. The hypothesis and the empirical test

**Hypothesis (maintainer):** spliced (mature) fragments are epistemically special — *stranded* (the strand is
fixed by the genomic splice motif, independent of library strand-specificity), *pure RNA*, and *never compete
with gDNA or the opposite strand*. Currently they are lumped into the same RNA density/precision as unspliced
fragments, so mature RNA can bleed into other components across message propagation. Preserving mature RNA is
the highest-value remaining calibration fix — **but this had not been proven on the new benchmark.**

We tested it by dissecting the **highest-error node** in `gdna300/ss0.99/capture-on` and decomposing its
messages into mature (`n_mat = SPs/SNs`, boundary spliced) vs nascent (`n_nasc = fbp·sm`, relayed unspliced
belief) via an env-gated per-edge trace in `bp_solver._scan`.

## 2. What the highest-error node proves

**Region 503 (single-strand `+` exon, node 1007), err = +21,386 — the #1 node.** `true f_g = 0.625`,
solved `f_g = 0.833` (a large gDNA **over-call**). Trace of its incoming RNA+ message (pass 2):

| quantity | value | meaning |
|---|---|---|
| `n_nasc / n_mat` | 243 / **805** | the message is **87 % MATURE** (spliced) |
| density nascent / mature | 1.22 / **8.05** | mature dominates the imputed density |
| target `f_pos` | **0.168** | vs strand-implied / true `f_pos ≈ 0.375` |
| **precision** | **1027.8** | vs ~1.4 on region-236-class AMBIG nodes |
| `s2_edge` | 0.000 | NOT silenced |
| `base_var` | 0.0010 | source boundary is *certain* (`src_var = 0.0000`) |

**The strand alone gives `f_g = 0.626` ≈ truth.** The mature-dominated message, at precision 1027, *overrides
the correct strand solve* and pulls `f_pos` down to 0.168 → the simplex complement inflates `f_g` to 0.833.
The junction spliced reads are **capture-depleted** (probes tile the exon body, not the junction), so the
mature *density under-estimates* the exon's RNA — and it is trusted at extreme precision.

**Cross-regime dissection** (AMBIG Σ|err| by signature) localizes it further:

| condition | mature-involving AMBIG (EX+EX−, EX−IN+, EX+IN−) | pure-nascent (IN+IN−) |
|---|---|---|
| gdna300, nascent, **capture ON** | ~104k (EX−IN+ net **+32k**, over-call) | 1.1k |
| gdna300, nascent, **capture OFF** | ~2.6k (solve fine) | 9.7k |
| gdna300, **no nascent**, capture ON | ~54k (EX+EX− net **−30k**, under-call) | 0.8k |
| **no gDNA**, nascent, capture ON | ~1.2k (low-count only) | 0.5k |

Three facts:
1. **The AMBIG error lives in MATURE-involving nodes.** Pure-nascent IN+IN− 2-D nodes solve well everywhere
   → the 2-D machinery is sound; the defect is specifically mature reliability.
2. **It is capture-driven.** Capture-OFF (mature ≈ unbiased) → mature-AMBIG solve fine; capture-ON → they blow
   up **bidirectionally** (over-call with nascent, under-call without). This confirms the maintainer's point
   that capture *enriches or depletes* junction reads unpredictably — there is no clean capture correction.
3. **Mature is strand-locked but bleeds to gDNA via the complement.** The message stays on the correct strand
   (`SPs` is motif-`+`), but by mis-setting `f_pos` its bias reaches `f_g` through `f_g = 1 − f_pos − f_neg`.

**Verdict: hypothesis CONFIRMED — mature RNA handling is the highest-value fix — with a nuance that reshapes
it.** The defect is *not* that mature is too weak (as region 236 in the old suite suggested). In the benchmark's
dominant regime mature is **too strong and wrong**: a count-rich, source-confident, capture-biased mature
density at precision 1027 overrides the intrinsic strand signal.

## 3. Why the precision is mis-calibrated (root cause)

In `bp_solver._scan` the per-strand RNA message is (RNA-pos shown; `bp_solver.py:1101-1113`):

```
n_nasc = fbp[lsrc]*sm ;  n_mat = SPs[lsrc]
rho    = n_nasc/er + n_mat/esp − rho_mat_dst          # (a) mature POOLED with nascent
mo     = log(max(rho,1/erd)/(md/erd))                 #     one mode
base_var = vbp[lsrc] + 1/(n_nasc+n_mat)               # (b) precision from SOURCE-belief + pooled count
s2_edge  = max((mo − lfp_loc[i])² − (base_var + vp_loc[i]), 0)   # (c) gated on the LOCAL belief
pr     = 1/(base_var + s2_edge)
fbp[i] = combine(local, (mo,pr))                      # (d) mature FOLDED into the relay belief
```

Four coupled defects:
- **(a) Pooling.** Mature (a gDNA-immune, motif-stranded *measurement*) and nascent (a relayed *belief*) share
  one density, one mode, one precision.
- **(b) Precision tracks source belief, not imputation reliability.** A count-rich mature from a *confident*
  boundary gives `base_var ≈ 0` → precision explodes (1027). The magnitude reflects the boundary's certainty
  about *its own fraction*, not whether its junction density predicts the neighbour exon's RNA fraction (the
  capture-geometry question — which is genuinely uncertain).
- **(c) The disagreement gate anchors too loosely.** `s2_edge` is measured against the node's *message-free
  local belief*, whose variance `vp_loc` is inflated by the (uncertain) global gDNA prior. So even a
  strand-**confident** node (503) has a loose `vp_loc` (≥ 0.44), `s2_edge = 0`, and the contradicting mature is
  **not silenced**. The intended behaviour ("a confident strand wins") fails because the anchor is not the
  strand.
- **(d) Relay leak.** Mature is folded into `fbp`, which relays downstream as *nascent* — but mature is spliced
  (skips introns) and does not flow contiguously, so its identity as pure stranded RNA is lost as it propagates.

## 4. Design — a separate, strand-deferring mature channel

The reconciling insight: **a spliced fragment is a high-confidence statement about RNA *presence and strand*,
but a low-confidence statement about the RNA *amount* (its imputed density is capture-biased).** The current
code treats the biased amount as a precise fraction. So the mature channel should establish presence/strand and
**defer to the intrinsic strand signal on the amount** — never override it.

Split the per-strand RNA message in `_scan` into two channels; apply both as independent Gaussians on
`log f_strand` in the solver (it already accepts per-strand messages, `simplex_logodds.py:346-354`):

1. **Nascent channel** — the contiguous unspliced field, unchanged: `rho_n = n_nasc/er − (nascent part of the
   exon→B absorption)`, disagreement-aware precision `1/(vbp[lsrc] + 1/n_nasc + s2_edge)`. This is the only RNA
   that flows contiguously; it alone updates the relay belief `fbp`.

2. **Mature channel** — NEW, fired only when `n_mat > 0` (a junction facing the exon):
   - **Strand-locked** (already): `n_mat = SPs` → RNA-pos only; `SNs` → RNA-neg only. Never gDNA, never the
     opposite strand.
   - **Presence-oriented mode** `mo_m = log(max(n_mat/esp, 1/erd)/(md/erd))` — the mature-implied fraction.
   - **Precision from the count but STRAND-GATED (the crux):** anchor the disagreement on the **strand-only
     belief** `(lf_strand[i], v_strand[i])` — the intrinsic gDNA-vs-RNA signal, *not* the global-inflated local
     belief:
     ```
     base_var_m = 1/n_mat                                   # count sampling only (a measurement)
     s2_m = max((mo_m − lf_strand[i])² − (base_var_m + v_strand[i]), 0)
     pr_m = 1/(base_var_m + s2_m)
     ```
     A **strand-confident** node (503: small `v_strand`) whose strand contradicts the mature → large `s2_m` →
     mature **silenced** → the strand wins (`f_g → 0.626` ≈ truth). A **strand-blind** node (region-236 AMBIG
     plateau: large `v_strand`) → `s2_m ≈ 0` → mature **speaks**, breaking the tilt toward the strand it proves
     RNA on. This is exactly "a confident strand wins, a strand-blind exon snaps" — but anchored on the strand,
     which is where it belongs.
   - **Not relayed:** the mature message updates the destination's *solve* but is **excluded from `fbp`** (the
     nascent-relay belief), so mature RNA never propagates multi-hop as fictitious nascent (defect d).
   - The exon→B mature **absorption** (`rho_mat_dst`) stays on the nascent channel (it removes an exon's mature
     from what crosses into an intron as nascent — a separate, correct mechanism).

`v_strand`/`lf_strand` are already computable — `node_sweep` builds the strand-only local solve `fg_strand`
(`bp_solver.py:1159`); extend it to return the per-strand `f_pos_strand`/`f_neg_strand` and their variances.

### Why this handles both bias directions
- **Capture-DEPLETED mature** (503, over-call): mature under-estimates, contradicts the confident strand →
  silenced → no phantom gDNA.
- **Capture-ENRICHED mature** (none/ON, under-call): mature over-estimates, contradicts the confident strand →
  silenced → no RNA over-attribution.
- **Unbiased mature** (capture-OFF): mature agrees with the strand → passes at count precision, harmless.
- **Strand-blind AMBIG** (region-236 class): strand cannot arbitrate → mature speaks, establishing which strand
  the RNA is on (the intended "more precise imputation").

The design never *assumes* a bias direction or a capture model (both ill-posed) — it lets the intrinsic strand
arbitrate and uses mature only as strand-consistent presence evidence.

## 5. Validation plan (in-process A/B on the benchmark)
Ship only on a clean result across the regimes that isolate the mechanism:
- **capture-ON, mature-involving AMBIG + single-strand exons (503-class):** Σ|err| should collapse (the #1
  node → ~0; EX−IN+ +32k / EX+EX− −30k toward 0).
- **capture-OFF & pure-nascent IN+IN−:** unchanged (already good — the fix must not regress them).
- **no-gDNA guardrail:** no new phantom gDNA on high-count nodes.
- Full 24-condition net-flow + mature-transcript accuracy: flat-or-better.

## 6. Risks & open questions
- **Strand-blind + biased mature** is the irreducible case: on a genuine AMBIG plateau, mature is the only
  RNA-strand evidence, so a capture-biased mature there is unarbitrated. Cross-boundary agreement (the two
  flanking junctions) is a secondary reliability signal to explore if this proves material.
- **`v_strand` for AMBIG nodes** is the 2-D strand posterior's per-strand variance — verify it is large exactly
  when the strand is non-identifying (the plateau) so the gate opens there and closes on single-strand nodes.
- **Nascent channel anchor:** left on the local belief for now; if nascent over-ride surfaces, the same
  strand-anchoring may apply, but the evidence here is mature-specific.
- **Determinism / cost:** one extra per-strand Gaussian + the strand-only per-strand solve; bounded, no new
  grid. Keep the strand-only solve the same call already used for `fg_strand`.
