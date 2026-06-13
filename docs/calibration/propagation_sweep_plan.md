# The three-term propagation sweep — implementation plan (reference)

**Status:** implementation plan, 2026-06-13. The design for the calibration **propagation sweep** — the
foundational innovation of the simplex method. Builds on the per-node inverse-variance fusion (now at byte
parity with production, `5acb9c8`) and restores the full design we planned (three-term vector + per-
component uncertainty), which my interim "propagate gDNA density only" sketch had thinned. Write/review
this before building.

## 1. The principle (what we planned, restated precisely)

Every node resolves to a **three-term vector** `(f₊, f₋, f_g)` — `+`-strand RNA, `−`-strand RNA, gDNA
fractions of the node's unspliced mass `U`, summing to 1 — **each carrying its own uncertainty**. When a
node communicates with a neighbour, it propagates the **whole three-term belief *with* its per-component
uncertainty**, so the receiver can integrate information *and* uncertainty for all three components
gracefully against what it already knows. One integrated belief per node; messages carry value+uncertainty.

**Why all three terms (not just gDNA).** An AMBIG region (overlapping `+`/`−` transcripts) cannot resolve
its own strands (balanced → strand-blind for `f_g`), but its `+`-RNA is the *same transcript* as the
adjacent single-strand `+` region (often the same exon), and its `−`-RNA the same as the `−` neighbour. So:
propagate `+`-RNA density from the `+` side, `−`-RNA density from the `−` side → reconstruct `f₊, f₋` → the
**residual is gDNA** (`f_g = 1 − f₊ − f₋`). This is the three-signal deconvolution; gDNA-only propagation
cannot do it.

**Reconciling "only gDNA propagates" (the issue-A point):** not wrong, incomplete. The three components
propagate on **different length scales**, encoded by their uncertainty: **gDNA** is smooth (uniform by
enrichment class) → low propagation variance → carries **far**; **per-strand RNA** is valid only within a
transcript/gene body (coverage varies between transcripts, somewhat within) → high propagation variance →
**decays fast**, but enough to reach the adjacent AMBIG node. The per-component uncertainty unifies both;
it is not optional.

## 2. What is built vs. missing

| piece | status |
|---|---|
| per-node pie resolve `(f₊,f₋,f_g)` (`solve_node`) | built — but used only for `f_g`; **per-component uncertainty not exposed** |
| per-node inverse-variance fusion strand+count (`deconv_regions_simplex`) | built (f_g only); **must emit the 3-term + per-component precision** |
| scalar forward–backward smoother (`_rts_smooth`) | built + tested (order-independent) — **scalar; must become 3-term (run per component)** |
| enrichment-class chains (gDNA) | built (`_smooth_density` grouping) |
| **per-strand gene-body chains** (for `+`-RNA / `−`-RNA propagation) | **missing** |
| **3-term message** (densities + per-component precision) | **missing** |
| **component-wise integration** (received ⊕ own, with the imputation penalty) | **missing** |
| **per-component propagation variance** (`var~mean` for gDNA; a coverage/decay model for RNA) | gDNA curve built; **RNA model missing** |
| exclude-self (BP message) rule | **missing** |

So the *machinery* (pie solve, scalar FB, chains, the count `var~mean`) exists; the **three-term
generalization, per-component uncertainty, per-strand chains, and the integration rule** are the new build.

## 3. The design

### 3.1 Node belief
`belief_i = (ρ₊, ρ₋, ρ_g)_i` carried as **densities** (mass per effective bp), each with a **precision**
`(τ₊, τ₋, τ_g)_i`. Densities (not fractions) are the transferable quantity (issue A): gDNA density is
uniform by enrichment class; per-strand RNA density is ~constant within a transcript body. A node converts
a received density to its own fraction via its `U` and effective length, which **caps the smear** at the
node's observed mass.

### 3.2 Local resolve (initialization)
Each node first resolves locally from its **own** strand + count by the inverse-variance fusion
(`deconv_regions_simplex`, β count-penalty), emitting `(ρ₊, ρ₋, ρ_g)` + per-component precision:
- strand-observable POS/NEG: `f_g` from the fusion; the RNA mass goes to the node's own strand
  (`ρ_sense`), the off-strand `ρ` is ~0 with high precision (signature says it's absent).
- AMBIG: `f_g` from the count only (low precision on `f_g`); `ρ₊, ρ₋` **unknown** (low precision) — to be
  filled by propagation.
- intergenic: `ρ_g` only (`f_g=1`), RNA components 0 (certain).
Per-component precision comes from the pie posterior (the grid marginal variances) or, for the
inverse-variance path, the fused precision per component.

### 3.3 Propagation chains (per component)
- **`ρ_g`**: along the **enrichment-class** chain (on-target exon / off-target intron+intergenic) — the
  built grouping.
- **`ρ₊`**: along **contiguous `+`-strand runs** (a gene body proxy); **`ρ₋`**: along `−`-strand runs.
  A node receives a `+`-RNA message only if `+` is active in its signature (AMBIG and POS qualify; NEG /
  intergenic do not).
The three components thus sweep on three (overlapping) chain systems.

### 3.4 The message + integration
Forward and backward sweeps (the order-independent two-pass BP on each chain). A message
`i→j` for component `c` carries `(ρ_c, τ_c)` — value + precision — **excluding** the message `j→i` already
sent (the BP exclude-self rule; prevents double-counting and guarantees the two-pass fixed point). The
receiver integrates per component by **inverse-variance**: `τ' = τ_own + τ_msg`,
`ρ' = (τ_own·ρ_own + τ_msg·ρ_msg)/τ'`, then **renormalizes** `(f₊,f₋,f_g) = (ρ₊,ρ₋,ρ_g)·(L/U)` onto the
simplex (clip ≥0, sum→1). The propagated message's precision is **penalized**: the **imputation penalty**
(harsher than the local β — it is a communicated, possibly-misspecified signal) × a **distance decay**
(precision falls per hop; the per-component `var~mean` / coverage variance accumulates). A confident-strand
node's own precision dwarfs the incoming message → it barely moves (no smear); an AMBIG node's wide own
belief yields to the propagated per-strand RNA densities.

### 3.5 Output
Per-region `f_g` → `CalibrationResult` (unchanged interface). Boundary sides similarly (or via the
existing transport, kept until they move to the sweep).

## 4. ISSUES / INCONSISTENCIES / RISKS (review these)

1. **Uncertainty representation: per-component Gaussian vs the full grid posterior.** Per-component
   `(mean, precision)` is tractable and gives clean inverse-variance integration, but it **ignores the
   simplex coupling** (the three are not independent; `Σ=1`) and is non-Gaussian near edges (`f=0/1`). The
   exact alternative is propagating the **joint grid posterior** (what `solve_node` computes) — exact but
   heavy messages. **Lean:** per-component Gaussian on the densities + renormalize; flag the approximation;
   fall back to grid only if a toy case shows it matters.
2. **Is short-range RNA-density propagation sound?** RNA coverage is non-uniform within a transcript (3′
   bias, isoform structure), so propagating `ρ₊` even one hop is approximate. **Mitigation:** high RNA
   propagation variance (fast decay) so it only assists where the local belief is very wide (AMBIG); a
   confident single-strand node ignores it. **Risk:** if the RNA `var~mean`/decay is mis-set, RNA could
   smear or fail to reach AMBIG. Needs a principled RNA-density variance (a coverage-dispersion model) — we
   have the gDNA `var~mean`; the RNA analogue is **new** and a real open item.
3. **Per-strand gene-body chains.** A contiguous `+`-strand signature run is a *proxy* for a gene body; two
   abutting `+` genes would be one run (RNA density jumps at the gene boundary). **Mitigation:** the RNA
   propagation variance + decay limits cross-gene leakage; or split runs at TSS/TES if the index exposes
   them. Flag: chain definition is approximate.
4. **The AMBIG-to-neighbour transcript identity.** The design assumes the AMBIG region's `+`-RNA equals the
   adjacent `+` region's (same transcript/exon). True when they are the same exon (the overlap case), but
   not if the `+` neighbour is a *different* `+` gene. The strand-run chain + decay handles the common case;
   pathological adjacencies are a residual risk to measure on the toy + flagship.
5. **Exclude-self / convergence.** Two-pass FB is exact on a **tree**; a genomic chain is a tree, so this
   holds **if** the message excludes the target's contribution. With three coupled components +
   renormalization the system is no longer exactly linear-Gaussian → may need 1–2 extra sweeps or light
   damping. Define `tol`; verify order-independence empirically (we have the scalar test; extend to 3-term).
6. **Local vs imputation penalty (β vs imputation).** The local β (count penalty) is settled. The
   **imputation penalty** on propagated components is **new** and must be set so: (a) propagated count/gDNA
   cannot tarnish a good local strand; (b) propagated per-strand RNA *can* rescue an AMBIG node. These pull
   in opposite directions for different components → likely **per-component imputation penalties**. Open
   calibration item.
7. **No-regression guarantee.** Single-strand / simple / no-gDNA nodes must stay at the parity result (the
   local belief is sharp → propagation barely moves it). `gdna_none` must stay FP-free. This is the
   acceptance gate — the sweep must **not** regress the nodes the per-node fusion already solves, and should
   **improve** AMBIG / capture-on.
8. **Performance.** Three component sweeps × chains × two passes × nodes — linear, ~3× the scalar cost; fine.
9. **Boundary sides.** Keep `deconv_sides` (production) for the sides initially; fold them into the 3-term
   sweep later (avoid a hybrid inconsistency in the flux transport — flag).
10. **`solve_node`'s role.** The local resolve currently uses the inverse-variance fusion (not `solve_node`)
    for `f_g`. To emit per-component precision we either (a) compute it from the fusion (RNA precision from
    the strand, gDNA from the fused precision), or (b) revive `solve_node`'s grid posterior for the
    marginal variances. Decide; (a) is lighter and consistent with the parity fix.

## 5. Build phases (incremental, each validated)

1. **3-term local emit:** `deconv_regions_simplex` → `(ρ₊,ρ₋,ρ_g)` + per-component precision (no
   propagation yet; verify `f_g` still byte-parity).
2. **Per-strand chains + the RNA `var~mean`/decay** (the new uncertainty model; issue 2).
3. **3-term FB sweep** with exclude-self + the imputation penalty; order-independence test (3-term).
4. **Integrate + renormalize**; per-region `f_g` → result.
5. **Validate:** toy AMBIG (the `+`/`−`/gDNA reconstruction recovers oracle), flagship net flow (**beats or
   matches** production; the AMBIG/capture-on leak drops), all scenarios, no-regression + `gdna_none` FP.
6. Only then productionize + teardown (per the roadmap).

## 6. Open questions for review (before building)

- **Q1.** Per-component Gaussian messages (lean) vs the joint grid posterior (exact, heavy)? (issue 1)
- **Q2.** The RNA-density propagation variance / decay — how to set it principled (the RNA analogue of the
  gDNA `var~mean`)? This is the biggest new modelling piece. (issue 2)
- **Q3.** Per-component imputation penalties (gDNA vs RNA pull opposite ways)? (issue 6)
- **Q4.** Per-strand chain = contiguous strand run (proxy), or split at gene boundaries if available? (issue 3)
