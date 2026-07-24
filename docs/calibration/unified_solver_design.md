# The unified per-node solver — one mode, not two (design)

**Status:** design, owner-specified, for review before implementation. **Date:** 2026-07-23.
**Owner directive:** consolidate the shift mode and the density mode into ONE solver, on the message algebra
`solve(node, left_message, right_message, priors)`. **Companions:** `enrichment_ratio_generalization.md`
(the bounding lemma §2, the total-density pivot), `enrichment_frame.py` (the `k`-transport primitives),
`message_arithmetic_reconciliation.md` (the mode inventory this retires).

---

## 0. Why — the murk, measured

The solver picks between two message MODES per edge (`use_shift` / `use_shift_g` in `_scan`), plus the
`c_b`/`_var_mat` peel, the `graft`/`absorb` mature terms, the honest clamp, and the gate-equality logic. Six
gated mechanisms in one hot loop. The cost is not just readability — it routes corrections away from where
they are needed. Measured on the verystrong cliff (`e2_edge_trace.py`, one scenario, 6,792 gDNA edges):

```
   edge type     mode                       count
    BND -> exon  SHIFT                        1818     <- the cliff INTO the enriched exon
    BND -> reg   SHIFT                        1578
   exon -> BND   SHIFT                        1818
    reg -> BND   SHIFT                        1146
    reg -> BND   density                       432
   ------------------------------------------------
   SHIFT: 6360 (93.6%)   density: 432 (6.4%)
```

The enrichment reframe (E2) was wired only onto the density path, so it reached **6.4%** of edges and **0** of
the `BND→exon` edges — the exact cliff it exists to correct. The mode split is the bug.

---

## 1. The two modes, precisely — and why they are duals

Both answer one question: given a source's per-component densities `ρ_c` (`c ∈ {g, +, −}`) and the
destination's per-component effective lengths `E_c`, what is the destination's composition `f_c`? Each
component's imputed mass is `M_c = ρ_c·E_c`. **They differ only in the normalizer:**

| | `f_c =` | normalizer | enrichment | component sets |
|---|---|---|---|---|
| **SHIFT** | `M_c / Σ M_c'` (imputed) | the **imputed** total | **cancels** — every `M_c` carries `e(src)` and so does `ΣM`; it divides out. Invariant. | needs **matched sets** — a source missing a component the dst has makes `ΣM` incomplete → asserts that component is absent (`f=0`). Wrong at seam→exon. |
| **DENSITY** | `M_c / M_dst` (observed) | the dst's **observed** mass | **does NOT cancel** — `ρ_c` carries `e(src)`, `M_dst` carries `e(dst)`; wrong by `e(src)/e(dst)` across a cliff. | works with **any set** — `M_dst` absorbs the dst's other components. |

Opposite weaknesses: shift is enrichment-safe but set-rigid; density is set-flexible but enrichment-broken.
The two-mode code chooses per edge and pays for both failure modes at the boundaries (§0).

---

## 2. ⭐ The unifying operation — REFRAME, then DENSITY-normalize by your own mass

The message factor on the destination's composition `f_c` is: **reframe** the source's component density into
the destination's frame by `r = ρ_tot(node)/ρ_tot(src)`, then divide by the destination's **own observed mass**
`M_dst` (its total density) — the density mode, made correct by the reframe:

```
    f_c(node) = ρ_c(src) · r · E_c^node / M_dst          r = ρ_tot(node)/ρ_tot(src) = e(dst)/e(src)
              = e(dst)·ρ̃_c · E_c^node / (e(dst)·Σ ρ̃·E)                     ← e(dst) CANCELS
              = ρ̃_c · E_c^node / Σ_c' ρ̃_c'·E_c'^node       enrichment-free, in the node's own frame
```

**The reframe converts `e(src) → e(dst)`; dividing by `M_dst` (which carries `e(dst)`) cancels it.** This is
literally the owner's rule — *"scale each incoming message to match the node's own total density"* — because
`M_dst/E` **is** the node's own total density. Two consequences, and together they ARE the consolidation:

* **It handles ANY component set.** A partial message — a seam sending gDNA only — gives the correct
  `f_g = ρ_g·E_g/M_dst < 1` (the dst's own mass supplies the RNA it lacks), where the **shift's** imputed-total
  normalizer would assert `f_g = 1` (pure gDNA) and be wrong. The RNA factor from that seam simply carries
  precision 0 and is ignored; the RNA comes from the other message. **This is what retires the shift/density
  split** — the density mode ÷`M_dst` was always the flexible one; it was only *wrong* because it lacked the
  reframe, and the reframe is the fix.
* **It equals the shift where the shift was right.** For a complete, matched, continuous source the reframed
  imputed total equals `M_dst`, so ÷`M_dst` and the shift's ÷(imputed total) coincide (verified algebraically).
  Using the destination's `E_c^node` transfers the content ratio `k = ρ_g/ρ_R` and re-forms `f_g` in the
  destination frame (the `k`-transport of `enrichment_frame.py`) — correct across the region-contained ↔
  boundary-crossing frame change, **not** an `f_g`-copy.

So there is **one operation** — reframe by `r`, then ÷`M_dst` — subsuming both modes. **The mode branch
disappears.**

**`ρ_tot` is LAZY and composition-aware — never a pure-gDNA precompute (owner, 2026-07-23).** `ρ_tot` is the
**sum of component densities, each in its own FL frame**:

```
    ρ_tot = ρ_g + ρ_R⁺ + ρ_R⁻ ,   ρ_g = count_g/E_g^gDNA ,   ρ_R^s = count_s/E_R^RNA   (+ ρ_spliced, §6)
```

RNA-FL for RNA, gDNA-FL for gDNA — the FL distributions differ substantially and that difference **is
information**, not noise. `ρ_tot` is computed **at the last second**, from the best composition available at
solve time (the strand self-solve, both incoming messages, the intron gDNA peel, and the **measured** spliced
RNA at junctions) — **not** precomputed. Treating a node as pure gDNA (`ρ_tot = M/E_g`) is only the **fallback**
when composition is genuinely unknown (an unstranded AMBIG node with no message); the bounding lemma
(`ρ_tot = M·[f_g/E_g + (1−f_g)/E_r]` is a convex blend, so any composition pins it to within the
effective-length ratio — **1.04–1.5× for normal nodes, 4×+ for short regions**) guarantees that fallback is
*safe*, but we use the real composition whenever we have it, which is far more precise than the ≤1.5× bound.
So the sequence per node is **gather messages + self-solve → best composition → `ρ_tot` (lazy) → enrichment
ratios → reframe → solve**, and pass 2 (with the fitted gDNA hyperprior, R1) refines it. Short regions
(`L ≲ fl_mean`) carry no contained mass and are the data-free relays of §5 — they never source an `r`.

---

## 3. The solver — `solve(node, left, right, priors)`

A message is per-component `(ρ_c, p_c)` — a density and a precision — plus the source's `ρ_tot(src)`, all in
the **source's** frame. `p_c > 0` means "I am informed about component c"; `p_c = 0` means "I know nothing".

```
solve(node, msgs=[left, right], priors):

    # (0) REFRAME each message into the node's frame — the ONE operation, applied to EVERY component.
    #     ρ_tot is LAZY (§2): the SUM of component densities in their OWN FL frames, from the BEST composition
    #     available now (self-solve + messages + spliced measurement) — gDNA-FL only as the no-belief fallback.
    #     It is DIRECTIONAL at a junction (§6): WITH ρ_spliced toward the exon/acceptor, WITHOUT toward the intron.
    node.rho_tot[side] = ρ_g + ρ_ν(+) + ρ_ν(-) (+ ρ_spliced if side == acceptor)   # each ÷ its own FL eff-len
    for msg in msgs:
        r = node.rho_tot[msg.side] / msg.rho_tot_src           # enrichment ratio, per side (≤1.5× worst case, §2)
        for c in {g, +, -}:
            msg.rho_c *= r                                      # commensurate; cancels within a matched msg

    # (1) FILTER to the node's active components (owner step 1). gDNA always kept.
    for msg in msgs:
        for s in {+, -}:
            if not node.active[s]:                              # inactive strand
                msg.rho_s, msg.p_s = 0, 0                       # KILLED

    # (2)(3) each message becomes a per-component FACTOR on f_c (the density mode ÷ M_dst, §2), carried into
    # the existing per-node ψ solve. p_c > 0 ⇒ the factor constrains f_c; p_c == 0 ⇒ no factor (ignored).
    for c in node.active_components:                            # {g} ∪ {active strands}
        factor_c = ( log( rho_c(reframed)·E_c^node / M_dst ) , p_c )   # one factor per message
        # The UNION of the two messages (owner step 3) must give some p_c > 0 for c to be identified; if BOTH
        # give p_c == 0 the component is unseen → the ψ solve falls back to the PRIOR (the population baseline).

    # (4) COMBINE — the EXISTING ψ solve reconciles the strand BB likelihood + prior + these factors (kept per
    # the owner's "keep the ψ solve" choice; R4 — a high-precision strand node dominates by construction). The
    # forward/backward α⊗β skeleton and `_comb` are unchanged; only how each factor's (mode, precision) is
    # BUILT changes — reframe + ÷M_dst, replacing the shift/density/c_b/absorb machinery.
    f_c = psi_solve( strand_BB(node), prior, {factor_c from left}, {factor_c from right} )
```

**The node's own belief** (`node.own_c`) is the message-free local solve already computed: the strand
likelihood (the Beta-Binomial tilt — the only intrinsic gDNA/RNA signal, at overdispersed Fisher precision)
plus the node's own mass. It is one more per-component density estimate in the fusion, not a special case.

**The prior** is the population gDNA baseline (`ρ_bg` / the NPMLE composition arm), the architecture's third
information source (`CALIBRATION_ARCHITECTURE.md`): strand + imputation + population. It is what a component
with `p_c* = 0` defers to — the honest handling of an unidentified node, replacing the `f_g = 1` default.

---

## 4. Identifiability — the precision algebra (owner steps 2–3), stated

* `p_c = 0` on a message ⇔ "no information about c" — either the strand is inactive (killed in step 1) or the
  source was uninformed. The two are indistinguishable to the receiver, which is correct: both mean *don't
  learn c from this message*.
* A node **solves component c** iff the **union** of {left, right, own, prior} has `p_c* > 0`. A single
  message missing c is fine when the other supplies it (the seam gives gDNA, the junction gives RNA — the exon
  fuses both). This is why the message need not be a complete composition: **completeness is a property of the
  union, not of one message.**
* When only the prior speaks (`p_c*` from messages + own = 0), the node is honestly unidentified and takes the
  population baseline at the prior's precision — it does not fabricate certainty. This is the
  precision-0 handoff (`§P7` of the restructure plan), now structural rather than a default.

---

## 5. Data-free / tiny nodes — relay, by construction

A node with `L ≲ fl_mean` cannot contain a fragment (§2), so `M_node ≈ 0` ⇒ no `ρ_tot(node)` ⇒ **no frame**.
It cannot reframe, so it **passes each message through unchanged** and contributes `p_c = 0` of its own. The
next framed node reframes across the gap by `ρ_tot(framed)/ρ_tot(framed_across_the_gap)` — the tiny node is
transparent to the ratio, exactly as measured (`enrichment_relay_validation.py`, median residual 0.000 across
relays). This falls out of "no mass ⇒ no frame ⇒ own precision 0 ⇒ adopt what arrived" — **no special case.**

---

## 6. The mature peel — which `ρ_tot`, per side (owner, 2026-07-23)

The mature enters in **two distinct places**, and they must not be conflated (this corrects an earlier draft
that claimed the per-side `ρ_tot` retires `absorb`/`graft` — it does not):

1. **The ENRICHMENT RATIO** — handled inside `ρ_tot`, one-sidedly (the table below). This is new and replaces
   the shift-mode `c_b`.
2. **The RNA MESSAGE routing** — `absorb` (exon→boundary) / `graft` (boundary→exon), a genuine density
   operation via the boundary's measured spliced. This is NECESSARY and stays (⚠ box after the table).

**Structural facts.** A boundary is ONE genomic position, so it supports **at most one splice junction** (one
motif) — the spliced mass is **one-sided**, on the exon/acceptor flank. The boundary's **unspliced** crossing
is inherently **mature-free** (mature splices away and reappears in the spliced channel), so there is nothing
to subtract. `ρ_spliced = spliced_mass / E_R^RNA` is fixed, known, invariant.

**The boundary's `ρ_tot` is DIRECTIONAL**, chosen to component-set-match the neighbor it is compared against:

| enrichment ratio toward | matched set | boundary `ρ_tot` uses |
|---|---|---|
| the **exon / acceptor** (carries mature) | `{g, ν, μ}` | `ρ_unspliced + ρ_spliced`  → `r_exon = ρ_tot(exon)/ρ_tot(bnd, WITH spliced)` |
| the **intron / other side** (mature-free) | `{g, ν}` | `ρ_unspliced` only          → `r_intron = ρ_tot(intron)/ρ_tot(bnd, WITHOUT spliced)` |

`ρ_unspliced = ρ_g + ρ_ν` (gDNA + nascent, the boundary's mature-free crossing). Including the spliced toward
the mature-bearing neighbor and excluding it toward the mature-free one is exactly the `r₁`(TOTAL) / `r₂`
(UNSPLICED) frame pair of `enrichment_ratio_generalization.md` §3.3 — component-set matching, made concrete.

**⚠ The per-side `ρ_tot` fixes the ENRICHMENT RATIO — it does NOT replace the RNA message routing.** A boundary
solves its **unspliced crossing** pie `(f_g, f_ν)`, mature-free. But an exon is a REGION with no spliced
channel, so its RNA density is `nascent + within-exon-mature` **lumped**; the density mode ÷`M_bnd` (mature-free
crossing mass) then over-states `f_ν` because `μ` does not cross. So the mature **routing still stands** (the
§13.0 asymmetry, `absorb` / `graft`), done via the boundary's own **measured** spliced:

* **exon → boundary (PEEL / `absorb`):** the boundary subtracts its measured mature from the incoming lumped
  RNA — `ρ_ν = ρ_R(exon, reframed) − ρ_spliced(bnd)` — because the exon cannot decompose its own RNA and only
  the boundary measures the departing share. A genuine subtraction.
* **boundary → exon (GRAFT):** the boundary adds its measured spliced back as an RNA density — the exon's RNA
  is `ν(crossing) + μ(spliced)`.

Both use the boundary's one measured `ρ_spliced`. **The consolidation retires the SHIFT and `c_b`** (the
shift-mode gDNA mature correction, now subsumed by ÷`M_dst`); the mature **routing** (`absorb`/`graft`) is
NECESSARY physics and stays. `_var_mat` / the honest clamp ride with the routing and are re-derived only when
the precision rework (R2/R3, the next task) lands.

---

## 7. What this retires

| retired | replaced by |
|---|---|
| `use_shift` / `use_shift_g` (the mode branch) | the single reframe + density-mode ÷`M_dst` (§2) |
| `_mode_shift` (+ `Mg/Mp/Mn/_den/comp_fl`) | `_mode_density` on every component, with a reframed `ρ_c` |
| `c_b` / `_var_mat` (the shift-mode gDNA mature correction) | the per-side `ρ_tot` (± spliced) enrichment ratio (§6) |
| the geo-mean / gate-equality (`gate_unequal`) murk | step (1) filter — kill inactive strands at the destination |
| E2's mode-gated reframe (`RIGEL_E2`) | the reframe applied to ALL components, ALL edges |

**KEPT (necessary, not mode logic):** the mature **routing** `absorb`/`graft` via the boundary's measured
spliced (§6 — the RNA departs at the junction; the exon can't self-decompose its RNA); the strand
likelihood + prior in the ψ solve (R4). `_var_mat` / the honest clamp are re-derived with the precision rework
(R2/R3), not here.

`_scan`'s message construction becomes: **reframe → filter → (route mature) → density-mode ÷`M_dst`.** No
per-edge shift/density branch, no `c_b`, no gate murk.

---

## 8. Open questions / risks — to measure, not assume

* **R1 — the reframe uses the node's current `f_g`** (through `ρ_tot(node)`), so the sweep is a fixed point.
  The bounding lemma bounds the error per pass at ≤1.5×; does one forward-backward pass converge, or is a
  second pass needed? *Measure* the per-node `f_g` change between passes.
* **R2 — variance of `r`.** `r` is a ratio of two `ρ_tot` estimates; each carries counting + composition
  uncertainty (`enrichment_frame.composition_logvar`). This must enter the fused precision, or an uncertain
  reframe is applied as exact (the §13.6h mistake). Across a matched message it cancels; across messages it
  does not — so it matters exactly where `r` matters.
* **R3 — σ²_transfer.** Once the mode is enrichment-correct, σ²_transfer's cliff-height damping
  (`(μ_dst−μ_src)²`) is the double-count of R2 and should be retired *with* this change, not before (owner's
  transfer-variance rework). Keep it until the unified solver is measured, then A/B its removal.
* **R4 — stranded.** The current strand-solved nodes are accurate (mwae 0.021); the unified fusion must not
  dilute a strong strand likelihood with messages. The precision algebra should protect this (a strand-anchored
  node's own precision dominates), but it is the primary regression risk and must be A/B'd per condition,
  stranded and unstranded.

---

## 9. Build plan (owner chose: design doc first → then build)

1. **This doc, reviewed.**
2. **`solve_node(node, left, right, priors)`** — a new, self-contained function implementing §3, behind a
   flag, leaving `_scan` as the default. Reuse `enrichment_frame.py` for `ρ_tot`, `k`, the reframe.
3. **Unit tests** on the theorem (§2): a matched single message reproduces the shift bit-for-bit; two
   different-source messages commensurate correctly; a killed strand drops out; `p_c*=0` defers to the prior.
4. **A/B** — `gdna_none` delta, the 32-scenario battery per condition (stranded AND unstranded), the paired
   diff. Flip the default only when it wins, then delete the retired machinery (§7) and regenerate goldens last.

**Standing gates unchanged.** No new constants. The reframe basis is the total density `ρ_tot ≈ M/E_g` with
the `k`-transport re-forming `f_g` in the node's frame (owner's choice), robust by the bounding lemma.

---

## 10. Precise implementation — `_scan_unified` + the two-iteration combine (design turn, 2026-07-23)

Behind `RIGEL_UNIFIED` (default off ⇒ current `_scan`). The forward-backward α/β skeleton and the final ψ solve
stay; what changes is what the relay carries and how the combine builds the imputation factors.

### 10.1 Per-node statics (precompute once, composition-free)
`mass` (M), `E_g`, `E_r`, `E_spl` (node-level, both-face-pooled); `spliced_mass` and the **acceptor face**
(the face carrying `spliced_* > 0`, from `node_geometry` — the junction's exon flank); `free_pos`/`free_neg`
(active strands); `is_boundary`. **`ρ_tot` is NOT precomputed** — it is composition-lazy (§2).

### 10.2 The lazy `ρ_tot` (composition-aware, per side)
`node_total_density(chain, geometry, f_g)` → `(ρ_unspl, ρ_with_spl)` (already built + tested). `f_g` is the
best composition available NOW — the **input belief** at the pass start (init/strand in pass 1; the pass-1
result in pass 2), refined at the combine by the two-iteration below. gDNA-FL fallback only where `f_g` is the
signature default.

### 10.3 The message = per-component `(ρ_g, ρ_+, ρ_-, p_g, p_+, p_-, ρ_tot_sender)` in the SENDER's frame.

### 10.4 `transport(msg, node, side)` — the ONE edge operation
```
r        = ρ_tot(node, side) / msg.ρ_tot_sender          # per-side (WITH spliced toward the acceptor)
ρ_c      = msg.ρ_c · r            for c in {g, +, -}      # reframe ALL components by the same r
# mature routing (§6), boundary endpoints only, via the boundary's MEASURED ρ_spliced:
#   exon → boundary : ρ_ν ← max(ρ_R − ρ_spliced(bnd), 0)     PEEL (absorb)
#   boundary → exon : ρ_R ← ρ_ν(crossing) + ρ_spliced(bnd)   GRAFT
for s in {+, -}: if not node.active[s]: ρ_s, p_s = 0, 0   # FILTER (dst-only)
return (ρ_g, ρ_+, ρ_-, p_g, p_+, p_-)
```

### 10.5 The relay (lightweight — NO ψ solve per hop)
Forward L→R: `fwd_belief[j] = fuse( own_density[j], transport(fwd_belief[j-1], j, left) )`, a precision-weighted
per-component density fuse (`own_density[j]` = the self-solve composition as densities). Data-free node ⇒
`own` precision 0 ⇒ adopts the transported message (pass-through, no frame). Backward R→L symmetric →
`bwd_belief[j]`. Emit `(densities, ρ_tot_sender)`; the receiver reframes.

### 10.6 The combine (per node, the ψ solve — with the two-iteration `ρ_tot`)
```
α = transport(fwd_belief[left(j)],  j, left)             # the left message, in j's frame
β = transport(bwd_belief[right(j)], j, right)            # the right message
for iter in (1, 2):                                      # iter-1 uses the input-belief ρ_tot; iter-2 the both-message one
    imp_c(mode) = density_mode_logfrac(fuse(α_c, β_c).ρ, E_c^j, M_j)   # ÷ M_dst
    imp_c(prec) = fuse(α_c, β_c).p                       # + Var(log r) — DEFERRED to R2 (precision task)
    f_c = ψ_solve(strand_j, prior_j, {imp_g, imp_+, imp_-})           # existing _solve_nodes_logodds_all
    ρ_tot(node) ← node_total_density(f_c)                # recompute lazily from the both-message composition
```
Iter-2 is where `ρ_tot` "waits for both messages" (owner). `n_rho_iters` (=2) is a numerical config knob, not a
model constant — measure whether iter-2 moves anything before committing the cost.

### 10.7 Precision — UNCHANGED this task (R2/R3 next)
`fuse` sums precisions (`_comb`); each message precision keeps the current formula
(`1/(Var(log f)+1/n+σ²_transfer)`). `Var(log r)` (R2) and retiring σ²_transfer (R3) are the **next** task; here
the modes are corrected, the precision is carried as-is, so the A/B is expected ~neutral (the corrected modes
gain weight only after R2/R3). Success = theorem holds + no regression (esp. stranded, R4) + the reframe reaches
`BND→exon` (0→100%, the trace).
