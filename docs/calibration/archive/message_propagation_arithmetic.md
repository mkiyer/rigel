# Message-Propagation Arithmetic — the density-currency modes

**Status (2026-07-21):** the rigorous, Monte-Carlo–validated derivation of the **message MODE** — how a source
constructs a message and how a recipient interprets it, along the

```
        intron region  ⟷  intron|exon boundary  ⟷  exon region
```

chain, with gDNA + nascent RNA (unspliced) + mature RNA (spliced) and hybrid-capture enrichment. This is the
arithmetic half of the message system: **what value crosses an edge and how the recipient re-frames it.** It
supersedes the mode content scattered across [`cliff_message_derivation.md`](cliff_message_derivation.md),
[`intergenic_boundary_behavior.md`](intergenic_boundary_behavior.md) (the owner's mental model) and
[`message_system_derivation.md`](message_system_derivation.md) §5–§6, folding them into one validated account.

**Scope / what is DEFERRED.** This turn derives **modes only** — the density values and the frame conversions.
Two things are explicitly out of scope and deferred to the next turn:
- **Message precision** (how `pr_c` is computed, and how the spliced-measurement and deconvolution-prediction
  precisions merge). Placeholders are marked `[PRECISION — deferred]`.
- **Solve gating** (the "can I solve? / skip-at-pass-0" decision, [`message_system_derivation.md`](message_system_derivation.md) §6B).

Validation harness: `scripts/debug/chain_mode_mc.py` (promoted from the session prototype). All numbers below are
reproducible with it inside the `rigel` env, `OMP_NUM_THREADS=1`.

---

## 0. The result in one paragraph

**Density is the message currency because it is frame-invariant; the composition fraction is not.** A source emits
per-component **densities** `ρ_c`; the recipient re-anchors them in its own frame (its own effective lengths and
total mass) to get *its* fractions. Along the intron→boundary→exon chain, the two hard phenomena — the capture
**cliff** and the **mature-RNA** composition change — occupy **different edges** and are handled by different, exact
mechanisms: the FL-frame log-odds **shift** on the mature-free `intron⟷boundary` cliff, and the **spliced
measurement** (add going to the exon, subtract coming from it) on the (uniform-capture) cliff-free `boundary⟷exon`
edge. Both are validated to |Δf_g| < 0.003 across fragment-length, density, capture, and ±splice variation. Two
caveats survive adversarial stress (§6a): the composition (f_g) transfer holds even under realistic *proportional*
capture, but the **mature** add/subtract acquires a differential-capture bias there (∝ mature fraction), and short
regions (`L < fl_mean`) hit an information wall the solver already handles by deference.

---

## 1. The belief object (verified)

A node owns a **fixed total unspliced mass** `M` (a count). Its unknown is the **composition** of that mass —
`(f_g, f_+, f_−)`, `Σ f = 1`, **2 DOF**, carried in the two natural coordinates:

- **λ = logit(f_g)** — the gDNA-vs-RNA-total level (`f_g = σ(λ)`);
- **θ (tilt)** — the +/− strand split, `f_± = (1−f_g)(1 ± sin θ)/2`.

The two precisions are `(pr_λ, pr_θ)`; `θ` is locked (`pr_θ = ∞`) for single-strand nodes, free for AMBIG. The
**gate** `(g, +, −)` (gDNA always on; each RNA strand on/off by structure) sets the DOF: `|active|−1` free axes.
Belief is stored as fractions **or** counts (interconvertible, `C_c = f_c M`); **density is not stored** — it is
derived on demand (§2). This section restates [`message_system_derivation.md`](message_system_derivation.md) §1.

---

## 2. Density is the frame-invariant currency (Theorem 1)

**The three frames.** The same physical location is measured differently by different node types:

| frame | what it counts | component eff-length `E_c` (`effective_length.py`) |
|---|---|---|
| **region-contained** (R) | fragments wholly inside a length-`L` region | `region_eff_length(L, FL_c) = E_f[max(0, L−ℓ+1)]` |
| **boundary-crossing** (B) | fragments straddling a point | `boundary_eff_length(FL_c) = fl_mean(FL_c)` |
| **spliced one-sided** (S) | spliced fragments crossing a junction, exon flank only | `spliced_side_eff_length(FL_r, R) = E_f[min(ℓ,R)²/2ℓ]` |

**Density** for component `c` at a node: `ρ_c = C_c / E_c` (component mass over that component's frame eff-length).

> **Theorem 1 (frame-invariance).** For a component present at a node, `ρ_c = C_c / E_c` equals the local
> **molecular density × capture** `d_c·κ`, *independent of the measuring frame*.
>
> *Why.* The mass is `C_c = d_c · κ · E_c`, where `E_c` is precisely the geometric **opportunity** (the count of
> fragment start-positions that produce the counted event in that frame) and `κ` is the local capture enrichment.
> Dividing it back out leaves `ρ_c = d_c κ`, a frame-free quantity. The FL distribution enters only through `E_c`,
> so gDNA and RNA (different FLs) have different `E_c` but each recovers its own `d_c κ`.

**MC confirmation.** The *same* mature RNA measured three ways agrees to ≈0.3–0.4 %:

| regime | ρ_mature (spliced frame, boundary) | ρ_mature (contained frame, exon) | rel. Δ |
|---|---|---|---|
| gDNA=RNA gauss200 | 1003.4 | 999.8 | 0.36 % |
| gDNA300>RNA150 | 1003.4 | 999.8 | 0.36 % |
| capture E=300 | 301 014.9 | 299 932.6 | 0.36 % |

The residual is MC noise plus the discrete `+1` in `region_eff_length`; it does **not** grow with the capture
scale (E=1 → 30 → 300 all give 0.36 %), confirming the invariance is exact, not a small-signal coincidence.

**Corollary — why the fraction shifts.** `f_c = ρ_c E_c / Σ_c' ρ_c' E_c'`. The `E_c` weights differ by frame,
so the *same* densities give *different* fractions in the contained vs crossing frames. Density is conserved;
composition is a frame-dependent projection of it. This is the entire origin of the "cliff shift" (§4).

---

## 3. The message and the recipient rule

**Message = per-component `(ρ_c, pr_c, gate_c)`** (the owner's format). The source emits **densities** (what it
measured) + **honest precisions** `[PRECISION — deferred]`; the **gate** is a structural constant read from a
vector, not on the wire. The source does **not** compute the recipient's composition and does **not** gate
emission — a zero-precision message is ignored, so "not emitting" ≡ "emitting pr=0" (emission may be skipped only
as a bit-identical performance optimization).

**Recipient rule (equal gates).** Given source densities `{ρ_g, ρ_+, ρ_−}` and its **own** eff-lengths
`{E_g^dst, E_r^dst}`, the recipient imputes its per-component masses `M_c = ρ_c E_c^dst` and normalizes by the
**imputed** total:

```
   f_g^dst = ρ_g E_g^dst / D,       D = ρ_g E_g^dst + (ρ_+ + ρ_−) E_r^dst
```

The capture scale `k` (composition-invariance: `ρ_c^dst = k·ρ_c^src`) **cancels** in the ratio ⇒ the mode is
**cliff-invariant**. This is the composition SHIFT, done recipient-side from densities.

> **Corrects the owner's sketch.** The sketch wrote `dst_rho_g = src_f_g · dst_count_g / dst_eff_len`, but the
> dst's per-component **count** `dst_count_g` is the unknown. The correct recipient-side form uses the dst's
> **total** mass and per-component **eff-lengths** — the normalized form above. The composition (a ratio)
> transfers; the recipient re-anchors the scale to its own total. (Equivalent to the log-odds shift in §4.)

---

## 4. The chain decomposition (Theorem 2 — the load-bearing structure)

> **Theorem 2.** Along `intron ⟷ boundary ⟷ exon`, the capture **cliff** and the **mature** composition change
> occupy **different edges**:
>
> - **`intron ⟷ boundary`** — a capture cliff (deep intron un-captured; the boundary's crossing fragments
>   overlap the exon ⇒ captured), but **mature-free**: a fragment crossing the intron|exon junction contains
>   intronic sequence ⇒ it is **unspliced** ⇒ gDNA or **nascent**, *never mature*. So the only composition change
>   is the FL-frame + capture scale — handled **exactly** by the shift (§4a).
> - **`boundary ⟷ exon`** — **no capture cliff under uniform (binary-overlap) capture** (both are exon-side), but
>   **mature-bearing**: the exon contains mature RNA the boundary's unspliced crossing does not, and the
>   boundary's **spliced** mass measures exactly that mature. So the only composition change is mature
>   reconciliation — handled **exactly** by the spliced measurement (§4b).

This separation is what makes the chain tractable: neither edge carries both hard phenomena at once.

> **⚠ Capture-model caveat (found by the proportional-capture stress, §6a).** The `boundary⟷exon` "no cliff" is
> exact only when capture is a **per-node scalar** (binary overlap). Under realistic **length-proportional**
> capture, a boundary-crossing fragment lies only *partly* in the exon ⇒ it is **under-enriched** vs a fully
> exon-contained fragment (≈48.5 % of the exon's enrichment at E=30, → 50 % asymptotically; ≈6 % for a narrow
> central probe). **Crucially, the *composition* transfer survives** this cliff: gDNA and nascent are **both**
> unspliced/gene-body and suffer the *identical* per-fragment cliff, so it **cancels in every f_g ratio** —
> C1 stays robust (|Δ|≤0.0009) and C3/C4 are exact when there is no mature (d_mat=0) *even with a 50 % cliff
> present*. What does **not** cancel is the **differential** cliff between the unspliced channel (~48.5 %) and the
> one-sided spliced/mature channel (~32 %, different overlap geometry): the mature add/subtract (§4b) mixes two
> differently-captured frames, biasing C4 by up to 0.11 and C3 by up to 0.04, **monotone in the mature fraction,
> zero at d_mat=0**. The mature reconciliation therefore needs a **capture-frame scale correction** under
> proportional capture — flagged in §4b and deferred (it is capture-model-specific and precision-adjacent).

### 4a. `intron ⟷ boundary` — the FL-frame log-odds shift (mature-free cliff)

Both endpoints admit `(g, +)` (or `(g,−)`); both see only gDNA + nascent. The shift converts the intron's
contained-frame `f_g` into the boundary's crossing-frame `f_g`:

```
   λ_bnd = λ_intron + [ log(E_g^B / E_g^R) − log(E_r^B / E_r^R) ]
```

with `E_g^R = region_eff_length(L_intron, FL_g)`, `E_g^B = fl_mean(FL_g)`, similarly for RNA. The capture
enrichment cancels (the shift depends only on eff-length **ratios**), and the shift is nonzero exactly because
gDNA and RNA have different FLs. **MC:** |Δf_g| ≤ 0.0013 across all FL pairs, and identical under E = 1/30/300
(cliff-invariant). This reproduces [`cliff_message_derivation.md`](cliff_message_derivation.md).

### 4b. `boundary ⟷ exon` — mature reconciliation via the spliced measurement

The boundary's **spliced** mass gives `ρ_mature = S / E_+^S` directly (Theorem 1 — same capture scale as the
exon, since both are exon-side; Theorem 2 — no cliff to cross). Then:

- **`boundary → exon` (add).** The exon's RNA is nascent **+** mature, so the RNA density the boundary reports is
  `ρ_r = ρ_nasc + ρ_mature`. Shift with `(ρ_g, ρ_nasc + ρ_mature)`:
  ```
     f_g^exon = ρ_g E_g^exon / ( ρ_g E_g^exon + (ρ_nasc + ρ_mature) E_r^exon )
  ```
  **MC:** |Δf_g| ≤ 0.0010.

- **`exon → boundary` (subtract — "the hardest node").** The boundary's unspliced crossing is nascent-*only*
  (mature does not cross the intron|exon junction). It recovers the nascent density by subtracting its **own**
  spliced from the exon's RNA message, then shifts into the crossing frame:
  ```
     ρ_nasc = ρ_r^exon − ρ_mature^bnd            (exact: both at the same exon-side scale, Theorem 2)
     f_g^bnd = ρ_g^exon E_g^B / ( ρ_g^exon E_g^B + ρ_nasc E_r^B )
  ```
  **MC:** |Δf_g| ≤ 0.0024.

**Why the exon→boundary message is worth the trouble.** The boundary *already* gets its unspliced composition from
the intron side (§4a) — but the intron is un-captured, sparse, and (in a depleted library) a weak gDNA estimate.
The exon is capture-**enriched** ⇒ many more counts ⇒ a far more precise gDNA density. The absorption lets the
boundary borrow the exon's statistical power. `[PRECISION — deferred: how the intron-side and exon-side estimates
combine, and how the spliced-measurement vs nascent-prediction precisions merge.]`

> **Two caveats on the subtraction (found by the workflow, §6a), both real but out of this turn's scope:**
> 1. **Capture-frame mismatch (proportional capture).** `ρ_mature^bnd` (spliced, one-sided geometry) and
>    `ρ_r^exon` (exon-contained geometry) are captured at *different* per-fragment rates under proportional
>    capture, so the raw subtraction mixes frames (§4-caveat). A **capture-scale correction** on the mature term
>    is needed there; under uniform capture it is exact.
> 2. **Conditioning when mature ≫ nascent.** `ρ_nasc = ρ_r^exon − ρ_mature` is a difference of two near-equal
>    large densities, so per-density noise is amplified by ≈`(1 + ρ_mat/ρ_nasc)` into the small recovered
>    `ρ_nasc`. At `ρ_mat/ρ_nasc = 50` the recovered nascent is ill-conditioned (a **precision**, not a mode,
>    issue — the boundary should down-weight the exon→boundary nascent estimate accordingly). `[PRECISION — deferred]`

---

## 5. Unequal gates & the "gDNA-always" question (Theorem 4)

> **Theorem 4.** The gDNA density transfers by the **same** frame/scale rules as any component — it is not
> algebraically special in *how* it transfers. What is special: gDNA is **gated on at every node** (present
> everywhere genomically), so a gDNA density is *always* an available currency, and the *only* component
> guaranteed present across an **unequal-gate** edge.

Consequences:
- **Equal gates** (the whole intron/boundary/exon chain above): the shift uses **all** components; gDNA is simply
  one of them. There is nothing extra to prove — the "composition shift works for gDNA" because it works for
  every component.
- **Unequal gates** (the canonical `intergenic-exon seam (g,0,0) → exon (g,+,0)`): the RNA source densities are
  absent, so the normalizer `D` is **incomplete** ⇒ the recipient's **fraction** `f_g` *cannot* be computed from
  the message alone. But the gDNA **density** still transfers, as a **one-sided lower anchor**:
  ```
     ρ_g^dst ≥ ρ_g^src     (capture only ADDS gDNA to the exon that the un-captured crossing background misses)
     ⇒ f_g^dst ≥ ρ_g^src E_g^dst / M_dst      (a weak LOWER bound on λ, never a point estimate, never a down-pull)
  ```

**So the answer to "is it a new theorem that composition shift always works for gDNA?"** — it is **true for the
density, false for the fraction**. The gDNA density is always transferable (that is Theorem 1 applied to the
ever-present component); the gDNA *fraction* at the destination still needs the RNA normalizer, which is exactly
what an unequal gate withholds. The unequal-gate case is therefore a **lower-bound anchor**, and the destination
composition is genuinely **unidentifiable** without the gDNA prior — [`message_system_derivation.md`](message_system_derivation.md) §6/§9.
(This section is the mode; the lower-anchor's ψ form and precision are `[PRECISION — deferred]` and belong to the
solve-gating turn.)

---

## 6. MC validation summary

`scripts/debug/chain_mode_mc.py`, `Nbase = 8·10⁶`. Worst |Δf_g| per claim over the grid {FL pairs: gauss/gamma,
symmetric & asymmetric; capture E ∈ 1/30/300; densities d_g:d_nasc:d_mat spanning gDNA-rich…RNA-rich; ±splice}:

| claim | edge | mechanism | worst |Δf_g| |
|---|---|---|---|
| **C1** | intron → boundary | FL-frame shift (§4a) | **0.0013** |
| **C2** | boundary ↔ exon | density frame-invariance (§2, Theorem 1) | **0.36 %** (rel.) |
| **C3** | boundary → exon | + spliced (§4b) | **0.0010** |
| **C4** | exon → boundary | − spliced absorption (§4b) | **0.0024** |

Cliff-invariance is exact: E = 1/30/300 give bit-comparable fractions with densities scaled by E. The `d_mat = 0`
("no splice junction") cases collapse C3/C4 to the pure-nascent shift correctly (|Δ| ≤ 0.0001).

### 6a. Adversarial stress (workflow `chain-mode-arith-validation`, 4 independent agents)

| agent | verdict | finding |
|---|---|---|
| **harness audit** (skeptic) | HONEST / non-circular | Predictions use only observable message inputs; truths are independent raw-count fractions (re-derived, Δ=0); eff-lengths imported from production `effective_length.py`, no shadow reimpl.; C4 re-derived from the generative model and matches the harness line. The validation does not beg the question. |
| **extreme regimes** (125-cell d_g:d_nasc:d_mat × FL grid, Nbase→96M) | ROBUST | Every |Δ|>0.01 traced to finite-count MC noise (∝1/√N, no systematic floor). C1≤0.008, C3≤0.002 everywhere incl. strongly asymmetric / bimodal / uniform FL and near-pure-gDNA / near-pure-RNA. |
| **short flanks** (Li/Le/Pe ∈ {80…1200}) | BREAKS_AT_EDGE | Two nameable edges (below). |
| **proportional capture** (per-fragment overlap weight) | REFUTED (a real caveat) | The `boundary⟷exon` "no cliff" is binary-capture-specific (below). |

**Three validity boundaries** established (none is a derivation gap in the *composition* transfer):

1. **Identifiability WALL — `L < fl_mean`** (by design, not a gap). Push a region below the fragment-length scale
   and its **contained population → 0**, so its density is `0/0` and a point-mode prediction blows up (C1 intron
   wall at Li=80/fl_mean_g=300 → 0.615; C4 exon wall → 0.547). This is **not** fixable by any eff-length (`+1`
   and `+0` blow up identically) — it is exactly the *honest-imprecision* regime `effective_length.py` documents:
   `region_eff` returns a tiny eff ⇒ `Var(log ρ)=1/n` explodes ⇒ the node emits an uninformative message and
   **defers to its neighbours**. The point-prediction "breaks" but the solver *knows* it cannot measure that
   node (this is the solve-gating concern — deferred).

2. **Discretization-frame inconsistency — `region (+1)` vs `spliced/boundary (continuous)`** (fixable, code-relevant).
   `region_eff_length` uses the **discrete** `E[max(0, L−ℓ+1)]` (the `+1`, correct for integer-bp accumulator
   data) while `boundary_eff_length = fl_mean` and `spliced_side_eff_length` are **continuous** integrals. At
   short flanks (`L ≈ fl_mean`) the differential `+1` is a several-percent frame bias that does **not** cancel in
   the gDNA−RNA shift (C1 → 0.004), the crossing-vs-contained comparison (C2 ρ_g → 0.026), and worst in the **C4
   absorption** (→ 0.028 at fl_mean=200/Le=150), because it subtracts a mature density measured in the
   *continuous* spliced frame from one measured in the *discrete* region frame. Under a matched frame (both `+0`
   for the continuous-start MC, or both `+1`-consistent for real data) all collapse to MC noise. **Actionable:**
   ensure the spliced/boundary eff-lengths share `region_eff_length`'s discretization convention (connects to the
   standing note in `effective_length.py` that `spliced_side_eff_length` is "deferred, not cosmetic").

3. **Proportional capture reintroduces a cliff on `boundary⟷exon` — but only the MATURE channel is hurt.** See
   the §4 capture-model caveat: the *composition* transfer survives (the unspliced cliff self-cancels in f_g
   ratios; C1 |Δ|≤0.0009, C3/C4 exact at d_mat=0 with a 50 % cliff present), while the **differential** cliff
   between the unspliced (~48.5 %) and one-sided spliced (~32 %) channels biases the mature add/subtract — C4 up
   to 0.11 (0.17 under a narrow probe), C3 up to 0.04, **monotone in mature fraction**. This is the one place the
   arithmetic needs a follow-up (a mature capture-scale correction) — flagged, precision-adjacent, deferred.
   *Reproducible in the shipped harness:* `simulate(..., capture_mode="proportional")` — at E=30 it measures
   `ρ_g^bnd/ρ_g^exon = 0.515` (the cliff is real) yet C1 |Δ|=0.0003 and C3-at-`d_mat=0` |Δ|=0.0001 (the composition
   transfer survives it), while C3-with-mature |Δ|=0.030.

**Bottom line.** The density-currency composition math (C1–C4) is honest and holds to MC noise across the entire
mainstream and most of the short-flank grid. The three boundaries are (1) a *known, handled* information wall,
(2) a *fixable* eff-length discretization-consistency item, and (3) a *real refinement* — the mature reconciliation
needs a capture-scale correction under proportional capture. None disturbs the core result: **density is the
frame-invariant currency, and the composition transfer is exact.**

---

## 7. The send/receive algorithm (mode only)

```
SOURCE emits, per component c present at the node:
    ρ_g   = f_g · M / E_g^src                         # gDNA density (crossing or contained frame)
    ρ_+   = f_+ · M / E_r^src   (+ spliced: + S_+/E_+^S,  − absorbed source-exon mature)
    ρ_−   = f_− · M / E_r^src   (+ spliced: + S_−/E_−^S,  − absorbed source-exon mature)
    gate  = structural (g:on; ± per free_s continuity)
    pr_c  = [PRECISION — deferred]

RECIPIENT interprets, per incoming message:
    if gate_src == gate_dst:                           # EQUAL gates → composition shift (§3, §4)
        M_c   = ρ_c · E_c^dst                          #   impute dst masses in the dst frame
        f_g^dst from the normalized shift              #   (mature already folded into ρ_± by add/subtract, §4b)
    else:                                              # UNEQUAL gates → gDNA-only lower anchor (§5)
        transferable = { gDNA }                        #   RNA components unconstrained
        f_g^dst ≥ ρ_g · E_g^dst / M_dst                #   one-sided; never a point estimate / down-pull
    fold onto the running (λ, θ) belief                #   [PRECISION — deferred: the EP moment-match weights]
```

Note the split is on **gate-equality** (structural), not on exon-membership (the current `use_shift = not
exon@either` proxy in `bp_solver`). The FL-frame **shift** is the correct mode on *every* equal-gate edge — cliff
or not — with the mature add/subtract layered onto exon-facing flanks. There is **no separate "density mode"** in
the arithmetic: C3/C4 above are the *shift* (normalize by the imputed total), not the observed-`md` density anchor.
This is the mode-level footing for [`message_system_derivation.md`](message_system_derivation.md) §6A ("retire the
density mode"); whether shift-everywhere is also *robust* in the solver is a **precision** question (deferred).

---

## 8. Relation to the current code (`bp_solver._scan`)

- **Confirmed correct.** The `rho_pos = fp_s·sm/_er + SPs/esp − absorb_p` construction ([bp_solver.py:609](../../src/rigel/calibration/bp_solver.py))
  is exactly §4b's add/subtract: `+SPs/esp` adds the dst-exon mature (boundary→exon), `−absorb_p` removes the
  source-exon mature (exon→boundary). The densities `rho_g/rho_pos/rho_neg` are the Theorem-1 currency.
- **Confirmed correct.** `use_shift` on clean edges → the §4a FL-frame shift; the eff-length frame conversion
  (`Mg = rho_g·egd`, normalize by the imputed total) is §3/§4.
- **The arithmetic says the density mode is unnecessary.** The exon-facing `use_shift=False` density branch
  (`log(ρ_c / (md/E_c))`, [bp_solver.py:643](../../src/rigel/calibration/bp_solver.py)) is the observed-total
  anchor; §4b shows the *shift* (imputed-total) is exact on those edges once mature is reconciled. Removing it is
  behavioral (a **precision**, not a mode, question) — deferred, gated by an A/B per the implementation plan.

---

## 9. Deferred to the next turn

1. **Message precision** — `pr_c` for each mode: the deconvolution-prediction precision (τ / reference-free
   evidence, count sampling `1/M`, transfer variance `σ²_transfer`) and the spliced-measurement precision
   (`S_eff/(1+S_eff·σ²_transfer)`), and **how they merge** when a boundary carries both nascent-prediction and
   spliced-measurement RNA (the owner's open question — a precision-weighted combine, converging to the confident
   channel as the other's mass → 0).
2. **Solve gating** — the DOF "can I solve?" criterion and the pass-0 skip rule ([`message_system_derivation.md`](message_system_derivation.md) §6B).
3. **The unequal-gate lower-anchor ψ form** — half-quadratic vs soft-hinge vs one-sided Gaussian (§5;
   [`message_system_derivation.md`](message_system_derivation.md) open issue #2).
