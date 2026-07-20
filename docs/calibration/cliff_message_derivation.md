# Cross-node message propagation across the density cliff — the effective-length-frame log-odds shift

**Status:** derivation (2026-07-20). Derives, from first principles and validates by Monte-Carlo across many
fragment-length distributions, the **correct arithmetic for a message between a region and a boundary** — the
only kind of message in the calibration chain. Resolves the density-cliff propagation bug: the current message
mode fails catastrophically across a capture enrichment/depletion cliff (isolated intron→boundary `|Δf_g| = 0.65`
under capture). The result is a single clean rule — a **log-odds shift** — that is exact across FL distributions
and **invariant to the capture cliff**. This note covers the clean case (**intron ↔ intron-exon boundary**) in
full and lays out the framework for extending to every message type (the mature-bearing exon transitions).

Companions: [`gdna_intron_factory_design.md`](gdna_intron_factory_design.md) (the intron factory that supplies
accurate source beliefs), [`background_reference_derivation.md`](background_reference_derivation.md) §3 (SCALE ×
SHAPE: enrichment cancels in composition), [`effective_length.py`](../../src/rigel/calibration/effective_length.py)
(the per-component eff-lengths), [`density_cliff_and_mature_absorption.md`](density_cliff_and_mature_absorption.md)
(the mature removal, the harder extension). Prototype + validation: `scripts/debug/cliff_message_mc.py`.

---

## 0. TL;DR

Every message is between a **region** (measures *contained* fragments) and a **boundary** (measures *crossing*
flux). Write `λ = logit(f_g)` — the gDNA-vs-RNA-total log-odds. The correct message is a **shift**:

```
    λ_dst  =  λ_src  +  [ log(E_g^dst / E_g^src)  −  log(E_r^dst / E_r^src) ] ,                       (0.1)
```

where `E_c^node` is the node's **per-component** effective length (contained for a region, crossing for a
boundary), from `effective_length.py`. **The capture enrichment cancels identically** (0.1 has no `e`), so it is
cliff-invariant; only the composition propagates. The shift is nonzero precisely because gDNA and RNA have
**different FL distributions**, so their eff-length *ratio* differs between the contained and crossing frames.

Validated by Monte-Carlo (`cliff_message_mc.py`) against ground-truth simulated fragments across gaussian /
gamma / bimodal / uniform FL pairs, and for enrichment `e ∈ {1, 30, 300}` (boundary `f_g` invariant to the
digit). The naive "`f_g` is invariant" (shift = 0) is wrong whenever the FL distributions differ (e.g. the true
boundary `f_g` is 0.57 while the intron reads 0.38). The current per-node density message — divide the source's
gDNA density by the destination's **single-component observed total** `M_dst/E_dst` — conflates the frame and the
cliff and is the bug.

---

## 1. The two node kinds and why `f_g` is frame-dependent

The chain alternates **region ↔ boundary ↔ region**. The two kinds measure *different geometric events*:

- A **region** of physical length `L` measures fragments **wholly contained** in it. A length-`ℓ` fragment is
  contained from `E_f[max(0, L − ℓ + 1)]` start positions — the count effective length
  `region_eff_length(FL_c, L)`. Different components use different FL pmfs (`gdna_fl_pmf` / `rna_fl_pmf`).
- A **boundary** is a genomic *point*; it measures **crossing flux** — fragments that span the point. A
  length-`ℓ` fragment crosses from `min(ℓ, flank)` start positions — the crossing effective length
  (`boundary_side_*_eff_length`, `→ fl_mean` for large flanks).

A node's `f_g` is a **count (mass) fraction**: `f_g = M_g / (M_g + M_r)` with `M_c = ρ_c · E_c^node`. Because
`E_g^node` uses the **gDNA** FL and `E_r^node` the **RNA** FL, and these FL distributions differ, **the count
`f_g` depends on which frame (contained vs crossing) the node is in** — even when the underlying composition is
identical. This is the crux the owner identified: *"boundaries have an effective length that depends on their
composition; if 100% gDNA, density = counts/(DNA FL); if 100% RNA, density = counts/(RNA FL)."*

---

## 2. The physical model and the invariant

At a genomic location `x`, let `d_c(x)` be the **intrinsic** (pre-capture) density of component `c` (molecules
per bp), and `e(x)` the **capture enrichment** (nucleic-acid-agnostic — it scales gDNA and RNA alike,
`background_reference_derivation.md` §3). The observed rate is `ρ_c(x) = e(x)·d_c(x)`.

The **composition** — the ratio `d_g : d_r` — is set by biology and is **locally uniform** for adjacent nodes
that share the same molecules (an intron and its flanking boundary see the *same* gDNA and the *same* nascent
pre-mRNA). The **density** is not uniform: capture makes it jump ~10²–10³× across the intron→exon cliff. So:

> **The invariant across the cliff is the composition `d_g/d_r` (enrichment cancels). The absolute density and
> the count `f_g` are frame-dependent. A message must propagate the composition, not the density or the count.**

For **synthetic data the composition is exactly uniform**, so a boundary *must* inherit its intron's composition
— any residual is sparsity or an arithmetic error, never biology (owner, 2026-07-20). Empirically the *density*
composition matches (intron 0.643 vs boundary 0.627 at capOFF); the *count* `f_g` differs (0.648 vs 0.546) purely
by the frame — exactly §1.

---

## 3. The derivation — a log-odds shift

A node's gDNA-vs-RNA-total log-odds:

```
    λ  =  logit(f_g)  =  log( M_g / M_r )  =  log( ρ_g · E_g^node / (ρ_r · E_r^node) )
       =  log(ρ_g / ρ_r)  +  log( E_g^node / E_r^node ) .                                            (3.1)
```

Split the density ratio into intrinsic × enrichment: `ρ_g/ρ_r = (e·d_g)/(e·d_r) = d_g/d_r` — **the enrichment
`e` cancels** (3.1 depends on `e` nowhere). So `log(ρ_g/ρ_r) = log(d_g/d_r)` is the invariant composition
log-odds, shared by source and destination. Subtracting (3.1) for the destination and source:

```
    λ_dst − λ_src  =  log(E_g^dst/E_r^dst) − log(E_g^src/E_r^src)
                   =  [ log(E_g^dst/E_g^src) − log(E_r^dst/E_r^src) ]  ≡  Δ(src→dst) .                (3.2)
```

**The message is a pure additive shift in `λ` by the per-component eff-length frame change (3.2).** The mode is
`λ_dst = λ_src + Δ`; enrichment never appears. Equivalent count-space form (what the code computes):
`f_g^dst = ρ_g^src·E_g^dst / (ρ_g^src·E_g^dst + ρ_r^src·E_r^dst)` with `ρ_c^src = f_c^src·M_src/E_c^src` — i.e.
impute the source's per-component **densities** onto the destination's **per-component** eff-lengths and
renormalize. (`½` differences between count- and density-eff-lengths cancel in `Δ`: they appear in both
`log E_g^dst` and `log E_r^dst`.)

**Properties.**
- `Δ = 0` when gDNA and RNA share the FL distribution (`E_g/E_r` identical in both frames) — then `f_g` *is*
  frame-invariant, the naive assumption.
- `Δ ≠ 0` when the FL distributions differ — the general, load-bearing case.
- `Δ(dst→src) = −Δ(src→dst)` — reversible, as a proper coordinate change must be.
- No `e`: a message crossing a 1000× cliff is byte-identical to one with no cliff. **This is the whole fix.**

---

## 4. Monte-Carlo validation across FL distributions (`scripts/debug/cliff_message_mc.py`)

Simulate gDNA (FL = gdna_fl) and nascent-RNA (FL = rna_fl) fragments over an intron + adjacent exon; deposit
*contained-in-intron* vs *crossing the intron|exon point*; measure the **ground-truth** `f_g` in each frame;
compare the boundary `f_g` to the shift prediction (3.2). Fragments overlapping the exon are enriched by `e_exon`
(positional capture).

| gDNA FL / RNA FL | fg_intron | **fg_boundary (MC truth)** | pred (shift 3.2) | naive (= intron) | Δ |
|---|---|---|---|---|---|
| gauss 200 / gauss 200 | 0.400 | 0.3996 | **0.4000** | 0.400 | 0.000 |
| gauss 300 / gauss 150 | 0.380 | 0.5720 | **0.5714** | 0.380 ✗ | +0.778 |
| gauss 120 / gauss 350 | 0.432 | 0.1857 | **0.1860** | 0.432 ✗ | −1.201 |
| gamma 300 / gamma 180 | 0.384 | 0.5271 | **0.5262** | 0.384 ✗ | +0.579 |
| bimodal / gauss 220 | 0.398 | 0.4136 | **0.4142** | 0.398 | +0.066 |
| unif 250±100 / unif 150±60 | 0.387 | 0.5269 | **0.5263** | 0.387 ✗ | +0.566 |

The shift prediction matches ground truth to MC noise for **every** distribution shape; the naive assumption is
off by up to 0.25 in `f_g`. **Enrichment invariance:** `e_exon = 1 / 30 / 300` → boundary `f_g = 0.5720` at every
level (unchanged to 4 digits). Both claims of §3 confirmed.

**Isolated in-tool check** (gdna300 sim, factory-solved intron source, exon side removed): the DENSITY message
`|Δf_g|` vs oracle is **0.651** under the cliff (catastrophic) vs **0.059** without; the shift/COMPOSITION message
is **0.167 / 0.051** — the cliff bug is the density mode's single-component denominator, and the shift fixes it.

---

## 5. Per-strand — the tilt `θ` is preserved

RNA carries a strand (`f_pos`, `f_neg`); gDNA does not. Both RNA strands are the **same molecule population**
sampled by the **same** RNA FL, so they share `E_r`. Hence the frame shift (3.2) — which is `log(E_g/E_r)` — acts
**only on `λ`** (gDNA vs RNA-total) and leaves the tilt `θ = arcsin(f_pos − f_neg over f_r)` **unchanged**. This
matches the solver's coordinate split: the message is a shift on the `λ` axis, `θ` relayed as-is. (Strand
*continuity* — whether a strand's nascent crosses the boundary at all — is the separate `free_s` structural gate,
unchanged.)

---

## 6. The full message taxonomy — every region ↔ boundary type

Every message pairs one region with one boundary, in one of two directions. What changes across types is **which
components are present** and **whether mature RNA must be reconciled**:

| region | components (contained) | boundary crossing (unspliced) | reconciliation |
|---|---|---|---|
| **intergenic** | gDNA only | gDNA | none — pure shift (trivially `f_g=1`) |
| **intron** | gDNA + nascent | gDNA + nascent | **none — pure shift (§7)** |
| **exon** | gDNA + nascent + **mature** | gDNA + nascent (mature **splices out**) | mature **removed** outbound / **added** inbound (§8) |

The key physical fact: **mature RNA lives *contained* in exons but does NOT cross a junction as *unspliced*** — a
mature molecule spanning the junction is a **spliced** read, counted separately. So a boundary's *unspliced*
crossing is always gDNA + nascent, never mature. Therefore:

- **intron/intergenic ↔ boundary:** the source and destination share the same component set (gDNA + nascent) →
  the **pure log-odds shift (3.2)** is the entire message. *(This note's fully-derived case, §7.)*
- **exon → boundary:** the exon's contained composition includes mature; the boundary's does not. The message
  must **remove the mature** from the source RNA before the shift (it does not cross). *(§8.)*
- **boundary → exon:** the boundary carries no mature; the exon does. The message must **add the exon's mature**
  to the destination RNA after the shift. *(§8.)*

The mature quantity is measured by the junction's **spliced** count `S`, projected to the contained/crossing
frame by the spliced eff-length `spliced_side_eff_length` — the geometry of
[`density_cliff_and_mature_absorption.md`](density_cliff_and_mature_absorption.md), now layered *on top of* the
validated shift rather than mixed into a saturating density mode.

---

## 7. The clean case, in full: intron ↔ intron-exon boundary

Both nodes carry exactly {gDNA, nascent-RNA}; no mature; the intron is off-target (the factory gives an accurate
source belief). The message is the pure shift (3.2), in both directions:

- **intron → boundary** (`E^src` = contained, `E^dst` = crossing):
  `Δ = log(E_g^cross/E_g^cont) − log(E_r^cross/E_r^cont)`, evaluated with the boundary's actual flank size
  (`boundary_side_eff_length(FL_c, R)`), `→ log(fl_mean_g/(L−fl_mean_g)) − log(fl_mean_r/(L−fl_mean_r))` for a
  long intron.
- **boundary → intron**: `Δ` with source/dst swapped = the negative of the above.

The precision is unchanged in spirit (§9): the source's own belief precision + the count term + the belief-free
`σ²_transfer` enrichment-crossing damping. The mode is the shift; nothing divides by a depleted quantity.

**This is exactly the composition message mode** (`bp_solver._COMPOSITION_MODE`, **LANDED 2026-07-20**,
commit `d8fef478`, per-edge `use_shift` gate on clean non-exon transitions): its `Mg = ρ_g^src·E_g^dst`,
`Mp = ρ_r^src·E_r^dst`, `f_g = Mg/(Mg+Mp+Mn)` is (3.2) rearranged. Zero-regression on the 32-scenario sim; the
aggregate boundary effect is small because the intron-exon boundaries are still dominated by the
higher-precision (density-mode) exon→boundary message — which §8 fixes.

---

## 8. Exon region ↔ intron-exon boundary — the mature reconciliation (DERIVED + MC-VALIDATED)

**Status:** derived + Monte-Carlo validated across FL distributions (`scripts/debug/cliff_exon_boundary_mc.py`);
implementation next. Grounded in §3's density-composition invariance + shift.

### 8.1 The component structure

An **exon region** observes ONLY *unspliced contained* fragments, per strand `s`:

```
    RNA_s^exon (contained)  =  nascent_s  +  mature_s          (both unspliced, in the exon body)
    gDNA^exon (contained)
```

An **intron-exon boundary** observes:

```
    unspliced crossing:  gDNA,  nascent_+ ,  nascent_−     (nascent crosses unspliced)
    SPLICED:             mature on the junction's ONE motif strand   (mature crosses ONLY as spliced)
```

The physical crux: **mature RNA has the SAME density `ρ_m` measured in two frames** — *contained-unspliced* in
the exon body (eff-len `E_r^cont = region_eff_length(RNA_FL, L)`) and *spliced-crossing* at the junction (eff-len
`E_spl`, the one-sided half-triangle `spliced_side_eff_length`). MC-verified: `ρ_m` recovered as
`mature_contained/E_r^cont` and as `spliced/E_spl` agree to ~0.3% across FL shapes. This is the **link** that
lets a boundary's spliced measurement reconcile an exon's contained mature.

The boundary's *unspliced* crossing carries **no mature** (it splices), so the exon's `f_g` (which counts mature
as RNA) and the boundary's unspliced `f_g` differ by the mature — reconciled below. A single junction carries
spliced on **one strand** only (the splice motif is stranded), so the reconciliation is applied to that strand.

### 8.2 Density arithmetic (the invariant), both directions

Work in the invariant **densities** (`ρ_g`, `ρ_n±` nascent, `ρ_m±` mature), then convert with the destination's
per-component eff-lengths (the shift, §3). Enrichment cancels throughout (MC: `e = 1` vs `100` → `f_g` identical).

**Boundary → Exon (ADD mature).** The boundary supplies `ρ_g`, `ρ_n±` (from its unspliced crossing:
`ρ_n = n_cross/E_r^cross`) and `ρ_m` on the junction strand (from its spliced: `ρ_m = S/E_spl`). The exon's RNA
density on that strand is `ρ_RNA_s = ρ_n_s + ρ_m_s`. Apply the exon's contained eff-lengths:

```
    f_g^exon  =  ρ_g·E_g^cont  /  ( ρ_g·E_g^cont  +  Σ_s (ρ_n_s + ρ_m_s)·E_r^cont ) .                  (8.2)
```

**Exon → Boundary (SUBTRACT mature).** The exon supplies `ρ_g`, `ρ_RNA_s = ρ_n_s + ρ_m_s` (it cannot split
nascent from mature). The boundary removes the mature *it directly measures* on its strand, `ρ_m_s = S/E_spl`, to
recover the nascent that actually crosses unspliced, then applies its crossing eff-lengths:

```
    ρ_n_s  =  max( ρ_RNA_s^exon  −  ρ_m_s^bnd ,  0 ) ,
    f_g^bnd(unspliced)  =  ρ_g·E_g^cross  /  ( ρ_g·E_g^cross  +  Σ_s ρ_n_s·E_r^cross ) .                (8.3)
```

MC (`cliff_exon_boundary_mc.py`, `ρ_g=0.4, ρ_n=0.3, ρ_m=0.6`): the reconciliation `ρ_RNA^exon − ρ_m = ρ_n`
holds to ~1%, and both (8.2)/(8.3) reproduce the ground-truth `f_g` to MC noise for gaussian / gamma / uniform
FL pairs, at every gDNA/RNA FL ratio.

### 8.3 The spliced-absorption GATE — WITH vs WITHOUT (owner's insight)

The `± ρ_m` reconciliation fires **only when the boundary carries spliced mass on the exon's strand and side** —
i.e. the boundary is a genuine splice site of the exon's transcript. This is a signature/geometry flag, already
computed: the boundary's `junction_strand` + the exon-flank routing (`node_geometry._spliced_faces` routes
spliced mass to the same-strand exon flank). Two message types result:

- **WITH spliced absorption** (`S_s > 0` on the exon's strand/side): apply `± ρ_m_s` (8.2 add / 8.3 subtract) —
  the boundary is a **source** of that strand's mature (into the exon) / a **sink** for it (out of the exon).
- **WITHOUT spliced absorption** (`S_s = 0`): the exon's mature does not splice at THIS boundary (it is a TSS/TES
  or the exon's other end); the message is the **pure shift (§7)** — no mature term.

This is exactly the `±SPs / −absorb` structure already in `_scan`, but **gated** by the same-strand-spliced flag
(the missing piece), so the `ρ_m` used is the mature that actually reconciles at THIS junction — never a
mismatched subtraction. On a clean (intron/intergenic) edge `S = 0` identically ⇒ WITHOUT ⇒ §7 (self-consistent).

### 8.4 Both-sides-exon boundaries

Some boundaries flank an exon on BOTH sides (adjacent exons / an opposite-strand exon↔exon seam). These decompose
by the **same gate applied per side/strand**: each side is an exon↔boundary transition that is WITH spliced
absorption iff that side's exon strand carries junction spliced, else WITHOUT. No new arithmetic — the
per-side/per-strand spliced flag fully determines the reconciliation. (An exon↔exon seam with no splicing is
WITHOUT on both sides = pure shift; a mutually-spliced pair is WITH on the spliced strand.)

### 8.5 Precision (unchanged) and the earlier regression, resolved

Precision is the §9 machinery (source belief + `1/M` + `σ²_transfer`), plus the spliced count credits its own
measurement precision. The earlier composition-mode regression ([[composition_mode_regresses_post_tau]], the
"300× under-removal") is now understood: it applied `−absorb` to exon edges **ungated** and normalized by the
imputed total without the density-frame link — subtracting a mismatched mature. §8.3's `ρ_m = S/E_spl` is the
mature that reconciles *at this junction*, gated by §8.3's flag; where the exon's mature splices at its *other*
junction, THIS boundary's `S ≈ 0` ⇒ WITHOUT ⇒ no (wrong) subtraction. The MC confirms the gated arithmetic is
exact.

### 8.6 Implementation plan

1. **Gate** the `±SPs/−absorb` terms in `_scan` by the same-strand-spliced flag (§8.3) → extend `use_shift`'s
   scope from clean-only to *all* region↔boundary edges, with the reconciliation applied WITH/WITHOUT per the
   flag. On a clean edge this is identical to §7 (S=0).
2. **Densities from the correct eff-lengths** — `ρ_m = S/E_spl` (half-triangle) added/subtracted in *density*
   space before the mass conversion (8.2/8.3), not the current count-space `S·E_r/E_spl` mixed into `rho_pos`.
3. **Validate** — the exon↔boundary MC as the unit ground truth (already passing), then pass-0 vs oracle on all
   32 scenarios (boundary + exon `|Δf_g|` must improve — this removes the exon→boundary domination that masks
   §7), no regressions, stranded controls flat.
4. **Then** the gDNA hyperprior study on the now-healthy pass-0.

---

## 9. Precision (unchanged in spirit)

The shift changes only the message **mode**. The **precision** is the source's honest belief precision
`1/(Var(log f_c^src) + 1/M_src)` + the belief-free `σ²_transfer = var_proj[dst] + (μ_proj[dst]−μ_proj[src])²`
(the enrichment-crossing damping from the NPMLE projection — which is *where* the cliff enters the precision,
correctly, as reduced confidence, never as a biased mode). A boundary's crossing count is genuinely lower than an
intron's contained count, so its message is honestly less precise — the count term carries that.

---

## 10. Implementation plan

1. **Land the clean shift** (§7): make the composition mode (`_COMPOSITION_MODE`) the production message *mode*
   for intron/intergenic ↔ boundary transitions (no mature) — retiring the density mode's single-component
   denominator, which is the cliff bug. Gate + A/B as usual.
2. **Validate** on the cached ambig scenarios (pass-0 vs oracle, all 32), with the FL-general MC
   (`cliff_message_mc.py`) as the unit-level ground truth. Gate: intron↔boundary propagation improves across the
   cliff (isolated `|Δf_g|` 0.65 → ~0.17), no node-type regressions, stranded controls flat.
3. **Extend** to exon ↔ boundary (§8) — the mature reconciliation layered on the shift.
4. **Then** the gDNA hyperprior study, on the now-healthy pass-0 (owner sequencing: the hyperprior learns garbage
   from a broken pass-0).

---

## 9. The hybrid enrichment-corrected density mode — DERIVED, TESTED, and REJECTED as a unifier

**Status:** derived + MC + 32-scenario sim A/B (2026-07-20). The owner's directive was to unify the two landed
message modes (the composition **shift** on clean/intron edges, the **density + mature** mode on exon edges) into
one "hybrid enrichment-corrected density" mode that keeps the observed-total anchor `md` while conserving the
density composition across the cliff. This section records the derivation, the implementation (`bp._HYBRID_MODE`,
default OFF), and the **decisive finding that it regresses** — so it is not re-attempted.

### 9.1 The two candidate hybrid modes

A message provides two soft targets to the dst's λ-fold: a gDNA factor on `log f_g`, an RNA-total factor on
`log f_r`, at explicitly-set precisions (identical across all modes — only the targets change). With
`M_c^pred = ρ_c^src·E_c^dst`, `_den = ΣM_c^pred`, `md` the observed dst total, and
`dmu = logR = mu_proj[dst] − mu_proj[src]` (the belief-free denoised enrichment jump the `σ²_transfer` projection
already computes), the density and composition modes plus the **two** ways to build a hybrid enrichment-corrected
density:

```
    density         t_c = log(M_c^pred / md)               observed anchor; Σexp(t)=_den/md ≠ 1 at a cliff
    composition     t_c = log(M_c^pred / _den)             Σexp(t)=1; the cliff-invariant shift (§7)
    hybrid-ALL      t_c = log(M_c^pred / md) + dmu   (all c)      "keep md + correct the cliff"
    hybrid-gDNA     t_g = log(M_g^pred / md) + dmu ; t_{r±} = log(M_{r±}^pred / md)   (gDNA ONLY — owner option B)
```

`hybrid-gDNA` is the owner's `composition_mode_regresses_post_tau` option B: keep density's *decoupled*, error-prone-
mature-robust RNA on the pure density mode, correct **only** the gDNA numerator's cliff bias (`ρ_g^src·e_dst/e_src`
then ÷md). Both were implemented behind `bp._HYBRID_MODE` (the gDNA-only form is the one that shipped in the
toggle) and A/B'd. **Both regress.**

### 9.2 hybrid-ALL is not distinct from composition (the algebra; adversarially confirmed)

`t_c^{hyb-ALL} − t_c^comp = log(_den/md) + dmu ≡ C` — **the same constant for every component** `c` (an
independent adversarial derivation confirmed `C_g − C_r < 1e-12`). So hybrid-ALL = composition + a *uniform*
offset. Consequences: (1) the implied log-odds `t_g − t_r = log(M_g/M_r)` is **identical** to the composition
shift (the `+C` cancels), so it conserves the same composition; (2) `C = log(R·_den/md)`, and `_den = md/R`
**exactly when the source composition is accurate** ⇒ `C = 0` ⇒ hybrid-ALL = composition. The `+dmu` is precisely
what turns the density mode *into* the composition mode when R is right.

The observed `md` survives only in the residual `C ≠ 0` under source error — and it is **not** a "weak offset
always" (the adversary showed |C| ≳ 1 shifts folded f_g by 0.08–0.31, asymmetrically); it is that the *operative*
`C` on these edges is **small** (|C| ≤ ~0.3 at true R) because `_den` barely tracks the composition split. MC
(`scripts/debug/cliff_hybrid_mc.py`): with true R, hybrid-ALL is byte-identical to composition for an accurate
source (|Δf_g| 0.0045) and within 0.002 under injected error; with the *estimated* R (the total-density ratio,
which the frame shift + composition change contaminate) it is **worse**. No regime beats composition.

### 9.3 hybrid-gDNA IS distinct — and regresses hardest (the `mu_proj` contamination)

Adding `dmu` to the gDNA factor **only** does *not* cancel in `t_g − t_r`, so hybrid-gDNA is a genuinely different
mode: it keeps density's decoupled RNA + observed-`md` anchor and shifts only the gDNA log-fraction by the measured
enrichment jump. In principle it fixes density's `−logR` gDNA cliff bias while staying robust to the error-prone
mature. **In practice it is the worst arm on unstranded exons**, because `dmu = mu_proj[dst] − mu_proj[src]` is the
*total-density* jump, which at a composition-changing exon↔boundary edge **conflates enrichment with the mature/RNA
composition change** (RNA present in the exon, absent from the boundary's unspliced crossing; a pure-gDNA seam).
The gDNA correction then mis-fires exactly where density's decoupling was winning (the pure-gDNA seams), driving
confident-wrong mass up to **9.4 %** (gdna300_none_off CW 1.66 → 9.38). A clean gDNA-only enrichment ratio would
need the *gDNA* density jump, which is circular (it is what we are solving for).

### 9.4 The sim verdict — every hybrid regresses; baseline is best

4-way pass-0-vs-oracle A/B over all 32 cached `ambig_dense_10mb` scenarios (a session A/B harness that toggled the
message mode per arm — those mode toggles have since been collapsed to the single production path; the reusable
single-arm error diagnostic is `scripts/debug/pass0_error_diagnostic.py`), factory + τ ON in every arm,
mass-weighted |Δf_g| vs oracle:

| arm | exon (all 32) | exon (unstranded ss0.50) | boundary | CW% |
|---|---|---|---|---|
| **baseline** (shift on clean, density+mature on exon) | **0.249** | **0.392** | **0.274** | **0.6** |
| hybrid-gDNA (owner option B) | 0.272 | 0.428 | 0.281 | 1.3 |
| hybrid-ALL (density+dmu on all c) | 0.264 | 0.414 | 0.276 | 0.9 |
| composition (`_SHIFT_ALL_EDGES`) | 0.265 | 0.416 | 0.286 | 1.0 |

All hybrids regress (19–24 flagged conditions), all in unstranded `ss_0.50` (e.g. hybrid-gDNA gdna300_none_off
exon 0.270→0.460, CW 1.66→9.38). Stranded (`ss_0.99`) exons are **flat at 0.011 under every mode**. No hybrid
gives any aggregate benefit.

### 9.5 Why — unstranded IDENTIFIABILITY + `mu_proj` contamination, not mode arithmetic

1. **Identifiability.** Stranded exons solve to 0.005–0.025 under *every* mode; unstranded sit at 0.39. The entire
   gap is the **loss of strand information** — the only intrinsic gDNA/RNA signal (count-zero-info). No message
   *mode* can recover information the unstranded data does not contain.
2. **Composition-conservation transports source noise.** On an unstranded exon the source composition is noise;
   the shift/hybrid-ALL propagate it as **confident-wrong mass**. The density mode's target-inconsistency
   (`Σexp(t) ≠ 1`) keeps messages **fuzzy** → self-limits error propagation → defers to the observed anchor /
   population prior, which *is* the better estimate when the source is noise. So the baseline is right to use
   composition only where the source is trustworthy (**factory-anchored introns**) and density where it is not
   (**exons**) — **the current architecture already IS the enrichment-corrected hybrid**, gated by source
   reliability, keeping the observed `md` anchor exactly on the unreliable-source edges.
3. **The enrichment estimate is unusable at exon edges.** `mu_proj` is a total-density landscape; its difference
   is a clean enrichment ratio only on **composition-preserving** (clean intron) edges — which is why the shift
   works there — and is contaminated on composition-changing exon edges, breaking hybrid-gDNA.
4. **The trust knob is PRECISION, not the mode (adversary Claim 3, confirmed).** With the cliff-invariant log-odds
   fixed, a mode's only remaining freedom is one scalar (the common level `C`), which encodes total-mass/enrichment
   consistency — *orthogonal* to composition accuracy (silent on real source error, harmful under R bias). No
   fixed-precision mode can be a cliff-invariant, non-collapsing, genuine composition-*trust* regularizer; that
   role belongs to the message precision (`pr = 1/(Var(log f_c^src) + 1/M_src + σ²_transfer)`).

### 9.6 Where the real headroom is (NOT the mode)

To push unstranded exons below the 0.39 floor the lever is **not** the message mode. It is either (a) the message
**precision** — down-weight untrustworthy (composition-vacuous) sources so a composition-conserving mode becomes
safe (this is what the τ core reaches for; the residual regression means it is not yet fully gating the compounded
confidence along the chain), or (b) inject the one signal unstranded data *does* carry — the **motif-stranded
spliced mature** measurement (the `_RNA_ABSORB` clean win, already in the exon density path) and the
**population / intron-factory gDNA prior**. Both are precision/prior levers, out of scope for a mode change.

**Decision:** keep the baseline. `bp._HYBRID_MODE` and `_SHIFT_ALL_EDGES` stay default OFF (byte-identical;
retained as the A/B artifacts alongside this write-up, to be removed on owner sign-off). Related:
[`background_reference_derivation.md`], the `composition_mode_regresses_post_tau` and
`pass0_mature_subtraction_gap` memories.

---

## Appendix — the one line

`λ_dst = λ_src + log(E_g^dst/E_g^src) − log(E_r^dst/E_r^src)`. *Because composition is what's invariant across the
capture cliff, and the only thing that changes between a contained frame and a crossing frame is the
per-component effective length.* Everything else — enrichment, absolute density, the count `f_g` — is frame
detail the shift absorbs exactly.
