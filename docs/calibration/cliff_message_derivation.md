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

**This is exactly the composition message mode** (`bp_solver._COMPOSITION_MODE`, currently dormant): its
`Mg = ρ_g^src·E_g^dst`, `Mp = ρ_r^src·E_r^dst`, `f_g = Mg/(Mg+Mp+Mn)` is (3.2) rearranged. Landing the shift =
enabling that mode for the clean transitions, with the density mode's single-eff-len denominator retired.

---

## 8. The extension framework — mature-bearing (exon) transitions

Layer the mature reconciliation onto the shift (do **not** re-mix it into a density denominator — that was the
saturating bug). For a source exon → boundary, in count space:

```
    nascent_c^src  =  max( f_c^src · M_src  −  mat_c ,  0 ) ,    mat_c  =  S_c · (E_r^src / E_spl^bnd)          (8.1)
```

(the contained mature `mat_c`, from the junction spliced `S_c` scaled to the exon body by the spliced
eff-length; `density_cliff_and_mature_absorption.md` §4.2), then apply the shift (3.2) to `{f_g^src, nascent}`.
For boundary → exon, **add** the destination exon's mature (`+S_c·E_r^dst/E_spl^bnd`) after the shift. The known
open item is the estimator for `mat_c` at exons whose mature splices at the *other* junction
([[composition_mode_regresses_post_tau]]); the shift itself is exact and FL-general (§4). **Sequencing:** land the
clean shift (§7), validate no-regression, then the mature reconciliation.

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

## Appendix — the one line

`λ_dst = λ_src + log(E_g^dst/E_g^src) − log(E_r^dst/E_r^src)`. *Because composition is what's invariant across the
capture cliff, and the only thing that changes between a contained frame and a crossing frame is the
per-component effective length.* Everything else — enrichment, absolute density, the count `f_g` — is frame
detail the shift absorbs exactly.
