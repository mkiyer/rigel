# The pass-0 gDNA over-call: the REFRAME, and why it breaks at transcript ends

**Owner-directed dissection, 2026-07-27.** Written because `pin_derivation.md` §12.2 sized a defect it did
not diagnose: the gDNA density delivered into exons is **5.4× too large** at capture-OFF, and correcting the
frame is worth **20×** the P-2 residual it was found chasing.

Evidence: `scratchpad/p2r_{f,g,h,i}_*.py`. Figure: `figures/gdna_reframe_terminus.png`.
Condition throughout: `gdna100_ss_0.50_nrna_none_capture_off` unless noted (the worst of the stratum;
`gdna300_…` and the stranded twin behave identically, which is itself informative — **this is not a
strandedness effect**).

---

## 1. The decomposition is an identity, and it names the culprit

A message's delivered gDNA density is `tg = rg_src · r`, so against the truth at the destination

```
    log10( tg / ρ_g^true(dst) )  =  log10( rg_src / ρ_g^true(src) )          ← the SOURCE's own error
                                 +  log10( r )                               ← the REFRAME
                                 +  log10( ρ_g^true(src) / ρ_g^true(dst) )   ← TRUE spatial difference
```

Every term is separately observable. Measured over all 1,289 exon destinations with a live gDNA message
(mass-weighted; identity residual **4.4e-16**):

| term | mass-wt mean | median |
|---|---|---|
| **TOTAL** `log10(delivered/true)` | **+1.564 dec (37×)** | +0.134 |
| (1) the source's own gDNA error | **+0.110** | −0.018 |
| (2) **the REFRAME `log10 r`** | **+1.508** | +0.171 |
| (3) true spatial difference | **−0.054** | −0.008 |

**The sources are nearly right, gDNA really is uniform (term 3 ≈ 0, as the physics says), and the reframe is
96 % of the error.** There is nothing to fix in the imputation itself.

## 2. WHERE it breaks: the boundary carries a transcript TERMINUS

`r = ρ_tot(dst)/ρ_tot(src)` is a ratio of **total** densities and exists to cancel the hybrid-capture step.
Its faces come from `_rho_faces`, which adds the one-sided **spliced** density at the acceptor. So:

* at a **splice junction**, RNA crosses the seam, both faces carry it, and `r ≈ 1` — the gDNA transfers
  unchanged, which is correct;
* at a **transcript start or end**, *no RNA crosses*. The boundary face is **pure gDNA** while the exon it
  feeds is RNA-dominated, so `r` becomes the exon's entire RNA-to-gDNA ratio.

Classifying every source boundary by that single structural bit (does any transcript start or end there):

| source boundary | n | median `log10 r` | **median delivered gDNA error** |
|---|---|---|---|
| **TERMINUS** | 426 | **+0.836** | **+0.847 dec = 7.0× too big** |
| **junction-only** | 732 | **+0.021** | **+0.004 dec = 1.0× — EXACT** |
| neither | 131 | +0.314 | +0.407 (2.6×) |

Reproduced on `gdna300` (0.560 vs 0.062) and on the **stranded** twin (0.796 vs 0.017).
**33 % of edges are terminus edges and they carry 66–68 % of the error mass** at capture-OFF.

## 3. ⭐⭐ THE CLOSED FORM — the message degenerates into "you are 100 % gDNA"

If the source face is pure gDNA then `ρ_tot(src) = rg_src` identically, and the reframe collapses:

```
    tg  =  rg_src · r  =  rg_src · ρ_tot(dst)/ρ_tot(src)  =  ρ_tot(dst)
```

**The delivered gDNA density becomes the destination's own total density — the message carries ZERO source
information.** Verified, not derived-and-hoped: on **63–66 % of terminus edges** `tg == ρ_tot(dst)` to
within **1e-9**, on every condition tested including capture-ON. The composition it implies is

```
    f_g^claimed  =  tg·E_g/M  =  1.18 – 1.27
```

— the message tells an exon that is **99.5 % RNA** that it is **more than 100 % gDNA**.

> ⚠ **This is the same signature as the pin bug** (`pin_derivation.md`): a message whose content is a pure
> function of the destination's own mass. The pin manufactured it from the destination's own *belief*; this
> manufactures it from the destination's own *total density*. That is why P-2 exposed it — `_pin_v` was
> rescaling the claim back down by ~0.648 and hiding it.

## 4. Two worked nodes, end to end

**Node 2651** — `chr_syn:5,922,363-5,924,674`, region 1325, signature `0x1 [exon−]`, 2,311 bp.
Transcripts `G0276.1` / `G0276.2` (gene `G0276_anchor`, − strand, span 5,922,363-5,950,646); the region is
their **first exon**, and `5,922,363` is the transcript's outer end.

```
   RAW      u_pos=37,061  u_neg=37,371   M=74,432   E_g=2,212.0  E_r=2,112.1
   ORACLE   gDNA=401  RNA=74,031   true f_g=0.0054   true ρ_g=0.1813   total ρ=33.65
   OWN      og=16.49 op=0.00 on=17.62  f_g_own=0.483   τ_own=0  (NO intrinsic strand evidence)

   LEFT  ← boundary 2650 @ 5,922,363  = the transcript's OUTER END
         source truth: gDNA=15, RNA=0, f_g=1.000, ρ_g=0.1500
         relayed rg_src = 0.149961        (the source's own error: +0.000 — EXACT)
         r = ρ_face(dst)/ρ_face(src) = 34.4604 / 0.14996 = 229.8
         delivered tg = 34.4604   ← identical to ρ_tot(dst); truth is 0.1813  ⇒ 190× too big
         ORIGIN: source +0.000  +  reframe +2.361  +  spatial −0.082  =  +2.279 dec

   RIGHT ← boundary 2652 @ 5,924,674  = a splice junction
         r = 34.4604/35.2748 = 0.977      delivered tg = 0.1745 vs truth 0.1813  ⇒ 1.0×  ✅

   ψ RECEIVES  gDNA level f_g=0.427 (prec 12.8)  vs oracle 0.0054   ⇒ SOLVED f_g = 0.084
```

**Node 3087** — `chr_syn:6,943,828-6,946,101`, region 1543, `0x2 [exon+]`. `G0316.1` (+ strand) has exons
`(6,934,891-6,934,927)` and `(6,943,828-6,946,101)`: the region is its **last exon**, and `6,946,101` is the
transcript's end.

```
   RAW      M=51,118   E_g=2,174.0  E_r=2,074.1
   ORACLE   gDNA=410  RNA=50,708   true f_g=0.0080   true ρ_g=0.1886
   LEFT  ← boundary 3086 @ 6,943,828  = a splice junction (RNA crosses)   r = 6.03   → 5.0× too big
   RIGHT ← boundary 3088 @ 6,946,101  = the transcript END (no RNA)       r = 114.6  → 127.6× too big
         relayed rg_src = 0.209945 (source's own error −0.000, EXACT) ; delivered tg = 24.068 = ρ_tot(dst)
   ψ RECEIVES  gDNA level f_g=0.469 (prec 12.4)  vs oracle 0.0080  ⇒ SOLVED f_g = 0.357   (|err| 0.349)
```

**In both nodes the *same* boundary type behaves correctly and the terminus does not**, inside one node,
with everything else held fixed. That is as clean as an attribution gets.

## 5. What this is, structurally

**It is P1g.** The region/boundary map has no TSS/TES, so the solver cannot tell a splice junction from a
transcript terminus — the debt already on record for `ω_graft`, which splits **≥30×** on exactly this bit
(`ROADMAP`, `P1D_P1E_DEBTS.md`). This is the **second independent place** the same missing bit does damage,
and here it is not a mis-priced variance but a **mode** error of up to 190×.

⚠ **It also re-opens `ROADMAP` §11** ("REFUTED, do not build: per-channel enrichment ratios at the boundary
face"). That refutation rested on substituting the oracle capture step for `r` buying ≈ 0 — measured while
`_pin_v` was cancelling `r`. It does not bind on the P-2 tree.

## 6. Sizing, and why the obvious fix cannot ship

`r_g ≡ 1` on the plain reframe (prior-free — "gDNA density is uniform"):

| stratum | refit=0 | refit=1 |
|---|---|---|
| **unstranded × capOFF × gDNA** | **0.0495 → 0.0146 (−0.0349), 6/6 better** | **0.0313 → 0.0077 (−0.0237), 6/6** |
| capture OFF (all 14) | 0.0339 → **0.0161**, 9 b / 0 w | 0.0204 → **0.0080**, 9 b / 0 w |
| stranded × capOFF | −0.0028 | −0.0025 |
| `gdna_none` | −0.0031 | −0.0022 |
| ⛔ unstranded × capON | **+0.1728** | **+0.2069** |

At capture-ON the reframe carries the real enrichment step and dropping it is fatal (the un-reframed error
there is 1.608 dec against the reframed 0.639). **The error flips sign exactly at capture**, so this arm
bounds the prize; it is not a candidate.

## 7. The candidates, and what has to be true of any of them

Any fix must (i) leave `r` intact where it carries a capture step, (ii) use no destination *belief*, and
(iii) add no tuned constant.

1. ⭐ **Put TSS/TES in the region/boundary map (P1g), then treat the terminus face structurally.** At a
   terminus the correct statement is *"no RNA crosses this seam"*, which is exactly what the face is
   currently getting wrong. This is the fix the ROADMAP already has queued for `ω_graft`, and it would pay
   twice. **Recommended.**
2. **A gDNA-specific frame.** `r_g = ρ_g^init(dst)/ρ_g^init(src)` instead of the total-density ratio. At
   `_RHO_ITERS = 1` the frames come from `_init_belief()`, which is belief-free (P4b), so this is at no
   worse BP standing than today's `r`. ⚠ Untested, and the init composition is the uninformative ½ on
   unstranded data, so it may just rescale the same error.
3. ⛔ **Not a variance.** The delivered mode is wrong by up to 190×; damping it cannot move it toward truth.
   This is the P1e lesson (`variance_ledger.md` §6) and it applies verbatim.

**Nothing is implemented.** The tree is unchanged: `src/` remains bit-identical 32/32 to the P-2 baseline.
