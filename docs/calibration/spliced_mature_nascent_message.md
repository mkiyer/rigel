# Spliced → mature → nascent: the RNA imputation message, from first principles

**Status:** the **mature-absorption** term implemented + validated on the toy 2026-06-30
(`calib-disagreement-precision`). Supersedes the NEXT_SESSION task-A premise (the `eff_spl` "~2×"
bug — **disproven**, see §5). Companion: `dispersion_aware_message_precision.md` (the gDNA/nascent
*precision*); this doc is the RNA *mode* (mature vs nascent) geometry.

---

## 1. Three molecular species, three channel signatures

Calibration deconvolves each node's **unspliced** mass into a pie `(f₊, f₋, f_g)`. The RNA part is
two physically distinct species that the strand model cannot tell apart but the **splice** channel
can:

| species | extent | crosses an intron–exon boundary as | within an exon, as |
|---|---|---|---|
| **gDNA** ρ_g | whole genome (≈uniform) | UNSPLICED (contiguous) | contained unspliced |
| **nascent** ρ_n (pre-mRNA) | whole gene body | UNSPLICED (contiguous) | contained unspliced |
| **mature** ρ_m (spliced mRNA) | exons only (introns removed) | **SPLICED** (gapped) — *never* unspliced | contained unspliced |

The deposit fills the mate gap and subtracts the cut introns (`fragment_genomic_spans`,
`bam_scanner.cpp`), so a mature fragment spanning a junction deposits as a **spliced crossing**; a
contiguous (gDNA / nascent) fragment deposits as an **unspliced crossing**.

Consequences that drive everything below:
- A node's **unspliced crossing** mass (a boundary) is **gDNA + nascent only** — mature is absent
  (it skips the intron as spliced). So a boundary's pie RNA fraction is **pure nascent**.
- A node's **contained unspliced** mass (an exon) is **gDNA + nascent + mature** — mature fragments
  that fit inside the exon look unspliced. So an exon's pie RNA fraction is **nascent + mature**.
- The **spliced** mass lives only on a boundary's **exon flank** (one-sided) and directly measures
  the mature density.

## 2. The contiguous (imputable) field is NASCENT, not total RNA

gDNA and nascent are contiguous along the genome → they are *imputed* across edges (the smooth
field; precision per `dispersion_aware_message_precision.md`). **Mature is not contiguous across
introns** → it is **not** an imputation. Mature enters the model only at a junction boundary, two
ways:
- **into the adjacent exon** — the same mature molecule continues into the exon body (a MEASUREMENT);
- **out of the exon, absorbed by the spliced channel** — when the exon is imputed onto the boundary,
  the part of the exon's RNA that is mature crosses as spliced, not as the boundary's unspliced
  crossing, so it must be removed.

So the RNA imputation message communicates the **nascent** density. The error fixed here: the
region→boundary message was imputing the exon's **total** RNA (mature+nascent) as nascent → mature
leaked into introns (the wholesale nascent hallucination).

## 3. Effective-length arithmetic (the densities)

FL pmf `f(ℓ)`; `fl_mean = Σ ℓ f(ℓ)`. For a region of length `L` and a boundary side bounded by a
region of length `R` (`effective_length.py`):

- **contained** count eff-len `E_r(L) = E[max(0, L−ℓ)]` — recovers a contained density `ρ = M/E_r`.
- **two-sided crossing** eff-len `E_side(R) = E[min(ℓ, R)]` — for a contiguous (gDNA/nascent)
  crossing.
- **one-sided spliced** eff-len `E_spl(R) = E[min(ℓ,R)²/(2ℓ)] → fl_mean/2` for `R ≫ ℓ`.

**Derivation of `E_spl` (the half-triangle).** A mature fragment of length ℓ crossing a junction has
donor coverage `a ∈ {1,…,min(ℓ−1, R)}` on its exon flank and deposits mass `a/ℓ` there. Per unit
mature density ρ_m, the flank's spliced mass is
`ρ_m · Σ_ℓ f(ℓ) Σ_{a=1}^{min(ℓ−1,R)} a/ℓ = ρ_m · E[min(ℓ,R)²/(2ℓ)]`, hence

```
ρ_m  =  M_spliced(exon flank) / E_spl(R)
```

**Validated** (`scripts/debug/spliced_efflen_calib.py`, high depth, internal exons 300–5000 bp):
`E_spl_effective = M_spl / ρ_RNA_true` matches the coded `fl_mean/2` to 2–6% and is
**length-independent** (the residual is the few-percent read-anchor/mate-gap loss from the spliced
channel). The `eff_spl` geometry is correct — do not rescale it.

### The mature/nascent split at an exon
The exon's contained RNA density is `ρ_RNA = f_RNA·M / E_r(L) = ρ_m + ρ_n` (both species fill the
contained eff-len). The junction's spliced gives `ρ_m`. Therefore the exon's **nascent** is the
residual:

```
ρ_n(exon)  =  f_RNA·M_exon / E_r(L)  −  ρ_m,        ρ_m = M_spl(junction flank) / E_spl
```

This is "the spliced mass absorbs the incoming RNA; the remainder is nascent."

## 4. The unified message (no gate — it falls out of one-sided spliced)

For an RNA message `src → dst` (per strand), in **density** units, re-expressed in the dst frame:

```
ρ_msg = (src nascent)              fbp·M_src / E_r_src
      + (src-face mature  ADDED)   SP[sf][src] / E_spl_src     # B→exon MEASUREMENT
      − (dst-face mature  ABSORBED) SP[df][dst] / E_spl_dst    # exon→B  (leaves nascent)
```

Because spliced is **one-sided** (exon flank only), at most one of the two mature terms is non-zero
on any edge, and **both are zero on intron↔boundary edges** — so introns receive pure nascent with
**no gate**:

| edge | src mature (add) | dst mature (absorb) | message |
|---|---|---|---|
| junction-B → exon | `SP[sf]>0` | 0 | nascent + **mature MEASUREMENT** |
| exon → junction-B | 0 | `SP[df]>0` | RNA − mature = **nascent** |
| B → intron | 0 | 0 | nascent |
| intron → B | 0 | 0 | nascent |

**Precision** (unchanged classes): the MEASUREMENT (`src mature present`) uses count precision and is
not disagreement-silenced (a gapped read *is* mature RNA, fused with the strand likelihood so a
confident strand wins and a strand-blind exon snaps to it). Every other case is a nascent IMPUTATION
at the disagreement-aware precision. The mature absorption only corrects the imputation **mode**
(subtracts ρ_m), keeping its IMPUTATION class.

Implementation: `bp_solver.py::node_sweep._scan`, blocks `emit_p`/`emit_n`. `SPd/SNd/ESPd = SP[df]/
SN[df]/ESP[df]` (the dst-face spliced). One subtraction line per strand.

## 5. What this fixes (and what it does NOT)

**Fixes** the message-driven intron nascent hallucination. TA+ `(1000,2000),(5000,10000)`, nascent=0,
stranded no-capture: intron f_g 0.866 → **0.972** (truth 1.0); exons unchanged; no regressions on the
other 3 conditions or the 188 calibration unit tests. (`toy_prod.py`.)

**Does NOT fix** the *unstranded* / *capture* intron leak (TA+ S2 0.54, S4 0.70) — there the leak is
**local**, not message-driven: a flat strand likelihood (κ≈½) plus a global gDNA prior that is the
population **mean**, not the depleted **floor**, fails to pin a depleted intron to f_g≈1. That is the
empirical-prior / "global minimum-not-mean" arc (NEXT_SESSION task D). Removing the boundary spliced
floor (physically correct — the floor double-counts mature into the boundary's unspliced pie) helps
the no-capture introns but **unleashes the capture density-cliff** (a depleted boundary's gDNA
message drags an enriched exon down) because ê(z) is degenerate on these small toys — so it is
**deferred** to land with the empirical prior. See `spliced_efflen_not_2x_nascent_subtraction` memo.
