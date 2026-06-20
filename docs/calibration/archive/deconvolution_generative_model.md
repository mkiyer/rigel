# The gDNA / nascent / mature deconvolution — the complete generative model

**Status:** first-principles derivation, rev 2 (2026-06-11). Rev 1 wrongly modelled the unspliced pool as
"gDNA + nascent." Corrected here: **mature RNA also contributes unspliced fragments** (those contained
within a single exon, not crossing a junction). This is the complete three-species model and what is
identifiable from it. A foundation to react to, not yet an implementation plan.

## 1. Three species, and how each makes fragments

| species | spliced fragments? | unspliced fragments? | strand | where | fragment length |
|---|---|---|---|---|---|
| **gDNA** | no | **yes** (always — DNA is contiguous) | unstranded (½) | exon, intron, intergenic | genomic FL |
| **nascent** (pre-mRNA) | no | **yes** (introns retained ⇒ contiguous) | sense rate κ | exon + intron (gene body) | RNA FL, genomic span |
| **mature** (spliced mRNA) | **yes** (junction-spanning) | **yes** (within a single exon) | sense rate κ | exon only | RNA FL, *exon-limited* |

The correction is the **mature unspliced** cell. A mature molecule fragmented for sequencing yields a
*spliced* read when the fragment straddles an exon–exon junction (the intron shows as the alignment gap)
and an *unspliced* read when the fragment falls entirely inside one exon. So:

```
unspliced in an INTRON      = gDNA + nascent
unspliced in an EXON        = gDNA + nascent + mature-contained
unspliced in INTERGENIC     = gDNA
spliced (anywhere)          = mature only          (junction-spanning)
```

## 2. What you observe at a node, and the Poisson means

Per region/node, four counts (unspliced oriented to transcript sense; spliced is motif-oriented):
`U_sense, U_anti, S_sense, S_anti`. Three latent rates — `ρ` (gDNA density), `ν` (nascent rate), `m`
(mature rate) — and four **known-geometry** effective lengths: `E_g, E_ν` (genomic span; gDNA-FL /
RNA-FL), `E_mu` (mature **contained within the exon** — RNA-FL, *exon-limited*), `E_ms` (mature
**junction-spanning** — RNA-FL). The thinned-Poisson means:

```
E[U_sense] = ½·ρE_g + κ·νE_ν + κ·mE_mu          # gDNA splits ½/½ ; all RNA reads sense at κ
E[U_anti ] = ½·ρE_g + (1−κ)·νE_ν + (1−κ)·mE_mu
E[S_sense] = κ·mE_ms                            # spliced = mature only
E[S_anti ] = (1−κ)·mE_ms
```

## 3. Identification — four observables, three latents ⇒ all three are recoverable

With κ known, the system is **over-determined** (4 equations, 3 unknowns), and the solution is a clean
cascade:

1. **Mature** from the spliced reads: `m̂ = (S_sense + S_anti) / E_ms`.
2. **Mature-contained unspliced is then *predicted*** — `mE_mu = m̂·E_mu` — and **subtracted** from the
   unspliced pool: `U'_sense = U_sense − κ·m̂E_mu`, `U'_anti = U_anti − (1−κ)·m̂E_mu`. *This is the new,
   essential step the rev-1 model was missing — the spliced reads tell you how much mature is hiding in
   the unspliced, via the geometric ratio `E_mu/E_ms`.*
3. **Nascent** from the strand imbalance of the residual: `νE_ν = (U'_sense − U'_anti)/(2κ−1)`.
4. **gDNA** from the rest: `ρE_g = (U'_sense + U'_anti) − νE_ν`.

So all three are identifiable — mature from the junctions, nascent from the strand tilt of what's left,
gDNA from the unstranded remainder.

## 4. The gDNA prior is *robust* to the RNA composition — the key well-posedness result

The calibration's job is the per-locus **gDNA-vs-RNA** split, not the nascent/mature split. And for that:
**gDNA is the only unstranded species; everything else (nascent *and* mature) reads sense at κ.** So the
strand isolates gDNA without needing to know the nascent/mature mix:

```
gDNA = (U_sense + U_anti) − (U_sense − U_anti)/(2κ−1)   →  2·U_anti  as κ→1
```

— and crucially the mature-contained term **cancels out of this** exactly the same way nascent does (both
are κ-stranded). So discovering that the unspliced pool is a *three*-way mix does **not** change the gDNA
estimate: the strand answer `gDNA = 2·antisense` (high κ) holds with all three present. That is why the
prior is well-posed even though the unspliced is a three-way mixture — the question we ask of it
(unstranded vs stranded) is two-way. The mature subtraction (§3 step 2) matters for *splitting RNA into
nascent vs mature*, which is the EM's job downstream, **not** for the gDNA prior.

## 5. So why does the completeness matter? Three concrete reasons

1. **Effective lengths differ by species (the FL-consistency axis).** A gDNA or nascent fragment over an
   exon may spill into flanking introns (genomic span); a *mature-contained* fragment cannot (it would
   cross the junction → spliced). So `E_mu < E_ν ≈ E_g` for short exons. When converting the gDNA *count*
   to a *density* or *fraction*, the RNA denominator's effective length is a **mixture** (genomic nascent
   + exon-limited mature), not one number. Getting this wrong is the short-exon bias the FL-consistency
   work already chases — now correctly attributed to the species mix.
2. **The spliced reads are an active constraint, not a side channel.** §3 step 2 uses them to predict and
   remove the mature-contained unspliced. Without that, a highly-expressed short-exon gene's
   mature-contained unspliced would be mislabelled nascent (harmless for the gDNA prior, wrong for nRNA).
3. **It tells us exactly what's identifiable where** (§6) — and where it isn't (the floors).

## 6. The inference, and how every regime falls out

Poisson likelihood over the four counts; `m` pinned by the spliced reads, `ν` free per node, `ρ` smooth
across a gDNA-homogeneous segment. Profile `m` and `ν`; estimate `ρ` by ML over the segment (the weights
are the likelihood curvatures — derived, not invented, the point of the whole reframe).

- **Intron node:** no mature ⇒ two-species (gDNA + nascent); strand separates them.
- **Exon node:** all three; subtract mature (spliced), then strand-separate gDNA + nascent.
- **Intergenic:** gDNA only; the magnitude `U ~ Poisson(ρE_g)` pins `ρ` directly (no strand needed).
- **κ → ½:** `(2κ−1)→0` ⇒ the strand imbalance is uninformative ⇒ gDNA and nascent are confounded *after*
  mature removal ⇒ `ρ` rests on the intergenic magnitude (the floor). **Mature is still identified** (the
  spliced reads don't need strand). The degradation is continuous in κ.
- **AMBIG** (overlapping `+`/`−`): spliced is motif-oriented ⇒ `m₊, m₋` still identified ⇒ subtract both
  matures; the residual is `gDNA + nascent₊ + nascent₋` with two free nascent rates and only two counts ⇒
  `ρ` under-determined ⇒ inherit the segment, **unless** a neighbour pins `nascent₊` or `nascent₋` (the
  carry-over — now precisely "supply the missing equation").

## 7. Open / hard parts (unchanged by the species correction)

- **Capture** — `ρ` constant per segment fails under probe enrichment (exon `ρ` jumps). Orthogonal to the
  species deconvolution; the splice fraction / IPR effective length target it. Separate decision.
- **The `E_mu/E_ms` geometry** (§3 step 2) needs careful derivation — the contained-vs-spanning split of
  mature fragments at an exon as a function of exon length and the RNA FL. This is the new piece to get
  right; it is pure geometry (no new parameter).
- **Profile-ML vs full-Bayesian** for `ν`/`m`; the segment/smoothness model for `ρ`.

## 8. Bottom line

The model is now complete and correct: three species, four observables, all identifiable. The reassuring
part is that the **gDNA prior — the thing calibration must output — depends only on the unstranded-vs-
stranded split, which is invariant to the nascent/mature composition**, so the strand estimator is sound.
The completeness pays off in the *effective lengths* (the species have different geometry, esp. short
exons) and in enabling the nascent/mature split for nRNA. The spliced reads earn a new, active role:
predicting and subtracting the mature-contained unspliced before the strand separates gDNA from nascent.
