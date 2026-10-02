# Fragment length in Rigel's calibration: what is wrong, and what a real fix looks like

**Status:** an open design question, for review by other agents and people. Sandbox document (`docs/dev/`):
nothing here is authoritative until it lands in `DESIGN.md` / `EQUATIONS.md`.

**Date:** 2026-10-02.

**Owner's ruling of the same day:** fragment-length modelling, which had been deferred past the 0.8.0 release,
is back in scope. We want to solve it properly, not patch it.

**What we ask of you:** read §0–§2 for the problem and §9 for the candidate designs, then answer the questions in
§12. Simple and elegant beats clever. Every number below was measured on the current code against simulator truth;
§A says how to reproduce it.

---

## 0. Summary

Rigel separates a sequencing library into gDNA (contaminating genomic DNA) and RNA, then quantifies transcripts.
Its **calibration** stage estimates, for every genomic object (exon, intron, intergenic stretch, boundary), how
much of the object's fragment count is gDNA. That stage was designed and benchmarked on panels where **gDNA and
RNA fragments have the same length distribution** — deliberately, so that the transcript-level EM could not
separate the two origins on length alone and hide calibration bugs.

The consequence, found this week, is that calibration silently assumes equal lengths in several places. Real
libraries do not oblige: the VCaP DNA and RNA halves differ by 63 bp, and cfRNA is short against
nucleosomal cfDNA. When the lengths differ, two things break.

1. **A library with no hybrid capture is read as captured, and isoforms collapse.** Take a stranded, capture-OFF
   library whose RNA (75 bp) is shorter than its gDNA (250 bp). Transcript error is 42.1 %, against 3.8 % for
   the same library unstranded (§3). The cause is that short exons, which gDNA's long fragments can barely fit
   into, receive a sliver of gDNA from strand noise. That sliver, divided by a near-zero gDNA "opportunity",
   reads as a gDNA density 200× the library's. Calibration takes the resulting cluster as hybrid-capture
   enrichment and contracts every transcript's length, unevenly.
2. **Under capture, RNA's capture is priced at gDNA's efficiency.** Capture efficiency depends on fragment
   length: a molecule binds a probe through its overlap with it. Calibration measures one efficiency per object
   from gDNA and applies it to RNA. With RNA shorter than gDNA this misprices RNA's effective lengths (§6):
   transcript error is 28.5 % against 5.0 % when RNA is the long component.

A one-line guard fixes the worst case of failure 1: no object trains the gDNA density model unless its gDNA
opportunity admits at least one fragment position. Transcripts go from 42 % to 3.5 %, and the ladder is flat. It
is a band-aid. Calibration still over-attributes gDNA 3–7× wherever gDNA can only just fit (§4), and the guard
does nothing for failure 2. Separately, it exposes a fragility in how both-stranded objects read out their
composition (§5).

**The common root (§2):** calibration measures in each object's count frame (which fragments landed here) but
reasons in the density frame (fragments per admissible position). The bridge between the two frames is the
per-component **opportunity** `E_c` — how many places a fragment of component `c` could have sat. Under equal
lengths, `E_gDNA ≈ E_RNA` at every object, so the frames coincide and every mismatch between them is invisible.
On the benchmark ladder, the mass-weighted mean of `|log(E_g/E_r)|` is 0.004. On the length-gap panels it
reaches 6 nats at short exons.

---

## 1. Rigel in one page

| term | meaning |
|---|---|
| **object** / **slot** | A genomic unit calibration solves: a REGION (exon, intron or intergenic stretch, cut at every annotation boundary), a BOUNDARY between two regions, or a splice junction. About 35k regions on the panels' genome; millions genome-wide. |
| **contained count** `M` | Fragments lying entirely inside a region. The accumulator stores this per region and strand, **with no length information** (a deliberate choice under equal lengths: a length row "adds no tilt"). Boundaries also store `Σ 1/A(w)` (below). |
| **length law** `P_c(w)` | The fragment-length distribution of component `c ∈ {gDNA, RNA}`. Estimated once, **before calibration**, from structurally pure pools: RNA from spliced fragments, gDNA from intergenic and intronic contained fragments (§7). |
| **admissible positions** `A_o(w)` | The number of start positions at which a fragment of length `w` is contained in object `o`: `max(0, L − w + 1)` for a region of length `L`. |
| **opportunity** `E_c(o)` | `Σ_w P_c(w) · A_o(w)`: the expected number of placements of a `c` fragment in `o`. Layer 2 of calibration computes it exactly. A 97 bp exon has `E_g ≈ 0.09` for 250 ± 60 bp gDNA but `E_r ≈ 20` for 75 ± 20 bp RNA. |
| **density** `ρ_c(o)` | Fragments of `c` per admissible position: expected count `= ρ_c · E_c`. gDNA is genomically uniform before capture, so `ρ_g` is roughly constant across the genome (0.05 per bp on the panel below) except where hybrid capture enriches it. |
| **composition** `f_g` | gDNA's share of an object's contained count: `f_g = ρ_g E_g / (ρ_g E_g + ρ_r E_r)`. It lives in the **count frame**. |
| **strand channel** | In a stranded library RNA reads sense with probability `1 − κ` (κ ≈ 0.01) and gDNA reads either strand with probability ½. An object's strand counts are therefore evidence about `f_g`. |
| **ψ** | The per-object solver. It builds a posterior over `λ = logit f_g` on a grid from the strand likelihood, a **reference measure** (Jeffreys `½ log f + ½ log(1 − f)` on `f_g`), the gDNA **landscape** prior when fitted, likelihood rows from neighbours ("messages"), and, for introns only, a background-density likelihood (the "intron factory"). |
| **landscape** | A population prior on `log ρ_g`, fitted from objects whose composition is "located" (`Var(log f_g) ≤ 1 nat²`). It is refitted three times ("refits"); pass 0 has no landscape. |
| **capture reference** `ρ_ref` | The landscape's enriched mode, if one exists, read as "the fully captured gDNA density". Each object's **capture efficiency** is its gDNA density against `ρ_ref`, clipped at 1. No enriched mode means no capture, and every efficiency is 1. |
| **ruler** | Uses those efficiencies to contract each transcript's (and each gDNA component's) **effective length** inside the EM, so that capture-enriched fragments are not mistaken for abundance. |
| **EM** | Per-locus assignment of fragments to transcripts, synthetic nascent spans and a gDNA component. It is length-aware fragment by fragment (it scores each candidate's transcript-space length under `P_r`, and gDNA under `P_g`). Calibration hands it the per-locus gDNA pseudocount and the capture-contracted lengths. |

Constraints every proposal must respect:
- **Three populations only:** gDNA, RNA+, RNA− (Axiom 0).
- **No unexplained constants:** every number is derived and documented.
- **One mechanism per A/B.**
- **Judged per stratum:** unstranded × capture-OFF, stranded × capture-OFF, stranded × capture-ON; unstranded ×
  capture-ON is reported but deferred.
- **No harm to the equal-length ladder,** the 16-condition panel the tool is ranked on.

---

## 2. The diagnosis: calibration assumes equal lengths in four places

Write `r(o) = log(E_g(o) / E_r(o))`, the log ratio of gDNA's and RNA's opportunity at object `o`.

| panel | lengths, gDNA / RNA (bp) | `r` at the 1st, 10th, 50th, 90th, 99th percentile of regions | mass-weighted mean of `\|r\|` |
|---|---|---|---|
| ladder g50 ss.99 OFF | equal | −0.11, −0.10, −0.02, 0, 0 | **0.004** |
| fl-gap RNA-short g50 ss.99 OFF | 250 / 75 | −5.99, −5.56, −0.64, −0.03, 0 | 0.430 |
| fl-gap RNA-long g50 ss.99 OFF | 75 / 250 | 0, 0.03, 0.73, 5.05, 5.62 | 0.104 |

Everything below is invisible when `r ≡ 0`.

| name | where | what it assumes | what happens when `r ≠ 0` |
|---|---|---|---|
| coordinate | ψ's coordinate | Its docstring defines `λ = logit f_g = log ρ_g − log ρ_rna` | The true identity is `logit f_g = log(ρ_g/ρ_r) + r`. Every object's coordinate is shifted by its own `r`. |
| reference | ψ's reference measure | `½ log f_g + ½ log(1 − f_g)`, symmetric in the count frame | In density terms the reference is shifted by `r`, up to 6 nats. In a short exon, "neutral" in the count frame means "gDNA as dense as the exon's RNA" in the density frame. |
| read-out | the density read-out | `ρ_g = f_g · M / E_g`, with location judged on `Var(log f_g)` | A small, located `f_g` from strand noise near the RNA vertex is divided by a tiny `E_g`. The density precision is never asked. |
| capture transfer | capture efficiency | One scalar per object, measured on gDNA's lengths, applied to RNA | Capture is length-selective: probe overlap is bounded by fragment length, and so is reach (§6). |

A fifth, related assumption is that **each component has one library-wide length law.** Off capture this is
true and our estimates are accurate (§7). Under capture the realised law depends on placement: on-target gDNA is
longer than off-target.

---

## 3. Failure 1: a capture-OFF library reads as captured

**The library:** the genome-scale fl-gap arm `flgap_rna_short` — RNA 75 ± 20 bp, gDNA 250 ± 60 bp, 100 bp
reads, 50 % gDNA, capture OFF.

| | stranded (ss 0.99) | unstranded |
|---|---|---|
| transcript error | **42.1 %** | 3.8 % |
| gene error | 4.2 % | 0.3 % |
| mRNA pool error | −102,416 | −976 |
| synthetic nascent-span error | +95,350 | −5,081 |
| calibration error against the oracle | 1.9 % | 2.1 % |

The worst single case: HPS4's ENST00000699228.1 reads 0 fragments against a truth of 8,408.

### The mechanism, step by step

Each step is measured.

1. **Strand noise credits short exons with gDNA.** A short exon (median 97 bp) holds dozens to hundreds of short
   RNA fragments, about 1 % of them on the wrong strand (κ). Calibration's strand deconvolution credits a
   fraction of a fragment to gDNA. Near the RNA vertex the posterior median of `f_g` sits above zero by the
   width of the strand term.
2. **gDNA can barely fit there.** A 250 ± 60 bp gDNA fragment almost never fits in a 97 bp exon: the slots that
   end up near the false mode have a median `E_g` of **0.14 positions**, and 92 % have fewer than 1.
3. **The density is absurd.** 0.87 fragments over 0.09 positions is about **10 gDNA fragments per bp**, against
   the library's 0.05 (200×).
4. **These slots train the landscape.** They pass every training rule: they have a composition, and the
   composition is located (the read-out assumption).
5. **They form a false "enriched" mode:** 9.65 per bp, 370 members.
6. **The ruler reads that mode as capture.** Every object's capture efficiency becomes 0.0006–0.6 (median 0.0056)
   on a library with no capture.
7. **Every EM length contracts, unevenly.** Transcript lengths drop to 0.17–10 % of their plain value (median
   0.54 %). Isoforms of one gene contract up to 20-fold differently, and the most contracted absorb the shared
   fragments. The gDNA component's length contracts about 200×.

### Why the failure has its shape

- **Stranded only.** Unstranded, a short exon's composition has no own evidence, so it is not located and does not
  train the landscape.
- **Full depth only.** At 10 % depth too few slots are located to form a mode. Transcript error is then 9.3 % for
  both libraries.
- **Short RNA only.** When RNA is the long component, short exons hold gDNA's short fragments and not RNA's.

### The one-variable A/Bs that isolate it

All on the stranded condition at full depth, pinned and with fractional assignment.

| arm | transcripts | mRNA | nascent spans |
|---|---|---|---|
| shipped | 42.1 % | −102k | +95k |
| the EM's strand term set to ½ (calibration stays stranded) | 42.8 % | −118k | +107k |
| uniform EM warm start | 42.1 % | −102k | +95k |
| MAP instead of VBEM | 42.1 % | −100k | +93k |
| calibration and the EM both read the library as unstranded | 3.6 % | −1.1k | −5.2k |
| **no capture reference, nothing else changed** | **3.5 %** | −950 | −3.2k |

Also checked:
- the length laws are correct (RNA 78.5 ± 17.2 bp, gDNA 249.6 ± 59.9 bp, both matching truth);
- the reads align to their transcripts' exons and junctions;
- the tree before this week's strand changes gives the shipped row to the fragment.

### It is a general over-attribution, not one exon

Calibration's gDNA against certified truth, binned by gDNA opportunity, on the same condition:

| `E_g` (positions) | regions | median `E_r/E_g` | true gDNA | calibrated, shipped | calibrated, with the guard of §4 |
|---|---|---|---|---|---|
| < 1 | 5,317 | 229 | 94 | **1,640 (17×)** | 671 (7×) |
| 1–3 | 1,245 | 47 | 150 | 614 (4×) | 447 (3×) |
| 3–10 | 1,185 | 21 | 477 | 560 | 543 |
| 10–30 | 1,027 | 9 | 1,087 | 1,132 | 1,047 |
| 30–100 | 1,545 | 4 | 4,856 | 5,394 | 4,489 |
| ≥ 100 | 12,415 | ~1 | ~4.72 M | within 1 % | within 1 % |

The over-call scales with `E_r/E_g`, exactly as the reference and read-out assumptions predict. In fragments it is small (a few thousand). In
density, the currency the landscape and the ruler read, it is enormous.

---

## 4. The one-position guard, and why it is a band-aid

The prototype rule: **a counted slot trains the gDNA landscape only where `E_g ≥ 1`**, i.e. at least one admissible
position for a contained gDNA fragment. It follows the existing lesson that a region below one fragment length has
no measurable density (`TRAPS: density-below-one-fragment-length`).

| | shipped | with the guard |
|---|---|---|
| RNA-short stranded × OFF transcripts | 42.08 % | **3.54 %** |
| — genes | 4.20 % | 0.23 % |
| — mRNA pool error | 102,416 | 920 |
| — nascent-span error | 95,350 | 3,144 |
| ladder, every stratum | — | flat (calibration identical to two decimals, transcripts ±0.01 pt) |
| ladder zero-gDNA controls (false gDNA) | 126 / 335 / 127 / 183 | 120 / 334 / 120 / 176 |
| odg05 and four other fl panels | — | transcripts within ±0.15 pt |
| test-chromosome stranded × capture-ON gDNA pool | 57.2k | 63.3k (+10.7 %), see §5 |

The guard tells true capture apart from the false mode cleanly:

| row | slots near the reference | with ≥ 1 position |
|---|---|---|
| ladder g05 capture-ON | 4,528 | 4,214 |
| ladder g50 capture-ON | 4,864 | 4,474 |
| test g05 capture-ON | 352 | 352 |
| the false mode | 376 | 29 |

**Why it is not the fix.** gDNA's length is a distribution. A region with 1.5 positions is "allowed", yet its
density is still a sliver of strand noise divided by a small number: the table in §3 shows 3× over-attribution
at 1–3 positions with the guard on. The guard also does nothing about the composition itself (the coordinate, reference and read-out assumptions), which
still enters the EM's per-locus gDNA count, or about failure 2. It is a threshold standing in for information
calibration does not use.

---

## 5. Why the guard regresses one row (g98 ss0.70 capture-ON)

The test chromosome's 98 %-gDNA, weakly stranded, capture-ON row loses 5.7k gDNA under the guard (of 1.15 M true;
transcripts 11.1k → 16.5k). Equal lengths (206 ± 98 bp) hold on this row, so it is not the length mechanism.

- **The guard drops exactly one training slot of ~1,300.** It is a short probed exon whose one true captured
  fragment puts its kernel at ~12 per bp, above the 3.2 per bp enriched mode. The landscape's modes are unchanged
  to the third decimal, and so is the capture reference (3.1967 either way).
- **Dropping a random located slot instead moves calibration by ±30 fragments,** so calibration is not generally
  fragile.
- **The 5.7k sits in three ~10 kb objects (96 %), all both-stranded** (exons on both strands, so the strand
  channel cancels). Their true gDNA density is 0.31 per bp. Two of them flip from 0.30 to 0.035 per bp: their
  composition is decided by the landscape prior alone, read out as the median of a **multimodal** posterior, and
  one tail kernel tips the median from one mode to the other.

So the regression is a separate fragility: **a both-stranded object's composition is a median of a multimodal,
prior-dominated posterior, and can flip under a negligible perturbation.** It belongs on the list (§9, the both-stranded read-out), but
it is not caused by the guard's logic.

---

## 6. Failure 2: capture efficiency depends on fragment length

**The physics.** A fragment is pulled down by a probe through its overlap with it. In the simulator a fragment's
capture weight is `off_target_weight + binding_per_base · overlap`, where `overlap` is its best contiguous overlap
with a probe. Real hybridisation behaves similarly and saturates. Two consequences:
- **On target, overlap is bounded by fragment length.** A 75 bp fragment inside a 120 bp probe weighs at most 751;
  a 250 bp fragment spanning it weighs 1,201.
- **Near a probe, reach is bounded by fragment length.** A long fragment starting 150 bp from a probe still
  touches it; a short one does not. Partially probed and unprobed objects near probes are therefore captured
  far more for long fragments than for short ones.

**What Rigel does.** It measures one efficiency per object, from gDNA (§1), and applies it to every component's
effective length, so the ratio of RNA's capture to gDNA's is taken as 1. The true ratio at object `o` is
`E_{Q_g,o}[r] / E_{P_g,o}[r]` with `r(w) = P_r(w)/P_g(w)`, where `Q` is the captured and `P` the reference law.
That is exactly 1 only under equal lengths.

**Measured, capture-ON, stranded ss 0.99.** Transcripts: RNA-short 28.5 %, RNA-long 5.0 %.

`ruler_vs_truth.py` scores the EM's capture-contracted lengths against the simulator's own capture-aware lengths.
Values are the median log error per transcript class:

| class | RNA-short, shipped | RNA-short, ruler fed true gDNA | RNA-long, shipped |
|---|---|---|---|
| fully probed (≥ 90 %) mRNA | −0.01 | −0.00 | +0.00 |
| partially probed > ½ mRNA | +0.20 | +0.17 | +0.02 |
| partially probed ≤ ½ mRNA | +0.59 | +0.45 | −0.06 |
| unprobed mRNA | +2.96 | +1.62 | +2.13 (50 transcripts) |

RNA-short's lengths are wrong even when the ruler is given the true gDNA counts, so the formula is the problem,
not calibration's input. The pattern — exact on fully probed transcripts, worse the less probed — is the reach
mechanism. The equal-length ladder cannot show it: there `r ≡ 1`, and the ratio is 1.

A September experiment (the "uncaptured frame") fixed a related census-length inconsistency: on the ladder's
stranded × capture-ON rows, transcripts −4.5 % and genes −14.5 %. It lost on real gaps even at the simulator's own
laws (genes +92 % and +10.5 %), and the reason is this mispricing (the capture-transfer assumption). It was parked on 2026-09-30, and that parking
is now released.

---

## 7. How accurate are our length laws, and when could they be accurate?

**How they are made today.** Once, before calibration, from structurally pure pools:
- **RNA:** spliced fragments (certified RNA), corrected for the tilt that longer fragments cross junctions more
  often.
- **gDNA:** intergenic and intronic contained fragments, plus boundary-crossing pools with an inversion, shrunk
  toward the global law by an empirical-Bayes pseudocount of 1,000 (an underived constant).

**How accurate they are.** Estimated law against the simulator's realised lengths by origin. TV is the total
variation distance between the estimated and true laws.

| panel | capture | RNA, estimate | RNA, truth | TV | gDNA, estimate | gDNA, truth | TV |
|---|---|---|---|---|---|---|---|
| ladder | OFF | 204.1 ± 82.9 | 212.7 ± 86.0 | 0.046 | 216.7 ± 86.8 | 216.7 ± 86.7 | 0.003 |
| ladder | ON | 220.7 ± 80.2 | 228.9 ± 82.1 | 0.046 | 216.7 ± 86.7 | 240.8 ± 84.2 | **0.118** |
| RNA-short | OFF | 78.5 ± 17.2 | 78.4 ± 16.7 | 0.003 | 249.6 ± 59.9 | 249.6 ± 59.9 | 0.003 |
| RNA-short | ON | 81.5 ± 17.5 | 81.6 ± 17.1 | 0.003 | 249.8 ± 59.9 | 258.6 ± 58.0 | 0.060 |
| RNA-long | OFF | 242.1 ± 59.3 | 247.6 ± 59.8 | 0.040 | 78.6 ± 16.9 | 78.6 ± 16.8 | 0.002 |
| RNA-long | ON | 248.0 ± 57.5 | 252.8 ± 57.7 | 0.036 | 78.8 ± 18.2 | 81.7 ± 17.2 | 0.077 |
| test fl gDNA-long | OFF | 115.4 ± 44.0 | 113.3 ± 39.4 | 0.029 | 260.2 ± 111.3 | 260.7 ± 111.1 | 0.011 |
| test fl gDNA-long | ON | 126.3 ± 44.1 | 122.4 ± 39.9 | 0.038 | 265.8 ± 109.4 | 282.6 ± 106.0 | 0.093 |

Read it this way:
- **Off capture the laws are good.** gDNA is essentially exact. RNA reads 2–4 % short when RNA is long: a
  residual of the junction-crossing selection.
- **Under capture the gDNA law is the off-target law.** It misses the lengthening on target (ladder 216.7 against a
  realised 240.8), because there is no single gDNA law under capture: the realised law depends on placement.
  That is failure 2 again, seen from the length side.
- **On real data,** the spliced pool also carries splice artifacts (misaligned gDNA), which contaminate the RNA
  law. That is an open item.
- **At zero gDNA there is nothing to estimate the gDNA law from.** A law built from a handful of fragments must
  carry its uncertainty, or any mechanism reading it will "see" a gap that is not there. This is exactly how the
  2026-08-10 length channel failed (§8).

**"When do we finally have accurate laws — after calibration?"** Partly. Off capture, the pre-calibration laws are
already accurate enough. What is not accurate is a single global law where reality is per-placement (capture),
and the selection residuals (spliced tilt, artifacts). After calibration we know each object's composition, so a
fixed point becomes possible: re-estimate both laws from every fragment weighted by its calibrated origin,
on-target and off-target separately. But the fixed point only helps once the model has a place to put a
per-placement law (the capture opportunity and the fixed point). A better global estimate cannot fix failure 2.

---

## 8. What has been tried, and the bar a new design must clear

- **2026-08-10, a fragment-length composition channel.** A Gaussian row on per-slot moments (count, `Σ 1/w`,
  `Σ w`) in ψ. It was found inadmissible:
  - it reported **54–57 % gDNA on a zero-gDNA library**. With no gDNA to fit a law on, the two laws sat 1.2 bp
    apart, and the row still spoke at 100 % of slots;
  - its answer was not a function of the length gap. Shrinking the gap made it *worse*: a Gaussian
    log-likelihood is linear in the composition with slope ∝ gap, so its maximum is a grid endpoint;
  - its precision and participation were switches that fired at any nonzero gap, while its information faded
    continuously;
  - the region (node) length moments were off by 2–32 % against the opportunity model.

  The lesson recorded then: **whatever carries length information must shrink to exactly nothing as the gap
  closes, and must pass the zero-gDNA controls.** The channel was then deferred with the release; the deferral is
  now released.
- **A model-free local mean length** (`count / Σ(1/(w−1)) + 1 = E[w]` at a boundary) is exact, but inverting it
  for the composition is ill-conditioned.
- **2026-09-30, the uncaptured frame** (§6): wins with equal lengths, loses on gaps, because of the capture-transfer assumption.
- **2026-08-24, a reference *location*** (pulling ψ toward a fixed composition) was refuted and is banned. Where
  the strand carried no information, the location was the entire answer at any depth. Background information
  may enter only as a **likelihood whose precision scales with counts** (`DESIGN.md` §6b.1, §0c.0d).

---

## 9. Toward the real fix

### Principles

1. **One frame.** Every quantity calibration reasons about — prior, reference, location test, capture — is
   stated in the density frame, with each component's own opportunity as the bridge to the count frame. When
   `E_g = E_r` nothing changes, so the ladder is protected by construction.
2. **Length is information exactly where it is information.** A length-derived term must carry Fisher information
   ∝ (gap)², vanishing smoothly when the laws coincide and when a law is unmeasured. No switch, no threshold.
3. **Capture is a property of a fragment, not of an object.** Its length dependence is physics, and it belongs
   in the opportunity, where every component's length is already integrated.
4. **Smallest mechanism first, one at a time.**

### Candidate designs

**The background likelihood: give every object the opportunity-aware gDNA likelihood introns already have.
(Smallest; recommended first.)**

Introns get a background-density likelihood row: gDNA count `g ~ NegBinom(ρ_0 · E_g, α)` on the `f_g` grid,
with `ρ_0` the genome-wide gDNA density measured on intergenic DNA before calibration, and `α` its fitted
dispersion. Exons get nothing comparable in pass 0: only the strand term, the Jeffreys reference in the count
frame, and messages.

Extend the row to every object, in the shape already ruled for capture (`DESIGN.md` §0c.3):
- a **spike** at `ρ_0 · E_g`;
- a **slab** of capture enrichment, bounded by the neighbouring boundary's density below and the object's own
  total above;
- mixing weight ½ from the reference, so no new constant.

Off capture the slab collapses onto the spike (measured off-target-to-boundary ratio 0.98).

- **Why it fixes failure 1 without a threshold:** at `E_g = 0.09` the spike expects 0.004 gDNA fragments, so the
  data's own precision rules out a sliver of strand noise; at `E_g = 1.5` it expects 0.08. The fit is continuous
  in `E_g`, and `E_g` already integrates gDNA's whole length distribution, including the short tail that does
  fit.
- **Why it is allowed:** it is a likelihood whose precision scales with counts, which §6b.1 permits; it is not a
  location.
- **Why it is elegant:** one rule, `g ~ NegBinom(ρ · E_g)`, for every object; it is already built for introns
  (`native/solve_kernel.cpp`'s `factory_row`).
- **What it does not fix:** failure 2.
- **Risks:**
  - under capture the slab must not suppress true enrichment, which matters most where the strand channel is weak;
  - unstranded exons gain an own-evidence row for the first time, which moves the unstranded strata.
- **A cheaper alternative inside the same idea, the reference assumption alone:** shift the Jeffreys reference into the density frame
  by `r`. It is necessary for consistency but probably insufficient: a density-frame reference still calls "gDNA
  as dense as the exon's RNA" neutral, which for a highly expressed short exon is still 50× the library's gDNA.

**The length likelihood: length as a proper likelihood, built right.**

Store a coarse per-region length histogram (`B` bins), so the contained fragments of a region are Poisson by bin:
`n_{o,b} ~ Poisson(ρ_g E_{g,o,b} + ρ_r E_{r,o,b})`, with `E_{c,o,b} = Σ_{w∈b} P_c(w) A_o(w)`. This is the exact
likelihood, not Gaussian moments:
- its composition information is identically zero when `P_g = P_r` (the bins become proportional), and it grows
  with the gap squared;
- marginalising each law's own sampling uncertainty (the pool counts behind each `P_c`) makes the zero-gDNA case
  carry no information, which is the 2026-08-10 failure.

It also gives unstranded gap libraries an own composition channel, which today they lack.

- **Cost:** an accumulator schema change of `B` floats per region per strand (35k regions here; a few million
  genome-wide). The scan is single-pass, so the bins must be fixed at scan time (log-spaced, say), not after
  the laws are known.
- **Open:** does the region opportunity model hold up per bin? The 2026-08 node moments were 2–32 % off.

**The capture opportunity: capture as a length-resolved opportunity.**

Replace the per-object scalar efficiency with a capture-weighted opportunity per component:
`E^cap_c(o) = Σ_w P_c(w) · Σ_{x ∈ placements} κ(overlap(x, w))`, so capture's length dependence enters exactly
where the length laws already are. Three ways to get `κ` and the overlaps:
- **(a) With the probe BED.** Kits ship one, and Rigel could accept it as optional input, like the
  splice-artifact blacklist. Geometry is then exact, and only the response curve `κ(overlap)` (affine and
  saturating, 1–2 parameters) is fitted, from gDNA's on-target against off-target densities.
- **(b) Without a BED, separable.** `c_o(w) = c_o · h(w)`, with one library-wide response `h(w) ∝ Q_g(w)/P_g(w)`
  read from gDNA's captured law (on-target pools) against its off-target law. RNA's capture at `o` is then
  `c_o · E_{P_r}[h] / E_{P_g}[h]`. This needs no accumulator change, but misses reach at partially probed
  objects, which is where most of today's error sits.
- **(c) Without a BED, learned per object** from the length likelihood's histograms. This is the "missing measurement" the
  September review named: gDNA's length-resolved mass on RNA's placement objects.

**The fixed point: length laws re-estimated after calibration.**

After calibration, re-estimate `P_g` and `P_r` from every fragment weighted by its calibrated origin (and, under
capture, per placement class), then re-solve once. This only makes sense together with the length likelihood or the capture opportunity, which give the
per-placement laws somewhere to go. It also offers a principled replacement for the empirical-Bayes pseudocount
(1,000) of §7: weight each pool by its own precision.

**The both-stranded read-out (§5; independent of length).**

A both-stranded object's composition is the median of a multimodal, prior-dominated posterior, and it flips
under negligible perturbations. Options:
- read out the posterior mean;
- carry the posterior's width or modes downstream, so the EM's gDNA count reflects the uncertainty;
- demand that a flip be earned by evidence.

This needs its own A/B; we record it here because the guard exposed it.

### Recommended order

1. **The background likelihood**, off capture first. Measure whether it alone flattens §3's over-attribution table across every `E_g`
   bin, and whether the one-position guard then becomes inert. If it does, the guard is not needed.
2. **The capture opportunity, separable form (b)**, as the cheapest capture × length correction, scored on `ruler_vs_truth.py` by class on both gap arms
   and the ladder. Then form (a) if Rigel takes a BED, or form (c) if the length likelihood is built.
3. **The length likelihood and the fixed point**, if unstranded gap libraries or the residual composition need an own length channel. Pass the
   2026-08-10 bar first: no information at zero gap or zero gDNA.
4. **The both-stranded read-out** on its own.

---

## 10. How a fix will be judged

**Panels**, each read per stratum:
- the ladder (equal lengths; must stay flat; zero-gDNA controls must not rise);
- the test chromosome's three fl arms (gDNA long, RNA long, equal-200 control; 30 conditions each);
- the genome-scale gap arms (RNA short, RNA long);
- the test panel and odg05;
- the real VCaP mix, whose halves differ by 63 bp and which has per-fragment read-name truth.

**Gaps in the test substrate, to close first:**
- the genome-scale gap arms are all g50 and carry a retired nascent model. One gDNA level nearly let the 2026-08
  channel land on a false 87 % win, so add g00, g05 and g98 rungs and re-simulate on the current nascent model;
- add at least one gap that is short against the read length in **both** directions.

**Instruments:**
- `calibration_vs_oracle.py` for calibration, plus §3's over-attribution-by-opportunity table;
- `ruler_vs_truth.py` by probed class under capture;
- `quant_accuracy.py`, pinned and fractional, for transcripts, genes and the gDNA, nascent-span and mRNA pools;
- a **gap sweep**: interpolate `P_g → P_r`, and require every length-derived term to fade continuously to
  exactly zero.

**Gates:**
- zero-gDNA rows: no false gDNA;
- information vanishing at zero gap;
- no new unexplained constants;
- the two failures fixed without a threshold.

---

## 11. What is not in question

- **Three populations, always** (gDNA, RNA+, RNA−). Length does not create a fourth. "Mature" and "nascent" are
  not calibration populations: RNA inside an intron is RNA that has not spliced there.
- **The EM is already length-aware** fragment by fragment. This is about calibration and the ruler.
- **The ladder keeps equal lengths.** It remains the panel that prices calibration without length. Gap mechanisms
  are judged beside it.

---

## 12. Questions for reviewers

1. **Is the diagnosis complete?** Are there places beyond the four assumptions of §2 where calibration or the ruler assumes equal
   lengths?
2. **The background likelihood:** is extending it to every object (spike-and-slab under capture) the right
   minimal fix for failure 1? What should the slab be, so that true capture enrichment survives where the strand
   channel is weak?
3. **The background likelihood against the reference:** should the Jeffreys reference itself move to the density frame as well, or is the likelihood
   enough? Is there a principled reference in the density frame?
4. **The length likelihood:** what is the right per-object representation of length evidence for a single-pass scanner — fixed bins,
   a few moments with an exact likelihood, or something else? How should the uncertainty of an estimated law be
   propagated, so the information vanishes at zero gap and at zero gDNA?
5. **The capture opportunity:** should Rigel accept a probe BED? Without one, is the separable response `h(w)` good enough, or is the
   per-object reach essential? What capture response curve is defensible from first principles?
6. **Laws:** off capture our pre-calibration laws are accurate. Is a post-calibration fixed point worth its
   complexity, and how should placement-dependent laws under capture be represented?
7. **The both-stranded read-out:** for a prior-dominated, multimodal composition posterior, what should a calibration read out — median,
   mean, the full posterior?
8. **Simplicity:** is there one mechanism that fixes both failures at once? A length-resolved opportunity with
   capture inside it (the length likelihood plus the learned capture opportunity) is the candidate we see. Is there a simpler one?
9. **Testing:** what simulation panels or real-data checks would you add to make a length mechanism convincing?

---

## A. Reproducing the numbers

Everything is under `~/Downloads/rigel_runs/prototypes/2026-10-01_shortrna/`.

| file | what it shows |
|---|---|
| `VERDICT.md` | §3–§6 in brief |
| `quant_sub.py` | in-process, pinned run with read-name truth; one-variable arms: `arm=em_unstranded`, `all_unstranded`, `no_reference` |
| `score.py` | read-name scoring |
| `training_census.py` | who trains the landscape, and their opportunity |
| `opportunity_bins.py` | §3's over-attribution table; runs under the current tree or a snapshot |
| `g98_dissect.py`, `drop_one.py` | §5 |
| `fl_accuracy.py` | §7 |
| `snap_proto/` | the guard's snapshot tree; `run_snapshot.sh` runs any tree outside the editable install |
| `proto/` | the guard's panel A/B: `compare.py`, `compare.txt` |
| `ruler/` | §6's `ruler_vs_truth.py` tables |

**Panels:**
- `~/Downloads/rigel_runs/suite/{ladder,flgap_rna_short,flgap_rna_long}`;
- `~/Downloads/rigel_runs/test_reference/scenarios*`;
- each has scan and oracle caches, and `slot_truth.npz` per condition, holding per-object true counts and both
  opportunities.

## B. Symbols

| symbol | meaning |
|---|---|
| `P_c(w)` | the length law of component `c` |
| `A_o(w)` | admissible contained positions for a length-`w` fragment in object `o` |
| `E_c(o) = Σ_w P_c(w) A_o(w)` | the opportunity |
| `ρ_c` | the density, so that the expected count is `ρ_c E_c` |
| `f_g = ρ_g E_g / (ρ_g E_g + ρ_r E_r)` | the count-frame composition |
| `r = log(E_g/E_r)` | the frame offset |
| κ | RNA's wrong-strand rate |
| `ρ_0` | the off-target gDNA density |
| `ρ_ref` | the capture reference |
| `Q_c` | the captured length law; `P_c` is the reference law |
