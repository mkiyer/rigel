# Effective length under hybrid capture — the problem, stated (working note, 2026-09-15)

A working note in the sandbox: the problem the EM's ruler solves, why it exists, what the shipped
estimator does and why, and what changes if the ruler is read from total abundance instead of gDNA. It is
written to be argued with. Nothing here is a ruling; the rulings it leans on are `DESIGN.md` §7.2 and
`EQUATIONS.md` §11, and the defect that prompted it is `ISSUES: the-ruler-reference-on-sparse-real-libraries`.

## 1. The problem

Write `c(x)` for the capture efficiency at genome position `x`: the probability, up to a constant, that a
fragment starting at `x` survives capture. Off a probe it is `c_off`; on a probe it is `c_on`, two to
three decades higher; at a probe's edge it tapers over a fragment length. In the simulator this is
`CaptureSampler.fragment_weight`, `off_target_weight + binding_per_base · overlap`.

A library holds gDNA at a uniform start rate `g` per base and transcripts `t`, each with a start rate
`a_t` per base of its footprint `F_t` (its exons, or its span for an unspliced entity). The fragment
intensity at `x` is

    n(x) = c(x) · ( g + Σ_{t ∋ x} a_t ).

The EM turns the fragments it assigns to `t`, `N_t`, into an abundance by dividing by a length:

    â_t = N_t / L_t^eff,      L_t^eff = Σ_{n ∈ F_t} S_n · c̃_n,      c̃_n = c_n / c_on,

over the objects `n` the transcript occupies with supports `S_n`. `L_t^eff` is the length the transcript
actually offered to the library. Only the ratios `c̃_n` matter: a global scale on every effective length
is a global scale on every abundance and cancels in any normalisation. **The estimand is a property of
the panel and the protocol.** It is the same at 50 % gDNA and at zero gDNA. What changes with gDNA is
only what we can see of it.

## 2. Why the lengths must shrink

Without the correction, `c̃ ≡ 1` and

    â_t = N_t / L_t = a_t · c̄_t,      c̄_t = Σ S_n c_n / Σ S_n,

so every transcript is under-estimated by its own mean efficiency: an unprobed transcript by `c_off/c_on`
against a probed one, a transcript with a long unprobed UTR by the UTR's share of its length. That is the
bias the owner names — transcripts with long unprobed regions read low against transcripts with better
probe overlap.

The second reason lives inside the EM. The locus's gDNA component spans introns and exons; under
capture its fragments sit in the probed exons. A uniform component at full span cannot explain that
concentration, so the EM hands the exonic gDNA to the mRNA, which can. Contracting every component by
the same `c̃` — the gDNA component to its exonic footprint, the mRNA hardly at all — is what lets the two
compete on one scale (`capture_eff_length`'s module docstring; `assemble_priors` applies the same rule to
the locus prior).

## 3. The current approach

Per object, the calibration's deconvolved gDNA density `ρ^g_n = m_n / S_n`, and

    c̃_n = min( ρ^g_n / ρ_ref , 1 ),      factor_t = Σ S_n c̃_n / Σ S_n,      eff_t = fl_t · factor_t,

with `ρ_ref` the fully-captured gDNA level — since 2026-09-14 the located enriched mode of the fitted
gDNA landscape, `None` when there is none, and `None` means `factor ≡ 1`. A multimapper floor
`w = C/(C+1)` on the object's contained evidence `C` pulls the factor toward 1 where little was counted.
A spliced transcript stitches its junctions in at the flanking exons' density so its uncontracted length
equals its fragment-length marginal.

## 4. Why it is justified, and exactly where the justification ends

gDNA is one template at a uniform rate, so `ρ^g_n = c_n · g`: its density **is** the efficiency map,
with nothing of the expression dynamic range in it. With `ρ_ref = c_on · g`,

    c̃_n = min( c_n g / c_on g , 1 ) = c_n / c_on,

exactly, and `â_t = a_t`. The cap makes the fully-captured level the unit, keeps `eff ≤ fl`, and clips
Poisson noise above the reference (the plug-in on a noisy field is biased below 1 by Jensen plus the
clip, which is why a reference read per object contracted a capture-OFF field by 8 %; `EQUATIONS.md` §11).
The reference is a population statement so that no single object's noise sets the unit.

The justification needs two things the sample must supply. **Per object, a measurement of `ρ^g_n`**:
`E[k_n] = c_n g S_n` gDNA fragments, and a region below one fragment has no location
(`Var(log k) = 1/k`, the one-fragment floor `_LOCATED_VAR` is derived from). **As a population, a located
enriched mode**: about `√n_train` probed regions at or above one fragment, or the mode's kernels render
decades wide and the verdict is `None`. Both are statements about `g · c_on · S`, gDNA fragments per
probed region. On the ladder at 5 % gDNA a probed exon holds tens; on LBX0190, with 3,434 gDNA fragments
across the genome, nearly every probed exon holds none, 1,118 regions of a million hold one, and the
gDNA that is captured sits in a few dozen hyper-efficient regions (real efficiency is a spectrum, not two
levels). There the per-object measurement is gone and the population mode is a handful of kernels: the
witness has faded, and the shipped rule then reads walls as a mode (`ISSUES:
the-ruler-reference-on-sparse-real-libraries`). At `g = 0` there is nothing to read at all.

## 5. The zero-gDNA case: what RNA alone can say

Set `g = 0`. Then `n(x) = c(x) · a(x)`, two unknown non-negative fields and one observation.

**The depleted mode.** Intronic and intergenic regions carry no probes and, with no gDNA, no fragments
either. They read zero: the anchors' wall, not a level. There is no depleted level to measure at `g = 0`;
what there is, is a floor. (At any `g > 0` they measure `c_off · g` precisely, and that is why the
depleted mode has never been the problem.)

**Exonic abundance.** An exon's density is `c_n · a_n`. High density is high expression, probed or not;
low density is low expression, probed or not. On a tiled panel `c` is constant over a transcript's
footprint and so is `a_t`, and then **only the product `c_t · a_t` is identifiable per transcript.** No
estimator escapes this; it is not a matter of cleverness. The only RNA signal that separates the two is
within-transcript contrast, a footprint that spans probed and unprobed sequence at one `a_t` — absent on
tiled panels and in cell-free RNA, and set aside by the owner.

**The population.** Over exons, `log ρ^tot = log c + log a`. If capture and expression were
independent, the observed distribution would be the expression distribution copied twice, shifted by
`log(c_on/c_off)`, and the shift would be estimable from the two modes. Panels are designed, not drawn:
the probed genes are chosen for a reason and their expression differs from the rest by an amount `Δ`
of either sign, and the modes then sit `log(c_on/c_off) + Δ` apart. What a total-landscape mode measures
is the typical probed exon's density, `c_on` times a typical probed expression, and its distance from the
other mode is not the enrichment ratio. Enrichment and expression enter the same observation additively
in the log and nothing in one sample's RNA tells them apart.

## 6. Total abundance in place of gDNA: what changes

**As the ruler's density.** Put `ρ^tot_n = c_n (g + a_n)` where the rule has `ρ^g_n`, with `ρ_ref` the
enriched mode of the total landscape. For a transcript alone in its locus the EM's count is its own
density integrated over its footprint, `N_t = Σ ρ^tot_n S_n`, so

    â_t = N_t / (L_t · factor_t) = Σ ρ_n S_n  /  Σ S_n · min(ρ_n / ρ_ref, 1).

If every object of `t` lies below the reference, `â_t = ρ_ref`, a constant. If every object lies above
it, `â_t` is the uncorrected mean density. **A ruler read from the same signal it corrects is a clamp:
every abundance below the reference is raised to the reference, and every abundance above it is left
uncorrected.** The multimapper floor only moves a transcript toward the uncorrected value where little
was counted, so the clamp bites hardest on well-measured, lowly expressed transcripts — the opposite of
what a correction is for. The owner's sketch is this result exactly: low-density regions get decreased
support and are inflated, up to precisely the reference; high-density regions keep their full length.
Compare the gDNA ruler, whose density `c_n g` does not contain `a_t`, so the same algebra returns `a_t`.

**As geometry only.** Use the total landscape to class each object probed or unprobed and give the
classes `c̃ = 1` and `c̃ = r = c_off/c_on`. The bias moves from every transcript to the misclassified
ones — a highly expressed unprobed exon read as probed, a silent probed exon read as unprobed — and each
costs `1/r`, two to three decades. And `r` itself is not identifiable from RNA (§5): it must come from
gDNA, or be assumed, or be read off the modes' separation with `Δ` folded in.

**As a bound.** Since `ρ^g_n ≤ ρ^tot_n`, `c̃_n ≤ ρ^tot_n / ρ_ref` with the gDNA reference. Sound and
one-sided; vacuous for every probed object (`ρ^tot ≥ ρ_ref`) and at `g = 0` (`ρ_ref → 0`); informative
for an unprobed object whose RNA is not much denser than its gDNA. This is the one way total abundance
enters without confounding, and it is a refinement of the gDNA ruler, not a replacement.

## 7. What happens, by regime

| gDNA per probed region | gDNA ruler (shipped form) | total-density ruler |
|---|---|---|
| tens or more (the ladder at 5 % and above) | exact up to the efficiency spectrum; verified on both panels | clamps every abundance below the reference |
| about one | the reference is located but each object's plug-in density is noise; the cap and the floor are doing the work | clamps |
| well below one, `g > 0` (LBX0190, MO_3021) | the population mode is a few kernels; the shipped membership test reads walls; the honest verdict is `None` | clamps, and now the total landscape is all RNA |
| zero | nothing to read; `None`; the bias of §2 stands undeclared | clamps on pure expression: the correction is expression itself |

## 8. A better idea, from the theory rather than the table

The principle the algebra of §6 gives is short: **a ruler must be read from a signal that does not
contain the quantity it corrects.** gDNA is that signal. Total abundance is not, and cannot be made one
by a cleverer detector. So the design question is not which landscape to read but how far down the gDNA
axis the gDNA witness can be made to reach, and what the tool says when it cannot reach.

* **Read the object's efficiency as an expectation, not a plug-in.** The landscape is already a prior
  over `log ρ^g` fitted from the population with zero-native Poisson kernels and an E-step for the
  location-free ones. The ruler currently takes the point estimate `m_n / S_n`. The principled quantity
  is the posterior mean of the clipped density under that prior, `E[min(ρ/ρ_ref, 1) | k_n, S_n]`: at high
  depth it is the plug-in, at low depth it is the population's own mixture weighted by the object's
  likelihood, and no constant enters. This is where the low-gDNA regime's information actually sits — in
  the population, not in the object — and it extends the working range down to wherever the population
  still locates its enriched mode.
* **Locate the mode on located kernels.** The membership defect of the open issue: an anchor's wall is
  not a location, and `√n` located kernels is the population's own resolution. With that, the boundary
  of the working range is a derived statement — about `√n_train` probed regions at one fragment or more
  — that the tool can evaluate and report on every run.
* **Declare the regime.** Below the boundary the tool cannot correct, and it should say so in the result
  and the report — the reference, the number of located kernels it rests on, and that the abundances
  carry the bias of §2 — rather than fabricate a map from RNA. Uncorrected-and-declared is what every
  other capture pipeline delivers; it is the honest floor.
* **Use the total only as the bound of §6**, where it sharpens an unprobed object's efficiency without
  confounding.

Two experiments settle the theory before any of this is built, both on the test chromosome with the
simulator's own `CaptureSampler.partition_array` as the truth. The clamp is falsifiable in one run: the
total-density ruler on any capture-ON row must collapse the abundances of every transcript below the
reference onto one value. And the regimes of §7 are a ladder — gDNA at 0, 0.1, 1, 5, 25, 50 % and depth
at full, a tenth, a hundredth — on which the shipped ruler, the expectation ruler and the total ruler are
scored per transcript class against the truth, up to a global scale. Where the gDNA rulers fall off that
ladder is the boundary the tool must declare.

## 9. Measured (2026-09-15, the same day; instruments under `~/Downloads/rigel_runs/prototypes/2026-09-15_ruler/s4/`)

The truth is the simulator's `CaptureSampler.partition_array`, fl-marginal over the pre-capture length pmf;
every ruler is scored as log(estimate / truth) anchored on the fully probed transcripts, per class.

**The formula is right where its witness is.** Test chromosome, 5 % gDNA, stranded, capture-ON: probed
transcripts 99 % within ±0.1 nat, partial 60–100 %, for the shipped ruler, the ruler without its floor, the
expectation ruler and the ruler on the certified true gDNA counts alike.

**The floor is a 30–40× bias on unprobed transcripts.** Shipped +3.42 nat at 5 %, +3.48 at 25 %, +3.71 at 50 %:
`1/(C+1)` caps the contraction. Without the floor: +0.12 / +0.04 / +0.09. The expectation ruler: −0.001 /
+0.01 / +0.05. The multimapper floor was protecting against a factor of exactly 0; the posterior mean under
the landscape prior has none.

**Total density fails as §6 derives.** "Largest mass above the depleted" finds the unprobed-expression level
and corrects nothing (+7.0 on unprobed transcripts); the highest located basin puts the reference on the
probed level and the probed class then scatters ±1.5 nat with the error sloped +0.4 to +0.75 per nat of true
expression — the clamp's signature; on capture-OFF the total ruler contracts a field that must read 1.

**The depth ladder** (0, 0.1, 1, 5, 25, 50 % gDNA × full, tenth and hundredth depth, `scenarios_depth_*`): the
reference is `None` below ~120 gDNA fragments and the tool says so; probed transcripts read within ±0.1 nat
for 30 % at ~1,200 fragments, 80–90 % at ~6,000, 92–98 % at ~12,000, all of them from ~30,000 — the axis is
gDNA fragments per probed exon piece (3 / 17 / 35 / 85 here), the same curve at every depth. The expectation
ruler lifts the low end (91 % vs 65 % at 5,850 fragments at a hundredth of the depth; 98 % vs 92 % at 11,700
at full depth) and holds the unprobed class to +0.4 nat at 1 % gDNA where the shipped ruler reads +3.5.

**The ladder shows the witness's geometry.** On `g05 ss.99 ON` the probed transcripts scatter −0.5 to +0.3
nat under every gDNA ruler, the certified true counts included: the ladder's panel is designed in transcript
coordinates, its probes span junctions, and gDNA at a split probe is captured at a fifth of the cDNA's weight
in the simulator. Its unprobed transcripts hold no gDNA fragment yet read +2.0 nat without the floor: the
calibration assigns gDNA to exons that have none, a composition residual the ruler inherits.

**The junction rule is part of the repair.** A junction is imputed from the pieces beside it; on the ladder those
are often shorter than a fragment with no support, and the flank then reads 0, fully captured, or the prior's mean
(0.096), so forty junctions at the crossing support swamp four hundred bases of measured regions — the other half
of the ladder's probed-class scatter, on the estimator's side. Four rules measured (probed within ±0.1 / unprobed
median on `g05 ss.99` and `g50 ss.50`): flanks 36 % / +2.0, 24 % / +0.3; exon pools 41 % / +1.7, 25 % / +0.2;
flank-where-supported-else-pool 41 % / +1.9, 27 % / +0.2; the transcript's support-weighted pool 47 % / +0.8,
42 % / +0.3 (and under the expectation ruler 54 % / +1.5, 38 % / +0.4). The transcript pool survives the dense
annotation at a small cost on the test chromosome's designed half-probed transcripts (60 % → 33–40 % of 15 within
±0.1), whose junctions truly read the flanks' mean; the reconciling derivation — the flank's evidence shrunk
toward the transcript's pool by its own support — is the next session's.

**What to build.** The expectation ruler in place of the plug-in and the floor, with the junction rule, in
`capture_eff_length` and `assemble_priors` alike; membership of the reference's mode on located kernels; the regime declared on the
result. Gated by the truth instrument on both panels and the depth ladder, and by the four real libraries as
test inputs. Recorded as `ISSUES: ruler-multimapper-floor-caps-the-correction`,
`ISSUES: ruler-witness-geometry-on-transcript-panels`, `ISSUES: the-ruler-reference-on-sparse-real-libraries`.
