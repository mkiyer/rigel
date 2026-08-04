

# summary

the calibration module has a solver that utilizes belief propagation to impute the solution across nodes and edges along contiguous genomic intervals.

each node/edge solves using
- internal information (strand model, density model, eventually the new fragment length model)
- messages (beliefs communicated from neighbor nodes)
- gdna prior (does not exist for pass-0)



# the challenge of density nonuniformity

we assume that gDNA is generally uniform across a locus (a locus is a contiguous region of the genome).

the challenge of solving relates to density nonuniformity. let's start with capture-OFF (no hybrid capture).

## capture-OFF sources of unspliced density nonuniformity

- unspliced density = gdna density + rna density
- gdna is uniform (our base assumption)
- rna is NOT uniform

let's look at what happens at boundaries.

take a simple transcript:
- TA+ exons (1000, 2000), (3000, 4000)

intergenic node (0, 1000): 
- uniform gdna density
- no RNA present (structural ZERO)

intergenic-exon edge @ 1000:
- uniform gdna density
- structural start (TSS)
- RNA active on + strand (RNA+ active, RNA- NOT active)
- here, we have unspliced density nonuniformity
- the absolute density of the intergenic node can be impute onto the intergenic-exon edge (density of gdna is uniform)
- the density of the RNA+ transcript is unknown

exon+ node (1000, 2000):
- gDNA and RNA+ active

exon-intron edge @ 2000:
- RNA+ has TWO paths
- RNA+ can *SPLICE out* using the splice junction (intron 2000-3000). the splice junction density is measured. we have a splice graph with splice junction fragment counts and densities
- RNA+ can continue contiguously (unspliced RNA)
- gDNA is uniform

intron+ (2000,3000):
- gdna and unspliced rna can be present
- density model used to deconvolute introns

intron-exon edge @ 3000:
- RNA+ comes from TWO sources
- RNA can *splice in* from the intron (2000,3000). this is the splice graph, splice junction fragments measured. need to lookup the incoming splice junctions at edge
- RNA can continue contiguously (unspliced)

exon+ (3000-4000)
- gDNA and RNA+ active

exon-intergenic edge @ 4000:
- structural END (TES)
- RNA does not continue
- gDNA continues


So in SUMMARY:
- transcript start sites (TSS) and end sites (TES) MAY density nonuniformity if RNA components start and stop there. transcripts may have multiple TSS And TES sites (different isoforms can use different TSS/TES)
- RNA can *splice out* and *splice in* at splice junctions. Splice junctions are a subset of edges of the graph. Spliced RNA is directly measured (observed) from the data. Not solved for.
- Unspliced RNA competes with gDNA

Traps:
- Deriving and solving the capture-OFF case is somewhat trivial. The logic completely fails under capture-ON.
- A more generalizable model is needed that accounts for hybrid capture on AND off.




## capture-ON sources of unspliced density nonuniformity

Hybrid capture uses probes targeting specific sequences. Those sequences are enriched (pulled down). The off-target molecules are depleted.

**Capture-ON ADDS density nonuniformity at captured nodes.**

**Assumption: exons are captured, introns/intergenic regions are not captured**

We have a sophisticated model derived, tested, and partially optimized for hybrid capture that is based on enrichment ratios between adjacent nodes/edges.

Enrichment ratios are the ratio of densities between adjacent nodes/edges.

For example:
- An exon has density 10 <-> The adjacent exon-intron edge has density 5 <-> The adjacent intron node has density of 1.0
- Instead of absolute densities, the enrichment/depletion induced by hybrid capture transforms the problem into a composition (fractional density problem).
- For certain node-edge pairs, we prove that we composition (the fraction of RNA+, RNA-, and gDNA) is unchanged.

For example, a single-stranded RNA+ intron with density 1 and the adjacent intron-exon node with density 100.
- The SAME RNA strands are active (RNA+)
- The only change is the potential for capture at the neighboring EXON. If the exon is captured, the probe may pull down some of the intron-exon edge crossing fragments, enriching them over background levels.
- The intron is deconvoluted relative to intergenic levels, which are assumed to be at background.
- So between an intron and an intron-exon edge, we can impute the composition as the message sent. The intron sends its composition and its densities to the intron-exon edge.
- The solver 'reframes' densities based on the enrichment ratio. Once the source (intron) and destination (intron-exon edge) are normalized by enrichment ratios (the 'reframe'), the solver can proceed.

The complexities arise for:
- intron-exon edge <-> exon node
- transcript start and stop (TSS and TES)

We built the new accumulator to better handle these cases. Now it is a matter of finishing the implementation and debugging it.


# what is left to do

## re-derive and prove the solver using the new accumulator

- utilize TSS and TES information (previously absent)
- prove and improve belief propagation arithmetic (derivation) using new accumulator data
- dissect and address one small issue at a time (for example, handle intergenic-exon edges, then handle intron-exon edges, etc)

## the path to success

- start with single-stranded transcripts, capture OFF and ON, derive the belief propagation arithmetic, the solver equations
- this is a bipartite graph (node <-> edge <-> node <-> edge) with a new splice junction data structure (at edges)
- we need to be able to lookup the *splice in* and *spliced out* counts/densities at edges
- we need to prove to ourselves the arithmetic works in isolated conditions.


# possible helpful additionals

## the toy harness

- we don't yet have the capability to "piggyback" on a cached scenario and add addition "spike in" transcripts to test additional behaviors
- using a cached scenario means a lot of parameters are "locked", including the FL distributions, the strand specificity, hybrid capture behavior, gDNA level, etc.
- but under those existing parameter constraints, we could add additional test transcripts to a "toy chromosome", simulate fragments (using the same simulation configuration as the matched cached scenario) and study the behavior of the test transcripts.

---

# ⭐ MEASUREMENT ADDENDUM (session 2026-08-03) — what the ladder says about the above

Added below the owner's derivation, not into it. Every number is from the 36-condition gDNA ladder with
the shipped oracle cache, pass-0 (prior-free), scored per object against the origin-split oracle.
⚠ The dropped/confirmed verdicts and the re-ranking live in `ROADMAP.md` §2; this section is only the
part that bears on the **derivation**.

## 1. ⭐⭐ The reframe's premise fails on a COMPONENT SET, and that is a rule, not a case

The doc's model of the reframe — "normalize source and destination by enrichment ratios, then impute the
composition" — is right, and there is one condition on it that is currently unstated and violated:

> **`r = ρ_tot(dst)/ρ_tot(src)` is an enrichment ratio only if `ρ_tot` accounts for the SAME component
> set at both ends.**

At a line, `ρ_tot` includes the junction channel; at a NODE it structurally cannot (an exon's mature
molecules live inside its *contained* count, never in a junction bank). So every `node → line` hop where
a junction attaches is scaled by

    r = [ρ_contig + ρ_J] / ρ_contig · r_true   =   (1 + ρ_J/ρ_contig) · r_true

**Measured: median 11.5× off capture, 2.9× under it**, and the delivered `φ = c_g·E_g/M` tracks the
inflation linearly across five decades (1.04 → 0.93; 50.9 → 23.3). The intron's 97 %-gDNA composition is
being applied to a total that includes mature RNA the intron cannot contain.

⭐ **The general rule this implies, and it explains the graft and the peel rather than sitting beside
them:** the reframe carries only components **structurally admissible at both endpoints**. A component
admissible at one end and not the other must cross by a **measurement**, not by the imputation — and
that is exactly what the two routing operators already are, one per direction:

| direction | what crosses the admissibility boundary | operator |
|---|---|---|
| line → EXON | mature RNA is measured at the line (`sj_count`) and belongs to the exon | the **graft** |
| EXON → line | the exon's mature RNA cannot cross the seam contiguously | the **peel** |
| intron → line | *nothing* — yet the line's `ρ_J` is in the ratio anyway | ⛔ **unhandled** |

So the third row is the gap, and it is the one the doc's own intron → intron-exon-edge worked example
walks through. ⚠ It is stated here as a derivation result and **not implemented**, which is why it is not
in `EQUATIONS.md`.

## 2. ⭐⭐ …but a COMMON-MODE scale error is second order, and that is why fixing it buys ~0

`λ = log(f_g/f_R)` is invariant to a shared rescaling of both arms, and the inflation above *is* shared.
Measured consequence: **a 50× inflation moves `f_g` by 0.12** (0.768 → 0.651). Removing the mismatch
entirely (`ρ_tot` over the contiguous set only, one thing varied, 4 conditions) moves the broken class by
at most **−6.8 %** and is **worse on two of four conditions**.

⭐ **The design lesson, and it should shape where effort goes:** getting `r` *wrong by a scale* is cheap;
what costs is an error that **breaks the common-mode symmetry** — i.e. one that hits the two arms
differently. The graft does exactly that (it adds an absolute measured density to one arm only), which is
why arm C below moves `node/exon` by 25 % while the frame fix moves nothing.

## 3. ⭐ The junction channel is LOAD-BEARING — measured, not argued

Deleting it outright (`junction_count = 0`: frame, graft and peel together) is **worse overall**:
`node/exon` mwae 0.285 → **0.611** and library `f_gdna` 0.774 → **0.841** against a truth of 0.750 on
the unstranded capture-off row. ⚠ But under capture it *improves* `node/exon` 0.0431 → **0.0317** and
`exon|exon` 0.0461 → 0.0424. **So the graft is right off capture and over-states under it** — a
capture-dependent frame problem, matching the solver's own M8/P1d note (median φ = 2.45 under capture).
That is the junction channel's real defect, and it is at the **exon**, not at the intron|exon line.

## 4. ⛔ The panel cannot answer the intron↔exon question, and this is structural

Every ladder condition is `nrna_none`. So `edge/intron|exon` truth is exactly 1.0 **because there is no
nascent RNA to cross it**, not because RNA cannot. `TRAPS.md` F9 says only that *mature* RNA does not
cross such a seam contiguously. ⭐ The doc's plan — "start with single-stranded transcripts, capture OFF
and ON, derive the arithmetic" — needs a **nascent-bearing** condition to have a non-degenerate target on
that seam at all. It is `ROADMAP.md` §2 step 3.

## 5. On the toy harness

Injection of the *library-level* scalars already exists and is the falsification handle two instruments
use (`InjectedCalibrationPriors`; `composition_evidence_census.py --inject-kappa` reproduces a natural
experiment exactly). What does **not** exist is the thing the doc asks for: adding spike-in transcripts to
a cached scenario and re-simulating only those under the cached scenario's locked parameters. Nothing was
built for it this session.




---

# ✅ THAT DERIVATION IS DONE — 2026-08-04

The section that used to be here asked for **one gDNA scale rule for capture-OFF and capture-ON**, with
capture-OFF as the `r = 1` limit of the general rule and no branch on capture. It is answered, and the
answer is not a better scale factor:

    rho_c(src) · rho_tot(dst)/rho_tot(src)  ==  phi_c(src) · rho_tot(dst)

The reframe is a **pure composition imputation** — the source's density SHARE applied to the destination's
observed total — with no level transport in it at all. So it cannot be repaired by a correction factor: the
missing factor is `phi_c(dst)`, the destination's own belief, and multiplying that back in *is* `TRAPS.md`
D4. ⭐ **The fix is a LICENCE.** The reframe is allowed exactly where a composition imputation is allowed —
the source must have SUPPLIED both components of the pair, which is the λ-emission gate's own predicate,
already derived and previously applied to the τ stream only — and where it is not, the gDNA **level crosses
unscaled**, because gDNA is uniform along the genome before capture.

⭐⭐ **Capture is not a case and needs no branch:** the capture landscape is carried by the
structurally-pure-gDNA population's OWN measurements, which the relay's mass pin restores at every such
object, so the level an exon receives is its flanking EDGE's *enriched* measurement rather than the
off-target floor. A per-capture-class landscape ratio was built to do this explicitly and measured
**byte-identical off capture**; it is deleted, and `node_geometry` records why.

| where it now lives | |
|---|---|
| the derivation | ⭐ **`EQUATIONS.md` §3.5** |
| the lesson, and the two dead ends | ⭐ **`TRAPS.md` D4c** |
| the gates | `tests/calibration/test_gdna_scale_rule.py` — 6, perturbation-verified |
| the numbers, and what is next | `ROADMAP.md` §2 · `ISSUES.md` C1 (closed) / C11 / C6 (restated) |

⚠ **Two residuals, both separate mechanisms, both measured, neither this rule's** — see `ISSUES.md`:
**C11**, the relay's mass pin as a random walk on the gDNA level; and the EDGE being a *partial*-overlap
crossing while an exon interior is fully covered, which makes the delivered level a lower bound under
capture and needs a probe-geometry model the tool does not have.

⛔ **And the trap that was avoided, because it would have scored well:** the affine-in-overlap
extrapolation `e_g[exon] = 2·e_g[crossing] − e_g[off]` is *exactly* the simulator's own retention law
(`sim/capture/sampler.py`: `off_target_weight + binding_per_base · overlap`), so fitting it would be
scoring against the substrate that generated it. Real hybridisation is sigmoidal in overlap, and by
Jensen a convex law makes the affine form an under-statement — consistent with the measured
`2·224 − 1 = 446` against a true `596`.
