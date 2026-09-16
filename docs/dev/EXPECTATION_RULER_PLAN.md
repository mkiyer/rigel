# ⛔ EXECUTED 2026-09-16 (session 4). The rulings are `DESIGN.md` §7.2 and `EQUATIONS.md` §11, the closures are in `ISSUES.md`, the instrument is `scripts/design/ruler_vs_truth.py`; nothing below is current. Kept as the record of the plan the session ran.

# The expectation ruler — derivation and implementation plan (working note, 2026-09-15)

The plan for the ruler's repair: what is settled by measurement, what remains to derive and the decision each
derivation needs, the commit sequence, and the gates. The problem statement and the measurements are the
companion note on effective length under capture; the issues are `ISSUES: ruler-multimapper-floor-caps-the-
correction`, `ISSUES: the-ruler-reference-on-sparse-real-libraries`, `ISSUES: ruler-witness-geometry-on-
transcript-panels`. Nothing here is a ruling until the owner takes it into `DESIGN.md`.

## 0. The plan on one page

**What the ruler produces.** For the EM, one number per transcript: its effective length under capture. That
is the spliced effective length the EM already uses — how many places an RNA fragment can start on the
transcript, junction-spanning starts included, from the RNA fragment-length distribution — times a shrinkage
factor between 0 and 1: the transcript's average capture efficiency relative to a fully probed transcript. The
locus's gDNA component gets the same treatment over its genomic span. The ruler decides nothing about which
fragments are RNA or DNA; the calibration did that before the ruler runs. A transcript competes with other
transcripts on its spliced length and with gDNA on the unspliced fragments; the EM already knows which fragments
are compatible with which template, and the ruler only says how much of each template capture sampled.

**Two different sets, and the confusion between them was the problem.** The LENGTH is counted on the spliced
transcript: every RNA start position, whether or not the fragment spans a junction. The EVIDENCE is unspliced
gDNA over the same sequence: gDNA is uniform before capture, so its density after capture on a piece of the
genome is that piece's efficiency. A transcript of tiny exons has almost no place for an unspliced fragment to
sit INSIDE an exon, and that is a statement about evidence, not about length — its length is its spliced
bases as for any transcript. Its gDNA evidence sits at its exon edges: the accumulator counts a fragment at
every boundary it crosses, so a gDNA fragment overlapping a tiny exon is counted at both edges.

**The three parts, in the order they are built.**

1. **The unit.** `ρ_ref`, the fully captured gDNA level, is the located mode of the fitted gDNA landscape; the
   membership test is repaired so anchors' walls cannot vote, and `None` is declared when there is no located
   mode. Small, isolated, gated on the real libraries. Commit one.
2. **The estimate.** Each piece's efficiency `min(ρ/ρ_ref, 1)` is the posterior mean under the landscape prior
   given the gDNA evidence overlapping the piece: first its own contained count (prototyped and measured: the
   30× floor bias gone), then its crossings apportioned by the geometry, which is what a tiny exon needs. The
   floor is deleted. Commit three, in two steps.
3. **The sum.** A transcript's factor is the length-weighted mean of its pieces' efficiencies over its own
   spliced bases, with the fragment-length end taper; `eff = fl × factor`. No boundary or junction object in the
   length. The locus gDNA component is the same sum over its genomic pieces. Commit three.

Between one and three: the truth instrument that scores every step against the simulator's own capture-aware
effective length becomes a design instrument (commit two), and the test chromosome gains a tiny-exon block,
probed and unprobed, so step 2's second half is judged on the structure it exists for.

**A worked locus.** A single 3 kb exon, probed: fl ≈ 2.8 kb; evidence is its own contained gDNA, dense;
efficiency ≈ 1; factor 1. A five-exon transcript with 500 bp exons, the first two probed: fl ≈ 2.3 kb; each
exon's own count gives its efficiency, ≈ 1 or ≈ 0.001; factor ≈ 0.4, the probed share of its bases. A
ten-exon transcript with 40 bp exons, probed across junctions: fl ≈ 200 bp; no exon holds a fragment; each
exon's evidence is its two edge crossings read against the intron's level; factor ≈ 1 if the edges carry
captured gDNA, ≈ 0.001 if not; what the junction probes add for cDNA over gDNA is invisible and declared.
Today the third transcript reads 0, then the floor.

**Decided.** The gDNA landscape is the witness; total density is refuted. The floor goes. The junction objects
go. The reference's members must have a location. **The open questions** (§5) are settled the owner's way
(2026-09-16): try each way on the truth instrument and show that it works; the numbers decide.

The rest of this note is the detail behind each part.

## 1. What measurement has settled

* The ruler's geometry is right: on the test chromosome, probed transcripts read within ±0.1 nat of the
  simulator's own capture-aware effective length for 99–100 % under every gDNA ruler.
* The floor `w = C/(C+1)` is a +3.4 to +3.7 nat bias on unprobed transcripts at every gDNA level; without it the
  plug-in is within 10 %, the population posterior within 5 %.
* The population posterior also lifts the low-gDNA end of the operating curve (91 % against 65 % of probed
  transcripts within ±0.1 at 5,850 gDNA fragments in a library of 11,700).
* The junction objects must go with it: on a dense annotation the pieces beside a junction have no support,
  and forty junctions reading a guess swamp four hundred bases of measurement. The best of four rules measured,
  the transcript's pool, is the same as having no junction objects; the derivation that removes them, and puts
  a piece shorter than a fragment back into the length with its edge crossings as evidence, is §2.3.
* The reference's membership test admits anchors' walls on sparse real libraries; that repair is a precondition
  (§3, commit one).
* Total density as the ruler is refuted with numbers; the gDNA landscape stays the witness (owner, 2026-09-15).

## 2. The derivation, item by item, with the decision each needs

### 2.1 The per-object posterior — two candidates, decided by the truth instrument

**(A) The population posterior** (prototyped and measured). `c̃_o = E[min(ρ/ρ_ref, 1) | k_o, S_o]` on the
landscape's own grid: prior `P(log ρ)` from the last refit's `DensityLandscape`, likelihood `Poisson(k_o | ρ S_o)`
with `k_o` the deconvolved gDNA mass and `S_o` the object's support (`gdna_region_eff_len`,
`gdna_boundary_eff_len`). No constant: the grid and the prior exist, the likelihood is the counting rule, the clip
is the ruler's definition. Two blemishes: `k_o` is a mass, not a count, so the Poisson through the gamma function
is a continuation; and `k_o` already carries the same prior through the solve, so the prior enters twice in a
small way.

**(B) The solve's own posterior.** `RegionBelief` publishes, per chain slot (regions and boundaries alike), the
gDNA share `f_g` and `Var(log f_g)` over the λ lattice. With `log f_g` Normal at that mean and variance,
`ρ = f_g M/S` is log-Normal and `E[min(ρ/ρ_ref, 1)]` is a closed form (a truncated log-Normal mean, two Φ terms).
No grid, no continuation, no second prior, O(n). Two risks: a slot without a composition (`var = ∞`) has no
posterior of its own and needs the population anyway; and a wide skewed posterior has `E[ρ]` far above its
median, which could over-read an unprobed exon whose solve is wide.

**Decision.** Prototype (B) beside (A) in `ruler_variants.py` and score both on the truth instrument over both
panels, the ladder's two rows and the depth ladder. Ship the simpler one that is not worse; if (B) matches (A)
where a slot is solved and needs (A) where it is not, ship (A) alone — one rule beats two.

### 2.2 The clip

Keep `min(·, 1)` inside the expectation. Measured: probed transcripts within ±0.1 nat; the shortfall at the
reference is precision-dependent and small (0.98 at 40 fragments, 0.96 at 52, 1 in the limit) and its common
part cancels in the anchoring. Carry `E[ρ]/ρ_ref` clipped afterwards as an arm for the record; it is exact below
the reference and noisier above, and the instrument will say whether the difference is visible.

### 2.3 The length and the evidence are two different sets

**What the code does today.** A transcript's object set is the pieces its exons overlap plus the boundaries
strictly interior to an exon (a piece cut by another transcript's start, end or alternative splice site);
the exon|intron and exon|intergenic edges are excluded, and one junction object per adjacent exon pair is
imputed from the two pieces beside it. Each object supplies its own deconvolved gDNA count as its evidence,
and each object's WEIGHT in the transcript's mean is its gDNA-frame support: the contained support
`E_f[(ℓ − w + 1)+]` for a piece, `w − 1` for a crossing. Two things follow that the owner's question exposes.
A piece shorter than a fragment has zero contained support, so its bases carry no weight and it supplies no
evidence — the transcript of tiny exons reads 0 today, then the floor. And the crossing supports of a
transcript's interior boundaries and junctions each count `w − 1` starts, which tile the transcript's start
positions only while every piece is longer than a fragment; where pieces are short a start is counted at
every boundary its fragment crosses, and the boundaries take too much of the weight.

**The deposit rule, which decides what a boundary's count means.** The accumulator's specification
(`tests/native/_accumulator_reference.py`): regions count fragments that fit inside, boundaries count
fragments that cross, and EVERY crossed boundary receives the full weight — no partitioning, so a boundary
count is a density and not a share of a total. A gDNA fragment overlapping a tiny exon is counted at both
of its edges. Its evidence exists; the ruler has been excluding it.

**Role one: the length.** A transcript's effective length under capture is the sum over its own start
positions of the efficiency of the fragment starting there. With a fragment's efficiency the mean of the
per-base efficiency over the bases it covers, that sum is

    eff_t = Σ_{x ∈ t} c̃(x) · τ(x),      τ(x) = (starts whose fragment covers base x) / w,  averaged over f,

the transcript's own bases at their pieces' efficiencies, weighted by the fragment-length end taper, and
`Σ_x τ(x)` is the fragment-length marginal `fl_t` exactly. So

    factor_t = Σ_p ℓ_p^τ · c̃_p / Σ_p ℓ_p^τ,      eff_t = fl_t · factor_t,

with `ℓ_p^τ` the taper-weighted base count of piece `p`, computable in closed form from the pmf. No
boundary object, no junction object and no contained support enters the length: a junction-spanning start
is counted through the bases it covers, a tiny exon carries its bases' weight, and `span = fl` holds for
every structure by construction — the property the junction objects were built to secure. The locus gDNA
component in `assemble_priors` takes the same form over the locus's pieces, introns included, since gDNA's
template is the genome.

**Role two: the evidence.** A piece's efficiency `c̃_p` is informed by every unspliced gDNA object whose
fragments overlap the piece: its own contained count, and the crossing counts at its boundaries — interior
boundaries, whose fragments lie in the transcript on both sides, and the exon|intron and exon|intergenic
edges, whose fragments lie half in the piece and half in a neighbour that is not the transcript's. A
crossing count's expectation is the support times the mean efficiency of the fragments crossing there, and
those fragments' bases are shared among the pieces they cover in proportions the geometry fixes — the pmf
and the pieces' lengths on each side — so

    E[k_o] = ρ_ref · S_o · Σ_p a_op · c̃_p,      a_op = the expected share of object o's fragment bases in piece p,

a small linear system per locus, with the population landscape as each piece's prior. A contained object
has `a = 1` on its own piece. A crossing at an interior boundary of a long exon has its bases in two pieces
of equal efficiency and reads that efficiency directly, which is the current rule. A crossing at the edge of
a tiny exon has a share `ℓ/w` of its bases in the exon and the rest in the intron, whose efficiency is pinned
by its own long contained support, so the exon's efficiency is identified from the edge counts against the
intron's known level. The posterior mean of each `c̃_p` under this model is what the length sums.

**The worst case, worked.** A transcript of many exons each shorter than a fragment: no piece contains a
fragment, every gDNA fragment over an exon crosses one or both of its edges, and the annotation cuts no
interior boundary into an exon that short. The length is its bases, as above. The evidence is the edge
crossings: for a probed tiny exon they carry many captured fragments against the intron's depleted level and
the exon reads probed; for an unprobed one they carry the depleted level and the exon reads depleted. The
only information the gDNA cannot give is the probe designed across the junction, the declared limit. Today
this transcript reads 0, then the floor.

**Implementation order.** (1) The length in per-base form, and the junction objects, the flank imputation
and the borrowed crossing support deleted. (2) The evidence from each piece's own contained count, the
population posterior of §2.1 — the prototype as measured. (3) The crossings added with the exact
apportionment `a_op`, the neighbours' efficiencies taken from step 2, one pass. (4) Iterate to the joint
posterior only if the truth instrument on a tiny-exon structure says one pass is not enough. Each step is an
arm on the truth instrument; the test chromosome gains a designed block of tiny exons, probed and unprobed,
before step 3 is judged.

**What no rule can see.** A probe designed across a junction binds the cDNA fully and the gDNA at a
fraction (`gdna_split_penalty` 0.2 in the simulator, unknown in reality): the gDNA witness under-reads the
mRNA's capture there by a factor the panel's design sets — `ISSUES: ruler-witness-geometry-on-transcript-
panels`, declared, not repaired.

### 2.4 The floor's removal and multimapper blindness

The floor's stated purpose is the accumulator's blindness to multimapping fragments (`AF_MULTIMAPPER_DROP`):
a repetitive region counts fewer fragments than exist, and the support `E_f[(L − w + 1)+]` is pure geometry with
no mappability in it (per-region mappability was removed from the index in v0.5.0 because the calibration did
not consume it; the index can still take an alignable store, which is the path a mappable support would revive). Under the posterior a blind region
reads as depleted and its transcript contracts. The synthetic panels carry no multimappers, so nothing measured
so far sees this.

**Decision.** Not a floor: a floor hedges every object for the sake of the repetitive few and costs the unprobed
class 30×. The honest home is the opportunity model — a support that counts only mappable start positions —
which is its own item in the opportunity layer (open it as `ISSUES: multimapper-blind-support`, adjacent to
`ISSUES: capture-blind-gdna-divisor`). Before the ruler ships without the floor, measure on the four real
libraries the ruler's factors by repeat overlap of the transcript (RepeatMasker, or the index's own splice
blacklist as a proxy): if repeat-rich transcripts contract anomalously, the support item moves ahead of the
port; if not, it stays queued. Either way the finding is recorded with its number.

### 2.5 The reference and its members — a precondition

`ISSUES: the-ruler-reference-on-sparse-real-libraries`: a basin's members are the kernels with a location
(count ≥ 1 fragment), the knn widths are read among them, and the choice among basins above the depleted one is
re-derived (largest located mass or highest located population, decided on the panels, which must not move, and
the four real libraries, where LBX0190 and MO_3021 must read `None` or their probed level). The result gains
the number of located members behind the reference, and the report prints the regime.

### 2.6 Where the code goes

`calibrate` (layer 7) computes the per-object efficiencies from the last refit's landscape and the deconvolved
masses and publishes them on the result beside the reference: `gdna_capture_efficiency_region`,
`gdna_capture_efficiency_boundary`, `gdna_reference_density`, `gdna_reference_members`. `capture_eff_length`
(layer 2) becomes geometry only — the object sets, the junction windows, the sums — and reads nothing but the
result; `assemble_priors` reads the same efficiency arrays for the locus gDNA effective length. The floor is
deleted in both. The layering is untouched: the result is a value handed down, not an import up. The schema test
and `_check_axis_array` cover the new arrays; the golden outputs move and are read for magnitude first.

## 3. The commit sequence

**Commit one — the reference's members** (one mechanism). Falsification test first: a landscape whose only "members"
above the depleted basin are anchors' walls must read `None`, verified failing on the shipped code; the
perturbation (walls admitted again) fires it. Gates: both panels' standing numbers unchanged to the fragment
(`calibration_vs_oracle.py` ③, `policy_benchmark.py`); the four real libraries as test input (LBX0190 and
MO_3021 → `None` or the probed level; LBX0588 and VCaP unchanged); the identity references re-frozen with the
reason. Docs: `DESIGN.md` §7.2 amended, the issue CLOSED with its numbers.

**Commit two — the truth instrument promoted.** `scripts/design/ruler_vs_truth.py` from `s4/ruler_vs_sampler.py`: the
sampler's effective length as the truth, arms as config values or a `--module` prototype, per-class table, a
`--self-test` on a synthetic landscape with a known partition, an index row in `CLAUDE.md` and a `TESTING.md`
section; the depth-ladder configs under `scripts/sim/configs/` as a side panel (owner decision: three files or
one with a depth axis). The suite's collected count re-derived, never adjusted.

**Commit three — the expectation ruler on the per-base length** (one mechanism, one commit). Order inside it: DERIVE
§2.1–§2.3 into `EQUATIONS.md` §11 (rewritten) → the test chromosome gains a tiny-exon block, probed and unprobed
(`docs/TESTING.md` §0a's recipe; both panels re-cached) → PROTOTYPE the (A)/(B) arms, the per-base length and the
evidence steps of §2.3 in `ruler_variants.py` →
A/B on the truth instrument across both panels, the ladder's two rows and the depth ladder; record the table →
falsification tests first and verified failing on the shipped code: an unprobed exon with no gDNA fragment
contracts to the depleted level and not to `1/(C+1)`; a transcript of exons shorter than a fragment reads its edge
crossings and not 0, and its pieces carry their bases' weight → `src/`: the result's new fields, `calibrate`'s computation, the ruler
and the prior as sums, the floor deleted → the perturbations (the likelihood dropped; the junction objects restored
at the flanks; the floor restored) each fire their gate → the truth instrument's standing numbers into `DESIGN.md` §7 →
`calibration_vs_oracle.py` ③ re-recorded (P/O stays 1 by construction of the O arm) and `policy_benchmark.py`
identical → `quant_accuracy.py --arm base` on both panels, the thermometer, which must not worsen → the identity
references re-frozen with the reason → the four real libraries re-run: factor distributions before and after, the
regime printed.

**Commit four — the multimapper check** (read-only). The measurement of §2.4 on the four real libraries; an issue opened
or closed with its number.

Then the port's re-baseline.

## 4. The gates, as a checklist

* Falsification tests written first, verified failing on the shipped code, then each perturbation fires.
* The truth instrument on both panels: probed within ±0.1 nat ≥ 99 % on the test chromosome; unprobed median
  |error| ≤ 0.1 nat at 25 % gDNA and above; the operating curve re-recorded on the depth ladder.
* The capture-OFF contract: every factor exactly 1, bit-identical to the fragment-length marginal.
* The zero-gDNA rows: `None`, every factor 1, declared.
* `calibration_vs_oracle.py` per stratum and `policy_benchmark.py --by-class` on both panels: the composition
  numbers unchanged; the ruler column re-recorded.
* `quant_accuracy.py --arm base` on both panels: not worse than the standing thermometer.
* The suite re-derived; goldens read for magnitude before `--update-golden`.
* `rename_identity.py --check` against the three references after every stage; re-frozen only with the reason.
* The four real libraries as test input, never design input; the regime printed for each.

## 5. The open questions — the owner's ruling (2026-09-16): try each way, and show that it works

Each question is two arms on the truth instrument, run on both panels, the ladder's two rows, the depth
ladder and the tiny-exon block; the table decides, and the losing arm is recorded with its number.

1. **The evidence model's depth.** One pass — each piece's efficiency from its own count, then the crossings
   apportioned with the neighbours held at that — against the joint posterior iterated to convergence. The
   tiny-exon block is where they can differ.
2. **The posterior.** The population's (A: the landscape prior with the Poisson likelihood on the deconvolved
   mass) against the solve's own (B: the belief's `f_g` and `Var(log f_g)` in closed form). If they agree where
   a slot is solved, (A) alone ships, since it also covers the unsolved slots.
3. **The clip.** `E[min(ρ/ρ_ref, 1)]` against `min(E[ρ]/ρ_ref, 1)`.
4. **Multimapper blindness.** The ruler without the floor, measured on the four real libraries by repeat
   overlap; the number says whether a mappable support has to precede the port.
5. **The depth ladder's home and the regime's wording** are housekeeping; the implementer proposes, the owner
   takes it into `DESIGN.md`.
