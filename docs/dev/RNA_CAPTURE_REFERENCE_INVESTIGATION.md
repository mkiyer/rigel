# The capture reference without a mode detector: an investigation

*2026-10-10, after the owner's ruling that an on/off capture detector cannot ship: either the reference the
shipped reader needs is approximated without a mode, or the honest reader is made fast enough to ship.
Sandbox note; nothing here is authoritative. Nothing in `src/` changed. Every number comes from prototype arms
run through the production reader and consumer at production cost, scored by the existing truth instruments
and the pinned end-to-end harness of the capture-reader benchmark (`RNA_CAPTURE_EXTERNAL_REVIEW.md`, whose
archive this note extends). Cells marked [pending] are runs in flight when this text was written.*

## 1. What the reference is for, and what the detector was really doing

The shipped reader reads each object's deconvolved gDNA count on its own opportunity as a Poisson
observation under the fitted gDNA landscape and publishes `E[min(ρ/ρ_ref, 1)]`. The EM is invariant to a
common factor across a locus, so `ρ_ref` matters only through two things: the **clip** (objects denser than
the reference saturate at 1) and the **`None` path** (no located mode → every efficiency exactly 1). An
always-defined reference therefore changes nothing on the captured objects it sits among; what it changes
is which objects saturate, and whether a library with no mode contracts at all.

The lineage, verified in the tagged source:

| version | the reference | fails where |
|---|---|---|
| 0.7.1 | the rightmost significant peak of a **mass-weighted** kernel density of per-region gDNA density (bandwidth 0.4, prominence 0.05, snapped to a node); `None` below 5 nodes | outliers: one exon holding 8 % of LBX0190's gDNA mass becomes the reference (10,801× the median) |
| shipped | the landscape's located enriched mode: the basin above the depleted one with the most located kernels, a mode only if it has more than √n members with median width ≤ 1 nat² | sparse panels: LBX0190 and MO_3021 have no such basin, so they contract nothing; a single probe never forms a mode |

The second line is the owner's objection. The first line is why the shipped rule was made stricter: the
2026-09-15 finding that anchors' resolution walls were being counted as located members on sparse plasma
libraries (`ISSUES`, the ruler-reference finding), after which the census was made to fail safe.

## 2. Reference estimators that need no detector, scored against truth

Computed from the production inputs (`count_gdna_region / gdna_region_eff_len`) on staged conditions, as a
multiple of the simulator's own fully captured density (1.0 is right) or, on capture-off rows, of the
background (`dump_inputs.py`, `inputs/`):

| row | located mode | mass-weighted median | mass-weighted mean | 90 % of mass below | 0.7.1 kernel peak | landscape top 1 % mass |
|---|---|---|---|---|---|---|
| benign g05 / g25 / g50 / g98 captured | 1.00 / 1.05 / 1.02 / 1.04 | 1.00 / 1.03 / 1.02 / 1.02 | 0.94 / 0.96 / 0.96 / 0.96 | 1.12 / 1.10 / 1.07 / 1.06 | 1.02 / 1.04 / 1.03 / 1.03 | 1.40 / 1.28 / 1.24 / 1.19 |
| strength 100 / depth 1/10 / suite RNA-short | 1.08 / 1.15 / 1.06 | 1.09 / 1.15 / 1.08 | 1.06 / 1.09 / 1.02 | 1.15 / 1.43 / 1.17 | 1.08 / 1.14 / 1.08 | 1.32 / 2.06 / 2.48 |
| weak capture (strength 0.1) | 0.99 | **0.08** | 0.30 | 1.02 | 1.00 | 1.28 |
| sparse layout / single probe | 0.44 / none | 0.44 / 0.00 | 0.45 / 0.02 | 0.86 / 0.00 | 0.45 / 0.00 | 1.16 / 0.01 |
| capture-off g05 / g50 / g98 (× background) | none | 0.93 / 1.00 / 1.00 | 1.10 / 1.00 / 1.00 | 1.67 / 1.04 / 1.03 | 0.97 / 1.00 / 1.00 | 2.48 / 1.65 / 1.59 |

On the real libraries, as a multiple of the median region density:

| library | gDNA mass; top 1 % of regions hold | located mode | mass-weighted median | mass-weighted mean | 0.7.1 kernel peak |
|---|---|---|---|---|---|
| LBX0190 | 4,855 fragments; 34.5 % (one exon alone 7.7 %) | none | 18.5× | 1,722× | 10,801× |
| MO_3021 | 82,381; 8.6 % | none | 1.0× | 149× | 0.9× |
| LBX0588 | 741,693; 3.6 % | 2,699× (12,551 members) | 1,294× | 1,590× | 2,352× |
| VCaP | 1,008,422; 0.4 % | 98× (24,932 members) | 64× | 90× | 107× |

Reading: the **mass-weighted median** reproduces the located mode within 10 % on every synthetic row that
has one, needs no constant, always exists, and is the one statistic the outlier exon of LBX0190 cannot drag.
It fails only where the located mode also fails (the sparse layout, a single probe) and where the captured
objects hold less than half of the gDNA mass (weak capture, 0.08×; VCaP is already at 0.65× of its mode).
The 0.7.1 kernel peak handles weak capture but needs two constants and follows the outlier.

## 3. The decisive test: a reference that always exists, end to end

Production reader, production consumer, production cost; only `ρ_ref` changed (`mmed`: the mass-weighted
median; `kde071`: 0.7.1's peak). Transcripts % / genes %, the pinned fractional harness:

| row | shipped | mass-weighted median | 0.7.1 peak |
|---|---|---|---|
| g00 ss0.50 capture-off | 5.67 / 1.01 | 12.08 / 1.35 | 12.08 / 1.35 |
| g00 ss0.70 capture-off | 6.35 / 0.90 | 13.32 / 0.98 | 13.32 / 0.98 |
| g00 ss0.99 capture-off | 6.99 / 0.83 | 12.67 / 0.88 | 12.67 / 0.88 |
| g00 ss0.99 captured | 7.35 / 0.67 | 7.77 / 0.67 | 13.78 / 0.66 |
| g05 ss0.50 capture-off | 7.91 / 1.09 | 8.01 / 1.04 | 7.99 / 1.04 |
| the remaining benign rows, single probe, strength 0.1 and 100, both probe layouts, RNA-short arm, depth 1/10, the ladder's six stranded rows, the four real libraries | [pending] | | |

**A zero-gDNA capture-off library, which is ordinary RNA-seq, loses six transcript points under any
always-defined reference.** The mechanism is not the reference: on those rows the only gDNA is the
calibration's own false gDNA, a few hundred fragments the solve attributes to RNA-rich exons
(`ISSUES: the-gdna-prior-enters-psi-twice` measures 366 at g00), and a reference placed anywhere inside that
distribution turns it into a contrast between a gene's RNA-rich and RNA-poor exons. The census's `None` was
not only "no mode found": with fewer than √n located members it refuses to read false gDNA as capture, and
that refusal is load-bearing on exactly the most common library type.

## 4. Reading the production inputs without a reference at all

If the reference matters only through the clip, the clip can go: read each object's posterior MODE under the
landscape from the same Poisson observation (`ship_map`, with the touched junction cap), or from the solver's
own belief with its variance (`ship_lnmap`). No reference, no clip, no `None`, production cost.

| row | shipped | posterior mode of the inferred count | belief-aware variant |
|---|---|---|---|
| g00 ss0.50 / 0.70 / 0.99 capture-off | 5.67 / 1.01; 6.35 / 0.90; 6.99 / 0.83 | 6.68 / 1.88; 6.75 / 1.31; 7.71 / 0.83 | 6.47 / 1.88; 6.79 / 1.31; 7.68 / 0.83 |
| g00 ss0.99 captured | 7.35 / 0.67 | 14.65 / 0.67 | 13.36 / 0.67 |
| g05 ss0.70 / 0.99 captured | 6.56 / 0.77; 4.93 / 0.35 | 6.09 / 0.77; 6.39 / 0.37 | 5.45 / 0.70; 5.20 / 0.35 |
| single probe, g50 rows (genes) | 3.56; 3.87; 2.55 | 2.66; 1.54; 1.43 | 2.62; 1.41; 1.42 |

Across the test chromosome's panels the mode readout of the inferred count tracks the honest candidate
(stratum means, transcripts / genes, shipped → readout): benign stranded capture-off 14.86 / 3.96 →
15.25 / 4.17, stranded captured 22.02 / 10.05 → 22.74 / 10.22, unstranded capture-off 16.15 / 4.63 →
16.16 / 4.56; strength 100 10.39 / 3.31 → 11.35 / 3.51; RNA-shorter gap 24.54 / 11.99 → 25.86 / 12.53;
depth 1/10 19.11 / 1.35 → 19.82 / 1.43; junction-probed 39.20 / 7.53 → 31.86 / 9.13; sparse-probed
36.17 / 12.26 → 34.63 / 12.68; single probe 9.22 / 2.12 → 8.95 / 1.22. Its failures are the zero-gDNA rows
alone: captured 7.83 / 1.42 → 11.51 / 1.42 (one row, g00 ss0.99, carries it: 7.35 → 14.65) and capture-off
6.67 / 0.87 → 7.23 / 1.07.

On the real libraries these readouts cost 3–6 s of reader time per whole genome at 6 GB peak (against the
honest prototype's 4.8 h and 37 GB) and move the EM's gDNA fraction the way the honest median did: LBX0190
0.0850 → 0.0971, MO_3021 0.1587 → 0.1638, LBX0588 0.9265 → 0.9260, VCaP 0.2395 → 0.2400 against a truth of
0.2518 (neither reader reaches it; the cheap readout moves it by 0.0005).

The same false gDNA breaks these too (g00 captured 7.35 → 14.65; capture-off genes 1.01 → 1.88), because
the inferred count is a point estimate made under the same prior, and no summary of its posterior can know
that the RNA could have explained those reads. The honest candidate of the benchmark reads the raw strand
columns with the RNA amount integrated out and is immune: its zero controls read 7.95 / 1.42 → 7.87 / 1.40
captured and 6.67 / 0.87 → 6.76 / 0.87 capture-off.

## 5. The structural variants: a reference from located gDNA only, or from gDNA-majority objects only

`mmedloc`, the mass-weighted median over regions whose gDNA the solve can locate (`Var(log f_g) ≤ 1 nat²`,
the production floor that decides what trains the landscape), is **refuted**: the false gDNA at RNA-rich
exons is located, because those exons have thousands of reads and the solve is confident in its small,
wrong gDNA share. The zero-gDNA capture-off rows read 12.08 / 13.32 / 12.67 against 5.67 / 6.35 / 6.99,
identical to the unrestricted median, and the g05 capture-off rows cost +2.0 and +2.4 points (8.27 against
6.24; 8.77 against 6.42) where the unrestricted median costs [see §3]. The captured rows are unchanged
(g05 ss0.99 4.93 → 4.93; g50 10.36 → 10.27) and the single probe is not read (7.85 → 9.63 at g05 ss0.99).

**Weak capture is where a mass-based reference breaks, and it does not break gracefully.** On the
strength-0.1 panel (10× enrichment, the captured objects holding well under half of the gDNA mass) the
mass-weighted median sits at the background (0.08× the captured level; the located mode reads 0.99×), and
at 50 % gDNA the rows go 17.24 / 5.14 → 23.78 / 17.64, 11.47 / 1.48 → 25.02 / 16.82 and 9.23 / 1.01 →
15.89 / 9.77. A reference below the enriched level saturates every captured object at 1 while the
uncaptured ones read their ratio to a background that is now the unit, and `priors.assemble_priors` prices
the gDNA component's own length through the same efficiencies, so the locus's gDNA-against-RNA split is
wrong by the missing factor and the gene counts follow. The census's membership test is what kept the
shipped reader at 0.99× there.

`mmedmaj`, the mass-weighted median over objects at which gDNA is the **majority** component of the
deconvolved mass (`f_g > ½`: the object's density is a gDNA measurement and not a leak from its RNA), is a
different statement in principle: on a zero-gDNA library no object should be gDNA-majority, so no reference
would exist and every weight would be 1, a structural `None` about the library's gDNA and not about a mode.
**Refuted in practice**: on every zero-gDNA capture-off row exactly one object is gDNA-majority (a tiny
object carrying only false gDNA), one object is enough to define a median, and the rows read 12.08 / 13.32 /
12.67 as before; a minimum membership would be the census's own threshold again. Where real gDNA exists it
behaves like the unrestricted median (g05 capture-off +2.0…+2.2; g25–g50 within ±0.15; the single probe not
read; weak capture the same catastrophe, genes 1.48 → 16.82).

## 6. What it would cost to make the honest reader fast

The NumPy prototype's time is the both-strand tilt integral: on the suite's genome-scale rows 14 % of slots
are both-strand and take 93 % of the time (median 0.1 s each, 0.9 s at the 90th percentile); single-strand
slots take 2 ms. The work per slot under the prototype's own lattice rules with the window rule is bounded
(about 17 nodes per axis; tilt nodes `max(24, ⌈(π/2)·√(n+1)⌉)`), so the evaluation count per genome is a
derivation, not a guess. Counted on the ladder's genome-scale slot population (70,176 slots, 14 % both-strand,
the deepest both-strand slot 35,184 reads needing 295 tilt nodes) and scaled to the plasma library's
2,087,476 slots:

| lattices | evaluations per genome | at 20 ns each, one thread | on 8 cores |
|---|---|---|---|
| windowed (ψ's own rule: nodes only within the integrand's support) | 2.9 × 10⁹ | 58 s | about 8 s |
| the prototype's full lattices (260 density nodes × up to 150 RNA nodes × tilt) | 3.9 × 10¹¹ | 7,800 s | about 16 min |

Both-strand slots are 82 % of the evaluations either way. The NumPy prototype uses the full lattices and
pays Python per slot, which is how it reached 4.8 hours; the window rule alone is a 135× reduction and the
native loop removes the per-slot overhead. A native kernel streaming with the blocks, bit-equal to the
prototype's frozen objects, lands inside the tool's three-minute budget with room to spare; the prototype is
the reference it is gated against.

## 7. Conclusion

1. **The reference can be approximated without a mode.** The mass-weighted median of the gDNA density
   reproduces the located mode within 10 % wherever the captured objects hold at least half of the gDNA mass
   (every normal-strength row, LBX0588, nearly VCaP), needs no constant, always exists, and the outlier exon
   that captures 0.7.1's kernel peak cannot move it. It gives LBX0190 a reference (18.5×; 26.9× over
   gDNA-majority objects) where the census gives none.
2. **An always-defined reference is not a shipping option on its own**, for three measured reasons, each
   of which the census's thresholds were silently handling: a zero-gDNA capture-off library, which is ordinary
   RNA-seq, loses six transcript points because the calibration's few hundred false gDNA fragments become a
   contrast (§3); the g05 capture-off rows lose two points to Poisson noise at low gDNA; and weak capture
   (captured mass share below half) breaks the gDNA component's pricing at high gDNA (genes 1.5 → 16.8, §5).
   Restricting the reference to located gDNA changes none of this, and restricting it to gDNA-majority
   objects changes none of it either: a single false-gDNA object defines the median on a zero-gDNA library.
3. **Reading the production inputs without any reference** (the posterior mode of the inferred count)
   matches the honest candidate on every panel at production cost (3–6 s per genome), and fails only on the
   zero-gDNA rows, again through the false gDNA (§4).
4. **The honest likelihood is the only reader immune to the false gDNA**, and its cost is a derivation,
   not a hope: 2.9 × 10⁹ evaluations per genome under the prototype's own lattice rules with the window
   rule, about 8 s on 8 cores in a native kernel (§6). The prototype supplies the bit-level references its
   gates need. This is the owner's option B, and it is the shipping path: build the native kernel with the
   speed gate first, keep the census as the fallback during the build, and the OFF-row spread of the
   benchmark stays the one declared cost.
5. **The defect underneath every cheap route is the false gDNA** (`ISSUES: the-gdna-prior-enters-psi-twice`,
   366 fragments at g00). Its measured fix cuts it to 13 but loses on every captured stratum, so it is parked
   on an open ruling. If that ruling lands, the mass-weighted median becomes a legitimate replacement for the
   census and the cheap mode readout a legitimate reader; until then neither is.

## Receipts

Scripts and results in `.cache/rigel_runs/2026-10-10_capture_reader_benchmark/` (`ship_reader.py`,
`dump_inputs.py`, `inputs/`, `results_ship/`, `results_ref/`, `e2e_ship*/`, `e2e_ref/`, `e2e_loc/`,
`real/ship_*`, `real/ref_*`, `real/loc_*`, `sweep_ship.log`, `sweep_ref.log`).
