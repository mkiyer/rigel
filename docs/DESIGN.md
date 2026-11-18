# DESIGN — what is built, and the rulings that constrain it

**Purpose.** Rigel as the source builds it, in the order a run executes, with the standing rulings that bind each
stage. Everything here is settled and is not re-litigated. §0 is the binding vocabulary, §0b the 0.8.0 scope. A
module's docstring is its manual; this document says what each part is for and which rules bind it. Not here:
derivations (`EQUATIONS.md`), lessons (`TRAPS.md`), open problems and refusals (GitHub issues), the ranked next steps
(`ROADMAP.md`), the panels and gates (`TESTING.md`), and history (git).

---

## 0. Vocabulary

**Axiom 0 — RNA is RNA.** There are three populations and no fourth: gDNA, RNA+ and RNA−. "Mature" and "nascent" are
not populations and not a degree of freedom: RNA inside an intron is RNA that has not spliced at that position. The one
distinction the solver draws is whether a fragment **is spliced** (certified RNA, nothing to deconvolve) or **is not**
(the whole deconvolution problem). The populations at a slot are a function of two bits of the annotation:

```
T(slot) = {gDNA} ∪ {RNA+ if statics.free_pos[slot]} ∪ {RNA− if statics.free_neg[slot]}        |T| ∈ {1, 2, 3}
```

A composition question that names "mature", "nascent" or "a third component" has broken the axiom; re-ask it as "what
is this channel's opportunity for RNA at this object?", a geometry the index supplies.

The table binds prose, comments, docstrings and identifiers alike:

| the term | what it means | not |
|---|---|---|
| **region** | a contiguous genomic interval of the partition, with a length in bp | — |
| **boundary** | a single genomic position between two adjacent regions, 0 bp wide | `line`, `seam`, `crossing` (as an object) |
| **edge** | a graph edge, a boundary or an sj (`edges_df`, `EDGE_KIND_CONTIGUOUS`, `EDGE_KIND_SJ`) | a synonym for boundary |
| **splice junction** (`sj`) | directed `src → dst` in GENOMIC order, `src < dst` always | a synonym for boundary |
| **low / high end** | an sj's two boundaries in genomic order: `_lo` / `_hi`, `boundary_left` / `boundary_right` | `donor`, `acceptor` |
| **slot** · **step** | one entry of the chain (region, boundary, region …) · one move along it | `hop` for a step |
| **G1 object** | a structurally pure-gDNA slot: no RNA strand admissible (an intergenic region, an `intergenic\|exon` boundary); `region_geometry.g1_locked` | `anchor` for an object |
| **density** = **abundance** | counts per base; **counts** are integers | the simulator's molar `abundance=` |
| **crossing fragment** | a fragment that spans a boundary | — |
| **switched off** | an A/B with one code path disabled | `ablated` |
| **annotated** / **nascent** | OUTPUT words only: RNA on a reference transcript, single-exon included / on a synthetic span | `mature`, `mrna`, `nrna` in an output |
| **splice-out** / **splice-in** | a message crossing a boundary in the direction mature RNA departs / arrives through the sj there | `graft` |

A G1 object's certainty is the annotation's — real intergenic space carries unannotated transcription — which is safe
for a pooled level and no licence to treat one slot as ground truth. `donor` / `acceptor` are strand-dependent 5′/3′
names over structures stored in genomic order, so on a minus-strand sj a "donor" holds the acceptor; they survive only
in the index's flag bits (`FLAG_DONOR_s` marks the genomic-low end of an `s`-strand intron), which every reader reads in
genomic terms. Splice-out and splice-in name a step's semantics, not its arithmetic: a step always uses the low slot's
high-flank total against the high slot's low-flank total, so a predicate on the SIDE is the right shape, never one on
the direction. `deconvolve` is a different verb: gDNA is deconvolved off, RNA is the residual
(`TRAPS: two-masks-one-name`).

---

## 0b. The 0.8.0 scope

The version on disk is `0.7.1` (`pyproject.toml`). 0.8.0 is a release of the TOOL: validation
sets its supported uses and limitations, its known failures are reported beside its gains, and a claim of robustness
across inputs stays unsupported while a known instability persists (issue #1).
No preview ships: 0.8.0 waits for every `release-blocker` issue (`ROADMAP.md`, Next).

| the question | the number | the instrument |
|---|---|---|
| is the TOOL right? — what the release ships on | the transcript table against per-transcript truth | `quant_accuracy.py` |
| is CALIBRATION right? — ranks a calibration mechanism | the `CalibrationResult` against an oracle calibration | `calibration_vs_oracle.py`, `prior_vs_oracle.py` |

Neither number stands in for the other: the table is downstream of calibration and the EM, so ranking a calibration
mechanism on it is refused. An end-to-end A/B pair runs with the scan pinned (`--set scan.total_threads=1`) and
`--set em.assignment_mode=fractional`, so it is exactly reproducible; an effect is judged by its size, with genes and
pools read beside transcripts.

Three strata are in scope — unstranded × capture-OFF, stranded × capture-OFF, stranded × capture-ON. Unstranded ×
capture-ON is **deferred**: reported on every benchmark, since a panel blind to the tool's worst cell cannot tell a win
from a re-labelling, open to robustness work, never a release requirement and never tuned to. It is blind by algebra:
at κ = ½ the gDNA fraction cancels from the strand mean. Every score is read per stratum, never pooled.

- **Robustness and generality over aggregate simulated accuracy**, fragment-length distributions especially. Roughly
  5 % DNA/RNA separation, 5 % transcript and under 2 % gene error would be very good — guides, never thresholds or
  constants. A simpler model at roughly the existing accuracy is preferred, and a mechanism that wins only on a
  panel's structure is not a mechanism.
- **Conservative allocation under uncertainty**: where the evidence is ambiguous, overestimating gDNA is preferred to
  underestimating it. It adds no floor or constant, and it is not a confident measurement of capture.
- **Fragment length is in scope.** A length mechanism carries information that vanishes smoothly as the gDNA/RNA gap
  closes, and passes the zero-gDNA controls. The ladder gives both origins equal lengths on purpose: the scorer reads a
  table per origin, so a gap lets the EM split origins on length alone and masks calibration's bugs
  (`TRAPS: a-length-gap-bypasses-calibration`). A length mechanism is judged on the fl-gap side panels, in both
  directions, beside the equal-length control, against the simulator's own truth table.
- **The nascent scope ruling.** Real RNA-seq holds nascent RNA at very low levels; Rigel models it simply (§2.2) so the
  tool is robust to it, and nascent-capture experiments are out of scope. Absence is the default, and it is the
  likelihood's, never the prior's (§9.3). A number at the panel's nascent share is a stress reading, never a design
  driver; state its unit, since a span covers a gene while a mature transcript is spliced. A background pooled over
  intergenic and intron slots is safe on real data, not on the panel.
- **The release holds for the annotated pool.** At 5 % gDNA and above, 0.8.0 does not ship while an in-scope row's
  annotated error exceeds 0.7.1's on the same reads (`SUCCESS.md`, bar 5; issue #76, issue #9). This bar has no
  stated-limitation fallback.
- **One message policy.** `transfer` ships on every stratum. A row where it loses to `silent` is a defect in the
  messages, repaired there, never a reason to switch policy by stratum (issue #26, issue #77).
- **The splice blacklist is retired** (issue #49, issue #63, issue #64). The recommended STAR settings report a
  genomic fragment's contiguous reading beside its spliced one, which is the protection the blacklist gave; what they
  do not cover is stated (issue #7). A blacklist brought back later is regenerated at those settings first.
- **Ruled, not yet built.** A multimapper's candidates are pruned within each place, never across places, so no place
  is orphaned and the EM decides between places (issue #62; §9.4). An efficiency of exactly 1 gives RNA's own length,
  and only the excess over 1 is carried (issue #65; §9.2).

---

## 1. The shape of the tool

Three stages, in `pipeline.run_pipeline`:

1. **Scan** (`scan_and_buffer`): one C++ htslib pass pairs mates (§3), resolves fragments against the index, trains the
   strand model, buffers every resolved fragment for the EM (`FragmentBuffer`) and deposits per-object tallies into the
   accumulator (`AccumulatorPayload`), holding whole a fragment with several surviving explanations.
2. **Second pass, then calibration**: `drain_side_buffer` draws one explanation per held fragment and re-deposits it
   (§4); `calibration.calibrate` deconvolves every region and boundary of the drained tally into gDNA, RNA+ and RNA− by
   message passing over the region chain (§5–§7).
3. **Quantify** (`quant_from_buffer`): the length tables are decided (§8), the buffer is scored and routed, and a
   per-locus EM assigns fragments to the locus's transcripts and one gDNA component; calibration enters each locus as
   two scalars, a gDNA count and a gDNA effective length, never per transcript (§9).

An optional pass re-reads the BAM and stamps every record (§10). The native modules (`native.py`) are `_bam_impl`
(`native/scan/`, `native/calibration/`), `_resolve_impl` (`native/resolve/`), `_scoring_impl`, `_em_impl` (`em_solver.cpp` with `em_*.h`
and `connected_components.h`) and `_solve_impl` (`solve_kernel.cpp` with `psi_kernel.h`, `transfer_kernel.h`,
`transfer_rows.h`, `transfer_store.h`). Each carries a `testing` submodule (`native/testing/`), its seams bound for
the tests. The resolver's overlap index is the vendored `cgranges`, compiled into `_resolve_impl` and `_bam_impl`.

**Where a number lives.** `config` holds what a run may set (`PipelineConfig`; `IndexConfig` for `rigel index`) and
what the code fixes (`CONSTANTS`, each field documented with why it has its value). Code reads a constant by name and
never restates one; an experiment varies one on a copy (`dataclasses.replace`). Outside `config` live only identifiers
and formats, named mathematical definitions, and values the kernels export (`native.transfer_rows.ARM_EPS`, and
to the gates `_solve_impl.testing`'s `JEFFREYS_REF` and `OWN_EVIDENCE_EPS`); `tests/test_constants.py` refuses any other bare module constant.

---

## 2. The index

`rigel index` writes format 11 (`INDEX_FORMAT_VERSION`; a load refuses any other): the splice graph
(`regions.feather`, `edges.feather`; `calibration/splice_graph.py`), the transcript table and the unobservable mask, with `manifest.json` recording the sources. The build is deterministic.

> The genome is a graph. Regions are intervals; boundaries connect them. A fragment is a path. A region counts the
> fragments contained in it; a boundary counts the fragments crossing it.

- **Regions** tile each reference, cut at every exon endpoint of every annotated transcript, never merged; a 1 bp
  region is legal and common. **Contiguous boundaries** sit between adjacent regions, carry gDNA and RNA, and hold 8
  structural flag bits, TSS / TES / DONOR / ACCEPTOR × {+, −}, not mutually exclusive, which consumers read raw.
  **Junction boundaries** run `src → dst` with an explicit strand and carry spliced RNA only.
- Genomic order is a topological order, so nothing traverses a graph. Every sj closes an undirected loop, so an sj is
  a factor on its two endpoint regions, never a message channel (`TRAPS: splicing-makes-the-graph-cyclic`).
- `validate_graph` asserts invariants I1–I13 (no I10), two of them re-derived by a different algorithm
  (`TRAPS: self-checking-validator`). A region or boundary census belongs to an annotation, never to the tool. `reach`
  is in neither index hash, so a calibration-facing cache keys it (`scan_cache.reach_digest`;
  `TRAPS: a-hash-that-misses-its-artifact`).

### 2.1 The observable genome

> A base on which no fragment can be observed is no opportunity. The index records such bases as a mask
> (`unobservable.feather`, `rigel.observable`): with an `alignable` store, its zero set at the longest read-length bin;
> without one, the FASTA's N runs. Every opportunity is read over the observable genome.

`region_arrays` cuts every region at the mask into observable pieces, and every opportunity form has an observable twin
equal to it to the bit on an empty mask. The geometry, the length laws, the landscape, the capture efficiencies, the
second pass, the published lengths and the priors all read them; a template's fl-marginal length sums its unmasked
runs (`capture_eff_length.observable_fl_eff_lengths`). Counting unobservable bases reads them as the cleanest gDNA there
is and under-calls gDNA. `summary.json`'s `observable_genome` block records the mask.

### 2.2 Nascent spans

Multi-exon transcripts whose start and end sites lie within `IndexConfig.nascent_tolerance` bp share one span, their
outer envelope (`index.create_nrna_transcripts`). A span an annotated single-exon transcript covers takes it as its
entity; every other is one synthetic transcript, a gene-neutral row of `index.t_df` flagged `is_synthetic`, named by
its contributors' `span_t_index`. `is_single_exon` is the structural fact (true of every synthetic span), never
"manufactured"; the annotated filter is `~is_synthetic` alone (`TRAPS: nrna-does-not-mean-synthetic`).

---

## 3. The scan

### 3.1 The input

A paired-end BAM, name-sorted or collated, from STAR at the recommended settings (`MANUAL.md`: `--sjdbScore 0
--outFilterMultimapScoreRange 2`, which report a spliced and a contiguous reading where they score alike).

> **Rigel reads what the aligner reported.** A hit is a reported alignment, checked: the scan reads, validates,
> withholds and counts, and never chooses between a template's hits; the EM does.

Rigel moves no base, forms no alignment the aligner did not report and rejects no junction it did — a scope kept on
the tool's main defects, with genome information inside Rigel open to discussion on evidence. Both readings of a
genomic fragment the aligner splices reach calibration and the EM, which weigh them on the whole likelihood; one whose
contiguous reading the aligner does not report is read as spliced RNA (issue #7). A read pair's multiplicity is its
usable hits; NH only flags a truncated output.

### 3.2 A template's hits

The scan pairs mates itself, for any aligner, by geometry alone; STAR's HI is not read (`scan/template_hits.cpp`'s
`template_hits`, `hit_verdict`; `tests/test_hit_verdict.py`).

- **Locations and candidates.** The records of one alignment (mate, reference, strand, leftmost position, CIGAR with
  M/=/X as one operation) are one location; every R1 location with every R2 location on one reference with opposite
  strands is a candidate hit, and a mate with no location leaves one-mate hits.
- **The pair rule** (`FragmentResolver::validate_pair`). The 5′ ends, soft clips walked back, face each other. The
  molecule spans both mates (a read trimmed at its 5′ end starts inside its molecule), and its length less every skip
  either mate observed lies in 1..`max_frag_length`. Failing that, one annotated transcript overlapping both mates and
  holding every observed junction may witness a shorter path by cutting its introns in the gap (an implicit splice), if
  none holds more than K aligned bases (`intron_overhang_tolerance`).
- **The annotated-path rule** (`FragmentResolver::pair_holders`). A hit that observes a junction needs one annotated
  transcript holding every observed junction and overlapping both mates; overlap, not containment, since the EM prices
  the overhang. Rigel assembles no transcripts; a refused hit is an unannotated splice.
- **Off the index**: a hit on a reference the index does not hold has no partition, transcript or length there.
  **Chimeric**: Rigel does not quantify chimeras; a location that is one piece of an alignment in pieces (SA, a 5′ hard
  clip, an unnamed supplementary record) is unusable (issue #16).

A failing hit is dropped with its reason; a template left with none is withheld.

### 3.3 One outcome per read pair, one status per record

A read pair (records sharing a name, filtered ones included) has one outcome: **excluded** (every record unmapped,
QC-failed or a skipped duplicate); **withheld**, by reason `off_index` > `chimeric` > `unannotated_splice` > `unpaired`
(an unpaired pair also by the furthest pair check reached); **excluded by the multimap option**; or **placed**. Each
record is filtered, used, unused or unplaced. Outcomes sum to the read pairs and statuses to the records by
construction, and the scan stops on a difference (`tests/test_accounting.py`). The annotation pass re-judges every pair
with the scan's code and stops unless its counts and a digest of its pairs equal the scan's.

### 3.4 The accumulator

`tests/native/_accumulator_reference.py` is the executable specification; the C++ is gated on byte-identity to it, and
where it and a document disagree, the reference wins.

> **The deposit rule.** A region owns the fragments CONTAINED in it; a boundary owns the fragments that CROSS it; no
> object's mass is ever moved onto another.

A channel is stored where a named consumer reads it, and nowhere else:

```
region_contained    count[2]
region              start_count[2]                      the ledger: Σ start = deposited
boundary_unspliced  count[2]  inv_length_sum  mass
boundary_spliced    count[2]                  mass      certified RNA: nothing deconvolves it
sj                  count[2]  inv_length_sum  mass[2]   the strand split detects aligner artifacts
```

`count` is `+1` on every object a fragment touches; `inv_length_sum` (`Σ 1/A(w)`) cancels its own opportunity and is an
exact model-free density at a boundary or sj; `mass` is the conserved count, one unit per fragment shared over every
object it crosses, sj included. No spliced fragment is ever contained, and no length-sum bank exists
(`TRAPS: equal-lengths-carry-no-composition`).

- **The deposit weight is 1/opportunity**, never 1/length: model-free unconditionally at a boundary or sj, at a region
  only within `w ≤ ℓ`, so no region bank takes it (`TRAPS: a-cancellation-is-conditional-on-its-support`).
- **Numbers and strands**: a count is an integer and a fraction float64. Count banks reproduce exactly; fraction banks
  vary by about 1e-15 with the per-worker merge, and bit-identity there is not a goal
  (`TRAPS: integer-channels-reproduce`). Everything is stored by GENOME strand (`kNStrandColumns`, + then −); sense is
  derived, never stored; `align_strand` and `sj_strand` are independent.

### 3.5 Fragment length — one definition

`L` is the genomic span less every observed junction and every intron the hypothesis cuts in an unsequenced gap within
K of its ends. One function measures it everywhere (`path_length.h`'s `path_length`): the deposit, the second pass
(`Accumulator.length_under`) and every EM candidate on its own path. A covered intron is not cut, so one molecule has
one length on a candidate however its mates lie, and a contradicted candidate's covered bases are priced as overhang.
Every length histogram comes from the payload, so a mixed-frame length is unrepresentable.

### 3.6 Arbitration and the side buffer

A fragment is offered once, after its last hit (`AccumulatorSet::offer`), with every hypothesis about its unsequenced
gaps — the empty one is the unspliced, genomic explanation — and every resolved hit as a placement, intergenic ones
included, since gDNA sits anywhere. The strand check and the length clip exclude a placement alone. One surviving
(placement, hypothesis) deposits; two or more hold the fragment WHOLE in the side buffer
(`deferred_undetermined_gap`, `deferred_multiple_placements`), and `deposited + deferred_* + dropped_* == offered`
counts fragments. A held record keeps each placement's extent, strands, introns, NM, clipped bases and hypotheses, in
one bank sorted on content. The scan also reads the library's mismatch rate ε — NM over the read bases its unique
unspliced hits align — onto the payload, so a cached scan drains at its own rate.

---

## 4. The second pass

`pipeline.drain_side_buffer` → `second_pass.score_held_fragments` → `choose_hypotheses` → `drain`.

- **It runs before calibration, and that is the structural decision**: every input of the score is pass one's, so
  calibration runs once, on the complete tally. Fit once, score once, drain once: the confident set is biased short,
  and feeding the drained fit back would prefer shorter, more-spliced paths in a loop that can run away.
- **One tally path**: the drain re-enters `deposit` with the chosen placement and hypothesis alone; `drain` is pure,
  so a cached scan drains at any seed. Afterwards nothing is held; pass one's numbers live in `DrainQC`.

> **The score**: a held record's candidates are scored by the EM's own emission, with pass one's traffic as the prior,
> `score(p, h) = A_p · ρ(h) · f(L) · s · e(p)`, formed in logs and normalised within the record
> (`second_pass.combine_factors`; `EQUATIONS.md`, the second pass's score).

- `ρ(h)` is the bottleneck (minimum) of pass-one inverse-length density over the junctions a spliced path uses, or the
  boundaries the unspliced path crosses where the others jump. `f(L)` is the RNA law on a path the deposit books
  spliced and the unconditional anchor otherwise. `s` is κ or 1 − κ for RNA, ½ for the genomic hypothesis; a placement
  whose observed intron is an annotated junction is certified RNA under every path.
- On a multimapper's placements only (on one placement they are common factors): `A_p`, the pass-one unique traffic
  at the objects the placement's path deposits on (`placement_abundance`), and `e(p)`, the hit's sequence at the EM's
  prices (§9.5).
- **gDNA's floor**: a copy whose fragments all multimap reads no pass-one traffic. gDNA is one uniform template before
  capture and capture only enriches, so a placement gDNA can be at reads `max(A_p, bulk)`, the bulk being the unique
  contained fragments per observable start over intergenic regions (`pass_one_bulk`); a spliced placement keeps its own.
- A factor zero for every surviving candidate is dropped as uninformative; one zero for some stays decisive
  (`TRAPS: an-all-zero-factor-is-inert`). The prior is pass one's traffic, not the EM's θ or the capture efficiencies,
  because the draw is calibration's input.

**The draw** is one multinomial per record from one seeded stream, keyed on queue position and never on content (a
content key would collapse a 60/40 posterior over identical fragments to 100/0), landing only on a candidate with mass;
one hypothesis wins the whole fragment (`TRAPS: fractional-mass-is-the-problem`), and the EM's later allocation is not
reconciled with it. What the floor leaves is issue #18; gDNA drawn onto a spliced path is
issue #17. **Instruments read the drained tally**: they drain cached scans at the
production seed (`scan_cache.calibration_inputs`) and lift oracle partitions into that frame by replaying the whole's
drawn choices (the oracle in `scripts/design/_oracle.py`).

---

## 5. Calibration: from the tally to one slot

`calibration.calibrate` owns the order and computes no model of its own:

```
substrate → chain, geometry, statics → strand model (κ; binomial) → G1/G2/G3 initial belief
  → pass 0: solve_chain with no fitted prior
  → refits (calib_refit_iters): fit the gDNA landscape on the belief, reset the belief, solve_chain again
  → deconvolved mass per region and boundary → capture efficiencies → CalibrationResult
```

### 5.1 The layers

`calibration/_layers.py` is authoritative: an import points down a layer or sideways, never up, and every module is
placed (`tests/calibration/test_layering.py`).

| layer | modules | a change about … |
|---|---|---|
| 0 · vocabulary | `errors` `signature` `region_chain` | the region bitmask, the N E N … chain |
| 1 · payload view | `splice_graph` `substrate` `region_arrays` | what a fragment tally MEANS |
| 2 · opportunity | `effective_length` `capture_eff_length` `sj_opportunity` `gdna_opportunity` `gdna_density` `fl` `shared_length_law` | how many places a fragment COULD have sat |
| 3 · geometry + per-slot solve | `region_geometry` `simplex_logodds` | one slot's own numbers, and ψ |
| 4 · strand | `strand_balance` | which strand a fragment came from |
| 5 · density and prior | `landscape` `capture_efficiency` | how dense a component is, and the priors |
| 6 · the solve | `diagnostics` `sweep` `messages/` | what one neighbour tells another |
| 7 · assemble | `calibrate` `priors` `result` `track` | turning the solve into a result |

### 5.2 Calibration is message passing

> **An exon's gDNA level cannot be measured. It can only be imputed from its neighbours, and imputation across the
> chain is message passing.**

gDNA is directly measurable only where no mature transcript crosses: intergenic and intron regions and the boundaries
against them. An exon holds gDNA plus RNA; its `intron|exon` flank's gDNA density bounds it from below and the flank's
certified sj flux from above. That is a fact of the annotation and the deposit rule, not of an estimator, so the solver
is a sweep over the chain, not a per-object fit.

- **RNA is the residual and is never predicted**: gDNA is near-uniform and measurable before any solve, RNA spans six
  decades with no genomic autocorrelation. Under capture an exon's roles exchange — its RNA is measured locally by the
  boundary's flux — so the residual is whichever component the data do not measure there.
- **The boundary axis splits on whether mature RNA can cross**, not on whether a junction attaches: an `exon|intron`
  boundary is near-pure gDNA under sparse nascent RNA, an `exon|exon` boundary is crossed freely
  (`TRAPS: mature-rna-never-crosses-a-boundary`).
- **Pass 0 bootstraps the prior**: it must solve enough slots, not all, to train the landscape (§7.1). **An imputation
  is weak by definition**: every step charges an honest width, and a regression on stranded data is a precision defect
  until proven otherwise (`TRAPS: an-imputation-must-cost-something-every-hop`). Strand specificity is a spectrum, so
  nothing switches off at a threshold.

### 5.3 The geometry and the strand model

`region_geometry.build_region_geometry` builds every slot's gDNA and RNA opportunity before any count is read and owns
every divisor below it — contained placements at a region, crossing placements at a boundary, each on its component's
own law — and the published gDNA lengths are projected off it. Precision is per component, never one scalar. A known
bias: an unspliced crossing's RNA half takes an unbounded reach, over-calling gDNA before a polyA site. A TOTAL
abundance is a pre-solve instrument only. The initial belief is the signature class: G1 locked all-gDNA and never
solved, G2 one strand and resolved by it, G3 both strands (AMBIG).

- **κ, the RNA sense fraction, is the genuine junctions' rate.** A stranded library's annotated junction is genuine
  RNA (at the library's minority-strand rate), a splice artifact (misaligned gDNA, ½) or reversed; the three-class
  likelihood over the per-junction strand table is maximised jointly, and κ is the `Beta(1, 1)` posterior mean over the
  reads its maximum credits to RNA (`junction_fit.genuine_sense_fraction`, `strand_balance`). Every κ consumer reads
  this one count.
- **The strand overdispersion is 0 for both components — binomial — by policy** (`calibrate._sweep_inputs`): no
  population measures one cleanly, and an od set too high mutes the channel. It returns only when proven needed, first
  as an RNA-only value from the deepest junctions, developed on real data across libraries.
- **The strand channel is live iff the protocol preserves strand** (`strand_balance.strand_discriminability`): the Bayes
  factor of a free κ against κ = ½ on the spliced 2×2, `ln BF = N·ln 2 + ln B(κ̂(N + 2), (1 − κ̂)(N + 2)) > 0`, with no
  multiple of σ chosen and no gDNA term. A library with no spliced read is unstranded, never an error. gDNA's strand
  term is ½.

### 5.4 ψ, the per-slot solve

ψ (`simplex_logodds`; natively `psi_kernel.h`'s `solve_slot`) is the posterior over a slot's composition on the
2-simplex, in λ, the gDNA-vs-RNA log-odds, and θ = arcsin τ, the RNA tilt: the strand term, the gDNA arm (the symmetric
Jeffreys reference, plus on the refits the landscape at the density each cell implies), the RNA arm (the reference
alone) and the delivered message rows. It never reads the spliced mass.

- **ψ has no reference location**: where the strand channel is silent a location term would be the whole answer at
  any depth, so none returns as a flag, constant or weakened variant; background reaches a slot only through the
  landscape and the messages.
- **One λ lattice for every consumer** at `CalibrationConfig.sweep_logodds_step` over `[−L, L]`, `L` the larger of
  `sweep_logodds_window` and the landscape's derived demand (`DensityLandscape.required_logodds_window`); the step is
  fixed, so a wider bracket never coarsens it. **The θ nodes follow the strand term's peak** (`tilt_window`: a window
  within `T = −log ε₆₄` nats of the maximum, `⌈2T/π⌉ + 1 = 24` trapezoid nodes), so the marginal is exact at any depth.
- **The tilt atom**: an AMBIG slot's tilt is {pure +, pure −, mixed} at equal reference weight, two columns at τ = ±1
  beside the continuum (`slot_cube`); a held RNA level on a strand rules the other strand's atom out
  (issue #42).
- **The read-out** (`compose`): `f_g` the continuous ½-quantile of the λ posterior, interpolated on λ
  (`TRAPS: interpolate-on-the-axis-where-the-lattice-is-uniform`), with `Var(log f_g)` its moment and the RNA split
  `(1 − f_g)·w_+`, `(1 − f_g)·(1 − w_+)`, so the composition closes identically. The median, because it is closer to
  truth at both vertices, where most error lives. A slot with no counts publishes `(0, 0, 0)`.

---

## 6. Calibration: the solve

### 6.1 The backbone and the policy

| file | role |
|---|---|
| `sweep.py` | THE BACKBONE: the structure, the library, the locus blocks, one native call (`native.solve_blocks`), the assertions; it knows nothing of capture, splices, levels or lanes (`test_sweep_backbone.py`) |
| `diagnostics.py` | `SweepDiagnostics` and `CalibrationDiagnostics`, the instruments' records of a solve; production passes none |
| `messages/__init__.py` | the interface (`Policy`: a name, a strand model, a library) and `ChainView`, with no belief field, and `view_fields`, the per-slot arrays the kernel takes |
| `messages/silent.py` | `SilentPolicy`, which sends nothing: the measured floor (`message_policy = "silent"`) |
| `messages/transfer.py` | `TransferPolicy`, the shipped default; its builders, passes and solve are `native/transfer_kernel.h`'s |

- **A terminal receives nothing.** A region admitting no RNA strand (`is_region & g1_locked`) is fixed before a message
  exists and no step enters it, so the chain breaks into independent loci there and is solved a block at a time
  (`region_chain.locus_blocks`), each block reading the terminal beyond its end.
- **One native call, one path.** `native.solve_blocks` runs each block end to end on its own thread and arena and
  reduces nothing across blocks but integer counts, so the answer is bit-identical at every thread count and block
  size (`sweep_block_slots` and `n_threads` are resources). The gates read the one implementation through bindings
  (`_solve_impl.testing`'s `transfer_prepare`, `transfer_pass`, `transfer_solve`; `tests/calibration/_transfer_harness.py`).
- **Every sweep runs the whole layer**, which reads observations, geometry and the library, never a belief or the
  landscape; a refit's messages equal the last refit's and are recomputed, never cached.
- **The one contract**: a message may use the destination's constants and observations, never its belief
  (`TRAPS: a-message-from-the-destinations-belief`); a recipient may only widen an arriving claim by its own counts. The
  kernel builds each row from the source's claim and what the source holds, so this holds by construction.
- **The assertions** (`sweep.AssertionCounts`): the kernel sees only the two neighbour states, and every delivered row
  is one finite row per slot on the grid. A violated assertion raises; each assertion counts its eligible slots
  (`TRAPS: could-the-arm-have-fired`).

## 6b. The message layer

`messages/transfer.py`'s docstring states the transfer policy in full; these are the rulings that bind it.

**A step carries a LEVEL or a COMPOSITION.** A level (`ρ_g`, gDNA per placement) survives population changes —
termini, strand flips, splicing — since gDNA is genomically continuous, and dies at an enrichment change. A composition
(`f_g`, gDNA's share of the crossing population) survives enrichment, which lifts both components alike, and dies at a
population change. One always survives, so the licence question is the currency question, asked per face kind —
`object class × {splice site, terminus, both}` from the flag bits (`TRAPS: an-object-class-does-not-see-a-terminus`). A
message's population is whatever enters the destination: the splice-in includes the spliced fragments, the splice-out
excludes them. Level and composition transfer are two point hypotheses about one log enrichment; nothing switches.

**The recipient decides, by formal forward-backward.** A sender publishes its claim unchanged; the recipient forwards,
modifies or stops it. `pass_block` runs a forward then a backward pass, each recipient holding its neighbour's claim
composed with what that neighbour holds from its far side, so no node hears its own claim back; silence and no
neighbour are the received table's two states (`ReceivedView`). `solve_block` solves each node once from its evidence,
its two held messages and the gDNA prior. **The completion contract**: no slot is skipped; a message crosses a run of
`exon|exon` boundaries charging each step its own width, and where composition cannot cross the gDNA level still does.

**The composition rules** are a kind and its parameters at `(destination, side)` — FORWARD, TRANSPORT, SPLICE-OUT,
EDGE, LEVEL (`transfer_rows.h`) — at most one per face, optional rows beside a presence mask read first:

- **intron ⇄ its boundary: FORWARD.** One unspliced population, so the density odds are shared and the count log-odds
  shift by `log(Eg_i/Eg_b) − log(Er_i/Er_b)`, adding no width; into the intron only under a shared single strand
  (`boundary_shares_strand`) and a live strand channel.
- **boundary → exon at a licensed face: TRANSPORT, the splice-in map** (`splice_faces`, `face_map_lambda`,
  `transport_row`; licensed is no terminus and one strand set, `face_is_licensed`). The certified flux caps the exon's
  claimable gDNA share, the row widens by `trigamma(n_u + ½) + trigamma(n_s + ½)`, and two faces sum as independent
  witnesses. Registered by geometry, never gated on the strand channel, so it survives where that channel is dead.
- **exon → boundary: SPLICE-OUT, the map read backwards** (`splice_out_row`): the mature share is removed and the
  enrichment cancels. The exon sends only its own strand row where its channel is live, so an unstranded exon says
  nothing; the width is the marginal over the measured ratio's log, at `trigamma(n_s + ½) + trigamma(n_u + ½)`.
- **A terminus** (`terminus_rules`, `outside_flank`): the outside flank is read off the flag (TSS+ / TES− bodies extend
  right, so the outside is the left flank) and licensed by the splice-out law with the spliced crossing's RNA rate. The
  inside region takes LEVEL, the level rule: the boundary's own strand profile through the level-kept map, blurred by
  both totals' counting and the pair's own discrepancies. A level comes from the sender's measurement only.
- **The alternative splice site** (`junction_flanks`): each flank pair holds two witnesses of one quantity, and where
  they disagree beyond counting the pair's messages widen by `max(0, d² − v)` (`blur_row`); nothing is shifted or
  pooled. At an sj+terminus boundary the terminus decides the side and the flux sits on the junction's exon side
  (`junction_exon_side`); a flux is never a crossing.
- **The `intergenic|exon` edge → exon: EDGE, one-sided** (`edge_level`): at least the edge's gDNA density at Poisson
  width and nothing above, since no local witness prices capture's lift of a probed interior; a zero count is vacuous
  (issue #41).

A blur preserves relative likelihood in its tails: small sums use a local log-sum-exp, no floor flattens distinct
tails, and an exactly-zero window stays zero.

**The gDNA level lane: gDNA is always conveyed, as a lower bound** (`gdna_lane`). A level is a profile over
`u = log(ρ/ρ_ref)`, `ρ_ref` the library's structurally pure gDNA density, so it needs no map and crosses a node with no
total. It is the default rule of every face without a composition rule, so a stop by omission is impossible. An empty
node forwards it; a full node emits the INTERSECTION of its own lower side and the held level (never the product, which
ratchets); the recipient prices the step (`hop_price`) and takes the level as a lower bound only — darkness under
capture is not absence. A held composition is never re-issued as a level. A two-sided bound needs a witness that the
library's gDNA is enriched, which no message has (issue #43).

**The RNA level lanes** (`rna_lane`). A both-stranded slot's split pins only a line in `(f_g, f_+)`; since
`½ f_g + κ f_+ + (1 − κ) f_− = p_+`, a lower bound on one strand's RNA is an upper bound on the gDNA share, so the three
lower bounds and the node's counts bracket it from both sides.

- One lane per strand (`LevelsView`), one shared RNA origin. Strand `s` crosses a face iff the boundary carries none of
  `s`'s bits and both nodes admit `s`; across `s`'s junction it enters `s`'s intron, not its exon; a terminus of `s`
  stops it. Two-sided only between an intron of `s` and its own boundary.
- Every step pays the pair's price, both column counts' counting plus their disagreement on `count − other`, where gDNA
  cancels. Sources: a single-strand node's own claim (`rna_level_of_profile`) and the certified flux at an exon's
  junctions, lower-sided and one step (`flux_level`), priced on the exon's column count over the protocol's share of
  the opportunity (`read_column`).
- Delivered as one row over ψ's cube at AMBIG nodes (`simplex_logodds.CubeRows`), and as a ceiling at single-strand
  nodes (`rna_row_of_level`) only from a face that sent no composition. The tilt has no lane, and no RNA level is read
  from a belief (`TRAPS: one-hop-lifted-out-is-still-the-relay`).

**The library** (`TransferPolicy.library`): the two coordinate origins and whether the strand split is a live witness,
the only things a message knows outside its block. **Every message idea is judged twice, never pooled**: at pass 0
(`calib_refit_iters = 0`) and through the full pipeline, since the landscape trains on pass 0 and returns its errors
everywhere (issue #23). Unstranded rows must win against `silent` and
stranded rows do minimal harm (issue #26).

---

## 7. Calibration: the landscape, the capture efficiencies and the result

### 7.1 The landscape prior

The landscape (`landscape.fit_landscape`) is `P(log ρ_g)` over the slots, fitted on the previous solve's deconvolved
gDNA and read by ψ's gDNA arm on every refit; nothing else carries between refits. The estimator is
component-agnostic; who trains it is the caller's (`calibrate._fit_gdna_hyperprior`;
`tests/calibration/test_landscape_training_population.py`):

- **Regions train; boundaries and AMBIG regions do not**: a boundary's two-flank mixture fills the valley between the
  modes, and AMBIG is what the prior resolves. A live region with no unspliced mass is a zero-count ANCHOR and trains
  regardless: the strongest depletion evidence there is.
- **A node whose only evidence is a bound does not train**: `RegionBelief.has_composition` (an own channel, structural
  certainty, or a received composition row) selects; levels, ceilings and cube rows are bounds.
- **The location floor**: a slot trains only where `Var(log f_g) ≤ 1 nat²` (`CONSTANTS.landscape.located_var`, the
  one-fragment wall through `Var(log c) = 1/c`); above it precision is a continuous weight (`landscape._reliability`).
- **The grid spans every slot the prior is read at** (`fit_landscape(domain=…)`), so the training cut never moves
  the axis; a location-free kernel is placed by the previous refit's landscape (`landscape.estep_kernels`), a counted
  one keeps its own location.

A provisional allocation is not an observed origin: no first landscape trains on intronic observations as known gDNA.
Weak-strand nascent RNA booked as gDNA can become the bulk (issue #24); the
prior's double entry into ψ is issue #21.

### 7.2 The capture efficiencies

**Capture is a spectrum**: no capture detector, library-level test or on/off status is built, because a plasma panel
enriches scarce transcripts that stay a tiny fraction of the library. Every capture weight is per object or placement
and continuous, and no enrichment is the limit where every weight is 1. **Capture has a local footprint**: a probe
binds any fragment overlapping it, so enrichment is shared with adjacent regions, boundaries and the connected
junction, represented at the object resolution, never by reconstructing probes or by a fixed cutoff.

**Capture inference is probe-free** (owner reaffirmed 2026-11-11): Rigel learns enrichment from the reads and
annotation; it must not depend on a bait-sequence file, probe coordinates or panel-design metadata. Probe designs
used to construct validation truth do not become inputs to the production estimator.

> **The reader** (`capture_efficiency.capture_efficiencies`, after the last refit). Each object's efficiency is its own
> gDNA density read against the BULK, the peak of the landscape's largest basin, since gDNA off the probes always
> outnumbers gDNA on them. The prior is Efron's two groups — the bulk basin collapsed onto its peak, elsewhere the
> population's own density; the likelihood is the object's own Poisson count, a region's contained, a boundary's
> crossing; the read-out is the posterior mean of the log density, saturated at its enriched basin's peak.

No mode is located and nothing gated: a one-basin landscape, no landscape, and an object with no gDNA support read 1.
The crossings are read once, by their boundary. The unit is the bulk (`gdna_bulk_density`), so an enriched object reads
far above 1, and the fully captured level is the saturation.

> **The weight** (`capture_efficiency.visible_share`, `weighted`; the owner's go, 2026-11-15). The published reading is
> the reader's raised to one library-wide weight `w` (`gdna_capture_weight`), every efficiency and the saturated level
> alike: `c → c^w`. `w` is the share of the one-strand unspliced mass whose object would show a gDNA fragment of its
> own (EQUATIONS §11). Where gDNA is dense enough to be seen where the library's reads sit, `w` is near 1 and the
> reading stands; where the reads sit on objects holding none, the few objects whose gDNA shows cannot speak for them,
> and the claim fades toward no capture anywhere (`w = 0`, the no-capture line).

The weight carries the release ruling of 2026-11-14: right wherever gDNA is easy to see, 5 % and above; where it is
sparse the capture claim shrinks to 1, keeping genes and giving up transcripts. It is one weight because capture is a
library pattern, and because the hole between 0 and 1 % gDNA is a transcript read inconsistently: some of its objects
hold real gDNA evidence and read thousands-fold, its RNA-heavy neighbours read 1, and shrinking each object toward 1
from its own evidence leaves that in place. It is continuous and has no constant: no detector, switch or threshold. Refused on the screen: a weight scaled by the landscape's spread (one probe
reads as weak capture and fades), and one from ψ's own variance (an unstranded slot says nothing of its own gDNA, so
unstranded capture ON fades). What it costs is issue #2.

False gDNA on a near-gDNA-free library still reads as capture, faded by the weight but not removed
(issue #1); that off-target cDNA is depleted as gDNA is remains untested
(issue #46).

### 7.3 The result

`CalibrationResult` publishes per region and boundary the deconvolved gDNA and RNA masses (certified spliced crossings
join RNA whole), the composition, the gDNA lengths, gDNA's conserved share at every boundary
(`gdna_boundary_conserved_len`, unbounded where unmasked, since gDNA's template does not end), the efficiencies at the
capture weight and the weight itself, κ and the QC density `gdna_density_global`, a ratio of sums. A zero-gDNA library
and a pure-RNA object are valid outputs.

---

## 8. The EM's two length tables

### 8.1 The five length pools

Each pool is selected by structure, never by a fit, and a fragment enters when exactly one hypothesis survived —
determinacy, not provenance (`TRAPS: a-purity-filter-is-a-length-filter`):

| pool | rule | component |
|---|---|---|
| `DNA_INTERGENIC` · `DNA_INTRONIC` | contained in one intergenic · intronic region | gDNA, with unannotated · intronic RNA |
| `DNA_INTRON_EXON` · `DNA_INTERGENIC_EXON` | crosses one boundary, flanks {intron · intergenic, exon} | gDNA: mature RNA never crosses there |
| `RNA_SPLICED` | uses an annotated sj on its one surviving path | RNA |

An exonic contained fragment or a multi-boundary crossing is a mixture and has no pool. The contained pair dominates
off capture, the crossing pair under it.

### 8.2 The two measured laws

`fl.build_fl_models`, from one payload in one frame (`pipeline.library_fl_models`):

- **gDNA's law** (`gdna_pmf`, the uniform frame every opportunity reads): each gDNA pool divided by its own opportunity,
  then combined, never pooled raw (`TRAPS: opposite-tilts-must-not-pool`; `TRAPS: divide-by-a-probability`). A two-pool
  CONTRAST divides out the contained pair's contaminant, its purities from `gdna_density`'s one-sided rate
  (transcription only adds, so clean objects sit low); it fades to the pools' own mixture as their separation falls
  into noise. The law is shrunk toward the anchor (`CONSTANTS.fragment_length.pool_prior_ess`;
  issue #20).
- **RNA's law** (`rna_pmf`): the spliced pool's own counts at their sj opportunity (`fl.own_count_law`), Poisson at
  `N·π(l)·law(l)` with the law the anchor times a random-walk deviation at its evidence maximum. The anchor
  (`global_pmf`, every deposited fragment) is only what an empty pool returns, read as no RNA measurement.

Calibration and the priors read both laws refit on the drained tally; the second pass reads pass one's.

### 8.3 One shared law unless the data show a difference capture cannot make

The E-step scores an unspliced fragment's length by `f_R(l)/f_G(l)`; a placement's capture weight is common to both
origins and cancels, so capture enters only the normalisers. But gDNA's law is read off target and RNA's census on
target, so a scorer handed both reads what capture did to RNA's census as chemistry.

> **RNA's law is both origins' table by default, and gDNA's moves away from it only as far as the data show a
> difference capture cannot have made** (`shared_length_law`, decided after calibration in `pipeline.length_tables`;
> `EQUATIONS.md`, the EM's two length tables).

1. **The capture mixture** `G_c = (1 − c)·g + c·g·l/μ_g`, `c` the captured fraction that best explains RNA's census in
   `[0, κ_c]`, `κ_c` RNA's captured fraction under calibration's efficiencies
   (`capture_eff_length.rna_captured_fraction`); capture only lengthens, so `c = 0` where RNA's census is not longer.
2. **The decision** `w = clip(D/B − 1, 0, 1)`: `D` the Hellinger distance between the measured laws less sampling noise,
   `B` the most capture at `κ_c` can make over every response that rewards length at most proportionally; gDNA's table
   is the geometric path from `rna_pmf` (`w = 0`) to `G_c` (`w = 1`). Off capture `B = 0` and a resolved difference
   keeps each origin's law; an RNA law no spliced fragment measured gives one law (`rna_measured`).

A decision, not a blend: a real gap partly shrunk shifts the RNA/gDNA balance as surely as a false gap kept, so the
ramp is whole once `D − B ≥ B`; its form is the design's one choice. It reads no strand
(issue #27). RNA's normaliser carries each efficiency to the tables as
`c·share(g)/share(T_G)`, so at one table it is `Σ_o c_o·share_o(g)`, the null model exactly.

---

## 9. The locus prior and the EM

### 9.1 Loci, components and the scorer

`routing.FragmentRouter` scores every buffered fragment against its compatible transcripts (`scoring.cpp`): a
deterministic unique spliced fragment is assigned directly, every other becomes an EM unit. Loci are connected
components of transcripts linked by shared fragments (`locus.build_multi_loci`), each holding its transcripts,
annotated and synthetic alike, and one gDNA component. A template no hit of which reaches a transcript is intergenic
gDNA. A candidate's terms are its length table, strand, coverage weight and per-base overhang penalty
(`FragmentScoringConfig.overhang_alpha`), plus a place weight and sequence on a multimapper; one under
`pruning_min_posterior` of the unit's best is pruned (issue #32).

> **gDNA can explain an alignment it can produce: unspliced, and implicit (a pair whose gap could hold an annotated
> intron). A multimapper can be gDNA if any of its alignments can; a sequenced junction cannot. gDNA then COMPETES as
> an RNA candidate does** (`scoring.cpp`'s `gdna_can_explain`, `gdna_competes`).

gDNA's term is the unspliced term at the genomic footprint whatever the label, and a footprint its length law cannot
produce is pruned with no length cut of its own.

### 9.2 The capture-contracted length — one shared rule

> **Every EM component — the locus gDNA, every synthetic span, every annotated transcript — takes its contracted
> length by one rule: the sum over every object its fragments deposit on of its CONSERVED SHARE there times that
> object's CAPTURE EFFICIENCY, `L_T = Σ_pieces S_T(p)·c_p + Σ_cuts M_T(k)·c_k`** (`EQUATIONS.md`, the conserved frame).

Pieces are the regions a template covers and a cut is the boundary or junction between consecutive pieces; the shares
follow the deposit rule's own split (`effective_length.conserved_cut_shares`) and total the fl-marginal length exactly.
Transcripts scale their length by `Σ share·c / Σ share` (`capture_eff_length.transcript_capture_eff_lengths`); the
gDNA component sums its locus's region supports and boundary shares at their efficiencies (`priors.assemble_priors`).
With no landscape every length is its uncontracted span, bit-identically.

- **A junction is priced by conservation of bases**: `c_junction = c_lo + c_hi − ½(c_intron,lo + c_intron,hi)`, within
  `[0, max(the junction cap, the best object it touches)]`, since a fragment binds its best single probe part
  (`capture_eff_length._junction_prices`). **Capture is local**: only objects within a fragment of it enter.
- **The junction cap is RNA's** (`capture_eff_length.fit_junction_cap`, landed 2026-11-13). gDNA never crosses a
  junction, so the landscape's saturated level only bounds what gDNA shows. The cap starts there, at the capture weight
  (§7.2), and is raised as far as the library's junctions support, each read on its own footprint: its conserved spliced
  mass against the unspliced RNA of the exon regions beside it that its transcripts alone cover (`EQUATIONS.md` §11).
  The counts never lower it, and with no level above the bulk (a saturated level of exactly 1) it is not fitted: no
  capture is read, so there is none for a junction to exceed. The lengths and a multimapper's place weights read the one
  value (`summary.json`'s `junction_capture_cap`).
- **A junction's price is never pooled**: some are probed by a probe spanning the junction and some are not, so each
  keeps its own two boundaries and intron, no junction's price is evidence for another's, and any gain from pooling
  them is chance another panel reverses. The cap is the rule's one library level, as the saturated level was, and its
  fit reads every junction on its own footprint. **Capture efficiency stays outside the EM**.
- **The yield's two consumers**: the contracted length is a YIELD and enters only the E-step, where only ratios within
  a locus decide, and the locus prior's gDNA length. TPM uses the plain fl-marginal length, since Rigel never knows
  which transcripts a panel probes. The yield has no floor — a component with no start position cannot emit — and its
  uncertainty is not propagated through the EM.

### 9.3 The locus prior

> **The EM's gDNA share is gDNA's share of every fragment in the locus, and both sides of it are one tally**:
> `P_g = U·G_c/T_c`, `P_R = U − P_g` (`pipeline.em_pseudocounts`), with `G_c` calibration's gDNA count of the locus,
> `T_c` its count of every fragment it deposited on the same objects, and `U` the units with a gDNA candidate. A
> fragment is in the numerator only if it is in the denominator (the owner's ruling, 2026-11-14).

The deterministic spliced fragments skip the E-step but enter every M-step, so the prior is neutral iff
`P_g : P_R = G : (N − G)` (`EQUATIONS.md`, the pseudocounts' odds). The strength is one pseudo-fragment per eligible
unit, with no knob and no derivation (issue #30). The RNA pseudocount goes to
every RNA component in proportion to its evidence, so a component with none receives nothing.

**The count** (`priors.assemble_priors`) holds each gDNA candidate fragment of a multi-locus once: the regions it
overlaps and every boundary touching them, the outer two included, since a fragment crossing a locus's boundary is its
candidate; only the crossing term converts, by the boundary's own `mass / count`. The pooled `q` is
issue #28.

**The fragment count** (`LocusPriors.fragment_count`) is every deposit on the same objects by the same shares: each
region's whole contained count, each boundary's two banks at their own conversions, and each sj's conserved mass at
the locus of the exon regions it joins, which the result names (`CalibrationResult.sj_src_region`, `sj_dst_region`).
Until 2026-11-14 the share's divisor was the EM's count of the locus's unique fragments, while the gDNA count held
every multimapper at its drawn placement. On an aligner's input the multimappers join most of the library into one
locus (61–96 % of the located fragments on the STAR rung) and are 4–7 % of it, so the share read that much too high:
1.030 at STAR g98 stranded capture ON, clipped to one, an RNA pseudocount of exactly 0 on 7.66 M fragments. The
simulator's alignments hold no multimapper and could not show it. A gDNA count above the fragment count is now
refused (`EQUATIONS.md`, the pseudocounts' odds).

### 9.4 A multimapper's places, and the intergenic background

> **A fragment's likelihood under a component sums, over its places, each place's term times its capture weight; hits
> overlapping on one reference are one place, read at their best** (`scoring.cpp`'s `merge_places`). **A place's weight
> is the deposit rule's mass at the reader's efficiencies, `w = Σ_o m_o·c_o`, on every candidate there**
> (`capture_eff_length.place_log_weights`). **A hit with no transcript candidate is an intergenic place: uncaptured, at
> calibration's bulk gDNA rate, a per-unit background no locus carries, never pruned** (`em_locus.h`'s
> `IntergenicTerm`; `EQUATIONS.md`, a multimapper's places).

Intergenic is uncaptured (`priors.intergenic_gdna_rate`) because an aligner puts captured genic fragments on
intergenic copies, and their efficiency would turn that artifact into gDNA. The rate is fixed from calibration, never
re-estimated across loci, which would couple every locus to every other. The background's mass is intergenic gDNA in
`summary.json`, never in the loci table; a unit it wins is stamped `ZF 0x25` on every record and `ZH 1` on none.

### 9.5 A hit's sequence

> **At the library's mismatch rate ε, every candidate at a multimapper's hit carries the hit's sequence against a
> matching aligned base, `NM·log[ε/(3(1 − ε))] + clip·log[1/(4(1 − ε))]`. A one-hit fragment reads neither: its
> sequence is common to every candidate** (`scoring.sequence_log_penalties`, `scoring.cpp`'s `score_mm_alignment`;
> `EQUATIONS.md`, a hit's sequence).

No unique unspliced base means no price; `ε = 0` makes a mismatch impossible; `ε ≥ ¾` is refused. The draw reads the
same prices (§4). The clip price can only favour the contiguous reading, the one gDNA can explain: a short-anchor
contiguous alternative clips what it cannot place, and a free clip would read it as a false tie
(issue #7).

### 9.6 The solver and whole counts

- VBEM by default (`EMConfig.mode`), SQUAREM-accelerated, from the coverage-weighted warm start projected through the
  calibration prior. **The accelerator never decides the answer**: a jump that would carry a live component below the
  floor is shortened toward the plain step, never clamped (`backtracked_squarem_step`;
  `tests/test_em_start_independence.py`), because a clamp is permanent.
- **VBEM's start-dependence is an accepted limit**: the digamma penalises a small component, so components sharing
  ambiguous fragments race and the start picks the winner; the start is fixed, so one input gives one answer. The
  iteration budget can stop the EM short of convergence (issue #34).

> **Whole counts are the default (`assignment_mode = "sample"`), fractional an option; there is no argmax mode and no
> cutoff on one fragment's posterior.** `sample` rounds each component's fractional count within its locus by largest
> remainder, draws each fragment from its posterior re-weighted by count still owed over mass still to come, and
> repairs a stranded fragment by moving a chain of fragments onto the counts, or onto counts within one of every
> fractional count where those are unreachable (`em_assign.h`'s `assign_posteriors`; `EQUATIONS.md`, whole counts).

It is feasible, not optimal, and is kept over the exact optimum per read because that optimum gives up calibration:
each component gets a representative sample of the fragments it could have produced. A molecule belongs to one
component, so a whole count shows read by read in the annotated BAM, and whether a transcript is real shows in its
total. The counts are the rounded fractional counts whatever the seed (`EMConfig.seed`, fixed).

---

## 10. Outputs

`rigel quant` writes `quant.feather` (annotated transcripts), `gene_quant.feather`, `nascent_quant.feather` (synthetic
spans), `loci.feather`, `fragment_lengths.feather`, `calibration_track.feather`, `summary.json` and `config.yaml`;
`--tsv` mirrors the tables.

- **RNA is reported as annotated or nascent**, never as a solver population: annotated is every reference transcript,
  single-exon included (it does not splice, so its nascent and mature sequence are one); nascent is the gene-neutral
  synthetic spans. A transcript names its span (`nascent_id`); `tpm` is over annotated transcripts, `tpm_total_rna`
  over transcripts and spans. TPM's length is the observable fl-marginal `effective_length`, floored at 1 bp; the EM's
  yield is `em_effective_length`.
- `summary.json`'s pools `count_annotated`, `count_nascent` and `count_gdna` sum to one in fraction, `count_gdna` being
  the EM's genic gDNA (`count_gdna_em`) plus intergenic gDNA; it records the read-pair outcomes, the observable genome,
  `fragment_stats.mismatch_rate`, κ pooled and genuine, the length decision, and the bulk and saturation.
- **The annotated BAM** (`--annotated-bam`) re-judges every pair with the scan's code and stamps each record
  (`annotate.py` is the schema): the winner (ZT, ZG, ZR, ZI, ZJ), the outcome bitfield (ZF), the posterior (ZW), the
  ambiguity class or a withheld pair's reason (ZC), the winning placement (ZH), the candidates (ZN), the record's own
  splice type (ZS) and the locus (ZL).
