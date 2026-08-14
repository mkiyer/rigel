# DESIGN — what is built, and the rulings behind it

**Everything here is SETTLED.** These are decisions that have been made, measured, and are not to be
re-litigated. Where a decision was reached by refuting an alternative, the refutation is in `TRAPS.md`,
not here. Derivations are in `EQUATIONS.md`.

---

## 0. ⭐⭐ VOCABULARY — one word per concept, owner-set 2026-08-04

⛔ **This table is binding on prose, comments, docstrings and identifiers alike.** Every banned word below
was in use for a concept that already had a name, and the ambiguity cost a reader real time.

| ✅ the term | what it means | ⛔ the banned synonym |
|---|---|---|
| **REGION** | a **contiguous genomic interval** — a region of the partition. Has a length in bp | `region` is tolerated: it is the index's own word for the same thing |
| **BOUNDARY** | a **single genomic position** separating two adjacent REGIONs. Zero bp wide | `boundary` · `boundary` · `crossing` (as a noun for an object) |
| **BOUNDARY** | ⭐ an accepted synonym for BOUNDARY, kept because the code and the earlier docs are full of it | — |
| **SPLICE JUNCTION** (`sj`) | ⚠ a splice junction — a *different* object on a *different* axis, directed `src → dst` in **GENOMIC** order (`src < dst`, always). Never a synonym for BOUNDARY | — |
| **the sj's LOW end / HIGH end** — `_lo` / `_hi`, or `boundary_left` / `boundary_right` | the two BOUNDARIES a sj is anchored at, named in GENOMIC order, which is the order the code stores them in | ⛔ `donor` · `acceptor` — see below |
| **slot** | one entry of the chain, which alternates REGION, BOUNDARY, REGION, BOUNDARY … A slot is a REGION **or** a BOUNDARY | — |
| **step** | one adjacency move along the chain: REGION→BOUNDARY or BOUNDARY→REGION. So REGION→BOUNDARY→REGION is **two steps** | `hop` |
| **structurally pure-gDNA object** (**G1 object**) | a slot at which no RNA strand is admissible, so its composition is CERTAIN: an intergenic REGION, or an `intergenic\|exon` BOUNDARY. Its gDNA density is directly observed, with nothing to deconvolve. The predicate is `node_geometry.g1_locked` | `anchor` — ⛔ that word had two meanings at once and now has only the one below |
| **the mass pin** | the operator that rescales a message so that `Σ_c ρ_c·E_c = M` at the destination (`messages/head.py`'s `_pin_v` and its scalar twin inside the scan kernel). "Pin" because the function is named for it | `the mass anchor` |
| **counts** | discrete integer fragment counts | — |
| **density** = **abundance** | counts per base. The two words mean the same thing | ⛔ not the simulator's molar `abundance=` field, which is a per-transcript weight |
| **crossing fragment** | ⭐ a **fragment** that spans a BOUNDARY. Legitimate and necessary — `crossing_eff_length` is the opportunity for exactly this — and it stays | ⛔ only the *object* sense is banned: objects are REGIONs and BOUNDARIES, never "crossings" |
| **switched off** | an A/B in which one code path is disabled and the run repeated, to establish that it is the cause | `ablated` |
| **splice-out** | ⭐ owner, 2026-08-05. A message crossing a BOUNDARY **in the direction in which mature RNA departs** through the sj there. The fragments that splice away leave the contiguous population; the residual continues | ⛔ `peel` |
| **splice-in** | ⭐ owner, 2026-08-05. A message crossing a BOUNDARY **in the direction in which mature RNA arrives** through the sj there — the spliced flux joins the destination's population | ⛔ `graft` |

⛔⛔ **`donor` / `acceptor` ARE BANNED, AND WHAT BANS THEM IS A MEASUREMENT** (owner ruling, carried here
2026-08-14). They are 5′/3′ names and therefore **strand-dependent**, while every structure in this tree is
stored in GENOMIC order. Measured on the suite index: `src < dst` holds on all **13,482** sj including all
**6,527** minus-strand ones — the code is CORRECT — so an identifier such as `donor_cut` holds the
**ACCEPTOR** for **48.4 %** of sj. ⛔ That is worse than a vague name: **the name lies about correct code**,
and a reader who trusts it writes a sign error. Say `boundary_left` / `boundary_right`, or `_lo` / `_hi`.

⚠ **The index's structural flag bits are the one place the words survive, and they are the hazard itself
rather than an exemption.** `FLAG_DONOR_s` marks the genomic-LOW end of an `s`-strand intron on **both**
strands, so on `−` it sits at the transcript's biological ACCEPTOR. `region_geometry.py` therefore states
this at every consumer and names its arrays `_lo` / `_hi`; `test_splice_flux_reframe` gates it on a
`−`-strand sj specifically. ⚠ `TESTING.md` §0b records the same trap from the toy side, where the bits
decide the sign of a whole scenario.

⭐⭐ **`splice-out` / `splice-in` are DIRECTIONAL, and that is the whole reason for the rename.** "Peel"
and "graft" named two operators; the same BOUNDARY is a splice-out for a message travelling one way and a
splice-in for a message travelling the other, so the pair names one thing seen from two sides. ⚠ The old
words are everywhere in `src/` (`_peel_share`, `graft_frame_logvar`, `graft_premise_logvar`, `_gr`) and
are converged as each area is touched, per this section's standing policy — not in one sweep.

⭐⭐⭐ **AND THE DERIVATION SHARPENED THAT INTO A RULING (2026-08-05, `EQUATIONS.md` §3.6c).** The two words
name the SEMANTICS of a step and **not** its arithmetic. An BOUNDARY presents one total to its genomic-LOW
neighbour and another to its genomic-HIGH one, and a hop between adjacent slots always uses the low slot's
HIGH-flank total against the high slot's LOW-flank total — *whichever* of them is the source. So the
message direction decides only whether the step is called a splice-out or a splice-in; it never decides
which number is used. ⛔ **A predicate on the message direction is therefore the wrong shape for anything in
this family, and one on the SIDE is the right one.** ⚠ The corollary that matters when writing the plumbing:
this is still not expressible as one array per pass, because within a single forward pass a BOUNDARY presents
its LOW-flank total as a destination and its HIGH-flank total as the very next hop's source.

⚠ **`docs/TESTING.md` §0b carries the counts/density half of this table** for readers who arrive there
first; this section is the canonical one.

⚠⚠ **The banned words are still widespread in code written before this ruling** — ~1,500 occurrences of
`boundary` / `boundary` / `crossing` / `anchor` across `src/`, `tests/` and `scripts/`, almost all in comments and
docstrings. They are being converged as each area is touched, not in one sweep: a blanket regex over that
much prose mis-fires (`crossing_eff_length` and "crossing fragment" are correct; the fragment-length
"anchor" is a third, unrelated meaning), and a 1,500-line mechanical diff is not reviewable.

---

## 0b. ⭐⭐⭐ THE 0.8.0 RELEASE SCOPE — owner ruling, 2026-08-14

⛔ **This is the organising frame for the whole phase, and it is settled.** The shipped version is
`0.7.1` (`pyproject.toml`); the target is **0.8.0**. `ROADMAP.md` ranks the work *inside* this scope;
what is in it and what is out of it is decided here.

### Three strata are the optimisation target, and the fourth is DEFERRED

The panel crosses **strandedness × capture**, and the four cells are not equally healthy.

> ⭐ **IN SCOPE FOR OPTIMISATION — three strata:** unstranded × capture-OFF, stranded × capture-OFF,
> stranded × capture-ON.
>
> ⛔ **DEFERRED — unstranded × capture-ON.** It is not a development target until the other three are
> fully optimised.

⛔⛔ **DEFERRED IS NOT DROPPED, and that distinction carries the ruling.** The deferred stratum **stays in
every benchmark and every measurement, and must keep being reported**: a panel that cannot see the cell
the tool is worst at cannot tell a real win from a re-labelling, and it is the one cell where a mechanism
that helps everywhere else can hide a catastrophe. ⭐ If it improves as a side effect of work on the other
three, that is a free win — it is never the *justification* for a change. ⛔ It follows that every score
is read **PER STRATUM and never pooled**: a panel total that averages the deferred cell into the three
in-scope ones states neither quantity.

⭐ **Why this shape.** Measured 2026-08-13/14 on the rebuilt 16-condition ladder (noise floor
0.996–1.013): the error is CONCENTRATED — unstranded × capture-ON carries **64.5 % of transcript error
and 90 % of gene-level error**, and the three in-scope strata carry the rest. ⛔ And that cell is not a
gradient anyone can descend yet: on it the tool emits a near-ZERO gDNA fraction *regardless of truth*
(exon `f_g` 0.040 / 0.0016 / 0.0021 at `g05` / `g50` / `g98` against truths 0.054 / 0.518 / 0.982), so it
looks acceptable at low gDNA **by coincidence**. A stratum whose answer barely responds to the truth is
one an optimiser overfits rather than fixes.

### ⛔⛔ The LENGTH CHANNEL is deferred until after 0.8.0

> ⛔ **The fragment-length / "length likelihood" channel as a CALIBRATION COMPOSITION channel is future
> work. 0.8.0 ships WITHOUT it. Do not propose it, do not list it as a candidate, do not rank it.**

⚠ **It does not exist in `src/`.** It was A/B'd once, on 2026-08-10, and never shipped — so this is a
scope ruling and not a code removal, and any doc presenting it as "the next thing to try" is stale.

⛔ **`length_likelihood` in `src/rigel/second_pass.py` is a DIFFERENT thing and is NOT affected** by this
ruling: it is the per-fragment second-pass assignment factor of §4 below, the `f(L)` in
`score = ρ × f(L) × s`. Same words, other axis.

⚠ The two rulings interlock: the optimisation panel gives both origins the SAME length distribution (the
ladder ruling below), so it is structurally incapable of pricing a length composition channel
(TRAPS: equal-lengths-carry-no-composition). The panel and the scope agree.

### CALIBRATION is the focus, and the metric is calibration against ORACLE CALIBRATION

> ⭐ **The focus of development is CALIBRATION.** The metric is the **calibration result scored against
> ORACLE calibration** — not the end-to-end transcript number.

⭐ The end-to-end number stays a thermometer rather than a target (`SUCCESS.md` says why): it moves for
reasons that have nothing to do with calibration — the EM's own evidence, the fragment-length scorer,
depth — so a change to calibration must be scored where it acts. ⛔ The oracle is the origin-split truth,
and the instruments that already own that comparison are `solvability_audit.py`, `pass0_vs_oracle.py` and
`prior_vs_oracle.py`; a calibration change is not measured until one of them has spoken.

⭐ **Scenarios must be CACHED so calibration can be re-run extremely efficiently.** A scan is minutes and
a re-solve is seconds, and every mechanism this project has found came out of the fast loop —
`scripts/design/build_scan_cache.py`, `scripts/design/flgap_study_cache.py`, and the `cache` stage of
`scripts/sim/panel.py`, which builds the oracle cache beside the scan one.

### ⭐⭐⭐ WHY THE LADDER GIVES gDNA AND RNA EQUAL FRAGMENT LENGTHS

> ⭐⭐ **THE EM ALREADY USES THE FRAGMENT-LENGTH DISTRIBUTION.** A large gDNA-vs-RNA length difference lets
> the EM assign fragments on LENGTH ALONE — **bypassing calibration entirely and MASKING bugs in it**.
> Equal lengths FORCE the calibration phase to be exercised.

⭐ It is structural rather than a preference: `scoring.cpp` carries **two** fragment-length lookup tables,
`fl_log_prob_` for RNA and `gdna_fl_log_prob_` for gDNA, so where the two distributions are far apart a
fragment's own length separates the components one fragment at a time and it barely matters what
calibration supplied. The panel would then read green with calibration broken.

⚠ **This reason is STRONGER than the one previously written down**, which was "the length channel is
neutralised, so residual error is attributable to density and strand". That is true and it says what the
panel MEASURES; it does not say what the panel PREVENTS. Where a doc still carries only the weaker form
it is incomplete rather than wrong.

⭐ **What calibration is left with is exactly STRAND and DENSITY** — plus belief propagation across
objects, which is currently OFF (`SilentPolicy`, §6.1 below) — and on an unstranded library it is density
alone, since the strand λ-term is exactly 0 at κ = ½ (§7 below). That is the substrate the three in-scope
strata are optimised on.

⛔ **"Equal" is a CONFIGURATION and never a guarantee**: identical parameters still leave a small realised
gap, because an mRNA fragment must fit inside its transcript and gDNA need not. `TESTING.md` carries the
measured residual and what it was priced at; ⛔ score the length axis against the simulator's own truth
table, never against a nominal parameter.

---

## 1. The shape of the tool

Three stages, `pipeline.py`:

1. **Scan** (`scan_and_buffer`) — a C++ htslib single-pass reader. Resolves fragments against the index,
   trains the strand and fragment-length models from unique mappers, buffers fragments into a columnar
   `FragmentBuffer`, and deposits per-object tallies into the C++ **accumulator** → `AccumulatorPayload`.
2. **Second pass** (`_drain_side_buffer`) — drains the fragments pass 1 held, then
   **Calibrate** (`calibration.calibrate`) — deconvolves each object's unspliced mass into gDNA vs RNA by
   a belief-propagation sweep over the region chain, and fits the library-level parameters.
3. **Quantify** (`quant_from_buffer`) — scores fragments, builds loci by connected components, runs a
   per-locus EM with `n_transcripts + 1` components. The calibration prior enters as **two per-locus
   Dirichlet scalars**, never per-transcript.

⚠ **Calibration cost is depth-INDEPENDENT** — every region in the index is solved regardless of read depth,
so a 971 k-fragment targeted BAM pays full genome-scale cost. On one real run: index load ~7 s, BAM scan
~2 s, calibration ~66 s, per-locus EM ~24 s. **The scan is ~2 % of runtime — accumulator work is
essentially free and calibration is the whole budget.**

---

## 2. The index

`INDEX_FORMAT_VERSION 8`, shipped as `nodes.feather` + `edges.feather`, built and checked by
`calibration/splice_graph.py`.

> The genome is a graph. **Regions** are intervals; **boundaries** connect them. A fragment is a **path**.
> Regions count fragments *contained*; boundaries count fragments *crossing* (a 0-bp boundary, no width).

- **Regions** tile each reference, region bound at **every exon endpoint** of every non-synthetic transcript, with
  **no merging** (`EQUATIONS.md` §1.7). 1 bp regions are legal and common — nothing may assume length > 1.
- **Two boundary kinds.** *Contiguous* boundaries are the boundary between genomically adjacent regions: bidirectional,
  carrying gDNA + RNA, endpoints **implicit** (boundary `i` sits between region `i` and `i+1`). *Junction* boundaries
  are directed `src → dst` in GENOMIC order (§0 — never donor→acceptor), **pure mature RNA** by
  construction, need explicit `(src, dst, strand)`,
  and carry no unspliced channel and no structural flags.
- A splice jump deposits on its **sj boundary only** — never on the contiguous boundaries it splices over.
- Boundaries always run `src < dst`, so **genomic order is a topological order and there is no graph traversal
  anywhere.** The graph is a DAG but **not** a polytree (every sj boundary closes one undirected loop),
  so sj boundaries must be *factors on their endpoint regions*, never message channels (TRAPS: splicing-makes-the-graph-cyclic).
- **8 structural flag bits per contiguous boundary**: TSS / TES / DONOR / ACCEPTOR × {+,−}, not mutually
  exclusive. Carry the raw bits to the consumer; do not pre-derive predicates in the plumbing.
- Validated by invariants I1–I13, two of which re-derive the answer by a **different algorithm**
  (TRAPS: self-checking-validator).
- `manifest.json` records the sources, their sha256 and the build flags. The build is deterministic.

⛔ **Never quote a stored census.** Region and boundary counts are properties of an *annotation*, not of the
tool: run `scripts/design/index_census.py`.

⚠ **`reach` is covered by no existing hash.** A rebuild moved 38 % of human reaches with both
`partition_hash` and `graph_hash` byte-identical, so any calibration-facing cache needs a third key
(TRAPS: a-hash-that-misses-its-artifact).

---

## 3. The accumulator

`tests/native/_accumulator_reference.py` **is the executable specification.** The C++ is gated on
byte-identity to it; where it and a document disagree, the reference wins.

### 3.1 What every object stores

⭐⭐ **A CHANNEL IS STORED WHERE A NAMED CONSUMER READS IT, AND NOWHERE ELSE** (owner-set, 2026-08-08).
The populations therefore do NOT all carry the same channels, and that asymmetry is the design:

    region_contained    count  inv_opportunity_sum
    boundary_unspliced  count  inv_length_sum       mass
    boundary_spliced    count                       mass      certified RNA — nothing deconvolves it
    sj                  count  inv_length_sum       mass[2]   inv_length_sum is LIVE in second_pass

| | | |
|---|---|---|
| `count` | `Σ 1` | statistical power — a count is a count |
| `inv_length_sum` / `inv_opportunity_sum` | `Σ 1/placements` (float64) | an exact model-free density **at a boundary**, and *not* one at a region — which is why it is not called `density` |
| `mass` | `Σ (slice/L)/n_bounds` (float64) | ⭐ the CONSERVED fragment count — sums to **one per fragment**, where `count` is `+1` on each of `max(K,1)` objects. ⭐⭐ A SJ BOUNDARY is a boundary exactly like a contiguous one, so a spliced fragment shares its one unit across every object it crosses, sj included. `EQUATIONS.md` §3b |

⛔ **Six banks were REMOVED on that rule** (2026-08-08): three `region_spanning_*`, the two
spliced-boundary length moments, and `sj_length_sum`. Structs went `Region` 80 → **24 B**, `Boundary`
80 → **48 B**, `SpliceJunction` 40 → **16 B**. ⚠ Deleting `region_spanning` has one structural consequence
worth knowing: **no spliced fragment touches the region axis at all** now — a spliced fragment can never be
*contained*, because both endpoints of an annotated intron are region bounds.

⛔⛔ **AND TWO MORE WENT ON 2026-08-13 — `region_contained_length_sum` and `boundary_unspliced_length_sum`
— BECAUSE THEIR STATED JUSTIFICATION WAS FALSE, WHICH IS WHAT MADE THE DELETION CLEAN RATHER THAN LOSSY.**
This section used to carry a `length_sum  Σ L` row claiming it *"carries the only information about the
gDNA/RNA split when the two components share a mean length"*. It does not: at `mu_g = mu_r = mu` the third
row is `(mu, mu)`, proportional to `(1, 1)` exactly like the other two, so the 3×2 system is still rank
one and `length_sum` removes nothing. It is an independent tilt **only when the means already differ** —
precisely when the first two rows already suffice (`TRAPS: equal-lengths-carry-no-composition`).
⚠ **The claim survived because nothing consumed it**: the banks reached `PopulationView.length_sum` and
stopped, so no test could disagree — a justification with no consumer is a claim with no gate. ⭐ The
retraction is written where the claim was, in `scan_payload.py`'s module docstring, not merely deleted.

⛔ **THE COUNTS KEEP BOTH GENOME-STRAND COLUMNS; THE INV-LENGTH SUMS AND THE BOUNDARY MASSES KEEP ONE —
AND `sj_mass` IS THE ONE REVERSAL OF THAT, RECORDED BELOW.** Which strand
a read aligned to says nothing about whether the molecule was gDNA or RNA, and every consumer summed the
two columns — so a strand axis there is half the bank wasted, and worse, a claim that the split is
meaningful. The counts keep both because the strand model is a Beta-Binomial over them, per strand.

⛔⛔ **AND `sj_count` KEEPS BOTH FOR A REASON THAT IS NOT THE STRAND MODEL** (owner, 2026-08-08). A
sj is stranded by its genomic splicing MOTIF, so the *fragments'* strand looks redundant, and every
consumer today sums the two. It is retained for **ALIGNER-ARTIFACT DETECTION**: aligners emit
false-positive `N` CIGAR ops from plain genomic DNA, `rigel.splice_blacklist` catches only those
`alignable` has enumerated by coordinate, and the empirical detector is that in a **stranded** library a
real sj inherits the global strand specificity while an artifact deposits onto BOTH strands.
⚠ Unstranded data cannot use it (κ = ½ leaves nothing to deviate from) — a property of the detector, not
a reason to drop the column. ⭐ The discriminating information lives ONLY in the split: a clean sj
and an artifactual one at the same depth carry the same total. Gated by
`test_the_junction_STRAND_SPLIT_IS_RETAINED_FOR_ALIGNER_ARTIFACT_DETECTION`.

⛔⛔⛔ **AND THAT SAME PREMISE REVERSES THE ONE-VALUE MASS RULING — ON THE SJ AXIS ALONE** (owner,
2026-08-12; landed 2026-08-13). `sj_mass` is `mass[2]`, and it is the only mass with a strand column.
⭐⭐ **The reversal is admissible because the PREMISE CHANGED, and the premise is recorded WITH the ruling
so that neither direction is re-litigated.** The original ruling was *"nothing reads a mass per strand"*.
That is now false for sj and only for sj: an **artifactual** splice junction accumulates SYMMETRICALLY on
both strands, exactly as gDNA does, so the strand model the tool already has can detect one — but only if it
is given a per-strand OBSERVABLE, and the count alone is not enough, because a count cannot separate a sj
used by many short fragments from one used by few long ones. ⚠ The second reason is structural: without the
bank, artifact filtering needs TWO passes over the BAM (tally, filter, re-accumulate the mass), which is the
one thing the single-pass architecture exists to avoid.
⛔ **The column is `col` — the SAME genome-strand column the count is deposited at** — so `mass[c]/count[c]`
is a per-strand mean and not a ratio of two different populations. ⚠ Summed over columns it is
byte-*comparable* to the single accumulator it replaced but **not bit-identical**: float addition is not
associative and the per-column deposit order differs (~1e-15 relative).
⛔ **SCOPE IT DELIBERATELY: the reversal is for the SJ axis, because that is where the consumer is.**
Whether the BOUNDARY axis's `unspliced_mass` / `spliced_mass` should also go per-strand is a SEPARATE
question with **no consumer named**, and `TRAPS: one-thing-varied` says do not bundle it — the one-value
ruling still STANDS there.
⭐ **`substrate` folds the strand axis at its own boundary**, so `PopulationView.mass` stays
strand-agnostic and no consumer below it changed at all. The per-strand values are deliberately not
re-exported there: their consumer is artifact filtering, which reads the payload.
⭐ **The memory ledger of the bundle, and it goes DOWN.** `SpliceJunction` 16 → **32 B** for the second
mass, against `Region` 24 → **16 B** and `Boundary` 48 → **40 B** for the two `length_sum` banks deleted
with it — **net ~19 MB lower at genome scale**. `sj_mass` moved `SINGLE_COLUMN_AXES` → `BANK_AXES`, and
since the shape digest is DERIVED from those two tables it moved itself.

⛔ **THE FIXED POINT IS GONE** (owner ruling, 2026-08-11: one numeric convention, and mixing two in one
schema is worse than either). **A COUNT is an integer; a FRACTION is float64.** There is no scale
constant and nothing decodes a bank. Headroom is no longer a question anyone has to answer.
⭐ It is also the more ACCURATE choice, measured against exact rational arithmetic on the
reciprocal-opportunity theorem: the fixed point missed the answer by 7.0e-10 … 2.0e-07 where float64
misses by 5.8e-15 … 2.8e-13 — 1e5–7e5× closer. Memory is unchanged **by this ruling** and still flat —
~85 MB at human scale on the struct sizes of the day, no hash map. ⚠ The 2026-08-13 bundle above took a
further ~19 MB off that; the current sizes are `Region` **16 B**, `Boundary` **40 B**, `SpliceJunction`
**32 B**, and the `static_assert`s beside them are the only place worth reading them from.

⭐⭐ **AND ONE OF THEM IS ALREADY LOAD-BEARING IN A WAY WORTH STATING**: `edge_unspliced_inv_length_sum` is LIVE in `second_pass` (via `pipeline`), and being the one channel whose opportunity and deposit cancel identically — `E[inv_length_sum] = rho` exactly, at a boundary, for ANY length distribution — that rho term is **the only provably fragment-length-gap-robust density estimator in the tree**. `EQUATIONS.md` §3c.

⛔ **THE TALLY IS NOT BIT-REPRODUCIBLE ACROSS WORKER COUNTS, and the owner has signed that off**
(2026-08-11). Integer addition is associative, so every COUNT bank still reproduces exactly — measured
0.000e+00 spread over 8 shipped multi-threaded runs. Float addition is not, so the FRACTION banks are
re-associated by the per-worker merge and wander by ~1e-15; that reaches `posterior_mean` at **1.503e-15**
and the deliverable at ~1e-11, five orders below `EMConfig.convergence_delta`.
⭐ **Tests validate the float banks within a DERIVED tolerance, bracketed from both sides** — a nudge just
past it must fail and one just inside must pass, so it is a tolerance rather than a threshold that never
binds. ⚠ This is a SECOND source of irreproducibility, independent of `EMConfig.seed`
(TRAPS: the-deliverable-is-not-reproducible-by-default); pinning the seed alone is no longer sufficient
for byte-identical output. `TRAPS: integer-channels-reproduce` carries the numbers.

### 3.1b ⭐⭐⭐ WHO OWNS A FRAGMENT — and nothing is ever re-attributed

**Owner ruling, 2026-08-08.** It was already true of the accumulator, the `CalibrationResult` and the
sweep; it is stated here because ONE consumer had drifted from it for a year.

> **A REGION owns the fragments CONTAINED in it. An BOUNDARY owns the fragments that CROSS it. No object's
> mass is ever moved onto another object.**

⭐ **A locus therefore collects BOTH kinds of object**: its REGIONs by genomic overlap, and **its BOUNDARIES are
the boundaries that TOUCH its regions**. Every region contributes both of its boundaries, so a locus of `k` contiguous
regions carries `k + 1` boundaries — **its two outer ones included**, because a fragment crossing a locus's
boundary overlaps the locus and is therefore one of its EM candidates.

⭐ The outer boundaries are unambiguous, and structurally so: a locus is bounded by intergenic sequence and
intergenic regions carry no transcripts, so no two loci contend for a boundary. ⚠ Contention is
*rare rather than impossible* — 20–34 boundaries per flgap condition, carrying ~0.01 % of the mass — so
`priors.contended_edges` REPORTS it and an assertion would have died on real data.

⛔ **What this replaced, and why it must not come back.** `assemble_priors` folded each boundary's mass into
ONE flank region, because `_project_regions_to_loci` divides by `region_size_bp` and so cannot see a 0-bp
object. The fold then needed a second heuristic — an intergenic RE-KEY — to stop a locus's far-left boundary
vanishing into its dropped flank. ⚠ The fold predates the accumulator rewrite: shipped **v0.7.1** stored
the two boundary sides as separate banks and used them directly for RNA (no pooling, no owner) while
sending gDNA through `_pooled_seam_arrays`, which rejoined them and had to pick an owner. Its stated
reason is the CONTRACTION — right there, wrong for the COUNT. One object served two purposes.
`TRAPS: a-fold-grows-a-heuristic`.

⭐ **The prior's target, stated once**: `n_gdna` in `em_solver.cpp:apply_grouped_prior_update` is a soft
count of the gDNA fragments that are candidates in ONE multi-locus, each counted once — a multi-locus
being a connected component of transcripts linked by shared fragments. ⛔ **A first-base count of a
locus's fragments is therefore NOT this quantity**: it drops every fragment that overlaps the locus but
starts outside it.

### 3.1c ⭐⭐⭐ THE PRIOR ARBITRATES THE UNSPLICED POOL, AND ITS STRENGTH HAS NO KNOB

**Owner ruling, 2026-08-09/10.** Two statements, both settled, both measured.

> **A spliced fragment is pure RNA and never receives a gDNA candidate**
> (`em_solver.cpp`: `has_gdna = !is_spliced && isfinite(gdna_ll)`), so it does not compete with gDNA and
> **must not enter the pool-level prior**. An unspliced intergenic fragment has no transcript candidate
> and never becomes an EM unit at all. What the prior arbitrates is exactly the unspliced fragments
> inside transcript bounds.

⭐ **AND THE CROSSING→FRAGMENT CONVERSION IS POPULATION-BLIND, MEASURED AND BOUNDED** (2026-08-10). The
boundary term is converted by `q = mass/count`, measured per boundary on the POOLED population and applied to each
component. gDNA and RNA have different true `q` — not because of fragment length (the equal-length null
carries the LARGEST error) but because gDNA crosses boundaries in long intergenic regions where `q → 1` while RNA
sits among short exons. On the TOTAL prior the region term dilutes it to `Δphi` **+0.00013 … +0.00596**, at
or below calibration's own noise floor, and it is not repairable in production because the driver is
placement rather than length. ⛔ Recorded, not fixed: `TRAPS: a-pooled-conversion-applied-per-component`.
⚠ A DIFFERENT quantity from the `Δphi ≤ 5e-4` below, which is the shipped `f_g` against the unspliced
pool; this one isolates the `q` conversion with a perfect `f_g`.

⛔ Putting spliced RNA into `a_r` would penalise gDNA with fragments it could never have won. ⭐ Measured:
against the population it describes — gDNA units + UNSPLICED RNA units — the shipped claim
`phi = a_g/(a_g+a_r)` is exact to **≤ 5e-4** on all four flgap conditions, drained and undrained. ⚠ An
entry once reported a "+0.07…+0.10 tilt" here; it divided by ALL RNA units
(`TRAPS: score-the-consumers-own-count`). The injection is population-matched too, by algebra:
`n_rna = S_r + U_r`, so `R = n_rna + a_r = S_r + (U_r + a_r)` puts the pseudo-count on the unspliced RNA
and leaves the spliced mass untouched; and `out[i] = raw[i]·(1 + a_r/n_rna)` is a UNIFORM scale, so it
changes no transcript's share of RNA — which is exactly "calibration says nothing about which
transcript".

> **The prior's STRENGTH is exactly one pseudo-fragment per real unspliced fragment, by construction, and
> there is deliberately NO knob** (owner, 2026-08-10: a strength parameter is another tunable and we are
> not adding one).

⭐ It follows from the conservation identity — `a_g + a_r` IS the locus's conserved unspliced count, which
IS the pool — and is measured at **Σa/Σpool = 0.999–1.000**, mass-weighted median `w` 1.000, with 0–1 loci
of ~1,300 above 2. So the MAP posterior is a 50/50 blend of calibration and the EM's own evidence.
⚠ Consequence, recorded not acted on: where calibration is wrong by half a unit, half the answer is too.

### 3.2 One strand convention

Everything is stored by **genome** strand (`CHANNEL_PLUS` / `CHANNEL_MINUS`). *Sense* / *antisense* is the
transcript-relative notion and is **derived, never stored**. Two strands exist and they are independent:
`align_strand` and `sj_strand`. An inferred splice is `sj_implicit`, never "inferred".

### 3.3 Two components only

**gDNA and RNA.** "RNA is RNA" — no mature/nascent split in the accumulator. Owner ruling.

### 3.4 Fragment length — ONE definition

`L` = genomic span minus region bound introns. The scanner's rival histogram, `FragmentLengthModels` and the
transcript-space definition are **deleted**, and every histogram `build_fl_models` reads comes from the
payload, so a mixed-frame call is unrepresentable.

⭐ **A gap intron is region bound on EVERY fragment**, not only unspliced ones, with the gaps the CIGAR already
explained excluded by **exact `(start, end)` equality** — overlap would let a different nearby intron
answer for one and make `L` too short.

⚠ `FragmentLengthModel` **singular** is the scorer and stays.

### 3.5 The five length pools

Each is pure **by construction**, and purity is what removes the circularity: a model is fitted only from a
population known to be one component, so nothing is estimated from the fragments it will later explain.

| pool | rule | component |
|---|---|---|
| `DNA_INTERGENIC` | contained in exactly one intergenic region | gDNA |
| `DNA_INTRONIC` | contained in exactly one intronic region | gDNA |
| `DNA_INTRON_EXON` | crosses exactly one boundary, flanks {intron, exon} | gDNA |
| `DNA_INTERGENIC_EXON` | crosses exactly one boundary, flanks {intergenic, exon} | gDNA |
| `RNA_SPLICED` | used an annotated sj, splice **observed** | RNA |

⭐ There is deliberately **no pool** for an exonic contained fragment or a multi-line crossing — those are
mixtures, and an impure pool is worse than a missing one.

⭐ **The pool is keyed on DETERMINACY, not provenance**: a fragment enters when exactly **one** hypothesis
survived, however it got there (TRAPS: a-purity-filter-is-a-length-filter).

⭐ **The two exon-crossing pools are gDNA**, because mature RNA never crosses an exon↔intron boundary
(TRAPS: mature-rna-never-crosses-a-boundary). They were filed for two milestones as "not gDNA"; they are.

### 3.6 Each pool is divided by its OWN opportunity

The gDNA model is fitted from **all four** gDNA pools, each divided by its own opportunity and then
combined — `calibration/gdna_opportunity.py`, derivation in `EQUATIONS.md` §4.3. The RNA model is the
sj pool divided by its own — `calibration/junction_opportunity.py`, `EQUATIONS.md` §4.2.

⛔ **The four gDNA pools must never be pooled raw** (TRAPS: opposite-tilts-must-not-pool), and every divisor is a
**probability**, not a count (TRAPS: divide-by-a-probability).

⭐ **The contained pair dominates off capture; the crossing pair dominates under it.** Under hybrid capture
the surviving off-target gDNA sits beside a probe, and a fragment beside a probe *reaches* the exon
boundary — so it stops being contained and becomes crossing. Fitting from the contained pair alone measures
the short half of one population.

⚠ **The second pass's scorer reads the same de-tilted pools calibration reads**, so there is one definition
of each length distribution rather than two.

### 3.7 The deposit weight is 1/opportunity

Not `1/length`. `EQUATIONS.md` §2 has the derivation and the two places where it stops being model-free.

---

## 4. The two-pass structure

**The accumulator arbitrates.** A fragment arrives with its hypothesis *set*; exactly one survivor
deposits, two or more are held WHOLE in the **side buffer**. A fragment's unsequenced mate gap may hold no
intron, one, or several, and *which* cannot be observed — the bases are not there.

**The second pass drains it**, between the scan and calibration:

1. **score** from pass-1 evidence alone — `score = ρ × f(L) × s`, `EQUATIONS.md` §10;
2. **one multinomial draw** per fragment;
3. **re-deposit through the SAME `deposit`**, with the chosen hypothesis alone — a set of size one, so
   arbitration is degenerate and the ordinary rules decide.

⭐ **One tally path.** There is no second deposit implementation and no duplicated crossing logic, so
byte-identity with the executable specification holds for free rather than by argument.

⭐ **It runs before calibration, and that is the structural decision.** Every input the score needs comes
from pass 1, so **calibration runs exactly once, on the complete tally.**

⚠ **The models the SCORER uses are pass one's; the models CALIBRATION uses are the drained tally's.** Fit
once, score once, drain once, stop. The confident set is biased short, so feeding the drained anchor back
into the score would prefer the shorter — more-spliced — path, and that loop can run away.

⭐ **The FRAGMENT is stored, never its consequences.** Object ids are large, derived, and would have to be
kept consistent with the partition; the fragment is small and replays exactly.

⛔ **The side buffer is the one bank whose ORDER is observable** — it is a list, not a sum — so the C++
export sorts it on the record's own content before it crosses the ABI (TRAPS: integer-channels-reproduce).

**After the drain nothing is held**: the bank is empty and the deferred counters are 0, so "the counter and
the fragments it counts are the same population" stays absolute. Pass one's numbers live in `DrainQC`.

### 4.1 Settled sub-decisions

| | |
|---|---|
| ⭐⭐ **a composition may be imputed across a step iff the source SUPPLIED both components AND the two objects measure the same RNA POPULATION** | owner, 2026-08-04: *"is the source of the message measuring the same thing that I am measuring?"* — if yes, attribute the density discrepancy to capture enrichment; if no, you cannot tell enrichment from a population difference. Derivation + the genomic form of the predicate: `EQUATIONS.md` §3.5b. ⛔ **Termini only** — DONOR/ACCEPTOR change the population too, but their flux is *measured* and the graft and the peel route it |
| ⭐ **the population test is written in GENOMIC terms, never TSS/TES** | the strand flips which flank a terminus implicates; `node_geometry.terminus_flank_gain` is the one home, gated on two mirror-image annotations (TRAPS: substitute-the-definitions-first family) |
| ⭐ **the mass pin carries the same licence, plus the structural case** | it fires iff no BELIEF can reach its budget: the composition was supplied, **or** the destination is a structurally pure-gDNA object whose `f_g = 1` is structure and whose `M/E_g` is therefore an observation. `EQUATIONS.md` §3.5c, TRAPS: the-pin-had-a-fixed-point/TRAPS: no-belief-not-no-numbers |
| gDNA's strand term is **0.5** | double-stranded, no sense direction. A fitted mixture marginal was implemented and refuted (`EQUATIONS.md` §5.3) |
| a flat-zero factor is **skipped**, not multiplied | TRAPS: an-all-zero-factor-is-inert |
| the draw is keyed on **queue position**, never content | `EQUATIONS.md` §10 |
| ambiguous assignment stays **integral** | TRAPS: fractional-mass-is-the-problem |

### 4.2 Still open

| | |
|---|---|
| **D-4** — should a density carry its weight of evidence? | a density of 0.01 from 1000 fragments and from 1 are not the same statement |
| the traffic term's **Poisson likelihood** | `ρ` enters as a hard zero; `P(0|λ,E) = e^(−λE)` is not zero (`EQUATIONS.md` §10). Deferred on a measurement, not on a missing piece — the exposure it needs is now available |

---

## 5. nRNA components

Not per-transcript shadows: unique nascent spans keyed by `(ref, strand, start, end)` are shared across
transcripts and materialized as ordinary transcript rows in `index.t_df`, flagged `is_nrna` /
`is_synthetic`.

⚠ On a **non-synthetic** row `is_nrna` means "single-exon, so mature ≡ nascent" — **not** "manufactured
span". **The real-transcript filter is `~is_synthetic`, alone** (TRAPS: nrna-does-not-mean-synthetic).

---

## 6. Code layout

**Python.** Top level: `cli` `pipeline` `config` `index` `scan` `scoring` `buffer` `scan_payload`
`_accumulator` `locus` `locus_partition` `scored_fragments` `estimator` `strand_model` `frag_length_model`
`second_pass` `splice` `splice_blacklist` `native` `gtf` `transcript` `annotate` `stats` `types`.

`calibration/`: `calibrate` (orchestrator) · `splice_graph` (the v8 index) · ⭐ **`sweep` (THE BACKBONE)**
and **`messages/` (the POLICY: `silent` `head` `variance`)** · `node_chain` `node_geometry` `node_init` ·
`substrate` `region_arrays` `signature` · `effective_length` `capture_eff_length` `fl`
`junction_opportunity` `gdna_opportunity` · `strand_likelihood` `gdna_strand` `strand_balance`
`strand_deconv` `strand_summary` · `density_deconv` `density_model` `gdna_landscape` `npmle`
`background_reference` · `simplex` `simplex_logodds` `derive` `run_fill` · `priors` `result` `errors`
`diagnostics` `track`.

**C++** (`src/rigel/native/`, nanobind, C++17, `-O3`, LTO, OpenMP):

| module | source | purpose |
|---|---|---|
| `_bam_impl` | `bam_scanner.cpp`, `calibration/accumulator.cpp` | BAM parsing, fragment grouping, model training, the accumulator |
| `_em_impl` | `em_solver.cpp` | per-locus EM, connected components (Kahan summation, SIMD `fast_exp.h`) |
| `_scoring_impl` | `scoring.cpp` | fragment likelihood scoring (`-ffast-math`, no SIMD) |
| `_resolve_impl` | `resolve.cpp` | fragment→transcript resolution via cgranges |
| `_cgranges_impl` | vendored | interval overlap |


### 6.1 ⭐⭐⭐ THE SOLVER IS A BACKBONE PLUS A POLICY — settled 2026-08-07, gated on BYTE-IDENTITY

The belief-propagation solver was ONE 1,635-line function in which the shape of the solve and every
argument about what a message should say were interleaved. It is now two things, and the split is
**structural rather than stylistic**:

| | | |
|---|---|---|
| `sweep.py` | ⭐ **THE BACKBONE.** Two directional scans, one combine, one ψ solve, one write-back, **five assertions**. | It knows nothing about capture, grafts, reframes, pins or enrichment — ⛔ and `test_sweep_backbone.py` asserts those words appear in none of its IDENTIFIERS, read from the AST rather than grepped, because grepping matches the docstring that says they are absent |
| `messages/silent.py` | ⭐ `SilentPolicy` — sends nothing. **THE DEFAULT**, five boundaries. | A reader who holds `sweep.py` plus this holds the entire working system |
| `messages/head.py` | `HeadPolicy` — every operator the evolved solver carried, each behind a NAMED switch (**17** of them) | So the panel prices them ONE AT A TIME rather than as a block |
| `messages/variance.py` | was `enrichment_frame.py` — the policy's variance toolbox | ⚠ `count_logvar` is also imported by `node_init`; it has ONE home and this is it |

⭐⭐⭐ **HOW IT WAS ACCEPTED, and this is the part that generalises: a restructure is gated, a rewrite is
not.** Two TRAPS: byte-identity-gates, opposite in direction, both passed:

| | must be | result |
|---|---|---|
| `--arm backbone_head` (every switch ON) | **byte-identical to `base`** | ✅ **1,872 / 1,872 scored fields, 72 / 72 rows** |
| `--arm backbone` (`SilentPolicy`) | **byte-identical to `msgfree_all`** | ✅ **1,872 / 1,872 scored fields, 72 / 72 rows** |

and, per-array on one real 70,176-slot chain (`scripts/design/backbone_parity.py`): **421,056 output
elements and 18,245,830 diagnostic elements, zero differences.** ⭐ The second gate is not only plumbing —
one arm runs the whole relay and discards it while the other never runs it at all, so passing it **PROVES
the relay reaches the answer ONLY through ψ's four imputed channels**.

⛔ **The alternative was tried and is why this route was chosen.** A clean rebuild came out **+103 %** and
took two sessions to decompose into a correct derivation, a UNIT ERROR (a log-odds delivered into an
angle's grid) and a STRUCTURAL disconnection (a one-slot-step channel on a bipartite chain reaches 0
REGIONs). One of the three was a typo-class error no amount of derivation review would have caught. **A
refactor gated on byte-identity has exactly zero of that risk: any difference is a bug, full stop.**

#### The interface, and its one contract

```python
relay = policy.prepare(ctx)                  # one working object per sweep
step, publish = relay.scan(backward=False)   # the per-hop kernel; None ⇒ nothing to relay
msg = relay.deliver(left_state, right_state) # -> PsiMessage, from the NEIGHBOURS only
```

⛔⛔ **TRAPS: a-message-from-the-destinations-belief, and the backbone enforces the enforceable half BY CONSTRUCTION.** `deliver` is handed
:class:`NeighbourState`, whose relayed arrays are **already gathered at the source slot** — so a policy
holding one has values FOR THE SOURCE and no way to ask the same array about the destination. TRAPS: a-message-from-the-destinations-belief's nine
costumes were all a message built from the destination's own relayed belief, and none of them is
expressible through that type. ⭐ A gather is exact, so making the rule structural costs no bits.

`StepContext` splits its fields under three headings — **observations** (either end), **geometry /
structure** (either end), **beliefs** (SOURCE-SIDE ONLY) — and the heading is what turns TRAPS: a-message-from-the-destinations-belief from a
discipline into something a reader can check. ⚠ **One belief field is read at the destination by the
shipped policy and it is a named, measured DEBT rather than an oversight**: `belief_fg` reaches the
reframe's frame pair, so the frame at a hop is a function of the destination's belief, at slots carrying
**57–77 % of library mass**.

#### The five assertions, and why they live in the backbone

⭐⭐ **A future policy can be as wrong as it likes and still cannot commit any of these** — each has shipped
at least once.

| the backbone asserts | it would have caught | on HEAD, `g50 ss0.50 capture_on`, pass-0 |
|---|---|---|
| `deliver` sees only the two NEIGHBOUR states | **TRAPS: a-message-from-the-destinations-belief — nine recurrences in nine costumes** | structural — inexpressible, not checked |
| every message mode lies inside its coordinate's own grid | **TRAPS: off-grid-message-mode** — the tilt bug, 74 % of `g00`'s error | ✅ tilt **0 / 4,795** · ⛔ gDNA **15,240 / 50,984 (29.9 %)** · RNA+ 15.0 % · RNA− 14.0 % · λ 0.46 % |
| every delivered share is in `[0, 1]` | the over-unit certified-RNA claim | ⛔ **3,013 / 15,629 (19.3 %)**, and the SUM ⛔ **31,174 / 50,984 (61.1 %)** |
| `\|T\| ≤ 3` | **AXIOM 0**, made executable | ✅ 0 / 70,176, and **9,912** slots reach 3 so it is not vacuous |
| the write-back touches only `solvable` slots | the basis mismatch that made an TRAPS: byte-identity-gate read `max\|Δ\| = 1.0` | ✅ 0 / 32,817 |

⭐⭐⭐ **AND THE ASSERTIONS PAID FOR THEMSELVES ON THEIR FIRST RUN — six of the ten counts above are new
facts about the shipped message layer.** The two that matter most:

* ⛔⛔ **61.1 % of live message packets assert that the three components TOGETHER account for MORE fragments
  than the slot observed.** That is the identity `Σ_c ρ_c·E_c = M` the mass pin exists to restore, and the
  pin is licensed in only two states, so everywhere else the residual is *delivered* rather than fixed. ⭐ It
  is consistent with a number already in the tree — `messages/variance.py` records the over-claim on 52–71 %
  of regions — but nothing had surfaced it as a **checkable invariant**, so nothing could rank it.
* ⛔ **29.9 % of live gDNA level modes are outside ψ's own log-share grid**, whose domain is
  `[log σ(−L), log σ(+L)] = [−10.000045, −4.54e-5]` and **not** `(−inf, 0]`. The cause is the `_EPS = 1e-9`
  floor: `log(1e-9) = −20.723`, **10.72 nats** below the low end. ⚠ **This is not 29.9 % of the error, and
  the difference from TRAPS: off-grid-message-mode is the honest part**: TRAPS: off-grid-message-mode's tilt pinned at the *wrong* corner, whereas a low-side
  share pins at "as little of this component as the grid can express", which at a slot that genuinely holds
  almost none of it is the *right* answer. What is certain is TRAPS: off-grid-message-mode's mechanism — no interior minimum, so
  precision buys a CORNER rather than a location. Whether the corner is right is **unmeasured**.

⚠ **The low-side check is the one the first draft of this section did not have**, and it was found by an
adversarial read rather than by the panel: the gates were already passing byte-identically, so nothing about
the restructure would ever have surfaced it. ⭐ `share_sum_at_most_one` is a SUPERSET of
`share_in_unit_interval` by construction — report both, never add them.

⛔⛔ **AN ASSERTION THE SHIPPED POLICY VIOLATES IS WAIVED WITH ITS MEASUREMENT, NEVER WIDENED.** The two
that fire are recorded in `sweep._KNOWN_VIOLATIONS` with the number beside each, and `test_sweep_backbone.py`
asserts that every waiver carries a written reason. Widening the predicate to fit the defect is how a gate
becomes vacuous (TRAPS: perturb-every-gate/TRAPS: a-gate-that-reconstructs); a waiver keeps the count VISIBLE and rankable.

⭐ **Each assertion also reports how many slots were ELIGIBLE for it**, because TRAPS: could-the-arm-have-fired: an assertion reporting
zero violations where its predicate can never fire is not evidence of anything.

⚠ **One tolerance in the family is derived and not tuned** (TRAPS: no-magic-numbers): the grid-domain check allows
one **grid SPACING** beyond the outermost grid point, because within one spacing the penalty's minimum is
still that boundary CELL — for the tilt, the legitimate answer "all RNA on one strand". It is load-bearing:
the shipped tilt mode overshoots `π/2` by **2 ULP** at 63 of 4,795 live slots, being a convex mean of two
messages that AGREE on `τ = ±1`, and TRAPS: off-grid-message-mode's real overshoot is **57 spacings** — four orders of magnitude
apart, so no threshold is being bought.

---

## 7. Where the error is, structurally

Two measurements that shape every design choice above, both worth keeping because they are about the
*structure* rather than about a particular library:

⚠ **Both decompose the error by OBJECT CLASS.** The decomposition by STRATUM — which cell of
strandedness × capture carries it, and therefore which three strata are the 0.8.0 optimisation target —
is §0b above, and the two are read together: an object class carrying the error is a mechanism, a stratum carrying it
is a scope.

**Regions and boundaries measure different components.** The gDNA/RNA opportunity ratio is **0.25 at a crossing
point** (4× RNA-selective) against **115.7 at a 100 bp region**, 25.5 at 147 bp, 2.5 at 300 bp, 1.19 at
1000 bp. A 200 bp RNA fragment cannot fit in a 147 bp region, so a short region is a good gDNA measurement and
says nothing about RNA. **Carry per-component precision, not one scalar.**

**Most of the error is in objects with no evidence of their own.** Over a 32-condition sweep, objects
carried entirely by neighbours were 54.1 % of mass and **91.9 % of all error** (6.2× the rate); objects
with own evidence 29.6 % / 8.1 %; structurally-locked pure-gDNA objects 16.4 % / 0.0 %. And 80.5 % of
total error was honest under-determination — only 1.9 % was confidently wrong.

⭐ **Re-derived per object against the origin-split oracle** (`pass0_vs_oracle.py`, 4 contaminated
conditions, both axes): the relay carries **93.1 %** of pass-0's error, at 1.0–3.7× its mass share on
every arm, and structurally-locked objects carry **0.0 %** — the same shape, on a different panel, with a
different instrument. Two things the older sweep could not say:

* ⭐ **On an unstranded library the DENSITY model carries the entire own-evidence budget.** At κ = ½ the
  strand λ-term is exactly 0 (`EQUATIONS.md` §5.2), and the intron factory still reaches **100 % of
  intron-node mass — both-stranded as well as single-stranded**. Density, not strand, is what makes an
  unstranded library solvable at all.
* ⛔ **But the factory is wired to `(REGION, INTRON)` and to nothing else**, so on an unstranded library
  the relay-only set is exactly *exon regions plus contiguous boundaries*: 27.6 % + 20.6 % of mass off capture,
  53.3 % + 44.6 % on it. ⚠ **Boundaries are the smaller half** — the earlier framing of this as an
  edge-axis property was wrong.
* ⭐⭐ **100.0 % of the relay-only mass sits on slots that HAVE a count and a gDNA opportunity.**
  Relay-only never meant "no information"; it means no channel is wired to read the information present.
  ⛔ So the coverage of the density peel is a *frame and region-type restriction*, not a limit of what is
  knowable — and `EQUATIONS.md` §3.2 (density is frame-invariant to ~0.36 %) is the argument that ρ_bg
  should transport to the crossing frame at all. What it cannot survive is non-uniform gDNA placement,
  which is precisely what capture creates at exons (TRAPS: capture-is-1000x-on-exons).
⛔⛔ **AND ALL OF THE ABOVE SCORES HONEST IGNORANCE AS ERROR, WHICH IS THE WRONG QUESTION FOR PASS-0.**
The prior-free solve exists to produce a **substrate the gDNA hyperprior is fitted against**, not to be
accurate everywhere; an object with no own evidence reporting `f_g ≈ ½` at *zero precision* is stating a
true fact about itself. Excluding the undetermined population, pass-0's error on the objects it CAN
solve is **0.0456** where the all-objects figure said 0.3150 — 99.5 % of the difference is the
undetermined class. **The measurement that matters is solvable → right / wrong → confidently wrong**,
because a wrong value with a tight variance outvotes correct neighbours and anchors the prior
(`scripts/design/solvability_audit.py`; the retired `node_error_attribution.py` and
`confident_fp_trace.py` had this before the refactor and it was lost).

⚠ On real cfRNA most confident-gDNA regions have **zero** counts (64–94 % across libraries), and genome-wide
80.5 % of regions carry no fragments at all. A density-space estimator floors at `1/E` and discards most of
the evidence.
