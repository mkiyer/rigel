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
| **NODE** | a **contiguous genomic interval** — a region of the partition. Has a length in bp | `region` is tolerated: it is the index's own word for the same thing |
| **EDGE** | a **single genomic position** separating two adjacent NODEs. Zero bp wide | `line` · `seam` · `crossing` (as a noun for an object) |
| **BOUNDARY** | ⭐ an accepted synonym for EDGE, kept because the code and the earlier docs are full of it | — |
| **junction** | ⚠ a **splice** junction — a *different* object on a *different* axis, directed donor→acceptor. Never a synonym for EDGE | — |
| **slot** | one entry of the chain, which alternates NODE, EDGE, NODE, EDGE … A slot is a NODE **or** an EDGE | — |
| **step** | one adjacency move along the chain: NODE→EDGE or EDGE→NODE. So NODE→EDGE→NODE is **two steps** | `hop` |
| **structurally pure-gDNA object** (**G1 object**) | a slot at which no RNA strand is admissible, so its composition is CERTAIN: an intergenic NODE, or an `intergenic\|exon` EDGE. Its gDNA density is directly observed, with nothing to deconvolve. The predicate is `node_geometry.g1_locked` | `anchor` — ⛔ that word had two meanings at once and now has only the one below |
| **the mass pin** | the operator that rescales a message so that `Σ_c ρ_c·E_c = M` at the destination (`bp_solver._pin_v` and its scalar twin inside `_relay`). "Pin" because the function is named for it | `the mass anchor` |
| **counts** | discrete integer fragment counts | — |
| **density** = **abundance** | counts per base. The two words mean the same thing | ⛔ not the simulator's molar `abundance=` field, which is a per-transcript weight |
| **crossing fragment** | ⭐ a **fragment** that spans an EDGE. Legitimate and necessary — `crossing_eff_length` is the opportunity for exactly this — and it stays | ⛔ only the *object* sense is banned: objects are NODEs and EDGEs, never "crossings" |
| **switched off** | an A/B in which one code path is disabled and the run repeated, to establish that it is the cause | `ablated` |
| **splice-out** | ⭐ owner, 2026-08-05. A message crossing an EDGE **in the direction in which mature RNA departs** through the junction there. The fragments that splice away leave the contiguous population; the residual continues | ⛔ `peel` |
| **splice-in** | ⭐ owner, 2026-08-05. A message crossing an EDGE **in the direction in which mature RNA arrives** through the junction there — the spliced flux joins the destination's population | ⛔ `graft` |

⭐⭐ **`splice-out` / `splice-in` are DIRECTIONAL, and that is the whole reason for the rename.** "Peel"
and "graft" named two operators; the same EDGE is a splice-out for a message travelling one way and a
splice-in for a message travelling the other, so the pair names one thing seen from two sides. ⚠ The old
words are everywhere in `src/` (`_peel_share`, `graft_frame_logvar`, `graft_premise_logvar`, `_gr`) and
are converged as each area is touched, per this section's standing policy — not in one sweep.

⭐⭐⭐ **AND THE DERIVATION SHARPENED THAT INTO A RULING (2026-08-05, `EQUATIONS.md` §3.6c).** The two words
name the SEMANTICS of a step and **not** its arithmetic. An EDGE presents one total to its genomic-LOW
neighbour and another to its genomic-HIGH one, and a hop between adjacent slots always uses the low slot's
HIGH-flank total against the high slot's LOW-flank total — *whichever* of them is the source. So the
message direction decides only whether the step is called a splice-out or a splice-in; it never decides
which number is used. ⛔ **A predicate on the message direction is therefore the wrong shape for anything in
this family, and one on the SIDE is the right one.** ⚠ The corollary that matters when writing the plumbing:
this is still not expressible as one array per pass, because within a single forward pass an EDGE presents
its LOW-flank total as a destination and its HIGH-flank total as the very next hop's source.

⚠ **`docs/TESTING.md` §0b carries the counts/density half of this table** for readers who arrive there
first; this section is the canonical one.

⚠⚠ **The banned words are still widespread in code written before this ruling** — ~1,500 occurrences of
`line` / `seam` / `crossing` / `anchor` across `src/`, `tests/` and `scripts/`, almost all in comments and
docstrings. They are being converged as each area is touched, not in one sweep: a blanket regex over that
much prose mis-fires (`crossing_eff_length` and "crossing fragment" are correct; the fragment-length
"anchor" is a third, unrelated meaning), and a 1,500-line mechanical diff is not reviewable.

---

## 1. The shape of the tool

Three stages, `pipeline.py`:

1. **Scan** (`scan_and_buffer`) — a C++ htslib single-pass reader. Resolves fragments against the index,
   trains the strand and fragment-length models from unique mappers, buffers fragments into a columnar
   `FragmentBuffer`, and deposits per-object tallies into the C++ **accumulator** → `AccumulatorPayload`.
2. **Second pass** (`_drain_side_buffer`) — drains the fragments pass 1 held, then
   **Calibrate** (`calibration.calibrate`) — deconvolves each object's unspliced mass into gDNA vs RNA by
   a belief-propagation sweep over the node chain, and fits the library-level parameters.
3. **Quantify** (`quant_from_buffer`) — scores fragments, builds loci by connected components, runs a
   per-locus EM with `n_transcripts + 1` components. The calibration prior enters as **two per-locus
   Dirichlet scalars**, never per-transcript.

⚠ **Calibration cost is depth-INDEPENDENT** — every node in the index is solved regardless of read depth,
so a 971 k-fragment targeted BAM pays full genome-scale cost. On one real run: index load ~7 s, BAM scan
~2 s, calibration ~66 s, per-locus EM ~24 s. **The scan is ~2 % of runtime — accumulator work is
essentially free and calibration is the whole budget.**

---

## 2. The index

`INDEX_FORMAT_VERSION 8`, shipped as `nodes.feather` + `edges.feather`, built and checked by
`calibration/splice_graph.py`.

> The genome is a graph. **Nodes** are intervals; **edges** connect them. A fragment is a **path**.
> Nodes count fragments *contained*; edges count fragments *crossing* (a 0-bp line, no width).

- **Nodes** tile each reference, cut at **every exon endpoint** of every non-synthetic transcript, with
  **no merging** (`EQUATIONS.md` §1.7). 1 bp nodes are legal and common — nothing may assume length > 1.
- **Two edge kinds.** *Contiguous* edges are the seam between genomically adjacent nodes: bidirectional,
  carrying gDNA + RNA, endpoints **implicit** (edge `i` sits between node `i` and `i+1`). *Junction* edges
  are directed donor→acceptor, **pure mature RNA** by construction, need explicit `(src, dst, strand)`,
  and carry no unspliced channel and no structural flags.
- A splice jump deposits on its **junction edge only** — never on the contiguous edges it splices over.
- Edges always run `src < dst`, so **genomic order is a topological order and there is no graph traversal
  anywhere.** The graph is a DAG but **not** a polytree (every junction edge closes one undirected loop),
  so junction edges must be *factors on their endpoint nodes*, never message channels (`TRAPS.md` E5).
- **8 structural flag bits per contiguous edge**: TSS / TES / DONOR / ACCEPTOR × {+,−}, not mutually
  exclusive. Carry the raw bits to the consumer; do not pre-derive predicates in the plumbing.
- Validated by invariants I1–I13, two of which re-derive the answer by a **different algorithm**
  (`TRAPS.md` A1).
- `manifest.json` records the sources, their sha256 and the build flags. The build is deterministic.

⛔ **Never quote a stored census.** Node and edge counts are properties of an *annotation*, not of the
tool: run `scripts/design/index_census.py`.

⚠ **`reach` is covered by no existing hash.** A rebuild moved 38 % of human reaches with both
`partition_hash` and `graph_hash` byte-identical, so any calibration-facing cache needs a third key
(`TRAPS.md` E9).

---

## 3. The accumulator

`tests/native/_accumulator_reference.py` **is the executable specification.** The C++ is gated on
byte-identity to it; where it and a document disagree, the reference wins.

### 3.1 What every object stores

**Three integer sums** per object per channel:

| | | |
|---|---|---|
| `count` | `Σ 1` | statistical power — a count is a count |
| `inv_length_sum` | `Σ round(2³²/placements)` | an exact model-free density **at an edge**, and *not* one at a node — which is why it is not called `density` |
| `length_sum` | `Σ L` | carries the only information about the gDNA/RNA split when the two components share a mean length; the other two carry **zero** there |

Fixed-point headroom is ~800×: with `L ∈ [20,2000]` each `round(2³²/L) ≤ 2.1e8`, so 1e8 fragments sum to
≤ 2.1e16 against a uint64 ceiling of 1.8e19. Memory is flat and small — ~85 MB at human scale (node 24 B,
contiguous edge 48 B, junction edge 24 B), no hash map.

⭐ **Every channel is an integer, and that is what makes the tally reproducible across worker counts**
(`TRAPS.md` E10).

### 3.2 One strand convention

Everything is stored by **genome** strand (`CHANNEL_PLUS` / `CHANNEL_MINUS`). *Sense* / *antisense* is the
transcript-relative notion and is **derived, never stored**. Two strands exist and they are independent:
`align_strand` and `sj_strand`. An inferred splice is `sj_implicit`, never "inferred".

### 3.3 Two components only

**gDNA and RNA.** "RNA is RNA" — no mature/nascent split in the accumulator. Owner ruling.

### 3.4 Fragment length — ONE definition

`L` = genomic span minus cut introns. The scanner's rival histogram, `FragmentLengthModels` and the
transcript-space definition are **deleted**, and every histogram `build_fl_models` reads comes from the
payload, so a mixed-frame call is unrepresentable.

⭐ **A gap intron is cut on EVERY fragment**, not only unspliced ones, with the gaps the CIGAR already
explained excluded by **exact `(start, end)` equality** — overlap would let a different nearby intron
answer for one and make `L` too short.

⚠ `FragmentLengthModel` **singular** is the scorer and stays.

### 3.5 The five length pools

Each is pure **by construction**, and purity is what removes the circularity: a model is fitted only from a
population known to be one component, so nothing is estimated from the fragments it will later explain.

| pool | rule | component |
|---|---|---|
| `DNA_INTERGENIC` | contained in exactly one intergenic node | gDNA |
| `DNA_INTRONIC` | contained in exactly one intronic node | gDNA |
| `DNA_INTRON_EXON` | crosses exactly one line, flanks {intron, exon} | gDNA |
| `DNA_INTERGENIC_EXON` | crosses exactly one line, flanks {intergenic, exon} | gDNA |
| `RNA_SPLICED` | used an annotated junction, splice **observed** | RNA |

⭐ There is deliberately **no pool** for an exonic contained fragment or a multi-line crossing — those are
mixtures, and an impure pool is worse than a missing one.

⭐ **The pool is keyed on DETERMINACY, not provenance**: a fragment enters when exactly **one** hypothesis
survived, however it got there (`TRAPS.md` C1).

⭐ **The two exon-crossing pools are gDNA**, because mature RNA never crosses an exon↔intron seam
(`TRAPS.md` F9). They were filed for two milestones as "not gDNA"; they are.

### 3.6 Each pool is divided by its OWN opportunity

The gDNA model is fitted from **all four** gDNA pools, each divided by its own opportunity and then
combined — `calibration/gdna_opportunity.py`, derivation in `EQUATIONS.md` §4.3. The RNA model is the
junction pool divided by its own — `calibration/junction_opportunity.py`, `EQUATIONS.md` §4.2.

⛔ **The four gDNA pools must never be pooled raw** (`TRAPS.md` C4), and every divisor is a
**probability**, not a count (`TRAPS.md` C3).

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
export sorts it on the record's own content before it crosses the ABI (`TRAPS.md` E10).

**After the drain nothing is held**: the bank is empty and the deferred counters are 0, so "the counter and
the fragments it counts are the same population" stays absolute. Pass one's numbers live in `DrainQC`.

### 4.1 Settled sub-decisions

| | |
|---|---|
| ⭐⭐ **a composition may be imputed across a step iff the source SUPPLIED both components AND the two objects measure the same RNA POPULATION** | owner, 2026-08-04: *"is the source of the message measuring the same thing that I am measuring?"* — if yes, attribute the density discrepancy to capture enrichment; if no, you cannot tell enrichment from a population difference. Derivation + the genomic form of the predicate: `EQUATIONS.md` §3.5b. ⛔ **Termini only** — DONOR/ACCEPTOR change the population too, but their flux is *measured* and the graft and the peel route it |
| ⭐ **the population test is written in GENOMIC terms, never TSS/TES** | the strand flips which flank a terminus implicates; `node_geometry.terminus_flank_gain` is the one home, gated on two mirror-image annotations (`TRAPS.md` D4c family) |
| ⭐ **the mass pin carries the same licence, plus the structural case** | it fires iff no BELIEF can reach its budget: the composition was supplied, **or** the destination is a structurally pure-gDNA object whose `f_g = 1` is structure and whose `M/E_g` is therefore an observation. `EQUATIONS.md` §3.5c, `TRAPS.md` D4d/D4e |
| gDNA's strand term is **0.5** | double-stranded, no sense direction. A fitted mixture marginal was implemented and refuted (`EQUATIONS.md` §5.3) |
| a flat-zero factor is **skipped**, not multiplied | `TRAPS.md` D7 |
| the draw is keyed on **queue position**, never content | `EQUATIONS.md` §10 |
| ambiguous assignment stays **integral** | `TRAPS.md` C5 |

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
span". **The real-transcript filter is `~is_synthetic`, alone** (`TRAPS.md` E6).

---

## 6. Code layout

**Python.** Top level: `cli` `pipeline` `config` `index` `scan` `scoring` `buffer` `scan_payload`
`_accumulator` `locus` `locus_partition` `scored_fragments` `estimator` `strand_model` `frag_length_model`
`second_pass` `splice` `splice_blacklist` `native` `gtf` `transcript` `annotate` `stats` `types`.

`calibration/`: `calibrate` (orchestrator) · `splice_graph` (the v8 index) · `bp_solver` (the solver) ·
`node_chain` `node_geometry` `node_init` · `substrate` `region_arrays` `signature` · `effective_length`
`capture_eff_length` `fl` `junction_opportunity` `gdna_opportunity` · `strand_likelihood` `gdna_strand`
`strand_balance` `strand_deconv` `strand_summary` · `density_deconv` `density_model` `gdna_landscape`
`npmle` `background_reference` · `simplex` `simplex_logodds` `enrichment_frame` `derive` `run_fill` ·
`priors` `result` `errors` `diagnostics` `track`.

**C++** (`src/rigel/native/`, nanobind, C++17, `-O3`, LTO, OpenMP):

| module | source | purpose |
|---|---|---|
| `_bam_impl` | `bam_scanner.cpp`, `calibration/accumulator.cpp` | BAM parsing, fragment grouping, model training, the accumulator |
| `_em_impl` | `em_solver.cpp` | per-locus EM, connected components (Kahan summation, SIMD `fast_exp.h`) |
| `_scoring_impl` | `scoring.cpp` | fragment likelihood scoring (`-ffast-math`, no SIMD) |
| `_resolve_impl` | `resolve.cpp` | fragment→transcript resolution via cgranges |
| `_cgranges_impl` | vendored | interval overlap |

---

## 7. Where the error is, structurally

Two measurements that shape every design choice above, both worth keeping because they are about the
*structure* rather than about a particular library:

**Nodes and edges measure different components.** The gDNA/RNA opportunity ratio is **0.25 at a crossing
point** (4× RNA-selective) against **115.7 at a 100 bp node**, 25.5 at 147 bp, 2.5 at 300 bp, 1.19 at
1000 bp. A 200 bp RNA fragment cannot fit in a 147 bp node, so a short node is a good gDNA measurement and
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
* ⛔ **But the factory is wired to `(NODE, INTRON)` and to nothing else**, so on an unstranded library
  the relay-only set is exactly *exon nodes plus contiguous edges*: 27.6 % + 20.6 % of mass off capture,
  53.3 % + 44.6 % on it. ⚠ **Edges are the smaller half** — the earlier framing of this as an
  edge-axis property was wrong.
* ⭐⭐ **100.0 % of the relay-only mass sits on slots that HAVE a count and a gDNA opportunity.**
  Relay-only never meant "no information"; it means no channel is wired to read the information present.
  ⛔ So the coverage of the density peel is a *frame and region-type restriction*, not a limit of what is
  knowable — and `EQUATIONS.md` §3.2 (density is frame-invariant to ~0.36 %) is the argument that ρ_bg
  should transport to the crossing frame at all. What it cannot survive is non-uniform gDNA placement,
  which is precisely what capture creates at exons (`TRAPS.md` F4).
⛔⛔ **AND ALL OF THE ABOVE SCORES HONEST IGNORANCE AS ERROR, WHICH IS THE WRONG QUESTION FOR PASS-0.**
The prior-free solve exists to produce a **substrate the gDNA hyperprior is fitted against**, not to be
accurate everywhere; an object with no own evidence reporting `f_g ≈ ½` at *zero precision* is stating a
true fact about itself. Excluding the undetermined population, pass-0's error on the objects it CAN
solve is **0.0456** where the all-objects figure said 0.3150 — 99.5 % of the difference is the
undetermined class. **The measurement that matters is solvable → right / wrong → confidently wrong**,
because a wrong value with a tight variance outvotes correct neighbours and anchors the prior
(`scripts/design/solvability_audit.py`; the retired `node_error_attribution.py` and
`confident_fp_trace.py` had this before the refactor and it was lost).

⚠ On real cfRNA most confident-gDNA nodes have **zero** counts (64–94 % across libraries), and genome-wide
80.5 % of nodes carry no fragments at all. A density-space estimator floors at `1/E` and discards most of
the evidence.
