# DESIGN — what is built, and the rulings behind it

**Purpose.** This document records the design of Rigel as it is built and the rulings that fixed it:
what was decided, the one measurement that decided it, and the date. Everything here is settled and is
not re-litigated. §0 is the binding vocabulary; §0b carries the 0.8.0 scope, the equal-fragment-lengths
ruling and the nascent scope ruling; §0c settles that calibration is message passing; §3–§6 describe the
accumulator, the two-pass structure and the solver as shipped; §6b–§6c the message layer and ψ; §7 where
the error sits. What does not belong here: derivations (`EQUATIONS.md`), lessons from mistakes
(`TRAPS.md`, cited by name), open problems and refusals with their numbers (`ISSUES.md`), the ranked
next steps (`ROADMAP.md`), how the panels are built (`TESTING.md`), and history (git). Section numbers
are anchors cited from the source, the tests and the other docs; a gap in the numbering is a deleted
section.

---

## 0. Vocabulary — one word per concept (owner, 2026-08-04)

This table is binding on prose, comments, docstrings and identifiers alike. Every banned word was in use
for a concept that already had a name, and the ambiguity cost a reader real time.

| the term | what it means | the banned synonym |
|---|---|---|
| **REGION** | a **contiguous genomic interval** — a region of the partition. Has a length in bp | `region` is tolerated: it is the index's own word for the same thing |
| **BOUNDARY** | a **single genomic position** separating two adjacent REGIONs. Zero bp wide | `line` · `seam` · `crossing` (as a noun for an object) |
| **edge** | kept, as the sanctioned GRAPH-EDGE sense — a BOUNDARY and/or an SJ — which is what makes `edges_df`, `EDGE_KIND_CONTIGUOUS` and `EDGE_KIND_SJ` correct | never as a synonym for BOUNDARY in prose |
| **SPLICE JUNCTION** (`sj`) | a splice junction — a *different* object on a *different* axis, directed `src → dst` in **GENOMIC** order (`src < dst`, always). Never a synonym for BOUNDARY | — |
| **the sj's LOW end / HIGH end** — `_lo` / `_hi`, or `boundary_left` / `boundary_right` | the two BOUNDARIES a sj is anchored at, named in GENOMIC order, which is the order the code stores them in | `donor` · `acceptor` — see below |
| **slot** | one entry of the chain, which alternates REGION, BOUNDARY, REGION, BOUNDARY … A slot is a REGION **or** a BOUNDARY | — |
| **step** | one adjacency move along the chain: REGION→BOUNDARY or BOUNDARY→REGION. So REGION→BOUNDARY→REGION is **two steps** | `hop` |
| **structurally pure-gDNA object** (**G1 object**) | a slot at which no RNA strand is admissible ACCORDING TO THE ANNOTATION, so its composition is certain *given the annotation*: an intergenic REGION, or an `intergenic\|exon` BOUNDARY. Its gDNA density is directly observed, with nothing to deconvolve. The predicate is `region_geometry.g1_locked`. On real data the certainty is annotation-derived rather than physical — intergenic space carries unannotated transcription the tool does not model (§0b) — which is safe for a POOLED level and is not a licence to treat one such slot as ground truth | `anchor` as an object name. The word survives only for the landscape's pooled intergenic level and the aligner's splice-anchor tolerance |
| **counts** | discrete integer fragment counts | — |
| **density** = **abundance** | counts per base. The two words mean the same thing | not the simulator's molar `abundance=` field, which is a per-transcript weight |
| **crossing fragment** | a **fragment** that spans a BOUNDARY. Legitimate and necessary — `crossing_eff_length` is the opportunity for exactly this — and it stays | only the *object* sense is banned: objects are REGIONs and BOUNDARIES, never "crossings" |
| **switched off** | an A/B in which one code path is disabled and the run repeated, to establish that it is the cause | `ablated` |
| **splice-out** | owner, 2026-08-05. A message crossing a BOUNDARY **in the direction in which mature RNA departs** through the sj there. The fragments that splice away leave the contiguous population; the residual continues | `splice_out` as a synonym for the deconvolution verb |
| **splice-in** | owner, 2026-08-05. A message crossing a BOUNDARY **in the direction in which mature RNA arrives** through the sj there — the spliced flux joins the destination's population | `graft` |

⛔ **`donor` / `acceptor` are banned, and what bans them is a measurement** (owner ruling, 2026-08-14).
They are 5′/3′ names and therefore strand-dependent, while every structure in this tree is stored in
GENOMIC order. Measured on the suite index: `src < dst` holds on all 13,482 sj including all 6,527
minus-strand ones — the code is correct — so an identifier such as `donor_cut` holds the ACCEPTOR for
48.4 % of sj. The name lies about correct code, and a reader who trusts it writes a sign error. Say
`boundary_left` / `boundary_right`, or `_lo` / `_hi`. The index's structural flag bits are the one place
the words survive, and they are the hazard itself: `FLAG_DONOR_s` marks the genomic-LOW end of an
`s`-strand intron on both strands, so on `−` it sits at the transcript's biological ACCEPTOR;
`region_geometry.py`, `transfer_rows.junction_exon_side` and `outside_flank` read the bits in genomic
terms, gated on a `−`-strand sj specifically.

**`splice-out` / `splice-in` are directional, and that is the whole reason for the pair.** The same
BOUNDARY is a splice-out for a message travelling one way and a splice-in for a message travelling the
other. `deconvolve` is a different verb — *"deconvolve the gDNA off, RNA is the residual"* — and a
global replace of it corrupts that sense (`TRAPS: two-masks-one-name`; `rename_census.py --sense` exists
to catch it). **The two words name the semantics of a step and not its arithmetic** (2026-08-05): a
BOUNDARY presents one total to its genomic-LOW neighbour and another to its genomic-HIGH one, and a hop
between adjacent slots always uses the low slot's HIGH-flank total against the high slot's LOW-flank
total, whichever of them is the source. ⛔ A predicate on the message direction is therefore the wrong
shape for anything in this family, and one on the SIDE is the right one.

`docs/TESTING.md` §0b carries the counts/density half of this table for readers who arrive there first;
this section is the canonical one. Re-derive the state of any remaining banned word with
`scripts/design/rename_census.py`; never quote a stored count.

---

## 0b. The 0.8.0 release scope (owner ruling, 2026-08-14)

The shipped version is `0.7.1` (`pyproject.toml`); the target is 0.8.0. `ROADMAP.md` ranks the work
*inside* this scope; what is in it and what is out of it is decided here.

### Three strata are the optimisation target, and the fourth is deferred

| stratum | 0.8.0 |
|---|---|
| unstranded × capture-OFF | **in scope** |
| stranded × capture-OFF | **in scope** |
| stranded × capture-ON | **in scope** |
| unstranded × capture-ON | ⛔ **deferred** |

**Deferred is not dropped, and that distinction carries the ruling.** The deferred stratum stays in every
benchmark and every measurement and must keep being reported: a panel that cannot see the cell the tool
is worst at cannot tell a real win from a re-labelling. If it improves as a side effect of work on the
other three, that is a free win; it is never the justification for a change. Every score is read per
stratum and never pooled.

**Why this shape.** Measured 2026-08-13/14 on the rebuilt 16-condition ladder: unstranded × capture-ON
carried 64.5 % of transcript error and 90 % of gene-level error, and that cell is not a gradient anyone
can descend — the tool emits a near-zero gDNA fraction there regardless of truth (exon `f_g` 0.0016 at
`g50` against a truth of 0.518), so it looks acceptable at low gDNA by coincidence. The algebra behind the
blindness: the gDNA fraction cancels from the strand mean at κ = ½, so an unstranded AMBIG slot has no
channel at all (`solvability_audit.py`).

### The length channel is deferred until after 0.8.0

⛔ The fragment-length likelihood as a CALIBRATION COMPOSITION channel is future work. 0.8.0 ships without
it: do not propose it, list it or rank it. It does not exist in `src/` (A/B'd once, 2026-08-10, never
shipped), so this is a scope ruling, not a code removal. Three other things called "length" are not
affected: layer 2's `fl` / `effective_length` / `capture_eff_length` (the OPPORTUNITY model);
`length_likelihood` in `src/rigel/second_pass.py` (the per-fragment assignment factor of §4); and the fl
PMFs priced by `em_fl_ceiling.py`. The panel gives both origins the same length distribution (below), so
it cannot price a length composition channel (`TRAPS: equal-lengths-carry-no-composition`).

### Calibration is the focus, and the metric is calibration against oracle calibration

The metric is the calibration result scored against ORACLE calibration — not the end-to-end transcript
number, which stays a thermometer (`SUCCESS.md` says why). The oracle is the origin-split truth, and the
instruments that own the comparison are `calibration_vs_oracle.py`, `solvability_audit.py` and
`prior_vs_oracle.py`. Scenarios are cached so calibration re-runs in seconds off a scan that took minutes
(`scripts/design/build_scan_cache.py`; `scripts/sim/panel.py cache` builds the oracle cache beside the
scan one). The pattern any future cache copies is the KEY: the scan manifest plus a content hash of every
source file that PRODUCES what is stored, with the module under test left out — sound only while nothing
that module produces is stored (`TRAPS: a-hash-that-misses-its-artifact`).

### Why the ladder gives gDNA and RNA equal fragment lengths

**The EM already uses the fragment-length distribution.** `scoring.cpp` carries two fragment-length
lookup tables, `fl_log_prob_` for RNA and `gdna_fl_log_prob_` for gDNA, so where the two distributions
are far apart a fragment's own length separates the components one fragment at a time and it barely
matters what calibration supplied: a length gap lets the EM split the origins on LENGTH ALONE, bypassing
calibration and masking bugs in it. Equal lengths force calibration to be exercised
(`TRAPS: a-length-gap-bypasses-calibration`); the fl-gap SIDE panel is what exerts the length mechanism
deliberately (`TESTING.md`). What calibration is left with is exactly STRAND and DENSITY, plus the
messages between objects; on an unstranded library it is density alone, since the strand λ-term is
exactly 0 at κ = ½ (`EQUATIONS.md` §5). "Equal" is a configuration and never a guarantee: an mRNA fragment
must fit inside its transcript and gDNA need not, so score the length axis against the simulator's own
truth table, never against a nominal parameter.

### The nascent scope ruling (owner, 2026-08-22)

Real RNA-seq has nascent RNA at very low levels — sparse and rare. Rigel models it simply — one synthetic
transcript spanning every multi-exon transcript — so the tool is ROBUST to it. That is the whole of
nascent RNA's place in the design.

* The default assumption and behaviour: nascent RNA is absent unless there is an abundance of evidence
  otherwise. Experiments designed to capture nascent RNA explicitly are out of scope for Rigel.
* A number measured at the panel's nascent fragment share (20.2 % capture-OFF, 2.6 % capture-ON) is a
  STRESS reading, not an expected-case reading. Robustness may be judged at stress level; design
  decisions may not be driven from it.
* **State the unit.** A nascent entity spans a whole gene (mean 40,667 bp on the ladder's index) while a
  mature transcript is spliced (mean 1,708 bp), so a molecularly sparse population still supplies a large
  share of fragments: on the rebuilt panel the molar ratio is 0.00895 and the fragment share is 20.2 % —
  the same population described two ways. Say which unit a nascent number is in, every time.
* Licensed by this ruling: a global gDNA background pooled over intergenic + intron slots is safe on real
  data — sparse nascent cannot move a pooled megabase-scale level. ⛔ Not on this panel: on the rebuilt
  ladder the intergenic + intron pool is inflated 1.18× at `g50` and 4.49× at `g05` capture-OFF, while
  intergenic-only is inflated exactly 1.0000× on every condition. With strand-specific data most exons
  solve directly — exon-imputation machinery earns its keep on UNSTRANDED data.
* The symmetric honesty: real intergenic regions carry unannotated transcription the tool does not model
  at all. Neither background pool is "clean" on real data, and both are safe for the same reason.

What this ruling is not: a licence to break AXIOM 0 (unspliced RNA at an intron is still RNA), or a
deletion of robustness — the synthetic nascent entity, the `--nrna` harness arms and the zero controls
all stay. It re-ranks concerns; it does not remove the model.

### What "sparse" means as a model (owner, 2026-08-22; the simulator's `sparse` mode)

Sparsity is a per-gene-span ON/OFF pattern, not a low global level. The three rulings that define the
mode (`sim.whole_genome.apply_sparse_nrna`, gated by `tests/test_whole_genome_sim_config.py`): the unit
of sparsity is the GENE SPAN (the nascent entity), not the transcript — isoforms share a span, so a
per-contributor draw would give a 5-isoform gene `1 − 0.9⁵ = 41 %` chance of nascent at
`on_fraction = 0.1`; the level is LOG-UNIFORM over its range (a linear draw on (1, 1000) puts 90 % of its
mass in the top decade); and the level is drawn INDEPENDENTLY of the mature level, so `nascent > mature`
is a real case the tool must survive — the nascent:mature ratio is a stability parameter, not an
expression one. The fragment share is emergent and the length geometry dominates it (a 23.8× factor):
level ~ logU(1, 100) gives 4.2 % at `on_fraction = 0.10` and 20.2 % at 0.50, the ladder's development
stress level. Price the share from the parameters before simulating.

---

## 0c. Calibration is message passing, and the exon is why (owner, 2026-08-17)

This section exists because the derivation below was re-derived from scratch in several sessions, each
time arriving at the same endpoint. If a design you are considering ends at *"then the exon's gDNA level
has to come from its neighbours"*, you have arrived here again.

### The structural fact

gDNA is directly measurable only where no mature transcript crosses. Three places, and the predicate is
the solver's own `mrna_active` — the same one §6b's boundary-axis ruling turns on: intergenic REGIONs,
intron REGIONs, and the BOUNDARIES against them (`exon|intron`, `exon|intergenic`). At an EXON the
unspliced mass is `gDNA + unspliced RNA`, and that is the quantity being solved for.

> ⛔ **An exon's gDNA level cannot be measured. It can only be imputed from its neighbours, and imputation
> across the chain IS message passing. There is no fourth option.**

This is a statement about the annotation and the deposit rule, not about any estimator, so no better
estimator, extra channel or cleverer prior removes it. It is why the solver is a belief-propagation sweep
over the region chain (§1, stage 2) rather than a per-object fit.

### The object system — what each object contributes

```
  intron REGION  <->  intron|exon BOUNDARY  <->  exon REGION  [ <->  exon|intron BOUNDARY  <->  intron REGION ]
       measures                measures              THE TARGET
```

| object | what it contributes |
|---|---|
| **intron REGION** | deconvolve → **gDNA + unspliced RNA**. A direct gDNA observation, subject to the nascent level |
| **`intron\|exon` BOUNDARY** | its UNSPLICED mass is **gDNA + unspliced RNA**; its **SPLICED** mass is **certified RNA** — gDNA cannot splice, so that arm needs no deconvolution |
| **exon REGION** | **gDNA + RNA — the target**, and unmeasurable alone |

Both flanking objects measure something the exon needs, and each measures a different one: the flank's
gDNA density bounds the exon from BELOW, and the flank's sj FLUX measures the exon's RNA, which bounds it
from ABOVE (§0c.3). The exon is the only slot in that picture with nothing of its own.

### The three escape hatches, all closed (2026-08-17)

Hybrid-capture enrichment is per exon, nonlinear and arbitrary — it depends on which probes the panel
contains. Measured at `g98 ss0.99`, the four gDNA-measuring rungs (intergenic, intron, `exon|intergenic`,
`exon|intron`) are one number off capture and span 122× under it, ordered by probe proximity, the two
boundary rungs only partially enriched (under-reading a true exon by 2.6–3.6×). So: ① a pooled scalar
reference is 3.90× worse than `base` on stranded × capture-ON; ② a local "nearest measured rung" is
1.27–1.50× worse there (and better off capture — a good estimator whose premise capture destroys);
③ `capture_eff_length` takes a solved `CalibrationResult` as its first argument and cannot be the thing
that produces the solve.

### 0c.0 A hop carries a LEVEL or a COMPOSITION, and the two invariances are complementary

The rule the message layer is built around (owner, 2026-08-18/19). One reason sits under every hop type:

| currency | what it is | INVARIANT to | DESTROYED by |
|---|---|---|---|
| **LEVEL** `rho_g` | gDNA fragments per placement | **POPULATION** changes — TSS, TES, strand flips, splicing. gDNA is genomically continuous and knows nothing about transcripts | **ENRICHMENT** changes — a probe edge between the two objects |
| **COMPOSITION** `f_g` | the gDNA share of the crossing population | **ENRICHMENT** changes — a probe enriches everything overlapping it, both components alike, so the ratio survives | **POPULATION** changes — the denominator is a different set of molecules |

They are complementary, so at every hop at least one is intact: use that one. The licence question *is*
the currency question, asked once per hop type. The transfer policy (§6b.12) carries a composition only
across a licensed face and a LEVEL everywhere else. **The population of a message is direction-dependent**
(owner, 2026-08-18: *"what crosses INTO this region?"*): BOUNDARY → EXON, the SPLICE IN, includes the
spliced fragments, which splice in to this exon; EXON → BOUNDARY, the SPLICE OUT, excludes them. One rule
— the message's population is whatever physically enters the destination — evaluated in two directions.

**Measured on the whole panel, 2026-08-19** (the hop-currency instrument, since retired; `transport_dispersion.py`
is the surviving decomposition): a terminus
— a TSS/TES, at a gene edge or inside another transcript — is a POPULATION change and carries a LEVEL; a
splice site into an exon carries the SPLICE-IN COMPOSITION; an intron into its own boundary carries a
COMPOSITION (exact to the fragment where a LEVEL is off by 78–98 % under capture); a hop OUT of an exon
carries a LEVEL on every arm. ⛔ So a hop type is `object class × {splice site, terminus, both}` read off
the boundary flags — the class alone conflates a TSS/TES inside a gene with a splice site, and the two
have opposite currencies (`TRAPS: an-object-class-does-not-see-a-terminus`). At a TERMINUS INTO AN EXON
UNDER CAPTURE neither currency survives — the residual §0c.3's spike-and-slab exists for. Three things
are settled and are not re-asked: forward-backward stays; the three populations `{gDNA, RNA+, RNA−}` are
carried natively (AXIOM 0); and messages DO flow into gDNA-measuring objects — an exon sends into its
`exon|intron` boundary and the boundary does the splice-out arithmetic.

### 0c.0b The two strategies are one continuum (owner, 2026-08-20)

§0c.0's table is a statement about invariances, not an instruction to pick one of two mechanisms.
Abundance transfer (the enrichment is 1) and composition transfer (the enrichment is exactly what the
totals report) are two point hypotheses about one unknown, the log enrichment; where a policy needs the
point between them it is fitted per hop from the data, never set. No shipped mechanism implements a
strategy switch, and none may.

### 0c.0c An imputation is a weak predictor by definition, and its precision must say so

Owner ruling, 2026-08-20: message propagation is an imputation, not a measurement; the strand model is a
measurement. The consequence is a design constraint: a regression on strand-specific data is a PRECISION
defect until proven otherwise. ⛔ And the goal is not a switch. Strand specificity is a spectrum, so the
tool must not have cliff behaviour at a threshold, and the answer is never "turn message propagation off
when stranded"; it is honest precision. The ladder simulates specificity as a binary (0.50 / 0.99) as a
convenience; a mechanism that works only because the panel is coarse is not a mechanism. The measurement
that makes this checkable is a split — score destinations that HAVE their own composition evidence apart
from those that do not (`TRAPS: an-imputation-must-cost-something-every-hop`).

### 0c.0d What pass-0 is for — the circular bootstrap

Owner, 2026-08-20. Calibration is circular from the beginning: the tool does not know whether the library
is capture-enriched nor how much gDNA it holds. The goal of pass-0 is to solve or impute ENOUGH exons —
not all of them — that the gDNA landscape can be trained on them; once trained, that per-object prior
does the work. Two rulings follow: the data's background enters ψ as a LIKELIHOOD wherever the data allow
it (the density λ-factor, whose precision scales with counts — never as a reference location, §6b.1), and
message propagation's irreducible job is the deep chain — long runs of `exon|exon` boundaries where
imputation is the only information there is, which is why it has to work and has to be weak.

### 0c.0e The completion contract of the message policy (owner, 2026-09-02; fulfilled 2026-09-09)

The policy is finished when every case is handled, and not before. No message may stay nullified and no
boundary or region may be skipped: for every node type and every boundary case the policy passes
whatever CAN be passed, forward and backward, so that every node solves with two messages, each with an
honest precision. Two laws follow. (1) **Multi-hop.** A message travels a chain of `exon|exon`
boundaries, each hop charging its own dampening, never a constant. (2) **gDNA is always conveyed.** Where
composition cannot cross — a terminus, a strand change — the gDNA LEVEL still can, with its own honest
width; the part of the message that survives is passed, never the whole message dropped. The acceptance
criterion: on the test chromosome, the ladder and every panel, each case either improves accuracy or
leaves it stable with minimal harm — the two halves judged apart, never pooled, node-locally at the
destinations beside the whole-library number. That held on 2026-09-09; the default flipped and the tool
shrank around one policy (§6.1). The contract stays as the definition every future message change is
held to.

### 0c.1 The mechanism is built and ships — do not build it again

The hop the derivation above asks for is the transfer policy's SPLICE-IN FACE MAP
(`messages/transfer._splice_faces`, `transfer_rows.face_map_lambda`; §6b.9 and §6b.12): the BOUNDARY's
measured sj flux is a density at the source, which joins the RNA claim entering the destination EXON —
the certified flux caps the claimable gDNA share. Only an EXON receives it; the flux is a measurement
with its own counting width, never an imputation; and it is registered by geometry, never gated on the
strand channel's precision, so it survives exactly the stratum where the strand channel is dead.
`CalibrationConfig.message_propagation = False` installs `SilentPolicy` — the measured floor (§6.1).

### 0c.3 The shape of the reference under capture is spike-and-slab

A ruling on the SHAPE, stated here in full. Write `rho_g,i = rho_0 · eps_i`, with
`rho_0` the un-enriched density and `eps_i ≥ 1` the enrichment at object `i`. A probe panel makes the
distribution of `eps` across exons a SPIKE at `eps = 1` (unprobed) plus a SLAB at high `eps` (probed) —
the physical structure of a capture panel, and the 122× ladder of §0c seen as a marginal. So the reference
on `log rho_g` is a two-component mixture, every part of it set by something already measured or ruled:

| part | what sets it |
|---|---|
| **spike** | `rho_0` and its width, measured exactly from the off-target objects — intergenic + intron REGIONs |
| **slab — lower endpoint** | the monotone bound from the adjacent BOUNDARY: enrichment is monotone in probe proximity, so a flank's gDNA density is a floor for the exon beside it |
| **slab — upper endpoint** | the object's OWN total density — it cannot hold more gDNA than the mass it holds |
| **the mixing weight** | the unprobed fraction. At pass-0 nothing is observed about the probe indicator, so it is the reference's own symmetric Jeffreys exponents (`simplex_logodds._JEFFREYS_REF`) with no observation added — mean and median ½. No new constant (`TRAPS: no-magic-numbers`) |

It degenerates to the shipped form where its extra assumption is inert: capture-OFF, the off-target to
boundary ratio is 0.98, so the slab collapses onto the spike; `g00`, both endpoints are 0 and `f_g → 0`;
`g98` capture-ON, the slab's upper endpoint reads 0.9918 against a truth of 0.9817.

**The two bounds.** Both come from neighbours — §0c's structural fact in arithmetic:
`f_g ≥ rho_flank · E_g / M` (the adjacent `exon|intron` flank) and `f_g ≤ 1 − (rho_r · E_r − S) / M`
(`rho_r` from the adjacent BOUNDARY's sj flux, `S` the certified spliced count, off the conserved identity
`rho_r · E_r = unspliced RNA + S`, `EQUATIONS.md` §3b). ⛔ They are a support constraint, not an estimate
(`TRAPS: an-upper-bound-is-not-an-estimate`). The upper bound is loose where RNA is abundant (0.6039
against 0.0000 at `g00`) and tight where it is scarce (0.9918 against 0.9817 at `g98 ss0.99` ON); the
lower bound is tight at low `f_g` and broken by capture; the slab is bracketed by both.

**It escapes the vertex theorem.** `EQUATIONS.md` §9a proves a simplex vertex unreachable without
evidence for every proper prior with a DENSITY; a spike-and-slab has an ATOM, so its median sits exactly
at the spike whenever the spike carries half the mass. The theorem stands for priors with a density, and
the vertex shortfall it prices (`vertex_ceiling.py`) is a property of the prior FAMILY, not of the data.

**Reconciling with §6b's "RNA is the residual and is never predicted" — by scope.** §6b's ruling is
about the OFF-TARGET case, where gDNA is near-uniform and RNA has no genomic autocorrelation, so a pooled
RNA density is inadmissible. At an EXON UNDER CAPTURE the roles are exchanged: gDNA is the unpredictable
one, while RNA is measured locally by the adjacent BOUNDARY's sj flux — that exon's own neighbour's
count, never a pooled density. Which component is the residual is decided by which one the data
measures at that object, and capture moves that line.

---

## 1. The shape of the tool

Three stages, `pipeline.py`:

1. **Scan** (`scan_and_buffer`) — a C++ htslib single-pass reader. Resolves fragments against the index,
   trains the strand and fragment-length models from unique mappers, buffers fragments into a columnar
   `FragmentBuffer`, and deposits per-object tallies into the C++ accumulator → `AccumulatorPayload`.
2. **Second pass** (`_drain_side_buffer`) — drains the fragments pass 1 held, then **Calibrate**
   (`calibration.calibrate`) — deconvolves each object's unspliced mass into gDNA vs RNA by a
   belief-propagation sweep over the region chain, and fits the library-level parameters.
3. **Quantify** (`quant_from_buffer`) — scores fragments, builds loci by connected components, runs a
   per-locus EM with `n_transcripts + 1` components. The calibration prior enters as two per-locus
   Dirichlet scalars, never per-transcript.

Calibration cost is depth-independent — every region in the index is solved regardless of read depth
(one real run: index load ~7 s, BAM scan ~2 s, calibration ~66 s, EM ~24 s). Calibration is the budget.

---

## 2. The index

`INDEX_FORMAT_VERSION 8`, shipped as `regions.feather` + `edges.feather`, built and checked by
`calibration/splice_graph.py`.

> The genome is a graph. **Regions** are intervals; **boundaries** connect them. A fragment is a **path**.
> Regions count fragments *contained*; boundaries count fragments *crossing* (a 0-bp boundary, no width).

- **Regions** tile each reference, cut at every exon endpoint of every non-synthetic transcript, with no
  merging (`EQUATIONS.md` §1). 1 bp regions are legal and common — nothing may assume length > 1.
- **Two boundary kinds.** *Contiguous* boundaries sit between genomically adjacent regions:
  bidirectional, carrying gDNA + RNA, endpoints implicit (boundary `i` sits between region `i` and
  `i+1`). *Junction* boundaries are directed `src → dst` in GENOMIC order (§0 — never donor→acceptor),
  pure mature RNA by construction, need explicit `(src, dst, strand)`, and carry no unspliced channel and
  no structural flags. A splice jump deposits on its sj boundary only.
- Boundaries always run `src < dst`, so genomic order is a topological order and there is no graph
  traversal anywhere. The graph is a DAG but not a polytree (every sj boundary closes one undirected
  loop), so sj boundaries are *factors on their endpoint regions*, never message channels
  (`TRAPS: splicing-makes-the-graph-cyclic`).
- **8 structural flag bits per contiguous boundary**: TSS / TES / DONOR / ACCEPTOR × {+,−}, not mutually
  exclusive. Carry the raw bits to the consumer; do not pre-derive predicates in the plumbing.
- Validated by thirteen invariants, two of which re-derive the answer by a different algorithm
  (`TRAPS: self-checking-validator`); `manifest.json` records the sources and the build is deterministic.

⛔ Never quote a stored census: region and boundary counts are properties of an annotation, not of the
tool — run `scripts/design/index_census.py`. And `reach` is covered by no existing hash: a rebuild moved
38 % of human reaches with both `partition_hash` and `graph_hash` byte-identical, so any
calibration-facing cache needs a third key (`TRAPS: a-hash-that-misses-its-artifact`).

---

## 3. The accumulator

`tests/native/_accumulator_reference.py` is the executable specification. The C++ is gated on
byte-identity to it; where it and a document disagree, the reference wins.

### 3.1 What every object stores

**A channel is stored where a named consumer reads it, and nowhere else** (owner, 2026-08-08). The
populations therefore do not all carry the same channels, and that asymmetry is the design:

    region_contained    count  inv_opportunity_sum
    region (per path)   start_count[2]  end_count[2]  span_count[2]
    boundary_unspliced  count  inv_length_sum       mass
    boundary_spliced    count                       mass      certified RNA — nothing deconvolves it
    sj                  count  inv_length_sum       mass[2]   inv_length_sum is LIVE in second_pass

| channel | | |
|---|---|---|
| `count` | `Σ 1` | statistical power — a count is a count |
| `inv_length_sum` / `inv_opportunity_sum` | `Σ 1/A(w)` (float64) | two deposit rules under two names, each cancelling its own opportunity on its own support (`E[Σ] = ρ·P(A>0)`, `EQUATIONS.md` §2): the BOUNDARY/sj form `1/(w−1)` has `P(w≥2) = 1` on any real library — an exact model-free density — while the REGION form `1/(ℓ−w+1)` reads `ρ·P(w≤ℓ)`, a density truncated by a per-component pmf functional (`TRAPS: a-cancellation-is-conditional-on-its-support`) |
| `mass` | `Σ (slice/L)/n_bounds` (float64) | the conserved fragment count — sums to one per fragment, where `count` is `+1` on each of `max(K,1)` objects. A SJ BOUNDARY is a boundary exactly like a contiguous one, so a spliced fragment shares its one unit across every object it crosses, sj included (`EQUATIONS.md` §3b) |

**The start/end/span region banks** (owner's taxonomy, 2026-08-21). Every accepted path books its first
covered base (`region_start_count`), its last (`region_end_count`), and every region it strictly spans
(`region_span_count`), each by genome strand. START and END have opportunity `ℓ` for every fragment
length — the composition-free TOTAL at a REGION — and are wall-blind only at one end each, so a consumer
side-selects; SPAN's opportunity is `(w−ℓ−1)₊`, a per-component pmf functional, consumer-gated. The
ledger closes twice over: `ΣS = ΣE = qc.deposited`; per region `contained ≤ min(S, E)`.

**What the banks do not carry, and why.** No spliced fragment touches the region axis: both endpoints of
an annotated intron are region bounds, so a spliced fragment can never be *contained*. No length-sum
bank exists: at `mu_g = mu_r` a `Σ L` row is proportional to the count row and adds no tilt
(`TRAPS: equal-lengths-carry-no-composition`). The counts keep both genome-strand columns (the strand
model is a Beta-Binomial over them); the inv-length sums and the boundary masses keep one, because every
consumer sums the columns. `sj_count` and `sj_mass` keep both, for aligner-artifact detection (owner,
2026-08-08 / 2026-08-12): in a stranded library a real sj inherits the global strand specificity while a
false-positive `N` op from genomic DNA deposits onto both strands, and a count alone cannot separate a sj
used by many short fragments from one used by few long ones. `substrate` folds the strand axis at its own
boundary, so `PopulationView.mass` stays strand-agnostic
(`test_the_junction_STRAND_SPLIT_IS_RETAINED_FOR_ALIGNER_ARTIFACT_DETECTION`).

**One numeric convention** (owner, 2026-08-11): a COUNT is an integer; a FRACTION is float64 — no fixed
point, no scale constant, nothing decodes a bank (measured against exact rational arithmetic, float64 is
1e5–7e5× closer than the fixed point it replaced). Memory is flat (~85 MB at human scale); the
`static_assert`s beside the structs are the only place worth reading the struct sizes from. **The tally
is not bit-reproducible across worker counts, and the owner has signed that off** (2026-08-11): every
COUNT bank reproduces exactly; the FRACTION banks are re-associated by the per-worker merge and wander by
~1e-15, reaching the deliverable at ~1e-11, five orders below `EMConfig.convergence_delta`. Tests
validate the float banks within a derived tolerance, bracketed from both sides
(`TRAPS: integer-channels-reproduce`).

#### 3.1a-i The four forms an abundance can take — the classification rule (2026-08-21)

A census of 240 sites found the tree forms a density four ways, and only two of them are defects:

| form | expectation | verdict |
|---|---|---|
| the REGION contained reciprocal bank, `Σ 1/(ℓ−w+1)` | `ρ·P(w ≤ ℓ)` — TRUNCATED | model-free but biased; 11.6× at a 98 bp exon |
| `count / E_contained(ℓ, pmf)` | `ρ` — unbiased | correct, but the fragment-length pmf enters the divisor, so a distorted pmf distorts the level |
| the BOUNDARY/sj reciprocal banks, `Σ 1/(w−1)` | `ρ` exactly, every library | already model-free and unbiased — leave them alone |
| `mass / E_c` with numerator and divisor COMPONENT-MATCHED | `ρ_c` | not an abundance site at all: a deconvolved gDNA mass over the gDNA opportunity is correct |

**The scope rule that follows (owner, 2026-08-21): a TOTAL abundance is a PRE-solve instrument.** It is
needed in exactly three places — the density model, the total-abundance LANDSCAPE, and ENRICHMENT RATIOS.
Everything post-calibration is out of scope by construction: once beliefs exist, each component has its
own fragment-length distribution and its own opportunity, so a per-component estimate is the right
instrument. `capture_eff_length` and `priors` consume a solved `CalibrationResult`, so their `mass/E_g` is
component-matched and a total does not belong there.

#### 3.1a-ii The wall rule and the side selection — the consumer half (2026-08-21)

`rigel/calibration/total_abundance.py` turns those banks into a per-slot TOTAL. Four things are settled:
(1) a side is exact iff its template distance clears `w_max − 1`, and `w_max` is READ from the support
end of `deposited_lengths` — never a quantile — the exact start opportunity being
`A_start(w | d) = min(ℓ, (d + ℓ − w + 1)₊)`, which equals `ℓ` for every `w` iff the template continues
`w_max − 1` bases past the region's genomic-HIGH bound (the END bank mirrors at the LOW bound); (2) the
distance is the component minimum over the populations AXIOM 0 admits at that slot — gDNA's template is
the contig, RNA on strand `s` only where `free_s`, taking the SPLICED distance where an exon covers the
region; (3) a double-walled slot is honestly not model-free and reads NaN — measured coverage
(`1 − double-walled`, START-mass weighted, the ladder) 94.7 % at capture-OFF, 84.3 % under capture;
(4) the population for the mature distances excludes synthetic spans, and the filter is written — on
every shipped index the exon rows happen to contain no synthetic transcript, so this was once true by
accident (`TRAPS: state-the-population-rule-do-not-inherit-it-from-a-table`).

#### 3.1a-iii Which of the landscape's outputs a consumer may read (2026-08-21, from the grid sweep)

`abundance_landscape.AbundanceLandscape` publishes `rho_0`, `span_R`, `w_slot`, the mode list and an
anchor verdict, and they are not equally trustworthy (measured with a grid sweep of the landscape's bandwidth, 16
conditions, `_N_GRID` swept over a 16× range): (1) `rho_0` and the anchor verdict are consumable —
`rho_0` moves 8–25 % across the whole range and the anchor-consistency verdict holds 12/12 on every
contaminated row at every grid; (2) `span_R` is NOT consumable as it stands — on `g50 ss0.99` capture-OFF
it reads 58 → 77 → 95.6 → 94.7 → 1.9 as the grid refines, because `split_basins` selects the enriched
mode by basin mass and over-resolution fragments the bulk into sub-bumps; the mode COUNT may not be read
at all (`TRAPS: a-mode-count-is-not-a-well-posed-quantity`); (3) the estimand is what makes this
landscape right, not the estimator — a fit on `mass / eff_gdna`, a total over one component's
opportunity model, carries the divisor's per-region spread (offset IQR 0.12 nats off capture, 1.66 under
it, removable by no bandwidth), while the landscape's divisor is a geometry. This module is not the
refused drop-in of `ISSUES: the-truncation-free-region-bank`; that refusal's bar still stands for any
consumer swap.

### 3.1b Who owns a fragment — and nothing is ever re-attributed (owner, 2026-08-08)

> **A REGION owns the fragments CONTAINED in it. A BOUNDARY owns the fragments that CROSS it. No object's
> mass is ever moved onto another object.**

A locus therefore collects both kinds of object: its REGIONs by genomic overlap, and its BOUNDARIES are
the boundaries that TOUCH its regions — a locus of `k` contiguous regions carries `k + 1` boundaries, its
two outer ones included, because a fragment crossing a locus's boundary overlaps the locus and is one of
its EM candidates. Contention for a boundary is rare rather than impossible (~0.01 % of the mass), so
`priors.contended_boundaries` reports it. What this replaced: `assemble_priors` folded each boundary's
mass into one flank region because `_project_regions_to_loci` cannot see a 0-bp object, and the fold then
needed an intergenic re-key (`TRAPS: a-fold-grows-a-heuristic`). The prior's target, stated once:
`n_gdna` in `em_solver.cpp:apply_grouped_prior_update` is a soft count of the gDNA fragments that are
candidates in one multi-locus, each counted once; a first-base count of a locus's fragments is not this
quantity, because it drops every fragment that overlaps the locus but starts outside it.

### 3.1b-ii Compatibility is checked before a fragment is called a chimera (owner, 2026-08-19)

A fragment is a chimera **only if its mates are genomically incompatible**:

    a different reference   OR   not facing inward   OR   implied fragment length > max_frag_length

⛔ Transcript-set disjointness alone is not evidence of a rearrangement: gDNA is genomically contiguous
and routinely spans two transcripts that share nothing. The orientation test needed no new predicate:
`build_fragment` keys blocks by `(ref_id, ref_strand)` with R2's orientation flipped, so
`unique_strands.size() == 1` (`CHIMERA_CIS_STRAND_SAME`) *is* "facing inward". The length is the implied
FRAGMENT LENGTH — outermost start to outermost end — never `min_gap` (`test_resolution.py` gates the
distinction). What it cost while wrong: 4,087 gDNA fragments per condition, every one a crosser
(`TRAPS: a-transcript-predicate-must-not-silently-drop-a-molecule`).

### 3.1c The prior arbitrates the unspliced pool, and its strength has no knob (owner, 2026-08-09/10)

> **A spliced fragment is pure RNA and never receives a gDNA candidate** (`em_solver.cpp`:
> `has_gdna = !is_spliced && isfinite(gdna_ll)`), so it must not enter the pool-level prior. An unspliced
> intergenic fragment has no transcript candidate and never becomes an EM unit. What the prior arbitrates
> is exactly the unspliced fragments inside transcript bounds.

Putting spliced RNA into `a_r` would penalise gDNA with fragments it could never have won; measured
against the population it describes, the shipped claim `phi = a_g/(a_g+a_r)` is exact to ≤ 5e-4. The
injection is population-matched by algebra: `R = n_rna + a_r = S_r + (U_r + a_r)` puts the pseudo-count
on the unspliced RNA, and `out[i] = raw[i]·(1 + a_r/n_rna)` is a uniform scale that changes no
transcript's share of RNA. The crossing→fragment conversion `q = mass/count` is population-blind (gDNA
crosses boundaries in long intergenic regions where `q → 1`); it dilutes to `Δphi` ≤ +0.006 on the total
prior — recorded, not fixed (`TRAPS: a-pooled-conversion-applied-per-component`).

> **The prior's strength is exactly one pseudo-fragment per real unspliced fragment, by construction, and
> there is deliberately no knob** (owner, 2026-08-10).

It follows from the conservation identity — `a_g + a_r` IS the locus's conserved unspliced count — and is
measured at `Σa/Σpool = 0.999–1.000`. The MAP posterior is a 50/50 blend of calibration and the EM's own
evidence; where calibration is wrong by half a unit, half the answer is too.

### 3.2 One strand convention

Everything is stored by **genome** strand (`CHANNEL_PLUS` / `CHANNEL_MINUS`). *Sense* / *antisense* is
the transcript-relative notion and is derived, never stored. Two strands exist and they are independent:
`align_strand` and `sj_strand`. An inferred splice is `sj_implicit`, never "inferred".

### 3.3 Two components only

**gDNA and RNA.** "RNA is RNA" — no mature/nascent split in the accumulator. Owner ruling.

### 3.3a No object class is asserted pure gDNA (owner, 2026-08-29)

"Pure" is a property of the annotation and the sample, not of the genome: pervasive transcription is
real, the intergenic space is whatever the user's GTF leaves over, and most genes are OFF in any one
sample but nobody knows which. So the gDNA strand-overdispersion fit trusts no class: it is the away-half
moment (`EQUATIONS.md` §6) over every genic count- and strand-observable object — intron regions,
`exon|intron` and gene-edge boundaries — with intergenic and AMBIG objects out because they cannot be
oriented. Every purity-based fit was moved by unannotated transcription on real cfRNA and on the
blank-chromosome control; the away-half was not. The simulator carries a supplemental `shadow_gtf` the
index never sees, so this stays testable.

Three further rulings (owner, 2026-08-30). **No seed is worth its pair count**: second-moment pooling
weights a seed by `n(n−1)/2`, and on real data one seed carried 77.8 % of a library's numerator, so the
fit uses inverse-variance weights `w_s = 1/(½ + c_s·V∞(ρ))`; trimming or Winsorizing is refused (it
biases the upper tail of the distribution whose mean is the estimand) and the concentration is reported
instead (`GdnaStrandModel.effective_seeds`). **The `Beta(2,2)` ceiling stays**, and a value at it is a
clamp, not a fit (`clamped_at_ceiling`, `raw_overdispersion`): a genuine intra-class correlation is
depth-invariant, and on every real library the moment rises monotonically with seed depth, so above the
ceiling the Beta-Binomial is absorbing a process it has no parameter for — the one asserted constant
left in the strand module. **No component shrinks toward a constant — they shrink toward each other**:
each component reports the precision of its own fitted estimate and the weaker borrows its deficit from
the better measured one (`EQUATIONS.md` §6); on the capture-ON rows, where hybrid capture depletes the
genomic seeds, the reconciliation lands on the oracle value exactly.

### 3.4 Fragment length — one definition

`L` = genomic span minus cut introns. The scanner's rival histogram, `FragmentLengthModels` and the
transcript-space definition are deleted, and every histogram `build_fl_models` reads comes from the
payload, so a mixed-frame call is unrepresentable. A gap intron is cut on every fragment, not only
unspliced ones, with the gaps the CIGAR already explained excluded by exact `(start, end)` equality.
`FragmentLengthModel` singular is the scorer and stays.

### 3.5 The five length pools

Each is pure by construction, and purity is what removes the circularity: a model is fitted only from a
population known to be one component, so nothing is estimated from the fragments it will later explain.

| pool | rule | component |
|---|---|---|
| `DNA_INTERGENIC` | contained in exactly one intergenic region | gDNA |
| `DNA_INTRONIC` | contained in exactly one intronic region | gDNA |
| `DNA_INTRON_EXON` | crosses exactly one boundary, flanks {intron, exon} | gDNA |
| `DNA_INTERGENIC_EXON` | crosses exactly one boundary, flanks {intergenic, exon} | gDNA |
| `RNA_SPLICED` | used an annotated sj, splice observed | RNA |

There is deliberately no pool for an exonic contained fragment or a multi-boundary crossing — those are
mixtures. The pool is keyed on DETERMINACY, not provenance: a fragment enters when exactly one hypothesis
survived, however it got there (`TRAPS: a-purity-filter-is-a-length-filter`). The two exon-crossing pools
are gDNA, because mature RNA never crosses an exon↔intron boundary
(`TRAPS: mature-rna-never-crosses-a-boundary`).

### 3.6 Each pool is divided by its own opportunity

The gDNA model is fitted from all four gDNA pools, each divided by its own opportunity and then combined
— `calibration/gdna_opportunity.py`, derivation in `EQUATIONS.md` §4. The RNA model is the sj pool divided
by its own — `calibration/sj_opportunity.py`. ⛔ The four gDNA pools must never be pooled raw
(`TRAPS: opposite-tilts-must-not-pool`), and every divisor is a probability, not a count
(`TRAPS: divide-by-a-probability`). The contained pair dominates off capture; the crossing pair dominates
under it, because the surviving off-target gDNA sits beside a probe and reaches the exon boundary. The
second pass's scorer reads the same de-tilted pools calibration reads.

### 3.7 The deposit weight is 1/opportunity

Not `1/length`. `EQUATIONS.md` §2 has the derivation, including the support factor `P(A > 0)` that bounds
where each form is model-free: the BOUNDARY/sj forms unconditionally (`P(w ≥ 2) = 1`), the REGION form
only within `w ≤ ℓ` (`TRAPS: a-cancellation-is-conditional-on-its-support`).

---

## 4. The two-pass structure

**The accumulator arbitrates.** A fragment arrives with its hypothesis *set*; exactly one survivor
deposits, two or more are held whole in the **side buffer**. A fragment's unsequenced mate gap may hold
no intron, one, or several, and *which* cannot be observed. **The second pass drains it**, between the
scan and calibration: (1) **score** from pass-1 evidence alone — `score = ρ × f(L) × s`, `EQUATIONS.md`
§10; (2) **one multinomial draw** per fragment; (3) **re-deposit through the same `deposit`**, with the
chosen hypothesis alone — a set of size one, so arbitration is degenerate and the ordinary rules decide.

* **One tally path.** No second deposit implementation and no duplicated crossing logic, so byte-identity
  with the executable specification holds for free.
* **It runs before calibration, and that is the structural decision.** Every input the score needs comes
  from pass 1, so calibration runs exactly once, on the complete tally. The models the scorer uses are
  pass one's; the models calibration uses are the drained tally's. Fit once, score once, drain once,
  stop — the confident set is biased short, so feeding the drained fit back into the score would prefer
  the shorter, more-spliced path, and that loop can run away.
* **The fragment is stored, never its consequences.** Object ids are large and derived; the fragment is
  small and replays exactly. The side buffer is the one bank whose order is observable — a list, not a
  sum — so the C++ export sorts it on the record's own content before it crosses the ABI
  (`TRAPS: integer-channels-reproduce`).
* **After the drain nothing is held**: the bank is empty and the deferred counters are 0. Pass one's
  numbers live in `DrainQC`.

### 4.1 Settled sub-decisions

| | |
|---|---|
| **a composition may be imputed across a step iff the source supplied both components AND the two objects measure the same RNA population** | owner, 2026-08-04: *"is the source of the message measuring the same thing that I am measuring?"* — if yes, attribute the density discrepancy to capture enrichment; if no, enrichment and a population difference cannot be told apart. Termini only — DONOR/ACCEPTOR change the population too, but their flux is *measured* and the splice-in and splice-out route it |
| **the population test is written in genomic terms, never TSS/TES** | the strand flips which flank a terminus implicates; `transfer_rows.outside_flank` is the one home (TSS+ / TES− bodies extend genomic-right, so the outside is the LEFT flank; TES+ / TSS− the reverse), gated on mirror-image annotations |
| gDNA's strand term is **0.5** | double-stranded, no sense direction. A fitted mixture marginal was implemented and refuted (`EQUATIONS.md` §5) |
| a flat-zero factor is **skipped**, not multiplied | `TRAPS: an-all-zero-factor-is-inert` |
| the draw is keyed on **queue position**, never content | `EQUATIONS.md` §10 |
| ambiguous assignment stays **integral** | `TRAPS: fractional-mass-is-the-problem` |

### 4.3 The frame ruling — instruments measure the DRAINED tally (owner, 2026-08-31)

**Every calibration-facing instrument scores the frame production calibrates: the drained one.** The
scan caches store pass one by design (so the drain can re-run at any seed); an instrument loads, drains
at the production seed (`scan_cache.calibration_inputs`), and puts cached ORACLE PARTITIONS into the same
frame by replaying the whole's already-drawn choices (`second_pass.lift_choices`, via
`OracleTruth.from_cached_parts` — one lift call per exact partitioning, any shared member first in each
list so two partitionings stay consistent). Sum-to-full re-validated on the drained frame is the lift's
own end-to-end identity gate. The pricing that decided it (2026-08-31): the P-vs-O metric moved ≤ 1.8 %
on the in-scope contaminated strata, but the `g00` capture-ON zero control moved −22.8 % and the
certified spliced-channel truth was understated 17–20 % undrained. The lift's origin-attribution
ambiguity is bounded and reported (`OracleTruth.n_ambiguous`; ladder-wide ~0.11 % of the library).

**The exact-zeros gate is frame-aware, and that is a semantics ruling.** "gDNA never splices" is true of
molecules and of pass-one deposits, so in the pass-one frame a spliced deposit in the gdna partition is a
hard refusal. In the drained frame the same deposits are a fact of production's tally — the multinomial
draw genuinely assigns a spliced hypothesis to some true-gDNA held fragments — so the oracle records them
per bank (`OracleTruth.gdna_spliced_leak`) instead of refusing to describe the frame. The leak itself is
an open defect: `ISSUES: drain-contaminates-certified-rna`.

---

## 5. nRNA components

Not per-transcript shadows: unique nascent spans keyed by `(ref, strand, start, end)` are shared across
transcripts and materialized as ordinary transcript rows in `index.t_df`, flagged `is_nrna` /
`is_synthetic`. On a non-synthetic row `is_nrna` means "single-exon, so mature ≡ nascent" — not
"manufactured span". The real-transcript filter is `~is_synthetic`, alone
(`TRAPS: nrna-does-not-mean-synthetic`).

---

## 6. Code layout

**Python.** Top level: `cli` `pipeline` `config` `index` `scan` `scoring` `buffer` `scan_payload`
`scan_cache` `locus` `locus_partition` `scored_fragments` `estimator` `strand_model` `frag_length_model`
`second_pass` `splice` `splice_blacklist` `native` `gtf` `transcript` `annotate` `stats` `types`, plus the
`report/` and `sim/` subpackages. `calibration/`: `calibrate` (orchestrator) · `splice_graph` (the v8
index) · `sweep` (the backbone) and `messages/` (the policy: `silent` · `transfer` + `transfer_rows`) ·
`region_chain` `region_geometry` `region_init` `structural_claims` · `substrate` `region_arrays`
`signature` · `effective_length` `capture_eff_length` `fl` `sj_opportunity` `gdna_opportunity` ·
`strand_likelihood` `gdna_strand` `strand_balance` `strand_summary` · `density_deconv`
`density_model` `landscape` `abundance_landscape` `total_abundance` · `simplex_logodds` `derive` ·
`priors` `result` `errors` `diagnostics` `track` · `_layers` (the layering the imports already had).
Re-derive this list rather than trusting it: `scripts/design/module_census.py` reads it off the AST.

**C++** (`src/rigel/native/`, nanobind, C++17, `-O3`, LTO, OpenMP):

| module | source | purpose |
|---|---|---|
| `_bam_impl` | `bam_scanner.cpp`, `calibration/accumulator.cpp` | BAM parsing, fragment grouping, model training, the accumulator |
| `_em_impl` | `em_solver.cpp` | per-locus EM, connected components (Kahan summation, SIMD `fast_exp.h`) |
| `_scoring_impl` | `scoring.cpp` | fragment likelihood scoring (`-ffast-math`, no SIMD) |
| `_resolve_impl` | `resolve.cpp` | fragment→transcript resolution via cgranges |
| `_cgranges_impl` | vendored | interval overlap |

### 6.1 The solver is a backbone plus a policy (settled 2026-08-07, gated on byte-identity)

| | | |
|---|---|---|
| `sweep.py` | **The backbone.** The self-solve, two directional passes, one ψ solve, one write-back, four assertions | It knows nothing about capture, splices, levels, lanes or enrichment — `test_sweep_backbone.py` asserts those words appear in none of its identifiers, read from the AST |
| `blocks.py` | the block plumbing of the locus solve: `view_fields` (every per-slot array a policy may read), `block_slice` (one block cut out of the chain, its links re-based), `gather` (the blocks' captures as the chain's) | nothing about what a solve or a message is |
| `message_cache.py` | `MessageCache` — the message layer's output shared across the refit sweeps, keyed on a digest of every field of the block's context, the library and the policy (§6b.15) | a field added to the context cannot be left out of the key: the digest iterates the dataclass |
| `messages/silent.py` | `SilentPolicy` — sends nothing. **The measured floor**, what `message_propagation = False` installs | A reader who holds `sweep.py` plus this holds the entire working system |
| `messages/transfer.py` | `TransferPolicy` — **the shipped default** (2026-09-09): every node's own claim, one named builder per message, the two passes and the solve (§6b.4–§6b.14) | `prepare` is a table of contents: a reader finds a message by its builder's name |
| `messages/faces.py` | `Faces` — the composition rules as typed tables over `(destination, side)`, `Faces.apply` the one home of the rule arithmetic, and the three helpers every reader of a face needs (`side_of`, `norm`, `fuse`) | gate: `test_transfer_faces.py` |
| `messages/lanes.py` | `LevelLane` — one class for the three populations' levels — and its two builders, `gdna_lane` (every face left without a composition rule) and `rna_lanes` (one per strand, faces from the flag bits) | gate: `test_transfer_rna_lanes.py` |
| `messages/transfer_rows.py` | the pure row constructors — every map, level, price and coordinate change, each a function of one face's numbers | `count_logvar` is the one home of the counting term; every hop price reads it |
| `messages/__init__.py` | the interface (`Policy`, `Prepared`), the message (`Message`, its lanes, `SILENCE`, `NO_NEIGHBOUR`), what ψ receives (`PsiMessage`) and what a policy may read (`BlockContext`) | every field of `BlockContext` has a reader in the policy or the backbone |

**A restructure is gated, a rewrite is not.** The split out of the one 1,635-line function passed two
`TRAPS: byte-identity-gate` gates of opposite direction and, per array on one real 70,176-slot chain
(`scripts/design/backbone_parity.py`), 421,056 output elements with zero differences. The alternative —
a clean rebuild — came out +103 %; a refactor gated on byte-identity has exactly zero of that risk.

#### The interface, and its one contract

```python
prepared = policy.prepare(ctx)                     # one working object per sweep: every node's OWN claim
receive  = prepared.propagate(backward=False)      # phase 1: the recipient's kernel, or None ⇒ all SILENCE
held[i]  = receive(source, destination)            # ... the BACKBONE runs the pass, in chain order
evidence = prepared.solve(from_left, from_right)   # phase 2, the policy's half -> PsiMessage
```

`TRAPS: a-message-from-the-destinations-belief`, and the backbone enforces the enforceable half by
construction: the kernel is called with two indices and builds the message into the destination from the
SOURCE's claim and what the source holds; the backbone writes `held` and the policy never reaches past its
hop. `BlockContext` splits its fields under three headings — **observations** (either end), **geometry /
structure** (either end), **beliefs** (source-side only). The shipped policy reads `belief_fg` once, at
`prepare`, for the variance freeze of each node's own strand profile — a source-side read by construction.

#### The backbone's assertions, and why they live in the backbone

| the backbone asserts | it would have caught |
|---|---|
| the kernel sees only the two NEIGHBOUR states | `TRAPS: a-message-from-the-destinations-belief` — nine recurrences in nine costumes (structural: inexpressible, not checked) |
| every delivered row is one row per slot on the solve grid, finite | a row array off the grid, or a NaN reaching ψ |
| `\|T\| ≤ 3` | AXIOM 0, made executable (on the ladder's `g50 ss0.50 capture_on`: 0 / 70,176, and 9,912 slots reach 3 so it is not vacuous) |
| the write-back touches only `solvable` slots | the basis mismatch that made a byte-identity gate read `max\|Δ\| = 1.0` |

The transfer policy delivers max-normalised profiles on ψ's own grid (`PsiMessage.lam_rows`,
`cube_rows`), which cannot be off-grid and claim no share. ⛔ An assertion the shipped policy violates is
waived with its measurement, never widened: the waiver table (`sweep._KNOWN_VIOLATIONS`) is empty, and
`test_sweep_backbone.py` asserts that any future entry carries a written reason. Each assertion also
reports how many slots were eligible for it (`TRAPS: could-the-arm-have-fired`).

---

## 6b. The message layer and ψ's reference — the rulings behind them

**RNA is the residual and is never predicted** (owner, 2026-08-16). gDNA is near-uniform over the genome
and measurable before any solve; RNA spans six orders of magnitude with essentially no genomic
autocorrelation, so a pooled RNA density is not a population parameter and pooling splice-junction flux
across the genome is inadmissible. The gDNA an object's own density predicts is `ρ_g·E_g,i / M_i`, and
whatever is left is RNA. This is scoped to the off-target case and §0c.3 carries the scope.

**The boundary axis splits on whether mature RNA can cross, not on whether a splice junction attaches**
(owner, 2026-08-15). A BOUNDARY owns the fragments that cross it contiguously, and mature RNA can do that
only where the template is contiguous exon on both sides. So an `exon|intron` boundary is near-pure gDNA
under sparse nascent, while an `exon|exon` boundary is crossed freely by mature RNA. The predicate is the
solver's own `mrna_active`. A pool defined as "a splice junction attaches here" measured true `f_g`
0.0000 over 955,428 fragments at the zero-gDNA control, because it lumps `exon|exon` in.

**The reference and the gDNA landscape partition the object universe rather than competing.** Where the
annotation determines the answer — intergenic and intron REGIONs, `exon|intron`, `intron|intron`,
gene-edge and opposite-strand `exon|exon` BOUNDARIES, 47.5 % of slots — the reference carries it before
any solve exists. Exons, same-strand `exon|exon` and AMBIG have no structural claim and are the
population the landscape is fitted to serve. The reference is worth one pseudo-fragment, so it is
swamped by any evidence channel — but on an evidence-free object posterior = prior at any depth.

### 6b.1 The reference location — refuted and deleted (owner, 2026-08-24)

ψ has no reference location. The structural (`m = 0.75` at `¬mrna_active`) and measured (ss-intron
background) location terms and their config flags were deleted outright; the surviving form is ψ's
symmetric Jeffreys Beta reference, which asserts nothing, and the location concept may not come back as a
flag, a constant or a "weakened" variant. The one measurement: where the strand channel carried no
information (κ = ½, or near a vertex) the location was the entire answer at any depth — 0.7471 at every N
from 10 to 10⁶ — and at zero refits it decided 100 % of the claimed-boundary error. Background
information enters as LIKELIHOOD terms whose precision scales with counts (`density_lambda_factor`) and,
per object, as the landscape prior (§7.1).

### 6b.2 Record — the RNA-anchored evidence factor

Deleted with the relay policy on 2026-09-09. The ruling that survives: the certified flux is a MESSAGE,
one hop, measured at boundaries and never solved — now §6b.13's flux level, priced by the junction–exon
pair.

### 6b.3 Record — the certified flux is a message (owner, 2026-08-25)

The spliced-fragment observation is a one-hop imputation from the flank boundary into the exon and may
not exist as an appendage beside the message framework: the sender publishes it unchanged and the
recipient decides. The ruling now lives in §6b.13.

### 6b.4 The exon → boundary message is the splice-in map read backwards (owner, 2026-09-02)

An `intron|exon` boundary's unspliced crossing holds gDNA and the unspliced RNA of the transcripts that
span it; the exon holds those AND the mature RNA that arrived by the splice junction. So the exon →
boundary message removes the mature share — the SPLICE-OUT direction. Rescaling the boundary's spliced
density into the exon's frame by the enrichment ratio, subtracting, and rescaling back, the enrichment
CANCELS and only the face's own spliced-to-unspliced ratio survives: `f_b = f_E · (U_b + S_b) / U_b` —
the splice-in face map solved for the boundary, so the boundary evaluates the exon's likelihood row AT
the map (`transfer_rows.splice_out_row`). Three rulings the measurements forced: (1) the exon publishes
its OWN evidence only — its strand row — and only when the solver's derived strand deadband
(`region_init.strand_evidence`) declares the channel live, so an unstranded library's exon says nothing,
exactly, with no constant; (2) both components convert counts to densities with ONE opportunity
treatment, the capture-blind geometric opportunity for gDNA and RNA alike (a capture-aware opportunity
on one component alone re-introduces a level across locales — a 12 % harm on sparse probes); (3) the
width is the marginal over the measured ratio (`log ρ ~ N(log S/U, 1/S + 1/U)`), not a uniform blur in
log-odds, because noise in the ratio moves the boundary's log-odds by `σ/(1 − f_b)`, unbounded at the
pure-gDNA vertex. The premise — spliced and unspliced fragments at the same face share capture affinity
— reads a bias of `a ≈ 1.3` under benign capture, recorded and not corrected. Measured: ladder stranded
capture-ON 0.987–0.995×, unstranded byte-identical (`tests/calibration/test_transfer_policy.py`).

### 6b.5 The boundary → intron message is the boundary's own strand row, verbatim (2026-09-02)

An intron and both of its `intron|exon` boundaries hold ONE unspliced population: mature RNA reaches an
exon by the splice junction and crosses neither boundary (certified: 0.02–0.03 % of the crossing mass on
the ladder, every case a ≤ 100 bp intron inside a fragment's unsequenced mate gap — an accepted
accumulator residual). Under §6b.4's one-opportunity rule the map is the identity, so the row is
delivered verbatim: the boundary's OWN strand row (`simplex_logodds.strand_row_logodds`, the variance
frozen at the boundary's incoming belief — a source-side read), never its belief. The hop adds nothing:
stage 0 on certified truth reads zero excess variance over counting between an intron's composition and
its boundaries' on every panel, so no widening ships. The licence: the boundary must admit the intron's
single strand set (`transfer_rows.boundary_shares_strand`); a terminus flag does not refuse; the deadband
is `tau_lam > 0` at the boundary, so an unstranded library sends nothing. Measured node-locally at the
receiving introns: −26…−38 % on the ladder's stranded capture-ON rows; the reversed row multiplies the
destination error by 1.7–30×.

### 6b.6 The exon|exon terminus boundary and its outside exon — the licence counts the spliced crossing (2026-09-02)

A transcript terminus at an `exon|exon` boundary covers exactly one flank, the INSIDE; the OUTSIDE is
read off the flag alone (TSS+ and TES− bodies extend genomic-right, so the outside is the left flank;
TES+ and TSS− the reverse; termini pointing both ways give no side — `transfer_rows.outside_flank`).
Certified: the outside flank's composition matches the crossing (13,985 ladder pairs) while the inside
flank differs by +0.25…+0.51 nats. The licence is NOT verbatim: a mature fragment that crosses an
`exon|exon` boundary and splices within its own extent is counted in the boundary's SPLICED bank, not
its unspliced crossing, while a fragment contained in the outside piece has no junction by geometry (on
1,003 deep pairs the boundary's true gDNA share exceeded its outside exon's by +0.11, entirely mature
RNA). So the licence is §6b.4's map with the spliced crossing as S: `f_b = f_O (U_b + S_b) / U_b` — the
residual +0.004 after the map against +0.055 before (`g50 ss.99 ON`). Three messages, every one a
shipped constructor: the outside exon's own strand row to the boundary through `splice_out_row`; the
boundary's own strand row to the outside exon through the face map (`transport_row`); and the composed
transport — what the splice-in map and the edge delivered to the outside exon, carried one hop further
into the boundary, the half that reaches unstranded data. The no-echo law is structural: a boundary can
never hear its own row back through the exon. The inside flank is never a destination here (§6b.12's
level rule serves it).

### 6b.7 The abundance-discrepancy rule — superseded

Superseded 2026-09-04 by the level rule (§6b.12): the two hypotheses for a total-abundance discrepancy
across a face — enrichment or new RNA — cannot be told apart, so no hypothesis is chosen, the value is
the boundary's level and the precision is dampened by the pair's own discrepancies.

### 6b.8 The alternative splice site, and the hop premise — the per-pair discrepancy rule (2026-09-02/03)

An `exon|exon` boundary carrying a splice junction and no terminus: one isoform continues contiguously
across it, the other splices out there. Read off the flag alone (a DONOR bit marks the intron's LOW end
on either strand; an ACCEPTOR bit its high end), the two flanks are C, the flank on the junction's
intron side, and E, the flank where both isoforms are exonic (`transfer_rows.junction_flanks`). C shares
the boundary's full unspliced crossing, so it is §6b.6's law with the spliced crossing alone; E holds the
crossing plus the isoform that splices out at this face, measured as the face's route flux `F`, so it is
§6b.4's law with `S_b + F`. Certified on the ladder: C − pred +0.003 and E − pred −0.002 in f off
capture; the flanks swapped open ±0.10 OFF, ±0.25 ON. Each flank's own strand row travels to the boundary
through `splice_out_row` and the boundary's to each flank through `transport_row`, deadband-gated.

**The discrepancy rule, per pair, and nothing pooled** (owner, 2026-09-03). Delivered at counting width
the messages harmed the alt-ss boundaries on benign capture-ON rows (`g50 ss.99 ON`: 126 → 249
node-locally): at every probed pair the gDNA's crossing-to-contained ratio runs ~10 % above the mature
RNA's, by an amount that depends on the locus. So each pair holds two witnesses of one quantity — the
boundary's own strand mode and the flank's mapped to it through the licence — and where they disagree
beyond counting, that pair's messages are widened by the excess `max(0, d² − v)` (`blur_row`); no mode is
shifted, no library-level quantity exists. ⛔ A pooled step per hop kind was landed for a day and
refused: other pairs' behaviour does not predict this one's, and on the ladder the two forms are within
0.25 % of each other on every row (`ISSUES.md` CLOSED/REFUSED, `the-pooled-hop-step`). Owner decision
recorded, not re-litigated: the capture-ON rows at low gDNA on the benign and junction panels
(+1.5…+9 %), where a systematic offset of the licence hides below each pair's counting — a shift would
recover it and shifts are refused. The per-pair rule is the per-hop dampening the completion contract
asks for: a chain of k hops accumulates k pairs' widths, none a constant and none pooled.

### 6b.9 The rebuild's foundation — the intron's forward, the splice-in map and the edge (owner, 2026-09-01/02)

**The paradigm.** One node type, one message, one boundary case at a time; every message a
COMPOSITION carried across ONE face by a derived map, so no level ever crosses a capture cliff and no
constant exists anywhere; silence, never a zero-filled channel, where there is nothing to say. ⛔ The
founding refusal: a gDNA level carried between locales under capture is refuted by probe placement alone
(measured on the adversarial probe panels), so a cross-locale level was refused as the rebuild's basis.

**The intron's forward.** The intron and its boundary share their unspliced population, so the intron's
own factory row is delivered unchanged at `intron|exon` boundaries; stage 0 on certified truth priced the
hop's cost at zero beyond the row's own width, and the blur constant that survived the prototypes was
deleted on the ladder A/B.

**The splice-in map.** The intron row travels into the exon through the splice-in FACE MAP at every
licensed face (`face_is_licensed`: no terminus, the same strand set), `face_map_lambda` monotone with the
certified flux capping the claimable gDNA share, `transport_row` reading the preimage and widening by
the face's own measured-ingredient variance `trigamma(n_u+½) + trigamma(n_s+½)`. Two licensed faces sum
as independent witnesses. Rows ride every sweep — the final-sweep-only citizenship was built and refused,
because true information training the prior bootstrap is where the blind-row value compounds. Ladder:
stranded ≤ silent 7/8; unstranded `g50 ON` −68 %, `g98 ON` −74 %. The adversarial probe panels (the
junction-probed and sparse-probed twins of the test chromosome) were built for this rung.

**The edge and the terminus excursion.** The `intergenic|exon` edge's profile-likelihood lower bound was
superseded by rule 5's level (§6b.12); the enrichment-ceiling upper side is owner-refused as
over-engineering and the zero-gDNA edge residual is an accepted error. The three terminus mechanisms
prototyped together on 2026-09-02 re-entered one at a time as §6b.6, the two-phase backbone and the level
rule (§6b.12).

### 6b.10 The scan seam — superseded

Superseded 2026-09-04 by §6b.12: under the two-phase backbone the passes compose what a node holds, and
the seam's ledger, consumed set and hop budget are unnecessary by construction.

### 6b.11 The propagation is formal forward-backward, and the recipient decides (owner, 2026-09-04)

**The sender does not decide; it sends.** For almost every message the source has no choice to make —
it states its claim and publishes it. **The recipient decides, during the propagation phase**, what to do
with what it received: FORWARD it, MODIFY it (the face's map, the hop's width), or STOP it (the
composition cannot cross here). Every rule the policy has is a decision made at the node that received,
never a filter applied at the node that spoke. **And the propagation is formal forward-backward: when it
ends, every node has received two messages, one from each neighbour, except the two nodes at a chain's
ends, which receive one.** A hop that carries nothing must still ARRIVE, as an explicitly uninformative
message, so the solve can tell silence from ignorance and every node is solved from two honest messages
(the completion contract, §0c.0e).

### 6b.12 The two-phase skeleton, the message's lanes and the level rule (owner rulings, 2026-09-04)

**The backbone is the owner's two phases and nothing straddles them.** `propagate` — a forward pass then a
backward pass; at each hop the RECIPIENT receives what its neighbour sends (the sender's own claim
composed with what the sender holds from its far side), decides to STOP, FORWARD or MODIFY it, and holds
the result; beliefs do not change; when both passes end every node holds one message from each
neighbour it has, `SILENCE` being a message and a missing neighbour (`NO_NEIGHBOUR`) not one. `solve` —
every node once, from its own evidence, the two held messages and the gDNA hyperprior. The names are
`prepare / propagate(backward) → receive(source, destination) / solve(from_left, from_right)`; the
per-sweep object is `Prepared`. Measured before the ruling: the formal form with the same messages won
both halves of the ladder against the one-hop form, 7/8 and 7/8 (`policy_prototype.py`).

**The message's lanes.** A node's unknown is its COMPOSITION on the simplex — two degrees of freedom
where both strands are live — and, where composition cannot cross a face, the LEVELS of the three
populations. So the message carries optional lanes: the gDNA-versus-RNA profile, and a level claim each
for gDNA, RNA+ and RNA− (§6c, §6b.13). The level lanes are for the faces composition cannot cross AND for
AMBIG regions, whose two degrees of freedom a single-stranded neighbour can impute one at a time. The
tilt needs no lane: both strands' bounds constrain it through the shares.

**The level rule** (owner design, landed 2026-09-04; replaces §6b.7). At a face composition cannot cross
— the region INSIDE a terminus, `exon|exon` and `exon|intron` alike — the gDNA LEVEL crosses: the
boundary's own strand profile carried through the level-kept map (the inside's gDNA share is the
boundary's share times its crossing density, times the inside's opportunity, over the inside's OWN total
— an observation), its shape preserved, blurred by both totals' counting and by the pair's own
discrepancies: the excess of the totals' disagreement over counting and, where both strand channels are
live, the excess of the two strand modes' disagreement over counting. The VALUE is kept; nothing is
pooled; no hypothesis is chosen for the discrepancy. **A level is made from the sender's MEASUREMENT
only** — its own claim and its total — never from what it holds: an imputation is not re-issued as a
level (a Gaussian summary of a one-sided profile and a level made from the held profile were both
refuted). With no own claim the boundary still sends what every library measures, the crossing total's
upper bound on the inside's gDNA density.

**Rule 5 is a level, and it is one-sided** (the ladder's verdict, 2026-09-04). The `intergenic|exon`
edge's crossing is structurally pure gDNA, so its COUNT measures the gDNA level the exon continues; the
exon converts it through its own total — at least the edge's gDNA density, at the count's exact Poisson
width — and nothing above. The upper side has no honest form: capture enriches a probed interior over
its edge by an amount no local witness measures (2.3× on the ladder), and every upper side tried pulled
unstranded probed exons toward a centre below the truth (the sparse-probe panel's `g98 ss.99 ON` 36,645
against 6,981). A ZERO count is vacuous: a dark edge under capture is not an empty one. The zero-gDNA rows
are the landscape prior's to win (§7.1).

**Every message-policy idea is judged twice, and the two readings are never pooled.** The gDNA
landscape prior is fitted after the first pass, on that pass's solved gDNA, and comes back in the second
— so a first-pass error poisons the prior and the prior returns it everywhere
(`ISSUES: gdna-landscape-trains-on-false-positives`). A mechanism is read at PASS ZERO
(`calib_refit_iters = 0`) and with the full pipeline, side by side; a win or a refutation read only
through the full pipeline is not attributed. The first pass does not have to solve every node; it has to
solve ENOUGH nodes confidently to train a landscape that solves the rest.

### 6b.13 The RNA level lanes — the both-stranded locus (owner rulings 2026-09-08; landed 2026-09-08/09)

The last propagation case: REGIONS and BOUNDARIES that admit RNA on both strands. Their composition has
two degrees of freedom and their own strand split pins only a LINE in `(f_g, f_+)`, so the gDNA share
there comes from the messages and the prior alone. **The key fact: a message about one strand's RNA level
IS a message about the gDNA share** — `½ f_g + κ f_+ + (1 − κ) f_− = p_+` links them, so a LOWER bound on
RNA+ is an UPPER bound on the gDNA share: the side the gDNA lane cannot give. The bracket theorem (gated
on a hand-built node, both κ): three lower bounds — gDNA, RNA+, RNA− — plus the node's own strand counts
give a two-sided gDNA share; removing any one opens a side. The owner's four rulings: the three levels
travel TOGETHER in one message (components only where measured, empties forwarded); the certified flux
joins as an RNA source; ONE representation everywhere — profiles on the solve grid, no Gaussian summary
anywhere in the transfer policy; the bar is about one percent of a row.

* **The lanes.** `Message.level_rna_pos` / `level_rna_neg`: a strand's level as a profile over
  `u_s = log(ρ_s / ρ_ref,s)` (`ρ_ref,s` the library's strand-`s` unspliced density over its
  single-strand exons — a coordinate). Faces from the flag bits, per strand: strand `s`'s level crosses a
  face iff the boundary carries none of `s`'s four bits and both nodes admit `s`; across `s`'s own
  junction it enters `s`'s intron and not `s`'s exon; a terminus of `s` stops `s` both ways. The intron
  test is PER STRAND (`BlockContext.exon_pos` / `exon_neg`): a region that admits `s` and carries no exon
  of `s` is `s`'s intron whatever the other strand does there. Two-sided only between an intron of `s` and
  its own boundary; lower-only everywhere else. **Every hop pays the pair's price, and the witness of the
  strand's abundance is the column split's asymmetry** (`lanes.LevelLane.witness`, 2026-09-09): both column
  counts' counting plus the disagreement, beyond its own counting, between the two nodes' estimates of
  this strand's RNA — `count − other` on each (gDNA splits evenly and cancels; the protocol's contrast
  `|1 − 2κ|` drops out of the ratio) — carried with the level across empty nodes. A node with no
  asymmetry is DARK, and two dark nodes agree whatever their column densities do. The decisive
  measurement: counting alone carried a lit intron's sharp upper side across a probe cliff unpriced (the
  ladder's `g05 ss.99 ON` read a 93 % RNA junction as 86 % gDNA); the split witness keeps that row.
* **The sources.** A single-strand node's own claim read as its live strand's RNA level
  (`rna_level_of_profile`); and the certified flux at each of an exon's junctions as that strand's level
  at the exon (`flux_level`: the spliced count's Poisson likelihood on the route rate's own opportunity,
  lower side; one hop, boundary → exon). **The junction's rate is an ESTIMATE of the exon's abundance,
  priced by the node pair** (owner, 2026-09-08): the spliced count at its route rate against the exon's
  count of that strand per RNA opportunity, read on the column the strand's RNA reads on (`read_column`).
  ⛔ It stays lower-sided: the two-sided estimate over-claimed at the probe cliff (the sparse-probe zero
  control 94 → 448). **An empty exon piece beside a lit junction is a source too** (landed 2026-09-09):
  the level is built there, priced by `hop_price` on the piece's zero count, and emitted with the flux's
  own witness (`ISSUES: the-empty-flux-source-at-the-junctions-counting-alone`).
* **The delivery at AMBIG nodes** (`PsiMessage.cube_rows`, `cube_row`): the held levels per strand as ONE
  row over ψ's `(λ, θ)` cube — at each cell `f_s = (1 − σ)(1 ± τ)/2`, the profile read at
  `log(ρ_s / ρ_ref,s)`; a one-sided profile stays one-sided (gated). The backbone adds the row inside the
  AMBIG solve, final solve only; absent, byte-identical (gated).
* **The upper side at single-strand nodes** (`rna_row_of_level`, `_ceilings`; landed 2026-09-08). An RNA
  level of the node's live strand says "at most this much gDNA": the held level read as a λ row through
  `f_s = 1 − σ`, bounds intersected. **Read only from a face that sent no composition** — a licensed
  face's splice-in map already carries the flux as its cap; the naive form that reads every level and
  every flux counts them twice (the weak-κ zero control 42 → 1,540).
* **Measured at landing** (full pipeline, halves apart): the ladder wins every non-zero row of both halves
  (unstranded 6/6, worst 1.000×; stranded 6/6, worst 0.994×; `g98 ss.99 ON` 201,578 → 176,703).
* ⛔ Refused with numbers, recorded in `ISSUES.md`: the gDNA lane emitting on every face (a one-sided
  floor arriving at a node with no channel of its own is a tilt, not a floor: ladder `g05 ss.99 OFF`
  +7 %); and a both-stranded node emitting its gDNA level (`ISSUES: ambig-node-as-a-gdna-source`).

### 6b.14 The sj+terminus boundary — the terminus decides the side, the junction places its flux (2026-09-08)

One boundary carrying a splice junction AND a transcript terminus of the same strand: a transcript that
starts at an internal exon's edge (RUNX1) or ends at one (LARGE1); 386 on the ladder. Until this ruling
both rules that could serve the inside exon refused it — the terminus rule because a junction was
present, the junction rules because a terminus was. **The rule.** `outside_flank` reads the orientation
from the terminus bit alone (a junction does not change which flank the terminating transcripts cover);
the level rule serves the inside exon and the outside map the outside exon exactly as at a plain
terminus. The junction's measured flux is placed where the junction's exon is (`junction_exon_side`):
added to the outside flank's population in the map when that flank is the outside, and to the boundary's
total in the level rule's totals' disagreement when the junction's exon is the inside. ⛔ The flux is not
a crossing: the level rule's strand-mode prediction keeps the crossing's own scaling. **What it is
worth.** The case-specific part — the boundary plus the inside exon — is 0.6 % of an in-scope row at the
counting floor; the rule is neutral there and 14 % better locally on the unstranded zero row (the inside
exons 182 → 63); whole library every ladder row within 0.3 %. Under sparse probing the deferred
`g98 ss.50 ON` row loses 4.8 % — the level rule's total bound at a depleted face — reported, not a
target. The substrate is the test chromosome's sj+terminus block (`sjterm` · `capsjterm`).

### 6b.15 The locus is the unit of the solve — a terminal receives nothing, and ψ's read-out is chunk-exact (owner, 2026-09-11)

**The decomposition is the locus, exactly as the EM already does it.** An intergenic region — any REGION
admitting no RNA strand, the predicate the SOLVE gate already locks — is a TERMINAL: structurally pure
gDNA, solved and fixed before a message exists, so nothing needs to cross it and nothing may. Verified on
the human chain before it was made structural: of 1,206,202 composition faces and 4,621,302 lane faces
the shipped policy built, none delivered into an intergenic region, none of a terminal's 65,852 gDNA-lane
faces ever carried anything (a terminal has no own level), and the policy wrote no held message at one
in either pass. The backbone now REFUSES to ask a kernel for a hop into a terminal (`sweep._pass`) and the
lane no longer lists faces from one, so no future policy can move the boundary condition; the terminal
predicate is `is_region & g1_locked`, one definition (33,120 slots on the human chain, exactly the
intergenic set). The chain therefore breaks into 33,018 loci (median 19 slots, the largest 2,477 =
0.12 % of the chain) and `sweep.solve_chain` solves it a LOCUS BLOCK at a time
(`region_chain.locus_blocks`): each block on its own slice of every input, reading one slot beyond
itself where that slot is the terminal its last node receives from (the halo is load-bearing: without it
the last node hears an open side instead of silence and its own flux is not read).

**The only information that crosses a block boundary is the policy's LIBRARY.** `Policy.library(view)`
runs once over the whole chain on a `ChainView` — observations and geometry, no beliefs, so a cross-block
reduction over beliefs has no field to read — and `prepare(ctx, library)` sees one block. The transfer
policy's library is three reference densities (the gDNA lane's, each RNA lane's) and whether the strand
split is a live witness, which is now derived from the deadband (`region_init.strand_discriminability`)
rather than from a per-slot solve. The intron factory's rows travel on the context (`factory_rows`, the
very array ψ adds as its λ-factor), which retired the policy's grid-keyed row callback.

**ψ's read-out is chunk-exact, and that is what makes the block size a knob rather than a choice.** The
shipped single-strand read-out was not: a fancy index on the last axis in `_regrid_global` returned an
F-ordered ψ whose row reductions summed in a row-count-dependent order, and the BLAS matrix–vector
moments dispatched a different kernel at one row — splitting any real 255-row tile moved ~70 % of its rows
by ≤ 1e-15, and halving `_SOLVE_BLOCK_BYTES` already moved slots on the shipped path. The repair (a
contiguous ψ, per-row moment sums) moves the answer by ≤ 3.1e-15 per slot per sweep, does not amplify
through four sweeps and three refits (the final belief ≤ 3.1e-15, `has_composition` never flips), leaves TPM and
effective lengths bit-identical on a real library and every aggregate of `calibration_vs_oracle.py` at the
last ulp with `ruler_n_moved` identical on all 16 conditions; the owner accepted it as identical to a
tolerance (2026-09-11). With it, every block size gives the same bits (gated on the toy for six sizes and
on a real 2.09M-slot sweep for eight), so `CalibrationConfig.sweep_block_slots` sets only the working set.

**The message layer is refit-invariant, so the refit sweeps share it** (derived and measured
2026-09-11, the first step after the decomposition). Everything the layer reads is on the context —
observations, geometry, the factory rows, the incoming belief's ``belief_fg`` and the liveness bits
``has_own_composition`` (`tau_lam > 0`, the one bit of the self-solve a policy may know; the context no longer
carries the self-solve object) — plus the library and the grid, and never the prior; and `calibrate`
resets the belief before every sweep. So for one grid every refit sweep's messages are the same:
measured on the human chain, sweeps 1–3 deliver identical ψ rows and cube rows to the bit and every node
hears the same thing. `message_cache.MessageCache` holds one grid's delivered messages, content-keyed on a
digest of every input the layer reads (a changed belief, row, count, library, grid or policy misses —
each channel gated by perturbation), sparsely (0.17 GB of rows plus 0.39 GB of cube rows per grid on
the 876k library, against a dense 2 GB); a refit sweep pays its two ψ solves and is served the rest.
Diagnostics never read from it. Pass 0's grid is never reused, so it is not held. On the 18.6M-fragment
library the refit grid is stable (`n_grid` 138 for all three refits), refits 2 and 3 are served entirely
(38 s each against 176 s), the run reads 0.65 of its wall in two back-to-back pairs, and the cache holds
2.68 GB (peak 15.0 → 17.8 GB) — the cube rows as the float32 the AMBIG solve casts them to; 4.1 GB as
float64.

**The rules are typed tables, and a face is a side** (2026-09-11, the port's data layout). Every
directed face is one of a node's two sides — it hears from its left neighbour or its right — so the
recipient's composition rule is a KIND and its parameters at ``(destination, side)``:
`messages.faces.Faces` holds ``(n, 2)`` tables (the kind, the face's unspliced and spliced counts,
the boundary's and far region's gDNA opportunity, a blur width, the level rule's width) and indices into
a row store of the ``(K,)`` maps; five kinds cover every shipped message — FORWARD, TRANSPORT (boundary →
region through the face map), SPLICE-OUT (region → boundary, the map read backwards), EDGE (the
intergenic|exon edge's one-sided level) and LEVEL (the terminus's level rule) — and `Faces.apply` is
the one place their arithmetic lives. A face carries ONE rule: the builders' faces are disjoint by
construction (the splice faces serve intron|exon pairs, the edge rule gene edges, the terminus rules
unlicensed faces, the alternative splice site junctions with no terminus), so the table refuses a
second rule at a face as it refuses a rule at a face that does not exist — the earlier "a later builder
replaces an earlier one" precedence had no instance on the toy or the human chain and was a hidden
assumption, not a rule. The level lanes hold their faces as ``(n, 2)`` bits and their junction flux as
a row table. The 1.2M closures, the 4.6M-tuple face sets and the neighbour-pair enumeration are gone,
bit-identically; `_SolveSite` needs no neighbour arrays. A compiled pass reads these buffers directly.

**Threads are the wrong tool for this sweep, and the executor waits for the port.** Measured on the real
sweep: the locus-split passes at 8 threads 0.83–0.94× (GIL-bound Python), ψ's grid solves 2.06×, the same
passes in 8 forked processes 6.16×. The owner's decision: no parallelism until the C/C++ port of the block
solve, which parallelises there; this ruling delivers the structure the port lands on and the memory half
of the problem — measured end to end on the 18.6M-fragment library at 8 threads, two back-to-back pairs:
the run's peak RSS 33.2 → 14.9 GB and 32.6 → 15.0 GB, the wall 0.97–0.98, the sweeps 0.96–0.97, and the
peak now set outside the solve (`build_region_geometry`, `init_beliefs`).

## 6c. ψ's composition is a point on the simplex, and closure is structural (2026-08-17)

**The composition has two degrees of freedom, not three.** ψ solves a point on the 2-simplex,
parametrised by `λ` (the gDNA-vs-RNA LEVEL) and `θ` (the RNA-internal TILT, a share with no absolute
scale). The composition is their image — `simplex_logodds._compose`:

    f_g  = the ½-quantile of the λ posterior        RNA total := 1 − f_g
    f_pos = (1 − f_g)·w        f_neg = (1 − f_g)·(1 − w)

so `f_g + f_pos + f_neg = 1` identically; 100.00 % of published objects close on every annotation class
(against 74.7 % of REGIONs before, when the RNA fractions were independent posterior means and the
closure error was the posterior's skew). ⛔ Taking means everywhere also closes and is refused: it scores
1.352 / 1.573 / 3.756 on the three in-scope strata and 1.801 on the zero control
(`vertex_ceiling.py --arm psi_mean`) — the median is closer to truth at both simplex vertices, where
49–83 % of in-scope error lives. Nothing is rescaled at publication.

**The ½-quantile is continuous and is read on λ.** Snapping to a lattice point put up to half a grid step
into `f_g` (`TRAPS: deriving-one-coordinate-propagates-its-error`); the interpolation must be on `λ`,
where the lattice is uniform — in `f_g`-space a concentrated posterior comes back biased toward ½ by
2.71e-03 at `n_grid` 60, on λ it returns its own grid point to 2.2e-16
(`TRAPS: interpolate-on-the-axis-where-the-lattice-is-uniform`). Admissibility is enforced inside the
map: a slot with no counts publishes `(0, 0, 0)` — "no data", not a composition claim — and a slot with
no admissible strand has composition `f_g` alone, which is honest and not a closure failure.

**The level lane — gDNA is always conveyed, as a lower bound** (landed 2026-09-05). Half the pass-zero
error on the ladder sat at nodes no message reached (51 % on `g50 ss.50 OFF`, 59 % on `g50 ss.99 ON`),
and the mechanism is the EMPTY node: 52 % of the ladder's exon pieces have no total (a piece shorter than
a fragment, or a dark piece), every rule's licence asks its flank for a total, so the boundaries on both
sides of such a piece heard nothing. The ruling that closes it:

* **A level is absolute.** `Message.level_gdna` is a PROFILE over `u = log(rho / rho_ref)` on the solve
  grid (`rho_ref` the library's structurally pure gDNA density, a coordinate choice); it needs no map and
  no knowledge of its recipient, which is what lets it cross a node with no total, and it is a profile
  because what travels on it is one-sided. A node's own composition profile becomes a level through its
  own total (`level_of_profile`) and a held level becomes the recipient's composition profile through ITS
  total (`profile_of_level`) — one map read both ways.
* **The lane is the default rule of every directed face without a composition rule** — strand-change
  faces, termini pointing both ways, the AMBIG complex, and every face into or out of an empty node. STOP
  by omission is impossible by construction (gated); a face with a composition rule sends composition
  only, so no witness is counted twice.
* **An empty node is transparent** — it holds levels only and forwards them unchanged. **A full node emits
  the intersection** of its own level's lower side and the priced level it holds — the pointwise tighter
  of two bounds, never their product. ⛔ The product form RATCHETS: on the ladder's `g05 ss.99 OFF` a chain
  of nine terminus boundaries with one true gDNA fragment each moved 2 → 29 apiece. **The recipient
  prices the hop** — both totals' counting plus the abundance discrepancy beyond it, per hop, nothing
  pooled — **and takes the level as a lower bound.**
* ⛔ **A level that crosses a face says "at least this much gDNA" and nothing more.** Every upper side was
  measured and refused: two-sided by class, +33 % on `g50 ss.99 ON`; the upper side kept for an own
  measurement's first hop only, +6 % on `g98 ss.70 ON`. Each is darkness under capture read as absence.
* ⛔ **What a node holds as a composition is never re-issued as a level** (reading one as a level harmed
  three in-scope ladder rows by 1–3 %), and **rule 5 and the level rule stay as landed** (re-expressing
  them through the lane broke the zero-gDNA controls, `g00 ss.50 ON` +19 %).
* **The prize the law forgoes is priced**: a two-sided lane read −26 % at pass zero on `g50 ss.50 OFF`,
  because off capture the minimum total density of an exon complex truly bounds its gDNA from above. The
  message layer cannot see whether the library's gDNA is enriched; `ISSUES: two-sided-exon-row` owns it.

## 7. Where the error is, by class

Read with §0b: an object class carrying the error is a mechanism, a stratum carrying it is a scope.
Regions and boundaries measure different components — the gDNA/RNA opportunity ratio is 0.25 at a
crossing point against 115.7 at a 100 bp region and 1.19 at 1,000 bp, so a short region is a good gDNA
measurement and says nothing about RNA; carry per-component precision, not one scalar. On an unstranded
library the density model carries the entire own-evidence budget: at κ = ½ the strand λ-term is exactly 0
(`EQUATIONS.md` §5), and the intron factory is what makes such a library solvable at all. Pass-0 scores
honest ignorance as error, which is the wrong question: an object with no own evidence reporting
`f_g ≈ ½` at zero precision is stating a true fact, and the measurement that matters is solvable → right /
wrong → confidently wrong (`solvability_audit.py`). Where it stands (re-derived 2026-09-10,
`policy_benchmark.py --by-class`; `ROADMAP.md` carries the ranking): in scope the residual sits on the
intron's own solve on unstranded capture-OFF and on `exon|exon` boundaries and walled exons on stranded
capture-ON; the deferred stratum is blind because the gDNA fraction cancels from the strand mean, so an
unstranded AMBIG slot has no channel.

### 7.1 The landscape prior — who trains it, where its kernels go, and what axis it lives on (owner rulings 2026-09-06 and 2026-09-10; landed 2026-09-10)

**Three rulings, one gate file** (`tests/calibration/test_landscape_training_population.py`; the arms and
their numbers `ISSUES: the-landscape-training-population-arms`; the instrument
`landscape_training_census.py`):

1. **A node whose only evidence is a bound, or which has none, does not train the prior.**
   `RegionBelief.has_composition` — an own composition channel (`has_own_composition_evidence`), structural
   certainty (`g1_locked`), or a COMPOSITION row received from a neighbour — is published by
   `sweep.solve_chain` from the held messages and selected on by `calibrate._fit_gdna_hyperprior`; a level
   lane, a ceiling and a cube row are bounds; the zero-count anchor trains regardless. ⛔ "Any non-flat
   λ-row" is NOT the predicate: `PsiMessage.lam_rows` fuses compositions and bounds (1,476 own-flux
   ceilings at the unstranded zero control; 137k against 111k).
2. **The grid spans every region and boundary the prior is read at** (`fit_landscape(domain=…)`), so the
   training cut changes which kernels are summed and never the axis; the `_MIN_TRAIN` guard measures the
   annotation-admitted population. Found on the gDNA-free golden toys, where the cut left the anchors
   alone: the grid collapsed to the floor (14 → 95 invented fragments of 1,000) and the guard refused the
   refit (52 → 201). Free on the ladder (byte-identical on 14 rows; the full domain widens the step ≤ 10 %).
3. **The E-step on the kernels that have no location** (`landscape._estep_kernels`). A region trained at
   less than one fragment is centred at its wall (the estimator's `max(count, 1)`), its Poisson kernel is
   flat below the wall, and normalised to unit mass it spreads that mass uniformly under the wall — a
   likelihood used as a density. At the unstranded zero control that put 1.0 % of the landscape's mass
   above −1 decade on the 1,310 anchors shorter than 100 bp and 0.6 % on the zero-count exons, and a
   blind exon's median was decided by that tail against ψ's Beta(½,½) reference's ½-nat/λ slope. The refit
   loop's previous landscape now places those kernels (kernel × P_prev, renormalised); COUNTED kernels
   keep their own location, so an enriched minority cannot be competed away. Nothing chosen.

The anchor holds 46–94 % of the estimator's weight on every ladder row and the refit loop de-entrenches
(false gDNA trained 777k → 71k across three refits at `g00 ss.50 OFF`); the zero-row residual was the
landscape's tail, not entrenchment.

**Measured** (the ladder, whole-library |gDNA − truth| in fragments; `policy_benchmark.py`): the zero
controls 149,552 → 674 and 166,837 → 283 (unstranded OFF / ON), 14,324 → 529 and 15,365 → 262 (stranded),
the rule and the E-step together; every in-scope contaminated row 0.95–1.00× of the merge-day default;
the deferred rows 0.96–1.02×; the zero-RNA controls within ±4 %. The test chromosome: wins or ties all 30
rows. Under both, `silent`'s zero controls fell to a few hundred fragments as well (the prior alone solves
a gDNA-free library), so `policy_benchmark.py`'s "beats silence" count reads 6/8 per half with every
losing row a zero control at that floor; every contaminated row still favours `transfer`.

⛔ Refused with numbers, in `ISSUES: the-landscape-training-population-arms`: every reading of "a bound"
that reaches the DELIVERED rows (under capture the probed exons' one-sided rows ARE the enriched mode's
witness), the likelihood-kernel estimator, the all-kernel E-step and six refits.
