# Intergenic alignments as gDNA candidates (review packet, 2026-10-15)

> **Subsequent review:** `INTERGENIC_GDNA_REVIEW.md` confirms incomplete paired hypotheses in the baseline scanner
> and reproduces the length-law collapse. Arm 2 is not approved to build as specified. The original packet and its
> measurements are retained below; the complete-pair contract, normalized opportunity and revised work order are in
> `SPLICE_ARTIFACTS_IMPLEMENTATION.md` and `RELEASE_0_8_0_FIX_PLAN.md`.

This packet is for an external reviewer and the project owner. It covers `ISSUES: multimapper-intergenic-alignments`:
what Rigel should do when a fragment's aligner reports an alignment outside every gene. Part A is the reviewer's
prompt; Part B is the packet.

**Where things stand.** The simplest form of the repair (arm 1, §5) is prototyped in a worktree and measured on a
locally aligned panel with exact truth labels (§6–§8). It helps the two contaminated capture-OFF conditions, is flat
at low contamination, and hurts the capture-ON condition; the cause of the capture-ON loss is found and confirmed
(§8). A derived form (arm 2, §9) is proposed and not built. Nothing is in the production source. All
numbers come from `main` at 41899ae4 and the scripts in §11. This file is sandbox: nothing in it is authoritative
until it moves to `DESIGN.md`.

---

## Part A — the reviewer's prompt (copy everything between the rules)

---

You are an independent technical reviewer. Your background is mixture models and EM for RNA-seq quantification
(Salmon, RSEM, kallisto or similar), multimapping reads, and short-read alignment.

You are reviewing a **design** for **Rigel**, a Bayesian RNA-seq transcript quantifier that separates RNA from
genomic-DNA (gDNA) contamination. The design concerns fragments whose aligner reports several placements, at least
one of them **outside every annotated gene**. Today Rigel drops such a placement before its EM, so a gDNA fragment
from an unannotated copy of a gene (a processed pseudogene, a repeat) can reach the EM with only its spliced alignment
across the parent gene, and is then counted as certified RNA.

**The owner's principle** (2026-10-15): "whenever a fragment aligns to an intergenic region we MUST honor the gDNA
hypothesis. The challenge is: what is the gDNA candidate? We could model each intergenic region as its own locus. Or
each chromosome. Or the entire genomic gDNA as one component. This is an unknown research question."

**The decision state.** One form (arm 1) is prototyped and measured: it helps the two contaminated capture-OFF
conditions, is flat at low contamination, and hurts the capture-ON condition, for a reason the packet diagnoses and
confirms. A derived form (arm 2) is proposed and not built. The work sits on the path to the 0.8.0
release, which waits for validation on aligned, artifact-containing data.

Your job:

1. Find the weaknesses in the derivation (packet §4) and in both arms (§5, §9).
2. Answer the owner's research question: what should the gDNA candidate for an intergenic placement be?
3. Judge the diagnosis of the capture-ON regression (§8): is it right, and is it complete?
4. Improve the proposal, and say what evidence should gate building it.

### Materials

- **The packet:** Part B of this document. **Read it first and completely.** It covers Rigel in brief, today's
  behaviour read from the code, the likelihood derived from a generative model, the two arms, the aligned test panel,
  the measured results, the capture-ON diagnosis, the proposal, and the open questions.
- **The evidence bundle** `rigel_intergenic_review.zip`, if the owner attached it, described by its `HANDOFF.md`: the
  code the packet cites, the prototype's diff and tests, and the measurement scripts with their outputs.
- **The panel's truth.** Every read name carries its fragment's true origin (transcript or genomic interval), so
  every per-fragment statement in the packet is exact for the simulated library. The simulated panel is private and
  not in the bundle.

### Rules your proposals must respect

These are the owner's standing rulings. You may argue that a ruling costs accuracy, but only with a measured or
derivable cost. Do not design around one silently.

1. **Axiom 0.** Every genomic object holds at most three populations: gDNA, RNA on the + strand and RNA on the −
   strand. "Mature" and "nascent" are not populations.
2. **No new empirical constant.** Every number is derived, read from the library, or already defined in the code. A
   yes/no cut on a continuous quantity counts as a constant. If a proposal needs one, name it as a discussion item.
3. **Real data are test inputs, never design inputs.** Design from the mechanism.
4. **Judge both errors together:** gDNA counted as RNA, and RNA counted as gDNA, per stratum (stranded or not,
   capture on or off), never pooled.
5. **One mechanism at a time.** Every change states how it would be A/B'd alone against the current tree, on the same
   inputs, scored against truth.
6. **Settled machinery, not re-litigated here:**
   - a per-locus EM (a locus is a connected component of transcripts linked by shared fragments) with one gDNA
     component per locus and a prior from calibration;
   - calibration's choice among a multimapper's placements is a hard, seeded draw;
   - capture enrichment enters the EM through capture-contracted effective lengths (§1).

### What to deliver

A markdown report of about 3–5 pages, with these sections.

1. **Verdict** in at most 10 lines: is arm 2 sound to build? List every **blocker** (a defect that must be fixed
   before building) with the evidence that it is one.
2. **Findings**, numbered. Each has a severity (**blocker**, **fix before building**, **consider later**), the packet
   section it concerns, the reasoning labelled *derived*, *measured* (cite the number) or *conjecture*, and a
   concrete falsifiable test.
3. **The research question.** Rank the candidates for an intergenic placement's gDNA reading: the locus's own gDNA
   component (arm 1); a fixed rate at calibration's density for the placement's region (arm 2); each intergenic
   region as its own component in a coupled EM; one component per chromosome; one genome-wide component. Say what
   each assumes, where it fails, and what it costs.
4. **The derivation** (§4). In particular:
   - the sum over placements against the shipped mean, and the sum's measured cost off capture (§8, §10 item 2),
     where the derivation says it should be right;
   - the capture efficiency at each placement, which cancels for a unique fragment and does not for a multimapper;
   - holding a region's density fixed inside a locus's EM: when is the plug-in's error second order, and when not
     (a region with no unique fragments, a high-copy repeat)?
   - where the outside reading's mass should be booked, and whether the locus's prior must change with it.
5. **The capture-ON diagnosis** (§8): agree or not, and the experiment that would settle it.
6. **Aligner dependence.** The panel was aligned with minimap2, which reports secondaries within 80 % of the primary
   score, per mate and often soft-clipped. STAR's defaults report far fewer. How should a placement's alignment quality
   enter the likelihood, if at all, given that the tool reads NM today but not the clipped bases or the aligner's
   score?
7. **Validation.** The minimum evidence that should gate building arm 2, and what result would stop it.

Be specific and terse. A finding without a test is a conjecture: say so. Where you are uncertain, give the experiment
that would settle it rather than a judgement.

---

## Part B — the packet

## 1. Rigel in brief

Rigel is a Bayesian RNA-seq transcript quantifier that separates RNA from genomic-DNA (gDNA) contamination.

- **The scan.** A single-pass C++ BAM scanner resolves each fragment against the annotation. Each of the fragment's
  aligner placements (a **hit**) is resolved separately. A hit that overlaps no transcript is **intergenic**.
- **Calibration.** It deconvolves every genomic region and region boundary, genome-wide and intergenic regions
  included, into gDNA, RNA+ and RNA−, from strand, density and fragment length. It reports a per-object gDNA mass, a
  gDNA density, and a **capture efficiency** `c` (the object's gDNA yield per position against the bulk; `c = 1`
  everywhere when the library is not captured). It draws one placement for each held multimapper (a seeded draw).
- **The locus EM.** A locus is a connected component of transcripts linked by shared fragments. Its components are
  its transcripts, synthetic spans for unspliced RNA in introns, and **one gDNA component**. Each fragment's weight
  under component `k` is `(θ_k/ℓ_k)·f_k(L)·s_k(o)·e`: abundance over effective length, the fragment-length law at the
  candidate's own length, the strand term (κ or 1 − κ under RNA, ½ under gDNA), and a read-base term (today a mismatch
  penalty per NM). The **effective lengths are capture-contracted**: each component's length counts its positions at
  their capture efficiencies, so `θ/ℓ` is a density per bulk-unit position.
- **The prior.** Calibration hands each locus a gDNA count. The EM's pseudocounts take their odds from calibration's
  gDNA share of the locus's unique fragments and their strength from the number of units with a gDNA candidate
  (`pipeline.em_pseudocounts`).
- **Spliced fragments are certified RNA.** gDNA cannot splice, so a fragment whose every placement is spliced gets no
  gDNA candidate.
- **Intergenic fragments.** A fragment whose every placement is intergenic never reaches a locus EM. Unspliced, it is
  counted as gDNA; spliced, it is counted as RNA that no annotated transcript holds.

## 2. The defect

### 2.1 What the tree does today (read from the code)

- **The scanner** (`bam_scanner.cpp`, the per-hit loop). For each hit of a fragment:
  - a hit with transcript candidates is buffered for the EM;
  - an intergenic hit of a **unique** fragment is buffered (calibration labels it intergenic; the scorer skips it);
  - an intergenic hit of a **multimapper** is **not buffered**. It is offered to calibration as a placement and
    nothing else.
- **The scorer** (`scoring.cpp`, `score_mm_alignment` and `flush_mm_group`). For a multimapper:
  - each transcript keeps its best hit (fewest overhang bases, then highest log-likelihood);
  - gDNA's term is taken from every hit gDNA can explain **that has a transcript candidate** (`n_cand > 0`), and is
    the **mean** over those hits: `log Σ_h exp(ℓ_h) − log n_h`, with `ℓ_h = log f_G(L_h) + log ½ + NM_h·log α`;
  - the unit is certified RNA when no hit gives gDNA a term, or gDNA's term falls below the pruning bound against the
    best RNA candidate.
- **The EM** subtracts `log ℓ_gDNA` of the locus from gDNA's term, so the reading is priced at the locus's own gDNA
  density.

So a gDNA fragment from an intergenic copy of a gene, whose reads also align spliced across the parent gene, reaches
the EM with its spliced hit alone and is **certified RNA on the parent**. If its genic hit is unspliced instead, gDNA
keeps a term, but only at the gene.

### 2.2 How often, measured

On the aligned panel (§6), condition `g50 ss.99 OFF`: 11,672 gDNA fragments have a spliced primary alignment and at
least one secondary. Of those:

| the fragment's unspliced secondaries | fragments | share |
|---|---|---|
| all intergenic | 10,271 | 88 % |
| at least one in a gene | 1,396 | 12 % |
| none | 5 | 0 % |

The measurement classified the secondaries by position, not by overlap with the true origin; §8's join checks the
true origin directly.

## 3. The owner's principle and the open question

The owner's rule: an intergenic alignment is a place gDNA could have produced the fragment, so the gDNA hypothesis
must be honoured. The open question is **what the gDNA candidate is**: the fragment's own locus's gDNA component, a
component per intergenic region, per chromosome, or one for the whole genome.

## 4. The likelihood, derived

**The generative model.** gDNA fragments start at genomic position `x` at rate `λ_G(x) = ρ·c(x)`, with `ρ` the
library's bulk gDNA density per start position and `c(x)` the capture efficiency. RNA fragments of transcript `t`
start at a transcript position at rate `(θ_t/ℓ_t)·c(x)`. The aligner proposes placements `h = 1…H` of the fragment's
bases. At placement `h` the fragment has footprint length `L_h`, orientation `o_h` and read-base emission `e_h`.

**The fragment's likelihood under each origin** sums over the placements that origin can explain:

```
P(u, gDNA) = Σ_{h ∈ G(u)}   ρ·c(x_h) · f_G(L_h) · ½ · e_h
P(u, t)    = Σ_{h ∈ T(u,t)} (θ_t/ℓ_t)·c(x_h) · f_t(L_h) · s_t(o_h) · e_h
```

`G(u)` is the set of placements gDNA can explain (unspliced, implicit or artifact), and `T(u,t)` those compatible
with transcript `t`.

**What the locus EM represents.** Inside locus `L`, the gDNA component stands for `ρ·c(x)` through
`(θ_g/ℓ_g)·c(x)`, where `ℓ_g` counts the locus's positions at their efficiencies. For a **unique** fragment,
`c(x)` multiplies every candidate alike and cancels in the E-step, which is why the scorer omits it.

**Three consequences for a multimapper:**

1. **The sum, not the mean.** gDNA's term is a sum over `G(u)`. The shipped mean divides it by `|G(u)|`. The two agree
   when one placement dominates; with equal placements the mean understates gDNA `|G(u)|`-fold, and a poor extra
   placement (an orphan mate, a clipped secondary) dilutes it. `ISSUES: multimapper-intergenic-alignments` already
   records that uniform gDNA placement derives the sum.
2. **The efficiency does not cancel.** Placements at different positions carry different `c(x_h)`. Dropping it is
   exact only when the unit's placements share one efficiency, which capture breaks: a probed exon and an unprobed
   intergenic copy differ by the probe's enrichment.
3. **An intergenic placement is outside the locus.** Its gDNA rate is the genome's at that position,
   `ρ(x')·c(x')`, which the locus's `θ_g/ℓ_g` does not describe. Calibration already estimates it: the gDNA density of
   the region holding `x'`, shrunk to its landscape prior where the region has little evidence of its own.

## 5. The candidates

| | the intergenic placement's gDNA reading | its density | where its mass is booked | built |
|---|---|---|---|---|
| **arm 1** | joins the locus's gDNA term (the shipped mean) | the locus's `θ_g/ℓ_g`, at the genic placement's efficiency | the locus's gDNA | prototype |
| **arm 2** | a fixed-rate reading beside the locus's components | calibration's gDNA density of its own region, at its own efficiency | intergenic gDNA | proposed (§9) |
| a coupled EM | a component per intergenic region | estimated jointly with the locus | that region | no |
| per chromosome / genome | one shared gDNA component | one estimated density | the shared component | no |

**Arm 2 unifies the last three.** Calibration's per-region density, shrunk to a landscape prior estimated over the
whole library, already is "the genome's gDNA component, evaluated locally": per region where a region has its own
evidence, the library's level where it has none. Holding that density fixed inside one locus's EM is the first step
of the coupled EM; it is exact to first order when the region's own unique fragments dominate its density.

**Arm 1** borrows the locus's density for a place outside the locus. Without capture the two densities agree up to
calibration's noise. Under capture the locus is probed and the intergenic copy is not, so arm 1 prices the
intergenic placement at the probed exon's efficiency. Arm 1 also keeps the shipped mean, so a poor intergenic
placement dilutes a gDNA term the genic placement already had.

## 6. The evidence: a locally aligned panel

- **The panel.** Five conditions of the 16-condition simulated ladder: gDNA 0 % (`g00`), 5 % (`g05`) and 50 % (`g50`);
  strand specificity 0.5 (`ss.50`, unstranded) and 0.99 (`ss.99`, stranded); capture off and on. Ten million
  fragments each.
- **The alignment.** The ladder's own simulated reads re-aligned with minimap2 (`-x splice:sr --secondary=yes
  -j junctions.bed12`, the annotation's junctions supplied), then name-sorted. The oracle ladder's BAMs hold the
  simulator's true alignments and no multimappers; the aligned panel holds 10.5 % multimappers.
- **minimap2's secondaries** are reported per mate, within 80 % of the primary's score by default, and often
  soft-clipped (`66M34S` on a 100-base read). The tool reads NM, not the clipped bases or the aligner's score. Here
  21 % of the records are secondary (`g00`); in the real VCaP library, aligned by STAR, 2.2 % are. This panel
  stresses multimappers harder than a STAR library would.
- **The truth.** Every read name carries the fragment's true origin, so a fragment's fate is scored exactly. The
  transcript and gene tables are scored by `quant_accuracy.py` against the same per-transcript truth as the oracle
  ladder.
- **The baseline is far from the oracle panel's**, and mostly not because of this defect:

| condition | tx Σ\|err\| oracle → aligned | gene Σ\|err\| oracle → aligned | gDNA estimate oracle / aligned (truth) | nascent RNA oracle / aligned (truth) |
|---|---|---|---|---|
| g00 ss.99 OFF | 154 k → 5.08 M (33×) | 10.8 k → 705 k | 1.1 k / 37 k (0) | 2.02 M / 2.23 M (2.02 M) |
| g05 ss.50 OFF | 149 k → 541 k | 13.2 k → 195 k | 0.50 M / 0.22 M (0.50 M) | 1.92 M / 2.11 M (1.92 M) |
| g50 ss.50 OFF | 99 k → 261 k | 13.9 k → 117 k | 4.98 M / 1.84 M (5.00 M) | 1.04 M / 3.10 M (1.01 M) |
| g50 ss.99 OFF | 95 k → 309 k | 11.1 k → 113 k | 5.00 M / 2.35 M (5.00 M) | 1.02 M / 2.60 M (1.01 M) |
| g50 ss.99 ON | 216 k → 565 k | 103 k → 343 k | 4.96 M / 3.50 M (5.00 M) | 0.10 M / 1.43 M (0.15 M) |

Two things make up the aligned panel's gDNA deficit, and neither is this packet's defect:

- **A panel artifact.** The panel's genome (chr21, chr22 and the ERCC spikes, 97.6 Mb) is 18.7 % N, and the simulator
  draws gDNA from the N gaps too. Those reads are all N, so no aligner places them: 0.99 M of the 5 M gDNA fragments
  are unmapped in each capture-OFF `g50` condition (36 k under capture, where gDNA sits on the probes). The oracle BAM
  places them by truth, so the oracle panel counts them.
- **A real defect on aligned input.** Of the gDNA fragments the aligner did place, the tool calls 57 % RNA at
  `g50 ss.50 OFF`, 44 % at `g50 ss.99 OFF` and 36 % at `g50 ss.99 ON`, mostly as nascent RNA. The oracle panel calls
  almost none. Under capture the loss sits mainly in **uniquely mapped unspliced** gDNA fragments (1.44 M of 4.24 M),
  not in multimappers. It is under investigation as its own defect (§10).

This packet's repair moves about 50 k fragments. Every arm below is read as a difference against the same aligned
baseline.

## 7. Results: arm 1 on the aligned panel

**The transcript and gene tables** (`quant_accuracy.py`, fractional assignment, the scan pinned; arm 1 − base):

| condition | tx Σ\|err\| | gene Σ\|err\| | false-positive mass | gDNA estimate |
|---|---|---|---|---|
| g00 ss.99 OFF | +135 (of 5.08 M) | −71 | +59 | +71 |
| g05 ss.50 OFF | −35 (of 541 k) | +62 | −174 | +7,109 |
| g50 ss.50 OFF | **−1,665** (of 261 k) | **−1,431** | −1,110 | +58,689 |
| g50 ss.99 OFF | **−1,865** (of 309 k) | **−1,272** | −1,273 | +46,249 |
| g50 ss.99 ON | **+3,839** (of 565 k) | **+4,191** | +1,659 | −7,478 |

**Both errors by true origin** (every fragment, the shipped sampled assignment at a fixed seed; "called RNA" means
resolved to an RNA pool, "called gDNA" means resolved to gDNA):

| condition | gDNA fragments called RNA: base → arm 1 | RNA fragments called gDNA: base → arm 1 | recovered per RNA fragment lost |
|---|---|---|---|
| g00 ss.99 OFF | 0 → 0 | 37,309 → 37,381 (+72) | — |
| g05 ss.50 OFF | 231,516 → 225,961 (**−5,555**) | 49,877 → 51,433 (+1,556) | 3.6 |
| g50 ss.50 OFF | 2,265,389 → 2,210,340 (**−55,049**) | 97,937 → 101,578 (+3,641) | 15 |
| g50 ss.99 OFF | 1,752,530 → 1,708,250 (**−44,280**) | 88,910 → 90,879 (+1,969) | 22 |
| g50 ss.99 ON | 1,783,091 → 1,790,376 (**+7,285**) | 320,972 → 320,781 (−191) | — |

Without capture, arm 1 recovers 15–22 gDNA fragments per genuine RNA fragment it loses at `g50`, and 3.6 at `g05`.
Under capture it moves the wrong way: more gDNA is called RNA, and slightly less RNA is called gDNA.

## 8. What happens under capture

**The method.** Both arms ran on `g50 ss.99 ON` at the same seed with the annotated BAM kept, and the two BAMs were
joined fragment by fragment. Each fragment was classed by its own records: unique; a multimapper whose placements
are all genic, all intergenic, or **mixed** (both); for a mixed one, whether its genic placement is spliced, and which
placement holds the true origin. The sampled draws do not stay aligned between arms, so a single fragment's flip is
mostly sampling noise (the unique class flips about 420 k fragments each way). Only a class's net change is read; its
noise is about ±1.4 k for the unique class and ±0.4 k for the mixed class.

**The loss sits entirely in the class the repair touches.** The unique class moves by −1,787 gDNA calls (within
noise) and the all-genic multimappers by −281. The mixed class:

| mixed multimappers (a genic and an intergenic placement) | fragments | gDNA calls, base → arm 1 |
|---|---|---|
| gDNA; the genic placement spliced (the defect of §2) | 261 | 0 → 236 (**+236**) |
| gDNA; the genic placement unspliced; true origin the intergenic copy | 15,983 | 13,530 → 14,477 (**+947**) |
| gDNA; the genic placement unspliced; true origin the genic copy | 258,975 | 135,184 → 128,768 (**−6,416**) |
| RNA; the genic placement spliced | 4,003 | 0 → 41 |
| RNA; the genic placement unspliced | 150,048 | 10,092 → 9,558 (−534) |

Of the 275 k gDNA fragments with an unspliced genic placement, 269 k have it on a probe.

**The mechanism is the shipped mean, not the density.** Under capture, a gDNA fragment that also aligns to an
intergenic copy nearly always came from the probed genic copy (259 k against 16 k): capture enriched it there. That
fragment already had a gDNA term at its genic placement. Arm 1 adds the intergenic placement to gDNA's **mean**,
`(e_genic + e_ig)/2`. For the 241,891 such fragments whose genic placements are all unspliced and overlap the true
origin (`nm_gap.py`; shares of all 241,891):

| the intergenic placement against the genic one | a mate pair | one mate only | gDNA's term under the mean |
|---|---|---|---|
| more mismatches (≥ 1) | 33.4 % | 31.0 % | about half: `e_ig ≤ α·e_genic`, with α = 0.1 per mismatch |
| equal mismatches | 26.5 % | 8.5 % | unchanged (the sum would double it) |
| fewer mismatches | 0.0 % | 0.6 % | raised |

A one-mate placement also has a different length term (its footprint is one read), so its rows are approximate. About
64 % of these fragments lose about half of gDNA's term, and move toward RNA. Without capture the same dilution should
happen, while far more of the fragments truly come from the intergenic copy, so the gain outweighs it. That is
inferred: the capture-OFF join is not yet run.

**The density borrowing predicted in §5 is real, but masked.** Pricing the unprobed intergenic placement at the
probed locus's efficiency over-states gDNA, which would move fragments toward gDNA. A mean can never exceed its best
member, so the mean hides the over-statement. A sum would expose it: an intergenic placement as good as the genic one
would double gDNA's term, where the truth multiplies it by about `1 + c(x_ig)/c(x_genic)`, near 1 under capture.

**The confirming experiment.** A diagnostic arm, arm 1 with gDNA's term summed over the placements instead of
averaged (one line in `flush_mm_group`), on the same condition, seed and join (base → arm → arm with the sum):

| class | origin | base | arm 1 (mean) | arm 1 with the sum |
|---|---|---|---|---|
| mixed, genic unspliced, true origin genic | gDNA called gDNA | 135,184 | −6,416 | **+41,361** |
| mixed, genic unspliced, true origin intergenic | gDNA called gDNA | 13,530 | +947 | +1,688 |
| mixed, genic spliced | gDNA called gDNA | 0 | +236 | +236 |
| mixed, genic unspliced | RNA called gDNA | 10,092 | −534 | **+14,729** |
| all-genic multimappers | gDNA called gDNA | 215,368 | −281 | +68,180 |
| all-genic multimappers | RNA called gDNA | 16,063 | +164 | **+18,833** |
| unique (indirect, through the locus abundances) | gDNA called gDNA | 2,805,538 | −1,787 | +21,484 |
| unique (indirect) | RNA called gDNA | 294,816 | +138 | +4,591 |

Two conclusions:

1. **The mean caused arm 1's capture-ON loss.** With the sum, the genic-origin class gains 41 k instead of losing 6 k.
2. **The sum alone over-states gDNA under capture**, as §4 predicts once the efficiency is dropped: genuine RNA
   multimappers lose 33.6 k fragments to gDNA (14.7 k mixed and 18.8 k all-genic), where the mean lost about none.
   Whether the all-genic class's placements straddle a probe is not yet measured.

**The diagnostic arm's tables, on all five conditions** (`quant_accuracy.py`, fractional, the scan pinned; each arm
minus base):

| condition | tx Σ\|err\|: arm 1 / with the sum | gene Σ\|err\|: arm 1 / with the sum | gDNA estimate: arm 1 / with the sum |
|---|---|---|---|
| g00 ss.99 OFF | +135 / +1,520 | −71 / +3,411 | +71 / **+12,288** (truth 0) |
| g05 ss.50 OFF | −35 / +1,486 | +62 / +1,570 | +7,110 / +27,056 |
| g50 ss.50 OFF | −1,665 / +3,462 | −1,431 / +2,680 | +58,689 / +122,045 |
| g50 ss.99 OFF | −1,865 / −832 | −1,272 / +1,693 | +46,248 / +109,643 |
| g50 ss.99 ON | +3,839 / **−21,697** | +4,191 / **−22,910** | −7,478 / +171,243 |

The sum recovers far more gDNA everywhere and transforms the capture-ON row (transcript error −3.8 %, gene −6.7 %),
because that row's baseline leaks so much gDNA into RNA. **Off capture it costs the tables slightly, and that is not
explained.** There every efficiency is 1, so §4 says the sum is the right form, yet the tables get worse and `g00`
gains 12 k false gDNA. Candidates, none tested: the broken gDNA length law (§10, item 1) prices every placement's gDNA
term; minimap2's equal-score secondaries at paralogs, where the sum doubles gDNA's term against each transcript's
single placement; or a flaw in the derivation. This is the strongest reason to repair the length law before judging
any form here.

So neither the mean nor the sum is right alone under capture, and the sum's cost off capture is open.

## 9. The proposal: the placement-exact multimapper likelihood

The likelihood of §4, built in three steps. Each is a separate mechanism with its own falsification tests and its
own A/B, in this order.

1. **The sum over gDNA's placements** (`scoring.cpp`, `flush_mm_group`: drop `− log n_h`). Alone it over-states gDNA
   under capture (§8), so it is measured alone and lands only with step 2. Its small cost off capture (§8) must be
   explained first, on the repaired length law.
2. **Each placement at its own capture efficiency.** Every placement's term, gDNA's and each transcript's, carries
   calibration's efficiency `ĉ` for the object it sits on. The scorer runs after calibration, and every buffered hit
   carries its genomic start, so the lookup is a region index.
   - For a unique fragment the factor is common to every candidate: the change is a numeric no-op, provable with
     `rename_identity.py --check`.
   - On an uncaptured library every efficiency is exactly 1 under a one-basin landscape, so the change is again a
     no-op.
3. **The intergenic placement at a fixed rate (arm 2).** Its gDNA reading enters the unit as one more weight:
   `ρ̂·ĉ_r · f_G(L) · ½ · e`. Here `ρ̂` is calibration's bulk gDNA density and `ĉ_r` the efficiency of the placement's
   own region, so `ρ̂·ĉ_r` is fragments per start position, the unit of the in-locus terms `(θ/ℓ)·c`.
   - **The E-step** weighs it beside the locus's components.
   - **The M-step** leaves it out: it is not a component of this locus, and its density is calibration's, held
     fixed.
   - **The mass it takes is intergenic gDNA**, reported with the intergenic count. The fragment ledger counts it
     once.
   - **Pruning** treats it as any other candidate.
   - **The locus's prior** is unchanged: its odds are read over unique fragments, and its strength counts units with
     a gDNA candidate.

**Why a fixed rate answers the owner's question.** The region's density is calibration's estimate from the region's
own unique fragments, shrunk toward the library-wide landscape where the region has little evidence. Per region
where the data speak, the genome's level where they do not: the landscape decides between the owner's "each region,
each chromosome or the whole genome" from the data. Holding it fixed is the first step of the coupled EM. Its error
is second order when the region's unique fragments dominate its density, and its weakest case is a high-copy repeat
with no unique fragments (§10).

**It is also source B's EM side.** The splice-artifact design's source B (`SPLICE_ARTIFACTS_IMPLEMENTATION.md`)
finds, at index time, a remote copy of a junction's joined exons. That copy is the same object: a gDNA placement
outside the locus. Step 3 is its scoring.

**The A/B plan for each step:**

- **Fixtures first**, each verified failing on the current tree, then the fixed code broken to watch each gate fire:
  - two equal gDNA placements double gDNA's term;
  - a probed and an unprobed placement are weighed in the ratio of their efficiencies;
  - an intergenic placement is priced at its region's density and booked as intergenic gDNA.
- **The oracle ladder** has no multimappers, so every step must be bit-identical there.
- **The aligned panel**, per stratum: both errors by true origin, the transcript and gene tables.
- **VCaP:** the DNA library's multimappers with an intergenic placement, both source libraries apart.
- **The cluster's STAR-aligned ladder** before the release.

**Order against the bigger defect.** Every gDNA term above reads the gDNA length law, which is broken on aligned input
(§10, item 1). The length law is repaired first, and the three steps are measured on the repaired tree.

## 10. Known limits and open questions

1. **The aligned-input gDNA length law (the larger defect).** On the aligned `g50 ss.99 OFF` condition, calibration's
   gDNA fragment-length law puts **76 % of its mass at exactly 100 bp**, the read length (mean 99 bp). On the same
   reads' oracle BAM it has mean 217 bp. The observed pools carry a read-length spike of only 2–9 % (intronic 2.0 %,
   intergenic 8.9 %; the oracle's 0.24 %). Calibration's bulk gDNA density agrees (0.052 aligned, 0.051 oracle), but
   the EM's gDNA falls from 2.40 M to 0.96 M and nascent RNA rises from 1.02 M to 2.60 M (truth 1.01 M). VCaP's STAR
   BAM (almost every fragment a proper pair) shows no collapse: its gDNA law has mean 276 bp.
   - **Hypothesis:** single-mate fragments (minimap2 reports secondaries per mate) are taken as whole fragments whose
     length is the read's span. The length module's two-pool contrast (`calibration/fl.py`) assigns the spike to gDNA,
     and the clipping of negative mass concentrates the law on it.
   - **Derived direction:** a single mate's span is a lower bound on its fragment's length, a censored observation.
     It should not enter the length pools, and its likelihood in the scorer is the survival function
     `P(L ≥ span)`, not the pmf at the span.
   - It gets its own derivation and A/B before the steps of §9.
2. **The sum's cost off capture** (§8). The derivation says the sum is right when every efficiency is 1, yet off
   capture the sum makes the transcript and gene tables slightly worse and adds 12 k false gDNA at `g00`. Which is
   wrong: the derivation, the length law it reads, or the treatment of equal-score secondaries?
3. **The panel's N gaps.** The simulator draws gDNA from N gaps (18.7 % of the panel genome); those reads are
   unmappable. Truth is scored over all fragments, so no aligned arm can reach the gDNA truth there.
4. **High-copy repeats.** A region with no unique fragments reads the landscape's level. Is that adequate, or does a
   repeat family need the coupled EM?
5. **Alignment quality beyond NM.** Soft-clipped bases are unpriced, and the aligner's score is unread; the mismatch
   penalty α = 0.1 is underived (`ISSUES: scoring-penalties-are-underived-constants`). A worse secondary is priced only
   by its NM and its length.
6. **RNA's placements.** Each transcript keeps its best placement, not the sum. A transcript holding two copies of the
   fragment's sequence is rare, but the derivation sums there too.
7. **Calibration and the EM.** Calibration's second pass draws a held multimapper's placement by each placement's
   unique evidence. Step 3 prices the intergenic placement by calibration's own density, so the two stages would read
   the same quantity. Is the draw's hard choice compatible with the EM's soft weight here?
8. **All-intergenic multimappers** stay out of every locus and are counted as intergenic, gDNA if unspliced. Their
   spliced form (an RNA no annotated transcript holds, aligned at several places) is counted as RNA.

## 11. Files and scripts

All under `$RIGEL_RUNS/prototypes/` unless noted.

- **Arm 1:** worktree `<worktree>` at 41899ae4 plus the change (`scoring.cpp`, `bam_scanner.cpp`).
  Gates: `tests/test_pipeline_routing.py` (an intergenic hit gives gDNA its term) and
  `tests/test_multimapper_intergenic_hits.py` (the scanner buffers it beside a genic hit, and never an all-intergenic
  multimapper's); both fail on the current tree and pass on the prototype. The break test of the second gate is not
  yet run.
- **The diagnostic sum arm:** worktree `<worktree>` (arm 1 plus the sum); `2026-10-15_mm_sum/` (its
  capture-ON join, `qa_sum.jsonl`).
- **The aligned panel:** `2026-10-15_splice_a/aligned/` (`build.sh`, the suite, `qa_base.jsonl`, the base census).
- **Arm 1's scores:** `2026-10-15_mm_intergenic/` (`score.sh`, `qa_mm.jsonl`, `census/`).
- **The capture-ON join:** `2026-10-15_mm_intergenic/diag_on/` (`run.sh`, `join.py`, `read_join.py`, `nm_gap.py`, and
  `where.py`'s outcome table).
- **The length law:** `2026-10-15_mm_intergenic/diag_fl/` (`run.sh`, the oracle and aligned runs' `summary.json` and
  length tables).
