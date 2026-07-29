# Accumulator v5 — the implementation specification

**2026-07-28. Written from scratch; stands alone.** Supersedes `04_accumulator_v3.md`, which supersedes
`03_path_accumulator.md` and `02_redesign_derivation.md`. Do not cite the superseded drafts — their
retractions are folded in here.

Companion, unchanged and still correct: **`../index/00_splice_graph_design.md`** (the graph, its builder,
its invariants, its toy-GTF test matrix). This document specifies what runs *on* that graph.

Prototypes backing every number: `scratchpad/acc_proto_{c,d,e,f,g}.py`, `acc_seam_check.py`,
`acc_measure{,2,3,4}.py`, `acc_terminus_census.py`.

---

## 0. WHAT v5 IS, AND THE DECISIONS IT CLOSES

Rigel keeps the architecture it has — **nodes own contained fragments, edges own crossing fragments** —
and changes four things:

1. **TSS/TES become partition cuts.** A node becomes internally homogeneous in RNA density. 59.5 % of
   human transcript termini currently fall *inside* a node, so today's partition is wrong, not merely
   unlabelled.
2. **Splice junctions become real directed edges**, donor node → acceptor node, with their own count and
   their own exact divisor. Today a spliced fragment is split across two boundary objects with a
   region-size-dependent partial divisor chosen by ~40 lines of signature inference.
3. **Edges carry an integer COUNT, not fractional mass.** This is what dissolves the multi-region
   partitioning problem. Boundary two-sidedness, fractional mass and the `n_cross` halving all disappear.
4. **Every object stores `(count, Σ 1/L)`.** At an edge `Σ 1/L` is a divisor-free, FL-model-free estimate
   of total density; with the count it is a complete two-component deconvolution.

**Owner decisions, closed:**

| # | decision |
|---|---|
| D1 | **No taper for nascent RNA** — its support is the full pre-mRNA span. ✅ **Confirmed 2026-07-28: this applies to NASCENT ONLY. The mature taper and the junction reaches are KEPT.** |
| D2 | Isoform ambiguity in the RNA reach — **MEASURED and closed (§4.4): use MAXIMAL reach.** Worth 3.3 % of exonic bp vs the mean convention; 75.7 % of positions have no ambiguity. |
| D3 | **Unannotated junctions**: rare; deposit as unspliced if contained in one node, else as a crossing on the nodes the blocks occupy (§5.1). |
| D4 | **Junction edges carry 2 channels** (sense/antisense against the annotated motif). |
| D5 | **Scan performance may degrade; a dedicated optimization phase is budgeted** (§14, phase W6). |
| D6 | **The pooled-seam factor-2 (§10.4) is fixed as part of this work**, not chased separately. |
| D7 | FL sensitivity is a real risk; measured and tracked (§11). |
| D8 | `factor_t` movement from the finer partition is expected; watched on the gDNA benchmark (§10.5). |

**Not doing** (settled): no path/cell store · no loopy BP · no per-node FL histograms · no anchor rule ·
no fractional mass anywhere.

---

## 1. THE MODEL IN ONE PAGE

```
    NODE = a region (genomic interval).   Owns fragments CONTAINED in it.
    EDGE = a transition.                  Owns fragments CROSSING it.
             CONTIGUOUS(i, i+1)  — a genomic point; undirected in physics
             JUNCTION(donor, acceptor, strand) — an annotated intron; directed

    deposit(fragment):
        crosses 0 edges  ->  node:  count += 1,  recip += 1/L
        crosses k edges  ->  EACH:  count += 1,  recip += 1/L          (k = 0,1,2,…)

    estimator at object n, component c:      ρ_c = (its share of the count) / E_c(n)
```

Everything else in this document is (a) what `E_c(n)` is, per object class and component, and (b) why the
deposit rule needs no partitioning.

---

## 2. THE GRAPH

Specified in full in `../index/00_splice_graph_design.md`. Recap of what v5 depends on:

* Nodes tile each reference, numbered in genomic order; ids contiguous within a reference.
* Every distinct exon endpoint of every real transcript is a cut. **No signature merging** — that is what
  preserves TSS/TES. Human: **992,068 nodes** (median 145 bp), **404,168 junction edges**.
* Every edge has `src < dst`, so **genomic order is a topological order** and no BFS is needed.
* Contiguous edges carry 8 structural bits: `TSS_s`, `TES_s`, `DONOR_s`, `ACCEPTOR_s` for `s ∈ {+,−}`.
* ⭐ **New in v5 — junction edges also store `reach_donor`, `reach_acceptor`** (`int32` each): the
  transcript's cumulative exonic length either side of the junction, taken over the annotation (§4.4).
  They are needed for the junction divisor and are pure index-time arithmetic.

---

## 3. THE OBSERVATION MODEL

### 3.1 The deposit rule, at every node count

A fragment **crosses** an edge iff it spans that point. Deposit `(+1, +1/L)` on **every** edge it crosses;
if it crosses none, it is contained.

| fragment overlaps | edges crossed | deposit |
|---|---|---|
| 1 node | 0 | node: `+1, +1/L` |
| 2 nodes | 1 | that edge: `+1, +1/L` |
| 3 nodes | 2 | **each**: `+1, +1/L` |
| N nodes | N−1 | **each**: `+1, +1/L` |

`L` is the fragment's **path length** — exonic bases, introns excluded.

**Worked** (`scratchpad/acc_proto_e.py`), on `A=500 B=50 C=10 D=1 E=1000`:

| object | contained | recip | flux | recip |
|---|---|---|---|---|
| node A | 1 | 1/200 | — | — |
| nodes B–E | 0 | 0 | — | — |
| edge A\|B | — | — | 2 | 1/300 + 1/200 |
| edges B\|C, C\|D, D\|E | — | — | 1 each | 1/200 each |

### 3.2 Why no partitioning is needed

A fragment `[s, s+w)` spans point `p` iff `s ∈ [p−w, p−1]` — exactly `w` start positions. So

```
    E[ flux(p) ]  =  ρ · Σ_w f(w)·w  =  ρ · fl_mean          exactly, any FL pmf
```

**Integer observable, exact Poisson variance, independent of both flanking region sizes.** Nothing is
deposited, so nothing is partitioned: a fragment over four nodes spans three points, and each of those
points genuinely has one more fragment over it. Measured against a uniform truth over 395 seams:
`flux/fl_mean ÷ truth` = **0.997 / 1.001 / 0.987** at region lengths 2000 / 500 / 200.

Contrast today: that same 4-node fragment writes six fractional numbers whose values depend on the region
sizes. **That asymmetry is the partitioning problem, and it exists only because a mass is being
conserved.**

The **encompassed middle node** gets no contained count — correct, nothing was contained in it. It relays,
and its flanking edges' divisors do not depend on its length.

⚠ **Not double-counting, and not the "+1 per junction" bias retracted in earlier drafts.** Each edge
answers its own question with its own expectation; there is no total to conserve across edges. And this is
already what ships: `00_design.md` §6 invariant 3 — *"Flux is the integer per-side crossing count and is
not split by the 50/50 rule — only mass is."* v5 keeps flux, drops mass, adds `Σ 1/L`.

**The real consequence is CORRELATION.** `flux(e₁)` and `flux(e₂)` share fragments when the edges are
closer than a fragment length — a variance question for BP, computable as `P(w > |B|)`, present identically
today. §12.2.

### 3.3 ⭐ The reciprocal-length moment

At a point the crossing opportunity for a length-`w` fragment is *exactly* `w`, so depositing `1/w` cancels
it identically:

```
    E[ Σ 1/L ]  =  Σ_c Σ_w ρ_c·w·f_c(w)·(1/w)  =  Σ_c ρ_c  =  ρ_total
```

**Divisor-free, FL-model-free, and valid for a mixture of components with different fragment lengths.** The
count gives the second equation:

```
    N      =  ρ_g·E_g[w]  +  ρ_r·E_r[w]
    Σ 1/L  =  ρ_g         +  ρ_r
```

identified iff the two **mean** lengths differ. The design's second row is `[1, 1]`, which is what makes it
well conditioned.

Measured (`acc_proto_g.py`, gDNA FL 50, `f_g = 0.15`):

| | RNA 200 | RNA 150 | RNA 100 | RNA 70 |
|---|---|---|---|---|
| `(N, Σ L)` sd / cond | 0.028 / 1137 | 0.027 / 726 | 0.029 / 452 | 0.043 / 459 |
| **`(N, Σ 1/L)` sd / cond** | **0.017 / 283** | **0.019 / 250** | **0.023 / 250** | **0.041 / 371** |

1.6× more precise, 4× better conditioned at a 4× FL contrast; never worse. `ρ_total` from `Σ 1/L` alone,
no model: **1.0024 / 0.9989 / 1.0004 / 1.0002 ×** truth. With the fitted gDNA FL wrong by 10 %, `ρ_total`
via `Σ 1/L` stays at **1.0012× / 1.0000×** where a `Σ L`-based model gives 0.9915× / **0.9652×**.

⚠ The divisor-free identity is an **edge** property. At a node the opportunity is `max(0, B−w+1)`, which
`1/w` does not cancel; there `Σ 1/L` is simply the better-conditioned second moment (cond 242→126 at
`B` = 1000, 275→124 at 5000).

### 3.4 Conservation invariants

Fractional mass conservation is replaced by three exact integer identities:

```
    Σ_nodes contained_count  ==  n_contained_fragments
    Σ_edges flux             ==  Σ_f (# edges fragment f crosses)      computable at deposit
    Σ_nodes contained_recip  ==  Σ_{contained f} round(2^32 / L_f)
```

### 3.5 The consumption rule

By Poisson thinning, counts on distinct objects are independent Poissons. One rule binds the solver:

> **Consume a partition, never a partition and one of its own marginals.**

---

## 4. EFFECTIVE LENGTHS

All are pure functions of (graph geometry, FL pmf). None depends on the deposit rule — that is the point.
Each is unit-testable against brute-force enumeration.

### 4.1 The table

| object | component | first moment (÷ count) | second moment (÷ `Σ 1/L`) |
|---|---|---|---|
| node, contained | gDNA | `Σ_w f_g(w)·max(0, B−w+1)` | `Σ_w f_g(w)·max(0, B−w+1)/w` |
| node, contained | nascent | same with `f_r`, over `B ∩ pre-mRNA span` | same |
| node, contained | mature | same with `f_r`, over `B ∩ exons` | same |
| contiguous edge | gDNA | `fl_mean(gDNA)` | **1** |
| contiguous edge | nascent | `fl_mean(RNA)`, **0** if no strand-`s` transcript spans the seam | **1**, likewise |
| contiguous edge | mature | `fl_mean(RNA)` if mature crosses contiguously, else **0** | **1**, likewise |
| junction edge | mature | **`E_J`** (§4.2) | **1** |

The structural zeros are read from the graph's edge flags — not observed from `(SP+SN) > 0` as today.

### 4.2 ⭐ ONE CONCEPT, TWO OBJECTS — `reach`

**`reach` = how many bases of the molecule's own sequence remain, starting here.** It exists because a
molecule must fit inside the thing it came from. gDNA: the rest of the chromosome (unbounded). Mature RNA:
the rest of the transcript — hard-bounded at the polyA site. Nascent: unbounded within the gene body (D1).

It applies to both object classes, and it is the *same* quantity:

* at a **node**, `reach_M(s)` is the exonic length remaining downstream of position `s`; it makes `E_M`
  taper over the last FL bases of a transcript (§4.3);
* at a **junction edge**, the two reaches `R_d` / `R_a` are the exonic lengths remaining either side of
  the junction; they make `E_J` a trapezoid (below) instead of `fl_mean`.

*(Earlier drafts called the junction pair "runways". Same thing — the term is retired; it is `reach` in
both places, and the edge columns are `reach_donor` / `reach_acceptor`.)*

### 4.2b The junction divisor — a trapezoid

A mature fragment spans `J` with `a ≥ 1` exonic bases donor-side and `w−a ≥ 1` acceptor-side, **and both
must fit inside the transcript**:

```
    E_J = E_f[ min(w−1, R_d, R_a, R_d + R_a − w + 1)_+ ]
```

On a `N(200,50)` RNA FL:

| `R_d`/`R_a` | 5000 | 550 | 300 | 200 | 147 | 100 | 50 (vs 5000) |
|---|---|---|---|---|---|---|---|
| `E_J` | 199.0 | **199.0** | 198.2 | 160.1 | 87.8 | 19.6 | **50.0** |

`fl_mean` is exact for a typical internal junction (median human transcript 1099 bp ⇒ `R ≈ 550` each side).
It tapers hard at junctions near a transcript end. **In the owner's own worked example** — TSS at 500,
first exon `[500,550)`, junction at 550 — `R_d = 50` and `E_J = 50`, so using `fl_mean` blindly would
under-call `ρ_mature` by **4×**.

### 4.3 ⚠ D1 — no taper, and the one confirmation this leaves

**Decided: nascent RNA takes the full pre-mRNA span, no taper.** Justification: nascent is *partially
transcribed*, so it has no fixed 3′ end and there is no hard fitting constraint — only a soft 5′→3′
gradient, which is an *abundance* effect already carried by `ρ_V`.

⚠ **Mature RNA is a different case and the numbers differ.** A mature molecule genuinely ends at the polyA
site, so the constraint is hard. Measured cost of dropping it:

| | |
|---|---|
| mean `E_M/\|A\|` over ALL human exonic bp at FL 200 | **0.893** ⇒ a **10.7 % systematic gDNA over-call** |
| exonic bp in the taper zone | 16.9 % |
| `f_g` error in the last region before a TES | **+0.36** |
| `E_J` error at the owner's worked junction | **4×** |

**Recommendation: keep the taper for mature (`E_M`) and the junction exonic reaches (`E_J`); drop it only for
nascent, as decided.** One line of confirmation needed before §14 W2 starts — it changes the schema
(`exonic reach_*` on junction edges) and the calibration-time geometry.

### 4.4 ⭐ D2 — isoform ambiguity, explained

**The concept.** `reach` asks *"starting here, how many bases of this molecule's own sequence remain?"* —
because a molecule must fit inside the thing it came from. For gDNA the answer is "the rest of the
chromosome" (unbounded). For mature RNA it is "the rest of the transcript". **The problem is that a
genomic position usually belongs to more than one transcript, and they disagree.**

Worked example — one gene, two isoforms sharing an exon:

```
    T1:  exon [1000,1200)  ─ intron ─  exon [3000,3100)                 ends at 3100
    T2:  exon [1000,1200)  ─ intron ─  exon [3000,5000)                 ends at 5000

    at genomic position 1150, inside the SHARED exon:
        reach along T1 = 50 (rest of exon 1) + 100 (exon 2)        =   150
        reach along T2 = 50                  + 2000                =  2050
```

A 200 bp mature fragment starting at 1150 **cannot exist on T1** (only 150 exonic bases remain) but
**can on T2**. So the opportunity at that position is `Σ_t θ_t·F(reach_t)` — and `θ_t`, the isoform
abundance, is exactly what calibration does **not** know (the EM resolves it downstream).

**The options:**

| | rule | effect |
|---|---|---|
| **(a)** | **maximal reach** — use `max_t reach_t(s)`; a position is "open" if *any* isoform supports it | over-states opportunity where the long isoform is lowly expressed ⇒ under-states `ρ_mature` ⇒ over-states `f_g` |
| (b) | minimal reach | the mirror error, and much larger — one short isoform closes every position |
| (c) | unweighted mean over isoforms | assumes equal expression, which is wrong by ~10⁴-fold in practice |
| (d) | iterate with the EM's `θ_t` | correct, but couples calibration to quantification and breaks the one-way pipeline |

**Recommendation (a), maximal.** It is the same convention Rigel already uses everywhere else (the region
signature is a *union* over transcripts), it is monotone and stable, and its error is bounded and one-sided.
It only bites within one FL of a transcript end **and** where isoforms disagree about that end.

**MEASURED** (`scratchpad/acc_isoform_reach.py`, 1,500 sampled multi-isoform human genes, 7.82 M exonic
positions, FL `N(200,50)`):

| | |
|---|---|
| genes with ≥ 2 real isoforms | 23,240 of 39,288 |
| exonic positions where **all isoforms agree exactly** on `F(reach)` | **75.7 %** |
| positions inside *someone's* taper zone (`F < 0.99`) | 13.2 % |
| where they disagree: median / p90 / p99 of \|F_max − F_mean\| | 0.052 / 0.413 / 0.667 |
| ⭐ **mean effect of choosing MAXIMAL over MEAN, over all exonic bp** | **3.29 %** |

**Verdict: (a) MAXIMAL is safe enough to ship.** Three-quarters of exonic positions have no ambiguity at
all, and the convention is worth **3.3 %** of the mature opportunity on average — against the **10.7 %**
the taper itself is worth (§4.3). *Getting the taper in at all is the first-order decision; the convention
is second-order.*

⚠ It is not free — at the p90 of disagreeing positions the two conventions differ by 0.41 in `F`. If that
turns out to matter, there is a cheap upgrade that needs no EM coupling:

> **(e) use the MANE / canonical transcript's reach where one exists, falling back to maximal.**
> `transcripts.feather` already carries `is_mane` / `is_ccds` / `is_basic`, so this is a principled prior on
> `θ_t` — most genes have a dominant isoform — at zero extra cost. Testable against (a) with the same
> script.

---

## 5. SPECIAL CASES

### 5.1 D3 — unannotated junctions

A junction absent from the annotation has **no edge to traverse**. Owner rule, formalised:

```
    if all blocks fall in ONE node:    deposit CONTAINED, UNSPLICED channel
    else:                              deposit (+1, +1/L) on each contiguous edge the BLOCKS span,
                                       SPLICED channel; the novel junction itself is not recorded
```

Two consequences worth stating: the RNA evidence is retained at the objects the fragment actually touches
without mutating the graph; and depositing the contained case as **unspliced** preserves the 2-channel node
(§6). ⚠ Count them and report the rate as QC — a rising unannotated-junction rate means a stale annotation.

### 5.2 Ambiguous paths and multimappers

A fragment with several compatible paths — an implicit splice whose mate gap matches more than one
annotated intron, or a multimapper — goes to a **small side buffer** holding only
`(candidate object ids, channel, L)`. `FragmentBuffer` is untouched. A second pass over that buffer alone,
after the RNA FL model is fitted, assigns each fragment to its **maximum-likelihood** path.

⚠ **Assignment must stay integral.** A fractional split reintroduces exactly the non-integer observable
this design removes, and `count` would stop being a count.

### 5.3 Degenerate geometry

1 bp nodes are legal (human has 15,315) and nothing may assume `length > 1`. A fragment running off a
reference end is clipped as today. A reference with no transcripts is one node and zero edges.

---

## 6. DATA MODEL

```
NODE   (per region; 2 channels — unspliced +/−.  A spliced fragment can never be
        contained, because every intron endpoint is a cut, §5.1 excepted)
    contained_count[c]     uint32
    contained_recip[c]     uint64      Σ round(2^32 / L_f)

EDGE   (contiguous: 4 channels — a spliced fragment does cross them inside an exon
        junction:   2 channels — sense/antisense vs the annotated motif, D4)
    flux[c]                uint32
    flux_recip[c]          uint64      Σ round(2^32 / L_f)
```

**No mass. No `left`/`right`. No floats.**

**Fixed point.** `L ∈ [20, 2000]` ⇒ each term ≤ 2.1e8; at 1e8 fragments the sum is ≤ 2.1e16 against a
`uint64` ceiling of 1.8e19 — 800× headroom, ~21 significant bits at the longest `L`.
⭐ **Integer addition is associative, so the result is bit-exact at any thread count** — this *removes* the
known ~2.6 % cross-process nondeterminism rather than extending it.

**Memory** (human): node `2×12 = 24 B` × 992 K + contiguous `4×12 = 48 B` × 992 K + junction `2×12 = 24 B`
× 404 K ≈ **81 MB**, flat, no hash map, no data-dependent growth.

---

## 7. DEPOSIT

```
deposit(spans, channel):                       # spans from fragment_genomic_spans (SHIPPED, unchanged)
    L = Σ span lengths ;  R = round(2^32 / L)
    walk the spans across the node cut array (searchsorted, as region_of_pos does today):
        if all spans in ONE node n:
            contained_count[n][ch] += 1 ;  contained_recip[n][ch] += R ;  return
        for each contiguous interface the fragment spans:
            flux[e][ch] += 1 ;  flux_recip[e][ch] += R
        for each inter-span gap:
            e = junction_edge(donor, acceptor, strand)     # index-time hash
            if e is None:  -> §5.1
            flux[e][ch] += 1 ;  flux_recip[e][ch] += R
```

No slice list, no per-slice division, no `n_cross` arithmetic, no float stores.

---

## 8. WHAT THIS DELETES

`boundary_mass_left/right` · `boundary_flux_left/right` · the `mass`/`n` duality through `substrate.py` and
`node_geometry.py` · `boundary_side_eff_length`'s ½ · `spliced_side_eff_length` and its
continues/terminates selector (~40 lines reconstructing junction structure from signature bits) · the
one-sided spliced routing · `boundary_junction_strand` and its `0`-is-identity worker merge · terminal
boundary objects · `BoundarySubstrate`'s two-sided view · the `_EPS` floors that exist because divisors
collapse · `region_types`/`fl_pool_idx` plumbing through `AccumulatorSet`.

---

## 9. CONSUMERS — the complete blast radius

Grepped: 7 files import the payload, 14 read mass/flux fields. **Nothing in the resolver, scorer, EM or
locus builder.**

| file | refs | change |
|---|---|---|
| `node_geometry.py` | 28 | largest single change: per-face arrays shrink; spliced selector deleted; junction edges enter |
| `scan_payload.py` | 27 | 4 boundary matrices → node + edge tables + edge topology |
| `substrate.py` | 21 | `(counts, mass)` pair collapses; `BoundarySubstrate` two-sidedness goes |
| `bp_solver.py` | 9 | `accept_*` → structural; graft/peel move to junction edges |
| `calibrate.py` | 8 | wiring + the effective-length calls |
| `effective_length.py` | — | rewritten to §4 |
| `priors.py`, `capture_eff_length.py` | 4+3 | the seam model + §10.4 |
| `density_model.py`, `gdna_strand.py`, `strand_deconv.py`, `simplex_logodds.py`, `node_init.py`, `density_deconv.py`, `background_reference.py`, `signature.py` | 1–3 each | input shape only |
| `fl.py` | — | pools reworked, §11.4 |

**C++**: `accumulator.{h,cpp}` (563 lines) rewritten; `bam_scanner.cpp` deposit call site only —
`fragment_genomic_spans` **survives unchanged**.

**Unaffected by design**: `background_reference.py` pools intergenic regions, which gain no new cuts (an
intergenic span contains no exon endpoints).

---

## 10. THE CAPTURE EFFECTIVE-LENGTH SHRINKAGE

### 10.1 What it computes

```
                Σ_{n∈t} S_n · min(ρ_n/ρ_ref, 1)
    factor_t =  ───────────────────────────────  ,   eff_em_t = fl_t·factor_t ,
                        Σ_{n∈t} S_n                  then factor ← w·factor + (1−w),  w = C_t/(C_t+1)
```

### 10.2 Why v5 does not disturb it

* **`S_n` is pure geometry × FL** — `gdna_region_eff_len` and `gdna_boundary_len` are functions of a region
  length and an FL pmf. Identical under v5, bit for bit.
* **The seam is already ONE pooled node with ONE density** — `_pooled_seam_arrays` sums the two sides. The
  shrinkage never uses them separately, so there is no per-side quantity to lose.

### 10.3 The substitution is exact

```
    E[m_s]/S_s  =  ρ  =  E[flux_s]/fl_mean
```

Same estimand; flux counts the same fragments at unit weight and is the sufficient statistic. Contained
region densities, junction-seam imputation and `S_n` are all unchanged; the evidence term `C_t` in
`w = C/(C+1)` improves from a mass-used-as-a-count to an actual count.

Capture-off bit-identity survives trivially: under uniform gDNA every `ρ_n = ρ_ref` ⇒ `factor = 1`, **for
any `S_n`**.

### 10.4 ⚠ D6 — the pooled-seam factor-2, fixed as part of this work

`boundary_side_eff_length`'s own docstring gives `E[per-face mass] = ρ·gdna_boundary_len(L_face)`, and
`_pooled_seam_arrays` sums the two faces — so `E[seam_m] = ρ·(gbl_r + gbl_{r+1})`. But it divides by
`0.5·(gbl_r + gbl_{r+1})`, and `gbl` is **already halved** (`calibrate.py:226`). Hence `2ρ`.

Measured (`acc_seam_check.py`, uniform truth, 395 seams): shipped **1.994 / 2.002 / 1.981** at region
lengths 2000 / 500 / 200; the summed form gives 0.997 / 1.001 / 0.990.

`_gdna_region_node_arrays`'s own docstring has the **correct** formula — `S_s = ½·(E[min_r] + E[min_{r+1}])`,
which since `gbl = E[min]/2` is the **sum** of the two `gdna_boundary_len` values. The prose beside it says
"AVERAGE", and the code follows the prose.

**No test caught it** because the clip rescues the uniform case: `min(2·S_s, S_s) = S_s` ⇒ factor 1 ⇒ the
"factor-1-under-uniform bedrock" invariant passes. Under capture, a seam whose true density lies in
`(ρ_ref/2, ρ_ref)` clips and contributes **no contraction when it should contribute some**. Both consumers
are affected.

**Owner decision D6: fix it inside this work.** v5's `flux/fl_mean` reads 1.00×, so the seam density
changes by 2× at the same arm — the A/B must be scoped to expect it, and it must not be attributed to the
accumulator change.

### 10.5 What genuinely moves

`S_r = E_g[max(0, L_r − ℓ)]` collapses as regions shrink, and v8 takes the median human exon region
147 → 102 bp. Contained nodes carry less of `span_full`, seams carry more — **the weighting between classes
shifts and `factor_t` moves.** Real, measurable, belongs to the *index* arm, and must be watched on the
gDNA benchmark where the nascent siphon this module fixes is visible (D8).

---

## 11. COMPONENT ASYMMETRY — gDNA 50 vs RNA 200

### 11.1 The selection effect is real and exactly priced

`E[flux_c] = ρ_c·E_c[w]`, so an edge over-represents RNA by `200/50 = 4×` per molecule. With a truth of
`f_g = 0.150` *of molecules*, an edge sees gDNA at **4.2 % of its events**. Yet:

| object | events | gDNA share | `f_g` | bias |
|---|---|---|---|---|
| edge, thin | 178 | 0.042 | 0.148 ± 0.079 | −0.0017 |
| edge, typical | 1,775 | 0.042 | 0.148 ± 0.027 | −0.0017 |
| edge, deep | 17,751 | 0.042 | 0.150 ± 0.008 | +0.0001 |

**Unbiased at every depth. Sparsity costs variance, not bias.**

### 11.2 Node and edge measure different components

| `E_g/E_r` | edge | L=100 | L=147 | L=300 | L=1000 | L=5000 |
|---|---|---|---|---|---|---|
| | **0.25** | **115.7** | 25.5 | 2.5 | 1.19 | 1.03 |

A short node is 116× gDNA-selective (a 200 bp RNA fragment cannot fit in a 147 bp node); an edge is 4×
RNA-selective. `L = 100` shows the honest failure: 95 % gDNA events ⇒ the RNA column collapses ⇒ `f_g` bias
+0.205, sd 0.372. That node is a good `ρ_g` measurement and says nothing about `ρ_r` — **and the 2×2 Fisher
reports exactly that.**

### 11.3 What this demands

1. ⛔ **Normalisation must not move to deposit as a component correction.** `Σ 1/L` is component-blind by
   construction (§3.3) and that is why it works; a per-component divisor at deposit would require knowing
   the component, which is the estimand.
2. **The solver must carry per-component precision.** `NodeBelief` already separates
   `var_pos/var_neg/var_gdna`; `node_init` gains the 2×2 Fisher from the moment solve.
3. **FL accuracy matters asymmetrically** (D7) — where a component holds 4 % of an object's events, an
   error in the *other* component's `E_c` dominates. Re-measure at 50/200 with `f_g = 0.15` before wiring.

### 11.4 ⭐ v5 improves the FL substrate, and fixes an ordering constraint

Junction-edge fragments are **pure mature RNA by construction**; intergenic contained fragments are **pure
gDNA**. Each therefore supplies a clean, structurally-unbiased moment for its own component
(`flux / flux_recip` = the abundance-weighted mean length). This replaces the region-type × compartment
`fl_pool_mass` heuristic.

**Pipeline ordering, and it is not circular:**

```
    scan → accumulate (counts + recip; NO FL model needed)
         → fit f_g, f_r from the STRUCTURALLY PURE pools
         → compute effective lengths
         → calibrate
```

The FL fit reads only pools that are pure by structure, never by inference, so nothing is estimated from
the fragments it will later explain.

---

## 12. HONEST LIMITATIONS

1. **Crossing observables are windowed.** `E[flux(p)] = Σ_w f(w)·mean(ρ over [p−w, p−1])` — a left-looking
   average over one fragment length. Only contained counts are strictly local. This is physics, it applies
   to today's design equally, and BP is what reconciles it.
2. **Adjacent edges of a short node are correlated** — a fragment spanning a 50 bp node crosses both its
   edges. Computable as `P(w > |B|)`; must be priced, not ignored.
3. **Density cannot be resolved below the fragment length.** A 1 bp node has no independently measurable
   density and no accumulator design can give it one. Composition *is* resolvable at fine scale — it
   depends on what fragments are, not where they are.
4. **The moment solve is a Poisson design on overdispersed real counts** ⇒ over-confident. The synthetic
   suite is Poisson **by construction** and cannot test this. Gate wiring on real cfRNA.
5. **Isoform ambiguity** (§4.4) is a modelling approximation, not a solved problem.

---

## 13. TEST PLAN

Reference implementation `tests/native/_accumulator_reference.py`, byte-for-byte against the C++ — the
existing contract shape. `_fragment_spans_reference.py` and `test_fragment_spans_spec.py` are **unchanged**
(the span primitive is upstream and survives).

| # | test |
|---|---|
| A1 | contained fragment → node `+1`, `recip += R(L)`; nothing else |
| A2 | 2-node crossing → one edge `+1` |
| A3 | 4-node crossing → **three** edges, each `+1` — the §3.1 table |
| A4 | fragment fully spanning a 1 bp node → both flanking edges credited, node untouched |
| A5 | simple spliced fragment → junction edge only; **no** deposit on the contiguous edges it splices over |
| A6 | multi-junction fragment → one deposit per junction traversed |
| A7 | bookended junction (zero-length intron) → junction edge, no node skipped |
| A8 | ⭐ count conservation: `Σ contained == n_contained`, `Σ flux == Σ_f (#edges crossed)`, exact |
| A9 | ⭐ fixed-point determinism: N workers in any order → **bit-identical** `count` and `recip` |
| A10 | opposite-strand junctions at identical coordinates → distinct edges |
| A11 | unannotated junction, contained → unspliced contained (§5.1) |
| A12 | unannotated junction, spanning → contiguous-edge crossings, spliced channel |
| A13 | ambiguous fragment → side buffer, not the main table; pass 2 assigns exactly one |
| A14 | fragment running off a reference end → clipped, still deposited |
| E1–E4 | effective lengths vs brute-force enumeration: node contained, contiguous edge, junction trapezoid (§4.2 table reproduced), structural zeros |
| S1 | ⭐ uniform-density recovery: `flux/fl_mean ÷ truth ∈ [0.98, 1.02]` at region lengths 2000/500/200 |
| S2 | ⭐ `Σ 1/L ÷ ρ_total ∈ [0.98, 1.02]` with **no FL model supplied** |

---

## 14. IMPLEMENTATION PHASES

```
W0  DECIDE  §4.3 confirmation ONLY (keep the mature taper + junction exonic reaches?).
            §4.4 is CLOSED by measurement: maximal reach.
            This is the single remaining item that changes the schema.

W0b ⭐ BUILD THE TRUST INSTRUMENT FIRST — a v8→v7 node COARSENING MAP.
    v8 refines v7 (every v7 region is a union of v8 nodes), so aggregating v8 nodes back
    to v7 regions makes `z2` and the oracle bench comparable ACROSS the partition change.
    ⚠ Without this the standing gate "held-fixed z2 must not regress" is UNDEFINED from W1
    onward, because the node set itself moves. Build it before W1, not after.

W1  INDEX v8   ../index/00_splice_graph_design.md X0-X6, plus reach_donor/acceptor
               on junction edges (§2).  Gate: G1-G18, I1-I12, P1-P5; determinism.

W2  ⭐ P1g OFF THE INDEX ALONE — C1 (reframe at terminus seams), C2 (ω_graft per
    structural class), C3 (structural accept_l/accept_r).  Needs NO accumulator change.
    Gate: the pre-registered A/B of P1G_SCOPE §8 incl. its falsification test.
    ⚠ Record the baseline on the v8 partition BEFORE wiring C1, so the arm varies one thing.

W3  PYTHON REFERENCE ACCUMULATOR + tests A1-A14, E1-E4, S1-S2, written to target and failing.
    Run it on a REAL cfRNA payload through a shim before any C++ exists.

W4  C++ ACCUMULATOR + nanobind; per-worker maps; deterministic fixed-point merge.
    SHADOW MODE: emit BOTH payloads from one scan.
    Gate: byte-for-byte vs the reference; legacy payload BIT-IDENTICAL 32/32.

W5  CONSUMERS, one observable per A/B arm, each re-recording its own baseline:
      W5a  contained counts        W5b  contiguous-edge flux
      W5c  junction edges + E_J    W5d  the §4.1 effective lengths
      W5e  the (count, Σ1/L) moment into node_init   [gate: REAL cfRNA, overdispersion]
      W5f  the §10.4 seam factor-2 (D6) — its own arm, expect a 2× seam-density move

W6  ⭐ OPTIMIZATION (D5) — budgeted, not opportunistic. Profile on cfrna:LBX0190,
    NEVER the toy. Targets: the finer cut array (~8 MB, past L2), deposit throughput,
    the calibration solve at 992 K nodes (vs today's 1.5 M — it may get FASTER).

W7  DELETIONS §8, then GOLDENS LAST

Cross-cutting, from W2 onward:
  * ⭐ an `accumulator_ledger.md` appended at EVERY arm — the four entangled effects
    (partition, deposit, junction divisor, seam factor-2) are only attributable if each
    arm's measured delta is recorded as it lands.
  * ⭐ EXTRACT, do not mutate, `node_geometry.py` (28 refs, the largest change): put the
    v5 per-object geometry in a NEW module so both paths stay runnable during migration.
  * keep a v7 index alongside the v8 one until W7 — v8 indexes are unreadable by v7 code,
    so shadow mode alone does not give a rollback for the INDEX arm.
  * `pass0_oracle_bench.py` needs a v5 object mapping; build it in W3 with the reference.
```

**Standing gates on every arm:** ruff · full suite · held-fixed `z2` must not regress · one thing varied ·
a pre-registered falsification test · baseline re-recorded from the current tree in the same session at
both refits.

⚠ **Shadow mode (W4) is the key de-risking device**: both accumulators fed by the same fragments in one
scan means every A/B is exact and free, and any arm can be reverted without a rescan.

---

## 15. RISKS

| risk | mitigation |
|---|---|
| **Scan slows on the finer partition** | accepted (D5); W6 is budgeted, not hoped for |
| Four entangled effects (partition, deposit, junction divisor, seam factor-2) all move the numbers | W2/W5a–f isolate them into separate arms; shadow mode makes each exact |
| The moment solve is over-confident on overdispersed data | W5e gated on real cfRNA, not the toy |
| Isoform ambiguity biases `E_M` and `E_J` | §4.4 measurement in W0; maximal-reach error is bounded and one-sided |
| Index rebuild invalidates every cache | `/tmp/rigel_selfsolve` namespaced by index hash in W1 |
| The toy validates nothing here | it overstates contained efficiency **4.6×**, understates multi-crossing **3.8×**, has no alternative TSS/TES and is Poisson by construction. **Gate on real cfRNA.** |
