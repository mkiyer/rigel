# The accumulator — design

**Draft 2.** The specification for the replacement accumulator: the per-fragment tally built during the
single-pass BAM scan. Everything downstream reads only its output.

Draft 1 was reviewed by six independent agents against the real cfRNA data; seven findings were blocking
and are folded in here. Where a number in draft 1 was wrong it is corrected in place and the correction
is noted, because a design document that quietly changes its numbers is how this project accumulated
seven contradictions in the first place.

**Notation.** ⭐ = measured this session, script named. Companion: `docs/SESSION_HANDOFF.md` (measured facts,
equations, traps).

> ⚠ **How to read the measured numbers.** They come from four cached cfRNA samples, three of them
> unusually sparse and all hybrid-capture enriched. **They are observations about those samples, not
> design constraints.** Rigel targets the >2 million human RNA-seq datasets in the SRA: depth,
> composition, fragment-length distributions, capture state and multimapper rates all vary by orders of
> magnitude across them. A number here is evidence that a mechanism *exists* and that the arithmetic is
> right — never a parameter to fit to. **Nothing in the deposit rule depends on any of them.**

---

## 1. Why replace it

The current accumulator chops each fragment at partition borders into *slices*, then decides "this slice
crosses its right border if another slice follows it". That test is a property of the slice **list**, not
of the genome, so it **cannot distinguish a contiguous crossing from a splice jump** — both are just
"there is a next slice".

⭐ Consequence on real cfRNA: **18–22 % of all "spliced" mass sits at positions with no annotated donor
and no annotated acceptor**, and `boundary_junction_strand` — whose comment asserts "this boundary is a
junction" — is set at every position a spliced fragment merely touches. Reproduced bit-for-bit in both
the Python reference and the C++.

No flag can repair this: by the time a consumer sees the tally, the two populations are already merged,
so a structural predicate can only subtract. ⭐ Measured: the best available flag substitution made the
benchmark *worse* (0.046675 → 0.046887) while deleting 34–41 % of genuine mature-RNA signal on cfRNA.
The fix has to be in the tally.

## 2. The model, in one page

> **The genome is a graph. Nodes are intervals; edges connect them. A fragment is a path.
> Nodes count fragments *contained* and *spanning*; edges count fragments *crossing*.
> Every object stores an integer count and a density term.**

**An edge is a 0-bp line.** No width. A fragment crosses it iff the line lies **strictly inside one of
the path's segments** — bases on both sides, adjacent in the molecule.

**Two kinds of edge.**

| | contiguous | junction |
|---|---|---|
| what it is | the seam between genomically adjacent nodes | a directed donor→acceptor splice link |
| direction | bidirectional | directed, single-stranded |
| composition | gDNA + RNA | **pure RNA** — gDNA cannot be spliced |
| endpoints | **implicit**: edge `i` lies between node `i` and node `i+1` | explicit `(src, dst, strand)` |
| count (human) | 1,043,595 | 404,168 |

**A splice jump deposits on its junction edge only** — never on the contiguous edges it splices over.
That is the whole point: "continued into the intron" and "spliced away" become two counters, and the
question calibration exists to answer becomes a subtraction instead of an inference.

**Worked example.** Exons `(0,1000)` and `(2000,3000)`; nodes `[0,1000)`, `[1000,2000)`, `[2000,3000)`;
contiguous edges at 1000 and 2000.

| fragment | covers | result |
|---|---|---|
| `(800,1000)` | 800…999 | contained in node 0 — ends *at* the line, does not cross |
| `(801,1001)` | 801…1000 | crosses contiguous edge @1000 (has bases 999 and 1000) |
| `(1000,1200)` | 1000…1199 | contained in node 1 (intronic — gDNA or nascent RNA) |
| blocks `(900,1000)`+`(2000,2100)` | 900…999, 2000…2099 | **junction edge 1000→2000 only.** Crosses no contiguous edge. |

**Two components, always: gDNA and RNA.** There is no mature/nascent distinction anywhere in the
accumulator. Spliced fragments measure RNA directly; unspliced fragments are the mixture to deconvolve.
Nascent RNA is modelled exactly as a single-exon transcript — no special case, no separate channel.

## 3. What a fragment is

This section did not exist in draft 1, and its absence was the single largest defect.

### 3.1 The molecule and its length

A fragment is a **path**: its aligned blocks, joined across *mate gaps* (unsequenced but part of the
molecule) and broken by *introns* (excised, not part of the molecule).

```
    span  =  rightmost block end  −  leftmost block start
    L     =  span  −  Σ (intron lengths)
```

**One definition.** What varies is not the definition but how an intron is *detected*:

| intron detection | source | trust |
|---|---|---|
| **observed** | a CIGAR `N` operation | high |
| **implicit** | `SPLICE_IMPLICIT` — an intron sits in the unsequenced mate gap, and the fragment length is otherwise incompatible with the fitted length distribution | **lower — mark it** |

⚠ An implicit splice must be flagged (`sj_implicit`), and a fragment carrying one **must not enter the
pure-RNA pool** (§8): its splice is a model inference, not an observation, so certifying it as RNA would
make the pool depend on the thing it is used to estimate.

⭐ **ONE WORD FOR IT: `sj_implicit`** (owner ruling, 2026-07-29). The type is already `SPLICE_IMPLICIT`
and `sj_` is already the codebase's prefix for a splice-junction attribute (`sj_strand`, `sj.feather`,
`sj_map`, `SJKey`, `sj_lookup_into`), so a `FragmentPath` carries `sj_strand` and `sj_implicit` — two
attributes of one splice under one prefix. The earlier `introns_inferred` / `inferred_intron_fragments`
were a **second** word for a concept that already had one, which is the mistake `sj_strand` was renamed
to fix. *Inferred* is gone from the vocabulary.

⛔ **`sj_implicit` is NOT a doubt about the path.** It says the splice was not sequenced. A fragment
whose candidate transcripts imply **different intron sets** has no determined `L` at all — that is
`path_ambiguous` (§9), a rejection rather than a flag on a deposit.

⚠ `L` must **not** be the sum of aligned block lengths. That drops the mate gap, which is part of the
molecule, and collapses the length distribution to a spike at twice the read length.

**Whatever counts toward `L` must also count as coverage for crossing**, or the density estimator is
biased. Mate gaps therefore cross lines; introns never do. This is forced, not a preference.

### 3.2 ⭐ Bounded influence — the fragment length limit

Draft 1 removed the old design's per-fragment mass budget of 1.0 and put nothing in its place. On real
data a few pathological alignments then own the entire edge channel.

⭐ Measured on LBX0190 (155,352 read-name groups with a mapped block, `bam_spans.py`):

| | measured |
|---|---|
| median `L` | **117 bp**; p90 377; p99 17,722 |
| p99.9 `L` | **76,895,968 bp**; max **234,545,096** (a whole chromosome) |
| groups with span > 1 Mb | 505 (0.33 %) — **95 % of them contain a supplementary record** |
| groups with blocks on more than one reference | 9,926 |
| **top 1,000 groups' share of all edge crossings, unbounded** | **99.8 %** |

**The rule: drop any fragment whose `L` exceeds `--max-fragment-length` (default 1000).** This is an
existing user-facing parameter with a documented default, not a new internal constant.

⭐ Its effect, and the reason it must be applied to `L` and not to the span:

| limit applied to | fragments dropped | top-1,000 share of crossings | mean crossings / fragment |
|---|---|---|---|
| nothing | — | 99.8 % | 16.5 |
| **`L` ≤ 1000** | **5.45 %** | **4.16 %** | **0.93** |
| span ≤ 1000 | **37.96 %** ⛔ | — | — |

A 300 bp molecule spanning a 10 kb intron has a 10 kb span. Limiting the **span** would discard every
legitimately spliced fragment; limiting **`L`** discards 5.45 %, of which 25 % carry a supplementary
record. Same number, opposite outcome.

**Chimeras are the other half, and they belong in the scanner, not here.** 95 % of the multi-megabase
spans carry a supplementary record, i.e. a split alignment grouped into one molecule. Two shipped
mechanisms let them through: `group_records_by_hit` files a `BAM_FSUPPLEMENTARY` record with its primary
(`bam_scanner.cpp:853-858`), and `detect_chimera` returns `CHIMERA_NONE` unless at least two blocks carry
a non-empty transcript set. A split alignment is a chimera — trans-splicing, a structural variant, or a
long unannotated intron (STAR only calls introns below ~300 kb). It should be classified as such, or
treated as `SPLICE_IMPLICIT`, and in neither case should it deposit a 1 Mb path.

### 3.3 The artifact catalogue

Every case a real BAM contains, and what `L` is for it. Silence here is how a spec becomes folklore.

| case | rule |
|---|---|
| **soft clips** (`S`) | **dropped — never counted in `L`** (owner ruling). The molecule may extend past the alignment, so `L` is understated on clipped fragments; this is accepted, not corrected |
| **deletion** (`D`) | advances the reference and **is** part of `L` — it is a real base of the template |
| **CIGAR `N`** | an intron; removed from `L`. **Every `N` is an intron — there is no minimum length** (owner ruling, §12). A floor would be a magic number with no evidence behind it, and it is what ships today |
| **overlapping introns on one fragment** | ⭐ **a real case, not a hypothetical** — the scanner reads `XS`/`ts` once per *record*, so a pair whose mates disagree about an acceptor yields two overlapping introns for one molecule (MO_3021 `chr8:138206290-138206943`, 1 group in 875,670). **Merged.** `L` is then *defined* as the total length of the path's segments, so no second formula for it can exist — a naive `span − Σ intron` disagrees by up to 1.5× and goes **negative** on a wide overlap |
| **abutting introns** | also **merged**, for a different and stronger reason: they imply a zero-length exon, which is physically impossible (owner ruling, §12) — a transcript with one is molecularly identical to a transcript without it, so no single molecule can legitimately use both. The index cannot produce the pair either: a zero-length exon is dropped when the exon arrays are built, fusing its two flanking introns into one |
| **overlapping mates** | ⭐ **the dominant case by far — 84.5 % / 89.4 % of proper pairs overlap or abut** (median gap −94 / −108 bp with 150 bp reads). The union of blocks handles it; a `TLEN`-based span would not. Draft 1 said 44 %, which is what `\|TLEN\|` reports and is wrong |
| **mate on another reference** | chimera; does not deposit |
| **supplementary / split** | chimera (§3.2); does not deposit as one molecule |
| **secondary** (`NH > 1`) | multimapper; side-buffered, §9 |
| **PCR duplicates** | the model assumes independent counts, so duplicates must be removed **upstream**. The accumulator does not de-duplicate; it must **assert** and report |
| **fragment overhanging a reference end** | clipped to the reference; `L` is the clipped length, so the placement count stays consistent |
| **blacklisted junction** | does not deposit at all — its true span is unrecoverable |

## 4. The deposit rule

### 4.1 The rule

```
for each accepted fragment, with path P and length L:
    if P lies entirely within one node        ->  node.contained   += (1, 1/L)
    else:
        for each contiguous line strictly inside a segment of P:
                                                  edge.crossing    += (1, 1/(L-1))
        for each intron of P matching an annotated junction:
                                                  junction.crossing += (1, 1/(L-1))
        for each node whose BOTH lines are crossed by P:
                                                  node.spanning += (1, 1/L)
```

**No partition. Every crossed edge receives the full weight**, and a fragment crossing several junctions
credits **every one of them** — each edge owns its own expectation, and the strand model is fitted from a
separate scan output (`pipeline.py:325`), so nothing downstream is distorted by crediting all of them.

**Unannotated introns create no junction edge and deposit nothing across the gap.** The blocks deposit
on the contiguous lines they span, in the **unspliced** channel, so an unannotated junction competes with
gDNA rather than being certified as RNA — the safe direction, since unannotated junctions are
disproportionately artifactual. If the whole path fits inside one node it is a contained unspliced
fragment. ⭐ Rate: **2.3 % / 1.7 %** of `N` operations are unannotated, but **18.9 % / 10.8 % of
*distinct* junctions** are — the latter is the honest "novel junction" figure and is an order of
magnitude larger.

### 4.2 ⭐ Why there is no partition, and why that is not "unfair to short fragments"

In the old design the observable was a *mass* — a conserved 1.0 spread over objects — so splitting was
mandatory. Here it is not, because of the cancellation that makes the estimator work:

```
    P(a length-L fragment crosses a given line)  ∝  L
    what it deposits there                       =  1/L
    expected contribution                        =  constant, INDEPENDENT of L
```

Every fragment length contributes **equally** to each edge's estimate. A long fragment appears at more
edges because it genuinely measured more places — more *evidence*, not more *weight*.

Dividing by `K` (the number of edges crossed) destroys the cancellation: with node spacing `s`,
`K ≈ L/s`, so the contribution becomes `s/L` — dependent on both fragment length and node spacing.

⭐ Measured (`partition_test.py`; uniform density, exact truth, estimate ÷ truth at every edge):

| node spacing | **no partition** | `1/K` | overlap-weighted |
|---|---|---|---|
| 50 bp | **0.9948** | 0.2806 | 0.2333 |
| 100 bp | **0.9944** | 0.5440 | 0.3790 |
| 200 bp | **0.9938** | 0.9072 | 0.4862 |
| 400 bp | **0.9960** | 0.9960 | 0.4980 |

`1/K` reads up to **3.6× low**. It agrees with the correct rule only where nodes are coarser than
fragments — i.e. where there is nothing to partition. The human median node is 151 bp. Overlap-weighting
is a flat factor-2 error and never converges. *(The `1/K` error is `E[1/K]`, so it depends on the
library's length distribution as well as on node spacing — which only strengthens the conclusion.)*

**The honest framing:** an edge count is not a share of a conserved total. There *is* no total. Density
is intensive, like temperature, not extensive, like mass.

### 4.3 ⭐ The one-base correction

A 0-bp line means "bases on both sides", which is `L − 1` admissible placements, not `L`:

```
    E[ Σ 1/L      over fragments crossing a line ]  =  ρ · (1 − E[1/L])     ← 0.54 % low
    E[ Σ 1/(L−1)  over fragments crossing a line ]  =  ρ                    ← exact
```

⭐ Verified on identical draws so the noise cancels: `Σ1/L` = 0.020012 vs `Σ1/(L−1)` = 0.020126, a ratio
of exactly `1 − E[1/L]`, predicted to four digits. The 0.994 column in §4.2 is this correction, not
noise. (The 4000 bp row's 0.9830 *is* noise — 90 interior edges.)

So: **`1/(L−1)` at an edge, `1/L` at a node.** Strictly the edge identity is `ρ · P(L ≥ 2)`; a fragment
of length 1 cannot cross anything.

⚠ *Correction to draft 1:* it justified this by analogy with the contained divisor's `+1`, claiming the
`+1` makes the divisor zero when the node is one fragment long. That is inverted — at `ℓ = w`,
`(ℓ − w + 1)_+ = 1`, and *dropping* the `+1` gives 0. The `+1` is simply the discrete count of start
positions.

## 5. What each object stores

Every channel is a pair: an **integer count** and a **fixed-point density**. Both, always — the count
carries statistical power (a Beta-Binomial needs an integer), the density carries level. **No fractional
quantity ever enters a likelihood.**

### 5.1 ⭐ ONE STRAND CONVENTION — the channel axis

**Every count and density array has exactly two channels, indexed by the fragment's own GENOME strand.**

```
STRAND_COLUMNS = {Strand.POS: 0, Strand.NEG: 1}   # without exception, in every bank
```

**Sense and antisense are the transcript-relative notion, and they are DERIVED, never stored.** A
junction edge carries its own genomic strand, so a consumer that wants transcript-relative counts
computes `sense = (fragment strand == junction strand)`. Nothing is lost — including the
wrong-stranded-junction aligner diagnostic, which is the antisense column of that derivation.

⚠ *This replaces draft 2's original table, which stored some banks by genome strand and others by
sense.* One schema holding two conventions is what made B2 possible: `node_spanning` was indexed
sense-relative for spliced fragments and genome-relative for unspliced ones. ⭐ Measured on real cfRNA,
**65–69 % of all `node_spanning` deposits come from spliced fragments, and 40–44 % of those landed in the
opposite column from the unspliced fragments beside them.** The authority is
`tests/native/_accumulator_reference.py`'s `genome_channel`.

⚠ **What the genome strand of a fragment is, exactly.** `align_strand` is the OR of the fragment's exon
blocks' strands, where read 1 keeps its reference orientation and read 2 is flipped to agree
(`bam_scanner.cpp:661-686`). So the channel is the fragment's **read-1 genome orientation** — a fixed,
protocol-independent quantity. Whether it equals the *RNA molecule's* strand is a property of the library
protocol (on a dUTP first-strand library read 1 is antisense), which is why strandedness is a declared
header field and not a channel convention. See §5.2.

**Node** — two disjoint populations, each two channels.

| population | definition | channels |
|---|---|---|
| **contained** | the whole path lies inside the node | {genome +, genome −} |
| **spanning** | one segment of the path crosses **both** of the node's lines | {genome +, genome −} |

There is **no spliced/unspliced split on a node.** ⭐ Every annotated intron has both endpoints as
partition cuts, so an annotated-spliced fragment always crosses its junction edge and can never be
*contained*; and an *unannotated* junction is routed to the unspliced population by design (§4.1). But
a spliced fragment's blocks routinely **span** nodes, and those deposits land in the same two channels as
everyone else's — which is precisely why the node banks must be genome-strand.

**Why spanning is stored:** it brackets the length distribution from above where containment
brackets it from below (§6), which is what gives short nodes any leverage at all. *(It was originally
proposed to price the correlation between adjacent edges. That justification is **withdrawn** — the
solver does not use correlation, and ⭐ a per-node counter would price it badly anyway: it captures only
the adjacent pair, which is 44–55 % of the off-diagonal covariance mass, and the density channel's
covariance needs `Σ 1/(L−1)²`, not a count.)*

**Contiguous edge** — fragments whose path crosses the line.

| channel | why |
|---|---|
| **unspliced** × {genome +, genome −} | the mixture being deconvolved |
| **spliced** × {genome +, genome −} | a *spliced* fragment crossing a contiguous line is a **certified RNA crossing** — gDNA cannot be spliced. The cleanest RNA marker at a seam, and the current design merges it away |
| structural flags (8 bits, from the index) | TSS / TES / DONOR / ACCEPTOR × strand |
| reach (4 columns, from the index) | `lo` / `hi` per strand (§7) |

"Spliced" here is a property of the **fragment** — it used at least one annotated junction somewhere —
not of the line it crossed. A fragment deposits into one of the two populations, never both.

⚠ It is certified **RNA**, not certified **mature** RNA — a partially spliced nascent molecule carries
junctions too. That matters because §7 must not apply a mature-only taper to this channel.

**Junction edge** — fragments using that exact donor→acceptor jump.

| channel | why |
|---|---|
| {genome +, genome −} | the same convention as everywhere else. The edge knows its own strand, so the column that disagrees with it **is** the antisense count — a free aligner-error diagnostic, derived rather than stored |
| reach (from the index) | ⚠ *Correction to draft 1*, which said junction edges need no reach. §7's placement formula is one formula for both edge kinds and requires it. A large minority of junctions sit close to a transcript end, where ignoring the taper is a multi-fold divisor error |

No unspliced channel and no structural flags on a junction: it is spliced by definition and is not a
genomic position with terminus semantics.

⚠ **`sj_strand` no longer selects a channel.** Its only remaining job is to disambiguate the
(constructible, biologically impossible) case of two annotated junctions sharing a coordinate pair and
differing only in strand. The `spliced_primary` branch that used to route on it is gone.

### 5.2 The payload

Draft 1 had no schema, which made "decide before the schema freezes" vacuous. The payload is:

```
nodes        contained_count[n,2]   contained_density[n,2]
             spanning_count[n,2]    spanning_density[n,2]
             start_count[n]                          # §10.2, the invariant
contiguous   unspliced_count[e,2]   unspliced_density[e,2]
             spliced_count[e,2]     spliced_density[e,2]
junction     count[j,2]             density[j,2]
pools        length_histogram[pool, max_len]         # integer, §8
qc           the denominators of §10.3
header       format_version · strandedness · max_fragment_length · graph_provenance
```

The trailing `2` is the §5.1 channel axis — genome strand — in every array without exception. `count` is
`uint32`, `density` is `uint64` fixed point (§10.1); the two are never given one column name, because
their divisors differ (`L` at a node, `L − 1` at a line). The authoritative field list is
`docs/accumulator/IMPLEMENTATION_PLAN.md` §3.5.

**Provenance is not optional.** ⚠ The current cache key covers `nodes.feather` only; on 2026-07-29 an
edge fix rewrote every `edges.feather` while leaving every `nodes.feather` byte-identical, so a stale
cache would have verified clean. The new payload is **edge-keyed by construction**, so its provenance
must cover nodes **and** edges.

**Strandedness is a declared field**, and §5.1's single convention is what makes it *only* a field.

⭐ A real bug found by the review, in the shipped code: `primary = (align_strand == motif_strand)`
(`bam_scanner.cpp:1493`) with `align_strand` derived from read 1. On a dUTP first-strand library read 1 is
antisense, so the channel **labelled** *sense* holds **0.6 %** of spliced fragments. This is the
unexplained κ ≈ 0.002–0.012 recorded across all four cfRNA libraries — a convention bug, not biology. The
gDNA/RNA split was unaffected (the model is symmetric in the strand convention), but the ± assignment was
inverted, and nothing could see it.

⭐ **The new schema cannot express that bug**: no channel is labelled *sense*, so no label can be wrong.
What remains is a genuine question about the *library*, one level up — does read-1 orientation agree with
the RNA molecule's strand? — and it is answered once, in the header, not per fragment. **Measure it,
declare it, assert it in QC.** The estimator is free: agreement of the read-1 strand with the
splice-motif strand over spliced fragments.

## 6. What the numbers estimate

Let component `c ∈ {gDNA, RNA}` have start-density `ρ_c` and length distribution `f_c`. At any object,
both observables are linear in `ρ_c`, with coefficients depending only on `f_c` and geometry:

```
    count    =  Σ_c ρ_c · A_c        A_c = E_c[ placements ]
    density  =  Σ_c ρ_c · B_c        B_c = E_c[ placements / weight ]
```

**At an edge, far from any transcript end** (`placements = w − 1`, `weight = w − 1`):

```
    count    =  ρ_g·(μ_g − 1)  +  ρ_r·(μ_r − 1)
    density  =  ρ_g            +  ρ_r
```

A 2×2 system with determinant `μ_g − μ_r` — **solvable whenever the two mean lengths differ.**

> ⚠ **Correction to draft 1, and it is the most important one.** Draft 1 said "every edge yields a
> gDNA/RNA split from fragment length alone". **That is a library-scale identity sold as a per-object
> one**, and it is wrong by three to four orders of magnitude.
>
> ⭐ Measured, building the design's own path rule from the BAM (`ck/cross2.py`): LBX0190 has **135,704
> clean crossings** over 1,043,595 contiguous edges. **98.66 % of seams carry zero crossings** and the
> 99th percentile is **one**. To reach `sd(f_g) = 0.05` at a single edge needs **1,122 / 415 / 1,707 /
> 425** crossings (Monte-Carlo validated: empirical sd 0.0496–0.0501, unbiased) — a threshold met by
> **12 seams out of 1,043,595**. Even a 10× looser target is met by 0.02 %.
>
> ⭐ And there is a hard ceiling, which is the number that should govern expectations: LBX0190's
> **entire genome-wide crossing budget supports at most ~120 independent `sd = 0.05` estimates**, and
> nearer **~54** counting distinct fragments (53 % of clean fragments cross no seam at all). Pooling
> adjacent seams is inefficient because crossings are concentrated — 183 seams hold half of them.
> **The edge pair is a library- or chromosome-scale statistic, plus a per-object statistic on the few
> hundred highest-flux seams, and nothing in between.**
>
> An under-populated object must therefore **emit nothing at zero precision** — never a floored
> division. That is a standing trap in this project: a "no data" default of 100 % gDNA was actively
> seeding false gDNA into neighbouring exons.
>
> ⚠ Two caveats on the information numbers. They are **strongly dependent on the true composition**
> (LBX0588 ranges 0.95 → 7.79 across `f_g` 0.1 → 0.9) and are quoted at an **assumed** `f_g = 0.15`;
> and crossings per edge is a **depth** statistic, not a design property — 0.130 on LBX0190 but 0.504
> and 0.489 on LBX0588 and MO_3021.

**No assumption about which component is longer.** ⭐ Across the four real libraries, gDNA mean lengths
are 146 / 102 / 157 / 309 and RNA 227 / 270 / 206 / 246 — **gDNA is shorter on three and longer on one.**
Discriminability is a function of `|μ_g − μ_r|` and the object's geometry, with no direction. Draft 1's
worked example assumed gDNA long and tight, which inverts on real data, and that table is removed.

**At a node there are THREE equations**, because containment and spanning bracket the length
distribution from opposite sides:

```
    contained     placements = max(0, ℓ − w + 1)   weight = w    fits inside      w ≤ ℓ
    spanning  placements = max(0, w − ℓ − 1)   weight = w    spans both ends  w ≥ ℓ + 2
```

Both are exact start-position counts — a node `[a,b)` is spanned iff the start lies in `[b−w+1, a−1]`,
which is `w − ℓ − 1` positions. The populations are disjoint, and `w = ℓ + 1` is in neither.

⭐ Verified against brute force (`spanning.py`, RNA ~ N(200,60), ρ = 0.05):

| node length | contained obs / pred | spanning obs / pred |
|---|---|---|
| 25 bp | 0.0000 / 0.0007 | **8.7213 / 8.7055** |
| 100 bp | 0.0573 / 0.0594 | 5.0087 / 5.0141 |
| 200 bp | 1.2100 / 1.2177 | 1.1373 / 1.1725 |
| 600 bp | **20.0813 / 20.0452** | 0.0000 / 0.0000 |

**Short nodes are measured by what spans them, long nodes by what fits inside them, and the crossover is
the mean fragment length.** No node scale is blind — though as `ℓ → 0` the spanning count converges
to the flanking edges' crossing count and stops being *independent* information.

**⚠ `Σ1/L` is model-free only at an edge.** At a node the `1/w` does not cancel `(ℓ−w+1)_+`; it is a
better-conditioned second moment and nothing more. ⭐ This matters more than it sounds: median `L` is
117 bp against a median node of 151 bp, and mean crossings per fragment after the length limit is 0.93 —
so **on real cfRNA most fragments are contained, i.e. most of the data lands in the frame that still
needs the length model.**

## 7. Where the estimators are valid

**Transcript ends.** An RNA fragment must *fit* in what remains. Past that point the admissible
placements stop being `L−1` and become the distance to the end — the same for every length — so the
weight cancels nothing:

| distance to transcript end | 20 bp | 50 | 100 | 200 | 300 | interior |
|---|---|---|---|---|---|---|
| `Σ1/(L−1)` reads | **0.108 ρ** | 0.270 ρ | 0.533 ρ | 0.923 ρ | 0.999 ρ | **1.0000 ρ** |

⭐ Analytic, confirmed by simulation on a 1099 bp transcript (the human median) to three decimals.

One formula covers both edge kinds and both frames, with `R_lo`/`R_hi` the molecule's own remaining
sequence either side:

```
    placements(w, R_lo, R_hi)  =  max( 0,  min( w − 1,  R_lo,  R_hi,  R_lo + R_hi − w + 1 ) )
```

Mean fragment length is its large-reach limit, not a separate case.

**The taper is per-component and comes from the annotation, not from a fit.** gDNA's template is the
chromosome, so `taper_g = 1`. RNA's template ends where its transcript ends.

### 7.1 What "reach" means when RNA is RNA

Nascent RNA is a single-exon transcript spanning the gene, so it tapers at **both** ends like any other
transcript. Reach is a **maximum over transcripts**, and genomic distance is never shorter than exonic
distance, so:

**RNA reach at a position = the genomic distance to the ends of the widest transcript covering it, per
strand.**

Two consequences, the second of which is the one that matters:

* inside an **exon** the nascent span sets the reach, not the exonic offset;
* inside an **intron**, mature RNA has zero reach but nascent RNA reads straight through, so **RNA reach
  in an intron is not zero**. A mature-only reach would declare zero RNA opportunity across every intron
  in the genome — exactly backwards, since the intron is where nascent RNA lives.

⚠ *Correction to draft 1:* it claimed this required a filter divergence (cuts from real transcripts,
reach including synthetic nascent rows). **It does not.** Synthetic nascent spans are shadow transcripts
sharing their parent's exact span, so they add no boundary information and no new rows are needed —
reach is computed from the same real transcripts' genomic spans. **One filter, `~is_synthetic`,
everywhere.** Draft 1 introduced precisely the filter divergence that caused this project's worst index
bug.

**One length distribution for RNA.** Fragment length is set by fragmentation and size selection, not by
whether the molecule was spliced.

### 7.2 Density steps, i.e. hybrid capture

Hybrid capture creates a density step at an exon boundary. This must be modelled, not treated as an
edge case.

⚠ **First, a number draft 1 asserted and could not support.** It claimed capture is "~1000×". ⭐ Measured
from the cached payloads, the pooled annotated:intergenic density ratio is **16.5 / 6.7 / 5.1 / 44.5**;
1000× appears only at the 99th–99.9th percentile of annotated node density. And that ratio **confounds
capture with expression** — an exon is dense both because it was captured and because transcript is
there, which is the very thing the deconvolution exists to separate. **The enrichment magnitude is not
established, and no design number should be conditioned on 1000×.**

⚠ **Second, and this is a design-level gap: the sign and the side of the smear depend on an assumption
nobody has stated.** "Density step" can mean three different things, and they disagree:

| reading | where the bias lands |
|---|---|
| **start-density** — a fragment is enriched iff its *start* base is in the target | positive, on the **enriched** side, peaking within one fragment length |
| **overlap-any** — a fragment is pulled down if its span *touches* a probe (the actual physics of hybrid capture) | **exactly zero inside the enriched region**; the entire smear lands on the **depleted** side |
| span-mean | a third answer again |

⭐ All three agree far from the step, so **no uniform-density test can settle this** — it must be
declared. The overlap-any reading is the one that matches how capture works, and under it the enriched
side is clean and the depleted side carries all the distortion.

⚠ **This is a calibration question, not an accumulator one**, and it is listed as such in §12. Which
reading holds determines which side of a boundary a *solver* should trust; the accumulator counts
fragments identically either way, and no line of the deposit rule branches on it.

⭐ Under the start-density reading, at a 2× step the gDNA-fraction bias peaks at **+0.033 / +0.035 /
+0.014 / −0.010** for the four libraries at `f_g = 0.15` (draft 1's +0.039 came from a noisy
single-realisation simulation and is replaced). Note vcap's **sign flips**, because its gDNA is longer
than its RNA. The peak sits **within one fragment length** (57–205 bp, sample-specific), and the bias
scales roughly as `f_g(1−f_g)`, so **any quoted figure must carry its `f_g`.**

**The contained population does not smear** — a contained fragment lies entirely inside one node, and
exon boundaries are node boundaries here, so a capture step falls on a node edge and never inside one:

| observable | localisation | length model |
|---|---|---|
| **contained** (and spanning) | **sharp at node resolution** | **required** |
| **edge crossing** | smeared over ~1 fragment length | not required |

⚠ **But it is not free, and draft 1 implied it was.** ⭐ Two measured costs:

1. **The sharpness is entirely a property of the cut.** Move a node so it straddles the enrichment
   boundary by half its width and the contained bias jumps to **+0.285**. So exon-edge cuts are
   *load-bearing*, not merely tidy — and a missed cut is a first-class failure, which is precisely what
   the v8 partition exists to prevent.
2. **Contained is not length-model-free.** At a 150 bp node the containment probability differs **6.6×**
   between gDNA and RNA, and a 10 % error in the fitted length models costs **0.010–0.026** of
   composition (0.081 for vcap). §8's corrections are therefore load-bearing too.

The two observables are genuinely complementary — the crossing channel gives a model-free level where
density is flat, the contained channel gives a sharp answer at a capture edge — and storing both is what
lets the solver have it both ways.

Density below one fragment length remains unresolvable by any design. That is physics, not a defect —
but the benchmark suite must contain a **step**, or none of this is exercised.

## 8. Fragment-length models

**Five pools, each structurally pure** (settled; `docs/accumulator/IMPLEMENTATION_PLAN.md` §3.4 has the enum). Purity is
the whole point — it is what stops a length model from being fitted to the fragments it will later explain.

| pool | population | ⭐ measured mean `L`, LBX0190 |
|---|---|---|
| `DNA_INTERGENIC` | contained in an intergenic node | **88.0** |
| `DNA_INTRONIC` | contained in an intronic node | **88.5** |
| `DNA_INTRON_EXON` | crossing one line with {intron, exon} flanks — a "splash" read | 138.8 |
| `DNA_INTERGENIC_EXON` | crossing one line with {intergenic, exon} flanks — a splash read | 211.6 |
| `RNA_SPLICED` | using an annotated junction, splice **observed** | **220.7** |

⭐ The two pure gDNA pools agree with each other at ~88 bp and the pure RNA pool sits at ~221 bp — a
**2.5× separation, measured with no fitted model anywhere**, and `μ_g − μ_r` is the determinant of the 2×2
in §6. MO_3021 gives 118.4 / 124.1 against 203.3, a 1.7× separation: real, and sample-specific.

Why these five and no others:

* **junction fragments are pure RNA** — gDNA cannot be spliced. ⭐ Mature RNA crosses an exon↔intron seam
  0 times in 1,146 seams. ⚠ Excluding `sj_implicit` fragments (§3.1), whose splice was never sequenced —
  certifying those as RNA would make the pool depend on the model it is used to fit.
* **intergenic and intronic contained fragments are pure gDNA** — ~99 % accurate on real data, the
  standing assumption of the tool.
* **the two splash pools are the only *on-target* gDNA population** (§8.2), so they are named rather than
  folded in. ⭐ They land *between* the two extremes (139, 212), which is the expected signature of
  on-target gDNA plus leaked RNA.
* **there is deliberately no pool** for an exonic contained fragment or a multi-line crossing. Those are
  gDNA/RNA mixtures, and an impure pool is worse than a missing one.

⚠ **This is not a re-labelling of what ships.** ⭐ The shipped `gdna_fl_pmf` has mean **146.05** on
LBX0190 because it sums four differently-tilted pools (intergenic and intronic × contained and boundary);
the pure intergenic-contained pool says **88.0**, and an independent measurement during the code recon put
the corrected pure estimate at 86.49. The new axis moves the gDNA mean length by **~40 %**, and the old
value was biased long by pooling the boundary compartments.

```
    scan and accumulate  (needs no length model)
      → fit f_gDNA and f_RNA from the pure pools
        → placements / effective lengths
          → calibrate
```

Nothing is ever estimated from the fragments it will later explain.

### 8.1 ⚠ Two corrections, both of which draft 1 got wrong

**(a) The mean-length estimator is frame-specific, and at an edge it has an exact closed form.**
⭐ At an edge, `Σ 1/(L−1) = ρ` identically, so

```
    count / Σ(1/(L−1))  =  E[L] − 1      exactly.      →      E[L] = count / Σ(1/(L−1)) + 1
```

Verified to **+0.0001 %** on all four libraries. Draft 1 omitted the `+1` and was low by exactly
`1/E[L]`. At a **contained node** the opportunity is `(ℓ−w+1)_+` and the weight `1/w`, so the same ratio
converges to `1/E[1/w]` — the **harmonic** mean, biased low. And the only gDNA pool is intergenic
*contained* fragments, i.e. exactly the frame where the naive recipe fails.

**(b) Every pool histogram is length-biased by its own opportunity, and draft 1 had the direction
backwards.** ⭐ `rna_fl_pmf` **is** the raw junction pool (`fl.py` builds it from the spliced histogram,
one bin per fragment), so truth is obtained by **dividing** by the opportunity, not multiplying. Draft 1
described tilting it further, which corresponds to nothing in the pipeline.

⚠ And the opportunity for the junction pool is **not** `w−1`. A fragment enters that pool if it crosses
*at least one* junction anywhere in its transcript, which is a transcript-level count over all real
transcripts — a different and larger quantity.

⭐ Corrected magnitudes: the RNA contained-effective-length at the median 151 bp node is understated by
**1.36× / 1.41× / 1.35× / 3.02×** (production values 8.96 / 10.65 / 11.85 / 0.31 against corrected
12.17 / 14.97 / 15.93 / 0.93) — not the "~2×" draft 1 claimed, and **vcap is 3×**, not 2×.

**The fix is the design's own reciprocal trick, one level up: divide each pooled histogram by its own
opportunity count before normalising** — `placements(w)` for a crossing pool, `(ℓ−w+1)_+` for a
contained pool, the transcript-level count for the junction pool. Then the mean is read from a corrected
`f_c` rather than from a ratio, and every frame agrees. This feeds `μ_g − μ_r`, the determinant of §6,
and (per §7.2) a 10 % error here costs 0.010–0.026 of composition — so it is load-bearing, not hygiene.

⚠ ⭐ Finally, `gdna_fl_pmf` is **four differently-tilted pools summed** (intergenic and intronic ×
contained and boundary), whose raw per-pool means differ by ~9× on one library. It is not a truth to
compare against. The clean gDNA estimate is the **intergenic-contained pool corrected by its own
opportunity** — and because intergenic nodes are ~11.5 kb, that correction moves it only 0.2 %.

### 8.2 The capture caveat

Under hybrid capture the intergenic pool is **depleted, not impure** — composition stays clean, coverage
may not. ⭐ On-target gDNA fragments are ~42 bp shorter than off-target, so a model fitted off-target is
mis-centred for exactly the fragments that leak. This is why the two splash pools exist as **named pools**
(§8) rather than being folded in — the comparison against the intergenic pool becomes a QC output instead
of a guess.

⚠ The intergenic pool is also thinned by the multimapper gate (§9), non-randomly toward repeats. Both
effects push the same way; the QC comparison is what will surface them.

## 9. Deferred: ambiguous paths and multimappers

**Deliberately after the main accumulator is built and verified** — but the schema must leave room now.

A fragment compatible with several paths (ambiguous junction assignment) or several loci (`NH > 1`) goes
to a **side buffer** of `(candidate object ids, channel, L)`. After the first pass, when the length models
and per-path densities exist, each is assigned by **multinomial probability across the matching paths**.

### 9.1 ⭐ The third population: an implicit splice whose intron set is not determined (owner ruling, 2026-07-29)

`SPLICE_IMPLICIT` fragments **do deposit** — an implicit splice overlaps an annotated intron and matches
in every other way, so the only thing missing is that the splice motif was not sequenced. `sj_strand` is
then inferred from the implying transcript, which is reasonable **iff there is only one of them**.

⛔ **When several candidate transcripts imply DIFFERENT INTRON SETS the path is not determined, and
nothing about the fragment can be tallied in the first pass.** The implied intron set fixes `L`, both
density quanta, the pool bin, the segment list, and therefore *which lines the fragment crosses*. There is
no partial answer to fall back on:

* it cannot deposit **spliced**, because which junction it used is the unknown;
* it cannot deposit **unspliced** either — `L` "involves an intron", so it does not fit the fragment-length
  distribution unless one of the candidate introns is cut out, which is the very choice in doubt. Only
  breaking the molecule into its aligned blocks and depositing each separately would avoid committing,
  and that discards the mate gap `L` depends on;
* and it must not be *forced* — picking one compatible transcript is choosing an `L` at random.

⭐ **The test is on distinct intron SETS, per candidate transcript** — not per gap, and not on `L`. Per gap
is wrong because a union of per-gap matches can take gap A's intron from transcript T1 and gap B's from
T2, producing a path **no single molecule has**; and `L` is too weak a key, because two introns of equal
length at different positions give the same `L` but different segments and hence different crossings.
⚠ The **empty** set counts as a distinct hypothesis: if some candidates imply an intron and others do not
(a retained-intron isoform), the fragment is ambiguous between a spliced `L` and an unspliced one.

So the rule is: **deposit iff the candidate transcripts yield exactly one distinct implied intron set;
otherwise the fragment is `path_ambiguous`** — it deposits on no object, leaves `start_count` untouched,
and is counted on its own denominator until the side buffer drains it. The second pass can then resolve it
exactly as it resolves a multimapper, and with more to go on: **the fragment length and the strand both
discriminate between the candidate transcripts.**

⛔ **The candidates are the REAL transcripts — `~is_synthetic`, alone.** Every gene carries a synthetic
nascent shadow transcript spanning its locus as one exon, so that row implies no intron in any gap, ever.
Counting it as a dissenting candidate makes the test unsatisfiable: ⭐ measured, it deferred **100 % of
implicit fragments on all three real cfRNA libraries and deposited none.** It also conflates two different
questions — *which intron did this molecule splice?*, which is what the test is for, with *was this molecule
nascent?*, which is a **component** and one the accumulator deliberately does not decide. ⚠ The filter is
`~is_synthetic` and never `is_nrna`: on a non-synthetic row that flag means "single-exon, so mature ≡
nascent", and using it as a realness filter has already deleted the termini of 26,475 real transcripts once
(`docs/SESSION_HANDOFF.md` §3 trap 3).

⭐ **The set test is the WHOLE test — it already implies the candidates agree on strand**, so there is no
second condition and no code for one. If every candidate implies the same intron at every gap, then every
candidate implies the *whole* set, so the first candidate that matches a gap is the same transcript at
every gap, and the strand is that one transcript's. A separate strand check could therefore never fire.
(An earlier statement of this ruling carried one; it was redundant, and a condition that cannot fire is a
line a later reader has to disprove.)

⭐ And the ruling *deleted* a special case rather than adding one. A strand-coincident implied intron
defers with everything else, so `contradictory_sj_strand` keeps its original and narrow meaning — the
**mates' motif tags disagreed about an observed splice** — and is unreachable for an implicit one.

⚠ **The assignment must be integral** — one fragment, one path. A fractional `1/NH` split re-creates the
non-integer observable this redesign exists to delete, and is provably biased for any second-moment
statistic.

**Recommended form, and it needs no maximum-likelihood machinery.** Group the buffered fragments by
`(candidate set, channel)`, then apportion each group's integer count across its candidates in proportion
to their densities, rounding by **largest remainder**. That is deterministic (no seed, no RNG),
integral, single-pass, unbiased in aggregate, and free of the winner-take-all bias that argmax would
introduce at a 51/49 locus. Sampling would be unbiased per fragment but would make the accumulator
non-reproducible; iterated ML would add a layer this tool does not need.

Until this phase lands, multimappers and `path_ambiguous` fragments do not deposit — and the accumulator
must **emit both denominators**,
because the loss is not uniform. ⭐ Measured multimapper share of intergenic fragments: **59.3 %
(LBX0190) / 39.5 % (MO_3021) / 20.2 % (vcap) / 16.7 % (LBX0588)** — library-dependent, with LBX0190 (by
far the smallest library) the extreme. Intronic is **~53.6 %**. Either way the gate thins the pure gDNA
pool non-randomly, toward repeats.

## 10. Arithmetic, invariants, QC

### 10.1 Determinism

**Counts:** `uint32`. **Densities:** `uint64` fixed point, `+= round(2^32 / weight)`, where `weight` is
`L` at a node and `L − 1` at an edge (two different quantities — do not give them one column name).

**Why not `float32`, which would halve the density storage.** The choice is not about precision, it is
about **reproducibility** — and on the numbers `float32` loses on both counts:

| | `float32` | `uint64` fixed point |
|---|---|---|
| storage, all density cells at human scale | 36.6 MB | 73.3 MB |
| …as a share of the tool's measured 8.6 GB peak | — | **+0.43 %** |
| relative error per deposit (L = 40 … 1000) | eps 1.2e-7 | ⭐ **3.7e-9 … 6.9e-8** — *better* |
| serial vs 8-way parallel sum, 10⁶ deposits | ⭐ **8.8e-5 relative** | **0 — bit-identical** |

⭐ Integer addition is associative, so the result is identical at any thread count, on any machine. That
matters more than it looks: the *current* float accumulator differs by only ~3.7e-7 per cell across
worker counts, and that propagates to a **~2.6 % difference in the calibration output** — the same BAM
gives different answers on different machines. Fixed point removes it entirely, and 37 MB is not a price
worth paying to keep it.

⚠ `uint32` fixed point is not viable: the scale needed to keep the quantisation error negligible
(`2^32`) overflows a 32-bit accumulator well below realistic per-object depth.

⚠ *Correction to draft 1:* the headroom argument was backwards. `round(2^32/weight)` is largest at
**small** `L`, so the binding constraint is a **minimum** length, not the maximum. At `L = 2` the weight
is 2^32 and 10⁸ fragments reach 4.3e17 against a ceiling of 1.8e19 — fine, but the precondition is
`L ≥ some floor`, which the shipped minimum read length already provides. State it as a floor.

The rounding mode must be pinned exactly (`llround`, ties away from zero) or byte-identity is undefined.

### 10.2 ⚠ An invariant that can actually fire

Draft 1's three "conservation identities" are **tautologies** — each right-hand side can only be
evaluated by re-running the deposit. The review's replay satisfied all three exactly while 91 % of the
crossings were junk.

**The replacement: one `uint32` per node counting fragments whose START lies in it.**

```
    Σ_nodes start_count  ==  number of accepted fragments        ← exact, externally checkable
```

It is 4 MB at human scale, it is the same counter a start-anchored diagnostic would need, and unlike the
crossing identities it can be checked against a number the scanner knows independently.

*(Start-anchored **attribution** as the primary estimator was considered and rejected: it keeps one point
and discards the path, which removes crossing flux, junction usage, and the spliced-vs-unspliced contrast
at boundaries — the structural gDNA/RNA signal. The count is kept as an invariant, not as an estimator.)*

### 10.3 QC the accumulator must emit

Not optional, and not derivable afterwards: **dropped multimappers · chimeras · fragments over the length
limit · blacklisted junctions · strand-undefined fragments · side-buffered ambiguous fragments ·
unannotated-junction rate · implicit-splice (`sj_implicit`) rate · **ambiguous-path fragments** ·
measured strandedness · duplicate assertion**. Every conservation statement must name its denominator.

### 10.4 The test matrix

Written before the code, each a named case: contained → node only · 4-node crossing → exactly 3 edges ·
fragment spanning a 1 bp node → both flanking edges **and** that node's spanning counter, contained
untouched · **spliced → junction edge only, no deposit on the contiguous edges it splices over** ·
zero-length intron · paired-end mate gap crossing a line → *does* count · `sj_implicit` marked and
excluded from the RNA pool · **`path_ambiguous` → no deposit anywhere, its own QC denominator, and
`start_count` untouched** · unannotated junction → unspliced channel · fragment over the length limit →
no deposit anywhere, QC incremented · clipping at a reference end · single-node reference · opposite-strand
junctions at identical coordinates are distinct edges · N workers → bit-identical · uniform density
recovered within [0.98, 1.02] · `Σ1/(L−1)` within [0.98, 1.02] **with no length model supplied**.

**Validate by re-deriving with a different algorithm.** Both bugs found in the index work were caught
this way and by nothing else. A validator that calls the builder's own helper validates nothing.

**The oracle is the production accumulator itself, partitioned by true fragment origin** — split the
simulated BAM by read-name origin, run the *same* scanner on each part, assert the parts sum to the whole.
Nothing is reimplemented, so nothing can be subtly wrong.

## 11. What gets built first

1. **The Python reference accumulator** (~300 readable lines) and the §10.4 matrix, tests written first
   and failing. This is the executable specification; the C++ is gated on byte-identity to it.
2. **Run it on a real cfRNA payload through a shim** before any C++ exists. The 10 Mb synthetic suite
   cannot validate this design: ⭐ its fragment length variance is **zero**, so nothing length-dependent —
   which is now the core of the design — is exercised at all.
3. **The benchmark suite** (owner's plan): real human genes as a backbone, cache the whole BAM-scan
   intermediate, then piggyback a synthetic stress chromosome. It must contain a **density step**,
   **fragment-length variance**, **alternative TSS/TES inside exons**, and a **low-gDNA × strong-capture
   corner**, or the design's own failure modes are invisible to it.

## 12. Decisions, and what actually blocks implementation

### Settled

| | ruling |
|---|---|
| **objects** | nodes (contained + spanning), contiguous edges, junction edges |
| **`L`** | leftmost block start → rightmost block end, minus introns. Soft clips dropped |
| **crossing** | a 0-bp line strictly inside a segment of the path |
| **partition** | none — every crossed edge gets the full weight |
| **weights** | `1/(L−1)` at an edge, `1/L` at a node |
| **fragment length limit** | `--max-fragment-length`, default 1000, applied to `L` |
| **unannotated junctions** | deposit **unspliced**, and **not along the intron** |
| **junctions per fragment** | credit **every** one |
| **components** | two: gDNA and RNA. Nascent RNA is a single-exon transcript |
| **storage** | counts `uint32`; densities `uint64` fixed point `round(2^32/weight)` |
| **spanning** | store count and density |
| **multimappers / ambiguous paths** | side buffer, deterministic largest-remainder apportionment, **a later phase** |
| **implicit splices** | ⭐ **they deposit**, with `sj_strand` inferred from the implying transcript and `sj_implicit` set, **iff the candidates imply exactly one intron set and agree on strand**. Otherwise `path_ambiguous` → side buffer (§9.1) |
| **optimisation** | aggressive, after implementation and accuracy |
| **strand** | ⭐ **ONE convention: genome strand, everywhere** (§5.1). Sense/antisense is derived, never stored |
| **minimum intron length** | **none — every CIGAR `N` is an intron** (A, below) |
| **malformed introns** | overlapping **and abutting** introns both merge, and the merge is counted (A, below) |

### A. The two rulings that closed the last open item

**A.1 There is no minimum intron length. Every CIGAR `N` is an intron.** This was listed as the single
item blocking implementation — the deposit branches on whether a gap is an intron, so a `few-bp N` that is
not a splice would be mis-handled. **Owner ruling: no floor.** A fixed threshold would be a magic number
with no evidence behind it, and "shorter than the shortest annotated intron on this reference" makes the
deposit rule depend on the annotation in a second, hidden place. It is also what ships today, so the
ruling changes no line of code — the design was the only thing out of step.

**A.2 A zero-length exon is physically impossible, so two ABUTTING introns are malformed and merge.**
There is no molecular difference between a transcript carrying a zero-length exon and one without it, so
a single molecule can never legitimately use two introns that share an endpoint — the pair is an alignment
artifact. Overlapping introns merge for the plainer reason that they are contradictory observations of one
molecule (§3.3). Both are counted in `introns_absorbed` rather than dropped silently.

⚠ A.2 **reverses** a strict-overlap decision taken during the S2 review, which had merged only genuine
overlaps on the theory that abutting introns were two legitimate splices. `_normalise_introns` merges on
`start <= previous_end`.

**Nothing now blocks implementation.**

### These are NOT accumulator questions — they belong to calibration

Listed so they are not re-raised as blockers. **None of them changes the deposit rule, the schema, or
a single line of the accumulator.** They are all about what a *consumer* does with the tally.

* **The capture model** — whether enrichment selects on a fragment's start base or on span-overlap
  determines which side of a boundary a *solver* should trust. The accumulator counts fragments either
  way.
* **Enrichment magnitude.** A property of a library, not of the tool.
* **How much composition information a single edge carries, and at what pooling scale it becomes
  usable.** Entirely a function of depth and composition, both of which vary by orders of magnitude
  across the SRA. The accumulator's job is to record the counts; deciding when a count is enough is the
  solver's.
* **Whether the spanning density term earns its bytes.** Stored; measure later on real workloads.
* **QC thresholds.** The accumulator must *emit* the denominators (§10.3). What counts as alarming is a
  downstream judgement.
* **The junction pool's opportunity function** (§8) — a fragment-length-model question, downstream of
  the tally.
