# RIGEL — CARRY FORWARD

All that survives from 278 docs / 74,823 lines + 142 memory files.

> **⭐ THIS FILE IS NOT DEAD — §2 and §3 are the most-used reference in the project.** §3's traps caught
> real defects repeatedly during the 2026-07-30 session (trap 1 shaped every validator, trap 15 is now
> `scripts/design/suite_resolves.py`, trap 16 stopped a single premature profile, trap 20 is why the
> suite carries 92 ERCC references, trap 25 decided what the scan cache stores, trap 27 merged two
> implementations of the capture landscape). §2 is the derivation reference the code depends on.
>
> **Read order: §3 (traps) → §2 (equations) → §1 (facts, re-derive before quoting).**
>
> ⚠ **The parenthetical citations are PROVENANCE, not links.** `PATH_MARGINALIZATION.md`,
> `METHODS.md`, `CALIBRATION_ARCHITECTURE.md`, `SIMULATOR.md`, `05_v5.md`, `scratchpad/acc_proto_g.py`
> and the rest were deleted in the 2026-07-29 purge and the 2026-07-30 cleanup. They record *where a
> number was originally measured*; git history is the only place they still exist. Do not chase them —
> re-derive instead, which is what `scripts/design/index_census.py` is for.
>
> ⚠ **§0 is the least current section.** Most of C1–C7 were settled by `ACCUMULATOR_DESIGN.md`; they are
> kept because the *measurements* in them are real and were never re-recorded elsewhere. Treat §0 as
> "measurements attached to questions", not as a live decision queue — the live queue is `TODO.md`.
> ⛔ **§4 was deleted 2026-07-30** (see the note at the foot of this file).

**Dead vocabulary:** *pass-0* = the first solve, before any fitted prior. *psi solve* = picking a node's gDNA/sense/antisense split on a
grid. *relay / reframe* = passing a density estimate to a neighbour, rescaled. *graft / peel* = adding/subtracting a junction's RNA when
moving an estimate across a seam. *P1d / omega_graft / P1e* = two fudge variance terms, both patches for the missing "is this seam a gene
end?" bit; both die with it. *count-zero-information* = a fragment count says nothing about gDNA-vs-RNA, only how precisely you know it.
*pin bug* = a message computed from the destination's own answer. *signature* = a 4-bit exon±/intron± node tag. *mwae* = mass-weighted
mean absolute error of estimated gDNA fraction vs simulated truth.

---

## §0 CONTRADICTIONS — the owner must decide

**C1 (biggest). Is `+1/L` per fragment actually unbiased?** Two measurements disagree; the design hangs on it.
- FOR: at a point `E[Σ1/L] = ρ` exactly for any fragment-length mix, no divisor, no FL model. Measured 1.0024 / 0.9989 / 1.0004 / 1.0002 ×
  truth. (`05_accumulator_v5.md:144`, `scratchpad/acc_proto_g.py`)
- AGAINST 1: with *spatially varying* ρ, `Σ1/L` returns a left-looking length-weighted average of ρ over one fragment length. gDNA smooths
  over ~60 bp, RNA over ~200 bp, so at a density STEP the gDNA fraction is over-called by **+0.155** (100 bp nodes) to **+0.565** (25 bp
  nodes). Exon boundaries under capture ARE density steps. Anchoring each fragment to its first node and dividing by that node's
  admissible-start count has bias +0.0000 and costs 1.4–2.7× variance. (`PATH_MARGINALIZATION.md:47-113`)
- AGAINST 2: near a transcript end `Σ1/L` reads **0.103×** truth (0.992× mid-transcript) — the crossing count falls 90 % while mean observed
  length moves only −5.6 %. A terminus removes *placements*, not *lengths*. **`Σ1/L` fixes the length bias; the reach taper fixes the
  placement loss; neither substitutes for the other.** (`accumulator_ledger.md:196-220`)
- Decide: `1/L` **plus** a reach taper, or anchor attribution. Not equivalent.

**C2. `Σ1/L` is fractional; the only intrinsic gDNA/RNA signal needs an integer.** The strand likelihood is a Beta-Binomial whose power is
the integer count, so a fractional density term re-creates the count-vs-mass duality the redesign exists to delete
(`PATH_MARGINALIZATION.md:95-113`). Likely resolution: the design stores BOTH — count for statistical power, `Σ1/L` for level, and no
fractional quantity ever enters a likelihood. Say so explicitly.

**C3. Is intergenic a clean pure-gDNA fragment-length pool?** YES: intergenic contained fragments are pure gDNA, so `count / Σ(1/L)` is the
gDNA mean length with no circularity (`05_v5.md:539-556`). NO: under capture intergenic is not captured — composition-pure but
coverage-empty — and on-target gDNA fragments are ~42 bp SHORTER (on-exon mean 343 vs off-target 386), so a model trained off-target is
mis-centred exactly for the fragments that leak; the signal lives in exon↔intron "splash" reads (`06_implementation_plan.md:131-141`).

**C4. How much gDNA is in the real cfRNA libraries?** Memory + calibration docs say LBX0190 ≈15 %, vcap ≈20–25 %. Current code measures
LBX0190 **0.0562**, MO_3021 0.1264, LBX0588 0.7718. No ground truth exists. **← the gDNA half is still open.**
- ✅ **THE SENSE-FRACTION HALF IS ANSWERED (S5.f, 2026-07-30): a CONVENTION MIRROR, and the INFERENCE IS CORRECT.** The suspicion was that
  "sense fraction 0.002–0.012 on all four real libraries (nearly fully antisense)" was a read-orientation bug. The chr22 pilot has ground
  truth: a library simulated at `strand_specificity = 0.99` **calibrates to κ = 0.0101** — `1 − truth`, on both capture arms (at 0.50 it reads
  0.4990–0.5002, where a mirror is invisible by construction).
  ⭐ **But the mirror is CONSISTENT, so it cancels.** Measured on a **zero-gDNA** stranded condition (truth `f_gdna = 0` exactly) by injecting
  κ: the fitted 0.0100 gives `f_gdna` **0.0030**; forcing the nominal truth 0.99 gives **0.4992** — **166× worse**, and 6× worse than forcing
  the uninformative 0.50 (0.0792). κ and the per-node sense columns are mirrored the same way, the strand likelihood scores a mirrored
  observation against a mirrored `p`, and the deconvolution is right. **Only the exported scalar is mis-labelled.**
  ⚠ **So fact 17's κ column below reads backwards**: the four real libraries are ordinary **highly stranded** libraries, not near-purely
  antisense ones. ⛔ Any fix must preserve `f_gdna = 0.0030` while moving κ to 0.99 — a change that moves κ alone has broken the inference to
  fix a label. `TODO.md` §6. (Originally `SESSION_2026_07_30_HANDOFF.md:333-341` vs memory `cfrna_sample_characteristics`)

**C5. Is `+1` on every crossed edge double counting?** Raised and retracted — each edge has its own expectation, so it is unbiased per edge.
But 12.97 % of human exonic crossing fragments cross more than one interface, and a fragment spanning a short node crosses both its edges,
so the counts are *correlated*. Overlap = `P(w > node length)`; with median node 151 bp vs RNA ~200 bp this fires for the MAJORITY of nodes.
Formula stated, never priced — **the largest unpriced risk in the new design.** (`05_v5.md:129,565`)

**C6. "No length taper for nascent RNA" is asserted, never measured.** Its mature counterpart was measured (11.0 %); the nascent case rests
entirely on "partial transcription ⇒ no fixed 3′ end". It sets the schema.

**C7. "The new solve may get faster" compares incompatible counts.** Today: 752 K regions + 753 K seams = 1.5 M objects. New: 1,043,881
nodes + 1,447,763 edges ≈ 2.5 M. Like-for-like the object count **doubles**.

**Resolved, do not re-litigate:** "59.5 % of termini" is WRONG → **53.4 %**; "992,068 nodes / median 145 bp" is stale → **1,043,881 /
151 bp**; "positions that are both a gene end and a splice site are the majority" is wrong → **0.99 %**.

---

## §1 MEASURED FACTS

1. ✅ **RE-DERIVED 2026-07-30 from the rebuilt index and every number reproduced exactly** — run
   `scripts/design/index_census.py INDEX --gtf GTF` rather than quoting this. It describes an ANNOTATION
   (GENCODE v46 / Ensembl 112 with ERCC controls), not the tool.
   **The human splice graph (index format v8):** 1,043,881 nodes (median 151 bp, **mean 2,970** genome-wide,
   15,687 of length 1, **591,442 = 56.7 % shorter than one 200 bp RNA fragment**), 1,043,595 contiguous
   edges, 404,168 junction edges = 1,447,763 edges. Replaces 752,654 merged regions, +38.7 %.
   254,319 transcripts parsed, 26,475 single-exon. (`splice_graph.py:16-19`)
   ⚠ The often-quoted **mean node 2,552 bp is chr1 only**; genome-wide it is 2,970.
2. **53.4 % of real human transcript ends fall strictly inside an old merged region** and were invisible to
   it. This is the entire reason for the redesign. ✅ Re-derived: **232,451**, = **53.42 %**.
   ⚠ The denominator was quoted as 435,291; it re-derives as **435,107** terminus-flagged contiguous seams.
   The numerator and the percentage are exact; the denominator depends on how a terminus is counted.
3. **Seam census over 1,043,595 contiguous edges:** terminus-only 40.70 %, splice-site-only 58.31 %, BOTH 0.99 % (10,337), neither 0.00 %. ✅ Re-derived exactly.
4. **The crossing-count estimator is unbiased and region-size-independent:** `count / mean_FL ÷ truth` = 0.997 / 1.001 / 0.987 at region
   lengths 2000 / 500 / 200 over 395 seams. This single measurement justifies the whole deposit rule. (`scratchpad/acc_seam_check.py`)
5. ⛔ **CORRECTED 2026-07-30 — DOES NOT GENERALISE, and the original claim is withdrawn.** It was
   measured at ONE setting (gDNA 50 / RNA 200, a 4× mean separation, deep in the region where the
   `μ_g − μ_r` determinant is large). Swept over 756 gDNA×RNA length pairs × 3 shape families × 4
   compositions (`scripts/design/observable_efficiency.py`) it is **false at every node ≥ 250 bp and at
   the edge**, and exactly reversed when the two components share a mean length — where `(count, Σ1/L)`
   carries **zero** information at any depth. Storing BOTH is what S5.a shipped. The original text:
   **`(count, Σ1/L)` beats `(count, ΣL)` everywhere:** sd 0.017 vs 0.028, condition number 283 vs 1137 (gDNA FL 50, RNA 200, f_g 0.15);
   1.6× more precise, 4× better conditioned, never worse. With the fitted gDNA length wrong by 10 %, total density reads 1.0012× truth vs
   0.9652×. (`scratchpad/acc_proto_g.py`)
6. ⛔ **CORRECTED 2026-07-30 — THE GEOMETRY REPRODUCES; THE "11.0 % gDNA OVER-CALL" DOES NOT.** The original text is kept at the foot of
   this item. What it measured is real; what it was *labelled* is not, and the label was load-bearing (it is why A7 was deferred past S5.f
   as a step expected to buy 11 %).
   **The geometry, re-measured on the chr22 pilot index at the suite's own FL:** the position-weighted taper ratio is **0.8738**, against the
   original 0.8904 on a different index, annotation and FL. ✅ **That half reproduces.**
   ⛔ **But the estimator is FRAGMENT-weighted, not POSITION-weighted, and that is the whole story.** Weighting each line by the unspliced
   mass actually on it, the taper ratio is **0.9596** — a 4.0 % geometric effect, not 11.0 %. Fragments concentrate mid-transcript, which is
   exactly where the taper is inert: **89.2 % of edge mass sits on lines the taper does not touch at all**, and only 2.7 % on lines tapered
   below 0.25×. A bp-weighted mean gives every terminal position equal say with every mid-transcript position; the calibration does not.
   ⛔ **And a geometric bias in a divisor does not pass through to composition 1:1.** Measured end to end (A/B on all 8 pilot conditions,
   arm A bit-identical to the S5.f baseline): turning the taper on moves the library gDNA fraction by **≤ 0.0002**, and *toward* truth on all
   four zero-gDNA conditions by less than 1e-4. Node gDNA mass moves by 0.021 % in total, edge gDNA mass by 0.52 %.
   ⭐ **The "+0.36" IS real for an individual node — and irrelevant.** Max per-node |Δf_gdna| = **0.3961**, so the original claim is right
   about the worst node. But only **6 nodes** move by more than 0.30, and they carry **17 fragments** between them = **0.0002 %** of node
   mass. Both statements are true at once; only the second decides anything.
   ⚠ **The condition under which this null would FAIL, and the one-number test for it.** The null holds because mass is mid-transcript. A
   3′-biased or heavily degraded library — cfRNA is degraded — or transcripts short relative to the fragment length would put mass exactly
   where the taper bites. **The screening test is the mass-weighted taper ratio** (0.9596 here): compute it before assuming the taper is
   negligible on a new library. It is far cheaper than wiring the taper.
   ⚠ **Contiguous seams remain WORSE than junctions (0.750 vs 0.886)** — a junction has exons both sides by construction; a contiguous seam
   is often inside a terminal exon. That comparison was never the disputed part.
   *Original text, kept because the geometry in it is sound:* "Ignoring the fragment-fit taper near a transcript end over-calls gDNA by
   11.0 % — bp-weighted mean of (crossing effective length ÷ mean FL) = 0.8904 over 3,000 genes / 10,957,350 exonic positions, RNA N(200,50).
   … The gDNA fraction is off by +0.36 in the last region before a polyA site.
7. **Isoform disagreement about "how much transcript remains" is small enough to ignore:** 75.7 % of exonic positions have every isoform
   agreeing exactly; max-over-isoforms vs mean moves the answer 3.29 %; taking the max independently per side of a junction agrees with the
   best realisable single-isoform pair on 93.9 % of disagreeing junctions (mean ratio 0.9989). 23,240 of 39,288 genes have ≥2 isoforms.
8. **Nodes and edges measure different components.** gDNA/RNA opportunity ratio: **0.25 at a crossing point** (4× RNA-selective) vs
   **115.7 at a 100 bp node**, 25.5 at 147 bp, 2.5 at 300 bp, 1.19 at 1000 bp. A 200 bp RNA fragment cannot fit in a 147 bp node, so a short
   node is a good gDNA measurement and says nothing about RNA. Carry per-component precision, not one scalar.
9. **Contained effective length collapses to EXACTLY 0.0 on 12.4 % of fine-partition nodes** (1.6 % of old merged regions): with RNA
   N(200,50), E = 0.000 at length 1 and 10, 0.030 at 59. Effective lengths are **superadditive** — over 118,195 splitting regions
   `Σ E(children)/E(whole)` = 0.7652, and 0.0917 for 305 bp → 145+160. So densities, effective lengths and variances cannot be pooled across
   two partitions; only truth pools additively (`f_g = ΣG/(ΣG+ΣR)`).
10. **The finer partition costs +9.7 % (LBX0190) / +6.9 % (MO_3021) of BAM-scan deposit time** on real cfRNA; cut array 6.0→8.4 MB.
    Affordable, not a cache cliff.
11. **Integer channels are bit-identical at 1/2/4/8 workers; float channels are not** (17/28 and 20/28 cells differ, max relative 3.7e-7).
    That IS the documented ~2.6 % cross-process nondeterminism: parallel FP reduction over a data-dependent chunk→worker split. Integer
    addition is associative; that is the whole fix.
12. **On real cfRNA most confident-gDNA nodes have ZERO counts:** 64 % (vcap), 82 % (LBX0588), **94 %** (LBX0190); count ≤3 in 90 / 97 /
    99.6 %. A density-space estimator floors at 1/E and discards most of the evidence. Genome-wide, **80.5 % of partition nodes carry zero
    fragments** (1,212,140 of 1,505,594) — skipping them was bit-identical and 5.6× faster.
13. **Mature RNA never crosses an exon↔intron seam:** 0 of 1,146 seams, on 7/7 conditions. Exon↔exon seams DO (44,774 mass over 121 seams).
    Transcript-end seams carry 0 of 430 and their true gDNA fraction is exactly 1.00 at p10/p50/p90. **The hard empirical case that a
    contiguous seam and a splice junction are physically different objects.**
14. **But a seam with RNA crossing it need not be a junction.** 33 % of all boundary unspliced mass (44,765 vs 89,128 gDNA) is mature RNA,
    all of it at seams WITHOUT spliced mass. One position can be a splice donor for transcript A and plain contiguous exon for transcript B;
    zero-gDNA libraries show seams with 44–55 k unspliced fragments that are 100 % RNA.
15. **Hybrid capture is ~1000× and only gDNA reads it cleanly.** Enrichment ratio: median 1.0 at intergenic and intron, **1153× at
    single-strand exons, 1212× at ambiguous exons**, continuous over 1→7700; off capture every class is 1.0. gDNA spans 4.6–5.2 decades under
    capture. RNA's own 10⁴ expression range hides the probe pattern; gDNA's uniform baseline does not.
16. **Capture destroys the intron signal 75×:** intronic unspliced density with real nascent RNA and no gDNA goes 0.0979 (off) → 0.0013 (on)
    against a no-nascent baseline of 0. **Nascent vs gDNA is fundamentally unidentifiable under capture.**
17. **Real cfRNA strand parameters** (LBX0190 / LBX0588 / MO_3021 / vcap): sense fraction κ = 0.0023 / 0.0120 / 0.0020 / 0.000060
    — ⚠ **MIRRORED; read these as `1 − κ`, i.e. highly stranded libraries** (§0 C4, measured against ground truth at S5.f). The mirror is
    consistent and cancels inside the solve, so nothing fitted from them is wrong; the reported SCALAR is; RNA strand
    overdispersion 0.0086 / 0.0137 / 0.0074 / 0.0134 (reproduced independently from deep junctions at 0.001–0.016); gDNA strand
    overdispersion 0.2000 / 0.0031 / 0.2000 / 0.0923 — **saturating its 0.2 ceiling on 2 of 4**, because the "pure gDNA" seed weight is
    identically 1.0 on real data, i.e. it re-encodes "the annotation said intergenic" and pools in unannotated transcription. Open defect;
    the raw estimate 0.9313 sits 66–600 σ above the ceiling, so it is bias, not noise, and no robust estimator reaches it.
18. **Real junction-count overdispersion is ≤0.02–0.03 and not fittable** (no technical replicates); naive estimates of 0.08–0.15 decay with
    depth (0.41→0.08 as mean 5→500), the signature of a between-junction fixed effect. **The synthetic suites are Poisson by construction**
    (`sim/wgs_engine.py:473` draws at fixed abundance; measured ω < 5e-5), so nothing dispersion-dependent validates there.
19. **Splice junctions leak:** 1/(fraction of RNA continuing unspliced past a junction) has median 2.3, p75 5.3, p90 9.0 — only ~40 % of
    junctions keep more than half the RNA going. Real junction read counts are small: p10 1, p25 4, median 36, p75 677; 34.7 % below 10.
20. **Fragment length is a real fourth information source:** gDNA vs RNA lengths have total-variation distance 0.886 (gDNA 111±32, RNA
    200±5) and per-fragment Fisher information about the gDNA fraction of 1.19–3.98, against the strand's 1.26 at κ=0.99 and **exactly 0** at
    κ=½. ⚠ The synthetic RNA sd of 5 is unrealistically tight; real size-selected libraries overlap far more.
21. **Where the error is:** over 32 conditions, nodes carried entirely by neighbours are 54.1 % of mass and **91.9 % of all error** (6.2× the
    rate); nodes with own evidence 29.6 % / 8.1 %; structurally-locked pure-gDNA nodes 16.4 % / 0.0 %. And 80.5 % of total error is honest
    under-determination (unstranded, low gDNA); only 1.9 % is confidently wrong.
22. **Calibration cost is depth-INDEPENDENT** — every node in the index is solved regardless of read depth, so a 971 k-fragment targeted BAM
    pays full genome-scale cost. One real run: index load ~7 s, BAM scan ~2 s, calibration ~66 s, per-locus EM ~24 s, peak RSS 8.6 GB.
    **The scan is ~2 % of runtime — accumulator work is essentially free; calibration is the whole budget.**
23. **Memory for the new accumulator is flat and small: ~85 MB** at human scale (node 24 B, contiguous edge 48 B, junction edge 24 B), no
    hash map. The old one is ~89 B per node **replicated per worker** = ~130 MB/worker at 1.5 M nodes. Fixed-point headroom is ~800×: with
    L ∈ [20,2000] each `round(2^32/L)` ≤ 2.1e8, so 1e8 fragments sum to ≤2.1e16 against a uint64 ceiling of 1.8e19.
24. ⛔ **VOID — the old-suite baseline.** It was 0.079005 prior-free / 0.046675 after 3 prior iterations on `ambig_dense_10mb`,
    32/32 reproducible at `3c293038`. **That suite was DELETED on 2026-07-30 along with every other benchmark and every index.**
    The number now refers to nothing: do not quote it, compare against it, or try to reproduce it. Kept only so that a reader who
    finds it quoted in an older document knows it is dead. See `LEDGER.md`'s deletion entry.

---

## §2 EQUATIONS THE CODE DEPENDS ON

- **Crossing count, the core identity.** A fragment `[s, s+w)` spans point `p` iff `s ∈ [p−w, p−1]` — exactly `w` start positions. So
  `E[count at p] = ρ·Σ_w f(w)·w = ρ·mean_FL`, exactly, for any fragment-length distribution and independent of both flanking region sizes.
  **This is why the partitioning problem dissolves.**
- **Reciprocal-length density, at an EDGE only.** `E[Σ1/L] = Σ_c Σ_w ρ_c·w·f_c(w)·(1/w) = Σ_c ρ_c` — the opportunity factor `w` cancels
  identically, across a mixture with different lengths. **At a NODE the opportunity is `max(0, B−w+1)`, which `1/w` does not cancel**; there
  `Σ1/L` is only a better-conditioned second moment. Do not claim node-level model-freeness.
- **The two-equation deconvolution at one object:** `N = ρ_g·E_g[w] + ρ_r·E_r[w]` and `Σ1/L = ρ_g + ρ_r`. Identified iff the two MEAN
  lengths differ; the second row being literally `[1,1]` is what makes the 2×2 well conditioned. `N / Σ(1/L)` is the abundance-weighted mean
  fragment length.
- **General crossing divisor, ONE formula for both edge kinds.** With `R_lo`, `R_hi` the molecule's own remaining sequence either side:
  `E_J = E_f[ min(w−1, R_lo, R_hi, R_lo+R_hi−w+1)_+ ]`. Mean fragment length is its large-reach limit, not a separate case. On RNA N(200,50):
  199.0 at R=550 (median human transcript 1099 bp ⇒ mid-transcript junctions exact), 160.1 at 200, 87.8 at 147, 19.6 at 100, 50.0 at R=50 —
  a **4× error** if mean length is used blindly at a first exon.
- **Exactly why `Σ1/L` fails at a terminus:** `E[Σ1/L] = ρ·E_f[placements(w)/w]`, equal to ρ only where placements ∝ w. At a point 50 bases
  from an end, placements = 50 for EVERY w > 51, independent of w. **Reach:** at position `p` inside exon `e`,
  `reach_lo = exonic_bases_before(e) + (p − exon.start)`, `reach_hi = total_exonic − reach_lo`; maximised over transcripts independently per
  side AND per strand. gDNA unbounded; mature stops at the polyA site. A reach of 0 is meaningful, not a sentinel.
- **Partition rule:** `cuts = unique({exon.start, exon.end over all non-synthetic transcripts} ∪ {0, ref_length})`, node `i` =
  `[cuts[i], cuts[i+1])`, **no merging**. Introns and termini are already exon endpoints, so only the merge step was removed.
- **Contained-node divisor:** `E_f[max(0, L−w+1)] = (L+1)·F(L) − S(L)`, `F` the FL CDF, `S(L)=Σ_{w≤L} w f(w)`; beyond support, `L+1−mean_FL`.
  **The `+1` is the discrete count of start positions, not a fudge** — drop it and the divisor is exactly 0 when a node is one fragment long.
- **The three old divisors** (needed until the old path retires): the old rule gave slice `i` a share `(slice_len/L)/n_crossed_sides`;
  integrating over a uniform start gives per-face seam `E_f[min(w,R)]/2` **exactly for every R** (the ½ is derived, not tuned), a one-sided
  spliced face `E_f[min(w,R)²/(2w)]`, and a pooled two-face seam divisor that is the **SUM** of the two faces, never their average.
- **The strand likelihood — the only intrinsic gDNA/RNA signal.** `p = ½·f_g + κ·(1−f_g)`;
  `var = N·p(1−p) + (N f_g)²·¼·od_g + (N(1−f_g))²·κ(1−κ)·od_r`; `loglik = −½(sense − N p)²/var − ½ log(var)`.
- **What strand can and cannot say.** With RNA tilt `d = f₊ − f₋`, `p = ½ + (κ−½)·d` — **the gDNA fraction cancels identically** (verified to
  5.6e-17). Strand measures the tilt; it reaches gDNA only through the triangle bound `f_g ≤ 1 − |d|` — tight on a single-strand node (that
  IS stranded identifiability), slack on a both-strand node. `I(f_g) = N(2κ−1)²·[f_g(1−f_g)]²/(4p(1−p))`, **exactly zero at κ=½ for any count
  and any overdispersion**, saturating in N at `(½−κ)²/(p(1−p)·od)` — capped by dispersion, not depth.
- **Overdispersion:** symmetric `Beta(a,a)` gives `od = 1/(2a+1)` (a=2→0.200, a=14→0.0345); effective count `n_eff = n/[1+(n−1)·od]` — at
  od 0.2 a 1,523-fragment seed is worth **five coin flips**. Pooled moments:
  `od = Σ_s[(k_s − n_s μ_s)² − n_s μ_s(1−μ_s)] / Σ_s n_s(n_s−1)μ_s(1−μ_s)`. Its exact null information
  `I = (Σn(n−1)pq)² / Σ[2n²p²q² + npq − 6np²q²]` collapses to the pair count `Σn(n−1)/2` **only at mean ½**; at RNA's κ the measured ratio is
  0.05–0.14. Second-moment evidence is counted in PAIRS of fragments inside one object — a singleton carries exactly zero.
- **Why an integer count must be stored:** `Var(log ρ_c) = 1/(f_c·n) ≡ Var(log f_c) + 1/n`, exactly. Mass sums fractional per-fragment
  shares, so `1/mass` is not a counting variance.
- **A faint background rate is measurable only in aggregate:** `ρ_bg = Σg/ΣE`, `Var(log ρ_bg) ≈ 1/Σg`. A region of effective length E
  resolves a rate only above ~1/E (Fisher information = ρ·E), so no per-region estimator finds it and true zero is resolved sharpest.
  **One-sided:** `ρ_bg > 0` proves DNA present; `ρ_bg ≈ 0` does NOT prove absence (capture depletes the off-target floor). Never a
  denominator or a scale.
- **Deconvolving counts against a gDNA background with NO strand data:** `P(g|C) ∝ P_bg(g)·1[0 ≤ g ≤ C]`, `g ~ NegBinom(ρ_bg·E_g, α_eff)`,
  flat one-sided prior on the RNA excess, truncation at the observed total. `α = Σμ²/max(Σ(g−μ)² − Σμ, 0⁺)` (∞ ⇒ Poisson),
  `1/α_eff = 1/α + 1/(Σg + n₀)`. No tuned constant. Measured vs truth: 0.875/0.885, 0.93/0.959, ~0/0.000. **The one part of the solver that
  demonstrably works.**
- **Self-limiting, knob-free damping.** Treat a claim and the destination's own estimate as two studies of one quantity:
  `b̂² = max(0, G² − v_msg − v_own)` with `G` the log gap ⇒ `p_eff = 1/max(v_msg, G² − v_own)`. Exact safety property: **a claim can outweigh
  a node's own belief only if it agrees within `√2·σ_own`**; where a node has no evidence `v_own = ∞` and the term switches itself off.
- **Density is the frame-invariant currency; a fraction is not.** `ρ_c = C_c/E_c` agrees across the contained / crossing / spliced frames to
  ~0.36 % and does not degrade across 1×/30×/300× capture. The log-odds shifts between frames by exactly
  `log(E_g^dst/E_g^src) − log(E_r^dst/E_r^src)`, and capture cancels identically — **a ratio transports across a capture cliff, an absolute
  density does not.** ⚠ Under `Σ1/L` each fragment's contribution is FL-independent, so this shift should vanish — verify it.
- **Conservation with unequal effective lengths is `Σ_c ρ_c·E_c = M`, not `Σ_c ρ_c = M/E`** — and its sensitivity is bounded: a purely
  compositional error moves `M/Σρ_c E_c` by only ×1.04 on a contained node and ×1.50 at a crossing, so a large violation is accumulated
  drift, never one hop.
- **"No prior" does not exist on a grid — omitting a term lets the grid supply Haldane** (`p(ρ)∝1/ρ`, improper, an amplifier toward the
  vertices). `p(ρ_c) ∝ ρ_c^(c−1)` gives Beta(c,c); `c=½` (Jeffreys) is the only grid-width-stable choice: posterior median spread over
  half-widths 4–20 is **0.045 at c=½ vs 0.525 at Haldane**.
- **Sums are well conditioned, differences are not.** Subtracting across a junction gives `Var(log ρ) = u²σ_T² + (u−1)²σ_μ²`,
  `u = 1/(continuing share)` — at the real median u=2.3 already at the edge of validity, at p75 u=5.3 hopeless at any depth. Prefer shares.

---

## §3 TRAPS — do not repeat these

1. **A validator that calls the builder's own helper validates nothing.** Deleting the negative-strand transcript-end swap left the ENTIRE
   1,289-test suite green and the graph validator accepting; only re-deriving the events by a *different algorithm* caught it. Worse than no
   check, because it reads as one. Corollaries: a validator that only compares non-empty classes never sees a spurious flag (emit all classes
   always); prove every validator FIRES by deliberately corrupting the input.
2. **A `min()` clip hid an exact factor of 2 for months.** The pooled seam density read 1.994/2.002/1.981× truth because the code SUMMED two
   faces' mass but divided by their AVERAGE, when each face divisor was already halved. It survived 29 tests: all four fixtures stored an
   un-halved length AND deposited half the mass (cancelling exactly), and `min(2S,S)=S` under a uniform field so the "contraction is 1 under
   a uniform field" bedrock invariant — written to catch exactly this — passed with the bug present. **Fixed 2026-07-29, do not fix twice.**
   Lessons: a clip can hide a scale error from a bedrock test; repair fixtures, never relax assertions; **an exact algebraic 2 is never a
   modelling approximation.**
3. **`~is_synthetic & ~is_nrna` as the "real transcript" filter deleted the start/end positions of 26,475 real single-exon transcripts
   (52,104 positions)** — exactly the visibility the redesign buys. On a non-synthetic row `is_nrna` means "single-exon, so mature ≡
   nascent", NOT "manufactured span". **ONE filter: `~is_synthetic`.** Nothing noticed because the flags had no consumer yet.
4. **Fractional mass IS the partitioning problem.** A fragment spanning 4 nodes writes six fractional numbers whose values depend on region
   sizes, purely because a mass is conserved; the same fragment's three crossing counts depend on nothing. **Corollary: multimapper /
   ambiguous-path assignment must stay INTEGRAL**, or the non-integer observable returns and the count stops being a count.
5. **Mass conservation does not catch mis-attribution.** In the shipped accumulator one fragment can credit the same boundary side TWICE (a
   spliced-out intron lying wholly inside one node makes that node appear twice in the slice list) and can credit a boundary it never
   crossed — with total mass still exactly 1.0. Verified on both the Python reference and the native build.
6. **A splice deposit is attributed to the NODE's edge, not the junction's coordinate.** Invisible for annotated introns (their ends are
   cuts); for unannotated junctions the mass lands kilobases away.
7. **A splice junction CANNOT be detected from a gap between deposited slices.** A contiguous spliced read whose exon body straddles an
   internal cut with no intron in the gap has no gap at all. The junction's identity is the cut-intron coordinates `(start, end, strand)` —
   which the BAM scanner already has and throws away before depositing. **Pass them through.**
8. **Spliced mass is deposited at covered (aligned) length scale but divided by a fragment-length divisor,** so mature density reads ~2.3×
   low and **63 % of mature RNA is misread as gDNA** (projected/true = 0.370; 409,148 fragments become phantom gDNA), stable 0.35–0.37 across
   capture and strandedness. Same family: binning the gDNA FL pool at *covered* length collapses it to a spike at 2×read-length, so every
   long gDNA fragment scores as RNA.
9. **The path / cell store was killed by MEMORY, not taste.** Possible unspliced paths ≈ 1 M nodes × 3–6 reachable ends = 3–6 M at ~100 B
   each = 0.3–0.6 GB, plus spliced paths, and a deep library saturates a large share. Owner's shorter verdict: expensive, and no consumer
   needs it. Do not re-propose path or cell enumeration without that census.
10. **Alternative splicing makes a node↔junction graph CYCLIC, and cycles are the common case** — a cassette exon is a 4-cycle, and the human
    graph has `edges − nodes + refs ≈ 404,000` independent cycles, one per junction. Two-sweep forward-backward is exact only on a tree.
    **Never break a cycle by dropping a junction edge** — that re-isolates the exon the edge exists for.
11. **A message computed from the destination's own belief or own total density confirms the destination and carries zero information.** Hit
    twice independently: one delivered an RNA fraction of exactly `1/(1+f_g_own)` (2.1e-16), reserving 33.6 % of the budget for imaginary
    gDNA so a **zero-gDNA library read back 29.3 % gDNA**; and at a gene end no RNA crosses, so the seam is pure gDNA, the rescale ratio
    becomes `1/f_g(dst)` and the delivered gDNA density is exactly the destination's own total (measured 7× too big in the median, up to
    190×, and **96 % of the error on that path**). **Rule: a message may use the destination's CONSTANTS (geometry, lengths), never its
    BELIEFS.** Any "fix" that divides the destination's belief back out rebuilds the bug.
12. **You cannot fix a biased mode with a variance.** Established three times independently. Under capture the counting estimate is
    systematically ~2× low but PRECISE — both flanking seams sit at the same enriched exon edge and agree on the same biased-low density, so
    the bias is trusted. **A disagreement-based variance model structurally cannot fix a bias.**
13. **Never hand a solver two Gaussians built from one latent.** A message on `log f_g` and one on `log(1−f_g)` are rank-1 with correlation
    exactly −1, so adding their Fisher information is **exactly 2× over-confident**, rising to ~7× with deep spliced content.
14. **Never fit a variance on the current, not-yet-solved belief.** Adjacent WRONG nodes agree, so the variance collapses, the messages turn
    confident, and the error propagates. Ablation: prior-refit only 0.118 ✓ / "honest" measured variance only 0.189 ✗ (worst corner 0.52 →
    0.95). *Honesty measured against a wrong belief is not honesty.* Same family: any component trained on the solver's own output is
    self-confirming — refit iterations 1→5 gave 0.0840 → 0.1056, monotonically worse.
15. ⛔ **DELETED 2026-07-30, and this is WHY it needed replacing. The 32-condition `ambig_dense_10mb` suite could not judge what it was used to judge.** Its fine node set is row-for-row IDENTICAL to its
    merged region set (1,698 == 1,698), so it cannot see a partition change; `frag_std = 0` (every fragment exactly 200 bp), so it cannot
    test anything fragment-length dependent including `Σ1/L`; it is Poisson by construction; it overstates contained efficiency 4.6 ×,
    understates multi-node crossing 3.8×, and over-represents the terminus+splice-site seam 12× (11.8 % vs 1.0 %). **Before running a
    benchmark, prove it can resolve the axis you are changing.** ⭐ This list is the requirements document for the replacement.
16. **The toy also ranks performance hotspots BACKWARDS.** At 3.4 k vs 1.5 M nodes: message relay 34 % → 81 %, per-node solves 9 % → 29 %,
    the prior's EM **28 % → 0.7 %**. A whole careful analysis was wasted on the toy's #1 hotspot. Profile on cached real cfRNA payloads.
17. **A bit-identity gate lied twice, in opposite directions.** An arm with ZERO rows scored "32/32 IDENTICAL" because the comparison looped
    over the new arm's rows; and a stored baseline went stale so unmodified HEAD no longer reproduced it. **Re-record the baseline from a
    `git stash` of HEAD in the same session; if HEAD-vs-baseline is not 100 %, the baseline is what is broken.**
18. **Every "is the declared precision earned?" number written before 2026-07-28 compared a LOG-space variance against a LINEAR squared
    error.** Corrected, the suite total goes 0.046 → 1.007 and the per-class ranking INVERTS. Any pre-2026-07-28 claim of over-confidence is
    void.
19. **A zero-gDNA "phantom" guard is ONE-SIDED** — in a library with no gDNA, any change that lowers the gDNA fraction scores better. It
    reversed a published verdict once (reported −13.1 % win; the full 32-condition battery made it 0.1396 → 0.1416 worse, stranded arm
    +51 %). Related: hard-label metrics (`net_flow`, MAP fragment labels) are byte-identical under changes that move the soft pools by tens
    of thousands of fragments. Score per node and per condition on soft quantities, never pooled — one run had a −47.6 k strand error and a
    +75.4 k prior error reading as a "nearly perfect" −1.9 k pool.
20. **Single-reference synthetic indexes hide reference-id-space mismatches.** The C++ resolver assigned ref-ids by first-seen interval order
    rather than the index's own order, so **476,719 of 476,732 real fragments were silently dropped inside `deposit()`** while every golden
    test passed.
21. **A fragment crossing K splice junctions must credit exactly ONE** (the leftmost annotated). Crediting all K shifts the library sense
    fraction 21–34 % and creates between-side correlation that reads as overdispersion (a simulator with ZERO overdispersion fit 0.092). A
    `1/K` split is provably biased by 4–12 σ, because `Var(Σ wᵢXᵢ) = pq·Σw²` while the estimator subtracts `pq·Σw`.
22. **"Effective lengths cancel, so a node's marginal is just its length" is FALSE near any transcript end** — a mature fragment must fit in
    the remaining transcript, gDNA need not.
23. **Density below one fragment length is not resolvable by ANY accumulator design.** A 1 bp node has no independently measurable density
    and never will (*composition* still is, since it depends on what the fragments are, not where). Corollaries: an object with zero
    opportunity for a component must emit nothing at zero precision, not a floored division; and "no data" must be inert, never "100 % gDNA"
    — that default was actively seeding false gDNA into neighbouring exons.
24. **`(src, kind, dst)` is not a total order for junction edges** — two strand-coincident junctions differ only in strand, so ordering
    becomes input-order-dependent and the duplicate check reads them as duplicates. GENCODE has zero of them (splice motifs are
    non-palindromic, so they are biologically impossible), so only a synthetic stress test can find it. Sort on `(src, kind, dst, strand)`.
25. **Cache keys that do not cover the artifact they cache.** The partition hash covers only the node file; a flag fix rewrote every edge
    file while leaving every node file byte-identical, so a stale cache would verify CLEAN and feed every downstream comparison the pre-fix
    flags. `/tmp/rigel_selfsolve` is also shared and un-namespaced — already cost two full index rebuilds. Never store a derived hash beside
    the data it describes; compute it on demand (39 ms at human scale).
26. **Worktrees silently run the wrong code** — an editable install's meta-path finder beats `PYTHONPATH`, so an A/B inside a git worktree
    executes the main repo's source.
27. **The prose next to the code said "the AVERAGE", the code followed the prose,** while a sibling module's docstring had the correct
    formula the whole time — two docstrings disagreed about one quantity for months and nobody diffed them.
28. ⛔ **A SUPPORT CEILING THAT MATCHES THE CLAMP IS NOT A MATCH — it is the clamp.** C1 recorded the fragment-length anchor's ceiling moving
    from 713 bp to 1000 as "the mismatch closed ⭐ entirely"; **713 was the library's true maximum** (`truth_fragment_lengths.tsv`) and 1000 is
    `max_frag_length`. The narrower estimate was the right one, and the "fix" was an uncut intron. ⚠ Whenever a distribution's support agrees
    with a configured limit, that is evidence the limit is binding, not evidence of agreement. C2.6, 2026-08-01.
29. ⭐ **A PURITY FILTER ON A LENGTH POOL IS A LENGTH FILTER.** C2.6's D1 bars a fragment from the pure-RNA pool when part of its `L` was
    inferred rather than sequenced — but the fragments that qualify are exactly the ones whose mates sit far apart, so the pool the
    fragment-length model is FITTED FROM became length-selected **short**: **−9.58 % mean / −22.46 % sd** vs truth, where keeping them reads
    **+0.67 % / +2.40 %**. ⚠ Before excluding a population from a pool, ask what the exclusion criterion correlates with — here it correlates
    with the axis being measured, so purity and accuracy point in opposite directions.
30. ⚠ **`git checkout -- <file>` does not undo a perturbation when the work is uncommitted — it deletes the work.** A perturbation harness must
    restore from a copy of the WORKING TREE. Cost one full re-implementation of a landed change, 2026-08-01.

---

## §4 — DELETED 2026-07-30

It was a **pre-settlement proposal list**, distilled before the design was agreed. Most of it became
`ACCUMULATOR_DESIGN.md`; six of its bullets had already been measured, rejected and marked
`⛔ SUPERSEDED` in place. A proposal list that contradicts the shipped deposit rule is a hazard, not
history — `ACCUMULATOR_DESIGN.md` and `tests/native/_accumulator_reference.py` are the authority, and
git retains the text.

⭐ Two items from it survived into live documents rather than being lost: the benchmark-suite
requirements are now `BENCHMARK_SUITE.md` plus `TODO.md`, and the concrete accumulator test matrix is
`tests/native/test_accumulator_spec.py`.
