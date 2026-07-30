# RIGEL — CARRY FORWARD

All that survives from 278 docs / 74,823 lines + 142 memory files. Read §0 (open decisions) then §3 (traps — the most expensive items).

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
LBX0190 **0.0562**, MO_3021 0.1264, LBX0588 0.7718. No ground truth exists. Also unexplained: sense fraction is 0.002–0.012 on all of them
(nearly fully antisense) — possibly a read-orientation convention bug. (`SESSION_2026_07_30_HANDOFF.md:333-341` vs memory `cfrna_sample_characteristics`)

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

1. **The human splice graph (verified on disk, index format v8):** 1,043,881 nodes (median 151 bp, 15,687 of length 1), 1,043,595 contiguous
   edges, 404,168 junction edges = 1,447,763 edges. Replaces 752,654 merged regions, +38.7 %. (`splice_graph.py:16-19`)
2. **53.4 % of real human transcript ends (232,451 of 435,291) fall strictly inside an old merged region** and were invisible to it. This is
   the entire reason for the redesign.
3. **Seam census over 1,043,595 contiguous edges:** terminus-only 40.70 %, splice-site-only 58.31 %, BOTH 0.99 % (10,337), neither 0.00 %.
4. **The crossing-count estimator is unbiased and region-size-independent:** `count / mean_FL ÷ truth` = 0.997 / 1.001 / 0.987 at region
   lengths 2000 / 500 / 200 over 395 seams. This single measurement justifies the whole deposit rule. (`scratchpad/acc_seam_check.py`)
5. **`(count, Σ1/L)` beats `(count, ΣL)` everywhere:** sd 0.017 vs 0.028, condition number 283 vs 1137 (gDNA FL 50, RNA 200, f_g 0.15);
   1.6× more precise, 4× better conditioned, never worse. With the fitted gDNA length wrong by 10 %, total density reads 1.0012× truth vs
   0.9652×. (`scratchpad/acc_proto_g.py`)
6. **Ignoring the fragment-fit taper near a transcript end over-calls gDNA by 11.0 %** — bp-weighted mean of (crossing effective length ÷
   mean FL) = 0.8904 over 3,000 genes / 10,957,350 exonic positions, RNA N(200,50). **Contiguous seams are WORSE than junctions (0.750 vs
   0.886)** — a junction has exons both sides by construction; a contiguous seam is often inside a terminal exon. The gDNA fraction is off
   by **+0.36** in the last region before a polyA site.
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
17. **Real cfRNA strand parameters** (LBX0190 / LBX0588 / MO_3021 / vcap): sense fraction κ = 0.0023 / 0.0120 / 0.0020 / 0.000060; RNA strand
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
24. **Old-suite baseline** (only if `ambig_dense_10mb` survives): mass-weighted error 0.079005 prior-free / 0.046675 after 3 prior
    iterations, 32/32 exactly reproducible at commit `3c293038`.

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
15. **The 32-condition `ambig_dense_10mb` suite cannot judge what it is used to judge.** Its fine node set is row-for-row IDENTICAL to its
    merged region set (1,698 == 1,698), so it cannot see a partition change; `frag_std = 0` (every fragment exactly 200 bp), so it cannot
    test anything fragment-length dependent including `Σ1/L`; it is Poisson by construction; it overstates contained efficiency 4.6 ×,
    understates multi-node crossing 3.8×, and over-represents the terminus+splice-site seam 12× (11.8 % vs 1.0 %). **Before running a
    benchmark, prove it can resolve the axis you are changing.**
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

---

## §4 IDEAS FOR THE NEW DESIGN

> ⛔ **READ THIS FIRST. §4 is a PRE-SETTLEMENT proposal list**, distilled from the deleted documentation
> before the design was agreed. Most of it survived into `ACCUMULATOR_DESIGN.md` — but **several items were
> subsequently measured and rejected**, and they are marked `⛔ SUPERSEDED` in place with what replaced
> them. **`ACCUMULATOR_DESIGN.md` and `tests/native/_accumulator_reference.py` are the authority; §4 is
> history.** The markers were added 2026-07-29, when six of them were found to contradict the deposit rule
> that S3 implements — an unmarked §4 would have handed a C++ author four wrong instructions.

**Object model — the index already ships it** (format v8, on disk, validated, 56 graph tests passing)

- Node = half-open interval tiling a reference, ids contiguous in genomic order. Edge = directed, always `src < dst`, so **genomic order is
  a topological order — no graph traversal anywhere**. `kind ∈ {0 CONTIGUOUS, 1 JUNCTION}`. Row order is the on-disk contract, so a node's
  out-edges are contiguous and a CSR is one `searchsorted` (none is built at load today — three O(n) passes to write).
- **8 structural flag bits per contiguous edge** (uint16, deliberately NOT mutually exclusive): TSS / TES / DONOR / ACCEPTOR × {+,−}. This
  is exactly what the old 4-bit signature was blind to, and it is 96 % of the largest known calibration error. **Carry the raw bits to the
  consumer; do not pre-derive predicates in the plumbing** — the one predicate specified up front turned out to be nearly the complement of
  what it replaced.
- **4 reach columns, named `lo`/`hi`, stored PER STRAND.** Never `donor`/`acceptor`: `src` is genomically left whatever the strand, so on a
  negative-strand junction the biological donor is on the right and the names mislabel half of 404,168 junctions. Per strand because at a
  mixed-strand seam the two genuinely differ (a strand-agnostic max over-stated one by 10×).
- **1 bp nodes are legal (15,687 of them); nothing may assume length > 1.** A reference with no transcripts is one node and zero edges.

**Deposit**

- **One rule: a fragment deposits on EVERY edge it crosses, with no partitioning**, because each edge answers its own question and there is no
  total to conserve across edges. ⛔ **SUPERSEDED in three particulars** (design §4.1, §4.3, §3.1) — as originally written this bullet said
  `+1/L` at an edge, "crossing none ⇒ contained", and "`L` is exonic bases only":
  - the edge weight is **`1/(L−1)`**, not `1/L`. A 0-bp line means bases on *both* sides, so there are `L−1` admissible placements; `1/L`
    reads **0.54 % low**, by exactly `1 − E[1/L]`, verified to four digits on identical draws.
  - **"crossing none ⇒ contained" is a BUG**, caught by the S2 spec matrix. An unannotated intron can swallow every line between two blocks,
    leaving a fragment that crosses nothing yet **straddles two nodes**. Contained means *the whole path lies inside one node*.
  - **`L` is NOT "exonic bases only"** — it is the genomic span with introns cut out, so the **mate gap counts**. Dropping it collapses the
    length distribution to a spike at twice the read length, and the design calls this out explicitly (§3.1).
- **Never divide by L at deposit time**, and never move normalisation into the deposit as a per-component correction — that needs to know
  the component, which is the estimand. `Σ1/L` works *because* it is component-blind. Normalisation lives in the effective lengths.
- **Everything integer or fixed point.** This REMOVES the known ~2.6 % cross-process nondeterminism. ⛔ The quantum is
  `round(2^32/placements)` with **`placements = L` at a node and `L−1` at an edge** — two different quantities, deliberately never given one
  column name (design §10.1). This bullet said `round(2^32/L)` for both.
- **Keep depositing the molecule's genomic span** (mate gap filled, spliced introns cut), not the sequenced read blocks — it removes the
  paired-end over-count at source. Keep excluding blacklisted-junction fragments entirely (their true span is unrecoverable).
- ⛔ **SUPERSEDED — the "three exact integer conservation identities" are TAUTOLOGIES** and were rejected (design §10.2). Each right-hand side
  can only be evaluated by re-running the deposit, so they cannot fail: the S2 review's deliberately-broken replay satisfied all three exactly
  while 91 % of the crossings were junk. **Replaced by one `uint32` per node counting fragments whose START lies in it**, giving
  `Σ start_count == n_accepted` — checkable against a number the scanner knows independently. ⚠ The denominator warning still stands:
  multimappers, chimeras, blacklisted and strand-undefined fragments never deposit, and overhanging fragments are clipped.
- **Unannotated junctions: never mutate the graph**, and report the rate as QC — a rising rate means a stale annotation. ⛔ **SUPERSEDED on the
  channel: they deposit UNSPLICED, not spliced** (owner ruling; design §4.1). Unannotated junctions are disproportionately artifactual, so they
  must compete with gDNA rather than be certified as RNA — the safe direction. They also deposit **nothing across the gap**. ⭐ And the rate is
  now measured: **2.3 % / 1.7 %** of `N` operations, but **18.9 % / 10.8 % of *distinct* junctions**, which is the honest figure.
- **Ambiguous-path fragments go to a small side buffer** `(candidate ids, channel, L)`, resolved in a second pass by maximum likelihood,
  INTEGRALLY, after the RNA fragment-length model is fitted. Size unmeasured.

**Solver shape**

- **Treat junction edges as FACTORS on their two endpoint nodes, not as message channels.** The chain stays linear (node ↔ contiguous-edge ↔
  node) so the forward-backward sweep stays exact. ~404,000 loops would make it approximate, remove any convergence guarantee, and
  double-count evidence around cycles — precisely the over-confidence failure this project has already hit twice. You still get gene-end
  typing, path counts and structural effective lengths.
- **Terminus seams are new, clean, structural pure-gDNA anchors.** Where no transcript of strand `s` spans a seam, the RNA opportunity is
  ZERO *by annotation* — declared, not guessed, and no ratio against a neighbour's own belief is ever formed. Today the only clean gDNA
  anchors are intergenic regions.
- **Carry per-component precision (the 2×2 from the moment solve), not one scalar** — a 100 bp node is a good gDNA measurement and says
  nothing about RNA; tell the solver that rather than let it believe a collapsed column.
- **Source, sink and relay are one precision-weighted fusion,** `running = fuse(own, incoming)`. A high-count intergenic region automatically
  sources, a low-count seam relays, an exon blends. Measured over 32 conditions: error 0.1396 → 0.1361, correlation 0.671 → 0.688, stranded
  arm neutral. **And always emit** — a zero-precision message is ignored anyway, but a boolean emission gate once silenced 52 % of real
  splice-junction measurements.
- **`claim / observed`** — fragments the incoming claims assert vs fragments the node saw — **is a free, oracle-free, prior-free diagnostic
  of message error** (1.00–1.01× ⇒ error 0.0007–0.0056; 1.26–1.33× ⇒ 0.165–0.218). The solver computes it and erases it. Use it, and
  apportion the correction by how well each component is known instead of rescaling both blindly.
- **The prior re-solve must fully reset** — nothing survives except the prior itself, so an over-confident node cannot refuse to budge.

**Fragment-length models**

- **Structurally pure pools remove a circular dependency:** junction-edge fragments are pure RNA by construction, intergenic contained fragments
  are pure gDNA. Ordering: accumulate (no FL model needed) → fit both FL distributions from the pure pools → effective lengths → calibrate.
  Nothing is estimated from the fragments it will later explain. ⛔ **SUPERSEDED in three particulars** (design §8, §8.1):
  - **"pure mature RNA" → pure RNA.** There is no mature/nascent distinction anywhere in the accumulator (owner ruling: "RNA is RNA"), and a
    partially spliced nascent molecule carries junctions too.
  - **the mean-length estimator is `count/Σ(1/(L−1)) + 1`, not `count/Σ(1/L)`.** The `+1` is exact at an edge and verified to +0.0001 % on all
    four libraries; omitting it is low by exactly `1/E[L]`. And at a **contained** node the same ratio converges to the **harmonic** mean —
    which is a problem, because intergenic-contained is the *only* gDNA pool, i.e. exactly the frame where the naive recipe fails.
  - **five pools, not two**, and the splash pools are named rather than optional: DNA intergenic / intronic / intron-exon / intergenic-exon and
    RNA spliced. ⭐ Measured means 88.0 / 88.5 / 138.8 / 211.6 / 220.7 — the splash pools land *between* the two extremes, as expected.
- **Bin each fragment's length ONCE, as an integer, into a named pool.** Today a crossing fragment smears its length fractionally across
  every region type it touches, accumulated as non-associative `double +=`.

**Validation and benchmarking**

- **Write the toy Python accumulator FIRST and make it the executable spec**, with a byte-for-byte native parity harness. ⭐ **Done in S2** —
  `tests/native/_accumulator_reference.py` plus a 46-case spec matrix, which caught four blocking bugs no later gate could have.
- **The oracle IS the production accumulator, partitioned by true fragment origin.** Split the simulated BAM by read-name origin, run the
  SAME scanner and accumulator on each partition, assert the partitions sum to the full payload. Four hard gates: integer channels sum
  exactly; mass within float rounding; gDNA is never spliced; every read lands in exactly one partition. Nothing is reimplemented, so
  nothing can be subtly wrong. (An earlier oracle on a *different* deposit basis manufactured a spurious 157 k "fix" and confounded months.)
- **Keep the concrete test matrix:** contained → node only · 4-node crossing → exactly THREE edges · fragment spanning a 1 bp node → both
  flanking edges **and that node's SPANNING counter, with contained untouched** (⛔ this bullet said "node untouched", from before `spanning`
  was stored) · spliced → junction edge only, NO deposit on the contiguous edges it splices over · zero-length intron · N workers →
  bit-identical · opposite-strand junctions at identical coordinates are distinct edges · unannotated junctions · clipping at a reference end ·
  effective lengths vs brute-force enumeration · uniform-density recovery within [0.98, 1.02] · `Σ1/(L−1)` within [0.98, 1.02] **with no FL
  model supplied**. ⛔ "exact count conservation" is dropped — see the tautology note under **Deposit**.
- **An executable unbiasedness test already exists** and is exactly the property the new edge counter needs:
  `tests/calibration/test_accumulator_span_unbiased.py` asserts crossing density ÷ contained density ≈ 1.0. Keep it, swap the estimator.
- **Validate by re-deriving with a DIFFERENT algorithm**, sharing only the input definition (node signature by midpoint containment vs
  cumulative difference arrays; terminus events by per-transcript min/max vs cumulative exonic offsets — both caught real bugs).
- **The new benchmark suite (owner's own plan):** real human genes as backbone, so calibration gets real RNA and DNA fragment-length
  distributions and a real strand model; **cache the whole BAM-scan intermediate** (FL models, strand model, accumulated counts); then
  piggyback a fake chromosome of toy stress cases onto the cached objects. Caching took a 24-condition run from ~13 min to **~9 s**,
  bit-identically (`scripts/debug/calib_cache.py` does this for the old path — rebuild it day one).
- **Design the blind spots deliberately.** From five separate failures the suite MUST have: fragment-length variance; alternative TSS/TES
  inside exons; non-Poisson counts; ample single-stranded nodes in any both-strand stress test (the population prior trains on
  single-stranded solved nodes); and a **low-gDNA × strong-capture corner (1–10 % gDNA)** — real libraries live there and 0/100/300 %
  conditions cannot see it. Include a density STEP, not just a uniform background: over a run of flat nodes a relayed message decays
  geometrically per hop, so a uniform scenario cannot distinguish "the relay works" from "the global prior reached it".
- **Land ONE effect per arm and append a ledger row as each lands.** Four effects move the same numbers and are only attributable if each
  delta is recorded against a baseline re-recorded from the same tree in the same session.
- **Build the coarsening map (new nodes → old regions) BEFORE the first partition arm.** Once the partition moves, per-node error statistics
  are not comparable across node sets — and scoring granularity ALONE moves the error metric 26 %.
- ⛔ **SUPERSEDED — "shadow mode: emit BOTH payloads from one scan" was REJECTED.** Owner ruling: **clean cut** (plan §5). Two payloads is
  duplicate state kept for comparison, which the standing rules forbid, and it would put both halves of the rework in flight at once — which is
  exactly what destroys the oracle. What replaces the attribution it would have bought is the **Python reference**: a complete correct
  implementation before any C++, so "is the C++ right?" and "is calibration rewired right?" are answered independently. Same for "ship new index
  tables additively so old indexes stay loadable" — no backwards compatibility.
- **Take the terminus win off the INDEX ALONE first** — the terminus-aware solver fixes need no accumulator change, so that value lands
  before any of the risky work.
- **Prefer structural (annotation-derived) predicates to observational (coverage-derived) ones:** a structural spliced-face count is constant
  at 503 across all 32 conditions; the observational one collapses 526 → 392 → 263 as capture strengthens — blind where evidence is scarcest.
- **The debug loop that works** (owner directive): full suite → worst scenario ranked by error MASS (`Σ mass·|error|`, not mean error) →
  worst objects inside it → root cause with measurements → fix → re-measure → repeat. Never optimise a metric without a mechanism.
- **Verify refactors by full per-node state, not aggregate scores** — every solved fraction, variance and precision across 8 conditions takes
  ~8 s and is strictly stronger than a 22-minute aggregate benchmark. If a "no-op" moves a number, stop and think, don't relax the gate.

**Standing owner directives** (all verified still binding)

- **No magic numbers.** Stop and discuss before adding any constant, heuristic or tunable; the previous calibration was burned down at ~91
  constants / ~25 qualitative cliffs; target ≤8 numeric literals. Every divisor must be *derived from the deposit rule* and unit-tested
  against brute-force enumeration — both the ½ at a seam and the `+1` on a contained node fell out that way.
- **No Greek letters in code** (fine in maths docs). Keep the module and constant count small. Develop mechanisms on controlled toys,
  validate on real data. Simplification is a gate, not a nicety.
- **All build / test / lint runs inside the activated `rigel` conda environment**; re-run `pip install --no-build-isolation -e .` after any
  C++ change.
