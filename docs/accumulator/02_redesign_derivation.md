# Accumulator redesign — specification

**v2, 2026-07-28.** Supersedes v1 (2026-07-27) — see §0 for what was wrong in it.
Concept: `dag_design.md`. Prerequisite that already shipped: `accumulator_fragment_span_redesign.md`
(deposit works on the molecule's contiguous spans, not read blocks).
Prototype: `scratchpad/acc_proto.py`. Measurements: `scratchpad/acc_measure{,2,3,4}.py`,
`scratchpad/acc_terminus_census.py`.

---

## 0. RETRACTIONS FROM v1

| v1 said | correct |
|---|---|
| the primary object is the node marginal, `E_c(A) = \|A\|` for every component, so effective lengths **cancel** | **False near any transcript terminus.** A mature fragment must fit inside the remaining transcript; gDNA has no such constraint. Uncorrected this over-calls `f_g` by **+0.36** in the last region before a TES (§3.2), and costs **10.7 %** averaged over all human exonic bp (§10). Fixed by the *admissible-start* effective length (§2.2). |
| a spliced fragment deposits **+1 per junction traversed** | **Wrong — it is a length-dependent count bias**, exactly as the owner said: a 5-junction fragment would score 5 and a 1-junction fragment 1. One fragment = **one count on one path** (§2.1). Junction flux is a *derived marginal*, never an accumulator. |
| three accumulators (node / cell / junction) | **One table**, keyed by path. Everything else is a projection (§4, §7). |
| paths collapse to `(first, last)` | **Paths are stored in full.** For unspliced fragments `(A,Z)` *is* the path; for spliced ones the junction list is kept (owner decision D2). |
| defer deposit through `FragmentBuffer` | Rejected (owner decision D4). `FragmentBuffer` holds equivalence classes, not fragment structure; loading it up would be unwieldy. Instead: **deposit directly during the scan, and buffer only the genuinely ambiguous fragments** in a small separate side buffer (§8). |
| lead with the anchor | The anchor was presented as if it were the storage. **It is not.** The storage is the path; the anchor is the path's first element and exists only to define which density the path's count is evidence *about* (§2.1). |

---

## 1. ESTIMAND

Per node, per component `c ∈ {G = gDNA, V = nascent RNA, M = mature RNA}`:

```
    ρ_c(A)  =  molecules of component c per bp of node A          [molecules · bp⁻¹]
```

Coverage is `ρ_c · E_{f_c}[w]` — a known constant per component times `ρ`, so the two framings carry
identical information and "molecules per bp" is chosen because the downstream consumer
(`priors.assemble_priors` → the per-locus Dirichlet) counts **fragments**, not bases.

Composition is `f_g = ρ_G / (ρ_G + ρ_V + ρ_M)` over the node's unspliced mass, unchanged.

---

## 2. OBSERVATION MODEL

### 2.1 The atom is the PATH

A fragment occupies a set of contiguous genomic spans (already produced by
`fragment_genomic_spans`). Its **path** `p` is the ordered sequence of regions it occupies together with
the junctions it uses:

```
    unspliced:  p = (A, A+1, …, Z)          — determined entirely by (A, Z)
    spliced:    p = (A, …, J₁, …, J₂, …, Z) — the junction list is part of the key
```

**Deposit rule: one fragment increments `count[p][channel]` by exactly 1.** No fraction, no fragment
entering two objects, no float anywhere.

`a(p)` (the first region) is the path's **anchor**. Its only role is this: a molecule is generated once,
by the process at its start position, so the count of paths anchored in `A` is evidence about `ρ(A)` and
about nothing else. Crediting every region a fragment covers would count one molecule several times —
which is what the fractional split exists to paper over.

> **The anchor is not a summary of the fragment.** The path is the complete description — it *is* the
> coverage pattern, and a 60 bp gDNA fragment and a 200 bp RNA fragment starting at the same base land on
> **different paths**. Nothing about length is discarded; §3 shows it is the main new signal.

### 2.2 Admissible-start effective length ⭐ (the terminus correction)

Each component has a **support** `S_c` and a FL pmf `f_c` with CDF `F_c`:

| c | support `S_c` | reach from position `s` |
|---|---|---|
| G | the reference, contiguous | `reach_G(s) = ref_end − s` |
| V | the pre-mRNA span, contiguous | `reach_V(s) = span_end − s` |
| M | the transcript's exonic path | `reach_M(s) = remaining transcript length from s` |

A molecule `(s, w)` exists only if it fits: `w ≤ reach_c(s)`. Therefore

```
  E1  path      E_c(p) = Σ_w f_c(w) · #{ s ∈ S_c : molecule (s,w) occupies exactly path p }
  E2  node      E_c(A) = Σ_{p : a(p)=A} E_c(p) = Σ_{s ∈ A ∩ S_c} F_c( reach_c(s) )
```

E2 follows from E1 by exchanging the sums; it is the **admissible-start count** and it is a prefix-sum of
a CDF lookup — `O(|A|)` to build once per (index, FL), or closed form on a piecewise-linear reach.

* `E_c(A) = |A ∩ S_c|` **exactly** wherever every `s ∈ A` has `reach_c(s)` beyond the FL support.
* It **tapers to 0** across the last FL bases of the support. That taper is the whole of the terminus
  effect and it is why v1's cancellation claim was false.
* `E_G` never tapers inside a reference; `E_M` tapers at every TES. **This is exactly the structural bit
  P1g needs**, so the terminus correction and the reframe fix are computed from the same annotation.

⚠ `reach_M(s)` is per transcript. Where isoforms of different length share `s`, the honest quantity is
`Σ_t θ_t F(reach_t(s))` with unknown `θ_t`. Calibration does not resolve isoforms — **use the annotation's
maximal remaining path** and record the sensitivity as an open item (§11).

### 2.3 Conservation

`Σ_{p : a(p)=A} E_c(p) = E_c(A)` by E2 — an identity, integers on the observable side, no ½ factors, at
every region length down to 1 bp. This replaces v1's fractional mass conservation and v1's
`Σ_Z E(A,Z) = |A|` (which was the special case `reach = ∞`).

### 2.4 Independence

By Poisson thinning the counts `{n(p)}` over distinct paths are **independent** Poissons, so they are legal
independent factors. One rule follows, and it is the rule the multi-junction bug violated:

> **Consume a partition, never a partition and one of its own marginals.**
> `{n(p) : a(p) = A}` is a partition. Junction flux `φ(J) = Σ_{p ∋ J} n(p)` and boundary crossing counts
> are **marginals** — legal to use *instead of* the partition, never *alongside* it, and their
> cross-junction covariance is positive.

---

## 3. IDENTIFIABILITY — what the path key actually buys

Prototype `scratchpad/acc_proto.py`, direct simulation of §2's generative model, exact analytic design
matrix, non-negative Poisson MLE. **No strand information, no prior, no imputation.**

### 3.1 Two components are separable from path counts alone, and the signal IS the FL contrast

Centre region of a uniform chain, gDNA `ρ=0.02` + RNA `ρ=0.05`, true `f_g = 0.286`:

| region | gDNA FL | RNA FL | design corr | Fisher cond | `f_g` est | bias |
|---|---|---|---|---|---|---|
| 150 bp | 60 | **200** | +0.33 | **2.5** | 0.287 ± 0.21 | **+0.001** |
| 150 bp | 60 | 150 | +0.52 | 3.4 | 0.303 | +0.018 |
| 150 bp | 60 | 100 | +0.87 | 14.2 | 0.477 | +0.191 |
| 150 bp | 60 | 65 | +1.00 | 876 | 0.618 | +0.332 |
| 150 bp | 60 | **60** | +1.000 | **8e16 (singular)** | 0.500 | +0.214 |

**The owner's cfRNA regime (gDNA 60 / RNA 200) is identified with a condition number of 2.5** — that is a
well-posed 2-component deconvolution from geometry × fragment length, prior-free and strand-free. It
degrades smoothly to exactly singular when the two FL distributions coincide. That is the honest
statement of the mechanism: **the information is the FL contrast, and nothing else.**

### 3.2 It is strongest where the current contained channel is destroyed

| region length | paths used | Fisher cond | bias |
|---|---|---|---|
| 60 bp | 10 | **2.4** | +0.058 |
| 100 bp | 7 | **2.3** | +0.009 |
| 150 bp | 5 | **2.5** | −0.035 |
| 300 bp | 3 | 5.1 | −0.039 |
| 600 bp | 2 | 14.7 | −0.007 |
| 1500 bp | 2 | **43.9** | +0.107 |

Short regions are *better* conditioned, because a short region resolves length finely. The current
contained channel is the exact opposite: at the median human exon region (147 bp) it retains **7.5 %** of
the region's opportunity and 48.7 % of exon regions have fewer than 10 contained start positions (§10).
**The two channels are complementary, and the redesign supplies the one that works at seams.**

⚠ This is not a brand-new information source. The current contained/left/right triple already has
different `E_G/E_R` ratios, so the contrast exists in principle today — but it is fractional (so the
Poisson likelihood is wrong), it collapses at short regions, and **no solver term uses it**: `node_init`'s
four sources are measured/intergenic, intron-factory, strand-deconv, unsolved-default. The redesign makes
an existing-but-unusable channel clean, integral, and available exactly where it is strongest.

### 3.3 The terminus bias, and that the correction removes it

Transcript `[2000, 6000)`, TES at 6000, 100 bp regions, true `f_g = 0.286`:

| distance to TES | `E_M(A)/\|A\|` | `f_g` naive (`E_M = \|A\|`) | `f_g` with admissible-start `E` |
|---|---|---|---|
| ≥ 300 bp | 1.000 | 0.28–0.35 | 0.28–0.35 |
| 200 bp | 0.810 | 0.279 | 0.249 |
| 100 bp | 0.200 | **0.512** | 0.409 |
| 0 bp (last region) | 0.004 | **0.645** | 0.422 |

**The owner's objection is confirmed and quantified: the naive rule over-calls `f_g` by up to +0.36.** The
admissible-start correction removes most of it; the residual at the very last region is genuine
low-information (`E_M/|A| = 0.004` — almost no RNA molecule can anchor there), and the Fisher information
reports it as such rather than hiding it.

---

## 4. DATA MODEL — one table

```
    paths        :  path_id -> { anchor_region:int32, end_region:int32,
                                 n_junctions:uint8, junction_ids:[...] }     (interned, CSR)
    path_count   :  uint32[n_paths, 4]                                       (the ONLY accumulator)
```

`channel` is unchanged: `ch = (spliced?2:0) + (primary?0:1)`.

Interning: per-worker open-addressing map, merged at the end exactly as `Accumulator::merge_from` does
today. Key = `(anchor, end)` packed into 64 bits for unspliced paths (exact, collision-free); spliced
paths key on the junction-id vector with exact comparison.

**Everything else is a projection** (§7) computed on demand from `path_count` — there is no second
accumulator, no `mass_left/right`, no `flux_left/right`, no float.

Memory is the one open number: bounded by the count of *distinct observed* paths. **M1 gates it** (§10).

---

## 5. DEPOSIT

```
deposit(fragment):
    spans   = fragment_genomic_spans(blocks, cut_introns)     # already shipped
    path    = walk(spans)                                     # regions touched + junctions used
    if path is unique:
        path_count[intern(path)][channel] += 1
    else:                                                     # >1 compatible path
        ambiguous_buffer.append(candidate_path_ids, channel)  # §8
```

`walk` is a single monotone sweep over the sorted region interfaces — the same `region_of_pos` +
advance loop the current `deposit` already runs, minus the slice/mass arithmetic.

---

## 6. EFFECTIVE LENGTHS — the four things to implement

All are pure functions of (graph geometry, FL pmf); each is unit-testable against brute-force enumeration,
and **none depends on the deposit rule** — that is the point.

1. `reach_c(s)` per component — a prefix quantity over the region, from the graph's TSS/TES/junction typing.
2. `E_c(A) = Σ_{s∈A∩S_c} F_c(reach_c(s))` — E2, the node effective length.
3. `E_c(p)` — E1, per path. For unspliced `p = (A,Z)`:
   `E_c(A,Z) = Σ_w f_c(w) · |A ∩ (Z − w + 1)|`, closed form in the region offsets.
   For spliced, the same convolution in transcript coordinates along `p`'s junction list.
4. The reference/support-end truncation, which falls out of `reach` with no special case.

---

## 7. PROJECTIONS

| consumer view | from `path_count` | kind |
|---|---|---|
| node (anchored) counts `n(A)` | `Σ_{a(p)=A} n(p)` | **partition** ✅ |
| cell counts `n(A,Z)` | `Σ_{a(p)=A, z(p)=Z} n(p)` | **partition** of `n(A)` ✅ |
| junction flux `φ(J)` | `Σ_{p ∋ J} n(p)` | marginal ⚠ correlated |
| boundary crossing count | `Σ_{p crosses b} n(p)` | marginal ⚠ correlated |
| **today's `AccumulatorPayload`** | **not derivable from paths alone** — the fractional shares need base coordinates. Computed by **dual-write at deposit time**, from the same spans, so it is bit-identical by construction. | legacy |

⚠ On the legacy view: the current solver's boundary nodes consume crossing flux, which is a *marginal*,
so it reuses one fragment across several boundaries. Dual-write reproduces that exactly (that is the
migration gate), and moving the solver onto anchored partitions is a later, separately measured arm — not
something the accumulator can fix on its own.

---

## 8. AMBIGUITY — the small side buffer (owner decisions D4/D5)

Two sources of a fragment having several compatible paths: an **implicit splice** whose mate gap matches
more than one annotated intron, and a **multimapper**. Both are resolvable once the RNA FL model is
trained, i.e. at the end of the scan.

* Unambiguous fragments — the overwhelming majority — deposit **directly during the scan**. No change to
  the scan architecture.
* Ambiguous ones go to a small separate buffer holding only `(candidate path_ids, channel)`. It is not
  `FragmentBuffer` and it stores nothing `FragmentBuffer` stores.
* Pass 2 (over that buffer only) assigns each fragment using the fitted FL pmfs and the paths' own counts.

⚠ **Assignment must stay integral.** Fractional assignment reintroduces exactly the non-integer observable
this redesign removes. Ship **hard assignment to the maximum-likelihood path**; an EM with fractional
weights is a later measured arm, and if it is ever adopted the variance must stop being `1/n`.

**M2 gates this**: measure the ambiguous fraction on real cfRNA before sizing the buffer.

---

## 9. DELETIONS

`boundary_mass_left/right` (8 float32/boundary) · `boundary_flux_left/right` (8 uint32/boundary) · the
`mass`/`n` duality through `substrate.py` and `node_geometry.py` · `boundary_side_eff_length`'s ½ (derived
from the deposit rule) · the continues/terminates selector for the spliced divisor (~40 lines in
`node_geometry.py` reconstructing junction/terminus structure from signature bits) ·
`boundary_junction_strand` + its `0`-is-identity worker merge (becomes index-time) · the `_EPS` floors on
every effective length · `region_types`/`fl_pool_idx` plumbing through `AccumulatorSet`.

---

## 10. MEASURED FACTS (real human index `refs/rigel_index_v7`, 752,654 regions, 227,844 transcripts)

| | |
|---|---|
| median exon region | **147 bp** |
| contained-channel efficiency `E/L` at exon regions, cfRNA-like FL | **0.075** (gDNA-like FL: **0.011**) |
| exon regions with fewer than 10 contained start positions | **48.7 %** |
| exonic starts crossing ≥1 interface | **37.7 %** |
| … of those crossers, crossing more than one | **12.97 %** (toy: 3.38 % — **3.8× understated**) |
| transcript termini strictly **inside** a region | **59.5 %** (383,187 total) |
| regions added by making termini partition events | **+30.3 %** (→ 980,508) |
| **exonic bp where `F_M(reach) < 0.99`** (FL 200) | **16.9 %** |
| **mean `E_M/\|A\|` over all exonic bp** (FL 200) | **0.893** — a 10.7 % systematic gDNA over-call if uncorrected |
| transcripts wholly inside the taper zone (FL 200) | 2.2 % |

⚠ The 10 Mb synthetic overstates contained efficiency **4.6×** (median exon region 270 bp vs 147) and
understates multi-crossing **3.8×**. It is also Poisson by construction. **Gate on real cfRNA.**

**Open, to measure before locking the spec:**

```
M1  distinct-path census + memory on cfrna:LBX0190          (gates §4)
M2  ambiguous-fragment fraction on real cfRNA               (gates §8)
M3  deposit cost: path interning vs the current two array
    increments, at genome scale                             (gates §5; profile on cfRNA, NEVER the toy)
M4  isoform sensitivity of reach_M's maximal-path convention (gates §2.2's ⚠)
M5  per-node FL histogram value-add, given §3 already
    extracts the FL contrast geometrically                   (gates D3 — may be redundant)
```

⭐ **M5 is new and matters**: §3 shows the path key already performs an FL deconvolution. A per-node FL
histogram (D3) may be largely redundant with it at short regions and additive only at long ones (where
§3.2 shows conditioning degrades, cond 43.9 at 1500 bp). Measure the overlap before wiring it.

---

## 11. IMPLEMENTATION PLAN

Owner decision **D6: do the full rework** — so P1g is not a separate workstream; it is P1 below.
(For the record, the distinction v1 drew: *P1g-lite* = restore the v6 boolean flags onto today's region
interfaces, partition unchanged, sees 40.5 % of termini; *P1g-full* = termini become partition events,
sees 100 %. The full rework subsumes both.)

```
P0  M1–M5, then lock this spec
P1  INDEX v8 — TSS/TES as partition events; typed edges (contiguous / junction, per strand, per side)
    gate: partition invariants; every index and every _selfsolve_cache payload rebuilt
P2  PATH ACCUMULATOR (C++), dual-writing today's payload from the same spans
    gate: legacy payload BIT-IDENTICAL 32/32 — it adds a store nobody reads yet
P3  EFFECTIVE LENGTHS — the four functions of §6, incl. admissible reach
    gate: unit tests vs brute-force enumeration; the §3.3 terminus table reproduced
P4  CONSUMERS migrate, one observable per A/B arm:
    node partition -> cells -> full paths.  Each arm re-records its own baseline.
P5  AMBIGUITY buffer (§8), hard assignment
P6  per-node FL histograms (D3) — LAST, and only if M5 says they add anything
P7  deletions (§9), then goldens
```

**Standing gates on every arm:** ruff; full suite; held-fixed `z2` must not regress; one thing varied per
arm; a pre-registered falsification test; baseline re-recorded from the current tree in the same session
at both refits (current reference **r0 0.078786 / r1 0.052470**).

⚠ **Sequencing against the calibration branch.** `src/` currently carries the uncommitted W4 landscape,
the additive `_gdna_arm` and the P-2 pin fix, with goldens already regenerated. P1 invalidates every cached
payload and every golden. **Get that tree accepted and committed before starting P1**, or the two change
sets become impossible to attribute.

⚠ **`/tmp/rigel_selfsolve` is shared and non-namespaced** — namespace it by index hash as part of P1;
it has already cost two full rebuilds.
