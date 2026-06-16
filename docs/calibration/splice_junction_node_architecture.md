# Splice-junction node architecture (bipartite boundary nodes)

**Status:** FUTURE / optional. Not a bug fix. Not needed for any open calibration
goal today (including the imminent pooled-seam effective-length PR). This document
records (1) the decisive verdict that spliced boundary one-sidedness is **not** a
correctness bug, and (2) a complete, implementable design for the bipartite
splice-junction node graph, to be built only when transcript-space mature-RNA
propagation becomes the active work item.

All file:line citations were verified against code (not docstrings) at
`main` around the displaced-mass fix (`17d0209`).

---

## 0. The observation under review

The calibration node graph is **linear today**: regions in genomic order, with a
boundary between each pair of genomically-adjacent regions. Each boundary stores
per-SIDE, per-CHANNEL mass — `mass_left[c]` / `mass_right[c]` for 4 channels
(`c0` unspliced+, `c1` unspliced−, `c2` spliced sense, `c3` spliced antisense)
(`tests/native/_accumulator_reference.py:62-79`; C++ `Boundary`,
`src/rigel/native/calibration/accumulator.cpp`). The substrate collapses these to
per-region left/right views with `mass_unspliced = m[:,0]+m[:,1]` and
`mass_spliced = m[:,2]+m[:,3]` (`src/rigel/calibration/substrate.py:52-53`).

The user's observation: merging boundary sides is correct for **unspliced** reads
(gDNA + unspliced nascent sit on *both* sides of a boundary) but **spliced** reads
are **one-sided** — a mature RNA splicing at an exon→intron boundary deposits mass
on the LEFT (exon) side and ZERO on the RIGHT (intron) side. The accumulator
already deposits spliced one-sided.

**This observation is correct and live in shipped code.** Proof: spec test T4
(`tests/native/test_accumulator_spec.py:153-179`), a spliced intron-skip
`blocks=[(1800,1950),(5050,5950)]` over partition exon1=[1000,2000),
intron=[2000,5000), exon2=[5000,6000):

```
b1.mass_left[ch]  == 150/L   b1.mass_right[ch] == 0.0   # exon1→intron: LEFT only
b1.flux_left[ch]  == 1       b1.flux_right[ch] == 0     # intron-facing zero
b2.mass_right[ch] == 900/L   b2.mass_left[ch]  == 0.0   # intron→exon2: RIGHT only
```

The deposit comment states the invariant outright (`accumulator.cpp:186-189`):
"a spliced jump credits one side of each flanking boundary, leaving the
intron-facing sides at zero (no false exon-intron flux)."

One sharp edge, not a bug: "spliced ⇒ one-sided" holds only at the seams
**flanking a skipped intron**. A contiguous spliced read whose exon body straddles
an internal region partition with *no* intron between (spec test T5,
`test_accumulator_spec.py:194-216`, same `spliced=True` flag) still splits its
interior slice 50/50 across both sides of that internal seam. The true invariant
is "one-sided **at intron-flanking junction seams**," guaranteed by the slice gap,
not "one-sided everywhere." This distinction is load-bearing for the design (§4).

---

## 1. IS IT A BUG?

**Verdict: NO. Spliced one-sidedness is not a correctness bug in any consumer
today, and the pooled-seam gDNA effective-length design is structurally immune.**
The user's reasoning is confirmed, and is in fact *stronger* than stated: the
immunity does not rest only on "the IPR is gDNA-only" — *every* consumer of the
spliced channel either treats one-sidedness as the intended signal, or never reads
boundary spliced at all.

The single invariant that protects everything:

> **gDNA mass is computed from `mass_unspliced` (c0/c1) only; spliced mass (c2/c3)
> is added wholesale to RNA and is pooled side-agnostically per region before any
> prior. No gDNA quantity ever reads the spliced channel.**

Per-consumer:

### 1a. Effective length / gDNA IPR — IMMUNE
The IPR sums `mass_gdna_*` only (`src/rigel/calibration/priors.py:222-226`), and
every gDNA quantity is `gdna = f_g · mass_unspl` by construction
(`src/rigel/calibration/simplex_sweep.py:289`;
`src/rigel/calibration/strand_deconv.py:215-216` — spliced flows only into
`rna = (1−frac)·mass_unspl + mass_spliced`). The capture-aware EM eff-lengths take
`mass_gdna_*` only (`src/rigel/calibration/capture_eff_length.py:127-134`); the
geometric length is pure region geometry (`src/rigel/calibration/derive.py:49-55`).
gDNA is set `spliced=true` only on a CIGAR-N intron
(`src/rigel/native/bam_scanner.cpp:1446-1447`), which a contiguous genomic read
cannot produce, so `mass_gdna_spliced ≡ 0` both by deposition and by readout.

**Under the pooled-seam design:** the seam pools `mass_gdna_right[r] +
mass_gdna_left[r+1]` (`effective_length_redesign_plan.md` §7). Both operands are
unspliced-derived gDNA; the spliced channel is not in the sum. Pooling the two
gDNA sides correctly reunifies one physical contiguous crossing (gDNA straddles a
seam → deposits both sides). **The pooled-seam eff-len cannot be affected by
spliced one-sidedness regardless of how spliced is routed.** The user is right that
this PR is safe.

### 1b. RNA prior `rna_prior_count` — IMMUNE (exact cancellation)
`rna_region = mass_rna_contained + mass_rna_left + mass_rna_right −
mass_rna_spliced` (`priors.py:233-239`). Whatever one-sided spliced mass was added
into `mass_rna_left`/`mass_rna_right` is subtracted by `mass_rna_spliced`, which is
the sum of `mass_spliced` over the **same three views**
(`src/rigel/calibration/calibrate.py:298-302`). The add and the subtract use the
identical per-side quantity, so the result is **independent of which side the
spliced mass landed on** — it would net to the same value even if spliced were
two-sided. This is an algebraic identity, the most robust kind of non-bug.

### 1c. Per-node deconv `deconv_sides` / `_deconv_per_node` — IMMUNE
The gDNA fraction uses `n_unspliced_pos/neg` only
(`strand_deconv.py:265-266`); the count fraction uses `mass_unspliced` only
(`strand_deconv.py:308-319`). Spliced enters solely as the deterministic additive
RNA term (`strand_deconv.py:216`). An intron-facing side simply has
`mass_spliced=0` and adds zero RNA there — which is correct, there is no mature RNA
on the intron side.

### 1d. `splice_junction.region_splice_gdna_frac` — one-sidedness is LOAD-BEARING (and correct)
This is the one place where one-sidedness is the *intended mechanism*. The
eligibility predicate classifies an adjacent region pair as a junction and names
the **exon side** (`exon_side="R"` / `"L"`); the anchor reads spliced flux only
from that exon side (`src/rigel/calibration/splice_junction.py:246-275`), gated by
`_anchor`'s `spl > 0` requirement (`:250`). It *expects* the intron-facing side to
be zero and never queries it. If spliced were two-sided, the intron region would
falsely read a mature reference and mis-debias. **Merging the sides naively would
break this.** Not a bug today — and the bipartite design (§3) makes this dependency
*structural* rather than latent.

### 1e. `rna_variance` / `mature_density` — one-sidedness is LOAD-BEARING (and correct)
Same exon-side anchoring (`src/rigel/calibration/rna_variance.py:55-63`;
`src/rigel/calibration/mature_density.py:94-119`). The `(d_L, d_R)` var~mean pair
is two exon-side measurements of the same mature density; intron-facing zeros are
never read. Two-sided spliced would inject spurious mature signal at the wrong
locations and corrupt the variance fit.

### 1f. Simplex sweep sided spliced lower bound — NOT A CONSUMER of boundary spliced
The sweep's `spliced_pos`/`spliced_neg` come from `c.n_spliced_sense` where
`c = substrate.contained` (`src/rigel/calibration/simplex_sweep.py:211,217-220`) —
the **contained** spliced count, never the boundary sides. It is a per-region RNA
floor (`f_s·U ≥ spliced_s`, `:79-82`). Boundary one-sidedness is structurally
irrelevant here.

### 1g. `derive` — IMMUNE
Sums `gdna_mass` only (`derive.py:49-51`). Spliced never reaches it.

### Could a bug be obscured?
The single latent risk: a **future** consumer that pools the RNA/spliced channel
across a *genomic* seam — e.g. `mass_rna_right[r] + mass_rna_left[r+1]` for an RNA
seam IPR. That seam would receive RNA from only the exon-flanking region and zero
from the intron-flanking region, giving an asymmetric half-populated seam. It is
not a bug today only because **no consumer pools the RNA/spliced channel across a
genomic seam**; every spliced consumer is exon-side-anchored. This is a latent
constraint, not a free invariant — every future boundary-mass consumer must respect
it. The bipartite design enforces it structurally. **Recommended cheap insurance:
Phase 0 below pins this constraint with a regression test now, independent of
whether the full rebuild ever happens.**

---

## 2. RECOMMENDATION

**DEFER. Do not implement the bipartite redesign now. Do exactly one thing today:
Phase 0 (§5) — a ~30-minute regression test pinning the latent one-sidedness
constraint.**

Rationale:

- **Zero current correctness benefit.** §1 proves every consumer is immune or
  one-sidedness-dependent-and-correct. No calibration metric moves. The pooled-seam
  eff-len — the imminent PR — is provably immune (§1a) and needs none of this.
- **High, real cost.** The accumulator change is C++ + a byte-for-byte reference
  rewrite + a `deposit` signature change + a payload schema reshape + every
  spliced-bearing accumulator golden moving. That is a large regression surface for
  no current gain.
- **The headline justification (isoform-swap) is mis-attributed and unverified.**
  The "68.5% of aggregate error is isoform-swap" figure is *end-to-end* error.
  Isoform-swap — *which* transcript within a gene — is resolved in the **per-locus
  EM** distributing RNA mass among compatible transcripts, *after* the calibration
  prior sets the gDNA-vs-RNA split. The calibration prior is **two per-locus
  Dirichlet scalars** (`alpha_rna_add`, `alpha_gdna_add`); it is constitutionally
  incapable of expressing which isoform. So transcript-space propagation in
  *calibration* cannot directly touch isoform-swap. At best it sharpens the
  **mature↔nascent split** (spliced-RNA vs unspliced-RNA), a different axis. The
  prior analyses conflated the two. The flagship is currently **EM-bound, not
  calibration-bound, at high strand-specificity** (memory:
  `flagship_em_bound_not_calibration`). Building a calibration-side propagation
  graph to attack an EM-side error is, on current evidence, solving the wrong
  problem at the most expensive layer.

**What would change the call** (all must hold before building):

1. The dissect tool (`scripts/debug/dissect_loci.py`) decomposes the aggregate
   error and shows a material slice is **mature↔nascent split that is
   calibration-reachable** — not EM-side isoform-swap, not gDNA leak. Threshold:
   the cross-intron coupling must plausibly recover more than the gDNA-strand /
   eff-len EM fixes already ranked #1–2 on the roadmap.
2. A **toy experiment** (long intron, depleted internal exon — the
   `calibration_feature_dev_toy_scenarios` recipe) demonstrates the linear chain
   *fails* and a hand-wired junction edge *fixes* it, before any accumulator cost
   is paid.
3. The EM-side gDNA-strand fix and gDNA-component eff-len (roadmap #1–2) have
   **already landed**, so the flagship is no longer EM-bound and the residual is
   provably in the calibration prior.
4. (Strongest amortization argument) the junction-node substrate is **also** wanted
   by the EM's transcript-space resolution — making the C++ cost amortize across
   two initiatives rather than one.

If 1+2+3 hold, build Phase A→C1 as inert plumbing, then Phase D as a guarded
research item.

---

## 3. THE DESIGN

### 3.1 Two node families (bipartite)

| Node kind | Status | Carries | Adjacency | Sidedness |
|---|---|---|---|---|
| **Region** | reused | contained mass | — | — |
| **Unspliced boundary (seam)** | reused, narrowed to c0/c1 | gDNA + unspliced nascent | **genomic**-adjacent: `region r ↔ seam ↔ region r+1` | two-sided (a contiguous crossing straddles both flanks) |
| **Splice-junction (SJ)** | **NEW**, sparse | spliced c2/c3 | **transcript**-adjacent: `donor-exon region ↔ junction ↔ acceptor-exon region` (intron region NOT adjacent) | one-sided by physics (mature on the exon flank) |

"All boundaries become one-sided (store one kind of mass)" is realized as: **each
node family stores exactly one mass kind.** The unspliced seam loses c2/c3 (drops
from 4 channels to 2). The SJ node is the only carrier of c2/c3. The graph becomes
**bipartite**: region ⊥ {seam ∪ SJ}; travelling between two region nodes always
visits a boundary-family node.

Caveat (the T5 sharp edge, §0): a *contiguous* spliced exon body that straddles an
internal region partition with no intron in the gap is **not** a junction. Its mass
is genomic-contiguous mature mass. The honest invariant is "SJ ⟺ intron-flanking
seam"; contiguous-spliced crossing mass needs an explicit home (§4, ATTACK 2).

### 3.2 Accumulator + payload change (C++ — the load-bearing work)

**The junction CANNOT be keyed from the slice gap. The `deposit` signature must
change to pass the cut-intron coordinates.** This is the single most important
correction to the earlier design sketches, verified against code:

- `deposit(block_starts, block_ends, n_blocks, spliced, primary)`
  (`accumulator.cpp:84-88`) receives **only post-intron-subtraction genomic
  spans** — `fragment_genomic_spans` (`bam_scanner.cpp:732-776`) has already cut
  the introns and discarded their coordinates before the call
  (`bam_scanner.cpp:1471-1483`). `deposit` rebuilds slices purely from blocks
  clipped to region edges (`accumulator.cpp:102-119`).
- The proposed discriminant "spliced AND slice-gap" is **not** a junction detector.
  Counterexample T5 (`test_accumulator_spec.py:194-216`): `spliced=True`, blocks
  `[(0,100),(100,180),(200,320)]`, contiguous R0→R1→R2 — **no gap**, yet spliced.
  The slice-gap discriminant misclassifies it. The gap is an artifact of partition
  *granularity* (introns become regions only where some transcript declares them,
  and `regions.py` merges adjacent same-signature segments), not of the splice.
- The junction identity is the **SJ tag `(intron_start, intron_end, strand)`**,
  which the scanner has (`cut_introns`, `cr.sj_strand`,
  `bam_scanner.cpp:1471-1472`) and throws away at the deposit boundary.

**Required signature:**
```cpp
void deposit(const std::int64_t* block_starts,
             const std::int64_t* block_ends,
             std::size_t n_blocks,
             const std::int64_t* intron_starts,   // NEW: cut introns
             const std::int64_t* intron_ends,     // NEW
             std::size_t n_introns,                // NEW (0 for unspliced)
             bool spliced,
             bool primary);
```
The caller (`bam_scanner.cpp:1471-1483`) already has `cut_introns`; pass them via
two new scratch vectors alongside the existing span scratch. With intron coords in
hand: the SJ key is the literal `(start,end,strand)`; the donor/acceptor *regions*
are looked up by `searchsorted(side="right")` on the intron endpoints — reusing the
exact containment idiom of `capture_eff_length._exon_region_incidence` to avoid the
off-by-one that bit the nascent eff-len (memory: `nascent_efflen_invariant_bug`);
and T4 (true junction, intron coord present in the gap) is cleanly distinguished
from T5 (contiguous spliced body, no intron coord in the gap).

**New junction store** (sparse, sorted-vector keyed by packed
`(intron_start, intron_end, strand)` per `Accumulator`, analogous to the global
nRNA-span table; `merge_from` extends to merge two sparse junction maps by key):
```cpp
struct SpliceJunction {
    std::int64_t intron_start, intron_end;  // the SJ key
    std::int32_t donor_region, acceptor_region;
    float         mass_sense[1], mass_antisense[1];   // c2/c3
    std::uint32_t flux_sense[1], flux_antisense[1];
};
```

**Routing** in the crossing path:
- c0/c1 unspliced contiguous crossing → genomic `Boundary` (now 2-channel),
  two-sided, **bit-identical to today's c0/c1**.
- c2/c3 across a cut intron → SJ node keyed by that intron; donor-end slice credits
  `mass_sense/antisense` (donor facet), acceptor-end slice credits the acceptor
  facet. **The two flanking genomic boundaries' c2/c3 are no longer touched.**
- c2/c3 contiguous straddle (T5, no intron in the gap) → see §4 ATTACK 2: must have
  an explicit destination (a separate contiguous-spliced seam store, or attribution
  to the longer slice's contained-spliced). It is crossing mass, NOT contained.

**Payload schema** (`src/rigel/scan_payload.py:57-74`): add a junction CSR table —
`ref_junction_offsets[n_refs+1]`, `junction_donor_region[]`,
`junction_acceptor_region[]`, `junction_intron_start[]`, `junction_intron_end[]`,
`junction_mass_sense[]`, `junction_mass_antisense[]`, `junction_flux_*[]`. Narrow
`boundary_mass_*`/`boundary_flux_*` to 2 channels (or keep at 4 with c2/c3 ≡ 0 for
one transition release behind a flag). The
`R_total + nonempty_refs == B_obj_total` invariant
(`scan_payload.py:138-143`) is unchanged (boundary cardinality unchanged).

**Substrate** (`src/rigel/calibration/substrate.py`): boundary left/right views read
c0/c1 only — `mass_unspliced` is *already* `m[:,0]+m[:,1]` so `_make_view` is
unchanged; only `mass_spliced`/`n_spliced_*` on the boundary views drop to 0. Add a
fourth `junctions` view keyed by region-pair exposing
`n_spliced_sense/antisense` + `mass_spliced`. Contained spliced is untouched
(it never went through a boundary), so the sweep's lower bound (§1f) is unaffected.

### 3.3 Does belief propagation change?

**Yes — the prior framing "BP does not change" is false as written, recoverable in
spirit.** The grid sum-product *kernel* is fully reusable: `_edge_logphi`
(`simplex_sweep.py:130-140`, the `(P,P)` log-coupling), `_local_loglik` /
`_mixture_strand_loglik`, the `_fg_median`/`_fg_var` readout, and the
decoupled-`O(P)` short-circuit + `(q_pos,q_neg)`-keyed cache (`:160-172`). But the
**driver** is a literal chain: `_sweep_chain` is `for i in range(1, m)` over one
per-reference ascending index array with edges implicit as `(idx[k], idx[k+1])`
(`simplex_sweep.py:143-188, 275-276`); edges carry finite `q_rna` only when both
endpoints are exon on the same strand (`exon_pos[a] & exon_pos[b]`, `:277-280`),
which decouples across introns.

Adding SJ edges (transcript-adjacent, `Q=q_rna`) makes the graph non-linear. The
driver must become a **tree message schedule** (`_sweep_tree`): build a per-component
adjacency list of genomic seam edges + junction edges, then leaf→root collect /
root→leaf distribute. For a *single* transcript this is a path → still a tree →
two-sweep exact, identical guarantee to today. **Alternative splicing breaks
tree-ness** (cassette exons → cycles) → loopy BP or junction-tree; this is the real
PR-D cost (§4 ATTACK 3). The math is unchanged; the iteration schedule is not, and
on a cyclic graph chain-exact sum-product ≠ loopy BP.

Recommended first cut: **junction node as a pure pass-through coupling**
(`ψ_junction ≡ 0`), so it only carries the `Q=q_rna` edge between donor and acceptor
regions. This delivers the transcript-space coupling with zero conservation risk and
defers the junction-as-evidence question.

### 3.4 The maps

- **gDNA → node:** unchanged. gDNA is genomic; it lives on region-contained (c0/c1)
  and genomic seams. No gDNA ever touches an SJ node (gDNA never splices). The
  pooled-seam IPR pools genomic-seam gDNA only.
- **Transcript → node:** the new map. A transcript is a path through region nodes
  connected by its own SJ nodes; the SJ key `(donor_region, acceptor_region)` is
  exactly a transcript's exon-exon adjacency. Reconstructable from `index.t_df`
  exon structure intersected with the region partition via the same
  `searchsorted(side="right")` incidence helper.

### 3.5 Conservation proof

**Mass (Law M).** Today's invariant
(`_accumulator_reference.py:254-263`):
`T = Σ_r Σ_ch contained[r,ch] + Σ_b Σ_ch (mass_left[b,ch] + mass_right[b,ch])`.
`ch` is fixed per fragment (`accumulator.cpp:92`), so every deposited gram lives in
exactly one of {c0,c1} XOR {c2,c3}. Partition the boundary double-sum by channel:
the unspliced term is byte-identical to the new 2-channel `boundary_mass_*`; the
spliced term is re-keyed from `(boundary, side)` to `(SJ node, facet)` by a
**bijection on deposited mass** (each spliced quantum has exactly one destination,
because c2/c3 are physically *removed* from the boundary store — there is no slot to
double-write). Hence `T' ≡ T` term-by-term, and the new
`total_mass_deposited()` (the three-family sum) must still equal `n_fragments` —
this is the byte-for-byte gate.

**Caveat the proof must honor (§4 ATTACK 2):** for a read that is BOTH spliced AND a
contiguous crossing (T9-shaped: an interior exon slice bounded by a junction edge on
one side and a contiguous edge on the other), you cannot simultaneously have
(a) seams carry only c0/c1, (b) total conserved, and (c) per-region
`mass_rna_spliced` bit-identical. Pick the resolution explicitly: route
contiguous-spliced crossing mass to a dedicated contiguous-spliced store and
re-derive `mass_rna_spliced` from {SJ nodes + that store}. The conservation proof
must be re-derived with the contiguous-spliced term explicit and validated on a
scanned T9 payload, not asserted.

**Length (Law L).** SJ nodes are **transcript-space**: they represent an intron
skip, occupy **zero genomic bp**, are NOT inserted into `boundary_positions`, and do
not alter `ref_region_offsets`/`ref_boundary_offsets`. Every length quantity reads
region-partition geometry (`derive.py:49-55`; `result.py`) or gDNA mass
(`priors.py:222-226`) only — SJ nodes carry neither. So `Σ region genomic bp` is
invariant and no length quantity reads an SJ node. Law L is preserved by
construction.

---

## 4. RISKS (adversarial findings)

| # | Risk | Severity | Mitigation |
|---|---|---|---|
| 1 | **SJ key from slice-gap is wrong.** `deposit` has no intron coords (`accumulator.cpp:84-88`); the slice-gap discriminant misclassifies contiguous-spliced T5 (`test_accumulator_spec.py:194-216`) and merge-fused introns. | **FATAL as naively designed** | **Change the `deposit` signature** to pass cut-intron coords (§3.2). Key SJ by `(intron_start,intron_end,strand)`; look up donor/acceptor regions by `searchsorted`. This distinguishes T4 (true junction) from T5 (contiguous). Cost: signature change ripples to every caller + byte-for-byte reference + every accumulator test. The design under-scoped this; it is the expensive part. |
| 2 | **Double-count / loss for spliced-AND-crossing (T9).** Cannot have seams-one-kind + total-conserved + per-region-spliced-bit-identical at once. | **SERIOUS** | Route per slice-edge: an edge spanning a cut intron → SJ node; a contiguous edge (no intron in gap) → a dedicated contiguous-spliced store (NOT contained, it crosses a boundary). Re-derive conservation with the contiguous-spliced term explicit. Validate `test_mass_conserved_per_node` (`tests/calibration/test_calibrate.py:42-46`) on a scanned T9 payload, not by assertion. PR-A/B are inert; PR-C is NOT (see #4). |
| 3 | **BP on a non-tree.** Cassette exons → cycles (isoform A–B–D vs A–C–D ⇒ 4-cycle). Region partition unions SJ edges over all overlapping isoforms (`regions.py` sweeps all transcripts), so cycles are the **common** case in multi-isoform loci, not an edge case. Chain-exact sum-product ≠ loopy BP. | **SERIOUS (PR-D only)** | PR-D is a research item: either a loopy-BP convergence study on these Gaussian-coupling log-odds potentials (with damping/schedule + a convergence guarantee that does not exist for free), or a junction-tree (exponential in treewidth). **Never** break cycles by dropping SJ edges — that re-isolates the cassette exon the feature targets. Collides with the already-open iteration-convergence problem (CLAUDE.md "iteration tightening") + the 4-pass bootstrap; risks non-determinism. Defer PR-D explicitly. |
| 4 | **"Inert refactor" is false for PR-C.** Rewiring `splice_junction`/`rna_variance`/`mature_density` to read SJ nodes by region-pair makes **exon-skip junctions** (non-genomically-adjacent regions) newly visible — they cannot be seen today via `splice_junction_eligibility(sig[i-1],sig[i])`. This is a behavior change that moves calibration golden outputs and could regress the fragile weak-SS(0.65) phantom-gDNA tolerance (~7.6%, memory: `calibration_pipeline_convergence_state`). | **SERIOUS** | **Split PR-C.** C1: re-key, restrict SJ visibility to genomically-adjacent junctions → assert golden-identical (true inert refactor). C2 (separate, benchmarked): admit exon-skip junctions. Do NOT bundle the visibility change into the plumbing refactor. |
| 5 | **Spliced-contained vs sweep lower bound.** If T9 contiguous-spliced crossing mass is mis-routed into `contained_spl`, the sweep's lower bound `spliced/n` (`simplex_sweep.py:79-82,211,217-220`) inflates and over-forces RNA. | **MANAGEABLE** | Keep the contained/crossing partition (`accumulator.cpp:160-173` `all_same` test) authoritative and orthogonal to the SJ split. Contiguous-spliced crossing mass is crossing mass, never contained. |
| 6 | gDNA eff-len immunity | **UNFOUNDED to attack — confirmed** | No spliced quantum can enter any gDNA channel (`simplex_sweep.py:289`; `strand_deconv.py:215`; `bam_scanner.cpp:1446-1447`). Pooled-seam IPR immune regardless of routing. This is why the feature buys zero correctness today. |
| 7 | Junction table is data-dependent sizing | minor | Sorted-vector-per-ref + CSR offsets, matching the nRNA-span-table pattern; `merge_from` merges by key. |

---

## 5. IMPLEMENTATION PLAN (FUTURE / optional)

**Sequence is non-negotiable: ship pooled-seam eff-len (independent, immune) and
land the EM-side gDNA-strand / eff-len fixes (roadmap #1–2) FIRST.** Only then,
and only if §2's "what would change the call" conditions hold, build A→D.

### Phase 0 — guard the latent invariant (DO NOW, ~30 min, no architecture change)
- Add a regression test asserting `region_splice_gdna_frac`, `rna_variance`,
  `mature_density` read spliced flux **only** from the predicate-named exon side;
  and that `total_mass_deposited() == n_fragments` continues to hold.
- Converts the undocumented C++↔Python one-sidedness coupling into a guarded
  contract; protects against a future contributor "fixing" the accumulator to
  two-sided spliced and silently corrupting the splice-junction debias. This is the
  only part of the whole initiative with positive cost/benefit today.

### Phase A — accumulator junction routing (C++ + reference, byte-for-byte) — EXPENSIVE, flagged
- **Files:** `src/rigel/native/calibration/accumulator.cpp` (+ `accumulator.h`),
  `src/rigel/native/bam_scanner.cpp` (pass `cut_introns` to `deposit`),
  `tests/native/_accumulator_reference.py` (mirror, byte-for-byte).
- Change `deposit` signature to take cut-intron coords (§3.2, ATTACK 1). Narrow
  `Boundary` to c0/c1. Add `SpliceJunction` sorted-vector store + `merge_from`. Add
  a contiguous-spliced store (ATTACK 2). Route by intron-coord, not slice-gap.
- **Gate behind a flag** so default = byte-identical.
- **Tests:** new spec tests — T4′ (intron-skip → SJ node; flanking boundaries c2/c3
  == 0), T9′ (multi-junction + contiguous mix; conservation), T5′ (contiguous
  spliced → contiguous-spliced store, NOT SJ). Extend `total_mass_deposited` to the
  three-family sum and keep all `== n_fragments` assertions
  (`test_accumulator_spec.py:78,179,216,252,280,311`). T4/T9 **move** — land as an
  explicit `--update-golden` commit so a downstream benchmark shift bisects cleanly.

### Phase B — payload + substrate
- **Files:** `src/rigel/scan_payload.py` (junction CSR table), the scanner emit
  side, `src/rigel/calibration/substrate.py` (junctions view; boundary views c0/c1).
- **Tests:** extend `tests/calibration/test_substrate_conservation.py` to assert the
  re-keying bijection (`Σ boundary[0:2] + Σ sj_mass == Σ non-terminal old-boundary
  mass`).

### Phase C — rewire spliced consumers (split into C1/C2 per ATTACK 4)
- **C1 (inert):** `src/rigel/calibration/calibrate.py:298-302` re-sources
  `mass_rna_spliced` from SJ nodes (donor+acceptor attribution → bit-identical
  per-region totals); `splice_junction.py`, `rna_variance.py`, `mature_density.py`
  read SJ nodes by region-pair key, **restricted to genomically-adjacent
  junctions**. Gate: `tests/calibration/test_calibrate.py:42-46`
  (`test_mass_conserved_per_node`) and the calibration golden outputs must be
  **unchanged** — the proof of inertness. Add `test_sj_attribution_inert`.
- **C2 (behavior change, benchmarked):** admit exon-skip junctions. Validate via the
  `calibration-benchmark` skill (net fragment-flow); watch the weak-SS phantom
  tolerance.

### Phase D — transcript-space mature propagation (the payoff; SEPARATE research initiative)
- **Files:** `src/rigel/calibration/simplex_sweep.py` — replace `_sweep_chain` with
  `_sweep_tree`; adjacency = genomic seams + SJ edges (`Q=q_rna`). For cyclic
  components: loopy BP (with a convergence-study deliverable) or junction-tree;
  never drop SJ edges to break cycles.
- **Validation:** develop on a toy first (small genome, multi-exon transcript with a
  long intron + an overlapping opposite-strand AMBIG exon — memory:
  `calibration_feature_dev_toy_scenarios`). Then `calibration-benchmark` net
  fragment-flow targeting `net_gdna_to_{nrna,mrna}` and the mature↔nascent metric on
  the `ss_0.99_nrna_rnd_capture_on` flagship.

---

## 6. WHAT IT UNLOCKS

The current chain couples the per-strand RNA:gDNA log-odds (the only
enrichment-and-content-invariant signal) only along **genomically-adjacent**
same-strand exon stretches (`simplex_sweep.py:277-280`). It explicitly **decouples
across an intron**: a transcript's exon-1 and exon-2, separated by a long intron,
sit on a decoupled chain edge (the intron is not an exon) — their odds cannot inform
each other, and the accumulator never even records the `(exon1, exon2)` adjacency.

The SJ node adds exactly that missing edge: a finite-`Q` coupling between
transcript-adjacent exons across the skipped intron. This enables **transcript-space
mature-RNA propagation** — a depleted internal exon can borrow mature odds from its
spliced-connected neighbours.

**Honest scope of the payoff.** This sharpens the **mature↔nascent split**
(spliced-mature vs unspliced-nascent RNA), which is what calibration *can* express.
It does **not** directly resolve isoform-swap (which isoform within a gene) — that
lives in the per-locus EM, downstream of calibration's two Dirichlet scalars (§2).
The frequently-cited "68.5% isoform-swap" figure is therefore the wrong target for
this feature; the right target is whatever fraction of the mature↔nascent confusion
is calibration-reachable, which must be measured (dissect tool) before committing.
The substrate may *also* serve the EM's transcript-space resolution later — if so,
that shared use is the strongest argument for paying the accumulator cost once.

**Bottom line:** not a bug, not needed now, provably immune under the pooled-seam
eff-len. The design above is complete and implementable for the future, with the
three corrections the earlier sketches got wrong made explicit: the junction must be
keyed from cut-intron coordinates (not the slice gap), PR-C is a behavior change not
an inert refactor, and BP on the cyclic alt-splice graph is loopy/junction-tree (not
"unchanged"). Build it only when transcript-space mature propagation is greenlit,
behind the EM-side fixes already prioritized.
