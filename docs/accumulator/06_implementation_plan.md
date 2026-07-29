# Accumulator v5 — the implementation plan

> ## ⛔ AMENDED 2026-07-29 BY OWNER DIRECTIVE — READ THIS BEFORE THE BODY
>
> **No backwards compatibility. No legacy retention. Converge and delete.** The tool is alpha; the
> production release in the community is a separate artifact. There is therefore **no reason to keep
> an old path alive for comparison** beyond checking that a change has not introduced a catastrophic
> bug — and keeping one alive is itself the risk, because the two become hard to disambiguate later.
>
> Three parts of the body below are superseded:
>
> 1. ⛔ **§3 W4's SHADOW MODE IS CANCELLED.** Emitting both payloads from one scan (~1.6 GB peak RSS,
>    two payload classes, `scan_payload_v5.py` beside `scan_payload.py`, every consumer written twice)
>    existed only to give each W5 arm a legacy A/B. **W3 is unchanged** — Python reference + the
>    byte-for-byte A1–A14 / E1–E4 / S1–S2 tests, written first and failing. W4 then replaces the C++
>    accumulator **outright**, gated on the reference and on the oracle bench. The reference-oracle
>    gate is the stronger one anyway: it checks against truth, not against the previous implementation.
>    §4's "≈1,600 written" shrinks accordingly, and F7's legacy-vs-legacy gate is moot.
> 2. ✅ **§3 W7's partition deletion HAPPENED AT W1b-clean**, not at W7. `calibration/regions.py` is
>    gone whole, `regions.feather`/`boundaries.feather` are out of the index, and P3 was replaced by
>    the stronger **I3b** (the signature recomputed from each node's midpoint — no v7 needed, and
>    proven to fire). Ledger W1b-clean.
> 3. ⚠ **F11/Q3's "regenerate goldens twice" did not happen at W1b** — the goldens did not move. All
>    20 golden scenarios have zero mergeable adjacencies, so the partition is unchanged for them.
>    One regeneration remains expected, at W7.
>
> ⭐ **And the finding that reframed W1b:** the 32-condition bench **cannot see a partition change**.
> `ambig_dense_10mb`'s v8 node partition is byte-identical to its v7 region partition (1,698 == 1,698,
> row-for-row). The partition effect was measured on `quick_3to1_5mb` instead (+25.0 % nodes):
> **+0.0069 at pass-0, −0.0014 at refit=3**, and it splits on whether the library contains gDNA at
> all. Ledger W1b.

**2026-07-28.** Written after a full code + implementation review of the v5 blast radius at `3c293038`.
Consumes `05_accumulator_v5.md` (THE SPEC) and `../index/00_splice_graph_design.md` (the v8 graph). This
document does **not** re-derive the design — every closed decision (D1–D8, no path store, no loopy BP, no
fractional mass, Σ1/L not ΣL) is taken as given. It records **what is actually in the tree**, **where the
spec and the tree disagree**, and **the order in which to change it**.

Everything below marked ⭐ MEASURED was run in this session, in the `rigel` env, against the real human
index (`~/Downloads/rigel_runs/refs/rigel_index_v7`) or the shipped code — not quoted from a doc.

---

## 0. THE ONE-PARAGRAPH VERDICT

The rework is a **net deletion of ≈3,800 lines against ≈1,600 written**, and the deletion list is sound: the
fractional-mass model, boundary two-sidedness, and the signature merge are cleanly separable and their
supporting machinery is exactly where §8 says it is. **But the phase plan in v5 §14 does not survive contact
with the code in five places** (§2 below), and **four spec statements are wrong or under-specified** in ways
that block a phase gate (§1). None of these is a design problem; all are sequencing and closure problems, and
all are fixable before a line of production code is written. The largest single change is not
`node_geometry.py` — it is the **48 files that index by `ref_boundary_offsets` / `RegionArrays.from_region_df`**
(65 and 100 references respectively, measured), which the spec's blast-radius table does not count.

---

## 1. FINDINGS THAT AMEND THE SPEC

Each of these changes a number, a gate, or a schema. They need owner sign-off before W1.

### F1 ⭐ MEASURED — the v8 census is **1,043,881 nodes**, not 992,068

The graph doc §3.4 counted with the filter `~is_synthetic & ~is_nrna` (227,844 transcripts). **The shipped
builder filters on `is_synthetic` only** (`calibration/regions.py:179`), which admits **26,475 rows that are
`is_nrna=True` but `is_synthetic=False`**, and their exon spans do contribute events to today's partition.

| | doc §3.4 | ⭐ measured, exact v7 filters |
|---|---|---|
| transcripts in the event set | 227,844 | **254,319** |
| nodes | 992,068 | **1,043,881** |
| median / p25 / p75 | 145 / 59 / 730 | **151 / 61 / 803** |
| 1 bp nodes | 15,315 | **15,687** |
| < 10 bp / < 200 bp | 8.5 % / 58.1 % | **8.2 % / 56.7 %** |
| contiguous edges | ≈ 992 K | **1,043,595** |
| junction edges | 404,168 | **404,168** ✅ exact match |
| strand-coincident junctions | 0 | **0** ✅ (G18 never fires, still must work) |
| ratio to v7 (752,654 regions) | +32 % | **+38.7 %** |

*(The strand filter is inert: all 254,319 non-synthetic rows have strand ∈ {POS, NEG}.)*

⛔ **THE SECOND FILTER BELOW IS STRUCK — IT WAS A BUG (ledger W1c, 2026-07-29).** On a
**non-synthetic** row `is_nrna` does not mean "manufactured span"; it means **single-exon, so mature
≡ nascent**. Measured: all 26,475 `is_nrna & ~is_synthetic` rows have `n_exons == 1` and none is a
`RIGEL_NRNA_*` row, while every manufactured span is already `is_synthetic`. So `~is_synthetic`
alone does exactly what the reasoning below requires, and the extra clause deleted the TSS/TES of
26,475 real transcripts — **52,104 distinct terminus positions**, i.e. the visibility v8 exists to
buy. **ONE filter: `~is_synthetic`.** Invariant **I13** (the flags ARE the events, both directions)
now catches this class.

**Consequence.** ⚠ Using the doc's filter makes **P3 impossible**: v7's interfaces are then *not* a subset of
v8's cuts, so "merge v8's equal-signature neighbours to reproduce v7" cannot hold. The event set **must** use
`~is_synthetic`, exactly as shipped. This closes graph-doc open question #1 by measurement — but note the
**flags are a different question**: `TSS_s`/`TES_s` must come from `~is_synthetic & ~is_nrna`, because an
nRNA span's ends are not real transcript termini. That is precisely P1G_SCOPE §5's "211 false positives are
100 % nRNA-span edges". **Two filters, two purposes; write both down in the builder.**

### F2 ⭐ MEASURED — graph-doc §7.2's **P2 is a broken gate** and would always fail

P2 as written asserts "the cut set of the v8 builder equals the cut set of the v7 builder (`np.union1d` of
v7 region starts/ends)". But the merge *deletes* cuts — that is the entire point of the change. Measured:
**chr1 has 25,552 v8-only cuts, chr2 21,278, chr3 18,661; v7-only = 0 on every one of 286 references.**

Replace with two gates that are both true and both checkable:

* **P2′** — v7's interfaces are a **subset** of v8's cuts. ⭐ Verified true on all 286 references.
* **P3** — merging v8's adjacent equal-signature nodes reproduces `regions.feather` **exactly**. Unchanged,
  and it is the real migration gate.

§3.1's prose ("the cut set is exactly today's cut set") is right about the **event** set and wrong about the
**emitted interfaces**; fix the wording so nobody re-derives P2.

### F3 ⭐ REPRODUCED — D6 is real, **and there is a second un-spec'd instance**

`scratchpad/acc_seam_check.py` re-run this session: shipped pooled-seam density reads **1.994 / 2.002 /
1.981 ×** truth at region lengths 2000 / 500 / 200 (and 1.972 / 1.847 at 100 / 50). The arithmetic:
`gdna_boundary_len = E[min]/2` (`calibrate.py:226`), `seam_m ≈ ρ·(E[min_r]+E[min_{r+1}])/2`, divided by
`0.5·(gbl_r+gbl_{r+1}) = (E[min_r]+E[min_{r+1}])/4` ⇒ **2ρ**. CONFIRMED.

⚠ **v5 §10.4 names only `_pooled_seam_arrays` (capture_eff_length.py:214). The identical `0.5 *` halving
appears a second time, inline and not routed through that function, at capture_eff_length.py:318** (the
splice-junction seam `s_j`). Patching only the named site makes
`test_capture_eff_length.py::test_no_nascent_mature_inversion_under_capture` fail with a real inversion
(nascent 310.17 < mature 412.48). **Both lines must change in one commit.**

### F4 ⛔ §8 deletes `region_types`; §11.4 requires it — RESOLVED: **`region_types` SURVIVES**

§8 deletes "`region_types`/`fl_pool_idx` plumbing through `AccumulatorSet`". §11.4 requires "intergenic
contained fragments are pure gDNA" as a structurally-pure FL pool — and **binning an intergenic-contained
fragment at deposit requires knowing the node's coarse type at deposit.** That *is* the plumbing.

⭐ **OWNER, 2026-07-29 — two domain facts that decide this, and they refute §11.4's premise:**

1. **The FL pools are largely DIAGNOSTIC**, and nascent RNA is extremely sparse in these libraries, so the
   shipped intergenic+intronic pool is a perfectly good gDNA proxy.
2. ⛔ **Under hybrid capture the gDNA FL signal comes from exon↔intron "splash" reads at the boundaries**,
   because those are partially enriched by the capture probes. **Intergenic regions are not captured.**

So an "intergenic-contained only" gDNA pool — §11.4's proposal, and my first plan — would **destroy the gDNA
FL estimate in exactly the cfRNA/capture regime the pipeline is built for**. §11.4's structural-purity
argument is correct about *composition* and wrong about *coverage*.

**RESOLVED: strike "`region_types`/`fl_pool_idx` plumbing" from §8's delete list. Keep the pool concept and
today's pool semantics.** The three improvements that are still real (§5 Q2) are orthogonal to purity:
integer pools, bin-once-per-fragment, and re-cutting the pool definitions on v8's structural flags so
"splash" becomes an explicit class instead of a `type × compartment` proxy.

### F5b ⚠ MEASURED 2026-07-29 — `Σ1/L` is **NOT divisor-free near a terminus** (v5 §3.3 amendment)

v5 §3.3 derives `E[Σ1/L] = ρ_total` exactly, "divisor-free and FL-model-free". That holds only where the
crossing opportunity is **proportional to `w`**. Near a transcript terminus the opportunity saturates at the
reach (`trapezoid(w) = R_a` for every `w > R_a + 1`), so it stops depending on `w`, and the identity fails:

```
    E[Σ1/L]  =  ρ · E_f[ trapezoid(w) / w ]      →  ρ  only when trapezoid(w) ≈ w
```

Simulated on true physics (`scratchpad/reach_worth_it.py`): `Σ1/L ÷ ρ_M` reads **0.992** mid-transcript and
**0.103** at 20 bp from the TES. ⭐ The mechanism: away from a terminus the crossing sample is *length
biased* (measured mean 214.2 vs `E[w²]/E[w] = 212.5`) and `1/L` cancels exactly that; near a terminus the
bias vanishes (mean returns to 200) and `1/L` over-divides. **`Σ1/L` corrects the length bias, `reach`
corrects the placement loss — different corrections, neither substitutes for the other.**

Consequence for the plan: wherever `Σ1/L` is used as `ρ_total` (W5e's moment solve), it needs the same
structural gate as the count. Ledger W1a-Q.

### F5 ⚠ `(count, Σ1/L)` gives **moments, never a pmf** — the FL histogram must survive

§4.1's node-contained and crossing divisors are **integrals against `f(w)`**; §11.4 supplies only an
abundance-weighted **mean** (`flux/flux_recip`). The production consumer is `pipeline.py:862`
`gdna_counts=gdna_fl_mass(calibration_payload)`, plus six more in tests/scripts. With F4 resolved this stops
being an ordering hazard and becomes a simple statement: **a library-wide FL histogram is required and stays.**
§0's "no per-node FL histograms" is untouched — library-wide pools are a different object.

Positive note: the RNA half of §11.4 is *already shipped* — `pipeline.py:861` feeds `rna_counts` from the
scanner's SPLICED_ANNOT histogram, which is already structurally pure.

### F6 ⚠ §8's "delete the `_EPS` floors" is right for EDGES and **wrong for NODES**

At a node the contained divisor genuinely collapses to **exactly 0.0** under v8, for both components, on tens
of thousands of nodes (RNA FL N(200,50): `E_contained` = 0.000 at L=1 and L=10, 0.030 at L=59 = v8's p25;
12.4 % of v8 nodes vs 1.6 % of v7 regions). `fl_mean` at an edge never collapses. **Land the structural-zero
mask first** — an object with `E_c == 0` emits nothing and has zero precision for component `c` — port
`test_node_chain.py::test_zero_opportunity_nodes_carry_zero_mass` to the v5 object model, and only then
remove the floors, and only at edges.

### F7 ⛔ The **W4 gate is unsatisfiable as written**

Two independent reasons:

1. `index.py:1079-1088` is a strict `!=` version gate. Once W1 bumps 7→8, **a v7 index is unloadable by the
   v8 tree**, so no process can produce a v7-partition legacy payload for comparison.
2. The legacy payload's float channels (`boundary_mass_left/right`, `fl_pool_mass`) are fractional doubles
   merged with `+=` over a data-dependent worker partition — the known ~2.6 % cross-process
   nondeterminism. They are **not bit-reproducible at `n_workers > 1` today**, before any v5 change.

**Restate the gate:** *legacy-vs-legacy on ONE FIXED (v8) cut array, across the C++ deposit rewrite only,
exact on the integer channels and either run at `total_threads=1` or tolerance-bounded on the floats.* This
requires a **new legacy baseline recorded at the end of W1b** (32 conditions, v8 partition, both refits),
which becomes the W4 reference. The pre-v5 baseline retires at W1b and must never be used as a W4 comparator.

### F8 ⭐ MEASURED — the coarsening map **cannot reconstruct a v7 effective length**

`region_eff_length(L) = Σ_w f(w)·max(0, L−w+1)` is superadditive. Over the **118,195 v7 regions that actually
split** under v8, `Σ E(children) / E(whole)` = **0.7652**; for a 305 bp region splitting into 145 + 160 the
ratio is **0.0917**.

**W0b must be scoped honestly.** The map **can** (i) define comparable populations and (ii) pool the
**oracle truth** — per-object gDNA/RNA counts are additive, so `f_o = ΣG/(ΣG+ΣR)` is well defined on any
object set. It **cannot** reconstruct a v7 contained density, a v7 contained eff-length, or a v7
`mass_global`. Any plan step that assumes otherwise is wrong.

### F9 ⛔ C1, C2 and C3 are **welded through one 10-line function** — W2 cannot be one arm

`_rho_faces` (bp_solver.py:537-546) consumes `accept_l`/`accept_r` (C3's target, :383-386) and is the sole
producer of `rho_l0`/`rho_r0` (:549), which feed `_seam_pair` :1042 (C2's fit population), `_relay`
:1046-1047 and `_transport` :1329-1334 (the reframe C1 gates), and `_flank_dom` :993 (C2's fit).

**W2 becomes three sub-arms with a re-recorded baseline between each: C3 → C1 → C2** — in that order, because
C3 changes both C1's and C2's inputs while neither changes C3's. ⚠ **This inverts P1G_SCOPE §10's recommended
order**, which puts C3 last as "lowest priority". It is lowest *value* and first *in sequence*.

### F10 ⚠ C3's specified predicate is **inverted** relative to what `accept_*` means

`accept_l/accept_r` means *"this face carries the one-sided mature (spliced) density"*. The specified
replacement `rna_crosses_contiguously(p,s) = exon_s(src) AND exon_s(dst) AND NOT (TSS_s or TES_s)` is
*"mature crosses this seam without splicing"* — nearly the **complement** at a junction. ⭐ Measured: only
**13.7 %** of the 404,168 junction seams have `exon_s` on both flanks, so a literal substitution would drop
the mature term from `ρ_tot` at ~86 % of junction seams and **fail P1G_SCOPE's own falsification test**
("junction-only edges must be UNMOVED"). The correct structural analogue is **per-face**:
`accept_face(s) = DONOR_s or ACCEPTOR_s at this seam, on the flank carrying exon_s`. C3 is not a 4-line
substitution; it also needs `node_total_density` to return the two per-face spliced densities separately
(today it sums them, node_geometry.py:370-372 — safe only because `accept_*` is currently observational).

### F11 ⚠ Goldens move at **W1b**, not W7

`test_golden_output.py` builds its own index in-test and runs the full pipeline including calibration,
comparing `_LOCI_NUMERIC_COLS` (which includes `gdna_prior_count`, `rna_prior_count`, `gdna_eff_len_em`) at
`rtol=1e-6`. A finer partition moves all three. **The standing "full suite" gate is unachievable from W1b
onward.** Choose explicitly: regenerate goldens **twice** (W1b and W7, both recorded in the ledger), or scope
the standing gate to `--deselect tests/test_golden_output.py` for W1b..W6. Silently carrying 21 red goldens
through six phases is the failure mode to avoid.

### F12 ⚠ **There is no baseline on disk**, and two live docs quote different ones

`/tmp/pass0_oracle_bench.tsv`, `/tmp/rigel_ablate` and `accumulator_ledger.md` do not exist. `P1G_SCOPE.md:239`
quotes "r0 0.078786 / **r1** 0.052470"; `ROADMAP.md` and `HANDOFF_18` quote "refit=0 0.079005 / **refit=3**
0.046675". P1G_SCOPE's is stale (HEAD runs refit=3). Until W0.1 lands, **no arm has a valid comparator**.

### F13 ⚠ `Accumulator::merge_from` is **not bound to Python** — test A9 has zero infrastructure

Defined at accumulator.cpp:234, called only from bam_scanner.cpp:1358, absent from the nanobind
`nb::class_<Accumulator>` block. A9 ("N workers in any order → bit-identical") is one of v5's two starred
determinism claims and **cannot be written today**. Bind it in W0.7 and write A9's harness against *today's*
accumulator first, so the claim has a control.

### F14 ⭐ DERIVED — the mature taper lives **only at edges**, and it is **one function**

The spec leaves `E_M`'s closed form open (§14 W0 calls it "the single remaining item that changes the
schema"). It resolves cleanly:

* **At a node, the taper is vacuous.** Nodes are homogeneous under v8, so for a node inside a transcript's
  exon set, `reach_M(s) = (node_end − s) + (exonic bases beyond node_end) ≥ node_end − s` — the
  node-containment constraint always binds first. **Node-contained mature eff-length is exactly
  `node_contained_eff_length(B, f_r)`, untapered, gated to 0 where no mature covers the node.** §4.1's
  "over `B ∩ exons`" is therefore vacuous: it is `B` or `∅`.
* **At an edge the taper is the whole story, and a contiguous edge and a junction edge take the SAME
  formula.** A mature fragment crossing genomic point `p` needs `a ≥ 1` exonic bases on the donor side and
  `w−a ≥ 1` on the acceptor side, both fitting inside the transcript — which is §4.2b's trapezoid
  `E_f[min(w−1, R_d, R_a, R_d+R_a−w+1)_+]` with `R_d`/`R_a` the exonic reaches either side of `p`. This is
  exactly §4.2's own claim ("it applies to both object classes, and it is the *same* quantity").
  `fl_mean` is its large-reach limit (§4.2b's table: 199.0 at R = 5000 *and* at R = 550), **not a separate
  case**.

**So there is ONE new divisor function, `crossing_eff_length(fl_pmf, R_d, R_a)`, used at both edge classes**,
and §4.1's "contiguous edge | mature | `fl_mean`" is its asymptote. That is the elegant reading and it
removes a special case rather than adding one. It requires **reaches on contiguous edges too**, not only
junction edges — a schema consequence the spec does not state.

⭐ **DECIDED BY MEASUREMENT, 2026-07-29.** The owner's objection was that a rich annotation carries many
*potential* TSS/TES, so a reach computed against it is speculative — and that MAXIMAL reach (the furthest
terminus on that strand) is the fair upper bound. **That is exactly D2**, already closed in §4.4 by
measurement (75.7 % of exonic positions have zero isoform disagreement; maximal vs the isoform mean is worth
3.29 % of the mature opportunity, against the taper's own 10.7 %). ⚠ Its bias direction, stated plainly:
maximal reach assumes the longest isoform exists ⇒ **over**-states opportunity ⇒ **under**-states `ρ_mature`
⇒ **over**-states `f_g`. One-sided, bounded, conservative toward gDNA.

The open sub-question — *does maximal reach collapse the taper to `fl_mean` at ordinary positions, making a
contiguous-edge reach column dead weight?* — is now measured. Both reaches taken maximal **independently per
side**, 3,000 sampled human genes, 10,957,350 exonic positions, RNA FL N(200,50):

| | |
|---|---|
| bp-weighted mean `E_cross / fl_mean` | **0.8904** ⇒ **11.0 %** systematic gDNA over-call if ignored |
| positions differing from `fl_mean` by > 1 % / > 5 % | **25.9 % / 21.6 %** |
| positions differing by > 25 % / > 50 % | **15.4 % / 10.2 %** |
| p1 / p5 / p10 / p25 of the ratio | **0.045 / 0.240 / 0.488 / 0.986** |

This independently reproduces §4.3's `0.893` / 10.7 % — and establishes that **§4.3's figure was already
computed under maximal reach**, so the 10.7 % is not an artefact of an optimistic convention. **The taper
does NOT collapse: reaches go on both edge classes.**

### F15 ⭐ MEASURED — the blast radius is understated

| symbol | spec / survey says | ⭐ measured |
|---|---|---|
| `ref_boundary_offsets` | — | **30 files, 65 references** |
| `RegionArrays.from_region_df` | 44 call sites | **100 call sites** |
| `build_node_chain` (scripts/debug) | 10 | **17** |
| the spliced selector | "~40 lines" | **74 lines** (node_geometry.py:189-262) |

Budget a compatibility shim that synthesises `ref_boundary_offsets` from `ref_node_offsets`, and keep
`from_region_df` as a deprecated alias through W5. ⚠ `scratchpad/` is **not** in the ruff target set, so its
~40 gate scripts will not surface mechanically — enumerate them by hand.

### F16 — smaller confirmed items

* ⭐ **The merge is 6 lines** (regions.py:128-133) and its invariant is 8 (regions.py:515-522). Deleting those
  6 lines and nothing else converts today's partition into v8's node set.
* **`index.boundary_df` is 100 % dead at runtime** — assigned at index.py:1186, consumed at :1187 by a
  validator whose only job is proving boundaries still equal the region interfaces. ⚠ But it is the *only*
  load-time cross-check catching a stale partition; do not delete it before `validate_graph` is wired at load.
* **`pipeline.py:281-287` silently skips calibration** for an index without `region_df`. A v8 index would
  disable calibration with no error anywhere. One line, fix it first.
* **`CalibrationResult._check_region_array` hard-rejects any non-float64 array** (result.py:22-35) — it
  blocks every integer-count arm and is on no delete list.
* **Three neighbour-differs sites, not one**: `test_regions.py:64`, `test_region_index_alignment.py:31-37`,
  and the validator itself.
* **`node_chain.build_node_chain` hard-asserts `(b1−b0) == k+1`** with the message "rebuild the index" — it
  fires first on a v8 partition and points the implementer at the wrong thing.
* **`fragment_genomic_spans` survives verbatim, and more strongly than the spec claims**: it fills the mate
  gap and cuts only at cut introns, so **every inter-span gap it produces IS a spliced intron** ⇒ v5 §7's
  junction lookup needs no ABI change. ⚠ One edge case to instrument, not comment: bam_scanner.cpp:776-785
  clips a cut intron to the fragment extent, and a leading intron starting exactly at `lo` emits no span
  before it — if that fires, the lookup gets clipped coordinates and silently misroutes an annotated
  junction to §5.1's unannotated path.
* ⭐ **The FL pmfs are invariant under the v8 partition** — all 118,195 splitting v7 regions are
  `coarse_type == EXON`, and `_GDNA_POOLS` is intergenic+intronic only. So `gdna_fl_pmf`/`rna_fl_pmf` do not
  move at W1, and **W1b isolates the geometry change cleanly**. Rely on this.
* ⭐ **The builder is fast**: a vectorised cut-set + junction pass over the real human GTF runs in **0.53 s**
  against the < 10 s budget, peak RSS 0.98 GB (dominated by loading the feathers).

---

## 2. TARGET ARCHITECTURE

The standard is fewer modules, fewer constants, one concept = one name. The rework should *end* with a
smaller calibration package than it started with.

```
DIES                                          BORN
────────────────────────────────────────────  ──────────────────────────────────────────────
calibration/regions.py         (548)          calibration/splice_graph.py        (~380)
calibration/region_arrays.py   (192)            build_splice_graph, validate_graph, SpliceGraph
calibration/node_chain.py      (112)            (nodes+edges arrays + CSR + the node↔edge chain)
                               ─────
                               852 lines       one module, one concept: THE GRAPH

node_geometry.py: the 18 per-face arrays      calibration/object_geometry.py     (~260)
  and the 74-line spliced selector              ObjectGeometry: per-object (count, recip, E_g, E_r)
  (node_geometry.py keeps NodeBelief /          MatureIncidence: per node, its junction edges' (flux, E_J)
   NodeStatics / init_beliefs, ~380 lines)

effective_length.py's 3 fractional divisors   effective_length.py, rewritten     (~140)
  (87 lines + 12 of docstring)                  node_contained_eff_length  (= today's, byte-identical)
                                                node_contained_recip_eff_length      NEW
                                                crossing_eff_length(fl_pmf, R_d, R_a) NEW — F14, ONE fn
                                                fl_mean                    (promoted, its asymptote)

the 3x2 FL pool GRID (reworked, not deleted)  2-3 named pools, INTEGER, binned once per fragment
  fl_pool_idx(region_type, boundary)            gdna  = intergenic + intronic (unchanged proxy)
                                                splash = crosses a DONOR_s/ACCEPTOR_s edge  <- named, not proxied
                                                mature = crosses a junction edge (free cross-check)
                                                region_types SURVIVES (F4 -- capture needs splash)

substrate.py's BoundarySubstrate + the        substrate.py: NodeSubstrate + EdgeSubstrate
  (counts, mass) duality (~90 of 203)            one ObjectView(count[c], recip[c])

C++ struct Boundary + slices + n_cross        C++ struct Node{u32 count[2]; u64 recip[2]}
  + FL grid + junction_strand (~170)                 struct Edge{u32 flux[c]; u64 flux_recip[c]}
```

**Two structural points the spec leaves implicit, stated here as the plan's position:**

1. **Junction edges are FACTORS, not BP variables.** Every junction closes a cycle (graph doc G7), and
   "no loopy BP" is settled. So the linear chain is `node ↔ contiguous-edge ↔ node ↔ …` — the *same* shape as
   today's bipartite chain, minus the two per-reference terminals (v8 has `k−1` contiguous edges for `k`
   nodes, not `k+1` boundaries). A junction edge's `(flux, Σ1/L, E_J)` enters as a local measurement at its
   incident nodes, which is what "graft/peel move to junction edges" means. **`build_node_chain` shrinks
   rather than grows.**
2. **`AccumulatorPayload` becomes two classes during shadow mode, not one rewritten class.** §9's table reads
   as a rewrite; shadow mode needs both live simultaneously. `scan_payload_v5.py` alongside, merged at W7.

---

## 3. THE PHASE PLAN

Changes from v5 §14 are marked ⚠. Every arm: ruff · full suite (see F11) · one thing varied · a
pre-registered falsification test · baseline re-recorded from the current tree in the same session · a row
appended to `docs/accumulator/accumulator_ledger.md`.

### ⚠ W0 — INSTRUMENTS, AND THE TWO INDEPENDENT v7 ARMS *(new; none of it touches the graph)*

| | task | gate |
|---|---|---|
| W0.1 | Re-record the baseline from a `git stash` of HEAD, **both refits (0 and 3)**. Create `accumulator_ledger.md`. Fix P1G_SCOPE §8's stale r1 quote (F12). | 32/32 vs HEAD |
| W0.2 | Unify the `z2` denominator: factor `_lin_var` (pass0_error_table.py:132) into one shared helper; fix `ablate_report.py:73` and `v_p4_paired_z2.py:28-32`, which use a different denominator. Re-record held-fixed z2 on the fixed scale. | one definition, one number |
| W0.3 | Bench provenance: add `partition` (v7/v8) and `mass_kind` columns to `pass0_oracle_bench`'s TSV; make `--report` and `p2res_report.py` **refuse** to diff arms across partitions without `--coarsen`. | attempted cross-diff errors |
| W0.4 | `index_hash` into `manifest.json`; namespace `/tmp/rigel_selfsolve` **and** `calib_cache` by it; `calib_cache.load()` raises on mismatch (it stores `index_dir` today and never checks it). | stale cache raises |
| W0.5 | ⚠ **ARM — D6, moved here from W5f.** Pure v7 arithmetic, independent of index/accumulator/solver. **Both** halving sites (capture_eff_length.py:214 **and** :318, F3) in one commit, plus the three fixture builders and the ~10 tests. | pre-register a **2× seam-density move**; promote `acc_seam_check.py` to a unit test; add the falsification case where the clip does **not** rescue (true seam density strictly inside `(ρ_ref/2, ρ_ref)`) |
| W0.6 | ⛔ **CANCELLED — claim REFUTED by its own gate (2026-07-29).** The gDNA seed weight is **not** identically 1.0: 23–66 % of seeds on real cfRNA and up to 100 % on synthetic differ from it (ledger W0.6). The channel is LIVE; no deletion, and `boundary_side_eff_length` keeps a real non-geometry consumer until W5. | gate ran, failed, arm dropped — `src/` untouched |
| W0.7 | Bind `Accumulator::merge_from` through nanobind; write A9's determinism harness against **today's** accumulator (F13). | A9 harness green on HEAD |

*Rationale for moving D6 and the seed-weight arm to W0: both are independent of everything else, and landing
them first removes two confounders from every subsequent arm's attribution.*

### ⚠ W0b — THE COARSENING MAP, HONESTLY SCOPED

Build `CoarseningMap` from cut arrays alone (two `searchsorted` + one `isin`), with four asserted invariants
(v7 interfaces ⊆ v8 cuts; every v8 node maps to exactly one v7 region; the map is surjective; per-reference
counts agree).

⚠ **Scope it to what F8 permits:** it defines comparable **populations** and pools the **oracle truth**
(counts are additive ⇒ `f_o = ΣG/(ΣG+ΣR)` is well defined on any object set). It does **not** reconstruct v7
eff-lengths, densities, or variances. **Prove it on a synthetic refinement before v8 exists** — split every
v7 region at its midpoint, re-derive the sub-node oracle by proportional allocation, coarsen back, assert
recovery.

✅ **DONE 2026-07-29 — `scripts/debug/coarsen.py` (ledger W0b).** Validated on the REAL v8→v7 relationship
(1,043,881 → 752,654, all four invariants, 0.7 s at human scale) as well as the synthetic refinement.
Two measured limits now on record: **27.9 % of v8 boundary slots have no v7 counterpart** and are excluded
(reported, never silent), and the **Jensen gap is 0.739×** — a coarse mwae reads 26 % better than the fine
one with no solver change at all, which is why Q4 scores on the fine population.

⚠ **Then decide the standing gate explicitly** (F8 makes the v5 §14 wording undefined). Recommendation:
score `z2` on the **fine (v8) population** against a baseline freshly recorded at W1b, and use the coarsening
map only for the *cross-partition* mwae sanity check.

### ⚠ W1 — INDEX v8, split into three arms

| | arm | gate |
|---|---|---|
| **W1a** | Build + validate + write `nodes.feather` / `edges.feather`; `INDEX_FORMAT_VERSION` 8; **not wired to the scanner**. Includes the flags (two filters, F1), the junction edges, and `reach_donor`/`reach_acceptor` on **both** edge kinds, maximal per side (Q1). | ⭐ **bit-identical 32/32** — this is the real wiring gate, and it *is* achievable because nothing reads the new artifact yet. G1–G18, I1–I12, P1/**P2′**/P3/P4/P5, determinism, §8 budgets. Plus: `crossing_eff_length` reproduces §4.2b's table (199.0/199.0/198.2/160.1/87.8/19.6/50.0) and the F14 census (0.8904 bp-weighted). |
| **W1b** | Flip `build_region_partition_arrays` to the node cut array. | ⚠ **expected to move.** Record the delta in the ledger as **the PARTITION effect**. ⭐ Record the new legacy baseline here — it becomes W4's reference (F7) and the fine-population `z2` baseline (Q4). **Regenerate goldens (1 of 2)** (F11/Q3). |
| **W1c** | The flags → `calibrate()` plumbing that W2 needs: `index.edges_df`, `build_boundary_flags_array(index) -> uint16[B]` aligned to the payload offsets (assert the +1 offset and the two zeroed terminals), a new keyword on `calibrate()`, new bool fields on `NodeStatics`. | bit-identical (nothing reads them yet) |

⚠ **This resolves P1G_SCOPE §8 item 2**, which demands the index arm be bit-identical 32/32 — true for W1a,
impossible for W1b, and the two were conflated because P1G_SCOPE was written for restoring columns on the
*v7* partition.

Also in W1: fix `pipeline.py:281` to an explicit v7/v8 dispatch that **raises** (F16); relax
`CalibrationResult._check_region_array`'s dtype gate; rebuild every index and re-key every cache.

### ⚠ W2 — P1g OFF THE INDEX ALONE — three arms, order **inverted** vs P1G_SCOPE §10

Baseline re-recorded between each (F9).

* **W2a — C3.** ⚠ Not a 4-line substitution (F10). Use the **per-face** structural predicate
  (`DONOR_s`/`ACCEPTOR_s` on the flank carrying `exon_s`), and split `node_total_density`'s summed spliced
  density into its two faces. Pre-register the exon↔exon-junction case (donor on one flank, acceptor on the
  other) explicitly.
* **W2b — C1.** The terminus gate on the **gDNA leg only** of the reframe. ⛔ Do not divide `r` by `f_g(dst)` —
  that is the pin bug rebuilt. Pre-registration per P1G_SCOPE §8, with its falsification test
  (**junction-only edges must be UNMOVED**).
* **W2c — C2.** `ω_graft` per structural class. ⚠ Watch the fit-population bias: exons with one live seam
  never enter the fit and are the worse-behaved half, and terminus-flanked exons are disproportionately
  one-seam. Decide the minimum per-class pair count before falling back to the pooled scalar.

### W3 — PYTHON REFERENCE ACCUMULATOR + TESTS, written to target and FAILING

A1–A14, E1–E4, S1–S2. Lift `acc_proto_e.py` whole as S1 and the A3/A4 fixture (its
`A=500 B=50 C=10 D=1 E=1000` geometry already covers the 1 bp node and the encompassed middle);
`acc_proto_g.py` as S2; `acc_proto_d.py`'s brute-force enumerators as E1–E4's oracle. Add a multi-reference
`build_test_index` and a declarative toy-graph fixture. Run the reference on a **real cfRNA payload** through
a shim before any C++ exists. Build `pass0_oracle_bench`'s v5 object mapping here.

⚠ Three contract items must be **frozen in the header before W4**, not discovered: the exact rounding of
`R = round(2^32/L)` (`llround` vs banker's — they differ at ties and byte-identity is undefined without it);
whether `L` is pre- or post-clip (observable in A14); and the per-class channel API.

### W4 — C++ ACCUMULATOR + SHADOW MODE

`set_graph(...)` as a **second** method, not an overload (the one-shot guard at bam_scanner.cpp:1118-1122
throws today). ⚠ **Budget peak RSS in W4, not W6**: ~90 B/node × 1.04 M v8 nodes ≈ 94 MB/copy × 9 at 8
workers ≈ 845 MB for the legacy arm alone, plus v5's 81 MB/copy ⇒ **~1.6 GB peak in shadow mode**.

⚠ **Gate restated per F7**: byte-for-byte vs the Python reference; legacy-vs-legacy on one fixed v8 cut array
across the deposit rewrite, exact on integer channels, at `total_threads=1` (or tolerance-bounded on floats,
stated).

⚠ **The FL pools are reworked here, not deleted** (Q2): integer, binned once per fragment, with `splash`
named structurally off the `DONOR_s`/`ACCEPTOR_s` flags. Its own arm, with a recorded `gdna_pmf` mean delta —
this is the arm that must not silently degrade the capture-regime gDNA FL, so gate it on the four cached
cfRNA payloads, never the toy.

### ⚠ W5 — CONSUMERS, one observable per arm *(f is gone — D6 moved to W0.5)*

`W5a` contained counts · `W5b` contiguous-edge flux · `W5c` junction edges + `E_J` · `W5d` the §4 effective
lengths (⚠ structural-zero masks **before** the `_EPS` removal, F6) · `W5e` the `(count, Σ1/L)` moment into
`node_init` — gated on **real cfRNA**, where the honest failure case of §11.2 must reproduce (an `L=100` node
reports high `ρ_g` precision and near-zero `ρ_r` precision, not a confident `f_g`).

⚠ Two hazards to budget: `chain_boundary_side_deconv` must **not** be deleted in the same arm as anything
touching the seam model; and the `_relay`/`_transport` and `_peel_share`/`_peel_share_scalar` twin pairs are
deliberate, explicitly forbidden from merging (a measured 15.7×/op), and the face index threads through six
functions — **every face-removal edit must be hand-mirrored into both arms with nothing enforcing it.** This
is the single highest-risk mechanical change; re-measure on `cfrna:LBX0190` after each edit, never the toy.

### W6 — OPTIMIZATION (budgeted)

Profile on `cfrna:LBX0190`. ⚠ Note the object count: 1.04 M nodes + 1.04 M contiguous edges + 404 K junction
edges ≈ **2.49 M objects vs today's 1.5 M chain nodes** — v5 §14's "may get FASTER" holds only if junction
edges are factors, not ψ solves (§2 point 1). Confirm which objects are solved before assuming the win.

### W7 — DELETIONS, THEN GOLDENS (2 of 2)

Run `ruff check src/ tests/ scripts/` after the API change and treat undefined-name failures as the
authoritative delete list. ⚠ `scratchpad/` is not linted — enumerate its ~40 gate scripts by hand.

---

## 4. DELETION CENSUS (verified, with locations)

| area | lines | key items |
|---|---|---|
| C++ accumulator | ~110 | `struct Boundary` (h:70-76) · `struct Slice` + expansion (cpp:19-28, 102-126) · the 50/50 crossing block (cpp:181-220) · `inv_L` · `boundary_junction_strand` (5 sites) |
| C++ bam_scanner | ~110 | the four boundary matrices' alloc/copy/dict (1934-2006) · 5 nanobind strided views (2660-2744) |
| `effective_length.py` | 99 | `boundary_side_eff_length` (32) · `boundary_side_crossing_count_eff_length` (23, **zero live consumers**) · `spliced_side_eff_length` (32) · `boundary_eff_length` alias (3) · 12 docstring |
| ~~C++ FL pools~~ | ~~130~~ | ⛔ **struck — `region_types` and the pools SURVIVE (F4).** Reworked in place: integer, bin-once, splash pool. |
| `node_geometry.py` | ~150 of 644 | the spliced selector (189-262, **74**) · `_spliced_faces` + exon masks (11) · 6 mass fields + 6 `_right` twins · `node_global_geometry` (19) |
| `substrate.py` | ~90 of 203 | `BoundarySubstrate` (60) · `left`/`right` + the D1 projection (12) · `_make_view`'s mass half |
| `region_arrays.py` | ~79 | `region_boundary_indices` (26) · `boundary_region_indices` (30) · `RegionArrays.order` (**dead field, zero readers**) |
| `regions.py` | ~150 of 548 | the merge (**6**) · neighbour-differs validator (8) · `_build_reference_rows` (57) · the entire boundary half (~87) |
| `index.py` | ~120 | `boundary_df` load+validate (10) · the write (3) · constants · prose |
| `bp_solver.py` | ~50 | `accept_l/r` (4) · `chain_boundary_side_deconv` (31) · `MS`/`_face_spl` (9, **already dead**) |
| `capture_eff_length.py` + `priors.py` | ~85 | `_pooled_seam_arrays` (24) · the junction-seam imputation (16) · `rekey_right` (24) ⚠ *its spirit must survive — see below* |
| `density_model.py` + `run_fill.py` | ~145 | `node_gdna_density` (87) + docstring (35) · `runfill_bidirectional` (27, **one consumer, dies with it**) |
| tests | ~900 | 17 of 34 in `test_accumulator_spec.py` (~250) · `test_density_model.py` (174, whole) · `test_message_frames.py` (151, whole) · `test_spliced_boundary_onesidedness.py` (159, whole) |
| `scripts/debug` | ~1,695 | the 13-file boundary-two-sidedness diagnostic family |
| **total** | **≈ 3,800** | against ≈ 1,600 written (the FL-pool rework is a rewrite, not a deletion) |

### ⛔ DO NOT DELETE (each was proposed and is wrong)

* **`region_eff_length` or its `+1`.** It is bit-for-bit §4.1's node-contained first moment and carries a
  hard-won discrete-count fix; without it, `eff = 0` was a live division by zero on 12.4 % of v7 regions —
  and v8 makes that population much larger. Rename to `node_contained_eff_length`, normalise the argument
  order, assert `np.array_equal` at the rename. Do not "simplify" the `+1`.
* **`fl_mean`.** §8 says nothing about it; it is *promoted* to the divisor for every contiguous edge.
* **`test_substrate_conservation.py::test_contained_mass_equals_count`.** It is an exact integer identity —
  §3.4's identity #1 in disguise — not a fractional-mass test. Rewrite it, don't drop it.
* **`derive.gdna_density_global`'s divisor choice.** It divides by `gdna_boundary_len` un-halved and
  un-averaged, which is *correct*, and it is the independent witness that `_pooled_seam_arrays` — not
  `gdna_boundary_len` — is the defective party in D6. Never harmonise it toward the seam form.
* **`strand_likelihood.py`, `strand_balance.py`, `strand_summary.py`.** They touch neither payload nor
  substrate. Any diff in them during W5 is a mistake — use them as a cheap tripwire.
* **`gdna_strand.fit_rna_strand_from_sj_table`.** Do not "unify" it with the junction edges' 2 channels: the
  accumulator's spliced channel pools SPLICED_IMPLICIT (arbitrary sense bit ≈ 0.49) and fitting `od_r` from
  it saturated 4/4 real libraries. The v5 junction edge may carry the same contamination — measure first.
* **`priors._gdna_region_node_arrays`'s `rekey_right` spirit.** The heuristic dies with
  `_pooled_seam_arrays`, but the failure it guards (a locus's far-left outer seam silently dropped,
  under-counting its gDNA) returns unless an edge whose `src` node is intergenic gets an explicit locus
  owner. It was found by reasoning, not by a test — nothing will catch its return.

---

## 5. DECISIONS — RESOLVED 2026-07-29

* ✅ **Q1 — reaches on BOTH edge classes.** Convention: **MAXIMAL**, taken independently per side (D2,
  owner-confirmed). One function `crossing_eff_length(fl_pmf, R_d, R_a)`; `fl_mean` is its large-reach limit,
  not a branch. Settled by the F14 measurement, not by preference: the taper survives maximal reach
  (`E/fl_mean` = 0.8904 bp-weighted; 21.6 % of exonic positions off by > 5 %). Cost: `reach_donor`/
  `reach_acceptor` int32 on contiguous edges too, ≈ 8 MB on `edges.feather`.
  ⚠ Record the bias direction in the ledger: maximal reach over-states opportunity ⇒ over-states `f_g`.
* ✅ **Q2 — `region_types` and the FL pools SURVIVE; §8's delete entry is struck.** Owner: the pools are
  largely diagnostic, nascent is sparse enough that intergenic+intronic is a good gDNA proxy, and — decisive
  — **under hybrid capture the gDNA FL signal comes from exon↔intron "splash" reads, not from intergenic
  regions, which are not captured.** §11.4's "intergenic contained ⇒ pure gDNA" pool is right about
  composition and wrong about coverage. Keep today's semantics; take the three orthogonal improvements below.
* ✅ **Q3 — regenerate goldens twice** (W1b and W7), both recorded in the ledger. The full suite stays a real
  gate throughout.
* ✅ **Q4 — score `z2` on the fine (v8) population** against a baseline recorded at W1b. The coarsening map is
  used only for the cross-partition mwae sanity check, which is all F8 permits it to do.
* ✅ **Q5 — D6 and the seed-weight arm move to W0.** Both are pure-v7 and independent; landing them first
  removes two confounders from every later arm. v5 D6 says "fixed inside this work" — this moves it *earlier*
  inside this work, not out of it.

### Q2's three improvements (kept; none of them touches pool purity)

1. ⭐ **Integer pools, binned once per fragment.** Today a crossing fragment smears its FL *fractionally*
   across every region-type it touches — `fl_pool_mass_[...] += (sl.end − sl.start) · inv_L`
   (accumulator.cpp:224-231) — accumulated as `double` and merged with a non-associative `+=`
   (accumulator.cpp:263-268). That is one of the known contributors to the ~2.6 % cross-process
   nondeterminism. Under v5 the natural form is to bin each fragment **once**, as an integer. Simpler, and
   bit-exact at any thread count — the same property §6 claims everywhere else.
2. ⭐ **Make "splash" an explicit pool.** The owner named the physically meaningful class: the exon↔intron
   seam. Today it is approximated by `region_type ∈ {intron, exon} × compartment = boundary`. With v8's
   structural flags it can be named exactly — a contiguous edge carrying `DONOR_s` or `ACCEPTOR_s` **is** an
   exon↔intron seam. Strictly better than the proxy, for exactly the measurement the owner relies on under
   capture.
3. **A junction-edge FL pool, free.** Junction-crossing fragments are pure mature RNA by construction, so
   they give an independent RNA FL estimate to cross-check the SPLICED_ANNOT histogram. The edge object
   already exists; this is diagnostic value at no new machinery.
   ⚠ Not an improvement, just a flag: `POOL_EB_PRIOR_ESS = 1000.0` (fl.py:60) carries its own "revisit the
   default on real-data pool sizes" comment. Left alone — raising it is a constant change and needs its own
   discussion.

---

## 6. OWNER DOMAIN FACTS THAT CONSTRAIN THE WORK

Recorded here because neither the spec nor the code states them, and both would have been violated:

1. ⛔ **Under hybrid capture the gDNA fragment-length signal comes from exon↔intron "splash" reads at the
   boundaries** — they are partially enriched by the capture probes. **Intergenic regions are not captured.**
   ⇒ any "pure gDNA = intergenic" pool silently guts the FL fit in the cfRNA regime. Refutes §11.4's pool.
2. **Nascent RNA is extremely sparse in these libraries**, which is why the shipped intergenic+intronic pool
   is a sound gDNA proxy despite not being composition-pure. ⇒ do not "fix" the proxy without measuring.
3. **The FL pools are largely diagnostic.** They are not load-bearing for the deconvolution; they set
   divisors and are read as QC. ⇒ their arm is gated on not-degrading, not on improving.
4. **Which annotated TSS/TES are actually used is unknowable before the EM runs** — the annotation carries
   many merely-potential termini. ⇒ the reach convention must be prior-free and one-sided, which is what
   MAXIMAL is (D2). A `θ_t`-weighted reach would couple calibration to quantification and break the one-way
   pipeline. *(There is a cheap, prior-free upgrade if maximal turns out to bite: §4.4(e) — use the
   MANE/canonical transcript's reach where one exists, falling back to maximal. `transcripts.feather`
   already carries `is_mane`/`is_ccds`/`is_basic`, and `acc_isoform_reach.py` can A/B it directly. Not
   proposed now; recorded so it is not rediscovered.)*

---

## 7. WHAT I AM CONFIDENT ABOUT

The design is sound and the code is more separable than the spec claims: the merge is 6 lines, the fractional
model is ~110 C++ lines, `fragment_genomic_spans` genuinely survives untouched (and its inter-span gaps are
*always* spliced introns, so the junction lookup needs no ABI change), the FL pmfs do not move at W1, and the
builder runs in half a second. The risk is not in the physics — it is in **attribution**: four entangled
effects, a baseline that does not exist on disk, a coarsening map that cannot do what the phase plan assumes,
and twenty-one goldens that go red six phases early. W0 exists to fix all four before anything moves.
