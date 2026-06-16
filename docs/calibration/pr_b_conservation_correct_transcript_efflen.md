# PR-B: conservation-correct mature/nascent transcript effective length

**Status:** plan (for the author to implement next).
**Target:** `src/rigel/calibration/capture_eff_length.py` — `transcript_capture_eff_lengths` + `_exon_region_incidence`.
**Reference to mirror (shipped):** `src/rigel/calibration/priors.py::assemble_priors` (PR-A, lines 215–326).
**Predecessor:** PR-A migrated the **gDNA component** locus eff-len to the pooled-seam node model. PR-B migrates the **other consumer** — the per-transcript MATURE/NASCENT EM eff-len — to the same model.

---

## 1. GOAL & CONTEXT

### 1.1 What PR-B does

Rebuild `transcript_capture_eff_lengths` (`capture_eff_length.py:99-164`) on the **conservation-correct pooled-seam node model** that PR-A shipped for the gDNA component. The output is unchanged in role — the per-transcript EM effective length `effective_lengths_em`, consumed only by the EM (`pipeline.py:446-452`; the TPM-facing full FL-marginal `effective_lengths` is untouched). What changes is the **per-region gDNA accounting** that builds the contraction factor.

### 1.2 Why — the other consumer still violates conservation

The current target is on the **retired path** that PR-A replaced for the gDNA side:

- It calls `_transport_boundary_flux` (`capture_eff_length.py:31, 127-134`) to fold each boundary's crossing mass into adjacent regions and iteratively move mass — the OLD transport path.
- It divides per-region gDNA mass by `gdna_geom_len` (`capture_eff_length.py:135-136`) and sums `gdna_geom_len` into the transcript span (`:147`). `gdna_geom_len = region_eff_len + left_eff + right_eff` (`derive.py:55`) — region length (FL-subtracted) PLUS both FL-scale boundary sides. PR-1 measured this sums to ≈ `1.065×` the genomic span (≈2× on small exons): a **genomic-length-conservation violation**. The EM divisor over-counts genomic length on small exons.

PR-B replaces both with the pooled-seam model: region nodes carry **contained mass only** at **genomic** length (`region_size_bp`), plus **one pooled seam node** per internal boundary at **FL-scale** length contributing **0** to the genomic span. This makes the per-transcript geometric span exact (§4) and makes mature/nascent diverge through the node set (not the length divisor).

### 1.3 What PR-A established that PR-B inherits

PR-A is the reference implementation. PR-B reuses three established mechanisms verbatim (mirror `priors.py:215-326`):

1. **Pooled seam nodes (NOT split).** Per internal boundary, ONE node with mass `mass_gdna_right[r] + mass_gdna_left[r+1]` (both halves POOLED), length `min(gdna_boundary_len[r], gdna_boundary_len[r+1])` (FL-scale), keyed to left-flank region `r`, gated on same-reference adjacency (`ref_id[r]==ref_id[r+1]`). Pooling not splitting is load-bearing: by Cauchy-Schwarz, splitting the seam mass into two nodes halves the IPR support → halves the contraction → nearly doubles the leak (PR-A measured split 13.5% vs pooled 5.5%). See `priors.py:222-230, 239-244`.
2. **The Laplace-smoothed IPR closed form** `eff = min((G+1)² / [supp + (2G+1)/span], span)` (`priors.py:298`, identical to the existing `capture_eff_length.py:152`). Canonical Laplace `+1` — no tunable constant; `G=0 ⇒ eff=span` (no contraction).
3. **Contained-evidence smooth shrinkage + the gDNA-blindness reasoning.** Calibration's accumulator is fed by **unique mappers only** (`bam_scanner.cpp`: `deposit_to_accumulator` inside `if (is_unique_mapper)`); the EM buffer separately gets multimappers at `1/NH`. So a multimapper-dominated locus (identical paralogs) has little/no **contained** mass, and a pooled seam importing FL-scale gDNA from an adjacent intron can over-contract. The IPR's Laplace `+1` cannot guard this — it shrinks on TOTAL `G` (contained + seam), which the seam itself inflates. PR-A instead shrinks on the **contained mass** `C` toward no-contraction, in **log space**: `w = C/(C+1)`; `eff ← eff^w · span^(1−w)`. `C→0 ⇒ w→0 ⇒ no contraction` (the multimap-blind null, reached smoothly); `C≫1 ⇒ w→1 ⇒ the earned contraction`. PR-B folds this in, anchored at the per-transcript geometric length instead of the locus span (§3.3). See `priors.py:300-320`.

> The current `transcript_capture_eff_lengths` already has a *different* shrinkage — `w = G/(G+n_reg)`, linear in factor space, shrinking on region count (`capture_eff_length.py:154-163`). PR-B **replaces** it with PR-A's contained-mass log-space form. Contained mass is the principled multimap-blindness discriminator; `n_reg` was a region-count proxy.

---

## 2. THE TRANSCRIPT→NODE MAP (the core)

Replace `_exon_region_incidence` (`capture_eff_length.py:41-96`, returns `(inc_t, inc_r)`) with:

```
_transcript_node_incidence(index, region_arrays) -> (inc_t, node_id, weight)
```

over a **unified node space**:
- `node_id < R` → **region node** `r` (genomic length `region_size_bp[r]`), `weight = overlap/region_size_bp ∈ (0,1]`;
- `node_id ≥ R` → **pooled seam node** keyed to left-flank region `r` (FL-scale length), `weight = 1.0`.

The map needs only `index` + `region_arrays` (masses are applied later by `node_id`), so it is a pure annotation-only function — sample-independent, like today.

### 2.1 Rule 1 — region membership: reuse the 20d43f2 bound VERBATIM

Keep the off-by-one-fixed containment idiom exactly (`capture_eff_length.py:72-73`), per-ref CSR slice `[lo0, hi0) = ref_offsets[rid:rid+2]`:

```python
lo = lo0 + searchsorted(ends[lo0:hi0],   a, side="right")   # first region whose end_r > a (contains a)
hi = lo0 + searchsorted(starts[lo0:hi0], b, side="left")    # last region whose start_r < b
```

The lower bound is read from region **ENDS** (`side="right"`). Commit `20d43f2`: reading from `starts` skips the region containing `a` when `a` is interior, dropping fully-contained spans (the nascent sink dropped 313K→184K, −41%, on the fix). **MANDATORY: preserve this exactly.**

Footprint **sources** (unchanged from today — this is what makes mature/nascent diverge; see §2.4):
- **mature** → EXON intervals from `intervals.feather`, `interval_type==EXON & t_index>=0` (`capture_eff_length.py:78-82`).
- **nascent** → full span `[start, end)` from `t_df[is_synthetic]` (`capture_eff_length.py:84-92`). Synthetic nRNA emit **no** EXON intervals (`index.py:509-510`), so the full-span branch is their only source and covers introns too.

### 2.2 Rule 2 — partial-overlap weight (NEW)

For a region `r` whose `[start_r, end_r)` is clipped by the transcript interval `[a, b)`:

```
overlap_r = min(end_r, b) − max(start_r, a)
weight_r  = overlap_r / region_size_bp[r]      # ∈ (0, 1]
```

The weight applies to **both** `g` and `L`, so the per-node density `g/L` is **weight-invariant** (only the magnitude scales). A multi-exon transcript hitting the same merged region accumulates the weight once per exon (`np.add.at`); the per-exon overlaps are disjoint by construction, so summing is correct. A region fully inside one exon → `overlap = region_size_bp` → `weight = 1.0`. This is NEW logic — the current map maps to whole regions with no overhang weighting.

> The `_add` COO emitter must therefore emit a `weight` column; the old `arange(lo, hi)` regions no longer all get weight 1 — the leftmost/rightmost get the partial weight, interior regions get 1.0.

### 2.3 Rule 3 — pooled-seam inclusion (NEW)

Include the interior seam node `R+r` (between regions `r` and `r+1`) **iff ALL of:**
- `r` is in the footprint AND `r+1` is in the footprint, AND
- `ref_id[r] == ref_id[r+1]` (same-reference adjacency, mirror `priors.py:240`).

The **leftmost** included region drops its left-outer seam; the **rightmost** drops its right-outer seam (a crossing overhanging past the transcript's first-exon-start / last-exon-end is excluded — its other half belongs to a boundary outside the footprint). Consequences:
- **single-region transcript → 0 seam nodes** (mandatory test).
- **mature**: footprint = exon regions ⇒ seams only between **adjacent kept exon-regions**. A mature exon↔exon splice junction is **NOT a genomic seam** — the intronic region *between* two mature exons is simply absent from the footprint (mature does not span introns), so no seam exists there. Exon↔exon junctions are **deferred SJ nodes** (`splice_junction_node_architecture.md`), out of PR-B scope.
- **nascent**: footprint = all regions in `[lo, hi)` ⇒ ALL interior seams including intron↔exon and intron↔intron — the nascent-only seam signal.

Seam node geometry (mirror `priors.py:241-244`, **pin to the `min`-of-flanks form** for cross-consumer bit-identity, NOT the redesign §7.2 left-flank form — see §5/R-extract):
```
seam_mass[r] = side_right[r] + side_left[r+1]            # POOL both halves
seam_len[r]  = min(gdna_boundary_len[r], gdna_boundary_len[r+1])   # FL-scale, tighter flank
```

### 2.4 Rule 4 — strand decode + the ≥1-node assertion (MANDATORY)

`t_df.strand ∈ {1, 2}` is the `Strand` IntEnum (POS=1, NEG=2; `types.py:34-36, 73-74`), **NOT** signed `{+1,−1}` and **NOT** the int8 `strand_class` on `RegionArrays` (TS_NEG=−1, TS_AMBIG=2; `signature.py:74-77`). A naive `strand>0 → POS` decode is True for BOTH strands.

**In THIS design the strand bit is UNUSED.** Inclusion is overlap-driven (Rule 1); the seam is a single pooled node with no left/right side to choose. So do NOT write a `{1,2}` strand decode — it would be dead, untested code (Risk F5). What is mandatory is the **coverage assertion**:

```python
assert np.unique(inc_t).size == n_t   # every transcript maps to ≥1 node
```

This replaces the current silent empty-return (`capture_eff_length.py:94-95`) with a raised error. It catches off-by-one regressions, ref-name misses, and any future inclusion bug — NOT a strand bug (nothing reads strand). Do not frame it as a strand guard.

### 2.5 Inclusion is interval/span-driven, NEVER signature-driven (the 0xA trap)

**Critical clarification the prior sketch language invites getting wrong.** A region can have a mixed signature, e.g. `0xA = BIT_EXON_POS | BIT_INTRON_NEG` (`signature.py:31-34`): it is **exon for the + transcript and intron for the − transcript simultaneously**. `coarse_type_array` (`signature.py:147-160`) collapses this to ONE coarse type (exon-wins-over-intron), and `transcript_strand_class` collapses `0xA` to `TS_AMBIG`.

The design docs invoke `BIT_EXON_*/BIT_INTRON_*` and "include intron region-nodes for nascent" as if inclusion were a **signature lookup**. **It is NOT, and must NOT be.** Inclusion is driven by **which intervals exist**: mature from `intervals.feather` EXON rows, nascent from the `t_df` full span. "EXON regions for mature, EXON+INTRON for nascent" is a *consequence* of which intervals exist for the transcript, not a per-region `signature`/`coarse_type` decision. If an implementer switches to signature-based inclusion, a `0xA` region is mis-assigned for one of the two overlapping opposite-strand transcripts.

**MANDATE:** `signature` / `coarse_type` / `strand_class` are **never consulted** for transcript→node membership. Membership = the searchsorted overlap (Rule 1) against the transcript's own intervals/span. Add a regression with an overlapping opposite-strand pair sharing a `0xA` region (§2.6 fixture).

### 2.6 The mandatory regression battery + the NEW fixture

The existing fixture `_MISALIGNED_GTF` (`tests/calibration/test_capture_eff_length.py:26-31`) is **same-strand, no nRNA twin, no opposite-strand overlap** — it exercises NONE of the new rules. PR-B requires a NEW fixture extending it with: (a) a multi-exon POS transcript whose intron is a distinct partition region; (b) a multi-exon NEG transcript; (c) a synthetic nRNA twin (same `[start,end)`, ≥1 interior intron region); (d) an overlapping opposite-strand pair sharing a `0xA` region.

| Surface | Hazard | Test |
|---|---|---|
| region bound `[lo,hi)` | reusing `starts` for `lo` drops contained exons (20d43f2) | **port** `test_incidence_maps_every_transcript`, `test_merged_region_interior_exon_is_mapped` (`test_capture_eff_length.py:39-95`) to the node map |
| partial overlap | unweighted inclusion double-counts edge gDNA owned by a neighbor; conservation breaks | `test_partial_overlap_span_equals_exonic` (restate the broken `test_incidence_len_t_ge_exonic` **inequality** as an **equality**); `test_partial_overlap_g_over_l_invariant` (`(w·g)/(w·L)==g/L`); add the both-sides-overhang case (exon ⊋ region AND region ⊋ exon) |
| outer-drop | first region's left seam / last's right seam double-counts overhang | `test_single_region_zero_seams` (1 region ⇒ 0 seams); `test_multi_region_internal_seams_only` (k regions ⇒ k−1 interior seams, no outer) |
| encompassed | interior region fully traversed; pooled-seam mass unchanged | `test_encompassed_region_seam` (3-region footprint, middle fully encompassed ⇒ `w_middle==1.0`, two interior seams) |
| mature/nascent divergence | the headline correctness property | `test_mature_nascent_node_sets_differ` (`mature_nodes ⊊ nascent_nodes`; extra = intron region-nodes + intron-bearing seams; `factor_nascent < factor_mature`) |
| mixed-signature 0xA | signature-based inclusion mis-assigns one strand | `test_opposite_strand_shares_0xA_region` (both transcripts map correctly; membership is overlap-driven, not signature) |
| coverage | a transcript maps to 0 nodes | `test_every_transcript_ge_one_node` (`np.unique(inc_t).size == n_t`) |
| capture-off | factor ≈ 1, NOT bit-identical | `test_capture_off_factor_one` (§3.4) |

---

## 3. THE PER-TRANSCRIPT EFF-LEN

Per transcript `c`:

```
eff_em_c = min( GEOMETRIC_c · factor_c , GEOMETRIC_c )
factor_c = factor_raw_c ^ w_c                              # contained-evidence shrinkage, log space
```

### 3.1 GEOMETRIC_c — the FL-marginal length (DO NOT recompute)

`GEOMETRIC_c = fl_eff_lengths[c]`, the existing FL-marginal `effective_lengths` passed in at `pipeline.py:450-451`. It is built once for ALL transcripts in `pipeline.py:438-444`:
```
exonic_lengths    = index.t_df["length"]                              # pipeline.py:438
effective_lengths = rna_fl.compute_all_transcript_eff_lens(exonic_lengths)   # pipeline.py:440
```
`compute_all_transcript_eff_lens` (`frag_length_model.py:493-545`) is length-agnostic — it integrates the RNA FL pmf against whatever length array it is handed. **Mature/nascent already diverge here**, carried entirely by `t_df["length"]`:
- **mature**: `length = Σ(exon.end − exon.start)` = **exonic** (`transcript.py:107-111`) ⇒ `GEOMETRIC = exonic − RNA-FL-integral`.
- **nascent synthetic**: `length = e − s` = **genomic span** (single full-span exon, `index.py:386-397`) ⇒ `GEOMETRIC = span − RNA-FL-integral`.

PR-B multiplies this by the factor; it **never** recomputes the geometric length. The capture-off identity is checked against this `fl` (factor 1 ⇒ `eff_em == fl`).

### 3.2 factor_raw_c — the pooled-seam IPR (reuse PR-A's closed form)

Accumulate over the transcript's node edges (`np.add.at` over `(inc_t, node_id, weight)`):

```
G_c    = Σ_nodes        weight · g_node                       # g_region=mass_gdna_contained; g_seam=pooled seam mass
supp_c = Σ_region_nodes weight · g²/region_size_bp            # contained² / genomic length
       + Σ_seam_nodes   1.0    · seam_mass²/seam_len          # FL-scale participation (NOT weighted)
span_c = Σ_region_nodes weight · region_size_bp               # GENOMIC; seams add 0  (replaces Σ gdna_geom_len at :147)
```

Then, **verbatim** PR-A / the existing `:152`:
```
eff_ipr_c   = min( (G_c + 1)² / [ supp_c + (2·G_c + 1)/span_c ], span_c )
factor_raw_c = eff_ipr_c / span_c       # ∈ (0,1]; 1 ⇒ no contraction
```

NEW vs reused:
- **REUSE** the Laplace IPR (`priors.py:298` == `capture_eff_length.py:152`) and the cap.
- **NEW** vs the current code: region `g` = `mass_gdna_contained` (contained-only, NOT `_transport_boundary_flux`); region `L` = `region_size_bp` (genomic, NOT `gdna_geom_len`); seam term added; `span` from `region_size_bp` (not `gdna_geom_len`); weighting (§2.2).

### 3.3 w_c — contained-evidence smooth shrinkage (fold in PR-A, anchored at GEOMETRIC)

```
C_c = Σ_region_nodes weight · (mass_gdna_contained + mass_rna_contained)    # footprint contained mass (gDNA+RNA)
w_c = C_c / (C_c + 1)                                                       # one pseudo-obs, magic-free
factor_c = factor_raw_c ^ w_c            # log space: exp(w·log(factor_raw) + (1−w)·log(1.0))
```

This mirrors `priors.py:315-320` exactly with one anchor difference: PR-A interpolates `eff_len` toward the locus **span** in log space; PR-B carries a dimensionless **factor**, so it interpolates the factor toward **1.0** (≡ `eff_em` toward `fl`, the geometric anchor). Algebraically the same shrinkage, anchor swapped.

> Use **gDNA + RNA** contained mass for `C_c`, matching PR-A's `contained_ev = gdna_contained + rna_contained` (`priors.py:315`), for cross-consumer consistency. (The verify notes were ambiguous "contained mass"; pin to PR-A's exact form.)

Behavior:
- `C → 0` (multimap-blind / zero-gDNA library) ⇒ `w → 0` ⇒ `factor → 1` ⇒ `eff_em = fl`, smoothly (no cliff). This is the identical-paralog collapse fix ported to the transcript path, and the **weak-SS zero-gDNA-phantom** guard.
- `C ≫ 1` (real captured locus) ⇒ `w → 1` ⇒ the earned IPR contraction.
- `G = 0` ⇒ `factor_raw → 1` AND `w → 0` (no contained) ⇒ double-guarded to `fl`.

### 3.4 Mature vs nascent divergence + capture-off identity

**Divergence.** The footprint `F` is the only thing differing between a mature transcript and its nascent twin (geometric `fl` already differs via §3.1). Nascent `F ⊋` mature `F`: nascent adds the **intron region-nodes** (depleted: large `L`, tiny `g` ⇒ inflate `span` without adding to `G`/`supp` ⇒ **stronger** contraction) plus the **intron-bearing interior seams**. The divergence is **dominated by the intron region-nodes**, not the seams.

**Capture-off identity.** Uniform gDNA enrichment ⇒ per-node density `g = ρ·L` ⇒ `(ρ·span+1)²/(ρ²·span + 2ρ·span + 1) → 1` ⇒ `factor → 1` ⇒ `eff_em == fl`. Seam nodes at uniform density sit at the same `ρ` and enter consistently, so the cap pins uniform transcripts to `fl`. **NOT bit-identical** on real data (gDNA is never perfectly uniform; PR-1 measured whole-genome factor 0.998822) — so real capture-off goldens move ~0.1% (§6 `--update-golden` discipline). The exact-uniform synthetic unit test IS bit-identical (`factor == 1.0` to `1e-9`).

---

## 4. CONSERVATION (proven, per transcript)

For one transcript with footprint regions `r ∈ F` (a contiguous run `[lo, hi)`):

### 4.1 Law (L) — genomic length, EXACT

```
span_c = Σ_{r∈F} weight_r · region_size_bp[r] + Σ_{seams} 0
       = Σ_{r∈F} (overlap_r / region_size_bp[r]) · region_size_bp[r]
       = Σ_{r∈F} overlap_r
       = exonic_length  (mature)   /   genomic_span  (nascent)        EXACTLY
```

Seams contribute **0** genomic bp — their FL-scale `seam_len` lives only inside the dimensionless `supp` term. This is exactly the conservation that `gdna_geom_len` (= region_eff + both sides ≈ 1.065× span) violated. **Test:** restate `test_incidence_len_t_ge_exonic` (inequality, `:50-69`) as an **equality** `span_c == exonic` / `== span` (no `±FL`).

### 4.2 Law (M) — gDNA mass, no double-count

```
G_c = Σ_{r∈F} weight_r · mass_gdna_contained[r]                       (region nodes — CONTAINED only)
    + Σ_{(r,r+1)∈F internal} (mass_gdna_right[r] + mass_gdna_left[r+1]) (pooled seams)
```

No double-count, by construction:
- region nodes carry `mass_gdna_contained` **only** — never `+left+right` (the `+left+right` fold was `_transport_boundary_flux`, deleted).
- the seam family visits each internal boundary's two sides **exactly once** (`mass_gdna_right[r]` and `mass_gdna_left[r+1]` are distinct array entries; the substrate `left`/`right` views partition every internal side bijectively, `result.py:44-47`). The seam is keyed to its left-flank region so it is emitted at most once per included adjacent pair.
- **outer-drop** prevents overhang double-count: a seam enters iff BOTH flanks ∈ F, so the leftmost region's left-outer half (belonging to a boundary outside F) and the rightmost's right-outer half never enter.

This is the locus-level PR-A mass law (`priors.py:222-223`, proven `isclose=True`) restricted to the transcript footprint. **Test:** `test_transcript_gdna_mass_no_double_count` — distinct known scalars per side; assert `G_c` equals the explicit sum; single-region transcript ⇒ `G_c == contained[r]` exactly (zero seam contribution).

---

## 5. RISKS (adversarial findings, each with mitigation) + the honest ROI call

### R-ROI (SERIOUS — the headline must be honest): PR-B does NOT meaningfully fix the nascent sink

The 3-pool decomposition at the flagship-with-nascent condition (`effective_length_redesign_plan.md:28-39`): `net_gdna_to_nrna = +125,907`, `net_gdna_to_mrna = +107,952`, **`net_nrna_to_mrna = −9,361`**. The nascent surplus is ~93% **gDNA-sourced** (PR-A territory — the locus `gdna_eff_len`). The **mature↔nascent face** — the only face PR-B's mature/nascent divergence can move — is `−9,361`, **under 4%** of the gDNA→nascent flux.

Worse, the *direction* is a hazard, not a safety property: per-transcript `effective_lengths_em` (PR-B) and the locus `gdna_eff_len` (PR-A) are sibling components in the SAME EM softmax (`em_solver.cpp:1838-1854`: `log_weight = log θ − log_eff_len`). Contracting a nascent transcript's eff-len raises its per-position rate ⇒ it competes **harder** ⇒ pulls mass TOWARD nascent. On a library that is already **over**-assigning nascent (`nrna_net_surplus = +135,268`), `factor_nascent < factor_mature` (§3.4) points toward **worsening** nascent over-assignment. The `nascent_efflen_invariant_bug` memory documents exactly this failure mode (contraction below mature inflated the sink 172K→1K causal proof). **Do NOT assert `factor_nascent < factor_mature` as a safety property** — it is the hazard direction; it is *correct geometry* but must be net-flow-gated, not assumed beneficial.

**Mitigation / reframe:** PR-B's defensible value is **conservation + consistency + dead-code enablement**, NOT "fix the nascent sink." It removes a genuine correctness defect (the `gdna_geom_len` genomic-length-conservation violation on the EM divisor) and unblocks deleting the retired transport path (§6). Frame the PR that way; gate it as net-flow-**neutral** (§6 gate), with the explicit expectation that the headline numbers barely move (≈ the −9,361 scale). The actual nascent-sink fix is **EM-side** (per-fragment nascent-vs-mature competition, isoform-swap) — out of PR-B scope.

### R1 — mixed-signature 0xA region (SERIOUS, latent doc-trap): mitigated by §2.5

If an implementer follows the docs' `BIT_EXON_*/BIT_INTRON_*` language and switches to signature-based inclusion, a `0xA` region is mis-assigned for one of two overlapping opposite-strand transcripts. **Mitigation:** §2.5 mandates inclusion stays interval/span-driven exactly as today, `signature`/`coarse_type` never consulted; `test_opposite_strand_shares_0xA_region` (§2.6) exercises it. The existing fixture cannot — it is same-strand.

### R2 — strand `{1,2}` landmine (UNFOUNDED for this design): no strand decode

`strand>0→POS` would empty ~42% of nascent footprints — but ONLY if strand gates inclusion. This design never reads strand (§2.4). **Mitigation:** do not write the `{1,2}` decode (dead code, Risk F5). Keep the ≥1-node assertion as off-by-one/ref-miss hygiene, NOT framed as a strand guard.

### R3 — off-by-one re-opened by partial-overlap weighting (MANAGEABLE): mitigated by §2.2/§4.1

The 20d43f2 fix is preserved verbatim (§2.1), but partial-overlap weighting is NEW logic on top. The conservation equality (§4.1) holds only if each exon's overlap is clipped to the region. **Mitigation:** the equality test (§2.6) + the both-sides-overhang test pin it.

### R4 — risk to PR-A's 5.62% leak via shared state (UNFOUNDED): fully decoupled

PR-A writes per-LOCUS `gdna_eff_len` (`priors.py:325` → `em_solver.cpp:1888, 2165`); PR-B writes per-TRANSCRIPT `effective_lengths_em` (`pipeline.py:450` → `em_solver.cpp:1891, 1839`). Different arrays, different EM parameters, no shared mutable output. The 5.62% leak is a function of `gdna_eff_len`, which PR-B does not touch.

### R-extract — risk to PR-A's leak via shared CODE (MANAGEABLE if sequenced): the only real coupling

The pooled-seam math is duplicated logic between `priors.py:239-244` and the new map. **Decision: for the first cut, DUPLICATE the ~6 lines in `capture_eff_length.py` — do NOT extract.** Refactoring a shipped, golden-pinned hot path on PR-B's schedule couples PR-B to re-baselining PR-A, for ~6 trivial lines. **If** extracting later (when a third consumer appears): (a) do it as a standalone behavior-preserving refactor commit with `tests/calibration/test_priors.py` + calibration goldens green BEFORE any PR-B logic; (b) **pin the shared helper to PR-A's `min`-of-flanks seam length** (`priors.py:242`), NOT the redesign §7.2 left-flank form — adopting left-flank would perturb PR-A's goldens for no measured benefit. PR-B uses `min`-of-flanks (§2.3) regardless.

---

## 6. IMPLEMENTATION SEQUENCE

PR-sized steps. PR-B is its own commit on top of shipped PR-A; the net-flow bisects cleanly to it.

**Step 0 — record baselines (post-PR-A).** Run the `calibration-benchmark` skill; record per condition `net_gdna_to_nrna`, `net_gdna_to_mrna`, `net_nrna_to_mrna`, `surplus_gdna`, `surplus_nrna` (metric names at `scripts/debug/qsweep_3pool.py:70-72`), the nRNA sink, and the flagship leak_frac.

**Step 1 — the NEW fixture.** Extend `_MISALIGNED_GTF` (§2.6): multi-exon POS with a distinct intron partition region, multi-exon NEG, a synthetic nRNA twin, an opposite-strand pair sharing a `0xA` region. This fixture is NOT optional — it is the only thing that exercises every new rule.

**Step 2 — `_transcript_node_incidence`** (replace `_exon_region_incidence`). Emit `(inc_t, node_id, weight)`; region membership reuses §2.1 verbatim; add partial-overlap weight (§2.2); add pooled-seam inclusion (§2.3); add the ≥1-node assertion (§2.4). Port + add the §2.6 battery; green.

**Step 3 — rewrite `transcript_capture_eff_lengths`.** Delete the `_transport_boundary_flux` import (`:31`) and call (`:127-134`) and the `gdna_geom_len` divisor (`:135-136`). Build region nodes (`mass_gdna_contained`, `region_size_bp`) + pooled seam nodes (DUPLICATE the `priors.py:239-244` math, `min`-of-flanks). Weighted `np.add.at` over the node edges → `G, supp, span, C` (§3.2-3.3). REUSE the IPR (`:152`) + cap (`:164`). REPLACE the `w=G/(G+n_reg)` shrinkage with the contained-evidence log-space form (§3.3).

**Step 4 — local tests + full suite.**
```
pytest tests/calibration/test_capture_eff_length.py -v   # the bug-surface battery (§2.6)
pytest tests/                                             # full suite (1113 baseline) green
```

**Step 5 — the net-flow gate (PRIMARY).** Run `calibration-benchmark`, at minimum `gdna300 × ss{0.99,0.50} × cap_on × {nrna_none, nrna_rnd}` plus `ss0.65 × gdna_none` (the zero-gDNA-phantom regression). **Proceed iff:**
1. `net_gdna_to_{nrna,mrna}` per condition do NOT worsen (≥ baseline, within net-flow noise).
2. No high-gDNA-tail inflation: flag any locus with `NEW/OLD eff_em > 1.5×`; abort if many loci inflate.
3. `net_nrna_to_mrna` moves only second-order (≈ −9,361 scale). A LARGE swing is a RED FLAG of a footprint-map regression (strand/off-by-one), NOT a win.
4. nRNA sink moves toward truth or holds; nascent must NOT collapse below mature beyond what net-flow tolerates.
5. The `ss0.65` zero-gDNA-phantom tolerance holds (the contained-evidence shrinkage → factor 1 when `C→0`).

**Step 6 — capture-off golden, `--update-golden` discipline.** The exact-uniform unit test (`test_capture_off_factor_one`) is bit-identical (`1e-9`). The REAL capture-off benchmark goldens move ~0.1% (PR-1's 0.998822). Land the golden update as its OWN `--update-golden` commit AFTER Step 5 confirms direction, so a golden move bisects to PR-B alone, never conflated with a logic regression. Document the expected ~0.1% magnitude.

**Step 7 — cleanup (SEPARATE follow-up commit, after PR-B green).** After PR-B removes the `capture_eff_length.py:31, 127-135` uses, both retired symbols have ZERO `src/` consumers (grep `_transport_boundary_flux|gdna_geom_len` over `src/ tests/ scripts/`):
- DELETE `_transport_boundary_flux` (`priors.py:103-159`) — only remaining reference is its own def + doc comments.
- DELETE `gdna_geom_len`: the `derive.py:30, 55-56` computation, the `DerivedDensity`/`CalibrationResult.gdna_geom_len` field (`result.py:61, 91`), the `calibrate.py:312` assignment.
- UPDATE breaking references: `tests/calibration/test_priors.py` (imports `_transport_boundary_flux` `:8`, `gdna_geom_len` kwarg `:16-282`, transport test block `:288-328`); `test_result_schema.py:23,43,56,68`; `test_calibrate.py:75-76`; `test_nrna_double_counting.py:248` (comment); debug scripts `scripts/debug/{dissect_efflen,enrichment_bedgraph,phase3_gdna_distribution,diag_imputation_truth}.py`.
- KEEP `gdna_boundary_len` (`result.py:66`) — the seam-node length in BOTH consumers.

Pure dead-code removal in its own commit ⇒ bisects independently from the logic change.

---

## 7. SCOPE CALL

**Ship a NARROWED PR-B now: conservation hygiene + dead-code enablement. Drop the "fixes the nascent sink" framing.**

Reasoning. The plan as originally sketched rested on a thesis the code refutes — that the transcript eff-len move materially fixes the nascent sink. It does not: `net_nrna_to_mrna = −9,361` is <4% of the gDNA→nascent flux (`+125,907`, PR-A's domain), and PR-B's mechanism on nascent (contract more → compete harder) points toward worsening the documented over-assignment, not fixing it (R-ROI). The genuine, defensible win is real and worth shipping: the transcript EM eff-len is currently on the retired `_transport_boundary_flux` + `gdna_geom_len` path that violates genomic-length conservation (≈1.065× span on the EM divisor) — PR-B fixes that, makes mature/nascent diverge through the (correct) node set, and unblocks deleting the retired path.

So:
- **DO PR-B now**, scoped as conservation hygiene, gated **net-flow-neutral** (Step 5), with duplicated (not extracted) pooled-seam math (R-extract), and WITHOUT the `factor_nascent < factor_mature` "safety" claim.
- **DEFER** any nascent-sink improvement to the **EM-side** work (per-fragment nascent-vs-mature competition / isoform-swap) — that is where the −9,361 face lives.
- **DO NOT** ship the original framing (sink-fix + shared-helper extraction on PR-B's schedule + the hazard-as-safety claim). It conflates a hygiene migration with a quant fix and refactors a golden-pinned path on the wrong schedule.

**Acceptable alternative if schedule-constrained: defer PR-B entirely behind the EM-side nascent work.** The mature↔nascent sink is EM-side; fixing the eff-len first risks a confounded benchmark and a wrong-sign nascent regression. The genomic-length-conservation defect on the EM divisor is real but second-order on quant — it can wait. Lowest risk, highest information. Recommended order of preference: (1) narrowed PR-B as hygiene now (if the dead-code deletion is wanted); (2) defer behind EM-side nascent work; (3) NOT the original plan.

---

### Files & anchors (all verified)
- Target: `src/rigel/calibration/capture_eff_length.py` — `:31` import (DELETE), `:41-96` `_exon_region_incidence` (REPLACE → `_transcript_node_incidence`), `:127-136` transport+geom (DELETE), `:138-148` accumulation (REWRITE weighted over nodes), `:150-153` IPR (KEEP), `:154-163` shrinkage (REPLACE w/ `C/(C+1)` log-space), `:164` cap (KEEP).
- Reference to mirror: `src/rigel/calibration/priors.py` — `:215-230` design comment, `:231-245` pooled-seam build, `:259-263` region+seam supp, `:298` Laplace IPR, `:315-320` contained-evidence log-shrinkage.
- Geometric length (DO NOT recompute): `src/rigel/pipeline.py:438-452`; `src/rigel/frag_length_model.py:493-545` (`compute_all_transcript_eff_lens`); nRNA span length `src/rigel/index.py:386-397`; mature exonic length `src/rigel/transcript.py:107-111`.
- Node geometry: `src/rigel/calibration/region_arrays.py:55` (`region_size_bp`), `:128-185` (adjacency); signature bits `src/rigel/calibration/signature.py:31-34, 147-177`.
- Strand enum: `src/rigel/types.py:34-36, 73-74` (POS=1, NEG=2).
- Synthetic nRNA emit NO exon intervals: `src/rigel/index.py:509-510`.
- SJ-node deferral: `docs/calibration/splice_junction_node_architecture.md`.
- EM softmax (PR-A/PR-B sibling components): `src/rigel/native/em_solver.cpp:1838-1854, 1888, 1891`.
- To retire after PR-B: `src/rigel/calibration/priors.py:103-159` (`_transport_boundary_flux`), `src/rigel/calibration/derive.py:30, 55-56`, `src/rigel/calibration/result.py:61, 91`, `src/rigel/calibration/calibrate.py:312`.
- Tests to port/restate: `tests/calibration/test_capture_eff_length.py:39-95`.
- Prior sketch refined here: `docs/calibration/effective_length_redesign_plan.md` §7.4, §7.7.
- Validation: `calibration-benchmark` skill (3-pool net-flow; metric names `scripts/debug/qsweep_3pool.py:70-72`; suite ~37s); rigel env `source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel`.
