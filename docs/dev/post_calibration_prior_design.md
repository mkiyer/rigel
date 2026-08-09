# Post-calibration → the EM prior — restoring the object model

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When an item settles, MOVE it to its permanent
    home (a ruling → `DESIGN.md`, a derivation → `EQUATIONS.md`, a lesson → `TRAPS.md`) and delete it
    here in the same edit.

    Opened 2026-08-08, from the owner's ruling: **an edge owns the fragments that cross it. A node owns
    only the fragments contained in it. Nothing is re-attributed.**

---

## ✅ STATUS — BOTH PHASES LANDED 2026-08-08

| | |
|---|---|
| suite | ✅ **0 failed / 3,253 passed** / 2 skipped / 9 xfail · lint clean |
| deleted | `edge_owner_nodes`, the intergenic re-key, `_component_node_arrays`, `_mass_where_there_is_opportunity`, `_left_keyed_edge_arrays` |
| added | `_region_locus_shares`, `_edge_locus_shares`, `contended_edges`, `_sum_by_locus` |
| schema | ⛔ **`payload_schema_digest` UNCHANGED** — no accumulator change, no cache rebuild |
| flgap | ⭐ **numerically identical**: `S−F` 0.0302 / 0.0084 / 0.0058 / 0.0035 and `O−F` 0.0406 / 0.0146 / 0.0081 / 0.0041, matching the pre-refactor values on all four conditions |
| goldens | 2 of 21 scenarios moved, by ≤ **4.8e-4** — the per-object contraction where a node used to own two lines |
| falsification | ⭐ **8/8 injected defects caught**, after the first pass found **2 gate holes** (below) |

⛔ **Two things the design predicted that measurement corrected:**

1. **Contended edges are NOT impossible — there are 20–34 per flgap condition**, carrying 7–252 units of
   mass against ~2.3 M total (≈0.01 %). The structural argument (loci are separated by intergenic
   sequence) holds for the overwhelming majority and not universally. ⭐ This is exactly why
   `contended_edges` reports rather than asserts; had it been an assert, the run would have died.
2. **`node_right_edge[r] == r` on a single reference**, so every existing fixture was blind to the
   difference between an edge index and a left-node index — the Phase-2 conversion was *untested* and a
   defect there would have been silent on real multi-chromosome data. Found by perturbation; closed by
   `test_the_boundary_incidence_is_an_EDGE_index_not_a_left_node_index` on a two-reference fixture.
   ⭐ The other hole: the per-object contraction clamp had no gate at all, because a uniform-density
   fixture cannot separate it from a folded one. Closed by
   `test_the_contraction_is_applied_PER_OBJECT_not_over_a_folded_total`.

---

## 0. ⭐⭐⭐ THE ONE THING TO READ

`assemble_priors` is the only place in `src/` that **folds a line's mass into a node's total**. Every
other consumer of the calibration keeps a line as its own object. That fold is why the prior needs an
owner map, why the owner map needs an intergenic re-key, and why the whole step reads as a heuristic.

**Delete the fold and all three go away.** A locus collects nodes *and* edges; each object keeps its own
mass; nothing moves. ⛔ **No accumulator schema change is required** — this is a change to attribution,
not to storage.

---

## 1. THE INVARIANT

    a NODE owns the fragments CONTAINED in it          node_contained_count, mass_*_node
    an EDGE owns the fragments that CROSS it           edge_unspliced_mass, mass_*_edge
    nothing is ever re-attributed from one to the other

`density_model.py` already states it: *"A contiguous edge is a first-class object with its own axis."*
The accumulator, `CalibrationResult`, the sweep, `strand_deconv` and `density_model` all honour it.

⭐ **`transcript_capture_eff_lengths` honours it too, and that is the model to copy.** A transcript
enumerates its own nodes *and* its own lines, and contracts each object on that object's own density:

```python
np.add.at(num,       rt, np.minimum(contained_m[rr] * inv, contained_S[rr]))   # its nodes
np.add.at(num,       bt, np.minimum(seam_m[br]    * inv, seam_S[br]))          # its lines
```

⚠ `_left_keyed_edge_arrays` looks like an attribution and is not: it builds a node-*indexed* array only
because `_transcript_node_incidence` emits a **left-node index** for each boundary. It is an indexing
convention. No line's mass is ever added into a node's total.

---

## 2. WHAT THE EM ACTUALLY CONSUMES — the target semantics

    em_solver.cpp:698    n_gdna = raw_counts[gdna_idx]        soft count, ONE multi-locus
                         G = n_gdna + a_g                     ← a_g is gdna_prior_count

A **multi-locus is a connected component of transcripts linked by shared fragments**, so a fragment is a
candidate in **exactly one** multi-locus. Combined with the owner's ruling that *any* overlapping
fragment is a candidate:

⭐⭐ **`gdna_prior_count[L]` = the number of gDNA fragments that are candidates at L, each counted ONCE.**
A fragment that starts intergenic and runs into the gene is a candidate, and must count.

⛔⛔ **THEREFORE `F` — `prior_vs_oracle`'s first-base fragment truth — IS THE WRONG YARDSTICK FOR THIS
ARM,** and a measured `S−F` is not by itself an error. `F` credits a fragment to whichever locus holds
its first base, so it drops exactly the straddlers the EM counts. Measured on `flgap_long` capture ON:
`net Σ(S−F) = +141,783 = +3.01 %` of `F`, with `net/abs = 0.998–1.000` — a pure one-way excess, which is
the signature of a population `F` cannot see rather than of a bias. ⭐ The boundary-crossing mass
predicts it to **1.002 / 1.000** off capture from geometry alone, which is now *confirmation of the
semantics gap*, not a bug report. ⚠ `prior_vs_oracle.py`'s docstring must be corrected to say so.

---

## 3. THE RULE — a locus's objects

**Nodes.** Unchanged: `_project_regions_to_loci`, genomic-overlap share.

**Edges.** ⭐ **A locus's edges are the edges that touch its nodes** (owner, 2026-08-08). Every interior
node contributes both of its lines; a locus of `k` nodes therefore has `k+1` lines, its two outer ones
included.

⭐ **The outer lines are unambiguous, and that is structural**: a locus is bounded by intergenic
sequence, intergenic regions carry no transcripts, so no two loci can contend for a boundary line.

    share(e, L) = max( share(lo(e), L), share(hi(e), L) )

⛔ **Gate it rather than trust it**: assert `Σ_L share(e, L) ≤ 1 + eps` for every edge and **report** any
violation with a count. The owner's argument says it cannot fire; a node split across two loci
(`test_node_split_between_two_loci`) is the only shape that could, and if it ever does we want the number,
not a silent renormalisation.

⭐ **This reproduces today's edge set exactly**, so the change is a *refactor* on real data: the current
code already reaches all `k+1` lines — the left outer one via the intergenic re-key, the right outer one
via default left-keying. What changes is that the rule becomes stated instead of emergent.

---

## 4. THE TWO OUTPUTS

```
gdna_prior_count[L] = Σ_{nodes r} share(r,L) · mass_gdna_node[r]
                    + Σ_{edges e} share(e,L) · mass_gdna_edge[e] · q[e]

rna_prior_count[L]  = Σ_{nodes r} share(r,L) · mass_rna_node[r]
                    + Σ_{edges e} share(e,L) · (mass_rna_edge[e] − mass_rna_spliced_edge[e]) · q[e]

gdna_eff_len[L]     = clamp( w·elen_L + (1−w)·span_L ),   w = C/(C+1)
    elen_L = Σ_{nodes} share·min(m_r/ρ_ref, S_r) + Σ_{edges} share·min(m_e/ρ_ref, S_e)
    span_L = Σ_{nodes} share·S_r                 + Σ_{edges} share·S_e
```

with `q = edge_mass_per_crossing` (the accumulator's own `mass/count`, converting an object-incidence
total to a fragment count) and `C` the projected contained mass.

⭐⭐ **`elen_L` is now the SAME operation `transcript_capture_eff_lengths` performs, over a different
object set.** That is the alignment the owner asked for, and it argues for one shared helper:

    contract(node_ids, node_share, edge_ids, edge_share, ρ_ref) -> (num, span)

called by the transcript path over a transcript's objects and by the locus path over a locus's. ⚠ The
transcript path additionally imputes a **splice-junction** object, which has no entry on the edge axis;
that stays where it is and is passed in.

---

## 5. WHAT CHANGES

| file | change |
|---|---|
| `calibration/priors.py` | ⛔ **delete `edge_owner_nodes`** and its re-key. ⛔ **delete `_component_node_arrays`** — its whole job was the fold. Replace `_conserved_count` with a per-locus edge sum. `assemble_priors` gains an edge→locus projection beside the node one |
| `calibration/priors.py` | `elen` / `span` are built from per-node **and** per-edge terms rather than from folded node arrays |
| `calibration/capture_eff_length.py` | ⚠ **phase 2, optional**: `_transcript_node_incidence` emits a left-node index per boundary, which is why `_left_keyed_edge_arrays` exists. Emitting an **edge index** deletes that helper and puts both paths on one axis. Not required for correctness; do it only if the shared `contract()` helper lands |
| `scripts/design/prior_vs_oracle.py` | ⛔ its docstring asserts `F` is the gDNA arm's exact truth. Correct it: `F` is a **first-base** count and the EM counts any overlapping candidate, so `S−F` carries a known one-way term. Consider adding an `F_overlap` arm that counts a gDNA fragment in every locus it overlaps |

**Deleted, net:** `edge_owner_nodes`, `_component_node_arrays`, the intergenic re-key, and (phase 2)
`_left_keyed_edge_arrays`.

---

## 6. WHAT DOES **NOT** CHANGE

⛔ Do not re-open any of these while doing the above.

* `_project_regions_to_loci` — the node projection and its overlap math.
* `ρ_ref` / `_global_reference_density` — the shared enriched-mode reference.
* Per-object `min(m/ρ_ref, S)` contraction, and the ruling that it is applied **per object**, never over
  a folded total (that would under-contract a captured exon whose line runs into a depleted intron).
* The contained-evidence shrinkage `w = C/(C+1)` — one pseudo-observation, no tunable.
* The eff-len clamp into `[floor, span]`.
* Withholding spliced mass from `rna_prior_count`, and excluding the junction flux (owner, 2026-07-30).
* The accumulator schema. **`payload_schema_digest` must not move.**
* `q = edge_mass_per_crossing` as the incidence→fragment converter, and the ruling that it is geometry
  and must never join `prior_vs_oracle.OVERRIDE_FIELDS`.

---

## 7. GATES — write them first, then perturb each

1. **Edge-set identity.** The new edge→locus set equals the current `edge_owner_nodes`-derived set on
   every flgap and ladder condition. ⭐ This is the refactor gate; a difference is a finding.
2. **Uniqueness.** `Σ_L share(e,L) ≤ 1 + eps` for every edge, count reported.
3. **Outer lines are present.** A locus of `k` nodes carries `k+1` edges — falsified by dropping the
   first or last and watching the count fall.
4. **Conservation.** `Σ_L (gdna + rna) prior` equals the projected conserved mass over all
   locus-touching objects, to the fragment.
5. **Factor 1 under uniform gDNA** — `gdna_eff_len == span` exactly. Unchanged, must stay green.
6. **`q ≠ 1` discrimination** — already in `test_priors.py`; must stay green.
7. **Perturbations that must each fire**: drop the outer edges; use `lo` always (no re-key, old bug);
   fold edge mass into nodes (the current behaviour); rescale the contained term by `q`.

---

## 8. OPEN

1. ⚠ **Does `S` become the *right* number once `F` is corrected?** The mass at a locus's edges sums to
   ≤ 1 per straddling fragment — the deficit is mass on lines wholly outside the locus, which is dropped.
   For a fragment whose excursion crosses no further line the sum is exactly 1; for a deep excursion it is
   less. **Unmeasured.** An `F_overlap` arm (§5) would settle it.
2. ⚠ The `+3.01 %` and the `1.002 / 1.000` geometric prediction were measured against `F`. They need
   re-reading against `F_overlap` before either is quoted as an error.
3. ⚠ `prior_vs_oracle`'s `S` arm feeds the truth per-component share, which is not available in
   production. It remains a diagnostic ceiling, not a candidate implementation.

---

## 9. HOW THE FOLD GOT THERE — so it does not come back

* **Shipped `v0.7.1` stored the split and then undid it.** Its accumulator carried `boundary_mass_left`
  and `boundary_mass_right` as separate banks, and `derive.py` described a region as *"its contained mass
  and its TWO BOUNDARY SIDES"*. ⭐ Its **RNA** prior used them directly — `mass_rna_contained +
  mass_rna_left + mass_rna_right − mass_rna_spliced`, no pooling, no owner, no re-key. ⛔ Its **gDNA**
  prior sent them through `_pooled_seam_arrays`, which re-joined the two sides into one seam and then had
  to choose an owner — reintroducing the problem the split had already solved, and growing the
  `rekey_right` patch to cover the intergenic case.
* **Why shipped pooled**: its stated reason is the *contraction* — *"the two halves are one physical
  crossing event."* That is right for `min(m/ρ_ref, S)` and wrong for the count. One object was used for
  both. **The fix is to separate the two uses**, which §4 does.
* **The current tree inherited the gDNA pattern and extended it to RNA**, because the accumulator no
  longer stores sides — `edge_unspliced_mass` is one number per line, and there was nothing else it could
  do. That removal was correct on its own terms (`DESIGN.md` §3.1: the sum-then-halve shape hid an exact
  factor of 2 for months) and it is not being reversed: **§3 needs no sides**, because a line is never
  split between nodes — it is assigned whole to the loci that touch it.

⭐ **The lesson, once this lands** (candidate for `TRAPS.md`): *a quantity defined on one axis, folded
onto another so a downstream consumer's interface can reach it, will grow a heuristic to repair the
fold — and the heuristic will outlive everyone's memory of why the fold was there.* The consumer's
interface was `_project_regions_to_loci`, which divides by `region_size_bp` and so structurally cannot
see a 0-bp object.
