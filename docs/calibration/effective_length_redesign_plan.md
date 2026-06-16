# Boundary-aware capture effective-length: diagnosis, redesign, and staged plan

**Status:** design / decision input — 2026-06-15. Companion to `new_eff_len_plan.md` (the author's
proposal). Targets the dominant capture-driven gDNA→RNA error on the single production path (the iterative
simplex sweep, post Phase-4 teardown).

**Provenance.** Multi-agent design pass (verify → design → adversarial critique → synthesize) followed by
**lead verification**: the diagnosis dissection was re-run independently and reproduced exactly, and the two
load-bearing new claims were checked against the code/data — one confirmed, one corrected (see §0). Every
`file:line` is current `main`; quantities are from the freshly re-quanted `quick_3to1_5mb` suite.

---

## 0. Lead verification & corrections (read first)

**Diagnosis independently reproduced.** Re-running the transport dissection on
`gdna_gdna300_ss_0.99_nrna_rnd_capture_on` (432 mature↔nascent twins) gave, identically to the design pass:
median `|log(nascent/mature eff)|` **1.489 (transport OFF) → 1.388 (ON)**; transport contracts nascent eff-len
~10× harder than mature (median Δeff **−714 vs −71**); near-collapsed twins (eff-len advantage < 0.5 nat)
**35 → 107 / 432**; eff-lens bit-identical to the shipped `em_effective_length`. **The author's mechanism is
real**: transport shrinks the nascent-vs-mature gap by concentrating gDNA into fewer high-density regions.

**Reframing CONFIRMED (with exact numbers): the dominant nascent sink is gDNA→nascent, not mature↔nascent.**
3-pool net flow at the flagship-with-nascent condition (`net_flow_3pool_per_condition.tsv`):

| quantity | value |
|---|---|
| `gdna_net_surplus` | **−233,859** (gDNA leaking out) |
| `net_gdna_to_nrna` | **+125,907** |
| `net_gdna_to_mrna` | **+107,952** |
| `net_nrna_to_mrna` | **−9,361** (mature↔nascent net flux is tiny) |
| `nrna_net_surplus` | +135,268 (true 190,588 → assigned 325,856) |

The nascent surplus is ~93% gDNA-sourced; the gDNA deficit (−233,859) is entirely the sum of the two
gDNA→RNA fluxes. The zero-gDNA control (`gdna_none…rnd`) shows nascent *under*-assigned (−35,845) — so the
nascent sink **appears only when gDNA is present**. This refines the earlier "mature↔nascent is the dominant
face" reading: **the eff-len lever must target the gDNA component first.** (The per-transcript isoform-swap
metric still reads 68.5% because nascent twins *absorbing* gDNA-sourced mass register as RNA reallocation;
the 3-pool decomposition reveals the source is gDNA.)

**Correction to the proposed first move (PR-1).** The design pass claimed the gDNA-component length is
"FL-subtracted → 0.1 for a 142 bp exon (1400× too small)." That conflates two quantities. `region_eff_length`
*is* ~0.7 for a 142 bp region (FL-subtracted) — but the gDNA component uses `gdna_geom_len = region_eff_len +
left_side + right_side` ([derive.py:55](../../src/rigel/calibration/derive.py)), which for that region is
**283** (≈ 2× the raw 142), confirmed numerically:

| region L | `region_eff_length` (FL-sub) | side `E[min(ℓ,L)]` | `gdna_geom_len` = r+2s | raw L |
|---|---|---|---|---|
| 142 (small exon) | 0.7 | 141.3 | **283.3** | 142 |
| 500 | 195.4 | 304.6 | 804.6 | 500 |
| 9119 (intron) | 8813.5 | 305.5 | 9424.5 | 9119 |

So the gDNA component is **not** FL-subtracted to ~0 — it is *inflated* on small exons by the boundary
crossings (and double-counts shared boundaries across adjacent regions, since each region adds both its own
sides). The author's plan (gDNA = raw locus span, no FL subtraction) is still the right direction — it would
*shrink* small-exon gDNA length and remove the boundary double-count — but the expected effect is the
opposite sign of the synthesis's reasoning, and its magnitude on the leak **must be measured, not assumed**.
This does not change the plan's structure, only PR-1's rationale and the need to gate it on the net-flow.

---

## 1. Diagnosis

`_transport_boundary_flux` ([priors.py:103-159](../../src/rigel/calibration/priors.py)) forms
`g = contained + left + right` and re-attributes each internal boundary's pooled mass onto the two flanking
**regions** ∝ `density_ratio·boundary_capacity`, iterated. That transported `gdna_region` drives **both** the
gDNA-component IPR (`priors.assemble_priors`) **and** the per-transcript capture contraction
([capture_eff_length.py:127](../../src/rigel/calibration/capture_eff_length.py)) — the same function called
twice. The EM competes on `theta_c = count_c / eff_len_c` ([em_solver.cpp:785,1299]), so eff-len is the only
lever; collapsing the nascent/mature eff-len gap hands the split to the prior + per-fragment likelihood.

**Mechanism, precisely.** Per-region tracing shows transport moves mass from the **depleted side to the
enriched side within each boundary** — sometimes intron→exon (one top twin: Δexon +13,700, Δintron −12,188),
sometimes into intronic capture *peaks* (Δexon −126, Δintron +121). Net effect is identical: fewer
high-density regions ⇒ harder IPR contraction of the nascent footprint. So the redesign should target
**footprint concentration**, not literally "intron→exon" (a refinement of the author's framing).

**Magnitude honesty.** The transport effect on the nascent/mature gap is real but **modest** at ss0.99
(1.49→1.39 nats; nascent still ~4× mature with transport ON). Per §0, the binding error is gDNA→{nascent,
mature} misclassification — the gDNA-component eff-len — which transport's mass concentration also drives.
Capture is the whole problem: capture-off has essentially no gDNA leak (`gdna_net_surplus` +2,838 →
−233,859 capture on).

---

## 2. The redesign

> **⚠ Superseded by §7.** §2.2–§2.5 below describe an earlier *split-side* node model that a real-pipeline
> adversarial test **falsified** (splitting each seam's mass into two side nodes nearly doubles the leak, by
> Cauchy-Schwarz halving the IPR support). The authoritative design is **§7 (pooled boundary seams)**, which
> recovers the leak below the old transport baseline (5.5% vs 7.97%) conservation-correctly. §2.1 (geometric
> lengths) and §2.4 (the transcript→node map rules) carry forward; read §7 for the node model + the PR plan.

### 2.1 Geometric (capture-independent) length per component

| Component | Geometric eff-len | Source | Status |
|---|---|---|---|
| Mature mRNA | `exonic_L − RNA-FL-integral(clip)` | `compute_all_transcript_eff_lens(t_df.length)` (exonic L) | **exists, reuse** |
| Nascent nRNA | `genomic_span − RNA-FL-integral(clip)` | same call; nRNA `length = end−start` | **exists, reuse** (mature/nascent already diverge here) |
| gDNA (per locus) | **raw locus span, no FL subtraction** | `Σ_blk (end−start)` from `MultiLocus` | **NEW** (today: `Σ gdna_geom_len`, region+sides — see §0) |

Mature and nascent geometric lengths **already diverge** via the existing FL-marginal code (nascent uses the
full genomic span). The gDNA change (raw locus span, replacing the region+sides sum that over-counts small
exons and double-counts boundaries) is the isolated PR-1 lever.

### 2.2 Node model — boundaries as first-class mass owners

Node set `𝒩 = regions ∪ boundary-left-sides ∪ boundary-right-sides`; each carries `(gDNA mass, geometric
length)`, all already on `CalibrationResult`:

| Node | gDNA mass | length | source |
|---|---|---|---|
| region-contained | `mass_gdna_contained[r]` | `region_eff_length(L_r, gdna_fl) = E[max(0,L−ℓ)]` | result.py:42 |
| boundary left-side | `mass_gdna_left[r]` | `boundary_side_eff_length(gdna_fl, L_r) = E[min(ℓ,L_r)]` | result.py:44 |
| boundary right-side | `mass_gdna_right[r]` | `boundary_side_eff_length(gdna_fl, L_r)` | result.py:46 |

**Correction to the author's plan:** the per-side length is `E[min(ℓ, L_side)]` (→ `fl_mean` for a long
region, → `L_side` for a short one), which is exactly `gdna_boundary_len` on the result today — **not**
`frag_len − 1` (a two-sided quantity). No new geometry; we just stop transporting it.

**Conservation (preserved).** A fragment lands in exactly one node (contained XOR the two crossing sides);
`Σ contained + Σ sides = n_fragments`. The redesign **never re-allocates mass between nodes** — it changes
only which length each mass divides by and which nodes a component sums. The one thing that moved mass —
`_transport_boundary_flux` — is deleted. The per-node `gdna+rna=total` result invariant and the C++
accumulator ledger are untouched.

### 2.3 The IPR / enrichment contraction (same closed form, wider node set)

Keep the Laplace-smoothed IPR bit-for-bit, changing only the node set and removing transport. For component
`c` over node set `𝒩_c`:

```
G_c     = Σ g_n ;  supp_c = Σ g_n²/L_n ;  span_c = Σ L_n
eff_ipr_c     = min( (G_c+1)² / (supp_c + (2G_c+1)/span_c), span_c )
contraction_c = eff_ipr_c / span_c              # ∈ (0,1]
w_c           = G_c / (G_c + |𝒩_c|)             # evidence shrinkage (no tunable constant)
eff_em_c      = min( L_geom_c · (w_c·contraction_c + (1−w_c)), L_geom_c )
```

Per-node enrichment `E_n = (g_n/L_n)/ρ_global` is exposed for QC; the contraction uses the mass-unit IPR so
the capture-off identity holds (§3).

### 2.4 The transcript→node map (the careful, mostly-NEW part)

Replaces `_exon_region_incidence`. Emits `(inc_t, node_id, weight)` over a unified node space
(`node_id < R` = region, `≥ R` = boundary side):

- **Region membership — reuse the off-by-one-fixed bound** (`lo = searchsorted(ends, a, "right")`,
  `hi = searchsorted(starts, b, "left")`). mRNA → exon intervals; nRNA → full span. Stays **overlap-driven,
  never strand-gated** (Risk R3).
- **Boundary-side inclusion (NEW):** include sides touching the component's included regions; the leftmost
  region drops its outer (left) side and the rightmost drops its outer (right) side ⇒ a single-region
  transcript gets **zero** side nodes.
- **Exon-side-only at exon↔intron seams (NEW):** include a seam side iff the touched region is exon **for the
  transcript's strand** — test `signature & BIT_EXON_{strand}` (the bit, not the coarse "exon-wins" type;
  mixed-signature regions exist, e.g. `0xA`).
- **Start/end partial-overlap weight (NEW):** weight a region node by `overlap/region_len` (applied to both
  `g` and `L`, so `g/L` is invariant); side nodes weight 1.0. Material on ~23% of mappings.

**Strand-encoding landmine (must guard):** `t_df.strand ∈ {1,2}`, not `{+1,−1}`. A naive `strand>0 →
BIT_EXON_POS` empties 42% of nascent footprints (silent factor=1, no contraction). **Mandatory assertion:
every transcript maps to ≥1 node.**

### 2.5 Global density & iteration

`gdna_density_global` (the enrichment denominator) is already computed over all nodes and consumes the
converged sweep — no new outer loop. Keep the simplex-prior `rho_global` region-observable-only (uncaptured
baseline) to avoid locking a captured phantom; captured boundary seeds enter only via the var~mean and the
enrichment denominator (a read-only post-hoc denominator — assert no re-solve).

### 2.6 Files: NEW vs reused

Reuse: FL primitives (`effective_length.py:40-94`), node masses (`result.py:42-46`), the Laplace IPR closed
form + `w=G/(G+n)` + `min(fl·factor,fl)`, `_project_regions_to_loci`, adjacency, signature bits, the
off-by-one bound, the mature/nascent FL-marginal lengths. **NEW:** `region_eff_len` surfaced on the result;
gDNA geom = raw locus span; per-node `Σ g²/L` over {regions ∪ sides}; the `(inc_t, node_id, weight)` map +
rules + strand guard. **DELETE:** `_transport_boundary_flux` + both call sites. No C++ change.

---

## 3. Why it works — and where it is honestly inert

**Mature/nascent divergence (restored signal).** Mature footprint = exon regions (enriched) + exon↔exon
sides → contracts modestly. Nascent footprint = exon regions + **intron regions** (depleted: large `L`, tiny
`g`) + the intron↔exon **exon-side** nodes → the depleted intron nodes inflate `span_c` without inflating
`G_c`/`supp_c` → stronger contraction toward the enriched footprint, but the kept enriched sides stop it
collapsing all the way to mature. The intron regions + intron↔exon sides are nascent-only — the signal
transport emptied.

**Why introns don't re-leak (the author's central worry, structurally answered).** Naively deleting transport
would leave crossing mass as contained-like mass on the *intron region at its full length* → near-uniform
`g/L` → introns stop contracting → leak. First-class boundaries **quarantine** that crossing mass on the side
node at `E[min(ℓ,L)] ≈ fl_mean` (~300 bp, not the intron's thousands) → the intron region keeps only its
sparse genuinely-contained gDNA at its long length → stays depleted → still contracts.

**Preserved invariants.** capture-off → factor 1 (uniform density ⇒ `eff_ipr=span`); sparse → `w→0`, no
contraction; `G=0` → geometric; conservation (no mass moved). *Caveat (R1):* capture-off identity is
**near**-identical for the transcript path, not bit-identical, and the gDNA raw-span change moves every
gDNA-bearing golden — use an explicit tolerance and `--update-golden`.

**Where it is inert (honest).** Zero-gDNA libraries (`G≈0`, factor=1) — the redesign does nothing; that
phantom is prior/EM-decided. Annotated single-exon `is_nrna` rows (no intron structure) — structurally inert.
The redesign's reach is the **gDNA-present, multi-exon-nascent** sink.

---

## 4. Risks

| # | Risk | Severity | Mitigation |
|---|---|---|---|
| **R2** | Removing transport may **inflate** gDNA-component eff-len on high-gDNA loci (transport does genuine concentration for the gDNA component) → could **worsen** the gDNA→RNA leak. Note §0: PR-1 (raw span) pushes the *other* way on small exons — net effect is uncertain and **must be measured**. | SERIOUS | **Decouple the two consumers.** Validate the gDNA-component variant against `net_gdna_to_{nrna,mrna}` directly, per-locus. Flag any locus with NEW/OLD gDNA eff-len > 1.5×. |
| R1 | "capture-off bit-identical" is false (gDNA raw-span moves goldens; transcript path only near-identical). | SERIOUS | Drop "bit-identical"; explicit tolerance; land gDNA raw-span as its own `--update-golden` commit so a benchmark move bisects cleanly. |
| R3 | Strand-bit *region* gating would re-drop regions (the 20d43f2 bug in disguise); `signature=0` fixtures silently disable everything. | SERIOUS | Region inclusion stays overlap-driven; strand bit only for the seam-side decision; guard all-`signature==0`; decode `{1,2}`. |
| R4 | `test_incidence_len_t_ge_exonic` breaks under partial-overlap weights. | SERIOUS | Re-state as `Σ weight·overlap == exonic ± fl` + an explicit "no region dropped" test. |
| R5 | An adversarial impl found NEW ≈ OLD on the nascent/mature ratio (intron region-nodes dilute the boundary exon-sides in the pooled IPR). The ratio may not be the right success metric (§0: it's second-order). | SERIOUS (open) | Decide empirically after the gDNA-component fix; if residual, give the boundary-exon-side enrichment its own term (PR-4). |
| R6 | `n_reg → n_node` silently re-tunes the "no-constant" shrinkage. | MANAGEABLE | Make the pseudocount explicit; verify net contraction still increases. |
| R7 | Every `CalibrationResult(...)` site needs the new `region_eff_len` kwarg. | MANAGEABLE | grep all sites (~5 test constructors + harness) before landing. |
| R8/R9 | encompassed 50/50 + side IPR non-linearity; C++ accumulator. | UNFOUNDED | Document; no C++ touched. |

---

## 5. Implementation plan (ordered, PR-sized, net-flow-gated)

Gate every step: `pytest tests/` (37 s) green **+** the net-flow benchmark with the **3-pool
`net_gdna_to_{nrna,mrna}` / `gdna_net_surplus` as the primary metric** (not twin eff-len ratios), at minimum
`gdna300 × ss{0.99,0.50} × cap_on × {nrna_none, nrna_rnd}`.

- **PR-0 — instrument the 3-pool net-flow** per condition in the benchmark output. Record baseline
  (`ss0.99 cap_on rnd`: gdna_net −233,859, g→n +125,907, g→m +107,952, nrna surplus +135,268).
- **PR-1 — gDNA component length = raw locus span** (replace `Σ gdna_geom_len` with the `MultiLocus` span;
  no FL subtraction, no boundary double-count). Isolated; touches no node structure. **Highest-ROI, lowest-
  risk — measure it ALONE** (its sign on the leak is not assumable, per §0). Regen gDNA-bearing goldens after
  validating direction. STOP and reassess if it recovers most of the −233K deficit.
- **PR-2 — delete `_transport_boundary_flux`; untransported per-node IPR over {regions ∪ sides}.** Gate is
  **R2**: if gDNA eff-len inflates on the high-gDNA tail and `net_gdna_to_*` worsens, do not proceed — keep
  transport-OFF for the transcript path only and retain a concentration-preserving gDNA-component formula.
- **PR-3 — the transcript→node map** (region membership + side inclusion + rules b/c/d + strand `{1,2}`
  guard). The most bug-prone surface — mandatory regression tests: single-exon → 0 sides; multi-exon →
  internal seam sides, outer excluded; exon↔intron → exon-side only (strand bit); mixed-signature `0xA`;
  partial overlap → `g/L` invariant; **every transcript ≥1 node**; −strand mirror; re-stated off-by-one
  guard.
- **PR-4 (conditional) — boundary-exon-side enrichment as its own discriminator term**, only if PR-1..3 leave
  residual gDNA→nascent flux that bisects to intron-density inflation (R5).

Measure at every gate: (a) `gdna_net_surplus`, `net_gdna_to_{nrna,mrna}` (primary); (b) nascent assigned vs
true 190,588 (recovery, not just suppression); (c) per-locus gDNA eff-len NEW/OLD tail (R2); (d) capture-off
transcript-path near-identity within tolerance.

---

## 6. Scope call

**Minimum viable = PR-1, measured.** PR-1 is one isolated change directly attacking the dominant gDNA→{n,m}
flux; its effect is uncertain in sign (§0) so it must be measured first and alone. If it recovers most of the
deficit, the bug-prone PR-3 node-map may be unnecessary.

**Full redesign (PR-2 + PR-3) only if PR-1 leaves residual** that bisects to transport concentration /
intron-density inflation. PR-3 re-litigates the leftmost/rightmost off-by-one history; the adversarial impl
showed it does **not** by itself improve the nascent/mature ratio — consistent with §0 (the sink is
gDNA-sourced, so the gDNA-component fix, not the nascent/mature gap, is the lever).

**Order vs the mature↔nascent EM prior:** the two are largely **orthogonal** here (`net_nrna_to_mrna` is only
−9,361). The eff-len fix targets the gDNA→{n,m} leak (capture-on conditions, the suite's definition); the EM
prior targets the `G≈0` / single-exon phantom that eff-len is structurally inert on. **Ship PR-1 first
regardless** — it is the dominant, independent lever; the prior fix can proceed in parallel.

---

# 7. DEFINITIVE first-class-boundary-node design + plan (supersedes §2.2–§5's PR sketch)

**Status:** decided — 2026-06-15. This section is the authoritative spec for the author to implement. It
supersedes the §2.2/§2.3 *split-side* node model and the §5 PR-1/PR-2/PR-3 ordering above, both of which were
falsified by a real-pipeline adversarial test (the *split* form nearly **doubles** the leak). It keeps §0–§1
(the diagnosis), §2.1 (geometric lengths), §2.4 (the transcript→node map rules), §3, §4, §6. Every `file:line`
is current `main@27c59cd`; numbers are from the freshly re-quanted `quick_3to1_5mb` flagship
(`gdna_gdna300_ss_0.99_nrna_rnd_capture_on`), measured end-to-end through the real `quant_from_buffer`.

## 7.1 What PR-1 taught us (the measured regression + why first-class nodes recover it)

PR-1 made the gDNA-component length conservation-correct (`gdna_region = contained+left+right` over bare
`region_size_bp`, dropping `_transport_boundary_flux` and the length-non-conserving `gdna_geom_len`), and in
doing so **regressed the gDNA→RNA leak from 7.97% (transport) to 9.00%** (measured, `/tmp/dissect_run.py`
variant A vs B; nRNA sink 231,049 → 269,524). The mechanism: the old transport **consolidated** each
intron↔exon seam's two pooled halves onto one high-density node, producing a strong inverse-participation-ratio
(IPR) contraction of the gDNA footprint (small `gdna_eff_len` ⇒ gDNA competes hard ⇒ little leak); PR-1's
conservation-correct split returns the seam mass to its two flanking regions, where the intron-side half is
diluted over the whole (up to ~50 kb) intron — the captured boundary signal is lost, the footprint
under-contracts, and gDNA leaks into nascent/mature. First-class boundary nodes recover the concentration
**without** transport's two sins (it moved mass between genomic positions, and it divided by the
length-inflating `gdna_geom_len` whose `Σ = 1.065×` the genomic span, 2.0× on small exons): a boundary is given
its own mass and its own **effective** length `≈ fl_mean`, so the captured crossing mass is *quarantined at
FL-scale density* (mass / ~300 bp) instead of diluted over the intron, while its **genomic** length stays ≈ 0
so the conserved span is untouched.

## 7.2 The node model — regions + **pooled boundary seams** (CRITICAL correction to §2.2)

> **The decisive correction.** §2.2 / `new_eff_len_plan.md` modelled each boundary as **two** per-region side
> nodes (`mass_gdna_left[r]` at one node, `mass_gdna_right[r]` at another). A real-pipeline test
> (`/tmp/dissect_run.py`, `/tmp/dissect_split.py`) proved this **catastrophic**: it makes the flagship leak
> **13.52%** — worse than both PR-1 (9.00%) and transport (7.97%) — because splitting a seam's mass `s` into
> halves `s_L, s_R` and entering them as two nodes contributes support `s_L²/ε + s_R²/ε ≈ 2·(s/2)²/ε = s²/(2ε)`,
> whereas the physical crossing is one event and **pooling** it contributes `(s_L+s_R)²/ε = s²/ε`. By
> Cauchy-Schwarz the split **halves the IPR support** — i.e. halves the contraction — by construction.
> **The first-class node is the SEAM, not the two sides.** Pool the two halves into one node per internal
> boundary. Measured: this recovers the leak to **5.50%** (below the pre-PR-1 7.97%) and lands nRNA on truth
> (190,323 assigned vs true 190,588), conservation-correct (`/tmp/dissect_run.py` variant D).

The node set `𝒩` for a reference is the disjoint union of two families:

| Node | gDNA mass | RNA mass | genomic length `Lg` | effective length `Le` | source (verified) |
|---|---|---|---|---|---|
| **region** `r` (one per region) | `mass_gdna_contained[r]` | `mass_rna_contained[r]` | `region_size_bp[r]` | `region_size_bp[r]` *(genomic bp; same as PR-1 — region density `g/L` is per-genomic-bp)* | `result.py:42-43`; `region_arrays.py:55,109` |
| **seam** `(r,r+1)` (one per **internal** boundary, `ref_id[r]==ref_id[r+1]`) | `mass_gdna_right[r] + mass_gdna_left[r+1]` (the **pooled** crossing mass — both halves of the one seam) | `mass_rna_right[r] + mass_rna_left[r+1]` | `0` (≈ 1 bp crossing position) | `gdna_boundary_len[r]` = `boundary_side_eff_length(gdna_fl, region_size_bp)[r]` = `E[min(ℓ,L)]` ≈ `fl_mean` | `result.py:44-47,66`; `effective_length.py:71-94`; `calibrate.py:91,313` |

There are exactly `R + n_internal_seams` nodes (`n_internal_seams = Σ_f (k_f − 1)` = `R − n_refs`), **not**
`R + 2R`. The seam is keyed to its **left-flank region `r`** (mass deposited at index `r` of a length-`R`
array, masked by `same_ref_left_right(ref_id)` / `ref_id[:-1]==ref_id[1:]`). Region edges coincide with
exon edges on real data, so the two halves `mass_gdna_right[r]` and `mass_gdna_left[r+1]` are the substrate's
own seam ownership (`substrate.py:113-114`: `right[r]` and `left[r+1]` are the two sides of boundary
`(r,r+1)`), pooled at the seam — exactly the substrate bijection, attributed to one node.

**Seam effective length:** use `gdna_boundary_len[r]` (the left flank's `E[min(ℓ,L_r)]`). On large introns this
→ `fl_mean` (≈ 300; the FL-scale quarantine the recovery needs); on a short exon flank it → the short length.
The two flanks of a seam share the same `R` only when both regions are the same size; in the measured fix we
took `min(gdna_boundary_len[r], gdna_boundary_len[r+1])` (the tighter cap) — either flank's value is acceptable
since both ≈ `fl_mean` on the dominant intron↔exon seams; pin to **the left flank's** `gdna_boundary_len[r]`
for determinism and document it.

### Conservation, stated and PROVEN on this node set

**Law (M) — mass.** `Σ_{n∈𝒩} g_n = Σ_r mass_gdna_contained[r] + Σ_{internal (r,r+1)} (mass_gdna_right[r] +
mass_gdna_left[r+1])`. The seam family visits each internal boundary's two sides exactly once
(`right[r]` for its left half, `left[r+1]` for its right half — distinct array entries; the substrate `left`/
`right` views partition every internal side bijectively, `substrate.py:105-114`, `boundary_region_indices`
`region_arrays.py:156-185`, terminal sides carry zero mass). So `Σ g_n = Σ contained + Σ left + Σ right` = the
total deconvolved gDNA, identical to PR-1's folded total (verified `isclose=True`, `/tmp/adv_conservation.py`)
and to the accumulator ledger `Σ contained + Σ sides` (`_accumulator_reference.py:254-263`). RNA identical via
`mass_rna_*`. **No mass moved (no transport), no double-count** (region nodes carry `contained` only; the
`+left+right` fold is removed and replaced by the pooled seam node). ∎

> **Honest restatement (Obj 2):** `Σ contained + Σ sides = total *deposited* mass`, **not** `= n_fragments`.
> On real data 0.2% of fragments are unresolved/outside all regions (3,992,068 vs 4,000,000). The conservation
> the design needs is *internal* — node sum = total deconvolved mass — which holds exactly. Strike "= n_fragments"
> wherever it appears (it is only exact on the toy 4-fragment accumulator spec).

**Law (L) — genomic length.** `Σ_{n∈𝒩} Lg_n = Σ_r region_size_bp[r] + Σ_seams 0 = Σ_r region_size_bp[r]` =
the covered genomic span (regions tile: verified `Σ = 5,000,000`, zero gaps/overlaps). Seam nodes contribute
**0** genomic bp — a boundary is a single crossing position; its `gdna_boundary_len ≈ fl_mean` is a *density
normalizer*, used only inside the dimensionless IPR factor, **never** added to the genomic span. This is exactly
the conservation that `gdna_geom_len = region_eff + 2·side_eff` violated (`Σ = 1.065×` span, 2.0× small exons).
∎

## 7.3 The gDNA effective length — geometric span × enrichment contraction

For a locus `ℓ` with overlap-weighted region set, the gDNA node set `𝒩_ℓ = {region nodes of ℓ} ∪ {internal
seam nodes within ℓ}`. Accumulate (the seam terms are pre-divided per region *before* the existing
overlap-weighted `_project_regions_to_loci`, so no new projection machinery — the seam is keyed to its
left-flank region `r` and rides that region's locus share):

```
G_ℓ      = Σ_{r}  share·g_region[r]   +  Σ_{seams}  share·g_seam            # total locus gDNA (region + pooled seam)
span_ℓ   = Σ_{r}  share·region_size_bp[r]                                   # CONSERVED genomic span (seams add 0)
supp_ℓ   = Σ_{r}  share·g_region[r]²/region_size_bp[r]                      # region density participation
         +  Σ_{seams}  share·g_seam²/gdna_boundary_len[seam]               # seam density participation (at FL)
gdna_eff_len = min( (G_ℓ+1)² / [ supp_ℓ + (2·G_ℓ+1)/span_ℓ ],  span_ℓ )    # Laplace +1, capped at genomic span
```

This is the **identical Laplace-smoothed IPR closed form** as `priors.py:272-275` (reused bit-for-bit): the
only change is `g_region[r]` is **contained-only** (drop the `+left+right` fold at `priors.py:222-226`) and the
two seam terms are added into `G_ℓ` and `supp_ℓ`. The cap and the Laplace reference stay purely genomic
(`span_ℓ` = `Σ region_size_bp`), so all the existing invariants are preserved (below). It is the
**`geometric genomic span × dimensionless contraction`** structure the brief names: `gdna_eff_len = span_ℓ ×
factor`, `factor = (G+1)²/[span·supp + (2G+1)] ∈ (0,1]`.

### The leak-recovery proof (toy + measured)

The captured crossing mass `s` at an intron↔exon seam enters `supp` at `gdna_boundary_len ≈ fl_mean ≈ 300`:
`s²/300`. Under PR-1 the same `s` was folded into the **intron region's** contained mass and divided by the
intron's `region_size_bp` (up to 495,958 bp): `s²/495958` — **~1600× smaller** support. The pooled seam node's
large `supp` contribution shrinks `(G+1)²/supp` ⇒ **smaller** `gdna_eff_len` ⇒ higher gDNA per-position rate ⇒
gDNA competes harder ⇒ **less leak** — *without moving mass* (transport's sin) and *without inflating the span*
(`gdna_geom_len`'s sin). Measured on the dominant leaker (locus 21, 40% of all leak — a 25-region multi-exon
gene whose captured gDNA is almost entirely at the seams): transport `eff_len=3,718`, PR-1 `12,933`, **split
(C) 15,565** (worse — the Cauchy-Schwarz halving), **pooled seam (D) 2,396** (below transport). Aggregate:

| variant | gDNA eff-len recipe | flagship leak_frac | nRNA assigned (true 190,588) |
|---|---|---|---|
| A transport (pre-PR-1) | `_transport_boundary_flux` + `gdna_geom_len` divisor | 0.0797 | 231,049 |
| B PR-1 (committed) | `contained+left+right` over `region_size_bp`, region-only IPR | 0.0900 | 269,524 |
| **C split-sides (the §2.2 design — REJECTED)** | region-contained + **2** side nodes at `gdna_boundary_len` | **0.1352** | 375,747 |
| **D pooled-seam (THIS design)** | region-contained + **1 seam node** `right[r]+left[r+1]` at `gdna_boundary_len` | **0.0550** | **190,323** |

(`/tmp/dissect_run.py`, real `quant_from_buffer`, all four swap **only** the gDNA eff-len in `assemble_priors`.)
D is the only variant that recovers below transport while staying conservation-correct.

### Preserved invariants (verified)

- **capture-off / uniform gDNA ⇒ factor → 1 ⇒ `eff_len = span`** (the cap pins it): uniform `g_n = ρ·L_n` gives
  `(ρ·span+1)²/(ρ²·span·span/span + 2ρ·span+1) → 1`; seam nodes under uniform gDNA also sit at density `ρ`
  (`g_seam ≈ ρ·fl_mean`), entering numerator and the (capped) inner participation identically, so the cap pins
  uniform loci to `span`. **NOT bit-identical** (Obj 4): measured capture-off whole-genome factor is 0.998822
  for PR-1 (gDNA is never perfectly uniform), so capture-off **goldens move** — strike "bit-identical", use
  tolerance + `--update-golden`.
- **sparse-but-concentrated ⇒ ≈ span:** the Laplace `+1` dominates at `G ≪ 1` ⇒ factor ≈ 1; the EM cannot
  amplify thin gDNA past the calibration's call.
- **G = 0 ⇒ span:** `factor = 1/(span·0+1)·…→ 1`. Verified on the zero-gDNA control (`gdna_none`): PR-1 and the
  node forms are identical to 4 sig figs (no phantom contraction from seam nodes, `/tmp/adv_perlocus2.py`).
- **EM contract unchanged:** the node IPR collapses to the single per-locus scalar `gdna_eff_len`
  (`em_solver.cpp:1846-1854` consumes one `log_eff_len[gdna_idx]`; `theta_c ∝ count/eff_len`). **No C++ change.**

## 7.4 The transcript→node & gDNA→node maps

### gDNA→node (for `priors.assemble_priors`)

Per locus, via the existing genomic-overlap `_project_regions_to_loci` (`priors.py:40-100`). Region nodes =
the locus's overlapped regions (existing). Seam nodes = the **internal seams between those regions**, keyed to
the left-flank region `r` (so the seam rides region `r`'s locus share). Build two extra length-`R` arrays
*before* projection: `g_seam[r] = (ref_id[r]==ref_id[r+1]) ? mass_gdna_right[r]+mass_gdna_left[r+1] : 0`, and
`supp_seam[r] = g_seam[r]²/gdna_boundary_len[r]`; project `{g: contained+g_seam, supp: contained²/L+supp_seam,
span: region_size_bp, rna: …}`. The "exon-side-only" / outer-drop rules below do **not** apply to the gDNA
component — gDNA is unstranded and its footprint is the whole locus; those rules are transcript-only.

### transcript→node (for `capture_eff_length`)

Replace `_exon_region_incidence` → `_transcript_node_incidence(index, region_arrays, calibration)` emitting
`(inc_t, node_id, weight)` over a unified node space (`node_id < R` = region; `node_id ≥ R` = seam, keyed to
its left-flank region). Because the node is now the **seam** (not two sides), the "exon-side-only at exon↔intron
seams" rule (§2.4) **simplifies to "include the interior seams of the footprint"** — a seam is a single node;
there is no left/right side to choose. The rules:

1. **Region membership — reuse the off-by-one-fixed bound verbatim** (the 20d43f2 fix, `capture_eff_length.py:72-73`):
   `lo = lo0 + searchsorted(ends, a, "right")`, `hi = lo0 + searchsorted(starts, b, "left")`. mRNA → exon
   intervals (`intervals.feather`, `interval_type==EXON`); nRNA synthetic → full span (`t_df`, absent from
   intervals). **Region inclusion stays overlap-driven, NEVER strand-gated** (Risk R3 — strand-gating inclusion
   re-introduces the 20d43f2 drop; and on this index it is a structural no-op anyway: 0/2705 exon→region
   overlaps touch a non-exon-for-strand region, so a strand-gated filter is dead code that ships untested).
2. **Seam-node inclusion (NEW, simplified):** include the interior seam `(r,r+1)` iff **both** `r` and `r+1`
   are in the footprint. The **leftmost** included region contributes no left-outer seam, the **rightmost** no
   right-outer seam — i.e. a crossing overhanging past the transcript's first-exon-start / last-exon-end is
   excluded. ⇒ **single-region transcript → 0 seam nodes** (mandatory regression case). mature: footprint =
   exon regions ⇒ seams between adjacent kept exon-regions; nascent: footprint = all regions in `[lo,hi)` ⇒ all
   interior seams (incl. intron↔exon and intron↔intron) — the nascent-only seam signal.
3. **Partial-overlap weight (NEW):** weight a region node by `w = overlap/region_size_bp` applied to **both**
   `g` and `L` (so the per-node density `g/L` is invariant; only the contribution magnitude scales). Seam nodes
   weight `1.0` (a crossing is whole or excluded). **Pin the IPR weighting to ONE form across the codebase:**
   `G = Σ w·g`, `supp = Σ w·g²/L`, `span = Σ w·L` — this is the form that preserves the capture-off identity
   (`(Σwg)²/(Σwg²/L)=ΣwL` under `g=ρL`). Material on ~24% of mappings (603/2466 exons start/end interior to a
   region).
4. **Strand `{1,2}` guard + ≥1-node assertion (NEW, mandatory):** `t_df.strand ∈ {1,2}` = the `Strand` enum
   (POS=1, NEG=2), **not** signed `{+1,−1}` and **not** the int8 `strand_class` (TS_NEG=−1, TS_AMBIG=2). The
   strand bit is needed **only** if a future PR-C re-introduces an exon-side discriminator; in this design it is
   unused (region inclusion is overlap-driven, the seam is a single node). Keep the decode + the runtime
   assertion `len(nodes) ≥ 1` for every transcript regardless — a naive `strand>0 → POS` empties 42% of
   footprints (silent factor=1). The IPR then runs over the emitted region+seam rows with the §7.3 formula,
   capped at the transcript's geometric `fl` (`min(fl·factor, fl)`, `capture_eff_length.py:164`, reused).

#### The bug surface (each with its mandatory regression test)

| surface | hazard | test |
|---|---|---|
| region bound `[lo,hi)` | reusing starts for `lo` drops fully-contained exons (the 20d43f2 bug) | `test_no_region_dropped` (keep the existing guard) |
| outer-drop | including the first region's left seam / last region's right seam double-counts overhang | `test_single_region_zero_seams`; `test_multi_region_internal_seams_only` |
| partial overlap | unweighted inclusion double-counts edge gDNA owned by a neighbor; `len_t ≥ exonic` test breaks | `test_partial_overlap_g_over_l_invariant` (`(w·g)/(w·L)==g/L`; `Σ w·overlap == exonic` **exactly**, no `±FL` — Obj 4); restate the broken `test_incidence_len_t_ge_exonic` |
| strand `{1,2}` | `strand>0→POS` empties 42% of footprints | `test_strand_decode_one_two`; `test_every_transcript_ge_one_node` |
| mature/nascent divergence | nascent must include intron region-nodes + intron-bearing interior seams | `test_mature_nascent_node_sets_differ` (live: 265/265 nascent ⊋ mature) |
| encompassed 50/50 | interior region fully traversed credits both its sides; as a seam node the pooled mass is the same | `test_encompassed_region_seam` (3-region footprint, middle encompassed) |
| capture-off | factor ≈ 1, NOT bit-identical | `test_capture_off_factor_one` (within tolerance) |

> **Note (Obj 1/3 honesty):** the mature/nascent divergence is dominated by the **intron region-nodes**
> (depleted: large `L`, tiny `g` ⇒ inflate `span` without `G`/`supp` ⇒ stronger contraction), not by the seam
> nodes. The adversarial finding (R5) that NEW ≈ OLD on the *twin nascent/mature eff-len ratio* still stands —
> the transcript-map change is **second-order** on `net_nrna_to_mrna` (only −9,361). The dominant lever is the
> **gDNA-component** seam fix (§7.3); validate the transcript map on **`net_gdna_to_{nrna,mrna}` per condition**,
> not twin ratios.

## 7.5 Enrichment ratios + global density + iteration

- **Per-node enrichment** `E_n = (g_n / Le_n) / ρ̄_g` (dimensionless; exposed for QC only — the IPR consumes the
  mass-unit sums `G, Σg²/L` directly for the capture-off identity, never the explicit ratio). A captured seam
  has `g_seam` large, `Le ≈ fl_mean ≈ 300` ⇒ `E_seam ≫ 1` (the quarantine); a deep-intron region has `E ≈ 0`.
- **Per-transcript / per-gDNA enrichment** = the IPR `factor_c ∈ (0,1]` of §7.3/§7.4 (the contraction); 1 ⇒ no
  enrichment (uniform / sparse), small ⇒ concentrated capture footprint.
- **Global gDNA density `ρ̄_g`** = the seed-region floor: `Σ_{count-observable region nodes} g / Σ Le`
  (intergenic + intron regions = `region_count_observable`, `density_model.py:88`; `signature` non-exon). This is
  exactly the sweep's existing `rho_global` (`calibrate.py:244-246`), restricted to observable regions to avoid
  locking an AMBIG/RNA phantom (`calibrate.py:239-243`). The **captured seam seeds** (intron↔exon crossings) are
  *consumers* of `ρ̄_g` (their `E_seam ≫ 1` is the signal), **not** averaged into it — averaging them inflates the
  reference and suppresses every enrichment.
- **Iteration** is the existing per-pass sweep (`calibrate.py:238-287`): pass 0 all-gDNA init; each pass re-fits
  `ρ̄_g` + the gDNA `var~mean` (`MonotoneVarMean`, already fed the boundary sides `gdna_left/right`,
  `calibrate.py:253-256`) on the running gDNA estimate, re-solves the per-node simplex, converges on `f_g`
  (tol 1e-3, `sweep_max_passes=4`). **No new outer loop.** The enrichment ratio (`E_n`) and the var~mean
  precision (`τ_count`) are orthogonal: `τ_count` weighs how much to trust the count magnitude in the simplex
  solve; `E_n` (via the IPR) sets how concentrated the resulting gDNA footprint is for the EM eff-len. The seam
  node's first-class status changes only the enrichment-IPR consumers + the global-density node set, not the
  per-node simplex solve (the `gdna_density_global` denominator is a read-only post-hoc — assert no re-solve).

## 7.6 Risks (adversarial findings, each with mitigation)

| # | finding | severity | mitigation |
|---|---|---|---|
| **F1 (FATAL, resolved)** | The **split-side** node model (§2.2) makes the leak **worse** (13.52% vs 9.00% vs 7.97%) — Cauchy-Schwarz halves the IPR support. | was FATAL | **Use the POOLED seam node** (`right[r]+left[r+1]`, one per internal boundary). Measured 5.50% (below transport). This design is built on D, not C. |
| F2 | "Σ contained + Σ sides = n_fragments" is false on real data (0.2% unresolved). | SERIOUS (doc) | State "= total *deposited* mass"; the internal node-sum = deconvolved-total conservation holds exactly. |
| F3 | "capture-off bit-identical" is false (PR-1 factor 0.9988; goldens move). | SERIOUS | Strike "bit-identical"; explicit tolerance; land each consumer change as its own `--update-golden` commit so a benchmark move bisects. |
| F4 (R2) | Removing transport from the **transcript** path has uncertain sign (may inflate gDNA-adjacent transcript eff-len). | SERIOUS | **Decouple the two consumers into separate PRs** (§7.7). Gate the transcript PR on `net_gdna_to_{nrna,mrna}`; flag loci with NEW/OLD eff-len > 1.5×; do not proceed if it worsens. |
| F5 (R3) | Strand-gating *region inclusion* re-drops regions; the strand-bit apparatus is a structural no-op (0/2705) ⇒ ships untested. | SERIOUS | Region inclusion overlap-driven; the seam-as-single-node makes the exon-side rule moot ⇒ the strand bit is unused in this design. Keep only the `{1,2}` decode + ≥1-node assertion. If a future PR-C needs the bit, build a synthetic opposite-strand-exon fixture first. |
| F6 (R4) | Partial-overlap breaks `test_incidence_len_t_ge_exonic`; three docs specified three different `supp` weightings. | SERIOUS | Pin `G=Σw·g, supp=Σw·g²/L, span=Σw·L`; restate the test as `Σ w·overlap == exonic` exactly (no `±FL`); add the explicit no-drop guard. |
| F7 (R5) | NEW ≈ OLD on the twin nascent/mature eff-len **ratio**; the divergence is intron region-nodes (dilutive), not seam nodes. | SERIOUS (open) | The dominant lever is the gDNA-component seam fix (PR-A), not the twin ratio. Validate on `net_gdna_to_{n,m}`. PR-C (boundary-exon-side own term) only if residual flux bisects to intron-density inflation. |
| F8 (R6) | `n_reg → n_node` re-tunes the `w=G/(G+n)` shrinkage. | MANAGEABLE | Make the pseudocount explicit; verify net contraction still increases. |
| F9 (R7) | `region_eff_len` schema field touches every `CalibrationResult(...)` site. | MANAGEABLE | **Not needed by PR-A** (priors uses `region_size_bp` + `gdna_boundary_len`, both already on the result). Defer to PR-B only if the transcript path needs `region_eff_length` separately; grep all sites first. |
| F10 (R8/R9) | encompassed 50/50 + seam IPR non-linearity; C++ accumulator. | UNFOUNDED | The pooled seam re-unifies the two halves the accumulator split, so the 50/50 asymmetry is moot. Document; no C++ touched (accumulator spec untouched, verified). |

## 7.7 Implementation sequence (PR-sized, net-flow-gated)

Gate every step: `pytest tests/` green **+** the `calibration-benchmark` net-flow with the **3-pool
`net_gdna_to_{nrna,mrna}` / `gdna_net_surplus`** as the primary metric (NOT twin eff-len ratios), at minimum
`gdna300 × ss{0.99,0.50} × cap_on × {nrna_none, nrna_rnd}`. Plus the weak-SS (ss0.65) zero-gDNA-phantom
tolerance regression (the brief's tol-250 case — pooling is more aggressive than PR-1, confirm it does not
over-contract a zero-gDNA library).

- **PR-0 — instrument the 3-pool net-flow** per condition. Record baselines: flagship `ss0.99 cap_on rnd`
  leak_frac 0.0900 (committed B), nRNA 269,524; pre-PR-1 transport 0.0797 / 231,049.

- **PR-A — gDNA-component pooled-seam nodes in `priors.assemble_priors`** (the isolated recovery lever).
  - Files: `src/rigel/calibration/priors.py` (`:222-226` drop the `+left+right` fold → `g_region = mass_gdna_contained`;
    build `g_seam[r] = same ? right[r]+left[r+1] : 0`, `supp_seam[r] = g_seam²/gdna_boundary_len[r]`; project
    `{g: contained+g_seam, support: contained²/L + supp_seam, span: region_size_bp, rna}`; IPR closed form
    `:272-275` unchanged). `rna_region` mirrors with `mass_rna_*` (keep the `−mass_rna_spliced` withholding,
    `:233-239`).
  - Conservation tests: `test_gdna_node_mass_conserved` (`Σ g_region+g_seam == Σ contained+left+right`,
    isclose); `test_gdna_uniform_eff_len_eq_span` (uniform gDNA ⇒ `gdna_eff_len == locus span` to tolerance —
    the cap pins it, proves seam nodes don't break capture-off); `test_gdna_G0_eq_span`;
    `test_single_region_locus_reduces_to_PR1` (no internal seam ⇒ identical to region-only).
  - **Measure ALONE.** Success: flagship leak_frac drops from 0.0900 toward **≤ 0.0797** (ideally ~0.055 as
    measured), nRNA assigned toward the true 190,588, with no ss0.65 zero-gDNA over-contraction. Regen
    gDNA-bearing goldens as a dedicated `--update-golden` commit *after* the net-flow confirms direction.
  - **STOP and reassess** if PR-A recovers most of the −233K gDNA deficit (it should — D landed nRNA on truth);
    PR-B/PR-C may then be unnecessary or purely for the "all components conserved" aesthetic.

- **PR-B — `capture_eff_length` pooled-seam nodes (mature/nascent transcript eff-len).** Only after PR-A is
  measured.
  - Files: `src/rigel/calibration/capture_eff_length.py` — delete `from .priors import _transport_boundary_flux`
    (`:31`) + the transport call (`:127-134`) + the `geom = gdna_geom_len` divisor (`:135`); replace
    `_exon_region_incidence` (`:41-96`) with `_transcript_node_incidence` emitting `(inc_t, node_id, weight)`
    over {region nodes (`region_size_bp`) ∪ interior pooled-seam nodes (`gdna_boundary_len`)} per §7.4; keep the
    Laplace IPR + `w=G/(G+n_node)` + `min(fl·factor, fl)` (`:150-164`).
  - Tests: the full §7.4 bug-surface battery on a NEW fixture extending `_MISALIGNED_GTF` (add a multi-exon POS,
    a multi-exon NEG, a synthetic nRNA twin; the existing fixture has no multi-region-with-interior-intron and
    no nRNA twin — it cannot exercise any new rule).
  - Gate F4: if `net_gdna_to_*` worsens or the high-gDNA tail inflates (NEW/OLD eff-len > 1.5× on many loci),
    **do not proceed** — keep transport-OFF gDNA-only and retain the existing transcript path. After both
    consumers migrate, `_transport_boundary_flux` (`priors.py:103-159`) and `derive.gdna_geom_len`
    (`derive.py:55`) have zero consumers → delete; `gdna_boundary_len` is repurposed as the seam-node length.

- **PR-C (conditional) — boundary-exon-side as its own discriminator term**, only if PR-A+B leave residual
  gDNA→nascent flux that bisects to intron region-node density dilution (F7). This is where the strand `{1,2}` /
  signature-bit apparatus would finally be load-bearing — build the synthetic opposite-strand-exon fixture
  first (F5).

**What to measure / success criterion.** Primary: flagship `gdna300 ss0.99 cap_on rnd` leak_frac **recovers
from 0.0900 to ≤ 0.0797 (ideally ~0.055, measured)** with nRNA assigned ≈ true 190,588, **while staying
conservation-correct** (`Σ node gDNA mass = total deconvolved`, `Σ node genomic length = covered span`).
Secondary: `gdna_net_surplus`, `net_gdna_to_{nrna,mrna}` per condition (recovery of the −233,859 deficit, not
mere suppression); ss0.65 zero-gDNA phantom tolerance held; capture-off transcript path within tolerance.
