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
