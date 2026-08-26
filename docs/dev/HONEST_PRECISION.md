# HONEST PRECISION — the policy arbitration and the message-layer redesign frame (2026-08-25)

The owner's mandate: the transport-dispersion decomposition is the key next step; the approach
must be DERIVATION-based; determine how variance should be computed for message propagation —
what the message/claim types are, how precision is computed for each, how it should be handled —
and tear down the existing infrastructure into a simple, clean, elegant message-propagation
interface with honest precision. This doc is the frame; every number in it was measured on the
current tree (post reference-deletion, post anchor-as-message) unless marked otherwise.

## 1. The policy arbitration — relay vs currency, measured clean for the first time

Both prior verdicts about CurrencyPolicy were confounded: "worst on every in-scope contaminated
stratum" predates the reference deletion (82–95 % of scored unstranded pass-0 error sat on slots
the deleted reference decided) and predates the anchor. On the CURRENT tree, all 16 conditions,
whole-library |est − truth| in fragments, plus the missing cell — currency carrying the
certified-flux stream (delivered via the same `PsiMessage.lam_rows` contract, patched in-process):

| condition (in-scope) | silent (control) | currency | currency+flux | relay+flux (shipped) |
|---|---|---|---|---|
| g00-OFF unstranded | 683,696 | **7,892** | 10,788 | 20,711 |
| g00-OFF stranded | 26,673 | **1,730** | — | 5,800 |
| g00-ON stranded | 14,524 | **4,549** | 5,350 | 19,538 |
| g05-OFF unstranded | **28,197** | 34,893 | 36,530 | 44,657 |
| g05-OFF stranded | **26,992** | 32,396 | — | 39,731 |
| g05-ON stranded | **38,518** | 52,099 | — | 67,664 |
| g50-OFF unstranded | **88,053** | 109,191 | — | 118,489 |
| g50-OFF stranded | **75,829** | 122,530 | — | 111,006 |
| g50-ON stranded | **74,936** | 106,302 | 103,185 | 164,477 |
| g98-OFF unstranded | **87,018** | 219,321 | — | 135,506 |
| g98-OFF stranded | **66,011** | 279,553 | 276,590 | 135,835 |
| g98-ON stranded | **80,685** | 130,282 | 111,263 | 158,469 |
| deferred g00/g05/g50/g98 | 227k/206k/2,785k/5,458k | 3.7k/199k/1,261k/2,079k | —/—/**388k**/— | 61k/134k/910k/1,294k |

**The readings, in order of importance:**

1. ⛔ **Silent still beats every message policy on every contaminated in-scope row.** Both
   composition schemes are net-harmful where the local solve has evidence. Whatever ships, the
   bar is silence, and the recorded attribution stands: the damage concentrates on destinations
   that HAVE own evidence (96 % of currency's at `g50` ss.99 OFF) — a PRECISION defect, not a
   value defect. This is the single strongest argument for the honest-precision rebuild: the
   problem is not which policy, it is that neither prices its hops honestly.
2. **Head-to-head, currency's scheme now beats relay's on ~11 of 16 rows** — all four zero
   controls decisively, the stranded capture-ON rows, most clean rows. Relay's scheme wins
   high-gDNA capture-OFF (`g98`-OFF: 136k vs 277k) and the heavy deferred rows.
3. **Currency + flux beats relay + flux on 6 of 7 cells measured** — g00-ON control 5.4k vs
   19.5k, g50-ON 103k vs 164k, g98-ON 111k vs 158k, deferred g50 388k vs 910k (2.3×; and 3.2×
   over currency alone — the stream and the clean scheme COMPOUND). The flux stream's tax on
   currency's zero controls is small (7.9k → 10.8k) where relay's is large.
4. The freeze verdict on currency is OVERTURNED as a ranking; what stands is the per-mechanism
   record: its enrichment-ratio input carries the `1/P(w ≤ ell_region)` factor (11.6× at a 98 bp
   exon) and a 58.6 % zero-bank accidental mute, and REPAIRING the input made the deliverable
   WORSE — more evidence that delivered PRECISION, not input value, is the binding problem.

**The verdict (proposed): neither file is the production policy.** The production policy is the
rebuild of §3, on **currency's chassis principles** — belief-free frames, one shrinkage knob
instead of licences-on-a-default, a premise charged on every hop, non-amplifying by construction,
no mass rescale, one shared hop table for both scan directions — **carrying relay's channel
completeness** — the certified-flux stream, the λ/θ delivery coordinates with the emission gate,
the per-strand population licence (currency's single largest structural gap; the 2026-08-18
dissection repair), the three recorded reception laws (mismatch deflation as reception-only). The
harvest table:

| keep from CURRENCY | keep from RELAY | keep from the ANCHOR work |
|---|---|---|
| the knob `w = (log r)²/((log r)² + v)` — enrichment as ONE unknown, no switch | the certified-flux stream + `lam_rows` contract | route-sum pooling (fix relay's splice-in at the geometry source) |
| precision-weighted premise by moments, every hop, floor 0 | per-strand population licence (three-case rule) | estimator refusal discipline (guards; refuse ⇒ flat, never confident) |
| residual disagreement `((1−w)·log r)²` — knob and damping one statement | λ/θ explicit coordinates + emission gate | NB counting law for count claims (`size = flux + ½`) |
| sampling variance as an idempotent CAP | mismatch deflation (reception: lowers precision, never moves a mode) | the citizenship rule: imputations never enter own-evidence precision |
| belief-free face totals (kills the `belief_fg` frame debt, 57–77 % of mass) | the one-sided/fused splice-in bound record | max-combine only where derivation is missing (and it is scaffolding to be replaced by §4) |
| structural absence of mass rescale | the switch-per-operator PRICING method (as instrument, not architecture) | |

## 2. The claim taxonomy — what actually flows, and what precision each carries

Everything that informs a slot reduces to FOUR claim kinds. (The full formula-level inventory
with file:line for every law was compiled 2026-08-25 and lives in the session record; this is the
load-bearing shape.)

| kind | examples | precision law today | honest-precision rule |
|---|---|---|---|
| **OWN observation** | strand counts; own densities; intron-factory factor | counting (`ψ'(n+½)`) + curvature | counting only — no imputation may enter (`has_own_evidence` is the arbiter the whole layer keys on) |
| **CERTIFIED count claim** | spliced flux (route-summed); the flux stream's rows | NB size + fitted scatter (the 3-estimator max) | counting + the DERIVED per-hop-class premise of §4 — the estimator stack is scaffolding |
| **TRANSPORTED belief** | relay/currency hop claims (abundance or composition) | source precision − counting hop costs − assorted dampers (20 relay operators / 4 currency terms) | source precision, then EXACTLY three costs: counting (shrinks with depth), reframe `(log r)²`-family (geometry), premise (never shrinks) — then reception |
| **RECEPTION transforms** | mismatch deflation; conservation surprise; the λ-emission gate | DL forms; some read destination belief legally (reception may) | reception may only LOWER precision, never move a mode; biases may not be priced as variances (P1e must shrink as its strata are diagnosed) |

The consumer contract (ψ) stays as-is and is already clean: four Gaussian channels in declared
coordinates plus ONE row-array channel for non-Gaussian evidence — the redesign changes what is
SENT, not what ψ accepts.

**The invariants the backbone should ASSERT (several already hold; make all structural):**
1. A single-source transform only reduces a precision (`p → p/(1+p·v)`, v ≥ 0); rises are
   additive fusions of independent witnesses only. (Verified true of every shipped transform
   today — the fan-out died partly for violating it. Assert it, don't trust it.)
2. A refused claim loses its precision in the same statement as its value.
3. One counting law: `count_logvar = ψ'(n+½)` — three raw `1/n` stragglers survive in the relay
   (`v_mu` at zero-count = ∞ annihilation; two conservation terms) and must die.
4. Claims live in declared coordinates, in-domain (off-grid mode = pin at the boundary).
5. The premise cost appears on EVERY hop and does not shrink with counting depth.
6. Frames are belief-free (currency proves it possible; the relay's `belief_fg` frame debt
   covers 57–77 % of library mass).
7. Imputations never enter phase-A or own-evidence precision (the anchor's citizenship rule,
   generalized: it is what `has_own_evidence`, the refits, and the deflation all key on).

## 3. The teardown target

ONE hop law replaces the operator pile. A hop transports a claim as:

```
value:      v_dst = v_src · r^w · (population licence)        # currency's knob, per strand
precision:  p_dst = reduce(p_src; counting, reframe, premise) # the three costs, nothing else
reception:  p_dst = deflate(p_dst; own evidence)              # may only lower
```

What dies when this lands (each with its measured epitaph already recorded): the relay's 20-switch
operator pile (kept only as the PRICING instrument until parity), the P1d pooled scalar AND
currency's global premise scalar AND the anchor's three-estimator max — all three replaced by the
ONE per-hop-class premise law §4 derives; the mass rescale (already structurally absent in
currency); `ONE_SIDED_RNA` + `rna_one_sided` (no shipped writer since the fan-out deletion — the
NB count claim is the two-sided upgrade); the relay's Gaussian splice-in claim (superseded by the
route-sum-correct flux stream); `splice_out_rna_logvar` (already unreferenced). The two premise
estimators' recorded arguments genuinely conflict (P1d's "no common effect ⇒ unweighted mean" vs
currency's "unweighted floors at zero on sparse hops ⇒ precision-weighted") — they estimate
different targets on different populations, and §4's per-class fit is the reconciliation: weighted
moments WITHIN a class, refusal below a population minimum, no pooling across classes that split
30× (the terminus-vs-sj split P1d itself records).

## 4. The transport-dispersion decomposition — the derivation this all waits on

The one number nobody has explained: two complete flanks of one exon disagree about the same RNA
by sd ≈ 0.3–0.65 in log at depths where counting noise is negligible. Everything the owner called
overengineered is scar tissue around that unexplained number, and the premise law of §3 is
un-derivable until it is decomposed. The program — measurement first, on certified truth, no
solver anywhere:

Enumerate the candidate sources of flank↔exon (and hop-ratio) disagreement, and classify each as
GEOMETRY (calculable from the opportunity model — stops being "dispersion" at all), STRUCTURE
(annotation), or STATISTICS (residual randomness the premise must carry):

| candidate source | class if confirmed | how to measure its share |
|---|---|---|
| unannotated junctions (missing routes at one flank) | STRUCTURE | truth: per-flank route census vs oracle molecule paths |
| transcript termini inside the exon (licensing leaks) | STRUCTURE | truth: molecules entering one flank only, per exon |
| fragment-length edge effects at short exons (`P(w ≤ ell)` truncation — the SAME factor already indicted in currency's ratio, 11.6× at 98 bp) | GEOMETRY | re-derive opportunity with the exact pmf per exon length decile |
| shared fragments between flanks of short exons (correlated, violates independent-disagreement) | STATISTICS (correlation term) | count fragments crossing both junctions |
| sj strand-column mismatch (summed vs matched) | GEOMETRY | re-pool per strand, measure the residual move |
| within-exon capture profile (probe placement) | GEOMETRY under capture | OFF-vs-ON share of the residual |
| isoform mix (different routes, same exon, different true rates) | STATISTICS (the true premise) | what remains after all of the above |

The deliverable: per-class variance shares, capture-OFF and capture-ON separately. Whatever is
GEOMETRY gets computed and removed from the dispersion. Whatever is STRUCTURE becomes a licence.
What remains IS the premise, fitted per hop class by weighted moments with refusal — and at that
point the three-estimator max, both global premise scalars, and the scatter-node "law" all reduce
to one derived quantity with one estimator.

### 4a. RUNG 1 RESULTS (measured 2026-08-26; `scripts/design/transport_dispersion.py`, six ladder conditions,
two-complete-flank exon pairs, deep = flank flux ≥ 100)

⭐⭐⭐ **THE HEADLINE: the "certified scatter ~0.44" was ~95 % the TRUTH ESTIMATE'S OWN COUNTING
NOISE, not transport dispersion.** The deep pairs' median certified-mature count is only 11–19
fragments, so `ψ'(n_mrna+½) ≈ 0.19–0.23` — sd 0.44 — was riding every certified-vs-truth reading.
On well-measured exons (n_mrna ≥ 100) the common-mode residual collapses to Var ≈ 0.009–0.013 at
capture-OFF. ⛔ Every prior "intrinsic dispersion ~0.4–0.65 sd" quote is retired; instruments
comparing against certified per-object truth MUST subtract the truth count's trigamma.

**The decomposition (deep pairs, per-flank log scale):**

| term | capture-OFF | capture-ON | class / verdict |
|---|---|---|---|
| flank counting | 0.008–0.011 | 0.007–0.009 | statistics, derived (trigamma) — exact |
| truth-side counting | 0.19–0.23 (of Var(c)) | 0.18–0.22 | an INSTRUMENT artifact, not a model term |
| length-curve center (common-mode) | −0.08 (ℓ<96) … +0.10..0.15 (ℓ>96) | −0.29..−0.14 short, ~0 long | **GEOMETRY, calculable** — the `P(w ≤ ℓ)`-family truncation; both flanks move together, so the pair estimator cannot see it; the exact transcript-coordinate opportunity derivation is rung 2's first job |
| residual common-mode | **0.009–0.013** | **0.069–0.093** | OFF: small residue; ON: PER-EXON capture-enrichment scatter (sd ~0.26–0.31) — the owner's "enrichment is a spectrum", measured as a variance, shared by both flanks |
| asymmetric flank excess | **0.015–0.019** | **0.033–0.054** | OFF: the true isoform/route premise (sd ~0.12–0.14); ON: + per-flank probe geometry |
| termini inside the exon body | exactly 0 exons on this panel | 0 | STRUCTURE — absent from the panel by construction; ⚠ a panel-validity note (real data has them; the panel cannot price this licence) |
| credit-leftmost deposit rule | no L/R center split at any length | same | REFUTED as a source — the reach-capped opportunity prices the rule correctly |
| sj strand-column matching | mis-keyed experiment (naive per-sj column matching exploded Var 0.04 → 0.31 at ss.99, ≈ no-op at ss.50) | — | OPEN — the substrate's count columns are genome-strand, not transcript-strand; needs the accumulator's column semantics before re-measuring |

**What this settles for the premise law (rung 2):**

1. ⭐ **The route-pair estimator was RIGHT at capture-OFF all along**: honest asym excess
   0.015–0.019 ≈ its readings (½·MAD² ≈ 0.011–0.016). The "it under-reads 0.44" charge is
   withdrawn — 0.44 was the contaminated target, not the truth.
2. ⭐⭐ **The common-mode blindness is now a NUMBER, and it is capture-ON-specific**: the pair
   estimator misses 0.069–0.093 of per-exon enrichment-scatter variance at ON (plus the
   length curve everywhere). The ON premise needs a second, common-mode term; the OFF premise is
   essentially the pair estimator alone.
3. The single-flank population reads ~0.05 BELOW the pairs' center (measured −0.05 relative
   offset OFF) — the recorded single-flank candidate, now sized.
4. Shallow pairs (flux < 100) are 100 % counting — no premise may be fitted on them (the
   population-minimum refusals were right).

Raw rows per condition: the instrument's `--out-dir` npz files (regenerable in ~4 min).

### 4b. RUNG 2 — THE DERIVATIONS AND THE TRANSFER-VARIANCE MODEL (2026-08-26)

**Derivation A — the length curve attributed, in absolute fragments** (`rung2_geometry.py`
pattern: per-transcript placement intensity λ_t = observed mature fragments / L_eff(spliced
length, TRUE pmf) — the panel writes both — then the crossing side `y_J = log(flux/(Σλ·A_J))`
and the contained side `z_e = log(n_mrna/(Σλ·eff_r))` are absolute, no free constant):

* The CROSSING model `A_J` carries a flank-exon-length mis-shape: y-medians −0.07/−0.08 at the
  shortest src/dst exon quintile → +0.05 at the longest, both capture states. Mechanism
  candidates, both calculable: the credit-leftmost cap absent from `reach_lo` (placements
  crossing the previous junction are in A but never in flux) and the max-over-transcripts reach
  (a short-continuation isoform's placements are over-counted). Per-junction spread at depth:
  MAD 0.099 OFF / 0.134 ON.
* The CONTAINED model carries a FITTED-PMF error: z = −0.143 with the fitted pmf vs −0.070 with
  the TRUE pmf at OFF — half the contained-side deficit is the fl model's pmf, alone worth ~7 %.
* ⭐ At capture-ON the contained side is ENRICHED relative to crossing, +0.31 at the shortest
  exon quintile → +0.02 at the longest: contained placements sit fully ON-probe, crossing
  placements half-off. This is the mechanism of the ON length curve AND the ON common-mode
  scatter — a per-exon contained-vs-crossing capture ratio, i.e. exactly the LOCAL enrichment
  quantity the owner ruled must be learned locally (and the currency knob's estimand).

**Derivation B — the premise's structure** (excess by depth band × route structure, g00):

* OFF single-route pairs: 0.000–0.016 — near ZERO. A single-route exon's two junctions abut the
  same exon and share their geometry error; it cancels in D.
* OFF multi-route pairs: 0.019–0.070 — 3–10× the single-route rows. Route mixtures surface the
  per-route A-error differentials. ⭐ So the OFF "premise" is largely GEOMETRY passing through
  route structure, and its honest statistical key is the ROUTE-STRUCTURE CLASS, not depth.
* ON: both classes lift (single 0.02–0.12, multi 0.03–0.28) — the per-junction probe geometry.
* ⛔ The naive precision-weighted moments over ALL pairs under-read the deep target 2–3×
  (weighting concentrates on ultra-deep pairs, which sit in the cleanest class) — the weighting
  is right WITHIN a class and wrong across classes. This resolves the recorded P1d-vs-currency
  estimator tension: weighted moments PER CLASS, never pooled across a 3–10× structural split.

**THE MODEL.** Scoring machinery unchanged (NB `size = flux + ½` ⊗ scatter nodes). The scatter:

    V_transport(exon) = V_route(class(exon)) + V_common
    V_route  = max(0, Var_w(D) − mean_w(v_cnt))/2   weighted moments, w ∝ 1/v_cnt, fitted PER
               ROUTE-STRUCTURE CLASS (single-route pairs | multi-route pairs), refusal below a
               population minimum, floor 0
    V_common = max(0, V_tail − V_route_pooled)      the pair-blind half: V_tail = the SAME
               moments law on the NEGATIVE obs-vs-prediction residuals (gDNA-safe side), with
               obs-counting subtracted — measures asym + common; the subtraction isolates common

* ⭐ ONE estimator form everywhere — `max(0, Var_w(residual) − mean_w(v_count))`, precision-
  weighted within its population, counting-subtracted, refusing below a minimum — applied to two
  residual sources (pair D's per class; negative obs residuals). It replaces `route_pair`'s MAD,
  the left-tail MAD, the guarded fit's width role, currency's global premise scalar and relay's
  P1d in one law.
* ⭐ The shipped `max(V_pair, V_tail)` is DERIVED, not a hack: V_tail sees asym+common, V_route
  sees asym; total = V_route + (V_tail − V_route)+ ≡ max of the two, each covering the other's
  blind spot (tail starves at high gDNA → falls back to V_route; pair is common-blind → topped
  up by tail). The guarded fit keeps its CENTER role only; its sd leaves the width.
* Validation against certified truth (§4a): predicted V_route(single) ≈ 0.00–0.02,
  V_route(multi) ≈ 0.02–0.07, V_common OFF ≈ 0.01, ON ≈ 0.07–0.09 — the model's terms match the
  decomposition's certified rows term by term.
* ⛔ The current shipped left-tail reads 0.08–0.20 at OFF where the certified transport total is
  ~0.03 — it carries un-subtracted obs counting and the (calculable) length curve. The new
  V_tail's counting subtraction fixes the first; rung 3's geometry corrections remove the second.

**Recorded for rung 3 (calculable geometry, each with its measured size — src changes, priced
on the panel before landing):** ① the per-route A_J corrections (credit-leftmost cap on
`reach_lo`; the reach MIXTURE over transcripts instead of the max) — worth ±5–8 % by flank
length and most of the multi-route premise; ② the fl model's pmf error — ~7 % on every
contained opportunity; ③ the contained-vs-crossing capture ratio — the ON length curve and
common mode, learnable LOCALLY per exon (the unified knob's estimand). After ①–③ the honest
premise shrinks toward the single-route floor (~0.01) and the model's fitted share shrinks with
it — the derivation-based teardown the owner asked for.

## 5. The campaign rungs

0. ✅ **Repair the pricing instruments** (landed 2026-08-25, same day as the audit that found
   them): `ladder_arm_ab`'s msgfree/msgscale arms now mute/temper `lam_rows` at delivery (the
   shared `_install_flux_hook`; the `backbone` ≡ `msgfree_all` identity is restored),
   `backbone_relay` asserts the relay CARRIES the flux evidence, and a new pure-config
   `anchor_off` arm prices the stream on the ladder; `backbone_parity` reuses the captured
   shipped policy's flux so `no_certified_flux` is a real arm (validated: 60,403/421,056 output
   elements differ vs head) and refuses it when vacuous; `toy_trace_error` replays
   `intron_prior ⊕ lam_rows` (fidelity 0.00e+00) with separate −factory / −stream arms;
   `psi_channel_ablation` gained the `flux` ablation (identity gate bit-for-bit; −5.0 % on its
   default condition, 11.6k slots reachable); `relay_pool_ab` / `calibration_walk` carry the
   changed-meaning notes. Two new gates: the `certified_flux` switch silences the stream with
   evidence present, and currency silently drops it — both perturbation-verified. In place of a
   flux-aware `RelayPolicy.name` (golden churn), `backbone_relay`'s flux assertion and
   `backbone_parity`'s captured-flux arms close the audit's identity hole.
1. ✅ **The decomposition of §4** — landed 2026-08-26, results in §4a. The premise's honest
   targets: OFF asym 0.015–0.019 (the pair estimator confirmed), ON + common 0.069–0.093
   (the quantified blindness); the length curve goes to geometry (rung 2's derivation).
2. ✅ **Derive the premise law** — §4b (2026-08-26): the one moments law, class-keyed, with
   V_common by tail-subtraction; the max() derived; installation is rung 3's first commit.
3. **The unified hop law** (§3) as the rebuilt policy, gated stream-by-stream against silent AND
   against both existing policies per stratum, zero controls, pure-gDNA/RNA-bearing split. Acceptance:
   ≥ silent on every in-scope stratum (the bar neither policy meets today) — that is what
   "honest precision" must buy, or the layer's honest size is smaller than what ships.
4. **Converge and delete** per §3's teardown list, verdicts recorded, the losers' epitaphs kept.


## 6. THE FOUNDATION SPEC and the code interrogation (owner architecture, ratified 2026-08-26)

⭐⭐⭐ **The spec is CODE**: `src/rigel/calibration/messages/foundation.py` — read its module
docstring as the architecture document. One Message type (three unspliced lanes = AXIOM 0's
populations, two spliced lanes = certified RNA per strand; provenance is the field name, a
spliced gDNA lane is inexpressible), one propagation-time `PropagationModel` (the pass-through /
integrate rule is the base template; `attenuate` is the propagation-variance override point), one
solve-time `SolveModel` (`solve_unspliced` supplies ψ's Gaussian channels, `solve_spliced`
the row-factor likelihood — the certified-flux treatment is an implementation of THAT). The base
classes enforce what no implementation may break: sender publishes unchanged, no single-source
amplification (refused at runtime), lane isolation in transit, pass-through at empty nodes.
Gate: `tests/calibration/test_message_foundation.py`, every rule perturbation-verified.

**The interrogation — where every existing piece stands relative to the spec:**

| CONFORMS (the spec found, not invented) | | 
|---|---|
| `sweep.solve_chain` | IS the two-timepoint machine: scans = propagation time, deliver + final ψ = solve time; the backbone needs no change |
| phase-A `build_region_init` | the own belief |
| currency's 6-number hop state | the unspliced lanes' shape, exactly; its per-hop own-fuse is the reconcile template |
| the flux stream's arithmetic (`rna_anchor`) | a `solve_spliced` implementation, one hop deep |
| `inv_sj_lo/hi`, per-face sj counts | the spliced lanes' raw observations, already in the context |
| `NeighbourState` source-indexing | enforces sender-publishes-unchanged by construction |
| mismatch deflation (as a concept) | recipient-side weighing — legal, belongs in `SolveModel` |

| VIOLATES (repair = re-expression against the spec) | |
|---|---|
| relay's 10-array state, 3 precision streams, 20 operators | no committed message type; propagation and credit smeared together |
| the frame pair reading `belief_fg` at both ends | the sender reading the recipient's belief — the recorded 57–77 % debt |
| hop-time dampers (P1d, conservation, transfer var) deciding the recipient's weighing at the sender | that weighing belongs to the solve; transit costs belong to `attenuate` |
| `PsiMessage.lam_rows` as a side lane | transitional: the spliced claim outside the message; folds into the spliced lanes when the unified policy is built |

| DISSOLVES (no counterpart in the spec) | |
|---|---|
| splice-in / splice-out operators | the spliced lanes carry flux natively; population bookkeeping becomes lane arithmetic |
| the certified-RNA Gaussian channel + `ONE_SIDED_RNA` / `rna_one_sided` | subsumed by `solve_spliced` (a count likelihood, two-sided, honest) |
| mass rescale | no place in the spec (currency already proved its structural absence) |
| the three-estimator max as an unexplained vote | becomes one `solve_spliced` variance model (the §4b law) |
| mode-fusion vs measurement precision streams | one precision per claim per lane |

**The build order under the spec** (each an open problem solved separately, in the owner's
sequence): ① the foundation — DONE, gated; ② the propagation variance model = ONE `attenuate`
implementation, derived (the §4b moments law is the candidate for the unspliced lanes; hop
distance and population licences enter here); ③ the solve = ONE `SolveModel`, built only once
the architecture has settled — the owner's warning stands: the solve weighs structure, two
messages and own belief, and is the hardest part. The unified policy then implements the
backbone's existing `Policy` protocol by delegating scan steps to the `PropagationModel` and
deliver to the `SolveModel` — the shipped policies stay untouched until it beats them on the
panel, per stratum, zero controls held.


## 7. THE UNIFIED-POLICY ROADMAP (owner request, 2026-08-26) — one policy from three strengths

The goal: ONE propagation policy and ONE solve, built against the foundation spec, integrating
what each existing piece is best at. What each contributes:

| from CURRENCY | from RELAY | from the SPLICED ANCHOR |
|---|---|---|
| the message shape (already the spec's unspliced lanes) | the per-strand population licence (the three-case rule) | the spliced lanes and their count-likelihood solve |
| the knob: enrichment as ONE unknown, no switch | mismatch-yielding at the recipient (a solve-time behavior) | the §4b variance law (route classes + the tail top-up) |
| the every-hop premise, precision-weighted | the operator-pricing METHOD (arms per mechanism) | the estimator discipline: guards, refusal ⇒ silence |
| belief-free frames | the structural licences (terminus flags, admissibility) | route-sum pooling and the geometry corrections of §4b |

**Step 0 — terminology.** ✅ 2026-08-26: propagate (was reconcile), solve (was credit), attenuate
(was discount), pure-gDNA / RNA-bearing (was DELIVER/REFUTE — which also collided with the
backbone's `deliver()` method). Instrument keys renamed (`pass0_claimed_ab` score keys are
`pure_gdna` / `rna_bearing`), self-test green.

**Step 1 — THE BRIDGE (next).** One spec-conforming policy runner: a Policy implementation that
plugs a `(PropagationModel, SolveModel)` pair into the backbone's existing protocol, and builds
each node's OWN Message from the banks the context already carries (a boundary's sj flux enters
its spliced lanes — which retires the `lam_rows` side lane architecturally). Two identity anchors
falsify it before anything new is claimed: ⓐ PassThroughPropagation + a silent SolveModel must be
byte-identical to SilentPolicy on the panel; ⓑ a SolveModel whose `solve_spliced` reproduces the
shipped anchor arithmetic from the messages' spliced lanes must be byte-identical to today's
relay flux rows (the parity-gate pattern exists). ⚠ Step 1 also settles the spec's one known
gap: a claim crossing an enrichment frame needs a REFRAME in the value before the fuse — where
that lives (inside `attenuate`'s contract or as its own extension point) is decided here, with
currency's belief-free face totals as the frame source.

**Step 2 — PROPAGATE (the unspliced attenuate).** The first real propagation variance model,
derived: currency's knob for the frame decision, the per-hop premise by weighted moments PER HOP
CLASS (the §4b estimator form), relay's population licences as structural zeroes. Priced per hop
first (`hop_currency`/`backbone_parity` style), then the ladder, per stratum.

**Step 3 — SOLVE.** The hardest part, built only when the architecture has settled and designed
WITH the owner: how a node weighs structure, two travelled messages, and its own belief.
Integrates the §4b spliced law (`solve_spliced`), mismatch-yielding as recipient behavior, and
the refit-vs-message arbitration (the ROADMAP's standing open question — two imputations at one
slot with nothing arbitrating them ends HERE, because the solve is the arbiter). Design doc and
owner review BEFORE code.

**Step 4 — THE CONTEST.** Unified vs silent vs relay vs currency, all 16 conditions, per
stratum, both zero controls, the pure-gDNA/RNA-bearing split never pooled. Acceptance: ≥ silent
on every in-scope stratum (the bar nothing meets today). Converge-and-delete the losers, verdicts
recorded, per the standing rule.
