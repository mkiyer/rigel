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
1. **The decomposition of §4** — measurement only, certified truth, no src changes.
2. **Derive the premise law** from the surviving STATISTICS terms; install behind a flag;
   backbone-parity per slot before any panel number.
3. **The unified hop law** (§3) as the rebuilt policy, gated stream-by-stream against silent AND
   against both existing policies per stratum, zero controls, DELIVER/REFUTE split. Acceptance:
   ≥ silent on every in-scope stratum (the bar neither policy meets today) — that is what
   "honest precision" must buy, or the layer's honest size is smaller than what ships.
4. **Converge and delete** per §3's teardown list, verdicts recorded, the losers' epitaphs kept.
