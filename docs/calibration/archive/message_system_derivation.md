# The Calibration Message-Passing System — Review, Derivation, Reframe & Plan

**Status (2026-07-21):** the formal, code-grounded companion to the owner's mental model in
[`intergenic_boundary_behavior.md`](intergenic_boundary_behavior.md). It (1) verifies what a node's belief stores,
(2) states precisely what the current message system does, with `bp_solver` line refs, (3) documents the verified
**τ-gag bug**, (4) derives the **per-component density message** and its equal-/unequal-gate interpretation from
first principles, (5) reconciles the gDNA prior ("Band-Aid") as the *designed* resolver, (6) gives a phased
**implementation plan** and (7) the **open issues & questions**. Sits under `CALIBRATION_ARCHITECTURE.md` (the
count-zero-info invariant, §10 here) and `CALIBRATION_MASTER.md`.

The nucleic-acid model is **gDNA vs RNA-total**, with the RNA split (+/−) a nuisance tilt — calibration models
only RNA-vs-gDNA (the per-locus EM separates nascent from mature downstream).

---

## 1. What a node's belief stores (VERIFIED)

A node owns a **fixed** total unspliced fragment mass `M` (count). Its unknown is the **composition** of that
mass — the mass fractions `(f_g, f_+, f_−)`, `f_g + f_+ + f_− = 1`. That is **2 degrees of freedom**, and the
solver already carries them as the two natural coordinates (`simplex_logodds`, `NodeDeconv.lam_*`/`theta_*`):

- **λ = logit(f_g) = log(f_g / (1−f_g))** — the **gDNA-vs-RNA-total** axis (the "gDNA level"). `f_g = σ(λ)`.
- **θ (tilt), τ = sin θ** — the **+/− strand balance**: `f_+ = (1−f_g)(1+τ)/2`, `f_− = (1−f_g)(1−τ)/2`.

**The two precisions the owner names ARE `(pr_λ, pr_θ)`** — verified: `NodeDeconv` carries `lam_var` and
`theta_var`, and `theta_var = 0` (θ locked at ±π/2) for single-strand nodes (`strand_deconv.py:60-61`). So:

| gate `(g, +, −)` | DOF | free precisions |
|---|---|---|
| intergenic `(1,0,0)` | 0 | none — `f_g=1` locked by structure |
| single-strand `(1,1,0)` / `(1,0,1)` | 1 | `pr_λ` only (θ locked to the live strand) |
| AMBIG `(1,1,1)` | 2 | `pr_λ` and `pr_θ` |

The `NodeBelief` **storage** currently keeps three component variances `(var_+, var_−, var_g)` — a *redundant
projection* of `(pr_λ, pr_θ)`. **The owner's model is correct; the fundamental object is `(λ, θ)` + two
precisions.** (Belief may be stored as fractions or counts — interconvertible via `C_c = f_c·M`; we do **not**
store density — see §2.)

**Density is derived, not stored.** Per-component density `ρ_c = C_c / E_c = f_c·M / E_c`, where `E_c` is the
component's **effective length** (`E_g` for gDNA, `E_r` for RNA — different, because gDNA and RNA have different
fragment-length distributions; `effective_length.py`). Total density
`ρ_tot = f_g·M/E_g + (1−f_g)·M/E_r` **depends on the composition** ⇒ *total density is not conserved as the
composition changes.* This is why density lives outside the belief and is recomputed on demand.

---

## 2. What the current message system actually does (code-grounded)

Per edge (source → dst), `bp_solver._scan` emits **three per-component messages** — gDNA `(amg, apg)`, RNA+
`(amp, app)`, RNA− `(amn, apn)` — each a `(mode, precision)` on a **log-fraction** axis (the recipient's ψ applies
`−½·pr·(log f_c(λ) − mode)²`, `simplex_logodds._local_loglik_logodds:297-310`). Two things determine each message:

### 2.1 The MODE — a composition/density **hybrid**, gated by edge type
`use_shift = (not exon@src) and (not exon@dst)` (`bp_solver.py:555`):

- **CLEAN edges** (`use_shift=True`; intron/intergenic ↔ boundary, no exon endpoint): the **composition
  log-odds shift** (`bp_solver.py:583-598`) — impute the dst's *fractions* from the source composition, with an
  eff-length frame correction so the capture enrichment cancels (cliff-invariant). Derivation in §5.
- **EXON-facing edges** (`use_shift=False`; *every edge into the enriched exons that under-call*): the **density
  mode** `log(ρ_c^src / (M_dst/E_c^dst))` (`bp_solver.py:605`) — the source's per-component density over the
  **recipient's observed total density**.

So the "composition vs density" question the owner raised is currently answered *both ways, frozen per edge* — the
unresolved "back-and-forth." **On the edges that matter (into exons), we send a density-ratio anchored to the
recipient's observed total**, which is exactly what mis-signs the TSS/TES seam (§3, and
`enrichment_sensitivity_worklog.md` §8c).

### 2.2 The PRECISION — honest, but behind an EMISSION GATE
gDNA message precision `pr = 1/(Var(log f_g) + 1/M_src + σ²_transfer)` (`bp_solver.py:606`): composition variance
+ count sampling + the enrichment-crossing transfer variance. **But emission is GATED**
(`bp_solver.py:546-548`): `emit_g = (M_src>0) and lam_ev`, where `lam_ev = struct_lock[src] or (τ[src] > 0)`, and
the reference-free evidence `τ = I_strand + I_struct + relayed`. The spliced-fragment precision credit
(`pr += S_eff/(1+S_eff·σ²_transfer)`, `bp_solver.py:630`) sits **inside** the `emit_p` block — *after* the gate.

### 2.3 The cliff (non-uniformity) handling
Two mechanisms stacked: **σ²_transfer** damps any message across an enrichment gap (penalty = mode-gap²,
`transfer_variance_formal_derivation.md`), and the **mode choice** (the composition shift *cancels* enrichment
identically; the density mode *anchors* to the recipient's observed total instead of crossing it).

---

## 3. The τ-gag BUG (VERIFIED — high priority)

**On unstranded data, `I_strand = 0` (the strand Fisher info vanishes at κ=½), and the τ-seed does NOT include
spliced-fragment evidence.** So a splice junction whose only composition evidence is its motif-stranded spliced
reads has `τ = 0 ⇒ lam_ev = False ⇒ emit_g = emit_p = False` — it emits **nothing**, and the spliced-precision
credit (which lives inside the gated block) never runs. *The spliced measurement cannot open the gate it would
then credit.*

**Measured** (`scratchpad/msg_gate_verify.py`, `gdna300_ss0.50_capON_nrna`):
- The specific edge: junction **1054** (96 spliced fragments, 120 unspliced crossing mass, real gDNA) → enriched
  exon 1055: `τ_src=0 ⇒ EMITTED pr_gDNA=0.0000, pr_RNA=0.0000`. Silent.
- Aggregate: of **764 spliced-carrying boundaries** (1528 emitting edges), **52% emit zero gDNA precision, 66%
  emit zero RNA precision, 52% emit nothing at all.** Half of our splice-junction messages — our single best
  RNA signal on unstranded data — are silenced.

The τ-gate's original purpose was real (kill the phantom cascade: using the *pooled-reference belief variance*
`σ²_λ` manufactured confidence on vacuous unstranded chains — `phantom_and_tau_precision_verified`). But it was
implemented as a **binary emission gate** that discards genuine evidence, which violates the principle in §7. The
fix is to make τ a **precision magnitude** (never a gate) and to **credit spliced evidence** in it (§7).

---

## 4. The reframe — a message is per-component `(density, precision, gate)`

Adopting the owner's format. A message carries, for each component `c ∈ {gDNA, +RNA, −RNA}`:

```
  msg[c] = (ρ_c,  pr_c,  gate_c)      gate_gDNA ≡ True (structural constant, need not be on the wire)
```

- **Content = DENSITY** `ρ_c` (absolute level = the common currency across nodes; a count is eff-length-dependent,
  a density is not). The source sends what *it* measured.
- **Per-component PRECISION** `pr_c` — the source may be confident about gDNA and know nothing about RNA (or vice
  versa). A single composition fraction *cannot* carry this; per-component precision degrades gracefully to
  partial information.
- **GATE** `gate_c` = the source's **signature** (what it can measure): intergenic → `(1,0,0)`; single-strand +
  → `(1,1,0)`; AMBIG → `(1,1,1)`. Read from a constant vector, not stored per message.

**The source does NOT compute composition, and does NOT gate emission** (§7). It sends densities + honest
precisions. **The RECIPIENT interprets** in its own frame (its own `M_dst`, `E_g^dst`, `E_r^dst`, and its own
gate), by one of two rules (§5, §6) chosen by comparing gates. This is the clean separation the current
source-computes-the-mode design lacks.

---

## 5. Derivation A — EQUAL gates: the composition shift (recipient-side)

**When `gate_src == gate_dst`** (both endpoints admit the same components), the **density composition** is assumed
conserved across the cliff — only the overall scale (the capture enrichment `k`) changes (`intergenic_boundary
_behavior.md` "principle of composition"; this is the only assumption available and it is nucleic-acid-agnostic —
capture pulls gDNA and RNA proportionally):

```
   ρ_c^dst = k · ρ_c^src     for every active component c,   common unknown k > 0.
```

The recipient knows its **fixed total mass** `M_dst` and its per-component eff-lengths. The scale `k` is fixed by
the recipient's own mass budget:

```
   M_dst = Σ_c C_c^dst = Σ_c ρ_c^dst · E_c^dst = k · ( ρ_g^src E_g^dst + (ρ_+^src + ρ_−^src) E_r^dst )
   ⇒  k = M_dst / D,     D ≡ ρ_g^src E_g^dst + (ρ_+^src + ρ_−^src) E_r^dst
```

giving the recipient's mass fractions (k **cancels** ⇒ cliff-invariant):

```
   f_g^dst = ρ_g^src E_g^dst / D
   f_+^dst = ρ_+^src E_r^dst / D
   f_−^dst = ρ_−^src E_r^dst / D
```

In log-odds this is the cliff shift `λ_dst = log(ρ_g^src/ρ_r^src) + log(E_g^dst/E_r^dst) = λ_src +
log(E_g^dst/E_g^src) − log(E_r^dst/E_r^src)` (`cliff_message_derivation.md`, MC-validated). **This matches the
current `use_shift` code** (`Mg=ρ_g·E_g^dst`, `Mp=ρ_+·E_r^dst`, `mode=log(Mg/ΣM)`, `bp_solver.py:594-603`) — the
arithmetic is correct; what is new is *doing it recipient-side from the message densities* and *only when gates
are equal*.

> **Corrects the owner's sketch.** The sketch wrote `dst_ρ_g = src_f_g · dst_count_g / dst_eff_len`, but
> `dst_count_g` (the dst's per-component gDNA count) is the unknown. The correct recipient-side form uses the dst's
> **total** mass `M_dst` and the dst's per-component **eff-lengths** — the fractions above. The composition (a
> ratio) transfers; the recipient re-anchors the scale to its own observed total.

---

## 6. Derivation B — UNEQUAL gates: the gDNA-only lower bound (intergenic↔exon)

**When `gate_src ≠ gate_dst`**, composition does **not** transfer. The canonical case is the **TSS/TES seam**:
source intergenic `(1,0,0)` → dst exon `(1,1,0)` (say + strand). The source can speak only to gDNA.

**Why the composition shift is invalid here.** Plugging `ρ_+^src = ρ_−^src = 0` into §5 gives `f_g^dst = 1` —
asserting the exon is *pure gDNA*. That is an over-claim: the source's RNA gate is **off** ("I cannot measure RNA
here"), which is **not** the same as "there is no RNA at the destination." The exon has RNA; the source simply
knows nothing about it.

**What the message legitimately says.** It conveys a gDNA **density** `ρ_g^src`. Because density is not conserved
(hybrid capture enriches the exon's *contained* gDNA far above the *crossing* background — the crossing is not
efficiently captured), `ρ_g^src` is a **weak lower bound** on the exon's gDNA density, never a point estimate:

```
   ρ_g^dst  ≥  ρ_g^src         (capture only ADDS gDNA to the exon; the un-captured crossing is the floor)
   ⇒  f_g^dst  ≥  ρ_g^src E_g^dst / M_dst      (a weak LOWER bound on f_g — usually small)
```

The RNA components are **unconstrained** by this message. Therefore the recipient must apply it as a **one-sided
(lower) anchor on the gDNA level λ, never a two-sided point estimate, and never a downward pull** (a gDNA-only
message cannot claim "little gDNA" — it can only establish gDNA *presence*).

**This is exactly where the current code fails**: the density mode uses `ρ_g^src` as a *point estimate*
(`log(ρ_g^src / (M_dst/E_g^dst))`), which for a depleted crossing (mass ~5) into an enriched exon (mass ~50k)
evaluates to `log(tiny) ≈ −7` ⇒ crushes `f_g → 0` (a pure-gDNA seam telling the exon "no gDNA"). The fix is the
lower-bound reading.

**The enrichment is genuinely unidentifiable from this message.** Exon total density 200 with a crossing floor of
2: the exon could be `gDNA=2, RNA=198` or `gDNA=200, RNA=0` — the 100× is either RNA or captured gDNA, and the
seam cannot tell (its crossing is un-captured background). **Resolution requires the gDNA prior** (§9). Until then
the honest belief is *"gDNA present at ≥ background, amount unknown"* — a wide λ posterior.

> **This is the source-physics ORIGIN of the "asymmetric-upward" projection** (`gdna_projection.py`, the
> enrichment-sensitivity "Band-Aid"). That projection empirically re-discovered that gDNA-only boundary evidence
> is a **lower bracket** (look UP, never down). Fixing the message at the source makes the asymmetry *emerge*
> instead of being imposed downstream.

**Precision here is very weak** and must be: `pr ≈ 1/(1/M_crossing + σ²_transfer + …)` with `M_crossing` tiny
(capture-depleted) and `σ²_transfer` huge (multi-decade cliff) ⇒ `pr → 0`. Per §11's summary, assigning
*unwarranted* precision to an unidentifiable node is dangerous — it cascades and blocks the prior from moving the
node later.

---

## 6A. The two modes are ONE — the shift + honest precision (retire the density mode)

Today `_scan` carries **two** message modes (§2.1): the **composition SHIFT** on clean edges and the **density
mode** on exon edges. They differ in exactly one thing — the **normalizer**. Both impute the dst's per-component
masses `M_c = ρ_c^src · E_c^dst`; then:
- **SHIFT** normalizes by the *imputed* total `ΣM_c` ⇒ `f_c^dst = M_c/ΣM_c`. Trusts the source composition; the
  enrichment scale `k` cancels ⇒ **cliff-invariant** (`λ_dst = λ_src + log(E_g^dst/E_g^src) − log(E_r^dst/E_r^src)`).
- **DENSITY** normalizes by the dst's *observed* total `md` ⇒ `f_c^dst = ρ_c^src·E_c^dst / md`. Trusts the dst's
  own total; the enrichment is baked in ⇒ **fails across a cliff** (source depleted ÷ dst enriched → `f_g→0`, the
  TSS/TES crush, `enrichment_sensitivity_worklog.md` §8c).

**Why the split exists (and why it is ad-hoc).** The shift is cliff-invariant but *trusts the source's
composition*; on an **unstranded exon** source that composition is degenerate (~0.5, unresolved), so the shift
propagated garbage confidently — the regression `cliff_message_derivation.md` §9 records. The density mode was
kept on exon edges as a **robustness crutch** ("distrust the exon source; anchor to observed data"), gated by the
`use_shift = not exon@either-end` proxy — a *reliability* heuristic, not a principled rule.

**The reconciliation — that crutch is now the PRECISION's job.** After the τ-fix the message precision is
`1/(Var(log f_c^src) + 1/M_src + σ²_transfer)`, where `Var(log f_c^src)` is sourced from the source's own
reference-free evidence `τ`. An unstranded exon source has `τ→0 ⇒ Var→∞ ⇒ precision→0` — its shift message is
**already ignored**. So we no longer need to *switch modes* to distrust an unreliable source; we use the
**composition shift everywhere and let honest precision decide trust.** The density mode is **retired** — it was
the shift's crutch from before the precision was honest. The only branch that remains is **structural**, not
ad-hoc (§5 vs §6):

- **gates match** ⇒ normalize over the transferred components = the shift (cliff-invariant).
- **gates differ** ⇒ the un-transferable components are unconstrained; the transferable gDNA imputes a one-sided
  **lower anchor** on its density (§6), never a two-sided point estimate.

So the *single* mode function is: **impute per-component masses via the eff-length frame; branch only on
gate-equality; trust via honest precision.** ⚠ **BEHAVIORAL, needs an A/B:** the memory's shift-on-exon regression
was measured at *full* precision — the honest-precision version must be benchmarked (`pass0_bench.py` +
`gdna_none` guard) before landing. Unifying the mode is therefore a *behavioral* phase, NOT part of the
byte-identical Phase-A refactor (which only *extracts* the two modes as-is).

---

## 6B. Solvability — "can I solve?" (the degrees-of-freedom criterion)

**The gate sets the DOF.** The composition `(f_g, f_+, f_−)` lives on the simplex, but the structural gate fixes
some components. With `|A| = 1 + [gate_+] + [gate_−]` active components, the composition has **`|A|−1` free DOF**,
in the two natural coordinates `λ` (the gDNA-vs-RNA level) and `θ` (the +/− tilt):

| gate `(g,+,−)` | class | DOF | free axes |
|---|---|---|---|
| (1,0,0) | intergenic / seam | **0** | none — `f_g=1` forced by structure |
| (1,1,0) / (1,0,1) | single-strand | **1** | `λ` (θ locked to the live strand) |
| (1,1,1) | AMBIG | **2** | `λ` and `θ` |

**A node is solvable ⟺ the Fisher information over its free axes is full-rank** — equivalently, **every free axis
has ≥1 nonzero-precision source.** The three information sources (the count-zero-info sources) each carry rank on
specific axes:

| source | axis it informs | rank |
|---|---|---|
| **strand tilt** (the Beta-Binomial) | the plus-fraction `p = ½f_g + κf_+ + (1−κ)f_−` | **1 if κ≠½, 0 if unstranded** |
| **a gDNA-bearing message** | `λ` | 1 |
| **a per-strand RNA message** (spliced / nascent) | `θ` | 1 |
| **the global gDNA prior** (Phase-2) | `λ` | 1 |

**Every case falls out** (deriving the owner's premise — "nonzero precision for each active component"):
- **Intergenic (0-DOF):** no free axis ⇒ self-solved with NO information (`f_g=1` by structure).
- **Single-strand STRANDED (1-DOF):** strand alone pins `λ` ⇒ solvable. (*"tilt solves 1-DOF."*)
- **Single-strand UNSTRANDED (strand rank 0):** `λ` has no pass-0 source (a TSS/TES neighbour gives only a weak
  gDNA *lower bound*, §6) ⇒ **unsolvable in pass-0 ⇒ skip; defer to the Phase-2 prior.** (the single-exon case, §8)
- **AMBIG (2-DOF):** strand gives *one* constraint (`p`) — and for a **balanced** node `p=½` regardless of `f_g`,
  so `λ` is unconstrained ⇒ **needs a second source** (a gDNA message or the prior). (*"strand needs additional
  info for 2-DOF"* — precisely the AMBIG two-root ambiguity, `CALIBRATION_MASTER.md` §4.)

**⟹ THE SKIP RULE (Phase C) — VALID for regions; boundary rule + metric pending. See [`solve_gate_design.md`](solve_gate_design.md).**
A node skips its pass-0 solve iff a free axis (`λ`; and `θ` for AMBIG) has zero total precision from {strand,
messages, prior}. **Node-type correlation test (2026-07-21) confirms the criterion for REGIONS** (77 % of mass):
the withheld single-strand + AMBIG regions are genuine *coin-flips* — their forced solve is uncorrelated with the
oracle (corr −0.04 / −0.07) — while the solved ones correlate 0.63 / 0.69. So withholding them is right; their
arbitrary default `f_g` must simply not be counted as error (an earlier mass-weighted-error "refutation" was a
metric artifact — `solve_gate_design.md` §2). **Open before it ships:** the criterion is *inverted for boundaries*
(a real bug — re-derive), and the DOF gate must count the fitted prior as a source so the Phase-2 refit *resolves*
the deferred nodes rather than re-skipping them.

---

## 7. Emission — do NOT gate; use honest precision (the τ-gag fix)

**Principle (owner): the source always emits; it never gates.** A message that carries no information carries
**zero precision**, and a zero-precision message is *ignored by the recipient* — so "not emitting" and "emitting
`pr=0`" are bit-identical. Emission may be skipped **only** as a performance optimization when the result is
provably bit-identical (all components `pr≈0`).

The τ-gate violates this: it is a **binary emission gate keyed on `τ>0`** that discards real spliced evidence
(§3). The fix keeps τ's *purpose* (a **reference-free precision** that is naturally weak on vacuous unstranded
chains — this is what prevents the phantom cascade) but changes its *use*:

1. **τ is a precision magnitude, never a gate.** A vacuous source gets `τ→0 ⇒ pr→0` (harmless: ignored by the
   recipient), so removing the *gate* does not re-open the phantom — the phantom came from the *pooled-reference
   belief variance*, not from emission. A real measurement gets `τ>0 ⇒ pr>0`.
2. **Credit spliced-fragment evidence in τ.** The motif-stranded spliced count is genuine, reference-free RNA
   evidence: it must contribute to `τ` (opening the door its own precision credit already funds). Then a spliced
   junction on unstranded data emits its DNA+RNA message — the powerful message §3 is silencing.

**This respects count-zero-info (§10):** the count enters only as *statistical power* (the spliced sampling
precision), never as a composition vote.

**Must re-verify:** the `gdna_none` false-positive guard (the phantom cascade, `τ→12` historically). Hypothesis:
it does not return, because vacuous sources still carry `pr≈0`. This is the key implementation risk (§11–12).

---

## 8. Single-exon / TSS / TES — the correct pass-0 default

On unstranded data a **single-exon transcript** has no strand and no internal junction ⇒ **no RNA measurement at
all** ⇒ `f_g` is *fundamentally unidentifiable*. Its only inputs are the two flanking TSS/TES gDNA-only messages
(§6). The correct pass-0 behavior is therefore: **a weak gDNA lower anchor, a wide λ posterior, low precision** —
*not* a confident value. Such nodes must be **withheld from gDNA-prior training** (they carry no honest
composition to teach the landscape). This generalizes to the **TSS/TES flanks of multi-exon genes** — the *outer*
boundary of a first/last exon is intergenic↔exon (gDNA-only); only the *inner* junction (exon↔intron / exon↔exon
with spliced) carries the RNA evidence that makes the gene identifiable (once §3 is fixed).

---

## 9. The gDNA prior as the resolver (pass-1) — reconciling the "Band-Aid"

The unidentifiable nodes of §6/§8 are resolved by the **global gDNA hyperprior** — `CALIBRATION_ARCHITECTURE.md`
§1.3's third information source, not a Band-Aid. The owner's worked example: an exon at total density 30, prior
with a depleted mode (2) and an enriched mode (20). The messages establish "gDNA present, weak precision"; the
prior, projected at density 30, pulls to its enriched mode 20 — *"you see 30, but gDNA can only account for ~20"*
⇒ `gDNA≈20, RNA≈10`. Correct.

**The reconciliation:** the enrichment-sensitivity work (`enrichment_sensitivity_worklog.md` §5 — the unified
landscape + the **asymmetric-upward projection**) is precisely this resolver, and its asymmetry is the §6
lower-bound physics. **Its prerequisite is that pass-0 leaves these nodes at honest *weak* precision** — not
crushed toward 0 (today's density-mode bug) and not over-committed. Fix the messages (§6–§7) so pass-0 hands the
prior an honest weak belief, and the prior does the identifiable-from-population lifting. The two halves are one
design.

---

## 10. Count-zero-information compliance (`CALIBRATION_ARCHITECTURE.md` §0)

Every element above keeps the count as **statistical power only**, never a composition vote:
- **Density message content**: the composition transferred (§5) is the *source's* composition, itself set by the
  source's strand/spliced *evidence*, not by any count. The count enters `pr_c` as sampling precision only.
- **gDNA-only lower bound** (§6): a **structural** (signature) fact — intergenic is gDNA by signature — plus a
  density floor; no count votes composition.
- **De-gating + spliced credit** (§7): the spliced *count* enters as its sampling precision (a §0 channel), never
  as "this is RNA."

---

## 11. Implementation plan (phased; each gated on the metrics below)

**Gates for every phase:** (a) the single-strand enriched under-call (`enrichment_sensitivity_worklog.md` §1,
`scratchpad/pass0_ss_dissect.py`), (b) the **`gdna_none` false-positive guard** (must not regress — the phantom),
(c) net fragment-flow (`oracle_and_benchmarking.md`), on the toy suites. Develop on toys (owner directive).

- **Phase 0 — the τ-gag fix (high priority, smallest, highest-value).** Credit spliced-fragment evidence in the
  `τ` seed and make τ a precision, not a binary emission gate (`bp_solver._scan`): a spliced junction on
  unstranded data emits its DNA+RNA message; a vacuous source still emits `pr≈0`. *Measure the 52%→ and whether
  `gdna_none` holds.* This alone should recover a chunk of the multi-exon bulk under-call.
- **Phase 1 — unequal-gate gDNA-only messages as one-sided lower anchors** (§6). Replace the exon-facing density
  *point* mode with the lower-bound reading for gDNA-only sources (never a downward pull; weak precision). Fixes
  the TSS/TES sign (the 38% wrong-signed seam messages, §3).
- **Phase 2 — the per-component density message + recipient-side interpretation** (§4–§6). Refactor the message
  to `(ρ_c, pr_c, gate_c)`; move the equal-gate composition shift and the unequal-gate lower bound to the
  recipient; select by comparing gates. Clean architecture; removes the `use_shift` edge-type proxy.
- **Phase 3 — the gDNA prior resolver + single-exon holdout** (§8–§9). Wire the unified landscape +
  asymmetric-upward projection as the pass-1 resolver on the now-honestly-weak unidentifiable nodes; hold single-
  exon transcripts out of training.

---

## 12. Open issues & questions

1. **Does de-gating re-open the phantom cascade in `gdna_none`?** (§7) The central risk. Hypothesis: no, because
   vacuous sources keep `pr≈0`. **Must measure** before committing Phase 0.
2. **Exact ψ form of the one-sided lower anchor** (§6). Is it a half-quadratic (`−½·pr·max(0, mode−λ)²`), a soft
   hinge, or a one-sided Gaussian? How does it compose with the 1-DOF λ solve and the reference arm?
3. **Where does a gDNA-only message land f_g in pass-0 — near the reference, or leaning gDNA?** The owner's note
   says "100% gDNA, weak precision"; the §6 derivation says "lower-bound + wide posterior ≈ near reference." These
   differ in the pass-0 point estimate (the *precision* agreement is what matters most). **Decision needed** — and
   it interacts with what the prior then sees.
4. **The spliced/mature vs nascent subtlety** (§3's powerful message). Spliced fragments measure **mature** RNA at
   the junction; the exon's *unspliced* mass competes gDNA vs **nascent**. How does a spliced-anchored RNA density
   transfer onto the exon's unspliced composition? (Mature presence bounds RNA *activity* and strand, but is a
   different pool than nascent.) **Needs its own derivation** — this is the crux of how much the τ-gag fix
   actually recovers.
5. **σ²_transfer at unequal gates** (§6). Is the enrichment-NPMLE projection variance still the right cliff
   damping when the source is gDNA-only and the dst is RNA-on?
6. **Belief storage: migrate to `(λ, θ)` + two precisions, or keep the 3-variance projection?** (§1) The relay
   already operates in `(λ, θ)`; the 3-variance `NodeBelief` is redundant. Refactor scope vs. churn.
7. **Recipient-side vs source-side computation** (§4). Moving the shift to the recipient is cleaner but touches
   the `_fold_lambda` EP moment-match and the message plumbing. Scope it.
8. **"Wait until the prior is fit before solving" vs. weak-precision pass-0.** The owner's note says defer solving
   unidentifiable nodes; the §9 reconciliation argues weak-precision pass-0 + prior pass-1 achieves the same with
   no separate deferral mechanism. **Confirm** they are equivalent (they should be, if precision is honest).
9. **Does Phase 0 alone move the benchmark**, or is it inert without Phase 1 (because the wrong-signed TSS/TES
   messages still dominate)? Sequencing/measurement question.
