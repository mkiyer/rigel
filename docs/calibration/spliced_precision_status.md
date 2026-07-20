# The unstranded RNA-evidence problem — status, challenges, and open questions

**Status:** active working notes (2026-07-18). This document records the current state of the calibration
message-precision work: what is derived and landed, what the remaining challenge is (the spliced-precision term, the spliced RNA
evidence), and the precise questions we are prototyping to answer. It is the companion "where we are and what we
are stuck on" note to the reviewer-facing derivations
([`message_precision_derivation.md`](message_precision_derivation.md), the full architecture; and
[`spliced_message_precision_synopsis.md`](spliced_message_precision_synopsis.md), the external-review synopsis).

---

## 1. The problem in one paragraph

Rigel's calibration deconvolves each genomic node's **unspliced** fragment mass into gDNA vs RNA on the axis
`λ = logit(f_g)`, by a forward-backward **belief-propagation (BP)** pass over a region↔boundary chain. The hard
case is the **unstranded, hybrid-capture** library: the strand tilt (the only intrinsic composition signal) is
**identically zero at κ=½**, so a node has *no intrinsic evidence* about its own composition. Composition can
then only come from three places: **structure** (intergenic = pure gDNA), **cross-node messages**, and the
**Phase-2 gDNA hyperprior**. The recurring failure mode is that the message layer *manufactures* confidence out
of the shared prior and propagates it, producing confident-wrong nodes that then corrupt the hyperprior fit.

---

## 2. What is derived and landed (the τ evidence compiler + the frame-fix)

**The reference-free evidence compiler `τ`.** Message precision must not be sourced from a node's *belief*
variance (which contains the shared Beta(½,½) reference and all pooled inflow) — doing so lets a
composition-vacuous chain pool the reference into confidence (a 35× phantom variance collapse, node 2931:
`var 2.8 → 0.08`, confident-wrong at oracle `f_g=0`). Instead we maintain a **reference-free evidence
precision `τ`** per node, seeded only from genuine composition evidence, and source message precision from `τ`.
`τ = 0` (⇒ emitted precision 0) exactly when there is no evidence. Implemented behind `_TAU_PRECISION` in
`bp_solver._scan`.

**The frame-fix (the debugging breakthrough).** The τ cavity accumulation `τ_dst += Σ pr` was adding message
precision in the **log-fraction frame** to the **λ frame** τ — a coordinate mismatch that amplified `1/(1−f_g)²
≈ 4×` per hop, so a tiny seed self-bootstrapped and saturated at `1/σ²_transfer` (measured: `τ → 12`). The fix
is a **Jacobian conversion** (`τ += pr·(1−f_g)²` for gDNA factors, `·f_g²` for RNA factors) — no tuned
constant. A toy (`scratchpad/tau_toy.py`) and the real pass both confirm: with the fix, a seed *decays* along a
vacuous chain instead of exploding.

**Measured effect (τ-core, frame-fixed, the spliced-precision term still off).**
- Node 2931 phantom collapses to the honest reference (`var 0.08 → 2.36`, `f_g 0.61 → 0.49`).
- Median belief-variance sharpening collapses across the suite (gdna300: `11.6× → 1.9×`).
- The **informative** low-gDNA hyperprior fits improve: gdna1% hybrid-NPMLE L1-to-oracle `0.85 → 0.37`;
  gdna5% `0.97 → 0.59`. (The zero-gDNA L1 stays ~2.0 — a degenerate-oracle artifact, not a phantom measure.)

This is the correct architectural foundation, validated. **But it is only half the machine.**

---

## 3. The current challenge — the spliced-precision term, a **mature-RNA measurement** into exons

**What the spliced-precision term is (the owner's corrected model — supersedes the earlier "RNA-vote / nascent" framing).**
Spliced fragments are **mature** RNA — the *antithesis* of nascent. Mature RNA lives in **exons**, never in
introns. So the spliced-precision term is **not** a source of nascent RNA and has nothing to do with the intron chain. It is a
**direct measurement of mature RNA that a splice-junction boundary emits to its transcript-continuous exon**:

- **Direction:** boundary → region (the **boundary is the source**; the exon is the sink). One-directional,
  routed only to the flank carrying the junction's motif strand.
- **Measurement, not prediction.** A spliced read is *guaranteed* pure RNA (only a processed mRNA can splice),
  so the boundary does not impute — it **measures** `ρ_mature = S / E_spl` with honest **count power `S_eff`**.
  This is a legitimate exception to count-zero-information: for the spliced pool the count *is* the confidence.
  The current code throws this away — spliced fragments are lumped with the unspliced mass and the message
  degrades to a low-confidence *prediction* ("I see this much predicted RNA") instead of a high-confidence
  *measurement* ("I **measured** this much RNA"). Restoring that measurement power is the whole point.
- **Mature is a COMPONENT, not a vote on `f_g`.** The exon's unspliced pool is `N = M_body + Nn + G`
  (mature-in-body + nascent + gDNA); `f_g = G/N`. The measurement pins the **mature** part; gDNA
  (neighbours + Phase-2 hyperprior) sits **on top**. It therefore *cannot* over-crush a gDNA-rich exon — it
  only removes the measured mature from the ambiguous residual. (This dissolves the earlier "`S ⊥ f_g`
  over-push" worry, which came from wrongly modeling mature as a competing vote.)

**Why it is load-bearing.** An exon has **no direct evidence of splicing of its own** — it measures only
unspliced fragments, which are ambiguous (nascent / mature / gDNA). The adjacent junction is the **discriminator**:
a fragment that splices *into* the region is known RNA; one that crosses unspliced is very likely gDNA. Without
feeding the spliced measurement in, a transcribed exon has a gDNA channel and no RNA counterweight and is
phantom-crushed to gDNA (toy §3.1: `gdna300 normal exon` → `f_g = 0.998` with no the spliced-precision term, vs truth `0.850`).

**The strength (settled by prototype).** The strength is **`ρ_mature · E_exon`** — the mature the *recipient
exon* actually holds — **not** the raw crossing count `S`. "The boundary reports the density; the recipient
handles it with its own effective length." Prototype `scratchpad/spliced_toy3.py`:

| scenario | truth `f_g` | no the spliced-precision term | `w = S` (raw) | **`w = ρ_mature·E_exon`** |
|---|---|---|---|---|
| nascent-rich normal exon (M=900, Nn=55, G=7) | 0.007 | 0.850 | 0.014 | **0.009 ✓** |
| gdna_none normal exon (M=900, Nn=100, G=0) | 0.000 | 0.500 | 0.003 | **0.002 ✓** |
| gdna300 normal exon (M=150, G=850) | 0.850 | 0.998 | 0.893 | **0.849 ✓** |
| fast-splice gDNA-rich SHORT exon (M=10, G=900) | 0.989 | 0.998 | **0.899 ✗** | **0.987 ✓** |

The `fast-splice gDNA-rich SHORT` row is the discriminator: its mature all splices *out* (tiny `E_exon`), so
its body is genuinely gDNA. Raw `S` over-pushes it to RNA (`0.899`); `ρ_mature · E_exon = 10` correctly keeps
it gDNA (`0.987`). The recipient's own `E_exon` bounds how much mature it can hold — exactly the
"recipient handles the density" principle.

---

## 4. The questions we are trying to answer

**Q-A (the strength) — SETTLED (prototype §3).** The strength is `ρ_mature · E_exon` (measured density × the
recipient's own effective length), with confidence from the spliced count power `S_eff`. Not raw `S` (over-pushes
short exons), not a saturating constant (would be a magic number), not a mode on `f_g` (mature is a component,
not a vote). No tuned constant.

**Q-B (mode vs precision split).** Mode = `ρ_mature = S/E_spl` (the measured density); precision = `S_eff`
attenuated by `σ²_transfer(b→e)`. The recipient converts the density to a mature pseudo-count via its own
`E_exon`, but the *confidence* is capped by the measuring count `S_eff` (you cannot be more certain than the
`S` reads that measured `ρ_mature`). Prototype the exact cap: `w_eff = min(ρ_mature·E_exon, S_eff)` vs folding
a Gaussian on `log ρ_r` at mode `log ρ_mature`, precision `S_eff`.

**Q-C (the τ-arbitration interaction).** Because mature is a *component* (gDNA sits on top), a confident mature
measurement no longer locks the node against the hyperprior: in `gdna300` the exon still lands at `f_g=0.849`
(the gDNA in `g_ev` wins on top of the pinned mature). So the earlier Q-C risk (a high-τ RNA junction that
Phase 2 cannot pull back) is *resolved by the component structure* — verify this holds in the real chain, that
the spliced-precision term raises the RNA-side precision without freezing the gDNA-side `τ` the hyperprior needs.

**Q-D (the clean relay, secondary).** Through a fragment-free tiny exon, the mature message **relays onward**
undamped to the next region that has exon body (the tiny exon has no `E_exon` to hold mature, so `M_body → 0`
locally and the density passes through). Adopt the reviewer's evidence-only relay `pr_out = 1/(1/τ_cav +
σ²_transfer)`, reconciled with the frame-fix. Mechanical once Q-B is fixed.

---

## 5. How we are answering them

**Prototype-first, on toys with known truth** (the owner's method), before any wiring into `_scan`:
- `scratchpad/tau_toy.py` — the τ-bootstrap and its frame-fix (done; §2).
- `scratchpad/spliced_toy.py`, `spliced_toy2.py` — the vote-on-`f_g` forms (done; **refuted** — mature is a
  component, not a vote).
- `scratchpad/spliced_toy3.py` — the corrected mature-measurement model with strength `ρ_mature·E_exon` (done;
  §3). Passes all regimes: nascent-rich → RNA, gdna300 → gDNA (no over-crush), fast-splice-short → gDNA,
  gdna_none → RNA — no tuned constant.

**LANDED (behind `_ISPLICED_PRECISION`, default off).** The change is **precision-only** on the boundary→region
RNA message; the mode already carries the spliced density (`rho_pos += SPs/esp`). Each of the three RNA
precision terms — the two per-strand emissions and the RNA-total λ-factor — gains the spliced **direct-count**
precision `S_eff/(1 + S_eff·σ²_transfer)` (no deconvolution variance; monotone in spliced count; `S=0 ⇒ no
change`). A/B on the pure-mature expressed-exon chain (`scratchpad/ab_ispliced.py`):

| regime | exon `f_g` OFF→ON | exon confidence (1/var) |
|---|---|---|
| stranded (κ=0.95) | 0.3000 → 0.3000 | **2.11×** |
| unstranded (κ=0.5) | 0.3000 → 0.3000 | **2.34×** |
| unstranded, no spliced | 0.3000 → 0.3000 | 1.00× (identical) |

Mode unchanged everywhere; the exon's RNA confidence rises most in the unstranded regime (spliced is the only
RNA evidence there); the no-spliced control is byte-identical; introns (pure gDNA) do not regress. Calibration
suite: 227 pass / 1 pre-existing known-red, unchanged.

### 5a. The two modeling choices — SETTLED (2026-07-19)

**Precision combination = precision-add** (`pr += S_eff/(1+S_eff·σ²_transfer)`). Combining two independent
Gaussian evidence sources sums their precisions — the standard Gaussian-BP update. The density-weighted
variance alternative injects the *running* density modes (ρ_u, ρ_s) into the precision layer, a
belief→precision feedback loop — exactly the phantom-confidence class the τ compiler exists to kill. `σ²_transfer`
is belief-free (the NPMLE projection, fit once), so the spliced term has no feedback. External review concurred.

**`S_eff` = raw spliced mass `SPs` (the accumulator's Σw), not the integer flux and not Kish — for now.** The
payload carries Σw (`boundary_mass`, float) and the integer event count (`boundary_flux`), but **not** Σw² —
so exact Kish `(Σw)²/Σw²` would need a new C++ accumulator channel (+ byte-for-byte Python-reference match).
Crucially, for weights ≤ 1 the ordering is provable: **raw mass ≤ Kish ≤ flux**. On dev data spliced reads are
heavily fractional (mass/flux ≈ **0.39**), so `boundary_mass` is already the ambiguity-discounted count — the
*conservative floor*, the safest against manufacturing certainty. Kish sits **above** raw mass, so it would
*raise* precision and worsen the gdna300 over-crediting (§5b). Exact Kish is therefore a precision-*recovery*
option to revisit only if the solve looks *starved* of spliced evidence — not a safety fix. (An external review
argued raw mass "inflates"; that inverts the ordering — the naive over-count is the integer flux, which we do
not use.)

### 5b. Real benchmark A/B — the crux tension (`scratchpad/ab_ispliced_bench.py`, quick_3to1_5mb)

Calibration-pool deconvolution error vs oracle (pre-EM, un-confounded), Δ = ON − OFF (negative = better):

| condition | Δ mwae_fg | gDNA→RNA leak | verdict |
|---|---|---|---|
| **none ss0.50 capture_on** (zero gDNA, unstranded) | **−0.035** | — | **WIN** (phantom removed) |
| **none ss0.50 capture_off** | **−0.024** | — | **WIN** |
| gdna300 ss0.99 / none ss0.99 (stranded) | ≈0 | flat | neutral ✓ |
| **gdna300 ss0.50 capture_on** (fast-splice gDNA-rich) | **+0.039** | **+86k** | **REGRESSION** |

Isolation (sub-flag A/B): the **per-strand emission** drives both the win and the regression; the RNA-total
λ-factor boost is **inert** on the region metric (`chain_region_deconv` reads the per-node ψ-solve, not the λ
belief). The regression is **not** a flaw in the message: `spliced_toy3.py` shows the identical mature
measurement lands the *correct* gdna300 `f_g=0.849` **when a confident gDNA counterweight is present**. Real
pass-1 runs with only the weak pass-0 `DensityNPMLE`, so the confident mature message has no gDNA counterweight
and pulls co-located enriched gDNA to RNA. This is the Pass-1/hyperprior division: **the mature measurement pins
the mature; the gDNA *amount* on top is the gDNA-hyperprior's job.** The honest evaluation of this term is
therefore *with the gDNA hyperprior in the loop* — not pass-1 alone. The zero-gDNA win + stranded-neutral is the
term doing exactly what it exists to do; the gdna300 number is the identifiability ceiling showing through a
missing counterweight, to be re-measured once the hyperprior fit is online.

**Remaining:** (Q-D) the clean evidence-only relay; re-run this A/B with the gDNA hyperprior online (the real
test); the pass-0-vs-oracle residual-source diagnostic (which nodes dominate the error).

---

## 6. The through-line

Every past calibration failure has been the count illegitimately voting composition. The τ compiler fixed the
*belief-pooling* version; the frame-fix fixed the *coordinate* version. the spliced-precision term is a different animal: it is
the one place a count *legitimately* carries evidence, because a **spliced read is a direct measurement of
mature RNA**, not an ambiguous unspliced crossing. The discipline is: take the spliced measurement at honest
count power (`S_eff`), let the **recipient exon** scale it by its own `E_exon` (so short exons that splice
their mature *out* are not over-credited), and treat the mature as a **component** the gDNA sits on top of —
never a vote on `f_g`. Nascent stays entirely out of this: it is the separate intron-vs-intergenic density
factory. That keeps the unstranded solve afloat without re-introducing a count-votes-composition bug.
