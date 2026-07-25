# Reconciling the ABSOLUTE-DENSITY and COMPOSITION regimes — the design brief

**Status: OPEN, and it is the next substantive task.** Owner-directed, 2026-07-25. Read
`ROADMAP.md` first, then `SESSION_2026_07_25_HANDOFF_6.md`, then this. Do NOT read `archive/`.

This document exists so the work below is not lost between sessions — it was surfaced while dissecting other
bugs, and every digression since has confirmed it is the real remaining structural problem in pass-0. It
collects (§1) the problem, (§2) the measured evidence, (§3) the derivation to do, (§4) two adjacent modelling
gaps the owner identified that bound how far §3 can go, and (§5) the record of what was tried and rejected so
nobody re-runs it.

---

## 1. The problem, in one page

Pass-0 has moved from comparing **absolute densities** between nodes to comparing **compositions**, normalising
between nodes by an enrichment ratio `r = ρ_tot(dst)/ρ_tot(src)` that cancels hybrid capture. That move is
correct and is not in question — it is what made the capture arm solvable at all. But it discarded something
the density regime had for free.

**The composition regime is scale-blind.** It reasons entirely in shares: "my neighbour is 95 % RNA, and after
scaling by the enrichment ratio between us, so am I." Nothing in that chain of reasoning is aware of how much
nucleic acid was actually *sequenced* at the destination. So a claim can be internally consistent, propagate
through many hops, and still be **physically impossible** — asserting more reads of one component than the node
observed in total.

**The density regime has the missing constraint but cannot handle capture.** In absolute terms a node that
sequenced `M` fragments cannot contain more than `M` fragments of anything, and a claim that crosses a 100×
depletion cliff has essentially no reads left to support it. That is exactly the information the composition
regime throws away.

**The synthesis we want.** The composition model sets the **shares**; the density model bounds the **scale**
and, crucially, tells us **how much confidence a claim has actually earned**. The two are not competitors —
the density regime is where the composition regime's premises get audited against data.

> **The owner's framing (keep this):** *"If we go off a giant enrichment cliff and lose a hundred x of our
> signal… there's no more reads left there. You can retain this fraction that theoretically there's still RNA
> signal, but you're dividing one or two reads between DNA and RNA at that point, and there really is no count
> left, and there's no precision left for you to maintain that signal. If it does go down to very low levels,
> the precision has to go down to extremely low levels because there are no more counts left to hold onto that
> belief."*

---

## 2. The measured evidence

### 2.1 The anchor: a phantom RNA claim that percolates four hops

Condition `gdna_gdna300_ss_0.99_nrna_none_capture_on` (a library with **no nascent RNA at all**, so every
intron's true RNA is exactly 0). Chain segment, exon 2201 → exon 2197:

| node | class | observed mass `M` | **claim (reads)** | **% of M** | msg precision (= eff. reads) | ORACLE RNA |
|---|---|---|---|---|---|---|
| exon 2201 | exon | 31,720 | 27,432 | 87 % | 7,350 | 19,709 ✔ |
| bnd 2200 | boundary | **120** | **222** | **185 %** | 151 | 0 |
| intron 2199 | intron | 466 | 434 | 93 % | 66 | 0 |
| bnd 2198 | boundary | **558** | **1,027** | **184 %** | 54 | 0 |
| exon 2197 | exon | 37,761 | 32,350 | 86 % | 64 | **91** |

Exon 2197's true RNA content is **0.24 %**; it solves at `f_g = 0.822` (self-solve 0.993, oracle 0.998) because
it receives an RNA− message at `f_neg = 0.479` with precision 63.6 — a **199× over-claim**, which overrides a
*correct* gDNA message (mode 0.991 at precision 321).

Hop-by-hop, the RNA density is attenuated correctly on the way in and manufactured back on the way out:

| hop | `r` | effect |
|---|---|---|
| exon 2201 → bnd 2200 | 0.01 | ÷100 |
| bnd 2200 → intron 2199 | 0.05 | ÷20 ✔ the intron correctly kills it (85.9 → 1.11 → 0.050) |
| intron 2199 → bnd 2198 | **103** | ×103 |
| bnd 2198 → exon 2197 | 5.77 | ×5.8 |

Net ×0.30. Note what **is** and **is not** wrong here: the reframe telescopes correctly for a pure
pass-through, and `_pin_v` makes `k ∝ 1/r` so `r` cancels from the *delivered* composition. The defect is that
the relay **fuses** each node's own belief in and then re-scales the fused value out of that node's frame — so
the intron's correct, confident *"there is no RNA here"* is expressed as a density, and re-multiplying by the
outgoing `r` undoes it. At intron 2199 the honest own belief `ρ_R = 0.00091` at precision **0.242** is
overwritten by an incoming `ρ_R = 0.0502` at precision **66.06** — a **273× weight ratio**.

### 2.2 The mass bound is violated on a majority of nodes

`Σ_c ρ_c·E_c / M` for the relay's running belief (it must be ≤ 1):

| condition | median | p90 | p99 | max | % over 1 | % over 2 |
|---|---|---|---|---|---|---|
| gdna300 ss0.99 capON | 1.07 | 1.71 | 57.9 | 343 | **71 %** | 7.4 % |
| gdna300 ss0.99 present capON | 1.06 | 1.59 | 43.2 | 343 | 69 % | 5.8 % |
| gdna100 ss0.50 present capON | 1.00 | 1.60 | **287.6** | 519 | 52 % | 6.5 % |
| gdna300 ss0.99 present capOFF | 1.01 | 1.84 | 31.0 | 287 | 57 % | 8.7 % |

**`_pin_v` — the operator that enforces this bound — exists, and its own docstring says the fix is "the
design's own ÷`M_dst`, applied at EVERY node rather than only at the final combine". It is applied only in
`_transport` (the combine). The relay is unpinned.** That is the single most concrete lead in this document.

### 2.3 Where the error actually lives (suite-wide, 74,494 nodes, refit = 0)

| stratum | nodes | err-mass | share | self-solve | solved |
|---|---|---|---|---|---|
| FULL RANK (`τ_own > 0`) | 26,813 | 1.05 M | **7.8 %** | 0.011 | 0.025 |
| NOT full rank (`τ_own = 0`) | 35,877 | 12.40 M | **92.2 %** | 0.333 | 0.164 |

So messages *halve* the error on the population that cannot self-solve — they are doing their job. The
tractable bug set is **7,462 full-rank nodes the messages made worse**: 707 k err-mass, self 0.021 → solved
0.090, 81.7 % exons, 92.4 % stranded × capture-ON.

But the same **mode** defect hits the `τ_own = 0` population harder, which is why this work is worth more than
the 5.3 % bug set suggests. RNA claim delivered vs oracle RNA, on RNA-**poor** exons (oracle RNA < 10 %):

| set | RNA claim | oracle RNA | over-claim | err-mass |
|---|---|---|---|---|
| full-rank RNA-poor exons | 0.285 | 0.026 | **11×** | 357 k |
| `τ_own = 0` RNA-poor exons | 0.276 | 0.027 | **10.4×** | **1,867 k (13.9 % of suite)** |
| RNA-**rich** exons (control) | 0.632 | 0.455 | 1.4× | — |

≈ **2.2 M error mass, 16.5 % of the suite**, and it is a MODE defect — correctable in pass-0, it does not wait
for the hyperprior.

---

## 2b. DERIVATION RESULTS (2026-07-25) — steps 1 and 2 are done

### 2b.1 The mass bound is an IDENTITY under the imputation premise
Under the premise (source and destination share composition), reframing gives
`ρ_c^msg = a_c·ρ_tot(dst) = ρ_c^dst,true`, hence `Σ_c ρ_c E_c = M` **exactly**. So its violation *is* the
premise error, measured against a hard observable, with no prior needed. But the sensitivity is limited:

```
    k = M / Σ_c ρ_c^msg E_c = [ Σ_c a_c^dst E_c ] / [ Σ_c a_c^src E_c ]
```

`k` sees a share mismatch **only through the eff-length-weighted average**, so when `E_g ≈ E_r` it is blind to
composition error entirely. Max `|log k|` from composition alone: **×1.04** on a contained region (E=2701/2801),
**×1.50** at a boundary crossing (300/200), ×3.3 on a pathological short region. **Therefore the measured p99 of
31–288× and max 519× cannot be single-hop composition error — they are ACCUMULATED drift** (the multiplicative
random walk `_relay`'s own comment describes). The mass bound's value is as a per-hop **anchor that stops the
drift compounding**, NOT as a per-hop composition detector. This materially changes §3.2 — see §3.2 below.

### 2b.2 Is the single `r` mis-specified, or is the premise wrong? — THE PREMISE (settled)
The orchestrator flagged that §4.2 might subsume §3. It does not. gDNA is uniformly present along the genome,
so the **true** capture step between two nodes is `[G(dst)/E_g(dst)] / [G(src)/E_g(src)]` from the oracle
(`scratchpad/derive_1_ratio_check.py`):

| condition | model `r` vs true capture step | RNA vs the capture premise |
|---|---|---|
| capture **off** | ×1.03 (corr 0.12 — no variation to track) | **×3.22** |
| capture on | ×1.07 (corr 0.923, slope 0.979) | ×2.67 |
| capture on + nRNA | ×1.13 (corr 0.917, slope 0.991) | ×2.23 |
| verystrong | ×1.43 (corr 0.955, slope 1.026) | ×1.97 |

The reframe estimates the true capture step **well** (slope ≈ 1, corr 0.92–0.96). The channel mixing of §4.2 is
real and sits exactly where predicted — `exon↔boundary` edges ×1.3 (×2.4 at verystrong) vs ×1.0 on
intron/intergenic edges — but it is the SMALLER effect. **The clean discriminator is capture-OFF: there is no
enrichment to mis-specify (reframe error ×1.03) and RNA still deviates ×3.22.** That is pure expression
difference, i.e. the imputation premise itself. So: fix the scale anchor (§3.1) first; per-channel `r` (§4.2) is
a separate, smaller, additive improvement.

### 2b.3 WHICH HOPS break the identity — graft and peel, not plain reframes
Bit-exact offline replay of `_relay` (`scratchpad/derive_2_relay_pin.py`, validated `max|Δρ| = 0.000e+00`
against the shipped arrays), per-hop `|log k|`:

| edge kind | median | p90 | fraction > 1.5× |
|---|---|---|---|
| **plain** reframe | ×1.05–1.11 | ×1.67–2.50 | 15–18 % |
| **graft** | ×1.08–1.12 | ×11.6–**84.4** | 22–29 % |
| **peel** | ×**1.31–1.58** | ×4.9–10.8 | **42–53 %** |

Plain reframes preserve the identity to within the §2b.1 analytic bound — **the theory is confirmed exactly
where it should hold**. The violations are concentrated at the ROUTING operations, which add (graft) or
subtract (peel) an **absolute measured density** into a **relative** claim. This is the same exon↔boundary face
§4.2 is about, reached independently.

### 2b.4 Pinning the relay: `_pin_v` semantics matter, and it works
Three variants replayed (RNA density vs oracle, mass-weighted `|Δ|`):

| | shipped | `pin=fused` (scale all three blindly) | **`pin=context`** (`_pin_v` semantics) |
|---|---|---|---|
| gdna100 ss0.50 present capON | 35.29 | 18.90 | 29.12 |
| gdna300 present capON | 27.23 | 22.14 | 25.01 |
| gdna300 present capture **OFF** | 0.176 | **0.641** ✗ | **0.068** ✓ |
| overshoot p99 / max | 31–288 / 519 | 1.00 / 1.00 | **1.47–1.58 / 1.9** |

`pin=fused` rescales all components blindly and **regresses capture-OFF 3.6×** — it forces a partial context to
account for the node's whole mass. `pin=context` applies the operator's real semantics (substitute the node's
own density for any component the context does not supply, so a partial claim stays partial), improves both
regimes, and **bounds the accumulated overshoot to the §2b.1 analytic bound** rather than to 1 — which is the
correct target, since the eff-length slack is legitimate.

## 3. The derivation

### 3.1 ✅ LANDED — the ceiling: anchor the relay context to the node's observed mass
`Σ_c ρ_c·E_c = M` at **every** node, not only at the combine — the identity of §2b.1, enforced with `_pin_v`'s
own semantics (a component the context does not supply is filled from the node's own density, so a PARTIAL
claim stays partial). It leaves the **share** untouched (composition regime intact) and forces the **scale** to
the node's observed mass (density regime enforced). No constant.

**A/B (32 conditions, `OMP_NUM_THREADS=1`):**

| | HEAD | + relay pin |
|---|---|---|
| refit=0 aggregate | 0.0964 | **0.0926** |
| refit=0 unstranded × capON | 0.1813 | **0.1720** |
| refit=0 capture OFF | 0.0474 | **0.0439** |
| refit=0 stranded ss_0.99 | 0.0390 | **0.0377** |
| refit=0 verystrong | 0.1865 | 0.1902 ⬆ |
| refit=0 gdna_none | 0.1065 | 0.1084 ⬆ |
| refit=1 aggregate | 0.0819 | **0.0779** |
| refit=1 unstranded × capON | 0.1681 | **0.1585** |

**14 better / 5 worse / 13 flat** (refit=0), **15 / 4 / 13** (refit=1). No zero-gDNA regression (unlike the
graft-frame fix of §5.1). The only losses are `verystrong` on the lowest-gDNA scenarios (gdna1/gdna5), where
capture is most extreme and the eff-length slack the pin permits is loosest — carry this into §4.2, which is
the per-channel enrichment work that should relieve it.

### 3.2 ⚠ RE-SCOPED BY §2b.1 — the pin correction as a variance term
**The original proposal (charge `(log k)²` as a scale-mismatch variance, the third sibling of `σ²_transfer` and
`b̂²`) is weaker than it first appeared, and must not be written as stated.** §2b.1 proves `k` sees composition
error only through the eff-length-weighted average, so on a contained region (`E_g ≈ E_r`) it is **blind** —
a maximally wrong composition produces `|log k| ≤ 0.036`. Charging `(log k)²` would therefore price the ROUTING
residual (graft/peel, §2b.3) and the accumulated drift, not composition mismatch, and after §3.1 the drift is
gone by construction. Before implementing anything here, settle: (a) what is left in `log k` once the relay is
pinned — measure it, it is the graft/peel identity break; (b) whether that residual is better priced as a
variance or removed by §4.2's per-channel ratios; (c) the null `Var(log k) = 1/n_dst + …` and the legitimate
`k ≠ 1` case. Original motivation retained below.

#### (original framing, for the record)
`k = M / Σ_c ρ_c E_c` is measured against a **hard observable**: the node's fragment count. When `k = 1/227`
the message's absolute claim was wrong by 227×, and today the pin silently renormalises that into a
confident-looking composition. That factor is free, already computed, and is the "no reads left to justify
this" evidence the owner describes.

**This is NOT the retired `σ²_cliff = (log r)²`.** That charged the enrichment cliff itself, which is legitimate
physics, and therefore over-damped pure-enrichment cliffs. `log k` charges only the residual the reframe
**failed to explain** — it is identically zero when the model is self-consistent, and it fired on 52–71 % of
nodes precisely because the model is not. Structurally it is the third sibling of the two terms already in the
law: `σ²_transfer` (scale sampling, M5) and `b̂²` (composition mismatch vs the self-solve, M7); `σ²_pin` would
be *scale mismatch vs the observed mass*. **Its decisive advantage over `b̂²` is that it works where
`τ_own = 0`** — i.e. on the 92 % — because it is anchored on counts, not on a self-solve the node does not have.

*Derivation needed:* the form (`(log k)²`? the DL treatment of it?), and the guard for the case where `k ≠ 1`
**legitimately** — a partial (gDNA-only) message, where `_pin_v` deliberately fills the missing components with
the node's own density so that "a seam sending gDNA only still gives `f_g < 1`". That case must not be charged.
MC-validate as we did M7 (`scripts/debug/message_variance_mc.py`), then implement, then per-condition A/B.

### 3.3 The floor — a node must be able to say "my mass is ALL RNA"
Discovered the hard way (§5.1). In a zero-gDNA library every fragment at every node is RNA, but **nothing in
prior-free pass-0 can establish that**: with no strand signal (`ss_0.50`) and no gDNA anywhere, a node has no
intrinsic evidence and ψ's uninformative Jeffreys reference drifts toward gDNA. A wrong 82× RNA over-claim in
the graft was silently doing that job; remove it and the solver hallucinates `f_g = 0.32` against an oracle of
`0.0000` over 72 % of the mass.

The exon's own observed mass **is** the evidence — but a pure composition model can never see it, because it
only ever compares shares between nodes. **So the absolute-density anchor is not only the brake on percolating
phantom RNA (§3.1/§3.2); it is also the honest source of the RNA authority that the over-claim was faking.**
The ceiling and the floor are the same missing mechanism, and §3.1–3.2 must be derived with §3.3 in view or the
fix will regress zero-gDNA libraries.

### 3.4 Bounds on justifiable depletion/enrichment
Owner's proposal, still to be worked: there may be a derivable bound on how much depletion a message can cross
and still carry confidence. At 100× depletion there are one or two reads left to split between DNA and RNA;
the precision must collapse to match, and today it does not (`σ²_transfer = Var(log r)` is *tiny* when both
nodes are count-rich — measured 0.0085 on the hop that does the damage). Candidate framing: express the
delivered claim as an implied **count** at the destination and require its precision to be consistent with the
counts that actually exist along the path (a data-processing/bottleneck bound). Measured on the anchor path,
the bottleneck is boundary 2200 at 240 unspliced fragments.

---

## 4. Two adjacent modelling gaps that bound how far §3 can go

### 4.1 ⚠ STRUCTURAL: the region/boundary map has no TSS/TES (deferred — big, and reaches to the index)
**Owner, 2026-07-25.** The region + boundary partition built from the GTF **does not annotate transcription
start / end sites**. A transcript END is therefore not represented anywhere in the region↔boundary map, so when
the solver meets one it has no way to know the RNA simply *stops* — it models the resulting density drop as a
**capture cliff**, which is wrong. It is truly the end of a transcript, not an enrichment gradient.

This is a genuine structural defect and fixing it reaches all the way back to **index building and the
accumulator**, so it is explicitly **NOT being fixed now**.

*Why the impact is expected to be small in practice:* it became visible only through a degenerate artifact of
the simulator, which generated a transcript whose exon interval coincided **exactly** with the interval of a
multi-exonic transcript — putting a transcript end and a splice junction at the identical position. That
essentially never occurs in real biology, so the sharp version of this failure is an artifact of simulated
data. It remains a real modelling gap; it is simply not the highest-yield fix.

*Where it showed up here:* exon 2201 carries a 6,701-fragment junction on one flank and only 8.4 on the other —
it is a **terminal** exon. Its RNA does not splice out on the terminal side, it stops. The peel only knows how
to subtract what splices out, so the residual is attributed to contiguous (nascent) RNA in a library that has
none. Note there is already an unused structural bool for the mature channel — `mrna_active_pos/neg`, computed
in `build_node_statics` — while the measurement stream is gated only by `fp_a/fn_a`, the *nascent* gate. That
is a partial mitigation available without the index work.

### 4.2 ⚠ THE PHYSICS IS CONFIRMED, BUT THE FIX IS REFUTED — do not build per-channel ratios for the relay

**Derivation result, 2026-07-25** (`scratchpad/derive_3_channel_ratios.py`, `derive_2_relay_pin.py`). The
owner's slope model is **true and now measured** (tables below). But the proposed remedy — replace the single
`r` with per-channel enrichment ratios — was tested at its **upper bound** and buys essentially nothing.

**The measurement.** gDNA is uniformly present along the genome, so any ratio of true gDNA densities IS a ratio
of capture efficiencies. Substituting that ORACLE capture step for the model's `r` in the (now pinned) relay is
the best any per-channel model could possibly achieve for the unspliced channel:

| condition | pinned, model `r` | pinned, **ORACLE** `r` | gain |
|---|---|---|---|
| gdna300 ss0.99 capON | 16.06 | 16.22 | **−0.16** |
| gdna300 present capON | 25.01 | 25.18 | **−0.17** |
| gdna100 ss0.50 present capON | 29.12 | 29.45 | **−0.33** |
| gdna300 present capture OFF | 0.0676 | 0.0543 | +0.013 |
| gdna100 verystrong | 43.08 | 41.26 | +1.83 (≈4 %) |
| gdna5 / gdna1 / none verystrong | 10.50 / 9.45 / 10.17 | 10.23 / 9.81 / 10.17 | +0.27 / −0.36 / 0.00 |

(mass-weighted `|Δ|` of the relayed RNA density vs oracle; replay fidelity against the shipped pinned relay is
`max|Δρ| = 0.00e+00`.)

**Why a perfect `r` buys nothing.** The reframe's only job is to convert density **scale** between frames, and
**composition is invariant under a common scale**. `r` already cancels in the pin (`k ∝ 1/r`). Once §3.1 pins
the relay at every hop, the scale is re-derived from the node's *observed mass* at each node, so `r` only has
to be roughly right — the pin absorbs the rest. **§3.1 retired §4.2's relay motivation.** Attempting a
per-channel `r` earlier, before the pin, would have looked productive and been superseded.

**What the channel distinction is STILL needed for.** Not ratio *accuracy* — frame *assignment* of the MEASURED
channels: whether the junction's spliced flux should be lifted by `r` at all. That is §5.1 (the graft-frame
fix), an 82× error at node 1909, and it remains live and blocked on §3.3. The spliced-channel half of the
owner's model is therefore real work; the unspliced-channel half is not.

**The physics, measured (retained — it is correct, and it is why §5.1 matters):**
**Owner, 2026-07-25.** Probe placement determines what happens across a boundary, and the current model has a
single `r` built from total density. Reality:

* A probe may **span an exon–exon junction**, enriching *spliced* fragments strongly — to the point that
  **spliced fragments can be denser than unspliced ones** across that boundary.
* A probe placed **mid-exon** leaves the exon's ends depleted, so by the time you reach the splice site the
  density is well below the exon's interior, and spliced fragments there are *less* dense than the exon.

Therefore **we cannot assume the boundary's density equals the exon's density on the spliced side, nor that it
equals the intron's density on the unspliced side.** The boundary genuinely sits on an enrichment/depletion
*face*, and the slope depends on probe placement. At minimum three ratios are needed:

| ratio | channel | current model |
|---|---|---|
| boundary ↔ exon, **spliced** channel | junction-crossing fragments, both blocks inside exons | assumed via the single `r` |
| boundary ↔ exon, **unspliced** channel | fragments straddling the exon↔intron seam, partly unbaited | assumed via the single `r` |
| boundary ↔ intron, unspliced channel | — | assumed via the single `r` |

**Measured support (oracle only, no solver in the loop).** Using oracle gDNA density ratios as pure capture
measurements, and oracle RNA for the spliced channel:

| condition | e(B unspliced)/e(EXON) | e(B unspliced)/e(INTRON) | e(B spliced)/e(EXON) |
|---|---|---|---|
| capture **off** (control) | 0.980 | 0.985 | 0.975 |
| capture on | **0.549** | **1.282** (p75 49.6) | 0.863 |
| capture on + nRNA | 0.467 | 1.158 (p75 45.8) | 0.543 |
| **verystrong** | **0.125** | **2113** | 0.560 |

The control is clean: with capture off all three ratios are ≈1 — no capture, no slope. Under capture the
boundary sits genuinely INTERMEDIATE and the face steepens sharply: at verystrong the exon:intron capture ratio
is ~17,000×, and the boundary sits at 0.125× the exon but 2113× the intron — i.e. **pulled strongly toward the
exon**, as expected for a seam-straddling fragment anchored on the baited side. The model's single face ratio is
roughly unbiased in the median but scattered per-edge on exon↔boundary (p25–p75 0.55–2.9 at verystrong) while
intron↔boundary is tight (0.88–0.99).

Two mechanisms were tested and did NOT explain the exon-edge error: excluding the spliced from the face
`ρ_tot` is WORSE at normal capture (−0.199, −0.106 in log) and better only at verystrong (+0.276); and
substituting the oracle capture step buys ~nothing (table above).

So the spliced channel's ratio to the exon interior is ≈1 with no capture and degrades smoothly as capture
strengthens — a **slope, graded by capture strength**, exactly as described, and measurable. (The simulator
appears not to model junction-spanning probes, so the ratio > 1 case is not represented in this suite; real
data will have it.)

Note a concrete place the single-`r` model already mixes the channels: `_rho_faces` includes the one-sided
spliced density in the acceptor face's `ρ_tot`, so that face's density is part exon-frame (the spliced part)
and part boundary-frame (the unspliced part), and one `r` is derived from the mixture.

---

## 5. What was tried and REJECTED — do not re-run

### 5.1 The graft-frame fix ("F1") — correct in isolation, NOT landable alone
The measured mature `ρ_μ = S/E_spl` is added to the source boundary's RNA **before** `·r`, so the destination
exon's enrichment step is applied to a measurement that already carries it. At node 1909: `ρ_μ = 0.1741`,
`r = 154.5`, delivered `f₊ = 0.7701` against an oracle RNA fraction of `0.0094` — an **82× over-claim**. Moving
the graft outside the reframe lands it at `f₊ = 0.0155`. §4.2's oracle table proves the premise (ratio ≈ 1, not
≈ `1/r`).

**A/B (32 conditions):** aggregate 0.0964 → **0.0789** (refit=0) and 0.0819 → **0.0687** (refit=1);
unstranded × capON 0.1813 → **0.1232**. It also flips
`test_mature_measurement_disagreement_silenced` from xfail to XPASS — a test written earlier to document
exactly this spliced↔gDNA channel coupling.

**Why it was reverted:** it regresses zero-gDNA libraries badly — `gdna_none_ss_0.50_nrna_none_capture_on`
0.0929 → 0.2693, with solved `f_g = 0.3228` against an oracle of `0.0000` over 72 % of the mass. Per-condition
16 better / 12 worse (r0), 12 / **13** (r1). The over-claim was **load-bearing**: it saturated `f_R → 1` and
produced the right answer for the wrong reason (§3.3). **Land this together with §3.3, not before.**

### 5.2 Moving the graft inside `_pin_v` ("F1b") — WORSE, hypothesis refuted
Hypothesis: the zero-gDNA regression came from adding the graft *after* the pin, removing the components'
competition for the node's observed mass. Tested: aggregate 0.0891 (vs F1's 0.0789), **9 better / 17 worse**
vs HEAD, and `gdna_none` still regressed (0.1385 vs HEAD 0.1065). The pin placement is not the mechanism.

### 5.3 The peel's frame ordering — a ≤7 % effect, not the driver
Subtracting the mature *before* vs *after* the reframe on real peel edges: 0.0000/0.0000, 1.0306/1.1135,
0.0625/0.0625. Either the mature density is small relative to the exon's RNA (so the side barely matters) or
large enough to fully consume and clip to zero either way. The peel's real defect is §4.1 — at a **terminal**
exon there is nothing to subtract, because the RNA does not splice out, it stops.

### 5.4 Rejected by the root-cause workflow (`wf_1ba425be`), with measurements
| proposal | verdict |
|---|---|
| force `r = 1` on the graft edge | worse (0.0743 → 0.0879 pooled); `r` carries real information |
| graft precision = spliced COUNT not mass | makes the message ~2.2× **stronger** — forbidden direction; +0.0020 |
| wire `graft_rna_logvar` (M2) as derived | +0.0152; `v_ν = ∞` makes the delta method return ∞ and kills every graft |
| charge `σ²_transfer` on the graft | 0.0474 → 0.0470 only, and it double-counts a common-mode `r` (M5 is proven) |
| remove the graft entirely / peel before reframe | all measured worse |
| `k = 1` (disable the pin) | big aggregate win but −0.0000 in-scope, +0.0070 on boundaries, 11/32 regress — correct diagnosis, wrong instrument (this is what §3.1/§3.2 must replace properly) |

---

## 6. Invariants any solution must preserve

* **The composition regime stays.** The enrichment-ratio model is correct and is what made the capture arm
  solvable; the density regime enters as a bound and an audit, never as a replacement.
* **Pass-0 must be WEAK and correctable.** An over-confident message that PINS a node wrong is worse than a
  weak one slightly off. A node lacking the data to solve itself must end up near-zero precision or unsolved —
  never moderate precision.
* **The `√2·σ_own` pin-safety inequality** (`p_eff = 1/max(v_msg, G²−v_own)`): a message can out-weigh a node's
  own belief only if it agrees to within `√2·σ_own`. Any change that breaks this is a regression even if the
  A/B improves.
* **No magic numbers.** Structural presence tests, derived quantities and exact limits only.
* `N` enters only as power, never as a composition vote. Variances in log-odds, never on simplex fractions.
* Counts are Poisson; the synthetic suite is Poisson by construction.

## 7. Method

Derive → MC-validate (`scripts/debug/message_variance_mc.py`, the M1–M7 arbiter, currently 0 failures) →
implement → per-condition A/B (`scripts/debug/pass0_oracle_bench.py`, refit=0 AND refit=1,
`OMP_NUM_THREADS=1`) → dissect any regression rather than assuming a theory flaw. Suite-wide node table:
`scratchpad/suite_dissect.py` → `/tmp/suite_nodes.npz`. Per-scenario: `scratchpad/dl_dissect.py`. Per-node:
`scratchpad/dump_node.py`, `scratchpad/seam_lambda_audit.py`. Goldens LAST.
