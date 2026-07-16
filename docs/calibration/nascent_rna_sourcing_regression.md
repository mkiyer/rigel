# The nascent-RNA sourcing regression — introns lose RNA support where the strand can't deconvolve

**Status:** KNOWN, ACCEPTED regression. Deferred to a dedicated future session. Written 2026-07-16, branch
`calib-ambig-init-wip`.
**Cause:** the mature-crossing gate ([`mature_crossing_gate.md`](mature_crossing_gate.md)) — correctly — stops an
exon from sending its RNA into a flanking intron. Where that (wrong) message was an intron's *only* RNA support,
the intron now over-calls gDNA.
**Owner decision (2026-07-16):** *"allow gDNA to siphon some of the nascent RNA for now; come back to fix it with
a dedicated session geared toward correctly sourcing nascent RNA from introns and exon–intron boundaries."*

---

## 1. The regression, measured

Production `calibrate` (fitted prior) vs oracle over the 24-condition `ambig_dense_10mb` cache, gate ON vs OFF:

* The gate is **net-positive** — gDNA→RNA leak **−7.0%** (16.86M → 15.68M fragments), mean `mwae_fg`
  0.175 → 0.168, 12/24 conditions improved.
* But **8/24 conditions regress**, and the regression is concentrated in one structural place: **`nrna_present`
  × unstranded (ss 0.50)** introns, which the gate pushes *toward gDNA*.

Per-node (prior-free, oracle-truthed, `scripts/debug/node_error_report.py`), the gate-worsened introns look like
this — true composition ~⅓ gDNA, believed correctly before the gate, over-called gDNA after:

| true f_g | f_g before gate | f_g after gate | eff_rna | class |
|---|---|---|---|---|
| 0.314 | 0.339 | **0.568** | large | intron (nrna_present, capOFF) |
| 0.384 | 0.432 | **0.727** | large | intron |
| 0.332 | 0.375 | **0.568** | large | intron |

The recurring `0.568` is the intron falling back to a near-prior value once its RNA message is removed.

## 2. Why it happens — a real information gap, not a bug

A `nrna_present` intron genuinely contains nascent RNA. To be called RNA it needs a *source of information* that
says so. Today there are only three, and in this regime all three are silent:

1. **Its own strand likelihood.** Under an **unstranded** library (κ ≈ ½) the per-strand counts carry **zero**
   gDNA-vs-RNA information (the count-zero-information principle — `CALIBRATION_ARCHITECTURE.md`). So the intron
   cannot tell from its own reads that it is RNA.
2. **A neighbouring intron / boundary nascent relay.** This is allowed by the gate, but on these AMBIG loci the
   adjacent nodes are equally strand-blind, so the relay carries little.
3. **The exon's message** — which, before the gate, leaked the exon's (mature-dominated) RNA density into the
   intron. It was **wrong** (an exon's mature does not cross into an intron — the whole point of the gate) but it
   happened to push the intron in the *right direction* (toward RNA). The gate removes it.

So the gate correctly deletes a **wrong** source, and the intron is left with **no correct source** for its
nascent in the one regime (unstranded) where it most needs one. The composition is honest ("I have no evidence I
am RNA, so I default toward the gDNA prior") — the deconvolution is simply under-determined. **This is the
predicted §4.4 muting from [`boundary_spliced_channel_design.md`](boundary_spliced_channel_design.md), now
observed end-to-end.**

## 3. What the fix requires — a nascent-RNA "factory"

The intron needs a *correct* source of nascent evidence. Two exist in principle; only the first is implemented.

### 3.1 Strand-specific libraries — implemented, sufficient where present
When the library is stranded (ss ≈ 0.99), the per-strand counts DO carry gDNA-vs-RNA information, and the
intron's own strand likelihood resolves its nascent. **The regression is small-to-absent on the stranded
conditions** — it is the *unstranded* cells that fail. So for stranded data the machinery is already adequate;
the gap is specifically the unstranded case.

### 3.2 Density-discrepancy nascent sourcing — NOT implemented, and nontrivial
On an unstranded library the only remaining signal is **density**: an intron whose total unspliced density
(reads/bp) exceeds the local gДNA baseline has *excess* mass that must be nascent RNA — even though the strand
cannot say which strand it is on. Concretely: `nascent_density ≈ max(0, ρ_unspliced − ρ_gdna_baseline)`, sourced
from:
* **introns** (all unspliced excess over the gDNA baseline is nascent — introns hold no mature), and
* **exon–intron boundaries** (the crossing unspliced excess, once mature is accounted for).

**This is the missing piece: a "nascent factory" that manufactures a nascent-RNA density message from a density
discrepancy alone, with no strand information.** It is nontrivial because:

* **The gDNA baseline is itself uncertain**, especially under capture (enrichment inflates the local baseline —
  the very confound that breaks the prior; see `calibration_roadmap.md`). Subtracting an uncertain baseline from
  a noisy density to estimate a small residual is exactly the ill-conditioned differencing this project keeps
  running into.
* **It must respect transcript structure.** Nascent is present wherever a transcript is (exon OR intron), so the
  factory sources from introns and exon–intron seams — but it must NOT double-count the mature already handled by
  the spliced channel, and it must route the resulting nascent so it propagates intron→exon (the direction the
  mature-crossing gate deliberately keeps open), never exon→intron.
* **Strand attribution is genuinely ambiguous** on AMBIG loci — the density excess could be nascent on either
  overlapping strand. The factory produces a *strand-agnostic* nascent density; how it splits across the two
  admissible strands (or stays pooled) is an open design question tied to the AMBIG 2-DOF solve.
* **It interacts with the pie / DOF model.** A nascent-density source is another message competing on the same
  `(λ, θ)` axes — its precision must be set coherently with everything else, which is exactly what
  [`dof_pie_model_fix.md`](dof_pie_model_fix.md) is about. The factory likely has to land *after* the pie model.

## 4. Scope of the dedicated session

The future "nascent sourcing" session should:
1. **Design the density-discrepancy nascent factory** — the estimator `max(0, ρ_unspliced − ρ_gdna_baseline)` at
   introns and exon–intron boundaries, its uncertainty (the baseline is the hard part), and its precision in the
   message layer.
2. **Decide the strand split** on AMBIG loci (pooled vs. tilt-split), in coordination with the 2-DOF solve.
3. **Guarantee non-double-counting** with the spliced/mature measurement channel (priority #3,
   `boundary_spliced_channel_design.md`) and with the mature-crossing gate's propagation rules.
4. **Re-measure the 24-condition regression** — the target is the `nrna_present × unstranded` introns recovering
   toward their true ~⅓ gDNA without re-introducing the exon→intron leak the gate removed.
5. **Sequence it after** the pie/DOF model and the reference/prior roadmap work, since the factory's precision
   and strand split depend on both.

## 5. Why it is safe to defer

* The gate is **net-positive today** (−7% leak); the regression is a **subset** of conditions and is bounded (the
  intron falls back to a prior-ish value, it does not blow up).
* On **stranded** libraries the intron's own strand likelihood already sources its nascent — the gap is specific
  to unstranded AMBIG data.
* The honest failure mode (gDNA siphons some nascent where nothing can adjudicate) is **strictly better than the
  pre-gate behaviour**, where an exon manufactured *wholesale* phantom nascent into every flanking intron.

**In one line:** the mature-crossing gate closed a wrong door; the nascent factory is the right door, and it is
not yet built.

## 6. Pointers

* The regression measurement + tool: `scripts/debug/node_error_report.py` (per-node before/after/oracle) and the
  24-condition production A/B (see `mature_crossing_gate.md` §Phase 6).
* The predicted muting: `boundary_spliced_channel_design.md` §4.4.
* The information principle: `CALIBRATION_ARCHITECTURE.md` (count-zero-information).
* The precision/DOF dependency: `dof_pie_model_fix.md`.
