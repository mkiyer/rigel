# Two-component {RNA, gDNA} message-passing — clean redesign

**Status:** design (2026-07-08). Replaces the 5-term {mature±, nascent±, gDNA} 3-channel message layer with
the original 3-term node vector **{RNA+, RNA−, gDNA}**, plus a boundary spliced source/sink and an
intron-density-vs-background-gDNA likelihood. Goal: one elegant, streamlined algorithm — not the current
spaghetti of layered ideas. Validated against `scratchpad/toy_msg.py` (background-DNA knob + synthetic σ²) and
the 16/24 benchmarks (alignment-based oracle).

## 1. Model

Each node carries a **3-term belief `(f_pos, f_neg, f_g)`** — the deconvolution of its **unspliced** fragment
mass into sense-RNA / antisense-RNA / gDNA. RNA here is UNDIFFERENTIATED nascent-vs-mature (both are RNA vs
gDNA); the mature/nascent split is reconstructed downstream by the per-locus EM, never by calibration. This is
already what `NodeBelief` is — the redesign REMOVES the message-layer mature/nascent overlay (`nasc_p/nasc_n`
running beliefs, `mat_elig` gate, the mature sub-message split), not the belief.

**Axioms:**
- **A1. Spliced fragments are pure stranded RNA (mature).** They exist only at splice-junction boundaries, are
  a FIXED measurement, and NEVER compete with gDNA / are never partitioned. A junction is a pure-RNA **source**
  (emits mature onto the abutting exon) and **sink** (absorbs the exon's mature so it does not leak into introns).
- **A2. A fragment count carries zero intrinsic gDNA/RNA information** (`CALIBRATION_ARCHITECTURE.md`). Unspliced
  mass is deconvolved by three sources only: the strand likelihood, the background-gDNA likelihood, and cross-node
  messages.
- **A3. gDNA flows genomically-contiguously; RNA flows only where its strand is continuous** (`free_s`). Nascent
  RNA (unspliced) extends contiguously into introns; mature RNA is confined to exons (A1 keeps it there).
- **A4. An exon's unspliced mass is not locally solvable when unstranded** — gDNA and nascent are indistinguishable
  there. The exon DEFERS to messages: nascent identified at the intron (background likelihood) and propagated back;
  mature identified by the junction spliced measurement (A1).

## 2. The per-node solve (3-term ψ), by regime

The node integrates, over the 2-simplex `(f_pos, f_neg, f_g)`:
- **Strand likelihood** — the Beta-Binomial tilt of the per-strand counts (the only INTRINSIC signal; count enters
  only as overdispersed Fisher information). **Strand-specific data (κ≠½): introns and exon-intron boundaries
  deconvolve unspliced → RNA + gDNA from strand ALONE, no message needed.**
- **Background-gDNA likelihood** (gDNA-eligible nodes = introns + intergenic; see §4). The likelihood that the
  node's unspliced density is drawn from the gDNA-background distribution alone; a density in EXCESS of that
  background is RNA (nascent). This is the intron-density deconvolution to RETAIN.
- **Cross-node messages** (§3) — the imputation that carries the identified split to nodes that cannot self-solve
  (unstranded exons). **Unstranded data (κ=½): introns deconvolve via message + background; exons via message.**

## 3. Message propagation — forward + backward

The chain is a forest of linear region↔boundary paths → BP is exact in one **forward + one backward** pass. Each
node solves once both directional messages have arrived (the terminal node of a path solves immediately in its
finishing direction; interior nodes combine both).

**What flows:** the unspliced (gDNA + nascent-RNA) density, per component, at the total-density disagreement
message precision `σ²_msg = σ²_imp + 1/n_src` (the settled total-density basis — NOT per-component, NOT per-edge).
gDNA flows across all seams; per-strand RNA flows only on `free_s` edges.

**The junction spliced source/sink (A1) — the replacement for the mature gate.** At a splice-junction boundary
with measured spliced mass `S` (pure mature RNA):
- **Source (boundary → exon):** add `S` to the abutting exon's RNA (the mature the exon expresses).
- **Sink (exon → intron):** the exon's mature must NOT transfer into the intron. The mature at the exon is the
  spliced measurement scaled to the exon footprint, `M_mat = S · (E_exon-body / E_spliced-junction)` (the
  junction-crossing fraction ↔ exon-body fraction by effective length). Subtract `M_mat` from the exon→intron
  transferred density; only the **nascent** remainder flows contiguously into the intron.
  - This is the crux fix. The current lumped path subtracts only the *spliced* mass `S` (`rho_mat_dst = SPd/ESPd`),
    NOT the exon-body mature it implies — so the exon-body mature relays through the intron (the confirmed
    through-relay bleed, `g0321_msg_trace.py`). Scaling `S` to the full exon mature closes it WITHOUT a mature
    channel.

**The flow, end to end** (one transcript, forward direction — backward mirrors):
TSS exon (nonspecific unspliced mass, A4) → junction boundary (spliced sink removes `M_mat`, rest transfers) →
intron (partition transferred + own unspliced into gDNA vs nascent via the background likelihood §4) → next
junction (spliced source adds mature; partitions its own unspliced from the incoming message) → terminal exon
(combines incoming nascent + sourced mature as total RNA, plus the gDNA estimate; solves completely).

## 4. Background-gDNA distribution (the intron deconvolution)

Model the genomic-DNA density distribution over gDNA-eligible regions; the intron partition asks *"is this
node's unspliced density improbably high under gDNA-alone sampling?"* — the excess is nascent RNA.

**Training set (confirmed):**
- **Fit from INTERGENIC regions alone first** (gDNA-clean; nascent-free). Theoretically correct.
- **Iteratively add intronic regions after the first nascent partition** — introns deconvolved as low-nascent
  join the background to grow sample size; nascent-heavy introns are excluded (they would self-contaminate).
- Rationale: in real data nascent is sparse, so full intergenic+intron training is *usually* safe (the current
  behavior), but it FAILS for toys/simulations (and any high-global-nascent regime) — `toy_msg.py` shows a 5×
  intron density excess (true f_g 0.20, 80% nascent) called ~gDNA because the nascent-heavy introns inflate the
  background floor. Intergenic-first + iterative fixes this while retaining sample size where nascent is sparse.

## 5. Cleanup — keep / remove

**KEEP (validated foundation):**
- The 3-term `(f_pos, f_neg, f_g)` belief + the 2-simplex per-node ψ solve.
- Total-density disagreement message precision (`σ²_imp + 1/n_src`).
- The strand Beta-Binomial (symmetric gDNA/RNA overdispersion; unstranded-uninformative).
- The intron-density-vs-background-gDNA likelihood (retooled to intergenic-first + iterative, §4).
- The junction spliced mass (`SP/SN`) as the pure-RNA source/sink — extended to scale `S → M_mat` (§3).
- The structural index (region types, splice-junction boundaries) — `mature_eligible` may be reusable for the
  spliced source/sink placement; re-audit whether `boundaries.feather`/v6 is still needed or reverts.

**REMOVE (the 5-term overlay + the A/B scaffolding):**
- The mature/nascent MESSAGE split: `nasc_p/nasc_n` running beliefs, `mat_elig` destination gate, the
  mature sub-message, the strand-gated nascent source (`pnasc_*_loc`, `w_str`).
- The shelved per-edge `local` shrinkage (`adjacent_disagreement_local` + per-edge threading) and the
  `sweep_disagreement_pass2` A/B knob (collapse to the single total-density basis).
- The `sweep_mature_nascent_split` A/B flag (the redesign IS split-off).
- Collapse to ONE production path (retire the `sweep_disagreement_shrinkage` legacy `σ²_edge` branch too, if the
  redesign subsumes it) — no more layered flags.

**Output schema + EM: UNCHANGED** — calibration still emits per-region gDNA-vs-RNA mass; the EM reconstructs
mature/nascent from reads (n_t+1 components). Confirmed calibration-only (design workflow, 2026-07-08).

## 6. Validation

- **Toy (`scratchpad/toy_msg.py`):** background-DNA knob (`gdna_fraction` genome-wide + background-rich genome)
  + synthetic σ² (monkeypatched fit). Sweep mature × nascent × gDNA × strand × capture; the redesign must (a)
  recover intron nascent (density excess → nascent, the current failure), (b) not bleed exon mature into introns
  (the spliced sink), (c) not orphan exon nascent to gDNA (A4 propagation).
- **Benchmarks:** 16 (quick_3to1_5mb) + 24-AMBIG (ambig_dense_10mb), alignment-based calib-vs-oracle
  (`scratchpad/calib_oracle_accuracy.py` corrected) + intron-bleed + AMBIG-exon + nascent-recovery. Must beat or
  match the current 3-channel on both, with the flat-gate regime (nascent > mature) fixed.
- Regenerate goldens once landed (production behavior changes).

## 7. Open design details (to resolve during implementation)

- The `E_exon-body / E_spliced-junction` eff-length ratio for `M_mat` (§3) under capture — spliced and unspliced
  are captured differently; the ratio may need a capture correction.
- Iterative background (§4) convergence: how many intron-addition rounds; the density-excess threshold for
  "low-nascent" intron admission.
- Whether the strand-specific local deconvolution (§2, κ≠½) and the unstranded message path unify into one solve
  or branch on a κ-derived weight (prefer unified: strand precision (2κ−1)² naturally → 0 unstranded).
