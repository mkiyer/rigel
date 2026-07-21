# Session handoff — message propagation arithmetic: the density-composition-conserving hybrid

**Branch:** `calib-ambig-init-wip`. **You are picking up the calibration message-propagation work.** Read this
prompt, then the linked docs, then the code. The north star below is the whole point of the effort.

---

## The north star

> **The message propagation system must CONSERVE DENSITY COMPOSITION across every region↔boundary transition,
> including hybrid-capture enrichment/depletion cliffs.**

The invariant is the per-component **density composition** `ρ_g : ρ_nascent± : ρ_mature±` — the capture enrichment
`e(x)` cancels in the ratios (`background_reference_derivation.md` §3: DNA rate = SCALE × SHAPE). A node's *count*
`f_g` is **frame-dependent** — a region measures *contained* fragments, a boundary measures *crossing* flux, and
gDNA vs RNA have different fragment-length distributions, so their per-component effective lengths differ between
the two frames. A message must convert frames correctly **while conserving the composition**, and it must do so
across a total-density cliff where the two nodes sit at very different enrichment.

---

## The clean wins — KEEP THESE (landed / validated)

1. **gDNA intron factory** (committed `ab7fbe69`; `gdna_intron_factory_design.md`). Peels confident gDNA from
   introns against the intergenic background NegBinom — a λ-factor (`config.intron_factory`, default ON via
   `_build_intron_prior`). Introns are off-target ⇒ `ρ_bg` is their true gDNA. Fixes the zero-gDNA false-positive
   and the gDNA under-call at unstranded introns; zero-regression on all 32 scenarios. This supplies the
   **accurate source beliefs** the message system needs at introns.

2. **The clean eff-length-frame log-odds SHIFT** (committed `d8fef478`; `cliff_message_derivation.md` §3, §7).
   `λ_dst = λ_src + log(E_g^dst/E_g^src) − log(E_r^dst/E_r^src)`. MC-validated across FL distributions
   (`scripts/debug/cliff_message_mc.py`); enrichment cancels. LANDED for CLEAN transitions only
   (intron/intergenic ↔ boundary, no mature — `bp_solver._COMPOSITION_MODE=True`, per-edge `use_shift` gated by
   `is_exon_node`).

3. **★ THE SPLICED-DENSITY SOURCE/SINK ABSORPTION — the clean win to STUDY and GENERALIZE** (`_RNA_ABSORB`,
   prior-session work, owner-confirmed CORRECT). In `bp_solver._scan`, `rho_pos`/`rho_neg` carry
   `+ SPs/esp` (a boundary SOURCE of mature: add the spliced density) and `− absorb = − SPd/ESPd` (a boundary
   SINK of mature: subtract the spliced density). This **keeps the mature (spliced) density composition the SAME
   across the splice-junction boundary ↔ exon region** — `ρ_m = spliced/E_spl` (density space, the per-face
   half-triangle `eff_spl`) equals the exon's contained-mature density (sim-verified: single-junction ratio
   ≈ 0.965; MC in `scripts/debug/cliff_exon_boundary_mc.py` reconstructs both directions exactly). This is the
   MODEL for how the whole system should conserve density composition across a cliff.

---

## The method to build — the HYBRID enrichment-corrected density (the PREFERRED, VALID mode)

This is option (3) from the prior session, and the owner's decision: **it is NOT a backdoor — it is a valid
computation for the message MODE.**

- **MODE** = an **enrichment-corrected density** that conserves density composition across the cliff. Keep the
  destination's **observed** total as the anchor (robust — it is real data, not derived from a possibly-inaccurate
  source belief), and correct the composition so the density ratio is conserved across the enrichment jump —
  exactly as the spliced absorption already does for the mature (it keeps `ρ_m` the same across the junction). The
  enrichment ratio, if needed explicitly, is already measured by the σ²_transfer projection
  (`mu_proj[dst] − mu_proj[src]`, the NPMLE enrichment landscape).
- **PRECISION** = **σ²_transfer (the enrichment-crossing damping) + the count variance** — UNCHANGED. The cliff
  enters the PRECISION honestly (reduced confidence when crossing an enrichment discontinuity), never as a biased
  mode.
- This is the preferred method precisely because it explicitly handles the "cliff" by conserving density
  composition while keeping the observed-total anchor.

---

## Why NOT the composition-normalization (÷ΣM) shift on exon edges — the finding that forces the hybrid

`_COMPOSITION_MODE` as first written normalizes by the **imputed total** `ΣM = Mg+Mp+Mn` (source-derived), which
is cliff-invariant BUT **discards the observed-total anchor `md`**. Consequence (validated this session):

- On CLEAN intron↔boundary edges it works — because the intron **factory** makes the source accurate (win #2).
- Extending it to EXON↔boundary edges (`_SHIFT_ALL_EDGES=True`, `scratchpad/validate_shift_all.py`) **REGRESSES**:
  34 flagged regressions, exon `|Δf_g|` 0.270→0.368, confidently-wrong mass 1.66→5.57, boundary 0.206→0.246 —
  because exons/boundaries have **no accurate source anchor** (no factory, unstranded), so the composition
  amplifies their inaccurate beliefs into confident-wrong mass; and even with accurate stranded sources it is
  slightly worse at capOFF because it throws away the robust `md`.
- The mature **reconciliation** (`±ρ_m`) is CORRECT (MC + the 0.965 sim check). It is the composition
  **normalization** on exon edges that fails. ⇒ keep the observed anchor: the **hybrid**.

See `[[composition_mode_regresses_post_tau]]` — this is the same architectural wall, now understood.

---

## Your task, in order

1. **Study the clean win** (win #3): read `bp_solver._scan`'s `rho_pos`/`rho_neg` (`±SPs/−absorb`), the density
   mode (`log(rho_pos/(md/erd))`), and `cliff_exon_boundary_mc.py`. Understand *exactly* how the spliced
   source/sink conserves the mature density composition across the junction while keeping the observed `md`
   anchor — this is the template.
2. **Derive the hybrid enrichment-corrected density MODE** for every message type — region↔boundary, both
   directions, all region types (intergenic / intron / exon) — conserving `ρ_g : ρ_nascent : ρ_mature` across the
   cliff, keeping the observed `md` anchor, precision = σ²_transfer + count variance. Extend the spliced
   source/sink principle to gDNA and nascent (not just mature). Be rigorous; write it up in
   `cliff_message_derivation.md` (extend §8 / add a §9).
3. **MC-validate** the hybrid across many gDNA/RNA FL distributions + enrichment levels (reuse / extend
   `cliff_message_mc.py` and `cliff_exon_boundary_mc.py`). Ensure the prototype is PERFECT before touching the
   sweep.
4. **Implement** in `bp_solver._scan` as the message MODE (the density mode with the composition-conserving
   correction; keep `_RNA_ABSORB` / the spliced source-sink), and **validate pass-0 vs oracle on all 32 cached
   `ambig_dense_10mb` scenarios** (factory + τ ON; `scratchpad/validate_cliff_shift.py` /
   `validate_intron_full.py` are the harnesses). Gate: boundary AND exon `|Δf_g|` improve, no node-type/condition
   regressions, stranded controls flat, confidently-wrong mass does not rise.
5. **Then** the gDNA hyperprior study — a broken pass-0 makes the hyperprior learn garbage, so it comes *after*
   a healthy message system.

---

## State, files, tools

- **Code:** `src/rigel/calibration/bp_solver.py` — `_scan` is the message loop. Toggles: `_COMPOSITION_MODE`
  (shift, default True), `_SHIFT_ALL_EDGES` (exon extension, default False — REGRESSES, keep OFF), `_RNA_ABSORB`
  (spliced source/sink, default True — the clean win), `_TAU_PRECISION` (default True). `use_shift` per-edge
  gate; `is_exon_node`; `rho_pos`/`rho_neg` carry the mature terms; `Mg/Mp/Mn/_den` are the composition
  normalization; the density mode is the `else` branch (`log(rho.../(md/e))`).
- **Docs (read in order):** `cliff_message_derivation.md` (the shift + §8 mature reconciliation, MC-validated) →
  `gdna_intron_factory_design.md` → `background_reference_derivation.md` §3/§8 (SCALE×SHAPE, the cliff, why
  `ρ_bg` is a floor for exons but a two-sided estimate for introns) → `CALIBRATION_ARCHITECTURE.md` §0–1
  (count-zero-information) → `CALIBRATION_STATUS.md`.
- **MC prototypes:** `scripts/debug/cliff_message_mc.py` (intron↔boundary shift), `scripts/debug/cliff_exon_boundary_mc.py`
  (exon↔boundary mature reconciliation) — both passing, the ground-truth arbiters.
- **A/B harnesses (scratchpad):** `validate_cliff_shift.py` (shift OFF/ON), `validate_shift_all.py` (exon
  extension OFF/ON), `validate_intron_full.py` (factory OFF/ON), all over the 32 cached scenarios with a
  regression detector.
- **Memories:** `cliff_message_derivation`, `composition_mode_regresses_post_tau`, `gdna_intron_factory_design`,
  `background_reference_derivation`, `phantom_and_tau_precision_verified`.

## Standing constraints (hard)

- **NO new magic numbers / constants / heuristics without pausing to discuss.** No clamps, cliffs, or binary
  cutoffs — smooth, honest Bayesian behaviour.
- Build/test/lint only inside the activated `rigel` conda env; always `OMP_NUM_THREADS=1`.
- Develop/validate on the **cached** `ambig_dense_10mb` scenarios, pass-0 vs oracle (no full pipeline needed).
- Always test AMBIG with ample single-strand nodes present.
- **The owner drives commits and golden regeneration.** Land behind a default-OFF toggle, A/B, then flip on
  owner sign-off.
