> ## ⛔ SUPERSEDED — THIS IS NOT THE LIVE HANDOFF
> The live handoff is **`SESSION_2026_07_26_HANDOFF_14.md`**; the entry point is **`ROADMAP.md`**.
> This file is kept as HISTORY — what was tried, measured and refuted. Its numbers describe the
> code as it was on its own date and many are now superseded. **Do not act on it without checking
> the live handoff first.** Its own DO-NOT-RE-RUN findings (if it has any) remain HERE and still
> stand — HANDOFF_14 §3 indexes which files carry them.

# Session handoff — the composition peel LANDED, the MESSAGE PACKET landed, and P1d is one measurement from shipping

*(this file's own original header — superseded)* Date: 2026-07-26. Branch `calib-ambig-init-wip`, HEAD `fae10a3f`.
Supersedes `SESSION_2026_07_26_HANDOFF_12.md` (still the composition-peel record).

Suite **0.0855 (refit=0) / 0.0671 (refit=1)**. Gates all green.

---

## 0. Read this first — where the project actually is

**The error MASS is close to its prior-free ceiling. The CALIBRATION of that error is not.** That distinction
is the whole frame for what remains, and it is measured, not asserted:

| | session start (`g5`) | **now (`pk2`)** |
|---|---|---|
| suite mwae, refit=0 / refit=1 | 0.0885 / 0.0678 | **0.0855 / 0.0671** |
| suite error, reads | 12,344,845 | **11,928,101** |
| **confidently-wrong, reads** | 1,777,658 (14.4 %) | **1,186,552 (9.9 %)** |
| exon AMBIG confidently-wrong | 693,557 (20.1 %) | **378,285 (11.8 %)** |

**But `z2` inside the confident quartile — *is the confidence earned?* — is still 8.4 on single-strand exons,
92.1 on AMBIG exons, 15.5 overall.** Introns (1.5) are the only honest class. A node at `z2 = 9` will be
weighted **nine times too heavily** by a hyperprior fit that trusts declared precision. That is the poison
mechanism, and it is what remains.

`SESSION_2026_07_25_HANDOFF_10.md` §3 already measured 67–75 % of the error in the dominant arms as
premise-limited — reachable only by a better *mode*, with ×10 on every precision moving nothing. **So chasing
mwae is nearly out of road; chasing `z2` is not.** Everything on the list below is chosen for its effect on
`z2`.

## 1. Documents, in reading order

1. **`ROADMAP.md`** — the entry point.
2. **`PASS0_FINISH_PLAN.md`** — ⭐ **the prioritized work list (P0–P6) with live status.** Read its §0b.
3. **`variance_ledger.md`** — ⭐ **the standing audit of every variance term**: what each prices, the proof
   that none of them compound, and P1d's full derivation + measurement + estimator. **Update it whenever a
   variance term is added, removed or re-scoped.**
4. **THIS FILE** — state, what landed, how to run everything.
5. `weighted_rescale_design.md` — the M12 conservation rescale and, in §9, the MESSAGE PACKET that subsumed it.
6. `SESSION_2026_07_26_HANDOFF_12.md` — the composition peel (M11 + the three-way level fuse).
7. `SESSION_2026_07_25_HANDOFF_10.md` — the boundary evidence base; §3 the premise-limited ceiling, §12 the
   grid/precision answers.
8. `CALIBRATION_ARCHITECTURE.md` — count-zero-information, the three sources.
9. `docs/calibration/boundary_work_notes.md` — **the owner's own notes file. Do not overwrite.**

**Never reference `docs/calibration/archive/`.**

## 2. What landed this session (11 commits, `d06f8055..fae10a3f`)

**A. The composition peel — steps 1–5 of `peel_by_composition_plan.md`, and the `mrna_active` gate is GONE.**
The standing blocker ("there is no LEVEL at the seams that matter") is resolved by **M11 `residual_level`** —
the seam's RNA level from the node's own observed mass closed against an imputed gDNA density. Two structural
corrections made it work: **there is no level PRECEDENCE** (three estimators of one quantity FUSE by inverse
variance, so the "no evidence ⇒ silence" arm that wrecked low-gDNA libraries does not exist), and **the fuse
is done in LINEAR density space** (a log fuse of positive modes cannot reach zero, so it cannot hear a
factory-solved intron saying "essentially no RNA here"; ablating that estimator costs 0.0010).

**B. M12 `conservation_rescale`** — derived, MC-free but 7-closed-form-tested, landed at the combine, then
**superseded**: once λ is fused correctly (C) it is completely inert (0/0/32). The law stays in
`enrichment_frame`; the default path no longer calls it.

**C. ⭐ THE MESSAGE PACKET.** A message must let the destination do three independent things — set each
component's **LEVEL**, set the **SPLIT** (λ, precision τ), set the **TILT** (θ). The solver was reading the
split and the tilt back off the *fused densities*, so both arrived weighted by the **level** precisions
instead of their own. Because λ and θ are scale-free, a message can state them directly. Fixing the weighting
is worth **0.0889 → 0.0855 / 0.0686 → 0.0671** and is what drove the 33 % drop in confidently-wrong mass.

**D. P0 — the debug loop's loss function.** `pass0_error_table.py` now orders by **confidently-wrong mass**
and prints `CWRONG`, **`z2`** (the calibration) and `%nodeQ1` (the selection confound).

**E. P2 REFUTED.** Introns' 91 % `errQ1conf` is a *selection artifact* — 48.4 % of intron nodes are in the
confident quartile because the factory legitimately makes them confident, and their `z2|Q1` is **1.73**, the
best in the solver. No fix to make. (The `intron_prior`/`tau_lam` two-path inconsistency survives as a
correctness nit only.)

**F. P1 root-caused (twice — the first diagnosis is corrected in place).** Not the spliced measurement's mode:
that is *exact* (36.6453 claimed vs 36.6734 true). The gDNA claim is 47× too big, the two together claim
1.75× the fragments the node observed, and the rescale divided **both** by 1.75 — taking 43 % off a
measurement to make room for a guess.

## 3. ▶ THE NEXT TASK — P1d, and it is one measurement from shipping

Full derivation, measurement and estimator: **`variance_ledger.md` §3–5.** In brief:

**The gap.** The graft hands an exon a spliced density measured over a ~100 bp junction window as a claim
about a ~2,100 bp exon — a 12–21× extrapolation that nothing in the ledger prices. M8 comes closest but
assumes capture is the only reason a junction and an exon differ; **off-capture `r = 1` and M8 charges exactly
zero**, while a point and a region still differ. Measured at capture-OFF: the graft share φ has median 1.11
(small bias) but **sd(log φ) = 0.76 ⇒ variance 0.58**, against a count variance the graft is charged of as
little as 0.0005. **~1000× over-confident.**

**The prize** (probe `RIGEL_XVAR=0.3`, oracle-free, default off): confidently-wrong **1,186,552 → 689,639
(−42 %)** and **exon-single `z2|Q1` 8.4 → 1.2 — honest** — for 2.4 % more error mass (mwae 0.0855 → 0.0875).

**The estimator EXISTS** (§5.1). An exon's two flanking junctions are essentially independent, and
`Var(gap)/2` recovers `Var(log φ)` to within 9–13 %. The method-of-moments form the calibrator already uses
twice gives `ω_graft = max(0, Var(gap) − E[1/n_l + 1/n_r])/2` = **0.594 / 0.601** against a true
**0.583 / 0.574** — within 2 %, prior-free, no constant.
⚠ An earlier note in this project called that estimator refuted. **That was a statistical error** (median |gap|
compared against a full sd on a tail-dominated distribution) and is corrected in §5.1.

**What is OPEN: the SHAPE, not the magnitude.** `Var(log φ)` spans **40×** across the junction-count range
(excess 1.674 below 30 counts, 0.041 above 1000) and Poisson explains almost none of it. A pooled ω would
over-charge well-counted junctions 14× and under-charge sparse ones 3×.

> **Working hypothesis, and the owner is the right person to confirm or kill it: the junction's count is a
> proxy for its SHARE of the exon's RNA.** A junction carrying most of the exon's transcripts is
> representative (φ → 1, tight); a minor-isoform junction represents little (φ ≪ 1, and `log φ` of a small
> share is volatile). If so, the per-junction observable is the junction's flux **relative to the exon's total
> RNA flux**, for which the exon's *other* junction's flux is the prior-free stand-in.

**Three questions put to the owner and not yet answered:** (1) is "low count ⇒ minor isoform ⇒
unrepresentative" the causal story, or is count merely tracking expression (a confound)? (2) is there a better
per-junction reliability signal — junction usage fraction, distance to a transcript end, isoform count? (3)
does ~2× irreducible scatter match real data? ⚠ The suite is Poisson by construction
(`synthetic_suite_is_poisson_omega_zero`), so anything overdispersion-shaped is **understated** here.

## 4. The rest of the list (`PASS0_FINISH_PLAN.md` has the detail)

| # | item | why |
|---|---|---|
| **P1d** | the graft's extrapolation variance | **NEXT** — derived, sized, estimator found; the shape is open |
| **P1e** | the conservation **SURPRISE** `z² = δ²/(αᵀΣα)` | owner-directed follow-up. The only damping term that is **never inert** — M7's DL needs `τ_own > 0` and so switches off on all unstranded data (84 % of error mass); the node's own mass is always known |
| P3 | AMBIG exon over-confidence | `z2\|Q1` still **92.1**, 378 k confidently-wrong. Down 45 % this session as a side-effect of the packet |
| P4 | the FAR level estimator is a **LOOKUP, not a message** | a BP correctness violation: the left message's share depends on the RIGHT neighbour's belief, so the two messages share information and the fused belief is *structurally* over-confident. Load-bearing (ablating costs 0.0010) — replace, do not delete |
| P5 | `r` from the gDNA channel | the "moot" verdict is stale — M11 is the first consumer for which `r` does not cancel |
| P6 | fragment length as a FOURTH source | **gated**: verify no circularity in `fl.py`'s pmfs, and get a realistic separability estimate first |
| — | then the **gDNA hyperprior refit**, then the re-solve | blocked on the above |

Also open: `CalibrationConfig`'s docstring still describes the **retired** NPMLE disagreement-variance
(`DensityNPMLE.project` is called nowhere). Fold the fix into whatever touches it next.

## 5. ⛔ DO NOT RE-RUN — settled this session

| item | verdict |
|---|---|
| introns as a trust defect (P2) | **selection artifact.** `z2\|Q1` = 1.73, the best class in the solver |
| ablating the RNA measurement channel (`RIGEL_RNAMEAS_OFF`) | 0.0855-era arm: 0.0895 → 0.1033, 4 better / 17 worse. It is the only thing that lets a zero-gDNA library say "my mass is all RNA" |
| the conservation rescale at the combine, **common-factor limit** | 0.0971 / 0.0761 — WORSE. A common rescale cannot change a composition, but it *does* move the level channels, and that is what breaks `none_ss0.50_nrna_none` (0.3624 → 0.5118) |
| the conservation rescale on all three fused densities (weighted) | same level damage, identical zero-gDNA numbers — proving the harm is touching the levels, not the apportionment |
| the conservation rescale **after** the packet | completely inert, 0 better / 0 worse / 32 flat |
| the per-message conservation rescale for λ | contributes nothing once λ is τ-fused (`pkt` ≡ `pktL`) |
| the per-message rescale on the LEVELS (`RIGEL_M12_MSG`) | net loss, but **halves capture-OFF at refit=1** (0.0289 → 0.0144, 8/2) and wrecks capture-ON. Needs two claims carried per message — its own derivation |
| the "node-frame" correction of the message's gDNA claim before M11 | conceptually wrong (the face selection makes the RATIO a clean capture step, not the units) and measured worse |
| validating a message variance against a posterior's own spread | wrong object. Validate against the **estimator's error vs truth**, as M1–M8 are |

Everything in `SESSION_2026_07_25_HANDOFF_11.md` §8 and HANDOFF_10 §4/§9.2 still holds.

## 6. HOW TO RUN EVERYTHING

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1                      # REQUIRED for A/B determinism
cd /Users/mkiyer/proj/rigel
```

### 6.1 Gates — all must stay green

```bash
python -m pytest tests/calibration tests/native -q   # fast loop (~7 s)
python -m pytest tests/ -q                           # 1248 pass, 2 xfail, 2 xpass
ruff check src/ tests/ scripts/                      # clean
python scripts/debug/message_variance_mc.py          # 0 failures over M1–M11
python -m pytest tests/ --update-golden -q           # goldens — LAST, after the solver is final
```

### 6.2 The A/B — 32 conditions, cached, ~1 s/scenario

```bash
P0_REFIT=0 python scripts/debug/pass0_oracle_bench.py --arm NAME_r0
P0_REFIT=1 python scripts/debug/pass0_oracle_bench.py --arm NAME_r1
python scratchpad/ab.py NEW BASE [extra arms...]     # strata + per-condition, both refits
```

Rows accumulate in `/tmp/pass0_oracle_bench.tsv`. **⚠ If `/tmp` has been cleared, the baseline arms are gone
— re-establish them by running the current HEAD as `pk2_r0`/`pk2_r1` before comparing anything.** Key arms:
`g5` = the pre-session HEAD (0.0885/0.0678), **`pk2` = HEAD now (0.0855/0.0671)**.

`scratchpad/ab.py` calls out the strata the owner asks for: `stranded ss_0.99`, `unstranded × capON`,
`verystrong`, `gdna_none`, `low gDNA`, `nrna_none`, capture OFF/ON.

### 6.3 The TRUST view — this is the metric that decides pass-0

```bash
python scripts/debug/pass0_error_table.py --refit 0        # writes /tmp/pass0_state.npz
```

Ordered by **confidently-wrong mass**. Columns: `CWRONG` (absolute reads), `errQ1conf` (share), **`z2`** (the
calibration: 1.0 honest, >1 over-confident), `%nodeQ1` (the selection confound — always read it beside
`errQ1conf`). The sharp per-class test (`z2` restricted to the confident quartile) is a few lines over the
saved `.npz`; see §0 for the current values.

### 6.4 Diagnostics (all default-OFF, all verified bit-identical when unset)

| flag | what it ablates |
|---|---|
| `RIGEL_S2T_OFF=1` | both cliff terms (M5 + M8) |
| `RIGEL_RNAMEAS_OFF=1` | the RNA measurement ψ factor |
| `RIGEL_M12_MSG=1` | the per-message pin → the weighted rescale (halves capture-OFF at refit=1, wrecks capture-ON) |
| `RIGEL_XVAR=<v>` | **P1d's probe** — a flat extra variance on the grafted spliced component |

### 6.5 `scratchpad/` (untracked — **keep it**)

`ab.py` (the multi-arm A/B reader) · `w4_provenance.py` (which level estimator sets each seam, from the
solver's own inert `_capture["_lvl"]`) · `p1_worked.py` (**one node in fragments and base pairs — start here
when something looks wrong**) · `p1_capoff.py` (the paired capture-OFF/ON dissection) · `wp1_prototype.py`
(the offline rescale prototype, on the real pre-pin state) · **`x1_graft.py` (P1d's graft share — the script
the next task extends)** · `e1_enrichment.py` / `e2_edges.py` (the enrichment-ratio interrogation, incl. the
fragment-length Fisher measurement) · `w1_level.py` / `w2_fuse.py` (the level-estimator studies) ·
`m11_level_fuse.patch`, `steps_1_to_4.patch` (restore points).

## 7. Invariants — preserve these

* the **`√2·σ_own`** DL pin-safety inequality (`mismatch_deflate` untouched all session);
* `Σ_c ρ_c·E_c = M` enforced with `_pin_v`'s partial-claim semantics (a component the message does not
  **supply** is filled from the node's OWN density and does not move — "supplied" is a statement about
  PRECISION, never about the density's value);
* `N` enters only as power — M11's `k` is literally an effective fragment count;
* the composition is ONE dof (λ); the tilt θ is a SEPARATE dof — **and each is now fused by its own
  precision**, which is what the packet fixed;
* **no magic numbers** — structural presence tests, derived quantities, exact limits, or parameters *fitted
  from the data* by method of moments (κ and both strand overdispersions are the precedent; P1d's `ω_graft`
  would be the third). A fitted parameter is not a tuned constant;
* **pass-0 must be WEAK and CORRECTABLE** — and this is now a tested property: M11's MC asserts honesty where
  the level is a measurement and **conservatism where it is marginal**, as two separate assertions;
* **before adding any variance term, name the physical event it prices and show which existing term would
  otherwise have absorbed it** (`variance_ledger.md` §2).

## 8. Vocabulary (owner's, use it)

RNA is **one species**. "Mature" and "nascent" are bookkeeping — only **SPLICED vs UNSPLICED** is observable.
A boundary can be an exon↔exon boundary that is *also* a splice junction, with RNA contiguous across it while
other RNA splices in or out. Both at once.
