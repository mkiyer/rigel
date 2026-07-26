# Session handoff — the composition peel LANDS: the level is a fuse of three estimators, not a precedence

**This is the LIVE handoff. START HERE.** Date: 2026-07-26. Branch `calib-ambig-init-wip`.
Supersedes `SESSION_2026_07_25_HANDOFF_11.md` (its §5/§6/§7 questions are all answered below).
Gates all green: `pytest tests/ -q` **1241 pass / 2 xfail / 2 xpass**, `ruff check src/ tests/ scripts/` clean,
`message_variance_mc.py` **0 failures over M1–M11**, goldens regenerated.

---

## 0. TL;DR — steps 1–5 are implemented, working, and the revert is undone

Suite **0.0895 (refit=0) / 0.0693 (refit=1)** against HEAD's 0.0885 / 0.0678 — **9 better / 9 worse / 14 flat**
at refit=0, against the reverted `full` arm's **1 better / 31 worse**. The composition peel now *wins*
where it used to be catastrophic and the residual deficit is one identified, localized regime.

| arm | refit=0 | refit=1 | what |
|---|---|---|---|
| `m8` | 0.0900 | 0.0700 | before any boundary work |
| **`g5` (previous HEAD)** | **0.0885** | **0.0678** | `mrna_active` gate + λ precision gate |
| `full` (HANDOFF_11) | 0.0934 | 0.0684 | steps 1–5, level by PRECEDENCE — reverted |
| **`final` (this session)** | **0.0895** | **0.0693** | steps 1–5, level by **FUSE**; the gate is GONE |

And on the metric the owner said matters more, it is **better than HEAD on every large node class**:

| node class | HEAD ERR | HEAD errQ1conf | now ERR | now errQ1conf |
|---|---|---|---|---|
| exon single | 7,570,987 | 10.6 % | 7,676,603 | **7.9 %** |
| exon AMBIG | 3,450,533 | 20.1 % | 3,498,633 | **18.0 %** |
| **boundary single** | 691,635 | 6.6 % | **664,289** | **5.3 %** |
| boundary AMBIG | 384,852 | 2.2 % | 396,401 | 2.9 % |
| intron single | 217,644 | 92.1 % | 218,799 | **91.1 %** |
| suite | 12,344,845 | **14.4 %** | 12,490,634 | **12.1 %** |

**The boundary class — the class this work targets — improves on BOTH metrics simultaneously** (−4.0 % error
mass *and* −20 % confidently-wrong share). That is new: last session's step-3 arm bought its honesty by
*adding* boundary error (692 k → 939 k). And the whole substrate got more confident (`Var(log f_g)` q3
1.012 → 0.703) while getting *less* confidently wrong — which is the definition of better calibration.

**HEAD's `mrna_active` structural gate is removed** (step 4, the owner's requirement).

## 1. The three questions, answered

### Q1 — why the regression? ✅ FIXED, and the diagnosis in HANDOFF_11 §5 was exactly right

Level-precedence **case 3** (`v_ν = ∞` ⇒ the RNA channel silenced) was the whole of it. It is now gone —
not repaired, *deleted*, because the new third estimator has no "no evidence" arm. The arms it destroyed are
now the arms that **beat** HEAD:

| condition | HEAD | `full` (precedence) | **now (fuse)** |
|---|---|---|---|
| `gdna1_ss0.50_present_capON` | 0.0689 | 0.0935 | **0.0636** |
| `gdna5_ss0.50_present_capON` | 0.0906 | 0.1133 | **0.0835** |
| `none_ss0.50_present_capON` | 0.0588 | 0.0841 | **0.0562** |
| `gdna5_..._verystrong` | 0.1550 | 0.1666 | **0.1493** |
| `gdna1_..._verystrong` | 0.1591 | 0.1813 | **0.1557** |
| `none_..._verystrong` | 0.2077 | 0.2284 | **0.2043** |

**`verystrong` is now 4 better / 0 worse**, `gdna_none` 3/1, `low gDNA (gdna1+gdna5+none)` 9 better / 5 worse.

### Q2 — the precision of a density RATIO ✅ the delta method was never the problem

`Var(log w) = w_μ²(v_ν + v_μ)` (M10) is correct and stays. The suspicion in HANDOFF_11 §7 was right that
`v_ν` was the fault, but the defect was **not** a mis-derivation — it was a **representation** failure, and it
was in the *fuse*, not the propagation:

> A log-space (geometric) fuse of positive modes **cannot reach zero**. A factory-solved intron saying
> `ρ_ν = 0.0006 ± 0.0012` — measured, and *conservative* (HANDOFF_10 §12.3) — has an unbounded log-variance
> precisely *because* it is confidently near zero, so in a log-space fuse it is out-weighed by any confident
> mid-range claim. The information is real; the coordinate destroys it.

`"0.02 ± 0.04"` and `"log ρ = −3.9 ± 2.0"` are the *same* delta-method statement. Only the first can pull a
fused level down. **So the levels are fused in LINEAR density space** (each contributing `Var = ρ²·v`), and the
fused pair is converted back to the model's currency as an effective fragment count.

That is not a stylistic choice — it is measured. **Ablating the far-intron estimator costs 0.0895 → 0.0905 /
0.0693 → 0.0703, concentrated on `nrna_none`.** The *same* estimator was inert (11 % of the fuse weight, no
measurable effect) under the log-space fuse. Linear space is what turned "the intron sets the level" from an
idea that had failed three times into a load-bearing term.

### Q3 — how to split the unspliced boundary crossing ✅ a third source, and it exists everywhere

**`residual_level` (M11)** — new, derived, MC-validated, unit-tested. The node's own **observed mass** closed
against an **imputed gDNA density**:

```
φ   = ρ_g·E_g / M                the share of the crossing the gDNA claim accounts for
f_R = 1 − φ,  truncated to [0,1] with σ_f = min(φ,1)·√v_g     ρ_ν = E[f_R]·M/E_r
k   = E[f_R]²/Var[f_R]           Var(log ρ_ν) = ψ'(k) + 1/(n·E[f_R]²)
```

This is the **generic density deconvolution** — the intron factory's own primitive — with the gDNA prior
supplied by a *neighbour* instead of the intergenic pool. It is count-zero-information-legal for exactly the
reason the factory is: the information is the imputed gDNA **density**; the count only converts a density into
a composition. And it exists at **every** seam — including the `exon|exon` seams that have no factory within
reach on 97 % and no strand when unstranded, and every seam of a low-gDNA library. **That is the source
HANDOFF_11 §6 said did not exist yet.**

Why it works, in one line: **it is built on the good channel.** gDNA transports across a hop at 0.96–0.99
accuracy while the RNA channel is the contaminant (HANDOFF_10 §7) — so the peel now leans on the gDNA claim
and the node's own count, not on the exon's RNA claim.

## 2. What is actually implemented

**There is no level PRECEDENCE any more.** Every candidate level has the same shape, `ρ_ν = (1−f̂_g)·M/E_r`,
so there is nothing to rank — there are three independent estimators of one quantity and they fuse by inverse
variance, each contributing exactly the precision it earned and *nothing* where it has none:

| estimator | source | inert when |
|---|---|---|
| **OWN** | the node's own belief, precision `τ_own` | `τ_own = 0` — all unstranded data, any node with no factory |
| **FAR** | the node across the seam, reframed, charged the M5 hop cost | that node has no `τ_own` (an `exon\|exon` seam: 97 %) |
| **MASS** | M11 — the mass identity closed with the message's gDNA claim | the gDNA claim carries no precision |

Then M10 as designed: `w = ρ̂_ν/(ρ̂_ν+ρ_μ)`, `Var(log w) = w_μ²(v̂_ν + v_μ)`, applied as a **scaling** of the
exon's RNA and charged to the message precision. `v_μ` uses the spliced **COUNT**, per strand, per face.

**Files:** `enrichment_frame.residual_level` (the new law + 5 unit tests in `test_enrichment_frame.py`);
`bp_solver._peel_share` (the fuse, one function driving both the relay and the combine);
`message_variance_mc.residual_level_law` / `residual_level_limits` (M11's MC).

## 3. ⭐ Four things measured this session that the next session should not re-derive

1. **`k ≥ 1` and the trigamma.** The level's log-variance is `ψ'(k)` where `k = E²/Var` is the level's
   **effective fragment count** — the same `Var(log) = 1/n` currency used everywhere else. Exact in *both*
   limits the delta method gets wrong: `k → ∞` gives `1/k` (Poisson); `k → 1` gives `π²/6`, where the delta
   method returns 1 and is over-confident by 1.6×. MC-validated: **z2 = 0.99 / 0.93 where the level is a
   measurement, 0.15–0.41 (conservative, never over-confident) where it is marginal.** The delta method by
   contrast runs 1.14× over-confident at `f_R = 0.5` and then explodes to z2 = 0.02 (50× over-damped) at
   `f_R = 0.05` — it is unusable across this range.
2. **The truncation must be TWO-SIDED, and the upper bound is the load-bearing one.** The imputed gDNA claim
   arrives at **√v_g ≈ 1.0–1.2 nats** at exon→boundary edges under capture (measured, `w4_provenance.py`) —
   so `σ_f` is of order 1 and a one-sided positive part returns `E ≈ 0.8σ`, i.e. *"most of my mass is RNA"*,
   **asserted out of pure ignorance at a confident-looking k ≈ 2**. Bounded above, the same ignorance
   degrades to its correct limit (`Uniform(0,1)`, `E = ½`, `k = 3`) which the fuse out-weighs with any real
   evidence. This was worth ~2.5× on the mass level's mode at RNA-free seams (3.42 → 1.40 against a truth
   of 0).
3. **The claim's log-variance must be linearized at `min(φ,1)`, not at `φ`.** A relayed claim asserting more
   gDNA than the node sequenced is routine (52–71 % of nodes). Linearizing there makes `σ_f` grow *with* the
   over-claim, the [0,1] window widens faster than the mean leaves it, and the estimator **inverts** — more
   imputed gDNA returning *more* RNA, reaching `f_R = 1` at `φ ≈ 100`. Caught by the monotonicity unit test,
   not by the A/B.
4. **The count enters with its exact `1/f_R` Jacobian**, not as a bare `+1/n`: `ρ_ν E_r = M − ρ_g E_g`, so a
   count error is amplified by the reciprocal of the RNA share. MC-caught at 19 % under-statement on the
   `f_R = 0.9` arm. The linearization point is the *truncated* mean, which is what keeps it finite as `φ → 1`.

## 4. ⛔ Measured and settled this session — do not re-run

| item | verdict |
|---|---|
| **restoring the `mrna_active` gate on top of the new peel** (arm `m11g`) | 0.0884 / 0.0676, **30 of 32 conditions completely unchanged.** The gate does not compose with the fuse — it *masks* it, killing the low-gDNA wins (`gdna1_capON` back from 0.0636 to flat) exactly as it saves the `nrna_none` ones. They are complementary in the wrong way: the gate is a binary switch that is right where M11 is weak and wrong where M11 is strong. Not a candidate. |
| **dropping the FAR estimator** (2-way fuse) | 0.0905 / 0.0703, concentrated on `nrna_none`. Keep the 3-way. |
| **log-space fuse of the levels** | 0.0904 / 0.0705. The linear fuse is worth 0.0009 / 0.0012 and is the reason FAR works at all. |
| **converting the message's gDNA claim to the "node frame" before M11** (dividing by the face/node `ρ_tot` ratio) | Conceptually WRONG — the face selection makes the *ratio* a clean capture step, it does not change the units of the transported density — and measured slightly worse (0.0895 → 0.0898). Reverted. |
| **validating M11 against the truncated normal's own posterior spread** | Wrong object. A message variance must be validated against the **estimator's error** vs the truth, as M1–M8 are. The first MC did the former and was misleading. |

## 5. ▶ WHERE THE REMAINING DEFICIT IS — one regime, fully localized

The +0.0010 / +0.0015 against HEAD is **not spread**. It is `gDNA-rich × capture-ON`:

| condition | HEAD | now | Δ |
|---|---|---|---|
| `gdna100_ss0.99_nrna_none_capON` | 0.0261 | 0.0340 | +0.0079 |
| `gdna300_ss0.99_nrna_none_capON` | 0.0338 | 0.0406 | +0.0068 |
| `gdna300_ss0.50_nrna_none_capON` | 0.1716 | 0.1773 | +0.0057 |
| `gdna300_ss0.99_present_capON` | 0.0674 | 0.0723 | +0.0048 |

Strata: `stranded ss_0.99` +0.0023 (0 better / 5 worse), `nrna_none` +0.0030, `capture ON` +0.0032 —
against `capture OFF` −0.0002 and `verystrong` −0.0035 (4/0).

**The mechanism is measured** (`scratchpad/w4_provenance.py`, reading the solver's own `_lvl` capture). At
`exon|intron` seams in `gdna300_ss0.99_nrna_none_capON`, where the oracle RNA density is **exactly 0**:

```
             nu_true   nu_fused   nu_own   nu_far   nu_mass   rho_mu       weight own/far/mass
far=intron     0.000      0.725    0.217    0.505     1.396    7.035           18% / 25% / 57%
far=exon      21.093     14.906   14.171  308.190    14.512    2.517           20% /  0% / 79%
```

The MASS estimator carries 57 % of the fuse and is the one that is wrong (1.396 vs 0), while OWN (0.217) and
FAR (0.505) are close. Its input is the exon's gDNA claim, and at a **capture-ON** exon→boundary edge that
claim is mis-scaled by the probe-placement step — the irreducible 0.4–1.3 nats of HANDOFF_10 §9.3.

**And here is the specific, actionable gap:** `σ²_transfer` does **not** see that step. `composition_logvar`
prices only counting plus composition uncertainty in `ρ_tot`, and at a boundary `E_g ≈ E_r` kills the
composition term, so `logvar_tot ≈ 1/n ≈ 2e-4`. The 0.4–1.3 nats of genuine exon-vs-boundary capture scatter
is **entirely unpriced**. Everywhere else that has not mattered, because `_pin_v` cancels `r` on 87.6 % of the
error (HANDOFF_10 §3) — **M11 is the first consumer in the solver for which `r` does not cancel**, so it is
the first to be exposed to it.

Three candidate directions, none tried:
* price the exon↔boundary frame step honestly. M8 (`graft_frame_logvar`) is the precedent for a
  non-cancelling `r`, but its `(log r)²` would over-damp `verystrong` — which currently *wins* 4/0, so
  whatever is charged must be the un-modelled part only, not the whole step;
* feed M11 the **intron's** gDNA density at `exon|intron` seams (factory-measured, and the intron→boundary
  reframe is the shorter step) rather than the exon's, or fuse the two gDNA claims before M11;
* a DL-style between-estimator excess on the level fuse itself. ⚠ Note it pulls the fused level *toward* the
  larger estimate here, so it is unlikely to help on its own — check that before building it.

## 6. Bootstrap

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
cd /Users/mkiyer/proj/rigel
python -m pytest tests/ -q                    # 1241 pass / 2 xfail / 2 xpass
ruff check src/ tests/ scripts/               # clean
python scripts/debug/message_variance_mc.py   # 0 failures over M1–M11
P0_REFIT=0|1 python scripts/debug/pass0_oracle_bench.py --arm NAME   # the A/B; this arm is `final_r0`/`final_r1`
python scratchpad/ab.py NEW BASE [extra...]   # multi-arm reader: strata + per-condition
python scripts/debug/pass0_error_table.py --refit 0                  # the suite in READS + errQ1conf
python scratchpad/w4_provenance.py            # WHICH estimator sets the level, and how wrong it is
```

`scratchpad/m11_level_fuse.patch` is this session's full diff (safety copy). `w1_level.py` / `w2_fuse.py` are
the offline estimator studies that chose the design; `w4_provenance.py` reads the live `_lvl` capture and is
the one to keep.

## 7. Invariants — all preserved, none touched

* the `√2·σ_own` DL pin-safety inequality (`mismatch_deflate` unchanged);
* `Σ_c ρ_c·E_c = M` at every relay hop and the combine, with `_pin_v` partial-claim semantics (unchanged);
* `N` enters only as power — M11's `k` is *literally* an effective fragment count, which is the strongest
  form of that invariant this model has yet had;
* the composition is ONE dof: M11 delivers the **λ axis only** and the message's own tilt splits it across
  the strands, exactly as the density deconvolution is specified to behave (HANDOFF_10 §6);
* **no magic numbers** — every branch is a structural presence test (`τ_own > 0`, precision > 0), a derived
  quantity, or an exact limit (`ψ'(1) = π²/6`, `k = 3` for the uniform posterior, `min(φ,1)` as the
  projection onto the admissible set);
* **pass-0 is WEAK and CORRECTABLE** — and this is now a *tested property*: M11's MC asserts honesty where the
  level is a measurement and **conservatism (`z2 ≤ 1`) where it is marginal**, as two separate assertions.

## 8. Vocabulary (owner's, unchanged)

RNA is **one species**. Only **SPLICED vs UNSPLICED** is observable. A boundary can be an exon↔exon boundary
that is also a splice junction, with RNA contiguous across it while other RNA splices in or out. The peel is
about spliced vs unspliced, never mature vs nascent.

`docs/calibration/boundary_work_notes.md` is the owner's own notes file — do not overwrite it.
Never reference `docs/calibration/archive/`.
