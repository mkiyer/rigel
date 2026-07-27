> ## ⛔ SUPERSEDED — THIS IS NOT THE LIVE HANDOFF
> The live handoff is **`SESSION_2026_07_26_HANDOFF_14.md`**; the entry point is **`ROADMAP.md`**.
> This file is kept as HISTORY — what was tried, measured and refuted. Its numbers describe the
> code as it was on its own date and many are now superseded. **Do not act on it without checking
> the live handoff first.** Its own DO-NOT-RE-RUN findings (if it has any) remain HERE and still
> stand — HANDOFF_14 §3 indexes which files carry them.

# Session handoff — the single-strand × capture study is DONE; M8 landed; AMBIG is next

*(this file's own original header — superseded)* It supersedes
`SESSION_2026_07_25_HANDOFF_7.md` (whose §4–§5 "next study" is now COMPLETE — this file is its result).
Do NOT read `docs/calibration/archive/`. Date: 2026-07-25. Branch `calib-ambig-init-wip`.

---

## 0. TL;DR

HANDOFF_7 §4 asked one question: under capture, does single-strand pass-0 break at (a) the node's own
self-solve, (b) the message MODES, or (c) the message PRECISIONS? **Answered: (b), and the mechanism is now
identified exactly.** The degradation is 77–92 % on EXONS and it is entirely a message defect. The delivered
exon composition is `ρ_g^src·E_g : (ρ_ν^src + ρ_μ)·E_r` — the reframe `r` and the pin factor `k` cancel to
`1.8e-15`. `ρ_μ` (the grafted junction flux) **is** the RNA claim (`ρ_μ/ρ_ν` = 400–37,000), and it is a
spliced measurement whose fragments are anchored in the flanking EXONS, so it already sits in the
destination's frame — while `ρ_g^src` sits in the boundary's, 6× lower under capture. Since `r` cancels, the
graft edge **never reframes the gDNA at all**. That is the whole 10× degradation.

The fix that shipped (**M8**) prices the un-cancelled frame step as a variance, `(log r)²`, on the grafted
component. Derived, MC-validated (rel ≤ 0.24 %), and A/B-won: **0.0926 → 0.0900 (refit=0)**, **0.0779 →
0.0700 (refit=1)**. **Next study: AMBIG**, the other ~50 % of the error mass (HANDOFF_7 §3).

## 1. What landed

| file | what |
|---|---|
| `enrichment_frame.graft_frame_logvar` | the M8 pure law + 3 closed-form unit tests |
| `bp_solver` (`_relay`, `_transport`) | M8 wired on the graft's spliced measurement, both places |
| `message_variance_mc.py` | the M8 MC arbiter case (`r = 1 / 2.6 / 6.1`) |
| `message_variance_derivation.md` §8 | the derivation, the measurements, and the mode-vs-variance decision |

**Gates:** `pytest tests/ -q` = **1232 pass, 2 xfail, 2 xpass** (1229 + 3 new); `ruff check src/ tests/
scripts/` clean; `message_variance_mc.py` **0 failures** over M1–M8. Goldens regenerated LAST.

## 2. The measurement chain (each step is a script in `scratchpad/`, all re-runnable)

| step | script | result |
|---|---|---|
| 1 | `ss_cap_1_decompose.py` | exons carry 77–92 % of capture-ON single-strand error; Δmsg +0.0005 → +0.052. In the `nrna_none` pair the exon SELF-solve is unchanged (0.0049 → 0.0049) while SOLVED degrades 4× — a pure message defect on a clean control |
| 2 | `ss_cap_2_exons.py` | on the worst exons `self ≈ oracle` and solved is pulled DOWN; the RNA measurement channel out-precisions the (correct) gDNA one |
| 3 | `ss_cap_3_replay.py` | **bit-exact ψ replay** (`max|Δ| = 0`) + channel ablation. No single channel owns it; `-gdna` makes it WORSE (the gDNA message is right) |
| 4 | `ss_cap_4_calib.py` | the reliability test. Exon RNA channel: **z2 = 1.0 / 1.1 / 1.0 off-capture** (perfectly calibrated) → **5.0 / 4.9 / 7.0 on**, and the error is almost pure BIAS (+1.42 / +1.05 / +0.53 nats) |
| 5 | `ss_cap_5_bias_source.py` | `r` and `k` cancel from the delivered share (fid `1.8e-15`); `ρ_μ/ρ_ν` = 400–37,000 ⇒ the graft IS the claim |
| 6 | `ss_cap_6_graft_frame.py` | the required per-edge frame factor: **median `log c` = +0.009 / −0.008 / +0.054 off-capture** (shipped model EXACT) and −1.26…−1.98 on; **`log r` does not predict it** (corr ≤ 0.35) |
| 7 | `ss_cap_7_oracle_channels.py` | **the root cause**, oracle-only: `ρ_R(exon)/ρ_spl(bnd)` = 1.02–1.86 and capture-INVARIANT, while the gDNA step exon↔boundary goes 1.03 → **6.1–6.8** |
| 8 | `ss_cap_8_term_check.py` | the term is right-sized: graft-edge `z2` **58–310 → 2.1–3.8**, uniformly across capture off and on |

⚠ **A diagnostic-side trap, now documented:** the oracle pools (`_oracle_per_node`) use the OPPOSITE ±
convention from `statics.free_pos`/`free_neg`, so a `free_pos` node carries its oracle RNA in `Rn`. Scoring a
per-strand channel against `Rp`/`Rn` produces a spurious constant ≈ +4.6-nat bias. **Score on the composition
axis (gDNA vs RNA-total); the tilt is a separate DOF anyway.**

## 3. The M8 law (see `message_variance_derivation.md` §8 for the full derivation)

```
    Var(log ρ_μ^msg)  +=  (log r)²          on a GRAFT edge only, 0 at r = 1
```

M5's graft exemption (`σ²_transfer = 0`, "`r` is common-mode across `{g, R}`") is correct for the imputed
`ρ_ν` and **false for the measured `ρ_μ`**, which has no matched gDNA partner to cancel `r` against — exactly
M5's peel/partial case. No tuned constant: it is the method-of-moments second moment of a single observation
of the un-cancelled step, the same logic as M7's `b̂²`. Applied to the spliced measurement's own precision, so
M2's `w_μ²` share weighting arises implicitly from the inverse-variance fusion.

## 4. ⭐ THE ABLATION — why VARIANCE and not the mode fix (this is the load-bearing decision)

Four arms, all A/B'd at both refits. `pin` = HEAD.

| arm | what | r0 agg | r0 b/w/f | r1 agg | r1 b/w/f | `gdna_none` |
|---|---|---|---|---|---|---|
| `pin` | HEAD | 0.0926 | — | 0.0779 | — | — |
| **`m8`** | **M8 variance only — SHIPPED** | **0.0900** | 17/15/0 | **0.0700** | **22/7/3** | **improves** |
| `m9` | M8 + graft in the dst frame | 0.0850 | 14/18/0 | 0.0591 | 20/9/3 | regresses |
| `m10` | m9 + drop M5's graft exemption | 0.0845 | 12/18/2 | 0.0652 | 12/19/1 | wrecks |
| `m11` | mode fix ONLY (= §5.1 today) | 0.0853 | 9/23/0 | 0.0666 | 10/20/2 | wrecks |

(counts at a 1e-4 threshold; `pass0_oracle_bench.py --report` quotes 15/9/8 and 20/6/6 for `m8` on its own
looser threshold. Rows are in `/tmp/pass0_oracle_bench.tsv`.)

**Three findings:**
1. **`m11` reproduces §5.1's documented revert exactly** on today's post-relay-pin baseline: `gdna_none`
   capOFF 0.3668 → **0.4369**, capON 0.3780 → **0.4917**. §5.1 stays REVERTED, and now for a re-measured
   reason, not a remembered one.
2. **`m10` confirms M5's graft exemption is real** — removing it loses everywhere (0.0652 r1, 12/19).
3. **`m8` is the only arm that improves `gdna_none`** (7 better / 0 worse / 1 flat at refit=1). It can only
   *remove* confidence; it never moves a mode toward a wrong answer. `m9`/`m11` buy their larger aggregate
   with exactly the regression `density_composition_reconciliation.md` §3.3 identified as forbidden.

**The governing principle decided this against the aggregate.** Worth remembering as precedent.

## 5. Known cost of M8, and the open refinement

`r` conflates the capture step with genuine density structure: off-capture `r ≈ 1.15–1.3` boundary→exon
(the exon really does carry more RNA per base) while the true frame step is 1.03. So M8 over-damps
off-capture and under-damps on-capture (`(log 2.6)² = 0.91` vs a true `(log 6.4)² = 3.4`). The cost lands on
**capture-OFF unstranded**, where the graft is the only RNA information:

| condition | pin | m8 |
|---|---|---|
| `gdna100_ss_0.50_nrna_none_capture_off` | 0.0866 | 0.1129 |
| `gdna300_ss_0.50_nrna_none_capture_off` | 0.0504 | 0.0631 |
| `gdna100_ss_0.50_nrna_present_capture_off` | 0.0737 | 0.0837 |

**The refinement: isolate the ENRICHMENT part of `r` prior-free.** gDNA is genomically uniform in content, so
a ratio of true gDNA densities IS a capture ratio — but the solver's own `ρ_g` at an unsolved node is the
100 %-gDNA default, so the naive `og/rg[src]` is not usable as-is. This is the next thing to derive if you
return to single-strand.

## 6. ▶ NEXT — AMBIG (HANDOFF_7 §3)

Single-strand is now materially better and the remaining half of the error mass is AMBIG, concentrated on
**`exon+intron(x-strand)`** (25.8 % of capture-ON error) — a region that is one strand's exon and the other's
intron. Prior-free, AMBIG nodes have `τ_own = 0` (no composition evidence at all: the strand likelihood
constrains only the tilt, never `f_g`), so M7's DL term is inert there by design and they are carried by
messages alone. Re-run the census (`scripts/debug/derive_4_boundary_classes.py`) on the new HEAD before
starting — M8 will have moved it.

⚠ **Do not read this as "AMBIG needs the hyperprior, so go do Phase 2."** The owner's directive stands: the
pass-0 solver must be correct first. The M8 result is the precedent — three sessions assumed the remaining
error needed the prior, and it turned out to be a message-frame defect measurable with oracle densities and
no prior at all.

## 7. ⛔ DO NOT RE-RUN — this session's additions to HANDOFF_7 §7

| item | verdict |
|---|---|
| **Graft outside the reframe (§5.1) on the post-pin baseline** | Re-measured as `m11`: 0.0853 / 0.0666, **9 better / 23 worse** at r0, `gdna_none` 0.3668 → 0.4369. Stays reverted. |
| **Graft in the dst frame + M8 (`m9`)** | Better aggregate (0.0850 / 0.0591) but carries §5.1's `gdna_none` regression. Rejected on the governing principle. |
| **Dropping M5's graft exemption from `σ²_transfer` (`m10`)** | Worse everywhere (0.0652 r1, 12/19). M5's exemption is confirmed correct for `ρ_ν`. |
| **`1/r` as the graft's frame factor** | REFUTED by direct measurement: required `log c` = 0 off-capture (so `1/r` over-corrects, `|log c|` 0.49 → `|log c + log r|` 0.77) and `log r` does not predict `log c` on-capture (corr ≤ 0.35). |
| **Per-strand oracle scoring of message channels** | Invalid — the `±` convention flip (§2). Score the composition axis. |

Everything in HANDOFF_7 §7 remains in force.

## 8. Tools and environment

Unchanged from HANDOFF_7 §6. The new per-node ψ replay (`scratchpad/ss_cap_3_replay.py`) is **bit-exact**
(`max|Δ| = 0`) for single-strand nodes and supersedes `pass0_node_dissect.py` for channel ablation (whose
replay fidelity is stale — HANDOFF_7 §8.5). It needs only `cap["_uni"][-1]` + `fg_init`/`intron_prior`
/`global_lp`; `theta_imp` is recomputed from `cp`/`cn` and is inert on single-strand nodes.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
python scripts/debug/suite_dissect.py --out /tmp/suite_nodes.npz    # ~4 min, 74,494 node rows
P0_REFIT=0 python scripts/debug/pass0_oracle_bench.py --arm NAME    # HEAD's arms are now m8_r0 / m8_r1
```
