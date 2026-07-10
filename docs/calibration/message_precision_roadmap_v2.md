# Message-precision roadmap v2 — resume with a clean signal

**Scope:** `src/rigel/calibration/bp_solver.py` (the BP-sweep message precision).
**Supersedes:** `message_precision_regression_and_fix_plan.md` (the paused bookmark).
**WIP code:** `calib-msgprec-eb-maxform-wip` (EB shrinkage + env-gated max-form, both implemented).
**Status:** ⏸️ **PAUSED 2026-07-10 — valuable, not abandoned.** Resumed 2026-07-09 (confounds fixed);
Phases 0 & 1 ran (results below). The easy suite has little message-precision headroom (prod is already the
best message rule where messages help), so the message-precision thread is paused with a clear open
candidate (the `f_g`-resolved EB, to be tested on the AMBIG suite — its real arena). Work pivoted to the
**capON local-solve residual** (the message-independent exonic-gDNA under-call) that Phases 0 & 1 isolate
as the dominant error. **What survives as valuable:** (a) the `RIGEL_MSG_MODE` harness + the soft-3-pool
bucketed measurement; (b) `off` as a diagnostic floor; (c) the `f_g`-resolved EB candidate (Change 3, plumbing
ported, untested); (d) the confirmed map of where messages help/hurt (below). What is dead is only the
*specific forms* the tests refuted (max, total-density EB) — not the idea of disagreement-aware precision.

## Why now (the precondition is finally met)

The paused doc was explicit: *"Message precision cannot fix a gDNA-independent error. Fix the siphon
first … then resume here with a clean signal."* Both masking errors are now gone:

- **Nascent siphon** — the eff-length inversion + per-transcript reference density that manufactured a
  gDNA-independent mature→nascent leak. Fixed and shipped: `8add036a` (global KDE-mode reference density),
  `48133062` (unified per-node eff-length contraction).
- **Fragment length** — interrogated and **ruled out** this session. Two controlled tests (no-RNA gDNA-only;
  flipped DNA=100/RNA=200) show capture does **not** size-select gDNA fragment length (captured FL == input
  FL); the earlier "on-target shorter" was a midpoint-classification artifact. The machinery is unbiased
  under equal FL. FL is not a siphon cause.

So for the first time the value of message propagation can be measured against a **clean** target.

## The measurement contract (non-negotiable)

- **Metric = soft 3-pool net surplus** (assigned − true, fragments, for gDNA / nascent / mature). The
  hard-label `net_flow` is *insensitive* to a soft-prior shift (a real change moves the 3-pool counts by tens
  of thousands of fragments while `net_flow` is byte-identical) — never use it to gate this work.
- **Siphon = SYNTHETIC nascent EM mass** (`estimator.nrna_em_count` = `t_counts[_synthetic_mask]`).
  Single-exon `is_nrna` transcripts are **non-synthetic** (both mature+nascent, molecularly
  indistinguishable) and count in the **mature** pool — never in the siphon. This is audited in the harness:
  `quantification.nrna_total` is exactly `nrna_em_count`, and `mrna_total` (= `get_counts_df` over
  `~is_synthetic`) keeps single-exon `is_nrna`. The three measurement traps: `scripts/debug/_metrics.py`.
- **Bucket every result** by capture × strand × gDNA level × nascent level. A single scalar verdict hid the
  max-form's stranded-win / unstranded-loss last time.
- **Baseline = messages OFF** (`RIGEL_MSG_MODE=off` forces `pr=0` in all three channels — the code's own
  designed no-message state). Judge every candidate against OFF, not against production.
- **Suite progression — easy → hard.** Start on `quick_3to1_5mb` (16 conditions: gDNA {none, 300} × ss
  {0.99, 0.50} × nascent {none, rnd} × capture {off, on}) — mostly stranded / single-strand loci, where
  strand carries most of the solve. Then build up to the **AMBIG** suite (`ambig_dense_10mb`, 305 AMBIG
  regions) where strand cannot solve the node and **message propagation is required** — that is the true
  test of a propagation scheme, not the easy quick suite.

## Design principles (the frame for this work)

1. **Adjacent-node pair disagreement informs message precision** — the cornerstone; keep it.
2. **Precision must be controlled** — two nodes can agree by chance; model that likelihood and temper
   runaway high precision. A single agreeing residual is *not* evidence of high precision.
3. **Precision can safely go to zero** — no harm; a muted message is just an ignored one.
4. **Measure benefit vs harm against a NO-MESSAGE baseline** — messages OFF is the floor everything is
   judged against.
5. **Cover the full matrix** — capture on/off, unstranded/stranded, gDNA levels, nascent levels.
6. **Locus complexity governs propagation** — messages are not always needed; a node that can solve itself
   should not be overridden by a neighbour.
7. **Baseline + status quo, then roadmap.**

## Disposition of prior work

| Idea | Verdict | Disposition |
|---|---|---|
| Adjacent-node disagreement → precision (#1) | Cornerstone, sound | **Foundation — keep** |
| Non-circular anchor (resid = message − dst *message-free* belief) | Settled correct; avoids the echo chamber | **Keep** (implemented on WIP) |
| Source-**belief** base_var (`vbg[lsrc]`) | Collapses → precision runaway ~1e9 | **Dead — never again** |
| Global count-only scalar `pr = n_src/(n_src·σ²_imp+1)` (production) | No self-silencing → capture-ON seam bleed | The incumbent to beat |
| Symmetric EB shrinkage, 1 residual/edge | `w=0` on flagship — can't distinguish edges, collapses to the scalar | **Dead *as-is* → needs residual pooling** |
| Asymmetric max-form self-silencing | Stranded win / unstranded loss (net +16% gDNA / −16% nascent) — siphon-confounded | **Best lever — re-baseline clean** |
| `net_flow` as the ship metric | Insensitive to soft-prior shifts | **Dead — soft 3-pool only** |

## The incumbent

Production message precision (bp_solver.py, all three channels — gDNA / RNA⁺ / RNA⁻):

```python
pr = n_src / (n_src * sig_imp + 1.0)   # σ²_msg = σ²_imp + 1/n_src ; σ²_imp fit once, one global scalar
```

No per-edge residual term, no self-silencing. Cures the v0.6.4 belief-derived runaway but drops the seam
term — the capture-ON regression the paused doc characterised.

## Roadmap

### Phase 0 — the no-message baseline (this doc's first experiment)

Add the `RIGEL_MSG_MODE=off` switch (done: bp_solver `pr=0`). Run the 16-condition soft-3-pool matrix,
bucketed, for `{prod, off}`. **Question:** does the current global-scalar message layer net-help or net-hurt
vs nothing, per bucket, now that the siphon is gone? *Expectation (to be confirmed): stranded ≈ OFF (strand
carries it); unstranded worse under OFF (relay is load-bearing).* Driver:
`scripts/debug/benchmark_ab_report.py --arms prod:0 off:0::RIGEL_MSG_MODE=off`.

#### Phase 0 RESULTS (2026-07-09, quick_3to1_5mb, soft 3-pool, |off|−|prod| = message value)

Message value on the **gDNA/RNA split** (mean |gdna surplus|, POSITIVE ⇒ messages help):

| bucket | \|gdna\| prod | \|gdna\| off | message value |
|---|--:|--:|--:|
| ss99 capOFF (stranded, no capture) | 6,624 | 10,674 | **+4,050 (help)** |
| ss99 capON  (stranded, capture) | 57,363 | 53,911 | −3,452 (~neutral) |
| ss50 capOFF (unstranded, no capture) | 18,061 | 62,798 | **+44,736 (strong help)** |
| ss50 capON  (unstranded, capture) | 164,374 | 56,941 | **−107,433 (catastrophic HURT)** |

**Findings.** (1) STRANDED ≈ neutral — strand carries the solve (confirms the prediction). (2) UNSTRANDED
capOFF: messages strongly HELP (+44k) and without them RNA→gDNA false positives explode (gNONE nrRND
capOFF: off invents 161k false gDNA) — the relay is load-bearing (confirms the prediction). (3) UNSTRANDED
capON is the OPPOSITE: production messages leak gDNA→RNA by up to **−347k** fragments (the capture
**seam bleed** the paused doc predicted) — turning messages OFF fixes it. Strand normally vetoes this
message; unstranded cannot. **The harm is seam-localized (capON); the help is seam-free (capOFF).**

**Second, decoupled finding — the capON nascent siphon is NOT primarily message precision.** On gdna300
capON the pure siphon (nrNONE) is ~90–120k under BOTH prod and off and BOTH strands (ss99: 104k/113k;
ss50: 121k/92k) — a message-independent, strand-independent residual (exonic-gDNA under-call). Message
precision is necessary-not-sufficient for the flagship: it can recover the −330k unstranded-capON split
leak, but a separate ~90k capON calibration residual remains (parallel target — reference-density /
exonic-gDNA under-call, not this doc).

**Revised Phase-2 priority (data-driven):** the **asymmetric max-form** (self-silence across density seams)
is now the PRIMARY lever — it directly targets the −107k capON seam harm while preserving the +44k capOFF
relay (no seams there to silence). Strand-information gating is the secondary safety valve, most valuable on
the AMBIG suite where the relay must stay alive. Phase-1 acceptance test: recover the unstranded-capON leak
WITHOUT regressing the unstranded-capOFF relay.

### Phase 1 — re-baseline the two implemented candidates (cheap; code exists)

Cherry-pick the WIP max-form + EB onto main behind `RIGEL_MSG_MODE` (`max`, `eb`); run them in the same
matrix vs `{off, prod}`. **Question:** does the max-form's stranded-win / unstranded-loss survive the siphon
fix, or was part of it the siphon? Fold the WIP env gate (`RIGEL_MSG_MAXFORM`) into `RIGEL_MSG_MODE`.

#### Phase 1 RESULTS (2026-07-09, quick_3to1_5mb, mean |gdna split err| per bucket)

| bucket | prod | off | max | eb |
|---|--:|--:|--:|--:|
| ss99 capOFF (stranded, no cap) | **6,299** | 10,473 | 8,002 | 6,394 |
| ss99 capON  (stranded, cap) | 57,599 | **53,950** | 61,544 | 57,827 |
| ss50 capOFF (unstranded, no cap) | **18,154** | 62,734 | 38,765 | 18,116 |
| ss50 capON  (unstranded, cap) | 164,468 | **57,196** | 196,849 | 164,568 |

**`max` is WORSE than both prod and off in every bucket → rejected.** Self-silencing anchored on the dst's
*local* belief entrenches the corrupted capON local belief (the exonic-gDNA under-call): it mutes the
*corrective* gDNA-neighbour messages (high disagreement with the wrong local) while keeping the
*leak-causing* RNA messages (low disagreement). The old "stranded win" was a **siphon-confound artifact** —
it vanished once the siphon was fixed. **`eb` ≡ prod** (`w=0` measured on all conditions — the total-density
frame has no between-edge signal). Inert.

**The deeper lesson.** No single message rule threads the suite: **prod wins where the relay is load-bearing
(capOFF-unstranded, 18k vs off's 63k), off wins where the seam bleeds (capON)**. The capON failure is NOT
fixable by a message rule — both `max` and `eb` anchor on the local belief, which at capON exons is *itself*
wrong (message-independent). So on the easy suite there is **little message-precision headroom: prod is
already the best message rule where messages are the right tool**, and the capON leak is a LOCAL-SOLVE
(reference-density / exonic-gDNA under-call) problem — the separate parallel target.

**Revised next step.** The one untested, principled candidate is the **`f_g`-resolved EB (Change 3)**:
measure the shrinkage weight from the *deconvolved* gDNA count `f_g·mass` after pass 1 (where capture's
coherent gDNA enrichment creates real between-edge signal ⇒ `w>0`), not per-edge disagreement with a
corrupted local belief. This is "self-silence on high disagreement" done right — disagreement of the
reliable gDNA *field*, population-shrinkage-bounded. **Test it on the AMBIG suite** (`ambig_dense_10mb`),
where the relay is load-bearing and a better rule can actually earn its keep — the easy suite cannot
distinguish it.

### Phase 2 — the improvement: **asymmetric self-silencing + strand-information gating**

Headline recommendation, motivated directly by the data:

- **Asymmetric max-form answers principle #2 by construction.** Precision is *capped* at the population
  average `σ²_imp + 1/n`; a lucky single agreement can never manufacture high precision — only genuine
  disagreement (a reliable signal at a seam) is allowed to *lower* it. Parameter-free runaway control.
- **Strand-information gating answers principles #5 and #6.** The stranded-win / unstranded-loss result *is*
  the design: self-silence the relay where the node can already solve itself (high local strand Fisher
  information — stranded, count-rich, simple loci), keep the relay load-bearing where strand cannot
  (low strand info — unstranded, AMBIG, complex loci). Scale the self-silencing strength by the node's local
  strand information. Each node asks *"do I need help?"* before listening — principle #6 made mechanical.
- **Pooled-residual EB is the fallback refinement.** Where a stratum shows real between-edge signal
  (`τ²_between > ψ'(½)`), aggregate residuals so `w` rises above 0 and below-average-disagreement edges earn
  *above*-average precision. The EB idea was starved of residuals, not wrong.

The AMBIG suite is the acceptance test for Phase 2: the easy quick suite can look fine with weak messages;
only the AMBIG loci prove a propagation scheme actually works.

## Dead — do not revisit

Belief-derived precision (`vbg`) → runaway. `net_flow` as the metric. Symmetric EB with one residual per
edge (`w=0`). "On-target gDNA is shorter" (a classification artifact; FL does not size-select).

## Code pointers

- Incumbent + switch: `src/rigel/calibration/bp_solver.py` (`node_sweep._scan`, `RIGEL_MSG_MODE`,
  `adjacent_disagreement_variance`).
- Candidates (WIP): `calib-msgprec-eb-maxform-wip` — `_eb_edge_precision`, `_max_edge_precision`,
  `_adjacent_log_density_residuals`, `adjacent_disagreement_shrinkage_weight`.
- Harness: `scripts/debug/benchmark_ab_report.py` (soft 3-pool, env-arm axis) + the `calibration-benchmark`
  skill. Suites: `quick_3to1_5mb` (easy), `ambig_dense_10mb` (hard).
- Metric guardrail: `scripts/debug/_metrics.py`.
