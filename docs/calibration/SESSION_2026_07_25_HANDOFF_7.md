# Session handoff — pass-0 CORRECTNESS first; single-strand × capture is the next study

**This is the LIVE handoff. Read `docs/calibration/ROADMAP.md` first, then this.** It supersedes
`SESSION_2026_07_25_HANDOFF_6.md` (whose §3 "next task is Phase 2" is **withdrawn** — see §1).
Do NOT read `docs/calibration/archive/`. Date: 2026-07-25. Branch `calib-ambig-init-wip`, HEAD `d88a7822`.

---

## 0. TL;DR — what to do, in one paragraph

Pass-0's message model is finished and winning (aggregate 0.1267 → **0.0926** refit=0, 0.1234 → **0.0779**
refit=1). A suite-wide census then relocated everything that is left, and it is **not** where the last three
sessions assumed. **Off-capture, single-strand pass-0 is essentially SOLVED — mwae 0.0051** (exons 0.0046,
retained introns 0.0051, intergenic exactly 0.0000). Two axes remain: **capture degrades single-strand 10×**
(0.0051 → 0.0476), and **AMBIG carries ~50 % of all error mass** at 3.7–16× the single-strand rate.
**The owner's directive: the pass-0 solver must be CORRECT before the gDNA hyperprior fit, and single-strand
must be debugged before AMBIG.** So: **the next study is single-strand × capture** (§4) — a clean controlled
comparison, because the same nodes with the same annotations score 0.0046 without capture and 0.0499 with it.
§5 is the exact experiment. §7 lists what is already REFUTED so you do not re-run it.

## 1. The direction, and one withdrawn recommendation

Calibration deconvolves each node's **unspliced** fragment mass into `(f_rna₊, f_rna₋, f_g)`. Pass-0 is the
**prior-free first pass**: it solves from only the two *intrinsic* sources (strand likelihood + cross-node
imputation), then its output trains a **gDNA hyperprior**, which is required to re-solve.

**The order is now fixed by the owner: pass-0 must be correct BEFORE the hyperprior fit.** An earlier
recommendation to move to Phase 2 (because three items in a row resolved as "the physics is real but the
pass-0 fix is not") was **withdrawn**. The census in §3 is why: it shows a large, structured, *correctable*
residual that has nothing to do with the missing prior.

**Governing principle (do not violate).** Pass-0 must be **WEAK and CORRECTABLE**. An over-confident message
that PINS a node wrong is far worse than a weak one that is slightly off. A node that lacks the data to solve
itself must end up at near-zero precision or unsolved — **never** moderate precision. Prefer the
under-confident option when unsure. **No magic numbers**: structural presence tests, derived quantities and
exact limits only. Counts are Poisson (the synthetic suite is Poisson by construction, so nothing
overdispersion-dependent can be validated on it).

**The owner's model vocabulary (use it).** RNA is just RNA — "mature" and "nascent" are **not** different
species, they are RNA in two places. The only distinction the data supports is **SPLICED vs UNSPLICED**. A
boundary can be an exon↔exon boundary that is *also* a splice junction: RNA can be contiguous across it while
other RNA splices in or out. **Both happen at once.**

## 2. State: what landed this session (all on `calib-ambig-init-wip`)

| commit | what |
|---|---|
| `f7c02c8e` | **DerSimonian–Laird composition-mismatch cliff term `b̂²`** — replaces the `(log r)²` proxy |
| `1a3e0a89` | retire the dead NPMLE σ²_transfer plumbing; struct_lock resolution + interior-anchor test |
| `64bb8ec3` | docs + goldens: the message-variance model is COMPLETE |
| `b81e926e` | **the λ-emission gate** — a source with only one component has no composition claim to make |
| `aa920e59` | the density↔composition reconciliation brief |
| `106616cd` | **anchor the relay context to the node's observed mass** (the relay pin) |
| `87ddb9f8` | §4.2 per-channel enrichment ratios — physics CONFIRMED, fix **REFUTED** |
| `364701d9` | §3.3 the RNA floor — **NOT achievable prior-free** |
| `d88a7822` | the boundary-class census (§3 below) |

**A/B, `scripts/debug/pass0_oracle_bench.py`, `OMP_NUM_THREADS=1`:**

| | session start | **HEAD** | | | session start | **HEAD** |
|---|---|---|---|---|---|---|
| refit=0 ALL 32 | 0.1267 | **0.0926** | | refit=1 ALL 32 | 0.1234 | **0.0779** |
| stranded ss_0.99 | 0.0262 | 0.0377 | | stranded ss_0.99 | 0.0293 | 0.0358 |
| unstranded ss_0.50 | 0.1981 | **0.1315** | | unstranded ss_0.50 | 0.1902 | **0.1077** |
| unstr × capON | 0.2543 | **0.1720** | | unstr × capON | 0.3874 | **0.1585** |
| verystrong | 0.2985 | **0.1902** | | verystrong | 0.1514 | 0.1350 |
| capture OFF | 0.0604 | **0.0439** | | capture OFF | 0.0059 | 0.0311 |

Gates: `pytest tests/` = **1229 pass, 2 xfail, 2 xpass**; `ruff check src/ tests/ scripts/` clean;
`scripts/debug/message_variance_mc.py` = **0 failures** over M1–M7. Goldens regenerated. Working tree clean
except untracked `scratchpad/` (keep it — it holds the derivation scripts).

**The message-precision law, as it now stands** (`message_variance_derivation.md`, M1–M7):

```
    p = 1 / ( Var(log f_c^src) + 1/n_src  +  σ²_transfer  +  b̂² )
             \__ strand ___/   \_count_/    \_ SCALE _/    \_ COMPOSITION _/
```

No prior of any kind enters message precision. Its safety property is exact and is an INVARIANT:
`p_eff = 1/max(v_msg, G²−v_own)` ⇒ **a message can out-weigh a node's own belief only if it agrees to within
`√2·σ_own`.** Any change that breaks this is a regression even if the A/B improves.

## 3. ⭐ THE CENSUS — where the error actually is (this is the map; trust it)

Regions classified by their **RAW 4-bit signature bits**, not the coarse type
(`scripts/debug/derive_4_boundary_classes.py`). Stranded conditions, refit=0:

| region class | DOF | capture OFF | capture ON | share of err (capON) |
|---|---|---|---|---|
| exon | single | **0.0046** | 0.0499 | 43.8 % |
| RETAINED intron | single | **0.0051** | 0.0746 | 1.8 % |
| intron | single | 0.0136 | 0.0378 | 1.0 % |
| intergenic | single | **0.0000** | **0.0000** | 0.0 % |
| exon | AMBIG | 0.1759 | 0.1879 | 17.3 % |
| **exon+intron(x-strand)** | AMBIG | 0.1191 | **0.1684** | **25.8 %** |
| RETAINED intron | AMBIG | 0.0806 | 0.1881 | 0.9 % |
| **TOTAL single-strand** | | **0.0051** | **0.0476** | 50.7 % |
| **TOTAL AMBIG** | | **0.0845** | **0.1771** | 49.3 % |

Boundaries by flanking pair: `exon | exon` **0.0305** (no junction) / **0.0477** (with junction) — both BETTER
than the suite average. `exon | intron` with junction 0.0530. The worst all touch the cross-strand region:
`exon | exon+intron(x-strand)` 0.2468 / 0.2708, `RETAINED | exon+intron(x-strand)` 0.2691.

**Read this as three facts:**
1. **Structure is handled.** Retained introns, exon|exon boundaries (with and without a junction), and
   intergenic anchors are all correct on single strand. The owner's `TA/TB` retained-intron case is
   represented exactly by the signature and scores identically to a plain exon (0.0051 vs 0.0046).
2. **Capture is a 10× degradation on single-strand** (0.0051 → 0.0476) with everything else held constant.
   **This is the next study.**
3. **AMBIG is half the error**, concentrated on `exon+intron(x-strand)` — a region that is one strand's EXON
   and the other strand's INTRON. That is the study after.

Suite-wide (all 32 conditions, `scripts/debug/suite_dissect.py` → `/tmp/suite_nodes.npz`): 92.2 % of error
mass sits on `τ_own = 0` nodes (no composition evidence of their own) and messages HALVE their error
(0.333 → 0.164); full-rank nodes self-solve at 0.011 and are degraded to 0.025 by messages.

## 4. ▶ THE NEXT STUDY — single-strand × capture

**Why this one.** It is the only place with a **clean control**: the identical nodes, annotations and
structures score **0.0046 off-capture and 0.0499 on-capture**. Nothing varies but the enrichment. Any
mechanism you find is therefore attributable, and no strand ambiguity confounds it.

**The question to answer, in this order:**
1. Does capture break the node's **own self-solve**, the message **modes**, or the message **precisions**?
   The suite table has `self` (message-free self-solve), `solved`, and every channel's mode + precision, so
   this is a decomposition, not a guess.
2. Whichever it is, localize it to a **class** (region class × boundary class from §3) and then to
   **individual nodes**, and trace the provenance the way node 2197 was traced.
3. Only then propose a fix — derive it, MC-validate it if it is a law, A/B it per-condition.

**Matched condition pairs to use** (identical but for capture):

```
gdna_gdna300_ss_0.99_nrna_none_capture_off      vs  gdna_gdna300_ss_0.99_nrna_none_capture_on
gdna_gdna300_ss_0.99_nrna_present_capture_off   vs  gdna_gdna300_ss_0.99_nrna_present_capture_on
gdna_gdna100_ss_0.99_nrna_present_capture_off   vs  gdna_gdna100_ss_0.99_nrna_present_capture_on
```

Restrict to **single-strand** nodes (`free_pos ^ free_neg`) throughout. Keep AMBIG out of it entirely.

**A hypothesis worth testing early, but do not assume it.** `σ²_transfer = Var(log r)` is *tiny* when both
nodes are count-rich (measured 0.0085 on the hop that did the most damage in the last dissection), yet capture
is exactly when a message crosses a real enrichment cliff. The owner's standing intuition applies: *"if we go
off a giant enrichment cliff and lose a hundred x of our signal there's no more reads left there… the
precision has to go down to extremely low levels because there are no more counts left to hold onto that
belief."* Whether the residual is a precision failure of that kind, or a mode failure, is question 1 above —
**measure it, do not assume it.**

## 5. The exact first experiment (copy-paste ready)

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
cd /Users/mkiyer/proj/rigel

# 1. Build the suite-wide node table (all 32 conditions, ~3-5 min, cached substrate).
python scripts/debug/suite_dissect.py --out /tmp/suite_nodes.npz
```

Then, in numpy, on `/tmp/suite_nodes.npz` (columns: `cond node cls dof solvable mass oracle self solved
tau_own c_tau lam_fg cm_g mo_g cm_p mo_p cm_n mo_n nl_cls nl_oracle nl_mass nr_cls nr_oracle nr_mass
n_unspl`), for each matched pair above, restricted to `dof != "ambig"`:

* **`|self − oracle|` off vs on** — if the self-solve degrades, the defect is in `node_init`, not the messages.
* **`|solved − oracle|` off vs on**, and `Δ = solved_err − self_err` — how much do messages ADD under capture?
* Per channel (`cm_g/mo_g`, `cm_p/mo_p`, `cm_n/mo_n`, `c_tau/lam_fg`): does the **mode** move away from the
  oracle under capture, or does the **precision** rise on an already-wrong mode? Those are different bugs.
* Split by region class (`scripts/debug/derive_4_boundary_classes.py` has the classifier) to find which class
  carries the degradation.

Then dissect the worst individual nodes with `scratchpad/dl_dissect.py <COND> --only protected` and
`scratchpad/dump_node.py`, and trace the message provenance hop by hop as in
`density_composition_reconciliation.md` §2.1.

## 6. Tools, files and environment

**Environment.** Everything runs inside the activated conda env; the C++ extension needs `$CONDA_PREFIX`:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1          # REQUIRED for A/B determinism
```

After any C++ change (`src/rigel/native/`): `pip install --no-build-isolation -e .`

**Gates (all must stay green).**

```bash
python -m pytest tests/calibration tests/native -q     # fast loop: 392 pass, 2 xfail, 2 xpass
python -m pytest tests/ -q                             # full: 1229 pass, 2 xfail, 2 xpass
python -m pytest tests/ --update-golden -q             # regenerate goldens — LAST, after the solver is final
ruff check src/ tests/ scripts/                        # must be clean
python scripts/debug/message_variance_mc.py            # the M1-M7 MC arbiter: 0 failures
```

**The benchmark suite** (cached, ~1 s/scenario): `~/Downloads/rigel_runs/ambig_dense_10mb`, 32 conditions,
`_selfsolve_cache`. The A/B:

```bash
P0_REFIT=0 python scripts/debug/pass0_oracle_bench.py --arm myarm_r0
P0_REFIT=1 python scripts/debug/pass0_oracle_bench.py --arm myarm_r1
python scripts/debug/pass0_oracle_bench.py --report --vs pin_r0 --new myarm_r0
# rows accumulate in /tmp/pass0_oracle_bench.tsv; HEAD's arms are `pin_r0` / `pin_r1`
```

**Diagnostics** (`scripts/debug/`): `suite_dissect.py` (suite-wide node table),
`derive_4_boundary_classes.py` (region/boundary census by raw signature bits),
`pass0_oracle_bench.py` (the A/B), `message_variance_mc.py` (the MC arbiter),
`pass0_node_dissect.py` (ψ channel ablation — ⚠ its replay fidelity is stale, see §8),
`seam_lambda_audit.py`. In `scratchpad/`: `dl_dissect.py` (per-scenario, splits error by DL-protection state
and prints message provenance), `dump_node.py`, `derive_1_ratio_check.py`, `derive_2_relay_pin.py` (bit-exact
offline relay replay — very useful, validated `max|Δρ| = 0`), `derive_3_channel_ratios.py`.

**Solver map.** `bp_solver.py` → `node_sweep` → `_unified_solve`: `_relay` (forward/backward, now pinned per
hop), `_transport` (the combine: reframe → route → `_pin_v` → λ gate → DL deflation), `_fuse`, and the three
streams (mode `pg/pp/pn`, measurement `mg/mp/mn`, composition `τ`). Pure laws in `enrichment_frame.py`
(M1–M5 + `mismatch_gap`/`mismatch_deflate`). Self-solve in `node_init.py`. ψ in `simplex_logodds.py`.
Signature bits in `signature.py`. `_capture["_dl"]` publishes the per-message DL gaps.

**Docs (live only).** `ROADMAP.md` (entry point) → this file →
`density_composition_reconciliation.md` (the density↔composition brief, and the full record of what has been
tried and rejected) → `message_variance_derivation.md` (M1–M7) → `CALIBRATION_ARCHITECTURE.md` (the
count-zero-information theory). Handoffs 4/5/6 are the arc. **Never reference `archive/`.**

## 7. ⛔ DO NOT RE-RUN — measured and rejected

| item | verdict |
|---|---|
| **Per-channel enrichment ratios at the boundary face** | Physics real (boundary is 0.125× the exon and 2113× the intron at verystrong) but the fix is **REFUTED**: substituting the ORACLE capture step buys −0.16/−0.17/−0.33 and ≈4 % at verystrong. The relay pin retired it. |
| **The graft-frame fix** (measured mature must not be reframed) | Premise PROVEN (82× over-claim at node 1909; aggregate 0.0964→0.0789) but **REVERTED**: it makes zero-gDNA libraries *more confidently wrong* (mode 0.378→0.492, `Var<0.02` 26.6 %→30.4 %). Needs the hyperprior first. |
| **Moving the graft inside `_pin_v`** | Worse: 9 better / 17 worse. |
| **An RNA "floor" from the spliced measurement** | Not achievable prior-free: exons have `mass_spliced = 0` ALWAYS; at boundaries `corr(log,log)` = 0.405/**0.014**/0.421. |
| **Excluding the spliced from the face `ρ_tot`** | Worse at normal capture (−0.199, −0.106), better only at verystrong. |
| **`σ²_pin` = `(log k)²` as a variance term** | Re-scoped, do NOT write as proposed: `k` sees composition error only through the eff-length-weighted average, so where `E_g ≈ E_r` it is **blind** (max `|log k|` = 0.036). |
| **`k = 1` (disable the pin)**, forcing `r = 1`, graft precision as COUNT, wiring M2's `graft_rna_logvar`, charging σ²_transfer on the graft, peel-before-reframe | All measured worse — see `density_composition_reconciliation.md` §5.4. |
| **Peel frame ordering** | ≤7 % effect, not the driver. |

## 8. Known defects, deliberately NOT fixed (each has numbers)

1. **`f_cur` at unsolvable nodes** (⚠ comment in `bp_solver`). ψ's output at a node with no free RNA strand is
   discarded by the write-back gate but still feeds the NEXT iteration's reframe, and the solver returns 0 for
   a node it never solved — so every gDNA anchor is re-framed as 100 % RNA (`ρ_tot` off by up to 1.8×). The
   one-line fix `np.where(solvable, ..., f_g)` scores slightly **worse** (aggregate +0.0004,
   unstranded-capON 0.1702→0.1740). Something downstream compensates; investigate the pair, do not flip blind.
2. **No TSS/TES in the region/boundary map** (owner). A transcript END is not represented, so the solver models
   the density drop as a capture cliff when the RNA simply stops. Reaches back to index building and the
   accumulator ⇒ deferred. Low expected real-data impact (it surfaced via a simulator artifact: an exon
   interval coinciding exactly with a multi-exonic transcript). NOTE: the `gdna_none` 100 %-RNA boundaries were
   traced to `exon|exon` boundaries, **not** transcript ends — so that evidence does NOT bear on this item.
3. **`coarse_type_array` collapses the retained-intron case** ("exon wins over intron"), and
   `bp_solver.is_exon_node` uses the coarse type — so a region that is one transcript's intron and another's
   exon receives the graft of junction flux that splices OVER it. Latent: RETAINED single-strand currently
   scores 0.0051, identical to a plain exon. Watch it when touching the graft.
4. **`simplex_logodds`'s docstring claim** that a boundary's unspliced crossing is "gDNA + nascent, disjoint
   from the spliced" is **wrong** — it assumes every boundary is an exon→intron junction. At an exon|exon
   boundary RNA is contiguous. Docstring only; no code depends on it.
5. **`pass0_node_dissect.py` replay fidelity is stale.** `_uni_msg` publishes the *mode-fusion* precisions
   (`cpg/cpp/cpn`) while ψ is handed `cm_g/cm_p/cm_n`, and neither `lam_imp` nor `theta_imp` is published.
   Fix it before relying on that tool (publish `cm_*`, `(lam_msg, c_tau)`, `(th_msg, th_prec)`).
6. **The peel's zero-truncation.** A fully-consumed peel emits `t_p = 0` — "no RNA continues past here" — at a
   live precision. `message_variance_derivation.md` §4 says `ρ_ν < 0` should be a PRIOR TRUNCATION, not an
   emission at zero. Peel hops break the mass identity worst (median ×1.31–1.58, 42–53 % over 1.5×).

## 9. INVARIANTS — preserve these

* The `√2·σ_own` pin-safety inequality (§2). Precision is discrete counts over a discrete length; enrichment
  scales the MODE, never the PRECISION.
* `Σ_c ρ_c·E_c = M` is an **identity** under the imputation premise — now enforced at every relay hop AND the
  combine, with `_pin_v` semantics (unsupplied components filled from the node's own density so a PARTIAL
  claim stays partial; rescaling all three blindly regresses capture-OFF 3.6×).
* `N` enters only as power (`τ_λ` Fisher or `1/n` sampling), never as a composition vote.
* Variances in log-odds `Var(λ)` / `Var(log f_c)`, NEVER on simplex fractions.
* The composition is ONE DOF (single-λ message); the tilt θ is a SEPARATE DOF (AMBIG only).
* `v_own` for the DL term is composition-ONLY (no `1/n_dst`): `M_dst` cancels from both sides of the gap.
* Do NOT add DL to `_relay` (no per-hop `_pin_v` there ⇒ it would charge the reframe residual as composition
  mismatch — verbatim the retired `(log r)²` defect).
* Do NOT delete `rna_*_frac_var`, `strand_likelihood.py`, `lam_var`, `var_gdna` — all live.

## 10. Method

Derive → MC-validate (if it is a law) → implement → **per-condition A/B at refit=0 AND refit=1** → dissect any
regression rather than assuming a theory flaw → repeat. Report per-condition better/worse/flat counts, not
just the aggregate; call out `stranded ss_0.99`, `unstranded × capON`, `verystrong` and `gdna_none` explicitly
(each has caught a different failure). Goldens LAST. Multi-agent derivation/audit workflows have worked very
well here (`wf_7d1708d4` caught two real defects; `wf_1ba425be` produced the RNA-over-claim root cause) — use
them for anything with more than one plausible mechanism.

---

## 11. ▶ KICKOFF PROMPT (copy-paste to start the next session)

> We are developing Rigel calibration — the prior-free first pass, "pass-0". Read
> `docs/calibration/ROADMAP.md`, then `docs/calibration/SESSION_2026_07_25_HANDOFF_7.md` (START HERE — it has
> the state, the census, the exact next experiment, and the do-not-re-run list), then
> `docs/calibration/density_composition_reconciliation.md` and `message_variance_derivation.md`. Do NOT read
> `docs/calibration/archive/`.
>
> Setup: `source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel`, always
> `OMP_NUM_THREADS=1`, repo `/Users/mkiyer/proj/rigel`, branch `calib-ambig-init-wip` at `d88a7822`. Gates:
> `pytest tests/ -q` = 1229 pass / 2 xfail / 2 xpass, `ruff check src/ tests/ scripts/` clean,
> `python scripts/debug/message_variance_mc.py` = 0 failures. Benchmark suite (cached, ~1 s/scenario):
> `~/Downloads/rigel_runs/ambig_dense_10mb`, 32 conditions; the A/B is
> `P0_REFIT=0|1 python scripts/debug/pass0_oracle_bench.py --arm NAME` and HEAD's arms are `pin_r0`/`pin_r1`
> (aggregate 0.0926 refit=0, 0.0779 refit=1). Untracked `scratchpad/` holds the derivation scripts — keep it.
>
> **THE TASK: the pass-0 SOLVER must be CORRECT before we fit the gDNA hyperprior** (owner's directive — do not
> propose moving to the hyperprior). A suite-wide census (HANDOFF_7 §3) found that off-capture, single-strand
> pass-0 is essentially solved (mwae **0.0051**; exons 0.0046, retained introns 0.0051, intergenic 0.0000, and
> exon|exon boundaries 0.0305/0.0477 — the structural cases are all handled). Two axes remain: **capture
> degrades single-strand 10× (0.0051 → 0.0476)** and **AMBIG is ~50 % of all error mass** (concentrated on
> `exon+intron(x-strand)`, a region that is one strand's exon and the other's intron). **Debug single-strand
> first, then AMBIG.**
>
> **So: run the single-strand × capture study (HANDOFF_7 §4–§5).** It has a clean control — the same nodes and
> annotations score 0.0046 off-capture and 0.0499 on-capture, so anything you find is attributable. Build the
> node table with `python scripts/debug/suite_dissect.py`, restrict to single-strand (`dof != "ambig"`), and on
> the matched capture-off/on condition pairs decompose *in this order*: does capture break (a) the node's own
> message-free self-solve, (b) the message MODES, or (c) the message PRECISIONS? Then localize to a region /
> boundary class (`scripts/debug/derive_4_boundary_classes.py`) and then to individual nodes
> (`scratchpad/dl_dissect.py`, `scratchpad/dump_node.py`), tracing message provenance hop by hop. Only then
> propose a fix: derive it, MC-validate it if it is a law, and A/B it per-condition at refit=0 AND refit=1.
>
> Vocabulary: RNA is just RNA — "mature"/"nascent" are not different species, only SPLICED vs UNSPLICED is
> observable; a boundary can be an exon↔exon boundary that is ALSO a splice junction, with RNA contiguous
> across it while other RNA splices in or out. GOVERNING PRINCIPLE: pass-0 must be WEAK and CORRECTABLE — an
> over-confident message that pins a node wrong is worse than a weak one slightly off; a node without the data
> to solve itself must end at near-zero precision or unsolved, never moderate precision; prefer the
> under-confident option when unsure. No magic numbers. Counts are Poisson. **Check HANDOFF_7 §7 before trying
> anything — per-channel enrichment ratios, the graft-frame fix, the RNA floor, `σ²_pin` and six other ideas
> are already measured and rejected.** Regenerate goldens LAST.
