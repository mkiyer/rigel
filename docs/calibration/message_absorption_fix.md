# Message-propagation fix — the spliced-absorption bug (detailed next steps)

**Status:** DIAGNOSED + owner-confirmed framework; NOT yet implemented. Written 2026-07-18. This is the immediate
next work item — a set of message-layer bugs whose fix is expected to substantially improve the solver.
**Read with:** `CALIBRATION_MASTER.md` (§3 sources, §5 the two NPMLE roles), `CALIBRATION_ARCHITECTURE.md` §0.

---

## 1. The bug, in one paragraph

At a **splice-junction boundary** (intron↔exon), the per-strand RNA imputation message
(`bp_solver._scan`, the `rho_pos`/`rho_neg` block, ~L430–455) subtracts the **destination boundary's spliced
count** from the **source's unspliced-RNA density**:

```
rho_pos = fp_s·sm/er   +   SPs[lsrc]/esp   −   SPd[i]/ESPd[i]        (the "− dst absorption" term)
mo      = log( max(rho_pos, 1/erd) / (md/erd) )                      (the clamp)
```

These two quantities are **disjoint deposits at different places and scales** — the source's *facing unspliced*
RNA vs. the boundary's *spliced* count — so the subtraction has no reason to stay positive. When it goes
negative, the `max(·, 1/erd)` **clamp launders it into "essentially no RNA,"** and — worst of all — sends that at
**full precision** (`pr` is computed from counts + belief variance, *not* from the fact that the content
cancelled). A strand-flat unstranded node reads "no RNA" as "all gDNA" and snaps to f_g≈0.97, then propagates a
confident gDNA message to its neighbours — the **boundary-dominated confident-false-positive cascade** (420 of
666 confident-wrong nodes in gdna5).

---

## 2. The worked example — the real counts (gdna5, chain nodes 2929/2930/2931)

Structure: **INTRON+ (2929) — boundary/splice-acceptor (2930) — EXON+ (2931)** of expressed gene G0299.

Oracle counts (from the accumulator; `mat_uns`/`nas_uns`/`mat_spl` = mature-unspliced / nascent-unspliced /
mature-spliced):

| node | unspliced total | gDNA | mature-unspl | **nascent-unspl** | **spliced (mature)** | correct unspliced f_g |
|---|---|---|---|---|---|---|
| INTRON 2929 | 78 | 0 | 0 | 78 | 0 | 0.000 |
| EXON 2931 | 119,021 | 493 | 73,639 | 44,889 | 0 | 0.004 |
| **BOUNDARY 2930** | **962** | **7** | **0** | **955** | **7,502** | **0.007** |

**What is physically true at the boundary:** two distinct RNA pools coexist — 7,502 **mature that splices out**
(the spliced reads, the mature leaving via the junction) and 955 **nascent that crosses unspliced** (pre-mRNA
still being transcribed). The unspliced crossing is 955 nascent + 7 gDNA ⇒ **correct f_g = 0.007** (essentially
pure RNA). Note the exon's unspliced is 73,639 mature + 44,889 nascent — **indistinguishable** (a read inside an
exon looks unspliced whether mature or nascent). We do **not** know the mature/nascent split in the exon.

**What the code does** (exon→boundary "backward" message, verified with the real geometry):
```
rho_pos =  exon unspliced-RNA density 41.58   +  0   −  boundary spliced density (7502/100) 75.02  =  −33.44
        →  clamp → mode −6.3 → "no RNA" → f_g = 0.97   (truth 0.007)
```
The eff-lengths are **consistent** (er = erd = ESPd ≈ 100) — so this is **not** a units mismatch. It is a
**disjoint-pool subtraction**: the boundary's *entire* spliced pool (75) exceeds the source's *facing-window*
unspliced RNA (41.58), because the spliced piles up at the junction while the facing unspliced is one
fragment-length window. Subtracting a total pool from a window density is meaningless and goes negative.

---

## 3. The correct model (owner-confirmed — the framework to implement)

> **The exon region can only report what it HAS; the boundary owns the splicing information and does the
> subtraction from its own side.**

1. **Region → boundary message = the region's current densities `(ρ_RNA+, ρ_RNA−, ρ_gDNA)`, NO spliced
   subtraction.** The exon says "this is the RNA and DNA density I have." It does not know what splices next.
   (Crucially: we do **not** infer nascent from abundance or from the spliced count — there is no biological
   correlation between expression level and nascent fraction; splicing efficiency and nascent half-life are
   independently regulated. **We can only measure.**)
2. **The boundary receives that density and does the partition using its OWN effective length and its OWN
   measured spliced count.** It converts the incoming RNA density to a count over its eff-length, subtracts the
   **spliced** it directly measured (the mature that left via the junction), and whatever remains **continues as
   nascent** unspliced RNA. The boundary's unspliced mass then partitions into **nascent (the residual RNA) vs
   gDNA** — the composition we solve. The spliced count is a *positive, direct RNA measurement* the boundary
   owns; it is never used to reduce the RNA imputation of a *neighbour*.
3. **Boundary → onward message propagates the NASCENT** (the unspliced RNA that continues), i.e. the spliced
   subtraction lives on the boundary's **outgoing** side, because the mature that spliced does not continue.
4. **Honest clamp (independent, land first):** a cancelled / non-positive imputation carries NO information and
   must be emitted as **no message (precision 0)** — never a confident value. The clamp must never manufacture
   confidence. (This is the "no clamps/cliffs" principle applied to messages.)

**The acceptance test for the fix, on this node:** the boundary 2930 must solve to **f_g ≈ 0.007** (nascent ≈
955, gDNA ≈ 7 of the 962 unspliced), NOT 0.97.

---

## 4. Implementation plan (specific, ordered)

**Step 0 — honest clamp (smallest, unambiguous, land + measure first).** In `bp_solver._scan`, when the
imputed `rho_pos`/`rho_neg` is ≤ 0 (or ≤ the count floor `1/erd`), **do not emit the RNA factor** — set its
precision to 0 and skip it (do not append to `lam_factors`, leave `amp/app = 0`). This removes the laundering
regardless of the deeper fix. Re-run the confident-FP A/B (`scripts/debug/confident_fp_trace.py`,
`hyperprior_fit_options.py`) and record how much of the cascade it removes on its own.

**Step 1 — remove the cross-node spliced subtraction from the INCOMING message.** Delete the
`− SPd[i]/ESPd[i]` (`absorb_p`) and `− SNd[i]/ESPd[i]` (`absorb_n`) terms — the destination's spliced must not
subtract from a neighbour's RNA imputation. The `_RNA_ABSORB` module toggle (added this session) already gates
these for the A/B; make the removal permanent once verified. Also re-examine the `+ SPs[lsrc]/esp` source-spliced
*addition*: a boundary source's spliced is the mature that left; whether it belongs in the *unspliced* RNA
message it emits is exactly the Step-2 question.

**Step 2 — move the spliced partition to the BOUNDARY's own side, in consistent units.** The boundary owns its
directly-measured **spliced** and **unspliced** counts. Implement the partition there: the RNA that continues as
unspliced (nascent) = (incoming RNA density × the boundary's RNA eff-length) − (the boundary's spliced count),
all in **counts on the boundary's eff-length basis**, floored at 0. The residual unspliced mass is gDNA vs
nascent. The boundary's *outgoing* RNA message then carries the nascent (post-subtraction). **The critical detail
is the eff-length basis** (see §5) — the incoming density and the spliced count must be reconciled to the same
support so the subtraction is physically meaningful (the acceptance test in §3 is the check).

**Step 3 — re-verify** on the node set (§2), the confident-FP A/B, and the full gradient
(`hyperprior_fit_options.py`) — the hyperprior should stop hallucinating gDNA at expressed splice junctions.

---

## 5. The critical eff-length reconciliation (the piece to get exactly right)

The A/B in §2 shows the subtraction goes negative because the incoming RNA density (facing window, ρ≈41.58) and
the spliced count (total junction pile-up, 7,502 over ESPd≈100 ⇒ ρ≈75) are **not on the same support**. The fix
in Step 2 must convert both to counts on the **boundary's own eff-length** before subtracting, so that on this
node the arithmetic yields nascent ≈ 955 and gDNA ≈ 7 (f_g ≈ 0.007). Concretely, the next session must:
- confirm what `mass_left/right`, `eff_rna_left/right`, `eff_spl_left/right`, `spliced_pos/neg_left/right`
  (the `NodeGeometry` faces) actually represent in counts vs density, using the oracle counts of §2 as ground
  truth (e.g. the boundary's spliced eff-length and the RNA eff-length must be such that the measured spliced
  =7,502 and the imputed total RNA reconcile to nascent≈955);
- derive the subtraction so it is exact on this node **and** degenerate-safe (no spliced ⇒ no change; pure intron
  ⇒ all unspliced RNA is nascent). Do NOT introduce a tuned constant — this is geometry, not a knob.

---

## 6. Tooling / reproduce

- `scripts/debug/confident_fp_trace.py --condition gdna_gdna5_ss_0.50_nrna_present_capture_on` — the
  confidently-wrong trace (per-node strand→loc→final, gDNA/RNA messages, neighbours). `--nodes 2930,2929,2928`
  traces the specific seed. `--dir over|under|both`.
- `scripts/debug/npmle_smoothness_diag.py` — the false-positive smear + bandwidth panels.
- `scripts/debug/hyperprior_fit_options.py` — the gradient plot (oracle-vs-fit); the ultimate downstream check.
- `bp_solver._RNA_ABSORB` (module toggle, default True) — the A/B switch for the absorption term.
- Oracle counts: `_scan_and_truth(...)['region_pools' / 'boundary_pools']` keyed by region/boundary id
  (= chain `ref_idx`); keys `gdna_pos/neg`, `mat_uns_pos/neg`, `nas_uns_pos/neg`, `mat_spl`, `nas_spl`.
- Env: `rigel` conda env, `OMP_NUM_THREADS=1`; caches at
  `~/Downloads/rigel_runs/ambig_dense_10mb/_selfsolve_cache`.

## 7. Why this matters (the expected impact)

The confident false positives are boundary-dominated and cascade; they are what corrupt the deconvolved-gDNA
hyperprior fit (they smear ~40% false positives across the density axis — `npmle_smoothness_diag.py`). Fixing
the RNA imputation at splice junctions removes the *seed* of the cascade, which should (a) collapse the
confident-FP population, (b) let the Phase-2 gDNA hyperprior be fit on clean deconvolved gDNA, and (c) improve
the unstranded low/zero-DNA corners without any bandwidth or strength tuning. It is a correctness fix, not a
tuning fix — exactly where the leverage is.
