# Roadmap item 1 — splice-junction absorption: the RNA sink

> ## ⚠ SUPERSEDED — do not implement §3.2 or §5 step 1
>
> **→ [`mature_crossing_gate.md`](mature_crossing_gate.md) is the live document.** Written 2026-07-16 after an
> adversarial review measured this one against the cached 7-condition oracle substrate.
>
> **The diagnosis (§1–§2) STANDS and is credited** — exon↔intron seams carry exactly 0 crossing mature (0 of
> 1,146, reproduced on 7/7 conditions) and the RNA message fires on 100% of them. That finding is the basis of
> the replacement.
>
> **The prescription is retired:**
> * **§3.2 is an exact algebraic no-op.** `rho_mature_est = SPd/ESPd ≡ rho_mat_dst`, so `rho_new ≡ rho_cur`.
>   Measured: `max|Δmode| = 4.263e-14`; **0 of 579** floored modes differ (the shipped `max(rho, 1.0/erd)` floor
>   dominates the proposed `max(·,0)` clamp). Three independent end-to-end implementations gave intron mwae
>   `0.1573 → ~0.1517` against §4's predicted **~0.0098**. §5 step 2's acceptance criterion is unachievable by
>   construction (`Δmode = 0.00e+00`).
> * **Its only real change (`pr` from the nascent count) is a regression** — it drops `n_mat` from the precision
>   on all 679 B→exon edges, degrading the only strand-free evidence in the system.
> * **The estimator §3.1 rests on is noise**: given an *oracle-perfect* source belief, corr = 0.311, additive
>   bias +837 (4.4× the true signal), SNR 0.13.
> * **§5 step 1's `coarse_type_array` taxonomy is strand-blind** — wrong on 181/909 coarse-EXON regions on the
>   very suite it validates against. Use the per-strand `mrna_active_{pos,neg}`, which already exists.
> * **Two headline numbers do not reproduce**: `rho_spliced/rho_mature p50 = 1.185 (n=614)` measures
>   **0.962 (n=599)** across nine density-definition variants (provenance unknown); §1's `0.0098` disagrees with
>   its own cited tool (`selfsolve_diag`: 0.006535 — different estimators; `0.0098` is the correct basis).
> * **The 20× is 1 of 7 conditions** (19.71 vs 1.20–8.95). Under capture the attribution **inverts** — gDNA does
>   the damage. *Not* cherry-picked: on the honest increment metric this cell ranks 4 of 7.
> * **B740 is not a "p5 outlier" to defer to item 2** — 7 of the 10 worst-damaged introns have `SPd = 0`, where
>   §3.2 has no effect at all. And `boundary_spliced_channel_design.md` §4.2 already names the mechanism
>   (alternative splicing: `n_mat/esp` measures only the isoforms using *that* junction).
>
> Retained below unchanged for provenance.

**Status:** ⚠ **SUPERSEDED** (see above). Written 2026-07-16, branch `calib-ambig-init-wip`.
**Companion:** [`dof_pie_model_fix.md`](dof_pie_model_fix.md) — item 2. **Read its §1 before A/B-ing this one:**
item 2 is what *bounds* this one's residual, so item 1 alone will look better than it is (§7).
**Reproduce everything here:** `scripts/debug/intron_message_trace.py`, `scripts/debug/selfsolve_diag.py`.

---

## 1. The symptom

`selfsolve_diag`, `gdna300_ss0.99_nrna_none_capture_off`, mass-weighted `|Δf_g|` over intron regions:

| | intron mwae |
|---|---|
| message-free self-solve | **0.0098** |
| + gDNA message only | **0.0093** |
| + **RNA message only** | **0.1932** — 20× worse |
| both (shipped) | 0.1573 |

**The gDNA channel is innocent; the RNA channel does 20× damage.** gDNA even partially *rescues* it
(0.1932 → 0.1573) — the correct message fighting the wrong one. Ablation is by replaying the final solve with
each channel's `(mode, prec)` suppressed; everything else identical.

Per-node, the 10 worst-damaged introns are **all** `B->exon|B->exon`, and all show the same shape:

| node | true f_g | f_loc | f_fin | RNA mode | RNA prec | gDNA mode | gDNA prec |
|---|---|---|---|---|---|---|---|
| 741 | 1.00 | 0.9937 | **0.0319** | **+3.096** | **1565.8** | −0.083 | 46.4 |
| 2367 | 1.00 | 0.9899 | **0.0296** | **+3.634** | **1860.7** | 0.007 | 46.4 |

`mode = log(target f_RNA)`, so `+3.4` demands `f_RNA = 30` — **3000%**, at 1600 pseudo-observations. The local
solve had these at 0.98–0.995 against a truth of 1.00.

---

## 2. THE DIAGNOSIS — the taxonomy, not the frames

⚠ **A previous diagnosis in this conversation ("the absorption subtracts in the wrong frame, 188× off") is
RETRACTED.** It was built on node B740, which measurement shows is a **p5 outlier**. The frames are correct.

Every exon-flanked boundary, classified by flanks × junction, against the oracle:

| flanks | junction? | n | RNA msg gated ON | **mature ACTUALLY crossing unspliced** |
|---|---|---|---|---|
| exon–intron | **YES** | 913 | 100% | **0** |
| exon–intergenic | no | 430 | **0%** | 0 |
| exon–intron | **no** | 233 | 100% | **0** |
| exon–exon | YES | 89 | 100% | 20,943 |
| exon–exon | no | 32 | 100% | 23,831 |

**What is already right:**
* **TSS/TES** (exon–intergenic, 430 seams): gated **0%**. Mature ends there and the message correctly never
  fires. `free_s = nrna_active_s(left) & nrna_active_s(right)` does its job.
* **exon–exon** (121 seams): gated ON, and mature genuinely **does** cross unspliced (44,774 mass). Correct —
  contiguous exon body, nothing splices.
* **The eff-length frames**: `rho_spliced(junction) / rho_mature(exon)` over 614 genuine junction faces is
  **p25=0.985, p50=1.185, p75=1.506**. If `esp` were wrong by 188× this ratio could not be ≈1. The junction's
  spliced density **does** recover the adjacent exon's mature density ~1:1.

**What is wrong — one row, 1,146 seams:**

> **exon–intron seams carry EXACTLY ZERO crossing mature (0 of 1,146), and the RNA message fires on 100% of
> them, carrying the exon's full mature density.**

Mature RNA never crosses an exon–intron seam. It **splices** (913 junction seams) or the transcript is simply
unspliced there (233). Either way it does not enter the boundary's unspliced channel. The message asserts it
does.

### 2.1 The arithmetic of the failure (B740, the traced case)

```
exon 739:  975 gDNA | 29,969 mature | 0 nascent   -> rho_unspl = 17.794   (fg_loc = 0.0296, correct)
B740    :   67 gDNA |      0 mature | 0 nascent   -> rho_unspl =  0.670   (fp_loc = 0.1071, correct)
                        9.7 SPLICED
intron741: 227 gDNA |      0 mature | 0 nascent   -> rho_unspl =  0.563   (fg_loc = 0.9937, correct)
```

Every node self-solves correctly. Then the message from exon 739 says *"your RNA density is 18.3"* to a node
whose **entire** unspliced density is **0.670**:

```
f_RNA x 0.670 = 18.3   =>   f_RNA = 27.3      IMPOSSIBLE
```

`f_pos` saturates at 0.984, `f_g` → 0.016. B740 — now believing it is 98% RNA — then relays that into intron
741 at precision 1565. **The traced 740→741 message was one hop downstream of the real event.**

The true gDNA density is **uniform** across this chain (0.561 / 0.670 / 0.563 — exon / boundary / intron), so
the geometry is right and the gDNA channel correctly sees no discrepancy at all. Only the RNA channel does,
and only because it is carrying a molecule that isn't there.

### 2.2 Why the existing absorption does not save it

```python
rho = fbp[lsrc]*sm/er  +  SPs[lsrc]/esp  -  SPd[i]/ESPd[i]
```

* On an **exon→B** edge (the sink), `SPs[lsrc] = 0` — `spliced_*` is **identically zero on every region**
  (`node_geometry`: `pick(zeros_R, b_spl_*)`; verified: total spliced mass over all region nodes = 0.0000 on
  both faces). So the mature *source* term never fires here.
* `SPd[i]` is the **destination boundary's** spliced mass, and it **does** fire — B740: `rho_abs = 0.0974`
  against `rho_msg = 18.30`. **0.53% absorbed.**
* On a **B→region** edge, `SPd[region] = 0` always ⇒ `rho_abs ≡ 0`. **The absorption is structurally
  unreachable into a region** — which is *correct* (the mature already left at the region→B edge), but it
  means nothing downstream can undo a boundary that has already been corrupted.

So the absorption fires on the right edge and subtracts a commensurable quantity — it simply cannot cover a
message that should never have been sent. **On 83.6% of genuine junction faces `rho_spl` covers ≥80% of
`rho_mat`, and after subtraction the required `f_RNA` is p25=0.000 / p50=0.000 / p75=0.003.** The subtraction
is the correct arithmetic. The 1,146 exon–intron seams are not an arithmetic problem; the message content is
wrong at the source.

---

## 3. THE FIX

**An exon's RNA message must carry its NASCENT only. Its mature must never enter the message.**

Three cases, decided by the seam's own structure — no new constant, no tuning:

| seam | mature's fate | RNA message should carry |
|---|---|---|
| **exon–intron** (913 junction + 233 not) | splices out, or simply does not cross | **nascent only** |
| **exon–exon** | crosses unspliced, contiguously | **nascent + mature** (unchanged — already correct) |
| **exon–intergenic** (TSS/TES) | ends | **nothing** (already gated off) |

Only the first row changes.

### 3.1 The mature estimate, and where it is available

At a **junction** seam the boundary's own spliced density measures the adjacent exon's mature density —
measured, median ratio **1.185** (§2). So the exon's mature is *estimable at the seam*, and subtracting it in
the exon's own frame leaves the nascent:

```
rho_msg_nascent = max( fbp[lsrc]*sm/er  -  rho_mature_est ,  0 )
rho_mature_est  = SP_face(B)/ESP_face(B)          # the seam's own spliced density, SAME frame (ratio ~1.19)
```

This is what the current `- SPd[i]/ESPd[i]` term already computes. **It is already in the right place and the
right frame** — see §4 for why it is nevertheless not doing the job.

At a **non-junction exon–intron** seam (233) there is no spliced count, hence no mature estimate — but there
is also nothing to estimate: the transcript is unspliced across that seam, so any RNA crossing IS nascent, and
the message is already correct. *(Measured: 0 crossing mature on all 233.)*

### 3.2 What actually has to change

The failure is not the formula, it is **which node's belief the message is built from**. `fbp[lsrc]` is the
source's **RUNNING** belief, and by the time the forward scan reaches boundary 740 it has already absorbed
exon 739's message and been flipped from 11% RNA to ~100% RNA. It then relays *that* onward. The absorption
subtracts the seam's mature *once*, at the exon→B edge; it cannot undo a belief that the same message already
corrupted.

**⇒ The exon→B absorption must reduce the SOURCE's contribution before it enters `fbg/fbp`, not merely the
message that leaves.** Concretely, in `bp_solver._scan`'s `emit_p` / `emit_n` blocks:

```python
# CURRENT (the absorption is applied to the outgoing rho only):
rho = n_nasc/er + n_mat/esp - rho_mat_dst
mo  = log(max(rho, 1.0/erd) / (md/erd))
pr  = n_src/(n_src*vbp[lsrc] + 1.0)
# ... then fbp[i] is updated from (mo, pr)

# REQUIRED: on an exon->B seam the mature must be removed from the message BEFORE the dst folds it in,
# and n_src must be the NASCENT count -- not the total RNA count -- or the precision inherits the mature.
n_mat_est = rho_mature_est * er            # the seam's mature, expressed in the SOURCE's frame
n_nasc    = max(fbp[lsrc]*sm - n_mat_est, 0.0)
rho       = n_nasc/er + n_mat/esp
mo        = log(max(rho, 1.0/erd) / (md/erd))
pr        = n_nasc/(n_nasc*vbp[lsrc] + 1.0)   # <-- precision from the NASCENT count
```

The second line of the fix matters as much as the first: `pr = n_src/(...)` with `n_src = fbp*sm` = the
**total** RNA count gives the exon a 27,209-pseudo-observation voice (measured at B740) built on mature it
should not be sending. Removing the mature from `rho` but leaving it in `n_src` would leave a quiet-but-still-
overconfident message.

---

## 4. Why this will not fully work on its own — and the number to expect

**9.4% of genuine junction faces demand an IMPOSSIBLE `f_RNA > 1` even after a correct subtraction
(max 41.5).** A composition cannot exceed 1. That the arithmetic can *ask* for 41.5 is item 2's defect — the
message is not constrained to be a composition, so nothing bounds it.

And **59.0% of junction faces have `rho_mat` alone exceeding the boundary's TOTAL capacity** (median
`rho_mat/rho_boundary_total` = **1.99**). A *fraction* of an impossible number is usually still impossible —
which is why a proportional-allocation rule (allocate the incoming RNA by the seam's own spliced:unspliced
split) does not rescue the failures either: at B740 it gives `18.3 × (3.4/13.1) = 4.76` ⇒ `f_+ = 14.9`, still
impossible.

**Predict, before measuring:** item 1 fixes the 1,146 exon–intron seams and the intron mwae should fall from
**0.1932 toward ~0.0098** (the message-free floor). A residual tail will remain on the ~4–16% of seams whose
junction does not measure the adjacent exon's mature (`rho_spl/rho_mat < 0.10` on 4.4%; B740 is one). **Do not
tune anything against that tail — it is item 2's.**

---

## 5. Implementation plan

1. **Seam taxonomy helper** (`node_geometry`): per boundary face, classify `exon–exon` / `exon–intron` /
   `exon–intergenic` from the flank signatures (`coarse_type_array`) — the same source `_boundary_strand_stats`
   already uses for `free_s`. No new data, no constant. *Verify: reproduce the §2 table exactly (913/430/233/89/32).*
2. **Mature removal at the source** (`bp_solver._scan`, `emit_p`/`emit_n`): on an **exon–intron** seam,
   subtract the seam's mature estimate from `n_nasc` **and** use the nascent count for `pr` (§3.2). Leave
   exon–exon and TSS/TES untouched. *Verify: `intron_message_trace` — the RNA mode on the 10 traced introns
   must fall from +3.1…+3.6 to ≤0, and `prec_p` from ~1600 to O(n_nascent).*
3. **Guard the invariant**: a unit test asserting **zero** crossing mature on exon–intron seams in the oracle
   (`0 of 1,146` — the fact the whole fix rests on), and that the RNA message on such a seam is 0 when the
   source exon has no nascent. *This is the test that would have caught the bug.*
4. **Re-run `selfsolve_diag --stage both`** on capOFF and capON. Expect intron mwae → ~0.01 on capOFF; expect a
   residual on the 4–16% tail.

**Do NOT** ship a proportional-allocation rule or tune `esp`. Both were considered and are refuted above
(§3, §4): the frames are right (median 1.185) and scaling an impossible message leaves it impossible.

---

## 6. What is established vs assumed

| | |
|---|---|
| gDNA channel innocent; RNA channel does 20× damage | **measured** (ablation, §1) |
| exon–intron seams carry 0 crossing mature; message fires on 100% | **measured** (oracle, 1,146 seams) |
| TSS/TES already gated off (0%); exon–exon correctly on | **measured** |
| junction spliced density recovers exon mature density ~1:1 | **measured** (median 1.185, n=614) |
| the eff-length frames are correct | **measured** (follows from the above) |
| the source's RUNNING belief is what carries the mature onward | **traced** (B740: fp_loc 0.107 → f_pos 0.984) |
| removing mature from `n_src` fixes the precision inflation | **assumed** — verify at step 2 |
| the residual tail is item 2's, not item 1's | **assumed** — the 9.4%-impossible figure supports it, but it is not proof |

---

## 7. The one thing to remember when A/B-ing this

**Item 2 bounds item 1's residual.** Item 1 alone will fix the bulk and leave a tail that still saturates,
because the representation permits `f_RNA = 41.5`. If we ship item 1 and measure, we will see a large win and
an unexplained remainder — that remainder is *predicted*, and it is not a reason to tune item 1.
