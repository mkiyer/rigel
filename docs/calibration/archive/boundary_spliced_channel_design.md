# Priority #3 — the boundary spliced channel: "mature is a MEASUREMENT, nascent is an IMPUTATION"

**Status:** design, **REVIEWED — §4.1/§4.2 decided, §4.3 CLOSED. One dependency remains (§4.2).**
Written 2026-07-15, branch `calib-ambig-init-wip`. Review outcome + what changed: **§7**.
**Prereqs landed:** #1 the derived reference (`reference_prior_derivation.md` §10); #2 honest message
precision (`pr = n_src/(n_src·vb_src + 1)`, σ²_transfer zeroed).
**Roadmap position:** step 3's gate. `selfsolve_diag.py` shows every solvable node at the no-information
value when unstranded; this is the only channel that can change that.

> **Do not implement §3 as written.** The diagnosis (§1) is correct and is the highest-value lever available.
> The *estimator* I first proposed is wrong in two ways, one of which I stated to the user as a headline
> result and which is **retracted** here (§2.3). §5 lists what is safe to land today.

---

## 1. The diagnosis (verified, not disputed by any reviewer)

`bp_solver._scan`'s RNA message already carries spliced mass. It is **lumped**:

```python
n_nasc = fbp[lsrc] * sm      # IMPUTATION: the source's BELIEVED unspliced RNA count
n_mat  = SPs[lsrc]           # MEASUREMENT: the source-face spliced mass (>0 only B->exon)
rho    = n_nasc/er + n_mat/esp - rho_mat_dst          # <-- ONE density
n_src  = n_nasc + n_mat
pr     = n_src / (n_src * vbp[lsrc] + 1.0)            # <-- ONE precision, from the NASCENT belief
```

So a **direct count of motif-stranded fragments** is delivered at the *nascent guess's* precision. Measured on
`gdna300_ss_0.50_nrna_none_capture_off`, boundaries with spliced mass (1015 / 1699):

| quantile | spliced mass | `pr` today (lumped) | `pr` if the mature channel carried its own count | ratio |
|---|---|---|---|---|
| p10 | 2 | 0.354 | 2 | 5× |
| p50 | **62** | **0.354** | **62** | **176×** |
| p90 | 1,141 | 0.355 | 1,141 | **3,215×** |

The precise, **strand-free** channel is muted by the imprecise one, at κ=½ where it is the *only* evidence in
the system. (`pr` is pinned at `1/2.804 = 0.357` — the reference's own `Var(log f_g)`.)

**Why this is strand-free:** `bam_scanner.cpp:1451,1483` — `deposit_strand = spliced ? motif_strand :
cr.align_strand`. The spliced channel's strand comes from the **splice motif**, not the library protocol.
Nothing consults the strand model.

### 1.1 What survives review (do not weaken these under pressure)

* **The frame is right.** `ρ` is an *abundance*, not a fragment density — the eff-length divisor converts a
  face's fragment mass to copies. `n_mat/esp` is `A_mature`; `n_nasc/er` is `A_nascent`; the destination
  exon's unspliced RNA *is* `A_mature(body) + A_nascent`. Two divisors for two geometries
  (`effective_length.py:135-161` half-triangle vs `:57-73` contained) is exactly what makes them the same `ρ`.
  **The "a junction's spliced count can't measure the exon body" objection dissolves** — they are the same
  molecules sampled by two windows.
* **Count-zero-information is NOT violated.** The evidence is the *classification* (CIGAR-N ⇒ not gDNA),
  structurally identical to a strand label; the count enters only as that label's Fisher information — the
  principle's own allowed role.
* **`n_nasc` and `n_mat` share no fragments.** ψ has no spliced term (`simplex_logodds` docstring), so
  `fbp[lsrc]` is spliced-free. *Corollary, worth saying out loud:* because ψ ignores a node's own spliced
  mass, **this message is the only route by which spliced evidence reaches an exon.** That raises the stakes
  on its variance; it does not lower them.

---

## 2. Why the estimator I proposed does not ship — four blockers

| # | Blocker | Class |
|---|---|---|
| **B1** | `Var(ρ_m) = ρ_m²/n_mat` → **NaN** at `n_mat = 0`, which is **≥50% of all edges** | mechanical (§3.1) |
| **B2** | `n_mat` is fractional **mass**, not a Poisson count | mechanical (§4.1) |
| **B3** | **The delta method is anti-conservative in the target regime** | real (§2.3, §3.2) |
| **B4** | The highest-precision channel in the system gets **no transfer term** | **decision** (§4.2) |

### 2.1 B1 — the NaN is the modal edge, and it is silent

The chain is `B0 R0 B1 R1 … Bk` (`node_chain.py:93-101`), so exactly half of all edges have a REGION source;
`SP` is identically zero on regions (`node_geometry.py:174-178`). Verified on numpy 2.4.3:

```
rho_m**2 / n_mat   with both 0.0   ->  nan       # silent; no exception; propagates into vbp
m_mat / esp**2     with m_mat=0    ->  0.0       # exactly 0 — the correct limit
```

Dissolves under the §3.1 rewrite. Not a repair — the *same algebra*, written so it can't divide by zero.

### 2.2 B2 — `n_mat` is mass, not a count

`node_geometry.py:125-126` reads `bsub.left.mass_spliced` — fractional mass (`accumulator.cpp:204`,
`share = (sl.end-sl.start)*inv_L/n_cross`). So `1/n_mat` is **not** a Poisson variance. Measured mass/flux =
**0.39 mean, 0.078 at p10**.

**But the reviewers disagreed on severity and the dissent is right on the math:** for compound-Poisson mass
`Σwᵢ` with `w ∈ (0,1]`, `CV² = 1/n_eff` with the Kish `n_eff = (Σw)²/Σw² ≥ Σw = mass`. So **`1/mass`
*overstates* the variance — a valid conservative bound, not a correctness bug.** It is 2.6× conservative on
average, 13× at p10. See §4.1 for the options; an integer flux **does** exist
(`boundary_flux_{left,right}[b, 2:4]`, uint32) and is already used as a count elsewhere in the package.

*Note: the mass is **correct** for the MODE either way — `ρ_m = m_mat/esp` must use the frame `esp` is defined
in. Only the variance wants the count.*

### 2.3 B3 — RETRACTION: my headline claim was wrong

I told the user: *"when `ρ_m ≫ ρ_n` the combined precision → `n_mat`, full count precision."* **That is false
as stated.** The delta method (`Var(log ρ) = Var(ρ)/ρ²`) linearizes at the mean and discards the upper tail of
a lognormal with log-sd `√2.804 = 1.67` — a fat tail. Even a 1% nascent addend moves the sum far more than the
linearization sees. `ρ_m ≫ ρ_n` is not enough; the claim needs `w = ρ_m/ρ ≥ 0.999`, i.e. **nascent below 0.1%
of the density**.

**⚠ The magnitude depends on a convention the code has not fixed** (§4.3). Both conventions agree on the
*direction*; they differ ~4× on how bad it is. Monte Carlo, 2M draws, `v_m = 1/500`, `v_n = 2.804`:

| convention | w | exact | **delta (mine)** | **FW** | delta/exact |
|---|---|---|---|---|---|
| mean-matched | 0.90 | 32.62 | 33.72 | 6.87 | **1.03×** |
| mean-matched | 0.99 | 357.83 | 446.31 | 285.14 | **1.25×** |
| median-matched | 0.90 | 6.87 | 33.72 | 1.09 | **4.90×** |
| median-matched | 0.99 | 107.66 | 446.31 | 39.09 | **4.15×** |

Delta is anti-conservative under both. FW is conservative under both.

---

## 3. The corrected design

### 3.1 Count space (kills B1; no behaviour change)

Never form a per-channel `Var(log ρ_c)`. Algebraically identical, continuous, and **exactly 0 at n=0**:

```
Var(ρ_m) = m_mat / esp²                        #  (m/esp)² / m  ==  m / esp²
Var(ρ_n) = n_nasc / er²  +  ρ_n² · vbp[lsrc]   #  (n/er)²  / n  ==  n / er²
```

`ρ_m = 0, Var = 0` at a region source is the analytically correct limit — a **structural absence**, not an
observed zero, so no Haldane question arises.

### 3.2 Fenton–Wilkinson instead of the delta method (fixes B3)

FW moment-matches the *sum of lognormals* to a lognormal. It is closed-form, adds **no constant**, and is
conservative under both conventions. **The decisive property is convention-independent:**

> **On a single live channel, FW reduces to today's message EXACTLY.** For one lognormal `E = e^{μ+σ²/2}`,
> `CV² = e^{σ²}−1` ⇒ `log1p(CV²) = σ²` and `μ_FW = μ`. Verified over 20k random `(μ, σ²)`:
> `max|σ²_FW − σ²| = 4.4e-16`, `max|μ_FW − μ| = 8.9e-16`.

So intron↔intron, intergenic, and **every region-source edge (≥50% of all edges)** are bit-identical no-ops
**by construction**, not approximately. That alone justifies FW over delta.

### 3.3 The per-edge message (RNA⁺; RNA⁻ mirrors)

```python
n_nasc = fbp[lsrc] * sm                             # IMPUTATION
rho_n0 = n_nasc / er
v_n    = vbp[lsrc]                                  # source RUNNING belief (log frame)
E_n    = rho_n0 * exp(v_n / 2)                      # <-- convention: see §4.3
V_n    = E_n**2 * (exp(v_n) - 1) + n_nasc / er**2   # belief spread + Poisson count

m_mat  = SPs[lsrc]                                  # MEASUREMENT: spliced mass
rho_m0 = m_mat / esp
v_m    = 1.0 / n_eff_mat + s2_transfer_mat          # §4.1 (n_eff), §4.2 (transfer)
E_m    = rho_m0 * exp(v_m / 2)
V_m    = E_m**2 * (exp(v_m) - 1)

E_a    = SPd[i] / espd;  V_a = SPd[i] / espd**2     # dst-face mature ABSORBED; 0 at SPd=0

E = E_m + E_n - E_a
V = V_m + V_n + V_a
if E <= 1.0 / erd:                                  # THE GATE — see §3.4
    pass                                            # emit NOTHING (pr = 0)
else:
    s2 = log1p(V / E**2)
    mo = (log(E) - 0.5 * s2) - log(md / erd)        # -> dst log-f_pos frame
    pr = 1.0 / s2
```

### 3.4 The `E ≤ 1/erd` gate replaces the floor — and this is load-bearing

Today `rho = … − rho_mat_dst` has **no non-negativity guard**; only the *mode* is floored
(`max(rho, 1.0/erd)`), and today's `pr` never touches `rho`, so an overshoot degrades gracefully at low
precision. **Under the design that stops being true:** at a full-absorption sink `pr → n_mat_dst`, i.e. the
floored *"your RNA is at min-observable"* message delivered at **full spliced-count precision** — firing
exactly on the capture-enriched-exon-next-to-a-busy-junction case that is the recorded flagship failure.

Neither "floor the ρ in the denominator" nor "use the raw ρ" is correct: **evaluating an expansion outside the
parameter space is not a repair.** `E ≤ floor` means *"the absorption explains everything I can see"* — the
honest message is **silence**, and `pr = 0` says exactly that. It also makes FW well-posed (no lognormal fits
a non-positive mean), and it is exactly right for cfRNA (no nascent ⇒ `E → 0` ⇒ the intron correctly receives
nothing, where today it receives a fabricated floor-mode message).

**This ships independently of the two-channel split and is defensible on its own.**

---

## 4. Open issues — RANKED

### 4.1 `n_eff` for the mature count: mass, flux, or Kish? **[USER DECISION]**

| option | what | verdict |
|---|---|---|
| **K1 — mass** (`n_eff = m_mat`) | zero plumbing; strictly **conservative** (§2.2) | **recommended now**: if the design wins under a conservative bound it wins *a fortiori*, and the agreed scope is theory, not benchmark chasing |
| K2 — integer flux | `boundary_flux_{left,right}[b,2:4]`; ~10 lines; precedent in-file (`node_geometry.py:346-347` already uses integer flux for BB power + float mass for the mode) | **anti-conservative** by the same 0.078–0.45 factor |
| K3 — exact Kish | needs a new `Σw²` accumulator channel: C++ + `_accumulator_reference.py` + golden regen | correct; later |

### 4.2 Does the mature channel get a transfer variance? **[USER DECISION — the main question]**

Dropping `vbp` from the mature channel is **correct** (spliced needs no deconvolution). Dropping
**σ²_transfer** is a *different claim*, and all three reviewers converged here independently.

σ²_transfer asks *"does the source's density equal the destination's?"* — it applies to **every** channel, and
arguably **more** to mature than to nascent. At `n_mat ≈ 500` the message asserts junction-face mature density
== exon-body mature density to **±4.5%**. Structural falsifiers:

* **Alternative splicing.** `n_mat/esp` is the abundance of isoforms using *that* junction; the body carries
  every isoform containing the exon. At Ψ=0.5 the message is low by `log 2 = −0.69` nats at `pr≈500` ⇒
  likelihood ratio `e^{−119}` ⇒ **exon RNA crushed, surplus to gDNA** — which is the *known dominant* residual
  error (EXON gDNA over-call, ~87%).
* 5′/3′ coverage bias; capture probe placement; incomplete splicing.
* **The project's own validation of this exact estimator** (`unstranded_spliced_derived_rho`) measured
  **corr 0.895 vs oracle at ss0.5** ⇒ residual ≥10% — **twice the granted error bar.**

With σ²_transfer zeroed (#2) and the local reference at `0.357`, an unbounded mature message wins **~1400:1** —
a hard pin under a ≥10%-wrong assumption. And it does not stay local: `vbp[i] = 1/(pp_loc+500) ≈ 0.002` makes
every downstream `pr → n_src`, so **one junction's spliced count unmutes the whole downstream chain** and
re-asserts the bias at full confidence at every hop.

Options: **(a)** ship at `s2_transfer_mat = 0` and accept the pin; **(b)** evaluate against the NPMLE
σ²_transfer bound it will live under; **(c)** fit a mature-specific transfer variance from adjacent-junction
disagreement (same machinery, no constant).

**Recommend (b).** Otherwise this gets A/B'd into production on a number scheduled to disappear, and mature
becomes the only message path in the system with **no model-mismatch guard** — and the highest-precision one.
*"No new constant" is the defect here, not the feature.*

### 4.3 The belief carries a MISMATCHED pair **[USER DECISION — cheap, but it must be made first]**

`_scan` carries `fbp` = the posterior **mean** of `f_pos` (`simplex_logodds`: `f_pos = sum(post * f_pos_g)`)
together with `vbp` = **`Var(log f_pos)`**. That is `E[f]` paired with `Var(log f)` — two different frames. A
lognormal needs one convention, and **the choice moves B3's severity by ~4×** (§2.3 table). Candidates:

* **mean-matched** (`μ = log(fbp) − v/2`): faithful to what `fbp` *is*. Delta only 1.25× anti-conservative.
* **median-matched** (`μ = log(fbp)`): what the *current mode* `mo = log(ρ)` implicitly assumes. Delta 4.15×.

The current code is *already* internally inconsistent here — it centres the message at `log(E[f])` while
calling the width `Var(log f)`. **This is pre-existing, not introduced by #3**, but #3 is the first thing to
depend on it quantitatively. Cleanest fix: carry `E[log f]` alongside (the solver already computes it —
`mLp`/`mLn` in `_solve_ambig_logodds`) and make the message a genuine log-frame Gaussian.

### 4.4 Predicted regression nobody asked for: **this MUTES nascent where mature dominates** *(flag, don't fix)*

At an exon→B edge `SPs = 0`, so `ρ = ρ_n − ρ_abs` where `ρ_n` is the exon's **total** unspliced RNA carrying
the full `vbp`, minus a large `ρ_abs`. Classic differencing amplification: with `r = ρ_M/ρ_N = 20`,
`vbp = 2.8` ⇒ `Var(log ρ) ≈ (1+r)²·vbp = 1235` ⇒ `pr = 8.1e-4` vs today's 0.357 — a **~440× muting**. Introns
beside well-expressed genes lose their nascent message ⇒ **intron gDNA over-call in nascent-rich libraries**.
The variance is *honest*; the behaviour change is real. **Named in advance so it is not misread as a
regression when it appears.**

### 4.5 "spliced ⇒ mature" is a property of the SIMULATOR, not the data **[decides validation strategy]**

`wgs_engine.py:594-596`: `if is_nrna: seq_bytes = self._get_premrna_bytes(t_idx)` — nascent is drawn from
**fully unspliced pre-mRNA**, so a nascent read can never carry a CIGAR-N. Measured across 3 conditions: the
boundary spliced channel is **gDNA 0.000% / nascent 0.000% / mature 100.000%** (n=684,168). `oracle.py:179`
`nas_spl` is identically zero **by construction**. And `oracle.py:132-137` asserts spliced purity **only on
`region_contained`** — there is **no assertion on the boundary spliced channel**, the very channel we propose
to trust at count precision.

On real data the premise leaks: `resolve_context.h:1379-1384` promotes to `SPLICE_IMPLICIT` whenever an
annotated intron lies wholly inside a PE mate gap — **origin-blind** — and `bam_scanner.cpp:1445-1447` folds
implicit into `spliced` → ch2/ch3. Human introns are frequently <300 bp; a 400 bp **gDNA** fragment spanning
one lands in the "gDNA-impossible" channel.

⇒ **No synthetic suite can falsify the premise underpinning the highest-precision channel in the system.** The
suite will show this design winning at any precision, with no signal of contamination. Actions: add the
boundary-spliced purity assertion to `oracle.py`; measure the `SPLICE_IMPLICIT` fraction of ch2/ch3 on
LBX0190 / MO_3005 **before** granting count precision; validate on cfRNA (nascent-free ⇒ the premise holds and
the `pr = 0` gate is a genuine win).

### 4.6 Smaller, resolvable *(no decision needed)*

* **Independence of `ρ_n` and `ρ_abs`.** `Var(ρ) = ΣVar` asserts independence; sampling is independent but
  `fbp` is informed by the sided spliced floor ⇒ positive unmodeled covariance ⇒ `Var(ρ)` overstated exactly
  where the nascent message matters. Compounds §4.4. **State the assumption explicitly.**
* **Silent cliff on tag-less BAMs.** `bam_scanner.cpp:1459-1463`: an unannotated spliced fragment, or any
  fragment in a BAM without XS/ts, is **dropped entirely** — not even deposited as unspliced. Since #3 makes
  mature the only strand-free evidence at κ=½, its silent absence is a cliff. **Emit a QC scalar: fraction of
  boundaries with `n_mat > 0`.**
* **`boundary_junction_strand` is last-write-wins.** `accumulator.cpp:211` overwrites; `:259-260` merge only
  fills zeros ⇒ **worker-order dependent** where opposite-strand junctions share a coordinate. It routes the
  entire mature floor. Today that's a coin flip at 0.357 pseudo-obs; under #3 it is a coin flip at `pr = n_mat`
  — on AMBIG loci, where the evidence matters most. **Log it; revisit if the QC scalar shows collisions.**
* **Stale docstrings that will mislead re-derivation.** `effective_length.py:143-145` and
  `node_geometry.py:240-243` claim `spliced_side_eff_length` is "exactly HALF" `boundary_side_eff_length`.
  **False since the ½ moved into `effective_length.py:132`.** Measured (Gaussian FL μ=200 σ=50): at
  R=400/1000/5000 the **ratio is 1.0000**; they diverge only at R ≲ FL support (0.531 at R=100). Keeping the
  divisors separate is still right — real at short flanks — **but not for the stated reason**, and anyone
  deriving `Var(ρ_m)` from these docstrings lands 2× off. Also `bp_solver.py:422-423` still says "the
  belief-free σ²_imp + 1/n_src (`sig_imp`)" — stale since #2.

### 4.7 A finding worth recording: not every boundary is a splice junction

Oracle, `gdna300_ss_0.50_nrna_none_capture_off`: boundaries carry **44,765** units of **unspliced MATURE**
mass against 89,128 of unspliced gDNA — i.e. **33% of boundary unspliced mass is mature**, which appears to
contradict "at a junction mature splices, so unspliced crossing = gDNA + nascent."

It does not: on boundaries **with** spliced mass (n=1015), true `f_g` is **1.000 at p10/p50/p90**. The mature
crossing lives at boundaries **without** spliced mass — exon-interior region transitions, which are not splice
junctions. `SPs[lsrc] > 0` gates exactly the real junctions. **The architectural claim is conditional on being
at a splice junction, and the docs state it unconditionally.** Worth a one-line fix.

*(And: `true f_g = 1.000` at junctions here is `nrna_none` **by construction** — do not let it validate
nascent≈0 in general.)*

---

## 5. Implementation plan

Steps 1–3 are **safe to land now** — no decision required, each independently verifiable.

1. **Docstring truth pass** — `effective_length.py:143-145`, `node_geometry.py:240-243` (the false 2×),
   `bp_solver.py:422-423` (stale σ²_imp), + §4.7's conditional. No behaviour. *Verify: read.*
2. **Count-space rewrite of today's precision** (`_scan`): `pr = n_src/(n_src·vbp+1)` ≡ `1/(vbp + 1/n_src)`.
   Pure refactor establishing the algebraic frame. *Verify: byte-exact `CalibrationResult` vs HEAD.*
3. **The `E ≤ 1/erd` gate** (§3.4) — ships independently; removes today's fabricated floor-mode message.
   *Verify: `selfsolve_diag` + cfRNA — intron messages should vanish, not floor.*
4. **Plumb the spliced integer flux** into `NodeGeometry` (~10 lines, no constant). Inert until step 6; needed
   for K3's diagnostic regardless. *Verify: unit test — `spliced_n_* == 0` on regions, `≥` mass on boundaries.*
5. **Resolve §4.3** (the convention), then **§4.1** (`n_eff`) and **§4.2** (transfer).
6. **FW two-channel combine** behind a default-off `CalibrationConfig` flag. Off ⇒ bit-identical (guaranteed by
   §3.2's single-channel exactness — *assert it*). *Verify: flag-off byte-exact; flag-on diff confined to
   boundary-adjacent exon nodes.*
7. **QC scalar** (§4.6) + the boundary-spliced purity assertion in `oracle.py` (§4.5).

---

## 6. Validation — theory only, per the agreed scope

**T1 — the NaN test, first.** Run `_scan` under `np.seterr(all='raise')` on a cached payload: today passes,
the design-as-written raises. Settles B1 empirically rather than by argument.

**T2 — `selfsolve_diag.py` should NOT change at all.** It measures the *message-free* self-solve; #3 touches
only the message layer. **If it moves, the change leaked into ψ and the count-zero-information boundary was
crossed — the cleanest single falsifier available.** What *should* move is `msg_trace`'s `pr` column: on
B→exon edges at κ=½, up from ~0.357 toward `n_eff_mat` (FW-capped — expect ~285 at w=0.99 under mean-matching,
**not** 446); on exon→B edges, *down* by up to ~440× (§4.4).

**T3 — the exactness ladder** (each byte-level): FW with one live channel ≡ today's message to 1e-12 (already
verified: 4.4e-16); intron↔intron / intergenic edges full-pipeline byte-exact vs HEAD; and assert the
at-most-one-of-`SPs`/`SPd` property directly rather than trusting the comment.

**T4 — the estimator's own honesty, by simulation.** Sweep `w ∈ [0.5, 0.9999]` × `vbp ∈ [0.1, 2.8]` × both
conventions. **Falsifier: if FW is ever anti-conservative (`fw_pr > exact_pr`) anywhere in the grid, the
combine is wrong.** Adjudicates delta-vs-FW without touching a BAM.

**T5 — the premise test, which no synthetic suite can run** (§4.5). Measure the `SPLICE_IMPLICIT` fraction of
ch2/ch3 on LBX0190 / MO_3005. If non-trivial, "spliced ⇒ mature" carries a gDNA path on real data and
`s2_transfer_mat = 0` is indefensible regardless of what the suite shows.

---

## 7. REVIEW OUTCOME *(added 2026-07-15, post-review)*

Verification for every claim below: `scripts/debug/boundary_spliced_check.py`.

### 7.1 ACCEPTED — the expectation bug. §4.3 is CLOSED (it was never a choice)

**`fbp` is the arithmetic mean** — `simplex_logodds.py:350`, `f_pos = np.sum(post * f_pos_g, axis=1)`. So
`ρ_n0 = fbp·sm/er` **is already** `E[ρ_n]`, and §3.3's `E_n = ρ_n0 · exp(v_n/2)` inflates the expected nascent
abundance by `e^{2.804/2} = 4.06×` **purely from grid uncertainty**. The reviewer is right, and right that this
**dissolves §4.3**: the physical identity of `fbp` forces the mapping. There is no convention to pick.

```python
E_n = rho_n0                                     # FIXED: rho_n0 IS the arithmetic mean
V_n = E_n**2 * (exp(v_n) - 1.0) + n_nasc/er**2
# and the log-frame location is structurally forced:  mu = log(rho_n0) - v_n/2
```

**Knock-on the reviewer did not follow through — and it corrects ME twice:**

1. **My "FW is decisively better" argument was wrong.** I justified FW by single-channel exactness. **Both
   estimators have it**: delta gives `ρ²v/ρ² = v`; FW gives `log1p(e^v−1) = v`. The ≥50%-of-edges no-op
   property does **not** discriminate them.
2. **The severity of B3 drops.** Under the correct convention (4M draws, `v_n = 2.804`):

| w = ρ_m/ρ | exact | delta | FW | delta/exact | FW/exact |
|---|---|---|---|---|---|
| 0.50 | 3.31 | 1.43 | 0.63 | 0.43× | 0.19× |
| 0.90 | 32.52 | 33.72 | 6.87 | **1.04×** | **0.21×** |
| 0.99 | 356.74 | 446.31 | 285.14 | **1.25×** | **0.80×** |
| 0.999 | 497.12 | 500.30 | 497.15 | 1.01× | 1.00× |

Delta is only **1.0–1.25×** anti-conservative — not the ~5× the median convention implied. FW is
**0.19–0.80×**, i.e. *conservative but wasteful* — it discards up to 79% of the precision. Neither is right
near w≈0.9, because `S` is genuinely **not lognormal** there.

**FW still stands, for a different reason than I gave: it is conservative, and this project's failure mode is
over-confident messages, not under-confident ones.** And the distinction is second-order anyway —
`0.357 → 285` (FW) vs `0.357 → 446` (delta) are both ~800–1250× wins. **Do not over-engineer the combine.**
*(Gauss–Hermite was tested and rejected: quadrature over the nascent channel alone is exact at w≤0.9 but
**152× anti-conservative at w=0.999**, where the mature channel's own variance dominates and the scheme misses
it. A correct GH needs both channels — not worth it for a second-order gain.)*

### 7.2 REJECTED — the differencing fix. The premise and the fix are both wrong

The reviewer proposes `m_net = max(0.0, SPs[lsrc] - SPd[i])`, arguing `E_a` is "a structural subset of the
exact same local physical pool" as `E_m`.

**The premise is factually wrong: `E_m` and `E_a` never coexist.** Measured over every forward edge:
`SPs>0 only = 343`, `SPd>0 only = 342`, **`BOTH = 0`**. `SP` is identically zero on regions
(`node_geometry.py:174-178`), so a B→exon edge has `SPs>0, SPd=0` and an exon→B edge has `SPs=0, SPd>0`. There
is no `Cov(E_m, E_a)` to omit — they are never in the same sum.

**The fix is worse: it silently deletes the absorption.** On an exon→B edge `SPs = 0`, so
`max(0, 0 − SPd) = 0`. Measured: the absorption is discarded on **342/342 edges (100%)**, dropping **114,141
units** of mature mass. The absorption exists so that only *nascent* crosses into an intron; deleting it
re-introduces mature bleed into introns — **regressing a shipped fix** (`spliced_efflen_not_2x_nascent_subtraction`,
−3.5%).

**The real pair is `(ρ_n, ρ_abs)` at exon→B edges, and their estimator errors ARE independent.** `ρ̂_n =
fbp·sm/er` draws its error from the exon's strand solve on its **unspliced** counts; `ρ̂_abs = SPd/espd` draws
its error from Poisson sampling of the junction's **spliced** counts — different fragments, different channels.
The *true values* are correlated (same transcript abundance); the *estimation errors* are not, and `Var(X−Y)`
takes the covariance of the **errors**. So `Var = ΣVar` is right.

⇒ **§4.4's ~440× muting is honest, not a covariance artifact.** If you don't know the exon's RNA fraction
(κ=½), you genuinely don't know how much nascent remains after subtracting a precisely-measured mature. The
message *should* go quiet. Suppressing that with `max(0,·)` would manufacture confidence from an omission —
the same error class as the original lumping. *(The one real correlation path — `fbp` informed by the sided
spliced floor — is §4.6, second-order, and stays flagged.)*

### 7.3 ACCEPTED — §4.1 **K1 (mass)** and §4.2 **option (b)**

Both match the document's own recommendations; the reviewer's reasoning is the reasoning in §4.1/§4.2.
Recorded as **decided**. Two amendments to their blueprint:

* `v_m = 1/max(1.0, m_net)` clamps `n_eff` at 1 — **a magic number**, and anti-conservative (a boundary at
  `m_net=0.1` would get `v_m=1` instead of 10). Under K1 use `v_m = 1/m_mat` with the **§3.1 count-space form**,
  which is `0` at `m_mat=0` by construction and needs no clamp.
* `regularizers.get_npmle_transfer_variance()` **does not exist** — σ²_transfer was zeroed in #2 and its NPMLE
  successor is unwritten. See §7.4.

### 7.4 THE REMAINING DEPENDENCY — option (b) sequences #3 behind σ²_transfer

Choosing (b) means the mature channel is only safe once bounded, and **the bound does not exist yet**. So #3
**cannot fully land** until the NPMLE σ²_transfer is designed. This is the correct order — it is exactly the
"do not A/B this against a transient zero" argument — but it must be stated: **#3's flag stays default-OFF
until σ²_transfer lands.** `s2_transfer_mat` ships as a single symbol wired to `0.0`, unreachable in
production, with the flag off.

### 7.5 Revised status of the blockers

| # | | status |
|---|---|---|
| B1 | NaN at `n_mat = 0` | **resolved** — §3.1 count-space rewrite |
| B2 | `n_mat` is mass not count | **decided** — K1 (conservative); flux plumbed but inert |
| B3 | delta anti-conservative | **resolved** — FW, and less severe than believed (§7.1) |
| B4 | no transfer term | **decided** — option (b); creates the §7.4 dependency |
| §4.3 | mismatched convention | **CLOSED** — not a choice (§7.1) |

**Steps 1–4 and 6 are now landable** (6 behind a default-OFF flag). Step 5 reduces to §7.4's dependency.
