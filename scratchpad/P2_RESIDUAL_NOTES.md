# The P-2 residual — diagnosis and the derived fix (working notes, 2026-07-27)

## 0. The stated mechanism was WRONG, and the measurement says so

`pin_derivation.md` §10 records the residual as *"with nothing to limit it, an RNA-only message now asserts
too much RNA"* — i.e. a **partial** claim, **under**-calling `f_g`. Measured (`p2r_a`, `p2r_c`, `p2r_d`), at
refit=0, on the regressing stratum:

| what §10 predicted | what is measured |
|---|---|
| partial (RNA-only) claims | **0.2 %** of the harmed mass is partial |
| `f_g` UNDER-called | **96–99 %** of the regression mass is **OVER**-calling (`f_g` went UP) |
| the message asserts too much RNA | the message asserts too much **gDNA**: `e^moG` 0.42–0.77 against an oracle of 0.008–0.043 |

The off-simplex reading is falsified too: `mo_R > 0` carries **3.3 %** of the stratum's error mass, on
boundaries, in the *over*-calling direction.

## 1. What the residual actually is

Split every node by whether the pin's budget `S` borrowed the DESTINATION's own density for a component the
message did not supply — the BP violation `pin_derivation.md` §6 ranks as the bug:

```
    contaminated := ∃c :  prec_c == 0  and  own_c > 0        (a structurally dead strand lends nothing)
```

P-2's move, split (refit=0, `p2r_d`):

| condition | Δ | CLEAN | CONTAMINATED |
|---|---|---|---|
| gdna100 ss0.50 nrna_none capOFF | **+0.0148** | **+0.0129** | +0.0019 |
| gdna300 ss0.50 nrna_none capOFF | +0.0076 | **+0.0067** | +0.0009 |
| gdna100 ss0.50 nrna_present capOFF | +0.0031 | +0.0022 | +0.0009 |
| gdna300 ss0.50 nrna_present capOFF | +0.0025 | +0.0019 | +0.0006 |
| **none ss0.50 nrna_none capOFF** (P-2's win) | **−0.1219** | +0.0000 (**0 clean nodes**) | **−0.1219** |
| gdna100 ss0.50 nrna_none capON | −0.0012 | −0.0028 | +0.0016 |
| gdna100 ss0.99 nrna_none capOFF | +0.0000 | −0.0000 | +0.0001 |

**P-2 fixed the contaminated branch and over-corrected the clean one.** Its entire `gdna_none` win is on
the contaminated branch; ~85 % of the residual is on the clean branch, and on the clean branch it is
**entirely EXONS** (net +0.0042 of the stratum; intron/boundary/intergenic ≈ 0), whose delivered level
claims over-account: median `S/M` = **1.108** on those exons.

This is exactly what `pin_derivation.md` §4/§8 already says — *"the pin is CORRECT for a full claim … the
entire defect is the partial-claim branch"* — and P-2 removed it from both.

## 2. The pin's two jobs are ONE operator, and they conflict

For a complete claim `t_c = ρ̂_c^src·r`, the pin's factor is `k = M/S = M/(r·Σ_c ρ̂_c^src E_c)`, so

```
    delivered  ρ_c^i = ρ̂_c^src · M_i / Σ_{c'} ρ̂_{c'}^src E_{c'}^i          ← r CANCELS
```

**Conservation IS r-cancellation** on a complete claim: you cannot enforce `Σ_c ρ_c E_c = M` and keep the
reframe. That is the whole tension, and it is why the two ends of the suite disagree:

* **capture OFF**: `r ≈ 1`, so cancelling it costs nothing and the budget bound is pure gain → the pin wins.
* **capture ON**: `r` IS the enrichment information (§11's lesson), so cancelling it is fatal → P-2 wins
  (clean-branch −0.0028 on the capON control).

So neither endpoint is right, and "restore the pin where it is BP-legal" trades one stratum for the other.

## 3. ⭐ THE DERIVED OPERATOR — the conditional mean, and it needs NO new constant

Write the delivered densities as `t_c = ρ̂_c^src·e^s` with a common log-scale `s`. Two estimators of `s`:

* the **reframe**: `s = log r`, with the common-mode variance the solver already computes,
  `σ_cm² = Var(log r) + 1/n_src` (M5's scale term plus the source's count) — this is P1e's own `_s2c`;
* **conservation**: `s = log r + δ`, `δ = log(M/S)`, exact if the composition transfers.

`M = Σ_c ρ_c^* E_c` is an identity, so `δ = αᵀε` with `ε_c = log(ρ_c^*/t_c)` and `α` the budget shares —
i.e. **δ is the share-weighted mean of the message's own log-errors.** Under P1e's error model
`Σ = σ_cm²·11ᵀ + diag(w_c)` plus the destination's Poisson `1/n_dst`,

```
    Var(δ) = αᵀΣα + 1/n_dst           (P1e's `_den`, already computed)
    E[ε_cm | δ] = δ · σ_cm² / (αᵀΣα + 1/n_dst)
```

so the message's common scale should be corrected by **the fraction of the observed conservation violation
that its own declared common-mode uncertainty can account for**:

```
    t_c ← t_c · exp(w·δ),      w = σ_cm² / (αᵀΣα + 1/n_dst) ∈ [0,1]
```

`w = 1` is exactly the pin (`e^δ = M/S`); `w = 0` is exactly P-2. **Both shipped behaviours are the two
endpoints of one derived shrinkage, and the weight is built from quantities P1e already computes** — no new
constant, no tuned knob, no gate on capture.

And the weight moves the right way on its own: a well-counted source across a well-counted hop has a small
`σ_cm²`, so the reframe is kept (capture ON, where counts are high); a poorly-determined scale gives
`w → 1`, so the mass frame takes over (capture OFF).

### Scope, and why it is not (★)

* Applied **only where the budget borrows nothing from the destination** (`_clean`), so no destination
  belief reaches the mode — `pin_derivation.md` §1's sharp line, satisfied by construction.
* It **does not discard the reframe**: `r` is applied first and survives with weight `1−w`. This is the
  requirement §11 leaves open (*"whatever bounds a partial RNA-only claim must do so without discarding the
  reframe"*), and it is why it is not the refuted (★) — (★) normalises by the SOURCE's mass, drops `r`
  unconditionally and applies to partial claims. This keeps `r`, uses the DESTINATION's mass, and is silent
  on partial claims.
* It touches **only the level modes**. `tlam`/`tth` are scale-free, and `_pin_v`'s output (the DL comparison
  frame) normalises to `M`, so a common factor cancels there identically. Nothing else in the packet moves.
* It is the **MODE fix for the bias half of `δ`** that `variance_ledger.md` §6 and `ROADMAP` §P1e say must
  exist — *"the right fix for the bias half is a MODE fix, not a variance"*.

## 4. Pre-registered predictions (before the A/B)

1. ⭐ **`gdna_none` does not move** — 0 clean nodes on the diagnosed condition. **FALSIFIER**: if
   `gdna_none` regresses, the branch split is wrong and the whole diagnosis fails.
2. The regressing stratum (unstranded × capOFF × gDNA-bearing) improves.
3. Stranded capture-OFF stays flat (its clean-branch delta is −0.0000).
4. capture-ON pays **less** than the crude restore (`w = 1`) would, because `w < 1` where counts are high.
5. Suite mwae improves at refit=0 and refit=1.

**Control arm (`w ≡ 1` on the clean branch)** = "restore the pin exactly where it is BP-legal". Run it too:
if it matches the derived `w`, the derivation adds nothing and the crude form ships (the lesson P1e's own
record insists on).
