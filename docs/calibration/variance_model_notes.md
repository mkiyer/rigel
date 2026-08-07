# THE TRANSFER DERIVATION — what crosses an EDGE, and what it costs

    Working design doc. Opened and rewritten 2026-08-06 after the first draft's §4 was found WRONG and
    its §5 was refuted by measurement. ⚠ **NOT a permanent fixture.** The six permanent docs are
    `CLAUDE.md`, `docs/SUCCESS.md`, `docs/ROADMAP.md`, `docs/TRAPS.md`, `docs/EQUATIONS.md`,
    `docs/DESIGN.md`, `docs/TESTING.md`.
    ⛔ **When this lands**: the derivation → `EQUATIONS.md`, the rulings → `DESIGN.md`, the lessons →
    `TRAPS.md`, the numbers → `ROADMAP.md`, and this file is DELETED.

⛔⛔ **READ §2 AND §6 BEFORE PROPOSING ANYTHING.** **Eleven** candidate mechanisms have been priced on the
full panel and all eleven are refuted — including **this document's own §2** and the **η rebuild itself**.
§6 says what each was and why it died. The body contains only what survived scrutiny; a previous draft
advocated a mechanism §6 now buries, and that self-contradiction is the worst failure mode for a design doc.

⛔ **THE CURRENT STATE OF THE WORK IS `SESSION_HANDOFF.md`, NOT §13 BELOW.** This file is the DERIVATION;
the handoff is what was measured against it. §13 has been pruned to avoid two homes for one account.

---

## §0 THE AXIOMS

**A1. THREE POPULATIONS: `gDNA`, `RNA+`, `RNA−`. THERE IS NO FOURTH.** "Mature" and "nascent" are not
populations, not species, not a degree of freedom. RNA in an intron is RNA that has not spliced *at that
position*. The only distinction is **spliced** (certified RNA, gDNA cannot splice, no deconvolution needed)
versus **unspliced** (the whole problem). ⛔ **The tell that this is being violated: a population set with
more than three members.** The first draft of this very document opened with five and every table built on
it came out wrong. `CLAUDE.md` **AXIOM 0** is the canonical statement.

⚠ **Three populations is not three composition components.** The composition is **two** — gDNA vs
RNA-total — one degree of freedom. The strand split is a **nuisance tilt**. Keep the three words apart.

**A2. PASS-0 IS PRIOR-FREE, AND IT IS A TRAINING SET, NOT A DELIVERABLE.** The gDNA hyperprior is fitted
*on pass-0's output*, so pass-0 may use a slot's own counts, its own lengths, and incoming messages — and
**no population statistic** (that is the circularity; owner's ruling 2026-08-06).

**A3. D4.** A message may use the destination's CONSTANTS and its OBSERVATIONS, never its BELIEFS.

**A4. NO TUNED CONSTANTS.** Every number below is a Jeffreys ½, a digamma/trigamma of it, or a geometric
quantity read off the index.

---

## §1 THE POPULATION SET

For slot `i`, the populations that can deposit into its **unspliced** bank:

$$T_i \;=\; \{\mathrm{g}\ \text{if}\ E_g^i>0\} \;\cup\; \{\mathrm{R^+}\ \text{if}\ \texttt{free\_pos}_i \wedge E_r^i>0\} \;\cup\; \{\mathrm{R^-}\ \text{if}\ \texttt{free\_neg}_i \wedge E_r^i>0\}$$

- gDNA is genomically continuous, so it is admissible wherever it has **opportunity** — there is no gDNA
  analogue of a forbidden strand.
- Each RNA strand is admissible iff the annotation frees it there **and** there is RNA opportunity.

⭐ `|T_i| ∈ {0,1,2,3}`, and these are exactly HEAD's own `live` predicates (`eff_global > 0`;
`free_s & (eff_r > 0)`), so the predicate is a **read**, not a new model. ⚠ The opportunity gates matter:
a NODE shorter than a fragment has `E = 0` and holds no population at all.

---

## §2 ⛔⛔ THE MEAN-LOCATION FORM IS REFUTED, AND IT IS NOT PART OF THE REBUILD

⛔ **Measured 2026-08-06 inside the `eta` prototype: +96,299 % on the `g00` zero control.** Two reasons it
comes out, and the second is the one that matters:

1. **It cannot say ZERO.** `exp(ψ₀(a+½))/E > 0` at every count, so a structurally pure-gDNA slot over 50 Mb
   of empty genome asserts a small POSITIVE gDNA density — and at `g00` the truth is exactly 0, so every
   one is a false positive with nothing to cancel it. `node_init.rho_g` is an EXACT zero at **60,544 of
   70,176** ladder slots, and that is the statement earning the −98 % at `g00`. §2 replaced all of them.
2. ⭐⭐ **It was credited with a benefit that belongs to §4.** The claim below — "C0d dissolves rather than
   being guarded against" — was **wrong**. C0d mattered because a *multiplicative* transport cannot carry
   zero. §4/§5 remove the multiplicative transport altogether: the conversion is an **additive constant in
   log space** and the density crossing is the identity, so nothing ever divides by a density and zero
   transports perfectly. §2 was never needed for it.

⚠ **And it is two things varied.** The rebuild is a STRUCTURAL simplification; §2 is an ESTIMATOR change
riding along inside it. It is now a separate later experiment with its own arm.

### ⭐ If a length-aware location is revisited, the derived answer is the MODE, not the mean

The same posterior already supplies the one-sided claim a zero count honestly is. For `Gamma(α, β)` the
mode is `(α−1)/β` for `α ≥ 1` and `0` otherwise, so with `α = a+½`, `β = E`:

$$\mathrm{mode}(\rho) = \frac{\max(a - \tfrac12,\,0)}{E}$$

- **exactly 0 at a zero count** — the statement that earns the −98 % at `g00`;
- same posterior, same Jeffreys ½, **no new constant**;
- the length-awareness moves into the PRECISION, where it belongs: linear variance `(a+½)/E²`, hence
  precision `2E²` at `a = 0` — **220,000×** apart between a 38 kb gap and an 82 b scrap;
- and it matches the interface, whose variable is literally named `mo_g` — a **mode**.

⚠ **Caveat to keep visible:** the mode has a kink at `a = ½` and `Gamma(½, E)` is badly non-Gaussian, so a
`(mode, precision)` summary is a poor approximation there. Worth a Monte-Carlo check before trusting it —
not a blocker.

### The original derivation, kept because the MODE form is built on it

⛔ Everything below is retained for the posterior it derives. **The mean-location conclusion it reaches is
refuted above; do not re-implement it.**

## §2a THE RATE POSTERIOR — and length lives in the LOCATION

`a` fragments (integer) over opportunity `E` (integer). Poisson likelihood, Jeffreys prior for a rate
(`ρ^{-1/2}` — the same `_JEFFREYS_REF = ½` ψ is built on). Exactly:

$$\rho \sim \mathrm{Gamma}(a+\tfrac12,\,E)
\quad\Longrightarrow\quad
\mathbb{E}[\log\rho]=\psi_0(a+\tfrac12)-\log E,
\qquad \mathrm{Var}(\log\rho)=\psi_1(a+\tfrac12)$$

⛔ **HEAD pairs this posterior's VARIANCE with a different posterior's LOCATION.** It uses `log(a/E)`,
which is not `E[log ρ]`; at `a = 0` that is a finite variance attached to a location of `−∞`. The repair is
one term, and it is where length enters:

$$\log\rho=\underbrace{\psi_0(a+\tfrac12)}_{\text{counts}}-\underbrace{\log E}_{\text{length}}$$

`ψ₀(½) = −γ − 2log2 = −1.9635`, so a zero count says `ρ = 0.1404/E`:

| empty stretch | claims | log-variance |
|---|---|---|
| 38,400 b — a typical inter-gene gap, gDNA-free library | `3.65e−6` | 4.9348 |
| 82 b — a scrap, 75 % gDNA under capture | `1.72e−3` | 4.9348 |

**470× apart in the claim, identical in relative confidence** — which is right: a zero count is always "a
factor of ~9 either way", and the length sets what it is a factor of nine *around*.

⛔⛔ **AND "no exact zero is ever formed, so `TRAPS.md` C0d dissolves" WAS THE ERROR.** That property is a
COST, not a benefit — it is exactly what removes the ability to say "there is none here". C0d is dissolved
by **§4**, which deletes the multiplicative transport that made zero unrepresentable. See §2 above.

⛔ **MEASURED TWICE, AND BOTH TIMES REFUSED:** alone this is the `zc_logmean` arm (**+6,264 %** on the
zero-gDNA control); wired structurally as the rebuild's level channel it is **+96,299 %** (§6).

---

## §3 WHEN COMPOSITION CAN CROSS

A composition claim may cross `s → d` only if the two slots draw on the same populations:

$$\boxed{T_s = T_d}$$

Otherwise only `T_s \cap T_d` may cross. Three regimes, three different problems:

| | | |
|---|---|---|
| `T_s = T_d` | same populations | §4 — the invariant crosses exactly |
| `T_s \supsetneq T_d` | source knows strands the destination cannot hold | **drop** them (the existing peel) |
| `T_s \subsetneq T_d`, or incomparable | destination holds strands the message is silent about | §5 |

**The table (A1 respected):**

| pair | `T_s` | `T_d` | |
|---|---|---|---|
| intergenic NODE ↔ intergenic\|exon EDGE | {g} | {g} | ✅ §4 (degenerate: `\|T\|=1`, no composition to carry) |
| intron NODE ↔ intron\|exon EDGE | {g, R_s} | {g, R_s} | ✅ §4 — RNA in an intron is RNA |
| intron\|exon EDGE ↔ exon NODE | {g, R_s} | {g, R_s} | ✅ §4 |
| exon\|exon EDGE ↔ exon NODE, same transcripts | {g, R_s} | {g, R_s} | ✅ §4 |
| **intergenic\|exon EDGE ↔ exon NODE** | **{g}** | **{g, R⁺/R⁻}** | ⛔ §5 |
| **across a TSS/TES bit** | the strand set changes | ⛔ §5 |
| **+strand exon ↔ −strand exon** | {g,R⁺} vs {g,R⁻} | ⛔ §5, intersection {g} |

⚠ **The "two faces" of an `intron|exon` EDGE is a DIFFERENT question.** Both faces have the same `T`; they
differ in which junction flux belongs on which side. §9 shows that question disappears.

---

## §4 ⭐⭐⭐ THE FRAME-FREE COORDINATE — the corrected centrepiece

⛔ **The first draft of this document claimed "composition transfer is the identity on λ". THAT IS FALSE**,
and the error matters because it is the difference between a belief-free exact operation and a wrong one.

**What actually transfers is the DENSITY, not the share.** "gDNA is uniform" and "the same transcripts run
through both slots" are statements about *densities*. Shares are derived, and they move whenever the two
slots' opportunity ratio differs — which between a NODE and an EDGE it always does. Verified:

| `E_g,E_r` at source → destination | `φ_g` | `λ` | `η` |
|---|---|---|---|
| (100, 200) → (250, 100) | 0.3333 → 0.7143 | −0.6931 → +0.9163 | **0 → 0** |
| (3000, 1500) → (254, 254) | 0.6667 → 0.5000 | +0.6931 → +0.0000 | **0 → 0** |
| (50, 400) → (900, 300) | 0.1111 → 0.7500 | −2.0794 → +1.0986 | **0 → 0** |

Since `λ = log(φ_g/φ_R) = log(ρ_g/ρ_R) + log(E_g/E_r)`, define

$$\boxed{\;\eta \;\equiv\; \lambda - \log\frac{E_g}{E_r} \;=\; \log\frac{\rho_g}{\rho_R}\;}$$

**`η` is the frame-free composition coordinate, and it transfers as the IDENTITY across any edge with
`T_s = T_d`.** The tilt `θ` — built from `(ρ_p − ρ_n)/ρ_R`, in which the common `E_r` cancels — likewise
transfers as the identity.

**Four consequences, and together they are the whole simplification:**

1. **The conversion is a known geometric constant.** `λ_d = λ_s + [log(E_g/E_r)_d − log(E_g/E_r)_s]`. No
   belief, no count, no unknown — read off the index.
2. **The mass identity holds by construction.** Rebuild at the destination as
   `ρ_c(d) = φ_c(d)·M_d/E_c(d)` with `φ` from `λ_d`; then `Σ_c ρ_c E_c = M_d Σφ_c = M_d`, exactly, for any
   `η` (verified to 1e-10). ⭐⭐ **So there is nothing for the mass pin to restore.**
3. **No ratio of totals is ever formed** — so no reframe `r`, no `framed` guard, no undefined ratio, no
   `rho > _EPS` guards.
3b. ⭐⭐ **AND THIS — NOT §2 — IS WHAT DISSOLVES `TRAPS.md` C0d.** C0d says a MULTIPLICATIVE transport
   cannot carry zero, so "there is none here" is unrepresentable by construction. Under §4/§5 the
   conversion is an **additive constant in log space** and the shared density crosses as the **identity**:
   `rho_g = 0` needs no `rho > _EPS` guard, is never a denominator, and arrives as `f_g = 0` intact. ⛔ The
   first draft credited this to §2's Jeffreys location, which is the opposite of true — §2's positive floor
   is what *destroys* the zero, and it cost +96,299 % on the control to find out (§6).
4. ⭐ **It explains what the reframe was for.** `r = ρ_tot(d)/ρ_tot(s)` is the code compensating, with
   *beliefs*, for working in `λ` where it should have worked in `η`. The correct compensation is exact and
   belief-free.

⚠ **The premise `η` inherits, stated plainly:** `η` is invariant only where gDNA density and RNA density are
each unchanged — i.e. within one capture stratum, with the same transcripts. `T_s = T_d` is a proxy for the
second condition; the first is the same uniform-gDNA assumption the whole tool rests on. This is not a new
weakness, but §4 does not remove it.

⚠ The **spliced** bank is untouched: a spliced fragment is certified RNA, entering as a *measurement* of the
RNA side. The identity above is over the **unspliced** bank, which is what needs deconvolving.

---

## §5 THE MISMATCHED EDGE — under the null there is NO UNKNOWN

`T_s ⊊ T_d`: the source (say an `intergenic|exon` seam, `T_s = {g}`) knows only `ρ_g`.

⭐ **A density is already frame-free.** So under the null — gDNA uniform, no capture enrichment —
`ρ_g` crosses **unchanged**, and the destination's own observations convert it:

$$f_g(d) \;=\; \frac{\rho_g(s)\cdot E_g(d)}{M_d}, \qquad \text{RNA mass} = M_d - \rho_g(s) E_g(d)$$

**Fully determined. No interval, no free parameter, no scale to integrate out.** The newly-active
population is not "unconstrained by the message" — the null *implies* it as the residual of the
destination's own count. This is the owner's formulation and it is a strict improvement on the first draft's
interval framing.

⭐⭐ **And this is exactly what HEAD already computes** — `mo_g = log(c_g·E_g/M) = log f_g`.

### §5.1 The general case — TSS/TES and strand changes, in one rule

Let `C = T_s \cap T_d` (shared) and `N = T_d \setminus T_s` (newly active). Under the null:

1. **Within `C`, the density RATIOS are frame-free and transfer as identities** (§4) — so the relative
   composition *among the shared populations* crosses exactly.
2. **Each shared density crosses unchanged**, so the destination-frame mass they account for is determined:
   $$A \;=\; \sum_{c \in C} \rho_c(s)\,E_c(d)$$
3. **`N` takes the residual** `M_d − A`, and how it splits depends only on `|N|`:

| `\|N\|` | example | the split |
|---|---|---|
| 0 | `T_s = T_d` | §4 — nothing new; the message is complete |
| 1 | a `+`-strand TSS: `{g} → {g, R⁺}`; or `{g,R⁺} → {g,R⁺,R⁻}` | **determined** — one unknown, one equation |
| 2 | `{g} → {g, R⁺, R⁻}` (an intergenic seam into a both-strand exon) | the RNA **total** is determined; the **tilt** `θ` is not, and falls to ψ's reference |

⭐⭐⭐ **So the only thing a component mismatch ever leaves undetermined is the TILT, and only when both
strands are newly active.** The gDNA-vs-RNA split — the quantity the tool exists to estimate — is
determined in every case. That is a much smaller residual unknown than the first draft supposed, and it is
`θ`, a nuisance parameter, not `λ`.

⚠ `A > M_d` is possible from sampling noise (the shared populations "account for" more than the destination
observed). Then the null is inconsistent with the destination's own count and the honest reading is `f_N = 0`
with the excess absorbed by the shared claim's own uncertainty. ⛔ Do **not** introduce a shift here; §6.

⛔⛔ **So the open problem is not "what crosses". It is: HOW MUCH SHOULD THE DESTINATION DOUBT THE NULL?**
Capture makes `ρ_g(d) = γ·ρ_g(s)` with `γ ≥ 1` possible; nothing at pass-0 is evidence that it *happened*.
§6 shows that every rule that introduces doubt is refused by the zero control, and §8 says where the doubt
must come from instead.

---

## §6 ⛔⛔⛔ THE GRAVEYARD — eleven candidates, all refuted, with the reason each died

⚠ Two scoring conventions, and mixing them is how "−7.7 %" gets read as "most of it": the **5-condition
localisation set** (`g75/g98 st ON`, `g75 un ON`, `g75 st OFF`, `g00 st ON`) and the **full 36**. All values
`Σ|err|` node+edge unless marked. `g00` is the owner-required zero-gDNA control: truth is 0, so every
fragment is a false positive with nothing to cancel it.

| candidate | what it did | `g00` | target | verdict |
|---|---|---|---|---|
| `zc_jeffreys_mean` | `ρ_g = ½/E_g` at zero mass | ⛔ +7,269 % | −13.9 % | mode → up |
| `zc_logmean` | `ρ_g = e^{ψ₀(½)}/E_g` (**= §2 alone**) | ⛔ +6,264 % | −11.3 % | mode → up |
| `zc_anchor_mute` | no `prec_g` at empty locked slots | ⛔ +5,554 % | −7.7 % | kills the win |
| `zc_struct_lock_g1` | scope `struct_lock` to `g1_locked ∧ NODE` | ⛔ +3,207 % | −1.2 % | the mis-scoped mask is load-bearing |
| `zc_reference_var` | `Var(f_g)=⅛` where `τ=0` | ✅ **+0.0 %** | −0.3 % | ⭐ passes the control, **inert** |
| `zc_discrepancy` (36) | `+½log D` shift, `(log D)²/12` | ⛔ +982 % | panel **+4.5 %** | mode → up |
| `zc_disc_var` (36) | variance only, mode untouched | ⛔ +255 % | panel **+0.9 %** | damping can't bite |
| `zc_ref_prior` (36) | own belief = ψ's reference, `τ+1/π²` | ⛔ +3,792 % | stON **−14.9 %** | mode → up |
| `zc_ref_prior_damp` (36) | the two above, paired | ⛔ +3,809 % | stON **−15.5 %** | ditto |

Plus four **decomposition reverts** used to attribute the regression, not proposed as fixes:
`zc_own_count` (+5,554 % / −1.2 %), `zc_live_count` (+5,554 % / −1.1 %), `zc_total_n` (inert),
`zc_transfer` (+5,820 % / −28.8 % — **this one reproduces the pre-fix tree**).

### ⛔⛔ §2 WIRED AS A STRUCTURAL FLOOR — the graveyard at 15× the magnitude (2026-08-06, `eta` arm)

The `eta` prototype put §2's mean location in the LEVEL channel, so every structurally pure-gDNA slot
originated a positive gDNA density instead of an exact zero. `zc_logmean` applied the same lever only at
zero-mass slots and cost +6,264 %; applied structurally it costs **15× more**. `abs_err_all`, node axis:

| condition | base | `eta` (with §2) | |
|---|---|---|---|
| ⛔ `g00 ss0.99 capture_on` | 913 | 885,962 | **+96,299 %** |
| ⛔ `g00 ss0.50 capture_off` | 36,651 | 1,696,332 | **+4,528 %** |

`library_f_gdna` at `g00 ss0.99 capture_on`: truth **0.0000**, base **0.0002**, `eta` **0.1514**.

⭐ **And removing §2 did NOT fix the control** — the same rows are unchanged with HEAD's location in place.
So §2 is refuted on its own merits (§2) but it is *not* this control's mechanism. Three hypotheses have
been eliminated by measurement, and the mechanism is still open:

| hypothesis | test | verdict |
|---|---|---|
| the LEVEL channel carries the off-probe floor everywhere | `--arm eta_nolevel` mutes it | ⛔ **inert** — +96,839 % vs +96,299 % |
| §2's positive floor removes the ability to say zero | HEAD's location restored | ⛔ **inert on `g00`** — byte-identical rows |
| a G1 **EDGE** is a phantom-gDNA emitter (`TRAPS.md` D4j) | emission scoped to NODE slots | ⛔ **inert** — every emitter was already a NODE |

⚠ ⭐ **The own-belief census refutes the "many weak claims fuse into one confident one" candidate too.**
At `g00 ss0.99 capture_on`, mass-weighted own `f_g` BEFORE any message is **0.0171 in BOTH arms** (the
self-solve is shared), and **0 of 70,176** slots have any own composition evidence (`tau = 0` — a zero-gDNA
library gives the strand channel nothing to say about gDNA). So own-belief formation is identical and the
entire difference is in the MESSAGE layer, which is where the next dissection must start.

⭐ **One real defect the census DID pin:** the gDNA measurement channel delivered a *share* of `f_g = 5.4e9`
(median `log f_g` = **+22.41** across 6,128 receiving slots, carried with real precision), because it
applied a neighbour's density to the destination's geometry without §5's own consistency cap. A share is
in [0,1] by definition. Fixed by capping at the destination's own mass — and it is **not** sufficient.

### ⛔⛔ THE η REBUILD ITSELF, ON THE FULL PANEL — candidates 10 and 11 (2026-08-06)

The prototype (`scripts/design/eta_node_sweep.py`, `--arm eta`) implements §1/§3/§4/§5 with §2 REMOVED and
HEAD's location. `abs_err_all_final`, 32 conditions, `g00` excluded because its +136,787 % swamps 35 others:

| stratum (node axis) | base | `eta` | |
|---|---|---|---|
| **ALL** | 11,699,312 | 23,730,212 | ⛔ **+103 %**, better on **6/32** |
| unstranded × capture **OFF** | 883,530 | 7,384,192 | ⛔ **+736 %** |
| stranded × capture OFF | 701,042 | 2,476,993 | ⛔ +253 % |
| stranded × capture ON | 1,408,188 | 2,631,181 | +87 % |
| unstranded × capture ON | 8,706,554 | 11,237,846 | +29 % |

⭐ **The one real win is `g75 ss0.50 capture_ON` at −51 % on the deliverable** (`library_f_gdna`
0.4981 → **0.7058** against a truth of 0.7942) — the largest single-condition improvement this campaign,
on the stratum §8's read flags as the one where the fitted prior is DEGRADED. ⚠ But it is ONE condition
inside the least-bad stratum, and capture-**OFF** is where the rebuild collapses. **`TRAPS.md` B18 for the
fourth time: no capture-OFF condition except `g00` had been tested before the panel ran.**

⛔ **So the composition transfer is not the problem and is not yet the answer either.** Its sign is gated
(0/20,000 mismatches on asymmetric fixtures), the mass identity holds, and `η` crosses exactly. What is
missing is whatever HEAD does on capture-OFF that the rebuild does not — and the `g00` mechanism (§6a) is
still open after five eliminated candidates.

### ⛔ THE `g00` CANDIDATES — five eliminated, mechanism OPEN

See `SESSION_HANDOFF.md` §6 for the full table and the census numbers. In short: the level channel (inert),
§2's floor (inert on `g00`), a G1-EDGE phantom emitter (inert — every emitter was already a NODE), "many
weak own-claims fusing" (refuted — own `f_g` is 0.0171 in BOTH arms and 0/70,176 slots have own composition
evidence), and a sign inversion (refuted — 0/20,000 mismatches). ⭐ The `+22.41` delivered share was the
DESTINATION-MASS epsilon (100 % of receivers have `M = 0` exactly), and ⛔ that channel cannot have caused
`g00` because those slots carry no mass — A14 a third time.

### ⛔⛔⛔ The pattern, and it is the campaign's real result

| candidate | how it moved `f_g` | `g00` |
|---|---|---|
| geometric midpoint | up by `√D` | +982 % |
| the ceiling, via the mass pin | up to 1 | +6,264 % |
| the reference own belief | up toward ½ | +2,118 % (final) |

Three independent derivations, three different mathematical routes, **one identical footprint**. They are
the same lever in different clothes: **"when in doubt, more gDNA."** Every one helps high-gDNA
stranded × capture-ON; every one is refused by the zero control, because at `g00` the doubt must resolve to
*no* gDNA.

$$\boxed{\text{No prior-free, library-agnostic rule separates "the discrepancy is RNA" from}}$$
$$\boxed{\text{"the discrepancy is capture-enriched gDNA" — they are observationally identical.}}$$

⛔ **THE PASS-0 ROUTE IS CLOSED. Do not spend another arm on it.**

### Two structural facts that explain *why* nothing worked

1. **Damping is inert where it is needed.** `_fuse(own, 0, msg, p) = msg` for **any** `p > 0`, and the
   dissection found **100 % of the top error on `relay_only` objects** — own precision exactly zero. So on
   precisely the objects carrying the regression, no amount of damping can have any effect. ⚠ *A reasoning
   error to avoid repeating:* having derived this, the first draft then predicted damping would fix the
   regression because "82 % of stranded mass has own evidence". Those are different populations — the 82 %
   is where the **mass** is, the `relay_only` remainder is where the **error** is. **Check which population
   carries the error before predicting a lever reaches it.**
2. **Giving the destination an own belief fixes the target and destroys the rest.** `zc_ref_prior` cuts the
   regression 15–24 % at *every* rung — more than anything else — because the reference `f_g = ½` is a pull
   toward more gDNA. At `g00` and on unstranded data the messages were already right, so the same pull is
   catastrophic (+2,118 % final, +28 % unstranded).

---

## §7 THE TWO TIMEPOINTS

**Propagation.** Every slot processes and relays, belief or no belief.

| | |
|---|---|
| `T_s = T_d` | `η ← η_s`, `θ ← θ_s`; convert to the slot's own `λ` by adding its known `log(E_g/E_r)`. No ratio, no shift |
| `T_s ⊊ T_d` | `ρ_g` crosses unchanged (§5); the destination's own `M_d`, `E_g` give `f_g`; the residual is RNA |
| `T_s ⊋ T_d` | drop `T_s \setminus T_d` first, then as above |

⚠ **Do not add a discrepancy shift or a discrepancy damping here.** Both were measured (§6) and both are
refused by the zero control.

**Solve.** Claims arrive on one axis with variances and ψ combines them with its reference. A `T_s = T_d`
message that carries certified RNA (a junction's spliced count) constrains the split directly; a `T_s ⊊ T_d`
message constrains `f_g` under the null and is silent on the tilt.

---

## §8 ⛔ ANSWERED — THE PRIOR'S CAPACITY IS NOT THE PROBLEM, AND NEITHER IS ITS SIGNAL

⛔⛔ **BOTH BRANCHES OF THE QUESTION BELOW ARE CLOSED BY MEASUREMENT (2026-08-06).** The fitted prior
already renders the bimodal landscape — **2.98 decades** of mode separation at `g75 ss0.99 capture_ON`,
30× more enriched mass ON than OFF, a single pile at the resolution wall at `g00`. And a prior fitted from
ORACLE truth is the SAME prior on the regressing strata (**0.04 dec** apart), so better pass-0 cannot help
it either. ⭐ The successor question — *why is prior fidelity ANTI-CORRELATED with deliverable quality?* —
is `ROADMAP.md` §4.1. `SESSION_HANDOFF.md` §2/§3.

### The original argument, kept for the reasoning it records

## §8a WHERE THE DOUBT WAS THOUGHT TO LIVE — pass-1, and the prior's CAPACITY

§6 closes pass-0. The quantity that separates the two worlds is **the library's gDNA density landscape**,
and that is exactly what the refit fits.

$$\text{pass-0: trust the null } (\gamma = 1)
\qquad\longrightarrow\qquad
\text{pass-1: } \pi(\gamma) \text{ from the fitted gDNA density prior}$$

⚠ **But this does not close by itself, and the measurement says so.** On the deliverable the regression
**grows** through the refit, +14.9 % → +30.0 % (`TRAPS.md` A15), while the win attenuates −37.2 % → −3.9 %. A prior
fitted on a pass-0 that is capture-blind cannot supply capture-awareness.

⭐⭐ **So the sharp open question, and it is the best next experiment after the rebuild:** the fitted prior
is described as an NPMLE-ish mixture with a smooth low-density background component. **Can it represent a
BIMODAL landscape — an off-target floor and an on-target enriched mode?** If yes, the doubt is
representable and the work is to give it enough signal. If it is structurally unimodal, no amount of pass-0
quality will help and the prior itself is the build. ⛔ **Answer this before designing another mechanism** —
it is a property of existing code, so it is a read, not an experiment.

---

## §9 WHAT GETS DELETED — the payoff

| deleted | because |
|---|---|
| the **mass pin**, both licence cases (`_lend or g1_locked`), and `_pin_v`'s budget | §4 consequence **2** — the identity holds by construction. ⭐ This operator has a measured failure in **both** directions: it caused this campaign's regression (a precise zero read as "supplied" ⇒ all mass to RNA) and the historical 29.3 % phantom |
| the reframe `r`, `framed`, and the flank pair `rho_lo`/`rho_hi` with its role-pairing in both twins | §4 consequences **3–4** — nothing forms a ratio of totals; the conversion is a known constant |
| four `rho > _EPS` guards | §2 — `ψ₀(a+½) − log E` is finite at every `a ≥ 0`, so no exact zero exists |
| `pop`, `lend`, graft/peel **as separate gating concepts** | §3 — one `T` comparison and one set intersection |
| `transfer_logvar`'s non-graft branch and `graft_frame_logvar` | §4 — a belief-free exact conversion has no scale noise to charge |
| ~~`struct_lock` as a flag distinct from `g1_locked`~~ | ⛔ **WITHDRAWN.** `\|T\| = 1` is structural certainty at a **NODE** and NOT at an **EDGE**: an EDGE's `T` is the INTERSECTION of its endpoints', so `T = {g}` says only "nothing but gDNA crosses here contiguously" — never that the mass on the seam IS gDNA. The two predicates differ for a real reason and the derivation REPRODUCES it rather than deleting it. ⚠ The mis-scoped mask is also load-bearing today (§6) |

⭐ A dozen hand-written special cases become **three reads and one conversion**: the population set, the
rate posterior, the opportunity ratio, and `η`.

---

## §10 THE BUILD — a standalone sweep, with a NEUTRALITY target

⭐⭐ **The acceptance criterion is accuracy-neutrality on the DELIVERABLE with the operators deleted — not
an accuracy win.** §6 forbids expecting this to fix the regression, and making that its bar would repeat
the campaign's central mistake. Evidence it is the right bar: on the deliverable the pre-fix tree and HEAD
are within 3.9 % of each other while every candidate is 16–20 % worse. **The shipped solve is remarkably
insensitive to pass-0 detail.**

**Why standalone rather than an in-place refactor.** Under §4 the new sweep is *small* — no pin, no reframe,
no flank pair, no graft/peel gating, no `framed`/`lend`/`pop`: a forward and a backward pass over `(η, θ)`.
On the order of 150 lines against ~1,200. Writing it whole in `scripts/design/` and injecting it through the
arm harness gives a clean measurement of the framework **as a whole** instead of a sequence of partial edits
each fighting the operators that remain — which is the obstruction that made this campaign hard to read.

⚠ **This is not a CLAUDE.md violation.** The no-legacy rule forbids keeping the old version *in `src/`* for
comparison. A prototype in `scripts/` is the pattern this project already uses; when it wins it moves in and
the old core is deleted in one commit.

**Order of work — ⛔ ALL FOUR ATTEMPTED; the outcome of each is recorded, not predicted:**

| | | outcome |
|---|---|---|
| 1 | §8's read — can the fitted prior represent a bimodal landscape? | ✅ **DONE. YES**, and both branches of the question closed (§8) |
| 2 | Write the standalone `(η, θ)` sweep implementing **§1, §3, §4, §5** | ✅ **DONE** — `scripts/design/eta_node_sweep.py`, ~230 lines, 24 gates, 11 firing perturbations. ⛔ **§2 is NOT in it** (§2, §4b of the handoff) |
| 3 | A/B it whole against HEAD on all 36 | ⛔ **DONE — NEGATIVE.** +103 % on the deliverable, 6/32 better, worst on capture-OFF (§6) |
| 4 | Move into `src/`, delete the old core | ⛔ **BLOCKED** — it has not earned it. The four STRICT xfails remain xfail |

⛔⛔ **THE OPEN BLOCKER IS THE `g00` MECHANISM**, +96,923 % with FIVE candidates eliminated by measurement.
`SESSION_HANDOFF.md` §6 carries the elimination table; the next step is a `worst_objects.py` dissection,
not a sixth hypothesis.

## §11 GATES

⭐ **DONE — these are now 24 real gates in `tests/calibration/test_eta_transfer.py`**, over the kernel in
`tests/calibration/_eta_reference.py` (ONE home, shared with the `scripts/design/eta_node_sweep.py`
prototype). `TRAPS.md` A2's second half was run: **11 perturbations, each a specific wrong derivation, all
11 fire**, control green.

⚠ **One of them fired NOTHING on the first pass, and that is A14 in miniature:** the perturbation scaled by
`mass_dst / mass_dst.mean()` and every transfer fixture was a **single-element array**, so the factor was
exactly 1.0 and the arm could not have fired at all. Closed with a whole-chain gate, because the sweep
calls the transfer vectorised over every slot.

**The original list, for the record:**

| | |
|---|---|
| `η` invariance | exact to 1e-12 across three opportunity pairs (§4) |
| the mass identity rebuilt from `η` | exact to 1e-10 |
| `Var(λ)` under `Beta(½,½)` | `2ψ₁(½) = π² = 9.86960440108936`, to 15 digits |
| `ψ₀(½)`, `ψ₁(½)` | `−1.9635`, `π²/2 = 4.9348` |

**Gates the build needs:**

- ⭐ **The η-invariance gate**: on a fixture where the two slots' `E_g/E_r` differ, the delivered `λ` must
  differ from the source's by exactly the geometric constant — and the shares must **not** be equal. A gate
  asserting share equality would enshrine the first draft's error.
- ⭐ **The mirror gate** (A13's shape — an invariance the code cannot fake): on a palindromic fixture a
  gDNA-only message and an RNA-only message must deliver mirror-image answers.
- **`Σ_c ρ_c E_c = M`** exactly, at every slot, on every arm.
- **Zero controls** on every arm: both `g00` captures and the `g98` end.
- **`|T| ≤ 3` asserted** — AXIOM 0 made executable.

---

## §12 OPEN QUESTIONS, ranked

1. ⛔ ~~**Can the fitted gDNA prior represent a bimodal landscape?**~~ **ANSWERED: YES**, and it already
   does (§8). Its LOCATION is right too. The successor is `ROADMAP.md` §4.1 — *why is prior fidelity
   anti-correlated with deliverable quality?*
2. ⭐⭐ **Why is the shipped solve so insensitive to pass-0?** A −37.2 % pass-0 change is −3.9 % shipped
   (`TRAPS.md` A15). Until this is understood, every pass-0 ranking in `ROADMAP.md` is suspect.
3. ⭐ **Should the 2026-08-06 change stay at all?** On the deliverable it is −3.9 % net with a +30 %
   stranded regression. An owner call, now well-measured.
4. ⚠ **`γ ≥ 1`** — "off-target gDNA density never exceeds on-target" is assumed by §5's null and is
   *measurable against the origin-split oracle with no solver*. Cheap; not yet done.
5. ⚠ **Seam contamination** — ⛔ DEFERRED BY RULING (owner, 2026-08-06): the simulated panel has none;
   address fragments that disagree with the annotation when real data is reached. Do not let it block a stage.

**What would refute §4:** if the η-invariance gate fails on a real fixture, the premise (density transfer
within a stratum) is wrong and the whole conversion is unfounded.

---

# §13 ⛔ PRUNED — THE HANDOFF IS `SESSION_HANDOFF.md`

This section carried a full session account (the measured state, the proven defects, the tree state, the
traps earned). It is now **one home away** in `docs/calibration/SESSION_HANDOFF.md`, which supersedes it in
full. ⚠ Keeping both was `TRAPS.md` A11's "two homes for one predicate" in documentation form — the same
mistake that deleted the per-issue REGISTER on 2026-08-05.

⭐ **What is kept below is the only part that is a property of the CODE rather than of a session:** the
code-site map, and the traps earned building against it.

## §13.4 The code sites, with line numbers

| what | where |
|---|---|
| ψ's solve entry, and it **already accepts** `lam_imp=(mode,prec)` and `theta_imp=(mode,prec)` | `bp_solver.py:213` → `simplex_logodds.py:592` |
| `mo_g = log(cg·E_g/M)` — ⭐ already `log φ_g`, message-only (D4-clean) | `bp_solver.py:1511` |
| the mass pin (scalar twin) — the `_k` block | `bp_solver.py:1096–1103` |
| `_pin_v` (vector twin) | `bp_solver.py:567` |
| `framed` (vector) · `r` (vector `:1198`, scalar `:982`) | `bp_solver.py:1197` · `:1198`, `:982` |
| `graft` / `lend` (vector) · `_gr` / `_lend` (scalar) · `_pop_l_a` | `bp_solver.py:1203` / `:1229` · `:990` / `:1008` · `:1180` |
| `_relay` / `_transport` — ⛔ **DELIBERATE TWINS, do-not-merge note at** `:929` | `:953` / `:1191` |
| `count_logvar` (the one home) · `composition_logvar` · `transfer_logvar` | `enrichment_frame.py:104` · `:126` · `:479` |
| `own_composition_logvar` · `own_precision` · `build_node_init` | `node_init.py:115` · `:141` · `:205` |
| `g1_locked` · `node_total_density` · `node_global_geometry` | `node_geometry.py:442` · `:342` · `:327` |
| `_JEFFREYS_REF = 0.5` | `simplex_logodds.py:74` |

⛔ **Four names are bound as module globals in MORE THAN ONE module** — `composition_logvar`,
`own_composition_logvar`, `build_node_init`, `_solve_nodes_logodds_all`. Patch every binding and assert the
patch fired (A10). `bp_solver` re-imports all four.

## §13.5 Traps earned in this campaign

1. ⛔ **Never edit `src/` while an arm is running** — each arm is a fresh process importing `src` at start-up.
   Adding a *new* arm to `scripts/` is safe.
2. ⛔ **zsh does not word-split unquoted variables.** Use an array and `"${CONDS[@]}"`. Cost two launches.
3. ⛔ **A wait-loop whose `pgrep` pattern matches its own wrapper deadlocks.** Wait on a log marker.
4. ⚠ **`P0.DEFAULT_SUITE` is the PILOT, not the ladder.** Always pass `--suite .../suite/ladder`.
5. ⚠ **Node-axis and node+edge figures differ by ~2×.** State which one, every time.
6. ⚠ **A composite arm fires only its components' names** — the A10 guard tripped on `zc_ref_prior_damp`
   *after* a complete, valid run. Guard narrowed; check *why* a guard fired before distrusting the data.
