# The capture weight as the posterior expectation of the contraction, not the plug-in mode — derivation (2026-10-11, sandbox, not authoritative)

**Status: BOTH arms REFUTED by their A/Bs (2026-10-11): the expectation of the density E[e^x] (§0) and the posterior
mean of log density E[x] (§0c). The mode is the location-faithful readout; a posterior average moves the location
instead of reducing the variance. The record is `ISSUES: the-capture-weights-are-noisy-on-a-capture-off-library`.**

## 0. The A/B of E[e^x] (prototype site1, both arms from one build, pinned, fractional)

| panel · stratum | mode arm → expectation arm, transcripts % (genes %) |
|---|---|
| ladder · stranded × OFF | 7.52 (0.53) → 15.57 (1.01) |
| ladder · stranded × ON | 9.71 (1.87) → 14.11 (2.30) |
| ladder · unstranded × OFF | 2.07 (0.24) → 11.26 (0.74) |
| ladder · unstranded × ON (deferred) | 40.42 (15.40) → 22.25 (5.16) |
| RNA-long · stranded × OFF / ON | 3.52 → 4.66 / 6.09 → 7.57 |
| RNA-short · stranded × OFF / ON | 4.99 → 7.07 / 14.24 → 15.39 |
| both gap panels · unstranded × ON (deferred) | 34.11 → 16.00 / 68.06 → 30.93 |
| VCaP gDNA fraction (truth 0.2518) | 0.2528 → 0.2489; weights q10/q50/q90 0.366/1/425 → 0.0091/1/6.6 |

Every in-scope row is worse; g00 unstranded OFF, which was fine, goes 1.59 → 23.39 %. The weights' spread says why:
`E[e^x]` under a posterior as broad as the prior is dominated by the prior's HIGH tail (Jensen: `E[e^x] ≥ e^{E[x]}`,
by orders of magnitude over a 20-nat support), so every shallow object was pulled ABOVE the deep ones, not toward the
typical level; the published unit (the median) then sat among the inflated shallow objects and the deep ones fell to
a hundredth of it. §2's claim that a shallow object shrinks "toward the prior's mean contraction, the one level" was
wrong in its scale: the prior's typical level is its mode or median in LOG density, not its mean of `e^x`. The
deferred stratum's large gain (and the same on both gap panels) is recorded, not pursued here.

## 0c. The A/B of E[x] (prototype site2, the same protocol)

| panel · stratum | mode arm → E[x] arm, transcripts % (genes %) |
|---|---|
| ladder · stranded × OFF | 7.52 (0.53) → 7.86 (0.52); the target row g00 13.61 → 9.50, but g05 2.72 → 7.38, g50 4.13 → 5.17 |
| ladder · stranded × ON | 9.71 (1.87) → 13.00 (2.66) |
| ladder · unstranded × OFF | 2.07 (0.24) → 4.16 (0.37); g00 1.59 → 4.05 |
| ladder · unstranded × ON (deferred) | 40.42 → 46.68; g50's gDNA fraction 0.343 → 0.164 |
| RNA-long · stranded × OFF / ON | 3.52 → 4.45 / 6.09 → 12.42 |
| RNA-short · stranded × OFF / ON | 4.99 → 7.28 / 14.24 → 23.83 |
| VCaP gDNA fraction (truth 0.2518) | 0.2528 → 0.2412; weights q10/q50/q90 0.366/1/425 → 0.091/1/311 |

The mirror image of §0: the mean of log density under a broad, one-sided posterior sits BELOW the mode (the "why the
mode" objection, 0.14–0.79 nats at a zero-read object, here on every shallow object), so shallow objects lose
enrichment they have, and every capture-ON stratum pays for it; the one row that improves is the one whose noise was
symmetric. Together the two arms say the same thing: the mode is the location-faithful readout of a one-sided
posterior, and replacing it by any posterior average moves the LOCATION of every shallow object (up under E[e^x],
down under E[x]) when what is wrong on a capture-OFF field is the VARIANCE of the modes across shallow objects. A
repair must reduce that variance without moving the location — shrinkage of the mode toward a LOCAL reference by the
posterior's width, which is the owner-approved local capture prior / local footprint direction
(`DESIGN.md` §7.2, 2026-10-07 and 2026-10-09) — and is a design of its own, not a readout change.

## 0b. The arm as run: the posterior mean of log density

The same uncertainty-shrinkage on the scale where the prior's typical level lives: publish
`w_o = exp(E_o[x] − median_o E[x])`, `E_o[x] = Σ_i x_i P_o(x_i)`. A variational EM maximises `E_P[log Pois(k; θ L(c))]`,
which carries `E[log c]` where one object dominates a length, so this is the derived form under that view; the
marginal-likelihood view that gave `E[c]` trusts the prior's tail literally and is what the A/B refused. The "why the
mode" ruling's objection (the mean under-reads a one-sided zero-read posterior by 0.14–0.79 nats) is measured on the
zero-read objects as before. Gates unchanged in form, the reference `Σ x_i P_i`.
 The owner chose this as option 1 on 2026-10-11, after the ladder's re-score and the
attribution below. Nothing enters the source before a prototype outside the tree is A/B'd per stratum against what ships.

## 1. What is wrong, measured

The reader publishes every object's capture weight as `w_o = exp(x*_o − median x*)`, `x*_o` the posterior MODE of its
log gDNA density under the landscape (`DESIGN.md` §7.2, `EQUATIONS.md` §11). On a capture-OFF library the true
weights are one level; the published ones are each object's own noise around it, largest where objects are shallow.
The EM's effective length of a transcript is the sum of its conserved shares at those weights
(`L_T = Σ S_T(p)·w_p + Σ M_T(k)·w_k`), so two isoforms sharing most objects but differing in a few shallow ones get
random length ratios, and mass moves to silent isoforms. Attribution (pinned, fractional, one condition, the
weights forced to one in a prototype build): the ladder's g00 stranded OFF row, transcripts 13.61 → 1.93 %, genes
0.83 → 0.14 %, the mass on silent isoforms 251,216 → 17,058, calibration's gDNA estimate unchanged (1,126 vs 1,112).
The weights own the whole of that row's regression; the gene level holds because the error is isoform allocation.

## 2. The derivation

At every object the reader already has the full posterior over the landscape's grid, `P_o(x_i) ∝ exp(log L_o(x_i) +
log P(x_i))`, before it takes the argmax. The weight is a function of that density, `c = ρ/ρ_ref = e^{x − x_ref}`, and
the EM uses it linearly: a transcript's fragment rate is `θ_T · L_T` with `L_T` linear in the `c`'s. Under
uncertainty about `c_o`, the rate's expectation is `θ_T · Σ S_T(o)·E[c_o]`, so the quantity the EM should read is the
posterior EXPECTATION of each object's contraction,

    E[c_o] = Σ_i e^{x_i − x_ref} · P_o(x_i),        published as  w_o = E[c_o] / median_o E[c_o],

the same unit convention as today (the typical read object is 1; the EM is invariant to the unit). No new constant, no
detector, no reference density: the landscape is the prior it already is, and the expectation is taken under the
object's own likelihood. What it does, by the shape of `P_o`:

- **a deep object** (many reads): `P_o` is a spike at the mode; `E[c] = e^{x*}` to within the spike's width — the
  weight is unchanged where the reader is right today;
- **a shallow object**: `P_o` is broad, close to the prior `P`; `E[c]` shrinks toward the prior's own mean
  contraction, which on a capture-OFF field is the one level — the noise collapses exactly where it comes from;
- **a zero-read object**: `P_o ∝ P(x)·e^{−ρEg}` is one-sided; the ruling for the MODE (`EQUATIONS.md` §11, "why the
  mode") was that the mean of `x` under-reads the background by 0.14–0.79 nats as it thins. `E[e^x]` is a different
  quantity (Jensen: `E[e^x] ≥ e^{E[x]}`), and the question is open: it is measured, not argued, on the ladder's zero-read
  objects against the simulator's own yield (`ruler_vs_truth.py --scale`, the unprobed class).

What it costs: nothing at a read slot — the curve is evaluated at every grid node already, and the expectation is one
pass over it; at a zero-read slot the coarse-to-fine search evaluates a subset, so either the expectation is taken
over the evaluated nodes with the dropped ones' mass bounded (they lie more than the stride from every local
maximum), or those slots take the full grid (≈ 250 CPU-s on VCaP at today's node cost, +15 s of wall). The prototype
measures both.

## 3. What it may cost on a captured field, to be measured

On a capture-ON library the landscape is bimodal (the background and the enriched mode). A shallow PROBED object's
posterior straddles both, so its `E[c]` is pulled toward the background by its own uncertainty, where the mode would
have picked the enriched side. That is a real loss where it happens and the honest statement of what a shallow object
knows. The A/B reads both ON strata beside both OFF strata, per stratum, and the gap panels beside the ladder; a
mechanism that buys the OFF rows by selling the ON rows is not taken on this derivation alone.

## 4. Gates and the A/B

| gate | holds |
|---|---|
| `expectation-is-the-reference` | the native `E[c]` against a NumPy expectation over the same curve on the 48-slot and the ladder fixtures, 1e-9; a deliberate defect (the prior left out) fires |
| `deep-object-identity` | on every slot whose posterior sd is below a tenth of a nat, `log E[c]` equals the mode to the spike's width |
| `unit-invariance` | the EM's counts and TPM invariant to the unit (the existing capture-reader gate) |
| the A/B | pinned, fractional: the ladder (16 rows, per stratum), the test chromosome, both gap panels, the VCaP mix against its read-name truth; `ruler_vs_truth.py --scale` for the probed and unprobed classes apart |

## 5. Stages

0. The prototype in a worktree: the kernel returns the expectation beside the mode (the standalone binding and the
   block solve), `calibrate.capture_weights` publishes it behind a prototype switch; the gates.
1. The A/B, every panel, per stratum; the zero-read question read on the unprobed class.
2. The ruling, with the numbers; then the source, by the move rule — `EQUATIONS.md` §11's "why the mode" paragraph
   is the one that changes.
