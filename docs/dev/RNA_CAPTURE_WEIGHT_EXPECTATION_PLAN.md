# The capture weight as the posterior expectation of the contraction, not the plug-in mode — derivation (2026-10-11, sandbox, not authoritative)

**Status: DERIVED, nothing built.** The owner chose this as option 1 on 2026-10-11, after the ladder's re-score and the
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
