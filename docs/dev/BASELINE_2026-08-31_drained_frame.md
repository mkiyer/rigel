# BASELINE — 2026-08-31, the drained frame

    ⚠ **A DEV DOC — a dated measurement snapshot, not a source of current numbers.** The way to get a
    current number is to run the instrument; this snapshot exists so the next session can see what
    moved without re-running everything. Substrate: the sparse-nascent 16-condition ladder +
    30-condition test chromosome, slot_truth certified in the DRAINED frame the same day. Tree: the
    frame-migration commits of 2026-08-31 (the ruling now lives at `DESIGN.md` §4.3). Every
    number below is drained-frame; ⛔ nothing measured before 2026-08-31 is comparable to these.

## ① The 0.8.0 metric — `calibration_vs_oracle.py`, all 16 ladder conditions (512 s)

Σ|Δ| fragments (region + boundary) against oracle calibration, per stratum (contaminated rows):

| stratum | Σ|Δ| frags | under / EXACT / over (objects) |
|---|---|---|
| stranded × capture OFF (IN SCOPE) | 152,695 | 69,580 / 7,059 / 24,235 |
| stranded × capture ON (IN SCOPE) | 624,999 | 56,987 / 6,215 / 9,931 |
| unstranded × capture OFF (IN SCOPE) | 179,990 | 67,559 / 7,087 / 26,100 |
| unstranded × capture ON (DEFERRED) | 2,428,329 | 58,476 / 6,256 / 8,249 |
| g00 zero controls (all strata) | 152,898 | 0 / 2,745 / 67,571 |

- ⭐ **Worst IN-SCOPE condition: `g98 ss.99 capture-ON`, Σ|Δ| = 456,838 — 450,025 UNDER vs 6,813
  over.** The systematic shape on every stratum: gDNA is UNDER-called at near-pure objects
  (shortfall +0.15…+0.29 at true f_g ∈ [0.999, 1], both axes) and mildly over-called at [0, 0.2].
  That is rank-1 dissection's target.
- Library `f_gdna` P vs O: in-scope within ~0.002–0.034; the known deferred blindness intact
  (`g98 ss.50 ON`: 0.796 vs 0.983).
- Zero controls: over-call only, 7.4 k (`g00 ss.99 OFF`) to 58.8 k (`g00 ss.50 OFF`) fragments.
- Frame extras now on every row: `gdna_spliced_leak`, `lift_n_ambiguous`.

## ② The prior the EM reads — `prior_vs_oracle.py`, all 16 (1,301 s)

| stratum | mwae φ (P vs O) | gdna_eff_len w-rel-err |
|---|---|---|
| stranded × OFF | 0.0150 | 0.037 |
| stranded × ON | 0.0254 | 0.031 |
| unstranded × OFF | 0.0166 | 0.030 |
| unstranded × ON (DEFERRED) | 0.1424 | 0.397 |
| g00 controls | 0.0070 | — |

- Prior gDNA totals P/O: in-scope 0.954–1.053 except `g50 ss.50 ON` 0.741 and `g98 ss.50 ON` 0.811
  (both DEFERRED); zero-control false-positive priors 5,132 → 112,706 fragments (ss.99 OFF → ss.50 ON).
- The drained-frame report, whole ladder: 3,150,429 held fragments; lift ambiguity 110,439 (3.51 % of
  held, ≈0.11 % of the library); the drain deposited 13,513 spliced records into the gdna partition
  across 12/16 conditions (`ISSUES: drain-contaminates-certified-rna`).

## ③ The policy benchmark — `policy_benchmark.py`, test (30) + ladder (16)

- **Rung-0 identity holds in the drained frame**: `message` ≡ `silent` on all 46 conditions.
- The two bars (ladder): **unstranded — relay beats silence 5/8, worst loss 1.49×**; **stranded —
  relay beats silence 2/8, worst harm 1.65×** (`g98 ss.99 OFF`, 126,467 vs 208,306). Relay's value
  stays concentrated where the local solve is blind (deferred rows 0.20× at `g98 ss.50 ON`; all four
  zero controls win).
- ⚠ The stranded-harm bar reads WORSE than the pre-rebuild records — re-derived, not comparable to
  them; this is the standing number for the message-policy charter now.

## ④ The thermometer + the attribution floor — `quant_accuracy.py` base / base_reseed, all 16

- Reseed floor per condition (|base − reseed| transcript Σ|Δ|): ~30 (`g98 ss.99 OFF`) to ~2,700
  (`g05 ss.50 ON`); library-fraction floor ~1e-4. ⛔ No `quant_accuracy` delta below the same-session
  floor is attributable (JSONs: session scratchpad `baseline_qa_base/reseed.json`).
- Transcript Σ|Δ| spans 90,204 (`g98 ss.99 OFF`) to 4.63 M (`g00 ss.99 ON`); the zero controls carry
  the largest misassignment (3.1–4.6 M) — the standing EM-assignment problem, not calibration.
- Library thermometer: in-scope within 0.002–0.06 of truth; deferred `g98 ss.50 ON` 0.757 vs 0.980.

## ⑤ Solvability — `solvability_audit.py`, whole panel

Read `mwae_all` / `Σ|err|` only (the other columns share a solver-moved denominator):

- In-scope contaminated rows: `mwae_all` 0.0075 (`g05 ss.99 OFF`) → 0.0413 (`g50 ss.50 OFF`);
  `Σ|err|` 41,016 → 304,662. Deferred rows 0.128–0.298 / up to 1.76 M, blind as recorded.
- Zero controls (false-positive checks): `Σ|err|` 6,389 (`g00 ss.99 OFF`) → 159,015 (`g00 ss.50 ON`).
- ⭐ **The relay HURTS the solvable set on 9/12 contaminated conditions** (mean Δ +0.0012) — the
  policy-charter reading, consistent with ③'s bars; and the declared precision is NOT earned
  (ratio > 1) on 11/12 rows.
- ⚠ This instrument's worst in-scope row (`g50 ss.50 OFF` by `Σ|err|`) differs from ①'s
  (`g98 ss.99 ON` by oracle-mass Σ|Δ|) — the two weigh objects differently; the dissection loop
  should look at both before choosing its scenario.
