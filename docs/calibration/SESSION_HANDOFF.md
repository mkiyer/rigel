# SESSION HANDOFF — ⭐⭐⭐ THE STRANDED × CAPTURE-ON RESIDUAL

    Written 2026-08-06. ⚠ **WORKING DOC, NOT A PERMANENT FIXTURE.** The six permanent docs are
      `CLAUDE.md`, `docs/SUCCESS.md`, `docs/ROADMAP.md`, `docs/TRAPS.md`, `docs/EQUATIONS.md`,
      `docs/DESIGN.md`, `docs/TESTING.md`.
    ⛔ **DELETE THIS FILE when its task lands**, promoting the numbers into `ROADMAP.md`, the lessons
      into `TRAPS.md`, the derivation into `EQUATIONS.md` and any ruling into `DESIGN.md`.

---

## §1 WHERE THE TREE IS

Branch `fragment-length-gold-standard`. Suite **22 failed / 2290 passed / 3 xfailed** — the 22 are the
standing set (21 `test_golden_output` + the paralog row). A 23rd failure or any other name is a
regression. `ruff check src/ tests/ scripts/` clean.

**A 39 % panel-wide error reduction landed** (`ROADMAP.md` §0). Two bugs, one root cause: the Poisson
log-count variance `1/n` was written out **twice**, and both copies diverge at zero counts.

| where | what it broke |
|---|---|
| `own_precision` | a zero-count object could not FORM its claim — 1,298 intergenic nodes silent |
| `composition_logvar` | `Var(log ρ_tot) = ∞` ⇒ `σ²_transfer = ∞` ⇒ `1/(1/p + ∞) = 0` annihilated **every** message that object sent, on all three streams |

Both now call **one** `count_logvar(n) = trigamma(n + ½)` living in `enrichment_frame` (the leaf module).
`TRAPS.md` C0c/C0d.

---

## §2 ⛔ THE RESIDUAL — one coherent stratum, 20–34 % WORSE

Everything else improved. **`ss_0.99` × `capture_on` at `g10` and above regressed**, monotonically and
without exception:

| condition | Σ\|err\| before | after | |
|---|---|---|---|
| **`g75 ss0.99 capture_on`** | 295,453 | **396,052** | **+34.0 %** ⭐ the worst, start here |
| `g50 ss0.99 capture_on` | 243,215 | 312,390 | +28.4 % |
| `g90 ss0.99 capture_on` | 306,731 | 393,446 | +28.3 % |
| `g25 ss0.99 capture_on` | 145,206 | 183,626 | +26.5 % |
| `g10 ss0.99 capture_on` | 74,407 | 91,383 | +22.8 % |
| `g98 ss0.99 capture_on` | 304,815 | 366,842 | +20.3 % |

Stranded overall: **2,393,419 → 2,606,346 (+8.9 %)**, 7/18 better. That is **10 % of the panel giving
back 213 K against 9.85 M gained elsewhere** — keep the fix, fix the residual.

⭐ **The shape is the diagnosis-in-waiting.** It is absent at `g00`–`g05` (which improve), absent
capture-OFF, and absent unstranded. It needs *stranded* **and** *capture-ON* **and** *real gDNA*
together.

---

## §3 ⚠ THE HYPOTHESIS — UNPROVEN, and it must be measured before it is built on

> **Under capture, an empty intergenic node does not mean "no gDNA here". It means "no probe here".**

The new zero-count anchor says `ρ_g = 0 ± trigamma(½)` wherever a structurally pure-gDNA slot has no
counts. Off capture that is true and it is the strongest statement in the library. **Under capture the
intergenic population is depleted ~24×** (`TRAPS.md` B19, `EQUATIONS.md` §3.6), so an intergenic node can
be empty while the library is 75 % gDNA — and `eff_gdna` does not know about probe depletion, so
`count/E_g` reads a coverage hole as a density of zero.

**Why that would hit exactly this stratum and no other:**

* **unstranded** — objects have no evidence of their own, so even a biased anchor beats ψ's ½ reference.
  It wins there (−44.2 %).
* **stranded** — objects DO have their own strand evidence. The biased anchor now *competes* with a
  correct local answer and drags it down.
* **`g00`–`g05`** — the anchor's claim of ~zero is approximately TRUE however it got there, so no harm.
* **capture-OFF** — intergenic is not depleted, so an empty node really is empty.

⛔ **This is a hypothesis with a mechanism, not a finding.** It has not been localized to a single object
and the alternative explanations below have not been excluded.

**Alternatives that must be ruled out, not assumed away:**

1. The regression is in `composition_logvar`'s **other** term, not the anchor — `σ²_transfer` is now
   finite everywhere, so *every* hop's damping changed, not just zero-count ones. A stranded library has
   more live messages, so it has more to lose from a mis-set transfer variance.
2. `trigamma(n + ½)` vs `1/n` at **small but nonzero** counts (n = 1–3, where they differ 7–9 %) —
   capture-ON intergenic nodes sit exactly there (`B19`: median 0, max 3).
3. A real anchor competing with a real strand solve is a **fusion-weight** problem, not an anchor
   problem — the two might both be right and be combined at the wrong ratio.

---

## §4 ⭐⭐⭐ THE EXPERIMENTS, in order. This is the owner's loop: dissect → localize → cause → fix → test → REPEAT

⛔ **Do not start with a ceiling.** A ceiling prices something that may be unreachable; the loop below
prices the actual defect. (Owner's ruling, 2026-08-05.)

### E1 — dissect the worst scenario

```bash
python scripts/design/worst_objects.py --condition gdna_g75_ss_0.99_nrna_none_capture_on \
  --suite ~/Downloads/rigel_runs/suite/ladder \
  --oracle-cache ~/Downloads/rigel_runs/suite/ladder/oracle_cache --axis both --top 25
```

⭐ **Read the CONCENTRATION curve first.** Concentrated ⇒ a mechanism and a handful of objects that
demonstrate it. Diffuse ⇒ a systematic bias and individual rows are noise. Then `fg_loc` vs `pred_fg`:
if they agree the messages are innocent and the local solve is the defect; if they diverge it is the
relay. ⚠ Rank by error **MASS** (fragments), never by mean error.

### E2 — the decisive A/B, and it tests the hypothesis directly

Restrict the zero-count anchor to slots whose opportunity is **real**, and re-measure. The cleanest
one-thing-varied arm is to withhold the anchor at capture-depleted intergenic slots and see whether the
stranded regression disappears **while the unstranded win survives**. If it does, the hypothesis is
confirmed and the fix is an opportunity model. If the regression persists, §3's alternatives 1–3 are
where it lives.

⭐ Cheap pre-flight, no solver: census the intergenic counts at `g75 ss0.99 capture_on` the way
`anchor_probe` did at `g00` — how many of the 1,298 are empty, what is `Σ E_g`, and what does the TRUE
gDNA density say the anchor *should* be? If the true density is far from zero at empty nodes, §3 is
confirmed before a single solve runs.

### E3 — the fix, once the cause is localized

The likely shape, if §3 holds: **`eff_gdna` must account for probe retention before a zero count
licenses a zero density.** ⛔ Derive it; do not tune it. And note the standing refusal — `EQUATIONS.md`
§3.5c records that the affine-in-overlap extrapolation `e_g[exon] = 2·e_g[EDGE] − e_g[off-probe]` is
**exactly the simulator's own retention law**, so fitting it would be scoring against the substrate that
generated it. Whatever is built must be derivable from the probe annotation, not from the panel.

### E4 — gates, then the panel, then land

`TRAPS.md` A2: falsification first, verified failing; then break the fixed code and watch **each** gate
fire. B14: name the observable for *each place* the change is made. ⛔ B18: run the 36-condition arm
before writing into `src/` — and score `mwae_all` + raw `Σ|err|`, never `solv%`/`mwae`/`conf-wrong`
(A12/A12b). Zero controls on every arm (`g00` and the `g98` end), and check the arm could have fired
(A14).

### E5 — REPEAT

Re-run the panel, take the new worst stratum, and go again. Two iterations of exactly this loop produced
the 39 %.

---

## §5 THE INSTRUMENTS YOU NEED

| | |
|---|---|
| `solvability_audit.py --suite` | the 36-condition table. ⭐ **Read `mwae_all` and `Σ\|err\|`** — the last two columns, fixed-denominator. ~15 min with the oracle cache |
| `worst_objects.py --condition` | step 2: one condition dissected to individual objects, ranked by error MASS |
| `certified_q_census.py` | the pattern for a truth-only census: reads the oracle cache, no solver, whole panel in seconds |
| `zero_controls.py` | ⛔ owner requires both arms on every experiment |
| `toy_harness.py` / `toy_trace_error.py` | isolate a mechanism once you have one |

⚠ The oracle cache at `~/Downloads/rigel_runs/suite/ladder/oracle_cache` is **VALID — do not delete.** It
depends only on the accumulator and the index, never on calibration.

---

## §6 ⛔ WHAT IS SETTLED — do not re-open

| | verdict |
|---|---|
| **the certified-RNA λ term** (`(½+S)·log(1−f_g)`) | ⛔ **REFUTED.** The splice-visibility `q` has median 0.19–0.71, so the term §2d dropped is the same size as the one kept; with `q` free the profile likelihood is exactly flat and the count carries ONE BIT. Worse than the uninformative reference on 12/36. `TRAPS.md` C0b, `certified_q_census.py`, `test_certified_rna_licence.py`. ⭐ The channel is *blocked on a missing exonic-reach opportunity*, not impossible |
| **the vertex problem** | ⛔ closed — the value of missing information, not headroom. `ROADMAP.md` §2 |
| **`Var(κ̂)` / the strand deadband** | ⛔ closed as an accuracy item (a wash); it was a REPORTING defect and is fixed instrument-side. Three solver designs refuted. `TRAPS.md` A12b |
| **a resolving-power floor on `τ`** | ⛔ derived, implemented, refuted — `τ` is continuous, so any floor is a tuned constant |
| **`length_likelihood`** | ⛔ stays OFF; the ladder keeps equal lengths and κ = ½ deliberately |
