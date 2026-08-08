# Validation campaign — handoff

    ⚠ **A DEV DOC. Provisional, not authoritative, and nothing may cite it.** When a finding here settles,
    MOVE it to its permanent home (`ROADMAP.md` for a number, `TRAPS.md` for a lesson, `DESIGN.md` for a
    ruling, `EQUATIONS.md` for a derivation) and delete it here in the same edit. See `docs/dev/README.md`.

    Opened 2026-08-07. ⭐ The settled version of everything is `ROADMAP.md` §1.

## The goal, in one sentence

**Characterise the whole tool end to end before adding anything to it** — with message propagation OFF and
`length_likelihood` OFF, so every number is attributable to the tool as it stands.

⛔ Both switches stay down for the entire campaign. `length_likelihood`'s panels and oracle caches are built
and ready and it is still not time: turning it on mid-campaign makes every other number unattributable.

## Status

| | item | state |
|---|---|---|
| 1 | **calibration-prior-vs-oracle** | ✅ **DONE.** `scripts/design/prior_vs_oracle.py` + 14 gates. Numbers MOVED to `ROADMAP.md` §1.1; lessons MOVED to `TRAPS.md` (`an-equal-length-panel-defeats-the-lift`) |
| 2 | **tool-absolute-accuracy** | 🔧 instrument built — `scripts/design/quant_accuracy.py --arm base` + 9 gates. Panel running |
| 3 | **error-downstream-of-calibration** | 🔧 same instrument, `--arm oracle`. Panel running |
| 4 | **performance** | ⏸ not started. The brief is `ROADMAP.md` §1.0 — ⛔ profile on real cfRNA, not on this panel |

## What is still genuinely open — do not re-derive, these are questions

- ⛔ **The assembler over-calls gDNA by 15.1 % under capture with PERFECT masses in** (`ROADMAP.md` §1.1
  ④). Diagnosed and not fixed. The suspect named by `assemble_priors`' own docstring is the
  support-weighted pooling `Σm/ΣS`, which is exact only where `rho_c` is uniform *inside* the locus —
  and capture puts a strong gradient there. **Nobody has tested that it is the cause.** A cheap probe:
  score `O − F` against a per-locus measure of internal gDNA density spread; if the residual tracks it,
  that is the mechanism.
- ⚠ **`nrna_est` is a THIRD false-positive channel and nothing was watching it.** On `nrna_none` panels
  the true nascent count is exactly 0, and at `g25 ss0.50 capture_on` the EM parked **~1.66 M** fragments
  on synthetic nRNA entities. ⛔ `get_counts_df` DROPS the synthetics, so this is invisible in every
  transcript-level table — `quant_accuracy.py` reports it on the `library` row for that reason. Whether
  it is a defect or the expected behaviour of a nascent candidate with no nascent RNA present has not
  been ruled on.
- **Is `nrna_parent_count` (44,030) meant to differ from `nrna_total` (30,125)?** Not obviously a defect,
  but nobody has stated the relation.
- **`tpm_total_rna` + the nRNA table's `tpm` sums to 1,001,068, not 1,000,000** — 0.107 %, too large for
  rounding, so the two tables normalise over slightly different denominators
  (`estimator.py:504-510`). Small, real, and unexplained.
- ⚠ **`CLAUDE.md`'s script table still cites numbered rule labels** — five distinct ones, in the
  `arm_identity.py` and `backbone_parity.py` rows and in the xfail paragraph. The naming rule bans them
  and `tests/test_no_jargon_labels.py` allowlists that one file, so the gate does not fire there.
  ⛔ Do not paste the labels into this file to list them: the gate DOES cover `docs/dev/`, and it fired
  on exactly that when this bullet was first written. Run
  `python -m pytest tests/test_no_jargon_labels.py -q` with `CLAUDE.md` removed from the allowlist to
  see them. They should be resolved to names, but that means knowing which rule each one WAS, and
  guessing would be worse than leaving them — someone who has the mapping should do it.

## Substrate, all cached and verified

| panel | conditions | oracle cache | what it resolves |
|---|---|---|---|
| `suite/ladder` | 36 | 341 M ✅ | gDNA 0→98 % × stranded/unstranded × capture on/off |
| `suite/flgap_short` | 4 | 32 M ✅ | gDNA ≪ RNA (realised gap −41.0 %) |
| `suite/flgap_long` | 4 | 32 M ✅ | gDNA ≫ RNA (+40.4 %) |
| `suite/pilot` | 8 | scan cache only | chr22 + ERCC, 0 %/100 % gDNA |

⚠ **`quant_accuracy.py` needs `--jobs 2`, not 4, on this machine.** `run_pipeline` holds **7–8.5 GB**
resident per 10 M-fragment condition; four concurrent shards drove the machine into swap and were slower
than two.
