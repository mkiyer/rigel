# Full-panel validation: frozen count candidate

*2026-10-07. Paired validation of the isolated count prototype; no new production source change.*

*Integration note, 2026-10-09: this count model is now in main's working tree. The
[landing audit](RNA_COUNT_READINESS_CHECKPOINT.md#final-engineering-landing) verifies the
cleanup chain and reproduces one condition per stratum through the ordinary main package.
The tables below remain the original receipts; the capture reader is unchanged.*

The candidate removes the intron factory, restores intron strand messages and makes those messages exact conditional strand likelihoods. The final psi observation model, point-count landscape and current capture reader are unchanged. The new density-evidence reference is a separate experiment.

All **54** end-to-end comparisons are complete: 30 test-chromosome conditions, 16 ladder conditions and four conditions in each length-gap panel. The earlier calibration-only screen covers all 24 ladder/gap conditions, with object classes and admitted strands separate. Those count findings are in [ISSUES.md](../ISSUES.md#the-gdna-prior-enters-psi-twice).

Cells are **transcript / gene absolute count error (%)**. Each contaminated summary combines only DNA levels within the same strand/capture stratum; controls, strand fidelities and capture states remain separate. Every individual condition is also shown. Unstranded capture-ON is deferred for 0.8.0 and remains a robustness concern.

Current and candidate runs use the same frozen source baseline, pinned scan threads, fractional assignment, seed zero and default calibration/EM threads. The candidate binary is selected and verified in every child. The 0.7.1 CLI receipts used fractional assignment and one thread; their truth fingerprints and RNA denominators match these runs. Timings are not compared across those different thread settings.

These are count-candidate results with the old reader. They do not establish a detector-free release or remove the remaining evidence, capture-prior and opportunity work. Individual high-DNA test conditions still fail 0.7.1, and all individual rows remain visible below. Real libraries and production landing checks are separate requirements.

## test_base

| Stratum | 0.7.1 | Current | Candidate |
|---|---:|---:|---:|
| ss0.50 OFF — contaminated | 9.25 / 1.76 | 8.66 / 1.33 | 8.81 / 1.34 |
| ss0.50 OFF — zero gDNA | 9.38 / 1.49 | 5.71 / 1.01 | 5.67 / 1.01 |
| ss0.50 ON — contaminated (deferred) | 18.94 / 10.31 | 15.97 / 5.21 | 11.96 / 3.22 |
| ss0.50 ON — zero gDNA (deferred) | 17.93 / 7.39 | 9.23 / 2.44 | 9.24 / 2.44 |
| ss0.70 OFF — contaminated | 9.44 / 1.54 | 8.62 / 1.25 | 8.00 / 1.25 |
| ss0.70 OFF — zero gDNA | 10.06 / 1.47 | 6.59 / 0.90 | 6.35 / 0.90 |
| ss0.70 ON — contaminated | 10.28 / 3.32 | 7.69 / 1.69 | 7.60 / 1.67 |
| ss0.70 ON — zero gDNA | 13.51 / 6.20 | 8.27 / 2.16 | 8.55 / 2.16 |
| ss0.99 OFF — contaminated | 8.20 / 1.20 | 7.67 / 1.08 | 7.62 / 1.08 |
| ss0.99 OFF — zero gDNA | 11.35 / 1.26 | 7.71 / 0.83 | 6.99 / 0.83 |
| ss0.99 ON — contaminated | 8.58 / 2.45 | 7.58 / 0.96 | 7.19 / 0.93 |
| ss0.99 ON — zero gDNA | 14.22 / 5.83 | 7.35 / 0.67 | 7.35 / 0.67 |

### Every condition

| Condition | 0.7.1 | Current | Candidate |
|---|---:|---:|---:|
| g00 ss0.50 OFF | 9.38 / 1.49 | 5.71 / 1.01 | 5.67 / 1.01 |
| g00 ss0.50 ON | 17.93 / 7.39 | 9.23 / 2.44 | 9.24 / 2.44 |
| g00 ss0.70 OFF | 10.06 / 1.47 | 6.59 / 0.90 | 6.35 / 0.90 |
| g00 ss0.70 ON | 13.51 / 6.20 | 8.27 / 2.16 | 8.55 / 2.16 |
| g00 ss0.99 OFF | 11.35 / 1.26 | 7.71 / 0.83 | 6.99 / 0.83 |
| g00 ss0.99 ON | 14.22 / 5.83 | 7.35 / 0.67 | 7.35 / 0.67 |
| g05 ss0.50 OFF | 8.37 / 1.28 | 7.88 / 1.09 | 7.93 / 1.09 |
| g05 ss0.50 ON | 14.04 / 4.39 | 14.80 / 2.60 | 10.80 / 1.14 |
| g05 ss0.70 OFF | 6.92 / 1.15 | 6.86 / 0.98 | 6.24 / 0.98 |
| g05 ss0.70 ON | 7.67 / 1.51 | 6.32 / 0.75 | 5.73 / 0.77 |
| g05 ss0.99 OFF | 6.96 / 0.97 | 6.55 / 0.89 | 6.45 / 0.89 |
| g05 ss0.99 ON | 6.03 / 1.15 | 5.23 / 0.35 | 4.91 / 0.35 |
| g25 ss0.50 OFF | 8.78 / 1.64 | 7.99 / 1.16 | 8.14 / 1.16 |
| g25 ss0.50 ON | 14.54 / 7.67 | 14.69 / 5.84 | 11.12 / 3.50 |
| g25 ss0.70 OFF | 10.40 / 1.44 | 8.65 / 1.14 | 7.57 / 1.15 |
| g25 ss0.70 ON | 10.33 / 3.69 | 6.49 / 1.34 | 6.75 / 1.26 |
| g25 ss0.99 OFF | 7.65 / 1.03 | 7.22 / 0.98 | 7.27 / 0.98 |
| g25 ss0.99 ON | 7.77 / 2.38 | 7.38 / 0.80 | 6.69 / 0.77 |
| g50 ss0.50 OFF | 10.23 / 2.14 | 9.99 / 1.53 | 10.33 / 1.54 |
| g50 ss0.50 ON | 27.38 / 18.98 | 17.25 / 6.91 | 12.45 / 4.42 |
| g50 ss0.70 OFF | 11.47 / 1.84 | 10.80 / 1.40 | 10.84 / 1.40 |
| g50 ss0.70 ON | 13.35 / 4.76 | 9.64 / 2.39 | 9.85 / 2.24 |
| g50 ss0.99 OFF | 10.27 / 1.43 | 9.42 / 1.17 | 9.28 / 1.16 |
| g50 ss0.99 ON | 12.97 / 4.07 | 10.21 / 1.23 | 10.17 / 1.20 |
| g98 ss0.50 OFF | 44.35 / 19.27 | 38.02 / 14.22 | 38.13 / 14.74 |
| g98 ss0.50 ON | 205.06 / 173.74 | 87.49 / 62.31 | 86.66 / 61.50 |
| g98 ss0.70 OFF | 42.05 / 16.09 | 36.83 / 13.94 | 36.59 / 13.47 |
| g98 ss0.70 ON | 55.98 / 39.60 | 69.50 / 42.44 | 71.94 / 45.54 |
| g98 ss0.99 OFF | 35.65 / 12.88 | 34.51 / 11.52 | 34.63 / 11.62 |
| g98 ss0.99 ON | 50.71 / 26.07 | 60.39 / 29.04 | 59.15 / 28.26 |

## ladder

| Stratum | 0.7.1 | Current | Candidate |
|---|---:|---:|---:|
| ss0.50 OFF — contaminated | 10.29 / 1.49 | 2.29 / 0.31 | 2.30 / 0.31 |
| ss0.50 OFF — zero gDNA | 20.09 / 2.53 | 1.60 / 0.14 | 1.59 / 0.14 |
| ss0.50 ON — contaminated (deferred) | 33.21 / 14.25 | 12.40 / 4.47 | 7.64 / 3.75 |
| ss0.50 ON — zero gDNA (deferred) | 10.22 / 3.18 | 7.27 / 3.22 | 7.27 / 3.22 |
| ss0.99 OFF — contaminated | 3.95 / 0.44 | 2.02 / 0.25 | 2.03 / 0.26 |
| ss0.99 OFF — zero gDNA | 3.70 / 0.47 | 1.93 / 0.14 | 1.93 / 0.14 |
| ss0.99 ON — contaminated | 9.91 / 2.52 | 3.88 / 1.38 | 3.88 / 1.38 |
| ss0.99 ON — zero gDNA | 20.22 / 3.47 | 6.47 / 2.36 | 6.46 / 2.36 |

### Every condition

| Condition | 0.7.1 | Current | Candidate |
|---|---:|---:|---:|
| g00 ss0.50 OFF | 20.09 / 2.53 | 1.60 / 0.14 | 1.59 / 0.14 |
| g00 ss0.50 ON | 10.22 / 3.18 | 7.27 / 3.22 | 7.27 / 3.22 |
| g00 ss0.99 OFF | 3.70 / 0.47 | 1.93 / 0.14 | 1.93 / 0.14 |
| g00 ss0.99 ON | 20.22 / 3.47 | 6.47 / 2.36 | 6.46 / 2.36 |
| g05 ss0.50 OFF | 12.03 / 1.34 | 1.87 / 0.17 | 1.88 / 0.17 |
| g05 ss0.50 ON | 28.37 / 5.60 | 10.90 / 1.40 | 4.95 / 1.41 |
| g05 ss0.99 OFF | 3.51 / 0.23 | 1.56 / 0.15 | 1.57 / 0.15 |
| g05 ss0.99 ON | 9.34 / 2.03 | 3.06 / 0.57 | 3.03 / 0.58 |
| g50 ss0.50 OFF | 5.43 / 0.73 | 2.49 / 0.35 | 2.50 / 0.35 |
| g50 ss0.50 ON | 29.03 / 17.78 | 10.26 / 5.50 | 9.55 / 5.26 |
| g50 ss0.99 OFF | 3.82 / 0.49 | 2.35 / 0.28 | 2.37 / 0.28 |
| g50 ss0.99 ON | 9.91 / 2.72 | 4.46 / 2.15 | 4.49 / 2.14 |
| g98 ss0.50 OFF | 49.48 / 27.30 | 16.87 / 5.65 | 17.51 / 5.96 |
| g98 ss0.50 ON | 367.97 / 337.04 | 137.01 / 125.19 | 88.07 / 77.79 |
| g98 ss0.99 OFF | 27.97 / 8.73 | 15.12 / 4.41 | 15.17 / 4.45 |
| g98 ss0.99 ON | 37.32 / 20.96 | 28.60 / 20.51 | 28.54 / 20.51 |

## flgap_rna_short

| Stratum | 0.7.1 | Current | Candidate |
|---|---:|---:|---:|
| ss0.50 OFF — contaminated | 6.00 / 0.61 | 3.79 / 0.27 | 3.81 / 0.27 |
| ss0.50 ON — contaminated (deferred) | 29.76 / 6.19 | 17.57 / 2.23 | 18.47 / 2.32 |
| ss0.99 OFF — contaminated | 6.31 / 0.44 | 3.53 / 0.23 | 3.54 / 0.23 |
| ss0.99 ON — contaminated | 20.33 / 3.52 | 13.33 / 1.40 | 13.31 / 1.42 |

### Every condition

| Condition | 0.7.1 | Current | Candidate |
|---|---:|---:|---:|
| g50 ss0.50 OFF | 6.00 / 0.61 | 3.79 / 0.27 | 3.81 / 0.27 |
| g50 ss0.50 ON | 29.76 / 6.19 | 17.57 / 2.23 | 18.47 / 2.32 |
| g50 ss0.99 OFF | 6.31 / 0.44 | 3.53 / 0.23 | 3.54 / 0.23 |
| g50 ss0.99 ON | 20.33 / 3.52 | 13.33 / 1.40 | 13.31 / 1.42 |

## flgap_rna_long

| Stratum | 0.7.1 | Current | Candidate |
|---|---:|---:|---:|
| ss0.50 OFF — contaminated | 5.29 / 0.59 | 2.16 / 0.22 | 2.14 / 0.23 |
| ss0.50 ON — contaminated (deferred) | 28.99 / 5.84 | 6.55 / 0.94 | 43.22 / 11.26 |
| ss0.99 OFF — contaminated | 5.23 / 0.46 | 2.27 / 0.19 | 2.26 / 0.19 |
| ss0.99 ON — contaminated | 7.38 / 1.62 | 4.23 / 0.66 | 4.26 / 0.68 |

### Every condition

| Condition | 0.7.1 | Current | Candidate |
|---|---:|---:|---:|
| g50 ss0.50 OFF | 5.29 / 0.59 | 2.16 / 0.22 | 2.14 / 0.23 |
| g50 ss0.50 ON | 28.99 / 5.84 | 6.55 / 0.94 | 43.22 / 11.26 |
| g50 ss0.99 OFF | 5.23 / 0.46 | 2.27 / 0.19 | 2.26 / 0.19 |
| g50 ss0.99 ON | 7.38 / 1.62 | 4.23 / 0.66 | 4.26 / 0.68 |

## Reproduction and remaining work

Initial receipts, source/native identity, calibration tables and the first 36 comparisons are in `.cache/rigel_runs/2026-10-07_full_panel_counts/`. The continuation is in `.cache/rigel_runs/2026-10-07_release_followup/`: commands, logs, per-run wall times/load, calibration-hook receipts, batch manifests and `panel_comparison.json`. `report_all.py` checks completeness, uniqueness, assignment settings, truth fingerprints and denominators before rendering this report.

The implementation and falsification checks for the count candidate are in [RNA_FACTORY_REMOVAL_CHECKPOINT.md](RNA_FACTORY_REMOVAL_CHECKPOINT.md). The separate evidence-interface implementation is in [RNA_MESSAGE_EVIDENCE_CHECKPOINT.md](RNA_MESSAGE_EVIDENCE_CHECKPOINT.md). Do not conflate their acceptance gates or accuracy claims.

All four serial real-library pairs are now complete; see [RNA_REAL_LIBRARY_VALIDATION.md](RNA_REAL_LIBRARY_VALIDATION.md). The subsequent [production-readiness audit](RNA_COUNT_READINESS_CHECKPOINT.md) ran the full unchanged suite and preflight on the frozen candidate. It identifies the required reference migration and a weak-strand stress-case count regression; actual source-landing checks are complete in the landing audit above. The detector-free reader and prior-once partner retain their specified owner checkpoints. No commit or push has been made.
