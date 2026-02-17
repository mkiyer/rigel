# FGFR2 Benchmark Report (zero gDNA, strand-specificity=1.0)

Date: 2026-02-16  
Region: chr10:121476640-121601584  
Fragments per run: 100000  
Seeds: 101, 202, 303, 404, 505

## Tag handling confirmation

- BAM SJ tag auto-detection returns: `('ts',)` on rerun FGFR2 alignments (minimap2).
- Current pipeline auto-detects both `XS` (STAR) and `ts` (minimap2), and falls back to strand-agnostic SJ matching when intron strand is unknown.

## Per-seed mean absolute error (MAE)

| Seed | hulkrna | salmon | kallisto | hulkrna/salmon | hulkrna/kallisto |
|---:|---:|---:|---:|---:|---:|
| 101 | 1130.215 | 60.610 | 69.567 | 18.65x | 16.25x |
| 202 | 2501.759 | 264.627 | 156.463 | 9.45x | 15.99x |
| 303 | 2585.889 | 668.373 | 99.254 | 3.87x | 26.05x |
| 404 | 1455.445 | 216.736 | 142.156 | 6.72x | 10.24x |
| 505 | 3192.075 | 303.384 | 430.125 | 10.52x | 7.42x |

## Aggregate metrics across 5 seeds

| Tool | Mean MAE | Median MAE | Mean RMSE | Mean Pearson | Mean Spearman |
|---|---:|---:|---:|---:|---:|
| hulkrna | 2173.077 | 2501.759 | 5326.365 | 0.70815 | 0.55862 |
| salmon | 302.746 | 264.627 | 743.457 | 0.99324 | 0.79776 |
| kallisto | 179.513 | 142.156 | 366.346 | 0.99835 | 0.86361 |

## Summary

- In this rerun cohort, hulkrna remains substantially behind both salmon and kallisto on FGFR2.
- Average MAE gap:
  - hulkrna vs salmon: 9.84x worse
  - hulkrna vs kallisto: 15.19x worse
- Relative ranking is stable: kallisto best, salmon second, hulkrna third.

## Output locations

- Rerun root: /Users/mkiyer/Downloads/hulkrna_runs/fgfr2_depth_runs_rerun_xsts
- Per-seed summaries:
  - /Users/mkiyer/Downloads/hulkrna_runs/fgfr2_depth_runs_rerun_xsts/seed_101/summary.json
  - /Users/mkiyer/Downloads/hulkrna_runs/fgfr2_depth_runs_rerun_xsts/seed_202/summary.json
  - /Users/mkiyer/Downloads/hulkrna_runs/fgfr2_depth_runs_rerun_xsts/seed_303/summary.json
  - /Users/mkiyer/Downloads/hulkrna_runs/fgfr2_depth_runs_rerun_xsts/seed_404/summary.json
  - /Users/mkiyer/Downloads/hulkrna_runs/fgfr2_depth_runs_rerun_xsts/seed_505/summary.json
