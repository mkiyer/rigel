# PR06 Scan Memory Cap Sweep

Date: 2026-05-13

Workload: VCAP RNA20M + gDNA20M

```bash
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --stages --threads 8 --scan-bgzf-threads 2 \
  --scan-buffer-size <cap> --scan-read-name-batch-size 512 \
  --memory-interval 250
```

Output directories:

- `/Users/mkiyer/Downloads/rigel_runs/profile_pr06_cap_4gib`
- `/Users/mkiyer/Downloads/rigel_runs/profile_pr06_cap_2gib`
- `/Users/mkiyer/Downloads/rigel_runs/profile_pr06_cap_1gib`

## Results

| Scan cap | Wall | Peak RSS | After scan RSS | Scan wall | Router wall | Spilled chunks | Buffer peak | CSR bytes | Partition bytes |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 GiB | 78.52 s | 16.05 GB | 9.05 GB | 35.47 s | 17.04 s | 14 | 4.13 GiB | 5.99 GiB | 5.99 GiB |
| 2 GiB | 80.61 s | 15.60 GB | 7.22 GB | 35.86 s | 18.08 s | 23 | 2.07 GiB | 5.99 GiB | 5.99 GiB |
| 1 GiB | 77.38 s | 15.49 GB | 6.50 GB | 36.16 s | 15.80 s | 28 | 1.08 GiB | 5.99 GiB | 5.99 GiB |

## Interpretation

The 2 GiB cap has no meaningful wall-time penalty in the full pipeline:
`80.61 s` versus `78.52 s` at 4 GiB, a 2.7% difference and inside the
ordinary run-to-run noise seen in this profiling setup.

The scan-stage RSS reduction is the main win. RSS immediately after scan
drops from 9.05 GB to 7.22 GB at 2 GiB, a 1.83 GB reduction. The new
PR05 `buffer_summary.memory_bytes_peak` metric shows the resident buffer
high-water dropping from 4.13 GiB to 2.07 GiB.

Total process peak RSS drops only 0.45 GB because the peak now occurs
after `fragment_router_scan`, where the global scoring CSR and partition
arrays dominate. PR07 (float32 scored-fragment payloads) and PR08 (fused
partition scatter) are therefore the right next memory targets.

## Decision

Use 2 GiB as the default scan buffer cap. It reduces the scan plateau
substantially, does not materially slow the full run, and keeps a
conservative margin versus the more aggressive 1 GiB low-memory setting.