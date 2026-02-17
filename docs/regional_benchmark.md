# Regional Benchmarking: hulkrna vs salmon vs kallisto

Use `scripts/benchmark_region_competition.py` to run focused-region, head-to-head
quantification benchmarks with real genome and GTF annotations.

## What it does

For each target region (`chr:start-end`):

1. Extracts region genome FASTA + fully-contained transcripts from the input GTF
2. Assigns transcript abundances (`random`, `uniform`, or from external file)
3. Simulates paired-end reads with configurable strand specificity and gDNA contamination
4. Runs quantification with:
   - `hulkrna` (alignment + index + pipeline)
   - `salmon`
   - `kallisto`
5. Writes per-region and aggregate comparison reports

## Required external tools

- `minimap2`
- `samtools`
- `salmon`
- `kallisto`

These are checked automatically at startup.

## Example command

```bash
PYTHONPATH=src conda run -n hulkrna python scripts/benchmark_region_competition.py \
  --genome /path/to/genome.fa \
  --gtf /path/to/annotation.gtf.gz \
  --region-file scripts/regions.example.tsv \
  --outdir bench_regions \
  --n-fragments 30000 \
  --strand-specificity 0.85 \
  --gdna-abundance 20 \
  --abundance-mode random \
  --threads 4
```

## Key outputs

- `bench_regions/summary.json`
- `bench_regions/summary.csv`
- `bench_regions/summary.md`
- `bench_regions/diagnostics.md`
- `bench_regions/<region>/per_transcript_counts.csv`

## Notes

- Regions keep transcripts only when all exons are fully contained in the target window.
- gDNA truth counts are parsed from read names in the simulated FASTQ.
- `kallisto` uses `--rf-stranded` when `--strand-specificity >= 0.9`.
- For external abundances, use `--abundance-mode file --abundance-file abundances.tsv`
  where file has columns: `transcript_id,abundance`.
