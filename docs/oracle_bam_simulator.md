# Oracle BAM Simulator

The **Oracle BAM Simulator** generates a name-sorted BAM file with
perfectly aligned paired-end reads, completely bypassing the
FASTQ → alignment → BAM pipeline.  This isolates hulkrna's
abundance-estimation errors from alignment errors in benchmarking.

## Motivation

When benchmarking hulkrna against salmon/kallisto, alignment errors
(misalignments, missed spliced reads, paralog cross-mapping) confound
the comparison.  The oracle BAM simulator lets you measure hulkrna's
**pure estimation accuracy** — how well the EM algorithm distributes
fragments among isoforms — without any alignment noise.

## Architecture

`OracleBamSimulator` wraps the existing `ReadSimulator` and reuses its
fragment-selection logic:

- Pool splitting (mature RNA / nascent RNA / gDNA)
- Abundance-weighted transcript selection
- Fragment-length sampling from truncated normal distribution
- Strand-specificity simulation (R1↔R2 swaps)

Instead of writing FASTQ, it projects each fragment's transcript-space
coordinates back to genomic coordinates and writes BAM records with:

| Fragment type | Coordinate mapping | CIGAR |
|---|---|---|
| **Mature RNA** | Transcript position → exon blocks | `M...N...M` (spliced) |
| **Nascent RNA** | Pre-mRNA position → genomic span | `M` (contiguous) |
| **gDNA** | Already genomic | `M` (contiguous) |

### BAM conventions

- **Sorting**: name-sorted (`SO:queryname`) for `run_pipeline`
- **Paired-end flags**: `0x1`, `0x2`, `0x40`/`0x80`, `0x10`/`0x20`
- **Mapping quality**: 255 (perfect)
- **`NH` tag**: 1 (unique mapping)
- **`XS` tag**: splice-junction strand for spliced reads (STAR/minimap2 convention)
- **Read names**: same format as `ReadSimulator` — ground truth encoded

## Usage

### Python API

```python
from hulkrna.sim import OracleBamSimulator, SimConfig, GDNAConfig
from hulkrna.sim.genome import MutableGenome

# Set up genome and transcripts (same as ReadSimulator)
genome = MutableGenome(length=10000, seed=42)
# ... build transcripts with GeneBuilder ...

config = SimConfig(
    frag_mean=250, frag_std=50,
    read_length=150,
    strand_specificity=0.95,
    seed=42,
)
gdna_config = GDNAConfig(abundance=10.0)

sim = OracleBamSimulator(
    genome, transcripts,
    config=config,
    gdna_config=gdna_config,
    ref_name="chr1",
)

# Write name-sorted BAM (default)
sim.write_bam("output.bam", n_fragments=50000)

# Write coordinate-sorted BAM with index
sim.write_bam("output.bam", n_fragments=50000, name_sorted=False)
```

### Benchmark integration

The oracle simulator is integrated into `benchmark_region_competition.py`
as the `--aligner oracle` option:

```bash
# Run benchmark with oracle BAM (no alignment)
python scripts/benchmarking/benchmark_region_competition.py \
    --genome /path/to/genome.fa \
    --gtf /path/to/annotation.gtf \
    --region chr11:5225464-5250625 \
    --aligner oracle \
    --n-fragments 50000 \
    --output-dir /path/to/output

# Compare with minimap2 alignment
python scripts/benchmarking/benchmark_region_competition.py \
    --genome /path/to/genome.fa \
    --gtf /path/to/annotation.gtf \
    --region chr11:5225464-5250625 \
    --aligner minimap2 \
    --n-fragments 50000 \
    --output-dir /path/to/output
```

When `--aligner oracle` is used:
1. `OracleBamSimulator` writes a perfectly-aligned BAM for hulkrna
2. A separate `ReadSimulator` (same seed) writes FASTQ for
   salmon/kallisto comparison
3. Truth counts are parsed from BAM read names
4. All tools are scored against identical ground truth

This lets you run benchmarks with alignment (`minimap2`, `hisat2`)
or without (`oracle`) and compare how much of each tool's error
budget comes from alignment vs. estimation.

## Coordinate projection details

### Mature RNA (spliced transcripts)

Fragment positions from `ReadSimulator` are in **transcript space**
(exon-concatenated, 5′→3′).  For multi-exon transcripts, the
simulator maps these back through the exon coordinate table:

```
Transcript space:   [0 ........... t_len]
                       exon1   exon2   exon3
Genomic space:      [--===---====---===--]
                         ↑         ↑
                         N         N   (intron skips in CIGAR)
```

For negative-strand transcripts, transcript coordinates are mirrored
(position 0 = rightmost genomic exon) before projection.

### Nascent RNA (pre-mRNA)

Pre-mRNA fragments span the full genomic extent including introns.
The `ReadSimulator` reverse-complements NEG-strand pre-mRNA, so
coordinates are mirrored before converting to genomic space.  The
resulting alignment is always a simple contiguous `M` block.

### gDNA

gDNA fragments are already in genomic coordinates.  The strand is
chosen randomly (50/50), and FR library convention determines R1/R2
orientation.

## Read orientation (FR library)

The simulator follows the standard Illumina dUTP (FR) convention:

| Transcript strand | R1 maps to | R2 maps to |
|---|---|---|
| + | − (reverse) | + (forward) |
| − | + (forward) | − (reverse) |

When `strand_specificity < 1.0`, R1↔R2 are swapped with probability
`1 − strand_specificity` to simulate imperfect dUTP incorporation.
