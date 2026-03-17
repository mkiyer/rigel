#!/usr/bin/env python
"""Generate a large calibration validation YAML config.

Creates a 2Mb genome with 60 expressed + 40 unexpressed transcripts
(100 total), producing ~400 annotation-partitioned regions.  Designed
to stress-test the purity-score-based gDNA calibration.

Usage:
    python scripts/gen_calibration_config.py > scripts/calibration_large.yaml
"""

import random
import yaml

GENOME = 2_000_000
SEED = 42
READ_LEN = 150
N_EXPRESSED = 60
N_UNEXPRESSED = 40
N_TOTAL = N_EXPRESSED + N_UNEXPRESSED

# Region spacing: ~18 kb per transcript unit → 100 × 18 kb ≈ 1.8 Mb
UNIT = GENOME // (N_TOTAL + 2)  # leave margin at both ends


def _make_exons(start, strand, rng):
    """Generate 2–4 exons within a ~3–8 kb span."""
    n_exons = rng.choice([2, 2, 3, 3, 4])
    exon_len_range = (200, 1200)
    intron_len_range = (800, 3000)

    exons = []
    pos = start
    for i in range(n_exons):
        elen = rng.randint(*exon_len_range)
        exons.append([pos, pos + elen])
        pos += elen
        if i < n_exons - 1:
            pos += rng.randint(*intron_len_range)
    return exons


def main():
    rng = random.Random(SEED)

    transcripts = {}
    nrnas = {}
    expressed_ids = []
    unexpressed_ids = []

    # Abundance profiles (exponential range for expressed transcripts)
    # Low: 10–50, Medium: 100–500, High: 1000–5000
    abundance_pool = (
        [rng.randint(10, 50) for _ in range(20)]
        + [rng.randint(100, 500) for _ in range(20)]
        + [rng.randint(1000, 5000) for _ in range(20)]
    )
    rng.shuffle(abundance_pool)

    # nRNA is ~30–60% of mRNA
    nrna_fracs = [rng.uniform(0.3, 0.6) for _ in range(N_EXPRESSED)]

    cursor = UNIT  # start with a gap
    strands = ["+", "-"]

    for i in range(N_TOTAL):
        idx = i + 1
        strand = strands[i % 2]
        exons = _make_exons(cursor, strand, rng)

        if i < N_EXPRESSED:
            tid = f"T{idx:02d}"
            expressed_ids.append(tid)
            # nRNA spans from first exon start to last exon end
            nrna_label = f"N{idx:02d}"
            nrnas[nrna_label] = [strand, exons[0][0], exons[-1][1]]
        else:
            tid = f"U{idx - N_EXPRESSED:02d}"
            unexpressed_ids.append(tid)

        transcripts[tid] = {"strand": strand, "exons": exons}

        # Advance cursor past the transcript span + intergenic gap
        span_end = exons[-1][1]
        gap = rng.randint(8000, 14000)
        cursor = span_end + gap

    # Build expression pattern
    pattern = {}
    for i, tid in enumerate(expressed_ids):
        ab = abundance_pool[i]
        pattern[tid] = ab
        nrna_label = f"N{i + 1:02d}"
        pattern[nrna_label] = int(ab * nrna_fracs[i])

    # Build sweep config
    sweep = {
        "strand_specificity": [0.5, 0.75, 0.9, 0.95],
        "n_rna_fragments": 50000,
        "gdna_fraction": [0.0, 0.1, 0.3, 0.5, 1.0],
        "gdna_frag_mean": [200, 350],
        "gdna_frag_std": 80,
        "gdna_strand_kappa": [5, 20, 100],
    }
    # Unexpressed transcripts always 0
    for uid in unexpressed_ids:
        sweep[uid] = 0

    config = {
        "genome_length": GENOME,
        "seed": SEED,
        "read_length": READ_LEN,
        "rna": {
            "frag_mean": 250,
            "frag_std": 50,
            "frag_min": 80,
            "frag_max": 600,
        },
        "gdna": {
            "frag_mean": 350,
            "frag_std": 80,
            "frag_min": 100,
            "frag_max": 1000,
        },
        "transcripts": transcripts,
        "nrnas": nrnas,
        "sweep": sweep,
        "patterns": [pattern],
    }

    header = (
        "# calibration_large.yaml — Large calibration validation config\n"
        "#\n"
        f"# Genome: {GENOME / 1e6:.0f} Mb, {N_EXPRESSED} expressed + "
        f"{N_UNEXPRESSED} unexpressed transcripts\n"
        "# Grid: 1 pattern × 5 gdna_frac × 4 SS × 3 kappa × 2 gdna_fl "
        "= 120 runs\n"
        "# (minus degenerate: ~116 effective runs)\n"
        "#\n"
        "# Validates:\n"
        "#   - Purity-score seed selection across region types\n"
        "#   - Frozen density estimation\n"
        "#   - κ recovery at various overdispersion levels\n"
        "#   - FL model recovery with similar vs different gDNA FL\n"
        "#   - Robustness from unstranded (0.5) to stranded (0.95)\n"
        "#\n"
    )

    print(header)
    yaml.dump(config, __import__("sys").stdout,
              default_flow_style=None, sort_keys=False, width=120)


if __name__ == "__main__":
    main()
