#!/usr/bin/env python3
"""Generate a 1M-fragment oracle BAM for the PVT1/MYC locus and profile hulkrna.

Usage:
    cd /Users/mkiyer/proj/hulkrna
    PYTHONPATH=src python scripts/profile_pvt1_1m.py
"""
import cProfile
import csv
import pstats
import sys
import time
from dataclasses import dataclass
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths — reuse extracted region from existing benchmark
# ---------------------------------------------------------------------------
BASE = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC")
REGION_FA = BASE / "region.fa"
REGION_GTF = BASE / "region.gtf"
INDEX_DIR = BASE / "hulkrna_index"
ABUNDANCE_FILE = Path(
    "/Users/mkiyer/proj/hulkrna/scripts/benchmarking/"
    "chr8_pvt1_myc_abundance.tsv"
)

N_FRAGMENTS = 200_000
BAM_OUT = BASE / "profile_200k" / "reads_namesort.bam"


# ---------------------------------------------------------------------------
# Minimal genome wrapper (same as benchmark's StringGenome)
# ---------------------------------------------------------------------------
@dataclass
class StringGenome:
    seq: str
    name: str

    def __getitem__(self, key):
        return self.seq[key]

    def __len__(self):
        return len(self.seq)

    def fetch(self, chrom, start, end):
        return self.seq[start:end]

    def get_reference_length(self, chrom):
        return len(self.seq)


def load_region_fasta(path: Path) -> tuple[str, str]:
    """Return (name, sequence) from a single-contig FASTA."""
    lines = path.read_text().strip().split("\n")
    name = lines[0].lstrip(">").split()[0]
    seq = "".join(lines[1:]).upper()
    return name, seq


def parse_region_transcripts(gtf_path: Path, region_label: str):
    """Parse transcripts from the region GTF."""
    from hulkrna.transcript import Transcript
    from hulkrna.types import Interval, Strand

    tx_exons: dict[str, list[tuple[int, int]]] = {}
    tx_meta: dict[str, dict] = {}

    with open(gtf_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature != "exon":
                continue

            seqname = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand_str = parts[6]
            attrs_str = parts[8]

            # Parse attributes
            attrs = {}
            for token in attrs_str.split(";"):
                token = token.strip()
                if not token:
                    continue
                kv = token.split(None, 1)
                if len(kv) == 2:
                    attrs[kv[0]] = kv[1].strip().strip('"')

            tid = attrs.get("transcript_id", "")
            gid = attrs.get("gene_id", "")
            gname = attrs.get("gene_name", gid)
            gtype = attrs.get("gene_type", "protein_coding")

            tx_exons.setdefault(tid, []).append((start, end))
            tx_meta[tid] = {
                "g_id": gid,
                "g_name": gname,
                "g_type": gtype,
                "strand": Strand.from_str(strand_str),
            }

    transcripts = []
    for t_index, tid in enumerate(sorted(tx_exons)):
        meta = tx_meta[tid]
        exons = sorted(tx_exons[tid], key=lambda x: x[0])
        rel_exons = [Interval(s, e) for s, e in exons]
        t = Transcript(
            t_id=tid,
            t_index=t_index,
            g_id=meta["g_id"],
            g_name=meta["g_name"],
            g_type=meta["g_type"],
            ref=region_label,
            strand=meta["strand"],
            exons=rel_exons,
        )
        transcripts.append(t)

    return transcripts


def assign_abundances_from_file(transcripts, abundance_file: Path):
    """Load abundances from TSV file."""
    abund_map = {}
    sep = "\t" if str(abundance_file).endswith(".tsv") else ","
    with open(abundance_file) as fh:
        reader = csv.DictReader(fh, delimiter=sep)
        for row in reader:
            abund_map[row["transcript_id"]] = float(row["abundance"])

    for t in transcripts:
        t.abundance = max(1e-9, abund_map.get(t.t_id, 1.0))


def main():
    from hulkrna.index import TranscriptIndex
    from hulkrna.pipeline import run_pipeline
    from hulkrna.sim import OracleBamSimulator, SimConfig

    # --- Load region ---
    region_label, region_seq = load_region_fasta(REGION_FA)
    print(f"Region: {region_label} ({len(region_seq):,} bp)")

    # --- Parse transcripts ---
    transcripts = parse_region_transcripts(REGION_GTF, region_label)
    print(f"Transcripts: {len(transcripts)}")

    # --- Assign abundances ---
    assign_abundances_from_file(transcripts, ABUNDANCE_FILE)
    n_expressed = sum(1 for t in transcripts if t.abundance > 1.0)
    print(f"Expressed transcripts: {n_expressed}")

    # --- Generate 1M oracle BAM ---
    sim_cfg = SimConfig(
        frag_mean=250.0,
        frag_std=50.0,
        frag_min=50,
        frag_max=1000,
        read_length=150,
        error_rate=0.0,
        seed=42,
    )
    genome = StringGenome(seq=region_seq, name=region_label)

    BAM_OUT.parent.mkdir(parents=True, exist_ok=True)

    if BAM_OUT.exists():
        print(f"\nOracle BAM already exists: {BAM_OUT}")
        print("  Delete it to regenerate.")
    else:
        print(f"\nGenerating {N_FRAGMENTS:,}-fragment oracle BAM...")
        t0 = time.perf_counter()
        oracle_sim = OracleBamSimulator(
            genome, transcripts,
            config=sim_cfg,
            ref_name=region_label,
        )
        oracle_sim.write_bam(BAM_OUT, N_FRAGMENTS)
        t1 = time.perf_counter()
        print(f"BAM generated in {t1 - t0:.1f}s: {BAM_OUT}")

    # --- Load index ---
    index = TranscriptIndex.load(str(INDEX_DIR))
    print(f"\nIndex: {index.num_transcripts} transcripts, {index.num_genes} genes")

    # --- Profile pipeline ---
    print(f"\n{'='*70}")
    print(f"Profiling hulkrna pipeline on {N_FRAGMENTS:,} fragments...")
    print(f"{'='*70}\n")

    t0 = time.perf_counter()
    profiler = cProfile.Profile()
    profiler.enable()
    result = run_pipeline(str(BAM_OUT), index)
    profiler.disable()
    wall = time.perf_counter() - t0

    print(f"\n{'='*70}")
    print(f"Pipeline wall time: {wall:.2f}s")
    print(f"{'='*70}\n")

    # --- Stats ---
    stats = result.stats
    print(f"Fragments processed: {stats.total}")
    frags_per_sec = stats.total / wall if wall > 0 else 0
    print(f"Throughput: {frags_per_sec:,.0f} fragments/sec")
    print()

    # --- cProfile output ---
    ps = pstats.Stats(profiler)

    ps.sort_stats("cumulative")
    print("=== Top 50 by cumulative time ===")
    ps.print_stats(50)

    ps.sort_stats("tottime")
    print("\n=== Top 50 by self time ===")
    ps.print_stats(50)

    # --- Function call count ---
    total_calls = sum(v[0] for v in ps.stats.values())
    print(f"\nTotal function calls: {total_calls:,}")

    # Save profile
    prof_path = "/tmp/hulkrna_pvt1_200k_profile.prof"
    profiler.dump_stats(prof_path)
    print(f"Profile saved to {prof_path}")


if __name__ == "__main__":
    main()
