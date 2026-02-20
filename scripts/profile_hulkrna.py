#!/usr/bin/env python3
"""Performance profiling benchmark for hulkrna.

Simulates reads in a complex region (many overlapping transcripts), aligns
with minimap2, then profiles hulkrna's pipeline to identify bottlenecks.

Usage:
    # Quick run (default region = EGFR):
    PYTHONPATH=src conda run -n hulkrna python scripts/profile_hulkrna.py

    # Custom region:
    PYTHONPATH=src conda run -n hulkrna python scripts/profile_hulkrna.py \
        --region chr7:55019017-55211628 --n-fragments 100000

    # Profile all benchmark regions:
    PYTHONPATH=src conda run -n hulkrna python scripts/profile_hulkrna.py \
        --region-file scripts/regions_benchmark.tsv --n-fragments 100000

Output per region:
    <outdir>/<region>/
        profile_stats.txt    – cProfile text summary (top 60 by cumtime)
        profile_stats.prof   – binary cProfile output (for snakeviz/pstats)
        profile_results.json – structured timing breakdown
        pipeline_breakdown.md – human-readable report
"""
from __future__ import annotations

import argparse
import cProfile
import io
import json
import logging
import os
import pstats
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

# Add src to path if needed
src = Path(__file__).resolve().parent.parent / "src"
if str(src) not in sys.path:
    sys.path.insert(0, str(src))

from hulkrna.gtf import GTF
from hulkrna.index import HulkIndex, write_bed12
from hulkrna.pipeline import run_pipeline, scan_and_buffer, count_from_buffer
from hulkrna.sim import GDNAConfig, ReadSimulator, SimConfig, reverse_complement
from hulkrna.transcript import Transcript
from hulkrna.types import Interval, Strand

logger = logging.getLogger(__name__)

# ── Defaults ─────────────────────────────────────────────────────────

DEFAULT_GENOME = Path(
    "/Users/mkiyer/Downloads/hulkrna_runs/refs/human/genome_controls.fasta.bgz"
)
DEFAULT_GTF = Path(
    "/Users/mkiyer/Downloads/hulkrna_runs/refs/human/genes_controls.gtf.gz"
)
DEFAULT_OUTDIR = Path("/Users/mkiyer/Downloads/hulkrna_runs/profiling")
DEFAULT_REGION = "chr7:55019017-55211628"  # EGFR
DEFAULT_N_FRAGMENTS = 100_000
DEFAULT_GDNA_ABUNDANCE_PCTL = 20.0  # moderate gDNA
DEFAULT_STRAND_SPECIFICITY = 0.95
DEFAULT_SEED = 42


# ── Region extraction (reused from benchmark_region_competition) ─────

import gzip
import re
from collections import Counter


def open_textmaybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_region(region_str: str):
    m = re.match(r"(.+):(\d+)-(\d+)$", region_str)
    if not m:
        raise ValueError(f"Invalid region format: {region_str}")
    chrom = m.group(1)
    start1 = int(m.group(2))
    end1 = int(m.group(3))
    label = f"{chrom}_{start1}_{end1}"
    return label, chrom, start1 - 1, end1


def parse_regions(args):
    regions = []
    for r in (args.region or []):
        label, chrom, s0, e0 = parse_region(r)
        regions.append((label, chrom, s0, e0))
    if args.region_file:
        with open(args.region_file) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                chrom = parts[0]
                s1 = int(parts[1])
                e1 = int(parts[2])
                label = parts[3] if len(parts) > 3 else f"{chrom}_{s1}_{e1}"
                regions.append((label, chrom, s1 - 1, e1))
    if not regions:
        # Default to EGFR
        label, chrom, s0, e0 = parse_region(DEFAULT_REGION)
        regions.append((label, chrom, s0, e0))
    return regions


def extract_region(genome_fa, gtf_path, chrom, start0, end0, label, out_dir):
    """Extract region FASTA + GTF + transcripts."""
    out_dir.mkdir(parents=True, exist_ok=True)
    out_fa = out_dir / "region.fa"
    out_gtf = out_dir / "region.gtf"

    # Extract FASTA
    result = subprocess.run(
        ["samtools", "faidx", str(genome_fa), f"{chrom}:{start0+1}-{end0}"],
        capture_output=True, text=True, check=True,
    )
    fasta_lines = result.stdout.strip().split("\n")
    region_seq = "".join(fasta_lines[1:]).upper()
    with open(out_fa, "w") as fh:
        fh.write(f">{label}\n")
        for i in range(0, len(region_seq), 80):
            fh.write(region_seq[i:i+80] + "\n")
    subprocess.run(["samtools", "faidx", str(out_fa)], check=True, capture_output=True)

    # Parse GTF
    tx_exons = defaultdict(list)
    tx_meta = {}
    with open_textmaybe_gzip(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[0] != chrom or parts[2] != "exon":
                continue
            s0_exon = int(parts[3]) - 1
            e0_exon = int(parts[4])
            strand_str = parts[6]
            attrs = {}
            for a in parts[8].split(";"):
                a = a.strip()
                if not a:
                    continue
                m = re.match(r'(\w+)\s+"([^"]*)"', a)
                if m:
                    attrs[m.group(1)] = m.group(2)
            tid = attrs.get("transcript_id", "")
            gid = attrs.get("gene_id", "")
            if not tid or not gid:
                continue
            tx_exons[tid].append((s0_exon, e0_exon))
            tx_meta[tid] = {
                "g_id": gid,
                "g_name": attrs.get("gene_name", gid),
                "g_type": attrs.get("gene_type", attrs.get("gene_biotype", "unknown")),
                "strand": Strand.from_str(strand_str),
            }

    # Filter to fully-contained transcripts
    transcripts = []
    t_index = 0
    for tid in sorted(tx_exons):
        exons = tx_exons[tid]
        if not all(s >= start0 and e <= end0 for s, e in exons):
            continue
        meta = tx_meta[tid]
        rel_exons = sorted(
            [Interval(s - start0, e - start0) for s, e in exons],
            key=lambda x: x.start,
        )
        t = Transcript(
            t_id=tid, t_index=t_index, g_id=meta["g_id"],
            g_name=meta["g_name"], g_type=meta["g_type"],
            ref=label, strand=meta["strand"], exons=rel_exons,
        )
        transcripts.append(t)
        t_index += 1

    # Write GTF
    with open(out_gtf, "w") as out:
        for t in transcripts:
            for exon in t.exons:
                gtf_obj = GTF(
                    seqname=label, source="hulkrna_profile",
                    feature="exon", start=exon.start, end=exon.end,
                    score=1.0, strand=Strand.to_str(t.strand), phase=".",
                    attrs={"gene_id": t.g_id, "transcript_id": t.t_id,
                           "gene_name": t.g_name, "gene_type": t.g_type},
                    tags={"basic"},
                )
                out.write(str(gtf_obj) + "\n")

    return out_fa, out_gtf, transcripts, region_seq


class StringGenome:
    def __init__(self, seq, name):
        self.seq = seq
        self.name = name
    def __getitem__(self, key):
        return self.seq[key]
    def __len__(self):
        return len(self.seq)
    def fetch(self, chrom, start, end):
        assert chrom == self.name
        return self.seq[start:end]
    def get_reference_length(self, chrom):
        assert chrom == self.name
        return len(self.seq)
    def references(self):
        return [self.name]
    @property
    def lengths(self):
        return {self.name: len(self.seq)}


def assign_abundances(transcripts, rng, min_ab=0.01, max_ab=100000.0):
    """Log-uniform abundance assignment."""
    log_ab = rng.uniform(np.log(max(min_ab, 1e-9)), np.log(max(max_ab, 1.0)),
                         size=len(transcripts))
    raw = np.exp(log_ab)
    med = float(np.median(raw))
    if med > 0:
        raw *= 100.0 / med
    for t, ab in zip(transcripts, raw):
        t.abundance = max(1e-9, float(ab))


def write_transcript_fasta(transcripts, region_seq, out_fa):
    with open(out_fa, "w") as out:
        for t in transcripts:
            parts = [region_seq[e.start:e.end] for e in t.exons]
            seq = "".join(parts)
            if t.strand == Strand.NEG:
                seq = reverse_complement(seq)
            out.write(f">{t.t_id}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")


def align_minimap2(region_fa, fq_r1, fq_r2, out_bam, bed_path=None):
    """Align PE reads with minimap2 → name-sorted BAM."""
    cmd = ["minimap2", "-ax", "splice:sr", "--secondary=yes"]
    if bed_path:
        cmd.extend(["-j", str(bed_path)])
    cmd.extend([str(region_fa), str(fq_r1), str(fq_r2)])
    sort_cmd = ["samtools", "sort", "-n", "-o", str(out_bam), "-"]
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(sort_cmd, stdin=p1.stdout, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    p1.stdout.close()
    p2.communicate()
    p1.wait()


# ── Profiling ────────────────────────────────────────────────────────

def profile_pipeline(bam_path, index, seed, region_label):
    """Run the full pipeline under cProfile and return (profiler, result, wall_time)."""
    import pysam

    profiler = cProfile.Profile()

    t0 = time.perf_counter()
    profiler.enable()
    pipe = run_pipeline(
        bam_path, index, seed=seed,
        sj_strand_tag="auto",
        include_multimap=True,
    )
    profiler.disable()
    wall_time = time.perf_counter() - t0

    return profiler, pipe, wall_time


def profile_stages_separately(bam_path, index, seed):
    """Profile scan_and_buffer and count_from_buffer separately for breakdown."""
    import pysam

    timings = {}

    # Stage 1: scan_and_buffer
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    t0 = time.perf_counter()
    stats, strand_models, insert_models, buf = scan_and_buffer(
        bam.fetch(until_eof=True), index,
        include_multimap=True,
        sj_strand_tag="auto",
    )
    timings["scan_and_buffer"] = time.perf_counter() - t0
    bam.close()

    timings["n_fragments"] = stats.n_fragments
    timings["n_buffered"] = buf.total_fragments
    timings["buffer_mb"] = buf.memory_bytes / 1024**2

    # Finalize models (cache derived values for fast scoring)
    strand_models.finalize()
    insert_models.finalize()

    # Stage 2: count_from_buffer
    t0 = time.perf_counter()
    counter = count_from_buffer(
        buf, index, strand_models, insert_models, stats,
        seed=seed,
    )
    timings["count_from_buffer"] = time.perf_counter() - t0

    timings["total"] = timings["scan_and_buffer"] + timings["count_from_buffer"]

    return timings, stats


def profile_run_pipeline_cprofile(bam_path, index, seed):
    """Full pipeline under cProfile with detailed function breakdown."""
    profiler = cProfile.Profile()

    t0 = time.perf_counter()
    profiler.enable()
    pipe = run_pipeline(
        bam_path, index, seed=seed,
        sj_strand_tag="auto",
        include_multimap=True,
    )
    profiler.disable()
    wall_time = time.perf_counter() - t0

    return profiler, pipe, wall_time


def extract_top_functions(profiler, n=60):
    """Extract top N functions by cumulative time from cProfile."""
    stream = io.StringIO()
    ps = pstats.Stats(profiler, stream=stream)
    ps.sort_stats("cumulative")
    ps.print_stats(n)
    return stream.getvalue()


def extract_callers(profiler, n=30):
    """Extract caller info for hot functions."""
    stream = io.StringIO()
    ps = pstats.Stats(profiler, stream=stream)
    ps.sort_stats("cumulative")
    ps.print_callers(n)
    return stream.getvalue()


def write_report(region_label, n_fragments, n_transcripts, n_genes,
                 timings, profile_text, callers_text, out_dir):
    """Write pipeline_breakdown.md."""
    lines = [
        f"# hulkrna Performance Profile: {region_label}",
        "",
        "## Setup",
        "",
        f"- **Fragments:** {n_fragments:,}",
        f"- **Transcripts:** {n_transcripts}",
        f"- **Genes:** {n_genes}",
        f"- **Buffered:** {timings.get('n_buffered', '?'):,}",
        f"- **Buffer memory:** {timings.get('buffer_mb', 0):.1f} MB",
        "",
        "## Stage Breakdown",
        "",
        f"| Stage | Time (s) | % of Total |",
        f"| --- | ---: | ---: |",
    ]
    total = timings.get("total", 1)
    for stage in ["scan_and_buffer", "count_from_buffer"]:
        t = timings.get(stage, 0)
        pct = t / total * 100 if total > 0 else 0
        lines.append(f"| {stage} | {t:.3f} | {pct:.1f}% |")
    lines.append(f"| **Total** | **{total:.3f}** | **100%** |")
    lines.append("")

    # Throughput
    nf = timings.get("n_fragments", n_fragments)
    if total > 0:
        lines.append("## Throughput")
        lines.append("")
        lines.append(f"- **{nf / total:,.0f} fragments/sec** (full pipeline)")
        lines.append(f"- **{nf / timings.get('scan_and_buffer', 1):,.0f} fragments/sec** (scan_and_buffer only)")
        lines.append(f"- **{nf / timings.get('count_from_buffer', 1):,.0f} fragments/sec** (count_from_buffer only)")
        lines.append("")

    lines.append("## cProfile: Top Functions by Cumulative Time")
    lines.append("")
    lines.append("```")
    lines.append(profile_text)
    lines.append("```")
    lines.append("")

    lines.append("## cProfile: Callers")
    lines.append("")
    lines.append("```")
    lines.append(callers_text)
    lines.append("```")

    report_path = out_dir / "pipeline_breakdown.md"
    report_path.write_text("\n".join(lines) + "\n")
    return report_path


# ── Main ─────────────────────────────────────────────────────────────

def run_region_profile(label, chrom, start0, end0, args, out_root):
    """Profile hulkrna on one region."""
    reg_dir = out_root / label
    reg_dir.mkdir(parents=True, exist_ok=True)

    logger.info("="*60)
    logger.info("Profiling region: %s (%s:%d-%d)", label, chrom, start0+1, end0)
    logger.info("="*60)

    # Extract region
    region_fa, region_gtf, transcripts, region_seq = extract_region(
        args.genome, args.gtf, chrom, start0, end0, label, reg_dir,
    )
    if not transcripts:
        logger.warning("No transcripts in region %s, skipping", label)
        return None

    n_transcripts = len(transcripts)
    n_genes = len({t.g_id for t in transcripts})
    logger.info("Region %s: %d transcripts, %d genes", label, n_transcripts, n_genes)

    # Assign abundances
    rng = np.random.default_rng(args.seed + abs(hash(label)) % 10000)
    assign_abundances(transcripts, rng, min_ab=args.abundance_min, max_ab=args.abundance_max)

    # Rewrite GTF with abundances
    with open(region_gtf, "w") as out:
        for t in transcripts:
            for exon in t.exons:
                gtf_obj = GTF(
                    seqname=label, source="hulkrna_profile",
                    feature="exon", start=exon.start, end=exon.end,
                    score=t.abundance, strand=Strand.to_str(t.strand), phase=".",
                    attrs={"gene_id": t.g_id, "transcript_id": t.t_id,
                           "gene_name": t.g_name, "gene_type": t.g_type},
                    tags={"basic"},
                )
                out.write(str(gtf_obj) + "\n")

    # Write transcript FASTA
    tx_fa = reg_dir / "transcripts.fa"
    write_transcript_fasta(transcripts, region_seq, tx_fa)

    # BED12 for minimap2
    bed_path = reg_dir / "annotation.bed"
    write_bed12(transcripts, bed_path)

    # Build hulkrna index
    index_dir = reg_dir / "hulkrna_index"
    HulkIndex.build(region_fa, region_gtf, index_dir, write_tsv=False)
    index = HulkIndex.load(index_dir)

    # Compute gDNA abundance (moderate level = 20th percentile)
    abundances = np.array([t.abundance for t in transcripts])
    gdna_abundance = float(np.percentile(abundances, args.gdna_percentile))
    logger.info("gDNA abundance (pctl %.0f): %.2f", args.gdna_percentile, gdna_abundance)

    # Simulate reads
    sim_cfg = SimConfig(
        frag_mean=250.0, frag_std=50.0, frag_min=50, frag_max=1000,
        read_length=150, error_rate=0.0,
        strand_specificity=args.strand_specificity,
        seed=args.seed,
    )
    gdna_cfg = GDNAConfig(
        abundance=gdna_abundance,
        frag_mean=350.0, frag_std=100.0, frag_min=100, frag_max=1000,
    ) if gdna_abundance > 0 else None

    sim_genome = StringGenome(region_seq, label)
    simulator = ReadSimulator(sim_genome, transcripts, config=sim_cfg, gdna_config=gdna_cfg)
    fastq_dir = reg_dir / "reads"
    logger.info("Simulating %d fragments...", args.n_fragments)
    fq_r1, fq_r2 = simulator.write_fastq(fastq_dir, args.n_fragments, prefix="sim")

    # Align
    bam_ns = reg_dir / "align" / "reads_namesort.bam"
    bam_ns.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Aligning with minimap2...")
    t0 = time.perf_counter()
    align_minimap2(region_fa, fq_r1, fq_r2, bam_ns, bed_path=bed_path)
    align_time = time.perf_counter() - t0
    logger.info("Alignment: %.2fs", align_time)

    # ── Profile: stage-level breakdown ──────────────────────────────

    logger.info("Profiling stage-level breakdown...")
    timings, stats = profile_stages_separately(bam_ns, index, args.seed)
    logger.info("scan_and_buffer: %.3fs", timings["scan_and_buffer"])
    logger.info("count_from_buffer: %.3fs", timings["count_from_buffer"])
    logger.info("Total: %.3fs", timings["total"])

    # ── Profile: cProfile function-level ────────────────────────────

    logger.info("Profiling full pipeline with cProfile...")
    profiler, pipe, wall_time = profile_run_pipeline_cprofile(bam_ns, index, args.seed)
    logger.info("cProfile wall time: %.3fs", wall_time)

    # Save profile data
    prof_path = reg_dir / "profile_stats.prof"
    profiler.dump_stats(str(prof_path))

    profile_text = extract_top_functions(profiler, n=60)
    callers_text = extract_callers(profiler, n=30)

    prof_txt_path = reg_dir / "profile_stats.txt"
    prof_txt_path.write_text(profile_text)

    callers_path = reg_dir / "profile_callers.txt"
    callers_path.write_text(callers_text)

    # Save timings JSON
    results = {
        "region": label,
        "n_transcripts": n_transcripts,
        "n_genes": n_genes,
        "n_fragments": args.n_fragments,
        "n_buffered": timings.get("n_buffered", 0),
        "gdna_pctl": args.gdna_percentile,
        "gdna_abundance": gdna_abundance,
        "strand_specificity": args.strand_specificity,
        "align_time_sec": align_time,
        "scan_and_buffer_sec": timings["scan_and_buffer"],
        "count_from_buffer_sec": timings["count_from_buffer"],
        "total_pipeline_sec": timings["total"],
        "cprofile_wall_sec": wall_time,
        "throughput_frags_per_sec": args.n_fragments / timings["total"] if timings["total"] > 0 else 0,
        "buffer_mb": timings.get("buffer_mb", 0),
    }
    json_path = reg_dir / "profile_results.json"
    json_path.write_text(json.dumps(results, indent=2))

    # Write report
    report_path = write_report(
        label, args.n_fragments, n_transcripts, n_genes,
        timings, profile_text, callers_text, reg_dir,
    )
    logger.info("Report: %s", report_path)

    return results


def build_parser():
    p = argparse.ArgumentParser(
        description="Profile hulkrna pipeline performance",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--genome", type=Path, default=DEFAULT_GENOME)
    p.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    p.add_argument("--region", action="append", default=None)
    p.add_argument("--region-file", type=Path, default=None)
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    p.add_argument("--n-fragments", type=int, default=DEFAULT_N_FRAGMENTS)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p.add_argument("--strand-specificity", type=float, default=DEFAULT_STRAND_SPECIFICITY)
    p.add_argument("--gdna-percentile", type=float, default=DEFAULT_GDNA_ABUNDANCE_PCTL)
    p.add_argument("--abundance-min", type=float, default=0.01)
    p.add_argument("--abundance-max", type=float, default=100000.0)
    p.add_argument("--verbose", action="store_true")
    return p


def main():
    parser = build_parser()
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    regions = parse_regions(args)
    args.outdir.mkdir(parents=True, exist_ok=True)

    all_results = []
    for label, chrom, s0, e0 in regions:
        result = run_region_profile(label, chrom, s0, e0, args, args.outdir)
        if result:
            all_results.append(result)

    # Write combined summary
    if all_results:
        summary_path = args.outdir / "profile_summary.json"
        summary_path.write_text(json.dumps(all_results, indent=2))

        # Print summary table
        print("\n" + "="*72)
        print("PERFORMANCE SUMMARY")
        print("="*72)
        print(f"{'Region':<30} {'Frags':>8} {'Tx':>4} {'Genes':>5} {'Scan(s)':>8} {'Count(s)':>9} {'Total(s)':>9} {'frags/s':>10}")
        print("-"*72)
        for r in all_results:
            print(
                f"{r['region']:<30} "
                f"{r['n_fragments']:>8,} "
                f"{r['n_transcripts']:>4} "
                f"{r['n_genes']:>5} "
                f"{r['scan_and_buffer_sec']:>8.2f} "
                f"{r['count_from_buffer_sec']:>9.2f} "
                f"{r['total_pipeline_sec']:>9.2f} "
                f"{r['throughput_frags_per_sec']:>10,.0f}"
            )
        print("="*72)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
