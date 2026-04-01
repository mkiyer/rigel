#!/usr/bin/env python3
"""Verify gDNA fragment generation in the whole-genome simulator.

Creates a small synthetic genome + annotation, runs WholeGenomeSimulator
with known gDNA configs, and validates:

1. Fragment lengths match the configured distribution (mean, std, bounds)
2. NH=1 tag is present on all gDNA BAM records
3. Read-name encoded coordinates match BAM TLEN
4. CIGAR is a single M block with correct length
5. Manifest stores the gDNA config
6. Multiple gDNA configs produce distinct distributions

Usage:
    conda activate rigel
    python scripts/debug/verify_gdna_sim.py
"""
from __future__ import annotations

import gzip
import json
import sys
import tempfile
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pysam

# Ensure the project root is importable
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from rigel.transcript import Transcript
from rigel.types import Interval, Strand
from scripts.sim.sim import (
    GDNASimConfig,
    SimConfig,
    SimulationParams,
    WholeGenomeSimulator,
    write_manifest,
)

# ── Helpers ──────────────────────────────────────────────────────────


def _make_test_fasta(path: Path, chrom_name: str = "chr1", length: int = 50_000):
    """Write a small random FASTA + index for testing."""
    rng = np.random.default_rng(42)
    seq = "".join(rng.choice(list("ACGT"), size=length))
    with open(path, "w") as f:
        f.write(f">{chrom_name}\n")
        for i in range(0, length, 80):
            f.write(seq[i : i + 80] + "\n")
    pysam.faidx(str(path))
    return chrom_name, length


def _make_test_transcripts(chrom: str) -> list[Transcript]:
    """Build a small set of transcripts on the test chromosome."""
    txs = []
    # Transcript at start of chrom (single exon, 500 bp)
    t1 = Transcript(
        ref=chrom, strand=Strand.POS,
        exons=[Interval(1000, 1500)],
        t_id="TX1", g_id="G1", g_name="Gene1",
        abundance=100.0, nrna_abundance=0.0,
    )
    t1.compute_length()
    t1.t_index = 0
    txs.append(t1)

    # Transcript with two exons (spliced, 400 bp)
    t2 = Transcript(
        ref=chrom, strand=Strand.NEG,
        exons=[Interval(5000, 5200), Interval(6000, 6200)],
        t_id="TX2", g_id="G2", g_name="Gene2",
        abundance=200.0, nrna_abundance=0.0,
    )
    t2.compute_length()
    t2.t_index = 1
    txs.append(t2)
    return txs


def _parse_gdna_bam(bam_path: str) -> list[dict]:
    """Parse gDNA R1 records from an oracle BAM."""
    records = []
    with pysam.AlignmentFile(bam_path, "rb") as f:
        for read in f:
            if not read.query_name.startswith("gdna:"):
                continue
            if not read.is_read1:
                continue
            tags = dict(read.tags)

            # Parse coordinates from read name: gdna:chr:start-end:strand:idx
            parts = read.query_name.split(":")
            name_start, name_end = parts[2].split("-")
            frag_len_from_name = int(name_end) - int(name_start)

            records.append({
                "qname": read.query_name,
                "frag_len_from_name": frag_len_from_name,
                "tlen": abs(read.template_length),
                "cigar": read.cigartuples,
                "read_len": read.query_length,
                "nh": tags.get("NH"),
                "ref_start": read.reference_start,
            })
    return records


# ── Test functions ───────────────────────────────────────────────────


def test_gdna_fragment_distribution(
    fasta_path: Path,
    transcripts: list[Transcript],
    gdna_cfg: GDNASimConfig,
    n_gdna: int,
    label: str,
) -> bool:
    """Run simulation with given gDNA config and validate fragment lengths."""
    sim_params = SimulationParams(
        n_rna_fragments=100,
        sim_seed=42,
        frag_mean=250.0, frag_std=50.0, frag_min=50, frag_max=1000,
        read_length=101, error_rate=0.0,
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        outdir = Path(tmpdir) / "out"
        sim = WholeGenomeSimulator(
            fasta_path, transcripts, sim_params, gdna_cfg,
            strand_specificity=0.9, seed=99,
        )
        _, _, bam_path = sim.simulate_and_write(outdir, n_rna=100, n_gdna=n_gdna)
        sim.close()

        records = _parse_gdna_bam(str(bam_path))

    if len(records) == 0:
        print(f"  [{label}] FAIL: No gDNA records found in BAM")
        return False

    frag_lens = np.array([r["frag_len_from_name"] for r in records])
    tlens = np.array([r["tlen"] for r in records])

    ok = True

    # Check 1: Fragment lengths within bounds
    if frag_lens.min() < gdna_cfg.frag_min:
        print(f"  [{label}] FAIL: min frag_len={frag_lens.min()} < frag_min={gdna_cfg.frag_min}")
        ok = False
    if frag_lens.max() > gdna_cfg.frag_max:
        print(f"  [{label}] FAIL: max frag_len={frag_lens.max()} > frag_max={gdna_cfg.frag_max}")
        ok = False

    # Check 2: Mean is close to configured mean (within 3 sigma / sqrt(n))
    expected_std_of_mean = gdna_cfg.frag_std / np.sqrt(len(frag_lens))
    mean_tolerance = max(3 * expected_std_of_mean, 10.0)  # at least 10 for small samples
    actual_mean = frag_lens.mean()
    if abs(actual_mean - gdna_cfg.frag_mean) > mean_tolerance:
        print(
            f"  [{label}] FAIL: mean frag_len={actual_mean:.1f} "
            f"vs expected ~{gdna_cfg.frag_mean} (tolerance={mean_tolerance:.1f})"
        )
        ok = False

    # Check 3: TLEN matches fragment length from read name
    tlen_mismatches = np.sum(tlens != frag_lens)
    if tlen_mismatches > 0:
        print(f"  [{label}] FAIL: {tlen_mismatches}/{len(records)} TLEN != frag_len from name")
        ok = False

    # Check 4: NH tag is present and equals 1
    nh_missing = sum(1 for r in records if r["nh"] is None)
    nh_wrong = sum(1 for r in records if r["nh"] is not None and r["nh"] != 1)
    if nh_missing > 0:
        print(f"  [{label}] FAIL: {nh_missing}/{len(records)} gDNA records missing NH tag")
        ok = False
    if nh_wrong > 0:
        print(f"  [{label}] FAIL: {nh_wrong}/{len(records)} gDNA records have NH != 1")
        ok = False

    # Check 5: CIGAR is a single M block
    for r in records:
        cigar = r["cigar"]
        if len(cigar) != 1 or cigar[0][0] != pysam.CMATCH:
            print(f"  [{label}] FAIL: unexpected CIGAR {cigar} for {r['qname']}")
            ok = False
            break

    # Check 6: Read length = min(read_len_cfg, frag_len)
    for r in records:
        expected_read_len = min(sim_params.read_length, r["frag_len_from_name"])
        if r["read_len"] != expected_read_len:
            print(
                f"  [{label}] FAIL: read_len={r['read_len']} "
                f"vs expected min({sim_params.read_length}, {r['frag_len_from_name']})={expected_read_len}"
            )
            ok = False
            break

    if ok:
        print(
            f"  [{label}] PASS: {len(records)} gDNA fragments, "
            f"mean={actual_mean:.1f} (expected ~{gdna_cfg.frag_mean}), "
            f"min={frag_lens.min()}, max={frag_lens.max()}, "
            f"all NH=1, CIGAR/TLEN correct"
        )
    return ok


def test_manifest_stores_gdna_config() -> bool:
    """Verify that write_manifest includes the gDNA config."""
    gdna_cfg = GDNASimConfig(
        frag_mean=400.0, frag_std=80.0, frag_min=120, frag_max=900,
    )
    cfg = SimConfig(
        genome="/fake/genome.fa",
        gtf="/fake/genes.gtf",
        gdna=gdna_cfg,
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        outdir = Path(tmpdir)
        write_manifest(outdir, cfg, [])

        with open(outdir / "manifest.json") as f:
            manifest = json.load(f)

    ok = True

    # Check gDNA section exists
    if "gdna" not in manifest:
        print("  [manifest] FAIL: 'gdna' section missing from manifest")
        return False

    gdna = manifest["gdna"]
    expected = asdict(gdna_cfg)
    for key, val in expected.items():
        if key not in gdna:
            print(f"  [manifest] FAIL: key '{key}' missing from manifest.gdna")
            ok = False
        elif gdna[key] != val:
            print(f"  [manifest] FAIL: gdna.{key}={gdna[key]} != expected {val}")
            ok = False

    # Check nrna and strand_specificities are also stored
    if "nrna" not in manifest:
        print("  [manifest] FAIL: 'nrna' section missing from manifest")
        ok = False
    if "strand_specificities" not in manifest:
        print("  [manifest] FAIL: 'strand_specificities' missing from manifest")
        ok = False

    if ok:
        print("  [manifest] PASS: gDNA config correctly stored in manifest")
    return ok


def test_distinct_configs_produce_different_distributions(
    fasta_path: Path, transcripts: list[Transcript],
) -> bool:
    """Two different gDNA configs should produce clearly different FL distributions."""
    sim_params = SimulationParams(
        n_rna_fragments=100, sim_seed=42,
        frag_mean=250.0, frag_std=50.0, frag_min=50, frag_max=1000,
        read_length=101, error_rate=0.0,
    )

    configs = [
        ("short", GDNASimConfig(frag_mean=150, frag_std=20, frag_min=100, frag_max=300)),
        ("long", GDNASimConfig(frag_mean=500, frag_std=50, frag_min=300, frag_max=800)),
    ]

    means = {}
    for label, gdna_cfg in configs:
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "out"
            sim = WholeGenomeSimulator(
                fasta_path, transcripts, sim_params, gdna_cfg,
                strand_specificity=0.9, seed=42,
            )
            _, _, bam_path = sim.simulate_and_write(outdir, n_rna=100, n_gdna=2000)
            sim.close()

            records = _parse_gdna_bam(str(bam_path))
            frag_lens = np.array([r["frag_len_from_name"] for r in records])
            means[label] = frag_lens.mean()

    diff = abs(means["long"] - means["short"])
    if diff < 100:
        print(
            f"  [distinct] FAIL: short mean={means['short']:.1f}, "
            f"long mean={means['long']:.1f}, diff={diff:.1f} (expected >100)"
        )
        return False

    print(
        f"  [distinct] PASS: short mean={means['short']:.1f}, "
        f"long mean={means['long']:.1f}, diff={diff:.1f}"
    )
    return True


# ── Main ─────────────────────────────────────────────────────────────


def main():
    print("=" * 60)
    print("gDNA Fragment Simulation Verification")
    print("=" * 60)

    # Setup: create temp genome + transcripts
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = Path(tmpdir) / "genome.fa"
        chrom, _ = _make_test_fasta(fasta_path)
        transcripts = _make_test_transcripts(chrom)

        all_ok = True

        # Test 1: Default gDNA config (mean=350, std=100, min=100, max=1000)
        print("\n1. Default GDNASimConfig (mean=350, std=100, min=100)")
        cfg1 = GDNASimConfig()  # defaults
        ok = test_gdna_fragment_distribution(
            fasta_path, transcripts, cfg1, n_gdna=5000, label="default",
        )
        all_ok &= ok

        # Test 2: Custom config (mean=200, std=30, min=100, max=400)
        print("\n2. Custom GDNASimConfig (mean=200, std=30, min=100, max=400)")
        cfg2 = GDNASimConfig(frag_mean=200, frag_std=30, frag_min=100, frag_max=400)
        ok = test_gdna_fragment_distribution(
            fasta_path, transcripts, cfg2, n_gdna=5000, label="custom",
        )
        all_ok &= ok

        # Test 3: Tight config — should produce narrow distribution
        print("\n3. Tight GDNASimConfig (mean=300, std=5, min=280, max=320)")
        cfg3 = GDNASimConfig(frag_mean=300, frag_std=5, frag_min=280, frag_max=320)
        ok = test_gdna_fragment_distribution(
            fasta_path, transcripts, cfg3, n_gdna=2000, label="tight",
        )
        all_ok &= ok

        # Test 4: Manifest stores gDNA config
        print("\n4. Manifest stores gDNA config")
        ok = test_manifest_stores_gdna_config()
        all_ok &= ok

        # Test 5: Different configs → different distributions
        print("\n5. Distinct configs produce different distributions")
        ok = test_distinct_configs_produce_different_distributions(
            fasta_path, transcripts,
        )
        all_ok &= ok

    print("\n" + "=" * 60)
    if all_ok:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print("=" * 60)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
