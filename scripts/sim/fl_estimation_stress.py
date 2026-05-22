#!/usr/bin/env python3
"""Stress-test Rigel RNA/gDNA fragment-length estimation on oracle mini-genomes.

This script intentionally stops after scan + calibration because the Phase 4
locus-prior/EM path is still under the fractional cutover gate. It still
exercises the production FL path:

1. `rigel.sim.Scenario.build_oracle` creates a mini-genome, arbitrary spliced
   transcripts, and a name-sorted perfect BAM.
2. `pipeline.scan_and_buffer` runs the native BAM scanner and fractional region
   accumulator.
3. `calibration.calibrate` builds finalized RNA/gDNA/global FL models.
4. The fitted FL means are compared with true sampled fragment lengths encoded
   in oracle BAM read names.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pysam

from rigel.calibration import calibrate
from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.sim import GDNAConfig, Scenario, SimConfig
from rigel.sim.truth import parse_origin

REPO_ROOT = Path(__file__).resolve().parents[2]
LOGGER = logging.getLogger("fl_estimation_stress")


@dataclass(frozen=True, slots=True)
class FLCase:
    name: str
    rna_mean: float
    rna_std: float
    gdna_mean: float
    gdna_std: float


DEFAULT_CASES: tuple[FLCase, ...] = (
    FLCase("rna150_gdna70", 150.0, 20.0, 70.0, 20.0),
    FLCase("rna70_gdna150", 70.0, 20.0, 150.0, 20.0),
    FLCase("rna150_gdna150", 150.0, 20.0, 150.0, 20.0),
    FLCase("rna250_gdna80", 250.0, 35.0, 80.0, 20.0),
    FLCase("rna80_gdna250", 80.0, 20.0, 250.0, 35.0),
)


@dataclass(slots=True)
class FLMetrics:
    case: str
    seed: int
    n_rna_requested: int
    gdna_fraction: float
    n_truth_mrna: int
    n_truth_gdna: int
    n_scan_fragments: int
    n_observed: int
    n_unannotated_ref: int
    n_fl_unavailable: int
    n_fl_rna: float
    n_fl_gdna: float
    rna_quality: str
    gdna_quality: str
    requested_rna_mean: float
    requested_gdna_mean: float
    truth_rna_mean: float
    truth_gdna_mean: float
    fitted_rna_mean: float
    fitted_gdna_mean: float
    scoring_rna_mean: float
    scoring_gdna_mean: float
    truth_rna_std: float
    truth_gdna_std: float
    fitted_rna_std: float
    fitted_gdna_std: float
    rna_mean_error: float
    gdna_mean_error: float
    rna_mean_abs_error: float
    gdna_mean_abs_error: float
    rna_mean_rel_error: float
    gdna_mean_rel_error: float
    rna_w1_bp: float
    gdna_w1_bp: float
    rna_l1: float
    gdna_l1: float
    rna_training_recovery: float
    gdna_pool_recovery: float


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def _nonnegative_float(value: str) -> float:
    parsed = float(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be non-negative")
    return parsed


def _frag_min(mean: float, std: float, lower_floor: int) -> int:
    return max(lower_floor, int(math.floor(mean - 4.0 * std)))


def _add_arbitrary_transcripts(
    scenario: Scenario,
    *,
    n_genes: int,
    genome_length: int,
) -> None:
    """Populate a mini-genome with non-overlapping multi-exon transcripts.

    Exons are intentionally short so even short RNA fragments frequently cross
    annotated splice junctions and contribute to the RNA FL training channel.
    """
    gene_spacing = max(1_500, genome_length // max(n_genes + 1, 1))
    for gene_idx in range(n_genes):
        base = 500 + gene_idx * gene_spacing
        exon_len = 45 + 5 * (gene_idx % 3)
        intron_len = 70 + 15 * (gene_idx % 4)
        n_exons = 7 + (gene_idx % 3)
        exons: list[tuple[int, int]] = []
        cursor = base
        for _ in range(n_exons):
            exons.append((cursor, cursor + exon_len))
            cursor += exon_len + intron_len
        if cursor + 500 >= genome_length:
            raise ValueError(
                "Mini-genome too short for requested transcript layout; "
                f"increase --genome-length above {genome_length}."
            )
        strand = "+" if gene_idx % 2 == 0 else "-"
        scenario.add_gene(
            f"gene{gene_idx}",
            strand,
            [
                {
                    "t_id": f"tx{gene_idx}",
                    "exons": exons,
                    "abundance": 100.0 + 10.0 * (gene_idx % 5),
                    "nrna_abundance": 0.0,
                }
            ],
            gene_name=f"Gene{gene_idx}",
        )


def _truth_lengths_from_bam(bam_path: Path) -> dict[str, np.ndarray]:
    seen: set[str] = set()
    lengths: dict[str, list[int]] = {"mrna": [], "gdna": [], "nrna": []}
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            qname = read.query_name
            if qname in seen:
                continue
            seen.add(qname)
            origin = parse_origin(qname)
            if origin.start is None or origin.end is None:
                continue
            lengths[origin.kind].append(origin.end - origin.start)
    return {kind: np.asarray(vals, dtype=np.int64) for kind, vals in lengths.items()}


def _truth_pmf(lengths: np.ndarray, max_size: int) -> np.ndarray:
    counts = np.zeros(max_size + 1, dtype=np.float64)
    if lengths.size == 0:
        return counts
    valid = lengths[(lengths >= 0) & (lengths <= max_size)]
    if valid.size == 0:
        return counts
    bincount = np.bincount(valid, minlength=max_size + 1).astype(np.float64)
    counts[: max_size + 1] = bincount[: max_size + 1]
    return counts / counts.sum()


def _model_pmf(model, max_size: int) -> np.ndarray:
    pmf = np.asarray(model.pmf, dtype=np.float64)
    out = np.zeros(max_size + 1, dtype=np.float64)
    n = min(out.size, pmf.size)
    out[:n] = pmf[:n]
    total = out.sum()
    return out / total if total > 0 else out


def _pmf_mean(pmf: np.ndarray) -> float:
    if pmf.sum() <= 0:
        return 0.0
    return float(np.dot(np.arange(pmf.size, dtype=np.float64), pmf))


def _l1(p: np.ndarray, q: np.ndarray) -> float:
    return float(np.abs(p - q).sum())


def _w1_bp(p: np.ndarray, q: np.ndarray) -> float:
    """1D Wasserstein distance for unit-spaced FL bins, in bp."""
    return float(np.abs(np.cumsum(p) - np.cumsum(q)).sum())


def _mean(values: np.ndarray) -> float:
    return float(values.mean()) if values.size else 0.0


def _std(values: np.ndarray) -> float:
    return float(values.std(ddof=0)) if values.size else 0.0


def _rel_error(error: float, truth: float) -> float:
    return float(error / truth) if truth else 0.0


def _run_case(
    case: FLCase,
    *,
    case_dir: Path,
    seed: int,
    n_rna: int,
    gdna_fraction: float,
    genome_length: int,
    n_genes: int,
    read_length: int,
    max_frag_length: int,
    prior_ess: float,
    pool_quality_good: int,
    pool_quality_weak: int,
    threads: int,
) -> FLMetrics:
    LOGGER.info("[%s] building oracle mini-genome", case.name)
    scenario = Scenario(
        case.name,
        genome_length=genome_length,
        seed=seed,
        work_dir=case_dir,
        ref_name="chr1",
        gdna_config=GDNAConfig(
            abundance=1.0,
            frag_mean=case.gdna_mean,
            frag_std=case.gdna_std,
            frag_min=_frag_min(case.gdna_mean, case.gdna_std, lower_floor=30),
            frag_max=max_frag_length,
        ),
    )
    _add_arbitrary_transcripts(
        scenario,
        n_genes=n_genes,
        genome_length=genome_length,
    )

    sim_config = SimConfig(
        frag_mean=case.rna_mean,
        frag_std=case.rna_std,
        frag_min=_frag_min(case.rna_mean, case.rna_std, lower_floor=30),
        frag_max=max_frag_length,
        read_length=read_length,
        error_rate=0.0,
        strand_specificity=1.0,
        seed=seed,
    )
    result = scenario.build_oracle(
        n_rna_fragments=n_rna,
        gdna_fraction=gdna_fraction,
        sim_config=sim_config,
    )
    truth_lengths = _truth_lengths_from_bam(result.bam_path)

    scan_config = BamScanConfig(
        sj_strand_tag="XS",
        include_multimap=True,
        max_frag_length=max_frag_length,
        total_threads=threads,
        bgzf_threads=min(2, max(1, threads)),
        fragments_per_chunk=250_000,
        read_name_batch_size=512,
        buffer_size_bytes=0,
    )
    stats, strand_models, scan_trained, buffer, payload = scan_and_buffer(
        str(result.bam_path),
        result.index,
        scan_config,
    )
    try:
        strand_models.finalize()
        if payload is None:
            raise RuntimeError("scan did not produce a calibration payload")
        calibration = calibrate(
            index=result.index,
            payload=payload,
            scan_trained=scan_trained,
            fl_prior_ess=prior_ess,
            pool_quality_good=pool_quality_good,
            pool_quality_weak=pool_quality_weak,
            resolver_splicing_anchor_tolerance=int(payload.resolver_splicing_anchor_tolerance),
        )
    finally:
        buffer.cleanup()

    fl = calibration.fl_models
    rna_truth = truth_lengths["mrna"]
    gdna_truth = truth_lengths["gdna"]
    rna_truth_pmf = _truth_pmf(rna_truth, max_frag_length)
    gdna_truth_pmf = _truth_pmf(gdna_truth, max_frag_length)
    rna_model_pmf = _model_pmf(fl.rna, max_frag_length)
    gdna_model_pmf = _model_pmf(fl.gdna, max_frag_length)

    truth_rna_mean = _mean(rna_truth)
    truth_gdna_mean = _mean(gdna_truth)
    rna_error = float(fl.rna.mean - truth_rna_mean)
    gdna_error = float(fl.gdna.mean - truth_gdna_mean)

    metrics = FLMetrics(
        case=case.name,
        seed=seed,
        n_rna_requested=n_rna,
        gdna_fraction=gdna_fraction,
        n_truth_mrna=int(rna_truth.size),
        n_truth_gdna=int(gdna_truth.size),
        n_scan_fragments=int(stats.n_fragments),
        n_observed=int(payload.n_observed),
        n_unannotated_ref=int(payload.n_unannotated_ref),
        n_fl_unavailable=int(payload.n_fl_unavailable),
        n_fl_rna=float(fl.n_rna),
        n_fl_gdna=float(fl.n_gdna),
        rna_quality=fl.rna_quality,
        gdna_quality=fl.gdna_quality,
        requested_rna_mean=case.rna_mean,
        requested_gdna_mean=case.gdna_mean,
        truth_rna_mean=truth_rna_mean,
        truth_gdna_mean=truth_gdna_mean,
        fitted_rna_mean=float(fl.rna.mean),
        fitted_gdna_mean=float(fl.gdna.mean),
        scoring_rna_mean=_pmf_mean(rna_model_pmf),
        scoring_gdna_mean=_pmf_mean(gdna_model_pmf),
        truth_rna_std=_std(rna_truth),
        truth_gdna_std=_std(gdna_truth),
        fitted_rna_std=float(fl.rna.std),
        fitted_gdna_std=float(fl.gdna.std),
        rna_mean_error=rna_error,
        gdna_mean_error=gdna_error,
        rna_mean_abs_error=abs(rna_error),
        gdna_mean_abs_error=abs(gdna_error),
        rna_mean_rel_error=_rel_error(rna_error, truth_rna_mean),
        gdna_mean_rel_error=_rel_error(gdna_error, truth_gdna_mean),
        rna_w1_bp=_w1_bp(rna_model_pmf, rna_truth_pmf),
        gdna_w1_bp=_w1_bp(gdna_model_pmf, gdna_truth_pmf),
        rna_l1=_l1(rna_model_pmf, rna_truth_pmf),
        gdna_l1=_l1(gdna_model_pmf, gdna_truth_pmf),
        rna_training_recovery=(float(fl.n_rna) / float(rna_truth.size) if rna_truth.size else 0.0),
        gdna_pool_recovery=(float(fl.n_gdna) / float(gdna_truth.size) if gdna_truth.size else 0.0),
    )
    LOGGER.info(
        "[%s] RNA mean truth=%.2f fitted=%.2f err=%+.2f; gDNA truth=%.2f fitted=%.2f err=%+.2f",
        case.name,
        metrics.truth_rna_mean,
        metrics.fitted_rna_mean,
        metrics.rna_mean_error,
        metrics.truth_gdna_mean,
        metrics.fitted_gdna_mean,
        metrics.gdna_mean_error,
    )
    return metrics


def _write_tsv(path: Path, rows: Iterable[FLMetrics]) -> None:
    rows = list(rows)
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def _write_json(path: Path, rows: list[FLMetrics], args: argparse.Namespace) -> None:
    args_dict = {
        key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()
    }
    payload = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "args": args_dict,
        "cases": [asdict(row) for row in rows],
        "summary": {
            "max_rna_abs_mean_error": max((r.rna_mean_abs_error for r in rows), default=0.0),
            "max_gdna_abs_mean_error": max((r.gdna_mean_abs_error for r in rows), default=0.0),
            "mean_rna_abs_mean_error": float(np.mean([r.rna_mean_abs_error for r in rows]))
            if rows
            else 0.0,
            "mean_gdna_abs_mean_error": float(np.mean([r.gdna_mean_abs_error for r in rows]))
            if rows
            else 0.0,
        },
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _parse_cases(selected: list[str] | None) -> tuple[FLCase, ...]:
    if not selected:
        return DEFAULT_CASES
    by_name = {case.name: case for case in DEFAULT_CASES}
    missing = [name for name in selected if name not in by_name]
    if missing:
        raise SystemExit(
            "Unknown case(s): " + ", ".join(missing) + "; valid cases: " + ", ".join(by_name)
        )
    return tuple(by_name[name] for name in selected)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Stress-test RNA and gDNA FL estimation on oracle mini-genomes."
    )
    default_out = (
        REPO_ROOT / "results" / ("fl_estimation_stress_" + datetime.now().strftime("%Y%m%d_%H%M%S"))
    )
    parser.add_argument("--out-dir", type=Path, default=default_out)
    parser.add_argument("--case", action="append", choices=[c.name for c in DEFAULT_CASES])
    parser.add_argument("--n-rna", type=_positive_int, default=6_000)
    parser.add_argument("--gdna-fraction", type=_nonnegative_float, default=1.0)
    parser.add_argument("--genome-length", type=_positive_int, default=80_000)
    parser.add_argument("--n-genes", type=_positive_int, default=12)
    parser.add_argument("--read-length", type=_positive_int, default=75)
    parser.add_argument("--max-frag-length", type=_positive_int, default=600)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--threads", type=_positive_int, default=2)
    parser.add_argument("--prior-ess", type=float, default=1_000.0)
    parser.add_argument("--pool-quality-good", type=_positive_int, default=5_000)
    parser.add_argument("--pool-quality-weak", type=_positive_int, default=1)
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING"])
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    cases = _parse_cases(args.case)
    LOGGER.info("Writing outputs to %s", args.out_dir)

    rows: list[FLMetrics] = []
    for offset, case in enumerate(cases):
        case_dir = args.out_dir / case.name
        case_seed = int(args.seed + offset * 1_003)
        rows.append(
            _run_case(
                case,
                case_dir=case_dir,
                seed=case_seed,
                n_rna=args.n_rna,
                gdna_fraction=args.gdna_fraction,
                genome_length=args.genome_length,
                n_genes=args.n_genes,
                read_length=args.read_length,
                max_frag_length=args.max_frag_length,
                prior_ess=args.prior_ess,
                pool_quality_good=args.pool_quality_good,
                pool_quality_weak=args.pool_quality_weak,
                threads=args.threads,
            )
        )

    _write_tsv(args.out_dir / "fl_estimation_metrics.tsv", rows)
    _write_json(args.out_dir / "fl_estimation_metrics.json", rows, args)

    print("\nFL estimation stress results")
    print(f"Output: {args.out_dir}")
    print(
        "case\trna_truth\trna_fit\trna_err\tgdna_truth\tgdna_fit\tgdna_err\t"
        "rna_q\tgdna_q\trna_recovery\tgdna_recovery"
    )
    for row in rows:
        print(
            f"{row.case}\t{row.truth_rna_mean:.2f}\t{row.fitted_rna_mean:.2f}\t"
            f"{row.rna_mean_error:+.2f}\t{row.truth_gdna_mean:.2f}\t"
            f"{row.fitted_gdna_mean:.2f}\t{row.gdna_mean_error:+.2f}\t"
            f"{row.rna_quality}\t{row.gdna_quality}\t"
            f"{row.rna_training_recovery:.3f}\t{row.gdna_pool_recovery:.3f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
