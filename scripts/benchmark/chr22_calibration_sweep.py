"""Calibration stress sweep on human chr22 (multimapper-aware).

Counterpart to ``calibration_sweep.py`` (the mini-genome sweep), but
uses the *real* human chr22 sequence + GENCODE annotation so that the
pipeline is exposed to pseudogenes, segmental duplications, and other
repetitive regions that produce multimapping reads.  The pre-built
production rigel index (``rigel_index``) is reused directly — no
re-indexing step.

Sweep axes::

    strand specificity (SS)  ∈ {0.5, 0.9, 0.99, 1.0}
    gDNA fraction-of-RNA     ∈ {0.0, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0}
    RNG seed                 ∈ {1, 2, 3}

For each (SS, gdna_fraction, seed) cell:

  1. Build an oracle BAM of chr22 reads with controlled gDNA contamination.
  2. Run ``run_pipeline`` against the full production index.
  3. Extract calibration estimates + truth λ_G.

Results TSV is written to ``--out`` and is consumed by the existing
``analyze_calibration_sweep.py`` analyser.

Run::

    python scripts/benchmark/chr22_calibration_sweep.py \\
        --out scripts/benchmark/results/calibration_v4_chr22 \\
        [--quick]
"""
from __future__ import annotations

import argparse
import json
import logging
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from rigel.calibration import calibrate_gdna
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, OracleBamSimulator, SimConfig
from rigel.transcript import Transcript

logger = logging.getLogger("chr22_calibration_sweep")


REF_NAME = "chr22"
DEFAULT_FASTA = Path(
    "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/chr22_only/chr22.fasta.bgz"
)
DEFAULT_INDEX = Path(
    "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/chr22_only/rigel_index"
)


# ------------------------------------------------------------------ genome

@dataclass
class _StringGenome:
    """Minimal genome shim for OracleBamSimulator/ReadSimulator.

    Only ``__len__``, ``__getitem__`` (int/slice) and the ``name``
    attribute are consumed by the simulators.
    """

    seq: str
    name: str

    def __len__(self) -> int:
        return len(self.seq)

    def __getitem__(self, key):
        if isinstance(key, slice):
            return self.seq[key]
        return self.seq[key]


def load_chr22_genome(fasta_path: Path) -> _StringGenome:
    logger.info("Loading %s sequence from %s ...", REF_NAME, fasta_path)
    with pysam.FastaFile(str(fasta_path)) as fa:
        if REF_NAME not in fa.references:
            raise ValueError(
                f"Expected {REF_NAME!r} in FASTA references; "
                f"found {list(fa.references)[:5]} ..."
            )
        seq = fa.fetch(REF_NAME).upper()
    logger.info("  %s length = %d bp", REF_NAME, len(seq))
    return _StringGenome(seq=seq, name=REF_NAME)


# ------------------------------------------------------------------ transcripts

def select_transcripts(
    index: TranscriptIndex,
    *,
    seed: int,
    max_transcripts: int,
    min_exons: int = 2,
) -> list[Transcript]:
    """Pick a subset of annotated multi-exon chr22 transcripts with abundances.

    Only annotated (non-synthetic, non-nRNA) transcripts are selected so
    we have an exon-based simulation substrate.  Abundances are drawn
    log-uniform from [5, 500] so the locus mix spans ~2 orders of
    magnitude.  ``max_transcripts`` caps the number sampled.
    """
    rng = np.random.default_rng(seed)
    t_df = index.t_df
    mask = (
        (t_df["ref"] == REF_NAME)
        & (~t_df["is_synthetic"])
        & (~t_df["is_nrna"])
        & (t_df["n_exons"] >= min_exons)
    )
    df = t_df.loc[mask].reset_index(drop=True)
    if len(df) == 0:
        raise RuntimeError(f"No annotated chr22 transcripts in index.")
    if len(df) > max_transcripts:
        keep = rng.choice(len(df), size=max_transcripts, replace=False)
        df = df.iloc[keep].reset_index(drop=True)

    # Pull exon structure from the index.  ``get_exon_intervals`` returns
    # an (n_exons, 2) int array per t_index.
    abundances = np.exp(rng.uniform(np.log(5.0), np.log(500.0), size=len(df)))
    out: list[Transcript] = []
    skipped = 0
    from rigel.types import Interval, Strand

    for row, abund in zip(df.itertuples(index=False), abundances):
        ex = index.get_exon_intervals(int(row.t_index))
        if ex is None or len(ex) < min_exons:
            skipped += 1
            continue
        strand = row.strand
        if isinstance(strand, (int, np.integer)):
            strand = Strand(int(strand))
        elif isinstance(strand, str):
            strand = Strand.from_str(strand)
        exons = [Interval(int(s), int(e)) for s, e in ex]
        out.append(Transcript(
            ref=REF_NAME,
            strand=strand,
            exons=exons,
            length=int(row.length) if not np.isnan(row.length) else None,
            t_id=str(row.t_id),
            g_id=str(row.g_id),
            g_name=str(row.g_name) if row.g_name else None,
            g_type=str(row.g_type) if row.g_type else None,
            t_index=int(row.t_index),
            g_index=int(row.g_index),
            is_basic=bool(row.is_basic),
            is_mane=bool(row.is_mane),
            is_ccds=bool(row.is_ccds),
            is_nrna=False,
            is_synthetic=False,
            abundance=float(abund),
            nrna_abundance=0.0,
        ))
    if skipped:
        logger.info("Skipped %d transcripts without cached exon intervals", skipped)
    logger.info("Selected %d chr22 transcripts (requested ≤%d)",
                len(out), max_transcripts)
    return out


# ------------------------------------------------------------------ truth

def truth_lambda_gdna(index: TranscriptIndex, n_gdna_truth: int,
                      ref: str = REF_NAME,
                      genome_len: int | None = None) -> float:
    """Truth λ_G on the pipeline's calibration domain.

    The pipeline's λ̂ is a density with units of fragments per
    *mappable* bp (``λ̂ ≈ a_i / E_i`` with ``E_i`` =
    ``mappable_effective_length``).

    ``OracleBamSimulator`` scatters gDNA fragments uniformly across
    the *whole genome* (here, all of ``ref``).  The expected number
    of fragment starts landing at a mappable position in region ``i``
    is ``n_gdna · E_i / L_genome``.  Summing over ``i`` and dividing
    by ``ΣE_i``:

        λ_true = n_gdna / L_genome

    NOT ``n_gdna / ΣE_i``.  The distinction matters whenever a
    non-trivial fraction of the simulated genome is unmappable:
    using ``ΣE_i`` in the denominator inflates truth by
    ``L_genome / ΣE_i`` (≈ 1.41× for chr22 at 100 bp reads).
    """
    if genome_len is None:
        r = index.region_df
        # Calibration regions tile the reference fully, so Σ region length
        # equals the genome length for that ref.
        genome_len = int(r.loc[r["ref"] == ref, "length"].sum())
    if genome_len <= 0:
        return float("nan")
    return float(n_gdna_truth) / float(genome_len)


# ------------------------------------------------------------------ runner

def run_one(
    out_dir: Path,
    *,
    genome: _StringGenome,
    transcripts: list[Transcript],
    index: TranscriptIndex,
    ss: float,
    gdna_fraction: float,
    seed: int,
    n_rna_fragments: int,
    rec_idx: int,
    keep_bam: bool = False,
) -> dict:
    label = f"ss{ss:.2f}_gdna{gdna_fraction:.2f}_seed{seed}"
    scen_dir = out_dir / "scenarios" / label
    scen_dir.mkdir(parents=True, exist_ok=True)
    bam_path = scen_dir / "oracle.bam"

    sim_cfg = SimConfig(
        frag_mean=200, frag_std=40, frag_min=80, frag_max=500,
        read_length=100, strand_specificity=ss, seed=seed,
    )
    gdna_cfg = (
        GDNAConfig(abundance=10.0, frag_mean=200, frag_std=40,
                   frag_min=80, frag_max=500)
        if gdna_fraction > 0 else None
    )

    oracle = OracleBamSimulator(
        genome, transcripts,
        config=sim_cfg, gdna_config=gdna_cfg, ref_name=REF_NAME,
    )

    n_mrna, n_nrna = oracle._sim.compute_rna_split(n_rna_fragments)
    n_gdna = int(round(n_rna_fragments * gdna_fraction)) if gdna_cfg else 0
    pool_split = (n_mrna, n_nrna, n_gdna)

    t0 = time.perf_counter()
    oracle.write_bam(bam_path, n_fragments=sum(pool_split),
                     name_sorted=True, pool_split=pool_split)
    t_sim = time.perf_counter() - t0

    n_rna_truth = n_mrna + n_nrna
    n_gdna_truth = n_gdna
    lam_truth = truth_lambda_gdna(index, n_gdna_truth)

    # Oracle BAM writes STAR-convention XS; "auto" also works.
    cfg = PipelineConfig(
        em=EMConfig(seed=seed),
        scan=BamScanConfig(sj_strand_tag="XS"),
    )
    t1 = time.perf_counter()
    try:
        pr = run_pipeline(str(bam_path), index, config=cfg)
        ok = True
        err = ""
    except Exception as exc:  # noqa: BLE001
        ok = False
        err = repr(exc)
        pr = None
    t_pipe = time.perf_counter() - t1

    row: dict = {
        "label": label,
        "ss": ss,
        "gdna_fraction": gdna_fraction,
        "seed": seed,
        "n_rna_truth": n_rna_truth,
        "n_gdna_truth": n_gdna_truth,
        "lambda_truth": lam_truth,
        "t_sim_s": t_sim,
        "t_pipe_s": t_pipe,
        "ok": ok,
        "err": err,
        "rec_idx": rec_idx,
    }

    if pr is not None:
        cal = pr.calibration
        est = pr.estimator
        row.update({
            "ss_estimated": (None if cal is None
                             else float(cal.strand_specificity)),
            "lambda_pool": None if cal is None else cal.lambda_gdna,
            "mu_R": None if cal is None else cal.mu_R,
            "sigma_R": None if cal is None else cal.sigma_R,
            "mixing_pi": None if cal is None else cal.mixing_pi,
            "mixing_pi_soft": (None if cal is None
                               else cal.mixing_pi_soft),
            "strand_used": None if cal is None else cal.strand_used,
            "strand_z": None if cal is None else cal.strand_z,
            "em_n_iter": None if cal is None else cal.em_n_iter,
            "em_converged": None if cal is None else cal.em_converged,
            "n_eligible": None if cal is None else cal.n_eligible,
            "n_soft": None if cal is None else cal.n_soft,
            "n_spliced_hard": (None if cal is None
                               else cal.n_spliced_hard),
            "gdna_em_count": float(est.gdna_em_count),
            "gdna_intergenic": int(pr.stats.n_intergenic),
            "gdna_total": int(pr.stats.n_intergenic) + float(est.gdna_em_count),
            "nrna_em_count": float(est.nrna_em_count),
            "gdna_contamination_rate": float(est.gdna_contamination_rate),
            "n_total_fragments": (None if cal is None
                                  else float(cal.region_n_total.sum())),
        })
        import math
        lp = row["lambda_pool"]
        lt = row["lambda_truth"]
        if lp is not None and lt is not None and lt > 0:
            row["lambda_relative_error"] = (lp - lt) / lt
            row["lambda_log2_ratio"] = math.log2(
                max(lp, 1e-30) / max(lt, 1e-30)
            )
        else:
            row["lambda_relative_error"] = None
            row["lambda_log2_ratio"] = None
        if n_gdna_truth > 0:
            row["gdna_em_relative_error"] = (
                (row["gdna_em_count"] - n_gdna_truth) / n_gdna_truth
            )
            row["gdna_total_relative_error"] = (
                (row["gdna_total"] - n_gdna_truth) / n_gdna_truth
            )
        else:
            row["gdna_em_relative_error"] = None
            row["gdna_total_relative_error"] = None

    if not keep_bam:
        try:
            bam_path.unlink(missing_ok=True)
        except OSError:
            pass

    return row


# ------------------------------------------------------------------ driver

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--n-rna-fragments", type=int, default=200_000,
                    help="RNA fragment count per cell (default 200k).")
    ap.add_argument("--max-transcripts", type=int, default=2000,
                    help="Cap on expressed chr22 transcripts per seed.")
    ap.add_argument("--quick", action="store_true",
                    help="Run a fast subset for smoke testing.")
    ap.add_argument("--match-star-grid", action="store_true",
                    help="Use the same 3×5×2 grid as chr22_star_calibration_sweep.")
    ap.add_argument("--keep-bam", action="store_true",
                    help="Retain per-cell oracle BAMs (default: delete).")
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args()

    logging.basicConfig(
        level=args.log_level,
        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    )

    args.out.mkdir(parents=True, exist_ok=True)

    if args.quick:
        ss_grid = [0.5, 0.9, 1.0]
        gdna_grid = [0.0, 0.5, 1.0, 2.0]
        seeds = [1]
        n_rna = max(50_000, args.n_rna_fragments // 4)
    elif args.match_star_grid:
        # Grid matched to chr22_star_calibration_sweep.py (3×5×2 = 30 cells)
        ss_grid = [0.5, 0.9, 1.0]
        gdna_grid = [0.0, 0.05, 0.25, 1.0, 2.0]
        seeds = [1, 2]
        n_rna = args.n_rna_fragments
    else:
        ss_grid = [0.5, 0.9, 0.99, 1.0]
        gdna_grid = [0.0, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0]
        seeds = [1, 2, 3]
        n_rna = args.n_rna_fragments

    logger.info("chr22 calibration sweep — %d SS × %d gDNA × %d seeds = %d cells",
                len(ss_grid), len(gdna_grid), len(seeds),
                len(ss_grid) * len(gdna_grid) * len(seeds))
    logger.info("  n_rna_fragments/cell = %d", n_rna)

    logger.info("Loading rigel index from %s ...", args.index)
    index = TranscriptIndex.load(args.index)

    genome = load_chr22_genome(args.fasta)

    # Pick transcripts once per seed — keep the expressed panel stable
    # across the SS×gDNA grid for that seed.
    panels: dict[int, list[Transcript]] = {}
    for s in seeds:
        panels[s] = select_transcripts(
            index, seed=s, max_transcripts=args.max_transcripts,
        )

    rows: list[dict] = []
    rec_idx = 0
    for ss in ss_grid:
        for gd in gdna_grid:
            for s in seeds:
                logger.info(
                    "[%d] ss=%.2f gdna=%.2f seed=%d", rec_idx, ss, gd, s,
                )
                row = run_one(
                    args.out,
                    genome=genome,
                    transcripts=panels[s],
                    index=index,
                    ss=ss, gdna_fraction=gd, seed=s,
                    n_rna_fragments=n_rna,
                    rec_idx=rec_idx,
                    keep_bam=args.keep_bam,
                )
                rows.append(row)
                rec_idx += 1
                # Flush progressively so partial runs are analysable.
                pd.DataFrame(rows).to_csv(
                    args.out / "results.tsv", sep="\t", index=False,
                )

    df = pd.DataFrame(rows)
    df.to_csv(args.out / "results.tsv", sep="\t", index=False)
    (args.out / "results.json").write_text(
        json.dumps(rows, indent=2, default=str)
    )
    logger.info("Wrote %d rows → %s", len(df), args.out / "results.tsv")


if __name__ == "__main__":
    main()
