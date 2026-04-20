"""Run the calibration stress sweep on the mini genome.

For each (strand_specificity, gdna_fraction) cell:

1. Generate an oracle BAM using :class:`OracleBamSimulator` with the
   per-transcript abundances baked into the mini GTF.
2. Run :func:`rigel.pipeline.scan_and_buffer` on that BAM to obtain
   region counts, fragment-length table, and trained strand/FL models.
3. Call :func:`rigel.calibration.calibrate_gdna` to produce the v4
   EM mixture fit (lam_G, pi, gamma, mu_R, sigma_R, …).
4. Compute *truth* λ_G = n_gdna / L_genome (uniform scatter) and
   per-region truth counts (mrna / nrna / gdna) from the oracle
   read names.
5. Emit:
     - ``sweep_results.tsv`` — one row per cell with summary metrics.
     - ``per_region/<label>.feather`` — per-region table: region_id,
       kind, rep_id, mappable_bp, n_u, n_s, n_gdna_true, gamma, E_gdna,
       posterior_class, …
     - ``per_read/<label>.feather`` — optional per-read truth vs.
       assignment trace (enabled with ``--per-read``).

Usage
-----
    python -m scripts.calibration.stress_mini.run_sweep \\
        --reference /tmp/mini_stress --out /tmp/mini_stress/sweep
"""
from __future__ import annotations

import argparse
import logging
import re
import sys
import time
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from rigel.calibration import calibrate_gdna
from rigel.config import BamScanConfig, PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer
from rigel.sim import GDNAConfig, SimConfig
from rigel.sim.oracle_bam import OracleBamSimulator
from rigel.sim.genome import MutableGenome
from rigel.types import Interval, Strand
from rigel.transcript import Transcript

logger = logging.getLogger("run_sweep")


# ------------------------------------------------------------------ config

DEFAULT_SS_LEVELS = (0.5, 0.75, 0.9, 1.0)
DEFAULT_GDNA_LEVELS = (0.0, 0.5, 1.0, 1.5, 2.0)
DEFAULT_N_RNA = 200_000          # baseline mature RNA fragments per cell
READ_LENGTH = 100
FRAG_MEAN = 200
FRAG_SD = 40
FRAG_MIN = 80
FRAG_MAX = 500
GDNA_FRAG_MEAN = 350
GDNA_FRAG_SD = 120
GDNA_FRAG_MIN = 100
GDNA_FRAG_MAX = 1000


# ------------------------------------------------------------------ helpers


class _StringGenome:
    """Minimal shim that satisfies :class:`OracleBamSimulator`'s genome API."""

    def __init__(self, seq: str, name: str) -> None:
        self.seq = seq
        self.name = name

    def __len__(self) -> int:
        return len(self.seq)

    def __getitem__(self, key):
        return self.seq[key]


def _load_genome(fasta: Path) -> _StringGenome:
    """Load the mini genome FASTA back into memory as a shim."""
    with pysam.FastaFile(str(fasta)) as fa:
        refs = fa.references
        if len(refs) != 1:
            raise RuntimeError(f"Expected single-ref genome, got {refs}")
        seq = fa.fetch(refs[0])
    return _StringGenome(seq=seq, name=refs[0])


def _transcripts_from_index(index: TranscriptIndex) -> list[Transcript]:
    """Reconstruct :class:`Transcript` objects for the simulator."""
    tdf = index.t_df
    out: list[Transcript] = []
    for row in tdf.itertuples(index=False):
        strand = row.strand
        if isinstance(strand, (int, np.integer)):
            strand = Strand(int(strand))
        elif isinstance(strand, str):
            strand = Strand.from_str(strand)
        ex = index.get_exon_intervals(int(row.t_index))
        if ex is None or len(ex) == 0:
            continue
        exons = [Interval(int(s), int(e)) for s, e in ex]
        out.append(Transcript(
            ref=str(row.ref),
            strand=strand,
            exons=exons,
            length=int(row.length),
            t_id=str(row.t_id),
            g_id=str(row.g_id),
            g_name=str(row.g_name) if row.g_name else None,
            g_type=str(row.g_type) if row.g_type else None,
            t_index=int(row.t_index),
            g_index=int(row.g_index),
            is_basic=bool(row.is_basic),
            is_mane=bool(row.is_mane),
            is_ccds=bool(row.is_ccds),
            is_nrna=bool(row.is_nrna),
            is_synthetic=bool(row.is_synthetic),
            abundance=float(row.abundance) if not np.isnan(row.abundance) else 100.0,
            nrna_abundance=0.0,
        ))
    return out


_READ_NAME_RE = re.compile(r"^([^:]+):(\d+)-(\d+):([^:]+):\d+$")


def _read_origin(qname: str) -> tuple[str, int, int]:
    """Parse origin label from oracle-BAM read name.

    Format: ``{label}:{frag_start}-{frag_end}:{strand}:{idx}``
    where ``label`` is ``<t_id>`` for mRNA, ``nrna_<t_id>`` for nRNA,
    and ``gdna`` for gDNA.
    Returns (kind, start, end) with kind ∈ {mrna, nrna, gdna}.
    """
    m = _READ_NAME_RE.match(qname)
    if m is None:
        return ("unknown", -1, -1)
    label, s, e, _strand = m.groups()
    s, e = int(s), int(e)
    if label == "gdna":
        kind = "gdna"
    elif label.startswith("nrna_"):
        kind = "nrna"
    else:
        kind = "mrna"
    return kind, s, e


def _truth_per_region(bam_path: Path, region_df: pd.DataFrame) -> pd.DataFrame:
    """Count true (mrna, nrna, gdna) fragments per calibration region.

    A fragment is assigned to whichever region contains the mid-point
    of the alignment on the reference — a simple, stable rule that
    matches how region counts are bucketed in the scanner.
    """
    # Build a position→region-id lookup.  Regions tile the reference
    # with non-overlapping intervals.
    if "ref" not in region_df.columns:
        raise ValueError("region_df missing 'ref' column")
    # For a single-ref genome we assume one ref.
    ref = region_df["ref"].iloc[0]
    starts = region_df["start"].to_numpy()
    ends = region_df["end"].to_numpy()
    order = np.argsort(starts)
    starts = starts[order]
    ends = ends[order]
    rid_by_sort = np.arange(len(region_df))[order]

    n_regions = len(region_df)
    cnt = np.zeros((n_regions, 3), dtype=np.int64)  # columns: mrna, nrna, gdna
    kind_col = {"mrna": 0, "nrna": 1, "gdna": 2}

    seen = set()       # dedupe by qname (paired-end pairs appear twice)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped:
                continue
            if rec.query_name in seen:
                continue
            seen.add(rec.query_name)
            kind, s, e = _read_origin(rec.query_name)
            if kind not in kind_col:
                continue
            if rec.reference_name != ref:
                continue
            mid = (s + e) // 2 if s >= 0 else rec.reference_start
            # binary search on `starts`
            idx = int(np.searchsorted(starts, mid, side="right")) - 1
            if idx < 0 or idx >= n_regions:
                continue
            if mid < starts[idx] or mid >= ends[idx]:
                continue
            cnt[rid_by_sort[idx], kind_col[kind]] += 1

    out = pd.DataFrame({
        "n_mrna_true": cnt[:, 0],
        "n_nrna_true": cnt[:, 1],
        "n_gdna_true": cnt[:, 2],
    })
    return out


# ------------------------------------------------------------------ one cell


def run_cell(
    *,
    reference_dir: Path,
    out_dir: Path,
    ss: float,
    gdna_fraction: float,
    n_rna_fragments: int,
    seed: int,
    emit_per_read: bool,
) -> dict:
    label = f"ss{ss:.2f}_gdna{gdna_fraction:.2f}"
    cell_dir = out_dir / "cells" / label
    cell_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()

    # -- simulator inputs ---------------------------------------------------
    index = TranscriptIndex.load(reference_dir / "rigel_index")
    genome = _load_genome(reference_dir / "genome.fa")
    transcripts = _transcripts_from_index(index)

    sim_cfg = SimConfig(
        frag_mean=FRAG_MEAN, frag_std=FRAG_SD,
        frag_min=FRAG_MIN, frag_max=FRAG_MAX,
        read_length=READ_LENGTH,
        strand_specificity=ss, seed=seed,
    )
    gdna_cfg = None
    n_gdna = 0
    if gdna_fraction > 0:
        gdna_cfg = GDNAConfig(
            abundance=1.0,  # placeholder; pool_split is explicit below
            frag_mean=GDNA_FRAG_MEAN, frag_std=GDNA_FRAG_SD,
            frag_min=GDNA_FRAG_MIN, frag_max=GDNA_FRAG_MAX,
        )
        n_gdna = int(round(n_rna_fragments * gdna_fraction))

    sim = OracleBamSimulator(
        genome=genome,
        transcripts=transcripts,
        config=sim_cfg,
        gdna_config=gdna_cfg,
        ref_name=genome.name,
    )

    bam_path = cell_dir / "oracle.bam"
    sim.write_bam(
        bam_path,
        n_fragments=n_rna_fragments + n_gdna,
        name_sorted=True,
        pool_split=(n_rna_fragments, 0, n_gdna),
    )
    t_sim = time.perf_counter() - t0

    # -- scan + calibrate ---------------------------------------------------
    cfg = PipelineConfig(
        scan=BamScanConfig(sj_strand_tag="XS"),
    )
    t1 = time.perf_counter()
    stats, strand_models, fl_models, buffer, region_counts, fl_table = scan_and_buffer(
        str(bam_path), index, cfg.scan,
    )
    strand_models.finalize()
    fl_models.build_scoring_models()
    fl_models.finalize(prior_ess=cfg.calibration.fl_prior_ess)

    cal = calibrate_gdna(
        region_counts,
        fl_table,
        index.region_df,
        strand_models.strand_specificity,
        mean_frag_len=fl_models.global_model.mean,
        intergenic_fl_model=fl_models.intergenic,
        fl_prior_ess=cfg.calibration.fl_prior_ess,
        diagnostics=True,
    )
    t_cal = time.perf_counter() - t1

    # -- truth --------------------------------------------------------------
    genome_len = len(genome)
    lam_truth = n_gdna / genome_len if genome_len > 0 else float("nan")
    truth_regions = _truth_per_region(bam_path, index.region_df)

    # -- per-region frame ---------------------------------------------------
    rdf = index.region_df.reset_index(drop=True).copy()
    per_region = pd.concat([rdf, truth_regions], axis=1)
    per_region["n_total"] = cal.region_n_total.astype(np.int64)
    per_region["gamma"] = (cal.region_gamma
                           if cal.region_gamma is not None
                           else np.full(len(rdf), np.nan))
    per_region["E_gdna"] = cal.region_e_gdna
    per_region["E_rna_expected"] = np.maximum(
        cal.region_n_total - cal.region_e_gdna, 0.0,
    )
    # region class from block truth
    blocks_df = pd.read_csv(reference_dir / "blocks.tsv", sep="\t")
    per_region = _annotate_region_class(per_region, blocks_df)
    per_region.to_feather(cell_dir / "per_region.feather")

    # -- optional per-read dump --------------------------------------------
    if emit_per_read:
        _dump_per_read(bam_path, cell_dir / "per_read.feather",
                       index.region_df, per_region)

    t_total = time.perf_counter() - t0

    row = {
        "label": label,
        "ss_true": ss,
        "gdna_fraction": gdna_fraction,
        "n_rna": n_rna_fragments,
        "n_gdna_truth": n_gdna,
        "genome_length": genome_len,
        "lambda_truth": lam_truth,
        "lambda_est": cal.lambda_gdna,
        "lambda_ratio": (cal.lambda_gdna / lam_truth
                         if lam_truth > 0 else float("nan")),
        "ss_est": cal.strand_specificity,
        "mu_R": cal.mu_R,
        "sigma_R": cal.sigma_R,
        "mixing_pi": cal.mixing_pi,
        "mixing_pi_soft": cal.mixing_pi_soft,
        "em_n_iter": cal.em_n_iter,
        "em_converged": bool(cal.em_converged),
        "n_eligible": cal.n_eligible,
        "strand_used": bool(cal.strand_used),
        "strand_z": cal.strand_z,
        "total_E_gdna": float(cal.region_e_gdna.sum()),
        "total_n_gdna_truth_regions": int(truth_regions["n_gdna_true"].sum()),
        "fl_gdna_mean": (cal.gdna_fl_model.mean
                         if cal.gdna_fl_model is not None else float("nan")),
        "t_sim_s": t_sim,
        "t_cal_s": t_cal,
        "t_total_s": t_total,
    }
    return row


def _annotate_region_class(per_region: pd.DataFrame,
                           blocks_df: pd.DataFrame) -> pd.DataFrame:
    """Attach block kind (U-tx / U-int / R2 / … / R50) to each region.

    The mapping is done by midpoint containment: the region's midpoint
    falls into exactly one block.
    """
    b_start = blocks_df["start"].to_numpy()
    b_end = blocks_df["end"].to_numpy()
    b_kind = blocks_df["kind"].to_numpy()
    b_rep_id = blocks_df["rep_id"].to_numpy()

    order = np.argsort(b_start)
    b_start_s = b_start[order]
    b_end_s = b_end[order]
    b_kind_s = b_kind[order]
    b_rep_s = b_rep_id[order]

    mid = ((per_region["start"].to_numpy() + per_region["end"].to_numpy()) // 2)
    idx = np.searchsorted(b_start_s, mid, side="right") - 1
    ok = (idx >= 0) & (idx < len(b_start_s))
    kinds = np.full(len(per_region), "unknown", dtype=object)
    reps = np.full(len(per_region), np.nan)
    ix = np.where(ok)[0]
    in_block = (mid[ix] >= b_start_s[idx[ix]]) & (mid[ix] < b_end_s[idx[ix]])
    kept = ix[in_block]
    kinds[kept] = b_kind_s[idx[kept]]
    reps[kept] = b_rep_s[idx[kept]]
    per_region["block_kind"] = kinds
    per_region["rep_id"] = reps
    return per_region


def _dump_per_read(bam_path: Path, out_feather: Path,
                   region_df: pd.DataFrame,
                   per_region: pd.DataFrame) -> None:
    """Emit one row per fragment: (truth_kind, mapped_region, gamma)."""
    starts = region_df["start"].to_numpy()
    ends = region_df["end"].to_numpy()
    order = np.argsort(starts)
    starts_s = starts[order]
    ends_s = ends[order]
    rid_by_sort = np.arange(len(region_df))[order]

    gamma = per_region["gamma"].to_numpy()
    kinds = per_region["block_kind"].to_numpy()

    rows: list[tuple] = []
    seen: set[str] = set()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.query_name in seen:
                continue
            seen.add(rec.query_name)
            truth, s, e = _read_origin(rec.query_name)
            mid = (s + e) // 2 if s >= 0 else rec.reference_start
            j = int(np.searchsorted(starts_s, mid, side="right")) - 1
            if j < 0 or j >= len(starts_s) or mid < starts_s[j] or mid >= ends_s[j]:
                rid = -1
                g = float("nan")
                k = "none"
            else:
                rid = int(rid_by_sort[j])
                g = float(gamma[rid])
                k = str(kinds[rid])
            rows.append((rec.query_name, truth, rid, k, g))

    df = pd.DataFrame(rows, columns=["qname", "truth", "region_id",
                                      "block_kind", "gamma"])
    df.to_feather(out_feather)


# ------------------------------------------------------------------ main


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--reference", type=Path, required=True,
                   help="Output of build_reference.py (contains rigel_index/, blocks.tsv)")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--ss", type=float, nargs="+", default=list(DEFAULT_SS_LEVELS))
    p.add_argument("--gdna", type=float, nargs="+", default=list(DEFAULT_GDNA_LEVELS))
    p.add_argument("--n-rna", type=int, default=DEFAULT_N_RNA)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--per-read", action="store_true",
                   help="Emit per-fragment truth/assignment tables per cell.")
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=args.log_level,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )

    args.out.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    for ss, gf in product(args.ss, args.gdna):
        logger.info("Cell ss=%.2f gdna=%.2f", ss, gf)
        try:
            row = run_cell(
                reference_dir=args.reference,
                out_dir=args.out,
                ss=ss, gdna_fraction=gf,
                n_rna_fragments=args.n_rna,
                seed=args.seed,
                emit_per_read=args.per_read,
            )
            rows.append(row)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Cell ss=%.2f gdna=%.2f failed: %s", ss, gf, exc)
            rows.append({
                "label": f"ss{ss:.2f}_gdna{gf:.2f}",
                "ss_true": ss, "gdna_fraction": gf,
                "error": repr(exc),
            })

    df = pd.DataFrame(rows)
    tsv_path = args.out / "sweep_results.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    logger.info("Wrote %s (%d rows)", tsv_path, len(df))
    return 0


if __name__ == "__main__":
    sys.exit(main())
