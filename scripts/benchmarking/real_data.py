"""Real-data 3-way benchmark: rigel vs salmon vs kallisto.

Unlike the simulation-based benchmarks in ``analysis.py``, real-data
libraries have no ground truth. This module compares tool outputs
against each other (pairwise correlations, MA-style residuals) and
exposes rigel-only diagnostics (mRNA / nRNA / gDNA fractional
decomposition, calibration outputs).

Designed for a small set of libraries where each has the layout::

    <library>/
      rigel/                  (legacy run, optional)
      rigel_v05/              (current rigel run, used here)
        quant.feather
        gene_quant.feather
        summary.json
        loci.feather
      salmon/quant.sf.gz
      kallisto/abundance.tsv.gz
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ----------------------------------------------------------------------
# Loaders
# ----------------------------------------------------------------------


def _load_rigel(rigel_dir: Path) -> pd.DataFrame:
    df = pd.read_feather(rigel_dir / "quant.feather")
    keep = ["transcript_id", "gene_id", "gene_name", "effective_length",
            "is_nrna", "count", "tpm"]
    return df[keep].rename(columns={"count": "rigel_count", "tpm": "rigel_tpm",
                                    "effective_length": "rigel_eff_len"})


def _load_salmon(salmon_dir: Path) -> pd.DataFrame:
    df = pd.read_csv(salmon_dir / "quant.sf.gz", sep="\t")
    return df.rename(columns={
        "Name": "transcript_id",
        "EffectiveLength": "salmon_eff_len",
        "NumReads": "salmon_count",
        "TPM": "salmon_tpm",
    })[["transcript_id", "salmon_eff_len", "salmon_count", "salmon_tpm"]]


def _load_kallisto(kallisto_dir: Path) -> pd.DataFrame:
    p = kallisto_dir / "abundance.tsv.gz"
    if not p.exists():
        p = kallisto_dir / "abundance.tsv"
    df = pd.read_csv(p, sep="\t")
    return df.rename(columns={
        "target_id": "transcript_id",
        "eff_length": "kallisto_eff_len",
        "est_counts": "kallisto_count",
        "tpm": "kallisto_tpm",
    })[["transcript_id", "kallisto_eff_len", "kallisto_count", "kallisto_tpm"]]


def _load_rigel_summary(rigel_dir: Path) -> dict:
    p = rigel_dir / "summary.json"
    if not p.exists():
        return {}
    with open(p) as f:
        return json.load(f)


# ----------------------------------------------------------------------
# Metrics (no scipy)
# ----------------------------------------------------------------------


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    rx = pd.Series(x).rank(method="average").to_numpy()
    ry = pd.Series(y).rank(method="average").to_numpy()
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2 or x.std() == 0 or y.std() == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def pairwise_metrics(a: np.ndarray, b: np.ndarray, *, eps: float = 1.0) -> dict:
    """Tool-vs-tool agreement metrics."""
    pa = _pearson(a, b)
    sa = _spearman(a, b)
    la = np.log2(a + eps)
    lb = np.log2(b + eps)
    lp = _pearson(la, lb)
    ls = _spearman(la, lb)
    # Concordance / disagreement at expressed cutoff
    expr_a = a > 0
    expr_b = b > 0
    both = int((expr_a & expr_b).sum())
    only_a = int((expr_a & ~expr_b).sum())
    only_b = int((~expr_a & expr_b).sum())
    return {
        "pearson_r": pa,
        "spearman_r": sa,
        "log2_pearson_r": lp,
        "log2_spearman_r": ls,
        "n": int(len(a)),
        "expressed_both": both,
        "expressed_only_a": only_a,
        "expressed_only_b": only_b,
        "sum_a": float(a.sum()),
        "sum_b": float(b.sum()),
    }


# ----------------------------------------------------------------------
# Per-library analysis
# ----------------------------------------------------------------------


@dataclass
class LibraryResult:
    label: str
    library_dir: Path
    rigel_dir: Path
    joined: pd.DataFrame
    rigel_summary: dict
    pairwise: dict[tuple[str, str], dict]
    rigel_pool: dict
    coverage: dict


def analyze_library(label: str, library_dir: Path, rigel_subdir: str = "rigel_v05",
                    *, top_n_outliers: int = 25) -> LibraryResult:
    library_dir = Path(library_dir)
    rigel_dir = library_dir / rigel_subdir
    salmon_dir = library_dir / "salmon"
    kallisto_dir = library_dir / "kallisto"

    logger.info(f"[{label}] loading tool outputs from {library_dir}")
    rigel = _load_rigel(rigel_dir)
    salmon = _load_salmon(salmon_dir)
    kallisto = _load_kallisto(kallisto_dir)
    rigel_summary = _load_rigel_summary(rigel_dir)

    # Coverage stats before joining
    coverage = {
        "rigel_n_tx": len(rigel),
        "salmon_n_tx": len(salmon),
        "kallisto_n_tx": len(kallisto),
        "rigel_total_count": float(rigel["rigel_count"].sum()),
        "salmon_total_count": float(salmon["salmon_count"].sum()),
        "kallisto_total_count": float(kallisto["kallisto_count"].sum()),
        "rigel_n_synthetic_nrna": int(rigel["is_nrna"].sum()) if "is_nrna" in rigel.columns else 0,
    }

    # Join: outer merge on transcript_id, keep only annotated transcripts
    # (drop rigel synthetic nRNA rows from agreement metrics — salmon/kallisto don't have them)
    rigel_annot = rigel[~rigel["is_nrna"].astype(bool)].copy() if "is_nrna" in rigel.columns else rigel.copy()
    j = rigel_annot.merge(salmon, on="transcript_id", how="inner") \
                   .merge(kallisto, on="transcript_id", how="inner")
    coverage["common_tx"] = len(j)

    # Pairwise metrics on common annotated transcripts
    pairwise = {
        ("rigel", "salmon"):    pairwise_metrics(j["rigel_count"].to_numpy(),
                                                  j["salmon_count"].to_numpy()),
        ("rigel", "kallisto"):  pairwise_metrics(j["rigel_count"].to_numpy(),
                                                  j["kallisto_count"].to_numpy()),
        ("salmon", "kallisto"): pairwise_metrics(j["salmon_count"].to_numpy(),
                                                  j["kallisto_count"].to_numpy()),
    }

    # rigel-only fragment-pool decomposition
    pool = {}
    if rigel_summary:
        cal = rigel_summary.get("calibration", {}) or {}
        em = rigel_summary.get("em_summary", {}) or rigel_summary.get("locus_em", {}) or {}
        # Try to surface useful keys without assuming full schema
        for key in ("pi_pool", "gdna_quality", "gdna_fl_obs", "rna_fl_obs",
                    "global_fl_obs", "category_counts", "ss_inferred",
                    "n_pool", "iter_count", "converged"):
            if key in cal:
                pool[f"cal.{key}"] = cal[key]
        # Emergent totals
        if "is_nrna" in rigel.columns:
            mrna_total = float(rigel.loc[~rigel["is_nrna"].astype(bool), "rigel_count"].sum())
            nrna_total = float(rigel.loc[rigel["is_nrna"].astype(bool), "rigel_count"].sum())
        else:
            mrna_total = float(rigel["rigel_count"].sum())
            nrna_total = 0.0
        gdna_em = rigel_summary.get("gdna_em_total")
        if gdna_em is None:
            gdna_em = rigel_summary.get("locus_em", {}).get("gdna_em_total")
        pool["rigel_mrna_total"] = mrna_total
        pool["rigel_nrna_total"] = nrna_total
        pool["rigel_gdna_em_total"] = float(gdna_em) if gdna_em is not None else None
        align = rigel_summary.get("alignment_stats", {})
        pool["rigel_total_reads"] = align.get("total_reads")
        pool["rigel_mapped_reads"] = align.get("mapped_reads")
        # Implied total accounted-for by rigel
        accounted = mrna_total + nrna_total + (float(gdna_em) if gdna_em else 0.0)
        pool["rigel_accounted"] = accounted

    return LibraryResult(
        label=label,
        library_dir=library_dir,
        rigel_dir=rigel_dir,
        joined=j,
        rigel_summary=rigel_summary,
        pairwise=pairwise,
        rigel_pool=pool,
        coverage=coverage,
    )


# ----------------------------------------------------------------------
# Outlier identification
# ----------------------------------------------------------------------


def top_outliers(j: pd.DataFrame, *, eps: float = 1.0, n: int = 25) -> dict[str, pd.DataFrame]:
    """For each tool pair, compute log2 fold-change residual and surface
    largest disagreements (where total signal >= 50)."""
    out = {}
    floor = 50.0
    for a, b in (("rigel", "salmon"), ("rigel", "kallisto"), ("salmon", "kallisto")):
        ca = j[f"{a}_count"].to_numpy()
        cb = j[f"{b}_count"].to_numpy()
        m = (ca + cb) >= floor
        sub = j.loc[m, ["transcript_id", "gene_name", f"{a}_count", f"{b}_count",
                        f"{a}_eff_len", f"{b}_eff_len"]].copy()
        sub[f"log2fc_{a}_vs_{b}"] = np.log2((ca[m] + eps) / (cb[m] + eps))
        sub["mean_count"] = (ca[m] + cb[m]) / 2.0
        sub = sub.reindex(sub[f"log2fc_{a}_vs_{b}"].abs().sort_values(ascending=False).index)
        out[f"{a}_vs_{b}"] = sub.head(n).reset_index(drop=True)
    return out


# ----------------------------------------------------------------------
# Markdown report writer
# ----------------------------------------------------------------------


def _fmt(v, prec=4):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "—"
    if isinstance(v, float):
        if abs(v) >= 1000:
            return f"{v:,.0f}"
        return f"{v:.{prec}f}"
    if isinstance(v, int):
        return f"{v:,}"
    return str(v)


def write_report(results: list[LibraryResult], out_path: Path, *, title: str) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = [f"# {title}", ""]
    lines.append("Real-data 3-way comparison of rigel, salmon, and kallisto on")
    lines.append("hybrid-capture libraries with varying gDNA spike-in. **No ground")
    lines.append("truth available** — agreement is measured pairwise.")
    lines.append("")

    # ------------------------------------------------------------------
    # §1 Per-library coverage
    # ------------------------------------------------------------------
    lines.append("## §1 Tool coverage and total counts")
    lines.append("")
    cov_rows = []
    for r in results:
        c = r.coverage
        cov_rows.append({
            "library": r.label,
            "rigel_tx": c["rigel_n_tx"],
            "salmon_tx": c["salmon_n_tx"],
            "kallisto_tx": c["kallisto_n_tx"],
            "common_tx": c["common_tx"],
            "rigel_count": c["rigel_total_count"],
            "salmon_count": c["salmon_total_count"],
            "kallisto_count": c["kallisto_total_count"],
            "rigel_synth_nrna": c["rigel_n_synthetic_nrna"],
        })
    lines.append("| library | rigel tx | salmon tx | kallisto tx | common | rigel Σcnt | salmon Σcnt | kallisto Σcnt | rigel synth nRNA |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in cov_rows:
        lines.append("| " + " | ".join([
            row["library"],
            _fmt(row["rigel_tx"]),
            _fmt(row["salmon_tx"]),
            _fmt(row["kallisto_tx"]),
            _fmt(row["common_tx"]),
            _fmt(row["rigel_count"]),
            _fmt(row["salmon_count"]),
            _fmt(row["kallisto_count"]),
            _fmt(row["rigel_synth_nrna"]),
        ]) + " |")
    lines.append("")

    # ------------------------------------------------------------------
    # §2 Pairwise agreement (count- and log-space)
    # ------------------------------------------------------------------
    lines.append("## §2 Pairwise tool agreement (annotated transcripts only)")
    lines.append("")
    for r in results:
        lines.append(f"### {r.label}")
        lines.append("")
        lines.append("| pair | n | Pearson r | log2 Pearson r | Spearman r | log2 Spearman r | Σa | Σb | both expressed | only a | only b |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for (a, b), m in r.pairwise.items():
            lines.append("| " + " | ".join([
                f"{a} vs {b}",
                _fmt(m["n"]),
                _fmt(m["pearson_r"]),
                _fmt(m["log2_pearson_r"]),
                _fmt(m["spearman_r"]),
                _fmt(m["log2_spearman_r"]),
                _fmt(m["sum_a"]),
                _fmt(m["sum_b"]),
                _fmt(m["expressed_both"]),
                _fmt(m["expressed_only_a"]),
                _fmt(m["expressed_only_b"]),
            ]) + " |")
        lines.append("")

    # ------------------------------------------------------------------
    # §3 rigel calibration / pool decomposition
    # ------------------------------------------------------------------
    lines.append("## §3 Rigel pool decomposition (rigel-only diagnostic)")
    lines.append("")
    for r in results:
        lines.append(f"### {r.label}")
        lines.append("")
        p = r.rigel_pool
        cal = r.rigel_summary.get("calibration", {})
        align = r.rigel_summary.get("alignment_stats", {})
        rows = [
            ("BAM total reads", align.get("total_reads")),
            ("BAM mapped reads", align.get("mapped_reads")),
            ("rigel mRNA total (annotated count)", p.get("rigel_mrna_total")),
            ("rigel synthetic-nRNA total", p.get("rigel_nrna_total")),
            ("rigel gDNA EM total", p.get("rigel_gdna_em_total")),
            ("rigel accounted", p.get("rigel_accounted")),
            ("calibration: pi_pool (gDNA fraction)", cal.get("pi_pool")),
            ("calibration: quality", cal.get("gdna_quality") or cal.get("quality")),
            ("calibration: ss_inferred", cal.get("ss_inferred")),
            ("calibration: pool n_obs", cal.get("n_pool")),
            ("calibration: iter / converged", f"{cal.get('iter_count')} / {cal.get('converged')}"),
        ]
        lines.append("| metric | value |")
        lines.append("|---|---:|")
        for k, v in rows:
            lines.append(f"| {k} | {_fmt(v)} |")
        lines.append("")

    # ------------------------------------------------------------------
    # §4 Effective-length comparison (rigel vs salmon)
    # ------------------------------------------------------------------
    lines.append("## §4 Effective-length agreement (rigel vs salmon)")
    lines.append("")
    for r in results:
        j = r.joined
        if "rigel_eff_len" not in j or "salmon_eff_len" not in j:
            continue
        ratio = j["rigel_eff_len"] / j["salmon_eff_len"].replace(0, np.nan)
        ratio = ratio.dropna()
        ratio = ratio[(ratio > 0) & np.isfinite(ratio)]
        lines.append(f"### {r.label}")
        lines.append("")
        lines.append(f"- n transcripts compared: **{len(ratio):,}**")
        lines.append(f"- median(rigel/salmon eff_len): **{ratio.median():.3f}**")
        lines.append(f"- mean(rigel/salmon eff_len):   **{ratio.mean():.3f}**")
        lines.append(f"- 5th / 95th pct ratio:         **{ratio.quantile(0.05):.3f} / {ratio.quantile(0.95):.3f}**")
        lines.append("")

    # ------------------------------------------------------------------
    # §5 Top outliers
    # ------------------------------------------------------------------
    lines.append("## §5 Largest pairwise disagreements (top 25 by |log2 fold change|, mean count ≥ 50)")
    lines.append("")
    for r in results:
        outl = top_outliers(r.joined)
        lines.append(f"### {r.label}")
        lines.append("")
        for pair_name, df in outl.items():
            if df.empty:
                continue
            lines.append(f"**{pair_name}**")
            lines.append("")
            cols = list(df.columns)
            lines.append("| " + " | ".join(cols) + " |")
            lines.append("|" + "|".join(["---"] * len(cols)) + "|")
            for _, row in df.iterrows():
                vals = []
                for c in cols:
                    v = row[c]
                    if isinstance(v, float):
                        vals.append(f"{v:.2f}" if abs(v) < 1000 else f"{v:,.0f}")
                    else:
                        vals.append(str(v))
                lines.append("| " + " | ".join(vals) + " |")
            lines.append("")
        lines.append("")

    out_path.write_text("\n".join(lines))
    logger.info(f"Report written: {out_path}")


# ----------------------------------------------------------------------
# CLI entry
# ----------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    import argparse

    p = argparse.ArgumentParser(description="Real-data 3-way rigel/salmon/kallisto benchmark")
    p.add_argument("--library", action="append", required=True,
                   help="label=path, can be repeated. e.g. dna00m=/path/to/lib")
    p.add_argument("--rigel-subdir", default="rigel_v05",
                   help="subdirectory under each library that holds rigel quant output")
    p.add_argument("-o", "--output", required=True, help="output markdown report path")
    p.add_argument("--title", default="Rigel vs salmon vs kallisto on hybrid-capture VCaP")
    args = p.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s",
                        datefmt="%H:%M:%S")

    libs = []
    for spec in args.library:
        if "=" not in spec:
            raise SystemExit(f"--library must be label=path, got {spec!r}")
        label, path = spec.split("=", 1)
        libs.append((label, Path(path)))

    results = [analyze_library(label, path, rigel_subdir=args.rigel_subdir)
               for label, path in libs]
    write_report(results, Path(args.output), title=args.title)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
