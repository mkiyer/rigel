#!/usr/bin/env python
"""Compare vcap mixture titration results: old pre-Option-B vs new Option-B.

For each library, sum mRNA counts, nRNA counts, and gDNA counts from:
  OLD: /scratch/.../runs/human/<lib>/rigel/{quant,nrna_quant}.feather + summary.json
  NEW: /scratch/.../runs/human_optionb/<lib>/{quant,nrna_quant}.feather + summary.json

Expected behavior with input = 20M RNA + {0,1,2,5,10,20,40,80}M gDNA:
  - mRNA total: flat across titration
  - nRNA total: flat across titration
  - gDNA total: monotonically increasing, ideally ~linear in input gDNA

The key Option-B improvement hypothesis: in mega-loci (where many tx share
genomic coordinates) the old sum-of-spans effective length artificially
suppresses gDNA likelihood. Option B replaces it with a per-fragment
harmonic-mean of hit-local spans, which should shift counts from
mRNA/nRNA → gDNA in contaminated libraries.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

OLD_BASE = Path("/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human")
NEW_BASE = Path("/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human_optionb")

LIBS = [
    "mctp_vcap_rna20m_dna00m",
    "mctp_vcap_rna20m_dna01m",
    "mctp_vcap_rna20m_dna02m",
    "mctp_vcap_rna20m_dna05m",
    "mctp_vcap_rna20m_dna10m",
    "mctp_vcap_rna20m_dna20m",
    "mctp_vcap_rna20m_dna40m",
    "mctp_vcap_rna20m_dna80m",
]


def _dna_m(lib: str) -> int:
    m = re.search(r"dna(\d+)m", lib)
    return int(m.group(1)) if m else -1


def summarise(run_dir: Path) -> dict:
    """Read summary.json; return totals (mRNA/nRNA/gDNA, calibration)."""
    out: dict = {
        "mRNA_count": float("nan"),
        "nRNA_count": float("nan"),
        "gDNA_count": float("nan"),
        "intergenic": float("nan"),
        "cal_E_gdna": float("nan"),
        "cal_lambda_gdna": float("nan"),
        "cal_SS": float("nan"),
        "cal_pi": float("nan"),
    }
    sj = run_dir / "summary.json"
    if not sj.exists():
        return out
    s = json.loads(sj.read_text())
    q = s.get("quantification", {}) or {}
    out["mRNA_count"] = float(q.get("mrna_total", float("nan")))
    out["nRNA_count"] = float(q.get("nrna_total", float("nan")))
    out["gDNA_count"] = float(q.get("gdna_total", float("nan")))
    out["intergenic"] = float(q.get("intergenic_total", float("nan")))
    out["n_transcripts"] = int(q.get("n_transcripts", 0))
    out["n_loci"] = int(q.get("n_loci", 0))
    cal = s.get("calibration", {}) or {}
    out["cal_E_gdna"] = float(cal.get("total_expected_gdna", float("nan")))
    out["cal_lambda_gdna"] = float(cal.get("lambda_gdna", float("nan")))
    out["cal_SS"] = float(cal.get("strand_specificity", float("nan")))
    out["cal_pi"] = float(cal.get("mixing_pi", float("nan")))
    return out


def main() -> None:
    rows = []
    for lib in LIBS:
        dna_m = _dna_m(lib)
        old = summarise(OLD_BASE / lib / "rigel")
        new = summarise(NEW_BASE / lib)
        row: dict = {"lib": lib, "dna_M_input": dna_m}
        for k, v in old.items():
            row[f"old_{k}"] = v
        for k, v in new.items():
            row[f"new_{k}"] = v
        rows.append(row)

    df = pd.DataFrame(rows)

    # Headline comparison
    def fmt(x):
        if pd.isna(x):
            return "  NA   "
        if abs(x) >= 1e6:
            return f"{x / 1e6:6.2f}M"
        if abs(x) >= 1e3:
            return f"{x / 1e3:6.2f}K"
        return f"{x:7.3f}"

    print(f"\nvcap mixture titration: old (pre-Option-B) vs new (Option-B)\n")
    print(
        f"{'lib':24s}  {'dna_in':>6s}  "
        f"{'old_mRNA':>8s} {'new_mRNA':>8s} {'Δ_mRNA_%':>9s}  "
        f"{'old_nRNA':>8s} {'new_nRNA':>8s} {'Δ_nRNA_%':>9s}  "
        f"{'old_gDNA':>8s} {'new_gDNA':>8s} {'Δ_gDNA_%':>9s}"
    )
    print("-" * 135)
    for r in rows:
        if pd.isna(r["new_mRNA_count"]):
            print(f"{r['lib']:24s}  {r['dna_M_input']:>5d}M  (new run not yet complete)")
            continue
        dm = (r["new_mRNA_count"] - r["old_mRNA_count"]) / r["old_mRNA_count"] * 100
        dn = (r["new_nRNA_count"] - r["old_nRNA_count"]) / max(r["old_nRNA_count"], 1) * 100
        dg = (r["new_gDNA_count"] - r["old_gDNA_count"]) / max(r["old_gDNA_count"], 1) * 100
        print(
            f"{r['lib']:24s}  {r['dna_M_input']:>5d}M  "
            f"{fmt(r['old_mRNA_count']):>8s} {fmt(r['new_mRNA_count']):>8s} {dm:+8.2f}%  "
            f"{fmt(r['old_nRNA_count']):>8s} {fmt(r['new_nRNA_count']):>8s} {dn:+8.2f}%  "
            f"{fmt(r['old_gDNA_count']):>8s} {fmt(r['new_gDNA_count']):>8s} {dg:+8.2f}%"
        )

    # Save
    out_tsv = NEW_BASE / "comparison.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\n→ {out_tsv}")


if __name__ == "__main__":
    main()
