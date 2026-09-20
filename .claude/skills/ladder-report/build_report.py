#!/usr/bin/env python
"""Build the Rigel ladder accuracy report — the markdown and the Artifact page — from arm jsonl.

ONE SET OF LABEL RULES, TWO RENDERINGS. The markdown is written by
``scripts/design/quant_accuracy.py --markdown``, which is gated by
``tests/calibration/test_quant_accuracy.py``; this script imports that same module and reuses its
``_load`` / ``stratum`` / ``_STRATA`` so the page cannot disagree with it about what a field means.
Neither renderer runs the pipeline: both read the arm files ``quant_accuracy`` already wrote, so the
report is a rendering of a measurement and never a second measurement.

⛔ The composition figures (how many truth-table rows are annotated transcripts, how many gene rows are
real genes) are DERIVED from a condition's own truth table, never typed in: they change with the index.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
SKILL = Path(__file__).resolve().parent
DEFAULT_ARMS = Path.home() / "Downloads/rigel_runs/suite/ladder/arms"

#: the per-axis fields the page renders. ⛔ `fn_mass`/`fn_n` are deliberately absent: a false negative
#: needs an estimate of exactly zero, which a fractional posterior never is, so the column reads 0 on
#: every condition and would look like a perfect score for something unmeasured
#: (TRAPS: could-the-arm-have-fired). `count_over`/`count_under` decompose the total honestly instead.
AXIS_FIELDS = ("n_expressed", "n_detected", "count_true", "count_abs_err", "count_net_err",
               "count_over", "count_under", "fp_mass", "fp_n", "mard", "spearman")


def load_instrument():
    """Import ``quant_accuracy.py`` by path so its label rules are shared, not re-implemented."""
    path = REPO / "scripts" / "design" / "quant_accuracy.py"
    spec = importlib.util.spec_from_file_location("quant_accuracy", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["quant_accuracy"] = mod
    spec.loader.exec_module(mod)
    return mod


def index_composition(suite: Path, condition: str) -> dict:
    """How much of the truth table is scoreable, derived from the table itself.

    The simulator draws from Rigel's own index, so the truth table carries one row per SYNTHETIC
    nascent entity as well as per annotated transcript. Those rows are zero on BOTH sides — zero truth,
    and zero estimate because the transcript table drops them — so they enter no figure. Reporting the
    raw row count as a transcript count claims thousands of perfectly scored transcripts that were
    never scored, which is why these four numbers are carried to the page.
    """
    import pandas as pd

    t = pd.read_csv(suite / condition / "truth_abundances.tsv", sep="\t")
    syn = t["transcript_id"].astype(str).str.startswith("RIGEL_NRNA_")
    return {
        "truth_rows": int(len(t)),
        "annotated_tx": int((~syn).sum()),
        "gene_rows": int(t["gene_id"].nunique()),
        "real_genes": int(t.loc[~syn, "gene_id"].nunique()),
    }


def build(arms: Path, arm: str, floor_arm: str, suite: Path | None) -> dict:
    qa = load_instrument()
    primary = qa._load(arms / f"{arm}.jsonl")
    floor_path = arms / f"{floor_arm}.jsonl"
    floor = qa._load(floor_path) if floor_path.is_file() else {}
    conds = sorted({c for c, _ax in primary})

    scenarios = []
    for c in conds:
        lib, tx, gn = primary[(c, "library")], primary[(c, "transcript")], primary[(c, "gene")]
        strand, capture = qa.stratum(c)
        cap = "ON" if capture.endswith("ON") else "OFF"
        f = floor.get((c, "transcript"))
        scenarios.append({
            "id": c,
            "level": next(k for k in ("g00", "g05", "g50", "g98") if f"_{k}_" in c),
            "strand": strand, "capture": cap,
            "deferred": strand == "unstranded" and cap == "ON",
            "pools": {
                # ⛔ EM + intergenic. Intergenic fragments reach no locus and never enter the EM, but
                # they ARE gDNA and the truth counts them; off capture they are more than half of it.
                "gdna": {"est": lib["gdna_est"] + lib["n_intergenic"], "true": lib["gdna_true"]},
                #: SYNTHETIC entities only — the split is `is_synthetic`, never `is_nrna`
                "nascent": {"est": lib["nrna_est"], "true": lib["nrna_true"]},
                "annotated": {"est": lib["mrna_est"], "true": lib["mrna_true"]},
            },
            "gdna_frac": {"est": lib["gdna_frac_est"], "true": lib["gdna_frac_true"]},
            "tx": {**{k: tx[k] for k in AXIS_FIELDS}, "median_rel_err": tx["median_rel_err"]},
            "gene": {k: gn[k] for k in AXIS_FIELDS},
            "floor": abs(tx["count_abs_err"] - f["count_abs_err"]) if f else 0.0,
        })

    modes = sorted({str(r.get("assignment_mode")) for r in primary.values()})
    if len(modes) > 1:
        raise SystemExit(f"⛔ the arm mixes assignment modes {modes} — one mode per report")
    #: the library depth, DERIVED from the truth's own three pools rather than assumed to be the
    #: panel's nominal 10 M — a panel rebuilt at another depth must not silently rescale every bar
    p0 = scenarios[0]["pools"]
    meta = {
        "arm": arm, "n": len(conds), "assignment": modes[0], "policy": "transfer",
        "depth": round(sum(p0[k]["true"] for k in p0)),
        "scored": time.strftime("%Y-%m-%d", time.localtime((arms / f"{arm}.jsonl").stat().st_mtime)),
    }
    if suite is not None:
        meta.update(index_composition(suite, conds[0]))
    return {"meta": meta, "scenarios": scenarios}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--arms", type=Path, default=DEFAULT_ARMS, help="directory of arm jsonl files")
    ap.add_argument("--arm", default="qa_ladder_base")
    ap.add_argument("--floor-arm", default="qa_ladder_base_reseed")
    ap.add_argument("--suite", type=Path, default=Path.home() / "Downloads/rigel_runs/suite/ladder",
                    help="panel directory, read ONLY to derive the truth table's composition")
    ap.add_argument("--html", type=Path, required=True)
    ap.add_argument("--markdown", type=Path, default=None,
                    help="also write the markdown via quant_accuracy.py --markdown")
    args = ap.parse_args()

    payload = build(args.arms, args.arm, args.floor_arm,
                    args.suite if args.suite and args.suite.is_dir() else None)
    tpl = (SKILL / "report_template.html").read_text()
    if "__DATA__" not in tpl:
        raise SystemExit("⛔ report_template.html has no __DATA__ placeholder")
    args.html.parent.mkdir(parents=True, exist_ok=True)
    args.html.write_text(tpl.replace("__DATA__", json.dumps(payload, separators=(",", ":"))))
    print(f"  ⭐ page -> {args.html}  ({args.html.stat().st_size:,} bytes, "
          f"{payload['meta']['n']} conditions)")

    if args.markdown:
        qa = load_instrument()
        paths = [args.arms / f"{args.arm}.jsonl"]
        fp = args.arms / f"{args.floor_arm}.jsonl"
        if fp.is_file():
            paths.append(fp)
        qa.markdown_report(paths, args.markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
