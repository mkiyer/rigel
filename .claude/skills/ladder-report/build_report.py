#!/usr/bin/env python
"""Build the Rigel report card — the Artifact page, and the markdown — from the scored arm jsonl.

The page renders, per ladder condition, the base arm on the simulator's own alignments beside the tagged 0.7.1
release and the oracle-prior ceiling; the same base arm on STAR alignments of the same reads (the STAR rung,
`panel.py align`), with what the aligner left unplaced; the two fragment-length gap panels; calibration against
the oracle calibration; and the VCaP mix's gDNA fraction. The markdown is written by
``scripts/design/quant_accuracy.py --markdown``, which is gated by ``tests/calibration/test_quant_accuracy.py``;
this script imports that same module and reuses its ``_load`` / ``stratum`` / ``gdna_estimate``.
Neither renderer runs the pipeline: both read files the instruments already wrote, so the report is a rendering of
a measurement and never a second measurement.

⛔ The composition figures (how many truth-table rows are annotated transcripts, how many gene rows are real genes)
are DERIVED from a condition's own truth table, never typed in: they change with the index.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
SKILL = Path(__file__).resolve().parent
# the panel as every instrument spells it: `scripts/design/_shared.py`, under the data root `$RIGEL_RUNS`
sys.path.insert(0, str(REPO / "scripts" / "design"))
from _shared import DEFAULT_SUITE, RUNS  # noqa: E402

DEFAULT_ARMS = DEFAULT_SUITE / "arms"
#: the gap panels beside the ladder, each with its own arms
FLGAP = ("flgap_rna_long", "flgap_rna_short")
LEVELS = ("g00", "g05", "g50", "g98")

#: the per-axis fields the page renders. ⛔ `fn_mass`/`fn_n` are deliberately absent: a false negative needs an
#: estimate of exactly zero, which a fractional posterior never is, so the column reads 0 on every condition and
#: would look like a perfect score for something unmeasured (TRAPS: could-the-arm-have-fired).
#: `count_over`/`count_under` decompose the total honestly instead.
AXIS_FIELDS = ("n_expressed", "n_detected", "count_true", "count_abs_err", "count_net_err",
               "count_over", "count_under", "fp_mass", "fp_n", "mard", "spearman")


def load_instrument():
    """Import ``quant_accuracy.py`` by path."""
    path = REPO / "scripts" / "design" / "quant_accuracy.py"
    spec = importlib.util.spec_from_file_location("quant_accuracy", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["quant_accuracy"] = mod
    spec.loader.exec_module(mod)
    return mod


def index_composition(suite: Path, condition: str) -> dict:
    """How much of the truth table is scoreable, derived from the table itself.

    The simulator draws from Rigel's own index, so the truth table carries one row per SYNTHETIC nascent entity
    as well as per annotated transcript. Those rows are zero on BOTH sides — zero truth, and zero estimate because
    the transcript table drops them — so they enter no figure. Reporting the raw row count as a transcript count
    claims thousands of perfectly scored transcripts that were never scored, which is why these four numbers are
    carried to the page.
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


def scenario(qa, rows: dict, c: str) -> dict:
    """One condition of one arm file: its pools, its gDNA fraction and its two scored axes."""
    lib, tx, gn = rows[(c, "library")], rows[(c, "transcript")], rows[(c, "gene")]
    strand, capture = qa.stratum(c)
    cap = "ON" if capture.endswith("ON") else "OFF"
    s = {
        "id": c,
        "level": next(k for k in LEVELS if f"_{k}_" in c),
        "strand": strand, "capture": cap,
        "deferred": strand == "unstranded" and cap == "ON",
        "pools": {
            # ⛔ EM + intergenic: the instrument's own definition, never restated here
            "gdna": {"est": qa.gdna_estimate(lib), "true": lib["gdna_true"]},
            #: SYNTHETIC entities only — the split is `is_synthetic`, never `is_single_exon`
            "nascent": {"est": lib["nrna_est"], "true": lib["nrna_true"]},
            "annotated": {"est": lib["mrna_est"], "true": lib["mrna_true"]},
        },
        "gdna_frac": {"est": lib["gdna_frac_est"], "true": lib["gdna_frac_true"]},
        "tx": {**{k: tx[k] for k in AXIS_FIELDS}, "median_rel_err": tx["median_rel_err"]},
        "gene": {k: gn[k] for k in AXIS_FIELDS},
    }
    if "unaligned_gdna" in lib:
        # an aligned input: what the aligner left unplaced, and the gDNA fraction's truth over what it placed
        s["unaligned"] = {k: lib[f"unaligned_{k}"] for k in ("gdna", "mrna", "nrna")}
        s["gdna_frac"]["true_aligned"] = lib["gdna_frac_true_aligned"]
    return s


def v071_rows(path: Path | None, panel: str) -> dict:
    """The tagged 0.7.1 release, scored the same way on the same libraries, by condition."""
    out: dict = {}
    if path is None or not path.is_file():
        return out
    for line in path.read_text().splitlines():
        r = json.loads(line)
        if r.get("panel_id") == panel and "axis" in r:
            out.setdefault(r["condition"], {})[r["axis"]] = r
    return out


def with_v071(qa, s: dict, v: dict | None) -> dict:
    if v:
        lib = v["library"]
        s["v071"] = {
            "tx_abs_err": v["transcript"]["count_abs_err"], "tx_true": v["transcript"]["count_true"],
            "gene_abs_err": v["gene"]["count_abs_err"], "gene_true": v["gene"]["count_true"],
            # 0.7.1's rows recorded one `n_intergenic`, booked as gDNA as today's `n_intergenic_gdna` is
            "pools": {"gdna": lib["gdna_est"] + lib.get("n_intergenic_gdna", lib.get("n_intergenic", 0.0)),
                      "nascent": lib["nrna_est"], "annotated": lib["mrna_est"]},
            "gdna_frac": lib["gdna_frac_est"],
        }
    return s


def calibration_rows(path: Path) -> list[dict]:
    """``calibration_vs_oracle.py --json``'s rows, one per condition."""
    out = []
    if not path.is_file():
        return out
    for r in json.loads(path.read_text()):
        c = r["condition"]
        strand = "stranded" if "_ss_0.99_" in c else "unstranded"
        cap = "ON" if c.endswith("capture_on") else "OFF"
        ax = r["axes"]
        out.append({"id": c, "level": next(k for k in LEVELS if f"_{k}_" in c), "strand": strand, "capture": cap,
                    "deferred": strand == "unstranded" and cap == "ON",
                    "f_P": r["library_f_gdna_P"], "f_O": r["library_f_gdna_O"],
                    "region_mwae": ax["region"]["mwae"], "boundary_mwae": ax["boundary"]["mwae"],
                    "region_net": ax["region"]["net_err"], "boundary_net": ax["boundary"]["net_err"]})
    return out


def head() -> str:
    try:
        return subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def build(args) -> dict:
    qa = load_instrument()
    arms, sources = args.arms, []

    def arm(path: Path):
        sources.append(str(path))
        return qa._load(path)

    base = arm(arms / "qa_ladder_base.jsonl")
    conds = sorted({c for c, _ax in base})
    v071 = v071_rows(args.v071, "ladder")
    extra = {k: arm(arms / f"qa_ladder_{k}.jsonl") for k in ("oracle",) if (arms / f"qa_ladder_{k}.jsonl").is_file()}
    scenarios = []
    for c in conds:
        s = with_v071(qa, scenario(qa, base, c), v071.get(c))
        for k, rows in extra.items():
            if (c, "transcript") in rows:
                s[k] = {"tx_abs_err": rows[(c, "transcript")]["count_abs_err"],
                        "gene_abs_err": rows[(c, "gene")]["count_abs_err"]}
        scenarios.append(s)

    star, star_meta = [], {}
    star_file = arms / f"qa_ladder_base__{args.star}.jsonl"
    if star_file.is_file():
        rows = arm(star_file)
        star = [scenario(qa, rows, c) for c in sorted({c for c, _ax in rows})]
        record = next((args.suite / c / f"{args.star}.align.json" for c in conds
                       if (args.suite / c / f"{args.star}.align.json").is_file()), None)
        star_meta = {"name": args.star}
        if record is not None:
            r = json.loads(record.read_text())
            star_meta.update({"version": r.get("version", "").splitlines()[0] if r.get("version") else "",
                              "aligner": Path(r.get("aligner", "")).name, "params": r.get("params", "")})

    flgap = []
    for panel in FLGAP:
        path = args.suite.parent / panel / "arms" / f"qa_{panel}_base.jsonl"
        if not path.is_file():
            continue
        rows, v = arm(path), v071_rows(args.v071, panel)
        for c in sorted({c for c, _ax in rows}):
            s = with_v071(qa, scenario(qa, rows, c), v.get(c))
            s["panel"] = panel.removeprefix("flgap_")
            flgap.append(s)

    calibration = {}
    for panel in ("ladder", *FLGAP):
        path = args.cvo / f"cvo_{panel}.json"
        if path.is_file():
            sources.append(str(path))
            calibration[panel] = calibration_rows(path)

    vcap = {}
    if args.vcap is not None and args.vcap.is_file():
        q = json.loads(args.vcap.read_text())["quantification"]
        vcap = {"truth": args.vcap_truth, "v071": args.vcap_v071, "rc": q["fraction_gdna"]}
        sources.append(str(args.vcap))
    if args.v071 is not None and args.v071.is_file():
        sources.append(str(args.v071))

    modes = sorted({str(r.get("assignment_mode")) for r in base.values()})
    if len(modes) > 1:
        raise SystemExit(f"⛔ the arm mixes assignment modes {modes} — one mode per report")
    #: the library depth, DERIVED from the truth's own three pools rather than assumed to be the panel's nominal
    #: 10 M — a panel rebuilt at another depth must not silently rescale every bar
    p0 = scenarios[0]["pools"]
    meta = {
        "arm": "qa_ladder_base", "head": args.head or head(), "n": len(conds), "assignment": modes[0],
        "policy": "transfer", "depth": round(sum(p0[k]["true"] for k in p0)),
        "scored": args.scored or time.strftime(
            "%Y-%m-%d", time.localtime((arms / "qa_ladder_base.jsonl").stat().st_mtime)),
        "baseline": args.baseline, "root": str(RUNS), "sources": sources, "star": star_meta,
    }
    if args.suite.is_dir():
        meta.update(index_composition(args.suite, conds[0]))
    return {"meta": meta, "scenarios": scenarios, "star": star, "flgap": flgap, "calibration": calibration,
            "vcap": vcap}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--arms", type=Path, default=DEFAULT_ARMS, help="the ladder's arm jsonl directory")
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE,
                    help="the ladder's directory: its truth table's composition and the aligned input's records")
    ap.add_argument("--star", default="rustar_recommended", help="the STAR rung's aligned input (panel.py align)")
    ap.add_argument("--cvo", type=Path, default=None, help="the directory of calibration_vs_oracle's cvo_<panel>.json")
    ap.add_argument("--vcap", type=Path, default=None, help="the VCaP mix's summary.json")
    ap.add_argument("--vcap-truth", type=float, default=None, help="the VCaP mix's gDNA fraction by read name")
    ap.add_argument("--vcap-v071", type=float, default=None, help="0.7.1's gDNA fraction on the VCaP mix")
    ap.add_argument("--v071", type=Path, default=None, help="the tagged 0.7.1 release's scored rows (jsonl)")
    ap.add_argument("--baseline", default="v0.7.1, scored 2026-10-05", help="the baseline's label on the page")
    ap.add_argument("--head", default=None, help="the commit the numbers come from (default: HEAD)")
    ap.add_argument("--scored", default=None, help="the date the arms were scored (default: the base arm's mtime)")
    ap.add_argument("--html", type=Path, required=True)
    ap.add_argument("--markdown", type=Path, default=None,
                    help="also write the markdown via quant_accuracy.py --markdown")
    args = ap.parse_args()
    if args.cvo is None:
        args.cvo = args.arms

    payload = build(args)
    tpl = (SKILL / "report_template.html").read_text()
    if "__DATA__" not in tpl:
        raise SystemExit("⛔ report_template.html has no __DATA__ placeholder")
    args.html.parent.mkdir(parents=True, exist_ok=True)
    args.html.write_text(tpl.replace("__DATA__", json.dumps(payload, separators=(",", ":"))))
    print(f"  ⭐ page -> {args.html}  ({args.html.stat().st_size:,} bytes; {payload['meta']['n']} ladder conditions, "
          f"{len(payload['star'])} on {args.star}, {len(payload['flgap'])} gap-panel rows, "
          f"{sum(len(v) for v in payload['calibration'].values())} calibration rows, vcap {payload['vcap'] or 'none'})")

    if args.markdown:
        load_instrument().markdown_report(args.arms / "qa_ladder_base.jsonl", args.markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
