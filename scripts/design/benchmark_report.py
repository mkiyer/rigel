#!/usr/bin/env python3
"""THE BENCHMARK AS ONE HTML PAGE — every scenario, in COUNTS, never pooled until the end.

⭐⭐⭐ **WHAT THIS IS FOR.** Owner, 2026-08-19: *"report the results across every scenario. Don't cherry
pick scenarios. DON'T POOL SCENARIOS. WE NEED TO SEE EVERY SCENARIO. ONLY POOL AT THE END."* — and
*"you may want to make a script to produce an HTML with the results, helpful to look at results for the
full benchmark and for individual transcripts on the test chromosome."* This renders exactly that: one
row per scenario per arm, the pooled total LAST and clearly marked as the summary of rows already shown.

⛔ **IT ADDS NO MEASUREMENT AND SCORES NOTHING.** The numbers come from `relay_pool_ab.py --out`, which
is the one scorer for this table; duplicating a scorer is how a baseline and a ceiling drift apart. This
file is a RENDERER — hand it a TSV and it draws it.

**THE COLUMNS, and every one of them is a COUNT of FRAGMENTS** (owner: *"the table needs to report
COUNTS (not ratios or fractions)"*)::

    gdna_true / rna_true      the origin-split oracle's own tally
    gdna_est  / rna_est       what calibration says
    net                       signed: est − true. Near zero with a large abs = mass is MISPLACED,
                              not mis-sized
    abs                       summed over OBJECTS — the misplaced mass, the one that cannot cancel
    true_nrna / true_mrna     the RNA truth split, so a gDNA over-call names which RNA it ate

⚠ **`abs` can exceed the library**, because it is summed per object: a fragment counted in the wrong
object contributes twice. That is the point of it, and `net` is printed beside it.

Usage::

    python scripts/design/relay_pool_ab.py --suite S --index I --oracle-cache O --out rows.tsv
    python scripts/design/benchmark_report.py rows.tsv -o benchmark.html
    python scripts/design/benchmark_report.py rows.tsv --transcripts T/scenarios -o benchmark.html
    python scripts/design/benchmark_report.py --self-test
"""

from __future__ import annotations

import argparse
import csv
import html
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

#: the axes a condition name encodes. ⛔ Parsed, never guessed at: a scenario whose name does not carry
#: them is reported under its own name rather than silently bucketed with something else.
AXES = (("gdna", ("g00", "g05", "g50", "g98")), ("ss", ("ss_0.50", "ss_0.99")),
        ("capture", ("capture_off", "capture_on")))


def parse_axes(condition: str) -> dict:
    out = {}
    for key, values in AXES:
        hit = [v for v in values if v in condition]
        out[key] = hit[0] if len(hit) == 1 else "?"
    return out


def read_rows(path: Path) -> list[dict]:
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        for k, v in list(r.items()):
            if k in ("condition", "arm", "axis"):
                continue
            try:
                r[k] = float(v)
            except (TypeError, ValueError):
                pass
    return rows


def pooled(rows: list[dict]) -> dict:
    """⭐ THE ONLY POOLING IN THIS FILE, and it is the LAST row of the table by construction."""
    out = {"condition": "ALL SCENARIOS (pooled — the rows above are the measurement)"}
    for k in ("gdna_est", "gdna_true", "gdna_net", "gdna_abs", "rna_est", "rna_true", "rna_net",
              "rna_abs", "true_gdna", "true_nrna", "true_mrna", "mass"):
        if rows and k in rows[0]:
            out[k] = sum(float(r.get(k, 0.0)) for r in rows)
    return out


def transcript_truth(scenarios: Path) -> dict[str, dict[str, float]]:
    """Per scenario, each transcript's OBSERVED fragment count from the simulator's own truth file.

    ⭐ This is what makes the page useful per TRANSCRIPT rather than only per scenario — the owner's
    *"results for individual transcripts on the test chromosome"*. ⛔ Read, never recomputed.
    """
    out: dict[str, dict[str, float]] = {}
    for cond in sorted(p.name for p in Path(scenarios).iterdir() if p.is_dir()):
        tsv = Path(scenarios) / cond / "truth_abundances.tsv"
        if not tsv.is_file():
            continue
        with open(tsv) as fh:
            per = {}
            for r in csv.DictReader(fh, delimiter="\t"):
                mature = float(r.get("observed_mrna_fragments", 0) or 0)
                nascent = float(r.get("observed_nrna_fragments", 0) or 0)
                per[r["transcript_id"]] = mature + nascent
        out[cond] = per
    return out


def _num(x, digits: int = 0) -> str:
    try:
        return f"{float(x):,.{digits}f}"
    except (TypeError, ValueError):
        return html.escape(str(x))


def _table(rows: list[dict], columns: list[tuple[str, str]], *, pool_row: dict | None = None) -> str:
    head = "".join(f"<th>{html.escape(label)}</th>" for _k, label in columns)
    body = []
    for r in rows:
        cells = []
        for key, _label in columns:
            v = r.get(key, "")
            cls = "num" if isinstance(v, float) else "txt"
            cells.append(f'<td class="{cls}">{_num(v) if isinstance(v, float) else html.escape(str(v))}</td>')
        body.append("<tr>" + "".join(cells) + "</tr>")
    if pool_row is not None:
        cells = []
        for key, _label in columns:
            v = pool_row.get(key, "")
            cells.append(f'<td class="num pool">{_num(v) if isinstance(v, float) else html.escape(str(v))}</td>')
        body.append('<tr class="pool">' + "".join(cells) + "</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"


_CSS = """
:root { --bg:#fff; --fg:#1a1a1a; --line:#d8dbe0; --head:#f3f5f7; --pool:#fff8e1; --warn:#b3261e; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg:#14161a; --fg:#e8eaed; --line:#333842; --head:#1d2026; --pool:#2a2417; --warn:#ff8a80; } }
:root[data-theme="dark"] { --bg:#14161a; --fg:#e8eaed; --line:#333842; --head:#1d2026;
  --pool:#2a2417; --warn:#ff8a80; }
body { background:var(--bg); color:var(--fg); font:14px/1.5 -apple-system,BlinkMacSystemFont,
  "Segoe UI",Helvetica,Arial,sans-serif; margin:0; padding:2rem 1.5rem 4rem; }
h1 { font-size:1.4rem; margin:0 0 .3rem; } h2 { font-size:1.05rem; margin:2.2rem 0 .5rem; }
p.note { color:#6b7280; margin:.2rem 0 1rem; max-width:60rem; }
.scroll { overflow-x:auto; border:1px solid var(--line); border-radius:6px; }
table { border-collapse:collapse; width:100%; font-variant-numeric:tabular-nums; }
th,td { padding:.32rem .6rem; border-bottom:1px solid var(--line); white-space:nowrap; }
th { background:var(--head); text-align:right; font-weight:600; position:sticky; top:0; }
th:first-child,td.txt { text-align:left; }
td.num { text-align:right; }
tr.pool td { background:var(--pool); font-weight:600; border-top:2px solid var(--line); }
code { background:var(--head); padding:.05rem .3rem; border-radius:3px; }
.warn { color:var(--warn); font-weight:600; }
"""


def render(rows: list[dict], per_transcript: dict, axis: str = "ALL") -> str:
    sel = [r for r in rows if r.get("axis", axis) == axis]
    conds = sorted({r["condition"] for r in sel})
    arms = sorted({r["arm"] for r in sel}) if sel and "arm" in sel[0] else [""]
    cols = [("condition", "scenario"), ("arm", "policy"),
            ("gdna_true", "gDNA true"), ("gdna_est", "gDNA est"),
            ("gdna_net", "gDNA net"), ("gdna_abs", "gDNA |err|"),
            ("rna_true", "RNA true"), ("rna_est", "RNA est"),
            ("rna_net", "RNA net"), ("rna_abs", "RNA |err|"),
            ("true_nrna", "nascent true"), ("true_mrna", "mature true")]
    parts = [f"<h1>Calibration benchmark — {len(conds)} scenarios, counts in fragments</h1>",
             '<p class="note">Every scenario is shown. The pooled row is <strong>last</strong> and is a '
             'summary of the rows above it, never a substitute for them. <code>net</code> is signed '
             '(est − true); <code>abs</code> is summed over objects — misplaced mass — and can exceed '
             'the library, which is what makes it the one that cannot cancel.</p>']

    for arm in arms:
        arm_rows = [r for r in sel if r.get("arm", "") == arm]
        if not arm_rows:
            continue
        label = {"off": "SilentPolicy (no message propagation)",
                 "on": "RelayPolicy (message propagation)",
                 "currency": "CurrencyPolicy (the Stage-3 rebuild, under development)",
                 }.get(arm, arm or "single arm")
        parts.append(f"<h2>{html.escape(label)}</h2>")
        parts.append('<div class="scroll">'
                     + _table(sorted(arm_rows, key=lambda r: r["condition"]), cols,
                              pool_row=pooled(arm_rows))
                     + "</div>")

    if per_transcript:
        t_ids = sorted({t for per in per_transcript.values() for t in per})
        parts.append("<h2>Per transcript — observed fragments in each scenario</h2>")
        parts.append('<p class="note">The simulator\'s own truth, read from each scenario\'s '
                     '<code>truth_abundances.tsv</code>. A transcript that gets zero fragments in a '
                     'scenario is shown as 0 rather than omitted — an absent row and a measured zero '
                     'are different facts.</p>')
        t_cols = [("condition", "scenario")] + [(t, t) for t in t_ids]
        t_rows = []
        for cond in sorted(per_transcript):
            row = {"condition": cond}
            row.update({t: float(per_transcript[cond].get(t, 0.0)) for t in t_ids})
            t_rows.append(row)
        parts.append('<div class="scroll">' + _table(t_rows, t_cols) + "</div>")
    return f"<title>Calibration benchmark</title><style>{_CSS}</style>" + "".join(parts)


def self_test() -> int:
    ok = fail = 0

    def check(name, cond):
        nonlocal ok, fail
        if cond:
            ok += 1
        else:
            fail += 1
            print(f"   ⛔ {name}")

    rows = [
        {"condition": "gdna_g00_ss_0.50_nrna_x_capture_off", "arm": "off", "axis": "ALL",
         "gdna_est": 100.0, "gdna_true": 0.0, "gdna_abs": 100.0, "rna_abs": 100.0},
        {"condition": "gdna_g50_ss_0.99_nrna_x_capture_on", "arm": "off", "axis": "ALL",
         "gdna_est": 10.0, "gdna_true": 8.0, "gdna_abs": 2.0, "rna_abs": 2.0},
        {"condition": "gdna_g00_ss_0.50_nrna_x_capture_off", "arm": "on", "axis": "REGION",
         "gdna_est": 1.0, "gdna_true": 0.0, "gdna_abs": 1.0, "rna_abs": 1.0},
    ]
    check("axes are parsed from the name",
          parse_axes("gdna_g50_ss_0.99_nrna_x_capture_on") == {"gdna": "g50", "ss": "ss_0.99",
                                                               "capture": "capture_on"})
    check("an unparseable name is marked, not bucketed", parse_axes("nonsense")["gdna"] == "?")

    page = render(rows, {})
    check("every scenario of the selected axis appears", page.count("gdna_g00_ss_0.50") >= 1
          and "gdna_g50_ss_0.99" in page)
    check("a row on a DIFFERENT axis is not mixed in", page.count('class="pool"') >= 1)
    # ⛔ the pooled row must be the LAST row of its table and must equal the sum of the rows shown
    p = pooled([r for r in rows if r["axis"] == "ALL" and r["arm"] == "off"])
    check("pooling sums exactly the rows above it", p["gdna_abs"] == 102.0)
    body = page[page.index("<tbody>"):]
    check("the pooled row comes last", body.rindex("pool") > body.index("gdna_g50"))
    check("counts are rendered with no decimals and thousands separators", ">100<" in page)

    per = {"c1": {"T1": 5.0}, "c2": {"T2": 7.0}}
    page2 = render(rows, per)
    check("a transcript missing from a scenario renders as 0, not blank",
          "T1" in page2 and "T2" in page2 and page2.count(">0<") >= 2)
    check("the page carries a title and its own styling", "<title>" in page2 and "--bg" in page2)
    check("the page is theme-aware in both directions",
          'prefers-color-scheme: dark' in page2 and '[data-theme="dark"]' in page2)

    print(f"\n   self-test: {ok} passed, {fail} failed")
    return 1 if fail else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rows", nargs="?", type=Path, help="the TSV from `relay_pool_ab.py --out`")
    ap.add_argument("-o", "--out", type=Path, default=Path("benchmark.html"))
    ap.add_argument("--axis", default="ALL", choices=("ALL", "REGION", "BOUNDARY"))
    ap.add_argument("--transcripts", type=Path, default=None,
                    help="a scenarios dir, to add the per-transcript observed-fragment table")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()
    if args.rows is None:
        ap.error("a rows TSV is required (from `relay_pool_ab.py --out`)")
    rows = read_rows(args.rows)
    per = transcript_truth(args.transcripts) if args.transcripts else {}
    args.out.write_text(render(rows, per, args.axis))
    n = len({r["condition"] for r in rows if r.get("axis") == args.axis})
    print(f"→ {args.out}   {n} scenarios on axis {args.axis}"
          + (f", {len(per)} scenarios of per-transcript truth" if per else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
