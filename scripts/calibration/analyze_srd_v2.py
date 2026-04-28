#!/usr/bin/env python
"""SRD v2 Phase 4 — analyse VCaP mixture sweep results.

Reads summary.json from each library run under
``--root srd_v2_vcap_mixture/default`` and (optionally) the SRD v1
baseline at ``--baseline srd_v1_vcap_mixture/default``, then prints a
side-by-side table of the acceptance metrics from
``docs/calibration/srd_v2_phase2plus_handoff.md`` §5a.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

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


def load(root: Path, lib: str) -> dict | None:
    p = root / lib / "summary.json"
    if not p.exists():
        return None
    return json.loads(p.read_text())


def fl_stats(fl_counts: list[float]) -> tuple[float, float, int]:
    """Return (bin_0_frac, bin_max_frac, mode)."""
    counts = list(fl_counts)
    n = len(counts)
    total = sum(counts)
    if total <= 0:
        return (float("nan"), float("nan"), -1)
    b0 = counts[0] / total
    bm = counts[-1] / total
    mode = max(range(n), key=lambda i: counts[i])
    return (b0, bm, mode)


def fmt(x: float | int | None, kind: str = "f") -> str:
    if x is None:
        return "—"
    if isinstance(x, float):
        if kind == "pct":
            return f"{100 * x:6.2f}%"
        if kind == "f4":
            return f"{x:.4f}"
        if kind == "f2":
            return f"{x:.2f}"
        return f"{x:.3f}"
    return str(x)


def row(lib: str, root: Path, baseline: Path | None) -> dict:
    cur = load(root, lib)
    base = load(baseline, lib) if baseline else None
    out: dict = {"lib": lib}
    if cur is None:
        out["status"] = "missing"
        return out

    cal = cur["calibration"]
    fl = cur["fragment_length"]
    quant = cur["quantification"]

    out["pi_pool"] = cal.get("pi_pool")
    out["quality"] = cal.get("gdna_fl_quality")
    out["n_pool"] = cal.get("n_pool")
    out["n_dropped_oor"] = cal.get("n_pool_dropped_out_of_range", 0)
    out["intronic_pos"] = cal.get("n_pool_intronic_strand_pos", -1)
    out["intronic_neg"] = cal.get("n_pool_intronic_strand_neg", -1)
    out["n_spliced"] = cal.get("n_spliced", -1)
    out["mixture_iters"] = cal.get("mixture_iterations")

    g = fl.get("gdna", {})
    g_hist = g.get("histogram") or g.get("counts")
    if isinstance(g_hist, dict):
        g_hist = g_hist.get("values")
    if g_hist is not None:
        b0, bm, mode = fl_stats(g_hist)
        out["gdna_fl_bin0"] = b0
        out["gdna_fl_binmax"] = bm
        out["gdna_fl_mode"] = mode
    out["gdna_fl_mean"] = cal.get("gdna_fl_mean")

    out["gdna_fraction"] = quant.get("gdna_fraction")
    out["mrna_fraction"] = quant.get("mrna_fraction")
    out["nrna_fraction"] = quant.get("nrna_fraction")

    if base is not None:
        bcal = base["calibration"]
        bfl = base["fragment_length"]
        bq = base["quantification"]
        bg = bfl.get("gdna", {})
        bg_hist = bg.get("histogram") or bg.get("counts")
        if isinstance(bg_hist, dict):
            bg_hist = bg_hist.get("values")
        if bg_hist is not None:
            b0, bm, mode = fl_stats(bg_hist)
            out["b_gdna_fl_bin0"] = b0
            out["b_gdna_fl_binmax"] = bm
            out["b_gdna_fl_mode"] = mode
        out["b_pi_pool"] = bcal.get("pi_pool")
        out["b_gdna_fl_mean"] = bcal.get("gdna_fl_mean")
        out["b_gdna_fraction"] = bq.get("gdna_fraction")
        out["b_mrna_fraction"] = bq.get("mrna_fraction")

    return out


def print_table(rows: list[dict]) -> None:
    cols = [
        ("lib", 32, "s"),
        ("quality", 8, "s"),
        ("pi_pool", 9, "f4"),
        ("gdna_fl_mean", 7, "f2"),
        ("gdna_fl_mode", 5, "d"),
        ("gdna_fl_bin0", 8, "pct"),
        ("gdna_fl_binmax", 8, "pct"),
        ("n_pool", 10, "d"),
        ("n_dropped_oor", 8, "d"),
        ("mrna_fraction", 8, "pct"),
        ("nrna_fraction", 8, "pct"),
        ("gdna_fraction", 8, "pct"),
    ]
    hdr = " ".join(f"{c[0]:>{c[1]}}" for c in cols)
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        if r.get("status") == "missing":
            print(f"{r['lib']:>32}  MISSING")
            continue
        cells = []
        for name, w, k in cols:
            v = r.get(name)
            if v is None:
                cells.append(f"{'—':>{w}}")
                continue
            if k == "s":
                cells.append(f"{v:>{w}}")
            elif k == "d":
                cells.append(f"{int(v):>{w}d}")
            elif k == "pct":
                cells.append(f"{100*float(v):>{w-1}.2f}%")
            elif k == "f4":
                cells.append(f"{float(v):>{w}.4f}")
            elif k == "f2":
                cells.append(f"{float(v):>{w}.2f}")
        print(" ".join(cells))


def print_compare(rows: list[dict]) -> None:
    cols = [
        ("lib", 32),
        ("v1 pi", 8),
        ("v2 pi", 8),
        ("v1 mode", 8),
        ("v2 mode", 8),
        ("v1 mean", 8),
        ("v2 mean", 8),
        ("v1 b0%", 7),
        ("v2 b0%", 7),
        ("v1 bM%", 7),
        ("v2 bM%", 7),
        ("v1 gdnaF", 9),
        ("v2 gdnaF", 9),
    ]
    hdr = " ".join(f"{c[0]:>{c[1]}}" for c in cols)
    print()
    print("==== v2 vs v1 baseline comparison ====")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        if r.get("status") == "missing":
            continue
        def f(v, kind="f"):
            if v is None or (isinstance(v, float) and v != v):
                return "—"
            if kind == "pct":
                return f"{100*v:.2f}"
            if kind == "f4":
                return f"{v:.4f}"
            if kind == "d":
                return str(int(v))
            return f"{v:.2f}"
        cells = [
            f"{r['lib']:>32}",
            f"{f(r.get('b_pi_pool'),'f4'):>8}",
            f"{f(r.get('pi_pool'),'f4'):>8}",
            f"{f(r.get('b_gdna_fl_mode'),'d'):>8}",
            f"{f(r.get('gdna_fl_mode'),'d'):>8}",
            f"{f(r.get('b_gdna_fl_mean')):>8}",
            f"{f(r.get('gdna_fl_mean')):>8}",
            f"{f(r.get('b_gdna_fl_bin0'),'pct'):>7}",
            f"{f(r.get('gdna_fl_bin0'),'pct'):>7}",
            f"{f(r.get('b_gdna_fl_binmax'),'pct'):>7}",
            f"{f(r.get('gdna_fl_binmax'),'pct'):>7}",
            f"{f(r.get('b_gdna_fraction'),'pct'):>9}",
            f"{f(r.get('gdna_fraction'),'pct'):>9}",
        ]
        print(" ".join(cells))


def print_acceptance(rows: list[dict]) -> None:
    print()
    print("==== Acceptance criteria (handoff §5a) ====")
    fails = []
    for r in rows:
        if r.get("status") == "missing":
            continue
        lib = r["lib"]
        b0 = r.get("gdna_fl_bin0")
        bm = r.get("gdna_fl_binmax")
        mode = r.get("gdna_fl_mode")
        n_pool = r.get("n_pool", 0)
        oor = r.get("n_dropped_oor", 0)
        if b0 is not None and b0 >= 0.01:
            fails.append((lib, f"FL bin-0 mass {100*b0:.2f}% >= 1%"))
        if bm is not None and bm >= 0.01:
            fails.append((lib, f"FL bin-max mass {100*bm:.2f}% >= 1%"))
        if mode is not None and not (100 <= mode <= 400):
            fails.append((lib, f"FL mode {mode} outside [100, 400]"))
        if n_pool > 0 and oor / n_pool >= 0.001:
            fails.append((lib, f"OOR fraction {100*oor/n_pool:.3f}% >= 0.1%"))

    # Monotonicity dna00m -> dna80m
    fracs = [(r["lib"], r.get("gdna_fraction")) for r in rows if r.get("status") != "missing"]
    spike_libs = [(L, f) for L, f in fracs if "dna00m" not in L]
    spike_fracs = [f for _, f in spike_libs if f is not None]
    if len(spike_fracs) >= 2:
        non_mono = sum(1 for i in range(1, len(spike_fracs))
                       if spike_fracs[i] < spike_fracs[i-1])
        if non_mono > 0:
            fails.append(("(spike series)", f"gdna_fraction non-monotone in {non_mono} steps"))

    # dna80m headline
    last = next((r for r in rows if r["lib"] == "mctp_vcap_rna20m_dna80m"), None)
    if last and last.get("gdna_fraction") is not None and last["gdna_fraction"] < 0.75:
        fails.append((last["lib"], f"dna80m gdna_fraction {100*last['gdna_fraction']:.2f}% < 75%"))

    if not fails:
        print("ALL PASS")
    else:
        for lib, msg in fails:
            print(f"FAIL  {lib}  {msg}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--baseline", type=Path, default=None)
    ns = ap.parse_args()
    rows = [row(L, ns.root, ns.baseline) for L in LIBS]
    print_table(rows)
    if ns.baseline is not None:
        print_compare(rows)
    print_acceptance(rows)


if __name__ == "__main__":
    main()
