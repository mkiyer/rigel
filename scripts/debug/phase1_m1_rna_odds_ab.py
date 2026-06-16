"""Phase-1 Measurement 1 — RNA-odds A/B on the complex-loci battery.

Runs the production complex_loci_benchmark TWICE: once unpatched (RNA-odds ON, the production
default q_rna=0.25) and once patched (RNA-odds OFF — every same-strand-exon edge decoupled by forcing
q_rna=inf, which makes _edge_logphi non-finite ⇒ the O(P) decoupled short-circuit fires everywhere).

calibrate.py does `from .simplex_sweep import deconv_regions_sweep`, so we patch the bound name
`rigel.calibration.calibrate.deconv_regions_sweep` to wrap the original and inject q_rna=inf.

Captures the per-locus table the benchmark prints (stdout) and the TOTAL, for ON vs OFF, plus a
family-prefix rollup. Reproducible:

    python scripts/debug/phase1_m1_rna_odds_ab.py [K=60] [gdna=120] [nrna=25]
"""
import importlib
import importlib.util
import io
import os
import re
import sys
from collections import defaultdict
from contextlib import redirect_stdout

# NOTE: `import rigel.calibration.calibrate as X` binds the *function* calibrate (re-exported in the
# package __init__ — it SHADOWS the submodule). Use importlib.import_module to get the real MODULE so
# patching its `deconv_regions_sweep` global actually reaches calibrate()'s `from .simplex_sweep import`.
cal_mod = importlib.import_module("rigel.calibration.calibrate")
assert hasattr(cal_mod, "deconv_regions_sweep"), "patch target missing — wrong module object"
_orig_sweep = cal_mod.deconv_regions_sweep

_BENCH_PATH = os.path.join(os.path.dirname(__file__), "complex_loci_benchmark.py")
_spec = importlib.util.spec_from_file_location("complex_loci_benchmark", _BENCH_PATH)
bench = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bench)


_FIRE = {"n": 0}


def _patched_off(*args, **kw):
    _FIRE["n"] += 1
    kw["q_rna"] = float("inf")
    return _orig_sweep(*args, **kw)


_ROW = re.compile(r"^\s*([A-Za-z0-9_]+)\s+(\d+)\s+([\d,]+)\s*$")
_TOTAL = re.compile(r"^\s*TOTAL\s+([\d,]+)\s*$")


def _parse(out):
    """Return (per_locus dict name->(n_amb, err), total)."""
    per = {}
    total = None
    for line in out.splitlines():
        mt = _TOTAL.match(line)
        if mt:
            total = float(mt.group(1).replace(",", ""))
            continue
        m = _ROW.match(line)
        if m and m.group(1) not in ("locus",):
            per[m.group(1)] = (int(m.group(2)), float(m.group(3).replace(",", "")))
    return per, total


def _run(argv):
    sys.argv = argv
    buf = io.StringIO()
    with redirect_stdout(buf):
        bench.main()
    out = buf.getvalue()
    return out, _parse(out)


def _family(name):
    # roll up by the leading token before the first underscore (cross/nested/.../multistrand/anchor/...)
    return name.split("_")[0]


def main():
    K = sys.argv[1] if len(sys.argv) > 1 else "60"
    gdna = sys.argv[2] if len(sys.argv) > 2 else "120"
    nrna = sys.argv[3] if len(sys.argv) > 3 else "25"
    argv = ["complex_loci_benchmark.py", K, gdna, nrna]

    print(f"=== M1: RNA-odds A/B (K={K} gDNA={gdna} nRNA={nrna}) ===\n")

    print(">>> RUN A: RNA-odds ON (production default q_rna=0.25)")
    out_on, (per_on, tot_on) = _run(list(argv))

    print(">>> RUN B: RNA-odds OFF (q_rna=inf — all edges decoupled)")
    cal_mod.deconv_regions_sweep = _patched_off
    try:
        out_off, (per_off, tot_off) = _run(list(argv))
    finally:
        cal_mod.deconv_regions_sweep = _orig_sweep
    print(f"    (OFF patch fired {_FIRE['n']} sweep calls — must be > 0)")

    # per-locus table (only loci with measurable AMBIG mass in either run)
    print("\n=== PER-LOCUS leak (mass-weighted |f_g-oracle|), ON vs OFF ===")
    print(f"{'locus':>34} {'#amb':>4} {'ON':>11} {'OFF':>11} {'OFF-ON':>11}")
    names = sorted(set(per_on) | set(per_off))
    for nm in names:
        non, eon = per_on.get(nm, (0, 0.0))
        noff, eoff = per_off.get(nm, (0, 0.0))
        nshow = non or noff
        delta = eoff - eon
        flag = ""
        if abs(delta) > 1.0:
            flag = "  ON wins" if delta > 0 else "  OFF wins"
        print(f"{nm:>34} {nshow:>4} {eon:>11,.0f} {eoff:>11,.0f} {delta:>11,.0f}{flag}")

    # family rollup
    fam_on = defaultdict(float)
    fam_off = defaultdict(float)
    for nm, (_, e) in per_on.items():
        fam_on[_family(nm)] += e
    for nm, (_, e) in per_off.items():
        fam_off[_family(nm)] += e
    print("\n=== FAMILY rollup (sum of per-locus leak) ===")
    print(f"{'family':>16} {'ON':>11} {'OFF':>11} {'OFF-ON':>11}")
    for fam in sorted(set(fam_on) | set(fam_off)):
        a, b = fam_on[fam], fam_off[fam]
        print(f"{fam:>16} {a:>11,.0f} {b:>11,.0f} {b - a:>11,.0f}")

    print("\n=== TOTAL ===")
    print(f"  ON  TOTAL = {tot_on:,.0f}")
    print(f"  OFF TOTAL = {tot_off:,.0f}")
    if tot_on is not None and tot_off is not None:
        d = tot_off - tot_on
        verdict = "OFF matches-or-beats ON" if d <= 0 else "ON wins"
        rel = (d / tot_on * 100.0) if tot_on else float("nan")
        print(f"  OFF-ON   = {d:,.0f}  ({rel:+.1f}% vs ON)  ->  {verdict}")


if __name__ == "__main__":
    main()
