"""Exhaustive single-strand-node stress battery (NO AMBIG). Varied exon sizes; sweep every knob.

Scenario: 6 non-overlapping (+) genes (3 expressed, 3 unexpressed), each with exons of sizes
[300 small, 1200 mid, 5000 large] separated by 3 kb introns. Sweeps κ, capture strength, gDNA amount,
RNA depth, and nascent RNA — one knob at a time off a base config. Reports per-class single-strand
accuracy (small/mid/large exon, intron) so we can see any remaining room for optimization and decide
whether to lock the solver. (Prior-free env: RIGEL_PRIOR_FREE=1.)
"""
from __future__ import annotations
import io
import contextlib
import numpy as np
from scripts.debug.toy_prod import run
from rigel.calibration.signature import BIT_EXON_POS, BIT_EXON_NEG, BIT_INTRON_POS, BIT_INTRON_NEG

SIZES = [300, 1200, 5000]  # small, mid, large exons
INTRON = 3000
GAP = 5000


def build_genes():
    genes = []
    x = GAP
    for k in range(6):
        s = x
        exons = []
        for sz in SIZES:
            exons.append((s, s + sz)); s = s + sz + INTRON
        ab = 100.0 if k % 2 == 0 else 0.0
        genes.append((f"G{k}", "+", exons, ab))
        x = s + GAP
    return genes, x + GAP


GENES, GLEN = build_genes()


def _run1(kappa, cap, cstr, gd, n_rna, nasc):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rdf, solved, truth, tg, tr, sig = run(
            "battery", GENES, kappa=kappa, n_rna=n_rna, gdna_fraction=gd, capture=cap,
            capture_strength=cstr, nascent=nasc, genome_length=GLEN, seed=13)
    is_exon = (sig & (BIT_EXON_POS | BIT_EXON_NEG)) != 0
    is_intron = (sig & (BIT_INTRON_POS | BIT_INTRON_NEG)) != 0
    L = (rdf["end"].to_numpy() - rdf["start"].to_numpy())
    obs = (tg + tr) >= 8
    ok = obs & np.isfinite(solved) & np.isfinite(truth)
    classes = {
        "exS": is_exon & (L < 700),
        "exM": is_exon & (L >= 700) & (L < 2500),
        "exL": is_exon & (L >= 2500),
        "intr": is_intron,
    }
    out = {}
    allerr = []
    for name, m in classes.items():
        mm = m & ok
        if mm.sum() == 0:
            out[name] = None
        else:
            e = solved[mm] - truth[mm]
            out[name] = (int(mm.sum()), float(np.median(np.abs(e))), float(np.mean(e)))
            allerr.extend(np.abs(e).tolist())
    out["ALL"] = float(np.median(allerr)) if allerr else None
    return out


def _row(label, o):
    def f(x):
        return "  -      " if x is None else f"{x[2]:+.2f}({x[0]})"
    allv = "  -  " if o["ALL"] is None else f"{o['ALL']:.3f}"
    print(f"{label:22s} | ALL|e|={allv} | exS {f(o['exS'])} exM {f(o['exM'])} "
          f"exL {f(o['exL'])} intr {f(o['intr'])}")


if __name__ == "__main__":
    import os
    print(f"prior_free={os.environ.get('RIGEL_PRIOR_FREE','1')}  (bias(n) per class; ALL=median|err|)")
    base = dict(kappa=0.99, cap=True, cstr=20.0, gd=0.3, n_rna=80000, nasc=0.0)
    print("== base: κ0.99 cap20 gd0.3 nrna80k nasc0 ==")
    _row("base", _run1(**base))
    print("== sweep κ (cap on) ==")
    for k in (0.5, 0.7, 0.9, 0.99):
        _row(f"  κ={k}", _run1(**{**base, "kappa": k}))
    print("== sweep capture strength ==")
    for cap, cs in ((False, 0.0), (True, 5.0), (True, 20.0), (True, 100.0)):
        _row(f"  cap={'off' if not cap else cs}", _run1(**{**base, "cap": cap, "cstr": cs}))
    print("== sweep gDNA amount ==")
    for gd in (0.1, 0.3, 0.6):
        _row(f"  gd={gd}", _run1(**{**base, "gd": gd}))
    print("== sweep RNA depth ==")
    for nr in (20000, 80000, 200000):
        _row(f"  nrna={nr}", _run1(**{**base, "n_rna": nr}))
    print("== sweep nascent RNA (κ0.99 cap on) ==")
    for na in (0.0, 5.0, 20.0):
        _row(f"  nasc={na}", _run1(**{**base, "nasc": na}))
    print("== UNSTRANDED slice (κ0.5) ==")
    for cap, cs in ((False, 0.0), (True, 20.0)):
        _row(f"  κ0.5 cap={'off' if not cap else cs}", _run1(**{**base, "kappa": 0.5, "cap": cap, "cstr": cs}))
