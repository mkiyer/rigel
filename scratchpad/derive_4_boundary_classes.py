"""DERIVATION STEP 4 — do we handle EXON|EXON boundaries and RETAINED INTRONS correctly?

OWNER'S POINT (2026-07-25). RNA is just RNA: there is SPLICED and UNSPLICED, not "mature" and "nascent" as
different species. And a boundary can be an exon-exon boundary that is ALSO a splice junction — RNA can be
contiguous across it WHILE other RNA splices in or out. Both happen at once. Example, both on +:

    TA+ exons (1000,2000), (9000,10000)      <- splices 2000 -> 9000
    TB+ exons (1000,10000)                   <- reads straight through (retained intron)

At position 2000 TA splices out while TB is contiguous. Region 2000-9000 is TA's INTRON and TB's EXON
simultaneously. The 4-bit signature represents this exactly (BIT_EXON_POS | BIT_INTRON_POS both set), but
`coarse_type_array` collapses it ("exon wins over intron") and `bp_solver.is_exon_node` uses that coarse type
- so such a region is treated as a plain exon and RECEIVES THE GRAFT of junction flux that splices OVER it.

This script classifies every region and boundary by its ACTUAL signature bits and measures the error per
class, so we know which of these cases are real in the suite and which are mis-handled.

    OMP_NUM_THREADS=1 python scratchpad/derive_4_boundary_classes.py
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import (  # noqa: E402
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
)
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9

CONDS = [
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_none_ss_0.50_nrna_none_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_index(index)


def region_class(sig):
    """Per-region class from the RAW bits, keeping the retained-intron case the coarse type destroys."""
    ep = (sig & BIT_EXON_POS) != 0
    en = (sig & BIT_EXON_NEG) != 0
    ip = (sig & BIT_INTRON_POS) != 0
    inn = (sig & BIT_INTRON_NEG) != 0
    # RETAINED INTRON: the same strand carries BOTH an exon and an intron annotation
    retained = (ep & ip) | (en & inn)
    ex, it = ep | en, ip | inn
    out = np.where(retained, "RETAINED", np.where(ex & it, "exon+intron(x-strand)",
                   np.where(ex, "exon", np.where(it, "intron", "intergenic"))))
    return out, retained


for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    fg, mass = np.asarray(cap["f_g"]), np.asarray(cap["mass_global"])
    ms = np.asarray(st.mass_spliced, float)
    kind = np.asarray(chain.kind)
    ridx = np.clip(np.asarray(chain.ref_idx, np.int64), 0, len(ra.signature) - 1)
    sig_node = np.asarray(ra.signature)[ridx]
    rcls, retained = region_class(sig_node)
    left, right = np.asarray(chain.left), np.asarray(chain.right)
    ok = np.isfinite(fo) & (mass > _EPS)
    err = np.abs(fg - fo)
    em = err * mass
    tot = em[ok].sum()

    print(f"\n{'=' * 104}\n{cond}   suite mwae {np.average(err[ok], weights=mass[ok]):.4f}\n{'=' * 104}")
    print("  REGIONS by TRUE signature (the coarse type collapses RETAINED into 'exon'):")
    print(f"    {'class':<26}{'nodes':>7}{'mass':>13}{'err-mass':>12}{'share':>8}{'mwae':>9}")
    isreg = ok & (kind == REGION)
    for c in ("exon", "RETAINED", "intron", "exon+intron(x-strand)", "intergenic"):
        m = isreg & (rcls == c)
        if not m.any():
            continue
        print(f"    {c:<26}{int(m.sum()):>7}{mass[m].sum():>13,.0f}{em[m].sum():>12,.0f}"
              f"{em[m].sum() / tot:>7.1%}{np.average(err[m], weights=mass[m]):>9.4f}")

    # BOUNDARIES by the pair of flanking signatures + whether a junction is present
    print("\n  BOUNDARIES by flanking pair  (JUNC = the boundary itself carries spliced mass):")
    print(f"    {'flank pair':<30}{'junc?':>7}{'nodes':>7}{'mass':>12}{'err-mass':>11}{'share':>8}{'mwae':>9}")
    isb = ok & (kind != REGION)
    pair = np.array(["-"] * len(fg), dtype=object)
    for i in np.flatnonzero(isb):
        lo, ro = int(left[i]), int(right[i])
        a = rcls[lo] if lo >= 0 else "edge"
        b = rcls[ro] if ro >= 0 else "edge"
        pair[i] = " | ".join(sorted([str(a), str(b)]))
    rows = []
    for p in sorted(set(pair[isb])):
        for jl, jm in (("yes", ms > 0), ("no", ms <= 0)):
            m = isb & (pair == p) & jm
            if m.sum() < 1:
                continue
            rows.append((em[m].sum(), p, jl, int(m.sum()), mass[m].sum(),
                         np.average(err[m], weights=mass[m]) if mass[m].sum() > 0 else np.nan))
    for e, p, jl, n, ma, mw in sorted(rows, reverse=True)[:10]:
        print(f"    {p:<30}{jl:>7}{n:>7}{ma:>12,.0f}{e:>11,.0f}{e / tot:>7.1%}{mw:>9.4f}")
