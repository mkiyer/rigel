"""P1e (1) + (2) — the VIOLATION delta and the SURPRISE z2, on the real solver, per condition.

    OMP_NUM_THREADS=1 python scratchpad/p1e_1_violation.py
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import p1e_lib as L  # noqa: E402

ORDER = [
    "gdna300_ss0.99_present_capOFF",
    "gdna300_ss0.50_present_capOFF",
    "gdna300_ss0.99_present_capON",
    "gdna100_ss0.50_present_VERYSTRONG",
    "none_ss0.50_present_capOFF",
    "gdna100_ss0.50_none_capOFF",
    "gdna300_ss0.50_none_capOFF",
]


def line(tag, d, frac):
    print(
        f"    {tag:<26}{d['n']:>7}{d['med']:>10.4f}{d['p25']:>10.4f}{d['p75']:>10.4f}"
        f"{d['p90']:>10.4f}{d['p99']:>11.4f}{100 * frac:>10.1f}%"
    )


def z2line(tag, d, gt1, gt10, gt100):
    print(
        f"    {tag:<26}{d['n']:>7}{d['med']:>10.4f}{d['p25']:>10.4f}{d['p75']:>10.3f}"
        f"{d['p90']:>11.2f}{d['p99']:>12.1f}{100 * gt1:>9.1f}%{100 * gt10:>9.1f}%{100 * gt100:>9.1f}%"
    )


for name in ORDER:
    cond = L.CONDS[name]
    inp, dbg = L.solve(cond)
    t, nf = L.message_table(inp, dbg)
    live = t["nsup"] > 0  # a message that supplies at least one component
    print(f"\n{'=' * 132}")
    print(f"{name}   ({int(live.sum())} live messages of {t['delta'].size} valid edges; "
          f"{int((~live).sum())} supply nothing and conserve exactly by construction)")
    print("=" * 132)

    print("\n  (1) THE VIOLATION  delta = log(M/S)      [S = the claim's asserted fragment count]")
    print(f"    {'stratum':<26}{'n':>7}{'median':>10}{'p25':>10}{'p75':>10}{'p90':>10}{'p99':>11}"
          f"{'|d|>0.1':>11}")
    dl = t["delta"][live]
    line("ALL live", L.q(dl), np.mean(np.abs(dl) > 0.1))
    for cls in ("exon", "boundary", "intron", "intergenic"):
        m = live & (t["cls"] == cls)
        if m.sum() == 0:
            continue
        line(cls, L.q(t["delta"][m]), np.mean(np.abs(t["delta"][m]) > 0.1))
    for tag, m in (("GRAFT edges", live & t["graft"]), ("non-graft", live & ~t["graft"])):
        if m.sum():
            line(tag, L.q(t["delta"][m]), np.mean(np.abs(t["delta"][m]) > 0.1))
    for tag, m in (
        ("exon x GRAFT", live & t["graft"] & (t["cls"] == "exon")),
        ("exon x non-graft", live & ~t["graft"] & (t["cls"] == "exon")),
    ):
        if m.sum():
            line(tag, L.q(t["delta"][m]), np.mean(np.abs(t["delta"][m]) > 0.1))
    # sign: does the claim OVER-state (S > M ⇒ delta < 0)?
    print(f"    over-claim (S>M) share: {100 * np.mean(dl < 0):.1f}%   "
          f"mass-wtd mean claim/obs = "
          f"{np.average(t['S'][live] / np.maximum(t['M'][live], 1e-9), weights=t['M'][live]):.3f}")

    print("\n  (2) THE SURPRISE  z2 = delta^2 / (alpha^T Sigma alpha + 1/n_dst)")
    print(f"    {'stratum':<26}{'n':>7}{'median':>10}{'p25':>10}{'p75':>10}{'p90':>11}{'p99':>12}"
          f"{'z2>1':>10}{'z2>10':>9}{'z2>100':>9}")
    zz = t["z2"][live]
    z2line("ALL live", L.q(zz), np.mean(zz > 1), np.mean(zz > 10), np.mean(zz > 100))
    for cls in ("exon", "boundary", "intron", "intergenic"):
        m = live & (t["cls"] == cls)
        if m.sum() == 0:
            continue
        z = t["z2"][m]
        z2line(cls, L.q(z), np.mean(z > 1), np.mean(z > 10), np.mean(z > 100))
    for tag, m in (("GRAFT edges", live & t["graft"]), ("non-graft", live & ~t["graft"])):
        if m.sum():
            z = t["z2"][m]
            z2line(tag, L.q(z), np.mean(z > 1), np.mean(z > 10), np.mean(z > 100))
    # mass-weighted: what the downstream actually inherits
    w = t["M"][live]
    print(f"    MASS-WEIGHTED mean z2 = {np.average(zz, weights=w):.1f}   "
          f"median-by-mass z2 = "
          f"{np.interp(0.5, np.cumsum(np.sort(w)[np.argsort(np.argsort(zz))]) / w.sum(), np.sort(zz)):.2f}"
          if w.sum() > 0 else "")
    print(f"    mean bhat2_cons (the DL excess) = {np.mean(t['bhat2'][live]):.4f}, "
          f"share with bhat2>0: {100 * np.mean(t['bhat2'][live] > 0):.1f}%")
