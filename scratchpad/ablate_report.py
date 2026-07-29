"""Score the ablation campaign: accuracy + honest precision + the INTERACTION term.

Protocol per SESSION_2026_07_26_HANDOFF_15 §4.2-4.4:
  * ACCURACY suite-wide (mass-weighted |f_g − oracle|) AND on the FIT SUBSTRATE (REGION & single-strand &
    live), because the substrate is a subset and the two can move in opposite directions;
  * HONEST PRECISION as z2 on a node set held FIXED to the BASELINE arm's confident quartile — never a
    per-arm self-quartile, which re-selects the population and hides the very effect being measured.

    python scratchpad/ablate_report.py base noP1d noM5c ...
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import sys as _sys
_sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from z2 import lin_var  # noqa: E402  -- THE single z2 denominator (log Var -> linear Var)

OUT = Path("/tmp/rigel_ablate")
_EPS = 1e-9
arms = sys.argv[1:] or ["base"]
BASE = arms[0]


def load(arm, refit):
    p = OUT / f"{arm}_r{refit}.npz"
    return np.load(p, allow_pickle=True) if p.exists() else None


def percond(d, m, n_cond=32):
    """Per-condition mass-weighted mwae — the owner's `N better / N worse / N flat` scoring at +-0.002.
    All conditions share one index, so the concatenated arrays reshape cleanly."""
    w, fg, fo = (x.reshape(n_cond, -1) for x in (d["mass"], d["fg"], d["fo"]))
    k = m.reshape(n_cond, -1) & np.isfinite(fo) & (w > _EPS)
    return np.array([np.average(np.abs(fg[i][k[i]] - fo[i][k[i]]), weights=w[i][k[i]])
                     if k[i].any() else np.nan for i in range(n_cond)])


def mwae(d, m):
    w, fg, fo = d["mass"], d["fg"], d["fo"]
    k = m & np.isfinite(fo) & (w > _EPS)
    return float(np.average(np.abs(fg[k] - fo[k]), weights=w[k])) if k.sum() else float("nan")


for refit in (0, 1):
    b = load(BASE, refit)
    if b is None:
        continue
    live = np.isfinite(b["fo"]) & (b["mass"] > _EPS)
    substrate = live & b["reg"].astype(bool) & b["ss"].astype(bool) & b["solv"].astype(bool)
    vb = b["var"]
    q1 = float(np.quantile(vb[np.isfinite(vb) & live], 0.25))
    frozen = live & np.isfinite(vb) & (vb <= q1)  # HELD FIXED to the baseline's confident quartile
    cls, amb = b["cls"], b["amb"].astype(bool)
    pops = [("ALL", frozen), ("exon single", frozen & (cls == "exon") & ~amb),
            ("exon AMBIG", frozen & (cls == "exon") & amb),
            ("boundary single", frozen & (cls == "boundary") & ~amb),
            ("boundary AMBIG", frozen & (cls == "boundary") & amb)]

    print(f"\n{'=' * 108}\nrefit={refit}   nodes live={live.sum():,}  fit-substrate={substrate.sum():,}  "
          f"frozen confident quartile={frozen.sum():,} (var <= {q1:.4g})")
    hdr = (f"  {'arm':<16}{'mwae suite':>12}{'Δ':>9}{'mwae substr':>13}{'Δ':>9}"
           + "".join(f"{p[0][:9]:>11}" for p in pops))
    print(hdr + "\n  " + "-" * (len(hdr) - 2))
    m0 = ms0 = None
    for arm in arms:
        d = load(arm, refit)
        if d is None:
            continue
        e, es = mwae(d, live), mwae(d, substrate)
        if m0 is None:
            m0, ms0 = e, es
        z = []
        for _, m in pops:
            k = m & np.isfinite(d["var"])
            num = float(np.sum(d["mass"][k] * (d["fg"][k] - d["fo"][k]) ** 2))
            den = float(np.sum(d["mass"][k] * lin_var(d["var"][k], d["fg"][k])))
            z.append(num / den if den > 0 else float("nan"))
        bw = ""
        if arm != BASE:
            db = percond(b, live)
            da = percond(d, live)
            better = int((da < db - 0.002).sum())
            worse = int((da > db + 0.002).sum())
            bw = f"   {better} better / {worse} worse / {32 - better - worse} flat"
        print(f"  {arm:<16}{e:>12.4f}{e - m0:>+9.4f}{es:>13.4f}{es - ms0:>+9.4f}"
              + "".join(f"{v:>11.2f}" for v in z) + bw)

    # ── THE INTERACTION CELL ──────────────────────────────────────────────────────────────────────────
    if all(load(x, refit) is not None for x in ("base", "noP1e", "noP1d", "noP1e_noP1d")):
        for nm, sel in (("suite", live), ("fit substrate", substrate)):
            m_0 = mwae(load("base", refit), sel)
            m_a = mwae(load("noP1e", refit), sel)
            m_b = mwae(load("noP1d", refit), sel)
            m_ab = mwae(load("noP1e_noP1d", refit), sel)
            inter = m_ab - m_a - m_b + m_0
            verdict = ("REDUNDANT — they were double-paying; removing both costs LESS than the sum"
                       if inter < -5e-5 else
                       "SYNERGISTIC — each needs the other; removing both costs MORE than the sum"
                       if inter > 5e-5 else "INDEPENDENT — effects add")
            print(f"\n  P1e x P1d interaction on {nm}:  cost(P1e off) {m_a - m_0:+.4f}   "
                  f"cost(P1d off) {m_b - m_0:+.4f}   cost(both off) {m_ab - m_0:+.4f}")
            print(f"    interaction = both − A − B + base = {inter:+.5f}   ⇒ {verdict}")
