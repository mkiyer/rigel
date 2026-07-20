"""Monte-Carlo validation of the cross-node message arithmetic (docs/calibration/cliff_message_derivation.md).

Simulates gDNA and nascent-RNA fragments over an intron + adjacent exon, deposits contained-in-intron vs
crossing the intron|exon boundary, measures the TRUE (ground-truth) f_g at each frame, and checks the derived
LOG-ODDS SHIFT prediction  λ_bnd = λ_intron + [log(E_g^B/E_g^I) − log(E_r^B/E_r^I)]  across many gDNA/RNA
fragment-length distributions, plus the capture-enrichment (cliff) invariance. Run inside the `rigel` env.

The point: gDNA and RNA have different FL distributions ⇒ different per-component effective lengths ⇒ the COUNT
composition f_g is FRAME-dependent (region-contained vs boundary-crossing). The naive "f_g is invariant" is wrong
whenever the FL distributions differ; the eff-length-frame log-odds shift is exact; and the enrichment cancels.
"""

from __future__ import annotations

import numpy as np

from rigel.calibration.effective_length import fl_mean, region_eff_length


def make_pmf(kind: str, m: float, s: float, L: int = 800) -> np.ndarray:
    """A fragment-length pmf over [1, L): gaussian / gamma / bimodal / uniform, mean≈m, spread≈s."""
    x = np.arange(1, L)
    if kind == "gauss":
        p = np.exp(-0.5 * ((x - m) / s) ** 2)
    elif kind == "gamma":
        k = (m / s) ** 2
        th = s * s / m
        p = x ** (k - 1) * np.exp(-x / th)
    elif kind == "bimodal":
        p = np.exp(-0.5 * ((x - m) / s) ** 2) + 0.7 * np.exp(-0.5 * ((x - (m + 180)) / (s * 0.6)) ** 2)
    elif kind == "unif":
        p = ((x >= m - s) & (x <= m + s)).astype(float)
    else:
        raise ValueError(kind)
    p = np.maximum(p, 0.0)
    pmf = np.zeros(L)
    pmf[1:] = p / p.sum()
    return pmf


def mc(gpmf, rpmf, Li=2000, dg=0.4, dr=0.6, e_exon=1.0, Nbase=6_000_000, seed=1):
    """Ground-truth (fg_intron, fg_boundary) from deposited fragments. ``e_exon`` = capture enrichment applied
    to fragments overlapping the exon (positional capture; nucleic-acid-agnostic)."""
    rng = np.random.default_rng(seed)
    maxfl = max(len(gpmf), len(rpmf))

    def run(pmf, dens):
        n = int(Nbase * dens)
        s = rng.uniform(-maxfl, Li, n)
        ell = rng.choice(len(pmf), size=n, p=pmf).astype(float)
        end = s + ell
        cont = (s >= 0) & (end <= Li)  # wholly contained in the intron
        cross = (s <= Li) & (end >= Li)  # crosses the intron|exon boundary point
        w = np.where(end > Li, e_exon, 1.0)  # overlaps the (captured) exon ⇒ enriched, gDNA & RNA alike
        return float((cont * 1.0).sum()), float((cross * w).sum())

    gc, gx = run(gpmf, dg)
    rc, rx = run(rpmf, dr)
    return gc / (gc + rc), gx / (gx + rx)


def predict_boundary_fg(gpmf, rpmf, fg_intron, Li=2000):
    """The derived log-odds shift: λ_bnd = λ_intron + log(E_g^B/E_g^I) − log(E_r^B/E_r^I). Large flanks ⇒ the
    crossing eff-len is ``fl_mean``; the contained eff-len is ``region_eff_length``."""
    eg_i = region_eff_length(np.array([float(Li)]), gpmf)[0]
    er_i = region_eff_length(np.array([float(Li)]), rpmf)[0]
    eg_b, er_b = fl_mean(gpmf), fl_mean(rpmf)
    shift = np.log(eg_b / eg_i) - np.log(er_b / er_i)
    odds = fg_intron / (1.0 - fg_intron)
    return 1.0 / (1.0 + 1.0 / (odds * np.exp(shift))), shift


def main():
    cases = [
        ("gDNA=RNA (gauss 200/40)", make_pmf("gauss", 200, 40), make_pmf("gauss", 200, 40)),
        ("gDNA 300 > RNA 150 (gauss)", make_pmf("gauss", 300, 50), make_pmf("gauss", 150, 30)),
        ("gDNA 120 < RNA 350 (gauss)", make_pmf("gauss", 120, 25), make_pmf("gauss", 350, 60)),
        ("gDNA gamma300 / RNA gamma180", make_pmf("gamma", 300, 90), make_pmf("gamma", 180, 60)),
        ("gDNA bimodal / RNA gauss220", make_pmf("bimodal", 180, 40), make_pmf("gauss", 220, 45)),
        ("gDNA unif 250±100 / RNA unif 150±60", make_pmf("unif", 250, 100), make_pmf("unif", 150, 60)),
    ]
    print(
        f"{'FL distributions':36s} {'fg_intron':>9s} {'fg_bnd(MC)':>10s} "
        f"{'pred(shift)':>11s} {'naive':>7s} {'shift':>7s}"
    )
    for name, g, r in cases:
        fg_i, fg_b = mc(g, r)
        pred, sh = predict_boundary_fg(g, r, fg_i)
        print(f"{name:36s} {fg_i:9.4f} {fg_b:10.4f} {pred:11.4f} {fg_i:7.4f} {sh:7.3f}")
    print("\nenrichment (cliff) invariance — boundary f_g must not move with e_exon:")
    g, r = make_pmf("gauss", 300, 50), make_pmf("gauss", 150, 30)
    for e in (1.0, 30.0, 300.0):
        fg_i, fg_b = mc(g, r, e_exon=e)
        print(f"  e_exon={e:6.0f}:  fg_intron={fg_i:.4f}  fg_boundary(MC)={fg_b:.4f}")


if __name__ == "__main__":
    main()
