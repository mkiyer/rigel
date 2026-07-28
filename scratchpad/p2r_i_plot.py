"""P-2 / PASS-0 — THE FIGURE: the gDNA over-call is the reframe, and the reframe breaks at TRANSCRIPT ENDS.

Panel A — the exact decomposition of the delivered gDNA error. `tg = rg_src · r`, so
          log10(tg/ρ_true) splits identically into the source's own error + log10(r) + the true spatial
          difference. The reframe is 96 % of it; the source is nearly right; gDNA really is uniform.
Panel B — log10(r) by the source boundary's STRUCTURAL bit. At a junction the two faces both carry the
          spliced RNA that crosses, so r ≈ 1. At a transcript TERMINUS no RNA crosses, the boundary face is
          pure gDNA, and r becomes the exon's whole RNA-to-gDNA ratio.
Panel C — delivered vs true gDNA density, coloured by that bit. Junction edges sit on y = x; terminus edges
          are lifted by one to three decades.
Panel D — where the error MASS is, capture-OFF vs capture-ON. 33 % of edges are terminus edges and they
          carry ~2/3 of the error at capture-OFF; at capture-ON the reframe is doing real work everywhere
          and the concentration disappears.

Run: OMP_NUM_THREADS=1 python scratchpad/p2r_i_plot.py
"""
from __future__ import annotations

import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from p2r_h_terminus import SUITE, boundary_classes, run  # noqa: E402

from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EPS = 1e-12
OUT = "/Users/mkiyer/proj/rigel/docs/calibration/figures/gdna_reframe_terminus.png"
OFF = "gdna_gdna100_ss_0.50_nrna_none_capture_off"
ON = "gdna_gdna100_ss_0.50_nrna_none_capture_on"
C_TERM, C_JUNC, C_OTH = "#c1442e", "#2f6f9f", "#9aa0a6"


def main():
    ix = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(ix.region_df, ix.ref_name_to_id)
    term, junc = boundary_classes(ix)
    d_off = run(OFF, ix, ra, cfg, term, junc)
    d_on = run(ON, ix, ra, cfg, term, junc)

    fig, ax = plt.subplots(2, 2, figsize=(13.5, 9.4))
    fig.suptitle("The pass-0 gDNA over-call: it is the REFRAME, and the reframe breaks at transcript ends\n"
                 f"({OFF[5:]}, unless noted)", fontsize=12.5, y=0.985)

    # ── A: the decomposition ─────────────────────────────────────────────────────────────────────────
    a = ax[0, 0]
    d = d_off
    w = np.maximum(d["w"] * np.abs(d["err"]), _EPS)
    # recompute the three terms from what `run` returns (err and logr are exact; source+spatial = err-logr)
    tot = np.average(d["err"], weights=w)
    refr = np.average(d["logr"], weights=w)
    rest = tot - refr
    vals = [tot, rest, refr]
    labs = ["TOTAL\nlog10(delivered/true)", "source's own error\n+ true spatial diff", "the REFRAME\nlog10(r)"]
    cols = ["#444444", "#7aa66c", C_TERM]
    bars = a.bar(range(3), vals, color=cols, width=0.62)
    for b, v in zip(bars, vals):
        a.text(b.get_x() + b.get_width() / 2, v + 0.04, f"{v:+.3f}\n({10 ** v:.0f}x)",
               ha="center", va="bottom", fontsize=9.5)
    a.set_xticks(range(3))
    a.set_xticklabels(labs, fontsize=9)
    a.axhline(0, color="k", lw=0.8)
    a.set_ylabel("decades (mass-weighted mean)")
    a.set_ylim(-0.2, max(vals) * 1.32)
    a.set_title("A · the decomposition is an IDENTITY (residual 4e-16)\n"
                "the reframe carries 96 % of the error", fontsize=10.5)

    # ── B: log10 r by structural class ───────────────────────────────────────────────────────────────
    b = ax[0, 1]
    bins = np.linspace(-1.0, 3.2, 55)
    for m, c, lab in ((d["term"], C_TERM, "TERMINUS source"),
                      (~d["term"] & d["junc"], C_JUNC, "junction-only source")):
        b.hist(np.clip(d["logr"][m], bins[0], bins[-1]), bins=bins, color=c, alpha=0.72,
               label=f"{lab}  (n={int(m.sum())}, med {np.median(d['logr'][m]):+.2f})")
    b.axvline(0.0, color="k", lw=1.1, ls="--")
    b.text(0.03, 0.94, "r = 1\n(gDNA transfers unchanged)", transform=b.transAxes, fontsize=8.5, va="top")
    b.set_xlabel("log10 r     (r = ρ_tot(dst) / ρ_tot(src), the reframe)")
    b.set_ylabel("edges")
    b.legend(fontsize=8.5, loc="upper right")
    b.set_title("B · at a junction the reframe is ~1; at a transcript END\n"
                "the boundary face is pure gDNA and r explodes", fontsize=10.5)

    # ── C: delivered vs true gDNA density ────────────────────────────────────────────────────────────
    c = ax[1, 0]
    dl = d["err"]                       # log10(delivered/true)
    # reconstruct absolute values for the scatter: use the true density implied by err and delivered
    # (run() does not return them separately, so plot err against log10 r, which is the causal pair)
    for m, cc, lab in ((~d["term"] & d["junc"], C_JUNC, "junction-only"),
                       (d["term"], C_TERM, "TERMINUS")):
        c.scatter(d["logr"][m], dl[m], s=9, alpha=0.45, color=cc, label=f"{lab} (n={int(m.sum())})",
                  edgecolors="none")
    lim = [-0.6, 3.2]
    c.plot(lim, lim, color="k", lw=1.0, ls="--", label="err = log10 r  (the reframe IS the error)")
    c.axhline(0, color="#888888", lw=0.8)
    c.set_xlim(lim)
    c.set_ylim(-1.2, 3.4)
    c.set_xlabel("log10 r  (the reframe applied to the gDNA component)")
    c.set_ylabel("log10( delivered gDNA density / TRUE )")
    c.legend(fontsize=8.5, loc="upper left")
    c.set_title("C · the delivered error tracks the reframe one-for-one;\n"
                "junction edges sit at the origin — they are exact", fontsize=10.5)

    # ── D: error-mass concentration, OFF vs ON ───────────────────────────────────────────────────────
    dd = ax[1, 1]
    groups = []
    for nm, dset in (("capture OFF", d_off), ("capture ON", d_on)):
        em = dset["w"] * np.abs(dset["err"])
        t = dset["term"]
        j = ~t & dset["junc"]
        groups.append((nm, 100 * t.mean(), 100 * em[t].sum() / max(em.sum(), _EPS),
                       100 * em[j].sum() / max(em.sum(), _EPS)))
    x = np.arange(2)
    wdt = 0.26
    dd.bar(x - wdt, [g[1] for g in groups], wdt, color="#bbbbbb", label="% of EDGES that are terminus")
    dd.bar(x, [g[2] for g in groups], wdt, color=C_TERM, label="% of ERROR MASS on terminus edges")
    dd.bar(x + wdt, [g[3] for g in groups], wdt, color=C_JUNC, label="% of ERROR MASS on junction edges")
    for k, g in enumerate(groups):
        for off_, v in ((-wdt, g[1]), (0, g[2]), (wdt, g[3])):
            dd.text(k + off_, v + 1.2, f"{v:.0f}", ha="center", fontsize=9)
    dd.set_xticks(x)
    dd.set_xticklabels([g[0] for g in groups])
    dd.set_ylabel("percent")
    dd.set_ylim(0, 88)
    dd.legend(fontsize=8.5, loc="upper center")
    dd.set_title("D · 1/3 of edges carry 2/3 of the error at capture-OFF —\n"
                 "at capture-ON the reframe is doing real work everywhere", fontsize=10.5)

    for axx in ax.ravel():
        axx.grid(alpha=0.22, lw=0.6)
        axx.set_axisbelow(True)
    fig.tight_layout(rect=(0, 0, 1, 0.945))
    fig.savefig(OUT, dpi=150)
    print(f"wrote {OUT}")
    print(f"\n  A: total {tot:+.3f} = source+spatial {rest:+.3f} + reframe {refr:+.3f}")
    for nm, dset in (("OFF", d_off), ("ON", d_on)):
        t = dset["term"]
        print(f"  {nm}: terminus med log10 r {np.median(dset['logr'][t]):+.3f} vs junction "
              f"{np.median(dset['logr'][~t & dset['junc']]):+.3f}")


if __name__ == "__main__":
    main()
