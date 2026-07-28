"""N1 — the figure: where the enriched mass goes, decomposed by the population that supplies it.

Top row  — three representative conditions. The fitted landscape is EXACTLY a weighted sum of per-node
kernels, so it is drawn as a stacked decomposition by what each training node IS (region/boundary x truly
enriched/depleted), against the oracle. The oracle's split (the truth's own valley) is marked.

Bottom row — the three summary views: the multiplicative attribution of the 0.40 recovery, where each
population's training weight lands, and the arms scoreboard including the control that refutes the
kernel-width idea.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_n1_plot.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_n1_massflow import conds, kernels, truth  # noqa: E402

_EPS = 1e-12
OUT = "/Users/mkiyer/proj/rigel/docs/calibration/figures/gdna_n1_massflow.png"

# Okabe-Ito, in fixed order; region vs boundary additionally carries a hatch, so identity never rests on
# colour alone (the two "depleted" and two "enriched" members are deliberately same-family).
POPS = [
    ("region x depleted", "#0072B2", None),
    ("region x ENRICHED", "#D55E00", None),
    ("boundary x depleted", "#56B4E9", "///"),
    ("boundary x ENRICHED", "#E69F00", "///"),
]
PICKS = ["gdna_gdna300_ss_0.50_nrna_none_capture_on",
         "gdna_gdna100_ss_0.99_nrna_none_capture_on",
         "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong"]


def decompose(s):
    """(split, oracle density, [per-population density contributions]) — the exact stacked identity."""
    orc = L.oracle_landscape(s)
    split = L.two_component(orc)["split"]
    mk = L.masks(s)
    sel = L.recipe_substrate(s, mk)
    tr = truth(s)
    w = np.zeros(sel.shape)
    w[sel] = L.recipe_weights(s, sel, mk)
    w /= max(w.sum(), _EPS)
    enr = tr > split
    parts = []
    for m in (sel & mk["region"] & ~enr, sel & mk["region"] & enr,
              sel & mk["boundary"] & ~enr, sel & mk["boundary"] & enr):
        parts.append((w[m][:, None] * kernels(s["g_hat"][m], s["eff"][m])).sum(0)
                     if m.any() else np.zeros_like(L.GRID))
    return split, orc, parts


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    scen = {s["cond"]: s for s in conds()}
    fig, axes = plt.subplots(2, 3, figsize=(17.5, 9.2))

    # ── top row: the stacked landscape decomposition ────────────────────────────────────────────────
    for ax, cond in zip(axes[0], PICKS):
        s = scen[cond]
        split, orc, parts = decompose(s)
        ax.stackplot(L.GRID, *parts, colors=[p[1] for p in POPS],
                     labels=[p[0] for p in POPS], edgecolor="white", linewidth=0.3)
        for coll, (_, _, h) in zip(ax.collections, POPS):
            if h:
                coll.set_hatch(h)
        ax.plot(L.GRID, orc, color="0.15", lw=2.2, label="ORACLE (region-only)")
        ax.axvline(split, color="0.15", ls=":", lw=1.4)
        # log y: the valley is three orders of magnitude below the depleted spike, and it is the whole point
        ax.set_yscale("log")
        ax.set_ylim(2e-6, 0.4)
        ax.annotate("oracle's split", (split, 0.25), xytext=(4, 0), textcoords="offset points",
                    fontsize=7.5, color="0.3", rotation=90, va="top")
        ax.set_title(f"{' · '.join(s['group'])}", fontsize=9.5)
        ax.set_xlim(-5, 2.5)
        ax.tick_params(labelsize=8)
        ax.grid(alpha=0.25, lw=0.6)
        ax.set_axisbelow(True)
    axes[0][0].legend(fontsize=7.5, loc="lower left", framealpha=0.92)
    axes[0][0].set_ylabel("landscape density (log)", fontsize=9)
    for ax in axes[0]:
        ax.set_xlabel("log$_{10}$ $\\rho_g$", fontsize=9)

    # ── bottom-left: the multiplicative attribution ─────────────────────────────────────────────────
    ax = axes[1][0]
    steps = [("oracle\n(region)", 1.000), ("AMBIG\nexcluded", 0.802), ("boundaries\nincluded", 0.682),
             ("reliability\nweight", 0.687), ("pass-0\ncounts", 1.014)]
    run = 1.0
    xs, bots, hts, cols, tops = [], [], [], [], []
    for i, (name, f) in enumerate(steps):
        if i == 0:
            xs, bots, hts, cols, tops = [name], [0.0], [1.0], ["0.55"], [1.0]
            continue
        new = run * f
        xs.append(name)
        bots.append(min(run, new))
        hts.append(abs(run - new))
        cols.append("#D55E00" if new < run else "#009E73")
        tops.append(new)
        run = new
    xs.append(f"FIT\n{run:.3f}")
    bots.append(0.0)
    hts.append(run)
    cols.append("#0072B2")
    tops.append(run)
    ax.bar(xs, hts, bottom=bots, color=cols, width=0.62)
    for x, b, h, t in zip(xs, bots, hts, tops):
        ax.text(x, b + h + 0.025, f"{t:.3f}", ha="center", fontsize=8)
    ax.set_ylim(0, 1.14)
    ax.set_ylabel("enriched-mass recovery vs oracle", fontsize=9)
    ax.set_title("attribution: 0.802 (AMBIG) × 0.682 (boundaries) × 0.687 (weight) × 1.014 (counts)\n"
                 "= 0.381, against a directly measured 0.399 — the rungs very nearly compose",
                 fontsize=9)
    ax.tick_params(labelsize=7.5)
    ax.grid(axis="y", alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)

    # ── bottom-middle: where each population's weight lands ─────────────────────────────────────────
    ax = axes[1][1]
    lands = np.array([[0.0003, 0.0041, 0.0372], [0.5030, 0.0373, 0.0005],
                      [0.0000, 0.0041, 0.0155], [0.2641, 0.1311, 0.0028]])
    order = [1, 0, 3, 2]  # region-dep, region-enr, boundary-dep, boundary-enr
    labels = ["region\n× depleted", "region\n× ENRICHED", "boundary\n× depleted", "boundary\n× ENRICHED"]
    band_c = ["#0072B2", "#999999", "#D55E00"]
    band_n = ["→ depleted band", "→ the VALLEY", "→ enriched band"]
    bot = np.zeros(4)
    for j in range(3):
        v = lands[order][:, j]
        ax.bar(labels, v, bottom=bot, color=band_c[j], label=band_n[j], width=0.62)
        bot += v
    ax.annotate("0.131 — 74 % of ALL\nvalley mass, from nodes\nthat are truly depleted",
                xy=(2.32, 0.33), xytext=(0.55, 0.40), fontsize=8.5,
                arrowprops=dict(arrowstyle="->", lw=1.1, color="0.25"), color="0.15")
    ax.set_ylabel("share of total training weight", fontsize=9)
    ax.set_title("destination of the training weight (mass is conserved)", fontsize=9.5)
    ax.legend(fontsize=8, loc="upper right")
    ax.tick_params(labelsize=7.5)
    ax.grid(axis="y", alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)

    # ── bottom-right: the arms, on both the headline set and the fabrication guard ──────────────────
    ax = axes[1][2]
    arms = ["BASE", "BND-OUT", "FLAT", "WIDTH", "WIDTH\n+BND-OUT", "ctrl:\nCONST-τ"]
    emd_gdna = [0.3045, 0.4356, 0.2977, 0.3061, 0.1704, 0.1777]
    emd_zero = [0.2059, 0.1190, 0.7392, 0.7777, 0.5491, 0.5354]
    x = np.arange(len(arms))
    ax.bar(x - 0.19, emd_gdna, 0.36, color="#0072B2", label="gDNA-bearing capture-ON (n=13)")
    ax.bar(x + 0.19, emd_zero, 0.36, color="#D55E00", label="ZERO-gDNA guard (n=9)")
    ax.axhline(0.3045, color="#0072B2", ls="--", lw=1.0)
    ax.axhline(0.2059, color="#D55E00", ls="--", lw=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(arms, fontsize=7.5)
    ax.set_ylabel("EMD to the region-only oracle (lower better)", fontsize=9)
    ax.set_title("the arms — and the control that refutes kernel-width", fontsize=9.5)
    ax.legend(fontsize=8, loc="upper left")
    ax.tick_params(labelsize=8)
    ax.grid(axis="y", alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)

    fig.suptitle("N1 — the missing enriched mass is a SUBSTRATE census effect, not a placement error. "
                 "Boundary nodes fill the valley between the two true modes.", fontsize=13)
    fig.tight_layout(rect=(0.005, 0.005, 1, 0.955))
    fig.savefig(OUT, dpi=110)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
