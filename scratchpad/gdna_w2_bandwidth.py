"""W2 — THE KERNEL: sweep the POPULATION RESOLUTION `h_pop`, empirically and with pictures.

The measured defect: the landscape's resolution comes ENTIRELY from each node's own Poisson count, so a
50k-fragment node contributes a ~0.002-decade near-delta. Mode count is therefore a function of TRAINING-SET
SIZE (3.0 -> 5.6 as n_train 1217 -> 121, oracle = 3) — the definition of overfitting. `h_i^2 =
h_within,i^2 + h_pop^2`; `h_pop` is the missing term.

No assumed rule (owner): sweep it. Four questions, four outputs:

  1. WHAT DOES IT LOOK LIKE?      --plot ladder  : the fit at each h_pop vs the oracle, per scenario
  2. WHAT DOES IT COST/BUY?       the metric table + --plot curves
  3. DOES IT FIX THE OVERFITTING? the dose-response gate (mode count vs n_train must go FLAT)
  4. CAN WE CHOOSE IT WITHOUT AN ORACLE?  held-out predictive likelihood — the criterion that also
     exists on real data, where EMD/mode-recall do not.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_w2_bandwidth.py [--plot]
"""
from __future__ import annotations

import sys

import numpy as np
from scipy.special import gammaln

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402

GRID = L.GRID
H_GRID = (0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.45, 0.60)
MODE_TOL = 0.30
ENR_THR = 0.0
_LN10 = np.log(10.0)


# ---------------- oracle-free criterion: held-out predictive likelihood ----------------
def heldout_loglik(s, h_pop, k_folds=5, seed=0):
    """Mean log predictive probability of a HELD-OUT node's count under a landscape fitted without it.

    P(g_j) = SUM_grid P(rho) * Poisson(g_j | rho * E_j)  — a mixture predictive, no constants, and it
    penalises over-smoothing and spurious modes alike. Defined identically on real data (no oracle)."""
    mk = L.masks(s)
    sel = np.flatnonzero(L.recipe_substrate(s, mk))
    if sel.size < 5 * k_folds:
        return np.nan
    rng = np.random.default_rng(seed)
    fold = rng.permutation(sel.size) % k_folds
    g_all, e_all = s["g_hat"][sel], s["eff"][sel]
    out = []
    for f in range(k_folds):
        tr = np.zeros(len(s["g_hat"]), bool)
        tr[sel[fold != f]] = True
        P = L.recipe(s, sel=tr, h_pop=h_pop)
        if P is None:
            continue
        te = fold == f
        lam = np.exp(GRID * _LN10)[None, :] * e_all[te][:, None]          # (n_te, G)
        gg = g_all[te][:, None]
        ll = gg * np.log(np.maximum(lam, 1e-300)) - lam - gammaln(gg + 1.0)
        m = ll.max(1, keepdims=True)
        out.append(np.log(np.maximum((np.exp(ll - m) * P[None, :]).sum(1), 1e-300)) + m[:, 0])
    return float(np.mean(np.concatenate(out))) if out else np.nan


# ---------------- metrics ----------------
def enriched_loc_err(P, orc):
    """|fit - oracle| location of the ENRICHED mode (argmax above ENR_THR). The owner's question: does
    smoothing consolidate fragmented spikes back into one recoverable enriched mode?"""
    hi = GRID > ENR_THR
    if P is None or orc is None or float(orc[hi].sum()) < 1e-3:
        return np.nan
    return abs(float(GRID[hi][np.argmax(P[hi])]) - float(GRID[hi][np.argmax(orc[hi])]))


def row(s, h):
    orc = L.oracle_landscape(s)
    P = L.recipe(s, h_pop=h)
    if P is None or orc is None:
        return None
    mo, mf = L.modes(orc), L.modes(P)
    rec = np.mean([any(abs(a - b) <= MODE_TOL for b, _ in mf) for a, _ in mo]) if mo else np.nan
    spu = np.mean([not any(abs(b - a) <= MODE_TOL for a, _ in mo) for b, _ in mf]) if mf else np.nan
    hi = GRID > ENR_THR
    om = float(orc[hi].sum())
    return (L.emd(P, orc), rec, spu, (float(P[hi].sum()) / om if om > 1e-4 else np.nan),
            enriched_loc_err(P, orc), L.spread(P) / max(L.spread(orc), 1e-9), len(mf))


def sweep():
    for su in ("ambig", "quick"):
        scen = L.load_scenarios(su)
        print(f"=== {su} — population resolution sweep ===")
        print(f"{'h_pop':>6s} {'EMD':>7s} {'recall':>7s} {'spurious':>9s} {'enr mass':>9s}"
              f" {'enr loc err':>12s} {'sd/oracle':>10s} {'n_modes':>8s} {'heldout LL':>11s}")
        for h in H_GRID:
            rows = [r for r in (row(s, h) for s in scen) if r is not None]
            a = np.array(rows, dtype=float)
            hl = np.nanmean([heldout_loglik(s, h) for s in scen])
            print(f"{h:6.2f} {np.nanmean(a[:, 0]):7.3f} {np.nanmean(a[:, 1]):7.2f} "
                  f"{np.nanmean(a[:, 2]):9.2f} {np.nanmean(a[:, 3]):9.2f} {np.nanmean(a[:, 4]):12.3f} "
                  f"{np.nanmean(a[:, 5]):10.2f} {np.nanmean(a[:, 6]):8.1f} {hl:11.4f}")
        print()


def dose_response():
    """THE OVERFITTING GATE: mode count must stop depending on training-set size."""
    print("=== the overfitting gate: mode count vs TRAINING-SET SIZE (want FLAT) ===")
    rng = np.random.default_rng(0)
    scen = {s["cond"]: s for s in L.load_scenarios("ambig")}
    conds = ["gdna_gdna300_ss_0.99_nrna_none_capture_on", "gdna_none_ss_0.50_nrna_none_capture_on"]
    print(f"{'condition':44s} {'h_pop':>6s} " + " ".join(f"{f'n={f:.2f}':>8s}" for f in (1.0, .5, .25, .1))
          + f" {'slope':>7s}")
    for cond in conds:
        s = scen[cond]
        mk = L.masks(s)
        sel = np.flatnonzero(L.recipe_substrate(s, mk))
        n_or = len(L.modes(L.oracle_landscape(s)))
        for h in (0.0, 0.15, 0.30):
            means = []
            for frac in (1.0, 0.5, 0.25, 0.10):
                c = []
                for _ in range(5):
                    idx = rng.choice(sel, max(int(sel.size * frac), 10), replace=False)
                    m = np.zeros(len(s["g_hat"]), bool)
                    m[idx] = True
                    c.append(len(L.modes(L.recipe(s, sel=m, h_pop=h))))
                means.append(float(np.mean(c)))
            slope = means[-1] - means[0]     # modes gained by shrinking the sample 10x
            print(f"{cond:44s} {h:6.2f} " + " ".join(f"{m:8.1f}" for m in means)
                  + f" {slope:+7.1f}   (oracle={n_or})")
        print()


def plot_ladder(out):
    """WHAT BANDWIDTH DOES, visually: the fit at each h_pop against the oracle."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    scen = {s["cond"]: s for s in L.load_scenarios("ambig")}
    picks = [("gdna_gdna300_ss_0.50_nrna_none_capture_on", "enriched mode, unstranded"),
             ("gdna_gdna100_ss_0.50_nrna_present_capture_verystrong", "verystrong: 11 fitted vs 8 oracle"),
             ("gdna_gdna1_ss_0.50_nrna_present_capture_verystrong", "the -0.90 fragmentation case"),
             ("gdna_gdna300_ss_0.50_nrna_none_capture_off", "sharp depleted spike (capOFF)"),
             ("gdna_none_ss_0.50_nrna_none_capture_on", "zero-gDNA: must NOT invent a mode"),
             ("gdna_gdna300_ss_0.99_nrna_none_capture_on", "stranded control")]
    arms = [("today (no resolution term)", dict(), "C3"),
            ("GLOBAL h_pop=0.30  (rejected)", dict(h_pop=0.30), "C7"),
            ("adaptive kNN 0.5", dict(knn_scale=0.5), "C0"),
            ("adaptive kNN 1.0", dict(knn_scale=1.0), "C2"),
            ("adaptive kNN 2.0", dict(knn_scale=2.0), "C1")]
    fig, axes = plt.subplots(2, 3, figsize=(21, 9))
    for ax, (cond, label) in zip(axes.ravel(), picks):
        s = scen[cond]
        orc = L.oracle_landscape(s)
        ax.fill_between(GRID, 0, orc, color="0.78", lw=0, label="ORACLE")
        for name, kw, c in arms:
            ax.plot(GRID, L.recipe(s, **kw), color=c, lw=1.5,
                    ls="--" if "GLOBAL" in name else "-", label=name)
        for loc, _h in L.modes(orc):
            ax.axvline(loc, color="0.45", ls=":", lw=0.9)
        ax.set_title(f"{label}\n{cond}", fontsize=9)
        ax.set_xlim(-5, 2.0)
        ax.tick_params(labelsize=8)
    axes.ravel()[0].legend(fontsize=8)
    fig.suptitle("W2 — the POPULATION RESOLUTION term: GLOBAL (rejected) vs ADAPTIVE nearest-neighbour "
                 "spacing (grey = oracle, dotted = oracle modes)", fontsize=13)
    fig.supxlabel("log10 rho_g   (gDNA density)")
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")


def plot_curves(out):
    """The metric-vs-h_pop curves, per suite — where the optimum is, and whether the criteria agree."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 4, figsize=(22, 8.5))
    for r, su in enumerate(("ambig", "quick")):
        scen = L.load_scenarios(su)
        M = {h: np.array([x for x in (row(s, h) for s in scen) if x is not None], float) for h in H_GRID}
        hl = [np.nanmean([heldout_loglik(s, h) for s in scen]) for h in H_GRID]
        series = [("EMD to oracle (lower better)", 0), ("spurious-mode rate (lower better)", 2),
                  ("enriched-mode location err", 4)]
        for c, (title, col) in enumerate(series):
            y = [np.nanmean(M[h][:, col]) for h in H_GRID]
            ax = axes[r, c]
            ax.plot(H_GRID, y, "o-", color="C0")
            ax.axvline(H_GRID[int(np.nanargmin(y))], color="C3", ls="--", lw=1,
                       label=f"best h={H_GRID[int(np.nanargmin(y))]:.2f}")
            ax.set_title(f"{su}: {title}", fontsize=9)
            ax.set_xlabel("h_pop (decades)")
            ax.legend(fontsize=8)
        ax = axes[r, 3]
        ax.plot(H_GRID, hl, "o-", color="C2")
        ax.axvline(H_GRID[int(np.nanargmax(hl))], color="C3", ls="--", lw=1,
                   label=f"best h={H_GRID[int(np.nanargmax(hl))]:.2f}")
        ax.set_title(f"{su}: HELD-OUT log-likelihood (higher better)\nthe ORACLE-FREE criterion", fontsize=9)
        ax.set_xlabel("h_pop (decades)")
        ax.legend(fontsize=8)
    fig.suptitle("W2 — choosing h_pop: oracle metrics vs the oracle-free held-out criterion", fontsize=13)
    fig.tight_layout(rect=(0, 0.02, 1, 0.94))
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")


if __name__ == "__main__":
    sweep()
    dose_response()
    if "--plot" in sys.argv:
        plot_ladder("/Users/mkiyer/proj/rigel/docs/calibration/figures/gdna_w2_ladder.png")
        plot_curves("/Users/mkiyer/proj/rigel/docs/calibration/figures/gdna_w2_curves.png")
