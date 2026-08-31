"""HOW GOOD IS THE TOTAL-DENSITY LANDSCAPE, AND WHAT IS ITS BANDWIDTH REALLY?
⭐⭐⭐ **A MEASUREMENT: no solver runs, no EM, and the only thing patched in `src/` is `_N_GRID`, per
arm, restored in a ``finally``.**

⛔⛔ **THIS INSTRUMENT USED TO CARRY AN NPMLE ARM AND IT WAS REMOVED WITH THE NPMLE ITSELF
(2026-08-21). You cannot A/B a deleted thing, and its verdict is not lost — it was MOVED.** The
head-to-head measured, on all 16 ladder conditions: the landscape's depleted level lands 4.8–21×
closer to certified gDNA truth (0.0056 / 0.0264 nats against 0.1191 / 0.1266 on the two capture-OFF
strata); the NPMLE's ``mass / eff_gdna`` axis carries an irreducible per-region offset spread (IQR
0.12 nats off capture, **1.66 under it**, comparable to the mode separation itself) that no bandwidth
removes; and on a generic held-out predictive likelihood the two TIE off capture, which is the finding
worth keeping — the NPMLE was not a bad density estimate, it was estimating the wrong quantity. The
permanent homes are `ISSUES: measured-prior-rung-4` and `DESIGN.md` §3.1a-iii; the rendered comparison is the
Abundance Landscape Atlas.

⭐⭐ **WHY THE REMAINING ARMS ARE STILL THE RIGHT ONES.** The three origin partitions' START banks SUM
TO THE FULL PAYLOAD's (an existing gate, re-asserted on every row here), so the measured total IS the
true total up to Poisson noise — this landscape's inputs ARE the truth. What the oracle partitions add
is the DECOMPOSITION: whether the depleted level the fit reports is the gDNA background it is read as,
and how much RNA contaminates the anchor pool that defines it.

**ⓐ HELD-OUT PREDICTIVE LOG-LIKELIHOOD — the yardstick that needs no chosen constant.**
Split the regions by index parity (deterministic, and it interleaves neighbours so both halves are
the same population), fit on one half, and score the mean ``log ∫ P(c | ρ·e)·P(log ρ) dρ`` of the
OTHER half's observations. ⭐ This is a proper score of a density estimate: it answers "which prior
would have predicted the data it did not see" with no reference distribution to choose and no
smoothing applied to a truth. ⛔ EMD is monotone in smoothing (`landscape`'s own header records it)
and must not select anything; a predictive likelihood is not, which is what makes the sweeps below a
SELECTION rather than a taste adjustment.

**ⓑ THE DEPLETED LEVEL AGAINST CERTIFIED gDNA TRUTH.** The census (the shipped ``split_basins``, one
rule with one home) gives a depleted mode; the truth is the ``gdna`` partition's pooled intergenic
rate through the SAME side selection. Error in nats. ⛔ At a zero-gDNA condition the pool's truth is
exactly 0, so this arm is STAMPED VACUOUS rather than reported — a ratio against zero is not a
measurement.

**② THE BANDWIDTH QUESTION.** ``--sweep`` prices ``_KNN_SCALE`` by ⓑ, per stratum, plus SPLIT-HALF
MODE STABILITY (how many of the modes found on one half reproduce on the other, within their own
fitted widths). ⛔ `TRAPS: no-magic-numbers` — the shipped 0.5 was SHAPE-selected on gDNA-shaped
data, and if a different scale wins here it lands with this evidence attached or not at all.

Usage::

    python scripts/design/landscape_head_to_head.py                     # every cached condition
    python scripts/design/landscape_head_to_head.py --condition NAME    # one condition
    python scripts/design/landscape_head_to_head.py --sweep             # + the _KNN_SCALE sweep
    python scripts/design/landscape_head_to_head.py --json out.json     # rows + both curves
    python scripts/design/landscape_head_to_head.py --self-test         # perturbed, no I/O
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402
from scipy.special import gammaln  # noqa: E402


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, Path(__file__).resolve().parent / name)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


OC = _sibling("object_composition.py")
PVO = OC.PVO

from rigel.calibration.abundance_landscape import (  # noqa: E402
    _census,
    fit_abundance_landscape,
    split_basins,
)
from rigel.calibration.landscape import DensityLandscape  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import RegionType, coarse_type_array  # noqa: E402
from rigel.calibration.splice_graph import (  # noqa: E402
    build_contiguous_boundary_reach_arrays,
    build_mature_wall_distances,
    build_sj_geometry_arrays,
)
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.calibration.total_abundance import (  # noqa: E402
    build_region_wall_mask,
    region_counts_and_exposure,
    w_max_from_deposited_lengths,
)
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import read_scan_cache  # noqa: E402

DEFAULT_SUITE = OC.DEFAULT_SUITE
DEFAULT_INDEX = OC.DEFAULT_INDEX

#: the three origin partitions. The START banks must sum to the full payload's — asserted, not assumed.
ORIGINS = ("gdna", "mrna", "nrna")

#: the `_KNN_SCALE` values `--sweep` prices. ⛔ NOT a constant in any estimator: a sweep axis, and the
#: shipped 0.5 is one of its points so the incumbent is always scored beside the challengers.
_SWEEP_SCALES = (0.125, 0.25, 0.5, 1.0, 2.0)

#: the `_N_GRID` values `--grid-sweep` prices, the shipped 260 among them.
#:
#: ⭐⭐⭐ **WHY THIS SWEEP EXISTS, AND IT IS THE ANSWER TO "IS THE BANDWIDTH OPTIMAL?"** Measured on the
#: ladder: at the shipped `_KNN_SCALE = 0.5`, **98.9 % of kernels sit AT the one-grid-step floor**
#: `knn_widths` applies (0.989 capture-OFF and capture-ON alike; still 0.92/0.76 at scale 2.0). The
#: population is DENSE on the log axis — ~30,700 regions, so the √n-th neighbour is nearer than one
#: grid cell for almost every region — and `knn_widths` is therefore doing exactly what it was designed
#: to do: reporting that the sample supports resolution finer than the axis can represent. ⛔ **The
#: consequence is that `_KNN_SCALE` is very nearly INERT at panel scale and the bandwidth ACTUALLY IN
#: FORCE is the grid step** — `(span in decades)/(_N_GRID − 1)`, ≈0.0336 dec over a data-derived
#: 8.7-decade span. So a bandwidth question asked of `_KNN_SCALE` cannot be answered on this substrate,
#: and this is the axis that can answer it.
_SWEEP_GRIDS = (65, 130, 260, 520, 1040)

_CLASS_NAMES = {
    int(RegionType.INTERGENIC): "intergenic",
    int(RegionType.INTRON): "intron",
    int(RegionType.EXON): "exon",
}

_EPS = 1e-12


# ---------------------------------------------------------------------------
# the scorers
# ---------------------------------------------------------------------------


def predictive_loglik(log_rho: np.ndarray, logP: np.ndarray, count, exposure) -> float:
    """Mean held-out ``log ∫ P(c | ρ·e)·P(log ρ) dρ`` per observation, over a fitted grid prior.

    ⭐ The prior is a discrete distribution on the grid, so the integral is a log-sum-exp over grid
    points of ``logP_k + log Poisson(c | exp(log_rho_k)·e)``. ⛔ The Poisson term is NOT row-normalised
    (``landscape._poisson_kernels`` normalises, which is right for a responsibility and wrong for a
    predictive — a normalised row cannot tell a well-predicted observation from a badly-predicted one).

    ⚠ An observation whose rate sits above the fitted grid's top scores LOW rather than being clamped
    to the edge value: the honest statement is that the prior did not cover it, and clamping would
    hide exactly the failure a held-out score exists to find.
    """
    c = np.asarray(count, dtype=np.float64)
    e = np.maximum(np.asarray(exposure, dtype=np.float64), _EPS)
    lp = np.asarray(logP, dtype=np.float64)
    lp = lp - np.log(np.sum(np.exp(lp - lp.max())) + _EPS) - lp.max()  # renormalise defensively
    grid = np.asarray(log_rho, dtype=np.float64)
    out = np.empty(c.shape[0], dtype=np.float64)
    chunk = 4096
    for lo in range(0, c.shape[0], chunk):
        hi = min(lo + chunk, c.shape[0])
        lam = np.exp(grid[None, :]) * e[lo:hi, None]
        ll = c[lo:hi, None] * np.log(np.maximum(lam, _EPS)) - lam - gammaln(c[lo:hi, None] + 1.0)
        a = ll + lp[None, :]
        mx = a.max(axis=1, keepdims=True)
        out[lo:hi] = (mx + np.log(np.sum(np.exp(a - mx), axis=1, keepdims=True))).ravel()
    return float(np.mean(out))


def census_of_curve(log_rho, logP, anchor_log_rho):
    """The SHIPPED census applied to any curve — `_census` for the basins, `split_basins` for the
    depleted/enriched selection. ⛔ Both are imported, never restated: an instrument that reimplements
    the rule it is scoring cannot detect a change in it."""
    ls = DensityLandscape(
        log_rho=np.asarray(log_rho, np.float64), logP=np.asarray(logP, np.float64), n_train=0
    )
    modes = _census(ls)
    depleted, enriched = split_basins(modes, float(anchor_log_rho))
    return modes, depleted, enriched


def modes_reproduce(modes_a, modes_b) -> float:
    """The share of half-A's modes that half-B also finds, each within its OWN fitted width.

    ⭐ The tolerance is the mode's own width — the fit's statement of how precisely it located itself
    — so no constant is chosen, and a mode that is genuinely sharp is held to a tight standard while a
    broad one is not. Mass-weighted by basin mass, so a wiggle owning no mass cannot dominate the
    verdict either way (the census's own threshold-free discipline)."""
    if not modes_a:
        return float("nan")
    wsum = sum(m.basin_mass for m in modes_a)
    if wsum <= 0.0:
        return float("nan")
    hit = 0.0
    for m in modes_a:
        tol = max(m.width, _EPS)
        if any(abs(n.log_rho - m.log_rho) <= tol for n in modes_b):
            hit += m.basin_mass
    return float(hit / wsum)


def restrict(mask, keep: np.ndarray):
    """A wall mask restricted to ``keep`` — how a HALF-SAMPLE fit is taken with no `src/` change.

    ⭐⭐ The module already has exactly one selection mechanism: a region that is not wall-exact on
    either side is not model-free, contributes zero count and zero exposure, and is excluded from the
    fit and from the anchor pool alike. Clearing both flags outside ``keep`` therefore expresses "fit
    on this half" in the estimator's OWN vocabulary — no new argument, and no second definition of
    the training population that could drift from the shipped one.
    """
    import dataclasses

    k = np.asarray(keep, bool)
    return dataclasses.replace(
        mask,
        start_exact=np.asarray(mask.start_exact, bool) & k,
        end_exact=np.asarray(mask.end_exact, bool) & k,
    )



# ---------------------------------------------------------------------------
# one condition
# ---------------------------------------------------------------------------


def floor_share(counts, exposure, sel, scale: float) -> float:
    """The share of kernels whose knn width is CLAMPED to one grid step — i.e. the share for which
    ``_KNN_SCALE`` has no effect at all.

    ⭐ This is the diagnostic that says whether a `_KNN_SCALE` sweep can mean anything on a given
    substrate: at 0.989 the knob is inert and the grid step is the bandwidth. Computed by asking
    ``knn_widths`` for the UNFLOORED widths (``grid_step = 0``) and comparing to the real step, so it
    reads the shipped kernel rather than a restatement of it."""
    from rigel.calibration.landscape import _grid, knn_widths

    s = np.asarray(sel, bool)
    c, e = np.asarray(counts, np.float64)[s], np.asarray(exposure, np.float64)[s]
    if c.size < 2:
        return float("nan")
    grid = _grid(c, e)
    step = float(grid[1] - grid[0])
    centres = np.clip(np.log10(np.maximum(c, 1.0)) - np.log10(e), grid[0], grid[-1])
    return float(np.mean(knn_widths(centres, 0.0, scale) <= step))


def grid_sweep(bundle, counts, exposure, keep_a, keep_b, ha, hb, grids) -> list[dict]:
    """Price the bandwidth ACTUALLY in force — the grid step — by the same held-out likelihood.

    ⛔⛔ **``_N_GRID`` IS PATCHED FOR THE ARM AND RESTORED, AND THE PATCH IS THE MEASUREMENT'S POINT.**
    It is a module constant documented as a computational budget ("discretization, not modelling");
    the floor share above says it is the binding bandwidth at panel scale, so it has to be swept to
    answer the question. The patch is a measurement arm, never a proposal, and it is undone in a
    ``finally`` so a raising fit cannot leak a changed constant into the rest of the run.
    """
    from rigel.calibration import landscape as _ls

    out = []
    original = _ls._N_GRID
    try:
        for g in grids:
            _ls._N_GRID = int(g)
            f_ab = fit_abundance_landscape(
                bundle["substrate"], bundle["region_arrays"], restrict(bundle["wall_mask"], keep_a)
            )
            f_ba = fit_abundance_landscape(
                bundle["substrate"], bundle["region_arrays"], restrict(bundle["wall_mask"], keep_b)
            )
            full = fit_abundance_landscape(
                bundle["substrate"], bundle["region_arrays"], bundle["wall_mask"]
            )
            if f_ab is None or f_ba is None or full is None:
                continue
            step = float(full.landscape.log_rho[1] - full.landscape.log_rho[0]) / np.log(10.0)
            out.append(
                {
                    "n_grid": int(g),
                    "grid_step_dec": step,
                    "ll_AB": predictive_loglik(
                        f_ab.landscape.log_rho, f_ab.landscape.logP, counts[hb], exposure[hb]
                    ),
                    "ll_BA": predictive_loglik(
                        f_ba.landscape.log_rho, f_ba.landscape.logP, counts[ha], exposure[ha]
                    ),
                    "n_modes_full": len(full.modes),
                    "mode_reproducibility": modes_reproduce(f_ab.modes, f_ba.modes),
                    "rho_0": full.rho_0,
                    "span_R": full.span_R,
                    "anchor_gap_nats": full.anchor_gap_nats,
                    "anchor_consistent": bool(full.anchor_consistent),
                }
            )
    finally:
        _ls._N_GRID = original
    for r in out:
        r["ll_mean"] = 0.5 * (r["ll_AB"] + r["ll_BA"])
    return out


def measure(bundle: dict, condition: str, *, scales=(), grids=()) -> dict:
    """One condition's head-to-head. PURE (arrays in, dict out) so the self-test feeds it synthetics.

    ``bundle`` carries: ``counts``/``exposure``/``model_free`` (the measured-total triple),
    ``mass``/``eff_gdna`` on the REGION axis and on the FULL chain axis, ``coarse``, and
    ``truth_gdna_counts`` (the gdna partition's counts through the same side selection).
    """
    counts = np.asarray(bundle["counts"], np.float64)
    exposure = np.asarray(bundle["exposure"], np.float64)
    sel = np.asarray(bundle["model_free"], bool) & (exposure > 0.0)
    coarse = np.asarray(bundle["coarse"], np.int64)
    ss, cap = PVO.stratum(condition)

    gates: list[dict] = []
    if "partition_sum" in bundle:
        s_parts, s_full = bundle["partition_sum"]
        gates.append(
            {
                "gate": "gate partitions-sum-to-full[start]",
                "ok": bool(np.array_equal(s_parts, s_full)),
                "detail": f"max|delta| = {float(np.max(np.abs(s_parts - s_full))):.0f}",
            }
        )

    # ── the two fits, on everything selected (the reported curves)
    al = fit_abundance_landscape(bundle["substrate"], bundle["region_arrays"], bundle["wall_mask"])
    if al is None:
        return {"condition": condition, "skipped": "fewer than two model-free regions"}

    # ── the anchor pool, and the CERTIFIED truth for the depleted level
    anchor_pool = sel & (coarse == int(RegionType.INTERGENIC))
    anchor_log_rho = al.anchor_log_rho
    if anchor_pool.any() and float(exposure[anchor_pool].sum()) > 0.0:
        gates.append(
            {
                "gate": "gate anchor-pool-matches-the-shipped-census",
                "ok": bool(
                    counts[anchor_pool].sum() <= 0.0
                    or abs(
                        float(
                            np.log(counts[anchor_pool].sum()) - np.log(exposure[anchor_pool].sum())
                        )
                        - float(anchor_log_rho)
                    )
                    < 1e-9
                ),
                "detail": "the instrument's intergenic pool IS the census's",
            }
        )
    tg = np.asarray(bundle["truth_gdna_counts"], np.float64)
    t_num = float(tg[anchor_pool].sum())
    t_den = float(exposure[anchor_pool].sum())
    truth_log_rho = float(np.log(t_num) - np.log(t_den)) if (t_num > 0 and t_den > 0) else float("nan")
    rna_share = (
        float((counts[anchor_pool].sum() - t_num) / counts[anchor_pool].sum())
        if counts[anchor_pool].sum() > 0
        else float("nan")
    )

    a_modes, a_dep, a_enr = al.modes, al.depleted, al.enriched

    # ── ⓑ held-out predictive likelihood, both arms on the SAME held-out pairs
    idx = np.flatnonzero(sel)
    keep_a = np.zeros(counts.shape[0], bool)
    keep_b = np.zeros(counts.shape[0], bool)
    keep_a[idx[0::2]] = True
    keep_b[idx[1::2]] = True
    ha, hb = idx[0::2], idx[1::2]
    ll = {}
    fits = {}
    if ha.size >= 2 and hb.size >= 2:
        for keep, sc_h, tag in ((keep_a, hb, "AB"), (keep_b, ha, "BA")):
            a_fit = fit_abundance_landscape(
                bundle["substrate"], bundle["region_arrays"], restrict(bundle["wall_mask"], keep)
            )
            fits[tag] = a_fit
            ll[f"landscape_{tag}"] = (
                predictive_loglik(
                    a_fit.landscape.log_rho, a_fit.landscape.logP, counts[sc_h], exposure[sc_h]
                )
                if a_fit is not None
                else float("nan")
            )
        if fits.get("AB") is not None and fits.get("BA") is not None:
            ll["landscape_mode_reproducibility"] = modes_reproduce(
                fits["AB"].modes, fits["BA"].modes
            )

    # ── the _KNN_SCALE sweep, scored by the SAME held-out likelihood
    sweep = []
    for s in scales:
        if ha.size < 2 or hb.size < 2:
            break
        f_ab = fit_abundance_landscape(
            bundle["substrate"],
            bundle["region_arrays"],
            restrict(bundle["wall_mask"], keep_a),
            knn_scale=s,
        )
        f_ba = fit_abundance_landscape(
            bundle["substrate"],
            bundle["region_arrays"],
            restrict(bundle["wall_mask"], keep_b),
            knn_scale=s,
        )
        if f_ab is None or f_ba is None:
            continue
        full = fit_abundance_landscape(
            bundle["substrate"], bundle["region_arrays"], bundle["wall_mask"], knn_scale=s
        )
        sweep.append(
            {
                "knn_scale": float(s),
                "ll_AB": predictive_loglik(
                    f_ab.landscape.log_rho, f_ab.landscape.logP, counts[hb], exposure[hb]
                ),
                "ll_BA": predictive_loglik(
                    f_ba.landscape.log_rho, f_ba.landscape.logP, counts[ha], exposure[ha]
                ),
                "n_modes_full": len(full.modes),
                "mode_reproducibility": modes_reproduce(f_ab.modes, f_ba.modes),
                "rho_0": full.rho_0,
                "span_R": full.span_R,
                "anchor_gap_nats": full.anchor_gap_nats,
            }
        )
    for row in sweep:
        row["ll_mean"] = 0.5 * (row["ll_AB"] + row["ll_BA"])
    for row in sweep:
        row["floor_share"] = floor_share(counts, exposure, sel, row["knn_scale"])

    gsweep = (
        grid_sweep(bundle, counts, exposure, keep_a, keep_b, ha, hb, grids)
        if grids and ha.size >= 2 and hb.size >= 2
        else []
    )

    return {
        "condition": condition,
        "strand": ss,
        "capture": cap,
        "scope": OC._scope(condition),
        "zero_gdna": bool(PVO.is_zero_gdna(condition)),
        "gates": gates,
        "n_selected": int(sel.sum()),
        "landscape": {
            "n_modes": len(a_modes),
            "rho_0": al.rho_0,
            "span_R": al.span_R,
            "depleted_log_rho": a_dep.log_rho,
            "depleted_width": a_dep.width,
            "enriched_log_rho": None if a_enr is None else a_enr.log_rho,
            "enriched_mass": None if a_enr is None else a_enr.basin_mass,
            "anchor_gap_nats": al.anchor_gap_nats,
            "curve": {
                "log_rho": al.landscape.log_rho.tolist(),
                "logP": al.landscape.logP.tolist(),
            },
        },
        "truth": {
            "anchor_log_rho_total": float(anchor_log_rho),
            "anchor_log_rho_gdna": truth_log_rho,
            "rna_share_of_anchor_pool": rna_share,
            "vacuous": not np.isfinite(truth_log_rho),
            "landscape_err_nats": float(abs(a_dep.log_rho - truth_log_rho))
            if np.isfinite(truth_log_rho)
            else float("nan"),
        },
        "loglik": ll,
        "sweep": sweep,
        "grid_sweep": gsweep,
        "floor_share_shipped": floor_share(counts, exposure, sel, 0.5),
    }



# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------


def load_condition(suite: Path, condition: str, index, cached: dict) -> dict:
    cache = read_scan_cache(Path(suite) / "scan_cache" / condition, index)
    payload = cache.payload
    ra = cached["region_arrays"]
    substrate = CalibrationSubstrate.from_payload(payload, ra)
    # ⭐ NO geometry is built here any more. It existed only to supply the NPMLE arm's
    # `(mass, eff_gdna)` pair; with that arm gone (see the header) this loader needs counts, lengths
    # and the wall mask alone — which is the whole point of the estimator that replaced it.
    mask = build_region_wall_mask(
        ra,
        cached["mature"],
        cached["reach"][0],
        cached["reach"][1],
        w_max=w_max_from_deposited_lengths(payload.deposited_lengths),
    )
    counts, exposure, model_free = region_counts_and_exposure(substrate, ra, mask)

    # the CERTIFIED truth: the gdna partition's counts through the SAME side selection, and the
    # sum-to-full gate on the START banks
    root = Path(suite) / "oracle_cache" / condition
    parts = {}
    for k in ORIGINS:
        d = root / k
        if not d.is_dir():
            raise FileNotFoundError(
                f"{d} is missing — run `panel.py cache` first; this instrument scores the depleted "
                "level against the origin partitions only."
            )
        parts[k] = read_scan_cache(d, index).payload
    s_parts = sum(
        np.asarray(p.region_start_count, np.float64).sum(1) for p in parts.values()
    )
    s_full = np.asarray(payload.region_start_count, np.float64).sum(1)
    truth_counts, _te, _tm = region_counts_and_exposure(
        CalibrationSubstrate.from_payload(parts["gdna"], ra), ra, mask
    )

    return {
        "counts": counts,
        "exposure": exposure,
        "model_free": model_free,
        "coarse": coarse_type_array(np.asarray(ra.signature, np.int64)),
        "truth_gdna_counts": truth_counts,
        "partition_sum": (s_parts, s_full),
        "substrate": substrate,
        "region_arrays": ra,
        "wall_mask": mask,
    }


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------


def _f(x, w=9, p=4):
    return f"{'—':>{w}}" if x is None or not np.isfinite(x) else f"{x:>{w}.{p}f}"


def report(d: dict) -> int:
    if "skipped" in d:
        print(f"\n── {d['condition']}: SKIPPED ({d['skipped']})")
        return 0
    bad = 0
    print(f"\n── {d['condition']}  [{d['strand']} x {d['capture']}]  {d['scope']}")
    for g in d["gates"]:
        flag = "✔" if g["ok"] else "⛔"
        if not g["ok"]:
            bad += 1
        print(f"   {flag} {g['gate']}: {g['detail']}")
    a = d["landscape"]
    print(
        f"   CENSUS  landscape: {a['n_modes']} modes  rho_0={a['rho_0']:.4g}  R={a['span_R']:.1f}  "
        f"anchor gap {_f(a['anchor_gap_nats'], 6, 3)} nats"
    )
    t = d["truth"]
    if t["vacuous"]:
        print(
            "   ⓒ DEPLETED vs CERTIFIED gDNA TRUTH: ⛔ STAMPED VACUOUS — the anchor pool's true "
            "gDNA rate is exactly 0 (a zero-gDNA condition), and a nats error against zero is not a "
            "measurement"
        )
    else:
        print(
            f"   ⓒ DEPLETED vs CERTIFIED gDNA TRUTH (true log rho = {t['anchor_log_rho_gdna']:+.4f}, "
            f"RNA share of the anchor pool {t['rna_share_of_anchor_pool']:.2e}):"
        )
        print(f"        landscape {_f(t['landscape_err_nats'], 8, 4)} nats")
    ll = d["loglik"]
    if ll:
        lsc = 0.5 * (ll.get("landscape_AB", np.nan) + ll.get("landscape_BA", np.nan))
        print(
            f"   ⓐ HELD-OUT PREDICTIVE LOGLIK per region (higher is better): {_f(lsc, 10, 3)}"
        )
        if "landscape_mode_reproducibility" in ll:
            print(
                f"        split-half mode reproducibility (mass-weighted, own-width tolerance): "
                f"{_f(ll['landscape_mode_reproducibility'], 6, 3)}"
            )
    print(
        f"   ⚠ KERNELS AT THE ONE-GRID-STEP FLOOR at the shipped _KNN_SCALE = 0.5: "
        f"{d['floor_share_shipped']:.3f} — above ~0.9 the knob is INERT and the GRID STEP is the "
        f"bandwidth in force"
    )
    if d["sweep"]:
        print("   ② _KNN_SCALE SWEEP — held-out loglik selects, mode count and stability describe")
        best = max(d["sweep"], key=lambda r: r["ll_mean"])
        for r in d["sweep"]:
            star = " ⭐ best" if r is best else ""
            ship = " (SHIPPED)" if r["knn_scale"] == 0.5 else ""
            print(
                f"        scale {r['knn_scale']:>6.3f}{ship:<10} ll {_f(r['ll_mean'], 10, 3)}  "
                f"modes {r['n_modes_full']:>3}  reproducibility {_f(r['mode_reproducibility'], 6, 3)}  "
                f"floor {_f(r.get('floor_share'), 6, 3)}  R {r['span_R']:>7.1f}{star}"
            )
    if d["grid_sweep"]:
        print("   ② _N_GRID SWEEP — the bandwidth ACTUALLY in force (the grid step)")
        best = max(d["grid_sweep"], key=lambda r: r["ll_mean"])
        for r in d["grid_sweep"]:
            star = " ⭐ best" if r is best else ""
            ship = " (SHIPPED)" if r["n_grid"] == 260 else ""
            cons = "✔" if r["anchor_consistent"] else "⛔"
            print(
                f"        n_grid {r['n_grid']:>5}{ship:<10} step {r['grid_step_dec']:.4f} dec  "
                f"ll {_f(r['ll_mean'], 10, 3)}  modes {r['n_modes_full']:>3}  "
                f"reproducibility {_f(r['mode_reproducibility'], 6, 3)}  {cons} anchor gap "
                f"{_f(r['anchor_gap_nats'], 6, 3)}  R {r['span_R']:>7.1f}{star}"
            )
    return bad


def summarise(rows: list[dict]) -> None:
    rows = [r for r in rows if "skipped" not in r]
    if not rows:
        return
    print(f"\n{'=' * 118}")
    print("PER STRATUM — never pooled across strata; the zero controls are their own row")
    print(f"{'=' * 118}")
    hdr = (
        f"{'selection':<44}{'n':>3}  {'held-out ll':>11}  {'nats vs truth':>13}"
    )
    print(hdr)
    print("-" * len(hdr))
    for name, pred in OC._SELECTIONS:
        sub = [r for r in rows if pred(r["condition"])]
        if not sub:
            continue

        def _m(f):
            v = [f(r) for r in sub]
            v = [x for x in v if x is not None and np.isfinite(x)]
            return float(np.mean(v)) if v else float("nan")

        la = _m(lambda r: 0.5 * (r["loglik"].get("landscape_AB", np.nan) + r["loglik"].get("landscape_BA", np.nan)))
        na = _m(lambda r: r["truth"]["landscape_err_nats"])
        print(f"{name:<44}{len(sub):>3}  {_f(la, 11, 3)}  {_f(na, 13, 4)}")
    sweeps = [r for r in rows if r["sweep"]]
    if sweeps:
        print(f"\n{'=' * 118}")
        print("② THE _KNN_SCALE SWEEP, PER STRATUM — the held-out likelihood is the selector")
        print(f"{'=' * 118}")
        for name, pred in OC._SELECTIONS:
            sub = [r for r in sweeps if pred(r["condition"])]
            if not sub:
                continue
            print(f"\n{name}  (n={len(sub)})")
            for s in _SWEEP_SCALES:
                vals = [
                    next((x for x in r["sweep"] if x["knn_scale"] == s), None) for r in sub
                ]
                vals = [v for v in vals if v is not None]
                if not vals:
                    continue
                lm = float(np.mean([v["ll_mean"] for v in vals]))
                nm = float(np.mean([v["n_modes_full"] for v in vals]))
                rp = float(np.nanmean([v["mode_reproducibility"] for v in vals]))
                wins = sum(
                    1
                    for r in sub
                    if max(r["sweep"], key=lambda x: x["ll_mean"])["knn_scale"] == s
                )
                ship = " (SHIPPED)" if s == 0.5 else ""
                print(
                    f"   scale {s:>6.3f}{ship:<10} mean ll {lm:>10.3f}   mean modes {nm:>5.1f}   "
                    f"reproducibility {rp:>5.3f}   best on {wins}/{len(sub)} conditions"
                )
    gs = [r for r in rows if r.get("grid_sweep")]
    if gs:
        print(f"\n{'=' * 118}")
        print("② THE _N_GRID SWEEP, PER STRATUM — the bandwidth actually in force")
        print(f"{'=' * 118}")
        for name, pred in OC._SELECTIONS:
            sub = [r for r in gs if pred(r["condition"])]
            if not sub:
                continue
            print(f"\n{name}  (n={len(sub)})")
            for g in _SWEEP_GRIDS:
                vals = [
                    next((x for x in r["grid_sweep"] if x["n_grid"] == g), None) for r in sub
                ]
                vals = [v for v in vals if v is not None]
                if not vals:
                    continue
                wins = sum(
                    1
                    for r in sub
                    if max(r["grid_sweep"], key=lambda x: x["ll_mean"])["n_grid"] == g
                )
                ship = " (SHIPPED)" if g == 260 else ""
                print(
                    f"   n_grid {g:>5}{ship:<10} step {np.mean([v['grid_step_dec'] for v in vals]):.4f} dec  "
                    f"mean ll {np.mean([v['ll_mean'] for v in vals]):>10.3f}   "
                    f"mean modes {np.mean([v['n_modes_full'] for v in vals]):>5.1f}   "
                    f"reproducibility {np.nanmean([v['mode_reproducibility'] for v in vals]):>5.3f}   "
                    f"anchor-consistent {sum(v['anchor_consistent'] for v in vals)}/{len(vals)}   "
                    f"best on {wins}/{len(sub)}"
                )
    print(
        "\n⛔ A `_KNN_SCALE` sweep cannot answer the bandwidth question on this substrate: ~99 % of "
        "kernels are clamped to one grid step, so the knob is inert and the GRID STEP is the "
        "bandwidth — which is what the `_N_GRID` sweep prices."
    )


# ---------------------------------------------------------------------------
# self-test — perturbed, no I/O
# ---------------------------------------------------------------------------


def self_test() -> int:  # noqa: C901
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests" / "calibration"))
    from test_abundance_landscape import bimodal_parts, parts, unimodal_parts

    ok = fail = 0

    def check(name, cond):
        nonlocal ok, fail
        if cond:
            ok += 1
        else:
            fail += 1
            print(f"⛔ {name}")

    # ── the predictive likelihood: a prior AT the truth must beat one three decades off
    rng = np.random.default_rng(3)
    e = np.full(400, 1000.0)
    c = rng.poisson(0.05 * 1000.0, 400).astype(np.float64)
    grid = np.linspace(np.log(1e-6), np.log(10.0), 240)

    def _spike(at):
        lp = -0.5 * ((grid - np.log(at)) / 0.05) ** 2
        return lp - np.log(np.sum(np.exp(lp)))

    right = predictive_loglik(grid, _spike(0.05), c, e)
    wrong = predictive_loglik(grid, _spike(50.0), c, e)
    check("a prior at the truth beats one three decades off", right > wrong)
    flat = predictive_loglik(grid, np.full_like(grid, -np.log(grid.size)), c, e)
    check("...and the truth also beats a flat prior", right > flat)
    check("a wrong spike is worse than being uninformative", flat > wrong)
    check(
        "the score is a MEAN per observation (doubling the data does not move it)",
        abs(
            predictive_loglik(grid, _spike(0.05), np.concatenate([c, c]), np.concatenate([e, e]))
            - right
        )
        < 1e-9,
    )
    check("it is deterministic", predictive_loglik(grid, _spike(0.05), c, e) == right)
    # ⛔ the un-normalised Poisson is the point: a row-normalised kernel cannot rank two priors
    check(
        "an observation ABOVE the grid top is not clamped to the edge",
        predictive_loglik(grid, _spike(0.05), np.full(5, 1e6), np.full(5, 1000.0))
        < predictive_loglik(grid, _spike(0.05), np.full(5, 50.0), np.full(5, 1000.0)),
    )

    # ── the census reuse: the shipped rule applied to a curve reproduces the shipped answer
    counts, lengths, sig, rho_lo, rho_hi = bimodal_parts()
    sub, ra, mask = parts(counts, lengths, sig)
    al = fit_abundance_landscape(sub, ra, mask)
    m2, d2, e2 = census_of_curve(al.landscape.log_rho, al.landscape.logP, al.anchor_log_rho)
    check("census_of_curve reproduces the shipped census", len(m2) == len(al.modes))
    check("...and the shipped depleted mode", abs(d2.log_rho - al.depleted.log_rho) < 1e-12)
    check(
        "...and the shipped enriched mode",
        (e2 is None) == (al.enriched is None)
        and (e2 is None or abs(e2.log_rho - al.enriched.log_rho) < 1e-12),
    )

    # ── mode reproducibility: identical fits reproduce, a shifted one does not
    check("identical censuses reproduce fully", modes_reproduce(al.modes, al.modes) == 1.0)
    from rigel.calibration.abundance_landscape import AbundanceMode

    shifted = tuple(
        AbundanceMode(m.log_rho + 50.0, m.basin_mass, m.width, m.lo, m.hi) for m in al.modes
    )
    check("a 50-nat shift reproduces nothing", modes_reproduce(al.modes, shifted) == 0.0)
    check(
        "the tolerance is the mode's OWN width",
        modes_reproduce(
            (AbundanceMode(0.0, 1.0, 0.5, -1.0, 1.0),), (AbundanceMode(0.4, 1.0, 0.5, -1.0, 1.0),)
        )
        == 1.0
        and modes_reproduce(
            (AbundanceMode(0.0, 1.0, 0.5, -1.0, 1.0),), (AbundanceMode(0.6, 1.0, 0.5, -1.0, 1.0),)
        )
        == 0.0,
    )
    check("an empty census is not a pass", np.isnan(modes_reproduce((), al.modes)))

    # ── ⛔ `restrict` MUST actually restrict: an inert one silently turns every held-out score into an
    # in-sample one, which is the failure `TRAPS: an-ablation-that-never-ran` names. The two halves are
    # disjoint, so their fits cannot be identical.
    idx_st = np.flatnonzero(region_counts_and_exposure(sub, ra, mask)[2])
    ka = np.zeros(mask.n_regions, bool)
    kb = np.zeros(mask.n_regions, bool)
    ka[idx_st[0::2]] = True
    kb[idx_st[1::2]] = True
    fa = fit_abundance_landscape(sub, ra, restrict(mask, ka))
    fb = fit_abundance_landscape(sub, ra, restrict(mask, kb))
    check(
        "restrict is NOT inert — disjoint halves give different fits",
        not np.array_equal(fa.landscape.logP, fb.landscape.logP),
    )
    check(
        "...and each half trains on about half the population",
        abs(fa.n_train - fb.n_train) <= 1 and fa.n_train < al.n_train,
    )
    check(
        "restrict never ADDS a region the shipped mask excluded",
        int(np.count_nonzero(region_counts_and_exposure(sub, ra, restrict(mask, ka))[2]))
        <= int(np.count_nonzero(region_counts_and_exposure(sub, ra, mask)[2])),
    )

    # ── the whole condition path on synthetics, both shapes
    def _bundle(counts_a, lengths_a, sig_a, *, truth_frac=1.0):
        s, r, m = parts(counts_a, lengths_a, sig_a)
        cts, exp, mf = region_counts_and_exposure(s, r, m)
        return {
            "counts": cts,
            "exposure": exp,
            "model_free": mf,
            "coarse": coarse_type_array(np.asarray(r.signature, np.int64)),
            "truth_gdna_counts": cts * truth_frac,
            "substrate": s,
            "region_arrays": r,
            "wall_mask": m,
        }

    b = _bundle(counts, lengths, sig)
    d = measure(b, "gdna_g50_ss_0.99_nrna_mid_capture_on")
    check("a bimodal condition reads bimodal", d["landscape"]["span_R"] > 10.0)
    check("the stratum parses", (d["strand"], d["capture"]) == ("stranded", "capture ON"))
    check("the held-out likelihood is finite", np.isfinite(d["loglik"]["landscape_AB"]))
    check("the depleted level is scored against gDNA truth", np.isfinite(d["truth"]["landscape_err_nats"]))
    check("the truth arm sees a pure-gDNA anchor pool", d["truth"]["rna_share_of_anchor_pool"] == 0.0)
    check("json round-trips", json.dumps(d) is not None)

    # ⛔ a zero-gDNA anchor pool must be STAMPED, never scored
    b0 = _bundle(counts, lengths, sig, truth_frac=0.0)
    d0 = measure(b0, "gdna_g00_ss_0.99_nrna_mid_capture_off")
    check("a zero-gDNA truth is stamped vacuous", d0["truth"]["vacuous"])
    check("...and reports no nats error", not np.isfinite(d0["truth"]["landscape_err_nats"]))

    c2, l2, s2 = unimodal_parts()
    d2u = measure(_bundle(c2, l2, s2), "gdna_g50_ss_0.50_nrna_mid_capture_off")
    check("a unimodal condition has span exactly 1", d2u["landscape"]["span_R"] == 1.0)

    # ── the sweep: it must MOVE, and the shipped scale must be one of the points scored
    ds = measure(b, "gdna_g50_ss_0.99_nrna_mid_capture_on", scales=_SWEEP_SCALES)
    check("the sweep scores every scale", len(ds["sweep"]) == len(_SWEEP_SCALES))
    check("the shipped 0.5 is scored", any(r["knn_scale"] == 0.5 for r in ds["sweep"]))
    check(
        "the sweep is not inert — the likelihood MOVES with the scale",
        len({round(r["ll_mean"], 6) for r in ds["sweep"]}) > 1,
    )
    check(
        "a wider kernel finds no more modes than a narrower one",
        min(r["n_modes_full"] for r in ds["sweep"] if r["knn_scale"] >= 2.0)
        <= max(r["n_modes_full"] for r in ds["sweep"] if r["knn_scale"] <= 0.25),
    )
    check(
        "the shipped scale's census matches the unswept fit",
        next(r for r in ds["sweep"] if r["knn_scale"] == 0.5)["n_modes_full"]
        == d["landscape"]["n_modes"],
    )

    # ── the FLOOR SHARE: it must be 0 on a deliberately sparse population and 1 on a dense one, or it
    # cannot support the claim that `_KNN_SCALE` is inert at panel scale
    sparse_c = np.array([1.0, 10.0, 1e3, 1e5], np.float64)
    sparse_e = np.full(4, 1e4)
    fs_sparse = floor_share(sparse_c, sparse_e, np.ones(4, bool), 0.5)
    check("a SPARSE population is not at the floor", fs_sparse < 0.5)
    # ⛔⛔ The dense fixture must be dense but NOT DEGENERATE. A first draft used 4,000 IDENTICAL
    # counts, so every kernel centre coincided, ``d_k`` was exactly 0, and the gate passed under a
    # perturbation that removed the floor comparison entirely — `TRAPS.md`'s "a zero arm can be
    # degenerate and then it tests nothing", found by perturbation 2026-08-21. Jittered counts give
    # distinct-but-closely-spaced centres, which is the real panel's shape.
    dense = 500.0 + (np.arange(4000) % 97).astype(np.float64)
    fs_dense = floor_share(dense, np.full(4000, 1e4), np.ones(4000, bool), 0.5)
    check("a DENSE population is entirely at the floor", fs_dense > 0.99)
    check(
        "...and its centres are DISTINCT, so the floor is doing the work rather than a tie",
        np.unique(np.log10(dense) - np.log10(1e4)).size > 50,
    )
    check(
        "the floor share RISES as the scale falls (a smaller scale clamps more)",
        floor_share(sparse_c, sparse_e, np.ones(4, bool), 0.01) >= fs_sparse,
    )
    check("fewer than two regions is not a number", np.isnan(floor_share(
        np.array([1.0]), np.array([1.0]), np.ones(1, bool), 0.5
    )))

    # ── the _N_GRID sweep: it must MOVE, the step must FALL as the grid grows, and the constant must
    # be RESTORED even when a fit raises
    from rigel.calibration import landscape as _ls

    before = _ls._N_GRID
    dg = measure(b, "gdna_g50_ss_0.99_nrna_mid_capture_on", grids=_SWEEP_GRIDS)
    check("_N_GRID is restored after the sweep", _ls._N_GRID == before)
    check("the grid sweep scores every point", len(dg["grid_sweep"]) == len(_SWEEP_GRIDS))
    steps = [r["grid_step_dec"] for r in dg["grid_sweep"]]
    check("a finer grid has a smaller step", steps == sorted(steps, reverse=True))
    check(
        "the sweep is not inert — the likelihood MOVES with the grid",
        len({round(r["ll_mean"], 6) for r in dg["grid_sweep"]}) > 1,
    )
    check(
        "a finer grid finds at least as many modes as a coarser one",
        max(r["n_modes_full"] for r in dg["grid_sweep"] if r["n_grid"] >= 520)
        >= min(r["n_modes_full"] for r in dg["grid_sweep"] if r["n_grid"] <= 130),
    )
    check(
        "the shipped 260 reproduces the unswept census exactly",
        next(r for r in dg["grid_sweep"] if r["n_grid"] == 260)["n_modes_full"]
        == d["landscape"]["n_modes"],
    )
    try:
        grid_sweep(
            {"substrate": None, "region_arrays": None, "wall_mask": None},
            np.zeros(2), np.zeros(2), np.ones(2, bool), np.ones(2, bool),
            np.array([0]), np.array([1]), (64,),
        )
    except Exception:
        pass
    check("_N_GRID is restored even when the arm RAISES", _ls._N_GRID == before)

    print(f"self-test: {ok} passed, {fail} failed")
    return 1 if fail else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--condition", default=None)
    ap.add_argument("--sweep", action="store_true", help="price _KNN_SCALE by held-out likelihood")
    ap.add_argument(
        "--grid-sweep",
        action="store_true",
        help="price _N_GRID — the bandwidth actually in force once the knn floor binds",
    )
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)
    if args.self_test:
        return self_test()

    index = TranscriptIndex.load(args.index)
    ra = RegionArrays.from_index(index)
    cached = {
        "region_arrays": ra,
        "sj": build_sj_geometry_arrays(index),
        "mature": build_mature_wall_distances(index, ra),
        "reach": build_contiguous_boundary_reach_arrays(index),
    }
    conds = (
        [args.condition]
        if args.condition
        else sorted(p.name for p in (args.suite / "scan_cache").iterdir() if p.is_dir())
    )
    scales = _SWEEP_SCALES if args.sweep else ()
    grids = _SWEEP_GRIDS if args.grid_sweep else ()
    rows, bad = [], 0
    for c in conds:
        bundle = load_condition(args.suite, c, index, cached)
        d = measure(bundle, c, scales=scales, grids=grids)
        bad += report(d)
        rows.append(d)
    summarise(rows)
    if args.json:
        args.json.write_text(json.dumps({"suite": str(args.suite), "rows": rows}))
        print(f"wrote {args.json}")
    if bad:
        print(f"\n⛔ {bad} gate(s) FAILED — the substrate is not validated, so no score above stands.")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
