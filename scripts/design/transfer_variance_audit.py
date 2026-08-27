"""WHERE DOES THE MESSAGE TRANSFER VARIANCE COME FROM, AND COULD THE AbundanceLandscape SUPPLY A
BETTER ONE? ⭐⭐⭐ **AN AUDIT: no solver runs, no EM, nothing in `src/` is patched, and it proposes no
change** — it establishes what is live, what is dead, and what a TOTAL's variance can and cannot say.

⭐⭐ **THE HISTORY, VERIFIED RATHER THAN QUOTED.** ``DensityNPMLE.project`` once supplied
``sigma^2_transfer`` as ``var_proj[dst] + (mu_proj[dst] - mu_proj[src])^2`` — a population quantity read
off a fitted landscape: a between-mode ambiguity plus a within-mode floor ``h^2``. ⭐⭐ **The module is
now DELETED outright (2026-08-21), so gate ① checks the strongest available fact — the file is gone and
nothing under ``src/`` calls ``.project()``** — rather than trusting a comment that says it is unused.
⚠ The scanner's own falsification moved with the deletion: it used to prove it could still SEE a call by
scanning ``test_npmle.py``, which no longer exists, so it now proves it on a synthetic source and also
proves a PROSE mention does not count.

⭐⭐⭐ **WHAT REPLACED IT, AND THE ONE STRUCTURAL DIFFERENCE THAT MATTERS.** The shipped law is
``messages.variance.transfer_logvar`` = ``Var(log rho_tot^dst) + Var(log rho_tot^src)``, each from
``composition_logvar`` = ``trigamma(M + 1/2)`` (counting) + a composition term that is exactly 0 wherever
``Var(f_g) = 0``; and it is identically 0 on a SPLICE IN, where the reframe is common-mode. ⛔ **Every
term in it SHRINKS AS EITHER SLOT IS COUNTED MORE DEEPLY.** The retired proxy's between-mode term did
not. So the retirement removed a POPULATION-HETEROGENEITY term and replaced it with sampling terms
only, and arm ⓑ measures what that costs: on the deeply-counted hops a transport is nearly FREE.

⭐ The Stage-3 ``CurrencyPolicy`` fills that hole by a different and better route —
``variance.premise_logvar``, a method-of-moments fit of the leftover spread of ``log r`` after each
hop's own counting variance is removed, floored at 0. ⭐⭐ **So the role the npmle held is not vacant
in the tool; it is vacant in the SHIPPED policy** (``message_policy`` defaults to ``relay``), and that
is this audit's central finding.

**ⓐ THE PER-SLOT POSTERIOR VARIANCE THE LANDSCAPE OFFERS, AGAINST THE COUNTING ONE IT WOULD REPLACE.**
The landscape gives ``Var(log rho | count, exposure)`` on its own grid in one line — the machinery the
npmle's ``project`` had, fitted on wall-exact totals instead of the forbidden ``mass/eff_gdna`` pair. The
audit compares it, per region, to ``count_logvar(count)``. ⛔⛔ **A RATIO BELOW 1 IS THE WRONG
DIRECTION AND IS THE POINT OF THIS ARM**: a fitted population prior SHRINKS a slot's rate toward a
mode, so its posterior is TIGHTER than the slot's own count warrants. Substituting it would make every
message MORE confident, which is the "transports for free" failure stated above, arriving through a new
door.

**ⓑ WHAT THE SHIPPED LAW ACTUALLY IS AT PASS-0 INIT, MEASURED ON REAL HOPS.** At init every slot
carries the structural default (``Var(f_g) = 0``), so ``composition_logvar`` reduces EXACTLY to
``count_logvar(M)`` — the shipped ``sigma^2_transfer`` is then the two counting terms and nothing else.
Reported as a distribution over the chain's real adjacent pairs, with the SPLICE-IN-style zero rows
counted separately, and beside it the deep-hop share: the fraction of hops where both slots are counted
deeply enough that the whole damping term is below one part in ten.

**ⓒ THE HETEROGENEITY THE LANDSCAPE SEES — the quantity a premise would be fitted from.** The
between-basin spread of the fitted density (its mass-weighted variance of ``log rho``, and the same
variance WITHIN the depleted basin alone) is a population number that does not shrink with any slot's
count. Reported beside ``premise_logvar`` computed from the same condition's adjacent-pair log ratios,
so the two candidate estimators of one quantity can be read against each other.

⛔⛔ **AND THE LIMIT IS STATED WITH THE NUMBERS, NEVER INFERRED FROM THEM: A TOTAL IS
COMPOSITION-VACUOUS.** ``rho_tot`` is what a message must be damped against only when the transported
currency is a LEVEL; a COMPOSITION hop needs the uncertainty of a gDNA-vs-RNA split, and no function of
the total can supply it (that vacuity is why the npmle was never allowed near the composition arm). So
the honest scope of any landscape-supplied variance is the LEVEL hops and the population premise — not
the composition ones. Arm ⓒ is reported as a candidate for that scope alone.

Usage::

    python scripts/design/transfer_variance_audit.py                     # every cached condition
    python scripts/design/transfer_variance_audit.py --condition NAME    # one condition
    python scripts/design/transfer_variance_audit.py --self-test         # perturbed, no I/O
"""

from __future__ import annotations

import argparse
import ast
import importlib.util
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, Path(__file__).resolve().parent / name)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


HTH = _sibling("landscape_head_to_head.py")
OC = _sibling("object_composition.py")
PVO = OC.PVO

from rigel.calibration.abundance_landscape import fit_abundance_landscape  # noqa: E402
from rigel.calibration.landscape import _poisson_kernels  # noqa: E402
from rigel.calibration.messages.variance import premise_logvar  # noqa: E402
from rigel.calibration.messages.variance import count_logvar, transfer_logvar  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import REGION, build_region_chain  # noqa: E402
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

#: "the damping is negligible" — a hop whose whole transfer term is below this share of one nat of
#: log-density leaves the message essentially undamped. ⛔ NOT a model constant and nothing branches on
#: it in `src/`: it is a REPORTING cut for arm ⓑ's descriptive share, and the distribution is printed
#: beside it so the reader never depends on the cut.
_NEGLIGIBLE_LOGVAR = 0.1

_CLASS_NAMES = {
    int(RegionType.INTERGENIC): "intergenic",
    int(RegionType.INTRON): "intron",
    int(RegionType.EXON): "exon",
}

_EPS = 1e-12


# ---------------------------------------------------------------------------
# ① the retirement gate — AST, not a grep, and it names what it found
# ---------------------------------------------------------------------------


def project_callers(root: Path) -> list[str]:
    """Every site under ``root`` that CALLS ``.project(...)`` — the retired role's only entry point.

    ⭐ Parsed from the AST rather than grepped, so a mention in a docstring or a comment (there is one
    in ``messages/variance.py``, deliberately, as the formula's provenance) is NOT a false positive.
    That distinction is the whole point: the comment SHOULD stay, and the call must not come back."""
    hits = []
    for f in sorted(Path(root).rglob("*.py")):
        try:
            tree = ast.parse(f.read_text())
        except SyntaxError:  # pragma: no cover - a broken tree is not this gate's business
            continue
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "project"
            ):
                hits.append(f"{f}:{node.lineno}")
    return hits


def npmle_module_present(src_root: Path) -> bool:
    """Is the NPMLE module still on disk at all?

    ⭐ Since 2026-08-21 the answer is NO — it was converge-and-deleted with its last consumer (a
    QC-only landscape fitted on the forbidden ``mass / eff_gdna`` pair). The retirement is therefore
    complete rather than merely unwired, and this is the strongest form the check can take."""
    return (Path(src_root) / "calibration" / "npmle.py").is_file()


# ---------------------------------------------------------------------------
# ⓐ the landscape's per-slot posterior variance
# ---------------------------------------------------------------------------


def posterior_logvar(landscape, count, exposure) -> np.ndarray:
    """``Var(log rho | count, exposure)`` per slot on the fitted grid — the one line the npmle's
    ``project`` needed a mixture for.

    The posterior is the slot's own Poisson kernel times the fitted population density, normalised on
    the grid (⭐ EXACTLY the product ``fit_abundance_landscape`` already integrates to get ``w``, so
    this is the same object read for its second moment rather than a new model). Returned in nats².
    """
    grid10 = np.asarray(landscape.log_rho, np.float64) / np.log(10.0)
    kern = _poisson_kernels(
        np.asarray(count, np.float64), np.maximum(np.asarray(exposure, np.float64), _EPS), grid10
    )
    post = kern * np.exp(np.asarray(landscape.logP, np.float64))[None, :]
    post /= np.maximum(post.sum(axis=1, keepdims=True), _EPS)
    x = np.asarray(landscape.log_rho, np.float64)[None, :]
    mu = (post * x).sum(axis=1, keepdims=True)
    return ((post * (x - mu) ** 2).sum(axis=1)).astype(np.float64)


def mass_weighted_logvar(landscape, lo=None, hi=None) -> float:
    """The fitted density's own variance of ``log rho``, optionally restricted to one basin — a
    POPULATION spread that no slot's count can shrink. This is the shape of quantity a premise is."""
    x = np.asarray(landscape.log_rho, np.float64)
    p = np.exp(np.asarray(landscape.logP, np.float64))
    if lo is not None and hi is not None:
        m = (x >= lo) & (x <= hi)
        x, p = x[m], p[m]
    s = float(p.sum())
    if s <= 0.0:
        return float("nan")
    p = p / s
    mu = float((p * x).sum())
    return float((p * (x - mu) ** 2).sum())


# ---------------------------------------------------------------------------
# one condition
# ---------------------------------------------------------------------------


def measure(bundle: dict, condition: str) -> dict:
    """PURE (arrays in, dict out) so the self-test feeds it synthetics."""
    counts = np.asarray(bundle["counts"], np.float64)
    exposure = np.asarray(bundle["exposure"], np.float64)
    sel = np.asarray(bundle["model_free"], bool) & (exposure > 0.0)
    coarse = np.asarray(bundle["coarse"], np.int64)
    ss, cap = PVO.stratum(condition)

    al = fit_abundance_landscape(bundle["substrate"], bundle["region_arrays"], bundle["wall_mask"])
    if al is None:
        return {"condition": condition, "skipped": "fewer than two model-free regions"}

    # ── ⓐ the landscape's posterior variance against the counting variance it would replace
    c, e = counts[sel], exposure[sel]
    v_post = posterior_logvar(al.landscape, c, e)
    v_count = count_logvar(c)
    ratio = v_post / np.maximum(v_count, _EPS)
    cc = coarse[sel]
    by_class = {}
    for cid, cname in _CLASS_NAMES.items():
        m = cc == cid
        if m.any():
            by_class[cname] = {
                "n": int(m.sum()),
                "median_ratio": float(np.median(ratio[m])),
                "median_v_post": float(np.median(v_post[m])),
                "median_v_count": float(np.median(v_count[m])),
            }

    # ── ⓑ the SHIPPED law at pass-0 init, on the chain's real adjacent pairs
    mass_slot = np.asarray(bundle["mass_slot"], np.float64)
    left = np.asarray(bundle["left"], np.int64)
    splice_in = np.asarray(bundle["splice_in"], bool)
    has = left >= 0
    dst = np.flatnonzero(has)
    src = left[has]
    # at init Var(f_g) = 0 everywhere, so composition_logvar reduces EXACTLY to count_logvar(M).
    # ⛔ ASSERTED, not assumed: the reduction is what makes this arm a reading of the shipped law
    # rather than a lookalike.
    from rigel.calibration.messages.variance import composition_logvar

    lv_init = np.asarray(
        composition_logvar(
            np.ones_like(mass_slot), np.full_like(mass_slot, 300.0),
            np.full_like(mass_slot, 300.0), np.zeros_like(mass_slot), mass_slot,
        ),
        np.float64,
    )
    reduction_holds = bool(np.allclose(lv_init, count_logvar(mass_slot), rtol=0, atol=1e-12))
    s2 = np.asarray(transfer_logvar(lv_init[dst], lv_init[src], splice_in[dst]), np.float64)
    live_hop = np.isfinite(s2)
    non_zero = live_hop & (s2 > 0.0)
    shipped = {
        "n_hops": int(dst.size),
        "n_splice_in_zero": int(np.count_nonzero(splice_in[dst])),
        "reduction_to_counting_holds": reduction_holds,
        "median_s2": float(np.median(s2[non_zero])) if non_zero.any() else float("nan"),
        "p10_s2": float(np.percentile(s2[non_zero], 10)) if non_zero.any() else float("nan"),
        "p90_s2": float(np.percentile(s2[non_zero], 90)) if non_zero.any() else float("nan"),
        "negligible_share": float(np.mean(s2[non_zero] < _NEGLIGIBLE_LOGVAR))
        if non_zero.any()
        else float("nan"),
    }

    # ── ⓒ the heterogeneity: the landscape's population spread vs a premise fitted from the hops
    rho_slot = np.where(
        np.asarray(bundle["exposure_slot"], np.float64) > 0.0,
        mass_slot / np.maximum(np.asarray(bundle["exposure_slot"], np.float64), _EPS),
        np.nan,
    )
    lr = np.log(np.maximum(rho_slot[dst], _EPS)) - np.log(np.maximum(rho_slot[src], _EPS))
    vr = count_logvar(mass_slot[dst]) + count_logvar(mass_slot[src])
    good = np.isfinite(lr) & (rho_slot[dst] > 0) & (rho_slot[src] > 0)
    prem = premise_logvar(lr[good], vr[good]) if int(np.count_nonzero(good)) >= 2 else float("nan")
    hetero = {
        "landscape_population_logvar": mass_weighted_logvar(al.landscape),
        "landscape_depleted_basin_logvar": mass_weighted_logvar(
            al.landscape, al.depleted.lo, al.depleted.hi
        ),
        "premise_from_hops": float(prem),
        "n_hops_fitted": int(np.count_nonzero(good)),
        "median_hop_counting_var": float(np.median(vr[good])) if good.any() else float("nan"),
    }

    return {
        "condition": condition,
        "strand": ss,
        "capture": cap,
        "scope": OC._scope(condition),
        "zero_gdna": bool(PVO.is_zero_gdna(condition)),
        "posterior": {
            "n": int(sel.sum()),
            "median_ratio": float(np.median(ratio)),
            "share_below_one": float(np.mean(ratio < 1.0)),
            "by_class": by_class,
        },
        "shipped": shipped,
        "heterogeneity": hetero,
    }


def load_condition(suite: Path, condition: str, index, cached: dict) -> dict:
    """⚠ **THE SELF-TEST CANNOT REACH THIS FUNCTION** — it does I/O by design, and `measure` is fed
    synthetic arrays instead. So a bank renamed under it is caught only by a REAL run, which is
    `TRAPS: a-green-suite-hid-five-dead-instruments` in its instrument form; its first real run here
    duly failed on `substrate.boundary_unspliced_count` (the banks live behind a `PopulationView`, so
    it is `substrate.boundary_unspliced.count`). Run the instrument, never just its self-test."""
    cache = read_scan_cache(Path(suite) / "scan_cache" / condition, index)
    payload = cache.payload
    ra = cached["region_arrays"]
    substrate = CalibrationSubstrate.from_payload(payload, ra)
    chain = build_region_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    mask = build_region_wall_mask(
        ra,
        cached["mature"],
        cached["reach"][0],
        cached["reach"][1],
        w_max=w_max_from_deposited_lengths(payload.deposited_lengths),
    )
    counts, exposure, model_free = region_counts_and_exposure(substrate, ra, mask)

    # ⭐⭐ THE CHAIN-AXIS PAIR IS THE RELAY'S OWN, TAKEN FROM THE GEOMETRY RATHER THAN REBUILT.
    # `composition_logvar`'s counting term reads the object's total `n` and its opportunity is
    # `eff_gdna` — CONTAINED placements at a REGION and CROSSING placements at a BOUNDARY, one rule and
    # one array (`region_gdna_geometry`). ⛔ A region-only exposure was the first draft and it dropped
    # EVERY hop: the chain alternates region/boundary, so a boundary with no exposure kills all of
    # them, and the premise fit read "0 hops". A boundary has zero genomic measure but a perfectly
    # well-defined crossing opportunity, which is exactly what the geometry supplies.
    from rigel.calibration.region_geometry import build_region_geometry, region_gdna_geometry
    from rigel.scan_cache import calibration_inputs

    kw = calibration_inputs(cache, index)
    geometry = build_region_geometry(
        chain, substrate, ra, cached["sj"], kw["gdna_fl_pmf"], kw["rna_fl_pmf"], None
    )
    mass_slot, exposure_slot = (np.asarray(a, np.float64) for a in region_gdna_geometry(geometry))

    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    is_r = kind == REGION
    n_slots = chain.n_slots
    b_spl = np.asarray(substrate.boundary_spliced.count, np.float64)
    if b_spl.ndim > 1:
        b_spl = b_spl.sum(axis=1)
    splice_in = np.zeros(n_slots, bool)
    splice_in[~is_r] = b_spl[obj[~is_r]] > 0.0

    return {
        "counts": counts,
        "exposure": exposure,
        "model_free": model_free,
        "coarse": coarse_type_array(np.asarray(ra.signature, np.int64)),
        "mass_slot": mass_slot,
        "exposure_slot": exposure_slot,
        "left": np.asarray(chain.left, np.int64),
        "splice_in": splice_in,
        "substrate": substrate,
        "region_arrays": ra,
        "wall_mask": mask,
    }


def _f(x, w=9, p=4):
    return f"{'—':>{w}}" if x is None or not np.isfinite(x) else f"{x:>{w}.{p}f}"


def report(d: dict) -> int:
    if "skipped" in d:
        print(f"\n── {d['condition']}: SKIPPED ({d['skipped']})")
        return 0
    bad = 0
    print(f"\n── {d['condition']}  [{d['strand']} x {d['capture']}]  {d['scope']}")
    p = d["posterior"]
    print(
        f"   ⓐ LANDSCAPE POSTERIOR vs COUNTING variance: median ratio {_f(p['median_ratio'], 7, 4)}"
        f"   share tighter than the count warrants {p['share_below_one']:.3f}   (n={p['n']:,})"
    )
    for cname, s in p["by_class"].items():
        print(
            f"        {cname:>10}: ratio {_f(s['median_ratio'], 7, 4)}   "
            f"v_post {_f(s['median_v_post'], 8, 4)}  v_count {_f(s['median_v_count'], 8, 4)}  "
            f"(n={s['n']:,})"
        )
    s = d["shipped"]
    flag = "✔" if s["reduction_to_counting_holds"] else "⛔"
    if not s["reduction_to_counting_holds"]:
        bad += 1
    print(
        f"   {flag} gate composition_logvar-reduces-to-counting-at-init "
        f"(Var(f_g) = 0 ⇒ the arm below IS the shipped law)"
    )
    print(
        f"   ⓑ SHIPPED sigma^2_transfer at pass-0 init over {s['n_hops']:,} real hops "
        f"({s['n_splice_in_zero']:,} SPLICE-IN rows are identically 0):"
    )
    print(
        f"        p10 {_f(s['p10_s2'], 8, 4)}  median {_f(s['median_s2'], 8, 4)}  "
        f"p90 {_f(s['p90_s2'], 8, 4)} nats^2   "
        f"share below {_NEGLIGIBLE_LOGVAR} nats^2 (essentially undamped): "
        f"{s['negligible_share']:.3f}"
    )
    h = d["heterogeneity"]
    print(
        f"   ⓒ HETEROGENEITY — a variance no count can shrink: landscape population "
        f"{_f(h['landscape_population_logvar'], 8, 4)}   depleted basin "
        f"{_f(h['landscape_depleted_basin_logvar'], 8, 4)}   "
        f"premise fitted from {h['n_hops_fitted']:,} hops {_f(h['premise_from_hops'], 8, 4)} nats^2"
    )
    return bad


def summarise(rows: list[dict]) -> None:
    rows = [r for r in rows if "skipped" not in r]
    if not rows:
        return
    print(f"\n{'=' * 120}")
    print("PER STRATUM — never pooled across strata; the zero controls are their own row")
    print(f"{'=' * 120}")
    hdr = (
        f"{'selection':<44}{'n':>3}  {'post/count':>10} {'undamped':>9}  "
        f"{'ship_med':>9}  {'land_het':>9} {'premise':>8}"
    )
    print(hdr)
    print("-" * len(hdr))
    for name, pred in OC._SELECTIONS:
        sub = [r for r in rows if pred(r["condition"])]
        if not sub:
            continue

        def _m(f, sub=sub):
            v = [f(r) for r in sub]
            v = [x for x in v if x is not None and np.isfinite(x)]
            return float(np.mean(v)) if v else float("nan")

        print(
            f"{name:<44}{len(sub):>3}  "
            f"{_f(_m(lambda r: r['posterior']['median_ratio']), 10, 4)} "
            f"{_f(_m(lambda r: r['shipped']['negligible_share']), 9, 3)}  "
            f"{_f(_m(lambda r: r['shipped']['median_s2']), 9, 4)}  "
            f"{_f(_m(lambda r: r['heterogeneity']['landscape_population_logvar']), 9, 4)} "
            f"{_f(_m(lambda r: r['heterogeneity']['premise_from_hops']), 8, 4)}"
        )
    print(
        "\n⛔ ⓐ's ratio BELOW 1 is the wrong direction for a transfer variance: a fitted population "
        "prior shrinks a slot toward a mode, so substituting it would make messages MORE confident."
        "\n⛔ A TOTAL is composition-vacuous, so nothing measured here licenses a COMPOSITION hop's "
        "variance — only a LEVEL hop's, and the population premise."
    )


def self_test() -> int:  # noqa: C901
    ok = fail = 0

    def check(name, cond):
        nonlocal ok, fail
        if cond:
            ok += 1
        else:
            fail += 1
            print(f"⛔ {name}")

    # ── ① the retirement gate, on the real tree
    src_root = Path(__file__).resolve().parents[2] / "src" / "rigel"
    check("the npmle module is GONE from src/", not npmle_module_present(src_root))
    check("...so nothing under src/ can call `.project()`", project_callers(src_root) == [])
    # ⛔⛔ **THE SCANNER MUST STILL BE ABLE TO SEE A CALL, OR THE GATE ABOVE IS VACUOUS.** It used to
    # prove that by scanning `tests/calibration/`, where `test_npmle.py` called `.project()` 7 times —
    # and that file was DELETED with the module, which would have left this check passing on an empty
    # tree forever (`TRAPS: an-ablation-that-never-ran`). It is proven on a synthetic source instead,
    # so it depends on nothing that is supposed to be gone.
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        (Path(td) / "m.py").write_text("def f(p, a, b):\n    return p.project(a, b)\n")
        seen = project_callers(Path(td))
    check("the AST scanner DOES detect a real `.project()` call", len(seen) == 1)
    with tempfile.TemporaryDirectory() as td:
        # a docstring mention must NOT count — the provenance comment in variance.py is exactly this
        (Path(td) / "m.py").write_text('"""calls p.project(a, b) historically."""\nX = 1\n')
        check("...and does NOT count a mention in prose", project_callers(Path(td)) == [])

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests" / "calibration"))
    from test_abundance_landscape import bimodal_parts, parts

    counts, lengths, sig, *_ = bimodal_parts()
    sub, ra, mask = parts(counts, lengths, sig)
    al = fit_abundance_landscape(sub, ra, mask)
    cts, exp, mf = region_counts_and_exposure(sub, ra, mask)
    m = mf & (exp > 0)

    # ── ⓐ the posterior variance
    v_post = posterior_logvar(al.landscape, cts[m], exp[m])
    v_count = count_logvar(cts[m])
    check("the posterior variance is positive and finite", np.all(np.isfinite(v_post)) and np.all(v_post > 0))
    check(
        "a fitted population prior TIGHTENS the typical slot (the arm's whole point)",
        float(np.median(v_post / v_count)) < 1.0,
    )
    check(
        "a FLAT prior does not tighten below the count's own variance",
        float(
            np.median(
                posterior_logvar(
                    type(al.landscape)(
                        log_rho=al.landscape.log_rho,
                        logP=np.full_like(al.landscape.logP, -np.log(al.landscape.logP.size)),
                        n_train=0,
                    ),
                    cts[m],
                    exp[m],
                )
                / v_count
            )
        )
        > float(np.median(v_post / v_count)),
    )
    check("it is deterministic", np.array_equal(v_post, posterior_logvar(al.landscape, cts[m], exp[m])))

    # ── the population spread
    full = mass_weighted_logvar(al.landscape)
    basin = mass_weighted_logvar(al.landscape, al.depleted.lo, al.depleted.hi)
    check("the population spread is positive", full > 0)
    check("one basin is TIGHTER than the whole bimodal population", basin < full)
    check(
        "restricting to the whole grid reproduces the unrestricted number",
        abs(
            mass_weighted_logvar(al.landscape, al.landscape.log_rho[0], al.landscape.log_rho[-1])
            - full
        )
        < 1e-12,
    )

    # ── ⓑ/ⓒ the whole condition path on a synthetic chain
    n = cts.shape[0]
    mass_slot = cts.astype(np.float64)
    bundle = {
        "counts": cts,
        "exposure": exp,
        "model_free": mf,
        "coarse": coarse_type_array(np.asarray(ra.signature, np.int64)),
        "mass_slot": mass_slot,
        "exposure_slot": exp,
        "left": np.concatenate([[-1], np.arange(n - 1)]).astype(np.int64),
        "splice_in": np.zeros(n, bool),
        "substrate": sub,
        "region_arrays": ra,
        "wall_mask": mask,
    }
    d = measure(bundle, "gdna_g50_ss_0.99_nrna_mid_capture_on")
    check("the init reduction to the counting law HOLDS", d["shipped"]["reduction_to_counting_holds"])
    check("the shipped law is measured over every real hop", d["shipped"]["n_hops"] == n - 1)
    check("the stratum parses", (d["strand"], d["capture"]) == ("stranded", "capture ON"))
    check("the premise is finite and non-negative", d["heterogeneity"]["premise_from_hops"] >= 0.0)
    check("json round-trips", json.dumps(d) is not None)

    # ⛔ a SPLICE IN hop must be identically 0 — the direction-dependence, not an average
    b2 = dict(bundle)
    b2["splice_in"] = np.ones(n, bool)
    d2 = measure(b2, "gdna_g50_ss_0.99_nrna_mid_capture_on")
    check(
        "every SPLICE-IN hop is identically zero", d2["shipped"]["n_splice_in_zero"] == n - 1
    )
    check("...and no non-zero hop survives", not np.isfinite(d2["shipped"]["median_s2"]))

    # ⛔ the deep-count property: doubling every count must LOWER the shipped transfer variance
    b3 = dict(bundle)
    b3["mass_slot"] = mass_slot * 100.0
    d3 = measure(b3, "gdna_g50_ss_0.99_nrna_mid_capture_on")
    check(
        "the shipped law SHRINKS with depth (the structural gap this audit names)",
        d3["shipped"]["median_s2"] < d["shipped"]["median_s2"],
    )
    check(
        "...while the landscape's population spread does NOT shrink with depth",
        abs(
            d3["heterogeneity"]["landscape_population_logvar"]
            - d["heterogeneity"]["landscape_population_logvar"]
        )
        < 1e-12,
    )

    print(f"self-test: {ok} passed, {fail} failed")
    return 1 if fail else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--condition", default=None)
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)
    if args.self_test:
        return self_test()

    src_root = Path(__file__).resolve().parents[2] / "src" / "rigel"
    callers = project_callers(src_root)
    present = npmle_module_present(src_root)
    print("① THE RETIREMENT, VERIFIED FROM THE AST rather than from the comment that claims it")
    if present or callers:
        print(f"   ⛔ npmle module on disk: {present}; `.project()` callers under src/: {callers}")
    else:
        print(
            "   ✔ the npmle module is DELETED and nothing under src/ calls `.project()` — the "
            "sigma^2_transfer role is not merely unwired, it is gone"
        )
        print(
            "   ⭐ and the role it filled is NOT vacant in the tool: `variance.premise_logvar` is the "
            "population term, fitted. It IS vacant in the shipped `relay` policy — see ⓑ/ⓒ."
        )

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
    rows, bad = [], 0
    for c in conds:
        d = measure(load_condition(args.suite, c, index, cached), c)
        bad += report(d)
        rows.append(d)
    summarise(rows)
    if args.json:
        args.json.write_text(json.dumps({"suite": str(args.suite), "rows": rows}))
        print(f"wrote {args.json}")
    return 1 if (bad or callers or present) else 0


if __name__ == "__main__":
    raise SystemExit(main())
