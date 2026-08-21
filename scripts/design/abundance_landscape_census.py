"""WHAT DOES THE TOTAL-DENSITY FIELD LOOK LIKE, PER CONDITION — the AbundanceLandscape census.

⭐⭐⭐ **A MEASUREMENT: no solver runs, no EM, nothing in `src/` is patched.** Per condition it fits
`calibration.abundance_landscape.fit_abundance_landscape` on the cached payload's wall-exact measured
totals and prints the CENSUS: every mode with its basin mass and width, `rho_0` (the depleted mode),
the span `R`, the pooled-anchor consistency verdict with its gap in nats, and the per-region
enrichment responsibility `w` summarised BY REGION CLASS (intergenic / intron / exon — mean and the
exposure-weighted mean, never a thresholded count). This is the instrument that says whether the
landscape can serve as the pass-0 measured reference: capture-OFF rows must read unimodal, capture-ON
rows bimodal with the enriched basin claiming the exon class, and the depleted mode must agree with
the pooled intergenic anchors on EVERY row.

⛔ **`w` CONFLATES ENRICHMENT WITH EXPRESSION, and every reading here carries that**: a hot unprobed
exon and a probed cold one can share a basin. The failure direction for any consumer is therefore
PERMISSIVE; the pricing of that confound is the plan's rung-5 read, not this instrument's.

⚠ The knn/width constants inherited from `landscape` have only ever been validated on gDNA-shaped
data — that module's own warning, repeated on every page this instrument writes.

`--json` dumps everything the census printed PLUS the fitted curve and a deterministic rug sample per
condition, so a report can be RENDERED elsewhere without refitting (rendering is not this
instrument's job — `benchmark_report.py`'s division of labour).

Usage::

    python scripts/design/abundance_landscape_census.py                      # every cached condition
    python scripts/design/abundance_landscape_census.py --condition NAME     # one condition
    python scripts/design/abundance_landscape_census.py --suite ... --index ...   # another panel
    python scripts/design/abundance_landscape_census.py --json out.json      # census + curves + rug
    python scripts/design/abundance_landscape_census.py --self-test          # perturbed, no I/O
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

from rigel.calibration.abundance_landscape import fit_abundance_landscape  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import RegionType, coarse_type_array  # noqa: E402
from rigel.calibration.splice_graph import (  # noqa: E402
    build_contiguous_boundary_reach_arrays,
    build_mature_wall_distances,
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

#: rug sample size per condition — a DISPLAY budget (the census itself uses every region).
#: Deterministic: an even stride over the class-sorted centres, never a random draw.
_RUG_MAX = 400

_CLASS_NAMES = {
    int(RegionType.INTERGENIC): "intergenic",
    int(RegionType.INTRON): "intron",
    int(RegionType.EXON): "exon",
}


def measure(substrate, region_arrays, wall_mask, condition: str) -> dict | None:
    """One condition's census. PURE (arrays in, dict out) so the self-test feeds it synthetics."""
    al = fit_abundance_landscape(substrate, region_arrays, wall_mask)
    if al is None:
        return None
    counts, exposure, model_free = region_counts_and_exposure(substrate, region_arrays, wall_mask)
    sel = np.asarray(model_free, bool) & (exposure > 0.0)
    coarse = coarse_type_array(np.asarray(region_arrays.signature, np.int64))
    ss, cap = PVO.stratum(condition)

    by_class = {}
    for cid, cname in _CLASS_NAMES.items():
        m = sel & (coarse == cid)
        if not m.any():
            continue
        w = al.w_slot[m]
        by_class[cname] = {
            "n": int(m.sum()),
            "mean_w": float(np.mean(w)),
            "exposure_weighted_w": float(np.average(w, weights=exposure[m])),
            "share_of_selected_exposure": float(exposure[m].sum() / exposure[sel].sum()),
        }

    centres = np.log(np.maximum(counts[sel], 1.0)) - np.log(exposure[sel])
    order = np.argsort(centres, kind="stable")
    stride = max(1, int(np.ceil(order.size / _RUG_MAX)))
    pick = order[::stride]
    sel_idx = np.flatnonzero(sel)

    return {
        "condition": condition,
        "strand": ss,
        "capture": cap,
        "scope": OC._scope(condition),
        "zero_gdna": bool(PVO.is_zero_gdna(condition)),
        "n_train": al.n_train,
        "n_modes": len(al.modes),
        "modes": [
            {
                "log_rho": m.log_rho,
                "basin_mass": m.basin_mass,
                "width": m.width,
                "lo": m.lo,
                "hi": m.hi,
            }
            for m in al.modes
        ],
        "rho_0": al.rho_0,
        "span_R": al.span_R,
        "unimodal": al.enriched is None,
        "anchor_log_rho": al.anchor_log_rho,
        "anchor_gap_nats": al.anchor_gap_nats,
        "anchor_consistent": bool(al.anchor_consistent),
        "depleted_log_rho": al.depleted.log_rho,
        "enriched_log_rho": None if al.enriched is None else al.enriched.log_rho,
        "enriched_basin": None
        if al.enriched is None
        else {"lo": al.enriched.lo, "hi": al.enriched.hi, "mass": al.enriched.basin_mass},
        "w_by_class": by_class,
        "curve": {
            "log_rho": [float(v) for v in al.landscape.log_rho],
            "logP": [float(v) for v in al.landscape.logP],
        },
        "rug": [
            {
                "log_rho": float(centres[i]),
                "region_class": _CLASS_NAMES.get(int(coarse[sel_idx[i]]), "other"),
                "w": float(al.w_slot[sel_idx[i]]),
            }
            for i in pick
        ],
    }


def report(d: dict) -> None:
    flag = "✔" if d["anchor_consistent"] else "⛔"
    shape = "UNIMODAL" if d["unimodal"] else f"bimodal R={d['span_R']:.1f}"
    print(f"\n── {d['condition']}  [{d['strand']} x {d['capture']}]  {d['scope']}")
    print(
        f"   {shape}  rho_0={d['rho_0']:.4g}  n_train={d['n_train']:,}  "
        f"{flag} anchor gap {d['anchor_gap_nats']:.3f} nats"
    )
    for m in d["modes"]:
        print(
            f"      mode log_rho={m['log_rho']:+8.3f}  basin_mass={m['basin_mass']:.4f}  "
            f"width={m['width']:.3f}"
        )
    for cname, s in d["w_by_class"].items():
        print(
            f"      w[{cname:>10}]  mean {s['mean_w']:.4f}  exposure-weighted "
            f"{s['exposure_weighted_w']:.4f}  (n={s['n']:,})"
        )


def load_condition(suite: Path, condition: str, index, cached) -> tuple:
    cache = read_scan_cache(Path(suite) / "scan_cache" / condition, index)
    payload = cache.payload
    ra = cached["region_arrays"]
    substrate = CalibrationSubstrate.from_payload(payload, ra)
    mask = build_region_wall_mask(
        ra,
        cached["mature"],
        cached["reach"][0],
        cached["reach"][1],
        w_max=w_max_from_deposited_lengths(payload.deposited_lengths),
    )
    return substrate, ra, mask


# ---------------------------------------------------------------------------
# self-test — perturbed, no I/O
# ---------------------------------------------------------------------------


def self_test() -> int:
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

    counts, lengths, sig, rho_lo, rho_hi = bimodal_parts()
    sub, ra, mask = parts(counts, lengths, sig)
    d = measure(sub, ra, mask, "gdna_g50_ss_0.99_nrna_mid_capture_on")
    check("a bimodal fixture reads bimodal", d is not None and not d["unimodal"])
    check(
        "the span is the mode ratio", abs(d["span_R"] - rho_hi / rho_lo) / (rho_hi / rho_lo) < 0.5
    )
    check("the stratum parses", (d["strand"], d["capture"]) == ("stranded", "capture ON"))
    check("basin masses sum to one", abs(sum(m["basin_mass"] for m in d["modes"]) - 1.0) < 1e-9)
    check(
        "the exon class carries the enriched responsibility",
        d["w_by_class"]["exon"]["mean_w"] > 0.5 > d["w_by_class"]["intergenic"]["mean_w"],
    )
    check("the rug is bounded and non-empty", 0 < len(d["rug"]) <= _RUG_MAX)
    check("the curve is the full grid", len(d["curve"]["log_rho"]) == len(d["curve"]["logP"]) > 0)
    check("anchors agree on the clean fixture", d["anchor_consistent"])

    c2, l2, s2 = unimodal_parts()
    sub2, ra2, mask2 = parts(c2, l2, s2)
    d2 = measure(sub2, ra2, mask2, "gdna_g50_ss_0.50_nrna_mid_capture_off")
    check(
        "a unimodal fixture reads unimodal with span exactly 1",
        d2["unimodal"] and d2["span_R"] == 1.0,
    )
    check(
        "unimodal w is 0 in every class", all(v["mean_w"] == 0.0 for v in d2["w_by_class"].values())
    )

    # ⛔ the census must MOVE when the field does — a reporter that cannot move reports nothing
    c3 = c2.copy()
    c3[-400:] = c3[-400:] * 300  # lift the WHOLE exon class (the fixture's last 400 regions) 300x
    # ⚠ a first draft lifted only the last 200 — HALF the class — and read mean_w 0.498, which is the
    # census being right about a half-lifted class, not a defect
    sub3, ra3, mask3 = parts(c3, l2, s2)
    d3 = measure(sub3, ra3, mask3, "gdna_g50_ss_0.50_nrna_mid_capture_off")
    check("lifting the exon class creates the enriched mode", not d3["unimodal"])
    check("...and the exon class claims it", d3["w_by_class"]["exon"]["mean_w"] > 0.9)

    # the JSON round-trips
    check("the dict is json-serialisable", json.dumps(d) is not None)

    print(f"self-test: {ok} passed, {fail} failed")
    return 1 if fail else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--condition", default=None)
    ap.add_argument("--json", type=Path, default=None, help="census + curves + rug, one file")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)
    if args.self_test:
        return self_test()

    index = TranscriptIndex.load(args.index)
    ra = RegionArrays.from_index(index)
    cached = {
        "region_arrays": ra,
        "mature": build_mature_wall_distances(index, ra),
        "reach": build_contiguous_boundary_reach_arrays(index),
    }
    conds = (
        [args.condition]
        if args.condition
        else sorted(p.name for p in (args.suite / "scan_cache").iterdir() if p.is_dir())
    )
    rows = []
    for c in conds:
        substrate, ra_c, mask = load_condition(args.suite, c, index, cached)
        d = measure(substrate, ra_c, mask, c)
        if d is None:
            print(f"⛔ {c}: fewer than two model-free regions — nothing to fit")
            continue
        report(d)
        rows.append(d)

    print(
        "\n⚠ the knn/width constants are inherited from `landscape` and were validated on "
        "gDNA-shaped data only; ⚠ w conflates enrichment with expression (permissive direction)."
    )
    if args.json:
        args.json.write_text(json.dumps({"suite": str(args.suite), "rows": rows}))
        print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
