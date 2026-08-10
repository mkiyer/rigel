#!/usr/bin/env python
"""⭐⭐⭐ **IS THE PRIOR'S CROSSING→FRAGMENT CONVERSION POPULATION-BLIND? — DRAINED, no solver runs.**

**The defect under test.** ``assemble_priors`` turns each line's crossing INCIDENCE into a FRAGMENT count
with one scalar per line (``priors.py``)::

    q = calibration.edge_mass_per_crossing = mass / count      # measured, not modelled
    gdna_edge = mass_gdna_edge * q ;   rna_edge = (mass_rna_edge - spliced) * q

⭐ The logic is right: a fragment crossing ``K`` lines is ``+1`` on each, so an incidence total must be
divided by the mean number of lines crossed, and ``mass/count`` is exactly that.

⛔ **But ``q`` is measured on the POOLED population and applied to the gDNA and RNA parts SEPARATELY.**
Under a uniform field ``q = [min(w-1,a) + min(w-1,b)] / 2(w-1)`` — an explicit function of the fragment
length ``w``. A longer fragment crosses more lines, so it carries LESS mass per crossing. **So gDNA and
RNA have different true ``q`` whenever their length distributions differ, and the assembler gives them
the same one.** It vanishes exactly when the two distributions coincide, which is why the equal-length
ladder cannot see it and why nothing has measured it.

⭐ **THE MEASUREMENT ISOLATES THE DEFECT FROM CALIBRATION'S OWN ERROR.** Feed the assembler a PERFECT
``f_g`` — the origin-split oracle's own per-edge gDNA share — and compare what it then computes against
what the conserved mass says is true::

    A_c = SUM_e count_c[e] * q_pooled[e]      what the assembler computes with a perfect f_g
    T_c = SUM_e mass_c[e]                     the truth: the conserved fragment count, by construction

``A_c / T_c - 1`` is the whole of the defect and nothing else. ⭐ And the composition claim
``phi = A_g/(A_g+A_r)`` against ``phi_true = T_g/(T_g+T_r)`` is what the EM actually consumes.

⛔ **THE NULL IS LOAD-BEARING AND IS RUN FIRST.** The ladder is built with EQUAL configured fragment
lengths, so ``q_g`` must equal ``q_r`` there and the error must vanish. A panel that showed an error on
flgap and none on the ladder is the signature; an error on BOTH means the instrument is measuring
something else (`TRAPS: could-the-arm-have-fired`).

⛔ It DRAINS by default, because production drains and the drain moves the RNA length distribution by
6-8 bp — which is a change to the very quantity this defect is a function of.

Usage::

    python scripts/design/edge_q_population.py                # flgap pair + the ladder null
"""

from __future__ import annotations

import argparse
import dataclasses
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
for _p in (_REPO / "scripts" / "design", _REPO / "tests" / "calibration"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

RUNS = Path.home() / "Downloads" / "rigel_runs"
INDEX = RUNS / "suite" / "rigel_index"

#: ⛔ The flgap pair carries the two LARGEST realised length gaps on disk and in OPPOSITE directions, so a
#: sign flip between them is available as a check. The ladder condition is the NULL: equal configured
#: lengths, so the defect must be absent there or the instrument is measuring something else.
_DEFAULT = (
    ("flgap_short", "gdna_g50_ss_0.50_nrna_none_capture_off"),
    ("flgap_short", "gdna_g50_ss_0.50_nrna_none_capture_on"),
    ("flgap_long", "gdna_g50_ss_0.50_nrna_none_capture_off"),
    ("flgap_long", "gdna_g50_ss_0.50_nrna_none_capture_on"),
    ("ladder", "gdna_g50_ss_0.50_nrna_none_capture_off"),
    ("ladder", "gdna_g50_ss_0.50_nrna_none_capture_on"),
)

_FAILED: list[str] = []


def _gate(name: str, ok: bool, detail: str) -> None:
    print(f"    {'✅' if ok else '⛔'} {name:<50} {detail}")
    if not ok:
        _FAILED.append(name)


def _edge_banks(payload):
    """``(count, mass)`` per contiguous edge, both strand columns summed and the mass descaled.

    ⚠ The descale is ``substrate.INV_LENGTH_SCALE`` because the mass is a fixed point at the SAME scale
    as ``inv_length_sum`` — decoding it any other way here would be a second decoder
    (`TRAPS: a-test-that-redefines`).
    """
    from rigel.calibration.substrate import INV_LENGTH_SCALE

    count = np.asarray(payload.edge_unspliced_count, np.float64)
    if count.ndim == 2:
        count = count.sum(axis=1)
    mass = np.asarray(payload.edge_unspliced_mass, np.float64) / INV_LENGTH_SCALE
    return count, mass


def _q(count, mass):
    """``mass / count``, and 1.0 where nothing crossed — the SHIPPED contract
    (`substrate.PopulationView.mass_per_crossing`): there is no mass at such a line to rescale."""
    out = np.ones(count.shape, np.float64)
    np.divide(mass, count, out=out, where=count > 0.0)
    return out


def load_drained(panel, cond, cfg, index):
    """The drained payload and its three drained origin partitions, exactly as production builds them."""
    from rigel.pipeline import _drain_side_buffer, _native_detect_sj_tag, scan_and_buffer
    from rigel.scan_cache import read_scan_cache
    from rigel.second_pass import drain as sp_drain, lift_choices

    from _oracle import ORIGINS

    suite = RUNS / "suite" / panel
    bam = str(suite / cond / "sim_oracle.bam")
    cache = suite / "oracle_cache" / cond
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    payload = read_scan_cache(cache / "_main", index, scan).payload
    parts = {k: read_scan_cache(cache / k, index, scan).payload for k in ORIGINS}

    _st, sm, _buf, _p = scan_and_buffer(bam, index, scan)
    lift: dict = {}
    payload = _drain_side_buffer(payload, index, sm, seed=cfg.second_pass_seed, _lift=lift)
    lifted, n_amb = lift_choices(lift["undrained"], [parts[k] for k in ORIGINS], lift["choices"])
    parts = {
        k: sp_drain(parts[k], ch, node_types=lift["node_types"], junctions=lift["junctions"])
        for k, ch in zip(ORIGINS, lifted)
    }
    return payload, parts, int(n_amb)


def main() -> int:  # noqa: C901
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.junction_opportunity import crossing_probability_from_index
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.substrate import CalibrationSubstrate
    from rigel.config import PipelineConfig
    from rigel.index import TranscriptIndex

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--conditions", nargs="*", default=None, help="panel:condition pairs")
    args = ap.parse_args()

    todo = [tuple(c.split(":", 1)) for c in args.conditions] if args.conditions else list(_DEFAULT)
    index = TranscriptIndex.load(str(INDEX))
    ra = RegionArrays.from_index(index)
    cfg = PipelineConfig()

    print()
    print("=" * 116)
    print("  ⭐⭐⭐ IS THE PRIOR'S CROSSING→FRAGMENT CONVERSION POPULATION-BLIND?   (DRAINED, no solver)")
    print("=" * 116)
    print(f"\n  {'condition':<40}{'mu_g-mu_r':>10}{'q_g':>7}{'q_r':>7}{'edge%':>8}"
          f"{'EDGE err':>11}{'TOTAL err':>11}{'d phi edge':>10}{'d phi TOT':>10}")
    print("  " + "-" * 114)

    rows = []
    for panel, cond in todo:
        if not (RUNS / "suite" / panel / "oracle_cache" / cond / "_main" / "payload.npz").exists():
            print(f"  ⚠ SKIP {panel}/{cond} — no cached payload")
            continue
        payload, parts, _n_amb = load_drained(panel, cond, cfg, index)

        c_all, m_all = _edge_banks(payload)
        c_g, m_g = _edge_banks(parts["gdna"])
        c_r = np.zeros_like(c_all)
        m_r = np.zeros_like(m_all)
        for k in ("mrna", "nrna"):
            a, b = _edge_banks(parts[k])
            c_r += a
            m_r += b

        q_pool = _q(c_all, m_all)

        # ── GATES, before any number is read ───────────────────────────────────────────────────────
        if not rows:
            print()
        sub = CalibrationSubstrate.from_payload(payload, ra)
        shipped_q = np.asarray(sub.edge_unspliced.mass_per_crossing, np.float64)
        _gate(f"Ⓖ1 pooled q reproduces the shipped one — {cond[:22]}",
              float(np.abs(q_pool - shipped_q).max()) == 0.0,
              f"max|Δ| = {float(np.abs(q_pool - shipped_q).max()):.3e}")
        _gate(f"Ⓖ2 the partition sums to the whole — {cond[:26]}",
              float(np.abs((c_g + c_r) - c_all).max()) == 0.0
              and float(np.abs((m_g + m_r) - m_all).max()) < 1e-6,
              f"count max|Δ| {float(np.abs((c_g + c_r) - c_all).max()):.3e}  "
              f"mass max|Δ| {float(np.abs((m_g + m_r) - m_all).max()):.3e}")

        # ── the defect, with a PERFECT f_g ────────────────────────────────────────────────────────
        # ⭐ The NODE term needs no conversion (a contained fragment deposits on exactly one node), so it
        # enters both arms identically and DILUTES the edge error. The EM sees the total, so the total is
        # what must be reported beside the edge-only figure — an edge-only number overstates the defect
        # by exactly the node share.
        def _node_count(p):
            v = np.asarray(p.node_contained_count, np.float64)
            return float(v.sum(axis=1).sum() if v.ndim == 2 else v.sum())

        n_g = _node_count(parts["gdna"])
        n_r = sum(_node_count(parts[k]) for k in ("mrna", "nrna"))

        a_g, t_g = float((c_g * q_pool).sum()), float(m_g.sum())
        a_r, t_r = float((c_r * q_pool).sum()), float(m_r.sum())
        tot_a_g, tot_t_g = n_g + a_g, n_g + t_g
        tot_a_r, tot_t_r = n_r + a_r, n_r + t_r
        edge_share = t_g / tot_t_g if tot_t_g > 0 else float("nan")
        err_tot_g = tot_a_g / tot_t_g - 1.0 if tot_t_g > 0 else float("nan")
        phi_tot_a = tot_a_g / (tot_a_g + tot_a_r) if (tot_a_g + tot_a_r) > 0 else float("nan")
        phi_tot_t = tot_t_g / (tot_t_g + tot_t_r) if (tot_t_g + tot_t_r) > 0 else float("nan")
        err_g = a_g / t_g - 1.0 if t_g > 0 else float("nan")
        err_r = a_r / t_r - 1.0 if t_r > 0 else float("nan")
        phi_a = a_g / (a_g + a_r) if (a_g + a_r) > 0 else float("nan")
        phi_t = t_g / (t_g + t_r) if (t_g + t_r) > 0 else float("nan")

        # count-weighted mean q per component — the mechanism's own variable
        qg = t_g / float(c_g.sum()) if c_g.sum() > 0 else float("nan")
        qr = t_r / float(c_r.sum()) if c_r.sum() > 0 else float("nan")
        qp = float(m_all.sum()) / float(c_all.sum()) if c_all.sum() > 0 else float("nan")

        max_w = int(payload.max_length)
        fl = build_fl_models(
            payload,
            junction_opportunity=crossing_probability_from_index(index, max_w),
            gdna_opportunity=gdna_opportunity_from_index(index, max_w),
        )
        w = np.arange(fl.gdna_pmf.shape[0], dtype=np.float64)
        gap = float((w * fl.gdna_pmf).sum()) - float((w * fl.rna_pmf).sum())

        rows.append((f"{panel}/{cond}", gap, qg, qr, qp, err_g, err_r, phi_a - phi_t))
        print(f"  {panel + '/' + cond.replace('nrna_none_', ''):<40}{gap:>+10.2f}{qg:>7.4f}{qr:>7.4f}"
              f"{100 * edge_share:>8.1f}%{100 * err_g:>10.3f}%{100 * err_tot_g:>10.3f}%"
              f"{phi_a - phi_t:>+10.5f}{phi_tot_a - phi_tot_t:>+10.5f}")

    print("  " + "-" * 114)
    print("  ⭐ A_c = SUM_e count_c*q_pooled (what the assembler computes with a PERFECT f_g);")
    print("     T_c = SUM_e mass_c (the conserved truth). `d phi` is the composition error the EM sees.")
    print("  ⛔ THE NULL: the ladder rows have equal configured lengths, so q_g must equal q_r and every")
    print("     error must vanish there. An error on the ladder too means this measures something else.")

    # ── Ⓖ3 the PERTURBATION: swapping the two components' q must move the answer ──────────────────
    if rows:
        print()
        biggest = max(rows, key=lambda r: abs(r[1]))
        _gate("Ⓖ3 the flgap arm HAS a length gap to be sensitive to",
              abs(biggest[1]) > 1.0, f"largest |mu_g-mu_r| = {abs(biggest[1]):.2f} bp on {biggest[0]}")
        # ⛔ The null is the ladder CAPTURE-OFF arm only. The ladder's CONFIGURED lengths are equal, but
        # hybrid capture SELECTS on length and manufactures a realised gap of its own (+15.33 bp measured,
        # consistent with the recorded +13.6 bp capture defect in the gDNA pmf) — so ladder capture_on is
        # not an equal-length condition and must not be read as the null.
        nulls = [r for r in rows if "ladder" in r[0] and "capture_off" in r[0]]
        if nulls:
            _gate("Ⓖ4 the NULL (ladder capture-OFF) is equal-length",
                  max(abs(r[1]) for r in nulls) < 8.0,
                  f"|mu_g-mu_r| = {max(abs(r[1]) for r in nulls):.2f} bp; "
                  f"⚠ capture_on is NOT a null, capture manufactures a gap")

    print()
    print("=" * 116)
    if _FAILED:
        print(f"  ⛔ {len(_FAILED)} GATE(S) FAILED — every number above is void until they pass:")
        for f in _FAILED:
            print(f"      {f}")
    else:
        print("  ✅ every gate passed")
    print("=" * 116)
    return 1 if _FAILED else 0


if __name__ == "__main__":
    raise SystemExit(main())
