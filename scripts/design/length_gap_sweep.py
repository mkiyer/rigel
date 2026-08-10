#!/usr/bin/env python
"""⭐⭐⭐ **IS THE LENGTH CHANNEL'S ANSWER A FUNCTION OF THE PMF GAP AT ALL? — DRAINED, no solver runs.**

**The question this exists to settle.** At ladder ``g00`` there is no gDNA, the four gDNA pools are
exactly empty, and ``gdna_pmf`` falls back to the EB anchor — yet the channel reports 54-57 % gDNA. The
natural reading is "the residual gap ``Delta = mu_g - mu_r`` is small but nonzero, and the channel reads it
as composition". ⛔ **If that reading is right, shrinking ``Delta`` must shrink the error.** This instrument
shrinks it continuously, over twelve orders of magnitude, and reports what the channel does.

The arm is one line::

    p_alpha = normalise( (1 - alpha) * rna_pmf + alpha * gdna_pmf )

``alpha = 1`` is the shipped channel exactly; ``alpha = 0`` feeds ONE pmf as both components, which is the
channel's established null and must be EXACTLY inert. Everything between is the same channel at a smaller
gap and nothing else changed.

⛔⛔ **IT DRAINS BY DEFAULT, AND THAT IS NOT A DETAIL.** Production drains the side buffer before
calibrating, the buffer holds fragments whose junction fell in the unsequenced mate gap — which requires a
LONG fragment — so draining adds tail mass to the very pool ``rna_pmf`` is fitted on. Measured: it moves
``mu_r`` by +6 to +8 bp and FLIPS the sign of ``Delta`` at ``g00 capture_on`` (+1.21 bp undrained,
-2.57 bp drained). ⚠ A predecessor of this measurement was run undrained and every magnitude in it is
void; ``--undrained`` is retained only so the two can be compared, and it prints as a diagnostic arm.

⭐ **FOUR THINGS ARE REPORTED PER ``alpha`` AND THEY ANSWER DIFFERENT QUESTIONS.** A channel that fails
gracefully has all four fading together; the whole finding is that they do not.

===================  ==========================================================================
 ``med ptp``          the row's AMPLITUDE. Should fade linearly in ``Delta``.
 ``live`` / ``mass``  PARTICIPATION — how many slots the row speaks at. The discrimination guard
                      is an EXACT inequality on five moments, so this may be BINARY in ``Delta``.
 ``impl f_gdna``      the ANSWER, mass-weighted from the per-slot argmax, against a known truth.
 ``med tau``          the DECLARED PRECISION, and the share of it sitting on the grid's own floor
                      ``1/Var(lambda grid)``, which is not evidence.
===================  ==========================================================================

⭐ ``argmax@end`` is the mechanism made visible: a Gaussian log-likelihood is asymptotically LINEAR in
``pi`` with slope proportional to ``Delta``, so as the gap closes the row degenerates to a straight tilt
and a straight tilt's maximum is a grid ENDPOINT. A high value there means the channel is emitting
``sign(Delta' V^-1 (x - N*mu))`` — the DIRECTION of the residual with its magnitude divided out — rather
than a location.

Usage::

    python scripts/design/length_gap_sweep.py                     # the three-condition default
    python scripts/design/length_gap_sweep.py --undrained         # ⛔ diagnostic only
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

#: ⭐ Twelve orders of magnitude, because the claim under test is that the answer is INVARIANT to the gap.
#: A three-point sweep could not distinguish "invariant" from "insensitive"; the endpoints are the two
#: things that must hold exactly (0 -> the established null, 1 -> the shipped channel).
_ALPHA = (0.0, 1e-12, 1e-9, 1e-6, 1e-3, 1e-2, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0)

#: ⛔ g00 is the FAILURE and flgap is the CONTRAST — a real 110 bp gap where the channel demonstrably
#: works. Reading either alone is one-sided.
_DEFAULT = (
    ("ladder", "gdna_g00_ss_0.50_nrna_none_capture_off"),
    ("ladder", "gdna_g00_ss_0.50_nrna_none_capture_on"),
    ("flgap_short", "gdna_g50_ss_0.50_nrna_none_capture_off"),
)

_FAILED: list[str] = []


def _gate(name: str, ok: bool, detail: str) -> None:
    print(f"    {'✅' if ok else '⛔'} {name:<48} {detail}")
    if not ok:
        _FAILED.append(name)


def _mean_bp(pmf) -> float:
    p = np.asarray(pmf, np.float64)
    return float((np.arange(p.shape[0], dtype=np.float64) * p).sum())


def _blend(rna_pmf, gdna_pmf, alpha):
    """``(1-alpha)*rna + alpha*gdna``, renormalised. ⛔ ``alpha = 0`` must return ``rna_pmf`` BIT-FOR-BIT
    and ``alpha = 1`` must return ``gdna_pmf`` bit-for-bit, or the two endpoint gates are testing a
    perturbed channel rather than the shipped one."""
    if alpha == 0.0:
        return np.asarray(rna_pmf, np.float64)
    if alpha == 1.0:
        return np.asarray(gdna_pmf, np.float64)
    b = (1.0 - alpha) * np.asarray(rna_pmf, np.float64) + alpha * np.asarray(gdna_pmf, np.float64)
    s = float(b.sum())
    return b / s if s > 0.0 else b


def load_condition(panel, cond, cfg, index, ra, drain):
    """The condition's payload and its three origin partitions, DRAINED unless told otherwise.

    ⛔ The parts are drained by REPLAYING the whole's choices (``lift_choices``), never independently —
    an independently drained partition resolves held fragments its own way and stops being a partition of
    the same library (``TRAPS: draining-breaks-the-oracle``).
    """
    from rigel.pipeline import _drain_side_buffer, _native_detect_sj_tag
    from rigel.scan_cache import read_scan_cache

    from _oracle import ORIGINS, OracleTruth

    suite = RUNS / "suite" / panel
    bam = str(suite / cond / "sim_oracle.bam")
    cache = suite / "oracle_cache" / cond
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    payload = read_scan_cache(cache / "_main", index, scan).payload
    parts = {k: read_scan_cache(cache / k, index, scan).payload for k in ORIGINS}

    if not drain:
        return payload, OracleTruth.from_parts(payload, parts), "UNDRAINED ⛔ DIAGNOSTIC ARM ONLY", {}

    from rigel.pipeline import scan_and_buffer
    from rigel.second_pass import drain as sp_drain, lift_choices

    n_before = float(np.asarray(payload.deposited_lengths, np.float64).sum())
    _st, sm, _buf, _p = scan_and_buffer(bam, index, scan)
    lift: dict = {}
    payload_d = _drain_side_buffer(payload, index, sm, seed=cfg.second_pass_seed, _lift=lift)
    lifted, n_amb = lift_choices(lift["undrained"], [parts[k] for k in ORIGINS], lift["choices"])
    parts = {
        k: sp_drain(parts[k], ch, node_types=lift["node_types"], junctions=lift["junctions"])
        for k, ch in zip(ORIGINS, lifted)
    }
    n_after = float(np.asarray(payload_d.deposited_lengths, np.float64).sum())
    leak = sum(
        int(np.asarray(getattr(parts["gdna"], b), np.int64).sum())
        for b in ("edge_spliced_count", "sj_count")
    )
    note = (f"DRAINED (production) — deposited {n_before:,.0f} → {n_after:,.0f} "
            f"(+{n_after - n_before:,.0f}), lift_ambiguous {int(n_amb):,}, gdna spliced leak {leak:,}")
    oracle = OracleTruth(full=payload_d, parts=parts,
                         read_counts={k: -1 for k in ORIGINS}, n_ambiguous=int(n_amb))
    return payload_d, oracle, note, {"added": n_after - n_before, "leak": leak}


def main() -> int:  # noqa: C901
    from rigel.calibration.density_deconv import density_factor_precision
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.junction_opportunity import crossing_probability_from_index
    from rigel.calibration.length_likelihood import build_slot_moments, length_loglik
    from rigel.calibration.node_chain import build_node_chain
    from rigel.calibration.node_geometry import build_node_geometry
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.simplex_logodds import _logodds_grid
    from rigel.calibration.splice_graph import build_junction_geometry_arrays
    from rigel.calibration.substrate import CalibrationSubstrate
    from rigel.config import PipelineConfig
    from rigel.index import TranscriptIndex

    from length_channel_census import slot_truth

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--undrained", action="store_true",
                    help="⛔ read the cached UNDRAINED payload. Production DRAINS and the drain moves "
                         "mu_r by 6-8 bp; this arm exists only to compare the two.")
    ap.add_argument("--conditions", nargs="*", default=None, help="panel:condition pairs")
    args = ap.parse_args()

    todo = list(_DEFAULT)
    if args.conditions:
        todo = [tuple(c.split(":", 1)) for c in args.conditions]

    index = TranscriptIndex.load(str(INDEX))
    ra = RegionArrays.from_index(index)
    cfg = PipelineConfig()

    print()
    print("=" * 122)
    print("  ⭐⭐⭐ IS THE LENGTH CHANNEL'S ANSWER A FUNCTION OF THE PMF GAP AT ALL?")
    print("=" * 122)

    for panel, cond in todo:
        if not (RUNS / "suite" / panel / "oracle_cache" / cond / "_main" / "payload.npz").exists():
            print(f"\n  ⚠ SKIP {panel}/{cond} — no cached payload")
            continue
        payload, oracle, note, _info = load_condition(panel, cond, cfg, index, ra, not args.undrained)
        max_w = int(payload.max_length)
        fl = build_fl_models(
            payload,
            junction_opportunity=crossing_probability_from_index(index, max_w),
            gdna_opportunity=gdna_opportunity_from_index(index, max_w),
        )
        chain = build_node_chain(payload.ref_node_offsets, payload.ref_edge_offsets)
        geom = build_node_geometry(
            chain, CalibrationSubstrate.from_payload(payload, ra), ra,
            build_junction_geometry_arrays(index), fl.gdna_pmf, fl.rna_pmf, None,
        )
        lam_grid, fg_grid = _logodds_grid(int(cfg.calibration.sweep_n_grid),
                                          float(cfg.calibration.sweep_logodds_window))
        cnt = np.asarray(geom.unspliced_count, np.float64).sum(axis=1)
        d_obs = np.asarray(geom.unspliced_inv_length_sum, np.float64)
        s_obs = np.asarray(geom.unspliced_length_sum, np.float64)
        fg_true, tot_true = slot_truth(oracle, chain, ra)
        m_r = build_slot_moments(chain, ra, fl.rna_pmf)
        floor = 1.0 / float(np.var(lam_grid))
        mu_g, mu_r = _mean_bp(fl.gdna_pmf), _mean_bp(fl.rna_pmf)
        truth = float(np.nansum(np.where(np.isfinite(fg_true), fg_true * tot_true, 0.0))
                      / max(float(tot_true[np.isfinite(fg_true)].sum()), 1e-9))

        print(f"\n{'─' * 122}")
        print(f"  {panel}/{cond}")
        print(f"  {note}")
        print(f"  mu_g {mu_g:8.3f}   mu_r {mu_r:8.3f}   Delta {mu_g - mu_r:+8.3f} bp   "
              f"n_gdna {fl.n_gdna:,.0f}   n_rna {fl.n_rna:,.0f}   "
              f"slot-truth f_g {truth:.4f}   grid floor {floor:.6f}")

        rows, ll_at = [], {}
        for a in _ALPHA:
            m_g = build_slot_moments(chain, ra, _blend(fl.rna_pmf, fl.gdna_pmf, a))
            ll = length_loglik(m_g, m_r, cnt, d_obs, s_obs, fg_grid)
            ll_at[a] = ll
            ptp = np.ptp(ll, axis=1)
            live = ptp > 0.0
            am = fg_grid[np.argmax(ll, axis=1)]
            tau = density_factor_precision(ll, lam_grid)
            w = tot_true.copy()
            w[~np.isfinite(fg_true)] = 0.0
            lw = w * live
            tot_lw = float(lw.sum())
            rows.append({
                "a": a,
                "delta": _mean_bp(_blend(fl.rna_pmf, fl.gdna_pmf, a)) - mu_r,
                "ptp": float(np.median(ptp[live])) if live.any() else 0.0,
                "n": int(live.sum()),
                "mass": 100.0 * tot_lw / max(float(w.sum()), 1e-9),
                "fg": float((am * lw).sum() / tot_lw) if tot_lw > 0 else float("nan"),
                "tau": float(np.median(tau[live])) if live.any() else 0.0,
                "atfloor": 100.0 * float(np.mean(np.abs(tau[live] - floor) < 1e-4)) if live.any() else 0.0,
                "end": (100.0 * float(lw[(am <= fg_grid[0]) | (am >= fg_grid[-1])].sum()) / tot_lw)
                if tot_lw > 0 else 0.0,
            })

        print("\n    ── GATES ──")
        p0 = np.ptp(ll_at[0.0], axis=1)
        _gate("Ⓖ1 alpha = 0 (one pmf both arms) EXACTLY inert",
              float(np.abs(p0).max()) == 0.0, f"max|ptp| = {float(np.abs(p0).max()):.3e}")
        shipped = length_loglik(build_slot_moments(chain, ra, fl.gdna_pmf), m_r,
                                cnt, d_obs, s_obs, fg_grid)
        _gate("Ⓖ2 alpha = 1 reproduces the shipped channel",
              float(np.abs(ll_at[1.0] - shipped).max()) == 0.0,
              f"max|Δ| = {float(np.abs(ll_at[1.0] - shipped).max()):.3e}")
        _gate("Ⓖ3 the arm could have fired", rows[-1]["n"] > 0,
              f"alpha = 1 live slots {rows[-1]['n']:,}  live mass {rows[-1]['mass']:.2f} %")
        # ⛔ PERTURBATION: blending toward the WRONG pmf must break Ⓖ2. A gate nobody has seen fail is
        # a gate nobody has tested (`perturb-every-gate`).
        bad = length_loglik(build_slot_moments(chain, ra, _blend(fl.gdna_pmf, fl.rna_pmf, 1.0)), m_r,
                            cnt, d_obs, s_obs, fg_grid)
        _gate("Ⓖ4 perturbation: swapped endpoints BREAKS Ⓖ2",
              float(np.abs(bad - shipped).max()) > 0.0,
              f"max|Δ| = {float(np.abs(bad - shipped).max()):.3e}  (must be > 0)")

        print(f"\n    {'alpha':>9}{'Delta bp':>11}{'med ptp':>12}{'live':>9}{'live mass':>11}"
              f"{'impl f_gdna':>13}{'med tau':>10}{'tau@floor':>11}{'argmax@end':>12}")
        print("    " + "-" * 108)
        for r in rows:
            print(f"    {r['a']:>9.0e}{r['delta']:>+11.4f}{r['ptp']:>12.3e}{r['n']:>9,}"
                  f"{r['mass']:>10.2f}%{r['fg']:>13.4f}{r['tau']:>10.4f}{r['atfloor']:>10.1f}%"
                  f"{r['end']:>11.1f}%")
        print(f"    ⭐ truth at these slots is {truth:.4f}. If 'impl f_gdna' is flat down the column, the")
        print("      channel's ANSWER is not a function of the gap and closing the gap cannot fix it.")

    print()
    print("=" * 122)
    if _FAILED:
        print(f"  ⛔ {len(_FAILED)} GATE(S) FAILED — every number above is void until they pass:")
        for f in _FAILED:
            print(f"      {f}")
    else:
        print("  ✅ every gate passed")
    print("=" * 122)
    return 1 if _FAILED else 0


if __name__ == "__main__":
    raise SystemExit(main())
