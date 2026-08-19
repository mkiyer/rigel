#!/usr/bin/env python
"""⭐⭐⭐ **IS THE PRIOR'S CROSSING→FRAGMENT CONVERSION POPULATION-BLIND? — DRAINED, no solver runs.**

**The defect under test.** ``assemble_priors`` turns each boundary's crossing INCIDENCE into a FRAGMENT count
with one scalar per boundary (``priors.py``)::

    q = calibration.boundary_mass_per_crossing = mass / count      # measured, not modelled
    gdna_boundary = mass_gdna_boundary * q ;   rna_boundary = (mass_rna_boundary - spliced) * q

⭐ The logic is right: a fragment crossing ``K`` boundaries is ``+1`` on each, so an incidence total must be
divided by the mean number of boundaries crossed, and ``mass/count`` is exactly that.

⛔ **But ``q`` is measured on the POOLED population and applied to the gDNA and RNA parts SEPARATELY.**
Under a uniform field ``q = [min(w-1,a) + min(w-1,b)] / 2(w-1)`` — an explicit function of the fragment
length ``w``. A longer fragment crosses more boundaries, so it carries LESS mass per crossing. **So gDNA and
RNA have different true ``q`` whenever their length distributions differ, and the assembler gives them
the same one.** It vanishes exactly when the two distributions coincide, which is why the equal-length
ladder cannot see it and why nothing has measured it.

⭐ **THE MEASUREMENT ISOLATES THE DEFECT FROM CALIBRATION'S OWN ERROR.** Feed the assembler a PERFECT
``f_g`` — the origin-split oracle's own per-boundary gDNA share — and compare what it then computes against
what the conserved mass says is true::

    A_c = SUM_e count_c[e] * q_pooled[e]      what the assembler computes with a perfect f_g
    T_c = SUM_e mass_c[e]                     the truth: the conserved fragment count, by construction

``A_c / T_c - 1`` is the whole of the defect and nothing else. ⭐ And the composition claim
``phi = A_g/(A_g+A_r)`` against ``phi_true = T_g/(T_g+T_r)`` is what the EM actually consumes.

⛔ **THE NULL IS LOAD-BEARING AND IS RUN FIRST.** The ladder is built with EQUAL configured fragment
lengths, so ``q_g`` must equal ``q_r`` there and the error must vanish. A panel that showed an error at a
length GAP and none at the null is the signature; an error on BOTH means the instrument is measuring
something else (`TRAPS: could-the-arm-have-fired`).

⛔⛔ **HALF THE RECORDED VERDICT IS NO LONGER ESTABLISHED BY ANYTHING ON DISK, AND NO ARM HERE CAN
RESTORE IT (2026-08-17).** The published verdict is *"the defect is real, is NOT driven by fragment
length (the equal-length null shows ``q_g`` 0.633 vs ``q_r`` 0.523 …), and is bounded at ≤ 0.6 pp of
composition"*. ⭐ **The equal-length null half REPRODUCES EXACTLY** — re-run 2026-08-17, ``q_g`` 0.6330 vs
``q_r`` 0.5233 — so the mechanism claim stands. ⛔ **The ≤ 0.6 pp BOUND does not.** It was established
across two OPPOSITE-SIGN length gaps carried by panels deleted on 2026-08-13; the two rows that survive
read ``d phi TOTAL`` **+0.00095** (null) and **+0.00244** (capture gap), which are consistent with the
bound but are not evidence FOR it, because the conditions that would stress it are the ones that are
missing. ⚠ **Do not re-quote the bound and do not invent a replacement.** ⛔ Nor is the answer to keep a
disabled arm here waiting for a substrate: this file carried four hard-coded panel rows that could only
ever SKIP, and they were deleted on 2026-08-17 with the panels and their configs. Restoring the
opposite-sign cross-check means designing a length-gap panel, not re-adding a name.

⭐ **What DOES run is the null plus ONE weaker gap arm, and the difference is stated on every row.**
The ladder's CONFIGURED lengths are equal, but hybrid capture SELECTS on length and manufactures a
realised gap of its own (`TRAPS: configured-lengths-are-not-realised`) — so ladder capture-OFF is the
null and ladder capture-ON is a genuine, single-signed, capture-sized length gap. The opposite-sign
cross-check is what was lost, not the arm.

⛔ It DRAINS by default, because production drains and the drain moves the RNA length distribution by
6-8 bp — which is a change to the very quantity this defect is a function of.

Usage::

    python scripts/design/boundary_q_population.py   # the ladder null + the capture gap arm

⛔⛔ **A SINGLE CONDITION IS NOT A RUN OF THIS INSTRUMENT, AND BOTH GATES NOW SAY SO.** The null alone
cannot see the defect and the gap arm alone cannot show the defect is not the instrument, so Ⓖ3 and Ⓖ4
each FAIL when their own arm is missing — measured 2026-08-17, `--conditions ladder:…capture_off` exits
1 on Ⓖ3 and `--conditions ladder:…capture_on` exits 1 on Ⓖ4. Pass the PAIR::

    python scripts/design/boundary_q_population.py --conditions \\
        ladder:gdna_g50_ss_0.50_nrna_none_capture_off ladder:gdna_g50_ss_0.50_nrna_none_capture_on
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

#: ⭐ Capture-OFF is the NULL (equal configured lengths, so the defect must vanish); capture-ON is the
#: only length-gap arm there is — capture SELECTS on length and manufactures a realised gap
#: (`TRAPS: configured-lengths-are-not-realised`). ⛔ Both are REQUIRED: Ⓖ3 and Ⓖ4 below each fail when
#: their own arm is absent, because neither alone is a run of this instrument.
_DEFAULT = (
    ("ladder", "gdna_g50_ss_0.50_nrna_mid_capture_off"),
    ("ladder", "gdna_g50_ss_0.50_nrna_mid_capture_on"),
)

_FAILED: list[str] = []


def _gate(name: str, ok: bool, detail: str) -> None:
    print(f"    {'✅' if ok else '⛔'} {name:<50} {detail}")
    if not ok:
        _FAILED.append(name)


def _boundary_banks(payload):
    """``(count, mass)`` per contiguous boundary, both strand columns summed.

    ⛔ **THERE IS NOTHING TO DECODE.** This helper used to divide the mass by
    ``substrate.INV_LENGTH_SCALE`` (2^32) because the bank was FIXED POINT. The whole fixed-point layer
    went at `94d283c0` under ONE NUMERIC CONVENTION — a count is an integer, a fraction is float64, and
    nothing decodes a bank — so the constant no longer exists and importing it is what killed this
    instrument (`TRAPS: a-green-suite-hid-five-dead-instruments`). ⚠ The repair is to DELETE the
    division, not to re-derive the scale: ``boundary_unspliced_mass`` is declared float64 in
    `scan_payload.py` and is already in fragment units. Ⓖ1 below is what proves it — the pooled ``q``
    built here must reproduce the SHIPPED ``mass_per_crossing`` bit for bit, and a stray 2^32 could not.
    """
    count = np.asarray(payload.boundary_unspliced_count, np.float64)
    if count.ndim == 2:
        count = count.sum(axis=1)
    mass = np.asarray(payload.boundary_unspliced_mass, np.float64)
    return count, mass


def _q(count, mass):
    """``mass / count``, and 1.0 where nothing crossed — the SHIPPED contract
    (`substrate.PopulationView.mass_per_crossing`): there is no mass at such a boundary to rescale."""
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
        k: sp_drain(parts[k], ch, region_types=lift["region_types"], sj=lift["sj"])
        for k, ch in zip(ORIGINS, lifted)
    }
    return payload, parts, int(n_amb)


def main() -> int:  # noqa: C901
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.sj_opportunity import crossing_probability_from_index
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
    print("=" * 131)
    print("  ⭐⭐⭐ IS THE PRIOR'S CROSSING→FRAGMENT CONVERSION POPULATION-BLIND?   (DRAINED, no solver)")
    print("=" * 131)
    # ⚠ EVERY WIDTH HERE MUST MATCH THE ROW BELOW. They did not: `q_r`, `boundary%` and `BOUNDARY err`
    # were given widths NARROWER than their own labels, so the header printed as `q_rboundary%BOUNDARY
    # err` while the numbers underneath were correctly spaced — a header that cannot be read off the
    # column it names.
    print(f"\n  {'condition':<40}{'mu_g-mu_r':>10}{'q_g':>8}{'q_r':>8}{'boundary%':>12}"
          f"{'BOUNDARY err':>13}{'TOTAL err':>13}{'d phi BOUNDARY':>15}{'d phi TOTAL':>12}")
    print("  " + "-" * 131)

    rows = []
    skipped: list[str] = []
    for panel, cond in todo:
        if not (RUNS / "suite" / panel / "oracle_cache" / cond / "_main" / "payload.npz").exists():
            print(f"  ⚠ SKIP {panel}/{cond} — no cached payload (`panel.py cache` builds it)")
            skipped.append(f"{panel}/{cond}")
            continue
        payload, parts, _n_amb = load_drained(panel, cond, cfg, index)

        c_all, m_all = _boundary_banks(payload)
        c_g, m_g = _boundary_banks(parts["gdna"])
        c_r = np.zeros_like(c_all)
        m_r = np.zeros_like(m_all)
        for k in ("mrna", "nrna"):
            a, b = _boundary_banks(parts[k])
            c_r += a
            m_r += b

        q_pool = _q(c_all, m_all)

        # ── GATES, before any number is read ───────────────────────────────────────────────────────
        if not rows:
            print()
        sub = CalibrationSubstrate.from_payload(payload, ra)
        shipped_q = np.asarray(sub.boundary_unspliced.mass_per_crossing, np.float64)
        _gate(f"Ⓖ1 pooled q reproduces the shipped one — {cond[:22]}",
              float(np.abs(q_pool - shipped_q).max()) == 0.0,
              f"max|Δ| = {float(np.abs(q_pool - shipped_q).max()):.3e}")
        _gate(f"Ⓖ2 the partition sums to the whole — {cond[:26]}",
              float(np.abs((c_g + c_r) - c_all).max()) == 0.0
              and float(np.abs((m_g + m_r) - m_all).max()) < 1e-6,
              f"count max|Δ| {float(np.abs((c_g + c_r) - c_all).max()):.3e}  "
              f"mass max|Δ| {float(np.abs((m_g + m_r) - m_all).max()):.3e}")

        # ── the defect, with a PERFECT f_g ────────────────────────────────────────────────────────
        # ⭐ The REGION term needs no conversion (a contained fragment deposits on exactly one region), so it
        # enters both arms identically and DILUTES the boundary error. The EM sees the total, so the total is
        # what must be reported beside the boundary-only figure — an boundary-only number overstates the defect
        # by exactly the region share.
        def _region_count(p):
            v = np.asarray(p.region_contained_count, np.float64)
            return float(v.sum(axis=1).sum() if v.ndim == 2 else v.sum())

        n_g = _region_count(parts["gdna"])
        n_r = sum(_region_count(parts[k]) for k in ("mrna", "nrna"))

        a_g, t_g = float((c_g * q_pool).sum()), float(m_g.sum())
        a_r, t_r = float((c_r * q_pool).sum()), float(m_r.sum())
        tot_a_g, tot_t_g = n_g + a_g, n_g + t_g
        tot_a_r, tot_t_r = n_r + a_r, n_r + t_r
        boundary_share = t_g / tot_t_g if tot_t_g > 0 else float("nan")
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
            sj_opportunity=crossing_probability_from_index(index, max_w),
            gdna_opportunity=gdna_opportunity_from_index(index, max_w),
        )
        w = np.arange(fl.gdna_pmf.shape[0], dtype=np.float64)
        gap = float((w * fl.gdna_pmf).sum()) - float((w * fl.rna_pmf).sum())

        rows.append((f"{panel}/{cond}", gap, qg, qr, qp, err_g, err_r, phi_a - phi_t))
        print(f"  {panel + '/' + cond.replace('nrna_none_', ''):<40}{gap:>+10.2f}{qg:>8.4f}{qr:>8.4f}"
              f"{100 * boundary_share:>11.1f}%{100 * err_g:>12.3f}%{100 * err_tot_g:>12.3f}%"
              f"{phi_a - phi_t:>+15.5f}{phi_tot_a - phi_tot_t:>+12.5f}")

    print("  " + "-" * 131)
    print("  ⭐ A_c = SUM_e count_c*q_pooled (what the assembler computes with a PERFECT f_g);")
    print("     T_c = SUM_e mass_c (the conserved truth). `d phi` is the composition error the EM sees.")
    print("  ⛔ THE NULL is ladder CAPTURE-OFF: equal configured lengths, so q_g must equal q_r and every")
    print("     error must vanish there. An error on the null too means this measures something else.")
    print("  ⚠ ladder capture-ON is NOT a null — capture SELECTS on length and manufactures a realised")
    print("     gap — but it is single-signed, so no row here can show a SIGN FLIP.")

    # ── Ⓖ3 the SENSITIVITY arm: without a realised length gap there is nothing for this to see ────
    if rows:
        print()
        # ⛔ The null is the ladder CAPTURE-OFF arm only. The ladder's CONFIGURED lengths are equal, but
        # hybrid capture SELECTS on length and manufactures a realised gap of its own (+15.33 bp measured,
        # consistent with the recorded +13.6 bp capture defect in the gDNA pmf) — so ladder capture_on is
        # not an equal-length condition and must not be read as the null.
        nulls = [r for r in rows if "ladder" in r[0] and "capture_off" in r[0]]
        # ⛔⛔ THIS GATE USED TO NAME A PANEL THAT WAS NOT ON DISK, and it kept passing on rows from
        # somewhere else entirely. ⛔ Worse, its candidate set included the NULL: Ⓖ3 asks for > 1.0 bp
        # and Ⓖ4 allows the null up to 8.0, so a run containing ONLY the null satisfied BOTH — the
        # sensitivity control was answered by the very row whose gap is supposed to be zero
        # (`TRAPS: could-the-arm-have-fired`). The candidates are now the NON-null rows, and a run with
        # none of them FAILS rather than reporting an arm that could not have fired.
        gap_arms = [r for r in rows if r not in nulls]
        biggest = max(gap_arms, key=lambda r: abs(r[1])) if gap_arms else None
        _gate("Ⓖ3 a NON-null arm has a realised length gap",
              biggest is not None and abs(biggest[1]) > 1.0,
              f"largest |mu_g-mu_r| = {abs(biggest[1]):.2f} bp on {biggest[0]}" if biggest
              else "NO non-null arm in this run — nothing here can see the defect, only the null")
        # ⛔⛔ AND RESCOPING Ⓖ3 LEFT THE MIRROR HOLE OPEN UNTIL 2026-08-17: Ⓖ4 was guarded by
        # `if nulls:`, so a run carrying ONLY the gap arm SKIPPED it and printed a clean
        # "✅ every gate passed" with no caveat, while the LOAD-BEARING null silently did not run —
        # measured, `--conditions ladder:…capture_on` alone. That is the same defect Ⓖ3 was rescoped to
        # remove, pointing the other way (`TRAPS: could-the-arm-have-fired`). The two arms are
        # SYMMETRIC — the gap arm alone cannot show the defect is not the instrument, exactly as the
        # null alone cannot see the defect — so neither may be skipped, and a missing arm FAILS.
        _gate("Ⓖ4 the NULL (ladder capture-OFF) is equal-length",
              bool(nulls) and max(abs(r[1]) for r in nulls) < 8.0,
              f"|mu_g-mu_r| = {max(abs(r[1]) for r in nulls):.2f} bp; "
              "⚠ capture_on is NOT a null, capture manufactures a gap" if nulls
              else "NO null arm in this run — nothing here can say the defect is not the instrument")

    # ⛔ WHAT DID NOT RUN IS PART OF THE RESULT, and the two gates above are now the whole of that
    # statement: a run missing the null fails Ⓖ4 and a run missing the gap arm fails Ⓖ3, so an
    # incomplete panel can no longer print a clean verdict. ⚠ `skipped` is echoed anyway, because a
    # SKIP line scrolls away above the table and the reader's eye is here.
    if skipped:
        print()
        print(f"  ⚠ {len(skipped)} condition(s) had no cached payload and did not run: "
              f"{', '.join(skipped)}")

    print()
    print("=" * 131)
    if _FAILED:
        print(f"  ⛔ {len(_FAILED)} GATE(S) FAILED — every number above is void until they pass:")
        for f in _FAILED:
            print(f"      {f}")
    else:
        print("  ✅ every gate passed")
    print("  ⛔ This is the NULL plus ONE single-signed capture gap. The published ≤ 0.6 pp bound was")
    print("     measured across two OPPOSITE-SIGN gaps that no panel on disk supplies — do not re-quote")
    print("     it from these rows.")
    print("=" * 131)
    return 1 if _FAILED else 0


if __name__ == "__main__":
    raise SystemExit(main())
