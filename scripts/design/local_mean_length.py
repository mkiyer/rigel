#!/usr/bin/env python
"""⭐⭐⭐ **IS THE MODEL-FREE LOCAL MEAN LENGTH TRUE, AND WHAT DOES IT BUY? — no solver runs.**

At a contiguous EDGE the crossing opportunity is ``w − 1`` and the accumulator's deposit is
``1/(w − 1)``, so the two cancel identically for ANY length distribution::

    E[count]          = rho * E_f[w − 1]
    E[inv_length_sum] = rho * 1
    =>  mhat := count / inv_length_sum + 1 = E_f[w]        <- a LOCAL MEAN LENGTH, no pmf anywhere

⭐ **This instrument asks four questions about that one line, in order, and answers a fifth it did not
promise.**

===  ===================================================================================
 ①   Is the identity TRUE per component, at an EDGE? gDNA is the CONTROL — it cannot splice
 ②   Does it FAIL at a NODE, where ``(ell − w + 1)+`` does not cancel? (the instrument's own
     falsification: an identity that "held" at a node would be measuring something else)
 ③   THE PRECONDITION — natural placement. Capture ON vs OFF, per component; and RNA at an
     EDGE OFF capture, where the only remaining selection is SPLICING
 ④   THE IMPOSSIBILITY SIGNATURE — what share of slots report an ``mhat`` OUTSIDE the two
     fitted component means, where no mixture mean can live
 ⑤   THE INVERSION and its CONDITION NUMBER ``1/(mu_g − mu_r)``
===  ===================================================================================

⛔⛔ **THE SCORING TARGET IS NOT ``length_sum / count``, AND THIS COST A TABLE TO LEARN.** A slot's
``length_sum / count`` is the mean over the fragments that CROSSED the line, i.e. the LANDED mean
``E[w(w−1)]/E[w−1] = E[w] + Var/E[w−1]`` — length-biased by construction, ~+46 bp on this panel. The
identity estimates the NATURAL mean ``E_f[w]``, so scoring ``mhat`` against the landed mean reads a
+20 % "failure" on a channel that is exact. The honest per-component truth is the origin partition's own
``deposited_lengths`` — one increment per ACCEPTED FRAGMENT, beside ``node_start_count`` — and BOTH are
printed here so the reader can see the length bias rather than absorb it.

⛔ Every truth source is gated before it is read (``Ⓖ``): sum-to-full on all three banks this reads —
``_oracle``'s own ``_BANKS`` tuple does NOT cover ``inv_length_sum`` or ``length_sum``, so that identity
had never been checked on the two channels the length model is built from — plus ``deposited_lengths``,
gDNA-cannot-splice, and an ANALYTIC NULL that must be exact to float.

Usage: ``python scripts/design/local_mean_length.py [--drain] [--specs <panel>/<condition> ...]``
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

#: ⭐ The gDNA axis is the axis that matters — one gDNA level cannot validate a composition estimator
#: (`TRAPS: a-single-level-panel-cannot-see-a-constant`), so g00 (truth exactly 0) and g98 are in by
#: default, and both flgap panels supply the LARGE-gap regime the ladder cannot.
_SPECS = tuple(
    f"{panel}/gdna_{lvl}_ss_0.50_nrna_none_capture_{cap}"
    for panel, lvl in (("ladder", "g00"), ("ladder", "g50"), ("ladder", "g98"),
                       ("flgap_short", "g50"), ("flgap_long", "g50"))
    for cap in ("off", "on")
)

#: Count strata. Imported from the census rather than restated — they bracket the recorded low-N bias
#: ladder (0.32 at N=1, 0.05 at N=5, 0.004 at N=50).
from length_channel_census import _BINS, slot_channels  # noqa: E402


def _ell_bins(mu: float):
    """Node-length strata as MULTIPLES OF THE FITTED MEAN FRAGMENT LENGTH — the only length scale the
    tilt ``(ell − w + 1)+`` has. ⭐ Derived, not chosen: below ``mu`` a fragment barely fits, at ``2 mu``
    the tilt is a mild slope, and by ``8 mu`` it is flat and the model reduces to "contained fragments
    have the library distribution".
    """
    return ((0.0, mu), (mu, 2 * mu), (2 * mu, 4 * mu), (4 * mu, 8 * mu), (8 * mu, np.inf))


def _pmf_mean(counts) -> float:
    p = np.asarray(counts, np.float64)
    tot = p.sum()
    return float((np.arange(p.shape[0]) * p).sum() / tot) if tot > 0 else float("nan")


def _landed_mean(counts) -> float:
    """The mean length of fragments CROSSING a line, predicted from a natural pmf: ``E[w(w−1)]/E[w−1]``.

    ⭐ This is what ``length_sum / count`` measures, and it is NOT ``E[w]``. Printed beside the identity
    so the length bias is visible as a number instead of being mistaken for the identity's error.
    """
    p = np.asarray(counts, np.float64)
    w = np.arange(p.shape[0], dtype=np.float64)
    a = np.maximum(w - 1.0, 0.0)
    den = float((p * a).sum())
    return float((p * a * w).sum() / den) if den > 0 else float("nan")


def _local_mean(count_sum, inv_sum):
    """⭐⭐ **THE ESTIMATOR, IN ONE PLACE** — ``count / inv_length_sum + 1``.

    ⛔ Every table AND the analytic null go through this function, so the null gates the FORMULA and not
    merely the data: dropping the ``+ 1``, or inverting the ratio, makes ``Ⓖ6`` fail. The first version
    of this script had the expression written out three times and the null could not see any of them —
    and one of those three copies then divided by a zero density and killed the run before the gate
    block had been printed.

    ⛔ NaN, never an exception, where there is no density to divide by: a degenerate input must be
    REPORTED by the gates, and a gate that raises instead of printing cannot be read.
    """
    num, den = np.asarray(count_sum, np.float64), np.asarray(inv_sum, np.float64)
    out = np.full(np.broadcast(num, den).shape, np.nan)
    np.divide(num, den, out=out, where=den > 0)
    out = out + 1.0
    return out if out.ndim else float(out)


def _ratio(num, den) -> float:
    """A plain safe mean — ⛔ deliberately NOT :func:`_local_mean`, which is the estimator and carries
    the ``+ 1`` the null gates. Two names because they are two quantities."""
    d = float(den)
    return float(num) / d if d > 0 else float("nan")


def _pcts(x, qs=(5, 25, 50, 75, 95)):
    x = np.asarray(x, np.float64)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return (float("nan"),) * len(qs)
    return tuple(float(np.percentile(x, q)) for q in qs)


# ──────────────────────────────────────────────────────────────────────────────────────────────────
#  Ⓖ  THE GATES.  Each one is a claim about a TRUTH SOURCE, and each returns (ok, detail).
# ──────────────────────────────────────────────────────────────────────────────────────────────────
_READ_BANKS = (
    "edge_unspliced_count",
    "edge_unspliced_inv_length_sum",
    "edge_unspliced_length_sum",
    "node_contained_count",
    "node_contained_inv_length_sum",
    "node_contained_length_sum",
)


def analytic_null(pmf, deposit_exponent: int = 1):
    """⭐ THE NULL, and it is ANALYTIC so it cannot be excused. Under natural placement the estimator is
    EXACT for any pmf: ``E[count]/E[inv_length_sum] + 1 = E_f[w]`` identically.

    ``deposit_exponent`` is the falsification handle — the accumulator deposits ``1/(w−1)^1``, and with
    the exponent at 2 the cancellation is destroyed and the null MUST break. A null that cannot be made
    to fail is not a null (`TRAPS: could-the-arm-have-fired`).
    """
    p = np.asarray(pmf, np.float64)
    p = p / p.sum()
    w = np.arange(p.shape[0], dtype=np.float64)
    opp = np.maximum(w - 1.0, 0.0)
    dep = np.zeros_like(w)
    np.divide(1.0, opp**deposit_exponent, out=dep, where=opp > 0)
    e_count = float((p * opp).sum())
    e_inv = float((p * opp * dep).sum())
    # ⭐ through the SAME function every table uses, so the null gates the FORMULA too.
    return _local_mean(e_count, e_inv), float((p * w).sum())


def run_gates(oracle, origins_map, pmf_for_null, drained: bool):
    gates = []

    def raw(payload, bank):
        return np.asarray(getattr(payload, bank), np.int64)

    # ── Ⓖ1 sum-to-full on the THREE banks this instrument reads. ⛔ `_oracle._BANKS` covers the count
    # banks only, so `inv_length_sum` and `length_sum` — the two channels the whole length model is
    # built out of — had never been validated on the origin split at all.
    worst = ("", 0)
    for bank in _READ_BANKS:
        whole = raw(oracle.full, bank)
        parts = sum(raw(oracle.parts[k], bank) for k in ("gdna", "mrna", "nrna"))
        d = int(np.abs(whole.astype(object) - parts.astype(object)).max()) if whole.size else 0
        if d > worst[1]:
            worst = (bank, d)
    gates.append(("Ⓖ1 sum-to-full on count+inv_length_sum+length_sum (EXACT integers)",
                  worst[1] == 0, f"worst |whole − Σparts| = {worst[1]} ({worst[0] or 'all banks'})"))

    # ── Ⓖ2 the per-component NATURAL length truth is itself a partition of the library's
    dl_whole = raw(oracle.full, "deposited_lengths")
    dl_parts = sum(raw(oracle.parts[k], "deposited_lengths") for k in ("gdna", "mrna", "nrna"))
    d2 = int(np.abs(dl_whole - dl_parts).max())
    gates.append(("Ⓖ2 sum-to-full on deposited_lengths (the natural-mean truth source)",
                  d2 == 0, f"worst |whole − Σparts| = {d2}; "
                           f"Σ = {int(dl_whole.sum()):,} accepted fragments"))

    # ── Ⓖ3 gDNA CANNOT SPLICE. If it did, "gDNA is the clean placement control" is void.
    # ⛔ MEASURED BY PERTURBATION: on the UNDRAINED arm this gate is unreachable — `OracleTruth.from_parts`
    # raises on a gDNA spliced deposit before it is ever evaluated, which is a strictly stronger gate. It
    # becomes load-bearing on the DRAINED arm, where `OracleTruth` is constructed directly and the lift
    # is ambiguous BY CONSTRUCTION, so a leak is expected there and is a BOUND rather than a failure —
    # the same reading `flgap_study_cache` records (1 record on flgap_short, 8,641 on flgap_long).
    leak = sum(int(raw(oracle.parts["gdna"], b).sum()) for b in ("edge_spliced_count", "sj_count"))
    g_edge = int(raw(oracle.parts["gdna"], "edge_unspliced_count").sum())
    gates.append(("Ⓖ3 gDNA partition has ZERO spliced incidence (it is the placement control)",
                  leak == 0 or drained,
                  f"edge_spliced_count + sj_count = {leak:,} against {g_edge:,} unspliced gDNA "
                  f"crossings ({100 * _ratio(leak, max(g_edge, 1)):.4f} %)"
                  + ("  ⚠ DRAINED arm: nonzero is the lift's own bound, so the per-component split is a "
                     "CROSS-CHECK, not a verdict" if (leak and drained) else "")))

    # ── Ⓖ4 co-support: a bank with a count must have both length channels, or `mhat` is a 0/0.
    bad = 0
    for axis in ("edge_unspliced", "node_contained"):
        for k in ("gdna", "mrna", "nrna"):
            c = raw(oracle.parts[k], f"{axis}_count").sum(axis=1)
            for ch in ("inv_length_sum", "length_sum"):
                v = raw(oracle.parts[k], f"{axis}_{ch}")
                bad += int(((c > 0) & (v <= 0)).sum())
    gates.append(("Ⓖ4 count > 0 => inv_length_sum > 0 and length_sum > 0 on every bank read",
                  bad == 0, f"{bad:,} violating (object, component) cells"))

    # ── Ⓖ5 COULD THE ARM HAVE FIRED? (`TRAPS: could-the-arm-have-fired`) Every component's accepted
    # fragment total is PRINTED, because a component with zero fragments scores an empty set and its
    # mean is NaN — which is exactly the g00 gDNA arm, and it must be visible rather than inferred.
    lo, hi = 1.0, float(dl_whole.shape[0] - 1)
    detail, ok5, n_live_comp = [], True, 0
    for c, ks in origins_map.items():
        dl = sum(raw(oracle.parts[k], "deposited_lengths") for k in ks)
        n, m = int(dl.sum()), _pmf_mean(dl)
        detail.append(f"{c}: {n:,} fragments, natural E[w] = "
                      f"{'EMPTY ARM — nothing scored' if n == 0 else f'{m:.2f} bp'}")
        if n > 0:
            n_live_comp += 1
            ok5 &= bool(lo <= m <= hi)
    gates.append(("Ⓖ5 could each component's arm have fired, and is its mean in [1, max_length]?",
                  ok5 and n_live_comp > 0, "   ".join(detail) + f"   support ≤ {hi:.0f} bp"))

    # ── Ⓖ6 THE ANALYTIC NULL, both arms.
    got, want = analytic_null(pmf_for_null, 1)
    rel = abs(got - want) / max(want, 1e-12)
    broke, _ = analytic_null(pmf_for_null, 2)
    fired = abs(broke - want) / max(want, 1e-12) > 1e-6
    gates.append(("Ⓖ6 ANALYTIC NULL: exact under natural placement, and it CAN break", rel < 1e-12
                  and fired, f"deposit 1/(w−1): mhat {got:.10f} vs E[w] {want:.10f} (rel {rel:.2e}); "
                             f"deposit 1/(w−1)² breaks it to {broke:.4f} ({'fires' if fired else 'DEAD'})"))
    return gates


# ──────────────────────────────────────────────────────────────────────────────────────────────────
def load_condition(panel: str, condition: str, index, cfg, ra, drain: bool):
    """The cached origin-split oracle for one condition, UNDRAINED by default."""
    from _oracle import ORIGINS, OracleTruth
    from rigel.pipeline import _drain_side_buffer, _native_detect_sj_tag
    from rigel.scan_cache import read_scan_cache

    suite = RUNS / "suite" / panel
    bam = str(suite / condition / "sim_oracle.bam")
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    cache = suite / "oracle_cache" / condition
    payload = read_scan_cache(cache / "_main", index, scan).payload
    parts = {k: read_scan_cache(cache / k, index, scan).payload for k in ORIGINS}
    note = "UNDRAINED (cached pass-1 tally — ⛔ PRODUCTION DRAINS)"
    if drain:
        # ⛔ The parts must be drained by REPLAYING the whole's choices, never independently
        # (`TRAPS: draining-breaks-the-oracle`).
        from rigel.pipeline import scan_and_buffer
        from rigel.second_pass import drain as sp_drain, lift_choices

        _st, sm, _buf, _p = scan_and_buffer(bam, index, scan)
        lift: dict = {}
        payload = _drain_side_buffer(payload, index, sm, seed=cfg.second_pass_seed, _lift=lift)
        lifted, n_amb = lift_choices(lift["undrained"], [parts[k] for k in ORIGINS], lift["choices"])
        parts = {k: sp_drain(parts[k], ch, node_types=lift["node_types"],
                             junctions=lift["junctions"]) for k, ch in zip(ORIGINS, lifted)}
        oracle = OracleTruth(full=payload, parts=parts,
                             read_counts={k: -1 for k in ORIGINS}, n_ambiguous=int(n_amb))
        note = f"DRAINED (production) — lift_ambiguous {int(n_amb):,}"
    else:
        oracle = OracleTruth.from_parts(payload, parts)
    return oracle, payload, scan, note


def measure(panel, condition, index, cfg, ra, drain):
    from length_channel_census import build_channel
    from rigel.calibration.node_chain import EDGE, NODE

    oracle, payload, _scan, note = load_condition(panel, condition, index, cfg, ra, drain)
    chain, geometry, fl, _ll = build_channel(payload, index, ra, cfg)
    kind = np.asarray(chain.kind)
    is_edge, is_node = kind == EDGE, kind == NODE

    origins_map = {"gdna": ("gdna",), "rna": ("mrna", "nrna")}
    gates = run_gates(oracle, origins_map, fl.global_pmf, bool(drain))

    # ── the mixture, exactly as the channel reads it, gated against the parts ──
    mix = slot_channels(oracle, chain, ("gdna", "mrna", "nrna"))
    g_c = np.asarray(geometry.unspliced_count, np.float64).sum(axis=1)
    g_u = np.asarray(geometry.unspliced_inv_length_sum, np.float64)
    dmax = max(float(np.abs(mix[0] - g_c).max()), float(np.abs(mix[1] - g_u).max()))
    gates.append(("Ⓖ7 mixture rebuilt from the parts == the SHIPPED geometry arrays",
                  dmax < 1e-6, f"max |Σparts − geometry| = {dmax:.3e}"))

    # ── Ⓖ8 the components must PARTITION the mixture, and the axis masks must partition the chain.
    # ⛔ Neither Ⓖ1 nor Ⓖ7 sees a component built from the wrong origin set: Ⓖ1 compares raw banks
    # component by component and Ⓖ7 compares the mixture, so double-counting mRNA into `rna` — the
    # nested-aggregate mistake `truth_pmfs` records — would pass both.
    comp_sum = [np.zeros(int(chain.n_slots), np.float64) for _ in range(3)]
    for origins in origins_map.values():
        for i, arr in enumerate(slot_channels(oracle, chain, origins)):
            comp_sum[i] += arr
    d8 = max(float(np.abs(comp_sum[i] - mix[i]).max()) for i in range(3))
    part_ok = bool(not (is_edge & is_node).any() and (is_edge | is_node).all())
    gates.append(("Ⓖ8 Σcomponents == mixture at every slot, and NODE/EDGE partition the chain",
                  d8 < 1e-6 and part_ok,
                  f"max |Σcomp − mix| = {d8:.3e}; masks partition = {part_ok}"))

    # ⛔⛔ PRINTED HERE, THE MOMENT THEY ARE KNOWN, AND NOT BY `report`. Measured by perturbation: with
    # the gate block printed downstream, zeroing `edge_unspliced_inv_length_sum` raised inside table ⑤
    # and the run died with the gate verdict already computed and never shown — a gate that cannot be
    # read is not a gate. Two perturbations were invisible for exactly this reason.
    print()
    print("=" * 112)
    print(f"  {panel}/{condition}")
    print(f"  ARM: {note}")
    print(f"  fitted mu_g = {_pmf_mean(fl.gdna_pmf):.2f} bp   mu_r = {_pmf_mean(fl.rna_pmf):.2f} bp   "
          f"gap = {_pmf_mean(fl.gdna_pmf) - _pmf_mean(fl.rna_pmf):+.2f} bp")
    print("=" * 112)
    print("\n  Ⓖ GATES — every truth source, before it is read")
    for gname, ok, detail in gates:
        print(f"    {'✅' if ok else '⛔ FAIL'}  {gname}")
        print(f"            {detail}")
    sys.stdout.flush()

    node_len = np.zeros(int(chain.n_slots), np.float64)
    nl = np.asarray(ra.region_size_bp, np.float64)
    node_len[is_node] = nl[np.asarray(chain.obj_idx, np.int64)[is_node]]

    mu = {"gdna": _pmf_mean(fl.gdna_pmf), "rna": _pmf_mean(fl.rna_pmf)}

    per = {}
    for comp, origins in origins_map.items():
        c, u, w = slot_channels(oracle, chain, origins)
        dl = sum(np.asarray(oracle.parts[k].deposited_lengths, np.float64) for k in origins)
        nat = _pmf_mean(dl)
        rec = {"nat": nat, "landed_pred": _landed_mean(dl), "n_frag": float(dl.sum()),
               "axes": {}, "bins": {}, "ell": {}}
        for aname, sel in (("EDGE", is_edge), ("NODE", is_node)):
            s = sel & (c > 0) & (u > 0)
            if not s.any():
                continue
            agg_m = _local_mean(c[s].sum(), u[s].sum())
            rec["axes"][aname] = {
                "slots": int(s.sum()), "frags": float(c[s].sum()), "mhat": agg_m,
                "harm": _ratio(c[s].sum(), u[s].sum()),
                "landed": _ratio(w[s].sum(), c[s].sum()),
                "ratio": nat / agg_m if agg_m > 0 else float("nan"),
            }
        s_e = is_edge & (c > 0) & (u > 0)
        for lo, hi in _BINS:
            b = s_e & (c >= lo) & (c <= hi)
            if b.any():
                rec["bins"][(lo, hi)] = (int(b.sum()), float(c[b].sum()),
                                         _pcts(_local_mean(c[b], u[b]) / nat))
        s_n = is_node & (c > 0) & (u > 0)
        for lo, hi in _ell_bins(mu["rna"]):
            b = s_n & (node_len >= lo) & (node_len < hi)
            if b.any():
                m = _local_mean(c[b].sum(), u[b].sum())
                rec["ell"][(lo, hi)] = (int(b.sum()), float(c[b].sum()), m,
                                        _ratio(c[b].sum(), u[b].sum()), nat / m if m > 0 else np.nan)
        per[comp] = rec

    # ── ④ the impossibility signature, on the MIXTURE and the SHIPPED fitted means ──
    mlo, mhi = min(mu["gdna"], mu["rna"]), max(mu["gdna"], mu["rna"])
    c_m, u_m, _w_m = mix
    live = is_edge & (c_m > 0) & (u_m > 0)
    mh = _local_mean(c_m, u_m)
    outside = live & ((mh < mlo) | (mh > mhi))
    tot_true = c_m
    imposs = {
        "live": int(live.sum()), "frac_slots": _ratio(outside.sum(), live.sum()),
        "frac_mass": _ratio(tot_true[outside].sum(), tot_true[live].sum()),
        "interval": (mlo, mhi), "width": mhi - mlo, "pcts": _pcts(mh[live]),
    }

    # ── ⑤ the inversion ──
    ew = _local_mean(c_m[live].sum(), u_m[live].sum())
    gap = mu["gdna"] - mu["rna"]
    phi_hat = (ew - mu["rna"]) / gap if gap != 0 else float("nan")
    # ⭐ THE ESTIMATOR'S OWN TARGET is the DENSITY share, and E[inv_length_sum] = rho exactly, so the
    # truth needs no model either: Σu_gdna / Σu_all over the same slots.
    cg, ug, _ = slot_channels(oracle, chain, ("gdna",))
    cr, ur, _ = slot_channels(oracle, chain, ("mrna", "nrna"))
    ns = {k: float(np.asarray(oracle.parts[k].node_start_count, np.float64).sum())
          for k in ("gdna", "mrna", "nrna")}
    inv = {
        "ew": ew, "mu_g": mu["gdna"], "mu_r": mu["rna"], "gap": gap,
        "cond": abs(1.0 / gap) if gap != 0 else float("inf"), "phi": phi_hat,
        "true_rho": _ratio(ug[live].sum(), ug[live].sum() + ur[live].sum()),
        "true_cnt": _ratio(cg[live].sum(), cg[live].sum() + cr[live].sum()),
        "true_lib": _ratio(ns["gdna"], sum(ns.values())),
    }
    return {"spec": f"{panel}/{condition}", "note": note, "gates": gates,
            "per": per, "imposs": imposs, "inv": inv, "mu": mu}


# ──────────────────────────────────────────────────────────────────────────────────────────────────
def report(r):
    """The measurement tables. ⛔ The gate block is printed by :func:`measure`, not here."""
    print("""
  ①② THE IDENTITY, PER COMPONENT AND PER AXIS.  mhat = count/inv_length_sum + 1
     ⭐ 'natural' is the origin partition's OWN deposited_lengths mean — one increment per ACCEPTED
     FRAGMENT. The 'landed' columns are Σw/N, length-BIASED by construction (E[w(w−1)]/E[w−1]).""")
    print(f"    {'comp':<6}{'axis':<6}{'slots':>9}{'fragments':>13}{'mhat':>9}"
          f"{'natural':>9}{'nat/mhat':>10}{'err %':>8}{'landed':>9}{'land/mhat':>11}"
          f"{'pred':>9}{'l/pred':>8}")
    print("    " + "-" * 113)
    for comp, rec in r["per"].items():
        for axis in ("EDGE", "NODE"):
            a = rec["axes"].get(axis)
            if a is None:
                continue
            print(f"    {comp:<6}{axis:<6}{a['slots']:>9,}{a['frags']:>13,.0f}{a['mhat']:>9.2f}"
                  f"{rec['nat']:>9.2f}{a['ratio']:>10.4f}{100 * (a['ratio'] - 1):>+8.2f}"
                  f"{a['landed']:>9.2f}{a['landed'] / max(a['mhat'], 1e-9):>11.4f}"
                  f"{rec['landed_pred']:>9.2f}"
                  f"{a['landed'] / max(rec['landed_pred'], 1e-9):>8.3f}")
    print("""    ⛔ 'land/mhat' is the ratio the brief asked for and it is NOT 1 even when the identity is
      EXACT — the landed mean is biased by +Var/E[w−1]. Read 'nat/mhat'.
    ⭐ EDGE nat/mhat = 1.000 means the identity holds. ⛔ A NODE row near 1.000 would mean this
      instrument is not measuring the identity: the node deposit is 1/L against an opportunity
      (ell−w+1)+, and nothing cancels.

  ① PER-SLOT DISTRIBUTION at EDGE — mhat / (that component's natural E[w]), by count N""")
    print(f"    {'comp':<6}{'N':<9}{'slots':>9}{'fragments':>14}{'p05':>8}{'p25':>8}{'p50':>8}"
          f"{'p75':>8}{'p95':>8}")
    print("    " + "-" * 78)
    for comp, rec in r["per"].items():
        for (lo, hi), (n, f, p) in rec["bins"].items():
            label = f"{lo}" if lo == hi else (f"{lo}+" if hi > 10**8 else f"{lo}-{hi}")
            print(f"    {comp:<6}{label:<9}{n:>9,}{f:>14,.0f}" + "".join(f"{v:>8.3f}" for v in p))

    print("\n  ② THE NODE CONTRAST, by node length ell (bins = multiples of the fitted RNA mean)")
    print(f"    {'comp':<6}{'ell bp':<14}{'slots':>9}{'fragments':>14}{'mhat':>9}{'harm mean':>11}"
          f"{'natural':>9}{'nat/mhat':>10}")
    print("    " + "-" * 84)
    for comp, rec in r["per"].items():
        for (lo, hi), (n, f, m, h, ratio) in rec["ell"].items():
            label = f"{lo:.0f}+" if not np.isfinite(hi) else f"{lo:.0f}-{hi:.0f}"
            print(f"    {comp:<6}{label:<14}{n:>9,}{f:>14,.0f}{m:>9.2f}{h:>11.2f}"
                  f"{rec['nat']:>9.2f}{ratio:>10.4f}")

    im, v = r["imposs"], r["inv"]
    print(f"""
  ④ THE IMPOSSIBILITY SIGNATURE — a mixture mean must lie between the component means
     interval [{im['interval'][0]:.2f}, {im['interval'][1]:.2f}] bp, width {im['width']:.2f} bp\
     live EDGE slots {im['live']:,}
     OUTSIDE:  {100 * im['frac_slots']:.2f} % of live slots   \
{100 * im['frac_mass']:.2f} % of crossing fragments
     per-slot mhat percentiles: p05 {im['pcts'][0]:.1f}  p25 {im['pcts'][1]:.1f}  \
p50 {im['pcts'][2]:.1f}  p75 {im['pcts'][3]:.1f}  p95 {im['pcts'][4]:.1f}

  ⑤ THE INVERSION — phi = (E[w] − mu_r)/(mu_g − mu_r), and its CONDITION NUMBER
     mass-weighted E[w] over live EDGE slots = {v['ew']:.3f} bp
     phi_hat = {v['phi']:+.4f}    |1/(mu_g − mu_r)| = {v['cond']:.4f} per bp
     TRUE density share Σu_g/Σu (the estimator's own target) = {v['true_rho']:.4f}
     TRUE edge unspliced COUNT share                        = {v['true_cnt']:.4f}
     TRUE library f_gdna (node_start_count)                 = {v['true_lib']:.4f}""")


def summary(rows):
    print("\n" + "=" * 112)
    print("""  ③ THE PRECONDITION, MEASURED — natural placement, capture OFF vs ON, per component
     err % = 100 * (natural E[w] / mhat − 1) at EDGE. gDNA cannot splice, so off capture it is the
     clean control; under capture the weight is a joint (start, length) function.""")
    print("=" * 112)
    print(f"    {'panel/level':<26}{'comp':<6}{'axis':<6}{'OFF err %':>11}{'ON err %':>11}"
          f"{'ON − OFF':>11}")
    print("    " + "-" * 72)
    by = {}
    for r in rows:
        panel, cond = r["spec"].split("/")
        key = f"{panel}/{cond.split('_')[1]}"
        by.setdefault(key, {})["on" if cond.endswith("_on") else "off"] = r
    for key, pair in by.items():
        for comp in ("gdna", "rna"):
            for axis in ("EDGE", "NODE"):
                e = {}
                for arm in ("off", "on"):
                    a = pair.get(arm, {}).get("per", {}).get(comp, {}).get("axes", {}).get(axis)
                    e[arm] = 100 * (a["ratio"] - 1) if a else float("nan")
                print(f"    {key:<26}{comp:<6}{axis:<6}{e['off']:>+11.2f}{e['on']:>+11.2f}"
                      f"{e['on'] - e['off']:>+11.2f}")

    # ⭐ THE SHARPER PRECONDITION QUESTION: at a LINE the bank holds only the UNSPLICED crossings, so
    # mature RNA's placement is ALREADY non-natural before any capture selection exists. OFF capture
    # there is no other selection left, so the residual there IS splicing — and it is printed beside
    # EQUATIONS 3e's recorded 7.5 % line deficit, in 3e's OWN statistic (the landed second moment),
    # so the two can be compared as numbers rather than as stories.
    print("\n  ③b RNA AT AN EDGE, OFF CAPTURE — the only selection left is SPLICING")
    print(f"    {'panel/level':<22}{'mhat deficit %':>16}{'landed m2 deficit %':>22}"
          f"{'EQUATIONS 3e':>15}")
    print("    " + "-" * 76)
    for key, pair in by.items():
        r = pair.get("off")
        a = (r or {}).get("per", {}).get("rna", {}).get("axes", {}).get("EDGE")
        if a is None:
            continue
        rec = r["per"]["rna"]
        print(f"    {key:<22}{100 * (a['mhat'] / rec['nat'] - 1):>+16.2f}"
              f"{100 * (a['landed'] / rec['landed_pred'] - 1):>+22.2f}{'−7.5 %':>15}")
    print("""    ⭐ Both are NEGATIVE and both are ONE mechanism: a spliced fragment's crossings go to
      `edge_spliced`, so `edge_unspliced` keeps only paths contained in one exonic block — a joint
      length x placement selection that truncates the long tail. The 'landed m2' column is 3e's own
      statistic, so if it lands near −7.5 % the two are the SAME number seen twice.""")
    print("\n" + "=" * 112)
    print("  ④⑤ THE BRIDGE TO STUDY 2 — the interval, who is outside it, and the condition number")
    print("=" * 112)
    print(f"    {'spec':<52}{'gap bp':>9}{'1/gap':>9}{'out slot%':>11}{'out mass%':>11}"
          f"{'phi_hat':>9}{'true rho':>9}")
    print("    " + "-" * 110)
    for r in rows:
        im, v = r["imposs"], r["inv"]
        print(f"    {r['spec']:<52}{v['gap']:>+9.2f}{v['cond']:>9.4f}"
              f"{100 * im['frac_slots']:>11.2f}{100 * im['frac_mass']:>11.2f}"
              f"{v['phi']:>+9.3f}{v['true_rho']:>9.3f}")

    print("""
  ⑥ THE HONEST NEGATIVE, from the numbers above:
     • The identity measures E[w] with NO pmf, and at an EDGE it is exact to a fraction of a percent
       on the gDNA control — that part is real and needs no length model.
     • It does NOT measure composition. Turning E[w] into phi needs BOTH global component means, so
       the model-freeness ends at the inversion.
     • The sensitivity is 1/(mu_g − mu_r) — the same 2x2 determinant every other channel in the
       deconvolution depends on. Read '1/gap' beside 'out mass%': where the gap is small the
       estimator is exact and USELESS, not biased.""")


def main() -> int:
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.config import PipelineConfig
    from rigel.index import TranscriptIndex

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--specs", nargs="*", default=list(_SPECS), help="panel/condition")
    ap.add_argument("--drain", action="store_true",
                    help="DRAIN the payloads (production). Slow: it re-scans the BAM.")
    args = ap.parse_args()

    index = TranscriptIndex.load(str(INDEX))
    cfg = PipelineConfig()
    ra = RegionArrays.from_index(index)

    print("\n" + "#" * 112)
    print("#  ⭐⭐⭐ THE MODEL-FREE LOCAL MEAN LENGTH — mhat = count/inv_length_sum + 1\n"
          f"#  ARM: {'DRAINED (production)' if args.drain else 'UNDRAINED (cached pass-1 tally)'}\n"
          "#  No solver runs. Every truth is the origin-split oracle cache.")
    print("#" * 112)

    rows, failed = [], 0
    for spec in args.specs:
        panel, condition = spec.split("/", 1)
        r = measure(panel, condition, index, cfg, ra, args.drain)
        report(r)
        failed += sum(1 for _n, ok, _d in r["gates"] if not ok)
        rows.append(r)
    if len(rows) > 1:
        summary(rows)
    print(f"\n  GATES FAILED: {failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
