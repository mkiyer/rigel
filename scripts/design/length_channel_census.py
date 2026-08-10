#!/usr/bin/env python
"""⭐⭐⭐ **DOES THE FRAGMENT-LENGTH CHANNEL KNOW ANYTHING TRUE? — no solver runs.**

**Why this exists, and why it must run before any A/B.** ``length_likelihood`` has never been switched
on. It is a written idea, not a measured feature: nothing has established that it fires, that what it
claims is true, or that the precision it declares is earned. Turning it on and reading a panel number
first would produce a figure attributable to the channel, the relay, the prior and the solver at once —
and an untrusted channel is the worst of those to have to infer.

⭐ **So this isolates it completely.** It builds the SHIPPED moments off the SHIPPED geometry, evaluates
the SHIPPED ``length_loglik`` on the real per-slot ``(count, inv_length_sum, length_sum)``, and scores
its argmax against the origin-split oracle's TRUE ``f_g`` at the same slot. No sweep, no messages, no
prior. Whatever this says is the channel's own opinion and nobody else's.

⛔⛔ **THE RISK IT IS AIMED AT IS ALREADY ON THE RECORD AND IS NOT SMALL.** The likelihood is a bivariate
Gaussian on ``(Σ 1/w, Σ w)`` — two statistics of the SAME draw, deterministically related — summed over
``N`` fragments, so it is asymptotic in ``N``. Measured in
``test_the_log_det_term_is_present_and_dominates_at_low_count``, the heteroscedastic ``−½ log det`` term
displaces the peak by **0.32 at N = 1**, 0.05 at 5, 0.004 at 50. ⚠ And the node axis got THINNER since
that was written: the three ``node_spanning`` banks were deleted 2026-08-08, so a NODE slot's evidence
is now the contained population alone. ⛔ ``tau_len`` is deliberately UNGATED (unlike ``i_strand``), so a
biased low-``N`` mode enters the solve carrying a precision — the shape that outvotes correct
neighbours.

⭐ **THE THREE QUESTIONS, IN ORDER. A "no" to an earlier one makes the later ones unreadable.**

===  ==================================================================  ==========================
 ①   Does it FIRE, and on how much of the library?                       live slots, live mass
 ②   Is the covariance CONDITIONED well enough to trust the location?    1 − rho², by slot type
 ③   Is what it claims TRUE, and is its precision EARNED?                |argmax − truth|, by N
===  ==================================================================  ==========================

⛔ **③ IS STRATIFIED BY ``N`` AND THAT IS THE POINT, NOT A FORMALITY.** A mass-weighted average over all
slots is dominated by the few deep ones, where the Gaussian is valid and the answer is unsurprising. The
question is what happens at the slots the genome is actually made of.

⭐ **TWO NULLS, BOTH REQUIRED** (`TRAPS: could-the-arm-have-fired`): feeding ONE pmf as both components
must make every row exactly inert — if it does not, the channel is reading something other than the
length gap — and the reported live fraction must be non-zero, or ③ is scoring an empty set.

Usage::

    python scripts/design/length_channel_census.py --panel flgap_short --condition ...capture_on
    python scripts/design/length_channel_census.py --suite ladder --condition gdna_g00_ss_0.50_...
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

#: The count bins ③ is stratified by. ⭐ Chosen to bracket the measured bias ladder (0.32 at N=1,
#: 0.05 at N=5, 0.004 at N=50), not by round numbers — so a row's bias is predictable from theory and
#: the table either confirms or refutes it.
_BINS = ((1, 1), (2, 4), (5, 19), (20, 49), (50, 10**9))


def build_channel(payload, index, region_arrays, config):
    """The SHIPPED construction, in the SHIPPED order — chain, geometry, moments, log-likelihood.

    ⛔ Every object here is built by the same call ``calibrate`` makes, in the same order, from the same
    inputs. A census that rebuilt the geometry its own way would be measuring a channel the tool does
    not ship (`TRAPS: a-test-that-redefines`).
    """
    from rigel.calibration.calibrate import _build_length_loglik
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.junction_opportunity import crossing_probability_from_index
    from rigel.calibration.node_chain import build_node_chain
    from rigel.calibration.node_geometry import build_node_geometry
    from rigel.calibration.splice_graph import build_junction_geometry_arrays
    from rigel.calibration.substrate import CalibrationSubstrate

    max_size = int(payload.max_length)
    fl = build_fl_models(
        payload,
        junction_opportunity=crossing_probability_from_index(index, max_size),
        gdna_opportunity=gdna_opportunity_from_index(index, max_size),
    )
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)
    chain = build_node_chain(payload.ref_node_offsets, payload.ref_edge_offsets)
    geometry = build_node_geometry(
        chain, substrate, region_arrays, build_junction_geometry_arrays(index),
        fl.gdna_pmf, fl.rna_pmf, None,
    )
    on = dataclasses.replace(config.calibration, length_likelihood=True)
    ll = _build_length_loglik(chain, geometry, region_arrays, fl.gdna_pmf, fl.rna_pmf, on)
    return chain, geometry, fl, ll





def truth_pmfs(suite, condition, max_size):
    """The condition's OWN realised per-component length histogram, from the simulator's truth table.

    ⭐ This is what separates the two possible repairs, and they are different work:

    * predicted-from-FITTED vs predicted-from-TRUE  →  the pmf FIT's error
    * predicted-from-TRUE   vs realised             →  a missing SELECTION term in the channel's model
      of its own population, which no amount of better fitting would remove

    ⚠ ``truth_fragment_lengths.tsv`` is the GENERATED library, not the accepted one — fragments dropped
    as too long, unmapped or chimeric are in it and not in the banks. Small, and stated rather than
    assumed.
    """
    import csv

    # ⛔⛔ THE `kind` COLUMN CARRIES NESTED AGGREGATES, NOT DISJOINT POPULATIONS: `mrna`, `nrna`, `gdna`
    # are the populations; `rna` is (mrna+nrna) again and `all` is the whole library again. The first
    # version of this function did `"gdna" if kind == "gdna" else "rna"`, which summed mrna + rna + all
    # into one bucket — mRNA counted twice plus the library once — and reported a "true RNA library" of
    # 191.06 bp where mRNA is 212.20. gDNA was unaffected, so the gDNA control read 1.000 and the RNA
    # arm looked broken: a self-inflicted "the control passes, the treatment fails" that survived two
    # rounds of interpretation. ⭐ Hence the explicit membership sets and the raise: an unrecognised
    # kind must stop the run, not pick a bucket.
    populations = {"mrna": "rna", "nrna": "rna", "gdna": "gdna"}
    aggregates = {"rna", "all", "total"}

    path = Path(suite) / condition / "truth_fragment_lengths.tsv"
    acc = {"gdna": np.zeros(max_size + 1), "rna": np.zeros(max_size + 1)}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            kind = row["kind"]
            if kind in aggregates:
                continue
            if kind not in populations:
                raise RuntimeError(
                    f"{path}: unrecognised kind {kind!r}. It is either a new population (add it) or a "
                    "new AGGREGATE (skip it) — and guessing is how mRNA came to be counted twice."
                )
            w = int(row["fragment_length"])
            if 0 <= w <= max_size:
                acc[populations[kind]][w] += float(row["count"])
    return {k: (v / v.sum() if v.sum() > 0 else v) for k, v in acc.items()}


def slot_channels(oracle, chain, origins):
    """Per-slot ``(count, sum_u, sum_w)`` for ONE component, off the origin split.

    ⭐ The EXACT three numbers ``length_loglik`` conditions on, but for a single component instead of the
    mixture — so the moments the channel ASSUMES can be compared with the moments the population
    actually has. ⛔ ``inv_length_sum`` is descaled here by the same ``INV_LENGTH_SCALE`` the shipped
    ``substrate.PopulationView`` applies; the raw payload is fixed point.
    """
    from rigel.calibration.node_chain import EDGE, NODE
    from rigel.calibration.substrate import INV_LENGTH_SCALE

    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    n = int(chain.n_slots)
    is_node, is_edge = kind == NODE, kind == EDGE

    def bank(name, axis2=False):
        a = sum(np.asarray(getattr(oracle.parts[k], name), np.float64) for k in origins)
        return a.sum(axis=1) if axis2 else a

    out = []
    for node_bank, edge_bank, scale, two in (
        ("node_contained_count", "edge_unspliced_count", 1.0, True),
        ("node_contained_inv_length_sum", "edge_unspliced_inv_length_sum", INV_LENGTH_SCALE, False),
        ("node_contained_length_sum", "edge_unspliced_length_sum", 1.0, False),
    ):
        v = np.zeros(n, np.float64)
        v[is_node] = bank(node_bank, two)[obj[is_node]] / scale
        v[is_edge] = bank(edge_bank, two)[obj[is_edge]] / scale
        out.append(v)
    return tuple(out)



def exonic_block_distances(chain, region_arrays, pure: bool = False):
    """Per EDGE slot, the distance to each end of the contiguous EXONIC block containing the line.

    ⭐ An unspliced mature-RNA fragment crossing a line cannot leave the exon it is in — leaving means
    using a junction, which puts it in `edge_spliced`, a different bank. So the two distances that
    bound it are the ends of the contiguous run of EXON regions around the line. Within an exon,
    transcript-space distance IS genomic distance (no introns inside), so this is a plain chain walk.

    ⛔ Zero on the side where the adjacent region is not exonic: the line is then already a junction (or
    an exon/intergenic boundary) and an unspliced RNA crossing would have to enter an intron, which is
    nascent, not mature.
    """
    from rigel.calibration.node_chain import EDGE, NODE
    from rigel.calibration.signature import coarse_type_array

    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    n = int(chain.n_slots)
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    rlen = np.asarray(region_arrays.region_size_bp, np.float64)

    is_node = kind == NODE
    node_slots = np.flatnonzero(is_node)
    if pure:
        # ⭐ PURELY exonic: the exon bit set AND the intron bit clear. `coarse_type_array` lets exon
        # WIN over intron, so a region that is exonic for transcript B and intronic for transcript A
        # counts as exon and the block runs straight through a junction that is real for A. That
        # over-extends d_left/d_right and over-predicts the mean — which is the sign of the residual.
        from rigel.calibration.signature import (
            BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS,
        )

        sig = np.asarray(region_arrays.signature).astype(np.int64)[obj[node_slots]]
        node_is_exon = ((sig & (BIT_EXON_POS | BIT_EXON_NEG)) != 0) & (
            (sig & (BIT_INTRON_POS | BIT_INTRON_NEG)) == 0
        )
    else:
        node_is_exon = rtype[obj[node_slots]] == 2  # RegionType.EXON
    node_len = rlen[obj[node_slots]]

    # cumulative exonic length running LEFT from each node, and RIGHT — reset at every non-exon node.
    # ⚠ A reference boundary shows up as a break in the chain's node/edge alternation, and a non-exon
    # node (intergenic) always sits at one, so the reset handles it.
    m = node_slots.size
    left = np.zeros(m)
    for i in range(m):
        left[i] = (left[i - 1] + node_len[i]) if (i > 0 and node_is_exon[i]) else (
            node_len[i] if node_is_exon[i] else 0.0
        )
    right = np.zeros(m)
    for i in range(m - 1, -1, -1):
        right[i] = (right[i + 1] + node_len[i]) if (i < m - 1 and node_is_exon[i]) else (
            node_len[i] if node_is_exon[i] else 0.0
        )

    d_left = np.zeros(n)
    d_right = np.zeros(n)
    pos = np.full(n, -1, np.int64)
    pos[node_slots] = np.arange(m)
    edge_slots = np.flatnonzero(kind == EDGE)
    li, ri = edge_slots - 1, edge_slots + 1
    ok = (li >= 0) & (ri < n) & (pos[np.clip(li, 0, n - 1)] >= 0) & (pos[np.clip(ri, 0, n - 1)] >= 0)
    e_ok = edge_slots[ok]
    pl, pr = pos[e_ok - 1], pos[e_ok + 1]
    d_left[e_ok] = np.where(node_is_exon[pl], left[pl], 0.0)
    d_right[e_ok] = np.where(node_is_exon[pr], right[pr], 0.0)
    return d_left, d_right


def bounded_crossing_moments(fl_pmf, d_left, d_right):
    """Crossing moments with the UNSPLICED opportunity ``A(w) = #{u : max(1,w−dr) <= u <= min(w−1,dl)}``.

    ⭐ One expression that does both jobs the measurement asked for: it drops the placements whose
    fragment would have run past a junction (and so landed in `edge_spliced`, not this bank), and it
    shrinks the opportunity where the exon is short. ⛔ ``w > d_left + d_right`` gives exactly 0 — a line
    inside an exon shorter than the fragment admits NO unspliced crossing, which is the owner's
    tiny-exon observation as a limiting case.

    ⚠ The deposit weight ``u(w) = 1/(w−1)`` is UNCHANGED — that is what the accumulator deposits, and
    this function may only change the OPPORTUNITY (`TRAPS: a-test-that-redefines`).
    """
    from rigel.calibration.length_likelihood import LandedMoments

    p = np.asarray(fl_pmf, np.float64)
    p = p / p.sum() if p.sum() > 0 else p
    w = np.arange(p.shape[0], dtype=np.float64)
    inv = np.zeros_like(w)
    np.divide(1.0, w - 1.0, out=inv, where=w >= 2.0)

    dl = np.asarray(d_left, np.float64)[:, None]
    dr = np.asarray(d_right, np.float64)[:, None]
    upper = np.minimum(w[None, :] - 1.0, dl)
    lower = np.maximum(1.0, w[None, :] - dr)
    A = np.maximum(0.0, upper - lower + 1.0)

    pw = p[None, :]
    eff = (A * pw).sum(axis=1)
    live = eff > 0.0

    def d(x):
        return np.divide(x, eff, out=np.zeros_like(eff), where=live)

    return LandedMoments(
        m1=d((A * pw * inv[None, :]).sum(axis=1)),
        m2=d((A * pw * w[None, :]).sum(axis=1)),
        q1=d((A * pw * inv[None, :] ** 2).sum(axis=1)),
        q2=d((A * pw * w[None, :] ** 2).sum(axis=1)),
        q12=d((A * pw * inv[None, :] * w[None, :]).sum(axis=1)),
        eff=eff,
    )


def fisher_information(m_g, m_r, count, fg_grid):
    """⭐⭐ **THE CHANNEL'S OWN FISHER INFORMATION ABOUT ``lambda``, ANALYTIC — the derived alternative
    to a discrimination GATE.**

    The likelihood is Gaussian with a mean that moves with ``pi`` and a covariance that also does. The
    mean-driven information is the whole of the asymptotically trustworthy part::

        Delta   = (m1_g - m1_r, m2_g - m2_r)        the two components' per-fragment moment gap
        V(pi)   = per-draw mixture covariance
        I_pi    = N * Delta' V(pi)^-1 Delta         <- the squared MAHALANOBIS distance between the
                                                       two components' length signatures, times N
        I_lam   = I_pi * [pi(1-pi)]^2               push forward to the log-odds axis

    ⭐ **Why this replaces the gate rather than tightening it.** ``I -> 0`` smoothly and QUADRATICALLY as
    the two pmfs converge, and ``Delta == 0`` exactly is just its limiting case — so the exact-inequality
    gate becomes redundant instead of needing a threshold nobody can derive. It also scales linearly with
    ``N``, and it never touches the solve grid, so it carries no ``1/Var(grid)`` floor.

    ⛔ It deliberately OMITS the ``1/2 tr[(V^-1 dV/dpi)^2]`` term. That is the heteroscedastic
    contribution — the same ``-1/2 log det`` that displaces the peak by 0.32 at ``N = 1`` — and it is
    ``O(1)`` in ``N``, i.e. it does not vanish at low depth where it is least trustworthy. Reporting the
    mean-driven part alone is the conservative reading.

    Returns ``(n_slots,)`` — the information at each slot's own maximum over the grid.
    """
    n = np.asarray(count, np.float64)[:, None]
    fg = np.asarray(fg_grid, np.float64)[None, :]

    def mix(a, b):
        return fg * np.asarray(a, np.float64)[:, None] + (1.0 - fg) * np.asarray(b, np.float64)[:, None]

    mu_d, mu_s = mix(m_g.m1, m_r.m1), mix(m_g.m2, m_r.m2)
    v_dd = mix(m_g.q1, m_r.q1) - mu_d * mu_d          # PER-DRAW, not times N
    v_ss = mix(m_g.q2, m_r.q2) - mu_s * mu_s
    v_ds = mix(m_g.q12, m_r.q12) - mu_d * mu_s
    det = v_dd * v_ss - v_ds * v_ds

    d_d = (np.asarray(m_g.m1, np.float64) - np.asarray(m_r.m1, np.float64))[:, None]
    d_s = (np.asarray(m_g.m2, np.float64) - np.asarray(m_r.m2, np.float64))[:, None]
    quad = v_ss * d_d * d_d - 2.0 * v_ds * d_d * d_s + v_dd * d_s * d_s
    ok = (det > 0.0) & (n > 0.0)
    i_pi = np.zeros_like(quad)
    np.divide(n * quad, det, out=i_pi, where=ok)
    i_lam = i_pi * (fg * (1.0 - fg)) ** 2
    return i_lam.max(axis=1)


def slot_truth(oracle, chain, region_arrays):
    """TRUE ``f_g`` per chain slot, on the SAME population the channel conditions on.

    ⛔ The channel reads ``node_contained`` at a NODE and ``edge_unspliced`` at an EDGE, so the truth
    must be the gDNA share of exactly those two banks — not of the locus, not of the mass. A truth on a
    different population is the error this campaign already made twice
    (`TRAPS: score-the-consumers-own-count`).
    """
    from rigel.calibration.node_chain import EDGE, NODE

    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    n = int(chain.n_slots)

    def bank(name, origins):
        return sum(
            np.asarray(getattr(oracle.parts[k], name), np.float64).sum(axis=1) for k in origins
        )

    node_g = bank("node_contained_count", ("gdna",))
    node_r = bank("node_contained_count", ("mrna", "nrna"))
    edge_g = bank("edge_unspliced_count", ("gdna",))
    edge_r = bank("edge_unspliced_count", ("mrna", "nrna"))

    g = np.zeros(n, np.float64)
    r = np.zeros(n, np.float64)
    is_node, is_edge = kind == NODE, kind == EDGE
    g[is_node], r[is_node] = node_g[obj[is_node]], node_r[obj[is_node]]
    g[is_edge], r[is_edge] = edge_g[obj[is_edge]], edge_r[obj[is_edge]]
    tot = g + r
    fg = np.full(n, np.nan)
    np.divide(g, tot, out=fg, where=tot > 0)
    return fg, tot


def main() -> int:
    from rigel.calibration.node_chain import EDGE, NODE
    from rigel.calibration.density_deconv import density_factor_precision
    from rigel.calibration.simplex_logodds import _logodds_grid
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.config import PipelineConfig
    from rigel.index import TranscriptIndex
    from rigel.pipeline import _drain_side_buffer, _native_detect_sj_tag

    from _oracle import ORIGINS, OracleTruth

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--panel", default="flgap_short")
    ap.add_argument("--condition", default="gdna_g50_ss_0.50_nrna_none_capture_on")
    ap.add_argument("--undrained", action="store_true",
                    help="read the cached UNDRAINED payload. ⛔ Production DRAINS, so this is a "
                         "diagnostic arm and not the tool's own substrate.")
    args = ap.parse_args()

    suite = RUNS / "suite" / args.panel
    bam = str(suite / args.condition / "sim_oracle.bam")
    index = TranscriptIndex.load(str(INDEX))
    cfg = PipelineConfig()
    ra = RegionArrays.from_index(index)
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))

    from rigel.scan_cache import read_scan_cache

    # ⛔⛔ PRODUCTION DRAINS BEFORE CALIBRATING, AND THE DRAIN IS NOT SMALL ON THIS CHANNEL.
    # The first version of this census read the cached (undrained) payload and understated the RNA
    # pmf's error by 8 bp: the side buffer holds fragments whose junction fell in the unsequenced mate
    # gap, which requires a LONG fragment, so draining adds tail mass to the very pool `rna_pmf` is
    # fitted on (217.93 -> 226.10 here). `TRAPS: draining-breaks-the-oracle` is why the parts must be
    # drained by REPLAYING the whole's choices rather than independently.
    cache = suite / "oracle_cache" / args.condition
    payload = read_scan_cache(cache / "_main", index, scan).payload
    parts = {k: read_scan_cache(cache / k, index, scan).payload for k in ORIGINS}
    drain_note = "UNDRAINED (diagnostic arm — production drains)"
    if not args.undrained:
        from rigel.pipeline import scan_and_buffer
        from rigel.second_pass import drain as sp_drain, lift_choices

        _st, sm, _buf, _p = scan_and_buffer(bam, index, scan)
        lift: dict = {}
        payload_d = _drain_side_buffer(payload, index, sm, seed=cfg.second_pass_seed, _lift=lift)
        lifted, n_amb = lift_choices(lift["undrained"], [parts[k] for k in ORIGINS], lift["choices"])
        parts = {
            k: sp_drain(parts[k], ch, node_types=lift["node_types"], junctions=lift["junctions"])
            for k, ch in zip(ORIGINS, lifted)
        }
        leak = sum(
            int(np.asarray(getattr(parts["gdna"], b), np.int64).sum())
            for b in ("edge_spliced_count", "sj_count")
        )
        payload = payload_d
        drain_note = (f"DRAINED (production) — lift_ambiguous {int(n_amb):,}, "
                      f"gdna spliced leak {leak:,} ⛔ gDNA cannot splice; a large leak makes the "
                      f"per-component split a cross-check, not a verdict")
        oracle = OracleTruth(full=payload, parts=parts,
                             read_counts={k: -1 for k in ORIGINS}, n_ambiguous=int(n_amb))
    else:
        oracle = OracleTruth.from_parts(payload, parts)

    chain, geometry, fl, ll = build_channel(payload, index, ra, cfg)
    kind = np.asarray(chain.kind)
    count = np.asarray(geometry.unspliced_count, np.float64).sum(axis=1)
    fg_true, tot_true = slot_truth(oracle, chain, ra)
    lam_grid, fg_grid = _logodds_grid(int(cfg.calibration.sweep_n_grid),
                                      float(cfg.calibration.sweep_logodds_window))

    live = np.ptp(ll, axis=1) > 0.0
    argmax = fg_grid[np.argmax(ll, axis=1)]
    tau = density_factor_precision(ll, lam_grid)

    # ── the REPAIRED channel: the RNA crossing opportunity bounded by the exonic block ──
    # ⭐ ASYMMETRIC BY PHYSICS, and that is the whole point: gDNA cannot splice, so a gDNA fragment
    # crossing a line never leaves via a junction and its opportunity stays `w − 1`. Only the RNA arm's
    # EDGE moments change; both NODE frames are untouched (containment already excludes splicing).
    from rigel.calibration.length_likelihood import build_slot_moments as _bsm, length_loglik as _ll

    _dl, _dr = exonic_block_distances(chain, ra)
    _bound = bounded_crossing_moments(fl.rna_pmf, _dl, _dr)
    _mg2, _mr2 = _bsm(chain, ra, fl.gdna_pmf), _bsm(chain, ra, fl.rna_pmf)
    _is_edge = kind == EDGE
    fixed_r = type(_mr2)(**{
        f: np.where(_is_edge, np.asarray(getattr(_bound, f), np.float64),
                    np.asarray(getattr(_mr2, f), np.float64))
        for f in ("m1", "m2", "q1", "q2", "q12", "eff")
    })
    ll_fix = _ll(_mg2, fixed_r, count,
                 np.asarray(geometry.unspliced_inv_length_sum, np.float64),
                 np.asarray(geometry.unspliced_length_sum, np.float64), fg_grid)
    live_fix = np.ptp(ll_fix, axis=1) > 0.0
    argmax_fix = fg_grid[np.argmax(ll_fix, axis=1)]
    tau_fix = density_factor_precision(ll_fix, lam_grid)

    print()
    print("=" * 104)
    print(f"  ⭐⭐⭐ THE LENGTH CHANNEL, ALONE — {args.panel}/{args.condition}")
    print(f"  {drain_note}")
    print(f"  mu_g = {float((np.arange(fl.gdna_pmf.shape[0]) * fl.gdna_pmf).sum()):.1f} bp   "
          f"mu_r = {float((np.arange(fl.rna_pmf.shape[0]) * fl.rna_pmf).sum()):.1f} bp   "
          f"grid K = {fg_grid.shape[0]}")
    print("=" * 104)

    # ── the NULL, first. A channel that speaks with ONE pmf on both arms is reading noise. ──
    from rigel.calibration.length_likelihood import build_slot_moments, length_loglik
    m_same = build_slot_moments(chain, ra, fl.gdna_pmf)
    null = length_loglik(m_same, m_same, count,
                         np.asarray(geometry.unspliced_inv_length_sum, np.float64),
                         np.asarray(geometry.unspliced_length_sum, np.float64), fg_grid)
    print(f"\n  ⭐ NULL (one pmf on both arms): max |ptp| over every slot = {float(np.abs(np.ptp(null, axis=1)).max()):.3e}"
          f"   {'✅ exactly inert' if float(np.abs(np.ptp(null, axis=1)).max()) == 0.0 else '⛔ SPEAKS WITH NO GAP'}")

    # ── ① does it fire, and on how much? ──
    print("\n  ① DOES IT FIRE?  (a slot is live iff its row is not flat)")
    print(f"    {'slot type':<14} {'slots':>12} {'live':>12} {'live %':>8} "
          f"{'fragments':>14} {'live frag':>14} {'live frag %':>12}")
    print("    " + "-" * 92)
    for name, sel in (("NODE", kind == NODE), ("EDGE", kind == EDGE), ("all", np.ones_like(kind, bool))):
        f_tot = float(tot_true[sel].sum())
        f_live = float(tot_true[sel & live].sum())
        print(f"    {name:<14} {int(sel.sum()):>12,} {int((sel & live).sum()):>12,} "
              f"{100 * (sel & live).sum() / max(sel.sum(), 1):>7.1f}% {f_tot:>14,.0f} "
              f"{f_live:>14,.0f} {100 * f_live / max(f_tot, 1):>11.1f}%")

    # ── ② the conditioning ──
    print("\n  ② IS THE COVARIANCE CONDITIONED?  u = 1/w and w are the SAME draw, so rho -> -1 by")
    print("     construction and det = v_dd*v_ss*(1-rho²) is a cancellation. 1-rho² is the headroom.")
    m_g = build_slot_moments(chain, ra, fl.gdna_pmf)
    m_r = build_slot_moments(chain, ra, fl.rna_pmf)
    print(f"    {'slot type':<14} {'min 1-rho²':>14} {'median 1-rho²':>16} {'p01':>14}")
    print("    " + "-" * 62)
    for name, sel in (("NODE", kind == NODE), ("EDGE", kind == EDGE)):
        s = sel & live
        if not s.any():
            print(f"    {name:<14} {'(no live slots)':>14}")
            continue
        pi = 0.5
        mu_d = pi * m_g.m1[s] + (1 - pi) * m_r.m1[s]
        mu_s = pi * m_g.m2[s] + (1 - pi) * m_r.m2[s]
        v_dd = pi * m_g.q1[s] + (1 - pi) * m_r.q1[s] - mu_d * mu_d
        v_ss = pi * m_g.q2[s] + (1 - pi) * m_r.q2[s] - mu_s * mu_s
        v_ds = pi * m_g.q12[s] + (1 - pi) * m_r.q12[s] - mu_d * mu_s
        denom = np.maximum(v_dd * v_ss, 1e-300)
        headroom = 1.0 - (v_ds * v_ds) / denom
        print(f"    {name:<14} {float(headroom.min()):>14.3e} {float(np.median(headroom)):>16.3e} "
              f"{float(np.percentile(headroom, 1)):>14.3e}")

    # ── ③ is the claim TRUE, stratified by N ──
    print("\n  ③ IS THE CLAIM TRUE?  |argmax − truth|, and whether tau is EARNED, by the count N")
    print("     ⚠ the recorded bias ladder predicts 0.32 at N=1, 0.05 at N=5, 0.004 at N=50")
    print("     ⭐ the last three columns are the SAME statistics with the RNA edge opportunity REPAIRED")
    print(f"    {'slot type':<8} {'N':<10} {'slots':>10} {'fragments':>13} {'mean|Δ|':>9} "
          f"{'med|Δ|':>9} {'bias':>9} {'med tau':>9} {'conf-wrong':>11}  |{'bias':>9} "
          f"{'mean|Δ|':>8} {'cf-wrong':>10}")
    print("    " + "-" * 126)
    scored = live & np.isfinite(fg_true)
    for name, tsel in (("NODE", kind == NODE), ("EDGE", kind == EDGE)):
        for lo, hi in _BINS:
            s = scored & tsel & (count >= lo) & (count <= hi)
            if not s.any():
                continue
            err = argmax[s] - fg_true[s]
            w = tot_true[s]
            # ⭐ "confidently wrong": the channel's own sd says the truth is >2 sd away. A tight wrong
            # value is what outvotes correct neighbours, so it is counted, not averaged away.
            sd = 1.0 / np.sqrt(np.maximum(tau[s], 1e-12))
            conf_wrong = float(w[np.abs(err) > 2.0 * sd].sum()) / max(float(w.sum()), 1.0)
            label = f"{lo}" if lo == hi else (f"{lo}+" if hi > 10**8 else f"{lo}-{hi}")
            # ⭐ the SAME statistics on the repaired channel, so the two are read side by side
            sf = live_fix & np.isfinite(fg_true) & tsel & (count >= lo) & (count <= hi)
            ef = argmax_fix[sf] - fg_true[sf]
            wf = tot_true[sf]
            sdf = 1.0 / np.sqrt(np.maximum(tau_fix[sf], 1e-12))
            cwf = float(wf[np.abs(ef) > 2.0 * sdf].sum()) / max(float(wf.sum()), 1.0)
            print(f"    {name:<8} {label:<10} {int(s.sum()):>10,} {float(w.sum()):>13,.0f} "
                  f"{float(np.abs(err).mean()):>9.4f} {float(np.median(np.abs(err))):>9.4f} "
                  f"{float(err.mean()):>+9.4f} {float(np.median(tau[s])):>9.3f} "
                  f"{100 * conf_wrong:>10.1f}%  |{float(ef.mean()):>+9.4f} "
                  f"{float(np.abs(ef).mean()):>8.4f} {100 * cwf:>9.1f}%")

    # ── ④ the DERIVED smooth shrinkage against the tabulated precision ──
    print("\n  ④ ⭐⭐ SMOOTH SHRINKAGE — the analytic Fisher information vs the tabulated tau")
    print("     I_lam = N * Delta' V^-1 Delta * [pi(1-pi)]^2. It vanishes QUADRATICALLY as the two pmfs")
    print(f"     converge and has no grid floor. ⛔ tau's floor is 1/Var(lam grid) = "
          f"{1.0 / float(np.var(lam_grid)):.4f} — anything at that value is the GRID, not evidence.")
    info = fisher_information(m_g, m_r, count, fg_grid)
    print(f"    {'slot type':<10} {'N':<10} {'slots':>10} {'med tau':>10} {'med I_lam':>12} "
          f"{'I/tau':>9} {'tau at floor':>13}")
    print("    " + "-" * 80)
    floor = 1.0 / float(np.var(lam_grid))
    for name, tsel in (("NODE", kind == NODE), ("EDGE", kind == EDGE)):
        for lo, hi in _BINS:
            s_ = live & tsel & (count >= lo) & (count <= hi)
            if not s_.any():
                continue
            mt, mi = float(np.median(tau[s_])), float(np.median(info[s_]))
            at_floor = float(np.mean(np.abs(tau[s_] - floor) < 1e-4))
            label = f"{lo}" if lo == hi else (f"{lo}+" if hi > 10**8 else f"{lo}-{hi}")
            print(f"    {name:<10} {label:<10} {int(s_.sum()):>10,} {mt:>10.4f} {mi:>12.4f} "
                  f"{mi / max(mt, 1e-12):>9.2f} {100 * at_floor:>12.1f}%")

    # ── ⑤ the model against the population it describes ──
    print("\n  ⑤ ⭐⭐ DOES EACH COMPONENT'S POPULATION HAVE THE MOMENTS THE CHANNEL ASSUMES?")
    print("     realised = the origin-split bank's own Sum_u/N and Sum_w/N; predicted = the shipped")
    print("     moments on the fitted pmf, count-weighted over the same slots.")
    print("     ⭐ gDNA CANNOT SPLICE, so it is the control: if RNA is off at EDGE and gDNA is not, the")
    print("     mechanism is the spliced crossings this bank excludes — not the pmf fit.")
    print(f"    {'component':<10} {'axis':<6} {'fragments':>13} {'m2 real':>9} {'m2 pred':>9} "
          f"{'ratio':>7} {'m1 real':>10} {'m1 pred':>10} {'ratio':>7}")
    print("    " + "-" * 92)
    for comp, origins, mom in (("gdna", ("gdna",), m_g), ("rna", ("mrna", "nrna"), m_r)):
        c_c, u_c, w_c = slot_channels(oracle, chain, origins)
        for name, tsel in (("NODE", kind == NODE), ("EDGE", kind == EDGE)):
            s_ = tsel & (c_c > 0)
            if not s_.any():
                continue
            tot = float(c_c[s_].sum())
            m2_real, m1_real = float(w_c[s_].sum()) / tot, float(u_c[s_].sum()) / tot
            m2_pred = float((c_c[s_] * np.asarray(mom.m2)[s_]).sum()) / tot
            m1_pred = float((c_c[s_] * np.asarray(mom.m1)[s_]).sum()) / tot
            print(f"    {comp:<10} {name:<6} {tot:>13,.0f} {m2_real:>9.2f} {m2_pred:>9.2f} "
                  f"{m2_real / max(m2_pred, 1e-12):>7.3f} {m1_real:>10.5f} {m1_pred:>10.5f} "
                  f"{m1_real / max(m1_pred, 1e-12):>7.3f}")
    print("    ⚠ ratio 1.000 = the channel's model of that population is exact. m2 is the mean length,")
    print("      m1 the mean reciprocal weight; they are the two coordinates the likelihood scores on.")

    # ── ⑥ WHICH REPAIR? the fitted pmf, or a missing selection term ──
    print("\n  ⑥ ⭐⭐⭐ WHICH REPAIR IS NEEDED — a better FIT, or a missing SELECTION term?")
    print("     Same comparison as ⑤, with the predicted moments rebuilt on the condition's TRUE")
    print("     per-component length histogram instead of the fitted pmf.")
    print("     ⭐ ratio -> 1.000 here means the FIT is the whole problem. A residual means the channel's")
    print("     model of its own population is missing a term, and no better fit removes it.")
    if "capture_on" in args.condition:
        print("    ⛔⛔ NOT INTERPRETABLE ON THIS CONDITION. `truth_fragment_lengths.tsv` is the GENERATED")
        print("    library and hybrid capture LENGTH-SELECTS what is then sequenced, so under capture the")
        print("    generated pmf is not the landed pmf and this arm conflates capture with everything")
        print("    else. Measured: the gDNA control reads 1.000 with BOTH pmfs off capture and the RNA")
        print("    arm goes 0.969 -> 0.860 under it, which is the selection, not a model error.")
        print("    ⭐ Read table ⑤ under capture; read ⑥ only off it.")
    tp = truth_pmfs(suite, args.condition, int(payload.max_length))
    print(f"    {'component':<10} {'axis':<6} {'m2 real':>9} {'m2 FITTED':>10} {'ratio':>7} "
          f"{'m2 TRUE':>9} {'ratio':>7} {'verdict':>26}")
    print("    " + "-" * 92)
    for comp, origins, mom in (("gdna", ("gdna",), m_g), ("rna", ("mrna", "nrna"), m_r)):
        true_m = build_slot_moments(chain, ra, tp[comp])
        c_c, _u_c, w_c = slot_channels(oracle, chain, origins)
        for name, tsel in (("NODE", kind == NODE), ("EDGE", kind == EDGE)):
            s_ = tsel & (c_c > 0)
            if not s_.any():
                continue
            tot = float(c_c[s_].sum())
            real = float(w_c[s_].sum()) / tot
            fit = float((c_c[s_] * np.asarray(mom.m2)[s_]).sum()) / tot
            tru = float((c_c[s_] * np.asarray(true_m.m2)[s_]).sum()) / tot
            r_fit, r_tru = real / max(fit, 1e-12), real / max(tru, 1e-12)
            closed = abs(r_tru - 1.0) < 0.25 * abs(r_fit - 1.0)
            verdict = "the FIT" if closed else ("a SELECTION term" if abs(r_tru - 1.0) > 0.01 else "both, fit dominant")
            print(f"    {comp:<10} {name:<6} {real:>9.2f} {fit:>10.2f} {r_fit:>7.3f} "
                  f"{tru:>9.2f} {r_tru:>7.3f} {verdict:>26}")

    # ── ⑦ WHERE ON THE NODE AXIS DOES THE OPPORTUNITY MODEL DEPART? ──
    print("\n  ⑦ ⭐⭐⭐ THE NODE-AXIS RESIDUAL, STRATIFIED — the tilt is (ell − w + 1), so it FLATTENS as")
    print("     ell grows. At long nodes the model reduces to 'contained fragments have the library")
    print("     distribution', so a residual that SURVIVES there is the population, and one that only")
    print("     appears at short ell is the tilt/placement.")
    print("     ⭐ gDNA is the control for the TILT ITSELF: it is uniformly placed, so if (ell−w+1) were")
    print("     the wrong shape gDNA would be off too.")
    from rigel.calibration.signature import coarse_type_array

    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    node_len_all = np.asarray(ra.region_size_bp, np.float64)
    slot_ell = np.zeros(int(chain.n_slots), np.float64)
    slot_rt = np.full(int(chain.n_slots), -1, np.int64)
    obj_all = np.asarray(chain.obj_idx, np.int64)
    slot_ell[kind == NODE] = node_len_all[obj_all[kind == NODE]]
    slot_rt[kind == NODE] = rtype[obj_all[kind == NODE]]
    _ELL = ((0, 149), (150, 299), (300, 599), (600, 1499), (1500, 10**9))
    print(f"    {'component':<10} {'ell':<12} {'slots':>9} {'fragments':>13} {'m2 real':>9} "
          f"{'m2 pred':>9} {'ratio':>7} {'m1 ratio':>9}")
    print("    " + "-" * 84)
    for comp, origins, mom in (("gdna", ("gdna",), m_g), ("rna", ("mrna", "nrna"), m_r)):
        c_c, u_c, w_c = slot_channels(oracle, chain, origins)
        for lo, hi in _ELL:
            s_ = (kind == NODE) & (c_c > 0) & (slot_ell >= lo) & (slot_ell <= hi)
            if not s_.any():
                continue
            tot = float(c_c[s_].sum())
            m2r, m1r = float(w_c[s_].sum()) / tot, float(u_c[s_].sum()) / tot
            m2p = float((c_c[s_] * np.asarray(mom.m2)[s_]).sum()) / tot
            m1p = float((c_c[s_] * np.asarray(mom.m1)[s_]).sum()) / tot
            label = f"{lo}+" if hi > 10**8 else f"{lo}-{hi}"
            print(f"    {comp:<10} {label:<12} {int(s_.sum()):>9,} {tot:>13,.0f} {m2r:>9.2f} "
                  f"{m2p:>9.2f} {m2r / max(m2p, 1e-12):>7.3f} {m1r / max(m1p, 1e-12):>9.3f}")
    print()
    print(f"    {'component':<10} {'node type':<12} {'slots':>9} {'fragments':>13} {'m2 real':>9} "
          f"{'m2 pred':>9} {'ratio':>7}")
    print("    " + "-" * 74)
    for comp, origins, mom in (("gdna", ("gdna",), m_g), ("rna", ("mrna", "nrna"), m_r)):
        c_c, _u, w_c = slot_channels(oracle, chain, origins)
        for rt, rname in ((0, "intergenic"), (1, "intron"), (2, "exon")):
            s_ = (kind == NODE) & (c_c > 0) & (slot_rt == rt)
            if not s_.any():
                continue
            tot = float(c_c[s_].sum())
            m2r = float(w_c[s_].sum()) / tot
            m2p = float((c_c[s_] * np.asarray(mom.m2)[s_]).sum()) / tot
            print(f"    {comp:<10} {rname:<12} {int(s_.sum()):>9,} {tot:>13,.0f} {m2r:>9.2f} "
                  f"{m2p:>9.2f} {m2r / max(m2p, 1e-12):>7.3f}")

    # ── ⑨ THE 2x2: the RNA opportunity repair x a gDNA pmf that is right under capture ──
    print("\n  ⑨ ⭐⭐⭐ THE 2x2 — do the two repairs CANCEL? Bias and mean|Δ| against truth at N >= 50.")
    print("     ⛔ A factorial, because §8 showed the RNA repair alone helps off capture and HURTS under")
    print("     it, and two effects that might cancel cannot be read one at a time.")
    print("     The gDNA column emulates a fixed capture defect by feeding the TRUE gDNA pmf.")
    tp2 = truth_pmfs(suite, args.condition, int(payload.max_length))
    _dl2, _dr2 = exonic_block_distances(chain, ra)
    rna_arms = {
        "A = w−1": _bsm(chain, ra, fl.rna_pmf),
        "A = unspliced": fixed_r,
    }
    gdna_arms = {
        "fitted": _bsm(chain, ra, fl.gdna_pmf),
        "TRUE pmf": _bsm(chain, ra, tp2["gdna"]),
    }
    print(f"\n    {'gDNA pmf':<12}{'RNA opportunity':<18}{'slot':<7}{'bias':>10}{'mean|Δ|':>10}"
          f"{'conf-wrong':>12}")
    print("    " + "-" * 70)
    for gname, gm in gdna_arms.items():
        for rname, rm in rna_arms.items():
            L = _ll(gm, rm, count,
                    np.asarray(geometry.unspliced_inv_length_sum, np.float64),
                    np.asarray(geometry.unspliced_length_sum, np.float64), fg_grid)
            lv = np.ptp(L, axis=1) > 0.0
            am = fg_grid[np.argmax(L, axis=1)]
            tt = density_factor_precision(L, lam_grid)
            for sname, tsel in (("EDGE", kind == EDGE), ("NODE", kind == NODE)):
                sel = lv & np.isfinite(fg_true) & tsel & (count >= 50)
                if not sel.any():
                    continue
                e = am[sel] - fg_true[sel]
                ww = tot_true[sel]
                sd = 1.0 / np.sqrt(np.maximum(tt[sel], 1e-12))
                cw = float(ww[np.abs(e) > 2.0 * sd].sum()) / max(float(ww.sum()), 1.0)
                print(f"    {gname:<12}{rname:<18}{sname:<7}{float(e.mean()):>+10.4f}"
                      f"{float(np.abs(e).mean()):>10.4f}{100 * cw:>11.1f}%")
    print("    ⭐ read the EDGE rows: if 'TRUE pmf x A = unspliced' is the best cell, the two repairs")
    print("      belong together and the gDNA capture defect must land FIRST.")

    # ── ⑧ THE DERIVED UNSPLICED OPPORTUNITY, as an arm ──
    print("\n  ⑧ ⭐⭐⭐ THE UNSPLICED CROSSING OPPORTUNITY — the derived repair, as an arm")
    print("     A(w) = #{u : max(1,w−d_right) <= u <= min(w−1,d_left)}, with d_* the distance to the")
    print("     ends of the contiguous EXONIC block. It drops placements whose fragment would have run")
    print("     past a junction (those land in `edge_spliced`) and shrinks where the exon is short.")
    print("     ⛔ NO free parameter: the formula is fully determined, so it lands on 1.000 or it is wrong.")
    c_r, _u_r, w_r = slot_channels(oracle, chain, ("mrna", "nrna"))
    s_ = (kind == EDGE) & (c_r > 0)
    tot = float(c_r[s_].sum())
    real = float(w_r[s_].sum()) / tot
    old = float((c_r[s_] * np.asarray(m_r.m2)[s_]).sum()) / tot
    print(f"\n    {'RNA, EDGE slots':<40}{'m2':>10}{'ratio real/pred':>18}")
    print("    " + "-" * 70)
    print(f"    {'realised (the bank)':<40}{real:>10.2f}{'':>18}")
    print(f"    {'predicted, A = w−1 (today)':<40}{old:>10.2f}{real / max(old, 1e-9):>18.3f}")
    arms = {}
    for label, pure in (("A = unspliced, COARSE block", False),
                        ("A = unspliced, PURELY-exonic block", True)):
        dl, dr = exonic_block_distances(chain, ra, pure=pure)
        bm = bounded_crossing_moments(fl.rna_pmf, dl, dr)
        arms[pure] = bm
        v = float((c_r[s_] * np.asarray(bm.m2)[s_]).sum()) / tot
        print(f"    {'predicted, ' + label:<40}{v:>10.2f}{real / max(v, 1e-9):>18.3f}")
    bm_r = arms[True]
    dead = float(c_r[(kind == EDGE) & (c_r > 0) & (np.asarray(bm_r.eff) <= 0)].sum())
    print(f"\n    ⚠ RNA at edges the PURELY-exonic opportunity calls IMPOSSIBLE: {dead:,.0f} "
          f"({100 * dead / max(tot, 1):.2f} %)")
    print("      A nonzero share means mature RNA is crossing a line its exonic block should forbid —")
    print("      nascent RNA, an annotation gap, or the block walk is wrong.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
