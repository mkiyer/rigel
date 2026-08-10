#!/usr/bin/env python
"""⭐⭐⭐ **CAN A SHORT EXON GIVE AN UNBIASED RNA LENGTH SAMPLE? — the measurement before the schema change.**

**The problem.** ``rna_pmf`` is fitted on ``POOL_RNA_SPLICED``, the only composition-pure RNA pool, and
that pool is length-SELECTED: a longer fragment crosses a junction more often, so the pool
over-represents it. The correction is ``pi(w)``, the modelled probability that a uniformly placed
length-``w`` window crosses a junction. Measured 2026-08-09 on the production (DRAINED) substrate,
``pi(w)`` is about **2× too flat** and the corrected pmf comes out **+10.7 %** long — 211.60 bp against
a true library of 191.06. Two recorded explanations were eliminated: it is not the uniform ``theta``
(true molar abundances change nothing), and it is not the "splice must be SEEN" requirement (that makes
selection *decrease* with length, the wrong sign).

⭐⭐ **THE IDEA (owner, 2026-08-09), AND IT REMOVES THE MODEL RATHER THAN IMPROVING IT.** Selection is
per-transcript geometry, and one geometry makes it trivial: a fragment of length ``w`` touching an
INTERNAL exon of length ``e < w`` **cannot be contained** — it must extend past a boundary, and both
boundaries of an internal exon are junctions. So it is spliced with probability exactly **1**. The
length bias vanishes by construction instead of being corrected.

⭐ Both quantities are known per fragment, so the population needs **no tunable**: *spliced fragments
for which containment was geometrically impossible*. And short exons are exactly where the current
correction is worst — ``pi`` reads 0.244 at 50–109 bp, the short end.

⛔ **WHY THIS SCRIPT EXISTS RATHER THAN A PATCH.** Reading that population in production needs a bank
the accumulator does not have (``sj_length_sum`` was deleted 2026-08-08 for want of a consumer). A
schema change costs a digest bump plus a ~2 h rebuild of 176 caches, so it is justified by a
measurement or not at all. Everything here is computed OFFLINE from the simulator's read names — which
carry ``{t_id}:{frag_start}-{frag_end}`` in TRANSCRIPT coordinates (``sim/reads.py``) — and the index's
own exon structure. No accumulator, no calibration, no solver.

⭐ **THREE ARMS, and the middle one is the incumbent:**

===  ==========================================================  ==========================
 T   the TRUE library — every mRNA fragment, straight from       the target
     the read names
 P   the SPLICED pool de-tilted by the shipped ``pi(w)``          what production does today
 S   the SHORT-EXON population de-tilted by ``w − 1 + e_bar``     the proposal
===  ==========================================================  ==========================

⛔ **GATED ON ITS OWN PREMISE FIRST.** If the read-name spans do not reproduce the condition's
``truth_fragment_lengths.tsv``, then they are not fragment spans and every number below is meaningless;
the script refuses rather than reporting.

Usage::

    python scripts/design/short_exon_fl_probe.py --panel flgap_short --condition ..._capture_off
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pysam  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
for _p in (_REPO / "scripts" / "design", _REPO / "tests" / "calibration"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

RUNS = Path.home() / "Downloads" / "rigel_runs"
INDEX = RUNS / "suite" / "rigel_index"


def read_fragments(bam: str):
    """``(t_id, start, w)`` for every mRNA fragment, from the read names. One entry per qname group."""
    from rigel.sim.read_name import parse_origin

    tids: list[str] = []
    starts: list[int] = []
    lens: list[int] = []
    current = None
    with pysam.AlignmentFile(bam, "rb") as fh:
        for rec in fh:
            q = rec.query_name
            if q == current:
                continue
            current = q
            o = parse_origin(q)
            if o.kind != "mrna":
                continue
            tids.append(o.transcript_id)
            starts.append(int(o.start))
            lens.append(int(o.end) - int(o.start))
    return tids, np.asarray(starts, np.int64), np.asarray(lens, np.int64)


def transcript_exons(index):
    """Per transcript, exon lengths in TRANSCRIPT order (5'→3').

    ⛔ A minus-strand transcript's transcript coordinates run against the genome, so its exons must be
    reversed. Getting this wrong would mis-place every junction on half the annotation, and it would
    look like a modelling result rather than a bug.
    """
    offsets, es, ee, _ = index.build_exon_csr()
    offsets = np.asarray(offsets)
    length = np.asarray(ee) - np.asarray(es)
    strand = index.t_df["strand"].to_numpy()
    out = {}
    for t, tid in enumerate(index.t_df["t_id"].to_numpy()):
        a, b = int(offsets[t]), int(offsets[t + 1])
        e = length[a:b]
        out[str(tid)] = e[::-1].copy() if str(strand[t]) == "-" else e.copy()
    return out


def main() -> int:
    from rigel.calibration.junction_opportunity import crossing_probability_from_index
    from rigel.index import TranscriptIndex

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--panel", default="flgap_short")
    ap.add_argument("--condition", default="gdna_g50_ss_0.50_nrna_none_capture_off")
    ap.add_argument("--max-exon", type=int, default=0,
                    help="0 = the derived rule (e < w, per fragment). A positive value additionally "
                         "restricts the qualifying exons to e <= this, which FIXES the exon set so it "
                         "cannot vary with w.")
    args = ap.parse_args()

    suite = RUNS / "suite" / args.panel
    bam = str(suite / args.condition / "sim_oracle.bam")
    index = TranscriptIndex.load(str(INDEX))

    tids, starts, w = read_fragments(bam)
    max_size = int(w.max())
    print()
    print("=" * 100)
    print(f"  ⭐⭐⭐ THE SHORT-EXON RNA LENGTH SAMPLE — {args.panel}/{args.condition}")
    print(f"  {len(tids):,} mRNA fragments from the read names")
    print("=" * 100)

    # ── GATE: are these really fragment spans? ──
    truth = np.zeros(max_size + 2)
    with (suite / args.condition / "truth_fragment_lengths.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["kind"] == "mrna":
                L = int(r["fragment_length"])
                if 0 <= L <= max_size:
                    truth[L] += float(r["count"])
    obs = np.bincount(w, minlength=max_size + 2).astype(float)
    m_obs = float((np.arange(obs.size) * obs).sum() / obs.sum())
    m_tru = float((np.arange(truth.size) * truth).sum() / max(truth.sum(), 1))
    print("\n  ⛔ PREMISE GATE — read-name spans vs the condition's own truth table")
    print(f"     mean {m_obs:.2f} bp observed against {m_tru:.2f} bp declared, "
          f"n {obs.sum():,.0f} vs {truth.sum():,.0f}")
    if abs(m_obs - m_tru) > 1.0:
        raise SystemExit(
            "  ⛔⛔ REFUSING: the read-name spans are not the fragment spans. Every arm below would be "
            "measuring something else."
        )
    print("     ✅ they are fragment spans; the arms below are readable")

    # ── classify every fragment against its own transcript's exon layout ──
    exons = transcript_exons(index)
    by_t = defaultdict(list)
    for i, t in enumerate(tids):
        by_t[t].append(i)

    all_h = np.zeros(max_size + 2)
    spl_h = np.zeros(max_size + 2)
    inc_h = np.zeros(max_size + 2)  # (fragment, qualifying exon) incidences
    e_sum = 0.0
    e_n = 0.0
    n_missing = 0
    for t, idx in by_t.items():
        e = exons.get(t)
        ii = np.asarray(idx, np.int64)
        if e is None or e.size == 0:
            n_missing += ii.size
            continue
        p, ww = starts[ii], w[ii]
        np.add.at(all_h, ww, 1.0)
        cum = np.concatenate(([0], np.cumsum(e)))          # exon i spans [cum[i], cum[i+1])
        lo, hi = cum[:-1][None, :], cum[1:][None, :]        # (1,K)
        P, W = p[:, None], ww[:, None]
        overlaps = (P < hi) & (P + W > lo)
        contained = (P >= lo) & (P + W <= hi)
        np.add.at(spl_h, ww[~contained.any(axis=1)], 1.0)
        internal = np.zeros(e.size, bool)
        if e.size >= 3:
            internal[1:-1] = True
        elen = e[None, :]
        qualifies = overlaps & internal[None, :] & (elen < W)
        if args.max_exon > 0:
            qualifies &= elen <= args.max_exon
        k = qualifies.sum(axis=1)
        np.add.at(inc_h, ww, k.astype(float))
        e_sum += float((np.broadcast_to(elen, qualifies.shape)[qualifies]).sum())
        e_n += float(k.sum())
    if n_missing:
        print(f"     ⚠ {n_missing:,} fragments on transcripts absent from the index — dropped")

    e_bar = e_sum / max(e_n, 1.0)
    wv = np.arange(all_h.size, dtype=float)
    pi = np.asarray(crossing_probability_from_index(index, all_h.size - 1), float)

    def norm(h):
        s = h.sum()
        return h / s if s > 0 else h

    def mean(h):
        s = h.sum()
        return float((wv * h).sum() / s) if s > 0 else float("nan")

    arm_T = norm(all_h)
    arm_P = norm(np.divide(spl_h, pi, out=np.zeros_like(spl_h), where=pi > 0))
    arm_S = norm(np.divide(inc_h, np.maximum(wv - 1.0 + e_bar, 1.0)))

    print(f"\n  qualifying (fragment, internal exon with e < w) incidences: {e_n:,.0f}   "
          f"mean qualifying exon e_bar = {e_bar:.1f} bp")
    print(f"  spliced fragments: {spl_h.sum():,.0f} of {all_h.sum():,.0f} "
          f"({100 * spl_h.sum() / max(all_h.sum(), 1):.1f} %)")
    print()
    print(f"  {'arm':<46}{'mean bp':>10}{'err vs T':>11}")
    print("  " + "-" * 68)
    print(f"  {'T  the TRUE library (all mRNA fragments)':<46}{mean(arm_T):>10.2f}{'':>11}")
    print(f"  {'P  spliced pool / pi(w)   [production today]':<46}{mean(arm_P):>10.2f}"
          f"{mean(arm_P) - mean(arm_T):>+11.2f}")
    print(f"  {'S  short-exon incidences / (w-1+e_bar)':<46}{mean(arm_S):>10.2f}"
          f"{mean(arm_S) - mean(arm_T):>+11.2f}")

    # ── THE EDGE MECHANISM: is the unspliced crossing population shorter, and by how much? ──
    print("\n  ⭐⭐⭐ THE EDGE MECHANISM — does spliced EXCLUSION explain the crossing deficit?")
    print("     `crossing_moments` tilts the library pmf by the crossing opportunity (w−1) and predicts")
    print("     m2 = 246.57 at EDGE slots. The bank `edge_unspliced` holds 228.01 — a −7.5 % deficit.")
    print("     ⛔ The bank excludes fragments that spliced ELSEWHERE (they land in `edge_spliced`), and")
    print("     longer RNA splices more. So the SAME tilt over the UNSPLICED subset should land on 228.")
    tilt = np.maximum(wv - 1.0, 0.0)
    uns_h = all_h - spl_h
    def tilted(h):
        num = float((tilt * wv * h).sum())
        den = float((tilt * h).sum())
        return num / den if den > 0 else float("nan")
    print(f"\n    {'population':<44}{'(w−1)-tilted mean':>19}{'target':>10}")
    print("    " + "-" * 74)
    print(f"    {'ALL RNA fragments  (what the model assumes)':<44}{tilted(all_h):>19.2f}{246.57:>10.2f}")
    print(f"    {'UNSPLICED only     (what the bank holds)':<44}{tilted(uns_h):>19.2f}{228.01:>10.2f}")
    print(f"    {'SPLICED only       (the excluded population)':<44}{tilted(spl_h):>19.2f}{'':>10}")
    print(f"\n    {'w':>10}{'P(spliced | w)':>16}{'share of RNA':>14}")
    print("    " + "-" * 42)
    for lo, hi in ((50, 79), (80, 109), (110, 149), (150, 199), (200, 249), (250, 299), (300, 399),
                   (400, 599)):
        a = all_h[lo:hi + 1].sum()
        if a <= 0:
            continue
        print(f"    {str(lo) + '-' + str(hi):>10}{spl_h[lo:hi + 1].sum() / a:>16.3f}"
              f"{a / all_h.sum():>14.4f}")
    print("    ⭐ a RISING P(spliced|w) is the mechanism: it removes the long tail from the bank.")

    print("\n  ⭐ PER-BIN SHAPE — arm / T, so 1.000 is exact. The short end is what pi(w) gets wrong.")
    print(f"    {'w':>10}{'T share':>10}{'P / T':>9}{'S / T':>9}")
    print("    " + "-" * 40)
    for lo, hi in ((50, 79), (80, 109), (110, 149), (150, 199), (200, 249), (250, 299), (300, 399),
                   (400, 599)):
        tT = arm_T[lo:hi + 1].sum()
        if tT <= 0:
            continue
        print(f"    {str(lo) + '-' + str(hi):>10}{tT:>10.4f}"
              f"{arm_P[lo:hi + 1].sum() / tT:>9.3f}{arm_S[lo:hi + 1].sum() / tT:>9.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
