#!/usr/bin/env python
"""⭐⭐⭐ THE REFRAME, WALKED — every count, both flank totals, and every hop in BOTH directions.

⛔ **This is not a summary statistic.** It prints one toy transcript's whole chain and then walks the
message layer hop by hop, so the include/exclude decision at each `(BOUNDARY, side)` is visible as a number
rather than inferred from an aggregate.

**Five sections:**

1. **THE COUNTS** — every REGION and every BOUNDARY: the unspliced crossing per genome strand, the spliced
   crossing, the sj flux per transcript strand, and the three divisors. Plus the sj axis.
2. ⭐⭐ **THE FLANK PAIR** — per slot, ``rho_unspliced``, ``rho_lo``, ``rho_hi``, and which side the
   sj flux landed on. ⛔ This is the whole change: a slot where the two differ is a slot the
   predecessor gave ONE number to.
3. ⭐⭐⭐ **THE FORWARD RELAY (low → high)** and **4. THE BACKWARD RELAY (high → low)**, hop by hop:
   which total each end presented, the reframe ``r``, what ``r`` would have been under the shipped
   single-total, what TRUTH says the same ratio is, whether the composition licence fired, and the gDNA
   level arriving and leaving. ⭐ A hop where ``r_new`` tracks truth and ``r_old`` does not is the
   derivation working, in one boundary.
5. **THE ANSWER** per slot against per-object truth.

⚠ **The geometry is REBUILT here, not read out of the solver** — ``build_region_geometry`` is a pure
function of (chain, substrate, region_arrays, sj, two pmfs, reach) and this calls it with the same
arguments `calibrate` did. ⛔ That re-derivation is then GATED against the frames the solver published
(``_uni_static['rho_lo'/'rho_hi']``) to 1e-12, so it cannot silently drift (TRAPS: self-checking-validator).

⚠ **`r_true` is built with the SAME flank rule from the ORACLE's counts**, so it is the ratio of the
populations the hop actually compares — not a different quantity. ``r_true_gdna`` beside it is the true
gDNA-density ratio, which is what a gDNA LEVEL should cross at (TRAPS: a-total-density-ratio) and is a DIFFERENT
number; both are printed because conflating them is TRAPS: recompute-from-the-oracle.

Usage:

    python scripts/design/reframe_walk.py --n-rna 200000            # all three gDNA arms
    python scripts/design/reframe_walk.py --n-rna 200000 --arms high
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import math
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

DESIGN = Path(__file__).resolve().parent
_s = importlib.util.spec_from_file_location("toy_harness", DESIGN / "toy_harness.py")
TH = importlib.util.module_from_spec(_s)
sys.modules["toy_harness"] = TH
_s.loader.exec_module(TH)

from rigel.calibration.region_chain import REGION  # noqa: E402
from rigel.calibration.region_geometry import (  # noqa: E402
    build_region_geometry,
    region_total_density,
)
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.calibration.splice_graph import build_sj_geometry_arrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.types import Strand  # noqa: E402

SUITE = Path.home() / "Downloads/rigel_runs/suite/ladder"
INDEX = Path.home() / "Downloads/rigel_runs/suite/rigel_index"
TYPES = {0: "intergenic", 1: "intron", 2: "exon"}
EPS = 1e-12

#: ⭐ The gDNA level is DERIVED from the donor, never set here (`TESTING.md` §0b), so "how much gDNA"
#: is chosen by picking a donor. All three are capture-OFF and UNSTRANDED, which is the simplest
#: possible regime: no enrichment landscape and exactly zero strand information (`EQUATIONS.md` §5.2).
ARMS = {
    # ⚠ `high` was `g75` until 2026-08-13. The ladder was rebuilt with four rungs (g00/g05/g50/g98)
    # and `g75` no longer exists, so this is the surviving middle — three DISTINCT gDNA levels is
    # what the three arms need, and 0 / 50 / 98 % still spans the range.
    "zero": "gdna_g00_ss_0.50_nrna_none_capture_off",
    "high": "gdna_g50_ss_0.50_nrna_none_capture_off",
    "veryhigh": "gdna_g98_ss_0.50_nrna_none_capture_off",
}


def labels(chain, ra):
    """A readable name per slot: the coarse type at a REGION, the flanking pair at an BOUNDARY."""
    kind, obj = np.asarray(chain.kind), np.asarray(chain.obj_idx, np.int64)
    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    starts, sizes = np.asarray(ra.start, np.int64), np.asarray(ra.region_size_bp, np.int64)
    out = []
    for s in range(int(chain.n_slots)):
        i = int(obj[s])
        if kind[s] == REGION:
            out.append(f"{TYPES[int(rtype[i])]} [{starts[i]:,}–{starts[i] + sizes[i]:,})")
        else:
            hi, lo = s + 1, s - 1
            b = int(rtype[obj[hi]]) if hi < int(chain.n_slots) and kind[hi] == REGION else -1
            a = int(rtype[obj[lo]]) if lo >= 0 and kind[lo] == REGION else -1
            pair = "|".join(TYPES.get(x, "?") for x in sorted((a, b)) if x >= 0)
            out.append(f"{pair} @{starts[obj[hi]]:,}" if b >= 0 else f"boundary #{i}")
    return out


def rebuild_geometry(r):
    """`build_region_geometry` with the arguments `calibrate` used — see the module docstring's gate."""
    ra = r.region_arrays
    sub = CalibrationSubstrate.from_payload(r.payload, ra)
    return build_region_geometry(
        r.chain,
        sub,
        ra,
        build_sj_geometry_arrays(r.index),
        r.donor.gdna_fl_pmf,
        r.donor.rna_fl_pmf,
    )


def truth_slot_arrays(r):
    """Per SLOT, from the oracle's origin split: gDNA / RNA unspliced counts, and the sj flux.

    ⚠ Built from ``truth.parts`` directly rather than from ``override_masses``, because that helper folds
    the SPLICED crossing into the RNA boundary mass and `region_total_density` deliberately excludes it from
    both flank totals — so using it here would compare two different populations."""
    ra = r.region_arrays
    subs = {k: CalibrationSubstrate.from_payload(r.truth.parts[k], ra) for k in ("gdna", "mrna", "nrna")}
    full = CalibrationSubstrate.from_payload(r.truth.full, ra)
    kind, obj = np.asarray(r.chain.kind), np.asarray(r.chain.obj_idx, np.int64)
    n = int(r.chain.n_slots)

    def per_slot(sub, region_pop, boundary_pop):
        nd = np.asarray(getattr(sub, region_pop).count, np.float64).sum(1)
        eg = np.asarray(getattr(sub, boundary_pop).count, np.float64).sum(1)
        out = np.zeros(n)
        for s in range(n):
            i = int(obj[s])
            out[s] = nd[i] if kind[s] == REGION else eg[i]
        return out

    g = per_slot(subs["gdna"], "region_contained", "boundary_unspliced")
    rna = per_slot(subs["mrna"], "region_contained", "boundary_unspliced") + per_slot(
        subs["nrna"], "region_contained", "boundary_unspliced"
    )
    return g, rna, full


def _flux_density(count, eff):
    """Σ_s count_s / eff_s over the two TRANSCRIPT-strand columns — 0 where a strand has no flux."""
    c, e = np.asarray(count, np.float64), np.asarray(eff, np.float64)
    out = np.zeros(c.shape[0])
    for k in (0, 1):
        out += np.where(e[:, k] > EPS, c[:, k] / np.where(e[:, k] > EPS, e[:, k], 1.0), 0.0)
    return out


def truth_flank_pair(r, geom, g_cnt, rna_cnt):
    """TRUTH's own ``(rho_lo, rho_hi)``, built with the SAME flank rule from the oracle's counts.

    ⭐ The sj flux is a MEASUREMENT — gDNA cannot splice — so its density needs no deconvolution
    and truth and solver share it exactly. What differs is the unspliced part, where truth knows the
    split and the solver has to derive it."""
    E_g = np.asarray(geom.eff_gdna, np.float64)
    E_r = np.asarray(geom.eff_rna, np.float64)
    unspl = np.where(E_g > EPS, g_cnt / np.where(E_g > EPS, E_g, 1.0), 0.0) + np.where(
        E_r > EPS, rna_cnt / np.where(E_r > EPS, E_r, 1.0), 0.0
    )

    lo = unspl + _flux_density(geom.sj_count_lo, geom.eff_sj_lo)
    hi = unspl + _flux_density(geom.sj_count_hi, geom.eff_sj_hi)
    rho_g = np.where(E_g > EPS, g_cnt / np.where(E_g > EPS, E_g, 1.0), 0.0)
    return unspl, lo, hi, rho_g


def section_counts(r, geom, lab, g_cnt, rna_cnt):
    print("\n── 1. THE COUNTS, every REGION and every BOUNDARY ─────────────────────────────────────────────")
    ra = r.region_arrays
    kind, obj = np.asarray(r.chain.kind), np.asarray(r.chain.obj_idx, np.int64)
    sizes = np.asarray(ra.region_size_bp, np.int64)
    U = np.asarray(geom.unspliced_count, np.float64)
    SP = np.asarray(geom.spliced_count, np.float64)
    JL = np.asarray(geom.sj_count_lo, np.float64)
    JH = np.asarray(geom.sj_count_hi, np.float64)
    E_g = np.asarray(geom.eff_gdna, np.float64)
    E_r = np.asarray(geom.eff_rna, np.float64)
    EJ = np.asarray(geom.eff_sj_lo, np.float64) + np.asarray(geom.eff_sj_hi, np.float64)
    print(f"   {'slot':>4} {'kind':<5} {'what':<26} {'bp':>6} {'unspl+':>8} {'unspl-':>8} "
          f"{'spliced':>8} {'juncLO':>7} {'juncHI':>7} {'E_g':>8} {'E_r':>8} {'E_J':>8}")
    print("   " + "-" * 118)
    for s in range(int(r.chain.n_slots)):
        i = int(obj[s])
        bp = int(sizes[i]) if kind[s] == REGION else 0
        print(f"   {s:>4} {'region' if kind[s] == REGION else 'boundary':<5} {lab[s]:<26} {bp:>6,} "
              f"{U[s, 0]:>8,.0f} {U[s, 1]:>8,.0f} {SP[s].sum():>8,.0f} "
              f"{JL[s].sum():>7,.0f} {JH[s].sum():>7,.0f} "
              f"{E_g[s]:>8,.1f} {E_r[s]:>8,.1f} {EJ[s].sum():>8,.1f}")
    jg = build_sj_geometry_arrays(r.index)
    starts = np.asarray(ra.start, np.int64)
    print()
    for k in range(int(getattr(jg, "n_sj", 0))):
        src, dst = int(np.asarray(jg.src_region)[k]), int(np.asarray(jg.dst_region)[k])
        st = "+" if int(np.asarray(jg.strand)[k]) == int(Strand.POS) else "-"
        flux = float(np.asarray(r.payload.sj_count, np.float64)[k].sum())
        lo = int(starts[src] + sizes[src])
        hi = int(starts[dst])
        print(f"   ⭐ SJ #{k}: {lo:,} → {hi:,}  strand {st}   flux = {flux:,.0f} fragments")
        print(f"      its genomic-LOW end is the BOUNDARY @{lo:,} (index bit DONOR{st.upper()}) — its exon "
              f"is on the LOW side there")
        print(f"      its genomic-HIGH end is the BOUNDARY @{hi:,} (index bit ACCEPTOR{st.upper()}) — its "
              f"exon is on the HIGH side there")
    print("\n   ⭐ TRUTH's origin split, per slot (the oracle BAM's per-fragment read names):")
    print(f"   {'slot':>4} {'what':<26} {'gDNA':>10} {'RNA unspl':>10} {'true f_g':>9}")
    for s in range(int(r.chain.n_slots)):
        tot = g_cnt[s] + rna_cnt[s]
        fg = f"{g_cnt[s] / tot:.4f}" if tot > 0 else "    —"
        print(f"   {s:>4} {lab[s]:<26} {g_cnt[s]:>10,.0f} {rna_cnt[s]:>10,.0f} {fg:>9}")


def section_flank_pair(r, geom, lab, st_cap, t_unspl, t_lo, t_hi):
    print("\n── 2. ⭐⭐ THE FLANK PAIR — which side of each BOUNDARY the sj flux belongs to ──────────")
    print("   ⛔ `rho_lo` is what this slot presents to its genomic-LOW neighbour, `rho_hi` to its")
    print("      genomic-HIGH one. They differ ONLY where sj flux attaches, and the predecessor")
    print("      used ONE number — the sj-INCLUSIVE total — on both sides of every such slot.")
    rho_lo, rho_hi = region_total_density(geom, st_cap["_fg_in"])
    # ⛔ the re-derivation gate: the rebuilt geometry must reproduce the solver's published frames
    np.testing.assert_allclose(rho_lo, np.asarray(st_cap["rho_lo"], float), rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(rho_hi, np.asarray(st_cap["rho_hi"], float), rtol=1e-12, atol=0.0)
    print("   ✅ the rebuilt geometry reproduces the solver's published frames to 1e-12\n")
    # ⚠ back out the shared unspliced part from the geometry's OWN flux banks, never from
    # ``|rho_lo − rho_hi|``: that only works when one bank is empty, and at a slot carrying both it is
    # silently wrong. Here it is exact by construction.
    unspl = rho_lo - _flux_density(geom.sj_count_lo, geom.eff_sj_lo)
    shipped_all = unspl + _flux_density(geom.sj_count, geom.eff_sj)
    print(f"   {'slot':>4} {'what':<26} {'rho_unspl':>10} {'rho_lo':>10} {'rho_hi':>10} "
          f"{'shipped':>10} | {'TRUE lo':>9} {'TRUE hi':>9} | which side got the flux")
    print("   " + "-" * 128)
    for s in range(int(r.chain.n_slots)):
        shipped = shipped_all[s]
        d = ""
        if abs(rho_lo[s] - rho_hi[s]) > EPS:
            d = "⭐ LOW flank INCLUDES it" if rho_lo[s] > rho_hi[s] else "⭐ HIGH flank INCLUDES it"
        print(f"   {s:>4} {lab[s]:<26} {unspl[s]:>10.5f} {rho_lo[s]:>10.5f} {rho_hi[s]:>10.5f} "
              f"{shipped:>10.5f} | {t_lo[s]:>9.5f} {t_hi[s]:>9.5f} | {d}")
    print("\n   ⚠ `shipped` is what the predecessor used at BOTH ends of every hop touching that slot.")
    return rho_lo, rho_hi, unspl, shipped_all


def section_relay(r, direction, rho_lo, rho_hi, shipped_all, t_lo, t_hi, t_rho_g, lab, st_cap, pin):
    """One direction of the relay, hop by hop, with the four ratios side by side."""
    fwd = direction == "forward"
    arrow = "low → high  (L→R genomic)" if fwd else "high → low  (R→L genomic)"
    what = "a SPLICE-OUT at a sj's low end" if fwd else "a SPLICE-IN at a sj's high end"
    n = 3 if fwd else 4
    print(f"\n── {n}. ⭐⭐⭐ THE {'FORWARD' if fwd else 'BACKWARD'} RELAY — {arrow} ─────────────────")
    print(f"   Each hop's message travels {arrow.split('(')[0].strip()}, which is {what}.")
    if fwd:
        print("   ⭐ The source is the destination's genomic-LOW neighbour, so the DESTINATION presents")
        print("      `rho_lo` and the SOURCE presents `rho_hi`.   r = rho_lo[dst] / rho_hi[src]")
    else:
        print("   ⭐ The source is the destination's genomic-HIGH neighbour, so the DESTINATION presents")
        print("      `rho_hi` and the SOURCE presents `rho_lo`.   r = rho_hi[dst] / rho_lo[src]")
    nbr = np.asarray(r.chain.left if fwd else r.chain.right, np.int64)
    lvl = np.asarray(st_cap["fwd_g" if fwd else "bwd_g"], float)
    prec = np.asarray(st_cap["fwd_pg" if fwd else "bwd_pg"], float)
    # ⭐⭐ THE THREE CONJUNCTS OF THE LICENCE, so nobody has to take `lend` on trust. The solver's own
    # expression is `_lend = pop[dst] & (pg[src] > 0) & ((pp + pn)[src] > 0)` — a POPULATION test on the
    # pair, and the SOURCE's own running gDNA and RNA precisions. All three are the relay's own state at
    # the moment the hop is taken (the forward pass writes slot s before it reads it at s < i, so the
    # published arrays ARE the values used).
    pg_s = np.asarray(st_cap["fwd_pg" if fwd else "bwd_pg"], float)
    pp_s = np.asarray(st_cap["fwd_pp" if fwd else "bwd_pp"], float)
    pn_s = np.asarray(st_cap["fwd_pn" if fwd else "bwd_pn"], float)
    lend = np.asarray(pin["lend"], bool)
    r_comb = np.asarray(pin["r"], float)
    num, den = (rho_lo, rho_hi) if fwd else (rho_hi, rho_lo)
    tnum, tden = (t_lo, t_hi) if fwd else (t_hi, t_lo)
    print(f"\n   {'hop':<44} {'r NEW':>8} {'r SHIPPED':>10} {'r TRUE':>8} {'r TRUE_g':>9} "
          f"| {'framed?':>7} {'pg[src]':>9} {'pR[src]':>9} {'pop':>5} {'lend':>5} {'r_g':>7} "
          f"{'level in':>9} {'level out':>10}")
    print("   " + "-" * 172)
    order = range(int(r.chain.n_slots)) if fwd else range(int(r.chain.n_slots) - 1, -1, -1)
    for i in order:
        s = int(nbr[i])
        if s < 0:
            continue
        rn = num[i] / den[s] if (den[s] > EPS and num[i] > EPS) else 1.0
        ro = (
            shipped_all[i] / shipped_all[s]
            if (shipped_all[s] > EPS and shipped_all[i] > EPS)
            else 1.0
        )
        rt = tnum[i] / tden[s] if (tden[s] > EPS and tnum[i] > EPS) else float("nan")
        rtg = t_rho_g[i] / t_rho_g[s] if t_rho_g[s] > EPS else float("nan")
        rg = rn if lend[i] else 1.0
        framed = bool(den[s] > EPS and num[i] > EPS)
        pR = float(pp_s[s] + pn_s[s])
        # pop is not published; it is RECOVERABLE exactly when the other two conjuncts hold
        pop = "True" if lend[i] else ("False" if (pg_s[s] > 0.0 and pR > 0.0) else "—")
        mark = ""
        if abs(ro - rn) > 1e-9:
            # the shipped ratio differs from the new one -> this is a hop the change touches
            e_new = abs(rn / rt - 1.0) if rt == rt and rt > 0 else float("nan")
            e_old = abs(ro / rt - 1.0) if rt == rt and rt > 0 else float("nan")
            mark = "  ⭐ CHANGED" + (
                f", new {e_new:.1%} vs shipped {e_old:.1%} off truth" if e_new == e_new else ""
            )
        print(f"   {lab[s][:20]:<20} → {lab[i][:20]:<20} {rn:>8.4f} {ro:>10.4f} {rt:>8.4f} "
              f"{rtg:>9.4f} | {('YES' if framed else 'NO-FRAME'):>7} {pg_s[s]:>9.3g} {pR:>9.3g} "
              f"{pop:>5} {str(bool(lend[i]))[:5]:>5} {rg:>7.4f} {lvl[s]:>9.5f} {lvl[i]:>10.5f}"
              f"{mark}")
        if abs(float(r_comb[i]) - rn) > 1e-9 and prec[s] >= 0:
            print(f"      ⚠ the combine published r = {float(r_comb[i]):.6f} for this hop; the relay's "
                  f"pairing gives {rn:.6f}")
    print("\n   ⭐ `r TRUE` is the ratio of the SAME two populations, from the oracle's own counts with")
    print("      the same flank rule. `r TRUE_g` is the true gDNA-DENSITY ratio — what a gDNA LEVEL")
    print("      should cross at, and a DIFFERENT number.")
    print("   ⛔ `framed?` = NO-FRAME means one END had ZERO density, so the reframe was SKIPPED and")
    print("      r was forced to 1 — a DEGENERATE hop, not a decision. A row that is NO-FRAME tells you")
    print("      nothing about the reframe logic either way.")
    print("   ⛔ `lend` is the solver's own `pop[dst] & pg[src]>0 & pR[src]>0`, printed beside all three")
    print("      conjuncts so it can be checked rather than trusted. `pop` shows `—` where it is not")
    print("      recoverable (one of the other two already failed).")


def section_decompose(r, geom, lab, g_cnt, rna_cnt):
    """⭐⭐⭐ WHY THE TOTAL-DENSITY RATIO IS NOT THE gDNA-DENSITY RATIO — the exact decomposition.

    ⛔ It is NOT a precision problem, and this section exists to show that with numbers. The identity is
    algebraic and holds to floating point::

        r_tot  =  phi_g(src) . r_g  +  phi_R(src) . r_R

    — the ratio of TOTALS is the share-weighted average of the two components' OWN density ratios,
    weighted by each component's share of the SOURCE's mass. So:

    * at a source that is 100 % gDNA (an intron, an intergenic REGION) ``phi_g = 1`` and ``r_tot`` IS
      ``r_g`` exactly — the scaling is legitimate and carries no foreign component;
    * at a source that is 99.94 % RNA (an expressed exon) ``r_tot`` is 0.9994 . r_R + 0.0006 . r_g, so it
      is a measurement of the RNA ratio and says ⭐ **almost nothing** about the gDNA ratio. Scaling the
      gDNA arm by it imports the RNA's ratio.

    ⭐⭐ **AND THE TWO COMPONENTS' RATIOS FAIL DIFFERENTLY, which is why more counts cannot settle this.**
    The last column is how many Poisson standard errors each ratio sits from **1.0**, which is what every
    ratio here SHOULD be: capture is off, so the true enrichment is 1 and gDNA is uniform along the genome
    while mature RNA is uniform along its transcript. Then:

    * the **gDNA** arm at an BOUNDARY is NOISE-limited — an BOUNDARY is 0 bp, so its opportunity is ~one mean
      fragment length whatever the chromosome does, and it holds tens of counts against a REGION's hundreds.
      Its ratio sits ~1-2 se from 1.0. ⭐ More depth DOES shrink this.
    * the **RNA** arm across a sj is BIAS-limited — it sits >13 se from 1.0 and more depth makes it
      MORE significant, not smaller. That is the sj-vs-exon frame gap: the sj's divisor
      ``E_J = E[w] - 1`` RISES with the fitted mean fragment length while the exon's ``E_r = e - E[w] + 1``
      FALLS, so a length-model error is amplified with opposite sign (`EQUATIONS.md` §3.6b, 0.62 %/bp).
      ⛔ TRAPS: a-variance-cannot-fix-a-bias: a variance cannot fix a biased mode.
    """
    print("\n── 3b. ⭐⭐⭐ WHY r_tot != r_gDNA — the exact decomposition, per hop ─────────────────────")
    print("   r_tot = phi_g(src)·r_g + phi_R(src)·r_R   — a SHARE-WEIGHTED AVERAGE of the two")
    print("   components' own density ratios. `se from 1` is how many Poisson standard errors the ratio")
    print("   sits from 1.0, which is what EVERY ratio here should be with capture off.")
    E_g = np.asarray(geom.eff_gdna, np.float64)
    E_r = np.asarray(geom.eff_rna, np.float64)
    fl_lo = np.asarray(geom.sj_count_lo, np.float64).sum(1)
    fl_hi = np.asarray(geom.sj_count_hi, np.float64).sum(1)
    ej_lo = _flux_density(geom.sj_count_lo, geom.eff_sj_lo)
    ej_hi = _flux_density(geom.sj_count_hi, geom.eff_sj_hi)

    def side(s, which):
        """(rho_g, rho_R, count_g, count_R) at slot s on flank ``which``."""
        rg = g_cnt[s] / E_g[s] if E_g[s] > EPS else 0.0
        rr = (rna_cnt[s] / E_r[s] if E_r[s] > EPS else 0.0) + (ej_lo[s] if which == "lo" else ej_hi[s])
        cr = rna_cnt[s] + (fl_lo[s] if which == "lo" else fl_hi[s])
        return rg, rr, g_cnt[s], cr

    for fwd in (True, False):
        nbr = np.asarray(r.chain.left if fwd else r.chain.right, np.int64)
        dside, sside = ("lo", "hi") if fwd else ("hi", "lo")
        print(f"\n   {'FORWARD' if fwd else 'BACKWARD'} — dst presents its {dside.upper()} flank, "
              f"src its {sside.upper()}")
        print(f"   {'hop':<44} {'phi_g(src)':>10} {'r_g':>9} {'se from 1':>10} {'r_R':>9} "
              f"{'se from 1':>10} {'r_tot':>9} {'= mix?':>9}")
        print("   " + "-" * 128)
        order = range(int(r.chain.n_slots)) if fwd else range(int(r.chain.n_slots) - 1, -1, -1)
        for i in order:
            s = int(nbr[i])
            if s < 0:
                continue
            gd, rd, cgd, crd = side(i, dside)
            gs, rs, cgs, crs = side(s, sside)
            tot_s = gs + rs
            if tot_s <= EPS or (gs <= EPS and rs <= EPS):
                continue
            phi_g = gs / tot_s
            r_g = gd / gs if gs > EPS else float("nan")
            r_R = rd / rs if rs > EPS else float("nan")
            r_tot = (gd + rd) / tot_s
            mix = phi_g * (r_g if r_g == r_g else 0.0) + (1.0 - phi_g) * (r_R if r_R == r_R else 0.0)

            def se(a, b, ratio):
                if not (a > 0 and b > 0 and ratio == ratio and ratio > 0):
                    return float("nan")
                return abs(math.log(ratio)) / math.sqrt(1.0 / a + 1.0 / b)

            se_g = se(cgd, cgs, r_g)
            se_R = se(crd, crs, r_R)

            def fmt(x, w=9, p=4):
                return f"{x:>{w}.{p}f}" if x == x else f"{'—':>{w}}"

            # ⛔⛔ THE ONE CASE THE IDENTITY CANNOT BE WRITTEN AT ALL, and it is the extreme of the
            # whole problem: the SOURCE has ZERO of a component the DESTINATION has plenty of. Then
            # ``r_R`` is undefined (a division by zero), ``phi_R(src)`` is 0, and the product is 0.inf —
            # so ``r_tot`` contains a term with NO counterpart at the source. It is a ratio of densities
            # only dimensionally; as a scale factor it is meaningless. ⭐ This is exactly what the
            # composition licence exists to catch, via its "did the source SUPPLY both components?"
            # conjunct — check the `lend` column in the next section: it is False on both these rows and
            # the gDNA level crosses UNSCALED. Marked NO-r_R rather than as a failure.
            if rs <= EPS or gs <= EPS:
                ok = "NO-r_R" if rs <= EPS else "NO-r_g"
            else:
                ok = "✅" if abs(mix - r_tot) < 1e-9 * max(1.0, abs(r_tot)) else "⛔"
            print(f"   {lab[s][:20]:<20} → {lab[i][:20]:<20} {phi_g:>10.6f} {fmt(r_g)} "
                  f"{fmt(se_g, 10, 1)} {fmt(r_R)} {fmt(se_R, 10, 1)} {fmt(r_tot)} {ok:>9}")
    print("\n   ⭐ Read the `phi_g(src)` column first. Where it is 1.000000 the source is pure gDNA and")
    print("      `r_tot` IS `r_g` — nothing foreign enters. Where it is ~0.0006 the source is an expressed")
    print("      exon and `r_tot` is the RNA ratio wearing a total's name.")
    print("   ⛔ And compare the two `se from 1` columns: the gDNA arm is NOISE (a couple of se, shrinks")
    print("      with depth), the RNA arm across a sj is BIAS (>13 se, GROWS with depth).")
    print("   ⛔⛔ `NO-r_R` marks a hop whose SOURCE has ZERO of a component the DESTINATION has plenty")
    print("      of, so no per-component ratio exists for it and `r_tot` contains a term with no")
    print("      counterpart at the source at all. That is the extreme of this problem and it is the")
    print("      case the composition licence already withholds — `lend` is False on those rows.")


def section_answer(r, lab, g_cnt, rna_cnt):
    print("\n── 5. THE ANSWER, per slot, against per-object truth ────────────────────────────────────")
    # ⚠ the SAME arrays `toy_harness.object_rows` reads. ``_uni['fg_out']`` is the LAST rho-iteration's
    # deconvolution and is NOT the shipped per-slot answer — reading it put 0.0000 at every
    # structurally pure-gDNA slot whose truth is 1.0000, which is how this was caught.
    cap = r.capture
    fg = np.asarray(cap["f_g"], float)
    loc = np.asarray(cap["fg_loc"], float)
    tau = np.asarray(cap["_tau0_lam"], float)
    print(f"   {'slot':>4} {'what':<26} {'n':>9} {'true f_g':>9} {'f_g loc':>9} {'f_g pred':>9} "
          f"{'Δ':>8} {'tau_own':>9} {'err frags':>10}")
    print("   " + "-" * 108)
    tot_err = 0.0
    for s in range(int(r.chain.n_slots)):
        n = g_cnt[s] + rna_cnt[s]
        if n <= 0:
            continue
        true = g_cnt[s] / n
        d = fg[s] - true
        err = abs(d) * n
        tot_err += err
        print(f"   {s:>4} {lab[s]:<26} {n:>9,.0f} {true:>9.4f} {loc[s]:>9.4f} {fg[s]:>9.4f} "
              f"{d:>+8.4f} {tau[s]:>9.3g} {err:>10.1f}")
    print(f"   {'':4} {'TOTAL error mass':<26} {'':>9} {'':>9} {'':>9} {'':>9} {'':>8} {'':>9} "
          f"{tot_err:>10.1f}")


def run_arm(arm, donor_name, spec, index, config, work_dir):
    print("\n" + "=" * 132)
    print(f"════ ARM: {arm.upper()} gDNA   ·   donor {donor_name}   ·   capture OFF · UNSTRANDED")
    print(TH.messages_stamp(bool(config.message_propagation)))
    print("=" * 132)
    donor = TH.harvest(SUITE / donor_name, index, config=config)
    r = TH.run_toy(spec, donor, work_dir, config=config)
    print(f"   RNA fragments {r.spec.n_rna_fragments:,}   ·   gDNA fragments {r.n_gdna_target:,} "
          f"(DERIVED from the donor's {donor.gdna_rate_per_base:.6g}/bp)   ·   "
          f"{r.spec.genome_length:,} bp   ·   {r.seconds:.1f} s")
    print(f"   kappa = {donor.priors.rna_sense_frac:.6f}  (½ ⇒ EXACTLY zero strand information, `EQUATIONS.md` §5.2)")
    geom = rebuild_geometry(r)
    lab = labels(r.chain, r.region_arrays)
    # ⛔⛔ EVERY SECTION BELOW SECTION 1 IS A HOP, and a hop exists only under `HeadPolicy`. The frames
    # (`rho_lo`/`rho_hi`), the relay state (`fwd_*`/`bwd_*`) and the mass pins (`_pin`) are all published
    # by the message layer, so this file died with `KeyError: 'rho_lo'` — the same class of defect as the
    # four `KeyError: '_uni'` instruments, at a different key. Refuse rather than walk an empty relay.
    TH.require_relay(r.capture, what="the reframe walk (every hop, both directions, with its frames)")
    st_cap = TH.relay_static(r.capture)
    # the belief the frames were built at — the sweep's INPUT f_g, which it publishes as fg_init
    st_cap["_fg_in"] = np.asarray(r.capture["fg_init"], float)
    g_cnt, rna_cnt, _full = truth_slot_arrays(r)
    t_unspl, t_lo, t_hi, t_rho_g = truth_flank_pair(r, geom, g_cnt, rna_cnt)
    section_counts(r, geom, lab, g_cnt, rna_cnt)
    rho_lo, rho_hi, unspl, shipped_all = section_flank_pair(r, geom, lab, st_cap, t_unspl, t_lo, t_hi)
    pins = r.capture["_pin"]
    section_decompose(r, geom, lab, g_cnt, rna_cnt)
    section_relay(r, "forward", rho_lo, rho_hi, shipped_all, t_lo, t_hi, t_rho_g, lab, st_cap, pins[0])
    section_relay(r, "backward", rho_lo, rho_hi, shipped_all, t_lo, t_hi, t_rho_g, lab, st_cap, pins[1])
    section_answer(r, lab, g_cnt, rna_cnt)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arms", nargs="*", default=list(ARMS), choices=list(ARMS))
    ap.add_argument("--n-rna", type=int, default=200_000,
                    help="HIGH by default: enough counts everywhere that nothing is sparse")
    ap.add_argument("--spec", default="spliced_exons")
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/rigel_reframe_walk"))
    # ⭐ DEFAULTS TO `on`: this instrument IS the relay — "every hop in BOTH directions". Under the
    # shipped mute there are no hops, so it refuses rather than printing an empty walk.
    TH.add_messages_flag(ap, default=True)
    args = ap.parse_args()

    index = TranscriptIndex.load(str(INDEX))
    config = TH.with_messages(
        dataclasses.replace(CalibrationConfig(), calib_refit_iters=0), TH.messages_on(args)
    )
    spec = dataclasses.replace(TH.SPECS[args.spec], n_rna_fragments=int(args.n_rna))
    print("=" * 132)
    print("⭐⭐⭐ THE REFRAME, WALKED — one two-exon transcript, every count, every hop, both directions")
    print("=" * 132)
    print(f"   spec {spec.name}: " + " ".join(
        f"{t['t_id']}{g['strand']} {t['exons']}" for g in spec.genes for t in g["transcripts"]))
    print("   PRIOR-FREE pass-0 (calib_refit_iters=0) — the substrate the gDNA hyperprior is "
          "later fitted against")
    for arm in args.arms:
        run_arm(arm, ARMS[arm], spec, index, config, args.work_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
