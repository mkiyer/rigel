"""The QUINTUPLE grid (owner interrogation 2) — the pass-0 nascent≫mature over-call harness.

Geometry: a 3-exon gene whose MIDDLE exon is flanked by two introns —
    near-intron (R) — IE boundary (B) — exon (R) — EI boundary (B) — far-intron (R)
The exon is a RELAY, not a source: at pass-0 it knows nothing about its own gDNA; that must be imputed from
the introns (the 3-hop game of telephone). So a signal that the exon's unspliced mass is NASCENT (not gDNA)
can only arrive from the introns via the boundaries.

This instruments the FULL quintuple (all 5 nodes, not just the exon) across a
    capture-strength × gDNA-fraction × nascent × mature
sweep and reports, per node: oracle f_g, solved f_g, the belief variance var_g (confidence — is a wrong
solve CONFIDENT or honest-weak?), the incoming gDNA message precision + mode (the relay strength/target),
and the strand-only / local (strand+prior, no-message) fractions.

PART A — the exon over-call MAP: (solved−oracle) f_g on a nascent×mature grid per (capture, gDNA).
PART B — the relay DISSECTION at the nascent≫mature EDGE: the 5-node table, answering
    (1) do the introns resolve as RNA (low solved f_g) but fail to RELAY it to the exon?
    (2) is the exon over-call CONFIDENT (low var_g — bad) or honest-weak (high var_g — acceptable)?

Diagnostic only (no production knobs). Run in the `rigel` env, OMP_NUM_THREADS=1.
"""
from __future__ import annotations
import contextlib
import io
import sys
import numpy as np
import pysam

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
import toy_prod
from rigel.calibration.node_chain import REGION

# 3-exon gene: exons 1000-2000, 4000-5000, 7000-8000; introns 2000-4000 and 5000-7000 (genome 10000).
# The MIDDLE exon (4000-5000) is the quintuple centre.
EXONS = [(1000, 2000), (4000, 5000), (7000, 8000)]
_MID_START = 4000

# sweep axes.  capture: (label, on, binding_per_base).  gDNA fraction; nascent (nrna_abundance); mature (abundance).
CAPTURE = [("off", False, 0.0), ("cap20", True, 20.0), ("cap100", True, 100.0)]
GDNA = [0.05, 0.5, 0.9]
NASC = [20, 200, 400]
MAT = [20, 200, 400]
_EDGE = (400, 20)  # the nascent≫mature corner Part B dissects
_EPS = 1e-9


def _crossing_fg(bam_path, coords):
    """Oracle boundary f_g = contiguously-crossing gDNA / (gDNA + RNA) at each genomic coord. Contiguous
    (per get_blocks) EXCLUDES splice-gap crossings, so this is the mature-free crossing composition — exactly
    what a boundary node sees (gDNA + nascent)."""
    g = {c: 0.0 for c in coords}
    r = {c: 0.0 for c in coords}
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for rd in bam.fetch(until_eof=True):
        if rd.is_secondary or rd.is_supplementary or rd.is_unmapped:
            continue
        isg = "gdna" in rd.query_name.lower()
        for bs, be in rd.get_blocks():
            for c in coords:
                if bs < c < be:
                    (g if isg else r)[c] += 1.0
    bam.close()
    return {c: (g[c] / (g[c] + r[c]) if (g[c] + r[c]) > 0 else np.nan) for c in coords}


def _quintuple_nodes(chain, mid_reg):
    """Return [near_intron, IE_bnd, exon, EI_bnd, far_intron] chain-node ids for the middle exon region."""
    kind = np.asarray(chain.kind)
    ref_idx = np.asarray(chain.ref_idx)
    left = np.asarray(chain.left)
    right = np.asarray(chain.right)
    reg2node = {int(ref_idx[n]): n for n in range(chain.n_nodes) if kind[n] == REGION}
    exon = reg2node[mid_reg]
    ie = int(left[exon])
    near = int(left[ie])
    ei = int(right[exon])
    far = int(right[ei])
    return [near, ie, exon, ei, far]


def _run_cell(cap_lab, cap_on, cap_str, gf, nasc, mat):
    """One grid cell → per-quintuple-node records. Returns (labels, dict-of-arrays) or (None, err)."""
    dbg = {}
    name = f"q_{cap_lab}_{gf}_{nasc}_{mat}"
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            rdf, solved, truth, _tg, _tr, _sig = toy_prod.run(
                name, [("G", "+", EXONS, float(mat))], kappa=0.5, n_rna=4000, genome_length=10000,
                gdna_fraction=gf, nascent=float(nasc), capture=cap_on, capture_strength=cap_str,
                seed=7, _debug=dbg)
    except Exception as e:  # noqa: BLE001 — diagnostic: a corner may hit a sim/calibrate guard
        return None, f"{type(e).__name__}: {e}"
    cap = dbg["capture"]
    chain = dbg["chain"]
    starts = rdf["start"].to_numpy()
    hits = np.where(starts == _MID_START)[0]
    if hits.size == 0:
        return None, "no-mid-exon"
    mid_reg = int(hits[0])
    quint = _quintuple_nodes(chain, mid_reg)
    ref_idx = np.asarray(chain.ref_idx)
    kind = np.asarray(chain.kind)
    fg = np.asarray(cap["f_g"])
    vg = np.asarray(cap["var_g"])
    pg = np.asarray(cap["prec_g"])
    mg = np.asarray(cap["mode_g"])
    pp = np.asarray(cap["prec_p"])
    pn = np.asarray(cap["prec_n"])
    fgs = np.asarray(cap["fg_strand"])
    fgl = np.asarray(cap["fg_loc"])
    # boundary oracle at the two crossing coords
    en = int(rdf["end"].to_numpy()[mid_reg])
    cross = _crossing_fg(dbg["bam_path"], [_MID_START, en])
    rec = {k: [] for k in ("oracle", "solved", "var_g", "prec_g", "mode_g", "prec_r", "fg_str", "fg_loc")}
    for n in quint:
        if kind[n] == REGION:
            orac = truth[int(ref_idx[n])]
        else:
            orac = cross[_MID_START] if n == quint[1] else cross[en]
        rec["oracle"].append(float(orac))
        rec["solved"].append(float(fg[n]))
        rec["var_g"].append(float(vg[n]))
        rec["prec_g"].append(float(pg[n]))
        rec["mode_g"].append(float(mg[n]))
        rec["prec_r"].append(float(pp[n] + pn[n]))
        rec["fg_str"].append(float(fgs[n]))
        rec["fg_loc"].append(float(fgl[n]))
    return ["near-intron", "IE-bnd", "EXON", "EI-bnd", "far-intron"], {k: np.array(v) for k, v in rec.items()}


def main():
    print("QUINTUPLE grid — capture × gDNA × nascent × mature; near-intron—IE—EXON—EI—far-intron\n")
    cache = {}
    errs = []
    for cap_lab, cap_on, cap_str in CAPTURE:
        for gf in GDNA:
            for na in NASC:
                for ma in MAT:
                    labels, res = _run_cell(cap_lab, cap_on, cap_str, gf, na, ma)
                    if labels is None:
                        errs.append((cap_lab, gf, na, ma, res))
                    else:
                        cache[(cap_lab, gf, na, ma)] = res
    if errs:
        print(f"[{len(errs)} cells failed] e.g. {errs[0]}\n")

    # ---- PART A: exon (solved − oracle) map, per (capture, gDNA), nascent(rows) × mature(cols) ----
    print("=" * 78)
    print("PART A — EXON over-call:  (solved − oracle) f_g   [+ = over-calls gDNA]")
    print("         rows = nascent, cols = mature; the top-right (nascent≫mature) is the edge\n")
    for cap_lab, _o, _s in CAPTURE:
        for gf in GDNA:
            print(f"  capture={cap_lab:>6}  gDNA_frac={gf:<5}   mat→ {'  '.join(f'{m:>6}' for m in MAT)}")
            for na in NASC:
                cells = []
                for ma in MAT:
                    res = cache.get((cap_lab, gf, na, ma))
                    if res is None:
                        cells.append("   nan")
                        continue
                    d = res["solved"][2] - res["oracle"][2]  # EXON = index 2
                    cells.append(f"{d:>+6.2f}")
                print(f"    nasc={na:<5}                {'  '.join(cells)}")
            print()

    # ---- PART B: the relay dissection at the nascent≫mature EDGE ----
    na_e, ma_e = _EDGE
    print("=" * 78)
    print(f"PART B — relay DISSECTION at the EDGE  nascent={na_e} ≫ mature={ma_e}")
    print("  Q1: do the introns solve LOW f_g (resolve as RNA) but the EXON stay HIGH (relay fails)?")
    print("  Q2: is a wrong EXON solve CONFIDENT (low var_g) or honest-weak (high var_g)?")
    print("  (mode_g = incoming gDNA message log-fraction target; fg_str = strand-only; fg_loc = +prior,no-msg)\n")
    for cap_lab, _o, _s in CAPTURE:
        for gf in GDNA:
            res = cache.get((cap_lab, gf, na_e, ma_e))
            if res is None:
                continue
            labels = ["near-intron", "IE-bnd", "EXON", "EI-bnd", "far-intron"]
            print(f"  --- capture={cap_lab:>6}  gDNA_frac={gf} ---")
            print(f"    {'node':>11} | {'orac':>6} {'solv':>6} {'err':>6} | {'fg_str':>6} {'fg_loc':>6} | "
                  f"{'var_g':>6} {'precG':>6} {'modeG':>6} {'precR':>6}")
            for i, lab in enumerate(labels):
                o, s = res["oracle"][i], res["solved"][i]
                err = "" if np.isnan(o) else f"{s - o:+.2f}"
                os_ = "  nan" if np.isnan(o) else f"{o:.3f}"
                print(f"    {lab:>11} | {os_:>6} {s:>6.3f} {err:>6} | {res['fg_str'][i]:>6.3f} "
                      f"{res['fg_loc'][i]:>6.3f} | {res['var_g'][i]:>6.3f} {res['prec_g'][i]:>6.2f} "
                      f"{res['mode_g'][i]:>+6.2f} {res['prec_r'][i]:>6.2f}")
            print()


if __name__ == "__main__":
    main()
