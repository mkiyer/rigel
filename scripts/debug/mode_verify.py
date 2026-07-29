"""Stage 1 — node-level verification of the message-MODE arithmetic in the REAL code.

Runs the production ``calibrate()`` on toy scenarios (via ``toy_prod.run(..., _debug=...)``) and, for every
message edge captured by ``bp_solver._scan`` (the inert ``_capture`` hook), checks:

  (1) CODE == DERIVED — the code's emitted mode equals the derivation's formula for the branch it took
      (shift = normalize by the imputed total; density = ÷ the dst's observed ``md``). Confirms
      ``docs/CARRY_FORWARD.md`` §3/§4 is what the code computes.
  (2) DENSITY vs SHIFT on exon-facing edges — how much the Stage-4/5 predicate flip would move each message,
      expressed as the target gDNA fraction ``exp(mode_g)``, compared against the destination's ORACLE f_g.
  (3) FLIP-SET CLASSIFICATION — every edge → {clean, B-safe, B-src, A, seam} (the plan §0 classes).

This is a read-only probe: the ``_capture`` hook is inert in production (``_capture is None``); goldens unaffected.
Run inside the ``rigel`` env, ``OMP_NUM_THREADS=1``.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

import toy_prod

_EPS = 1.0e-9  # must match bp_solver._EPS
SCRATCH = Path("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/4f7a248b-0c78-4b40-9030-462373aefb19/scratchpad")
toy_prod.SCRATCH = SCRATCH
REGION = 0  # chain.kind


# --------------------------------------------------------------------------- derived mode formulas (match code)
def _shift_gdna(e):
    cont_p = e["fp_s"] and e["fp_d"]
    cont_n = e["fn_s"] and e["fn_d"]
    Mg = e["rho_g"] * e["egd"]
    Mp = (max(e["rho_pos"], 0.0) * e["erd"]) if cont_p else 0.0
    Mn = (max(e["rho_neg"], 0.0) * e["erd"]) if cont_n else 0.0
    den = Mg + Mp + Mn
    den = den if den > _EPS else _EPS
    comp_fl = 1.0 / e["md"]
    return math.log(max(Mg / den, comp_fl))


def _density_gdna(e):
    return math.log(max(e["rho_g"], 1.0 / e["egd"]) / (e["md"] / e["egd"]))


def _classify(e):
    ge = (e["fp_s"] == e["fp_d"]) and (e["fn_s"] == e["fn_d"])
    if e["exon_d"]:  # dst is an exon region ⇒ src is a boundary
        return ("B-safe" if ge else "seam"), ge
    if e["exon_s"]:  # src is an exon region ⇒ dst is a boundary
        return ("B-src" if ge else "seam"), ge
    return ("clean" if ge else "A"), ge  # both non-exon


def verify(name, genes, **kw):
    debug: dict = {}
    rdf, solved, truth, tg, tr, sig = toy_prod.run(name, genes, _debug=debug, **kw)
    edges = debug.get("capture", {}).get("_edge_modes", [])

    # (1) code == derived, over emitted gDNA messages
    max_code_derived = 0.0
    n_checked = 0
    by_class: dict[str, list] = {}
    for e in edges:
        if e["prec_g"] <= 0.0:  # not emitted ⇒ mode_g==0, skip
            continue
        derived = _shift_gdna(e) if e["use_shift"] else _density_gdna(e)
        d = abs(derived - e["mode_g"])
        max_code_derived = max(max_code_derived, d)
        n_checked += 1
        cls, _ = _classify(e)
        by_class.setdefault(cls, []).append(e)

    # (2)+(3) per-class summary; on exon edges compute density-vs-shift target f_g and vs oracle
    print(f"\n===== {name} =====")
    print(f"  edges: {len(edges)} total, {n_checked} with emitted gDNA;  CODE==DERIVED max|Δmode_g| = {max_code_derived:.2e}")
    print(f"  {'class':>8} {'n':>4} | {'exp(mode) density→shift (target f_g)':>38} | {'oracle f_g (dst region)':>24}")
    for cls in ("clean", "A", "B-safe", "B-src", "seam"):
        es = by_class.get(cls, [])
        if not es:
            continue
        d_tf, s_tf, orac = [], [], []
        for e in es:
            d_tf.append(math.exp(_density_gdna(e)))
            s_tf.append(math.exp(_shift_gdna(e)))
            if e["kind_d"] == REGION and not np.isnan(truth[e["ref_d"]]):
                orac.append(truth[e["ref_d"]])
        d_m = np.mean(d_tf); s_m = np.mean(s_tf)
        o_m = np.mean(orac) if orac else float("nan")
        flip = "→shift" if cls in ("B-safe", "B-src") else ("→density" if cls == "A" else "unchanged")
        print(f"  {cls:>8} {len(es):>4} | density {d_m:6.3f}  shift {s_m:6.3f}  ({flip:>9}) | {o_m:>10.3f}"
              f"   {'closer:shift' if abs(s_m-o_m) < abs(d_m-o_m) else 'closer:density' if not np.isnan(o_m) else ''}")
    return max_code_derived


def _emits(e):
    return (e["prec_g"] > 0.0) or (e["prec_p"] > 0.0) or (e["prec_n"] > 0.0)


def flip_set_report(scenarios):
    """Stage 3 — the blast radius of ``use_shift = gates_equal``. For every DIRECTED edge, classify into
    {clean, seam (unchanged) | B-safe, B-src, A (CHANGE)} and count total vs emitting. Cross-checks that the
    3 changing classes are EXACTLY the edges where ``use_shift != gates_equal`` (self-consistency)."""
    CHANGING = ("B-safe", "B-src", "A")
    ALL = ("clean", "seam", "B-safe", "B-src", "A")
    agg = {c: [0, 0] for c in ALL}  # [total, emitting]
    print("\n" + "=" * 96)
    print("STAGE 3 — FLIP-SET DIAGNOSIS  (change = use_shift != gates_equal;  'emit' = sends any message)")
    print(f"  {'scenario':>22} | " + " | ".join(f"{c:>12}" for c in CHANGING) + " | clean seam | chk")
    for name, kw in scenarios:
        debug: dict = {}
        toy_prod.run(name, kw.pop("genes"), _debug=debug, **kw)
        edges = debug.get("capture", {}).get("_edge_modes", [])
        per = {c: [0, 0] for c in ALL}
        mismatch = 0
        for e in edges:
            cls, ge = _classify(e)
            per[cls][0] += 1
            per[cls][1] += int(_emits(e))
            # self-consistency: changing ⟺ (use_shift != gates_equal)
            if (cls in CHANGING) != (bool(e["use_shift"]) != ge):
                mismatch += 1
        for c in ALL:
            agg[c][0] += per[c][0]; agg[c][1] += per[c][1]
        cells = " | ".join(f"{per[c][0]:>4} ({per[c][1]:>3})" for c in CHANGING)
        chk = "OK" if mismatch == 0 else f"FAIL×{mismatch}"
        print(f"  {name:>22} | {cells} | {per['clean'][0]:>4} {per['seam'][0]:>4} | {chk}")
    print("  " + "-" * 94)
    cells = " | ".join(f"{agg[c][0]:>4} ({agg[c][1]:>3})" for c in CHANGING)
    print(f"  {'AGGREGATE (tot / emit)':>22} | {cells} | {agg['clean'][0]:>4} {agg['seam'][0]:>4} |")
    print("\n  Reading: 'N (M)' = N edges in the class, M of which emit a message. Only EMITTING edges feel the")
    print("  mode flip. B-safe = Stage 4 (boundary→exon), B-src = Stage 5 (exon→boundary), A = Stage 6 (gate-")
    print("  unequal non-exon → density). clean/seam are UNCHANGED by the flip.")


def main():
    TA = [("TA", "+", [(1000, 2000), (5000, 10000)])]
    # two genes, opposite strands, with an intergenic gap ⇒ exercises strand-transition (class-A) edges
    TWO = [("TA", "+", [(1000, 2000), (5000, 9000)]), ("TB", "-", [(14000, 16000), (19000, 23000)])]
    m = 0.0
    m = max(m, verify("verify_S4_unstr_cap", TA, kappa=0.5, n_rna=4000, gdna_fraction=0.5, capture=True, capture_strength=20.0))
    m = max(m, verify("verify_S2_unstr", TA, kappa=0.5, n_rna=4000, gdna_fraction=0.5, capture=False))
    m = max(m, verify("verify_S3_str_cap", TA, kappa=1.0, n_rna=4000, gdna_fraction=0.5, capture=True, capture_strength=20.0))
    print(f"\n=== OVERALL CODE==DERIVED max|Δmode_g| across scenarios: {m:.2e} "
          f"({'PASS — code implements the derivation' if m < 1e-9 else 'INVESTIGATE'}) ===")

    flip_set_report([
        ("1gene str  cap", dict(genes=[g[:] for g in TA], kappa=1.0, capture=True, capture_strength=20.0)),
        ("1gene unstr cap", dict(genes=[g[:] for g in TA], kappa=0.5, capture=True, capture_strength=20.0)),
        ("1gene unstr off", dict(genes=[g[:] for g in TA], kappa=0.5, capture=False)),
        ("2gene str  off", dict(genes=[g[:] for g in TWO], kappa=1.0, capture=False, genome_length=25000)),
        ("2gene unstr cap", dict(genes=[g[:] for g in TWO], kappa=0.5, capture=True, capture_strength=20.0, genome_length=25000)),
    ])


if __name__ == "__main__":
    main()
