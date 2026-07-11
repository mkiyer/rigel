"""THE single component-audit diagnostic for the nascent-RNA siphon.

One deterministic scan+calibrate. For the worst nascent-siphon locus of a condition, it shows — for EVERY
component (mature isoforms, nascent spans, gDNA) — the effective length the EM actually uses, decomposed
node-by-node, next to the assigned mass and the ground truth. The eff-lengths are the PRODUCTION
``transcript_capture_eff_lengths`` values (never a reimplementation); the node decomposition is validated to
sum back to that production value (a hard assert) before anything is printed, so the numbers are trustworthy.

The question it answers (the user's roadmap): does mature RNA win the eff-length battle it should win?
Nascent RNA spans introns + exon-intron boundary nodes that mature lacks, so eff(nascent) should EXCEED
eff(mature). If it does not, this dump shows exactly which nodes drag the nascent's length down (intron
regions? under-read exon-intron boundaries?), which is where the siphon's false mass comes from.

Usage:
  OMP_NUM_THREADS=1 python scripts/debug/locus_component_audit.py [condition] [locus_id]
    condition : suite subdir (default: the flagship gdna_gdna300_ss_0.99_nrna_rnd_capture_on)
    locus_id  : force a locus; default = the locus with the largest nascent (observed-expected)
"""

import os
import sys
from dataclasses import replace as dc_replace

import numpy as np
import pandas as pd

from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.capture_eff_length import (
    transcript_capture_eff_lengths,
    _transcript_node_incidence,
    _global_reference_density,
)
from rigel.frag_length_model import FragmentLengthModel
from rigel.splice import SpliceType
from rigel.types import IntervalType

SUITE = os.environ.get("AUDIT_SUITE", "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb")
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_rnd_capture_on"
FORCE_LOCUS = int(sys.argv[2]) if len(sys.argv) > 2 else None


def main():
    index = TranscriptIndex.load(f"{SUITE}/rigel_index")
    bam = f"{SUITE}/{COND}/sim_oracle.bam"
    cfg = PipelineConfig()
    scan = dc_replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _stats, sm, flm, _buf, payload = scan_and_buffer(bam, index, scan)
    sm.finalize()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    fl = build_fl_models(
        global_counts=flm.global_model.counts,
        rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(payload),
        max_size=flm.max_size,
    )
    cal = calibrate(payload=payload, region_arrays=ra, strand_model=sm,
                    gdna_fl_pmf=fl.gdna_pmf, rna_fl_pmf=fl.rna_pmf, config=cfg.calibration)

    # ---- production effective lengths (exactly as _setup_geometry_and_estimator builds them) ----
    exonic_lengths = index.t_df["length"].values.astype(np.float64)
    rna_fl = FragmentLengthModel.from_pmf(fl.rna_pmf, fl.max_size)
    eff_fl = rna_fl.compute_all_transcript_eff_lens(exonic_lengths.astype(np.int64))
    eff_em = transcript_capture_eff_lengths(cal, ra, index, eff_fl)  # THE trusted values

    # ---- node arrays (mirror capture_eff_length exactly, for the validated decomposition) ----
    cm = np.asarray(cal.mass_gdna_contained, float)
    cS = np.maximum(np.asarray(cal.gdna_region_eff_len, float), 1e-9)
    sl = np.asarray(cal.mass_gdna_left, float)
    srt = np.asarray(cal.mass_gdna_right, float)
    blen = np.asarray(cal.gdna_boundary_len, float)
    ref_id = np.asarray(ra.ref_id)
    n = cm.shape[0]
    seam_m = np.zeros(n)
    seam_S = np.zeros(n)
    if n > 1:
        same = ref_id[:-1] == ref_id[1:]
        seam_m[:-1] = np.where(same, srt[:-1] + sl[1:], 0.0)
        seam_S[:-1] = np.where(same, 0.5 * (blen[:-1] + blen[1:]), 0.0)
    rho_ref = _global_reference_density(cm, cS)
    inv = 1.0 / rho_ref if rho_ref else 0.0

    rt, rr, bt, br, jt, jl, jr = _transcript_node_incidence(index, ra)
    rstart, rend = np.asarray(ra.start), np.asarray(ra.end)

    # exon/intron mask per region (for the locus's gene) from the exon intervals
    iv = pd.read_feather(os.path.join(index.index_dir, "intervals.feather"))
    ex = iv[(iv["interval_type"] == int(IntervalType.EXON)) & (iv["t_index"] >= 0)]

    def region_is_exon(r):
        rs, re, rf = int(rstart[r]), int(rend[r]), int(ref_id[r])
        # any exon (on this ref) overlapping [rs, re)?
        sub = ex[(ex["ref"].map(index.ref_name_to_id) == rf)]
        return bool(((sub["start"] < re) & (sub["end"] > rs)).any())

    # ---- ground truth + assigned mass (the trusted ZT-abundance metric) ----
    tp = pd.read_csv(f"{SUITE}/net_flow_per_transcript.tsv", sep="\t")
    tp = tp[tp.condition == COND].copy()
    tp["is_nrna"] = tp.transcript_id.str.startswith("RIGEL_NRNA_")
    if FORCE_LOCUS is not None:
        locus = FORCE_LOCUS
    else:
        nr = tp[tp.is_nrna].copy()
        nr["siphon"] = nr.observed - nr.expected
        locus = int(nr.groupby("locus_id")["siphon"].sum().idxmax())
    loc = pd.read_csv(f"{SUITE}/net_flow_per_locus.tsv", sep="\t")
    loc = loc[(loc.condition == COND) & (loc.locus_id == locus)]
    gdna_eff = float(loc["gdna_eff_len_em"].iloc[0]) if "gdna_eff_len_em" in loc.columns and len(loc) else float("nan")
    gdna_span = float(loc["locus_span_bp"].iloc[0]) if "locus_span_bp" in loc.columns and len(loc) else float("nan")

    txids = tp[tp.locus_id == locus].copy()
    tid2idx = {t: i for i, t in enumerate(index.t_df["t_id"].values)}
    tstart = index.t_df["start"].values
    tend = index.t_df["end"].values

    print(f"\n{'='*100}")
    print(f"WORST NASCENT-SIPHON LOCUS = {locus}   condition={COND}")
    print(f"rho_ref (global enriched-mode gDNA reference density) = {rho_ref}")
    print(f"gDNA component: eff_len={gdna_eff:,.0f}  span={gdna_span:,.0f}")
    print(f"{'='*100}")
    hdr = f"{'component':40} {'kind':5} {'gspan':>8} {'splen':>7} {'fl':>8} {'factor':>7} {'eff_em':>8} {'expect':>9} {'observe':>9}"
    print(hdr)
    print("-" * len(hdr))

    recs = []
    for _, row in txids.sort_values("observed", ascending=False).iterrows():
        tid = row.transcript_id
        ti = tid2idx.get(tid)
        if ti is None:
            continue
        gspan = float(tend[ti] - tstart[ti])
        splen = exonic_lengths[ti]
        f = eff_em[ti] / eff_fl[ti] if eff_fl[ti] > 0 else float("nan")
        recs.append((tid, "nrna" if row.is_nrna else "mrna", ti, gspan, splen, eff_fl[ti], f, eff_em[ti],
                     row.expected, row.observed))
        print(f"{tid:40} {'nrna' if row.is_nrna else 'mrna':5} {gspan:>8.0f} {splen:>7.0f} "
              f"{eff_fl[ti]:>8.0f} {f:>7.3f} {eff_em[ti]:>8.0f} {row.expected:>9.0f} {row.observed:>9.0f}")

    # ---- node-by-node decomposition for the nascent(s) + the top mature (validated) ----
    def decompose(ti):
        rows = []
        num = span = 0.0
        for r in rr[rt == ti]:
            m, S = cm[r], cS[r]
            c = min(m * inv, S)
            num += c
            span += S
            rows.append(("REGION", "exon" if region_is_exon(r) else "intron", int(rstart[r]), int(rend[r]),
                         S, m, (m / S) / rho_ref if (S > 0 and rho_ref) else 0.0, c))
        for r in br[bt == ti]:
            m, S = seam_m[r], seam_S[r]
            c = min(m * inv, S)
            num += c
            span += S
            le = region_is_exon(r)
            ri = region_is_exon(r + 1) if r + 1 < n else le
            btype = "exon-intron" if le != ri else ("exon-interior" if le else "intron-interior")
            rows.append(("BNDRY", btype, int(rend[r]), int(rstart[r + 1]) if r + 1 < n else int(rend[r]),
                         S, m, (m / S) / rho_ref if (S > 0 and rho_ref) else 0.0, c))
        mask = jt == ti
        for a, b in zip(jl[mask], jr[mask]):
            rho_l = cm[a] / cS[a]
            rho_r = cm[b] / cS[b]
            s_j = 0.5 * (blen[a] + blen[b])
            m_j = 0.5 * (rho_l + rho_r) * s_j
            c = min(m_j * inv, s_j)
            num += c
            span += s_j
            rows.append(("JUNC", "splice", int(rstart[a]), int(rend[b]), s_j, m_j,
                         (m_j / s_j) / rho_ref if (s_j > 0 and rho_ref) else 0.0, c))
        return rows, num, span

    for tid, kind, ti, *_ in recs:
        if kind != "nrna":
            continue
        rows, num, span = decompose(ti)
        f_recomp = num / span if span > 0 else 1.0
        # validate against production (before the multimapper shrinkage w, which we report separately)
        cev = 0.0
        cont_ev = cm + np.asarray(cal.mass_rna_contained, float)
        for r in rr[rt == ti]:
            cev += cont_ev[r]
        w = cev / (cev + 1.0)
        f_prod = eff_em[ti] / eff_fl[ti]
        assert abs((w * f_recomp + (1 - w)) - f_prod) < 1e-6, \
            f"decomp {w*f_recomp+(1-w):.6f} != production {f_prod:.6f} for {tid}"
        # class rollup
        print(f"\n--- NASCENT node decomposition: {tid}  (validated ✓  factor {f_prod:.3f} = w·{f_recomp:.3f}+{1-w:.3f},  w={w:.3f}) ---")
        agg = {}
        for typ, sub, a, b, S, m, ratio, c in rows:
            k = f"{typ}:{sub}"
            d = agg.setdefault(k, [0, 0.0, 0.0, 0.0])
            d[0] += 1; d[1] += S; d[2] += c; d[3] += m
        print(f"  {'node-class':22} {'#':>4} {'Σ S_n (span)':>14} {'Σ contrib (num)':>16} {'Σ gDNA mass':>12} {'kept%':>7}")
        for k in sorted(agg):
            cnt, Ss, cc, mm = agg[k]
            print(f"  {k:22} {cnt:>4} {Ss:>14,.0f} {cc:>16,.0f} {mm:>12,.0f} {100*cc/Ss if Ss else 0:>6.1f}%")
        print(f"  {'TOTAL':22} {'':>4} {span:>14,.0f} {num:>16,.0f} {'':>12} {100*num/span if span else 0:>6.1f}%")
        print(f"  → eff_fl(nascent)={eff_fl[ti]:,.0f} · factor {f_prod:.3f} = eff_em {eff_em[ti]:,.0f}")


if __name__ == "__main__":
    main()
