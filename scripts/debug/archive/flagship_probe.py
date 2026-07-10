"""Flagship single-strand-node accuracy: capture-ON, strand-specific (ss0.99), ±gDNA.

The measure of the Phase-1 single-strand solver (must work WITHOUT a driving prior). Per single-strand
REGION node: oracle f_g (from the BAM, contained gDNA/RNA via mate+TLEN) vs solved f_g (calibrate's
contained masses). Reports median |err| + signed bias for single-strand EXONS and INTRONS separately
(AMBIG excluded — that's Phase 3). Partial summary controls for the gDNA/expression co-location trap by
splitting exons into expressed (RNA>0) vs unexpressed.

Toggle the solver config by editing bp_solver (weak floor / ê(z) off) and re-running.
"""
from __future__ import annotations
import sys
import collections
from pathlib import Path
import numpy as np
import pysam

from rigel.index import TranscriptIndex
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.pipeline import scan_and_buffer
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.calibrate import calibrate
from rigel.calibration.signature import (
    BIT_EXON_POS, BIT_EXON_NEG, BIT_INTRON_POS, BIT_INTRON_NEG, TS_POS, TS_NEG,
)
from rigel.splice import SpliceType

SUITE = Path.home() / "Downloads/rigel_runs/quick_3to1_5mb"
IDX = SUITE / "rigel_index"


def _oracle_fg(bam_path, ra, ref_names):
    """Per-region oracle (gDNA, RNA) contained counts via read1 mate+TLEN span fully inside the region."""
    starts = ra.start; ends = ra.end; rid = ra.ref_id; R = starts.shape[0]
    by_ref = collections.defaultdict(list)
    for i in range(R):
        by_ref[int(rid[i])].append(i)
    for k in list(by_ref):
        arr = sorted(by_ref[k], key=lambda i: starts[i])
        by_ref[k] = (np.array([starts[i] for i in arr]), np.array([ends[i] for i in arr]),
                     np.array(arr))
    name2ref = {n: i for i, n in enumerate(ref_names)}
    g = np.zeros(R); r = np.zeros(R)
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for rd in bam.fetch(until_eof=True):
        if (rd.is_secondary or rd.is_supplementary or rd.is_unmapped or not rd.is_read1
                or rd.mate_is_unmapped):
            continue
        is_g = "gdna" in rd.query_name.lower()
        rf = name2ref.get(bam.get_reference_name(rd.reference_id))
        if rf is None or rf not in by_ref:
            continue
        tl = abs(rd.template_length)
        if tl == 0:
            continue
        lo = min(rd.reference_start, rd.next_reference_start); hi = lo + tl
        s, e, ii = by_ref[rf]; j = np.searchsorted(s, lo, side="right") - 1
        if 0 <= j < R and s[j] <= lo and hi <= e[j]:
            (g if is_g else r)[ii[j]] += 1.0
    bam.close()
    return g, r


def probe(cond, label):
    bam = SUITE / cond / "sim_oracle.bam"
    idx = TranscriptIndex.load(IDX)
    _st, sm, flm, _buf, pl = scan_and_buffer(str(bam), idx, BamScanConfig())
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    cal = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, CalibrationConfig())
    gc = np.asarray(cal.mass_gdna_contained); rc = np.asarray(cal.mass_rna_contained)
    solved = np.where((gc + rc) > 0, gc / np.maximum(gc + rc, 1e-9), np.nan)
    g, r = _oracle_fg(bam, ra, idx.ref_names)
    oracle = np.where((g + r) > 0, g / np.maximum(g + r, 1e-9), np.nan)
    sig = np.asarray(ra.signature).astype(int)
    sc = np.asarray(ra.strand_class)
    ss = (sc == TS_POS) | (sc == TS_NEG)  # single-strand
    is_exon = (sig & (BIT_EXON_POS | BIT_EXON_NEG)) != 0
    is_intron = (sig & (BIT_INTRON_POS | BIT_INTRON_NEG)) != 0
    obs = (g + r) >= 5  # only score regions with enough oracle fragments

    def rep(mask, name):
        m = mask & ss & obs & np.isfinite(solved) & np.isfinite(oracle)
        if m.sum() == 0:
            print(f"    {name:28s} n=0"); return
        err = solved[m] - oracle[m]
        print(f"    {name:28s} n={int(m.sum()):4d}  median|err|={np.median(np.abs(err)):.3f}  "
              f"bias={np.mean(err):+.3f}  oracle_fg med={np.median(oracle[m]):.3f}")

    print(f"\n=== {label}: {cond} ===")
    rep(is_exon, "single-strand EXON (all)")
    rep(is_exon & (r > 0) & (g > 0), "  exon expressed+gDNA")
    rep(is_exon & (r == 0), "  exon UNexpressed (pure gDNA)")
    rep(is_intron, "single-strand INTRON (all)")


def dive(cond):
    """Decompose the pure-gDNA single-strand under-call: strand-only vs local(prior) vs final, + the
    message that drags them. Uses the node_sweep capture hook."""
    bam = SUITE / cond / "sim_oracle.bam"
    idx = TranscriptIndex.load(IDX)
    _st, sm, flm, _buf, pl = scan_and_buffer(str(bam), idx, BamScanConfig())
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    cap = {}
    cm = sys.modules["rigel.calibration.calibrate"]
    orig = cm.node_sweep
    cm.node_sweep = lambda *a, **k: orig(*a, _capture=cap, **k)
    try:
        cm.calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, CalibrationConfig())
    finally:
        cm.node_sweep = orig
    from rigel.calibration.node_chain import build_node_chain, REGION
    ch = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    kind = np.asarray(ch.kind); ridx = np.asarray(ch.ref_idx)
    reg_node = {int(ridx[n]): n for n in range(ch.n_nodes) if kind[n] == REGION}
    g, r = _oracle_fg(bam, ra, idx.ref_names)
    sig = np.asarray(ra.signature).astype(int); sc = np.asarray(ra.strand_class)
    ss = (sc == TS_POS) | (sc == TS_NEG)
    is_exon = (sig & (BIT_EXON_POS | BIT_EXON_NEG)) != 0
    is_intron = (sig & (BIT_INTRON_POS | BIT_INTRON_NEG)) != 0
    fgs = np.asarray(cap["fg_strand"]); fgl = np.asarray(cap["fg_loc"]); fgf = np.asarray(cap["f_g"])
    pg = np.asarray(cap["prec_g"]); pp = np.asarray(cap["prec_p"]); pn = np.asarray(cap["prec_n"])
    mp = np.asarray(cap["mode_p"]); mn = np.asarray(cap["mode_n"])

    def cls(mask, name):
        regs = [i for i in np.where(mask & ss & ((g + r) >= 20))[0] if i in reg_node]
        if not regs:
            print(f"  {name:26s} n=0"); return
        nodes = [reg_node[i] for i in regs]
        orc = np.array([g[i] / max(g[i] + r[i], 1e-9) for i in regs])
        print(f"  {name:26s} n={len(regs):3d} | oracle_fg={np.median(orc):.3f} | "
              f"fg_strand={np.median(fgs[nodes]):.3f} fg_loc={np.median(fgl[nodes]):.3f} "
              f"fg_final={np.median(fgf[nodes]):.3f} | gMSGpr={np.median(pg[nodes]):.2f} "
              f"pMSGpr={np.median(pp[nodes]+pn[nodes]):.2f} pMSGmode={np.median(mp[nodes]+mn[nodes]):+.2f}")

    print(f"\n=== DIVE {cond}: strand vs local vs final (pure-gDNA nodes; fg should be ~1.0) ===")
    cls(is_exon & (r == 0), "exon UNexpressed")
    cls(is_exon & (r > 0) & (g > 0), "exon expressed+gDNA")
    cls(is_intron, "intron (gDNA)")


if __name__ == "__main__":
    args = sys.argv[1:]
    if args and args[0] == "dive":
        dive("gdna_gdna300_ss_0.99_nrna_none_capture_on")
    else:
        conds = [("gdna_gdna300_ss_0.99_nrna_none_capture_on", "FLAGSHIP +gDNA (cap,ss0.99)"),
                 ("gdna_none_ss_0.99_nrna_none_capture_on", "FLAGSHIP −gDNA (cap,ss0.99)")]
        for c, lab in conds:
            probe(c, lab)
