"""Production-faithful toy driver: full simulator → oracle BAM → REAL calibrate().

No hand-rolled harness — this IS the production calibration path, on tiny hand-specified transcripts,
so there is no harness-vs-production gap. Reports per-region (oracle truth f_g, calibrate()-solved f_g).

API:  run(name, genes, kappa=, n_rna=, gdna_fraction=, capture=, capture_strength=, nascent=)
      genes = [("TA","+",[(1000,2000),(5000,10000)]), ...]
"""
from __future__ import annotations

import collections
from pathlib import Path

import numpy as np
import pysam

from rigel.sim import Scenario, ReadSimConfig, GDNAConfig, CaptureConfig
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.pipeline import scan_and_buffer
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.calibrate import calibrate
from rigel.calibration.signature import BIT_EXON_POS, BIT_EXON_NEG, BIT_INTRON_POS, BIT_INTRON_NEG
from rigel.splice import SpliceType

SCRATCH = Path("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/e1a42517-d9da-4b58-b12a-0b84f986ddef/scratchpad")


def _write_probes(path, ref, genes):
    with open(path, "w") as f:
        for gid, _strand, exons, *_rest in genes:
            for i, (s, e) in enumerate(exons):
                f.write(f"{ref}\t{s}\t{e}\t{gid}:p{i}\t0\t+\t{s}\t{e}\t0\t1\t{e - s}\t0\n")


def _truth_fg(bam_path, ra, ref_names):
    """Per-region oracle f_g = contained gDNA / (contained gDNA + contained RNA), via mate+TLEN span."""
    starts = ra.start; ends = ra.end; rid = ra.ref_id; R = starts.shape[0]
    by_ref = collections.defaultdict(list)
    for i in range(R):
        by_ref[int(rid[i])].append(i)
    for k in list(by_ref):
        arr = sorted(by_ref[k], key=lambda i: starts[i])
        by_ref[k] = (np.array([starts[i] for i in arr]), np.array([ends[i] for i in arr]), np.array(arr))
    name2ref = {n: i for i, n in enumerate(ref_names)}
    g = np.zeros(R); r = np.zeros(R)
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for rd in bam.fetch(until_eof=True):
        if rd.is_secondary or rd.is_supplementary or rd.is_unmapped or not rd.is_read1 or rd.mate_is_unmapped:
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
    return np.where((g + r) > 0, g / np.maximum(g + r, 1e-9), np.nan), g, r


def run(name, genes, *, kappa=1.0, n_rna=4000, gdna_fraction=0.5, capture=False,
        capture_strength=20.0, nascent=0.0, genome_length=12000, seed=7, instrument=False,
        config=None, _debug=None):
    wd = SCRATCH / f"prod_{name}"
    sc = Scenario(name, genome_length=genome_length, seed=seed, work_dir=wd, ref_name="chr1")
    for gid, strand, exons, *rest in genes:
        ab = float(rest[0]) if rest else 100.0  # optional per-gene abundance (0 = unexpressed, pure gDNA)
        sc.add_gene(gid, strand, [{"t_id": gid, "exons": exons, "abundance": ab}])
    cap_cfg = None
    if capture:
        probes = wd / "probes.bed"; wd.mkdir(parents=True, exist_ok=True)
        _write_probes(probes, "chr1", genes)
        cap_cfg = CaptureConfig(probes=str(probes), binding_per_base=capture_strength)
    result = sc.build_oracle(
        n_rna_fragments=n_rna, gdna_fraction=gdna_fraction,
        sim_config=ReadSimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                                 read_length=120, strand_specificity=kappa, seed=seed),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=250, frag_std=50),
        capture_config=cap_cfg, nrna_abundance=nascent,
    )
    idx = result.index
    if _debug is not None:
        _debug["bam_path"] = str(result.bam_path)  # oracle BAM, for downstream per-node oracle probes
    _st, sm, flm, _buf, pl = scan_and_buffer(str(result.bam_path), idx, BamScanConfig())
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    if _debug is not None:
        from rigel.calibration.effective_length import boundary_eff_length as _bel
        _debug["gdna_fl_mean"] = float(_bel(fl.gdna_pmf))  # gDNA mean fragment length (tiny-node threshold)
    ra = RegionArrays.from_index(idx)
    from rigel.calibration.strand_balance import fit_strand_balance
    kfit = float(fit_strand_balance(sm).rna_sense_frac)
    cap = {}
    if instrument:
        import sys as _sys
        calmod = _sys.modules["rigel.calibration.calibrate"]
        _orig = calmod.node_sweep
        calmod.node_sweep = lambda *a, **k: _orig(*a, _capture=cap, **k)
        try:
            cal = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, config or CalibrationConfig(), _debug=_debug)
        finally:
            calmod.node_sweep = _orig
    else:
        cal = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, config or CalibrationConfig(), _debug=_debug)
    if _debug is not None:  # expose the node chain to callers (calibrate already stashes _debug["capture"])
        from rigel.calibration.node_chain import build_node_chain
        _debug["chain"] = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    gc = np.asarray(cal.mass_gdna_contained); rc = np.asarray(cal.mass_rna_contained)
    solved = np.where((gc + rc) > 0, gc / np.maximum(gc + rc, 1e-9), np.nan)
    truth, tg, tr = _truth_fg(result.bam_path, ra, idx.ref_names)

    sig = np.asarray(ra.signature).astype(int)

    def rtype(s):
        if s & (BIT_EXON_POS | BIT_EXON_NEG):
            return "exon"
        if s & (BIT_INTRON_POS | BIT_INTRON_NEG):
            return "intron"
        return "intergenic"

    rdf = idx.nodes_df
    print(f"\n===== {name}: {kappa=} capture={'ON x'+str(capture_strength) if capture else 'OFF'} "
          f"n_rna={n_rna} gdna_frac={gdna_fraction} nascent={nascent}  (fit κ={kfit:.3f}) =====")
    print(f"  {'reg':>3} {'span':>13} {'type':>10} {'sig':>3} | {'truth':>6} {'solved':>6} {'err':>7}")
    for i in range(len(rdf)):
        st = rdf['start'].iloc[i]; en = rdf['end'].iloc[i]
        tt = float('nan') if np.isnan(truth[i]) else truth[i]
        ts = 'nan' if np.isnan(tt) else f"{tt:.3f}"
        err = '' if np.isnan(tt) else f"{solved[i] - tt:+.3f}"
        print(f"  {i:>3} {f'{st}-{en}':>13} {rtype(sig[i]):>10} {sig[i]:>3} | {ts:>6} {solved[i]:>6.3f} {err:>7}")
    if instrument and cap:
        from rigel.calibration.node_chain import build_node_chain, REGION
        ch = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
        kind = np.asarray(ch.kind); ridx = np.asarray(ch.ref_idx)
        fgl = np.asarray(cap["fg_loc"]); fgs = np.asarray(cap["fg_strand"])
        mg = np.asarray(cap["mode_g"]); pg = np.asarray(cap["prec_g"])
        mp = np.asarray(cap["mode_p"]); pp = np.asarray(cap["prec_p"])
        mn = np.asarray(cap["mode_n"]); pn = np.asarray(cap["prec_n"])
        eff_gl = np.asarray(cap["eff_gdna_l"]); eff_gr = np.asarray(cap["eff_gdna_r"])
        ml = np.asarray(cap["mass_l"]); mr = np.asarray(cap["mass_r"])
        fgn = np.asarray(cap["f_g"])
        print(f"  rho_global={cap['rho_global']:.4f}  (anchor=fg_loc, msg modes are log-fraction targets)")
        print(f"  {'node':>10} {'kind':>4} | {'fg_str':>7} {'fg_loc':>7} {'fg_fin':>7} | {'gMSGmo':>7} {'gMSGpr':>7} {'pMSGmo':>7} {'pMSGpr':>7}")
        for n in range(ch.n_nodes):
            isreg = kind[n] == REGION
            tag = (f"R{ridx[n]}" if isreg else f"B{ridx[n]}")
            # gDNA density this node presents (left face) for context
            print(f"  {tag:>10} {'reg' if isreg else 'bnd':>4} | {fgs[n]:7.3f} {fgl[n]:7.3f} {fgn[n]:7.3f}"
                  f" | {mg[n]:7.2f} {pg[n]:7.2f} {mp[n]:7.2f} {pp[n]:7.2f}")
    return rdf, solved, truth, tg, tr, sig


if __name__ == "__main__":
    TA = [("TA", "+", [(1000, 2000), (5000, 10000)])]
    for lab, kap, cap in [("S1", 1.0, False), ("S2", 0.5, False), ("S3", 1.0, True), ("S4", 0.5, True)]:
        run(f"TA_{lab}", TA, kappa=kap, n_rna=4000, gdna_fraction=0.5, capture=cap, capture_strength=20.0)
