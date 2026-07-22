"""Toolkit: seed POPULATION calibration priors from a cached genome-scale scenario, then inject them into a
small FULL-TRANSCRIPT toy for controlled local testing of the pass-0 message propagation.

WHY: a tiny toy cannot fit rna_sense_frac (needs genome-scale spliced), the enrichment-density NPMLE (needs a
population to see enriched-vs-depleted modes), or the intergenic ρ_bg — all directly-observable population
priors that GROUND the solver (owner, 2026-07-22). So we (1) extract them from an ambig_dense_10mb condition's
cache via calibrate(_debug=...)["calibration_priors"], and (2) inject them into calibrate() on the toy.

THE TOY is a FULL 3-exon transcript with intergenic ENDS (owner's design) — a complete, self-grounded problem:
  intergenic — B — exon1(TSS) — B — intron — B — exon2 — B — intron — B — exon3(TES) — B — intergenic
The intergenic ends seed the baseline gDNA level (ρ_bg) that propagates through the messages. The toy is
generated at ambig-matched FRAGMENT LENGTHS (RNA=200, gDNA=100) + capture (binding_per_base=10) so its densities
land on the injected absolute-log-density NPMLE landscape; `--n-rna` tunes DEPTH to match ρ_bg.

Usage:
  python toy_inject.py                       # inject capture-on priors into a matched toy; compare vs toy-fit
  python toy_inject.py --cond <condition>    # choose the reference condition
"""
from __future__ import annotations
import argparse
import collections
import contextlib
import io
import pickle
from pathlib import Path

import numpy as np
import pysam

from rigel.sim import Scenario, ReadSimConfig, GDNAConfig, CaptureConfig
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.pipeline import scan_and_buffer
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.calibrate import calibrate
from rigel.calibration.node_chain import build_node_chain, REGION
from rigel.calibration.signature import coarse_type_array
from rigel.index import TranscriptIndex
from rigel.splice import SpliceType

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
SCRATCH = Path("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/4f7a248b-0c78-4b40-9030-462373aefb19/scratchpad")
_EPS = 1e-9
# ambig_dense_10mb generative params (scripts/sim/configs/ambig_dense_10mb.yaml) — MATCH these in the toy:
RNA_FL = dict(frag_mean=200, frag_std=0, frag_min=100, frag_max=300, read_length=100)
GDNA_FL = dict(frag_mean=100, frag_std=0)
CAP = dict(binding_per_base=10.0, gdna_split_penalty=0.2, off_target_weight=1.0)


def extract_priors(cond: str):
    """Extract the population priors fitted by calibrate on a cached ambig_dense_10mb condition."""
    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    cfg = CalibrationConfig()
    with open(SUITE / "_selfsolve_cache" / f"{cond}.pkl", "rb") as fh:
        c = pickle.load(fh)
    dbg = {}
    calibrate(c["payload"], ra, c["strand_model"], np.asarray(c["gdna_fl_pmf"]),
              np.asarray(c["rna_fl_pmf"]), cfg, _debug=dbg)
    return dbg["calibration_priors"]


def _write_probes(path, ref, exons, fraction=1.0):
    """Probe BED: cover the CENTRAL `fraction` of each exon.

    NOTE (measured 2026-07-22): `fraction` is a LOAD-BEARING sweep dimension, not cosmetic. Splice-junction
    reads live at exon EDGES, so a centred partial probe (fraction=0.5) DEPLETES them ~300× (331 → 1 spliced
    alignments in the same library) and with them the boundary's spliced measurement ``S_B``. Default to
    full-exon probes so the junction channel survives; sweep `fraction` to test graceful degradation."""
    with open(path, "w") as f:
        for i, (s, e) in enumerate(exons):
            w = int((e - s) * fraction)
            ps = s + ((e - s) - w) // 2
            pe = ps + w
            f.write(f"{ref}\t{ps}\t{pe}\tp{i}\t0\t+\t{ps}\t{pe}\t0\t1\t{pe - ps}\t0\n")


def build_toy(name, *, exons, strand="+", gdna_fraction, nascent, mature, capture_on,
              n_rna, genome_length, kappa=0.5, seed=7):
    """Generate a full-transcript toy at ambig-matched FL + capture; scan → calibrate inputs + oracle + chain."""
    wd = SCRATCH / f"inj_{name}"
    sc = Scenario(name, genome_length=genome_length, seed=seed, work_dir=wd, ref_name="chr1")
    sc.add_gene("G", strand, [{"t_id": "G", "exons": exons, "abundance": float(mature)}])
    cap_cfg = None
    if capture_on:
        wd.mkdir(parents=True, exist_ok=True)
        probes = wd / "probes.bed"
        _write_probes(probes, "chr1", exons)
        cap_cfg = CaptureConfig(probes=str(probes), **CAP)
    result = sc.build_oracle(
        n_rna_fragments=n_rna, gdna_fraction=gdna_fraction,
        sim_config=ReadSimConfig(strand_specificity=kappa, seed=seed, **RNA_FL),
        gdna_config=GDNAConfig(abundance=0.0, **GDNA_FL),
        capture_config=cap_cfg, nrna_abundance=float(nascent),
    )
    idx = result.index
    _st, sm, flm, _buf, pl = scan_and_buffer(str(result.bam_path), idx, BamScanConfig())
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    truth = _truth_fg(result.bam_path, ra, idx.ref_names)
    return dict(payload=pl, ra=ra, strand_model=sm, gdna_fl_pmf=np.asarray(fl.gdna_pmf),
                rna_fl_pmf=np.asarray(fl.rna_pmf), chain=chain, truth=truth, rdf=idx.region_df)


def _truth_fg(bam_path, ra, ref_names):
    starts = ra.start; ends = ra.end; rid = ra.ref_id; R = starts.shape[0]
    by_ref = collections.defaultdict(list)
    for i in range(R):
        by_ref[int(rid[i])].append(i)
    for k in list(by_ref):
        arr = sorted(by_ref[k], key=lambda i: starts[i])
        by_ref[k] = (np.array([starts[j] for j in arr]), np.array([ends[j] for j in arr]), np.array(arr))
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
        s, e, ii = by_ref[rf]; j = int(np.searchsorted(s, lo, side="right") - 1)
        if 0 <= j < R and s[j] <= lo and hi <= e[j]:
            (g if is_g else r)[ii[j]] += 1.0
    bam.close()
    return np.where((g + r) > 0, g / np.maximum(g + r, _EPS), np.nan)


def solve(toy, priors=None, refit=0):
    cfg = CalibrationConfig() if refit else __import__("dataclasses").replace(CalibrationConfig(), calib_refit_iters=0)
    dbg = {}
    calibrate(toy["payload"], toy["ra"], toy["strand_model"], toy["gdna_fl_pmf"], toy["rna_fl_pmf"],
              cfg, _debug=dbg, injected_priors=priors)
    return dbg


def _report(toy, dbg, label):
    ra = toy["ra"]; chain = toy["chain"]; truth = toy["truth"]; rdf = toy["rdf"]
    cap = dbg["capture"]
    ctype = coarse_type_array(np.asarray(ra.signature)).astype(int)
    kind = np.asarray(chain.kind); ref_idx = np.asarray(chain.ref_idx)
    fg = np.asarray(cap["f_g"])
    reg2node = {int(ref_idx[n]): n for n in range(chain.n_nodes) if kind[n] == REGION}
    tn = {0: "intergenic", 1: "intron", 2: "exon"}
    print(f"\n  [{label}]  per-region solve vs oracle:")
    print(f"    {'reg':>3} {'type':>10} {'span':>15} | {'oracle':>6} {'solved':>6} {'err':>6}")
    for i in range(len(rdf)):
        if i not in reg2node:
            continue
        o = truth[i]
        os_ = "  nan" if np.isnan(o) else f"{o:.3f}"
        err = "" if np.isnan(o) else f"{fg[reg2node[i]] - o:+.2f}"
        span = f"{int(rdf['start'].iloc[i])}-{int(rdf['end'].iloc[i])}"
        print(f"    {i:>3} {tn.get(int(ctype[i]), '?'):>10} {span:>15} | {os_:>6} {fg[reg2node[i]]:>6.3f} {err:>6}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cond", default="gdna_gdna300_ss_0.50_nrna_present_capture_on")
    ap.add_argument("--n-rna", type=int, default=20000)
    ap.add_argument("--gdna", type=float, default=3.0)
    ap.add_argument("--nascent", type=float, default=200.0)
    ap.add_argument("--mature", type=float, default=100.0)
    args = ap.parse_args()

    capture_on = "capture_on" in args.cond
    print(f"reference condition (priors source): {args.cond}   capture_on={capture_on}")
    priors = extract_priors(args.cond)
    print(f"  injected: κ={priors.rna_sense_frac:.4f}  od_g={priors.gdna_strand_overdispersion:.4f}  "
          f"od_r={priors.rna_strand_overdispersion:.4f}  NPMLE_cells={priors.enrichment_prior.n_cells}  "
          f"n_gdna_obs={priors.n_gdna_obs:.0f}  bg={'yes' if priors.background else 'no'}")

    # full 3-exon transcript with intergenic ends (genome 0-30000; gene exons in the middle)
    exons = [(6000, 7000), (11000, 12000), (16000, 17000)]
    with contextlib.redirect_stdout(io.StringIO()):
        toy = build_toy("full_tx", exons=exons, gdna_fraction=args.gdna, nascent=args.nascent,
                        mature=args.mature, capture_on=capture_on, n_rna=args.n_rna, genome_length=30000)

    # scale diagnostic: the toy's intergenic gDNA density should ≈ the reference ρ_bg
    from rigel.calibration.bp_solver import build_node_geometry, node_global_geometry
    from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
    sub = CalibrationSubstrate.from_payload(toy["payload"], toy["ra"])
    bsub = BoundarySubstrate.from_payload(toy["payload"])
    geo = build_node_geometry(toy["chain"], sub, bsub, toy["ra"], toy["gdna_fl_pmf"], toy["rna_fl_pmf"])
    mass, eff = node_global_geometry(toy["chain"], geo)
    ctype = coarse_type_array(np.asarray(toy["ra"].signature)).astype(int)
    kind = np.asarray(toy["chain"].kind); ref_idx = np.asarray(toy["chain"].ref_idx)
    dens = np.where(eff > _EPS, mass / np.maximum(eff, _EPS), 0.0)
    ig = dens[(kind == REGION) & (ctype[np.clip(ref_idx, 0, len(ctype) - 1)] == 0)]
    ex = dens[(kind == REGION) & (ctype[np.clip(ref_idx, 0, len(ctype) - 1)] == 2)]
    print(f"  toy scale: intergenic density med={np.median(ig[ig > 0]) if (ig > 0).any() else 0:.4f}  "
          f"exon density med={np.median(ex[ex > 0]) if (ex > 0).any() else 0:.4f}  "
          f"(match ref intergenic ρ_bg — tune --n-rna)")

    _report(toy, solve(toy, priors=None), "TOY-FIT priors, Phase-1")
    _report(toy, solve(toy, priors=priors), "INJECTED priors, Phase-1 (pass-0)")
    _report(toy, solve(toy, priors=priors, refit=1), "INJECTED priors, +Phase-2 hyperprior (ρ_bg peel)")


if __name__ == "__main__":
    main()
