"""Deterministic per-region / per-boundary calibration dissection vs an oracle.

Splits a benchmark condition's ``sim_oracle.bam`` by TRUE origin (gdna-only / rna-only, from the
simulator read name), scans each separately to get the per-region oracle composition (using the
production accumulator deposit), runs the production ``calibrate()`` on the full scan, and prints
production-vs-oracle per region (and per boundary side), ranked by gDNA-mass error.

Determinism: invoke with ``OMP_NUM_THREADS=1`` so the scan accumulator (and thus the whole pipeline)
is bit-reproducible across reruns on the same BAM. The three scans are cached (pickle) so fix-testing
re-runs only the fast ``calibrate()`` on the cached payloads.

    OMP_NUM_THREADS=1 python scripts/debug/dissect_regions.py [condition] [--top-n N] [--rebuild] [--locus L]

The default condition is the gDNA-3x / capture-on / ss-0.99 / no-nascent case.
"""
from __future__ import annotations

import argparse
import os
import pickle
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SUITE = Path.home() / "Downloads/rigel_runs/quick_3to1_5mb"
CACHE = Path("/tmp/dissect_cache")
_EPS = 1e-9


def split_bam_by_origin(src: Path, gdna_out: Path, rna_out: Path) -> None:
    """Split a name-sorted BAM into gDNA-only and RNA-only by the simulator read name.

    gDNA reads are named ``gdna:...``; everything else (``t_id:...`` mRNA, ``nrna_...`` nascent) is RNA.
    awk keeps the @ header and both mates (mates share the qname), preserving name-sort + NH tags.
    """
    if gdna_out.exists() and rna_out.exists():
        return
    subprocess.run(
        f"samtools view -h '{src}' | awk '/^@/||$1~/^gdna:/' | samtools view -b -o '{gdna_out}' -",
        shell=True, check=True,
    )
    subprocess.run(
        f"samtools view -h '{src}' | awk '/^@/||$1!~/^gdna:/' | samtools view -b -o '{rna_out}' -",
        shell=True, check=True,
    )


def scan_payload(bam: Path, index, scan_cfg):
    """Scan one BAM → (AccumulatorPayload, StrandModels, FragmentLengthModels)."""
    from rigel.pipeline import scan_and_buffer

    _stats, strand, fl, buffer, payload = scan_and_buffer(str(bam), index, scan_cfg)
    del buffer  # not needed for calibration
    return payload, strand, fl


def build_or_load_cache(cond: str, rebuild: bool):
    from rigel.index import TranscriptIndex
    from rigel.config import BamScanConfig
    from rigel.calibration.fl import build_fl_models, gdna_fl_mass
    from rigel.splice import SpliceType

    index = TranscriptIndex.load(SUITE / "rigel_index")
    scan_cfg = BamScanConfig()
    ck = CACHE / f"{cond}.pkl"
    if ck.exists() and not rebuild:
        with open(ck, "rb") as fh:
            blob = pickle.load(fh)
        return index, blob

    CACHE.mkdir(parents=True, exist_ok=True)
    cdir = SUITE / cond
    gdna_bam, rna_bam = CACHE / f"{cond}.gdna.bam", CACHE / f"{cond}.rna.bam"
    print(f"[scan] splitting {cdir/'sim_oracle.bam'} by origin ...", file=sys.stderr)
    split_bam_by_origin(cdir / "sim_oracle.bam", gdna_bam, rna_bam)

    print("[scan] full ...", file=sys.stderr)
    payload_full, strand_full, fl_full = scan_payload(cdir / "sim_oracle.bam", index, scan_cfg)
    print("[scan] gdna-only ...", file=sys.stderr)
    payload_gdna, _, _ = scan_payload(gdna_bam, index, scan_cfg)
    print("[scan] rna-only ...", file=sys.stderr)
    payload_rna, _, _ = scan_payload(rna_bam, index, scan_cfg)

    fl_models = build_fl_models(
        global_counts=fl_full.global_model.counts,
        rna_counts=fl_full.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(payload_full),
        max_size=fl_full.max_size,
    )
    blob = dict(
        payload_full=payload_full, payload_gdna=payload_gdna, payload_rna=payload_rna,
        strand_full=strand_full, gdna_pmf=fl_models.gdna_pmf, rna_pmf=fl_models.rna_pmf,
    )
    try:
        with open(ck, "wb") as fh:
            pickle.dump(blob, fh)
        print(f"[scan] cached → {ck}", file=sys.stderr)
    except Exception as e:  # noqa: BLE001
        print(f"[scan] cache pickle failed ({e}); will re-scan next run", file=sys.stderr)
    return index, blob


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("condition", nargs="?", default="gdna_gdna300_ss_0.99_nrna_none_capture_on")
    ap.add_argument("--top-n", type=int, default=25)
    ap.add_argument("--rebuild", action="store_true")
    ap.add_argument("--kill-rna-prec", action="store_true",
                    help="EXPERIMENT: monkeypatch fit_rna_varmean to emit ~0 precision (honest-σ²_bio "
                         "direction) — tests whether over-confident RNA messages cause the AMBIG f_g→0.")
    ap.add_argument("--soften-global", type=float, default=0.0,
                    help="EXPERIMENT: scale σ_global UP by this factor (↓ global precision by factor²) — "
                         "tests whether the over-precise global prior pins AMBIG f_g→0.")
    ap.add_argument("--global-all-regions", action="store_true",
                    help="EXPERIMENT: compute rho_global over ALL regions with data (not count-observable "
                         "only) — the agreed capture-aware baseline.")
    args = ap.parse_args()

    from rigel.calibration import calibrate
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.substrate import CalibrationSubstrate
    from rigel.config import CalibrationConfig

    index, blob = build_or_load_cache(args.condition, args.rebuild)
    region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

    if args.kill_rna_prec:
        import rigel.calibration.bp_solver as bp

        class _LowPrec:  # σ²_bio→∞ ⇒ message τ→0 (RNA imputation effectively ignored)
            def predict(self, mu):
                return np.full(np.asarray(mu, float).shape, 1e6)

        bp.fit_rna_varmean = lambda *a, **k: _LowPrec()
        print("  [EXPERIMENT] RNA-message precision killed (fit_rna_varmean σ²_bio→∞)")

    if args.soften_global > 0:
        import rigel.calibration.bp_solver as bp
        _orig_g = bp.global_gdna_prior

        def _soft(*a, **k):
            rho, sig = _orig_g(*a, **k)
            return rho, sig * args.soften_global  # ↓ precision by factor²

        bp.global_gdna_prior = _soft
        print(f"  [EXPERIMENT] global σ scaled ×{args.soften_global} (precision ÷{args.soften_global**2:.0f})")

    if args.global_all_regions:
        import rigel.calibration.bp_solver as bp
        _orig_ga = bp.global_gdna_prior

        def _allreg(region_f_g, mass, eff, obs):  # rho over ALL regions with data, not count-observable only
            return _orig_ga(region_f_g, mass, eff, np.asarray(mass, float) > 0.0)

        bp.global_gdna_prior = _allreg
        print("  [EXPERIMENT] rho_global computed over ALL regions with data (capture-aware baseline)")

    cal = calibrate(
        payload=blob["payload_full"], region_arrays=region_arrays,
        strand_model=blob["strand_full"], gdna_fl_pmf=blob["gdna_pmf"],
        rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig(),
    )

    sub_g = CalibrationSubstrate.from_payload(blob["payload_gdna"], region_arrays)
    sub_r = CalibrationSubstrate.from_payload(blob["payload_rna"], region_arrays)
    sub_f = CalibrationSubstrate.from_payload(blob["payload_full"], region_arrays)

    rdf = index.region_df.reset_index(drop=True)

    # contained-node oracle vs production
    g_or = np.asarray(sub_g.contained.mass_unspliced, float)
    r_or = np.asarray(sub_r.contained.mass_unspliced, float)
    g_pr = np.asarray(cal.mass_gdna_contained, float)
    r_pr = np.asarray(cal.mass_rna_contained, float)
    fg_or = g_or / (g_or + r_or + _EPS)
    fg_pr = g_pr / (g_pr + r_pr + _EPS)
    gerr = g_pr - g_or  # + = production OVER-calls gDNA; − = UNDER-calls

    # raw strand counts (full scan) — the strand balance the deconv sees
    upos = np.asarray(sub_f.contained.n_unspliced_pos, float)
    uneg = np.asarray(sub_f.contained.n_unspliced_neg, float)
    spl = np.asarray(sub_f.contained.n_spliced, float)

    from rigel.calibration.signature import coarse_strand_from_signature, coarse_type_from_signature
    sig = rdf["signature"].to_numpy()
    strand_cls = np.array([coarse_strand_from_signature(int(s)) for s in sig])
    type_cls = np.array([coarse_type_from_signature(int(s)) for s in sig])
    SC = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}
    TC = {0: "intergenic", 1: "intron", 2: "exon"}

    print(f"=== {args.condition} : contained-node gDNA deconv (production − oracle) ===")
    print(f"  totals: oracle gDNA(contained)={g_or.sum():,.0f}  prod gDNA(contained)={g_pr.sum():,.0f}  "
          f"net err={gerr.sum():,.0f}")
    print(f"  gdna_density_global={cal.gdna_density_global:.4g}  rna_sense_frac={cal.rna_sense_frac:.3f}")

    # ── error aggregated by strand class × region type ──
    print(f"\n  gDNA-mass error by (strand_class × region_type):")
    print(f"  {'class':>6} {'type':>11} {'nreg':>5} {'g_oracle':>11} {'g_prod':>11} {'net_err':>11} {'|err|':>11}")
    for sc in (3, 1, 2, 0):
        for tc in (2, 1, 0):
            m = (strand_cls == sc) & (type_cls == tc)
            if not m.any():
                continue
            print(f"  {SC[sc]:>6} {TC[tc]:>11} {int(m.sum()):>5} {g_or[m].sum():>11,.0f} "
                  f"{g_pr[m].sum():>11,.0f} {gerr[m].sum():>+11,.0f} {np.abs(gerr[m]).sum():>11,.0f}")

    order = np.argsort(-np.abs(gerr))[: args.top_n]
    print(f"\n  top-{args.top_n} regions by |gDNA err|  (upos/uneg = raw unspliced strand counts):")
    print(f"  {'reg':>5} {'sig':>3} {'class':>5} {'start':>9} {'mass_u':>9} {'g_orac':>9} "
          f"{'g_prod':>9} {'gERR':>9} {'fg_or':>6} {'fg_pr':>6} {'upos':>7} {'uneg':>7} {'spl':>6}")
    for i in order:
        mu = g_or[i] + r_or[i]
        print(f"  {i:>5} {int(sig[i]):>3} {SC[int(strand_cls[i])]:>5} {int(rdf.start.iloc[i]):>9} "
              f"{mu:>9,.0f} {g_or[i]:>9,.0f} {g_pr[i]:>9,.0f} {gerr[i]:>+9,.0f} "
              f"{fg_or[i]:>6.3f} {fg_pr[i]:>6.3f} {upos[i]:>7,.0f} {uneg[i]:>7,.0f} {spl[i]:>6,.0f}")

    over = gerr[gerr > 0].sum()
    under = gerr[gerr < 0].sum()
    print(f"\n  TOTAL contained gDNA over-call (+)={over:,.0f}  under-call (−)={under:,.0f}")


if __name__ == "__main__":
    main()
