"""Suite-level CALIBRATION-basis pool deconvolution benchmark (pre-EM).

Answers one question across every condition of a sim suite: **how well does production calibration
deconvolve each region's unspliced fragment mass into the 2-simplex ``(RNA₊, RNA₋, gDNA)``?** This is
the calibration layer ALONE — it does not run the EM, so it is not confounded by downstream
assignment / eff-length effects (use the ``net_flow`` 3-pool TSVs for the post-EM view).

It is built entirely on two VALIDATED primitives, so there is nothing new to get subtly wrong:

* truth  — :class:`oracle.OracleTruth` (the production accumulator, split by true fragment origin;
  ``region_pools()`` gives the per-region gDNA / mature / nascent × genome-strand contained mass).
* estimate — the production :func:`rigel.calibration.calibrate`, read through its ``_debug`` hook +
  :func:`chain_region_deconv` for the per-region ``(f_g, f₊, f₋)`` pie (the CalibrationResult itself
  does not carry the per-strand split).

We scan the full BAM ONCE and feed that payload to BOTH calibrate and the oracle (``full_payload=``),
so the estimate and the truth are measured on byte-identical input.

Metric (pool = summed over all region-contained nodes, the competition basis):
  * per-pool net surplus  ``Σ(cal − true)``  and dispersion ``Σ_r |cal_r − true_r|`` (net cancels
    per-region over/under-call; report both) — for gDNA, RNA₊, RNA₋.
  * mass-weighted ``|Δf_g|`` — the headline per-region gDNA-fraction error (0 = perfect).
  * mass-weighted simplex TV ``Σ_r M_r·½(|Δf_g|+|Δf₊|+|Δf₋|) / Σ_r M_r`` — the single 3-way
    deconvolution scalar (0 = perfect, 1 = worst), the right target for the AMBIG stress.
  * directional gDNA under-call (leaks to RNA) vs over-call.

Calibration models RNA-vs-gDNA only; mature vs nascent is the EM's job. ``mature`` / ``nascent`` are
therefore reported as the TRUE composition of the RNA the calibration lumps together (the RNA-pool
error is on the total). Spliced mass is guaranteed-RNA (no gDNA rival) and is excluded from the
unspliced competition basis; its true total is reported for context.

    OMP_NUM_THREADS=1 python scripts/debug/calib_pool_benchmark.py --suite DIR [--conditions a,b] \
        --out calib_pool_benchmark.tsv
"""
from __future__ import annotations

import argparse
import os
import pickle
from dataclasses import replace as dc
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd

from rigel.calibration import calibrate
from rigel.calibration.bp_solver import chain_region_deconv
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.splice import SpliceType

from oracle import OracleTruth  # noqa: E402  (same directory)


def _parse_condition(name: str) -> dict:
    """Decode the ``gdna_{gdna}_ss_{ss}_nrna_{nrna}_capture_{capture}`` condition name into its factor
    levels. Position-based on the fixed naming scheme (``ss``/``nrna``/``capture`` are found by their
    label; the gDNA level is the token after the leading ``gdna`` — NOT matched against a value list,
    since ``none`` also appears in ``nrna_none``)."""
    toks = name.split("_")
    get = lambda key, default="?": (  # noqa: E731
        toks[toks.index(key) + 1] if key in toks and toks.index(key) + 1 < len(toks) else default
    )
    gdna = toks[1] if len(toks) > 1 and toks[0] == "gdna" else "?"
    return dict(gdna=gdna, ss=get("ss"), nrna=get("nrna"), capture=get("capture"))


def _scan_and_truth(suite: Path, cond: str, index, cfg, work_dir: Path, cache_dir: Path | None) -> dict:
    """The EXPENSIVE, calibration-INDEPENDENT half: the production scan + the validated oracle truth.
    Cached to ``cache_dir`` (one pickle per condition) so iterating on the calibration code re-runs only
    ``calibrate`` + the metrics (seconds) instead of re-scanning + re-oracling every BAM (minutes).

    The cache holds exactly the calibrate inputs (payload / strand model / FL pmfs) and the oracle's
    per-region truth pools — none of which depend on the calibration being benchmarked. Delete the dir to
    invalidate (e.g. after a scanner/accumulator/oracle change)."""
    cache = (cache_dir / f"{cond}.pkl") if cache_dir is not None else None
    if cache is not None and cache.exists():
        with open(cache, "rb") as fh:
            return pickle.load(fh)

    bam = str(suite / cond / "sim_oracle.bam")
    # ---- ONE production scan; feed the same payload to calibrate AND the oracle ----
    sc = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _stats, sm, flm, _buf, payload = scan_and_buffer(bam, index, sc)
    sm.finalize()
    fl = build_fl_models(
        global_counts=flm.global_model.counts,
        rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(payload),
        max_size=flm.max_size,
    )
    # ---- validated per-origin truth on the SAME payload ----
    # boundary_mass is float32; over ~8M fragments its sum-to-full rounding (~1e-2) can exceed the
    # 1e-2 default (tuned for the 5Mb suite). The INTEGER region_contained check — the real proof the
    # partition is the production split, and the ONLY basis this metric uses — still holds exactly, so
    # a generous float tolerance (~1e-7 relative here) is safe.
    orc = OracleTruth.from_bam(bam, index, cfg, work_dir, cond, full_payload=payload,
                               boundary_mass_tol=0.5)
    out = dict(
        payload=payload, strand_model=sm,
        gdna_fl_pmf=np.asarray(fl.gdna_pmf), rna_fl_pmf=np.asarray(fl.rna_pmf),
        pools={k: np.asarray(v) for k, v in orc.region_pools().items()},
    )
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        with open(cache, "wb") as fh:
            pickle.dump(out, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return out


def evaluate_condition(suite: Path, cond: str, index, cfg, work_dir: Path,
                       cache_dir: Path | None = None) -> dict:
    """Run production calibration + the validated oracle on one condition → a pool-error row."""
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    inp = _scan_and_truth(suite, cond, index, cfg, work_dir, cache_dir)
    payload = inp["payload"]

    dbg: dict = {}
    calibrate(payload=payload, region_arrays=ra, strand_model=inp["strand_model"],
              gdna_fl_pmf=inp["gdna_fl_pmf"], rna_fl_pmf=inp["rna_fl_pmf"],
              config=cfg.calibration, _debug=dbg)
    # per-region solved pie (f_g, f₊, f₋) over the contained UNSPLICED mass
    reg = chain_region_deconv(dbg["chain"], dbg["belief"], dbg["substrate"])
    M = np.asarray(dbg["substrate"].contained.mass_unspliced, np.float64)  # contained unspliced total
    cal_g = reg.gdna_frac * M
    cal_p = reg.rna_pos_frac * M
    cal_n = reg.rna_neg_frac * M

    p = inp["pools"]
    true_g = p["gdna_pos"] + p["gdna_neg"]
    true_p = p["mat_uns_pos"] + p["nas_uns_pos"]
    true_n = p["mat_uns_neg"] + p["nas_uns_neg"]
    true_mat_uns = p["mat_uns_pos"] + p["mat_uns_neg"]
    true_nas_uns = p["nas_uns_pos"] + p["nas_uns_neg"]
    true_spl = p["mat_spl"] + p["nas_spl"]

    # conservation: contained unspliced total == gDNA + RNA(±) truth (validated sum-to-full)
    resid = float(np.abs(M - (true_g + true_p + true_n)).max())
    if resid > 1e-6:
        raise AssertionError(f"{cond}: unspliced conservation off by {resid:.3e}")

    # ---- strand-orientation alignment ----
    # calibration's f_pos/f_neg are TRANSCRIPT-strand (sense/antisense; free_pos gates on the region's
    # ± transcript bits), while the accumulator channels ch0/ch1 (and thus true_p/true_n) are GENOME
    # strand. For a stranded protocol transcript-'+' RNA maps to a fixed genome strand, so the two
    # frames differ by ONE global swap (verified: a clean pool-level swap, not per-region noise). We
    # detect that single library-level bit by whichever labelling fits the RNA pool better and align
    # the truth to the calibration frame, so the per-strand error is the REAL deconvolution error
    # (not the convention offset). f_g and the RNA total are orientation-invariant.
    ident = abs(cal_p.sum() - true_p.sum()) + abs(cal_n.sum() - true_n.sum())
    swap = abs(cal_p.sum() - true_n.sum()) + abs(cal_n.sum() - true_p.sum())
    swapped = swap < ident
    if swapped:
        true_p, true_n = true_n, true_p

    # ---- per-region fractions on the unspliced simplex (only where there is mass) ----
    m = M > 0
    w = M[m]
    tfg = true_g[m] / M[m]
    tfp = true_p[m] / M[m]
    tfn = true_n[m] / M[m]
    dfg = np.abs(reg.gdna_frac[m] - tfg)
    dfp = np.abs(reg.rna_pos_frac[m] - tfp)
    dfn = np.abs(reg.rna_neg_frac[m] - tfn)
    wsum = max(float(w.sum()), 1.0)
    mwae_fg = float(np.sum(w * dfg) / wsum)
    tv = float(np.sum(w * 0.5 * (dfg + dfp + dfn)) / wsum)  # aligned mass-wt simplex TV

    row = dict(condition=cond, **_parse_condition(cond))
    row.update(
        # pool totals (contained unspliced competition basis)
        true_gdna=true_g.sum(), true_rna=true_p.sum() + true_n.sum(),
        true_mature=true_mat_uns.sum(), true_nascent=true_nas_uns.sum(),
        true_rna_pos=true_p.sum(), true_rna_neg=true_n.sum(),
        true_spliced=true_spl.sum(),
        cal_gdna=cal_g.sum(), cal_rna=cal_p.sum() + cal_n.sum(),
        cal_rna_pos=cal_p.sum(), cal_rna_neg=cal_n.sum(),
        # pool-level calibration error
        gdna_net=cal_g.sum() - true_g.sum(),
        gdna_abs=float(np.abs(cal_g - true_g).sum()),
        rna_pos_net=cal_p.sum() - true_p.sum(),
        rna_neg_net=cal_n.sum() - true_n.sum(),
        rna_pos_abs=float(np.abs(cal_p - true_p).sum()),
        rna_neg_abs=float(np.abs(cal_n - true_n).sum()),
        gdna_undercall=float(np.maximum(true_g - cal_g, 0.0).sum()),  # gDNA→RNA leak
        gdna_overcall=float(np.maximum(cal_g - true_g, 0.0).sum()),   # RNA→gDNA siphon
        mwae_fg=mwae_fg, simplex_tv=tv,
        unspliced_total=float(M.sum()),
        strand_orient=("swap" if swapped else "identity"),
    )
    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--conditions", default=None,
                    help="comma-separated subset (default: all sim_oracle.bam dirs)")
    ap.add_argument("--out", default=None, help="output TSV (default: <suite>/calib_pool_benchmark.tsv)")
    ap.add_argument(
        "--cache-dir", default=None,
        help="cache the per-condition scan + oracle truth here (calibration-independent), so re-running "
        "after a calibration change costs only calibrate + metrics (seconds, not minutes). Delete the dir "
        "to invalidate after a scanner/accumulator/oracle change.",
    )
    args = ap.parse_args()

    suite = Path(args.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work_dir = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_calib_pool_split"
    cache_dir = Path(args.cache_dir) if args.cache_dir else None

    if args.conditions:
        conds = args.conditions.split(",")
    else:
        conds = sorted(d.name for d in suite.iterdir()
                       if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_"))

    rows = []
    for i, cond in enumerate(conds, 1):
        print(f"[{i}/{len(conds)}] {cond} ...", flush=True)
        rows.append(evaluate_condition(suite, cond, index, cfg, work_dir, cache_dir))

    df = pd.DataFrame(rows)
    out = Path(args.out) if args.out else suite / "calib_pool_benchmark.tsv"
    df.to_csv(out, sep="\t", index=False)

    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 40)
    show = ["condition", "true_gdna", "cal_gdna", "gdna_net", "gdna_abs",
            "true_rna_pos", "cal_rna_pos", "rna_pos_net",
            "true_rna_neg", "cal_rna_neg", "rna_neg_net", "mwae_fg", "simplex_tv"]
    print("\n=== CALIBRATION-basis pool error (contained unspliced) ===")
    print(df[show].to_string(index=False, float_format=lambda x: f"{x:,.4g}"))
    print(f"\nfull table -> {out}")


if __name__ == "__main__":
    main()
