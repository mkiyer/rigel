"""Full end-to-end A/B benchmark: baseline vs a calibration change, across all 16 conditions.

For each arm (baseline / candidate) and condition, runs `rigel quant` with a chosen
`--gdna-prior-mixture-bridge` epsilon, then collects the FULL set of tool-output metrics the
benchmark should report before shipping a change:

  * 3-POOL NET SURPLUS (assigned − true), fragments, for gDNA / nascent / mature — the pool-level
    deconvolution result (soft EM counts; sensitive to a calibration-prior change, unlike the
    hard-label net-flow). Assigned from summary.json `quantification`; true from truth_summary.json
    `origin_counts`.
  * ABSOLUTE mature error Σ_tx |measured − true| (positive/negative per-transcript flow cancels in
    the net, so the absolute sum is reported alongside).
  * TRANSCRIPT-level Spearman (truth vs measured), MARD, n_FP, n_FN — the vetted abundance-accuracy
    metrics from `rigel.sim.analysis` (abundance_per_condition.tsv).

Writes a structured JSON (both arms, all conditions) consumed by the HTML report renderer.

    python scripts/debug/benchmark_ab_report.py <suite> --out ab_report.json \
        --arms baseline:0 fix1:0.01 [--threads 4]
"""
import argparse, json, subprocess, sys, time
from pathlib import Path
import numpy as np, pandas as pd


def run_quant(cond_dir: Path, index: Path, eps: float, threads: int, n_grid_ss: int | None = None) -> None:
    out = cond_dir / "rigel_out"
    cmd = ["rigel", "quant", "--bam", str(cond_dir / "sim_oracle.bam"), "--index", str(index),
           "-o", str(out), "--annotated-bam", str(cond_dir / "annotated.bam"),
           "--sj-strand-tag", "auto", "--emit-locus-stats", "--tsv",
           "--gdna-prior-mixture-bridge", str(eps), "--threads", str(threads)]
    if n_grid_ss is not None:
        cmd += ["--sweep-n-grid-single-strand", str(n_grid_ss)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"quant failed for {cond_dir.name} (eps={eps}):\n{r.stderr[-1500:]}")


def _observed_truth(cond_dir: Path) -> dict:
    """OBSERVED fragment counts per origin from manifest.json (`n_*_observed`) — the correct pool truth.
    origin_counts in truth_summary.json is the MOLECULAR count (same for capture on/off), which over-states
    the observed gDNA on capture-off (most gDNA is never sequenced), so it must NOT be used as the pool truth."""
    manifest = json.load(open(cond_dir.parent / "manifest.json"))
    for cd in manifest.get("conditions", []):
        if cd.get("name") == cond_dir.name:
            return {"gdna": float(cd.get("n_gdna_observed", 0)), "nrna": float(cd.get("n_nrna_observed", 0)),
                    "mrna": float(cd.get("n_mrna_observed", 0))}
    return {"gdna": 0.0, "nrna": 0.0, "mrna": 0.0}


def pool_metrics(cond_dir: Path) -> dict:
    """3-pool assigned (summary EM counts) vs OBSERVED true (manifest n_*_observed) + surplus."""
    q = json.load(open(cond_dir / "rigel_out" / "summary.json")).get("quantification", {})
    # gDNA pool = in-locus gDNA + intergenic gDNA (the summary tracks them separately; on capture-off most
    # gDNA is intergenic, so it MUST be included or the gDNA pool looks 60% under-called).
    assigned = {"gdna": float(q.get("gdna_total", 0)) + float(q.get("intergenic_total", 0)),
                "nrna": float(q.get("nrna_total", 0)), "mrna": float(q.get("mrna_total", 0))}
    true = _observed_truth(cond_dir)
    return {p: {"assigned": assigned[p], "true": true[p], "surplus": assigned[p] - true[p]}
            for p in ("gdna", "nrna", "mrna")}


def transcript_metrics(cond_dir: Path) -> dict:
    """MATURE-transcript accuracy vs truth, abundance-aware. Mature RNA is the tool's priority, so we
    isolate it (drop the `is_nrna` rows quant.feather carries) and report accuracy on the transcripts that
    matter — the abundant ones — separately from the low-abundance detection-limit noise.

    Measured = quant.feather `count_em` on MATURE rows (`is_nrna == False`); truth =
    `observed_mrna_fragments` (post-capture observed mature per transcript). Metrics:
      * spearman_expr / mard_expr — over EXPRESSED transcripts (true > 0).
      * spearman_abund / mard_abund — over the top-abundance quartile (the transcripts that matter).
      * n_fp / abs_mature_err — false positives (measured but true=0) and Σ|measured − true|.
    """
    from scipy.stats import spearmanr
    # quant.feather IS the mature-transcript table (synthetic nascent shadows live separately in
    # nrna_quant.feather). It includes annotated single-exon transcripts flagged is_nrna — those are mature
    # AND nascent at once (unspliced, indistinguishable), so they are KEPT as mature transcripts.
    q = pd.read_feather(cond_dir / "rigel_out" / "quant.feather")
    mcol = "count_em" if "count_em" in q.columns else "count"
    tr = pd.read_csv(cond_dir / "truth_abundances.tsv", sep="\t")
    tcol = "observed_mrna_fragments" if "observed_mrna_fragments" in tr.columns else "post_capture_mrna_fragments"
    j = q[["transcript_id", mcol]].merge(tr[["transcript_id", tcol]], on="transcript_id", how="left").fillna(0.0)
    m = j[mcol].to_numpy(float); t = j[tcol].to_numpy(float)

    def _sp(mm, tt):
        return float(spearmanr(mm, tt).statistic) if len(tt) > 2 and np.ptp(tt) > 0 else float("nan")
    expr = t > 0
    me, te = m[expr], t[expr]
    sp_expr = _sp(me, te)
    mard_expr = float(np.mean(np.abs(me - te) / te)) if expr.any() else float("nan")
    # top-abundance quartile (the transcripts that matter)
    abund = te >= np.quantile(te, 0.75) if expr.sum() >= 4 else np.ones(expr.sum(), bool)
    sp_abund = _sp(me[abund], te[abund])
    mard_abund = float(np.mean(np.abs(me[abund] - te[abund]) / te[abund])) if abund.any() else float("nan")
    n_fp = int(np.sum((m > 0.5) & (t == 0)))
    abs_err = float(np.sum(np.abs(m - t)))
    net_err = float(np.sum(m - t))
    return {"spearman_expr": sp_expr, "mard_expr": mard_expr,
            "spearman_abund": sp_abund, "mard_abund": mard_abund,
            "n_expr": int(expr.sum()), "n_fp": n_fp,
            "abs_mature_err": abs_err, "net_mature_err": net_err}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("suite", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--arms", nargs="+", default=["baseline:0", "fix1:0.01"],
                    help="arm spec name:eps[:n_grid_ss] (n_grid_ss optional → config default)")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    index = args.suite / "rigel_index"
    conds = sorted(p for p in args.suite.glob("gdna_*") if (p / "sim_oracle.bam").exists())

    def _parse(a):
        parts = a.split(":")
        return (parts[0], float(parts[1]), int(parts[2]) if len(parts) > 2 else None)
    arms = [_parse(a) for a in args.arms]
    report: dict = {"conditions": [c.name for c in conds], "arms": [a for a, _, _ in arms], "data": {}}
    for arm, eps, ngss in arms:
        report["data"][arm] = {}
        for c in conds:
            t0 = time.time()
            run_quant(c, index, eps, args.threads, n_grid_ss=ngss)
            report["data"][arm][c.name] = {"pools": pool_metrics(c), "tx": transcript_metrics(c)}
            print(f"[{arm} eps={eps}] {c.name}  ({time.time()-t0:.0f}s)", flush=True)
    json.dump(report, open(args.out, "w"), indent=2)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
