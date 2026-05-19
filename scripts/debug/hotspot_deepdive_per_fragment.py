"""Deep-dive per-fragment analysis of gDNA→RNA mispredictions in 5 hotspot regions.

For each of 5 hotspot windows, we:
  1. Identify the locus_id and the top assigned RNA target.
  2. Pull locus-level posterior mass (theta_m, theta_n, theta_g) from loci.feather.
  3. Pull target effective length from quant.feather / nrna_quant.feather.
  4. Compute the per-fragment posterior ratio (target_k vs gDNA) implied by the EM.
  5. Pull every gDNA-source FP fragment that lies in the window from annotated.bam
     and emit a per-fragment likelihood walkthrough.

Sources:
- truth labels: query-name flowcell  (H7MFFDSXY = gDNA, C6EL5ANXX = RNA)
- annotated BAM tags: ZW (winner posterior), ZF (winner txp idx), ZC (category),
  ZS (splice type), ZL (locus_id), ZT (target txp id), ZN (n_candidates).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pysam

RUN = Path("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/no_mm")
BAM = RUN / "annotated.bam"

# (chrom, win_start_1based, win_end, label) — picked for diversity
HOTSPOTS = [
    ("chr11", 65_500_001, 65_510_000, "MALAT1-adjacent nRNA, high RNA expression"),
    ("chr9", 76_700_001, 76_710_000, "nRNA chr9:76.7M, high RNA expression"),
    ("chr2", 178_560_001, 178_570_000, "giant 278-kb nRNA chr2:178M, low RNA expr"),
    ("chr1", 152_350_001, 152_360_000, "FLG2 mRNA, near-zero local RNA expr"),
    ("chr11", 64_260_001, 64_270_000, "nRNA chr11:64.26M, opp-strand FPs"),
]

LOG_FL_RATIO_DEFAULT = 0.0  # will be filled from summary


def load_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    loci = pd.read_feather(RUN / "loci.feather").set_index("locus_id")
    quant = pd.read_feather(RUN / "quant.feather")
    nrna = pd.read_feather(RUN / "nrna_quant.feather").set_index("nrna_id")
    return loci, quant, nrna


def gather_all_fp(hotspots: list[tuple[str, int, int, str]]) -> dict[tuple[str, int, int], pd.DataFrame]:
    """One sequential pass over the (unsorted) BAM; collect FP rows per window."""
    by_win: dict[tuple[str, int, int], list[dict]] = {(c, s, e): [] for c, s, e, _ in hotspots}
    bam = pysam.AlignmentFile(BAM, "rb", check_sq=False)
    n = 0
    for rec in bam:
        n += 1
        if n % 5_000_000 == 0:
            print(f"  ...scanned {n:,} records")
        if rec.is_secondary or rec.is_supplementary or rec.is_read2:
            continue
        if rec.is_unmapped:
            continue
        qn = rec.query_name
        parts = qn.split(":")
        fc = parts[2] if len(parts) >= 3 else qn
        if fc != "H7MFFDSXY":
            continue
        try:
            zf = rec.get_tag("ZF")
        except KeyError:
            continue
        # Predicted: gDNA if ZF<0; else read ZT to decide nRNA vs mRNA
        zt = rec.get_tag("ZT") if rec.has_tag("ZT") else ""
        pred = "gdna" if zf < 0 else ("nrna" if zt.startswith("RIGEL_NRNA") else "mrna")
        if pred == "gdna":
            continue  # not a FP
        chrom = rec.reference_name
        pos1 = rec.reference_start + 1
        # Find matching window
        for key in by_win:
            c, ws, we = key
            if c == chrom and ws <= pos1 <= we:
                by_win[key].append(dict(
                    qname=qn, chrom=chrom, pos=pos1, end=rec.reference_end,
                    is_reverse=rec.is_reverse, tlen=rec.template_length,
                    cigar=rec.cigarstring, mapq=rec.mapping_quality,
                    nm=rec.get_tag("NM") if rec.has_tag("NM") else 0,
                    zf=zf,
                    zw=rec.get_tag("ZW") if rec.has_tag("ZW") else np.nan,
                    zc=rec.get_tag("ZC") if rec.has_tag("ZC") else "",
                    zs=rec.get_tag("ZS") if rec.has_tag("ZS") else "",
                    zl=rec.get_tag("ZL") if rec.has_tag("ZL") else -1,
                    zn=rec.get_tag("ZN") if rec.has_tag("ZN") else 0,
                    zt=zt,
                    zr=rec.get_tag("ZR") if rec.has_tag("ZR") else "",
                    zg=rec.get_tag("ZG") if rec.has_tag("ZG") else "",
                    pred=pred,
                ))
                break
    bam.close()
    return {k: pd.DataFrame(v) for k, v in by_win.items()}


def locus_summary(loci: pd.DataFrame, lid: int) -> dict:
    if lid not in loci.index:
        return {}
    row = loci.loc[lid]
    return {
        "locus_id": lid,
        "locus_span_bp": int(row["locus_span_bp"]),
        "n_transcripts": int(row["n_transcripts"]),
        "n_nrna_entities": int(row["n_nrna_entities"]),
        "n_em_fragments": int(row["n_em_fragments"]),
        "mrna": float(row["mrna"]),
        "nrna": float(row["nrna"]),
        "gdna": float(row["gdna"]),
        "gdna_eff_len": float(row["gdna_eff_len"]),
        "gdna_prior": float(row["gdna_prior"]),
    }


def target_eff_len(quant: pd.DataFrame, nrna: pd.DataFrame, zt: str) -> tuple[float, float, str]:
    """Return (target_eff_len, target_count_em, pool)."""
    if not zt:
        return (np.nan, np.nan, "?")
    if zt.startswith("RIGEL_NRNA"):
        if zt in nrna.index:
            row = nrna.loc[zt]
            return (float(row["effective_length"]), float(row["count"]), "nRNA")
    # transcript ID
    m = quant[quant["transcript_id"] == zt]
    if len(m):
        row = m.iloc[0]
        pool = "nRNA" if bool(row["is_nrna"]) else "mRNA"
        return (float(row["effective_length"]), float(row["count"]), pool)
    return (np.nan, np.nan, "?")


def compute_ratio(theta_target: float, theta_g: float, target_eff_len_bp: float,
                  gdna_eff_len_bp: float, fl_rna_mean: float, fl_gdna_mean: float,
                  fl_p_rna: float, fl_p_gdna: float) -> dict:
    """
    Per-fragment posterior ratio for target vs gDNA, assuming an unspliced
    fragment fully contained in both candidate spans.

    P(frag | target) = (1 / target_eff_len) * fl_pdf_rna(L)
    P(frag | gDNA)  = (1 / (2 * gdna_eff_len)) * fl_pdf_gdna(L)
    posterior_target / posterior_gdna
       = (theta_target / theta_g) * (2 * gdna_eff_len / target_eff_len)
         * (fl_pdf_rna(L) / fl_pdf_gdna(L))
    """
    if theta_g <= 0 or target_eff_len_bp <= 0:
        return {"ratio": float("inf")}
    geom = (2.0 * gdna_eff_len_bp) / target_eff_len_bp
    theta_ratio = theta_target / theta_g
    fl_ratio = fl_p_rna / fl_p_gdna if fl_p_gdna > 0 else 1.0
    ratio = theta_ratio * geom * fl_ratio
    return {
        "theta_ratio": theta_ratio,
        "geom_ratio_2gL_over_tL": geom,
        "fl_ratio": fl_ratio,
        "posterior_ratio_target_vs_gdna": ratio,
        "posterior_target_share": ratio / (1.0 + ratio),
    }


def main() -> None:
    loci, quant, nrna = load_tables()
    # FL means from summary
    fl_rna_mean = 246.17
    fl_gdna_mean = 310.51

    out_root = Path("/Users/mkiyer/proj/rigel/results/vcap_hotspot_deepdive_2026-05-17")
    out_root.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    print("Scanning annotated BAM once for all hotspots...")
    by_win = gather_all_fp(HOTSPOTS)
    for chrom, ws, we, label in HOTSPOTS:
        print(f"\n{'='*100}\n>>> {chrom}:{ws:,}-{we:,}  [{label}]\n{'='*100}")
        fps = by_win[(chrom, ws, we)]
        # Total gDNA-source in window: from the previous analyzer output, we know
        # it; here we only count FPs (RNA-predicted gDNA-source).
        print(f"  gDNA-source called RNA (false pos):    {len(fps):,}")
        if len(fps) == 0:
            continue
        # locus aggregate (skip locus_id=-1 / unassigned)
        zl_valid = fps.loc[fps["zl"] >= 0, "zl"]
        if not len(zl_valid):
            print("  No valid locus_id among FPs — skipping locus math.")
            continue
        lid_top = int(zl_valid.mode().iloc[0])
        ls = locus_summary(loci, lid_top)
        if not ls:
            print(f"  locus_id {lid_top} not in loci.feather — skipping")
            continue
        # dominant target (skip empty ZT)
        zt_nonempty = fps.loc[fps["zt"].astype(str).str.len() > 1, "zt"]
        top_target = zt_nonempty.value_counts().idxmax() if len(zt_nonempty) else ""
        target_eff_len_bp, target_count, target_pool = target_eff_len(quant, nrna, top_target)
        # theta
        total_theta = ls["mrna"] + ls["nrna"] + ls["gdna"]
        print(f"  Locus_id={lid_top}: span={ls['locus_span_bp']:,} bp  n_tx={ls['n_transcripts']}  "
              f"n_nrna={ls['n_nrna_entities']}  n_em_frag={ls['n_em_fragments']:,}")
        print(f"  EM mass: mRNA={ls['mrna']:.0f}  nRNA={ls['nrna']:.0f}  gDNA={ls['gdna']:.0f}  "
              f"(theta_g={ls['gdna']/total_theta:.3f})")
        print(f"  Locus gdna_eff_len = {ls['gdna_eff_len']:,.0f} bp  "
              f"(prior eta_g={ls['gdna_prior']:.4f})")
        print(f"  Top target:  {top_target}  pool={target_pool}  "
              f"eff_len={target_eff_len_bp:,.0f} bp  count_em={target_count:,.1f}")
        # Geom ratio (worst case: assume all theta is in the top target → upper bound on FP pull)
        # Conservative estimate: distribute nrna mass equally over n_nrna entities
        if target_pool == "nRNA":
            theta_target_local = ls["nrna"] / max(ls["n_nrna_entities"], 1)
        else:
            theta_target_local = ls["mrna"] / max(ls["n_transcripts"] - ls["n_nrna_entities"], 1)
        # FL ratio for FL=200 (typical paired-end short read insert): approximate
        # using Gaussian pdf w/ std=88
        L = 250  # typical insert
        fl_p_rna = np.exp(-((L - fl_rna_mean) / 88.0) ** 2 / 2) / (88.0 * np.sqrt(2 * np.pi))
        fl_p_gdna = np.exp(-((L - fl_gdna_mean) / 88.0) ** 2 / 2) / (88.0 * np.sqrt(2 * np.pi))
        r = compute_ratio(theta_target_local, ls["gdna"], target_eff_len_bp,
                          ls["gdna_eff_len"], fl_rna_mean, fl_gdna_mean, fl_p_rna, fl_p_gdna)
        print(f"  --- per-fragment posterior (insert L=250, target=dominant nRNA/mRNA) ---")
        print(f"  theta_target/theta_g                  = {r['theta_ratio']:.4f}")
        print(f"  2 * gdna_eff_len / target_eff_len     = {r['geom_ratio_2gL_over_tL']:.2f}")
        print(f"  FL_pdf_rna / FL_pdf_gdna at L=250     = {r['fl_ratio']:.3f}")
        print(f"  --> P(target) / P(gDNA)               = {r['posterior_ratio_target_vs_gdna']:.2f}")
        print(f"  --> P(target | frag)                  = {r['posterior_target_share']:.4f}")
        # ZW distribution
        print(f"  Observed ZW for FP fragments: median={fps['zw'].median():.3f}  "
              f"mean={fps['zw'].mean():.3f}  frac>=0.9={(fps['zw']>=0.9).mean():.3f}")
        # Sample 5 fragments
        print(f"\n  Sample 5 FP fragments:")
        for _, r in fps.head(5).iterrows():
            print(f"    {r['chrom']}:{r['pos']:,}-{r['end']:,}  "
                  f"{'rev' if r['is_reverse'] else 'fwd'}  tlen={r['tlen']}  "
                  f"ZW={r['zw']:.3f}  ZC={r['zc']}  ZS={r['zs']}  pred={r['pred']}  "
                  f"target={r['zt']}")
        # Save
        out = out_root / f"{chrom}_{ws}_{we}_fp.tsv"
        fps.to_csv(out, sep="\t", index=False)
        summary_rows.append({
            "chrom": chrom, "window": f"{ws}-{we}", "label": label,
            "fps": len(fps), "locus_id": lid_top, "locus_span_bp": ls["locus_span_bp"],
            "n_nrna_entities": ls["n_nrna_entities"], "theta_g": ls["gdna"]/total_theta,
            "gdna_eff_len_locus": ls["gdna_eff_len"], "top_target": top_target,
            "target_eff_len": target_eff_len_bp, "target_count_em": target_count,
            "theta_target_per_entity": theta_target_local,
            "geom_ratio_2gL_over_tL": r.get("geom_ratio_2gL_over_tL", np.nan),
            "posterior_target_share_at_L250": r.get("posterior_target_share", np.nan),
            "median_ZW": fps["zw"].median(),
        })
    pd.DataFrame(summary_rows).to_csv(out_root / "summary.tsv", sep="\t", index=False)
    print(f"\nWrote summary + per-window TSVs to {out_root}")


if __name__ == "__main__":
    main()
