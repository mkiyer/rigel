"""Analyze the capture binding-energy x DNA:RNA sweep: leak (by pool) + effective lengths per condition.

Consumes the evaluate_suite outputs of efflen_binding_sweep. For each condition (parsed binding + DNA:RNA)
reports, END-TO-END:
  - net gDNA->RNA leak, split by pool (-> mature mRNA, -> nascent nRNA) via net_flow_per_transcript
  - gDNA effective length (median gdna_eff_len_em) + locus span  -> contraction gE/span
  - transcript effective length (median, mature) -> the exon-competition footprint
  - the KEY ratio gDNA_eff / transcript_eff  (should -> 1 as binding -> infinity if gDNA correctly
    contracts to the exon footprint; if it stays > 1, the gDNA is under-weighted -> leak)

    python scripts/debug/efflen_sweep_analysis.py
"""
import re
from pathlib import Path
import numpy as np, pandas as pd
pd.set_option("display.width", 240)

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/efflen_binding_sweep")
BIND = {"off": 0.0, "b2": 2.0, "b10": 10.0, "b50": 50.0, "b200": 200.0}
GDNA = {"gdna05": 0.05, "gdna100": 1.0, "gdna300": 3.0}


def parse(cond):
    b = re.search(r"capture_(\w+)$", cond); g = re.search(r"gdna_(gdna\d+)_", cond)
    return (GDNA.get(g.group(1)) if g else None, BIND.get(b.group(1)) if b else None)


def main():
    locf = SUITE / "net_flow_per_locus.tsv"
    txf = SUITE / "net_flow_per_transcript.tsv"
    if not locf.exists():
        print("no net_flow yet — run evaluate_suite first"); return
    loc = pd.read_csv(locf, sep="\t")
    tx = pd.read_csv(txf, sep="\t")
    tx["is_nrna"] = tx.transcript_id.str.startswith("RIGEL_NRNA_")
    rows = []
    for cond in sorted(loc.condition.unique()):
        gr, bind = parse(cond)
        if gr is None or bind is None:
            continue
        net = loc[loc.condition == cond].net_gdna_to_rna.sum()
        t = tx[tx.condition == cond]
        leak_mat = t[~t.is_nrna].net_from_gdna.sum()
        leak_nasc = t[t.is_nrna].net_from_gdna.sum()
        # eff-lengths from the per-condition rigel_out
        od = SUITE / cond / "rigel_out"
        gE = span = tE = np.nan
        lf = od / "loci.feather"
        if lf.exists():
            L = pd.read_feather(lf)
            gE = float(L.gdna_eff_len_em.median()) if "gdna_eff_len_em" in L else np.nan
            span = float(L.locus_span_bp.median()) if "locus_span_bp" in L else np.nan
        qf = od / "quant.feather"
        if qf.exists():
            Q = pd.read_feather(qf)
            effcol = next((c for c in ["em_effective_length", "effective_length_em", "effective_length"] if c in Q), None)
            if effcol and "is_nrna" in Q:
                tE = float(Q.loc[~Q.is_nrna.astype(bool), effcol].median())
            elif effcol:
                tE = float(Q[effcol].median())
        rows.append(dict(bind=bind, dna_rna=gr, net_leak=net, leak_mat=leak_mat, leak_nasc=leak_nasc,
                         gdna_eff=gE, tx_eff=tE, span=span,
                         gE_over_span=gE / span if span else np.nan,
                         gE_over_tE=gE / tE if tE else np.nan))
    df = pd.DataFrame(rows).sort_values(["dna_rna", "bind"])
    print(df.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
    print("\n=== pivot: net_leak by DNA:RNA (rows) x binding (cols) ===")
    print(df.pivot(index="dna_rna", columns="bind", values="net_leak").to_string(float_format=lambda x: f"{x:,.0f}"))
    print("\n=== pivot: gDNA_eff/tx_eff ratio (should ->1 at high binding if gDNA contracts to exon footprint) ===")
    print(df.pivot(index="dna_rna", columns="bind", values="gE_over_tE").to_string(float_format=lambda x: f"{x:,.2f}"))


if __name__ == "__main__":
    main()
