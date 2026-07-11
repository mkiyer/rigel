"""Dissect the worst-error loci in a benchmark condition (3-pool: gDNA / nRNA / mRNA).

Ranks loci by total absolute per-transcript error |observed - expected| (the true error vs
ground truth) and decomposes each locus's error into the directed pool flows:

    gDNA -> nRNA   (gDNA mis-called as nascent — the gene-body unspliced twin)
    gDNA -> mRNA   (gDNA mis-called as a mature isoform)
    mRNA -> nRNA   (mature stolen by the nascent component — RNA->RNA, INVISIBLE to net_gdna_to_rna)
    isoform-swap   (mature mass redistributed among mature isoforms within the locus)

The mRNA->nRNA flow is the leak the gDNA-vs-RNA metric hides; it appears even in a pure-RNA
(gdna_none) library, so the nascent sink is an EM/prior artifact independent of gDNA.

Usage: python scripts/debug/dissect_loci.py [condition] [top_n]
"""
import sys

import pandas as pd

S = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
top_n = int(sys.argv[2]) if len(sys.argv) > 2 else 12

loc = pd.read_csv(f"{S}/net_flow_per_locus.tsv", sep="\t")
loc = loc[loc.condition == cond].set_index("locus_id")
tx = pd.read_csv(f"{S}/net_flow_per_transcript.tsv", sep="\t")
tx = tx[tx.condition == cond].copy()
tx["is_nrna"] = tx.transcript_id.str.startswith("RIGEL_NRNA_")


def locus_decomp(g):
    nr, mr = g[g.is_nrna], g[~g.is_nrna]
    gdna_to_nrna = nr.net_from_gdna.sum()
    nrna_absorbed = nr.observed.sum()  # truth nrna==0 in nrna_none ⇒ all artifact
    mrna_to_nrna = nrna_absorbed - gdna_to_nrna  # nascent mass that came from mature isoforms
    gdna_to_mrna = mr.net_from_gdna.sum()
    # mature gained from gdna (gdna_to_mrna) and lost to nascent; the rest of mature's churn is swap.
    iso_swap = mr.net_from_rna_isoforms.clip(lower=0).sum()
    return pd.Series({
        "abs_err": g.delta.abs().sum(),
        "gdna_to_nrna": gdna_to_nrna, "gdna_to_mrna": gdna_to_mrna,
        "mrna_to_nrna": mrna_to_nrna, "iso_swap": iso_swap,
        "nrna_absorbed": nrna_absorbed, "n_nrna": int(g.is_nrna.sum()),
    })


d = tx.groupby("locus_id").apply(locus_decomp, include_groups=False)
d = d.join(loc[["gdna_prior_count", "rna_prior_count", "gdna_eff_len_em",
                "locus_span_bp", "n_em_fragments", "n_transcripts"]])

T = d.sum()
print(f"=== {cond} ===")
print(f"TOTAL  abs_err={T.abs_err:,.0f} | gDNA->nRNA {T.gdna_to_nrna:,.0f}  "
      f"gDNA->mRNA {T.gdna_to_mrna:,.0f}  mRNA->nRNA {T.mrna_to_nrna:,.0f}  "
      f"iso_swap {T.iso_swap:,.0f} | nrna_absorbed(artifact) {T.nrna_absorbed:,.0f}")
print(f"loci={len(d)}  top-{top_n} hold {100*d.nlargest(top_n,'abs_err').abs_err.sum()/T.abs_err:.0f}% of abs_err\n")

hdr = (f"{'locus':>5} {'ntx':>3} {'abs_err':>8} {'g->nr':>7} {'g->mr':>7} {'mr->nr':>7} "
       f"{'swap':>7} {'nr_abs':>7} | {'gdna_pri':>9} {'rna_pri':>9} {'efflen':>7} {'span':>7}")
print(hdr)
for lid, r in d.nlargest(top_n, "abs_err").iterrows():
    print(f"{int(lid):>5} {int(r.n_transcripts):>3} {r.abs_err:>8,.0f} {r.gdna_to_nrna:>7,.0f} "
          f"{r.gdna_to_mrna:>7,.0f} {r.mrna_to_nrna:>7,.0f} {r.iso_swap:>7,.0f} {r.nrna_absorbed:>7,.0f} | "
          f"{r.gdna_prior_count:>9,.0f} {r.rna_prior_count:>9,.0f} {r.gdna_eff_len_em:>7,.0f} "
          f"{r.locus_span_bp:>7,.0f}")

print("\n\n===== per-transcript drill-down (worst loci) =====")
for lid in d.nlargest(min(4, top_n), "abs_err").index:
    g = tx[tx.locus_id == lid].sort_values("observed", ascending=False)
    r = d.loc[lid]
    print(f"\n--- locus {int(lid)}: abs_err={r.abs_err:,.0f}  g->nr={r.gdna_to_nrna:,.0f} "
          f"g->mr={r.gdna_to_mrna:,.0f} mr->nr={r.mrna_to_nrna:,.0f} | "
          f"gdna_pri={r.gdna_prior_count:,.0f} rna_pri={r.rna_prior_count:,.0f} "
          f"efflen={r.gdna_eff_len_em:,.0f} span={r.locus_span_bp:,.0f} ---")
    print(f"  {'transcript':>40} {'kind':>5} {'nex':>3} {'splen':>6} "
          f"{'expect':>8} {'observe':>8} {'<-gdna':>7} {'<-iso':>7}")
    for _, t in g.iterrows():
        kind = "nrna" if t.is_nrna else ("mono" if t.single_exon else "mrna")
        nm = t.transcript_id if len(t.transcript_id) <= 40 else t.transcript_id[-40:]
        print(f"  {nm:>40} {kind:>5} {t.n_exons:>3.0f} {t.spliced_length:>6.0f} "
              f"{t.expected:>8,.0f} {t.observed:>8,.0f} {t.net_from_gdna:>7,.0f} "
              f"{t.net_from_rna_isoforms:>7,.0f}")
