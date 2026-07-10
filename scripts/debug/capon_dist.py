"""Quantify the capON directional gDNA under-call (RNA over-attribution = the leak driver) distribution:
by strand class (AMBIG both-strand vs single-strand POS/NEG vs NONE) and mass concentration, at the
strand-only stage vs the FINAL calibration. Directional (non-cancelling) so it reflects the leak, not the
net (which cancels). Uses the pass_trace NPZ."""
import sys
import numpy as np

NPZ = sys.argv[1] if len(sys.argv) > 1 else \
    "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb/gdna_gdna300_ss_0.99_nrna_none_capture_on/pass_trace.npz"
d = np.load(NPZ, allow_pickle=True)
raw = d["reg_raw"].astype(float)
true_g = d["reg_true_g"].astype(float)
cls = d["reg_cls"]
reg_node = d["reg_node"].astype(int)           # chain-node index per region
fp = d["free_pos"].astype(bool)[reg_node]      # per-region (free_pos/neg are per chain node)
fn = d["free_neg"].astype(bool)[reg_node]
p1 = d["reg_p1_str"].astype(float)   # strand-only f_g
fin = d["reg_p2_fin"].astype(float)  # final f_g
cls = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cls])

# strand class
sc = np.where(fp & fn, "AMBIG", np.where(fp & ~fn, "POS", np.where(~fp & fn, "NEG", "NONE")))
exon = cls == "exon"


def directional(fg):
    """gDNA under-called (=RNA over-attributed, leaks to RNA) and over-called, per region, in mass units."""
    cal_g = fg * raw
    under = np.maximum(true_g - cal_g, 0.0)  # gDNA under -> RNA over (LEAK driver)
    over = np.maximum(cal_g - true_g, 0.0)
    return under, over


for stage, fg in [("strand-only", p1), ("FINAL(p2)", fin)]:
    under, over = directional(fg)
    print(f"\n===== {stage}: directional gDNA under-call (RNA over-attribution) =====")
    print(f"  ALL regions: under={under.sum():>10,.0f}  over={over.sum():>10,.0f}  net={under.sum()-over.sum():>+10,.0f}")
    print(f"  EXON only  : under={under[exon].sum():>10,.0f}  over={over[exon].sum():>10,.0f}  "
          f"net={under[exon].sum()-over[exon].sum():>+10,.0f}")
    print("  EXON under-call by strand class (share of exon under-call):")
    tot = under[exon].sum()
    for k in ["AMBIG", "POS", "NEG", "NONE"]:
        m = exon & (sc == k)
        if under[m].sum() > 0 or m.sum():
            print(f"    {k:6} n={m.sum():>4}  under={under[m].sum():>10,.0f} ({100*under[m].sum()/max(tot,1):4.0f}%)"
                  f"  mean_raw={raw[m].mean() if m.sum() else 0:>8,.0f}")

# mass concentration of the FINAL exon under-call
under, _ = directional(fin)
ue = under * exon
order = np.argsort(ue)[::-1]
cum = np.cumsum(ue[order]) / max(ue.sum(), 1)
print(f"\n===== FINAL exon under-call mass concentration (total={ue.sum():,.0f}) =====")
for topn in [5, 10, 20, 50, 100]:
    if topn <= len(order):
        print(f"  top {topn:>3} exons: {100*cum[topn-1]:4.0f}% of the under-call")
nnz = int((ue > 0).sum())
print(f"  ({nnz} exon regions have any under-call)")

# how much did messages+KDE correct, by strand class (strand-only under -> final under)
u1, _ = directional(p1)
uf, _ = directional(fin)
print("\n===== correction (strand-only under -> FINAL under) by exon strand class =====")
for k in ["AMBIG", "POS", "NEG", "NONE"]:
    m = exon & (sc == k)
    if m.sum():
        print(f"  {k:6} n={m.sum():>4}  strand-only under={u1[m].sum():>10,.0f} -> FINAL under={uf[m].sum():>10,.0f}"
              f"  (corrected {100*(1-uf[m].sum()/max(u1[m].sum(),1)):4.0f}%)")
