import numpy as np
D = np.load("/Users/mkiyer/proj/rigel/scratchpad/x2_poisson_edges.npz")
CONDS = ["gdna_gdna300_ss_0.99_nrna_present_capture_off",
         "gdna_gdna300_ss_0.50_nrna_present_capture_off",
         "gdna_gdna100_ss_0.50_nrna_present_capture_off",
         "gdna_gdna300_ss_0.50_nrna_none_capture_off",
         "gdna_gdna300_ss_0.99_nrna_present_capture_on",
         "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong"]
def g(c,k): return D[f"{c}__{k}"]
print("A. eff_spl: is it constant?  (report claims 99.91 at EVERY edge)")
for c in CONDS:
    e = g(c,"eff_spl")
    print(f"  {c[5:]:<48} n={e.size:5d} min={e.min():.4f} med={np.median(e):.4f} max={e.max():.4f} nuniq={np.unique(np.round(e,6)).size}")
print()
print("B. is rho_mu binning IDENTICAL to n_spl binning?  rho_mu*eff_spl - n_spl:")
for c in CONDS[:1]:
    print("   max|rho_mu*eff_spl - n_spl| =", np.abs(g(c,"rho_mu")*g(c,"eff_spl")-g(c,"n_spl")).max())
print()
print("C. n_spl is the accumulator MASS (SP_l/SP_r), not the COUNT. distribution:")
for c in CONDS:
    x=g(c,"n_spl"); print(f"  {c[5:]:<48} min={x.min():.4f} p1={np.percentile(x,1):.3f} p25={np.percentile(x,25):.2f} med={np.median(x):.2f} frac<1={np.mean(x<1):.3f}")
print()
print("D. n_Rs_exon identically zero?")
for c in CONDS:
    x=g(c,"n_Rs_exon"); print(f"  {c[5:]:<48} min={x.min()} max={x.max()} nonzero={np.count_nonzero(x)}")
print()
print("E. min of the other four counts (the 1/n and psi' arguments):")
for c in CONDS:
    row=[f"{k}:{g(c,k).min():.4g}" for k in ("n_unspl_R_bnd","n_R_exon","n_g_bnd","n_g_exon")]
    print(f"  {c[5:]:<48} "+"  ".join(row))
print()
print("F. v_num_ex tail: does 1/n on a small MASS blow up?")
for c in CONDS[:4]:
    v=g(c,"v_num_ex"); print(f"  {c[5:]:<48} max={v.max():.3f} p99={np.percentile(v,99):.3f} mean={v.mean():.4f}  frac of E[v_num] from top 1%: {np.sort(v)[-max(1,v.size//100):].sum()/v.sum():.3f}")
