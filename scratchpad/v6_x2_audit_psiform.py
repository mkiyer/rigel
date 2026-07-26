import numpy as np
from scipy.special import polygamma
_T=1e-12
A=np.load("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/df34908c-f6e7-45d0-827a-31f2a818e9d9/scratchpad/a2.npz")
CONDS=["gdna_gdna300_ss_0.99_nrna_present_capture_off","gdna_gdna300_ss_0.50_nrna_present_capture_off",
       "gdna_gdna100_ss_0.50_nrna_present_capture_off","gdna_gdna300_ss_0.50_nrna_none_capture_off"]
psi1=lambda n: polygamma(1,np.maximum(np.asarray(n,float),_T))
D={c:{k:A[f"{c}__{k}"] for k in ("n_spl_cnt","n_spl_mass","n_unspl_R_bnd","n_R_exon","n_g_bnd","n_g_exon","w_nu","w_mu","phi","share")} for c in CONDS}
def bud(d,ns,form):
    if form=="lin": v=np.where(d["w_nu"]>0,d["w_nu"]**2/np.maximum(d["n_unspl_R_bnd"],_T),0.0)+d["w_mu"]**2/np.maximum(ns,_T)
    else:           v=np.where(d["w_nu"]>0,d["w_nu"]**2*psi1(d["n_unspl_R_bnd"]),0.0)+d["w_mu"]**2*psi1(ns)
    return v+psi1(d["n_R_exon"])+psi1(d["n_g_bnd"])+psi1(d["n_g_exon"])
print("1) finding #2 re-tested with the TRUE spliced COUNT: is w^2*psi' still 'incoherent'?")
for form in ("lin","psi"):
    rs=[]
    for c in CONDS:
        d=D[c]; rs.append(np.var(np.log(d["phi"]))-bud(d,d["n_spl_cnt"],form).mean())
    rs=np.array(rs); print(f"   {form:>4}: {np.round(rs,4)}  spread max/min = {rs.max()/rs.min():.2f}x  mean={rs.mean():.4f}")
print("   (x2 reported, on the MASS: psi' -> 0.1527/0.3887/0.2691/0.3957 = 2.55x spread; lin -> 1.12x)")
for form in ("lin","psi"):
    rs=[]
    for c in CONDS:
        d=D[c]; rs.append(np.var(np.log(d["phi"]))-bud(d,d["n_spl_mass"],form).mean())
    rs=np.array(rs); print(f"   MASS {form:>4}: {np.round(rs,4)} spread {rs.max()/rs.min():.2f}x")

print("\n2) 'well-counted' threshold, CORRECTED counts (spliced=count, pools as-is)")
LP=np.concatenate([np.log(D[c]["phi"]) for c in CONDS])
BC=np.concatenate([bud(D[c],D[c]["n_spl_cnt"],"lin") for c in CONDS])
F={k:np.concatenate([D[c][k] for c in CONDS]) for k in ("n_spl_cnt","n_unspl_R_bnd","n_R_exon","n_g_bnd","n_g_exon")}
for thr in (5,10,20,30,50,100,200):
    ok=np.ones(LP.size,bool)
    for k in F: ok&=F[k]>=thr
    if ok.sum()<10: print(f"   thr={thr:4d}: n={ok.sum()} too few"); continue
    print(f"   thr={thr:4d}: n={ok.sum():5d} Var={np.var(LP[ok]):.4f} E[vp]={BC[ok].mean():.4f} RESID={np.var(LP[ok])-BC[ok].mean():.4f}")
print("   (x2's reported 'all five >= 30' -> 0.1119, but its n_spl is the MASS => really >=60 in count)")
