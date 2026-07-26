"""FINAL correction: the ONLY count/mass error is the SPLICED arm (SP/SN are per-face MASS; the true
per-face COUNT is ~2.0x larger and is what bp_solver's v_mu uses).  Oracle pools are already counts."""
import numpy as np
from scipy.special import polygamma
_T=1e-12
A=np.load("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/df34908c-f6e7-45d0-827a-31f2a818e9d9/scratchpad/a2.npz")
CONDS=["gdna_gdna300_ss_0.99_nrna_present_capture_off","gdna_gdna300_ss_0.50_nrna_present_capture_off",
       "gdna_gdna100_ss_0.50_nrna_present_capture_off","gdna_gdna300_ss_0.50_nrna_none_capture_off"]
psi1=lambda n: polygamma(1,np.maximum(np.asarray(n,float),_T))
def bud(d,ns):
    return (np.where(d["w_nu"]>0,d["w_nu"]**2/np.maximum(d["n_unspl_R_bnd"],_T),0.0)
            + d["w_mu"]**2/np.maximum(ns,_T) + psi1(d["n_R_exon"])+psi1(d["n_g_bnd"])+psi1(d["n_g_exon"]))
D={c:{k:A[f"{c}__{k}"] for k in ("n_spl_cnt","n_spl_mass","n_unspl_R_bnd","n_R_exon","n_g_bnd","n_g_exon","w_nu","w_mu","phi","share")} for c in CONDS}
print("A) per-condition, EXACT form, spliced arm on the MASS (x2) vs on the true COUNT")
print(f"{'condition':<48}{'Var':>9}{'vp MASS':>9}{'RES MASS':>10}{'%exp':>7}{'vp CNT':>9}{'RES CNT':>10}{'%exp':>7}")
rm=[];rc=[]
for c in CONDS:
    d=D[c]; lp=np.log(d["phi"]); V=np.var(lp)
    bm=bud(d,d["n_spl_mass"]).mean(); bc=bud(d,d["n_spl_cnt"]).mean()
    rm.append(V-bm); rc.append(V-bc)
    print(f"{c[5:]:<48}{V:>9.4f}{bm:>9.4f}{V-bm:>10.4f}{100*bm/V:>6.1f}%{bc:>9.4f}{V-bc:>10.4f}{100*bc/V:>6.1f}%")
rm,rc=np.array(rm),np.array(rc)
print(f"   x2 (mass):   mean={rm.mean():.4f} sd={rm.std(ddof=1):.4f} ({100*rm.std(ddof=1)/rm.mean():.1f}%)   [report: 0.436 +- 0.021, 22-33% expl]")
print(f"   CORRECTED:   mean={rc.mean():.4f} sd={rc.std(ddof=1):.4f} ({100*rc.std(ddof=1)/rc.mean():.1f}%)")
LP=np.concatenate([np.log(D[c]["phi"]) for c in CONDS])
BM=np.concatenate([bud(D[c],D[c]["n_spl_mass"]) for c in CONDS]); BC=np.concatenate([bud(D[c],D[c]["n_spl_cnt"]) for c in CONDS])
print(f"   POOLED n={LP.size}: Var={np.var(LP):.4f}  RES mass={np.var(LP)-BM.mean():.4f} ({100*BM.mean()/np.var(LP):.1f}% expl)"
      f"   RES cnt={np.var(LP)-BC.mean():.4f} ({100*BC.mean()/np.var(LP):.1f}% expl)")

print("\nB) shape by junction spliced COUNT, corrected budget, with bootstrap SE")
NS=np.concatenate([D[c]["n_spl_cnt"] for c in CONDS]); rng=np.random.default_rng(0)
def boot(m,R=2000):
    a,b=LP[m],BC[m]; n=a.size
    s=np.array([np.var(a[i])-np.mean(b[i]) for i in rng.integers(0,n,(R,n))]); return np.var(a)-np.mean(b),s.std()
prev=None
for a,b in ((-np.inf,30),(30,100),(100,300),(300,1000),(1000,np.inf)):
    m=(NS>=a)&(NS<b); r,se=boot(m)
    print(f"   {a:g}-{b:g}: n={m.sum():5d} med={np.median(NS[m]):8.1f} Var={np.var(LP[m]):8.4f} E[vp]={BC[m].mean():7.4f} RESID={r:8.4f} +- {se:.4f}")
print("   -> report's table (2) claims a MONOTONE 23.5x fall (1.0447/0.4788/0.1905/0.1178/0.0444)")

print("\nC) MC NULL (true phi == 1 everywhere), correct Poisson means, both budgets")
print(f"{'condition':<48}{'RES (x2 mass budget)':>22}{'RES (count budget)':>20}")
for c in CONDS:
    d=D[c]
    lm=np.maximum(d["n_spl_cnt"],1e-6); ln=np.maximum(d["n_unspl_R_bnd"],0.0)
    lR=np.maximum(d["n_R_exon"],1e-6); lgb=np.maximum(d["n_g_bnd"],1e-6); lge=np.maximum(d["n_g_exon"],1e-6)
    sm=np.maximum(d["n_spl_cnt"]/np.maximum(d["n_spl_mass"],_T),1.0)
    knu=np.where(ln>0,d["w_nu"]/np.maximum(ln,_T),0.0); kmu=d["w_mu"]/lm
    numt=knu*ln+kmu*lm   # == 1
    rng2=np.random.default_rng(3); am=[];ac=[]
    for _ in range(400):
        Nm=rng2.poisson(lm);Nn=rng2.poisson(ln);NR=rng2.poisson(lR);Ngb=rng2.poisson(lgb);Nge=rng2.poisson(lge)
        num=knu*Nn+kmu*Nm; ok=(Nm>0)&(NR>0)&(Ngb>0)&(Nge>0)&(num>0)
        # phi_hat / phi_true , phi_true = numt/(lR^-1 ... ) built so that the ratio is exactly 1 at the means
        lp=np.log((num[ok]/numt[ok])/(NR[ok]/lR[ok])*(Nge[ok]/lge[ok])/(Ngb[ok]/lgb[ok]))
        wn=knu[ok]*Nn[ok]/num[ok]; wm=1.0-wn
        def B(ns): return np.where(wn>0,wn**2/np.maximum(Nn[ok],_T),0.0)+wm**2/np.maximum(ns,_T)+psi1(NR[ok])+psi1(Ngb[ok])+psi1(Nge[ok])
        am.append(np.var(lp)-np.mean(B(Nm[ok]/sm[ok]))); ac.append(np.var(lp)-np.mean(B(Nm[ok])))
    print(f"{c[5:]:<48}{np.mean(am):>22.4f}{np.mean(ac):>20.4f}")
