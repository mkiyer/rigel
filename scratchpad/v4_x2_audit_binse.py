import numpy as np
from scipy.special import polygamma
_T=1e-12
A=np.load("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/df34908c-f6e7-45d0-827a-31f2a818e9d9/scratchpad/a2.npz")
X=np.load("/Users/mkiyer/proj/rigel/scratchpad/x2_poisson_edges.npz")
CONDS=["gdna_gdna300_ss_0.99_nrna_present_capture_off","gdna_gdna300_ss_0.50_nrna_present_capture_off",
       "gdna_gdna100_ss_0.50_nrna_present_capture_off","gdna_gdna300_ss_0.50_nrna_none_capture_off"]
psi1=lambda n: polygamma(1,np.maximum(np.asarray(n,float),_T))

# ---- 1. MC null run through x2's OWN (mass-as-count) budget ---------------------------------------
print("1. MC NULL through x2's OWN budget (mass used as the count). Unbiased => 0.")
for c in CONDS:
    d={k:A[f"{c}__{k}"] for k in ("n_spl_cnt","n_unspl_R_bnd","c_bnd","n_R_exon","c_exon","n_g_bnd","n_g_exon","w_nu","w_mu")}
    lam_mu=np.maximum(d["n_spl_cnt"],1e-6); lam_nu=np.maximum(d["n_unspl_R_bnd"]*d["c_bnd"],0.0)
    lam_R=np.maximum(d["n_R_exon"]*d["c_exon"],1e-6); lam_gb=np.maximum(d["n_g_bnd"]*d["c_bnd"],1e-6)
    lam_ge=np.maximum(d["n_g_exon"]*d["c_exon"],1e-6)
    sm=np.maximum(d["n_spl_cnt"]/np.maximum(A[f"{c}__n_spl_mass"],_T),1.0)  # count-per-mass at the junction
    k_nu=np.where(lam_nu>0,d["w_nu"]/np.maximum(lam_nu,_T),0.0); k_mu=d["w_mu"]/lam_mu
    rng=np.random.default_rng(1); res=[]
    for _ in range(300):
        Nmu=rng.poisson(lam_mu);Nnu=rng.poisson(lam_nu);NR=rng.poisson(lam_R)
        Ngb=rng.poisson(lam_gb);Nge=rng.poisson(lam_ge)
        num=k_nu*Nnu+k_mu*Nmu; ok=(Nmu>0)&(NR>0)&(Ngb>0)&(Nge>0)&(num>0)
        lp=np.log(num[ok]/NR[ok]*Nge[ok]/Ngb[ok])
        # x2's budget: every count replaced by its MASS  (spl: N/sm ; unspl/gdna: N/c)
        mmu=Nmu[ok]/sm[ok]; mnu=Nnu[ok]/d["c_bnd"][ok]; mR=NR[ok]/d["c_exon"][ok]
        mgb=Ngb[ok]/d["c_bnd"][ok]; mge=Nge[ok]/d["c_exon"][ok]
        wn=np.where(num[ok]>0,k_nu[ok]*Nnu[ok]/num[ok],0.0); wm=1.0-wn
        vp=np.where(wn>0,wn**2/np.maximum(mnu,_T),0.0)+wm**2/np.maximum(mmu,_T)+psi1(mR)+psi1(mgb)+psi1(mge)
        res.append(np.var(lp)-np.mean(vp))
    res=np.array(res); print(f"   {c[5:]:<48} mean={res.mean():>8.4f} sd={res.std():.4f}  (0 = unbiased)")

# ---- 2. bin-level bootstrap SE for the 'FLAT in extrapolation' claim ------------------------------
def cat(k,src=X): return np.concatenate([src[f"{c}__{k}"] for c in CONDS])
lp=cat("log_phi"); vp=cat("v_tot_ex"); Lx=cat("len_exon"); es=cat("eff_spl")
extrap=Lx/es
rng=np.random.default_rng(0)
def boot(m,R=2000):
    a,b=lp[m],vp[m]; n=a.size
    s=np.array([np.var(a[i])-np.mean(b[i]) for i in (rng.integers(0,n,(R,n)))])
    return np.var(a)-np.mean(b), s.std()
print("\n2. bin-level RESIDUAL +- bootstrap SE  (report calls extrapolation 'FLAT, 1.13x range')")
print("   by len(exon)/eff_spl:")
for a,b in ((-np.inf,4),(4,8),(8,16),(16,np.inf)):
    m=(extrap>=a)&(extrap<b); r,se=boot(m)
    print(f"     {a:g}-{b:g}: n={m.sum():5d} RESID={r:.4f} +- {se:.4f}")
print("   by exon length bp:")
for a,b in ((-np.inf,400),(400,700),(700,1200),(1200,np.inf)):
    m=(Lx>=a)&(Lx<b); r,se=boot(m)
    print(f"     {a:g}-{b:g}: n={m.sum():5d} RESID={r:.4f} +- {se:.4f}")
print("   by junction spliced MASS (x2's table 2):")
ns=cat("n_spl")
for a,b in ((-np.inf,30),(30,100),(100,300),(300,1000),(1000,np.inf)):
    m=(ns>=a)&(ns<b); r,se=boot(m)
    print(f"     {a:g}-{b:g}: n={m.sum():5d} RESID={r:.4f} +- {se:.4f}")

# ---- 3. 'well-counted' threshold sensitivity ------------------------------------------------------
print("\n3. the 'all five counts >= 30' threshold (report's cleanest estimate 0.1119). Sensitivity:")
five=["n_spl","n_unspl_R_bnd","n_R_exon","n_g_bnd","n_g_exon"]
FF={k:cat(k) for k in five}
for thr in (5,10,20,30,50,100,200):
    ok=np.ones(lp.size,bool)
    for k in five: ok&=FF[k]>=thr
    if ok.sum()<10: print(f"   thr={thr:4d}: n={ok.sum()} too few"); continue
    print(f"   thr={thr:4d}: n={ok.sum():5d}  Var={np.var(lp[ok]):.4f}  E[vp]={np.mean(vp[ok]):.4f}  RESID={np.var(lp[ok])-np.mean(vp[ok]):.4f}")

# ---- 4. tail concentration reproduce -------------------------------------------------------------
print("\n4. tail concentration (report: top1% 32-39%, top5% 68-77%)")
for c in CONDS:
    a=X[f"{c}__log_phi"]; dev=(a-a.mean())**2; o=np.sort(dev)[::-1]
    print(f"   {c[5:]:<48} top1%={o[:max(1,o.size//100)].sum()/o.sum():.3f} top5%={o[:max(1,o.size//20)].sum()/o.sum():.3f}")

# ---- 5. is the rho_mu cut a relabel of the n_spl cut? ---------------------------------------------
print("\n5. rho_mu bins vs n_spl bins: identical membership?")
rm=cat("rho_mu")
for (a,b),(a2,b2) in ((( -np.inf,0.3),(-np.inf,30)),((0.3,1.5),(30,150))):
    m1=(rm>=a)&(rm<b); m2=(ns>=a2)&(ns<b2)
    print(f"     rho_mu[{a:g},{b:g}) n={m1.sum()}   n_spl[{a2:g},{b2:g}) n={m2.sum()}   identical={np.array_equal(m1,m2)}")
