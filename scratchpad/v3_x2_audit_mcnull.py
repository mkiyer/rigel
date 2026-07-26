"""MC NULL TEST of x2's estimator: if the graft premise holds EXACTLY (true phi == 1 at every edge),
does RESIDUAL = Var(log phi) - E[v_pois] return 0?   And if a known premise variance s2 is injected,
does it return s2?   Real per-edge count geometry, five independent Poissons, x2's own filters."""
import numpy as np
from scipy.special import polygamma
_T=1e-12
D=np.load("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/df34908c-f6e7-45d0-827a-31f2a818e9d9/scratchpad/a2.npz")
CONDS=["gdna_gdna300_ss_0.99_nrna_present_capture_off","gdna_gdna300_ss_0.50_nrna_present_capture_off",
       "gdna_gdna100_ss_0.50_nrna_present_capture_off","gdna_gdna300_ss_0.50_nrna_none_capture_off"]
psi1=lambda n: polygamma(1,np.maximum(np.asarray(n,float),_T))

def run(cond, s2_true, R=400, use_counts=True, seed=0):
    d={k:D[f"{cond}__{k}"] for k in ("n_spl_cnt","n_spl_mass","n_unspl_R_bnd","c_bnd","n_R_exon","c_exon",
                                     "n_g_bnd","n_g_exon","w_nu","w_mu")}
    # TRUE Poisson means = the observed counts (counts, not masses)
    lam_mu=np.maximum(d["n_spl_cnt"],1e-6)
    lam_nu=np.maximum(d["n_unspl_R_bnd"]*d["c_bnd"],0.0)
    lam_R =np.maximum(d["n_R_exon"]*d["c_exon"],1e-6)
    lam_gb=np.maximum(d["n_g_bnd"]*d["c_bnd"],1e-6)
    lam_ge=np.maximum(d["n_g_exon"]*d["c_exon"],1e-6)
    # constants k so that the TRUE phi == 1 exactly:  phi = (k_nu N_nu + k_mu N_mu)/(k_R N_R) * N_ge/N_gb * C
    # choose k_nu,k_mu to reproduce the observed numerator shares w_nu,w_mu; then set C to force phi_true=1.
    k_nu=np.where(lam_nu>0, d["w_nu"]/np.maximum(lam_nu,_T), 0.0); k_mu=d["w_mu"]/lam_mu
    num_t=k_nu*lam_nu+k_mu*lam_mu                       # == 1 by construction
    C=num_t/ (1.0/1.0) * lam_R*lam_gb/lam_ge            # so that phi_true = num_t/(C^-1 ...) = 1
    rng=np.random.default_rng(seed); out=[]
    for _ in range(R):
        Nmu=rng.poisson(lam_mu); Nnu=rng.poisson(lam_nu); NR=rng.poisson(lam_R)
        Ngb=rng.poisson(lam_gb); Nge=rng.poisson(lam_ge)
        ok=(Nmu>0)&(NR>0)&(Ngb>0)&(Nge>0)
        num=k_nu*Nnu+k_mu*Nmu
        ok&= num>0
        lt=np.where(s2_true>0, rng.normal(0.0,np.sqrt(s2_true),size=lam_mu.size),0.0)  # injected premise error
        phi=num/np.maximum(NR,1)*np.maximum(Nge,1)/np.maximum(Ngb,1)*C/np.maximum(num_t,_T)*np.exp(lt)
        lp=np.log(np.maximum(phi[ok],_T))
        wn=np.where(num[ok]>0,k_nu[ok]*Nnu[ok]/num[ok],0.0); wm=1.0-wn
        vnum=np.where(wn>0,wn**2/np.maximum(Nnu[ok],_T),0.0)+wm**2/np.maximum(Nmu[ok],_T)
        vp=vnum+psi1(NR[ok])+psi1(Ngb[ok])+psi1(Nge[ok])
        out.append(np.var(lp)-np.mean(vp))
    return np.array(out)

print("MC NULL — true phi == 1 at every edge (premise EXACT). An unbiased estimator returns 0.")
print(f"{'condition':<48}{'mean RESID':>12}{'sd':>9}{'p2.5':>9}{'p97.5':>9}")
for c in CONDS:
    r=run(c,0.0); print(f"{c[5:]:<48}{r.mean():>12.4f}{r.std():>9.4f}{np.percentile(r,2.5):>9.4f}{np.percentile(r,97.5):>9.4f}")
print("\nMC with an INJECTED premise variance s2 = 0.50 (the corrected all-edge residual).")
print(f"{'condition':<48}{'mean RESID':>12}{'sd':>9}{'bias':>9}")
for c in CONDS:
    r=run(c,0.50); print(f"{c[5:]:<48}{r.mean():>12.4f}{r.std():>9.4f}{r.mean()-0.50:>9.4f}")
print("\nMC with s2 = 0.10 (the 'well-counted' residual):")
for c in CONDS:
    r=run(c,0.10); print(f"{c[5:]:<48}{r.mean():>12.4f}{r.std():>9.4f}{r.mean()-0.10:>9.4f}")
