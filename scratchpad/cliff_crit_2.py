"""CRITIC #2 — completeness: does a PRECISION-only cross-cliff fix recover node 1909,
or is a MODE fix also needed? The message MODE is stretched-wrong (f_pos=0.718 onto a
~1.5%-RNA node). We test whether crushing the message precision alone lets the node's own
weak-but-correct belief (fg_loc=0.953, tau=1.6) win, and which derivation's delivered
precision actually lands f_g ~ 0.95.

Faithful 2-simplex psi-style fusion in (lambda=logit f_g, tilt s):
  f_g = sigmoid(lambda), f_R = 1-f_g, f_pos = f_R*s, f_neg = f_R*(1-s).
Objective = own belief on lambda (the strand solve, mode logit(0.953), prec tau_own=1.6)
          + gDNA measurement on log f_g (mode mo_g, prec cm_g)
          + RNA+ measurement on log f_pos (mode log 0.718, prec cm_p)   <-- WRONG mode
          + RNA- measurement on log f_neg (mode mo_n, prec cm_n)
          + composition message on lambda (mode lam_msg, prec c_tau).
Posterior mean f_g by grid integration over (lambda, s).
"""
import numpy as np

EPS = 1e-12

# ---- anchor node 1909 numbers (from HANDOFF_4 sec6 + the derivations) ----
FG_OWN   = 0.953        # own message-free solve (CORRECT, weak)
TAU_OWN  = 1.6          # own precision on lambda=logit(f_g)  (the strand tilt)
FG_ORACLE= 0.985
LAM_OWN  = np.log(FG_OWN/(1-FG_OWN))     # logit(0.953)=3.006

# corrupting RNA+ measurement (the reported bug channel)
CM_P_BUG = 26.0
MO_P_BUG = np.log(0.718)                 # f_pos=0.718 stretched-wrong
# gDNA measurement (roughly correct, continues from the mostly-gDNA exon)
CM_G_BUG = 25.0
MO_G_BUG = np.log(0.90)                  # gDNA message ~ high f_g (deriv3: ~1 nat off truth 0.985)
# RNA- measurement: negligible (sense-stranded), keep small precision, low mode
CM_N_BUG = 2.0
MO_N_BUG = np.log(0.01)
# composition message c_tau: RNA-heavy neighbour lambda (also corrupting per ablation)
C_TAU_BUG = 12.0
LAM_MSG   = np.log(0.36/(1-0.36))        # logit(0.36) src RNA-heavy ~ -0.575

# cliff quantities
R      = 407.0
LOGR2  = np.log(R)**2                     # 36.6
N_SPL  = 47.0

# grid
LAM = np.linspace(-6.0, 10.0, 801)
S   = np.linspace(1e-3, 1-1e-3, 401)
LG, SG = np.meshgrid(LAM, S, indexing='ij')
FG = 1.0/(1.0+np.exp(-LG))
FR = 1.0-FG
FP = FR*SG
FN = FR*(1.0-SG)

def fused_fg(cm_g, mo_g, cm_p, mo_p, cm_n, mo_n, c_tau, lam_msg,
             tau_own=TAU_OWN, lam_own=LAM_OWN):
    L  = -0.5*tau_own*(LG-lam_own)**2
    L += -0.5*cm_g*(np.log(np.maximum(FG,EPS))-mo_g)**2
    L += -0.5*cm_p*(np.log(np.maximum(FP,EPS))-mo_p)**2
    L += -0.5*cm_n*(np.log(np.maximum(FN,EPS))-mo_n)**2
    L += -0.5*c_tau*(LG-lam_msg)**2
    L -= L.max()
    W = np.exp(L)
    Z = W.sum()
    return float((W*FG).sum()/Z)

# ---- 0. sanity: with NO messages, own belief alone -> ~0.953 ----
fg_own_only = fused_fg(0,0,0,0,0,0,0,0)
print(f"[0] own belief alone (no messages)            f_g = {fg_own_only:.3f}  (target ~0.953)")

# ---- 1. reproduce the BUG (full-precision wrong messages) ----
fg_bug = fused_fg(CM_G_BUG, MO_G_BUG, CM_P_BUG, MO_P_BUG, CM_N_BUG, MO_N_BUG, C_TAU_BUG, LAM_MSG)
print(f"[1] BUG (full-precision wrong cm_p+c_tau)      f_g = {fg_bug:.3f}  (reported collapse ~0.51)")

print()
print("--- sweep: DAMP BOTH corrupting streams (cm_p AND c_tau) by a common factor ---")
print(" the ablation says BOTH must be damped. sweep delivered precision, mode LEFT WRONG:")
for cmp_new, ctau_new, tag in [
    (26.0, 12.0, "none (bug)"),
    (5.0,  5.0,  "mild"),
    (1.0,  1.0,  "1.0"),
    (0.5,  0.5,  "0.5"),
    (0.1,  0.1,  "0.1"),
    (0.03, 0.03, "0.03"),
    (0.01, 0.01, "0.01"),
    (0.0,  0.0,  "0 (drop both)"),
]:
    fg = fused_fg(CM_G_BUG, MO_G_BUG, cmp_new, MO_P_BUG, CM_N_BUG, MO_N_BUG, ctau_new, LAM_MSG)
    print(f"   cm_p={cmp_new:6.3f} c_tau={ctau_new:6.3f} [{tag:13s}]  f_g = {fg:.3f}")

print()
print("--- each DERIVATION's delivered precision on the corrupting streams (mode LEFT WRONG) ---")

# Deriv #1: tau_msg = 1/(1/tau_src + (log r)^2), fed to ALL streams
cmp_d1  = 1.0/(1.0/CM_P_BUG + LOGR2)
ctau_d1 = 1.0/(1.0/C_TAU_BUG + LOGR2)
fg_d1 = fused_fg(CM_G_BUG, MO_G_BUG, cmp_d1, MO_P_BUG, CM_N_BUG, MO_N_BUG, ctau_d1, LAM_MSG)
print(f"[D1] (log r)^2 to all: cm_p={cmp_d1:.4f} c_tau={ctau_d1:.4f}  f_g = {fg_d1:.3f}")

# Deriv #2: cov = 1/max(1,r) on the MEASUREMENT stream (cm_p) only (minimal change); c_tau kept
cov = 1.0/max(1.0, R)
cmp_d2  = CM_P_BUG*cov
fg_d2_measonly = fused_fg(CM_G_BUG, MO_G_BUG, cmp_d2, MO_P_BUG, CM_N_BUG, MO_N_BUG, C_TAU_BUG, LAM_MSG)
print(f"[D2a] cov=1/r on cm_p ONLY (c_tau UNTOUCHED=12): cm_p={cmp_d2:.4f}  f_g = {fg_d2_measonly:.3f}")
# Deriv #2 extended: cov also on c_tau (its A/B fallback)
ctau_d2 = C_TAU_BUG*cov
fg_d2_both = fused_fg(CM_G_BUG, MO_G_BUG, cmp_d2, MO_P_BUG, CM_N_BUG, MO_N_BUG, ctau_d2, LAM_MSG)
print(f"[D2b] cov=1/r on BOTH cm_p and c_tau:            cm_p={cmp_d2:.4f} c_tau={ctau_d2:.4f}  f_g = {fg_d2_both:.3f}")

# Deriv #3: DerSimonian-Laird per component. Anchor gap on own self-solve.
# RNA+ arm: G = mo_p_msg - mo_p_own; v_src=1/cm_p; v_own,R=(f_g_own)^2/tau_own
mo_p_own = np.log(max((1-FG_OWN)*0.99, EPS))   # own f_pos ~ f_R*sense ~ 0.047*0.99
G_R = MO_P_BUG - mo_p_own
v_src_R  = 1.0/CM_P_BUG
v_own_R  = (FG_OWN**2)/TAU_OWN
b2_R = max(0.0, G_R**2 - v_src_R - v_own_R)
cmp_d3 = 1.0/(v_src_R + b2_R)
# composition arm: lambda gap
G_lam = LAM_MSG - LAM_OWN
v_src_lam = 1.0/C_TAU_BUG
v_own_lam = 1.0/TAU_OWN
b2_lam = max(0.0, G_lam**2 - v_src_lam - v_own_lam)
ctau_d3 = 1.0/(v_src_lam + b2_lam)
fg_d3 = fused_fg(CM_G_BUG, MO_G_BUG, cmp_d3, MO_P_BUG, CM_N_BUG, MO_N_BUG, ctau_d3, LAM_MSG)
print(f"[D3] DL: cm_p={cmp_d3:.4f} (G_R={G_R:.2f}) c_tau={ctau_d3:.4f} (G_lam={G_lam:.2f})  f_g = {fg_d3:.3f}")

print()
print("--- MODE-FIX comparison: correct the RNA+ mode to truth (f_pos~0.015) at FULL precision ---")
mo_p_true = np.log(0.015)
fg_modefix = fused_fg(CM_G_BUG, MO_G_BUG, CM_P_BUG, mo_p_true, CM_N_BUG, MO_N_BUG, C_TAU_BUG, np.log(0.985/0.015))
print(f"[MF] correct modes, full precision:            f_g = {fg_modefix:.3f}")

print()
print("==== SENSITIVITY: vary the gDNA-message mode/precision + own-belief precision ====")
for mo_g, cm_g in [(np.log(0.985),25.0),(np.log(0.90),25.0),(np.log(0.80),25.0),(np.log(0.90),10.0)]:
    fg_bug2 = fused_fg(cm_g, mo_g, CM_P_BUG, MO_P_BUG, CM_N_BUG, MO_N_BUG, C_TAU_BUG, LAM_MSG)
    # D2a measurement-only vs D1 all-stream
    fg_d2a = fused_fg(cm_g, mo_g, CM_P_BUG/R, MO_P_BUG, CM_N_BUG, MO_N_BUG, C_TAU_BUG, LAM_MSG)
    fg_d1x = fused_fg(cm_g, mo_g, 1/(1/26+LOGR2), MO_P_BUG, CM_N_BUG, MO_N_BUG, 1/(1/12+LOGR2), LAM_MSG)
    print(f" mo_g=log{np.exp(mo_g):.3f} cm_g={cm_g:4.1f}: bug={fg_bug2:.3f}  D2a(meas-only)={fg_d2a:.3f}  D1(both)={fg_d1x:.3f}")

print()
print("==== MATCHED-COMPOSITION BIG CLIFF (deriv3 case C): mode is CORRECT, r large ====")
print(" a message that AGREES with own (same composition) at full precision, does damping hurt?")
# message mode matches own: f_pos ~ own f_pos ~ 0.015; lambda ~ own lambda 3.0
mo_p_match = np.log(0.015); lam_match = LAM_OWN
fg_match_full  = fused_fg(CM_G_BUG, MO_G_BUG, 26.0, mo_p_match, CM_N_BUG, MO_N_BUG, 12.0, lam_match)
fg_match_damp1 = fused_fg(CM_G_BUG, MO_G_BUG, 1/(1/26+LOGR2), mo_p_match, CM_N_BUG, MO_N_BUG, 1/(1/12+LOGR2), lam_match)
print(f"   matched msg full precision  f_g = {fg_match_full:.3f}")
print(f"   matched msg (log r)^2 damp  f_g = {fg_match_damp1:.3f}   (damping a CORRECT msg -> harmless, own agrees)")

print()
print("==== AMBIG SOLE-INFO: tau_own=0, message is the ONLY info, r large mismatched ====")
print(" here precision-only 'recovery' = fall to no-info; there is NO own belief to defer to:")
for cmp_v,tag in [(26.0,'full (trust bad msg)'),(1/(1/26+LOGR2),'(log r)^2 damp')]:
    fg_amb = fused_fg(0.0, 0.0, cmp_v, MO_P_BUG, 0.0, 0.0, 0.0, 0.0, tau_own=1e-6, lam_own=0.0)
    print(f"   cm_p={cmp_v:.4f} [{tag:22s}]  f_g = {fg_amb:.3f}")
