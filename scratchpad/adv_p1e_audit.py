"""ADVERSARIAL audit of P1e's DL null.

Three questions the headline does not answer:

 Q1  The conservation identity ``Σ_c ρ_c E_c = M`` holds only for a claim that supplies EVERY component.
     A PARTIAL claim knows only the INEQUALITY ``Σ_sup ρ_c E_c ≤ M`` — "I did not explain all your mass"
     is not a contradiction. But ``b̂² = max(0, δ² − ...)`` is SYMMETRIC in δ. How much of the term's
     firing mass sits on the UNDER-claim side (δ > 0), where there is nothing to contradict?

 Q2  ``αᵀΣα`` charges the unsupplied components ZERO variance (``s_c = 0``) while they contribute
     ``1 − α_sup`` of the mass budget S from the node's OWN self-solve. That own density has its own
     error — on unstranded nodes it is the ``τ_own = 0`` default, i.e. INFINITE variance — and the MoM
     null does not subtract it. The project's own rule: "a MoM fit must subtract EVERY noise source the
     model already knows."  How much of S is own-fill, on the stratum that carries the effect?

 Q3  is ``αᵀΣα`` even the quadratic form of the Σ the code writes down?  With the unsupplied rows zeroed
     the restricted quadratic form is ``σ_cm²·α_sup² + Σ α_c² w_c``; the code computes
     ``σ_cm²·α_sup + Σ α_c² w_c`` (it dots ``Σα`` built with ``Σ_all α = 1`` against ``α`` over the
     supplied set only). Size the gap.

    OMP_NUM_THREADS=1 python scratchpad/adv_p1e_audit.py
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import p1e_lib as L  # noqa: E402

ORDER = [
    "gdna300_ss0.99_present_capOFF",
    "gdna300_ss0.50_present_capOFF",
    "gdna300_ss0.99_present_capON",
    "gdna100_ss0.50_present_VERYSTRONG",
    "none_ss0.50_present_capOFF",
    "gdna100_ss0.50_none_capOFF",
    "gdna300_ss0.50_none_capOFF",
]

rows = []
for name in ORDER:
    inp, dbg = L.solve(L.CONDS[name])
    t, nf = L.message_table(inp, dbg)
    live = t["nsup"] > 0
    T = {k: v[live] for k, v in t.items() if isinstance(v, np.ndarray)}
    T["cls"] = t["cls"][live]
    n = T["delta"].size

    # α_sup: the share of the claim's mass budget that the message actually SUPPLIED.
    a_sup = np.where(T["sup_g"], T["alpha_g"], 0.0) + np.where(
        T["sup_r"], T["alpha_p"] + T["alpha_n"], 0.0
    )
    own_fill = 1.0 - a_sup
    # supplied-only admissibility: does the supplied part ALONE exceed the node's mass?
    a_sup_over = a_sup * (T["M"] / np.maximum(T["M"], 1e-9)) * (T["S"] / np.maximum(T["M"], 1e-9))

    fires = T["bhat2"] > 0.0
    under = T["delta"] > 0.0  # S < M  → the claim did not explain all the mass
    # the term's total "firing mass": b̂² weighted by the destination mass it damps
    fm = T["bhat2"] * T["M"]
    tot = fm.sum()

    one = T["nsup"] == 1
    print(f"\n{'=' * 118}\n{name}   live={n}\n{'=' * 118}")
    print(f"  fires (b̂²>0): {fires.mean() * 100:5.1f}%   of which UNDER-claim (δ>0): "
          f"{(fires & under).sum() / max(fires.sum(), 1) * 100:5.1f}%")
    print(f"  b̂²·M share on UNDER-claim messages: {fm[under].sum() / max(tot, 1e-9) * 100:5.1f}%")
    print(f"  b̂²·M share on ONE-COMPONENT msgs  : {fm[one].sum() / max(tot, 1e-9) * 100:5.1f}%")
    print(f"  b̂²·M share on 1-comp × UNDER      : {fm[one & under].sum() / max(tot, 1e-9) * 100:5.1f}%")
    print(f"  supplied-alone inadmissible (Σ_sup m_c > M): all {100 * np.mean(a_sup_over > 1):.1f}%  "
          f"| fires {100 * np.mean(a_sup_over[fires] > 1):.1f}%  "
          f"| 1-comp&fires {100 * np.mean(a_sup_over[one & fires] > 1) if (one & fires).any() else float('nan'):.1f}%")
    print(f"  own-fill share of S (1-α_sup): median {np.median(own_fill):.3f}  "
          f"mean {own_fill.mean():.3f}  | on fires: median {np.median(own_fill[fires]) if fires.any() else float('nan'):.3f}"
          f"  | 1-comp: median {np.median(own_fill[one]) if one.any() else float('nan'):.3f}")
    # Q3: the two quadratic forms
    w_g = T["w_g"]
    aSa_code = T["aSa"]
    aSa_restr = T["s2cm"] * a_sup * a_sup + (
        np.where(T["sup_g"], T["alpha_g"] ** 2 * w_g, 0.0)
        + np.where(T["sup_r"], (T["alpha_p"] + T["alpha_n"]) ** 2 * T["w_r"], 0.0)
    )
    rr = aSa_code / np.maximum(aSa_restr, 1e-30)
    print(f"  Q3  αᵀΣα(code)/αᵀΣα(restricted Σ): median {np.median(rr):.3f}  p90 {np.percentile(rr, 90):.3f}"
          f"   (>1 ⇒ the code is the CONSERVATIVE one)")
    rows.append((name, fm[under].sum() / max(tot, 1e-9), fm[one].sum() / max(tot, 1e-9),
                 np.median(own_fill[fires]) if fires.any() else np.nan))

print(f"\n{'=' * 118}\nSUMMARY   share of the DAMPING MASS (Σ b̂²·M) by stratum\n{'=' * 118}")
print(f"  {'condition':<36}{'UNDER-claim':>14}{'one-component':>16}{'med own-fill|fires':>22}")
for r in rows:
    print(f"  {r[0]:<36}{100 * r[1]:>13.1f}%{100 * r[2]:>15.1f}%{r[3]:>22.3f}")
