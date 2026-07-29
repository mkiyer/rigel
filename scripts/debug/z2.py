"""THE z2 calibration statistic — one definition, one denominator.

``z2 = E[(f_g − oracle)²] / E[Var(f_g)]``, mass-weighted. **1.0 = honest, > 1 = over-confident.**

⚠ **Why this module exists.** ``NodeBelief.var_gdna`` is a grid moment of ``log f_g``
(:func:`rigel.calibration.simplex_logodds._solve_nodes_logodds`) — a **log-space** variance, unbounded above
(measured: it exceeds ¼ on 33 % of scored nodes and reaches 4.48). The numerator is a **linear** squared
fraction error. Dividing one by the other compares two different units, and **every z2 recorded before
2026-07-28 is on that mixed scale** — the suite total read 0.046 ("20× conservative") when it is really ≈1,
and the per-class ranking inverts once corrected.

That was fixed in ``pass0_error_table.py`` and **nowhere else**: ``ablate_report.py``,
``v_p4_paired_z2.py``, ``ho_subz2.py`` and ``v_ho_subcheck.py`` all kept the raw log-space denominator, so
the standing gate "held-fixed z2 must not regress" was being scored on two incompatible scales depending on
which tool ran. This module is the single source of truth; every consumer imports it.

    import sys; sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
    from z2 import lin_var, z2
"""

from __future__ import annotations

import numpy as np

__all__ = ["lin_var", "z2"]

_EPS = 1e-9


def lin_var(var_log, f_g):
    """Convert the belief's log-space ``Var(log f_g)`` to the linear ``Var(f_g)`` the error is measured in.

    Two derived steps, no tuned constant:

    * ``f_g = exp(log f_g)`` ⇒ ``Var(f_g) = f_g²·(e^v − 1)`` **exactly** under a lognormal, reducing to the
      first-order delta method ``f_g²·v`` as ``v → 0``. The exact form matters because ``v`` reaches ~4.5,
      where the first-order term is not a usable approximation.
    * capped at ``f_g·(1 − f_g)`` — the greatest variance ANY ``[0, 1]`` variable with mean ``f_g`` can have.
      This is the same bound ``bp_solver.node_sweep`` applies to ``_var_fg`` before ``composition_logvar``.

    The 700 in the clip is a float64 overflow guard on ``expm1``, not a modelling constant: the cap binds for
    any ``v`` beyond a few units, so clipping the exponent cannot change the result.
    """
    v = np.asarray(var_log, dtype=np.float64)
    f = np.asarray(f_g, dtype=np.float64)
    return np.minimum(f * f * np.expm1(np.clip(v, 0.0, 700.0)), f * (1.0 - f))


def z2(mass, err_lin, var_log, f_g, mask=None):
    """Mass-weighted ``E[err²] / E[Var(f_g)]`` over ``mask``, with both sides in LINEAR fraction units.

    ``err_lin`` is the per-object ``|f_g − oracle|`` (linear). Objects whose ``var_log`` is non-finite carry
    no opinion and are excluded — a node that declares infinite variance cannot be over- or under-confident.

    The numerator is deliberately left in linear units: a log-space numerator would need ``log(oracle)``, and
    the oracle ``f_g`` is exactly 0 on 23.5 % of scored nodes.
    """
    m = np.ones(np.shape(mass), dtype=bool) if mask is None else np.asarray(mask, dtype=bool).copy()
    v = np.asarray(var_log, dtype=np.float64)
    m &= np.isfinite(v)
    if not m.any():
        return float("nan")
    w = np.asarray(mass, dtype=np.float64)[m]
    num = float(np.sum(w * np.asarray(err_lin, dtype=np.float64)[m] ** 2))
    den = float(np.sum(w * lin_var(v[m], np.asarray(f_g, dtype=np.float64)[m])))
    return num / den if den > 0 else float("nan")
