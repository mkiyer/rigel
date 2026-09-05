"""The ROW CONSTRUCTORS of the composition transfer — pure functions of one face's numbers, each a
max-normalised log-likelihood over the solve grid ``lam`` (``f_g = sigma(lam)``) that `TransferPolicy`
delivers as ``PsiMessage.lam_rows``. Nothing here reads the chain or the context; everything here is
unit-gated in ``tests/calibration/test_transfer_policy.py``.

The laws every constructor keeps: every ratio is formed WITHIN one locale, so no level ever crosses a
capture cliff; a claim below its own evidence is silence (an all-zero row), never a near-zero row; the
width is derived from counts, never a constant.
"""

from __future__ import annotations

import numpy as np
from scipy.special import polygamma
from scipy.stats import norm

from ..splice_graph import (
    FLAG_ACCEPTOR_NEG as _ACC_NEG,
    FLAG_ACCEPTOR_POS as _ACC_POS,
    FLAG_DONOR_NEG as _DON_NEG,
    FLAG_DONOR_POS as _DON_POS,
    FLAG_TES_NEG as _TES_NEG,
    FLAG_TES_POS as _TES_POS,
    FLAG_TSS_NEG as _TSS_NEG,
    FLAG_TSS_POS as _TSS_POS,
)

__all__ = [
    "level_bound_row",
    "level_map_lambda",
    "level_row",
    "blur_row",
    "boundary_shares_strand",
    "edge_bound_row",
    "face_is_licensed",
    "face_map_lambda",
    "junction_flanks",
    "outside_flank",
    "splice_out_row",
    "transport_row",
]

EPS = 1.0e-9
TERMINUS = _TSS_POS | _TSS_NEG | _TES_POS | _TES_NEG
SJ_FLAGS = _DON_POS | _DON_NEG | _ACC_POS | _ACC_NEG
#: a transcript body that extends genomic-RIGHT from its terminus (a + start, a − end) leaves the
#: OUTSIDE flank on the left; one that extends LEFT (a + end, a − start) leaves it on the right
_BODY_RIGHT = _TSS_POS | _TES_NEG
_BODY_LEFT = _TES_POS | _TSS_NEG
#: the marginal over ``log rho`` is taken on equal-probability nodes of the standard normal —
#: quadrature resolution, like ``n_grid``, not a model constant
_MARGINAL_NODES = norm.ppf((np.arange(9) + 0.5) / 9.0)


def blur_row(row, lam, v):
    """The delta-method counting width: a Gaussian blur of variance ``v`` along ``lam`` applied to a
    max-normalised log-row (the kernel `transport_row`, `level_row` and item 7's pair width share)."""
    out = np.asarray(row, np.float64) - np.max(row)
    if v > 0.0 and lam.shape[0] > 1:
        dlam = float(lam[1] - lam[0])
        half = max(int(np.ceil(4.0 * float(np.sqrt(v)) / dlam)), 1)
        x = np.arange(-half, half + 1) * dlam
        kern = np.exp(-0.5 * x * x / v)
        kern /= kern.sum()
        pr = np.convolve(np.pad(np.exp(out), half, mode="edge"), kern, mode="valid")
        out = np.log(np.maximum(pr, 1.0e-300))
    return out - out.max()


def face_is_licensed(flags_b, fp_e, fn_e, fp_i, fn_i) -> bool:
    """A face is licensed for composition transfer when NO transcript terminus sits on it and both
    flanks admit the SAME strand set. An UNMEASURED population change (a terminus admits uncounted
    molecules; a strand flip changes membership) refuses the hop; the measured spliced route does
    not — it joins the message instead."""
    return (not (int(flags_b) & TERMINUS)) and bool(fp_e) == bool(fp_i) and bool(fn_e) == bool(fn_i)


def boundary_shares_strand(fp_b, fn_b, fp_i, fn_i) -> bool:
    """An intron|exon boundary may send its own strand row to its intron when the two admit the SAME
    strand set and that set is a SINGLE strand: the row is a statement about one live strand, and an
    AMBIG boundary's strand split constrains only the tilt, never the gDNA level (the Schur complement
    the local solve already applies). A terminus flag does not refuse: a true intron flank carries no
    exon bit, so no transcript terminating at the boundary covers it — the intron is always the
    OUTSIDE flank, whose composition the crossing shares."""
    return bool(fp_b) == bool(fp_i) and bool(fn_b) == bool(fn_i) and bool(fp_b) != bool(fn_b)


def outside_flank(flags_b, left, right):
    """The flank of a TERMINUS boundary that the terminating transcripts do NOT cover — the one whose
    population is exactly what crosses the boundary, spliced crossing included — and the covered
    (INSIDE) flank, as ``(outside, inside)``; ``(None, None)`` when the termini point both ways, when
    no terminus sits here, or when a splice junction shares the boundary (the sj+terminus case is its
    own item). The orientation is read off the flag alone: TSS+ and TES− bodies extend genomic-right,
    so the outside is the LEFT flank; TES+ and TSS− extend left, so it is the RIGHT."""
    f = int(flags_b)
    if not (f & TERMINUS) or (f & SJ_FLAGS):
        return None, None
    to_right, to_left = bool(f & _BODY_RIGHT), bool(f & _BODY_LEFT)
    if to_right and not to_left:
        return left, right
    if to_left and not to_right:
        return right, left
    return None, None


#: a DONOR bit marks the intron's LOW end on either strand (the flags are genomic-order), so the
#: intron lies to the boundary's right; an ACCEPTOR bit marks its HIGH end, the intron to the left
_INTRON_RIGHT = _DON_POS | _DON_NEG
_INTRON_LEFT = _ACC_POS | _ACC_NEG


def junction_flanks(flags_b, left, right):
    """At an exon|exon boundary carrying a splice junction and no terminus: ``(C, E)`` — the flank on
    the junction's INTRON side, which shares the boundary's full unspliced crossing (the outside-flank
    law, `splice_out_row` with the spliced crossing alone), and the flank where both isoforms are
    exonic, which holds the crossing plus the isoform that splices out here, MEASURED at the face as
    the route flux (the splice-in law with the spliced crossing plus the flux). ``(None, None)`` when a
    terminus shares the boundary (its own case), when no junction sits here, or when junctions leave
    both ways."""
    f = int(flags_b)
    if (f & TERMINUS) or not (f & SJ_FLAGS):
        return None, None
    to_right, to_left = bool(f & _INTRON_RIGHT), bool(f & _INTRON_LEFT)
    if to_right and not to_left:
        return right, left
    if to_left and not to_right:
        return left, right
    return None, None


def face_map_lambda(lam, n_u, a_g_b, a_r_b, e_g_e, e_r_e, s):
    """The SPLICE-IN face map ``lam_e(lam_u)``: the exon's gDNA log-odds as a function of the
    crossing's, given the face's unspliced count ``n_u``, the face (``a``) and exon (``e``)
    opportunities and the measured spliced density ``s`` that joins on the exon's side —

        lam_e = log(n_u*sigma(lam_u)*e_g/a_g) − log((n_u*(1−sigma(lam_u))/a_r + s)*e_r)

    monotone nondecreasing, saturating at the certified-flux ceiling when ``s > 0`` (a measured spliced
    density CAPS the claimable gDNA share) and a pure opportunity shift at ``s = 0``."""
    lam = np.asarray(lam, np.float64)
    sig = 1.0 / (1.0 + np.exp(-lam))
    g_arm = float(n_u) * sig / float(a_g_b) * float(e_g_e)
    r_arm = (float(n_u) * (1.0 - sig) / float(a_r_b) + float(s)) * float(e_r_e)
    return np.log(np.maximum(g_arm, 1e-300)) - np.log(np.maximum(r_arm, 1e-300))


def transport_row(row, lam, lam_e_of_u, n_u, n_s):
    """A source row carried INTO an exon through the face map: the row read at the map's preimage
    (a likelihood evaluated at the corresponding point — no Jacobian), widened by the face's counting
    variance ``trigamma(n_u+1/2) + trigamma(n_s+1/2)`` (the delta-method price of the map's measured
    position), max-normalised. Beyond the map's ceiling the row takes the flat limit — a soft cap."""
    lam = np.asarray(lam, np.float64)
    r = np.asarray(row, np.float64)
    lam_u_of_x = np.interp(lam, np.asarray(lam_e_of_u, np.float64), lam, left=lam[0], right=lam[-1])
    out = np.interp(lam_u_of_x, lam, r - r.max())
    v = float(polygamma(1, float(n_u) + 0.5) + polygamma(1, float(n_s) + 0.5))
    return blur_row(out, lam, v)


def splice_out_row(row_e, lam, n_u, n_s, a_g_b, a_g_e):
    """An exon's own row carried OUT to its intron|exon boundary: the splice-in map read backwards
    (``f_b = f_E (n_u + n_s) / n_u`` — the enrichment ratio cancels, only the face's own
    spliced-to-unspliced ratio survives), marginalised over that ratio's counting uncertainty
    ``log rho ~ N(log n_s/n_u, trigamma(n_s+1/2) + trigamma(n_u+1/2))`` on equal-probability nodes.

    Both components convert counts to densities on the SAME, capture-blind opportunity ``a_g`` (the
    spliced density is ``n_s / a_g_b``): a capture-aware opportunity on one component alone
    re-introduces a level across locales. Vacuous (all zeros) at a depleted face or a flat row."""
    lam = np.asarray(lam, np.float64)
    r = np.asarray(row_e, np.float64)
    if not (n_u > 0.0 and a_g_b > 0.0 and a_g_e > 0.0) or np.ptp(r) <= EPS:
        return np.zeros_like(lam)
    r = r - r.max()
    sd = float(np.sqrt(polygamma(1, float(n_s) + 0.5) + polygamma(1, float(n_u) + 0.5)))
    acc = np.zeros_like(lam)
    for z in _MARGINAL_NODES:
        s = float(n_s) / float(a_g_b) * float(np.exp(z * sd))
        m = face_map_lambda(lam, n_u, a_g_b, a_g_b, a_g_e, a_g_e, s)
        acc += np.exp(np.interp(m, lam, r, left=r[0], right=r[-1]))
    out = np.log(np.maximum(acc / _MARGINAL_NODES.size, 1.0e-300))
    out -= out.max()
    return out if np.ptp(out) > EPS else np.zeros_like(lam)


def level_map_lambda(lam, density_b, opportunity_i, total_i):
    """THE LEVEL-KEPT MAP ``lam_i(lam_b)``: the gDNA level crosses a face composition cannot — the
    inside's gDNA share is the sender's gDNA share times the sender's crossing density
    ``density_b`` (its total per base of gDNA opportunity), times the inside's gDNA opportunity, over
    the inside's OWN total (an observation). Monotone nondecreasing; clipped at the grid's ends."""
    lam = np.asarray(lam, np.float64)
    sig = 1.0 / (1.0 + np.exp(-lam))
    f_i = np.clip(sig * float(density_b) * float(opportunity_i) / float(total_i), 1e-9, 1.0 - 1e-9)
    return np.log(f_i / (1.0 - f_i))


def level_row(row_b, lam, lam_i_of_b, v):
    """A sender's profile carried across a level face: read at the map's preimage (its shape —
    and its sidedness — preserved; no value is summarised), then blurred by ``v``: the two totals'
    counting plus the pair's own dampening. Vacuous at a flat profile."""
    r = np.asarray(row_b, np.float64)
    if np.ptp(r) <= EPS:
        return np.zeros_like(np.asarray(lam, np.float64))
    lam = np.asarray(lam, np.float64)
    pre = np.interp(lam, np.asarray(lam_i_of_b, np.float64), lam, left=lam[0], right=lam[-1])
    out = blur_row(np.interp(pre, lam, r - r.max()), lam, v)
    return out if np.ptp(out) > EPS else np.zeros_like(lam)


def level_bound_row(lam, density_b, opportunity_i, total_i, v):
    """The level every library can send: the crossing's TOTAL density bounds the gDNA density
    above, so the inside's gDNA share times its own density may not exceed it — one-sided, at the
    level face's width ``v``."""
    lam = np.asarray(lam, np.float64)
    sig = 1.0 / (1.0 + np.exp(-lam))
    x = np.log(np.maximum(sig * float(total_i) / float(opportunity_i), 1e-300)) - np.log(
        float(density_b)
    )
    row = -0.5 * np.maximum(x, 0.0) ** 2 / max(float(v), 1e-12)
    return row - row.max()


def edge_bound_row(lam, n_b, n_e, a_g_b, a_g_e):
    """The intergenic|exon EDGE's claim on its exon, LOWER BOUND ONLY. The edge's crossing is
    structurally pure gDNA and an exon is at least as enriched as its own edge (probe panels target
    exons), so the honest claim is the profile likelihood over the nuisance enrichment ``s >= 1``,
    ``sup_s Pois(n_b; c(lam)/s)`` with ``c = sigma(lam)*n_e*a_g_b/a_g_e``: exactly 0 wherever
    ``c >= n_b`` (some enrichment explains any excess), the edge count's own one-sided Poisson tail
    below it, and identically zero at ``n_b = 0`` — a zero edge is vacuous, never a claim."""
    lam = np.asarray(lam, np.float64)
    n_b = float(n_b)
    if not n_b > 0.0:
        return np.zeros_like(lam)
    sig = 1.0 / (1.0 + np.exp(-lam))
    c = sig * float(n_e) * float(a_g_b) / float(a_g_e)
    return np.where(c >= n_b, 0.0, n_b * np.log(np.maximum(c, 1e-300) / n_b) - (c - n_b))
