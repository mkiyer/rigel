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
    "edge_level_row",
    "face_is_licensed",
    "face_map_lambda",
    "junction_exon_side",
    "junction_flanks",
    "outside_flank",
    "splice_out_row",
    "transport_row",
    "hop_price",
    "intersect",
    "level_of_profile",
    "lower_side",
    "poisson_level",
    "profile_of_level",
    "cube_row",
    "flux_level",
    "read_column",
    "rna_level_of_profile",
    "rna_row_of_level",
    "strand_bits",
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
    max-normalised log-row (the kernel `transport_row`, `level_row` and the alternative splice site's
    pair width share)."""
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
    (INSIDE) flank, as ``(outside, inside)``; ``(None, None)`` when the termini point both ways or when
    no terminus sits here. The orientation is read off the terminus flag alone: TSS+ and TES− bodies
    extend genomic-right, so the outside is the LEFT flank; TES+ and TSS− extend left, so it is the
    RIGHT. A splice junction sharing the boundary does not change which flank is inside — the
    sj+terminus case, where a transcript starts or ends exactly at another isoform's exon edge; that
    junction's flux is placed by `junction_exon_side` instead."""
    f = int(flags_b)
    if not (f & TERMINUS):
        return None, None
    to_right, to_left = bool(f & _BODY_RIGHT), bool(f & _BODY_LEFT)
    if to_right and not to_left:
        return left, right
    if to_left and not to_right:
        return right, left
    return None, None


def junction_exon_side(flags_b, left, right):
    """The flank on a junction's EXON side: a DONOR bit marks the intron's low end (the intron lies right,
    the exon left), an ACCEPTOR bit its high end; ``None`` when junctions leave both ways or none is
    present. At an sj+terminus boundary the junction's measured flux belongs to the population of this
    flank — the RNA that splices in or out here — and is placed there: in the terminus's outside map when
    this flank is the outside, in the terminus rule's totals' disagreement when it is the inside."""
    f = int(flags_b)
    don, acc = bool(f & (_DON_POS | _DON_NEG)), bool(f & (_ACC_POS | _ACC_NEG))
    if don and not acc:
        return left
    if acc and not don:
        return right
    return None


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


def edge_level_row(lam, n_b, n_e, a_g_b, a_g_e):
    """THE EDGE'S LEVEL, ONE-SIDED. The intergenic|exon edge's crossing is structurally pure gDNA, so
    its COUNT measures the gDNA level the exon continues; for each hypothesised exon share the implied
    edge count is ``c = sigma(lam) * n_e * a_g_b / a_g_e``. BELOW the edge's level the row is the
    count's exact Poisson likelihood ``n_b log(c/n_b) − (c − n_b)`` — the exon has at least the edge's
    gDNA density, at counting width — and NOTHING above it.

    ⛔ The upper side has no honest form and must not be added: capture enriches a probed interior over
    its edge by an amount no local witness measures, so a two-sided or dampened upper side pulls an
    unstranded exon toward a centre below the truth. A ZERO count is vacuous for the same reason —
    under capture a dark edge beside a probed exon is not an empty one — and the zero controls such a
    claim would win belong to the landscape prior."""
    lam = np.asarray(lam, np.float64)
    n_b = float(n_b)
    if not n_b > 0.0:
        return np.zeros_like(lam)
    sig = 1.0 / (1.0 + np.exp(-lam))
    c = sig * float(n_e) * float(a_g_b) / float(a_g_e)
    return np.where(c >= n_b, 0.0, n_b * np.log(np.maximum(c, 1e-300) / n_b) - (c - n_b))


# ── THE LEVEL LANE: a gDNA level as an absolute profile over u = log(rho / rho_ref) ──────────────


def poisson_level(u, n, a, rho_ref):
    """A structurally pure gDNA count ``n`` on gDNA opportunity ``a`` — the gene edge's crossing — as a
    level profile: the Poisson log-likelihood of the count at each density, max-normalised. A zero
    count is a profile falling with the density (nothing below zero is claimed)."""
    c = float(rho_ref) * np.exp(np.asarray(u, np.float64)) * float(a)
    n = float(n)
    row = n * np.log(np.maximum(c, 1e-300)) - c if n > 0.0 else -c
    return row - row.max()


def level_of_profile(row, lam, u, n, a, rho_ref):
    """A node's own composition profile (over ``lam``) read as a LEVEL profile (over ``u``) through the
    node's own total: the density ``rho`` implies the gDNA share ``rho a / n``, so the level profile is
    the composition profile read at that share; above the total the share is impossible and the level
    falls as the total's Poisson tail (the total bounds the level)."""
    row = np.asarray(row, np.float64)
    lam = np.asarray(lam, np.float64)
    c = float(rho_ref) * np.exp(np.asarray(u, np.float64)) * float(a)
    n = float(n)
    f = np.clip(c / n, EPS, 1.0 - EPS)
    out = np.interp(np.log(f / (1.0 - f)), lam, row - row.max())
    out = out + np.where(c >= n, n * np.log(np.maximum(c, 1e-300) / n) - (c - n), 0.0)
    return out - out.max()


def lower_side(profile):
    """A profile's LOWER SIDE: non-decreasing in ``u``, max-normalised — the claim "at least this much",
    which is all a level that crosses a face may say (the interior may be enriched over the source,
    never depleted; every upper side measured harmful on the stranded capture-ON rows)."""
    p = np.maximum.accumulate(np.asarray(profile, np.float64))
    return p - p.max()


def intersect(bounds):
    """Two or more bounds on ONE density combine by INTERSECTION: the pointwise minimum of their
    log-profiles (the tighter bound wins at each density), max-normalised. Bounds intersect, they do
    not multiply — a product of one-sided claims sharpens where nothing was measured."""
    out = None
    for b in bounds:
        b = np.asarray(b, np.float64)
        out = b if out is None else np.minimum(out, b)
    return None if out is None else out - out.max()


def profile_of_level(profile, u, lam, n, a, rho_ref):
    """A held level read as THIS node's composition profile through its own total — `level_of_profile`'s
    map read backwards, a pure coordinate change: ``u(lam) = log(sigma(lam) n / (a rho_ref))``."""
    p = np.asarray(profile, np.float64)
    lam = np.asarray(lam, np.float64)
    u_of_lam = np.log(1.0 / (1.0 + np.exp(-lam)) * float(n) / (float(a) * float(rho_ref)))
    out = np.interp(u_of_lam, np.asarray(u, np.float64), p, left=p[0], right=p[-1])
    return out - out.max()


def count_logvar(count) -> np.ndarray:
    """``Var(log rho)`` for a Poisson rate seen as ``count`` events over an opportunity — exactly, at
    every count including zero: under the Jeffreys prior the rate's posterior is ``Gamma(count + 1/2,
    E)``, whose log has variance ``trigamma(count + 1/2)``, independent of the opportunity ``E`` (it
    moves the location and cannot sharpen the claim). It agrees with the ``1/n`` it replaces at
    moderate counts; the whole difference is at small ones, and at ``n = 0`` it is ``pi^2/2`` rather
    than infinity — a zero count is a measurement, not an absence. THE ONE HOME of the counting term:
    every hop price here and `region_init`'s own precision read it, so there is one definition and
    nothing to keep in step."""
    return polygamma(1, np.asarray(count, np.float64) + 0.5)


def hop_price(n_s, a_s, n_x, a_x):
    """One hop's price on the lane's WITNESS counts, per hop and nothing pooled: both
    counts' counting (`count_logvar` each) plus the discrepancy of the two count densities beyond
    what counting explains, ``max(0, log(r)^2 - (1/n_s + 1/n_x))`` with ``r`` the ratio of the
    densities ``n_x / a_x`` and ``n_s / a_s``. A discrepancy is never attributed (capture, new
    transcription, noise): it widens. The gDNA lane prices on the two nodes' totals; an RNA lane on
    the strand's own counts, because the totals are the wrong witness for it (entering an overlap the
    total jumps because the OTHER strand joins while this strand's density is unchanged). A zero count
    on either side has no density ratio, so counting is the whole price (``trigamma(1/2)`` is
    finite)."""
    n_s, n_x = float(n_s), float(n_x)
    v = float(count_logvar(n_s) + count_logvar(n_x))
    if n_s > 0.0 and n_x > 0.0:
        r = (n_x / float(a_x)) / (n_s / float(a_s))
        v += max(0.0, float(np.log(r)) ** 2 - (1.0 / n_s + 1.0 / n_x))
    return v


# ── THE RNA LEVEL LANES (the both-stranded locus) ─────────────────────────────────────────────────

#: a strand's four boundary bits, its junction bits and its terminus bits, by strand key
strand_bits = {
    "pos": (
        _TSS_POS | _TES_POS | _DON_POS | _ACC_POS,
        _DON_POS | _ACC_POS,
        _TSS_POS | _TES_POS,
    ),
    "neg": (
        _TSS_NEG | _TES_NEG | _DON_NEG | _ACC_NEG,
        _DON_NEG | _ACC_NEG,
        _TSS_NEG | _TES_NEG,
    ),
}


def rna_level_of_profile(row, lam, u, n, a_r, rho_ref):
    """A single-strand node's own composition profile (over ``lam``) read as its live strand's RNA
    LEVEL (over ``u = log(rho_s / rho_ref_s)``) through the node's own total: the RNA density
    ``rho_s`` implies the RNA share ``rho_s a_r / n``, hence ``f_g = 1 − rho_s a_r / n`` and the ``lam``
    the profile is read at; above the total the share is impossible and the level falls as the
    total's Poisson tail. `level_of_profile` with ``1 − sigma`` in place of ``sigma``. Near
    ``f_r → 1`` the coordinate saturates at the total and several ``lam`` cells share one ``u`` cell —
    the same grid limit the gDNA level has at ``f_g → 1``."""
    row = np.asarray(row, np.float64)
    lam = np.asarray(lam, np.float64)
    c = float(rho_ref) * np.exp(np.asarray(u, np.float64)) * float(a_r)
    n = float(n)
    f_r = np.clip(c / n, EPS, 1.0 - EPS)
    out = np.interp(np.log((1.0 - f_r) / f_r), lam, row - row.max())
    out = out + np.where(c >= n, n * np.log(np.maximum(c, 1e-300) / n) - (c - n), 0.0)
    return out - out.max()


def rna_row_of_level(profile, u, lam, n, a_r, rho_ref):
    """A held RNA level read as THIS single-strand node's composition row — `rna_level_of_profile`'s map
    read backwards, a pure coordinate change: ``u_s(lam) = log((1 − sigma(lam)) n / (a_r rho_ref))``. A
    lower-only level (non-decreasing in u) is NON-INCREASING in lam: "at least this much RNA" is "at most
    this much gDNA", the upper side of the gDNA share, the partner a gDNA floor needs."""
    p = np.asarray(profile, np.float64)
    lam = np.asarray(lam, np.float64)
    f_r = 1.0 / (1.0 + np.exp(lam))
    u_of_lam = np.log(f_r * float(n) / (float(a_r) * float(rho_ref)))
    out = np.interp(u_of_lam, np.asarray(u, np.float64), p, left=p[0], right=p[-1])
    return out - out.max()


def flux_level(u, count, rate, rho_ref, v=0.0):
    """THE CERTIFIED FLUX at one of an exon's junctions as that strand's RNA level at the exon.
    Spliced fragments are RNA of a KNOWN strand: measured, never solved, and strictly one hop
    (boundary → exon). The junction's rate is an ESTIMATE of the exon's abundance, priced by the node
    pair's disagreement — the spliced count's Poisson profile at each hypothesised density on the route
    rate's own opportunity ``count / rate``, widened by the hop's price ``v``, then its LOWER SIDE.
    ⛔ Lower-only, because a two-sided estimate over-claims at the probe cliff: the exon takes "at least
    the RNA its junction's isoforms carry" and nothing above. A zero count claims nothing (``None``)."""
    count, rate = float(count), float(rate)
    if not (count > 0.0 and rate > 0.0):
        return None
    pl = poisson_level(u, count, count / rate, rho_ref)
    return lower_side(blur_row(pl, u, v) if v > 0.0 else pl)


def read_column(col, kappa):
    """The genome-strand column strand ``col``'s RNA READS on: its own when the library reads sense
    (``kappa >= 1/2``, or no fitted strand model), the other under an antisense protocol. A junction's
    route rate is in transcript-strand terms, so the exon count it is priced against must be the
    count of the reads that strand's RNA produces; reading the other column inverts a node's strand
    share and blurs its floor to nothing."""
    return int(col) if (kappa is None or float(kappa) >= 0.5) else 1 - int(col)


def cube_row(profiles, u, lam, theta, n, a_r, rho_refs):
    """The held RNA levels of an AMBIG node as ONE row over ψ's ``(lam, theta)`` cube: at each cell the
    strand's share ``f_s = (1 − sigma)(1 ± tau)/2`` implies the density ``f_s n / a_r``, and the held
    profile is read at ``log(rho_s / rho_ref_s)`` — the map `profile_of_level` applies on the λ axis,
    with the tilt inside. A one-sided profile stays one-sided (the map is monotone in each share), so
    "at least this much RNA+" arrives as a wall in the cube and no parametric summary is made.
    ``profiles`` is ``{"pos": profile, "neg": profile}`` (either may be absent), ``rho_refs`` the two
    lanes' coordinates."""
    lam = np.asarray(lam, np.float64)
    theta = np.asarray(theta, np.float64)
    sig = 1.0 / (1.0 + np.exp(-lam))
    tau = np.sin(theta)
    f_act = (1.0 - sig)[:, None]
    shares = {
        "pos": f_act * (1.0 + tau)[None, :] / 2.0,
        "neg": f_act * (1.0 - tau)[None, :] / 2.0,
    }
    out = np.zeros((lam.shape[0], theta.shape[0]))
    for name, prof in profiles.items():
        if prof is None:
            continue
        prof = np.asarray(prof, np.float64)
        with np.errstate(divide="ignore"):
            u_s = np.log(shares[name] * float(n) / float(a_r)) - np.log(float(rho_refs[name]))
        out += np.interp(u_s, np.asarray(u, np.float64), prof, left=prof[0], right=prof[-1])
    return out - out.max()
