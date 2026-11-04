"""The composition-transfer message policy: its name, its strand model and its library reduction.

A message is either a composition profile (a max-normalised log-likelihood over the solve grid of the destination's
gDNA share, carried across one face by a derived map) or a population level (an absolute profile over the log
density, carried where composition cannot cross, as a bound). Each node states its own claim from data only; each
directed face carries the recipient's rule (stop, forward, or a map); and the two passes and the solve do the rest. A
message is built from the source's claim and the recipient's constants and observations, never from the recipient's
belief; a flat profile is no claim, never a zero-filled one; and no node hears its own claim back. The builders
(``splice_faces``, ``edge_level``, ``terminus_rules``, ``alternative_splice_site``, ``gdna_lane``, ``rna_lane``), the
passes and the solve are native: ``native/transfer_kernel.h``, with the row constructors in ``native/transfer_rows.h``,
run per block by ``native/solve_kernel.cpp``.

What stays in Python is the library (:meth:`TransferPolicy.library`): the level lanes' two coordinates and whether the
strand split is a live witness, reduced once per sweep over the whole chain from observations and geometry. They are
the only things a message knows about slots outside its own block.

Gate: ``tests/calibration/test_transfer_policy.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ...native import transfer_rows as _rows
from . import ChainView

__all__ = ["TransferPolicy"]


@dataclass(frozen=True, slots=True)
class _Library:
    """What the whole library says, reduced once per sweep and read by every block.

    ``rho_gdna``: the structurally pure gDNA density, the gDNA lane's coordinate (``0.0``: no positive density
    anywhere, so no level lane can be built). ``rho_rna``: the RNA lanes' one coordinate, the unspliced RNA density
    over the single-strand exons with both strands pooled, or over every exon when none has RNA opportunity; when that
    is zero, the available RNA count/opportunity pairs supply the scale. A level is absolute and the coordinate only
    its origin, so one coordinate serves both strands. ``split_live``: whether the library's protocol verdict permits
    local strand-difference witnesses."""

    rho_gdna: float
    rho_rna: float
    split_live: bool


class TransferPolicy:
    """The composition-transfer policy. Observed strand claims use the fitted RNA sense fraction ``kappa``: eligible
    regions and boundaries supply the conditional-binomial likelihood of their observed strand counts."""

    name = "transfer"

    def __init__(self, kappa: float):
        self.kappa = float(kappa)

    # ── the library: the lanes' coordinates and the strand witness's liveness, once ─────────────────
    def library(self, view: ChainView) -> _Library:
        is_exon = np.asarray(view.is_exon_region, bool)
        fp, fn = np.asarray(view.free_pos, bool), np.asarray(view.free_neg, bool)
        n_u = np.asarray(view.n_slot, np.float64)
        a_g = np.asarray(view.eff_gdna, np.float64)
        a_r = np.asarray(view.eff_rna, np.float64)
        cnt = np.asarray(view.unspliced_count, np.float64)
        # the gDNA lane's coordinate: the structurally pure gDNA density over the intergenic regions
        # (``ChainView.terminal``). It is only the reference of a log axis, so any positive density serves; with no
        # intergenic reads the library-wide density stands in, because a zero coordinate would build no gDNA lane.
        pure = view.terminal & (a_g > 0.0)
        rho = float(n_u[pure].sum() / a_g[pure].sum()) if a_g[pure].sum() > 0.0 else 0.0
        if not rho > 0.0:
            den = float(a_g[a_g > 0.0].sum())
            rho = float(n_u[a_g > 0.0].sum()) / den if den > 0.0 else 0.0
        kappa = self.kappa
        # the split is a witness of a strand's RNA only where the library's strand channel is live
        # (``ChainView.strand_live``); each hop reads its own observations.
        single = fp != fn
        split_live = bool(view.strand_live)
        # the RNA lanes' coordinate: the unspliced RNA density over the single-strand exons, each strand read on its
        # own column (the kernel's ``read_column``: its own under a sense protocol, the other under an antisense one),
        # both strands pooled; every exon when none has RNA opportunity
        num = den = 0.0
        for free, col in ((fp, 0), (fn, 1)):
            sel = is_exon & free & single & (a_r > 0.0)
            col_read = int(_rows.read_column(col, kappa))
            num += float(cnt[sel, col_read].sum())
            den += float(a_r[sel].sum())
        if den <= 0.0:
            sel = is_exon & (a_r > 0.0)
            num, den = float(cnt[sel].sum()), float(a_r[sel].sum())
        rho_rna = num / den if den > 0.0 else 0.0
        if rho_rna <= 0.0:
            # a log coordinate must stay positive wherever local RNA evidence can exist; these pairs set only its
            # units, and the native builders still select the sources
            available = (fp | fn) & (n_u > 0.0) & (a_r > 0.0)
            num, den = float(n_u[available].sum()), float(a_r[available].sum())
            for flux, rate in (
                (view.sj_count_lo, view.route_rate_lo),
                (view.sj_count_hi, view.route_rate_hi),
            ):
                available = (flux > 0.0) & (rate > 0.0)
                num += float(flux[available].sum())
                den += float((flux[available] / rate[available]).sum())
            if num > 0.0:
                rho_rna = num / den
        return _Library(rho, rho_rna, split_live)
