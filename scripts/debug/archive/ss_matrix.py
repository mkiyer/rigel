"""Comprehensive single-strand-node optimization matrix (NO AMBIG nodes).

A single-strand-only scenario (all + genes, non-overlapping → zero AMBIG), mixing EXPRESSED and
UNEXPRESSED multi-exon genes, run across κ × capture × ±gDNA. Reports per-node-type accuracy
(expressed exon / unexpressed exon=pure gDNA / intron=gDNA) so we can optimize EVERY regime and catch
regressions (esp. unstranded). Toggle solver config via env: RIGEL_BSF (boundary spliced floor),
RIGEL_PRIOR_FREE. Usage:  RIGEL_BSF=1 python -m scripts.debug.ss_matrix
"""
from __future__ import annotations
import io
import contextlib
import os
import numpy as np

from scripts.debug.toy_prod import run
from rigel.calibration.signature import BIT_EXON_POS, BIT_EXON_NEG, BIT_INTRON_POS, BIT_INTRON_NEG

# 6 single-strand (+) genes, non-overlapping, 3 exons each (→ 2 introns); 3 expressed, 3 unexpressed.
def _gene(g0, ab):
    return [(g0[0], "+", g0[1], ab)]
GENES = []
base = 2000
for k in range(6):
    s = base + k * 12000
    exons = [(s, s + 1000), (s + 3000, s + 5000), (s + 8000, s + 9000)]
    ab = 100.0 if k % 2 == 0 else 0.0  # alternate expressed / unexpressed
    GENES.append((f"G{k}", "+", exons, ab))
GENOME = 2000 + 6 * 12000 + 3000


def _summary(solved, truth, tg, tr, sig):
    is_exon = (sig & (BIT_EXON_POS | BIT_EXON_NEG)) != 0
    is_intron = (sig & (BIT_INTRON_POS | BIT_INTRON_NEG)) != 0
    obs = (tg + tr) >= 5
    out = {}
    for name, m in [("exon_expr", is_exon & (tr > 0) & (tg > 0)),
                    ("exon_unexpr", is_exon & (tr == 0) & (tg > 0)),
                    ("exon_noGDNA", is_exon & (tr > 0) & (tg == 0)),  # FP test: oracle f_g=0, want ~0
                    ("intron", is_intron & (tg > 0))]:
        mm = m & obs & np.isfinite(solved) & np.isfinite(truth)
        if mm.sum() == 0:
            out[name] = None
        else:
            err = solved[mm] - truth[mm]
            out[name] = (int(mm.sum()), float(np.median(np.abs(err))), float(np.mean(err)))
    return out


def matrix():
    print(f"config: RIGEL_BSF={os.environ.get('RIGEL_BSF','1')} "
          f"RIGEL_PRIOR_FREE={os.environ.get('RIGEL_PRIOR_FREE','1')}")
    print(f"{'cap':>4} {'kappa':>6} {'gdna':>5} | "
          f"{'exon_expr':>22} {'exon_unexpr':>22} {'exon_noGDNA':>22} {'intron':>22}")
    for cap in (False, True):
        for kappa in (0.5, 0.7, 0.99):
            for gd, gdlab in ((0.5, "yes"), (0.0, "no")):
                buf = io.StringIO()
                with contextlib.redirect_stdout(buf):
                    rdf, solved, truth, tg, tr, sig = run(
                        f"ssm_{int(cap)}_{kappa}_{gdlab}", GENES, kappa=kappa, n_rna=60000,
                        gdna_fraction=gd, capture=cap, capture_strength=20.0,
                        genome_length=GENOME, seed=11)
                s = _summary(solved, truth, tg, tr, sig)
                def fmt(x):
                    return "  n=0 " if x is None else f"n={x[0]:2d} |e|={x[1]:.3f} b={x[2]:+.3f}"
                print(f"{('ON' if cap else 'off'):>4} {kappa:>6} {gdlab:>5} | "
                      f"{fmt(s['exon_expr']):>22} {fmt(s['exon_unexpr']):>22} "
                      f"{fmt(s['exon_noGDNA']):>22} {fmt(s['intron']):>22}")


if __name__ == "__main__":
    matrix()
