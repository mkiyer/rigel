"""Is the κ-band gDNA under-call dataset-SIZE dependent? (od_g is a GLOBAL fit shrunk to a heavy prior.)

Build single-strand-only toys of increasing gene count at κ=0.7 (the worst band); report the fitted
od_g (global) and the unexpressed-gene (pure-gDNA) exon bias. If more seeds → od_g→0 → bias→0, the
toy anomaly is the prior dominating a small seed set (the global fit needs enough seeds), and on a huge
real dataset od_g is accurate. If the bias persists at large size, it's a deeper bug.
"""
from __future__ import annotations
import sys
import io
import contextlib
import numpy as np
from scripts.debug.toy_prod import run
from rigel.calibration.signature import BIT_EXON_POS, BIT_EXON_NEG
import rigel.calibration.gdna_strand as gs

SIZES_BP = [300, 1200, 5000]
INTRON = 3000
GAP = 5000


def genes_for(n_genes):
    genes = []
    x = GAP
    for k in range(n_genes):
        s = x
        ex = []
        for sz in SIZES_BP:
            ex.append((s, s + sz)); s = s + sz + INTRON
        genes.append((f"G{k}", "+", ex, 100.0 if k % 2 == 0 else 0.0))
        x = s + GAP
    return genes, x + GAP


_odbox = {}
_orig = gs._fit_overdispersion
def _wrap(*a, **k):
    r = _orig(*a, **k)
    _odbox.setdefault("hits", []).append((r[0], r[1]))  # (od, n_seed_nodes), first call = gDNA
    return r
gs._fit_overdispersion = _wrap


def one(n_genes, kappa=0.7):
    genes, glen = genes_for(n_genes)
    _odbox["hits"] = []
    with contextlib.redirect_stdout(io.StringIO()):
        rdf, solved, truth, tg, tr, sig = run(
            f"sz{n_genes}", genes, kappa=kappa, n_rna=12000 * n_genes, gdna_fraction=0.3,
            capture=True, capture_strength=20.0, genome_length=glen, seed=17)
    od_g, n_seed = _odbox["hits"][0]
    is_exon = (sig & (BIT_EXON_POS | BIT_EXON_NEG)) != 0
    m = is_exon & (tr == 0) & (tg > 0) & ((tg + tr) >= 8) & np.isfinite(solved) & np.isfinite(truth)
    bias = float(np.mean(solved[m] - truth[m])) if m.sum() else float("nan")
    print(f"  n_genes={n_genes:3d} | od_g_fit={od_g:.4f} (n_seed={n_seed:4d}) | "
          f"unexpressed-exon bias={bias:+.3f} (n={int(m.sum())})")


if __name__ == "__main__":
    print("κ=0.7 capture: dataset-size dependence of od_g over-fit + pure-gDNA-exon under-call")
    for n in (6, 12, 24, 48):
        one(n)
