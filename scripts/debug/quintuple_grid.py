"""The QUINTUPLE grid (owner's interrogation 2): intron R - IE B - exon R - EI B - intron R, a 3-exon gene's
MIDDLE exon flanked by two introns. The exon is a RELAY, not a source (it knows nothing at pass-0; its gDNA is
imputed from the introns - the 3-hop game of telephone). Sweep gDNA x nascent x mature and ask: does the geo-mode
message propagation resolve the middle exon across the WHOLE grid, especially the nascent>>mature corner where
change-1 (suppress exon unspliced RNA) is most suspect?

Reports the middle-exon region solved vs oracle f_g (toy_prod's production-faithful scan+calibrate). Run with
RIGEL_BND_GEO_GDNA unset (baseline) and =1 (geo-mode) to A/B."""
from __future__ import annotations
import contextlib
import io
import os
import sys
import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
import toy_prod
from rigel.calibration.signature import BIT_EXON_POS, BIT_EXON_NEG

# 3-exon gene: exons 1000-2000, 4000-5000, 7000-8000; introns 2000-4000 and 5000-7000 (genome 10000).
# middle exon (4000-5000) = quintuple centre: intron(2000-4000) - IE - exon - EI - intron(5000-7000).
EXONS = [(1000, 2000), (4000, 5000), (7000, 8000)]

# (label, gdna_fraction, nascent, mature). mature >= 20 always (calibrate rejects a zero-spliced library).
GRID = [
    ("mature>>nasc (panel)", 0.5, 20, 400),
    ("nascent>>mat (EDGE)", 0.5, 400, 20),
    ("balanced", 0.5, 200, 200),
    ("gDNA>> nascent>>mat", 0.9, 400, 20),
    ("gDNA>> mature>>nasc", 0.9, 20, 400),
    ("gDNA~0 nascent>>mat", 0.05, 400, 20),
    ("gDNA~0 mature>>nasc", 0.05, 20, 400),
]


def _mid_exon_fg(name, gf, nasc, mat):
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            rdf, solved, truth, tg, tr, sig = toy_prod.run(
                name, [("G", "+", EXONS, float(mat))], kappa=0.5, n_rna=4000, genome_length=10000,
                gdna_fraction=gf, nascent=float(nasc), capture=True, capture_strength=20.0, seed=7)
    except Exception as e:
        return np.nan, np.nan, type(e).__name__
    sig = np.asarray(sig).astype(int)
    exon = [i for i in range(len(sig)) if sig[i] & (BIT_EXON_POS | BIT_EXON_NEG)]
    mid = exon[len(exon) // 2] if exon else -1
    return (solved[mid], truth[mid], "") if mid >= 0 else (np.nan, np.nan, "no-exon")


mode = "GEO" if os.environ.get("RIGEL_BND_GEO_GDNA") == "1" else "baseline"
print(f"[mode={mode}]  middle-exon f_g across the quintuple grid (capture ON, ss=0.5)\n")
print(f"{'regime':>22} | {'gDNA':>5} {'nasc':>5} {'mat':>5} | {'oracle':>6} {'solved':>6} {'err':>6} {'note':>10}")
for lab, gf, nasc, mat in GRID:
    s, t, note = _mid_exon_fg(f"q_{gf}_{nasc}_{mat}_{mode}", gf, nasc, mat)
    err = "" if np.isnan(t) or np.isnan(s) else f"{s - t:+.2f}"
    ts = "nan" if np.isnan(t) else f"{t:.2f}"
    ss = "nan" if np.isnan(s) else f"{s:.2f}"
    print(f"{lab:>22} | {gf:>5.2f} {nasc:>5} {mat:>5} | {ts:>6} {ss:>6} {err:>6} {note:>10}")
