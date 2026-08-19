#!/usr/bin/env python
"""⭐⭐⭐ IS THE CERTIFIED COUNT COMPOSITION EVIDENCE? — the splice-visibility `q`, measured against TRUTH.

⛔ **THIS INSTRUMENT'S VERDICT IS A NEGATIVE, AND IT IS THE POINT OF THE FILE.** A spliced fragment cannot
be gDNA, so ``boundary_spliced`` is certified RNA — and the obvious next step is to make that certified count
speak about the slot's own gDNA-vs-RNA split. It cannot, and this measures why.

**The derivation being tested.** At one boundary the unspliced bank holds ``M`` crossings split ``f_g``
gDNA / ``1−f_g`` RNA, and the certified bank holds ``S`` crossings that spliced elsewhere. The two are ONE
population of contiguous crossings split by whether a splice was VISIBLE, so with ``q`` that visibility
probability::

    E[S]  =  (q / (1−q)) · (1−f_g) · M

    log P(S | f_g)  =  S·log(1 − f_g)   −   (q/(1−q))·(1−f_g)·M   +   const
                       \\___ retained ___/    \\_______ dropped _______/

⭐ The retained term's coefficient is the RAW COUNT and is **frame-free** — any opportunity ratio between
the two banks' divisors multiplies ``E[S]``, so in log space it is an additive constant and cannot reach
the coefficient of ``log(1−f_g)``. TRAPS: two-divisors-opposite-sign's opposite-sign trap is structurally absent here: neither
``eff_rna`` nor ``eff_sj`` appears in the retained term at all. **That half of the derivation is
sound and `tests/calibration/test_certified_rna_licence.py` gates it.**

⛔⛔ **WHAT IS NOT SOUND IS DROPPING THE SECOND TERM.** Keeping only the first makes the claim one-sided —
"there is at least this much RNA here" — which is honest ONLY if the dropped term is small, i.e. only if
``q → 1``. This file measures ``q`` directly from the origin-split oracle, where it is an identity::

    q  =  S / (S + C_R)          C_R = the TRUE RNA in the unspliced bank (mrna + nrna)

---

⭐⭐⭐ **WHAT IT MEASURED — RE-MEASURED 2026-08-17 ON THE 16-CONDITION LADDER. No solver runs, ~4 s.**

⛔ **The numbers below replace a table measured on the 36-condition ladder retired 2026-08-13, and the
replacement was a RE-RUN rather than an edit** (`TRAPS: re-record-the-baseline`). Every row that still
has a condition to stand on reproduced to the printed digit; the one row that did not is named under
point 2, because a rung that no longer exists cannot be re-measured and must not be quietly renumbered.
The command is the Usage line: ``python scripts/design/certified_q_census.py``.

**1. ``q`` is nowhere near 1.** Mass-weighted median ``q`` = **0.192–0.706**, and **60.4–97.9 % of the
mass sits at q < 0.9** in every one of the 16 conditions. So the dropped term is not a correction, it is
comparable to the term being kept.

**2. The raw-count term is therefore not a floor — it is a PRIOR TOWARD RNA, and it is panel-negative on
more than a third of the ladder.** Scoring the median of Beta(½, ½+S) — which is exactly ψ's reference plus
the term with nothing else speaking — against per-BOUNDARY truth, mass-weighted over the BOUNDARIES the
term fires on:

=========================  ==========  ==========  ==========  ==========
condition                  q med       mwae term   mwae ref    Δ
=========================  ==========  ==========  ==========  ==========
g00 ss0.50 capture_off     0.200       0.0107      0.5000      **−0.4893**
g00 ss0.50 capture_on      0.195       0.0077      0.5000      **−0.4923**
g50 ss0.50 capture_off     0.234       0.0695      0.4540      −0.3845
g50 ss0.50 capture_on      0.356       0.5245      0.3285      **+0.1960**
g98 ss0.50 capture_off     0.706       0.6854      0.3347      +0.3507
g98 ss0.50 capture_on      0.694       0.9230      0.4779      +0.4451
g98 ss0.99 capture_on      0.700       0.9244      0.4785      **+0.4459**
=========================  ==========  ==========  ==========  ==========

**WORSE on 6 of 16 conditions**; worst **+0.4459** at ``g98 ss0.99 capture_on``, where **100.0 % of the
mass** sits on BOUNDARIES whose truth is gDNA-rich (mass-weighted mean true ``f_g`` = **0.9785**) and the
term answers **0.0541**.

⛔ **ONE RECORDED NUMBER IS UNREPRODUCIBLE AND IS NOT REPLACED BY A NEARBY ONE.** The retired verdict's
worst row was ``+0.4578`` at ``g90 ss0.50 capture_on``. The rebuilt ladder's gDNA rungs are
``g00/g05/g50/g98`` — **there is no g90**, and ``g98`` is a different mixture, not a re-measurement of
the same one. What would reproduce it is a panel carrying a ``g90`` rung, and nothing on disk does.

⭐⭐ **3. And the two ends of the ladder ARE the two zero controls, which is what makes the diagnosis
certain rather than suggestive.** At ``g00`` the truth is all-RNA and the term is worth **−0.4893 …
−0.4923** — the best results on the panel. At ``g98 capture_on`` the truth is nearly all-gDNA and the term
is worth **+0.4451 … +0.4459** — the worst. The sign of the effect is the sign of the truth. ⛔ **A channel
whose benefit tracks the answer rather than the evidence has not added information; it has added a prior**
(TRAPS: honesty-metrics-reward-ignorance's shape).

**4. And ``q`` is a property of the RNA GEOMETRY, not of the composition.** Per-BOUNDARY ``q`` at ``g00`` vs
``g50`` agrees to a Spearman **ρ = 0.9257** on 5,241 shared BOUNDARIES — a 50-point swing in the gDNA fraction
barely moves it, because ``q`` asks only "of the RNA crossing this boundary, what fraction shows a splice?".
⚠ That pair is the one the code below reaches for by name, so it is re-derived on every run rather than quoted.
⭐ So the missing quantity is
computable in principle: ``1−q`` is the share of the crossing opportunity that fits inside the unbroken
EXONIC reach either side of the boundary. ⛔ But that array does not exist —
``splice_graph.build_contiguous_boundary_reach_arrays`` is deliberately GENOMIC (a contiguous boundary is also
crossed by nascent RNA, which does not splice), and the exonic version is per-TRANSCRIPT and
abundance-weighted, which calibration does not know at pass-0. That is the build this channel is blocked
on, and it needs its own brute-force enumeration gate.

Usage::

    python scripts/design/certified_q_census.py                 # the whole ladder
    python scripts/design/certified_q_census.py --conditions gdna_g50_ss_0.50_nrna_mid_capture_on
"""

from __future__ import annotations

import argparse
import dataclasses
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402
from scipy.stats import beta as _Beta  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests" / "calibration"))

from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import _native_detect_sj_tag  # noqa: E402
from rigel.scan_cache import read_scan_cache  # noqa: E402

DEFAULT_SUITE = Path.home() / "Downloads/rigel_runs/suite/ladder"
DEFAULT_INDEX = Path.home() / "Downloads/rigel_runs/suite/rigel_index"
ORIGINS = ("gdna", "mrna", "nrna")


def _load(cache_root: Path, suite: Path, tag: str, index, cfg):
    """The three origin partitions from the oracle cache. ⚠ Loaded through the SHIPPED
    ``read_scan_cache``, which refuses a payload whose ``graph_hash`` / ``reach_digest`` /
    ``payload_schema_digest`` does not describe this index — so a stale cache is refused loudly rather
    than feeding a truth number (`memory: reach_is_covered_by_no_hash`)."""
    bam = str(suite / tag / "sim_oracle.bam")
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    return {k: read_scan_cache(cache_root / tag / k, index, scan).payload for k in ORIGINS}


def measure(parts) -> dict:
    """Per-BOUNDARY truth on the CERTIFIED axis. Everything here is an identity on the oracle, not a fit.

    ⚠ ``true_fg`` is the gDNA fraction of the **unspliced** bank — the solver's own definition. Scoring
    against a spliced-INCLUSIVE fraction reads a small mass defect as a half-unit error, which is how the
    first run of `certified_rna_audit.py` mis-reported this channel.
    """
    es = sum(np.asarray(parts[k].boundary_spliced_count, np.float64).sum(1) for k in ORIGINS)
    eu = {k: np.asarray(parts[k].boundary_unspliced_count, np.float64).sum(1) for k in ORIGINS}
    c_g, c_r = eu["gdna"], eu["mrna"] + eu["nrna"]
    mass = c_g + c_r
    live = (es > 0) & (mass > 0)
    S, Cg, Cr, M = es[live], c_g[live], c_r[live], mass[live]
    return dict(
        live=live,
        S=S,
        C_g=Cg,
        C_R=Cr,
        M=M,
        # the splice VISIBILITY: of the RNA that crossed this boundary, the share that showed a splice.
        q=S / np.maximum(S + Cr, 1e-12),
        true_fg=Cg / M,
        # ⭐ what the RAW-COUNT term alone would answer: psi's reference (Beta(½,½)) updated by S
        # certified RNA observations is EXACTLY Beta(½, ½+S) on f_g — the term IS the Bayes update the
        # reference was standing in for, which is why it needs no separate normalisation.
        term_fg=_Beta.ppf(0.5, 0.5, 0.5 + S),
    )


def _mass_quantile(x, w, f):
    o = np.argsort(x)
    return float(x[o][np.searchsorted(np.cumsum(w[o]) / w.sum(), f)])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default=str(DEFAULT_SUITE))
    ap.add_argument("--index", default=str(DEFAULT_INDEX))
    ap.add_argument("--cache", default=None, help="default: <suite>/oracle_cache")
    ap.add_argument("--conditions", nargs="*", default=None)
    args = ap.parse_args()

    suite = Path(args.suite)
    cache = Path(args.cache) if args.cache else suite / "oracle_cache"
    index = TranscriptIndex.load(args.index)
    cfg = PipelineConfig()
    conds = args.conditions or sorted(p.name for p in cache.iterdir() if p.is_dir())

    print(f"{'condition':<42} {'nBOUNDARY':>7} {'mass':>10} {'q med':>6} {'q<.9':>6} "
          f"{'mwae term':>9} {'mwae ref':>8} {'Δ':>8} {'harm n':>7} {'harm%':>6}")
    print("-" * 120)
    agg, per_boundary = [], {}
    for tag in conds:
        try:
            m = measure(_load(cache, suite, tag, index, cfg))
        except Exception as e:  # noqa: BLE001
            print(f"{tag:<42} SKIP  {type(e).__name__}: {e}")
            continue
        if m["S"].size == 0:
            print(f"{tag:<42} — no BOUNDARY carries boundary_spliced")
            continue
        M, t = m["M"], m["true_fg"]
        mw = lambda x: float(np.sum(M * np.abs(x - t)) / M.sum())  # noqa: E731
        term, ref = mw(m["term_fg"]), mw(np.full_like(t, 0.5))
        harm = t > 0.5
        print(f"{tag:<42} {m['S'].size:>7,} {M.sum():>10,.0f} "
              f"{_mass_quantile(m['q'], M, 0.5):>6.3f} "
              f"{float(M[m['q'] < 0.9].sum() / M.sum()):>6.3f} "
              f"{term:>9.4f} {ref:>8.4f} {term - ref:>+8.4f} {int(harm.sum()):>7,} "
              f"{100 * M[harm].sum() / M.sum():>5.1f}%")
        agg.append((tag, term - ref))
        per_boundary[tag] = (m["live"], m["q"])

    if not agg:
        print("\nnothing measured — is the oracle cache present?")
        return 1

    print("\n" + "=" * 120)
    worse = [a for a in agg if a[1] > 0]
    print(f"⛔ the raw-count term is WORSE than the uninformative reference on "
          f"{len(worse)} of {len(agg)} conditions")
    print(f"   worst {max(agg, key=lambda a: a[1])[0]}  {max(a[1] for a in agg):+.4f}")
    print(f"   best  {min(agg, key=lambda a: a[1])[0]}  {min(a[1] for a in agg):+.4f}")
    for axis in ("capture_off", "capture_on", "ss_0.50", "ss_0.99"):
        sub = [a for a in agg if axis in a[0]]
        if sub:
            print(f"   {axis:<12} mean Δ {np.mean([a[1] for a in sub]):+.4f}   "
                  f"worse on {sum(1 for a in sub if a[1] > 0)}/{len(sub)}")

    # ⭐ IS q A PROPERTY OF THE RNA GEOMETRY RATHER THAN OF THE COMPOSITION? Two conditions differing
    #   ONLY in gDNA share the same transcripts, so a q that moves with gDNA would mean q is not
    #   geometry — and then no annotation-derived divisor could ever recover it.
    pair = [t for t in ("gdna_g00_ss_0.50_nrna_mid_capture_off",
                        "gdna_g50_ss_0.50_nrna_mid_capture_off") if t in per_boundary]
    if len(pair) == 2:
        (la, qa), (lb, qb) = per_boundary[pair[0]], per_boundary[pair[1]]
        both = la & lb
        xa, xb = qa[both[la]], qb[both[lb]]
        if xa.size > 2:
            from scipy.stats import spearmanr

            r = float(spearmanr(xa, xb).statistic)
            print(f"\n⭐ q is GEOMETRY, not composition: q(g00) vs q(g50) on {xa.size:,} shared BOUNDARIES "
                  f"— Spearman ρ = {r:.4f}")
            print("   ⇒ the missing quantity is the unbroken EXONIC reach either side of the boundary, which "
                  "is\n     annotation-derivable in principle and does NOT exist as an array today.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
