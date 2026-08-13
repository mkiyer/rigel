#!/usr/bin/env python
"""⭐⭐ IS HYBRID CAPTURE DOING WHAT IT CLAIMS? — the same toy, probes ON and OFF, nothing else moved.

Capture is a physical claim with a direction: **probes on the exons enrich the sequence that binds them
and thereby DEPLETE, relatively, the sequence that does not.** That is testable without a solver and
without re-implementing the sampler — simulate the identical chromosome, identical seed, identical gDNA
rate and identical RNA budget twice, once with probes and once without, and read the three consequences
off the ground truth:

1. ⭐ **gDNA per base, per region.** gDNA is uniform along the genome before capture, so its post-capture
   density profile IS the capture landscape, measured directly. Exon regions must go UP, intergenic and
   intronic regions must go DOWN — relative to the same library's own mean.
2. ⭐ **the mRNA length marginal.** A longer fragment presents more sequence to a probe, so capture
   selects for length until the overlap saturates at the probe length. The realised mean must RISE.
3. ⭐ **the junction-crossing share.** ``_toy_probes`` tiles probes WITHIN each exon so none spans the
   junction, so a junction-crossing fragment's best single probe covers only its longer overhang. Below
   about ``2 x probe_length`` that is less overlap than a contained fragment gets, so crossing fragments
   must be relatively DEPLETED.

⛔ Each is a direction the knobs predict, not a number this file fits. The magnitudes are reported
because they are the useful part; the gates are on sign.

⭐⭐ **MEASURED 2026-08-05** on `spliced_exons`, 120 kb, gDNA rate raised to 0.05/bp for power (both arms):

* probed exon regions **65-66x enriched**; the off-probe INTERIOR (116 kb) **0.046x**, i.e. 22x depleted;
* ⛔ but a whole 7 kb intron REGION reads **0.80x**, NOT depleted — because the +-500 bp COLLAR abutting a
  probed exon is itself **5.4x ENRICHED**. A region beside a probe is not off-probe, and an intron's
  measured gDNA density under capture is a MIXTURE of a depleted interior and two enriched ends;
* the mRNA length mean rises **+6.6 bp** (z = 14.8) — capture selects for length, as the engine says;
* the junction depletion has the hard onset the weight law predicts: crossing/uniform is **0.62x below
  120 bp** (one probe length), 0.97 in 120-240, and 0.97 above 240 where both saturate. ⚠ The residual
  3 % above 240 bp is the probe TILING PHASE: 1,000 bp of exon at 120 bp per probe leaves a 40 bp runt
  at the exon end, and that runt is what a junction-crossing fragment's overhang lands on.
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import os
import sys
from collections import Counter
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import math  # noqa: E402
import numpy as np  # noqa: E402

SCR = Path(__file__).resolve().parent
_s = importlib.util.spec_from_file_location("vts", SCR / "verify_toy_substrate.py")
V = importlib.util.module_from_spec(_s)
sys.modules["vts"] = V
_s.loader.exec_module(V)

from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

TYPE_NAMES = {0: "intergenic", 1: "intron", 2: "exon"}


def run(spec, donor, work_dir, tag):
    sub = V._simulate(dataclasses.replace(spec, name=f"{spec.name}_{tag}"), donor, work_dir)
    frags = V._fragments(sub["bam"])
    return sub, frags


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spec", default="spliced_exons")
    ap.add_argument("--donor", default="gdna_g50_ss_0.50_nrna_none_capture_on")
    ap.add_argument("--n-rna", type=int, default=40000)
    ap.add_argument("--genome-length", type=int, default=120000)
    ap.add_argument("--gdna-rate", type=float, default=None,
                    help="⚠ override the donor's gDNA density PURELY FOR STATISTICAL POWER. Both arms "
                         "get the same value, so it cannot bias the ON/OFF contrast; at the donor's "
                         "own rate a 7 kb intron holds ~10 fragments and no ratio is resolvable.")
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/rigel_verify_capture"))
    args = ap.parse_args()

    index = TranscriptIndex.load(str(V.INDEX))
    cfg = dataclasses.replace(CalibrationConfig(), calib_refit_iters=0)
    donor_on = V.TH.harvest(V.SUITE / args.donor, index, config=cfg)
    if not donor_on.capture_on:
        print("pick a capture_on donor", file=sys.stderr)
        return 2
    # ⛔ ONE thing varied. Same gDNA rate, same lengths, same strand, same seed — only the probes go.
    if args.gdna_rate is not None:
        donor_on = dataclasses.replace(donor_on, gdna_rate_per_base=float(args.gdna_rate))
    donor_off = dataclasses.replace(donor_on, capture_on=False)

    spec = V.TH.SPECS[args.spec]
    spec = dataclasses.replace(spec, genome_length=args.genome_length,
                               n_rna_fragments=args.n_rna, name=f"cap_{args.spec}")
    tset = [t for g in spec.genes for t in g["transcripts"]]
    geom = V.Geom(tuple(tuple(e) for e in tset[0]["exons"]), spec.genes[0]["strand"])

    print("=" * 104)
    print(f"⭐⭐ HYBRID CAPTURE — {spec.name}, probes ON vs OFF, everything else identical")
    print("=" * 104)
    k = donor_on.capture_knobs
    print(f"   knobs: off_target_weight {k['off_target_weight']}  binding_per_base "
          f"{k['binding_per_base']}  probe_length {k['probe_length']}  "
          f"gdna_split_penalty {k['gdna_split_penalty']}")
    print(f"   chromosome {spec.genome_length:,} bp   gDNA rate {donor_on.gdna_rate_per_base:.6g}/bp"
          f"   RNA budget {spec.n_rna_fragments:,}")

    sub_on, fr_on = run(spec, donor_on, args.work_dir, "on")
    sub_off, fr_off = run(spec, donor_off, args.work_dir, "off")

    ra = sub_on["ra"]
    starts = np.asarray(ra.start, np.int64)
    sizes = np.asarray(ra.region_size_bp, np.int64)
    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)

    def gdna_per_region(frags):
        c = Counter()
        for f in frags:
            if f["kind"] != "gdna":
                continue
            mid = (f["start"] + f["end"]) // 2
            i = int(np.searchsorted(starts, mid, side="right")) - 1
            if 0 <= i < starts.size:
                c[i] += 1
        return c

    g_on, g_off = gdna_per_region(fr_on), gdna_per_region(fr_off)
    n_on = sum(g_on.values()) or 1
    n_off = sum(g_off.values()) or 1
    print("\n── 1. gDNA DENSITY PER REGION (fragments per kb, normalised to each library's own total) ──")
    print(f"   total gDNA fragments: probes ON {n_on:,}   OFF {n_off:,}")
    print(f"\n   {'region':<34} {'bp':>9} {'n OFF':>7} {'n ON':>7} {'OFF /kb':>9} {'ON /kb':>9} "
          f"{'ratio':>8} {'±':>7}")
    fails = []
    for i in range(starts.size):
        bp = int(sizes[i])
        if bp <= 0:
            continue
        d_off = g_off.get(i, 0) / n_off / bp * 1000
        d_on = g_on.get(i, 0) / n_on / bp * 1000
        label = f"{TYPE_NAMES[int(rtype[i])]} [{starts[i]:,},{starts[i] + bp:,})"
        k_on, k_off = g_on.get(i, 0), g_off.get(i, 0)
        ratio = d_on / d_off if d_off > 0 else float("inf")
        # ⛔ a ratio of two Poisson counts: log-sd is sqrt(1/k_on + 1/k_off). With ten fragments a side
        # that is 45 %, so a 1.3x reading means nothing and must not be gated on.
        rel = math.sqrt(1.0 / max(k_on, 1) + 1.0 / max(k_off, 1))
        resolvable = k_on >= 25 and k_off >= 25
        print(f"   {label:<34} {bp:>9,} {k_off:>7,} {k_on:>7,} {d_off:>9.4f} {d_on:>9.4f} "
              f"{ratio:>8.2f} {'±' + f'{100 * rel:.0f}%':>7} {'' if resolvable else ' (too few)'}")
        if not resolvable:
            continue
        if TYPE_NAMES[int(rtype[i])] == "exon" and not (math.log(ratio) > 2 * rel):
            fails.append(f"exon {label} not ENRICHED (ratio {ratio:.2f} ± {100 * rel:.0f}%)")
        # ⛔ NOT "every non-exon region is depleted". A region ABUTTING a probed exon is not off-probe: a
        # fragment lying in its first or last ~fragment-length can still overlap the neighbour's probe
        # and be captured. Only the INTERIOR is off-probe, and it is split out below.
    V.check(not fails, "⭐ every probed exon region is ENRICHED", "; ".join(fails))

    # ── the interior/edge split, which is where "off-probe is depleted" is actually testable ──
    print("\n   ⭐ A region beside a probed exon is NOT off-probe. Split each long region into the")
    print("      collar within one max-fragment-length of a probed exon, and the interior beyond it:")
    collar = int(donor_on.frag_max)
    exon_bounds = [(s, e) for s, e in geom.exons]

    def near_probe(pos):
        return any(s - collar < pos < e + collar for s, e in exon_bounds)

    def split_counts(frags):
        edge = interior = 0
        for f in frags:
            if f["kind"] != "gdna":
                continue
            mid = (f["start"] + f["end"]) // 2
            if any(s <= mid < e for s, e in exon_bounds):
                continue
            if near_probe(mid):
                edge += 1
            else:
                interior += 1
        return edge, interior

    e_on, i_on = split_counts(fr_on)
    e_off, i_off = split_counts(fr_off)
    # the bp behind each population, so these are densities and not raw counts
    total_bp = int(spec.genome_length)
    exon_bp = sum(e - s for s, e in exon_bounds)
    collar_bp = sum(min(e + collar, total_bp) - max(s - collar, 0) for s, e in exon_bounds) - exon_bp
    interior_bp = total_bp - exon_bp - collar_bp
    for label, k_on, k_off, bp in (("probe COLLAR (±%d bp)" % collar, e_on, e_off, collar_bp),
                                   ("off-probe INTERIOR", i_on, i_off, interior_bp)):
        d_on = k_on / n_on / max(bp, 1) * 1000
        d_off = k_off / n_off / max(bp, 1) * 1000
        r = d_on / d_off if d_off > 0 else float("inf")
        rel = math.sqrt(1.0 / max(k_on, 1) + 1.0 / max(k_off, 1))
        print(f"   {label:<34} {bp:>9,} {k_off:>7,} {k_on:>7,} {d_off:>9.4f} {d_on:>9.4f} "
              f"{r:>8.2f} {'±' + f'{100 * rel:.0f}%':>7}")
        if label.startswith("off-probe"):
            V.check(math.log(r) < -2 * rel,
                    "⭐ the OFF-PROBE INTERIOR is depleted by capture",
                    f"ratio {r:.3f} ± {100 * rel:.0f}%")
        else:
            V.check(math.log(r) > 2 * rel,
                    "⭐ and the COLLAR beside a probe is ENRICHED — which is why a whole intron REGION "
                    "reads ~1.0",
                    f"ratio {r:.2f} ± {100 * rel:.0f}%")

    print("\n── 2. THE mRNA LENGTH MARGINAL — capture selects for length ──────────────────────────")
    w_on = np.array([f["end"] - f["start"] for f in fr_on if f["kind"] == "mrna"])
    w_off = np.array([f["end"] - f["start"] for f in fr_off if f["kind"] == "mrna"])
    se = float(np.sqrt(w_on.var() / w_on.size + w_off.var() / w_off.size))
    print(f"   mean fragment length   OFF {w_off.mean():.2f}   ON {w_on.mean():.2f}   "
          f"Δ {w_on.mean() - w_off.mean():+.2f} bp   (se of the difference {se:.2f})")
    V.check(w_on.mean() - w_off.mean() > 4 * se,
            "⭐ capture SELECTS FOR LENGTH — the realised mean rises",
            f"z = {(w_on.mean() - w_off.mean()) / max(se, 1e-9):+.1f}")

    print("\n── 3. THE JUNCTION-CROSSING SHARE — probes tile PER EXON, so none spans the junction ──")
    e1 = geom.exon_lengths[0]

    def share(frags):
        m = [f for f in frags if f["kind"] == "mrna"]
        c = sum(1 for f in m if f["start"] < e1 < f["end"])
        return c, len(m), c / max(len(m), 1)

    c_on, n_m_on, s_on = share(fr_on)
    c_off, n_m_off, s_off = share(fr_off)
    # ⚠ the two arms have different LENGTH marginals (that is result 2), and the crossing share depends
    # on length — so compare against each arm's OWN uniform-placement expectation, which removes it.
    def expected_share(frags):
        m = [f for f in frags if f["kind"] == "mrna"]
        by = Counter(f["end"] - f["start"] for f in m)
        exp = sum(n * geom.n_placements(w)[1] / max(geom.n_placements(w)[0], 1)
                  for w, n in by.items())
        return exp / max(len(m), 1)

    e_on, e_off = expected_share(fr_on), expected_share(fr_off)
    print(f"   probes OFF   crossing {c_off:,}/{n_m_off:,} = {s_off:.4f}   "
          f"uniform-placement expectation {e_off:.4f}   ratio {s_off / e_off:.3f}")
    print(f"   probes ON    crossing {c_on:,}/{n_m_on:,} = {s_on:.4f}   "
          f"uniform-placement expectation {e_on:.4f}   ratio {s_on / e_on:.3f}")
    V.check(abs(s_off / e_off - 1.0) < 0.05,
            "with probes OFF the crossing share matches uniform placement")

    # ⭐⭐ THE EFFECT IS CONFINED, AND THE WEIGHT LAW SAYS EXACTLY WHERE. A contained fragment's best
    # single-probe overlap is ``min(w, probe_length)``; a crossing one's is
    # ``min(max(a, w−a), probe_length)`` with ``a`` its exon-1 overhang. For ``w >= 2·probe_length`` the
    # longer overhang already exceeds the probe, so BOTH saturate at ``probe_length`` and there is NO
    # depletion at all. Below that the crossing fragment loses overlap, and below ``probe_length`` it
    # loses the most. ⛔ So a flat "crossing fragments are depleted" is the WRONG prediction — the right
    # one is a gradient with a hard onset at 2·probe_length.
    plen = int(k["probe_length"])
    print(f"\n   {'length band':<22} {'OFF obs/exp':>12} {'ON obs/exp':>12} {'ON/OFF':>9} "
          f"{'n(ON)':>8}   predicted")
    bands = [(0, plen, "STRONG depletion"), (plen, 2 * plen, "mild depletion"),
             (2 * plen, 10 ** 9, "NONE — both saturate")]
    for lo, hi, pred in bands:
        def band_ratio(frags):
            m = [f for f in frags if f["kind"] == "mrna" and lo <= f["end"] - f["start"] < hi]
            if not m:
                return float("nan"), 0
            obs = sum(1 for f in m if f["start"] < e1 < f["end"]) / len(m)
            by = Counter(f["end"] - f["start"] for f in m)
            exp = sum(n * geom.n_placements(w)[1] / max(geom.n_placements(w)[0], 1)
                      for w, n in by.items()) / len(m)
            return (obs / exp if exp > 0 else float("nan")), len(m)
        r_on, n1 = band_ratio(fr_on)
        r_off, _ = band_ratio(fr_off)
        rr = r_on / r_off if r_off and np.isfinite(r_off) else float("nan")
        print(f"   {f'{lo}-{min(hi, 500)} bp':<22} {r_off:>12.3f} {r_on:>12.3f} {rr:>9.3f} "
              f"{n1:>8,}   {pred}")
    def band(frags, lo, hi):
        m = [f for f in frags if f["kind"] == "mrna" and lo <= f["end"] - f["start"] < hi]
        obs = sum(1 for f in m if f["start"] < e1 < f["end"]) / max(len(m), 1)
        by = Counter(f["end"] - f["start"] for f in m)
        exp = sum(n * geom.n_placements(w)[1] / max(geom.n_placements(w)[0], 1)
                  for w, n in by.items()) / max(len(m), 1)
        return obs / exp if exp > 0 else float("nan")
    V.check(band(fr_on, 0, plen) < band(fr_off, 0, plen),
            f"⭐ below one probe length ({plen} bp) crossing fragments ARE depleted",
            f"{band(fr_on, 0, plen):.3f} vs {band(fr_off, 0, plen):.3f}")
    V.check(abs(band(fr_on, 2 * plen, 10 ** 9) - band(fr_off, 2 * plen, 10 ** 9)) < 0.03,
            f"⭐ and above two probe lengths ({2 * plen} bp) there is NO effect, as the law predicts",
            f"{band(fr_on, 2 * plen, 10 ** 9):.3f} vs {band(fr_off, 2 * plen, 10 ** 9):.3f}")

    print()
    print("=" * 104)
    if V.FAIL:
        print(f"⛔ {len(V.FAIL)} GATE(S) FAILED:")
        for f in V.FAIL:
            print(f"     - {f}")
    else:
        print("✅ EVERY CAPTURE GATE PASSED")
    print("=" * 104)
    return 1 if V.FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
