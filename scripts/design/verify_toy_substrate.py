#!/usr/bin/env python
"""⭐⭐⭐ IS THE INPUT TO THE SOLVER CORRECT? — the toy substrate, end to end, against first principles.

⛔ **NO SOLVER RUNS HERE.** Every number below is either (a) derived by hand from the transcript's own
geometry, (b) read off the simulator's per-fragment ground truth in the BAM read names, or (c) read off
the accumulator payload. Nothing is compared against a previous run and nothing is tolerated.

**The generative model, read out of the simulator rather than assumed** (`wgs_engine`, `sampling`):

1. each fragment's LENGTH ``w`` is drawn from ``Normal(mean, sd)`` truncated to ``[min, max]``;
2. its START is drawn over the ``seq_len - w + 1`` legal starts of its own template — uniformly with
   capture off, and from the probe landscape with capture on (`capture.sampler.sample_starts`);
3. the read name records ``{t_id}:{start}-{start+w}:{strand}:{i}`` in **TEMPLATE** coordinates —
   spliced-transcript space for mRNA, pre-mRNA space for nascent, genomic for gDNA.

So for a two-exon transcript with exons of ``e1`` and ``e2`` bases and spliced length ``L = e1 + e2``,
a length-``w`` fragment has exactly ``L - w + 1`` legal starts, of which

    w - 1                  cross the sj        (start in [e1-w+1, e1-1])
    (e1 - w + 1)+          lie inside exon 1
    (e2 - w + 1)+          lie inside exon 2

and those three add back to ``L - w + 1``. ⭐ That identity is the spine of the whole file: it is checked
per length, so a placement bug cannot hide in an aggregate.

**The deposit rule, transcribed from `tests/native/_accumulator_reference.py`** (the executable
specification, which the C++ is gated byte-identical to) and re-implemented here from the TRUTH rather
than from the BAM — so this compares *the molecule the simulator made* with *the tally the accumulator
holds*, which no existing gate does:

* the path is the fragment's contiguous genomic SEGMENTS; a sj splits it in two;
* ``region_start_count`` gets +1 at the region containing the path's FIRST covered base — one per accepted
  fragment, so its sum is the deposited count;
* a region_bound boundary strictly inside a segment is CROSSED: ``boundary_spliced`` if the path used any annotated
  sj, else ``boundary_unspliced``;
* a region between two consecutively-crossed boundaries OF ONE SEGMENT is SPANNED;
* every annotated sj used gets ``sj_count`` +1;
* CONTAINED — ``region_contained`` — iff the path used no sj AND its first and last bases are in
  the same region. ⛔ A spliced fragment is never contained.

⭐⭐ **What that predicts for a ONE-transcript toy, and it is worth stating before measuring it:** on a
pure-mRNA two-exon transcript whose exons ARE the regions, every fragment deposits either one
``region_contained`` (in exon 1 or exon 2) or one ``sj_count`` — and ``boundary_unspliced`` at the two
``intron|exon`` boundaries is **exactly zero**, because mature RNA cannot cross an exon↔intron boundary
contiguously. ⛔ **That last clause is NOT "zero everywhere"**, and the difference is the whole reason the
structural gate is now a SET: as soon as another transcript's exon spans the intron, mature RNA crosses
that boundary legitimately. On `splice_both_strands` TB+ and TC− both do, at all six interior BOUNDARIES.

⭐⭐⭐ **ANY NUMBER OF TRANSCRIPTS, ON EITHER STRAND (2026-08-05).** The predecessor refused a ``−``
transcript and refused more than one, which made the rung the splice-flux reframe must be derived against
unverifiable input (TRAPS: self-checking-validator). Both restrictions are gone, and the mechanism is one idea: work in
``u``-space, the transcript coordinate RE-ORIENTED to ascend with the genome. See :class:`Geom`. The
reflections are transcribed from `rigel.sim.bam.transcript_to_genomic_blocks` and
``premrna_to_genomic_interval`` rather than assumed, per TRAPS: two-divisors-opposite-sign's rule about reading the
simulator's code.

⛔⛔ **WHAT IT MEASURED, 2026-08-05.** Every bank matched the truth-derived deposit **exactly** on
`spliced_exons` × 4 arms (pure mRNA, + gDNA, + nascent, capture-ON at 120 kb) and on
**`splice_both_strands` × 2 arms** — pure mRNA at 100 k, and + gDNA + nascent(60) — with **zero**
fragments held on any of them. On `splice_both_strands`, 4 transcripts / 2 strands / 2 sj: 22
nonzero banks all Δ = 0, ``sj_count`` 174 and 190 against truth 174 and 190, and both sj'
crossing counts within 0.5 z of the exact placement prediction.

⭐⭐ **The finding worth carrying:** ``sj_count`` equals the TRUE sj-crossing count exactly, and it
still does with FOUR overlapping transcripts on both strands, where a fragment's sj set has to
resolve against its own transcript. The implicit-splice path is not lossy here — which is what makes this
toy a clean substrate, and what says a pool shortfall on the real panel is about ANNOTATION AMBIGUITY,
not read length.

⚠ **Perturbations — TRAPS: self-checking-validator's second half, and they are a ``--perturb`` FLAG rather than a hand
edit so every claim here is reproducible.** All six fire — ⛔⛔ **ON `splice_both_strands`, AND THE SPEC
IS PART OF THE CLAIM.** Re-verified 2026-08-13. On the owner's `spliced_exons` rung only `drop_sj` fires
and the other five are SILENT — not a hole, but not a pass either: that spec declares ONE `+`-strand
transcript, so `pos_blocks` has no `−` fragment to mirror, `single_geom`/`transcript_order` have no
second geometry to confuse, and the structural-set gate `unspliced_zero_everywhere` flips is VACUOUS
there (no exon spans a BOUNDARY, so both its sets are empty). ⭐ Reading those five silences as a pass
is `TRAPS: could-the-arm-have-fired`, and it was read that way once. `pos_premrna` additionally needs
`--nrna`; the arms below are `--spec splice_both_strands --no-gdna [--nrna 60]`.

| perturbation | what fired |
|---|---|
| `pos_blocks` | ⭐ the deposit gate — a ``−`` mature fragment mapped as ``+`` lands mirror-image |
| `pos_premrna` | ⭐ the deposit gate — same for nascent, against the transcript's genomic SPAN |
| `transcript_order` | the deposit gate, and 2,586 fragments stop mapping at all |
| `single_geom` | the deposit gate — the old one-transcript assumption loses HALF the deposits |
| `unspliced_zero_everywhere` | ⭐⭐ the structural-set gate, naming all six BOUNDARIES. ⚠ Only on the
  pure-mRNA arm, because that is the only arm the claim is about — a scope statement, not a hole |
| `drop_sj` | the sj-axis equality AND the deposit gate |

⚠ **Two gates are RESOLUTION-AWARE and say so out loud** rather than passing quietly. The
drawn-vs-realised length discrimination needs the templates to be short enough to tilt the marginal: on
`splice_both_strands` the four templates (2/8/10/4 kb) separate the two laws by only 1.10 bp, so the arm
needs ≳ 83 k mRNA fragments before that check has teeth, and it prints the required depth when it does
not (TRAPS: key-on-a-realised-quantity). ⭐ The start-uniformity test is POOLED across lengths by standardising each
fragment against its own ``(eff−1)/2 ± √((eff²−1)/12)``, which is exact and has teeth at any depth — the
per-length form it replaces needed n ≥ 500 in ONE length and therefore reported "not checked" on every
toy, which reads as a pass.
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pysam  # noqa: E402

DESIGN = Path(__file__).resolve().parent
_s = importlib.util.spec_from_file_location("toy_harness", DESIGN / "toy_harness.py")
TH = importlib.util.module_from_spec(_s)
sys.modules["toy_harness"] = TH
_s.loader.exec_module(TH)

from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.sim.read_name import parse_origin  # noqa: E402

SUITE = Path.home() / "Downloads/rigel_runs/suite/ladder"
INDEX = Path.home() / "Downloads/rigel_runs/suite/rigel_index"
TYPE_NAMES = {0: "intergenic", 1: "intron", 2: "exon"}

FAIL: list[str] = []

#: ⭐⭐ THE PERTURBATION UNDER TEST — TRAPS: self-checking-validator's second half, as a flag rather than a hand edit so
#: every claim below is reproducible by anyone. ``--perturb <name>`` corrupts exactly one thing and the
#: run must then FAIL; a perturbation that fires nothing is a hole in the gate set, not a pass.
#: The verdicts are in this module's docstring.
PERTURB = "none"
PERTURBATIONS = {
    "none": "the honest run",
    "pos_blocks": "map a − transcript's mature interval as if it were + (drop the t → L−t reflection)",
    "pos_premrna": "map a − transcript's NASCENT interval as if it were + (drop the reflection)",
    "transcript_order": "take a − transcript's exons in TRANSCRIPT order, not genomic — so its "
    "sj u-region_bound and its exon lengths swap",
    "single_geom": "route every RNA fragment through the FIRST transcript's geometry — the "
    "one-transcript assumption this file used to make",
    "unspliced_zero_everywhere": "assert the OLD structural claim (boundary_unspliced == 0 everywhere on "
    "a pure-mRNA library), which is true only of a one-transcript two-exon toy",
    "drop_sj": "hide one of the spec's introns from the sj map, so its crossings read as "
    "contiguous",
}


def check(ok: bool, label: str, detail: str = "") -> bool:
    print(f"   {'PASS' if ok else '⛔ FAIL'}  {label}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAIL.append(label)
    return ok


# ──────────────────────────────────────────────────────────────────────────────────────────────────
# the transcript's own geometry — everything predicted below comes from these three numbers
# ──────────────────────────────────────────────────────────────────────────────────────────────────


@dataclasses.dataclass(frozen=True)
class Geom:
    """One transcript's geometry, on EITHER strand.

    ⭐⭐ **Everything below works in ``u``-space — the transcript coordinate RE-ORIENTED so it ascends
    with the genome** — and that single change is what makes one implementation serve both strands.
    ``rigel.sim.bam.transcript_to_genomic_blocks`` maps a ``−``-strand transcript interval by reflecting
    it (``t → L − t``) and then walking the exons in ascending genomic order, so::

        u = t                on ``+``          u_interval([a,b)) = [a, b)
        u = L − t            on ``−``          u_interval([a,b)) = [L−b, L−a)

    and in ``u``-space the exons are the genomic ones in genomic order, the sj region_bounds are their
    cumulative lengths, and a uniform draw over ``t`` is a uniform draw over ``u``. So placement
    combinatorics, containment and per-sj crossing are all strand-FREE once expressed here.

    ⛔ The reflection is transcribed from the simulator's own code, not assumed — TRAPS: two-divisors-opposite-sign's rule.
    """

    exons: tuple[tuple[int, int], ...]  #: genomic [start, end) per exon, ascending
    strand: str

    @property
    def exon_lengths(self):
        """Genomic-ascending exon lengths — i.e. ``u``-space, on both strands."""
        out = tuple(e - s for s, e in self.exons)
        if PERTURB == "transcript_order" and self.strand == "-":
            return out[::-1]
        return out

    @property
    def spliced_length(self) -> int:
        return sum(self.exon_lengths)

    @property
    def offsets(self):
        """Cumulative ``u``-space offset of each exon's first base (genomic order)."""
        out, acc = [], 0
        for length in self.exon_lengths:
            out.append(acc)
            acc += length
        return tuple(out)

    @property
    def region_bounds_u(self) -> tuple[int, ...]:
        """The sj' ``u``-space offsets, ascending — one per intron, genomic order."""
        return self.offsets[1:]

    @property
    def introns(self) -> tuple[tuple[int, int], ...]:
        """Genomic ``(intron_start, intron_end)`` per intron, ascending — the sj axis's own key."""
        return tuple(
            (self.exons[i][1], self.exons[i + 1][0]) for i in range(len(self.exons) - 1)
        )

    @property
    def span(self) -> tuple[int, int]:
        """The transcript's genomic span — the pre-mRNA template nascent RNA is drawn from."""
        return self.exons[0][0], self.exons[-1][1]

    def u_interval(self, t_start: int, t_end: int) -> tuple[int, int]:
        """A TRANSCRIPT-space interval re-oriented to ascend with the genome (see the class docstring)."""
        if self.strand == "-" and PERTURB != "pos_blocks":
            L = self.spliced_length
            return L - t_end, L - t_start
        return t_start, t_end

    def blocks(self, t_start: int, t_end: int) -> list[tuple[int, int]]:
        """Map a TRANSCRIPT-space interval to its genomic segments, ascending, on either strand.

        ⚠ Byte-for-byte the same arithmetic as ``rigel.sim.bam.transcript_to_genomic_blocks``: reflect
        on ``−``, then walk the exons in genomic order."""
        u0, u1 = self.u_interval(t_start, t_end)
        out = []
        for (gs, ge), off in zip(self.exons, self.offsets):
            length = ge - gs
            lo, hi = max(u0, off), min(u1, off + length)
            if hi > lo:
                out.append((gs + (lo - off), gs + (hi - off)))
        return out

    def premrna_block(self, t_start: int, t_end: int) -> tuple[int, int]:
        """A NASCENT (pre-mRNA-space) interval mapped to its single genomic segment, either strand.

        ⚠ The same reflection, against the transcript's genomic SPAN rather than its spliced length —
        ``rigel.sim.bam.premrna_to_genomic_interval``. Getting this wrong on ``−`` puts every nascent
        fragment at the mirror-image position inside the gene, which is invisible in a total and
        visible per region."""
        g0, g1 = self.span
        n = g1 - g0
        if self.strand == "-" and PERTURB != "pos_premrna":
            t_start, t_end = n - t_end, n - t_start
        return g0 + t_start, g0 + t_end

    def n_placements(self, w: int) -> tuple[int, int, list[int]]:
        """``(total, crossing-any, per-exon contained)`` legal starts for a length-``w`` fragment."""
        total = max(self.spliced_length - w + 1, 0)
        inside = [max(e - w + 1, 0) for e in self.exon_lengths]
        return total, total - sum(inside), inside

    def crossings_per_sj(self, w: int) -> list[int]:
        """Legal starts that cross EACH sj, in ``u``-space — one entry per intron.

        A length-``w`` fragment crosses the region_bound at ``c`` iff its ``u``-start lies in
        ``[c − w + 1, c − 1]``, intersected with the ``[0, L − w]`` legal range. ⭐ Per sj, not
        pooled: on a 3-exon transcript a long fragment can cross two, and only the per-sj form
        can tell a mis-placed sj from a mis-counted one."""
        total = max(self.spliced_length - w + 1, 0)
        if total <= 0:
            return [0] * len(self.region_bounds_u)
        return [
            max(min(c - 1, total - 1) - max(0, c - w + 1) + 1, 0) for c in self.region_bounds_u
        ]

    def exon_of_u(self, u0: int, u1: int) -> int | None:
        """The exon index wholly containing ``[u0, u1)``, or ``None`` if it crosses a sj."""
        for k, off in enumerate(self.offsets):
            if u0 >= off and u1 <= off + self.exon_lengths[k]:
                return k
        return None


# ──────────────────────────────────────────────────────────────────────────────────────────────────
# an INDEPENDENT implementation of the deposit rule, driven by TRUTH
# ──────────────────────────────────────────────────────────────────────────────────────────────────


class TruthTally:
    """Re-derive every accumulator bank from the simulator's own per-fragment truth.

    ⚠ Deliberately NOT a copy of the production code path: it takes the TRUE genomic segments (from the
    read name, mapped through the annotation) rather than the aligned records, so a disagreement with
    the payload is a real statement about fidelity and not a tautology (TRAPS: self-checking-validator)."""

    def __init__(
        self, region_bounds: np.ndarray, n_regions: int, n_boundaries: int, sj_id_by_intron: dict
    ):
        self.region_bounds = np.asarray(region_bounds, np.int64)
        self.region_contained = np.zeros(n_regions, np.int64)
        self.region_start = np.zeros(n_regions, np.int64)
        self.boundary_unspliced = np.zeros(n_boundaries, np.int64)
        self.boundary_spliced = np.zeros(n_boundaries, np.int64)
        #: ⛔ TWO DIFFERENT THINGS, AND THEY ONCE SHARED THE NAME ``sj`` — the crossing TALLY below
        #: was overwritten by the lookup map, so ``deposit`` raised on the first spliced fragment.
        #: This file's own `TRAPS: two-masks-one-name`; keep the tally and the map named apart.
        self.sj = Counter()  #: jid -> crossings, the bank ``sj_count`` is scored against
        self.sj_id_by_intron = sj_id_by_intron  #: (intron_start, intron_end) -> jid
        self.n_deposited = 0

    def _region_of(self, pos: int) -> int:
        """⚠ ``region_bounds`` carries BOTH reference boundaries — ``n_regions = n_region_bounds − 1`` — so the
        region index is one less than the insertion point. Transcribed from the reference's own
        ``_local_region``; getting it wrong shifts every deposit by one region, which is exactly what a
        first run of this file did."""
        return min(max(int(np.searchsorted(self.region_bounds, pos, side="right")) - 1, 0), self.region_bounds.size - 2)

    def deposit(self, segments: list[tuple[int, int]], introns: list[tuple[int, int]]):
        if not segments:
            return
        sj_ids = [self.sj_id_by_intron[i] for i in introns if i in self.sj_id_by_intron]
        spliced = bool(sj_ids)
        first_base, last_base = segments[0][0], segments[-1][1] - 1
        self.region_start[self._region_of(first_base)] += 1
        self.n_deposited += 1
        # ⛔ THE ARRAY AND THE LOOP INDEX MUST NOT SHARE A NAME. They did — ``boundary[boundary - 1]``
        # raised ``TypeError: 'int' object is not subscriptable`` the moment a fragment crossed one.
        # The bank carries the ``_count`` suffix exactly as `tests/native/_accumulator_reference.py`
        # names it, which is the file this one transcribes.
        boundary_count = self.boundary_spliced if spliced else self.boundary_unspliced
        for seg_start, seg_end in segments:
            first = int(np.searchsorted(self.region_bounds, seg_start, side="right"))
            last = int(np.searchsorted(self.region_bounds, seg_end, side="left"))
            for boundary in range(first, last):
                boundary_count[boundary - 1] += 1
        # ⚠ ``boundary`` indexes ``region_bounds``; boundary ``boundary-1`` follows the reference exactly.
        for jid in sj_ids:
            self.sj[jid] += 1
        if not sj_ids and self._region_of(first_base) == self._region_of(last_base):
            self.region_contained[self._region_of(first_base)] += 1


# ──────────────────────────────────────────────────────────────────────────────────────────────────


def _drawn_pmf(donor):
    """The pre-capture length marginal the sampler actually draws, as ``(lengths, pmf)``.

    ⚠ The engine rounds the config to INTEGERS and its sampler does ``.astype(int)``, which FLOORS a
    continuous normal — so ``P(w) = Phi((w+1−μ)/σ) − Phi((w−μ)/σ)``, not the density at ``w``."""
    lo, hi = donor.frag_min, donor.frag_max
    mu_i, sd_i = float(round(donor.frag_mean)), float(round(donor.frag_std))
    xs = np.arange(0, hi + 2)
    cdf = 0.5 * (1.0 + np.vectorize(math.erf)((xs - mu_i) / (sd_i * math.sqrt(2.0))))
    p = np.diff(cdf, append=1.0)[: xs.size]
    p[(xs < lo) | (xs > hi)] = 0.0
    return xs, p / p.sum()


def _pool_length_gate(label, widths, templates, donor, tag):
    """The realised length marginal of ONE pool, against ``f_post(w) ∝ f_pre(w)·Σ_t a_t·eff_t(w)``.

    ⭐⭐ **THE REALISED MARGINAL IS NOT THE DRAWN ONE**, and with several templates it is not one
    template's either. `wgs_engine._post_capture_length_allocation` reweights the pre-capture draw by
    the pool's TOTAL opportunity summed over templates — ``Σ_t a_t·(L_t − w + 1)+`` off capture — because
    a library cannot yield more fragments of a length than its templates have placements for. Each pool
    (mature / nascent) has its own template lengths, so each is gated separately: pooling them would let
    a long template's tilt cancel a short one's.

    ``templates`` is ``[(abundance, template_length), …]``. Returns nothing; checks in place."""
    w = np.asarray(widths, float)
    if not w.size or not templates:
        return
    lo, hi = donor.frag_min, donor.frag_max
    check(int(w.min()) >= lo and int(w.max()) <= hi,
          f"every {label} length is inside the configured [{lo}, {hi}]",
          f"observed [{int(w.min())}, {int(w.max())}]")
    xs, p = _drawn_pmf(donor)
    mu_drawn = float((xs * p).sum())
    eff = np.zeros(xs.size, float)
    for a, tlen in templates:
        eff += float(a) * np.maximum(tlen - xs + 1, 0).astype(float)
    if donor.capture_on:
        # ⭐ Under capture ``total_eff(w)`` is the PROBE-WEIGHTED opportunity, and the engine's own note
        # says capture SELECTS FOR LENGTH: a longer fragment presents more sequence, overlaps more of a
        # probe, and is captured better — until the overlap saturates at the probe length. With
        # `off_target_weight` a and `binding_per_base` b and probes tiling the exons, a fragment's best
        # single-probe overlap is ``min(w, probe_length)``, so
        #     total_eff(w) ≈ Σ_t a_t (L_t − w + 1)+ · (a + b·min(w, probe_length))
        # ⚠ APPROXIMATE — it ignores the sj split and the per-probe tiling phase — so the gate
        # below only asserts the DIRECTION under capture, never the value.
        k = donor.capture_knobs
        a0, b0 = float(k["off_target_weight"]), float(k["binding_per_base"])
        plen = float(k["probe_length"])
        eff = eff * (a0 + b0 * np.minimum(xs.astype(float), plen))
    q = p * eff
    if q.sum() <= 0:
        return
    q = q / q.sum()
    mu_pred = float((xs * q).sum())
    sd_pred = float(np.sqrt((q * (xs - mu_pred) ** 2).sum()))
    se = sd_pred / math.sqrt(w.size)
    print(f"      {label}: drawn mean {mu_drawn:.2f} · realised prediction {mu_pred:.2f} "
          f"· observed {w.mean():.2f}   (n = {w.size:,}, se = {se:.2f})")
    print(f"         templates {tag}")
    if donor.capture_on:
        check(w.mean() > mu_drawn,
              f"⭐ capture SELECTS FOR LENGTH: the realised {label} mean exceeds the drawn mean",
              f"observed {w.mean():.2f} vs drawn {mu_drawn:.2f}; the approximate probe-weighted "
              f"law predicts {mu_pred:.2f}")
        return
    check(abs(w.mean() - mu_pred) < 4 * se,
          f"the realised {label} length mean matches f_pre(w)·Σ_t a_t·eff_t(w)",
          f"observed {w.mean():.2f} vs predicted {mu_pred:.2f}  (4 se = {4 * se:.2f})")
    # ⭐⭐ DOES THIS SPEC RESOLVE THE TWO LAWS AT ALL? The separation the geometry offers is
    # ``|mu_pred − mu_drawn|``, and it is a property of the TEMPLATE LENGTHS, not of the depth: a
    # 10,000 bp template has a nearly flat ``(L − w + 1)`` over the fragment range, so a mixture
    # dominated by long templates tilts the marginal by ~1 bp where a 2 kb one tilts it by ~5.
    # ⛔ So an unconditional "it must be 4 se from the drawn mean" FAILS on a long-template spec for a
    # reason that has nothing to do with fidelity. The honest form is `suite_resolves.py`'s: assert the
    # discrimination where the design supplies it, and where it does not, say so and say what depth
    # would (TRAPS: key-on-a-realised-quantity — key the claim on a realised quantity, and print it).
    gap = abs(mu_pred - mu_drawn)
    if gap > 4 * se:
        check(abs(w.mean() - mu_drawn) > 4 * se,
              f"⭐ and the {label} pool is DISTINGUISHABLE from the drawn marginal, so it has teeth",
              f"drawn {mu_drawn:.2f} is {abs(w.mean() - mu_drawn) / se:.1f} se away")
    else:
        need = int(math.ceil((4.0 * sd_pred / max(gap, 1e-9)) ** 2)) if gap > 0 else -1
        print(f"      ⚠ NOT RESOLVED: the two laws differ by only {gap:.2f} bp here against 4 se = "
              f"{4 * se:.2f}, so this arm cannot tell them apart.")
        print(f"        ⭐ the accuracy gate above still holds; the discrimination needs "
              f"n ≳ {need:,} {label} fragments (template lengths, not depth, set the gap).")
    print(f"         length sd: observed {w.std():.2f} vs predicted {sd_pred:.2f}")


def gate_simulator(frags, geoms, abund, nrna_abund, spec, donor):
    """S1 — did the simulator produce what it was asked for, in the coordinates it claims?"""
    print("\n── GATE S1: THE SIMULATOR ────────────────────────────────────────────────────────────")
    mrna = [f for f in frags if f["kind"] == "mrna"]
    gdna = [f for f in frags if f["kind"] == "gdna"]
    nrna = [f for f in frags if f["kind"] == "nrna"]
    print(f"   fragments in the oracle BAM: {len(frags):,}   "
          f"mRNA {len(mrna):,} · nascent {len(nrna):,} · gDNA {len(gdna):,}")
    # ⚠ ``n_rna_fragments`` is the RNA BUDGET, not the mature count: with ``nrna_abundance > 0`` the
    # nascent pool is drawn from the same budget. Checking it against the mature count alone reads as a
    # simulator failure and is a reader failure.
    check(len(mrna) + len(nrna) == spec.n_rna_fragments,
          "mRNA + nascent count == the RNA budget requested",
          f"{len(mrna):,} + {len(nrna):,} vs {spec.n_rna_fragments:,}")
    if spec.nrna_abundance:
        print(f"      nascent share of the RNA budget: {len(nrna) / max(len(mrna) + len(nrna), 1):.1%}")

    # ⛔ every t_id in the BAM must be one this spec declared — a typo'd id would otherwise be silently
    # skipped by every per-transcript prediction below and read as "nothing to check".
    unknown = sorted({f["t_id"] for f in mrna + nrna} - set(geoms))
    check(not unknown, "every RNA fragment's transcript id is one the spec declared", f"{unknown}")

    for tid, geom in sorted(geoms.items()):
        L = geom.spliced_length
        own = [f for f in mrna if f["t_id"] == tid]
        bad = [f for f in own if not (0 <= f["start"] and f["end"] <= L and f["end"] > f["start"])]
        check(not bad, f"{tid}{geom.strand}: every mRNA interval lies inside [0, {L:,}) of "
                       "TRANSCRIPT space", f"{len(bad)} violation(s) of {len(own):,}")
    if nrna:
        for tid, geom in sorted(geoms.items()):
            g0, g1 = geom.span
            n = g1 - g0
            own = [f for f in nrna if f["t_id"] == tid]
            bad = [f for f in own if not (0 <= f["start"] and f["end"] <= n and f["end"] > f["start"])]
            check(not bad, f"{tid}{geom.strand}: every nascent interval lies inside [0, {n:,}) of "
                           "PRE-mRNA space", f"{len(bad)} violation(s) of {len(own):,}")

    # ⭐ the two RNA pools are drawn SEPARATELY (`wgs_engine._accumulate_pool` per pool), so each gets
    # its own opportunity-reweighted prediction against its own template lengths.
    _pool_length_gate(
        "mRNA", [f["end"] - f["start"] for f in mrna],
        [(abund[t], g.spliced_length) for t, g in geoms.items() if abund.get(t, 0.0) > 0],
        donor,
        " · ".join(f"{t}:{g.spliced_length:,}bp×{abund[t]:g}"
                   for t, g in sorted(geoms.items()) if abund.get(t, 0.0) > 0),
    )
    if nrna:
        _pool_length_gate(
            "nascent", [f["end"] - f["start"] for f in nrna],
            [(nrna_abund[t], g.span[1] - g.span[0])
             for t, g in geoms.items() if nrna_abund.get(t, 0.0) > 0],
            donor,
            " · ".join(f"{t}:{g.span[1] - g.span[0]:,}bp×{nrna_abund[t]:g}"
                       for t, g in sorted(geoms.items()) if nrna_abund.get(t, 0.0) > 0),
        )
    return mrna, nrna, gdna


def gate_splice_combinatorics(mrna, geoms, capture_on):
    """S2 — the sj-crossing counts, against the exact placement count, PER TRANSCRIPT and PER
    SJ.

    ⭐ Per sj rather than pooled: on a multi-exon transcript a long fragment can cross two, and
    only the per-sj form separates "a sj is in the wrong place" from "the total is off".
    Returns ``{(t_id, sj_index): observed crossings}``."""
    print("\n── GATE S2: THE SPLICE, AGAINST FIRST-PRINCIPLES COMBINATORICS ───────────────────────")
    out: dict[tuple[str, int], int] = {}
    spliced = {t: g for t, g in geoms.items() if len(g.exons) > 1}
    if not spliced:
        print("   no multi-exon transcript on this spec — nothing to check")
        return out
    for tid, geom in sorted(spliced.items()):
        own = [f for f in mrna if f["t_id"] == tid]
        by_len = Counter(f["end"] - f["start"] for f in own)
        us = [geom.u_interval(f["start"], f["end"]) for f in own]
        print(f"\n   {tid}{geom.strand}   L = {geom.spliced_length:,}   exon lengths "
              f"{geom.exon_lengths} (genomic order)   n = {len(own):,}")
        for j, region_bound in enumerate(geom.region_bounds_u):
            obs = sum(1 for u0, u1 in us if u0 < region_bound < u1)
            out[(tid, j)] = obs
            exp = var = 0.0
            for w, n in by_len.items():
                total = max(geom.spliced_length - w + 1, 0)
                if total <= 0:
                    continue
                p = geom.crossings_per_sj(w)[j] / total
                exp += n * p
                var += n * p * (1 - p)
            sd = math.sqrt(var)
            gpos = geom.introns[j]
            print(f"      sj {j} @ genomic {gpos[0]:,}→{gpos[1]:,}  (u-region_bound {region_bound:,}): "
                  f"observed {obs:,} · uniform prediction {exp:,.1f} ± {sd:,.1f}")
            if capture_on:
                continue
            check(abs(obs - exp) < 4 * sd,
                  f"{tid}: sj {j}'s crossing count matches uniform placement",
                  f"z = {(obs - exp) / max(sd, 1e-9):+.2f}")
        if capture_on:
            print("      ⚠ capture is ON, so uniform placement is NOT the null — the probe landscape")
            print("        reweights the starts. Reported for scale only; GATE S3 tests the law.")
        # ⛔ THE PARTITION IDENTITY, per length, so a placement bug cannot hide in an aggregate: a
        # fragment lies wholly inside exactly one exon, or it crosses at least one sj.
        worst = (0.0, None)
        for w, n in sorted(by_len.items()):
            if n < 200:
                continue
            total, crossing, inside = geom.n_placements(w)
            if total <= 0:
                continue
            got = Counter()
            for u0, u1 in us:
                if u1 - u0 != w:
                    continue
                k = geom.exon_of_u(u0, u1)
                got["cross" if k is None else f"e{k + 1}"] += 1
            expect = [("cross", crossing / total)] + [
                (f"e{k + 1}", inside[k] / total) for k in range(len(inside))
            ]
            for label, exp_p in expect:
                sd_l = math.sqrt(n * exp_p * (1 - exp_p))
                z = (got[label] - n * exp_p) / max(sd_l, 1e-9)
                if abs(z) > abs(worst[0]):
                    worst = (z, f"w={w} {label}: {got[label]} vs {n * exp_p:.1f}")
        print(f"      worst per-length z over lengths with n ≥ 200: {worst[0]:+.2f}   ({worst[1]})")
        if not capture_on:
            check(abs(worst[0]) < 4.5, f"{tid}: no single fragment length is mis-placed")
        n_cross = sum(1 for u0, u1 in us if geom.exon_of_u(u0, u1) is None)
        n_in = len(own) - n_cross
        check(n_cross + n_in == len(own),
              f"{tid}: every mRNA fragment is EITHER inside one exon OR crossing ≥1 sj",
              f"{n_in:,} contained + {n_cross:,} crossing = {len(own):,}")
    return out


def _uniform_start_z(starts, widths, template_len):
    """The POOLED z for "starts are uniform over each fragment's own legal range", exactly.

    ⭐ A length-``w`` fragment has ``eff = template_len − w + 1`` legal starts, so under uniform
    placement ``E[start] = (eff−1)/2`` and ``Var[start] = (eff²−1)/12`` — both known in closed form and
    both length-dependent. Standardising EACH fragment by its OWN mean and sd gives i.i.d. mean-0
    variance-1 terms, so ``Σz/√n`` is a single N(0,1) test over the whole pool.

    ⛔ This replaces a per-length test that needed n ≥ 500 IN ONE LENGTH: a realistic draw spreads
    ~450 distinct lengths over the pool, so that form reported "not checked" on every toy and the
    absence read as a pass. Pooling costs nothing and has teeth at any depth.

    Returns ``(z, n_used)``; ``n_used`` excludes degenerate fragments (``eff ≤ 1``, no freedom)."""
    zs = []
    for s, w in zip(starts, widths):
        eff = template_len - w + 1
        if eff <= 1:
            continue
        var = (eff * eff - 1) / 12.0
        zs.append((s - (eff - 1) / 2.0) / math.sqrt(var))
    if not zs:
        return 0.0, 0
    arr = np.asarray(zs, float)
    return float(arr.mean() * math.sqrt(arr.size)), arr.size


def gate_capture(mrna, nrna, gdna, geoms, res, donor, spec):
    """S3 — does the realised start distribution match the sampler's OWN weight law?"""
    print("\n── GATE S3: HYBRID CAPTURE ───────────────────────────────────────────────────────────")
    if not donor.capture_on:
        print("   capture is OFF on this donor — the null is uniform placement, tested in GATE S2.")
        # ⛔ and that null must actually hold, or 'off' is not off. Per transcript: a mixture of
        # templates of different lengths is NOT uniform on any one of them.
        for tid, geom in sorted(geoms.items()):
            own = [f for f in mrna if f["t_id"] == tid]
            z, n = _uniform_start_z([f["start"] for f in own],
                                    [f["end"] - f["start"] for f in own], geom.spliced_length)
            print(f"   {tid}{geom.strand}: pooled mRNA start-uniformity z = {z:+.2f}  (n = {n:,})")
            if n:
                check(abs(z) < 4.0, f"{tid}: mRNA starts are uniform on the transcript, capture off")
            nas = [f for f in nrna if f["t_id"] == tid]
            if nas:
                g0, g1 = geom.span
                z, n = _uniform_start_z([f["start"] for f in nas],
                                        [f["end"] - f["start"] for f in nas], g1 - g0)
                print(f"   {tid}{geom.strand}: pooled nascent start-uniformity z = {z:+.2f} "
                      f"(n = {n:,}, pre-mRNA template {g1 - g0:,} bp)")
                check(abs(z) < 4.0, f"{tid}: nascent starts are uniform on the pre-mRNA, capture off")
        # ⭐⭐ AND THE ONE THE REFRAME DERIVATION LEANS ON: gDNA is uniform along the chromosome before
        # capture. the truth's own gDNA-density RATIO is the comparand for a gDNA
        # level precisely because of this, so it is worth a gate rather than an assumption.
        if gdna:
            z, n = _uniform_start_z([f["start"] for f in gdna],
                                    [f["end"] - f["start"] for f in gdna], spec.genome_length)
            print(f"   gDNA: pooled start-uniformity z = {z:+.2f}  (n = {n:,} over "
                  f"{spec.genome_length:,} bp)")
            check(abs(z) < 4.0,
                  "⭐⭐ gDNA starts are UNIFORM along the chromosome with capture off — the premise "
                  "the gDNA-level comparand rests on")
        return
    print(f"   capture ON.  knobs: {donor.capture_knobs}")
    print("   ⭐ the ON-vs-OFF contrast, which is where the capture LAW is actually tested, lives in")
    print("      `verify_capture.py`. What follows is the shape this run alone can show.")
    # ⭐ The probes tile each exon, so on this spec EVERY transcript base is under some probe and the
    # discriminating quantity is not "on probe vs off" but the SPLIT penalty at the sj: a probe
    # never spans the sj (they tile per exon), so a sj-crossing fragment's best single
    # probe group covers at most the longer of its two overhangs.
    for tid, geom in sorted(geoms.items()):
        if len(geom.exons) < 2:
            continue
        by_len = defaultdict(lambda: [0, 0])
        for f in mrna:
            if f["t_id"] != tid:
                continue
            u0, u1 = geom.u_interval(f["start"], f["end"])
            by_len[u1 - u0][0] += 1
            by_len[u1 - u0][1] += 1 if geom.exon_of_u(u0, u1) is None else 0
        print(f"\n   {tid}{geom.strand}"
              f"\n   {'w':>6} {'n':>8} {'crossing':>9} {'observed p':>11} {'uniform p':>10} {'ratio':>7}")
        rows = 0
        for w in sorted(by_len):
            n, c = by_len[w]
            if n < 300:
                continue
            total, crossing, _ = geom.n_placements(w)
            if total <= 0:
                continue
            pu = crossing / total
            po = c / n
            print(f"   {w:>6} {n:>8,} {c:>9,} {po:>11.4f} {pu:>10.4f} {po / max(pu, 1e-9):>7.3f}")
            rows += 1
            if rows >= 8:
                break
    print("   ⭐ ratio < 1 means capture DEPLETES sj-crossing fragments relative to uniform,")
    print("      which is what per-exon probe tiling predicts: no probe spans the sj, so a")
    print("      crossing fragment has less overlap with its best single probe than a contained one.")


def gate_accumulator(frags, geoms, res, payload, ra, spec):
    """A1 — every bank, re-derived from the TRUTH, against the payload. The main event."""
    print("\n── GATE A1: THE ACCUMULATOR, AGAINST TRUTH-DERIVED DEPOSITS ──────────────────────────")
    index = res.index
    region_bounds = np.asarray(payload.region_bounds, np.int64)
    n_regions, n_boundaries = int(payload.n_regions), int(payload.n_boundaries)
    # the sj map: (intron_start, intron_end) -> jid, from the index's own sj axis
    from rigel.calibration.splice_graph import build_sj_geometry_arrays

    jg = build_sj_geometry_arrays(index)
    starts = np.asarray(ra.start, np.int64)
    sizes = np.asarray(ra.region_size_bp, np.int64)
    sj = {}
    for k in range(int(getattr(jg, "n_sj", 0))):
        src = int(np.asarray(jg.src_region)[k])
        dst = int(np.asarray(jg.dst_region)[k])
        sj[(int(starts[src] + sizes[src]), int(starts[dst]))] = k
    # ⛔ EVERY spec-declared intron must appear on the index's sj axis, and nothing else may. A
    # sj the map misses would make its crossings read as contiguous and silently deposit on the
    # wrong bank; the equality catches both directions.
    if PERTURB == "drop_sj" and sj:
        sj.pop(sorted(sj)[0])
    declared = {i for g in geoms.values() for i in g.introns}
    check(declared == set(sj),
          "⭐ the index's sj axis is EXACTLY the set of introns the spec declares",
          f"spec {sorted(declared)} vs index {sorted(sj)}")

    tt = TruthTally(region_bounds, n_regions, n_boundaries, sj)

    first_tid = sorted(geoms)[0]
    for f in frags:
        # the one-transcript assumption, as a perturbation: route everything through one geometry
        tid = first_tid if PERTURB == "single_geom" else f.get("t_id")
        if f["kind"] == "mrna":
            segs = geoms[tid].blocks(f["start"], f["end"])
            introns = [(segs[i][1], segs[i + 1][0]) for i in range(len(segs) - 1)]
        elif f["kind"] == "nrna":
            segs, introns = [geoms[tid].premrna_block(f["start"], f["end"])], []
        else:
            segs, introns = [(f["start"], f["end"])], []
        segs = [(a, b) for a, b in segs if b > a]
        tt.deposit(segs, introns)

    def col(bank):
        return np.asarray(bank, np.int64).sum(axis=1)

    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    print(f"\n   {'object':<30} {'bank':<16} {'TRUTH':>9} {'accumulator':>12} {'Δ':>8}")
    print("   " + "-" * 80)
    ok_all = True
    pay = {
        # ⛔ NO `region_spanning` ROW, AND THE TRUTH-SIDE TALLY IS GONE WITH IT. The bank was DELETED
        # from the payload on 2026-08-08 (`5591cc01`, "six dead banks gone"); this table went on
        # naming it, so the gate raised `AttributeError` before it compared anything. A `src/`
        # deletion is the mechanism that kills an instrument — there is nothing left to score it
        # against, so keeping the tally would be scoring a bank that no longer exists.
        "region_contained": (col(payload.region_contained_count), tt.region_contained, "region"),
        "region_start": (np.asarray(payload.region_start_count, np.int64), tt.region_start, "region"),
        "boundary_unspliced": (col(payload.boundary_unspliced_count), tt.boundary_unspliced, "boundary"),
        "boundary_spliced": (col(payload.boundary_spliced_count), tt.boundary_spliced, "boundary"),
    }
    for i in range(n_regions):
        label = f"REGION {TYPE_NAMES[int(rtype[i])]} [{starts[i]:,},{starts[i] + sizes[i]:,})"
        for bank, (got, exp, axis) in pay.items():
            if axis != "region" or (got[i] == 0 and exp[i] == 0):
                continue
            d = int(got[i]) - int(exp[i])
            ok_all &= d == 0
            print(f"   {label:<30} {bank:<16} {int(exp[i]):>9,} {int(got[i]):>12,} {d:>+8,}")
    for i in range(n_boundaries):
        label = f"BOUNDARY @{starts[i + 1]:,}" if i + 1 < n_regions else f"BOUNDARY #{i}"
        for bank, (got, exp, axis) in pay.items():
            if axis != "boundary" or (got[i] == 0 and exp[i] == 0):
                continue
            d = int(got[i]) - int(exp[i])
            ok_all &= d == 0
            print(f"   {label:<30} {bank:<16} {int(exp[i]):>9,} {int(got[i]):>12,} {d:>+8,}")
    sj_got = col(payload.sj_count)
    for k in range(sj_got.shape[0]):
        d = int(sj_got[k]) - int(tt.sj.get(k, 0))
        ok_all &= d == 0
        print(f"   {'SJ boundary #' + str(k):<30} {'sj_count':<16} "
              f"{tt.sj.get(k, 0):>9,} {int(sj_got[k]):>12,} {d:>+8,}")
    print()
    check(ok_all, "⭐⭐ EVERY BANK MATCHES THE TRUTH-DERIVED DEPOSIT EXACTLY")

    # ── ⭐⭐ THE ZERO THAT IS STRUCTURAL, PREDICTED FROM THE ANNOTATION ALONE ──────────────────────
    # On a pure-mRNA library an BOUNDARY can carry unspliced flux only where some transcript's EXON spans
    # its position — mature RNA cannot cross an exon↔intron boundary contiguously (TRAPS: mature-rna-never-crosses-a-boundary). ⛔ The
    # predecessor asserted `boundary_unspliced == 0` EVERYWHERE, which is true only of a one-transcript
    # two-exon toy and false the moment another transcript's exon spans the intron: on
    # `splice_both_strands` TB+ and TC− both do. So the prediction is now a SET, derived from the exon
    # intervals with no reference to the fragments, and the gate is an EQUALITY in both directions —
    # a spurious nonzero and a missing one are both failures (TRAPS: self-checking-validator's "emit all classes").
    e_un = col(payload.boundary_unspliced_count)
    mrna_only = all(f["kind"] == "mrna" for f in frags)
    if mrna_only:
        exonic = {
            int(starts[i + 1])
            for i in range(n_boundaries)
            if i + 1 < n_regions
            and any(s < int(starts[i + 1]) < e for g in geoms.values() for s, e in g.exons)
        }
        if PERTURB == "unspliced_zero_everywhere":
            exonic = set()
        got = {int(starts[i + 1]) for i in range(n_boundaries)
               if i + 1 < n_regions and int(e_un[i]) > 0}
        print(f"   BOUNDARIES an exon spans contiguously: {sorted(exonic)}")
        print(f"   BOUNDARIES with nonzero boundary_unspliced: {sorted(got)}")
        check(got == exonic,
              "⛔ pure-mRNA library: boundary_unspliced is nonzero EXACTLY where an exon spans the BOUNDARY",
              f"only-observed {sorted(got - exonic)} · only-predicted {sorted(exonic - got)}")
    return tt


def gate_invariants(payload, tt, frags):
    """A2 — the accumulator's own externally-checkable identities."""
    print("\n── GATE A2: THE ACCUMULATOR'S OWN INVARIANTS ─────────────────────────────────────────")
    qc = {f.name: getattr(payload.qc, f.name) for f in dataclasses.fields(payload.qc)}
    start_sum = int(np.asarray(payload.region_start_count, np.int64).sum())
    dep_len_sum = int(np.asarray(payload.deposited_lengths, np.int64).sum())
    print(f"   Σ region_start_count      {start_sum:,}")
    print(f"   Σ deposited_lengths     {dep_len_sum:,}")
    print(f"   fragments in the BAM    {len(frags):,}")
    print(f"   TRUTH-derived deposits  {tt.n_deposited:,}")
    check(start_sum == dep_len_sum,
          "Σ region_start_count == Σ deposited_lengths (one per accepted fragment)")
    held = len(frags) - start_sum
    print(f"   ⭐ fragments NOT deposited: {held:,}  ({held / max(len(frags), 1):.2%})")
    for k, v in sorted(qc.items()):
        if isinstance(v, (int, float)) and v:
            print(f"      qc[{k}] = {v:,}")
    return held


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spec", default="spliced_exons")
    ap.add_argument("--donor", default="gdna_g50_ss_0.50_nrna_none_capture_off")
    ap.add_argument("--n-rna", type=int, default=None)
    ap.add_argument("--nrna", type=float, default=None)
    ap.add_argument("--genome-length", type=int, default=0)
    ap.add_argument("--no-gdna", action="store_true",
                    help="⭐ the cleanest arm: a pure-mRNA library, where every prediction is exact")
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/rigel_verify_substrate"))
    ap.add_argument("--perturb", default="none", choices=sorted(PERTURBATIONS),
                    help="⭐⭐ corrupt ONE thing and confirm a gate fires (TRAPS: self-checking-validator). "
                         + " · ".join(f"{k}: {v}" for k, v in PERTURBATIONS.items()))
    args = ap.parse_args()
    global PERTURB
    PERTURB = args.perturb
    if PERTURB != "none":
        print(f"⚠⚠ PERTURBATION ACTIVE — {PERTURB}: {PERTURBATIONS[PERTURB]}")
        print("   this run MUST fail. A perturbation that fires nothing is a HOLE in the gate set.")

    index = TranscriptIndex.load(str(INDEX))
    config = dataclasses.replace(CalibrationConfig(), calib_refit_iters=0)
    donor = TH.harvest(SUITE / args.donor, index, config=config)
    if args.no_gdna:
        donor = dataclasses.replace(donor, gdna_rate_per_base=0.0)
    spec = TH.SPECS[args.spec]
    if args.genome_length:
        spec = dataclasses.replace(spec, genome_length=int(args.genome_length))
    if args.n_rna:
        spec = dataclasses.replace(spec, n_rna_fragments=int(args.n_rna))
    if args.nrna is not None:
        spec = dataclasses.replace(spec, nrna_abundance=float(args.nrna))
    spec = dataclasses.replace(spec, name=f"verify_{spec.name}")

    # ⭐ ANY number of transcripts, on EITHER strand — `Geom` works in genome-ascending ``u``-space.
    geoms, abund, nrna_abund = {}, {}, {}
    for gene in spec.genes:
        for t in gene["transcripts"]:
            tid = str(t["t_id"])
            if tid in geoms:
                print(f"duplicate transcript id {tid!r} in the spec", file=sys.stderr)
                return 2
            exons = tuple(tuple(int(x) for x in e) for e in sorted(t["exons"]))
            geoms[tid] = Geom(exons, str(gene["strand"]))
            abund[tid] = float(t.get("abundance", 0.0))
            # ⚠ ``spec.nrna_abundance`` OVERRIDES the per-transcript value when > 0 — that is
            # `Scenario.build_oracle`'s own rule, and reading the per-transcript field instead would
            # predict a zero nascent pool on every toy that has one.
            nrna_abund[tid] = (
                float(spec.nrna_abundance) if spec.nrna_abundance > 0
                else float(t.get("nrna_abundance", 0.0))
            )

    print("=" * 100)
    print(f"⭐⭐⭐ SUBSTRATE VERIFICATION — {spec.name}")
    print("=" * 100)
    for tid, geom in sorted(geoms.items()):
        print(f"   transcript {tid}{geom.strand}  exons {list(geom.exons)}  "
              f"spliced length {geom.spliced_length:,}  abundance {abund[tid]:g}"
              + (f"  introns {list(geom.introns)}" if geom.introns else ""))
    print(f"   chromosome {spec.genome_length:,} bp   nascent {spec.nrna_abundance:g}   "
          f"donor {args.donor}")
    print(f"   gDNA rate {donor.gdna_rate_per_base:.6g}/bp"
          + ("   ⭐ FORCED TO ZERO (pure-mRNA arm)" if args.no_gdna else ""))

    # simulate + scan, reusing the harness so this is the SAME substrate every other instrument sees
    sub = _simulate(spec, donor, args.work_dir)
    frags = _fragments(sub["bam"])
    mrna, nrna, gdna = gate_simulator(frags, geoms, abund, nrna_abund, spec, donor)
    gate_splice_combinatorics(mrna, geoms, donor.capture_on)
    gate_capture(mrna, nrna, gdna, geoms, sub["res"], donor, spec)
    tt = gate_accumulator(frags, geoms, sub["res"], sub["payload"], sub["ra"], spec)
    gate_invariants(sub["payload"], tt, frags)

    print()
    print("=" * 100)
    if FAIL:
        print(f"⛔ {len(FAIL)} GATE(S) FAILED:")
        for f in FAIL:
            print(f"     - {f}")
    else:
        print("✅ EVERY GATE PASSED")
    print("=" * 100)
    return 1 if FAIL else 0


def _simulate(spec, donor, work_dir):
    """`toy_harness.run_toy`'s first half, kept verbatim so this verifies the SHIPPED path."""
    import rigel.pipeline as PL
    from rigel.sim import CaptureConfig, GDNAConfig, ReadSimConfig, Scenario

    pc = PipelineConfig()
    wd = Path(work_dir) / spec.name
    wd.mkdir(parents=True, exist_ok=True)
    n_gdna = int(round(donor.gdna_rate_per_base * spec.genome_length))
    sim_cfg = ReadSimConfig(
        frag_mean=int(round(donor.frag_mean)), frag_std=int(round(donor.frag_std)),
        frag_min=donor.frag_min, frag_max=donor.frag_max, read_length=donor.read_length,
        strand_specificity=donor.strand_specificity, seed=spec.seed,
    )
    gdna_cfg = GDNAConfig(abundance=0.0, frag_mean=int(round(donor.frag_mean)),
                          frag_std=int(round(donor.frag_std)))
    sc = Scenario(spec.name, genome_length=spec.genome_length, seed=spec.seed, work_dir=wd / "sim")
    for gene in spec.genes:
        sc.add_gene(gene["gene_id"], gene["strand"], gene["transcripts"])
    capture_cfg = None
    if donor.capture_on:
        probes = TH._toy_probes(spec, wd / "probes.tsv", donor.capture_knobs)
        k = donor.capture_knobs
        capture_cfg = CaptureConfig(
            probes=probes, probe_format="transcript",
            off_target_weight=float(k["off_target_weight"]),
            binding_per_base=float(k["binding_per_base"]),
            gdna_split_penalty=float(k["gdna_split_penalty"]),
            min_overlap=int(k["min_overlap"]),
        )
    res = sc.build_oracle(
        n_rna_fragments=int(spec.n_rna_fragments),
        gdna_fraction=n_gdna / max(spec.n_rna_fragments, 1),
        nrna_abundance=float(spec.nrna_abundance),
        sim_config=sim_cfg, gdna_config=gdna_cfg, capture_config=capture_cfg,
    )
    bam = str(res.bam_path)
    scan = dataclasses.replace(pc.scan, sj_strand_tag=PL._native_detect_sj_tag(bam))
    _st, _sm, _b, payload = PL.scan_and_buffer(bam, res.index, scan)
    return {"res": res, "bam": bam, "payload": payload, "ra": RegionArrays.from_index(res.index)}


def _fragments(bam: str) -> list[dict]:
    """One record per QNAME, with the simulator's own ground truth decoded from the name."""
    out, seen = [], set()
    with pysam.AlignmentFile(bam, "rb") as fh:
        for rec in fh.fetch(until_eof=True):
            q = rec.query_name
            if q in seen:
                continue
            seen.add(q)
            o = parse_origin(q)
            out.append({"qname": q, "kind": o.kind, "t_id": o.transcript_id,
                        "start": o.start, "end": o.end, "strand": o.strand})
    return out


if __name__ == "__main__":
    raise SystemExit(main())
