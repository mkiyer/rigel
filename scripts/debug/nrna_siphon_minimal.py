"""nRNA siphon root-cause diagnostic — direct per-locus EM probe.

Bypasses BAM I/O / scoring and feeds hand-crafted log-likelihoods to the
batch locus-EM C++ kernel.  The goal is to isolate the dynamics that
cause pure-gDNA fragments to be siphoned into the nRNA component.

Scenario design
---------------
A single locus on chr1 with two transcripts:
  t0 — annotated mRNA, 5 exons, exonic_length L_exon = 1000 bp,
       genomic span S = 10 000 bp (introns total 9000 bp)
  t1 — synthetic nRNA, single exon spanning the full pre-mRNA span,
       length S = 10 000 bp, ``is_nrna=True``, ``is_synthetic=True``

We then synthesise per-fragment "units" with hand-set log-likelihoods
(no real scoring).  Each unit represents one intronic gDNA-origin
fragment, with two candidates:
  * nRNA (t1) at ``log_lik_nrna``
  * gDNA   at ``gdna_log_lik`` (already pre-corrected per-fragment)

By varying the relative likelihood and the locus-prior
(``alpha_gdna`` / ``alpha_rna``) we can tease apart whether the siphon
is driven by:
  (A) likelihood degeneracy (nRNA ≈ gDNA per-fragment),
  (B) prior misassignment (alpha_gdna under-counts contamination),
  (C) strand-asymmetric exploitation (sense fragments preferentially
      assigned to nRNA when gDNA's strand prior is uniform).

Usage::

    python scripts/debug/nrna_siphon_minimal.py
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

# Make the rigel package importable when running as a script.
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from rigel.config import EMConfig  # noqa: E402
from rigel.estimator import AbundanceEstimator  # noqa: E402
from rigel.locus_partition import partition_and_free  # noqa: E402
from rigel.scored_fragments import Locus, ScoredFragments  # noqa: E402
from rigel.splice import SpliceStrandCol  # noqa: E402

UNSPLICED_SENSE = int(SpliceStrandCol.UNSPLICED_SENSE)
UNSPLICED_ANTI = int(SpliceStrandCol.UNSPLICED_ANTISENSE)


# ---------------------------------------------------------------------------
# Mock index that mimics what TranscriptIndex provides to the EM solver
# ---------------------------------------------------------------------------


class _MockIndex:
    """Two-transcript locus: t0 = mRNA, t1 = synthetic nRNA."""

    def __init__(self, span: int = 10_000, exonic_len: int = 1000):
        self.num_transcripts = 2
        # `length` is the value used by the C++ apply_bias_correction_uniform
        # for each transcript-component:
        #   mRNA → exonic length (effective sampling space)
        #   nRNA → genomic span  (synthetic single-exon spanning pre-mRNA)
        self.t_df = pd.DataFrame(
            {
                "t_id": ["t0", "t1_nrna"],
                "ref": ["chr1", "chr1"],
                "start": [0, 0],
                "end": [span, span],
                "length": np.array([exonic_len, span], dtype=np.int64),
                "is_nrna": np.array([False, True]),
                "is_synthetic": np.array([False, True]),
            }
        )


# ---------------------------------------------------------------------------
# Build a ScoredFragments object with N intronic "gDNA-origin" units
# ---------------------------------------------------------------------------


@dataclass
class FragmentSpec:
    """Per-fragment hand-set state."""

    nrna_log_lik: float            # log P(frag | nRNA component)
    gdna_log_lik: float            # already-corrected log P(frag | gDNA)
    count_col: int = UNSPLICED_SENSE   # SpliceStrandCol


def build_em_data(specs: list[FragmentSpec]) -> tuple[ScoredFragments, list[Locus]]:
    """Construct a ScoredFragments + single Locus from per-frag specs.

    Each fragment has exactly one RNA candidate (the nRNA component, t=1)
    plus one gDNA candidate appended by the C++ extractor when
    ``is_spliced=False`` and ``gdna_log_liks`` is finite.
    """
    n = len(specs)
    offsets = np.zeros(n + 1, dtype=np.int64)
    flat_t = np.empty(n, dtype=np.int32)
    flat_lk = np.empty(n, dtype=np.float64)
    flat_cc = np.empty(n, dtype=np.uint8)

    for i, spec in enumerate(specs):
        flat_t[i] = 1                                   # nRNA component (t1)
        flat_lk[i] = float(spec.nrna_log_lik)
        flat_cc[i] = np.uint8(spec.count_col)
        offsets[i + 1] = offsets[i] + 1

    locus_t = flat_t.copy()
    locus_cc = flat_cc.copy()

    is_spliced = np.zeros(n, dtype=bool)                # all unspliced
    gdna_lk = np.array([s.gdna_log_lik for s in specs], dtype=np.float64)

    em_data = ScoredFragments(
        offsets=offsets,
        t_indices=flat_t,
        log_liks=flat_lk,
        count_cols=flat_cc,
        coverage_weights=np.ones(n, dtype=np.float64),
        tx_starts=np.zeros(n, dtype=np.int32),
        tx_ends=np.full(n, 200, dtype=np.int32),        # frag length 200
        locus_t_indices=locus_t,
        locus_count_cols=locus_cc,
        is_spliced=is_spliced,
        gdna_log_liks=gdna_lk,
        frag_ids=np.arange(n, dtype=np.int64),
        frag_class=np.zeros(n, dtype=np.int8),
        splice_type=np.zeros(n, dtype=np.uint8),
        n_units=n,
        n_candidates=n,
    )

    locus = Locus(
        locus_id=0,
        transcript_indices=np.array([0, 1], dtype=np.int32),
        unit_indices=np.arange(n, dtype=np.int32),
        gdna_span=10_000,
        merged_intervals=[("chr1", 0, 10_000)],
    )
    return em_data, [locus]


# ---------------------------------------------------------------------------
# Run the EM and return the per-component pool counts
# ---------------------------------------------------------------------------


def run_em(
    specs: list[FragmentSpec],
    *,
    alpha_gdna: float,
    alpha_rna: float,
    em_mode: str = "vbem",
    em_iterations: int = 2000,
    span: int = 10_000,
    exonic_len: int = 1000,
) -> dict:
    """Build inputs and run batch locus EM. Returns split into mRNA/nRNA/gDNA."""
    index = _MockIndex(span=span, exonic_len=exonic_len)
    em_data, loci = build_em_data(specs)

    rc = AbundanceEstimator(
        num_transcripts=2,
        em_config=EMConfig(
            seed=42,
            mode=em_mode,
            iterations=em_iterations,
            assignment_mode="fractional",
            n_threads=1,
        ),
        is_nrna=np.array([False, True]),
        is_synthetic=np.array([False, True]),
    )
    rc._transcript_spans = np.array([span, span], dtype=np.float64)
    rc._exonic_lengths = np.array([exonic_len, span], dtype=np.float64)

    partitions = partition_and_free(em_data, loci)
    partition_tuples = [
        (
            p.offsets, p.t_indices, p.log_liks, p.coverage_weights,
            p.tx_starts, p.tx_ends, p.count_cols, p.is_spliced,
            p.gdna_log_liks, p.locus_t_indices, p.locus_count_cols,
        )
        for p in [partitions[i] for i in range(len(loci))]
    ]
    locus_t_lists = [loc.transcript_indices for loc in loci]

    total_gdna, locus_mrna, locus_gdna = rc.run_batch_locus_em_partitioned(
        partition_tuples,
        locus_t_lists,
        np.array([alpha_gdna], dtype=np.float64),
        np.array([alpha_rna], dtype=np.float64),
        index,
        em_iterations=em_iterations,
    )
    rc._gdna_em_total += total_gdna

    # mRNA component (t0) gets nothing here because no fragment lists it.
    em_per_t = rc.em_counts.sum(axis=1)
    return {
        "mrna_t0": float(em_per_t[0]),
        "nrna_t1": float(em_per_t[1]),
        "gdna": float(total_gdna),
        "n_input": len(specs),
    }


# ---------------------------------------------------------------------------
# Helper: bias-corrected log-likelihoods for a synthetic fragment
# ---------------------------------------------------------------------------


def raw_rna_lik(strand_log: float, log_fl: float) -> float:
    """Per-fragment RNA-candidate log-lik *before* bias correction.

    The C++ ``apply_bias_correction_uniform`` adds ``-log(t_length - frag + 1)``
    to every RNA candidate at EM time, so we must NOT bake that in here.
    """
    return strand_log + log_fl


def gdna_lik(log_fl: float, span: int, gdna_flank: int, frag_len: int,
             splice_pen: float = 0.0, log_strand: float = 0.0) -> float:
    """Mimic the per-hit gDNA scoring (NH=1, no NM penalty).

    Unlike the RNA path, gDNA log-liks are *pre-corrected* by the scorer
    (Option B harmonic-mean baked in) and the EM bias correction skips
    the gDNA component (`local_comp >= n_t`).
    """
    L_h = span + gdna_flank
    e_h = max(L_h - frag_len + 1, 1)
    LOG_HALF = math.log(0.5)
    return log_fl + splice_pen + LOG_HALF + log_strand - math.log(e_h)


def effective_log_lik_pair(
    *, log_fl_rna: float, log_fl_gdna: float, log_p_strand_rna: float,
    span: int, frag_len: int, gdna_flank: int = 200,
) -> tuple[float, float, float]:
    """Compute the *effective* log-liks (post C++ bias correction) for
    diagnostic display only — useful to compare what the EM actually sees.
    """
    eff_rna = max(span - frag_len + 1, 1)
    raw_rna = raw_rna_lik(log_p_strand_rna, log_fl_rna)
    em_rna = raw_rna - math.log(eff_rna)
    em_gdna = gdna_lik(log_fl_gdna, span, gdna_flank, frag_len)
    return em_rna, em_gdna, em_rna - em_gdna


# ---------------------------------------------------------------------------
# Experiments
# ---------------------------------------------------------------------------


def _print_header(title: str) -> None:
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")


def _print_row(name: str, res: dict, truth_n: int | None = None) -> None:
    n = res["n_input"]
    g = res["gdna"]
    nrna = res["nrna_t1"]
    mrna = res["mrna_t0"]
    tot = g + nrna + mrna
    extra = ""
    if truth_n is not None:
        extra = f"  | gdna_err={g - truth_n:+.0f}  nrna_siphon={nrna:.0f}/{truth_n}"
    print(
        f"  {name:<40} N={n:>6}  mRNA={mrna:>8.1f}  nRNA={nrna:>8.1f}  gDNA={g:>8.1f}"
        f"  total={tot:>8.1f}{extra}"
    )


def experiment_unstranded_pure_gdna(N: int = 10_000) -> None:
    """Pure-gDNA fragments, unstranded (uniform strand prior)."""
    _print_header("E1 — UNSTRANDED, pure gDNA (intronic, 1 transcript)")
    span = 10_000
    exonic = 1000
    frag_len = 200
    log_fl = -5.0
    LOG_HALF = math.log(0.5)

    raw_nrna = raw_rna_lik(LOG_HALF, log_fl)
    g_lik = gdna_lik(log_fl, span, gdna_flank=200, frag_len=frag_len)
    em_n, em_g, delta = effective_log_lik_pair(
        log_fl_rna=log_fl, log_fl_gdna=log_fl, log_p_strand_rna=LOG_HALF,
        span=span, frag_len=frag_len,
    )
    print(
        f"  raw nRNA={raw_nrna:+.3f}  gDNA={g_lik:+.3f}    "
        f"effective: nRNA(EM)={em_n:+.3f}  gDNA(EM)={em_g:+.3f}  Δ={delta:+.3f}"
    )

    specs = [FragmentSpec(raw_nrna, g_lik, UNSPLICED_SENSE)] * N

    for label, ag, ar in [
        ("alpha_gdna≫alpha_rna  (truth-aware)", 0.99 * N, 0.01 * N),
        ("alpha_gdna=alpha_rna   (uninformative)", 0.5 * N, 0.5 * N),
        ("alpha_gdna≪alpha_rna  (anti-truth)", 0.01 * N, 0.99 * N),
    ]:
        res = run_em(specs, alpha_gdna=ag, alpha_rna=ar, em_mode="vbem",
                     span=span, exonic_len=exonic)
        _print_row(label, res, truth_n=N)


def experiment_stranded_pure_gdna(N: int = 10_000, ss: float = 0.99) -> None:
    """Pure-gDNA fragments under stranded RNA-Seq (50/50 sense:antisense)."""
    _print_header(f"E2 — STRANDED ss={ss}, pure gDNA (intronic, 1 transcript)")
    span = 10_000
    exonic = 1000
    frag_len = 200
    log_fl = -5.0
    log_p_sense = math.log(ss)
    log_p_anti = math.log(1.0 - ss + 1e-300)

    raw_n_sense = raw_rna_lik(log_p_sense, log_fl)
    raw_n_anti = raw_rna_lik(log_p_anti, log_fl)
    g_lik = gdna_lik(log_fl, span, gdna_flank=200, frag_len=frag_len)
    em_ns, em_g, ds = effective_log_lik_pair(
        log_fl_rna=log_fl, log_fl_gdna=log_fl, log_p_strand_rna=log_p_sense,
        span=span, frag_len=frag_len,
    )
    em_na, _, da = effective_log_lik_pair(
        log_fl_rna=log_fl, log_fl_gdna=log_fl, log_p_strand_rna=log_p_anti,
        span=span, frag_len=frag_len,
    )
    print(f"  per-fragment raw: nRNA(sense)={raw_n_sense:+.3f} nRNA(anti)={raw_n_anti:+.3f}  gDNA={g_lik:+.3f}")
    print(f"  effective at EM:  nRNA(sense)={em_ns:+.3f}  nRNA(anti)={em_na:+.3f}  gDNA={em_g:+.3f}")
    print(f"  per-frag effective Δ (nRNA−gDNA): sense={ds:+.3f}  anti={da:+.3f}")

    half = N // 2
    specs = (
        [FragmentSpec(raw_n_sense, g_lik, UNSPLICED_SENSE)] * half
        + [FragmentSpec(raw_n_anti, g_lik, UNSPLICED_ANTI)] * (N - half)
    )

    for label, ag, ar in [
        ("alpha_gdna≫alpha_rna", 0.99 * N, 0.01 * N),
        ("alpha_gdna=alpha_rna", 0.5 * N, 0.5 * N),
        ("alpha_gdna≪alpha_rna", 0.01 * N, 0.99 * N),
    ]:
        res = run_em(specs, alpha_gdna=ag, alpha_rna=ar, em_mode="vbem",
                     span=span, exonic_len=exonic)
        _print_row(label, res, truth_n=N)


def experiment_lik_offset_sweep(N: int = 10_000) -> None:
    """Hold prior 50/50, sweep nRNA−gDNA log-lik offset (post-correction)."""
    _print_header("E3 — EFFECTIVE LIK OFFSET SWEEP (unstranded, 50/50 prior)")
    span = 10_000
    exonic = 1000
    frag_len = 200
    eff_rna = max(span - frag_len + 1, 1)
    g_lik = gdna_lik(-5.0, span, gdna_flank=200, frag_len=frag_len)  # fixed
    # Choose raw_rna so that effective_rna - effective_gdna = delta
    for delta in [-3.0, -1.0, -0.5, 0.0, 0.5, 1.0, 3.0]:
        target_eff_rna = g_lik + delta
        raw_rna = target_eff_rna + math.log(eff_rna)
        specs = [FragmentSpec(raw_rna, g_lik, UNSPLICED_SENSE)] * N
        res = run_em(specs, alpha_gdna=0.5 * N, alpha_rna=0.5 * N,
                     span=span, exonic_len=exonic)
        gdna_frac = res["gdna"] / N
        nrna_frac = res["nrna_t1"] / N
        print(f"  effective Δ={delta:+.2f}nat (nRNA−gDNA)   gdna_frac={gdna_frac:.3f}  "
              f"nrna_frac={nrna_frac:.3f}")


def experiment_prior_sweep(N: int = 10_000) -> None:
    """Hold likelihoods equal (effective), sweep alpha_gdna fraction."""
    _print_header("E4 — PRIOR SWEEP (equal effective likelihoods)")
    span = 10_000
    exonic = 1000
    frag_len = 200
    eff_rna = max(span - frag_len + 1, 1)
    g_lik = gdna_lik(-5.0, span, gdna_flank=200, frag_len=frag_len)
    raw_rna = g_lik + math.log(eff_rna)   # so effective nRNA == gDNA
    specs = [FragmentSpec(raw_rna, g_lik, UNSPLICED_SENSE)] * N
    for ag_frac in [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99]:
        ag = ag_frac * N
        ar = (1.0 - ag_frac) * N
        res = run_em(specs, alpha_gdna=ag, alpha_rna=ar,
                     span=span, exonic_len=exonic)
        gdna_frac = res["gdna"] / N
        print(f"  alpha_gdna_frac={ag_frac:.2f}   gdna_frac={gdna_frac:.3f}  "
              f"nrna_frac={res['nrna_t1'] / N:.3f}")


def experiment_mixed_locus(N_mrna: int = 9000, N_gdna: int = 1000) -> None:
    """The realistic case: locus dominated by exonic mRNA + a small minority
    of intronic gDNA-origin fragments.

    Each "mRNA" fragment has only the mRNA (t0) candidate (no nRNA, no gDNA
    — it's spliced or geometrically incompatible with gDNA).
    Each "gDNA-intronic" fragment has nRNA (t1) + gDNA candidates as in E1.

    Truth: gDNA component should hold N_gdna; nRNA should hold 0.
    """
    _print_header(
        f"E5 — MIXED LOCUS: {N_mrna} exonic-mRNA + {N_gdna} intronic-gDNA "
        "(unstranded, single transcript)"
    )
    span = 10_000
    exonic = 1000
    frag_len = 200
    log_fl = -5.0
    LOG_HALF = math.log(0.5)

    raw_n = raw_rna_lik(LOG_HALF, log_fl)
    g_lik = gdna_lik(log_fl, span, gdna_flank=200, frag_len=frag_len)
    raw_m = raw_rna_lik(LOG_HALF, log_fl)
    em_m_eff = raw_m - math.log(max(exonic - frag_len + 1, 1))
    em_n_eff = raw_n - math.log(max(span - frag_len + 1, 1))
    print(f"  effective at EM: mRNA={em_m_eff:+.3f}  nRNA={em_n_eff:+.3f}  gDNA={g_lik:+.3f}")

    # Build mixed specs by directly constructing ScoredFragments here
    # (the helper above only supports nRNA+gDNA per fragment).
    from rigel.scored_fragments import Locus, ScoredFragments
    from rigel.locus_partition import partition_and_free
    n_total = N_mrna + N_gdna
    offsets = np.zeros(n_total + 1, dtype=np.int64)
    flat_t = np.empty(n_total, dtype=np.int32)
    flat_lk = np.empty(n_total, dtype=np.float64)
    flat_cc = np.empty(n_total, dtype=np.uint8)

    # First N_mrna fragments → mRNA candidate only (t0), assignable
    for i in range(N_mrna):
        flat_t[i] = 0
        flat_lk[i] = raw_m
        flat_cc[i] = UNSPLICED_SENSE
        offsets[i + 1] = i + 1
    # Next N_gdna fragments → nRNA candidate (t1), with gDNA appended by C++
    for j in range(N_gdna):
        i = N_mrna + j
        flat_t[i] = 1
        flat_lk[i] = raw_n
        flat_cc[i] = UNSPLICED_SENSE
        offsets[i + 1] = i + 1

    is_spliced = np.zeros(n_total, dtype=bool)
    gdna_lk_arr = np.full(n_total, -np.inf, dtype=np.float64)
    gdna_lk_arr[N_mrna:] = g_lik   # only intronic units get a finite gDNA lik

    em_data = ScoredFragments(
        offsets=offsets, t_indices=flat_t, log_liks=flat_lk,
        count_cols=flat_cc,
        coverage_weights=np.ones(n_total, dtype=np.float64),
        tx_starts=np.zeros(n_total, dtype=np.int32),
        tx_ends=np.full(n_total, frag_len, dtype=np.int32),
        locus_t_indices=flat_t.copy(),
        locus_count_cols=flat_cc.copy(),
        is_spliced=is_spliced,
        gdna_log_liks=gdna_lk_arr,
        frag_ids=np.arange(n_total, dtype=np.int64),
        frag_class=np.zeros(n_total, dtype=np.int8),
        splice_type=np.zeros(n_total, dtype=np.uint8),
        n_units=n_total,
        n_candidates=n_total,
    )
    locus = Locus(
        locus_id=0,
        transcript_indices=np.array([0, 1], dtype=np.int32),
        unit_indices=np.arange(n_total, dtype=np.int32),
        gdna_span=span,
        merged_intervals=[("chr1", 0, span)],
    )

    index = _MockIndex(span=span, exonic_len=exonic)
    rc = AbundanceEstimator(
        num_transcripts=2,
        em_config=EMConfig(
            seed=42, mode="vbem", iterations=2000,
            assignment_mode="fractional", n_threads=1,
        ),
        is_nrna=np.array([False, True]),
        is_synthetic=np.array([False, True]),
    )
    rc._transcript_spans = np.array([span, span], dtype=np.float64)
    rc._exonic_lengths = np.array([exonic, span], dtype=np.float64)

    partitions = partition_and_free(em_data, [locus])
    p = partitions[0]
    pt = (p.offsets, p.t_indices, p.log_liks, p.coverage_weights,
          p.tx_starts, p.tx_ends, p.count_cols, p.is_spliced,
          p.gdna_log_liks, p.locus_t_indices, p.locus_count_cols)

    # Sweep alpha_gdna fraction WITHIN this locus (truth ratio = N_gdna/n_total)
    truth_gdna_frac = N_gdna / n_total
    print(f"  truth_gdna_frac = {truth_gdna_frac:.3f}")
    c_base = float(n_total)
    for ag_frac in [0.001, 0.01, truth_gdna_frac, 0.1, 0.3, 0.5]:
        ag = c_base * ag_frac
        ar = c_base * (1.0 - ag_frac)
        rc.em_counts.fill(0.0)
        rc._gdna_em_total = 0.0
        total_gdna, _, _ = rc.run_batch_locus_em_partitioned(
            [pt], [np.array([0, 1], dtype=np.int32)],
            np.array([ag], dtype=np.float64),
            np.array([ar], dtype=np.float64),
            index,
            em_iterations=2000,
        )
        em_per_t = rc.em_counts.sum(axis=1)
        nrna = float(em_per_t[1])
        mrna = float(em_per_t[0])
        gdna = float(total_gdna)
        siphon = nrna  # nRNA truth = 0
        print(f"  alpha_gdna_frac={ag_frac:.4f}  mRNA={mrna:>7.1f}  "
              f"nRNA={nrna:>6.1f}  gDNA={gdna:>6.1f}  "
              f"siphon={siphon:>5.1f}/{N_gdna}  "
              f"({100*siphon/max(N_gdna,1):.1f}%)")


def main() -> None:
    experiment_unstranded_pure_gdna()
    experiment_stranded_pure_gdna(ss=0.99)
    experiment_stranded_pure_gdna(ss=0.50)
    experiment_lik_offset_sweep()
    experiment_prior_sweep()
    experiment_mixed_locus(N_mrna=9000, N_gdna=1000)
    experiment_mixed_locus(N_mrna=99000, N_gdna=1000)
    experiment_mixed_locus(N_mrna=999000, N_gdna=1000)
    experiment_mixed_locus_stranded(ss=0.99, N_mrna=9000, N_gdna=1000)
    experiment_mixed_locus_stranded(ss=0.99, N_mrna=99000, N_gdna=1000)
    experiment_proposed_fix()


def experiment_proposed_fix() -> None:
    """Compare current locus-prior derivation vs the PROPOSED FIX
    (restrict prior derivation to gDNA-eligible units; scale c_base to
    eligible count) on a battery of mixed-locus scenarios.
    """
    _print_header("E7 — PROPOSED FIX vs CURRENT (per-locus prior derivation)")
    span = 10_000
    exonic = 1000
    frag_len = 200
    log_fl = -5.0
    LOG_HALF = math.log(0.5)
    raw_n = raw_rna_lik(LOG_HALF, log_fl)
    raw_m = raw_rna_lik(LOG_HALF, log_fl)
    g_lik = gdna_lik(log_fl, span, gdna_flank=200, frag_len=frag_len)

    # Compute the per-fragment gamma values via the same formula
    # as compute_locus_priors_from_partitions (z = log_prior + gdna_ll - max_rna)
    # using the EFFECTIVE rna lik (post bias correction).
    em_rna = raw_n - math.log(span - frag_len + 1)

    scenarios = [
        # (label, N_mrna, N_gdna, pi_pool)
        ("100% intronic-pool gDNA (truth pi_pool=1.0)", 9000, 1000, 0.999),
        ("100% intronic-pool gDNA (truth pi_pool=1.0)", 99000, 1000, 0.999),
        ("50% intronic-pool gDNA (mixed truth)", 9000, 1000, 0.50),
        ("calibrated underestimate (pi_pool=0.45 vs truth 1.0)", 9000, 1000, 0.45),
    ]

    from rigel.locus_partition import partition_and_free
    for label, N_mrna, N_gdna, pi_pool in scenarios:
        n_total = N_mrna + N_gdna
        offsets = np.zeros(n_total + 1, dtype=np.int64)
        flat_t = np.empty(n_total, dtype=np.int32)
        flat_lk = np.empty(n_total, dtype=np.float64)
        flat_cc = np.empty(n_total, dtype=np.uint8)
        for i in range(N_mrna):
            flat_t[i] = 0; flat_lk[i] = raw_m; flat_cc[i] = UNSPLICED_SENSE
            offsets[i + 1] = i + 1
        for j in range(N_gdna):
            i = N_mrna + j
            flat_t[i] = 1; flat_lk[i] = raw_n; flat_cc[i] = UNSPLICED_SENSE
            offsets[i + 1] = i + 1
        is_spliced = np.zeros(n_total, dtype=bool)
        gdna_lk_arr = np.full(n_total, -np.inf, dtype=np.float64)
        gdna_lk_arr[N_mrna:] = g_lik

        locus = Locus(
            locus_id=0, transcript_indices=np.array([0, 1], dtype=np.int32),
            unit_indices=np.arange(n_total, dtype=np.int32),
            gdna_span=span, merged_intervals=[("chr1", 0, span)],
        )

        # Compute priors two ways:
        log_prior = math.log(pi_pool / (1.0 - pi_pool))
        max_rna_per_unit = np.empty(n_total)
        max_rna_per_unit[:N_mrna] = (raw_m - math.log(exonic - frag_len + 1))
        max_rna_per_unit[N_mrna:] = em_rna
        z = log_prior + gdna_lk_arr - max_rna_per_unit
        z = np.where(np.isfinite(z), np.clip(z, -50, 50), -50)
        gamma_f = 1.0 / (1.0 + np.exp(-z))

        c_base_current = 5.0
        ag_cur = c_base_current * float(gamma_f.mean())
        ar_cur = c_base_current * (1.0 - float(gamma_f.mean()))

        eligible = np.isfinite(gdna_lk_arr)
        gamma_eligible = gamma_f[eligible]
        # FIXED: eligible-only mean, c_base=5.0 (matches new locus.py)
        gamma_l_e = float(gamma_eligible.mean())
        ag_fix = c_base_current * gamma_l_e
        ar_fix = c_base_current * (1.0 - gamma_l_e)
        # FIXED-STRONG: c_base=10 (current production default)
        ag_strong = 10.0 * gamma_l_e
        ar_strong = 10.0 * (1.0 - gamma_l_e)

        index = _MockIndex(span=span, exonic_len=exonic)

        def _build_em_data():
            return ScoredFragments(
                offsets=offsets.copy(), t_indices=flat_t.copy(),
                log_liks=flat_lk.copy(),
                count_cols=flat_cc.copy(),
                coverage_weights=np.ones(n_total, dtype=np.float64),
                tx_starts=np.zeros(n_total, dtype=np.int32),
                tx_ends=np.full(n_total, frag_len, dtype=np.int32),
                locus_t_indices=flat_t.copy(), locus_count_cols=flat_cc.copy(),
                is_spliced=is_spliced.copy(),
                gdna_log_liks=gdna_lk_arr.copy(),
                frag_ids=np.arange(n_total, dtype=np.int64),
                frag_class=np.zeros(n_total, dtype=np.int8),
                splice_type=np.zeros(n_total, dtype=np.uint8),
                n_units=n_total, n_candidates=n_total,
            )

        def _run(ag, ar):
            rc = AbundanceEstimator(
                num_transcripts=2,
                em_config=EMConfig(seed=42, mode="vbem", iterations=2000,
                                   assignment_mode="fractional", n_threads=1),
                is_nrna=np.array([False, True]),
                is_synthetic=np.array([False, True]),
            )
            rc._transcript_spans = np.array([span, span], dtype=np.float64)
            rc._exonic_lengths = np.array([exonic, span], dtype=np.float64)
            partitions = partition_and_free(_build_em_data(), [locus])
            p = partitions[0]
            pt = (p.offsets, p.t_indices, p.log_liks, p.coverage_weights,
                  p.tx_starts, p.tx_ends, p.count_cols, p.is_spliced,
                  p.gdna_log_liks, p.locus_t_indices, p.locus_count_cols)
            total_gdna, _, _ = rc.run_batch_locus_em_partitioned(
                [pt], [np.array([0, 1], dtype=np.int32)],
                np.array([ag], dtype=np.float64),
                np.array([ar], dtype=np.float64),
                index, em_iterations=2000,
            )
            em_per_t = rc.em_counts.sum(axis=1)
            return float(em_per_t[0]), float(em_per_t[1]), float(total_gdna)

        m_cur, n_cur, g_cur = _run(ag_cur, ar_cur)
        m_fix, n_fix, g_fix = _run(ag_fix, ar_fix)
        m_str, n_str, g_str = _run(ag_strong, ar_strong)

        print(f"\n  scenario: {label}  (N_mRNA={N_mrna}, N_gDNA={N_gdna}, pi_pool={pi_pool})")
        print(f"    CURRENT  alpha=({ag_cur:.3f},{ar_cur:.3f})  →  "
              f"mRNA={m_cur:.0f} nRNA={n_cur:.0f} gDNA={g_cur:.0f}  "
              f"siphon={n_cur:.0f}/{N_gdna} ({100*n_cur/max(N_gdna,1):.1f}%)")
        print(f"    FIXED    alpha=({ag_fix:.3f},{ar_fix:.3f})  →  "
              f"mRNA={m_fix:.0f} nRNA={n_fix:.0f} gDNA={g_fix:.0f}  "
              f"siphon={n_fix:.0f}/{N_gdna} ({100*n_fix/max(N_gdna,1):.1f}%)")
        print(f"    STRONG   alpha=({ag_strong:.3f},{ar_strong:.3f})  →  "
              f"mRNA={m_str:.0f} nRNA={n_str:.0f} gDNA={g_str:.0f}  "
              f"siphon={n_str:.0f}/{N_gdna} ({100*n_str/max(N_gdna,1):.1f}%)")


def experiment_mixed_locus_stranded(ss: float, N_mrna: int, N_gdna: int) -> None:
    """Same as E5 but stranded — half sense, half antisense gDNA fragments."""
    _print_header(
        f"E6 — MIXED LOCUS STRANDED ss={ss}: {N_mrna} exonic-mRNA + "
        f"{N_gdna} intronic-gDNA (50/50 sense:anti)"
    )
    span = 10_000
    exonic = 1000
    frag_len = 200
    log_fl = -5.0
    log_p_sense = math.log(ss)
    log_p_anti = math.log(1.0 - ss + 1e-300)

    raw_m = raw_rna_lik(log_p_sense, log_fl)   # mRNA assumed sense
    raw_n_sense = raw_rna_lik(log_p_sense, log_fl)
    raw_n_anti = raw_rna_lik(log_p_anti, log_fl)
    g_lik = gdna_lik(log_fl, span, gdna_flank=200, frag_len=frag_len)
    em_ns = raw_n_sense - math.log(span - frag_len + 1)
    em_na = raw_n_anti - math.log(span - frag_len + 1)
    print(f"  effective intronic Δ (nRNA−gDNA): "
          f"sense={em_ns - g_lik:+.3f}  anti={em_na - g_lik:+.3f}")

    from rigel.locus_partition import partition_and_free
    n_total = N_mrna + N_gdna
    half_g = N_gdna // 2
    offsets = np.zeros(n_total + 1, dtype=np.int64)
    flat_t = np.empty(n_total, dtype=np.int32)
    flat_lk = np.empty(n_total, dtype=np.float64)
    flat_cc = np.empty(n_total, dtype=np.uint8)
    for i in range(N_mrna):
        flat_t[i] = 0
        flat_lk[i] = raw_m
        flat_cc[i] = UNSPLICED_SENSE
        offsets[i + 1] = i + 1
    for j in range(N_gdna):
        i = N_mrna + j
        flat_t[i] = 1
        if j < half_g:
            flat_lk[i] = raw_n_sense
            flat_cc[i] = UNSPLICED_SENSE
        else:
            flat_lk[i] = raw_n_anti
            flat_cc[i] = UNSPLICED_ANTI
        offsets[i + 1] = i + 1

    is_spliced = np.zeros(n_total, dtype=bool)
    gdna_lk_arr = np.full(n_total, -np.inf, dtype=np.float64)
    gdna_lk_arr[N_mrna:] = g_lik

    em_data = ScoredFragments(
        offsets=offsets, t_indices=flat_t, log_liks=flat_lk,
        count_cols=flat_cc,
        coverage_weights=np.ones(n_total, dtype=np.float64),
        tx_starts=np.zeros(n_total, dtype=np.int32),
        tx_ends=np.full(n_total, frag_len, dtype=np.int32),
        locus_t_indices=flat_t.copy(),
        locus_count_cols=flat_cc.copy(),
        is_spliced=is_spliced,
        gdna_log_liks=gdna_lk_arr,
        frag_ids=np.arange(n_total, dtype=np.int64),
        frag_class=np.zeros(n_total, dtype=np.int8),
        splice_type=np.zeros(n_total, dtype=np.uint8),
        n_units=n_total,
        n_candidates=n_total,
    )
    locus = Locus(
        locus_id=0, transcript_indices=np.array([0, 1], dtype=np.int32),
        unit_indices=np.arange(n_total, dtype=np.int32),
        gdna_span=span, merged_intervals=[("chr1", 0, span)],
    )
    index = _MockIndex(span=span, exonic_len=exonic)
    rc = AbundanceEstimator(
        num_transcripts=2,
        em_config=EMConfig(seed=42, mode="vbem", iterations=2000,
                           assignment_mode="fractional", n_threads=1),
        is_nrna=np.array([False, True]),
        is_synthetic=np.array([False, True]),
    )
    rc._transcript_spans = np.array([span, span], dtype=np.float64)
    rc._exonic_lengths = np.array([exonic, span], dtype=np.float64)
    partitions = partition_and_free(em_data, [locus])
    p = partitions[0]
    pt = (p.offsets, p.t_indices, p.log_liks, p.coverage_weights,
          p.tx_starts, p.tx_ends, p.count_cols, p.is_spliced,
          p.gdna_log_liks, p.locus_t_indices, p.locus_count_cols)

    # Two prior regimes:
    #   (a) "real-data" — alpha_gdna derived from WHOLE locus (mostly mRNA, dilutes gDNA)
    #   (b) "intronic-pool" — alpha_gdna scaled to the intronic-only fraction
    c_base_whole = float(n_total)
    c_base_intronic = float(N_gdna)
    truth_intronic_gdna_frac = 1.0  # intronic pool is 100% gDNA in this scenario
    print(f"  truth: locus has {N_mrna} mRNA + {N_gdna} gDNA; intronic pool is 100% gDNA")
    print("  --- regime (a): alpha_gdna derived from WHOLE-LOCUS share (current behaviour) ---")
    for ag_frac in [0.001, N_gdna / n_total, 0.1, 0.5]:
        ag = c_base_whole * ag_frac
        ar = c_base_whole * (1.0 - ag_frac)
        rc.em_counts.fill(0.0); rc._gdna_em_total = 0.0
        total_gdna, _, _ = rc.run_batch_locus_em_partitioned(
            [pt], [np.array([0, 1], dtype=np.int32)],
            np.array([ag], dtype=np.float64),
            np.array([ar], dtype=np.float64),
            index, em_iterations=2000,
        )
        em_per_t = rc.em_counts.sum(axis=1)
        nrna = float(em_per_t[1]); mrna = float(em_per_t[0]); gdna = float(total_gdna)
        print(f"    alpha_gdna_frac={ag_frac:.4f}  mRNA={mrna:>7.1f}  "
              f"nRNA={nrna:>6.1f}  gDNA={gdna:>6.1f}  "
              f"siphon={nrna:>5.1f}/{N_gdna} ({100*nrna/max(N_gdna,1):.1f}%)")
    print("  --- regime (b): alpha_gdna scaled to INTRONIC POOL only ---")
    ag = c_base_intronic * 1.0
    ar = c_base_intronic * 0.001
    rc.em_counts.fill(0.0); rc._gdna_em_total = 0.0
    total_gdna, _, _ = rc.run_batch_locus_em_partitioned(
        [pt], [np.array([0, 1], dtype=np.int32)],
        np.array([ag], dtype=np.float64),
        np.array([ar], dtype=np.float64),
        index, em_iterations=2000,
    )
    em_per_t = rc.em_counts.sum(axis=1)
    nrna = float(em_per_t[1]); mrna = float(em_per_t[0]); gdna = float(total_gdna)
    print(f"    alpha_gdna={ag:.0f} (intronic-only c_base)   mRNA={mrna:>7.1f}  "
          f"nRNA={nrna:>6.1f}  gDNA={gdna:>6.1f}  "
          f"siphon={nrna:>5.1f}/{N_gdna} ({100*nrna/max(N_gdna,1):.1f}%)")


if __name__ == "__main__":
    main()
