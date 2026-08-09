#!/usr/bin/env python
"""⭐⭐⭐ **CALIBRATION'S ENDPOINT AGAINST TRUTH** — ``LocusPriors`` vs the origin-split oracle.

**The question, and why nothing answered it before.** Calibration does not ship a number; it ships a
PRIOR. Three float64 arrays indexed by ``multi_locus_id`` (``calibration/priors.py``)::

    gdna_prior_count   the gDNA component's Dirichlet pseudocount
    rna_prior_count    the RNA group's pseudocount (the EM splits it by evidence)
    gdna_eff_len       the IPR-contracted effective length of the gDNA component

The first two are **FRAGMENT COUNTS** — ``em_solver.cpp:apply_grouped_prior_update`` adds them straight
to the EM's own soft counts (``G = n_gdna + a_g``, ``R = n_rna + a_r``) — and that is what makes an
oracle version well-defined rather than a modelling opinion. Every other instrument in this repo scores
an INTERMEDIATE: ``solvability_audit.py`` and ``pass0_vs_oracle.py`` score per-object ``f_g``,
``calibration_truth_ab.py`` scores the library figure. Nothing scored the thing the EM actually reads.

⭐ **THREE ARMS, AND THE TWO DIFFERENCES BETWEEN THEM ARE THE WHOLE DIAGNOSTIC.**

===  ==========================================================  ===================================
 P   the SHIPPED prior                                           ``assemble_priors(cal, ra, loci)``
 O   the prior a PERFECT DECONVOLUTION would produce             truth masses, same assembler
 F   the DIRECT per-locus fragment truth                         ``node_start_count``, per origin
===  ==========================================================  ===================================

``P − O`` is calibration's own error — everything upstream of :func:`assemble_priors`.
``O − F`` is the ASSEMBLER's error — the mass → density → fragment-count conversion, the
overlap projection, and the support-weighted pooling ``assemble_priors``' own docstring flags as an
approximation. They are different repairs in different files and pooling them would hide which.

⛔ **O IS NOT AN ESTIMATOR AND MUST NEVER BECOME ONE** (the shape ``pass0_vs_oracle.py`` refuses too).
It is the shipped :func:`assemble_priors` fed the ONE lever that already exists —
``OracleTruth.override_masses``, the true per-object gDNA/RNA split — with every other input untouched.
Writing a "best possible prior" here would be the magic-number failure with an estimator in place of a
constant.

⭐ **F IS EXACT ON THE gDNA ARM AND A BOUND ON THE RNA ARM, and that asymmetry is physics.** gDNA does
not splice, so the gdna partition's ``node_start_count`` — one deposit per accepted fragment, at the
node holding its first base — projected through the SHIPPED
``priors._project_regions_to_loci`` IS the true number of gDNA fragments starting in the locus, with
nothing to subtract. The RNA arm's target is the UNSPLICED count (``assemble_priors`` withholds spliced
mass because certified RNA has no gDNA competitor in the EM), and ``node_start_count`` carries no
splice bit — so ``F_rna`` here is spliced-INCLUSIVE and is reported as an **upper bound**, labelled as
one. Separating it needs a five-way origin×splice BAM split and a second oracle cache; it is not done.

⛔⛔ **UNDRAINED ON EVERY ARM, AND THAT IS FORCED — WITH A MEASUREMENT, NOT AN ASSUMPTION.**
Production drains the side buffer before calibration. ``lift_choices`` exists to make a DRAINED oracle
valid, and on this panel it is **not**: measured at ``g25 ss0.50 capture_on``, 6,450 of 209,658 held
fragments (3.08 %) tie on the deferred bank's canonical key across origins, and replaying the whole's
choices then deposits 1,026 ``edge_spliced`` + 526 ``sj`` records inside the **gdna** partition —
which ``OracleTruth._validate`` refuses, correctly, because gDNA cannot splice. The ladder is built
with EQUAL configured fragment lengths, so span collisions are the common case there by construction.
⭐ So both sides run undrained, and the caveat is PRICED rather than waved at: the ``drained`` arm
recomputes the shipped prior on the production (drained) payload, so ``P_drained − P`` is the
drain's own contribution as a number. Pass ``--no-drained-arm`` to skip it.

⛔ **THE UNITS THE ERROR IS REPORTED IN, and both are needed.**

* **fragments** — ``Σ|ΔA|``, additive, comparable across strata, and the unit the EM adds in.
* **the composition claim** ``phi = a_g / (a_g + a_r)`` — mass-weighted ``|Δphi|``, the analogue of
  every other instrument's ``mwae``. A prior can carry the right RATIO at the wrong SCALE (it then
  pulls the EM the right way but too weakly) or the right scale at the wrong ratio (it pulls hard in
  the wrong direction), and one number cannot distinguish those. ``scale = a_g + a_r`` is printed too.

Gates: ``tests/calibration/test_prior_vs_oracle.py``, each carrying its own perturbation.

Usage::

    python scripts/design/prior_vs_oracle.py --suite ~/Downloads/rigel_runs/suite/ladder \\
        --oracle-cache ~/Downloads/rigel_runs/suite/ladder/oracle_cache --jobs 6
    python scripts/design/prior_vs_oracle.py --conditions gdna_g25_ss_0.50_nrna_none_capture_on
    python scripts/design/prior_vs_oracle.py --json out.json --emit-oracle-masses DIR
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "tests" / "calibration"))

from _oracle import ORIGINS, OracleTruth, _split_bam  # noqa: E402

import rigel.calibration.priors as PRIORS  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import (  # noqa: E402
    _drain_side_buffer,
    _native_detect_sj_tag,
    quant_from_buffer,
    scan_and_buffer,
)
from rigel.scan_cache import ScanCacheKeyError, read_scan_cache, write_scan_cache  # noqa: E402

_RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_SUITE = _RUNS / "suite" / "ladder"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"

#: The three ``LocusPriors`` fields, in the order every table prints them.
PRIOR_FIELDS = ("gdna_prior_count", "rna_prior_count", "gdna_eff_len")

#: The mass arrays ``OracleTruth.override_masses`` replaces. Named here so the ``noop`` gate can
#: re-inject exactly this set from the SHIPPED result and demand byte-identity — an override applied
#: to a field nothing reads is an override that never ran (TRAPS: an-ablation-that-never-ran).
OVERRIDE_FIELDS = (
    "mass_gdna_node",
    "mass_rna_node",
    "mass_gdna_edge",
    "mass_rna_edge",
    "mass_rna_spliced_edge",
    "mass_rna_junction",
)


# ── the scoring ──────────────────────────────────────────────────────────────────────────────────


def composition(gdna: np.ndarray, rna: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``(phi, scale)`` per locus — ``phi = a_g/(a_g+a_r)`` is **NaN** where the prior is empty.

    ⛔ NaN, never 0, and for the same reason ``pass0_vs_oracle.object_fractions`` gives: a locus whose
    prior is ``(0, 0)`` makes no composition claim at all, and a floored 0 reads as the confident claim
    "no gDNA here". The mass-weighted mean is blind to the difference — a zero-scale locus carries zero
    weight — so it shows up only in the COUNT of loci scored, which is exactly where an inflated
    denominator is invisible.
    """
    g = np.asarray(gdna, np.float64)
    scale = g + np.asarray(rna, np.float64)
    phi = np.full(scale.shape, np.nan, dtype=np.float64)
    np.divide(g, scale, out=phi, where=scale > 0.0)
    return phi, scale


@dataclass(frozen=True, slots=True)
class ArmScore:
    """One arm scored against one reference, on ONE prior arm (gDNA or RNA), over one locus selection.

    Every field is a fragment total except ``rel`` and ``mwae_phi`` — so the totals ADD across a
    partition of the loci and the two rates do not.
    """

    n_loci: int  #: loci in the selection. ⚠ never the panel's locus count
    n_claiming: int  #: loci where EITHER arm or reference puts a nonzero count on this arm
    total_arm: float  #: Σ a, the arm's own fragment total
    total_ref: float  #: Σ a*, the reference's
    net_err: float  #: Σ (a − a*). What a library-level figure would see
    abs_err: float  #: Σ |a − a*|. What the per-locus answer is
    over_call: float  #: Σ (a − a*)+
    under_call: float  #: Σ (a* − a)+

    @property
    def rel(self) -> float:
        """``Σ|a − a*| / Σ a*`` — the absolute error as a fraction of the reference total."""
        return self.abs_err / self.total_ref if self.total_ref > 0 else float("nan")

    @property
    def cancellation(self) -> float:
        """``Σ|err| / |net|``. ⭐ Large means a big under-call is sitting next to a big over-call and
        any library-level summary of this arm is flattering it."""
        return self.abs_err / abs(self.net_err) if self.net_err != 0.0 else float("inf")


def score_arm(arm: np.ndarray, ref: np.ndarray, select: np.ndarray | None = None) -> ArmScore:
    """Score one prior arm (a per-locus fragment count) against a reference arm."""
    a = np.asarray(arm, np.float64)
    r = np.asarray(ref, np.float64)
    if a.shape != r.shape:
        raise ValueError(
            f"arm has {a.shape[0]} loci, reference has {r.shape[0]} — these are different locus "
            "partitions. A prior scored against another run's loci is not a small error."
        )
    if select is not None:
        a, r = a[select], r[select]
    d = a - r
    return ArmScore(
        n_loci=int(a.shape[0]),
        n_claiming=int(np.count_nonzero((a > 0.0) | (r > 0.0))),
        total_arm=float(a.sum()),
        total_ref=float(r.sum()),
        net_err=float(d.sum()),
        abs_err=float(np.abs(d).sum()),
        over_call=float(np.maximum(d, 0.0).sum()),
        under_call=float(np.maximum(-d, 0.0).sum()),
    )


def score_composition(arm, ref, select=None) -> dict:
    """The ``phi`` / ``scale`` view of one arm-vs-reference pair — what the EM actually feels.

    ⭐ **Weighted by the REFERENCE's own scale**, never by the arm's. Weighting by the arm lets a
    mechanism improve the number by shrinking its own prior to nothing at the loci it gets wrong, which
    is TRAPS: honesty-metrics-reward-ignorance with the weight in place of the denominator. The reference's scale is a
    property of the condition, so it is fixed across every arm compared here.
    """
    phi_a, sc_a = composition(arm.gdna_prior_count, arm.rna_prior_count)
    phi_r, sc_r = composition(ref.gdna_prior_count, ref.rna_prior_count)
    live = np.isfinite(phi_a) & np.isfinite(phi_r)
    if select is not None:
        live = live & select
    w = sc_r[live]
    d_phi = np.abs(phi_a[live] - phi_r[live])
    tot_w = float(w.sum())
    # scale: a RATIO of totals, reported in log10 so a 10x under-prior and a 10x over-prior are
    # symmetric. ⚠ Only over the loci both arms make a claim at, matching mwae_phi's population.
    sa, sr = float(sc_a[live].sum()), float(sc_r[live].sum())
    return {
        "n_scored": int(live.sum()),
        "weight": tot_w,
        "mwae_phi": float((w * d_phi).sum() / tot_w) if tot_w > 0 else float("nan"),
        "median_phi_err": float(np.median(d_phi)) if d_phi.size else float("nan"),
        "scale_arm": sa,
        "scale_ref": sr,
        "scale_log10_ratio": float(np.log10(sa / sr)) if sa > 0 and sr > 0 else float("nan"),
    }


def score_eff_len(arm, ref, select=None) -> dict:
    """``gdna_eff_len``, weighted by the reference's gDNA prior count.

    ⚠ Weighted by the **gDNA** count and not by the total: ``gdna_eff_len`` divides the gDNA
    component's abundance alone, so a locus with no gDNA is a locus where this array does nothing, and
    including it at full weight would report the error of a number nothing reads
    (TRAPS: weight-it-like-the-consumer).
    """
    a = np.asarray(arm.gdna_eff_len, np.float64)
    r = np.asarray(ref.gdna_eff_len, np.float64)
    w = np.asarray(ref.gdna_prior_count, np.float64)
    live = np.isfinite(a) & np.isfinite(r) & (r > 0.0)
    if select is not None:
        live = live & select
    w = np.where(live, w, 0.0)
    tot = float(w.sum())
    rel = np.zeros_like(a)
    np.divide(np.abs(a - r), r, out=rel, where=live)
    return {
        "n_scored": int(live.sum()),
        "w_rel_err": float((w * rel).sum() / tot) if tot > 0 else float("nan"),
        "median_rel_err": float(np.median(rel[live])) if live.any() else float("nan"),
        "total_arm": float(a[live].sum()),
        "total_ref": float(r[live].sum()),
    }


# ── the arms ─────────────────────────────────────────────────────────────────────────────────────


def capture_priors(buffer, index, strand_models, fl, region_arrays, stats, calibration, payload,
                   pipeline_config):
    """Run the production quant path far enough to get ``(multi_loci, LocusPriors)``, then STOP.

    ⭐ **The loci are the production loci, not a re-derivation.** ``build_multi_loci`` builds
    connected components of transcripts linked by SCORED fragments, so the locus partition is a
    function of the scoring stage and cannot be reconstructed from the index alone. Wrapping
    :func:`~rigel.calibration.priors.assemble_priors` — which ``quant_from_buffer`` imports
    function-locally, so patching the module attribute is picked up at call time — takes both objects
    from the call production itself makes (TRAPS: a-test-that-redefines).

    ⛔ The sentinel exception is what makes this affordable: the per-locus EM is the single most
    expensive stage (47 % of a condition) and this instrument does not read its output. ⚠ Item 3 —
    injecting the oracle prior and re-quantifying — needs the EM and must NOT reuse this path.
    """

    class _StopAfterPriors(Exception):
        pass

    captured: dict = {}
    original = PRIORS.assemble_priors

    def _wrapper(cal, ra, multi_loci):
        captured["multi_loci"] = multi_loci
        captured["priors"] = original(cal, ra, multi_loci)
        raise _StopAfterPriors

    PRIORS.assemble_priors = _wrapper
    try:
        quant_from_buffer(
            buffer, index, strand_models, fl, region_arrays, stats, calibration, payload,
            em_config=pipeline_config.em, scoring=pipeline_config.scoring,
        )
    except _StopAfterPriors:
        pass
    finally:
        PRIORS.assemble_priors = original
    if "priors" not in captured:
        # ⛔ TRAPS: an-ablation-that-never-ran. ``quant_from_buffer`` returns early on an empty unit set, and a
        # silently-absent capture would read here as a condition with no loci rather than as a
        # harness that never fired.
        raise RuntimeError(
            "assemble_priors was never called — quant_from_buffer returned before the prior stage "
            "(no EM units, or no multi-loci). This is not a condition with zero error."
        )
    return captured["multi_loci"], captured["priors"]


def oracle_priors(oracle: OracleTruth, calibration, region_arrays, multi_loci):
    """**O** — the prior a perfect deconvolution would produce, and the ``noop`` gate beside it.

    ⭐ The ONE lever: ``override_masses`` swaps the six mass arrays for the origin-split truth and
    changes nothing else. ⛔ The ``noop`` arm re-injects the SHIPPED masses through the identical
    ``dataclasses.replace`` and must come back byte-identical (TRAPS: byte-identity-gate) — an override that silently
    landed on the wrong field, or a ``replace`` that dropped a field, would otherwise look like a
    result. Returns ``(O, noop)``.
    """
    override = oracle.override_masses(region_arrays)
    missing = set(OVERRIDE_FIELDS) - set(override)
    extra = set(override) - set(OVERRIDE_FIELDS)
    if missing or extra:
        raise RuntimeError(
            f"override_masses no longer writes OVERRIDE_FIELDS: missing={sorted(missing)} "
            f"extra={sorted(extra)}. The noop gate would be testing a different set than the arm."
        )
    o = PRIORS.assemble_priors(
        dataclasses.replace(calibration, **override), region_arrays, multi_loci
    )
    noop = PRIORS.assemble_priors(
        dataclasses.replace(calibration, **{f: getattr(calibration, f) for f in OVERRIDE_FIELDS}),
        region_arrays,
        multi_loci,
    )
    return o, noop


def share_priors(oracle: OracleTruth, calibration, region_arrays, multi_loci):
    """**S** — the O arm, plus each component rescaled by its OWN true per-line share.

    ⭐⭐ **WHY THIS ARM EXISTS.** ``assemble_priors`` rescales BOTH components at a line by ONE pooled
    share, ``mass / count`` off the mixture. That is exact when the two components share a length
    distribution and biased when they do not — by exactly ``share_r / share_g``, independent of the true
    mixing ratio. ⛔ The bias is **purely compositional**: the locus total is conserved to the last
    fragment, so no conservation gate can see it. Only a per-component comparison can.

    ``O − S`` is therefore the pooled share's own contribution, isolated, and ``S − F`` is everything
    else the assembler does wrong. Until this arm existed the two were summed inside ``O − F`` and there
    was no way to tell which was which.

    ⚠ **The assembler takes ONE share, so the arm is built by calling it TWICE** — once with the gDNA
    truth share (keeping its gDNA arm) and once with the RNA truth share (keeping its RNA arm). That is
    not a re-implementation: it is the shipped function, run twice with one input varied.

    ⛔ **The shares are MEASURED off the origin split** (``OracleTruth.component_shares``), never derived
    from a pmf — see that method for why an analytic share would make this a model arm.
    """
    override = oracle.override_masses(region_arrays)
    shares = oracle.component_shares()
    truth_cal = dataclasses.replace(calibration, **override)
    gdna = PRIORS.assemble_priors(
        dataclasses.replace(truth_cal, edge_mass_per_crossing=shares["gdna"]),
        region_arrays, multi_loci,
    )
    rna = PRIORS.assemble_priors(
        dataclasses.replace(truth_cal, edge_mass_per_crossing=shares["rna"]),
        region_arrays, multi_loci,
    )
    return PRIORS.LocusPriors(
        gdna_prior_count=gdna.gdna_prior_count,
        rna_prior_count=rna.rna_prior_count,
        gdna_eff_len=gdna.gdna_eff_len,
    ), shares


def eff_len_inflation(calibration, region_arrays, multi_loci) -> dict:
    """⭐ Is ``gdna_eff_len`` clamped by an INCIDENCE-support sum rather than the genomic span?

    ``assemble_priors`` clamps ``gdna_eff_len`` to ``span = Σ share·(S_node + S_edge)``. ``S_edge`` is
    ``E_g[w − 1] ≈ mu_g − 1`` PER LINE, so every interior line adds most of a fragment length to a
    locus whose nodes may be a few hundred bases — an incidence-like sum, not a genomic extent. The EM
    divides the gDNA component's abundance by this array, so an inflation here is a direct scale error
    on one of the three numbers calibration ships.

    ⚠ Reports the ratio to the locus's GENOMIC span, mass-weighted by the gDNA prior, so the number is
    what the consumer feels rather than what an unweighted locus average would say
    (``TRAPS: weight-it-like-the-consumer``).
    """
    # ⭐ Nodes and edges are projected on their OWN axes, exactly as `assemble_priors` does — no line is
    # folded onto a flank node. Re-deriving the fold here would measure a span the assembler no longer
    # builds (`TRAPS: a-test-that-redefines`).
    n_loci = len(multi_loci)
    node_support = np.maximum(np.asarray(calibration.gdna_node_eff_len, np.float64), 0.0)
    edge_support = np.maximum(np.asarray(calibration.gdna_edge_eff_len, np.float64), 0.0)
    proj = PRIORS._project_regions_to_loci(
        region_arrays, multi_loci, n_loci,
        {
            "node_only": node_support,
            "genomic": np.asarray(region_arrays.region_size_bp, np.float64),
        },
    )
    e_idx, e_lid, e_w = PRIORS._edge_locus_shares(region_arrays, multi_loci, n_loci)
    proj["support"] = proj["node_only"] + PRIORS._sum_by_locus(
        e_idx, e_lid, e_w, edge_support, n_loci
    )
    live = proj["genomic"] > 0
    ratio = np.divide(proj["support"], proj["genomic"], out=np.ones_like(proj["support"]), where=live)
    node_ratio = np.divide(proj["node_only"], proj["genomic"],
                           out=np.ones_like(proj["support"]), where=live)
    return {
        "median_support_over_genomic": float(np.median(ratio[live])) if live.any() else float("nan"),
        "median_node_over_genomic": float(np.median(node_ratio[live])) if live.any() else float("nan"),
        "total_support": float(proj["support"].sum()),
        "total_genomic": float(proj["genomic"].sum()),
    }


def fragment_truth(oracle: OracleTruth, region_arrays, multi_loci):
    """**F** — the DIRECT per-locus true fragment count, from ``node_start_count`` per origin.

    ⭐ **Why this bank and no other.** ``node_start_count`` is the accumulator's one real invariant —
    ``Σ node_start_count == qc.deposited``, one increment per accepted fragment at the node holding its
    first base. Every other bank counts OBJECT INCIDENCES (a fragment deposits on ``max(K, 1)``
    objects), which is exactly the unit error ``assemble_priors``' own docstring exists to correct. So
    this is the only bank on which "the locus's true fragment count" is a count rather than a model.

    ⭐ **And it is the same quantity ``assemble_priors`` targets, by construction.** That function
    computes ``rho_c · span_bp`` with ``span_bp`` the locus's GENOMIC span; under a component of
    density ``rho`` fragment-starts-per-bp, the expected number of starts inside the span is exactly
    ``rho · span_bp``. The comparison is therefore an identity check, not an analogy.

    ⛔ **Projected through the SHIPPED ``_project_regions_to_loci``.** A second overlap-share
    implementation here would drift from the one under test and the difference would read as
    assembler error (TRAPS: a-test-that-redefines).

    ⚠ **The RNA arm is spliced-INCLUSIVE and is a BOUND, not the target.** ``rna_prior_count``
    withholds spliced mass; ``node_start_count`` has no splice bit. Returns
    ``(f_gdna, f_rna_upper, dropped)`` where ``dropped`` is the per-origin fragment count whose start
    node overlaps no locus — intergenic, correctly outside every prior, and reported so that
    ``Σ F + dropped == the library total`` is checkable rather than assumed.
    """
    starts = {k: np.asarray(oracle.parts[k].node_start_count, np.float64) for k in ORIGINS}
    rna = starts["mrna"] + starts["nrna"]
    proj = PRIORS._project_regions_to_loci(
        region_arrays, multi_loci, len(multi_loci), {"gdna": starts["gdna"], "rna": rna}
    )
    dropped = {
        "gdna": float(starts["gdna"].sum() - proj["gdna"].sum()),
        "rna": float(rna.sum() - proj["rna"].sum()),
    }
    return proj["gdna"], proj["rna"], dropped


# ── one condition ────────────────────────────────────────────────────────────────────────────────


@dataclass
class ConditionResult:
    """Everything one condition produced. The arrays are kept so a gate can re-score them."""

    condition: str
    n_loci: int
    priors: dict  #: arm name -> LocusPriors
    f_gdna: np.ndarray
    f_rna_upper: np.ndarray
    f_dropped: dict
    noop_identical: dict  #: field -> bool
    #: ⭐ is gdna_eff_len clamped by an incidence sum? See :func:`eff_len_inflation`
    eff_len: dict
    drain: dict  #: the measured drain caveat
    library: dict  #: condition-level totals, for the header row
    seconds: float
    #: ⭐ The three inputs a gate needs to RE-DERIVE any arm above rather than trust the arrays it was
    #: handed. Kept deliberately: a falsification test that can only read the outputs can check they
    #: are self-consistent and never that they are the outputs of the shipped code.
    oracle: OracleTruth = None
    region_arrays: object = None
    multi_loci: list = None
    calibration: object = None


def _oracle_parts(bam, index, scan, pipeline_config, work_dir, tag, cache_root):
    """The three origin partitions, from the cache when it is valid — otherwise split and scan.

    ⛔ Keyed by the SHIPPED ``read_scan_cache``, never a home-made key: ``reach`` is covered by no other
    hash, so a rebuilt index would verify clean against one. Same argument as
    ``pass0_vs_oracle.load_or_build_oracle``, which this deliberately mirrors.
    """
    if cache_root is not None:
        dirs = {k: Path(cache_root) / tag / k for k in ORIGINS}
        try:
            return {k: read_scan_cache(dirs[k], index, scan).payload for k in ORIGINS}
        except (FileNotFoundError, KeyError, ScanCacheKeyError):
            pass
    paths, _counts = _split_bam(bam, Path(work_dir), tag)
    parts = {}
    for origin in ORIGINS:
        _s, strand_model, _b, payload = scan_and_buffer(paths[origin], index, scan)
        parts[origin] = payload
        if cache_root is not None:
            d = Path(cache_root) / tag / origin
            write_scan_cache(d, payload=payload, strand_model=strand_model, index=index,
                             bam=paths[origin], scan_config=scan)
    return parts


def _calibrate_and_prior(payload, strand_model, buffer, stats, index, ra, pipeline_config):
    """calibrate → score → ``(calibration, fl, multi_loci, LocusPriors)`` on ONE payload."""
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.junction_opportunity import crossing_probability_from_index
    from rigel.calibration.splice_graph import (
        build_edge_flags_array,
        build_junction_geometry_arrays,
    )

    max_size = int(payload.max_length)
    fl = build_fl_models(
        payload,
        junction_opportunity=crossing_probability_from_index(index, max_size),
        gdna_opportunity=gdna_opportunity_from_index(index, max_size),
    )
    cal = calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=strand_model,
        gdna_fl_pmf=fl.gdna_pmf,
        rna_fl_pmf=fl.rna_pmf,
        config=pipeline_config.calibration,
        junctions=build_junction_geometry_arrays(index),
        edge_flags=build_edge_flags_array(index),
    )
    multi_loci, priors = capture_priors(
        buffer, index, strand_model, fl, ra, stats, cal, payload, pipeline_config
    )
    return cal, fl, multi_loci, priors


def measure_condition(bam, index, pipeline_config, work_dir, tag, *, oracle_cache=None,
                      drained_arm=True, emit_masses=None) -> ConditionResult:
    """Scan once, build P / O / F, and price the drain caveat.

    ⚠ **Two scans of the same BAM when ``drained_arm`` is on, and they are not interchangeable.** The
    fragment BUFFER is consumed by ``quant_from_buffer`` (and freed by it), so a second prior needs a
    second scan; the payload could be cached but the buffer cannot.
    """
    start = time.perf_counter()
    scan = dataclasses.replace(pipeline_config.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    ra = RegionArrays.from_index(index)

    stats, strand_model, buffer, payload = scan_and_buffer(bam, index, scan)
    n_held = int(payload.deferred.n_fragments)
    parts = _oracle_parts(bam, index, scan, pipeline_config, work_dir, tag, oracle_cache)
    oracle = OracleTruth.from_parts(payload, parts)  # sum-to-full, UNDRAINED, raises if invalid

    cal, _fl, multi_loci, p_arm = _calibrate_and_prior(
        payload, strand_model, buffer, stats, index, ra, pipeline_config
    )
    o_arm, noop = oracle_priors(oracle, cal, ra, multi_loci)
    s_arm, _shares = share_priors(oracle, cal, ra, multi_loci)
    f_gdna, f_rna_upper, f_dropped = fragment_truth(oracle, ra, multi_loci)
    eff_len = eff_len_inflation(cal, ra, multi_loci)

    noop_identical = {
        f: bool(np.array_equal(getattr(noop, f), getattr(p_arm, f))) for f in PRIOR_FIELDS
    }

    if emit_masses is not None:
        # locus-FREE, on purpose: item 3 must rebuild the oracle prior on the loci its OWN run
        # produced. The locus partition is a function of the scoring stage, so an array keyed by
        # multi_locus_id is not portable between runs.
        d = Path(emit_masses) / tag
        d.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(d / "oracle_masses.npz", **oracle.override_masses(ra))

    # ⭐ THE DRAIN CAVEAT, PRICED. Two numbers, and they answer different questions: whether a DRAINED
    # oracle is admissible at all (it is not, and the leak says why), and how far the shipped prior
    # moves when the drain runs (which is what "undrained" costs this table).
    drain: dict = {"n_held": n_held, "n_ambiguous": None, "gdna_spliced_leak": None,
                   "p_drained_vs_p": None}
    if drained_arm and n_held:
        from rigel.second_pass import drain as sp_drain, lift_choices

        lift: dict = {}
        # ⭐ ``drain`` is pure — payload in, payload out — so ``payload`` is still the undrained
        # tally the arms above were scored on, and ``payload_d`` is what production would calibrate.
        payload_d = _drain_side_buffer(
            payload, index, strand_model, seed=pipeline_config.second_pass_seed, _lift=lift
        )
        lifted, n_amb = lift_choices(
            lift["undrained"], [parts[k] for k in ORIGINS], lift["choices"]
        )
        parts_d = {
            k: sp_drain(parts[k], ch, node_types=lift["node_types"], junctions=lift["junctions"])
            for k, ch in zip(ORIGINS, lifted)
        }
        drain["n_ambiguous"] = int(n_amb)
        drain["gdna_spliced_leak"] = int(
            np.asarray(parts_d["gdna"].edge_spliced_count, np.int64).sum()
            + np.asarray(parts_d["gdna"].sj_count, np.int64).sum()
        )
        # ⚠ A SECOND SCAN, and it buys one thing: a fresh fragment BUFFER. The first was handed to
        # ``quant_from_buffer`` and a buffer is not re-scannable; the payload is cacheable and the
        # buffer is not. The scan is a function of the BAM, the index and the scan config alone, so
        # this reproduces the first one exactly.
        _s2, sm2, buf2, _p2 = scan_and_buffer(bam, index, scan)
        _c2, _f2, ml2, p_drained = _calibrate_and_prior(
            payload_d, sm2, buf2, _s2, index, ra, pipeline_config
        )
        # ⚠ The drained run builds its OWN loci, and they need not match: build_multi_loci depends on
        # the scored fragments, which depend on the calibration, which depends on the drain. A
        # mismatch is REPORTED, never silently index-aligned.
        drain["n_loci_drained"] = len(ml2)
        drain["loci_match"] = len(ml2) == len(multi_loci)
        if drain["loci_match"]:
            drain["p_drained_vs_p"] = {
                "gdna": dataclasses.asdict(
                    score_arm(p_drained.gdna_prior_count, p_arm.gdna_prior_count)
                ),
                "rna": dataclasses.asdict(
                    score_arm(p_drained.rna_prior_count, p_arm.rna_prior_count)
                ),
            }

    g = float(np.asarray(oracle.parts["gdna"].node_start_count, np.int64).sum())
    r = float(
        np.asarray(oracle.parts["mrna"].node_start_count, np.int64).sum()
        + np.asarray(oracle.parts["nrna"].node_start_count, np.int64).sum()
    )
    return ConditionResult(
        condition=tag,
        n_loci=len(multi_loci),
        priors={"P": p_arm, "O": o_arm, "S": s_arm},
        f_gdna=f_gdna,
        f_rna_upper=f_rna_upper,
        f_dropped=f_dropped,
        noop_identical=noop_identical,
        eff_len=eff_len,
        drain=drain,
        library={
            "true_gdna_fragments": g,
            "true_rna_fragments": r,
            "true_f_gdna": g / (g + r) if (g + r) > 0 else 0.0,
        },
        seconds=time.perf_counter() - start,
        oracle=oracle,
        region_arrays=ra,
        multi_loci=multi_loci,
        calibration=cal,
    )


# ── reporting ────────────────────────────────────────────────────────────────────────────────────


def stratum(cond: str) -> tuple[str, str]:
    """The panel's two binary axes. Same definition as ``arm_score.stratum`` — one home would be
    better, but that module is a scorer for a different file format and importing it here would drag
    its ``$RIGEL_ARMS`` path resolution in."""
    return (
        "stranded" if "ss_0.99" in cond else "unstranded",
        "capture ON" if "capture_on" in cond else "capture OFF",
    )


def is_zero_gdna(cond: str) -> bool:
    """``g00`` — the owner-required ZERO-gDNA control. Truth is exactly 0, so every gDNA fragment in
    the prior there is a false positive with nothing to cancel it, and a relative change is unbounded.
    Reported on its own row, never inside ALL."""
    return "_g00_" in cond


def _agg(scores):
    """Sum a list of ``ArmScore`` into one. ⛔ The rates are RE-DERIVED from the summed totals, never
    averaged: a mean of ratios over conditions of wildly different depth is a number with no
    consumer (TRAPS: never-pool-the-strata's third way)."""
    scores = [s for s in scores if s is not None]
    if not scores:
        return None
    return ArmScore(
        n_loci=sum(s.n_loci for s in scores),
        n_claiming=sum(s.n_claiming for s in scores),
        total_arm=sum(s.total_arm for s in scores),
        total_ref=sum(s.total_ref for s in scores),
        net_err=sum(s.net_err for s in scores),
        abs_err=sum(s.abs_err for s in scores),
        over_call=sum(s.over_call for s in scores),
        under_call=sum(s.under_call for s in scores),
    )


_STRATA = (
    ("stranded", "capture OFF"),
    ("stranded", "capture ON"),
    ("unstranded", "capture OFF"),
    ("unstranded", "capture ON"),
)

#: Every selection every table prints, in order. ⭐ One list, so a stratum added here appears on every
#: table at once and cannot appear on some of them (which is how two tables come to disagree about
#: what "ALL" means).
_SELECTIONS = (
    ("ALL (g00 excluded)", lambda c: not is_zero_gdna(c)),
    *((" x ".join(st), (lambda c, st=st: stratum(c) == st and not is_zero_gdna(c)))
      for st in _STRATA),
    (None, None),  # a rule line
    ("⛔ g00 ZERO-gDNA control", is_zero_gdna),
)


def report(rows: list[dict]) -> None:
    """The whole report, from the per-condition JSON — **the only report path there is.**

    ⭐⭐ **It reads JSON in the serial case too, and that is the point.** The first version had a rich
    in-process report and a reduced one for the sharded path, printing the same tables from two
    different sources. That is `TRAPS: a-test-that-redefines` with a report in place of a test: the
    numbers a `--jobs 1` run printed and the numbers a `--jobs 6` run printed were produced by
    different code and nothing compared them. Everything either table needed was already in
    :func:`to_json`, so there is one table builder and the shard merge is not a special case.
    """
    rows = sorted(rows, key=lambda r: r["condition"])
    print()
    print("=" * 104)
    print("  ⭐⭐⭐ CALIBRATION'S ENDPOINT vs THE ORACLE — LocusPriors, in FRAGMENTS")
    print(f"  {len(rows)} conditions   messages OFF   length_likelihood OFF   UNDRAINED "
          "(the drain caveat is priced below)")
    print("=" * 104)

    # ── the gate first. A table read before its gate is a table nobody checked. ──
    bad = [(r["condition"], f) for r in rows for f, ok in r["noop_identical"].items() if not ok]
    if bad:
        print(f"\n  ⛔⛔ NOOP GATE FAILED on {len(bad)} (condition, field) pairs — the override "
              "plumbing is NOT inert and every number below is unattributable:")
        for c, f in bad[:10]:
            print(f"       {c}  {f}")
        raise SystemExit(2)
    print(f"\n  ✅ noop gate: re-injecting the shipped masses reproduces P byte-identically on all "
          f"{len(rows)} x {len(PRIOR_FIELDS)} arrays (TRAPS: byte-identity-gate)")

    def arm_table(title: str, key: str, arm: str, note: str = "") -> None:
        print()
        print(f"  {title}")
        if note:
            print(f"  {note}")
        print(f"    {'stratum':<26} {'ref total':>14} {'arm total':>14} {'Σ|Δ|':>14} "
              f"{'rel':>8} {'net':>14} {'canc':>7}")
        print("    " + "-" * 102)
        for label, sel in _SELECTIONS:
            if label is None:
                print("    " + "-" * 102)
                continue
            s = _agg([ArmScore(**r[key][arm]) for r in rows if sel(r["condition"])])
            if s is None:
                print(f"    {label:<26} {'(empty)':>14}")
                continue
            print(f"    {label:<26} {s.total_ref:>14,.0f} {s.total_arm:>14,.0f} "
                  f"{s.abs_err:>14,.0f} {s.rel:>8.3f} {s.net_err:>+14,.0f} "
                  f"{s.cancellation:>7.1f}")

    _RNA_BOUND = ("⚠ F IS AN UPPER BOUND: node_start_count carries no splice bit, so it counts the "
                  "spliced fragments rna_prior_count deliberately withholds. A NEGATIVE net is "
                  "expected here and is not an error.")

    arm_table("① P vs O — CALIBRATION'S OWN ERROR (a perfect deconvolution, same assembler) · gDNA arm",
              "P_vs_O", "gdna")
    arm_table("   … RNA arm", "P_vs_O", "rna")
    arm_table("④ O vs F — THE ASSEMBLER'S OWN ERROR (truth masses in) · gDNA arm",
              "O_vs_F", "gdna",
              "⭐ F is EXACT here — gDNA does not splice. This prices the mass→density→fragment-count "
              "conversion, the projection and the pooling, alone.")
    arm_table("   … RNA arm", "O_vs_F", "rna", _RNA_BOUND)
    arm_table("⑤ P vs F — THE TOTAL PRIOR ERROR (what the EM is handed, vs the true counts) · gDNA arm",
              "P_vs_F", "gdna")
    arm_table("   … RNA arm", "P_vs_F", "rna", _RNA_BOUND)

    # ── ⑧ the two halves of ④, separated ──
    arm_table("⑧ S vs F — THE ASSEMBLER WITH PERFECT PER-COMPONENT SHARES · gDNA arm",
              "S_vs_F", "gdna",
              "⭐ Everything ④ measures EXCEPT the pooled share. The gap between ④ and ⑧ is what one "
              "share for two components costs.")
    arm_table("⑨ O vs S — THE POOLED SHARE'S OWN CONTRIBUTION, ISOLATED · gDNA arm",
              "O_vs_S", "gdna",
              "⛔ The locus TOTAL is conserved here to the fragment, so this error is purely "
              "compositional and no conservation gate can see it. Equal component lengths ⇒ identically "
              "zero, which is why the ladder is structurally blind to it.")

    # ── ⑩ is gdna_eff_len clamped by an incidence sum? ──
    print()
    print("  ⑩ ⭐ gdna_eff_len's CLAMP — is `span` a genomic extent or an INCIDENCE sum?")
    print(f"    {'stratum':<26} {'support/genomic':>16} {'nodes only':>12} {'Σ support':>16} "
          f"{'Σ genomic':>16}")
    print("    " + "-" * 92)
    for label, sel in _SELECTIONS:
        if label is None:
            print("    " + "-" * 92)
            continue
        sub_rows = [r["eff_len"] for r in rows if sel(r["condition"]) and r.get("eff_len")]
        if not sub_rows:
            continue
        print(f"    {label:<26} "
              f"{float(np.median([x['median_support_over_genomic'] for x in sub_rows])):>16.2f} "
              f"{float(np.median([x['median_node_over_genomic'] for x in sub_rows])):>12.2f} "
              f"{sum(x['total_support'] for x in sub_rows):>16,.0f} "
              f"{sum(x['total_genomic'] for x in sub_rows):>16,.0f}")
    print("    ⚠ `support/genomic` well above 1 means every interior line is adding ~mu_g − 1 to the "
          "locus's clamp.")
    print("    The EM divides the gDNA component's abundance by gdna_eff_len, so this is a direct "
          "scale error on a shipped number.")

    # ── ② the composition claim, which is what the EM feels ──
    print()
    print("  ② THE COMPOSITION CLAIM  phi = a_g/(a_g+a_r)  — P against O")
    print(f"    {'stratum':<26} {'n loci':>8} {'mwae_phi':>10} {'median':>10} "
          f"{'scale P':>14} {'scale O':>14} {'log10 P/O':>10}")
    print("    " + "-" * 96)
    for label, sel in _SELECTIONS:
        if label is None:
            print("    " + "-" * 96)
            continue
        sub = [r["composition_P_vs_O"] for r in rows if sel(r["condition"])]
        if not sub:
            continue
        w = np.array([x["weight"] for x in sub])
        tot = float(w.sum())
        mwae = float((w * np.array([x["mwae_phi"] for x in sub])).sum() / tot) if tot else np.nan
        med = float(np.median([x["median_phi_err"] for x in sub]))
        sp, so = sum(x["scale_arm"] for x in sub), sum(x["scale_ref"] for x in sub)
        lr = np.log10(sp / so) if sp > 0 and so > 0 else float("nan")
        print(f"    {label:<26} {sum(x['n_scored'] for x in sub):>8,} {mwae:>10.4f} {med:>10.4f} "
              f"{sp:>14,.0f} {so:>14,.0f} {lr:>+10.3f}")

    # ── ③ the third array ──
    print()
    print("  ③ gdna_eff_len — P against O, weighted by O's own gDNA prior count")
    print(f"    {'stratum':<26} {'n loci':>8} {'w rel err':>11} {'median rel':>11}")
    print("    " + "-" * 60)
    for label, sel in _SELECTIONS:
        if label is None:
            print("    " + "-" * 60)
            continue
        sub = [r["eff_len_P_vs_O"] for r in rows if sel(r["condition"])]
        if not sub:
            continue
        print(f"    {label:<26} {sum(x['n_scored'] for x in sub):>8,} "
              f"{float(np.mean([x['w_rel_err'] for x in sub])):>11.4f} "
              f"{float(np.median([x['median_rel_err'] for x in sub])):>11.4f}")

    # ── ⑥ the drain caveat, as a number ──
    print()
    print("  ⑥ ⛔ THE DRAIN CAVEAT — measured, not assumed")
    d = [r["drain"] for r in rows if r["drain"].get("n_ambiguous") is not None]
    if not d:
        print("     (not measured — --no-drained-arm)")
    else:
        amb, held = sum(x["n_ambiguous"] for x in d), sum(x["n_held"] for x in d)
        leak = sum(x["gdna_spliced_leak"] for x in d)
        n_leak = sum(1 for x in d if x["gdna_spliced_leak"] > 0)
        print(f"     a DRAINED oracle is INADMISSIBLE here: {amb:,} of {held:,} held fragments "
              f"({amb / max(held, 1):.2%}) tie on the deferred bank's key across origins,")
        print(f"     depositing {leak:,} spliced records into the gdna partition on "
              f"{n_leak}/{len(d)} conditions. gDNA cannot splice; OracleTruth refuses it.")
        mism = [x for x in d if not x.get("loci_match", True)]
        if mism:
            print(f"     ⚠ {len(mism)} conditions build a DIFFERENT number of loci once drained "
                  "(the partition is a function of the scoring stage) — not compared.")
        ok = [x for x in d if x.get("p_drained_vs_p")]
        for arm in ("gdna", "rna"):
            a = _agg([ArmScore(**x["p_drained_vs_p"][arm]) for x in ok])
            if a:
                print(f"     what the drain costs this table, {arm:>4} arm: "
                      f"Σ|P_drained − P| = {a.abs_err:,.0f} on a P total of {a.total_ref:,.0f} "
                      f"({a.rel:.3%}), net {a.net_err:+,.0f}")

    # ── ⑦ the per-condition ladder, because a stratum total hides the shape ──
    print()
    print("  ⑦ PER CONDITION — the gDNA arm, P and O as fragment totals")
    print(f"    {'condition':<44} {'true f_g':>9} {'O_g':>13} {'P_g':>13} {'P/O':>8} "
          f"{'Σ|Δ|':>13} {'mwae_phi':>9}")
    print("    " + "-" * 116)
    for r in rows:
        g = ArmScore(**r["P_vs_O"]["gdna"])
        ratio = g.total_arm / g.total_ref if g.total_ref > 0 else float("nan")
        print(f"    {r['condition']:<44} {r['library']['true_f_gdna']:>9.4f} {g.total_ref:>13,.0f} "
              f"{g.total_arm:>13,.0f} {ratio:>8.3f} {g.abs_err:>13,.0f} "
              f"{r['composition_P_vs_O']['mwae_phi']:>9.4f}")

    print()
    print(f"  total wall clock {sum(r['seconds'] for r in rows):,.0f} s over {len(rows)} conditions")


def to_json(results: list[ConditionResult]) -> list[dict]:
    out = []
    for r in results:
        row: dict = {
            "condition": r.condition,
            "stratum": list(stratum(r.condition)),
            "n_loci": r.n_loci,
            "noop_identical": r.noop_identical,
            "drain": r.drain,
            "eff_len": r.eff_len,
            "library": r.library,
            "f_dropped": r.f_dropped,
            "seconds": r.seconds,
        }
        p, o = r.priors["P"], r.priors["O"]
        sa = r.priors["S"]
        for ref_name, gref, rref, arm in (
            ("P_vs_O", o.gdna_prior_count, o.rna_prior_count, p),
            ("O_vs_F", r.f_gdna, r.f_rna_upper, o),
            ("P_vs_F", r.f_gdna, r.f_rna_upper, p),
            # ⭐ the two halves of O_vs_F, separated
            ("S_vs_F", r.f_gdna, r.f_rna_upper, sa),
            ("O_vs_S", sa.gdna_prior_count, sa.rna_prior_count, o),
        ):
            row[ref_name] = {
                "gdna": dataclasses.asdict(score_arm(arm.gdna_prior_count, gref)),
                "rna": dataclasses.asdict(score_arm(arm.rna_prior_count, rref)),
            }
        row["composition_P_vs_O"] = score_composition(p, o)
        row["eff_len_P_vs_O"] = score_eff_len(p, o)
        out.append(row)
    return out


# ── main ─────────────────────────────────────────────────────────────────────────────────────────


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--oracle-cache", type=Path, default=None,
                    help="defaults to <suite>/oracle_cache when that directory exists")
    ap.add_argument("--work-dir", type=Path,
                    default=Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_prior_oracle")
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--emit-oracle-masses", type=Path, default=None,
                    help="write the oracle mass override per condition (locus-FREE) for item 3")
    ap.add_argument("--no-drained-arm", action="store_true",
                    help="skip the drain caveat measurement (halves the wall clock)")
    ap.add_argument("--jobs", type=int, default=1,
                    help="run this many conditions CONCURRENTLY by re-invoking on shards. The "
                         "conditions are independent, so this changes no number.")
    args = ap.parse_args()

    names = args.conditions or sorted(
        p.name for p in args.suite.iterdir() if (p / "sim_oracle.bam").is_file()
    )
    if not names:
        raise SystemExit(f"no conditions with a sim_oracle.bam under {args.suite}")
    cache = args.oracle_cache
    if cache is None and (args.suite / "oracle_cache").is_dir():
        cache = args.suite / "oracle_cache"

    if args.jobs > 1 and len(names) > 1:
        # ⭐ Shards, not threads: conditions share nothing but a read-only index and cache, and
        # re-invoking the single-process path keeps the measured code byte-for-byte the serial one.
        # ⚠ OMP_NUM_THREADS=1 is forced at import so the workers do not fight.
        shards = [s for s in (names[i:: args.jobs] for i in range(args.jobs)) if s]
        tmp = args.work_dir / "_shards"
        tmp.mkdir(parents=True, exist_ok=True)
        procs, outs = [], []
        for i, sh in enumerate(shards):
            o = tmp / f"{i}.json"
            outs.append(o)
            cmd = [sys.executable, str(Path(__file__).resolve()),
                   "--suite", str(args.suite), "--index", str(args.index),
                   "--work-dir", str(args.work_dir / f"shard{i}"), "--json", str(o),
                   "--conditions", *sh]
            if cache is not None:
                cmd += ["--oracle-cache", str(cache)]
            if args.no_drained_arm:
                cmd.append("--no-drained-arm")
            if args.emit_oracle_masses is not None:
                cmd += ["--emit-oracle-masses", str(args.emit_oracle_masses)]
            procs.append(subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                          stderr=subprocess.STDOUT, text=True))
        rc = 0
        for i, pr in enumerate(procs):
            out, _ = pr.communicate()
            if pr.returncode != 0:
                rc = pr.returncode
                print(f"  ⛔ shard {i} FAILED (rc={pr.returncode}):\n{out}", flush=True)
            else:
                print(f"  shard {i}: {len(shards[i])} conditions ok", flush=True)
        if rc:
            # ⛔ TRAPS: an-ablation-that-never-ran's shape — a short output file reads as a complete panel.
            raise SystemExit("a shard failed; refusing to report a partial panel")
        merged: list[dict] = []
        for o in outs:
            merged += json.loads(o.read_text())
        if args.json is not None:
            args.json.write_text(json.dumps(merged, indent=1))
        report(merged)
        return 0

    index = TranscriptIndex.load(str(args.index))
    pipeline_config = PipelineConfig()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for name in names:
        bam = str(args.suite / name / "sim_oracle.bam")
        print(f"  … {name}", flush=True)
        results.append(measure_condition(
            bam, index, pipeline_config, args.work_dir, name,
            oracle_cache=cache, drained_arm=not args.no_drained_arm,
            emit_masses=args.emit_oracle_masses,
        ))
    payload = to_json(results)
    if args.json is not None:
        args.json.write_text(json.dumps(payload, indent=1))
    report(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
