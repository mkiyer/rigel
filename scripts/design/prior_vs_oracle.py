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

⭐ **FIVE ARMS, AND THE DIFFERENCES BETWEEN THEM ARE THE WHOLE DIAGNOSTIC.**

===  ==========================================================  ===================================
 P   the SHIPPED prior                                           ``assemble_priors(cal, ra, loci)``
 O   the prior a PERFECT DECONVOLUTION would produce             truth masses, same assembler
 S   O, plus each component rescaled by its OWN true share        :func:`share_priors`
 Fo  ⭐ the EM's OWN candidate count, per origin                  :func:`overlap_truth`
 F   the per-locus START-BASE count                              ``region_start_count``, per origin
===  ==========================================================  ===================================

``P − O`` is calibration's own error — everything upstream of :func:`assemble_priors`.
``O − Fo`` is the ASSEMBLER's error — the mass → fragment-count conversion, the overlap projection
and the one pooled per-component share. ``O − S`` splits that share out of it and ``S − Fo`` is
everything else. They are different repairs in different files and pooling them would hide which.

⛔⛔ **``F`` IS NOT THE EM's TARGET, AND EVERY ``O − F`` / ``S − F`` NUMBER THIS FILE ONCE PRINTED WAS
SCORED AGAINST A QUANTITY NOTHING CONSUMES** (found and corrected 2026-08-08; the paragraph this
replaced claimed ``F`` was "EXACT on the gDNA arm ... with nothing to subtract"). ``region_start_count``
deposits at the region holding a fragment's FIRST BASE, so ``F`` credits a locus only for the fragments
that *start* inside it. The EM counts differently: ``n_gdna`` in
``em_solver.cpp:apply_grouped_prior_update`` is the soft count over that multi-locus's own EM UNITS,
and a fragment becomes one by being a scored CANDIDATE — which a fragment that starts in the
intergenic flank and reaches into the locus is. ``F`` drops exactly that straddling population.
⭐ ``F`` is retained, labelled, and priced on its own table, because the size of the correction is
itself a result.

⭐⭐ **``Fo`` (``F_overlap``) IS THAT COUNT, TAKEN FROM THE RUN'S OWN SCORING STAGE.** Every EM unit of
a multi-locus (``MultiLocus.unit_indices``, the array ``locus_partition`` scatters by) labelled with
the TRUE origin of its fragment (``_oracle.frag_id_origins``, joined on ``frag_id``). Nothing is
re-implemented: the unit→locus map is ``build_multi_loci``'s own and the origin is the simulator's read
name (TRAPS: a-test-that-redefines). ⛔ It is not circular — ``assemble_priors`` never sees a unit
count; grep it for ``unit_indices``.

⭐ **AND ``Fo`` GIVES THE RNA ARM AN EXACT TARGET FOR THE FIRST TIME.** ``rna_prior_count`` withholds
spliced mass, and ``region_start_count`` carries no splice bit — which is why ``F_rna`` could only ever
be an upper BOUND, "separating it needs a five-way origin×splice BAM split". A scored unit carries
``is_spliced``, so the split is free: ``Fo`` reports the RNA arm's target (**unspliced** RNA units)
and, separately, every RNA unit.

⛔⛔ **AND ``rna_unspliced`` IS THE ONE TO SCORE ``a_r`` AGAINST — ``rna_all`` IS A TRAP THIS FILE FELL
INTO ON THE DAY IT WAS WRITTEN.** A spliced unit never receives a gDNA candidate
(``em_solver.cpp``: ``has_gdna = !is_spliced && isfinite(gdna_ll)``), so spliced RNA does not compete
with gDNA and must not enter a prior that arbitrates that competition — putting it in would penalise
gDNA with fragments it could never have won. The prior's population is therefore **gDNA units +
UNSPLICED RNA units**, and against it the composition claim is exact to ≤ 5e-4. Scored against
``rna_all`` instead it reads a phantom +0.07…+0.10 tilt, which is the denominator and not the prior
(``TRAPS: score-the-consumers-own-count``, committed and then immediately repeated).

⛔ **O IS NOT AN ESTIMATOR AND MUST NEVER BECOME ONE** (the shape ``pass0_vs_oracle.py`` refuses too).
It is the shipped :func:`assemble_priors` fed the ONE lever that already exists —
``OracleTruth.override_masses``, the true per-object gDNA/RNA split — with every other input untouched.
Writing a "best possible prior" here would be the magic-number failure with an estimator in place of a
constant.

⚠ **WHAT ``F`` IS STILL GOOD FOR, stated so it is not read as the target again.** It is the only arm
that needs no scoring stage — one deposit per accepted fragment at the region holding its first base,
projected through the SHIPPED ``priors._project_regions_to_loci`` — so it is the arm that keeps working
when the question is about the *projection* rather than about the EM. ``Fo − F`` is the straddling
population, per locus, and table ⑪ is it.

⛔⛔ **UNDRAINED ON EVERY ARM, AND THAT IS FORCED — WITH A MEASUREMENT, NOT AN ASSUMPTION.**
Production drains the side buffer before calibration. ``lift_choices`` exists to make a DRAINED oracle
valid, and on this panel it is **not**: measured at ``g25 ss0.50 capture_on``, 6,450 of 209,658 held
fragments (3.08 %) tie on the deferred bank's canonical key across origins, and replaying the whole's
choices then deposits 1,026 ``boundary_spliced`` + 526 ``sj`` records inside the **gdna** partition —
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
    python scripts/design/prior_vs_oracle.py --conditions gdna_g50_ss_0.50_nrna_mid_capture_on
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

from _oracle import (  # noqa: E402
    ORIGIN_CODE,
    ORIGINS,
    OracleTruth,
    _split_bam,
    check_walk_alignment,
    frag_id_origins,
)

import rigel.calibration.priors as PRIORS  # noqa: E402
import rigel.locus as LOCUS  # noqa: E402
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
    "mass_gdna_region",
    "mass_rna_region",
    "mass_gdna_boundary",
    "mass_rna_boundary",
    "mass_rna_spliced_boundary",
    "count_rna_sj",
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
    """Run the production quant path far enough to get ``(multi_loci, LocusPriors, units)``, then STOP.

    ⭐ **The loci are the production loci, not a re-derivation.** ``build_multi_loci`` builds
    connected components of transcripts linked by SCORED fragments, so the locus partition is a
    function of the scoring stage and cannot be reconstructed from the index alone. Wrapping
    :func:`~rigel.calibration.priors.assemble_priors` — which ``quant_from_buffer`` imports
    function-locally, so patching the module attribute is picked up at call time — takes both objects
    from the call production itself makes (TRAPS: a-test-that-redefines).

    ⭐ **``units`` comes from the same run, through the same trick.** ``build_multi_loci`` is the one
    call that sees the scored CSR, so wrapping it yields the two per-unit arrays :func:`overlap_truth`
    needs — ``frag_ids`` (the identity that joins to origin truth) and ``is_spliced`` (the bit that
    separates the assembler's RNA target from the EM's RNA population). ⛔ Re-scoring the buffer a
    second time to obtain them would be a different scoring stage than the one that built the loci.

    ⛔ The sentinel exception is what makes this affordable: the per-locus EM is the single most
    expensive stage (47 % of a condition) and this instrument does not read its output. ⚠ Item 3 —
    injecting the oracle prior and re-quantifying — needs the EM and must NOT reuse this path.
    """

    class _StopAfterPriors(Exception):
        pass

    captured: dict = {}
    original = PRIORS.assemble_priors
    original_ml = LOCUS.build_multi_loci

    def _wrapper(cal, ra, multi_loci):
        captured["multi_loci"] = multi_loci
        captured["priors"] = original(cal, ra, multi_loci)
        raise _StopAfterPriors

    def _ml_wrapper(em_data, idx):
        captured["units"] = {
            "frag_ids": np.array(em_data.frag_ids, dtype=np.int64, copy=True),
            "is_spliced": np.array(em_data.is_spliced, dtype=bool, copy=True),
            "n_units": int(em_data.n_units),
        }
        return original_ml(em_data, idx)

    PRIORS.assemble_priors = _wrapper
    LOCUS.build_multi_loci = _ml_wrapper
    try:
        quant_from_buffer(
            buffer, index, strand_models, fl, region_arrays, stats, calibration, payload,
            em_config=pipeline_config.em, scoring=pipeline_config.scoring,
        )
    except _StopAfterPriors:
        pass
    finally:
        PRIORS.assemble_priors = original
        LOCUS.build_multi_loci = original_ml
    if "priors" not in captured or "units" not in captured:
        # ⛔ TRAPS: an-ablation-that-never-ran. ``quant_from_buffer`` returns early on an empty unit set, and a
        # silently-absent capture would read here as a condition with no loci rather than as a
        # harness that never fired.
        raise RuntimeError(
            "assemble_priors was never called — quant_from_buffer returned before the prior stage "
            "(no EM units, or no multi-loci). This is not a condition with zero error."
        )
    return captured["multi_loci"], captured["priors"], captured["units"]


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
    """**S** — the O arm, plus each component rescaled by its OWN true per-boundary share.

    ⭐⭐ **WHY THIS ARM EXISTS.** ``assemble_priors`` rescales BOTH components at a boundary by ONE pooled
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
    shares = oracle.component_shares()
    truth_cal = dataclasses.replace(calibration, **oracle.override_masses(region_arrays))
    gdna = PRIORS.assemble_priors(
        dataclasses.replace(truth_cal, boundary_mass_per_crossing=shares["gdna"]),
        region_arrays, multi_loci,
    )
    rna = PRIORS.assemble_priors(
        dataclasses.replace(truth_cal, boundary_mass_per_crossing=shares["rna"]),
        region_arrays, multi_loci,
    )
    return PRIORS.LocusPriors(
        gdna_prior_count=gdna.gdna_prior_count,
        rna_prior_count=rna.rna_prior_count,
        gdna_eff_len=gdna.gdna_eff_len,
    ), shares


def eff_len_inflation(calibration, region_arrays, multi_loci) -> dict:
    """⭐ Is ``gdna_eff_len`` clamped by an INCIDENCE-support sum rather than the genomic span?

    ``assemble_priors`` clamps ``gdna_eff_len`` to ``span = Σ share·(S_region + S_boundary)``. ``S_boundary`` is
    ``E_g[w − 1] ≈ mu_g − 1`` PER BOUNDARY, so every interior boundary adds most of a fragment length to a
    locus whose regions may be a few hundred bases — an incidence-like sum, not a genomic extent. The EM
    divides the gDNA component's abundance by this array, so an inflation here is a direct scale error
    on one of the three numbers calibration ships.

    ⚠ Reports the ratio to the locus's GENOMIC span, mass-weighted by the gDNA prior, so the number is
    what the consumer feels rather than what an unweighted locus average would say
    (``TRAPS: weight-it-like-the-consumer``).
    """
    # ⭐ Regions and boundaries are projected on their OWN axes, exactly as `assemble_priors` does — no boundary is
    # folded onto a flank region. Re-deriving the fold here would measure a span the assembler no longer
    # builds (`TRAPS: a-test-that-redefines`).
    n_loci = len(multi_loci)
    region_support = np.maximum(np.asarray(calibration.gdna_region_eff_len, np.float64), 0.0)
    boundary_support = np.maximum(np.asarray(calibration.gdna_boundary_eff_len, np.float64), 0.0)
    proj = PRIORS._project_regions_to_loci(
        region_arrays, multi_loci, n_loci,
        {
            "region_only": region_support,
            "genomic": np.asarray(region_arrays.region_size_bp, np.float64),
        },
    )
    e_idx, e_lid, e_w = PRIORS._boundary_locus_shares(region_arrays, multi_loci, n_loci)
    proj["support"] = proj["region_only"] + PRIORS._sum_by_locus(
        e_idx, e_lid, e_w, boundary_support, n_loci
    )
    live = proj["genomic"] > 0
    ratio = np.divide(proj["support"], proj["genomic"], out=np.ones_like(proj["support"]), where=live)
    region_ratio = np.divide(proj["region_only"], proj["genomic"],
                           out=np.ones_like(proj["support"]), where=live)
    return {
        "median_support_over_genomic": float(np.median(ratio[live])) if live.any() else float("nan"),
        "median_region_over_genomic": float(np.median(region_ratio[live])) if live.any() else float("nan"),
        "total_support": float(proj["support"].sum()),
        "total_genomic": float(proj["genomic"].sum()),
    }


@dataclass(frozen=True, slots=True)
class OverlapTruth:
    """**Fo** — the EM's OWN per-locus candidate count, by TRUE origin. See :func:`overlap_truth`."""

    gdna: np.ndarray  #: float64[n_loci] — the target of ``gdna_prior_count``
    rna_unspliced: np.ndarray  #: float64[n_loci] — the target of ``rna_prior_count``
    rna_all: np.ndarray  #: float64[n_loci] — every RNA unit, the EM's own RNA evidence
    diag: dict  #: the accounting a caller must report — see :func:`overlap_truth`


def unit_origins(unit_frag_ids: np.ndarray, frag_origin: np.ndarray) -> np.ndarray:
    """``int8[n_units]`` — each EM unit's TRUE origin, joined on ``frag_id``.

    ⛔ **The range check is the whole gate and it must ABORT.** ``frag_origin`` is indexed BY
    ``frag_id``, so a walk that grouped the BAM differently than the scanner did produces an array of
    the wrong length — and numpy would happily wrap a negative index or raise a bare ``IndexError``
    that reads as a bug in this file rather than as a broken join.
    """
    ids = np.asarray(unit_frag_ids, np.int64)
    n = int(np.asarray(frag_origin).shape[0])
    if ids.size and (int(ids.min()) < 0 or int(ids.max()) >= n):
        raise RuntimeError(
            f"unit frag_ids span [{ids.min()}, {ids.max()}] but the origin walk found {n} fragments. "
            "The walk did not reproduce the scanner's frag_id — every origin label below would be "
            "shifted, and a shifted label is a plausible number."
        )
    return np.asarray(frag_origin, np.int8)[ids]


def overlap_truth(multi_loci, unit_origin: np.ndarray, unit_is_spliced: np.ndarray,
                  n_units: int, walk: dict) -> OverlapTruth:
    """**Fo** — how many of a multi-locus's EM CANDIDATES were truly gDNA, and truly RNA.

    ⭐⭐ **This is the quantity the EM's prior is added to**, and it is a count of UNITS, not of
    genomic overlaps and not of start positions: ``em_solver.cpp:apply_grouped_prior_update`` forms
    ``G = n_gdna + a_g`` where ``n_gdna`` is the soft gDNA count over the units
    ``locus_partition`` handed that locus — i.e. over ``MultiLocus.unit_indices`` exactly.

    ⭐ **Every input is somebody else's output.** ``unit_indices`` is ``build_multi_loci``'s,
    ``unit_is_spliced`` is the scoring stage's, and ``unit_origin`` is the simulator's read name joined
    on ``frag_id`` by :func:`unit_origins`. Nothing here re-derives a locus, a candidate set or a
    splice call (TRAPS: a-test-that-redefines).

    ⚠ ``unit_origin`` arrives ALREADY JOINED rather than as ``(frag_ids, frag_origin)``, so a cache can
    store one int8 per unit instead of an int64 ``frag_id`` per unit plus the whole per-fragment walk.
    The join's own gate lives in :func:`unit_origins`, where its inputs are.

    ⛔ **THREE ARRAYS BECAUSE THERE ARE THREE POPULATIONS AND THEY ARE NOT INTERCHANGEABLE.**

    * ``gdna`` — the target of ``gdna_prior_count``. Exact: a spliced unit cannot be gDNA (its
      ``gdna_log_lik`` is ``-inf``) and gDNA cannot splice, so nothing has to be withheld.
    * ``rna_unspliced`` — the target of ``rna_prior_count``, which withholds spliced mass.
    * ``rna_all`` — every RNA unit, which is what the EM's ``n_rna`` sees from the unit axis.

    ⚠ ``diag`` carries the accounting, and a caller must report it: ``spliced_gdna_units`` is the
    join's SECONDARY diagnostic (gDNA cannot splice; the hard gate is ``_oracle.check_walk_alignment``
    and it is a count identity, because this one is blind to a slip smaller than a population block),
    ``orphan_units`` counts units no multi-locus claimed, and
    ``nonunit_fragments`` is per origin the fragments that never became a candidate at all —
    intergenic gDNA, gated fragments, and the deterministic spliced-unambig RNA that bypasses the EM.
    ``Σ Fo[origin] + nonunit_fragments[origin] == walk["totals"][origin]`` is then checkable rather
    than assumed.
    """
    n_loci = len(multi_loci)
    n_units = int(n_units)
    origin = np.asarray(unit_origin, np.int8)
    spliced = np.asarray(unit_is_spliced, bool)
    if origin.shape != (n_units,) or spliced.shape != (n_units,):
        raise RuntimeError(
            f"per-unit arrays are {origin.shape} / {spliced.shape} against n_units={n_units} — these "
            "are not the same scoring stage's units."
        )

    lid = np.full(n_units, -1, dtype=np.int64)
    for m in multi_loci:
        lid[np.asarray(m.unit_indices, np.int64)] = m.multi_locus_id
    claimed = lid >= 0

    def per_locus(sel: np.ndarray) -> np.ndarray:
        return np.bincount(lid[sel & claimed], minlength=n_loci).astype(np.float64)

    is_gdna = origin == ORIGIN_CODE["gdna"]
    is_rna = ~is_gdna
    gdna = per_locus(is_gdna)
    rna_all = per_locus(is_rna)
    rna_unspliced = per_locus(is_rna & ~spliced)

    counted = {"gdna": float(gdna.sum()), "rna": float(rna_all.sum())}
    diag = {
        "n_units": n_units,
        "orphan_units": int((~claimed).sum()),
        "spliced_gdna_units": int((is_gdna & spliced).sum()),
        "spliced_rna_units": int((is_rna & spliced).sum()),
        "unit_totals": {"gdna": counted["gdna"], "rna": counted["rna"]},
        "nonunit_fragments": {
            "gdna": float(walk["totals"]["gdna"]) - counted["gdna"],
            "rna": float(walk["totals"]["mrna"] + walk["totals"]["nrna"]) - counted["rna"],
        },
        "walk": walk,
    }
    return OverlapTruth(gdna=gdna, rna_unspliced=rna_unspliced, rna_all=rna_all, diag=diag)


def fragment_truth(oracle: OracleTruth, region_arrays, multi_loci):
    """**F** — the per-locus true count of fragments whose FIRST BASE lands in the locus.

    ⛔⛔ **THIS IS NOT THE PRIOR'S TARGET — :func:`overlap_truth` IS.** ``region_start_count`` deposits at
    the region holding a fragment's first base, so a fragment that starts in the intergenic flank and
    reaches into the locus is a candidate the EM counts and a fragment ``F`` does not. Until
    2026-08-08 this docstring asserted the opposite ("the same quantity ``assemble_priors`` targets, by
    construction"), on an argument about ``rho_c · span_bp`` — a rule the assembler no longer computes.
    ``Fo − F`` is the straddling population and is reported as such.

    ⭐ **What it is still the right instrument for.** ``region_start_count`` is the accumulator's one real
    invariant — ``Σ region_start_count == qc.deposited``, one increment per accepted fragment — so ``F``
    is the only per-locus truth that needs no scoring stage, no candidate set and no EM. That makes it
    the arm to reach for when the question is about the *projection*, and it is why it survives here
    rather than being deleted.

    ⛔ **Projected through the SHIPPED ``_project_regions_to_loci``.** A second overlap-share
    implementation here would drift from the one under test and the difference would read as
    assembler error (TRAPS: a-test-that-redefines).

    ⚠ **The RNA arm is spliced-INCLUSIVE and is a BOUND.** ``rna_prior_count`` withholds spliced mass;
    ``region_start_count`` has no splice bit. ⭐ ``overlap_truth`` does have one, so the bound is no
    longer the best available RNA target — this one is kept only as ``F``'s own RNA companion. Returns
    ``(f_gdna, f_rna_upper, dropped)`` where ``dropped`` is the per-origin fragment count whose start
    region overlaps no locus — intergenic, correctly outside every prior, and reported so that
    ``Σ F + dropped == the library total`` is checkable rather than assumed.
    """
    parts = {k: np.asarray(oracle.parts[k].region_start_count, np.float64) for k in ORIGINS}
    g = parts["gdna"]
    r = parts["mrna"] + parts["nrna"]
    proj = PRIORS._project_regions_to_loci(
        region_arrays, multi_loci, len(multi_loci), {"gdna": g, "rna": r}
    )
    dropped = {
        "gdna": float(g.sum() - proj["gdna"].sum()),
        "rna": float(r.sum() - proj["rna"].sum()),
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
    #: ⭐ **Fo** — the EM's own candidate count by true origin. THE reference for every arm.
    overlap: OverlapTruth
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
    #: the Fo arm's two raw inputs, for the same reason — a gate that cannot re-tally them cannot
    #: perturb the join, and the join is where a silent one-fragment shift would live
    units: dict = None
    frag_origin: np.ndarray = None


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
    from rigel.calibration.sj_opportunity import crossing_probability_from_index
    from rigel.calibration.splice_graph import (
        build_boundary_flags_array,
        build_sj_geometry_arrays,
    )

    max_size = int(payload.max_length)
    fl = build_fl_models(
        payload,
        sj_opportunity=crossing_probability_from_index(index, max_size),
        gdna_opportunity=gdna_opportunity_from_index(index, max_size),
    )
    cal = calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=strand_model,
        gdna_fl_pmf=fl.gdna_pmf,
        rna_fl_pmf=fl.rna_pmf,
        config=pipeline_config.calibration,
        sj=build_sj_geometry_arrays(index),
        boundary_flags=build_boundary_flags_array(index),
    )
    multi_loci, priors, units = capture_priors(
        buffer, index, strand_model, fl, ra, stats, cal, payload, pipeline_config
    )
    return cal, fl, multi_loci, priors, units


def measure_condition(bam, index, pipeline_config, work_dir, tag, *, oracle_cache=None,
                      drained_arm=True, emit_masses=None) -> ConditionResult:
    """Scan once, build P / O / S / Fo / F, and price the drain caveat.

    ⚠ **Two scans of the same BAM when ``drained_arm`` is on, and they are not interchangeable.** The
    fragment BUFFER is consumed by ``quant_from_buffer`` (and freed by it), so a second prior needs a
    second scan; the payload could be cached but the buffer cannot.

    ⚠ **Plus one pysam WALK of the same BAM**, for the ``frag_id → true origin`` key the ``Fo`` arm
    joins on. Deliberately uncached: it is ~30 s at 10 M fragments against the two scans it sits beside,
    and a cache keyed on anything weaker than the scan cache's own manifest is how a stale truth array
    gets read as a result.
    """
    start = time.perf_counter()
    scan = dataclasses.replace(pipeline_config.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    ra = RegionArrays.from_index(index)

    stats, strand_model, buffer, payload = scan_and_buffer(bam, index, scan)
    n_held = int(payload.deferred.n_fragments)
    parts = _oracle_parts(bam, index, scan, pipeline_config, work_dir, tag, oracle_cache)
    oracle = OracleTruth.from_parts(payload, parts)  # sum-to-full, UNDRAINED, raises if invalid
    frag_origin, walk = frag_id_origins(bam, scan)
    # ⛔ Before anything is scored: the walk must have issued the same frag_ids the scan did.
    check_walk_alignment(walk, stats)

    cal, _fl, multi_loci, p_arm, units = _calibrate_and_prior(
        payload, strand_model, buffer, stats, index, ra, pipeline_config
    )
    o_arm, noop = oracle_priors(oracle, cal, ra, multi_loci)
    s_arm, _shares = share_priors(oracle, cal, ra, multi_loci)
    f_gdna, f_rna_upper, f_dropped = fragment_truth(oracle, ra, multi_loci)
    overlap = overlap_truth(
        multi_loci, unit_origins(units["frag_ids"], frag_origin),
        units["is_spliced"], units["n_units"], walk,
    )
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
            k: sp_drain(parts[k], ch, region_types=lift["region_types"], sj=lift["sj"])
            for k, ch in zip(ORIGINS, lifted)
        }
        drain["n_ambiguous"] = int(n_amb)
        drain["gdna_spliced_leak"] = int(
            np.asarray(parts_d["gdna"].boundary_spliced_count, np.int64).sum()
            + np.asarray(parts_d["gdna"].sj_count, np.int64).sum()
        )
        # ⚠ A SECOND SCAN, and it buys one thing: a fresh fragment BUFFER. The first was handed to
        # ``quant_from_buffer`` and a buffer is not re-scannable; the payload is cacheable and the
        # buffer is not. The scan is a function of the BAM, the index and the scan config alone, so
        # this reproduces the first one exactly.
        _s2, sm2, buf2, _p2 = scan_and_buffer(bam, index, scan)
        _c2, _f2, ml2, p_drained, _u2 = _calibrate_and_prior(
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

    g = float(np.asarray(oracle.parts["gdna"].region_start_count, np.int64).sum())
    r = float(
        np.asarray(oracle.parts["mrna"].region_start_count, np.int64).sum()
        + np.asarray(oracle.parts["nrna"].region_start_count, np.int64).sum()
    )
    return ConditionResult(
        condition=tag,
        n_loci=len(multi_loci),
        priors={"P": p_arm, "O": o_arm, "S": s_arm},
        f_gdna=f_gdna,
        f_rna_upper=f_rna_upper,
        f_dropped=f_dropped,
        overlap=overlap,
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
        units=units,
        frag_origin=frag_origin,
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
    (None, None),  # a rule boundary
    ("⛔ g00 ZERO-gDNA control", is_zero_gdna),
)


def _rel(x: float) -> str:
    """``rel`` in 8 columns. ⛔ Switches to scientific below 1e-3 rather than rounding to ``0.000``:
    the S arm's residual against Fo is ~4e-5 and the whole point of that row is that it is small but
    NOT zero — a fixed 3-decimal format would have printed the verdict as an exact zero."""
    if not np.isfinite(x):
        return f"{'nan':>8}"
    return f"{x:>8.1e}" if 0.0 < abs(x) < 1e-3 else f"{x:>8.4f}"


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

    # ── ⛔ THE Fo JOIN, AND ITS TWO CHECKS ARE NOT EQUALLY STRONG.
    # The HARD one already ran per condition (`_oracle.check_walk_alignment`: the walk's record and
    # group counts against the scanner's own `stats.total` / `stats.n_read_names`) and raised if it
    # failed, so reaching here means the join is exact. What is printed is the SECONDARY diagnostic —
    # gDNA cannot splice — plus the number that says how much it is worth on this substrate: the
    # simulator writes each population in a block, so a 10 M-fragment condition has ~15 origin
    # transitions and a one-fragment slip would mislabel ~15 fragments. Loud for a big slip, blind to a
    # small one, which is exactly why it is not the gate.
    bad_align = [(r["condition"], r["overlap"]["spliced_gdna_units"], r["overlap"]["n_units"])
                 for r in rows if r["overlap"]["spliced_gdna_units"] > r["overlap"]["n_units"] // 1000]
    if bad_align:
        print(f"\n  ⛔⛔ gDNA IS SPLICED on {len(bad_align)} conditions — over 0.1 % of the "
              "gDNA-labelled units, which is impossible physics. The origin labels are wrong:")
        for c, n, tot in bad_align[:10]:
            print(f"       {c}  {n:,} spliced gdna units of {tot:,}")
        raise SystemExit(2)
    worst = max((r["overlap"]["spliced_gdna_units"] for r in rows), default=0)
    orphans = sum(r["overlap"]["orphan_units"] for r in rows)
    trans = sum(r["overlap"]["walk"]["n_transitions"] for r in rows)
    print(f"  ✅ frag_id join: record+group counts match the scanner on all {len(rows)} conditions "
          f"(hard gate). Secondary: at most {worst:,} of "
          f"{sum(r['overlap']['n_units'] for r in rows):,} units are spliced-and-gdna, "
          f"{orphans:,} claimed by no locus")
    print(f"     ⚠ {trans:,} origin transitions in BAM order across the panel — the secondary "
          "diagnostic is only sensitive to a slip larger than a population block.")

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
                  f"{s.abs_err:>14,.0f} {_rel(s.rel)} {s.net_err:>+14,.0f} "
                  f"{s.cancellation:>7.1f}")

    _RNA_TARGET = ("⭐ The reference is the UNSPLICED RNA units — the population rna_prior_count "
                   "actually targets. Fo carries the splice bit per unit, so this arm is no longer a "
                   "bound (the F_rna column was, and table ⑫ prices the difference).")

    arm_table("① P vs O — CALIBRATION'S OWN ERROR (a perfect deconvolution, same assembler) · gDNA arm",
              "P_vs_O", "gdna")
    arm_table("   … RNA arm", "P_vs_O", "rna")
    arm_table("④ O vs Fo — THE ASSEMBLER'S OWN ERROR (truth masses in) · gDNA arm",
              "O_vs_FO", "gdna",
              "⭐ Fo is the EM's own candidate count — every unit of the locus, labelled by true "
              "origin. This prices the mass→fragment-count conversion, the projection and the pooled "
              "share, alone.")
    arm_table("   … RNA arm", "O_vs_FO", "rna", _RNA_TARGET)
    arm_table("⑤ P vs Fo — THE TOTAL PRIOR ERROR (what the EM is handed, vs the true counts) · gDNA arm",
              "P_vs_FO", "gdna")
    arm_table("   … RNA arm", "P_vs_FO", "rna", _RNA_TARGET)

    # ── ⑧ the two halves of ④, separated ──
    arm_table("⑧ S vs Fo — THE ASSEMBLER WITH PERFECT PER-COMPONENT SHARES · gDNA arm",
              "S_vs_FO", "gdna",
              "⭐ Everything ④ measures EXCEPT the pooled share. The gap between ④ and ⑧ is what one "
              "share for two components costs.")
    arm_table("   … RNA arm", "S_vs_FO", "rna", _RNA_TARGET)
    arm_table("⑨ O vs S — THE POOLED SHARE'S OWN CONTRIBUTION, ISOLATED · gDNA arm",
              "O_vs_S", "gdna",
              "⛔ The locus TOTAL is conserved here to the fragment, so this error is purely "
              "compositional and no conservation gate can see it. Equal component lengths ⇒ identically "
              "zero, which is why the ladder is structurally blind to it.")

    # ── ⑪ what the yardstick correction is worth ──
    print()
    print("  ⑪ ⛔ THE YARDSTICK ITSELF — Fo (the EM's candidates) against F (first-base starts)")
    print("  ⭐ Fo − F is the STRADDLING population: fragments that overlap a locus but start outside "
          "it.")
    print("  ⛔ Every O−F / S−F number printed before 2026-08-08 was scored against F. The `rel` "
          "columns say what that cost.")
    print(f"    {'stratum':<26} {'Σ Fo':>13} {'Σ F':>13} {'Σ|Fo−F|':>11} {'rel':>8} "
          f"{'O−F rel':>8} {'O−Fo rel':>8} {'S−F rel':>8} {'S−Fo rel':>8}")
    print("    " + "-" * 112)
    for label, sel in _SELECTIONS:
        if label is None:
            print("    " + "-" * 112)
            continue
        sub = [r for r in rows if sel(r["condition"])]
        if not sub:
            continue
        y = _agg([ArmScore(**r["FO_vs_F"]["gdna"]) for r in sub])
        cells = {
            k: _agg([ArmScore(**r[k]["gdna"]) for r in sub])
            for k in ("O_vs_F", "O_vs_FO", "S_vs_F", "S_vs_FO")
        }
        print(f"    {label:<26} {y.total_arm:>13,.0f} {y.total_ref:>13,.0f} {y.abs_err:>11,.0f} "
              f"{_rel(y.rel)} {_rel(cells['O_vs_F'].rel)} {_rel(cells['O_vs_FO'].rel)} "
              f"{_rel(cells['S_vs_F'].rel)} {_rel(cells['S_vs_FO'].rel)}")

    # ── ⑫ the prior's POPULATION, and its STRENGTH ──
    print()
    print("  ⑫ ⭐ THE PRIOR DESCRIBES THE UNSPLICED POOL — is its claim right, and how strong is it?")
    print("  A spliced unit never gets a gDNA candidate (`em_solver.cpp`: has_gdna = !is_spliced && …),")
    print("  so the population `a_g : a_r` describes is gDNA units + UNSPLICED RNA units. ⛔ Scoring it")
    print("  against ALL RNA units reads a phantom +0.07…+0.10 tilt that is the denominator, not the")
    print("  prior (`TRAPS: score-the-consumers-own-count` — committed, then repeated, 2026-08-08).")
    print(f"    {'stratum':<26} {'pool gDNA':>13} {'pool RNA':>13} {'phi true':>9} {'phi S':>8} "
          f"{'Δ':>8} {'strength':>9} {'spliced RNA':>13}")
    print("    " + "-" * 116)
    for label, sel in _SELECTIONS:
        if label is None:
            print("    " + "-" * 116)
            continue
        sub = [r["overlap"] for r in rows if sel(r["condition"])]
        srows = [r for r in rows if sel(r["condition"])]
        if not sub:
            continue
        rg = sum(x["unit_totals"]["gdna"] for x in sub)
        rus = sum(x["rna_unspliced_total"] for x in sub)
        spliced = sum(x["unit_totals"]["rna"] for x in sub) - rus
        s_g = sum(ArmScore(**r["S_vs_FO"]["gdna"]).total_arm for r in srows)
        s_r = sum(ArmScore(**r["S_vs_FO"]["rna"]).total_arm for r in srows)
        pool = rg + rus
        phi_true = rg / pool if pool > 0 else float("nan")
        phi_s = s_g / (s_g + s_r) if (s_g + s_r) > 0 else float("nan")
        strength = (s_g + s_r) / pool if pool > 0 else float("nan")
        print(f"    {label:<26} {rg:>13,.0f} {rus:>13,.0f} {phi_true:>9.4f} "
              f"{phi_s:>8.4f} {phi_s - phi_true:>+8.4f} {strength:>9.3f} {spliced:>13,.0f}")
    print("    ⭐ `strength` is Σ(a_g+a_r) / the unspliced pool — pseudo-fragments per real fragment.")
    print("    It is ~1.000 BY CONSTRUCTION (the conserved count is that pool), so the posterior is a")
    print("    50/50 blend of calibration and the EM's own evidence and there is no knob. Not a defect;")
    print("    a design fact nothing had priced.")

    # ── ⑩ is gdna_eff_len clamped by an incidence sum? ──
    print()
    print("  ⑩ ⭐ gdna_eff_len's CLAMP — is `span` a genomic extent or an INCIDENCE sum?")
    print(f"    {'stratum':<26} {'support/genomic':>16} {'regions only':>12} {'Σ support':>16} "
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
              f"{float(np.median([x['median_region_over_genomic'] for x in sub_rows])):>12.2f} "
              f"{sum(x['total_support'] for x in sub_rows):>16,.0f} "
              f"{sum(x['total_genomic'] for x in sub_rows):>16,.0f}")
    print("    ⚠ `support/genomic` well above 1 means every interior boundary is adding ~mu_g − 1 to the "
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
            "overlap": {
                **r.overlap.diag,
                "rna_unspliced_total": float(r.overlap.rna_unspliced.sum()),
            },
            "seconds": r.seconds,
        }
        p, o = r.priors["P"], r.priors["O"]
        sa = r.priors["S"]
        fo_g, fo_r = r.overlap.gdna, r.overlap.rna_unspliced
        for ref_name, gref, rref, arm in (
            ("P_vs_O", o.gdna_prior_count, o.rna_prior_count, p),
            # ⭐ Fo is the reference every assembler arm is scored against — the EM's own candidate
            # count. F is kept beside it on the same arms so table ⑪ can price the correction.
            ("O_vs_FO", fo_g, fo_r, o),
            ("P_vs_FO", fo_g, fo_r, p),
            ("S_vs_FO", fo_g, fo_r, sa),
            ("O_vs_F", r.f_gdna, r.f_rna_upper, o),
            ("P_vs_F", r.f_gdna, r.f_rna_upper, p),
            ("S_vs_F", r.f_gdna, r.f_rna_upper, sa),
            ("O_vs_S", sa.gdna_prior_count, sa.rna_prior_count, o),
            # ⛔ the yardstick itself, as an arm: Fo scored against F
            ("FO_vs_F", r.f_gdna, r.f_rna_upper,
             PRIORS.LocusPriors(gdna_prior_count=fo_g, rna_prior_count=fo_r,
                                gdna_eff_len=np.zeros_like(fo_g))),
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
