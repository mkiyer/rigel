"""The toy harness — a mini chromosome a test defines, calibrated in seconds with every object's answer beside
its truth. A TEST SUBSTRATE, loaded by ``test_toy_harness.py`` (its gates) and ``test_encompassing_locus.py``.

A toy spec is a handful of genes with exactly the structure under interrogation, simulated, scanned,
drained by the second pass exactly as production runs it, split by origin, and calibrated, with the
per-object result read beside the per-object truth from the simulator's read names. The library-level
quantities a toy cannot fit for itself — the fitted fragment-length pmfs, the strand balance with its
overdispersions and noise-floor sample sizes, the enrichment landscape, the intron background, the capture
knobs and the read-simulation settings — are harvested from one real cached condition, the donor, and
injected through `InjectedCalibrationPriors`; the toy supplies only the controlled geometry. Depth matching
is not optional: the landscape is an absolute log-density model, so the toy's gDNA count is derived from the
donor's own gDNA density per base (measured on its intergenic slots) times the chromosome length, never
chosen. There is no gDNA knob; the RNA side is the experimental variable. A region's stored counts are
contained counts (``density x effective_length``, not ``density x bp``). Magnitudes do not transfer
between donors (direction is preserved, size is not); a toy cannot rank defects; it is a correctness
substrate, never a profiling one (`TRAPS: toys-rank-hotspots-backwards`).
"""

from __future__ import annotations

import dataclasses
import json
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _oracle import OracleTruth  # noqa: E402

from rigel.calibration.calibrate import InjectedCalibrationPriors, calibrate  # noqa: E402
from rigel.calibration.region_chain import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import (  # noqa: E402
    _drain_side_buffer,
    _native_detect_sj_tag,
    scan_and_buffer,
)
from rigel.scan_cache import index_derived_inputs  # noqa: E402
from rigel.sim import CaptureConfig, GDNAConfig, ReadSimConfig, Scenario  # noqa: E402

#: Coarse region type names, shared with every other instrument (`signature.coarse_type_array`).
TYPE_NAMES = {0: "intergenic", 1: "intron", 2: "exon"}


# ──────────────────────────────────────────────────────────────────────────────────────────────────
# THE DONOR — everything the toy cannot fit for itself
# ──────────────────────────────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class DonorGlobals:
    """The library-level conditions harvested from one real cached condition.

    Deliberately not cached to disk. The bundle is a function of the calibration code that fit it, so
    a stored copy goes stale on exactly the changes this harness exists to test, and a stale global is
    invisible — it does not crash, it just quietly answers a different question. Harvesting costs one
    scan plus one calibrate; harvest once per session and run many toys against it.
    """

    condition: str
    priors: InjectedCalibrationPriors
    gdna_fl_pmf: np.ndarray
    rna_fl_pmf: np.ndarray
    #: gDNA molecules per base, measured on the donor's own structurally-pure-gDNA population.
    gdna_rate_per_base: float
    #: the donor's read-simulation settings, so the toy draws from the same library
    frag_mean: float
    frag_std: float
    frag_min: int
    frag_max: int
    read_length: int
    strand_specificity: float
    #: capture: whether the donor had it on, and the numeric knobs to reproduce it on toy probes
    capture_on: bool
    capture_knobs: dict = field(default_factory=dict)

    def describe(self) -> str:
        return (
            f"donor={self.condition}\n"
            f"  kappa={self.priors.rna_sense_frac:.6f}  n_rna_obs={self.priors.n_rna_obs:,.0f}\n"
            f"  od_rna={self.priors.rna_strand_overdispersion:.4g}  "
            f"od_gdna={self.priors.gdna_strand_overdispersion:.4g}\n"
            f"  gDNA rate = {self.gdna_rate_per_base:.6g} molecules/base   "
            f"FL {self.frag_mean:.0f}+-{self.frag_std:.0f} [{self.frag_min},{self.frag_max}] "
            f"read {self.read_length}\n"
            f"  strand_specificity={self.strand_specificity:.4g}  capture={'ON' if self.capture_on else 'off'}\n"
            f"  abundance_landscape="
            f"{'yes' if self.priors.abundance_landscape is not None else 'NONE'}  "
            # A `describe()` string is executed by no test, so a field deleted from
            # `InjectedCalibrationPriors` leaves this line raising while the suite stays green
            # (`TRAPS: a-green-suite-hid-five-dead-instruments`). Grep this function when a field
            # leaves that class.
            f"intron_background={'yes' if self.priors.intron_background is not None else 'NONE'}"
        )


def harvest(
    donor_dir: Path,
    index: TranscriptIndex,
    *,
    config: CalibrationConfig | None = None,
    pipeline_config: PipelineConfig | None = None,
    bam: str | None = None,
    name: str | None = None,
) -> DonorGlobals:
    """Calibrate one cached condition and keep everything the toy cannot fit for itself.

    ``bam`` / ``name`` override the defaults (``donor_dir/sim_oracle.bam`` and ``donor_dir.name``) so a
    gate can harvest from a scenario it built itself rather than from a 10 M-fragment panel condition
    that only exists on one machine.
    """
    config = config or CalibrationConfig()
    pipeline_config = pipeline_config or PipelineConfig()
    bam = bam or str(donor_dir / "sim_oracle.bam")
    scan = dataclasses.replace(pipeline_config.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _stats, strand_model, _buf, payload = scan_and_buffer(bam, index, scan)

    ra = RegionArrays.from_index(index)
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.sj_opportunity import crossing_probability_from_index

    max_size = int(payload.max_length)
    fl = build_fl_models(
        payload,
        sj_opportunity=crossing_probability_from_index(index, max_size),
        gdna_opportunity=gdna_opportunity_from_index(index, max_size),
    )
    debug: dict = {}
    calibrate(
        payload=payload,
        strand_model=strand_model,
        gdna_fl_pmf=fl.gdna_pmf,
        rna_fl_pmf=fl.rna_pmf,
        config=config,
        _debug=debug,
        **index_derived_inputs(index),
    )

    # the gDNA rate, from the pure population, in the SAME frame the toy will be measured in
    rate = _rate_from_capture(debug["capture"], debug["chain"], ra)

    name = name or donor_dir.name
    sim = _donor_sim_params(donor_dir, name)
    return DonorGlobals(
        condition=name,
        priors=debug["calibration_priors"],
        gdna_fl_pmf=fl.gdna_pmf,
        rna_fl_pmf=fl.rna_pmf,
        gdna_rate_per_base=rate,
        capture_on="capture_on" in name,
        **sim,
    )


def _rate_from_capture(capture, chain, region_arrays) -> float:
    """``sum(count) / sum(eff_gdna)`` over the donor's INTERGENIC region slots.

    Both arrays are the solver's own (``capture.count`` and ``capture.eff_gdna``), so the rate
    is in exactly the frame the toy's own regions will be measured in — no second implementation of an
    effective length (`TRAPS: two-docstrings-one-quantity`).
    """
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    is_region = kind == REGION
    pure = is_region.copy()
    pure[is_region] = rtype[obj[is_region]] == 0
    count = np.asarray(capture.count, np.float64).sum(axis=1)
    eff = np.asarray(capture.eff_gdna, np.float64)
    if not pure.any() or eff[pure].sum() <= 0.0:
        raise ValueError("no intergenic region slots with gDNA opportunity in the donor")
    return float(count[pure].sum() / eff[pure].sum())


def _donor_sim_params(donor_dir: Path, name: str) -> dict:
    """The donor's read-simulation settings.

    Fragment lengths come from ``truth_summary.json``'s post-capture measurement where it exists,
    never from a configured ``frag_mean`` — capture selects for length, so the configured parameters
    describe a library that was never sequenced (`TRAPS: capture-selects-for-length`). Strand
    specificity is read off the condition name, which is where the panel encodes it.
    """
    ss = 0.5
    for part in name.split("_"):
        try:
            v = float(part)
        except ValueError:
            continue
        if 0.0 <= v <= 1.0:
            ss = v
            break
    summary = donor_dir / "truth_summary.json"
    frag_mean, frag_std, frag_min, frag_max = 206.0, 98.0, 50, 500
    if summary.is_file():
        d = json.loads(summary.read_text())
        allrow = d.get("fragment_lengths", {}).get("all", {})
        if allrow.get("mean") is not None:
            frag_mean = float(allrow["mean"])
            frag_std = float(allrow["std"])
            frag_min = int(allrow["min"])
            frag_max = int(allrow["max"])
    return {
        "frag_mean": frag_mean,
        "frag_std": frag_std,
        "frag_min": frag_min,
        "frag_max": frag_max,
        "read_length": 100,
        "strand_specificity": ss,
        "capture_knobs": {
            "off_target_weight": 1.0,
            "binding_per_base": 10.0,
            "min_overlap": 1,
            "probe_length": 120,
            "capture_fraction": 1.0,
        },
    }


# ──────────────────────────────────────────────────────────────────────────────────────────────────
# THE TOY
# ──────────────────────────────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class ToySpec:
    """A mini chromosome: what to put on it, and how much RNA to express.

    There is no gDNA knob. The gDNA level is not a free parameter — it is set to match the donor's
    measured density per base, because the injected enrichment landscape is absolute. The RNA side is
    the experimental knob, and ``abundance`` on each transcript is what you vary.
    """

    name: str
    what_it_probes: str
    genome_length: int
    genes: list[dict]
    n_rna_fragments: int = 40_000
    nrna_abundance: float = 0.0
    seed: int = 7
    #: transcript ids to put a probe on when the donor is capture-ON. ``None`` = all of them.
    captured: tuple[str, ...] | None = None


@dataclass
class ToyResult:
    spec: ToySpec
    donor: DonorGlobals
    result: object  #: the CalibrationResult
    truth: OracleTruth
    payload: object
    region_arrays: RegionArrays
    chain: object
    capture: dict
    n_gdna_target: int
    seconds: float
    #: the toy's own index, so a caller can rebuild anything index-derived — in particular
    #: `splice_graph.build_boundary_flags_array`, the TSS/TES/DONOR/ACCEPTOR bits per boundary.
    index: object = None


def _toy_probes(spec: "ToySpec", out: Path, knobs: dict) -> str:
    """Probes over the toy's transcripts, written from the spec's own exon coordinates.

    Deliberately not ``write_random_capture_probes``: that draws a random subset of genes, which is
    right for a panel and wrong for a toy, where you decide what is captured and a random draw would
    make the condition depend on a seed. Every transcript named in ``spec.captured`` (default: all of
    them) is tiled with probes per exon, which makes probe placement part of the controlled geometry.

    The file is the transcript-coordinate TSV the sampler already reads: ``transcript_id start end``,
    coordinates in transcript space (`sampler._load_transcript_probes`).
    """
    plen = int(knobs.get("probe_length", 120))
    lines = ["transcript_id\tstart\tend"]
    for gene in spec.genes:
        for t in gene["transcripts"]:
            if spec.captured is not None and t["t_id"] not in spec.captured:
                continue
            # Tiled, not one centred probe: a single central probe leaves a transcript's ends
            # uncovered, so its `intergenic|exon` boundaries stay at off-target density, and a 0-bp
            # line's counts are `density x mean_FL` however long the chromosome is. The donor panels
            # are tiled (`design_suite_probes.py` at probe_density 1.0), which is what enriches a
            # first/last exon's boundary in the first place. Match that.
            #
            # And tiled per exon, so every probe abuts the intron|exon boundaries and none straddles a
            # sj. Probes are written in transcript space, so a probe spanning an internal sj offset has
            # a genomic footprint in two blocks, and a gDNA fragment then binds only the better of the
            # two parts — the population that spans an intron|exon boundary loses the other half.
            # Per-exon tiling leaves every probe inside one exon, unsplit, ending exactly on the
            # boundary, so a boundary-crossing fragment takes the full binding weight for its exon-side
            # overlap. It is also the honest geometry: a
            # probe boundary at an exon end is what a real panel produces, and the split-probe case is
            # a separate population worth its own rung.
            off = 0
            for s, e in t["exons"]:
                elen = e - s
                for start in range(0, elen, plen):
                    lines.append(f"{t['t_id']}\t{off + start}\t{off + min(start + plen, elen)}")
                off += elen
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n")
    return str(out)


def run_toy(
    spec: ToySpec,
    donor: DonorGlobals,
    work_dir: Path,
    *,
    config: CalibrationConfig | None = None,
    pipeline_config: PipelineConfig | None = None,
) -> ToyResult:
    """Simulate → scan → split by origin → calibrate with the donor's globals injected."""
    import time

    t0 = time.perf_counter()
    config = config or CalibrationConfig()
    pipeline_config = pipeline_config or PipelineConfig()
    wd = Path(work_dir) / spec.name
    wd.mkdir(parents=True, exist_ok=True)

    # ── the gDNA level is DERIVED, not chosen: match the donor's molecules per base ──────────────
    n_gdna = int(round(donor.gdna_rate_per_base * spec.genome_length))
    gdna_fraction = n_gdna / max(spec.n_rna_fragments, 1)

    sim_cfg = ReadSimConfig(
        frag_mean=int(round(donor.frag_mean)),
        frag_std=int(round(donor.frag_std)),
        frag_min=donor.frag_min,
        frag_max=donor.frag_max,
        read_length=donor.read_length,
        strand_specificity=donor.strand_specificity,
        seed=spec.seed,
    )
    gdna_cfg = GDNAConfig(
        abundance=0.0,
        frag_mean=int(round(donor.frag_mean)),
        frag_std=int(round(donor.frag_std)),
    )

    sc = Scenario(
        spec.name,
        genome_length=spec.genome_length,
        seed=spec.seed,
        work_dir=wd / "sim",
    )
    for gene in spec.genes:
        sc.add_gene(gene["gene_id"], gene["strand"], gene["transcripts"])

    capture_cfg = None
    if donor.capture_on:
        probes = _toy_probes(spec, wd / "probes.tsv", donor.capture_knobs)
        k = donor.capture_knobs
        capture_cfg = CaptureConfig(
            probes=probes,
            probe_format="transcript",
            off_target_weight=float(k["off_target_weight"]),
            binding_per_base=float(k["binding_per_base"]),
            min_overlap=int(k["min_overlap"]),
        )

    res = sc.build_oracle(
        n_rna_fragments=int(spec.n_rna_fragments),
        gdna_fraction=gdna_fraction,
        nrna_abundance=float(spec.nrna_abundance),
        sim_config=sim_cfg,
        gdna_config=gdna_cfg,
        capture_config=capture_cfg,
    )

    bam = str(res.bam_path)
    scan = dataclasses.replace(pipeline_config.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _stats, strand_model, _buf, pass_one = scan_and_buffer(bam, res.index, scan)
    ra = RegionArrays.from_index(res.index)

    # ── the second pass, on the whole, exactly as production runs it ─────────────────────────────
    # An undrained tally understates the spliced population: a held fragment is held precisely because
    # its unsequenced gap admits more than one intron path. `_lift` is what makes the oracle valid
    # afterwards: score and draw once on the whole, then replay each fragment's chosen hypothesis
    # inside whichever origin partition holds it (`TRAPS: draining-breaks-the-oracle`).
    lift: dict = {}
    payload = _drain_side_buffer(
        pass_one, res.index, strand_model, seed=pipeline_config.second_pass_seed, _lift=lift
    )
    truth = OracleTruth.from_bam(
        bam,
        res.index,
        pipeline_config,
        wd / "split",
        spec.name,
        # Drained whole here (it is what calibration reads and what sum-to-full is asserted against);
        # undrained whole inside `drain_with` (the drained bank holds nothing, so it has no key pool).
        full_payload=payload,
        drain_with=(
            (lift["undrained"], lift["choices"], lift["region_types"], lift["sj"]) if lift else None
        ),
    )

    debug: dict = {}
    out = calibrate(
        payload=payload,
        strand_model=strand_model,
        gdna_fl_pmf=donor.gdna_fl_pmf,
        rna_fl_pmf=donor.rna_fl_pmf,
        config=config,
        injected_priors=donor.priors,
        _debug=debug,
        **index_derived_inputs(res.index),
    )
    return ToyResult(
        spec=spec,
        donor=donor,
        result=out,
        truth=truth,
        payload=payload,
        region_arrays=ra,
        chain=debug["chain"],
        capture=debug["capture"],
        n_gdna_target=n_gdna,
        seconds=time.perf_counter() - t0,
        index=res.index,
    )


# ──────────────────────────────────────────────────────────────────────────────────────────────────
# THE REPORT — every object, because there are few enough to read
# ──────────────────────────────────────────────────────────────────────────────────────────────────


def object_rows(r: ToyResult) -> list[dict]:
    """One row per region and per contiguous boundary, in genomic order along the chain."""
    chain, ra = r.chain, r.region_arrays
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    start = np.asarray(ra.start, np.int64)
    size = np.asarray(ra.region_size_bp, np.float64)

    cap = r.capture
    fg_loc = np.asarray(cap.fg_loc, np.float64)
    fg = np.asarray(cap.f_g, np.float64)
    var_g = np.asarray(cap.var_g, np.float64)
    tau = np.asarray(cap.tau_lam, np.float64)
    count = np.asarray(cap.count, np.float64).sum(axis=1)
    mature = np.asarray(cap.mature, np.float64)
    spliced = np.asarray(cap.spliced, np.float64)

    ov = r.truth.override_masses(ra)
    tg = {
        "region": np.asarray(ov["count_gdna_region"], np.float64),
        "boundary": np.asarray(ov["count_gdna_boundary"], np.float64),
    }
    tr = {
        "region": np.asarray(ov["count_rna_region"], np.float64),
        "boundary": np.asarray(ov["count_rna_boundary"], np.float64),
    }
    pg = {
        "region": np.asarray(r.result.count_gdna_region, np.float64),
        "boundary": np.asarray(r.result.count_gdna_boundary, np.float64),
    }
    pr = {
        "region": np.asarray(r.result.count_rna_region, np.float64),
        "boundary": np.asarray(r.result.count_rna_boundary, np.float64),
    }

    rows = []
    for s in range(int(chain.n_slots)):
        i = int(obj[s])
        axis = "region" if kind[s] == REGION else "boundary"
        if axis == "region":
            label = TYPE_NAMES[int(rtype[i])]
            where = f"{int(start[i]):,}–{int(start[i] + size[i]):,}"
            bp = int(size[i])
        else:
            lo, hi = s - 1, s + 1
            a = int(rtype[obj[lo]]) if lo >= 0 and kind[lo] == REGION else -1
            b = int(rtype[obj[hi]]) if hi < int(chain.n_slots) and kind[hi] == REGION else -1
            pair = "|".join(TYPE_NAMES.get(x, "?") for x in sorted((a, b)) if x >= 0)
            label = pair or "boundary"
            where = f"@{int(start[obj[hi]]):,}" if b >= 0 else f"#{i}"
            bp = 0
        t_tot = tg[axis][i] + tr[axis][i]
        p_tot = pg[axis][i] + pr[axis][i]
        rows.append(
            {
                "slot": s,
                "axis": axis,
                "type": label,
                "where": where,
                "bp": bp,
                "n": float(count[s]),
                "spliced": float(spliced[s]),
                "sj": float(mature[s]),
                "true_fg": float(tg[axis][i] / t_tot) if t_tot > 0 else float("nan"),
                "fg_loc": float(fg_loc[s]),
                "pred_fg": float(fg[s]),
                "sd_fg": float(np.sqrt(max(var_g[s], 0.0))),
                "tau": float(tau[s]),
                "err": float(abs(pg[axis][i] - tg[axis][i])),
                "mass": float(t_tot if t_tot > 0 else p_tot),
            }
        )
    return rows


def _gene(gid, strand, exons, abundance, t_id=None):
    return {
        "gene_id": gid,
        "strand": strand,
        "transcripts": [{"t_id": t_id or f"{gid}_t1", "exons": exons, "abundance": abundance}],
    }
