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
contained counts (``density x effective_length``, not ``density x bp``). `SPECS` is an ordered ladder, each
rung adding one structure to the one before it. Magnitudes do not transfer between donors (direction is
preserved, size is not); a toy cannot rank defects; it is a correctness substrate, never a profiling one
(`TRAPS: toys-rank-hotspots-backwards`).
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


#: Ordered simplest-first, and each one adds exactly one structure to the one before it. That is the
#: point: when a row goes wrong, the thing that changed is the thing to look at.
SPECS: dict[str, ToySpec] = {
    # There is no gene-free rung: `TranscriptIndex` requires at least one transcript, so a
    # chromosome with no annotation cannot be indexed at all. A silent gene is the right first rung
    # anyway — it makes every object in the toy structurally pure gDNA, so any deviation from
    # f_g = 1 is a pure false positive with nothing to trade off against.
    "silent": ToySpec(
        name="silent",
        what_it_probes="⭐ ALL objects are pure gDNA (one silent gene) — every deviation from f_g = 1 "
        "is a false positive, with nothing to cancel against it",
        genome_length=60_000,
        genes=[_gene("g1", "+", [(20_000, 23_000), (28_000, 31_000)], 0.0)],
        n_rna_fragments=1,
    ),
    "TA_single_exon": ToySpec(
        name="TA_single_exon",
        what_it_probes="⭐⭐ OWNER'S SPEC. 5 kb chromosome, ONE single-exon transcript TA+ (1000,3000). "
        "NO introns and NO sj, so the exon can ONLY be solved through the two "
        "intergenic|exon BOUNDARIES: intergenic -> boundary -> exon -> boundary -> intergenic. It therefore "
        "tests exactly one thing — can an accurate intergenic gDNA level reach a single-stranded "
        "exon by message passing?",
        genome_length=5_000,
        genes=[_gene("TA", "+", [(1_000, 3_000)], 100.0, t_id="TA")],
        n_rna_fragments=1_000,
    ),
    "one_exon": ToySpec(
        name="one_exon",
        what_it_probes="a single-exon gene: one exon region, two intergenic|exon boundaries, no sj",
        genome_length=60_000,
        genes=[_gene("g1", "+", [(20_000, 23_000)], 400.0)],
    ),
    "two_exon": ToySpec(
        name="two_exon",
        what_it_probes="ONE intron between two exons — the intron|exon boundaries and the sj flux",
        genome_length=60_000,
        genes=[_gene("g1", "+", [(20_000, 23_000), (28_000, 31_000)], 400.0)],
    ),
    "spliced_exons": ToySpec(
        name="spliced_exons",
        what_it_probes="⭐⭐ OWNER'S SPEC. ONE two-exon transcript TA+ (1,000, 2,000) (9,000, 10,000) "
        "— so this is `nested_exons`'s TWIN at the same gene boundaries on the same 12 kb "
        "chromosome, with an INTRON and a SJ where the nesting was. FIVE REGIONS, FOUR "
        "contiguous BOUNDARIES and ⭐ the ladder's first SJ BOUNDARY:\n"
        "          REGION intergenic [0, 1000)        BOUNDARY @1,000   intergenic|exon, pure gDNA (TSS+)\n"
        "          REGION exon  [1000, 2000)   TA e1  BOUNDARY @2,000   intron|exon, the DONOR+ side\n"
        "          REGION intron [2000, 9000)  TA i1  BOUNDARY @9,000   intron|exon, the ACCEPTOR+ side\n"
        "          REGION exon  [9000, 10000)  TA e2  BOUNDARY @10,000  intergenic|exon, pure gDNA (TES+)\n"
        "          REGION intergenic [10000, 12000)\n"
        "          SJ BOUNDARY 2,000 → 9,000 (+), pure mature RNA, NOT a chain slot\n"
        "        ⭐⭐ What it adds over every rung before it, and why it is the hard one: the two "
        "exon↔intron BOUNDARIES. Mature RNA cannot cross an exon↔intron boundary contiguously, so their truth "
        "is pure gDNA — but the solver's own continuity gate says a strand IS admissible there "
        "(nascent RNA could cross), so they are NOT G1 and the solver must *derive* what the structure "
        "already implies. ⛔ On an `nrna_none` donor that is the maximally-violated case of the "
        "intron↔exon imputation premise, so a nascent rung is the control this one needs.\n"
        "        ⛔⛔ CAPTURE-ON NEEDS `--genome-length 120000` ON THIS RUNG. At 12 kb the whole "
        "chromosome gets ~39 gDNA fragments and the two exon↔intron BOUNDARIES carry 2 and 5 — not a solve, an "
        "empty chromosome. At 120 kb they carry 20 and 36 at the EXON's own capture stratum (density "
        "0.079–0.142 against the exon interior's 0.158–0.162 and the intron interior's 0.00015), which is "
        "the regime the object actually matters in. ⭐ That is capture working as intended: the gDNA "
        "signal LEAVES the intergenic and intronic REGIONS and arrives at the BOUNDARIES abutting the exon.\n"
        "        ⭐ And unlike `nested_exons` there IS own evidence inside the gene: the 7,000 bp intron "
        "REGION is where the intron factory lives, so the gDNA level does not have to travel from the "
        "gene ends. The two exons each sit between a G1 gene-boundary BOUNDARY and an exon|intron BOUNDARY, and "
        "the sj's flux is the only measurement of their mature RNA.",
        genome_length=12_000,
        genes=[_gene("g1", "+", [(1_000, 2_000), (9_000, 10_000)], 300.0, t_id="TA")],
    ),
    "alt_splice": ToySpec(
        name="alt_splice",
        what_it_probes="⭐⭐⭐ OWNER'S SPEC, 2026-08-05 — ALTERNATIVE SPLICING: several sj meeting "
        "at ONE boundary, which is the case a per-BOUNDARY sj total cannot represent.\n"
        "          TA+ (1,000, 2,000) (5,000, 6,000) (9,000, 10,000)   3 exons — the INCLUSION isoform\n"
        "          TB+ (1,000, 2,000) (9,000, 10,000)                  2 exons — the SKIPPING isoform\n"
        "        THREE sj over TWO shared sites:\n"
        "          j 2,000 -> 5,000   (TA's first intron)\n"
        "          j 6,000 -> 9,000   (TA's second intron)\n"
        "          j 2,000 -> 9,000   (TB's only intron — the exon-skipping jump)\n"
        "        ⭐⭐ SO THE SITES ARE SHARED, AND THAT IS THE POINT: the BOUNDARY @2,000 is the genomic-LOW "
        "end of TWO sj and the BOUNDARY @9,000 is the genomic-HIGH end of TWO. Both of @2,000's fluxes "
        "belong to its LOW flank and both of @9,000's to its HIGH flank, so each bank must POOL them as "
        "`Sum(count)/Sum(E)` — the ratio of sums, never the mean of ratios. ⛔ A single sj-inclusive "
        "total per BOUNDARY cannot express this at all, and neither can a per-sj rule that forgets the "
        "two share a line.\n"
        "        ⭐ It also adds a region that is exon AND intron on the SAME strand: [5,000, 6,000) is TA's "
        "middle exon and lies inside TB's intron. `splice_both_strands` had that contrast only ACROSS "
        "strands; here it is within one, so no strand bit can separate them and `coarse_type_array` calls "
        "it `exon`.\n"
        "        ⚠ What it does NOT cover: a BOUNDARY that is one sj's LOW end and another's HIGH end "
        "at once. That needs one transcript's intron to END where another's BEGINS.",
        genome_length=12_000,
        genes=[
            {
                "gene_id": "gA",
                "strand": "+",
                "transcripts": [
                    {
                        "t_id": "TA",
                        "exons": [(1_000, 2_000), (5_000, 6_000), (9_000, 10_000)],
                        "abundance": 300.0,
                    },
                    {"t_id": "TB", "exons": [(1_000, 2_000), (9_000, 10_000)], "abundance": 300.0},
                ],
            },
        ],
        n_rna_fragments=4_000,
    ),
    "tes_readthrough": ToySpec(
        name="tes_readthrough",
        what_it_probes="⭐⭐⭐ OWNER'S SPEC, 2026-08-05 — the CERTIFIED-RNA CHANNEL AT A TERMINUS BOUNDARY, "
        "which is the case no other rung can produce at all.\n"
        "          TA+ (1,050, 2,000) (9,000,  9,100)\n"
        "          TB+ (1,000, 2,000) (9,050, 11,000)\n"
        "        Two sj from ONE shared donor: j 2,000 -> 9,000 (TA) and j 2,000 -> 9,050 (TB).\n"
        "        ⭐⭐⭐ **BOUNDARY @9,100 IS THE POINT.** It is TA's TES and NO sj touches it — yet "
        "transcription CONTINUES past it, because TB's second exon runs to 11,000. A TB fragment that "
        "USED TB's sj and reaches >50 bp past 9,050 crosses 9,100 **contiguously having spliced "
        "elsewhere**, so it lands in `boundary_spliced` at a line with no sj to price it against. That "
        "is exactly the population the new TSS/TES boundaries create, and if it is not binned as spliced it "
        "falls into the UNSPLICED pool and gets deconvolved — certified RNA fed to the gDNA solver.\n"
        "        ⛔ **No previous toy can make this fragment.** On every earlier rung the exons ARE the "
        "regions, so a spliced molecule never crosses an interior line contiguously and `boundary_spliced` is "
        "structurally zero everywhere (measured 0 on `alt_splice`, including at exons holding 68,000 RNA "
        "fragments).\n"
        "        ⭐ The other three structures, each a separate stress:\n"
        "          BOUNDARY @9,050 — TB's sj ACCEPTOR **and** a plain contiguity line for TA, whose "
        "exon 2 spans 9,000-9,100 unbroken. So one line carries sj flux for one transcript and an "
        "unspliced RNA crossing for another.\n"
        "          BOUNDARY @1,050 — TA's TSS, with TB already transcribing through it.\n"
        "          REGION [9,000, 9,050) — TA exon AND TB intron on the SAME strand, 50 bp wide, so it is "
        "also below one fragment length and has no resolvable density of its own (TRAPS: density-below-one-fragment-length).\n"
        "        ⚠ Both `abundance` values are meant to be SWEPT: the certified channel's strength at "
        "@9,100 is TB's alone, while the unspliced crossing there is gDNA + TB, so the TA/TB ratio moves "
        "the two independently. A single abundance pair tests one corner of that.",
        genome_length=13_000,
        genes=[
            {
                "gene_id": "gA",
                "strand": "+",
                "transcripts": [
                    {"t_id": "TA", "exons": [(1_050, 2_000), (9_000, 9_100)], "abundance": 300.0},
                    {"t_id": "TB", "exons": [(1_000, 2_000), (9_050, 11_000)], "abundance": 300.0},
                ],
            },
        ],
        n_rna_fragments=4_000,
    ),
    "splice_both_strands": ToySpec(
        name="splice_both_strands",
        what_it_probes="⭐⭐⭐ OWNER'S SPEC, 2026-08-05 — the rung the SPLICE-FLUX REFRAME logic must be "
        "derived against. FOUR transcripts, BOTH strands, overlapping exons AND overlapping introns, "
        "and TWO splice junctions pointing opposite ways:\n"
        "          TA+ (2,000, 3,000) (9,000, 10,000)      2 exons, + strand, intron 3,000-9,000\n"
        "          TB+ (2,000, 10,000)                     1 exon,  + strand, spans TA's intron\n"
        "          TC− (1,000, 11,000)                     1 exon,  − strand, spans everything\n"
        "          TD− (1,000, 2,500) (8,500, 11,000)      2 exons, − strand, intron 2,500-8,500\n"
        "        ⭐⭐ WHY THIS ONE. Every previous rung let a BOUNDARY answer 'is my neighbour an exon?' "
        "with a yes or a no. Here it cannot: THREE regions are simultaneously an INTRON on one strand and "
        "an EXON on the other — [2,500, 3,000), [3,000, 8,500) and [8,500, 9,000) — so 'exon' is not a "
        "property of a region at all, it is a property of (region, strand). And the two sj are on "
        "OPPOSITE strands, so a boundary can be the DONOR of one and sit beside the ACCEPTOR of the other.\n"
        "        ⛔ The question it exists to answer is per (BOUNDARY, side, strand, donor-or-acceptor, "
        "message direction): when this boundary reframes against that neighbour, does its splice flux belong "
        "in the total or not? The derivation is open.\n"
        "        ⚠ MY READING OF THE OWNER'S SPEC, flagged rather than assumed: the owner wrote the last "
        "two transcripts both as `TC-`. Two transcripts cannot share an id, so they are TC− and TD− "
        "here. If the intent was one transcript with two isoforms the ids change and nothing else does.",
        genome_length=12_000,
        genes=[
            {
                "gene_id": "gP",
                "strand": "+",
                "transcripts": [
                    {"t_id": "TA", "exons": [(2_000, 3_000), (9_000, 10_000)], "abundance": 300.0},
                    {"t_id": "TB", "exons": [(2_000, 10_000)], "abundance": 300.0},
                ],
            },
            {
                "gene_id": "gM",
                "strand": "-",
                "transcripts": [
                    {"t_id": "TC", "exons": [(1_000, 11_000)], "abundance": 300.0},
                    {"t_id": "TD", "exons": [(1_000, 2_500), (8_500, 11_000)], "abundance": 300.0},
                ],
            },
        ],
        n_rna_fragments=4_000,
    ),
    "nested_exons": ToySpec(
        name="nested_exons",
        what_it_probes="⭐⭐ NESTED EXONS (owner's spec). THREE nested single-exon transcripts on one gene:\n"
        "                     TA+ (1,000, 10,000)   TB+ (2,000, 9,000)   TC+ (3,000, 8,000)\n"
        "        The partition is therefore SEVEN REGIONS and SIX BOUNDARIES, and there is NO INTRON "
        "anywhere:\n"
        "          REGION intergenic [0, 1000)        BOUNDARY @1,000   intergenic|exon, pure gDNA\n"
        "          REGION exon [1000, 2000)   TA      BOUNDARY @2,000   exon|exon\n"
        "          REGION exon [2000, 3000)   TA+TB   BOUNDARY @3,000   exon|exon\n"
        "          REGION exon [3000, 8000)   TA+B+C  BOUNDARY @8,000   exon|exon\n"
        "          REGION exon [8000, 9000)   TA+TB   BOUNDARY @9,000   exon|exon\n"
        "          REGION exon [9000, 10000)  TA      BOUNDARY @10,000  intergenic|exon, pure gDNA\n"
        "          REGION intergenic [10000, 12000)\n"
        "        ⭐ Every one of the five exon REGIONS carries RNA and the library is unstranded, so NONE "
        "of them has composition evidence of its own; and with no intron there is no object in the "
        "middle of the gene that can re-derive the gDNA level for itself. So the gDNA level has to "
        "travel from the two pure-gDNA BOUNDARIES at the gene ends through five evidence-free objects. "
        "⭐⭐ Because the transcripts are NESTED, the RNA density is a symmetric staircase (1, 2, 3, 2, 1 "
        "transcripts) while the gDNA density is UNIFORM — so the truth is known per object and any "
        "systematic drift in the delivered gDNA level is visible against it.",
        genome_length=12_000,
        genes=[
            {
                "gene_id": "g1",
                "strand": "+",
                "transcripts": [
                    {"t_id": "TA", "exons": [(1_000, 10_000)], "abundance": 300.0},
                    {"t_id": "TB", "exons": [(2_000, 9_000)], "abundance": 300.0},
                    {"t_id": "TC", "exons": [(3_000, 8_000)], "abundance": 300.0},
                ],
            }
        ],
        n_rna_fragments=2_000,
    ),
    "nested_exons_neg": ToySpec(
        name="nested_exons_neg",
        what_it_probes="⭐ THE STRAND MIRROR of `nested_exons` — the same three nested transcripts on the "
        "MINUS strand. A − transcript runs right-to-left, so its TSS is at its HIGH coordinate and its "
        "TES at its low one: the FLAG_TSS_NEG / FLAG_TES_NEG bits must therefore implicate the OPPOSITE "
        "flank from their + counterparts. ⛔ That is a convention, not a derivation, and a convention "
        "that is assumed rather than pinned is how a sign gets flipped silently — so this rung exists to "
        "pin it against oracle truth.",
        genome_length=12_000,
        genes=[
            {
                "gene_id": "g1",
                "strand": "-",
                "transcripts": [
                    {"t_id": "TA", "exons": [(1_000, 10_000)], "abundance": 300.0},
                    {"t_id": "TB", "exons": [(2_000, 9_000)], "abundance": 300.0},
                    {"t_id": "TC", "exons": [(3_000, 8_000)], "abundance": 300.0},
                ],
            }
        ],
        n_rna_fragments=2_000,
    ),
    "deep_exon": ToySpec(
        name="deep_exon",
        what_it_probes="⭐ TRAPS: a-purity-filter-is-a-length-filter: a LARGE, RNA-rich exon beside an intron",
        genome_length=80_000,
        genes=[_gene("g1", "+", [(20_000, 34_000), (40_000, 43_000)], 3000.0)],
        n_rna_fragments=120_000,
    ),
    "nascent": ToySpec(
        name="nascent",
        what_it_probes="two_exon plus NASCENT RNA — the only way an intron|exon boundary legitimately "
        "carries RNA, which the whole ladder panel cannot express",
        genome_length=60_000,
        genes=[_gene("g1", "+", [(20_000, 23_000), (28_000, 31_000)], 400.0)],
        nrna_abundance=40.0,
    ),
    "two_genes": ToySpec(
        name="two_genes",
        what_it_probes="one expressed gene and one SILENT gene — the silent one's exons are pure gDNA",
        genome_length=90_000,
        genes=[
            _gene("g1", "+", [(20_000, 23_000), (28_000, 31_000)], 400.0),
            _gene("g2", "-", [(50_000, 53_000), (58_000, 61_000)], 0.0),
        ],
    ),
    "opposite_strands": ToySpec(
        name="opposite_strands",
        what_it_probes="two genes on OPPOSITE strands, well separated — the strand channel's control",
        genome_length=90_000,
        genes=[
            _gene("g1", "+", [(20_000, 23_000), (28_000, 31_000)], 400.0),
            _gene("g2", "-", [(50_000, 53_000), (58_000, 61_000)], 400.0),
        ],
    ),
    "encompassing": ToySpec(
        name="encompassing",
        what_it_probes="⭐⭐ OWNER'S SPEC, 2026-09-13 — a single-exon TB− (10,000–30,000) ENCOMPASSING a two-exon "
        "TA+ (11,000–12,000, 19,000–20,000): every slot from 11,000 to 20,000 is both-stranded, TA+'s exons "
        "are exons of both genes and its intron is TB−'s exon, and TB−'s level — measured in its single-strand "
        "flanks — must cross every one of TA+'s boundaries (they carry only TA+'s bits) while TA+'s junction "
        "flux is the + source at its exons although TA+ has no single-strand exon anywhere. Gated at four "
        "abundance regimes by `tests/calibration/test_encompassing_locus.py`; this rung is TA+ ≫ TB−.",
        genome_length=40_000,
        genes=[
            _gene("gA", "+", [(11_000, 12_000), (19_000, 20_000)], 1000.0),
            _gene("gB", "-", [(10_000, 30_000)], 30.0),
        ],
        n_rna_fragments=20_000,
    ),
}
