"""The whole-genome condition grid behind `simulate_reads.py`: :func:`run_simulation` and every step it takes. A helper
module, not a script.

:func:`run_simulation` loads the transcript set from the rigel index (``rigel.sim.transcriptome``: annotated
transcripts plus the synthetic nascent entities, ``t_index`` equal to the index row), appends any shadow transcripts,
assigns the base abundances once (``_sim_abundance``) and runs the grid (:func:`run_condition_grid`): nascent x gDNA
rate x gDNA strand overdispersion x strand specificity x capture. It writes ``manifest.json`` last.

Per condition the grid copies the transcripts at their base abundances, applies the nascent mode, splits the fragment
budget (:func:`resolve_depths`), simulates on the capture sampler its panel shares, writes the truth (``_sim_truth``)
and returns a manifest entry. A condition whose output already exists is skipped. The seed keys on ``(gdna, strand
specificity)`` only (:func:`capture_paired_condition_seed`), so the capture, overdispersion and nascent variants of
one base condition are paired.
"""

from __future__ import annotations

import copy
import hashlib
import json
import logging
import time
from dataclasses import asdict, dataclass, is_dataclass, replace
from itertools import product
from pathlib import Path
from typing import Any

from rigel.index import TranscriptIndex, load_manifest, load_reference_lengths, source_record
from rigel.sim.capture import CaptureSampler
from rigel.sim.transcriptome import load_transcripts_from_index
from rigel.sim.wgs_config import GDNASimConfig
from rigel.sim.wgs_engine import WholeGenomeSimulator
from rigel.transcript import Transcript

from _sim_abundance import apply_nrna_ratio, apply_sparse_nrna, assign_base_abundances
from _sim_config import CaptureScenario, NRNAConfig, WholeGenomeSimConfig, refuse_a_run_that_keeps_no_reads
from _sim_truth import write_molar_truth, write_post_capture_truth

logger = logging.getLogger(__name__)


# ── the transcript set ───────────────────────────────────────────────────────────────────────────


def load_transcripts(gtf_path: str | Path) -> list[Transcript]:
    """The transcripts of a GTF (may be gzipped), ``t_index`` assigned sequentially."""
    logger.info("Loading transcripts from %s", gtf_path)
    transcripts = Transcript.read_gtf(str(gtf_path), parse_mode="warn-skip")
    logger.info("Read %d transcripts from GTF", len(transcripts))
    for i, t in enumerate(transcripts):
        t.t_index = i
    logger.info(
        "Final: %d transcripts from %d genes", len(transcripts), len({t.g_id for t in transcripts})
    )
    return transcripts


def merge_shadow_transcripts(
    transcripts: list[Transcript], shadow_gtf: str | Path
) -> list[Transcript]:
    """Append shadow transcripts, unannotated transcription the simulator draws from and the index never sees.

    The shadow GTF's rows follow the index's (``transcripts``) in ``t_index``, with ``is_synthetic = False``,
    ``span_t_index = -1`` and no nascent. A shadow whose ``t_id`` the index knows is refused, since the tool could
    know it. Their fragments are named like any RNA fragment, so the truth shows RNA where the annotation has none.
    """
    shadows = load_transcripts(shadow_gtf)
    known = {t.t_id for t in transcripts}
    clash = sorted(t.t_id for t in shadows if t.t_id in known)
    if clash:
        raise ValueError(
            f"{len(clash)} shadow transcript(s) are KNOWN to the index and are therefore not shadows: "
            f"{clash[:5]}{'…' if len(clash) > 5 else ''} — remove them from the shadow GTF or the annotation"
        )
    base = max((t.t_index for t in transcripts), default=-1) + 1
    for k, t in enumerate(shadows):
        t.t_index = base + k
        t.is_synthetic = False
        t.span_t_index = -1
        t.nrna_abundance = 0.0
    logger.info(
        "Shadow transcripts: %d appended from %s (unknown to the index by construction)",
        len(shadows),
        shadow_gtf,
    )
    return transcripts + shadows


def ensure_index(cfg: WholeGenomeSimConfig) -> Path:
    """The index the transcriptome comes from: ``cfg.index`` if set, else one built from ``genome`` + ``gtf`` into
    ``<outdir>/rigel_index`` (duplicates collapsed) and reused on the next run."""
    if cfg.index:
        index_dir = Path(cfg.index)
        if not index_dir.is_dir():
            raise FileNotFoundError(f"index not found: {index_dir}")
        return index_dir
    index_dir = Path(cfg.outdir) / "rigel_index"
    if not index_dir.is_dir():
        logger.info(
            "No index configured — building one from %s + %s into %s", cfg.genome, cfg.gtf, index_dir
        )
        TranscriptIndex.build(
            cfg.genome, cfg.gtf, index_dir, write_tsv=False, collapse_duplicate_transcripts=True
        )
    return index_dir


def _index_provenance(index_dir: Path, gtf: str | Path) -> dict:
    """The index record ``manifest.json`` carries: its resolved path and its own manifest verbatim.

    Refuses a ``gtf`` whose sha256 is not the one the index recorded, so a config cannot name one annotation and
    simulate another; a byte-identical copy elsewhere passes. Gate: ``tests/test_whole_genome_sim_config.py``."""
    index_manifest = load_manifest(index_dir)
    built_from = index_manifest["sources"]["gtf"]
    given = source_record(gtf)
    if given["sha256"] != built_from["sha256"]:
        raise ValueError(
            f"gtf {given['path']} is not the GTF the index {Path(index_dir).resolve()} was built from "
            f"({built_from['path']}): sha256 {given['sha256']} against {built_from['sha256']}. The "
            "simulated transcriptome is the index's, so gtf: must name that GTF; or build the index "
            "from this one."
        )
    return {"path": str(Path(index_dir).resolve()), "manifest": index_manifest}


# ── the axes ─────────────────────────────────────────────────────────────────────────────────────


def gdna_label_for_rate(rate: float, labels: list[str] | None, index: int) -> str:
    """Return a stable gDNA label for a configured rate."""
    if labels is not None and index < len(labels):
        return labels[index]
    return f"r{rate:g}"


def _gdna_axis(gdna: GDNASimConfig) -> list[tuple[str, float]]:
    return [(gdna_label_for_rate(rate, gdna.rate_labels, i), rate) for i, rate in enumerate(gdna.rates)]


def _genomic_refs(gdna: GDNASimConfig) -> list[str]:
    """The references gDNA is drawn from. Which references are genomic is an input, never inferred: an annotated
    reference can be an RNA-only spike-in."""
    if gdna.genomic_refs is not None:
        return gdna.genomic_refs
    if any(rate > 0 for rate in gdna.rates):
        raise ValueError(
            "gdna.genomic_refs must list the references genomic DNA is drawn from whenever any "
            "gdna.rates entry is non-zero. Every reference is genomic or RNA-only and the "
            "classification is an input, not an inference from the annotation."
        )
    return []


def _overdispersion_axis(gdna: GDNASimConfig) -> list[tuple[str, float]]:
    """The strand-overdispersion axis: one condition per value, or one at ``strand_overdispersion``."""
    values = (
        gdna.strand_overdispersions
        if gdna.strand_overdispersions is not None
        else [gdna.strand_overdispersion]
    )
    labels = gdna.strand_overdispersion_labels
    axis: list[tuple[str, float]] = []
    for i, od in enumerate(values):
        od = float(od)
        if not (0.0 <= od < 1.0):
            raise ValueError(f"gdna strand_overdispersion must be in [0, 1); got {od}")
        if labels is not None and i < len(labels):
            label = labels[i]
        else:
            label = "od00" if od == 0.0 else f"od{od:g}".replace(".", "p")
        axis.append((label, od))
    return axis


def nrna_label_for_ratio(ratio: float | tuple[float, float], labels: list[str] | None, idx: int) -> str:
    if labels and idx < len(labels):
        return labels[idx].strip()
    if isinstance(ratio, tuple):
        lo, hi = ratio
        return f"range_{lo:g}_{hi:g}"
    if ratio == 0:
        return "none"
    return f"ratio_{ratio:g}"


def _nrna_axis(
    nrna: NRNAConfig, has_file_nrna: bool
) -> list[tuple[str, str, float | tuple[float, float] | None, int]]:
    """One ``(label, mode, value, index)`` per nascent condition: per ``ratios`` entry for ``additive_ratio``, per
    ``abundance_ranges`` pair for ``sparse`` (which folds ``index`` into its seed), or the single
    ``("file", "file", None, 0)`` when the abundance file supplied nRNA values.
    """
    if has_file_nrna:
        return [("file", "file", None, 0)]
    values = nrna.ratios if nrna.mode == "additive_ratio" else nrna.abundance_ranges
    return [
        (nrna_label_for_ratio(value, nrna.ratio_labels, i), nrna.mode, value, i)
        for i, value in enumerate(values)
    ]


def capture_axis(cfg: WholeGenomeSimConfig) -> tuple[list[CaptureScenario], bool]:
    """The capture scenarios, and whether condition names carry their labels: only a ``capture.configs`` list
    names them."""
    if cfg.capture_configs:
        return cfg.capture_configs, True
    label = "on" if cfg.capture.probes else "off"
    return [CaptureScenario(label=label, config=cfg.capture)], False


# ── a condition's depth, seed and name ───────────────────────────────────────────────────────────


@dataclass(frozen=True, slots=True)
class ConditionDepths:
    """How one condition's fragment budget splits between RNA and gDNA. The nascent share of ``n_rna`` is not
    imposed: nascent entities are rows of the one RNA multinomial, so it follows from molecules and lengths."""

    n_rna: int
    n_gdna: int

    @property
    def total(self) -> int:
        return self.n_rna + self.n_gdna


def resolve_depths(sim, *, gdna_rate: float) -> ConditionDepths:
    """Split one condition's fragment budget between RNA and gDNA:

    * ``n_total_fragments is None`` (the default): the RNA depth is fixed and gDNA is added at ``rate x n_rna``,
      so a ``rate`` of 1.0 is a 50 % gDNA library at twice the depth of ``rate`` 0.
    * ``n_total_fragments`` set: the total is fixed and ``rate`` decides only the split,
      ``n_rna = total/(1+rate)``, ``n_gdna = total − n_rna``. This mirrors a sequencing run's fixed budget and
      reaches high gDNA fractions at a bounded depth; the RNA side thins as gDNA rises.

    Nascent comes out of ``n_rna``, never on top of it. Raises when the split leaves no RNA.

    Gates: ``tests/test_sim_fixed_total_depth.py``.
    """
    if sim.n_total_fragments is None:
        # additive: RNA depth fixed, gDNA added on top
        n_rna = int(sim.n_rna_fragments)
        return ConditionDepths(n_rna, round(gdna_rate * n_rna))

    # fixed total: `rate` decides the split only
    total = int(sim.n_total_fragments)
    n_rna = round(total / (1.0 + float(gdna_rate)))
    # gDNA must stay the remainder: two shares computed independently drift from the total by a fragment.
    n_gdna = total - n_rna
    if n_rna <= 0:
        raise ValueError(
            f"n_total_fragments={total} at gdna_rate={gdna_rate} leaves {n_rna} RNA fragments. "
            "A condition with no RNA is not an RNA-seq library; raise the total or lower the rate."
        )
    return ConditionDepths(n_rna, n_gdna)


def stable_seed(base_seed: int, *parts: object) -> int:
    """Derive a reproducible 32-bit seed from a base seed and string parts."""
    text = "\0".join([str(base_seed), *(str(part) for part in parts)])
    digest = hashlib.blake2b(text.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "big") & 0xFFFF_FFFF


def capture_paired_condition_seed(base_seed: int, gdna_label: str, strand_specificity: float) -> int:
    """Seed shared by the capture and nascent variants of one ``(gdna, ss)`` base condition.

    A nascent-on cell and its nascent-off twin share the seed but not an identical mature stream, since nascent
    rows are drawn in the same multinomial and re-allocate the RNA budget.
    """
    return stable_seed(base_seed, condition_dir_name(gdna_label, strand_specificity, "_paired"))


def condition_dir_name(
    gdna_label: str,
    strand_specificity: float,
    nrna_label: str,
    capture_label: str | None = None,
    gdna_strand_overdispersion: float = 0.0,
) -> str:
    """Return a condition's directory name, the one definition every consumer shares. The ``_od_<value>`` slug
    appears only when ``gdna_strand_overdispersion > 0``, so a suite that does not sweep it keeps short names."""
    name = f"gdna_{gdna_label}"
    if gdna_strand_overdispersion > 0.0:
        name = f"{name}_od_{gdna_strand_overdispersion:g}".replace(".", "p")
    name = f"{name}_ss_{strand_specificity:.2f}_nrna_{nrna_label}"
    if capture_label is not None:
        name = f"{name}_capture_{capture_label}"
    return name


# ── one condition ────────────────────────────────────────────────────────────────────────────────


def _nascent_transcripts(
    transcripts: list[Transcript], nrna: NRNAConfig, point: tuple
) -> tuple[list[Transcript], dict]:
    """A copy of the base transcripts at one nascent setting, and the manifest fields that record it.

    ``additive_ratio`` pools each entity's molecules from its contributors (``apply_nrna_ratio``), ``sparse`` sets
    each entity's absolute abundance (``apply_sparse_nrna``), and ``file`` keeps the loaded abundances.
    """
    label, mode, value, index = point
    copied = copy.deepcopy(transcripts)
    ratio: float | None = None
    abundance_range: tuple[float, float] | None = None
    if mode == "additive_ratio":
        ratio = float(value or 0.0)
        apply_nrna_ratio(copied, ratio)
    elif mode == "sparse":
        abundance_range = tuple(value)
        ratio = apply_sparse_nrna(
            copied, abundance_range, on_fraction=nrna.on_fraction, seed=nrna.seed + index
        )
    return copied, {
        "nrna_label": label,
        "nrna_mode": mode,
        "nrna_ratio": ratio,
        "nrna_abundance_range": abundance_range,
        "nrna_on_fraction": nrna.on_fraction if mode == "sparse" else None,
    }


def _condition_entry(
    cfg: WholeGenomeSimConfig,
    nascent: dict,
    molar_truth: str,
    gdna: tuple[str, float],
    overdispersion: tuple[str, float],
    strand_specificity: float,
    capture: CaptureScenario,
    capture_in_names: bool,
) -> dict:
    """One condition's manifest entry before it is simulated: its point on every axis, its depths, its seed and the
    paths of its files, each truth table under the key its ``truth_summary.json`` uses."""
    gdna_label, gdna_rate = gdna
    od_label, od = overdispersion
    name = condition_dir_name(
        gdna_label,
        strand_specificity,
        nascent["nrna_label"],
        capture.label if capture_in_names else None,
        gdna_strand_overdispersion=od,
    )
    depths = resolve_depths(cfg.simulation, gdna_rate=gdna_rate)
    return {
        "name": name,
        "gdna_label": gdna_label,
        "gdna_rate": gdna_rate,
        "gdna_strand_overdispersion": od,
        "gdna_strand_overdispersion_label": od_label,
        "strand_specificity": strand_specificity,
        **nascent,
        "capture_label": capture.label,
        "capture_enabled": bool(capture.config.probes),
        "capture_config": capture.config,
        "n_rna": depths.n_rna,
        "n_gdna": depths.n_gdna,
        "n_total": depths.total,
        "seed": capture_paired_condition_seed(cfg.simulation.sim_seed, gdna_label, strand_specificity),
        "truth_kind": "post_capture_empirical",
        "pre_capture_abundances": molar_truth,
        "post_capture_abundances": f"{name}/truth_abundances.tsv",
        "post_capture_fragment_lengths": f"{name}/truth_fragment_lengths.tsv",
        "truth_summary": f"{name}/truth_summary.json",
        "fastq_r1": f"{name}/sim_R1.fq.gz",
        "fastq_r2": f"{name}/sim_R2.fq.gz",
    }


def _simulate_condition(
    cfg: WholeGenomeSimConfig,
    entry: dict,
    transcripts: list[Transcript],
    genomic_refs: list[str],
    sampler: CaptureSampler,
) -> str | None:
    """Simulate one condition into its directory unless its output exists; return its oracle BAM's path in the
    panel, ``None`` when it has none.

    The resume key is the oracle BAM whenever one is requested, so a condition whose BAM exists is skipped even if
    its FASTQs are gone; delete the BAM to regenerate it.
    """
    cond_dir = Path(cfg.outdir) / entry["name"]
    bam = f"{entry['name']}/sim_oracle.bam"
    if (cond_dir / ("sim_oracle.bam" if cfg.oracle_bam else "sim_R1.fq.gz")).exists():
        print("    Output exists, skipping", flush=True)
        return bam if cfg.oracle_bam else None
    print("    Simulating...", end="", flush=True)
    t0 = time.monotonic()
    simulator = WholeGenomeSimulator(
        Path(cfg.genome),
        transcripts,
        replace(cfg.simulation, sim_seed=entry["seed"]),
        replace(cfg.gdna, strand_overdispersion=entry["gdna_strand_overdispersion"]),
        genomic_refs=genomic_refs,
        strand_specificity=entry["strand_specificity"],
        capture_config=entry["capture_config"],
        capture_sampler=sampler,
    )
    _, _, bam_path = simulator.simulate_and_write(
        cond_dir,
        entry["n_rna"],
        entry["n_gdna"],
        oracle_bam=cfg.oracle_bam,
        prefix="sim",
        n_workers=cfg.simulation.n_workers,
    )
    simulator.close()
    print(f" done ({time.monotonic() - t0:.1f}s)", flush=True)
    return bam if bam_path else None


def _write_condition_truth(
    outdir: Path, entry: dict, transcripts: list[Transcript], *, emit_fastq: bool
) -> dict:
    """Write one condition's truth from its reads' names, drop its FASTQs when they are not kept, and return the
    observed origin counts its entry records."""
    cond_dir = outdir / entry["name"]
    bam = cond_dir / "sim_oracle.bam"
    summary = write_post_capture_truth(
        transcripts,
        outdir / entry["post_capture_abundances"],
        outdir / entry["post_capture_fragment_lengths"],
        outdir / entry["truth_summary"],
        bam_path=bam if bam.exists() else None,
        fastq_path=cond_dir / "sim_R1.fq.gz",
        condition=entry["name"],
        pre_capture_abundances=entry["pre_capture_abundances"],
        gdna_strand_overdispersion=entry["gdna_strand_overdispersion"],
    )
    # The FASTQs are dropped only after the truth is written, since it falls back to them.
    if not emit_fastq:
        for fastq in (cond_dir / "sim_R1.fq.gz", cond_dir / "sim_R2.fq.gz"):
            fastq.unlink(missing_ok=True)
    counts = summary["origin_counts"]
    return {f"n_{kind}_observed": int(counts.get(kind, 0)) for kind in ("mrna", "nrna", "gdna")}


# ── the grid ─────────────────────────────────────────────────────────────────────────────────────


def run_condition_grid(
    cfg: WholeGenomeSimConfig, transcripts: list[Transcript], *, has_file_nrna: bool = False
) -> list[dict]:
    """Run every condition *cfg* lays out over *transcripts*, which carry their base abundances, and return the
    per-condition manifest entries for the caller to write. *has_file_nrna* replaces the nascent sweep with the
    abundance file's own levels. Every axis is checked before anything is written."""
    gdna_axis = _gdna_axis(cfg.gdna)
    genomic_refs = _genomic_refs(cfg.gdna)
    od_axis = _overdispersion_axis(cfg.gdna)
    nrna_axis = _nrna_axis(cfg.nrna, has_file_nrna)
    captures, capture_in_names = capture_axis(cfg)

    outdir = Path(cfg.outdir)
    # One sampler per capture scenario, shared by its conditions: the probe layout and the per-width partition
    # depend only on the panel and the templates.
    ref_lengths = load_reference_lengths(Path(cfg.genome))
    samplers = {
        c.label: CaptureSampler.from_config(c.config, transcripts, ref_lengths) for c in captures
    }
    cells = list(product(gdna_axis, od_axis, cfg.strand_specificities, captures))
    total = len(nrna_axis) * len(cells)
    conditions: list[dict] = []
    for point in nrna_axis:
        nascent_transcripts, nascent = _nascent_transcripts(transcripts, cfg.nrna, point)
        molar_truth = f"truth_abundances_nrna_{nascent['nrna_label']}.tsv"
        write_molar_truth(nascent_transcripts, outdir / molar_truth)
        for gdna, od, strand_specificity, capture in cells:
            entry = _condition_entry(
                cfg, nascent, molar_truth, gdna, od, strand_specificity, capture, capture_in_names
            )
            print(
                f"\n  [{len(conditions) + 1}/{total}] {entry['name']}: RNA={entry['n_rna']:,} "
                f"gDNA={entry['n_gdna']:,} SS={strand_specificity:.2f} nRNA={entry['nrna_label']} "
                f"capture={capture.label}",
                flush=True,
            )
            entry["oracle_bam"] = _simulate_condition(
                cfg, entry, nascent_transcripts, genomic_refs, samplers[capture.label]
            )
            entry.update(
                _write_condition_truth(outdir, entry, nascent_transcripts, emit_fastq=cfg.emit_fastq)
            )
            conditions.append(entry)
    return conditions


# ── the manifest ─────────────────────────────────────────────────────────────────────────────────


def _jsonable(value: Any) -> Any:
    if is_dataclass(value):
        return _jsonable(asdict(value))
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def write_manifest(
    outdir: Path,
    cfg: WholeGenomeSimConfig,
    conditions: list[dict[str, Any]],
    *,
    index: dict[str, Any],
    n_shadow_transcripts: int = 0,
) -> Path:
    """Write ``manifest.json`` under ``outdir`` and return its path: the config, the conditions, the absolute input
    paths and the index the transcriptome comes from. ``index`` is ``{"path", "manifest"}`` as
    :func:`_index_provenance` returns it; each truth file's path is in its condition's entry."""
    manifest = {
        "version": 1,
        # the index's resolved path and its own manifest: its FASTA and GTF with their sha256, and its build flags
        "index": _jsonable(index),
        # unannotated transcripts simulated from `shadow_gtf`, so a scored panel says it carried them
        "n_shadow_transcripts": int(n_shadow_transcripts),
        "simulation": _jsonable(cfg.simulation),
        "gdna": _jsonable(cfg.gdna),
        "nrna": _jsonable(cfg.nrna),
        "capture": _jsonable(cfg.capture),
        "capture_configs": _jsonable(cfg.capture_configs),
        "abundance": _jsonable(cfg.abundance),
        "conditions": _jsonable(conditions),
    }
    for key in ("genome", "gtf", "shadow_gtf"):
        value = getattr(cfg, key)
        manifest[key] = str(Path(value).resolve()) if value else value
    manifest["strand_specificities"] = _jsonable(cfg.strand_specificities)
    path = outdir / "manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        json.dump(manifest, handle, indent=2)
    return path


# ── the run ──────────────────────────────────────────────────────────────────────────────────────


def run_simulation(cfg: WholeGenomeSimConfig) -> list[dict]:
    """Run the whole simulation: load the index's transcripts and any shadows, assign the base abundances once, run
    the condition grid and write the manifest. Returns the per-condition manifest entries.
    """
    # checked again: a caller can set `oracle_bam` after the config is read (--no-oracle)
    refuse_a_run_that_keeps_no_reads(cfg)
    outdir = Path(cfg.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if not Path(cfg.genome).exists():
        raise FileNotFoundError(f"Genome not found: {cfg.genome}")
    if not Path(cfg.gtf).exists():
        raise FileNotFoundError(f"GTF not found: {cfg.gtf}")

    index_dir = ensure_index(cfg)
    transcripts = load_transcripts_from_index(index_dir)
    if not transcripts:
        raise RuntimeError(f"No transcripts loaded from index {index_dir}")
    index_record = _index_provenance(index_dir, Path(cfg.gtf))
    n_annotated = len(transcripts)
    if cfg.shadow_gtf:
        if not Path(cfg.shadow_gtf).exists():
            raise FileNotFoundError(f"shadow GTF not found: {cfg.shadow_gtf}")
        transcripts = merge_shadow_transcripts(transcripts, cfg.shadow_gtf)

    has_file_nrna = assign_base_abundances(transcripts, cfg.abundance)
    conditions = run_condition_grid(cfg, transcripts, has_file_nrna=has_file_nrna)
    manifest_path = write_manifest(
        outdir, cfg, conditions, index=index_record, n_shadow_transcripts=len(transcripts) - n_annotated
    )
    logger.info("Wrote manifest to %s", manifest_path)
    return conditions
