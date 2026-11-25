"""The whole-genome condition grid behind `simulate_reads.py`: :func:`run_simulation` and every step it takes. A helper
module, not a script.

:func:`run_simulation` loads the transcript set from the rigel index (``rigel.sim.transcriptome``: annotated
transcripts plus the synthetic nascent entities, ``t_index`` equal to the index row), appends any shadow transcripts,
assigns the base abundances once (``_sim_abundance``) and runs the grid (:func:`run_condition_grid`): nascent x gDNA
rate x gDNA strand overdispersion x strand specificity x capture. It writes ``manifest.json`` last.

Per condition the grid copies the transcripts at their base abundances, applies the nascent mode, splits the fragment
budget (:func:`resolve_depths`), simulates on the capture sampler its panel shares, writes the truth (``_sim_truth``)
and returns a manifest entry. A condition whose output already exists is skipped only when the panel recorded it
under the inputs this config would simulate it from (:func:`refuse_a_stale_resume`); otherwise the run refuses. The
seed keys on ``(gdna, strand specificity)`` only (:func:`capture_paired_condition_seed`), so the capture,
overdispersion and nascent variants of one base condition are paired.
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

import pysam

from rigel.index import TranscriptIndex, load_manifest, load_reference_lengths, source_record
from rigel.sim.capture import CaptureSampler
from rigel.sim.transcriptome import load_transcripts_from_index
from rigel.sim.wgs_config import GDNASimConfig
from rigel.sim.wgs_engine import WholeGenomeSimulator
from rigel.transcript import Transcript

from _sim_abundance import apply_nascent_ratio, apply_sparse_nascent, assign_base_abundances
from _sim_config import CaptureScenario, NascentConfig, WholeGenomeSimConfig, refuse_a_run_that_keeps_no_reads
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
        t.nascent_abundance = 0.0
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
        od = float(od)  # in [0, 1): GDNASimConfig refuses any other
        if labels is not None and i < len(labels):
            label = labels[i]
        else:
            label = "od00" if od == 0.0 else f"od{od:g}".replace(".", "p")
        axis.append((label, od))
    return axis


def nascent_label_for_ratio(ratio: float | tuple[float, float], labels: list[str] | None, idx: int) -> str:
    if labels and idx < len(labels):
        return labels[idx].strip()
    if isinstance(ratio, tuple):
        lo, hi = ratio
        return f"range_{lo:g}_{hi:g}"
    if ratio == 0:
        return "none"
    return f"ratio_{ratio:g}"


def _nascent_axis(
    nascent: NascentConfig, has_file_nascent: bool
) -> list[tuple[str, str, float | tuple[float, float] | None, int]]:
    """One ``(label, mode, value, index)`` per nascent condition: per ``ratios`` entry for ``additive_ratio``, per
    ``abundance_ranges`` pair for ``sparse`` (which folds ``index`` into its seed), or the single
    ``("file", "file", None, 0)`` when the abundance file supplied nascent RNA values.
    """
    if has_file_nascent:
        return [("file", "file", None, 0)]
    values = nascent.ratios if nascent.mode == "additive_ratio" else nascent.abundance_ranges
    return [
        (nascent_label_for_ratio(value, nascent.ratio_labels, i), nascent.mode, value, i)
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
    """Seed shared by the capture and nascent variants of one ``(gdna, ss)`` base condition, keyed on those two
    values, never on a condition's name, so renaming a label or the name's format moves no draw.

    A nascent-on cell and its nascent-off twin share the seed but not an identical mature stream, since nascent
    rows are drawn in the same multinomial and re-allocate the RNA budget.
    """
    return stable_seed(base_seed, gdna_label, f"{strand_specificity:.2f}")


def condition_dir_name(
    gdna_label: str,
    strand_specificity: float,
    nascent_label: str,
    capture_label: str | None = None,
    gdna_strand_overdispersion: float = 0.0,
) -> str:
    """Return a condition's directory name, the one definition every consumer shares. The ``_od_<value>`` slug
    appears only when ``gdna_strand_overdispersion > 0``, so a suite that does not sweep it keeps short names."""
    name = f"gdna_{gdna_label}"
    if gdna_strand_overdispersion > 0.0:
        name = f"{name}_od_{gdna_strand_overdispersion:g}".replace(".", "p")
    name = f"{name}_ss_{strand_specificity:.2f}_nascent_{nascent_label}"
    if capture_label is not None:
        name = f"{name}_capture_{capture_label}"
    return name


# ── one condition ────────────────────────────────────────────────────────────────────────────────


def _nascent_transcripts(
    transcripts: list[Transcript], nascent: NascentConfig, point: tuple
) -> tuple[list[Transcript], dict]:
    """A copy of the base transcripts at one nascent setting, and the manifest fields that record it.

    ``additive_ratio`` pools each entity's molecules from its contributors (``apply_nascent_ratio``), ``sparse`` sets
    each entity's absolute abundance (``apply_sparse_nascent``), and ``file`` keeps the loaded abundances.
    """
    label, mode, value, index = point
    copied = copy.deepcopy(transcripts)
    ratio: float | None = None
    abundance_range: tuple[float, float] | None = None
    if mode == "additive_ratio":
        ratio = float(value or 0.0)
        apply_nascent_ratio(copied, ratio)
    elif mode == "sparse":
        abundance_range = tuple(value)
        ratio = apply_sparse_nascent(
            copied, abundance_range, on_fraction=nascent.on_fraction, seed=nascent.seed + index
        )
    return copied, {
        "nascent_label": label,
        "nascent_mode": mode,
        "nascent_ratio": ratio,
        "nascent_abundance_range": abundance_range,
        "nascent_on_fraction": nascent.on_fraction if mode == "sparse" else None,
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
        nascent["nascent_label"],
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


# ── resuming ─────────────────────────────────────────────────────────────────────────────────────

#: Each finished condition's record of what it was simulated from, written in its directory once its truth is, so a
#: run that dies before the panel's manifest is written still resumes every condition it finished.
CONDITION_RECORD = "condition.json"

#: A condition entry's fields that record what its simulation produced, not what it was asked to simulate.
_OUTCOME_KEYS = frozenset({"oracle_bam", "n_annotated_observed", "n_nascent_observed", "n_gdna_observed"})

#: Manifest fields that lay out the grid (which conditions exist; each entry records its own point on every axis), and
#: so say nothing about what one condition is. `gtf` is a path: the index's record carries the GTF the transcriptome
#: comes from by digest, and the run refuses a `gtf` that is not it. The worker count is identity, not layout: the
#: serial writer and the sharded one draw different reads.
_NOT_IDENTITY = {
    "": frozenset({"conditions", "strand_specificities", "capture", "capture_configs", "gtf"}),
    "gdna": frozenset(
        {"rates", "rate_labels", "strand_overdispersion", "strand_overdispersions", "strand_overdispersion_labels"}
    ),
}


def condition_identity(manifest: dict, name: str) -> dict | None:
    """What condition *name* of a manifest-shaped record was simulated from: the panel's inputs (the index, the
    simulation and gDNA laws, the nascent and abundance models, the genome, the shadow transcripts) and its own entry,
    less what lays out the grid or records the outcome. ``None`` when the record lists no such condition."""
    entry = next((c for c in manifest.get("conditions", []) if c.get("name") == name), None)
    if entry is None:
        return None
    panel = {k: v for k, v in manifest.items() if k not in _NOT_IDENTITY[""]}
    for section, keys in _NOT_IDENTITY.items():
        if section and isinstance(panel.get(section), dict):
            panel[section] = {k: v for k, v in panel[section].items() if k not in keys}
    identity = {"panel": panel, "condition": {k: v for k, v in entry.items() if k not in _OUTCOME_KEYS}}
    # the record is compared as JSON reads it back
    return json.loads(json.dumps(identity))


def _recorded_identity(outdir: Path, name: str, previous: dict | None) -> dict | None:
    """The condition's identity as the panel recorded it: its own record when it has one, else the manifest the panel
    held before this run, else ``None``."""
    record = outdir / name / CONDITION_RECORD
    if record.is_file():
        return json.loads(record.read_text())
    return condition_identity(previous, name) if previous is not None else None


def refuse_a_stale_resume(outdir: Path, name: str, current: dict, previous: dict | None) -> None:
    """Resume a condition whose output exists only when the panel's record of it is what this config would simulate:
    otherwise the old reads would be kept under a manifest that names the new config."""
    recorded = _recorded_identity(outdir, name, previous)
    if recorded is None:
        raise SystemExit(
            f"⛔ {outdir / name} has output but no record of the config it was simulated from (no {CONDITION_RECORD}, "
            f"and the panel's manifest does not list it). Delete the directory to simulate it again."
        )
    if recorded != current:
        moved = sorted(
            f"{part}.{k}"
            for part in ("panel", "condition")
            for k in set(recorded.get(part, {})) | set(current[part])
            if recorded.get(part, {}).get(k) != current[part].get(k)
        )
        raise SystemExit(
            f"⛔ {outdir / name} was simulated under a different config (differs in {moved}). Delete the directory "
            f"to simulate it under this one, or run the config it was simulated from."
        )


def _output(cfg: WholeGenomeSimConfig, entry: dict) -> Path:
    """The file a condition's output is resumed on: its oracle BAM whenever one is requested, else its FASTQ."""
    return Path(cfg.outdir) / entry["name"] / ("sim_oracle.bam" if cfg.oracle_bam else "sim_R1.fq.gz")


def _simulate_condition(
    cfg: WholeGenomeSimConfig,
    entry: dict,
    transcripts: list[Transcript],
    genomic_refs: list[str],
    sampler: CaptureSampler,
) -> str | None:
    """Simulate one condition into its directory unless its output exists; return its oracle BAM's path in the
    panel, ``None`` when it has none. The caller resumes only a condition the panel's record describes
    (:func:`refuse_a_stale_resume`).

    The resume key is the oracle BAM whenever one is requested, so a condition whose BAM exists is skipped even if
    its FASTQs are gone; delete the BAM to regenerate it.
    """
    cond_dir = Path(cfg.outdir) / entry["name"]
    bam = f"{entry['name']}/sim_oracle.bam"
    if _output(cfg, entry).exists():
        print("    Output exists and its record matches this config, resuming", flush=True)
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
    return {f"n_{kind}_observed": int(counts.get(kind, 0)) for kind in ("annotated", "nascent", "gdna")}


# ── the grid ─────────────────────────────────────────────────────────────────────────────────────


def run_condition_grid(
    cfg: WholeGenomeSimConfig,
    transcripts: list[Transcript],
    *,
    has_file_nascent: bool = False,
    index: dict[str, Any] | None = None,
    n_shadow_transcripts: int = 0,
) -> list[dict]:
    """Run every condition *cfg* lays out over *transcripts*, which carry their base abundances, and return the
    per-condition manifest entries for the caller to write. *has_file_nascent* replaces the nascent sweep with the
    abundance file's own levels. Every axis is checked before anything is written. *index* and
    *n_shadow_transcripts* are the manifest's (:func:`write_manifest`): a condition whose output exists resumes only
    when the panel recorded it under the same inputs, and each finished condition records them
    (:data:`CONDITION_RECORD`)."""
    previous_path = Path(cfg.outdir) / "manifest.json"
    previous = json.loads(previous_path.read_text()) if previous_path.is_file() else None
    gdna_axis = _gdna_axis(cfg.gdna)
    genomic_refs = _genomic_refs(cfg.gdna)
    od_axis = _overdispersion_axis(cfg.gdna)
    nascent_axis = _nascent_axis(cfg.nascent, has_file_nascent)
    captures, capture_in_names = capture_axis(cfg)

    outdir = Path(cfg.outdir)
    # One sampler per capture scenario, shared by its conditions: the probe layout, each probe's strength (read from
    # its sequence) and the per-width partition depend only on the panel, the genome and the templates.
    ref_lengths = load_reference_lengths(Path(cfg.genome))
    with pysam.FastaFile(str(cfg.genome)) as fasta:
        samplers = {
            c.label: CaptureSampler.from_config(c.config, transcripts, ref_lengths, fasta=fasta)
            for c in captures
        }
    cells = list(product(gdna_axis, od_axis, cfg.strand_specificities, captures))
    total = len(nascent_axis) * len(cells)
    conditions: list[dict] = []
    for point in nascent_axis:
        nascent_transcripts, nascent = _nascent_transcripts(transcripts, cfg.nascent, point)
        molar_truth = f"truth_abundances_nascent_{nascent['nascent_label']}.tsv"
        write_molar_truth(nascent_transcripts, outdir / molar_truth)
        for gdna, od, strand_specificity, capture in cells:
            entry = _condition_entry(
                cfg, nascent, molar_truth, gdna, od, strand_specificity, capture, capture_in_names
            )
            print(
                f"\n  [{len(conditions) + 1}/{total}] {entry['name']}: RNA={entry['n_rna']:,} "
                f"gDNA={entry['n_gdna']:,} SS={strand_specificity:.2f} nascent RNA={entry['nascent_label']} "
                f"capture={capture.label}",
                flush=True,
            )
            identity = condition_identity(
                manifest_record(cfg, [entry], index=index or {}, n_shadow_transcripts=n_shadow_transcripts),
                entry["name"],
            )
            if _output(cfg, entry).exists():
                refuse_a_stale_resume(outdir, entry["name"], identity, previous)
            entry["oracle_bam"] = _simulate_condition(
                cfg, entry, nascent_transcripts, genomic_refs, samplers[capture.label]
            )
            entry.update(
                _write_condition_truth(outdir, entry, nascent_transcripts, emit_fastq=cfg.emit_fastq)
            )
            (outdir / entry["name"] / CONDITION_RECORD).write_text(json.dumps(identity, indent=2))
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
    """Write ``manifest.json`` under ``outdir`` and return its path: :func:`manifest_record`."""
    manifest = manifest_record(cfg, conditions, index=index, n_shadow_transcripts=n_shadow_transcripts)
    path = outdir / "manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        json.dump(manifest, handle, indent=2)
    return path


def manifest_record(
    cfg: WholeGenomeSimConfig,
    conditions: list[dict[str, Any]],
    *,
    index: dict[str, Any],
    n_shadow_transcripts: int = 0,
) -> dict:
    """The manifest: the config, the conditions, the absolute input paths and the index the transcriptome comes from.
    ``index`` is ``{"path", "manifest"}`` as :func:`_index_provenance` returns it; each truth file's path is in its
    condition's entry."""
    manifest = {
        "version": 1,
        # the index's resolved path and its own manifest: its FASTA and GTF with their sha256, and its build flags
        "index": _jsonable(index),
        # unannotated transcripts simulated from `shadow_gtf`, so a scored panel says it carried them
        "n_shadow_transcripts": int(n_shadow_transcripts),
        "simulation": _jsonable(cfg.simulation),
        "gdna": _jsonable(cfg.gdna),
        "nascent": _jsonable(cfg.nascent),
        "capture": _jsonable(cfg.capture),
        "capture_configs": _jsonable(cfg.capture_configs),
        "abundance": _jsonable(cfg.abundance),
        "conditions": _jsonable(conditions),
    }
    for key in ("genome", "gtf", "shadow_gtf"):
        value = getattr(cfg, key)
        manifest[key] = str(Path(value).resolve()) if value else value
    manifest["strand_specificities"] = _jsonable(cfg.strand_specificities)
    return manifest


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

    has_file_nascent = assign_base_abundances(transcripts, cfg.abundance)
    n_shadow = len(transcripts) - n_annotated
    conditions = run_condition_grid(
        cfg, transcripts, has_file_nascent=has_file_nascent, index=index_record, n_shadow_transcripts=n_shadow
    )
    manifest_path = write_manifest(outdir, cfg, conditions, index=index_record, n_shadow_transcripts=n_shadow)
    logger.info("Wrote manifest to %s", manifest_path)
    return conditions
