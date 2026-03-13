#!/usr/bin/env python
"""synthetic_sim_sweep.py — Combinatorial simulation sweep for Rigel testing.

Entities
--------
- **Transcripts**: defined by strand and exon coordinates.
- **Nascent RNAs (nRNAs)**: defined by genomic coordinates
  ``[strand, start, end]``, matching Rigel's internal nRNA model.
  Transcripts whose span falls within the nRNA coordinates are
  automatically assigned to that group.  Each group has a single
  abundance value; the widest-span multi-exon transcript carries the
  ``nrna_abundance`` in the simulator.
- **gDNA**: single global contamination entity.

Sweep grid
----------
Each entity gets an abundance dimension: constant, list, or
``{start, stop, step}``.  All dimensions form a Cartesian product.

Optional ``patterns`` define correlated parameter bundles — each entry
is one point in a pattern dimension, cross-producted with remaining
sweep dims.  Entities in patterns are excluded from independent sweep.

Config structure (YAML)
-----------------------
::

    genome_length: 50000
    seed: 42
    read_length: 150

    rna: { frag_mean: 200, frag_std: 30, frag_min: 80, frag_max: 450 }
    gdna: { frag_mean: 350, frag_std: 100, frag_min: 100, frag_max: 1000 }

    transcripts:
      TA1:
        strand: "+"
        exons: [[1000, 2000], [5000, 5500], [7000, 7500], [9000, 10000]]
      ...

    nrna_groups:           # optional; auto-derived from overlapping spans
      nrna_TA: [TA1, TA2, TA3, TA4]

    nrnas:                 # coordinate-defined nRNAs (preferred)
      NTA: ["+", 1000, 10000]

    sweep:                 # each entity → constant / list / {start,stop,step}
      nrna_TA: [0, 100]
      gdna: [0, 100]
      strand_specificity: 1.0
      n_fragments: 2000

    params:                # optional; pipeline hyperparameter overrides
      em:                  # EMConfig fields
        prior_alpha: [0.001, 0.01, 0.1]
        strand_symmetry_kappa: 6.0
        mode: "map"
      scan:                # BamScanConfig fields
        strand_prior_kappa: [1.0, 2.0, 4.0]
      scoring:             # FragmentScoringConfig fields
        overhang_log_penalty: -4.605

    patterns:              # optional; correlated bundles
      - {TA1: 100, TA2: 0, TA3: 0, TA4: 3}
      - ...

Usage
-----
    python scripts/synthetic_sim_sweep.py -c scripts/nrna_sweep_config.yaml
    python scripts/synthetic_sim_sweep.py -c config.yaml -o outdir/ -v
    python scripts/synthetic_sim_sweep.py -c config.yaml --gtf genes.gtf -L 50000
    python scripts/synthetic_sim_sweep.py -c config.yaml --dry-run
"""

import argparse
import csv
import dataclasses as _dc
import itertools
import json
import logging
import sys
import tempfile
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path

import yaml

# Ensure src/ is importable when running from the scripts/ directory.
_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "src"))

from rigel.config import (
    EMConfig, PipelineConfig, BamScanConfig, FragmentScoringConfig,
)
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig, run_benchmark

# ---------------------------------------------------------------------------
# Config-field registry (auto-derived from dataclasses)
# ---------------------------------------------------------------------------

_EM_FIELDS = frozenset(f.name for f in _dc.fields(EMConfig))
_SCAN_FIELDS = frozenset(f.name for f in _dc.fields(BamScanConfig))
_SCORING_FIELDS = frozenset(f.name for f in _dc.fields(FragmentScoringConfig))
_ALL_PARAM_FIELDS = _EM_FIELDS | _SCAN_FIELDS | _SCORING_FIELDS

# Mapping from params YAML sub-key to the field set
_PARAM_SECTIONS = {
    "em": _EM_FIELDS,
    "scan": _SCAN_FIELDS,
    "scoring": _SCORING_FIELDS,
}

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Value resolution
# ---------------------------------------------------------------------------

def resolve_values(spec):
    """Resolve a sweep value spec to a concrete list.

    Accepts:
        scalar:      100                        → [100]
        list:        [0, 50, 100]               → [0, 50, 100]
        sweep dict:  {start: 0, stop: 100, step: 50} → [0, 50, 100]
    """
    if isinstance(spec, (int, float)):
        return [spec]
    if isinstance(spec, list):
        return list(spec)
    if isinstance(spec, dict) and "start" in spec and "stop" in spec:
        start = spec["start"]
        stop = spec["stop"]
        step = spec.get("step", 1)
        vals = []
        v = start
        while v <= stop + 1e-9:
            vals.append(round(v, 10))
            v += step
        return vals
    raise ValueError(f"Cannot resolve sweep spec: {spec!r}")


# ---------------------------------------------------------------------------
# Transcript / GTF parsing
# ---------------------------------------------------------------------------

def parse_gtf(gtf_path):
    """Parse a GTF file into transcript definitions.

    Returns
    -------
    transcripts : dict
        t_id → {strand: str, exons: [[start, end], ...]}
    gene_groups : dict
        gene_id → [t_id, ...]  (for nRNA group inference)
    """
    transcripts: dict[str, dict] = {}
    gene_groups: dict[str, list[str]] = {}

    with open(gtf_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            strand = fields[6]
            attrs: dict[str, str] = {}
            for attr in fields[8].rstrip(";").split(";"):
                attr = attr.strip()
                if not attr:
                    continue
                if " " in attr:
                    k, v = attr.split(" ", 1)
                    attrs[k] = v.strip('"')
                elif "=" in attr:
                    k, v = attr.split("=", 1)
                    attrs[k] = v.strip('"')
            gene_id = attrs.get("gene_id", "")
            t_id = attrs.get("transcript_id", "")
            if not t_id:
                continue
            transcripts.setdefault(t_id, {"strand": strand, "exons": []})
            # GTF is 1-based inclusive → 0-based half-open
            transcripts[t_id]["exons"].append(
                [int(fields[3]) - 1, int(fields[4])]
            )
            if gene_id:
                gene_groups.setdefault(gene_id, [])
                if t_id not in gene_groups[gene_id]:
                    gene_groups[gene_id].append(t_id)

    # Sort exons within each transcript
    for t_def in transcripts.values():
        t_def["exons"].sort()

    return transcripts, gene_groups


# ---------------------------------------------------------------------------
# nRNA group computation
# ---------------------------------------------------------------------------

@dataclass
class NrnaGroup:
    """One nascent RNA entity shared by a group of transcripts."""
    label: str
    t_ids: list[str]
    strand: str
    start: int
    end: int
    carrier: str  # transcript ID that carries the nrna_abundance


def _transcript_span(t_def):
    """Return (strand, start, end, n_exons) for a transcript definition."""
    exons = t_def["exons"]
    return (
        t_def["strand"],
        min(e[0] for e in exons),
        max(e[1] for e in exons),
        len(exons),
    )


def _pick_carrier(t_ids, t_spans):
    """Pick the transcript best suited to carry nRNA abundance.

    Prefers the widest-span multi-exon transcript.  Falls back to the
    widest single-exon transcript (the simulator zeros single-exon
    nRNA automatically since it's indistinguishable from mRNA).
    """
    best_multi = (None, -1)
    best_any = (None, -1)
    for t_id in t_ids:
        strand, start, end, n_exons = t_spans[t_id]
        span = end - start
        if span > best_any[1]:
            best_any = (t_id, span)
        if n_exons > 1 and span > best_multi[1]:
            best_multi = (t_id, span)
    return best_multi[0] if best_multi[0] is not None else best_any[0]


def _group_overlapping(t_ids, t_spans):
    """Group transcript IDs by overlapping genomic spans on the same strand.

    Returns list of (t_id_list, strand, merged_start, merged_end).
    """
    by_strand: dict[str, list[tuple[int, int, str]]] = {}
    for t_id in t_ids:
        strand, start, end, _ = t_spans[t_id]
        by_strand.setdefault(strand, []).append((start, end, t_id))

    groups = []
    for strand in sorted(by_strand):
        intervals = sorted(by_strand[strand])
        merged: list[tuple[int, int, list[str]]] = []
        for start, end, t_id in intervals:
            if merged and start <= merged[-1][1]:
                # Overlaps with current group
                prev_start, prev_end, prev_ids = merged[-1]
                merged[-1] = (prev_start, max(end, prev_end),
                              prev_ids + [t_id])
            else:
                merged.append((start, end, [t_id]))
        for start, end, grp_ids in merged:
            groups.append((grp_ids, strand, start, end))
    return groups


def compute_nrna_groups(transcripts_config, user_groups=None):
    """Compute nRNA groups from transcript definitions.

    Parameters
    ----------
    transcripts_config : dict
        t_id → {strand, exons}
    user_groups : dict or None
        label → [t_id, ...]  (from ``nrna_groups`` in YAML)

    Returns
    -------
    dict[str, NrnaGroup]
    """
    t_spans = {t_id: _transcript_span(t_def)
               for t_id, t_def in transcripts_config.items()}

    assigned: set[str] = set()
    groups: dict[str, NrnaGroup] = {}

    # 1. User-defined groups
    if user_groups:
        for label, t_ids in user_groups.items():
            for t_id in t_ids:
                if t_id not in transcripts_config:
                    raise ValueError(
                        f"nrna_groups: transcript '{t_id}' in group "
                        f"'{label}' not found in transcripts")
            strands = set(t_spans[t][0] for t in t_ids)
            if len(strands) > 1:
                raise ValueError(
                    f"nrna_groups: group '{label}' contains transcripts "
                    f"on different strands: {strands}")
            strand = strands.pop()
            start = min(t_spans[t][1] for t in t_ids)
            end = max(t_spans[t][2] for t in t_ids)
            carrier = _pick_carrier(t_ids, t_spans)
            groups[label] = NrnaGroup(label, list(t_ids), strand,
                                      start, end, carrier)
            assigned.update(t_ids)

    # 2. Auto-derive remaining transcripts from overlapping spans
    remaining = [t_id for t_id in transcripts_config if t_id not in assigned]
    if remaining:
        for grp_ids, strand, start, end in _group_overlapping(remaining,
                                                               t_spans):
            if len(grp_ids) == 1:
                label = f"nrna_{grp_ids[0]}"
            else:
                # Find a label that doesn't conflict
                label = f"nrna_{grp_ids[0]}"
                if label in groups:
                    label = f"nrna_{strand}_{start}_{end}"
            carrier = _pick_carrier(grp_ids, t_spans)
            groups[label] = NrnaGroup(label, grp_ids, strand,
                                      start, end, carrier)

    return groups


def parse_nrna_coords(nrnas_config, transcripts_config):
    """Build nRNA groups from coordinate definitions.

    Each nRNA is defined as ``label → [strand, start, end]``.
    Transcripts whose genomic span is contained within the nRNA's
    coordinates (same strand) are automatically assigned to that group.
    Unassigned transcripts get auto-derived singleton groups.

    Raises
    ------
    ValueError
        If coordinates have the wrong shape or no transcript falls
        within the specified region.
    """
    t_spans = {t_id: _transcript_span(t_def)
               for t_id, t_def in transcripts_config.items()}

    assigned: set[str] = set()
    groups: dict[str, NrnaGroup] = {}

    # Sort nRNA entries by span size (ascending) so that more specific
    # (smaller) groups take priority over broader ones when transcripts
    # fall within multiple overlapping nRNA coordinate regions.
    sorted_nrnas = sorted(
        nrnas_config.items(),
        key=lambda item: int(item[1][2]) - int(item[1][1]),
    )

    for label, coords in sorted_nrnas:
        if not isinstance(coords, list) or len(coords) != 3:
            raise ValueError(
                f"nrna '{label}': expected [strand, start, end], "
                f"got {coords!r}")
        strand = str(coords[0])
        start, end = int(coords[1]), int(coords[2])

        # Match unassigned transcripts within (strand, start, end)
        matching: list[str] = []
        for t_id, (t_strand, t_start, t_end, _) in t_spans.items():
            if t_id in assigned:
                continue
            if t_strand == strand and t_start >= start and t_end <= end:
                matching.append(t_id)

        if not matching:
            logger.warning(
                "nrna '%s' at (%s, %d, %d) has no unassigned transcripts "
                "(all already claimed by more specific nRNA groups)",
                label, strand, start, end)
            continue

        carrier = _pick_carrier(matching, t_spans)
        groups[label] = NrnaGroup(
            label=label, t_ids=matching, strand=strand,
            start=start, end=end, carrier=carrier)
        assigned.update(matching)

    # Auto-derive groups for unassigned transcripts
    remaining = [t_id for t_id in transcripts_config if t_id not in assigned]
    if remaining:
        for grp_ids, strand, start, end in _group_overlapping(remaining,
                                                               t_spans):
            lbl = f"nrna_{grp_ids[0]}"
            if lbl in groups:
                lbl = f"nrna_{strand}_{start}_{end}"
            carrier = _pick_carrier(grp_ids, t_spans)
            groups[lbl] = NrnaGroup(
                label=lbl, t_ids=grp_ids, strand=strand,
                start=start, end=end, carrier=carrier)

    return groups


# ---------------------------------------------------------------------------
# Sweep grid builder
# ---------------------------------------------------------------------------

def _parse_params_section(params_config):
    """Parse the ``params:`` YAML section into sweepable dimensions.

    Supports nested sub-keys (``em:``, ``scan:``, ``scoring:``) and
    validated against known config-dataclass fields.  Each field value
    uses the same scalar / list / ``{start, stop, step}`` syntax as
    sweep dimensions.

    Returns
    -------
    param_dims : dict[str, list]
        field_name → resolved value list
    """
    if not params_config:
        return {}

    param_dims: dict[str, list] = {}

    for section_key, field_set in _PARAM_SECTIONS.items():
        section = params_config.get(section_key, {})
        if not isinstance(section, dict):
            continue
        for k, v in section.items():
            if k not in field_set:
                logger.warning(
                    "Unknown param '%s' in params.%s — ignoring", k, section_key)
                continue
            param_dims[k] = resolve_values(v)

    # Also accept flat (un-nested) keys for convenience
    for k, v in params_config.items():
        if k in _PARAM_SECTIONS:
            continue
        if k in _ALL_PARAM_FIELDS:
            if k not in param_dims:  # nested takes precedence
                param_dims[k] = resolve_values(v)
        else:
            logger.warning(
                "Unknown param '%s' in params section — ignoring", k)

    return param_dims


def build_sweep_grid(sweep_config, patterns_config, all_t_ids, nrna_labels,
                     param_dims=None):
    """Build combinatorial dimensions from sweep + patterns + params.

    Entities in patterns are removed from the Cartesian product of
    sweep dims and replaced by a single pattern-index dimension.
    Hyperparameters from the ``params:`` section are added as additional
    sweep dimensions.

    Returns (dims, linked_names) where dims is an OrderedDict.

    Raises
    ------
    ValueError
        If any entity appears in both ``sweep`` and ``patterns``.
    """
    # Validate no entity in both sweep and patterns
    linked_names: set[str] = set()
    if patterns_config:
        for entry in patterns_config:
            linked_names.update(entry.keys())
        conflict = linked_names & set(sweep_config.keys())
        if conflict:
            raise ValueError(
                f"Entities cannot appear in both 'sweep' and 'patterns': "
                f"{sorted(conflict)}")

    dims: OrderedDict[str, list] = OrderedDict()

    # Patterns become a single index dimension
    if patterns_config:
        dims["_pattern_idx"] = list(range(len(patterns_config)))

    # Independent sweep params (skip pattern-linked entities and
    # config-field names — those belong in params:, but we accept
    # them here for backward compatibility)
    for key, spec in sweep_config.items():
        if key in linked_names:
            continue
        dims[key] = resolve_values(spec)

    # Default: transcript IDs not in sweep or patterns → [0]
    for t_id in all_t_ids:
        if t_id not in dims and t_id not in linked_names:
            dims[t_id] = [0]

    # Default: nRNA groups not in sweep or patterns → [0]
    for label in nrna_labels:
        if label not in dims and label not in linked_names:
            dims[label] = [0]

    # Core params with defaults
    dims.setdefault("gdna", [0])
    dims.setdefault("strand_specificity", [1.0])

    # Fragment count: mode (a) n_fragments=fixed total (default),
    # or mode (b) n_rna_fragments + gdna_fraction
    has_rna_mode = "n_rna_fragments" in dims
    if not has_rna_mode:
        dims.setdefault("n_fragments", [2000])
    if has_rna_mode:
        dims.setdefault("gdna_fraction", [0.0])

    # Merge params dimensions (params: section takes precedence over
    # legacy config keys that happen to appear in sweep:)
    if param_dims:
        for k, v in param_dims.items():
            dims[k] = v  # overwrite if already in dims from sweep:

    return dims, linked_names


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def run_sweep(config, output_dir, *, gtf_path=None,
              genome_length_override=None, dry_run=False):
    """Run the full combinatorial sweep."""
    genome_length = genome_length_override or config.get("genome_length", 50000)
    seed = config.get("seed", 42)

    # -- Fragment distribution parameters --
    rna_frag_params = config.get("rna", {})
    gdna_frag_params = config.get("gdna", {})
    read_length = config.get("read_length", 150)

    # -- Parse transcripts --
    gtf_gene_groups = None
    if gtf_path:
        transcripts_config, gtf_gene_groups = parse_gtf(gtf_path)
    elif "transcripts" in config:
        transcripts_config = config["transcripts"]
    else:
        raise SystemExit(
            "No transcript definitions: provide 'transcripts' in YAML "
            "or --gtf on the command line")

    all_t_ids = list(transcripts_config.keys())

    # -- Compute nRNA groups --
    nrnas_config = config.get("nrnas")
    if nrnas_config:
        nrna_groups = parse_nrna_coords(nrnas_config, transcripts_config)
    else:
        user_groups = config.get("nrna_groups")
        if user_groups is None and gtf_gene_groups:
            user_groups = {f"nrna_{gid}": t_ids
                           for gid, t_ids in gtf_gene_groups.items()}
        nrna_groups = compute_nrna_groups(transcripts_config, user_groups)
    nrna_labels = list(nrna_groups.keys())

    # Log nRNA groups
    logger.info("nRNA groups (%d):", len(nrna_groups))
    for label, grp in nrna_groups.items():
        logger.info("  %s: %s (%s-strand, %d–%d, carrier=%s)",
                     label, grp.t_ids, grp.strand, grp.start, grp.end,
                     grp.carrier)

    # -- Build sweep grid --
    sweep_config = config.get("sweep", {})
    patterns_config = config.get("patterns")
    param_dims = _parse_params_section(config.get("params", {}))
    dims, linked_names = build_sweep_grid(
        sweep_config, patterns_config, all_t_ids, nrna_labels,
        param_dims=param_dims)

    # Cartesian product
    dim_names = list(dims.keys())
    combinations = list(itertools.product(*(dims[k] for k in dim_names)))

    # Log grid
    logger.info("Sweep grid:")
    for name in dim_names:
        if name == "_pattern_idx":
            logger.info("  patterns: %d configurations", len(patterns_config))
            for i, p in enumerate(patterns_config):
                logger.info("    [%d] %s", i, p)
        else:
            logger.info("  %s: %s", name, dims[name])
    logger.info("Total runs: %d", len(combinations))

    if dry_run:
        print(f"\nDry run: {len(combinations)} combinations would be executed.")
        return

    # -- Resolve CSV columns --
    param_cols: list[str] = []
    if patterns_config:
        seen: set[str] = set()
        for entry in patterns_config:
            for k in entry:
                if k not in seen:
                    param_cols.append(k)
                    seen.add(k)
    for name in dim_names:
        if name != "_pattern_idx" and name not in param_cols:
            param_cols.append(name)

    result_cols = (
        [f"{t}_expected" for t in all_t_ids]
        + [f"{t}_observed" for t in all_t_ids]
        + [f"{t}_abs_diff" for t in all_t_ids]
        + [f"{t}_rel_err" for t in all_t_ids]
        + ["nrna_expected", "nrna_observed", "nrna_abs_diff", "nrna_rel_err"]
        + ["gdna_expected", "gdna_observed", "gdna_abs_diff"]
        + ["total_mrna_expected", "total_mrna_observed", "total_mrna_rel_err"]
        + ["total_rna_expected", "total_rna_observed", "total_rna_rel_err"]
        + ["n_fragments_actual", "n_intergenic", "n_chimeric"]
        # Fragment length distribution recovery metrics
        + ["rna_fl_true_mean", "rna_fl_est_mean", "rna_fl_mean_err"]
        + ["rna_fl_true_std", "rna_fl_est_std"]
        + ["gdna_fl_true_mean", "gdna_fl_est_mean", "gdna_fl_mean_err"]
        + ["gdna_fl_true_std", "gdna_fl_est_std"]
        + ["rna_fl_n_obs", "gdna_fl_n_obs"]
    )
    fieldnames = param_cols + result_cols

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    tsv_path = out / "sweep_results.tsv"
    json_path = out / "sweep_results.json"

    summary_rows: list[dict] = []
    all_rows: list[dict] = []

    with open(tsv_path, "w", newline="") as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames,
                                delimiter="\t")
        writer.writeheader()

        for run_idx, combo in enumerate(combinations):
            raw = dict(zip(dim_names, combo))

            # Merge pattern + individual parameters
            params: dict[str, float] = {}
            if patterns_config and "_pattern_idx" in raw:
                params.update(patterns_config[int(raw["_pattern_idx"])])
            for k, v in raw.items():
                if k != "_pattern_idx":
                    params[k] = v

            # Human-readable label
            parts: list[str] = []
            for t in all_t_ids:
                v = params.get(t, 0)
                if v:
                    parts.append(f"{t}={v}")
            for label in nrna_labels:
                if params.get(label, 0) > 0:
                    parts.append(f"{label}={params[label]}")
            if params.get("gdna", 0) > 0:
                parts.append(f"gdna={params['gdna']}")
            ss = float(params.get("strand_specificity", 1.0))
            if ss < 1.0:
                parts.append(f"ss={ss}")
            for ek in sorted(_ALL_PARAM_FIELDS & params.keys()):
                parts.append(f"{ek}={params[ek]}")
            run_label = " ".join(parts) or "baseline"

            logger.info("=" * 60)
            logger.info("Run %d/%d: %s", run_idx + 1,
                        len(combinations), run_label)

            n_frags = int(params.get("n_fragments", 2000))
            n_rna_frags = (int(params["n_rna_fragments"])
                           if "n_rna_fragments" in params else None)
            gdna_frac = (float(params["gdna_fraction"])
                         if "gdna_fraction" in params else None)
            gdna_ab = float(params.get("gdna", 0))

            with tempfile.TemporaryDirectory(prefix="rigel_sweep_") as tmpdir:
                sc = Scenario(
                    "sweep",
                    genome_length=genome_length,
                    seed=seed,
                    work_dir=Path(tmpdir) / "work",
                )

                # Add transcripts grouped by nRNA group (gene_id = group label)
                for grp_label, grp in nrna_groups.items():
                    nrna_ab = float(params.get(grp_label, 0))
                    t_list = []
                    for t_id in grp.t_ids:
                        t_def = transcripts_config[t_id]
                        ab = float(params.get(t_id, 0))
                        # Only the carrier gets the nRNA abundance
                        nrna = nrna_ab if t_id == grp.carrier else 0.0
                        t_list.append({
                            "t_id": t_id,
                            "exons": [tuple(e) for e in t_def["exons"]],
                            "abundance": ab,
                            "nrna_abundance": nrna,
                        })
                    sc.add_gene(grp_label, grp.strand, t_list)

                # gDNA config — needed when gdna abundance > 0 or
                # gdna_fraction > 0 (mode b)
                needs_gdna = gdna_ab > 0 or (gdna_frac and gdna_frac > 0)
                gdna_cfg = None
                if needs_gdna:
                    gdna_cfg = GDNAConfig(
                        abundance=gdna_ab if gdna_ab > 0 else 1.0,
                        frag_mean=gdna_frag_params.get("frag_mean", 350),
                        frag_std=gdna_frag_params.get("frag_std", 100),
                        frag_min=gdna_frag_params.get("frag_min", 100),
                        frag_max=gdna_frag_params.get("frag_max", 1000),
                    )

                sim_cfg = SimConfig(
                    frag_mean=rna_frag_params.get("frag_mean", 200),
                    frag_std=rna_frag_params.get("frag_std", 30),
                    frag_min=rna_frag_params.get("frag_min", 80),
                    frag_max=rna_frag_params.get("frag_max", 450),
                    read_length=read_length,
                    strand_specificity=ss,
                    seed=seed,
                )

                # nrna_abundance=0 preserves per-transcript values from add_gene
                result = sc.build_oracle(
                    n_fragments=n_frags,
                    sim_config=sim_cfg,
                    gdna_config=gdna_cfg,
                    nrna_abundance=0.0,
                    n_rna_fragments=n_rna_frags,
                    gdna_fraction=gdna_frac,
                )

                # -- Build pipeline config from params --
                em_kwargs = {k: params[k] for k in _EM_FIELDS
                             if k in params and k != "seed"}
                scan_kwargs: dict = {"sj_strand_tag": "auto"}
                scan_kwargs.update(
                    {k: params[k] for k in _SCAN_FIELDS if k in params})
                scoring_kwargs = {k: params[k] for k in _SCORING_FIELDS
                                  if k in params}

                pipe_cfg = PipelineConfig(
                    em=EMConfig(seed=seed, **em_kwargs),
                    scan=BamScanConfig(**scan_kwargs),
                    scoring=FragmentScoringConfig(**scoring_kwargs),
                )
                pr = run_pipeline(
                    result.bam_path, result.index, config=pipe_cfg,
                )
                bench = run_benchmark(
                    result, pr, scenario_name=f"sweep_{run_idx}",
                )

                logger.info("\n%s", bench.summary())

                # -- Collect CSV row --
                row: dict = {col: params.get(col, 0) for col in param_cols}

                t_map = {ta.t_id: ta for ta in bench.transcripts}
                for t_id in all_t_ids:
                    ta = t_map.get(t_id)
                    if ta:
                        row[f"{t_id}_expected"] = ta.expected
                        row[f"{t_id}_observed"] = round(ta.observed, 2)
                        row[f"{t_id}_abs_diff"] = round(ta.abs_diff, 2)
                        row[f"{t_id}_rel_err"] = round(ta.rel_error, 4)
                    else:
                        for suf in ("_expected", "_observed",
                                    "_abs_diff", "_rel_err"):
                            row[f"{t_id}{suf}"] = 0

                # nRNA
                row["nrna_expected"] = bench.n_nrna_expected
                row["nrna_observed"] = round(bench.n_nrna_pipeline, 2)
                row["nrna_abs_diff"] = round(bench.nrna_abs_diff, 2)
                nrna_rel = (bench.nrna_abs_diff / bench.n_nrna_expected
                            if bench.n_nrna_expected > 0 else 0.0)
                row["nrna_rel_err"] = round(nrna_rel, 4)

                # gDNA
                row["gdna_expected"] = bench.n_gdna_expected
                row["gdna_observed"] = round(bench.n_gdna_pipeline, 2)
                row["gdna_abs_diff"] = round(bench.gdna_abs_diff, 2)

                # Aggregated RNA totals
                row["total_mrna_expected"] = bench.total_expected
                row["total_mrna_observed"] = round(bench.total_observed, 2)
                mrna_rel = (abs(bench.total_observed - bench.total_expected)
                            / bench.total_expected
                            if bench.total_expected > 0 else 0.0)
                row["total_mrna_rel_err"] = round(mrna_rel, 4)

                trna_exp = bench.total_expected + bench.n_nrna_expected
                trna_obs = bench.total_observed + bench.n_nrna_pipeline
                rna_rel = (abs(trna_obs - trna_exp) / trna_exp
                           if trna_exp > 0 else 0.0)
                row["total_rna_expected"] = trna_exp
                row["total_rna_observed"] = round(trna_obs, 2)
                row["total_rna_rel_err"] = round(rna_rel, 4)

                row["n_fragments_actual"] = bench.n_fragments
                row["n_intergenic"] = bench.n_intergenic
                row["n_chimeric"] = bench.n_chimeric

                # Fragment length distribution recovery metrics
                fl = pr.frag_length_models
                rna_m = fl.rna_model
                gdna_m = fl.gdna_model

                rna_true_mean = rna_frag_params.get("frag_mean", 200)
                rna_true_std = rna_frag_params.get("frag_std", 30)
                row["rna_fl_true_mean"] = rna_true_mean
                row["rna_fl_true_std"] = rna_true_std
                row["rna_fl_n_obs"] = rna_m.n_observations
                if rna_m.n_observations > 0:
                    row["rna_fl_est_mean"] = round(rna_m.mean, 2)
                    row["rna_fl_est_std"] = round(rna_m.std, 2)
                    row["rna_fl_mean_err"] = round(
                        rna_m.mean - rna_true_mean, 2)
                else:
                    row["rna_fl_est_mean"] = ""
                    row["rna_fl_est_std"] = ""
                    row["rna_fl_mean_err"] = ""

                gdna_true_mean = gdna_frag_params.get("frag_mean", 350)
                gdna_true_std = gdna_frag_params.get("frag_std", 100)
                row["gdna_fl_true_mean"] = gdna_true_mean
                row["gdna_fl_true_std"] = gdna_true_std
                row["gdna_fl_n_obs"] = gdna_m.n_observations
                if gdna_m.n_observations > 0:
                    row["gdna_fl_est_mean"] = round(gdna_m.mean, 2)
                    row["gdna_fl_est_std"] = round(gdna_m.std, 2)
                    row["gdna_fl_mean_err"] = round(
                        gdna_m.mean - gdna_true_mean, 2)
                else:
                    row["gdna_fl_est_mean"] = ""
                    row["gdna_fl_est_std"] = ""
                    row["gdna_fl_mean_err"] = ""

                writer.writerow(row)
                tsvfile.flush()
                all_rows.append(row)

                summary_rows.append({
                    "run": run_idx + 1,
                    "label": run_label,
                    "mrna_rel": mrna_rel,
                    "nrna_exp": bench.n_nrna_expected,
                    "nrna_diff": bench.nrna_abs_diff,
                    "gdna_exp": bench.n_gdna_expected,
                    "gdna_diff": bench.gdna_abs_diff,
                })

    # Write JSON
    with open(json_path, "w") as jf:
        json.dump(all_rows, jf, indent=2)

    _print_summary(summary_rows, tsv_path, json_path)


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def _print_summary(rows, tsv_path, json_path):
    sep = "=" * 74
    print(f"\n{sep}")
    print(f"SWEEP COMPLETE — {len(rows)} runs")
    print(f"Results: {tsv_path}")
    print(f"         {json_path}")
    print(sep)

    hdr = f"{'#':>3}  {'Configuration':<40}  {'mRNA%':>6}  {'nRNA':>6}  {'gDNA':>6}"
    print(hdr)
    print("-" * len(hdr))

    for r in rows:
        lbl = r["label"][:40]
        m = f"{r['mrna_rel']*100:.1f}%" if r["mrna_rel"] > 0.001 else "ok"
        n = f"{r['nrna_diff']:.0f}" if r["nrna_exp"] > 0 else "."
        g = f"{r['gdna_diff']:.0f}" if r["gdna_exp"] > 0 else "."
        print(f"{r['run']:>3}  {lbl:<40}  {m:>6}  {n:>6}  {g:>6}")
    print()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Combinatorial simulation sweep for Rigel nRNA/gDNA testing.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Sweep value specifications (YAML):
    scalar:   100                             → [100]
    list:     [0, 50, 100]                    → [0, 50, 100]
    sweep:    {start: 0, stop: 100, step: 50} → [0, 50, 100]

Examples:
  python scripts/synthetic_sim_sweep.py -c scripts/nrna_sweep_config.yaml
  python scripts/synthetic_sim_sweep.py -c config.yaml -o outdir/ -v
  python scripts/synthetic_sim_sweep.py -c config.yaml --gtf genes.gtf -L 50000
  python scripts/synthetic_sim_sweep.py -c config.yaml --dry-run
""",
    )
    parser.add_argument("-c", "--config", required=True,
                        help="YAML sweep configuration file")
    parser.add_argument("-g", "--gtf", default=None,
                        help="GTF file for transcript structures "
                             "(overrides 'transcripts' in YAML)")
    parser.add_argument("-L", "--genome-length", type=int, default=None,
                        help="Genome length in bp (overrides YAML)")
    parser.add_argument("-o", "--outdir", default="sweep_results",
                        help="Output directory (default: sweep_results)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show sweep grid without running")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable debug-level logging")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-5s %(name)s — %(message)s",
    )

    with open(args.config) as f:
        config = yaml.safe_load(f)

    run_sweep(
        config,
        args.outdir,
        gtf_path=args.gtf,
        genome_length_override=args.genome_length,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()
