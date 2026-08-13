"""⭐⭐ NET FRAGMENT FLOW — WHERE DID EACH MISASSIGNED FRAGMENT GO?

    Gate: ``tests/test_net_flow.py``

⭐ **This answers a question no other instrument does.** ``quant_accuracy.py`` measures the MAGNITUDE of
transcript error; this decomposes its DIRECTION — each transcript's surplus or deficit split into
gDNA-sourced and RNA-isoform-sourced flow, over a confusion matrix of true origin against assigned
destination. "6 M gDNA fragments left the gDNA pool at ``g98 ss0.50 capture_on``" is a `quant_accuracy`
number; "and here is the transcript they landed on" is this one.

⛔ **THIS MODULE IS WHAT SURVIVED ``sim/analysis.py``** (retired 2026-08-11, owner). That file was a
1,589-boundary SECOND SCORER: it ran the tool and rendered its own accuracy tables beside
``quant_accuracy.py``'s, against its own definition of truth. Two scorers is how a baseline and a ceiling
drift apart (`TRAPS: score-the-consumers-own-count`), so the scoring and report-rendering halves went and
the flow decomposition — the part with no duplicate — moved here with its tests.

⚠ The truth loaders below came along because the flow analysis needs them; they are panel plumbing, not
a second definition of truth.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .manifest import condition_manifest_map, load_manifest
from .read_name import parse_origin



_POOLS_3 = ("gdna", "nrna", "mrna")


def load_truth(sim_base: Path, truth_name: str | None = None) -> pd.DataFrame:
    """Load a ground-truth abundance table."""
    manifest = load_manifest(sim_base)
    if truth_name is None:
        truth_name = manifest.get("truth_abundances", "truth_abundances.tsv")
    truth_path = sim_base / truth_name
    if not truth_path.exists():
        condition_map = condition_manifest_map(manifest)
        truth_names = [
            str(row.get("truth_abundances"))
            for row in condition_map.values()
            if row.get("truth_abundances")
        ]
        if truth_names:
            truth_path = sim_base / truth_names[0]
    return pd.read_csv(truth_path, sep="\t")


def load_condition_truth(
    sim_base: Path,
    condition: str,
    condition_meta: dict[str, dict],
    fallback_truth: pd.DataFrame,
) -> pd.DataFrame:
    """Load condition-specific truth, falling back to the supplied table."""
    truth_name = condition_meta.get(condition, {}).get("truth_abundances")
    if not truth_name:
        return fallback_truth
    truth_path = sim_base / str(truth_name)
    if not truth_path.exists():
        return fallback_truth
    return pd.read_csv(truth_path, sep="\t")


def parse_condition(cond_name: str) -> dict:
    """Parse a ``gdna_<lbl>_ss_<val>_nrna_<lbl>[...]`` condition name into its parts.

    Structural only — the gDNA *rate* comes from the manifest (``gdna_rate``), never the name.
    """
    parts = cond_name.split("_")
    return {
        "condition": cond_name,
        "strand_specificity": float(parts[3]),
        "gdna_label": f"{parts[0]}_{parts[1]}",
        "nrna_label": parts[5] if len(parts) > 5 else "none",
    }


def load_loci(out_dir: Path) -> pd.DataFrame:
    """Load loci output."""
    feather_path = out_dir / "loci.feather"
    tsv_path = out_dir / "loci.tsv"
    if feather_path.exists():
        return pd.read_feather(feather_path)
    elif tsv_path.exists():
        return pd.read_csv(tsv_path, sep="\t")
    else:
        return pd.DataFrame()


def collect_fragment_flows(
    sim_base: Path, conditions: list[str]
) -> tuple[dict[str, FlowData], list[dict]]:
    """Single pass over each condition's annotated BAM.

    Returns ``(flows_by_condition, overview_rows)``. ``overview_rows`` reproduces the
    legacy gross-confusion schema (consumed by the condition report + acceptance checks);
    ``flows_by_condition`` carries the sparse per-locus flow matrix for the net analysis.
    """
    import pysam

    flows_by_cond: dict[str, FlowData] = {}
    overview_rows: list[dict] = []

    for cond in conditions:
        cond_dir = sim_base / cond
        annotated_bam = cond_dir / "annotated.bam"
        if not annotated_bam.exists():
            continue

        # Component registry (lazy integer ids).
        key2cid: dict[tuple, int] = {}
        comp_name: dict[int, str] = {}
        comp_kind: dict[int, str] = {}
        comp_locus: dict[int, int] = {}
        # RNA components' home locus is the modal ZL across fragments touching them
        # (BAM-derived, self-consistent with the ZL-keyed gDNA components). This puts
        # nRNA-shadow spans — which carry locus_id=-1 in quant.feather — into their
        # true locus, so mRNA→nRNA flow stays within-locus instead of leaking to cross_locus.
        comp_zl: dict[int, Counter] = defaultdict(Counter)

        def cid_for_transcript(tid: str) -> int:
            key = ("t", tid)
            cid = key2cid.get(key)
            if cid is None:
                cid = len(key2cid)
                key2cid[key] = cid
                comp_name[cid] = tid
                comp_kind[cid] = "rna"
                comp_locus[cid] = -1  # finalized to modal ZL after the pass
            return cid

        def cid_for_gdna(locus: int) -> int:
            key = ("g", locus)
            cid = key2cid.get(key)
            if cid is None:
                cid = len(key2cid)
                key2cid[key] = cid
                comp_name[cid] = f"gdna@{locus}"
                comp_kind[cid] = "gdna"
                comp_locus[cid] = locus
            return cid

        def cid_unassigned() -> int:
            key = ("u",)
            cid = key2cid.get(key)
            if cid is None:
                cid = len(key2cid)
                key2cid[key] = cid
                comp_name[cid] = "unassigned"
                comp_kind[cid] = "unassigned"
                comp_locus[cid] = -1
            return cid

        flow: Counter = Counter()
        pool_flow: Counter = Counter()  # (true_pool, assigned_pool) -> count, pools gdna/nrna/mrna

        # Legacy gross-confusion counters (back-compat for downstream consumers).
        correct_tx = correct_gene = wrong_tx = 0
        rna_as_gdna = gdna_as_rna = gdna_correct = 0
        total_rna = total_mrna = total_nrna = total_gdna = 0
        nrna_as_rna = nrna_as_gdna = mrna_as_gdna = mrna_as_nrna = 0

        with pysam.AlignmentFile(str(annotated_bam), "rb") as bam:
            for read in bam:
                if read.is_read2 or read.is_secondary or read.is_supplementary:
                    continue

                origin = parse_origin(read.query_name)

                def _tag(name: str, default):
                    try:
                        return read.get_tag(name)
                    except KeyError:
                        return default

                zt = str(_tag("ZT", "") or "")
                zg = str(_tag("ZG", "") or "")
                category = str(_tag("ZC", "") or "")
                zf = int(_tag("ZF", 0) or 0)
                zl = int(_tag("ZL", -1))
                # ZF bit 2 (0x04) = gDNA pool; intergenic also carries it but guard on ZC too.
                is_assigned_gdna = bool(zf & 0x04) or category == "intergenic"
                is_assigned_nrna = bool(zf & 0x08)

                # Assigned component (the tool's hard MAP destination).
                if is_assigned_gdna:
                    a_cid = cid_for_gdna(zl)
                elif zt and zt != ".":
                    a_cid = cid_for_transcript(zt)
                else:
                    a_cid = cid_unassigned()

                # True component (oracle) + legacy gross counters.
                if origin.kind == "gdna":
                    total_gdna += 1
                    t_cid = cid_for_gdna(zl)
                    if is_assigned_gdna:
                        gdna_correct += 1
                    else:
                        gdna_as_rna += 1
                else:
                    true_tx = origin.transcript_id or ""
                    t_cid = cid_for_transcript(true_tx)
                    total_rna += 1
                    if origin.kind == "nrna":
                        total_nrna += 1
                        if is_assigned_gdna:
                            rna_as_gdna += 1
                            nrna_as_gdna += 1
                        else:
                            nrna_as_rna += 1
                    else:
                        total_mrna += 1
                        if is_assigned_gdna:
                            rna_as_gdna += 1
                            mrna_as_gdna += 1
                        elif is_assigned_nrna:
                            mrna_as_nrna += 1
                        elif zt == true_tx:
                            correct_tx += 1
                        elif zg and true_tx.rsplit(".", 1)[0] == zg:
                            correct_gene += 1
                        else:
                            wrong_tx += 1

                flow[(t_cid, a_cid)] += 1
                assigned_pool = (
                    "gdna" if is_assigned_gdna else ("nrna" if is_assigned_nrna else "mrna")
                )
                pool_flow[(origin.kind, assigned_pool)] += 1
                # Vote the fragment's genomic locus (ZL) for any RNA endpoint.
                if comp_kind[t_cid] == "rna":
                    comp_zl[t_cid][zl] += 1
                if comp_kind[a_cid] == "rna":
                    comp_zl[a_cid][zl] += 1

        # Finalize RNA components' home locus to their modal ZL.
        for cid, counter in comp_zl.items():
            if counter:
                comp_locus[cid] = counter.most_common(1)[0][0]

        flows_by_cond[cond] = FlowData(
            condition=cond,
            flow=dict(flow),
            comp_name=comp_name,
            comp_kind=comp_kind,
            comp_locus=comp_locus,
            total_gdna_true=total_gdna,
            pool_flow=dict(pool_flow),
        )
        overview_rows.append(
            {
                "condition": cond,
                "total": total_rna + total_gdna,
                "total_rna": total_rna,
                "total_mrna": total_mrna,
                "total_nrna": total_nrna,
                "total_gdna": total_gdna,
                "correct_tx": correct_tx,
                "correct_gene": correct_gene,
                "wrong_tx": wrong_tx,
                "rna_as_gdna": rna_as_gdna,
                "mrna_as_gdna": mrna_as_gdna,
                "mrna_as_nrna": mrna_as_nrna,
                "nrna_as_rna": nrna_as_rna,
                "nrna_as_gdna": nrna_as_gdna,
                "gdna_as_rna": gdna_as_rna,
                "gdna_correct": gdna_correct,
            }
        )

    return flows_by_cond, overview_rows


def _spearman(x: "pd.Series", y: "pd.Series") -> float:
    """Robust Spearman that returns nan on degenerate input."""
    mask = x.notna() & y.notna()
    if mask.sum() < 5 or x[mask].nunique() < 2 or y[mask].nunique() < 2:
        return float("nan")
    from scipy.stats import spearmanr

    r, _ = spearmanr(x[mask], y[mask])
    return float(r)


@dataclass
class FlowData:
    """Per-condition fragment-flow matrix over integer-keyed components."""

    condition: str
    flow: dict[tuple[int, int], int]  # (true_cid, assigned_cid) -> fragment count
    comp_name: dict[int, str]  # cid -> transcript_id | "gdna@<locus>" | "unassigned"
    comp_kind: dict[int, str]  # cid -> "rna" | "gdna" | "unassigned"
    comp_locus: dict[int, int]  # cid -> locus_id (gdna@L -> L; unknown -> -1)
    total_gdna_true: int = 0  # true gDNA fragments (oracle), all loci + intergenic
    # 3-pool fragment flow (true_pool, assigned_pool) -> count, pools in {gdna, nrna, mrna}. The
    # net reduction (net(a→b)=flow[a,b]-flow[b,a]) cancels sequence-identical unrecoverable
    # misassignment; per-pool net surplus = assigned-true is the un-conflated leak/FP measure.
    pool_flow: dict[tuple[str, str], int] = field(default_factory=dict)


def _flow_marginals(flow: dict[tuple[int, int], int]) -> tuple[dict[int, int], dict[int, int]]:
    """Return (expected[c]=row sum from c, observed[c]=col sum into c)."""
    expected: Counter = Counter()
    observed: Counter = Counter()
    for (a, b), n in flow.items():
        expected[a] += n
        observed[b] += n
    return dict(expected), dict(observed)


def _pool_flow_3way_row(fd: FlowData) -> dict:
    """One per-condition 3-pool (gDNA/nRNA/mRNA) net-flow row.

    Reports each pool's true/assigned totals + **net surplus** (assigned − true; + = the pool is
    net-inflated, − = net-deficit) and the three **net fluxes** between pool pairs. Net cancels the
    sequence-identical, truly-unidentifiable misassignment, so the per-pool net surplus is the
    un-conflated leak / false-positive measure (a gross sum over-counts unrecoverable exchange).
    """
    pf = fd.pool_flow

    def cnt(a: str, b: str) -> int:
        return pf.get((a, b), 0)

    true = {p: sum(cnt(p, b) for b in _POOLS_3) for p in _POOLS_3}
    asg = {p: sum(cnt(a, p) for a in _POOLS_3) for p in _POOLS_3}
    row = {"condition": fd.condition}
    for p in _POOLS_3:
        row[f"{p}_true"] = true[p]
        row[f"{p}_assigned"] = asg[p]
        row[f"{p}_net_surplus"] = asg[p] - true[p]
    row["net_gdna_to_nrna"] = cnt("gdna", "nrna") - cnt("nrna", "gdna")
    row["net_gdna_to_mrna"] = cnt("gdna", "mrna") - cnt("mrna", "gdna")
    row["net_nrna_to_mrna"] = cnt("nrna", "mrna") - cnt("mrna", "nrna")
    return row


def _net_flow_rows(fd: FlowData) -> tuple[list[dict], list[dict]]:
    """Reduce one condition's flow matrix to per-transcript and per-locus net-flow rows.

    Per transcript T (home locus L):
      delta = observed - expected
            = net_from_gdna(gdna@L → T) + net_from_rna_isoforms(Σ T'≠T in L) + cross_locus
    """
    flow = fd.flow
    expected, observed = _flow_marginals(flow)

    def net(a: int, b: int) -> int:
        return flow.get((a, b), 0) - flow.get((b, a), 0)

    # Index components by home locus.
    rna_by_locus: dict[int, list[int]] = {}
    gdna_by_locus: dict[int, int] = {}
    for cid, kind in fd.comp_kind.items():
        loc = fd.comp_locus[cid]
        if kind == "rna":
            rna_by_locus.setdefault(loc, []).append(cid)
        elif kind == "gdna":
            gdna_by_locus[loc] = cid

    tx_rows: list[dict] = []
    for loc, rna_cids in rna_by_locus.items():
        gdna_cid = gdna_by_locus.get(loc)
        for cid in rna_cids:
            exp = expected.get(cid, 0)
            obs = observed.get(cid, 0)
            delta = obs - exp
            net_gdna = net(gdna_cid, cid) if gdna_cid is not None else 0
            net_iso = sum(net(other, cid) for other in rna_cids if other != cid)
            tx_rows.append(
                {
                    "condition": fd.condition,
                    "transcript_id": fd.comp_name[cid],
                    "locus_id": loc,
                    "expected": exp,
                    "observed": obs,
                    "delta": delta,
                    "net_from_gdna": net_gdna,
                    "net_from_rna_isoforms": net_iso,
                    "cross_locus": delta - net_gdna - net_iso,
                }
            )

    locus_rows: list[dict] = []
    all_loci = set(rna_by_locus) | set(gdna_by_locus)
    for loc in all_loci:
        rna_cids = rna_by_locus.get(loc, [])
        gdna_cid = gdna_by_locus.get(loc)
        rna_exp = sum(expected.get(c, 0) for c in rna_cids)
        rna_obs = sum(observed.get(c, 0) for c in rna_cids)
        gdna_exp = expected.get(gdna_cid, 0) if gdna_cid is not None else 0
        gdna_obs = observed.get(gdna_cid, 0) if gdna_cid is not None else 0
        locus_rows.append(
            {
                "condition": fd.condition,
                "locus_id": loc,
                "n_transcripts": len(rna_cids),
                "rna_expected": rna_exp,
                "rna_observed": rna_obs,
                "rna_delta": rna_obs - rna_exp,
                "gdna_expected": gdna_exp,
                "gdna_observed": gdna_obs,
                # + => gDNA mass leaked to RNA; - => RNA siphoned into gDNA.
                "net_gdna_to_rna": gdna_exp - gdna_obs,
            }
        )

    return tx_rows, locus_rows


def analyze_net_flow(
    sim_base: Path,
    conditions: list[str],
    flows: dict[str, FlowData] | None = None,
) -> str:
    """Primary gDNA-vs-RNA accuracy: net fragment-flow deconvolution.

    Writes ``net_flow_per_transcript.tsv`` and ``net_flow_per_locus.tsv`` (covariate-joined)
    and returns a text report: pool net-leak & direction, per-transcript Δ distribution,
    covariate root-cause ranking, and the identifiability diagnostic.
    """
    if flows is None:
        flows, _ = collect_fragment_flows(sim_base, conditions)

    cmeta = condition_manifest_map(load_manifest(sim_base))
    fallback_truth = load_truth(sim_base)

    tx_frames: list[pd.DataFrame] = []
    locus_frames: list[pd.DataFrame] = []
    for cond in conditions:
        fd = flows.get(cond)
        if fd is None:
            continue
        tx_rows, locus_rows = _net_flow_rows(fd)
        meta = cmeta.get(cond, {})
        tag = {
            "gdna_label": meta.get("gdna_label", parse_condition(cond)["gdna_label"]),
            "ss": meta.get("strand_specificity", parse_condition(cond)["strand_specificity"]),
            "capture": meta.get("capture_label", "off"),
        }

        if tx_rows:
            tdf = pd.DataFrame(tx_rows)
            # Transcript covariates from truth (n_exons, spliced_length, gene_id, strand).
            ctruth = load_condition_truth(sim_base, cond, cmeta, fallback_truth)
            cov_cols = [
                c
                for c in (
                    "transcript_id",
                    "gene_id",
                    "n_exons",
                    "spliced_length",
                    "strand",
                    "mrna_abundance",
                )
                if c in ctruth.columns
            ]
            tdf = tdf.merge(ctruth[cov_cols], on="transcript_id", how="left")
            if "n_exons" in tdf.columns:
                tdf["single_exon"] = tdf["n_exons"].fillna(0) <= 1
            for k, v in tag.items():
                tdf[k] = v
            tx_frames.append(tdf)

        if locus_rows:
            ldf = pd.DataFrame(locus_rows)
            loci_out = load_loci(sim_base / cond / "rigel_out")
            cov_cols = [
                c
                for c in (
                    "locus_id",
                    "locus_span_bp",
                    "n_em_fragments",
                    "gdna_prior_count",
                    "rna_prior_count",
                    "gdna_eff_len_em",
                )
                if c in loci_out.columns
            ]
            if cov_cols and not loci_out.empty:
                ldf = ldf.merge(loci_out[cov_cols], on="locus_id", how="left")
            for k, v in tag.items():
                ldf[k] = v
            locus_frames.append(ldf)

    boundaries: list[str] = []
    hr = "═" * 100
    boundaries.append(f"\n{hr}")
    boundaries.append("  NET FRAGMENT-FLOW DECONVOLUTION  (primary gDNA↔RNA accuracy)")
    boundaries.append(hr)
    boundaries.append(
        "  Net flow cancels symmetric (sequence-identical, unrecoverable) misassignment;\n"
        "  only systematic bias remains. + net_gdna_to_rna ⇒ gDNA leaked into RNA;\n"
        "  - ⇒ RNA siphoned into gDNA. Per transcript, Δ = net(gDNA→T) + net(RNA isoforms→T).\n"
    )

    if not tx_frames:
        boundaries.append("  [no annotated BAMs / loci available — run quant with --annotated-bam]")
        return "\n".join(boundaries)

    tx_all = pd.concat(tx_frames, ignore_index=True)
    locus_all = pd.concat(locus_frames, ignore_index=True) if locus_frames else pd.DataFrame()

    tx_path = sim_base / "net_flow_per_transcript.tsv"
    locus_path = sim_base / "net_flow_per_locus.tsv"
    tx_all.to_csv(tx_path, sep="\t", index=False)
    if not locus_all.empty:
        locus_all.to_csv(locus_path, sep="\t", index=False)

    # ── Pool net-leak & direction, per condition ──
    boundaries.append("  POOL NET LEAK (signed; fraction of true gDNA), by condition:")
    boundaries.append(
        f"  {'gdna':9}{'cap':5}{'ss':6} | {'true_gDNA':>10} {'net→RNA':>9} "
        f"{'leak%':>7} | {'in_locus':>9} {'intergenic':>11}"
    )
    boundaries.append("  " + "-" * 74)
    pool_rows = []
    for cond in conditions:
        fd = flows.get(cond)
        if fd is None:
            continue
        meta = cmeta.get(cond, {})
        sub = locus_all[locus_all["condition"] == cond] if not locus_all.empty else pd.DataFrame()
        in_locus_net = (
            int(sub[sub["locus_id"] >= 0]["net_gdna_to_rna"].sum()) if not sub.empty else 0
        )
        intergenic_net = (
            int(sub[sub["locus_id"] < 0]["net_gdna_to_rna"].sum()) if not sub.empty else 0
        )
        net_total = in_locus_net + intergenic_net
        true_gdna = fd.total_gdna_true
        leak_frac = net_total / true_gdna if true_gdna else 0.0
        pool_rows.append(
            {
                "gdna_label": meta.get("gdna_label", "?"),
                "gdna_rate": float(meta.get("gdna_rate", 0.0)),
                "capture": meta.get("capture_label", "off"),
                "ss": meta.get("strand_specificity", float("nan")),
                "true_gdna": true_gdna,
                "net_to_rna": net_total,
                "leak_frac": leak_frac,
                "in_locus": in_locus_net,
                "intergenic": intergenic_net,
            }
        )
    pool_rows.sort(key=lambda r: (r["gdna_rate"], r["capture"], -r["ss"]))
    for r in pool_rows:
        boundaries.append(
            f"  {r['gdna_label']:9}{r['capture']:5}{r['ss']:<6.2f} | {r['true_gdna']:>10,} "
            f"{r['net_to_rna']:>+9,} {r['leak_frac'] * 100:>6.2f}% | "
            f"{r['in_locus']:>+9,} {r['intergenic']:>+11,}"
        )
    cond_path = sim_base / "net_flow_per_condition.tsv"
    pd.DataFrame(pool_rows).to_csv(cond_path, sep="\t", index=False)

    # ── 3-pool net flow: gDNA ↔ nRNA ↔ mRNA (the resurrected nRNA pool needs 3 pools) ──
    boundaries.append(
        "\n  3-POOL NET FLOW (gDNA/nRNA/mRNA): per-pool net surplus (assigned−true) + net pair fluxes."
    )
    boundaries.append(
        "  + surplus ⇒ pool net-inflated (false gain); − ⇒ net-deficit. Watch mRNA surplus (mature FP)."
    )
    boundaries.append(
        f"  {'gdna':9}{'cap':4}{'ss':5}{'nrna':5} | {'gDNA_surp':>10}{'nRNA_surp':>10}{'mRNA_surp':>10}"
        f" | {'g→n':>8}{'g→m':>8}{'n→m':>8}"
    )
    boundaries.append("  " + "-" * 86)
    pool3_rows = []
    for cond in conditions:
        fd = flows.get(cond)
        if fd is None or not fd.pool_flow:
            continue
        meta = cmeta.get(cond, {})
        row = _pool_flow_3way_row(fd)
        row.update(
            {
                "gdna_label": meta.get("gdna_label", "?"),
                "gdna_rate": float(meta.get("gdna_rate", 0.0)),
                "capture": meta.get("capture_label", "off"),
                "ss": meta.get("strand_specificity", float("nan")),
                "nrna": meta.get("nrna_label", "rnd" if "nrna_rnd" in cond else "none"),
            }
        )
        pool3_rows.append(row)
    pool3_rows.sort(key=lambda r: (r["gdna_rate"], r["capture"], -r["ss"], r["nrna"]))
    for r in pool3_rows:
        boundaries.append(
            f"  {r['gdna_label']:9}{r['capture']:4}{r['ss']:<5.2f}{r['nrna']:5} | "
            f"{r['gdna_net_surplus']:>+10,}{r['nrna_net_surplus']:>+10,}{r['mrna_net_surplus']:>+10,}"
            f" | {r['net_gdna_to_nrna']:>+8,}{r['net_gdna_to_mrna']:>+8,}{r['net_nrna_to_mrna']:>+8,}"
        )
    if pool3_rows:
        pd.DataFrame(pool3_rows).to_csv(
            sim_base / "net_flow_3pool_per_condition.tsv", sep="\t", index=False
        )

    # ── Per-transcript Δ distribution + source decomposition ──
    boundaries.append("\n  PER-TRANSCRIPT Δ (observed − expected) and its decomposition:")
    boundaries.append(
        f"  {'gdna':9}{'cap':5}{'ss':6} | {'n_tx':>5} {'meanΔ':>7} {'sdΔ':>7} "
        f"{'mean|Δ|':>8} | {'fromGDNA':>9} {'fromISO':>8} | {'|Δ|>10':>6}"
    )
    boundaries.append("  " + "-" * 78)
    for (gl, cap, ss), grp in tx_all.groupby(["gdna_label", "capture", "ss"], dropna=False):
        boundaries.append(
            f"  {str(gl):9}{str(cap):5}{float(ss):<6.2f} | {len(grp):>5} "
            f"{grp['delta'].mean():>7.2f} {grp['delta'].std():>7.2f} "
            f"{grp['delta'].abs().mean():>8.2f} | "
            f"{grp['net_from_gdna'].mean():>+9.2f} {grp['net_from_rna_isoforms'].mean():>+8.2f} | "
            f"{int((grp['delta'].abs() > 10).sum()):>6}"
        )

    # ── Root cause: covariate ranking against gDNA contamination inflow ──
    boundaries.append(
        "\n  ROOT CAUSE — Spearman(net_from_gdna, covariate) over transcripts in gDNA>0 conditions:"
    )
    contam = tx_all[tx_all["gdna_label"] != "none"].copy()
    if not contam.empty:
        cov_candidates = {
            "true_RNA_depth(expected)": contam.get("expected"),
            "n_exons": contam.get("n_exons"),
            "spliced_length": contam.get("spliced_length"),
            "single_exon": contam.get("single_exon").astype(float)
            if "single_exon" in contam.columns
            else None,
            "mrna_abundance": contam.get("mrna_abundance"),
        }
        ranked = []
        for name, series in cov_candidates.items():
            if series is None:
                continue
            ranked.append((name, _spearman(contam["net_from_gdna"], series)))
        ranked = [r for r in ranked if not np.isnan(r[1])]
        ranked.sort(key=lambda r: abs(r[1]), reverse=True)
        for name, r in ranked:
            boundaries.append(f"    {name:<28} ρ = {r:+.3f}")
        if not ranked:
            boundaries.append("    (insufficient variation to rank covariates)")

    # ── Identifiability diagnostic: gross confusion vs net (expected-unrecoverable vs bias) ──
    boundaries.append("\n  IDENTIFIABILITY — single-exon (gDNA-identical) vs multi-exon transcripts:")
    if "single_exon" in contam.columns and not contam.empty:
        boundaries.append(
            f"    {'class':<12} {'n_tx':>6} {'mean|Δ|':>8} {'meanΔ':>8} {'mean net_from_gdna':>18}"
        )
        for label, mask in (
            ("single-exon", contam["single_exon"]),
            ("multi-exon", ~contam["single_exon"]),
        ):
            grp = contam[mask]
            if grp.empty:
                continue
            boundaries.append(
                f"    {label:<12} {len(grp):>6} {grp['delta'].abs().mean():>8.2f} "
                f"{grp['delta'].mean():>8.2f} {grp['net_from_gdna'].mean():>18.2f}"
            )
        boundaries.append(
            "    (single-exon transcripts are sequence-identical to gDNA and have no isoforms, so\n"
            "     their Δ ≈ net_from_gdna; compare gDNA inflow across classes to localize where\n"
            "     contamination concentrates. meanΔ near 0 with large spread ⇒ unbiased-but-noisy;\n"
            "     meanΔ ≈ mean net_from_gdna ⇒ systematic gDNA→RNA leak on that class.)"
        )

    boundaries.append(f"\n  Wrote {cond_path}")
    boundaries.append(f"  Wrote {tx_path}")
    if not locus_all.empty:
        boundaries.append(f"  Wrote {locus_path}")
    return "\n".join(boundaries)
