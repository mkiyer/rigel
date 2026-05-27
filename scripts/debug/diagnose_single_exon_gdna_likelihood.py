from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import digamma

from rigel.calibration import calibrate
from rigel.calibration.latent_states import STATE_NAMES
from rigel.calibration.signature import (
    CHANNEL_STRAND_NEG,
    CHANNEL_STRAND_POS,
    COMPARTMENT_BOUNDARY_LEFT,
    COMPARTMENT_BOUNDARY_RIGHT,
    COMPARTMENT_CONTAINED,
    SPLICE_SPLICED,
    SPLICE_UNSPLICED,
    channel_index,
)
from rigel.calibration.prior import assemble_priors
from rigel.config import BamScanConfig, EMConfig, FragmentScoringConfig, PipelineConfig
from rigel.locus import build_multi_loci
from rigel.locus_partition import partition_and_free
from rigel.pipeline import (
    _assign_locus_ids,
    _calibration_strand_summary,
    _native_detect_sj_tag,
    _score_fragments,
    _setup_geometry_and_estimator,
    scan_and_buffer,
)
from rigel.scoring import FragmentScorer, GDNA_SPLICE_PENALTIES, LOG_HALF, LOG_SAFE_FLOOR
from rigel.sim import ReadSimConfig, Scenario, run_benchmark
from rigel.splice import SPLICE_UNSPLICED, SpliceType
from rigel.types import Strand

SIM_SEED = 42
PIPELINE_SEED = 42
WORK_DIR = Path("/tmp/rigel_diagnose_single_exon_gdna_likelihood")
FRAGMENT_TSV = Path("/tmp/rigel_single_exon_gdna_likelihood_fragments.tsv")


def sim_config() -> ReadSimConfig:
    return ReadSimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=1.0,
        seed=SIM_SEED,
    )


def build_scenario():
    sc = Scenario("single_exon", genome_length=8000, seed=SIM_SEED, work_dir=WORK_DIR)
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(500, 1500)], "abundance": 100}])
    sc.add_gene(
        "g_helper",
        "+",
        [{"t_id": "t_helper", "exons": [(2500, 3000), (3500, 4000)], "abundance": 50}],
    )
    sc.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(6500, 6800)], "abundance": 0}])
    return sc.build_oracle(
        n_fragments=500,
        sim_config=sim_config(),
        gdna_config=None,
        nrna_abundance=0,
    )


def transcript_label(index, t_idx: int) -> str:
    row = index.t_df.iloc[int(t_idx)]
    for col in ("transcript_id", "t_id", "target_id"):
        if col in row.index:
            return str(row[col])
    return str(t_idx)


def transcript_index_by_id(index, transcript_id: str) -> int:
    for col in ("transcript_id", "t_id", "target_id"):
        if col in index.t_df.columns:
            hits = np.flatnonzero(index.t_df[col].astype(str).to_numpy() == transcript_id)
            if hits.size:
                return int(hits[0])
    raise RuntimeError(f"could not find transcript id {transcript_id!r}")


def snapshot_buffer_features(buffer) -> dict[int, dict[str, object]]:
    rows: dict[int, dict[str, object]] = {}
    for chunk in buffer.iter_chunks():
        for i in range(chunk.size):
            start = int(chunk.t_offsets[i])
            end = int(chunk.t_offsets[i + 1])
            t_indices = chunk.t_indices[start:end].astype(np.int32, copy=True)
            rows[int(chunk.frag_id[i])] = {
                "frag_id": int(chunk.frag_id[i]),
                "t_indices": t_indices,
                "frag_lengths": chunk.frag_lengths[start:end].astype(np.int32, copy=True),
                "exon_bp": chunk.exon_bp[start:end].astype(np.int32, copy=True),
                "splice_type": int(chunk.splice_type[i]),
                "exon_strand": int(chunk.exon_strand[i]),
                "sj_strand": int(chunk.sj_strand[i]),
                "fragment_class": int(chunk.fragment_classes[i]),
                "num_hits": int(chunk.num_hits[i]),
                "read_length": int(chunk.read_length[i]),
                "genomic_footprint": int(chunk.genomic_footprint[i]),
                "genomic_start": int(chunk.genomic_start[i]),
                "nm": int(chunk.nm[i]),
            }
    return rows


def log_likelihood(model, length: int) -> float:
    try:
        return float(model.log_likelihood(int(length)))
    except Exception:
        return float("nan")


def summarize(name: str, values: np.ndarray) -> str:
    arr = np.asarray(values, dtype=np.float64)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return f"{name}: no finite values"
    qs = np.percentile(finite, [0, 5, 25, 50, 75, 95, 100])
    return (
        f"{name}: n={finite.size} min={qs[0]:.4f} p05={qs[1]:.4f} "
        f"p25={qs[2]:.4f} p50={qs[3]:.4f} p75={qs[4]:.4f} "
        f"p95={qs[5]:.4f} max={qs[6]:.4f} mean={finite.mean():.4f}"
    )


def grouped_prior_update(raw_counts, carried_state, gdna_idx, alpha_gdna, alpha_rna, has_gdna):
    raw = np.maximum(np.nan_to_num(raw_counts, nan=0.0, posinf=0.0, neginf=0.0), 0.0)
    carried = None if carried_state is None else np.maximum(carried_state, 0.0)
    out = np.zeros_like(raw, dtype=np.float64)
    rna_mask = np.ones(raw.shape[0], dtype=bool)
    rna_mask[gdna_idx] = False
    n_rna = float(raw[rna_mask].sum())
    carried_rna = 0.0 if carried is None else float(carried[rna_mask].sum())
    a_g = float(alpha_gdna) if has_gdna else 0.0
    a_r = float(alpha_rna) if has_gdna else 0.0
    if n_rna <= 1.0e-12 and carried_rna <= 1.0e-12:
        a_r = 0.0
    n_g = float(raw[gdna_idx]) if has_gdna else 0.0
    out[gdna_idx] = n_g + a_g if has_gdna else 0.0
    rna_total = n_rna + a_r
    if n_rna > 1.0e-12:
        out[rna_mask] = rna_total * raw[rna_mask] / n_rna
    elif carried is not None and carried_rna > 1.0e-12:
        out[rna_mask] = rna_total * carried[rna_mask] / carried_rna
    return out


def build_local_arrays(part, locus, estimator, gdna_eff_len, alpha_gdna, alpha_rna):
    local_t = [int(x) for x in locus.transcript_indices]
    global_to_local = {t: i for i, t in enumerate(local_t)}
    gdna_idx = len(local_t)
    n_components = gdna_idx + 1
    unit_components: list[np.ndarray] = []
    unit_log_liks: list[np.ndarray] = []
    unit_cov_wts: list[np.ndarray] = []
    has_gdna = False
    for ui in range(part.n_units):
        comps = []
        lls = []
        wts = []
        start = int(part.offsets[ui])
        end = int(part.offsets[ui + 1])
        for j in range(start, end):
            t = int(part.t_indices[j])
            if t in global_to_local:
                comps.append(global_to_local[t])
                lls.append(float(part.log_liks[j]))
                wts.append(float(part.coverage_weights[j]))
        if int(part.is_spliced[ui]) == 0 and np.isfinite(float(part.gdna_log_liks[ui])):
            comps.append(gdna_idx)
            lls.append(float(part.gdna_log_liks[ui]))
            wts.append(1.0)
            has_gdna = True
        order = np.argsort(np.asarray(comps, dtype=np.int32))
        unit_components.append(np.asarray(comps, dtype=np.int32)[order])
        unit_log_liks.append(np.asarray(lls, dtype=np.float64)[order])
        unit_cov_wts.append(np.asarray(wts, dtype=np.float64)[order])
    log_eff = np.zeros(n_components, dtype=np.float64)
    for t, local in global_to_local.items():
        log_eff[local] = math.log(max(float(estimator.em_effective_lengths[t]), 1.0))
    log_eff[gdna_idx] = math.log(max(float(gdna_eff_len), 1.0))
    unambig = np.zeros(n_components, dtype=np.float64)
    for t, local in global_to_local.items():
        unambig[local] = float(estimator.unambig_counts[t].sum())
    warm_raw = unambig.copy()
    for comps, wts in zip(unit_components, unit_cov_wts, strict=True):
        total = float(wts.sum())
        if total <= 0.0:
            total = 1.0
        for comp, wt in zip(comps, wts, strict=True):
            warm_raw[int(comp)] += float(wt) / total
    alpha0 = grouped_prior_update(
        warm_raw,
        None,
        gdna_idx,
        alpha_gdna,
        alpha_rna,
        has_gdna,
    )
    return local_t, gdna_idx, unit_components, unit_log_liks, log_eff, unambig, alpha0, has_gdna


def vbem_fixed_point(
    unit_components,
    unit_log_liks,
    log_eff,
    unambig,
    alpha0,
    gdna_idx,
    alpha_gdna,
    alpha_rna,
    has_gdna,
    max_iter=10000,
    tol=1.0e-10,
):
    alpha = np.maximum(np.asarray(alpha0, dtype=np.float64), 1.0e-12)
    n_components = alpha.shape[0]
    post_by_unit: list[np.ndarray] = []
    for it in range(max_iter):
        dg_sum = float(digamma(max(float(alpha.sum()), 1.0e-12)))
        log_weights = digamma(np.maximum(alpha, 1.0e-12)) - dg_sum - log_eff
        em_totals = np.zeros(n_components, dtype=np.float64)
        post_by_unit = []
        for comps, lls in zip(unit_components, unit_log_liks, strict=True):
            vals = lls + log_weights[comps]
            maxv = float(np.max(vals))
            expv = np.exp(vals - maxv)
            denom = float(expv.sum())
            post = expv / denom if denom > 0.0 and np.isfinite(denom) else np.zeros_like(expv)
            post_by_unit.append(post)
            for comp, p in zip(comps, post, strict=True):
                em_totals[int(comp)] += float(p)
        raw_counts = unambig + em_totals
        new_alpha = grouped_prior_update(
            raw_counts,
            alpha,
            gdna_idx,
            alpha_gdna,
            alpha_rna,
            has_gdna,
        )
        new_alpha = np.maximum(new_alpha, 1.0e-12)
        old_theta = alpha / max(float(alpha.sum()), 1.0e-300)
        new_theta = new_alpha / max(float(new_alpha.sum()), 1.0e-300)
        delta = float(np.abs(new_theta - old_theta).sum())
        alpha = new_alpha
        if delta < tol:
            break
    theta = alpha / max(float(alpha.sum()), 1.0e-300)
    final_log_weights = np.log(theta + 1.0e-300) - log_eff
    final_post = []
    for comps, lls in zip(unit_components, unit_log_liks, strict=True):
        vals = lls + final_log_weights[comps]
        maxv = float(np.max(vals))
        expv = np.exp(vals - maxv)
        final_post.append(expv / float(expv.sum()))
    return theta, alpha, final_post, it + 1


def row_for_unit(
    index,
    scorer: FragmentScorer,
    calibration,
    feature_by_frag_id,
    part,
    local_t,
    gdna_idx,
    unit_components,
    unit_log_liks,
    log_eff,
    theta,
    final_post,
    ui: int,
):
    fid = int(part.frag_ids[ui])
    feature = feature_by_frag_id[fid]
    comps = unit_components[ui]
    lls = unit_log_liks[ui]
    posts = final_post[ui]
    rna_mask = comps != gdna_idx
    rna_positions = np.flatnonzero(rna_mask)
    best_pos = int(rna_positions[np.argmax(lls[rna_mask])])
    best_comp = int(comps[best_pos])
    best_t = int(local_t[best_comp])
    gdna_pos_arr = np.flatnonzero(comps == gdna_idx)
    gdna_pos = int(gdna_pos_arr[0]) if gdna_pos_arr.size else -1
    best_rna_ll = float(lls[best_pos])
    gdna_ll = float(lls[gdna_pos]) if gdna_pos >= 0 else float("-inf")
    best_rna_post = float(posts[best_pos])
    gdna_post = float(posts[gdna_pos]) if gdna_pos >= 0 else 0.0

    t_candidates = feature["t_indices"]
    t_match = np.flatnonzero(t_candidates == best_t)
    feature_pos = int(t_match[0]) if t_match.size else -1
    flen = int(feature["frag_lengths"][feature_pos]) if feature_pos >= 0 else -1
    exon_bp = int(feature["exon_bp"][feature_pos]) if feature_pos >= 0 else -1
    read_length = int(feature["read_length"])
    overhang = max(read_length - exon_bp, 0) if exon_bp >= 0 else -1
    stype = int(feature["splice_type"])
    exon_strand = int(feature["exon_strand"])
    nm = int(feature["nm"])
    same_strand = False
    if exon_strand in (int(Strand.POS), int(Strand.NEG)):
        same_strand = exon_strand == int(index.t_to_strand_arr[best_t])
        strand_term = scorer.log_p_sense if same_strand else scorer.log_p_antisense
    else:
        strand_term = LOG_HALF
    mismatch_term = nm * scorer.mismatch_log_penalty if nm > 0 else 0.0
    overhang_term = overhang * scorer.overhang_log_penalty if overhang >= 0 else 0.0
    rna_fl = log_likelihood(calibration.fl_models.rna, flen)
    global_fl = log_likelihood(calibration.fl_models.global_, flen)
    gfp = int(feature["genomic_footprint"])
    gdna_fl = log_likelihood(calibration.fl_models.gdna, gfp)
    gdna_splice_prob = GDNA_SPLICE_PENALTIES.get(stype, 1.0)
    gdna_splice_term = math.log(max(float(gdna_splice_prob), LOG_SAFE_FLOOR))
    calc_rna_ll = strand_term + rna_fl + overhang_term + mismatch_term
    calc_gdna_ll = gdna_fl + gdna_splice_term + LOG_HALF + mismatch_term
    raw_delta = gdna_ll - best_rna_ll
    fl_delta = gdna_fl - rna_fl
    strand_delta = LOG_HALF - strand_term
    overhang_delta = -overhang_term
    eff_delta = -float(log_eff[gdna_idx]) + float(log_eff[best_comp])
    theta_delta = math.log(theta[gdna_idx] + 1.0e-300) - math.log(theta[best_comp] + 1.0e-300)
    final_delta = raw_delta + eff_delta + theta_delta
    winner = "gDNA" if gdna_post > best_rna_post else transcript_label(index, best_t)
    return {
        "unit": int(ui),
        "frag_id": fid,
        "winner": winner,
        "gdna_post": gdna_post,
        "best_rna_post": best_rna_post,
        "raw_delta_gdna_minus_rna": raw_delta,
        "final_delta_with_theta_and_eff": final_delta,
        "theta_delta": theta_delta,
        "eff_delta": eff_delta,
        "fl_delta": fl_delta,
        "fl_delta_if_rna_used_global_fl": gdna_fl - global_fl,
        "raw_delta_if_rna_used_global_fl": gdna_fl - global_fl + strand_delta + overhang_delta,
        "strand_delta": strand_delta,
        "overhang_delta": overhang_delta,
        "best_t": best_t,
        "best_t_label": transcript_label(index, best_t),
        "best_rna_ll": best_rna_ll,
        "gdna_ll": gdna_ll,
        "calc_rna_ll": calc_rna_ll,
        "calc_gdna_ll": calc_gdna_ll,
        "calc_resid_rna": calc_rna_ll - best_rna_ll,
        "calc_resid_gdna": calc_gdna_ll - gdna_ll,
        "strand_term_rna": strand_term,
        "gdna_log_half": LOG_HALF,
        "rna_fl": rna_fl,
        "global_fl": global_fl,
        "gdna_fl": gdna_fl,
        "overhang": overhang,
        "overhang_term": overhang_term,
        "read_length": read_length,
        "exon_bp": exon_bp,
        "frag_length_rna": flen,
        "genomic_footprint": gfp,
        "genomic_start": int(feature["genomic_start"]),
        "genomic_midpoint": int(feature["genomic_start"]) + gfp // 2 if int(feature["genomic_start"]) >= 0 else -1,
        "exon_strand": exon_strand,
        "t_strand": int(index.t_to_strand_arr[best_t]),
        "same_strand": bool(same_strand),
        "splice_type": stype,
        "nm": nm,
        "count_col": int(part.locus_count_cols[ui]),
        "coverage_weight_best": float(part.coverage_weights[int(part.offsets[ui]) + best_pos])
        if int(part.offsets[ui]) + best_pos < int(part.offsets[ui + 1])
        else float("nan"),
    }


def unit_log_liks_with_global_rna_fl(
    unit_components,
    unit_log_liks,
    local_t,
    gdna_idx,
    part,
    feature_by_frag_id,
    calibration,
):
    adjusted = []
    for ui, (comps, lls) in enumerate(zip(unit_components, unit_log_liks, strict=True)):
        out = lls.copy()
        feature = feature_by_frag_id[int(part.frag_ids[ui])]
        t_candidates = feature["t_indices"]
        for pos, comp in enumerate(comps):
            if int(comp) == gdna_idx:
                continue
            t = int(local_t[int(comp)])
            t_match = np.flatnonzero(t_candidates == t)
            if not t_match.size:
                continue
            feature_pos = int(t_match[0])
            flen = int(feature["frag_lengths"][feature_pos])
            out[pos] += log_likelihood(calibration.fl_models.global_, flen) - log_likelihood(
                calibration.fl_models.rna,
                flen,
            )
        adjusted.append(out)
    return adjusted


def print_region_audit(index, payload, calibration, target_regions):
    region_df = index.region_df.copy()
    rc = calibration.region_calibration
    local = calibration.boundary_local
    sweep = calibration.boundary_sweep
    strand = calibration.strand_channels

    def channel(compartment, splice, strand_idx):
        return payload.region_counts[:, channel_index(compartment, splice, strand_idx)]

    contained_unspliced = channel(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS) + channel(
        COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG
    )
    contained_spliced = channel(COMPARTMENT_CONTAINED, SPLICE_SPLICED, CHANNEL_STRAND_POS) + channel(
        COMPARTMENT_CONTAINED, SPLICE_SPLICED, CHANNEL_STRAND_NEG
    )
    left_unspliced = channel(COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED, CHANNEL_STRAND_POS) + channel(
        COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG
    )
    right_unspliced = channel(
        COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED, CHANNEL_STRAND_POS
    ) + channel(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)
    left_spliced = channel(COMPARTMENT_BOUNDARY_LEFT, SPLICE_SPLICED, CHANNEL_STRAND_POS) + channel(
        COMPARTMENT_BOUNDARY_LEFT, SPLICE_SPLICED, CHANNEL_STRAND_NEG
    )
    right_spliced = channel(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_SPLICED, CHANNEL_STRAND_POS) + channel(
        COMPARTMENT_BOUNDARY_RIGHT, SPLICE_SPLICED, CHANNEL_STRAND_NEG
    )

    rows = []
    for rid in target_regions:
        row = region_df.iloc[int(rid)]
        rows.append(
            {
                "region_id": int(rid),
                "span": f"{int(row.start)}-{int(row.end)}",
                "signature": int(row.signature),
                "boundary_kind_left": int(row.boundary_kind_left),
                "boundary_kind_right": int(row.boundary_kind_right),
                "contained_unspliced": float(contained_unspliced[rid]),
                "contained_spliced": float(contained_spliced[rid]),
                "boundary_left_unspliced": float(left_unspliced[rid]),
                "boundary_right_unspliced": float(right_unspliced[rid]),
                "boundary_left_spliced": float(left_spliced[rid]),
                "boundary_right_spliced": float(right_spliced[rid]),
                "prior_unspliced_total": float(rc.prior_mass.unspliced_total[rid]),
                "prior_gdna_mean": float(rc.prior_mass.gdna_unspliced_mean[rid]),
                "prior_rna_mean": float(rc.prior_mass.rna_unspliced_mean[rid]),
                "mu_gdna": float(rc.mu_gdna[rid]),
                "upper_gdna": float(rc.upper_gdna[rid]),
                "rna_lower": float(rc.rna_lower[rid]),
                "A_r": float(rc.A_r[rid]),
                "gamma_r": float(rc.gamma_r[rid]),
                "local_alpha_excess": float(local.alpha_excess[rid]),
                "local_mu": float(local.mu_local[rid]),
                "sweep_alpha_excess": float(sweep.alpha_excess[rid]),
                "sweep_mu": float(sweep.mu_sweep[rid]),
                "p_background": float(rc.p_states[rid, 0]),
                "p_gdna_only_capture": float(rc.p_states[rid, 1]),
                "p_expressed_capture": float(rc.p_states[rid, 2]),
                "p_expressed_offtarget": float(rc.p_states[rid, 3]),
                "strand_contained_mean": float(strand.contained_mean[rid]) if strand is not None else float("nan"),
                "strand_boundary_left_mean": float(strand.boundary_left_mean[rid])
                if strand is not None
                else float("nan"),
                "strand_boundary_right_mean": float(strand.boundary_right_mean[rid])
                if strand is not None
                else float("nan"),
            }
        )
    audit = pd.DataFrame(rows)
    print("\n=== Region / boundary audit around t1 ===")
    print(audit.to_string(index=False))
    print("\nState names:", STATE_NAMES)


def main() -> None:
    result = build_scenario()
    index = result.index
    sample_config = PipelineConfig(
        em=EMConfig(seed=PIPELINE_SEED, assignment_mode="sample", rna_call_bias=0.5),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    from rigel.pipeline import run_pipeline

    sample_pr = run_pipeline(result.bam_path, index, config=sample_config)
    bench = run_benchmark(result, sample_pr, scenario_name="diagnose_single_exon_sample")
    print("=== Reproduction with default sample assignment ===")
    print(bench.summary())
    print(sample_pr.estimator.get_loci_df(index)[["locus_id", "n_em_fragments", "mrna", "nrna", "gdna", "alpha_gdna_add", "alpha_rna_add", "prior_rna_share_final", "enable_gdna"]].to_string(index=False))

    em_config = EMConfig(seed=PIPELINE_SEED, assignment_mode="fractional", rna_call_bias=0.5)
    scan = BamScanConfig(sj_strand_tag="auto")
    if scan.sj_strand_tag == "auto":
        scan = replace(scan, sj_strand_tag=_native_detect_sj_tag(str(result.bam_path)))

    stats, strand_models, frag_length_models, buffer, payload = scan_and_buffer(
        str(result.bam_path), index, scan
    )
    strand_models.finalize()
    feature_by_frag_id = snapshot_buffer_features(buffer)
    strand_summary = _calibration_strand_summary(strand_models)
    calibration = calibrate(
        index=index,
        payload=payload,
        scan_trained=frag_length_models,
        strand_summary=strand_summary,
    )

    geometry, estimator = _setup_geometry_and_estimator(index, calibration.fl_models.rna, em_config)
    scoring = FragmentScoringConfig()
    scorer = FragmentScorer.from_models(
        strand_models,
        calibration.fl_models.rna,
        calibration.fl_models.gdna,
        index,
        estimator,
        overhang_log_penalty=scoring.overhang_log_penalty,
        mismatch_log_penalty=scoring.mismatch_log_penalty,
        gdna_splice_penalties=scoring.gdna_splice_penalties,
        pruning_min_posterior=scoring.pruning_min_posterior,
    )
    em_data = _score_fragments(
        buffer,
        index,
        strand_models,
        calibration.fl_models.rna,
        calibration.fl_models.gdna,
        stats,
        estimator,
        scoring,
        scan.log_every,
        None,
    )
    multi_loci = build_multi_loci(em_data, index)
    _assign_locus_ids(estimator, multi_loci)
    prior_table = assemble_priors(
        multi_loci=multi_loci,
        em_data=em_data,
        index=index,
        calibration=calibration,
        em_config=em_config,
    )
    t1_idx = transcript_index_by_id(index, "t1")
    t1_locus = next(locus for locus in multi_loci if t1_idx in set(map(int, locus.transcript_indices)))
    lid = int(t1_locus.multi_locus_id)
    print("\n=== t1 locus identity ===")
    print("t1_idx", t1_idx, "locus_id", lid)
    print(
        "locus transcript indices:",
        [(int(t), transcript_label(index, int(t))) for t in t1_locus.transcript_indices],
    )
    print(
        "prior:",
        {
            "alpha_gdna_add": float(prior_table.alpha_gdna_add[lid]),
            "alpha_rna_add": float(prior_table.alpha_rna_add[lid]),
            "prior_rna_share_final": float(prior_table.prior_rna_share_final[lid]),
            "gdna_eff_len": float(prior_table.gdna_eff_len[lid]),
            "gdna_eff_len_unweighted": float(prior_table.gdna_eff_len_unweighted[lid]),
            "gdna_em_exposure_weight": float(prior_table.gdna_em_exposure_weight[lid]),
            "enable_gdna": int(prior_table.enable_gdna[lid]),
            "prior_n_local_gdna": float(prior_table.prior_n_local_gdna[lid]),
            "prior_n_local_rna": float(prior_table.prior_n_local_rna[lid]),
        },
    )

    global_midpoint = em_data.genomic_midpoint.copy()
    global_offsets = em_data.offsets.copy()
    global_t_indices = em_data.t_indices.copy()
    global_log_liks = em_data.log_liks.copy()
    global_gdna_liks = em_data.gdna_log_liks.copy()
    global_frag_ids = em_data.frag_ids.copy()
    print(
        "global raw score sanity:",
        {
            "n_units": int(em_data.n_units),
            "finite_gdna_units": int(np.isfinite(global_gdna_liks).sum()),
            "gDNA_raw_gt_best_RNA_all_units": int(
                sum(
                    float(global_gdna_liks[u])
                    > float(np.max(global_log_liks[global_offsets[u] : global_offsets[u + 1]]))
                    for u in range(em_data.n_units)
                    if np.isfinite(float(global_gdna_liks[u]))
                )
            ),
        },
    )

    partitions = partition_and_free(em_data, multi_loci)
    part = partitions[lid]
    partition_tuple = (
        part.offsets,
        part.t_indices,
        part.log_liks,
        part.coverage_weights,
        part.count_cols,
        part.is_spliced,
        part.gdna_log_liks,
        part.locus_t_indices,
        part.locus_count_cols,
    )
    native_result = estimator.run_batch_locus_em_partitioned(
        [partition_tuple],
        [t1_locus.transcript_indices],
        np.array([prior_table.alpha_gdna_add[lid]], dtype=np.float64),
        index,
        rna_prior_count=np.array([prior_table.alpha_rna_add[lid]], dtype=np.float64),
        gdna_eff_len=np.array([prior_table.gdna_eff_len[lid]], dtype=np.float64),
        enable_gdna=np.array([prior_table.enable_gdna[lid]], dtype=np.uint8),
        em_iterations=em_config.iterations,
        em_convergence_delta=em_config.convergence_delta,
        emit_assignments=True,
    )
    print("\n=== Native fractional EM for t1 locus ===")
    print(
        {
            "native_locus_rna": float(native_result[1][0]),
            "native_locus_gdna": float(native_result[2][0]),
            "native_gdna_fraction": float(native_result[2][0] / part.n_units),
            "native_map_winner_gdna_count": int(np.count_nonzero(native_result[3] == -2)),
            "native_mean_winner_posterior": float(np.asarray(native_result[4], dtype=np.float64).mean()),
        }
    )

    (
        local_t,
        gdna_idx,
        unit_components,
        unit_log_liks,
        log_eff,
        unambig,
        alpha0,
        has_gdna,
    ) = build_local_arrays(
        part,
        t1_locus,
        estimator,
        prior_table.gdna_eff_len[lid],
        prior_table.alpha_gdna_add[lid],
        prior_table.alpha_rna_add[lid],
    )
    theta, alpha, final_post, n_iter = vbem_fixed_point(
        unit_components,
        unit_log_liks,
        log_eff,
        unambig,
        alpha0,
        gdna_idx,
        prior_table.alpha_gdna_add[lid],
        prior_table.alpha_rna_add[lid],
        has_gdna,
    )
    global_rna_log_liks = unit_log_liks_with_global_rna_fl(
        unit_components,
        unit_log_liks,
        local_t,
        gdna_idx,
        part,
        feature_by_frag_id,
        calibration,
    )
    theta_global, alpha_global, final_post_global, n_iter_global = vbem_fixed_point(
        unit_components,
        global_rna_log_liks,
        log_eff,
        unambig,
        alpha0,
        gdna_idx,
        prior_table.alpha_gdna_add[lid],
        prior_table.alpha_rna_add[lid],
        has_gdna,
    )
    global_gdna_post_sum = float(
        sum(
            post[np.flatnonzero(comps == gdna_idx)[0]]
            for comps, post in zip(unit_components, final_post_global, strict=True)
            if np.any(comps == gdna_idx)
        )
    )
    rows = [
        row_for_unit(
            index,
            scorer,
            calibration,
            feature_by_frag_id,
            part,
            local_t,
            gdna_idx,
            unit_components,
            unit_log_liks,
            log_eff,
            theta,
            final_post,
            ui,
        )
        for ui in range(part.n_units)
    ]
    df = pd.DataFrame(rows)
    df.to_csv(FRAGMENT_TSV, sep="\t", index=False)

    print("\n=== Python VBEM check for t1 locus ===")
    print("iterations", n_iter)
    print("theta", {"gDNA": float(theta[gdna_idx]), **{transcript_label(index, t): float(theta[i]) for i, t in enumerate(local_t)}})
    print("alpha", {"gDNA": float(alpha[gdna_idx]), **{transcript_label(index, t): float(alpha[i]) for i, t in enumerate(local_t)}})
    print("sum p_gdna", float(df["gdna_post"].sum()), "mean", float(df["gdna_post"].mean()))
    print(
        "counterfactual RNA=global FL",
        {
            "iterations": int(n_iter_global),
            "theta_gdna": float(theta_global[gdna_idx]),
            "alpha_gdna": float(alpha_global[gdna_idx]),
            "sum_p_gdna": global_gdna_post_sum,
            "mean_p_gdna": global_gdna_post_sum / float(part.n_units),
        },
    )
    print("fragment table", str(FRAGMENT_TSV))

    print("\n=== Score deltas: gDNA minus best RNA ===")
    for col in (
        "raw_delta_gdna_minus_rna",
        "final_delta_with_theta_and_eff",
        "fl_delta",
        "fl_delta_if_rna_used_global_fl",
        "raw_delta_if_rna_used_global_fl",
        "strand_delta",
        "overhang_delta",
        "eff_delta",
        "theta_delta",
        "gdna_post",
    ):
        print(summarize(col, df[col].to_numpy()))
    print(
        "counts",
        {
            "n_units": int(len(df)),
            "raw_gdna_gt_rna": int((df["raw_delta_gdna_minus_rna"] > 0.0).sum()),
            "final_gdna_posterior_gt_0.5": int((df["gdna_post"] > 0.5).sum()),
            "overhang_positive": int((df["overhang"] > 0).sum()),
            "same_strand": int(df["same_strand"].sum()),
            "opposite_or_unknown_strand": int((~df["same_strand"]).sum()),
            "calc_rna_abs_resid_max": float(df["calc_resid_rna"].abs().max()),
            "calc_gdna_abs_resid_max": float(df["calc_resid_gdna"].abs().max()),
        },
    )
    show_cols = [
        "unit",
        "frag_id",
        "winner",
        "gdna_post",
        "raw_delta_gdna_minus_rna",
        "fl_delta",
        "strand_delta",
        "overhang_delta",
        "eff_delta",
        "theta_delta",
        "frag_length_rna",
        "genomic_footprint",
        "overhang",
        "same_strand",
        "genomic_start",
    ]
    print("\nTop 20 gDNA-favored fragments by raw score:")
    print(df.sort_values("raw_delta_gdna_minus_rna", ascending=False)[show_cols].head(20).to_string(index=False))
    print("\nTop 20 RNA-favored fragments by raw score:")
    print(df.sort_values("raw_delta_gdna_minus_rna", ascending=True)[show_cols].head(20).to_string(index=False))

    t1_midpoints = global_midpoint[t1_locus.unit_indices]
    print(
        "\nmidpoint sanity",
        {
            "min": int(np.min(t1_midpoints)),
            "max": int(np.max(t1_midpoints)),
            "outside_t1_exon": int(np.count_nonzero((t1_midpoints < 500) | (t1_midpoints >= 1500))),
        },
    )
    target_regions = [0, 1, 2]
    print_region_audit(index, payload, calibration, target_regions)
    print("\nFL model summary:")
    print(calibration.fl_models.to_summary_dict())


if __name__ == "__main__":
    main()
