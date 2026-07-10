"""Canonical calibration/EM measurement helpers — the ONE place the nascent-siphon metric is defined.

This module exists because the nascent-RNA siphon was mis-measured THREE different ways during diagnosis,
each producing a different wrong answer and eroding trust in the results (see
``docs/calibration/siphon_measurement.md``). Every diagnostic that reports the siphon MUST import from here
rather than re-deriving it.

THE THREE TRAPS (all give WRONG siphon numbers):
  1. ``estimator.get_counts_df(index)["count"]`` is ``np.where(is_synthetic, 0.0, t_total)`` — it ZEROES the
     synthetic nascent shadow rows, so summing it over nascent rows returns 0 ("the toy never siphons").
  2. Filtering by ``is_nrna`` CONFLATES real annotated SINGLE-EXON transcripts (legitimately
     ``is_nrna=True`` — a single-exon transcript is both mature and nascent) with the synthetic shadow
     spans. Their mass is REAL, not siphon; including it over-counts (the false "structural residual").
  3. ``net_flow`` ZF-bit (0x08) and ZT-tag hard labels (``rigel.sim.analysis``) are per-fragment hard MAP
     labels, NOT the EM abundance — insensitive to the soft posterior and not comparable across runs.

THE TRUTH: the siphon is the EM-assigned mass on the SYNTHETIC nascent shadow spans (the unspliced twins of
MULTI-exon transcripts, which have no independent molecular existence). On an ``nrna_none`` scenario (true
nascent = 0) this equals the pure siphon. It is read from ``estimator.em_counts`` (indexed over ALL
transcripts, incl. the synthetic rows that ``get_counts_df`` drops), NEVER from a display/hard-label column.
"""
import numpy as np


def em_counts_per_transcript(estimator):
    """EM-assigned fragment count per transcript (summed over channels), indexed over ALL transcripts
    including the synthetic nascent shadows. The faithful basis for every RNA-component measurement."""
    return np.asarray(estimator.em_counts, dtype=np.float64).sum(axis=1)


def nascent_siphon(estimator, index=None):
    """The nascent-RNA siphon = total count on the SYNTHETIC nascent shadow spans.

    Delegates to the estimator's own canonical scalar ``estimator.nrna_em_count`` (``t_counts`` summed over
    the ``_synthetic_mask``) — the PRODUCTION code always had the right metric; the traps were in tools that
    re-derived it. This is the ONLY correct siphon measure (module docstring lists the three traps). On an
    ``nrna_none`` scenario (true nascent abundance = 0) any nonzero value is pure siphon. ``index`` is
    accepted for signature symmetry but unused when the property is available."""
    prop = getattr(estimator, "nrna_em_count", None)
    if prop is not None:
        return float(prop)
    em = em_counts_per_transcript(estimator)  # fallback: synthetic mask over EM counts
    syn = index.t_df["is_synthetic"].to_numpy(bool)
    return float(em[syn].sum())


def rna_component_breakdown(estimator, index):
    """(siphon, single_exon_nrna, multiexon_mature) EM mass — the three RNA-component buckets, so a report
    never again conflates the legitimate single-exon ``is_nrna`` mass with the synthetic-shadow siphon."""
    em = em_counts_per_transcript(estimator)
    syn = index.t_df["is_synthetic"].to_numpy(bool)
    is_nrna = index.t_df["is_nrna"].to_numpy(bool)
    single_exon_nrna = is_nrna & ~syn  # real single-exon transcripts (both mature+nascent) — NOT siphon
    mature = ~is_nrna & ~syn
    return float(em[syn].sum()), float(em[single_exon_nrna].sum()), float(em[mature].sum())


def ontarget_gdna_fl_pmf(bam_path, region_arrays, index, max_size):
    """The TRUE on-target (captured, exonic-region) gDNA fragment-length pmf from a sim oracle BAM — the
    length distribution of the gDNA that actually competes with RNA in loci. The production gDNA FL model
    trains on the OFF-target (intergenic+intronic) accumulator pools (see calibration.fl.gdna_fl_mass),
    which capture leaves LONGER than the size-selected on-target gDNA; feeding this on-target pmf to
    calibrate()+scorer tests the capture-aware fix. Returns a normalized pmf aligned to max_size+1."""
    import pysam
    from rigel.sim.read_name import parse_origin
    from rigel.calibration.signature import BIT_EXON_POS, BIT_EXON_NEG

    starts = np.asarray(region_arrays.start, np.int64)
    ends = np.asarray(region_arrays.end, np.int64)
    ref_off = np.asarray(region_arrays.ref_offsets, np.int64)
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    name2id = index.ref_name_to_id
    exon_bit = BIT_EXON_POS | BIT_EXON_NEG
    counts = np.zeros(int(max_size) + 1, dtype=np.float64)
    with pysam.AlignmentFile(str(bam_path), "rb") as f:
        default = f.references[0] if f.references else None
        seen = set()
        for r in f:
            q = r.query_name
            if q in seen:
                continue
            seen.add(q)
            o = parse_origin(q)
            if o.kind != "gdna" or o.start is None:
                continue
            rid = name2id.get(str(o.ref if o.ref is not None else default))
            if rid is None:
                continue
            lo0, hi0 = int(ref_off[rid]), int(ref_off[rid + 1])
            mid = (int(o.start) + int(o.end)) // 2
            j = lo0 + int(np.searchsorted(ends[lo0:hi0], mid, side="right"))
            if not (lo0 <= j < hi0) or not (sig[j] & exon_bit):
                continue
            counts[min(int(o.end) - int(o.start), int(max_size))] += 1.0
    total = counts.sum()
    return counts / total if total > 0 else counts


def oracle_node_masses(bam_path, region_arrays, index):
    """Build the TRUE per-node gDNA/RNA masses from a sim oracle BAM's read-name origins — the "oracle
    calibration" used to separate calibration error from structural EM effects (feed to
    ``dataclasses.replace(cal, **oracle_node_masses(...))``). Contained fragments go to their region;
    crossing fragments split 0.5/side across the interior boundaries they span; spliced RNA is tallied
    separately. Returns the dict of mass arrays keyed exactly as ``CalibrationResult`` expects."""
    import pysam
    from rigel.sim.read_name import parse_origin

    starts = np.asarray(region_arrays.start, np.int64)
    ends = np.asarray(region_arrays.end, np.int64)
    ref_off = np.asarray(region_arrays.ref_offsets, np.int64)
    ref_id = np.asarray(region_arrays.ref_id)
    name_to_id = index.ref_name_to_id
    n = region_arrays.n_regions
    g_c, r_c = np.zeros(n), np.zeros(n)
    g_l, r_l = np.zeros(n), np.zeros(n)
    g_r, r_rt = np.zeros(n), np.zeros(n)
    r_spl = np.zeros(n)
    spliced_q = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as f:
        for rd in f:
            q = rd.query_name
            sp = rd.cigartuples is not None and any(op == 3 for op, _ in rd.cigartuples)
            spliced_q[q] = spliced_q.get(q, False) or sp
    with pysam.AlignmentFile(str(bam_path), "rb") as f:
        default_ref = f.references[0] if f.references else None
        seen = set()
        for rd in f:
            q = rd.query_name
            if q in seen:
                continue
            seen.add(q)
            o = parse_origin(q)
            if o.start is None:
                continue
            rid = name_to_id.get(str(o.ref if o.ref is not None else default_ref))
            if rid is None:
                continue
            lo0, hi0 = int(ref_off[rid]), int(ref_off[rid + 1])
            a, b = int(o.start), int(o.end)
            is_g = o.kind == "gdna"
            lo = lo0 + int(np.searchsorted(ends[lo0:hi0], a, side="right"))
            hi = lo0 + int(np.searchsorted(starts[lo0:hi0], b, side="left"))
            if hi <= lo:
                continue
            if hi - lo == 1:
                (g_c if is_g else r_c)[lo] += 1.0
                if not is_g and spliced_q.get(q, False):
                    r_spl[lo] += 1.0
            else:
                nb = hi - 1 - lo
                w = 0.5 / nb
                for rr in range(lo, hi - 1):
                    if ref_id[rr] != ref_id[rr + 1]:
                        continue
                    (g_r if is_g else r_rt)[rr] += w
                    (g_l if is_g else r_l)[rr + 1] += w
                if not is_g and spliced_q.get(q, False):
                    r_spl[lo] += 1.0
    return dict(mass_gdna_contained=g_c, mass_rna_contained=r_c, mass_gdna_left=g_l,
                mass_rna_left=r_l, mass_gdna_right=g_r, mass_rna_right=r_rt, mass_rna_spliced=r_spl)
