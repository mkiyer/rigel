"""THE calibration oracle — the single, validated ground-truth source for debugging/testing/benchmarking.

Principle: the oracle IS the production accumulator, partitioned by TRUE fragment origin. We split the sim
BAM into gdna / mrna / nrna by read-name origin (:func:`rigel.sim.read_name.parse_origin`), run the SAME
production scanner+accumulator on each partition, and assert the partitions sum to the full payload
(byte-exact on the integer channels; float32-rounding tolerance on boundary mass). Because the accumulator
deposits each fragment independently, this sum-to-full identity PROVES the partition is the production
payload split by origin — no reimplementation, nothing to get subtly wrong.

This replaces the retired ``oracle_node_masses`` (in the deleted ``_metrics``/``oracle_*`` scripts), which
deposited WHOLE fragments by SPAN with no intron-cutting — an INCOMPATIBLE basis with the accumulator the
calibration actually consumes (per-base coverage, introns cut). That mismatch (e.g. it reported 0 RNA in
high-expression exons where the accumulator has the real unspliced exon-body mRNA) confounded earlier
"calibration error" conclusions. See docs/calibration/oracle_and_benchmarking.md.

Accumulator channels (region_contained[R,4] and boundary_mass_{left,right}[B,4]):
  ch0 = unspliced genome+   ch1 = unspliced genome−   ch2 = spliced sense   ch3 = spliced antisense
The gDNA-vs-RNA deconvolution (calibration) is over the UNSPLICED channels (0,1); spliced (2,3) is mature
RNA that never competes with gDNA (result.py withholds it from the RNA prior). gDNA fragments are never
spliced (validated: ch2/3 == 0 in the gdna partition).

    OMP_NUM_THREADS=1 python scripts/debug/oracle.py [condition] [--suite DIR]
"""
from __future__ import annotations
import argparse
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pysam

from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.sim.read_name import parse_origin

ORIGINS = ("gdna", "mrna", "nrna")
_UNSPL = (0, 1)   # unspliced channels (the gDNA-vs-RNA competition basis)
_SPL = (2, 3)     # spliced channels (mature RNA; never competes with gDNA)


def _split_bam(bam: str, out_dir: Path, tag: str) -> tuple[dict[str, str], dict[str, int]]:
    """Split a name-sorted BAM into per-origin BAMs. Both mates share a qname ⇒ same origin ⇒ same file;
    iteration order (hence name-sort) is preserved. EVERY read is written to exactly one partition — the
    'account for every fragment' guarantee — and the total is asserted to reconcile."""
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {k: str(out_dir / f"{tag}.{k}.bam") for k in ORIGINS}
    counts = {k: 0 for k in ORIGINS}
    n_in = 0
    with pysam.AlignmentFile(bam, "rb") as fin:
        w = {k: pysam.AlignmentFile(paths[k], "wb", template=fin) for k in ORIGINS}
        for r in fin:
            n_in += 1
            k = parse_origin(r.query_name).kind  # raises on any unclassifiable read (no silent drop)
            w[k].write(r)
            counts[k] += 1
        for x in w.values():
            x.close()
    if sum(counts.values()) != n_in:
        raise AssertionError(f"oracle split dropped reads: in={n_in} out={sum(counts.values())}")
    return paths, counts


def _scan_payload(bam: str, index, cfg):
    from dataclasses import replace as dc
    sc = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _stats, _sm, _flm, _buf, payload = scan_and_buffer(bam, index, sc)
    return payload


@dataclass
class OracleTruth:
    """Validated per-origin accumulator payloads for one condition. Construct via :meth:`from_bam`, which
    runs the sum-to-full validation as a HARD gate (raises if the partition does not reconstruct the full
    payload — i.e. if the oracle is ever not trustworthy)."""
    full: object
    parts: dict          # origin -> payload
    read_counts: dict    # origin -> reads written (every input read accounted for)
    boundary_mass_tol: float

    @classmethod
    def from_bam(cls, bam: str, index, cfg, work_dir: Path, tag: str,
                 boundary_mass_tol: float = 1e-2, full_payload=None) -> "OracleTruth":
        """Split the BAM by origin, scan each partition, and validate sum-to-full.

        ``full_payload`` lets a caller that has ALREADY scanned the full BAM (e.g. to run
        ``calibrate`` on it) hand that payload in, skipping a redundant full re-scan. It must be
        the production scan of ``bam`` with the same ``cfg`` — sum-to-full then also PROVES the
        oracle partitions reconstruct the exact payload the calibration consumed."""
        paths, read_counts = _split_bam(bam, work_dir, tag)
        full = full_payload if full_payload is not None else _scan_payload(bam, index, cfg)
        parts = {k: _scan_payload(paths[k], index, cfg) for k in ORIGINS}
        self = cls(full=full, parts=parts, read_counts=read_counts, boundary_mass_tol=boundary_mass_tol)
        self._validate()
        return self

    def _validate(self) -> None:
        rc_full = np.asarray(self.full.region_contained, np.int64)
        rc_sum = sum(np.asarray(self.parts[k].region_contained, np.int64) for k in ORIGINS)
        if not np.array_equal(rc_sum, rc_full):
            raise AssertionError(
                f"oracle INVALID: region_contained partitions do not sum to full "
                f"(max|diff|={np.abs(rc_sum - rc_full).max()}). The partition is not the production split.")
        for arr in ("boundary_flux_left", "boundary_flux_right"):
            af = np.asarray(getattr(self.full, arr), np.int64)
            asum = sum(np.asarray(getattr(self.parts[k], arr), np.int64) for k in ORIGINS)
            if not np.array_equal(asum, af):
                raise AssertionError(f"oracle INVALID: {arr} partitions do not sum to full.")
        for arr in ("boundary_mass_left", "boundary_mass_right"):
            af = np.asarray(getattr(self.full, arr), np.float64)
            asum = sum(np.asarray(getattr(self.parts[k], arr), np.float64) for k in ORIGINS)
            md = float(np.abs(asum - af).max())
            if md > self.boundary_mass_tol:
                raise AssertionError(f"oracle INVALID: {arr} sum-to-full maxdiff {md:.3e} > tol.")
        # gDNA is never spliced (physical): the gdna partition must have zero spliced contained mass.
        g_spl = np.asarray(self.parts["gdna"].region_contained, np.int64)[:, _SPL].sum()
        if g_spl != 0:
            raise AssertionError(f"oracle INVALID: gdna partition has {g_spl} spliced contained reads (>0).")

    # ---- per-region TRUE masses on the accumulator basis ----
    def region_unspliced(self):
        """(G, R) per region: TRUE unspliced gDNA vs unspliced RNA contained mass — the gDNA-vs-RNA
        competition basis the calibration deconvolves. R = mrna+nrna unspliced (exon-body + nascent)."""
        rc = lambda k: np.asarray(self.parts[k].region_contained, np.float64)  # noqa: E731
        G = rc("gdna")[:, _UNSPL].sum(1)
        R = (rc("mrna") + rc("nrna"))[:, _UNSPL].sum(1)
        return G, R

    def region_true_fg(self):
        """Per-region TRUE gDNA fraction of the unspliced contained mass (NaN where no unspliced mass)."""
        G, R = self.region_unspliced()
        tot = G + R
        return np.where(tot > 0, G / np.maximum(tot, 1e-12), np.nan), tot

    def region_pools(self) -> dict:
        """Per-region TRUE contained mass on the accumulator basis, split by ORIGIN × genome STRAND.

        The gDNA-vs-RNA calibration deconvolves the **unspliced** channels into the 2-simplex
        ``(RNA₊, RNA₋, gDNA)`` (genome strand: ch0=+, ch1=−). Spliced (ch2 sense, ch3 antisense) is
        guaranteed-mature RNA that never competes with gDNA. Calibration cannot split mature from
        nascent — that is the downstream EM's job — so ``mature`` / ``nascent`` here are the TRUE
        composition of the RNA the calibration lumps together. Every array is float64[R]; all eight
        components sum (over origins × channels) to the full per-region contained mass (the validated
        sum-to-full identity). Keys:
          gdna_pos/gdna_neg          — unspliced gDNA by genome strand (should be ~50/50)
          mat_uns_pos/mat_uns_neg    — unspliced mature (exon-body) RNA by genome strand
          nas_uns_pos/nas_uns_neg    — unspliced nascent RNA by genome strand
          mat_spl/nas_spl            — spliced RNA (mature / nascent), guaranteed-RNA, no gDNA rival
        """
        rc = lambda k: np.asarray(self.parts[k].region_contained, np.float64)  # noqa: E731
        g, m, n = rc("gdna"), rc("mrna"), rc("nrna")
        return dict(
            gdna_pos=g[:, 0], gdna_neg=g[:, 1],
            mat_uns_pos=m[:, 0], mat_uns_neg=m[:, 1],
            nas_uns_pos=n[:, 0], nas_uns_neg=n[:, 1],
            mat_spl=m[:, 2] + m[:, 3], nas_spl=n[:, 2] + n[:, 3],
        )

    def override_masses(self, region_arrays) -> dict:
        """The TRUE per-region CalibrationResult mass arrays, built DIRECTLY from the per-origin substrates
        (the exact schema calibrate assembles). gDNA = the gdna partition's unspliced mass; RNA = the
        (mrna+nrna) unspliced mass + ALL spliced (spliced-inclusive, matching cal's
        ``rna = (1−f_g)·M_unspliced + M_spliced``); ``mass_rna_spliced`` = the full spliced (identical to
        cal — spliced is never deconvolved). Conservation ``gdna+rna = total node mass`` holds because the
        partitions sum to the full payload (the validated identity). Feed via
        ``dataclasses.replace(cal, **override_masses(ra))`` for the perfect-calibration lever."""
        from rigel.calibration.substrate import CalibrationSubstrate
        g = CalibrationSubstrate.from_payload(self.parts["gdna"], region_arrays)
        m = CalibrationSubstrate.from_payload(self.parts["mrna"], region_arrays)
        n = CalibrationSubstrate.from_payload(self.parts["nrna"], region_arrays)
        f = CalibrationSubstrate.from_payload(self.full, region_arrays)

        def U(sub, side):
            return np.asarray(getattr(sub, side).mass_unspliced, np.float64)

        def Sp(sub, side):
            return np.asarray(getattr(sub, side).mass_spliced, np.float64)

        rna_u = lambda side: U(m, side) + U(n, side)  # noqa: E731 (unspliced RNA = mrna+nrna)
        spl = lambda side: Sp(f, side)                # noqa: E731 (all spliced is RNA; == cal's)
        return dict(
            mass_gdna_contained=U(g, "contained"),
            mass_rna_contained=rna_u("contained") + spl("contained"),
            mass_gdna_left=U(g, "left"),
            mass_rna_left=rna_u("left") + spl("left"),
            mass_gdna_right=U(g, "right"),
            mass_rna_right=rna_u("right") + spl("right"),
            mass_rna_spliced=spl("contained") + spl("left") + spl("right"),
        )


def _main():
    ap = argparse.ArgumentParser()
    ap.add_argument("condition", nargs="?", default="gdna_gdna300_ss_0.99_nrna_none_capture_on")
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb")
    args = ap.parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    wd = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_oracle_split"
    index = TranscriptIndex.load(f"{args.suite}/rigel_index")
    cfg = PipelineConfig()
    bam = f"{args.suite}/{args.condition}/sim_oracle.bam"
    print(f"=== ORACLE {args.condition} ===")
    orc = OracleTruth.from_bam(bam, index, cfg, wd, args.condition)
    print("VALIDATION PASSED: per-origin partitions sum to the full production payload.")

    G, R = orc.region_unspliced()
    print(f"\nTRUE unspliced contained mass: gDNA={G.sum():,.0f}  RNA={R.sum():,.0f}  "
          f"(RNA = exon-body mRNA + nascent)")

    # calibration accuracy on the CORRECT basis: compare per-region f_g
    from rigel.calibration import calibrate
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.fl import build_fl_models, gdna_fl_mass
    from rigel.splice import SpliceType
    from dataclasses import replace as dc
    sc = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    stats, sm, flm, buffer, payload = scan_and_buffer(bam, index, sc)
    sm.finalize()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(payload), max_size=flm.max_size)
    cal = calibrate(payload=payload, region_arrays=ra, strand_model=sm,
                    gdna_fl_pmf=fl.gdna_pmf, rna_fl_pmf=fl.rna_pmf, config=cfg.calibration)
    cal_g = np.asarray(cal.mass_gdna_contained, np.float64)
    cal_r = np.asarray(cal.mass_rna_contained, np.float64)
    # cal contained total vs payload unspliced contained (should match — spliced contained ~0)
    print(f"\ncal contained total (g+r)={ (cal_g+cal_r).sum():,.0f}  vs TRUE unspliced contained={ (G+R).sum():,.0f}")
    true_fg, tot = orc.region_true_fg()
    cal_fg = np.where((cal_g + cal_r) > 0, cal_g / np.maximum(cal_g + cal_r, 1e-12), np.nan)
    ok = np.isfinite(true_fg) & np.isfinite(cal_fg)
    w = tot[ok]
    err_g_mass = (cal_g - G)  # per-region gDNA contained mass error
    print("\n=== CALIBRATION ACCURACY on the correct (accumulator) basis ===")
    print(f"  contained gDNA mass: cal={cal_g.sum():,.0f}  true={G.sum():,.0f}  "
          f"net err={ (cal_g-G).sum():+,.0f}  Σ|err|={np.abs(cal_g-G).sum():,.0f}")
    mwae = float(np.sum(w * np.abs(cal_fg[ok] - true_fg[ok])) / max(w.sum(), 1))
    print(f"  mass-weighted |Δf_g| = {mwae:.4f}   (0 = perfect per-region gDNA fraction)")
    dirn = np.maximum(G - cal_g, 0.0)  # gDNA under-called (leaks to RNA)
    print(f"  directional gDNA under-call (RNA over-attribution) = {dirn.sum():,.0f}  "
          f"over-call = {np.maximum(cal_g-G,0).sum():,.0f}")


if __name__ == "__main__":
    _main()
