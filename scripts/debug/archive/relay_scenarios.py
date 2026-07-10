"""Phase 0 — relay-exposing benchmark scenarios + oracle baseline measurement.

Builds the two user-specified scenarios that EXPOSE the message-relay gap (which the existing
sim suite hides), runs the CURRENT calibration on each, and compares per-region gDNA mass to a
by-origin oracle (split the BAM into gdna-only / rna-only by the simulator read name, scan each).

  Scenario A — tiny middle exon (RNA relay across a gap):
    T1(+) exons (1000,2000),(10000,11000);  T2(-) exons (5000,6000),(7000,7000+L),(8000,9000)
    Sweep the middle-exon length L down through the fragment length → it zeroes out and KILLS the
    -strand relay between T2's outer exons. Run stranded AND unstranded.

  Scenario B — gDNA-signal killers (gDNA relay across zero-count regions):
    Many small single-exon transcripts (length L_k) tiled along the genome on (+); one (-) transcript
    with a giant intron overlapping them; uniform gDNA background. Sweep L_k down through the gDNA
    fragment length → killers go to zero counts and break the gDNA-density relay → underestimate.

Determinism: OMP_NUM_THREADS=1. Usage:
    OMP_NUM_THREADS=1 python scripts/debug/relay_scenarios.py [A|B|both]
"""
from __future__ import annotations

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from rigel.config import BamScanConfig, CalibrationConfig
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario
from rigel.pipeline import scan_and_buffer
from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.signature import coarse_strand_from_signature, coarse_type_from_signature
from rigel.splice import SpliceType

_EPS = 1e-9
SEED = 42
SC = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}
TC = {0: "intergenic", 1: "intron", 2: "exon"}


def _sim_cfg(ss: float) -> ReadSimConfig:
    return ReadSimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                         read_length=100, strand_specificity=ss, seed=SEED)


def _gdna_cfg(abundance: float) -> GDNAConfig | None:
    if abundance <= 0:
        return None
    return GDNAConfig(abundance=abundance, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000)


# ── scenario builders ────────────────────────────────────────────────────────

def build_scenario_A(L: int, work_dir: Path) -> Scenario:
    """Tiny middle exon (7000, 7000+L) on T2(-); long-intron T1(+) overlaps all of T2."""
    sc = Scenario(f"relayA_L{L}", genome_length=60000, seed=SEED, work_dir=work_dir)
    sc.add_gene("g1", "+", [{"t_id": "T1", "exons": [(1000, 2000), (10000, 11000)], "abundance": 100}])
    sc.add_gene("g2", "-", [{"t_id": "T2",
                             "exons": [(5000, 6000), (7000, 7000 + L), (8500, 9500)], "abundance": 100}])
    # training transcripts (multi-exon, both strands) so FL + strand + var~mean models can fit
    for k, base in enumerate(range(20000, 50000, 6000)):
        strand = "+" if k % 2 == 0 else "-"
        sc.add_gene(f"tr{k}", strand, [{"t_id": f"TR{k}",
                    "exons": [(base, base + 800), (base + 1500, base + 2300), (base + 3000, base + 3800)],
                    "abundance": 80}])
    return sc


def build_scenario_B(L_k: int, work_dir: Path) -> Scenario:
    """Many small single-exon (+) killers tiled 5000..45000; one (-) giant-intron transcript over them."""
    sc = Scenario(f"relayB_Lk{L_k}", genome_length=60000, seed=SEED, work_dir=work_dir)
    killers = []
    for j, pos in enumerate(range(5000, 45001, 2500)):
        killers.append({"t_id": f"K{j}", "exons": [(pos, pos + L_k)], "abundance": 60})
    sc.add_gene("killers", "+", killers)
    # the (-) giant-intron transcript: two far-apart exons, intron spans the whole killer field
    sc.add_gene("giant", "-", [{"t_id": "GIANT", "exons": [(2000, 3000), (47000, 48000)], "abundance": 100}])
    # training transcripts far from the killer field
    for k, base in enumerate(range(50000, 58000, 3500)):
        sc.add_gene(f"tr{k}", "+" if k % 2 == 0 else "-",
                    [{"t_id": f"TR{k}", "exons": [(base, base + 600), (base + 1200, base + 1800)],
                      "abundance": 80}])
    return sc


# ── calibrate + oracle ───────────────────────────────────────────────────────

def _scan(bam: Path, index, scan_cfg):
    stats, strand, fl, buffer, payload = scan_and_buffer(str(bam), index, scan_cfg)
    del buffer
    return strand, fl, payload


def _split_by_origin(src: Path, gdna_out: Path, rna_out: Path) -> None:
    subprocess.run(f"samtools view -h '{src}' | awk '/^@/||$1~/^gdna:/' | samtools view -b -o '{gdna_out}' -",
                   shell=True, check=True)
    subprocess.run(f"samtools view -h '{src}' | awk '/^@/||$1!~/^gdna:/' | samtools view -b -o '{rna_out}' -",
                   shell=True, check=True)


def calibrate_and_oracle(result, scan_cfg):
    """Run the production calibrate on the full BAM + a by-origin oracle. Returns (index, cal, g_or, r_or)."""
    index = result.index
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    strand_full, fl_full, payload_full = _scan(result.bam_path, index, scan_cfg)

    fl_models = build_fl_models(
        global_counts=fl_full.global_model.counts,
        rna_counts=fl_full.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(payload_full),
        max_size=fl_full.max_size,
    )
    cal = calibrate(payload=payload_full, region_arrays=ra, strand_model=strand_full,
                    gdna_fl_pmf=fl_models.gdna_pmf, rna_fl_pmf=fl_models.rna_pmf,
                    config=CalibrationConfig())

    tmp = Path(tempfile.mkdtemp())
    gbam, rbam = tmp / "g.bam", tmp / "r.bam"
    _split_by_origin(result.bam_path, gbam, rbam)
    g_or = np.asarray(CalibrationSubstrate.from_payload(_scan(gbam, index, scan_cfg)[2], ra)
                      .contained.mass_unspliced, float)
    r_or = np.asarray(CalibrationSubstrate.from_payload(_scan(rbam, index, scan_cfg)[2], ra)
                      .contained.mass_unspliced, float)
    return index, cal, g_or, r_or


def region_table(index, cal, g_or, r_or, lo: int, hi: int):
    """Per-region production-vs-oracle gDNA for regions overlapping [lo, hi]."""
    rdf = index.region_df.reset_index(drop=True)
    sig = rdf["signature"].to_numpy()
    g_pr = np.asarray(cal.mass_gdna_contained, float)
    r_pr = np.asarray(cal.mass_rna_contained, float)
    start = rdf.start.to_numpy(); end = rdf.end.to_numpy()
    sel = np.where((end > lo) & (start < hi))[0]
    rows = []
    for i in sel:
        rows.append(dict(reg=int(i), start=int(start[i]), end=int(end[i]),
                         cls=SC[coarse_strand_from_signature(int(sig[i]))],
                         typ=TC[coarse_type_from_signature(int(sig[i]))],
                         g_or=g_or[i], g_pr=g_pr[i], r_or=r_or[i], r_pr=r_pr[i]))
    return rows


# ── drivers ──────────────────────────────────────────────────────────────────

def run_A():
    print("\n" + "=" * 96)
    print("SCENARIO A — tiny middle exon (RNA relay). gDNA=0 ⇒ oracle gDNA=0; any prod gDNA is PHANTOM.")
    print("=" * 96)
    scan_cfg = BamScanConfig(sj_strand_tag="auto")
    for ss in (0.99, 0.50):
        print(f"\n  --- strand_specificity = {ss} ({'stranded' if ss > 0.9 else 'UNSTRANDED'}) ---")
        print(f"  {'L':>5} {'tot_gPHANTOM':>13} {'gDNA@T2regions':>15} {'#T2 regs':>9}")
        for L in (1000, 300, 100, 50, 10):
            with tempfile.TemporaryDirectory() as wd:
                sc = build_scenario_A(L, Path(wd))
                result = sc.build_oracle(n_fragments=4000, sim_config=_sim_cfg(ss), gdna_config=None)
                index, cal, g_or, r_or = calibrate_and_oracle(result, scan_cfg)
                g_pr = np.asarray(cal.mass_gdna_contained, float)
                rows = region_table(index, cal, g_or, r_or, 4500, 10000)  # the T2 span
                t2_g = sum(r["g_pr"] for r in rows)
                print(f"  {L:>5} {g_pr.sum():>13,.0f} {t2_g:>15,.1f} {len(rows):>9}")


def run_B():
    print("\n" + "=" * 96)
    print("SCENARIO B — gDNA killers (gDNA relay). Uniform gDNA background ⇒ killers SHOULD read background.")
    print("=" * 96)
    scan_cfg = BamScanConfig(sj_strand_tag="auto")
    print(f"  {'L_k':>5} {'tot_g_oracle':>13} {'tot_g_prod':>12} {'recov%':>7} "
          f"{'giant-intron g_or':>18} {'giant-intron g_pr':>18}")
    for L_k in (800, 400, 200, 100, 50):
        with tempfile.TemporaryDirectory() as wd:
            sc = build_scenario_B(L_k, Path(wd))
            result = sc.build_oracle(n_fragments=8000, sim_config=_sim_cfg(0.99),
                                     gdna_config=_gdna_cfg(50))
            index, cal, g_or, r_or = calibrate_and_oracle(result, scan_cfg)
            g_pr = np.asarray(cal.mass_gdna_contained, float)
            # giant-intron region(s): the (-) transcript's intron 3000..47000 (intergenic/intron-class, big)
            rows = region_table(index, cal, g_or, r_or, 3000, 47000)
            big = [r for r in rows if (r["end"] - r["start"]) > 5000]
            gi_or = sum(r["g_or"] for r in big); gi_pr = sum(r["g_pr"] for r in big)
            recov = 100 * g_pr.sum() / max(g_or.sum(), _EPS)
            print(f"  {L_k:>5} {g_or.sum():>13,.0f} {g_pr.sum():>12,.0f} {recov:>6.0f}% "
                  f"{gi_or:>18,.0f} {gi_pr:>18,.0f}")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    if which in ("A", "both"):
        run_A()
    if which in ("B", "both"):
        run_B()
