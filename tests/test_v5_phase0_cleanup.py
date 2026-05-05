"""cleanup tests for calibration.

Validates the destructive cleanup contract:

1. ``rigel.calibration`` package imports cleanly (no v1 SRD code left).
2. ``calibrate`` raises :class:`NotImplementedError` until Phase 9.
3. ``bootstrap_fl_calibration`` returns a :class:`CalibrationStub`
 with a populated ``gdna_fl_model``.
4. End-to-end pipeline run yields ``RunResult.calibration.version == "stub"``.
5. No file under ``src/rigel/`` matches the regex ``\\bSRD v[123]\\b``
 (regression guard for any reintroduction of v1 SRD references).
6. ``summary.json``-style summary dict has ``calibration.active is False``.
7. ``compute_default_locus_priors`` returns equal gDNA/RNA pseudocounts.

See ``docs/calibration/calibration.md`` §3 for the
authoritative contract.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

import rigel.calibration as cal_pkg
from rigel.calibration import (
 CalibrationStub,
 bootstrap_fl_calibration,
 calibrate,
)
from rigel.frag_length_model import FragmentLengthModels
from rigel.locus import compute_default_locus_priors
from rigel.scored_fragments import MultiLocus


_SRD_VERSIONED_RE = re.compile(r"\bSRD v[123]\b")
_RIGEL_SRC = Path(__file__).resolve().parent.parent / "src" / "rigel"


def _make_finalized_fl_models(n_obs: int = 200) -> FragmentLengthModels:
 """Build a finalized FL container with synthetic observations."""
 flm = FragmentLengthModels(max_size=600)
 rng = np.random.default_rng(42)
 lengths = rng.normal(loc=200.0, scale=40.0, size=n_obs).astype(np.int64)
 lengths = np.clip(lengths, 50, 500)
 flm.global_model.observe_batch(lengths)
 flm.build_scoring_models()
 flm.finalize()
 return flm


def test_package_imports_cleanly():
 """The new package surface exports exactly three public names."""
 assert hasattr(cal_pkg, "CalibrationStub")
 assert hasattr(cal_pkg, "bootstrap_fl_calibration")
 assert hasattr(cal_pkg, "calibrate")
 assert set(cal_pkg.__all__) == {
 "CalibrationStub",
 "bootstrap_fl_calibration",
 "calibrate",
 # global gDNA-density estimation surface (§2.6, §5.1.3).
 "GlobalGdnaDensity",
 "GlobalDensityTable",
 "compute_global_densities",
 "KappaEstimate",
 "estimate_kappa",
 "KAPPA_DEFAULT",
 "KAPPA_MIN",
 "KAPPA_MAX",
 # locoregional EB shrinkage + per-Locus gDNA prior (§5.2.2).
 "shrink_to_loco",
 "C_BASE_DEFAULT",
 "LocusGdnaEstimate",
 "MultiLocusPrior",
 "PriorTable",
 "estimate_locus_gdna",
 "assemble_multilocus_prior",
 "assemble_priors",
 "build_prior_weight_rna",
 "FLAG_INTERGENIC_ZERO_LEFF",
 "FLAG_INTRON_ZERO_LEFF",
 "FLAG_EXON_INTRON_NO_ELIGIBLE",
 "FLAG_PI_CLIPPED",
 "FLAG_NO_REGIONS",
 # pool fragment-length models (§5.3.1).
 "PoolFLModels",
 "compute_pool_fl_models",
 "DEFAULT_PRIOR_ESS",
 "DEFAULT_QUALITY_THRESHOLD_GOOD",
 "DEFAULT_QUALITY_THRESHOLD_WEAK",
 # CalibrationResult schema + builder (§5.3.2).
 "CalibrationResult",
 "build_calibration_result",
 }


def test_calibrate_v5_raises_until_phase9():
 """``calibrate`` must fail loudly until lands."""
 with pytest.raises(NotImplementedError, match="Phase 9"):
 calibrate()


def test_bootstrap_returns_stub_with_gdna_fl_model():
 """The bootstrap returns a CalibrationStub with a populated FL model."""
 flm = _make_finalized_fl_models()
 stub = bootstrap_fl_calibration(flm)

 assert isinstance(stub, CalibrationStub)
 assert stub.version == "stub"
 assert stub.active is False
 assert stub.gdna_fl_model is not None
 # Must be a real, non-degenerate distribution (mean roughly tracks
 # the global histogram mean we fed in).
 assert 100.0 < stub.gdna_fl_model.mean < 300.0


def test_pipeline_run_uses_stub(tmp_path):
 """End-to-end run leaves ``RunResult.calibration.version == 'stub'``."""
 pytest.importorskip("rigel.sim")
 from rigel.config import BamScanConfig, EMConfig, PipelineConfig
 from rigel.index import TranscriptIndex
 from rigel.pipeline import run_pipeline
 from rigel.sim import Scenario, SimConfig

 seed = 42
 sc = Scenario("phase0_stub", genome_length=5000, seed=seed,
 work_dir=tmp_path / "phase0")
 sc.add_gene("g1", "+", [
 {"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 50},
 ])
 sim_cfg = SimConfig(
 frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
 read_length=100, strand_specificity=1.0, seed=seed,
 )
 res = sc.build_oracle(n_fragments=200, sim_config=sim_cfg)
 index = TranscriptIndex.load(res.index_dir)

 pcfg = PipelineConfig(
 em=EMConfig(seed=seed),
 scan=BamScanConfig(sj_strand_tag="auto"),
 )
 pr = run_pipeline(res.bam_path, index, config=pcfg)

 assert pr.calibration is not None
 assert pr.calibration.version == "stub"
 assert pr.calibration.active is False


def test_no_srd_versioned_strings_in_src():
 """Regression guard: no ``SRD v1/v2/v3`` references survive in src/."""
 offenders: list[tuple[Path, int, str]] = []
 for path in _RIGEL_SRC.rglob("*"):
 if not path.is_file():
 continue
 if path.suffix not in {".py", ".cpp", ".h", ".hpp", ".cc"}:
 continue
 try:
 text = path.read_text(encoding="utf-8")
 except UnicodeDecodeError:
 continue
 for lineno, line in enumerate(text.splitlines(), start=1):
 if _SRD_VERSIONED_RE.search(line):
 offenders.append((path, lineno, line.strip()))
 assert not offenders, (
 "Found SRD v[123] references in src/:\n"
 + "\n".join(f" {p}:{n} {ln}" for p, n, ln in offenders)
 )


def test_summary_dict_marks_calibration_inactive():
 """``CalibrationStub.to_summary_dict()`` reports ``active=False``."""
 flm = _make_finalized_fl_models()
 stub = bootstrap_fl_calibration(flm)
 d = stub.to_summary_dict()
 assert d["version"] == "stub"
 assert d["active"] is False
 assert "gdna_fl_mean" in d


def test_compute_default_locus_priors_is_symmetric():
 """Default per-locus prior is split equally between gDNA and RNA."""
 loci = [
 MultiLocus(locus_id=i, transcript_indices=np.array([i], dtype=np.int32),
 unit_indices=np.array([], dtype=np.int32), gdna_span=1000,
 merged_intervals=[("chr1", 0, 1000)])
 for i in range(5)
 ]
 alpha_g, alpha_r = compute_default_locus_priors(loci)
 assert alpha_g.shape == (5,)
 assert alpha_r.shape == (5,)
 np.testing.assert_allclose(alpha_g, alpha_r)
 # Total mass per locus equals the documented c_base default.
 np.testing.assert_allclose(alpha_g + alpha_r, 10.0)
