"""Minimal driver for py-spy profiling of the rigel pipeline."""
import sys
from rigel.index import TranscriptIndex
from rigel.pipeline import run_pipeline
from rigel.config import PipelineConfig, EMConfig, BamScanConfig

bam = sys.argv[1]
index_dir = sys.argv[2]

idx = TranscriptIndex.load(index_dir)
cfg = PipelineConfig(
    em=EMConfig(n_threads=8),
    scan=BamScanConfig(n_scan_threads=8),
)
result = run_pipeline(bam, idx, cfg)
