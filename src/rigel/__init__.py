"""rigel — Bayesian RNA-seq transcript quantification that separates RNA from genomic-DNA contamination.

A single-pass C++ BAM scan tallies the fragments, calibration deconvolves the library into gDNA and RNA, and a
per-locus EM assigns the RNA to transcripts.
"""

import os as _os

# Keep numpy's OpenMP from spawning idle thread pools that compete with rigel's own C++ threads.
_os.environ.setdefault("OMP_NUM_THREADS", "1")

from importlib.metadata import version as _version, PackageNotFoundError

try:
    __version__ = _version("rigel-rnaseq")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"
