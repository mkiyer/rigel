"""Hybrid-capture probe weighting for the simulator.

Re-exports the public API: configuration (:mod:`capture.config`) + runtime sampler
(:mod:`capture.sampler`). Probe design lives in :mod:`capture.design`.
"""

from .config import CaptureConfig, CaptureScenario
from .sampler import CaptureSampler, WeightedInterval

__all__ = ["CaptureConfig", "CaptureScenario", "CaptureSampler", "WeightedInterval"]
