"""
hulkrna.core — Backward-compatible re-export shim.

All foundational types have moved to ``hulkrna.types`` and
``hulkrna.categories``.  This module re-exports everything so that
existing ``from .core import ...`` statements continue to work.
"""

# Re-export foundational types
from .types import (  # noqa: F401
    Strand,
    Interval,
    GenomicInterval,
    IntervalType,
    RefInterval,
    MergeCriteria,
    MergeResult,
    EMPTY_MERGE,
)

# Re-export count enums
from .categories import (  # noqa: F401
    CountCategory,
    CountStrand,
    CountType,
    NUM_COUNT_TYPES,
)

# merge_sets_with_criteria lives in hulkrna.resolution now.
# Re-export it here for backward compatibility, but import lazily
# to avoid circular imports (resolution → fragment → core).
def merge_sets_with_criteria(*args, **kwargs):
    from .resolution import merge_sets_with_criteria as _fn
    return _fn(*args, **kwargs)
