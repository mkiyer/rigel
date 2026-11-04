"""Length models from a scanned payload alone, for tests that need the second pass's inputs and no index.

Production builds its models with the index's opportunities and partition (``pipeline.library_fl_models``). These
read the payload's own pools through the models' kernel (``fl.fl_models_from_histograms``) with no selection: RNA's
law from the spliced pool's raw counts, gDNA's from the two contained pools. A fixture whose expectations are set by
its own length table (the scorer's ballast) needs exactly that, with no annotation-derived divisor reshaping it.
"""

from __future__ import annotations

import numpy as np

from rigel.calibration.fl import fl_models_from_histograms, rna_fl_mass
from rigel.scan_payload import POOL_DNA_INTERGENIC, POOL_DNA_INTRONIC


def payload_fl_models(payload):
    """The length models over ``payload``'s own pools, with no opportunity divisor."""
    pools = np.asarray(payload.pool_lengths, dtype=np.float64)
    return fl_models_from_histograms(
        global_counts=payload.deposited_lengths,
        rna_counts=rna_fl_mass(payload),
        gdna_counts=pools[[POOL_DNA_INTERGENIC, POOL_DNA_INTRONIC]].sum(axis=0),
        max_size=int(payload.max_length),
        pool_counts=payload.pool_lengths,
    )
