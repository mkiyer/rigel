"""Hand-built synthetic substrate inputs for calibration unit tests.

A plain helper module (not ``conftest.py``) to avoid colliding with the
top-level ``conftest`` that several tests import by name.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS
from rigel.scan_payload import AccumulatorPayload


def make_synthetic_payload() -> tuple[AccumulatorPayload, RegionArrays]:
    """A 1-ref, 3-region payload + aligned RegionArrays with known channels.

    Regions (chr1): r0 +exon, r1 −exon, r2 intergenic.
    ``region_contained`` rows ``[ch0,ch1,ch2,ch3]`` = ``[unspl+, unspl−, spl+, spl−]``::

        r0: [10, 2, 3, 0]    r1: [1, 20, 0, 5]    r2: [7, 8, 0, 0]

    4 boundaries b0..b3 (b0, b3 terminal). b1 is shared by r0's right view and
    r1's left view::

        b1: flux [4,1,0,0]  mass_left [1.5,0.5,0,0]  mass_right [2.5,0.5,0,0]
        b2: flux [2,3,1,0]  mass_left [3.0,1.0,0.5,0] mass_right [0.5,1.0,0,0]
        b0, b3: all zero (reference terminals)
    """
    region_contained = np.array([[10, 2, 3, 0], [1, 20, 0, 5], [7, 8, 0, 0]], dtype=np.uint32)
    flux = np.zeros((4, 4), dtype=np.uint32)
    mass_left = np.zeros((4, 4), dtype=np.float32)
    mass_right = np.zeros((4, 4), dtype=np.float32)
    flux[1] = [4, 1, 0, 0]
    mass_left[1] = [1.5, 0.5, 0, 0]
    mass_right[1] = [2.5, 0.5, 0, 0]
    flux[2] = [2, 3, 1, 0]
    mass_left[2] = [3.0, 1.0, 0.5, 0]
    mass_right[2] = [0.5, 1.0, 0, 0]

    payload = AccumulatorPayload(
        boundaries=np.array([0, 100, 200, 300], dtype=np.int64),
        ref_pos_offsets=np.array([0, 4], dtype=np.int64),
        ref_region_offsets=np.array([0, 3], dtype=np.int64),
        ref_boundary_offsets=np.array([0, 4], dtype=np.int64),
        region_contained=region_contained,
        boundary_mass_left=mass_left,
        boundary_mass_right=mass_right,
        boundary_flux=flux,
        n_refs=1,
    )
    region_df = pd.DataFrame(
        {
            "region_id": np.arange(3, dtype=np.int64),
            "ref_name": pd.array(["chr1"] * 3, dtype="string"),
            "start": np.array([0, 100, 200], dtype=np.int64),
            "end": np.array([100, 200, 300], dtype=np.int64),
            "length": np.array([100, 100, 100], dtype=np.int64),
            "signature": np.array([BIT_EXON_POS, BIT_EXON_NEG, 0], dtype=np.uint8),
        }
    )
    region_arrays = RegionArrays.from_region_df(region_df, {"chr1": 0})
    return payload, region_arrays
