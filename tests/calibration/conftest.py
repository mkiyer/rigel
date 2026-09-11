"""Shared fixtures for the calibration tests.

`sweep_inputs` is the one every transfer-policy test needs: a real chain driven through the
two-phase backbone by `_transfer_harness.capture_sweep_inputs`, captured once per module so a file
of row-level assertions pays for it once. It lives here rather than in each file because eight
copies of the same two lines is eight places for the scope or the builder to drift.
"""

from __future__ import annotations

import pytest

from _transfer_harness import capture_sweep_inputs


@pytest.fixture(scope="module")
def sweep_inputs(tmp_path_factory):
    """The captured backbone inputs for one module's worth of transfer-policy assertions."""
    return capture_sweep_inputs(tmp_path_factory)
