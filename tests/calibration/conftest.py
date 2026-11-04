"""Shared fixtures for the calibration tests.

`sweep_inputs` is the fixture every transfer-policy test needs: a real chain driven through the two-phase backbone by
`_transfer_harness.capture_sweep_inputs`, captured once per module so a module of row-level assertions pays for it once,
and defined here so its scope and builder have one home.
"""

from __future__ import annotations

import pytest

from _transfer_harness import capture_sweep_inputs


@pytest.fixture(scope="module")
def sweep_inputs(tmp_path_factory):
    """The captured backbone inputs for one module's worth of transfer-policy assertions."""
    return capture_sweep_inputs(tmp_path_factory)
