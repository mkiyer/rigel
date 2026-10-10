"""An intergenic background estimate must not become a second intron count prior."""

import importlib


from rigel.config import CalibrationConfig

calibrate_module = importlib.import_module("rigel.calibration.calibrate")


def _calibrate(inputs):
    return calibrate_module.calibrate(
        payload=inputs["payload"],
        **inputs["calibrate_kw"],
        config=CalibrationConfig(),
    )


def test_initial_and_refitted_sweeps_do_not_receive_an_intron_prior(sweep_inputs, monkeypatch):
    original = calibrate_module.solve_chain
    seen = []

    def observe(*args, **kwargs):
        seen.append((kwargs["gdna_prior"] is None, kwargs.get("intron_prior")))
        return original(*args, **kwargs)

    monkeypatch.setattr(calibrate_module, "solve_chain", observe)
    _calibrate(sweep_inputs)
    assert seen and seen[0][0], "The initial solve must be exercised"
    assert any(not initial for initial, _ in seen), "A fitted-landscape solve must be exercised"
    assert all(prior is None for _, prior in seen), "An intron factor still enters a count solve"


def test_count_inference_does_not_evaluate_the_unused_background_fit(sweep_inputs, monkeypatch):
    def unwanted_fit(*args, **kwargs):
        raise AssertionError("The removed intron prior still causes a background fit")

    # Restoring the old imported fitter makes this gate fail. No unused fitter is
    # retained in production merely to provide the injection seam.
    monkeypatch.setattr(calibrate_module, "fit_intron_background", unwanted_fit, raising=False)
    _calibrate(sweep_inputs)
