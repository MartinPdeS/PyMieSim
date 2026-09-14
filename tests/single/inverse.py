"""Tests for the dependency-free inverse-fitting API."""

import numpy as np
import pytest

from PyMieSim import FitResult, Observation, Parameter, fit_parameters, ureg


def test_fit_recovers_scalar_parameter_without_scipy():
    observation = Observation(values=np.array([3.0, 5.0, 7.0]), uncertainty=0.1)

    def model(parameters):
        return parameters["offset"] + np.array([1.0, 3.0, 5.0])

    result = fit_parameters(
        model, observation,
        [Parameter("offset", initial=0.0, bounds=(-10.0, 10.0))],
        max_iterations=100,
    )
    assert isinstance(result, FitResult)
    assert result.parameters["offset"] == pytest.approx(2.0, abs=2e-5)
    assert result.objective == pytest.approx(0.0, abs=1e-8)
    assert result.success


def test_fit_passes_units_to_forward_model():
    observation = Observation(values=np.array([1.0, 2.0]))
    seen = []

    def model(parameters):
        seen.append(parameters["diameter"])
        return np.full(2, parameters["diameter"].to("nanometer").magnitude / 100)

    result = fit_parameters(
        model, observation,
        [Parameter("diameter", 100 * ureg.nanometer, (50 * ureg.nanometer, 250 * ureg.nanometer))],
    )
    assert result.parameters["diameter"].units == ureg.nanometer
    assert seen


def test_observation_and_model_shapes_are_checked():
    with pytest.raises(ValueError):
        Observation(values=[])
    with pytest.raises(ValueError):
        Observation(values=[1, 2], uncertainty=[1, 2, 3])
    with pytest.raises(ValueError, match="shape"):
        fit_parameters(lambda _: [1.0], Observation([1.0, 2.0]), [Parameter("x", 0, (-1, 1))])


def test_fit_validates_parameters_and_bounds():
    with pytest.raises(ValueError):
        Parameter("x", 2, (0, 1))
    with pytest.raises(ValueError):
        fit_parameters(lambda _: [1.0], Observation([1.0]), [], max_iterations=1)
    with pytest.raises(ValueError):
        fit_parameters(lambda _: [1.0], Observation([1.0]), [Parameter("x", 0, (-1, 1))], initial_step=0)


def test_fit_progress_is_optional(capsys):
    model = lambda parameters: [parameters["x"]]
    observation = Observation(values=[1.0])
    parameters = [Parameter(name="x", initial=0.0, bounds=(-2.0, 2.0))]

    fit_parameters(model=model, observation=observation, parameters=parameters)
    assert capsys.readouterr().out == ""

    fit_parameters(
        model=model,
        observation=observation,
        parameters=parameters,
        show_progress=True,
    )
    output = capsys.readouterr().out
    assert "iteration" in output
    assert "objective" in output
    assert "converged" in output
