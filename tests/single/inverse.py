"""Tests for the dependency-free inverse-fitting API."""

import numpy as np
import pytest
from pint.errors import DimensionalityError

from PyMieSim import FitParameter, MeasuredData, ParameterFitResult, fit_parameters, ureg


@pytest.mark.parametrize(
    ("name", "legacy_name"),
    [("FitParameter", "Parameter"), ("MeasuredData", "Observation"), ("ParameterFitResult", "FitResult")],
)
def test_descriptive_class_names_and_compatibility_aliases(name, legacy_name):
    import PyMieSim
    import PyMieSim.inverse as inverse

    cls = getattr(inverse, name)
    assert cls.__name__ == name
    assert getattr(PyMieSim, name) is cls
    assert getattr(inverse, legacy_name) is cls
    assert getattr(PyMieSim, legacy_name) is cls
    assert name in PyMieSim.__all__


def test_fit_recovers_scalar_parameter_without_scipy():
    observation = MeasuredData(values=np.array([3.0, 5.0, 7.0]), uncertainty=0.1)

    def model(parameters):
        return parameters["offset"] + np.array([1.0, 3.0, 5.0])

    result = fit_parameters(
        model, observation,
        [FitParameter("offset", initial=0.0, bounds=(-10.0, 10.0))],
        max_iterations=100,
    )
    assert isinstance(result, ParameterFitResult)
    assert repr(result).startswith("ParameterFitResult(")
    assert result.parameters["offset"] == pytest.approx(2.0, abs=2e-5)
    assert result.objective == pytest.approx(0.0, abs=1e-8)
    assert result.success


def test_fit_passes_units_to_forward_model():
    observation = MeasuredData(values=np.array([1.0, 2.0]))
    seen = []

    def model(parameters):
        seen.append(parameters["diameter"])
        return np.full(2, parameters["diameter"].to("nanometer").magnitude / 100)

    result = fit_parameters(
        model, observation,
        [FitParameter("diameter", 100 * ureg.nanometer, (50 * ureg.nanometer, 250 * ureg.nanometer))],
    )
    assert result.parameters["diameter"].units == ureg.nanometer
    assert seen


def test_fit_converts_mixed_units_and_handles_zero_initial_quantity():
    seen = []

    def model(parameters):
        diameter = parameters["diameter"]
        seen.append(diameter)
        return (2 * diameter).to("micrometer")

    parameter = FitParameter(
        name="diameter",
        initial=0 * ureg.nanometer,
        bounds=(-0.1 * ureg.micrometer, 0.0003 * ureg.millimeter),
    )
    result = fit_parameters(
        model=model,
        observation=MeasuredData(values=240 * ureg.nanometer, uncertainty=0.01 * ureg.micrometer),
        parameters=[parameter],
    )
    assert result.success
    assert result.parameters["diameter"].units == ureg.nanometer
    assert result.parameters["diameter"].magnitude == pytest.approx(120, abs=0.001)
    assert result.initial_parameters["diameter"].magnitude == 0
    assert all(value.units == ureg.nanometer and np.isfinite(value.magnitude) for value in seen)
    assert result.prediction == pytest.approx(240, abs=0.002)
    assert result.objective == pytest.approx(0, abs=1e-7)


def test_residuals_convert_predictions_and_uncertainty_to_observation_units():
    result = fit_parameters(
        model=lambda _: np.array([1.0, 3.0]) * ureg.micrometer,
        observation=MeasuredData(
            values=np.array([500.0, 2000.0]) * ureg.nanometer,
            uncertainty=np.array([0.1, 0.2]) * ureg.micrometer,
        ),
        parameters=[FitParameter(name="x", initial=0, bounds=(-1, 1))],
    )
    np.testing.assert_allclose(result.prediction, [1000, 3000])
    np.testing.assert_allclose(result.observed, [500, 2000])
    np.testing.assert_allclose(result.residuals, [5, 5])
    assert result.objective == pytest.approx(50)


def test_parameter_validates_bounds_after_unit_conversion():
    FitParameter(
        name="diameter",
        initial=150 * ureg.nanometer,
        bounds=(0.05 * ureg.micrometer, 0.0003 * ureg.millimeter),
    )
    with pytest.raises(ValueError, match="initial value must lie within bounds"):
        FitParameter(
            name="diameter",
            initial=150 * ureg.nanometer,
            bounds=(0.2 * ureg.micrometer, 300 * ureg.nanometer),
        )
    with pytest.raises(ValueError, match="lower bound < upper bound"):
        FitParameter(
            name="diameter",
            initial=150 * ureg.nanometer,
            bounds=(0.2 * ureg.micrometer, 100 * ureg.nanometer),
        )


@pytest.mark.parametrize("initial", [1.0, 1 * ureg.nanometer])
def test_parameter_rejects_incompatible_bound_dimensions(initial):
    with pytest.raises(DimensionalityError):
        FitParameter(name="x", initial=initial, bounds=(0 * ureg.second, 2 * ureg.second))


@pytest.mark.parametrize("values", [1.0, 1 * ureg.nanometer])
def test_observation_rejects_incompatible_uncertainty_dimensions(values):
    with pytest.raises(DimensionalityError):
        MeasuredData(values=values, uncertainty=1 * ureg.second)


@pytest.mark.parametrize("values", [1.0, 1 * ureg.nanometer])
def test_fit_rejects_incompatible_prediction_dimensions(values):
    with pytest.raises(DimensionalityError):
        fit_parameters(
            model=lambda _: 1 * ureg.second,
            observation=MeasuredData(values=values),
            parameters=[FitParameter(name="x", initial=0, bounds=(-1, 1))],
        )


def test_bare_magnitudes_use_reference_units():
    result = fit_parameters(
        model=lambda parameters: parameters["length"].magnitude,
        observation=MeasuredData(values=2 * ureg.nanometer, uncertainty=0.1),
        parameters=[FitParameter(name="length", initial=0 * ureg.nanometer, bounds=(-5, 5))],
    )
    assert result.parameters["length"].magnitude == pytest.approx(2, abs=2e-5)


def test_dimensionless_quantities_are_converted_for_bare_references():
    result = fit_parameters(
        model=lambda parameters: parameters["fraction"] * 100 * ureg.percent,
        observation=MeasuredData(values=0.5, uncertainty=10 * ureg.percent),
        parameters=[FitParameter(name="fraction", initial=0, bounds=(0 * ureg.percent, 100 * ureg.percent))],
    )
    assert result.parameters["fraction"] == pytest.approx(0.5, abs=2e-5)
    assert result.objective == pytest.approx(0, abs=1e-8)


def test_observation_and_model_shapes_are_checked():
    with pytest.raises(ValueError):
        MeasuredData(values=[])
    with pytest.raises(ValueError):
        MeasuredData(values=[1, 2], uncertainty=[1, 2, 3])
    with pytest.raises(ValueError, match="shape"):
        fit_parameters(lambda _: [1.0], MeasuredData([1.0, 2.0]), [FitParameter("x", 0, (-1, 1))])


def test_fit_validates_parameters_and_bounds():
    with pytest.raises(ValueError):
        FitParameter("x", 2, (0, 1))
    with pytest.raises(ValueError):
        fit_parameters(lambda _: [1.0], MeasuredData([1.0]), [], max_iterations=1)
    with pytest.raises(ValueError):
        fit_parameters(lambda _: [1.0], MeasuredData([1.0]), [FitParameter("x", 0, (-1, 1))], initial_step=0)


def test_fit_progress_is_optional(capsys):
    model = lambda parameters: [parameters["x"]]
    observation = MeasuredData(values=[1.0])
    parameters = [FitParameter(name="x", initial=0.0, bounds=(-2.0, 2.0))]

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
