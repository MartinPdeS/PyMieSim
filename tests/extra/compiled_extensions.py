"""Direct tests for the distributions and inverse compiled extensions."""

import numpy as np
import pytest

import PyMieSim.distributions as distributions
import PyMieSim.inverse as inverse


@pytest.mark.parametrize("order", [1, 2, 4, 8, 16, 32])
def test_compiled_legendre_matches_numpy(order):
    nodes, weights = distributions.legendre(order=order)
    expected_nodes, expected_weights = np.polynomial.legendre.leggauss(order)

    np.testing.assert_allclose(nodes, expected_nodes, rtol=1e-14, atol=1e-15)
    np.testing.assert_allclose(weights, expected_weights, rtol=1e-13, atol=1e-15)


@pytest.mark.parametrize("order", [1, 2, 4, 8, 16, 32])
def test_compiled_hermite_matches_numpy(order):
    nodes, weights = distributions.hermite(order=order)
    expected_nodes, expected_weights = np.polynomial.hermite.hermgauss(order)

    np.testing.assert_allclose(nodes, expected_nodes, rtol=1e-13, atol=1e-14)
    np.testing.assert_allclose(weights, expected_weights, rtol=1e-12, atol=1e-14)


def test_compiled_inverse_api_and_progress(capsys):
    observation = inverse.Observation(values=[3.0])
    parameter = inverse.Parameter(name="offset", initial=0.0, bounds=(-5.0, 5.0))

    result = inverse.fit_parameters(
        model=lambda values: [values["offset"]],
        observation=observation,
        parameters=[parameter],
        show_progress=True,
    )

    assert result.success
    assert result.parameters["offset"] == pytest.approx(3.0, abs=1e-5)
    assert "objective" in capsys.readouterr().out
