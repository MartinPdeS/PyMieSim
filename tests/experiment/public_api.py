"""Tests for the public parameter-sweep API."""

from PyMieSim import (
    CoreShellSet,
    Experiment,
    GaussianSet,
    PhotodiodeSet,
    PolarizationSet,
    SellmeierMaterial,
    SellmeierMedium,
    TabulatedMaterial,
)


def test_experiment_components_are_public():
    assert Experiment.__name__ == "Setup"
    assert all(
        component is not None
        for component in (
            CoreShellSet,
            GaussianSet,
            PhotodiodeSet,
            PolarizationSet,
            SellmeierMaterial,
            TabulatedMaterial,
            SellmeierMedium,
        )
    )
