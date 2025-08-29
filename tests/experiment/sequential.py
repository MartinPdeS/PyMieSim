#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from unittest.mock import patch
from TypedUnit import ureg

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyOptik import Material

TOTAL_SIZE = 5


@patch("matplotlib.pyplot.show")
def test_valid_experiment(mock_show):
    """
    Test a valid experiment configuration with proper broadcasting.

    The test constructs a valid source, scatterer, and detector using scalar and array
    parameters, then runs the experiment. It is expected that no error is raised.
    """
    source = Gaussian.build_sequential(
        wavelength=np.linspace(600, 1000, TOTAL_SIZE) * ureg.nanometer,
        polarization=0 * ureg.degree,
        optical_power=1e-3 * ureg.watt,
        NA=0.2 * ureg.AU,
        total_size=TOTAL_SIZE,
    )

    scatterer = Sphere.build_sequential(
        source=source,
        diameter=np.linspace(400, 1400, TOTAL_SIZE) * ureg.nanometer,
        property=1.4 * ureg.RIU,
        medium_property=1.0 * ureg.RIU,
        total_size=TOTAL_SIZE,
    )

    detector = CoherentMode.build_sequential(
        mode_number="LP01",
        rotation=0 * ureg.degree,
        NA=0.2 * ureg.AU,
        polarization_filter=np.nan * ureg.degree,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
        sampling=100 * ureg.AU,
        total_size=TOTAL_SIZE,
    )

    experiment = Setup(scatterer=scatterer, source=source, detector=detector)

    # This call should complete without raising an error.
    experiment.get_sequential(Sphere.available_measure_list[0])


def test_invalid_medium_property():
    """
    Test that using an invalid medium property (Material.water) raises a ValueError.

    When medium_property is defined using Material.water, the configuration should be rejected.
    """
    source = Gaussian.build_sequential(
        wavelength=np.linspace(600, 1000, TOTAL_SIZE) * ureg.nanometer,
        polarization=0 * ureg.degree,
        optical_power=1e-3 * ureg.watt,
        NA=0.2 * ureg.AU,
        total_size=TOTAL_SIZE,
    )

    with pytest.raises(AssertionError):
        Sphere.build_sequential(
            source=source,
            diameter=np.linspace(400, 1400, TOTAL_SIZE) * ureg.nanometer,
            property=1.4 * ureg.RIU,
            medium_property=Material.water,  # This should trigger an error
            total_size=TOTAL_SIZE,
        )


def test_invalid_property():
    """
    Test that using an invalid property (Material.water) for the scatterer raises a ValueError.

    When property is defined using Material.water, the configuration should be rejected.
    """
    source = Gaussian.build_sequential(
        wavelength=np.linspace(600, 1000, TOTAL_SIZE) * ureg.nanometer,
        polarization=0 * ureg.degree,
        optical_power=1e-3 * ureg.watt,
        NA=0.2 * ureg.AU,
        total_size=TOTAL_SIZE,
    )

    with pytest.raises(AssertionError):
        Sphere.build_sequential(
            source=source,
            diameter=np.linspace(400, 1400, TOTAL_SIZE) * ureg.nanometer,
            property=Material.water,  # This should trigger an error
            medium_property=1.0 * ureg.RIU,
            total_size=TOTAL_SIZE,
        )


def test_parameter_broadcasting():
    """
    Test that scalar parameters are correctly broadcast to arrays of the specified total size.

    This test uses scalar values for diameter, property, and medium_property, and verifies that
    the resulting arrays all have a size equal to TOTAL_SIZE.
    """
    source = Gaussian.build_sequential(
        wavelength=800 * ureg.nanometer,  # Scalar value
        polarization=45 * ureg.degree,
        optical_power=2e-3 * ureg.watt,
        NA=0.3 * ureg.AU,
        total_size=TOTAL_SIZE,
    )

    scatterer = Sphere.build_sequential(
        source=source,
        diameter=500 * ureg.nanometer,  # Scalar value
        property=1.5 * ureg.RIU,  # Scalar value
        medium_property=1.2 * ureg.RIU,  # Scalar value
        total_size=TOTAL_SIZE,
    )

    # Check that each broadcasted attribute is an array with the proper length.
    assert scatterer.diameter.size == TOTAL_SIZE
    assert scatterer.property.size == TOTAL_SIZE
    assert scatterer.medium_property.size == TOTAL_SIZE


def test_broadcasted_values():
    """
    Test that broadcasted arrays contain the expected repeated scalar values.

    This test verifies that when scalar parameters are provided, the resulting broadcasted arrays
    consist entirely of the original scalar value.
    """
    value = 700 * ureg.nanometer
    source = Gaussian.build_sequential(
        wavelength=value,
        polarization=0 * ureg.degree,
        optical_power=1e-3 * ureg.watt,
        NA=0.2 * ureg.AU,
        total_size=TOTAL_SIZE,
    )
    scatterer = Sphere.build_sequential(
        source=source,
        diameter=value,
        property=1.6 * ureg.RIU,
        medium_property=1.0 * ureg.RIU,
        total_size=TOTAL_SIZE,
    )

    # Create an expected array filled with the broadcasted value.
    expected = np.full(TOTAL_SIZE, value.magnitude, dtype=object) * value.units
    np.testing.assert_array_equal(scatterer.diameter.magnitude, expected.magnitude)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
