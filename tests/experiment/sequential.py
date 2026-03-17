#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
from unittest.mock import patch

from PyMieSim.units import ureg

from PyMieSim.material import SellmeierMaterial, SellmeierMedium
from PyMieSim.experiment.detector_set import CoherentModeSet
from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup


TOTAL_SIZE = 5


@patch("matplotlib.pyplot.show")
def test_valid_experiment(mock_show):
    """
    Test a valid experiment configuration with proper broadcasting.

    The test constructs a valid source, scatterer, and detector using scalar and array
    parameters, then runs the experiment. It is expected that no error is raised.
    """
    polarization_set = PolarizationSet(angles=[0] * TOTAL_SIZE * ureg.degree)
    source = GaussianSet.build_sequential(
        target_size=TOTAL_SIZE,
        wavelength=np.linspace(600, 1000, TOTAL_SIZE) * ureg.nanometer,
        polarization=polarization_set,
        optical_power=1e-3 * ureg.watt,
        numerical_aperture=0.2,
    )

    scatterer = SphereSet.build_sequential(
        diameter=np.linspace(400, 1400, TOTAL_SIZE) * ureg.nanometer,
        material=1.4,
        medium=1.0,
        target_size=TOTAL_SIZE,
    )

    detector = CoherentModeSet.build_sequential(
        mode_number="LP01",
        rotation=0 * ureg.degree,
        numerical_aperture=[0.2],
        polarization_filter=np.nan * ureg.degree,
        medium=1.0,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
        sampling=100,
        target_size=TOTAL_SIZE,
    )

    experiment = Setup(
        scatterer_set=scatterer,
        source_set=source,
        detector_set=detector
    )

    # This call should complete without raising an error.
    experiment.get_sequential(SphereSet.available_measure_list[0])


def test_invalid_medium_refractive_index():
    """
    Test that using an invalid medium refractive_index (Material.water) raises a ValueError.

    When medium_refractive_index is defined using Material.water, the configuration should be rejected.
    """
    with pytest.raises((AttributeError, RuntimeError)):
        SphereSet.build_sequential(
            diameter=np.linspace(400, 1400, TOTAL_SIZE) * ureg.nanometer,
            material=1.4,
            medium=SellmeierMedium("water"),  # This should trigger an error
            target_size=TOTAL_SIZE,
        )


def test_invalid_refractive_index():
    """
    Test that using an invalid refractive_index (Material.water) for the scatterer raises a ValueError.

    When refractive_index is defined using Material.water, the configuration should be rejected.
    """
    with pytest.raises((AttributeError, RuntimeError)):
        SphereSet.build_sequential(
            diameter=np.linspace(400, 1400, TOTAL_SIZE) * ureg.nanometer,
            material=SellmeierMaterial("water"),  # This should trigger an error
            medium=1.0,
            target_size=TOTAL_SIZE,
        )


def test_parameter_broadcasting():
    """
    Test that scalar parameters are correctly broadcast to arrays of the specified total size.

    This test uses scalar values for diameter, refractive_index, and medium_refractive_index, and verifies that
    the resulting arrays all have a size equal to TOTAL_SIZE.
    """
    scatterer = SphereSet.build_sequential(
        diameter=500 * ureg.nanometer,  # Scalar value
        material=1.5,  # Scalar value
        medium=1.2,  # Scalar value
        target_size=TOTAL_SIZE,
    )

    # Check that each broadcasted attribute is an array with the proper length.

    assert scatterer.diameter.size == TOTAL_SIZE
    assert len(scatterer.material) == TOTAL_SIZE
    assert len(scatterer.medium) == TOTAL_SIZE


def test_broadcasted_values():
    """
    Test that broadcasted arrays contain the expected repeated scalar values.

    This test verifies that when scalar parameters are provided, the resulting broadcasted arrays
    consist entirely of the original scalar value.
    """
    value = 700 * ureg.nanometer

    scatterer = SphereSet.build_sequential(
        diameter=value,
        material=1.6,
        medium=1.0,
        target_size=TOTAL_SIZE,
    )

    # Create an expected array filled with the broadcasted value.
    expected = np.full(TOTAL_SIZE, value.magnitude) * value.units
    np.testing.assert_allclose(
        scatterer.diameter.to("nm").magnitude,
        expected.to("nm").magnitude
    )


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
