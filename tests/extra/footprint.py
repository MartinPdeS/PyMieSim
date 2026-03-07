#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import matplotlib.pyplot as plt
from unittest.mock import patch
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.polarization import PolarizationState
from PyMieSim.single.detector import Photodiode
from PyMieSim.single.representations import Footprint

@patch("matplotlib.pyplot.show")
def test_sphere_plottings(mock_show_plt):
    """
    Test the plotting functionality for the footprint data of a spherical scatterer.
    The 'matplotlib.pyplot.show' function is patched to avoid displaying the plot during testing.

    Parameters:
        mock_show_plt (Mock): Mock object for 'matplotlib.pyplot.show' to suppress plot display.
    """
    # Create a Gaussian light source
    source = Gaussian(
        wavelength=750 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        optical_power=1 * ureg.watt,
        numerical_aperture=0.3 * ureg.AU,
    )

    # Create a spherical scatterer
    scatterer = Sphere(
        diameter=100 * ureg.nanometer,
        source=source,
        material=1.4 * ureg.RIU,
        medium=1.1 * ureg.RIU,
    )

    # Create a photodiode detector
    detector = Photodiode(
        numerical_aperture=0.1 * ureg.AU,
        phi_offset=0 * ureg.degree,
        gamma_offset=0 * ureg.degree,
        polarization_filter=0 * ureg.degree,
        medium=1.0 * ureg.RIU
    )

    # Retrieve the footprint data for the scatterer
    data = Footprint(scatterer=scatterer, detector=detector)

    # Plot the data (mocked to avoid showing the plot during tests)
    data.plot()

    plt.close()

    # Ensure the data object is not None
    assert data is not None, "Footprint data should not be None."


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
