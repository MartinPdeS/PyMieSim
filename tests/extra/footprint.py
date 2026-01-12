#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import matplotlib.pyplot as plt
from unittest.mock import patch
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode


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
        wavelength=750 * ureg.nanometer,  # Wavelength in meters (e.g., 750 nm)
        polarization=0 * ureg.degree,  # Polarization angle
        optical_power=1 * ureg.watt,  # Optical power in ureg.watts
        NA=0.3 * ureg.AU,  # Numerical aperture
    )

    # Create a spherical scatterer
    scatterer = Sphere(
        diameter=100 * ureg.nanometer,  # Diameter in meters (e.g., 100 nm)
        source=source,  # Associated light source
        property=1.4 * ureg.RIU,  # Refractive index of the sphere
        medium_property=1.1 * ureg.RIU,  # Refractive index of the surrounding medium
    )

    # Create a photodiode detector
    detector = Photodiode(
        NA=0.1 * ureg.AU,  # Numerical aperture
        phi_offset=0 * ureg.degree,  # Azimuthal angle offset
        gamma_offset=0 * ureg.degree,  # Polar angle offset
        polarization_filter=0 * ureg.degree,  # Polarization filter angle
    )

    # Retrieve the footprint data for the scatterer
    data = scatterer.get_footprint(detector=detector)

    # Plot the data (mocked to avoid showing the plot during tests)
    data.plot()

    plt.close()

    # Ensure the data object is not None
    assert data is not None, "Footprint data should not be None."


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
