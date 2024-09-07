#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from unittest.mock import patch
import matplotlib.pyplot as plt

@patch('matplotlib.pyplot.show')
def test_sphere_plottings(mock_show_plt):
    """
    Test the plotting functionality for the footprint data of a spherical scatterer.
    The 'matplotlib.pyplot.show' function is patched to avoid displaying the plot during testing.

    Parameters:
        mock_show_plt (Mock): Mock object for 'matplotlib.pyplot.show' to suppress plot display.
    """
    # Create a Gaussian light source
    source = Gaussian(
        wavelength=750e-9,   # Wavelength in meters (e.g., 750 nm)
        polarization=0,      # Polarization angle
        optical_power=1,     # Optical power in watts
        NA=0.3               # Numerical aperture
    )

    # Create a spherical scatterer
    scatterer = Sphere(
        diameter=100e-9,     # Diameter in meters (e.g., 100 nm)
        source=source,       # Associated light source
        index=1.4,           # Refractive index of the sphere
        medium_index=1.1     # Refractive index of the surrounding medium
    )

    # Create a photodiode detector
    detector = Photodiode(
        NA=0.1,              # Numerical aperture
        phi_offset=0,        # Azimuthal angle offset
        gamma_offset=0,      # Polar angle offset
        polarization_filter=0 # Polarization filter angle
    )

    # Retrieve the footprint data for the scatterer
    data = scatterer.get_footprint(detector=detector)

    # Plot the data (mocked to avoid showing the plot during tests)
    data.plot()

    plt.close()

    # Ensure the data object is not None
    assert data is not None, "Footprint data should not be None."


if __name__ == "__main__":
    pytest.main([__file__])
