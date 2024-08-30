import pytest
from unittest.mock import patch
from PyMieSim.single.scatterer import Cylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from PyMieSim.single import plot_system
from math import sqrt

@pytest.fixture
def source():
    """
    Fixture to create a Gaussian source with predefined parameters.

    Returns:
        Gaussian: A Gaussian light source object.
    """
    return Gaussian(
        wavelength=1550e-9,  # 1550 nm wavelength
        polarization=0,      # Linear polarization angle in radians
        optical_power=1,     # Optical power in arbitrary units
        NA=0.3               # Numerical Aperture
    )

@pytest.fixture
def scatterer(source):
    """
    Fixture to create a Cylinder scatterer with predefined parameters.

    Args:
        source (Gaussian): The Gaussian light source used for scattering.

    Returns:
        Cylinder: A cylindrical scatterer object.
    """
    return Cylinder(
        diameter=7.8e-6,     # 7.8 micrometers diameter
        source=source,
        medium_index=1.0,    # Refractive index of the surrounding medium
        index=sqrt(1.5)      # Refractive index of the scatterer
    )

@pytest.fixture
def detector():
    """
    Fixture to create a Photodiode detector with predefined parameters.

    Returns:
        Photodiode: A photodiode detector object.
    """
    return Photodiode(
        NA=0.1,                 # Numerical Aperture
        gamma_offset=90,        # Gamma offset in degrees
        phi_offset=0,           # Phi offset in degrees
        polarization_filter=0   # Polarization filter angle in degrees
    )

@patch('pyvista.Plotter.show')
def test_plot_system(mock_show, source, scatterer, detector):
    """
    Test the plot_system function to ensure it can plot a source, scatterer, and detector
    without calling the actual show method from PyVista.

    Args:
        mock_show (MagicMock): Mock object for pyvista.Plotter.show.
        source (Gaussian): Fixture providing a Gaussian light source.
        scatterer (Cylinder): Fixture providing a Cylinder scatterer.
        detector (Photodiode): Fixture providing a Photodiode detector.
    """
    spf = scatterer.get_spf()
    plot_system(source, detector, spf)

    # Ensure the show method was called exactly once
    mock_show.assert_called_once()

if __name__ == "__main__":
    pytest.main([__file__])