import pytest
from unittest.mock import patch
from PyMieSim.single.scatterer import Cylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.detector import Photodiode
from PyMieSim.single import plot_system
from math import sqrt
import matplotlib.pyplot as plt
from PyMieSim.units import nanometer, degree, watt, AU, RIU


@pytest.fixture
def source():
    """
    Fixture to create a Gaussian source with predefined parameters.

    Returns:
        Gaussian: A Gaussian light source object.
    """
    return Gaussian(
        wavelength=1550 * nanometer,  # 1550 nm wavelength
        polarization=0 * degree,      # Linear polarization angle in radians
        optical_power=1 * watt,     # Optical power in arbitrary units
        NA=0.3 * AU               # Numerical Aperture
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
        diameter=780 * nanometer,     # 7.8 micrometers diameter
        source=source,
        medium_property=1.0 * RIU,    # Refractive index of the surrounding medium
        property=sqrt(1.5) * RIU      # Refractive index of the scatterer
    )


@pytest.fixture
def detector():
    """
    Fixture to create a Photodiode detector with predefined parameters.

    Returns:
        Photodiode: A photodiode detector object.
    """
    return Photodiode(
        NA=0.1 * AU,                     # Numerical Aperture
        gamma_offset=90 * degree,        # Gamma offset in degrees
        phi_offset=0 * degree,           # Phi offset in degrees
        polarization_filter=0 * degree   # Polarization filter angle in degrees
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

    plt.close()

    # Ensure the show method was called exactly once
    mock_show.assert_called_once()


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
