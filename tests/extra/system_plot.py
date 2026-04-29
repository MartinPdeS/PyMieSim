import pytest
from unittest.mock import patch
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import InfiniteCylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.polarization import PolarizationState
from PyMieSim.single.detector import Photodiode
import matplotlib.pyplot as plt

@patch("matplotlib.pyplot.show")
def test_plot_system(mock_show):
    """
    Test the plot_system function to ensure it can plot a source, scatterer, and detector
    without calling the actual show method from PyVista.

    """
    source = Gaussian(
        wavelength=1550 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        optical_power=1 * ureg.watt,
        numerical_aperture=0.3,
    )

    scatterer = InfiniteCylinder(
        diameter=780 * ureg.nanometer,
        medium=1.0,
        material=1.5,
    )

    scatterer.material.initialize(source.wavelength)
    scatterer.medium.initialize(source.wavelength)

    detector = Photodiode(
        numerical_aperture=0.1,
        gamma_offset=90 * ureg.degree,
        phi_offset=0 * ureg.degree,
        polarization_filter=0 * ureg.degree,
        medium=1.0
    )

    detector.medium.initialize(source.wavelength)
    detector.initialize_mesh(scatterer)

    figure = plt.figure()
    ax = figure.add_subplot(111, projection="3d")

    scatterer._add_to_ax(ax)
    source._add_to_ax(ax)
    detector._add_to_ax(ax)

    plt.show()
    mock_show.assert_called_once()


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
