import pytest
from unittest.mock import patch
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import InfiniteCylinder
from PyMieSim.single.source import Gaussian
from PyMieSim.single.polarization import PolarizationState
from PyMieSim.single.detector import Photodiode
from PyMieSim.single import SystemPlotter
from PyMieSim.single.representations import SPF
from math import sqrt
import matplotlib.pyplot as plt


@patch("pyvista.Plotter.show")
def test_plot_system(mock_show):
    """
    Test the plot_system function to ensure it can plot a source, scatterer, and detector
    without calling the actual show method from PyVista.

    """
    source = Gaussian(
        wavelength=1550 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        optical_power=1 * ureg.watt,
        numerical_aperture=0.3 * ureg.AU,
    )

    scatterer = InfiniteCylinder(
        diameter=780 * ureg.nanometer,
        source=source,
        medium=1.0 * ureg.RIU,
        material=sqrt(1.5) * ureg.RIU,
    )

    detector = Photodiode(
        numerical_aperture=0.1 * ureg.AU,
        gamma_offset=90 * ureg.degree,
        phi_offset=0 * ureg.degree,
        polarization_filter=0 * ureg.degree,
        medium=1.0 * ureg.RIU
    )
    spf = SPF(scatterer=scatterer)

    plotter = SystemPlotter(show_axis_label=False)

    plotter.plot(source, detector, data=spf)

    plt.close()

    mock_show.assert_called_once()


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
