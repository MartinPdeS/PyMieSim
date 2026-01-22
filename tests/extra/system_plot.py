import pytest
from unittest.mock import patch
from PyMieSim.units import ureg

from PyMieSim.single.scatterer import InfiniteCylinder
from PyMieSim.single.source import Gaussian
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
        wavelength=1550 * ureg.nanometer,  # 1550 nm wavelength
        polarization=0 * ureg.degree,  # Linear polarization angle in radians
        optical_power=1 * ureg.watt,  # Optical power in arbitrary units
        NA=0.3 * ureg.AU,  # Numerical Aperture
    )

    scatterer = InfiniteCylinder(
        diameter=780 * ureg.nanometer,  # 7.8 micrometers diameter
        source=source,
        medium_refractive_index=1.0 * ureg.RIU,  # Refractive index of the surrounding medium
        refractive_index=sqrt(1.5) * ureg.RIU,  # Refractive index of the scatterer
    )

    detector = Photodiode(
        NA=0.1 * ureg.AU,  # Numerical Aperture
        gamma_offset=90 * ureg.degree,  # Gamma offset in ureg.degrees
        phi_offset=0 * ureg.degree,  # Phi offset in ureg.degrees
        polarization_filter=0* ureg.degree,  # Polarization filter angle in ureg.degrees
        medium_refractive_index=1.0 * ureg.RIU
    )
    spf = SPF(scatterer=scatterer)

    plotter = SystemPlotter(show_axis_label=False)

    plotter.plot(source, detector, data=spf)

    plt.close()

    mock_show.assert_called_once()


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
