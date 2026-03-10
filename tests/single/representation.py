#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.units import ureg
from unittest.mock import patch
import matplotlib.pyplot as plt


from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.polarization import PolarizationState
from PyMieSim.single.detector import Photodiode
from PyMieSim.single import Setup

representation_list = ["farfields", "stokes", "spf", "s1s2", "footprint"]



# Parametrized test for plotting functions
@patch("pyvista.Plotter.show")
@patch("matplotlib.pyplot.show")
@pytest.mark.parametrize("representation", representation_list)
def test_plottings(mock_show_plt, mock_show_pyvista, representation):

    source = Gaussian(
        wavelength=750 * ureg.nanometer,
        polarization=PolarizationState(angle=0 * ureg.degree),
        optical_power=1 * ureg.watt,
        numerical_aperture=0.3 * ureg.AU,
    )

    scatterer = Sphere(
        diameter=100 * ureg.nanometer,
        medium=1.0 * ureg.RIU,
        material=1.4 * ureg.RIU,
    )

    detector = Photodiode(
        sampling=100,
        numerical_aperture=0.2 * ureg.AU,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
        medium=1.0 * ureg.RIU
    )

    setup = Setup(
        scatterer=scatterer,
        source=source,
        detector=detector,
    )

    data = setup.get_representation(representation)

    assert data is not None
    data.plot()
    plt.close()


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
