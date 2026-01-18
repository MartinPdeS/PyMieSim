#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.units import ureg
from unittest.mock import patch
import matplotlib.pyplot as plt


from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single.representations import FarField, Stokes, SPF, S1S2

representations = [FarField, Stokes, SPF, S1S2]



# Parametrized test for plotting functions
@pytest.mark.parametrize("representation", representations)
@patch("pyvista.Plotter.show")
@patch("matplotlib.pyplot.show")
def test_plottings(mock_show_plt, mock_show_pyvista, representation):

    source = Gaussian(
        wavelength=750 * ureg.nanometer,
        polarization=0 * ureg.degree,
        optical_power=1 * ureg.watt,
        NA=0.3 * ureg.AU,
    )

    scatterer = Sphere(
        diameter=100 * ureg.nanometer,
        source=source,
        medium_refractive_index=1.0 * ureg.RIU,
        refractive_index=1.4 * ureg.RIU,
    )

    data = representation(scatterer=scatterer)

    assert data is not None
    data.plot()
    plt.close()


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
