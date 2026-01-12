#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from PyMieSim.units import ureg
from unittest.mock import patch
import matplotlib.pyplot as plt


from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian
from PyMieSim.single import plot_system

representations = ["get_farfield", "get_stokes", "get_spf", "get_s1s2"]


@pytest.fixture()
def source():
    return Gaussian(
        wavelength=750 * ureg.nanometer,
        polarization=0 * ureg.degree,
        optical_power=1 * ureg.watt,
        NA=0.3 * ureg.AU,
    )


@pytest.fixture()
def scatterer(source):
    return Sphere(
        diameter=100 * ureg.nanometer,
        source=source,
        medium_property=1.0 * ureg.RIU,
        property=1.4 * ureg.RIU,
    )


# Parametrized test for plotting functions
@pytest.mark.parametrize("representation", representations)
@patch("pyvista.Plotter.show")
@patch("matplotlib.pyplot.show")
def test_plottings(mock_show_plt, mock_show_pyvista, representation, scatterer):
    data = getattr(scatterer, representation)()
    assert data is not None
    data.plot()
    plt.close()


@patch("pyvista.Plotter.show")
def test_plot_system(mock_show, scatterer, source):
    plot_system(scatterer, source)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
