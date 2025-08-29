#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
import matplotlib.pyplot as plt
from TypedUnit import ureg

from PyMieSim.experiment.detector import CoherentMode
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup


def get_experiment_dataframe():
    """
    Generate experiment data in the same way as in your example.
    Returns a PyMieSimDataFrame instance.
    """
    source = Gaussian(
        wavelength=np.linspace(600, 1000, 50) * ureg.nanometer,
        polarization=0 * ureg.degree,
        optical_power=1e-3 * ureg.watt,
        NA=0.2 * ureg.AU,
    )
    scatterer = Sphere(
        diameter=100 * ureg.nanometer,
        source=source,
        property=1.4 * ureg.RIU,
        medium_property=1.1 * ureg.RIU,
    )
    detector = CoherentMode(
        mode_number="LP01",
        rotation=0 * ureg.degree,
        NA=[0.1, 0.2] * ureg.AU,
        polarization_filter=None,
        gamma_offset=0 * ureg.degree,
        phi_offset=0 * ureg.degree,
        sampling=100 * ureg.AU,
    )
    experiment = Setup(scatterer=scatterer, source=source, detector=detector)
    # Get the measurement dataframe (assumed to be a PyMieSimDataFrame subclass)
    dataframe = experiment.get("coupling", drop_unique_level=True, scale_unit=True)
    return dataframe


def test_plot_valid_with_std():
    """
    Test that plotting with a valid x index and std parameter returns a valid Axes.
    """
    df = get_experiment_dataframe()
    # We assume that the coupling dataframe has a MultiIndex with
    # levels including 'source:wavelength' and 'detector:NA'.
    ax = df.plot(x="source:wavelength", std="detector:NA", show=False)
    assert isinstance(ax, plt.Axes)

    plt.close()


def test_plot_valid_without_std():
    """
    Test that plotting without std (and with a valid y parameter) returns a valid Axes.
    """
    df = get_experiment_dataframe()
    # We assume that the dataframe has a data column named 'coupling'
    # (or the first column should be the one to plot).
    ax = df.plot(x="source:wavelength", show=False)
    assert isinstance(ax, plt.Axes)
    # Check that at least one line was drawn
    assert len(ax.get_lines()) > 0, "No line was plotted for coupling data."
    plt.close()


def test_plot_invalid_x():
    """
    Test that plotting with an invalid x index level raises a ValueError.
    """
    df = get_experiment_dataframe()
    with pytest.raises(
        ValueError, match="x parameter 'invalid_x' is not in the DataFrame index"
    ):
        df.plot(x="invalid_x", std="detector:NA", show=False)
    plt.close()


if __name__ == "__main__":
    pytest.main(["-W", "error", __file__])
