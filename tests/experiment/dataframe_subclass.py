#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
import matplotlib.pyplot as plt
from PyMieSim.units import ureg

from PyMieSim.experiment.detector_set import CoherentModeSet
from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup


def get_experiment_dataframe():
    """
    Generate experiment data in the same way as in your example.
    Returns a PyMieSimDataFrame instance.
    """
    source = GaussianSet(
        wavelength=np.linspace(600, 1000, 50) * ureg.nanometer,
        polarization=PolarizationSet(0 * ureg.degree),
        optical_power=[1e-3] * ureg.watt,
        numerical_aperture=[0.2],
    )
    scatterer = SphereSet(
        diameter=[100] * ureg.nanometer,
        material=[1.4 ],
        medium=[1.1],
    )
    detector = CoherentModeSet(
        mode_number=["LP01"],
        rotation=[0] * ureg.degree,
        numerical_aperture=[0.1, 0.2],
        gamma_offset=[0] * ureg.degree,
        phi_offset=[0] * ureg.degree,
        sampling=[100],
    )
    experiment = Setup(scatterer_set=scatterer, source_set=source, detector_set=detector)
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
    assert isinstance(ax, plt.Figure)

    plt.close()


def test_plot_valid_without_std():
    """
    Test that plotting without std (and with a valid y parameter) returns a valid Axes.
    """
    df = get_experiment_dataframe()
    # We assume that the dataframe has a data column named 'coupling'
    # (or the first column should be the one to plot).
    figure = df.plot(x="source:wavelength", show=False)
    assert isinstance(figure, plt.Figure)
    # Check that at least one line was drawn
    assert len(figure.axes[0].get_lines()) > 0, "No line was plotted for coupling data."
    plt.close()


def test_plot_options_are_applied():
    """Custom plotting options should be applied without the helper decorator."""
    df = get_experiment_dataframe()

    figure = df.plot(
        x="source:wavelength",
        show=False,
        title="Coupling overview",
        figure_size=(5, 3),
        xscale="linear",
        xlim=(600, 1000),
        ylim=(-1, 1),
    )

    axis = figure.axes[0]

    assert figure.get_size_inches()[0] == pytest.approx(5)
    assert figure.get_size_inches()[1] == pytest.approx(3)
    assert axis.get_title() == "Coupling overview"
    assert axis.get_xscale() == "linear"
    assert axis.get_xlim() == pytest.approx((600, 1000))
    assert axis.get_ylim() == pytest.approx((-1, 1))

    plt.close()


def test_plot_supports_polar_projection_alias():
    """The convenience alias ``project`` should enable polar plots."""
    df = get_experiment_dataframe()

    figure = df.plot(x="source:wavelength", show=False, project="polar")

    assert figure.axes[0].name == "polar"

    plt.close()


def test_plot_invalid_x():
    """
    Test that plotting with an invalid x index level raises a ValueError.
    """
    df = get_experiment_dataframe()
    with pytest.raises(ValueError):
        df.plot(x="invalid_x", std="detector:NA", show=False)
    plt.close()


if __name__ == "__main__":
    pytest.main(["-W", "error", __file__])
