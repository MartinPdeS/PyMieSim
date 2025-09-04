#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from pydantic.dataclasses import dataclass
from MPSPlots import helper

from PyMieSim.single.representations.base import BaseRepresentation
from PyMieSim.utils import config_dict


@dataclass(config=config_dict, kw_only=True)
class S1S2(BaseRepresentation):
    """
    Represents the S1 and S2 scattering functions, which are components of the scattering matrix.

    Parameters
    ----------
    scatterer : BaseScatterer
        The scatterer object.
    sampling : int
        Number of points for evaluating the S1 and S2 functions.
    distance : Quantity
        Distance at which the fields are evaluated.

    Methods:
        compute_components: Computes the S1 and S2 functions based on the scatterer's properties.
        plot: Visualizes the S1 and S2 functions on a polar plot.
    """

    def __post_init__(self):
        self.phi = numpy.linspace(-180, 180, self.sampling)

        self.compute_components()

    def compute_components(self) -> None:
        """
        Computes the S1 and S2 scattering parameters based on the scatterer's properties and the scattering angle phi.

        S1 and S2 are integral parts of the scattering matrix describing the change in polarization state of light upon scattering.

        The method calculates these parameters for a range of phi angles and stores them as the S1 and S2 attributes of the instance.
        """
        self.S1, self.S2 = self.scatterer._cpp_get_s1s2(
            phi=numpy.deg2rad(self.phi) + numpy.pi / 2
        )

    @helper.pre_plot(nrows=1, ncols=2, subplot_kw={"polar": True})
    def plot(self, axes) -> None:
        """
        Plots the S1 and S2 Stokes parameters on polar plots.

        The method generates two polar plots: one for the absolute values of the S1 parameter and another
        for the S2 parameter, filling the area between the radial axis and the parameter values.

        Returns
        -------
        None
            This method does not return a value. It displays the polar plots.
        """
        # Plot for S1 parameter
        axes[0].set(title=r"S$_1$ parameter")
        axes[0].fill_between(
            numpy.deg2rad(self.phi),
            y1=0,
            y2=numpy.abs(self.S1),
            color="C0",
            alpha=0.7,
            # edgecolor="black",
        )

        # Plot for S2 parameter
        axes[1].set(title=r"S$_2$ parameter")
        axes[1].fill_between(
            numpy.deg2rad(self.phi),
            y1=0,
            y2=numpy.abs(self.S2),
            color="C1",
            alpha=0.7,
            # edgecolor="black",
        )
