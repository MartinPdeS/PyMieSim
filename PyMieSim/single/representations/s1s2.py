#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from pydantic.dataclasses import dataclass
from MPSPlots import helper

from PyMieSim.utils import config_dict
from PyMieSim.units import ureg

@dataclass(config=config_dict, kw_only=True)
class S1S2():
    r"""
    Compute the S1 and S2 scattering amplitude functions for a spherical scatterer.

    The S1 and S2 parameters represent the scattering amplitudes for perpendicular and parallel polarizations of light, respectively. These parameters are fundamental in Mie theory, which describes the scattering of electromagnetic waves by spherical particles.

    The formulas for \( S_1 \) and \( S_2 \) are:

    .. math::
        S_1 = \sum\limits_{n=1}^{n_{\text{max}}} \frac{2n+1}{n(n+1)} \left( a_n \pi_n + b_n \tau_n \right) \\
        S_2 = \sum\limits_{n=1}^{n_{\text{max}}} \frac{2n+1}{n(n+1)} \left( a_n \tau_n + b_n \pi_n \right)

    Where:

    - :math:`a_n` and :math:`b_n`: Mie coefficients, which depend on the size, shape, and refractive index of the scatterer.
    - :math:`\pi_n` and :math:`\tau_n` \)`: Angular functions related to the angular components of the incident and scattered fields.
    - :math:`n_{\text{max}}`: Maximum number of terms in the series, determined by the size parameter of the scatterer.

    These scattering amplitude functions are essential for calculating properties such as scattering phase functions, efficiencies, and angular distribution of scattered light.

    Parameters
    ----------
    sampling : int
        The number of angular points used to sample the S1 and S2 functions. Higher sampling improves the resolution of the scattering pattern but increases computation time.
    distance : Length, optional
        The distance from the scatterer at which the S1 and S2 parameters are evaluated. This is typically set to 1 meter by default, but can be adjusted for specific setups.

    Returns
    -------
    representations.S1S2
        An object containing the computed S1 and S2 parameters, representing the scattering amplitudes for the two polarization components.

    Notes
    -----
    - The S1 and S2 parameters are central to Mie scattering theory and are used to derive many important scattering properties, such as intensity distributions and polarization effects.
    - The `sampling` parameter controls how finely the angular distribution is resolved. A higher value of `sampling` provides more detailed scattering information, which can be critical for accurately modeling the far-field pattern.

    Example
    -------
    You can use this method to compute the scattering properties of spherical particles, particularly in experiments where the polarization and scattering pattern of the light are important.

    Example usage:

    >>> s1s2 = scatterer.get_s1s2(sampling=500)
    >>> print(s1s2.S1, s1s2.S2)

    """
    scatterer: object
    sampling: int = 200

    def __post_init__(self):
        self.phi = numpy.linspace(-180, 180, self.sampling) * ureg.degree

        self.S1, self.S2 = self.scatterer.get_s1s2(
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
