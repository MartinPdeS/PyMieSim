#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import matplotlib.pyplot as plt
from typing import Optional, Sequence, Tuple

from PyMieSim.units import ureg


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
    def __init__(self, setup, sampling: int = 200):
        """
        Initialize the S1/S2 far-field representation.

        Parameters
        ----------
        setup
            Scattering setup used to compute the amplitude functions.
        sampling
            Number of angular samples spanning $[-\\pi, \\pi]$.
        """
        self.setup = setup
        self.sampling = sampling
        self.phi = numpy.linspace(-numpy.pi, numpy.pi, self.sampling) * ureg.radian

        self.S1, self.S2 = self.setup.get_s1s2(
            angles=self.phi
        )

    @staticmethod
    def _finalize_figure(
        figure: plt.Figure,
        *,
        tight_layout: bool = True,
        save_as: Optional[str] = None,
        show: bool = True,
        ylim: Optional[Tuple[float, float]] = None,
    ) -> plt.Figure:
        """
        Apply common figure finalization steps for S1/S2 plots.

        Parameters
        ----------
        figure
            Figure produced by the plotting routine.
        tight_layout
            If True, call ``figure.tight_layout()`` before returning.
        save_as
            Optional output path used to save the figure.
        show
            If True, display the figure with ``plt.show()``.
        ylim
            Optional radial limits applied to every polar subplot.

        Returns
        -------
        plt.Figure
            Finalized matplotlib figure.
        """
        if ylim is not None:
            for ax in figure.axes:
                ax.set_ylim(*ylim)

        if tight_layout:
            figure.tight_layout()

        if save_as is not None:
            figure.savefig(save_as, dpi=300)

        if show:
            plt.show()

        return figure

    def plot(
        self,
        figure_size: Tuple[float, float] = (12, 5),
        titles: Sequence[str] = (r"S$_1$ parameter", r"S$_2$ parameter"),
        tight_layout: bool = True,
        save_as: Optional[str] = None,
        show: bool = True,
        ylim: Optional[Tuple[float, float]] = None,
        style=None,
    ) -> plt.Figure:
        """
        Plot the S1 and S2 scattering amplitudes on polar subplots.

        Parameters
        ----------
        figure_size
            Figure size passed to ``plt.subplots``.
        titles
            Two subplot titles applied to the S1 and S2 panels.
        tight_layout
            If True, call ``figure.tight_layout()`` before returning.
        save_as
            Optional output path used to save the figure.
        show
            If True, display the figure with ``plt.show()``.
        ylim
            Optional radial limits applied to both polar subplots.
        style
            Matplotlib style applied while creating the figure. If omitted, the
            runtime-loaded ``MPSPlots.styles.scientific`` style is used.

        Returns
        -------
        plt.Figure
            Figure containing the two polar amplitude plots.
        """
        if len(titles) != 2:
            raise ValueError("titles must contain exactly two entries for the S1 and S2 subplots")

        if style is None:
            import MPSPlots
            style = MPSPlots.styles.scientific

        with plt.style.context(style):
            figure, axes = plt.subplots(1, 2, figsize=figure_size, subplot_kw={"polar": True})

            axes[0].set(title=titles[0])
            axes[0].fill_between(
                self.phi,
                y1=0,
                y2=numpy.abs(self.S1),
                color="C0",
                alpha=0.7,
                edgecolor="black",
            )

            axes[1].set(title=titles[1])
            axes[1].fill_between(
                self.phi,
                y1=0,
                y2=numpy.abs(self.S2),
                color="C1",
                alpha=0.7,
            )

            for ax in axes:
                ax.set_yticklabels([])
                ax.set_xlabel("")
                ax.set_ylabel("")

            return self._finalize_figure(
                figure,
                tight_layout=tight_layout,
                save_as=save_as,
                show=show,
                ylim=ylim,
            )
