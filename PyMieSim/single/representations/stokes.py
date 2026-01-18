#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from pydantic.dataclasses import dataclass
from typing import List
import pyvista
from MPSPlots.colormaps import blue_black_red

from PyMieSim.units import ureg
from PyMieSim.single.representations.base import BaseRepresentation
from PyMieSim.utils import config_dict, spherical_to_cartesian


@dataclass(config=config_dict, kw_only=True)
class Stokes():
    r"""
    Compute the four Stokes parameters: I, Q, U, and V, which describe the polarization state of scattered light.

    The Stokes parameters provide a complete description of the polarization state of electromagnetic radiation, and are defined as follows:

    .. math::
        I &= \big| E_x \big|^2 + \big| E_y \big|^2 \quad &\text{(total intensity)} \\
        Q &= \big| E_x \big|^2 - \big| E_y \big|^2 \quad &\text{(linear polarization along x and y)} \\
        U &= 2 \mathcal{Re} \big\{ E_x E_y^* \big\} \quad &\text{(linear polarization at +45° and -45°)} \\
        V &= 2 \mathcal{Im} \big\{ E_x E_y^* \big\} \quad &\text{(circular polarization)}

    Where:

    - :math:`I`: Total intensity of the scattered light.
    - :math:`Q`: The degree of linear polarization along the x and y axes.
    - :math:`U`: The degree of linear polarization at +45° and -45° to the axes.
    - :math:`V`: The degree of circular polarization.

    These parameters provide a powerful way to represent the full polarization state of the scattered light.

    Parameters
    ----------
    sampling : int
        The number of angular points used to sample the Stokes parameters. A higher sampling value increases the resolution of the computed polarization pattern.
    distance : Length, optional
        The distance from the scatterer at which the Stokes parameters are evaluated. The default is 1 meter, but this can be adjusted based on the experimental setup.

    Returns
    -------
    Stokes
        An object containing the computed Stokes parameters (I, Q, U, V), which describe the polarization state of the scattered light.

    Notes
    -----
    - The Stokes parameters are essential for understanding and characterizing the polarization of scattered light. They are used in a wide variety of applications, including atmospheric optics, remote sensing, and optical communication.
    - The `sampling` parameter controls the angular resolution of the calculated Stokes parameters. Higher values provide finer detail, which can be important in accurately modeling polarization effects in scattering experiments.

    Example
    -------
    You can use this method to compute the Stokes parameters for a spherical scatterer, helping you analyze how the scatterer affects the polarization state of light.

    Example usage:

    >>> stokes_params = scatterer.get_stokes(sampling=500)
    >>> print(stokes_params.I, stokes_params.Q, stokes_params.U, stokes_params.V)

    """
    scatterer: object
    sampling: int

    def __post_init__(self):
        fields = self.scatterer.get_full_farfields(
            sampling=self.sampling, distance=1 * ureg.meter
        )
        self.E_phi, self.E_theta, self.theta, self.phi = fields

        self.compute_components()

    def get_colormap_limits(self, scalar: numpy.ndarray, symmetric: bool = False):
        if symmetric:
            max_abs = numpy.abs(scalar).max()
            return [-max_abs, max_abs]
        else:
            return None

    def compute_components(self) -> None:
        r"""
        Computes the Stokes parameters (I, Q, U, V) based on the electric field components (E_phi and E_theta).

        The method calculates the normalized intensity (I), linear polarizations (Q and U), and circular polarization (V) of the light
        scattered by the particle, using the electric field components in spherical coordinates.

        The Stokes parameters are calculated using the following formulas:

        .. math:
            - I = |E_phi|^2 + |E_theta|^2
            - Q = |E_phi|^2 - |E_theta|^2
            - U = 2 * Re{E_phi * E_theta*}
            - V = -2 * Im{E_phi * E_theta*}

        The results are stored as attributes of the instance: I, Q, U, and V.

        """
        intensity = numpy.abs(self.E_phi) ** 2 + numpy.abs(self.E_theta) ** 2

        self.I = intensity / numpy.max(intensity)  # noqa: E741
        self.Q = (numpy.abs(self.E_phi) ** 2 - numpy.abs(self.E_theta) ** 2) / intensity
        self.U = (+2 * numpy.real(self.E_phi * self.E_theta.conjugate())) / intensity
        self.V = (-2 * numpy.imag(self.E_phi * self.E_theta.conjugate())) / intensity

    def plot(
        self,
        unit_size: List[float] = (400, 400),
        background_color: str = "white",
        show_edges: bool = False,
        colormap: str = blue_black_red,
        opacity: float = 1.0,
        symmetric_colormap: bool = False,
        show_axis_label: bool = False,
    ) -> None:
        """
        Visualizes the Stokes parameters (I, Q, U, V) on a 3D plot.

        Parameters
        ----------
        unit_size : List[float]
            The size of each subplot in pixels (width, height). Default is (400, 400).
        background_color : str
            The background color of the plot. Default is 'white'.
        show_edges : bool
            If True, displays the edges of the mesh. Default is False.
        colormap : str
            The colormap to use for scalar mapping. Default is 'blue_black_red'.
        opacity : float
            The opacity of the mesh. Default is 1.0.
        symmetric_colormap : bool
            If True, the colormap will be symmetric around zero. Default is False.
        show_axis_label : bool
            If True, shows the axis labels. Default is False.
        """
        phi_mesh, theta_mesh = numpy.meshgrid(self.phi, self.theta)
        x, y, z = spherical_to_cartesian(
            r=numpy.full_like(phi_mesh, 0.5), phi=phi_mesh, theta=theta_mesh
        )

        window_size = (unit_size[1] * 4, unit_size[0])  # Four subplots horizontally

        scene = pyvista.Plotter(
            theme=pyvista.themes.DocumentTheme(), window_size=window_size, shape=(1, 4)
        )
        scene.set_background(background_color)

        for idx, (name, field) in enumerate(
            zip(["I", "Q", "U", "V"], [self.I, self.Q, self.U, self.V])
        ):
            field = field.flatten(order="F")
            mesh = pyvista.StructuredGrid(x, y, z)
            scene.subplot(0, idx)

            colormap_limits = self.get_colormap_limits(
                scalar=field, symmetric=symmetric_colormap
            )

            mapping = scene.add_mesh(
                mesh,
                cmap=colormap,
                scalars=field,
                opacity=opacity,
                style="surface",
                show_edges=show_edges,
                clim=colormap_limits,
                show_scalar_bar=False,
            )

            scene.add_axes_at_origin(labels_off=not show_axis_label)
            scene.add_scalar_bar(mapper=mapping.mapper, title=f"{name} field")

        scene.show()
