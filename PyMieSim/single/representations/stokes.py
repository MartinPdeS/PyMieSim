#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import List
import pyvista
from MPSPlots.colormaps import blue_black_red

from PyMieSim.units import ureg
from PyMieSim.single.full_mesh import FullMesh # Necessary for loading the class, even if not directly used in this file


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

    def __init__(self, setup: object, sampling: int = 200, distance: ureg.Quantity = 1.0 * ureg.meter):
        self.setup = setup
        self.sampling = sampling
        self.distance = distance

        self.I, self.Q, self.U, self.V, self.mesh = self.setup.get_stokes(
            sampling=self.sampling, distance=self.distance
        )

    def plot(
        self,
        unit_size: List[float] = (400, 400),
        background_color: str = "white",
        show_edges: bool = False,
        colormap: str = blue_black_red,
        opacity: float = 1.0,
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
        show_axis_label : bool
            If True, shows the axis labels. Default is False.
        """
        cartesian = self.mesh.spherical_mesh.to_cartesian()

        window_size = (unit_size[1] * 4, unit_size[0])  # Four subplots horizontally

        scene = pyvista.Plotter(
            theme=pyvista.themes.DocumentTheme(), window_size=window_size, shape=(1, 4)
        )
        scene.set_background(background_color)

        for idx, (name, field) in enumerate(
            zip(["I", "Q", "U", "V"], [self.I, self.Q, self.U, self.V])
        ):
            field = field.flatten(order="F").magnitude
            mesh = pyvista.StructuredGrid(
                cartesian.x.to("meter").magnitude,
                cartesian.y.to("meter").magnitude,
                cartesian.z.to("meter").magnitude
            )
            scene.subplot(0, idx)

            mapping = scene.add_mesh(
                mesh,
                cmap=colormap,
                scalars=field,
                opacity=opacity,
                style="surface",
                show_edges=show_edges,
                show_scalar_bar=False,
            )

            scene.add_axes_at_origin(labels_off=not show_axis_label)
            scene.add_scalar_bar(mapper=mapping.mapper, title=f"{name} field")

        scene.show()
