#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from pydantic.dataclasses import dataclass
from typing import List
import pyvista

from PyMieSim.single.representations.base import BaseRepresentation
from PyMieSim.utils import config_dict, spherical_to_cartesian


@dataclass(config=config_dict, kw_only=True)
class SPF(BaseRepresentation):
    r"""
    Represents the Scattering Phase Function (SPF).

    Inherits from BaseRepresentation and computes the SPF, which is a measure of how light is scattered by a particle at different angles.

    .. math::
        \text{SPF} = E_{\parallel}(\phi,\theta)^2 + E_{\perp}(\phi,\theta)^2

    Methods:
        compute_components: Computes the SPF based on the electric field components.
        plot: Visualizes the SPF on a 3D plot.

    """

    def compute_components(self) -> None:
        """
        Computes the Scattering Phase Function (SPF) based on the electric field components (E_phi and E_theta).

        The SPF is calculated as the square root of the sum of the squared magnitudes of the electric field components, representing
        the total scattering intensity distribution as a function of angles.

        The result is stored as the SPF attribute of the instance.
        """
        self.SPF = numpy.sqrt(numpy.abs(self.E_phi) ** 2 + numpy.abs(self.E_theta) ** 2)

    def plot(
        self,
        unit_size: List[float] = (400, 400),
        background_color: str = "white",
        show_edges: bool = False,
        colormap: str = "viridis",
        opacity: float = 1.0,
        set_surface: bool = True,
        show_axis_label: bool = False,
    ) -> None:
        """
        Visualizes the scattering phase function on a 3D plot.

        This method creates a 3D visualization of the scattering phase function (SPF). It allows customization
        of the plot's appearance, including the colormap, mesh opacity, and whether or not to display mesh edges
        and axis labels.

        Parameters
        ----------
        unit_size : List[float]
            The size of the plot window in pixels (width, height). Default is (400, 400).
        background_color : str
            The background color of the plot. Default is 'white'.
        show_edges : bool
            If True, displays the edges of the mesh. Default is False.
        colormap : str
            The colormap to use for scalar mapping. Default is 'viridis'.
        opacity : float
            The opacity of the mesh. Default is 1.0.
        set_surface : bool
            If True, the surface represents the scaled SPF; if False, a unit sphere is used. Default is True.
        show_axis_label : bool
            If True, shows the axis labels. Default is False.
        """
        # Define the window size based on the unit size provided
        window_size = (unit_size[1], unit_size[0])  # One subplot

        # Create a PyVista plotting scene with the specified theme and window size
        scene = pyvista.Plotter(
            theme=pyvista.themes.DocumentTheme(), window_size=window_size
        )

        # Set the background color of the scene
        scene.set_background(background_color)

        # Add the 3D axis-aligned plot to the scene using the specified settings
        mapping = self._add_to_3d_ax(
            scene=scene,
            colormap=colormap,
            opacity=opacity,
            show_edges=show_edges,
            set_surface=set_surface,
        )

        # Optionally add axis labels
        scene.add_axes_at_origin(labels_off=not show_axis_label)

        # Add a scalar bar to the scene to represent the scattering phase function
        scene.add_scalar_bar(mapper=mapping.mapper, title="Scattering Phase Function")

        # Display the scene
        scene.show()

    def _add_to_3d_ax(
        self,
        scene: pyvista.Plotter,
        set_surface: bool = False,
        show_edges: bool = False,
        colormap: str = "viridis",
        opacity: float = 1.0,
    ) -> None:
        """
        Adds a 3D surface plot to the provided PyVista scene based on the scattering phase function (SPF).

        This method generates a 3D surface plot of the SPF using spherical coordinates, and adds it to the scene.
        The surface can either represent the actual SPF or a normalized unit sphere, depending on the `set_surface` flag.
        The appearance of the surface can be customized using various parameters.

        Parameters
        ----------
        scene : pyvista.Plotter
            The PyVista plotting scene where the surface will be added.
        set_surface : bool
            If True, the surface will represent the scaled SPF; if False, a unit sphere is used. Default is True.
        show_edges : bool
            If True, edges of the mesh will be displayed. Default is False.
        colormap : str
            The colormap to use for visualizing the scalar field. Default is 'viridis'.
        opacity : float
            The opacity of the surface mesh. Default is 1.0.
        """
        # Create mesh grids for phi and theta
        phi_mesh, theta_mesh = numpy.meshgrid(self.phi, self.theta)

        # Normalize the scattering phase function (SPF) for visualization
        scalar = self.SPF / self.SPF.max() * 2

        # Determine the coordinates based on whether the surface represents the SPF or a unit sphere
        if set_surface:
            x, y, z = spherical_to_cartesian(r=scalar, phi=phi_mesh, theta=theta_mesh)
        else:
            x, y, z = spherical_to_cartesian(
                r=numpy.ones(phi_mesh.shape) * 0.5, phi=phi_mesh, theta=theta_mesh
            )

        # Create a structured grid from the calculated coordinates
        mesh = pyvista.StructuredGrid(x, y, z)

        # Add the surface mesh to the scene
        mapping = scene.add_mesh(
            mesh,
            cmap=colormap,
            scalars=scalar.flatten(order="F"),
            opacity=opacity,
            style="surface",
            show_edges=show_edges,
            show_scalar_bar=False,
        )

        return mapping
