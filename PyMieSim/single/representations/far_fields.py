#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import pyvista
from pydantic.dataclasses import dataclass
from typing import List

from MPSPlots.colormaps import blue_black_red
from PyMieSim.single.representations.base import BaseRepresentation
from PyMieSim.utils import config_dict, spherical_to_cartesian


@dataclass(config=config_dict, kw_only=True)
class FarField(BaseRepresentation):
    r"""
    Represents the scattering far-field in spherical coordinates.

    Inherits from BaseRepresentation and visualizes the far-field pattern characterized by the perpendicular and parallel components
    of the electric field in spherical coordinates.

    .. math::
        \text{Fields} = E_{||}(\phi,\theta)^2, E_{\perp}(\phi,\theta)^2

    Methods:
        compute_components: Calculates the field components. This implementation is a placeholder, as the components are precomputed.
        plot: Visualizes the far-field pattern in a 3D plot.

    """

    def compute_components(self) -> None:
        """
        Placeholder method in FarField class. Does not perform any computation as field components are precomputed.

        This method is intended to be consistent with the structure of BaseRepresentation but does not need to modify or compute
        any attributes for FarField instances.
        """
        return

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
        Visualizes the Far field (in phi and theta vector projections) on a 3D plot.

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

        window_size = (unit_size[1] * 4, unit_size[0])  # Two subplots horizontally

        scene = pyvista.Plotter(
            theme=pyvista.themes.DocumentTheme(), window_size=window_size, shape=(1, 4)
        )
        scene.set_background(background_color)

        repr = [self.E_phi.real, self.E_phi.imag, self.E_theta.real, self.E_theta.imag]
        repr_label = ["phi real", "phi imag", "theta real", "theta imag"]

        for idx, (label, field) in enumerate(zip(repr_label, repr)):
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
            if "theta" in label:
                self.add_theta_vector_to_3d_plot(scene=scene, radius=0.6)

            if "phi" in label:
                self.add_phi_vector_to_3d_plot(scene=scene, radius=0.6)

            scene.add_axes_at_origin(labels_off=not show_axis_label)
            scene.add_scalar_bar(mapper=mapping.mapper, title=f"{label} field")

        scene.show()
