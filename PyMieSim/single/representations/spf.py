#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import List
import pyvista

from PyMieSim.units import ureg
from PyMieSim.utils import spherical_to_cartesian
from PyMieSim.single.full_mesh import FullMesh # Necessary for loading the class, even if not directly used in this file


class SPF():
    r"""
    Compute the scattering phase function (SPF) for the scatterer.

    The scattering phase function describes how the intensity of scattered light varies as a function of the scattering angles \( \theta \) (polar) and \( \phi \) (azimuthal). It is a key quantity in light scattering, as it characterizes the angular distribution of scattered light.

    The scattering phase function is computed as:

    .. math::
        \text{SPF} = \sqrt{ E_{\parallel}(\phi, \theta)^2 + E_{\perp}(\phi, \theta)^2 }

    Where:

    - :math:`E_{\parallel}`: The parallel component of the scattered electric field.
    - :math:`E_{\perp}`: The perpendicular component of the scattered electric field.
    - :math:`\phi` and :math:`\theta`: The azimuthal and polar angles, respectively, which describe the angular position of the scattered light.

    The SPF combines the intensity contributions from both polarization components to describe the total scattered intensity in different directions.

    Parameters
    ----------
    sampling : int
        The number of angular points used to sample the scattering phase function. A higher sampling value increases the angular resolution of the computed SPF.
    distance : Length, optional
        The distance from the scatterer at which the scattering phase function is evaluated. By default, this is set to 1 meter, but it can be adjusted depending on the specific experimental configuration.

    Returns
    -------
    representations.SPF
        An object containing the computed scattering phase function (SPF), representing the angular distribution of scattered light intensity.

    Notes
    -----
    - The scattering phase function is a critical quantity in the study of light scattering, as it provides insight into the directionality of the scattered light. It is commonly used in fields like atmospheric optics, biomedical imaging, and remote sensing.
    - The `sampling` parameter controls the angular resolution of the SPF. Higher sampling values yield more detailed scattering patterns, especially important when capturing small angular features.

    Example
    -------
    You can use this method to compute the scattering phase function of a spherical scatterer, helping to understand the angular distribution of light scattered by the particle.

    Example usage:

    >>> spf = scatterer.get_spf(sampling=500)
    >>> print(spf)

    """
    def __init__(self, setup, sampling: int = 200):
        self.setup = setup
        self.sampling = sampling
        self.SPF, self.mesh = self.setup.get_structured_spf(
            sampling=self.sampling,
            distance=1.0 * ureg.meter
        )

    def plot(
        self,
        unit_size: List[float] = (800, 800),
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
            The size of the plot window in pixels (width, height). Default is (800, 800).
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
        window_size = (unit_size[1], unit_size[0])

        scene = pyvista.Plotter(
            theme=pyvista.themes.DocumentTheme(), window_size=window_size
        )

        scene.set_background(background_color)

        mapping = self._add_to_3d_ax(
            scene=scene,
            colormap=colormap,
            opacity=opacity,
            show_edges=show_edges,
            set_surface=set_surface,
        )

        scene.add_axes_at_origin(labels_off=not show_axis_label)

        scene.add_scalar_bar(mapper=mapping.mapper, title="Scattering Phase Function")

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
        x, y, z = spherical_to_cartesian(
            r=self.SPF,
            phi=self.mesh.spherical_mesh.phi.to("radian").magnitude,
            theta=self.mesh.spherical_mesh.theta.to("radian").magnitude
        )

        mesh = pyvista.StructuredGrid(x, y, z)

        mapping = scene.add_mesh(
            mesh,
            cmap=colormap,
            scalars=self.SPF.T.flatten(),
            opacity=1.0,
            style="surface",
            show_edges=show_edges,
            show_scalar_bar=False,
        )

        return mapping
