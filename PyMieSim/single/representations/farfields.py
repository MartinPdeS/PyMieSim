#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import pyvista
from typing import List
from MPSPlots.colormaps import blue_black_red

from PyMieSim.units import ureg, Length


class FarFields():
    r"""
    Compute the far-field scattering pattern for the scatterer.

    The far fields describe the behavior of the scattered electromagnetic waves at a large distance from the scatterer, where the waves can be approximated as planar.
    The computed far fields represent the scattered electric field components in the directions parallel and perpendicular to the plane of incidence.

    The far fields are computed as:

    .. math::
        \text{Fields} = E_{||}(\phi, \theta)^2, \, E_{\perp}(\phi, \theta)^2

    These components represent the intensities of the scattered electric field in the parallel (\( E_{||} \)) and perpendicular (\( E_{\perp} \)) directions as functions of the spherical angles \( \phi \) (azimuthal) and \( \theta \) (polar).

    The fields are expressed up to a constant phase factor:

    .. math::
        \exp{\left(-i k r \right)}

    Where:

    - :math:`k`: The wave number, related to the wavelength of the incident light.
    - :math:`r`: The distance from the scatterer (assumed to be large in the far-field approximation).

    Parameters
    ----------
    sampling : int
        The number of angular points used to sample the far-field pattern. A higher sampling value increases the angular resolution of the computed far-field intensities.
    distance : Length, optional
        The distance from the scatterer at which the far fields are evaluated. By default, this is set to 1 meter, but it can be adjusted to suit the specific experimental configuration.

    Returns
    -------
    representations.FarField
        An object containing the computed far-field components (parallel and perpendicular intensities), which represent the angular distribution of the scattered light in the far field.

    Notes
    -----
    - The far-field approximation assumes that the distance from the scatterer is large enough that the curvature of the wavefronts can be neglected, and the fields can be treated as planar.
    - This method is useful for determining the angular distribution of the scattered light, which is critical in applications such as radar cross-section analysis, optical scattering experiments, and remote sensing.


    """
    def __init__(self,  setup: object, sampling: int = 200, distance: Length = 1.0 * ureg.meter):
        self.setup = setup
        self.sampling = sampling
        self.distance = distance

        self.E_phi, self.E_theta, self.mesh = self.setup.get_farfields(
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
        show_axis_label : bool
            If True, shows the axis labels. Default is False.
        """
        cartesian = self.mesh.spherical_mesh.to_cartesian()

        window_size = (unit_size[1] * 4, unit_size[0])  # Two subplots horizontally

        scene = pyvista.Plotter(
            theme=pyvista.themes.DocumentTheme(), window_size=window_size, shape=(1, 4)
        )
        scene.set_background(background_color)

        repr = [self.E_phi.real, self.E_phi.imag, self.E_theta.real, self.E_theta.imag]
        repr_label = ["phi real", "phi imag", "theta real", "theta imag"]

        for idx, (label, field) in enumerate(zip(repr_label, repr)):
            field = field.flatten(order="F").magnitude
            mesh = pyvista.StructuredGrid(
                cartesian.x.to("meter").magnitude,
                cartesian.y.to("meter").magnitude,
                cartesian.z.to("meter").magnitude
            )
            scene.subplot(0, idx)

            max_abs = numpy.abs(field).max()
            colormap_limits = [-max_abs, max_abs]

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
                self.add_theta_vector_to_3d_plot(scene=scene, radius=1.1, color="White")

            if "phi" in label:
                self.add_phi_vector_to_3d_plot(scene=scene, radius=1.1, color="White")

            scene.add_axes_at_origin(labels_off=not show_axis_label)
            scene.add_scalar_bar(mapper=mapping.mapper, title=f"{label} field")

        scene.show()

    def add_theta_vector_to_3d_plot(
        self,
        scene: pyvista.Plotter,
        n_points: int = 20,
        opacity: float = 1.0,
        radius: float = 1.0,
        color: str = "black",
    ) -> None:
        """
        Adds a vector field to the 3D plot, representing vectors in the theta direction.

        Parameters
        ----------
        scene : pyvista.Plotter
            The 3D plotting scene to which the vectors will be added.
        n_points : int
            Number of points to generate along the theta and phi directions. Default is 100.
        opacity : float
            Opacity of the vectors. Default is 1.0.
        radius : float
            Radius at which to place the vectors. Default is 1.0.
        color : str
            Color of the vectors. Default is 'black'.
        """
        theta = numpy.linspace(0, 360, n_points)
        phi = numpy.linspace(180, 0, n_points)

        # Convert spherical coordinates to Cartesian coordinates
        x, y, z = pyvista.transform_vectors_sph_to_cart(theta=theta, phi=phi, r=radius, u=1, v=0, w=0)

        # Combine the Cartesian coordinates into a vector array
        vector_field = numpy.c_[x.ravel(), y.ravel(), z.ravel()]

        # Create a structured grid from spherical coordinates
        spherical_grid = pyvista.grid_from_sph_coords(theta, phi, radius)
        spherical_grid.point_data["component"] = vector_field * 0.1

        # Generate glyphs (arrows) for the vectors
        glyphs = spherical_grid.glyph(
            orient="component", scale="component", tolerance=0.005
        )

        # Add the vector glyphs to the scene
        scene.add_mesh(glyphs, color=color, opacity=opacity)

    def add_phi_vector_to_3d_plot(
        self,
        scene: pyvista.Plotter,
        n_points: int = 20,
        opacity: float = 1.0,
        radius: float = 1.0,
        color: str = "black",
    ) -> None:
        """
        Adds a vector field to the 3D plot, representing vectors in the phi direction.

        Parameters
        ----------
        scene : pyvista.Plotter
            The 3D plotting scene to which the vectors will be added.
        n_points : int
            Number of points to generate along the theta and phi directions. Default is 100.
        opacity : float
            Opacity of the vectors. Default is 1.0.
        radius : float
            Radius at which to place the vectors. Default is 1.0.
        color : str
            Color of the vectors. Default is 'black'.
        """
        theta = numpy.linspace(0, 360, n_points)
        phi = numpy.linspace(180, 0, n_points)

        # Convert spherical coordinates to Cartesian coordinates
        x, y, z = pyvista.transform_vectors_sph_to_cart(theta=theta, phi=phi, r=radius, u=0, v=1, w=0)

        # Combine the Cartesian coordinates into a vector array
        vector_field = numpy.c_[x.ravel(), y.ravel(), z.ravel()]

        # Create a structured grid from spherical coordinates
        spherical_grid = pyvista.grid_from_sph_coords(theta, phi, radius)
        spherical_grid.point_data["component"] = vector_field * 0.1

        # Generate glyphs (arrows) for the vectors
        glyphs = spherical_grid.glyph(
            orient="component", scale="component", tolerance=0.005
        )

        # Add the vector glyphs to the scene
        scene.add_mesh(glyphs, color=color, opacity=opacity)
