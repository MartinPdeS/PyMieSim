#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from pydantic.dataclasses import dataclass
import pyvista

from TypedUnit import Length


config_dict = dict(
    arbitrary_types_allowed=True, kw_only=True, slots=True, extra="forbid"
)


@dataclass(config=config_dict, kw_only=True)
class BaseRepresentation:
    """
    Base class for scattering representations.

    Parameters
    ----------
    scatterer : BaseScatterer
        The scatterer object, representing the physical scatterer in the simulation.
    sampling : int
        The number of points used for evaluating the Stokes parameters in spherical coordinates (default is 100).
    distance : float
        The distance from the scatterer at which fields are evaluated (default is 1.0).

    Methods:
        compute_components: A placeholder method intended to be overridden by subclasses for computing specific scattering components.
    """

    scatterer: object
    sampling: int
    distance: Length

    def __post_init__(self):
        fields = self.scatterer._cpp_get_full_farfields(
            sampling=self.sampling, distance=self.distance.to_base_units().magnitude
        )
        self.E_phi, self.E_theta, self.theta, self.phi = fields

        self.compute_components()

    def compute_components(self) -> None:
        """
        Placeholder method for computing scattering components. Intended to be overridden by subclasses.
        """
        raise NotImplementedError("This method should be implemented by subclasses.")

    def get_colormap_limits(self, scalar: numpy.ndarray, symmetric: bool = False):
        if symmetric:
            max_abs = numpy.abs(scalar).max()
            return [-max_abs, max_abs]
        else:
            return None

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

        # Define the vector direction (unit vector along x-axis)
        vector = numpy.array([1, 0, 0])

        # Convert spherical coordinates to Cartesian coordinates
        x, y, z = pyvista.transform_vectors_sph_to_cart(theta, phi, radius, *vector)

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

        # Define the vector direction (unit vector along y-axis)
        vector = numpy.array([0, 1, 0])

        # Convert spherical coordinates to Cartesian coordinates
        x, y, z = pyvista.transform_vectors_sph_to_cart(theta, phi, radius, *vector)

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
