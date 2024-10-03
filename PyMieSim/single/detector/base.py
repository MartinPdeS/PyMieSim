#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from typing import Optional
from pydantic.dataclasses import dataclass
from pydantic import field_validator
from typing import Union, Optional, Tuple

from PyMieSim.single.representations import Footprint
from PyMieSim.binary.Fibonacci import FibonacciMesh as CPPFibonacciMesh  # has to be imported as extension  # noqa: F401
from MPSPlots.colormaps import blue_black_red
from PyMieSim.units import Quantity, degree, AU

from PyMieSim.single.scatterer.base import BaseScatterer
import pyvista

c = 299792458.0  #: Speed of light in vacuum (m/s).
epsilon0 = 8.854187817620389e-12  #: Vacuum permittivity (F/m).

config_dict = dict(
    kw_only=True,
    slots=True,
    extra='forbid',
    arbitrary_types_allowed=True
)

@dataclass(config=config_dict)
class BaseDetector():
    """
    Base class for different types of detectors with methods for setup, rotation,
    calculating light coupling, and generating footprint and 3D plots.
    """
    mode_number: str
    NA: Quantity
    gamma_offset: Quantity
    phi_offset: Quantity
    mean_coupling: bool
    coherent: bool
    sampling: Optional[Quantity] = 200 * AU
    polarization_filter: Optional[Quantity] = numpy.nan * degree
    rotation: Optional[Quantity] = 90 * degree

    @field_validator('polarization_filter', mode='before')
    def _validate_polarization(cls, value):
        """
        Ensures that gamma_offset, phi_offset, polarization_filter, and rotation are Quantity objects with angle units.
        Converts them to numpy arrays after validation.
        """
        if value is None:
            value = numpy.nan * degree

        if not value.check(degree):
            raise ValueError(f"{value} must have angle units (degree or radian).")

        return value

    @field_validator('gamma_offset', 'phi_offset', 'rotation', mode='before')
    def _validate_angle_quantity(cls, value):
        """
        Ensures that gamma_offset, phi_offset, and rotation are Quantity objects with angle units.
        Converts them to numpy arrays after validation.
        """
        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with angle units.")

        if not value.check(degree):
            raise ValueError(f"{value} must have angle units (degree or radian).")

        return value


    def __post_init__(self):
        self.initialize()

    def coupling(self, scatterer: BaseScatterer) -> float:
        r"""
        Calculate the light coupling between the detector and a scatterer.

        .. math::
            |\iint_{\Omega}  \Phi_{det} \,\, \Psi_{scat}^* \,  d \Omega|^2

        | Where:
        |   :math:`\Phi_{det}` is the capturing field of the detector and
        |   :math:`\Psi_{scat}` is the scattered field.

        Args:
            scatterer (BaseScatterer): The scatterer object.

        Returns:
            float: The coupling in watts.
        """
        return getattr(self.binding, "Coupling" + type(scatterer).__name__)(scatterer.binding)

    def get_footprint(self, scatterer: BaseScatterer) -> Footprint:
        r"""
        Generate the footprint of the scattered light coupling with the detector.

        .. math::
            \big| \mathscr{F}^{-1} \big\{ \tilde{ \psi } (\xi, \nu),\
                   \tilde{ \phi}_{l,m}(\xi, \nu)  \big\}
                   (\delta_x, \delta_y) \big|^2

        | Where:
        |   :math:`\Phi_{det}` is the capturing field of the detector and
        |   :math:`\Psi_{scat}` is the scattered field.

        Args:
            scatterer (BaseScatterer): The scatterer object.

        Returns:
            Footprint: The scatterer footprint with this detector.
        """
        return Footprint(scatterer=scatterer, detector=self)

    def plot(self,
            unit_size: Tuple[float, float] = (600, 600),
            background_color: str = 'black',
            show_axis_label: bool = False,
            colormap: str = blue_black_red) -> None:
        """
        Plot the real and imaginary parts of the scattered fields in a 3D scene.

        This method creates a 3D plot of the scattered fields, displaying the real part of the field
        as colored points in the scene. It also adds a translucent sphere and a directional cone
        to represent additional geometrical information.

        Parameters
        ----------
        unit_size : Tuple[float, float]
            The size of the plot window (width, height). Default is (600, 600).
        background_color : str
            The background color of the plot. Default is 'black'.
        show_axis_label : bool
            If True, axis labels will be shown. Default is False.
        colormap : str
            The colormap to use for the scalar field. Default is 'blue_black_red'.

        Returns:
            None: This method does not return a value. It displays the 3D plot.
        """
        # Stack the mesh coordinates
        coordinates = numpy.row_stack((
            self.binding.mesh.x,
            self.binding.mesh.y,
            self.binding.mesh.z
        ))

        # Create a 3D plotting scene with the specified window size and theme
        window_size = (unit_size[1], unit_size[0])

        scene = pyvista.Plotter(
            theme=pyvista.themes.DarkTheme(),
            window_size=window_size
        )

        # Set the background color
        scene.set_background(background_color)

        self._add_to_3d_ax(scene=scene, colormap=colormap)

        # Optionally add axis labels
        scene.add_axes_at_origin(labels_off=not show_axis_label)

        # Add a translucent sphere to the scene
        sphere = pyvista.Sphere(radius=1)
        scene.add_mesh(sphere, opacity=0.3)

        # Display the scene
        scene.show()

    def _add_to_3d_ax(self, scene: pyvista.Plotter, colormap: str = blue_black_red) -> None:
        """
        Adds the scalar field and a directional cone to a 3D PyVista plotting scene.

        This method adds points representing the real part of the scalar field to the given 3D scene,
        along with a cone mesh to indicate directional information. It also includes a scalar bar
        to display the values of the scalar field.

        Parameters
        ----------
        scene : pyvista.Plotter
            The PyVista plotting scene where the elements will be added.
        colormap : str
            The colormap to use for the scalar field visualization.

        Returns:
            None: This method does not return a value. It modifies the provided scene.
        """
        # Stack the mesh coordinates into a single array
        coordinates = numpy.row_stack((
            self.binding.mesh.x,
            self.binding.mesh.y,
            self.binding.mesh.z
        ))

        # Wrap the coordinates for PyVista visualization
        points = pyvista.wrap(coordinates.T)

        # Add the points to the scene, representing the real part of the scalar field
        mapping = scene.add_points(
            points,
            scalars=numpy.real(self.binding.scalar_field),
            point_size=30,
            render_points_as_spheres=True,
            cmap=colormap,
            show_scalar_bar=False
        )

        # Create a cone mesh to indicate directional information
        cone_mesh = pyvista.Cone(
            center=coordinates.mean(axis=1) / 2,
            direction=-coordinates.mean(axis=1),
            height=numpy.cos(self.max_angle),
            resolution=100,
            angle=numpy.rad2deg(self.max_angle)
        )

        # Add the cone mesh to the scene with specified color and opacity
        scene.add_mesh(cone_mesh, color='blue', opacity=0.3)

        # Add a scalar bar to the scene for the real part of the field
        scene.add_scalar_bar(mapper=mapping.mapper, title='Collecting Field Real Part')


    def get_poynting_vector(self, scatterer: BaseScatterer) -> float:
        r"""

        Method return the Poynting vector norm defined as:

        .. math::
            \vec{S} = \epsilon c^2 \vec{E} \times \vec{B}

        Parameters :
            Mesh : Number of voxel in the 4 pi space to compute energy flow.

        """
        Ephi, Etheta = scatterer.get_farfields_array(phi=self.binding.mesh.phi, theta=self.binding.mesh.theta, r=1.)

        E_norm = numpy.sqrt(numpy.abs(Ephi)**2 + numpy.abs(Etheta)**2)

        B_norm = E_norm / c

        poynting = epsilon0 * c**2 * E_norm * B_norm

        return poynting

    def get_energy_flow(self, scatterer: BaseScatterer) -> float:
        r"""
        Returns energy flow defined as:

        .. math::
            W_a &= \sigma_{sca} * I_{inc} \\[10pt]
            P &= \int_{A} I dA \\[10pt]
            I &= \frac{c n \epsilon_0}{2} |E|^2 \\[10pt]

        | With:
        |     I : Energy density
        |     n  : Refractive index of the medium
        |     :math:`\epsilon_0` : Vaccum permitivity
        |     E  : Electric field
        |     \sigma_{sca}: Scattering cross section.

        More info on wikipedia link (see ref[6]).

        :param      Mesh:  The mesh
        :type       Mesh:  FibonacciMesh

        :returns:   The energy flow.
        :rtype:     float
        """

        poynting = self.get_poynting_vector(scatterer=scatterer)

        total_power = 0.5 * numpy.sum(poynting) * self.binding.mesh.d_omega

        return total_power