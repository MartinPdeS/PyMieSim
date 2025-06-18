#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from typing import Optional

from PyMieSim.single.representations import Footprint
from MPSPlots.colormaps import blue_black_red
from PyMieSim.units import Quantity, watt, meter, second, farad, volt
from PyMieSim import units
from PyMieSim.single.scatterer.base import BaseScatterer
import pyvista

c = 299792458.0 * meter / second  #: Speed of light in vacuum (m/s).
epsilon0 = 8.854187817620389e-12 * farad / meter  #: Vacuum permittivity (F/m).


class BaseDetector(units.UnitsValidation):
    medium_refractive_index: Optional[Quantity] = 1.0 * units.RIU

    @property
    def max_angle(self) -> Quantity:
        """
        Returns the maximum angle of the detector in radians.
        This is used to determine the angular coverage of the detector.
        """
        return self._cpp_max_angle * units.radian

    @property
    def min_angle(self) -> Quantity:
        """
        Returns the minimum angle of the detector in radians.
        This is used to determine the angular coverage of the detector.
        """
        return self._cpp_min_angle * units.radian

    def _validate_angle_units(cls, value):
        """
        Ensures that gamma_offset, phi_offset, and rotation are Quantity objects with angle units.
        Converts them to numpy arrays after validation.
        """
        if not isinstance(value, Quantity) or not value.check(units.degree):
            raise ValueError(f"{value} must be a Quantity with angle units [degree or radian].")

        return value

    def get_coupling(self, scatterer: BaseScatterer) -> Quantity:
        r"""
        Compute the light coupling between the detector and a scatterer.

        The coupling quantifies the interaction between the field captured by the detector and the scattered field produced by the scatterer. Mathematically, the coupling is calculated as:

        .. math::
            |\iint_{\Omega}  \Phi_{det} \, \Psi_{scat}^* \,  d \Omega|^2

        Where:

        - \( \Phi_{det} \): The capturing field of the detector, representing the sensitivity of the detector to the incoming scattered field.
        - \( \Psi_{scat} \): The scattered field produced by the scatterer.
        - \( \Omega \): The solid angle over which the integration is performed, typically covering the full \( 4\pi \) steradians around the scatterer.
        - \( d\Omega \): The differential solid angle element.

        This integral computes the overlap between the detector's sensitivity pattern and the scattered field, which is then squared to represent the power coupled into the detector.

        Parameters
        ----------
        scatterer : BaseScatterer
            The scatterer object that interacts with the incident light, producing the scattered field.

        Returns
        -------
        Quantity
            The power coupling between the detector and the scatterer, expressed in watts (W). This value represents the amount of scattered power that is captured by the detector.

        Notes
        -----
        - The method internally invokes the appropriate binding method based on the type of scatterer (e.g., Sphere, Cylinder) to calculate the coupling.
        - The coupling depends on both the geometry of the detector and the nature of the scattered field, making it essential for evaluating the efficiency of light collection in scattering experiments.

        Example
        -------
        A common use case is to evaluate how much of the scattered light from a nanoparticle is captured by a photodiode or integrating sphere. The result can be used to estimate the efficiency of light collection for scattering measurements.

        """
        return self._cpp_get_coupling(scatterer) * units.watt

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
        coordinates = numpy.vstack((
            self._cpp_mesh.cartesian.x,
            self._cpp_mesh.cartesian.y,
            self._cpp_mesh.cartesian.z
        ))

        # Wrap the coordinates for PyVista visualization
        points = pyvista.wrap(coordinates.T)

        scalar_field = numpy.asarray(self._cpp_scalar_field).real

        abs_max = abs(scalar_field).max()

        # Add the points to the scene, representing the real part of the scalar field
        mapping = scene.add_points(
            points,
            scalars=scalar_field,
            point_size=20,
            render_points_as_spheres=True,
            cmap=colormap,
            show_scalar_bar=False,
            clim=[-abs_max, abs_max]
        )

        # Create a cone mesh to indicate directional information
        cone_mesh = pyvista.Cone(
            center=coordinates.mean(axis=1) / 2,
            direction=-coordinates.mean(axis=1),
            height=numpy.cos(self.max_angle),
            resolution=100,
            angle=self.max_angle.to('degree').magnitude,
        )

        # Add the cone mesh to the scene with specified color and opacity
        scene.add_mesh(cone_mesh, color='blue', opacity=0.6)

        # Add a scalar bar to the scene for the real part of the field
        scene.add_scalar_bar(mapper=mapping.mapper, title='Collecting Field Real Part')

    def get_poynting_vector(self, scatterer: BaseScatterer, distance: Quantity = 1 * meter) -> float:
        r"""
        Compute the Poynting vector norm, representing the energy flux density of the electromagnetic field.

        The Poynting vector describes the directional energy transfer per unit area for an electromagnetic wave. It is defined as:

        .. math::
            \vec{S} = \epsilon_0 c^2 \, \vec{E} \times \vec{B}

        Where:

        - \( \vec{S} \): Poynting vector (W/m²)
        - \( \epsilon_0 \): Permittivity of free space (F/m)
        - \( c \): Speed of light in vacuum (m/s)
        - \( \vec{E} \): Electric field vector (V/m)
        - \( \vec{B} \): Magnetic field vector (T)

        The cross product of the electric and magnetic field vectors results in the Poynting vector, which represents the flow of electromagnetic energy in space.

        Parameters
        ----------
        scatterer : BaseScatterer
            The scatterer object that interacts with the incident electromagnetic wave, affecting the fields and energy flow.
        distance : Quantity
            The distance at which the Poynting vector is computed.

        Returns
        -------
        Quantity
            The norm of the Poynting vector, which gives the magnitude of the energy flux density in watts per square meter (W/m²).

        Notes
        -----
        The Poynting vector is computed over a 3D mesh of voxels that cover the entire solid angle of \( 4\pi \) steradians. This method calculates the local energy flux at each voxel and returns the norm, which represents the magnitude of energy flow at each point in space around the scatterer.

        The Poynting vector is fundamental in understanding how energy is transmitted through space in the form of electromagnetic waves.

        Example
        -------
        This method is used to assess the distribution of energy around a scatterer. The total energy flow can be obtained by integrating the Poynting vector over the surface enclosing the scatterer.

        """
        Ephi, Etheta = scatterer.get_farfields_array(
            phi=self._cpp_mesh.spherical.phi,
            theta=self._cpp_mesh.spherical.theta,
            r=distance
        )

        E_norm = numpy.sqrt(numpy.abs(Ephi)**2 + numpy.abs(Etheta)**2) * volt / meter

        B_norm = (E_norm / c).to('tesla')

        poynting = epsilon0 * c**2 * E_norm * B_norm

        return poynting.to(watt/meter ** 2)

    def get_energy_flow(self, scatterer: BaseScatterer, distance: Quantity = 1 * meter) -> Quantity:
        r"""
        Calculate the total energy flow (or radiated power) from the scatterer based on the Poynting vector.

        The energy flow is computed using the following relationship between the scattered energy and the incident intensity:

        .. math::
            W_a &= \sigma_{sca} \cdot I_{inc} \\[10pt]
            P &= \int_{A} I \, dA \\[10pt]
            I &= \frac{c n \epsilon_0}{2} \, |E|^2

        Where:

        - \( W_a \): Energy flow (W)
        - \( \sigma_{sca} \): Scattering cross section (m²)
        - \( I_{inc} \): Incident intensity (W/m²)
        - \( P \): Radiated power (W)
        - \( I \): Energy density (W/m²)
        - \( c \): Speed of light in vacuum (m/s)
        - \( n \): Refractive index of the surrounding medium
        - \( \epsilon_0 \): Permittivity of free space (F/m)
        - \( E \): Electric field (V/m)

        The total power is computed by integrating the intensity over the surface area of the scatterer.

        Parameters
        ----------
        scatterer : BaseScatterer
            The scatterer object, which contains information about the scattering properties of the particle, such as geometry and material.
        distance : Quantity
            The distance at which the Poynting vector is computed. It should change the computed total power.

        Returns
        -------
        Quantity
            The total energy flow (radiated power) from the scatterer, expressed in watts.

        Notes
        -----
        This method computes the energy flow by calculating the Poynting vector (which represents the directional energy flux) and summing it over the surface mesh of the scatterer. The final result is the total radiated power.

        """
        poynting_vector = self.get_poynting_vector(scatterer=scatterer, distance=distance)

        dA = self._cpp_mesh._cpp_d_omega * distance ** 2

        total_power = 0.5 * numpy.trapezoid(y=poynting_vector, dx=dA)

        return total_power

