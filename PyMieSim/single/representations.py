#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from PyMieSim.single.scatterer import Sphere, CoreShell, Cylinder
    from PyMieSim.single.detector import Photodiode, LPmode, IntegratingSphere


import numpy
from MPSPlots.render3D import SceneList as SceneList3D
from MPSPlots.render2D import SceneList
from dataclasses import dataclass
from PyMieSim.tools.special_functions import spherical_to_cartesian, rotate_on_x


@dataclass
class BaseRepresentation():
    """
    Base class for scattering representations.

    Attributes:
        scatterer (Union[Sphere, CoreShell, Cylinder]): The scatterer object, representing the physical scatterer in the simulation.
        sampling (int): The number of points used for evaluating the Stokes parameters in spherical coordinates (default is 100).
        distance (float): The distance from the scatterer at which fields are evaluated (default is 1.0).

    Methods:
        compute_components: A placeholder method intended to be overridden by subclasses for computing specific scattering components.
    """
    scatterer: Sphere | CoreShell | Cylinder
    sampling: int = 100
    distance: float = 1.0

    def __post_init__(self):
        fields = self.scatterer.Bind.get_full_fields(
            sampling=self.sampling,
            r=self.distance
        )

        self.E_phi, self.E_theta, self.theta, self.phi = fields

        self.compute_components()

    def compute_components(self) -> None:
        """
        Placeholder method for computing scattering components. Intended to be overridden by subclasses.
        """
        raise NotImplementedError("This method should be implemented by subclasses.")


class Stokes(BaseRepresentation):
    r"""
    Represents the scattering far-field in the Stokes representation.

    Inherits from BaseRepresentation and calculates the Stokes parameters which describe the polarization state of light.

    The Stokes parameters (I, Q, U, V) are defined according to their conventional definitions, representing the total intensity,
    difference in intensities between horizontal and vertical polarizations, difference in intensities between two diagonal polarizations,
    and the right and left circular polarizations, respectively.

    | The stokes parameters are:
    |     I : intensity of the fields
    |     Q : linear polarization parallel to incident polarization
    |     U : linear polarization 45 degree to incident polarization
    |     V : Circular polarization

    .. math:
        I &= \big| E_x \big|^2 + \big| E_y \big|^2 \\[10pt]

        Q &= \big| E_x \big|^2 - \big| E_y \big|^2 \\[10pt]

        U &= 2 \mathcal{Re} \big\{ E_x E_y^* \big\} \\[10pt]

        V &= 2 \mathcal{Im} \big\{ E_x E_y^* \big\} \\[10pt]

    Methods:
        compute_components: Computes the Stokes parameters based on the electric field components.
        plot: Visualizes the Stokes parameters on a 3D plot.

    """

    def compute_components(self) -> None:
        """
        Computes the Stokes parameters (I, Q, U, V) based on the electric field components (E_phi and E_theta).

        The method calculates the normalized intensity (I), linear polarizations (Q and U), and circular polarization (V) of the light
        scattered by the particle, using the electric field components in spherical coordinates.

        The Stokes parameters are calculated using the following formulas:
        I = |E_phi|^2 + |E_theta|^2
        Q = |E_phi|^2 - |E_theta|^2
        U = 2 * Re{E_phi * E_theta*}
        V = -2 * Im{E_phi * E_theta*}

        The results are stored as attributes of the instance: I, Q, U, and V.
        """
        intensity = numpy.abs(self.E_phi)**2 + numpy.abs(self.E_theta)**2

        self.I = intensity / numpy.max(intensity)
        self.Q = (numpy.abs(self.E_phi)**2 - numpy.abs(self.E_theta)**2) / intensity
        self.U = (+2 * numpy.real(self.E_phi * self.E_theta.conjugate())) / intensity
        self.V = (-2 * numpy.imag(self.E_phi * self.E_theta.conjugate())) / intensity

    def plot(self) -> SceneList3D:
        """
        Visualizes the Stokes parameters (I, Q, U, V) on a 3D plot.

        Returns:
            - SceneList3D: An object containing the 3D visualization of the Stokes parameters.
        """

        phi_mesh, theta_mesh = numpy.meshgrid(self.phi, self.theta)

        x, y, z = spherical_to_cartesian(
            r=phi_mesh * 0 + 0.5,
            phi=phi_mesh,
            theta=theta_mesh
        )

        figure = SceneList3D()

        for field_name in ['I', 'Q', 'U', 'V']:
            field = getattr(self, field_name)

            ax = figure.append_ax()

            artist = ax.add_mesh(
                x=x,
                y=y,
                z=z,
                scalar_coloring=field,
                colormap='seismic',
                show_edges=False
            )

            ax.add_unit_axis(show_label=False)
            ax.add_colorbar(artist=artist, title=f'{field_name} field')

        return figure


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

    def plot(self) -> SceneList3D:
        """
        Visualizes the far-field pattern in a 3D plot, representing the squared magnitudes of the parallel and perpendicular components
        of the electric field in spherical coordinates.

        Returns:
        - SceneList3D: An object containing the 3D visualization of the far-field pattern.
        """
        phi_mesh, theta_mesh = numpy.meshgrid(self.phi, self.theta)

        x, y, z = spherical_to_cartesian(r=phi_mesh * 0 + 0.5, phi=phi_mesh, theta=theta_mesh)

        figure = SceneList3D()

        for field_name in ['phi', 'theta']:
            field_array = getattr(self, f"E_{field_name}")
            ax = figure.append_ax()

            artist = ax.add_mesh(
                x=x,
                y=y,
                z=z,
                scalar_coloring=field_array.real,
                colormap='seismic',
                show_edges=False,
            )

            ax.add_unit_axis(show_label=False)
            ax.add_colorbar(artist=artist, title=f'{field_name} field [real]')

            if 'phi' in field_name:
                ax.add_unit_phi_vector(radius=1 / 2)
            elif 'theta' in field_name:
                ax.add_unit_theta_vector(radius=1 / 2)

            field_array = getattr(self, f"E_{field_name}")
            ax = figure.append_ax()
            artist = ax.add_mesh(
                x=x,
                y=y,
                z=z,
                scalar_coloring=field_array.imag,
                colormap='seismic',
                show_edges=False,
            )

            ax.add_unit_axis(show_label=False)
            ax.add_colorbar(artist=artist, title=f'{field_name} field [imag]')

            if 'phi' in field_name:
                ax.add_unit_phi_vector(radius=1 / 2)
            elif 'theta' in field_name:
                ax.add_unit_theta_vector(radius=1 / 2)

        return figure


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
        self.SPF = numpy.sqrt(numpy.abs(self.E_phi)**2 + numpy.abs(self.E_theta)**2)

    def plot(self) -> SceneList3D:
        """
        Visualizes the Scattering Phase Function (SPF) on a 3D plot.

        The method maps the SPF values to radial distances from the origin in spherical coordinates, providing a visual representation
        of how light is scattered by the particle in different directions.

        Returns:
        - SceneList3D: An object containing the 3D visualization of the SPF.
        """
        scalar_mesh = self.SPF / self.SPF.max() * 2

        phi_mesh, theta_mesh = numpy.meshgrid(self.phi, self.theta)

        x, y, z = spherical_to_cartesian(
            r=scalar_mesh,
            phi=phi_mesh,
            theta=theta_mesh
        )

        figure = SceneList3D()

        ax = figure.append_ax()
        artist = ax.add_mesh(
            x=x,
            y=y,
            z=z,
            scalar_coloring=scalar_mesh,
            show_edges=False
        )

        ax.add_unit_axis(show_label=False)
        ax.add_colorbar(artist=artist, title='Scattering phase function')

        return figure


@dataclass
class S1S2():
    """
    Represents the S1 and S2 scattering functions, which are components of the scattering matrix.

    Attributes:
        scatterer (Union[Sphere, CoreShell, Cylinder]): The scatterer object.
        sampling (int): Number of points for evaluating the S1 and S2 functions.
        distance (float): Distance at which the fields are evaluated.

    Methods:
        compute_components: Computes the S1 and S2 functions based on the scatterer's properties.
        plot: Visualizes the S1 and S2 functions on a polar plot.
    """
    scatterer: Sphere | CoreShell | Cylinder
    sampling: int = 300
    distance: float = 1.0

    def __post_init__(self):
        self.phi = numpy.linspace(-180, 180, self.sampling)

        self.compute_components()

    def compute_components(self) -> None:
        """
        Computes the S1 and S2 scattering parameters based on the scatterer's properties and the scattering angle phi.

        S1 and S2 are integral parts of the scattering matrix describing the change in polarization state of light upon scattering.

        The method calculates these parameters for a range of phi angles and stores them as the S1 and S2 attributes of the instance.
        """
        self.S1, self.S2 = self.scatterer.Bind.get_s1s2(phi=numpy.deg2rad(self.phi) + numpy.pi / 2)

    def plot(self) -> SceneList:
        """
        Visualizes the S1 and S2 scattering parameters on a polar plot.

        The method creates two polar plots representing the absolute values of S1 and S2 as functions of the scattering angle phi.

        Returns:
        - SceneList: An object containing the polar plots of the S1 and S2 parameters.
        """
        figure = SceneList(unit_size=(3, 3))

        ax_s1 = figure.append_ax(projection='polar', title=r'S_1 parameter')
        ax_s2 = figure.append_ax(projection='polar', title=r'S_2 parameter')

        ax_s1.add_fill_line(
            x=numpy.deg2rad(self.phi),
            y0=numpy.zeros(self.phi.shape),
            y1=numpy.abs(self.S1),
            color='C0',
            line_style='-'
        )

        ax_s2.add_fill_line(
            x=numpy.deg2rad(self.phi),
            y0=numpy.zeros(self.phi.shape),
            y1=numpy.abs(self.S2),
            color='C1',
            line_style='-'
        )

        return figure


@dataclass
class Footprint():
    r"""
    Represents the footprint of the scatterer as detected by various detectors.

    .. math::
        \text{Footprint} = \big| \mathscr{F}^{-1} \big\{ \tilde{ \psi }\
        (\xi, \nu), \tilde{ \phi}_{l,m}(\xi, \nu)  \big\} \
        (\delta_x, \delta_y) \big|^2

    For more information see references in the
    `documentation <https://pymiesim.readthedocs.io/en/latest>`_
    The footprint is defined as:

    Attributes:
        detector (Union[Photodiode, LPmode, IntegratingSphere]): The detector object.
        scatterer (Union[Sphere, CoreShell, Cylinder]): The scatterer object.
        sampling (int): Number of points to evaluate the Stokes parameters in spherical coordinates (default is 500).
        padding_factor (int): Padding factor for the Fourier transform (default is 20).

    Methods:
        compute_footprint: Computes the footprint based on the far-field patterns and detector characteristics.
        plot: Visualizes the computed footprint.
    """
    detector: Photodiode | LPmode | IntegratingSphere
    scatterer: Sphere | CoreShell | Cylinder
    sampling: int = 500
    padding_factor: int = 20

    def __post_init__(self):
        assert self.sampling == 500, "Only a sampling of 500 is available for the moment."

        self.compute_footprint()

    def compute_footprint(self):
        """
        Computes the footprint of the scatterer as detected by the specified detector.

        The footprint is calculated based on the far-field scattering patterns and the characteristics of the detector,
        using a Fourier transform to project the far-field onto the detector plane.

        The computed footprint and the corresponding spatial coordinates are stored as attributes of the instance.
        """
        max_angle = self.detector.max_angle
        n_point = complex(self.sampling)

        phi, theta = numpy.mgrid[
            -max_angle: max_angle: n_point, 0: numpy.pi: n_point
        ]

        max_distance_direct_space = 1 / (numpy.sin(self.detector.max_angle) * self.scatterer.source.k / (2 * numpy.pi))

        x = y = numpy.linspace(-1, 1, self.sampling) * self.sampling / 2 * max_distance_direct_space / self.padding_factor

        _, phi, theta = rotate_on_x(phi + numpy.pi / 2, theta, numpy.pi / 2)

        far_field_para, far_field_perp = self.scatterer.get_farfields_array(
            phi=phi.flatten() + numpy.pi / 2,
            theta=theta.flatten(),
            r=1.0,
            structured=False
        )

        detector_structured_farfield = self.detector.get_structured_scalarfield(sampling=self.sampling)

        perpendicular_projection = detector_structured_farfield * far_field_perp.reshape(theta.shape)

        parallel_projection = detector_structured_farfield * far_field_para.reshape(theta.shape)

        fourier_parallel = self.get_fourier_component(parallel_projection)
        fourier_perpendicular = self.get_fourier_component(perpendicular_projection)

        self.mapping = (fourier_parallel + fourier_perpendicular)
        self.direct_x = x
        self.direct_y = y

    def get_fourier_component(self, scalar: numpy.ndarray) -> numpy.ndarray:
        """
        Computes the Fourier component of a given scalar field.

        This method performs a two-dimensional inverse Fourier transform on the input scalar field, which represents
        a projection (either parallel or perpendicular) of the far-field pattern. It then extracts a central portion
        of the result, effectively applying a padding factor to increase the resolution of the Fourier transform.

        Parameters:
        - scalar (numpy.ndarray): A two-dimensional numpy array representing the scalar field of which the Fourier component
          is to be computed. This field could represent either the parallel or perpendicular projection of the far-field
          pattern onto the detector plane.

        Returns:
        - numpy.ndarray: A two-dimensional numpy array representing the computed Fourier component. This array is a square
          section, extracted from the center of the full Fourier transform, with dimensions determined by the original
          sampling rate and the padding factor of the instance. The values in the array represent the intensity distribution
          of the light in the detector plane, providing insights into the spatial characteristics of the scattering pattern.

        The method uses numpy's fft module to perform the Fourier transform, applying a padding factor to the input to
        achieve a higher resolution in the Fourier domain. The resulting Fourier transform is then squared and fftshifted
        to center the zero-frequency component, and a central portion is extracted to match the intended output size.
        """
        # Calculate the target size based on the sampling and padding factor, and the indices for the central portion extraction.
        total_size = self.sampling * self.padding_factor
        offset = (total_size - self.sampling) // 2

        # Apply zero-padding to the scalar field to increase the resolution of the Fourier transform.
        padded_scalar = numpy.pad(scalar, pad_width=((offset, offset), (offset, offset)), mode='constant', constant_values=0)

        # Perform the two-dimensional inverse Fourier transform on the padded scalar field.
        fourier_transformed = numpy.fft.ifft2(padded_scalar)

        # Compute the squared magnitude and center the zero-frequency component.
        fourier_magnitude_squared = numpy.abs(numpy.fft.fftshift(fourier_transformed))**2

        # Extract the central portion corresponding to the original sampling rate adjusted by the padding factor.
        central_portion = fourier_magnitude_squared[offset:-offset, offset:-offset]

        return central_portion

    def plot(self) -> SceneList:
        """
        Visualizes the footprint of the scatterer on the detector plane.

        Creates a 2D plot representing the intensity distribution of the scattered light as it would be detected by the specified detector,
        providing insight into the spatial characteristics of the scattering pattern.

        Returns:
        - SceneList: An object containing the 2D plot of the scatterer's footprint.
        """
        figure = SceneList(unit_size=(6, 6))

        ax = figure.append_ax(
            title='Scatterer Footprint',
            x_label=r'Offset distance in X-axis [$\mu$m]',
            y_label=r'Offset distance in Y-axis [$\mu$m]',
        )

        artist = ax.add_mesh(
            x=self.direct_y * 1e6,
            y=self.direct_x * 1e6,
            scalar=self.mapping,
        )

        ax.add_colorbar(artist=artist, colormap='gray')

        return figure


# -
