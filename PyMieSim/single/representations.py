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
    scatterer: Sphere | CoreShell | Cylinder
    """ The scatterer instance """
    sampling: int = 100
    """ Number of point to evaluate the Stokes parameters in spherical coord """
    distance: float = 1.0
    """ Distance at which evaluate the fields """

    def __post_init__(self):
        fields = self.scatterer.Bind.get_full_fields(
            sampling=self.sampling,
            r=self.distance
        )

        self.E_phi, self.E_theta, self.theta, self.phi = fields

        self.compute_components()


class Stokes(BaseRepresentation):
    r"""
    Class representing scattering Far-field in the Stokes representation.
    The parameters are defined as in this wikipedia article: https://en.wikipedia.org/wiki/Stokes_parameters

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

    """

    def compute_components(self) -> None:
        """
        Calculates the Stokes components.

        :returns:   No returns.
        :rtype:     None
        """
        intensity = numpy.abs(self.E_phi)**2 + numpy.abs(self.E_theta)**2

        self.I = intensity / numpy.max(intensity)
        self.Q = (numpy.abs(self.E_phi)**2 - numpy.abs(self.E_theta)**2) / intensity
        self.U = (+2 * numpy.real(self.E_phi * self.E_theta.conjugate())) / intensity
        self.V = (-2 * numpy.imag(self.E_phi * self.E_theta.conjugate())) / intensity

    def plot(self, show_ource=True, show_axes=True):
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
    Class representing scattering Far-field in a spherical coordinate representation.
    The Far-fields are defined as:

    .. math::
        \text{Fields} = E_{||}(\phi,\theta)^2, E_{\perp}(\phi,\theta)^2

    """

    def compute_components(self) -> None:
        """
        Calculates the fields components. They are already computed.

        :returns:   No returns.
        :rtype:     None
        """
        return

    def plot(self) -> SceneList3D:
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
    Class representing scattering phase function of SPF in short.
    The SPF is defined as:

    .. math::
        \text{SPF} = E_{\parallel}(\phi,\theta)^2 + E_{\perp}(\phi,\theta)^2

    """

    def compute_components(self) -> None:
        """
        Calculates the Stokes components.

        :returns:   No returns.
        :rtype:     None
        """
        self.SPF = numpy.sqrt(numpy.abs(self.E_phi)**2 + numpy.abs(self.E_theta)**2)

    def plot(self) -> SceneList3D:
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
    """ Class representing S1 and S2 function. S1 and S2 are defined as: """
    scatterer: Sphere | CoreShell | Cylinder
    """ The scatterer instance """
    sampling: int = 300
    """ Number of point to evaluate the Stokes parameters in spherical coord """
    distance: float = 1.0
    """ Distance at which evaluate the fields """

    def __post_init__(self):
        self.phi = numpy.linspace(-180, 180, self.sampling)

        self.compute_components()

    def compute_components(self) -> None:
        self.S1, self.S2 = self.scatterer.Bind.get_s1s2(phi=numpy.deg2rad(self.phi) + numpy.pi / 2)

    def plot(self) -> SceneList:
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
    detector: Photodiode | LPmode | IntegratingSphere
    """ The detector instance """
    scatterer: Sphere | CoreShell | Cylinder
    """ The scatterer instance """
    sampling: int = 500
    """ Number of point to evaluate the Stokes parameters in spherical coord. At the moment only 500 is available"""
    padding_factor: int = 20
    """ Padding factor for the fourier transform """

    r"""
    Class representing footprint of the scatterer.
    The footprint usually depend on the scatterer and the detector.
    For more information see references in the
    `documentation <https://pymiesim.readthedocs.io/en/latest>`_
    The footprint is defined as:

    .. math::
        \text{Footprint} = \big| \mathscr{F}^{-1} \big\{ \tilde{ \psi }\
        (\xi, \nu), \tilde{ \phi}_{l,m}(\xi, \nu)  \big\} \
        (\delta_x, \delta_y) \big|^2

    """

    def __post_init__(self):
        assert self.sampling == 500, "Only a sampling of 500 is available for the moment."

        self.compute_footprint()

    def compute_footprint(self):
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
        total_size = self.sampling * self.padding_factor

        start = int(total_size / 2 - numpy.floor(self.sampling / 2))
        end = int(total_size / 2 + numpy.ceil(self.sampling / 2))

        fourier = numpy.fft.ifft2(scalar, s=[total_size, total_size])

        fourier = numpy.abs(numpy.fft.fftshift(fourier))**2

        return fourier[start: end, start: end]

    def plot(self) -> SceneList:
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
