#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy


from MPSPlots.render3D import Scene3D
from MPSPlots.render2D import SceneList
import PyMieSim
from PyMieSim.tools.special_functions import spherical_to_cartesian, rotate_on_x


class Stokes():
    # https://en.wikipedia.org/wiki/Stokes_parameters
    r"""
    Class representing scattering Far-field in the Stokes
    representation.
    | The stokes parameters are:
    |     I : intensity of the fields
    |     Q : linear polarization parallel to incident polarization
    |     U : linear polarization 45 degree to incident polarization
    |     V : Circular polarization

    .. math:
        I &= \\big| E_x \big|^2 + \\big| E_y \\big|^2

        Q &= \\big| E_x \big|^2 - \\big| E_y \\big|^2

        U &= 2 \\mathcal{Re} \\big\{ E_x E_y^* \\big\}

        V &= 2 \\mathcal{Im} \\big\{ E_x E_y^* \\big\}

    Parameters
    ----------
    parent_scatterer : :class:`Scatterer`
        The scatterer parent.
    sampling : :class:`int`
        samplingber of point to evaluate the Stokes parameters in spherical coord.
    distance : :class:`float`
        distance at which we evaluate the Stokes parameters.

    """

    def __init__(self, parent_scatterer, sampling: int = 100, distance: float = 1.):
        self.parent_scatterer = parent_scatterer

        self.E_phi, self.E_theta, self.theta, self.phi = parent_scatterer.Bind.get_full_fields(sampling=sampling, r=1)

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

        figure = Scene3D(shape=(1, 4), window_size=[2500, 800])

        enumeration = zip(
            [(0, 0), (0, 1), (0, 2), (0, 3)],
            [self.I, self.Q, self.U, self.V],
            ['I', 'Q', 'U', 'V']
        )

        for plot_number, field, name in enumeration:
            figure.add_mesh(
                plot_number=plot_number,
                x=x,
                y=y,
                z=z,
                scalars=field.T,
                color_map='seismic',
                scalar_bar_args={'title': f'{name} Amplitude'}
            )

            figure.add_unit_axes_to_ax(plot_number=plot_number)
            figure.add_text_to_axes(plot_number=plot_number, text=f'{name} field')

        return figure


class SPF():
    r"""
    Class representing scattering phase function of SPF in short.
    The SPF is defined as:
    .. math::
        \\text{SPF} = E_{\\parallel}(\\phi,\\theta)^2 + E_{\\perp}(\\phi,\\theta)^2

    Parameters
    ----------
    parent_scatterer : :class:`Scatterer`
        The scatterer parent.
    sampling : :class:`int`
        samplingber of point to evaluate the SPF in spherical coord.
    distance : :class:`float`
        distance at which we evaluate the SPF.

    """

    def __init__(self, parent_scatterer, sampling: int = 100, distance: float = 1.):

        self.parent_scatterer = parent_scatterer

        self.E_phi, self.E_theta, self.theta, self.phi = parent_scatterer.Bind.get_full_fields(sampling=sampling, r=1)

        self.SPF = numpy.sqrt(numpy.abs(self.E_phi)**2 + numpy.abs(self.E_theta)**2)

    def plot(self, show_source=True, show_axes=True):

        scalar_mesh = self.SPF / self.SPF.max() * 2

        phi_mesh, theta_mesh = numpy.meshgrid(self.phi, self.theta)

        x, y, z = spherical_to_cartesian(
            r=scalar_mesh,
            phi=phi_mesh,
            theta=theta_mesh
        )

        figure = Scene3D(shape=(1, 1), window_size=[1800, 1000])

        plot_number = (0, 0)

        figure.add_mesh(
            plot_number=plot_number,
            x=x,
            y=y,
            z=z,
            scalars=scalar_mesh.T,
            scalar_bar_args={'title': 'intensity'}
        )

        figure.add_unit_axes_to_ax(plot_number=plot_number)
        figure.add_text_to_axes(plot_number=plot_number, text='Scattering phase function')

        return figure


class S1S2():
    r"""Dict subclass representing S1 and S2 function.
    S1 and S2 are defined as:

    Parameters
    ----------
    parent_scatterer : :class:`Scatterer`
        The scatterer parent.
    sampling : :class:`int`
        samplingber of point to evaluate the S1 and S2 in spherical coord.

    """
    def __init__(self, parent_scatterer, Phi: numpy.ndarray = None, sampling: int = None):
        self.parent_scatterer = parent_scatterer

        if sampling is None:
            sampling = 200

        if Phi is None:
            Phi = numpy.linspace(-180, 180, sampling)

        self.S1, self.S2 = parent_scatterer.Bind.get_s1s2(phi=numpy.deg2rad(Phi) + numpy.pi / 2)
        self.phi = Phi

    def plot(self) -> SceneList:
        figure = SceneList(unit_size=(3, 3))

        zero = 0 * numpy.abs(self.S1)

        enumeration = zip(
            [0, 1],
            [self.S1, self.S2],
            ['S1', 'S2'],
            ['C0', 'C1']
        )

        for col, s_param, name, color in enumeration:
            ax = figure.append_ax(projection='polar', title=f'{name} parameter')

            ax.add_fill_line(
                x=numpy.deg2rad(self.phi),
                y0=zero,
                y1=numpy.abs(s_param),
                color=color,
                line_style='-'
            )

        return figure


class FarField():
    r"""
    Class representing scattering Far-field in a spherical
    coordinate representation.
    The Far-fields are defined as:

    .. math::
        \\text{Fields} = E_{||}(\\phi,\\theta)^2, E_{\\perp}(\\phi,\\theta)^2

    Parameters
    ----------
    parent_scatterer : :class:`Scatterer`
        The scatterer parent.
    sampling : :class:`int`
        samplingber of point to evaluate the far-fields in spherical coord.
    distance : :class:`float`
        distance at which we evaluate the far-fields.
    """

    def __init__(self, sampling: int = 200, parent_scatterer=None, distance: float = 1.):
        self.parent_scatterer = parent_scatterer

        self.E_phi, self.E_theta, self.theta, self.phi = parent_scatterer.Bind.get_full_fields(sampling=sampling, r=1)

    def plot(self, show_ource=True, show_axes=True):
        phi_mesh, theta_mesh = numpy.meshgrid(self.phi, self.theta)

        x, y, z = spherical_to_cartesian(r=phi_mesh * 0 + 0.5, phi=phi_mesh, theta=theta_mesh)

        figure = Scene3D(shape=(1, 4), window_size=[2500, 800])

        enumeration = zip(
            [(0, 0), (0, 1), (0, 2), (0, 3)],
            [self.E_phi.real, self.E_phi.imag, self.E_theta.real, self.E_theta.imag],
            ['Phi real', 'Phi imaginary', 'Theta real', 'Theta imaginary']
        )

        for plot_number, field, name in enumeration:
            figure.add_mesh(
                plot_number=plot_number,
                x=x,
                y=y,
                z=z,
                scalars=field.T,
                color_map='seismic',
                scalar_bar_args={'title': f'{name} Amplitude'}
            )

            if 'Phi' in name:
                figure.add_phi_vector_field(plot_number)
            elif 'Theta' in name:
                figure.add_theta_vector_field(plot_number)

            figure.add_unit_axes_to_ax(plot_number=plot_number)
            figure.add_text_to_axes(plot_number=plot_number, text=f'{name} field')

        return figure


class Footprint():
    r"""
    Class representing footprint of the scatterer.
    The footprint usually depend on the scatterer and the detector.
    For more information see references in the
    `documentation <https://pymiesim.readthedocs.io/en/latest>`_
    The footprint is defined as:

    .. math::
        \\text{Footprint} = \\big| \\mathscr{F}^{-1} \\big\\{ \\tilde{ \\psi }\
        (\\xi, \\nu), \\tilde{ \\phi}_{l,m}(\\xi, \\nu)  \\big\\} \
        (\\delta_x, \\delta_y) \\big|^2


    Parameters
    ----------
    scatterer : :class:`Scatterer`
        The scatterer.
    detector : :class:`Detector`
        The detector.
    sampling : :class:`int`
        samplingber of point to evaluate the footprint in cartesian coord.

    """

    def __init__(self, scatterer, detector):
        self.detector = detector
        self.scatterer = scatterer
        self.padding_factor = 10

        self.sampling = 500 if isinstance(detector, PyMieSim.detector.LPmode) else detector.sampling

        self._compute_footprint_()

    def _compute_footprint_(self):
        phi, theta = numpy.mgrid[-self.detector.max_angle:self.detector.max_angle:complex(self.sampling),
                                 0:numpy.pi:complex(self.sampling)]

        max_direct = 1 / (numpy.sin(self.detector.max_angle) * self.scatterer.source.k / (2 * numpy.pi))

        x = y = numpy.linspace(-1, 1, self.sampling) * self.sampling / 2 * max_direct / self.padding_factor

        _, phi, theta = rotate_on_x(phi + numpy.pi / 2, theta, numpy.pi / 2)

        far_field_para, far_field_perp = self.scatterer._FarField(
            phi=phi.flatten() + numpy.pi / 2,
            theta=theta.flatten(),
            r=1.0,
            structured=False
        )

        scalarfield = self.detector.get_structured_scalarfield()[0]

        perp = scalarfield * far_field_perp.reshape(theta.shape)

        para = scalarfield * far_field_para.reshape(theta.shape)

        fourier_para = self.get_fourier_component(para)
        fourier_perp = self.get_fourier_component(perp)

        self.mapping = (fourier_para + fourier_perp)
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

        ax.add_mesh(
            x=self.direct_y * 1e6,
            y=self.direct_x * 1e6,
            scalar=self.mapping,
            colormap='gray'
        )

        return figure


# -
