#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import logging
from dataclasses import dataclass

from PyMieSim.representations import Footprint
from PyMieSim.mesh import FibonacciMesh
from PyMieSim.binary.DetectorInterface import BindedDetector
from PyMieSim import load_lp_mode
from PyMieSim.tools.special_functions import NA_to_angle

from MPSPlots.render3D import SceneList as SceneList3D


@dataclass
class GenericDetector():
    r"""
    .. note::
        Detector type class representing a photodiode, light coupling is
        thus independant of the phase of the latter.
    """
    scalar_field: numpy.ndarray
    """ Array representing the detection field distribution. """
    NA: float
    """ Numerical aperture of imaging system. """
    gamma_offset: float
    """ Angle [Degree] offset of detector in the direction perpendicular to polarization. """
    phi_offset: float
    """ Angle [Degree] offset of detector in the direction parallel to polarization. """
    polarization_filter: float
    """ Angle [Degree] of polarization filter in front of detector. """
    coupling_mode: str = 'Point'
    """ Method for computing mode coupling. Either Point or Mean. """
    coherent: bool = False
    """ Describe the detection scheme coherent or uncoherent. """

    def __post_init__(self):
        self.scalar_field = self.scalar_field.astype(complex)
        self.sampling = self.scalar_field.size
        self.max_angle = NA_to_angle(self.NA)
        self.polarization_filter = numpy.float64(self.polarization_filter)

        self.Mesh = FibonacciMesh(
            max_angle=self.max_angle,
            sampling=self.sampling,
            phi_offset=self.phi_offset,
            gamma_offset=self.gamma_offset
        )

        self._get_binding_()

    def _get_binding_(self):
        point_coupling = True if self.coupling_mode.lower() == 'point' else False

        self.cpp_binding = BindedDetector(
            scalar_field=self.scalar_field,
            NA=self.NA,
            phi_offset=numpy.deg2rad(self.phi_offset),
            gamma_offset=numpy.deg2rad(self.gamma_offset),
            polarization_filter=numpy.deg2rad(self.polarization_filter),
            coherent=self.coherent,
            point_coupling=point_coupling
        )

    def get_structured_scalarfield(self):
        return numpy.ones([self.sampling, self.sampling])

    def coupling(self, scatterer):
        r"""
        .. note::
            Return the value of the scattererd light coupling as computed as:

            .. math::
                |\iint_{\Omega}  \Phi_{det} \,\, \Psi_{scat}^* \,  d \Omega|^2

            | Where:
            |   :math:`\Phi_{det}` is the capturing field of the detector and
            |   :math:`\Psi_{scat}` is the scattered field.

        Parameters
        ----------
        Scatterer : :class:`Scatterer`
            Scatterer instance (sphere, cylinder, ...).

        Returns
        -------
        :class:`float`
            Value of the coupling.

        """

        return getattr(self.cpp_binding, "Coupling" + type(scatterer).__name__)(scatterer.Bind)

    def get_footprint(self, scatterer) -> Footprint:
        r"""
        .. note::
            Return the footprint of the scattererd light coupling with the
            detector as computed as:

            .. math::
                \big| \mathscr{F}^{-1} \big\{ \tilde{ \psi } (\xi, \nu),\
                       \tilde{ \phi}_{l,m}(\xi, \nu)  \big\}
                       (\delta_x, \delta_y) \big|^2

            | Where:
            |   :math:`\Phi_{det}` is the capturing field of the detector and
            |   :math:`\Psi_{scat}` is the scattered field.

        Parameters
        ----------
        Scatterer : :class:`Scatterer`.
            Scatterer instance (sphere, cylinder, ...).

        Returns
        -------
        :class:`Footprint`.
            Dictionnary subclass with all pertienent information.

        """
        return Footprint(scatterer=scatterer, detector=self)

    def plot(self) -> SceneList3D:
        r"""
        .. note::
            Method that plot the real part of the scattered field
            (:math:`E_{\theta}` and :math:`E_{\phi}`).

        """
        coordinate = numpy.c_[self.Mesh.X, self.Mesh.Y, self.Mesh.Z].T

        figure = SceneList3D()

        for scalar_type in ['real', 'imag']:
            scalar = getattr(self.scalar_field, scalar_type)

            ax = figure.append_ax()
            artist = ax.add_unstructured_mesh(
                coordinates=coordinate,
                scalar_coloring=scalar,
                symmetric_map=True,
                symmetric_colormap=True
            )

            ax.add_unit_sphere()
            ax.add_unit_axis(show_label=False)
            ax.add_colorbar(artist=artist, title=f'field [{scalar_type}]')

        return figure


class Photodiode(GenericDetector):
    def __init__(self,
            NA: float,
            sampling: int,
            gamma_offset: float,
            phi_offset: float,
            polarization_filter: float = None):

        scalar_field = numpy.ones(sampling)

        super().__init__(
            scalar_field=scalar_field,
            NA=NA,
            phi_offset=phi_offset,
            gamma_offset=gamma_offset,
            polarization_filter=polarization_filter,
            coherent=False,
            coupling_mode='Point'
        )


class IntegratingSphere(GenericDetector):
    def __init__(self, sampling: int, polarization_filter: float = None):
        scalar_field = numpy.ones(sampling)

        super().__init__(
            scalar_field=scalar_field,
            NA=2,
            phi_offset=0,
            gamma_offset=0,
            polarization_filter=polarization_filter,
            coherent=False,
            coupling_mode='Point'
        )


class LPmode(GenericDetector):
    def __init__(self,
            mode_number: str,
            NA: float,
            gamma_offset: float,
            phi_offset: float,
            sampling: int = 200,
            rotation: float = 0,
            polarization_filter: float = None,
            coupling_mode: str = 'Point'):

        if NA > 0.3 or NA < 0:
            logging.warning("High values of NA do not comply with paraxial approximation. Value under 0.3 are prefered.")

        self.mode_number = mode_number

        scalar_field = load_lp_mode(
            mode_number=self.mode_number,
            structure_type='unstructured',
            sampling=sampling
        )

        super().__init__(
            scalar_field=scalar_field,
            NA=NA,
            phi_offset=phi_offset,
            gamma_offset=gamma_offset,
            polarization_filter=polarization_filter,
            coherent=True,
            coupling_mode=coupling_mode
        )

    def get_structured_scalarfield(self):
        return load_lp_mode(
            mode_number=self.mode_number,
            structure_type='structured',
        ),


# -
