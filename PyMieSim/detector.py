#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from PyMieSim.scatterer import Sphere, CoreShell, Cylinder


import numpy
import logging
from dataclasses import dataclass, field

from PyMieSim.representations import Footprint
from PyMieSim.mesh import FibonacciMesh
from PyMieSim.binary.DetectorInterface import BindedDetector
from PyMieSim import load_lp_mode
from PyMieSim.tools.special_functions import NA_to_angle

from MPSPlots.render3D import SceneList as SceneList3D


class GenericDetector():

    def initialize(self):
        self.max_angle = NA_to_angle(NA=self.NA)
        self.polarization_filter = numpy.float64(self.polarization_filter)

        self.fibonacci_mesh = FibonacciMesh(
            max_angle=self.max_angle,
            sampling=self.sampling,
            phi_offset=self.phi_offset,
            gamma_offset=self.gamma_offset
        )

        self.set_cpp_binding()

    def set_cpp_binding(self) -> None:
        point_coupling = True if self.coupling_mode.lower() == 'point' else False

        self.cpp_binding = BindedDetector(
            scalar_field=self.unstructured_farfield,
            NA=self.NA,
            phi_offset=numpy.deg2rad(self.phi_offset),
            gamma_offset=numpy.deg2rad(self.gamma_offset),
            polarization_filter=numpy.deg2rad(self.polarization_filter),
            coherent=self.coherent,
            point_coupling=point_coupling
        )

    def coupling(self, scatterer: Sphere | CoreShell | Cylinder) -> float:
        r"""
        Return the value of the scattererd light coupling as computed as:

        .. math::
            |\iint_{\Omega}  \Phi_{det} \,\, \Psi_{scat}^* \,  d \Omega|^2

        | Where:
        |   :math:`\Phi_{det}` is the capturing field of the detector and
        |   :math:`\Psi_{scat}` is the scattered field.

        :param      scatterer:  The scatterer
        :type       scatterer:  Sphere | CoreShell | Cylinder

        :returns:   The coupling in watt
        :rtype:     float
        """
        return getattr(self.cpp_binding, "Coupling" + type(scatterer).__name__)(scatterer.Bind)

    def get_footprint(self, scatterer: Sphere | CoreShell | Cylinder) -> Footprint:
        r"""
        Return the footprint of the scattererd light coupling with the
        detector as computed as:

        .. math::
            \big| \mathscr{F}^{-1} \big\{ \tilde{ \psi } (\xi, \nu),\
                   \tilde{ \phi}_{l,m}(\xi, \nu)  \big\}
                   (\delta_x, \delta_y) \big|^2

        | Where:
        |   :math:`\Phi_{det}` is the capturing field of the detector and
        |   :math:`\Psi_{scat}` is the scattered field.

        :param      detector:  The detector
        :type       detector:  Sphere | CoreShell | Cylinder

        :returns:   The scatterer footprint.
        :rtype:     Footprint
        """
        return Footprint(scatterer=scatterer, detector=self)

    def plot(self) -> SceneList3D:
        r"""
        Method that plot the real part of the scattered fields (:math:`E_{\theta}` and :math:`E_{\phi}`).

        :returns:   3D plotting scene
        :rtype:     SceneList3D
        """
        coordinate = numpy.c_[self.fibonacci_mesh.X, self.fibonacci_mesh.Y, self.fibonacci_mesh.Z].T

        figure = SceneList3D()

        for scalar_type in ['real', 'imag']:
            scalar = getattr(self.unstructured_farfield, scalar_type)

            ax = figure.append_ax()
            artist = ax.add_unstructured_mesh(
                coordinates=coordinate,
                scalar_coloring=scalar,
                symmetric_map=True,
                symmetric_colormap=True
            )

            ax.add_unit_sphere(opacity=0.3)
            ax.add_unit_axis(show_label=False)
            ax.add_colorbar(artist=artist, title=f'field [{scalar_type}]')

        return figure


@dataclass
class Photodiode(GenericDetector):
    """
    Detector type class representing a photodiode, light coupling mechanism is non-coherent and thus
    independent of the phase of the impinging scattered light field.
    """
    NA: float
    """ Numerical aperture of imaging system. """
    gamma_offset: float
    """ Angle [Degree] offset of detector in the direction perpendicular to polarization. """
    phi_offset: float
    """ Angle [Degree] offset of detector in the direction parallel to polarization. """
    sampling: int = 200
    """ Sampling of the farfield distribution """
    polarization_filter: float = None
    """ Angle [Degree] of polarization filter in front of detector. """
    coherent: bool = field(default=False, init=False)
    """ Indicate if the coupling mechanism is coherent or not """
    coupling_mode: str = field(default='point', init=False)
    """ indicate if the coupling mechanism is point-wise or mean-wise. Value is either point or mean. """

    def __post_init__(self):
        self.unstructured_farfield = numpy.ones(self.sampling, dtype=complex)

        super().initialize()

    def get_structured_scalarfield(self, sampling: int = 500) -> numpy.ndarray:
        return numpy.ones([sampling, sampling])


@dataclass
class IntegratingSphere(GenericDetector):
    """
    Detector type class representing a photodiode, light coupling mechanism is non-coherent and thus
    independent of the phase of the impinging scattered light field.
    """
    sampling: int = 200
    """ Sampling of the farfield distribution """
    polarization_filter: float = None
    """ Angle [Degree] of polarization filter in front of detector. """
    NA: float = field(default=2, init=False)
    """ Numerical aperture of imaging system. """
    gamma_offset: float = field(default=0, init=False)
    """ Angle [Degree] offset of detector in the direction perpendicular to polarization. """
    phi_offset: float = field(default=0, init=False)
    """ Angle [Degree] offset of detector in the direction parallel to polarization. """
    coherent: bool = field(default=False, init=False)
    """ Indicate if the coupling mechanism is coherent or not """
    coupling_mode: str = field(default='point', init=False)
    """ indicate if the coupling mechanism is point-wise or mean-wise. Value is either point or mean. """

    def __post_init__(self):
        self.unstructured_farfield = numpy.ones(self.sampling, dtype=complex)

        super().initialize()

    def get_structured_scalarfield(self, sampling: int = 500) -> numpy.ndarray:
        return numpy.ones([sampling, sampling])


@dataclass
class LPmode(GenericDetector):
    """
    Detector type class representing a fiber LP mode, light coupling mechanism is coherent
    and thus, dependent of the phase of the impinging scattered light field.
    """
    mode_number: str
    """ String representing the LP mode to be initialized (e.g. 'LP01', 'LP11' etc)"""
    NA: float
    """ Numerical aperture of imaging system. """
    gamma_offset: float
    """ Angle [Degree] offset of detector in the direction perpendicular to polarization. """
    phi_offset: float
    """ Angle [Degree] offset of detector in the direction parallel to polarization. """
    sampling: int = 200
    """ Sampling of the farfield distribution """
    rotation: float = 0
    """ Rotation along the NA axis of the farfield distribution [degree]"""
    polarization_filter: float = None
    """ Angle [Degree] of polarization filter in front of detector. """
    coupling_mode: str = 'Point'
    """ indicate if the coupling mechanism is point-wise or mean-wise. Value is either point or mean. """
    coherent: bool = field(default=True, init=False)
    """ Indicate if the coupling mechanism is coherent or not. """

    def __post_init__(self):
        if self.NA > 0.3 or self.NA < 0:
            logging.warning(f"High values of NA: {self.NA} do not comply with the paraxial approximation. Value under 0.3 are prefered.")

        self.unstructured_farfield = load_lp_mode(
            mode_number=self.mode_number,
            structure_type='unstructured',
            sampling=self.sampling
        )

        super().initialize()

    def get_structured_scalarfield(self, sampling: int = 500) -> numpy.ndarray:
        far_field = load_lp_mode(
            mode_number=self.mode_number,
            structure_type='structured',
            sampling=sampling,
        )

        return far_field


# -
