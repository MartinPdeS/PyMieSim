#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from PyMieSim.single.scatterer import Sphere, CoreShell, Cylinder


import numpy
import logging
from dataclasses import dataclass, field

from PyMieSim.single.representations import Footprint
from PyMieSim.binary.Fibonacci import FibonacciMesh as CPPFibonacciMesh  # has to be imported as extension  # noqa: F401
from PyMieSim.binary.DetectorInterface import BindedDetector
from PyMieSim.binary import ModeField
from PyMieSim.tools.special_functions import NA_to_angle
from MPSPlots.render3D import SceneList as SceneList3D


class GenericDetector():
    """
    Base class for different types of detectors with methods for setup, rotation,
    calculating light coupling, and generating footprint and 3D plots.
    """

    def __init__(self, **kwargs):
        """
        Initialize detector attributes from keyword arguments.
        """
        for key, value in kwargs.items():
            setattr(self, key, value)

    def coupling(self, scatterer: Sphere | CoreShell | Cylinder) -> float:
        r"""
        Calculate the light coupling between the detector and a scatterer.

        .. math::
            |\iint_{\Omega}  \Phi_{det} \,\, \Psi_{scat}^* \,  d \Omega|^2

        | Where:
        |   :math:`\Phi_{det}` is the capturing field of the detector and
        |   :math:`\Psi_{scat}` is the scattered field.

        Args:
            scatterer (Sphere|CoreShell|Cylinder): The scatterer object.

        Returns:
            float: The coupling in watts.
        """
        return getattr(self.cpp_binding, "Coupling" + type(scatterer).__name__)(scatterer.binding)

    def get_footprint(self, scatterer: Sphere | CoreShell | Cylinder) -> Footprint:
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
            scatterer (Sphere|CoreShell|Cylinder): The scatterer object.

        Returns:
            Footprint: The scatterer footprint with this detector.
        """
        return Footprint(scatterer=scatterer, detector=self)

    def plot(self) -> SceneList3D:
        """
        Plot the real and imaginary parts of the scattered fields.

        Returns:
            SceneList3D: The 3D plotting scene containing the field plots.
        """
        coordinate = numpy.row_stack((
            self.cpp_binding.mesh.x,
            self.cpp_binding.mesh.y,
            self.cpp_binding.mesh.z
        ))

        figure = SceneList3D()

        for scalar_type in ['real', 'imag']:
            scalar = getattr(numpy.asarray(self.cpp_binding.scalar_field), scalar_type)

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
    mean_coupling: bool = field(default=False, init=False)
    """ Indicate if the coupling mechanism is point-wise or mean-wise. Value is either point or mean. """
    rotation: str = field(default=0, init=False)
    """ Indicate the rotation of the field in the axis of propagation. """

    def __post_init__(self):
        self.polarization_filter = numpy.float64(self.polarization_filter)

        self.max_angle = NA_to_angle(NA=self.NA)

        self.cpp_binding = BindedDetector(
            mode_number='NC00',
            sampling=self.sampling,
            NA=self.NA,
            phi_offset=numpy.deg2rad(self.phi_offset),
            gamma_offset=numpy.deg2rad(self.gamma_offset),
            polarization_filter=numpy.deg2rad(self.polarization_filter),
            rotation=numpy.deg2rad(self.rotation),
            coherent=self.coherent,
            mean_coupling=self.mean_coupling
        )

    def get_structured_scalarfield(self, sampling: int = 100) -> numpy.ndarray:
        return numpy.ones([sampling, sampling])


@dataclass
class IntegratingSphere(Photodiode):
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
    mean_coupling: bool = field(default=False, init=False)
    """ indicate if the coupling mechanism is point-wise or mean-wise. Value is either point or mean. """
    rotation: float = field(default=0, init=False)
    """ Indicate the rotation of the field in the axis of propagation. """

    def __post_init__(self):
        super().__post_init__()

    def get_structured_scalarfield(self, sampling: int = 100) -> numpy.ndarray:
        return numpy.ones([sampling, sampling])


@dataclass
class CoherentMode(GenericDetector):
    """
    Detector type class representing a laser Hermite-Gauss mode, light coupling mechanism is coherent
    and thus, dependent of the phase of the impinging scattered light field.
    """
    mode_number: str
    """ String representing the HG mode to be initialized (e.g. 'LP01', 'HG11', 'LG22' etc)"""
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
    mean_coupling: bool = False
    """ indicate if the coupling mechanism is point-wise (if setted True) or mean-wise (if setted False). """
    coherent: bool = field(default=True, init=False)
    """ Indicate if the coupling mechanism is coherent or not. """
    rotation: float = 90
    """ Indicate the rotation of the field in the axis of propagation. """

    def __post_init__(self):
        if self.NA > 0.3 or self.NA < 0:
            logging.warning(f"High values of NA: {self.NA} do not comply with the paraxial approximation. Value under 0.3 are prefered.")

        self.mode_family = self.mode_number[:2]

        if self.mode_family.lower() not in ['lp', 'lg', 'hg']:
            raise ValueError(f'Invalid mode family: {self.mode_family}. Options are ["LP", "LG", "HG"]')

        number_0, number_1 = self.mode_number[2:]
        self.number_0, self.number_1 = int(number_0), int(number_1)

        match self.mode_family.lower():
            case 'lp':
                self.azimuthal_number, self.radial_number = self.number_0, self.number_1
                self.cpp_mode_field_getter = ModeField.get_LP
            case 'lg':
                self.azimuthal_number, self.radial_number = self.number_0, self.number_1
                self.cpp_mode_field_getter = ModeField.get_LG
            case 'hg':
                self.x_number, self.y_number = self.number_0, self.number_1
                self.cpp_mode_field_getter = ModeField.get_HG

        self.max_angle = NA_to_angle(NA=self.NA)

        self.polarization_filter = numpy.float64(self.polarization_filter)

        self.cpp_binding = BindedDetector(
            mode_number=self.mode_number,
            sampling=self.sampling,
            NA=self.NA,
            phi_offset=numpy.deg2rad(self.phi_offset),
            gamma_offset=numpy.deg2rad(self.gamma_offset),
            polarization_filter=numpy.deg2rad(self.polarization_filter),
            rotation=numpy.deg2rad(self.rotation),
            coherent=self.coherent,
            mean_coupling=self.mean_coupling
        )

    def get_structured_scalarfield(self, sampling: int = 100) -> numpy.ndarray:
        x_mesh, y_mesh = numpy.mgrid[-100:100:complex(sampling), -100:100:complex(sampling)]

        coordinates = numpy.row_stack((
            x_mesh.ravel(),
            y_mesh.ravel(),
        ))

        norm = numpy.sqrt(numpy.square(coordinates).sum(axis=0)).max()

        coordinates /= norm

        field = self.cpp_mode_field_getter(
            coordinates[0],
            coordinates[1],
            self.number_0,
            self.number_1
        )

        return field.reshape([sampling, sampling])


# -
