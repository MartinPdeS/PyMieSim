#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from PyMieSim import single
    from PyMieSim.single.scatterer import Sphere, CoreShell, Cylinder


import numpy
import logging
from dataclasses import dataclass, field
from PyMieSim import load_single_lp_mode

from PyMieSim.single.representations import Footprint
from PyMieSim.binary.Fibonacci import FIBONACCI  # has to be imported as extension
from PyMieSim.binary.DetectorInterface import BindedDetector
from PyMieSim.tools.special_functions import NA_to_angle

from MPSPlots.render3D import SceneList as SceneList3D


class GenericDetector():
    """
    Base class for different types of detectors with methods for setup, rotation,
    calculating light coupling, and generating footprint and 3D plots.
    """

    def __init__(self, **kwargs):
        self.initialize_attributes(**kwargs)
        self.initialize()

    def initialize_attributes(self, **kwargs):
        """
        Initialize detector attributes from keyword arguments.
        """
        for key, value in kwargs.items():
            setattr(self, key, value)

    def initialize(self):
        """
        Perform initial setup by setting max angle based on NA, converting
        polarization_filter to float, and binding C++ backend.
        """
        self.max_angle = NA_to_angle(NA=self.NA)
        self.polarization_filter = numpy.float64(self.polarization_filter)
        self.set_cpp_binding()

    def set_cpp_binding(self):
        """
        Establish binding with the C++ backend for the detector.
        """
        point_coupling = (self.coupling_mode.lower() == 'point')

        self.cpp_binding = BindedDetector(
            scalar_field=self.unstructured_farfield,
            NA=self.NA,
            phi_offset=numpy.deg2rad(self.phi_offset),
            gamma_offset=numpy.deg2rad(self.gamma_offset),
            polarization_filter=numpy.deg2rad(self.polarization_filter),
            coherent=self.coherent,
            rotation_angle=0,
            point_coupling=point_coupling
        )

    def rotate_around_axis(self, rotation_angle: float):
        """
        Rotate the detector mesh around its principal axis by a given angle.

        Args:
            rotation_angle (float): Rotation angle in degrees.
        """
        self.cpp_binding.mesh.rotate_around_axis(rotation_angle)

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
        return getattr(self.cpp_binding, "Coupling" + type(scatterer).__name__)(scatterer.Bind)

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
class LPMode(GenericDetector):
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
    polarization_filter: float = None
    """ Angle [Degree] of polarization filter in front of detector. """
    coupling_mode: str = 'Point'
    """ indicate if the coupling mechanism is point-wise or mean-wise. Value is either point or mean. """
    coherent: bool = field(default=True, init=False)
    """ Indicate if the coupling mechanism is coherent or not. """

    def __post_init__(self):
        if self.NA > 0.3 or self.NA < 0:
            logging.warning(f"High values of NA: {self.NA} do not comply with the paraxial approximation. Value under 0.3 are prefered.")

        self.interpret_mode_name()

        self.unstructured_farfield = load_single_lp_mode(
            mode_number=self.mode_number,
            structure_type='unstructured',
            sampling=self.sampling
        )

        super().initialize()

        if self.rotation_angle != 0:
            self.rotate_around_axis(self.rotation_angle)

    def interpret_mode_name(self) -> None:
        """
        Intepret the mode_number parameter to check if there is a rotation associated

        :returns:   { description_of_the_return_value }
        :rtype:     None
        """
        assert isinstance(self.mode_number, str), "Mode number must be a string. Example 'LP01' or 'LP11:90 (rotation of 90 degree)'"
        if ':' not in self.mode_number:
            self.mode_number += ':0'

        split = self.mode_number.split(':')

        self.mode_number, rotation_angle = split

        self.rotation_angle = float(rotation_angle)

    def get_structured_scalarfield(self, sampling: int = 500) -> numpy.ndarray:
        far_field = load_single_lp_mode(
            mode_number=self.mode_number,
            structure_type='structured',
            sampling=sampling,
        )

        return far_field


# -
