#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import logging
from dataclasses import field
from pydantic.dataclasses import dataclass
from typing import Union, Optional

from PyMieSim.single.representations import Footprint
from PyMieSim.binary.DetectorInterface import BindedDetector
from PyMieSim.binary import ModeField
from PyMieSim.binary.Fibonacci import FibonacciMesh as CPPFibonacciMesh  # has to be imported as extension  # noqa: F401
from PyMieSim.special_functions import NA_to_angle
from MPSPlots.render3D import SceneList as SceneList3D

from PyMieSim.single import scatterer

c = 299792458.0  #: Speed of light in vacuum (m/s).
epsilon0 = 8.854187817620389e-12  #: Vacuum permittivity (F/m).


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

    def coupling(self, scatterer: Union[scatterer.Sphere, scatterer.CoreShell, scatterer.Cylinder]) -> float:
        r"""
        Calculate the light coupling between the detector and a scatterer.

        .. math::
            |\iint_{\Omega}  \Phi_{det} \,\, \Psi_{scat}^* \,  d \Omega|^2

        | Where:
        |   :math:`\Phi_{det}` is the capturing field of the detector and
        |   :math:`\Psi_{scat}` is the scattered field.

        Args:
            scatterer (Union[scatterer.Sphere, scatterer.CoreShell, scatterer.Cylinder]): The scatterer object.

        Returns:
            float: The coupling in watts.
        """
        return getattr(self.binding, "Coupling" + type(scatterer).__name__)(scatterer.binding)

    def get_footprint(self, scatterer: Union[scatterer.Sphere, scatterer.CoreShell, scatterer.Cylinder]) -> Footprint:
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
            scatterer (Union[scatterer.Sphere, scatterer.CoreShell, scatterer.Cylinder]): The scatterer object.

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
            self.binding.mesh.x,
            self.binding.mesh.y,
            self.binding.mesh.z
        ))

        figure = SceneList3D()

        for scalar_type in ['real', 'imag']:
            scalar = getattr(numpy.asarray(self.binding.scalar_field), scalar_type)

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

    def get_poynting_vector(self, scatterer: Union[scatterer.Sphere, scatterer.CoreShell, scatterer.Cylinder]) -> float:
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

    def get_energy_flow(self, scatterer: Union[scatterer.Sphere, scatterer.CoreShell, scatterer.Cylinder]) -> float:
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


@dataclass(kw_only=True, slots=True, config=dict(extra='forbid'))
class Photodiode(GenericDetector):
    """
    Detector class representing a photodiode with a non-coherent light coupling mechanism.
    This means it is independent of the phase of the impinging scattered light field.

    Attributes:
        NA (float): Numerical aperture of the imaging system.
        gamma_offset (float): Angle [Degree] offset of the detector in the direction perpendicular to polarization.
        phi_offset (float): Angle [Degree] offset of the detector in the direction parallel to polarization.
        sampling (int): Sampling rate of the far-field distribution. Default is 200.
        polarization_filter (Union[float, None]): Angle [Degree] of the polarization filter in front of the detector.
        coherent (bool): Indicates if the coupling mechanism is coherent. Default is False.
        mean_coupling (bool): Indicates if the coupling mechanism is point-wise or mean-wise. Default is False.
        rotation (float): Rotation angle of the field along the axis of propagation. Default is 0.
    """

    NA: float
    gamma_offset: float
    phi_offset: float
    sampling: int = 200
    polarization_filter: Union[float, None] = None
    coherent: bool = field(default=False, init=False)
    mean_coupling: bool = field(default=False, init=False)
    rotation: float = field(default=0, init=False)

    def __post_init__(self):
        self.polarization_filter = numpy.float64(self.polarization_filter)
        self.max_angle = NA_to_angle(NA=self.NA)
        self.binding = BindedDetector(
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

    def get_structured_scalarfield(self, sampling: Optional[int] = 100) -> numpy.ndarray:
        """
        Generate a structured scalar field as a numpy array.

        Args:
            sampling (int): The sampling rate for the scalar field. Default is 100.

        Returns:
            numpy.ndarray: A 2D array representing the structured scalar field.
        """
        return numpy.ones([sampling, sampling])


@dataclass(kw_only=True, slots=True, config=dict(extra='forbid'))
class IntegratingSphere(Photodiode):
    """
    Detector class representing a photodiode with a non-coherent light coupling mechanism.
    This implies independence from the phase of the impinging scattered light field.

    Attributes:
        sampling (int): Sampling rate of the far-field distribution. Default is 200.
        polarization_filter (Union[float, None]): Angle [Degree] of the polarization filter in front of the detector.
        NA (float): Numerical aperture of the imaging system. Default is 2.
        gamma_offset (float): Angle [Degree] offset of the detector in the direction perpendicular to polarization. Default is 0.
        phi_offset (float): Angle [Degree] offset of the detector in the direction parallel to polarization. Default is 0.
        coherent (bool): Indicates if the coupling mechanism is coherent. Default is False.
        mean_coupling (bool): Indicates if the coupling mechanism is point-wise or mean-wise. Default is False.
        rotation (float): Rotation angle of the field along the axis of propagation. Default is 0.
    """

    sampling: int = 200
    polarization_filter: Union[float, None] = None
    NA: float = field(default=2, init=False)
    gamma_offset: float = field(default=0, init=False)
    phi_offset: float = field(default=0, init=False)
    coherent: bool = field(default=False, init=False)
    mean_coupling: bool = field(default=False, init=False)
    rotation: float = field(default=0, init=False)

    def __post_init__(self):
        super(IntegratingSphere, self).__post_init__()

    def get_structured_scalarfield(self, sampling: Optional[int] = 100) -> numpy.ndarray:
        """
        Generate a structured scalar field as a numpy array.

        Args:
            sampling (int): The sampling rate for the scalar field. Default is 100.

        Returns:
            numpy.ndarray: A 2D array representing the structured scalar field.
        """
        return numpy.ones([sampling, sampling])


@dataclass(kw_only=True, slots=True, config=dict(extra='forbid'))
class CoherentMode(GenericDetector):
    """
    Detector class representing a laser Hermite-Gauss mode with a coherent light coupling mechanism.
    This means it depends on the phase of the impinging scattered light field.

    Attributes:
        mode_number (str): String representing the HG mode to be initialized (e.g., 'LP01', 'HG11', 'LG22').
        NA (float): Numerical aperture of the imaging system.
        gamma_offset (float): Angle [Degree] offset of the detector in the direction perpendicular to polarization.
        phi_offset (float): Angle [Degree] offset of the detector in the direction parallel to polarization.
        sampling (int): Sampling rate of the far-field distribution. Default is 200.
        polarization_filter (Union[float, None]): Angle [Degree] of the polarization filter in front of the detector.
        mean_coupling (bool): Indicates if the coupling mechanism is point-wise (True) or mean-wise (False). Default is False.
        coherent (bool): Indicates if the coupling mechanism is coherent. Default is True.
        rotation (float): Rotation angle of the field along the axis of propagation. Default is 90.
    """

    mode_number: str
    NA: float
    gamma_offset: float
    phi_offset: float
    sampling: int = 200
    polarization_filter: Union[float, None] = None
    mean_coupling: bool = False
    coherent: bool = field(default=True, init=False)
    rotation: float = 90

    def __post_init__(self):
        if self.NA > 0.3 or self.NA < 0:
            logging.warning(f"High values of NA: {self.NA} do not comply with the paraxial approximation. Values under 0.3 are preferred.")

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

        self.binding = BindedDetector(
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

    def get_structured_scalarfield(self, sampling: Optional[int] = 100) -> numpy.ndarray:
        """
        Generate a structured scalar field as a numpy array.

        Args:
            sampling (int): The sampling rate for the scalar field. Default is 100.

        Returns:
            numpy.ndarray: A 2D array representing the structured scalar field.
        """
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
