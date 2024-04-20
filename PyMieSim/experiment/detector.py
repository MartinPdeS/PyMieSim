#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from PyMieSim.experiment.setup import Setup
    from collections.abc import Iterable
    from typing import NoReturn

import numpy
from dataclasses import dataclass, field

from DataVisual import units
from PyMieSim.binary.Sets import CppDetectorSet
from PyMieSim.binary.Fibonacci import FibonacciMesh as CPPFibonacciMesh
from PyMieSim.modes import hermite_gauss, laguerre_gauss, linearly_polarized


@dataclass
class BaseDetector():
    """
    A base class for detectors in Mie scattering simulations, handling the common functionalities
    for different types of detectors. It is designed to format input parameters, compute necessary
    detector properties, and bind the detector to a simulation experiment setup.

    Attributes:
        NA (Iterable[float]): Numerical aperture of the imaging system. It defines the range of angles
                              over which the system can accept or emit light.
        gamma_offset (Iterable[float]): Angle offset of the detector in the direction perpendicular
                                        to polarization, in degrees.
        phi_offset (Iterable[float]): Angle offset of the detector in the direction parallel to
                                      polarization, in degrees.
        polarization_filter (Iterable[float]): Angle of the polarization filter in front of the
                                              detector, in degrees. This specifies the orientation
                                              of the polarization filter.
        sampling (int): Number of sampling points for field evaluation. It dictates the resolution
                        at which the field is sampled by the detector.

    Note:
        This class is not intended to be instantiated directly but serves as a base for specific
        detector types.
    """
    NA: Iterable
    gamma_offset: Iterable
    phi_offset: Iterable
    polarization_filter: Iterable
    sampling: int

    def __post_init__(self) -> NoReturn:
        """
        Initializes and prepares the detector instance by formatting inputs, computing field arrays,
        determining rotation angles, building parameters for visualization, and initializing C++ bindings.
        This method is automatically called after the class has been initialized.

        Returns:
            NoReturn
        """
        self.mapping = {
            'scalarfield': None,
            'NA': None,
            'phi_offset': None,
            'gamma_offset': None,
            'polarization_filter': None,
        }

        self.compute_scalar_fields()

        self.get_rotation_angle_from_mode_number()

        self.initialize_binding()

    def get_rotation_angle_from_mode_number(self) -> NoReturn:
        """
        Computes the rotation angle from the mode number for CoherentMode detectors; for others, sets it to zero.
        This method establishes the orientation of the detector based on its mode number, which is crucial
        for accurate simulation of light detection.

        Returns:
            NoReturn
        """
        if self.name.lower() == 'photodiode':
            self.rotation_angle = numpy.zeros(self.scalar_fields.size)
            return

        rotation_angle_list = []

        self.mode_number = [
            mode if ":" in mode else mode + ":00" for mode in self.mode_number
        ]

        for mode_number in self.mode_number:
            _, rotation_angle = mode_number.split(':')

            rotation_angle_list.append(rotation_angle)

        self.rotation_angle = numpy.asarray(rotation_angle_list).astype(float)

    def bind_to_experiment(self, experiment: Setup) -> NoReturn:
        """
        Binds this detector to a specified experimental setup, integrating it into the simulation workflow.

        Parameters:
            experiment (Setup): The experimental setup to which the detector will be bound.

        Returns:
            NoReturn
        """
        experiment.binding.set_detector(self.binding)

    def get_datavisual_table(self) -> NoReturn:
        """
        Appends the scatterer's properties to a given table for visualization purposes. This enables the
        representation of scatterer properties in graphical formats.

        Parameters:
            table (list): The table to which the scatterer's properties will be appended.

        Returns:
            list: The updated table with the scatterer's properties included.
        """

        self.mapping['scalarfield'] = units.Custom(
            long_label='Field',
            short_label='field',
            base_values=self.scalar_fields,
            value_representation=self.detector_name
        )

        self.mapping['NA'] = units.Index(
            long_label='Numerical aperture',
            short_label='NA',
            base_values=self.NA,
            use_prefix=False,
            string_format=""
        )

        self.mapping['phi_offset'] = units.Degree(
            long_label='Phi angle',
            short_label=r'$\phi_{offset}$',
            base_values=self.phi_offset,
            use_prefix=False,
            string_format='.1f'
        )

        self.mapping['gamma_offset'] = units.Degree(
            long_label='Gamma angle',
            short_label=r'$\phi_{offset}$',
            base_values=self.gamma_offset,
            use_prefix=False,
            string_format='.1f'
        )

        self.mapping['polarization_filter'] = units.Degree(
            long_label=r'Polarization filter',
            short_label=r'f$_{pol}$',
            base_values=self.polarization_filter,
            use_prefix=False,
            string_format='.1f'
        )

        return [v for k, v in self.mapping.items() if v is not None]

    def initialize_binding(self) -> NoReturn:
        """
        Initializes the C++ binding for the detector, configuring it with the necessary parameters for
        simulation. This step is essential for integrating the Python-defined detector with the underlying
        C++ simulation engine.

        Returns:
            NoReturn
        """
        point_coupling = True if self.coupling_mode == 'point' else False

        self.binding_kwargs = dict(
            scalar_field=self.scalar_fields.astype(complex),
            NA=numpy.atleast_1d(self.NA).astype(float),
            phi_offset=numpy.deg2rad(numpy.atleast_1d(self.phi_offset).astype(float)),
            gamma_offset=numpy.deg2rad(numpy.atleast_1d(self.gamma_offset).astype(float)),
            polarization_filter=numpy.deg2rad(numpy.atleast_1d(self.polarization_filter).astype(float)),
            point_coupling=point_coupling,
            coherent=self.coherent,
            rotation_angle=self.rotation_angle.astype(float)
        )

        self.binding = CppDetectorSet(**self.binding_kwargs)

    def interpret_mode_name(self, mode_name: str) -> tuple[int, int, float]:
        """
        Intepret the mode_number parameter to check if there is a rotation associated

        :returns:   The two mode numbers plus rotation angle
        :rtype:     tuple[int, int, float]
        """
        assert isinstance(mode_name, str), "Mode number must be a string. Example 'LP01' or 'HG11:90 (rotation of 90 degree)'"
        if ':' not in mode_name:
            mode_name += ':0'

        mode_family_name = mode_name[:2].lower()

        split = mode_name.split(':')

        mode_name, rotation_angle = split

        number_0, number_1 = mode_name[-2:]

        number_0, number_1 = int(number_0), int(number_1)

        rotation_angle = float(rotation_angle)

        return mode_family_name, number_0, number_1, rotation_angle


@dataclass
class Photodiode(BaseDetector):
    """
    Represents a photodiode detector in Mie scattering simulations, extending the BaseDetector
    with specific attributes and methods for photodiode-type detection.

    Attributes:
        coupling_mode (str): Defines the method for computing mode coupling, either 'point' or 'Mean'.
                             Default is 'point'.
        coherent (bool): Indicates whether the detection scheme is coherent (True) or incoherent (False).
                         This is initialized to False and not meant to be modified.
        name (str): The name of the detector, initialized to "Photodiode" and not intended for modification.

    Methods:
        compute_scalar_fields(): Overrides the BaseDetector method to provide the field arrays specific
                            to a photodiode detector. It generates a scalar field array representing
                            the detection scheme.

    Note:
        This class specifically models a photodiode detector and its interaction within a Mie scattering
        simulation experiment.
    """
    coupling_mode: str = 'point'
    coherent: bool = field(default=False, init=False)
    name: str = field(default="Photodiode", init=False)

    def compute_scalar_fields(self) -> numpy.ndarray:
        """
        Generates a scalar field array representing the photodiode detection scheme. This method overrides
        the BaseDetector's compute_scalar_fields method to provide functionality specific to photodiode detectors.

        Returns:
            numpy.ndarray: An array of scalar fields corresponding to the photodiode detection scheme.
        """
        self.scalar_fields = numpy.ones([1, self.sampling])

        self.detector_name = ['Photodiode']


@dataclass
class CoherentMode(BaseDetector):
    mode_number: Iterable
    """ List of mode to be used. """
    coupling_mode: str = 'point'
    """ Method for computing mode coupling. Either point or Mean. """
    coherent: bool = field(default=True, init=False)
    """ Describe the detection scheme coherent or uncoherent. """
    name: str = field(default="CoherentMode", init=False)
    """ name of the set """

    def compute_scalar_fields(self) -> numpy.ndarray:
        """
        Loads and generates complex scalar field arrays representing the specified CoherentMode modes. This method
        overrides the BaseDetector's compute_scalar_fields to cater specifically to CoherentMode detectors.

        Returns:
            numpy.ndarray: An array of complex scalar fields representing the CoherentMode modes involved in the detection.
        """
        self.mode_number = numpy.atleast_1d(self.mode_number).astype(str)

        proxy_fibonacci_mesh = CPPFibonacciMesh(
            max_angle=0.3,
            phi_offset=0,
            gamma_offset=0,
            sampling=self.sampling,
            rotation_angle=0
        )

        self.scalar_fields = numpy.zeros([len(self.mode_number), self.sampling]).astype(complex)

        self.detector_name = self.mode_number

        for idx, mode_name in enumerate(self.mode_number):
            mode_family_name, number_0, number_1, rotation_angle = self.interpret_mode_name(mode_name)

            match mode_family_name:
                case 'lp':
                    mode_module = linearly_polarized
                case 'hg':
                    mode_module = hermite_gauss
                case 'lg':
                    mode_module = laguerre_gauss
                case _:
                    raise ValueError('Invalid mode family name, it has to be in either: LP, HG or LG')

            self.scalar_fields[idx] = mode_module.interpolate_from_fibonacci_mesh(
                fibonacci_mesh=proxy_fibonacci_mesh,
                number_0=number_0,
                number_1=number_1,
            )

            print(self.scalar_fields[idx])
# -
