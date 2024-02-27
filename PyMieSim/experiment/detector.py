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

from DataVisual import Xparameter
from PyMieSim import load_lp_mode
import PyMieSim.datavisual_x_parameters as Kwargs
from PyMieSim.binary.Sets import CppDetectorSet


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

    parameter_str_list = [
        'NA',
        'phi_offset',
        'gamma_offset',
        'polarization_filter',
    ]

    def __post_init__(self) -> NoReturn:
        """
        Initializes and prepares the detector instance by formatting inputs, computing field arrays,
        determining rotation angles, building parameters for visualization, and initializing C++ bindings.
        This method is automatically called after the class has been initialized.

        Returns:
            NoReturn
        """
        self.format_inputs()

        self.get_fields_array()

        self.get_rotation_angle_from_mode_number()

        self.build_x_parameters()

        self.initialize_binding()

    def get_rotation_angle_from_mode_number(self) -> NoReturn:
        """
        Computes the rotation angle from the mode number for LPMode detectors; for others, sets it to zero.
        This method establishes the orientation of the detector based on its mode number, which is crucial
        for accurate simulation of light detection.

        Returns:
            NoReturn
        """
        if self.name.lower() == 'photodiode':
            self.rotation_angle = numpy.zeros(self.scalarfield.values.size)
            return

        rotation_angle_list = []

        self.mode_number = [
            mode if ":" in mode else mode + ":00" for mode in self.mode_number
        ]

        for mode_number in self.mode_number:
            _, rotation_angle = mode_number.split(':')

            rotation_angle_list.append(rotation_angle)

        self.rotation_angle = numpy.asarray(rotation_angle_list).astype(float)

    def format_inputs(self) -> NoReturn:
        """
        Formats the input attributes (NA, phi_offset, gamma_offset, polarization_filter) into numpy arrays
        for further processing. This ensures that all input parameters are correctly structured for the
        simulation and C++ binding processes.

        Returns:
            NoReturn
        """
        self.NA = numpy.atleast_1d(self.NA).astype(float)
        self.phi_offset = numpy.atleast_1d(self.phi_offset).astype(float)
        self.gamma_offset = numpy.atleast_1d(self.gamma_offset).astype(float)
        self.polarization_filter = numpy.atleast_1d(self.polarization_filter).astype(float)

    def bind_to_experiment(self, experiment: Setup) -> NoReturn:
        """
        Binds this detector to a specified experimental setup, integrating it into the simulation workflow.

        Parameters:
            experiment (Setup): The experimental setup to which the detector will be bound.

        Returns:
            NoReturn
        """
        experiment.binding.set_detector(self.binding)

    def build_x_parameters(self) -> NoReturn:
        """
        Constructs the Xparameters for visualization, translating the detector's attributes into a format
        suitable for inclusion in data visualization tools. This facilitates the graphical representation
        of simulation results.

        Returns:
            NoReturn
        """
        for parameter_str in self.parameter_str_list:
            parameter = getattr(self, parameter_str)

            parameter = numpy.array(parameter)

            kwargs_parameter = getattr(Kwargs, parameter_str)

            x_parameter = Xparameter(values=parameter, **kwargs_parameter)

            setattr(self, parameter_str, x_parameter)

    def append_to_table(self, table: list) -> list:
        """
        Appends the detector's attributes to a given visualization table.

        Parameters:
            table (list): The table to which the detector's attributes will be appended.

        Returns:
            list: The updated table with the detector's attributes included.
        """
        return [*table, self.scalarfield, self.NA, self.phi_offset, self.gamma_offset, self.polarization_filter]

    def initialize_binding(self) -> NoReturn:
        """
        Initializes the C++ binding for the detector, configuring it with the necessary parameters for
        simulation. This step is essential for integrating the Python-defined detector with the underlying
        C++ simulation engine.

        Returns:
            NoReturn
        """
        point_coupling = True if self.coupling_mode == 'point' else False

        phi_offset_rad = numpy.deg2rad(self.phi_offset.values)

        gamma_offset_rad = numpy.deg2rad(self.gamma_offset.values)

        polarization_filter_rad = numpy.deg2rad(self.polarization_filter.values)

        self.binding = CppDetectorSet(
            scalarfield=self.scalarfield.values.astype(complex),
            NA=self.NA.values,
            phi_offset=phi_offset_rad,
            gamma_offset=gamma_offset_rad,
            polarization_filter=polarization_filter_rad,
            point_coupling=point_coupling,
            coherent=self.coherent,
            rotation_angle=self.rotation_angle.astype(float)
        )


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
        get_fields_array(): Overrides the BaseDetector method to provide the field arrays specific
                            to a photodiode detector. It generates a scalar field array representing
                            the detection scheme.

    Note:
        This class specifically models a photodiode detector and its interaction within a Mie scattering
        simulation experiment.
    """
    coupling_mode: str = 'point'
    coherent: bool = field(default=False, init=False)
    name: str = field(default="Photodiode", init=False)

    def get_fields_array(self) -> numpy.ndarray:
        """
        Generates a scalar field array representing the photodiode detection scheme. This method overrides
        the BaseDetector's get_fields_array method to provide functionality specific to photodiode detectors.

        Returns:
            numpy.ndarray: An array of scalar fields corresponding to the photodiode detection scheme.
        """
        scalarfield = numpy.ones([1, self.sampling])

        self.scalarfield = Xparameter(
            values=scalarfield,
            representation=['Photodiode'],
            **Kwargs.scalarfield
        )


@dataclass
class LPMode(BaseDetector):
    mode_number: Iterable
    """ List of mode to be used. """
    coupling_mode: str = 'point'
    """ Method for computing mode coupling. Either point or Mean. """
    coherent: bool = field(default=True, init=False)
    """ Describe the detection scheme coherent or uncoherent. """
    name: str = field(default="LPMode", init=False)
    """ name of the set """

    def get_fields_array(self) -> numpy.ndarray:
        """
        Loads and generates complex scalar field arrays representing the specified LP modes. This method
        overrides the BaseDetector's get_fields_array to cater specifically to LPMode detectors.

        Returns:
            numpy.ndarray: An array of complex scalar fields representing the LP modes involved in the detection.
        """
        self.mode_number = numpy.atleast_1d(self.mode_number).astype(str)

        scalarfield = load_lp_mode(
            mode_number=self.mode_number,
            sampling=self.sampling,
            structure_type='unstructured'
        ).astype(complex)

        self.scalarfield = Xparameter(
            values=scalarfield,
            representation=self.mode_number,
            **Kwargs.scalarfield
        )
