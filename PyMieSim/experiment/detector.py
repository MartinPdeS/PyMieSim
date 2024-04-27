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
            'mode_number': None,
            'sampling': None,
            'rotation': None,
            'NA': None,
            'phi_offset': None,
            'gamma_offset': None,
            'polarization_filter': None,
        }

        self.initialize_binding()

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
            mode_number=numpy.atleast_1d(self.mode_number).astype(str),
            sampling=numpy.atleast_1d(self.sampling).astype(int),
            NA=numpy.atleast_1d(self.NA).astype(float),
            phi_offset=numpy.deg2rad(numpy.atleast_1d(self.phi_offset).astype(float)),
            gamma_offset=numpy.deg2rad(numpy.atleast_1d(self.gamma_offset).astype(float)),
            polarization_filter=numpy.deg2rad(numpy.atleast_1d(self.polarization_filter).astype(float)),
            rotation=numpy.deg2rad(numpy.atleast_1d(self.rotation)).astype(float),
            point_coupling=point_coupling,
            coherent=self.coherent
        )

        self.binding = CppDetectorSet(**self.binding_kwargs)

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

        self.mapping['mode_number'] = units.Custom(
            long_label='Mode number',
            short_label='mode',
            base_values=self.mode_number,
        )

        self.mapping['sampling'] = units.Custom(
            long_label='Sampling',
            short_label='sampling',
            base_values=self.sampling,
        )

        self.mapping['rotation'] = units.Degree(
            long_label='Rotation angle',
            short_label='rot',
            base_values=self.rotation,
            string_format='.1f'
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
            short_label=r'$\phi$',
            base_values=self.phi_offset,
            use_prefix=False,
            string_format='.1f'
        )

        self.mapping['gamma_offset'] = units.Degree(
            long_label='Gamma angle',
            short_label=r'$\gamma$',
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
    rotation: float = 0
    coupling_mode: str = 'point'
    coherent: bool = field(default=False, init=False)
    mode_number: str = field(default='NC00', init=False)


@dataclass
class CoherentMode(BaseDetector):
    mode_number: Iterable
    """ List of mode to be used. """
    coupling_mode: str = 'point'
    """ Method for computing mode coupling. Either point or Mean. """
    rotation: float = 0
    """ Rotation of the coherent mode field. (degree)"""
    coherent: bool = field(default=True, init=False)
    """ Describe the detection scheme coherent or uncoherent. """

    def __post_init__(self) -> numpy.ndarray:
        """
        Loads and generates complex scalar field arrays representing the specified CoherentMode modes. This method
        overrides the BaseDetector's compute_scalar_fields to cater specifically to CoherentMode detectors.

        Returns:
            numpy.ndarray: An array of complex scalar fields representing the CoherentMode modes involved in the detection.
        """
        self.mode_number = numpy.atleast_1d(self.mode_number).astype(str)

        for idx, mode_name in enumerate(self.mode_number):
            mode_family_name = mode_name[0: 2]

            if mode_family_name not in ['LP', 'HG', 'LG', 'NC']:
                raise ValueError('Invalid mode family name, it has to be in either: LP, HG or LG')

        super().__post_init__()

# -
