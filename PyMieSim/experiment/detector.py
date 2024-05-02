#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from PyMieSim.experiment.setup import Setup
    from typing import NoReturn, Iterable

import numpy
from dataclasses import dataclass, field

from DataVisual import units
from PyMieSim.binary.Sets import CppDetectorSet


@dataclass
class BaseDetector():
    """
    Base class for constructing detectors in Mie scattering simulations, managing common attributes and methods.

    This class is responsible for formatting input parameters, calculating detector characteristics, and
    binding detectors to a simulation experiment. It should be subclassed to create specific detector types.

    Attributes:
        NA (Iterable[float]): Defines the numerical aperture, which is the range of angles the system can accept light.
        gamma_offset (Iterable[float]): Specifies the angular offset perpendicular to polarization (in degrees).
        phi_offset (Iterable[float]): Specifies the angular offset parallel to polarization (in degrees).
        polarization_filter (Iterable[float]): Sets the angle of the polarization filter (in degrees).
        sampling (Iterable[int]): Dictates the resolution for field sampling.

    This class is not intended for direct instantiation.
    """

    def __post_init__(self) -> NoReturn:
        """
        Called automatically after the class initialization to prepare the detector by formatting inputs, computing
        field arrays, setting up rotation angles, and initializing visualization and C++ bindings.
        """
        self.mapping = self.setup_mapping()
        self.initialize_binding()

    def setup_mapping(self) -> dict:
        """
        Creates a dictionary to map detector settings to their respective values, serving as a base for data visualization
        and C++ integration setup.

        Returns:
            dict: A mapping of detector settings to values.
        """
        return {
            'mode_number': None,
            'sampling': None,
            'rotation': None,
            'NA': None,
            'phi_offset': None,
            'gamma_offset': None,
            'polarization_filter': None,
        }

    def initialize_binding(self) -> NoReturn:
        """
        Prepares and initializes the C++ bindings necessary for the simulation, configuring the detector with simulation
        parameters.
        """
        self.binding_kwargs = dict(
            mode_number=numpy.atleast_1d(self.mode_number).astype(str),
            sampling=numpy.atleast_1d(self.sampling).astype(int),
            NA=numpy.atleast_1d(self.NA).astype(float),
            phi_offset=numpy.deg2rad(numpy.atleast_1d(self.phi_offset).astype(float)),
            gamma_offset=numpy.deg2rad(numpy.atleast_1d(self.gamma_offset).astype(float)),
            polarization_filter=numpy.deg2rad(numpy.atleast_1d(self.polarization_filter).astype(float)),
            rotation=numpy.deg2rad(numpy.atleast_1d(self.rotation)).astype(float),
            mean_coupling=self.mean_coupling,
            coherent=self.coherent
        )

        self.binding = CppDetectorSet(**self.binding_kwargs)

    def bind_to_experiment(self, experiment: Setup) -> NoReturn:
        """
        Binds the detector to a specified experimental setup, ensuring integration into the simulation workflow.

        Parameters:
            experiment (Setup): The experiment setup to which the detector is bound.
        """
        experiment.binding.set_detector(self.binding)

    def get_datavisual_table(self) -> NoReturn:
        """
        Compiles the detector's properties into a table format for data visualization.

        Returns:
            list: A list of formatted data visual elements representing the detector's properties.
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
    Represents a photodiode detector tailored for Mie scattering simulations, enhancing the BaseDetector with specific features.

    Attributes:
        coupling_mode (str): Indicates the mode coupling method, either 'point' or 'Mean'.
        coherent (bool): Specifies if the detection is coherent (True) or incoherent (False). Defaults to False.
        name (str): The name of the detector, initialized to "Photodiode". Not intended to be modified.

    Note:
        This class is specifically configured to simulate a photodiode detector within a Mie scattering experiment.
    """
    NA: Iterable[float] | float
    gamma_offset: Iterable[float] | float
    phi_offset: Iterable[float] | float
    polarization_filter: Iterable[float] | float
    sampling: Iterable[int] | int
    mean_coupling: bool = True
    rotation: Iterable[float] | float = field(default=0, init=False)
    coherent: bool = field(default=False, init=False)
    mode_number: str = field(default='NC00', init=False)


@dataclass
class CoherentMode(BaseDetector):
    """
    Specialized for coherent detection modes in Mie scattering simulations, this class extends BaseDetector.

    It manages the initialization and representation of coherent detection modes, specifically addressing their unique requirements.

    Attributes:
        mode_number (Iterable[str] | str): Designates the mode numbers involved in the detection.
        mean_coupling (bool): Indicates whether to use average coupling for calculations. Defaults to False.
        coherent (bool): Specifies if the detection is inherently coherent. Defaults to True.

    Note:
        This class is specifically designed to handle and simulate coherent detection modes.
    """
    mode_number: Iterable[str] | str
    rotation: Iterable[float] | float
    NA: Iterable[float] | float
    gamma_offset: Iterable[float] | float
    phi_offset: Iterable[float] | float
    polarization_filter: Iterable[float] | float
    sampling: Iterable[int] | int
    mean_coupling: bool = False
    coherent: bool = field(default=True, init=False)

    def __post_init__(self):
        """
        Initializes complex scalar field arrays to represent CoherentMode modes, preparing the detector for simulation.

        Returns:
            numpy.ndarray: An array of complex scalar fields representing the modes used in detection.
        """
        self.mode_number = numpy.atleast_1d(self.mode_number).astype(str)

        for idx, mode_name in enumerate(self.mode_number):
            mode_family_name = mode_name[:2]
            if mode_family_name not in ['LP', 'HG', 'LG', 'NC']:
                raise ValueError('Invalid mode family name, must be one of: LP, HG, LG, NC')

        super().__post_init__()

# -
