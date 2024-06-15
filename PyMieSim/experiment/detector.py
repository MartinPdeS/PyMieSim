#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy
from dataclasses import field
from pydantic.dataclasses import dataclass
from pydantic import ConfigDict
from PyMieSim.binary.Sets import CppDetectorSet
from PyMieSim.experiment import parameters
from typing import List, Union, NoReturn


config_dict = ConfigDict(
    kw_only=True,
    slots=True,
    extra='forbid',
    arbitrary_types_allowed=True
)


class BaseDetector():
    """
    Base class for constructing detectors in Mie scattering simulations, managing common attributes and methods.

    This class is responsible for formatting input parameters, calculating detector characteristics, and
    binding detectors to a simulation experiment. It should be subclassed to create specific detector types.

    Attributes:
        NA (List[float]): Defines the numerical aperture, which is the range of angles the system can accept light.
        gamma_offset (List[float]): Specifies the angular offset perpendicular to polarization (in degrees).
        phi_offset (List[float]): Specifies the angular offset parallel to polarization (in degrees).
        polarization_filter (List[float]): Sets the angle of the polarization filter (in degrees).
        sampling (List[int]): Dictates the resolution for field sampling.

    This class is not intended for direct instantiation.
    """

    def __post_init__(self) -> NoReturn:
        """
        Called automatically after the class initialization to prepare the detector by formatting inputs, computing
        field arrays, setting up rotation angles, and initializing visualization and C++ bindings.
        """
        self.mode_number = numpy.atleast_1d(self.mode_number).astype(str)
        self.sampling = numpy.atleast_1d(self.sampling).astype(int)
        self.NA = numpy.atleast_1d(self.NA).astype(float)
        self.phi_offset = numpy.deg2rad(numpy.atleast_1d(self.phi_offset).astype(float))
        self.gamma_offset = numpy.deg2rad(numpy.atleast_1d(self.gamma_offset).astype(float))
        self.polarization_filter = numpy.deg2rad(numpy.atleast_1d(self.polarization_filter).astype(float))
        self.rotation = numpy.deg2rad(numpy.atleast_1d(self.rotation)).astype(float)

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
            mode_number=self.mode_number,
            sampling=self.sampling,
            NA=self.NA,
            phi_offset=self.phi_offset,
            gamma_offset=self.gamma_offset,
            polarization_filter=self.polarization_filter,
            rotation=self.rotation,
            mean_coupling=self.mean_coupling,
            coherent=self.coherent
        )

        self.binding = CppDetectorSet(**self.binding_kwargs)

    def get_datavisual_table(self) -> NoReturn:
        """
        Compiles the detector's properties into a table format for data visualization.

        Returns:
            list: A list of formatted data visual elements representing the detector's properties.
        """

        for parameter_str in ['mode_number', 'sampling', 'rotation', 'NA', 'phi_offset', 'gamma_offset', 'polarization_filter']:
            parameters_values = getattr(self, parameter_str)
            self.mapping[parameter_str] = getattr(parameters, parameter_str)
            self.mapping[parameter_str].set_base_values(parameters_values)

        return [v for k, v in self.mapping.items() if v is not None]


@dataclass(config=config_dict)
class Photodiode(BaseDetector):
    """
    A photodiode detector tailored for Mie scattering simulations, enhancing the BaseDetector with specific features.

    Attributes:
        NA (Union[List[float], float]): Numerical aperture(s) for the detector.
        gamma_offset (Union[List[float], float]): Gamma offset(s) for the detector.
        phi_offset (Union[List[float], float]): Phi offset(s) for the detector.
        polarization_filter (Union[List[Optional[float]], Optional[float]]): Polarization filter(s) for the detector.
        sampling (Union[List[int], int]): Sampling rate(s) for the detector.
        mean_coupling (bool): Specifies if mean coupling is used. Defaults to True.
        rotation (Union[List[float], float]): Rotation angle(s) for the detector. Initialized to 0.
        coherent (bool): Indicates if the detection is coherent. Initialized to False.
        mode_number (str): Mode number of the detector. Initialized to 'NC00'.

    Notes:
        This class is specifically configured to simulate a photodiode detector within a Mie scattering experiment.
    """
    NA: Union[numpy.ndarray, List[float], float]
    gamma_offset: Union[numpy.ndarray, List[float], float]
    phi_offset: Union[numpy.ndarray, List[float], float]
    polarization_filter: Union[numpy.ndarray, List[float | None], float | None]
    sampling: Union[numpy.ndarray, List[int], int]
    mean_coupling: bool = True
    rotation: Union[numpy.ndarray, List[float] | float] = field(default=0, init=False)
    coherent: bool = field(default=False, init=False)
    mode_number: str = field(default='NC00', init=False)


@dataclass(config=config_dict)
class CoherentMode(BaseDetector):
    """
    Specialized for coherent detection modes in Mie scattering simulations, this class extends BaseDetector.

    It manages the initialization and representation of coherent detection modes, specifically addressing their unique requirements.

    Attributes:
        mode_number (List[str] | str): Designates the mode numbers involved in the detection.
        mean_coupling (bool): Indicates whether to use average coupling for calculations. Defaults to False.
        coherent (bool): Specifies if the detection is inherently coherent. Defaults to True.

    Note:
        This class is specifically designed to handle and simulate coherent detection modes.
    """
    mode_number: Union[numpy.ndarray, List[str], str]
    rotation: Union[numpy.ndarray, List[float], float]
    NA: Union[numpy.ndarray, List[float], float]
    gamma_offset: Union[numpy.ndarray, List[float], float]
    phi_offset: Union[numpy.ndarray, List[float], float]
    polarization_filter: Union[numpy.ndarray, List[float | None], float | None]
    sampling: Union[numpy.ndarray, List[int], int]
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

        super(CoherentMode, self).__post_init__()

# -
