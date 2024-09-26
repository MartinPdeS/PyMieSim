#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from dataclasses import field
import pint_pandas
from pydantic.dataclasses import dataclass
from pydantic import ConfigDict, validator
from PyMieSim.binary.SetsInterface import CppDetectorSet
from PyMieSim.units import Quantity, radian, degree
from typing import List, Union, Optional, Any


config_dict = ConfigDict(
    kw_only=True,
    slots=True,
    extra='forbid',
    arbitrary_types_allowed=True
)



class BaseDetector:
    """
    Base class for detectors in Mie scattering simulations, handling common attributes and methods.

    This class handles the initialization and binding of a detector to a simulation, formatting inputs, and managing
    the C++ bindings needed for high-performance calculations. It should be subclassed to create specific types of detectors.

    Attributes:
        NA (List[float]): Numerical aperture(s) defining the range of angles the system can accept light.
        gamma_offset (Quantity): Angular offset perpendicular to polarization (in degrees).
        phi_offset (Quantity): Angular offset parallel to polarization (in degrees).
        polarization_filter (List[float]): Angle(s) of the polarization filter (in degrees).
        sampling (List[int]): Resolution for field sampling.

    This class is not intended to be instantiated directly.
    """

    def __post_init__(self) -> None:
        """Post-initialization: sets up mapping and binds the detector to the C++ engine."""
        self.mapping = dict.fromkeys(['mode_number', 'sampling', 'rotation', 'NA', 'phi_offset', 'gamma_offset', 'polarization_filter'])

        self._initialize_binding()

    def _initialize_binding(self) -> None:
        """
        Initializes C++ bindings for the simulation, configuring the detector with provided parameters.
        """
        self.binding_kwargs = {
            "mode_number": self.mode_number,
            "sampling": self.sampling,
            "NA": self.NA,
            "polarization_filter": self.polarization_filter.to(radian).magnitude if self.polarization_filter is not None else np.nan,
            "phi_offset": self.phi_offset.to(radian).magnitude,
            "gamma_offset": self.gamma_offset.to(radian).magnitude,
            "rotation": self.rotation.to(radian).magnitude,
        }

        # Ensure all values are at least 1D arrays
        self.binding_kwargs = {k: np.atleast_1d(v) for k, v in self.binding_kwargs.items()}

        # Additional detector settings
        self.binding_kwargs["mean_coupling"] = self.mean_coupling
        self.binding_kwargs["coherent"] = self.coherent

        self.binding = CppDetectorSet(**self.binding_kwargs)

    def _generate_mapping(self) -> list:
        """
        Appends the scatterer's properties to a given table for visualization purposes. This enables the
        representation of scatterer properties in graphical formats.

        Parameters:
            table (list): The table to which the scatterer's properties will be appended.

        Returns:
            list: The updated table with the scatterer's properties included.
        """
        self.mapping['mode_number'] = self.mode_number
        self.mapping['rotation'] = pint_pandas.PintArray(self.rotation, dtype=self.rotation.units)
        self.mapping['NA'] = pint_pandas.PintArray(self.NA, dtype=self.NA.units)
        self.mapping['phi_offset'] = pint_pandas.PintArray(self.phi_offset, dtype=self.phi_offset.units)
        self.mapping['gamma_offset'] = pint_pandas.PintArray(self.gamma_offset, dtype=self.gamma_offset.units)
        self.mapping['sampling'] = pint_pandas.PintArray(self.sampling, dtype=self.sampling.units)
        self.mapping['polarization_filter'] = self.polarization_filter

    @validator('polarization_filter', pre=True, always=True)
    def validate_polarization(cls, value):
        """
        Ensures that gamma_offset, phi_offset, polarization_filter, and rotation are Quantity objects with angle units.
        Converts them to numpy arrays after validation.
        """
        if value is None:
            value = np.nan * degree

        value = np.atleast_1d(value).astype(float)

        if not value.check(degree):
            raise ValueError(f"{value} must have angle units (degree or radian).")

        return value

    @validator('gamma_offset', 'phi_offset', 'rotation', pre=True, always=True)
    def validate_angle_quantity(cls, value):
        """
        Ensures that gamma_offset, phi_offset, and rotation are Quantity objects with angle units.
        Converts them to numpy arrays after validation.
        """
        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with angle units.")

        if not value.check(degree):
            raise ValueError(f"{value} must have angle units (degree or radian).")

        return np.atleast_1d(value).astype(float)

    @validator('NA', 'sampling', 'mode_number', pre=True, always=True)
    def validate_array(cls, value):
        """Ensure that arrays are properly converted to numpy arrays."""
        if not isinstance(value, np.ndarray):
            value = np.atleast_1d(value)

        return value


@dataclass(config=config_dict)
class Photodiode(BaseDetector):
    """
    Photodiode detector tailored for Mie scattering simulations, extending BaseDetector with specific features.

    Attributes:
        NA (Union[List[float], float]): Numerical aperture(s) of the detector.
        gamma_offset (Quantity): Gamma angular offset (in degrees).
        phi_offset (Quantity): Phi angular offset (in degrees).
        polarization_filter (Optional[Quantity]): Polarization filter angle (in degrees).
        sampling (Union[List[int], int]): Sampling rate(s) for the detector.
        mean_coupling (bool): Whether mean coupling is used. Defaults to True.
        rotation (Quantity): Rotation angle of the detector. Defaults to 0 degrees.
        coherent (bool): Indicates if the detection is coherent. Defaults to False.
        mode_number (str): Mode number of the detector. Defaults to 'NC00'.
    """
    NA: Quantity
    gamma_offset: Quantity
    phi_offset: Quantity
    sampling: Quantity
    polarization_filter: Optional[Quantity | None] = None
    mean_coupling: Optional[bool] = True

    coherent: bool = field(default=False, init=False)
    mode_number: Union[np.ndarray, List[str], str] = field(default='NC00', init=False)
    rotation: Quantity = field(default=[0] * degree, init=False)

@dataclass(config=config_dict)
class CoherentMode(BaseDetector):
    """
    Coherent mode detector for Mie scattering simulations, handling coherent detection modes.

    This detector is designed specifically for coherent modes, such as LP, HG, LG, and NC modes, which require specific handling
    in Mie scattering experiments.

    Attributes:
        mode_number (Union[List[str], str]): Mode number(s) involved in the detection.
        NA (Union[List[float], float]): Numerical aperture(s) of the detector.
        gamma_offset (Quantity): Gamma angular offset (in degrees).
        phi_offset (Quantity): Phi angular offset (in degrees).
        rotation (Quantity): Rotation angle of the detector.
        sampling (Union[List[int], int]): Sampling rate(s) for the detector.
        polarization_filter (Optional[Quantity]): Polarization filter angle (in degrees).
        mean_coupling (Optional[bool]): Whether mean coupling is used. Defaults to False.
        coherent (bool): Specifies if the detection is coherent. Defaults to True.
    """
    mode_number: Union[np.ndarray, List[str], str]
    NA: Quantity
    gamma_offset: Quantity
    phi_offset: Quantity
    rotation: Quantity
    sampling: Quantity
    polarization_filter: Optional[Quantity] = None
    mean_coupling: Optional[bool] = False

    coherent: bool = field(default=True, init=False)

    @validator('mode_number', pre=True)
    def validate_mode_number(cls, mode_number):
        """Ensure mode numbers are valid and belong to supported families."""
        mode_number = np.atleast_1d(mode_number).astype(str)
        for mode in mode_number:
            if mode[:2] not in ['LP', 'HG', 'LG', 'NC']:
                raise ValueError(f'Invalid mode family {mode[:2]}. Must be one of: LP, HG, LG, NC')
        return mode_number
