#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pint_pandas
from typing import Optional
from pydantic.dataclasses import dataclass
from pydantic import ConfigDict, field_validator
from PyMieSim.binary.SetsInterface import CppDetectorSet
from PyMieSim.units import Quantity, radian, degree, AU
from typing import Any


config_dict = ConfigDict(
    kw_only=True,
    slots=True,
    extra='forbid',
    arbitrary_types_allowed=True
)

@dataclass(config=config_dict)
class BaseDetector:
    """
    Base class for detectors in Mie scattering simulations, handling common attributes and methods.

    This class handles the initialization and binding of a detector to a simulation, formatting inputs, and managing
    the C++ bindings needed for high-performance calculations. It should be subclassed to create specific types of detectors.

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

    This class is not intended to be instantiated directly.
    """
    mode_number: str
    NA: Quantity
    gamma_offset: Quantity
    phi_offset: Quantity
    rotation: Quantity
    mean_coupling: bool
    coherent: bool
    sampling: Optional[Quantity] = 200 * AU
    polarization_filter: Optional[Quantity] = np.nan * degree


    @field_validator('polarization_filter', mode='before')
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

    @field_validator('gamma_offset', 'phi_offset', 'rotation', mode='before')
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

    @field_validator('NA', 'sampling', 'mode_number', mode='before')
    def validate_au_quantity(cls, value):
        """Ensure that arrays are properly converted to numpy arrays."""
        if not isinstance(value, np.ndarray):
            value = np.atleast_1d(value)

        return value

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