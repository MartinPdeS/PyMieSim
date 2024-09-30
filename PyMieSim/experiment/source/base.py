#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from pydantic.dataclasses import dataclass
from pydantic import ConfigDict, field_validator
import pint_pandas

from PyMieSim import polarization
from PyMieSim.units import degree
from PyMieSim.binary.SetsInterface import CppSourceSet
from PyMieSim.units import Quantity, meter

config_dict = ConfigDict(
    kw_only=True,
    slots=True,
    extra='forbid',
    arbitrary_types_allowed=True
)

@dataclass
class BaseSource:
    """
    Base class for light sources in PyMieSim experiments.
    """

    def __post_init__(self):
        self.mapping = {}

        if not isinstance(self.polarization, polarization.BasePolarization):
            self.polarization = polarization.Linear(self.polarization)

        self._generate_binding_kwargs()

        binding_kwargs = {
            k: v.to_base_units().magnitude if isinstance(v, Quantity) else v for k, v in self.binding_kwargs.items()
        }

        self.binding = CppSourceSet(**binding_kwargs)

    @field_validator('wavelength', mode='before')
    def validate_length_quantity(cls, value):
        """
        Ensures that diameter is Quantity objects with length units."""
        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with meters units.")

        if not value.check(meter):
            raise ValueError(f"{value} must have length units (meters).")

        return numpy.atleast_1d(value)

    def _generate_mapping(self) -> None:
        """
        Appends the scatterer's properties to a given table for visualization purposes. This enables the
        representation of scatterer properties in graphical formats.

        Parameters:
            table (list): The table to which the scatterer's properties will be appended.

        Returns:
            list: The updated table with the scatterer's properties included.
        """
        self.mapping['wavelength'] = pint_pandas.PintArray(self.wavelength, dtype=self.wavelength.units)
        self.mapping['source_NA'] = pint_pandas.PintArray(self.NA, dtype=self.NA.units)
        self.mapping['optical_power'] = pint_pandas.PintArray(self.optical_power, dtype=self.optical_power.units)

        if hasattr(self.polarization, 'angle'):
            self.mapping['polarization'] = pint_pandas.PintArray(self.polarization.angle, dtype=degree)
        else:
            self.mapping['polarization'] = [repr(e) for e in self.polarization.element]
