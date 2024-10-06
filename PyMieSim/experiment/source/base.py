#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from pydantic.dataclasses import dataclass
from pydantic import ConfigDict, field_validator
import pint_pandas
from typing import Union
from dataclasses import fields
from PyMieSim.polarization import BasePolarization, Linear
from PyMieSim.units import degree
from PyMieSim.binary.SetsInterface import CppSourceSet
from PyMieSim.units import Quantity, meter

config_dict = ConfigDict(
    kw_only=True,
    slots=True,
    extra='forbid',
    arbitrary_types_allowed=True
)

@dataclass(config=config_dict)
class BaseSource:
    """
    Base class for light sources in PyMieSim experiments.
    """
    wavelength: Quantity
    polarization: Union[BasePolarization, Quantity]

    def __post_init__(self):
        self.mapping = {}

        if not isinstance(self.polarization, BasePolarization):
            self.polarization = Linear(self.polarization)

        self._generate_binding_kwargs()

        binding_kwargs = {
            k: v.to_base_units().magnitude if isinstance(v, Quantity) else v for k, v in self.binding_kwargs.items()
        }

        self.binding = CppSourceSet(**binding_kwargs)

    @field_validator('wavelength', mode='before')
    def _validate_length_quantity(cls, value):
        """
        Ensures that diameter is Quantity objects with length units."""
        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with meters units.")

        if not value.check(meter):
            raise ValueError(f"{value} must have length units (meters).")

        return numpy.atleast_1d(value)

    def _generate_mapping(self) -> None:
        """
        Constructs a table of the scatterer's properties formatted for data visualization.
        This method populates the `mapping` dictionary with user-friendly descriptions and formats of the scatterer properties.

        Returns:
            list: A list of visual representations for each property in the `mapping` dictionary that has been populated.
        """

        for attr in [f.name for f in fields(self) if f.name != 'source']:
            values = getattr(self, attr)

            if values is None: continue

            # attr = attr.replace('_', ' ').capitalize()
            print(values)

            if hasattr(values, 'magnitude'):
                magnitude = values.magnitude
                units  = values.units
                self.mapping[attr] = pint_pandas.PintArray(magnitude, dtype=units)
            else:
                self.mapping[attr] = [repr(m) for m in values]

    # def _generate_mapping(self) -> None:
    #     """
    #     Appends the scatterer's properties to a given table for visualization purposes. This enables the
    #     representation of scatterer properties in graphical formats.

    #     Parameters:
    #         table (list): The table to which the scatterer's properties will be appended.

    #     Returns:
    #         list: The updated table with the scatterer's properties included.
    #     """
    #     self.mapping['wavelength'] = pint_pandas.PintArray(self.wavelength, dtype=self.wavelength.units)
    #     self.mapping['source_NA'] = pint_pandas.PintArray(self.NA, dtype=self.NA.units)
    #     self.mapping['optical_power'] = pint_pandas.PintArray(self.optical_power, dtype=self.optical_power.units)

    #     if hasattr(self.polarization, 'angle'):
    #         self.mapping['polarization'] = pint_pandas.PintArray(self.polarization.angle, dtype=degree)
    #     else:
    #         self.mapping['polarization'] = [repr(e) for e in self.polarization.element]
