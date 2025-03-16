#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy
from pydantic import field_validator
import pint_pandas
from dataclasses import fields
from PyMieSim.polarization import BasePolarization, Linear
from PyMieSim.binary.interface_sets import CppSourceSet
from PyMieSim.units import Quantity, meter, watt, AU, degree


class BaseSource:
    """
    Base class for light sources in PyMieSim experiments.
    """
    def _generate_binding(self):
        self.mapping = {}

        if not isinstance(self.polarization, BasePolarization):
            self.polarization = Linear(self.polarization)

        self._generate_binding_kwargs()

        binding_kwargs = {
            k: v.to_base_units().magnitude if isinstance(v, Quantity) else v for k, v in self.binding_kwargs.items()
        }

        self.binding = CppSourceSet(**binding_kwargs)

    @field_validator('wavelength', mode='plain')
    def _validate_length(cls, value):
        """
        Ensures that diameter is Quantity objects with length units.
        """
        if not isinstance(value, Quantity) or not value.check(meter):
            raise ValueError(f"{value} must be a Quantity with meters units [meter].")

        return numpy.atleast_1d(value)

    @field_validator('optical_power', mode='plain')
    def _validate_power(cls, value):
        """
        Ensure that arrays are properly converted to numpy arrays.
        """
        if not isinstance(value, Quantity) or not value.check(watt):
            raise ValueError(f"{value} must be a Quantity with power units [watt].")

        if not isinstance(value, numpy.ndarray):
            value = numpy.atleast_1d(value)

        return value

    @field_validator('NA', mode='plain')
    def _validate_arbitrary_units(cls, value):
        """
        Ensure that arrays are properly converted to numpy arrays.
        """
        if not isinstance(value, Quantity) or not value.check(AU):
            raise ValueError(f"{value} must be a Quantity with arbitrary units [AU].")

        if not isinstance(value, numpy.ndarray):
            value = numpy.atleast_1d(value)

        return value

    @field_validator('polarization', mode='before')
    def _validate_polarization(cls, value):
        """
        Ensures that polarization is well defined.
        """
        if (not isinstance(value, Quantity) or not value.check(degree)) and not isinstance(value, BasePolarization):
            raise ValueError(f"{value} must be a Quantity with degree units [degree] or a Polarization instance.")

        if isinstance(value, BasePolarization):
            return value

        if isinstance(value, Quantity):
            return Linear(value)

    def _generate_mapping(self) -> None:
        """
        Constructs a table of the scatterer's properties formatted for data visualization.
        This method populates the `mapping` dictionary with user-friendly descriptions and formats of the scatterer properties.

        Returns:
            list: A list of visual representations for each property in the `mapping` dictionary that has been populated.
        """
        for attr in [f.name for f in fields(self) if f.name != 'source']:
            values = getattr(self, attr)

            if values is None:
                continue

            if hasattr(values, 'magnitude'):
                magnitude = values.magnitude
                units = values.units
                self.mapping["source:" + attr] = pint_pandas.PintArray(magnitude, dtype=units)
            else:
                self.mapping["source:" + attr] = [repr(m) for m in values]
