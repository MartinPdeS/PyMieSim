#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy
from pydantic import field_validator
import pint_pandas
from dataclasses import fields
from PyMieSim.single.polarization import BasePolarization, Linear

from TypedUnit import Length, Angle, Power, Dimensionless


class BaseSource:
    """
    Base class for light sources in PyMieSim experiments.
    """

    @field_validator("wavelength", mode="plain")
    def _validate_length(cls, value):
        """
        Ensures that diameter is Quantity objects with length units.
        """
        Length.check(value)
        return numpy.atleast_1d(value)

    @field_validator("optical_power", mode="plain")
    def _validate_power(cls, value):
        """
        Ensure that arrays are properly converted to numpy arrays.
        """
        Power.check(value)
        return numpy.atleast_1d(value)

    @field_validator("NA", mode="plain")
    def _validate_arbitrary_units(cls, value):
        """
        Ensure that arrays are properly converted to numpy arrays.
        """
        Dimensionless.check(value)
        return numpy.atleast_1d(value)

    @field_validator("polarization", mode="before")
    def _validate_polarization(cls, value):
        """
        Ensures that polarization is well defined.
        """
        if isinstance(value, BasePolarization):
            return value

        Angle.check(value)

        return Linear(value)

    @field_validator("amplitude", mode="before")
    def _validate_array(cls, value):
        """Ensure that arrays are properly converted to numpy arrays."""
        if not isinstance(value, numpy.ndarray):
            value = numpy.atleast_1d(value)

        return value

    def _generate_mapping(self) -> None:
        """
        Constructs a table of the scatterer's properties formatted for data visualization.
        This method populates the `mapping` dictionary with user-friendly descriptions and formats of the scatterer properties.

        Returns:
            list: A list of visual representations for each property in the `mapping` dictionary that has been populated.
        """
        for attr in [f.name for f in fields(self) if f.name != "source"]:
            values = getattr(self, attr)

            if values is None:
                continue

            if hasattr(values, "magnitude"):
                magnitude = values.magnitude
                units = values.units
                self.mapping["source:" + attr] = pint_pandas.PintArray(
                    magnitude, dtype=units
                )
            else:
                self.mapping["source:" + attr] = [repr(m) for m in values]
