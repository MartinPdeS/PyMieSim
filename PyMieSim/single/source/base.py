#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydantic import validator
from PyMieSim.polarization import BasePolarization, Linear
from PyMieSim.units import Quantity, watt, meter, degree, watt, AU

config_dict = dict(
    kw_only=True,
    slots=True,
    extra='forbid',
    arbitrary_types_allowed=True
)


class BaseSource:
    @validator('NA', pre=True, always=True)
    def validate_NA(cls, value):
        """
        Ensures that diameter is Quantity objects with AU units."""
        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with AU units.")

        if not value.check(AU):
            raise ValueError(f"{value} must have AU units.")

        return value

    @validator('optical_power', pre=True, always=True)
    def validate_optical_power(cls, value):
        """
        Ensures that diameter is Quantity objects with power units."""
        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with power units.")

        if not value.check(watt):
            raise ValueError(f"{value} must have power units (watt).")

        return value

    @validator('wavelength', pre=True, always=True)
    def validate_wavelength(cls, value):
        """
        Ensures that diameter is Quantity objects with length units."""
        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with meters units.")

        if not value.check(meter):
            raise ValueError(f"{value} must have length units (meters).")

        return value

    @validator('polarization', pre=True, always=True)
    def validate_polarization(cls, value):
        """
        Ensures that polarization is well defined."""
        if isinstance(value, BasePolarization):
            return value

        if isinstance(value, Quantity):
            return Linear(value)

        raise ValueError(f"{value} must be a Linear or a Quantity with degree units.")


