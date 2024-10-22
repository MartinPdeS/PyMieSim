#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydantic import field_validator
from PyMieSim.polarization import BasePolarization, Linear
from PyMieSim.units import Quantity, watt, meter, AU

config_dict = dict(
    kw_only=True,
    slots=True,
    extra='forbid',
    arbitrary_types_allowed=True
)


class BaseSource:
    @field_validator('NA', mode='before')
    def _validate_NA(cls, value):
        """
        Ensures that diameter is Quantity objects with AU units."""
        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with AU units.")

        if not value.check(AU):
            raise ValueError(f"{value} must have AU units.")

        return value

    @field_validator('optical_power', mode='before')
    def _validate_optical_power(cls, value):
        """
        Ensures that diameter is Quantity objects with power units."""
        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with power units.")

        if not value.check(watt):
            raise ValueError(f"{value} must have power units (watt).")

        return value

    @field_validator('wavelength', mode='before')
    def _validate_wavelength(cls, value):
        """
        Ensures that diameter is Quantity objects with length units."""
        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with meters units.")

        if not value.check(meter):
            raise ValueError(f"{value} must have length units (meters).")

        return value

    @field_validator('polarization', mode='before')
    def _validate_polarization(cls, value):
        """
        Ensures that polarization is well defined."""
        if isinstance(value, BasePolarization):
            return value

        if isinstance(value, Quantity):
            return Linear(value)

        raise ValueError(f"{value} must be a Linear or a Quantity with degree units.")
