#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyMieSim.polarization import BasePolarization, Linear
from PyMieSim.units import Quantity, watt, meter, AU


class BaseSource:

    def _validate_AU(cls, value):
        """
        Ensures that diameter is Quantity objects with AU units.
        """
        if not isinstance(value, Quantity) or not value.check(AU):
            raise ValueError(f"{value} must be a Quantity with arbitrary units [AU].")

        return value

    def _validate_watt(cls, value):
        """
        Ensures that diameter is Quantity objects with power units.
        """
        if not isinstance(value, Quantity) or not value.check(watt):
            raise ValueError(f"{value} must be a Quantity with power units [watt].")

        return value

    def _validate_length(cls, value):
        """
        Ensures that diameter is Quantity objects with length units.
        """
        if not isinstance(value, Quantity) or not value.check(meter):
            raise ValueError(f"{value} must be a Quantity with meters units [meter].")

        return value

    def _validate_polarization(cls, value):
        """
        Ensures that polarization is well defined.
        """
        if isinstance(value, BasePolarization):
            return value

        if isinstance(value, Quantity):
            return Linear(value)

        raise ValueError(f"{value} must be a Linear or a Quantity with degree units.")
