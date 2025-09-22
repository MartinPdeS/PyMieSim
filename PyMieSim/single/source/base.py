#!/usr/bin/env python
# -*- coding: utf-8 -*-

from TypedUnit import Angle
from pydantic import field_validator


class BaseSource:
    @field_validator("polarization", mode="plain")
    def _validate_source_polarization(cls, value):
        """
        Ensures that polarization is well defined.
        """
        from PyMieSim.single.polarization import BasePolarization, Linear

        if isinstance(value, BasePolarization):
            return value

        Angle.check(value)

        return Linear(value)
