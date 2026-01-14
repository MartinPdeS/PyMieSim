#!/usr/bin/env python
# -*- coding: utf-8 -*-

from TypedUnit import Angle
from pydantic import field_validator


class BaseSource:
    @classmethod
    def check(cls, value):
        if not isinstance(value, cls):
            raise AssertionError(f"Expected instance of {cls.__name__} but got {type(value).__name__}")

        return value

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
