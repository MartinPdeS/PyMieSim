#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pint_pandas
from typing import Union
from PyMieSim.single.source import (
    Gaussian as _,
)  # noqa:F401  # Necessary for pybind11 binding initialization

from TypedUnit import Length, Dimensionless, Power, Angle

from PyMieSim.experiment.source.base import BaseSource
from PyMieSim.experiment.utils import Sequential
from PyMieSim.single.polarization import BasePolarization, Linear
from PyMieSim.binary.interface_experiment import GaussianSourceSet


class Gaussian(BaseSource, Sequential):
    """
    Represents a Gaussian light source with a specified numerical aperture and optical power.

    Inherits from BaseSource and adds specific attributes for Gaussian sources.

    Parameters
    ----------
    wavelength : Quantity
        The wavelength(s) of the light source.
    polarization : Union[UnitPolarizationAngle, Quantity]
        Polarization state of the light field, if float is given it is assumed Linear polarization of angle theta.
    NA : Quantity
        The numerical aperture(s) of the Gaussian source.
    optical_power : Quantity
        The optical power of the source, in Watts.
    """

    attributes = ["wavelength", "NA", "optical_power", "polarization"]

    def __init__(
        self,
        wavelength: Length,
        polarization: Union[BasePolarization, Angle],
        NA: Dimensionless,
        optical_power: Power,
    ):
        self.mapping = {}

        self.NA = np.atleast_1d(NA)
        self.wavelength = np.atleast_1d(wavelength)
        self.optical_power = np.atleast_1d(optical_power)

        if not isinstance(polarization, BasePolarization):
            Angle.check(polarization)
            polarization = Linear(polarization)

        self.polarization = polarization

        self.binding_kwargs = dict(
            wavelength=self.wavelength,
            jones_vector=self.polarization.element,
            NA=self.NA,
            optical_power=self.optical_power,
            is_sequential=self.is_sequential,
        )

        self.set = GaussianSourceSet(**self.binding_kwargs)

