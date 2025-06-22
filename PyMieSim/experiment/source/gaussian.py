#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import Union
from PyMieSim.single.source import Gaussian as _  # noqa:F401  # Necessary for pybind11 binding initialization
from pydantic.dataclasses import dataclass
from PyMieSim.units import Quantity
from PyMieSim.experiment.source.base import BaseSource
from PyMieSim.experiment.utils import config_dict, Sequential
from PyMieSim.polarization import BasePolarization, Linear
from PyMieSim.binary.interface_sets import CppGaussianSourceSet

@dataclass(config=config_dict)
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
    NA: Quantity
    optical_power: Quantity
    wavelength: Quantity
    polarization: Union[BasePolarization, Quantity]

    def _generate_binding(self) -> None:
        """
        Prepares the keyword arguments for the C++ binding based on the scatterer's properties. This
        involves evaluating material indices and organizing them into a dictionary for the C++ interface.

        """
        self.mapping = {}

        if not isinstance(self.polarization, BasePolarization):
            self.polarization = Linear(self.polarization)

        self.binding_kwargs = dict(
            wavelength=self.wavelength,
            jones_vector=self.polarization.element,
            NA=self.NA,
            optical_power=self.optical_power,
            is_sequential=self.is_sequential
        )

        self.binding = CppGaussianSourceSet(
            wavelength=self.wavelength.to('meter').magnitude,
            jones_vector=self.polarization.element,
            NA=self.NA.to('dimensionless').magnitude,
            optical_power=self.optical_power.to('watt').magnitude,
            is_sequential=self.is_sequential
        )