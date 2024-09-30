#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Union
import numpy
from pydantic.dataclasses import dataclass
from pydantic import field_validator
from dataclasses import field
from PyMieSim import polarization
from PyMieSim.units import Quantity, meter
from PyMieSim.experiment.source.base import BaseSource, config_dict

@dataclass(config=config_dict)
class Gaussian(BaseSource):
    """
    Represents a Gaussian light source with a specified numerical aperture and optical power.

    Inherits from BaseSource and adds specific attributes for Gaussian sources.

    Attributes:
        wavelength (Quantity): The wavelength(s) of the light source.
        polarization (Union[UnitPolarizationAngle, float]): Polarization state of the light field, if float is given it is assumed Linear polarization of angle theta.
        NA (List): The numerical aperture(s) of the Gaussian source.
        optical_power (float): The optical power of the source, in Watts.
    """
    wavelength: Quantity
    polarization: Union[polarization.BasePolarization, Quantity]
    NA: Quantity
    optical_power: Quantity

    name: str = field(default='PlaneWave', init=False)

    @field_validator('wavelength', mode='before')
    def validate_length_quantity(cls, value):
        """
        Ensures that diameter is Quantity objects with length units."""
        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with meters units.")

        if not value.check(meter):
            raise ValueError(f"{value} must have length units (meters).")

        return numpy.atleast_1d(value)

    @field_validator('NA', 'optical_power', mode='before')
    def validate_array(cls, value):
        """Ensure that arrays are properly converted to numpy arrays."""
        if not isinstance(value, numpy.ndarray):
            value = numpy.atleast_1d(value)

        return value

    def _generate_binding_kwargs(self) -> None:
        """
        Prepares the keyword arguments for the C++ binding based on the scatterer's properties. This
        involves evaluating material indices and organizing them into a dictionary for the C++ interface.

        Returns:
            None
        """
        self.binding_kwargs = dict(
            wavelength=self.wavelength,
            jones_vector=self.polarization.element,
            NA=self.NA,
            optical_power=self.optical_power
        )
