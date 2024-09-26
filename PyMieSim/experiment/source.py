#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import List, Union
import numpy
from pydantic.dataclasses import dataclass
from pydantic import ConfigDict, validator
from dataclasses import field
import pint_pandas

from PyMieSim import polarization
from PyMieSim.units import degree
from PyMieSim.binary.SetsInterface import CppSourceSet
from PyMieSim.units import Quantity, meter

config_dict = ConfigDict(
    kw_only=True,
    slots=True,
    extra='forbid',
    arbitrary_types_allowed=True
)


@dataclass
class BaseSource:
    """
    Base class for light sources in PyMieSim experiments.
    """

    def __post_init__(self):
        self.mapping = {}

        if not isinstance(self.polarization, polarization.BasePolarization):
            self.polarization = polarization.Linear(self.polarization)

        self._generate_binding_kwargs()

        binding_kwargs = {
            k: v.to_base_units().magnitude if isinstance(v, Quantity) else v for k, v in self.binding_kwargs.items()
        }

        self.binding = CppSourceSet(**binding_kwargs)

    @validator('wavelength', pre=True, always=True)
    def validate_length_quantity(cls, value):
        """
        Ensures that diameter is Quantity objects with length units."""
        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with meters units.")

        if not value.check(meter):
            raise ValueError(f"{value} must have length units (meters).")

        return numpy.atleast_1d(value)

    def _generate_mapping(self) -> None:
        """
        Appends the scatterer's properties to a given table for visualization purposes. This enables the
        representation of scatterer properties in graphical formats.

        Parameters:
            table (list): The table to which the scatterer's properties will be appended.

        Returns:
            list: The updated table with the scatterer's properties included.
        """
        self.mapping['wavelength'] = pint_pandas.PintArray(self.wavelength, dtype=self.wavelength.units)
        self.mapping['polarization'] = pint_pandas.PintArray(self.polarization.angles, dtype=degree)
        self.mapping['source_NA'] = pint_pandas.PintArray(self.NA, dtype=self.NA.units)
        self.mapping['optical_power'] = pint_pandas.PintArray(self.optical_power, dtype=self.optical_power.units)



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

    @validator('wavelength', pre=True, always=True)
    def validate_length_quantity(cls, value):
        """
        Ensures that diameter is Quantity objects with length units."""
        if not isinstance(value, Quantity):
            raise ValueError(f"{value} must be a Quantity with meters units.")

        if not value.check(meter):
            raise ValueError(f"{value} must have length units (meters).")

        return numpy.atleast_1d(value)

    @validator('polarization', 'NA', 'optical_power', pre=True)
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
            jones_vector=self.polarization.jones_vector,
            NA=self.NA,
            optical_power=self.optical_power
        )


@dataclass(config=config_dict)
class PlaneWave(BaseSource):
    """
    Represents a Plane Wave light source with a specified amplitude.

    Inherits from BaseSource and specifies amplitude directly.

    Attributes:
        wavelength (Quantity): The wavelength(s) of the light source.
        polarization (Union[UnitPolarizationAngle, float]): Polarization state of the light field, if float is given it is assumed Linear polarization of angle theta.
        amplitude (Quantity): The amplitude of the plane wave, in Watts.
    """
    wavelength: Quantity
    polarization: Union[polarization.JonesVector, Quantity]
    amplitude: Quantity

    name: str = field(default='PlaneWave', init=False)

    def _generate_binding_kwargs(self) -> None:
        """
        Prepares the keyword arguments for the C++ binding based on the scatterer's properties. This
        involves evaluating material indices and organizing them into a dictionary for the C++ interface.

        Returns:
            None
        """
        self.binding_kwargs = dict(
            wavelength=self.wavelength,
            jones_vector=self.polarization.jones_vector,
            amplitude=self.amplitude
        )

# -
