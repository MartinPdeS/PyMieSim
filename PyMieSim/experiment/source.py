#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import List, Union, NoReturn
import numpy
from pydantic.dataclasses import dataclass
from pydantic import ConfigDict
from dataclasses import field

from PyMieSim.experiment import parameters
from PyMieSim import polarization
from PyMieSim.binary.Sets import CppSourceSet

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
        self.mapping = {
            'wavelength': None,
            'polarization': None,
            'NA': None,
            'optical_power': None
        }

        self.format_inputs()
        self.generate_binding()

    def format_inputs(self):
        """Formats input wavelengths and polarization values into numpy arrays."""
        self.wavelength = numpy.atleast_1d(self.wavelength).astype(float)

        if not isinstance(self.polarization, polarization.JonesVector):
            self.polarization = polarization.Linear(self.polarization)

        self.NA = numpy.atleast_1d(self.NA).astype(float)
        self.optical_power = numpy.atleast_1d(self.optical_power).astype(float)

    def get_datavisual_table(self) -> NoReturn:
        """
        Appends the scatterer's properties to a given table for visualization purposes. This enables the
        representation of scatterer properties in graphical formats.

        Parameters:
            table (list): The table to which the scatterer's properties will be appended.

        Returns:
            list: The updated table with the scatterer's properties included.
        """
        self.mapping['wavelength'] = parameters.wavelength
        self.mapping['wavelength'].set_base_values(self.wavelength)

        self.mapping['polarization'] = parameters.polarization
        self.mapping['polarization'].set_base_values(self.polarization.elements)

        self.mapping['NA'] = parameters.NA_source
        self.mapping['NA'].set_base_values(self.NA)

        self.mapping['optical_power'] = parameters.optical_power
        self.mapping['optical_power'].set_base_values(self.optical_power)

        return [v for k, v in self.mapping.items() if v is not None]


@dataclass(config=config_dict)
class Gaussian(BaseSource):
    """
    Represents a Gaussian light source with a specified numerical aperture and optical power.

    Inherits from BaseSource and adds specific attributes for Gaussian sources.

    Attributes:
        wavelength (List): The wavelength(s) of the light source.
        polarization (List): The polarization values of the light source, in degrees.
        NA (List): The numerical aperture(s) of the Gaussian source.
        optical_power (float): The optical power of the source, in Watts.
    """
    wavelength: Union[numpy.ndarray, List[float], float]
    polarization: Union[polarization.JonesVector, numpy.ndarray, List[float], float]
    NA: Union[numpy.ndarray, List[float], float]
    optical_power: Union[numpy.ndarray, List[float], float]
    name: str = field(default='PlaneWave', init=False)

    def generate_binding(self) -> NoReturn:
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

        self.binding = CppSourceSet(**self.binding_kwargs)


@dataclass(config=config_dict)
class PlaneWave(BaseSource):
    """
    Represents a Plane Wave light source with a specified amplitude.

    Inherits from BaseSource and specifies amplitude directly.

    Attributes:
        wavelength (List): The wavelength(s) of the light source.
        polarization (List): The polarization values of the light source, in degrees.
        amplitude (float): The amplitude of the plane wave, in Watts.
    """
    wavelength: Union[numpy.ndarray, List[float], float]
    polarization: Union[polarization.JonesVector, numpy.ndarray, List[float], float]
    amplitude: Union[numpy.ndarray, List[float], float]
    name: str = field(default='PlaneWave', init=False)

    def generate_binding(self) -> NoReturn:
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

        self.binding = CppSourceSet(**self.binding_kwargs)

# -
