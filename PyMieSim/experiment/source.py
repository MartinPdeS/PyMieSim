#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING, Iterable
if TYPE_CHECKING:
    from PyMieSim.experiment.setup import Setup

import numpy
from dataclasses import dataclass, field

from PyMieSim.physics import power_to_amplitude
from PyMieSim import polarization
from PyMieSim.binary.Sets import CppSourceSet
from DataVisual import units


@dataclass
class BaseSource:
    """
    Base class for light sources in PyMieSim experiments.

    Attributes:
        wavelength (Iterable): The wavelength(s) of the light source.
        polarization_value (Iterable): The polarization values of the light source, in degrees.
        name (str): The name of the source set, defaults to 'PlaneWave'.
    """
    wavelength: Iterable
    polarization_value: Iterable
    name: str = field(default='PlaneWave', init=False)

    def __post_init__(self):
        self.format_inputs()
        self.generate_polarization_attribute()
        self.generate_amplitude()
        self.generate_binding()
        self.build_x_parameters()

    def format_inputs(self):
        """Formats input wavelengths and polarization values into numpy arrays."""
        self.polarization_value = numpy.atleast_1d(self.polarization_value).astype(float)
        self.wavelength = numpy.atleast_1d(self.wavelength).astype(float)

    def generate_polarization_attribute(self):
        """
        Generates the polarization attribute based on specified polarization values.
        Raises an error for invalid polarization types not handled in subclasses.
        """
        self.polarization = polarization.Linear(*self.polarization_value)
        self.jones_vector = numpy.atleast_1d(self.polarization.jones_vector).astype(numpy.complex128).T

    def generate_amplitude(self):
        """Abstract method for generating amplitude, to be implemented by subclasses."""
        raise NotImplementedError("Subclass must implement this method.")

    def build_x_parameters(self):
        """Builds Xparameters for wavelength and linear polarization for DataVisual representations."""
        self.wavelength = units.Length(
            long_label='Wavelength',
            short_label=r'$\lambda$',
            base_values=self.wavelength,
            string_format='.2f'
        )

        self.linear_polarization = units.Degree(
            long_label='Polarization',
            short_label=r'Pol.',
            base_values=self.polarization_value,
            string_format='.2f'
        )

    def generate_binding(self):
        """Generates the C++ binding for the source set."""
        self.binding = CppSourceSet(
            wavelength=self.wavelength,
            jones_vector=self.jones_vector,
            amplitude=self.amplitude
        )

    def bind_to_experiment(self, experiment: Setup):
        """
        Binds the source set to an experiment.

        Parameters:
            experiment (Setup): The experiment setup to bind the source to.
        """
        experiment.binding.set_source(self.binding)

    def append_to_table(self, table: list) -> list:
        """
        Appends the wavelength and linear polarization to a given table.

        Parameters:
            table (list): The table to append data to.

        Returns:
            list: The updated table with appended data.
        """
        return [*table, self.wavelength, self.linear_polarization]


@dataclass
class Gaussian(BaseSource):
    """
    Represents a Gaussian light source with a specified numerical aperture and optical power.

    Inherits from BaseSource and adds specific attributes for Gaussian sources.

    Attributes:
        NA (Iterable): The numerical aperture(s) of the Gaussian source.
        optical_power (float): The optical power of the source, in Watts.
        polarization_type (str): The type of polarization, defaults to 'linear'.
    """
    NA: Iterable
    optical_power: float
    polarization_type: str = 'linear'

    def generate_amplitude(self):
        """Generates the amplitude of the Gaussian source based on its optical power and numerical aperture."""
        self.amplitude = power_to_amplitude(
            wavelength=self.wavelength,
            optical_power=self.optical_power,
            NA=self.NA
        )


@dataclass
class PlaneWave(BaseSource):
    """
    Represents a Plane Wave light source with a specified amplitude.

    Inherits from BaseSource and specifies amplitude directly.

    Attributes:
        amplitude (float): The amplitude of the plane wave, in Watts.
        polarization_type (str): The type of polarization, defaults to 'linear'.
    """
    amplitude: float
    polarization_type: str = 'linear'

    def generate_amplitude(self):
        """Sets the amplitude of the plane wave as a numpy array."""
        self.amplitude = numpy.atleast_1d(self.amplitude)
