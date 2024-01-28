#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from PyMieSim.experiment.setup import Setup
    from collections.abc import Iterable

import numpy
from dataclasses import dataclass, field

from DataVisual import Xparameter

from PyMieSim.physics import power_to_amplitude
from PyMieSim.polarization import LinearPolarization
import PyMieSim.datavisual_x_parameters as Kwargs

from PyMieSim.binary.Sets import CppSourceSet


@dataclass
class BaseSource(object):
    wavelength: Iterable
    """ Wavelenght of the light field. """
    polarization_value: Iterable
    """ Polarization of the light field in degree. """
    name: str = field(default='PlaneWave', init=False)
    """ Name of the set """

    def __post_init__(self):
        self.format_inputs()

        self.generate_polarization_attribute()

        self.generate_amplitude()

        self.build_x_parameters()

        self.generate_binding()

    def generate_binding(self) -> None:
        """
        Generate the C++ binding

        :returns:   No return
        :rtype:     None
        """
        self.binding = CppSourceSet(
            wavelength=self.wavelength.values,
            jones_vector=self.jones_vector,
            amplitude=self.amplitude
        )

    def bind_to_experiment(self, experiment: Setup) -> None:
        """
        Bind this specific set to a Setup experiment.

        :param      experiment:  The experiment
        :type       experiment:  Setup

        :returns:   No return
        :rtype:     None
        """
        experiment.binding.set_source(self.binding)

    def append_to_table(self, table: list) -> list:
        """
        Append elements to the xTable from the DataVisual library for the plottings.

        :param      table:  The table
        :type       table:  list

        :returns:   The updated list
        :rtype:     list
        """
        return [*table, self.wavelength, self.linear_polarization]

    def build_x_parameters(self) -> None:
        """
        Builds the parameters that will be passed in XTable for DataVisual.

        :returns:   No Return
        :rtype:     None
        """
        self.wavelength = Xparameter(values=self.wavelength, **Kwargs.wavelength)

        self.linear_polarization = Xparameter(
            values=self.jones_vector,
            representation=self.polarization.angle_list,
            **Kwargs.linear_polarization
        )

    def format_inputs(self) -> None:
        """
        Format the inputs given by the user into numpy array. Those inputs are subsequently
        sent to the cpp binding.

        :returns:   No return
        :rtype:     None
        """
        self.polarization_value = numpy.atleast_1d(self.polarization_value).astype(float)

        self.wavelength = numpy.atleast_1d(self.wavelength).astype(float)

    def generate_polarization_attribute(self) -> None:
        match self.polarization_type.lower():
            case 'linear':
                self.polarization = LinearPolarization(*self.polarization_value)
            case _:
                raise f'Invalid polarization type: {self.polarization_type}. Valid options ["linear", "jones vector"]'

        self.jones_vector = numpy.atleast_1d(self.polarization.jones_vector).astype(numpy.complex128).T


@dataclass
class Gaussian(BaseSource):
    NA: Iterable
    """ Numerical aperture of the source """
    optical_power: float
    """ Optical power in unit of Watt """
    polarization_type: str = 'linear'
    """ How to interpret the polarization value """

    def generate_amplitude(self) -> None:
        self.amplitude = power_to_amplitude(
            wavelength=self.wavelength,
            optical_power=self.optical_power,
            NA=self.NA,
        )


@dataclass
class PlaneWave(BaseSource):
    amplitude: float
    """ Optical power in unit of Watt """
    polarization_type: str = 'linear'
    """ How to interpret the polarization value """

    def generate_amplitude(self) -> None:
        self.amplitude = numpy.atleast_1d(self.amplitude)
