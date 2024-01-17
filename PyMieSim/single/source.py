#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections.abc import Iterable

import numpy
from dataclasses import dataclass
from PyMieSim.polarization import LinearPolarization, JonesVector, RightCircularPolarization, LeftCircularPolarization
from PyMieSim.physics import power_to_amplitude
from PyMieSim.tools.special_functions import NA_to_angle

from MPSPlots.render3D import SceneList as SceneList3D


@dataclass(kw_only=True)
class PlaneWave():
    """
    Class representing plane wave beam as a light source for light scattering.
    """
    wavelength: float
    """ Wavelenght of the light field. """
    polarization_value: float
    """ Polarization of the light field in degree. """
    polarization_type: str = 'linear'
    """ How to interpret the polarization value """
    amplitude: float
    """ Amplitude of the electric field """

    def __post_init__(self):
        self.k = 2 * numpy.pi / self.wavelength

        self.amplitude = power_to_amplitude(
            wavelength=self.wavelength,
            optical_power=self.optical_power,
            NA=self.NA,
        )

        self.generate_polarization_attribute()

    def generate_polarization_attribute(self) -> None:
        match self.polarization_type.lower():
            case 'linear':
                self.polarization = self.interpret_linear_polarization(value=self.polarization_value)
            case 'jones vector':
                self.polarization = self.interpret_jones_vector(value=self.polarization_value)
            case 'circular':
                self.polarization = self.interpret_circular_polarization(value=self.polarization_value)
            case _:
                raise f'Invalid polarization type: {self.polarization_type}. Valid options ["linear", "jones vector"]'

    def interpret_jones_vector(self, value: Iterable) -> JonesVector:
        value = numpy.atleast_1d(value)

        assert value.size == 2, 'Jones vector is supposed to be a size 2 complex vector.'

        return JonesVector(value)

    def interpret_linear_polarization(self, value: float) -> LinearPolarization:
        return LinearPolarization(value)

    def interpret_circular_polarization(self, value: str) -> LinearPolarization:
        match value.lower():
            case 'right':
                return RightCircularPolarization()
            case 'left':
                return LeftCircularPolarization()

        return LinearPolarization(value)

    def plot(self) -> SceneList3D:
        """
        Plots the structure of the source.

        :returns:   3D plotting scene
        :rtype:     SceneList3D
        """
        max_angle = NA_to_angle(NA=self.NA)

        max_angle = numpy.rad2deg(max_angle)

        figure = SceneList3D()

        ax = figure.append_ax()

        ax.add_cone(
            center=(0.0, 0.0, 0.45),
            direction=(0.0, 0.0, -1.0),
            height=0.9,
            resolution=100,
            angle=max_angle,
            color='red',
            opacity=0.7
        )

        ax.add_unit_sphere(opacity=0.3)
        ax.add_unit_axis(show_label=False)

        return figure


@dataclass(kw_only=True)
class Gaussian():
    """
    Class representing plane wave beam as a light source for light scattering.
    """
    wavelength: float
    """ Wavelenght of the light field. """
    optical_power: float
    """ Optical power in unit of Watt """
    NA: float
    """ Numerical aperture of the source """
    polarization_value: float
    """ Polarization of the light field in degree. """
    polarization_type: str = 'linear'
    """ How to interpret the polarization value """

    def __post_init__(self):
        self.k = 2 * numpy.pi / self.wavelength

        self.amplitude = power_to_amplitude(
            wavelength=self.wavelength,
            optical_power=self.optical_power,
            NA=self.NA,
        )

        self.generate_polarization_attribute()

    def generate_polarization_attribute(self) -> None:
        match self.polarization_type.lower():
            case 'linear':
                self.polarization = self.interpret_linear_polarization(value=self.polarization_value)
            case 'jones vector':
                self.polarization = self.interpret_jones_vector(value=self.polarization_value)
            case 'circular':
                self.polarization = self.interpret_circular_polarization(value=self.polarization_value)
            case _:
                raise f'Invalid polarization type: {self.polarization_type}. Valid options ["linear", "jones vector"]'

    def interpret_jones_vector(self, value: Iterable) -> JonesVector:
        value = numpy.atleast_1d(value)

        assert value.size == 2, 'Jones vector is supposed to be a size 2 complex vector.'

        return JonesVector(value)

    def interpret_linear_polarization(self, value: float) -> LinearPolarization:
        return LinearPolarization(value)

    def interpret_circular_polarization(self, value: str) -> LinearPolarization:
        match value.lower():
            case 'right':
                return RightCircularPolarization()
            case 'left':
                return LeftCircularPolarization()

        return LinearPolarization(value)

    def plot(self) -> SceneList3D:
        """
        Plots the structure of the source.

        :returns:   3D plotting scene
        :rtype:     SceneList3D
        """
        max_angle = NA_to_angle(NA=self.NA)

        max_angle = numpy.rad2deg(max_angle)

        figure = SceneList3D()

        ax = figure.append_ax()

        ax.add_cone(
            center=(0.0, 0.0, 0.45),
            direction=(0.0, 0.0, -1.0),
            height=0.9,
            resolution=100,
            angle=max_angle,
            color='red',
            opacity=0.7
        )

        ax.add_unit_sphere(opacity=0.3)
        ax.add_unit_axis(show_label=False)

        return figure


# -
