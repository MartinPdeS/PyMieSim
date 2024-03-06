#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections.abc import Iterable
import numpy
from dataclasses import dataclass

from PyMieSim.physics import power_to_amplitude
from PyMieSim import polarization
from PyMieSim.tools.special_functions import NA_to_angle
from MPSPlots.render3D import SceneList as SceneList3D


@dataclass(kw_only=True)
class LightSource:
    """
    Abstract class for light sources in light scattering simulations.

    Attributes:
        wavelength (float): Wavelength of the light field in meters.
        polarization_value (float): Polarization state of the light field.
        polarization_type (str): Specifies how the polarization_value should be interpreted ('linear', 'jones vector', 'circular').
        amplitude (float): Amplitude of the electric field.
    """
    wavelength: float
    polarization_value: float
    polarization_type: str = 'linear'

    def __post_init__(self):
        self.k = 2 * numpy.pi / self.wavelength  # Wave number
        self.generate_polarization_attribute()

    def generate_polarization_attribute(self) -> None:
        """
        Generates the polarization attribute based on the specified polarization_type and polarization_value.
        """
        match self.polarization_type.lower():
            case 'linear':
                self.polarization = polarization.Linear(self.polarization_value)
            case 'jones vector':
                self.polarization = self.interpret_jones_vector(self.polarization_value)
            case 'circular':
                self.polarization = self.interpret_circular_polarization(self.polarization_value)
            case _:
                raise ValueError(f'Invalid polarization type: {self.polarization_type}. Supported options are "linear", "jones vector", "circular".')

    def interpret_jones_vector(self, value: Iterable) -> polarization.JonesVector:
        """
        Interprets the given value as a Jones vector.

        Parameters:
            value (Iterable): A size 2 iterable representing the Jones vector.

        Returns:
            polarization.JonesVector: The Jones vector representation of the polarization.
        """
        value = numpy.atleast_1d(value)
        assert value.size == 2, 'Jones vector must be a size 2 complex vector.'
        return polarization.JonesVector(value)

    def interpret_circular_polarization(self, value: str):
        """
        Interprets the given value as circular polarization.

        Parameters:
            value (str): 'right' or 'left' indicating the circular polarization direction.

        Returns:
            polarization.RightCircular or polarization.LeftCircular: The circular polarization object.
        """
        match value.lower():
            case 'right':
                return polarization.RightCircular()
            case 'left':
                return polarization.LeftCircular()
            case _:
                raise ValueError('Circular polarization value must be either "right" or "left".')

    def plot(self) -> SceneList3D:
        """
        Abstract method for plotting the structure of the source. To be implemented by subclasses.
        """
        raise NotImplementedError("Subclass must implement this method.")


@dataclass(kw_only=True)
class PlaneWave(LightSource):
    """
    Represents a plane wave light source for light scattering simulations.

    Inherits from LightSource and specifies amplitude directly.
    """

    amplitude: float | None = None

    def __post_init__(self):
        super().__post_init__()
        # Assuming amplitude is directly specified; no additional calculation required for PlaneWave.

    def plot(self) -> SceneList3D:
        """
        Plots the structure of the PlaneWave source.

        Returns:
            SceneList3D: A 3D plotting scene object.
        """
        # Example plotting logic; adjust as needed for your specific visualization requirements.
        figure = SceneList3D()
        ax = figure.append_ax()
        ax.add_unit_sphere(opacity=0.3)
        ax.add_unit_axis(show_label=False)
        return figure


@dataclass(kw_only=True)
class Gaussian(LightSource):
    """
    Represents a Gaussian light source for light scattering simulations, characterized by its optical power and numerical aperture.

    Attributes:
        optical_power (float): Optical power of the source in Watts.
        NA (float): Numerical aperture of the source.
    """
    optical_power: float
    NA: float

    def __post_init__(self):
        self.amplitude = power_to_amplitude(
            wavelength=self.wavelength,
            optical_power=self.optical_power,
            NA=self.NA
        )

        super().__post_init__()

    def plot(self) -> SceneList3D:
        """
        Plots the structure of the Gaussian source.

        Returns:
            SceneList3D: A 3D plotting scene object.
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
