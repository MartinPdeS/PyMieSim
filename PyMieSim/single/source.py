#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Union, NoReturn
import numpy
from pydantic.dataclasses import dataclass
from dataclasses import field

from PyMieSim.polarization import UnitPolarizationAngle
from PyMieSim.special_functions import NA_to_angle
from MPSPlots.render3D import SceneList as SceneList3D
from PyMieSim.binary.SourceInterface import BindedGaussian, BindedPlanewave, BindedBaseSource  # noqa: F401


config_dict = dict(
    kw_only=True,
    slots=True,
    extra='forbid'
)


@dataclass(config=config_dict)
class PlaneWave():
    """
    Represents a plane wave light source for light scattering simulations.

    Inherits from LightSource and specifies amplitude directly.

    Attributes:
        wavelength (float): Wavelength of the light field in meters.
        polarization (Union[UnitPolarizationAngle, float]): Polarization state of the light field, if float is given it is assumed Linear polarization of angle theta.
        amplitude (float): Amplitude of the electric field.
    """
    wavelength: float
    polarization: Union[UnitPolarizationAngle, float]
    amplitude: float = 1

    def __post_init__(self) -> NoReturn:
        if not isinstance(self.polarization, UnitPolarizationAngle):
            self.polarization = UnitPolarizationAngle(self.polarization)

        self.wavenumber = 2 * numpy.pi / self.wavelength

        self.binding = BindedPlanewave(
            wavelength=self.wavelength,
            amplitude=self.amplitude,
            jones_vector=self.polarization.jones_vector
        )

    def plot(self) -> SceneList3D:
        """
        Plots the structure of the PlaneWave source.

        Returns:
            SceneList3D: A 3D plotting scene object.
        """
        figure = SceneList3D()
        ax = figure.append_ax()
        ax.add_unit_sphere(opacity=0.3)
        ax.add_unit_axis(show_label=False)
        return figure


@dataclass(config=config_dict)
class Gaussian():
    """
    Represents a Gaussian light source for light scattering simulations, characterized by its optical power and numerical aperture.

    Attributes:
        wavelength (float): Wavelength of the light field in meters.
        polarization_value (float): Polarization state of the light field, if float is given it is assumed Linear polarization of angle theta.
        amplitude (float): Amplitude of the electric field.
        optical_power (float): Optical power of the source in Watts.
        NA (float): Numerical aperture of the source.
    """
    wavelength: float
    polarization: Union[UnitPolarizationAngle, float]
    optical_power: float
    NA: float
    amplitude: float = field(init=False, repr=False)

    def __post_init__(self) -> NoReturn:
        if not isinstance(self.polarization, UnitPolarizationAngle):
            self.polarization = UnitPolarizationAngle(self.polarization)

        self.wavenumber = 2 * numpy.pi / self.wavelength

        self.binding = BindedGaussian(
            wavelength=self.wavelength,
            NA=self.NA,
            optical_power=self.optical_power,
            jones_vector=self.polarization.jones_vector
        )

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
