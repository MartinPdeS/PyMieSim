#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from typing import Union
from PyMieSim.single.source import (
    PlaneWave as _,
)  # noqa:F401  # Necessary for pybind11 binding initialization

from TypedUnit import Length, Angle, ElectricField

from PyMieSim.experiment.source.base import BaseSource
from PyMieSim.single.polarization import BasePolarization, Linear
from PyMieSim.experiment.utils import Sequential
from PyMieSim.binary.interface_experiment import PlaneWaveSourceSet

class PlaneWave(BaseSource, Sequential):
    """
    Represents a Plane Wave light source with a specified amplitude.

    Inherits from BaseSource and specifies amplitude directly.

    Parameters
    ----------
    wavelength : Quantity
        The wavelength(s) of the light source.
    polarization : Union[UnitPolarizationAngle, Quantity]
        Polarization state of the light field, if float is given it is assumed Linear polarization of angle theta.
    amplitude : Quantity
        The amplitude of the plane wave, in Volt/Meter.
    """

    attributes = ["wavelength", "amplitude", "polarization"]

    def __init__(
        self,
        wavelength: Length,
        polarization: Union[BasePolarization, Angle],
        amplitude: ElectricField,
    ):
        self.wavelength = np.atleast_1d(wavelength)
        self.amplitude = np.atleast_1d(amplitude)
        self.polarization = polarization

        if not isinstance(self.polarization, BasePolarization):
            Angle.check(self.polarization)
            self.polarization = Linear(self.polarization)

        self.binding_kwargs = dict(
            wavelength=self.wavelength,
            jones_vector=self.polarization.element,
            amplitude=self.amplitude,
            is_sequential=self.is_sequential,
        )

        self.set = PlaneWaveSourceSet(**self.binding_kwargs)
