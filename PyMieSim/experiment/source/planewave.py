#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import Union
from PyMieSim.single.source import (
    PlaneWave as _,
)  # noqa:F401  # Necessary for pybind11 binding initialization
from pydantic.dataclasses import dataclass
from TypedUnit import Length, Angle, ElectricField, AnyUnit

from PyMieSim.experiment.source.base import BaseSource
from PyMieSim.single.polarization import BasePolarization
from PyMieSim.experiment.utils import Sequential
from PyMieSim.binary.interface_sets import CppPlaneWaveSourceSet
from PyMieSim.utils import config_dict


@dataclass(config=config_dict)
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

    amplitude: ElectricField
    wavelength: Length
    polarization: Union[BasePolarization, Angle]

    def _generate_binding(self) -> None:
        """
        Prepares the keyword arguments for the C++ binding based on the scatterer's properties. This
        involves evaluating material indices and organizing them into a dictionary for the C++ interface.

        Returns
        -------
        None
        """
        self.mapping = {}

        self.binding_kwargs = dict(
            wavelength=self.wavelength,
            jones_vector=self.polarization.element,
            amplitude=self.amplitude,
            is_sequential=self.is_sequential,
        )

        self.set = CppPlaneWaveSourceSet(
            **{
                k: v.to_base_units().magnitude if isinstance(v, AnyUnit) else v
                for k, v in self.binding_kwargs.items()
            }
        )
