#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Union
from pydantic.dataclasses import dataclass
from dataclasses import field

from PyMieSim import polarization
from PyMieSim.units import Quantity
from PyMieSim.experiment.source.base import BaseSource, config_dict

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
