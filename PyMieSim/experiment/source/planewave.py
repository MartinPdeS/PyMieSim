#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from pydantic.dataclasses import dataclass
from pydantic import field_validator
from PyMieSim.units import Quantity
from PyMieSim.experiment.source.base import BaseSource, config_dict


@dataclass(config=config_dict)
class PlaneWave(BaseSource):
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
        The amplitude of the plane wave, in Watts.
    """
    amplitude: Quantity

    @field_validator('amplitude', mode='before')
    def _validate_array(cls, value):
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
            jones_vector=self.polarization.element,
            amplitude=self.amplitude
        )
