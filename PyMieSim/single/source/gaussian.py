#!/usr/bin/env python
# -*- coding: utf-8 -*-
from TypedUnit import Length, Angle, Power, Dimensionless

from PyMieSim.single.polarization import BasePolarization
from PyMieSim.binary.interface_single import GAUSSIAN
from PyMieSim.single.source.base import BaseSource


class Gaussian(GAUSSIAN, BaseSource):
    """
    Represents a Gaussian light source for optical simulations.
    """
    def __init__(
        self,
        wavelength: Length,
        polarization: Angle | BasePolarization,
        optical_power: Power,
        NA: Dimensionless,
    ) -> None:
        """
        Represents a Gaussian light source for optical simulations.

        wavelength: Length
            Wavelength of the light field in meters.
        polarization : Angle | BasePolarization
            Polarization state of the light field, if float is given it is assumed Linear polarization of angle theta.
        optical_power: Power
            Optical power of the source in Watts.
        NA: Dimensionless
            Numerical aperture of the source.
        """
        NA = Dimensionless.check(NA)
        wavelength = Length.check(wavelength)
        optical_power = Power.check(optical_power)

        self.polarization = self._validate_source_polarization(polarization)

        super().__init__(
            wavelength=wavelength,
            jones_vector=self.polarization.element[0],
            optical_power=optical_power,
            NA=NA,
        )

# -
