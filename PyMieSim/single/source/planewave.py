#!/usr/bin/env python
# -*- coding: utf-8 -*-
from TypedUnit import Length, ElectricField, Angle

from PyMieSim.single.polarization import BasePolarization
from PyMieSim.binary.interface_single import PLANEWAVE
from PyMieSim.single.source.base import BaseSource


class PlaneWave(PLANEWAVE, BaseSource):
    """
    Represents a plane wave light source for optical simulations.

    Parameters
    ----------
    wavelength : Length
        Wavelength of the light field in meters.
    polarization : BasePolarization | Angle
        Polarization state of the light field, if float is given it is assumed Linear polarization of angle theta.
    amplitude : ElectricField
        Amplitude of the electric field.
    """
    def __init__(
        self,
        wavelength: Length,
        polarization: Angle | BasePolarization,
        amplitude: ElectricField,
    ) -> None:
        """
        Represents a plane wave light source for optical simulations.

        wavelength : Length
            Wavelength of the light field in meters.
        polarization : BasePolarization | Angle
            Polarization state of the light field, if float is given it is assumed Linear polarization of angle theta.
        amplitude : ElectricField
            Amplitude of the electric field.
        """
        wavelength = Length.check(wavelength)
        amplitude = ElectricField.check(amplitude)
        self.polarization = self._validate_source_polarization(polarization)

        super().__init__(
            wavelength=wavelength,
            jones_vector=self.polarization.element[0],
            amplitude=amplitude,
        )
