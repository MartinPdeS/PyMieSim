#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from dataclasses import dataclass
from PyMieSim.polarization import LinearPolarization, JonesVector


@dataclass
class PlaneWave():
    """
    .. note::
        Class representing plane wave beam as a light source for
        light scattering.
    """
    wavelength: float
    """ Wavelenght of the light field. """
    polarization: float = 0
    """ Polarization of the light field in degree. """
    amplitude: float = 1
    """ Maximal value of the electric field at focus point. """

    def __post_init__(self):
        self.k = 2 * numpy.pi / self.wavelength
        self.amplitude = numpy.atleast_1d(self.amplitude).astype(float)
        self.polarization = numpy.atleast_1d(self.polarization).astype(float)
        self.wavelength = numpy.atleast_1d(self.wavelength).astype(float)

        if isinstance(self.polarization, JonesVector):
            self.polarization = self.polarization
        else:
            self.polarization = LinearPolarization(*self.polarization)
