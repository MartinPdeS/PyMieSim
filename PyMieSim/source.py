#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from dataclasses import dataclass
from PyMieSim.polarization import LinearPolarization, JonesVector


@dataclass
class PlaneWave():
    """
    Class representing plane wave beam as a light source for light scattering.
    """
    wavelength: float
    """ Wavelenght of the light field. """
    linear_polarization: float = 0
    """ Polarization of the light field in degree. """
    amplitude: float = 1
    """ Maximal value of the electric field at focus point. """

    def __post_init__(self):
        self.k = 2 * numpy.pi / self.wavelength

        self.format_inputs()

        if isinstance(self.linear_polarization, JonesVector):
            self.linear_polarization = self.linear_polarization
        else:
            self.linear_polarization = LinearPolarization(*self.linear_polarization)

    def format_inputs(self) -> None:
        self.linear_polarization = numpy.atleast_1d(self.linear_polarization).astype(float)

# -
