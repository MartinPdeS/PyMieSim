import numpy as np
from typing import Tuple
from PyMieCoupling.functions.converts import deg2rad



class Source(object):

    def __init__(self,
                 Wavelength: float,
                 Polarization: float):

        self.Wavelength = Wavelength

        self.Polarization = Polarization

        self.k = 2 * np.pi / Wavelength
