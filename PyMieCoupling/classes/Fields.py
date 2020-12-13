import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple

from PyMieCoupling.functions.converts import rad2deg, deg2rad, NA2Angle
from PyMieCoupling.classes.Meshes import AngleMeshes, DirectMeshes



global Fontsize, pi, cmapPad
Fontsize, pi, cmapPad = 7, 3.141592, 0.2

class Source(object):

    def __init__(self,
                 Wavelength:   float,
                 Polarization: float,
                 Power:        float = 1):

        self.Wavelength = Wavelength

        self.k = 2 * pi / Wavelength

        self.Power = Power

        if Polarization != None:
            self.Polarization = Angle(Polarization)
        else:
            self.Polarization = None




class Angle(object):

    def __init__(self, input):

        if input != 'None':
            self.Degree = input
            self.Radian = deg2rad(input)
        else:
            self.Degree = 'None'
            self.Radian = 'None'



# -
