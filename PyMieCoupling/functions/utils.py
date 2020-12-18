from scipy.interpolate import griddata
import numpy as np
from PyMieCoupling.functions.converts import deg2rad, rad2deg

class Angle(object):

    def __init__(self, input, unit='Degree'):
        if unit == 'Degree':
            self.Degree = input
            self.Radian = deg2rad(input)
        if unit == 'Radian':
            self.Degree = rad2deg(input)
            self.Radian = input
