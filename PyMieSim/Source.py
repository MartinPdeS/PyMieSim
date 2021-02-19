import numpy as np




class PlaneWave(object):
    def __init__(self, Wavelength, Polarization):
        self.Wavelength = Wavelength
        self.k = 2 * np.pi / Wavelength

        self.Polarization = Polarization

    def expansion(self, order):
        return (-1j)**order/(self.k*1j) * (2*order+1) / (order*(order+1));
