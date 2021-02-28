from PyMieSim.Source import GaussianBeam, PlaneWave
import matplotlib.pyplot as plt
import numpy as np


""" Results shoudl be the same as ref[2] figure 2. """




beam = GaussianBeam(Wavelength   = 0.6328e-6,
                    NA           = 0.21,
                    Polarization = 0,
                    offset       = [0e-6,0e-6,0])


bsc = beam.Anm(n=2, m=-3)


print(bsc)
