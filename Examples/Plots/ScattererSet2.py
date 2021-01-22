
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Sets import ExperimentalSet, ScattererSet

LightSource = Source(Wavelength = 405e-9, Polarization = 0)

Detector0 = Photodiode(Name              = 'Detector',
                       NA                = 0.3,
                       Source            = LightSource,
                       Sampling          = 501,
                       GammaOffset       = 0,
                       PhiOffset         = 0,
                       CouplingMode      = 'Centered')


Detector1 = Photodiode(Name              = 'Detector',
                       NA                = 0.3,
                       Source            = LightSource,
                       Sampling          = 501,
                       GammaOffset       = 0,
                       PhiOffset         = 90   ,
                       CouplingMode      = 'Centered')


Set = ExperimentalSet(DiameterList  = np.linspace(100e-9, 1000e-9, 200),
                      RIList        = np.linspace(1.5, 1.5, 1).round(1),
                      Detectors     = [Detector0, Detector1],
                      Source        = LightSource)


ScatSet = ScattererSet(DiameterList  = np.linspace(100e-9, 500e-9, 5),
                       RIList        = np.linspace(1.5, 1.5, 1).round(1),
                       Source        = LightSource)

DataFrame = Set.DataFrame

DataFrame.Plot(y='Coupling', logy=True)

plt.show()

# -
