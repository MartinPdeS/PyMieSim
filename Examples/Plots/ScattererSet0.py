
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

LightSource = Source(Wavelength   = 950e-9,
                     Polarization = 0)

Detector0 = LPmode(NA                = 0.05,
                   Source            = LightSource,
                   Sampling          = 401,
                   GammaOffset       = 40,
                   PhiOffset         = 0,
                   Mode              = (0,1),
                   CouplingMode      = 'Mean')

Detector1 = LPmode(NA                = 0.05,
                   Source            = LightSource,
                   Sampling          = 401,
                   GammaOffset       = 40,
                   PhiOffset         = 0,
                   Mode              = (1,1),
                   CouplingMode      = 'Mean')


Set = ExperimentalSet(DiameterList  = np.linspace(100e-9, 5000e-9, 400),
                      RIList        = np.linspace(1.5, 1.5, 1).round(1),
                      Detectors     = [Detector0, Detector1],
                      Source        = LightSource)



ScatSet = ScattererSet(DiameterList  = np.linspace(100e-9, 500e-9, 5),
                       RIList        = np.linspace(1.5, 1.5, 1).round(1),
                       Source        = LightSource)



DataFrame = ScatSet.S1S2

DataFrame.Plot(y='S2')

Array = Set.Coupling

DataFrame = Set.DataFrame


DataFrame.Plot(y='Coupling')


plt.show()

# -
