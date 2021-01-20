
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Sets import ScattererSet

LightSource = Source(Wavelength   = 950e-9,
                     Polarization = 0)

Detector0 = Photodiode(NA                = 0.1,
                       Source            = LightSource,
                       Sampling          = 1001,
                       GammaOffset       = 0,
                       PhiOffset         = 0)

Detector1 = LPmode(NA                = 0.1,
                   Source            = LightSource,
                   Sampling          = 401,
                   GammaOffset       = 0,
                   PhiOffset         = 0,
                   Mode              = (1,1))


Set = ScattererSet(DiameterList  = np.linspace(100e-9, 20000e-9, 400),
                   RIList        = np.linspace(1.4, 1.4, 1).round(1),
                   Detectors     = [Detector0, Detector1],
                   Source        = LightSource,
                   Mode          = 'Mean')


Array = Set.GetCouplingArray()


print(Array)
#DF = Set.GetCouplingFrame()

#DF.Plot(y='Coupling')

#DF.Plot(y='STD')

#plt.show()

# -
