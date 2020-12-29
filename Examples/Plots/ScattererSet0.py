
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Scattering import ScattererSet

LightSource = Source(Wavelength   = 950e-9,
                     Polarization = 0)

Detector = Photodiode(NA                = 0.1,
                      Source            = LightSource,
                      Samples           = 1001,
                      GammaOffset       = 0,
                      PhiOffset         = 0)

Detector1 = LPmode(NA                = 0.1,
                   Source            = LightSource,
                   Samples           = 401,
                   GammaOffset       = 0,
                   PhiOffset         = 0,
                   Mode              = (1,1))


Detector1.Plot()

Set = ScattererSet(DiameterList  = np.linspace(100e-9, 20000e-9, 400),
                   RIList        = np.linspace(1.4, 1.4, 1).round(1),
                   Detector      = Detector1,
                   Source        = LightSource,
                   Mode          = 'Mean')


DF = Set.GetCouplingFrame()

DF.Plot(y='Coupling')

#DF.Plot(y='STD')

plt.show()

# -
