
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

LightSource = Source(Wavelength   = 450e-9,
                     Polarization = 0)

Detector = Photodiode(NA                = 0.9,
                      Source            = Source,
                      Samples           = 1001,
                      GammaOffset       = 0,
                      PhiOffset         = 0)




Set = ScattererSet(DiameterList  = np.linspace(100e-9, 1000e-9, 400),
                   RIList        = np.linspace(1.3, 1.6, 4).round(1),
                   Detector      = Detector,
                   Source        = LightSource
                   )


DF = Set.GetCouplingFrame()

DF.Plot(y='Coupling')

DF.Plot(y='STD')

plt.show()

# -
