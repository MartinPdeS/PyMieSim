
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.functions.Optimization import LoopRIDiameter
from PyMieCoupling.classes.Fields import Source
from PyMieCoupling.classes.DataFrame import Frame

LightSource = Source(Wavelength   = 1000e-9,
                     Polarization = 0)

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

Detector = LPmode(Fiber       = Fiber,
                  Mode        = (1, 1),
                  Source      = LightSource,
                  Npts        = 101,
                  ThetaOffset = 0,
                  PhiOffset   = 0,
                  Filter      = 0,
                  NA          = 0.1)

"""
Detector = Photodiode(NA                = 0.2,
                      Source            = LightSource,
                      Npts              = 201,
                      ThetaOffset       = 0,
                      PhiOffset         = 0)"""


RIList = np.linspace(1.3, 1.6, 1).round(4)

DiameterList = np.linspace(1000,2000,250).round(4) * 1e-9


DF = Frame(RIList        = RIList,
           DiameterList  = DiameterList,
           Detector      = Detector,
           Source        = LightSource)


DF.Plot(y='Coupling', Polarization='Perpendicular')

DF.Plot(y='Coupling', Polarization='Parallel')

DF.Plot(y='Coupling', Polarization='Filtered')

DF.Plot(y='STD', Polarization='Perpendicular')

DF.Plot(y='STD', Polarization='Parallel')

DF.Plot(y='STD', Polarization='Filtered')


plt.show()







# -
