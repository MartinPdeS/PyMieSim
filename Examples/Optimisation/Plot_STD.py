
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import matplotlib.pyplot as plt
import numpy as np
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.functions.Optimization import CouplingStat

npts=201

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)


Detector0 = LPmode(Fiber         = Fiber,
                   Mode          = (0, 1),
                   Wavelength    = 400e-9,
                   Npts          = npts,
                   ThetaOffset   = 0,
                   PhiOffset     = 35,
                   Magnification = 2.)


DiameterList = np.linspace(100,1000,20).round(2) * 1e-9

RIList = np.linspace(1.3, 2.0, 10).round(2)

DataFrame = CouplingStat(RIList,
                         DiameterList,
                         Detector = Detector0,
                         Wavelength = 400e-9,
                         Npts       = 101)

DataFrame.plot(y='Coupling')

DataFrame.plot(y='STD')

plt.show()















# -













# -
