
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
from PyMieCoupling.classes.Misc import Source

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

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


DataFrame = CouplingStat(RIList        = np.linspace(1.33, 1.65, 10).round(2),
                         DiameterList  = np.linspace(100,1000,20).round(2) * 1e-9,
                         Detector      = Detector0,
                         Source        = LightSource)

DataFrame.plot(y='Coupling')

DataFrame.plot(y='STD')

plt.show()















# -













# -
