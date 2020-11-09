
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import matplotlib.pyplot as plt
import numpy as np
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.classes.DataFrame import Frame
from PyMieCoupling.classes.Misc import Source

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

npts = 51
cuda = False

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)


Detector0 = LPmode(Fiber         = Fiber,
                   Mode          = (0, 1),
                   Source        = LightSource,
                   Npts          = npts,
                   ThetaOffset   = 0,
                   PhiOffset     = 35,
                   Magnification = 2.,
                   cuda          = cuda)


DataFrame = Frame(RIList        = np.linspace(1.33, 1.65, 6).round(2),
                  DiameterList  = np.linspace(100,1000,100).round(2) * 1e-9,
                  Detector      = Detector0,
                  Source        = LightSource,
                  cuda          = cuda)

#DataFrame.plot(y='Coupling')

#DataFrame.plot(y='STD')

#plt.show()















# -













# -
