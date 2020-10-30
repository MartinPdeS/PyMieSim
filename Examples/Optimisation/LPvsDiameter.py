
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face LP01 and LP11 Mode
For different scatterer diameters.
_________________________________________________________
"""

import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.Optimization import CouplingStat
from PyMieCoupling.classes.Misc import Source

npts = 151

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP11 = LPmode(Fiber         = Fiber,
              Mode          = (1, 1),
              Source        = LightSource,
              Npts          = npts,
              ThetaOffset   = 0,
              PhiOffset     = 0,
              Name          = 'LP11',
              Magnification = 2.)

LP01 = LPmode(Fiber         = Fiber,
              Mode          = (0, 1),
              Source        = LightSource,
              Npts          = npts,
              ThetaOffset   = 0,
              PhiOffset     = 0,
              Name          = 'LP01',
              Magnification = 2.)


LP01DataFrame = CouplingStat(RIList        = [1.4],
                             DiameterList  = np.linspace(100,1000,5).round(3) * 1e-9,
                             Detector      = LP01,
                             Source        = LightSource)

LP01DataFrame.plot(y = 'Coupling')

LP11DataFrame = CouplingStat(RIList        = [1.4],
                             DiameterList  = np.linspace(100,1000,5).round(3) * 1e-9,
                             Detector      = LP11,
                             Source        = LightSource)

LP11DataFrame.plot(y = 'Coupling')

plt.show()







# -
