
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face LP01 and LP11 Mode
For different scatterer diameters.
_________________________________________________________
"""



import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Modes import mode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.couplings import PointFieldCoupling
from PyMieCoupling.functions.Optimization import CouplingStat

npts=201

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP11 = mode(fiber       = Fiber,
            LPmode      = (1, 1),
            wavelength  = 400e-9,
            npts        = npts,
            ThetaOffset = 0,
            PhiOffset   = 0,
            Name        = 'LP11',
            Magnification=2.
            )

LP01 = mode(fiber       = Fiber,
            LPmode      = (0, 1),
            wavelength  = 400e-9,
            npts        = npts,
            ThetaOffset = 0,
            PhiOffset   = 0,
            Name        = 'LP01',
            Magnification=2.
            )

DiameterList = np.linspace(100,1000,5) * 1e-9

LP01DataFrame = CouplingStat(RIList        = [1.4],
                             DiameterList  = DiameterList,
                             Detector      = LP01,
                             wavelength    = 400e-9,
                             npts          = 101)

LP01DataFrame.plot(y = 'Coupling')

LP11DataFrame = CouplingStat(RIList       = [1.4],
                             DiameterList = DiameterList,
                             Detector     = LP11,
                             wavelength   = 400e-9,
                             npts         = 101)

LP11DataFrame.plot(y = 'Coupling')

plt.show()







# -
